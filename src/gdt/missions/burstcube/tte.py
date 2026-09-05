# Copyright 2024-2025 by the BurstCube Team.
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
"""BurstCube TTE (time-tagged event) data: a 1024-channel photon list for one
detector, built from the archive's ``events/*_tte_uf.evt.gz`` files.

Per archive caveat #4, the ``EVENTS`` extension's own ``TSTART``/``TSTOP``
keywords are unreliable (sometimes written as strings, sometimes with
``TSTOP < TSTART``, which drives ``TELAPSE``/``EXPOSURE`` negative). This
reader never uses those two keywords to build the event data or its GTI --
the GTI always comes from the ``STDGTI`` extension, and the time range comes
from the event times themselves, which is the fallback the spec calls for.
A `UserWarning` is raised if the broken condition is detected, but nothing
crashes and no negative exposure is ever produced.

There is no ``Rsp.to_tte()``: instead, :meth:`BurstCubeTte.to_64_channels`
applies the CALDB ``reb64`` grouping (1024 -> 64 channels) so that TTE data
can be folded against the native 64-channel detector response.
"""
import warnings

import numpy as np
import astropy.io.fits as fits

from gdt.core.tte import PhotonList
from gdt.core.data_primitives import Gti, EventList

from . import caldb
from .headers import TteHeaders

__all__ = ['BurstCubeTte']


class BurstCubeTte(PhotonList):
    """BurstCube time-tagged event (TTE) data for one detector: a 1024-channel
    photon list, with a full CALDB energy calibration (``eb1024``).
    """

    def slice_time(self, time_ranges, **kwargs):
        """Slice the photon list to one or more time ranges.

        This overrides :meth:`gdt.core.tte.PhotonList.slice_time` to avoid
        losing counts. The base implementation reassembles the sliced
        segments with ``EventList.merge(..., force_unique=True)``, which
        treats two events sharing a timestamp as duplicates and keeps only
        one. BurstCube records event times at a 10 microsecond granularity
        (``TIMEDEL = 1e-05``), coarse enough that many genuinely distinct
        photons share a timestamp: in ``bc240530cs0_tte_uf.evt.gz``, 6364
        events carry only 1644 distinct times, so the base method silently
        discards about 74% of the counts.

        Deduplication only exists to handle segments that overlap in time.
        This implementation rejects overlapping ranges outright and then
        merges with ``force_unique=False``, so every event in the requested
        ranges is kept.

        Args:
            time_ranges ([(float, float), ...]): The time ranges to slice to

        Returns:
            (:class:`BurstCubeTte`)

        Raises:
            ValueError: If any two of ``time_ranges`` overlap, since the
                result would double-count the overlap.
        """
        time_ranges = self._assert_range_list(time_ranges)
        ordered = sorted(self._assert_range(r) for r in time_ranges)
        for (_, prev_stop), (next_start, _) in zip(ordered, ordered[1:]):
            if next_start < prev_stop:
                raise ValueError(
                    'time_ranges must not overlap; got overlapping ranges '
                    f'ending at {prev_stop} and starting at {next_start}')

        segments = [self.data.time_slice(*r) for r in ordered]

        gti = Gti.from_list(self.gti.as_list())
        for segment in segments:
            gti = Gti.intersection(gti, Gti.from_list([segment.time_range]))

        data = EventList.merge(segments, sort=True, force_unique=False)

        headers = self._build_headers(self.trigtime, *data.time_range,
                                      data.num_chans)
        return self.from_data(data, gti=gti, trigger_time=self.trigtime,
                              headers=headers,
                              event_deadtime=self.event_deadtime,
                              overflow_deadtime=self.overflow_deadtime,
                              **kwargs)

    @property
    def detector(self):
        """(str): The detector name (e.g. ``'CS0'``), from ``INSTRUME``."""
        return self.headers[0]['INSTRUME']

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a BurstCube TTE FITS file.

        Args:
            file_path (str): The file path of the FITS file

        Returns:
            (:class:`BurstCubeTte`)
        """
        obj = super().open(file_path, **kwargs)

        hdrs = [hdu.header for hdu in obj.hdulist]
        headers = TteHeaders.from_headers(hdrs)

        events_idx = obj.hdu_index_from_name('EVENTS')
        gti_idx = obj.hdu_index_from_name('STDGTI')

        detector = headers['EVENTS']['INSTRUME']
        cls._warn_if_tstart_tstop_broken(headers['EVENTS'])

        ebounds = caldb.ebounds(detector, 1024)
        times = obj.column(events_idx, 'TIME')
        channels = obj.column(events_idx, 'PHA')
        data = EventList(times=times, channels=channels, ebounds=ebounds)

        # per archive caveat #4, EVENTS' own TSTART/TSTOP are unreliable;
        # the GTI always comes from STDGTI's own START/STOP columns instead.
        # the _cl/_uf CBD schema difference (round-1 finding) does not apply
        # here: TTE STDGTI is always the 6-column schema, but START/STOP are
        # common to both, so reading just those two columns is safe either way.
        gti = Gti.from_bounds(obj.column(gti_idx, 'START'),
                              obj.column(gti_idx, 'STOP'))

        obj.close()

        return cls.from_data(data, gti=gti, filename=obj.filename,
                             headers=headers)

    def to_64_channels(self):
        """Regroup the 1024-channel data to the native 64-channel scheme,
        using the CALDB ``reb64`` grouping for this detector, so that TTE can
        be folded against the 64-channel detector response
        (``BurstCubeRsp``, added in a later version of this plugin -- there
        is no ``Rsp.to_tte()``, this is the intended path in the other
        direction).

        Returns:
            (:class:`BurstCubeTte`)
        """
        reb = caldb.rebin(self.detector, 64)
        # reb64's groups are contiguous (chan_min[i] == chan_max[i-1] + 1),
        # so the channel edges are just the group starts plus one final edge.
        new_edges = np.append(reb.chan_min, reb.chan_max[-1] + 1)

        def _group_by_caldb_edges(counts, count_uncert, exposure, old_edges):
            return counts, count_uncert, exposure, new_edges

        return self.rebin_energy(_group_by_caldb_edges)

    @staticmethod
    def _warn_if_tstart_tstop_broken(events_header):
        """Warn (never raise) if this file exhibits the archive's known
        broken EVENTS TSTART/TSTOP (caveat #4): TSTART/TSTOP stored as
        strings and/or TSTOP < TSTART, which would otherwise propagate a
        negative TELAPSE/EXPOSURE.

        Args:
            events_header (:class:`~gdt.core.headers.Header`): The EVENTS
                header, after :class:`~.headers.TteHeaders` coercion (so
                TSTART/TSTOP are floats here even if stored as strings on
                disk).
        """
        tstart = events_header['TSTART']
        tstop = events_header['TSTOP']
        if tstop < tstart:
            warnings.warn(
                'This TTE file has TSTOP < TSTART in its EVENTS header '
                '(archive caveat #4); ignoring both and using the STDGTI '
                'extension and the event times instead.', UserWarning,
                stacklevel=3)

    def _build_hdulist(self):
        hdulist = fits.HDUList()

        primary_hdu = fits.PrimaryHDU(header=self.headers['PRIMARY'])
        for key, val in self.headers['PRIMARY'].items():
            primary_hdu.header[key] = val
        hdulist.append(primary_hdu)

        hdulist.append(self._events_table())
        hdulist.append(self._gti_table())
        return hdulist

    def _build_headers(self, trigtime, tstart, tstop, num_chans):
        headers = self.headers.copy()
        for hdu in headers:
            try:
                hdu['TSTART'] = tstart
                hdu['TSTOP'] = tstop
            except KeyError:
                pass
            try:
                hdu['DETCHANS'] = num_chans
            except KeyError:
                pass
        return headers

    def _events_table(self):
        time_col = fits.Column(name='TIME', format='1D', unit='s',
                               array=np.copy(self.data.times))
        pha_col = fits.Column(name='PHA', format='1I',
                              array=self.data.channels)
        hdu = fits.BinTableHDU.from_columns([time_col, pha_col],
                                            header=self.headers['EVENTS'])
        for key, val in self.headers['EVENTS'].items():
            hdu.header[key] = val
        return hdu

    def _gti_table(self):
        start_col = fits.Column(name='START', format='1D', unit='s',
                                array=np.array(self.gti.low_edges()))
        stop_col = fits.Column(name='STOP', format='1D', unit='s',
                               array=np.array(self.gti.high_edges()))
        hdu = fits.BinTableHDU.from_columns([start_col, stop_col],
                                            header=self.headers['STDGTI'])
        for key, val in self.headers['STDGTI'].items():
            hdu.header[key] = val
        return hdu
