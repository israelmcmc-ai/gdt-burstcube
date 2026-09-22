"""BurstCube TTE (time-tagged event) data: a 1024-channel photon list for one
detector, built from the archive's ``events/*_tte_uf.evt.gz`` files.

The ``EVENTS`` extension's own ``TSTART``/``TSTOP`` keywords are unusable in
**every** TTE file in the archive. All 28 store them as strings rather than
numbers, and in all 28 ``TSTOP - TSTART`` disagrees with the file's own
``ONTIME``. In 16 of them ``TSTOP < TSTART``, which drives ``TELAPSE`` and
``EXPOSURE`` negative and leaves the primary ``DATE-END`` earlier than
``DATE-OBS``. This is not described in the official archive caveats document
-- caveat #4 concerns how few TTE files were downlinked and caveat #5 the
2024-08-15 timeline alignment, neither the header keywords. See the README's
Caveats section for the per-file breakdown.

This reader therefore never uses those two keywords. The GTI comes from the
``STDGTI`` extension and the time range from the event times, which agree
with each other exactly in all 28 files and with ``ONTIME`` to better than a
millisecond. A `UserWarning` is raised when the keywords are self-
contradictory, but nothing crashes and no negative exposure is produced.

There is no ``Rsp.to_tte()``: instead, :meth:`BurstCubeTTE.to_64_channels`
applies the CALDB ``reb64`` grouping (1024 -> 64 channels) so that TTE data
can be folded against the native 64-channel detector response.
"""
import warnings

import numpy as np
import astropy.io.fits as fits

from gdt.core.tte import PhotonList
from gdt.core.data_primitives import Gti, EventList

from . import caldb
from .headers import TTEHeaders
from .time import check_met_epoch

__all__ = ['BurstCubeTTE', 'DEFAULT_GAP_THRESHOLD']

#: The default gap threshold for :meth:`BurstCubeTTE.recording_blocks`,
#: in seconds.
DEFAULT_GAP_THRESHOLD = 0.5


class BurstCubeTTE(PhotonList):
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
            (:class:`BurstCubeTTE`)

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

        # The result covers the original GTI intersected with the UNION of
        # the kept ranges. Accumulating by repeated intersection instead --
        # as gdt-core's PhotonList.slice_time does -- empties the GTI as soon
        # as two ranges are disjoint, and the next intersection then raises
        # on the empty result.
        #
        # Ranges containing no events are dropped: an empty EventList has a
        # time_range of None, which Gti.from_list cannot take, and
        # EventList.merge reduces over each segment's times and so fails on
        # an empty one. Dropping them is also what a user means -- TTE
        # coverage is clustered, so a reasonable set of windows can leave
        # some empty.
        populated = [segment for segment in segments if segment.size > 0]
        if not populated:
            raise ValueError(
                'None of the requested time ranges contain any events; the '
                f'data spans {self.time_range[0]} to {self.time_range[1]}.')

        gti = Gti.intersection(
            Gti.from_list(self.gti.as_list()),
            Gti.from_list([segment.time_range for segment in populated]))

        data = EventList.merge(populated, sort=True, force_unique=False)

        headers = self._build_headers(self.trigtime, *data.time_range,
                                      data.num_chans)
        return self.from_data(data, gti=gti, trigger_time=self.trigtime,
                              headers=headers,
                              event_deadtime=self.event_deadtime,
                              overflow_deadtime=self.overflow_deadtime,
                              **kwargs)

    def recording_blocks(self, gap_threshold=DEFAULT_GAP_THRESHOLD):
        """The intervals over which this file actually recorded events.

        TTE does not cover its nominal span continuously: event times arrive
        in short blocks separated by gaps of comparable length, and across a
        gap the event list is simply empty. CBD -- the same detector and the
        same event stream, binned on board rather than written out event by
        event -- counts at its normal rate through those gaps, which is how
        we know they are not quiet sky. See the README caveat *TTE gaps*;
        ``examples/tte_gaps_vs_cbd.py`` reproduces the measurement behind
        it.

        Every rate derived from TTE needs these blocks, because dividing a
        count by an elapsed duration that spans a gap dilutes the rate by
        the gap's share of that duration. On 2024-08-14 ``CS0`` the live
        time is 97 s of a 224 s span, a factor of 2.3.

        Each block runs from its first to its last event, so the returned
        total is live time and not elapsed time. The gaps themselves are the
        complement of the result over the file's time range, which
        :meth:`~gdt.missions.burstcube.gti.BurstCubeGTI.complement` will
        give you.

        Args:
            gap_threshold (float, optional): Start a new block wherever the
                wait for the next event exceeds this, in seconds. The
                default of 0.5 s sits far from both scales it separates --
                in-block waits are milliseconds and real gaps are around a
                second -- so the result is insensitive to the exact value.

        Returns:
            (:class:`~gdt.core.data_primitives.Gti`)

        Raises:
            ValueError: If there are no events, since there are no blocks to
                report, or if ``gap_threshold`` is not positive.
        """
        if gap_threshold <= 0.0:
            raise ValueError(f'gap_threshold must be positive, got {gap_threshold}')

        times = np.sort(self.data.times)
        if times.size == 0:
            raise ValueError('There are no events, so there are no recording '
                             'blocks.')

        gap_idx = np.flatnonzero(np.diff(times) > gap_threshold)
        starts = times[np.concatenate(([0], gap_idx + 1))]
        stops = times[np.concatenate((gap_idx, [times.size - 1]))]
        return Gti.from_bounds(starts, stops)

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
            (:class:`BurstCubeTTE`)
        """
        obj = super().open(file_path, **kwargs)

        hdrs = [hdu.header for hdu in obj.hdulist]
        headers = TTEHeaders.from_headers(hdrs)
        check_met_epoch(headers['EVENTS'])

        events_idx = obj.hdu_index_from_name('EVENTS')
        gti_idx = obj.hdu_index_from_name('STDGTI')

        detector = headers['EVENTS']['INSTRUME']
        cls._warn_if_tstart_tstop_broken(headers['EVENTS'])

        ebounds = caldb.ebounds(detector, 1024)
        times = obj.column(events_idx, 'TIME')
        channels = obj.column(events_idx, 'PHA')
        data = EventList(times=times, channels=channels, ebounds=ebounds)

        # EVENTS' own TSTART/TSTOP are unusable in every archive TTE file
        # (see the module docstring), so the GTI always comes from STDGTI's
        # own START/STOP columns instead.
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
            (:class:`BurstCubeTTE`)
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
        broken EVENTS TSTART/TSTOP: TSTART/TSTOP stored as
        strings and/or TSTOP < TSTART, which would otherwise propagate a
        negative TELAPSE/EXPOSURE.

        Args:
            events_header (:class:`~gdt.core.headers.Header`): The EVENTS
                header, after :class:`~.headers.TTEHeaders` coercion (so
                TSTART/TSTOP are floats here even if stored as strings on
                disk).
        """
        tstart = events_header['TSTART']
        tstop = events_header['TSTOP']
        if tstop < tstart:
            warnings.warn(
                'This TTE file has TSTOP < TSTART in its EVENTS header '
                '; ignoring both and using the STDGTI '
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
