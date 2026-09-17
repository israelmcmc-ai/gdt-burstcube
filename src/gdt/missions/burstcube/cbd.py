"""BurstCube CBD (continuous binned data): a 16-channel time history of
counts for one detector.

Two critical things this reader gets right that a naive port would not:

* ``TIMEPIXR=1`` means the single ``TIME`` column is the *end* of each bin,
  not the start. Bin edges are ``[TIME - TIMEDEL, TIME]``. Getting this
  backwards silently shifts every light curve by one bin width (0.256 s).
* The time grid is not uniform: the instrument was frequently off or
  unstable, leaving gaps of seconds to tens of minutes between bins (more so
  in ``_cl`` than ``_uf``). Bin edges are built from the actual ``TIME``
  values, and gdt-core's :class:`~gdt.core.data_primitives.TimeEnergyBins`
  already rebins/slices each contiguous run separately, so a gap is never
  bridged or filled.
"""
import warnings

import numpy as np
import astropy.io.fits as fits

from gdt.core.phaii import Phaii
from gdt.core.data_primitives import Gti, TimeEnergyBins

from . import caldb
from .headers import CBDHeaders, CBDUnfilteredHeaders

__all__ = ['BurstCubeCBD']

#: Tolerance, in seconds, separating float64 rounding artifacts in the bin
#: timestamps from real sub-TIMEDEL structure. The artifacts sit at one ULP
#: (~1.5e-08 s at BurstCube MET magnitudes) and the smallest real feature is
#: 1.0e-03 s, so anything in between works; this sits ~67x above the noise and
#: ~1000x below the feature. See the discussion in :meth:`BurstCubeCBD.open`.
_SNAP_TOL = 1e-6


class BurstCubeCBD(Phaii):
    """BurstCube continuous binned data (CBD) for one detector: a 16-channel
    time history of counts, built from the archive's ``monitor/*_3cbd_{uf,cl}.fits.gz``
    files.

    The auxiliary per-bin diagnostic columns (:attr:`original_time`,
    :attr:`time_syst_error`, :attr:`sumtot`, :attr:`sumch_2_15`,
    :attr:`blended`) are only populated when the object comes directly from
    :meth:`open`. They have no well-defined meaning after :meth:`slice_time`,
    :meth:`rebin_time`, or :meth:`merge` (which combine or drop bins), so
    those operations return ``None`` for all five.
    """

    def __init__(self):
        super().__init__()
        self._original_time = None
        self._time_syst_error = None
        self._sumtot = None
        self._sumch_2_15 = None
        self._blended = None

    @property
    def detector(self):
        """(str): The detector name (e.g. ``'CS0'``), from ``INSTRUME``."""
        return self.headers[0]['INSTRUME']

    @property
    def original_time(self):
        """(numpy.ndarray or None): The raw, uncorrected time from telemetry
        for each bin (the ``ORIGINAL_TIME`` column). See archive caveat #11.
        """
        return self._original_time

    @property
    def time_syst_error(self):
        """(numpy.ndarray or None): The systematic uncertainty on ``TIME``
        for each bin, in seconds (the ``TIME_SYST_ERROR`` column).
        """
        return self._time_syst_error

    @property
    def sumtot(self):
        """(numpy.ndarray or None): The sum of counts over all 16 channels
        for each bin, as recorded in the file (the ``SUMTOT`` column) --
        independent of, and a check against, summing :attr:`data`'s counts.

        ``None`` for unfiltered (``_uf``) files, which do not carry this
        column; only the cleaned (``_cl``) files do.
        """
        return self._sumtot

    @property
    def sumch_2_15(self):
        """(numpy.ndarray or None): The sum of counts over channels 2-15 for
        each bin, as recorded in the file (the ``SUMCH_2_15`` column).

        ``None`` for unfiltered (``_uf``) files, which do not carry this
        column; only the cleaned (``_cl``) files do.
        """
        return self._sumch_2_15

    @property
    def blended(self):
        """(numpy.ndarray of bool or None): True for a bin whose ``TIME`` is
        less than ``TIMEDEL`` after the previous bin's ``TIME``, i.e. the
        packet cadence appears faster than the nominal 0.256 s bin width.
        Per archive caveat #9, such bins may blend telemetry from different
        periods and should be treated with extreme caution; they are the
        bins filtered out of ``_cl``. The first bin is never flagged (there
        is no preceding bin to compare against).
        """
        return self._blended

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a BurstCube CBD FITS file.

        Args:
            file_path (str): The file path of the FITS file

        Returns:
            (:class:`BurstCubeCBD`)
        """
        obj = super().open(file_path, **kwargs)

        hdrs = [hdu.header for hdu in obj.hdulist]

        # The STDGTI extension comes in two schemas (see CBDGtiHeader): the
        # 2-column standard form and the 6-column form shared with TTE. Pick
        # by looking at the extension, not at the _uf/_cl filename, since
        # nothing in the files guarantees that correlation holds everywhere.
        gti_hdr = obj.hdulist['STDGTI'].header
        headers_cls = (CBDHeaders if 'HDUCLAS2' in gti_hdr
                       else CBDUnfilteredHeaders)
        headers = headers_cls.from_headers(hdrs)

        cbd_idx = obj.hdu_index_from_name('CBD')
        gti_idx = obj.hdu_index_from_name('STDGTI')

        detector = headers['CBD']['INSTRUME']
        timedel = headers['CBD']['TIMEDEL']

        # TIMEPIXR=1: TIME is the END of the bin.
        time = obj.column(cbd_idx, 'TIME')
        tstart = time - timedel
        tstop = time

        # TIME is ~1e8 (MET), where a float64 ULP is ~1.5e-8 s: a
        # TIME - TIMEDEL subtraction can therefore land a ULP away from the
        # *previous* row's TIME even when the two bins are genuinely
        # back-to-back. gdt-core's contiguous-segment detection (used by
        # rebin_time/slice_time) requires an *exact* match to recognize
        # that, so snap tstart to the previous tstop where they already
        # agree to within the noise floor.
        #
        # The tolerance has to clear the noise but stay well below the
        # smallest REAL feature. Measured on bc240530cs0_3cbd_uf.fits.gz,
        # the two populations are cleanly separated by five orders of
        # magnitude:
        #   * float64 artifacts: |tstart[i] - tstop[i-1]| == 1.49e-08 s
        #     (exactly 1 ULP), 17989 of 98187 boundaries;
        #   * genuine short bins: 1.000e-03 s short of TIMEDEL, 432 of them.
        # _SNAP_TOL sits between the two. A tolerance of 1e-3 would cut
        # straight through the second population, snapping the 136 whose
        # noise puts them just under 1 ms and leaving the other 296 alone --
        # treating one physical feature two different ways depending on
        # rounding.
        close = np.isclose(tstart[1:], tstop[:-1], atol=_SNAP_TOL, rtol=0.0)
        tstart[1:][close] = tstop[:-1][close]

        # The 1 ms short bins leave tstart[i] genuinely *before* tstop[i-1],
        # so consecutive bins overlap by a millisecond. A counting bin cannot
        # overlap its predecessor, and gdt-core refuses to merge overlapping
        # bins outright -- rebin_time raises ValueError on every real CBD file
        # (53 overlaps in one day of CS0 _cl, 432 in _uf) while slice_time
        # happens to succeed, so the failure only shows up on rebinning.
        #
        # Clamp the start to the previous stop. Those bins then carry a
        # 0.255 s exposure rather than the nominal TIMEDEL, which is the
        # honest reading: the packet arrived a millisecond early, and the bin
        # covers only the time since the previous one ended. Nothing is
        # hidden -- the `blended` flag below still marks every one of them.
        np.maximum(tstart[1:], tstop[:-1], out=tstart[1:])

        exposure = tstop - tstart
        if np.any(exposure <= 0.0):
            warnings.warn(
                f'{np.sum(exposure <= 0.0)} CBD bins have a non-positive '
                'exposure after clamping overlaps, meaning consecutive TIME '
                'values are not increasing by less than TIMEDEL. The bin '
                'edges for those rows should not be trusted.',
                RuntimeWarning, stacklevel=2)

        ebounds = caldb.ebounds(detector, 16)
        data = TimeEnergyBins(obj.column(cbd_idx, 'COUNTS'), tstart, tstop,
                              exposure, ebounds.low_edges(),
                              ebounds.high_edges())

        # Caveat #9: CBD recorded with a packet cadence < TIMEDEL may blend
        # telemetry from different periods and "should be treated with
        # extreme caution". The flag applies that threshold literally,
        # backing off only by the float64 noise floor -- not by enough to
        # also swallow the 1 ms short bins, which are a real feature of the
        # timestamps and are therefore reported rather than hidden. The
        # first bin has no predecessor to compare against.
        dt = np.diff(time)
        blended = np.concatenate(([False], dt < (timedel - _SNAP_TOL)))

        # the STDGTI schema differs between _cl (2 columns) and _uf (6
        # columns); START/STOP are common to both, so read only those.
        gti = Gti.from_bounds(obj.column(gti_idx, 'START'),
                              obj.column(gti_idx, 'STOP'))

        original_time = obj.column(cbd_idx, 'ORIGINAL_TIME')
        time_syst_error = obj.column(cbd_idx, 'TIME_SYST_ERROR')

        # SUMTOT and SUMCH_2_15 are present only in the cleaned (_cl) files;
        # the unfiltered (_uf) files carry TIME, COUNTS, ORIGINAL_TIME and
        # TIME_SYST_ERROR only. Verified across several days and detectors.
        # They stay None when absent -- the same sums are always available
        # from the COUNTS array anyway.
        present = obj.hdulist[cbd_idx].columns.names
        sumtot = obj.column(cbd_idx, 'SUMTOT') if 'SUMTOT' in present else None
        sumch_2_15 = (obj.column(cbd_idx, 'SUMCH_2_15')
                      if 'SUMCH_2_15' in present else None)

        obj.close()

        new_obj = cls.from_data(data, gti=gti, filename=obj.filename,
                                headers=headers)
        new_obj._original_time = original_time
        new_obj._time_syst_error = time_syst_error
        new_obj._sumtot = sumtot
        new_obj._sumch_2_15 = sumch_2_15
        new_obj._blended = blended
        return new_obj

    def _build_hdulist(self):
        hdulist = fits.HDUList()

        primary_hdu = fits.PrimaryHDU(header=self.headers['PRIMARY'])
        for key, val in self.headers['PRIMARY'].items():
            primary_hdu.header[key] = val
        hdulist.append(primary_hdu)

        hdulist.append(self._cbd_table())
        hdulist.append(self._gti_table())
        return hdulist

    def _build_headers(self, trigtime, tstart, tstop, num_chans):
        headers = self.headers.copy()
        for hdu in headers:
            # PRIMARY carries neither TSTART/TSTOP nor DETCHANS; only the
            # CBD and STDGTI extensions do.
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

    def _cbd_table(self):
        time = np.copy(self.data.tstop)
        counts_col = fits.Column(name='COUNTS', format='16I', unit='count',
                                 array=self.data.counts)
        time_col = fits.Column(name='TIME', format='1D', unit='s', array=time)
        hdu = fits.BinTableHDU.from_columns([time_col, counts_col],
                                            header=self.headers['CBD'])
        for key, val in self.headers['CBD'].items():
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
