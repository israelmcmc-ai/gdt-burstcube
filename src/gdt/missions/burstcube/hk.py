"""BurstCube detector housekeeping: the ``DETECTOR_HK1`` (1 s cadence) and
``DETECTOR_HK2`` (60 s cadence) extensions of ``auxil/bcYYMMDDcsa.hk.gz``,
each exposed as an :class:`~astropy.timeseries.TimeSeries`.

Per archive caveat #7 the per-detector energy thresholds were raised
mid-mission to suppress a non-Poissonian low-energy noise component.
:meth:`~BurstCubeHK.base_threshold` is how that shows up here.
:meth:`~BurstCubeHK.peak_threshold` reads the other threshold column, which
is a pulse-shape cut rather than an energy one and never moves; both
docstrings spell out the difference.

The keV figures caveat #7 quotes are on a different energy calibration from
CALDB's, and read ~20% high against the archive's own scale -- see the
README caveat *CALDB energy scale disagrees with caveat #7's threshold
table*.
"""
from astropy.timeseries import TimeSeries

from gdt.core.file import FitsFileContextManager

from .headers import DetectorHKHeaders
from .time import Time

__all__ = ['BurstCubeHK']

# the 4-element housekeeping arrays are indexed by detector number 0..3,
# matching BurstCubeDetectors.number
_ENABLE_FLAG_COLUMNS = ('TRIG_ENABLED', 'DET_ENABLED', 'CBD_ENABLED',
                        'TTE_ENABLED')
_THRESHOLD_COLUMNS = ('PEAK_THRES', 'BASE_THRES')


class BurstCubeHK(FitsFileContextManager):
    """Reader for a BurstCube detector housekeeping file."""

    @property
    def hk1(self) -> TimeSeries:
        """(astropy.timeseries.TimeSeries): The ``DETECTOR_HK1`` data
        (1 s cadence): bias voltages/currents, temperatures, rate counts,
        the enable flags, and more."""
        return self._series('DETECTOR_HK1')

    @property
    def hk2(self) -> TimeSeries:
        """(astropy.timeseries.TimeSeries): The ``DETECTOR_HK2`` data
        (60 s cadence): the energy thresholds and bias/temperature control
        settings."""
        return self._series('DETECTOR_HK2')

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a BurstCube detector housekeeping FITS file.

        Args:
            file_path (str): The file path of the FITS file

        Returns:
            (:class:`BurstCubeHK`)
        """
        obj = super().open(file_path, **kwargs)
        hdrs = [hdu.header for hdu in obj.hdulist]
        obj._headers = DetectorHKHeaders.from_headers(hdrs)
        return obj

    def trig_enabled(self, detector=None):
        """The ``TRIG_ENABLED`` flags from ``DETECTOR_HK1``: whether the
        onboard trigger algorithm was enabled for a detector at each time.

        Args:
            detector (int, str, or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
                The detector, e.g. ``0`` or ``BurstCubeDetectors.CS0``. If
                omitted, returns the full (`n`, 4) array for all detectors.

        Returns:
            (numpy.ndarray)
        """
        return self._enable_flag('TRIG_ENABLED', detector)

    def det_enabled(self, detector=None):
        """The ``DET_ENABLED`` flags from ``DETECTOR_HK1``: whether a
        detector itself was enabled at each time.

        Args:
            detector (int, str, or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
                The detector. If omitted, returns the full (`n`, 4) array.

        Returns:
            (numpy.ndarray)
        """
        return self._enable_flag('DET_ENABLED', detector)

    def cbd_enabled(self, detector=None):
        """The ``CBD_ENABLED`` flags from ``DETECTOR_HK1``: whether CBD data
        collection was enabled for a detector at each time.

        Args:
            detector (int, str, or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
                The detector. If omitted, returns the full (`n`, 4) array.

        Returns:
            (numpy.ndarray)
        """
        return self._enable_flag('CBD_ENABLED', detector)

    def tte_enabled(self, detector=None):
        """The ``TTE_ENABLED`` flags from ``DETECTOR_HK1``: whether TTE data
        collection was enabled for a detector at each time.

        Args:
            detector (int, str, or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
                The detector. If omitted, returns the full (`n`, 4) array.

        Returns:
            (numpy.ndarray)
        """
        return self._enable_flag('TTE_ENABLED', detector)

    def peak_threshold(self, detector=None):
        """The ``PEAK_THRES`` column of ``DETECTOR_HK2``, in mV: the
        pulse-shape criterion the IDAB front end uses to accept a sample as
        a genuine peak.

        The file's own column comment is "Min lvl past pre-post window x
        valid IDAB peak" -- how far a candidate sample must stand above the
        samples on either side of it, within the ``PHA_WIN`` window, to
        count as a peak rather than a shoulder or a slow drift. It is a
        pulse-shape cut, not an energy cut.

        This is therefore **not** the threshold archive caveat #7's
        mid-mission change moved -- that is :meth:`base_threshold`, which
        measures height above the running baseline instead. ``PEAK_THRES``
        reads 6.0 mV for every detector in every housekeeping file in the
        archive and never changes.

        Neither column is described in the archive caveats document or the
        archive guide; the ``TTYPE`` comments above are the only
        documentation, and there is no equivalent in any Fermi GBM data
        product.

        Args:
            detector (int, str, or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
                The detector. If omitted, returns the full (`n`, 4) array.

        Returns:
            (numpy.ndarray)
        """
        return self._threshold('PEAK_THRES', detector)

    def base_threshold(self, detector=None):
        """The ``BASE_THRES`` discriminator threshold from ``DETECTOR_HK2``,
        in mV: how far a pulse must rise above the running baseline average
        to be recorded. The file's own column comment is "Min lvl past
        baseline avg x valid IDAB peak".

        This is the energy threshold archive caveat #7's mid-mission change
        moved: 82 -> 238 mV on ``CS0``, and 98 -> 257, 74 -> 247, 74 -> 259
        on ``CS1``-``CS3``. The new values were tried temporarily before
        being made permanent, so a single day's file can hold both.

        Caveat #7's Table 1 quotes these in keV (26.93 -> 100.12 keV on
        detector 0). The turn-on measured in the TTE spectra sits ~20% below
        that, consistently across all four detectors and on both sides of
        the change, because the archive's PHA-to-energy mapping is CALDB's
        pre-launch ``eb1024`` and the table is on a different calibration.
        Use the CALDB scale for anything you compare against other archive
        energies; see notebook 1.

        Args:
            detector (int, str, or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
                The detector. If omitted, returns the full (`n`, 4) array.

        Returns:
            (numpy.ndarray)
        """
        return self._threshold('BASE_THRES', detector)

    def _series(self, ext_name: str) -> TimeSeries:
        idx = self.hdu_index_from_name(ext_name)
        columns = [c for c in self.get_column_names(idx) if c != 'TIME']
        data = {name: self.column(idx, name) for name in columns}
        time = Time(self.column(idx, 'TIME'), format='burstcube')
        return TimeSeries(time=time, data=data)

    def _enable_flag(self, column: str, detector):
        idx = self.hdu_index_from_name('DETECTOR_HK1')
        return self._detector_column(idx, column, detector)

    def _threshold(self, column: str, detector):
        idx = self.hdu_index_from_name('DETECTOR_HK2')
        return self._detector_column(idx, column, detector)

    def _detector_column(self, hdu_idx: int, column: str, detector):
        values = self.column(hdu_idx, column)
        if detector is None:
            return values
        if hasattr(detector, 'number'):
            number = detector.number
        elif isinstance(detector, str):
            from .detectors import BurstCubeDetectors
            number = BurstCubeDetectors.from_str(detector.upper()).number
        else:
            number = int(detector)
        return values[:, number]
