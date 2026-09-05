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
"""BurstCube detector housekeeping: the ``DETECTOR_HK1`` (1 s cadence) and
``DETECTOR_HK2`` (60 s cadence) extensions of ``auxil/bcYYMMDDcsa.hk.gz``,
each exposed as an :class:`~astropy.timeseries.TimeSeries`.

Per archive caveat #7, the per-detector energy thresholds were raised mid-mission
from ~21-30 keV to ~100 keV to suppress a non-Poissonian low-energy noise
component; :meth:`~BurstCubeDetectorHk.peak_threshold` and
:meth:`~BurstCubeDetectorHk.base_threshold` are how that change shows up in
the housekeeping data.
"""
from astropy.timeseries import TimeSeries

from gdt.core.file import FitsFileContextManager

from .headers import DetectorHkHeaders
from .time import Time

__all__ = ['BurstCubeDetectorHk']

# the 4-element housekeeping arrays are indexed by detector number 0..3,
# matching BurstCubeDetectors.number
_ENABLE_FLAG_COLUMNS = ('TRIG_ENABLED', 'DET_ENABLED', 'CBD_ENABLED',
                        'TTE_ENABLED')
_THRESHOLD_COLUMNS = ('PEAK_THRES', 'BASE_THRES')


class BurstCubeDetectorHk(FitsFileContextManager):
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
            (:class:`BurstCubeDetectorHk`)
        """
        obj = super().open(file_path, **kwargs)
        hdrs = [hdu.header for hdu in obj.hdulist]
        obj._headers = DetectorHkHeaders.from_headers(hdrs)
        return obj

    def trig_enabled(self, detector=None):
        """The ``TRIG_ENABLED`` flags from ``DETECTOR_HK1``: whether the
        onboard trigger algorithm was enabled for a detector at each time.

        Args:
            detector (int or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
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
            detector (int or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
                The detector. If omitted, returns the full (`n`, 4) array.

        Returns:
            (numpy.ndarray)
        """
        return self._enable_flag('DET_ENABLED', detector)

    def cbd_enabled(self, detector=None):
        """The ``CBD_ENABLED`` flags from ``DETECTOR_HK1``: whether CBD data
        collection was enabled for a detector at each time.

        Args:
            detector (int or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
                The detector. If omitted, returns the full (`n`, 4) array.

        Returns:
            (numpy.ndarray)
        """
        return self._enable_flag('CBD_ENABLED', detector)

    def tte_enabled(self, detector=None):
        """The ``TTE_ENABLED`` flags from ``DETECTOR_HK1``: whether TTE data
        collection was enabled for a detector at each time.

        Args:
            detector (int or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
                The detector. If omitted, returns the full (`n`, 4) array.

        Returns:
            (numpy.ndarray)
        """
        return self._enable_flag('TTE_ENABLED', detector)

    def peak_threshold(self, detector=None):
        """The ``PEAK_THRES`` energy threshold from ``DETECTOR_HK2``, in mV.
        See archive caveat #7 for the mid-mission threshold change.

        Args:
            detector (int or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
                The detector. If omitted, returns the full (`n`, 4) array.

        Returns:
            (numpy.ndarray)
        """
        return self._threshold('PEAK_THRES', detector)

    def base_threshold(self, detector=None):
        """The ``BASE_THRES`` energy threshold from ``DETECTOR_HK2``, in mV.
        See archive caveat #7 for the mid-mission threshold change.

        Args:
            detector (int or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`, optional):
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
        number = getattr(detector, 'number', detector)
        return values[:, number]
