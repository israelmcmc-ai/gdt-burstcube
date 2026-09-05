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
"""The BurstCube observation catalog, ``burcbmastr``, via HEASARC Browse.

Modelled on ``gdt.missions.fermi.gbm.catalogs``: subclassing
:class:`~gdt.core.heasarc.BrowseCatalog` is a ~5-line exercise once the table
name is known.

Per the catalog schema, ``num_events`` is 4 on exactly the days with TTE data
and 0 elsewhere, so ``num_events > 0`` is the supported way to find TTE days
(:meth:`~BurstCubeObsCatalog.tte_days`). Many rows have empty
exposure/integration-time/rate fields for detectors with no CBD that day;
those are represented as NaN in the underlying table and must be treated as
missing, not as 0.0 -- the per-detector accessors here never coerce a
missing value to zero.
"""
import os

import numpy as np

from gdt.core import cache_path
from gdt.core.heasarc import BrowseCatalog

__all__ = ['BurstCubeObsCatalog']

burstcube_cache_path = os.path.join(cache_path, 'burstcube')


class BurstCubeObsCatalog(BrowseCatalog):
    """Interfaces with the BurstCube observation catalog (``burcbmastr``)
    via HEASARC Browse.

    Parameters:
        cache_path (str): The path where the cached catalog will live.
        cached (bool, optional): Set to True to read from the cached file
                                 instead of querying HEASARC. Default is False.
        verbose (bool, optional): Default is True
    """

    def __init__(self, cache_path=burstcube_cache_path, **kwargs):
        super().__init__(cache_path, table='burcbmastr', **kwargs)

    def tte_days(self):
        """Return the subset of the catalog with TTE data, i.e. the
        observation days with ``num_events > 0``.

        Returns:
            (:class:`BurstCubeObsCatalog`)
        """
        return self.slice('num_events', lo=1)

    def exposure(self, detector):
        """The per-observation total exposure for one detector, in seconds.
        Rows with no CBD data for that detector (an empty field in the
        catalog) are NaN, not 0.0.

        Args:
            detector (str or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`):
                The detector, e.g. ``'CS0'``.

        Returns:
            (numpy.ndarray)
        """
        return self._missing_as_nan(self._detector_column('exposure', detector))

    def integration_time(self, detector):
        """The per-observation CBD time resolution for one detector, in
        seconds (nominally 0.256). Rows with no CBD data for that detector
        are NaN, not 0.0.

        Args:
            detector (str or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`):
                The detector, e.g. ``'CS0'``.

        Returns:
            (numpy.ndarray)
        """
        return self._missing_as_nan(
            self._detector_column('integration_time', detector))

    def count_rate(self, detector):
        """The per-observation average count rate for one detector, in
        counts/s. Rows with no CBD data for that detector are NaN, not 0.0.

        Args:
            detector (str or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`):
                The detector, e.g. ``'CS0'``.

        Returns:
            (numpy.ndarray)
        """
        return self._missing_as_nan(self._detector_column('count_rate', detector))

    def _detector_column(self, prefix, detector):
        name = getattr(detector, 'name', detector).lower()
        return f'{prefix}_{name}'

    def _missing_as_nan(self, column):
        """Return a catalog column as float64, preserving (never
        introducing) NaN for missing entries, rather than letting them read
        as 0.0.

        Args:
            column (str): The column name

        Returns:
            (numpy.ndarray)
        """
        data = np.ma.asarray(self._data[column])
        if np.ma.is_masked(data) or data.mask is not np.ma.nomask:
            return data.filled(np.nan).astype(float)
        return data.astype(float)
