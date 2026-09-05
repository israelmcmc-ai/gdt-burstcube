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
"""The BurstCube South Atlantic Anomaly (SAA) boundary polygon, from the
CALDB ``bcf/saa`` region file, modelled on
:mod:`gdt.missions.fermi.gbm.saa`.

Unlike GBM, which has a documented history of successive SAA polygons over
the mission (see ``GbmSaaPolygon1`` .. ``GbmSaaPolygon15``), BurstCube's
CALDB carries a single, unversioned SAA polygon, so there is no time-indexed
collection here.
"""
from pathlib import Path
from typing import Optional, Union

import numpy as np
from matplotlib.path import Path as MplPath

from gdt.core.geomagnetic import SouthAtlanticAnomaly

from . import caldb

__all__ = ['BurstCubeSaa']


class BurstCubeSaa(SouthAtlanticAnomaly):
    """The BurstCube SAA boundary polygon, read from the CALDB SAA region
    file (``bcf/saa/bccsa_saareg_20230101v001.fits``).

    Args:
        cache_dir (Path, optional): The local CALDB cache directory, passed
            to :func:`~gdt.missions.burstcube.caldb.saa_region`. Defaults to
            the standard CALDB cache location.
    """

    def __init__(self, cache_dir: Optional[Union[str, Path]] = None):
        region = caldb.saa_region(cache_dir=cache_dir)
        self._latitude = region.latitude
        self._longitude = region.longitude
        super().__init__()

    def contains(self, longitude, latitude):
        """Determine whether one or more (longitude, latitude) points fall
        within the SAA boundary polygon.

        Args:
            longitude (float or numpy.ndarray): East longitude(s), in degrees
            latitude (float or numpy.ndarray): Latitude(s), in degrees

        Returns:
            (bool or numpy.ndarray): A scalar bool if both inputs were
            scalar, otherwise a boolean array.
        """
        path = MplPath(np.column_stack((self.longitude, self.latitude)))
        lon = np.atleast_1d(longitude)
        lat = np.atleast_1d(latitude)
        mask = path.contains_points(np.column_stack((lon, lat)))

        if np.isscalar(longitude) and np.isscalar(latitude):
            return bool(mask[0])
        return mask
