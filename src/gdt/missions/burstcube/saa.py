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

    Two corrections are applied to what the file literally contains, both
    documented at :func:`~gdt.missions.burstcube.caldb.saa_region`:

    * its ``X``/``Y`` columns are latitude/longitude, not the
      longitude/latitude its own header comments claim;
    * its 19 vertices do not repeat the first point, so the polygon is open.
      A 20th vertex closing it is appended here, which makes
      :meth:`~gdt.core.geomagnetic.SouthAtlanticAnomaly.is_closed` true and
      means plotting ``longitude`` against ``latitude`` draws a closed
      outline rather than a broken one. It changes no containment result --
      :class:`~matplotlib.path.Path` closes an open polygon implicitly --
      only what you see and what ``is_closed`` reports.

    Args:
        cache_dir (Path, optional): The local CALDB cache directory, passed
            to :func:`~gdt.missions.burstcube.caldb.saa_region`. Defaults to
            the standard CALDB cache location.
    """

    def __init__(self, cache_dir: Optional[Union[str, Path]] = None):
        region = caldb.saa_region(cache_dir=cache_dir)
        latitude = np.asarray(region.latitude, dtype=float)
        longitude = np.asarray(region.longitude, dtype=float)
        if latitude[0] != latitude[-1] or longitude[0] != longitude[-1]:
            latitude = np.append(latitude, latitude[0])
            longitude = np.append(longitude, longitude[0])
        self._latitude = latitude
        self._longitude = longitude
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
