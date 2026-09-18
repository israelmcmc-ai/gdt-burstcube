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


def _crossings(latitude, longitude) -> int:
    """Count the pairs of non-adjacent edges that cross, treating the
    vertices as a closed polygon. Zero means the polygon is simple.

    Args:
        latitude (np.array): The vertex latitudes, in degrees
        longitude (np.array): The vertex East longitudes, in degrees

    Returns:
        (int)
    """
    points = list(zip(longitude, latitude))
    if points[0] != points[-1]:
        points.append(points[0])
    num_edges = len(points) - 1

    def side(origin, a, b):
        return ((a[0] - origin[0]) * (b[1] - origin[1])
                - (a[1] - origin[1]) * (b[0] - origin[0]))

    def crosses(p, q, r, s):
        return (((side(r, s, p) > 0) != (side(r, s, q) > 0))
                and ((side(p, q, r) > 0) != (side(p, q, s) > 0)))

    return sum(
        crosses(points[i], points[i + 1], points[j], points[j + 1])
        for i in range(num_edges) for j in range(i + 2, num_edges)
        # the first and last edges are adjacent in a closed polygon
        if not (i == 0 and j == num_edges - 1))


class BurstCubeSaa(SouthAtlanticAnomaly):
    """The BurstCube SAA boundary polygon, read from the CALDB SAA region
    file (``bcf/saa/bccsa_saareg_20230101v001.fits``).

    Three corrections are applied to what the file literally contains. The
    first is in :func:`~gdt.missions.burstcube.caldb.saa_region`, which
    reads the ``X``/``Y`` columns as latitude/longitude rather than the
    longitude/latitude its own header comments claim. The other two are
    here:

    * **Two pairs of vertices are listed out of order**, so the boundary
      doubles back on itself and the polygon self-intersects -- visibly, as
      a spur off the eastern and western corners. Sorting the vertices by
      angle about their centroid removes both crossings and leaves the other
      15 in exactly the order the file gives them (see
      :meth:`_ordered_by_angle`).
    * **The polygon is left open**: the 19 vertices do not repeat the first
      point, so drawing them gives a broken outline. A 20th vertex closing
      it is appended, which makes
      :meth:`~gdt.core.geomagnetic.SouthAtlanticAnomaly.is_closed` true.
      Closing changes no containment result on its own --
      :class:`~matplotlib.path.Path` closes an open polygon implicitly.

    Reordering *does* change containment, by about 1% of the area: the two
    self-intersections were small bow-ties that a point-in-polygon test
    counts with the opposite sign to the rest of the region.

    Args:
        cache_dir (Path, optional): The local CALDB cache directory, passed
            to :func:`~gdt.missions.burstcube.caldb.saa_region`. Defaults to
            the standard CALDB cache location.
    """

    def __init__(self, cache_dir: Optional[Union[str, Path]] = None):
        region = caldb.saa_region(cache_dir=cache_dir)
        latitude = np.asarray(region.latitude, dtype=float)
        longitude = np.asarray(region.longitude, dtype=float)

        latitude, longitude = self._ordered_by_angle(latitude, longitude)

        if latitude[0] != latitude[-1] or longitude[0] != longitude[-1]:
            latitude = np.append(latitude, latitude[0])
            longitude = np.append(longitude, longitude[0])

        self._latitude = latitude
        self._longitude = longitude
        super().__init__()

    @staticmethod
    def _ordered_by_angle(latitude, longitude):
        """Order the vertices counter-clockwise by angle about their
        centroid, but only if that removes a self-intersection the file's own
        order has.

        The SAA is a single blob containing its own centroid, so its boundary
        vertices are in angular order around it and sorting recovers that
        order. That is not true of every polygon -- a sufficiently concave
        one has vertices no angular sweep can reach in order -- so this
        checks that the file's order really is broken first, and that
        sorting really does fix it. If either check fails the file's order is
        kept untouched, which means a future CALDB revision cannot be
        silently rearranged into something else.

        Args:
            latitude (np.array): The vertex latitudes, in degrees
            longitude (np.array): The vertex East longitudes, in degrees

        Returns:
            (np.array, np.array): The reordered latitudes and longitudes
        """
        if _crossings(latitude, longitude) == 0:
            return latitude, longitude

        angle = np.arctan2(latitude - latitude.mean(),
                           longitude - longitude.mean())
        order = np.argsort(angle)
        # keep whichever vertex the file started from as the start
        order = np.roll(order, -int(np.flatnonzero(order == 0)[0]))

        if _crossings(latitude[order], longitude[order]) > 0:
            return latitude, longitude
        return latitude[order], longitude[order]

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
