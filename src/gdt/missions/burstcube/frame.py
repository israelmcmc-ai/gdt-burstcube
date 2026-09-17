"""The BurstCube spacecraft frame.

Attitude (spacecraft orientation) is reconstructed for only 3 brief epochs
across the entire BurstCube archive (see the ``trend/attitude/bc_csa_att.fits``
caveats). For every other time, a user who wants to convert a sky position
into BurstCube azimuth/zenith, or vice versa, must supply their own attitude
quaternion. :meth:`BurstCubeFrame.from_quaternion` is the documented path for
doing so.

Quaternion component order
---------------------------
``trend/attitude/bc_csa_att.fits`` stores each attitude quaternion as a
4-element ``QPARAM`` column, and separately gives the resulting spacecraft
boresight pointing (RA, Dec, roll) in a ``POINTING`` column, which lets the
component order be checked empirically rather than assumed. Its first row is
``QPARAM = [0.0844, 0.6314, 0.6416, 0.4273]`` and
``POINTING = [48.723068, 10.861862, 333.95041]``. Treating ``QPARAM`` as
**scalar-last** (i.e. passing it to :class:`~gdt.core.coords.Quaternion` with
the default ``scalar_first=False``, so ``x=0.0844, y=0.6314, z=0.6416,
w=0.4273``) and rotating the BurstCube frame's boresight (``az=0, zenith=0``)
into ICRS reproduces ``POINTING``'s RA/Dec to 6 decimal places. Treating it as
scalar-first instead reproduces neither the RA nor the Dec. **``QPARAM`` is
therefore scalar-last**, matching the :class:`~gdt.core.coords.Quaternion`
default, so callers passing an archive ``QPARAM`` row to
:meth:`BurstCubeFrame.from_quaternion` should leave ``scalar_first=False``.
This is only established for the frame transform itself; it says nothing
about a reader for the attitude file, which is added in a later version of
this plugin.
"""
from typing import Optional, Union

import astropy.coordinates.representation as r
from astropy.coordinates import FunctionTransform, ICRS, frame_transform_graph
from astropy.time import Time

from gdt.core.coords import Quaternion, SpacecraftFrame
from gdt.core.coords.spacecraft.frame import icrs_to_spacecraft, spacecraft_to_icrs

from .detectors import BurstCubeDetectors

__all__ = ['BurstCubeFrame', 'burstcube_to_icrs', 'icrs_to_burstcube']


class BurstCubeFrame(SpacecraftFrame):
    """The BurstCube spacecraft frame, in azimuth and zenith. The frame is
    defined by a quaternion that represents a rotation from the BurstCube
    frame to the ICRS frame. This class is a wholesale inheritance of
    :class:`~gdt.core.coords.SpacecraftFrame`, with a convenience constructor,
    :meth:`from_quaternion`, for the common case of a user-supplied attitude
    quaternion.

    Example use:

        >>> from astropy.time import Time
        >>> from gdt.missions.burstcube.frame import BurstCubeFrame
        >>> quat = [0.0844, 0.6314, 0.6416, 0.4273]  # scalar-last (x, y, z, w)
        >>> frame = BurstCubeFrame.from_quaternion(quat, Time('2024-06-29T16:53:28'))
    """
    @classmethod
    def from_quaternion(cls, quaternion: Union[Quaternion, list, tuple],
                        obstime: Time, *, scalar_first: bool = False,
                        obsgeoloc: Optional[r.CartesianRepresentation] = None,
                        obsgeovel: Optional[r.CartesianRepresentation] = None
                        ) -> 'BurstCubeFrame':
        """Build a BurstCubeFrame from a user-supplied attitude quaternion.

        This is the primary way to obtain a BurstCube frame with a working
        ICRS transform for the vast majority of the mission, for which no
        reconstructed attitude exists in the archive. The quaternion may come
        from, e.g., a scientist's own attitude reconstruction, or one of the
        3 rows of ``trend/attitude/bc_csa_att.fits`` (read by the standalone
        attitude reader added in a later version of this plugin).

        Args:
            quaternion (:class:`~gdt.core.coords.Quaternion`, list, or tuple):
                A 4-element attitude quaternion (or an (`n`, 4) array of
                quaternions), representing the rotation from the BurstCube
                frame to the ICRS frame. If not already a
                :class:`~gdt.core.coords.Quaternion`, it is converted using
                `scalar_first`.
            obstime (astropy.time.Time): The time(s) at which the quaternion
                applies.
            scalar_first (bool, optional): Set to True if `quaternion` is in
                scalar-first order, False if scalar-last. Ignored if
                `quaternion` is already a :class:`~gdt.core.coords.Quaternion`.
                Default is False.
            obsgeoloc (astropy.coordinates.representation.CartesianRepresentation, optional):
                The spacecraft position in Earth-centered Inertial
                coordinates, e.g. interpolated from an orbit file. If not
                given, defaults to the origin, and Earth-occultation
                properties (`sun_visible`, `location_visible`, etc.) will not
                be meaningful.
            obsgeovel (astropy.coordinates.representation.CartesianRepresentation, optional):
                The spacecraft orbital velocity in Earth-centered Inertial
                coordinates.

        Returns:
            (:class:`BurstCubeFrame`)
        """
        if not isinstance(quaternion, Quaternion):
            quaternion = Quaternion(quaternion, scalar_first=scalar_first)

        kwargs = {'quaternion': quaternion, 'obstime': obstime,
                 'detectors': BurstCubeDetectors}
        if obsgeoloc is not None:
            kwargs['obsgeoloc'] = obsgeoloc
        if obsgeovel is not None:
            kwargs['obsgeovel'] = obsgeovel
        return cls(**kwargs)


@frame_transform_graph.transform(FunctionTransform, BurstCubeFrame, ICRS)
def burstcube_to_icrs(burstcube_frame, icrs_frame):
    """Convert from the BurstCube frame to the ICRS frame.

    Args:
        burstcube_frame (:class:`BurstCubeFrame`): The BurstCube frame
        icrs_frame (astropy.coordinates.ICRS)

    Returns:
        (astropy.coordinates.ICRS)
    """
    return spacecraft_to_icrs(burstcube_frame, icrs_frame)


@frame_transform_graph.transform(FunctionTransform, ICRS, BurstCubeFrame)
def icrs_to_burstcube(icrs_frame, burstcube_frame):
    """Convert from the ICRS frame to the BurstCube frame.

    Args:
        icrs_frame (astropy.coordinates.ICRS)
        burstcube_frame (:class:`BurstCubeFrame`): The BurstCube frame

    Returns:
        (:class:`BurstCubeFrame`)
    """
    return icrs_to_spacecraft(icrs_frame, burstcube_frame)
