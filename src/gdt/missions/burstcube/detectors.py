"""The BurstCube detector definitions.

The azimuth and zenith of each detector normal (in spacecraft coordinates) are
derived numerically below from the bundled CALDB alignment files
(``bcf/align/bccsN_align_20210101v000.fits``), rather than hardcoded by eye.
Each alignment file is header-only and holds a 3x3 DET->SAT rotation matrix in
the keywords ``ALIGNM11`` .. ``ALIGNM33``. The CALDB header COMMENT block
(quoted here for reference) says:

    * the spacecraft z-axis is the instrument boresight
    * the spacecraft y-axis points toward the solar panels
    * the x-axis is right-handed
    * the 4 detectors sit in a square pyramid, each 45 deg from vertical, fixed
    * CS0 is in (-x,+y), CS1 in (-x,-y), CS2 in (+x,-y), CS3 in (+x,+y)

The detector normal in its own frame is +z = (0, 0, 1). Rotating that vector
into spacecraft coordinates with the alignment matrix, and then converting to
spherical coordinates (zenith measured from the spacecraft +z axis, azimuth
measured from +x towards +y), gives:

    CS0: zenith = 45.0 deg, azimuth = 135.0 deg   (quadrant -x,+y)
    CS1: zenith = 45.0 deg, azimuth = 225.0 deg   (quadrant -x,-y)
    CS2: zenith = 45.0 deg, azimuth = 315.0 deg   (quadrant +x,-y)
    CS3: zenith = 45.0 deg, azimuth =  45.0 deg   (quadrant +x,+y)

confirming both properties promised by the CALDB comment block: the zenith
angle is exactly 45 degrees for all four detectors, and the azimuths are
spaced exactly 90 degrees apart, in the quadrants the comment names. This
derivation is repeated verbatim by
``tests/missions/burstcube/offline/test_detectors.py``, which recomputes it
directly from the bundled files rather than trusting the constants below.
"""
from importlib.resources import files

import numpy as np
from astropy.io import fits
import astropy.units as u

from gdt.core.detector import Detectors

__all__ = ['BurstCubeDetectors']

# Detector normal in its own frame: +z (the boresight direction).
_DET_BORESIGHT = np.array([0.0, 0.0, 1.0])


def _alignment_matrix(det_num: int) -> np.ndarray:
    """Read the bundled CALDB alignment file for one detector and return its
    3x3 DET->SAT rotation matrix.

    Args:
        det_num (int): The detector number, 0-3

    Returns:
        (numpy.ndarray): A (3, 3) rotation matrix
    """
    filename = f'bccs{det_num}_align_20210101v000.fits'
    path = files('gdt.missions.burstcube.data').joinpath(filename)
    with fits.open(path) as hdulist:
        header = hdulist[0].header
        matrix = np.array([
            [header['ALIGNM11'], header['ALIGNM12'], header['ALIGNM13']],
            [header['ALIGNM21'], header['ALIGNM22'], header['ALIGNM23']],
            [header['ALIGNM31'], header['ALIGNM32'], header['ALIGNM33']],
        ])
    return matrix


def _pointing_from_alignment(det_num: int):
    """Derive the (azimuth, zenith) of a detector's normal in spacecraft
    coordinates from its bundled CALDB alignment matrix.

    Args:
        det_num (int): The detector number, 0-3

    Returns:
        (float, float): The azimuth and zenith, in degrees
    """
    matrix = _alignment_matrix(det_num)
    normal_sat = matrix @ _DET_BORESIGHT
    zenith = np.degrees(np.arccos(normal_sat[2]))
    azimuth = np.degrees(np.arctan2(normal_sat[1], normal_sat[0])) % 360.0
    return azimuth, zenith


_CS0_AZ, _CS0_ZEN = _pointing_from_alignment(0)
_CS1_AZ, _CS1_ZEN = _pointing_from_alignment(1)
_CS2_AZ, _CS2_ZEN = _pointing_from_alignment(2)
_CS3_AZ, _CS3_ZEN = _pointing_from_alignment(3)


# unfortunately Sphinx has a major bug that prevents the autodoc of Enums,
# so we have to define all of this in the docstring...

class BurstCubeDetectors(Detectors):
    """The BurstCube Detector name and orientation definitions.

    BurstCube has 4 CsI scintillator detectors, CS0 through CS3, arranged in a
    square pyramid, each 45 degrees from the spacecraft zenith and spaced 90
    degrees apart in azimuth. There is no 5th "CSA" detector: CSA in file and
    header names means "all detectors / spacecraft-level", not a detector.

    .. rubric:: Attributes Summary
    .. autosummary::

      azimuth
      elevation
      full_name
      number
      zenith

    .. rubric:: Methods Summary

    .. autosummary::

      from_full_name
      from_num
      from_str
      pointing
      skycoord

    .. rubric:: Attributes Documentation

    .. autoattribute:: azimuth
    .. autoattribute:: elevation
    .. autoattribute:: full_name
    .. autoattribute:: number
    .. autoattribute:: zenith

    .. rubric:: Methods Documentation

    .. automethod:: from_full_name
    .. automethod:: from_num
    .. automethod:: from_str
    .. automethod:: pointing
    .. automethod:: skycoord
    """
    CS0 = ('CS0', 0, _CS0_AZ * u.deg, _CS0_ZEN * u.deg)
    CS1 = ('CS1', 1, _CS1_AZ * u.deg, _CS1_ZEN * u.deg)
    CS2 = ('CS2', 2, _CS2_AZ * u.deg, _CS2_ZEN * u.deg)
    CS3 = ('CS3', 3, _CS3_AZ * u.deg, _CS3_ZEN * u.deg)
