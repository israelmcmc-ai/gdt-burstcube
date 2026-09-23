"""Offline tests for gdt.missions.burstcube.frame."""
import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from astropy.time import Time as AstropyTime

from gdt.missions.burstcube.detectors import BurstCubeDetectors
from gdt.missions.burstcube.frame import BurstCubeFrame

# Row 1 of trend/attitude/bc_csa_att.fits: QPARAM and its own POINTING (RA,
# Dec, roll), used to check the quaternion component-order convention.
_ROW1_QPARAM = [0.0844, 0.6314, 0.6416, 0.4273]
_ROW1_RA_DEC = (48.723068, 10.861862)
_ROW1_OBSTIME = AstropyTime('2024-06-29T16:53:28')


def test_from_quaternion_builds_a_frame_with_detectors_attached():
    """from_quaternion must attach BurstCubeDetectors so that
    detector_angle() and skycoord() work without further setup.
    """
    frame = BurstCubeFrame.from_quaternion(_ROW1_QPARAM, _ROW1_OBSTIME)
    assert frame.detectors is BurstCubeDetectors
    assert frame.quaternion is not None


def test_scalar_last_reproduces_archive_pointing():
    """QPARAM must be treated as scalar-last (the Quaternion default) to
    reproduce the attitude file's own POINTING column: rotating the
    BurstCube frame's boresight (az=0, zenith=0) into ICRS must reproduce
    POINTING's RA/Dec for row 1, to 6 decimal places.
    """
    frame = BurstCubeFrame.from_quaternion(_ROW1_QPARAM, _ROW1_OBSTIME)
    boresight = SkyCoord(0 * u.deg, 90 * u.deg, frame=frame)
    icrs = boresight.transform_to('icrs')
    assert icrs.ra.deg[0] == pytest.approx(_ROW1_RA_DEC[0], abs=1e-5)
    assert icrs.dec.deg[0] == pytest.approx(_ROW1_RA_DEC[1], abs=1e-5)


def test_scalar_first_does_not_reproduce_archive_pointing():
    """The opposite convention (scalar-first) must NOT reproduce POINTING,
    confirming that scalar-last is not merely one of two conventions that
    both happen to work.
    """
    frame = BurstCubeFrame.from_quaternion(_ROW1_QPARAM, _ROW1_OBSTIME,
                                           scalar_first=True)
    boresight = SkyCoord(0 * u.deg, 90 * u.deg, frame=frame)
    icrs = boresight.transform_to('icrs')
    assert icrs.ra.deg[0] != pytest.approx(_ROW1_RA_DEC[0], abs=1.0)


def test_round_trip_icrs_to_frame_and_back():
    """A sky coordinate transformed into the BurstCube frame and back to
    ICRS must reproduce the original RA/Dec.
    """
    frame = BurstCubeFrame.from_quaternion(_ROW1_QPARAM, _ROW1_OBSTIME)
    original = SkyCoord(10.0 * u.deg, 20.0 * u.deg)
    az_zen = original.transform_to(frame)
    back = az_zen.transform_to('icrs')
    assert back.ra.deg == pytest.approx(original.ra.deg, abs=1e-6)
    assert back.dec.deg == pytest.approx(original.dec.deg, abs=1e-6)


def test_from_quaternion_accepts_orbit_position():
    """from_quaternion must accept an optional spacecraft position, for
    combining a user quaternion with an interpolated orbit position.
    """
    import astropy.coordinates.representation as r

    position = r.CartesianRepresentation([7000.0, 0.0, 0.0], unit=u.km)
    frame = BurstCubeFrame.from_quaternion(_ROW1_QPARAM, _ROW1_OBSTIME,
                                           obsgeoloc=position)
    assert frame.obsgeoloc.x.to(u.km).value == pytest.approx(7000.0)
