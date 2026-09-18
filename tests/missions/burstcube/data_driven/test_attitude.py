"""Data-driven test against the real attitude file: spec section 17's
resolved quaternion convention, checked for all 3 real rows (not just the
synthetic fixture used offline).
"""
import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord

from gdt.missions.burstcube.attitude import BurstCubeAttitude

from .conftest import real_file

# spec section 17's verification table: row -> (RA, Dec) of the +z boresight
_EXPECTED_RA_DEC = {
    0: (48.723068, 10.861862),
    1: (350.98593, 49.458153),
    2: (182.68195, -46.037468),
}


def test_all_three_rows_reproduce_pointing_on_the_real_file():
    path = real_file('bc_csa_att.fits')
    att = BurstCubeAttitude.open(path)
    assert att.num_rows == 3

    for row, (expected_ra, expected_dec) in _EXPECTED_RA_DEC.items():
        frame = att.frame(row)
        boresight = SkyCoord(0 * u.deg, 90 * u.deg, frame=frame)
        icrs = boresight.transform_to('icrs')
        assert icrs.ra.deg == pytest.approx(expected_ra, abs=1e-3)
        assert icrs.dec.deg == pytest.approx(expected_dec, abs=1e-3)
