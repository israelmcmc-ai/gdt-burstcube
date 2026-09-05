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
"""Offline tests for gdt.missions.burstcube.attitude."""
import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord

from gdt.missions.burstcube.attitude import BurstCubeAttitude

from .conftest import make_attitude_fits


def test_exactly_three_rows(tmp_path):
    path = tmp_path / 'att.fits'
    make_attitude_fits(path)
    att = BurstCubeAttitude.open(path)
    assert att.num_rows == 3
    assert att.time.size == 3
    assert att.quaternion.scalar_last.shape == (3, 4)
    assert att.pointing.shape == (3, 3)


@pytest.mark.parametrize('row', [0, 1, 2])
def test_frame_boresight_reproduces_pointing(tmp_path, row):
    """Each row's own POINTING RA/Dec must be reproduced by rotating that
    row's frame boresight (az=0, zenith=0) into ICRS, per the resolved
    quaternion convention (spec section 17): agreement to <1e-3 deg.
    """
    path = tmp_path / 'att.fits'
    _, _, pointing = make_attitude_fits(path)
    att = BurstCubeAttitude.open(path)

    frame = att.frame(row)
    boresight = SkyCoord(0 * u.deg, 90 * u.deg, frame=frame)
    icrs = boresight.transform_to('icrs')

    assert icrs.ra.deg == pytest.approx(pointing[row, 0], abs=1e-3)
    assert icrs.dec.deg == pytest.approx(pointing[row, 1], abs=1e-3)


def test_no_interpolation_api_offered(tmp_path):
    """This is a standalone reader with no interpolation: there must be no
    method that accepts an arbitrary time and returns an attitude. The only
    per-epoch accessor is `frame(row)`, indexed by row, not by time.
    """
    path = tmp_path / 'att.fits'
    make_attitude_fits(path)
    att = BurstCubeAttitude.open(path)

    for name in ('at', 'interpolate', 'nearest'):
        assert not hasattr(att, name)
