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
"""Offline tests for gdt.missions.burstcube.orbit."""
import numpy as np
import pytest

from gdt.missions.burstcube.frame import BurstCubeFrame
from gdt.missions.burstcube.orbit import BurstCubeOrbit

from .conftest import make_orbit_fits


def test_get_spacecraft_frame_has_position_and_velocity_no_quaternion(tmp_path):
    path = tmp_path / 'orbit.fits'
    time = make_orbit_fits(path)
    orbit = BurstCubeOrbit.open(path)
    frame = orbit.get_spacecraft_frame()

    assert isinstance(frame, BurstCubeFrame)
    assert frame.quaternion is None
    assert frame.detectors is not None
    np.testing.assert_allclose(frame.obstime.burstcube, time)
    # obsgeoloc/obsgeovel are populated (non-default) for every row
    assert np.all(frame.obsgeoloc.x.value != 0.0)
    assert np.all(frame.obsgeovel.x.value != 0.0)


def test_position_converted_to_meters(tmp_path):
    """X/Y/Z are stored in km; the frame's obsgeoloc must be in (equivalent)
    SI meters, per SpacecraftFrame's own convention.
    """
    import astropy.units as u

    path = tmp_path / 'orbit.fits'
    make_orbit_fits(path)
    orbit = BurstCubeOrbit.open(path)
    frame = orbit.get_spacecraft_frame()
    assert frame.obsgeoloc.x[0].to(u.km).value == pytest.approx(7000.0)


def test_merge_combines_and_deduplicates(tmp_path):
    """merge() must combine two orbit files into one, time-sorted, with
    duplicate TIME entries removed.
    """
    path1 = tmp_path / 'orbit1.fits'
    path2 = tmp_path / 'orbit2.fits'
    make_orbit_fits(path1, n=5)
    make_orbit_fits(path2, n=5)  # deliberately identical/overlapping times

    orbit1 = BurstCubeOrbit.open(path1)
    orbit2 = BurstCubeOrbit.open(path2)
    merged = BurstCubeOrbit.merge(orbit1, orbit2)

    frame = merged.get_spacecraft_frame()
    # fully overlapping inputs must dedupe down to the original 5 rows
    assert frame.obstime.size == 5
    assert np.all(np.diff(frame.obstime.burstcube) > 0)
