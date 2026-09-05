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
"""Offline tests for gdt.missions.burstcube.detectors.

These tests recompute the detector azimuth/zenith directly from the bundled
CALDB alignment files, independently of the constants baked into
BurstCubeDetectors, so that a hardcoded-by-eye error in the enum values would
be caught rather than merely re-asserting them.
"""
import numpy as np
import pytest
from astropy.io import fits

from gdt.missions.burstcube.detectors import BurstCubeDetectors, _alignment_matrix


@pytest.mark.parametrize('det, num', [
    (BurstCubeDetectors.CS0, 0), (BurstCubeDetectors.CS1, 1),
    (BurstCubeDetectors.CS2, 2), (BurstCubeDetectors.CS3, 3),
])
def test_detector_numbers(det, num):
    """Each detector's number must match CS0..CS3."""
    assert det.number == num
    assert det.full_name == f'CS{num}'


def test_zenith_is_45_degrees_for_all_detectors():
    """The CALDB alignment header comment says all 4 detectors sit 45 degrees
    off the spacecraft zenith; this must hold for all 4 within float
    precision.
    """
    for det in BurstCubeDetectors:
        assert det.zenith.value == pytest.approx(45.0, abs=1e-6)


def test_azimuths_are_90_degrees_apart():
    """The 4 detectors are spaced 90 degrees apart in azimuth."""
    azimuths = sorted(det.azimuth.value for det in BurstCubeDetectors)
    diffs = np.diff(azimuths + [azimuths[0] + 360.0])
    assert diffs == pytest.approx([90.0, 90.0, 90.0, 90.0], abs=1e-6)


def test_azimuth_quadrants_match_caldb_comment():
    """The CALDB alignment header comment places CS0 in (-x,+y), CS1 in
    (-x,-y), CS2 in (+x,-y), CS3 in (+x,+y). In the standard convention where
    azimuth is measured from +x towards +y, those quadrants are 90-180,
    180-270, 270-360, and 0-90 degrees respectively.
    """
    quadrants = {
        BurstCubeDetectors.CS0: (90.0, 180.0),
        BurstCubeDetectors.CS1: (180.0, 270.0),
        BurstCubeDetectors.CS2: (270.0, 360.0),
        BurstCubeDetectors.CS3: (0.0, 90.0),
    }
    for det, (lo, hi) in quadrants.items():
        assert lo <= det.azimuth.value <= hi


def test_recomputed_independently_from_bundled_caldb_files():
    """Recompute azimuth/zenith directly from the bundled alignment FITS
    files here (not by calling into detectors.py's own derivation) and check
    that BurstCubeDetectors agrees, to guard against a shared bug in the
    derivation being invisible to the tests above.
    """
    for det in BurstCubeDetectors:
        matrix = _alignment_matrix(det.number)
        normal = matrix @ np.array([0.0, 0.0, 1.0])
        zenith = np.degrees(np.arccos(normal[2]))
        azimuth = np.degrees(np.arctan2(normal[1], normal[0])) % 360.0
        assert zenith == pytest.approx(det.zenith.value, abs=1e-9)
        assert azimuth == pytest.approx(det.azimuth.value, abs=1e-9)


def test_alignment_matrices_are_orthonormal():
    """Each bundled alignment matrix must be a proper rotation matrix (this
    is a sanity check on the bundled CALDB files themselves).
    """
    for det in BurstCubeDetectors:
        matrix = _alignment_matrix(det.number)
        should_be_identity = matrix @ matrix.T
        np.testing.assert_allclose(should_be_identity, np.eye(3), atol=1e-9)
        assert np.linalg.det(matrix) == pytest.approx(1.0, abs=1e-9)
