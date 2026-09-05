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
"""Offline tests for gdt.missions.burstcube.hk."""
import numpy as np
from astropy.timeseries import TimeSeries

from gdt.missions.burstcube.detectors import BurstCubeDetectors
from gdt.missions.burstcube.hk import BurstCubeDetectorHk

from .conftest import make_hk_fits


def test_hk1_and_hk2_are_timeseries(tmp_path):
    path = tmp_path / 'hk.fits'
    make_hk_fits(path)
    hk = BurstCubeDetectorHk.open(path)
    assert isinstance(hk.hk1, TimeSeries)
    assert isinstance(hk.hk2, TimeSeries)
    assert len(hk.hk1) == 10
    assert len(hk.hk2) == 3


def test_enable_flags_full_array_and_per_detector(tmp_path):
    path = tmp_path / 'hk.fits'
    make_hk_fits(path)
    hk = BurstCubeDetectorHk.open(path)

    full = hk.trig_enabled()
    assert full.shape == (10, 4)
    np.testing.assert_array_equal(hk.trig_enabled(1), full[:, 1])
    np.testing.assert_array_equal(hk.trig_enabled(0), full[:, 0])
    assert np.all(hk.trig_enabled(1) == 1)
    assert np.all(hk.trig_enabled(0) == 0)


def test_enable_flags_accept_detector_enum(tmp_path):
    path = tmp_path / 'hk.fits'
    make_hk_fits(path)
    hk = BurstCubeDetectorHk.open(path)
    np.testing.assert_array_equal(hk.det_enabled(BurstCubeDetectors.CS2),
                                  hk.det_enabled(2))


def test_enable_flags_and_thresholds_accept_detector_name_string(tmp_path):
    """A plain detector name string (e.g. 'CS0', as a user would naturally
    type) must work the same as the enum member or the bare int -- not just
    be documented as accepted.
    """
    path = tmp_path / 'hk.fits'
    make_hk_fits(path)
    hk = BurstCubeDetectorHk.open(path)
    np.testing.assert_array_equal(hk.det_enabled('CS2'), hk.det_enabled(2))
    np.testing.assert_array_equal(hk.peak_threshold('cs1'), hk.peak_threshold(1))


def test_thresholds_from_hk2(tmp_path):
    path = tmp_path / 'hk.fits'
    make_hk_fits(path)
    hk = BurstCubeDetectorHk.open(path)

    peak = hk.peak_threshold()
    assert peak.shape == (3, 4)
    np.testing.assert_allclose(peak[0], [26.93, 30.33, 21.40, 20.30], atol=1e-2)
    np.testing.assert_allclose(hk.peak_threshold(0), peak[:, 0])

    base = hk.base_threshold()
    assert base.shape == (3, 4)
