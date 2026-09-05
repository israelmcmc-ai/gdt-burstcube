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
"""Offline tests for gdt.missions.burstcube.saa. Uses the bundled CALDB SAA
region file (downloads blocked, forcing the bundled fallback), same as the
round-1 CALDB tests.
"""
import numpy as np
import pytest

from gdt.missions.burstcube import caldb
from gdt.missions.burstcube.saa import BurstCubeSaa


def test_polygon_matches_caldb_file(tmp_path, monkeypatch):
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))
    saa = BurstCubeSaa(cache_dir=tmp_path)
    assert saa.num_points == 19
    assert saa.longitude[0] == pytest.approx(-30.0, abs=1e-3)
    assert saa.latitude[0] == pytest.approx(33.9, abs=1e-3)


def test_contains_interior_and_exterior_points(tmp_path, monkeypatch):
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))
    saa = BurstCubeSaa(cache_dir=tmp_path)

    assert saa.contains(-20.0, -40.0) is True
    assert saa.contains(100.0, 50.0) is False


def test_contains_accepts_arrays(tmp_path, monkeypatch):
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))
    saa = BurstCubeSaa(cache_dir=tmp_path)

    lon = np.array([-20.0, 100.0])
    lat = np.array([-40.0, 50.0])
    result = saa.contains(lon, lat)
    np.testing.assert_array_equal(result, [True, False])
