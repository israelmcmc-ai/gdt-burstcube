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
"""Offline tests for gdt.missions.burstcube.gti."""
import numpy as np
import pytest
from gdt.core.data_primitives import Gti

from gdt.missions.burstcube.cbd import BurstCubeCbd
from gdt.missions.burstcube.gti import BurstCubeGti, apply_to, complement, intersect, union

from .conftest import make_cbd_fits, make_trend_gti_fits


def test_open_reads_gti_from_file(tmp_path):
    path = tmp_path / 'trend.gti'
    starts, stops = make_trend_gti_fits(path)
    reader = BurstCubeGti.open(path)
    np.testing.assert_allclose(reader.gti.low_edges(), starts)
    np.testing.assert_allclose(reader.gti.high_edges(), stops)


def test_intersect_and_union():
    g1 = Gti.from_list([(0, 10), (20, 30)])
    g2 = Gti.from_list([(5, 25)])
    assert intersect(g1, g2).as_list() == [(5.0, 10.0), (20.0, 25.0)]
    assert union(g1, g2).as_list() == [(0.0, 30.0)]


def test_complement_finds_the_gaps():
    g = Gti.from_list([(0, 10), (20, 30)])
    assert complement(g, 0, 30).as_list() == [(10.0, 20.0)]


def test_complement_of_full_coverage_is_none():
    g = Gti.from_list([(0, 30)])
    assert complement(g, 0, 30) is None


def test_complement_includes_leading_and_trailing_gaps():
    g = Gti.from_list([(10, 20)])
    assert complement(g, 0, 30).as_list() == [(0.0, 10.0), (20.0, 30.0)]


def test_apply_to_reproduces_cl_style_cut_from_uf(tmp_path):
    """apply_to() must cut a Phaii down to a GTI, e.g. reproducing a _cl
    file's cut given its _uf counterpart and the matching trend GTI.
    """
    n = 20
    time = 107629263.33 + (np.arange(n) + 1) * 0.256
    path = tmp_path / 'cbd_uf.fits'
    make_cbd_fits(path, time=time, gti=[(time[0] - 0.256, time[-1])])
    cbd = BurstCubeCbd.open(path)

    # a trend GTI that only covers the first half of the data
    cutoff = time[9]
    trimmed_gti = Gti.from_list([(time[0] - 0.256, cutoff)])
    trimmed = apply_to(trimmed_gti, cbd)

    assert trimmed.time_range[1] == pytest.approx(cutoff)
    assert trimmed.data.size[0] == 10
