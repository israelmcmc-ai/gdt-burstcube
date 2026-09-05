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
"""Offline tests for gdt.missions.burstcube.catalog. Uses `cached=True` to
read a synthetic catalog FITS table without any network access, exactly
as gdt.core.heasarc.BrowseCatalog supports.
"""
import numpy as np

from gdt.missions.burstcube.catalog import BurstCubeObsCatalog

from .conftest import make_catalog_fits


def _make_catalog(tmp_path):
    make_catalog_fits(
        tmp_path / 'burcbmastr.fit',
        obsid=['240530', '240531', '240602'],
        num_events=[4, 0, 4],
        exposure_by_detector={
            'cs0': [3219.4, np.nan, 1000.0],
            'cs1': [3219.4, np.nan, 1000.0],
            'cs2': [3219.4, np.nan, 1000.0],
            'cs3': [3219.4, np.nan, 1000.0],
        })
    return BurstCubeObsCatalog(cache_path=str(tmp_path), cached=True)


def test_tte_days_selects_num_events_greater_than_zero(tmp_path):
    cat = _make_catalog(tmp_path)
    tte = cat.tte_days()
    assert tte.num_rows == 2
    assert set(tte.get_table(columns=['obsid'])['OBSID']) == {'240530', '240602'}


def test_exposure_preserves_missing_as_nan_not_zero(tmp_path):
    cat = _make_catalog(tmp_path)
    exposure = cat.exposure('CS0')
    assert exposure[0] == 3219.4
    assert np.isnan(exposure[1])
    assert exposure[2] == 1000.0


def test_exposure_accepts_detector_enum(tmp_path):
    from gdt.missions.burstcube.detectors import BurstCubeDetectors
    cat = _make_catalog(tmp_path)
    np.testing.assert_allclose(
        cat.exposure(BurstCubeDetectors.CS0)[[0, 2]], [3219.4, 1000.0])
