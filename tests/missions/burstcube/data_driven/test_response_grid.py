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
"""Data-driven tests for gdt.missions.burstcube.response's
BurstCubeResponseGrid, against the real pixel 0/1 files for all 4 detectors.
These files are NOT split into per-detector subdirectories (they are
whatever `gdt-data download burstcube` flattens `burstcube.urls` into), so
this also exercises pixel->filename resolution correctly picking the right
detector out of a directory that holds all 4.
"""
import numpy as np
import pytest

from gdt.missions.burstcube.response import BurstCubeResponseGrid

from .conftest import BURSTCUBE_DATA_DIR, real_file


@pytest.mark.parametrize('detector', ['CS0', 'CS1', 'CS2', 'CS3'])
def test_get_drm_by_pixel_resolves_the_correct_detector(detector):
    """Each detector's own pixel-0/1 files must be found and opened as that
    detector, even though all 4 detectors' files sit in one flat directory.
    """
    real_file(f'bc{detector.lower()}_px0000t003p045_20240417v211.rsp.gz')
    real_file(f'bc{detector.lower()}_px0001t003p135_20240417v211.rsp.gz')

    grid = BurstCubeResponseGrid(detectors=detector, local_dir=str(BURSTCUBE_DATA_DIR))
    rsp0 = grid.get_drm(detector, pix=0)
    rsp1 = grid.get_drm(detector, pix=1)

    assert rsp0.detector == detector
    assert rsp1.detector == detector
    assert rsp0.pixel == 0
    assert rsp1.pixel == 1


def test_four_detectors_at_the_same_pixel_have_different_responses():
    """A sanity check that resolution truly distinguishes detectors: the
    four detectors' pixel-0 DRMs must not be identical to each other.
    """
    for det in ('CS0', 'CS1', 'CS2', 'CS3'):
        real_file(f'bc{det.lower()}_px0000t003p045_20240417v211.rsp.gz')

    grid = BurstCubeResponseGrid(local_dir=str(BURSTCUBE_DATA_DIR))
    matrices = {det: grid.get_drm(det, pix=0).drm.matrix for det in
               ('CS0', 'CS1', 'CS2', 'CS3')}
    dets = list(matrices)
    for i in range(len(dets)):
        for j in range(i + 1, len(dets)):
            assert not np.allclose(matrices[dets[i]], matrices[dets[j]])
