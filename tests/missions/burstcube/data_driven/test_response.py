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
"""Data-driven tests for gdt.missions.burstcube.response, against the real
CS0 pixel-0 and pixel-1 response files.
"""
import gzip
import shutil

import numpy as np
import pytest
from gdt.core.spectra.functions import PowerLaw

from gdt.missions.burstcube.response import (EDGE_INDICES_64_TO_16,
                                             BurstCubeRsp)

from .conftest import real_file


def _open_gz(gz_path, tmp_path):
    """gzip.open works directly with astropy.io.fits, but BurstCubeRsp.open
    goes through FitsFileContextManager.open(), which expects a plain path;
    decompress to tmp_path once per test.
    """
    dest = tmp_path / gz_path.stem
    with gzip.open(gz_path, 'rb') as src, open(dest, 'wb') as dst:
        shutil.copyfileobj(src, dst)
    return dest


def test_pixel0_matrix_matches_spec_regression_values(tmp_path):
    """Spec section 10: the pixel-0 CS0 matrix has max 18.918, sum 7762.86."""
    path = _open_gz(real_file('bccs0_px0000t003p045_20240417v211.rsp.gz'), tmp_path)
    rsp = BurstCubeRsp.open(path)

    assert rsp.detector == 'CS0'
    assert rsp.pixel == 0
    assert rsp.num_chans == 64
    assert rsp.drm.matrix.max() == pytest.approx(18.918, abs=1e-3)
    assert rsp.drm.matrix.sum() == pytest.approx(7762.86, abs=1e-2)


def test_folded_power_law_matches_spec_regression_value(tmp_path):
    """Spec section 21: folding PowerLaw(0.01, -2.0) through the pixel-0 CS0
    response gives a native 64-channel total of 128.6532926116.
    """
    path = _open_gz(real_file('bccs0_px0000t003p045_20240417v211.rsp.gz'), tmp_path)
    rsp = BurstCubeRsp.open(path)

    pl = PowerLaw()
    folded = rsp.fold_spectrum(pl.fit_eval, (0.01, -2.0))
    # rel=1e-9 rather than 1e-6: the photon edges are widened to float64 at
    # read time, so this reproduces the independently computed reference to
    # full double precision. A looser tolerance would hide a regression back
    # to float32 centroids, which shifts this by ~2.5e-9.
    assert folded.counts.sum() == pytest.approx(128.6532926116, rel=1e-9)


def test_to_cbd_gives_16_channels_and_preserves_folded_total(tmp_path):
    """to_cbd() must give exactly 16 channels, preserve the folded total,
    and its per-group sums must equal np.add.reduceat of the native 64-
    channel folded spectrum at the edge indices.
    """
    path = _open_gz(real_file('bccs0_px0000t003p045_20240417v211.rsp.gz'), tmp_path)
    rsp = BurstCubeRsp.open(path)
    pl = PowerLaw()

    folded_64 = rsp.fold_spectrum(pl.fit_eval, (0.01, -2.0)).counts

    cbd_rsp = rsp.to_cbd()
    assert cbd_rsp.num_chans == 16
    folded_16 = cbd_rsp.fold_spectrum(pl.fit_eval, (0.01, -2.0)).counts

    # regrouping the DRM first and then folding sums in a different order
    # than folding at 64 channels and summing the result, so the two totals
    # agree only to the ~1e-7 relative precision of the underlying float32
    # MATRIX column, not bit-for-bit -- this is expected floating-point
    # summation-order noise, not a regression.
    assert folded_16.sum() == pytest.approx(folded_64.sum(), rel=1e-6)
    expected_groups = np.add.reduceat(folded_64, EDGE_INDICES_64_TO_16[:-1])
    np.testing.assert_allclose(folded_16, expected_groups, rtol=1e-6)


def test_to_cbd_uses_caldb_eb16_not_rsp_internal_ebounds(tmp_path):
    """The regrouped channel energies must be CALDB's eb16 for CS0 (top
    channel 1298.79-1648.51 keV), not the .rsp's own internal (rounded,
    detector-independent) grid (top channel 1290-5000 keV).
    """
    path = _open_gz(real_file('bccs0_px0000t003p045_20240417v211.rsp.gz'), tmp_path)
    rsp = BurstCubeRsp.open(path)
    cbd_rsp = rsp.to_cbd()

    assert cbd_rsp.ebounds.low_edges()[1] == pytest.approx(28.48, abs=1e-2)
    assert cbd_rsp.ebounds.high_edges()[-1] == pytest.approx(1648.51, abs=1e-2)
    # and NOT the .rsp's own internal grid's top edge
    assert cbd_rsp.ebounds.high_edges()[-1] != pytest.approx(5000.0, abs=1.0)


def test_fold_is_linear(tmp_path):
    path = _open_gz(real_file('bccs0_px0000t003p045_20240417v211.rsp.gz'), tmp_path)
    rsp = BurstCubeRsp.open(path)
    pl = PowerLaw()

    f = rsp.fold_spectrum(pl.fit_eval, (0.01, -2.0)).counts
    g = rsp.fold_spectrum(pl.fit_eval, (0.02, -1.5)).counts

    scaled = rsp.drm.fold_spectrum(
        lambda params, e: 3.0 * pl.fit_eval(params, e), (0.01, -2.0))
    np.testing.assert_allclose(scaled, 3.0 * f, rtol=1e-9)

    def sum_model(params, e):
        return pl.fit_eval((0.01, -2.0), e) + pl.fit_eval((0.02, -1.5), e)
    summed = rsp.drm.fold_spectrum(sum_model, None)
    np.testing.assert_allclose(summed, f + g, rtol=1e-9)


def test_pixel0_and_pixel1_are_healpix_ring_neighbors(tmp_path):
    """Spec section 19: pixels 0 and 1 are HEALPix RING neighbours in the
    first ring (theta 2.924 deg, phi 45 and 135 deg) -- the two real pixel
    files exercise the 4-pixel interpolation for a direction between them.
    This is the forward-hemisphere pair the shipped data-driven files
    support; the forward-vs-rearward effective-area margin (which needs a
    rearward pixel not in the verified 17-file set) is covered by the
    offline synthetic-DRM test instead.
    """
    import healpy as hp

    theta0, phi0 = hp.pix2ang(16, 0)
    theta1, phi1 = hp.pix2ang(16, 1)
    assert np.degrees(theta0) == pytest.approx(2.924, abs=1e-2)
    assert np.degrees(theta1) == pytest.approx(2.924, abs=1e-2)
    assert np.degrees(phi0) == pytest.approx(45.0, abs=1e-2)
    assert np.degrees(phi1) == pytest.approx(135.0, abs=1e-2)

    path0 = _open_gz(real_file('bccs0_px0000t003p045_20240417v211.rsp.gz'), tmp_path)
    path1 = _open_gz(real_file('bccs0_px0001t003p135_20240417v211.rsp.gz'), tmp_path)
    rsp0 = BurstCubeRsp.open(path0)
    rsp1 = BurstCubeRsp.open(path1)
    assert rsp0.pixel == 0
    assert rsp1.pixel == 1
    # distinct (not identical) DRMs, since they're different sky directions
    assert not np.allclose(rsp0.drm.matrix, rsp1.drm.matrix)


def test_photon_edges_widened_to_float64(tmp_path):
    """The FITS energy columns are big-endian float32. ResponseMatrix keeps
    whatever dtype it is given, and photon_bin_centroids is sqrt(emin*emax),
    so leaving them narrow evaluates the geometric mean in float32 and loses
    ~2e-4 keV near 5000 keV. Widening at read time is free and removes it.
    """
    path = _open_gz(real_file('bccs0_px0000t003p045_20240417v211.rsp.gz'), tmp_path)
    drm = BurstCubeRsp.open(path).drm
    assert drm._emin.dtype == np.float64
    assert drm._emax.dtype == np.float64

    lo = np.asarray(drm.photon_bins.low_edges(), dtype=np.float64)
    hi = np.asarray(drm.photon_bins.high_edges(), dtype=np.float64)
    centroids = np.asarray(drm.photon_bin_centroids, dtype=np.float64)
    assert np.allclose(centroids, np.sqrt(lo * hi), rtol=0, atol=1e-12)
