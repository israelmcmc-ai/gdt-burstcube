"""Offline tests for gdt.missions.burstcube.response, using synthetic DRMs
(no network, no real .rsp files needed)."""
import healpy as hp
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.time import Time as AstropyTime
from gdt.core.data_primitives import ResponseMatrix
from gdt.core.spectra.functions import PowerLaw

from gdt.missions.burstcube import caldb
from gdt.missions.burstcube.frame import BurstCubeFrame
from gdt.missions.burstcube.response import (BurstCubeResponseGrid,
                                             BurstCubeRsp)


def _synthetic_drm(num_ebins=20, num_chans=64, seed=0, scale=1.0):
    """Build a synthetic 64-channel ResponseMatrix with reproducible random
    values, scaled by `scale` (used to fabricate a "brighter"/"dimmer"
    direction for the forward/rearward comparison).
    """
    rng = np.random.default_rng(seed)
    photon_lo = np.geomspace(10.0, 5000.0, num_ebins + 1)[:-1]
    photon_hi = np.geomspace(10.0, 5000.0, num_ebins + 1)[1:]
    chan_lo = np.linspace(0.0, 1000.0, num_chans + 1)[:-1]
    chan_hi = np.linspace(0.0, 1000.0, num_chans + 1)[1:]
    matrix = scale * rng.uniform(0.0, 5.0, size=(num_ebins, num_chans))
    return ResponseMatrix(matrix, photon_lo, photon_hi, chan_lo, chan_hi)


def _synthetic_rsp(detector='CS0', **kwargs):
    return BurstCubeRsp.from_data(_synthetic_drm(**kwargs), detector=detector)


def _block_caldb_downloads(monkeypatch):
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))


def test_fold_matches_independent_geometric_mean_formula():
    """Independently (without calling fold_spectrum's own internals)
    reimplement sum_i matrix[i,j] * f(centroid_i) * width_i, with the
    centroid as the geometric mean sqrt(elo*ehi) -- section 21 confirms
    gdt-core uses exactly this, not a trapezoid or midpoint rule.
    """
    rsp = _synthetic_rsp()
    pl = PowerLaw()
    params = (0.01, -2.0)

    folded = rsp.drm.fold_spectrum(pl.fit_eval, params)

    photon_lo = rsp.drm.photon_bins.low_edges()
    photon_hi = rsp.drm.photon_bins.high_edges()
    matrix = rsp.drm.matrix
    expected = np.zeros(rsp.num_chans)
    for i in range(rsp.num_ebins):
        centroid = np.sqrt(photon_lo[i] * photon_hi[i])
        width = photon_hi[i] - photon_lo[i]
        flux = pl.fit_eval(params, np.array([centroid]))[0]
        expected += matrix[i, :] * flux * width

    np.testing.assert_allclose(folded, expected, rtol=1e-12)


def test_to_cbd_gives_16_channels_preserves_total_and_group_sums(tmp_path, monkeypatch):
    _block_caldb_downloads(monkeypatch)
    rsp = _synthetic_rsp()
    pl = PowerLaw()
    params = (0.01, -2.0)

    folded_64 = rsp.drm.fold_spectrum(pl.fit_eval, params)
    cbd_rsp = rsp.to_cbd()

    assert cbd_rsp.num_chans == 16
    folded_16 = cbd_rsp.drm.fold_spectrum(pl.fit_eval, params)
    assert folded_16.sum() == pytest.approx(folded_64.sum(), rel=1e-9)

    edge_indices = caldb.regroup_edges('CS0', 64, 16)
    expected_groups = np.add.reduceat(folded_64, edge_indices[:-1])
    np.testing.assert_allclose(folded_16, expected_groups, rtol=1e-9)


def test_to_cbd_uses_caldb_eb16_energies_not_rsp_internal(monkeypatch):
    _block_caldb_downloads(monkeypatch)
    rsp = _synthetic_rsp(detector='CS0')
    cbd_rsp = rsp.to_cbd()

    eb16 = caldb.ebounds('CS0', 16, cache_dir=None)
    np.testing.assert_allclose(cbd_rsp.ebounds.low_edges(), eb16.low_edges())
    np.testing.assert_allclose(cbd_rsp.ebounds.high_edges(), eb16.high_edges())
    # and NOT the synthetic DRM's own internal 0-1000 keV channel grid
    assert cbd_rsp.ebounds.high_edges()[-1] != pytest.approx(1000.0)


def test_slice_channels_keeps_range_untouched_and_drops_the_rest(monkeypatch):
    """slice_channels selects channels; unlike to_cbd it does not regroup
    them, so each kept channel's own column of the DRM (and its own ebounds
    edges) should be untouched, not merged or reweighted.
    """
    _block_caldb_downloads(monkeypatch)
    rsp = _synthetic_rsp(num_chans=16)

    sliced = rsp.slice_channels(2, 15)

    assert sliced.num_chans == 14
    np.testing.assert_allclose(sliced.drm.matrix, rsp.drm.matrix[:, 2:16])
    np.testing.assert_allclose(sliced.ebounds.low_edges(),
                               rsp.ebounds.low_edges()[2:16])
    np.testing.assert_allclose(sliced.ebounds.high_edges(),
                               rsp.ebounds.high_edges()[2:16])


def test_slice_channels_out_of_range_raises(monkeypatch):
    _block_caldb_downloads(monkeypatch)
    rsp = _synthetic_rsp(num_chans=16)

    with pytest.raises(ValueError):
        rsp.slice_channels(15, 2)  # start after stop

    with pytest.raises(ValueError):
        rsp.slice_channels(0, 16)  # stop out of bounds (0-indexed, 16 chans)


def test_caldb_regroup_edges_have_17_entries_for_16_channels():
    """Regression for the spec-section-21 correction: a 16-entry list
    (omitting the trailing 63) silently gives 15 channels while still
    preserving the folded total, so it would NOT raise -- only checking the
    entry count and the actual regrouped num_chans catches it.
    """
    edge_indices = caldb.regroup_edges('CS0', 64, 16)
    assert len(edge_indices) == 17

    rsp = _synthetic_rsp()
    correct = rsp.drm.rebin(edge_indices=edge_indices)
    assert correct.num_chans == 16

    # the buggy earlier draft of section 10 (spec section 21) omitted the
    # value 63 specifically -- not simply the last entry -- merging the
    # (59,63) and (63,64) groups into one (59,64) group.
    buggy_earlier_draft = np.array([i for i in edge_indices if i != 63])
    assert len(buggy_earlier_draft) == 16
    wrong = rsp.drm.rebin(edge_indices=buggy_earlier_draft)
    assert wrong.num_chans == 15  # silently wrong, not a crash -- the whole danger

    pl = PowerLaw()
    folded_full = rsp.drm.fold_spectrum(pl.fit_eval, (0.01, -2.0))
    folded_truncated = wrong.fold_spectrum(pl.fit_eval, (0.01, -2.0))
    # the total is preserved even though a channel silently went missing
    assert folded_truncated.sum() == pytest.approx(folded_full.sum(), rel=1e-9)


def test_fold_is_linear():
    rsp = _synthetic_rsp()
    pl = PowerLaw()

    f = rsp.drm.fold_spectrum(pl.fit_eval, (0.01, -2.0))
    g = rsp.drm.fold_spectrum(pl.fit_eval, (0.02, -1.5))

    scaled = rsp.drm.fold_spectrum(lambda p, e: 3.0 * pl.fit_eval(p, e), (0.01, -2.0))
    np.testing.assert_allclose(scaled, 3.0 * f, rtol=1e-12)

    def summed_model(p, e):
        return pl.fit_eval((0.01, -2.0), e) + pl.fit_eval((0.02, -1.5), e)
    summed = rsp.drm.fold_spectrum(summed_model, None)
    np.testing.assert_allclose(summed, f + g, rtol=1e-12)


def test_get_interp_weights_sum_to_one():
    for theta_deg, phi_deg in [(2.0, 90.0), (45.0, 10.0), (170.0, 300.0), (0.0, 0.0)]:
        _, weights = hp.get_interp_weights(caldb.response_grid().nside, np.radians(theta_deg),
                                           np.radians(phi_deg))
        assert weights.sum() == pytest.approx(1.0, abs=1e-10)


def test_interp_at_exact_pixel_centre_reproduces_nearest_pixel_drm():
    """interp=True at an exact pixel centre must reproduce the nearest-
    pixel (interp=False) DRM exactly, since get_interp_weights returns
    [1., 0., 0., 0.] there (verified, spec section 21).
    """
    class _FakeGrid(BurstCubeResponseGrid):
        def _load_pixel_uncached(self, det_name, pixel):
            return _synthetic_rsp(detector=det_name, seed=pixel)

    grid = _FakeGrid(detectors='CS0')
    theta0, phi0 = hp.pix2ang(caldb.response_grid().nside, 0)
    az0, zen0 = np.degrees(phi0), np.degrees(theta0)

    interp_rsp = grid.get_drm('CS0', az=az0, zen=zen0, interp=True)
    nearest_rsp = grid.get_drm('CS0', pix=0)
    np.testing.assert_allclose(interp_rsp.drm.matrix, nearest_rsp.drm.matrix)


def test_forward_hemisphere_effective_area_exceeds_rearward_by_a_margin():
    """Effective area from a forward-hemisphere direction must exceed a
    rearward one by a generous margin. Does NOT assert the maximum is
    exactly at the boresight.
    """
    forward = _synthetic_rsp(seed=1, scale=10.0)
    rearward = _synthetic_rsp(seed=2, scale=1.0)

    forward_area = forward.drm.channel_effective_area().counts.sum()
    rearward_area = rearward.drm.channel_effective_area().counts.sum()
    assert forward_area > 3.0 * rearward_area


def test_skycoord_without_attitude_raises_documented_error():
    grid = BurstCubeResponseGrid(detectors='CS0')
    coord = SkyCoord(10, 20, unit='deg')

    with pytest.raises(ValueError, match='caveat #1'):
        grid.get_drm('CS0', skycoord=coord, frame=None)

    frame_no_quat = BurstCubeFrame(obstime=AstropyTime('2024-05-30T00:00:00'))
    with pytest.raises(ValueError, match='caveat #1'):
        grid.get_drm('CS0', skycoord=coord, frame=frame_no_quat)
