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
"""Offline tests for gdt.missions.burstcube.cbd."""
import warnings

import numpy as np
import pytest
from gdt.core.binning.binned import combine_by_factor

from gdt.missions.burstcube.cbd import BurstCubeCBD
from gdt.missions.burstcube.headers import CbdHeaders, CbdUnfilteredHeaders

from .conftest import TIMEDEL_CBD, make_cbd_fits


def test_timepixr_bin_edges(tmp_path):
    """TIMEPIXR=1 means TIME is the end of the bin: edges must be
    [TIME - TIMEDEL, TIME], not [TIME, TIME + TIMEDEL]. Getting this
    backwards shifts every light curve by one bin width.
    """
    path = tmp_path / 'cbd.fits'
    time, _ = make_cbd_fits(path)
    cbd = BurstCubeCBD.open(path)
    np.testing.assert_allclose(cbd.data.tstart, time - TIMEDEL_CBD)
    np.testing.assert_allclose(cbd.data.tstop, time)


def test_detector_and_ebounds(tmp_path):
    """The detector name comes from INSTRUME, and the 16-channel energy axis
    must be the CALDB eb16 edges for that specific detector (not eb16 for a
    different detector, and not the raw channel numbers).
    """
    path = tmp_path / 'cbd.fits'
    make_cbd_fits(path, detector='CS0')
    cbd = BurstCubeCBD.open(path)
    assert cbd.detector == 'CS0'
    assert cbd.num_chans == 16
    np.testing.assert_allclose(cbd.ebounds.low_edges()[1], 28.48, atol=1e-2)
    np.testing.assert_allclose(cbd.ebounds.high_edges()[-1], 1648.51, atol=1e-2)


def test_gap_is_preserved_not_filled(tmp_path):
    """A deliberate gap in the time grid must show up as two separate GTI
    segments and must never be bridged, filled, or interpolated across.
    """
    n = 20
    time = 107629263.33 + (np.arange(n) + 1) * TIMEDEL_CBD
    time[10:] += 100.0  # a 100 s gap after the 10th bin
    gti = [(time[0] - TIMEDEL_CBD, time[9]), (time[10] - TIMEDEL_CBD, time[-1])]

    path = tmp_path / 'cbd_gap.fits'
    make_cbd_fits(path, time=time, gti=gti)
    cbd = BurstCubeCBD.open(path)

    assert cbd.gti.num_intervals == 2
    np.testing.assert_allclose(cbd.gti.as_list(), gti)
    # the underlying bin edges themselves must also show the actual gap,
    # not a synthetic filled-in bin
    dt = np.diff(cbd.data.tstart)
    assert np.max(dt) > 50.0


def test_blended_flag_catches_short_dt(tmp_path):
    """A bin whose TIME arrives less than TIMEDEL after the previous bin's
    TIME must be flagged (archive caveat #9); the first bin is never
    flagged, and normal 0.256 s spacing is never flagged.
    """
    n = 10
    time = 107629263.33 + (np.arange(n) + 1) * TIMEDEL_CBD
    time[5:] += (0.1 - TIMEDEL_CBD)  # index 5 now arrives only 0.1 s after index 4
    path = tmp_path / 'cbd_blend.fits'
    make_cbd_fits(path, time=time)
    cbd = BurstCubeCBD.open(path)

    expected = np.zeros(n, dtype=bool)
    expected[5] = True
    np.testing.assert_array_equal(cbd.blended, expected)
    assert cbd.blended[0] == False


def test_original_time_and_sums_exposed(tmp_path):
    """ORIGINAL_TIME, TIME_SYST_ERROR, SUMTOT, and SUMCH_2_15 must be
    exposed and match the file's own columns.
    """
    path = tmp_path / 'cbd.fits'
    time, counts = make_cbd_fits(path)
    cbd = BurstCubeCBD.open(path)
    np.testing.assert_allclose(cbd.original_time, time - 0.01)
    np.testing.assert_allclose(cbd.time_syst_error, 0.05)
    np.testing.assert_array_equal(cbd.sumtot, counts.sum(axis=1))
    np.testing.assert_array_equal(cbd.sumch_2_15, counts[:, 2:16].sum(axis=1))


def test_auxiliary_columns_are_none_after_slice(tmp_path):
    """The per-bin diagnostic columns have no meaning after a slice/rebin
    (bins are combined or dropped), so they must be None on the result.
    """
    path = tmp_path / 'cbd.fits'
    make_cbd_fits(path)
    cbd = BurstCubeCBD.open(path)
    sliced = cbd.slice_time([cbd.time_range])
    assert sliced.original_time is None
    assert sliced.blended is None


def test_combine_by_factor_conserves_counts_within_a_contiguous_segment(tmp_path):
    """combine_by_factor must conserve the total counts and produce the
    expected halved number of bins, when the contiguous segment length is
    an exact multiple of the bin factor.
    """
    n = 20  # one contiguous segment, divisible by 2
    path = tmp_path / 'cbd.fits'
    time, counts = make_cbd_fits(path)
    cbd = BurstCubeCBD.open(path)

    rebinned = cbd.rebin_time(combine_by_factor, 2)
    assert rebinned.data.size == (n // 2, 16)
    assert rebinned.data.counts.sum() == counts.sum()
    # the new edges must be every-other original edge
    np.testing.assert_allclose(rebinned.data.tstop, time[1::2])


def test_slicing_across_a_gap_does_not_fill_it(tmp_path):
    """Slicing over a time range that spans a gap must not fill or
    interpolate counts into the gap; the resulting GTI must still show two
    segments.
    """
    n = 20
    time = 107629263.33 + (np.arange(n) + 1) * TIMEDEL_CBD
    time[10:] += 100.0
    gti = [(time[0] - TIMEDEL_CBD, time[9]), (time[10] - TIMEDEL_CBD, time[-1])]
    path = tmp_path / 'cbd_gap.fits'
    make_cbd_fits(path, time=time, gti=gti)
    cbd = BurstCubeCBD.open(path)

    sliced = cbd.slice_time([(time[0] - TIMEDEL_CBD, time[-1])])
    assert sliced.gti.num_intervals == 2
    assert sliced.data.size[0] == n  # no bins were fabricated to fill the gap


def test_cbd_gti_schema_uf_variant_has_extra_columns(tmp_path):
    """The unfiltered CBD STDGTI has 6 columns rather than 2 and lacks
    HDUCLAS2/HDUVERS/TIMEZERO. open() must detect that from the extension
    and use the matching header template, so the file reads cleanly with no
    warnings -- reading only START/STOP, which both schemas share.
    """
    path = tmp_path / 'cbd_uf.fits'
    make_cbd_fits(path, gti_schema='uf')
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        cbd = BurstCubeCBD.open(path)
    assert isinstance(cbd.headers, CbdUnfilteredHeaders)
    assert cbd.gti.num_intervals == 1


def test_cbd_gti_schema_cl_variant_uses_standard_headers(tmp_path):
    """The cleaned variant's 2-column STDGTI selects the standard template,
    likewise without warnings."""
    path = tmp_path / 'cbd_cl.fits'
    make_cbd_fits(path, gti_schema='cl')
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        cbd = BurstCubeCBD.open(path)
    assert isinstance(cbd.headers, CbdHeaders)
    assert not isinstance(cbd.headers, CbdUnfilteredHeaders)


def test_bin_edges_snap_only_float_noise_not_real_short_bins(tmp_path):
    """At BurstCube MET magnitudes (~1e8 s) a float64 ULP is ~1.5e-8 s, so
    deriving tstart = TIME - TIMEDEL lands a ULP off the previous TIME even
    for genuinely back-to-back bins; gdt-core needs exact equality to treat
    them as contiguous. The snap must close that gap without touching the
    real 1 ms short bins, which sit five orders of magnitude away.
    """
    t0 = 107629263.33
    time = t0 + (np.arange(6) + 1) * TIMEDEL_CBD
    time[4:] -= 1e-3          # a genuine 1 ms short bin at index 4
    path = tmp_path / 'cbd_noise.fits'
    make_cbd_fits(path, time=time)

    cbd = BurstCubeCBD.open(path)
    d = cbd.data
    resid = d.tstart[1:] - d.tstop[:-1]

    # float-noise boundaries are snapped to exactly contiguous
    assert np.all(resid[:3] == 0.0)

    # The real 1 ms short bin is NOT treated as noise. Its start would fall
    # before the previous bin's end, so it is clamped to that end -- bins may
    # not overlap, and gdt-core refuses to rebin them if they do. The feature
    # survives as a shortened exposure and in the `blended` flag, rather than
    # being absorbed into the snap.
    assert resid[3] == 0.0
    assert d.exposure[4] == pytest.approx(0.255, abs=1e-7)
    assert d.exposure[3] == pytest.approx(TIMEDEL_CBD, abs=1e-7)
    assert cbd.blended[4]
    assert np.all(d.exposure > 0)


def test_short_bins_flagged_consistently(tmp_path):
    """All bins short by the same amount must be flagged the same way. A
    tolerance set at the feature size itself would split them by rounding.
    """
    time = 107629263.33 + (np.arange(8) + 1) * TIMEDEL_CBD
    time[2:] -= 1e-3
    time[5:] -= 1e-3
    path = tmp_path / 'cbd_short.fits'
    make_cbd_fits(path, time=time)

    cbd = BurstCubeCBD.open(path)
    assert cbd.blended.sum() == 2
    assert cbd.blended[2] and cbd.blended[5]


def test_sumtot_absent_in_unfiltered_files(tmp_path):
    """SUMTOT/SUMCH_2_15 appear only in cleaned files; the unfiltered ones
    carry TIME, COUNTS, ORIGINAL_TIME and TIME_SYST_ERROR only."""
    path = tmp_path / 'cbd_nosums.fits'
    make_cbd_fits(path, with_sums=False)
    cbd = BurstCubeCBD.open(path)
    assert cbd.sumtot is None
    assert cbd.sumch_2_15 is None
    assert cbd.data.counts.sum() > 0


def test_sumtot_present_in_cleaned_files(tmp_path):
    path = tmp_path / 'cbd_sums.fits'
    make_cbd_fits(path, with_sums=True)
    cbd = BurstCubeCBD.open(path)
    assert cbd.sumtot is not None
    assert np.array_equal(cbd.sumtot, cbd.data.counts.sum(axis=1))


def test_rebin_works_when_short_bins_would_otherwise_overlap(tmp_path):
    """Regression: bins short of TIMEDEL leave tstart[i] before tstop[i-1],
    and gdt-core refuses to merge overlapping bins, so rebin_time raised
    ValueError on every real CBD file. slice_time succeeded regardless, which
    is why it went unnoticed until rebinning was exercised.
    """
    n = 16
    time = 107629263.33 + (np.arange(n) + 1) * TIMEDEL_CBD
    time[4:] -= 1e-3
    time[9:] -= 1e-3
    path = tmp_path / 'cbd_overlap.fits'
    make_cbd_fits(path, time=time)

    cbd = BurstCubeCBD.open(path)
    assert np.all(cbd.data.tstart[1:] >= cbd.data.tstop[:-1])

    rebinned = cbd.rebin_time(combine_by_factor, 2)
    assert rebinned.data.num_times == n // 2
    assert rebinned.data.counts.sum() == cbd.data.counts.sum()
