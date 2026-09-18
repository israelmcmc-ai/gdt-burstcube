"""Data-driven tests for CBD, TTE, orbit, GTI and housekeeping.

Every bug found in this plugin so far was found by running against the real
archive rather than the synthetic fixtures, because the fixtures were built
from the same assumptions as the code and so agreed with it. These tests pin
the quantities that were only ever established by hand: the timestamp
populations that set the CBD snap tolerance, the timestamp collisions that
break a deduplicating merge on TTE, the two STDGTI schemas, and the header
fidelity of a file whose own date keywords are internally inconsistent.

They need the files listed in ``src/gdt/data/burstcube.urls``
(``gdt-data download burstcube``) and skip cleanly without them.
"""
import warnings

import numpy as np
import pytest
from astropy.io import fits

from gdt.core.binning.binned import combine_by_factor
from gdt.core.binning.unbinned import bin_by_edges, bin_by_time
from gdt.core.phaii import Phaii

from gdt.missions.burstcube import caldb
from gdt.missions.burstcube.cbd import BurstCubeCBD
from gdt.missions.burstcube.gti import BurstCubeGTI
from gdt.missions.burstcube.hk import BurstCubeHK
from gdt.missions.burstcube.orbit import BurstCubeOrbit
from gdt.missions.burstcube.tte import BurstCubeTTE
from gdt.missions.burstcube.headers import CBDHeaders, CBDUnfilteredHeaders

from .conftest import real_file

CBD_CL = 'bc240530cs0_3cbd_cl.fits.gz'
CBD_UF = 'bc240530cs0_3cbd_uf.fits.gz'
TTE = 'bc240530cs0_tte_uf.evt.gz'
# 240814 CS0 is the pairing behind the README's TTE recording-gap caveat:
# it is the only archive day with both a quiet CBD baseline and TTE blocks
# long enough to match CBD bins against one for one.
CBD_CL_240814 = 'bc240814cs0_3cbd_cl.fits.gz'
TTE_240814 = 'bc240814cs0_tte_uf.evt.gz'
ORBIT = 'bc240530.hk.gz'
DET_HK = 'bc240530csa.hk.gz'
DET_HK_240814 = 'bc240814csa.hk.gz'
SAA_GTI = 'bc_csa_saa_in_cl.gti'
BIN_GTI = 'bc_cs0_bin0_256.gti'

TIMEDEL = 0.256


def _open_cbd(basename):
    with warnings.catch_warnings():
        warnings.simplefilter('error')       # reading must be warning-free
        return BurstCubeCBD.open(real_file(basename))


# --------------------------------------------------------------------------
# CBD
# --------------------------------------------------------------------------

@pytest.mark.parametrize('basename,expected_bins', [(CBD_UF, 98188), (CBD_CL, 12576)])
def test_cbd_row_counts(basename, expected_bins):
    """Row counts of the real files, so a parsing change that drops or
    duplicates rows is caught."""
    assert _open_cbd(basename).data.num_times == expected_bins


def test_cbd_timestamps_are_bin_ends():
    """TIMEPIXR=1: the TIME column is the END of each bin. Reading it as the
    start shifts every light curve by 0.256 s.
    """
    cbd = _open_cbd(CBD_CL)
    raw = fits.open(real_file(CBD_CL))['CBD'].data['TIME'].astype(float)
    d = cbd.data

    # every bin ENDS at its TIME value, exactly
    assert np.allclose(d.tstop, raw, rtol=0, atol=1e-9)

    # and starts one TIMEDEL earlier, except where that would overlap the
    # previous bin (the 1 ms short bins), which are clamped to the previous
    # stop and so carry a 0.255 s exposure
    nominal = raw - TIMEDEL
    clamped = np.concatenate(([False], nominal[1:] < raw[:-1] - 1e-6))
    assert clamped.sum() == 53
    assert np.allclose(d.tstart[~clamped], nominal[~clamped], rtol=0, atol=1e-6)
    assert np.allclose(d.tstart[1:][clamped[1:]], d.tstop[:-1][clamped[1:]],
                       rtol=0, atol=0)
    assert d.exposure[clamped] == pytest.approx(0.255, abs=1e-6)
    assert d.exposure[~clamped] == pytest.approx(TIMEDEL, abs=1e-6)


def test_cbd_snap_tolerance_sits_between_the_real_populations():
    """The snap tolerance exists to remove float64 artifacts without touching
    real sub-TIMEDEL structure. In this file the two populations are five
    decades apart: artifacts at one ULP (~1.5e-8 s) and genuine bins short by
    1.0e-3 s. Nothing may land in between, or the tolerance is splitting a
    real feature.
    """
    raw = fits.open(real_file(CBD_UF))['CBD'].data['TIME'].astype(float)
    resid = np.abs((raw - TIMEDEL)[1:] - raw[:-1])
    nonzero = resid[resid > 0]

    ulp_scale = np.spacing(raw[0])
    artifacts = nonzero[nonzero < 1e-6]
    assert artifacts.size > 0
    assert np.all(artifacts <= 10 * ulp_scale)

    between = nonzero[(nonzero >= 1e-6) & (nonzero < 1e-4)]
    assert between.size == 0, (
        f'{between.size} residuals lie in the gap the tolerance occupies; '
        'the two populations are no longer cleanly separated')


def test_cbd_short_bins_are_all_flagged_the_same_way():
    """The 1 ms short bins are one physical feature and must be treated
    uniformly. A tolerance set at the feature size itself would flag only the
    subset whose rounding fell on one side.
    """
    cbd = _open_cbd(CBD_UF)
    raw = fits.open(real_file(CBD_UF))['CBD'].data['TIME'].astype(float)
    dt = np.diff(raw)
    short = np.concatenate(([False], np.isclose(dt, 0.255, atol=1e-5)))
    assert short.sum() == 432
    assert np.array_equal(cbd.blended, short)


def test_cbd_gaps_are_preserved_not_filled():
    """The instrument was frequently off; the cleaned file has multi-minute
    gaps. Slicing across them must not fabricate bins.
    """
    cbd = _open_cbd(CBD_CL)
    d = cbd.data
    gaps = (d.tstart[1:] - d.tstop[:-1])
    assert (gaps > 1.0).sum() >= 4               # 120, 180, 840, 1260 s
    sliced = cbd.slice_time([(d.tstart[0], d.tstop[-1])])
    assert sliced.data.num_times == d.num_times


def test_cbd_rebin_runs_on_real_files():
    """Regression: the 1 ms short bins used to leave tstart[i] before
    tstop[i-1], so bins overlapped and gdt-core refused to merge them --
    rebin_time raised ValueError on every real CBD file. slice_time happened
    to succeed, which is why it went unnoticed.
    """
    cbd = _open_cbd(CBD_CL)
    rebinned = cbd.rebin_time(combine_by_factor, 4)
    assert rebinned.data.num_times > 0


def test_cbd_rebin_conserves_counts_within_a_contiguous_run():
    """Counts are conserved exactly when no trailing partial group is
    dropped. combine_by_factor discards the remainder at the end of each
    contiguous segment, so this uses one segment whose length divides the
    factor -- otherwise the test would be measuring gdt-core's remainder
    policy rather than our binning.
    """
    cbd = _open_cbd(CBD_CL)
    d = cbd.data
    boundaries = np.flatnonzero(d.tstart[1:] != d.tstop[:-1])
    first_len = int(boundaries[0]) + 1
    usable = first_len - (first_len % 4)

    run = cbd.slice_time([(d.tstart[0], d.tstop[usable - 1])])
    before = run.data.counts.sum()
    assert run.data.num_times % 4 == 0

    rebinned = run.rebin_time(combine_by_factor, 4)
    assert rebinned.data.counts.sum() == before
    assert rebinned.data.num_times == run.data.num_times // 4


def test_cbd_rebin_only_loses_trailing_partial_groups():
    """Across the whole gappy file the deficit must be exactly the trailing
    bins combine_by_factor drops per contiguous segment -- nothing else."""
    cbd = _open_cbd(CBD_CL)
    d = cbd.data
    boundaries = np.flatnonzero(d.tstart[1:] != d.tstop[:-1])
    seg_lengths = np.diff(np.concatenate(([0], boundaries + 1, [d.num_times])))
    expected_dropped = int((seg_lengths % 2).sum())

    rebinned = cbd.rebin_time(combine_by_factor, 2)
    assert d.num_times - 2 * rebinned.data.num_times == expected_dropped


def test_cbd_sum_columns_present_only_in_cleaned_files():
    """SUMTOT/SUMCH_2_15 exist only in _cl. Reading them unconditionally
    raised KeyError on every _uf file.
    """
    assert _open_cbd(CBD_UF).sumtot is None
    cl = _open_cbd(CBD_CL)
    assert cl.sumtot is not None
    assert np.array_equal(cl.sumtot, cl.data.counts.sum(axis=1))


def test_cbd_gti_schema_detected_from_the_extension():
    """The two variants carry different STDGTI schemas; the reader picks the
    header template by inspecting the extension, not the filename."""
    assert isinstance(_open_cbd(CBD_UF).headers, CBDUnfilteredHeaders)
    cl_headers = _open_cbd(CBD_CL).headers
    assert isinstance(cl_headers, CBDHeaders)
    assert not isinstance(cl_headers, CBDUnfilteredHeaders)


def test_cbd_energy_axis_is_per_detector_caldb():
    """The 16-channel axis comes from CALDB eb16 for this file's detector."""
    cbd = _open_cbd(CBD_CL)
    assert cbd.detector == 'CS0'
    lo = np.asarray(cbd.data.emin, dtype=float)
    assert lo[:3] == pytest.approx([0.0, 28.48, 62.23], abs=1e-2)
    assert np.asarray(cbd.data.emax, dtype=float)[-1] == pytest.approx(1648.51, abs=1e-2)


# --------------------------------------------------------------------------
# TTE
# --------------------------------------------------------------------------

def _open_tte():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')      # the broken-header warning is expected
        return BurstCubeTTE.open(real_file(TTE))


def test_tte_broken_time_keywords_are_recovered():
    """TSTART/TSTOP are strings with TSTOP < TSTART. The reader must warn,
    fall back to STDGTI and the event times, and never yield a negative span.
    The recovered span should match the file's own ONTIME.
    """
    with pytest.warns(UserWarning, match='TSTOP'):
        tte = BurstCubeTTE.open(real_file(TTE))

    start, stop = tte.time_range
    assert stop > start
    ontime = fits.open(real_file(TTE))['EVENTS'].header['ONTIME']
    assert (stop - start) == pytest.approx(ontime, abs=1e-3)


def test_tte_timestamp_collisions_exist_in_the_real_file():
    """The premise behind the slice_time override: BurstCube timestamps at
    10 us, so distinct photons share a time. If this stops being true the
    override is still correct but no longer load-bearing.
    """
    ev = _open_tte().data
    assert ev.size == 6364
    assert len(np.unique(ev.times)) == 1644


def test_tte_slice_time_keeps_every_event():
    """gdt-core's PhotonList.slice_time dedupes by timestamp via
    EventList.merge(force_unique=True), which discarded ~74% of this file.
    """
    tte = _open_tte()
    assert tte.slice_time([tte.time_range]).data.size == 6364


def test_tte_binning_matches_an_independent_histogram():
    """Bin the real events and compare against a plain NumPy histogram."""
    tte = _open_tte()
    start, stop = tte.time_range
    phaii = tte.to_phaii(bin_by_time, TIMEDEL, time_range=[(start, stop)])
    assert phaii.data.counts.sum() == 6364

    edges = np.arange(start, stop + TIMEDEL, TIMEDEL)
    reference, _ = np.histogram(tte.data.times, bins=edges)
    assert phaii.data.counts.sum() == reference.sum()


def test_tte_to_64_channels_preserves_events():
    """The CALDB reb64 regroup changes the channel axis, not the event count."""
    tte = _open_tte()
    regrouped = tte.to_64_channels()
    assert regrouped.data.size == tte.data.size
    assert regrouped.data.channels.max() <= 63


def test_tte_headers_match_the_file_on_disk():
    """This file's three extensions carry identical DATE-OBS/DATE-END but
    different TSTART/TSTOP. Reading must not resync the dates from the timing
    keywords, or the extensions end up disagreeing with each other and with
    the archive.
    """
    disk = {(hdu.name or 'PRIMARY'): (hdu.header.get('DATE-OBS'),
                                      hdu.header.get('DATE-END'))
            for hdu in fits.open(real_file(TTE))}
    tte = _open_tte()
    for name, expected in disk.items():
        header = tte.headers[name]
        assert (header['DATE-OBS'], header['DATE-END']) == expected


# --------------------------------------------------------------------------
# orbit, GTI, housekeeping
# --------------------------------------------------------------------------

def test_orbit_positions_are_physical():
    """Position and velocity magnitudes should match a low-Earth orbit."""
    orbit = BurstCubeOrbit.open(real_file(ORBIT))
    frame = orbit.get_spacecraft_frame()
    radius = np.linalg.norm(frame.obsgeoloc.xyz.to_value('km'), axis=0)
    speed = np.linalg.norm(frame.obsgeovel.xyz.to_value('km/s'), axis=0)
    assert np.all((radius > 6500) & (radius < 7000))
    assert np.all((speed > 7.0) & (speed < 8.0))


def test_orbit_merge_deduplicates_and_is_order_independent():
    """Merging a file with itself must not double its samples, and the result
    must not depend on argument order."""
    orbit = BurstCubeOrbit.open(real_file(ORBIT))
    n = len(orbit.get_spacecraft_frame())

    other = BurstCubeOrbit.open(real_file(ORBIT))
    forward = BurstCubeOrbit.merge(orbit, other)
    backward = BurstCubeOrbit.merge(other, orbit)

    assert len(forward.get_spacecraft_frame()) == n
    assert np.array_equal(forward.get_spacecraft_frame().obstime.value,
                          backward.get_spacecraft_frame().obstime.value)


def test_gti_files_read_and_intersect():
    """Trend GTI files read, and intersecting SAA-in with the binning GTI
    yields intervals contained in both."""
    saa = BurstCubeGTI.open(real_file(SAA_GTI))
    binning = BurstCubeGTI.open(real_file(BIN_GTI))
    assert saa.gti.num_intervals == 586
    assert binning.gti.num_intervals == 107

    combined = BurstCubeGTI.intersect(saa.gti, binning.gti)
    intervals = combined.as_list()
    assert len(intervals) > 0
    # every interval of an intersection must be non-empty and lie inside both
    # inputs; the total duration cannot exceed either input's
    for low, high in intervals:
        assert high > low
    total = sum(high - low for low, high in intervals)
    saa_total = sum(high - low for low, high in saa.gti.as_list())
    bin_total = sum(high - low for low, high in binning.gti.as_list())
    assert total <= min(saa_total, bin_total)


def test_detector_hk_exposes_thresholds_and_enable_flags():
    """HK1 and HK2 carry per-detector arrays indexed 0-3."""
    hk = BurstCubeHK.open(real_file(DET_HK))
    assert len(hk.hk1) > 0
    assert len(hk.hk2) > 0

    # without a detector, the full (n, 4) array; with one, that detector's column
    all_thresholds = hk.peak_threshold()
    assert all_thresholds.shape[1] == 4
    for det in range(4):
        assert np.array_equal(hk.peak_threshold(det), all_thresholds[:, det])

    flags = hk.cbd_enabled()
    assert flags.shape[1] == 4
    assert set(np.unique(flags)).issubset({0, 1})


def test_tte_recording_blocks_are_gaps_cbd_kept_counting_through():
    """The README caveat *TTE gaps*, as a test. On 240814 CS0, TTE arrives in
    ~94 short blocks, and across the gaps between them it records nothing at
    all over roughly 100 s while CBD -- the same detector and the same event
    stream, binned on board -- counts at its usual rate.

    TTE is binned on the CBD file's own bin edges, so each pair of bins
    covers exactly the same interval and no rate scaling is involved.

    This pins the measurement the caveat and
    ``examples/tte_gaps_vs_cbd.py`` report: if the archive is ever
    reprocessed to fill the gaps, this test is how we find out.
    """
    cbd = _open_cbd(CBD_CL_240814)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')      # the known broken TSTART/TSTOP
        tte = BurstCubeTTE.open(real_file(TTE_240814))

    tstart, tstop = tte.time_range
    blocks = tte.recording_blocks()
    gaps = BurstCubeGTI.complement(blocks, tstart, tstop)
    assert blocks.num_intervals == gaps.num_intervals + 1

    live = sum(stop - start for start, stop in blocks.as_list())
    assert live < 0.5 * (tstop - tstart)     # less than half the span is live

    data = cbd.slice_time([(tstart, tstop)]).data
    edges = np.append(data.tstart, data.tstop[-1])
    tte_binned = tte.to_phaii(bin_by_edges, edges, phaii_class=Phaii)
    assert tte_binned.data.counts.sum() == tte.data.size

    cbd_counts = data.counts.sum(axis=1)
    tte_counts = tte_binned.data.counts.sum(axis=1)
    # _cl drops intervals, so a couple of CBD bin boundaries are not
    # contiguous and the bins spanning them are not like-for-like pairs
    straddles_a_cbd_gap = np.append(data.tstart[1:] - data.tstop[:-1] > 1e-9,
                                    False)

    def bins_inside(gti):
        """Bins lying entirely inside one interval of `gti`, so the block and
        gap samples are disjoint and each bin's own exposure is exact."""
        mask = np.zeros(data.tstart.size, dtype=bool)
        for low, high in gti.as_list():
            mask |= (data.tstart >= low) & (data.tstop <= high)
        return mask & ~straddles_a_cbd_gap

    in_block, in_gap = bins_inside(blocks), bins_inside(gaps)
    block_exposure = data.exposure[in_block].sum()
    gap_exposure = data.exposure[in_gap].sum()
    assert block_exposure > 50.0         # enough of each to mean anything
    assert gap_exposure > 50.0

    # across the gaps TTE is empty -- not low, empty -- while CBD is
    # unchanged, and that is a lot of counts TTE never wrote
    assert tte_counts[in_gap].sum() == 0
    cbd_in_gap = cbd_counts[in_gap].sum() / gap_exposure
    cbd_in_block = cbd_counts[in_block].sum() / block_exposure
    assert cbd_in_gap / cbd_in_block == pytest.approx(1.0, abs=0.2)
    assert cbd_in_gap * gap_exposure > 1e4

    # inside the blocks the deficit is confined to the brightest bins: below
    # 200 ct/s the two agree to well under a percent, above it TTE is short
    cbd_rate = cbd_counts / data.exposure
    quiet = in_block & (cbd_rate < 200.0)
    bright = in_block & (cbd_rate >= 200.0)
    assert data.exposure[quiet].sum() > 50.0
    assert data.exposure[bright].sum() > 1.0

    assert (cbd_counts[quiet].sum() / tte_counts[quiet].sum()
            == pytest.approx(1.0, abs=0.01))
    assert cbd_counts[bright].sum() / tte_counts[bright].sum() > 1.03

    # and the whole in-block deficit lives in those bright bins
    deficit = cbd_counts[in_block].sum() - tte_counts[in_block].sum()
    bright_deficit = cbd_counts[bright].sum() - tte_counts[bright].sum()
    assert bright_deficit >= deficit


def test_base_thres_is_the_caveat_7_threshold_and_peak_thres_is_not():
    """``PEAK_THRES`` reads 6.0 mV for every detector on every day and never
    moves; ``BASE_THRES`` is the one archive caveat #7's mid-mission change
    raised, and 240814 catches a transition inside one file.

    Easy to get backwards -- the names suggest otherwise -- so pin it.
    """
    early = BurstCubeHK.open(real_file(DET_HK))
    late = BurstCubeHK.open(real_file(DET_HK_240814))

    for hk in (early, late):
        assert np.array_equal(np.unique(hk.peak_threshold()), [6.0])

    np.testing.assert_array_equal(
        [np.unique(early.base_threshold()[:, d]) for d in range(4)],
        [[82.0], [98.0], [74.0], [74.0]])
    assert [list(np.unique(late.base_threshold()[:, d])) for d in range(4)] == [
        [82.0, 238.0], [98.0, 257.0], [74.0, 247.0], [74.0, 259.0]]


def test_tte_threshold_turn_on_sits_below_caveat_7s_table_by_a_constant_factor():
    """Caveat #7 Table 1 puts CS0's threshold at 26.93 keV before the change
    and 100.12 keV after. The turn-on measured in the TTE spectra is ~20%
    lower -- on *both* sides, so it is a calibration scale offset between
    CALDB's pre-launch ``eb1024`` and whatever the table is on, not scatter.

    This is what notebook 1 reports, and the reason the plugin quotes the
    measured value rather than the table's.
    """
    def turn_on_kev(basename):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')      # broken TSTART/TSTOP
            tte = BurstCubeTTE.open(real_file(basename))
        ebounds = caldb.ebounds(tte.detector, 1024)
        low = np.asarray(ebounds.low_edges())
        width = np.asarray(ebounds.high_edges()) - low
        density = np.bincount(tte.data.channels, minlength=1024) / width
        return low[np.flatnonzero(density >= 0.5 * density.max())[0]]

    before, after = turn_on_kev(TTE), turn_on_kev(TTE_240814)
    table_before, table_after = 26.93, 100.12        # caveat #7, detector 0

    # the change itself is unmistakable: a factor of ~3.7 either way
    assert after / before > 3.0
    assert table_after / table_before > 3.0

    # but each sits ~20% below the table, by the same factor
    ratio_before = before / table_before
    ratio_after = after / table_after
    assert 0.7 < ratio_before < 0.9
    assert ratio_after == pytest.approx(ratio_before, abs=0.05)
