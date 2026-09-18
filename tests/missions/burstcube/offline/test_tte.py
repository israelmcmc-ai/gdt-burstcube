"""Offline tests for gdt.missions.burstcube.tte."""
import warnings

import numpy as np
import pytest
from gdt.core.binning.unbinned import bin_by_time
from gdt.core.phaii import Phaii

from gdt.missions.burstcube.gti import BurstCubeGTI
from gdt.missions.burstcube.tte import BurstCubeTTE

from .conftest import TIMEDEL_CBD, make_tte_fits


def test_open_reads_1024_channels_and_correct_detector(tmp_path):
    path = tmp_path / 'tte.fits'
    times, channels = make_tte_fits(path, detector='CS0', broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(path)
    assert tte.detector == 'CS0'
    assert tte.num_chans == 1024
    np.testing.assert_allclose(np.sort(tte.data.times), np.sort(times))


def test_broken_tstart_tstop_does_not_crash_and_warns(tmp_path):
    """The archive's broken EVENTS TSTART/TSTOP (strings, TSTOP < TSTART)
    must not crash open(), must warn, and the GTI must come from STDGTI
    (i.e. reflect the actual event time range), never a negative exposure.
    """
    path = tmp_path / 'tte_broken.fits'
    times, channels = make_tte_fits(path, broken_tstart_tstop=True)
    with pytest.warns(UserWarning, match='TSTOP < TSTART'):
        tte = BurstCubeTTE.open(path)

    assert tte.gti.range[0] == pytest.approx(times.min())
    assert tte.gti.range[1] == pytest.approx(times.max())
    assert tte.get_exposure() >= 0.0


def test_well_formed_tstart_tstop_does_not_warn(tmp_path):
    path = tmp_path / 'tte_ok.fits'
    make_tte_fits(path, broken_tstart_tstop=False)
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter('always')
        BurstCubeTTE.open(path)
    user_warnings = [w for w in record if issubclass(w.category, UserWarning)]
    assert not user_warnings


def test_to_64_channels_preserves_events_and_top_energy_edge(tmp_path):
    """to_64_channels must not lose or gain events, and the top-channel
    energy edge must be unchanged (it is the same CALDB detector, just
    regrouped).
    """
    path = tmp_path / 'tte.fits'
    make_tte_fits(path, broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(path)

    tte64 = tte.to_64_channels()
    assert tte64.num_chans == 64
    assert tte64.data.size == tte.data.size
    assert tte64.ebounds.high_edges()[-1] == pytest.approx(
        tte.ebounds.high_edges()[-1])
    assert tte64.ebounds.low_edges()[0] == pytest.approx(0.0)


def test_bin_by_time_conserves_total_counts(tmp_path):
    """TTE binned with bin_by_time -> Phaii must have a total count equal to
    the number of events in range.
    """
    path = tmp_path / 'tte.fits'
    times, channels = make_tte_fits(path, broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(path)

    phaii = tte.to_phaii(bin_by_time, TIMEDEL_CBD, phaii_class=Phaii)
    assert phaii.data.counts.sum() == tte.data.size


def test_event_at_known_time_lands_in_expected_bin(tmp_path):
    """A single event placed at a known time must land in the bin whose
    [tstart, tstop) contains it.
    """
    t0 = 107629263.33
    times = np.array([t0 + 0.1])
    channels = np.array([5], dtype=np.int16)
    path = tmp_path / 'tte_single.fits'
    make_tte_fits(path, times=times, channels=channels, broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(path)

    phaii = tte.to_phaii(bin_by_time, TIMEDEL_CBD, tstart=t0, tstop=t0 + 1.0,
                         phaii_class=Phaii)
    # the event at t0+0.1 must land in the first bin [t0, t0+0.256)
    assert phaii.data.counts[0].sum() == 1
    assert phaii.data.counts[1:].sum() == 0


def test_event_exactly_on_a_bin_edge(tmp_path):
    """An event exactly on a bin edge must land on the documented side.
    gdt-core's EventList.bin() uses numpy.histogram2d, whose bins are
    half-open [lo, hi) except the final bin, which is closed -- so an event
    exactly at an interior edge lands in the bin that edge starts, not the
    one it ends.
    """
    t0 = 107629263.33
    edge = t0 + TIMEDEL_CBD  # exactly the boundary between bin 0 and bin 1
    times = np.array([edge])
    channels = np.array([5], dtype=np.int16)
    path = tmp_path / 'tte_edge.fits'
    make_tte_fits(path, times=times, channels=channels, broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(path)

    phaii = tte.to_phaii(bin_by_time, TIMEDEL_CBD, tstart=t0, tstop=t0 + 1.0,
                         phaii_class=Phaii)
    assert phaii.data.counts[0].sum() == 0
    assert phaii.data.counts[1].sum() == 1


def test_tte_binned_at_0_256_reproduces_synthetic_cbd(tmp_path):
    """TTE binned at 0.256 s must reproduce a synthetic CBD built from the
    same events: same per-bin total counts.
    """
    t0 = 107629263.33
    n_bins = 8
    rng = np.random.default_rng(42)

    # place events uniformly within each of n_bins 0.256 s windows
    times = []
    for i in range(n_bins):
        lo, hi = t0 + i * TIMEDEL_CBD, t0 + (i + 1) * TIMEDEL_CBD
        num = rng.integers(1, 6)
        times.extend(rng.uniform(lo, hi - 1e-6, size=num))
    times = np.sort(np.array(times))
    channels = rng.integers(0, 1024, size=times.size).astype(np.int16)

    tte_path = tmp_path / 'tte_for_roundtrip.fits'
    make_tte_fits(tte_path, times=times, channels=channels, broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(tte_path)

    phaii = tte.to_phaii(bin_by_time, TIMEDEL_CBD, tstart=t0,
                         tstop=t0 + n_bins * TIMEDEL_CBD, phaii_class=Phaii)
    tte_counts_per_bin = phaii.data.counts.sum(axis=1)

    # build the "CBD" bin edges honoring TIMEPIXR=1 (end-of-bin TIME) and
    # histogram the same events by hand as the independent comparison
    cbd_time = t0 + (np.arange(n_bins) + 1) * TIMEDEL_CBD
    counts_by_hand = np.histogram(times, bins=np.concatenate(
        ([t0], cbd_time)))[0]

    np.testing.assert_array_equal(tte_counts_per_bin, counts_by_hand)
    assert tte_counts_per_bin.sum() == times.size


def test_slice_time_keeps_events_sharing_a_timestamp(tmp_path):
    """BurstCube timestamps events at 10 us, so many distinct photons share a
    time. gdt-core's PhotonList.slice_time reassembles segments with
    EventList.merge(force_unique=True), which treats those as duplicates and
    drops all but one -- about 74% of the counts in a real file. The override
    must keep every event.
    """
    n_per_stamp = 4
    stamps = 107629263.5 + np.arange(50) * 0.01
    times = np.repeat(stamps, n_per_stamp)
    channels = np.tile(np.arange(n_per_stamp) + 100, len(stamps))
    path = tmp_path / 'tte_dup_times.fits'
    make_tte_fits(path, times=times, channels=channels)

    tte = BurstCubeTTE.open(path)
    assert tte.data.size == len(times)
    assert len(np.unique(tte.data.times)) == len(stamps)  # the collision

    sliced = tte.slice_time([tte.time_range])
    assert sliced.data.size == len(times)


def test_to_phaii_with_time_range_keeps_all_counts(tmp_path):
    """The same loss reached users through to_phaii(time_range=...), which
    calls slice_time internally."""
    stamps = 107629263.5 + np.arange(40) * 0.01
    times = np.repeat(stamps, 3)
    channels = np.tile([10, 20, 30], len(stamps))
    path = tmp_path / 'tte_dup_phaii.fits'
    make_tte_fits(path, times=times, channels=channels)

    tte = BurstCubeTTE.open(path)
    binned = tte.to_phaii(bin_by_time, 0.256, time_range=[tte.time_range])
    assert binned.data.counts.sum() == len(times)

    unrestricted = tte.to_phaii(bin_by_time, 0.256)
    assert binned.data.counts.sum() == unrestricted.data.counts.sum()


def test_slice_time_rejects_overlapping_ranges(tmp_path):
    """Overlapping ranges would double-count, so they are refused rather than
    silently deduplicated."""
    path = tmp_path / 'tte_overlap.fits'
    make_tte_fits(path)
    tte = BurstCubeTTE.open(path)
    t0 = tte.time_range[0]
    with pytest.raises(ValueError, match='must not overlap'):
        tte.slice_time([(t0, t0 + 10.0), (t0 + 5.0, t0 + 15.0)])


def test_slice_time_accepts_disjoint_ranges(tmp_path):
    """Disjoint ranges are fine and keep every event in either range."""
    times = 107629263.5 + np.arange(300) * 0.1
    path = tmp_path / 'tte_disjoint.fits'
    make_tte_fits(path, times=times, channels=np.full(len(times), 100))
    tte = BurstCubeTTE.open(path)
    t0 = times[0]
    r1, r2 = (t0, t0 + 5.0), (t0 + 10.0, t0 + 15.0)
    expected = int((((times >= r1[0]) & (times <= r1[1])) |
                    ((times >= r2[0]) & (times <= r2[1]))).sum())
    assert tte.slice_time([r1, r2]).data.size == expected


def test_slice_time_handles_disjoint_ranges_where_some_are_empty(tmp_path):
    """TTE coverage is clustered, so a reasonable set of windows can leave
    some containing no events. An empty segment has a time_range of None,
    which Gti.from_list rejects, and EventList.merge reduces over each
    segment's times and so fails on an empty one. Both must be handled.
    """
    times = 107629263.5 + np.concatenate([np.arange(50) * 0.01,
                                          100.0 + np.arange(50) * 0.01])
    path = tmp_path / 'tte_clustered.fits'
    make_tte_fits(path, times=times, channels=np.full(len(times), 100))
    tte = BurstCubeTTE.open(path)

    t0 = times[0]
    ranges = [(t0, t0 + 1.0),          # populated
              (t0 + 50.0, t0 + 51.0),  # empty
              (t0 + 100.0, t0 + 101.0)]
    expected = int(sum(((times >= a) & (times <= b)).sum() for a, b in ranges))
    assert tte.slice_time(ranges).data.size == expected


def test_slice_time_raises_when_no_range_contains_events(tmp_path):
    path = tmp_path / 'tte_none.fits'
    make_tte_fits(path)
    tte = BurstCubeTTE.open(path)
    far = tte.time_range[1] + 1e6
    with pytest.raises(ValueError, match='contain any events'):
        tte.slice_time([(far, far + 10.0)])


def test_recording_blocks_splits_on_gaps_and_reports_live_time(tmp_path):
    """Three deliberate blocks separated by 2 s gaps, each block made of
    events 10 ms apart. The blocks must come back exactly, and their total
    must be live time -- strictly less than the elapsed span.
    """
    t0 = 107629263.0
    times = np.concatenate([t0 + offset + np.arange(5) * 0.01
                            for offset in (0.0, 2.0, 4.0)])
    path = tmp_path / 'tte_blocks.fits'
    make_tte_fits(path, times=times, broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(path)

    blocks = tte.recording_blocks()
    assert blocks.num_intervals == 3
    np.testing.assert_allclose(
        blocks.as_list(),
        [(t0, t0 + 0.04), (t0 + 2.0, t0 + 2.04), (t0 + 4.0, t0 + 4.04)])

    live = sum(stop - start for start, stop in blocks.as_list())
    assert live == pytest.approx(0.12)
    assert live < tte.time_range[1] - tte.time_range[0]


def test_recording_blocks_gap_threshold_changes_the_split(tmp_path):
    """A threshold above the inter-block spacing must coalesce everything
    into one block; below the in-block spacing, every event is its own.
    """
    t0 = 107629263.0
    times = np.concatenate([t0 + offset + np.arange(5) * 0.01
                            for offset in (0.0, 2.0, 4.0)])
    path = tmp_path / 'tte_blocks.fits'
    make_tte_fits(path, times=times, broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(path)

    assert tte.recording_blocks(gap_threshold=3.0).num_intervals == 1
    assert tte.recording_blocks(gap_threshold=0.001).num_intervals == times.size


def test_recording_blocks_complement_is_the_gaps(tmp_path):
    """The gaps are the complement of the blocks over the file's time range,
    which is how the README and the example script derive them.
    """
    t0 = 107629263.0
    times = np.concatenate([t0 + offset + np.arange(5) * 0.01
                            for offset in (0.0, 2.0, 4.0)])
    path = tmp_path / 'tte_blocks.fits'
    make_tte_fits(path, times=times, broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(path)

    gaps = BurstCubeGTI.complement(tte.recording_blocks(), *tte.time_range)
    assert gaps.num_intervals == 2
    np.testing.assert_allclose(gaps.as_list(),
                               [(t0 + 0.04, t0 + 2.0), (t0 + 2.04, t0 + 4.0)])


def test_recording_blocks_is_one_block_when_events_are_contiguous(tmp_path):
    path = tmp_path / 'tte.fits'
    times, _ = make_tte_fits(path, broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(path)

    blocks = tte.recording_blocks()
    assert blocks.num_intervals == 1
    assert blocks.range == pytest.approx((times.min(), times.max()))


def test_recording_blocks_rejects_a_nonpositive_threshold(tmp_path):
    path = tmp_path / 'tte.fits'
    make_tte_fits(path, broken_tstart_tstop=False)
    tte = BurstCubeTTE.open(path)

    with pytest.raises(ValueError, match='must be positive'):
        tte.recording_blocks(gap_threshold=0.0)


def test_recording_blocks_is_unaffected_by_event_ordering(tmp_path):
    """Real event times are sorted, but nothing in the reader guarantees it,
    so the method must sort rather than trust the column order.
    """
    t0 = 107629263.0
    ordered = np.concatenate([t0 + offset + np.arange(5) * 0.01
                              for offset in (0.0, 2.0, 4.0)])
    shuffled = np.random.default_rng(0).permutation(ordered)

    paths = []
    for name, times in (('ordered.fits', ordered), ('shuffled.fits', shuffled)):
        path = tmp_path / name
        make_tte_fits(path, times=times, broken_tstart_tstop=False)
        paths.append(path)

    blocks = [BurstCubeTTE.open(p).recording_blocks().as_list() for p in paths]
    np.testing.assert_allclose(blocks[0], blocks[1])
