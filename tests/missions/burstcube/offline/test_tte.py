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
"""Offline tests for gdt.missions.burstcube.tte."""
import warnings

import numpy as np
import pytest
from gdt.core.binning.unbinned import bin_by_time
from gdt.core.phaii import Phaii

from gdt.missions.burstcube.tte import BurstCubeTte

from .conftest import TIMEDEL_CBD, make_tte_fits


def test_open_reads_1024_channels_and_correct_detector(tmp_path):
    path = tmp_path / 'tte.fits'
    times, channels = make_tte_fits(path, detector='CS0', broken_tstart_tstop=False)
    tte = BurstCubeTte.open(path)
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
        tte = BurstCubeTte.open(path)

    assert tte.gti.range[0] == pytest.approx(times.min())
    assert tte.gti.range[1] == pytest.approx(times.max())
    assert tte.get_exposure() >= 0.0


def test_well_formed_tstart_tstop_does_not_warn(tmp_path):
    path = tmp_path / 'tte_ok.fits'
    make_tte_fits(path, broken_tstart_tstop=False)
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter('always')
        BurstCubeTte.open(path)
    user_warnings = [w for w in record if issubclass(w.category, UserWarning)]
    assert not user_warnings


def test_to_64_channels_preserves_events_and_top_energy_edge(tmp_path):
    """to_64_channels must not lose or gain events, and the top-channel
    energy edge must be unchanged (it is the same CALDB detector, just
    regrouped).
    """
    path = tmp_path / 'tte.fits'
    make_tte_fits(path, broken_tstart_tstop=False)
    tte = BurstCubeTte.open(path)

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
    tte = BurstCubeTte.open(path)

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
    tte = BurstCubeTte.open(path)

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
    tte = BurstCubeTte.open(path)

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
    tte = BurstCubeTte.open(tte_path)

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

    tte = BurstCubeTte.open(path)
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

    tte = BurstCubeTte.open(path)
    binned = tte.to_phaii(bin_by_time, 0.256, time_range=[tte.time_range])
    assert binned.data.counts.sum() == len(times)

    unrestricted = tte.to_phaii(bin_by_time, 0.256)
    assert binned.data.counts.sum() == unrestricted.data.counts.sum()


def test_slice_time_rejects_overlapping_ranges(tmp_path):
    """Overlapping ranges would double-count, so they are refused rather than
    silently deduplicated."""
    path = tmp_path / 'tte_overlap.fits'
    make_tte_fits(path)
    tte = BurstCubeTte.open(path)
    t0 = tte.time_range[0]
    with pytest.raises(ValueError, match='must not overlap'):
        tte.slice_time([(t0, t0 + 10.0), (t0 + 5.0, t0 + 15.0)])


def test_slice_time_accepts_disjoint_ranges(tmp_path):
    """Disjoint ranges are fine and keep every event in either range."""
    times = 107629263.5 + np.arange(300) * 0.1
    path = tmp_path / 'tte_disjoint.fits'
    make_tte_fits(path, times=times, channels=np.full(len(times), 100))
    tte = BurstCubeTte.open(path)
    t0 = times[0]
    r1, r2 = (t0, t0 + 5.0), (t0 + 10.0, t0 + 15.0)
    expected = int((((times >= r1[0]) & (times <= r1[1])) |
                    ((times >= r2[0]) & (times <= r2[1]))).sum())
    assert tte.slice_time([r1, r2]).data.size == expected
