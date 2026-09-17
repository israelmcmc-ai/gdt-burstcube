#!/usr/bin/env python
"""Minimal reproduction of the README caveat *TTE gaps*.

BurstCube TTE covers its nominal span in short recording blocks separated by
gaps of comparable length. Across a gap the event list is empty, so a
TTE-only light curve reads zero and nothing in the file says whether that
was a quiet sky or TTE not writing. CBD is the same detector and the same
event stream, just binned on board instead of written out event by event,
so it settles the question: it runs continuously through the gaps.

This downloads one day's cleaned CBD and TTE for one detector, finds the TTE
recording blocks, and plots the two rates against each other with the gaps
shaded. TTE is binned on CBD's own bin edges rather than a fresh uniform
grid, so each pair of bins covers exactly the same interval and the panels
can be read against each other bin for bin. It also prints the CBD and TTE
rates inside the blocks and inside the gaps, and the ratio between them as
a function of rate.

Only cleaned (``_cl``) CBD is used, so nothing here depends on unfiltered
data.

Usage::

    python examples/tte_gaps_vs_cbd.py [--obs-id 240814] [--detector CS0]
                                       [--outfile tte_gaps_vs_cbd.png]
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from gdt.core import data_path
from gdt.core.binning.unbinned import bin_by_edges
from gdt.core.phaii import Phaii

from gdt.missions.burstcube.cbd import BurstCubeCBD
from gdt.missions.burstcube.finders import BurstCubeObsFinder
from gdt.missions.burstcube.gti import complement
from gdt.missions.burstcube.tte import BurstCubeTTE


def bins_inside(data, intervals):
    """Which CBD bins lie entirely inside one of `intervals`.

    Bins straddling an edge are excluded rather than split, so the block and
    gap samples stay disjoint and each bin's own exposure stays exact.

    Args:
        data (:class:`~gdt.core.data_primitives.TimeEnergyBins`): The CBD data
        intervals ([(float, float), ...]): The intervals to test against

    Returns:
        (np.array): A boolean mask over the bins
    """
    mask = np.zeros(data.tstart.size, dtype=bool)
    for start, stop in intervals:
        mask |= (data.tstart >= start) & (data.tstop <= stop)
    return mask


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--obs-id', default='240814',
                        help='observation day, YYMMDD (default: 240814)')
    parser.add_argument('--detector', default='CS0',
                        help='detector name (default: CS0)')
    parser.add_argument('--outfile', default='tte_gaps_vs_cbd.png',
                        help='where to write the plot')
    args = parser.parse_args()

    download_dir = Path(data_path) / 'burstcube' / args.obs_id
    download_dir.mkdir(parents=True, exist_ok=True)

    finder = BurstCubeObsFinder()
    finder.cd(args.obs_id)
    cbd_path, = finder.get_cbd(download_dir, detectors=args.detector,
                               variant='cl', verbose=False)
    tte_path, = finder.get_tte(download_dir, detectors=args.detector,
                               verbose=False)

    cbd = BurstCubeCBD.open(cbd_path)
    tte = BurstCubeTTE.open(tte_path)

    t0, t1 = tte.time_range
    blocks = tte.recording_blocks()
    gaps = complement(blocks, t0, t1)
    live = sum(stop - start for start, stop in blocks.as_list())

    print(f'{args.obs_id} {args.detector}: {tte.data.size} TTE events in '
          f'{blocks.num_intervals} recording blocks')
    print(f'  {live:.1f} s of live time in a {t1 - t0:.1f} s span '
          f'(a factor of {(t1 - t0) / live:.2f})')

    # Bin TTE on CBD's own bin edges, not a fresh uniform grid: a grid
    # anchored on the first event lands each TTE bin across two CBD bins,
    # and then neither the plot nor the ratio compares like with like.
    data = cbd.slice_time([(t0, t1)]).data
    edges = np.append(data.tstart, data.tstop[-1])
    tte_binned = tte.to_phaii(bin_by_edges, edges, phaii_class=Phaii)

    cbd_counts = data.counts.sum(axis=1)
    tte_counts = tte_binned.data.counts.sum(axis=1)

    # Two of CBD's own bin boundaries are not contiguous -- _cl drops
    # intervals -- so two of these bins bridge a CBD gap and are wider than
    # the rest. They are not a like-for-like pair, so leave them out.
    straddles_a_cbd_gap = np.append(data.tstart[1:] - data.tstop[:-1] > 1e-9,
                                    False)

    print()
    for label, gti in [('inside a TTE block', blocks),
                       ('inside a TTE gap', gaps)]:
        mask = bins_inside(data, gti.as_list()) & ~straddles_a_cbd_gap
        exposure = data.exposure[mask].sum()
        print(f'  {label:20s} {mask.sum():4d} bins {exposure:7.2f} s   '
              f'CBD {cbd_counts[mask].sum() / exposure:6.1f} ct/s   '
              f'TTE {tte_counts[mask].sum() / exposure:6.1f} ct/s')

    # Inside the blocks the two do not quite agree in total, and the
    # residual is not spread evenly: it is confined to the brightest bins.
    print('\n  CBD/TTE inside a block, by CBD rate:')
    in_block = bins_inside(data, blocks.as_list()) & ~straddles_a_cbd_gap
    cbd_rate = cbd_counts / data.exposure
    for lo, hi in [(0, 100), (100, 200), (200, 300), (300, np.inf)]:
        mask = in_block & (cbd_rate >= lo) & (cbd_rate < hi)
        if not mask.any():
            continue
        cbd_total, tte_total = cbd_counts[mask].sum(), tte_counts[mask].sum()
        ratio = cbd_total / tte_total
        error = ratio * np.sqrt(1 / cbd_total + 1 / tte_total)
        print(f'    {lo:5.0f}-{hi:<5.0f} ct/s  {mask.sum():4d} bins '
              f'{data.exposure[mask].sum():6.2f} s   CBD {cbd_total:6d}   '
              f'TTE {tte_total:6.0f}   ratio {ratio:.4f} +/- {error:.4f}')

    # Blank the bridging bins rather than draw them, in both panels.
    plot_edges = edges.copy()
    cbd_rates = np.where(straddles_a_cbd_gap, np.nan, cbd_rate)
    tte_rates = np.where(straddles_a_cbd_gap, np.nan,
                         tte_counts / data.exposure)

    fig, axes = plt.subplots(3, 1, figsize=(11, 8))
    for ax, rates, color, title in [
            (axes[0], cbd_rates, 'C0', 'CBD _cl, native bins'),
            (axes[1], tte_rates, 'C1', "TTE binned on CBD's own bin edges")]:
        shade(ax, gaps, t0)
        ax.stairs(rates, plot_edges - t0, color=color)
        ax.set_title(f'{args.detector} {args.obs_id} -- {title}')
        ax.set_ylabel('Rate (ct/s)')
        ax.set_xlim(0, t1 - t0)
    axes[1].set_xlabel(f'Time since MET {t0:.3f} s '
                       '(shaded: TTE recording gaps)')

    # The full span packs ~90 gaps into one axis; overlay the first 20 s so
    # the correspondence is legible bin by bin.
    shade(axes[2], gaps, t0)
    axes[2].stairs(cbd_rates, plot_edges - t0, color='C0', label='CBD _cl')
    axes[2].stairs(tte_rates, plot_edges - t0, color='C1', label='TTE')
    axes[2].set_xlim(0, 20)
    axes[2].set_title('first 20 s, overlaid on identical bins: TTE is exactly '
                      'zero in every shaded gap, CBD is not')
    axes[2].set_xlabel(f'Time since MET {t0:.3f} s')
    axes[2].set_ylabel('Rate (ct/s)')
    axes[2].set_ylim(bottom=-12)   # so TTE's zero sits off the axis line
    axes[2].legend(loc='upper right')

    fig.tight_layout()
    fig.savefig(args.outfile, dpi=130)
    print(f'\nwrote {args.outfile}')


def shade(ax, gti, t0):
    """Shade each interval of `gti`, relative to `t0`."""
    for start, stop in gti.as_list():
        ax.axvspan(start - t0, stop - t0, color='0.85', zorder=0)


if __name__ == '__main__':
    main()
