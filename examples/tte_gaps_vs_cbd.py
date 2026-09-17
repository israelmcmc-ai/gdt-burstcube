#!/usr/bin/env python
"""Minimal reproduction of the README caveat *TTE stops writing while the
detector keeps counting*.

BurstCube TTE covers its nominal span in short recording blocks separated by
gaps of comparable length. Across a gap the event list is empty, so a
TTE-only light curve reads zero and nothing in the file says whether the sky
was quiet or TTE simply wrote nothing. CBD watches the same detector
continuously and settles it.

This downloads one day's cleaned CBD and TTE for one detector, finds the TTE
recording blocks, and plots the two rates against each other with the gaps
shaded. It also prints the CBD and TTE rates inside the blocks and inside
the gaps.

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
from gdt.core.binning.unbinned import bin_by_time
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

    # Count TTE events per CBD bin, so both instruments are measured over
    # exactly the same intervals and no rate scaling is needed. The bin is
    # half-open, [tstart, tstop): CBD bins are contiguous, so counting the
    # closing edge too would assign an event landing exactly on a shared
    # edge to both neighbours -- which happens here, for the one event that
    # opens a recording block on a bin boundary.
    data = cbd.data
    times = np.sort(tte.data.times)
    cbd_counts = data.counts.sum(axis=1)
    tte_counts = (np.searchsorted(times, data.tstop, 'left')
                  - np.searchsorted(times, data.tstart, 'left'))

    print()
    for label, gti in [('inside a TTE block', blocks),
                       ('inside a TTE gap', gaps)]:
        mask = bins_inside(data, gti.as_list())
        exposure = data.exposure[mask].sum()
        print(f'  {label:20s} {mask.sum():4d} bins {exposure:7.2f} s   '
              f'CBD {cbd_counts[mask].sum() / exposure:6.1f} ct/s   '
              f'TTE {tte_counts[mask].sum() / exposure:6.1f} ct/s')

    tte_lc = tte.to_phaii(bin_by_time, 0.256,
                          phaii_class=Phaii).to_lightcurve()
    cbd_lc = cbd.slice_time([(t0, t1)]).to_lightcurve()

    fig, axes = plt.subplots(3, 1, figsize=(11, 8))
    for ax, lc, color, title in [
            (axes[0], cbd_lc, 'C0', 'CBD _cl, 0.256 s bins'),
            (axes[1], tte_lc, 'C1', 'TTE binned to 0.256 s')]:
        shade(ax, gaps, t0)
        ax.step(lc.lo_edges - t0, lc.rates, where='post', color=color)
        ax.set_title(f'{args.detector} {args.obs_id} -- {title}')
        ax.set_ylabel('Rate (ct/s)')
        ax.set_xlim(0, t1 - t0)
    axes[1].set_xlabel(f'Time since MET {t0:.3f} s '
                       '(shaded: TTE recording gaps)')

    # The full span packs ~90 gaps into one axis; overlay the first 20 s so
    # the correspondence is legible bin by bin.
    shade(axes[2], gaps, t0)
    axes[2].step(cbd_lc.lo_edges - t0, cbd_lc.rates, where='post',
                 color='C0', label='CBD _cl')
    axes[2].step(tte_lc.lo_edges - t0, tte_lc.rates, where='post',
                 color='C1', label='TTE')
    axes[2].set_xlim(0, 20)
    axes[2].set_title('first 20 s, overlaid: TTE is exactly zero in every '
                      'shaded gap, CBD is not')
    axes[2].set_xlabel(f'Time since MET {t0:.3f} s')
    axes[2].set_ylabel('Rate (ct/s)')
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
