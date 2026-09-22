#!/usr/bin/env python
"""Minimal standalone plot of the BurstCube light curve around GRB 240629A,
with the archive's own quoted timing uncertainty on a shared time axis.

GBM triggered on this burst (bn240629704) at 2024-06-29T16:53:52.729 UTC.
That instant, converted onto the BurstCube clock using the CBD file's own
MJDREFI/MJDREFF/TIMESYS keywords, is t0 -- the zero of the x axis. The x axis
is BurstCube MET, just shifted so that t0 sits at 0.

That header-derived epoch is the **defective** one: read at face value (per
OGIP, MJDREF is expressed in the scale TIMESYS names), it resolves to
2021-01-01 00:00:00 UTC, which is 37.000 s later than the true MET epoch,
2021-01-01 00:00:00 TAI. That is exactly why the burst appears at
t0+``OBSERVED_OFFSET`` (34.14 s) rather than at t0 itself:
``OBSERVED_OFFSET`` is the offset measured under the archive's *stated*
epoch. Correcting to the true TAI epoch would place the burst at t0-2.86 s
instead -- agreement with the GBM trigger within the true timing
uncertainty, and itself evidence that the archive's quoted
``TIME_SYST_ERROR`` (the lower panel here) is underestimated in at least
some periods. See the README Caveats section and
``gdt.missions.burstcube.time`` for the full evidence; this script
deliberately keeps the archive's own (defective) epoch so that the
discrepancy it exists to show remains visible on the plot.

The lower panel is the TIME_SYST_ERROR column: the uncertainty the archive
itself quotes on each bin's TIME, which is a ground-reconstructed quantity
(see the archive caveats document, caveat #11 "Time Corrections", and the
ORIGINAL_TIME column alongside it).

The tallest native-resolution spike in the plotted window is at t0+91 s and
is a single 0.256 s bin in CS0 alone (a single-detector phosphorescence
event), not a burst -- it is left in rather than clipped, since nothing here
is filtered. (That t0+91 s is in this script's own, uncorrected frame;
against a corrected t0 the same spike sits at +53.98 s, i.e. 37 s earlier in
relative terms. It is not renumbered here, because this script deliberately
keeps the uncorrected frame.)

Uses no GDT package -- only astropy, numpy, matplotlib and the standard
library. URLs are hardcoded; the epoch is not, it comes from the downloaded
header -- deliberately the archive's own (defective) one; see above.

Usage::

    python examples/grb240629a_lightcurve.py [outfile.png]
"""
import gzip
import io
import sys
import urllib.request

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta

BASE = "https://heasarc.gsfc.nasa.gov/FTP/burstcube/data/obs/2024_06/240629/monitor"
CBD_URLS = (
    f"{BASE}/bc240629cs0_3cbd_cl.fits.gz",
    f"{BASE}/bc240629cs1_3cbd_cl.fits.gz",
    f"{BASE}/bc240629cs2_3cbd_cl.fits.gz",
    f"{BASE}/bc240629cs3_3cbd_cl.fits.gz",
)

# Fermi GBM trigger bn240629704, from the GBM burst catalog (fermigbrst).
T0_UTC = "2024-06-29T16:53:52.729"

HALF_WIDTH = 300.0      # seconds each side of t0 to plot
REBIN = 8               # 0.256 s native -> 2.048 s, for legibility

# Where the burst actually appears in this data, under the archive's own
# *stated* (defective) epoch used above -- NOT derived by this script:
# measured by fitting GBM's own burst profile to these counts. Drawn only so
# the offset can be read off the figure. Under the corrected TAI epoch this
# would instead be t0-2.86 s (see the module docstring).
OBSERVED_OFFSET = 34.14


def fetch(url):
    with urllib.request.urlopen(url) as response:
        return response.read()


def read_cbd(url):
    with fits.open(io.BytesIO(gzip.decompress(fetch(url)))) as hdul:
        header = hdul["CBD"].header
        data = hdul["CBD"].data
        return header, data["TIME"].copy(), data["SUMTOT"].copy(), \
            data["TIME_SYST_ERROR"].copy()


headers, times, counts, syst_errors = [], [], [], []
for url in CBD_URLS:
    print(f"downloading {url.rsplit('/', 1)[-1]}")
    header, time, sumtot, syst_err = read_cbd(url)
    headers.append(header)
    times.append(time)
    counts.append(sumtot)
    syst_errors.append(syst_err)

header = headers[0]
timedel = header["TIMEDEL"]

# All four detectors share identical bin edges, so the light curves can simply
# be summed; check rather than assume.
for other in times[1:]:
    assert np.array_equal(times[0], other), "detectors do not share bin edges"

# OGIP: MJDREF is expressed in the scale named by TIMESYS. This is the
# epoch as the header states it -- the defective one, 37 s later than the
# true MET epoch (2021-01-01 00:00:00 TAI); see the module docstring.
epoch = Time(header["MJDREFI"] + header["MJDREFF"], format="mjd",
             scale=header["TIMESYS"].lower())
t0_met = (Time(T0_UTC, scale="utc") - epoch).sec
print(f"\nepoch  = {epoch.utc.isot} UTC (from MJDREFI/MJDREFF/TIMESYS)")
print(f"t0     = {T0_UTC} UTC = MET {t0_met:.3f}")

# TIMEPIXR=1, so TIME is the END of each bin; plot bin centres.
centre = times[0] - 0.5 * timedel
rel = centre - t0_met
rate = np.sum(counts, axis=0) / timedel
syst_err = syst_errors[0]
for other in syst_errors[1:]:
    assert np.array_equal(syst_err, other), "TIME_SYST_ERROR differs by detector"

window = np.abs(rel) <= HALF_WIDTH
rel, rate, syst_err = rel[window], rate[window], syst_err[window]
print(f"{rel.size} bins within +/-{HALF_WIDTH:.0f} s of t0")
print(f"TIME_SYST_ERROR over this window: {syst_err.min():.3f} to "
      f"{syst_err.max():.3f} s")

trim = (rel.size // REBIN) * REBIN
rel_b = rel[:trim].reshape(-1, REBIN).mean(axis=1)
rate_b = rate[:trim].reshape(-1, REBIN).mean(axis=1)

fig, (ax_lc, ax_err) = plt.subplots(
    2, 1, figsize=(10, 6.5), sharex=True, height_ratios=(3, 1))

ax_lc.step(rel, rate, where="mid", lw=0.6, color="0.75",
           label=f"{timedel:g} s (native)")
ax_lc.step(rel_b, rate_b, where="mid", lw=1.6, color="#1f5fa9",
           label=f"{timedel * REBIN:g} s")
ax_lc.set_ylabel("Rate, CS0-CS3 summed (ct/s)")
ax_lc.set_title("BurstCube CBD around GRB 240629A "
                "(t0 = Fermi GBM trigger bn240629704)")
ax_lc.legend(loc="upper right", fontsize=9, frameon=False)

ax_err.step(rel, syst_err, where="mid", lw=1.6, color="#1f5fa9")
ax_err.set_ylabel("TIME_SYST_ERROR (s)")
ax_err.set_xlabel(f"MET - t0 (s)          [t0 = MET {t0_met:.3f}]")
ax_err.set_ylim(0, max(1.0, syst_err.max() * 1.6))

for ax in (ax_lc, ax_err):
    ax.axvline(0.0, color="#b2182b", lw=1.2)
    ax.axvline(OBSERVED_OFFSET, color="#b2182b", lw=1.2, ls="--")
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_xlim(-HALF_WIDTH, HALF_WIDTH)

ax_lc.annotate("GBM t0", xy=(0, 1), xytext=(-4, -10),
               textcoords="offset points", xycoords=("data", "axes fraction"),
               ha="right", va="top", color="#b2182b", fontsize=9)
ax_lc.annotate(f"burst seen here, t0{OBSERVED_OFFSET:+.1f} s",
               xy=(OBSERVED_OFFSET, 1), xytext=(6, -10),
               textcoords="offset points", xycoords=("data", "axes fraction"),
               ha="left", va="top", color="#b2182b", fontsize=9)

fig.tight_layout()
outfile = sys.argv[1] if len(sys.argv) > 1 else "grb240629a_lightcurve.png"
fig.savefig(outfile, dpi=130)
print(f"\nwrote {outfile}")
