#!/usr/bin/env python
"""Minimal reproduction of the README caveat *The archive's MET epoch is
stated wrong by 37 seconds*.

The BurstCube MET epoch is defined by ``MJDREFI``/``MJDREFF``/``TIMESYS`` in
every science file's header. Read per the OGIP convention -- ``MJDREF`` is
expressed in the scale named by ``TIMESYS`` -- those keywords place the epoch
at 2021-01-01 00:00:00 UTC, written as 2021-01-01 00:01:09.184 TT. **That is
the archive's defect**, not the truth: the true MET epoch is 2021-01-01
00:00:00 TAI, 37.000 s earlier (see the README Caveats section and
``gdt.missions.burstcube.time`` for the independent GRB 240629A evidence that
settles which side is correct).

Applying the header's stated (defective) epoch to a CBD file's own
``TSTART``/``TSTOP`` reproduces its ``DATE-OBS``/``DATE-END`` strings exactly
-- unsurprising, since those strings were generated from the same defective
epoch. Applying that *same* stated epoch to the trend timeline CSV's MET
column does *not* reproduce that file's own UTC column: every row is off by
exactly 37.000 s, which is TAI - UTC since 2017-01-01. The timeline CSV was
the one telling the truth all along; it is the FITS headers' stated epoch
that is wrong by that same 37 s.

Unlike the other scripts here this one deliberately does **not** use
``gdt-burstcube``, or any GDT package -- it needs only ``astropy`` and the
standard library, so it can be run by anyone holding the archive files
without installing this plugin first. Nothing is hardcoded: the epoch comes
from the downloaded file's own header keywords, taken at face value -- which
is exactly what makes it demonstrate the defect rather than the corrected
epoch this package's own ``BurstCubeSecTime`` actually uses.

Usage::

    python examples/timeline_utc_vs_met.py
"""
import gzip
import io
import urllib.request

from astropy.io import fits
from astropy.time import Time, TimeDelta

BASE = "https://heasarc.gsfc.nasa.gov/FTP/burstcube/data"
CBD = f"{BASE}/obs/2024_06/240629/monitor/bc240629cs0_3cbd_cl.fits.gz"
TIMELINE = f"{BASE}/trend/timeline/20250606_timeline_final.csv"


def fetch(url):
    with urllib.request.urlopen(url) as response:
        return response.read()


# ---------------------------------------------------------------- the epoch
with fits.open(io.BytesIO(gzip.decompress(fetch(CBD)))) as hdul:
    header = hdul["CBD"].header

mjdrefi = header["MJDREFI"]
mjdreff = header["MJDREFF"]
timesys = header["TIMESYS"]

# OGIP: MJDREF is expressed in the scale named by TIMESYS. This is the
# epoch the header states -- the defective one, 37 s later than the true
# MET epoch (2021-01-01 00:00:00 TAI); see the module docstring.
epoch = Time(mjdrefi + mjdreff, format="mjd", scale=timesys.lower())

print(f"MJDREFI = {mjdrefi}")
print(f"MJDREFF = {mjdreff!r}  ({mjdreff * 86400:.3f} s)")
print(f"TIMESYS = {timesys!r}")
print(f"epoch   = {epoch.utc.isot} UTC  (= {epoch.tt.isot} TT)")


def met_to_utc(met):
    return (epoch + TimeDelta(met, format="sec")).utc


# --------------------------------- 1. the FITS file agrees with itself
# (both DATE-OBS/DATE-END and this epoch were generated from the same
# defective epoch, so of course they agree)
print("\n1. CBD file, MET vs its own DATE-OBS/DATE-END")
for met_key, date_key in (("TSTART", "DATE-OBS"), ("TSTOP", "DATE-END")):
    met = header[met_key]
    written = Time(header[date_key], scale="utc")
    derived = met_to_utc(met)
    print(f"   {met_key} = {met:.3f}")
    print(f"      epoch + MET      : {derived.isot}")
    print(f"      {date_key} in file : {written.isot}")
    print(f"      difference       : {(derived - written).sec:+.3f} s")

# ---------------------------- 2. the timeline CSV does not (it is correct;
# the header's stated epoch, used above, is the defective one)
print("\n2. Trend timeline CSV, MET vs its own UTC column")
rows = fetch(TIMELINE).decode().splitlines()
for line in rows[:3]:
    met_str, utc_str, description = line.split(",", 2)
    derived = met_to_utc(float(met_str))
    written = Time(utc_str, scale="utc")
    print(f"   {description.strip()}")
    print(f"      MET              : {float(met_str):.3f}")
    print(f"      epoch + MET      : {derived.isot}")
    print(f"      UTC col in file  : {written.isot}")
    print(f"      difference       : {(derived - written).sec:+.3f} s")

offsets = []
for line in rows:
    if not line.strip():
        continue
    met_str, utc_str, _ = line.split(",", 2)
    offsets.append((met_to_utc(float(met_str)) - Time(utc_str, scale="utc")).sec)

print(f"\n   all {len(offsets)} rows: min {min(offsets):+.3f} s, "
      f"max {max(offsets):+.3f} s")
print("\n37.000 s is exactly TAI - UTC since 2017-01-01.")
