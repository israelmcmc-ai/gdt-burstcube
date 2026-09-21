#!/usr/bin/env python
"""Minimal reproduction of the README caveat *Timeline UTC column is 37
seconds off its own MET column*.

The BurstCube MET epoch is defined by ``MJDREFI``/``MJDREFF``/``TIMESYS`` in
every science file's header. Read per the OGIP convention -- ``MJDREF`` is
expressed in the scale named by ``TIMESYS`` -- those keywords place the epoch
at 2021-01-01 00:00:00 UTC, written as 2021-01-01 00:01:09.184 TT.

Applying that epoch to a CBD file's own ``TSTART``/``TSTOP`` reproduces its
``DATE-OBS``/``DATE-END`` strings exactly. Applying the *same* epoch to the
trend timeline CSV's MET column does not reproduce that file's own UTC
column: every row is off by exactly 37.000 s, which is TAI - UTC since
2017-01-01.

Note what this does and does not show. It establishes that two archive
products disagree by exactly a leap-second offset. It does **not** establish
which side is correct -- per the archive's own caveat #11 the MET in the
science files is a reconstructed quantity with documented errors of tens of
seconds, so neither 0 nor 37.000 s is necessarily the truth.

Unlike the other scripts here this one deliberately does **not** use
``gdt-burstcube``, or any GDT package -- it needs only ``astropy`` and the
standard library, so it can be run by anyone holding the archive files
without installing this plugin first. Nothing is hardcoded: the epoch comes
from the downloaded file's own header keywords.

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

# OGIP: MJDREF is expressed in the scale named by TIMESYS.
epoch = Time(mjdrefi + mjdreff, format="mjd", scale=timesys.lower())

print(f"MJDREFI = {mjdrefi}")
print(f"MJDREFF = {mjdreff!r}  ({mjdreff * 86400:.3f} s)")
print(f"TIMESYS = {timesys!r}")
print(f"epoch   = {epoch.utc.isot} UTC  (= {epoch.tt.isot} TT)")


def met_to_utc(met):
    return (epoch + TimeDelta(met, format="sec")).utc


# ------------------------------------------- 1. the FITS file agrees with itself
print("\n1. CBD file, MET vs its own DATE-OBS/DATE-END")
for met_key, date_key in (("TSTART", "DATE-OBS"), ("TSTOP", "DATE-END")):
    met = header[met_key]
    written = Time(header[date_key], scale="utc")
    derived = met_to_utc(met)
    print(f"   {met_key} = {met:.3f}")
    print(f"      epoch + MET      : {derived.isot}")
    print(f"      {date_key} in file : {written.isot}")
    print(f"      difference       : {(derived - written).sec:+.3f} s")

# ---------------------------------------------- 2. the timeline CSV does not
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
