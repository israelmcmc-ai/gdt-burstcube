# Release Notes for Gamma-ray Data Tools: BurstCube

## Version 0.1.0 (Unreleased)

Initial package foundation.

- Added `BurstCubeSecTime` (MET, continuous TT seconds since 2021-01-01) and
  `BurstCubeObsId` (`YYMMDD` observation-day) time formats.
- Added `BurstCubeDetectors`, with azimuth and zenith derived numerically from the
  bundled CALDB alignment matrices for CS0-CS3.
- Added `BurstCubeFrame`, including construction from a user-supplied quaternion for
  the many mission epochs with no reconstructed attitude.
- Added the CALDB resolution chain (`$CALDB` -> local cache -> download -> bundled
  package data) and accessors for the bundled alignment, ebounds, rebin, and SAA
  region files.
- Added FITS header definitions for the CBD, TTE, orbit, attitude, and detector
  housekeeping products.
- Bundled the 10 BurstCube CALDB files needed for offline operation (233,280 bytes).
- Added `BurstCubeTTE.recording_blocks()`, which returns the intervals a TTE
  file actually recorded over as a `Gti`. TTE covers its nominal span in short
  blocks separated by gaps of comparable length, and reads zero across a gap
  while the detector keeps counting, so any rate taken over an interval that
  spans a gap is diluted. See the README's *TTE gaps* caveat and
  `examples/tte_gaps_vs_cbd.py`.
- Licensed under MIT, matching the BurstCube team's `bctools`. Per-file license
  headers are gone; `LICENSE` and `NOTICE` carry the terms and the attribution
  owed to the GDT Core authors.
