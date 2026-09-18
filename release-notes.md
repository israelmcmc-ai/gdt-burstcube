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
- Corrected the CALDB SAA polygon's column mapping: the region file's `X` is
  latitude and `Y` is longitude, the opposite of its own `TTYPE` comments.
  `BurstCubeSaa` also closes the polygon, which the file leaves open. Both
  are cross-checked against Fermi GBM's `GbmSaaPolygon5`, whose first 11
  vertices BurstCube's polygon reproduces exactly.
- Renamed `BurstCubeGti` to `BurstCubeGTI` and moved `intersect`, `union`,
  `complement` and `apply_to` onto it as class methods.
- Added `caldb.regroup_edges()`, composing a coarse-to-fine channel regroup
  from CALDB's own rebinning tables instead of a hardcoded index list.
- Added `BurstCubeObsFinder.obs_id_from()` and an `obs_id` property.
- Documented that `DETECTOR_HK2`'s `BASE_THRES`, not `PEAK_THRES`, carries
  caveat #7's mid-mission energy threshold change, and that the energy it
  corresponds to in CALDB's scale is ~20% below caveat #7's table.
- The CALDB SAA polygon also lists two pairs of vertices out of order, so the
  boundary crosses itself twice. `BurstCubeSaa` sorts the vertices by angle
  about their centroid, but only when the file's own order self-intersects
  and sorting fixes it.
- Removed `response.NSIDE`/`NUM_PIXELS`. The grid's resolution now comes from
  `caldb.response_grid()`, which derives it from the CALDB simulation file's
  row count; call sites read `.nside`/`.num_pixels` off it directly. That
  file is bundled, so the grid stays available offline (bundled CALDB: 10
  files, 233 kB -> 11 files, 582 kB).
- Documented what `DETECTOR_HK2`'s `PEAK_THRES` actually is -- a pulse-shape
  cut, not an energy threshold -- from the file's own column comments, since
  nothing else in the archive documents either column.
