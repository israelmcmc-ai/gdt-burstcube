# gdt-burstcube

A Gamma-ray Data Tools (GDT) plugin for BurstCube, built on `astro-gdt`
(gdt-core 2.2.3) and modelled on `gdt-fermi`. Readers and finders for every
public BurstCube data product, plus the HEALPix detector response grid.

Package: `astro-gdt-burstcube`, namespace `gdt.missions.burstcube`. MIT
licensed (matching the BurstCube team's `bctools`); `NOTICE` carries the
attribution owed to the GDT Core authors and is kept regardless of license.

## Layout

```
src/gdt/missions/burstcube/    the plugin
    time.py detectors.py frame.py caldb.py headers.py   definitions
    cbd.py tte.py response.py                           science data
    orbit.py attitude.py gti.py saa.py hk.py timeline.py  ancillary
    finders.py catalog.py                               archive access
    data/                     11 bundled CALDB files (582 kB)
src/gdt/data/burstcube.urls   20 real archive files, for gdt-data download
tests/missions/burstcube/
    offline/                  157 tests, synthetic fixtures + bundled CALDB
    data_driven/              38 tests, real archive files, skip if absent
docs/notebooks/               4 tutorials, committed with real outputs
examples/                     3 standalone scripts: TTE gaps (uses the
                              plugin); timeline UTC vs MET and the GRB
                              240629A light curve (no GDT package at all)
```

## Commands

```sh
python -m pytest tests -q                              # all 195
python -m pytest tests/missions/burstcube/offline -q   # no network needed
gdt-data download burstcube                            # fetch data_driven files
python -m build --wheel                                # check package data ships
```

Notebooks are committed **executed, with outputs**. Re-run any you touch:

```sh
cd docs/notebooks && jupyter nbconvert --execute --inplace \
  --ExecutePreprocessor.kernel_name=burstcube-kernel <nb>.ipynb
```

Then verify zero `output_type == 'error'` cells and zero null
`execution_count`s before committing — `nbconvert` exits 0 on some failures.

CI (`.github/workflows/ci.yml`, Python 3.11/3.12) runs **only the offline
tier**: no network, no `gdt-data` step. The data-driven tier is for local use
and skips cleanly when the archive files are missing.

## Conventions

- **Acronyms are fully capitalised in class names**: `BurstCubeCBD`,
  `BurstCubeTTE`, `BurstCubeHK`, `BurstCubeGTI`, `BurstCubeSAA`, `SAARegion`,
  `GTIHeaders`, `CBDGTIHeader`, `TTEGTIHeader`, `GTITrend*Header`,
  `DetectorHK*`. Deliberate exceptions, left as-is by request: the `Rsp*`
  classes and `BurstCubeObsId`. `Obs`, `Aux`, `Det`, `Sec` are abbreviations,
  not acronyms, and stay mixed-case.
- **Nothing derivable from CALDB is hardcoded.** `caldb.py` is the single
  source for ebounds, rebinning, alignment, the SAA polygon and the response
  grid's `nside`. Call `caldb.response_grid().nside` at the point of use
  rather than caching it in a module constant (~3 ms per call, measured).
- **Readability over optimisation** unless the difference exceeds ~2x.
- **Tutorials and examples use cleaned (`_cl`) CBD only.** The two exceptions
  are the sections whose subject *is* the `_uf`/`_cl` difference (notebook 1)
  and `BurstCubeGTI.apply_to`, which reconstructs a cut from `_uf`.
- Archive quirks go in the README's **Caveats** section, with the measurement
  behind them; tutorials point at it rather than re-explaining.

## The lesson that cost the most

**Synthetic fixtures agreed with the code because both came from the same
assumptions.** Every one of the seven implementation bugs found in this repo
was found by running against real archive files, not by the offline suite.
The offline tier is for fast regression cover; it is not evidence that a
reader works. Before believing a new reader, run it over the real files in
`gdt-data`'s download directory, and prefer a `data_driven` test that pins
what the real file actually contains.

Corollary that also cost time: verify claims about the data before writing
them down. Several confident conclusions in this repo's history were wrong
(a "100x TTE/CBD disagreement" that was a duty-cycle artifact, a "saturation
at 508 ct/s" that was one day's coincidence, a `np.isclose` result ruined by
the default `rtol` at MET ~1e8). Measure, state the method, and say what the
data does not settle.

## BurstCube data facts worth not re-deriving

- **MET**: `MJDREFI=59215`, `MJDREFF=0.00080074074074074` (= 69.184 s =
  32.184 TT−TAI + 37 TAI−UTC), `TIMESYS='TT'`. Epoch is 2021-01-01 00:00:00
  UTC written in TT; MET counts TT seconds with no leap-second bookkeeping.
- **CBD**: `TIMEPIXR=1`, so `TIME` is the **end** of each bin; edges are
  `[TIME − TIMEDEL, TIME]`. Nominal 0.256 s. Float noise is ~1 ULP (1.5e-8 s)
  while genuine short bins are exactly 1.0e-3 s short — `_SNAP_TOL = 1e-6` in
  `cbd.py` sits between them deliberately. `SUMTOT`/`SUMCH_2_15` exist only
  in `_cl` files.
- **TTE**: 28 files total (7 days × 4 detectors), 1024 channels. Event times
  are quantised to 10 µs, so many distinct photons share a timestamp — hence
  `BurstCubeTTE.slice_time` overrides the base class, whose
  `EventList.merge(force_unique=True)` would discard ~74% of the counts.
- **Channel schemes**: CALDB describes `reb16` and `reb64` against the 1024
  native channels, never against each other; `caldb.regroup_edges` composes
  them. The 64→16 result is 17 edge indices — a 16-entry list silently gives
  15 channels while preserving the folded total.
- **Response grid**: HEALPix nside=16 RING, 3072 pixels, in spacecraft
  coordinates. Defined by `cpf/simloc/bccsa_sim_20221001v001.fits` (one row
  per pixel); the `.rsp` files record only their own `PIXEL` and `ORDERING`.
  `get_drm` interpolates over the 4 neighbouring pixels, matching
  `bctools.InstrumentResponse.get_drm`.
- **Archive layout is deterministic** — filenames follow from day, detector
  and product type — so the finders build paths directly and never list a
  remote directory. A missing file is a per-file 404, warned and skipped.

## Archive defects this plugin works around

Full write-ups with evidence are in README §Caveats. In brief:

| Defect | Where handled |
|---|---|
| TTE `TSTART`/`TSTOP` unusable in all 28 files (strings; 16 have `TSTOP < TSTART`) | `tte.py` uses `STDGTI` + event times, warns |
| TTE records in short blocks with gaps; reads zero while the detector counts | `BurstCubeTTE.recording_blocks()` |
| TTE also drops events inside bins above ~200 ct/s (agreement is 0.999 below that) | documented only |
| CALDB SAA polygon: `X` is latitude and `Y` longitude, contradicting its own `TTYPE` comments | `caldb.saa_region` |
| Same polygon: two vertex pairs transposed, so it self-intersects; and left open | `BurstCubeSAA` |
| CALDB `eb1024` and caveat #7's keV table differ by a constant ~0.8 gain factor | documented only |
| `PEAK_THRES` is a pulse-shape cut, not the energy threshold (`BASE_THRES` is) | docstrings, `hk.rst` |
| `.rsp` internal `EBOUNDS` is an older, detector-independent grid | `to_cbd()` substitutes CALDB `eb16` |
| Timeline CSV's UTC column is 37 s off its own MET column | `timeline.py` exposes it as `utc_as_written` only |
| Attitude exists for 3 epochs mission-wide, ~10-12° errors | no interpolation; `BurstCubeFrame.from_quaternion` |

Three of these are defects in the archive/CALDB files themselves, not in this
package: the SAA column swap, the SAA vertex order, and the `eb1024` scale.
If they are ever fixed upstream, the corresponding tests fail loudly, which
is intentional.

## gdt-core issues encountered

Tracked in [#3](https://github.com/israelmcmc-ai/gdt-burstcube/issues/3),
[#4](https://github.com/israelmcmc-ai/gdt-burstcube/issues/4) and
[#5](https://github.com/israelmcmc-ai/gdt-burstcube/issues/5), all labelled
`dependencies`:

- `ChannelEffectiveArea` plots `photon_effective_area()` under a "Channel
  Energy" axis; `channel_effective_area()` is an unweighted column sum, so
  not in cm².
- The DRM plot classes drop their artists if the returned object is
  garbage-collected before the figure is drawn — assign them to a variable.
- `Phaii.slice_time` raises on a real GTI two ways (empty segments give
  `time_range is None`; disjoint segments empty its cumulative accumulator).
  `BurstCubeGTI.apply_to` works around both; the comments in its body explain
  how.
- `Phaii.to_pha()` hardcodes the GBM `'SPECTRUM'` FITS extension name when
  copying header keywords, so it raises `KeyError` on any mission whose data
  extension is named something else. `BurstCubeCBD.to_pha()` overrides it,
  pointed at `'CBD'`.

## Workflow

Development happens on `claude/burstcube-gdt-plugin-ops82w`, reviewed in
[PR #1](https://github.com/israelmcmc-ai/gdt-burstcube/pull/1). The reviewer
pushes direct edits to the same branch, so **fetch and merge before
starting** — and preserve their wording when re-touching a file they edited.
Reply to review threads individually, with the commit SHA and what actually
changed; every GitHub comment ends with the Claude Code attribution footer.
