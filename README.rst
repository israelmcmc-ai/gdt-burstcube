================
GDT-BurstCube
================

GDT-BurstCube is an extension to the Gamma-ray Data Tools (GDT) that adds
readers, finders, and a detector response interface for the BurstCube
mission: continuous binned data (CBD), time-tagged events (TTE), the
detector response grid, orbit/attitude, detector housekeeping, the mission
timeline, and the ``burcbmastr`` observation catalog.

Normal Installation
--------------------

If you don't plan to contribute code to the project, the recommended install
method is installing from PyPI using:

.. code-block:: sh

   pip install astro-gdt-burstcube
   gdt-data init

The ``gdt-data init`` is required to initialize the library after installation
of astro-gdt. You do not need to perform the initialization again if
astro-gdt was already installed and initialized. There is no harm in running
it again "just in case".

Quickstart
----------

.. code-block:: python

   from gdt.missions.burstcube.finders import BurstCubeObsFinder
   from gdt.missions.burstcube.cbd import BurstCubeCbd

   # find and download the cleaned CBD file for one detector, one day
   finder = BurstCubeObsFinder('240530')
   paths = finder.get_cbd('./data', detectors='CS0', variant='cl')

   cbd = BurstCubeCbd.open(paths[0])
   lc = cbd.to_lightcurve()

See the Jupyter notebooks (inside docs/notebooks)
for three complete, executed worked examples against the live archive: data
types and binning, detector response, and ancillary data (orbit, attitude,
GTI/SAA, housekeeping, timeline, catalog).

**Read the Caveats section below before drawing conclusions from BurstCube
data** -- several of the mission's real quirks are not obvious from the API
alone.


================
Caveats
================

BurstCube public archive has some
rough edges. This toolkit works around what it safely can and documents the
rest here. This section covers the issues most likely
to affect analysis. It is not exhaustive, see the official caveats list
`BurstCube Archive Caveats <https://heasarc.gsfc.nasa.gov/docs/burstcube/archive/burstcube_archive_caveats.pdf>`_
document (2025-07-22).

Response EBOUNDS mismatch
--------------------------

A ``.rsp`` file's own internal ``EBOUNDS`` extension is an **older, rounded,
detector-independent** energy grid -- it is not the same grid CALDB uses for
CBD's 16-channel data, and not the same for every detector. For example, on
CS0:

.. code-block::

   .rsp EBOUNDS  top channel:  1290.00 -- 5000.00 keV   (same for all 4 detectors)
   CALDB eb16    top channel:  1298.79 -- 1648.51 keV   (CS0-specific)

The channel **grouping** (which of the native 64 channels combine into each
of CBD's 16) is identical between the two grids -- only the energies differ.
Because of this, ``BurstCubeRsp.to_cbd()``
regroups channels **by index** and then replaces the channel energies with
CALDB's per-detector ``eb16``, discarding the ``.rsp``'s own EBOUNDS
energies entirely. If you fold a spectrum through a 16-channel response and
report channel energies, make sure you are reporting CALDB's, not the
``.rsp`` file's own.

Attitude is unavailable for most of the mission
--------------------------------------------------

Real, reconstructed attitude information exists for only **3** brief epochs
across the entire archive (``trend/attitude/bc_csa_att.fits``), each with an
approximate 10-12 degree 1-sigma pointing error. The detector response grid
is defined in spacecraft coordinates, so **without attitude you cannot place
a source on the sky** -- there is no way to convert a sky position into a
response pixel, or vice versa, for the other ~86 days.

The documented, supported path for every other time is to supply your own
independently reconstructed attitude quaternion via
``BurstCubeFrame.from_quaternion()``. Asking
for a sky-position response without one raises a clear error rather than
silently assuming an orientation.

TTE header keywords are broken in at least one real file
-------------------------------------------------------------

at has an ``EVENTS`` extension whose own
``TSTART``/``TSTOP`` header keywords cannot be trusted: they are stored as
**strings** rather than numbers, ``TSTOP`` is less than ``TSTART``, and the
derived ``EXPOSURE`` is negative. ``BurstCubeTte`` ignores these keywords
entirely and derives the time range from the ``STDGTI`` extension and the
event times themselves, emitting a ``UserWarning`` when it does so. TTE
exists for only 7 of the roughly 89 observation days in the archive, so this
mainly matters if you read TTE headers directly rather than through this
reader.

Timeline UTC column is 37 seconds off its own MET column
-------------------------------------------------------------

``trend/timeline/*_timeline_final.csv`` has three columns and no header row:
a MET, a human-readable UTC string, and an event description.

.. code-block::

   106367555.28030825,2024-05-16T02:31:58.280,Spacecraft Reboot

**The UTC string is 37 seconds earlier than the MET on the same row.** The
MET is correct; the UTC string is not.

*Why there should be no offset at all.* BurstCube MET is defined by
``MJDREFI=59215``, ``MJDREFF=0.00080074074074074``, ``TIMESYS='TT'``. That
fraction is 69.184 s = 32.184 (TT-TAI) + 37 (TAI-UTC in 2021), so the epoch
is 2021-01-01 00:00:00 UTC written in TT, and MET counts TT seconds from it.
TT-UTC held constant at 69.184 s from 2017 through the whole mission -- no
leap second occurred -- so the conversion collapses to a plain addition with
no correction term:

.. code-block::

   UTC = 2021-01-01T00:00:00 + MET

*What the file actually contains.* All 505 rows reproduce exactly --
string-identical after rounding to milliseconds -- under:

.. code-block::

   UTC_csv = 2021-01-01T00:00:00 + MET - 37 s

That is the signature of treating the epoch as 2021-01-01 00:00:00 **TAI**
and the MET count as TAI seconds, then converting TAI to UTC by subtracting
the 37 leap seconds. The conversion step is right; the premise about what
MET counts is wrong. It is the same 37-second confusion that caveat #11 of
the official document describes for ground-commanded time jams.

The residual scatter after that model is at most half a millisecond, which
is entirely the CSV's own rounding to three decimals -- the underlying
offset is a clean constant, not a drift.

*Which one is authoritative.* The MET is. The epoch above reproduces
``DATE-OBS`` and ``DATE-END`` to the millisecond in the CBD (both ``_uf``
and ``_cl``) and orbit files, and the ``burcbmastr`` catalog's own UTC
strings agree with it too. The FITS products are self-consistent; *this CSV
is the outlier*.

``BurstCubeTimeline`` therefore parses the MET column and converts it with
``BurstCubeSecTime``. The CSV's own string is exposed unconverted as
``utc_as_written`` for provenance only -- **do not use it for analysis.**

----

For the full list of caveats -- including unphysical detector counts, CBD
sub-second clock drifts after PPS dropouts, CBD data blending below the
nominal 0.256 s cadence, and the detailed time-correction bookkeeping behind
``ORIGINAL_TIME``/``TIME``/``TIME_SYST_ERROR`` -- see the official
`BurstCube Archive Caveats document`__.

__ https://heasarc.gsfc.nasa.gov/docs/burstcube/archive/burstcube_archive_caveats.pdf


Setting up a Development Environment
--------------------------------------

If you do want to contribute code to this project (and astro-gdt), you can
use the following commands to quickly set up a development environment:

.. code-block:: sh

   mkdir gdt-devel
   cd gdt-devel
   python -m venv venv
   . venv/bin/activate
   pip install --upgrade pip setuptools wheel
   git clone git@github.com:USRA-STI/gdt-core.git
   git clone git@github.com:USRA-STI/gdt-burstcube.git
   pip install -e gdt-core/
   gdt-data init
   pip install -e gdt-burstcube/

This should result in git-devel having the following directory structure::

   .
   ├── venv
   ├── gdt-core
   └── gdt-burstcube

and both gdt-core and gdt-burstcube installed in the virtual environment named
venv.


Helping with Documentation
-----------------------------

You can contribute additions and changes to the documentation. In order to
use sphinx to compile the documentation source files, we recommend that you
install the packages contained within ``docs/requirements.txt``.

To compile the documentation, use the following commands:

.. code-block:: sh

   cd gdt-burstcube/docs
   pip install -r requirements.txt
   make html
