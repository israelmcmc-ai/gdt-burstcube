.. _notebooks:


Jupyter Notebook Tutorials
==========================

These notebooks were executed top-to-bottom against the live HEASARC archive
(observation day 240530, detector CS0 as the running example unless noted),
and are committed here with their real outputs -- nothing in them is
hand-written. Where something did not work as a naive port from another
mission might expect, the notebook shows the actual error and the real
workaround, rather than hiding it.

.. toctree::
   :maxdepth: 1

   notebooks/1_data_types_and_binning
   notebooks/2_response
   notebooks/3_ancillary
   notebooks/4_grb_240629a

----

**1. Data Types and Binning** -- CBD (continuous binned data) and TTE
(time-tagged events), the difference between ``_uf`` (unfiltered) and ``_cl``
(cleaned) variants, finding and downloading a full observation day, time
selection, rebinning, light curves, and count spectra -- including the
caveats that bite in practice: bin-end (``TIMEPIXR=1``) timestamps, the
``TIME_SYST_ERROR`` column, the mid-mission energy threshold change, and a
real rebinning failure caused by genuine 1 ms-short bins in the archive data.
It also compares TTE against CBD on day 240814 -- binned on CBD's own bin
edges, in time and in energy, plus the same events at TTE's native 1024
channels -- which is where the TTE recording gaps show up.

**2. Detector Response** -- the HEALPix response grid (nside=16, spacecraft
coordinates), downloading only the pixels needed for one direction,
interpolating a DRM, plotting the DRM and effective area, folding a spectrum,
regrouping to CBD's 16 channels, and comparing the four detectors at one
location.

**3. Ancillary Data** -- orbit/ephemeris and Earth position, the 3
reconstructed attitude epochs in the whole archive, supplying your own
attitude quaternion, GTI and SAA filtering, detector housekeeping (enable
flags and energy thresholds), the mission timeline (and its known 37 s UTC
offset), and the ``burcbmastr`` observation catalog.

**4. GRB 240629A: A Joint BurstCube/GBM Worked Example** -- a complete
worked analysis of one real burst, chaining the previous 3 notebooks'
pieces together the way a real analysis would: getting the trigger time,
T90, and sky position from GBM's own Trigger and Burst Catalogs (via
``astro-gdt-fermi``); downloading and plotting the corresponding BurstCube
light curve, individually and summed across detectors, against GBM's own
(binning GBM's unbinned TTE onto BurstCube's bin edges, since GBM's CTIME
uses adaptive rather than fixed time binning); confirming which of
BurstCube's 3 reconstructed attitude epochs actually covers the trigger and
using it to get the response; a background fit excluding GBM's T90 plus a
10 s buffer; and a joint 4-detector power-law spectral fit. The result is a
clean non-detection and a 90% upper limit -- not a burst BurstCube saw, but
a fully worked, honestly reported analysis of the case where it did not.

----

Standalone Example Scripts
--------------------------

``examples/`` holds short, self-contained scripts that each download what
they need and reproduce one specific result:

* ``examples/tte_gaps_vs_cbd.py`` -- plots TTE's recording gaps against the
  CBD rate over the same interval, the minimal reproduction of the README
  caveat *TTE gaps*.
