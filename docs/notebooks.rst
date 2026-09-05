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

----

**1. Data Types and Binning** -- CBD (continuous binned data) and TTE
(time-tagged events), the difference between ``_uf`` (unfiltered) and ``_cl``
(cleaned) variants, finding and downloading a full observation day, time
selection, rebinning, light curves, and count spectra -- including the
caveats that bite in practice: bin-end (``TIMEPIXR=1``) timestamps, the
``TIME_SYST_ERROR`` column, the mid-mission energy threshold change, and a
real rebinning failure caused by genuine 1 ms-short bins in the archive data.

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
