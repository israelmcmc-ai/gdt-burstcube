.. _gdt-burstcube:

********************************************************
Welcome to BurstCube Gamma-ray Data Tools Documentation!
********************************************************

The BurstCube Gamma-ray Data Tools (GDT) is a toolkit for BurstCube data
built on the :external:ref:`GDT Core Package<gdt-core>`.

BurstCube is a CubeSat gamma-ray burst mission with four scintillator
detectors, observing the unocculted sky. This toolkit provides readers for
every public BurstCube data product -- continuous binned data (CBD),
time-tagged events (TTE), the detector response grid, orbit/attitude,
housekeeping, the mission timeline, and the observation catalog -- plus
finders that download products directly from the HEASARC archive.

.. rubric:: Caveats

BurstCube is a young mission with real, documented rough edges in its public
data: broken header keywords in some files, mid-mission changes to detector
thresholds, an EBOUNDS table in the response files that disagrees with
CALDB, attitude reconstructed for only 3 brief epochs across the whole
mission, and frequent data gaps. This toolkit works around what it safely
can and surfaces the rest rather than hiding it. Read the **Caveats** section
of the top-level `README <https://github.com/USRA-STI/gdt-burstcube#caveats>`_
before drawing conclusions from BurstCube data.

.. rubric:: Additional Resources

For questions, bug reports, and comments, please visit the
`GDT-BurstCube GitHub repository <https://github.com/USRA-STI/gdt-burstcube>`_.

***************
Getting Started
***************
.. toctree::
   :maxdepth: 1

   install
   notebooks

******************
User Documentation
******************

BurstCube Definitions
=====================
.. toctree::
   :maxdepth: 1

   missions/burstcube/time
   missions/burstcube/detectors
   missions/burstcube/frame
   missions/burstcube/caldb
   missions/burstcube/headers

Data Types
==========
.. toctree::
   :maxdepth: 1

   missions/burstcube/cbd
   missions/burstcube/tte
   missions/burstcube/response

Ancillary Data
==============
.. toctree::
   :maxdepth: 1

   missions/burstcube/orbit
   missions/burstcube/attitude
   missions/burstcube/gti
   missions/burstcube/saa
   missions/burstcube/hk
   missions/burstcube/timeline

Data Finders and Catalogs
==========================
.. toctree::
   :maxdepth: 1

   missions/burstcube/finders
   missions/burstcube/catalog


Indices and tables
===================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
