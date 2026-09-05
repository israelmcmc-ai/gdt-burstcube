================
GDT-BurstCube
================

GDT-BurstCube is an extension to the Gamma-ray Data Tools (GDT) that adds
readers, finders, and a detector response interface for the BurstCube
mission: continuous binned data (CBD), time-tagged events (TTE), the
detector response grid, orbit/attitude, detector housekeeping, the mission
timeline, and the ``burcbmastr`` observation catalog.

This software is not subject to EAR.

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

Caveats
--------

BurstCube's public archive has real, known rough edges -- an EBOUNDS grid in
each response file that does not match CALDB's per-detector energies, CBD
timestamps that mark the end of each bin rather than the start, attitude
reconstructed for only 3 brief epochs across the whole mission, broken
header keywords in at least one real TTE file, a mission timeline whose own
UTC column is 37 seconds off its own MET column, and frequent data gaps.
This toolkit documents each of these rather than hiding them. See the full
**Caveats** section of the
`README on GitHub <https://github.com/USRA-STI/gdt-burstcube#caveats>`_
before drawing conclusions from BurstCube data.

Contributing Code or Documentation
--------------------------------------

If you plan to help with the development or documentation of astro-gdt, then
please visit our github site at https://github.com/USRA-STI/gdt-burstcube.
