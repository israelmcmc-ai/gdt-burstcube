.. _burstcube-caldb:

*************************************************************
BurstCube CALDB Access (:mod:`gdt.missions.burstcube.caldb`)
*************************************************************

A request for a CALDB file (an EBOUNDS table, a rebinning table, an
alignment file, the SAA region) is resolved through a chain of four steps,
tried in order: an ``$CALDB`` environment variable using the standard HEASoft
layout, a local cache directory, a fresh download from HEASARC, and finally
the 10 files bundled directly with this package, which keeps every reader in
this plugin usable fully offline.

This module is normally used indirectly, through
:func:`~gdt.missions.burstcube.detectors.BurstCubeDetectors` and the readers
in :mod:`~gdt.missions.burstcube.cbd`, :mod:`~gdt.missions.burstcube.saa`, and
:mod:`~gdt.missions.burstcube.response`, rather than called directly.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.caldb
   :inherited-members:
