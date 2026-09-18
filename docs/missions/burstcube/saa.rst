.. _burstcube-saa:
.. |BurstCubeSAA| replace:: :class:`~gdt.missions.burstcube.saa.BurstCubeSAA`

*******************************************************************
BurstCube SAA Boundary (:mod:`gdt.missions.burstcube.saa`)
*******************************************************************

|BurstCubeSAA| reads the South Atlantic Anomaly boundary polygon from the
CALDB SAA region file, modeled on ``gdt.missions.fermi.gbm.saa``. Unlike GBM,
which has a documented history of successive SAA polygons over the mission,
BurstCube's CALDB carries a single, unversioned polygon, so there is no
time-indexed collection here.

    >>> from gdt.missions.burstcube.saa import BurstCubeSAA
    >>> saa = BurstCubeSAA()
    >>> saa.contains(-40.0, -25.0)
    True

This is a static geographic (longitude, latitude) polygon, not a time
interval -- it needs a ground position (e.g. from
:mod:`~gdt.missions.burstcube.orbit`) to evaluate against, not a time. See
the :ref:`notebooks` (notebook 3) for a worked example overlaying a real
day's ground track on the boundary.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.saa
   :inherited-members:
