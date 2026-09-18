.. _burstcube-finders:
.. |BurstCubeObsFinder| replace:: :class:`~gdt.missions.burstcube.finders.BurstCubeObsFinder`
.. |BurstCubeTrendFinder| replace:: :class:`~gdt.missions.burstcube.finders.BurstCubeTrendFinder`

*******************************************************************
BurstCube Data Finders (:mod:`gdt.missions.burstcube.finders`)
*******************************************************************

|BurstCubeObsFinder| finds and downloads products for a single observation
day (``obs/YYYY_MM/YYMMDD/``); |BurstCubeTrendFinder| finds the mission-long
trend products (attitude, GTI, timeline), which are not organized by day.

Unlike GBM, every BurstCube archive filename is fully determined by the
observation day, detector, and product type, so neither finder needs to list
a remote directory to know what to download -- each builds the exact
relative path and downloads it directly. A missing file (including an entire
missing observation day; day directories are **not** contiguous) is a
per-file 404, handled by skipping with a warning rather than raising.

    >>> from gdt.missions.burstcube.finders import BurstCubeObsFinder
    >>> finder = BurstCubeObsFinder('240530')
    >>> paths = finder.get_cbd('./data', detectors='CS0', variant='cl')

|BurstCubeTrendFinder| has no per-instance navigation argument (the trend
directory is a single fixed location), so it is ready to download from
immediately on construction -- no separate ``cd()`` call is needed:

    >>> from gdt.missions.burstcube.finders import BurstCubeTrendFinder
    >>> trend_finder = BurstCubeTrendFinder()
    >>> att_path = trend_finder.get_attitude('./data')

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.finders
   :inherited-members:
