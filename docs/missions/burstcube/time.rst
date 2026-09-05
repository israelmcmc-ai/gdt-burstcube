.. _burstcube-time:

*******************************************************
BurstCube Time (:mod:`gdt.missions.burstcube.time`)
*******************************************************

BurstCube science files carry ``MJDREFI = 59215``, ``MJDREFF =
0.00080074074074074``, and ``TIMESYS = 'TT'``. That epoch is
2021-01-01 00:00:00 UTC expressed in TT, and BurstCube Mission Elapsed Time
(MET) is a continuous count of TT seconds from it -- with no leap-second
bookkeeping, unlike Fermi MET. :class:`~gdt.missions.burstcube.time.BurstCubeSecTime`
registers this as an :class:`~astropy.time.Time` format named ``'burstcube'``:

    >>> from gdt.missions.burstcube.time import Time
    >>> t = Time(107629263.202, format='burstcube')
    >>> t.iso
    '2024-05-30 17:02:12.386'

This module also registers the ``'burstcube_obsid'`` format, for the
``YYMMDD`` observation-day strings used throughout the archive and the
finders:

    >>> from gdt.missions.burstcube.time import Time
    >>> t = Time('240530', format='burstcube_obsid')
    >>> t.burstcube
    107568000.0

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.time
   :inherited-members:
