.. _burstcube-time:

*******************************************************
BurstCube Time (:mod:`gdt.missions.burstcube.time`)
*******************************************************

The BurstCube Mission Elapsed Time (MET) epoch is **2021-01-01 00:00:00
TAI**. This is a corrected value: every BurstCube science file carries
``MJDREFI = 59215``, ``MJDREFF = 0.00080074074074074``, and ``TIMESYS =
'TT'``, which read at face value (per the OGIP convention that ``MJDREF`` is
expressed in the scale ``TIMESYS`` names) resolve to 2021-01-01 00:00:00
**UTC** instead -- 37 s later, and a known defect in the archive's own
headers, not in this package. See the README's Caveats section for the
evidence, including GRB 240629A. MET is a continuous count of seconds from
the corrected epoch, with no leap-second bookkeeping, unlike Fermi MET.
:class:`~gdt.missions.burstcube.time.BurstCubeSecTime` registers this as an
:class:`~astropy.time.Time` format named ``'burstcube'``:

    >>> from gdt.missions.burstcube.time import Time
    >>> t = Time(107629263.202, format='burstcube')
    >>> t.iso
    '2024-05-30 17:01:03.202'

Opening a file with this package checks its own ``MJDREFI``/``MJDREFF``/
``TIMESYS`` against the corrected epoch, via
:func:`~gdt.missions.burstcube.time.check_met_epoch`, and warns if the file
still states the known-defective value (silently, if a future HEASARC
revision fixes it).

This module also registers the ``'burstcube_obsid'`` format, for the
``YYMMDD`` observation-day strings used throughout the archive and the
finders:

    >>> from gdt.missions.burstcube.time import Time
    >>> t = Time('240530', format='burstcube_obsid')
    >>> t.burstcube
    107568037.0

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.time
   :inherited-members:
