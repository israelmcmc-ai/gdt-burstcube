.. _burstcube-tte:
.. |BurstCubeTte| replace:: :class:`~gdt.missions.burstcube.tte.BurstCubeTte`

****************************************************************
BurstCube Time-Tagged Events (:mod:`gdt.missions.burstcube.tte`)
****************************************************************

|BurstCubeTte| reads a BurstCube TTE file (``events/*_tte_uf.evt.gz``):
unbinned event data for one detector, subclassing
:class:`~gdt.core.tte.PhotonList`. TTE exists for only 7 of the roughly 89
observation days in the archive.

Per archive caveat #4, the ``EVENTS`` extension's own ``TSTART``/``TSTOP``
header keywords are unreliable in at least one real file (stored as strings,
and with ``TSTOP < TSTART`` and a negative ``EXPOSURE``); this reader ignores
them and derives the time range from the ``STDGTI`` extension and the event
times themselves, warning when it does so.

    >>> import warnings
    >>> from gdt.missions.burstcube.tte import BurstCubeTte
    >>> with warnings.catch_warnings(record=True) as caught:
    ...     warnings.simplefilter('always')
    ...     tte = BurstCubeTte.open('bc240530cs0_tte_uf.evt.gz')
    >>> tte.data.size
    6364

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.tte
   :inherited-members:
