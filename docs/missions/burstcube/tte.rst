.. _burstcube-tte:
.. |BurstCubeTTE| replace:: :class:`~gdt.missions.burstcube.tte.BurstCubeTTE`

****************************************************************
BurstCube Time-Tagged Events (:mod:`gdt.missions.burstcube.tte`)
****************************************************************

|BurstCubeTTE| reads a BurstCube TTE file (``events/*_tte_uf.evt.gz``):
unbinned event data for one detector, subclassing
:class:`~gdt.core.tte.PhotonList`. TTE exists for only 7 of the roughly 89
observation days in the archive.

The ``EVENTS`` extension's own ``TSTART``/``TSTOP`` header keywords are
unusable in **every** archive TTE file -- all 28 store them as strings, and
16 have ``TSTOP < TSTART`` and so a negative ``EXPOSURE``. This reader
ignores them and derives the time range from the ``STDGTI`` extension and
the event times themselves, warning when it does so.

    >>> import warnings
    >>> from gdt.missions.burstcube.tte import BurstCubeTTE
    >>> with warnings.catch_warnings(record=True) as caught:
    ...     warnings.simplefilter('always')
    ...     tte = BurstCubeTTE.open('bc240530cs0_tte_uf.evt.gz')
    >>> tte.data.size
    6364

TTE also does not cover its span continuously: events arrive in short
recording blocks separated by gaps of comparable length, and across a gap
the event list is empty while the detector keeps counting at its normal
rate. Any rate taken over an interval that spans a gap is therefore diluted.
:meth:`~gdt.missions.burstcube.tte.BurstCubeTTE.recording_blocks` returns the
blocks as a :class:`~gdt.core.data_primitives.Gti`, and
:meth:`~gdt.missions.burstcube.gti.BurstCubeGTI.complement` turns them into
the gaps.

    >>> blocks = tte.recording_blocks()
    >>> blocks.num_intervals
    46
    >>> sum(stop - start for start, stop in blocks.as_list())  # live time
    9.051...
    >>> tte.time_range[1] - tte.time_range[0]                  # elapsed span
    295.07...

See the README caveat *TTE gaps*, and ``examples/tte_gaps_vs_cbd.py``, for
the measurement behind this.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.tte
   :inherited-members:
