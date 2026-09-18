.. _burstcube-hk:
.. |BurstCubeHK| replace:: :class:`~gdt.missions.burstcube.hk.BurstCubeHK`

*******************************************************************
BurstCube Detector Housekeeping (:mod:`gdt.missions.burstcube.hk`)
*******************************************************************

|BurstCubeHK| reads the ``DETECTOR_HK1`` (1 s cadence) and
``DETECTOR_HK2`` (60 s cadence) extensions of ``auxil/bcYYMMDDcsa.hk.gz``,
each exposed as an :class:`~astropy.timeseries.TimeSeries`. Per archive
caveat #7, the per-detector energy thresholds were raised mid-mission to
suppress a non-Poissonian low-energy noise component.
:meth:`~gdt.missions.burstcube.hk.BurstCubeHK.base_threshold` is where that
shows up: 82 -> 238 mV on ``CS0``, and similarly on the other three. The
new values were tried temporarily before being made permanent, so a single
day's file can hold both.

.. warning::
    :meth:`~gdt.missions.burstcube.hk.BurstCubeHK.peak_threshold` is *not*
    that threshold, despite the name. ``PEAK_THRES`` is a pulse-shape cut --
    how far a candidate sample must stand above its neighbours within the
    ``PHA_WIN`` window to count as a peak rather than a shoulder -- where
    ``BASE_THRES`` measures height above the running baseline. It reads
    6.0 mV for every detector in every file in the archive and never
    changes.

Neither column is described in the archive caveats document or the archive
guide. The ``DETECTOR_HK2`` ``TTYPE`` comments are the only documentation
there is:

.. code-block::

   BASE_THRES   mV   Min lvl past baseline avg x valid IDAB peak
   PEAK_THRES   mV   Min lvl past pre-post window x valid IDAB peak
   PHA_WIN           Samples pre-post potential peak for IDAB events

There is no equivalent in any Fermi GBM data product, so the GBM
convention is no guide here.

Per-detector accessors accept an int, a plain detector-name string (e.g.
``'CS0'``), or a :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`
member:

    >>> from gdt.missions.burstcube.hk import BurstCubeHK
    >>> hk = BurstCubeHK.open('bc240530csa.hk.gz')
    >>> hk.det_enabled('CS0').mean()
    np.float64(0.9966824644549763)

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.hk
   :inherited-members:
