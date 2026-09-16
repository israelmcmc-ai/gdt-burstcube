.. _burstcube-hk:
.. |BurstCubeHK| replace:: :class:`~gdt.missions.burstcube.hk.BurstCubeHK`

*******************************************************************
BurstCube Detector Housekeeping (:mod:`gdt.missions.burstcube.hk`)
*******************************************************************

|BurstCubeHK| reads the ``DETECTOR_HK1`` (1 s cadence) and
``DETECTOR_HK2`` (60 s cadence) extensions of ``auxil/bcYYMMDDcsa.hk.gz``,
each exposed as an :class:`~astropy.timeseries.TimeSeries`. Per archive
caveat #7, the per-detector energy thresholds were raised mid-mission from
~21-30 keV to ~100 keV to suppress a non-Poissonian low-energy noise
component; :meth:`~gdt.missions.burstcube.hk.BurstCubeHK.peak_threshold`
and :meth:`~gdt.missions.burstcube.hk.BurstCubeHK.base_threshold` are
how that shows up in the housekeeping data. A single day only samples one
side of that mission-long change.

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
