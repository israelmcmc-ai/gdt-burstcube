.. _burstcube-cbd:
.. |BurstCubeCbd| replace:: :class:`~gdt.missions.burstcube.cbd.BurstCubeCbd`

********************************************************************
BurstCube Continuous Binned Data (:mod:`gdt.missions.burstcube.cbd`)
********************************************************************

|BurstCubeCbd| reads a BurstCube CBD file (``monitor/*_3cbd_{uf,cl}.fits.gz``):
a 16-channel time history of counts for one detector, subclassing
:class:`~gdt.core.phaii.Phaii`.

Two things this reader gets right that a naive port from GBM would not:

* ``TIMEPIXR=1``: the file's single ``TIME`` column is the **end** of each
  bin, not the start. Getting this backwards silently shifts every light
  curve by one bin width (0.256 s).
* The time grid is not uniform -- the instrument was frequently off or
  unstable, leaving gaps from seconds to tens of minutes. Bin edges come from
  the real ``TIME`` values, and gdt-core's contiguous-segment detection
  ensures a gap is never bridged.

    >>> from gdt.missions.burstcube.cbd import BurstCubeCbd
    >>> cbd = BurstCubeCbd.open('bc240530cs0_3cbd_cl.fits.gz')
    >>> cbd.detector
    'CS0'
    >>> cbd.data.num_chans
    16

See the :ref:`notebooks` (notebook 1) for a worked example including a real
rebinning failure this file's own genuine 1 ms-short bins can trigger, and
the honest workaround.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.cbd
   :inherited-members:
