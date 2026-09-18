.. _burstcube-response:
.. |BurstCubeRsp| replace:: :class:`~gdt.missions.burstcube.response.BurstCubeRsp`
.. |BurstCubeResponseGrid| replace:: :class:`~gdt.missions.burstcube.response.BurstCubeResponseGrid`
.. |BurstCubeRspFinder| replace:: :class:`~gdt.missions.burstcube.response.BurstCubeRspFinder`

********************************************************************
BurstCube Detector Response (:mod:`gdt.missions.burstcube.response`)
********************************************************************

BurstCube's detector response is a grid of Detector Response Matrices (DRMs),
one per HEALPix pixel (``NSIDE=16``, RING ordering, ``NUM_PIXELS=3072``) per
detector, in **spacecraft** coordinates -- not one response computed on the
fly per source position.

* |BurstCubeRsp| (subclassing :class:`~gdt.core.response.Rsp`) reads a single
  ``.rsp`` file: 64 native channels. Its own internal ``EBOUNDS`` is an
  older, rounded, detector-independent grid; :meth:`~gdt.missions.burstcube.response.BurstCubeRsp.to_cbd`
  regroups to CBD's 16 channels by index but replaces the channel energies
  with CALDB's per-detector ``eb16`` -- the two disagree even in the top
  channel (see the top-level README's Caveats section for the exact
  numbers).
* |BurstCubeResponseGrid| downloads (and LRU-caches) only the grid pixels
  actually requested, never the whole ~1 GB grid, and interpolates a DRM for
  a given ``(az, zen)``, HEALPix pixel, or sky position (which needs a frame
  carrying attitude).
* |BurstCubeRspFinder| figures out which pixels a direction needs before
  downloading anything.

    >>> from gdt.missions.burstcube.response import BurstCubeResponseGrid
    >>> grid = BurstCubeResponseGrid(detectors='CS0')
    >>> drm = grid.get_drm('CS0', az=45.0, zen=30.0, interp=True)
    >>> drm.num_chans
    64

See the :ref:`notebooks` (notebook 2) for a full worked example: plotting the
DRM and effective area, folding a spectrum, regrouping to CBD, and comparing
all four detectors at one location.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.response
   :inherited-members:
