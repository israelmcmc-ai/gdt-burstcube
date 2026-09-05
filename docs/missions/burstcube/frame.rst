.. _burstcube-frame:
.. |BurstCubeFrame| replace:: :class:`~gdt.missions.burstcube.frame.BurstCubeFrame`
.. |Quaternion| replace:: :class:`~gdt.core.coords.Quaternion`

****************************************************************
BurstCube Spacecraft Frame (:mod:`gdt.missions.burstcube.frame`)
****************************************************************

|BurstCubeFrame| is the frame aligned with BurstCube spacecraft coordinates
(azimuth and elevation, or equivalently zenith), defined by a quaternion that
rotates from spacecraft coordinates to ICRS. Real attitude is reconstructed
for only 3 brief epochs across the entire mission (see
:mod:`~gdt.missions.burstcube.attitude`); for every other time, build the
frame directly from a quaternion of your own:

    >>> from astropy.time import Time
    >>> from gdt.missions.burstcube.frame import BurstCubeFrame
    >>> # QPARAM is stored (and expected here) scalar-last: (x, y, z, w)
    >>> quaternion = [0.0844, 0.6314, 0.6416, 0.4273]
    >>> frame = BurstCubeFrame.from_quaternion(quaternion, Time('2024-06-29T16:53:28'))

Now a sky position can be rotated into the spacecraft frame:

    >>> from astropy.coordinates import SkyCoord
    >>> coord = SkyCoord(ra=48.723068, dec=10.861862, unit='deg')
    >>> sc = coord.transform_to(frame)
    >>> (sc.az.deg, sc.el.deg)
    (array([283.63702267]), array([90.]))

This example is the file's own row 1 boresight, so ``el=90`` (directly along
the boresight) is expected; note that azimuth is not physically meaningful
exactly at that pole, the same way RA is undefined exactly at a celestial
pole.

The quaternion component order (scalar-last, matching the
:class:`~gdt.core.coords.Quaternion` default) was established empirically
against the real attitude file's own ``POINTING`` column -- see the module's
own docstring for the full derivation.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.frame
   :inherited-members:
