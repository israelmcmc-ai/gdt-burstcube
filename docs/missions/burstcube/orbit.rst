.. _burstcube-orbit:
.. |BurstCubeOrbit| replace:: :class:`~gdt.missions.burstcube.orbit.BurstCubeOrbit`

**********************************************************
BurstCube Orbit Data (:mod:`gdt.missions.burstcube.orbit`)
**********************************************************

|BurstCubeOrbit| reads the ``ORBIT`` extension of ``auxil/bcYYMMDD.hk.gz``:
the spacecraft position and velocity in J2000 ECI coordinates, reconstructed
from TLEs. This file spans the full on-orbit lifetime, including days with no
science data at all, and its time is independent of -- and more accurate
than -- the instrument clock.

There is no attitude information in this file. Calling
:meth:`~gdt.missions.burstcube.orbit.BurstCubeOrbit.get_spacecraft_frame`
gives a :class:`~gdt.missions.burstcube.frame.BurstCubeFrame` with
``obsgeoloc``/``obsgeovel`` set and no quaternion -- usable for
Earth-visibility and geocenter calculations, but not for converting sky
positions to/from BurstCube coordinates. Combine it with a user-supplied
quaternion (see :mod:`~gdt.missions.burstcube.frame`) for that.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.orbit
   :inherited-members:
