.. _burstcube-detectors:

************************************************************************
BurstCube Detector Definitions (:mod:`gdt.missions.burstcube.detectors`)
************************************************************************

BurstCube flies four scintillator detectors, ``CS0`` through ``CS3``, in a
square pyramid, each 45 degrees off the spacecraft +z (boresight) axis and
90 degrees apart in azimuth. Their pointings are derived numerically from the
bundled CALDB alignment files rather than hardcoded, and confirm that
geometry exactly: CS0 is at (zenith, azimuth) = (45, 135) deg, CS1 at
(45, 225), CS2 at (45, 315), and CS3 at (45, 45).

    >>> from gdt.missions.burstcube.detectors import BurstCubeDetectors
    >>> BurstCubeDetectors.CS0.number
    0
    >>> BurstCubeDetectors.from_str('CS2')
    BurstCubeDetectors.CS2

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.detectors
   :inherited-members:
