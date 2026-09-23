.. _burstcube-attitude:
.. |BurstCubeAttitude| replace:: :class:`~gdt.missions.burstcube.attitude.BurstCubeAttitude`

***********************************************************
BurstCube Attitude (:mod:`gdt.missions.burstcube.attitude`)
***********************************************************

|BurstCubeAttitude| reads ``trend/attitude/bc_csa_att.fits``, which holds
exactly **3** manually reconstructed attitude epochs across the entire
archive, each with an approximate 10-12 degree 1-sigma pointing error. This
reader deliberately does not interpolate or extrapolate between them -- with
only 3 widely-spaced epochs and no information about the attitude in
between, any interpolation would be fabricating data.

    >>> from gdt.missions.burstcube.attitude import BurstCubeAttitude
    >>> att = BurstCubeAttitude.open('bc_csa_att.fits')
    >>> att.num_rows
    3
    >>> att.pointing[0]  # (RA, Dec, roll) of the boresight, in degrees
    array([ 48.723068,  10.861862, 333.95041 ], dtype='>f8')

For every other time -- essentially all of the mission -- there is no
attitude file to open. The documented path is a user-supplied quaternion via
:meth:`~gdt.missions.burstcube.frame.BurstCubeFrame.from_quaternion`; see
:mod:`~gdt.missions.burstcube.frame` and the :ref:`notebooks` (notebook 3).

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.attitude
   :inherited-members:
