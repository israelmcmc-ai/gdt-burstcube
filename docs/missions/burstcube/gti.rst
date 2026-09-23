.. _burstcube-gti:
.. |BurstCubeGTI| replace:: :class:`~gdt.missions.burstcube.gti.BurstCubeGTI`

**************************************************
BurstCube GTI (:mod:`gdt.missions.burstcube.gti`)
**************************************************

|BurstCubeGTI| reads the standalone BurstCube trend GTI files
(``trend/gti_binning``, ``trend/gti_poscnt``, ``trend/gti_saa``), each a
single ``STDGTI`` extension of ``START``/``STOP`` pairs, and carries the set
operations gdt-core does not provide as class methods.

:class:`~gdt.core.data_primitives.Gti` already has ``intersection`` and a
gap-merging ``merge``, which |BurstCubeGTI| exposes as
:meth:`~gdt.missions.burstcube.gti.BurstCubeGTI.intersect` and
:meth:`~gdt.missions.burstcube.gti.BurstCubeGTI.union` for symmetry. The two
that are new are
:meth:`~gdt.missions.burstcube.gti.BurstCubeGTI.complement` (the gaps in a
GTI over a bounding range) and
:meth:`~gdt.missions.burstcube.gti.BurstCubeGTI.apply_to` (cutting a
:class:`~gdt.core.phaii.Phaii` or :class:`~gdt.core.tte.PhotonList` down to a
GTI). All four take and return plain
:class:`~gdt.core.data_primitives.Gti` objects, so they apply to any GTI, not
only to one read from a trend file.

    >>> from gdt.missions.burstcube.gti import BurstCubeGTI
    >>> saa = BurstCubeGTI.open('bc_csa_saa_in_cl.gti')
    >>> saa.gti.num_intervals
    586
    >>> clean = BurstCubeGTI.apply_to(saa.gti, cbd_uf)

Why ``apply_to`` exists
=======================

It is not a convenience wrapper: gdt-core's ``slice_time`` raises on a real
GTI, for two reasons that both come from the shapes these files have rather
than from anything BurstCube-specific about slicing. Most intervals of a
mission-long trend GTI select no data, and an empty segment's ``time_range``
is ``None``, which raises inside the primitive; and disjoint intervals empty
the cumulative GTI accumulator ``Phaii.slice_time`` builds. Both failures,
and what ``apply_to`` does instead, are commented in its body.

See the :ref:`notebooks`, notebook 3, for it in use.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.gti
   :inherited-members:
