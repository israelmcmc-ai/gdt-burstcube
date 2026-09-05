.. _burstcube-gti:
.. |BurstCubeGti| replace:: :class:`~gdt.missions.burstcube.gti.BurstCubeGti`

**************************************************
BurstCube GTI (:mod:`gdt.missions.burstcube.gti`)
**************************************************

Readers and set operations for the standalone BurstCube trend GTI files
(``trend/gti_binning``, ``trend/gti_poscnt``, ``trend/gti_saa``).
:class:`~gdt.core.data_primitives.Gti` already provides ``intersection`` and
a gap-merging ``merge``; this module adds :func:`~gdt.missions.burstcube.gti.complement`
(the gaps in a GTI over a bounding range) and
:func:`~gdt.missions.burstcube.gti.apply_to` (cutting a
:class:`~gdt.core.phaii.Phaii` or :class:`~gdt.core.tte.PhotonList` down to a
GTI).

Note:
    ``apply_to`` is a thin wrapper over gdt-core's own ``Phaii.slice_time``,
    which intersects multiple requested time ranges' own GTIs cumulatively
    rather than as a union. Passing a GTI with more than one disjoint segment
    in a single call can therefore raise ``TypeError`` from gdt-core itself;
    apply each segment separately and merge the results instead (see the
    :ref:`notebooks`, notebook 3, for a worked example of both the failure
    and the workaround).

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.gti
   :inherited-members:
