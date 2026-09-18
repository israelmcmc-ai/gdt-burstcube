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

It is not a convenience wrapper. Handing a real trend GTI straight to
gdt-core's ``slice_time`` fails, for two reasons that both come from the
shapes these files have rather than from anything BurstCube-specific about
slicing:

* **Most intervals select nothing.** A trend GTI spans the whole mission --
  the SAA one has 586 intervals over five months -- while a single CBD or
  TTE file covers minutes to hours. ``slice_time`` builds
  ``Gti.from_list([segment.time_range])`` per requested range, and an empty
  segment's ``time_range`` is ``None``, which raises ``TypeError`` from
  inside the primitive. ``apply_to`` drops those intervals first. Filtering
  on the file's overall time range is not enough: BurstCube data is gappy,
  so an interval can sit inside the span and still contain no bins.

* **Disjoint intervals empty the accumulator.** ``Phaii.slice_time``
  accumulates a GTI by *intersecting* each segment's own range into a
  running total, which empties as soon as two segments are disjoint and then
  raises on the empty result -- even though that accumulated value is
  discarded (``from_data`` is handed ``self.gti`` instead). ``apply_to``
  slices one interval per call, so that loop runs exactly once and never
  empties, then recombines with ``merge``.

See the :ref:`notebooks`, notebook 3, for a worked example of both failures
and the result.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.gti
   :inherited-members:
