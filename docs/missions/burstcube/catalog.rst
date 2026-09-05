.. _burstcube-catalog:
.. |BurstCubeObsCatalog| replace:: :class:`~gdt.missions.burstcube.catalog.BurstCubeObsCatalog`

*********************************************************************
BurstCube Observation Catalog (:mod:`gdt.missions.burstcube.catalog`)
*********************************************************************

|BurstCubeObsCatalog| interfaces with the BurstCube observation catalog
(``burcbmastr``) via HEASARC Browse, subclassing
:class:`~gdt.core.heasarc.BrowseCatalog`.

    >>> from gdt.missions.burstcube.catalog import BurstCubeObsCatalog
    >>> catalog = BurstCubeObsCatalog(cached=False, verbose=False)
    >>> catalog.tte_days().num_rows
    7

Per the catalog schema, ``num_events`` is 4 on exactly the days with TTE data
and 0 elsewhere, so ``num_events > 0`` (via
:meth:`~gdt.missions.burstcube.catalog.BurstCubeObsCatalog.tte_days`) is the
supported way to find TTE days. Many rows have empty exposure/rate/
integration-time fields for detectors with no CBD that day; the per-detector
accessors here return those as NaN, never coerced to 0.0.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.catalog
   :inherited-members:
