.. _burstcube-headers:

**************************************************************
BurstCube FITS Headers (:mod:`gdt.missions.burstcube.headers`)
**************************************************************

FITS header definitions for every BurstCube data product, modeled on
``gdt.missions.fermi.gbm.headers`` with the keywords BurstCube adds that GBM
does not need -- ``PROCVER``, ``CALDBVER``, ``SEQPNUM``, and the
``TIMEPIXR``/``TIMEDEL`` pair (GBM's PHAII bins are always start-of-bin, so it
has no need for ``TIMEPIXR``). Every keyword list here was checked directly
against real archive files rather than the mission documentation alone, since
the two do not always agree -- see :mod:`~gdt.missions.burstcube.cbd` for the
most consequential place that mattered (the ``STDGTI`` schema).

This module is used internally by the readers in this plugin; user code
rarely needs to import it directly.

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.headers
   :inherited-members:
