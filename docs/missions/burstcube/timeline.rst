.. _burstcube-timeline:
.. |BurstCubeTimeline| replace:: :class:`~gdt.missions.burstcube.timeline.BurstCubeTimeline`

*******************************************************************
BurstCube Mission Timeline (:mod:`gdt.missions.burstcube.timeline`)
*******************************************************************

|BurstCubeTimeline| reads ``trend/*_timeline_final.csv``, a reconstructed log
of spacecraft/instrument events: a 3-column file with no header row (MET,
UTC-string, event description).

**The CSV's own UTC column is wrong by exactly +37 s**, verified against its
own MET column, evidently from treating MET as TAI seconds since the epoch
and omitting the TAI-UTC leap-second offset. This reader parses and converts
the MET column (authoritative, exposed as
:attr:`~gdt.missions.burstcube.timeline.BurstCubeTimeline.time`) and exposes
the CSV's own UTC string separately, unconverted, as
:attr:`~gdt.missions.burstcube.timeline.BurstCubeTimeline.utc_as_written` --
**do not use that column for analysis.**

    >>> from gdt.missions.burstcube.timeline import BurstCubeTimeline
    >>> timeline = BurstCubeTimeline.open('20250606_timeline_final.csv')
    >>> timeline.matching('Spacecraft Reboot').sum()
    68

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.timeline
   :inherited-members:
