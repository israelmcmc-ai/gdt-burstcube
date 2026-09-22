.. _burstcube-timeline:
.. |BurstCubeTimeline| replace:: :class:`~gdt.missions.burstcube.timeline.BurstCubeTimeline`

*******************************************************************
BurstCube Mission Timeline (:mod:`gdt.missions.burstcube.timeline`)
*******************************************************************

|BurstCubeTimeline| reads ``trend/*_timeline_final.csv``, a reconstructed log
of spacecraft/instrument events: a 3-column file with no header row (MET,
UTC-string, event description).

**The CSV's own UTC column is correct.** It agrees with the MET column
converted through :class:`~gdt.missions.burstcube.time.BurstCubeSecTime` --
this package's corrected MET epoch, 2021-01-01 00:00:00 TAI -- to better
than a millisecond. It was the archive's own FITS headers that stated the
wrong epoch (2021-01-01 00:00:00 UTC, 37 s late); see the README's Caveats
section for the evidence, including GRB 240629A. This reader parses and
converts the MET column (exposed as
:attr:`~gdt.missions.burstcube.timeline.BurstCubeTimeline.time`) and also
exposes the CSV's own UTC string separately, unconverted, as
:attr:`~gdt.missions.burstcube.timeline.BurstCubeTimeline.utc_as_written`
(the name predates this finding -- it was originally believed to disagree by
37 s -- and is kept as-is to avoid churn).

    >>> from gdt.missions.burstcube.timeline import BurstCubeTimeline
    >>> timeline = BurstCubeTimeline.open('20250606_timeline_final.csv')
    >>> timeline.matching('Spacecraft Reboot').sum()
    68

Reference/API
=============

.. automodapi:: gdt.missions.burstcube.timeline
   :inherited-members:
