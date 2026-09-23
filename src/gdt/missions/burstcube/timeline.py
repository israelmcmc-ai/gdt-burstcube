"""The BurstCube mission timeline: a reconstructed log of spacecraft/instrument
events from ``trend/timeline/*_timeline_final.csv``, a 3-column file with no
header row: ``MET, UTC-string, event-description``.

**The CSV's own UTC column is correct.** It agrees with the MET column
converted through :class:`~gdt.missions.burstcube.time.BurstCubeSecTime` --
this package's corrected MET epoch, 2021-01-01 00:00:00 TAI -- to better
than a millisecond on every row checked. It was the archive's own FITS
headers that were wrong: read at face value, their ``MJDREFI``/``MJDREFF``/
``TIMESYS`` state an epoch of 2021-01-01 00:00:00 UTC, 37 s later than the
truth (see the ``gdt.missions.burstcube.time`` module docstring and the
README Caveats section for the evidence, including GRB 240629A). This
reader:

* parses the MET column and converts it with :class:`~gdt.missions.burstcube.time.BurstCubeSecTime`
  (exposed as :attr:`~BurstCubeTimeline.time`);
* exposes the CSV's own UTC column separately, unconverted, as
  :attr:`~BurstCubeTimeline.utc_as_written`, for reference against the raw
  file -- kept under this name (rather than renamed now that it is known to
  agree) to avoid unnecessary churn.

Per archive caveat #6, the timeline is reconstructed from mixed packet types
at different cadences, so events can appear out of order by seconds; this
reader does not re-sort them.
"""
import csv
import re
from pathlib import Path
from typing import Union

import numpy as np

from .time import Time

__all__ = ['BurstCubeTimeline']

# "CBD Data Timestamp Drift (GPS Reboot). Unnaccounted absolute time drift: X"
# -- note the source's own misspelling "Unnaccounted", reproduced verbatim.
_DRIFT_PATTERN = re.compile(r'Unnaccounted absolute time drift:\s*([+-]?\d+(?:\.\d+)?)')


class BurstCubeTimeline:
    """The BurstCube mission timeline.

    Args:
        met (numpy.ndarray): The MET (column 0) of each event
        utc_as_written (numpy.ndarray of str): The CSV's own UTC string
            (column 1) of each event -- see the module docstring for why
            this is untrustworthy
        description (numpy.ndarray of str): The event description (column 2),
            verbatim
    """

    def __init__(self, met, utc_as_written, description):
        self._met = np.asarray(met, dtype=float)
        self._utc_as_written = np.asarray(utc_as_written, dtype=object)
        self._description = np.asarray(description, dtype=object)
        self._drift_seconds = self._parse_drift(self._description)

    @property
    def num_rows(self):
        """(int): The number of timeline events."""
        return self._met.size

    @property
    def time(self):
        """(astropy.time.Time): The authoritative event times, from the MET
        column, converted with :class:`~gdt.missions.burstcube.time.BurstCubeSecTime`.
        """
        return Time(self._met, format='burstcube')

    @property
    def utc_as_written(self):
        """(numpy.ndarray of str): The CSV's own UTC string for each event.
        Agrees with :attr:`time` (converted to UTC) to better than a
        millisecond for every row this was verified against -- the name
        predates that finding (it was originally believed to disagree by
        37 s, which turned out to be a defect in the archive's FITS headers,
        not in this CSV) and is kept as-is to avoid churn. Kept as a
        separate, unconverted property for reference/debugging against the
        raw file; prefer :attr:`time` for analysis, since it carries full
        :class:`~astropy.time.Time` precision rather than millisecond text.
        """
        return self._utc_as_written

    @property
    def description(self):
        """(numpy.ndarray of str): The raw event description for each event,
        exactly as written in the CSV -- including the source's own
        misspellings ("Continous", "Unnaccounted"). Not normalized or
        corrected; use :meth:`matching` to search it.
        """
        return self._description

    @property
    def drift_seconds(self):
        """(numpy.ndarray of float): The drift magnitude, in seconds, parsed
        out of the "CBD Data Timestamp Drift (GPS Reboot)" rows' own
        description text. ``NaN`` for every other row.
        """
        return self._drift_seconds

    def matching(self, substring):
        """Select the events whose description contains a substring.

        Note:
            Matching is a plain, case-sensitive substring search against the
            raw description text. The source's own two misspellings are not
            normalized, so match "Continous" and "Unnaccounted" (not
            "Continuous"/"Unaccounted") if searching for those events.

        Args:
            substring (str): The substring to match, e.g.
                ``'Spacecraft Reboot'`` or ``'Unnaccounted'``.

        Returns:
            (numpy.ndarray of bool): A mask into this timeline's rows.
        """
        return np.array([substring in d for d in self._description])

    @classmethod
    def open(cls, file_path: Union[str, Path]):
        """Read a BurstCube timeline CSV file.

        Args:
            file_path (str): The path to the ``*_timeline_final.csv`` file

        Returns:
            (:class:`BurstCubeTimeline`)
        """
        met, utc, description = [], [], []
        with open(file_path, newline='') as fh:
            for row in csv.reader(fh):
                if not row:
                    continue
                met.append(float(row[0]))
                utc.append(row[1])
                description.append(row[2])
        return cls(met, utc, description)

    @staticmethod
    def _parse_drift(description):
        drift = np.full(len(description), np.nan)
        for i, text in enumerate(description):
            match = _DRIFT_PATTERN.search(text)
            if match:
                drift[i] = float(match.group(1))
        return drift
