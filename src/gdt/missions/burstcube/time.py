"""BurstCube time formats, registered with :class:`astropy.time.Time`.

BurstCube science files carry ``MJDREFI = 59215`` and ``MJDREFF =
0.00080074074074074``, ``TIMESYS = 'TT'``.  MJD 59215 is 2021-01-01, and
``0.00080074074074074 * 86400 = 69.184`` seconds, which is exactly
``32.184`` (TT-TAI) plus ``37`` (TAI-UTC in 2021).  In other words, the MET
epoch is 2021-01-01 00:00:00 UTC expressed in TT, and MET is a continuous
count of TT seconds with no leap-second bookkeeping -- simpler than Fermi,
which counts UTC seconds (including leap seconds) since its own epoch.
"""
import re

import erfa
import numpy as np
from astropy.time import Time, TimeFromEpoch, TimeUnique
from astropy.time.utils import day_frac

__all__ = ['BurstCubeSecTime', 'BurstCubeObsId', 'Time']


class BurstCubeSecTime(TimeFromEpoch):
    """Represents the number of elapsed TT seconds since 2021-01-01 00:00:00
    UTC (expressed as ``2021-01-01 00:01:09.184`` TT). This is the BurstCube
    Mission Elapsed Time (MET). Unlike Fermi MET, there is no leap-second
    bookkeeping: BurstCube MET is a continuous count of TT seconds.
    """
    name = 'burstcube'
    """(str): Name of the mission"""

    unit = 1.0 / 86400
    """(float): unit in days"""

    epoch_val = '2021-01-01 00:01:09.184'
    """(str): The epoch in Terrestrial Time"""

    epoch_val2 = None

    epoch_scale = 'tt'
    """(str): The scale of :attr:`epoch_val`"""

    epoch_format = 'iso'
    """(str): Format of :attr:`epoch_val`"""


class BurstCubeObsId(TimeUnique):
    """Represents the BurstCube observation-day directory name, ``YYMMDD``,
    used throughout the archive (e.g. ``obs/2024_05/240530/``). This allows
    a directory name to be used directly wherever a :class:`~astropy.time.Time`
    is expected, e.g. ``Time('240530', format='burstcube_obsid')``.

    The represented instant is midnight UTC at the start of that day.
    """
    name = 'burstcube_obsid'
    """(str): The name of the time format"""

    _obsid_pattern = re.compile(r'^(?P<year>\d\d)(?P<month>\d\d)(?P<day>\d\d)$')

    def _check_val_type(self, val1, val2):
        if not all(self._obsid_pattern.match(val) is not None for val in val1.flat):
            raise TypeError('Input values for {} class must be BurstCube '
                            'observation IDs (YYMMDD)'.format(self.name))
        if val2 is not None:
            raise ValueError(
                f'{self.name} objects do not accept a val2 but you provided {val2}')
        return val1, None

    def set_jds(self, val1, val2):
        """Convert the observation ID contained in val1 to jd1, jd2"""
        iterator = np.nditer([val1, None, None, None],
                             flags=['refs_ok', 'zerosize_ok'],
                             op_dtypes=[None] + 3 * [np.intc])

        for val, iy, im, iday in iterator:
            m = self._obsid_pattern.match(str(val))
            if not m:
                raise ValueError('observation ID was not the correct format.')

            obsid = m.groupdict()
            iy[...] = int(obsid['year']) + 2000
            im[...] = int(obsid['month'])
            iday[...] = int(obsid['day'])

        zeros = np.zeros(iterator.operands[1].shape, dtype=np.intc)
        zero_secs = np.zeros(iterator.operands[1].shape, dtype=np.double)
        jd1, jd2 = erfa.dtf2d(self.scale.upper().encode('ascii'),
                              iterator.operands[1], iterator.operands[2],
                              iterator.operands[3], zeros, zeros, zero_secs)
        self.jd1, self.jd2 = day_frac(jd1, jd2)

    def to_value(self, parent=None, out_subfmt=None):
        """Convert to string representation of the observation ID"""
        if out_subfmt is not None:
            self._select_subfmts(out_subfmt)

        scale = self.scale.upper().encode('ascii')
        iys, ims, ids, ihmsfs = erfa.d2dtf(scale, 0, self.jd1, self.jd2)
        iterator = np.nditer([iys, ims, ids, None],
                             flags=['refs_ok', 'zerosize_ok'],
                             op_dtypes=3 * [None] + [object])

        for iy, im, iday, out in iterator:
            out[...] = f'{iy - 2000:02d}{im:02d}{iday:02d}'

        return iterator.operands[-1]

    value = property(to_value)
