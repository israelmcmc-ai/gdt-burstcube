"""BurstCube time formats, registered with :class:`astropy.time.Time`.

The BurstCube MET epoch is **2021-01-01 00:00:00 TAI**. This is a corrected
value, not the one the archive's own FITS headers state.

Every BurstCube science file carries ``MJDREFI = 59215`` and ``MJDREFF =
0.00080074074074074``, ``TIMESYS = 'TT'``. Read at face value, per the OGIP
convention that ``MJDREF`` is expressed in the scale named by ``TIMESYS``,
that resolves to 2021-01-01 00:00:00 **UTC** written in TT (``0.00080074074074074
* 86400 = 69.184`` seconds is ``32.184`` TT-TAI plus ``37`` TAI-UTC in 2021).
That is wrong by exactly one leap second's worth of confusion: the true epoch
is 2021-01-01 00:00:00 TAI, 37.000 s earlier in absolute terms.

The evidence: GRB 240629A (Fermi GBM trigger bn240629704, 2024-06-29T16:53:52.729
UTC) appears in BurstCube CBD at t0+34.14 s under the archive's stated epoch;
correcting to the TAI epoch moves that to t0-2.86 s, bringing BurstCube into
agreement with the GBM trigger within the true timing uncertainty (the
residual is itself evidence that the archive's quoted ``TIME_SYST_ERROR`` is
underestimated in at least some periods -- see the README's Caveats section).
The trend timeline CSV's own UTC column, which a naive reading of the FITS
headers makes look wrong by 37 s, is in fact correct under this epoch.

``BurstCubeSecTime`` (registered as the ``'burstcube'`` format) uses the
corrected epoch, so every :class:`~astropy.time.Time` this package returns
is 37 s earlier than converting the same MET with the archive's own stated
epoch would give. MET remains a continuous count of TT-rate seconds with no
leap-second bookkeeping -- simpler than Fermi, which counts UTC seconds
(including leap seconds) since its own epoch -- only the instant the count
starts from has changed. See :func:`check_met_epoch` for how a file's own
header is checked against this on read, and ``headers.py`` for what this
package writes.
"""
import re
import warnings

import erfa
import numpy as np
from astropy.time import Time, TimeFromEpoch, TimeUnique
from astropy.time.utils import day_frac

__all__ = ['BurstCubeSecTime', 'BurstCubeObsId', 'Time', 'check_met_epoch',
          'UnrecognizedMETEpochWarning']


class BurstCubeSecTime(TimeFromEpoch):
    """Represents the number of elapsed seconds since 2021-01-01 00:00:00
    **TAI**. This is the BurstCube Mission Elapsed Time (MET) -- using the
    corrected epoch, not the one the archive's own FITS headers state (see
    the module docstring). Unlike Fermi MET, there is no leap-second
    bookkeeping: BurstCube MET is a continuous count of seconds.
    """
    name = 'burstcube'
    """(str): Name of the mission"""

    unit = 1.0 / 86400
    """(float): unit in days"""

    epoch_val = '2021-01-01 00:00:00'
    """(str): The epoch, in Terrestrial Time"""

    epoch_val2 = None

    epoch_scale = 'tai'
    """(str): The scale of :attr:`epoch_val` -- the corrected epoch scale;
    the archive's own headers state ``'utc'``, which is the defect this
    epoch corrects (see the module docstring)."""

    epoch_format = 'iso'
    """(str): Format of :attr:`epoch_val`"""


class UnrecognizedMETEpochWarning(UserWarning):
    """Warned by :func:`check_met_epoch` when a file's own
    ``MJDREFI``/``MJDREFF``/``TIMESYS`` resolve to an instant that is
    neither this package's corrected MET epoch (2021-01-01 00:00:00 TAI)
    nor the known archive defect (2021-01-01 00:00:00 UTC). A distinct
    class from the plain ``UserWarning`` :func:`check_met_epoch` raises for
    the known defect, since this case is unexpected and worth a louder
    signal: something about the file does not match either explanation this
    package knows about.
    """


#: The corrected BurstCube MET epoch, as an absolute instant.
_CORRECTED_MET_EPOCH = Time('2021-01-01 00:00:00', scale='tai')

#: The known archive defect: MJDREFI/MJDREFF/TIMESYS taken at face value (per
#: the OGIP convention) resolve to this instant instead, 37.000 s later.
_DEFECTIVE_MET_EPOCH = Time('2021-01-01 00:00:00', scale='utc')

#: Tolerance for matching a header's resolved epoch against the two known
#: instants above. Both known cases match to double-precision floating-point
#: noise (< 1 microsecond); 1 ms leaves an enormous margin while still being
#: far tighter than the 37 s gap between them.
_EPOCH_MATCH_TOL_SEC = 1e-3


def check_met_epoch(header):
    """Check what epoch a FITS header's ``MJDREFI``/``MJDREFF``/``TIMESYS``
    actually resolve to (per the OGIP convention that ``MJDREF`` is
    expressed in the scale ``TIMESYS`` names), and warn if it is not this
    package's corrected BurstCube MET epoch.

    This is driven entirely by the header's own keyword values, not
    hardcoded to always warn: a future HEASARC revision that writes the
    corrected epoch is accepted silently. See the module docstring for the
    evidence behind the correction.

    Args:
        header: A FITS header exposing ``header['MJDREFI']``,
            ``header['MJDREFF']`` and ``header['TIMESYS']`` (a raw
            :class:`astropy.io.fits.Header` or a
            :class:`~gdt.core.headers.Header`).
    """
    try:
        mjdrefi = header['MJDREFI']
        mjdreff = header['MJDREFF']
        timesys = str(header['TIMESYS']).strip().lower()
    except KeyError:
        return

    try:
        header_epoch = Time(mjdrefi, mjdreff, format='mjd', scale=timesys)
    except Exception:
        warnings.warn(
            f"This file's MJDREFI={mjdrefi!r}, MJDREFF={mjdreff!r}, "
            f"TIMESYS={timesys!r} could not be interpreted as a time at "
            'all, so the MET epoch it states could not be checked against '
            "this package's corrected epoch "
            f'({_CORRECTED_MET_EPOCH.tai.iso} TAI).',
            UnrecognizedMETEpochWarning, stacklevel=3)
        return

    if abs((header_epoch - _CORRECTED_MET_EPOCH).sec) < _EPOCH_MATCH_TOL_SEC:
        # Already the corrected epoch -- e.g. a future HEASARC revision that
        # fixed the headers. Nothing to warn about.
        return

    if abs((header_epoch - _DEFECTIVE_MET_EPOCH).sec) < _EPOCH_MATCH_TOL_SEC:
        warnings.warn(
            "This file's MJDREFI/MJDREFF/TIMESYS state an epoch of "
            '2021-01-01 00:00:00 UTC, which is 37 s later than the true '
            'MET epoch, 2021-01-01 00:00:00 TAI (see the gdt.missions.'
            'burstcube.time module docstring and the README Caveats '
            "section). This package uses the corrected TAI epoch for every "
            "time it returns, regardless of what this header states.",
            UserWarning, stacklevel=3)
        return

    warnings.warn(
        "This file's MJDREFI/MJDREFF/TIMESYS resolve to "
        f'{header_epoch.tai.iso} TAI, which is neither this package\'s '
        f'corrected MET epoch ({_CORRECTED_MET_EPOCH.tai.iso} TAI) nor the '
        f'known archive defect ({_DEFECTIVE_MET_EPOCH.tai.iso} TAI). This '
        'package still assumes the corrected epoch regardless -- verify '
        'this file before trusting any time it returns.',
        UnrecognizedMETEpochWarning, stacklevel=3)


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
