# Copyright 2024-2025 by the BurstCube Team.
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
"""FITS header definitions for the BurstCube data products.

These are modeled on ``gdt.missions.fermi.gbm.headers``, with the keywords
BurstCube carries that GBM does not: ``PROCVER``, ``CALDBVER``, ``SEQPNUM``,
and the ``TIMEPIXR``/``TIMEDEL`` pair (GBM does not need ``TIMEPIXR`` because
its PHAII bins are always start-of-bin). Every keyword list below was checked
directly against real archive files (``bc240530cs0_3cbd_cl.fits.gz``,
``bc240530cs0_tte_uf.evt.gz``, ``bc240530.hk.gz``, ``bc240530csa.hk.gz``,
``bc_csa_att.fits``), not just the mission documentation, since the two do not
always agree (see the CBD/TTE ``STDGTI`` note below).
"""
from gdt.core.headers import FileHeaders, Header

from .time import Time

__all__ = ['AttitudeHeaders', 'CbdHeaders', 'CbdUnfilteredHeaders', 'DetectorHkHeaders', 'GtiHeaders',
          'OrbitHeaders', 'TteHeaders']

# mission definitions
_telescope = 'BURSTCUBE'
_origin = 'GSFC'
_timesys = 'TT'
_timeunit = 's'
_timeref = 'LOCAL'
_tassign = 'SATELLITE'
_mjdrefi = 59215
_mjdreff = 0.00080074074074074

# common keyword cards
_telescope_card = ('TELESCOP', _telescope, 'Telescope mission name')
_instrument_card = ('INSTRUME', '', 'Instrument name')
_datamode_card = ('DATAMODE', '', 'Data mode')
_obs_id_card = ('OBS_ID', '', 'Observation ID')
_origin_card = ('ORIGIN', _origin, 'Origin of the FITS file')
_date_card = ('DATE', '', 'File creation date')
_date_obs_card = ('DATE-OBS', '', 'Start Date')
_date_end_card = ('DATE-END', '', 'Stop Date')
_extname_card = ('EXTNAME', '', 'Name of this binary table extension')
_mjdrefi_card = ('MJDREFI', _mjdrefi, 'MJD reference day 01 Jan 2021 00:00:00')
_mjdreff_card = ('MJDREFF', _mjdreff, 'MJD reference (fraction of day)')
_timeref_card = ('TIMEREF', _timeref, 'Reference Frame')
_tassign_card = ('TASSIGN', _tassign, 'Time assigned')
_timesys_card = ('TIMESYS', _timesys, 'Time system')
_timeunit_card = ('TIMEUNIT', _timeunit, 'Time unit for timing header keywords')
_clockapp_card = ('CLOCKAPP', False, 'If clock corrections are applied (T/F)')
_tstart_card = ('TSTART', 0.0, 'Start Time')
_tstop_card = ('TSTOP', 0.0, 'Stop Time')
_telapse_card = ('TELAPSE', 0.0, 'TSTOP-TSTART in sec')
_ontime_card = ('ONTIME', 0.0, 'Sum good time interval in sec')
_deadapp_card = ('DEADAPP', False, 'if deadtime applied')
_hduclass_card = ('HDUCLASS', 'OGIP', 'format conforms to the OGIP standard')
_procver_card = ('PROCVER', '', 'Processing version')
_caldbver_card = ('CALDBVER', '', 'CALDB version')
_seqpnum_card = ('SEQPNUM', 1, 'Number of times the dataset has been processed')

# cards for the processing trailer common to most extensions: ORIGIN,
# CREATOR, PROCVER, CALDBVER, SEQPNUM (BurstCube-specific: not in GBM)
_bc_trailer_cards = [_origin_card, Header.creator(), _procver_card,
                     _caldbver_card, _seqpnum_card]

# the common timing block on every science-data extension
_bc_time_cards = [_mjdrefi_card, _mjdreff_card, _timeref_card, _tassign_card,
                  _timesys_card, _timeunit_card]


class BurstCubeHeader(Header):
    """A :class:`~gdt.core.headers.Header` that keeps ``DATE-OBS``/``DATE-END``
    in sync with ``TSTART``/``TSTOP`` (interpreted as BurstCube MET). Per
    archive caveat #4, ``TSTART``/``TSTOP`` are sometimes written as strings
    with ``TSTOP < TSTART`` in TTE files; this sync is best-effort and quietly
    does nothing if the value cannot be interpreted as a MET, rather than
    raising, so that opening a file with malformed timing keywords does not
    crash.
    """

    def __setitem__(self, key, val):
        if not isinstance(key, tuple) and not isinstance(val, tuple):
            date_key = {'TSTART': 'DATE-OBS', 'TSTOP': 'DATE-END'}.get(key.upper())
            if date_key is not None and date_key in self:
                try:
                    self[date_key] = Time(float(val), format='burstcube').utc.isot
                except (TypeError, ValueError):
                    pass

        super().__setitem__(key, val)


class DataPrimaryHeader(BurstCubeHeader):
    """PRIMARY header for CBD and TTE data files."""
    name = 'PRIMARY'
    keywords = [_telescope_card, _instrument_card, _datamode_card,
                _obs_id_card, _date_obs_card, _date_end_card, _origin_card,
                _date_card, Header.creator()]


class AuxPrimaryHeader(BurstCubeHeader):
    """PRIMARY header for the orbit and detector housekeeping auxiliary
    files, which (unlike CBD/TTE) carry no ``DATAMODE``.
    """
    name = 'PRIMARY'
    keywords = [_telescope_card, _instrument_card, _obs_id_card,
                _date_obs_card, _date_end_card, _origin_card, _date_card,
                Header.creator()]


class AttitudePrimaryHeader(BurstCubeHeader):
    """PRIMARY header for the attitude file, which spans the whole mission
    and so carries no ``OBS_ID`` (and, in the real archive file, no ``DATE``).
    """
    name = 'PRIMARY'
    keywords = [_telescope_card, _instrument_card, _date_obs_card,
                _date_end_card, _origin_card, Header.creator()]


class CbdDataHeader(BurstCubeHeader):
    """Header for the ``CBD`` extension. ``TIMEPIXR=1`` means the ``TIME``
    column is the *end* of each 0.256 s bin, not the start; getting this
    backwards shifts every light curve by one bin width.
    """
    name = 'CBD'
    keywords = [_extname_card, _telescope_card, _instrument_card,
                _datamode_card, _obs_id_card] + _bc_time_cards + [
                ('TIMEPIXR', 1, 'Time stamps refer to the end of each bin'),
                ('TIMEDEL', 0.256, 'Integration time'), _clockapp_card,
                _date_obs_card, _date_end_card, _tstart_card, _tstop_card,
                _deadapp_card, _telapse_card, _ontime_card, _hduclass_card,
                ('HDUCLAS1', 'ARRAY', 'hduclass1'),
                ('HDUCLAS2', 'TOTAL', 'hduclass2'),
                # WCS-like descriptors for the COUNTS column's channel axis
                ('1CDLT2', 1, 'Pixel increment along CTYP axis'),
                ('1CDPX2', 1, 'reference pixel of the CTYP axis'),
                ('1CRVL2', 0, 'Value of the CTYP axis reference pixel'),
                ('1CUNI2', 'chan', 'Unit of the pixel along the CTYP axis'),
                ('1CTYP2', 'CHANNEL', 'CTYP axis name')] + _bc_trailer_cards \
               + [_date_card]


class CbdGtiHeader(BurstCubeHeader):
    """Header for the ``STDGTI`` extension of a *cleaned* (``_cl``) CBD file:
    the 2-column ``START``/``STOP`` OGIP standard form, carrying
    ``HDUCLAS2='STANDARD'``, ``HDUVERS`` and ``TIMEZERO``.

    Note:
        The two CBD variants do not share a ``STDGTI`` schema. Unfiltered
        (``_uf``) files instead carry the 6-column form used by TTE
        (``START``, ``STOP``, ``START_ORIGINAL_TIME``,
        ``START_TIME_SYST_ERROR``, ``STOP_ORIGINAL_TIME``,
        ``STOP_TIME_SYST_ERROR``) with ``HDUCLAS1='GTI'`` and none of the
        three keywords above -- see :class:`TteGtiHeader`. Because nothing in
        the files guarantees this correlation holds for every observation,
        :meth:`~gdt.missions.burstcube.cbd.BurstCubeCbd.open` selects between
        the two by inspecting the extension itself rather than by the
        ``_uf``/``_cl`` filename.
    """
    name = 'STDGTI'
    keywords = [_extname_card, _hduclass_card,
                ('HDUCLAS1', 'GTI', 'Contains good time intervals'),
                ('HDUCLAS2', 'STANDARD', 'Contains standard good time intervals'),
                ('HDUVERS', '1.0.0', 'Version of GTI header'),
                ('TIMEZERO', 0.0, 'Zero-point offset for TIME column'),
                _tstart_card, _tstop_card, _telescope_card, _instrument_card,
                _datamode_card, _obs_id_card] + _bc_time_cards + [
                _clockapp_card, _date_obs_card, _date_end_card] \
               + _bc_trailer_cards + [_date_card]


class EventsHeader(BurstCubeHeader):
    """Header for the ``EVENTS`` extension of a TTE file.

    Per archive caveat #4, ``TSTART``/``TSTOP`` are sometimes stored as
    strings with ``TSTOP < TSTART``, giving negative ``TELAPSE``/``EXPOSURE``
    while ``ONTIME`` stays positive; readers must fall back to the ``STDGTI``
    extension and the event times themselves rather than trusting these two
    keywords. The trigger keywords (``TRIG_ID`` etc.) are always zero/empty:
    BurstCube never had an onboard astrophysical trigger fire.
    """
    name = 'EVENTS'
    keywords = [_extname_card, _telescope_card, _instrument_card,
                _datamode_card, _obs_id_card,
                ('OBJECT', 'UNKNOWN', 'Object name'),
                ('RA_OBJ', -2000.0, 'R.A. Object (placeholder: -2000)'),
                ('DEC_OBJ', -2000.0, 'Dec Object (placeholder: -2000)'),
                ('POSFLAG', 0, 'Position FLAG 1=good, 0=best effort, -1=false'),
                ('EQUINOX', 2000.0, 'Equinox'),
                ('RADECSYS', 'FK5', 'Coordinates System')] + _bc_time_cards + [
                ('TIMEDEL', 1e-5, 'Integration time'), _clockapp_card,
                _date_obs_card, _date_end_card, _tstart_card, _tstop_card,
                _telapse_card, _ontime_card,
                ('EXPOSURE', 0.0, 'Exposure'), _deadapp_card,
                ('EVTTYPE', 'T3E', 'Event Type of T3E or RTTE'),
                ('TRIG_ID', '', 'trigger id'),
                ('TRIGGER', 0, 'Trigger number'),
                ('TRIGTIME', 0.0, 'Trigger Time sec since 01 Jan 2021 00:00:00'),
                ('TRIGUTC', 0.0, 'TIME of Trigger time in UTC'),
                ('TRIGDET', '', 'Trigger Detector 0=no 1=yes'),
                ('TRIGSCAL', 0, '[s] Trigger timescale'),
                ('TRIGSIGN', 0.0, 'Trigger significance'),
                ('TRIGEMIN', 0.0, '[keV] Trigger energy min'),
                ('TRIGEMAX', 0.0, '[keV] Trigger energy max'),
                _hduclass_card,
                ('HDUCLAS1', 'EVENT', 'hduclass1'),
                ('HDUCLAS2', 'ALL', 'hduclass2')] + _bc_trailer_cards \
               + [_date_card]


class TteGtiHeader(BurstCubeHeader):
    """Header for the ``STDGTI`` extension of a TTE file. Unlike
    :class:`CbdGtiHeader`, its columns are ``START``, ``STOP``,
    ``START_ORIGINAL_TIME``, ``START_TIME_SYST_ERROR``,
    ``STOP_ORIGINAL_TIME``, ``STOP_TIME_SYST_ERROR``, and its header carries
    only ``HDUCLAS1`` (no ``HDUCLAS2``, ``HDUVERS``, or ``TIMEZERO``).
    """
    name = 'STDGTI'
    keywords = [_extname_card, _telescope_card, _instrument_card,
                _datamode_card, _obs_id_card] + _bc_time_cards + [
                _clockapp_card, _date_obs_card, _date_end_card, _tstart_card,
                _tstop_card, _hduclass_card,
                ('HDUCLAS1', 'GTI', 'hduclass1')] + _bc_trailer_cards \
               + [_date_card]


class OrbitDataHeader(BurstCubeHeader):
    """Header for the ``ORBIT`` extension. ``INSTRUME`` is ``'CS'`` here
    (spacecraft-level), not a detector name or ``'CSA'``.
    """
    name = 'ORBIT'
    keywords = [_extname_card, _telescope_card, _instrument_card,
                _obs_id_card] + _bc_time_cards + [
                _clockapp_card, _date_obs_card, _date_end_card, _tstart_card,
                _tstop_card, _hduclass_card,
                ('HDUCLAS1', 'TEMPORALDATA', 'hduclass1'),
                ('HDUCLAS2', 'EPHEM', 'hduclass2')] + _bc_trailer_cards


class AttitudeDataHeader(BurstCubeHeader):
    """Header for the ``ATTITUDE`` extension. Exactly 3 rows exist in the
    entire archive; see :mod:`gdt.missions.burstcube.frame` for how the
    ``QPARAM`` quaternion component order was established against this
    extension's ``POINTING`` column.
    """
    name = 'ATTITUDE'
    keywords = [_extname_card, _telescope_card, _instrument_card] \
               + _bc_time_cards + [
                _clockapp_card, _date_obs_card, _date_end_card, _tstart_card,
                _tstop_card, _hduclass_card,
                ('HDUCLAS1', 'TEMPORALDATA', 'hduclass1'),
                ('HDUCLAS2', 'ATTITUDE', 'hduclas2'), _origin_card,
                Header.creator(), _procver_card, _caldbver_card,
                _seqpnum_card]


class DetectorHkHeader(BurstCubeHeader):
    """Shared header shape for the ``DETECTOR_HK1`` and ``DETECTOR_HK2``
    extensions, which differ only in ``EXTNAME`` and ``TIMEDEL`` (1 s for
    HK1, 60 s for HK2). ``INSTRUME`` is always ``'CSA'``.
    """
    name = 'DETECTOR_HK'
    keywords = [_extname_card, _telescope_card, _instrument_card,
                _obs_id_card] + _bc_time_cards + [
                ('TIMEDEL', 1.0, 'Integration time'), _clockapp_card,
                _date_obs_card, _date_end_card, _tstart_card, _tstop_card,
                _hduclass_card,
                ('HDUCLAS1', 'TEMPORALDATA', 'hduclass1'),
                ('HDUCLAS2', 'HKP', 'hduclass2')] + _bc_trailer_cards \
               + [_date_card]


class DetectorHk1Header(DetectorHkHeader):
    """Header for the ``DETECTOR_HK1`` extension (1 s cadence)."""
    name = 'DETECTOR_HK1'
    keywords = DetectorHkHeader.keywords.copy()


class DetectorHk2Header(DetectorHkHeader):
    """Header for the ``DETECTOR_HK2`` extension (60 s cadence)."""
    name = 'DETECTOR_HK2'
    keywords = [kw if kw[0] != 'TIMEDEL' else ('TIMEDEL', 60.0, 'Integration time')
               for kw in DetectorHkHeader.keywords]


class GtiTrendPrimaryHeader(BurstCubeHeader):
    """PRIMARY header for a standalone trend GTI file (``trend/gti_*/*.gti``).

    Note:
        Unlike the other BurstCube products, these files are produced by a
        chain of generic HEASoft tools (``fcalc``/``fmerge``/``maketime``),
        and different trend GTI files carry noticeably different amounts of
        carried-over PRIMARY content: e.g. ``trend/gti_binning/*.gti``
        PRIMARY headers are close copies of a CBD file's header (complete
        with ``TIMEPIXR``/``TSTART``/etc. and a long ``HISTORY`` trail),
        while ``trend/gti_saa/*.gti`` PRIMARY headers carry almost nothing.
        Only the keywords common to both are kept here; the ``STDGTI``
        extension is consistent across variants and is fully captured by
        :class:`GtiTrendDataHeader`.
    """
    name = 'PRIMARY'
    keywords = [_telescope_card, _instrument_card, _origin_card,
                Header.creator()]


class GtiTrendDataHeader(BurstCubeHeader):
    """Header for the ``STDGTI`` extension of a standalone trend GTI file
    (``trend/gti_binning``, ``trend/gti_poscnt``, ``trend/gti_saa``).

    Note:
        ``HDUVERS`` and ``TIMEZERO`` appear on some trend GTI variants (e.g.
        ``trend/gti_saa``, matching :class:`CbdGtiHeader`) but not others
        (e.g. ``trend/gti_binning``), so, unlike :class:`CbdGtiHeader`, they
        are not included here.
    """
    name = 'STDGTI'
    keywords = [_extname_card, _hduclass_card,
                ('HDUCLAS1', 'GTI', 'Contains good time intervals'),
                ('HDUCLAS2', 'STANDARD', 'Contains standard good time intervals'),
                _tstart_card, _tstop_card, _telescope_card, _instrument_card,
                _datamode_card] + _bc_time_cards + [
                _clockapp_card, _date_obs_card, _date_end_card] \
               + _bc_trailer_cards + [_date_card]


#-------------------------------------

class CbdHeaders(FileHeaders):
    """FITS headers for a CBD file whose ``STDGTI`` uses the 2-column
    standard schema (in practice, the cleaned ``_cl`` files)."""
    _header_templates = [DataPrimaryHeader(), CbdDataHeader(), CbdGtiHeader()]


class CbdUnfilteredHeaders(FileHeaders):
    """FITS headers for a CBD file whose ``STDGTI`` uses the 6-column schema
    shared with TTE (in practice, the unfiltered ``_uf`` files).

    Identical to :class:`CbdHeaders` apart from the GTI extension. See the
    note on :class:`CbdGtiHeader`.
    """
    _header_templates = [DataPrimaryHeader(), CbdDataHeader(), TteGtiHeader()]


class TteHeaders(FileHeaders):
    """FITS headers for TTE (time-tagged event) files."""
    _header_templates = [DataPrimaryHeader(), EventsHeader(), TteGtiHeader()]


class OrbitHeaders(FileHeaders):
    """FITS headers for the orbit/ephemeris file (``auxil/bcYYMMDD.hk.gz``)."""
    _header_templates = [AuxPrimaryHeader(), OrbitDataHeader()]


class AttitudeHeaders(FileHeaders):
    """FITS headers for the attitude file (``trend/attitude/bc_csa_att.fits``)."""
    _header_templates = [AttitudePrimaryHeader(), AttitudeDataHeader()]


class DetectorHkHeaders(FileHeaders):
    """FITS headers for the detector housekeeping file
    (``auxil/bcYYMMDDcsa.hk.gz``).
    """
    _header_templates = [AuxPrimaryHeader(), DetectorHk1Header(),
                         DetectorHk2Header()]


class GtiHeaders(FileHeaders):
    """FITS headers for a standalone trend GTI file (``trend/gti_*/*.gti``)."""
    _header_templates = [GtiTrendPrimaryHeader(), GtiTrendDataHeader()]
