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
from contextlib import contextmanager

from gdt.core.headers import FileHeaders, Header

from .time import Time

__all__ = ['AttitudeHeaders', 'BurstCubeFileHeaders', 'CBDHeaders',
           'CBDUnfilteredHeaders', 'DetectorHKHeaders', 'GTIHeaders',
           'OrbitHeaders', 'RspHeaders', 'TTEHeaders']

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


#: Whether setting TSTART/TSTOP should also refresh DATE-OBS/DATE-END. Off
#: while headers are being populated from a file on disk -- see
#: :func:`_dates_unsynced`.
_SYNC_DATES = True


@contextmanager
def _dates_unsynced():
    """Suspend the ``TSTART``/``TSTOP`` -> ``DATE-OBS``/``DATE-END`` sync.

    Reading a file must not silently rewrite its own header values. The
    templates list ``DATE-OBS``/``DATE-END`` before ``TSTART``/``TSTOP``, so
    without this the copied date would be overwritten moments later by one
    recomputed from the timing keywords. That is a real divergence, not a
    cosmetic one: in ``bc240530cs0_tte_uf.evt.gz`` all three extensions carry
    the same DATE-OBS/DATE-END on disk, but ``STDGTI`` has a different (and
    saner) TSTART/TSTOP than ``EVENTS``, so resyncing while loading left the
    three extensions reporting three different observation times -- and
    writing the object back out would have persisted values the archive never
    contained.

    The sync itself is still wanted when *we* set a new time range, e.g. in
    ``_build_headers`` after a slice or rebin.
    """
    global _SYNC_DATES
    previous = _SYNC_DATES
    _SYNC_DATES = False
    try:
        yield
    finally:
        _SYNC_DATES = previous


class BurstCubeHeader(Header):
    """A :class:`~gdt.core.headers.Header` that keeps ``DATE-OBS``/``DATE-END``
    in sync with ``TSTART``/``TSTOP`` (interpreted as BurstCube MET) whenever a
    new time range is set on it.

    The sync is suspended while headers are populated from a file, so that
    reading never alters what the archive wrote; see :func:`_dates_unsynced`.
    In every archive TTE file ``TSTART``/``TSTOP`` are written as
    strings with ``TSTOP < TSTART`` in TTE files, so the sync is best-effort
    and quietly does nothing if the value cannot be interpreted as a MET,
    rather than raising.
    """

    def __setitem__(self, key, val):
        if _SYNC_DATES and not isinstance(key, tuple) and not isinstance(val, tuple):
            date_key = {'TSTART': 'DATE-OBS', 'TSTOP': 'DATE-END'}.get(key.upper())
            if date_key is not None and date_key in self:
                try:
                    self[date_key] = Time(float(val), format='burstcube').utc.isot
                except (TypeError, ValueError):
                    pass

        super().__setitem__(key, val)


class BurstCubeFileHeaders(FileHeaders):
    """A :class:`~gdt.core.headers.FileHeaders` that preserves the date
    keywords a file actually carries.

    Populating from disk happens inside :func:`_dates_unsynced`, so a header
    read back from a file matches the file, byte for byte, even where the
    archive's own values are internally inconsistent.
    """

    @classmethod
    def from_headers(cls, headers):
        """Build from a list of headers read from a file, leaving their
        ``DATE-OBS``/``DATE-END`` values untouched.

        Args:
            headers (list of :class:`astropy.io.fits.Header`): The headers

        Returns:
            (:class:`BurstCubeFileHeaders`)
        """
        with _dates_unsynced():
            return super().from_headers(headers)


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


class CBDDataHeader(BurstCubeHeader):
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


class CBDGTIHeader(BurstCubeHeader):
    """Header for the ``STDGTI`` extension of a *cleaned* (``_cl``) CBD file:
    the 2-column ``START``/``STOP`` OGIP standard form, carrying
    ``HDUCLAS2='STANDARD'``, ``HDUVERS`` and ``TIMEZERO``.

    Note:
        The two CBD variants do not share a ``STDGTI`` schema. Unfiltered
        (``_uf``) files instead carry the 6-column form used by TTE
        (``START``, ``STOP``, ``START_ORIGINAL_TIME``,
        ``START_TIME_SYST_ERROR``, ``STOP_ORIGINAL_TIME``,
        ``STOP_TIME_SYST_ERROR``) with ``HDUCLAS1='GTI'`` and none of the
        three keywords above -- see :class:`TTEGTIHeader`. Because nothing in
        the files guarantees this correlation holds for every observation,
        :meth:`~gdt.missions.burstcube.cbd.BurstCubeCBD.open` selects between
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

    In every archive TTE file ``TSTART``/``TSTOP`` are stored as
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


class TTEGTIHeader(BurstCubeHeader):
    """Header for the ``STDGTI`` extension of a TTE file. Unlike
    :class:`CBDGTIHeader`, its columns are ``START``, ``STOP``,
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


class DetectorHKHeader(BurstCubeHeader):
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


class DetectorHK1Header(DetectorHKHeader):
    """Header for the ``DETECTOR_HK1`` extension (1 s cadence)."""
    name = 'DETECTOR_HK1'
    keywords = DetectorHKHeader.keywords.copy()


class DetectorHK2Header(DetectorHKHeader):
    """Header for the ``DETECTOR_HK2`` extension (60 s cadence)."""
    name = 'DETECTOR_HK2'
    keywords = [kw if kw[0] != 'TIMEDEL' else ('TIMEDEL', 60.0, 'Integration time')
               for kw in DetectorHKHeader.keywords]


class GTITrendPrimaryHeader(BurstCubeHeader):
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
        :class:`GTITrendDataHeader`.
    """
    name = 'PRIMARY'
    keywords = [_telescope_card, _instrument_card, _origin_card,
                Header.creator()]


class GTITrendDataHeader(BurstCubeHeader):
    """Header for the ``STDGTI`` extension of a standalone trend GTI file
    (``trend/gti_binning``, ``trend/gti_poscnt``, ``trend/gti_saa``).

    Note:
        ``HDUVERS`` and ``TIMEZERO`` appear on some trend GTI variants (e.g.
        ``trend/gti_saa``, matching :class:`CBDGTIHeader`) but not others
        (e.g. ``trend/gti_binning``), so, unlike :class:`CBDGTIHeader`, they
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

class CBDHeaders(BurstCubeFileHeaders):
    """FITS headers for a CBD file whose ``STDGTI`` uses the 2-column
    standard schema (in practice, the cleaned ``_cl`` files)."""
    _header_templates = [DataPrimaryHeader(), CBDDataHeader(), CBDGTIHeader()]


class CBDUnfilteredHeaders(BurstCubeFileHeaders):
    """FITS headers for a CBD file whose ``STDGTI`` uses the 6-column schema
    shared with TTE (in practice, the unfiltered ``_uf`` files).

    Identical to :class:`CBDHeaders` apart from the GTI extension. See the
    note on :class:`CBDGTIHeader`.
    """
    _header_templates = [DataPrimaryHeader(), CBDDataHeader(), TTEGTIHeader()]


class TTEHeaders(BurstCubeFileHeaders):
    """FITS headers for TTE (time-tagged event) files."""
    _header_templates = [DataPrimaryHeader(), EventsHeader(), TTEGTIHeader()]


class OrbitHeaders(BurstCubeFileHeaders):
    """FITS headers for the orbit/ephemeris file (``auxil/bcYYMMDD.hk.gz``)."""
    _header_templates = [AuxPrimaryHeader(), OrbitDataHeader()]


class AttitudeHeaders(BurstCubeFileHeaders):
    """FITS headers for the attitude file (``trend/attitude/bc_csa_att.fits``)."""
    _header_templates = [AttitudePrimaryHeader(), AttitudeDataHeader()]


class DetectorHKHeaders(BurstCubeFileHeaders):
    """FITS headers for the detector housekeeping file
    (``auxil/bcYYMMDDcsa.hk.gz``).
    """
    _header_templates = [AuxPrimaryHeader(), DetectorHK1Header(),
                         DetectorHK2Header()]


class GTIHeaders(BurstCubeFileHeaders):
    """FITS headers for a standalone trend GTI file (``trend/gti_*/*.gti``)."""
    _header_templates = [GTITrendPrimaryHeader(), GTITrendDataHeader()]


# ------------------------------------------------------------------------
# Response (.rsp) headers. These use a different keyword vocabulary from
# the science-data products above (CALDB CCLS0001-style boilerplate rather
# than OBS_ID/PROCVER/etc.), so they get their own small set of cards.

_filter_card = ('FILTER', 'NONE', 'Filter')
_ordering_card = ('ORDERING', 'RING', 'Pixel ordering scheme, either RING or NESTED')
_pixel_card = ('PIXEL', 0, 'Pixel number for HEAPIX representation')
_chantype_card = ('CHANTYPE', 'PHA', 'PHA or PI channel')
_detchans_card = ('DETCHANS', 64, 'Detector channel')

# the CALDB "boilerplate" block (CCLS0001 etc.) common to both extensions
_rsp_caldb_cards = [
    ('CCLS0001', 'CPF', 'Dataset is Basic Calibration File'),
    ('CDTP0001', 'DATA', 'Calibration file contains data'),
    ('CVSD0001', '2020-01-01', 'UTC date when calibration should first be used'),
    ('CVST0001', '00:01:00', 'UTC time when calibration should first be used'),
    ('CBD10001', 'DATAMODE(EVENT,CBD,ATD)', 'Applicable modes'),
    ('CBD20001', '', 'Pixel number for HEAPIX representation'),
    ('CBD30001', '', 'Distance from optical axis'),
    ('CBD40001', '', 'Azimuthal angle'),
    _hduclass_card, ('HDUCLAS1', 'RESPONSE', 'hduclass1'),
    ('HDUVERS', '1.2.0', 'Version of format (OGIP memo OGIP-92-007)'),
    ('HDUVERS1', '1.2.0', 'Obsolete included for back compatibility'),
]


class RspPrimaryHeader(BurstCubeHeader):
    """PRIMARY header for a ``.rsp`` file. Much sparser than the science-data
    PRIMARY headers: no ``OBS_ID``, no ``DATAMODE``, no ``DATE-OBS``/``DATE-END``.
    """
    name = 'PRIMARY'
    keywords = [_telescope_card, _instrument_card, _origin_card, _date_card,
               Header.creator()]


class RspEboundsHeader(BurstCubeHeader):
    """Header for the ``EBOUNDS`` extension of a ``.rsp`` file: the (older,
    rounded, detector-independent) 64-channel energy grid the response
    itself was computed on -- see :meth:`~gdt.missions.burstcube.response.BurstCubeRsp.to_cbd`
    for why this must not be confused with CALDB's own ``eb16``/``eb64``.
    """
    name = 'EBOUNDS'
    keywords = [_extname_card, _chantype_card, _detchans_card,
               _telescope_card, _instrument_card, _filter_card,
               _ordering_card, _pixel_card, ('CCNM0001', 'EBOUNDS',
               'Type of calibration data')] + _rsp_caldb_cards + [
               ('HDUCLAS2', 'EBOUNDS', 'hduclas2'),
               ('CDES0001', '', 'Description keyword')] + [_origin_card,
               _date_card]


class RspSpecrespHeader(BurstCubeHeader):
    """Header for the ``SPECRESP MATRIX`` extension of a ``.rsp`` file: the
    response matrix itself, for one HEALPix pixel (``PIXEL``) of one
    detector (``INSTRUME``).
    """
    name = 'SPECRESP MATRIX'
    keywords = [_extname_card, _chantype_card, _detchans_card,
               _telescope_card, _instrument_card, _filter_card,
               _ordering_card, _pixel_card, ('CCNM0001', 'MATRIX',
               'Type of calibration data')] + _rsp_caldb_cards + [
               ('HDUCLAS2', 'RSP_MATRIX', 'hduclas2'),
               ('HDUCLAS3', 'FULL', 'hduclas3'),
               ('CDES0001', '', 'Description keyword')] + [_origin_card,
               _date_card]


class RspHeaders(BurstCubeFileHeaders):
    """FITS headers for a single-DRM BurstCube response file (``.rsp``)."""
    _header_templates = [RspPrimaryHeader(), RspEboundsHeader(),
                         RspSpecrespHeader()]
