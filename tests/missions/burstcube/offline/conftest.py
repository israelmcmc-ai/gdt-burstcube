"""Synthetic FITS fixture builders shared across the round-2 offline tests.
Each function reproduces the real header structure verified against sample
archive files (see the plugin spec), so that opening these files exercises
the same code paths a real file would, without needing network access.
"""
import numpy as np
import astropy.io.fits as fits

TIMEDEL_CBD = 0.256


def _common_time_cards(instrume, datamode, obs_id, tstart, tstop, timedel=None,
                       clockapp=False):
    cards = dict(TELESCOP='BURSTCUBE', INSTRUME=instrume, DATAMODE=datamode,
                OBS_ID=obs_id, MJDREFI=59215, MJDREFF=0.00080074074074074,
                TIMEREF='LOCAL', TASSIGN='SATELLITE', TIMESYS='TT',
                TIMEUNIT='s', CLOCKAPP=clockapp, TSTART=tstart, TSTOP=tstop,
                ORIGIN='GSFC', CREATOR='BCTOOLS', PROCVER='00.00.00.00',
                CALDBVER='20250226', SEQPNUM=1, DATE='2025-07-17T07:21:09')
    if timedel is not None:
        cards['TIMEDEL'] = timedel
    return cards


def make_cbd_fits(path, detector='CS0', time=None, counts=None, gti=None,
                  gti_schema='cl', with_sums=None):
    """Build a synthetic CBD FITS file (PRIMARY, CBD, STDGTI).

    Args:
        path (Path): Where to write the file
        detector (str): The INSTRUME value
        time (numpy.ndarray, optional): The end-of-bin TIME values (TIMEPIXR=1).
            Defaults to 20 bins of exactly 0.256 s starting at MET 107629263.33.
        counts (numpy.ndarray, optional): (n, 16) counts array. Defaults to
            small deterministic values.
        gti (list of (float, float), optional): GTI intervals. Defaults to
            the full time range as one interval.
        gti_schema (str): 'cl' for the 2-column START/STOP schema, or 'uf'
            for the 6-column schema, matching the real files.
        with_sums (bool, optional): Whether to write the SUMTOT and
            SUMCH_2_15 columns. Real files carry them only in the cleaned
            variant. Defaults to True for `gti_schema='cl'` and False for
            `'uf'`, mirroring the archive.

    Returns:
        (numpy.ndarray, numpy.ndarray): the `time` and `counts` arrays used,
        for the caller to compare against.
    """
    if time is None:
        n = 20
        time = 107629263.33 + (np.arange(n) + 1) * TIMEDEL_CBD
    else:
        n = len(time)
    if counts is None:
        rng = np.random.default_rng(0)
        counts = rng.integers(0, 5, size=(n, 16)).astype(np.int16)
    if gti is None:
        gti = [(time[0] - TIMEDEL_CBD, time[-1])]

    original_time = time - 0.01
    time_syst_error = np.full(n, 0.05)
    sumtot = counts.sum(axis=1).astype(np.int32)
    sumch_2_15 = counts[:, 2:16].sum(axis=1).astype(np.int32)

    tstart, tstop = time[0] - TIMEDEL_CBD, time[-1]

    primary_hdr = fits.Header()
    for k, v in _common_time_cards(detector, 'CBD', '240530', tstart, tstop).items():
        if k in ('MJDREFI', 'MJDREFF', 'TIMEREF', 'TASSIGN', 'TIMESYS',
                 'TIMEUNIT', 'TSTART', 'TSTOP', 'PROCVER', 'CALDBVER',
                 'SEQPNUM', 'CLOCKAPP'):
            continue
        primary_hdr[k] = v
    primary_hdr['DATE-OBS'] = '2024-05-30T17:01:03.330'
    primary_hdr['DATE-END'] = '2024-05-30T18:34:42.719312'

    if with_sums is None:
        with_sums = (gti_schema == 'cl')

    cbd_cols = [
        fits.Column(name='TIME', format='1D', unit='s', array=time),
        fits.Column(name='COUNTS', format='16I', unit='count', array=counts),
        fits.Column(name='ORIGINAL_TIME', format='D', unit='s', array=original_time),
        fits.Column(name='TIME_SYST_ERROR', format='D', unit='s', array=time_syst_error),
    ]
    if with_sums:
        cbd_cols += [
            fits.Column(name='SUMTOT', format='1J', unit='count', array=sumtot),
            fits.Column(name='SUMCH_2_15', format='1J', unit='count',
                        array=sumch_2_15),
        ]
    cbd_hdu = fits.BinTableHDU.from_columns(cbd_cols, name='CBD')
    for k, v in _common_time_cards(detector, 'CBD', '240530', tstart, tstop,
                                   timedel=TIMEDEL_CBD, clockapp=True).items():
        cbd_hdu.header[k] = v
    cbd_hdu.header['TIMEPIXR'] = 1
    cbd_hdu.header['DATE-OBS'] = '2024-05-30T17:01:03.330'
    cbd_hdu.header['DATE-END'] = '2024-05-30T18:34:42.719312'
    cbd_hdu.header['DEADAPP'] = False
    cbd_hdu.header['TELAPSE'] = tstop - tstart
    cbd_hdu.header['ONTIME'] = sum(hi - lo for lo, hi in gti)
    cbd_hdu.header['HDUCLASS'] = 'OGIP'
    cbd_hdu.header['HDUCLAS1'] = 'ARRAY'
    cbd_hdu.header['HDUCLAS2'] = 'TOTAL'
    cbd_hdu.header['1CDLT2'] = 1
    cbd_hdu.header['1CDPX2'] = 1
    cbd_hdu.header['1CRVL2'] = 0
    cbd_hdu.header['1CUNI2'] = 'chan'
    cbd_hdu.header['1CTYP2'] = 'CHANNEL'

    gti_start = np.array([lo for lo, hi in gti])
    gti_stop = np.array([hi for lo, hi in gti])
    if gti_schema == 'cl':
        gti_cols = [
            fits.Column(name='START', format='1D', unit='s', array=gti_start),
            fits.Column(name='STOP', format='1D', unit='s', array=gti_stop),
        ]
    else:
        gti_cols = [
            fits.Column(name='START', format='1D', unit='s', array=gti_start),
            fits.Column(name='STOP', format='1D', unit='s', array=gti_stop),
            fits.Column(name='START_ORIGINAL_TIME', format='D', unit='s', array=gti_start),
            fits.Column(name='START_TIME_SYST_ERROR', format='D', unit='s',
                       array=np.full(gti_start.size, 0.05)),
            fits.Column(name='STOP_ORIGINAL_TIME', format='D', unit='s', array=gti_stop),
            fits.Column(name='STOP_TIME_SYST_ERROR', format='D', unit='s',
                       array=np.full(gti_stop.size, 0.05)),
        ]
    gti_hdu = fits.BinTableHDU.from_columns(gti_cols, name='STDGTI')
    for k, v in _common_time_cards(detector, 'CBD', '240530', gti_start[0],
                                   gti_stop[-1], clockapp=True).items():
        gti_hdu.header[k] = v
    gti_hdu.header['HDUCLASS'] = 'OGIP'
    gti_hdu.header['HDUCLAS1'] = 'GTI'
    if gti_schema == 'cl':
        gti_hdu.header['HDUCLAS2'] = 'STANDARD'
        gti_hdu.header['HDUVERS'] = '1.0.0'
        gti_hdu.header['TIMEZERO'] = 0.0
    gti_hdu.header['DATE-OBS'] = '2024-05-30T17:01:03.330'
    gti_hdu.header['DATE-END'] = '2024-05-30T18:34:42.719312'

    hdulist = fits.HDUList([fits.PrimaryHDU(header=primary_hdr), cbd_hdu, gti_hdu])
    hdulist.writeto(path, overwrite=True)
    return time, counts


def make_tte_fits(path, detector='CS0', times=None, channels=None,
                  broken_tstart_tstop=True):
    """Build a synthetic TTE FITS file (PRIMARY, EVENTS, STDGTI).

    Args:
        path (Path): Where to write the file
        detector (str): The INSTRUME value
        times (numpy.ndarray, optional): Event times. Defaults to 500
            randomly ordered events over ~5 s.
        channels (numpy.ndarray, optional): Event PHA channels (0-1023).
        broken_tstart_tstop (bool): If True (the default, matching the real
            archive), EVENTS' TSTART/TSTOP are written as strings with
            TSTOP < TSTART, as 16 of the 28 archive TTE files do.

    Returns:
        (numpy.ndarray, numpy.ndarray): the `times` and `channels` arrays used.
    """
    if times is None:
        rng = np.random.default_rng(1)
        times = np.sort(rng.uniform(107629263.33, 107629263.33 + 20 * TIMEDEL_CBD,
                                    size=500))
    if channels is None:
        rng = np.random.default_rng(2)
        channels = rng.integers(0, 1024, size=len(times)).astype(np.int16)

    primary_hdr = fits.Header()
    primary_hdr['TELESCOP'] = 'BURSTCUBE'
    primary_hdr['INSTRUME'] = detector
    primary_hdr['DATAMODE'] = 'EVENT'
    primary_hdr['OBS_ID'] = '240530'
    primary_hdr['DATE-OBS'] = '2024-05-30T17:05:34'
    primary_hdr['DATE-END'] = '2024-05-30T17:05:18'
    primary_hdr['ORIGIN'] = 'GSFC'
    primary_hdr['DATE'] = '2025-07-17T07:24:18'
    primary_hdr['CREATOR'] = 'BCTOOLS'

    ev_cols = [
        fits.Column(name='TIME', format='1D', unit='s', array=times),
        fits.Column(name='PHA', format='1I', unit='chan', array=channels),
        fits.Column(name='ORIGINAL_TIME', format='D', unit='s', array=times - 0.01),
        fits.Column(name='TIME_SYST_ERROR', format='D', unit='s',
                   array=np.full(len(times), 0.05)),
    ]
    ev_hdu = fits.BinTableHDU.from_columns(ev_cols, name='EVENTS')
    if broken_tstart_tstop:
        ev_tstart, ev_tstop = '107629534.257', '107629518.648'
    else:
        ev_tstart, ev_tstop = float(times[0]), float(times[-1])
    for k, v in dict(TELESCOP='BURSTCUBE', INSTRUME=detector, DATAMODE='EVENT',
                     OBS_ID='240530', OBJECT='UNKNOWN', RA_OBJ=-2000.0,
                     DEC_OBJ=-2000.0, POSFLAG=0, EQUINOX=2000.0, RADECSYS='FK5',
                     MJDREFI=59215, MJDREFF=0.00080074074074074, TIMEREF='LOCAL',
                     TASSIGN='SATELLITE', TIMESYS='TT', TIMEUNIT='s', TIMEDEL=1e-5,
                     CLOCKAPP=False, **{'DATE-OBS': '2024-05-30T17:05:34',
                                       'DATE-END': '2024-05-30T17:05:18'},
                     TSTART=ev_tstart, TSTOP=ev_tstop, TELAPSE=-15.6, ONTIME=295.07,
                     EXPOSURE=-15.6, DEADAPP=False, EVTTYPE='T3E', TRIG_ID='',
                     TRIGGER=0, TRIGTIME=0.0, TRIGUTC=0.0, TRIGDET='', TRIGSCAL=0,
                     TRIGSIGN=0.0, TRIGEMIN=0.0, TRIGEMAX=0.0, HDUCLASS='OGIP',
                     HDUCLAS1='EVENT', HDUCLAS2='ALL', ORIGIN='GSFC',
                     CREATOR='BCTOOLS', PROCVER='00.00.00.00', CALDBVER='20250226',
                     SEQPNUM=1, DATE='2025-07-17T07:24:18').items():
        ev_hdu.header[k] = v

    gti_start = np.array([times[0]])
    gti_stop = np.array([times[-1]])
    gti_cols = [
        fits.Column(name='START', format='1D', unit='s', array=gti_start),
        fits.Column(name='STOP', format='1D', unit='s', array=gti_stop),
        fits.Column(name='START_ORIGINAL_TIME', format='D', unit='s', array=gti_start),
        fits.Column(name='START_TIME_SYST_ERROR', format='D', unit='s', array=[0.05]),
        fits.Column(name='STOP_ORIGINAL_TIME', format='D', unit='s', array=gti_stop),
        fits.Column(name='STOP_TIME_SYST_ERROR', format='D', unit='s', array=[0.05]),
    ]
    gti_hdu = fits.BinTableHDU.from_columns(gti_cols, name='STDGTI')
    for k, v in dict(TELESCOP='BURSTCUBE', INSTRUME=detector, DATAMODE='EVENT',
                     OBS_ID='240530', MJDREFI=59215, MJDREFF=0.00080074074074074,
                     TIMEREF='LOCAL', TASSIGN='SATELLITE', TIMESYS='TT',
                     TIMEUNIT='s', CLOCKAPP=False,
                     **{'DATE-OBS': '2024-05-30T17:05:34',
                       'DATE-END': '2024-05-30T17:05:18'},
                     TSTART=float(gti_start[0]), TSTOP=float(gti_stop[0]),
                     HDUCLASS='OGIP', HDUCLAS1='GTI', CREATOR='BCTOOLS',
                     PROCVER='00.00.00.00', CALDBVER='20250226', SEQPNUM=1,
                     ORIGIN='GSFC', DATE='2025-07-17T07:24:18').items():
        gti_hdu.header[k] = v

    hdulist = fits.HDUList([fits.PrimaryHDU(header=primary_hdr), ev_hdu, gti_hdu])
    hdulist.writeto(path, overwrite=True)
    return times, channels


def make_orbit_fits(path, n=8):
    """Build a synthetic orbit FITS file (PRIMARY, ORBIT)."""
    time = 107568002.370496 + np.arange(n) * 30.0
    x = np.linspace(7000, 7010, n)
    y = np.linspace(0, 10, n)
    z = np.linspace(-100, -90, n)
    vx = np.full(n, 1.0)
    vy = np.full(n, 2.0)
    vz = np.full(n, 3.0)

    primary_hdr = fits.Header()
    primary_hdr['TELESCOP'] = 'BURSTCUBE'
    primary_hdr['INSTRUME'] = 'CS'
    primary_hdr['OBS_ID'] = '240530'
    primary_hdr['DATE-OBS'] = '2024-05-30T00:00:02.370'
    primary_hdr['DATE-END'] = '2024-05-30T23:59:32.924'
    primary_hdr['ORIGIN'] = 'GSFC'
    primary_hdr['DATE'] = '2025-07-17T13:42:08'
    primary_hdr['CREATOR'] = 'BCTOOLS'

    cols = [
        fits.Column(name='TIME', format='1D', unit='s', array=time),
        fits.Column(name='X', format='1D', unit='km', array=x),
        fits.Column(name='Y', format='1D', unit='km', array=y),
        fits.Column(name='Z', format='1D', unit='km', array=z),
        fits.Column(name='Vx', format='1D', unit='km/s', array=vx),
        fits.Column(name='Vy', format='1D', unit='km/s', array=vy),
        fits.Column(name='Vz', format='1D', unit='km/s', array=vz),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name='ORBIT')
    for k, v in dict(TELESCOP='BURSTCUBE', INSTRUME='CS', OBS_ID='240530',
                     MJDREFI=59215, MJDREFF=0.00080074074074074, TIMEREF='LOCAL',
                     TASSIGN='SATELLITE', TIMESYS='TT', TIMEUNIT='s', CLOCKAPP=False,
                     **{'DATE-OBS': '2024-05-30T00:00:02.370',
                       'DATE-END': '2024-05-30T23:59:32.924'},
                     TSTART=float(time[0]), TSTOP=float(time[-1]), HDUCLASS='OGIP',
                     HDUCLAS1='TEMPORALDATA', HDUCLAS2='EPHEM', CREATOR='BCTOOLS',
                     PROCVER='00.00.00.00', CALDBVER='20250226', SEQPNUM=1,
                     ORIGIN='GSFC').items():
        hdu.header[k] = v

    hdulist = fits.HDUList([fits.PrimaryHDU(header=primary_hdr), hdu])
    hdulist.writeto(path, overwrite=True)
    return time


def make_attitude_fits(path):
    """Build a synthetic attitude FITS file with the 3 real rows verified in
    the plugin spec.
    """
    time = np.array([110220808.0, 111192020.0, 111465000.0])
    qparam = np.array([
        [0.0844, 0.6314, 0.6416, 0.4273],
        [0.1083, 0.3291, 0.1500, 0.9260],
        [-0.1891, 0.9078, 0.0591, -0.3696],
    ])
    pointing = np.array([
        [48.723068, 10.861862, 333.95041],
        [350.98593, 49.458153, 0.0],
        [182.68195, -46.037468, 0.0],
    ])

    primary_hdr = fits.Header()
    primary_hdr['TELESCOP'] = 'BURSTCUBE'
    primary_hdr['INSTRUME'] = 'CSA'
    primary_hdr['DATE-OBS'] = '2024-06-29T16:53:28'
    primary_hdr['DATE-END'] = '2024-07-14T02:30:00'
    primary_hdr['ORIGIN'] = 'GSFC'
    primary_hdr['CREATOR'] = 'BCTOOLS'

    cols = [
        fits.Column(name='TIME', format='1D', unit='s', array=time),
        fits.Column(name='QPARAM', format='4D', array=qparam),
        fits.Column(name='POINTING', format='3D', unit='deg', array=pointing),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name='ATTITUDE')
    for k, v in dict(TELESCOP='BURSTCUBE', INSTRUME='CSA', MJDREFI=59215,
                     MJDREFF=0.00080074074074074, TIMEREF='LOCAL', TASSIGN='SATELLITE',
                     TIMESYS='TT', TIMEUNIT='s', CLOCKAPP=False,
                     **{'DATE-OBS': '2024-06-29T16:53:28',
                       'DATE-END': '2024-07-14T02:30:00'},
                     TSTART=float(time[0]), TSTOP=float(time[-1]), HDUCLASS='OGIP',
                     HDUCLAS1='TEMPORALDATA', HDUCLAS2='ATTITUDE', ORIGIN='GSFC',
                     CREATOR='BCTOOLS', PROCVER='00.00.00.00', CALDBVER='20250226',
                     SEQPNUM=0).items():
        hdu.header[k] = v

    hdulist = fits.HDUList([fits.PrimaryHDU(header=primary_hdr), hdu])
    hdulist.writeto(path, overwrite=True)
    return time, qparam, pointing


def make_hk_fits(path, n1=10, n2=3):
    """Build a synthetic detector housekeeping FITS file (PRIMARY,
    DETECTOR_HK1, DETECTOR_HK2).
    """
    time1 = 107629263.33 + np.arange(n1) * 1.0
    time2 = 107629263.33 + np.arange(n2) * 60.0

    trig_enabled = np.zeros((n1, 4), dtype=np.int16)
    trig_enabled[:, 1] = 1
    det_enabled = np.ones((n1, 4), dtype=np.int16)
    cbd_enabled = np.ones((n1, 4), dtype=np.int16)
    tte_enabled = np.zeros((n1, 4), dtype=np.int16)

    peak_thres = np.tile([26.93, 30.33, 21.40, 20.30], (n2, 1)).astype(np.float32)
    base_thres = np.tile([10.0, 10.0, 10.0, 10.0], (n2, 1)).astype(np.float32)

    primary_hdr = fits.Header()
    primary_hdr['TELESCOP'] = 'BURSTCUBE'
    primary_hdr['INSTRUME'] = 'CSA'
    primary_hdr['OBS_ID'] = '240530'
    primary_hdr['DATE-OBS'] = '2024-05-30T00:00:17'
    primary_hdr['DATE-END'] = '2024-05-30T23:59:51'
    primary_hdr['ORIGIN'] = 'GSFC'
    primary_hdr['DATE'] = '2025-07-17T07:20:45'
    primary_hdr['CREATOR'] = 'BCTOOLS'

    hk1_cols = [
        fits.Column(name='TIME', format='1D', unit='s', array=time1),
        fits.Column(name='TRIG_ENABLED', format='4I', array=trig_enabled),
        fits.Column(name='DET_ENABLED', format='4I', array=det_enabled),
        fits.Column(name='CBD_ENABLED', format='4I', array=cbd_enabled),
        fits.Column(name='TTE_ENABLED', format='4I', array=tte_enabled),
    ]
    hk1_hdu = fits.BinTableHDU.from_columns(hk1_cols, name='DETECTOR_HK1')
    for k, v in dict(TELESCOP='BURSTCUBE', INSTRUME='CSA', OBS_ID='240530',
                     MJDREFI=59215, MJDREFF=0.00080074074074074, TIMEREF='LOCAL',
                     TASSIGN='SATELLITE', TIMESYS='TT', TIMEUNIT='s', TIMEDEL=1.0,
                     CLOCKAPP=False, **{'DATE-OBS': '2024-05-30T00:00:17',
                                       'DATE-END': '2024-05-30T23:59:51'},
                     TSTART=float(time1[0]), TSTOP=float(time1[-1]), HDUCLASS='OGIP',
                     HDUCLAS1='TEMPORALDATA', HDUCLAS2='HKP', ORIGIN='GSFC',
                     CREATOR='BCTOOLS', PROCVER='00.00.00.00', CALDBVER='20250226',
                     SEQPNUM=1, DATE='2025-07-17T07:20:45').items():
        hk1_hdu.header[k] = v

    hk2_cols = [
        fits.Column(name='TIME', format='1D', unit='s', array=time2),
        fits.Column(name='PEAK_THRES', format='4E', unit='mV', array=peak_thres),
        fits.Column(name='BASE_THRES', format='4E', unit='mV', array=base_thres),
    ]
    hk2_hdu = fits.BinTableHDU.from_columns(hk2_cols, name='DETECTOR_HK2')
    for k, v in dict(TELESCOP='BURSTCUBE', INSTRUME='CSA', OBS_ID='240530',
                     MJDREFI=59215, MJDREFF=0.00080074074074074, TIMEREF='LOCAL',
                     TASSIGN='SATELLITE', TIMESYS='TT', TIMEUNIT='s', TIMEDEL=60.0,
                     CLOCKAPP=False, **{'DATE-OBS': '2024-05-30T00:00:17',
                                       'DATE-END': '2024-05-30T23:59:51'},
                     TSTART=float(time2[0]), TSTOP=float(time2[-1]), HDUCLASS='OGIP',
                     HDUCLAS1='TEMPORALDATA', HDUCLAS2='HKP', ORIGIN='GSFC',
                     CREATOR='BCTOOLS', PROCVER='00.00.00.00', CALDBVER='20250226',
                     SEQPNUM=1, DATE='2025-07-17T07:20:45').items():
        hk2_hdu.header[k] = v

    hdulist = fits.HDUList([fits.PrimaryHDU(header=primary_hdr), hk1_hdu, hk2_hdu])
    hdulist.writeto(path, overwrite=True)
    return time1, time2


def make_trend_gti_fits(path, instrume='CS0', starts=None, stops=None):
    """Build a synthetic standalone trend GTI file (PRIMARY, STDGTI)."""
    if starts is None:
        starts = np.array([107629263.33, 107629366.146])
        stops = np.array([107629266.146, 107629368.45])

    primary_hdr = fits.Header()
    primary_hdr['TELESCOP'] = 'BURSTCUBE'
    primary_hdr['INSTRUME'] = instrume
    primary_hdr['ORIGIN'] = 'GSFC'
    primary_hdr['CREATOR'] = 'BCTOOLS'

    cols = [
        fits.Column(name='START', format='1D', unit='s', array=starts),
        fits.Column(name='STOP', format='1D', unit='s', array=stops),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name='STDGTI')
    for k, v in dict(TELESCOP='BURSTCUBE', INSTRUME=instrume, DATAMODE='CBD',
                     MJDREFI=59215, MJDREFF=0.00080074074074074, TIMEREF='LOCAL',
                     TASSIGN='SATELLITE', TIMESYS='TT', TIMEUNIT='s', CLOCKAPP=True,
                     **{'DATE-OBS': '2024-05-30T17:01:03.330',
                       'DATE-END': '2024-05-30T18:34:42.719312'},
                     TSTART=float(starts[0]), TSTOP=float(stops[-1]), HDUCLASS='OGIP',
                     HDUCLAS1='GTI', HDUCLAS2='STANDARD', ORIGIN='GSFC',
                     CREATOR='BCTOOLS', PROCVER='00.00.00.00', CALDBVER='20250226',
                     SEQPNUM=1, DATE='2025-07-17T07:21:09').items():
        hdu.header[k] = v

    hdulist = fits.HDUList([fits.PrimaryHDU(header=primary_hdr), hdu])
    hdulist.writeto(path, overwrite=True)
    return starts, stops


def make_catalog_fits(path, obsid, num_events, exposure_by_detector):
    """Build a synthetic burcbmastr-like catalog FITS table.

    Args:
        path (Path): Where to write the file
        obsid (list of str): Observation IDs
        num_events (list of int): Number of TTE files per row
        exposure_by_detector (dict): e.g. {'cs0': [...], 'cs1': [...], ...},
            with NaN for missing values
    """
    cols = [
        fits.Column(name='obsid', format='6A', array=np.array(obsid)),
        fits.Column(name='num_events', format='1J',
                   array=np.array(num_events, dtype=np.int32)),
    ]
    for det, values in exposure_by_detector.items():
        cols.append(fits.Column(name=f'exposure_{det}', format='1D',
                                array=np.array(values, dtype=np.float64)))
    hdu = fits.BinTableHDU.from_columns(cols)
    hdulist = fits.HDUList([fits.PrimaryHDU(), hdu])
    hdulist.writeto(path, overwrite=True)
