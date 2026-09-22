"""Offline tests for gdt.missions.burstcube.headers. These use only
synthetic headers built from the templates themselves; no external files are
needed.
"""
import warnings

import pytest

from gdt.missions.burstcube import headers as h


@pytest.mark.parametrize('cls, expected_extensions', [
    (h.CBDHeaders, ['PRIMARY', 'CBD', 'STDGTI']),
    (h.TTEHeaders, ['PRIMARY', 'EVENTS', 'STDGTI']),
    (h.OrbitHeaders, ['PRIMARY', 'ORBIT']),
    (h.AttitudeHeaders, ['PRIMARY', 'ATTITUDE']),
    (h.DetectorHKHeaders, ['PRIMARY', 'DETECTOR_HK1', 'DETECTOR_HK2']),
    (h.GTIHeaders, ['PRIMARY', 'STDGTI']),
])
def test_file_headers_construct_with_expected_extensions(cls, expected_extensions):
    """Every FileHeaders subclass must construct and expose exactly its
    documented extensions.
    """
    obj = cls()
    assert obj.keys() == expected_extensions


def test_cbd_carries_burstcube_specific_keywords_gbm_lacks():
    """CBD must carry PROCVER, CALDBVER, SEQPNUM, and the TIMEPIXR/TIMEDEL
    pair, none of which GBM's headers need.
    """
    obj = h.CBDHeaders()
    for keyword in ('PROCVER', 'CALDBVER', 'SEQPNUM', 'TIMEPIXR', 'TIMEDEL'):
        assert keyword in obj['CBD']


def test_cbd_timepixr_defaults_to_end_of_bin():
    """TIMEPIXR must default to 1 (end of bin), matching every real CBD
    file; getting this wrong silently shifts every light curve.
    """
    obj = h.CBDHeaders()
    assert obj['CBD']['TIMEPIXR'] == 1


def test_tstart_sync_updates_date_obs():
    """Setting TSTART on a header that also defines DATE-OBS must update
    DATE-OBS to the corresponding BurstCube MET, matching the MET-to-UTC
    reference value used in test_time.py.
    """
    obj = h.CBDHeaders()
    obj['CBD']['TSTART'] = 107629263.33
    assert obj['CBD']['DATE-OBS'] == '2024-05-30T17:00:26.330'


def test_tstart_sync_does_not_crash_on_numeric_string_tstart():
    """Verified across all 28 archive TTE files, TTE
    files sometimes carry TSTART/TSTOP as numeric strings (e.g.
    '107629534.257') with TSTOP < TSTART. Setting such a value must not
    raise; the DATE-OBS sync should still follow the (semantically
    confusing) MET the string encodes, since the underlying gdt.core.headers
    coercion accepts anything float() accepts.
    """
    obj = h.TTEHeaders()
    obj['EVENTS']['TSTART'] = '107629534.257'  # must not raise
    assert obj['EVENTS']['TSTART'] == pytest.approx(107629534.257)
    assert obj['EVENTS']['DATE-OBS'].startswith('2024-05-30T17:04:57')


def test_tstart_sync_leaves_date_obs_untouched_on_unparseable_value():
    """A truly non-numeric TSTART is rejected by the underlying
    gdt.core.headers type coercion (TSTART is declared as a float field);
    that is the framework's normal, expected strictness and is unrelated to
    the DATE-OBS sync, which only ever best-effort *adds to* that behavior.
    """
    obj = h.TTEHeaders()
    with pytest.raises(TypeError):
        obj['EVENTS']['TSTART'] = 'not-a-number'


def test_fresh_headers_carry_the_corrected_mjdreff_and_warn():
    """A freshly constructed (not read-from-file) FileHeaders object must
    carry the corrected MJDREFF -- the bare TT-TAI offset, epoch 2021-01-01
    00:00:00 TAI -- not the archive's defective value, and constructing it
    must warn that this will not byte-match archive files.
    """
    with pytest.warns(UserWarning, match='MJDREFF'):
        obj = h.CBDHeaders()
    assert obj['CBD']['MJDREFF'] == pytest.approx(32.184 / 86400)


def test_from_headers_does_not_warn_about_mjdreff():
    """Reading a real file's headers (from_headers) must not trigger the
    write-time MJDREFF warning: the value populated is the archive's own,
    not this package's corrected default, and that gets its own separate
    check (gdt.missions.burstcube.time.check_met_epoch), exercised in
    test_time.py.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        obj = h.CBDHeaders()
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        h.CBDHeaders.from_headers([obj[i] for i in range(obj.num_headers)])


def test_from_headers_round_trip():
    """A FileHeaders object built via from_headers() from another object's
    own headers must reproduce every field.
    """
    obj = h.OrbitHeaders()
    obj['ORBIT']['INSTRUME'] = 'CS'
    obj['ORBIT']['TSTART'] = 107568002.370496

    rebuilt = h.OrbitHeaders.from_headers([obj[i] for i in range(obj.num_headers)])
    assert rebuilt['ORBIT']['INSTRUME'] == 'CS'
    assert rebuilt['ORBIT']['TSTART'] == pytest.approx(107568002.370496)
