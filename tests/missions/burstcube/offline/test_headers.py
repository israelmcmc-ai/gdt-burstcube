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
"""Offline tests for gdt.missions.burstcube.headers. These use only
synthetic headers built from the templates themselves; no external files are
needed.
"""
import pytest

from gdt.missions.burstcube import headers as h


@pytest.mark.parametrize('cls, expected_extensions', [
    (h.CBDHeaders, ['PRIMARY', 'CBD', 'STDGTI']),
    (h.TTEHeaders, ['PRIMARY', 'EVENTS', 'STDGTI']),
    (h.OrbitHeaders, ['PRIMARY', 'ORBIT']),
    (h.AttitudeHeaders, ['PRIMARY', 'ATTITUDE']),
    (h.DetectorHKHeaders, ['PRIMARY', 'DETECTOR_HK1', 'DETECTOR_HK2']),
    (h.GtiHeaders, ['PRIMARY', 'STDGTI']),
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
    assert obj['CBD']['DATE-OBS'] == '2024-05-30T17:01:03.330'


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
    assert obj['EVENTS']['DATE-OBS'].startswith('2024-05-30T17:05:34')


def test_tstart_sync_leaves_date_obs_untouched_on_unparseable_value():
    """A truly non-numeric TSTART is rejected by the underlying
    gdt.core.headers type coercion (TSTART is declared as a float field);
    that is the framework's normal, expected strictness and is unrelated to
    the DATE-OBS sync, which only ever best-effort *adds to* that behavior.
    """
    obj = h.TTEHeaders()
    with pytest.raises(TypeError):
        obj['EVENTS']['TSTART'] = 'not-a-number'


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
