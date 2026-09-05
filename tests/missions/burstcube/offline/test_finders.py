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
"""Offline tests for gdt.missions.burstcube.finders. No live network: `cd()`
is never called (it would require a real directory listing), so the finder's
`_args` is set directly to simulate having already changed to a day, and the
protocol-level download is monkeypatched.
"""
import warnings
from pathlib import Path

import pytest
from astropy.time import Time

from gdt.missions.burstcube.finders import (BurstCubeObsFinder,
                                            BurstCubeTrendFinder,
                                            _detector_numbers)


def test_construct_path_from_string():
    finder = BurstCubeObsFinder()
    assert finder._construct_path('240530') == 'burstcube/data/obs/2024_05/240530'


def test_construct_path_from_time():
    finder = BurstCubeObsFinder()
    t = Time('240530', format='burstcube_obsid')
    assert finder._construct_path(t) == 'burstcube/data/obs/2024_05/240530'


def test_construct_path_handles_a_day_at_month_end():
    finder = BurstCubeObsFinder()
    assert finder._construct_path('240831') == 'burstcube/data/obs/2024_08/240831'


def test_trend_finder_path_is_fixed():
    finder = BurstCubeTrendFinder()
    assert finder._construct_path() == 'burstcube/data/trend'


@pytest.mark.parametrize('value, expected', [
    (None, [0, 1, 2, 3]),
    ('CS2', [2]),
    (2, [2]),
    ([0, 'CS1'], [0, 1]),
])
def test_detector_numbers_normalization(value, expected):
    assert _detector_numbers(value) == expected


def _finder_at(obs_id):
    """Build a finder as if cd() had already succeeded, without any network
    access.
    """
    finder = BurstCubeObsFinder()
    finder._args = (obs_id,)
    return finder


def test_get_cbd_requests_the_expected_relative_paths(monkeypatch):
    finder = _finder_at('240530')
    requested = []
    monkeypatch.setattr(finder._protocol, 'download',
                        lambda file, dest, verbose: requested.append(file) or Path(dest) / file)

    finder.get_cbd('/tmp/out', detectors=['CS0', 'CS1'], variant='cl', verbose=False)
    assert requested == ['monitor/bc240530cs0_3cbd_cl.fits.gz',
                        'monitor/bc240530cs1_3cbd_cl.fits.gz']


def test_get_tte_requests_the_expected_relative_paths(monkeypatch):
    finder = _finder_at('240530')
    requested = []
    monkeypatch.setattr(finder._protocol, 'download',
                        lambda file, dest, verbose: requested.append(file) or Path(dest) / file)

    finder.get_tte('/tmp/out', detectors='CS3', verbose=False)
    assert requested == ['events/bc240530cs3_tte_uf.evt.gz']


def test_get_orbit_and_hk_paths(monkeypatch):
    finder = _finder_at('240530')
    requested = []
    monkeypatch.setattr(finder._protocol, 'download',
                        lambda file, dest, verbose: requested.append(file) or Path(dest) / file)

    finder.get_orbit('/tmp/out', verbose=False)
    finder.get_hk('/tmp/out', verbose=False)
    assert requested == ['auxil/bc240530.hk.gz', 'auxil/bc240530csa.hk.gz']


def test_missing_day_is_skipped_cleanly_not_raised(monkeypatch):
    """A missing observation day (or a missing product for an existing day)
    must surface as a per-file None with a warning, never an unhandled
    exception -- day directories are not contiguous.
    """
    finder = _finder_at('991301')  # a day that does not exist

    def _raise_not_found(file, dest, verbose):
        raise OSError('404: Not Found')
    monkeypatch.setattr(finder._protocol, 'download', _raise_not_found)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter('always')
        result = finder.get_cbd('/tmp/out', detectors=['CS0'], verbose=False)

    assert result == [None]
    assert any('Could not download' in str(w.message) for w in record)


def test_get_all_returns_partial_results_for_a_day_with_no_tte(monkeypatch):
    """get_all() must return every product's outcome even when some (e.g.
    TTE, present on only 7 of 89 days) are missing.
    """
    finder = _finder_at('240531')  # a day with CBD but (in this scenario) no TTE

    def _download(file, dest, verbose):
        if 'events/' in file:
            raise OSError('404: Not Found')
        return Path(dest) / file
    monkeypatch.setattr(finder._protocol, 'download', _download)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        result = finder.get_all('/tmp/out', detectors=['CS0'], verbose=False)

    assert result['tte'] == [None]
    assert result['cbd_cl'] != [None]
    assert result['orbit'] is not None


def test_trend_finder_is_ready_to_download_without_an_explicit_cd(monkeypatch):
    """Unlike BurstCubeObsFinder, BurstCubeTrendFinder has no per-instance
    navigation argument, so a bare `BurstCubeTrendFinder()` must be able to
    download immediately -- a user should not have to call a no-op `cd()`
    first just to satisfy BaseFinder's internal state, and this must not
    require a live directory listing to do so (this test has no network
    access at all).
    """
    finder = BurstCubeTrendFinder()
    requested = []
    monkeypatch.setattr(finder._protocol, 'download',
                        lambda file, dest, verbose: requested.append(file) or Path(dest) / file)

    finder.get_attitude('/tmp/out', verbose=False)
    assert requested == ['attitude/bc_csa_att.fits']
