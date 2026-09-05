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
"""Offline tests for gdt.missions.burstcube.caldb. No network access needed:
downloads are monkeypatched to fail, forcing the resolver through to the
bundled package data.
"""
import shutil

import numpy as np
import pytest

from gdt.missions.burstcube import caldb

# The CS0 eb16 edges from the CALDB, given in the plugin's spec for testing.
_CS0_EB16_EDGES = [0, 28.48, 62.23, 106.86, 165.56, 233.85, 316.55, 410.77,
                   494.52, 585.14, 685.96, 792.33, 907.6, 1030.33, 1160.67,
                   1298.79, 1648.51]


def _block_downloads(monkeypatch):
    """Make any download attempt raise, simulating no network access."""
    def _raise(*args, **kwargs):
        raise OSError('downloads are blocked in this test')
    monkeypatch.setattr(caldb, '_download', _raise)


@pytest.fixture(autouse=True)
def _no_caldb_env(monkeypatch):
    """Ensure $CALDB is unset by default; individual tests can set it."""
    monkeypatch.delenv('CALDB', raising=False)


class TestResolutionChain:
    """Tests for the $CALDB -> cache -> download -> bundled chain."""

    def test_bundled_fallback_when_caldb_unset_and_downloads_blocked(
            self, tmp_path, monkeypatch, caplog):
        """With $CALDB unset, an empty cache, and downloads blocked, the
        resolver must fall through to the bundled package data and log that
        this is what happened.
        """
        _block_downloads(monkeypatch)
        with caplog.at_level('INFO', logger='gdt.missions.burstcube.caldb'):
            path = caldb.resolve_caldb_file(
                'bcf/align/bccs0_align_20210101v000.fits',
                cache_dir=tmp_path / 'cache')
        assert path.is_file()
        assert 'gdt/missions/burstcube/data' in path.as_posix()
        assert any('bundled package data' in msg for msg in caplog.messages)

    def test_cache_is_used_when_present(self, tmp_path, monkeypatch, caplog):
        """If the file is already in the local cache, the resolver must use
        it directly rather than attempting a download.
        """
        relative = 'bcf/align/bccs0_align_20210101v000.fits'
        cache_dir = tmp_path / 'cache'
        cached_file = cache_dir / relative
        cached_file.parent.mkdir(parents=True)
        bundled = caldb.resolve_caldb_file(relative, cache_dir=tmp_path / 'unused')
        shutil.copyfile(bundled, cached_file)

        def _fail_if_called(*args, **kwargs):
            raise AssertionError('download should not have been attempted')
        monkeypatch.setattr(caldb, '_download', _fail_if_called)

        with caplog.at_level('INFO', logger='gdt.missions.burstcube.caldb'):
            path = caldb.resolve_caldb_file(relative, cache_dir=cache_dir)
        assert path == cached_file
        assert any('local cache' in msg for msg in caplog.messages)

    def test_caldb_env_var_takes_priority(self, tmp_path, monkeypatch, caplog):
        """With $CALDB set and containing the file (in the standard HEASoft
        layout), it must be used ahead of the cache, download, or bundled
        steps.
        """
        relative = 'bcf/align/bccs0_align_20210101v000.fits'
        fake_caldb_root = tmp_path / 'fake_caldb'
        caldb_file = fake_caldb_root / 'data' / 'burstcube' / 'csa' / relative
        caldb_file.parent.mkdir(parents=True)
        bundled = caldb.resolve_caldb_file(relative, cache_dir=tmp_path / 'unused')
        shutil.copyfile(bundled, caldb_file)

        monkeypatch.setenv('CALDB', str(fake_caldb_root))

        def _fail_if_called(*args, **kwargs):
            raise AssertionError('download should not have been attempted')
        monkeypatch.setattr(caldb, '_download', _fail_if_called)

        with caplog.at_level('INFO', logger='gdt.missions.burstcube.caldb'):
            path = caldb.resolve_caldb_file(relative, cache_dir=tmp_path / 'cache')
        assert path == caldb_file
        assert any('$CALDB' in msg for msg in caplog.messages)

    def test_download_populates_cache_and_is_reported(
            self, tmp_path, monkeypatch, caplog):
        """A successful download must be reported and land in the cache at
        the expected relative path.
        """
        relative = 'bcf/align/bccs0_align_20210101v000.fits'
        bundled = caldb.resolve_caldb_file(relative, cache_dir=tmp_path / 'unused')

        def _fake_download(url, dest):
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(bundled, dest)
        monkeypatch.setattr(caldb, '_download', _fake_download)

        cache_dir = tmp_path / 'cache'
        with caplog.at_level('INFO', logger='gdt.missions.burstcube.caldb'):
            path = caldb.resolve_caldb_file(relative, cache_dir=cache_dir)
        assert path == cache_dir / relative
        assert path.is_file()
        assert any('download' in msg for msg in caplog.messages)


class TestEboundsAndRebinSelection:
    """Tests that ebounds/rebin files are selected by INSTRUME, not
    extension index, and return the right per-detector values.
    """

    def test_cs0_eb16_edges_match_expected_values(self, tmp_path, monkeypatch):
        _block_downloads(monkeypatch)
        eb = caldb.ebounds('CS0', 16, cache_dir=tmp_path)
        low = eb.low_edges()
        high = eb.high_edges()
        np.testing.assert_allclose(low, _CS0_EB16_EDGES[:-1], atol=1e-2)
        np.testing.assert_allclose(high, _CS0_EB16_EDGES[1:], atol=1e-2)

    def test_ebounds_accepts_detector_enum_member(self, tmp_path, monkeypatch):
        """ebounds() must accept a BurstCubeDetectors member, not just a
        plain string, and get the same answer either way.
        """
        from gdt.missions.burstcube.detectors import BurstCubeDetectors

        _block_downloads(monkeypatch)
        by_string = caldb.ebounds('CS0', 16, cache_dir=tmp_path)
        by_enum = caldb.ebounds(BurstCubeDetectors.CS0, 16, cache_dir=tmp_path)
        np.testing.assert_allclose(by_string.low_edges(), by_enum.low_edges())

    def test_different_detectors_give_different_edges(self, tmp_path, monkeypatch):
        """CS0 and CS3 have distinct ebounds (per the plugin spec, their
        first channel differs: CS0 28.48-62.23 keV, CS3 29.22-62.30 keV), so
        selecting by INSTRUME must actually distinguish them.
        """
        _block_downloads(monkeypatch)
        cs0 = caldb.ebounds('CS0', 16, cache_dir=tmp_path)
        cs3 = caldb.ebounds('CS3', 16, cache_dir=tmp_path)
        assert cs0.low_edges()[1] != cs3.low_edges()[1]

    def test_ebounds_selection_is_not_by_extension_index(self, tmp_path, monkeypatch):
        """Directly confirm that all 4 ENERGY_EBOUNDS extensions in the
        bundled eb16 file are distinguishable only by INSTRUME (i.e. an
        index-based reader would silently get this wrong for at least one
        detector).
        """
        from astropy.io import fits

        _block_downloads(monkeypatch)
        path = caldb.resolve_caldb_file(
            'cpf/ebounds/bccsa_eb16_20221001v001.fits', cache_dir=tmp_path)
        with fits.open(path) as hdulist:
            instrumes = [hdu.header['INSTRUME'] for hdu in hdulist[1:]]
        assert instrumes == ['CS0', 'CS1', 'CS2', 'CS3']

    def test_reb16_first_and_last_rows(self, tmp_path, monkeypatch):
        """CS0 reb16 first rows are (0,23,1), (24,45,1), (46,74,1); last row
        is (817,1023,1), per the plugin spec.
        """
        _block_downloads(monkeypatch)
        reb = caldb.rebin('CS0', 16, cache_dir=tmp_path)
        np.testing.assert_array_equal(reb.chan_min[:3], [0, 24, 46])
        np.testing.assert_array_equal(reb.chan_max[:3], [23, 45, 74])
        assert reb.chan_min[-1] == 817
        assert reb.chan_max[-1] == 1023

    def test_reb64_first_rows(self, tmp_path, monkeypatch):
        """CS0 reb64 first rows are (0,23,1), (24,27,1), (28,31,1), per the
        plugin spec.
        """
        _block_downloads(monkeypatch)
        reb = caldb.rebin('CS0', 64, cache_dir=tmp_path)
        np.testing.assert_array_equal(reb.chan_min[:3], [0, 24, 28])
        np.testing.assert_array_equal(reb.chan_max[:3], [23, 27, 31])


class TestAlignmentAndSaa:
    """Tests for the alignment and SAA-region CALDB accessors."""

    def test_alignment_matches_bundled_header(self, tmp_path, monkeypatch):
        _block_downloads(monkeypatch)
        al = caldb.alignment('CS0', cache_dir=tmp_path)
        np.testing.assert_allclose(
            al.matrix[0], [0.7071067811865475, -0.5, -0.5])
        assert al.rollsign == 1
        assert al.rolloff == 0.0

    def test_saa_region_is_a_19_point_polygon(self, tmp_path, monkeypatch):
        _block_downloads(monkeypatch)
        saa = caldb.saa_region(cache_dir=tmp_path)
        assert saa.shape == 'POLYGON'
        assert len(saa.longitude) == 19
        assert len(saa.latitude) == 19
        assert saa.longitude[0] == pytest.approx(-30.0, abs=1e-3)
        assert saa.latitude[0] == pytest.approx(33.9, abs=1e-3)
