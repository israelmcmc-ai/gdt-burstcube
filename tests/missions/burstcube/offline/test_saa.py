"""Offline tests for gdt.missions.burstcube.saa. Uses the bundled CALDB SAA
region file (downloads blocked, forcing the bundled fallback), same as the
round-1 CALDB tests.
"""
import numpy as np
import pytest

from gdt.missions.burstcube import caldb
from gdt.missions.burstcube.saa import BurstCubeSAA, _crossings


def test_polygon_matches_caldb_file_and_is_closed(tmp_path, monkeypatch):
    """19 vertices in the file, plus a 20th repeating the first so the
    polygon is closed -- the file leaves it open, which draws a broken
    outline and makes gdt-core's ``is_closed`` false."""
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))
    saa = BurstCubeSAA(cache_dir=tmp_path)

    assert caldb.saa_region(cache_dir=tmp_path).latitude.size == 19
    assert saa.num_points == 20
    assert saa.is_closed()
    assert saa.latitude[0] == pytest.approx(-30.0, abs=1e-3)
    assert saa.longitude[0] == pytest.approx(33.9, abs=1e-3)


#: Fermi GBM's ``GbmSaaPolygon5`` (gdt.missions.fermi.gbm.saa), where these
#: same numbers are stored in explicitly named _latitude/_longitude lists.
#: BurstCube's CALDB polygon opens with exactly these 11 vertices and then
#: extends further west, which is the independent check that its X column is
#: latitude and its Y column longitude, not the other way round.
GBM_POLYGON_5_LATITUDE = [-30.000, -19.867, -9.733, 0.400, 2.000, 2.000,
                          -1.000, -6.155, -8.880, -14.220, -18.404]
GBM_POLYGON_5_LONGITUDE = [33.900, 12.398, -9.103, -30.605, -38.400, -45.000,
                           -65.000, -84.000, -89.200, -94.300, -94.300]


def test_polygon_opens_with_gbms_own_saa_vertices(tmp_path, monkeypatch):
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))
    saa = BurstCubeSAA(cache_dir=tmp_path)

    np.testing.assert_allclose(saa.latitude[:11], GBM_POLYGON_5_LATITUDE,
                               atol=1e-3)
    np.testing.assert_allclose(saa.longitude[:11], GBM_POLYGON_5_LONGITUDE,
                               atol=1e-3)


def test_polygon_covers_the_south_atlantic_anomaly(tmp_path, monkeypatch):
    """Read with the columns swapped the polygon is a narrow band off Brazil
    that the orbit crosses at the wrong times; read correctly it is the SAA.
    """
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))
    saa = BurstCubeSAA(cache_dir=tmp_path)

    assert saa.latitude.min() == pytest.approx(-53.6, abs=0.1)
    assert saa.latitude.max() == pytest.approx(2.0, abs=0.1)
    assert saa.longitude.min() == pytest.approx(-94.3, abs=0.1)
    assert saa.longitude.max() == pytest.approx(33.9, abs=0.1)

    # a point in the heart of the anomaly, and two well outside it
    assert saa.contains(-45.0, -30.0) is True      # 45W, 30S
    assert saa.contains(0.0, 0.0) is False         # Gulf of Guinea
    assert saa.contains(100.0, 40.0) is False      # central Asia


def test_contains_interior_and_exterior_points(tmp_path, monkeypatch):
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))
    saa = BurstCubeSAA(cache_dir=tmp_path)

    assert saa.contains(-20.0, -40.0) is True
    assert saa.contains(100.0, 50.0) is False


def test_contains_accepts_arrays(tmp_path, monkeypatch):
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))
    saa = BurstCubeSAA(cache_dir=tmp_path)

    lon = np.array([-20.0, 100.0])
    lat = np.array([-40.0, 50.0])
    result = saa.contains(lon, lat)
    np.testing.assert_array_equal(result, [True, False])


def test_polygon_is_simple_after_reordering(tmp_path, monkeypatch):
    """The file lists two pairs of vertices out of order, so its boundary
    doubles back and crosses itself twice -- visible as a spur off each of
    the eastern and western corners. Sorting by angle about the centroid
    removes both.
    """
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))
    region = caldb.saa_region(cache_dir=tmp_path)
    saa = BurstCubeSAA(cache_dir=tmp_path)

    assert _crossings(region.latitude, region.longitude) == 2
    assert _crossings(saa.latitude, saa.longitude) == 0


def test_reordering_moves_only_the_two_transposed_pairs(tmp_path, monkeypatch):
    """Everything except vertices 11/12 and 17/18 keeps the file's order, so
    this is a repair of two transpositions, not a wholesale rearrangement.
    """
    monkeypatch.setattr(caldb, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(OSError('blocked')))
    region = caldb.saa_region(cache_dir=tmp_path)
    saa = BurstCubeSAA(cache_dir=tmp_path)

    file_order = list(zip(np.round(region.latitude, 3),
                          np.round(region.longitude, 3)))
    fixed = list(zip(np.round(saa.latitude[:-1], 3),
                     np.round(saa.longitude[:-1], 3)))

    expected = list(file_order)
    expected[11], expected[12] = expected[12], expected[11]
    expected[17], expected[18] = expected[18], expected[17]
    assert fixed == expected


def test_already_simple_polygon_is_left_alone():
    """A future CALDB revision that is already correctly ordered must pass
    through untouched, including a concave one no angular sweep would order
    the same way."""
    latitude = np.array([0.0, 0.0, 5.0, 10.0, 10.0])
    longitude = np.array([0.0, 10.0, 5.0, 10.0, 0.0])   # concave, simple
    assert _crossings(latitude, longitude) == 0

    out_lat, out_lon = BurstCubeSAA._ordered_by_angle(latitude, longitude)
    np.testing.assert_array_equal(out_lat, latitude)
    np.testing.assert_array_equal(out_lon, longitude)
