"""Offline tests for gdt.missions.burstcube.timeline."""
import numpy as np
import pytest

from gdt.missions.burstcube.timeline import BurstCubeTimeline

_SAMPLE_ROWS = [
    # verified real row from the plugin spec section 18: MET, the CSV's own
    # UTC (37 s earlier than the truth), and the event description.
    (106367555.28030825, '2024-05-16T02:31:58.280', 'Spacecraft Reboot'),
    (106400100.0, '2024-05-16T11:33:23.000',
     'Continous Binned Data (CBD) Capture Turned On'),
    (106450100.0, '2024-05-17T01:46:03.000',
     'CBD Data Timestamp Drift (GPS Reboot). Unnaccounted absolute time drift: -5.14'),
    (106460100.0, '2024-05-17T04:32:43.000',
     'CBD Data Timestamp Drift (GPS Reboot). Unnaccounted absolute time drift: 1.28'),
]


def _write_csv(path):
    with open(path, 'w') as fh:
        for met, utc, desc in _SAMPLE_ROWS:
            fh.write(f'{met},{utc},{desc}\n')


def test_parses_three_columns_with_no_header(tmp_path):
    path = tmp_path / 'timeline.csv'
    _write_csv(path)
    tl = BurstCubeTimeline.open(path)
    assert tl.num_rows == 4
    np.testing.assert_allclose(tl._met, [row[0] for row in _SAMPLE_ROWS])


def test_time_property_uses_met_and_disagrees_with_csv_utc_by_37s(tmp_path):
    """The FITS-verified 37 s discrepancy (spec section 18): the CSV's own
    UTC string is exactly 37.000 s earlier than the MET column converted
    with BurstCubeSecTime.
    """
    path = tmp_path / 'timeline.csv'
    _write_csv(path)
    tl = BurstCubeTimeline.open(path)

    from astropy.time import Time as AstropyTime
    written_utc = AstropyTime(list(tl.utc_as_written), format='isot', scale='utc')
    discrepancy = (tl.time.utc - written_utc).sec
    np.testing.assert_allclose(discrepancy[0], 37.0, atol=1e-3)


def test_utc_as_written_is_kept_separate_and_unconverted(tmp_path):
    path = tmp_path / 'timeline.csv'
    _write_csv(path)
    tl = BurstCubeTimeline.open(path)
    assert tl.utc_as_written[0] == '2024-05-16T02:31:58.280'


def test_drift_seconds_parsed_from_gps_reboot_rows(tmp_path):
    path = tmp_path / 'timeline.csv'
    _write_csv(path)
    tl = BurstCubeTimeline.open(path)

    assert np.isnan(tl.drift_seconds[0])
    assert np.isnan(tl.drift_seconds[1])
    assert tl.drift_seconds[2] == pytest.approx(-5.14)
    assert tl.drift_seconds[3] == pytest.approx(1.28)


def test_matching_tolerates_source_misspellings_without_normalizing(tmp_path):
    path = tmp_path / 'timeline.csv'
    _write_csv(path)
    tl = BurstCubeTimeline.open(path)

    # the source's own misspelling must be matched verbatim...
    assert tl.matching('Continous').sum() == 1
    # ...and the "corrected" spelling must NOT match anything, proving the
    # description text is not being normalized under the hood
    assert tl.matching('Continuous').sum() == 0

    assert tl.matching('Unnaccounted').sum() == 2
    assert tl.matching('Unaccounted').sum() == 0

    assert tl.matching('Spacecraft Reboot').sum() == 1
