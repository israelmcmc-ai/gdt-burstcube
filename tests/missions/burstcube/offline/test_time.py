"""Offline tests for gdt.missions.burstcube.time. No network access needed."""
import numpy as np
import pytest

from gdt.missions.burstcube.time import Time


class TestBurstCubeSecTime:
    """Tests for the BurstCube MET time format."""

    def test_met_to_utc_reference_value(self):
        """MET 107629263.33 must convert to 2024-05-30T17:01:03.330 UTC. This
        is the DATE-OBS of the CBD extension of bc240530cs0_3cbd_cl.fits.gz,
        so it validates the whole epoch/scale definition against a real
        archive file, not just the arithmetic.
        """
        t = Time(107629263.33, format='burstcube')
        assert t.utc.isot == '2024-05-30T17:01:03.330'

    def test_epoch_is_2021_01_01_utc_in_tt(self):
        """MET 0 must be 2021-01-01T00:00:00.000 UTC, per MJDREFI=59215."""
        t = Time(0.0, format='burstcube')
        assert t.utc.isot == '2021-01-01T00:00:00.000'

    def test_round_trip(self):
        """Converting a MET to UTC and back must reproduce the original MET,
        to within floating-point precision.
        """
        met = 107629263.33
        t = Time(met, format='burstcube')
        recovered = Time(t.utc.isot, format='isot', scale='utc').burstcube
        assert recovered == pytest.approx(met, abs=1e-3)

    def test_array_input(self):
        """The format must work on an array of MET values, not just scalars."""
        mets = np.array([0.0, 107629263.33, 2 * 107629263.33])
        t = Time(mets, format='burstcube')
        assert t.utc.isot[1] == '2024-05-30T17:01:03.330'


class TestBurstCubeObsId:
    """Tests for the BurstCube observation-day (YYMMDD) time format."""

    def test_round_trip(self):
        """A YYMMDD string must round-trip through Time unchanged."""
        obsid = '240530'
        t = Time(obsid, format='burstcube_obsid')
        assert t.burstcube_obsid == obsid

    def test_represents_midnight_utc(self):
        """The represented instant is midnight UTC at the start of the day."""
        t = Time('240530', format='burstcube_obsid')
        assert t.utc.isot == '2024-05-30T00:00:00.000'

    def test_array_round_trip(self):
        """An array of YYMMDD strings must round-trip too."""
        obsids = ['240530', '240602', '240813']
        t = Time(obsids, format='burstcube_obsid')
        assert list(t.burstcube_obsid) == obsids

    def test_rejects_malformed_string(self):
        """A string that isn't YYMMDD must raise, not silently misparse.
        Astropy wraps the format's own TypeError in a ValueError when the
        format is given explicitly (one candidate format, so it reports
        immediately with the original error chained in).
        """
        with pytest.raises(ValueError, match='burstcube_obsid'):
            Time('2024-05-30', format='burstcube_obsid')

    def test_accepts_time_or_string_interchangeably(self):
        """A finder should be able to accept either a plain string or a
        Time object built from that string and get the same day.
        """
        from_string = Time('240530', format='burstcube_obsid')
        from_time = Time(from_string)
        assert from_time.burstcube_obsid == '240530'
