"""Offline tests for gdt.missions.burstcube.time. No network access needed."""
import warnings

import numpy as np
import pytest
from astropy.io import fits

from gdt.missions.burstcube.time import (Time, check_met_epoch,
                                         UnrecognizedMETEpochWarning)


class TestBurstCubeSecTime:
    """Tests for the BurstCube MET time format."""

    def test_met_to_utc_reference_value(self):
        """MET 107629263.33 must convert to 2024-05-30T17:00:26.330 UTC:
        37.000 s *before* 2024-05-30T17:01:03.330, the DATE-OBS the CBD
        extension of bc240530cs0_3cbd_cl.fits.gz carries on disk. That 37 s
        gap is expected, not a bug: DATE-OBS was itself generated under the
        archive's defective epoch (2021-01-01 00:00:00 UTC, per the file's
        own MJDREFI/MJDREFF/TIMESYS read at face value), 37 s later than the
        corrected epoch (2021-01-01 00:00:00 TAI) this package uses. See the
        gdt.missions.burstcube.time module docstring and the README Caveats
        section for the evidence.
        """
        t = Time(107629263.33, format='burstcube')
        assert t.utc.isot == '2024-05-30T17:00:26.330'

    def test_epoch_is_2021_01_01_tai(self):
        """MET 0 must be 2021-01-01T00:00:00.000 TAI (= 2020-12-31T23:59:23.000
        UTC), the corrected BurstCube MET epoch -- not the 2021-01-01
        00:00:00 UTC the archive's own FITS headers state at face value
        (see the module docstring).
        """
        t = Time(0.0, format='burstcube')
        assert t.tai.isot == '2021-01-01T00:00:00.000'
        assert t.utc.isot == '2020-12-31T23:59:23.000'

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
        assert t.utc.isot[1] == '2024-05-30T17:00:26.330'


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


def test_archive_mjdreff_is_37s_later_than_the_corrected_epoch():
    """Regression test for the epoch correction itself.

    Every archive FITS file carries MJDREFI=59215, MJDREFF=0.00080074074074074,
    TIMESYS='TT'. Read at face value (OGIP: MJDREF is expressed in the scale
    TIMESYS names), that resolves to exactly 37.000 s later than this
    package's corrected epoch, 2021-01-01 00:00:00 TAI. If HEASARC ever
    fixes the archive headers to state the corrected epoch directly, this
    test fails loudly -- the same treatment this repo already gives the SAA
    polygon and eb1024 defects.
    """
    archive_epoch = Time(59215, 0.00080074074074074, format='mjd', scale='tt')
    corrected_epoch = Time('2021-01-01 00:00:00', scale='tai')
    assert (archive_epoch - corrected_epoch).sec == pytest.approx(37.0, abs=1e-6)


class TestCheckMetEpoch:
    """Tests for check_met_epoch(), the data-driven epoch check on read."""

    @staticmethod
    def _header(mjdreff, timesys='TT', mjdrefi=59215):
        return fits.Header({'MJDREFI': mjdrefi, 'MJDREFF': mjdreff,
                            'TIMESYS': timesys})

    def test_silent_on_the_corrected_epoch(self):
        """A header already stating the corrected epoch (e.g. a future
        HEASARC fix) must not warn."""
        header = self._header(32.184 / 86400)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            check_met_epoch(header)  # must not raise/warn

    def test_warns_on_the_known_archive_defect(self):
        """A header stating the known defective epoch (2021-01-01 00:00:00
        UTC) must raise exactly one UserWarning, mentioning the 37 s gap."""
        header = self._header(0.00080074074074074)
        with pytest.warns(UserWarning, match='37 s') as record:
            check_met_epoch(header)
        assert len(record) == 1
        assert not any(isinstance(w.message, UnrecognizedMETEpochWarning)
                       for w in record)

    def test_warns_louder_on_an_unrecognized_epoch(self):
        """A header resolving to neither known epoch gets the distinct,
        louder warning class."""
        header = self._header(0.5)
        with pytest.warns(UnrecognizedMETEpochWarning):
            check_met_epoch(header)

    def test_missing_keywords_do_not_raise(self):
        """A header without MJDREFI/MJDREFF/TIMESYS at all (e.g. PRIMARY on
        some products) is simply skipped, not an error."""
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            check_met_epoch(fits.Header())  # must not raise
