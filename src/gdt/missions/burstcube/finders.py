"""Finders for the BurstCube HEASARC archive
(``https://heasarc.gsfc.nasa.gov/FTP/burstcube/data/``).

Unlike GBM, where a trigger's exact filenames are not fully predictable in
advance, every BurstCube archive filename is fully determined by the
observation day, detector, and product type (see the archive layout in the
plugin spec). So, unlike :class:`gdt.missions.fermi.gbm.finders.GbmFinder`,
these finders do not need to list a remote directory to know what to
download: they build the exact relative file path for each requested product
and download it directly, which also means a missing file (including an
entire missing observation day -- day directories are **not** contiguous)
is simply a per-file 404, handled by skipping that file with a warning
rather than raising.
"""
import warnings
from pathlib import Path
from typing import List, Optional, Union

from gdt.core.heasarc import BaseFinder

from .time import Time as AstropyTime

__all__ = ['BurstCubeObsFinder', 'BurstCubeTrendFinder']

_ALL_DETECTORS = (0, 1, 2, 3)


def _detector_numbers(detectors) -> List[int]:
    """Normalize a detector or list of detectors to a list of numbers 0-3.

    Args:
        detectors: None (all 4), or one or more of an int, a detector name
            string (e.g. ``'CS0'``), or a
            :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`
            member.

    Returns:
        (list of int)
    """
    if detectors is None:
        return list(_ALL_DETECTORS)
    if not isinstance(detectors, (list, tuple)):
        detectors = [detectors]
    numbers = []
    for det in detectors:
        if hasattr(det, 'number'):
            numbers.append(det.number)
        elif isinstance(det, str):
            numbers.append(int(det.upper().removeprefix('CS')))
        else:
            numbers.append(int(det))
    return numbers


class _BurstCubeFinderMixin(BaseFinder):
    """Shared download-resilience behavior for the BurstCube finders: since
    filenames are deterministic, downloads are attempted directly rather
    than via a directory listing, and a missing file is skipped (with a
    warning) instead of raising.
    """

    def _download_each(self, download_dir: Union[str, Path],
                       files: List[str], verbose: bool = True) -> List[Optional[Path]]:
        """Download each file independently, skipping (with a warning)
        any that fail -- e.g. a product that does not exist for this
        observation day.

        Args:
            download_dir (str or Path): The download directory
            files (list of str): The file paths, relative to this finder's
                current directory
            verbose (bool, optional): Passed through to the download call.
                Default is True.

        Returns:
            (list of Path or None): One entry per requested file; None
            where that file could not be downloaded.
        """
        paths = []
        for file in files:
            try:
                paths.append(self._protocol.download(file, download_dir, verbose))
            except Exception as err:
                warnings.warn(f'Could not download {file}: {err}',
                              RuntimeWarning, stacklevel=2)
                paths.append(None)
        return paths


class BurstCubeObsFinder(_BurstCubeFinderMixin):
    """Finder for a single BurstCube observation day
    (``obs/YYYY_MM/YYMMDD/``).

    Parameters:
        obs_id (str or astropy.time.Time, optional): The observation day,
            e.g. ``'240530'``, or a :class:`~astropy.time.Time` (any format
            that converts to ``burstcube_obsid``). If omitted, call
            :meth:`cd` before downloading.
        protocol (str, optional): The connection protocol. Default is HTTPS.
    """
    _root = 'burstcube/data'

    def cd(self, obs_id: Union[str, AstropyTime]):
        """Change to a new observation day.

        Args:
            obs_id (str or astropy.time.Time): The observation day
        """
        super().cd(obs_id)

    def get_cbd(self, download_dir, detectors=None, variant='cl', **kwargs):
        """Download CBD files for this observation day.

        Args:
            download_dir (str): The download directory
            detectors (optional): One or more detectors to download (see
                :func:`_detector_numbers`). If omitted, downloads all 4.
            variant (str, optional): ``'cl'`` (cleaned) or ``'uf'``
                (unfiltered). Default is ``'cl'``.
            verbose (bool, optional): If True, outputs the download status.

        Returns:
            (list of Path or None)
        """
        if variant not in ('cl', 'uf'):
            raise ValueError("variant must be 'cl' or 'uf'")
        obs_id = self.obs_id
        files = [f'monitor/bc{obs_id}cs{n}_3cbd_{variant}.fits.gz'
                for n in _detector_numbers(detectors)]
        return self._download_each(download_dir, files, **kwargs)

    def get_tte(self, download_dir, detectors=None, **kwargs):
        """Download TTE files for this observation day. Per the archive,
        TTE exists on only 7 days across the whole mission; on any other
        day, every file downloaded here will be reported as missing.

        Args:
            download_dir (str): The download directory
            detectors (optional): One or more detectors to download. If
                omitted, downloads all 4.
            verbose (bool, optional): If True, outputs the download status.

        Returns:
            (list of Path or None)
        """
        obs_id = self.obs_id
        files = [f'events/bc{obs_id}cs{n}_tte_uf.evt.gz'
                for n in _detector_numbers(detectors)]
        return self._download_each(download_dir, files, **kwargs)

    def get_orbit(self, download_dir, **kwargs):
        """Download the orbit/ephemeris file for this observation day.

        Args:
            download_dir (str): The download directory
            verbose (bool, optional): If True, outputs the download status.

        Returns:
            (Path or None)
        """
        obs_id = self.obs_id
        return self._download_each(download_dir, [f'auxil/bc{obs_id}.hk.gz'],
                                   **kwargs)[0]

    def get_hk(self, download_dir, **kwargs):
        """Download the detector housekeeping file for this observation day.

        Args:
            download_dir (str): The download directory
            verbose (bool, optional): If True, outputs the download status.

        Returns:
            (Path or None)
        """
        obs_id = self.obs_id
        return self._download_each(
            download_dir, [f'auxil/bc{obs_id}csa.hk.gz'], **kwargs)[0]

    def get_all(self, download_dir, detectors=None, **kwargs):
        """Download every product available for this observation day: both
        CBD variants, TTE, the orbit file, and the detector housekeeping
        file. Products that do not exist for this day (e.g. TTE on the 82
        of 89 days without it) are skipped with a warning, not raised.

        Args:
            download_dir (str): The download directory
            detectors (optional): One or more detectors to download. If
                omitted, downloads all 4.
            verbose (bool, optional): If True, outputs the download status.

        Returns:
            (dict): Keys ``'cbd_cl'``, ``'cbd_uf'``, ``'tte'`` (each a list
            of Path or None), and ``'orbit'``, ``'hk'`` (each a Path or None).
        """
        return {
            'cbd_cl': self.get_cbd(download_dir, detectors, 'cl', **kwargs),
            'cbd_uf': self.get_cbd(download_dir, detectors, 'uf', **kwargs),
            'tte': self.get_tte(download_dir, detectors, **kwargs),
            'orbit': self.get_orbit(download_dir, **kwargs),
            'hk': self.get_hk(download_dir, **kwargs),
        }

    def _construct_path(self, obs_id: Union[str, AstropyTime]) -> str:
        """Construct the remote path for an observation day.

        Args:
            obs_id (str or astropy.time.Time): The observation day, e.g.
                ``'240530'``

        Returns:
            (str): e.g. ``'burstcube/data/obs/2024_05/240530'``
        """
        obs_id = self.obs_id_from(obs_id)
        year, month = '20' + obs_id[0:2], obs_id[2:4]
        return f'{self._root}/obs/{year}_{month}/{obs_id}'

    @property
    def obs_id(self) -> Optional[str]:
        """(str or None): The observation day currently selected, as
        ``YYMMDD``, or None before the first :meth:`cd`."""
        args = getattr(self, '_args', None)
        return self.obs_id_from(args[0]) if args else None

    @staticmethod
    def obs_id_from(when: Union[str, AstropyTime]) -> str:
        """The ``YYMMDD`` observation day containing a given time.

        BurstCube organizes the archive by UTC calendar day, and the
        ``burstcube_obsid`` time format registered by
        :mod:`gdt.missions.burstcube.time` is the conversion, so this is
        also available directly as ``Time(...).burstcube_obsid``. Every
        method here that takes an observation day -- the constructor,
        :meth:`cd` -- accepts a :class:`~astropy.time.Time` too and runs it
        through this, so an explicit conversion is only needed when you want
        the string itself (to name a directory, say).

        Args:
            when (str or astropy.time.Time): A time, in any format
                :class:`~astropy.time.Time` accepts -- including BurstCube
                MET via ``Time(met, format='burstcube')`` -- or a
                ``YYMMDD`` string, which is returned unchanged.

        Returns:
            (str): The observation day, e.g. ``'240530'``

        Example:
            >>> from gdt.missions.burstcube.time import Time
            >>> BurstCubeObsFinder.obs_id_from(Time(114214208.4, format='burstcube'))
            '240814'
        """
        if isinstance(when, AstropyTime):
            return when.burstcube_obsid
        return str(when)


class BurstCubeTrendFinder(_BurstCubeFinderMixin):
    """Finder for the BurstCube mission-long trend products
    (``trend/attitude``, ``trend/gti_binning``, ``trend/gti_poscnt``,
    ``trend/gti_saa``, ``trend/timeline``), which are not organized by
    observation day.

    Unlike :class:`BurstCubeObsFinder`, there is no per-instance navigation
    parameter to pass (the trend directory is a single fixed location), so
    this finder is ready to download from immediately on construction --
    no separate, argument-less ``cd()`` call is needed first. As the module
    docstring notes, BurstCube filenames are fully deterministic, so (per
    :class:`_BurstCubeFinderMixin`) this never needs a real directory
    listing to know what exists; the fixed trend path is set directly
    rather than through :meth:`~gdt.core.heasarc.BaseFinder.cd`, which
    would otherwise perform one.

    Parameters:
        protocol (str, optional): The connection protocol. Default is HTTPS.
    """
    _root = 'burstcube/data/trend'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self._args is None:
            self._args = ()
            self._cwd = self._construct_path()
            self._protocol._cd(self._cwd)

    def get_attitude(self, download_dir, **kwargs):
        """Download the attitude file (3 rows, whole mission).

        Args:
            download_dir (str): The download directory
            verbose (bool, optional): If True, outputs the download status.

        Returns:
            (Path or None)
        """
        return self._download_each(
            download_dir, ['attitude/bc_csa_att.fits'], **kwargs)[0]

    def get_gti_binning(self, download_dir, detector, **kwargs):
        """Download the "valid 0.256 s binning" trend GTI for one detector.

        Args:
            download_dir (str): The download directory
            detector: The detector (see :func:`_detector_numbers`)
            verbose (bool, optional): If True, outputs the download status.

        Returns:
            (Path or None)
        """
        n = _detector_numbers(detector)[0]
        return self._download_each(
            download_dir, [f'gti_binning/bc_cs{n}_bin0_256.gti'],
            **kwargs)[0]

    def get_gti_poscnt(self, download_dir, detector, **kwargs):
        """Download the "positive counts" trend GTI for one detector.

        Args:
            download_dir (str): The download directory
            detector: The detector (see :func:`_detector_numbers`)
            verbose (bool, optional): If True, outputs the download status.

        Returns:
            (Path or None)
        """
        n = _detector_numbers(detector)[0]
        return self._download_each(
            download_dir, [f'gti_poscnt/bc_cs{n}_poscounts.gti'],
            **kwargs)[0]

    def get_gti_saa(self, download_dir, variant='in_cl', **kwargs):
        """Download an SAA trend GTI.

        Args:
            download_dir (str): The download directory
            variant (str, optional): One of ``'in'``, ``'out'``, ``'in_cl'``,
                ``'out_cl'`` (``_cl`` is the filtered version). Default is
                ``'in_cl'``.
            verbose (bool, optional): If True, outputs the download status.

        Returns:
            (Path or None)
        """
        valid = ('in', 'out', 'in_cl', 'out_cl')
        if variant not in valid:
            raise ValueError(f'variant must be one of {valid}')
        return self._download_each(
            download_dir, [f'gti_saa/bc_csa_saa_{variant}.gti'],
            **kwargs)[0]

    def get_timeline(self, download_dir, filename='20250606_timeline_final.csv',
                     **kwargs):
        """Download the mission timeline CSV.

        Note:
            The timeline filename is date-stamped and may change with a
            future archive release; ``filename`` defaults to the one
            verified against the live archive (see the plugin spec), but
            can be overridden.

        Args:
            download_dir (str): The download directory
            filename (str, optional): The remote filename.
            verbose (bool, optional): If True, outputs the download status.

        Returns:
            (Path or None)
        """
        return self._download_each(
            download_dir, [f'timeline/{filename}'], **kwargs)[0]

    def _construct_path(self) -> str:
        """Construct the remote path for the trend directory.

        Returns:
            (str): ``'burstcube/data/trend'``
        """
        return self._root
