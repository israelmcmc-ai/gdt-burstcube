"""Access to the BurstCube CALDB (calibration database) files.

A request for a CALDB file is resolved through a chain of four steps, in this
order, and the resolver logs which step actually satisfied the request:

1. ``$CALDB``, if set, using the standard HEASoft layout,
   ``$CALDB/data/burstcube/csa/<subdir>/<filename>``.
2. A local cache directory (by default ``~/.gdt/burstcube/caldb``).
3. A fresh download into that cache, from
   ``https://heasarc.gsfc.nasa.gov/FTP/caldb/data/burstcube/csa/``.
4. The bundled package data in ``gdt.missions.burstcube.data``, which is
   always present and lets everything in this plugin work fully offline.

Only the 10 files needed for the offline plugin (the 4 detector alignment
files, the SAA region, all 3 ebounds files, and both rebin files) are bundled;
the CALDB response files are not (see ``BurstCubeRspFinder``, added in a later
version of this plugin, which downloads only the response pixels actually
needed for a calculation).
"""
import logging
import os
import shutil
import tempfile
import urllib.request
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Optional, Union
from urllib.error import URLError

import numpy as np
from astropy.io import fits

from gdt.core.data_primitives import Ebounds

__all__ = ['alignment', 'ebounds', 'rebin', 'saa_region', 'resolve_caldb_file',
          'Alignment', 'Rebinning', 'SaaRegion', 'CALDB_REMOTE_ROOT',
          'DEFAULT_CACHE_DIR']

log = logging.getLogger(__name__)

#: The CALDB remote root for the BurstCube CSA instrument.
CALDB_REMOTE_ROOT = 'https://heasarc.gsfc.nasa.gov/FTP/caldb/data/burstcube/csa/'

#: The standard HEASoft CALDB layout, relative to $CALDB.
_CALDB_LOCAL_SUBPATH = Path('data') / 'burstcube' / 'csa'

#: The default local cache directory, used if $CALDB is unset or does not
#: have the file, and no `cache_dir` is given explicitly.
DEFAULT_CACHE_DIR = Path.home() / '.gdt' / 'burstcube' / 'caldb'

_BUNDLED_PACKAGE = 'gdt.missions.burstcube.data'

_ALIGN_FILES = {f'CS{n}': f'bcf/align/bccs{n}_align_20210101v000.fits'
               for n in range(4)}
_SAA_FILE = 'bcf/saa/bccsa_saareg_20230101v001.fits'
_EBOUNDS_FILES = {16: 'cpf/ebounds/bccsa_eb16_20221001v001.fits',
                  64: 'cpf/ebounds/bccsa_eb64_20221001v001.fits',
                  1024: 'cpf/ebounds/bccsa_eb1024_20221001v001.fits'}
_REBIN_FILES = {16: 'cpf/rebin/bccsa_reb16_20221001v001.fits',
                64: 'cpf/rebin/bccsa_reb64_20221001v001.fits'}


@dataclass(frozen=True)
class Alignment:
    """The detector alignment for one BurstCube detector, from
    ``bcf/align/bccsN_align_20210101v000.fits``.

    Attributes:
        matrix (numpy.ndarray): The (3, 3) DET->SAT rotation matrix
        rollsign (int): The ``ROLLSIGN`` keyword
        rolloff (float): The ``ROLLOFF`` keyword, in degrees
    """
    matrix: np.ndarray
    rollsign: int
    rolloff: float


@dataclass(frozen=True)
class Rebinning:
    """A CALDB channel-rebinning scheme for one BurstCube detector, from
    ``cpf/rebin/bccsa_rebNN_20221001v001.fits``.

    Attributes:
        chan_min (numpy.ndarray): The lowest native (1024-channel) channel in
            each rebinned group
        chan_max (numpy.ndarray): The highest native (1024-channel) channel in
            each rebinned group
        rebinning (numpy.ndarray): The number of native channels combined
            into each rebinned group (``chan_max - chan_min + 1``)
    """
    chan_min: np.ndarray
    chan_max: np.ndarray
    rebinning: np.ndarray


@dataclass(frozen=True)
class SaaRegion:
    """The South Atlantic Anomaly boundary polygon, from
    ``bcf/saa/bccsa_saareg_20230101v001.fits``.

    Attributes:
        shape (str): The region shape, e.g. ``'POLYGON'``
        longitude (numpy.ndarray): The vertex East longitudes, in degrees
        latitude (numpy.ndarray): The vertex latitudes, in degrees
        r (numpy.ndarray): The ``R`` column (unused for a polygon shape)
        rotang (numpy.ndarray): The ``ROTANG`` column (unused for a polygon
            shape)
        component (int): The ``COMPONENT`` number
    """
    shape: str
    longitude: np.ndarray
    latitude: np.ndarray
    r: np.ndarray
    rotang: np.ndarray
    component: int


def resolve_caldb_file(relative_path: Union[str, Path],
                       cache_dir: Optional[Path] = None) -> Path:
    """Resolve a single CALDB file through the resolution chain, logging
    which step satisfied the request.

    Args:
        relative_path (str or Path): The file's path relative to the
            BurstCube CSA CALDB root, e.g.
            ``'bcf/align/bccs0_align_20210101v000.fits'``.
        cache_dir (Path, optional): The local cache directory to check and
            download into. Defaults to :data:`DEFAULT_CACHE_DIR`.

    Returns:
        (Path): The resolved local file path.
    """
    relative_path = Path(relative_path)
    cache_dir = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR

    # 1. $CALDB, using the standard HEASoft layout.
    caldb_env = os.environ.get('CALDB')
    if caldb_env:
        candidate = Path(caldb_env) / _CALDB_LOCAL_SUBPATH / relative_path
        if candidate.is_file():
            log.info('Resolved %s via $CALDB (%s)', relative_path, candidate)
            return candidate

    # 2. Local cache.
    cached = cache_dir / relative_path
    if cached.is_file():
        log.info('Resolved %s via local cache (%s)', relative_path, cached)
        return cached

    # 3. Download into the local cache.
    try:
        _download(CALDB_REMOTE_ROOT + relative_path.as_posix(), cached)
    except (URLError, OSError) as err:
        log.info('Download of %s failed (%s); falling back to bundled data',
                relative_path, err)
    else:
        log.info('Resolved %s via download (%s)', relative_path, cached)
        return cached

    # 4. Bundled package data (always present; the offline fallback).
    bundled = files(_BUNDLED_PACKAGE).joinpath(relative_path.name)
    log.info('Resolved %s via bundled package data (%s)', relative_path,
            bundled)
    return Path(bundled)


def _download(url: str, dest: Path):
    """Download a single CALDB file to `dest`, creating parent directories as
    needed. Downloads to a temporary file first and moves it into place, so
    that a failed or interrupted download never leaves a corrupt file in the
    cache.

    Args:
        url (str): The remote URL
        dest (Path): The local destination path
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=30) as response:
        with tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as tmp:
            shutil.copyfileobj(response, tmp)
            tmp_path = Path(tmp.name)
    tmp_path.replace(dest)


def _select_by_instrume(hdulist: fits.HDUList, instrume: str) -> fits.BinTableHDU:
    """Select the extension of a multi-detector CALDB file whose ``INSTRUME``
    keyword matches. **Never select by extension index**: the ebounds and
    rebin files have 4 extensions, one per detector, and their order is not
    guaranteed by CALDB.

    Args:
        hdulist (astropy.io.fits.HDUList): The open FITS file
        instrume (str): The detector name to match, e.g. ``'CS0'``

    Returns:
        (astropy.io.fits.BinTableHDU)
    """
    for hdu in hdulist[1:]:
        if hdu.header.get('INSTRUME') == instrume:
            return hdu
    raise ValueError(f'No extension with INSTRUME={instrume!r} found in '
                     f'{hdulist.filename()}')


def _det_name(det) -> str:
    """Normalize a detector reference (a string or a
    :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors` member) to
    its short CALDB/INSTRUME name, e.g. ``'CS0'``.

    Args:
        det (str or Enum): The detector

    Returns:
        (str)
    """
    return str(getattr(det, 'name', det)).upper()


def alignment(det, cache_dir: Optional[Path] = None) -> Alignment:
    """Retrieve the CALDB alignment for one detector.

    Args:
        det (str or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`):
            The detector, e.g. ``'CS0'``
        cache_dir (Path, optional): The local CALDB cache directory. Defaults
            to :data:`DEFAULT_CACHE_DIR`.

    Returns:
        (:class:`Alignment`)
    """
    name = _det_name(det)
    path = resolve_caldb_file(_ALIGN_FILES[name], cache_dir=cache_dir)
    with fits.open(path) as hdulist:
        header = hdulist[0].header
        matrix = np.array([
            [header['ALIGNM11'], header['ALIGNM12'], header['ALIGNM13']],
            [header['ALIGNM21'], header['ALIGNM22'], header['ALIGNM23']],
            [header['ALIGNM31'], header['ALIGNM32'], header['ALIGNM33']],
        ])
        rollsign = header['ROLLSIGN']
        rolloff = header['ROLLOFF']
    return Alignment(matrix=matrix, rollsign=rollsign, rolloff=rolloff)


def ebounds(det, nchan: int = 16, cache_dir: Optional[Path] = None) -> Ebounds:
    """Retrieve the CALDB energy-channel bounds for one detector.

    Args:
        det (str or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`):
            The detector, e.g. ``'CS0'``
        nchan (int, optional): The number of channels: 16, 64, or 1024.
            Default is 16.
        cache_dir (Path, optional): The local CALDB cache directory. Defaults
            to :data:`DEFAULT_CACHE_DIR`.

    Returns:
        (:class:`~gdt.core.data_primitives.Ebounds`)
    """
    name = _det_name(det)
    path = resolve_caldb_file(_EBOUNDS_FILES[nchan], cache_dir=cache_dir)
    with fits.open(path) as hdulist:
        hdu = _select_by_instrume(hdulist, name)
        e_min = np.asarray(hdu.data['E_MIN'], dtype=float)
        e_max = np.asarray(hdu.data['E_MAX'], dtype=float)
    return Ebounds.from_bounds(e_min, e_max)


def rebin(det, nchan: int = 16, cache_dir: Optional[Path] = None) -> Rebinning:
    """Retrieve the CALDB channel-rebinning scheme (1024 native channels down
    to `nchan`) for one detector.

    Args:
        det (str or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`):
            The detector, e.g. ``'CS0'``
        nchan (int, optional): The rebinned number of channels: 16 or 64.
            Default is 16.
        cache_dir (Path, optional): The local CALDB cache directory. Defaults
            to :data:`DEFAULT_CACHE_DIR`.

    Returns:
        (:class:`Rebinning`)
    """
    name = _det_name(det)
    path = resolve_caldb_file(_REBIN_FILES[nchan], cache_dir=cache_dir)
    with fits.open(path) as hdulist:
        hdu = _select_by_instrume(hdulist, name)
        chan_min = np.asarray(hdu.data['CHANMIN'], dtype=int)
        chan_max = np.asarray(hdu.data['CHANMAX'], dtype=int)
        rebinning = np.asarray(hdu.data['REBINNING'], dtype=int)
    return Rebinning(chan_min=chan_min, chan_max=chan_max, rebinning=rebinning)


def saa_region(cache_dir: Optional[Path] = None) -> SaaRegion:
    """Retrieve the CALDB South Atlantic Anomaly boundary polygon.

    Args:
        cache_dir (Path, optional): The local CALDB cache directory. Defaults
            to :data:`DEFAULT_CACHE_DIR`.

    Returns:
        (:class:`SaaRegion`)
    """
    path = resolve_caldb_file(_SAA_FILE, cache_dir=cache_dir)
    with fits.open(path) as hdulist:
        row = hdulist['REGION'].data[0]
    return SaaRegion(shape=str(row['SHAPE']).strip(),
                     longitude=np.asarray(row['X'], dtype=float),
                     latitude=np.asarray(row['Y'], dtype=float),
                     r=np.asarray(row['R'], dtype=float),
                     rotang=np.asarray(row['ROTANG'], dtype=float),
                     component=int(row['COMPONENT']))
