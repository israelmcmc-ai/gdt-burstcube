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
"""BurstCube detector response: one DRM per HEALPix pixel per detector, from
CALDB ``cpf/response/csN/{th0_90,th91_180}/*.rsp.gz``.

The response grid is 3072 pixels (HEALPix nside=16, RING ordering) in
**spacecraft** coordinates, split across two directories by hemisphere
(``th0_90`` for pixels 0-1567, ``th91_180`` for 1568-3071). That is ~1 GB in
total, so nothing here ever bulk-downloads or holds the whole grid in memory:
:class:`BurstCubeResponseGrid` loads pixels lazily with an LRU cache, and
:class:`BurstCubeRspFinder` downloads only the (typically 4-per-detector)
pixels actually needed for one direction.

Per-pixel filenames (e.g. ``bccs0_px0000t003p045_20240417v211.rsp.gz``) embed
a *rounded* theta/phi label that is not reconstructible from the pixel number
for every pixel (82 of 3072 sit on a half-integer phi from float noise), and
the date/version fields change when CALDB republishes. So **the pixel number
is the only key in the public API here**; filenames are resolved privately,
never constructed from theta/phi, and never hardcoded into a static index.
"""
import re
from functools import lru_cache
from pathlib import Path
from typing import Optional, Union

import astropy.units as u
import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord

from gdt.core.data_primitives import ResponseMatrix
from gdt.core.heasarc import Http
from gdt.core.response import Rsp

from . import caldb
from .headers import RspHeaders

__all__ = ['BurstCubeResponseGrid', 'BurstCubeRsp', 'BurstCubeRspFinder']

#: The HEALPix resolution of the response grid (RING ordering).
NSIDE = 16

#: Total number of pixels in the grid (``12 * NSIDE**2``).
NUM_PIXELS = hp.nside2npix(NSIDE)

#: The 64->16 channel regroup, verified against the real CS0 pixel-0 file
#: (spec section 21): 17 edge indices for 16 channels. A 16-entry list
#: (omitting the trailing 63) silently yields 15 channels while still
#: preserving the folded total -- :meth:`BurstCubeRsp.to_cbd` asserts
#: ``num_chans == 16`` to catch that class of mistake immediately.
EDGE_INDICES_64_TO_16 = np.array(
    [0, 1, 6, 11, 16, 21, 26, 31, 35, 39, 43, 47, 51, 55, 59, 63, 64])

#: The CALDB root for response files, under the same burstcube/csa release
#: tree as the other CALDB products (see gdt.missions.burstcube.caldb).
RESPONSE_REMOTE_ROOT = 'https://heasarc.gsfc.nasa.gov/FTP/caldb/data/burstcube/csa/cpf/response/'

#: The default local cache for downloaded response files (never the 233 kB
#: bundled CALDB directory -- these are not bundled, per the plugin spec).
DEFAULT_RESPONSE_CACHE_DIR = Path.home() / '.gdt' / 'burstcube' / 'responses'

_PIXEL_FROM_FILENAME = re.compile(r'px(\d{4})')


def _hemisphere(pixel: int) -> str:
    """The response-grid hemisphere directory for a pixel number.

    Args:
        pixel (int): The HEALPix pixel number, 0-3071

    Returns:
        (str): ``'th0_90'`` for pixels 0-1567, ``'th91_180'`` for 1568-3071
    """
    return 'th0_90' if pixel < NUM_PIXELS // 2 else 'th91_180'


def _resolve_pixel_path(detector: str, pixel: int,
                        cache_dir: Optional[Union[str, Path]] = None,
                        local_dir: Optional[Union[str, Path]] = None) -> Path:
    """Resolve a (detector, pixel) pair to a local ``.rsp`` file path,
    downloading it if necessary. This is the one place that knows how to
    find a response file; nothing else in this module (or its public API)
    ever constructs or depends on the filename pattern.

    Resolution order:

    1. Glob ``local_dir`` (if given) for ``*px{pixel:04d}*.rsp*``.
    2. Glob the local cache directory the same way.
    3. Fetch and cache the remote directory listing for this
       detector/hemisphere (once per detector/hemisphere per process, via
       :func:`_remote_listing`), find the entry whose filename contains
       ``px{pixel:04d}``, and download it into the cache.

    Args:
        detector (str): The detector name, e.g. ``'CS0'``
        pixel (int): The HEALPix pixel number, 0-3071
        cache_dir (str or Path, optional): The local cache directory.
            Defaults to :data:`DEFAULT_RESPONSE_CACHE_DIR`.
        local_dir (str or Path, optional): An additional directory to check
            first, e.g. a directory of files a user already downloaded by
            hand.

    Returns:
        (Path): The local path to the resolved ``.rsp`` (or ``.rsp.gz``) file.
    """
    det_dir = detector.lower()
    hemisphere = _hemisphere(pixel)
    # the detector prefix ("bccs0_") is a stable, structural part of the
    # naming convention -- unlike the rounded theta/phi labels, it is not
    # something the spec warns us off relying on. Include it in the glob so
    # that a flat directory holding more than one detector's files (e.g. a
    # user's own download folder, or the data-driven test fixtures, which
    # are not organized into per-detector subdirectories) can't match the
    # wrong detector's file for the same pixel number.
    pattern = f'*{det_dir}_px{pixel:04d}*.rsp*'
    cache_dir = Path(cache_dir) if cache_dir is not None else DEFAULT_RESPONSE_CACHE_DIR
    pixel_cache_dir = cache_dir / det_dir / hemisphere

    search_dirs = [Path(local_dir)] if local_dir is not None else []
    search_dirs.append(pixel_cache_dir)
    for directory in search_dirs:
        if directory.is_dir():
            matches = sorted(directory.glob(pattern))
            if matches:
                return matches[0]

    # not found locally: fetch (and cache) the remote directory listing,
    # then download just the one file we need.
    listing = _remote_listing(det_dir, hemisphere)
    key = f'px{pixel:04d}'
    filename = next((name for name in listing if key in name), None)
    if filename is None:
        raise FileNotFoundError(
            f'No response file for {detector} pixel {pixel} found in the '
            f'remote {det_dir}/{hemisphere} listing.')

    pixel_cache_dir.mkdir(parents=True, exist_ok=True)
    dest = pixel_cache_dir / filename
    url = f'{RESPONSE_REMOTE_ROOT}{det_dir}/{hemisphere}/{filename}'
    _download_url(url, pixel_cache_dir)
    return dest


@lru_cache(maxsize=None)
def _remote_listing(det_dir: str, hemisphere: str) -> tuple:
    """Fetch and cache (once per detector/hemisphere, for the life of the
    process) the remote directory listing of ``.rsp.gz`` filenames.

    Args:
        det_dir (str): The lowercase detector directory name, e.g. ``'cs0'``
        hemisphere (str): ``'th0_90'`` or ``'th91_180'``

    Returns:
        (tuple of str): The filenames in that remote directory.
    """
    http = Http(url=RESPONSE_REMOTE_ROOT)
    return tuple(http.ls(f'{det_dir}/{hemisphere}'))


def _download_url(url: str, dest_dir: Path):
    """Download one response file. Split out as its own function so tests
    can monkeypatch just the network step.

    Args:
        url (str): The remote URL
        dest_dir (Path): The local directory to download into
    """
    Http(url=RESPONSE_REMOTE_ROOT).download_url(url, dest_dir, verbose=False)


class BurstCubeRsp(Rsp):
    """A single BurstCube detector response matrix (DRM), for one HEALPix
    pixel of one detector. **Native 64 channels by default** -- see
    :meth:`to_cbd` to regroup to CBD's 16 channels.
    """

    def __init__(self):
        super().__init__()
        self._pixel = None

    @property
    def pixel(self):
        """(int or None): The HEALPix pixel number (nside=16, RING) this DRM
        was computed for, in spacecraft coordinates. None for a DRM that was
        not read from (or explicitly tagged with) a single grid pixel, e.g.
        an interpolated result from :class:`BurstCubeResponseGrid`.
        """
        return self._pixel

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a single BurstCube ``.rsp`` file.

        Args:
            file_path (str): The file path of the FITS file

        Returns:
            (:class:`BurstCubeRsp`)
        """
        obj = super().open(file_path, **kwargs)

        hdrs = [hdu.header for hdu in obj.hdulist]
        headers = RspHeaders.from_headers(hdrs)

        ebounds_idx = obj.hdu_index_from_name('EBOUNDS')
        specresp_idx = obj.hdu_index_from_name('SPECRESP MATRIX')

        # unlike GBM's RspHeaders, BurstCube .rsp files carry no NUMEBINS
        # keyword, so the number of photon bins is simply the row count.
        num_ebins = obj.column(specresp_idx, 'ENERG_LO').size
        num_chans = headers['EBOUNDS']['DETCHANS']
        fchan = obj.column(specresp_idx, 'F_CHAN')
        nchan = obj.column(specresp_idx, 'N_CHAN')
        matrix = cls._decompress_drm(obj.column(specresp_idx, 'MATRIX'),
                                     num_ebins, num_chans, fchan, nchan)

        # The energy columns are big-endian float32 in the file. ResponseMatrix
        # keeps whatever dtype it is handed, and photon_bin_centroids is the
        # geometric mean sqrt(emin * emax) -- evaluated in float32 that loses
        # ~2e-4 keV at the top of the 10-5000 keV range, which carries through
        # fold_spectrum as a ~2.5e-9 relative error. Widening to native float64
        # here costs 256 values and removes the loss (and avoids propagating a
        # byte-swapped dtype through every downstream array).
        to_f64 = lambda col: np.asarray(obj.column(specresp_idx, col),
                                        dtype=np.float64)
        to_f64_eb = lambda col: np.asarray(obj.column(ebounds_idx, col),
                                           dtype=np.float64)
        drm = ResponseMatrix(matrix, to_f64('ENERG_LO'), to_f64('ENERG_HI'),
                             to_f64_eb('E_MIN'), to_f64_eb('E_MAX'))

        detector = headers['EBOUNDS']['INSTRUME']
        pixel = headers['EBOUNDS']['PIXEL']

        obj.close()

        rsp = cls.from_data(drm, filename=obj.filename, headers=headers,
                            detector=detector)
        rsp._fchan = fchan
        rsp._nchan = nchan
        rsp._pixel = pixel
        return rsp

    def to_cbd(self):
        """Regroup this DRM's 64 native channels to CBD's 16 channels.

        This does two things, both required (spec section 21):

        1. Regroups the channel axis using the verified 17-entry
           :data:`EDGE_INDICES_64_TO_16`, and asserts the result actually
           has 16 channels (a truncated 16-entry list would silently give
           15 while still preserving the folded total).
        2. Replaces the regrouped channel boundaries -- which
           :meth:`~gdt.core.data_primitives.ResponseMatrix.rebin` inherits
           from this DRM's own (older, rounded, detector-independent)
           internal ``EBOUNDS`` -- with this detector's CALDB ``eb16``
           energies. The two disagree meaningfully at the top channel (the
           .rsp says 1290-5000 keV; CALDB says ~1298.79-1648.51 keV for CS0),
           so skipping this step would silently mislabel the top channel's
           energy range.

        Returns:
            (:class:`BurstCubeRsp`)
        """
        regrouped = self.drm.rebin(edge_indices=EDGE_INDICES_64_TO_16)
        if regrouped.num_chans != 16:
            raise RuntimeError(
                f'Expected 16 channels after the CBD regroup, got '
                f'{regrouped.num_chans}; EDGE_INDICES_64_TO_16 may be wrong.')

        eb16 = caldb.ebounds(self.detector, 16)
        final_drm = ResponseMatrix(regrouped.matrix,
                                   regrouped.photon_bins.low_edges(),
                                   regrouped.photon_bins.high_edges(),
                                   eb16.low_edges(), eb16.high_edges())

        rsp = BurstCubeRsp.from_data(final_drm, filename=self.filename,
                                     start_time=self.tstart,
                                     stop_time=self.tstop,
                                     trigger_time=self.trigtime,
                                     headers=self.headers,
                                     detector=self.detector)
        rsp._pixel = self.pixel
        return rsp

    def _build_headers(self, num_chans, num_ebins):
        headers = self.headers.copy()
        headers['EBOUNDS']['DETCHANS'] = num_chans
        headers['SPECRESP MATRIX']['DETCHANS'] = num_chans
        return headers

    @staticmethod
    def _decompress_drm(matrix, num_photon_bins, num_channels, fchan, nchan):
        """Decompress a DRM using the OGIP ``F_CHAN``/``N_CHAN``/``N_GRP``
        scheme.

        Note:
            **BurstCube's ``F_CHAN`` is 0-indexed** (``TLMIN4=0`` in the
            header, and verified directly against the real CS0 pixel-0
            file), unlike GBM's 1-indexed convention -- do not subtract 1
            here. Every sampled BurstCube ``.rsp`` file has ``N_GRP=1`` for
            every row (a fixed, non-variable-length ``F_CHAN``/``N_CHAN``
            column format that cannot represent more than one group per
            row), so each row is decompressed as a single contiguous block.

        Args:
            matrix (numpy.ndarray): The (num_photon_bins, num_channels)
                ``MATRIX`` column
            num_photon_bins (int): The number of photon (row) bins
            num_channels (int): The number of detector channels
            fchan (numpy.ndarray): The ``F_CHAN`` column (0-indexed start
                channel of the single group in each row)
            nchan (numpy.ndarray): The ``N_CHAN`` column (number of channels
                in the single group in each row)

        Returns:
            (numpy.ndarray)
        """
        drm = np.zeros((num_photon_bins, num_channels))
        for i in range(num_photon_bins):
            start = int(np.ravel(fchan[i])[0])
            width = int(np.ravel(nchan[i])[0])
            drm[i, start:start + width] = np.ravel(matrix[i])[:width]
        return drm


class BurstCubeResponseGrid:
    """The nside=16 RING HEALPix response grid, in spacecraft coordinates,
    for one or more BurstCube detectors. Pixels are loaded lazily (with an
    LRU cache) as they are requested -- this never holds all 3072 pixels'
    matrices in memory at once.

    Args:
        detectors (optional): One or more detector names (e.g. ``'CS0'``) or
            :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`
            members this grid serves. If omitted, all 4 detectors.
        cache_dir (str or Path, optional): The local cache directory for
            downloaded response files. Defaults to
            :data:`DEFAULT_RESPONSE_CACHE_DIR`.
        local_dir (str or Path, optional): An additional directory to check
            for already-downloaded response files before any network access.
        max_cached_pixels (int, optional): The number of (detector, pixel)
            DRMs to keep in the LRU cache at once. Default is 64 (16 get_drm
            calls' worth of 4-pixel interpolation, times 4 detectors).
    """

    def __init__(self, detectors=None, cache_dir=None, local_dir=None,
                max_cached_pixels: int = 64):
        if detectors is None:
            from .detectors import BurstCubeDetectors
            detectors = list(BurstCubeDetectors)
        elif not isinstance(detectors, (list, tuple)):
            detectors = [detectors]
        self._detectors = [getattr(d, 'name', d).upper() for d in detectors]
        self._cache_dir = cache_dir
        self._local_dir = local_dir
        self._load_pixel = lru_cache(maxsize=max_cached_pixels)(self._load_pixel_uncached)

    @property
    def detectors(self):
        """(list of str): The detector names this grid serves."""
        return list(self._detectors)

    def get_drm(self, det, az=None, zen=None, pix=None, skycoord=None,
               frame=None, interp=True):
        """Get a DRM for one detector, by spacecraft azimuth/zenith, by
        HEALPix pixel number, or by sky position.

        Args:
            det (str or :class:`~gdt.missions.burstcube.detectors.BurstCubeDetectors`):
                The detector, e.g. ``'CS0'``.
            az (float, optional): Spacecraft azimuth, in degrees. Requires
                `zen`.
            zen (float, optional): Spacecraft zenith angle, in degrees.
                Requires `az`.
            pix (int, optional): A HEALPix pixel number (0 to
                :data:`NUM_PIXELS` - 1). Bypasses interpolation: `interp` is
                ignored.
            skycoord (astropy.coordinates.SkyCoord, optional): A sky
                position. Requires `frame` to carry a valid attitude
                quaternion (see archive caveat #1); raises otherwise.
            frame (:class:`~gdt.missions.burstcube.frame.BurstCubeFrame`, optional):
                The spacecraft frame to transform `skycoord` into. Required
                if `skycoord` is given.
            interp (bool, optional): If True (the default), interpolate
                over the 4 nearest pixels with
                :func:`healpy.get_interp_weights` (matching
                ``bctools.io.InstrumentResponse.get_drm``). If False, use
                the nearest pixel (:func:`healpy.ang2pix`). Ignored if `pix`
                is given.

        Returns:
            (:class:`BurstCubeRsp`)
        """
        det_name = getattr(det, 'name', det).upper()

        if pix is not None:
            return self._load_pixel(det_name, int(pix))

        if skycoord is not None:
            az, zen = self._az_zen_from_skycoord(skycoord, frame)
        elif az is None or zen is None:
            raise ValueError('Provide (az, zen), pix, or (skycoord, frame).')

        theta = np.radians(zen)
        phi = np.radians(az)

        if not interp:
            pixel = int(hp.ang2pix(NSIDE, theta, phi))
            return self._load_pixel(det_name, pixel)

        pixels, weights = hp.get_interp_weights(NSIDE, theta, phi)
        matrices = [self._load_pixel(det_name, int(p)).drm.matrix for p in pixels]
        weighted_matrix = sum(w * m for w, m in zip(weights, matrices))

        template = self._load_pixel(det_name, int(pixels[0]))
        drm = ResponseMatrix(weighted_matrix, template.drm.photon_bins.low_edges(),
                             template.drm.photon_bins.high_edges(),
                             template.ebounds.low_edges(),
                             template.ebounds.high_edges())
        return BurstCubeRsp.from_data(drm, detector=det_name)

    def _load_pixel_uncached(self, det_name: str, pixel: int) -> BurstCubeRsp:
        path = _resolve_pixel_path(det_name, pixel, cache_dir=self._cache_dir,
                                   local_dir=self._local_dir)
        return BurstCubeRsp.open(path)

    @staticmethod
    def _az_zen_from_skycoord(skycoord: SkyCoord, frame) -> tuple:
        """Convert a sky position to spacecraft (az, zen), requiring a frame
        with attitude.

        Args:
            skycoord (astropy.coordinates.SkyCoord): The sky position
            frame: The spacecraft frame to transform into

        Returns:
            (float, float): (az, zen) in degrees
        """
        if frame is None or frame.quaternion is None:
            raise ValueError(
                'A SkyCoord requires a BurstCubeFrame with a valid attitude '
                'quaternion to convert to spacecraft coordinates. Per '
                'archive caveat #1, attitude is unavailable for most of the '
                'BurstCube mission; supply one via '
                'BurstCubeFrame.from_quaternion() (see gdt.missions.burstcube.frame).')
        local = skycoord.transform_to(frame)
        az = np.asarray(local.az.to(u.deg).value)
        zen = np.asarray((90.0 * u.deg - local.el.to(u.deg)).value)
        if az.size == 1:
            az, zen = az.item(), zen.item()
        return az, zen


class BurstCubeRspFinder:
    """Downloads only the response-grid pixels actually needed for one sky
    or spacecraft direction -- for one location that is 4 pixels per
    detector (~16 files, ~320 kB total for all 4 detectors), rather than the
    ~1 GB full grid. This is the entire point of this class: it never bulk-
    downloads.

    Args:
        cache_dir (str or Path, optional): The local cache directory to
            download into. Defaults to :data:`DEFAULT_RESPONSE_CACHE_DIR`.
    """

    def __init__(self, cache_dir=None):
        self._cache_dir = cache_dir

    def needed_pixels(self, az, zen, interp=True):
        """The pixel numbers needed for one direction.

        Args:
            az (float): Spacecraft azimuth, in degrees
            zen (float): Spacecraft zenith angle, in degrees
            interp (bool, optional): If True (the default), the 4
                interpolation pixels; if False, just the 1 nearest pixel.

        Returns:
            (list of int)
        """
        theta, phi = np.radians(zen), np.radians(az)
        if not interp:
            return [int(hp.ang2pix(NSIDE, theta, phi))]
        pixels, _ = hp.get_interp_weights(NSIDE, theta, phi)
        return [int(p) for p in pixels]

    def get_pixels(self, detectors, pixels):
        """Download (or resolve from a local cache) a set of pixels for a
        set of detectors.

        Args:
            detectors (list): The detector names, e.g. ``['CS0', 'CS1']``
            pixels (list of int): The HEALPix pixel numbers to fetch

        Returns:
            (dict): ``{detector: {pixel: Path}}``
        """
        return {det: {pix: _resolve_pixel_path(det, pix, cache_dir=self._cache_dir)
                     for pix in pixels}
               for det in detectors}

    def get_drms(self, detectors, az, zen, interp=True):
        """Download only the pixels needed for one direction, across one or
        more detectors, and return the opened DRMs.

        Args:
            detectors (list): The detector names, e.g. ``['CS0', 'CS1']``
            az (float): Spacecraft azimuth, in degrees
            zen (float): Spacecraft zenith angle, in degrees
            interp (bool, optional): If True (the default), fetch the 4
                interpolation pixels; if False, just the nearest pixel.

        Returns:
            (dict): ``{detector: {pixel: BurstCubeRsp}}``
        """
        pixels = self.needed_pixels(az, zen, interp=interp)
        paths = self.get_pixels(detectors, pixels)
        return {det: {pix: BurstCubeRsp.open(path) for pix, path in by_pixel.items()}
               for det, by_pixel in paths.items()}
