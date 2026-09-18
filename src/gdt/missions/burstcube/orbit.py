"""BurstCube orbit/ephemeris data: the spacecraft position and velocity in
J2000 ECI coordinates, reconstructed from TLEs (NORAD 59562), from the
archive's ``auxil/bcYYMMDD.hk.gz`` files (single ``ORBIT`` extension).

Per archive caveat #2, the orbit file spans the full on-orbit lifetime
(including days with no science data at all), and its time is independent of
-- and more accurate than -- the instrument clock; there may be unknown
offsets between orbit time and CBD/HK/TTE MET.

There is no attitude information in this file (``INSTRUME='CS'``, not a
detector or ``'CSA'``): the resulting :class:`~gdt.missions.burstcube.frame.BurstCubeFrame`
carries ``obsgeoloc``/``obsgeovel`` only, with no quaternion. It is usable
for Earth-visibility and geocenter calculations, but not for converting sky
positions to/from BurstCube coordinates -- combine it with a user-supplied
quaternion via :meth:`~gdt.missions.burstcube.frame.BurstCubeFrame.from_quaternion`
for that (see archive caveat #1).
"""
import numpy as np
import astropy.coordinates.representation as r
import astropy.table as table
import astropy.io.fits as fits
import astropy.units as u

from gdt.core.coords.spacecraft import SpacecraftFrameModelMixin
from gdt.core.file import FitsFileContextManager

from .detectors import BurstCubeDetectors
from .frame import BurstCubeFrame
from .headers import OrbitHeaders
from .time import Time

__all__ = ['BurstCubeOrbit']


class BurstCubeOrbit(SpacecraftFrameModelMixin, FitsFileContextManager):
    """Reader for a BurstCube orbit/ephemeris file (the ``ORBIT`` extension
    of ``auxil/bcYYMMDD.hk.gz``).
    """

    def _reorder_bytes(self, arr):
        """Reorder bytes for compatibility across numpy's byte-order API,
        which changed between numpy 1.x and 2.x. FITS data is big-endian,
        and some downstream operations (e.g. scipy interpolation) require
        native byte order.

        Args:
            arr (numpy.ndarray): The big-endian array read from the FITS file

        Returns:
            (numpy.ndarray): The same data in native byte order
        """
        if np.__version__ >= '2.0.0':
            return arr.view(arr.dtype.newbyteorder()).byteswap()
        return arr.byteswap().newbyteorder()

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a BurstCube orbit FITS file.

        Args:
            file_path (str): The file path of the FITS file

        Returns:
            (:class:`BurstCubeOrbit`)
        """
        obj = super().open(file_path, **kwargs)
        hdrs = [hdu.header for hdu in obj.hdulist]
        obj._headers = OrbitHeaders.from_headers(hdrs)
        return obj

    def get_spacecraft_frame(self) -> BurstCubeFrame:
        """Retrieve the spacecraft position and velocity as a
        :class:`~gdt.missions.burstcube.frame.BurstCubeFrame`. No quaternion
        is set (this file carries no attitude information).

        Returns:
            (:class:`~gdt.missions.burstcube.frame.BurstCubeFrame`)
        """
        orbit_idx = self.hdu_index_from_name('ORBIT')
        frame = BurstCubeFrame(
            obstime=Time(self._reorder_bytes(self.column(orbit_idx, 'TIME')),
                        format='burstcube'),
            obsgeoloc=r.CartesianRepresentation(
                x=self._reorder_bytes(self.column(orbit_idx, 'X')),
                y=self._reorder_bytes(self.column(orbit_idx, 'Y')),
                z=self._reorder_bytes(self.column(orbit_idx, 'Z')),
                unit=u.km),
            obsgeovel=r.CartesianRepresentation(
                x=self._reorder_bytes(self.column(orbit_idx, 'Vx')) * u.km / u.s,
                y=self._reorder_bytes(self.column(orbit_idx, 'Vy')) * u.km / u.s,
                z=self._reorder_bytes(self.column(orbit_idx, 'Vz')) * u.km / u.s,
                unit=u.km / u.s),
            detectors=BurstCubeDetectors)
        return frame

    @classmethod
    def merge(cls, orbit1, orbit2):
        """Merge two BurstCubeOrbit objects into a single object, e.g. for
        consecutive days. The order of the inputs does not matter, and any
        duplicate time entries are removed.

        Args:
            orbit1 (:class:`BurstCubeOrbit`): The first orbit object
            orbit2 (:class:`BurstCubeOrbit`): The second orbit object

        Returns:
            (:class:`BurstCubeOrbit`)
        """
        idx1 = orbit1.hdu_index_from_name('ORBIT')
        idx2 = orbit2.hdu_index_from_name('ORBIT')
        table1 = table.Table(orbit1._hdulist[idx1].data)
        table2 = table.Table(orbit2._hdulist[idx2].data)

        if orbit1.headers['ORBIT']['TSTART'] < orbit2.headers['ORBIT']['TSTART']:
            new_table = table.vstack([table1, table2])
            ref_obj, other_obj = orbit1, orbit2
        else:
            new_table = table.vstack([table2, table1])
            ref_obj, other_obj = orbit2, orbit1

        new_table = table.unique(new_table, keys='TIME')
        new_table.sort('TIME')

        prihdr = ref_obj._hdulist[0].header.copy()
        prihdr['DATE-END'] = other_obj._hdulist[0].header['DATE-END']

        orbit_hdr = ref_obj._hdulist[idx1].header.copy()
        orbit_hdr['DATE-END'] = other_obj._hdulist[idx2].header['DATE-END']
        orbit_hdr['TSTOP'] = other_obj._hdulist[idx2].header['TSTOP']

        hdulist = fits.HDUList()
        hdulist.append(fits.PrimaryHDU(header=prihdr))
        hdulist.append(fits.BinTableHDU(data=new_table, header=orbit_hdr,
                                        name='ORBIT'))

        obj = cls()
        obj._hdulist = hdulist
        obj._headers = OrbitHeaders.from_headers(
            [hdu.header for hdu in obj.hdulist])
        return obj
