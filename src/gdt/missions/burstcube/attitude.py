"""BurstCube attitude data: a standalone reader for
``trend/attitude/bc_csa_att.fits``, which holds exactly 3 manually
reconstructed attitude epochs across the entire archive.

This is deliberately **not** folded into :class:`~gdt.missions.burstcube.frame.BurstCubeFrame`
and does **not** interpolate or extrapolate across the mission: with only 3
widely-spaced epochs (none of which coincide with a TTE day) and no
information about the attitude in between, any interpolation would be
fabricating data. Use :meth:`BurstCubeAttitude.frame` to get an exact
:class:`~gdt.missions.burstcube.frame.BurstCubeFrame` for one of the 3
epochs, and see archive caveat #1 for the many other times attitude simply
is not available -- for which the documented path is a user-supplied
quaternion via :meth:`~gdt.missions.burstcube.frame.BurstCubeFrame.from_quaternion`.
"""
from gdt.core.coords import Quaternion
from gdt.core.file import FitsFileContextManager

from .detectors import BurstCubeDetectors
from .frame import BurstCubeFrame
from .headers import AttitudeHeaders
from .time import Time, check_met_epoch

__all__ = ['BurstCubeAttitude']


class BurstCubeAttitude(FitsFileContextManager):
    """Reader for the BurstCube attitude file. Exactly 3 rows exist in the
    whole archive; this class exposes them directly rather than offering any
    interpolation between them.
    """

    @property
    def num_rows(self):
        """(int): The number of attitude epochs (always 3, in the real
        archive file)."""
        return self.time.size

    @property
    def time(self):
        """(astropy.time.Time): The MET of each attitude epoch."""
        attitude_idx = self.hdu_index_from_name('ATTITUDE')
        return Time(self.column(attitude_idx, 'TIME'), format='burstcube')

    @property
    def quaternion(self):
        """(:class:`~gdt.core.coords.Quaternion`): The attitude quaternion(s)
        for each epoch. ``QPARAM`` is stored scalar-last (x, y, z, w) -- see
        :mod:`gdt.missions.burstcube.frame` for how that was established --
        so no reordering is needed.
        """
        attitude_idx = self.hdu_index_from_name('ATTITUDE')
        return Quaternion(self.column(attitude_idx, 'QPARAM'))

    @property
    def pointing(self):
        """(numpy.ndarray): The (RA, Dec, roll) of the spacecraft +z axis
        (the instrument boresight) in ICRS, in degrees, for each epoch, as
        recorded in the file's own ``POINTING`` column.
        """
        attitude_idx = self.hdu_index_from_name('ATTITUDE')
        return self.column(attitude_idx, 'POINTING')

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open the BurstCube attitude FITS file.

        Args:
            file_path (str): The file path of the FITS file

        Returns:
            (:class:`BurstCubeAttitude`)
        """
        obj = super().open(file_path, **kwargs)
        hdrs = [hdu.header for hdu in obj.hdulist]
        obj._headers = AttitudeHeaders.from_headers(hdrs)
        check_met_epoch(obj._headers['ATTITUDE'])
        return obj

    def frame(self, row: int) -> BurstCubeFrame:
        """Build an exact :class:`~gdt.missions.burstcube.frame.BurstCubeFrame`
        for one of the reconstructed attitude epochs. This is only valid at
        the epoch itself -- there is no interpolation to any other time.

        Args:
            row (int): The row index of the attitude epoch, 0 to
                :attr:`num_rows` - 1.

        Returns:
            (:class:`~gdt.missions.burstcube.frame.BurstCubeFrame`)
        """
        return BurstCubeFrame(obstime=self.time[row],
                              quaternion=self.quaternion[row],
                              detectors=BurstCubeDetectors)
