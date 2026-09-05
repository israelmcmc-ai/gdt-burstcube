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
from .time import Time

__all__ = ['BurstCubeAttitude']

#: Approximate 1-sigma pointing reconstruction error, in degrees, for each of
#: the 3 rows of the attitude file. These are NOT present in the FITS file
#: itself: they come from ``trend/README`` and archive caveat #1, and are
#: recorded here only for reference/documentation. They are approximate
#: ("~10 deg" in the source), so treat them as order-of-magnitude only.
APPROXIMATE_POINTING_ERROR_DEG = (10.0, 12.0, 10.5)


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
