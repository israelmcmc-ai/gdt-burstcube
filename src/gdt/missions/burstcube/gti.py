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
"""Readers and set operations for the BurstCube trend GTI files
(``trend/gti_binning``, ``trend/gti_poscnt``, ``trend/gti_saa``).

:class:`~gdt.core.data_primitives.Gti` already provides ``intersection`` and
a gap-merging ``merge`` (its own ``insert`` coalesces overlapping intervals,
so ``merge`` is a proper union); this module adds the :func:`complement`
that gdt-core does not provide, plus :func:`apply_to` for cutting a
:class:`~gdt.core.phaii.Phaii` or :class:`~gdt.core.tte.PhotonList` down to
a GTI -- e.g. to reproduce a CBD ``_cl`` file's cut starting from its ``_uf``
counterpart and the matching ``trend/gti_*`` file.
"""
from gdt.core.data_primitives import Gti
from gdt.core.file import FitsFileContextManager

from .headers import GtiHeaders

__all__ = ['BurstCubeGti', 'apply_to', 'complement', 'intersect', 'union']


class BurstCubeGti(FitsFileContextManager):
    """Reader for one of the standalone BurstCube trend GTI files. Each file
    is a single ``STDGTI`` extension of ``START``/``STOP`` pairs.
    """

    @property
    def gti(self):
        """(:class:`~gdt.core.data_primitives.Gti`): The good time intervals
        in this file."""
        idx = self.hdu_index_from_name('STDGTI')
        return Gti.from_bounds(self.column(idx, 'START'),
                               self.column(idx, 'STOP'))

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a BurstCube trend GTI file.

        Args:
            file_path (str): The file path of the FITS file

        Returns:
            (:class:`BurstCubeGti`)
        """
        obj = super().open(file_path, **kwargs)
        hdrs = [hdu.header for hdu in obj.hdulist]
        obj._headers = GtiHeaders.from_headers(hdrs)
        return obj


def intersect(gti1: Gti, gti2: Gti) -> Gti:
    """The intersection of two GTIs: time covered by both.

    Args:
        gti1 (:class:`~gdt.core.data_primitives.Gti`): A GTI
        gti2 (:class:`~gdt.core.data_primitives.Gti`): Another GTI

    Returns:
        (:class:`~gdt.core.data_primitives.Gti`)
    """
    return Gti.intersection(gti1, gti2)


def union(gti1: Gti, gti2: Gti) -> Gti:
    """The union of two GTIs: time covered by either. Overlapping or
    touching intervals are coalesced into one.

    Args:
        gti1 (:class:`~gdt.core.data_primitives.Gti`): A GTI
        gti2 (:class:`~gdt.core.data_primitives.Gti`): Another GTI

    Returns:
        (:class:`~gdt.core.data_primitives.Gti`)
    """
    return Gti.merge(gti1, gti2)


def complement(gti: Gti, tstart: float, tstop: float) -> Gti:
    """The complement of a GTI within a bounding time range: the gaps.

    Args:
        gti (:class:`~gdt.core.data_primitives.Gti`): The GTI to complement
        tstart (float): The start of the bounding time range
        tstop (float): The end of the bounding time range

    Returns:
        (:class:`~gdt.core.data_primitives.Gti` or None): None if the GTI
        entirely covers ``[tstart, tstop]`` with no gaps.
    """
    gaps = []
    cursor = tstart
    for lo, hi in sorted(gti.as_list()):
        lo, hi = max(lo, tstart), min(hi, tstop)
        if lo > hi:
            continue
        if lo > cursor:
            gaps.append((cursor, lo))
        cursor = max(cursor, hi)
    if cursor < tstop:
        gaps.append((cursor, tstop))

    if not gaps:
        return None
    return Gti.from_list(gaps)


def apply_to(gti: Gti, data_obj):
    """Cut a :class:`~gdt.core.phaii.Phaii` or
    :class:`~gdt.core.tte.PhotonList` down to a GTI, e.g. to reproduce a
    ``_cl`` file's cut starting from its ``_uf`` counterpart.

    Args:
        gti (:class:`~gdt.core.data_primitives.Gti`): The GTI to apply
        data_obj (:class:`~gdt.core.phaii.Phaii` or :class:`~gdt.core.tte.PhotonList`):
            The data to cut down

    Returns:
        (:class:`~gdt.core.phaii.Phaii` or :class:`~gdt.core.tte.PhotonList`):
        A new object of the same type as ``data_obj``, sliced to the GTI.
    """
    return data_obj.slice_time(gti.as_list())
