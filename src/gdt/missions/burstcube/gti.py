"""Reader and set operations for the BurstCube trend GTI files
(``trend/gti_binning``, ``trend/gti_poscnt``, ``trend/gti_saa``).

:class:`~gdt.core.data_primitives.Gti` already provides ``intersection`` and
a gap-merging ``merge`` (its own ``insert`` coalesces overlapping intervals,
so ``merge`` is a proper union). :class:`BurstCubeGTI` wraps the file format
and adds, as class methods, the :meth:`~BurstCubeGTI.complement` that
gdt-core does not provide plus a :meth:`~BurstCubeGTI.apply_to` for cutting
a :class:`~gdt.core.phaii.Phaii` or :class:`~gdt.core.tte.PhotonList` down to
a GTI -- e.g. to reproduce a CBD ``_cl`` file's cut starting from its ``_uf``
counterpart and the matching ``trend/gti_*`` file.

All four operations take and return plain
:class:`~gdt.core.data_primitives.Gti` objects, so they work on any GTI --
one read from a trend file, a data product's own ``.gti``, or one built by
hand -- not only on the contents of a :class:`BurstCubeGTI` file.
"""
import numpy as np

from gdt.core.data_primitives import Gti
from gdt.core.file import FitsFileContextManager

from .headers import GTIHeaders

__all__ = ['BurstCubeGTI']


class BurstCubeGTI(FitsFileContextManager):
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
            (:class:`BurstCubeGTI`)
        """
        obj = super().open(file_path, **kwargs)
        hdrs = [hdu.header for hdu in obj.hdulist]
        obj._headers = GTIHeaders.from_headers(hdrs)
        return obj

    @classmethod
    def intersect(cls, gti1: Gti, gti2: Gti) -> Gti:
        """The intersection of two GTIs: time covered by both.

        Args:
            gti1 (:class:`~gdt.core.data_primitives.Gti`): A GTI
            gti2 (:class:`~gdt.core.data_primitives.Gti`): Another GTI

        Returns:
            (:class:`~gdt.core.data_primitives.Gti`)
        """
        return Gti.intersection(gti1, gti2)

    @classmethod
    def union(cls, gti1: Gti, gti2: Gti) -> Gti:
        """The union of two GTIs: time covered by either. Overlapping or
        touching intervals are coalesced into one.

        Args:
            gti1 (:class:`~gdt.core.data_primitives.Gti`): A GTI
            gti2 (:class:`~gdt.core.data_primitives.Gti`): Another GTI

        Returns:
            (:class:`~gdt.core.data_primitives.Gti`)
        """
        return Gti.merge(gti1, gti2)

    @classmethod
    def complement(cls, gti: Gti, tstart: float, tstop: float) -> Gti:
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
        for low, high in sorted(gti.as_list()):
            low, high = max(low, tstart), min(high, tstop)
            if low > high:
                continue
            if low > cursor:
                gaps.append((cursor, low))
            cursor = max(cursor, high)
        if cursor < tstop:
            gaps.append((cursor, tstop))

        if not gaps:
            return None
        return Gti.from_list(gaps)

    @classmethod
    def apply_to(cls, gti: Gti, data_obj):
        """Cut a :class:`~gdt.core.phaii.Phaii` or
        :class:`~gdt.core.tte.PhotonList` down to a GTI, e.g. to reproduce a
        ``_cl`` file's cut starting from its ``_uf`` counterpart.

        This exists because gdt-core's own ``slice_time`` cannot be handed a
        real GTI directly; see the comments in the body for what goes wrong
        and why.

        Args:
            gti (:class:`~gdt.core.data_primitives.Gti`): The GTI to apply
            data_obj (:class:`~gdt.core.phaii.Phaii` or :class:`~gdt.core.tte.PhotonList`):
                The data to cut down

        Returns:
            (:class:`~gdt.core.phaii.Phaii` or :class:`~gdt.core.tte.PhotonList`):
            A new object of the same type as ``data_obj``, sliced to the GTI.

        Raises:
            ValueError: If no interval of ``gti`` overlaps the data at all,
                since there is no meaningful empty object to return.
        """
        data = data_obj.data
        if hasattr(data, 'tstart'):
            # binned (Phaii): keep intervals that touch at least one bin
            def selects_data(low, high):
                return bool(np.any((data.tstop > low) & (data.tstart < high)))
        else:
            # unbinned (PhotonList): keep intervals containing >= 1 event
            def selects_data(low, high):
                return bool(np.any((data.times >= low) & (data.times <= high)))

        # Problem 1: most intervals of a real GTI select nothing. A trend GTI
        # spans the whole mission -- the SAA one has 586 intervals over five
        # months -- while a single CBD or TTE file covers minutes to hours.
        # slice_time builds Gti.from_list([segment.time_range]) per requested
        # range, and an empty segment's time_range is None, so it raises
        # TypeError from inside the primitive rather than returning an empty
        # result. Dropping those intervals here is the fix. Filtering on the
        # file's overall time span instead would not be enough: BurstCube data
        # is gappy, so an interval can sit inside the span and still contain
        # no bins.
        overlapping = [(low, high) for low, high in gti.as_list()
                       if selects_data(low, high)]

        if not overlapping:
            data_start, data_stop = data_obj.time_range
            intervals = gti.as_list()
            raise ValueError(
                f'No interval of the GTI contains any data. The data spans '
                f'{data_start} to {data_stop}; the GTI covers '
                f'{intervals[0][0]} to {intervals[-1][1]} in '
                f'{len(intervals)} intervals. Check that the GTI and the '
                'data are for the same observation and detector.')

        if len(overlapping) == 1:
            return data_obj.slice_time(overlapping)

        # Problem 2: disjoint intervals empty slice_time's GTI accumulator.
        # Phaii.slice_time intersects each sliced segment's own range into a
        # running total, which goes empty as soon as two segments are
        # disjoint, and the next iteration raises on the empty result -- even
        # though that accumulated value is then discarded, since from_data is
        # handed self.gti instead. Slicing one interval per call runs that
        # loop exactly once, so it never empties; merge recombines the
        # pieces. It fails on a CBD _uf file's own 4-segment STDGTI, not just
        # on the mission-long trend files.
        pieces = [data_obj.slice_time([interval]) for interval in overlapping]
        return type(data_obj).merge(pieces)
