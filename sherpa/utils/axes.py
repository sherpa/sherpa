#  Copyright (C) 2017-2025
#  Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

"""Represent axis and axes.

"""

# How much of sherpa.models.regrid.EvaluationSpace1D/2D would we want
# here? There are features of the model that could be useful, but
# for now leave as is.

from __future__ import annotations

from abc import ABCMeta, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, overload

import numpy as np

from .err import DataErr
from .types import ArrayType


@overload
def _to_readable_array(x: None) -> None:
    ...

@overload
def _to_readable_array(x: ArrayType) -> np.ndarray:
    ...

def _to_readable_array(x: ArrayType | None) -> np.ndarray | None:
    """Convert x into a ndarray that can not be edited (or is None)."""

    if x is None:
        return None

    x = np.asarray(x).copy()
    if not np.iterable(x):
        raise DataErr("notanarray")

    if TYPE_CHECKING:
        assert isinstance(x, np.ndarray)

    if x.ndim != 1:
        raise DataErr("not1darray")

    x.setflags(write=False)
    return x


class Axis(metaclass=ABCMeta):
    """Represent an axis of a N-D object.

    .. versionchanged:: 4.18.0
       This class is now defined in sherpa.utils.axes rather than
       sherpa.models.regrid.

    """

    # This is set when the data is set
    _is_ascending = None

    @property
    def is_empty(self) -> bool:
        """Is the axis empty?

        An empty axis is either set to `None` or a zero-element
        sequence.
        """
        return self.size == 0

    @property
    @abstractmethod
    def is_integrated(self) -> bool:
        """Is the axis integrated?"""
        pass

    @property
    def is_ascending(self) -> bool:
        """Is the axis ascending?

        The axis is ascending if the elements in `lo` are sorted in
        ascending order.  Only the first and last elements are
        checked, and it is assumed that the elements are sorted.
        """
        if self.is_empty:
            raise DataErr("Axis is empty or has a size of 0")

        if TYPE_CHECKING:
            assert self._is_ascending is not None

        return self._is_ascending

    @property
    @abstractmethod
    def start(self) -> float:
        """The starting point (lowest value) of the data axis."""
        pass

    @property
    @abstractmethod
    def end(self) -> float:
        """The ending point of the data axis

        If the data axis is ascending the end boundary is the last
        element of the `hi` array when the axis is integrated,
        otherwise it's the last element of `lo`. Conversely, for
        descending axes, the last element is either the first element
        of the `hi` array or of the `lo` array, depending on whether
        the axis is integrated or not, respectively.
        """
        pass

    @property
    @abstractmethod
    def size(self) -> int:
        """The size of the axis."""
        pass

    def overlaps(self, other: Axis) -> bool:
        """Check if this axis overlaps with another.

        Parameters
        ----------
        other : Axis instance
            The axis to compare to.

        Returns
        -------
        bool
            True if they overlap, False if not
        """
        # Could apply this check but this is not expected to be
        # called directly, so leave for now.
        #
        # if not isinstance(other, Axis):
        #     raise TypeError("other argument must be an axis")

        num = max(0, min(self.end, other.end) - max(self.start, other.start))
        return bool(num != 0)


class PointAxis(Axis):
    """Represent a point (not integrated) axis of a N-D object.

    The length can not be changed once set.

    .. versionchanged:: 4.18.0
       This class is now defined in sherpa.utils.axes rather than
       sherpa.models.regrid.

    Parameters
    ----------
    x : array_like or None
        The starting point of the axis. If `None` or `[]` then the
        data axis is said to be empty. The axis can be in ascending or
        descending order but this is not checked.

    """

    def __init__(self, x: ArrayType | None) -> None:
        self._x = _to_readable_array(x)
        if self._x is None:
            return

        nx = len(self._x)
        if nx > 0:
            self._is_ascending = self._x[-1] > self._x[0]

    @property
    def x(self) -> np.ndarray | None:
        """Read-only access to the axis data."""
        return self._x

    @property
    def is_integrated(self) -> bool:
        return False

    # TODO: what happens if the axis is empty?
    @property
    def start(self) -> float:
        if self.is_ascending:
            return self.x[0]

        return self.x[-1]

    # TODO: what happens if the axis is empty?
    @property
    def end(self) -> float:
        if self.is_ascending:
            return self.x[-1]

        return self.x[0]

    @property
    def size(self) -> int:
        if self.x is None:
            return 0

        return self.x.size


class IntegratedAxis(Axis):
    """Represent an integrated axis of a N-D object.

    .. versionchanged:: 4.18.0
       This class is now defined in sherpa.utils.axes rather than
       sherpa.models.regrid.

    Parameters
    ----------
    lo : array_like or None
        The starting point of the axis.  If `lo` is `None` or `[]`
        then the data axis is said to be empty. The axis can be in
        ascending or descending order.
    hi : array_like or None
        The ending point of the axis. The number of elements must match `lo` (either `None`
        or a sequence of the same size). Each element is expected to
        be larger than the corresponding element of the `lo` axis,
        even if the `lo` array is in descending order.

    """

    _lo = None
    _hi = None

    @overload
    def __init__(self,
                 lo: ArrayType,
                 hi: ArrayType
                 ) -> None:
        ...

    @overload
    def __init__(self,
                 lo: None,
                 hi: None
                 ) -> None:
        ...

    def __init__(self,
                 lo: ArrayType | None,
                 hi: ArrayType | None
                 ) -> None:
        self._lo = _to_readable_array(lo)
        self._hi = _to_readable_array(hi)

        if self._lo is None:
            if self._hi is None:
                return

            raise DataErr("mismatchn", "lo", "hi", "None", len(self._hi))

        nlo = len(self._lo)
        if self._hi is None:
            raise DataErr("mismatchn", "lo", "hi", nlo, "None")

        nhi = len(self._hi)
        if nlo != nhi:
            raise DataErr("mismatchn", "lo", "hi", nlo, nhi)

        if nlo > 0:
            self._is_ascending = self._lo[-1] > self._lo[0]

    @property
    def lo(self) -> np.ndarray | None:
        """Read-only access to the axis data (low edge)."""
        return self._lo

    @property
    def hi(self) -> np.ndarray | None:
        """Read-only access to the axis data (high edge)."""
        return self._hi

    @property
    def is_integrated(self) -> bool:
        return True

    # TODO: what happens if the axis is empty?
    @property
    def start(self) -> float:
        if self.is_ascending:
            return self.lo[0]

        return self.lo[-1]

    # TODO: what happens if the axis is empty?
    @property
    def end(self) -> float:
        if self.is_ascending:
            return self.hi[-1]

        return self.hi[0]

    @property
    def size(self) -> int:
        if self.lo is None:
            return 0

        return self.lo.size


# TODO: how best to send in metadata, such as axis information (e.g.
# grouping data, units, or WCS information)? Does it make sense to
# include any or all of this? How about a generic "metadata" field?
#
# This could include the dependent axis as well, but there are many
# times this would not be useful, so leave out for now.
#
@dataclass(kw_only=True)
class Axes:
    """Represent a set of axes."""

    axes: list[Axis]
    """The axis or axes."""
