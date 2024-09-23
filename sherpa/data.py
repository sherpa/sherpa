#
#  Copyright (C) 2008, 2015 - 2017, 2019 - 2024
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

"""Tools for creating, storing, inspecting, and manipulating data sets.

The main classes for representing data sets are `Data1D`,
`Data1DInt`, and `Data2D`, to handle (x, y), (xlo,
xhi, y), and (x1, x2, y) data, although there are also
more-specialized cases, such as `Data1DAsymmetricErrs`. These
classes build on the `Data` class, which supports dynamic
filtering of data - to select a subset of the data range - as well as
data access and model evaluation to match the data range.

The `Filter` class is used to handle data filtering - that
is, to combine filters such as selecting the range a to b (``notice``)
and hiding the range c to d (``ignore``). This is used with the
`DataSpace1D` and `DataSpace2D` classes to handle
evaluating models on different grids to the data, and then converting
back to the data space, whether by rebinning or interpolation.

Design
------

The design for the `Data` class assumes

- the basic fields are the independent axis - which is thought of as a
  tuple of arrays - and the dependent axis, where normally the
  independent axis is labelled starting with ``x``, and the dependent axis
  ``y``.

- there are a number of related fields - such as the statistical and
  systematic errors - that follow the same behavior as the
  independent and dependent axes.

- fields are converted to `ndarray` when read in.

- the independent axes can either be points (`~sherpa.models.regrid.PointAxis`)
  or integrated (`~sherpa.models.regrid.IntegratedAxis`). There is
  currently no support for an object with a combination of point and
  integrated axes but it could be added.

- the independent axis is mediated by a "Data Space" (`DataSpaceND`,
  `DataSpace1D`, `IntegratedDataSpace1D`, `DataSpace2D`, and
  `IntegratedDataSpace2D`).

- `Data` objects contain the `~Data.ndim` field to indicate the number of
  independent axes it supports.

- to support arbitrary dimensions, the data is treated as one-dimensional
  arrays.

- `Data` objects can be created with no data, but are fixed once an
  axis - expected to be the independent axis but it need not be - is
  set (see the `~Data.size` attribute). In general the data is
  assumed to be set when the object is created.

- there are checks to ensure that the data has the correct size, shape,
  and potentially data type, but this is not intended to catch all
  possible problems.

- Numpy masked arrays can be used to initialize the dependent variable
  and the mask is converted to the format of the `~Data.mask`
  attribute of the `Data` object, taking into account that for Sherpa
  a value of True indicates a valid quantity, while the opposite is
  true in numpy.

  In general, it is expected that any filtering has been applied
  before the `Data` object is created.

- the independent axes are marked as read-only, so the only way to change
  a value is to replace the whole axis, in which case any existing
  filter will be removed.

Notebook support
----------------

The Data objects support the rich display protocol of IPython, with
HTML display of a table of information highlighting the relevant data.
Examples can be found at [NoteBook]_.

References
----------

.. [NoteBook] https://sherpa.readthedocs.io/en/latest/NotebookSupport.html

Examples
--------

Create a data set representing the independent axis (``x``) and
dependent axis (``y``) then filter to select only those values between
500-520 and 530-700:

>>> import numpy as np
>>> x = np.arange(1000)
>>> y = np.random.normal(size=1000)
>>> d1 = Data1D('example', x, y)
>>> d1.notice(500, 700)
>>> d1.ignore(520, 530)

"""

from abc import ABCMeta
import logging
from typing import Any, Literal, Optional, Sequence, Union, \
    overload
import warnings

import numpy as np

from sherpa.models.regrid import EvaluationSpace1D, IntegratedAxis, PointAxis
from sherpa.utils import NoNewAttributesAfterInit, formatting, \
    print_fields, create_expr, create_expr_integrated, \
    calc_total_error, bool_cast, filter_bins
from sherpa.utils.err import DataErr
from sherpa.utils.numeric_types import SherpaFloat
from sherpa.utils.parallel import parallel_map_funcs
from sherpa.utils.types import ArrayType, ModelFunc, StatErrFunc

warning = logging.getLogger(__name__).warning


__all__ = ('Data', 'DataSimulFit', 'Data1D', 'Data1DInt',
           'Data1DAsymmetricErrs', 'Data2D', 'Data2DInt')


# The alias does not really save any characters, but it's used
# throughout the code so it seems useful.
#
FieldsType = tuple[str, ...]


@overload
def _check(array: None) -> None:
    ...

@overload
def _check(array: ArrayType) -> np.ndarray:
    ...

def _check(array):
    """Ensure the data is a 1D array, or can be converted to one."""

    # Special case support for None.
    if array is None:
        return None

    # Assume that having a shape attribute means "this is near-enough
    # an ndarray" that we don't need to convert it.
    #
    if hasattr(array, "shape"):
        # Special-case the 0-dim case
        if len(array.shape) == 0:
            raise DataErr("notanarray")
        if len(array.shape) != 1:
            raise DataErr("not1darray")

        return array

    return _check(np.asanyarray(array))


# Can we get away with assuming array is not None?
#
@overload
def _check_nomask(array: None) -> None:
    ...

@overload
def _check_nomask(array: ArrayType) -> np.ndarray:
    ...

def _check_nomask(array):
    if hasattr(array, 'mask'):
        warnings.warn(f'Input array {array} has a mask attribute. Because masks are supported for dependent variables only the mask attribute of the independent array is ignored and values `behind the mask` are used.')

    # Ensure we do NumPy conversions after checking for mask to
    # make sure we don't end up removing .mask in the following.
    #
    return _check(array)


@overload
def _check_dep(array: None) -> tuple[None, Literal[True]]:
    ...

@overload
def _check_dep(array: ArrayType) -> tuple[np.ndarray, bool]:
    ...

def _check_dep(array):
    if not hasattr(array, 'mask'):
        return _check(array), True

    # We know the mask convention is opposite to sherpa
    if isinstance(array, np.ma.MaskedArray):
        return _check(array), ~array.mask

    # We don't know what the mask convention is
    warnings.warn(f'Format of mask for array {array} not supported thus the mask is is ignored and values `behind the mask` are used. Set .mask attribute manually or use "set_filter" function.')
    return _check(array), True


class DataSpace1D(EvaluationSpace1D):
    """
    Class for representing 1-D Data Space. Data Spaces are spaces that describe the data domain. As models can be
    evaluated over data spaces, data spaces can be considered evaluation spaces themselves. However this "is-a"
    relationship is in the code mostly for convenience and could be removed in future versions.
    """
    def __init__(self, filter, x):
        """
        Parameters
        ----------
        filter : Filter
            a filter object that initialized this data space
        x : array_like
            the x axis of this data space
        """
        self.filter = filter
        super().__init__(_check_nomask(x))

    def get(self, filter=False):
        """
        Get a filtered representation of this data set. If `filter` is `False` this object is returned.

        Parameters
        ----------
        filter : bool
            whether the data set should be filtered before being returned

        Returns
        -------
        DataSpace1D
        """

        if not bool_cast(filter):
            return self

        data = self.grid[0]

        data = self.filter.apply(data)
        return DataSpace1D(self.filter, data)

    def for_model(self, model):
        """
        Models can be defined over arbitrary evaluation spaces. However,
        at evaluation time during a fit, the model's evaluation space shall
        be done at the user's request space only and set to 0 every where else.

        Parameters
        ----------
        model : The model whose evaluation space needs to be joined with the dataset's data space.

        Returns
        -------
        DataSpace1D
            A data space that joins this data space with the model's evaluation space. if the model does not have an
            evaluation space assigned to itself then `self` is returned.
        """
        evaluation_space = None

        if model is not None and hasattr(model, "evaluation_space") \
           and self not in model.evaluation_space:
            evaluation_space = self

        return self if evaluation_space is None else evaluation_space


class IntegratedDataSpace1D(EvaluationSpace1D):
    """
    Same as DataSpace1D, but for supporting integrated data sets.
    """
    def __init__(self, filter, xlo, xhi):
        """
        Parameters
        ----------
        filter : Filter
            a filter object that initialized this data space
        xlo : array_like
            the lower bounds array of this data space
        xhi : array_like
            the higher bounds array of this data space
        """
        xlo = _check_nomask(xlo)
        xhi = _check_nomask(xhi)

        # We only have to check if xhi is None, as if it's set then the
        # superclass will check the length. We do need to support both
        # values being None though.
        #
        if xlo is not None and xhi is None:
            raise DataErr("mismatchn", "lo", "hi", len(xlo), "None")

        self.filter = filter
        super().__init__(xlo, xhi)

    def get(self, filter=False):
        """
        Get a filtered representation of this data set. If `filter` is `False` this object is returned.

        Parameters
        ----------
        filter : bool
            whether the data set should be filtered before being returned

        Returns
        -------
        IntegratedDataSpace1D
        """

        if not bool_cast(filter):
            return self

        data = self.grid

        data = tuple(self.filter.apply(axis) for axis in data)
        return IntegratedDataSpace1D(self.filter, *data)

    def for_model(self, model):
        """
        Models can be defined over arbitrary evaluation spaces. However, at evaluation time during a fit, the model's
        evaluation space and the data space will be joined together and the model will be evaluated over the joined
        domain. This makes sure that when the models are rebinned back to the data space the evaluation does not have
        to be extrapolated from the model's evaluation space alone.

        Parameters
        ----------
        model : The model whose evaluation space needs to be joined with the dataset's data space.

        Returns
        -------
        IntegratedDataSpace1D
            A data space that joins this data space with the model's evaluation space. if the model does not have an
            evaluation space assigned to itself then `self` is returned.
        """
        evaluation_space = None

        if model is not None and hasattr(model, "evaluation_space") \
           and self not in model.evaluation_space:
            evaluation_space = self

        return self if evaluation_space is None else evaluation_space


# We do not inherit from EvaluationSpace2D because that would require
# that the x0 and x1 arrays form a grid (e.g. nx by ny) and we need
# to be able to support a non-grid set of points.
#
# The code is written to resemble sherpa.models.grid.EvaluationSpace2D
# though, so that the checks added for that code are run here.
#
class DataSpace2D:
    """Class for representing 2-D Data Spaces.

    Data Spaces are spaces that describe the data domain.

    Parameters
    ----------
    filter : Filter
        A filter object that initialized this data space.
    x0, x1 : array_like
        The first and second axes of this data space. The arrays are
        copied.

    """

    def __init__(self, filter, x0, x1):
        x0 = _check_nomask(x0)
        x1 = _check_nomask(x1)

        self.filter = filter
        self.x_axis = PointAxis(x0)
        self.y_axis = PointAxis(x1)
        if self.x_axis.size != self.y_axis.size:
            raise DataErr("mismatchn", "x0", "x1", self.x_axis.size, self.y_axis.size)

    def get(self, filter=False):
        """
        Get a filtered representation of this data set. If `filter` is `False` this object is returned.

        Parameters
        ----------
        filter : bool
            whether the data set should be filtered before being returned

        Returns
        -------
        DataSpace2D
        """

        if not bool_cast(filter):
            return self

        data = self.grid

        data = tuple(self.filter.apply(axis) for axis in data)
        return DataSpace2D(self.filter, *data)

    @property
    def grid(self):
        """The grid representation of this dataset.

        The x0 and x1 arrays in the grid are one-dimensional representations of the meshgrid obtained
        from the x and y axis arrays, as in `numpy.meshgrid(x, y)[0].ravel()`

        Returns
        -------
        (x0, x1) : (ndarray, ndarray)
            A tuple representing the x0 and x1 axes.
        """
        return self.x_axis.x, self.y_axis.x

    # Should these be deprecated now that the code attempts to mimic EvaluationSpace2D?
    #
    @property
    def x0(self):
        """Return the first axis."""
        return self.x_axis.x

    @property
    def x1(self):
        """Return the second axis."""
        return self.y_axis.x


class IntegratedDataSpace2D:
    """Same as DataSpace2D, but for supporting integrated data sets.

    Parameters
    ----------
    filter : Filter
        A filter object that initialized this data space.
    x0lo, x1lo, x0hi, x1hi : array_like
        The lower bounds array of the first and second axes, then the
        upper bounds.

    """

    def __init__(self, filter, x0lo, x1lo, x0hi, x1hi):
        x0lo = _check_nomask(x0lo)
        x1lo = _check_nomask(x1lo)
        x0hi = _check_nomask(x0hi)
        x1hi = _check_nomask(x1hi)

        self.filter = filter
        self.x_axis = IntegratedAxis(x0lo, x0hi)
        self.y_axis = IntegratedAxis(x1lo, x1hi)
        if self.x_axis.size != self.y_axis.size:
            raise DataErr("mismatchn", "x0", "x1", self.x_axis.size, self.y_axis.size)

    def get(self, filter=False):
        """
        Get a filtered representation of this data set. If `filter` is `False` this object is returned.

        Parameters
        ----------
        filter : bool
            whether the data set should be filtered before being returned

        Returns
        -------
        IntegratedDataSpace2D
        """

        if not bool_cast(filter):
            return self

        data = self.grid

        data = tuple(self.filter.apply(axis) for axis in data)
        return IntegratedDataSpace2D(self.filter, *data)

    @property
    def grid(self):
        """The grid representation of this dataset.

        The x0 and x1 arrays in the grid are one-dimensional representations of the meshgrid obtained
        from the x and y axis arrays, as in `numpy.meshgrid(x, y)[0].ravel()`

        Returns
        -------
        (x0lo, x1lo, x0hi, x1hi) : (ndarray, ndarray, ndarray, ndarray)
            The axis data.

        """
        return self.x_axis.lo, self.y_axis.lo, self.x_axis.hi, self.y_axis.hi

    # Should these be deprecated now that the code attempts to mimic EvaluationSpace2D?
    #
    @property
    def x0lo(self):
        """Return the first axis (low edge)"""
        return self.x_axis.lo

    @property
    def x0hi(self):
        """Return the first axis (high edge)."""
        return self.x_axis.hi

    @property
    def x1lo(self):
        """Return the second axis (low edge)."""
        return self.y_axis.lo

    @property
    def x1hi(self):
        """Return the second axis (high edge)."""
        return self.y_axis.hi


class DataSpaceND:
    """Class for representing arbitrary N-Dimensional data domains

    Parameters
    ----------
    filter : Filter
        A filter object that initialized this data space.
    indep : tuple of array_like
        The tuple of independent axes. The arrays are copied rather
        than storing a reference.

    """

    def __init__(self, filter, indep):
        self.filter = filter
        self.indep = tuple(_check_nomask(d) for d in indep)

    def get(self, filter=False):
        """
        Get a filtered representation of this data set. If `filter` is `False` this object is returned.

        Parameters
        ----------
        filter : bool
            whether the data set should be filtered before being returned

        Returns
        -------
        DataSpaceND
        """

        if not bool_cast(filter):
            return self

        data = tuple(self.filter.apply(axis) for axis in self.indep)
        return DataSpaceND(self.filter, data)

    @property
    def grid(self):
        """
        Return the grid representation of this dataset.

        The independent arrays are returned unchanged, i.e. unlike the DataSpace2D class they are not meshed

        Returns
        -------
        tuple
            A tuple representing the independent axes.
        """
        return self.indep


# NOTE:
#
# Although this is labelled as supporting N-D datasets, the
# implementation assumes that the data has been flattened to 1D arrays
# - in particular the notice method. It is likely that we can document
# this - i.e. that the mask is going to be 1D.
#
class Filter:
    """A class for representing filters of N-Dimensional datasets.

    The filter does not know the size of the dataset or the values of
    the independent axes.

    """
    def __init__(self) -> None:
        self._mask: Union[np.ndarray, bool] = True

    @property
    def mask(self) -> Union[np.ndarray, bool]:
        """Mask array for dependent variable

        Returns
        -------
        mask : bool or numpy.ndarray
        """
        return self._mask

    @mask.setter
    def mask(self, val: Union[ArrayType, bool]) -> None:
        if val is None:
            raise DataErr('ismask')

        # The code below has to deal with bool-like values, such as
        # numpy.ma.nomask (this evaluates to False but is not a bool),
        # and numpy.bool_ values. The following code is intended to
        # catch these "special" values. Note that we explicitly check
        # for boolean values ('val is True' and 'val is False') rather
        # than just whether val behaves like a boolean - for example
        # 'if val or not val: ...'.
        #
        if not np.isscalar(val):
            # Is it possible for the following to throw a conversion
            # error? If so, we could catch it and convert it into a
            # DataErr, but it does not seem worth it (as it's not
            # obvious what the error to catch would be).
            #
            self._mask = np.asarray(val, bool)

        elif val is np.ma.nomask:
            self._mask = True

        elif (val is True) or (val is False):
            self._mask = val

        elif isinstance(val, np.bool_):
            # are there other types we could be looking for here?
            self._mask = bool(val)

        else:
            # Note that setting mask to 2.0 will fail, but an array
            # of 2.0's will get converted to booleans.
            #
            raise DataErr('ismask')

    @overload
    def apply(self, array: None) -> None:
        ...

    @overload
    def apply(self, array: ArrayType) -> np.ndarray:
        ...

    def apply(self, array):
        """Apply this filter to an array

        Parameters
        ----------
        array : array_like or None
            Array to be filtered

        Returns
        -------
        array_like : ndarray or None

        Raises
        ------
        sherpa.utils.err.DataErr
            The filter has removed all elements or there is a
            mismatch between the `mask` and the ``array`` argument.

        See Also
        --------
        notice

        """
        if array is None:
            return None

        # Note that mask may not be a boolean but an array.
        if self.mask is False:
            raise DataErr('notmask')

        # Ensure we always return a ndarray.
        #
        array = _check(array)
        if self.mask is True:
            return array

        if array.size != self.mask.size:
            raise DataErr("mismatchn", "mask", "data array", self.mask.size, array.size)

        return array[self.mask]

    def notice(self,
               mins: ArrayType,
               maxes: ArrayType,
               axislist: Sequence[ArrayType],
               ignore: bool = False,
               integrated: bool = False
               ) -> None:
        """Select a range to notice or ignore (remove).

        The ``axislist`` argument is expected to be sent the
        independent axis of a `Data` object - so ``(x, )`` for
        one-dimensional data, ``(xlo, xhi)`` for integrated
        one-dimensional data, ``(x0, x1)`` for two-dimensional data,
        and ``(x0lo, x1lo, x0hi, x1hi)`` for integrated two-dimensinal
        data. The ``mins`` and ``maxes`` must then be set to match
        this ordering.

        Parameters
        ----------
        mins : sequence of values
           The minimum value of the valid range (elements may be None
           to indicate no lower bound). When not None, it is treated
           as an inclusive limit, so points >= min are included.
        maxes : sequence of values
           The maximum value of the valid range (elements may be None
           to indicate no upper bound). It is treated as an inclusive
           limit (points <= max) when integrated is False, and an
           exclusive limit (points < max) when integrated is True.
        axislist: sequence of arrays
           The axis to apply the range to. There must be the same
           number of elements in mins, maxes, and axislist.  The
           number of elements of each element of axislist must also
           agree (the cell values do not need to match).
        ignore : bool, optional
           If True the range is to be ignored, otherwise it is
           included.  The default is to include the range.
        integrated : bool, optional
           Is the data integrated (we have low and high bin edges)?  The
           default is False. When True it is expected that axislist
           contains a even number of rows, where the odd values are the
           low edges and the even values the upper edges, and that the
           mins and maxes only ever contain a single value, given in
           (None, hi) and (lo, None) ordering.

        See Also
        --------
        apply

        Examples
        --------

        Select points in xs which are in the range 1.5 <= x <= 4:

        >>> f = Filter()
        >>> f.mask
        True
        >>> xs = [1, 2, 3, 4, 5]
        >>> f.notice([1.5], [4], (xs, ))
        >>> f.mask
        array([False,  True,  True,  True, False])

        Filter the data to select all points with x0 >= 1.5 and x1 <= 4:

        >>> f = Filter()
        >>> x0 = [1, 1.4, 1.6, 2, 3]
        >>> x1 = [2, 2, 4, 4, 6]
        >>> f.notice([1.5, None], [None, 4], (x0, x1))
        >>> f.mask
        array([False, False,  True,  True, False])

        For integrated data sets the lower and upper edges should be
        sent separately with the max and min limits, along with
        setting the integrated flag. The following selects the bins
        that cover the range 2 to 4 and 1.5 to 3.5:

        >>> xlo = [1, 2, 3, 4, 5]
        >>> xhi = [2, 3, 4, 5, 6]
        >>> f = Filter()
        >>> f.notice([None, 2], [4, None], (xlo, xhi), integrated=True)
        >>> f.mask
        array([False,  True,  True,  False, False])
        >>> f.notice([None, 1.5], [3.5, None], (xlo, xhi), integrated=True)
        >>> f.mask
        array([ True,  True,  True, False, False])

        """

        # If integrated is True then we should have an even number
        # of axislist elements, but we do not require this.
        #
        ignore = bool_cast(ignore)
        for vals, label in zip([mins, maxes, axislist],
                               ['lower bound', 'upper bound', 'grid']):
            if any(isinstance(val, str) for val in vals):
                raise DataErr('typecheck', label)

        mask = filter_bins(mins, maxes, axislist, integrated=integrated)

        if mask is None:
            self.mask = not ignore
        elif not ignore:
            if self.mask is True:
                self.mask = mask
            else:
                self.mask |= mask
        else:
            mask = ~mask
            if self.mask is False:
                self.mask = mask
            else:
                self.mask &= mask


class BaseData(metaclass=ABCMeta):
    """
    Base class for all data classes. Left for compatibility with older versions.
    """
    pass


# DATA-NOTE: ND Data cannot be plotted
class Data(NoNewAttributesAfterInit, BaseData):
    """Generic, N-Dimensional data sets.

    A data class is the collection of a data space and a number of
    data arrays for the dependent variable and associated errors.

    Parameters
    ----------
    name : string
        name of this dataset
    indep: tuple of array_like
        the tuple of independent arrays.
    y : array_like
        The values of the dependent observable. If this is a numpy
        masked array, the mask will used to initialize a mask.
    staterror : array_like
        the statistical error associated with the data
    syserror : array_like
        the systematic error associated with the data

    Notes
    -----

    This class can be extended by classes defining data sets of
    specific dimensionality. Extending classes should override the
    `_init_data_space` method.

    This class provides most of the infrastructure for extending
    classes for free.

    Data classes contain a ``mask`` attribute, which can be used
    ignore certain values in the array when fitting or plotting that
    data. The convention in Sherpa is that ``True`` marks a values as
    *valid* and ``False`` as *invalid* (note that this is opposite to
    the numpy convention). When a `Data` instance is initialized with
    a dependent array that has a ``mask`` attribute (e.g. numpy masked
    array), it will attempt to convert that mask to the Sherpa
    convention and raise a warning otherwise. In any case, the user
    can set ``data.mask`` after initialization if that conversion does
    not yield the expected result.

    """

    _fields: FieldsType = ("name", "indep", "dep", "staterror", "syserror")
    """The main data values stored by the object (as a tuple).

    This is used to identify the column data - that is values that
    can be written out in tabular form - other than the "name"
    field. Other fields are listed in _extra_fields.
    """

    _extra_fields: FieldsType = ()
    """Any extra fields that should be displayed by str(object)."""

    _related_fields: FieldsType = ("y", "staterror", "syserror")
    """What fields must match the size of the independent axis.

    These fields are set to None whenever the independent axis size
    is set or changed.
    """

    _y: Optional[np.ndarray] = None
    _size: Optional[int] = None
    _staterror: Optional[np.ndarray] = None
    _syserror: Optional[np.ndarray] = None

    ndim: Optional[int] = None
    "The dimensionality of the dataset, if defined, or None."

    def __init__(self,
                 name: str,
                 indep: Union[Sequence[ArrayType],
                              Sequence[None]],
                 y: Optional[ArrayType],
                 staterror: Optional[ArrayType] = None,
                 syserror: Optional[ArrayType] = None
                 ) -> None:
        self.name = name
        self._data_space = self._init_data_space(Filter(), *indep)
        self.y, self.mask = _check_dep(y)
        self.staterror = staterror
        self.syserror = syserror

        self._ylabel = 'y'
        NoNewAttributesAfterInit.__init__(self)

    def _check_data_space(self, dataspace):
        """Check that the data space has the correct size.

        Note that this also sets the size of the data object if it has
        not been set.

        Parameters
        ----------
        dataspace
            The dataspace object for the data object.

        """

        indep0 = dataspace.get().grid[0]
        nold = self.size

        # If there is no data then we need to update the size, unless
        # the new data is also empty.
        #
        if nold is None:
            if indep0 is None:
                return

            self._size = len(indep0)
            return

        # We could allow the independent axis to be reset to None.
        #
        if indep0 is None:
            raise DataErr("independent axis can not be cleared")

        nnew = len(indep0)
        if nnew != nold:
            raise DataErr(f"independent axis can not change size: {nold} to {nnew}")

    def _init_data_space(self,
                         filter: Filter,
                         *data: ArrayType
                         ) -> DataSpaceND:
        """
        Extending classes should implement this method to provide the proper data space construction.

        Parameters
        ----------
        filter : Filter
            a filter object passed by the initializer upon initialization of extending classes.
        data : tuple of array_like
            the tuple of independent arrays used to build the data space.

        Returns
        -------
        object
            an instance of the dataspace associated with this data set.

        """

        ds = DataSpaceND(filter, data)
        self._check_data_space(ds)
        return ds

    def _set_related(self,
                     attr: str,
                     val: Optional[ArrayType],
                     check_mask: bool = True,
                     **kwargs) -> None:
        """Set a field that must match the independent axes size.

        The value can be None or something with the same length as the
        independent axis. This is intended to be used from the property
        setter. There is also a check to warn the user if the value
        contains a NumPy masked array which does not match the
        dependent axis.

        """
        if val is None:
            setattr(self, f"_{attr}", None)
            return

        # This will check whether val is an array or not.
        #
        vals = _check(val)
        nval = len(val)

        # Check the mask value of val, and not vals, since _check
        # could have called np.asarray and so lost the mask setting.
        #
        if check_mask and hasattr(val, "mask"):
            if not hasattr(self.y, "mask") or \
               len(self.y) != nval or \
                   not np.all(self.y.mask == val.mask):

                warnings.warn(f"The mask of {attr} differs from the dependent array, only the mask of the dependent array is used in Sherpa.")

        nelem = self.size
        if nelem is None:
            # We set the object size here
            self._size  = nval
            setattr(self, f"_{attr}", vals)
            return

        if nval != nelem:
            raise DataErr('mismatchn', 'independent axis', attr, nelem, nval)

        setattr(self, f"_{attr}", vals)

    @property
    def y(self) -> Optional[np.ndarray]:
        """The dependent axis.

        If set, it must match the size of the independent axes.
        """
        return self._y

    @y.setter
    def y(self, val: ArrayType) -> None:
        self._set_related("y", val, check_mask=False)

    @property
    def size(self) -> Optional[int]:
        """The number of elements in the data set.

        Returns
        -------
        size : int or None
            If the size has not been set then None is returned.
        """
        return self._size

    # Allow users to len() a data object. The idea is that this
    # represents the "number of elements" (e.g. the size of the
    # independent axis). Does this make sense for some of the data
    # classes (in particular RMF)?
    #
    def __len__(self) -> int:
        if self.size is None:
            return 0

        return self.size

    @property
    def dep(self) -> Optional[np.ndarray]:
        """
        Left for compatibility with older versions
        """
        return self.y

    @dep.setter
    def dep(self, val: Optional[ArrayType]) -> None:
        self.y = val

    @property
    def mask(self) -> Union[np.ndarray, bool]:
        """
        Mask array for dependent variable

        Returns
        -------
        mask : bool or numpy.ndarray
        """
        return self._data_space.filter.mask

    @mask.setter
    def mask(self, val: Union[ArrayType, bool]) -> None:

        # If we have a scalar then
        # - we do not check sizes (as it's possible to set even though
        #   the independent axis is not set)
        # - we do not want to convert it to a ndarray
        #
        # val can not be None but that is caught when setting the
        # filter.mask attribute.
        #
        if val is not None and np.ndim(val) != 0:
            if self.size is None:
                raise DataErr("The independent axis has not been set yet")

            nval = len(val)
            if nval != self.size:
                raise DataErr("mismatchn", "independent axis", "mask",
                              self.size, nval)

        self._data_space.filter.mask = val

    # This is overloaded by Data1D and Data2D so can it be made "virtual"?
    #
    def get_dims(self) -> tuple[int, ...]:
        """
        Return the dimensions of this data space as a tuple of tuples.
        The first element in the tuple is a tuple with the dimensions of the data space, while the second element
        provides the size of the dependent array.

        Returns
        -------
        tuple
        """
        # What should this do when the independent or dependent axes
        # are not set?
        indep_size = tuple(indep.size for indep in self.indep)
        return indep_size, self.dep.size

    def _clear_filter(self) -> None:
        """Clear out the existing filter.

        This is designed for use by @indep.setter.

        """

        # This is currently a no-op. It may be overridden by a
        # subclass.
        #
        pass

    @property
    def indep(self
              ) -> Union[tuple[np.ndarray, ...],
                         tuple[None, ...]]:
        """The grid of the data space associated with this data set.

        When set, the field must be set to a tuple, even for a
        one-dimensional data set. The "related" fields such as the
        dependent axis and the error fields are set to None if their
        size does not match.

        .. versionchanged:: 4.14.1
           The filter created by `notice` and `ignore` is now cleared
           when the independent axis is changed.

        Returns
        -------
        tuple of array_like or None

        """
        return self._data_space.get().grid

    @indep.setter
    def indep(self,
              val: Union[tuple[ArrayType, ...],
                         tuple[None, ...]]
              ) -> None:

        # This is a low-level check so raise a normal Python error
        # rather than DataErr. Do we want to allow a sequence here,
        # that is not forced to be a tuple? The concern then is that it gets
        # harder to distinguish between a user accidentally giving x,
        # where it's an array, rather than (x,), say.
        #
        if not isinstance(val, tuple):
            raise TypeError(f"independent axis must be sent a tuple, not {type(val).__name__}")

        ds = self._init_data_space(self._data_space.filter, *val)
        self._data_space = ds
        self._clear_filter()

    def get_indep(self,
                  filter: bool = False
                  ) -> Union[tuple[np.ndarray, ...],
                             tuple[None, ...]]:
        """Return the independent axes of a data set.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.

        Returns
        -------
        axis: tuple of arrays
           The independent axis values for the data set. This gives
           the coordinates of each point in the data set.

        See Also
        --------
        get_dep : Return the dependent axis of a data set.

        """
        return self._data_space.get(filter).grid

    def set_indep(self,
                  val: Union[tuple[ArrayType, ...],
                             tuple[None, ...]]
                  ) -> None:
        self.indep = val

    def get_dep(self,
                filter: bool = False
                ) -> Optional[np.ndarray]:
        """Return the dependent axis of a data set.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.

        Returns
        -------
        axis: array
           The dependent axis values for the data set. This gives
           the value of each point in the data set.

        See Also
        --------
        get_indep : Return the independent axis of a data set.
        get_error : Return the errors on the dependent axis of a data set.
        get_staterror : Return the statistical errors on the dependent axis of a data set.
        get_syserror : Return the systematic errors on the dependent axis of a data set.

        """
        dep = self.dep
        if bool_cast(filter):
            dep = self.apply_filter(dep)

        return dep

    # Do we want to allow val to be None?
    def set_dep(self, val: Union[ArrayType, float]) -> None:
        """
        Set the dependent variable values.

        Parameters
        ----------
        val : sequence or number
            If a number then it is used for each element.

        """
        if np.iterable(val):
            dep = np.asarray(val, SherpaFloat)
        else:
            nelem = self.size
            if nelem is None:
                raise DataErr("sizenotset", self.name)

            dep = SherpaFloat(val) * np.ones(nelem, dtype=SherpaFloat)

        self.y = dep

    # It is not clear when default values need to be given in overload
    # statements. The current choice is based on mypy 1.10.0.
    #
    @overload
    def get_y(self,
              filter: bool,
              yfunc: None,
              use_evaluation_space: bool = False
              ) -> np.ndarray:
        ...

    @overload
    def get_y(self,
              filter: bool,
              yfunc: ModelFunc,
              use_evaluation_space: bool = False
              ) -> tuple[np.ndarray, ArrayType]:
        ...

    def get_y(self, filter=False, yfunc=None, use_evaluation_space=False):
        """Return dependent axis in N-D view of dependent variable

        Parameters
        ----------
        filter
        yfunc
        use_evaluation_space

        Returns
        -------
        y: array or (array, array) or None
            If yfunc is not None and the dependent axis is set then
            the return value is (y, y2) where y2 is yfunc evaluated on
            the independent axis.

        """
        y = self.get_dep(filter)
        if yfunc is None:
            return y

        if filter:
            y2 = self.eval_model_to_fit(yfunc)
        else:
            y2 = self.eval_model(yfunc)

        return (y, y2)

    @property
    def staterror(self) -> Optional[np.ndarray]:
        """The statistical error on the dependent axis, if set.

        This must match the size of the independent axis.
        """
        return self._staterror

    @staterror.setter
    def staterror(self, val: ArrayType) -> None:
        self._set_related("staterror", val)

    @property
    def syserror(self) -> Optional[np.ndarray]:
        """The systematic error on the dependent axis, if set.

        This must match the size of the independent axis.
        """
        return self._syserror

    @syserror.setter
    def syserror(self, val: ArrayType) -> None:
        self._set_related("syserror", val)

    def get_staterror(self,
                      filter: bool = False,
                      staterrfunc: Optional[StatErrFunc] = None
                      ) -> Optional[ArrayType]:
        """Return the statistical error on the dependent axis of a data set.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.
        staterrfunc : function
           If no statistical error has been set, the errors will
           be calculated by applying this function to the
           dependent axis of the data set.

        Returns
        -------
        axis : array or `None`
           The statistical error for each data point. A value of
           `None` is returned if the data set has no statistical error
           array and `staterrfunc` is `None`.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        get_syserror : Return the systematic errors on the dependent axis of a data set.

        """
        staterror = getattr(self, 'staterror', None)
        if bool_cast(filter):
            staterror = self.apply_filter(staterror)

        if (staterror is None) and (staterrfunc is not None):
            dep = self.get_dep(filter)
            staterror = staterrfunc(dep)

        return staterror

    def get_syserror(self,
                     filter: bool = False
                     ) -> Optional[np.ndarray]:
        """Return the systematic error on the dependent axis of a data set.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.

        Returns
        -------
        axis : array or None
           The systematic error for each data point. A value of
           `None` is returned if the data set has no systematic
           errors.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        get_staterror : Return the statistical errors on the dependent axis of a data set.

        """
        syserr = getattr(self, 'syserror', None)
        if bool_cast(filter):
            syserr = self.apply_filter(syserr)

        return syserr

    def get_error(self, filter=False, staterrfunc=None):
        """Return the total error on the dependent variable.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.
        staterrfunc : function
           If no statistical error has been set, the errors will
           be calculated by applying this function to the
           dependent axis of the data set.

        Returns
        -------
        axis : array or `None`
           The error for each data point, formed by adding the
           statistical and systematic errors in quadrature.

        See Also
        --------
        get_dep : Return the independent axis of a data set.
        get_staterror : Return the statistical errors on the dependent axis of a data set.
        get_syserror : Return the systematic errors on the dependent axis of a data set.

        """
        return calc_total_error(self.get_staterror(filter, staterrfunc),
                                self.get_syserror(filter))

    def get_yerr(self, filter=False, staterrfunc=None):
        """Return errors in dependent axis in N-D view of dependent variable.

        Parameters
        ----------
        filter
        staterrfunc

        Returns
        -------

        """
        return self.get_error(filter, staterrfunc)

    # It looks like we never set yfunc when calling this method so we
    # could remove it.
    #
    def get_ylabel(self, yfunc=None) -> str:
        """Return label for dependent axis in N-D view of dependent variable.

        Parameters
        ----------
        yfunc
           Unused.

        Returns
        -------
        label: str
           The label.

        See Also
        --------
        set_ylabel

        """
        return self._ylabel

    def set_ylabel(self, label: str) -> None:
        """Set the label for the dependent axis.

        .. versionadded:: 4.17.0

        Parameters
        ----------
        label: str
            The new label.

        See Also
        --------
        get_ylabel

        """
        self._ylabel = label

    @overload
    def apply_filter(self, data: None) -> None:
        ...

    @overload
    def apply_filter(self, data: ArrayType) -> np.ndarray:
        ...

    def apply_filter(self, data):
        if data is None:
            return None

        nelem = self.size
        if nelem is None:
            raise DataErr("sizenotset", self.name)

        data = _check(data)
        ndata = len(data)

        if ndata != nelem:
            raise DataErr("mismatchn", "data", "array", nelem, ndata)

        return self._data_space.filter.apply(data)

    def notice(self, mins, maxes,
               ignore: bool = False,
               integrated: bool = False
               ) -> None:
        self._data_space.filter.notice(mins, maxes, self.get_indep(),
                                       ignore=ignore,
                                       integrated=integrated)

    def ignore(self, *args, **kwargs) -> None:
        kwargs['ignore'] = True
        self.notice(*args, **kwargs)

    def _can_apply_model(self, modelfunc: ModelFunc) -> None:
        """Check if model dimensions match data dimensions."""

        # Should this also check whether the independent axis is set?
        mdim = getattr(modelfunc, "ndim", None)
        ddim = getattr(self, "ndim", None)
        if None not in [mdim, ddim] and mdim != ddim:
            raise DataErr(f"Data and model dimensionality do not match: {ddim}D and {mdim}D")

    def eval_model(self, modelfunc: ModelFunc) -> ArrayType:
        """Evaluate the model on the independent axis."""

        self._can_apply_model(modelfunc)
        return modelfunc(*self.get_indep())

    def eval_model_to_fit(self,
                          modelfunc: ModelFunc) -> ArrayType:
        """Evaluate the model on the independent axis after filtering."""

        self._can_apply_model(modelfunc)
        return modelfunc(*self.get_indep(filter=True))

    def to_guess(self) -> tuple[Optional[np.ndarray], ...]:

        # Should this also check whether the independent and dependent
        # axes are set?
        #
        arrays: list[Optional[np.ndarray]]
        arrays = [self.get_y(filter=True, yfunc=None)]
        arrays.extend(self.get_indep(True))
        return tuple(arrays)

    def to_fit(self,
               staterrfunc: Optional[StatErrFunc] = None
               ) -> tuple[Optional[np.ndarray],
                          Optional[ArrayType],
                          Optional[np.ndarray]]:
        return (self.get_dep(True),
                self.get_staterror(True, staterrfunc),
                self.get_syserror(True))

    def __str__(self) -> str:
        """
        Return a listing of the attributes listed in self._fields and,
        if set, self._extra_fields.
        """

        fields = self._fields + self._extra_fields
        fdict = {f: getattr(self, f) for f in fields}
        return print_fields(fields, fdict)

    def __repr__(self) -> str:
        r = f'<{type(self).__name__} data set instance'
        if hasattr(self, 'name'):
            r += f" '{self.name}'"
        r += '>'
        return r


class DataSimulFit(NoNewAttributesAfterInit):
    """Store multiple data sets.

    This class lets multiple data sets be treated as a single
    dataset by concatenation. That is, if two data sets have lengths
    n1 and n2 then they can be considered as a single data set of
    length n1 + n2.

    Parameters
    ----------
    name : str
        The name for the collection of data.
    datasets : sequence of Data objects
        The datasets to be stored; there must be at least one. They are
        assumed to behave as `sherpa.data.Data` objects, but there is no
        check for this condition.

    Attributes
    ----------
    datasets : sequence of Data

    See Also
    --------
    sherpa.models.model.SimulFitModel

    Examples
    --------

    >>> d1 = Data1D('d1', [1, 2, 3], [10, 12, 15])
    >>> d2 = Data1D('d2', [1, 2, 5, 7], [4, 15, 9, 24])
    >>> dall = DataSimulFit('comb', (d1, d2))
    >>> yvals, _, _ = dall.to_fit()
    >>> print(yvals)
    [10 12 15  4 15  9 24]

    """

    def __init__(self,
                 name: str,
                 datasets: Sequence[Data],
                 numcores: int = 1) -> None:
        if len(datasets) == 0:
            raise DataErr('zerodatasimulfit', type(self).__name__)
        self.name = name
        self.datasets = tuple(datasets)
        self.numcores = numcores
        super().__init__()

    def eval_model_to_fit(self,
                          modelfuncs: Sequence[ModelFunc]
                          ) -> np.ndarray:
        if self.numcores == 1:
            total_model = []
            for func, data in zip(modelfuncs, self.datasets):
                tmp_model = data.eval_model_to_fit(func)
                total_model.append(tmp_model)

            return np.concatenate(total_model)

        # best to make this a different derived class
        funcs = []
        datasets = []
        for func, data in zip(modelfuncs, self.datasets):
            funcs.append(func)
            datasets.append(data.get_indep(filter=False))
        total_model = parallel_map_funcs(funcs, datasets, self.numcores)
        all_model = []
        for model, data in zip(total_model, self.datasets):
            all_model.append(data.apply_filter(model))
        return np.concatenate(all_model)

    def to_fit(self,
               staterrfunc: Optional[StatErrFunc] = None
               ) -> tuple[np.ndarray,
                          Optional[np.ndarray],
                          Optional[np.ndarray]]:
        total_dep = []
        total_staterror = []
        total_syserror = []

        no_staterror = True
        no_syserror = True

        for data in self.datasets:
            dep, staterror, syserror = data.to_fit(staterrfunc)

            total_dep.append(dep)

            if staterror is not None:
                no_staterror = False
            total_staterror.append(staterror)

            if syserror is not None:
                no_syserror = False
            else:
                syserror = np.zeros_like(dep)
            total_syserror.append(syserror)

        total_dep = np.concatenate(total_dep)

        if no_staterror:
            total_staterror = None
        elif np.any([np.equal(array, None).any()
                     for array in total_staterror]):
            raise DataErr('staterrsimulfit')
        else:
            total_staterror = np.concatenate(total_staterror)

        if no_syserror:
            total_syserror = None
        else:
            total_syserror = np.concatenate(total_syserror)

        return total_dep, total_staterror, total_syserror

    # DATA-NOTE: this implementation is weird. Is this even used?
    # Typing of this is a bit hard too since Data does not have
    # a to_plot method.
    #
    def to_plot(self,
                yfunc=None,
                staterrfunc: Optional[StatErrFunc] = None
                ):
        return self.datasets[0].to_plot(yfunc.parts[0], staterrfunc)


class Data1D(Data):
    '''Dateset for 1D data.

    Parameters
    ----------
    name : string
        name of this dataset
    x : array-like
        Coordinates of the independent variable
    y : array-like
        The values of the dependent observable. If this is a numpy masked
        array, the mask will used to initialize a mask.
    staterror : array-like
        the statistical error associated with the data
    syserror : array-like
        the systematic error associated with the data
    '''
    _fields: FieldsType = ("name", "x", "y", "staterror", "syserror")
    ndim = 1

    def __init__(self,
                 name: str,
                 x: Optional[ArrayType],
                 y: Optional[ArrayType],
                 staterror: Optional[ArrayType] = None,
                 syserror: Optional[ArrayType] = None
                 ) -> None:
        self._xlabel = 'x'
        super().__init__(name, (x, ), y, staterror, syserror)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the data
        """
        return html_data1d(self)

    def _init_data_space(self,
                         filter: Filter,
                         *data: ArrayType
                         ) -> DataSpace1D:
        ndata = len(data)
        if ndata != 1:
            raise DataErr("wrongaxiscount", self.name, 1, ndata)

        ds = DataSpace1D(filter, *data)
        self._check_data_space(ds)
        return ds

    def get_x(self,
              filter: bool = False,
              model: Optional[ModelFunc] = None,
              use_evaluation_space: bool = False
              ) -> Optional[np.ndarray]:

        if model is not None:
            mdim = getattr(model, "ndim", None)
            if mdim is not None and mdim != 1:
                raise DataErr(f"Data and model dimensionality do not match: 1D and {mdim}D")

        return self.get_evaluation_indep(filter, model, use_evaluation_space)[0]

    def get_xerr(self,
                 filter: bool = False,
                 yfunc=None
                 ) -> Optional[np.ndarray]:
        """Return linear view of bin size in independent axis/axes.

        Parameters
        ----------
        filter
        yfunc

        Returns
        -------
        xerr : array or None

        """
        return None

    def get_xlabel(self) -> str:
        """Return label for linear view of independent axis/axes

        Returns
        -------
        label : str

        See Also
        --------
        set_xlabel

        """
        return self._xlabel

    def set_xlabel(self, label: str) -> None:
        """Set the label for the independent axis.

        .. versionadded:: 4.17.0

        Parameters
        -------
        label : str
           The new label.

        See Also
        --------
        set_xlabel

        """
        self._xlabel = label

    # The superclass suggests returning both the independent and
    # dependent axis sizes.
    #
    def get_dims(self, filter: bool = False) -> tuple[int, ...]:
        if self.size is None:
            raise DataErr("sizenotset", self.name)

        return len(self.get_x(filter)),

    @overload
    def get_y(self,
              filter: bool,
              yfunc: None = None,
              use_evaluation_space: bool = False
              ) -> np.ndarray:
        ...

    @overload
    def get_y(self,
              filter: bool,
              yfunc: ModelFunc,
              use_evaluation_space: bool = False
              ) -> tuple[np.ndarray, ArrayType]:
        ...

    def get_y(self, filter=False, yfunc=None,
              use_evaluation_space=False):
        """Return the dependent axis.

        Parameters
        ----------
        filter
        yfunc
        use_evaluation_space

        Returns
        -------
        y : array or (array, array) or None
            If yfunc is not None and the dependent axis is set then
            the return value is (y, y2) where y2 is yfunc evaluated on
            the independent axis.

        """
        y = self.get_dep(filter)
        if yfunc is None:
            return y

        mdim = getattr(yfunc, "ndim", None)
        if mdim is not None and mdim != 1:
            raise DataErr(f"Data and model dimensionality do not match: 1D and {mdim}D")

        model_evaluation = yfunc(*self.get_evaluation_indep(filter, yfunc, use_evaluation_space))
        return (y, model_evaluation)

    @overload
    def get_bounding_mask(self) -> tuple[bool, None]:
        ...

    @overload
    def get_bounding_mask(self) -> tuple[np.ndarray, tuple[int]]:
        ...

    def get_bounding_mask(self):
        mask = self.mask
        size = None
        if np.iterable(self.mask):
            # create bounding box around noticed image regions
            mask = np.array(self.mask)
            size = (mask.size,)

        return mask, size

    def get_img(self, yfunc=None):
        """Return 1D dependent variable as a 1 x N image.

        Parameters
        ----------
        yfunc

        Returns
        -------

        """
        y_img = self.get_y(False, yfunc)
        if yfunc is not None:
            y_img = (y_img[0].reshape(1, y_img[0].size),
                     y_img[1].reshape(1, y_img[1].size))
        else:
            y_img = y_img.reshape(1, y_img.size)

        return y_img

    def get_imgerr(self):
        err = self.get_error()
        if err is None:
            return None

        return err.reshape(1, err.size)

    def get_filter(self,
                   format: str = '%.4f',
                   delim: str = ':'
                   ) -> str:
        """Return the data filter as a string.

        Parameters
        ----------
        format : str, optional
            The formatting of the numeric values.
        delim : str, optional
            The string used to mark the low-to-high range.

        Returns
        -------
        filter : str
           The filter, represented as a collection of single values or
           ranges, separated by commas.

        See Also
        --------
        get_filter_expr, ignore, notice

        Examples
        --------

        >>> import numpy as np
        >>> x = np.asarray([1, 2, 3, 5, 6])
        >>> y = np.ones(5)
        >>> d = Data1D('example', x, y)
        >>> d.get_filter()
        '1.0000:6.0000'
        >>> d.ignore(2.5, 4.5)
        >>> d.get_filter()
        '1.0000:2.0000,5.0000:6.0000'

        >>> d.get_filter(format='%i', delim='-')
        '1-2,5-6'

        """

        x = self.get_x(filter=True)
        if np.iterable(self.mask):
            mask = self.mask
        else:
            mask = np.ones(len(x), dtype=bool)

        return create_expr(x, mask=mask, format=format, delim=delim)

    def get_filter_expr(self) -> str:
        """Return the data filter as a string along with the units.

        This is a specialised version of get_filter which adds the
        axis units.

        Returns
        -------
        filter : str
           The filter, represented as a collection of single values or
           ranges, separated by commas.

        See Also
        --------
        get_filter

        Examples
        --------

        >>> d = Data1D('example', [1., 2., 3., 5., 6., 7.], [0, .4, .5, .6, .7, .8])
        >>> d.notice(1., 6.)
        >>> d.ignore(2.5, 4.)
        >>> d.get_filter_expr()
        '1.0000-2.0000,5.0000-6.0000 x'

        Note that the expression lists the valid data points. While we ignore
        only the range 2.5-4.0, there is no data point between 4. and 5., so
        the second part of the valid range is 5.0 to 6.0.
        """
        return self.get_filter(delim='-') + ' ' + self.get_xlabel()

    # The return value is hard to type, so it is skipped at the
    # moment.
    def to_plot(self,
                yfunc: Optional[ModelFunc] = None,
                staterrfunc: Optional[StatErrFunc] = None):
        # As we introduced models defined on arbitrary grids, the x array can also depend on the
        # model function, at least in principle.
        return (self.get_x(True, yfunc),
                self.get_y(True, yfunc),
                self.get_yerr(True, staterrfunc),
                self.get_xerr(True, yfunc),
                self.get_xlabel(),
                self.get_ylabel())

    # The return value is hard to type, so it is skipped at the
    # moment.
    def to_component_plot(self,
                          yfunc: Optional[ModelFunc] = None,
                          staterrfunc: Optional[StatErrFunc] = None):
        # As we introduced models defined on arbitrary grids, the x array can also depend on the
        # model function, at least in principle.
        return (self.get_x(True, yfunc, use_evaluation_space=True),
                self.get_y(True, yfunc, use_evaluation_space=True),
                self.get_yerr(True, staterrfunc),
                self.get_xerr(True, yfunc),
                self.get_xlabel(),
                self.get_ylabel())

    def get_evaluation_indep(self,
                             filter: bool = False,
                             model: Optional[ModelFunc] = None,
                             use_evaluation_space: bool = False
                             ) -> Optional[np.ndarray]:
        data_space = self._data_space.get(filter)
        if use_evaluation_space:
            return data_space.for_model(model).grid

        return data_space.grid

    def notice(self,
               xlo: Optional[float] = None,
               xhi: Optional[float] = None,
               ignore: bool = False) -> None:
        """Notice or ignore the given range.

        Ranges are inclusive for both the lower and upper limits.

        Parameters
        ----------
        xlo, xhi : number or None, optional
            The range to change. A value of None means the minimum or
            maximum permitted value.
        ignore : bool, optional
            Set to True if the range should be ignored. The default is
            to notice the range.

        See Also
        --------
        get_filter, get_filter_expr

        Notes
        -----
        If no ranges have been ignored then a call to `notice` with
        `ignore=False` will select just the `lo` to `hi` range, and
        exclude any bins outside this range. If there has been a
        filter applied then the range `lo` to `hi` will be added to
        the range of noticed data (when `ignore=False`).

        Examples
        --------

        >>> import numpy as np
        >>> x = np.arange(0.4, 2.6, 0.2)
        >>> y = np.ones_like(x)
        >>> d = Data1D('example', x, y)
        >>> d.x[0], d.x[-1]
        (0.4, 2.4000000000000004)
        >>> d.notice()
        >>> d.get_filter(format='%.1f')
        '0.4:2.4'
        >>> d.notice(0.8, 1.2)
        >>> d.get_filter(format='%.1f')
        '0.8:1.2'
        >>> d.notice(1.5, 2.1)
        >>> d.get_filter(format='%.1f')
        '0.8:1.2,1.6:2.0'

        """

        Data.notice(self, (xlo,), (xhi,), ignore)

    @property
    def x(self) -> Optional[np.ndarray]:
        """
        Used for compatibility, in particular for __str__ and __repr__
        """
        return self.get_x()


class Data1DAsymmetricErrs(Data1D):
    """1-D data set with asymmetric errors

    See Also
    --------
    sherpa.sim.ReSampleData

    Notes
    -----
    The elo and ehi fields are stored as delta values relative to y.

    Sherpa's data and statistic objects assume symmetric errors, and
    so the methods here sometimes return a single array for an error,
    and sometimes a pair or 2D array (giving the low and high values).
    The assumption is that the staterror field is set when the object
    is created.

    """

    _fields: FieldsType = ("name", "x", "y", "staterror", "syserror", "elo", "ehi")

    def __init__(self,
                 name: str,
                 x: Optional[ArrayType],
                 y: Optional[ArrayType],
                 elo: Optional[ArrayType],
                 ehi: Optional[ArrayType],
                 staterror=None,
                 syserror=None
                 ) -> None:
        self.elo = elo
        self.ehi = ehi
        super().__init__(name, x, y, staterror=staterror, syserror=syserror)

    # TODO: should we change get_error and get_staterror?
    #
    def get_yerr(self, filter=False, staterrfunc=None):
        """Return the y error.

        The staterrfunc argument is currently ignored.
        """

        if not filter:
            # A minor optimization
            return self.elo, self.ehi

        return self.apply_filter(self.elo), self.apply_filter(self.ehi)


class Data1DInt(Data1D):
    """1-D integrated data set.

    This class is designed to handle non-consecutive bins.

    Parameters
    ----------
    name : string
        name of this dataset
    xlo : array-like
        lower bounds of the bins of the independent coordinates
    xhi : array-like
        Upper bound of the bins of the independent coordinates
    y : array-like
        The values of the dependent observable. If this is a numpy masked
        array, the mask will used to initialize a mask.
    staterror : array-like
        the statistical error associated with the data
    syserror : array-like
        the systematic error associated with the data
    """
    _fields: FieldsType = ("name", "xlo", "xhi", "y", "staterror", "syserror")

    def __init__(self,
                 name: str,
                 xlo: Optional[ArrayType],
                 xhi: Optional[ArrayType],
                 y: Optional[ArrayType],
                 staterror: Optional[ArrayType] = None,
                 syserror: Optional[ArrayType] = None
                 ) -> None:

        # Note: we do not call the superclass here.
        self._xlabel = 'x'
        Data.__init__(self, name, (xlo, xhi), y, staterror, syserror)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the data
        """
        return html_data1dint(self)

    def _init_data_space(self,
                         filter: Filter,
                         *data: ArrayType
                         ) -> IntegratedDataSpace1D:
        ndata = len(data)
        if ndata != 2:
            raise DataErr("wrongaxiscount", self.name, 2, ndata)

        ds = IntegratedDataSpace1D(filter, *data)
        self._check_data_space(ds)
        return ds

    def get_x(self,
              filter: bool = False,
              model: Optional[ModelFunc] = None,
              use_evaluation_space: bool = False
              ) -> np.ndarray:
        indep = self.get_evaluation_indep(filter, model, use_evaluation_space)
        if len(indep) == 1:
            # assume all data has been filtered out
            return np.asarray([])

        return (indep[0] + indep[1]) / 2.0

    def get_xerr(self,
                 filter: bool = False,
                 model: Optional[ModelFunc] = None
                 ) -> Optional[np.ndarray]:
        """Returns an X "error".

        The error value for the independent axis is not well defined
        in Sherpa.

        .. versionchanged:: 4.16.1
           The return value is now half the bin width instead of the
           full bin width.

        Parameters
        ----------
        filter : bool, optional
           Should the values be filtered to the current notice range?
        model : Model or None, optional

        Returns
        -------
        xerr : ndarray
           The half-width of each bin. Although the types suggest it
           can be None, it will always be an array, even if empty.

        """
        indep = self.get_evaluation_indep(filter, model)
        if len(indep) == 1:
            # assume all data has been filtered out
            return np.asarray([])

        # Return the half-width: see #748 #1817
        xlo, xhi = indep
        return (xhi - xlo) / 2

    def get_filter(self, format='%.4f', delim=':') -> str:
        """Return the data filter as a string.

        For each noticed range the filter is reported as
        starting at the low edge of the first bin and ends
        at the upper edge of the last bin in the range.

        .. versionchanged:: 4.14.0
           Prior to 4.14.0 the filter used the mid-point of the bin,
           not its low or high value.

        Parameters
        ----------
        format : str, optional
            The formatting of the numeric values.
        delim : str, optional
            The string used to mark the low-to-high range.

        Returns
        -------
        filter : str
           The filter, represented as a collection of ranges separated
           by commas.

        See Also
        --------
        get_filter_expr, ignore, notice

        Examples
        --------

        >>> import numpy as np
        >>> xlo = np.asarray([1, 2, 3, 5, 6])
        >>> xhi = xlo + 1
        >>> y = np.ones(5)
        >>> d = Data1DInt('example', xlo, xhi, y)
        >>> d.get_filter()
        '1.0000:7.0000'
        >>> d.ignore(2.5, 4.5)
        >>> d.get_filter()
        '1.0000:2.0000,5.0000:7.0000'

        >>> d.get_filter(format='%i', delim='-')
        '1-2,5-7'

        """

        # We could just use the middle of each bin (which is what
        # would happen if we didn't over-ride Data1D.filter) but let's
        # use the start and end values for each selected bin.
        #
        indep = self.get_evaluation_indep(filter=True)
        if len(indep) == 1:
            # assume all data has been filtered out
            return ''

        if np.iterable(self.mask):
            mask = self.mask
        else:
            # Unlike create_expr we do not need to send
            # a mask in here
            mask = None

        return create_expr_integrated(indep[0], indep[1], mask=mask,
                                      format=format, delim=delim)

    def notice(self,
               xlo: Optional[float] = None,
               xhi: Optional[float] = None,
               ignore: bool = False
               ) -> None:
        """Notice or ignore the given range.

        Ranges are inclusive for the lower limit and exclusive
        for the upper limit.

        .. versionchanged:: 4.14.0
           Filtering Data1DInt datasets has been improved to fix a
           number of corner cases. As part of this the upper limit has
           been changed to be exclusive whereas previously it was not
           obvious what the filter was doing.

        Parameters
        ----------
        xlo, xhi : number or None, optional
            The range to change. A value of None means the minimum or
            maximum permitted value.
        ignore : bool, optional
            Set to True if the range should be ignored. The default is
            to notice the range.

        See Also
        --------
        get_filter, get_filter_expr

        Notes
        -----
        If no ranges have been ignored then a call to `notice` with
        `ignore=False` will select just the `lo` to `hi` range, and
        exclude any bins outside this range. If there has been a
        filter applied then the range `lo` to `hi` will be added to
        the range of noticed data (when `ignore=False`).

        Examples
        --------

        >>> import numpy as np
        >>> edges = np.arange(0.4, 2.6, 0.2)
        >>> xlo, xhi = edges[:-1], edges[1:]
        >>> y = np.ones_like(xlo)
        >>> d = Data1DInt('example', xlo, xhi, y)
        >>> d.xlo[0], d.xhi[-1]
        (0.4, 2.400000000000001)
        >>> d.get_filter(format='%.1f')
        '0.4:2.4'
        >>> d.notice(0.8, 1.9)
        >>> d.get_filter(format='%.1f')
        '0.8:2.0'

        >>> d.notice()
        >>> d.get_filter(format='%.1f')
        '0.4:2.4'

        """

        Data.notice(self, (None, xlo), (xhi, None),
                    ignore=ignore, integrated=True)

    @property
    def xlo(self) -> Optional[np.ndarray]:
        """
        Property kept for compatibility
        """
        return self._data_space.x_axis.lo

    @property
    def xhi(self) -> Optional[np.ndarray]:
        """
        Property kept for compatibility
        """
        return self._data_space.x_axis.hi


class Data2D(Data):
    '''2D data set

    This class represents a 2D data set. It is desigend to work with flattened
    arrays for coordinates and independent variables, which makes it easy to
    deal with filters and sparse or irregularly-placed grids.

    Of course, the same structure can also be used for regularly-gridded data,
    it just has to be passed in as a flattened array. However, in this case
    the more specialized `~sherpa.astro.data.DataIMG` class might be more appropriate.

    Parameters
    ----------
    name : string
        name of this dataset
    x0 : array-like
        Independent coordinate values for the first dimension
    x1 : array-like
        Independent coordinate values for the second dimension
    y : array-like
        The values of the dependent observable. If this is a numpy masked
        array, the mask will be used to initialize a mask.
    shape : tuple
        Shape of the data grid for regularly gridded data (optional). This is
        used to return the data as an image e.g. for display,
        but it not needed for fitting and modelling.
        For irregularly gridded data, shape must be `None`.
    staterror : array-like
        the statistical error associated with the data
    syserror : array-like
        the systematic error associated with the data

    Examples
    --------
    An irregularly-gridded 2D dataset, with points at (-200, -200),
    (-200, 0), (0, 0), (200, -100), and (200, 150) can be created
    with:

        >>> irreg2d = Data2D("irregular2d",
        ...                  [-200, -200, 0, 200, 200],
        ...                  [-200, 0, 0, -100, 150],
        ...                  [12, 15, 23, 45, -2])

    A regularly-gridded 2D dataset can be created, but the
    arguments must be flattened::

        >>> import numpy as np
        >>> x1, x0 = np.mgrid[20:30:2, 5:20:2]
        >>> datashape = x0.shape
        >>> y = np.sqrt((x0 - 10)**2 + (x1 - 31)**2)
        >>> x0 = x0.flatten()
        >>> x1 = x1.flatten()
        >>> y = y.flatten()
        >>> reg2d = Data2D("regular2d", x0, x1, y, shape=datashape)

    .. note::
        Sherpa provides the `~sherpa.astro.data.DataIMG` class to handle
        regularly-gridded data more easily.
    '''
    _fields: FieldsType = ("name", "x0", "x1", "y", "shape", "staterror", "syserror")
    ndim = 2

    # Why should we add shape to extra-fields instead? See #1359 to
    # fix #47.
    #
    # It would change the notebook output slightly so we don't change
    # for now.
    #
    # _extra_fields = ("shape", )

    def __init__(self,
                 name: str,
                 x0: Optional[ArrayType],
                 x1: Optional[ArrayType],
                 y: Optional[ArrayType],
                 shape: Optional[tuple[int, int]] = None,
                 staterror: Optional[ArrayType] = None,
                 syserror: Optional[ArrayType] = None
                 ) -> None:
        self.shape = shape

        self._x0label = 'x0'
        self._x1label = 'x1'

        super().__init__(name, (x0, x1), y, staterror, syserror)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the data
        """
        return html_data2d(self)

    def _init_data_space(self,
                         filter: Filter,
                         *data: ArrayType
                         ) -> DataSpace2D:
        ndata = len(data)
        if ndata != 2:
            raise DataErr("wrongaxiscount", self.name, 2, ndata)

        ds = DataSpace2D(filter, *data)
        self._check_data_space(ds)
        return ds

    def get_x0(self, filter: bool = False) -> Optional[np.ndarray]:
        return self._data_space.get(filter).x0

    def get_x1(self, filter: bool = False) -> Optional[np.ndarray]:
        return self._data_space.get(filter).x1

    def get_x0label(self) -> str:
        """Return label for first dimension in 2-D view of independent axis/axes.

        Returns
        -------
        label: str

        See Also
        --------
        get_x1label, set_x0label

        """
        return self._x0label

    def get_x1label(self) -> str:
        """Return label for second dimension in 2-D view of independent axis/axes.
        Returns
        -------
        label: str

        See Also
        --------
        get_x0label, set_x1label

        """
        return self._x1label

    def set_x0label(self, label: str) -> None:
        """Set the label for the first independent axis.

        .. versionadded:: 4.17.0

        Parameters
        -------
        label: str

        See Also
        --------
        get_x0label, set_x1label

        """
        self._x0label = label

    def set_x1label(self, label: str) -> None:
        """Set the label for the second independent axis.

        .. versionadded:: 4.17.0

        Parameters
        -------
        label: str

        See Also
        --------
        get_x1label, set_x0label

        """
        self._x1label = label

    def get_axes(self) -> tuple[np.ndarray, np.ndarray]:
        self._check_shape()
        # FIXME: how to filter an axis when self.mask is size of self.y?
        return (np.arange(self.shape[1]) + 1,
                np.arange(self.shape[0]) + 1)

    def get_dims(self, filter: bool = False) -> tuple[int, ...]:
        if self.size is None:
            raise DataErr("sizenotset", self.name)

        # self._check_shape()
        if self.shape is not None:
            return self.shape[::-1]

        return len(self.get_x0(filter)), len(self.get_x1(filter))

    def get_filter_expr(self) -> str:
        return ''

    def get_filter(self) -> str:
        return ''

    def get_max_pos(self,
                    dep: Optional[np.ndarray] = None
                    ) -> Union[tuple[float, float],
                               list[tuple[float, float]]]:
        """Return the coordinates of the maximum value.

        Parameters
        ----------
        dep : ndarray or None, optional
            The data to search and it must match the current data
            filter. If not given then the dependent axis is used.

        Returns
        -------
        coords : pair or list of pairs
            The coordinates of the maximum location. The data values
            match the values returned by `get_x0` and `get_x1`. If
            there is only one maximum pixel then a pair is returned
            otherwise a list of pairs is returned.

        See Also
        --------
        get_dep, get_x0, get_x1

        """
        if dep is None:
            dep = self.get_dep(True)

        x0 = self.get_x0(True)
        x1 = self.get_x1(True)

        # TODO: Should be able to just use numpy.argmax
        pos = np.asarray(np.where(dep == dep.max())).squeeze()
        if pos.ndim == 0:
            pos = int(pos)
            return x0[pos], x1[pos]

        return [(x0[index], x1[index]) for index in pos]

    # For images, only need y-array
    # Also, we do not filter, as imager needs M x N (or
    # L x M x N) array
    def get_img(self, yfunc=None):
        """Return the dependent axis as a 2D array.

        The data is not filtered.

        Parameters
        ----------
        yfunc : sherpa.models.model.Model instance or None, optional
            If set then it is a model that is evaluated on the data
            grid and returned along with the dependent axis.

        Returns
        -------
        img : ndarray or (ndarray, ndarray)
            The data as a 2D array or a pair of 2D arrays when yfunc
            is set.

        """

        self._check_shape()
        y_img = self.get_y(False, yfunc)
        if yfunc is None:
            return y_img.reshape(*self.shape)

        return (y_img[0].reshape(*self.shape),
                y_img[1].reshape(*self.shape))

    def get_imgerr(self):
        self._check_shape()
        err = self.get_error()
        if err is None:
            return None

        return err.reshape(*self.shape)

    def to_contour(self, yfunc=None):
        return (self.get_x0(True),
                self.get_x1(True),
                self.get_y(True, yfunc),
                self.get_x0label(),
                self.get_x1label())

    def _check_shape(self) -> None:
        if self.shape is None:
            raise DataErr('shape', self.name)

    @property
    def x0(self) -> Optional[np.ndarray]:
        """
        kept for compatibility
        """
        return self.get_x0()

    @property
    def x1(self) -> Optional[np.ndarray]:
        """
        kept for compatibility
        """
        return self.get_x1()

    def notice(self,
               x0lo: Optional[float] = None,
               x0hi: Optional[float] = None,
               x1lo: Optional[float] = None,
               x1hi: Optional[float] = None,
               ignore: bool = False
               ) -> None:
        Data.notice(self, (x0lo, x1lo), (x0hi, x1hi),
                    ignore=ignore)


class Data2DInt(Data2D):
    '''2-D integrated data set

    This class represents a 2D data set. It is desigend to work with flattened
    arrays for coordinates and independent variables, which makes it easy to
    deal with filters and sparse or irregularly-placed grids.

    Of course, the same structure can also be used for regularly-gridded data,
    it just has to be passed in as a flattened array.  However, in this case
    the more specialiazed `~sherpa.astro.data.DataIMGInt` class might be more appropriate.

    Parameters
    ----------
    name : string
        name of this dataset
    x0lo : array-like
        Lower bounds of the bins in the first dimension of the independent coordinate
    x1lo : array-like
        Lower bound of the bins in the second dimension of the independent coordinate
    x0hi : array-like
        Upper bound of the bins in the first dimension of the independent coordinate
    x1hi : array-like
        Upper bound of the bins in the second dimension of the independent coordinate
    y : array-like
        The values of the dependent observable. If this is a numpy masked
        array, the mask will be used to initialize a mask.
    shape : tuple
        Shape of the data grid for regularly gridded data (optional). This is
        used return the data as an image e.g. for display,
        but it not needed for fitting and modelling.
        For irregularly gridded data, shape must be `None`.
    staterror : array-like
        the statistical error associated with the data
    syserror : array-like
        the systematic error associated with the data

    Examples
    --------
    Create an irregularly-gridded 2D dataset that has measurements in three areas
    of different size. The first is in a box of size 10 centered on
    (-15, -20), the second is a small rectangle centered on (0, 0), and the third
    another rectangle. The measured data values in these areas are 12, 15 and 23.

        >>> x0lo = [-20, -1, 10]
        >>> x1lo = [-25, -1, -10]
        >>> x0hi = [-10, 1, 15]
        >>> x1hi = [-10, 1, -8]
        >>> y = [12, 15, 23]
        >>> three_rect = Data2DInt("irregular2d_binned", x0lo=x0lo, x1lo=x1lo,
        ...                        x0hi=x0hi, x1hi=x1hi, y=y)

    The regular version of this might represent a binned image, here with
    square pixels. In this example, we create ``x1, x0`` arrays that represent
    the pixel center, but need to pass in the pixel edges to make the data object::

        >>> import numpy as np
        >>> x1, x0 = np.mgrid[20:30, 5:20]
        >>> datashape = x0.shape
        >>> y = np.sqrt((x0 - 10)**2 + (x1 - 31)**2)
        >>> x0 = x0.flatten()
        >>> x1 = x1.flatten()
        >>> y = y.flatten()
        >>> image = Data2DInt("binned_image", x0lo=x0 - 0.5, x1lo=x1 - 0.5,
        ...                   x0hi=x0 + 0.5, x1hi=x1 + 0.5, y=y, shape=datashape)

    .. note::
        Sherpa provides the `~sherpa.astro.data.DataIMGInt` class to make it easier
        to work with regularly-gridded data.
    '''
    _fields: FieldsType = ("name", "x0lo", "x1lo", "x0hi", "x1hi", "y", "staterror", "syserror")
    _extra_fields: FieldsType = ("shape", )

    def __init__(self,
                 name: str,
                 x0lo: Optional[ArrayType],
                 x1lo: Optional[ArrayType],
                 x0hi: Optional[ArrayType],
                 x1hi: Optional[ArrayType],
                 y: Optional[ArrayType],
                 shape: Optional[Sequence[int]] = None,
                 staterror: Optional[ArrayType] = None,
                 syserror: Optional[ArrayType] = None
                 ) -> None:

        # Note: we do not call the superclass here.
        self.shape = shape
        self._x0label = 'x0'
        self._x1label = 'x1'
        self._ylabel = 'y'
        Data.__init__(self, name, (x0lo, x1lo, x0hi, x1hi), y, staterror, syserror)

    def _init_data_space(self,
                         filter: Filter,
                         *data: ArrayType
                         ) -> IntegratedDataSpace2D:
        ndata = len(data)
        if ndata != 4:
            raise DataErr("wrongaxiscount", self.name, 4, ndata)

        ds = IntegratedDataSpace2D(filter, *data)
        self._check_data_space(ds)
        return ds

    def get_x0(self, filter: bool = False) -> Optional[np.ndarray]:
        if self.size is None:
            return None
        indep = self._data_space.get(filter)
        return (indep.x0lo + indep.x0hi) / 2.0

    def get_x1(self, filter=False) -> Optional[np.ndarray]:
        if self.size is None:
            return None
        indep = self._data_space.get(filter)
        return (indep.x1lo + indep.x1hi) / 2.0

    def notice(self,
               x0lo: Optional[float] = None,
               x0hi: Optional[float] = None,
               x1lo: Optional[float] = None,
               x1hi: Optional[float] = None,
               ignore: bool = False
               ) -> None:
        Data.notice(self, (None, None, x0lo, x1lo), (x0hi, x1hi, None, None),
                    ignore=ignore, integrated=True)

    @property
    def x0lo(self) -> Optional[np.ndarray]:
        """
        Property kept for compatibility
        """
        return self._data_space.x0lo

    @property
    def x0hi(self) -> Optional[np.ndarray]:
        """
        Property kept for compatibility
        """
        return self._data_space.x0hi

    @property
    def x1lo(self) -> Optional[np.ndarray]:
        """
        Property kept for compatibility
        """
        return self._data_space.x1lo

    @property
    def x1hi(self) -> Optional[np.ndarray]:
        """
        Property kept for compatibility
        """
        return self._data_space.x1hi


# Notebook representations
#
def html_data1d(data: Data1D) -> str:
    """HTML representation: Data1D

    If have matplotlib then plot the data, otherwise summarize it.
    """

    from sherpa import plot

    dtype = type(data).__name__

    plotter = plot.DataPlot()
    plotter.prepare(data)

    summary = f'{dtype} Plot'
    try:
        out = plot.backend.as_html_plot(plotter, summary)
    except AttributeError:
        out = None

    if out is not None:
        return formatting.html_from_sections(data, [out])

    # Summary properties
    #
    meta: list[tuple[str, Any]] = []
    if data.name is not None and data.name != '':
        meta.append(('Identifier', data.name))

    meta.append(('Number of bins', len(data.x)))

    # Should this only be displayed if a filter has been applied?
    #
    fexpr = data.get_filter_expr()
    nbins = data.get_dep(filter=True).size
    meta.append(('Using', f'{fexpr} with {nbins} bins'))

    # Rely on the _fields ordering, ending at staterror. This
    # ignores _extra_fields.
    #
    for f in data._fields[1:]:
        if f == 'staterror':
            break

        meta.append((f.upper(), getattr(data, f)))

    if data.staterror is not None:
        meta.append(('Statistical error', data.staterror))

    if data.syserror is not None:
        meta.append(('Systematic error', data.syserror))

    ls = [formatting.html_section(meta, summary=f'{dtype} Summary',
                                  open_block=True)]
    return formatting.html_from_sections(data, ls)


def html_data1dint(data: Data1DInt) -> str:
    """HTML representation: Data1DInt

    If have matplotlib then plot the data, otherwise summarize it.
    """

    from sherpa import plot

    dtype = type(data).__name__

    plotter = plot.DataHistogramPlot()
    plotter.prepare(data)

    summary = f'{dtype} Plot'
    try:
        out = plot.backend.as_html_plot(plotter, summary)
    except AttributeError:
        out = None

    if out is not None:
        return formatting.html_from_sections(data, [out])

    # Summary properties
    #
    meta: list[tuple[str, Any]] = []
    if data.name is not None and data.name != '':
        meta.append(('Identifier', data.name))

    meta.append(('Number of bins', len(data.xlo)))

    # Should this only be displayed if a filter has been applied?
    #
    fexpr = data.get_filter_expr()
    nbins = data.get_dep(filter=True).size
    meta.append(('Using', f'{fexpr} with {nbins} bins'))

    # Rely on the _fields ordering, ending at staterror
    for f in data._fields[1:]:
        if f == 'staterror':
            break

        meta.append((f.upper(), getattr(data, f)))

    if data.staterror is not None:
        meta.append(('Statistical error', data.staterror))

    if data.syserror is not None:
        meta.append(('Systematic error', data.syserror))

    ls = [formatting.html_section(meta, summary=f'{dtype} Summary',
                                  open_block=True)]
    return formatting.html_from_sections(data, ls)


def html_data2d(data: Data2D) -> str:
    """HTML representation: Data2D and derived classes

    """

    dtype = type(data).__name__

    # It would be nice to plot the plot, but there are several questions to
    # resolve, such as:
    #
    #  - do we plot each point (okay for sparse data) or binned
    #  - simple binning, adaptive binning, hexagonal binning?
    #  - do we just pick a number, like 100, to bin the data to

    # Summary properties
    #
    meta: list[tuple[str, Any]] = []
    if data.name is not None and data.name != '':
        meta.append(('Identifier', data.name))

    # NOTE: shape is not well defined, is it x by y or
    # the inverse, so I am not going to display it at the
    # moment.
    # if data.shape != None:
    #     meta.append(('Shape', data.shape))

    meta.append(('Number of bins', len(data.y)))

    # Rely on the _fields ordering, ending at shape
    for f in data._fields[1:]:
        if f == 'shape':
            break

        meta.append((f.upper(), getattr(data, f)))

    if data.staterror is not None:
        meta.append(('Statistical error', data.staterror))

    if data.syserror is not None:
        meta.append(('Systematic error', data.syserror))

    ls = [formatting.html_section(meta, summary=f'{dtype} Summary',
                                  open_block=True)]
    return formatting.html_from_sections(data, ls)
