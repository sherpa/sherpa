# -*- coding: utf-8 -*-
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

"""
Evaluate a model on a different grid to the requested one.

This is intended to support convolution-style models, where the
convolved model should be evaluated on a different grid to the
data - e.g. a larger grid, since the convolution will account
for signal outside the data range - and then be regridded to
match the desired grid.
"""

from abc import ABCMeta, abstractmethod
import itertools
import logging
import warnings

import numpy as np
from sherpa.utils._utils import rebin  # type: ignore
from sherpa.utils.akima import akima

from sherpa.astro.utils import reshape_2d_arrays
from sherpa.utils.err import DataErr, ModelErr

warning = logging.getLogger(__name__).warning


PIXEL_RATIO_THRESHOLD = 0.1


def _to_readable_array(x):
    """Convert x into a ndarray that can not be edited (or is None)."""

    if x is None:
        return None

    x = np.asarray(x).copy()
    if not np.iterable(x):
        raise DataErr("notanarray")

    if x.ndim != 1:
        raise DataErr("not1darray")

    x.setflags(write=False)
    return x


class Axis(metaclass=ABCMeta):
    """Represent an axis of a N-D object."""

    # This is set when the data is set
    _is_ascending = None

    @property
    def is_empty(self):
        """Is the axis empty?

        An empty axis is either set to `None` or a zero-element
        sequence.
        """
        return self.size == 0

    @property
    @abstractmethod
    def is_integrated(self):
        """Is the axis integrated?"""
        pass

    @property
    def is_ascending(self):
        """Is the axis ascending?

        The axis is ascending if the elements in `lo` are sorted in
        ascending order.  Only the first and last elements are
        checked, and it is assumed that the elements are sorted.
        """
        if self.is_empty:
            raise DataErr("Axis is empty or has a size of 0")

        return self._is_ascending

    @property
    @abstractmethod
    def start(self):
        """The starting point (lowest value) of the data axis."""
        pass

    @property
    @abstractmethod
    def end(self):
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
    def size(self):
        """The size of the axis."""
        pass

    def overlaps(self, other):
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

    Parameters
    ----------
    x : array_like or None
        The starting point of the axis. If `None` or `[]` then the
        data axis is said to be empty. The axis can be in ascending or
        descending order but this is not checked.

    """

    def __init__(self, x):
        self._x = _to_readable_array(x)
        if self._x is None:
            return

        nx = len(self._x)
        if nx > 0:
            self._is_ascending = self._x[-1] > self._x[0]

    @property
    def x(self):
        return self._x

    @property
    def is_integrated(self):
        return False

    @property
    def start(self):
        if self.is_ascending:
            return self.x[0]

        return self.x[-1]

    @property
    def end(self):
        if self.is_ascending:
            return self.x[-1]

        return self.x[0]

    @property
    def size(self):
        if self.x is None:
            return 0

        return self.x.size


class IntegratedAxis(Axis):
    """Represent an integrated axis of a N-D object.

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

    def __init__(self, lo, hi):
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
            self._is_ascending = lo[-1] > lo[0]

    @property
    def lo(self):
        return self._lo

    @property
    def hi(self):
        return self._hi

    @property
    def is_integrated(self):
        return True

    @property
    def start(self):
        if self.is_ascending:
            return self.lo[0]

        return self.lo[-1]

    @property
    def end(self):
        if self.is_ascending:
            return self.hi[-1]

        return self.hi[0]

    @property
    def size(self):
        if self.lo is None:
            return 0

        return self.lo.size


# TO-DO: Why does this not inherit from NoNewAttributesAfterInit?
class EvaluationSpaceND():
    """Class for N-D Evaluation Spaces.

    Some methods return flat-True by default, others do not.
    Currently, I trying to set the default behavior such that it's
    backwards compatible, but some consistency would also be useful.

    Also, if we did ways with the ability to initialize with None,
    that would help a lot and simplify the code conceptually
    (and also reduce length and make typing more useful).


    """
    def __init__(self, axes, integrated=False):
        if integrated and not len(axes) % 2 == 0:
            raise ValueError("In integrated mode, independent variables must be pairs of lo, hi.")

        self.axes = [self._clean(arg) for arg in axes]
        if not integrated:
            self.axes = [PointAxis(arg) for arg in self.axes]
        else:
            self.axes = [IntegratedAxis(arg[0], arg[1]) for arg in zip(self.axes[:len(axes)//2],
                                                                       self.axes[len(axes)//2:])]

    @staticmethod
    def _clean(array):
        if array is None:
            return None

        # We need to take extra care not to change the order of the arrays, hence
        # the additional complexity
        array_unique, indexes = np.unique(array, return_index=True)
        return array_unique[indexes.argsort()]

    @property
    def is_empty(self):
        """Is the space empty (the x axis has no elements)?"""
        return any([a.is_empty for a in self.axes])

    @property
    def is_integrated(self):
        """Is the space integrated?"""
        # We currently only support all axes to be integrated or none of them
        # so it's sufficient to check the first one here.
        return self.axes[0].is_integrated

    @property
    def grid(self, flat=True):
        """The grid representation of the space.

        Returns
        -------
        tuple
            A tuple representing the x axis. The tuple will contain
            two arrays if the dataset is integrated, one otherwise.
        """
        if self.is_integrated:
            out = tuple(axis.lo for axis in self.axes)
            out2 = tuple(axis.hi for axis in self.axes)
        else:
            out = tuple(axis.x for axis in self.axes)
            out2 = ()

        # Again: Special treatment for None. Would be so much easier
        # if we did not allow initializing with empty axes.
        if flat and not any(o is None for o in out):
            meshgrid_axes = np.meshgrid(*out)
            meshgrid_axes2 = np.meshgrid(*out2)
            # meshgrid create a copy by default, so writing into
            # the result would not propagate back to the original arrays.
            # To avoid user confusion, make sure it's not writable.
            for m in meshgrid_axes + meshgrid_axes2:
                m.flags.writeable = False
            return tuple(m.ravel() for m in meshgrid_axes) + tuple(m.ravel() for m in meshgrid_axes2)

        return out + out2

    @property
    def is_ascending(self):
        """Is the space ascending?

        Returns
        -------
        (xflag, yflag) : (bool, bool)
            True if the axis is ascending, False otherwise, for the
            x and y axes respectively.
        """
        return tuple(a.is_ascending for a in self.axes)

    @property
    def start(self):
        """The start (lowest value) of the space.

        Returns
        -------
        (xstart, ystart) : (number, number)
            The start of the x and y axis arrays, respectively.
        """
        return tuple(a.start for a in self.axes)

    @property
    def end(self):
        """The end (highest value) of the space.

        Returns
        -------
        (xend, yend) : (number, number)
            The end of the x and y axis arrays, respectively.
        """
        return tuple(a.end for a in self.axes)

    @property
    def shape(self):
        """The sizes of the x and y axes."""
        return tuple(a.size for a in self.axes)

    def zeros_like(self, flat=True):
        """Returns zeroes for each element of the space.

        Returns
        -------
        array
            A one-dimensional array.
        """
        zeros = np.zeros(self.shape)
        if flat:
            return zeros.ravel()
        return zeros

    @property
    def midpoint_grid(self):
        """The mid-points of the space.

        For non-integrated spaces this returns the X axis.

        Returns
        -------
        array
            Return the average point of the bins of integrated axes,
            for each bin, or the non-integrated x axis array.

        """
        if self.is_integrated:
            return tuple((axis.lo + axis.hi)/2 for axis in self.axes)

        return tuple(axis.x for axis in self.axes)

    def overlaps(self, other):
        """Check if this evaluation space overlaps with another.

        Parameters
        ----------
        other : EvaluationSpaceND
            The space to compare to.

        Returns
        -------
        bool
            True if they overlap, False if not
        """
        if len(self.axes) != len(other.axes):
            raise ValueError("Can only compare spaces with the same number of axes")
        return any(a.overlaps(b) for a, b in zip(self.axes, other.axes))

    def __contains__(self, other):
        """Are all elements of other within the range (start, end) of this space?

        Parameters
        ----------
        other : EvaluationSpace1D
            The space to compare to.

        Returns
        -------
        boolean
        """
        if len(self.axes) != len(other.axes):
            raise ValueError("Can only compare spaces with the same number of axes")

        # OL: I have mixed feelings about overriding this method. On one hand it makes the
        # tests more expressive and natural, on the other this method is intended to check
        # if an element is in a collection, so it's a bit of a stretch semantically.
        return all(a.start <= b.start and a.end >= b.end for a, b in zip(self.axes, other.axes))

class EvaluationSpace1D(EvaluationSpaceND):
    """Class for 1D Evaluation Spaces.

    An Evaluation Space is a set of data axes representing the data
    space over which a model can be evaluated. A 1D Evaluation Space
    has only one axis.

    Parameters
    ----------
    x : array_like or None, optional
        The data array, or the low end of the data bins if the dataset is "integrated"
    xhi : array_like or None, optional
        The high end of the data bins for integrated datasets.
    """
    def __init__(self, x=None, xhi=None):
        """The input arrays are used to instantiate a single axis."""
        if xhi is None:
            super().__init__((x,), integrated=False)
        else:
            super().__init__((x, xhi), integrated=True)


    # Multi-D spaces try to reconstruct the original 1D input arrays by
    # removing duplicates. Doe 1D we've never done that and allowed
    # several points to have the same x-value.
    # Do we want to unify that?
    # For now, keep the old behavior.
    @staticmethod
    def _clean(array):
        return array

    @property
    def is_ascending(self):
        """Is the space ascending?"""
        return super().is_ascending[0]

    #@property
    #def grid(self):
    #    """The grid representation of the space.
    #
    #    Returns
    #    -------
    #    tuple
    #        A tuple representing the x axis. The tuple will contain
    #        two arrays if the dataset is integrated, one otherwise.
    #    """
        # We can not just rely on the is_integrated setting since
        # an integrated axis can have is_integrated set to False
        #
        # TODO: should we fix this? It only affects a few corner-case tests
        #       so maybe it's something we can address upstream? Or work out
        #       why we want is_integrated to be False when the axis size is 0?
        #
    #    if self.x_axis.is_integrated:
    #        return self.x_axis.lo, self.x_axis.hi
    #
    #    if isinstance(self.x_axis, IntegratedAxis):
    #        return self.x_axis.lo,
    #
    #    return self.x_axis.x,

    @property
    def midpoint_grid(self):
        """The mid-points of the space.

        For non-integrated spaces this returns the X axis.

        Returns
        -------
        array
            Return the average point of the bins of integrated axes,
            for each bin, or the non-integrated x axis array.

        """
        return super().midpoint_grid[0]

    @property
    def start(self):
        """The start (lowest value) of the space."""
        return super().start[0]

    @property
    def end(self):
        """The end (highest value) of the space."""
        return super().end[0]

    @property
    def x_axis(self):
        """The x axis of the space.

        Or not deprecated? Used in RMFs internally.
        Deprecated: For backwards compatibility
        """
        return self.axes[0]


class EvaluationSpace2D(EvaluationSpaceND):
    """Class for 2D Evaluation Spaces.

    An Evaluation Space is a set of data axes representing the data
    space over which a model can be evaluated. A 2D Evaluation Space
    has two axes, x and y.

    Parameters
    ----------
    x, y : array_like or None, optional
        The data array, or the low end of the x and y data bins if the dataset
        is "integrated". These are not required to be the same length.
    xhi, yhi : array_like or None, optional
        The high end of the x and y data bins for integrated datasets.
    """

    def __init__(self, x=None, y=None, xhi=None, yhi=None):
        # In the 2D case the arrays are redundant, as they are flattened from a meshgrid.
        # We need to clean them up first to have proper axes.
        # This may happen when an EvaluationSpace2D is instantiated using the arrays passed to
        # the calc method.
        #
        # This means that this class does not check that x and y (if set) have
        # the same length.
        #
        # The following cleaning is now done in the Evaluation space.
        #x_unique, y_unique, xhi_unique, yhi_unique = self._clean_arrays(x, y, xhi, yhi)

        if xhi is None and yhi is None:
            super().__init__((x, y), integrated=False)
        else:
            super().__init__((x, y, xhi, yhi), integrated=True)

    @property
    def x_axis(self):
        """The x axis of the space.

        Or not deprecated? Used in RMFs internally.
        Deprecated: For backwards compatibility
        """
        return self.axes[0]

    @property
    def y_axis(self):
        """The y axis of the space.

        Or not deprecated? Used in RMFs internally.
        Deprecated: For backwards compatibility
        """
        return self.axes[1]

    # TO-DO: This would more accurately be called "same-endpoints"
    # Why are we more stringent in 2D than in 1D?
    # Do we want to change that? For now, but leave like this for backwards
    # compatibility.
    def overlaps(self, other):
        """Check if this evaluation space overlaps with another.

        Note that this is more stringent for 2D, as the boundaries
        need to coincide in this case.

        Parameters
        ----------
        other : EvaluationSpace2D
            The space to compare to.

        Returns
        -------
        bool
            True if they overlap, False if not
        """
        return bool(self.x_axis.start == other.x_axis.start
                    and self.y_axis.start == other.y_axis.start
                    and self.x_axis.end == other.x_axis.end
                    and self.y_axis.end == other.y_axis.end)

class ModelDomainRegridder1D():
    """Allow 1D models to be evaluated on a different grid.

    This class is not used directly in a model expression;
    instead it creates an instance that is used to evaluate
    the model.

    Parameters
    ----------
    evaluation_space : object or None, optional
    name : str, optional
        The default name is 'regrid1d'.
    interp : callable, optional
        The interpolation function: it should accept three arrays: the
        output x values and the x, y values to interpolate, and return
        the interpolated y values. The default is to use
        `sherpa.utils.akima.akima`.

    Examples
    --------

    The "internal" model (gaussian plus constant) will be
    evaluated on the grid 0 to 10 (spacing of 0.5), and then
    linearly-interpolated onto the desired grid (1 to 8,
    spacing of 0.7). In this example there is no benefit to
    this approach - it is easier just to evaluate
    `internal_mdl` on the grid `x` - but it illustrates
    the approach.

    >>> from sherpa.models import Gauss1D, Const1D
    >>> internal_mdl = Gauss1D() + Const1D()
    >>> eval_space = EvaluationSpace1D(np.arange(0, 10, 0.5))
    >>> rmdl = ModelDomainRegridder1D(eval_space)
    >>> mdl = rmdl.apply_to(internal_mdl)
    >>> x = np.arange(1, 8, 0.7)
    >>> y = mdl(x)

    """

    def __init__(self, evaluation_space=None, name='regrid1d', **kwargs):
        self.name = name
        self.integrate = True
        self.evaluation_space = evaluation_space if evaluation_space is not None else EvaluationSpace1D()

        self.method = kwargs.get("interp", akima)

    @property
    def method(self):
        """Interpolate the data from the internal to requested grid.

        This is *only* used for point grids, as integrated grids use a
        simple rebinning scheme. The default is
        `sherpa.utils.akima.akima`.  The callable should accept `(xout,
        xin, yin)` arguments and interpolate the `(xin, yin)` data
        onto the `xout` grid, returning the interpolated data.
        """
        return self._method

    @method.setter
    def method(self, method):
        if not callable(method):
            raise TypeError(f"method argument '{repr(method)}' is not callable")
        self._method = method

    @property
    def grid(self):
        """The grid of the associated evaluation space."""
        return self.evaluation_space.grid

    @grid.setter
    def grid(self, value):
        try:  # value is an iterable (integrated models) to be unpacked
            self.evaluation_space = EvaluationSpace1D(*value)
        except TypeError:  # value is a single array (non-integrated models)
            self.evaluation_space = EvaluationSpace1D(value)

    def apply_to(self, model):
        """Evaluate a model on a different grid."""
        from sherpa.models.model import RegridWrappedModel
        return RegridWrappedModel(model, self)

    def calc(self, pars, modelfunc, *args, **kwargs):
        """Evaluate and regrid a model

        Evaluate the model on the internal grid and then convert it
        onto the desired grid either preserving flux (rebinning) or
        via interpolation.

        Parameters
        ----------
        pars : sequence of numbers
            The parameter values of the model.
        modelfunc
            The model to evaluate (the calc attribute of the model)
        args
            The grid to interpolate the model onto. This must match the
            format of the grid attribute of the model - i.e.
            non-integrate (single array) or integrated (a pair of
            equal-sized arrays).
        kwargs
            Keyword arguments for the model.

        Notes
        -----
        If the requested grid (i.e. that defined by args) does not overlap
        the stored grid (the grid attribute) then all values are set to 0.
        However, if the grids partially overlap then there will be
        extrapolation (depending on the method).

        It is not clear yet whether the restriction on grid type (i.e.
        must match between the requested grid and the internal grid
        whether it is integrated or non-integrated) is too restrictive.

        """

        if self.evaluation_space.is_empty:  # Simply pass through
            return modelfunc(pars, *args, **kwargs)

        requested_eval_space = self._make_and_validate_grid(args)

        return self._evaluate(requested_eval_space, pars, modelfunc, **kwargs)

    def _make_and_validate_grid(self, args_array):
        """Validate input grid and check whether it's point or integrated.

        Parameters
        ----------
        args_array : list
            The array or arguments passed to the `call` method

        Returns
        -------
        requested_eval_space : EvaluationSpace1D
        """
        nargs = len(args_array)
        if nargs == 0:
            raise ModelErr('nogrid')

        requested_eval_space = EvaluationSpace1D(*args_array)

        # Ensure the two grids match: integrated or non-integrated.
        if self.evaluation_space.is_integrated and not requested_eval_space.is_integrated:
            raise ModelErr('needsint')

        if requested_eval_space.is_integrated and not self.evaluation_space.is_integrated:
            raise ModelErr('needspoint')

        if self.evaluation_space.is_integrated and requested_eval_space.is_integrated:
            lo = self.evaluation_space.grid[0]
            hi = self.evaluation_space.grid[1]
            if np.any(lo[1:] < hi[:-1]) or np.any(lo == hi):
                raise ModelErr('needsint')

        return requested_eval_space

    def eval_non_integrated(self, pars, modelfunc, data_grid, eval_grid,
                            **kwargs):
        """Interpolate the model.

        Parameters
        ----------
        pars : list of numbers
            The parameters for the model.
        modelfunc : callable
            The model to evaluate. It is called as
            `modelfunc(pars, x, **kwargs)`
        data_grid : sequence of numbers
            The grid on which to return the values.
        eval_grid : sequence of numbers
            The grid on which to evaluate the model.
        kwargs
            Any arguments to be sent to modelfunc.

        Returns
        -------
        model : ndarray
            The model values matching the data_grid bins.

        """

        # eval_grid is out of data_grid range
        if eval_grid[-1] < data_grid[0] or eval_grid[0] > data_grid[-1]:
            return np.zeros(data_grid.size)

        #
        # join all elements of data_grid within
        # eval_spaee to minimize interpolation
        #
        indices = np.where((data_grid > eval_grid[0]) &
                           (data_grid < eval_grid[-1]))
        my_eval_grid = np.unique(np.append(eval_grid, data_grid[indices]))

        y_tmp = modelfunc(pars, my_eval_grid, **kwargs)
        y_interpolate = self.method(data_grid, my_eval_grid, y_tmp)

        if y_interpolate.size == data_grid.size and \
           eval_grid[0] < data_grid[0] and eval_grid[-1] > data_grid[-1]:
            # data space all within eval_grid
            return y_interpolate

        # find indices within data_grid
        indices = np.unique(data_grid.searchsorted(my_eval_grid))
        indices = indices[np.where(indices < data_grid.size)]

        y = np.zeros(data_grid.size)
        y[indices] = y_interpolate[indices]

        return y

    def eval_integrated(self, pars, modelfunc, data_grid, eval_grid,
                        **kwargs):
        """Rebin the model onto a grid with low and high bins.

        Parameters
        ----------
        pars : list of numbers
            The parameters for the model.
        modelfunc : callable
            The model to evaluate. It is called as
            `modelfunc(pars, lo, hi, **kwargs)`
        data_grid : (sequence of numbers, sequence of numbers)
            The grid on which to return the values, as low and
            high edges.
        eval_grid : (sequence of numbers, sequence of numbers)
            The grid on which to evaluate the model, as low and
            high edges.
        kwargs
            Any arguments to be sent to modelfunc.

        Returns
        -------
        model : ndarray
            The model values matching the data_grid bins.

        """

        y = modelfunc(pars, eval_grid[0], eval_grid[1],
                      **kwargs)
        return rebin(y, eval_grid[0], eval_grid[1],
                     data_grid[0], data_grid[1])

    def _evaluate(self, data_space, pars, modelfunc, **kwargs):
        """Evaluate the model and then convert to the requested grid.

        The model is evaluated using the evaluation_space attribute and
        then converted to match the data_space argument. If the data_space
        is integrated and the integrate attribute is set then the conversion
        is done by rebinning the data (and so preserving the signal)
        otherwise the method attribute is used to interpolate the data.

        Parameters
        ----------
        data_space : EvaluationSpace1D instance
            The output grid for the model.
        pars : list of numbers
            The parameters for the model.
        modelfunc : callable
            The model to evaluate. It is called as
            modelfunc(pars, *args, **kwargs) where args is
            either one or two arguments.
        kwargs
            Any arguments to be sent to modelfunc.

        Notes
        -----
        This is based on sherpa.models.TableModel but is simplified as
        we do not provide a fold method.

        """
        # Not really sure I need this, but let's be safe
        kwargs['integrate'] = self.integrate

        eval_space = self.evaluation_space
        if data_space.is_integrated and self.integrate:
            return self.eval_integrated(pars, modelfunc,
                                        data_space.grid, eval_space.grid,
                                        **kwargs)

        # We either have integrate=False or a non-integrated bin is given.
        return self.eval_non_integrated(pars, modelfunc,
                                        data_space.midpoint_grid,
                                        eval_space.midpoint_grid,
                                        **kwargs)


class ModelDomainRegridder2D():
    """Allow 2D models to be evaluated on a different grid.

    This class is not used directly in a model expression;
    instead it creates an instance that is used to evaluate
    the model.

    Parameters
    ----------
    evaluation_space : object or None, optional
    name : str, optional
        The default name is 'regrid2d'.

    Examples
    --------

    The "internal" model (gaussian plus constant) will be
    evaluated on the grid 0 to 10 (spacing of 0.5), and then
    linearly-interpolated onto the desired grid (1 to 8,
    spacing of 0.7). In this example there is no benefit to
    this approach - it is easier just to evaluate
    ``internal_mdl`` on the grid ``x, y`` - but it illustrates
    the approach.

    It is not obvious why this example appears to fail,
    but it is left in as documentation. See issue 840 at
    https://github.com/sherpa/sherpa/issues/840

    >>> from sherpa.models import Gauss2D, Const2D
    >>> internal_mdl = Gauss2D() + Const2D()
    >>> eval_space = EvaluationSpace2D(np.arange(0, 10, 0.5), np.arange(0, 10, 0.5))
    >>> rmdl = ModelDomainRegridder2D(eval_space)
    >>> mdl = rmdl.apply_to(internal_mdl)
    >>> x = np.arange(1, 8, 0.7)
    >>> y = np.arange(1, 8, 0.7)
    >>> x, y = reshape_2d_arrays(x, y)
    >>> z = mdl(x, y)  # doctest: +SHOW_WARNINGS
    UserWarning: requested space and evaluation space do not overlap, evaluating model to 0

    """

    def __init__(self, evaluation_space=None, name='regrid2d'):
        self.name = name
        self.evaluation_space = evaluation_space\
            if evaluation_space is not None else EvaluationSpace2D()

    @property
    def grid(self):
        return self.evaluation_space.grid

    @grid.setter
    def grid(self, value):
        self.evaluation_space = EvaluationSpace2D(*value)

    def apply_to(self, model):
        """Evaluate a model on a different grid."""
        from sherpa.models.model import RegridWrappedModel
        return RegridWrappedModel(model, self)

    def calc(self, pars, modelfunc, *args, **kwargs):
        """Evaluate and regrid a model

        Evaluate the model on the internal grid and then
        interpolate onto the desired grid.

        Parameters
        ----------
        pars : sequence of numbers
            The parameter values of the model.
        modelfunc
            The model to evaluate (the calc attribute of the model)
        args
            The grid to interpolate the model onto. This must match the
            format of the grid attribute of the model - i.e.
            non-integrate (x, y arrays) or integrated (xlo, ylo, xhi, yhi).
        kwargs
            Keyword arguments for the model.

        Notes
        -----
        If the requested grid (i.e. that defined by args) does not overlap
        the stored grid (the grid attribute) then all values are set to 0.
        However, if the grids partially overlap then there will be
        extrapolation (depending on the method).

        It is not clear yet whether the restriction on grid type (i.e.
        must match between the requested grid and the internal grid
        whether it is integrated or non-integrated) is too restrictive.
        """

        if self.evaluation_space.is_empty:  # Simply pass through
            return modelfunc(pars, *args, **kwargs)

        requested_eval_space = self._make_and_validate_grid(args)

        return self._evaluate(requested_eval_space, pars, modelfunc, **kwargs)

    def _make_and_validate_grid(self, args_array):
        """
        Validate input grid and check whether it's point or integrated.

        Parameters
        ----------
        args_array : list
            The array or arguments passed to the `call` method

        Returns
        -------
        requested_eval_space : EvaluationSpace2D
        """
        nargs = len(args_array)
        if nargs == 0:
            raise ModelErr('nogrid')

        requested_eval_space = EvaluationSpace2D(*args_array)

        # Ensure the two grids match: integrated or non-integrated.
        if self.evaluation_space.is_integrated and not requested_eval_space.is_integrated:
            raise ModelErr('needsint')
        if requested_eval_space.is_integrated and not self.evaluation_space.is_integrated:
            raise ModelErr('needspoint')

        return requested_eval_space

    def _evaluate(self, requested_space, pars, modelfunc, **kwargs):
        # Evaluate the model on the user-defined grid and then rebin
        # onto the desired grid.

        if not requested_space.overlaps(self.evaluation_space):
            warnings.warn("requested space and evaluation space do not overlap, evaluating model to 0")
            return requested_space.zeros_like()

        y = modelfunc(pars, *self.grid, **kwargs)
        return rebin_2d(y, self.evaluation_space, requested_space).ravel()


def rebin_2d(y, from_space, to_space):
    to_x_dim = to_space.x_axis.size
    to_y_dim = to_space.y_axis.size

    from_x_dim = from_space.x_axis.size
    from_y_dim = from_space.y_axis.size

    if hasattr(from_space, "data_2_psf_pixel_size_ratio"):
        ratio = from_space.data_2_psf_pixel_size_ratio
        scale_x, scale_y = 1/ratio[0], 1/ratio[1]
    else:
        scale_x = from_x_dim / to_x_dim
        scale_y = from_y_dim / to_y_dim

    scale = scale_x * scale_y

    if scale == 1:
        return y

    reshaped_y = y.reshape(from_x_dim, from_y_dim)
    reshaped_scaled_y = reshaped_y / scale

    if (abs(scale_x - round(scale_x)) > PIXEL_RATIO_THRESHOLD
            or abs(scale_y - round(scale_y)) > PIXEL_RATIO_THRESHOLD):
        return rebin_no_int(reshaped_scaled_y, dimensions=(to_x_dim, to_y_dim))

    return rebin_int(reshaped_scaled_y, int(round(scale_x)), int(round(scale_y)))


def rebin_int(array, scale_x, scale_y):
    """Rebin array by an integer scale on both x and y

    Parameters
    ----------
    array : array_like
        The array to be rebinned
    scale_x : int
        The pixel ratio on the x axis
    scale_y : int
        The pixel ratio on the y axis

    Returns
    -------
    array_like

    """
    xedge = np.shape(array)[0] % scale_x
    yedge = np.shape(array)[1] % scale_y
    sub_array = array[xedge:, yedge:]
    binned_x_shape = np.shape(sub_array)[0] // scale_x
    binned_y_shape = np.shape(sub_array)[1] // scale_y

    image = np.reshape(sub_array, (binned_x_shape, scale_x, binned_y_shape, scale_y))
    image = np.sum(image, axis=3)
    image = np.sum(image, axis=1)

    return image


def rebin_no_int(array, dimensions=None, scale=None):
    """Rebin the array, conserving flux.

    Return the array ``array`` to the new ``dimensions`` conserving flux,
    so that the sum of the output matches the sum of ``array``.

    Raises
    ------
    AssertionError
        If the totals of the input and result array don't agree, raise an error because computation may have gone wrong

    Notes
    -----
    This routine is based on the example at
    https://martynbristow.co.uk/
    which was released as GPL v3 © Martyn Bristow 2015. It has been
    slightly modified for Sherpa.

    Examples
    --------

    >>> ar = np.array([
    ...    [0,1,2],
    ...    [1,2,3],
    ...    [2,3,4],
    ...    ])
    >>> rebin_no_int(ar, (2,2))
    array([[1.5, 4.5],
           [4.5, 7.5]])

    """
    if dimensions is not None:
        if isinstance(dimensions, float):
            dimensions = [int(dimensions)] * len(array.shape)
        elif isinstance(dimensions, int):
            dimensions = [dimensions] * len(array.shape)
        elif len(dimensions) != len(array.shape):
            raise RuntimeError('')
    elif scale is not None:
        if isinstance(scale, (float, int)):
            dimensions = map(int, map(round, map(lambda x: x * scale, array.shape)))
        elif len(scale) != len(array.shape):
            raise RuntimeError('')
    else:
        raise RuntimeError('Incorrect parameters to rebin.\n\trebin(array, dimensions=(x,y))\n\trebin(array, scale=a')

    dY, dX = map(divmod, map(float, array.shape), dimensions)

    result = np.zeros(dimensions)
    for j, i in itertools.product(*map(range, array.shape)):
        (J, dj) = divmod(j * dimensions[0], array.shape[0])
        (I, di) = divmod(i * dimensions[1], array.shape[1])
        (J1, dj1) = divmod(j + 1, array.shape[0] / float(dimensions[0]))
        (I1, di1) = divmod(i + 1, array.shape[1] / float(dimensions[1]))

        # Moving to new bin
        # Is this a discrete bin?
        dx, dy = 0, 0
        if (I1 - I == 0) | ((I1 - I == 1) & (di1 == 0)):
            dx = 1
        else:
            dx = 1 - di1
        if (J1 - J == 0) | ((J1 - J == 1) & (dj1 == 0)):
            dy = 1
        else:
            dy = 1 - dj1
        # Prevent it from allocating outide the array
        I_ = min(dimensions[1] - 1, I + 1)
        J_ = min(dimensions[0] - 1, J + 1)
        result[J, I] += array[j, i] * dx * dy
        result[J_, I] += array[j, i] * (1 - dy) * dx
        result[J, I_] += array[j, i] * dy * (1 - dx)
        result[J_, I_] += array[j, i] * (1 - dx) * (1 - dy)

    allowError = 0.001
    assert array.sum() == 0 or \
           (array.sum() < result.sum() * (1 + allowError)) and \
           (array.sum() > result.sum() * (1 - allowError))
    return result
