# -*- coding: utf-8 -*-
#  Copyright (C) 2017, 2018  Smithsonian Astrophysical Observatory
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
import warnings

import numpy as np

from sherpa.astro.utils import reshape_2d_arrays
from sherpa.utils import interpolate, neville, rebin
from sherpa.utils.err import ModelErr


class Axis(object):
    def __init__(self, lo, hi):
        self.lo = np.asarray(lo) if lo is not None else None
        self.hi = np.asarray(hi) if hi is not None else None

    @property
    def is_empty(self):
        """Is the axis empty or None?"""
        return self.lo is None or not self.lo.size

    @property
    def is_integrated(self):
        return self.hi is not None and self.hi.size > 0

    @property
    def is_ascending(self):
        try:
            return self.lo[-1] > self.lo[0]
        except TypeError:
            raise ValueError("{} does not seem to be an array".format(self.lo))

    @property
    def start(self):
        if self.is_ascending:
            return self.lo[0]
        return self.lo[-1]

    @property
    def end(self):
        if self.is_ascending and self.is_integrated:
            return self.hi[-1]
        if self.is_ascending and not self.is_integrated:
            return self.lo[-1]
        if self.is_integrated:
            return self.hi[0]
        return self.lo[0]

    @property
    def size(self):
        return self.lo.size

    def overlaps(self, other):
        """
        Check if this axis overlaps with another
        Parameters
        ----------
        other : Axis

        Returns
        -------
        overlaps : bool
            True if they overlap, False if not
        """
        num = max(0, min(self.end, other.end) - max(self.start, other.start))
        return bool(num != 0)


class EvaluationSpace1D(object):
    def __init__(self, x=None, xhi=None):
        self.x_axis = Axis(x, xhi)

    @property
    def is_empty(self):
        """Is the grid empty or None?"""
        return self.x_axis.is_empty

    @property
    def is_integrated(self):
        """Is the grid integrated (True) or point (False)?"""
        return self.x_axis.is_integrated

    @property
    def is_ascending(self):
        """Is the grid in ascending (True) or descending (False) order?"""
        return self.x_axis.is_empty

    @property
    def grid(self):
        if self.x_axis.is_integrated:
            return self.x_axis.lo, self.x_axis.hi
        else:
            return self.x_axis.lo

    @property
    def start(self):
        return self.x_axis.start
    
    @property
    def end(self):
        return self.x_axis.end

    def zeros_like(self):
        return np.zeros(self.x_axis.lo.size)

    def overlaps(self, other):
        """
        Check if this evaluation space overlaps with another
        Parameters
        ----------
        other : EvaluationSpace1D

        Returns
        -------
        overlaps : bool
            True if they overlap, False if not
        """
        return self.x_axis.overlaps(other.x_axis)


# TODO How much validation should be done here?
class EvaluationSpace2D(object):
    def __init__(self, x=None, y=None, xhi=None, yhi=None):
        # In the 2D case the arrays are redundant, as they are flattened from a meshgrid.
        # We need to clean them up first to have proper axes.
        # This may happen when an EvaluationSpace2D is instantiated using the arrays passed to
        # the calc method.
        x_unique, y_unique, xhi_unique, yhi_unique = self._clean_arrays(x, y, xhi, yhi)
        self.x_axis = Axis(x_unique, xhi_unique)
        self.y_axis = Axis(y_unique, yhi_unique)

    def _clean_arrays(self, x, y, xhi, yhi):
        return self._clean(x), self._clean(y), self._clean(xhi), self._clean(yhi)

    @staticmethod
    def _clean(array):
        if array is not None:
            # We need to take extra care not to change the order of the arrays, hence
            # the additional complexity
            array_unique, indexes = np.unique(array, return_index=True)
            return array_unique[indexes.argsort()]

    @property
    def is_empty(self):
        return self.x_axis.is_empty or self.y_axis.is_empty

    @property
    def is_integrated(self):
        """Is the grid integrated (True) or point (False)?"""
        return (not self.is_empty)\
               and self.x_axis.is_integrated\
               and self.y_axis.is_integrated

    @property
    def is_ascending(self):
        """Is the grid in ascending (True) or descending (False) order?
        Return a tuple with (is_ascending(self.x), is_ascending(self.y))
        """
        return self.x_axis.is_ascending, self.y_axis.is_ascending

    @property
    def start(self):
        return self.x_axis.start, self.y_axis.start

    @property
    def end(self):
        return self.x_axis.end, self.y_axis.end

    def overlaps(self, other):
        """
        Check if this evaluation space overlaps with another
        Note that this is more stringent for 2D, as the boundaries
        need to coincide in this case.
        (#TODO come up with actual requirements on the overlap conditions)

        Parameters
        ----------
        other : EvaluationSpace2D

        Returns
        -------
        overlaps : bool
            True if they overlap, False if not
        """
        return bool(self.x_axis.start == other.x_axis.start\
               and self.y_axis.start == other.y_axis.start\
               and self.x_axis.end == other.x_axis.end\
               and self.y_axis.end == other.y_axis.end)

    @property
    def grid(self):
        x, y = reshape_2d_arrays(self.x_axis.lo, self.y_axis.lo)
        if self.x_axis.is_integrated:
            xhi, yhi = reshape_2d_arrays(self.x_axis.hi, self.y_axis.hi)
            return x, y, xhi, yhi
        else:
            return x, y

    def zeros_like(self):
        size = self.x_axis.lo.size * self.y_axis.lo.size
        return np.zeros(size)


class ModelDomainRegridder1D(object):
    """Allow 1D models to be evaluated on a different grid.

    This class is not used directly in a model expression;
    instead it creates an instance that is used to evaluate
    the model.

    Attributes
    ----------
    method
        The function that interpolates the data from the internal
        grid onto the requested grid. The default is
        sherpa.utils.neville. This is *only* used for point
        grids, as integrated grids use a simple rebinning scheme.

    Examples
    --------

    The "internal" model (gaussian plus constant) will be
    evaluated on the grid 0 to 10 (spacing of 0.5), and then
    linearly-interpolated onto the desired grid (1 to 8,
    spacing of 0.7). In this example there is no benefit to
    this approach - it is easier just to evaluate
    ``internal_mdl`` on the grid ``x`` - but it illustrates
    the approach.

    >>> from sherpa.models import Gauss1D, Const1D
    >>> internal_mdl = Gauss1D() + Const1D()
    >>> eval_space = EvaluationSpace1D(np.arange(0, 10, 0.5))
    >>> rmdl = ModelDomainRegridder1D(eval_space)
    >>> mdl = rmdl.apply_to(internal_mdl)
    >>> x = np.arange(1, 8, 0.7)
    >>> y = mdl(x)

    """

    def __init__(self, evaluation_space=None, name='regrid1d'):
        self.name = name
        self.evaluation_space = evaluation_space\
            if evaluation_space is not None else EvaluationSpace1D()

        # The tests show that neville (for simple interpolation-style
        # analysis) is much-more accurate than linear_interp, so use
        # that. If the user cares more about speed than accuracy
        # then they can switch to sherpa.utils.linear_interp.
        # Note that I have not tested the speed, so I am just assuming
        # that linear_interp is faster than neville.
        #
        self.method = neville

    @property
    def grid(self):
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

        Evaluate the model on the internal grid and then
        interpolate onto the desired grid.

        Parameters
        ----------
        pars : sequence of numbers
            The parameter values of the model.
        modelfunc
            The model to evaluate (the calc attribute of the model)
            TODO: should this be model as we call model.calc here?
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
        must match between the requested grid and the intenal grid
        whether it is integrated or non-integrated) is too restrictive.
        """

        if self.evaluation_space.is_empty:  # Simply pass through
            return modelfunc(pars, *args, **kwargs)

        requested_eval_space = self._make_and_validate_grid(args)

        return self._evaluate(requested_eval_space, pars, modelfunc)

    def _make_and_validate_grid(self, args_array):
        """
        Validate input grid and check whether it's point or integrated.

        Parameters
        ----------
        args_array : list
            The array or arguments passed to the `call` method

        Returns
        -------
        requested_eval_space : EvaluationSpace1D
        """
        # FIXME Isn't this fragile?
        nargs = len(args_array)
        if nargs == 0:
            raise ModelErr('nogrid')

        requested_eval_space = EvaluationSpace1D(*args_array)

        # Ensure the two grids match: integrated or non-integrated.
        if self.evaluation_space.is_integrated and not requested_eval_space.is_integrated:
            raise ModelErr('needsint')
        if requested_eval_space.is_integrated and not self.evaluation_space.is_integrated:
            raise ModelErr('needspoint')

        return requested_eval_space

    def _evaluate(self, requested_space, pars, modelfunc):
        # Evaluate the model on the user-defined grid and then interpolate
        # onto the desired grid. This is based on sherpa.models.TableModel
        # but is simplified as we do not provide a fold method.
        #
        # TODO: can we use _modelfcts.integrate1d at all here?
        #

        if not requested_space.overlaps(self.evaluation_space):
            warnings.warn("requested space and evaluation space do not overlap, evaluating model to 0")
            return requested_space.zeros_like()

        if requested_space.is_integrated:
            # TODO: should there be some check that the grid size
            #       is "compatible"? Note that test_regrid1d_int_flux
            #       appears to fail if the grid width used for modelfunc
            #       is larger than the output.
            #
            y = modelfunc(pars, *self.evaluation_space.grid)
            return rebin(y, self.grid[0], self.grid[1],
                         *requested_space.grid)
        else:
            y = modelfunc(pars, self.grid)
            return interpolate(requested_space.grid, self.grid, y,
                               function=self.method)


class ModelDomainRegridder2D(object):
    """Allow 2D models to be evaluated on a different grid.

    This class is not used directly in a model expression;
    instead it creates an instance that is used to evaluate
    the model.

    Examples
    --------

    The "internal" model (gaussian plus constant) will be
    evaluated on the grid 0 to 10 (spacing of 0.5), and then
    linearly-interpolated onto the desired grid (1 to 8,
    spacing of 0.7). In this example there is no benefit to
    this approach - it is easier just to evaluate
    ``internal_mdl`` on the grid ``x, y`` - but it illustrates
    the approach.

    >>> from sherpa.models import Gauss2D, Const2D
    >>> internal_mdl = Gauss2D() + Const2D()
    >>> eval_space = EvaluationSpace2D(np.arange(0, 10, 0.5), np.arange(0, 10, 0.5))
    >>> rmdl = ModelDomainRegridder2D(eval_space)
    >>> mdl = rmdl.apply_to(internal_mdl)
    >>> x = np.arange(1, 8, 0.7)
    >>> y = np.arange(1, 8, 0.7)
    >>> x, y = reshape_2d_arrays(x, y)
    >>> z = mdl(x, y)

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
            TODO: should this be model as we call model.calc here?
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
        must match between the requested grid and the intenal grid
        whether it is integrated or non-integrated) is too restrictive.
        """

        if self.evaluation_space.is_empty:  # Simply pass through
            return modelfunc(pars, *args, **kwargs)

        requested_eval_space = self._make_and_validate_grid(args)

        return self._evaluate(requested_eval_space, pars, modelfunc)

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
        # FIXME Isn't this fragile?
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

    def _evaluate(self, requested_space, pars, modelfunc):
        # Evaluate the model on the user-defined grid and then rebin
        # onto the desired grid.
        #
        # TODO: should there be some check that the grid size
        #       is "compatible"?

        if not requested_space.overlaps(self.evaluation_space):
            warnings.warn("requested space and evaluation space do not overlap, evaluating model to 0")
            return requested_space.zeros_like()

        y = modelfunc(pars, *self.grid)
        return rebin_2d(y, self.evaluation_space, requested_space).ravel()


def rebin_2d(y, custom_space, requested_space):
    # we assume that the spaces have an integer ratio, and other regularities, for now.
    requested_x_dim = requested_space.x_axis.size
    requested_y_dim = requested_space.y_axis.size

    custom_x_dim = custom_space.x_axis.size
    custom_y_dim = custom_space.y_axis.size

    reshaped_y = y.reshape(custom_x_dim, custom_y_dim)

    return rebin_flux(reshaped_y, requested_x_dim, requested_y_dim)


def rebin_flux(array, dimensions=None, scale=None):
    """ Return the array ``array`` to the new ``dimensions`` conserving flux the flux in the bins
    The sum of the array will remain the same. From http://martynbristow.co.uk/wordpress/blog/rebinning-data/,
    released as GPL v3 © Martyn Bristow 2015. Slightly modified for Sherpa.

    >>> ar = np.array([
    ...    [0,1,2],
    ...    [1,2,3],
    ...    [2,3,4],
    ...    ])
    >>> rebin_flux(ar, (2,2))
    array([[1.5, 4.5],
           [4.5, 7.5]])

    Raises
    ------

    AssertionError
        If the totals of the input and result array don't agree, raise an error because computation may have gone wrong

    Reference
    =========
    +-+-+-+
    |1|2|3|
    +-+-+-+
    |4|5|6|
    +-+-+-+
    |7|8|9|
    +-+-+-+
    """
    if dimensions is not None:
        if isinstance(dimensions, float):
            dimensions = [int(dimensions)] * len(array.shape)
        elif isinstance(dimensions, int):
            dimensions = [dimensions] * len(array.shape)
        elif len(dimensions) != len(array.shape):
            raise RuntimeError('')
    elif scale is not None:
        if isinstance(scale, float) or isinstance(scale, int):
            dimensions = map(int, map(round, map(lambda x: x * scale, array.shape)))
        elif len(scale) != len(array.shape):
            raise RuntimeError('')
    else:
        raise RuntimeError('Incorrect parameters to rebin.\n\trebin(array, dimensions=(x,y))\n\trebin(array, scale=a')
    # print(dimensions)
    # print("Rebinning to Dimensions: %s, %s" % tuple(dimensions))
    import itertools
    dY, dX = map(divmod, map(float, array.shape), dimensions)

    result = np.zeros(dimensions)
    for j, i in itertools.product(*map(range, array.shape)):
        (J, dj), (I, di) = divmod(j * dimensions[0], array.shape[0]), divmod(i * dimensions[1], array.shape[1])
        (J1, dj1), (I1, di1) = divmod(j + 1, array.shape[0] / float(dimensions[0])),\
                               divmod(i + 1, array.shape[1] / float(dimensions[1]))

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
    allowError = 0.1
    assert array.sum() == 0 or\
           (array.sum() < result.sum() * (1 + allowError)) and\
           (array.sum() > result.sum() * (1 - allowError))
    return result
