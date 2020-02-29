# -*- coding: utf-8 -*-
from __future__ import division
#  Copyright (C) 2017, 2018, 2019, 2020
#         Smithsonian Astrophysical Observatory
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
from sherpa.utils._utils import rebin

from sherpa.astro.utils import reshape_2d_arrays
from sherpa.utils import interpolate, neville
from sherpa.utils.err import ModelErr

import logging
warning = logging.getLogger(__name__).warning


PIXEL_RATIO_THRESHOLD = 0.1


class Axis():
    """
    Class for representing N-D axes objects, for both "integrated" and "non-integrated" datasets
    """
    def __init__(self, lo, hi):
        """
        In integrated datasets axes are defined by bins. In this case both `lo` and `hi` are
        not None. `lo` and `hi` will be converted to `numpy` arrays if they are not.

        If `lo` is `None` or empty then the data axis is said to be empty.

        Parameters
        ----------
        lo : array_like
            The starting point of the axis
        hi : array_like
            The ending point of the axis
        """
        self.lo = np.asarray(lo) if lo is not None else None
        self.hi = np.asarray(hi) if hi is not None else None

    @property
    def is_empty(self):
        """


        Returns
        -------
        bool
            Whether the axis is empty, i.e. if `lo` is `None` or an empty array.
        """
        return self.lo is None or not self.lo.size

    @property
    def is_integrated(self):
        """
        Is the axis integrated?

        Returns
        -------
        bool
            The axis is integrated is `hi` is not None and not empty.
        """
        return self.hi is not None and self.hi.size > 0

    @property
    def is_ascending(self):
        """
        Is the axis ascending?

        Returns
        -------
        bool
            The axis is ascending if the elements in `lo` are sorted in ascending order.
            Only the first and last elements are checked, and it is assumed that the
            elements are sorted.
        """
        try:
            return self.lo[-1] > self.lo[0]
        except TypeError:
            raise ValueError("{} does not seem to be an array".format(self.lo))

    @property
    def start(self):
        """
        Starting point of the data axis

        Returns
        -------
        number
            The first element in `lo` if the axis is ascending, or the last element otherwise.
        """
        if self.is_ascending:
            return self.lo[0]
        return self.lo[-1]

    @property
    def end(self):
        """
        Ending point of the data axis

        Returns
        -------
        number
            If the data axis is ascending the end boundary is the last element of the `hi` array when
            the axis is integrated, otherwise it's the last element of `lo`.

            Conversely, for descending axes, the last element is either the first element of the `hi`
            array or of the `lo` array, depending on whether the axis is integrated or not,
            respectively.
        """
        if self.is_ascending and self.is_integrated:
            return self.hi[-1]
        if self.is_ascending and not self.is_integrated:
            return self.lo[-1]
        if self.is_integrated:
            return self.hi[0]
        return self.lo[0]

    @property
    def size(self):
        """
        The size of the axis.

        Returns
        -------
        number
            The size of the axis.
        """
        return self.lo.size

    def overlaps(self, other):
        """
        Check if this axis overlaps with another
        Parameters
        ----------
        other : Axis

        Returns
        -------
        bool
            True if they overlap, False if not
        """
        num = max(0, min(self.end, other.end) - max(self.start, other.start))
        return bool(num != 0)


class EvaluationSpace1D():
    """
    Class for 1D Evaluation Spaces. An Evaluation Space is a set of data axes representing
    the data space over which a model can be evaluated.

    A 1D Evaluation Space has only one axis.
    """
    def __init__(self, x=None, xhi=None):
        """
        The input arrays are used to instantiate a single axis.

        Parameters
        ----------
        x : array_like
            The data array, or the low end of the data bins if the dataset is "integrated"
        xhi: array_like
            The high end of the data bins for integrated datasets.
        """
        self.x_axis = Axis(x, xhi)

    @property
    def is_empty(self):
        """
        Is the dataset empty?

        Returns
        -------
        bool
            True if the x axis is empty, False otherwise
        """
        return self.x_axis.is_empty

    @property
    def is_integrated(self):
        """
        Is the grid integrated?

        Returns
        -------
        bool
            True if the x axis is integrated, False otherwise.
        """
        return self.x_axis.is_integrated

    @property
    def is_ascending(self):
        """
        Is the dataset ascending?

        Returns
        -------
        bool
            True if the x axis is ascending, False otherwise.
        """
        return self.x_axis.is_empty

    @property
    def grid(self):
        """
        Return the grid representation of this dataset. The grid is always a tuple, even if the
        dataset is 1-D and not integrated. This is due to the existing architecture of Sherpa's
        model classes and the fact that there is no signature difference among 1-D and 2-D models,
        as 1-D models can receive 1 or 2 arrays and 2-D models can receive 2 or 4 arrays.

        Returns
        -------
        tuple
            A tuple representing the x axis. The tuple will contain two arrays if the dataset is
            integrated, one otherwise.
        """
        if self.x_axis.is_integrated:
            return self.x_axis.lo, self.x_axis.hi
        else:
            return self.x_axis.lo,

    @property
    def midpoint_grid(self):
        """
        Return a single array representing the dataset.

        Returns
        -------
        array
            Return the average point of the bins of integrated axes, for each bin, or the non-integrated
            x axis array.
        """
        if self.x_axis.is_integrated:
            return (self.x_axis.lo + self.x_axis.hi)/2
        else:
            return self.x_axis.lo

    @property
    def start(self):
        """
        The start of the dataset.

        Returns
        -------
        number
            The start of the x axis array
        """
        return self.x_axis.start

    @property
    def end(self):
        """
        The end of the dataset.

        Returns
        -------
        number
            The end of the x axis array
        """
        return self.x_axis.end

    def zeros_like(self):
        """
        Utility function that returns an array of zeros that has the same shape as the dataset.

        Returns
        -------
        array
        """
        return np.zeros(self.x_axis.lo.size)

    def overlaps(self, other):
        """
        Check if this evaluation space overlaps with another
        Parameters
        ----------
        other : EvaluationSpace1D

        Returns
        -------
        bool
            True if they overlap, False if not
        """
        return self.x_axis.overlaps(other.x_axis)

    def __contains__(self, other):
        """
        check if this space properly contains the `other` space, i.e. if the `other` space is contained
        within the boundaries of `self`.

        Parameters
        ----------
        other : EvaluationSpace1D

        Returns
        -------
        boolean
        """

        # OL: I have mixed feelings about overriding this method. On one hand it makes the
        # tests more expressive and natural, on the other this method is intended to check
        # if an element is in a collection, so it's a bit of a stretch semantically.
        return self.start <= other.start and self.end >= other.end


class EvaluationSpace2D():
    """
    Class for 2D Evaluation Spaces. An Evaluation Space is a set of data axes representing
    the data space over which a model can be evaluated.

    A 2D Evaluation Space has two axes, x and y.
    """
    def __init__(self, x=None, y=None, xhi=None, yhi=None):
        """
        The input arrays are used to instantiate the x and y axes.

        Parameters
        ----------
        x : array_like
            The data array, or the low end of the x data bins if the dataset is "integrated"
        xhi: array_like
            The high end of the x data bins for integrated datasets.
        y : array_like
            The data array, or the low end of the y data bins if the dataset is "integrated"
        yhi: array_like
            The high end of the y data bins for integrated datasets.
        """

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
        """
        Is the dataset empty?

        Returns
        -------
        bool
            True if the x axis or y axis are empty, False otherwise
        """
        return self.x_axis.is_empty or self.y_axis.is_empty

    @property
    def is_integrated(self):
        """
        Is the grid integrated?

        Returns
        -------
        bool
            True if the axes are integrated, False otherwise.
        """
        return (not self.is_empty)\
               and self.x_axis.is_integrated\
               and self.y_axis.is_integrated

    @property
    def is_ascending(self):
        """
        Is the dataset ascending?

        Returns
        -------
        tuple(bool)
            True if the axis is ascending, False otherwise, for the x and y axes respectively
        """
        return self.x_axis.is_ascending, self.y_axis.is_ascending

    @property
    def start(self):
        """
        The start of the dataset.

        Returns
        -------
        tuple
            The start of the x and y axis arrays, respectively
        """
        return self.x_axis.start, self.y_axis.start

    @property
    def end(self):
        """
        The enf of the dataset.

        Returns
        -------
        tuple
            The end of the x and y axis arrays, respectively
        """
        return self.x_axis.end, self.y_axis.end

    @property
    def shape(self):
        return self.x_axis.size, self.y_axis.size

    def overlaps(self, other):
        """
        Check if this evaluation space overlaps with another
        Note that this is more stringent for 2D, as the boundaries
        need to coincide in this case.

        Parameters
        ----------
        other : EvaluationSpace2D

        Returns
        -------
        bool
            True if they overlap, False if not
        """
        return bool(self.x_axis.start == other.x_axis.start\
               and self.y_axis.start == other.y_axis.start\
               and self.x_axis.end == other.x_axis.end\
               and self.y_axis.end == other.y_axis.end)

    @property
    def grid(self):
        """
        Return the grid representation of this dataset. The grid is always a tuple, even if the
        dataset is 1-D and not integrated. This is due to the existing architecture of Sherpa's
        model classes and the fact that there is no signature difference among 1-D and 2-D models,
        as 1-D models can receive 1 or 2 arrays and 2-D models can receive 2 or 4 arrays.

        The x and y arrays in the grid are one-dimentional representations of the meshgrid obtained
        from the x and y axis arrays, as in `numpy.meshgrid(x, y)[0].ravel()`

        Returns
        -------
        tuple
            A tuple representing the x and y axes. The tuple will contain four arrays if the dataset is
            integrated, two otherwise.
        """

        x, y = reshape_2d_arrays(self.x_axis.lo, self.y_axis.lo)
        if self.x_axis.is_integrated:
            xhi, yhi = reshape_2d_arrays(self.x_axis.hi, self.y_axis.hi)
            return x, y, xhi, yhi
        else:
            return x, y

    def zeros_like(self):
        """
        Utility function that returns an array of zeros that has the same shape as the dataset.

        Returns
        -------
        array
        """
        size = self.x_axis.lo.size * self.y_axis.lo.size
        return np.zeros(size)


class ModelDomainRegridder1D():
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
        self.integrate = True
        self.evaluation_space = evaluation_space if evaluation_space is not None else EvaluationSpace1D()

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
            if np.any(lo[1:] < hi[:-1]):
                raise ModelErr('needsint')

        return requested_eval_space

    def eval_non_integrated(self, pars, modelfunc, data_space, eval_space,
                            **kwargs):
        y = np.zeros(data_space.size)

        # eval_space is out of data_space range
        if eval_space[-1] < data_space[0] or eval_space[0] > data_space[-1]:
            return y

        indices = data_space.searchsorted(eval_space)
        indices_within_data_space = np.where(indices < len(data_space))
        my_eval_space = eval_space[indices_within_data_space]
        yy = modelfunc(pars, eval_space, **kwargs)
        y_interpolate = interpolate(data_space, eval_space, yy,
                                    function=self.method)
        if y_interpolate.size == data_space.size and \
           eval_space[0] < data_space[0] and eval_space[-1] > data_space[-1]:
            return y_interpolate
        indices_within_y = np.where(indices < len(y))
        y[indices[indices_within_y]] = y_interpolate[indices[indices_within_y]]
        return y

    def _evaluate(self, data_space, pars, modelfunc, **kwargs):
        """
        Evaluate the model on the user-defined grid and then interpolate/rebin
        onto the desired grid. This is based on sherpa.models.TableModel
        but is simplified as we do not provide a fold method.
        """
        kwargs['integrate'] = self.integrate  # Not really sure I need this, but let's be safe

        eval_space = self.evaluation_space
        if data_space.is_integrated:
            if self.integrate:
                # This should be the most common case
                y = modelfunc(pars, eval_space.grid[0], eval_space.grid[1],
                              **kwargs)
                return rebin(y, eval_space.grid[0], eval_space.grid[1],
                             data_space.grid[0], data_space.grid[1])
            else:
                # The integrate flag is set to false, so just evaluate the model
                # and then interpolate using the grids midpoints.
                y = modelfunc(pars, eval_space.midpoint_grid, **kwargs)
                return interpolate(data_space.midpoint_grid,
                                   eval_space.midpoint_grid, y,
                                   function=self.method)
        else:
            return self.eval_non_integrated(pars, modelfunc,
                                            data_space.midpoint_grid,
                                            eval_space.midpoint_grid,
                                            **kwargs)


class ModelDomainRegridder2D():
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

        if not requested_space.overlaps(self.evaluation_space):
            warnings.warn("requested space and evaluation space do not overlap, evaluating model to 0")
            return requested_space.zeros_like()

        y = modelfunc(pars, *self.grid)
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
    """
    Rebin array by an integer scale on both x and y

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
    http://martynbristow.co.uk/wordpress/blog/rebinning-data/
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
        if isinstance(scale, float) or isinstance(scale, int):
            dimensions = map(int, map(round, map(lambda x: x * scale, array.shape)))
        elif len(scale) != len(array.shape):
            raise RuntimeError('')
    else:
        raise RuntimeError('Incorrect parameters to rebin.\n\trebin(array, dimensions=(x,y))\n\trebin(array, scale=a')
    import itertools
    dY, dX = map(divmod, map(float, array.shape), dimensions)

    result = np.zeros(dimensions)
    for j, i in itertools.product(*map(range, array.shape)):
        (J, dj), (I, di) = divmod(j * dimensions[0], array.shape[0]), divmod(i * dimensions[1], array.shape[1])
        (J1, dj1), (I1, di1) = divmod(j + 1, array.shape[0] / float(dimensions[0])), \
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
    allowError = 0.001
    assert array.sum() == 0 or \
           (array.sum() < result.sum() * (1 + allowError)) and \
           (array.sum() > result.sum() * (1 - allowError))
    return result
