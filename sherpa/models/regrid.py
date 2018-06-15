#
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

.. note::

   This is an experimental module.

"""

import numpy as np

from sherpa.models.model import Model, ArithmeticModel, \
    ArithmeticFunctionModel, CompositeModel
from sherpa.utils import interpolate, neville, rebin
from sherpa.utils.err import ModelErr


__all__ = ('Regrid1D', 'RegridModel1D')


class EvaluationSpace(object):
    def __init__(self, x=None, xhi=None):
        self.xlo = np.asarray(x) if x is not None else None
        self.xhi = self.xhi = np.asarray(xhi) if xhi is not None else None

    @property
    def is_empty(self):
        """Is the grid empty or None?"""
        return self.xlo is None or not self.xlo.size

    @property
    def is_integrated(self):
        """Is the grid integrated (True) or point (False)?"""
        return self.xhi is not None

    @property
    def is_ascending(self):
        """Is the grid in ascending (True) or descending (False) order?"""
        return self.xlo[-1] > self.xlo[0]

    @property
    def grid(self):
        if self.xhi is not None:
            return self.xlo, self.xhi
        else:
            return self.xlo

    @property
    def start(self):
        if self.is_ascending:
            return self.xlo[0]
        return self.xlo[-1]
    
    @property
    def end(self):
        if self.is_ascending and self.is_integrated:
            return self.xhi[-1]
        if self.is_ascending and not self.is_integrated:
            return self.xlo[-1]
        if self.is_integrated:
            return self.xhi[0]
        return self.xlo[0]

    def zeros_like(self):
        return np.zeros(self.xlo.size)

    def overlaps(self, other):
        """
        Check if this evaluation space overlaps with another
        Parameters
        ----------
        other : EvaluationSpace

        Returns
        -------
        overlaps : bool
            True if they overlap, False if not
        """
        return max(0, min(self.end, other.end) - max(self.start, other.start))


class Regrid1D(object):
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

    >>> internal_mdl = Gauss1D() + Const1D()
    >>> eval_space = EvaluationSpace(np.arange(0, 10, 0.5))
    >>> rmdl = Regrid1D(eval_space)
    >>> mdl = rmdl(internal_mdl)
    >>> x = np.arange(1, 8, 0.7)
    >>> y = mdl(x)

    """

    def __init__(self, evaluation_space=None, name='regrid1d'):
        self.name = name
        self.evaluation_space = evaluation_space\
            if evaluation_space is not None else EvaluationSpace()

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
            self.evaluation_space = EvaluationSpace(*value)
        except TypeError:  # value is a single array (non-integrated models)
            self.evaluation_space = EvaluationSpace(value)

    def apply_to(self, model):
        """Evaluate a model on a different grid."""
        return RegridModel1D(model, self)

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

        if not requested_eval_space.overlaps(self.evaluation_space):
            return requested_eval_space.zeros_like()

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
        requested_eval_space : EvaluationSpace
            Whether the requested grid is point or integrated.
        """
        # FIXME Isn't this fragile?
        nargs = len(args_array)
        if nargs == 0:
            raise ModelErr('nogrid')

        requested_eval_space = EvaluationSpace(*args_array)

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


class RegridModel1D(CompositeModel, ArithmeticModel):

    @staticmethod
    def wrapobj(obj):
        if isinstance(obj, ArithmeticModel):
            return obj
        else:
            return ArithmeticFunctionModel(obj)

    def __init__(self, model, wrapper):
        self.model = self.wrapobj(model)
        self.wrapper = wrapper
        CompositeModel.__init__(self,
                                "{}({})".format(self.wrapper.name,
                                                self.model.name),
                                (self.model, ))

    def calc(self, p, *args, **kwargs):
        return self.wrapper.calc(p, self.model.calc, *args, **kwargs)
