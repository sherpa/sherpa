#
#  Copyright (C) 2017  Smithsonian Astrophysical Observatory
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
from sherpa.utils import interpolate, neville
from sherpa.utils.err import ModelErr


__all__ = ('Regrid1D', 'Regrid2D', 'RegridModel1D')


# TODO: how much of this can be in a Regrid class which can
#       be used by both 1D and 2D versions.

class Regrid1D(Model):
    """Allow 1D models to be evaluated on a different grid.

    This class is not used directly in a model expression;
    instead it creates an instance that is used to evaluate
    the model.

    Attributes
    ----------
    method
        The function that interpolates the data from the internal
        grid onto the requested grid. The default is
        sherpa.utils.neville.

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
    >>> rmdl = Regrid1D()
    >>> rmdl.grid = np.arange(0, 10, 0.5)
    >>> mdl = rmdl(internal_mdl)
    >>> x = np.arange(1, 8, 0.7)
    >>> y = mdl(x)

    """

    def __init__(self, name='regrid1d'):
        self.__grid = None
        self.__grid_is_integrated = None
        self.__grid_is_ascending = None

        # The tests show that neville (for simple interpolation-style
        # analysis) is much-more accurate than linear_interp, so use
        # that. If the user cares more about speed than accuracy
        # then they can switch to sherpa.utils.linear_interp.
        # Note that I have not tested the speed, so I am just assuming
        # that linear_interp is faster than neville.
        #
        self.method = neville

        Model.__init__(self, name, ())

    @property
    def grid_is_integrated(self):
        """Is the grid integrated (True) or point (False)?"""
        return self.__grid_is_integrated

    @property
    def grid_is_ascending(self):
        """Is the grid in ascending (True) or descending (False) order?"""
        return self.__grid_is_ascending

    @property
    def grid(self):
        """Return the grid on which the model is to be evaluate."""
        return self.__grid

    def set_grid(self, x, xhi=None):
        """Set the grid.

        Parameters
        ----------
        x, xhi : array or None
            The grid, which can either be a single array
            (xhi is None) or the low and high edges of each bin,
            in which case the two arrays must have the same number
            of elements. The arrays are assumed to be 1D and
            monotonic: that is, the first and last bins define
            the start and end of the grid. The arrays can be
            in descending order (or at least, that is the current
            plan). If x is None then the grid is removed.
        """

        if x is None:
            # does not matter whether xhi is set or not
            self.__grid_is_integrated = None
            self.__grid_is_ascending = None
            self.__grid = None
        elif xhi is None:
            x = np.asarray(x)
            self.__grid_is_integrated = False
            self.__grid_is_ascending = x[-1] > x[0]
            self.__grid = x
        else:
            x = np.asarray(x)
            xhi = np.asarray(xhi)
            self.__grid_is_integrated = True
            self.__grid_is_ascending = x[-1] > x[0]
            self.__grid = (x, xhi)

    def __call__(self, model):
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

        nargs = len(args)
        if nargs == 0:
            raise ModelErr('nogrid')

        if self.grid is None:
            # There is no need for a regrid, so just evaluate
            # the model with the input grid
            #
            return modelfunc(pars, *args, **kwargs)

        # Ensure the two grids match: integrated or non-integrated.
        #
        requested_grid_point = nargs == 1
        if requested_grid_point and self.grid_is_integrated:
            raise ModelErr('needsint')
        if not requested_grid_point and not self.grid_is_integrated:
            raise ModelErr('needspoint')

        # NOTE: for the moment focussing on both grids in ascending
        #       order. Will need to add tests to check if there is
        #       a mis-match.
        #

        # Check if there is any overlap. For now just rely on
        # the first array (even if integrated).
        #
        # Perhaps we should just let this through, and let the
        # interpolation or extrapolation happen.
        x1 = args[0][0]
        x2 = args[0][-1]
        if x2 > x1:
            x1, x2 = x2, x1

        if requested_grid_point:
            g = self.grid
        else:
            g = self.grid[0]

        g1 = g[0]
        g2 = g[-1]

        if not self.grid_is_ascending:
            g1, g2 = g2, g1

        if (g1 > x2) or (g2 < x1):
            return np.zeros(args[0].size)

        # Evaluate the model on the user-defined grid and then interpolate
        # onto the desired grid. This is based on sherpa.models.TableModel
        # but is simplified as we do not provide a fold method.
        #
        # TODO: can we use _modelfcts.integrate1d at all here?
        #
        if requested_grid_point:
            y = modelfunc(pars, self.grid)
            return interpolate(args[0], self.grid, y,
                               function=self.method)

        else:
            raise NotImplementedError()


class Regrid2D(Model):

    def __init__(self):
        raise NotImplementedError


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
                                (self.wrapper, self.model))

    def calc(self, p, *args, **kwargs):
        return self.wrapper.calc(p, self.model.calc, *args, **kwargs)
