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
from sherpa.utils import interpolate, neville, rebin
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
    >>> rmdl = Regrid1D()
    >>> rmdl.set_grid(np.arange(0, 10, 0.5))
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
        x = args[0]
        x1 = x[0]
        x2 = x[-1]
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
            return np.zeros(x.size)

        # Evaluate the model on the user-defined grid and then interpolate
        # onto the desired grid. This is based on sherpa.models.TableModel
        # but is simplified as we do not provide a fold method.
        #
        # TODO: can we use _modelfcts.integrate1d at all here?
        #
        if requested_grid_point:
            y = modelfunc(pars, self.grid)
            return interpolate(x, self.grid, y,
                               function=self.method)

        else:
            # TODO: should there be some check that the grid size
            #       is "compatible"? Note that test_regrid1d_int_flux
            #       appears to fail if the grid width used for modelfunc
            #       is larger than the output.
            #
            y = modelfunc(pars, self.grid[0], self.grid[1])
            return rebin(y, self.grid[0], self.grid[1],
                         args[0], args[1])


class Regrid2D(Model):
    """Allow 2D models to be evaluated on a different grid.

    Currently this is restricted to 2D point grids (i.e. the
    model is evaluated at a point (X0,X1) and not over a cell).
    This class is not used directly in a model expression;
    instead it creates an instance that is used to evaluate
    the model.

    Attributes
    ----------
    method
        The function that interpolates the data from the internal
        grid onto the requested grid. The default is
        sherpa.utils.neville. This is *only* used for point
        grids.

    Examples
    --------

    The "internal" model (gaussian plus constant) will be
    evaluated on the grid 0 to 10 (spacing of 0.5) for X0
    and -5 to 14 (spacing 0.4) for X1, and then interpolated
    onto the grid that the model is called with (as defined by
    the ``x0`` and ``x1`` arrays below). In this example there
    is no benefit to this approach - it is easier just to evaluate
    ``internal_mdl`` on the grid ``x0, x1`` - but it illustrates
    the approach.

    >>> internal_mdl = Gauss2D() + Const2D()
    >>> rmdl = Regrid2D()
    >>> g0 = np.arange(0, 10, 0.5)
    >>> g1 = np.arange(-5, 15, 0.4)
    >>> rmdl.set_grid)(g0, g1)
    >>> mdl = rmdl(internal_mdl)
    >>> x0 = np.arange(1, 8, 0.7)
    >>> x1 = np.arange(-2, 12, 0.7)
    >>> y = mdl(x0, x1)

    """

    def __init__(self, name='regrid2d'):
        self.__grid = None
        self.__grid_is_integrated = None
        self.__grid_is_ascending = None

        # TODO: what methods so we have for 2D interpolation?
        self.method = None

        Model.__init__(self, name, ())

    @property
    def grid_is_integrated(self):
        """Is the grid integrated (True) or point (False)?"""
        return self.__grid_is_integrated

    @property
    def grid_is_ascending(self):
        """Is the grid in ascending (True) or descending (False) order?

        The returned value is ``None`` or a tuple, for the X0 and X1
        axes.
        """
        return self.__grid_is_ascending

    @property
    def grid(self):
        """Return the grid on which the model is to be evaluate."""
        return self.__grid

    def set_grid(self, x0, x1):
        """Set the grid.

        At present only point grids are supported.

        Parameters
        ----------
        x0, x1 : array or None
            The grids for the X0 and X1 independent axes.
            Each array is assumed to be 1D and monotonic: that is,
            the first and last bins define the start and end of that
            grid. The arrays can be in descending order (or at least,
            that is the current plan). There is no requirement that
            the two grids have to have the same length or spacing.
            Either both arguments are set, or both are None.
        """

        if x0 is None and x1 is None:
            self.__grid_is_integrated = None
            self.__grid_is_ascending = None
            self.__grid = None
            return

        if x0 is None or x1 is None:
            # Do not want to add a ModelErr for this at the moment
            raise ValueError("Both x0 and x1 are None or both are set")

        x0 = np.asarray(x0)
        x1 = np.asarray(x1)
        self.__grid_is_integrated = False
        self.__grid_is_ascending = (x0[-1] > x0[0],
                                    x1[-1] > x1[0])
        self.__grid = (x0, x1)

    def __call__(self, model):
        """Evaluate a model on a different grid."""
        return RegridModel2D(model, self)

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
            format of the grid attribute of the model. At present this
            means a point-grid, so x0 and x1 arrays.
        kwargs
            Keyword arguments for the model.

        Notes
        -----
        If the requested grid (i.e. that defined by args) does not overlap
        the stored grid (the grid attribute) then all values are set to 0.
        However, if the grids partially overlap then there will be
        extrapolation (depending on the method).

        """

        nargs = len(args)
        if nargs == 0:
            raise ModelErr('nogrid')

        if self.grid is None:
            # There is no need for a regrid, so just evaluate
            # the model with the input grid
            #
            return modelfunc(pars, *args, **kwargs)

        if nargs < 2:
            raise ValueError("Expected a 2D grid")

        # For now we know that self.grid_is_integrated is False
        requested_grid_point = nargs == 2
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
        x0 = args[0]
        x1 = args[1]

        x01 = x0[0]
        x02 = x0[-1]
        if x02 > x01:
            x01, x02 = x02, x01

        x11 = x1[0]
        x12 = x1[-1]
        if x12 > x11:
            x11, x12 = x12, x11

        if requested_grid_point:
            g0, g1 = self.grid
        else:
            raise NotImplementedError("No support for 2D integrated grids")

        g01 = g0[0]
        g02 = g0[-1]
        if not self.grid_is_ascending[0]:
            g01, g02 = g02, g01

        g11 = g1[0]
        g12 = g1[-1]
        if not self.grid_is_ascending[1]:
            g11, g12 = g12, g11

        if ((g01 > x02) or (g02 < x01)) and ((g11 > x12) or (g12 < x11)):
            return np.zeros(x0.size, x1.size)

        # Evaluate the model on the user-defined grid and then interpolate
        # onto the desired grid. This is based on sherpa.models.TableModel
        # but is simplified as we do not provide a fold method.
        #
        #
        QUS: what does TableModel2D do? Do we have that? NOPE
        if requested_grid_point:
            y = modelfunc(pars, self.grid)
            return interpolate(args[0], self.grid, y,
                               function=self.method)

        else:
            # TODO: should there be some check that the grid size
            #       is "compatible"? Note that test_regrid1d_int_flux
            #       appears to fail if the grid width used for modelfunc
            #       is larger than the output.
            #
            y = modelfunc(pars, self.grid[0], self.grid[1])
            return rebin(y, self.grid[0], self.grid[1],
                         args[0], args[1])




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
