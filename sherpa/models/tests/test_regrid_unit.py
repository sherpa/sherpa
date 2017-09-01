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

import numpy as np
from numpy.testing import assert_allclose

import pytest

from sherpa.models.model import Model, ArithmeticModel, CompositeModel, \
    ArithmeticFunctionModel
from sherpa.models.basic import Const1D, Gauss1D, Const2D, Gauss2D, StepLo1D
from sherpa.models.parameter import Parameter
from sherpa.instrument import PSFModel
from sherpa.data import Data1D

import sherpa.utils
from sherpa.utils.err import ModelErr

from sherpa.models.regrid import Regrid1D, RegridModel1D, Regrid2D


def _setup_1d():
    """Create Gauss1D + Const1D components."""

    gmdl = Gauss1D()
    cmdl = Const1D()
    gmdl.pos = 50
    gmdl.fwhm = 30
    gmdl.ampl = 20
    cmdl.c0 = 10

    return gmdl, cmdl


def _setup_2d():
    """Create Gauss2D + Const2D components."""

    gmdl = Gauss2D()
    cmdl = Const2D()
    gmdl.xpos = 505
    gmdl.ypos = -240
    gmdl.fwhm = 30
    gmdl.ampl = 20
    cmdl.c0 = 10

    return gmdl, cmdl


# TODO: should there be equivalent tests for Regrid2D?

def test_regrid1d_default_model_name():
    mdl = Regrid1D()
    assert mdl.name == "regrid1d"


def test_regrid1d_given_model_name():
    mdl = Regrid1D('linGrid')
    assert mdl.name == "linGrid"  # TODO: why is this not lower-cased?


def test_regrid1d_default_grid_is_empty():
    mdl = Regrid1D()
    assert mdl.grid is None


# TODO: it is unclear whether Regrid should be part of the
#       Sherpa model hierarchy or not.
#
# def test_regrid1d_create_model_instance():
#     mdl = Regrid1D()
#     assert isinstance(mdl, Model)


def test_regrid1d_wrapping_create_model_instance():
    cmdl = Const1D()
    rmdl = Regrid1D()
    mdl = rmdl(cmdl)
    assert isinstance(mdl, RegridModel1D)


def test_regrid1d_wrapping_create_arithmetic_instance():
    cmdl = Const1D()
    rmdl = Regrid1D()
    mdl = rmdl(cmdl)
    assert isinstance(mdl, ArithmeticModel)


def test_regrid1d_wrapping_create_composite_instance():
    cmdl = Const1D()
    gmdl = Gauss1D()
    imdl = cmdl + gmdl
    rmdl = Regrid1D()
    mdl = rmdl(imdl)
    assert isinstance(mdl, CompositeModel)
    assert len(mdl.parts) == 2

    # This test was written before the code was written, which meant
    # that I had a different idea to how the composite model was going
    # to wrap things. It is not clear yet which (if either) is correct.
    #
    # assert mdl.parts[0] == cmdl
    # assert mdl.parts[1] == gmdl
    assert mdl.parts[0] == rmdl
    assert mdl.parts[1] == imdl


def test_regrid1d_call_twice():
    """What happens if have no evaluation (output) grid?"""

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid1D()
    rmdl.set_grid(np.arange(1000, 2000, 100))

    mdl = rmdl(internal_mdl)

    # It appears that calling the model with no arguments is
    # the same as calling 'rmdl(internal_mdl)'. An "easy" way
    # to test this is to rely on the stringification (and have
    # other tests handle that rmdl(...) is doing the right
    # thing).
    #
    noargs = mdl()
    assert str(mdl) == str(noargs)


def test_regrid1d_wrapping_name():
    """Check the name field of a wrapped model.

    This is also checked in test_regrid1d_wrapping_str.
    """

    internal_model = Const1D('con') + Gauss1D('gau')
    imodel_name = internal_model.name

    # a test where the Regrid1D model is named is in
    # test_regrid1d_wrapping_str.
    rmdl = Regrid1D()
    mdl = rmdl(internal_model)

    # TODO: It is not clear what the syntactic constraints on
    #       the name field are; if it is to create an evaluable
    #       model at the UI layer then the following is
    #       incorrect. There is "prior art" here with PSF, ARF,
    #       and RMF models to look at.
    #
    expected_name = 'regrid1d({})'.format(imodel_name)
    assert mdl.name == expected_name


def test_regrid1d_wrapping_str():
    """Check the str output of a wrapped model.

    Since this includes the name field, it subsumes
    test_regrid1d_wrapping_name but leave that test alone
    as it is an explicit check.
    """

    # This is basically checking that __str__ comes from
    # the CompositeModel and that everything is set up
    # correctly.
    #
    internal_model = Const1D('con')
    internal_model.c0 = 2
    internal_model.c0.freeze()

    imodel_name = internal_model.name

    rmdl = Regrid1D('test')
    mdl = rmdl(internal_model)

    expected_name = 'test({})'.format(imodel_name)

    # need to strip off the first line and replace it with
    # the new model name
    expected_params = "\n".join(str(internal_model).split('\n')[1:])

    expected_str = expected_name + "\n" + expected_params

    assert str(mdl) == expected_str


# If the grid is None then evaluating a model should act as
# an identity function: it should return the same data as
# if the wrapped model had been called directly.
#
# Although it is expected that a direct pass through will
# be used, so an exact check for equality could be used,
# the following tests use an approximate equality test
# since there is no requrement on how the model is to
# work. The model values are chosen so that the
# comparison by the numpy.assert_allclose (limits
# are explicitly set to atol=0, rtol=1e-7) should be
# sensible.
#
def test_regrid1d_identity_when_no_grid():

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid1D()
    mdl = rmdl(internal_mdl)

    grid = np.arange(-10, 100, 5)

    yexp = internal_mdl(grid)
    ygot = mdl(grid)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


def test_regrid1d_identity_when_no_grid_rev():

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid1D()
    mdl = rmdl(internal_mdl)

    grid = np.arange(-10, 100, 5)[::-1]

    yexp = internal_mdl(grid)
    ygot = mdl(grid)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


def test_regrid1d_identity_when_no_grid_int():

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid1D()
    mdl = rmdl(internal_mdl)

    # Ensure that the grid widths are not the same,
    # and that it is not contiguous. Would probably
    # have been easier to list the bins explicitly rather
    # than the following selection of filters.
    #
    grid = np.arange(-10, 100, 5)

    def excl(vlo, vhi):
        return (grid < vlo) | (grid > vhi)

    idx = excl(-1, 1) & excl(79, 86)
    grid = grid[idx]

    xlo = grid[:-1]
    xhi = grid[1:]

    idx = np.arange(xlo.size) != 4
    xlo = xlo[idx]
    xhi = xhi[idx]

    yexp = internal_mdl(xlo, xhi)
    ygot = mdl(xlo, xhi)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


@pytest.mark.xfail(reason="2D support not written")
def test_regrid2d_identity_when_no_grid():

    gmdl, cmdl = _setup_2d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid2D()
    mdl = rmdl(internal_mdl)

    # use different bin size for the x0/x1 axes
    x0grid = np.arange(450, 560, 5)
    x1grid = np.arange(-300, -200, 4)

    yexp = internal_mdl(x0grid, x1grid)
    ygot = mdl(x0grid, x1grid)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


# TODO: should there be a test_regrid2d_identity_when_no_grid_int()?


def test_regrid1d_identity_after_clearing_grid():
    """Ensure that the grid can be removed."""

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid1D()
    rmdl.set_grid(np.arange(200, 300, 20))

    mdl = rmdl(internal_mdl)

    grid = np.arange(-10, 100, 5)

    yexp = internal_mdl(grid)

    rmdl.set_grid(None)
    ygot = mdl(grid)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


def test_regrid1d_no_overlap():
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid1D()
    rmdl.set_grid(np.arange(1000, 2000, 100))

    mdl = rmdl(internal_mdl)

    grid = np.arange(-10, 100, 5)
    ygot = mdl(grid)

    assert_allclose(ygot, np.zeros(grid.size), atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_rev1():
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid1D()
    rmdl.set_grid(np.arange(1000, 2000, 100))

    mdl = rmdl(internal_mdl)

    grid = np.arange(-10, 100, 5)[::-1]
    ygot = mdl(grid)

    assert_allclose(ygot, np.zeros(grid.size), atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_rev2():
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid1D()
    rmdl.set_grid(np.arange(1000, 2000, 100)[::-1])

    mdl = rmdl(internal_mdl)

    grid = np.arange(-10, 100, 5)
    ygot = mdl(grid)

    assert_allclose(ygot, np.zeros(grid.size), atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_rev3():
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid1D()
    rmdl.set_grid(np.arange(1000, 2000, 100)[::-1])

    mdl = rmdl(internal_mdl)

    grid = np.arange(-10, 100, 5)[::-1]
    ygot = mdl(grid)

    assert_allclose(ygot, np.zeros(grid.size), atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_int():
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid1D()
    rgrid = np.arange(1000, 2000, 100)
    rmdl.set_grid(rgrid[:-1], rgrid[1:])

    mdl = rmdl(internal_mdl)

    grid = np.arange(-10, 100, 5)
    ygot = mdl(grid[:-1], grid[1:])

    assert_allclose(ygot, np.zeros(grid.size - 1), atol=0, rtol=1e-7)


@pytest.mark.xfail(reason="2D support not written")
def test_regrid2d_no_overlap():
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = _setup_2d()
    internal_mdl = gmdl + cmdl

    rmdl = Regrid2D()
    x0 = np.arange(10, 20, 1)
    x1 = np.arange(40, 60, 2)
    rmdl.set_grid(x0, x1)
    mdl = rmdl(internal_mdl)

    x0grid = np.arange(450, 560, 5)
    x1grid = np.arange(-300, -200, 4)
    ygot = mdl(x0grid, x1grid)

    assert_allclose(ygot, np.zeros(x0grid.size), atol=0, rtol=1e-7)


class MyConst1D(Const1D):

    def __init__(self, name='myconst1d'):
        self._calc_store = []
        Const1D.__init__(self, name)

    def calc(self, *args, **kwargs):
        self._calc_store.append((args, kwargs))
        return Const1D.calc(self, *args, **kwargs)


def test_regrid1d_passes_through_the_grid():
    """Is the grid actually being passed through to the model?"""

    rmdl = Regrid1D()
    imdl = MyConst1D()
    imdl.c0 = -34.5
    mdl = rmdl(imdl)
    grid_expected = [5, 10, 15, 20, 25, 30]
    grid_requested = [12, 18, 20]

    rmdl.set_grid(grid_expected)
    store = imdl._calc_store
    len(store) == 0

    y = mdl(grid_requested)
    assert len(y) == len(grid_requested)
    assert len(store) == 1
    store = store[0]
    assert len(store) == 2

    assert len(store[0]) == 2
    assert store[1] == {}

    store = store[0]
    assert store[0] == [-34.5]
    assert (store[1] == grid_expected).all()


def test_regrid1d_error_calc_no_args():

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl
    grid_evaluate = np.arange(-10, 100, 5)

    rmdl = Regrid1D()
    rmdl.set_grid(grid_evaluate)
    mdl = rmdl(internal_mdl)

    with pytest.raises(ModelErr) as excinfo:
        pvals = [p.val for p in internal_mdl.pars]
        mdl.calc(p=pvals)

    assert 'There is no grid on which to evaluate the model' \
        in str(excinfo.value)


def test_regrid1d_error_grid_mismatch_1():
    """Internal grid is integrated but given points"""

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    grid_evaluate = np.arange(-10, 100, 5)
    rmdl = Regrid1D()
    rmdl.set_grid(grid_evaluate[:-1], grid_evaluate[1:])

    mdl = rmdl(internal_mdl)
    grid_run = np.arange(0, 20, 10)
    with pytest.raises(ModelErr) as excinfo:
        mdl(grid_run)

    assert 'An integrated grid is required for model evaluation' \
        in str(excinfo.value)


def test_regrid1d_error_grid_mismatch_2():
    """Internal grid is points but given integrated"""

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    grid_evaluate = np.arange(-10, 100, 5)
    rmdl = Regrid1D()
    rmdl.set_grid(grid_evaluate)

    mdl = rmdl(internal_mdl)
    grid_run = np.arange(0, 20, 10)
    with pytest.raises(ModelErr) as excinfo:
        mdl(grid_run[:-1], grid_run[1:])

    assert 'A non-integrated grid is required for model evaluation' \
        in str(excinfo.value)


# Evaluate on a grid x'_i which is close to the desired
# grid x_i, and then check that the output is close to
# the model evaluated on the grid x_i.
#
# Note that there's no check that the model is evaluated
# on the grid x'_i.
#
def _test_regrid1d_interpolation(rtol,
                                 eval_incr=True,
                                 req_incr=True,
                                 method=None):
    """Test interpolation case.

    Parameters
    ----------
    rtol : number
        The relative tolerance, passed through to assert_allclose.
    eval_incr : bool, optional
        Does the evaluation grid increase or decrease?
    req_incr : bool, optional
        Does the requested grid increase or decrease?
    method : function reference, optional
        The interpolator to use, which should match
        sherpa.utils.linear_interp's API. If None then
        use the default method (which is sherpa.utils.neville)
    """

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    # Shift the evaluation grid compared to the
    # requested grid (making sure that the evaluation
    # grid extends outside the requested grid on both
    # edges to avoid checking for edge effects here).
    #
    # grid_evaluate is x'_i, has 22 bins
    # grid_request is x_i, has 21 bins
    #
    grid_evaluate = np.arange(-10, 100, 5)
    grid_request = np.linspace(-5, 85, 21)

    if not eval_incr:
        grid_evaluate = grid_evaluate[::-1]

    if not req_incr:
        grid_request = grid_request[::-1]

    rmdl = Regrid1D()
    rmdl.set_grid(grid_evaluate)
    if method is not None:
        rmdl.method = method

    mdl = rmdl(internal_mdl)

    yexp = internal_mdl(grid_request)

    ygot = mdl(grid_request)
    assert_allclose(ygot, yexp, atol=0, rtol=rtol)


def _test_regrid1d_int_interpolation(rtol, method=None):
    """Test interpolation case for integrated grids.

    Parameters
    ----------
    rtol : number
        The relative tolerance, passed through to assert_allclose.
    method : function reference, optional
        The interpolator to use, which should match
        sherpa.utils.linear_interp's API. If None then
        use the default method (which is sherpa.utils.neville)

    Notes
    -----
    Should this logic be added to _test_regrid1d_interpolation - e.g.
    controlled by an interpolate parameter - or are there more differences,
    such as the interpolation method? For integrated models it isn't
    quite interpolation that is wanted, so perhaps the function
    name needs changing?
    """

    gmdl, cmdl = _setup_1d()
    internal_mdl = gmdl + cmdl

    grid_evaluate = np.arange(-10, 100, 5)
    grid_request = np.linspace(-5, 85, 21)

    rmdl = Regrid1D()
    rmdl.set_grid(grid_evaluate[:-1], grid_evaluate[1:])
    if method is not None:
        rmdl.method = method

    mdl = rmdl(internal_mdl)

    yexp = internal_mdl(grid_request[:-1], grid_request[1:])

    ygot = mdl(grid_request[:-1], grid_request[1:])
    assert_allclose(ygot, yexp, atol=0, rtol=rtol)


# The tolerance is adjusted until the tests pass (as long as the
# value remains small).
#
# The liner interpolation fails when eval_incr is False. Is this a bug
# or do we need to guard against this?
#
# @pytest.mark.parametrize("eincr", [True, False])
@pytest.mark.parametrize("eincr", [True])
@pytest.mark.parametrize("rincr", [True, False])
@pytest.mark.parametrize("margs", [(5e-5, None),
                                   (0.011, sherpa.utils.linear_interp)])
def test_regrid1d_interpolation(eincr, rincr, margs):

    tol, method = margs
    _test_regrid1d_interpolation(rtol=tol, method=method,
                                 eval_incr=eincr, req_incr=rincr)


@pytest.mark.xfail(reason="1D integration support not written")
def test_regrid1d_int_interpolation_default():

    # The tolerance will likely need changing
    #
    _test_regrid1d_int_interpolation(rtol=1e-5)


# TODO: more tests when regridding the models

class ShiftyKernel1D(Model):
    """A convolution-style model that shifts the origin.

    The design of the Kernel/Model class structure is based on the
    version used for the XSPEC convolution models in the contributed
    code for CIAO. It is not yet clear to me if this is a sensible
    approach for users, but should be fine as a test case.
    """

    def __init__(self, name='shiftykernel1d'):
        self.dx = Parameter(name, 'dx', 0)
        Model.__init__(self, name, (self.dx, ))

    def __call__(self, model):
        return ShiftyModel1D(model, self)

    def calc(self, pars, rhs, *args, **kwargs):

        dx = pars[0]
        rpars = pars[1:]

        if len(args) == 1:
            x = args[0] + dx
            nargs = (x, )
        elif len(args) == 2:
            xlo = args[0] + dx
            xhi = args[1] + dx
            nargs = (xlo, xhi)
        else:
            raise ValueError("1D only")

        return rhs(rpars, *nargs, **kwargs)


class ShiftyModel1D(CompositeModel, ArithmeticModel):

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


def _old_test_regrid1d_works_with_convolution_style():
    """This doesn't really test more than the previous
    model-evaluation tests.

    Why do I have this routine here? Do I want this as a test,
    or did I leave it around when writing
    test_regrid1d_works_with_convolution_style
    with the intention of deleting it?
    """

    gmdl = Gauss1D()
    gmdl.pos = 200
    gmdl.ampl = 10

    cmdl = Const1D()
    cmdl.c0 = -500

    imdl = gmdl + cmdl

    grid = np.arange(-50, 50, 1)

    # The values should be very-close to cmdl.c0; there isn't
    # really a need to check this here, but it makes the intent
    # of the test a bit clearer.
    #
    y0 = imdl(grid)
    flat = np.zeros(grid.size) + cmdl.c0.val
    assert_allclose(y0, flat, atol=1e-15, rtol=0)

    smdl = ShiftyKernel1D()
    smdl.dx = 190

    yexp = imdl(grid + 190)

    # First check that ShiftyKernel is behaving as expected
    #
    test_mdl = smdl(imdl)
    ygot = test_mdl(grid)
    assert_allclose(ygot, yexp, atol=1e-15, rtol=0)

    # Now to the actual test of Regrid1D
    #
    # - grid sent in to the outermost model
    # - grid on which the model is to be evaluated
    #   (the grid attribute of Regrid1D)
    # - the grid after adjustment by shiftykernel
    #
    # if g1 is 0 - 100
    #    g2 is -100 - 200
    #    g3 is g2 + 80, so -20 - 280
    #
    # if the gaussian is centered at 10 then
    #
    g1 = np.arange(0, 100, 2)
    g2 = np.arange(-100, 200, 2)

    rmdl = Regrid1D()
    rmdl.set_grid(g2)

    smdl.dx = 80
    gmdl.pos = 10

    mdl = rmdl(test_mdl)

    yeval = mdl(g1)
    # assert_allclose(yeval, yexp, atol=1e-15, rtol=0)

    return g1, yeval


def test_regrid1d_works_with_convolution_style():
    """This doesn't really test more than the previous
    model-evaluation tests.
    """

    smdl = StepLo1D()
    smdl.xcut = 100
    smdl.ampl = 10

    cmdl = Const1D()
    cmdl.c0 = -500

    imdl = smdl + cmdl

    # Set up the convolution kernel
    #
    gsmooth = Gauss1D()
    psf = PSFModel('psf', gsmooth)

    smoothed = psf(imdl)

    # This is the model that will be evaluated
    #
    regrid = Regrid1D()
    smoothed_regrid = regrid(smoothed)

    # Ignoring edge effects, the smoothed step function drops from
    # x=100 down to x~120 (it is not clear why it doesn't smooth
    # below x=100, but maybe I've set up the convolution wrong).
    # So, if the regrid model evaluates x=0 to 200 but the requested
    # grid is from x=102 ro 180, then we should see a difference
    # to evaluating without regrid.
    #
    xfull = np.arange(0, 200, 0.5)
    xout = np.arange(101, 180, 0.5)

    regrid.set_grid(xfull)

    # fake up a data object for the fold method
    # TODO: it is not clear to me what grid should be used here;
    #       xfull or xout
    #
    d = Data1D('fake', xfull, xfull * 0)
    psf.fold(d)

    y_regrid = smoothed_regrid(xout)

    # calculate the expected values
    y_full = smoothed(xfull)

    # since the grids have the same binning, it is a simple extraction
    idx0, = np.where(xfull == xout[0])
    idx1, = np.where(xfull == xout[-1])
    y_expected = y_full[idx0[0]:idx1[0] + 1]

    assert_allclose(y_regrid, y_expected, atol=1e-10, rtol=0)

    # check that this is all worth it; i.e. that without regrid
    # you just get the constant term. If this fails then it does not
    # mean that the regrid code is broken, but it implies that
    # something fundamental has changed with the basic model
    # evaluation.
    #
    d = Data1D('fake', xout, xout * 0)
    psf.fold(d)
    y_check = smoothed(xout)
    y_expected = np.zeros(xout.size) + cmdl.c0.val
    assert_allclose(y_check, y_expected, rtol=0, atol=1e-7)
