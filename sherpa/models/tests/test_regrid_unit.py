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

import numpy as np
from numpy.testing import assert_allclose

import pytest

from sherpa.models.model import Model, ArithmeticModel, CompositeModel, \
    ArithmeticFunctionModel
from sherpa.models.basic import Const1D, Gauss1D, Const2D, Gauss2D, \
    PowLaw1D, StepLo1D
from sherpa.models.parameter import Parameter
from sherpa.instrument import PSFModel
from sherpa.data import Data1D

import sherpa.utils
from sherpa.utils.err import ModelErr

from sherpa.models.regrid import ModelDomainRegridder1D, RegridModel1D, EvaluationSpace1D


@pytest.fixture
def setup_1d():
    """Create Gauss1D + Const1D components."""

    gmdl = Gauss1D()
    cmdl = Const1D()
    gmdl.pos = 50
    gmdl.fwhm = 30
    gmdl.ampl = 20
    cmdl.c0 = 10

    return gmdl, cmdl


@pytest.mark.parametrize("cls,name",
                         [(ModelDomainRegridder1D, "regrid1d"),
                          ])
def test_default_model_name(cls, name):
    mdl = cls()
    assert mdl.name == name


@pytest.mark.parametrize("cls", [ModelDomainRegridder1D])
def test_given_model_name(cls):
    mdl = cls(name='linGrid')
    assert mdl.name == "linGrid"  # TODO: why is this not lower-cased?


@pytest.mark.parametrize("regrid_class,model_class,regrid_model_class",
                         [(ModelDomainRegridder1D, Const1D, RegridModel1D)])
def test_wrapping_create_model_instance(regrid_class,
                                        model_class,
                                        regrid_model_class):
    rmdl = regrid_class()
    cmdl = model_class()
    mdl = rmdl.apply_to(cmdl)
    assert isinstance(mdl, regrid_model_class)


@pytest.mark.parametrize("regrid_class,model_class",
                         [(ModelDomainRegridder1D, Const1D)])
def test_wrapping_create_arithmetic_instance(regrid_class, model_class):
    rmdl = regrid_class()
    cmdl = model_class()
    mdl = rmdl.apply_to(cmdl)
    assert isinstance(mdl, ArithmeticModel)


def test_regrid1d_wrapping_create_composite_instance():
    # This test depends on what we want the regridded model to look like, which is
    # somewhat arbitrary
    cmdl = Const1D()
    gmdl = Gauss1D()
    imdl = cmdl + gmdl
    rmdl = ModelDomainRegridder1D()
    mdl = rmdl.apply_to(imdl)
    assert isinstance(mdl, CompositeModel)
    assert len(mdl.parts) == 1
    assert mdl.parts[0] is imdl


def test_regrid1d_call_twice(setup_1d):
    """What happens if have no evaluation (output) grid?"""

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    eval_space = EvaluationSpace1D(np.arange(1000, 2000, 100))
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

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
    rmdl = ModelDomainRegridder1D()
    mdl = rmdl.apply_to(internal_model)

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

    rmdl = ModelDomainRegridder1D(name='test')
    mdl = rmdl.apply_to(internal_model)

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
def test_regrid1d_identity_when_no_grid(setup_1d):

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    rmdl = ModelDomainRegridder1D()
    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)

    yexp = internal_mdl(grid)
    ygot = mdl(grid)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


def test_regrid1d_identity_when_no_grid_rev(setup_1d):

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    rmdl = ModelDomainRegridder1D()
    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)[::-1]

    yexp = internal_mdl(grid)
    ygot = mdl(grid)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


def test_regrid1d_identity_when_no_grid_int(setup_1d):

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    rmdl = ModelDomainRegridder1D()
    mdl = rmdl.apply_to(internal_mdl)

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
# TODO: should there be a test_regrid2d_identity_when_no_grid_int()?


def test_regrid1d_identity_after_clearing_grid(setup_1d):
    """Ensure that the grid can be removed."""

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    eval_space = EvaluationSpace1D(np.arange(200, 300, 20))

    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)

    yexp = internal_mdl(grid)

    rmdl.grid = None
    ygot = mdl(grid)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


def test_regrid1d_no_overlap(setup_1d):
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    eval_space = EvaluationSpace1D(np.arange(1000, 2000, 100))
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)
    ygot = mdl(grid)

    assert_allclose(ygot, np.zeros(grid.size), atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_rev1(setup_1d):
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    eval_space = EvaluationSpace1D(np.arange(1000, 2000, 100))
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)[::-1]
    ygot = mdl(grid)

    assert_allclose(ygot, np.zeros(grid.size), atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_rev2(setup_1d):
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    eval_space = EvaluationSpace1D(np.arange(1000, 2000, 100)[::-1])
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)
    ygot = mdl(grid)

    assert_allclose(ygot, np.zeros(grid.size), atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_rev3(setup_1d):
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    eval_space = EvaluationSpace1D(np.arange(1000, 2000, 100)[::-1])
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)[::-1]
    ygot = mdl(grid)

    assert_allclose(ygot, np.zeros(grid.size), atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_int(setup_1d):
    """If the two grids have no overlap, return value is 0."""

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    array = np.arange(1000, 2000, 100)
    eval_space = EvaluationSpace1D(array[:-1], array[1:])
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)
    ygot = mdl(grid[:-1], grid[1:])

    assert_allclose(ygot, np.zeros(grid.size - 1), atol=0, rtol=1e-7)


class MyConst1D(Const1D):

    def __init__(self, name='myconst1d'):
        self._calc_store = []
        Const1D.__init__(self, name)

    def calc(self, *args, **kwargs):
        self._calc_store.append((args, kwargs))
        return Const1D.calc(self, *args, **kwargs)


def test_regrid1d_passes_through_the_grid():
    """Is the grid actually being passed through to the model?"""

    rmdl = ModelDomainRegridder1D()
    imdl = MyConst1D()
    imdl.c0 = -34.5
    mdl = rmdl.apply_to(imdl)
    grid_expected = [5, 10, 15, 20, 25, 30]
    grid_requested = [12, 18, 20]

    rmdl.grid = grid_expected
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


def test_regrid1d_error_calc_no_args(setup_1d):

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl
    grid_evaluate = EvaluationSpace1D(np.arange(-10, 100, 5))

    rmdl = ModelDomainRegridder1D(grid_evaluate)
    mdl = rmdl.apply_to(internal_mdl)

    with pytest.raises(ModelErr) as excinfo:
        pvals = [p.val for p in internal_mdl.pars]
        mdl.calc(p=pvals)

    assert 'There is no grid on which to evaluate the model' \
        in str(excinfo.value)


def test_regrid1d_error_grid_mismatch_1(setup_1d):
    """Internal grid is integrated but given points"""

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    grid_evaluate = np.arange(-10, 100, 5)
    eval_space = EvaluationSpace1D(grid_evaluate[:-1], grid_evaluate[1:])
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)
    grid_run = np.arange(0, 20, 10)
    with pytest.raises(ModelErr) as excinfo:
        mdl(grid_run)

    assert 'An integrated grid is required for model evaluation' \
        in str(excinfo.value)


def test_regrid1d_error_grid_mismatch_2(setup_1d):
    """Internal grid is points but given integrated"""

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    grid_evaluate = EvaluationSpace1D(np.arange(-10, 100, 5))
    rmdl = ModelDomainRegridder1D(grid_evaluate)

    mdl = rmdl.apply_to(internal_mdl)
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
                                 method=None,
                                 setup_1d=None):
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

    gmdl, cmdl = setup_1d
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

    rmdl = ModelDomainRegridder1D(EvaluationSpace1D(grid_evaluate))
    if method is not None:
        rmdl.method = method

    mdl = rmdl.apply_to(internal_mdl)

    yexp = internal_mdl(grid_request)

    ygot = mdl(grid_request)
    assert_allclose(ygot, yexp, atol=0, rtol=rtol)


def _test_regrid1d_int(rtol,
                       eval_incr=True,
                       req_incr=True,
                       setup_1d=None):
    """Test with for integrated grids.

    Parameters
    ----------
    rtol : number
        The relative tolerance, passed through to assert_allclose.
    eval_incr : bool, optional
        Does the evaluation grid increase or decrease?
        CURRENTLY UNSUPPORTED IF False.
    req_incr : bool, optional
        Does the requested grid increase or decrease?
        CURRENTLY UNSUPPORTED IF False.

    Notes
    -----
    This is very similar to  _test_regrid1d_interpolation except that
    it does not support the method parameter.
    """

    # have not coded this, since need xlo[0] > xlo[1] but xhi[0] > xlo[0]
    # (I think).
    assert eval_incr
    assert req_incr

    gmdl, cmdl = setup_1d
    internal_mdl = gmdl + cmdl

    grid_evaluate = np.arange(-10, 100, 5)
    grid_request = np.linspace(-5, 85, 21)

    eval_space = EvaluationSpace1D(grid_evaluate[:-1], grid_evaluate[1:])

    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

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
def test_regrid1d_interpolation(eincr, rincr, margs, setup_1d):

    tol, method = margs
    _test_regrid1d_interpolation(rtol=tol, method=method,
                                 eval_incr=eincr, req_incr=rincr, setup_1d=setup_1d)


def test_regrid1d_int(setup_1d):
    _test_regrid1d_int(rtol=0.015, setup_1d=setup_1d)


# Can use a "calculate the flux" style model (e.g. XSPEC's c[p]flux)
# to test out the 1D integrated case.
#
class ReNormalizerKernel1DInt(Model):
    """A convolution-style model which renormalizes supplied model.

    The signal between the lo and hi values is forced to equal
    the flux parameter.
    """
    def __init__(self, name='renormalizerkernel1d'):
        self.flux = Parameter(name, 'flux', 1.0)
        self.lo = Parameter(name, 'lo', 0, alwaysfrozen=True)
        self.hi = Parameter(name, 'hi', 100, alwaysfrozen=True)
        Model.__init__(self, name, (self.flux, self.lo, self.hi))

    def __call__(self, model):
        return ReNormalizerModel1DInt(model, self)

    def calc(self, pars, rhs, *args, **kwargs):

        flux, lo, hi = pars[0:3]
        rpars = pars[3:]

        if len(args) == 2:
            xlo = args[0]
            xhi = args[1]
            nargs = (xlo, xhi)
        else:
            raise ValueError("1D Int only")

        # In a real-world version edge effects would be an issue,
        # but here just worry about whether a bin is in or out.
        # Can also just worry about the grids used in the test!
        #
        # Assumes the grid is increasing, and that idx is not empty.
        #
        idx = (xhi > lo) & (xlo < hi)
        if not idx.any():
            return np.zeros(xlo.size)

        yorig = rhs(rpars, *nargs, **kwargs)
        yflux = yorig[idx].sum()
        return yorig * flux / yflux


class ReNormalizerModel1DInt(CompositeModel, ArithmeticModel):

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


# TODO: more tests when regridding the models

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
    regrid = ModelDomainRegridder1D()
    smoothed_regrid = regrid.apply_to(smoothed)

    # Ignoring edge effects, the smoothed step function drops from
    # x=100 down to x~120 (it is not clear why it doesn't smooth
    # below x=100, but maybe I've set up the convolution wrong).
    # So, if the regrid model evaluates x=0 to 200 but the requested
    # grid is from x=102 ro 180, then we should see a difference
    # to evaluating without regrid.
    #
    xfull = np.arange(0, 200, 0.5)
    xout = np.arange(101, 180, 0.5)

    regrid.grid = xfull

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


def test_regrid1d_int_flux():
    """Check integration models using an XSPEC-style c[p]flux model.
    """

    gamma = 1.7
    pmdl = PowLaw1D()
    pmdl.gamma = gamma

    fluxmdl = ReNormalizerKernel1DInt()
    fluxmdl.flux = 20.0
    fluxmdl.lo = 10.0
    fluxmdl.hi = 20

    grid = np.arange(1.0, 9.0, 0.01)
    glo = grid[:-1]
    ghi = grid[1:]

    mdl = fluxmdl(pmdl)

    # Without any regridding, the model should be zero
    yzero = mdl(glo, ghi)
    assert_allclose(yzero, np.zeros(glo.size), atol=1e-10, rtol=0)

    regrid = ModelDomainRegridder1D()
    rmdl = regrid.apply_to(mdl)

    # ensure it covers the 10 - 20 range as well as 1-9. Pick
    # a smaller grid size than the output grid.
    #
    xfull = np.arange(0, 22, 0.005)
    regrid.grid = xfull[:-1], xfull[1:]

    y_regrid = rmdl(glo, ghi)
    assert y_regrid.shape == glo.shape

    # Model flux in lo/hi range is int_lo^hi x^-gamma dx (since norm
    # and reference of the power law are both 1.0). This is, since
    # gamma is not 1, [x^(1-gamma)]^hi_lo / (1 - gamma).
    #
    term = 1.0 - gamma
    renorm = (20**term - 10**term) / term
    y_expected = pmdl(glo, ghi) * 20 / renorm

    assert_allclose(y_regrid, y_expected, atol=1e-10, rtol=0)
