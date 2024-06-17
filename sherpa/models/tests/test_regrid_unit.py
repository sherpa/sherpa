#
#  Copyright (C) 2017 - 2024
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

import re

import numpy as np
from numpy.testing import assert_allclose

import pytest

from sherpa.models.model import Model, ArithmeticModel, BinaryOpModel, \
    CompositeModel, ArithmeticConstantModel, ArithmeticFunctionModel, \
    RegriddableModel1D, RegridWrappedModel, UnaryOpModel
from sherpa.models.basic import Box1D, Const1D, Gauss1D, \
    PowLaw1D, StepLo1D
from sherpa.models.parameter import Parameter
from sherpa.instrument import PSFModel
from sherpa.data import Data1D, Data1DInt
from sherpa.astro import ui

from sherpa.utils import linear_interp
from sherpa.utils.err import DataErr, ModelErr
from sherpa.utils.numeric_types import SherpaFloat

from sherpa.models.regrid import ModelDomainRegridder1D, EvaluationSpace1D, \
    EvaluationSpace2D, PointAxis, IntegratedAxis


@pytest.fixture(params=[True, False])
def setup_1d(request):
    """Create Gauss1D + Const1D components."""

    gmdl = Gauss1D()
    cmdl = Const1D()
    gmdl.pos = 5000
    gmdl.fwhm = 30
    gmdl.ampl = 20
    cmdl.c0 = 10

    return_composite = request.param

    if return_composite:
        return gmdl + cmdl

    return cmdl


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
                         [(ModelDomainRegridder1D, Const1D, RegridWrappedModel)])
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

    internal_mdl = setup_1d

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
    expected_name = f'regrid1d({imodel_name})'
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

    expected_name = f'test({imodel_name})'

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
# since there is no requirement on how the model is to
# work. The model values are chosen so that the
# comparison by the numpy.assert_allclose (limits
# are explicitly set to atol=0, rtol=1e-7) should be
# sensible.
#
def test_regrid1d_identity_when_no_grid(setup_1d):

    internal_mdl = setup_1d

    rmdl = ModelDomainRegridder1D()
    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)

    yexp = internal_mdl(grid)
    ygot = mdl(grid)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


def test_regrid1d_identity_when_no_grid_rev(setup_1d):

    internal_mdl = setup_1d

    rmdl = ModelDomainRegridder1D()
    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)[::-1]

    yexp = internal_mdl(grid)
    ygot = mdl(grid)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


def test_regrid1d_identity_when_no_grid_int(setup_1d):

    internal_mdl = setup_1d

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


def test_regrid1d_identity_after_clearing_grid(setup_1d):
    """Ensure that the grid can be removed."""

    internal_mdl = setup_1d

    eval_space = EvaluationSpace1D(np.arange(200, 300, 20))

    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)

    yexp = internal_mdl(grid)

    rmdl.grid = None
    ygot = mdl(grid)

    assert_allclose(ygot, yexp, atol=0, rtol=1e-7)


def test_regrid1d_no_overlap(setup_1d):
    """If the two grids have no overlap, return value is the same as the model evaluated over the data space."""

    internal_mdl = setup_1d

    eval_space = EvaluationSpace1D(np.arange(1000, 2000, 100))
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)

    ygot = mdl(grid)

    assert_allclose(ygot, [0., ]*grid.size, atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_rev1(setup_1d):
    """If the two grids have no overlap, return value is the same as the model evaluated over the data space."""

    internal_mdl = setup_1d

    eval_space = EvaluationSpace1D(np.arange(1000, 2000, 100))
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)[::-1]

    ygot = mdl(grid)

    assert_allclose(ygot, [0., ]*grid.size, atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_rev2(setup_1d):
    """If the two grids have no overlap, return value is the same as the model evaluated over the data space."""

    internal_mdl = setup_1d

    eval_space = EvaluationSpace1D(np.arange(1000, 2000, 100)[::-1])
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)

    ygot = mdl(grid)

    assert_allclose(ygot, [0., ]*grid.size, atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_rev3(setup_1d):
    """If the two grids have no overlap, return value is the same as the model evaluated over the data space."""

    internal_mdl = setup_1d

    eval_space = EvaluationSpace1D(np.arange(1000, 2000, 100)[::-1])
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)[::-1]

    ygot = mdl(grid)

    assert_allclose(ygot, [0., ]*grid.size, atol=0, rtol=1e-7)


def test_regrid1d_no_overlap_int(setup_1d):
    """If the two grids have no overlap, return value is the same as the model evaluated over the data space."""

    internal_mdl = setup_1d

    array = np.arange(1000, 2000, 100)
    eval_space = EvaluationSpace1D(array[:-1], array[1:])
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)

    grid = np.arange(-10, 100, 5)

    ygot = mdl(grid[:-1], grid[1:])

    assert_allclose(ygot, [0., ]*(grid.size - 1), atol=0, rtol=1e-7)


class MyConst1D(Const1D):

    def __init__(self, name='myconst1d'):
        self._calc_store = []
        Const1D.__init__(self, name)

    def calc(self, p, *args, **kwargs):
        self._calc_store.append((p, args, kwargs))
        return Const1D.calc(self, p, *args, **kwargs)


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
    assert len(store) == 0

    y = mdl(grid_requested)
    assert len(y) == len(grid_requested)
    assert len(store) == 1

    p, args, kwargs = store[0]

    assert p == [-34.5]

    combine = np.unique(np.append(grid_expected, grid_requested))
    indices = combine.searchsorted(args)
    assert (args == combine[indices]).all()

    assert kwargs == {'integrate': True}


def test_regrid1d_error_calc_no_args(setup_1d):

    internal_mdl = setup_1d
    grid_evaluate = EvaluationSpace1D(np.arange(-10, 100, 5))

    rmdl = ModelDomainRegridder1D(grid_evaluate)
    mdl = rmdl.apply_to(internal_mdl)

    with pytest.raises(ModelErr,
                       match=ModelErr.dict['nogrid']):
        pvals = [p.val for p in internal_mdl.pars]
        mdl.calc(p=pvals)


def test_regrid1d_error_grid_mismatch_1(setup_1d):
    """Internal grid is integrated but given points"""

    internal_mdl = setup_1d

    grid_evaluate = np.arange(-10, 100, 5)
    eval_space = EvaluationSpace1D(grid_evaluate[:-1], grid_evaluate[1:])
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)
    grid_run = np.arange(0, 20, 10)
    with pytest.raises(ModelErr,
                       match=re.escape(ModelErr.dict['needsint'])):
        mdl(grid_run)


def test_ui_regrid1d_non_overlapping_not_allowed():
    """Integrated data space must not overlap"""

    ui.dataspace1d(1, 100, 2, dstype=Data1DInt)
    b1 = Box1D()
    ui.set_model(b1)
    b1.xlow = 10
    b1.xhi = 80
    b1.ampl.max = 100
    grid_hi = np.linspace(2, 101, 600)
    grid_lo = np.linspace(1, 100, 600)
    with pytest.raises(ModelErr,
                       match=re.escape(ModelErr.dict['needsint'])):
        b1.regrid(grid_lo, grid_hi)


def test_low_level_regrid1d_non_overlapping_not_allowed():
    """Integrated data space must not overlap"""

    c = Box1D()
    lo = np.linspace(1, 100, 600)
    hi = np.linspace(2, 101, 600)
    with pytest.raises(ModelErr,
                       match=re.escape(ModelErr.dict['needsint'])):
        c.regrid(lo, hi)


def test_regrid1d_error_grid_mismatch_2(setup_1d):
    """Internal grid is points but given integrated"""

    internal_mdl = setup_1d

    grid_evaluate = EvaluationSpace1D(np.arange(-10, 100, 5))
    rmdl = ModelDomainRegridder1D(grid_evaluate)

    mdl = rmdl.apply_to(internal_mdl)
    grid_run = np.arange(0, 20, 10)
    with pytest.raises(ModelErr,
                       match=ModelErr.dict['needspoint']):
        mdl(grid_run[:-1], grid_run[1:])


@pytest.mark.parametrize("requested",
                         [np.arange(1, 7, 0.1),
                          np.arange(1, 7, 0.05),
                          np.arange(1, 7, 0.2)])
def test_low_level_regrid1d_full_overlap(requested):
    """Base case of test_low_level_regrid1d_partial_overlap
    """

    # The range over which we want the model evaluated
    xgrid = np.arange(2, 6, 0.1)
    d = Data1D('tst', xgrid, np.ones_like(xgrid))

    mdl = Box1D()
    mdl.xlow = 3.1
    mdl.xhi = 4.2
    mdl.ampl = 0.4

    yexpected = d.eval_model(mdl)
    assert yexpected.min() == pytest.approx(0.0)
    assert yexpected.max() == pytest.approx(0.4)

    ygot = d.eval_model(mdl.regrid(requested))
    assert ygot == pytest.approx(yexpected)


@pytest.mark.parametrize("requested",
                         [np.arange(2.5, 7, 0.2),
                          np.arange(1, 5.1, 0.2),
                          np.arange(2.5, 5.1, 0.2),
                          np.arange(2.5, 7, 0.075),
                          np.arange(1, 5.1, 0.075),
                          np.arange(2.5, 5.1, 0.075)])
def test_low_level_regrid1d_partial_overlap(requested):
    """What happens if there is partial overlap of the grid?

    The question becomes how do we evaluate the model
    "outside" the regrid range. There's at least two
    options:

      a) set to 0
      b) use the original grid

    This test is chosen so that it holds with both
    possibilities: the model evaluates to 0 outside
    of x=3.1 - 4.2, and the partial overlaps are
    carefully chosen to always include this full
    range.

    See https://github.com/sherpa/sherpa/issues/722
    """

    # The range over which we want the model evaluated
    xgrid = np.arange(2, 6, 0.1)
    d = Data1D('tst', xgrid, np.ones_like(xgrid))

    mdl = Box1D()
    mdl.xlow = 3.1
    mdl.xhi = 4.2
    mdl.ampl = 0.4

    yexpected = d.eval_model(mdl)
    assert yexpected.min() == pytest.approx(0.0)
    assert yexpected.max() == pytest.approx(0.4)

    ygot = d.eval_model(mdl.regrid(requested))
    assert ygot == pytest.approx(yexpected)


@pytest.mark.parametrize("requested, tol",
                         [(np.arange(1, 7, 0.1), 1e-7),
                          (np.arange(1, 7, 0.05), 1e-7),
                          (np.arange(1, 7, 0.2), 0.02)])
def test_low_level_regrid1d_int_full_overlap(requested, tol):
    """Base case of test_low_level_regrid1d_int_partial_overlap
    """

    # The range over which we want the model evaluated
    dx = 0.1
    xgrid = np.arange(2, 6, dx)
    xlo = xgrid[:-1]
    xhi = xgrid[1:]
    d = Data1DInt('tst', xlo, xhi, np.ones_like(xlo))

    mdl = Box1D()
    mdl.xlow = 3.1
    mdl.xhi = 4.2
    mdl.ampl = 0.4

    yexpected = d.eval_model(mdl)
    assert yexpected.min() == pytest.approx(0.0)
    assert yexpected.max() == pytest.approx(0.4 * dx)

    ygot = d.eval_model(mdl.regrid(requested[:-1], requested[1:]))
    assert ygot == pytest.approx(yexpected, abs=tol)


@pytest.mark.parametrize("requested, tol",
                         [(np.arange(2.5, 7, 0.075), 0.007),
                          (np.arange(1, 5.1, 0.075), 0.007),
                          (np.arange(2.5, 5.1, 0.075), 0.007),
                          (np.arange(2.5, 7, 0.12), 0.007),
                          (np.arange(1, 5.1, 0.12), 0.012),
                          (np.arange(2.5, 5.1, 0.12), 0.007),
                          (np.arange(2.5, 7, 0.2), 0.02),
                          (np.arange(1, 5.1, 0.2), 0.02),
                          (np.arange(2.5, 5.1, 0.2), 0.02)])
def test_low_level_regrid1d_int_partial_overlap(requested, tol):
    """What happens if there is partial overlap of the grid?

    See test_low_level_regrid1d_partial_overlap.

    See the comments for when the tol value is used (only for
    edges, not "flat" sections of the model).
    """

    # The range over which we want the model evaluated
    dx = 0.1
    xgrid = np.arange(2, 6, dx)
    xlo = xgrid[:-1]
    xhi = xgrid[1:]
    d = Data1DInt('tst', xlo, xhi, np.ones_like(xlo))

    mdl = Box1D()
    mdl.xlow = 3.1
    mdl.xhi = 4.2
    mdl.ampl = 0.4

    yexpected = d.eval_model(mdl)
    assert yexpected.min() == pytest.approx(0.0)
    assert yexpected.max() == pytest.approx(0.4 * dx)

    rlo = requested[:-1]
    rhi = requested[1:]
    ygot = d.eval_model(mdl.regrid(rlo, rhi))

    # Note that due to the different bin widths of the input and
    # output bins, the start and end bins of the integrated box
    # model will not line up nicely, and so the rising and falling
    # edges may cause differences for more than a single bin. We
    # can think of there being five regions (going from left to
    # right along the grid):
    #
    #    before
    #    rising edge
    #    middle of box
    #    falling edge
    #    after
    #
    # The expected bin values for the middle are 0.4 * dx = 0.04
    # and before and after should be 0. The edges depend on the
    # bin widths, so the tolerance here had been adjusted until
    # the tests pass.
    #
    before = np.arange(0, 10)
    rising = np.arange(10, 12)
    middle = np.arange(12, 21)
    trailing = np.arange(21, 23)
    after = np.arange(23, 39)

    # This ensures that if xgrid is changed (in length at least) we
    # have a reminder to check the above ranges.
    #
    assert len(xlo) == 39

    for idxs in [before, middle, after]:
        assert ygot[idxs] == pytest.approx(yexpected[idxs])

    for idxs in [rising, trailing]:
        assert ygot[idxs] == pytest.approx(yexpected[idxs], abs=tol)


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

    internal_mdl = setup_1d

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

    internal_mdl = setup_1d

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
                                   (0.011, linear_interp)])
def test_regrid1d_interpolation(eincr, rincr, margs, setup_1d):

    tol, method = margs
    _test_regrid1d_interpolation(rtol=tol, method=method,
                                 eval_incr=eincr, req_incr=rincr,
                                 setup_1d=setup_1d)


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

        return ArithmeticFunctionModel(obj)

    def __init__(self, model, wrapper):
        self.model = self.wrapobj(model)
        self.wrapper = wrapper
        CompositeModel.__init__(self, f"{self.wrapper.name}({self.model.name})",
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

    # Ignoring edge effects, the smoothed step function drops from
    # x=100 down to x~120 (it is not clear why it doesn't smooth
    # below x=100, but maybe I've set up the convolution wrong).
    # So, if the regrid model evaluates x=0 to 200 but the requested
    # grid is from x=102 ro 180, then we should see a difference
    # to evaluating without regrid.
    #
    xfull = np.arange(0, 200, 0.5)
    xout = np.arange(101, 180, 0.5)

    # Set up the convolution kernel
    #
    gsmooth = Gauss1D()
    psf = PSFModel('psf', gsmooth)

    # fake up a data object for the fold method
    # TODO: it is not clear to me what grid should be used here;
    #       xfull or xout
    #
    d = Data1D('fake', xfull, xfull * 0)
    psf.fold(d)

    # calculate the expected values
    smoothed = psf(imdl)
    y_full = smoothed(xfull)

    # since the grids have the same binning, it is a simple extraction
    idx0, = np.where(xfull == xout[0])
    idx1, = np.where(xfull == xout[-1])
    assert idx0[0] == 202
    assert idx1[0] == 359
    y_expected = y_full[202:360]

    # This is the model that will be evaluated
    #
    regrid = ModelDomainRegridder1D()
    regrid.grid = xfull

    smoothed_regrid = regrid.apply_to(smoothed)
    y_regrid = smoothed_regrid(xout)

    assert y_regrid == pytest.approx(y_expected, abs=1e-10)

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

    assert y_check == pytest.approx(y_expected, abs=1e-7)


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


@pytest.mark.parametrize('x, y', [
    (None, None),
    (None, [2, 3, 4]),
    ([], [2, 3, 4]),
    ([1, 2], None),
    ([1, 2], [])
])
def test_evaluation_space2d_empty(x, y):
    assert EvaluationSpace2D(x=x, y=y).is_empty


def test_evaluation_space2d_empty_no_args():
    assert EvaluationSpace2D().is_empty


@pytest.mark.parametrize("xaxis", [True, False])
@pytest.mark.parametrize("lo,hi", [([], [2, 3]),
                                   (None, [2, 3]),
                                   ([1], [2, 3]),
                                   ([1, 2], []),
                                   ([1, 2], None),
                                   ([1, 2], [1])
                                   ])
def test_evaluation_space2d_check_axes_are_integrated(xaxis, lo, hi):
    """Check these calls fail.

    Prior to 4.14.1 you could set the integrated axes to
    different sizes.
    """

    olo = [1, 2]
    ohi = [2, 3]
    if xaxis:
        args = {"x": lo, "xhi": hi, "y": olo, "yhi": ohi}
    else:
        args = {"x": olo, "xhi": ohi, "y": lo, "yhi": hi}

    with pytest.raises(DataErr,
                       match=r"^size mismatch between lo and hi: (\d|None) vs (\d|None)$"):
        EvaluationSpace2D(**args)


@pytest.mark.parametrize('xlo, xhi, ylo, yhi, is_integrated', [
    ([1, 2], [2, 3], [1, 2, 3], [2, 3, 4], True),
    ([1, 2], None, [1, 2], None, False),
])
def test_evaluation_space2d_is_integrated(xlo, xhi, ylo, yhi, is_integrated):
    assert EvaluationSpace2D(x=xlo, xhi=xhi, y=ylo, yhi=yhi).is_integrated\
           is is_integrated


@pytest.mark.parametrize('xlo, xhi, ylo, yhi, is_ascending', [
    ([1, 2], [2, 3], [1, 2], [2, 3], (True, True)),
    ([1, 2], [2, 3], [2, 1], [3, 2], (True, False)),
    ([2, 1], [3, 2],  [1, 2], [2, 3], (False, True)),
    ([2, 1], [3, 2], [2, 1], [3, 2], (False, False)),
    ([1, 2], None, [1, 2], None, (True, True)),
    ([1, 2], None, [2, 1], None, (True, False)),
    ([2, 1], None, [1, 2], None, (False, True)),
    ([2, 1], None, [2, 1], None, (False, False)),
])
def test_evaluation_space2d_is_ascending(xlo, xhi, ylo, yhi, is_ascending):
    assert EvaluationSpace2D(x=xlo, xhi=xhi, y=ylo, yhi=yhi).is_ascending \
           == is_ascending


def test_evaluation_space2d_is_ascending_error():
    with pytest.raises(DataErr,
                       match="Axis is empty or has a size of 0"):
        EvaluationSpace2D(x=None, xhi=None, y=None, yhi=None).is_ascending


@pytest.mark.parametrize('xlo, xhi, ylo, yhi', [
    ([1, 2], [2, 3], [1, 2], [2, 3]),
    ([1, 2], [2, 3], [2, 1], [3, 2]),
    ([2, 1], [3, 2], [1, 2], [2, 3]),
    ([2, 1], [3, 2], [2, 1], [3, 2]),
    ([1, 2, 3], None, [1, 2, 3], None),
    ([3, 2, 1.5, 1], None, [1, 1.5, 2, 3], None)
])
def test_evaluation_space2d_start_end(xlo, xhi, ylo, yhi):
    assert EvaluationSpace2D(x=xlo, xhi=xhi, y=ylo, yhi=yhi).start == (1, 1)
    assert EvaluationSpace2D(x=xlo, xhi=xhi, y=ylo, yhi=yhi).end == (3, 3)


@pytest.mark.parametrize('x_overlaps, y_overlaps', [
    (True, True),
    (True, False),
    (False, False),
])
@pytest.mark.parametrize('integrated', [
    True, False
])
def test_evaluation_space2d_overlaps(x_overlaps, y_overlaps, integrated, setup_overlapping_spaces):
    x1, y1, x2, y2, overlaps = setup_overlapping_spaces

    space_one = EvaluationSpace2D(x=x1[0], xhi=x1[1], y=y1[0], yhi=y1[1])
    space_two = EvaluationSpace2D(x=x2[0], xhi=x2[1], y=y2[0], yhi=y2[1])

    assert space_one.overlaps(space_two) is overlaps


def test_evaluation_space2d_grid():
    xlo = [1, 2, 3]
    xhi = [2, 3, 4]
    ylo = [4, 5, 6]
    yhi = [5, 6, 7]

    expected_xlo, expected_ylo = [1, 2, 3, 1, 2, 3, 1, 2, 3], [4, 4, 4, 5, 5, 5, 6, 6, 6]
    expected_xhi, expected_yhi = [2, 3, 4, 2, 3, 4, 2, 3, 4], [5, 5, 5, 6, 6, 6, 7, 7, 7]

    eval_space = EvaluationSpace2D(x=xlo, xhi=xhi, y=ylo, yhi=yhi)

    np.testing.assert_array_equal(eval_space.grid, (expected_xlo, expected_ylo, expected_xhi, expected_yhi))
    # assert eval_space.grid == (expected_xlo, expected_ylo, expected_xhi, expected_yhi)


@pytest.fixture
def setup_overlapping_spaces(integrated, x_overlaps, y_overlaps):
    if integrated:
        x1 = [1.0, 2.0], [2.0, 3.0]
        y1 = [1.0, 2.0], [2.0, 3.0]
        if x_overlaps:
            x2 = [1, 1.5, 2], [2, 2.5, 3]
        else:
            x2 = [1, 2, 3], [1, 2, 3]
        if y_overlaps:
            y2 = [1, 1.5, 2], [2, 2.5, 3]
        else:
            y2 = [-1, 0], [0, 1]
        return x1, y1, x2, y2, (x_overlaps and y_overlaps)

    # To try and mix things up, below I use descending axes as well as a different number of elements
    # in the grid

    x1 = [2.0, 1.0], None
    y1 = [2.0, 1.0], None
    if x_overlaps:
        x2 = [2, 1.5, 1], None
    else:
        x2 = [3, 2, 1], None
    if y_overlaps:
        y2 = [2, 1.5, 1], None
    else:
        y2 = [-1, -2, -3], None
    return x1, y1, x2, y2, (x_overlaps and y_overlaps)


def test_wrong_kwargs():
    xgrid = np.arange(2, 6, 0.1)
    d = Data1D('tst', xgrid, np.ones_like(xgrid))
    mdl = Box1D()
    requested = np.arange(1, 7, 0.1)
    with pytest.raises(TypeError,
                       match="unknown keyword argument: 'fubar'"):
        d.eval_model(mdl.regrid(requested, fubar='wrong_kwargs'))


def test_interp_method_is_callable():
    """Check that method is callable."""
    rmdl = ModelDomainRegridder1D()
    with pytest.raises(TypeError,
                       match="method argument 'True' is not callable"):
        rmdl.method = True


def test_axis_check_not_a_scalar():
    """Mainly to ensure this error condition is checked"""

    with pytest.raises(DataErr,
                       match="Array must be a sequence or None"):
        PointAxis(23)


def test_axis_check_not_multidim():
    """Mainly to ensure this error condition is checked"""

    with pytest.raises(DataErr,
                       match="Array must be 1D"):
        PointAxis(np.arange(12).reshape(3, 4))


def test_pointaxis_not_integrated():
    """Check for the empty case"""

    assert not PointAxis([]).is_integrated


def test_integratedaxis_is_integrated():
    """Check for the empty case"""

    assert IntegratedAxis([], []).is_integrated


def test_pointaxis_size():
    """Pick a descending axis for fun"""
    assert PointAxis([4, 3, 2]).size == 3


def test_integratedaxis_size():
    assert IntegratedAxis([4, 3, 2], [4.5, 3.5, 3]).size == 3


def test_pointaxis_start():
    """Pick a descending axis for fun"""
    assert PointAxis([4, 3, 2]).start == 2


def test_pointaxis_end():
    """Pick a descending axis for fun"""
    assert PointAxis([4, 3, 2]).end == 4


def test_integratedaxis_start():
    """Pick a descending axis for fun"""
    assert IntegratedAxis([4, 3, 2], [4.5, 3.5, 3]).start == 2


def test_integratedaxis_end():
    """Pick a descending axis for fun"""
    assert IntegratedAxis([4, 3, 2], [4.5, 3.5, 3]).end == pytest.approx(4.5)


def test_evaluationspace1d_zeros_like_empty():

    assert EvaluationSpace1D().zeros_like() == pytest.approx([])


def test_evaluationspace1d_zeros_like_point():

    assert EvaluationSpace1D([3, 2, -1]).zeros_like() == pytest.approx([0, 0, 0])


def test_evaluationspace1d_zeros_like_integrated():

    assert EvaluationSpace1D([3, 2, -1], [2, 1, -3]).zeros_like() == pytest.approx([0, 0, 0])


@pytest.mark.parametrize("flag,mdl",
                         [(False, Box1D() * Const1D() * Gauss1D()),
                          (False, Box1D() * Const1D() * (1 + Gauss1D())),
                          (False, 1 / (1 + Const1D())),
                          (False, 1 / (Box1D() * Gauss1D())),
                          (False, 1 / (1 + Box1D() * Gauss1D())),
                          (False, ((Box1D() + Const1D()) / (Box1D() + Gauss1D()))),
                          (False, (Const1D() / (Box1D() + Gauss1D()) + Box1D())),
                          (False, (1 / (1 - (-Const1D())))),
                          (False, ((-Box1D()) / ((-Const1D()) - (-Gauss1D())))),
                          (True, Box1D() * Gauss1D()),
                          (True, Box1D() / (Const1D() + Gauss1D())),
                          (True, Const1D() / (1 + Box1D() * Gauss1D())),
                          (True, (Box1D() + Const1D() / (Box1D() + Gauss1D())))
                          ])
@pytest.mark.parametrize("args", [([1, 2, 3], ), ([1, 2, 4], [2, 3, 5])])
def test_recursion_1802(flag, mdl, args):
    """Can we create moderately-complex models to regrid?

    Test cases from #1802. It is unlikely that the axis choice (point
    or integrated) makes a difference, but add a test just in case.

    """

    # This test could be marked XFAIL but I want to make it obvious
    # once it gets fixed.
    #
    if flag:
        rmdl = mdl.regrid(*args)
        # basic check to see if we have got a regrid model
        assert isinstance(rmdl, RegridWrappedModel)
        assert rmdl.name == f"regrid1d({mdl.name})"
        return

    with pytest.raises(RecursionError):
        mdl.regrid(*args)


def test_unop_regrid():
    """Can we regrid a unary op model?"""

    mdl = Gauss1D()
    mdl.pos = 100
    mdl.ampl = 10

    egrid = [80, 90, 95, 107, 115]
    expected = -mdl(egrid)

    nmdl = -mdl
    with pytest.raises(AttributeError):
        rgrid = nmdl.regrid(np.arange(70, 120, 2))

    # assert rgrid(egrid) == pytest.approx(expected)


def test_unop_binop_combo():
    """What happens if a binop is given two unops: can we regrid?"""

    mdl1 = Box1D()
    mdl2 = Gauss1D()
    mdl = (-mdl1 - (-mdl2))
    with pytest.raises(ModelErr,
                       match="Neither component supports regrid method"):
        mdl.regrid([1, 2, 3])


def test_regrid_binop_arithmeticconstantmodel():
    """If you have managed to create binop(constants) then regridded, does it work"""

    grid = np.linspace(20, 30, 21)
    orig = BinaryOpModel(4, 6, np.add, '+')

    with pytest.raises(ModelErr,
                       match="Neither component supports regrid method"):
        orig.regrid(grid)


def test_box1d_point():
    """This model is used in several tests so check it out"""

    cpt = Box1D()
    cpt.xlow = 130
    cpt.xhi = 189

    mdl = cpt
    exp = np.asarray([0, 0, 1, 1, 1, 1, 1])

    xgrid = np.arange(105, 195, 13)
    assert mdl(xgrid) == pytest.approx(exp)

    xbase = np.arange(100, 200, 10)
    rmdl = mdl.regrid(xbase)
    assert rmdl(xgrid) == pytest.approx(exp)


def test_constant_box1d_point():
    """Check 3 + mdl"""

    cpt = Box1D()
    cpt.xlow = 130
    cpt.xhi = 189

    mdl = 3 + cpt

    exp = 3 + np.asarray([0, 0, 1, 1, 1, 1, 1])

    xgrid = np.arange(105, 195, 13)
    assert mdl(xgrid) == pytest.approx(exp)

    xbase = np.arange(100, 200, 10)
    rmdl = mdl.regrid(xbase)
    assert rmdl(xgrid) == pytest.approx(exp)


def test_box1d_bin():
    """This model is used in several tests so check it out"""

    cpt = Box1D()
    cpt.xlow = 130
    cpt.xhi = 189

    mdl = cpt
    exp = np.asarray([0, 1, 13, 13, 13, 13])

    # Ideally this would be exp, but it's not quite
    exp2 = np.asarray([0, 1, 13, 13, 13, 12.7])

    xgrid = np.arange(105, 195, 13)
    xg1 = xgrid[:-1]
    xg2 = xgrid[1:]
    assert mdl(xg1, xg2) == pytest.approx(exp)

    xbase = np.arange(100, 200, 10)
    rmdl = mdl.regrid(xbase[:-1], xbase[1:])
    assert rmdl(xg1, xg2) == pytest.approx(exp2)


def test_constant_box1d_bin():
    """Check 3 + mdl"""

    cpt = Box1D()
    cpt.xlow = 130
    cpt.xhi = 189

    mdl = 3 + cpt

    # What is the correct value here? That is, what do
    # we expect the "integrated" version of the constant
    # to be? It should just be 3 (i.e. it doesn't care about
    # the bin width), but once we start rebinning things it
    # gets complicated.
    #
    exp = 3 + np.asarray([0, 1, 13, 13, 13, 13])

    # Just report the current value
    exp2 = 3.9 + np.asarray([0, 1, 13, 13, 13, 12.7])

    xgrid = np.arange(105, 195, 13)
    xg1 = xgrid[:-1]
    xg2 = xgrid[1:]
    assert mdl(xg1, xg2) == pytest.approx(exp)

    xbase = np.arange(100, 200, 10)
    rmdl = mdl.regrid(xbase[:-1], xbase[1:])
    assert rmdl(xg1, xg2) == pytest.approx(exp2)


@pytest.mark.parametrize("integrate", [True, False])
def test_deep_binop_points(integrate):
    """Can we handle a "deep" binop tree?

    (2 * model) / (3 + model)

    We evaluate it on a set of points, not a grid,
    so the integration setting doesn't matter.
    """

    yexp = np.asarray([0, 2/3, 0.5, 0.5, 0.5, 0, 0])

    xbase = np.arange(100, 200, 10)
    xgrid = np.arange(105, 195, 13)

    m1 = Box1D('m1')
    m2 = Box1D('m2')
    m1.xlow = 110
    m1.xhi = 160
    m2.xlow = 130
    m2.xhi = 189

    m1.integrate = integrate
    m2.integrate = integrate

    mdl = (2 * m1) / (3 + m2)
    with pytest.raises(RecursionError):
        regrid = mdl.regrid(xbase)

    # ym = mdl(xgrid)
    # assert ym == pytest.approx(yexp)

    # yr = regrid(xgrid)
    # assert yr == pytest.approx(ym)


@pytest.mark.parametrize("integrate,yexp,yexp2",
                         [(True,
                           [16/3, 26/4, 26/16, 26/16, 6/16, 0, 0],
                           [40/7.5, 8.15384615, 2, 2, 0.46153846, 0, 0]),
                          (False,
                           [0, 2/3, 0.5, 0.5, 0.5, 0, 0],
                           [4/7.5, 0.85, 0.65, 0.65, 0.65, 0, 0])
                          ])
def test_deep_binop_bins(integrate, yexp, yexp2):
    """Can we handle a "deep" binop tree?

    (2 * model) / (3 + model)

    We evaluate it on bins (lo,hi) so the integration
    setting does matter here.

    We have to give separate "expected" results for the un-regridded
    and regridded models as they are not close enough to use the same
    values. The values for integrate=False (yexp, yexp2) may be wrong,
    but hard to diagnose when issue #1802 is not fixed.

    """

    yexp = np.asarray(yexp)
    yexp2 = np.asarray(yexp2)

    xbase = np.arange(100, 200 + 10, 10)
    xgrid = np.arange(105, 195 + 13, 13)

    m1 = Box1D('m1')
    m2 = Box1D('m2')
    m1.xlow = 110
    m1.xhi = 160
    m2.xlow = 130
    m2.xhi = 189

    m1.integrate = integrate
    m2.integrate = integrate

    mdl = (2 * m1) / (3 + m2)
    with pytest.raises(RecursionError):
        regrid = mdl.regrid(xbase[:-1], xbase[1:])

    # ym = mdl(xgrid[:-1], xgrid[1:])
    # assert ym == pytest.approx(yexp)

    # yr = regrid(xgrid[:-1], xgrid[1:])
    # assert yr == pytest.approx(yexp2)


class SimpleModel(RegriddableModel1D):
    """A model which uses the middle of integrated bins"""

    def __init__(self, name='simplemodel'):
        self.scale = Parameter(name, "scale", 0, min=-5, max=5)
        RegriddableModel1D.__init__(self, name, (self.scale, ))

    def calc(self, p, xlo, xhi=None, **kwargs):
        "This ignores the integrate setting"
        x = np.asarray(xlo)
        if xhi is not None:
            x = x + np.asarray(xhi)
            x /= 2

        return p[0] * x


def get_non_integrated_model_on_integrated_axis_1d(integrated):

    mdl = SimpleModel()
    mdl.scale = 2
    assert mdl.integrate is True  # regression test

    mdl.integrate = integrated

    eval_grid = [0.5, 1, 1.5, 2, 2.5, 3]
    rmdl = mdl.regrid(eval_grid[:-1], eval_grid[1:])

    # These bins do not fill the space (ie there are gaps between
    # bins).
    #
    req_lo = [0.6, 0.8, 1.0, 1.4, 2.0, 2.6]
    req_hi = [0.8, 1.0, 1.2, 2.0, 2.4, 2.8]

    # We expect 2 * (xlo + xhi) / 2, = xlo + xh, so
    #
    #   [ 1.4, 1.8, 2.2, 3.4, 4.4, 5.4]
    #
    # at least when "non integrated".
    #
    return rmdl(req_lo, req_hi)


def test_non_integrated_model_on_integrated_axis_1d_false():
    """Check what happens with a 'multiplicative' model given an integrated axis."""

    got = get_non_integrated_model_on_integrated_axis_1d(False)

    # The relative tolerance is guessed here as the test currently
    # fails quite badly for the first bin.
    #
    # expected = [1.4, 1.8, 2.2, 3.4, 4.4, 5.4]
    # assert got == pytest.approx(expected, rel=0.1)

    # For the moment this is a regression test, so we just check
    # the values we get from running the code on 4.15.1 era code.
    #
    expected = [0, 1.8, 2.2, 3.4, 4.4, 5.4]
    assert got == pytest.approx(expected)


def test_non_integrated_model_on_integrated_axis_1d_true():
    """Check what happens with a 'multiplicative' model given an integrated axis."""

    got = get_non_integrated_model_on_integrated_axis_1d(True)

    # The relative tolerance is guessed here as the test currently
    # fails quite badly (as there's uncertainty over the behavior).
    #
    # expected = [1.4, 1.8, 2.2, 3.4, 4.4, 5.4]
    # assert got == pytest.approx(expected, rel=0.1)

    # For the moment this is a regression test, so we just check
    # the values we get from running the code on 4.15.1 era code.
    #
    expected = [0.6, 0.6, 1, 4, 3.6, 2.2]
    assert got == pytest.approx(expected)


@pytest.mark.parametrize("cls", [EvaluationSpace1D, EvaluationSpace2D])
def test_evaluationspace_empty_is_empty(cls):
    """Simple check"""
    espace = cls()
    assert espace.is_empty is True


@pytest.mark.parametrize("cls", [EvaluationSpace1D, EvaluationSpace2D])
def test_evaluationspace_empty_is_integrated(cls):
    """Simple check"""
    espace = cls()
    assert espace.is_integrated is False


@pytest.mark.parametrize("cls", [EvaluationSpace1D, EvaluationSpace2D])
def test_evaluationspace_empty_is_ascending(cls):
    """Simple check"""
    espace = cls()
    with pytest.raises(DataErr, match="^Axis is empty or has a size of 0$"):
        # This is a property
        _ = espace.is_ascending


@pytest.mark.parametrize("cls", [EvaluationSpace1D, EvaluationSpace2D])
@pytest.mark.parametrize("meth", ["start", "end"])
def test_evaluationspace_empty_range(cls, meth):
    """Simple check"""
    espace = cls()
    with pytest.raises(DataErr, match="^Axis is empty or has a size of 0$"):
        # These are properties, so accessing the field causes the error
        _ = getattr(espace, meth)


@pytest.mark.parametrize("cls,expected",
                         [(EvaluationSpace1D, (None, )),
                          (EvaluationSpace2D, ([None], [None]))
                          ])
def test_evaluationspace_empty_grid(cls, expected):
    """Simple check.

    EvaluationSpace1D has midpoint_grid but the 2D version does
    not.
    """
    espace = cls()
    assert espace.grid == expected


@pytest.mark.parametrize("cls", [EvaluationSpace1D, EvaluationSpace2D])
def test_evaluationspace_empty_zeroes(cls):
    """Simple check"""
    espace = cls()
    assert espace.zeros_like().shape == (0, )


@pytest.mark.parametrize("cls", [EvaluationSpace1D, EvaluationSpace2D])
def test_evaluationspace_empty_overlaps(cls):
    """Simple check"""
    espace1 = cls()
    espace2 = cls()
    with pytest.raises(DataErr, match="^Axis is empty or has a size of 0$"):
        espace1.overlaps(espace2)


@pytest.mark.parametrize("cls", [EvaluationSpace1D])
def test_evaluationspace_empty_contains(cls):
    """This is a regression test.

    EvaluationSpace2D does not have a __contains__ method.
    """
    espace1 = cls()
    espace2 = cls()
    with pytest.raises(DataErr, match="^Axis is empty or has a size of 0$"):
        _ = espace1 in espace2
