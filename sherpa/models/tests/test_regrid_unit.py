#
#  Copyright (C) 2017, 2018, 2019, 2020, 2021
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

import numpy as np
from numpy.testing import assert_allclose

import pytest

from sherpa.models.model import Model, ArithmeticModel, CompositeModel, \
    ArithmeticFunctionModel, RegridWrappedModel
from sherpa.models.basic import Box1D, Const1D, Gauss1D, \
    PowLaw1D, StepLo1D
from sherpa.models.parameter import Parameter
from sherpa.instrument import PSFModel
from sherpa.data import Data1D, Data1DInt
import sherpa.astro.ui as ui

import sherpa.utils
from sherpa.utils.err import ModelErr

from sherpa.models.regrid import ModelDomainRegridder1D, EvaluationSpace1D, \
    EvaluationSpace2D


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
    assert store[1] == {'integrate': True}

    store = store[0]
    assert store[0] == [-34.5]
    combine = np.unique(np.append(grid_expected, grid_requested))
    indices = combine.searchsorted(store[1])
    assert (store[1] == combine[indices]).all()


def test_regrid1d_error_calc_no_args(setup_1d):

    internal_mdl = setup_1d
    grid_evaluate = EvaluationSpace1D(np.arange(-10, 100, 5))

    rmdl = ModelDomainRegridder1D(grid_evaluate)
    mdl = rmdl.apply_to(internal_mdl)

    with pytest.raises(ModelErr) as excinfo:
        pvals = [p.val for p in internal_mdl.pars]
        mdl.calc(p=pvals)

    assert ModelErr.dict['nogrid'] in str(excinfo.value)


def test_regrid1d_error_grid_mismatch_1(setup_1d):
    """Internal grid is integrated but given points"""

    internal_mdl = setup_1d

    grid_evaluate = np.arange(-10, 100, 5)
    eval_space = EvaluationSpace1D(grid_evaluate[:-1], grid_evaluate[1:])
    rmdl = ModelDomainRegridder1D(eval_space)

    mdl = rmdl.apply_to(internal_mdl)
    grid_run = np.arange(0, 20, 10)
    with pytest.raises(ModelErr) as excinfo:
        mdl(grid_run)

    assert ModelErr.dict['needsint'] in str(excinfo.value)


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
    with pytest.raises(ModelErr) as excinfo:
        b1.regrid(grid_lo, grid_hi)

    assert ModelErr.dict['needsint'] in str(excinfo.value)


def test_low_level_regrid1d_non_overlapping_not_allowed():
    """Integrated data space must not overlap"""

    c = Box1D()
    lo = np.linspace(1, 100, 600)
    hi = np.linspace(2, 101, 600)
    with pytest.raises(ModelErr) as excinfo:
        c.regrid(lo, hi)

    assert ModelErr.dict['needsint'] in str(excinfo.value)


def test_regrid1d_error_grid_mismatch_2(setup_1d):
    """Internal grid is points but given integrated"""

    internal_mdl = setup_1d

    grid_evaluate = EvaluationSpace1D(np.arange(-10, 100, 5))
    rmdl = ModelDomainRegridder1D(grid_evaluate)

    mdl = rmdl.apply_to(internal_mdl)
    grid_run = np.arange(0, 20, 10)
    with pytest.raises(ModelErr) as excinfo:
        mdl(grid_run[:-1], grid_run[1:])

    assert ModelErr.dict['needspoint'] in str(excinfo.value)


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
                                   (0.011, sherpa.utils.linear_interp)])
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

@pytest.mark.xfail  # see #1334
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
    (None, [1, 2]),
    ([1, 2], None),
    (None, None),
    ([], [1, 2]),
    ([1, 2], [])
])
def test_evaluation_space2d_empty(x, y):
    assert EvaluationSpace2D(x=x, y=y).is_empty


def test_evaluation_space2d_empty_no_args():
    assert EvaluationSpace2D().is_empty


@pytest.mark.parametrize('xlo, xhi, ylo, yhi, is_integrated', [
    ([1, 2], [2, 3], [1, 2], [2, 3], True),
    ([1, 2], [2, 3], None, None, False),
    ([1, 2], [2, 3],  [], [2, 3], False),
    ([1, 2], [2, 3], [1, 2], [], False),
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
    with pytest.raises(ValueError):
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
    with pytest.raises(TypeError) as excinfo:
        d.eval_model(mdl.regrid(requested, fubar='wrong_kwargs'))

    assert "unknown keyword argument: 'fubar'" in str(excinfo.value)


def test_interp_method_is_callable():
    """Check that method is callable."""
    rmdl = ModelDomainRegridder1D()
    with pytest.raises(TypeError) as exc:
        rmdl.method = True

    assert str(exc.value) == "method argument 'True' is not callable"
