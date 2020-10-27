# Copyright 2018, 2020 Smithsonian Astrophysical Observatory
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
# following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
# disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
# products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
import numpy as np
import pytest
from pytest import approx

from sherpa.astro.data import DataIMG, DataIMGInt
from sherpa.astro.ui.utils import Session
from sherpa.data import Data1DInt, Data1D
from sherpa.models.basic import Box1D
from sherpa.models import Const1D, RegriddableModel1D, Parameter, Const2D, \
    RegriddableModel2D, ArithmeticModel, Gauss2D, basic, model
from sherpa.utils.err import ModelErr
from sherpa.utils import neville, linear_interp
from sherpa.utils import akima

@pytest.fixture
def setup():
    const = Const1D("const")
    const.c0 = 0
    const.c0.freeze()

    my_model = MyModel("my_model")
    my_model.integrate = False

    return Session(), my_model, const


@pytest.fixture
def setup2d():
    const = Const2D("const")
    const.c0 = 0
    const.c0.freeze()

    x = [2, 3, 2, 3]
    y = [2, 2, 3, 3]
    xhi = [2.1, 3.5, 2.1, 3.5]
    yhi = [2.1, 2.1, 3, 3.5]

    # This is the result when rebinning [100, ] * 4
    z = [225, ] * 4

    my_model = MyModel2D("my_model")

    return Session(), my_model, const, (x, y, xhi, yhi, z)


def test_evaluate_model_on_arbitrary_grid_point_list(setup):
    """
    The idea of the test is that the model will evaluate differently depending on the grid it is evaluated on.
    This is really arbitrary, it just exercises the high level API for a common workflow while making sure the results
    are the expected ones.
    """
    ui, my_model, const = setup

    # Load data
    ui.load_arrays(1, [1, 2, 3], [100, 100, 100])

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid([1, 2, 2.5, 4, 5])

    # The model will usually be part of a complex model expression, so let's pretend we add another component,
    # although that component is muted.
    ui.set_source(regrid_model + const)

    # Fit and check the result
    assert_fit(ui, my_model, 1)

    # Now fit with a different grid.
    # This is also the important part.
    regrid_model.grid = [1, 2, 3, 4, 5]

    assert_fit(ui, my_model, 0)


def test_evaluate_model_on_arbitrary_grid_point_list_2d(setup2d):
    """
    The idea of the test is that the model will evaluate differently depending on the grid it is evaluated on.
    This is really arbitrary, it just exercises the high level API for a common workflow while making sure the results
    are the expected ones.
    """
    ui, my_model, const, data = setup2d
    x, y, _, _, z = data

    # Load data
    ui.load_arrays(1, x, y, z, DataIMG)

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid([2, 2.5, 3], [2, 2.5, 3])

    # The model will usually be part of a complex model expression, so let's pretend we add another component,
    # although that component is muted.
    ui.set_source(regrid_model + const)

    # Fit and check the result
    assert_fit(ui, my_model, (1, 1))

    # Now fit with a different grid.
    # This is also the important part.
    regrid_model.grid = [2, 3], [2, 3]

    assert_fit(ui, my_model, (0, 0))


def test_evaluate_model_on_arbitrary_grid_integrated_list(setup):
    """
    Same as above, but with integrated models.
    """
    ui, my_model, const = setup

    # Load data
    ui.load_arrays(1, [1.5, 2.5, 3.5], [2.5, 3.5, 4.5], [100, 100, 100], Data1DInt)

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid([0, 1, 2], [1, 2, 3])

    # The model will be part of a complex model expression, so let's pretend we add another component
    ui.set_source(regrid_model + const)

    # Fit and check the result
    assert_fit(ui, my_model, 1)

    # Now fit with a different grid.
    # This is also the important part.
    regrid_model.grid = [1.5, 2.5, 3.5], [2.5, 3.5, 4.5]

    assert_fit(ui, my_model, 0)


def test_evaluate_model_on_arbitrary_grid_integrated_list_2d(setup2d):
    """
    Same as above, but with integrated models
    """
    ui, my_model, const, data = setup2d
    x, y, xhi, yhi, z = data

    # Load data
    ui.load_arrays(1, x, y, xhi, yhi, z, DataIMGInt)

    regrid_lo = [2, 2.5, 3]
    regrid_hi = np.array([2, 2.5, 3.5])

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid(regrid_lo, regrid_lo, regrid_hi, regrid_hi)

    # The model will usually be part of a complex model expression, so let's pretend we add another component,
    # although that component is muted.
    ui.set_source(regrid_model + const)

    # Fit and check the result
    assert_fit(ui, my_model, (1, 1))

    # Now fit with a different grid.
    # This is also the important part.
    regrid_model.grid = x, y, xhi, yhi

    assert_fit(ui, my_model, (0, 0))


def test_evaluate_model_on_arbitrary_grid_point_ndarray(setup):
    """
    The idea of the test is that the model will evaluate differently depending on the grid it is evaluated on.
    This is really arbitrary, it just exercises the high level API for a common workflow while making sure the results
    are the expected ones.
    """
    ui, my_model, const = setup

    # Load data
    ui.load_arrays(1, [1, 2, 3], [100, 100, 100])

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid(np.array([1, 2, 2.5, 4, 5]))

    # The model will be part of a complex model expression, so let's pretend we add another component
    ui.set_source(regrid_model + const)

    # Fit and check the result
    assert_fit(ui, my_model, 1)

    # Now fit with a different regrid.
    # This is also the important part.
    regrid_model.grid = np.array([1, 2, 3, 4, 5])

    assert_fit(ui, my_model, 0)


def test_evaluate_model_on_arbitrary_grid_integrated_ndarray(setup):
    """
    Same as above, but with integrated models.
    """
    ui, my_model, const = setup

    # Load data
    ui.load_arrays(1, [1.5, 2.5, 3.5], [2.5, 3.5, 4.5], [100, 100, 100], Data1DInt)

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid(np.array([0, 1, 2]), [1, 2, 3])

    # The model will be part of a complex model expression, so let's pretend we add another component
    ui.set_source(regrid_model + const)

    # Fit and check the result
    assert_fit(ui, my_model, 1)

    # Now fit with a different grid.
    # This is also the important part.
    regrid_model.grid = [1.5, 2.5, 3.5], np.array([2.5, 3.5, 4.5])

    assert_fit(ui, my_model, 0)


def test_evaluate_model_on_arbitrary_grid_no_overlap(setup):
    """
    If grids do not overlap, issue a warning and return zeros
    """
    ui, my_model, _ = setup

    # Get a model that evaluates on a different grid
    # This is the important part. Note that there is overlap, but
    # the start and end p
    with pytest.raises(ModelErr) as excinfo:
        regrid_model = my_model.regrid([2, 2.5], [2, 2.5])
    assert ModelErr.dict['needsint'] in str(excinfo.value)


def test_evaluate_model_on_arbitrary_grid_no_overlap_2d(setup2d):
    """
    In the 2D case, the overlap is way more stringent than in the 1D case, due to the complexity of rebinning
    """
    ui, my_model, _, data = setup2d
    x, y, _, _, _ = data

    my_model.x_has_25 = 1  # To force the model to evaluate to something other than 0.

    # Get a model that evaluates on a different grid
    # This is the important part. Note that there is overlap, but
    # the start and end points are different.
    regrid_model = my_model.regrid([2, 2.5], [2, 2.5])

    with pytest.warns(UserWarning):
        np.testing.assert_array_equal(regrid_model(x, y), [0, 0, 0, 0])


def test_runtime_interp():
    def tst_runtime_interp(model, requested, interp):
        regrid_model = mdl.regrid(requested, interp=interp)
        yregrid = regrid_model(xgrid)
        return yregrid

    xgrid = np.arange(2, 6, 0.1)
    requested = np.arange(2.5, 5.1, 0.075)
    mdl = Box1D()
    mdl.xlow = 3.1
    mdl.xhi = 4.2
    mdl.ampl = 0.4
    yregrid = tst_runtime_interp(mdl, requested, akima.akima)
    assert 4.4 == approx(yregrid.sum())
    yregrid = tst_runtime_interp(mdl, requested, linear_interp)
    assert 4.4 == approx(yregrid.sum())
    yregrid = tst_runtime_interp(mdl, requested, neville)
    assert - 5.0e6 > yregrid.sum()

    d = Data1D('tst', xgrid, np.ones_like(xgrid))
    yexpected = d.eval_model(mdl)
    requested = np.arange(2.5, 7, 0.2)
    rmdl = mdl.regrid(requested)
    ygot = d.eval_model(rmdl)
    assert ygot == approx(yexpected)


def test_regrid_binaryop_1d():
    """issue #762, Cannot regrid a composite model (BinaryOpModel)"""
    from sherpa.stats import LeastSq
    from sherpa.fit import Fit
    from sherpa.optmethods import LevMar


    class MyConst1D(RegriddableModel1D):

        def __init__(self, name='myconst1d'):
            self.c0 = Parameter(name, 'c0', 3.1)
            self.counter = 0
            ArithmeticModel.__init__(self, name, (self.c0,))

        def calc(self, par, *args, **kwargs):
            x = args[0]
            self.counter += x.size
            return par[0]


    class MyGauss(RegriddableModel1D):

        def __init__(self, name='mygauss'):
            self.sigma = Parameter(name, 'sigma', 10, min=0, max=10)
            self.pos = Parameter(name, 'pos', 0, min=-10, max=10)
            self.ampl = Parameter(name, 'ampl', 5)
            self.counter = 0
            ArithmeticModel.__init__(self, name, (self.sigma, self.pos, self.ampl))

        def calc(self, par, *args, **kwargs):
            sigma, pos, ampl = par[0], par[1], par[2]
            x = args[0]
            self.counter += x.size
            return ampl * np.exp(-0.5 * (args[0] - pos)**2 / sigma**2)


    np.random.seed(0)
    leastsq = LeastSq()
    levmar = LevMar()
    mygauss = MyGauss()
    myconst = MyConst1D()
    mymodel = mygauss + myconst
    x = np.linspace(-5., 5., 5)
    err = 0.25
    y = mymodel(x) + np.random.normal(mygauss.pos.val, err, x.shape)
    mygauss.counter = 0
    myconst.counter = 0
    data = Data1D('one', x, y)
    fit = Fit(data, mymodel, leastsq, levmar)
    result = fit.fit()
    assert result.numpoints == x.size
    assert result.statval < 1.0
    assert mygauss.counter == myconst.counter
    assert (result.nfev + 4) * x.size == mygauss.counter

    mygauss.counter = 0
    myconst.counter = 0
    x_regrid = np.linspace(-5., 5., 25)
    mymodel_regrid = mymodel.regrid(x_regrid)
    fit = Fit(data, mymodel_regrid, leastsq, levmar)
    result = fit.fit()
    assert result.numpoints == x.size
    assert result.statval < 1.0
    assert mygauss.counter == myconst.counter
    assert (result.nfev + 4) * x_regrid.size == mygauss.counter

def test_regrid_binaryop_2d():
    y0, x0 = np.mgrid[20:29, 10:20]
    y0 = y0.flatten()
    x0 = x0.flatten()

    gmdl = Gauss2D()
    gmdl.fwhm = 14
    gmdl.xpos = 15
    gmdl.ypos = 24
    gmdl.ampl = 10

    cmdl = Const2D()
    cmdl.c0 = 4

    xr1 = np.arange(10, 20, 1)
    yr1 = np.arange(20, 29, 1)

    rmdlg = gmdl.regrid(xr1, yr1)
    rmdlc = cmdl.regrid(xr1, yr1)

    shape = y0.shape
    truthg = gmdl(x0, y0).reshape(shape)
    truthc = cmdl(x0, y0).reshape(shape)
    truth = truthg + truthc

    ans1 = rmdlg(x0, y0).reshape(shape)
    ans2 = rmdlc(x0, y0).reshape(shape)
    assert (ans1 == truthg).all()
    assert (ans2 == truthc).all()

    rmdl = (gmdl + cmdl).regrid(xr1, yr1)
    ans3 = rmdl(x0, y0).reshape(shape)
    assert (ans3 == truth).all()


def test_regrid_call_behavior():
    class Wrappable1D(model.RegriddableModel1D):

        def __init__(self, cls, name):
            self.ncalled = []  # record the number of elements
            self.baseclass = cls
            self.baseclass.__init__(self, name)

        def calc(self, pars, xlo, *args, **kwargs):
            xlo = np.asarray(xlo)
            self.ncalled.append((xlo[0], xlo[-1], xlo.size))
            return self.baseclass.calc(self, pars, xlo, *args, **kwargs)


    m1 = Wrappable1D(basic.Const1D, 'm1')
    m2 = Wrappable1D(basic.Gauss1D, 'm2')

    m2.pos = 5

    xregrid = np.arange(0, 20, 0.2)
    xdata = np.arange(1.5, 12.5, 0.5)

    morig = m1 + m2
    mwrap = morig.regrid(xregrid)

    y = mwrap(xdata)

    # Check both components were called with the same grid
    assert m1.ncalled == m2.ncalled

    # Check that m1 was called with the expected grid (ie that
    # it is larger than xdata).
    got = m1.ncalled
    assert len(got) == 1
    minval, maxval, nbins = m1.ncalled[0]
    assert minval == pytest.approx(0)
    assert maxval == pytest.approx(19.8)

    assert nbins > xdata.size
    assert nbins == 111


class MyModel(RegriddableModel1D):
    """
    A model that returns [100, ] * len(x) if 2.5 is in the input array x
    """
    def __init__(self, name):
        self.has_25 = Parameter(name, "has_25", 0, min=0, max=1)
        ArithmeticModel.__init__(self, name, (self.has_25,))

    def guess(self, dep, *args, **kwargs):
        raise NotImplementedError()

    def get_center(self):
        raise NotImplementedError()

    def set_center(self, *args, **kwargs):
        raise NotImplementedError()

    def calc(self, p, *args, **kwargs):
        x = args[0]

        if 2.5 not in x:
            if p[0] == 0:
                return [100, ] * len(x)
            else:
                return [100-p[0]*100,] * len(x)
        if p[0] == 1:
            return [100, ] * len(x)
        else:
            return [p[0]*100, ] * len(x)


class MyModel2D(RegriddableModel2D):
    """
    A 2D model that returns [100, ] * len(x) * len(y) if 2.5 is in the input arrays x and y
    """
    def __init__(self, name):
        self.x_has_25 = Parameter(name, "x_has_25", 0, min=0, max=1)
        self.y_has_25 = Parameter(name, "y_has_25", 0, min=0, max=1)
        RegriddableModel2D.__init__(self, name, (self.x_has_25, self.y_has_25))

    def guess(self, dep, *args, **kwargs):
        raise NotImplementedError()

    def get_center(self):
        raise NotImplementedError()

    def set_center(self, *args, **kwargs):
        raise NotImplementedError()

    def calc(self, p, *args, **kwargs):
        x, y, x_has_25, y_has_25 = args[0], args[1], p[0], p[1]

        x_eval = np.array(self._eval(x, x_has_25))
        y_eval = np.array(self._eval(y, y_has_25))

        return (x_eval + y_eval)/2

    def _eval(self, array, has_25):
        if 2.5 not in array:
            if has_25 == 0:
                return [100, ] * len(array)
            else:
                return [100-has_25*100,] * len(array)
        if has_25 == 1:
            return [100, ] * len(array)
        else:
            return [has_25*100, ] * len(array)


def assert_fit(ui, model, value):
    ui.fit()

    try:  # 2D, two values
        len(value)
        assert model.x_has_25.val == approx(value[0])
        assert model.y_has_25.val == approx(value[1])
    except TypeError:  # 1D, one value
        assert model.has_25.val == approx(value)
