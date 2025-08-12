#
#  Copyright (C) 2025
#  MIT
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
"""Test interface to Scipy optimization methods.

These tests are intentionally sparse, as the optimizers themselves are
already tested in the scipy library. Here, we only want to test that the
interface works and parameters are passed to scipy correctly.
"""

import pytest
pytest.importorskip("scipy", reason="scipy not installed")

import numpy as np

from sherpa import ui
from sherpa import optmethods
from sherpa.fit import Fit
from sherpa.data import Data1D
from sherpa.models import Const1D, Gauss1D, Delta1D
from sherpa.stats import CStat, Chi2
from sherpa.optmethods import optscipy
from sherpa.models.parameter import hugeval


data = Data1D('data1', x=np.arange(10),
              y=[34.5, 23.4, 22.3, 45.6, 56.7, 67.8, 58.9, 43.0, 30.1, 25.2],
              staterror=5 * np.ones(10))



def test_convert_bounds_to_scipy():
    """Test the conversion of bounds to scipy format."""

    # Test case 1: No bounds
    xmin = np.array([-hugeval, -hugeval])
    xmax = np.array([hugeval, hugeval])
    scipy_bounds = optscipy.convert_bounds_to_scipy(xmin, xmax, requires_finite_bounds=False)
    assert scipy_bounds is None

    # Test case 2: Finite bounds
    xmin = np.array([-10, 0])
    xmax = np.array([10, 5])
    scipy_bounds = optscipy.convert_bounds_to_scipy(xmin, xmax, requires_finite_bounds=False)
    assert scipy_bounds == [(-10, 10), (0, 5)]

    # Test case 3: Mixed bounds
    xmin = np.array([-10, -hugeval, 0])
    xmax = np.array([hugeval, 5, 10])
    scipy_bounds = optscipy.convert_bounds_to_scipy(xmin, xmax, requires_finite_bounds=False)
    assert scipy_bounds == [(-10, None), (None, 5), (0, 10)]


def test_convert_bounds_to_scipy_raises():
    """Test the conversion of bounds to scipy format raises ValueError."""
    # Test case 1: Finite bounds
    xmin = np.array([-10, 0])
    xmax = np.array([10, 5])
    scipy_bounds = optscipy.convert_bounds_to_scipy(xmin, xmax, requires_finite_bounds=True)
    assert scipy_bounds == [(-10, 10), (0, 5)]

    # Test case 2: Mixed bounds
    xmin = np.array([-10, 0])
    xmax = np.array([10, hugeval])
    with pytest.raises(ValueError, match="The scipy function requires finite bounds"):
        optscipy.convert_bounds_to_scipy(xmin, xmax, requires_finite_bounds=True)


def test_optimizers_listed_in_ui():
    """Check that scipy optimizers are listed in the UI."""
    listed_optimizers = ui.list_methods()

    for opt in ['scipy_basinhopping',
                'scipy_differentialevolution',
                'scipy_direct',
                'scipy_dualannealing',
                'scipy_minimize',
                'scipy_shgo',
                ]:
        assert opt in listed_optimizers


def test_bounds_correctly_converted(recwarn):
    """
    This test uses both a direct and an indirect way to check that
    bounds are set to `None` if they are really just numerical max/min
    of a floating point number.

    We perform a fit with the `CG` method, which does not support bounds and will
    raise a warning if any bounds are set. This tests that the bounds are actually
    passed in to scipy correctly. We also check the recorded bounds in the
    extra_output field.
    """
    dat = Data1D('dat', np.arange(10), np.array([2,3,2,3,1,0,1,3,5,2]))
    bkg = Const1D(name='bkg')
    bkg.c0 = 3.4
    line = Delta1D(name='line')
    line.pos = 5.
    line.pos.frozen = True
    line.ampl = -1.2
    model = bkg + line

    method = optmethods.Scipy_Minimize()
    method.method = 'CG'
    fitter = Fit(data=dat, model=model, stat=CStat(), method=method)
    fit_res = fitter.fit(record_steps=True)

    assert len(recwarn) == 0
    assert fit_res.extra_output['input_bounds'] is None

    # And while we are at it, check a few other things.
    # These are regression test and might change if scipy changes the
    # default behavior of the CG method.
    assert fit_res.succeeded
    assert fit_res.nfev == 34
    assert np.allclose(fit_res.parvals, (2.444444465983341, -3.0014523573812197))


def test_instance_name():
    opt = optmethods.Scipy_Minimize()
    assert opt.name == 'scipy.optimize.minimize'

    opt = optmethods.Scipy_Minimize(name='my_name')
    assert opt.name == 'my_name'


def test_scipy_default_config():
    """Test that the kwargs that a automatically filled are not
    in the list of default config options."""
    opt = optmethods.Scipy_Minimize()
    assert opt.default_config == {'method': None, 'tol': None,
                                  'callback': None, 'options': None}


@pytest.mark.parametrize(
    "optimizer,works_unconstrained",
    [(optmethods.Scipy_Minimize, True),
     (optmethods.Scipy_Basinhopping, True),
     (optmethods.Scipy_DifferentialEvolution, False),
     (optmethods.Scipy_DualAnnealing, False),
     (optmethods.Scipy_Shgo, False),
     (optmethods.Scipy_Direct, False),
     ],
)
def test_scipy_unconstrainted(optimizer, works_unconstrained):
    """Test the scipy minimize interface."""
    model = Const1D(name='const')

    opt = optimizer()
    fit = Fit(data=data, model=model, stat=Chi2(), method=opt)

    if works_unconstrained:
        # Perform the fit
        result = fit.fit()

        # Check the result
        assert result.succeeded
        assert result.parvals[0] == pytest.approx(40.75)
    else:
        with pytest.raises(ValueError, match="The scipy function requires finite bounds"):
            # Perform the fit
            fit.fit()



@pytest.mark.parametrize(
    "optimizer,rtol",
    [(optmethods.Scipy_Minimize, 1e-5),
     (optmethods.Scipy_Basinhopping, 1e-5),
     (optmethods.Scipy_DifferentialEvolution, 1e-5),
     (optmethods.Scipy_DualAnnealing, 1e-5),
     (optmethods.Scipy_Shgo, 1e-5),
     (optmethods.Scipy_Direct, 0.001),
     ],
)
def test_scipy_constrained(optimizer, rtol):
    """Test the scipy minimize interface when all parameters have limits.

    We fit a simple one-component model with just three parameters, where
    every fit can succeed with default parameters for the optimizer.
    """
    gauss = Gauss1D(name='gauss')

    gauss.pos.min = 0
    gauss.pos.max = 20
    gauss.fwhm.min = 0.1
    gauss.fwhm.max = 10
    gauss.ampl.min = 0
    gauss.ampl.max = 100

    opt = optimizer()
    fit = Fit(data=data, model=gauss, stat=Chi2(), method=opt)

    # Perform the fit
    result = fit.fit()

    # Regression tests, not from first principles
    # However, it's a good sign when all algorithms give the same result.

    assert gauss.fwhm.val == pytest.approx(7.05428, rel=rtol)
    assert gauss.pos.val == pytest.approx(4.92955, rel=rtol)
    assert gauss.ampl.val == pytest.approx(59.2477, rel=rtol)

    if optimizer == optmethods.Scipy_Direct:
        # Brute gets the right answer, but runs into the maxiterlimit
        assert result.succeeded is False
    else:
        assert result.succeeded is True
