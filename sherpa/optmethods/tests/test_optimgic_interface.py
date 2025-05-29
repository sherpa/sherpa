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
"""Test interface to Optimagic optimization methods.

These tests are intentionally sparse, as the optimizers themselves are
already tested in the optimagic library. Here, we only want to test that the
interface works and parameters are passed to optimagic correctly.

optimagic interfaces to optimizers in many different libraries, but of those
only scipy is a required dependency. Thus, all tests are run with scipy
optimizers and we leave it to the optimagic library to test that their interface
to the other backends works - after all the whole point of optimagic is to offer
the same interface to all of them.
"""

import pytest
pytest.importorskip("optimagic", reason="optimagic not installed")

import optimagic as om
import numpy as np

from sherpa import ui
from sherpa import optmethods
from sherpa.fit import Fit
from sherpa.data import Data1D
from sherpa.models import Const1D, Gauss1D
from sherpa.stats import Chi2
from sherpa.optmethods import Optimagic


data = Data1D('data1', x=np.arange(10),
              y=[34.5, 23.4, 22.3, 45.6, 56.7, 67.8, 58.9, 43.0, 30.1, 25.2],
              staterror=5 * np.ones(10))

def test_optimizers_listed_in_ui():
    """Check that scipy optimizers are listed in the UI."""
    listed_optimizers = ui.list_methods()

    assert 'optimagic' in listed_optimizers


def test_instance_name():
    opt = optmethods.Optimagic()
    assert opt.name == 'optimagic.minimize'

    opt = optmethods.Optimagic(name='my_name')
    assert opt.name == 'my_name'


@pytest.mark.parametrize(
    "optimizer,works_unconstrained",
    [('scipy_lbfgsb', True),
     ('scipy_differential_evolution', False),
     ('pounders', True),  # A Least-Squares algorithm
     ],
)
def test_optimagic_unconstrainted(optimizer, works_unconstrained):
    """Test the scipy minimize interface."""
    model = Const1D(name='const')

    opt = Optimagic()
    opt.algorithm = optimizer
    fit = Fit(data=data, model=model, stat=Chi2(), method=opt)

    if works_unconstrained:
        # Perform the fit
        result = fit.fit()

        # Check the result
        assert result.succeeded
        assert result.parvals[0] == pytest.approx(40.75)
    else:
        with pytest.raises(ValueError, match="The optimizer requires finite bounds"):
            fit.fit()


@pytest.mark.parametrize(
    "optimizer,rtol",
    [('scipy_lbfgsb', 1e-5),  # algorithm uses scalar stats values
     (om.algorithms.ScipyLBFGSB, 1e-5),  # call with class
     (om.algorithms.ScipyLBFGSB(), 1e-5),  # call with instance
     ('scipy_differential_evolution', 1e-5),  # a global optimizer
     ("pounders", 1e-5),  # algorithm uses least-squares structure
     ],
)
def test_optimagic_constrained(optimizer, rtol):
    """Test the scipy minimize interface when all parameters have limits.

    We fit a simple one-component model with just three parameters, where
    every fit can succeed with default parameters for the optimizer.
    """
    gauss = Gauss1D(name='gauss')

    gauss.pos.min = 0
    gauss.pos.max = 20
    gauss.fwhm.min = 0.1
    gauss.fwhm.max = 20
    gauss.ampl.min = 0
    gauss.ampl.max = 100

    opt = Optimagic(algorithm=optimizer)
    fit = Fit(data=data, model=gauss, stat=Chi2(), method=opt)

    # Perform the fit
    result = fit.fit()

    # Regression tests, not from first principles
    # However, it's a good sign when all algorithms give the same result.

    assert gauss.fwhm.val == pytest.approx(7.05428, rel=rtol)
    assert gauss.pos.val == pytest.approx(4.92955, rel=rtol)
    assert gauss.ampl.val == pytest.approx(59.2477, rel=rtol)

    assert result.succeeded is True


def test_optimagic_set_options():
    """Test that we can pass options to the optimagic optimizer.
    
    We need ot choose an option that is (a) not the default and (b)
    has a noticeable impact on results so that we can be sure it was
    actually used.

    To some extent this is already tested in the test_optimagic_constrained
    and test_optimagic_unconstrainted tests, where we use the "algorithm"
    option, but we want to test that we can set other options too.
    """

    gauss = Gauss1D(name='gauss')

    gauss.pos.min = 0
    gauss.pos.max = 20
    gauss.fwhm.min = 0.1
    gauss.fwhm.max = 10
    gauss.ampl.min = 0
    gauss.ampl.max = 100

    opt = Optimagic()
    opt.algorithm = "scipy_bfgs"
    fit = Fit(data=data, model=gauss, stat=Chi2(), method=opt)

    # Perform the fit
    result = fit.fit()

    # Regression tests, not from first principles
    assert gauss.fwhm.val == pytest.approx(7.05428)
    assert gauss.pos.val == pytest.approx(4.92955)
    assert gauss.ampl.val == pytest.approx(59.2477)

    assert result.succeeded is True

    gauss.reset()
    opt.algo_options = {'stopping.maxiter': 2, 'norm': 1000}
    result2 = fit.fit()
    # After just 2 steps, it's not converged yet, so the answers are different
    assert result2.nfev is None
    assert gauss.fwhm.val == 0.1
    assert gauss.pos.val == 0.

    assert result2.succeeded is False


def test_optimagic_invalid_algorithm():
    """Test that we get optimagic error for non-exiting algorithm.
    """

    gauss = Gauss1D(name='gauss')

    opt = Optimagic()
    opt.algorithm = "does_not_exist"
    fit = Fit(data=data, model=gauss, stat=Chi2(), method=opt)

    # Perform the fit
    with pytest.raises(AttributeError, match="optimagic does not have an algorithm named does_not_exist."):   
        result = fit.fit()


def test_optimagic_setattr():
    """Check the call-through-to-config option works"""

    opt = Optimagic()
    opt.algorithm = "scipy_lbfgsb"
    # Unlike test_optmethod_getattr, only check one config value
    oldval = opt.config["algorithm"]
    assert opt.algorithm == pytest.approx(oldval)

    newval = "scipy_differential_evolution"
    opt.algorithm= newval
    assert opt.config["algorithm"] == pytest.approx(newval)
    assert opt.algorithm == pytest.approx(newval)

    # just to check that ftol doesn't happen to match the new value,
    # which would mean changing newval
    #
    assert oldval != pytest.approx(newval)
