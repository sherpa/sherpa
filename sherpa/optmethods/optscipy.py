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
"""Interface to Scipy optimization methods.

This module contains classes that wrap the optimization functions in
`scipy.optimize` to match the calling signature and return values
to the Sherpa interface.

If `scipy <https://scipy.org>`_ is installed, classes are created automatically and
can be used in the same way as other optimizers in Sherpa.
The most versatile function is `scipy.optimize.minimize`, wrapped into the
`Scipy_Minimize` class.
`scipy.optimize.minimize` is itself a wrapper around several different
optimization algorithms. Which one is used by default depends on the bounds
places on the parameter values of the model to be fit.

`scipy.optimize` also contains several global optimizers that aim to explore
the parameter space more fully. Most of these will only work if meaningful
limits are placed in the parameters, see
`the scipy docs for global optimizers <https://docs.scipy.org/doc/scipy/tutorial/optimize.html#global-optimization>`_
for details.
"""
from collections.abc import Callable
import functools
import inspect
from typing import Any


import numpy as np
from sherpa.optmethods.buildin import OptMethod
from sherpa.utils.types import (ArrayType,
                                OptFunc, OptReturn,
                                StatFunc)
from sherpa.models.parameter import hugeval

# Will be filled out programmatically at the bottom of this file
__all__ = []


def convert_bounds_to_scipy(parmins: ArrayType,
                            parmaxes: ArrayType,
                            requires_finite_bounds: bool) -> list[tuple[float | None, float | None]] | None:
    """Convert parmins and parmaxes to format used by `scipy.optimize` for bounds.

    Parameters
    ----------
    parmins : ArrayType
        The minimum values for the parameters.
    parmaxes : ArrayType
        The maximum values for the parameters.
    requires_finite_bounds : bool
        If True, the scipy function requires finite bounds.

    Returns
    -------
    list[tuple[float | None, float | None]] | None
        The bounds in the format used by scipy.optimize, or None if
        there are no bounds.
    """
    if np.allclose(parmins, -hugeval) and np.allclose(parmaxes, hugeval):
        bounds = None
    else:
        bounds = [(None if pmin == -hugeval else pmin,
                 None if pmax == hugeval else pmax)
                for pmin, pmax in zip(parmins, parmaxes)]
    if requires_finite_bounds:
        if None in np.array(bounds):
            raise ValueError("The scipy function requires finite bounds, but "
                             "the Sherpa model has some bounds set to HUGEVAL.")
    return bounds


# At first sight, this looks like it should be a decorator, but instead
# we make it a higher-order function, because decorated functions cannot
# be pickled.
# https://pythonicthoughtssnippets.github.io/2020/08/09/PTS13-rethinking-python-decorators.html
def wrap_scipy_fcn(func: Callable,
                   requires_finite_bounds: bool,
                   stat: StatFunc,
                   x0: np.ndarray,
                   xmin: np.ndarray,
                   xmax: np.ndarray,
                   **kwargs) -> OptReturn:
    """Wrap a function in scipy.optimize to the Sherpa interface.
    """
    sig = inspect.signature(func)

    def stat_wrapper(x):
        # The function is called with the parameters
        # and returns the statistic and per-bin values
        return stat(x)[0]

    converted_args: dict = {}
    if 'x0' in sig.parameters.keys():
        converted_args['x0'] = x0
    if 'bounds' in sig.parameters.keys():
        converted_args['bounds'] = convert_bounds_to_scipy(xmin, xmax,
                                                            requires_finite_bounds)

    result = func(stat_wrapper, **converted_args, **kwargs)
    for arg in ['bounds', 'ranges']:
        if arg in converted_args:
            result[f'input_{arg}'] = converted_args[arg]
    return (result.success, result.x, result.fun, result.message, result)


SCIPY_KEYWORDS_NOT_APPLICABLE = ['fun', 'func',
                                 'args', 'x0', 'bounds', 'ranges',
                                 'jac', 'hess', 'hessp',
                                 'constraints',
                                 'integrality', 'vectorized',
                                 ]
'''List of keywords NOT to expose when wrapping scipy optimization functions.

This is a list of keywords in the signature of functions in scipy.optimize that we do not
want to expose, either because the interface converts Sherpa input to them
automatically, or because they are not applicable to the Sherpa interface.
'''


class ScipyBase(OptMethod):
    """Base class for wrapping scipy optimization functions.
    This class wraps that function to match the calling signature and
    return values to the Sherpa interface.
    """

    _scipy_func: Callable
    """Optimization function in scipy

    This class wraps that function to match the calling signature and
    return values to the Sherpa interface.
    """

    _requires_finite_bounds: bool
    """If True, the scipy function requires finite bounds on the parameters.

    This is usually the case for global optimizers.
    """
    def _get_default_config(self) -> dict[str, Any]:
        sig = inspect.signature(self._scipy_func)
        return {p.name: p.default for p in sig.parameters.values()
                if p.kind == p.POSITIONAL_OR_KEYWORD and
                 p.name not in SCIPY_KEYWORDS_NOT_APPLICABLE}

    default_config = property(_get_default_config,
                              doc='The default settings for the optimiser.')

    def __init__(self, name : str | None = None) -> None:
         super().__init__(name=f'scipy.optimize.{self._scipy_func.__name__}' if name is None else name,
                         optfunc=functools.partial(wrap_scipy_fcn,
                                                   self._scipy_func,
                                                   self._requires_finite_bounds),
                        )


try:
    from scipy import optimize

    class Scipy_Minimize(ScipyBase):
        """Optimizer using `scipy.optimize.minimize`.

        See the
        `Scipy User Guide for minimzation <https://docs.scipy.org/doc/scipy/tutorial/optimize.html#local-minimization-of-multivariate-scalar-functions-minimize>`_ for details.
        for details on how `scipy.optimize.minimize` works, a short summary
        is below.

        `scipy.optimize.minimize` implements several different algorithms
        than can be chosen with the `method` keyword. Each method has its own
        advantages and drawbacks, e.g. some only work with finite bounds or
        require the function to be smooth. Each method has its own set of parameters
        that controls the optimization process, e.g. for the step size or maximum
        number of iterations.

        See the `scipy.optimize.minimize` documentation for a full list of methods
        and their parameters.
        By default, `scipy.optimize.minimize` will choose a method that this
        appropriate for the given input, but since Sherpa does not supply
        analytic gradients or Hessians, some methods will fail if explicitly
        selected.

        Sherpa will automatically convert statistics functions, input values,
        parameter limits etc. to the format required by the scipy function.
        The following attributes of `scipy.optimize.minimize` can be set
        as attributes of this class.

        Attributes
        ----------
        method : str or callable
            The optimization method to use. See the `scipy.optimize.minimize`
            documentation for a list of available methods.
        tol : float
            Tolerance for termination.
        options : dict
            A dictionary of solver options.
        callback : callable
            A callable called after each iteration.

        Example
        -------

        We begin with a simple example of fitting a Gaussian and a constant
        to some data, using all the default options for the optimizer.

        >>> import numpy as np
        >>> from sherpa.data import Data1D
        >>> from sherpa.models import Const1D, Gauss1D
        >>> from sherpa.stats import Chi2
        >>> from sherpa.optmethods import Scipy_Minimize
        >>> from sherpa.fit import Fit
        >>> data = Data1D('data1', x=np.arange(10),
        ...              y=[34.5, 23.4, 22.3, 45.6, 56.7, 67.8, 58.9, 43.0, 30.1, 25.2],
        ...              staterror=5 * np.ones(10))
        >>> const = Const1D(name='const')
        >>> gauss = Gauss1D(name='gauss')
        >>> gauss.pos = 5.0  # start not too far from the data
        >>> model = const + gauss
        >>> optimizer = Scipy_Minimize()
        >>> fit = Fit(data=data, model=model, stat=Chi2(), method=optimizer)
        >>> result = fit.fit()
        >>> print(result)
        datasets       = None
        itermethodname = none
        methodname     = scipy_minimize
        statname       = chi2
        succeeded      = True
        parnames       = ('const.c0', 'gauss.fwhm', 'gauss.pos', 'gauss.ampl')
        parvals        = (26.06715425374017, 3.285588554348949, 5.0425169187544565, 42.00624527682233)
        statval        = 6.879954385987281
        istatval       = 700.2311623262916
        dstatval       = 693.3512079403043
        numpoints      = 10
        dof            = 6
        qval           = 0.33209185997701957
        rstat          = 1.1466590643312136
        message        = CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH
        nfev           = 195

        We can change the method used in the optimization and set some of the
        options for the optimizer, such as the goal for the precision and the
        maximum number of function evaluations before the fit stops.

        >>> optimizer.method = 'TNC'
        >>> optimizer.options = {'maxfun': 200, 'xtol': 1e-4}
        >>> model.reset()
        >>> print(fit.fit())
        datasets       = None
        ...
        rstat          = 1.1466659346262622
        message        = Converged (|x_n-x_(n-1)| ~= 0)
        nfev           = 725
        """


        _scipy_func = staticmethod(optimize.minimize)
        _requires_finite_bounds = False


    class Scipy_Basinhopping(ScipyBase):
        """Optimizer using `scipy.optimize.basinhopping`.

        Basin-hopping is a two-phase method that combines a global stepping
        algorithm with local minimization at each step. Designed to mimic the
        natural process of energy minimization of clusters of atoms, it works
        well for similar problems with “funnel-like, but rugged” energy
        landscapes.

        See the `scipy.optimize.basinhopping` documentation for details of
        all parameters.
        Sherpa will automatically convert statistics functions, input values,
        parameter limits etc. to the format required by the scipy function.
        The following attributes can be set as attributes of this class.

        Attributes
        ----------
        niter : int
            Number of iterations to perform.
        T : float
            The “temperature” parameter for the acceptance or rejection
            criterion. Higher “temperatures” mean that larger jumps in function
            value will be accepted. For best results T should be comparable to
            the separation (in function value) between local minima.
        step_size : float
            Maximum step size for use in the random displacement.
        minimizer_kwargs : dict
            A dictionary of options to pass to the local minimizer.
        take_step : callable
            Replace the default step-taking routine with this routine.
        accept_test : callable
            Define a test which will be used to judge whether to accept the
            step. This will be used in addition to the Metropolis test based
            on “temperature” T.
        callback : callable
            A callable called after each iteration.
        niter_success : int
            Stop the run if the global minimum candidate remains the same for this
            number of iterations.
        rng : {None, int, `numpy.random.Generator`}
            Random number generator instance or seed.
        target_acceptance_rate : float
            The target acceptance rate for the step acceptance test.
        stepwise_factor : float
            The stepsize is multiplied or divided by this stepwise factor upon each
            update. Range is (0, 1). Default is 0.9.
        """
        _scipy_func = staticmethod(optimize.basinhopping)
        _requires_finite_bounds = False


    class Scipy_DifferentialEvolution(ScipyBase):
        """Optimizer using `scipy.optimize.differential_evolution`.

        The differential evolution method is stochastic and can search
        large areas of candidate space at the cost of running a longer than
        typical gradient-based methods.

        See the `scipy.optimize.differential_evolution` documentation for details
        of all parameters.
        Sherpa will automatically convert statistics functions, input values,
        parameter limits etc. to the format required by the scipy function.
        The following attributes can be set as attributes of this class.

        Attributes
        ----------
        strategy : str
            The differential evolution strategy to use. See the `scipy.optimize.differential_evolution`
            documentation for a list of available strategies.
        maxiter : int
            The maximum number of generations to be evaluated.
        popsize : int
            A multiplier for setting the total population size.
        tol : float
            Relative tolerance for convergence.
        mutation : float or tuple
            The mutation constant.
        recombination : float
            The recombination constant.
        rng : {None, int, `numpy.random.Generator`}
            Random number generator instance or seed.
        disp : bool
            Set to True to print convergence messages.
        callback : callable
            A callable called after each iteration.
        polish : bool
            If True, the best solution is refined by a local optimizer.
        init : str or array_like
            The method used to initialize the population, see the
            `scipy.optimize.differential_evolution`
            documentation for a full list.
        atol : float
            Absolute tolerance for convergence.
        updating : {'immediate', 'deferred'}
            Whether to update the population immediately or only
            once per generation.
        workers : int or map-like callable
            The number of workers to use for parallelization.
        """
        _scipy_func = staticmethod(optimize.differential_evolution)
        _requires_finite_bounds = True


    class Scipy_DualAnnealing(ScipyBase):
        """Optimizer using `scipy.optimize.dual_annealing`.

        See the `scipy.optimize.dual_annealing` documentation for details
        of all parameters.
        Sherpa will automatically convert statistics functions, input values,
        parameter limits etc. to the format required by the scipy function.
        The following attributes can be set as attributes of this class.

        Attributes
        ----------
        maxiter : int
            The maximum number of iterations.
        minimizer_kwargs : dict
            A dictionary of options to pass to the local minimizer.
        initial_temp : float
            A higher value allows dual_annealing to escape
            local minima that it is trapped in.
        restart_temp_ratio : float
            During the annealing process, temperature is decreasing, when it
            reaches initial_temp * restart_temp_ratio, the reannealing process
            is triggered. Default value of the ratio is 2e-5. Range is (0, 1).
        visit : float
            Parameter for visiting distribution. Default value is 2.62.
            Higher values allow jumps to more distant regions.
            The value range is (1, 3].
        accept : float
            Control for the acceptance probability. Default value is -5.0
            with a range (-1e4, -5].
        maxfun : int
            Soft limit for the number of objective function calls.
            Default value is 1e7.
        rng : {None, int, `numpy.random.Generator`}
            Random number generator instance or seed.
        no_local_search : bool
            If True, the local search is not performed.
        callback : callable
            A callable called after each iteration.
        """
        _scipy_func = staticmethod(optimize.dual_annealing)
        _requires_finite_bounds = True


    class Scipy_Shgo(ScipyBase):
        """Optimizer using `scipy.optimize.shgo`.

        The `scipy.optimize.shgo` function implements the simplicial homology
        global optimization algorithm. It is a global optimization algorithm
        that is designed to find the global minimum.

        Sherpa will automatically convert statistics functions, input values,
        parameter limits etc. to the format required by the scipy function and those
        cannot be set by the user. Below is a list of attributes most likely to be
        be of interest for Sherpa users, see `scipy.optimize.shgo`
        documentation for details of all parameters.

        Attributes
        ----------
        n : int
            Number of sampling points used in the construction of the simplicial complex.
        iters : int
            Number of iterations used in the construction of the simplicial complex.
            Default is 1.
        callback : callable
            A callable called after each iteration.
        minimizer_kwargs : dict
            A dictionary of options to pass to the local minimizer.
        options : dict
            A dictionary of solver options. See the `scipy.optimize.shgo`
            documentation for a full list.
        sampling_method : str or function
            Current built in sampling method options are `halton`, `sobol`
            and `simplicial`.
            See the `scipy.optimize.shgo` documentation for details.
        workers : int or map-like callable
            The number of workers to use for parallelization.
        """
        _scipy_func = staticmethod(optimize.shgo)
        _requires_finite_bounds = True


    class Scipy_Direct(ScipyBase):
        """Optimizer using `scipy.optimize.direct`.

        The `scipy.optimize.direct` function implements the DIRECT algorithm
        for global optimization.

        Sherpa will automatically convert statistics functions, input values,
        parameter limits etc. to the format required by the scipy function and those
        cannot be set by the user. Below is a list of attributes most likely to be
        be of interest for Sherpa users, see `scipy.optimize.direct`
        documentation for details of all parameters.

        Attributes
        ----------
        eps : float
            Minimal required difference between steps. Default is 1e-4.
        maxfun : int or None
            Approximate upper bound on objective function evaluations.
        maxiter : int
            Maximum number of iterations.
        locally_biased : bool
            If True (default), use the locally biased variant of the algorithm
            known as DIRECT_L.
        vol_tol : float
            Terminate the optimization once the volume of the hyperrectangle
            containing the lowest function value is smaller than vol_tol of the
            complete search space. Must lie between 0 and 1. Default is 1e-16.
        len_tol : float
            Allowed range is 0 to 1. Default is 1e-6.
        """
        _scipy_func = staticmethod(optimize.direct)
        _requires_finite_bounds = True


    __all__.append('Scipy_Minimize')
    __all__.append('Scipy_Basinhopping')
    __all__.append('Scipy_DifferentialEvolution')
    __all__.append('Scipy_DualAnnealing')
    __all__.append('Scipy_Shgo')
    __all__.append('Scipy_Direct')


except ImportError:
    # scipy is not available, so we cannot create the classes
    # that wrap the scipy functions
    pass
