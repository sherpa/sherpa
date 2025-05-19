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
"""Interface to optimagic optimization methods.

This module contains classes that wrap the optimization functions in
`optimagic <https://estimagic.org>` to match the calling signature and 
return values to the Sherpa interface.

If `optimagic <https://estimagic.org>`_ is installed, classes are created automatically and
can be used in the same way as other optimizers in Sherpa.

Optimagic supports well over 50 different algorithms
----------------------------------------------------
Optimagic supports several different optimizers from about a dozen different
libraries, including scipy, nlopt, pygmo2, and others; most if these dependencies
need to be installed separately; they will only be available to you if the
respective library is installed in your Python environment.

While optimagic aims to provide a consistent interface to all of these libraries,
the options that can be set for those optimizers depends on the algorithm,
e.g. in some cases one may set the maximum number of iterations, while in others
one may define the direction of the optimization. The 
`optimagic documentation <https://optimagic.readthedocs.io/en/latest/algorithms.html>`_
lists all available algorithms and their options.

There is no default algorithm, so you always need to specify the 
algorithm to use before you can use the optimizer in a fit:

    >>> from sherpa.optmethods import Optimagic
    >>> opt = Optimagic(name='my_cool_optimizer')
    >>> opt.algorithm = 'scipy_ls_dogbox'

Setting additional algorithm options is supported, but not required
(all options have defaults):

    >>> opt.algo_options = {'stopping.maxfun': 1000}

Most inputs for optimagic are automatically filled by Sherpa
------------------------------------------------------------
This module wraps the `optimagic.minimize` function in such a way that
all the information that Sherpa already has about the model, the data, and the
statistic is passed to the optimizer. Sherpa will automatically pick the correct
format of the statistic function, convert the parameter minimum and maximum
to the format for optimagic bounds etc.

Caveats to watch out for
-------------------------
Optimagic contains both local and global optimizers (see :ref:`fit-strategies-setminmax`
for a discussion). Global optimizers need all parameters to be bounded
between finite values. 

Sherpa never provides analytical gradients for any model 
(because many Sherpa models are not differentiable, e.g. use tabulated values).
While optimagic can use numerical gradients, this means that much of the 
efficiency gain for typical gradient-based optimizers is lost.

Optimagic is designed to be fault tolerant if, e.g. some bins return NaN
values, but how exactly that impacts the optimization depends on the
algorithm used. If errors occur, the error from optimagic will be passed to
the user, which may differ from typical Sherpa error messages.

What statistic functions are supported?
---------------------------------------
Most statistic functions return a single scalar value that is minimized
by the optimizer. However, some optimizers in optimagic are designed to
work with explicitly with least-squares or likelihood functions and can converge faster
if the optimizer is given the full array of per-bin values.
Sherpa will automatically detect if the optimizer is one of those and return
the statistic values in the correct format. However, **it is up to user** to make
sure that statistic and the optimizer are compatible. For example, one *can* feed
a `~sherpa.stats.LeastSq` statistic into the `optimagic.algos.bhhh` optimizer
(which expects a likelihood), but the fit might take much longer or may not
converge at all and the results could be wrong.

Further diagnostics
-------------------
The last element of the optimizer output is a dictionary which holds the output
from optimagic in an `optimagic.optimization.optimize_result.OptimizeResult` object.
When used in a fit, this dictionary is stored in the `extra_output` field
of the fit result: `result.extra_output['optimagic result']`.
These output values can be used for example to visualize the `optimization history
<https://optimagic.readthedocs.io/en/latest/how_to/how_to_visualize_histories.html>`_
or to `benchmark different optimizers
<https://optimagic.readthedocs.io/en/latest/how_to/how_to_benchmarking.html>`_
as described in the optimagic documentation.




"""
import inspect
from typing import Any


import numpy as np
from sherpa.optmethods.buildin import OptMethod
from sherpa.utils.types import (OptReturn,
                                StatFunc)
from .optscipy import convert_bounds_to_scipy

# Will be filled out programmatically at the bottom of this file
__all__ = []


APPLICABLE_OPTIMAGIC_MINIMIZE_KEYWORDS = [
    'algorithm', 'algo_options', 'numdiff_options', 
    'error_handling', 'error_penalty', 'scaling',
    ]
'''Keyword arguments to `optimagic.minimize` that are exposed in Sherpa.'''


try:
    import optimagic as om
    from optimagic.typing import AggregationLevel

    # At first sight, this looks like it should be a decorator, but instead
    # we make it a higher-order function, because decorated functions cannot
    # be pickled.
    # https://pythonicthoughtssnippets.github.io/2020/08/09/PTS13-rethinking-python-decorators.html
    def wrap_optimagic_minimize(
                    stat: StatFunc,
                    x0: np.ndarray,
                    xmin: np.ndarray,
                    xmax: np.ndarray,
                    **kwargs) -> OptReturn:
        """Wrap `optimagic.minimize` to the Sherpa interface.
        """

        bounds = None
        if hasattr(kwargs['algorithm'], 'algo_info'):
            algo_info = kwargs['algorithm'].algo_info
        else:
            try:
                algorithm = getattr(om.algos, kwargs['algorithm'])
            except AttributeError:
                raise AttributeError('optimagic does not have an algorithm '
                                     f'named {kwargs["algorithm"]}. ')
            else:
                algo_info = algorithm.algo_info

        bounds = None
        if algo_info.supports_bounds:
            bounds = convert_bounds_to_scipy(xmin, xmax, algo_info.is_global)

        match algo_info.solver_type:
            case AggregationLevel.SCALAR:
                def stat_wrapper(x):
                    return stat(x)[0]
            case AggregationLevel.LEAST_SQUARES:
                @om.mark.least_squares
                def stat_wrapper(x):
                    return stat(x)[1]
            case AggregationLevel.LIKELIHOOD:
                @om.mark.likelihood
                def stat_wrapper(x):
                    return -stat(x)[1]
            case _:
                # As of May 2025, there is no other option in optimagic
                # but to future-proof, raise an error if new options are added.
                raise NotImplementedError
        

        result = om.minimize(stat_wrapper, 
                             x0,
                             bounds=bounds,
                             **kwargs)
        
        out = (# I don't think we'll ever get None as a return value, but Mypy complains
               # that the type declaration in optimagic is "bool | None", while we want "bool".
               # None certainly isn't a valid success.
               result.success if result.success else False, 
               result.params, 
               result.fun, 
               # The message in optimagic is "str | None", while Sherpa expects "str".
               result.message if result.message else '', 
               {'optimagic result': result, 'input_bounds': bounds})
        return out


    class Optimagic(OptMethod):
        """Optimizer using the optimagic library.

        This class is a wrapper around the optimagic optimization methods.
        It provides a common interface for all optimizers in the optimagic library.
        See `the very useful and extensive optimagic documentation 
        <https://optimagic.readthedocs.io/en/latest/>`_ for general information
        on optimizations and all the algorithms available.

        Sherpa will automatically convert statistics functions, input values,
        parameter limits etc. to the format required by the scipy function.
        The following attributes of `optimagic.minimize` can be set
        as attributes of this class.   
        
        Some optimizers in optimagic are designed tobwork with explicitly with 
        least-squares or likelihood functions and can converge faster
        if the optimizer is given the full array of per-bin values.
        Sherpa will automatically detect if the optimizer is one of those and return
        the statistic values in the correct format. However, **it is up to user** to make
        sure that statistic and the optimizer are compatible. For example, one *can* feed
        a `~sherpa.stats.LeastSq` statistic into the `optimagic.algos.bhhh` optimizer
        (which expects a likelihood), but the fit might take much longer or may not
        converge at all and the results could be wrong.

        Attributes
        ----------
        algorithm :
            The optimization algorithm to use, see 
            https://optimagic.readthedocs.io/en/latest/algorithms.html
            for list of optimizers and their dependencies.
            Can be a string, subclass of `optimagic.algorithm.Algorithm` or an instance
            of a subclass of `optimagic.algorithm.Algorithm`.
        algo_options : dict
            A dictionary of options for the algorithm. The options depend on
            the algorithm used. See https://optimagic.readthedocs.io/en/latest/algorithms.html
            for a list of algorithms and their options.
        numdiff_options : dict or `optimagic.differentiation.numdiff_options.NumdiffOptions`
            Options for numerical differentiation.
        error_handling : str
            Set what to do if an error (e.g. a NaN) occurs during the optimization.
            Default is "raise", which will stop the optimization. Alternatively,
            this can be set to "continue", which will employ a penalty function
            to guide the optimizer back to values that can be evaluated.
        error_penalty : dict
            Dictionary with keys "slope" and "constant" that influences the magnitude
            of the penalty values. Both number should be positive.
        scaling
            If None or False, the parameter space is not rescaled. 
            If True, a heuristic is used to improve the conditioning of the
            optimization problem. To choose which heuristic is used and to customize
            the scaling, provide a dictionary or an instance of 
            `optimagic.parameters.scaling.ScalingOptions`.
        
        Example
        -------
        In this example, we fit a Gaussian and a constant to some data:

        >>> import numpy as np
        >>> from sherpa.data import Data1D
        >>> from sherpa.models import Const1D, Gauss1D
        >>> from sherpa.stats import Chi2
        >>> from sherpa.optmethods import Optimagic
        >>> from sherpa.fit import Fit
        >>> data = Data1D('data1', x=np.arange(10),
        ...              y=[34.5, 23.4, 22.3, 45.6, 56.7, 67.8, 58.9, 43.0, 30.1, 25.2],
        ...              staterror=5 * np.ones(10))
        >>> const = Const1D(name='const')
        >>> gauss = Gauss1D(name='gauss')
        >>> gauss.pos = 5.0  # start not too far from the data
        >>> model = const + gauss

        When we create the optimizer, we always have to set the algorithm
        before we can use the optimizer in a fit. 
        Optimagic has numerous optional dependencies, but only scipy is required,
        so for this example we will use one of the scipy algorithms.

        >>> optimizer = Optimagic()
        >>> optimizer.algorithm = 'scipy_bfgs'

        >>> fit = Fit(data=data, model=model, stat=Chi2(), method=optimizer)
        >>> result = fit.fit()
        >>> print(result)
        datasets       = None
        itermethodname = none
        methodname     = optimagic
        statname       = chi2
        succeeded      = True
        parnames       = ('const.c0', 'gauss.fwhm', 'gauss.pos', 'gauss.ampl')
        parvals        = (26.067327420836456, 3.285600683850382, 5.042514578467134, 42.005725690386164)
        statval        = 6.879954370635533
        istatval       = 700.2311623262916
        dstatval       = 693.351207955656
        numpoints      = 10
        dof            = 6
        qval           = 0.33209186143330294
        rstat          = 1.1466590617725887
        message        = Optimization terminated successfully.
        nfev           = None

        The last element of the fit output is a dictionary which holds the output
        from optimagic in an `optimagic.optimization.optimize_result.OptimizeResult` object,
        which can provide further information about the optimization process:

        >>> result.extra_output['optimagic result'].convergence_report['five_steps']
        {'relative_criterion_change': 0.00044802722082809587,
        'relative_params_change': 0.0025520339599334517,
        'absolute_criterion_change': 0.003082406836099949,
        'absolute_params_change': 0.041909631395145426}
        """
        def _get_default_config(self) -> dict[str, Any]:
            sig = inspect.signature(om.minimize)
            return {p.name: p.default for p in sig.parameters.values()
                    if p.name in APPLICABLE_OPTIMAGIC_MINIMIZE_KEYWORDS}   
        
        default_config = property(_get_default_config,
                                doc='The default settings for the optimiser.') 
        
        def __init__(self, name : str | None = None) -> None:
            super().__init__(name='optimagic.minimize' if name is None else name,
                             optfunc=wrap_optimagic_minimize,
                             )
            
    __all__.append('Optimagic')

except ImportError:
    # optimagic is not installed
    pass