from __future__ import print_function
#
#  Copyright (C) 2009, 2015, 2016  Smithsonian Astrophysical Observatory
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

from six.moves import zip as izip
from six.moves import xrange
import logging
import os
import signal
from numpy import arange, array, abs, iterable, sqrt, where, \
    ones_like, isnan, isinf, float, float32, finfo, nan, any, zeros, append, \
    int, ones
from sherpa.utils import NoNewAttributesAfterInit, print_fields, erf, igamc, \
    bool_cast, is_in, is_iterable, list_to_open_interval, sao_fcmp
from sherpa.utils.err import FitErr, EstErr, SherpaErr
from sherpa.data import DataSimulFit
from sherpa.estmethods import Covariance, EstNewMin
from sherpa.models import SimulFitModel
from sherpa.optmethods import LevMar, NelderMead
from sherpa.stats import Chi2, Chi2Gehrels, Cash, CStat, Chi2ModVar, LeastSq, \
    Likelihood, WStat, UserStat

warning = logging.getLogger(__name__).warning
info = logging.getLogger(__name__).info

__all__ = ('FitResults', 'ErrorEstResults', 'Fit')


# def get_valid_args(func):
#     valid_args = func.func_code.co_varnames[:func.func_code.co_argcount]
#     kwargs_length = len(func.func_defaults)  # number of keyword arguments
#     valid_kwargs = valid_args[-kwargs_length:]  # because kwargs are last
#     return valid_kwargs


def evaluates_model(func):
    """
    Fit object decorator that runs model startup() and teardown()

    """
    def run(fit, *args, **kwargs):
        fit.model.startup()
        result = func(fit, *args, **kwargs)
        fit.model.teardown()
        return result

    run.__name__ = func.__name__
    run.__doc__ = func.__doc__
    return run


class StatInfoResults(NoNewAttributesAfterInit):

    """A summary of the current statistic value for
    one or more data sets.

    Attributes
    ----------
    name : str
       The name of the data set, or sets.
    ids : sequence of int or str
       The data set ids (it may be a tuple or array) included in the
       results.
    bkg_ids: sequence of int or str, or `None`
       The background data set ids (it may be a tuple or array)
       included in the results, if any.
    statname : str
       The name of the statistic function.
    statval : number
       The statistic value.
    numpoints : int
       The number of bins used in the fits.
    dof: int
       The number of degrees of freedom in the fit (the number of
       bins minus the number of free parameters).
    qval: number or `None`
       The Q-value (probability) that one would observe the reduced
       statistic value, or a larger value, if the assumed model is
       true and the current model parameters are the true parameter
       values. This will be `None` if the value can not be calculated
       with the current statistic (e.g. the Cash statistic).
    rstat: number of `None`
       The reduced statistic value (the `statval` field divided by
       `dof`). This is not calculated for all statistics.

    """

    _fields = ('name', 'ids', 'bkg_ids', 'statname', 'statval',
               'numpoints', 'dof', 'qval', 'rstat')

    def __init__(self, statname, statval, numpoints, model, dof,
                 qval=None, rstat=None):
        self.name = ''
        self.ids = None
        self.bkg_ids = None

        self.statname = statname
        self.statval = statval
        self.numpoints = numpoints
        self.model = model
        self.dof = dof
        self.qval = qval
        self.rstat = rstat

    def __repr__(self):
        return '<Statistic information results instance>'

    def __str__(self):
        return print_fields(self._fields, vars(self))

    def format(self):
        s = ''
        if (self.ids is not None and self.bkg_ids is None):
            if len(self.ids) == 1:
                s = 'Dataset               = %s\n' % str(self.ids[0])
            else:
                s = 'Datasets              = %s\n' % str(self.ids).strip("()")
        elif (self.ids is not None and self.bkg_ids is not None):
            s = 'Background %s in Dataset = %s\n' % (str(self.bkg_ids[0]),
                                                     str(self.ids[0]))
        s += 'Statistic             = %s\n' % self.statname
        s += 'Fit statistic value   = %g\n' % self.statval
        s += 'Data points           = %g\n' % self.numpoints
        s += 'Degrees of freedom    = %g' % self.dof
        if self.qval is not None:
            s += '\nProbability [Q-value] = %g' % self.qval
        if self.rstat is not None:
            s += '\nReduced statistic     = %g' % self.rstat
        return s


class FitResults(NoNewAttributesAfterInit):

    """A summary of the fit results.

    Attributes
    ----------
    datasets : sequence of int or str
       A sequence of the data set ids included in the results.
    itermethodname : str or `None`
       What iterated-fit scheme was used, if any.
    statname : str
       The name of the statistic function.
    succeeded : bool
       Was the fit successful (did it converge)?
    parnames : tuple of str
       the parameter names that were varied in the fit
       (the thawed parameters in the model expression).
    parvals : tuple of number
       The parameter values, in the same order as `parnames`.
    statval : number
       The statistic value after the fit.
    istatval : number
       The statistic value at the start of the fit.
    dstatval : number
       The change in the statistic value (`istatval - statval`).
    numpoints : int
       The number of bins used in the fits.
    dof : int
       The number of degrees of freedom in the fit (the number of
       bins minus the number of free parameters).
    qval : number or `None`
       The Q-value (probability) that one would observe the reduced
       statistic value, or a larger value, if the assumed model is
       true and the current model parameters are the true parameter
       values. This will be `None` if the value can not be calculated
       with the current statistic (e.g.  the Cash statistic).
    rstat : number or `None`
       The reduced statistic value (the `statval` field divided by
       `dof`). This is not calculated for all statistics.
    message : str
       A message about the results of the fit (e.g. if the fit was
       unable to converge). The format and contents depend on the
       optimisation method.
    nfev : int
       The number of model evaluations made during the fit.

    """

    _fields = ('datasets', 'itermethodname', 'methodname', 'statname',
               'succeeded', 'parnames', 'parvals', 'statval', 'istatval',
               'dstatval', 'numpoints', 'dof', 'qval', 'rstat', 'message',
               'nfev')

    def __init__(self, fit, results, init_stat, param_warnings):
        _vals = fit.data.eval_model_to_fit(fit.model)
        _dof = len(_vals) - len(tuple(results[1]))
        _qval = None
        _rstat = None
        _covarerr = results[4].get('covarerr')
        if (isinstance(fit.stat, (CStat, WStat, Chi2)) and
                not isinstance(fit.stat, LeastSq)):
            if _dof > 0 and results[2] >= 0.0:
                _qval = igamc(_dof / 2., results[2] / 2.)
                _rstat = results[2] / _dof
            else:
                _rstat = nan

        self.succeeded = results[0]
        self.parnames = tuple(p.fullname for p in fit.model.pars
                              if not p.frozen)
        self.parvals = tuple(results[1])
        self.istatval = init_stat
        self.statval = results[2]
        self.dstatval = abs(results[2] - init_stat)
        self.numpoints = len(_vals)
        self.dof = _dof
        self.qval = _qval
        self.rstat = _rstat
        self.message = results[3]
        if _covarerr is not None:
            self.covarerr = tuple(_covarerr)
        else:
            self.covarerr = None
        self.nfev = results[4].get('nfev')
        self.extra_output = results[4]
        self.modelvals = _vals
        self.methodname = type(fit.method).__name__.lower()
        self.itermethodname = fit._iterfit.itermethod_opts['name']
        statname = type(fit.stat).__name__.lower()
        if isinstance(fit.stat, Chi2) and not isinstance(fit.stat, LeastSq):
            isSimulFit = isinstance(fit.data, DataSimulFit)
            if isSimulFit:
                is_error_set = [
                    d.staterror is not None for d in fit.data.datasets]
                if all(is_error_set):
                    statname = 'chi2'
            elif fit.data.staterror is not None:
                statname = 'chi2'
        self.statname = statname
        self.datasets = None  # To be filled by calling function
        self.param_warnings = param_warnings
        NoNewAttributesAfterInit.__init__(self)

    def __setstate__(self, state):
        self.__dict__.update(state)

        if 'itermethodname' not in state:
            self.__dict__['itermethodname'] = 'none'

    def __nonzero__(self):
        return self.succeeded

    def __repr__(self):
        return '<Fit results instance>'

    def __str__(self):
        return print_fields(self._fields, vars(self))

    def format(self):
        s = ''
        if (self.datasets is not None):
            if len(self.datasets) == 1:
                s = 'Dataset               = %s\n' % str(self.datasets[0])
            else:
                s = 'Datasets              = %s\n' % str(
                    self.datasets).strip("()")
        if (self.itermethodname is not None and self.itermethodname != 'none'):
            s += 'Iterative Fit Method  = %s\n' % self.itermethodname.capitalize()
        s += 'Method                = %s\n' % self.methodname
        s += 'Statistic             = %s\n' % self.statname
        s += 'Initial fit statistic = %g\n' % self.istatval
        s += 'Final fit statistic   = %g' % self.statval
        if self.nfev is not None:
            s += ' at function evaluation %d' % self.nfev

        s += '\nData points           = %g' % self.numpoints
        s += '\nDegrees of freedom    = %g' % self.dof

        if self.qval is not None:
            s += '\nProbability [Q-value] = %g' % self.qval
        if self.rstat is not None:
            s += '\nReduced statistic     = %g' % self.rstat
        s += '\nChange in statistic   = %g' % self.dstatval

        if self.covarerr is None:
            for name, val in izip(self.parnames, self.parvals):
                s += '\n   %-12s   %-12g' % (name, val)
        else:
            for name, val, covarerr in izip(self.parnames, self.parvals, self.covarerr):
                s += '\n   %-12s   %-12g +/- %-12g' % (name, val, covarerr)

        if self.param_warnings != "":
            s += "\n" + self.param_warnings

        return s


class ErrorEstResults(NoNewAttributesAfterInit):

    _fields = ('datasets', 'methodname', 'iterfitname', 'fitname', 'statname',
               'sigma', 'percent', 'parnames', 'parvals', 'parmins',
               'parmaxes', 'nfits')

    def __init__(self, fit, results, parlist=None):
        if (parlist is None):
            parlist = [p for p in fit.model.pars if not p.frozen]

        from sherpa.estmethods import est_success, est_hardmin, est_hardmax, est_hardminmax

        warning_hmin = "hard minimum hit for parameter "
        warning_hmax = "hard maximum hit for parameter "
        self.datasets = None  # To be set by calling function
        self.methodname = type(fit.estmethod).__name__.lower()
        self.iterfitname = fit._iterfit.itermethod_opts['name']
        self.fitname = type(fit.method).__name__.lower()
        self.statname = type(fit.stat).__name__.lower()
        self.sigma = fit.estmethod.sigma
        self.percent = erf(self.sigma / sqrt(2.0)) * 100.0
        self.parnames = tuple(p.fullname for p in parlist if not p.frozen)
        self.parvals = tuple(p.val for p in parlist if not p.frozen)
        self.parmins = ()
        self.parmaxes = ()
        self.nfits = 0
        success = True
        for i in range(len(parlist)):
            if (results[2][i] != est_success):
                success = False
            if (results[2][i] == est_hardmin or
                    results[2][i] == est_hardminmax):
                self.parmins = self.parmins + (None,)
                warning(warning_hmin + self.parnames[i])
            else:
                self.parmins = self.parmins + (results[0][i],)

            if (results[2][i] == est_hardmax or
                    results[2][i] == est_hardminmax):
                self.parmaxes = self.parmaxes + (None,)
                warning(warning_hmax + self.parnames[i])
            else:
                self.parmaxes = self.parmaxes + (results[1][i],)

        self.nfits = results[3]
        self.extra_output = results[4]

        NoNewAttributesAfterInit.__init__(self)

    def __setstate__(self, state):
        self.__dict__.update(state)

        if 'iterfitname' not in state:
            self.__dict__['iterfitname'] = 'none'

    def __repr__(self):
        return '<%s results instance>' % self.methodname

    def __str__(self):
        return print_fields(self._fields, vars(self))

    def format(self):
        s = ""
        if (self.datasets is not None):
            if len(self.datasets) == 1:
                s = 'Dataset               = %s\n' % str(self.datasets[0])
            else:
                s = 'Datasets              = %s\n' % str(
                    self.datasets).strip("()")
        s += 'Confidence Method     = %s\n' % self.methodname
        if (self.iterfitname is not None or self.iterfitname != 'none'):
            s += 'Iterative Fit Method  = %s\n' % self.iterfitname.capitalize()
        s += 'Fitting Method        = %s\n' % self.fitname
        s += 'Statistic             = %s\n' % self.statname

        s += "%s %g-sigma (%2g%%) bounds:" % (self.methodname, self.sigma,
                                              self.percent)

        def myformat(hfmt, str, lowstr, lownum, highstr, highnum):
            str += hfmt % ('Param', 'Best-Fit', 'Lower Bound', 'Upper Bound')
            str += hfmt % ('-' * 5, '-' * 8, '-' * 11, '-' * 11)

            for name, val, lower, upper in izip(self.parnames, self.parvals,
                                                self.parmins, self.parmaxes):

                str += '\n   %-12s %12g ' % (name, val)
                if is_iterable(lower):
                    str += ' '
                    str += list_to_open_interval(lower)
                elif (lower is None):
                    str += lowstr % '-----'
                else:
                    str += lownum % lower
                if is_iterable(upper):
                    str += '  '
                    str += list_to_open_interval(upper)
                elif (upper is None):
                    str += highstr % '-----'
                else:
                    str += highnum % upper

            return str

        low = map(is_iterable, self.parmins)
        high = map(is_iterable, self.parmaxes)
        in_low = is_in(True, low)
        in_high = is_in(True, high)
        mymethod = self.methodname == 'confidence'

        lowstr = '%12s '
        lownum = '%12g '
        highstr = '%12s'
        highnum = '%12g'

        if in_low and in_high and mymethod:
            hfmt = '\n   %-12s %12s %29s %29s'
            lowstr = '%29s '
            lownum = '%29g '
            highstr = '%30s'
            highnum = '%30g'
        elif in_low and not in_high and mymethod:
            hfmt = '\n   %-12s %12s %29s %12s'
            lowstr = '%29s '
            lownum = '%29g '
            highstr = '%13s'
            highnum = '%13g'
        elif not in_low and in_high and mymethod:
            hfmt = '\n   %-12s %12s %12s %29s'
            highstr = '%29s'
            highnum = '%29g'
        else:
            hfmt = '\n   %-12s %12s %12s %12s'

        return myformat(hfmt, s, lowstr, lownum, highstr, highnum)


class IterFit(NoNewAttributesAfterInit):

    def __init__(self, data, model, stat, method, itermethod_opts=None):
        if itermethod_opts is None:
            itermethod_opts = {'name': 'none'}
        # Even if there is only a single data set, I will
        # want to treat the data and models I am given as
        # collections of data and models -- so, put data and
        # models into the objects needed for simultaneous fitting,
        # if they are not already in such objects.
        self.data = data
        if (type(data) is not DataSimulFit):
            self.data = DataSimulFit('simulfit data', (data,))
        self.model = model
        if (type(model) is not SimulFitModel):
            self.model = SimulFitModel('simulfit model', (model,))
        self.stat = stat
        self.method = method
        # Data set attributes needed to store fitting values between
        # calls to fit
        self._dep = None
        self.extra_args = None
        self._staterror = None
        self._syserror = None
        self._nfev = 0
        self._file = None
        # Options to send to iterative fitting method
        self.itermethod_opts = itermethod_opts
        self.iterate = False
        self.funcs = {'primini': self.primini, 'sigmarej': self.sigmarej}
        self.current_func = None
        if (itermethod_opts['name'] != 'none'):
            self.current_func = self.funcs[itermethod_opts['name']]
            self.iterate = True

    # SIGINT (i.e., typing ctrl-C) can dump the user to the Unix prompt,
    # when signal is sent from G95 compiled code.  What we want is to
    # get to the Sherpa prompt instead.  Typically the user only thinks
    # to interrupt during long fits or projection, so look for SIGINT
    # here, and if it happens, raise the KeyboardInterrupt exception
    # instead of aborting.
    def _sig_handler(self, signum, frame):
        raise KeyboardInterrupt()

    def get_extra_args(self, dep):
        """get the bkg data for wstat and the user defined statistics"""

        def vectorize_backscale_ratio(bkg_backscal, src_backscal, num):
            '''return the ratio of bkg_backscal/src_backscal as a numpy
            array of length n. bkg_backscal and/or src_backscal can be
            either a floating point number or a numpy array'''
            if hasattr(bkg_backscal, '__iter__') is False and \
                    hasattr(src_backscal, '__iter__') is False:
                # both backscal are floating point numbers, generate arrays
                bkg = bkg_backscal * ones(num)
                src = src_backscal
            else:
                # at least one of the backscals is a numpy array
                bkg = bkg_backscal
                src = src_backscal
            return bkg / src

        result = {'bkg': None, 'backscale_ratio': None, 'data_size': None,
                  'exposure_time': None}
        bkg_dep = []
        data_size = None
        exposure_time = None
        backscale_ratio = []

        len_datasets = len(self.data.datasets)
        data_size = zeros(len_datasets, dtype=int)
        exposure_time = zeros(2 * len_datasets)
        for index in xrange(len_datasets):
            mydata = self.data.datasets[index]
            data_size[index] = mydata.get_dep(True).size
        result['data_size'] = data_size
        for index in xrange(len_datasets):
            mydata = self.data.datasets[index]
            if hasattr(mydata, 'response_ids') and \
                    hasattr(mydata, 'background_ids') and \
                    len(mydata.response_ids) and len(mydata.background_ids):
                bkg = mydata.get_background(mydata.background_ids[0])
                tmp_bkg_dep = bkg.get_dep(True)
                # data_size[index] = tmp_bkg_dep.size
                bkg_dep = append(bkg_dep, tmp_bkg_dep)
                exposure_time[2 * index] = mydata.exposure
                exposure_time[2 * index + 1] = bkg.exposure
                ratio = vectorize_backscale_ratio(bkg.backscal,
                                                  mydata.backscal,
                                                  data_size[index])
                backscale_ratio = append(backscale_ratio, ratio)
            else:
                return result

        result['bkg'] = bkg_dep
        result['backscale_ratio'] = backscale_ratio
        # result['data_size'] = data_size
        result['exposure_time'] = exposure_time
        return result

    def _get_callback(self, outfile=None, clobber=False):
        if len(self.model.thawedpars) == 0:
            # raise FitError('model has no thawed parameters')
            raise FitErr('nothawedpar')

        # support Sherpa use with SAMP
        try:
            signal.signal(signal.SIGINT, self._sig_handler)
        except ValueError as e:
            warning(e)

        self._dep, self._staterror, self._syserror = self.data.to_fit(
            self.stat.calc_staterror)

        self.extra_args = self.get_extra_args(self._dep)
        self._nfev = 0
        if outfile is not None:
            if os.path.isfile(outfile) and not clobber:
                # raise FitError("'%s' exists, and clobber==False" % outfile)
                raise FitErr('noclobererr', outfile)
            self._file = open(outfile, 'w')
            names = ['# nfev statistic']
            names.extend(['%s' % par.fullname for par in self.model.pars
                          if not par.frozen])
            print(' '.join(names), file=self._file)

        def cb(pars):
            # We need to store the new parameter values in order to support
            # linked parameters

            self.model.thawedpars = pars
            model = self.data.eval_model_to_fit(self.model)
            stat = self.stat.calc_stat(self._dep, model, self._staterror,
                                       syserror=self._syserror,
                                       extra_args=self.extra_args)
            if self._file is not None:
                vals = ['%5e %5e' % (self._nfev, stat[0])]
                vals.extend(['%5e' % val for val in self.model.thawedpars])
                print(' '.join(vals), file=self._file)

            self._nfev += 1
            return stat

        return cb

    def primini(self, statfunc, pars, parmins, parmaxes, statargs=(),
                statkwargs=None):
        """An iterative scheme, where the variance is computed from
        the model amplitudes.

        This is a chi-square statistic where the variance is computed
        from model amplitudes derived in the previous iteration of the
        fit. This 'Iterative Weighting' ([1]_) attempts to remove
        biased estimates of model parameters which is inherent in
        chi-square2 statistics ([2]_).

        The variance in bin i is estimated to be:

        sigma^2_i^j = S(i, t_s^(j-1)) + (A_s/A_b)^2 B_off(i, t_b^(j-1))

        where j is the number of iterations that have been carried out
        in the fitting process, B_off is the background model
        amplitude in bin i of the off-source region, and t_s^(j-1) and
        t_b^(j-1) are the set of source and background model parameter
        values derived during the iteration previous to the current
        one. The variances are set to an array of ones on the first
        iteration.

        In addition to reducing parameter estimate bias, this
        statistic can be used even when the number of counts in each
        bin is small (< 5), although the user should proceed with
        caution.

        References
        ----------

        .. [1] "Multiparameter linear least-squares fitting to Poisson
               data one count at a time", Wheaton et al. 1995, ApJ 438,
               322
               http://adsabs.harvard.edu/abs/1995ApJ...438..322W

        .. [2] "Bias-Free Parameter Estimation with Few Counts, by
               Iterative Chi-Squared Minimization", Kearns, Primini, &
               Alexander, 1995, ADASS IV, 331
               http://adsabs.harvard.edu/abs/1995ASPC...77..331K

        """
        if statkwargs is None:
            statkwargs = {}
        # Primini's method can only be used with chi-squared;
        # raise exception if it is attempted with least-squares,
        # or maximum likelihood
        if (isinstance(self.stat, Chi2) and
                type(self.stat) is not LeastSq):
            pass
        else:
            raise FitErr('needchi2', 'Primini\'s')
        # Get tolerance, max number of iterations from the
        # dictionary for Primini's method
        tol = self.itermethod_opts['tol']
        if (type(tol) != int and
                type(tol) != float):
            raise SherpaErr(
                "'tol' value for Primini's method must be a number")
        maxiters = self.itermethod_opts['maxiters']
        if (type(maxiters) != int):
            raise SherpaErr(
                "'maxiters' value for Primini's method must be an integer")

        # Store original statistical errors for all data sets.
        # Then, set all statistical errors equal to one, to
        # prepare for Primini's method.
        staterror_original = []

        for d in self.data.datasets:
            st = d.get_staterror(filter=False)
            staterror_original.append(st)
            d.staterror = ones_like(st)

        # Keep record of current and previous statistics;
        # when these are within some tolerace, Primini's method
        # is done.
        previous_stat = float32(finfo(float32).max)
        current_stat = statfunc(pars)[0]
        nfev = 0
        iters = 0

        # This is Primini's method.  The essence of the method
        # is to fit; then call the statistic's staterror function
        # on model values (calculated with current parameter values);
        # then call fit *again*.  Do this until the final fit
        # statistic for the previous and current calls to fit
        # agree to within tolerance.
        final_fit_results = None
        try:
            while (sao_fcmp(previous_stat, current_stat, tol) != 0 and
                   iters < maxiters):
                final_fit_results = self.method.fit(statfunc,
                                                    self.model.thawedpars,
                                                    parmins, parmaxes,
                                                    statargs, statkwargs)
                previous_stat = current_stat
                current_stat = final_fit_results[2]
                nfev += final_fit_results[4].get('nfev')
                iters += 1

                # Call stat.staterror with *model*values*, not data
                # Model values calculated using best-fit parameter
                # values from the just-completed call to the fit
                # function.
                model_iterator = iter(self.model())
                for d in self.data.datasets:
                    d.staterror = self.stat.calc_staterror(
                        d.eval_model(next(model_iterator)))

            # Final number of function evaluations is the sum
            # of the numbers of function evaluations from all calls
            # to the fit function.
            final_fit_results[4]['nfev'] = nfev
        finally:
            # Clean up *no*matter*what* -- we must always
            # restore original statistical errors.
            staterror_original.reverse()
            for d in self.data.datasets:
                d.staterror = staterror_original.pop()

        # Return results from Primini's iterative fitting method
        return final_fit_results

    def sigmarej(self, statfunc, pars, parmins, parmaxes, statargs=(),
                 statkwargs=None):
        """Exclude points that are significately far away from the best fit.

        The `sigmarej` scheme is based on the IRAF `sfit` function
        [1]_, where after a fit data points are excluded if the value
        of `(data-model)/error)` exceeds a threshold, and the data
        re-fit. This removal of data points continues until the fit
        has converged. The error removal can be asymmetric, since
        there are separate parameters for the lower and upper limits.

        References
        ----------

        .. [1] http://iraf.net/irafhelp.php?val=sfit

        """
        if statkwargs is None:
            statkwargs = {}
        # Sigma-rejection can only be used with chi-squared;
        # raise exception if it is attempted with least-squares,
        # or maximum likelihood
        if (isinstance(self.stat, Chi2) and
                type(self.stat) is not LeastSq):
            pass
        else:
            raise FitErr('needchi2', 'Sigma-rejection')

        # Get maximum number of allowed iterations, high and low
        # sigma thresholds for jrection of data points, and
        # "grow" factor (i.e., how many surrounding data points
        # to include with rejected data point).

        maxiters = self.itermethod_opts['maxiters']
        if (type(maxiters) != int):
            raise SherpaErr(
                "'maxiters' value for sigma rejection method must be an integer")
        if (maxiters < 1):
            raise SherpaErr("'maxiters' must be one or greater")

        hrej = self.itermethod_opts['hrej']
        if (type(hrej) != int and
                type(hrej) != float):
            raise SherpaErr(
                "'hrej' value for sigma rejection method must be a number")
        if (not (hrej > 0)):
            raise SherpaErr("'hrej' must be greater than zero")

        lrej = self.itermethod_opts['lrej']
        # FIXME: [OL] There are more reliable ways of checking if an object
        # is (not) a number.
        if (type(lrej) != int and type(lrej) != float):
            raise SherpaErr(
                "'lrej' value for sigma rejection method must be a number")
        if lrej <= 0:
            raise SherpaErr("'lrej' must be greater than zero")

        grow = self.itermethod_opts['grow']
        if (type(grow) != int):
            raise SherpaErr(
                "'grow' value for sigma rejection method must be an integer")
        if (grow < 0):
            raise SherpaErr("'grow' factor must be zero or greater")

        # Keep record of current and previous statistics
        previous_stat = float32(finfo(float32).max)
        current_stat = statfunc(pars)[0]
        nfev = 0
        iters = 0

        # Store original masks (filters) for each data set.
        mask_original = []

        for d in self.data.datasets:
            # If there's no filter, create a filter that is
            # all True
            if (iterable(d.mask) != True):
                mask_original.append(d.mask)
                d.mask = ones_like(array(d.get_dep(False), dtype=bool))
            else:
                mask_original.append(array(d.mask))

        self.model.teardown()
        final_fit_results = None
        rejected = True
        try:
            while (rejected and iters < maxiters):
                # Update stored y, staterror and syserror values
                # from data, so callback function will work properly
                self._dep, self._staterror, self._syserror = self.data.to_fit(
                    self.stat.calc_staterror)
                self.model.startup()
                final_fit_results = self.method.fit(statfunc,
                                                    self.model.thawedpars,
                                                    parmins, parmaxes,
                                                    statargs, statkwargs)
                model_iterator = iter(self.model())
                rejected = False

                for d in self.data.datasets:
                    # For each data set, compute
                    # (data - model) / staterror
                    # over filtered data space
                    residuals = (d.get_dep(True) - d.eval_model_to_fit(
                        next(model_iterator))) / d.get_staterror(True, self.stat.calc_staterror)

                    # For each modeled value that exceeds
                    # sigma thresholds, set the corresponding
                    # filter value from True to False
                    ressize = len(residuals)
                    filsize = len(d.mask)
                    newmask = d.mask

                    j = 0
                    kmin = 0
                    for i in xrange(0, ressize):
                        while (newmask[j] == False and
                               j < filsize):
                            j = j + 1
                        if (j >= filsize):
                            break
                        if (residuals[i] <= -lrej or
                                residuals[i] >= hrej):
                            rejected = True
                            kmin = j - grow
                            if (kmin < 0):
                                kmin = 0
                            kmax = j + grow
                            if (kmax >= filsize):
                                kmax = filsize - 1
                            for k in xrange(kmin, kmax + 1):
                                newmask[k] = False
                        j = j + 1

                        # If we've masked out *all* data,
                        # immediately raise fit error, clean up
                        # on way out.
                        if (any(newmask) == False):
                            raise FitErr('nobins')
                        d.mask = newmask

                # For data sets with backgrounds, correct that
                # backgrounds have masks that match their sources
                for d in self.data.datasets:
                    if (hasattr(d, "background_ids") == True and
                            hasattr(d, "get_background") == True):
                        for bid in d.background_ids:
                            b = d.get_background(bid)
                            if (iterable(b.mask) == True and
                                    iterable(d.mask) == True):
                                if (len(b.mask) == len(d.mask)):
                                    b.mask = d.mask

                # teardown model, get ready for next iteration
                self.model.teardown()
                iters = iters + 1
                nfev += final_fit_results[4].get('nfev')
                final_fit_results[4]['nfev'] = nfev
        except:
            # Clean up if exception occurred
            mask_original.reverse()
            for d in self.data.datasets:
                d.mask = mask_original.pop()

            # Update stored y, staterror and syserror values
            # from data, so callback function will work properly
            self._dep, self._staterror, self._syserror = self.data.to_fit(
                self.stat.calc_staterror)
            self.model.startup()
            raise

        self._dep, self._staterror, self._syserror = self.data.to_fit(
            self.stat.calc_staterror)
        self.model.startup()

        # N.B. -- If sigma-rejection in Sherpa 3.4 succeeded,
        # it did *not* restore the filter to its state before
        # sigma-rejection was called.  If points are filtered
        # out, they stay out.  So we emulate the same behavior
        # if our version of sigma-rejection succeeds.

        # Mind you, if sigma-rejection *fails*, then we *do*
        # restore the filter, and re-raise the exception in
        # the above exception block.

        # Return results from sigma rejection
        return final_fit_results

    def fit(self, statfunc, pars, parmins, parmaxes, statargs=(), statkwargs=None):
        if statkwargs is None:
            statkwargs = {}
        if not self.iterate:
            return self.method.fit(statfunc, pars, parmins, parmaxes,
                                   statargs, statkwargs)
        else:
            return self.current_func(statfunc, pars, parmins, parmaxes,
                                     statargs, statkwargs)


class Fit(NoNewAttributesAfterInit):

    def __init__(self, data, model, stat=None, method=None, estmethod=None,
                 itermethod_opts=None):
        if itermethod_opts is None:
            itermethod_opts = {'name': 'none'}
        self.data = data
        self.model = model

        if stat is None:
            stat = Chi2Gehrels()
        if method is None:
            method = LevMar()
        if estmethod is None:
            estmethod = Covariance()

        self.stat = stat
        self.method = method
        self.estmethod = estmethod
        # Confidence limit code freezes one parameter
        # at a time.  Keep a record here of which one
        # that is, in case an exception is raised and
        # this parameter needs to be thawed in the
        # exception handler.
        self.thaw_indices = ()
        iter = 0
        for current_par in self.model.pars:
            if current_par.frozen is True:
                pass
            else:
                self.thaw_indices = self.thaw_indices + (iter,)
            iter = iter + 1
        self.current_frozen = -1

        # The number of times that reminimization has occurred
        # during an attempt to compute confidence limits.  If
        # that number equals self.estmethod.maxfits, cease all
        # further attempt to reminimize.
        self.refits = 0

        # Set up an IterFit object, so that the user can select
        # an iterative fitting option.
        self._iterfit = IterFit(self.data, self.model, self.stat, self.method,
                                itermethod_opts)
        NoNewAttributesAfterInit.__init__(self)

    def __setstate__(self, state):
        self.__dict__.update(state)

        if '_iterfit' not in state:
            self.__dict__['_iterfit'] = IterFit(self.data, self.model, self.stat, self.method,
                                                {'name': 'none'})

    def __str__(self):
        return (('data      = %s\n' +
                 'model     = %s\n' +
                 'stat      = %s\n' +
                 'method    = %s\n' +
                 'estmethod = %s') %
                (self.data.name,
                 self.model.name,
                 type(self.stat).__name__,
                 type(self.method).__name__,
                 type(self.estmethod).__name__))

    def guess(self, **kwargs):
        """
        kwargs = { 'limits' : True, 'values' : False }
        """
        self.model.guess(*self.data.to_guess(), **kwargs)

    # QUS: should this have an @evaluates_model decorator?
    #
    def _calc_stat(self):
        """Calculate the current statistic value.

        Returns
        -------
        statval, fvec : number, array of numbers
            The overall statistic value and the "per-bin" value.
        """

        # TODO: is there anything missing here that
        #       self._iterfit.get_extra_args calculates?
        return self.stat.calc_stat_from_data(self.data, self.model)

    def calc_stat(self):
        """Calculate the statistic value.

        Evaluate the statistic for the current model and data
        settings (e.g. parameter values and data filters).

        Returns
        -------
        stat : number
           The current statistic value.

        See Also
        --------
        calc_chisqr, calc_stat_info

        """

        return self._calc_stat()[0]

    def calc_chisqr(self):
        """Calculate the per-bin chi-squared statistic.

        Evaluate the per-bin statistic for the current model and data
        settings (e.g. parameter values and data filters).

        Returns
        -------
        chisq : array or `None`
           The chi-square value for each bin of the data, using the
           current statistic (as set by `set_stat`).  A value of
           `None` is returned if the statistic is not a chi-square
           distribution.

        See Also
        --------
        calc_stat, calc_stat_info

        """

        # Since there is some setup work needed before calling
        # this routine, and to avoid catching any AttributeErrors
        # thrown by the routine, use this un-pythonic check.
        #
        if not hasattr(self.stat, 'calc_chisqr'):
            return None

        return self.stat.calc_chisqr(self.data, self.model)

    def calc_stat_info(self):
        """Calculate the statistic value and related information.

        Evaluate the statistic for the current model and data
        settings (e.g. parameter values and data filters).

        Returns
        -------
        statinfo : StatInfoResults instance
           The current statistic value.

        See Also
        --------
        calc_chisqr, calc_stat

        """

        # TODO: does this need to be updated due to changes in _calc_stat?
        #       also, this logic would be better in the stat class
        #       than here
        #
        stat, fvec = self._calc_stat()
        model = self.data.eval_model_to_fit(self.model)

        qval = None
        rstat = None
        numpoints = len(model)
        dof = numpoints - len(self.model.thawedpars)
        if (isinstance(self.stat, (CStat, Chi2)) and
                not isinstance(self.stat, LeastSq)):
            if stat >= 0.0:
                qval = igamc(dof / 2., stat / 2.)
            rstat = stat / dof

        name = self.stat.name

        if isinstance(self.stat, Chi2) and not isinstance(self.stat, LeastSq):
            isSimulFit = isinstance(self.data, DataSimulFit)
            if isSimulFit:
                is_error_set = [
                    d.staterror is not None for d in self.data.datasets]
                if all(is_error_set):
                    name = 'chi2'
            elif self.data.staterror is not None:
                name = 'chi2'

        return StatInfoResults(name, stat, numpoints, model,
                               dof, qval, rstat)

    @evaluates_model
    def fit(self, outfile=None, clobber=False):
        dep, staterror, syserror = self.data.to_fit(self.stat.calc_staterror)

        # TODO: add tests to check this happens before it can be removed,
        #       as it should be handled by the calc_stat routine
        if not iterable(dep) or len(dep) == 0:
            # raise FitError('no noticed bins found in data set')
            raise FitErr('nobins')

        if ((iterable(staterror) and 0.0 in staterror) and
                isinstance(self.stat, Chi2) and
                type(self.stat) != Chi2 and
                type(self.stat) != Chi2ModVar):
            # raise FitError('zeros found in uncertainties, consider using' +
            #               ' calculated uncertainties')
            raise FitErr('binhas0')

        init_stat = self.calc_stat()
        # output = self.method.fit ...
        output = self._iterfit.fit(self._iterfit._get_callback(outfile, clobber),
                                   self.model.thawedpars,
                                   self.model.thawedparmins,
                                   self.model.thawedparmaxes)
        # LevMar always calculate chisquare, so call calc_stat
        # just in case statistics is something other then chisquare
        self.model.thawedpars = output[1]
        tmp = list(output)
        tmp[2] = self.calc_stat()
        output = tuple(tmp)
        # end of the gymnastics 'cause one cannot write to a tuple

        # check if any parameter values are at boundaries,
        # and warn user.
        tol = finfo(float32).eps
        param_warnings = ""
        for par in self.model.pars:
            if not par.frozen:
                if sao_fcmp(par.val, par.min, tol) == 0:
                    param_warnings += ("WARNING: parameter value %s is at its minimum boundary %s\n" %
                                       (par.fullname, str(par.min)))
                if sao_fcmp(par.val, par.max, tol) == 0:
                    param_warnings += ("WARNING: parameter value %s is at its maximum boundary %s\n" %
                                       (par.fullname, str(par.max)))

        if self._iterfit._file is not None:
            vals = ['%5e %5e' % (self._iterfit._nfev, tmp[2])]
            vals.extend(['%5e' % val for val in self.model.thawedpars])
            print(' '.join(vals), file=self._iterfit._file)
            self._iterfit._file.close()
            self._iterfit._file = None

        return FitResults(self, output, init_stat, param_warnings.strip("\n"))

    @evaluates_model
    def simulfit(self, *others):
        if len(others) == 0:
            return self.fit()

        fits = (self,) + others
        d = DataSimulFit('simulfit data', tuple(f.data for f in fits))
        m = SimulFitModel('simulfit model', tuple(f.model for f in fits))

        f = Fit(d, m, self.stat, self.method)
        return f.fit()

    @evaluates_model
    def est_errors(self, methoddict=None, parlist=None):
        # Define functions to freeze and thaw a parameter before
        # we call fit function -- projection can call fit several
        # times, for each parameter -- that parameter must be frozen
        # while the others freely vary.
        def freeze_par(pars, parmins, parmaxes, i):
            # Freeze the indicated parameter; return
            # its place in the list of all parameters,
            # and the current values of the parameters,
            # and the hard mins amd maxs of the parameters
            self.model.pars[self.thaw_indices[i]].val = pars[i]
            self.model.pars[self.thaw_indices[i]].frozen = True
            self.current_frozen = self.thaw_indices[i]
            keep_pars = ones_like(pars)
            keep_pars[i] = 0
            current_pars = pars[where(keep_pars)]
            current_parmins = parmins[where(keep_pars)]
            current_parmaxes = parmaxes[where(keep_pars)]
            return (current_pars, current_parmins, current_parmaxes)

        def thaw_par(i):
            if (i < 0):
                pass
            else:
                self.model.pars[self.thaw_indices[i]].frozen = False
                self.current_frozen = -1

        # confidence needs to know which parameter it is working on.
        def get_par_name(ii):
            return self.model.pars[self.thaw_indices[ii]].fullname

        # Call from a parameter estimation method, to report
        # that limits for a given parameter have been found
        def report_progress(i, lower, upper):
            if (i < 0):
                pass
            else:
                name = self.model.pars[self.thaw_indices[i]].fullname
                if isnan(lower) or isinf(lower):
                    info("%s \tlower bound: -----" % name)
                else:
                    info("%s \tlower bound: %g" % (name, lower))
                if isnan(upper) or isinf(upper):
                    info("%s \tupper bound: -----" % name)
                else:
                    info("%s \tupper bound: %g" % (name, upper))

        # If starting fit statistic is chi-squared or C-stat,
        # can calculate reduced fit statistic -- if it is
        # more than 3, don't bother calling method to estimate
        # parameter limits.

        if (type(self.stat) is LeastSq):
            # raise FitError('cannot estimate confidence limits with ' +
            #               type(self.stat).__name__)
            raise EstErr('noerr4least2', type(self.stat).__name__)

        if (type(self.stat) is not Cash):
            dep, staterror, syserror = self.data.to_fit(
                self.stat.calc_staterror)

            if not iterable(dep) or len(dep) == 0:
                #raise FitError('no noticed bins found in data set')
                raise FitErr('nobins')

            # For chi-squared and C-stat, reduced statistic is
            # statistic value divided by number of degrees of
            # freedom.

            # Degress of freedom are number of data bins included
            # in fit, minus the number of thawed parameters.
            dof = len(dep) - len(self.model.thawedpars)
            if (dof < 1):
                #raise FitError('degrees of freedom are zero or lower')
                raise EstErr('nodegfreedom')

            if (hasattr(self.estmethod, "max_rstat") and
                    (self.calc_stat() / dof) > self.estmethod.max_rstat):
                # raise FitError('reduced statistic larger than ' +
                #               str(self.estmethod.max_rstat))
                raise EstErr('rstat>max', str(self.estmethod.max_rstat))

        # If statistic is chi-squared, change fitting method to
        # Levenberg-Marquardt; else, switch to NelderMead.  (We
        # will do fitting during projection, and therefore don't
        # want to use LM with a stat other than chi-squared).

        # If current method is not LM or NM, warn it is not a good
        # method for estimating parameter limits.
        if (type(self.estmethod) is not Covariance and
                type(self.method) is not NelderMead and
                type(self.method) is not LevMar):
            warning(self.method.name + " is inappropriate for confidence " +
                    "limit estimation")

        oldmethod = self.method
        if (hasattr(self.estmethod, "fast") and
                bool_cast(self.estmethod.fast) and
                methoddict is not None):
            if (isinstance(self.stat, Likelihood)):
                if (type(self.method) is not NelderMead):
                    self.method = methoddict['neldermead']
                    warning("Setting optimization to " + self.method.name
                            + " for confidence limit search")
            else:
                if (type(self.method) is not LevMar):
                    self.method = methoddict['levmar']
                    warning("Setting optimization to " + self.method.name
                            + " for confidence limit search")

        # Now, set up before we call the confidence limit function
        # Keep track of starting values, will need to set parameters
        # back to starting values when we are done.
        startpars = self.model.thawedpars
        startsoftmins = self.model.thawedparmins
        startsoftmaxs = self.model.thawedparmaxes
        starthardmins = self.model.thawedparhardmins
        starthardmaxs = self.model.thawedparhardmaxes

        # If restricted to soft_limits, only send soft limits to
        # method, and do not reset model limits
        if (bool_cast(self.estmethod.soft_limits) is True):
            starthardmins = self.model.thawedparmins
            starthardmaxs = self.model.thawedparmaxes
        else:
            self.model.thawedparmins = starthardmins
            self.model.thawedparmaxes = starthardmaxs

        self.current_frozen = -1

        # parnums is the list of indices of the thawed parameters
        # we want to visit.  For example, if there are three thawed
        # parameters, and we want to derive limits for only the first
        # and third, then parnums = [0,2].  We construct the list by
        # comparing each parameter in parlist to the thawed model
        # parameters.  (In the default case, when parlist is None,
        # that means get limits for all thawed parameters, so parnums
        # is [0, ... , numpars - 1], if the number of thawed parameters
        # is numpars.)
        parnums = []
        if parlist is not None:
            allpars = [p for p in self.model.pars if not p.frozen]
            for p in parlist:
                count = 0
                match = False
                for par in allpars:
                    if p is par:
                        parnums.append(count)
                        match = True
                    count = count + 1
                if not match:
                    raise EstErr('noparameter', p.fullname)
            parnums = array(parnums)
        else:
            parlist = [p for p in self.model.pars if not p.frozen]
            parnums = arange(len(startpars))

        # If we are here, we are ready to try to derive confidence limits.
        # General rule:  if failure because a hard limit was hit, find
        # out which parameter it was so we can tell the user.
        # If a new minimum statistic was found, start over, with parameter
        # values that yielded new lower statistic as the new starting point.
        output = None
        results = None
        oldremin = -1.0
        if (hasattr(self.estmethod, "remin")):
            oldremin = self.estmethod.remin
        try:
            output = self.estmethod.compute(self._iterfit._get_callback(),
                                            self._iterfit.fit,
                                            self.model.thawedpars,
                                            startsoftmins,
                                            startsoftmaxs,
                                            starthardmins,
                                            starthardmaxs,
                                            parnums,
                                            freeze_par, thaw_par,
                                            report_progress, get_par_name)
        except EstNewMin as e:
            # If maximum number of refits has occurred, don't
            # try to reminimize again.
            if (hasattr(self.estmethod, "maxfits") and
                    not(self.refits < self.estmethod.maxfits - 1)):
                self.refits = 0
                thaw_par(self.current_frozen)
                self.model.thawedpars = startpars
                self.model.thawedparmins = startsoftmins
                self.model.thawedparmaxes = startsoftmaxs
                self.method = oldmethod
                if (hasattr(self.estmethod, "remin")):
                    self.estmethod.remin = -1.0
                warning("Maximum number of reminimizations reached")

            # First report results of new fit, then call
            # compute limits for those new best-fit parameters
            for p in parlist:
                p.frozen = False
            self.current_frozen = -1

            if e.args != ():
                self.model.thawedpars = e.args[0]

            self.model.thawedparmins = startsoftmins
            self.model.thawedparmaxes = startsoftmaxs
            results = self.fit()
            self.refits = self.refits + 1
            warning(
                "New minimum statistic found while computing confidence limits")
            warning("New best-fit parameters:\n" + results.format())

            # Now, recompute errors for new best-fit parameters
            results = self.est_errors(methoddict, parlist)
            self.model.thawedparmins = startsoftmins
            self.model.thawedparmaxes = startsoftmaxs
            self.method = oldmethod
            if (hasattr(self.estmethod, "remin")):
                self.estmethod.remin = oldremin
            return results
        except:
            for p in parlist:
                p.frozen = False
            self.current_frozen = -1
            self.model.thawedpars = startpars
            self.model.thawedparmins = startsoftmins
            self.model.thawedparmaxes = startsoftmaxs
            self.method = oldmethod
            if (hasattr(self.estmethod, "remin")):
                self.estmethod.remin = oldremin
            raise

        for p in parlist:
            p.frozen = False
        self.current_frozen = -1
        self.model.thawedpars = startpars
        self.model.thawedparmins = startsoftmins
        self.model.thawedparmaxes = startsoftmaxs
        results = ErrorEstResults(self, output, parlist)
        self.method = oldmethod
        if (hasattr(self.estmethod, "remin")):
            self.estmethod.remin = oldremin

        return results
