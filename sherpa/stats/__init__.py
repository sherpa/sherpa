# 
#  Copyright (C) 2009  Smithsonian Astrophysical Observatory
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

import numpy
from sherpa.utils import NoNewAttributesAfterInit
from sherpa.utils.err import StatErr
import sherpa.stats._statfcts


__all__ = ('Stat', 'Cash', 'CStat', 'LeastSq',
           'Chi2Gehrels', 'Chi2ConstVar', 'Chi2DataVar', 'Chi2ModVar',
           'Chi2XspecVar', 'Chi2')


from sherpa import get_config
from ConfigParser import ConfigParser

config = ConfigParser()
config.read(get_config())

# truncation_flag indicates whether or not model truncation
# should be performed.  If true, use the truncation_value from
# the config file.
truncation_flag = config.get('statistics','truncate').upper()
truncation_value = float(config.get('statistics','trunc_value'))
if (bool(truncation_flag) is False or truncation_flag == "FALSE" or
    truncation_flag == "NONE" or truncation_flag == "0"):
    truncation_value = -1.0


class Stat(NoNewAttributesAfterInit):

    def __init__(self, name):
        self.name = name
        NoNewAttributesAfterInit.__init__(self)

    def __repr__(self):
        if self.__doc__ is not None:
            return self.__doc__
        return ("<%s statistic instance '%s'>" %
                (type(self).__name__, self.name))

    def calc_staterror(self, data):
        raise NotImplementedError

    def calc_stat(self, data, model, staterror=None, syserror=None,
                  weight=None):
        raise NotImplementedError

class Likelihood(Stat):
    """Maximum likelihood function"""
    def __init__(self, name='likelihood'):
        Stat.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        # Likelihood stats do not have 'errors' associated with them.
        # return 1 to avoid dividing by 0 by some optimization methods.
        return numpy.ones_like(data)

### Ahelp ingest: 2015-05-04 DJB
## DOC-TODO: where is the truncate/trunc_value stored for objects
##           AHA: it appears to be taken straight from the config
##           file rather than associated with the Stat class
## DOC-TODO: where to talk about the .sherpa.rc config file?
class Cash(Likelihood):
    """Maximum likelihood function.

    Counts are sampled from the Poisson distribution, and so the best
    way to assess the quality of model fits is to use the product of
    individual Poisson probabilities computed in each bin i, or the
    likelihood L:

    L = (product)_i [ M(i)^D(i)/D(i)! ] * exp[-M(i)]

    where M(i) = S(i) + B(i) is the sum of source and background model
    amplitudes, and D(i) is the number of observed counts, in bin i.

    The Cash statistic [1]_ is derived by (1) taking the logarithm of
    the likelihood function, (2) changing its sign, (3) dropping the
    factorial term (which remains constant during fits to the same
    dataset), and (4) multiplying by two:

    C = 2 * (sum)_i [ M(i) - D(i) log M(i) ]

    The factor of two exists so that the change in cash statistic from
    one model fit to the next, (Delta)C, is distributed approximately
    as (Delta)chi-square when the number of counts in each bin is
    high.  One can then in principle use (Delta)C instead of
    (Delta)chi-square in certain model comparison tests. However,
    unlike chi-square, the cash statistic may be used regardless of
    the number of counts in each bin.

    The magnitude of the Cash statistic depends upon the number of
    bins included in the fit and the values of the data
    themselves.  Hence one cannot analytically assign a
    goodness-of-fit measure to a given value of the Cash statistic.
    Such a measure can, in principle, be computed by performing
    Monte Carlo simulations. One would repeatedly sample new
    datasets from the best-fit model, fit them, and note where the
    observed Cash statistic lies within the derived distribution
    of Cash statistics. Alternatively, the `cstat` statistic can
    be used.

    Notes
    -----
    The background should not be subtracted from the data when this
    statistic is used.  It should be modeled simultaneously with the
    source.

    The Cash statistic function evaluates the logarithm of each data
    point. If the number of counts is zero or negative, it's not
    possible to take the log of that number. The behavior in this case
    is controlled by the `truncate` and `trunc_value` settings in the
    .sherpa.rc file:

    - if `truncate` is `True` (the default value), then
      `log(trunc_value)` is used whenever the data value is <= 0.  The
      default is `trunc_value=1.0e-25`.

    - when `truncate` is `False` an error is raised.

    References
    ----------

    .. [1] "Parameter estimation in astronomy through application of
           the likelihood ratio", Cash, W. 1979, ApJ 228, 939
           http://adsabs.harvard.edu/abs/1979ApJ...228..939C

    """
    def __init__(self, name='cash'):
        Likelihood.__init__(self, name)

    @staticmethod
    def calc_stat(data, model, staterror=None, syserror=None, weight=None):
        return _statfcts.calc_cash_stat(data, model, staterror, syserror,
                                        weight, truncation_value)


### Ahelp ingest: 2015-05-04 DJB
class CStat(Likelihood):
    """Maximum likelihood function (XSPEC style).

    This is equivalent to the X-Spec implementation of the
    Cash statistic [1]_. It does *not* include the background
    modelling.

    Counts are sampled from the Poisson distribution, and so the best
    way to assess the quality of model fits is to use the product of
    individual Poisson probabilities computed in each bin i, or the
    likelihood L:

    L = (product)_i [ M(i)^D(i)/D(i)! ] * exp[-M(i)]

    where M(i) = S(i) + B(i) is the sum of source and background model
    amplitudes, and D(i) is the number of observed counts, in bin i.

    The cstat statistic is derived by (1) taking the logarithm of the
    likelihood function, (2) changing its sign, (3) dropping the
    factorial term (which remains constant during fits to the same
    dataset), (4) adding an extra data-dependent term (this is what
    makes it different to `Cash`, and (5) multiplying by two:

    C = 2 * (sum)_i [ M(i) - D(i) + D(i)*[log D(i) - log M(i)] ]

    The factor of two exists so that the change in the cstat statistic
    from one model fit to the next, (Delta)C, is distributed
    approximately as (Delta)chi-square when the number of counts in
    each bin is high.  One can then in principle use (Delta)C instead
    of (Delta)chi-square in certain model comparison tests. However,
    unlike chi-square, the cstat statistic may be used regardless of
    the number of counts in each bin.

    The inclusion of the data term in the expression means that,
    unlike the Cash statistic, one can assign an approximate
    goodness-of-fit measure to a given value of the cstat statistic,
    i.e. the observed statistic, divided by the number of degrees of
    freedom, should be of order 1 for good fits.

    Notes
    -----
    The background should not be subtracted from the data when this
    statistic is used.  It should be modeled simultaneously with the
    source.

    The cstat statistic function evaluates the logarithm of each data
    point. If the number of counts is zero or negative, it's not
    possible to take the log of that number. The behavior in this case
    is controlled by the `truncate` and `trunc_value` settings in the
    .sherpa.rc file:

    - if `truncate` is `True` (the default value), then
      `log(trunc_value)` is used whenever the data value is <= 0.  The
      default is `trunc_value=1.0e-25`.

    - when `truncate` is `False` an error is raised.

    References
    ----------

    .. [1] The description of the Cash statistic (`cstat`) in
           https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html

    """
    def __init__(self, name='cstat'):
        Likelihood.__init__(self, name)

    @staticmethod
    def calc_stat(data, model, staterror=None, syserror=None, weight=None):
        return _statfcts.calc_cstat_stat(data, model, staterror, syserror,
                                         weight, truncation_value)


class Chi2(Stat):
    """Chi Squared"""
    def __init__(self, name='chi2'):
        Stat.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
         raise StatErr('chi2noerr')

    @staticmethod
    def calc_stat(data, model, staterror, syserror=None, weight=None):
        return _statfcts.calc_chi2_stat(data, model, staterror,
                                        syserror, weight, truncation_value)

class LeastSq(Chi2):
    """Least Squared"""
    def __init__(self, name='leastsq'):
        Stat.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        return numpy.ones_like(data)        

    @staticmethod
    def calc_stat(data, model, staterror, syserror=None, weight=None):
        return _statfcts.calc_lsq_stat(data, model, staterror,
                                       syserror, weight, truncation_value)
    

class Chi2Gehrels(Chi2):
    """Chi Squared with Gehrels variance"""
    def __init__(self, name='chi2gehrels'):
        Chi2.__init__(self, name)

    calc_staterror = _statfcts.calc_chi2gehrels_errors


class Chi2ConstVar(Chi2):
    """Chi Squared with constant variance"""
    def __init__(self, name='chi2constvar'):
        Chi2.__init__(self, name)

    calc_staterror = _statfcts.calc_chi2constvar_errors


class Chi2DataVar(Chi2):
    """Chi Squared with data variance"""
    def __init__(self, name='chi2datavar'):
        Chi2.__init__(self, name)

    calc_staterror = _statfcts.calc_chi2datavar_errors


class Chi2ModVar(Chi2):
    """Chi Squared with model amplitude variance"""
    def __init__(self, name='chi2modvar'):
        Chi2.__init__(self, name)

    # Statistical errors are not used
    @staticmethod
    def calc_staterror(data):
        return numpy.zeros_like(data)

    @staticmethod
    def calc_stat(data, model, staterror, syserror=None, weight=None):
        return _statfcts.calc_chi2modvar_stat(data, model, staterror,
                                              syserror, weight,
                                              truncation_value)


class Chi2XspecVar(Chi2):
    """Chi Squared with data variance (XSPEC style)"""
    def __init__(self, name='chi2xspecvar'):
        Chi2.__init__(self, name)

    calc_staterror = _statfcts.calc_chi2xspecvar_errors


class UserStat(Stat):

    def __init__(self, statfunc=None, errfunc=None, name='userstat'):
        self._statfuncset = False
        self.statfunc = (lambda x: None)
        
        self._staterrfuncset = False
        self.errfunc = (lambda x: None)
        
        if statfunc is not None:
            self.statfunc = statfunc
            self._statfuncset = True

        if errfunc is not None:
            self.errfunc = errfunc
            self._staterrfuncset = True
            

        Stat.__init__(self, name)


    def __getstate__(self):
        state = self.__dict__.copy()
        # Function pointers to methods of the class
        # (of type 'instancemethod') are NOT picklable
        # remove them and restore later with a coord init
        del state['statfunc']
        del state['errfunc']

        return state

    def __setstate__(self, state):
        # Populate the function pointers we deleted at pickle time with
        # no-ops.
        self.__dict__['statfunc']=(lambda x: None)
        self.__dict__['errfunc']=(lambda x: None)
        self.__dict__.update(state)


    def set_statfunc(self, func):
        self.statfunc = func
        self._statfuncset = True


    def set_errfunc(self, func):
        self.errfunc = func
        self._staterrfuncset = True


    def calc_staterror(self, data):
        if not self._staterrfuncset:
            raise StatErr('nostat', self.name, 'calc_staterror()')
        return self.errfunc(data)


    def calc_stat(self, data, model, staterror=None, syserror=None,
                  weight=None):
        if not self._statfuncset:
            raise StatErr('nostat', self.name, 'calc_stat()')
        return self.statfunc(data, model, staterror, syserror, weight)
