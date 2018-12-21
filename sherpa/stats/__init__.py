#
#  Copyright (C) 2009, 2015, 2016, 2017, 2018
#             Smithsonian Astrophysical Observatory
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

from __future__ import absolute_import

import warnings

import numpy
from sherpa.utils import NoNewAttributesAfterInit, igamc
from sherpa.utils.err import FitErr, StatErr
from sherpa.data import DataSimulFit
from sherpa.models import SimulFitModel

from . import _statfcts

from sherpa import get_config
from six.moves.configparser import ConfigParser


__all__ = ('Stat', 'Cash', 'CStat', 'LeastSq',
           'Chi2Gehrels', 'Chi2ConstVar', 'Chi2DataVar', 'Chi2ModVar',
           'Chi2XspecVar', 'Chi2',
           'UserStat', 'WStat')


config = ConfigParser()
config.read(get_config())

# truncation_flag indicates whether or not model truncation
# should be performed.  If true, use the truncation_value from
# the config file.
truncation_flag = config.get('statistics', 'truncate').upper()
truncation_value = float(config.get('statistics', 'trunc_value'))
if (bool(truncation_flag) is False or truncation_flag == "FALSE" or
        truncation_flag == "NONE" or truncation_flag == "0"):
    truncation_value = 1.0e-25


class Stat(NoNewAttributesAfterInit):
    """The base class for calculating a statistic given data and model."""

    # Used by calc_stat
    #
    _calc = None

    # This should be overridden by derived classes and set to True if the rstat and qvalue
    # figures can be calculated for that class of statistics.
    _can_calculate_rstat = None

    def __init__(self, name):
        self.name = name
        NoNewAttributesAfterInit.__init__(self)

    def __repr__(self):
        if self.__doc__ is not None:
            return self.__doc__
        return ("<%s statistic instance '%s'>" %
                (type(self).__name__, self.name))

    @staticmethod
    def _bundle_inputs(data, model):
        """Convert input into SimulFit instances.

        Convert the inputs into DataSimulFit and SimulFitModel
        instances.

        Parameters
        ----------
        data : a Data or DataSimulFit instance
            The data set, or sets, to use.
        model : a Model or SimulFitModel instance
            The model expression, or expressions. If a SimulFitModel
            is given then it must match the number of data sets in the
            data parameter.

        Returns
        -------
        data, model : DataSimulFit instance, SimulFitModel instance
            If the input was a SimulFit object then this is just
            the input value.
        """

        if not isinstance(data, DataSimulFit):
            data = DataSimulFit('simulfit data', (data,))

        if not isinstance(model, SimulFitModel):
            model = SimulFitModel('simulfit model', (model,))

        return data, model

    @staticmethod
    def _check_has_bins(data):
        """Raise an error if there are no noticed bins in the dataset.

        Parameters
        ----------
        data : a DataSimulFit instance
            The data sets to use.

        Raises
        ------
        FitErr

        Notes
        -----
        It is unclear whether this should error out if any particular
        dataset has no noticed bins, or only if all the datasets have
        no noticed bins (the latter approach is taken).

        """

        for dset in data.datasets:
            # Assume that the error column does not need to be
            # calculated for this check, so staterrfunc can be
            # None.
            dep, _, _ = dset.to_fit(staterrfunc=None)
            if numpy.iterable(dep) and len(dep) > 0:
                return

        raise FitErr('nobins')

    @staticmethod
    def _check_sizes_match(data, model):
        """Raise an error if number of datasets and models do not match.

        Parameters
        ----------
        data : a DataSimulFit instance
            The data sets to use.
        model : a SimulFitModel instance
            The model expressions for each data set. It must match
            the data parameter (the models are in the same order
            as the data objects).

        Raises
        ------
        StatErr

        """

        ndata = len(data.datasets)
        nmdl = len(model.parts)
        if ndata != nmdl:
            raise StatErr('mismatch',
                          'number of data sets', ndata,
                          'model expressions', nmdl)

    def _validate_inputs(self, data, model):
        """Ensure that the inputs are correct for the statistic.

        The default behavior is to check that the data contains
        at least one bin and that the number of datasets matches
        the number of models. It also converts single values to
        simultaneous objects, if necessary.

        Parameters
        ----------
        data : a Data or DataSimulFit instance
            The data set, or sets, to use.
        model : a Model or SimulFitModel instance
            The model expressions for each data set. It must match
            the data parameter (the models are in the same order
            as the data objects).

        Returns
        -------
        data, model : DataSimulFit instance, SimulFitModel instance
            If the input was a SimulFit object then this is just
            the input value.
        """

        if self._calc is None:
            # This is a programmer error rather than a user error,
            # so use NotImplementedError rather than
            # StatErr('nostat', self.name, '_calc')
            #
            raise NotImplementedError("_calc method has not been set")

        data, model = self._bundle_inputs(data, model)
        self._check_has_bins(data)
        self._check_sizes_match(data, model)
        return data, model

    def _get_fit_model_data(self, data, model):
        data, model = self._validate_inputs(data, model)
        fitdata = data.to_fit(staterrfunc=self.calc_staterror)
        modeldata = data.eval_model_to_fit(model)

        return fitdata, modeldata

    # TODO:
    #  - should this accept sherpa.data.Data input instead of
    #    "raw" data (i.e. to match calc_stat)
    #  - should this be moved out of the base Stat class since
    #    isn't relevant for likelihood statistics?
    #
    def calc_staterror(self, data):
        """Return the statistic error values for the data.

        Parameters
        ----------
        data : scalar or 1D array of numbers
            The data values.

        Returns
        -------
        staterror : scalar or array of numbers
            The errors for the input data values (matches the data
            argument).

        """
        raise NotImplementedError

    # TODO: add *args, **kwargs?
    def calc_stat(self, data, model):
        """Return the statistic value for the data and model.

        Parameters
        ----------
        data : a Data or DataSimulFit instance
            The data set, or sets, to use.
        model : a Model or SimulFitModel instance
            The model expression, or expressions. If a SimulFitModel
            is given then it must match the number of data sets in the
            data parameter.

        Returns
        -------
        statval, fvec : number, array of numbers
            The statistic value and the per-bin "statistic" value.

        """

        raise NotImplementedError

    def goodness_of_fit(self, statval, dof):
        """Return the reduced statistic and q value.

        The reduced statisitc is conceptually simple, as it is just
        statistic / degrees-of-freedom, but it is not meaningful for
        all statistics, and it is only valid if there are any degrees
        of freedom.

        Parameters
        ----------
        statval : float
            The statistic value. It is assumed to be finite.
        dof : int
            The number of degrees of freedom, which may be 0 or negative.

        Returns
        -------
        rstat, qval : float or NaN or None, float or NaN or None
            The reduced statistic and q value. If the statistic does not support
            a goodness of fit then the return values are None. If it does then
            NaN is returned if either the number of degrees of freedom is 0
            (or less), or the statistic value is less than 0.

        """

        if not self._can_calculate_rstat:
            return None, None

        if dof > 0 and statval >= 0.0:
            qval = igamc(dof / 2.0, statval / 2.0)
            rstat = statval / dof
            return rstat, qval

        return numpy.nan, numpy.nan


class Likelihood(Stat):
    """Maximum likelihood function"""

    def __init__(self, name='likelihood'):
        Stat.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        # Likelihood stats do not have 'errors' associated with them.
        # return 1 to avoid dividing by 0 by some optimization methods.
        return numpy.ones_like(data)

    def _check_background_subtraction(self, data):
        """Raise an error if any dataset has been background subtracted.

        Parameters
        ----------
        data : a DataSimulFit instance
            The data sets to use.

        Raises
        ------
        FitErr
        """

        for dobj in data.datasets:
            if getattr(dobj, 'subtracted', False):
                # TODO:
                # Historically this has been a FitErr, but it
                # would make more sense to be a StatErr. This could
                # break people's code, so hold off for now
                # (October 2016).
                #
                raise FitErr('statnotforbackgsub', self.name)

    def _validate_inputs(self, data, model):
        data, model = Stat._validate_inputs(self, data, model)
        self._check_background_subtraction(data)
        return data, model

    def calc_stat(self, data, model):
        fitdata, modeldata = self._get_fit_model_data(data, model)

        return self._calc(fitdata[0], modeldata, None,
                          truncation_value)


# DOC-TODO: where is the truncate/trunc_value stored for objects
#           AHA: it appears to be taken straight from the config
#           file rather than associated with the Stat class
# DOC-TODO: where to talk about the .sherpa.rc config file?


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

    _calc = _statfcts.calc_cash_stat

    def __init__(self, name='cash'):
        Likelihood.__init__(self, name)


class CStat(Likelihood):
    """Maximum likelihood function (XSPEC style).

    This is equivalent to the XSPEC implementation of the
    Cash statistic [1]_ except that it requires a model to be fit
    to the background. To handle the background in the same manner
    as XSPEC, use the WStat statistic.

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

    _calc = _statfcts.calc_cstat_stat
    _can_calculate_rstat = True

    def __init__(self, name='cstat'):
        Likelihood.__init__(self, name)


class Chi2(Stat):
    """Chi Squared statistic.

    The chi-square statistic is:

        chi^2 = (sum)_i [ [ N(i,S) - B(i,x,pB) - S(i,x,pS) ]^2 / sigma(i)^2 ]

    where N(i,S) is the total number of observed counts in bin i of
    the on-source region; B(i,x,pB) is the number of predicted
    background model counts in bin i of the on-source region (zero for
    background-subtracted data), rescaled from bin i of the off-source
    region, and computed as a function of the model argument x(i)
    (e.g., energy or time) and set of background model parameter
    values pB; S(i,x,pS) is the number of predicted source model
    counts in bin i, as a function of the model argument x(i) and set
    of source model parameter values pS; and sigma(i) is the error in
    bin i.

    N(i,B) is the total number of observed counts in bin i of the
    off-source region; A(B) is the off-source "area", which could be
    the size of the region from which the background is extracted, or
    the length of a background time segment, or a product of the two,
    etc.; and A(S) is the on-source "area". These terms may be defined
    for a particular type of data: for example, PHA data sets A(B) to
    `BACKSCAL * EXPOSURE` from the background data set and A(S) to
    `BACKSCAL * EXPOSURE` from the source data set.

    There are different ways of defining the sigma(i) terms,
    supported by the sub-classes.

    Notes
    -----
    It is assumed that there is a one-to-one mapping between a given
    background region bin and a given source region bin.  For
    instance, in the analysis of PHA data, it is assumed that the
    input background counts spectrum is binned in exactly the same way
    as the input source counts spectrum, and any filter applied to the
    source spectrum automatically applied to the background spectrum.
    This means that the user cannot, for example, specify arbitrary
    background and source regions in two dimensions and get correct
    results. This limitation *only* applies to backgrounds included
    included as part of the data set - e.g. as with PHA files - and
    can be avoided by treating the background as a separate data set.

    """

    _calc = _statfcts.calc_chi2_stat
    _can_calculate_rstat = True

    def __init__(self, name='chi2'):
        Stat.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        raise StatErr('chi2noerr')

    def calc_stat(self, data, model):
        fitdata, modeldata = self._get_fit_model_data(data, model)

        return self._calc(fitdata[0], modeldata,
                          fitdata[1], fitdata[2],
                          None,  # TODO: weights
                          truncation_value)

    def calc_chisqr(self, data, model):
        """Return the chi-square value for each bin.

        Parameters
        ----------
        data : a Data or DataSimulFit instance
            The data set, or sets, to use.
        model : a Model or SimulFitModel instance
            The model expression, or expressions. If a SimulFitModel
            is given then it must match the number of data sets in the
            data parameter.

        Returns
        -------
        chisqr : array of numbers
            The per-bin chi-square values.

        """

        _, fvec = self.calc_stat(data, model)
        return fvec * fvec


class LeastSq(Chi2):
    """Least Squared Statistic.

    The least-square statistic is equivalent to a chi-square
    statistic where the error on each point - sigma(i) - is 1.

    """

    _calc = _statfcts.calc_lsq_stat
    _can_calculate_rstat = False

    def __init__(self, name='leastsq'):
        Stat.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        return numpy.ones_like(data)


class Chi2Gehrels(Chi2):
    """Chi Squared with Gehrels variance.

    The variance is estimated from the number of counts in each bin,
    but unlike `Chi2DataVar`, the Gaussian approximation is not
    used. This makes it more-suitable for use with low-count data.

    The standard deviation for each bin is calculated using the
    approximation from [1]_:

        sigma(i,S) = 1 + sqrt(N(i,s) + 0.75)

    where the higher-order terms have been dropped. This is accurate
    to approximately one percent. For data where the background has
    not been subtracted then the error term is:

        sigma(i) = sigma(i,S)

    whereas with background subtraction,

        sigma(i)^2 = sigma(i,S)^2 + [A(S)/A(B)]^2 sigma(i,B)^2

    See Also
    --------
    Chi2DataVar, Chi2ModVar, Chi2XspecVar

    Notes
    -----
    The accuracy of the error term when the background has been
    subtracted has not been determined. A preferable approach to
    background subtraction is to model the background as well as the
    source signal.

    References
    ----------

    .. [1] "Confidence limits for small numbers of events in
           astrophysical data", Gehrels, N. 1986, ApJ, vol 303,
           p. 336-346.
           http://adsabs.harvard.edu/abs/1986ApJ...303..336G

    """

    def __init__(self, name='chi2gehrels'):
        Chi2.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        return _statfcts.calc_chi2gehrels_errors(data)


class Chi2ConstVar(Chi2):
    """Chi Squared with constant variance.

    The variance is the same in each bin, and set to be the mean
    number of counts in the data:

        sigma(i)^2 = (1/N) * (sum)_(j=1)^N N(j,S) + [A(S)/A(B)]^2 N(j,B)

    where N is the number of on-source (and off-source) bins included
    in the fit. The background term appears only if an estimate of the
    background has been subtracted from the data.

    """

    def __init__(self, name='chi2constvar'):
        Chi2.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        return _statfcts.calc_chi2constvar_errors(data)


class Chi2DataVar(Chi2):
    """Chi Squared with data variance.

    The variance in each bin is estimated from the data value in that
    bin.

    If the number of counts in each bin is large, then the shape of
    the Poisson distribution from which the counts are sampled tends
    asymptotically towards that of a Gaussian distribution, with
    variance

        sigma(i)^2 = N(i,S) + [A(S)/A(B)]^2 N(i,B)

    where N is the number of on-source (and off-source) bins included
    in the fit. The background term appears only if an estimate of the
    background has been subtracted from the data.

    See Also
    --------
    Chi2Gehrels, Chi2ModVar, Chi2XspecVar

    """

    def __init__(self, name='chi2datavar'):
        Chi2.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        return _statfcts.calc_chi2datavar_errors(data)


class Chi2ModVar(Chi2):
    """Chi Squared with model amplitude variance.

    The variance in each bin is estimated from the *model* value in
    that bin. This contrasts with other Chi-squared statics - such
    as `Chi2DataVar`, Chi2XspecVar`, and `Chi2Gehrels` - which use
    the data values. The variance is

        sigma(i)^2 = S(i) + [A(S)/A(B)]^2 B(i,off)

    where B(i,off) is the background model amplitude in bin i of the
    off-source region.

    See Also
    --------
    Chi2DataVar, Chi2Gehrels, Chi2XspecVar

    Notes
    -----
    The background should not be subtracted from the data when this
    statistic is used, as it underestimates the variance when fitting
    background-subtracted data.

    """

    _calc = _statfcts.calc_chi2modvar_stat

    def __init__(self, name='chi2modvar'):
        Chi2.__init__(self, name)

    # Statistical errors are not used
    @staticmethod
    def calc_staterror(data):
        return numpy.zeros_like(data)


class Chi2XspecVar(Chi2):
    """Chi Squared with data variance (XSPEC style).

    The variance in each bin is estimated from the data value in that
    bin.

    The calculation of the variance is the same as `Chi2DataVar`
    except that if the number of counts in a bin is less than 1
    then the variance for that bin is set to 1.

    See Also
    --------
    Chi2DataVar, Chi2Gehrels, Chi2ModVar

    """

    def __init__(self, name='chi2xspecvar'):
        Chi2.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        return _statfcts.calc_chi2xspecvar_errors(data)


class UserStat(Stat):
    """Support simple user-supplied statistic calculations.

    Notes
    -----
    This class is used by the `sherpa.ui.load_user_stat`
    to provide a user-definable statistic calculation as
    a function. For more complicated cases it is suggested that
    users should write their own class instead of using
    this one.

    """

    def __init__(self, statfunc=None, errfunc=None, name='userstat'):
        self._statfuncset = False
        self.statfunc = (lambda x: None)

        self._staterrfuncset = False
        self.errfunc = (lambda x: None)

        if statfunc is not None:
            self.statfunc = statfunc
            self._calc = statfunc
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
        self.__dict__['statfunc'] = (lambda x: None)
        self.__dict__['errfunc'] = (lambda x: None)
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

    def calc_stat(self, data, model):
        if not self._statfuncset:
            raise StatErr('nostat', self.name, 'calc_stat()')

        fitdata, modeldata = self._get_fit_model_data(data, model)

        return self.statfunc(fitdata[0],
                             modeldata,
                             staterror=fitdata[1],
                             syserror=fitdata[2],
                             weight=None)  # TODO weights


class WStat(Likelihood):
    """Maximum likelihood function including background (XSPEC style).

    This is equivalent to the XSPEC implementation of the
    W statistic for CStat [1]_, and includes the background data in
    the fit statistic. If a model is being fit to the background then
    the CStat statistic should be used.

    The following description is taken from [1]_.

    Suppose that each bin in the background spectrum is given its own
    parameter so that the background model is b_i = f_i. A standard fit
    for all these parameters would be impractical; however there is an
    analytical solution for the best-fit f_i in terms of the other
    variables which can be derived by using the fact that the derivative
    of the likelihood (L) will be zero at the best fit. Solving for the
    f_i and substituting gives the profile likelihood::

        W = 2 sum_(i=1)^N t_s m_i + (t_s + t_b) f_i -
            S_i ln(t_s m_i + t_s f_i) - B_i ln(t_b f_i) -
            S_i (1- ln(S_i)) - B_i (1 - ln(B_i))

    where::

        f_i = (S_i + B_i - (t_s + t_b) m_i + d_i) / (2 (t_s + t_b))
        d_i = sqrt([(t_s + t_b) m_i - S_i - B_i]^2 +
                   4(t_s + t_b) B_i m_i)

    If any bin has S_i and/or B_i zero then its contribution to W (W_i)
    is calculated as a special case. So, if S_i is zero then::

        W_i = t_s m_i - B_i ln(t_b / (t_s + t_b))

    If B_i is zero then there are two special cases. If
    m_i < S_i / (t_s + t_b) then::

        W_i = - t_b m_i - S_i ln(t_s / (t_s + t_b))

    otherwise::

        W_i = t_s m_i + S_i (ln(S_i) - ln(t_s m_i) - 1)

    In practice, it works well for many cases but for weak sources can
    generate an obviously wrong best fit. It is not clear why this happens
    although binning to ensure that every bin contains at least one count
    often seems to fix the problem. In the limit of large numbers of counts
    per spectrum bin a second-order Taylor expansion shows that W tends to::

        sum_(i=1)^N ( [S_i - t_s m_i - t_s f_i]^2 / (t_s (m_i + f_i)) +
                      [B_i - t_b f_i]^2 / (t_b f_i) )

    which is distributed as chi^2 with N - M degrees of freedom, where
    the model m_i has M parameters (include the normalization).

    References
    ----------

    .. [1] The description of the W statistic (`wstat`) in
           https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html

    """

    _calc = _statfcts.calc_wstat_stat
    _can_calculate_rstat = True

    def __init__(self, name='wstat'):
        Likelihood.__init__(self, name)

    def calc_stat(self, data, model):

        data, model = self._validate_inputs(data, model)

        # Need access to backscal values and background data filtered
        # and grouped in the same manner as the data. There is no
        # easy access to this via the Data API (in part because the
        # Data class has no knowledge of grouping or backscale values).
        #
        # An alternative approach would be to just calculate the
        # statistic for each dataset individually and then
        # sum the statistics for the return value, but the
        # original code used this approach.
        #
        data_src = []
        data_model = data.eval_model_to_fit(model)
        data_bkg = []
        nelems = []
        exp_src = []
        exp_bkg = []
        backscales = []

        for dset, mexpr in zip(data.datasets, model.parts):

            y = dset.to_fit(staterrfunc=None)[0]
            data_src.append(y)
            nelems.append(y.size)

            try:
                bids = dset.background_ids
            except AttributeError:
                raise StatErr('usecstat')

            nbkg = len(bids)
            if nbkg == 0:
                raise StatErr('usecstat')

            elif nbkg > 1:
                # TODO: improve warning
                warnings.warn("Only using first background component for data set {}".format(dset.name))

            bid = bids[0]

            bset = dset.get_background(bid)

            # TODO: the following should be reviewed to see what
            # happens if optional information is missing (e.g. if
            # BACKSCAL is not set we should default to all 1's,
            # but does this code handle it?)
            #

            data_bkg.append(dset.apply_filter(bset.get_dep(False),
                                              groupfunc=numpy.sum))

            # The assumption is that the source and background datasets
            # have the same number of channels (before any grouping or
            # filtering is applied).
            #
            # Since the backscal values can be a scalar or array, it is
            # easiest just to convert everything to an array.
            #
            dummy = numpy.ones(dset.get_dep(False).size)

            # Combine the BACKSCAL values (use the default _middle
            # scheme as this is used elsewhere when combining
            # BACKSCAL values; perhaps there should be an API call
            # for this?).
            #
            src_backscal = dset.apply_filter(dset.backscal * dummy,
                                             groupfunc=dset._middle)
            bkg_backscal = dset.apply_filter(bset.backscal * dummy,
                                             groupfunc=dset._middle)

            backscales.append(bkg_backscal / src_backscal)

            # The AREASCAL values are applied to the exposure
            # times, since this is how XSPEC handles this (at
            # least that's my undertanding of a conversation with
            # Keith Arnaud, for XSPEC ~ version 12.9). This requires
            # turning an exposure into an array if there's no
            # AREASCAl value
            #
            # For now we follow the same approach as the BACKSCAL
            # values if the data is grouped.
            #
            #
            if dset.areascal is None:
                ascal = dummy[:dset.get_dep(True).size]
            else:
                ascal = dset.apply_filter(dset.areascal * dummy,
                                          groupfunc=dset._middle)

            exp_src.append(dset.exposure * ascal)

            if bset.areascal is None:
                ascal = dummy[:dset.get_dep(True).size]
            else:
                ascal = dset.apply_filter(bset.areascal * dummy,
                                          groupfunc=dset._middle)

            exp_bkg.append(bset.exposure * ascal)

        data_src = numpy.concatenate(data_src)
        exp_src = numpy.concatenate(exp_src)
        exp_bkg = numpy.concatenate(exp_bkg)
        data_bkg = numpy.concatenate(data_bkg)
        backscales = numpy.concatenate(backscales)

        return self._calc(data_src, data_model, nelems,
                          exp_src, exp_bkg,
                          data_bkg, backscales,
                          truncation_value)
