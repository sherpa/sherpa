#
#  Copyright (C) 2009, 2015, 2016, 2019, 2020, 2021
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

"""Calculate fluxes for Astronomical data sets.

The functions in this module are designed primarily for analysis
of X-ray data - that is, for use with DataPHA objects - but can
be used with other data classes.

"""

import logging

import numpy
import numpy.random

from sherpa.astro.utils import calc_energy_flux
from sherpa.utils import parallel_map
from sherpa.utils.err import ArgumentErr, FitErr, ModelErr
from sherpa.sim import NormalParameterSampleFromScaleMatrix, \
    NormalParameterSampleFromScaleVector
from sherpa.models.model import SimulFitModel

info = logging.getLogger(__name__).info

__all__ = ['calc_flux', 'sample_flux', 'calc_sample_flux']


class CalcFluxWorker():
    """Internal class for use by calc_flux.

    We always return the full sample, even when only
    a subset of them are needed to calculate the flux.
    """

    def __init__(self, method, data, src, lo, hi, subset=None):
        self.method = method
        self.data = data
        self.src = src
        self.lo = lo
        self.hi = hi
        self.subset = subset

    def __call__(self, sample):
        if self.subset is None:
            self.src.thawedpars = sample
        else:
            self.src.thawedpars = sample[self.subset]

        flux = self.method(self.data, self.src, self.lo, self.hi)
        return numpy.asarray([flux] + list(sample))


def calc_flux(data, src, samples, method=calc_energy_flux,
              lo=None, hi=None, numcores=None, subset=None):
    """Calculate model fluxes from a sample of parameter values.

    Given a set of parameter values, calculate the model flux for
    each set.

    .. versionchanged:: 4.12.2
       The subset parameter was added.

    .. versionchanged:: 4.12.1
       The fit parameter was removed.

    Parameters
    ----------
    data : sherpa.data.Data subclass
        The data object to use.
    src : sherpa.models.Arithmetic instance
        The source model (without instrument response for PHA data)
    samples : 2D array
        The rows indicate each set of sample, and the columns the
        parameter values to use. If there are n free parameters
        in the model then the array must have a size of num by m,
        where num is the number of fluxes to calculate and m >= n.
        If m > n then the subset argument must be set.
    method : function, optional
        How to calculate the flux: assumed to be one of calc_energy_flux
        or calc_photon_flux
    lo : number or None, optional
        The lower edge of the dataspace range for the flux calculation.
        If None then the lower edge of the data grid is used.
    hi : number or None, optional
        The upper edge of the dataspace range for the flux calculation.
        If None then the upper edge of the data grid is used.
    numcores : int or None, optional
        Should the analysis be split across multiple CPU cores?
        When set to None all available cores are used.
    subset : list of ints or None, optional
        This is only used when the samples array has more parameters
        in it than are free in src. In this case the subset array lists
        the column number of the free parameters in src. So, if the
        samples represented 'nh', 'gamma', and 'ampl' values for each
        row, but the src model only contained the 'gamma' and 'ampl'
        parameters then subset would be [1, 2].

    Returns
    -------
    vals : 2D NumPy array
        If the samples array has a shape of (num, nfree) then vals
        has the shame (num, nfree + 1). The first column is the flux
        for the row, and the remaining columns are copies of the input
        samples array.

    See Also
    --------
    sample_flux

    """

    old_vals = src.thawedpars
    worker = CalcFluxWorker(method, data, src, lo, hi, subset)
    try:
        fluxes = parallel_map(worker, samples, numcores)
    finally:
        src.thawedpars = old_vals

    return numpy.asarray(fluxes)


def _sample_flux_get_samples_with_scales(fit, src, correlated, scales,
                                         num, clip='hard'):
    """Return the parameter samples given the parameter scales.

    Parameters
    ----------
    fit : sherpa.fit.Fit instance
        The fit instance. The fit.model expression is assumed to include
        any necessary response information. The number of free parameters
        in fit.model is mfree.
    src : sherpa.models.ArithmeticModel instance
        The model for which the flux is being calculated. This must be
        a subset of the fit.model expression, and should not include the
        response information. There must be at least one thawed parameter
        in this model. The number of free parameters in src is sfree.
    correlated : bool
        Are the parameters assumed to be correlated or not? If correlated
        is True then scales must be 2D.
    scales : 1D or 2D array
        The parameter scales. When 1D they are the gaussian sigma
        values for the parameter, and when a 2D array they are
        the covariance matrix. The scales parameter must match the number
        of parameters in fit (mfree) and not in src (sfree) when they
        are different. For 1D the size is mfree and for 2D it is
        mfree by mfree.
    num : int
        Tne number of samples to return. This must be 1 or greater.
    clip : {'hard', 'soft', 'none'}, optional
        What clipping strategy should be applied to the sampled
        parameters. The default ('hard') is to fix values at their
        hard limits if they exceed them. A value of 'soft' uses the
        soft limits instead, and 'none' applies no clipping. The last
        column in the returned arrays indicates if the row had any
        clipped parameters (even when clip is set to 'none').

    Returns
    -------
    samples, clipped : 2D NumPy array, 1D NumPy array
        The dimensions are num by mfree. The ordering of the parameter
        values in each row matches that of the free parameters in
        fit.model.  The clipped array indicates whether a row had one
        or more clipped parameters.

    Raises
    ------
    ArgumentErr
        If the scales argument contains invalid (e.g. None or IEEE
        non-finite values) values, or is the wrong shape.
    ModelErr
        If the scales argument has the wrong size (that is, it
        does not represent mfree parameter values).

    Notes
    -----
    The support for src being a subset of the fit.model argument
    has not been tested for complex models, that is when fit.model
    is rmf(arf(source_model)) and src is a combination of components
    in source_model but not all the components of source_model.

    """

    npar = len(src.thawedpars)
    mpar = len(fit.model.thawedpars)
    assert mpar >= npar

    scales = numpy.asarray(scales)

    # A None value will cause scales to have a dtype of object,
    # which is not supported by isfinite, so check for this
    # first.
    #
    # Numpy circa 1.11 raises a FutureWarning with 'if None in scales:'
    # about this changing to element-wise comparison (which is what
    # we want). To avoid this warning I use the suggestion from
    # https://github.com/numpy/numpy/issues/1608#issuecomment-9618150
    #
    if numpy.equal(None, scales).any():
        raise ArgumentErr('bad', 'scales',
                          'must not contain None values')

    # We require that scales only has finite values in it.
    # The underlying sample routines are assumed to check other
    # constraints, or deal with negative values (for the 1D case
    # uncorrelated case the absolute value is used).
    #
    if not numpy.isfinite(scales).all():
        raise ArgumentErr('bad', 'scales',
                          'must only contain finite values')

    if scales.ndim == 2 and (scales.shape[0] != scales.shape[1]):
        raise ArgumentErr('bad', 'scales',
                          'scales must be square when 2D')

    # Ensure the scales array matches the correlated parameter:
    #  - when True it must be the covariance matrix (2D)
    #  - when False it can be either be a 1D array of sigmas or
    #    the covariance matrix, which we convert to an array of
    #    sigmas
    #
    if correlated:
        if scales.ndim != 2:
            raise ArgumentErr('bad', 'scales',
                              'when correlated=True, scales must be 2D')
    elif scales.ndim == 2:
        # convert from covariance matrix
        scales = numpy.sqrt(scales.diagonal())
    elif scales.ndim != 1:
        raise ArgumentErr('bad', 'scales',
                          'when correlated=False, scales must be 1D or 2D')

    # At this point either 1D or 2D square array. Now to check the
    # number of elements.
    #
    if scales.shape[0] != mpar:
        raise ModelErr('numthawed', mpar, scales.shape[0])

    if correlated:
        sampler = NormalParameterSampleFromScaleMatrix()
    else:
        sampler = NormalParameterSampleFromScaleVector()

    samples = sampler.get_sample(fit, scales, num=num)
    clipped = sampler.clip(fit, samples, clip=clip)
    return samples, clipped


def _sample_flux_get_samples(fit, src, correlated, num, clip='hard'):
    """Return the parameter samples, using fit to define the scales.

    The covariance method is used to estimate the errors for the
    parameter sampling.

    Parameters
    ----------
    fit : sherpa.fit.Fit instance
        The fit instance. The fit.model expression is assumed to include
        any necessary response information. The number of free parameters
        in the model expression is mfree.
    src : sherpa.models.ArithmeticModel instance
        The model for which the flux is being calculated. This must be
        a subset of the fit.model expression, and should not include the
        response information. There must be at least one thawed parameter
        in this model. The number of free parameters in src is sfree and
        must be <= mfree.
    correlated : bool
        Are the parameters assumed to be correlated or not?
    num : int
        Tne number of samples to return. This must be 1 or greater.
    clip : {'hard', 'soft', 'none'}, optional
        What clipping strategy should be applied to the sampled
        parameters. The default ('hard') is to fix values at their
        hard limits if they exceed them. A value of 'soft' uses the
        soft limits instead, and 'none' applies no clipping. The last
        column in the returned arrays indicates if the row had any
        clipped parameters (even when clip is set to 'none').

    Returns
    -------
    samples, clipped : 2D NumPy array, 1D NumPy array
        The dimensions are num by mfree. The ordering of the parameter
        values in each row matches that of the free parameters in
        fit.model. The clipped array indicates whether a row
        had one or more clipped parameters.

    Notes
    -----
    The support for src being a subset of the fit.model argument
    has not been tested for complex models, that is when fit.model
    is rmf(arf(source_model)) and src is a combination of components
    in source_model but not all the components of source_model.
    """

    npar = len(src.thawedpars)
    mpar = len(fit.model.thawedpars)
    assert mpar >= npar

    if correlated:
        sampler = NormalParameterSampleFromScaleMatrix()
    else:
        sampler = NormalParameterSampleFromScaleVector()

    samples = sampler.get_sample(fit, num=num)
    clipped = sampler.clip(fit, samples, clip=clip)
    return samples, clipped


def decompose(mdl):
    """Return the individual instances in the model expression.

    Parameters
    ----------
    mdl : sherpa.models.model.ArithmeticModel instance
        The model expression. This can be a single component
        or a combination, and can also be a SimulFitModel
        instance.

    Returns
    -------
    cpts : set of sherpa.models.model.Model instances
        The individual components in mdl

    """

    # The __contains__ behavior of the model class appears to
    # return the individual components, including sub-expressions,
    # so we just drop the sub-expressions.
    #
    out = set([])

    if isinstance(mdl, SimulFitModel):
        # This is messier than I'd like
        for mcpt in mdl:
            for cpt in mcpt:
                if hasattr(cpt, 'parts'):
                    continue

                out.add(cpt)

        return out

    for cpt in mdl:
        if hasattr(cpt, 'parts'):
            continue

        out.add(cpt)

    return out


def sample_flux(fit, data, src,
                method=calc_energy_flux, correlated=False,
                num=1, lo=None, hi=None, numcores=None, samples=None,
                clip='hard'):
    """Calculate model fluxes from a sample of parameter values.

    Draw parameter values from a normal distribution and then calculate
    the model flux for each set of parameter values. The values are
    drawn from normal distributions, and the distributions can either
    be independent or have correlations between the parameters.

    .. versionchanged:: 4.12.2
       The clip parameter was added and an extra column is added to
       the return to indicate if each row was clipped.

    Parameters
    ----------
    fit : sherpa.fit.Fit instance
        The fit object. The src parameter is assumed to be a subset of
        the fit.model expression (to allow for calculating the flux of
        a model component), and the fit.data object matches the data
        object. The fit.model argument is expected to include instrumental
        models (for PHA data sets). These objects can represent
        simultaneous fits (e.g. sherpa.data.DataSimulFit and
        sherpa.models.model.SimulFitModel).
    data : sherpa.data.Data subclass
        The data object to use. This is not a DataSimulFit instance.
    src : sherpa.models.Arithmetic instance
        The source model (without instrument response for PHA data) that
        is used for calculating the flux. This is not a SimulFitModel
        instance.
    method : function, optional
        How to calculate the flux: assumed to be one of calc_energy_flux
        or calc_photon_flux
    correlated : bool, optional
        Are the parameter draws independent of each other?
    num : int, optional
        The number of iterations.
    lo : number or None, optional
        The lower edge of the dataspace range for the flux calculation.
        If None then the lower edge of the data grid is used.
    hi : number or None, optional
        The upper edge of the dataspace range for the flux calculation.
        If None then the upper edge of the data grid is used.
    numcores : int or None, optional
        Should the analysis be split across multiple CPU cores?
        When set to None all available cores are used.
    samples : 1D or 2D array, optional
        What are the errors on the parameters? If set to None then
        the covariance method is used to estimate the parameter errors.
        If given and correlated is True then samples must be a 2D array,
        and contain the covariance matrix for the free parameters in
        fit.model, and the matrix must be positive definite. If correlated
        is False then samples can either be sent the covariance matrix or a
        1D array of the error values (i.e. the sigma of the normal
        distribution). If there are n free parameters then the 1D array has
        to have n elements and the 2D array n by n elements.
    clip : {'hard', 'soft', 'none'}, optional
        What clipping strategy should be applied to the sampled
        parameters. The default ('hard') is to fix values at their
        hard limits if they exceed them. A value of 'soft' uses the
        soft limits instead, and 'none' applies no clipping. The last
        column in the returned arrays indicates if the row had any
        clipped parameters (even when clip is set to 'none').

    Returns
    -------
    vals : 2D NumPy array
        The shape of samples is (num, nfree + 2), where nfree is the
        number of free parameters in fit.model. Each row contains one
        iteration, and the columns are the calculated flux, followed
        by the free parameters, and then a flag column indicating if
        the parameters were clipped (1) or not (0).

    See Also
    --------
    calc_flux, calc_sample_flux

    Notes
    -----
    The ordering of the samples array, and the columns in the output,
    matches that of the free parameters in the fit.model expression. That is::

        [p.fullname for p in fit.model.pars if not p.frozen]

    If src is a subset of the full source expression then samples,
    when not None, must still match the number of free parameters in
    the full source expression (that given by fit.model).

    """

    if num <= 0:
        raise ArgumentErr('bad', 'num', 'must be a positive integer')

    # Ensure we have free parameters. Note that this uses the
    # 'src' model, not 'fit.model'.
    #
    npar = len(src.thawedpars)
    if npar == 0:
        raise FitErr('nothawedpar')

    # Check that src is a "subset" of fit.model. The current approach is
    # probably not sufficient to capture the requirements, where
    # ordering of the components is important, but is a start.
    #
    cpts_full = decompose(fit.model)
    for cpt in decompose(src):
        if cpt in cpts_full:
            continue

        raise ArgumentErr('bad', 'src', 'model contains term not in fit')

    mpar = len(fit.model.thawedpars)
    if npar > mpar:
        # This should not be possible given the above check, but leave in.
        raise ArgumentErr('bad', 'src',
                          'more free parameters than expected')

    # The argument to sample_flux should really be called scales and
    # not samples.
    #
    scales = samples
    if scales is None:
        samples, clipped = _sample_flux_get_samples(fit, src, correlated,
                                                    num, clip=clip)
    else:
        samples, clipped = _sample_flux_get_samples_with_scales(fit, src, correlated,
                                                                scales, num, clip=clip)

    # When a subset of the full model is used we need to know how
    # to select which rows in the samples array refer to the
    # parameters of interest. We could compare on fullname,
    # but is not sufficient to guarantee the match.
    #
    if npar < mpar:
        full_pars = dict(map(reversed,
                             enumerate([p for p in fit.model.pars
                                        if not p.frozen])))
        cols = []
        for src_par in [p for p in src.pars if not p.frozen]:
            try:
                cols.append(full_pars[src_par])
            except KeyError:
                # This should not be possible at this point but the
                # decompose check above may be insufficient.
                raise ArgumentErr('bad', 'src',
                                  'unknown parameter "{}"'.format(src_par.fullname))

        cols = numpy.asarray(cols)
        assert cols.size == npar, 'We have lost a parameter somewhere'
    else:
        cols = None

    # Need to append the clipped array (it would be nice to retain
    # the boolean nature of this).
    #
    vals = calc_flux(data, src, samples, method, lo, hi, numcores,
                     subset=cols)
    return numpy.concatenate((vals, numpy.expand_dims(clipped, 1)), axis=1)


def calc_sample_flux(lo, hi, fit, data, samples, modelcomponent,
                     confidence):
    """Given a set of parameter samples, estimate the flux distribution.

    This is similar to sample_flux but returns values for both the full
    and component models (which can be the same). The set of parameter
    samples is searched to remove rows where the parameters lie outside
    the soft limits, and then this set of rows is used to calculate the
    flux of the modelcomponent expression (the flux for the full model
    is taken from the samples argument). The resulting set of fluxes is
    used to calculate the median and the confidence range.

    .. versionchanged:: 4.13.1
       The clipping now includes parameters at the soft limits;
       previously they were excluded which could cause problems with
       parameters for which we can only calculate an upper limit. The
       id and session arguments have been removed and the
       modelcomponent argument can now be set to None. The statistic
       value is now returned for each row, even those that were
       excluded from the flux calculation.

    Parameters
    ----------
    lo : number or None, optional
        The lower edge of the dataspace range for the flux calculation.
        If None then the lower edge of the data grid is used.
    hi : number or None, optional
        The upper edge of the dataspace range for the flux calculation.
        If None then the upper edge of the data grid is used.
    fit : sherpa.fit.Fit instance
        The fit object. The src parameter is assumed to be a subset of
        the fit.model expression (to allow for calculating the flux of
        a model component), and the fit.data object matches the data
        object. The fit.model argument is expected to include instrumental
        models (for PHA data sets). These objects can represent
        simultaneous fits (e.g. sherpa.data.DataSimulFit and
        sherpa.models.model.SimulFitModel).
    data : sherpa.data.Data subclass
        The data object to use. This is not a DataSimulFit instance.
    samples
        The output of sample_flux for the data (assumed to have been
        created with method=calc_energy_flux and clip='hard').
    modelcomponent : sherpa.models.Arithmetic instance or None
        The source model (without instrument response for PHA data)
        that is used for calculating the "unabsorbed" flux. This is
        not a SimulFitModel instance. If None then we just re-use the
        input flux values (first column of samples).
    confidence : number, optional
       The confidence level for the upper and lower values, as a
       percentage (0 to 100). A value of 68.3 will return the
       one-sigma range.

    Returns
    -------
    (fullflux, cptflux, vals)
        The fullflux and cptflux arrays contain the results for the
        full source model and the flux of the `modelcomponent`
        argument (they can be the same). They have three elements and
        give the median value, the value containing 100 - confidence/2
        of the data, and the fraction containing confidence/2 of the
        flux distribution. For the default confidence argument of 68
        this means the last two give the one-sigma upper and lower
        bounds. The vals array has a shape of ``(num+1, N+3)``, where
        ``N`` is the number of free parameters and num is the `num`
        parameter. The rows of this array contain the flux value for
        the iteration (for the full source model), the parameter
        values, a flag indicating whether any parameter in that row
        was clipped to the hard-limits of the parameter, and the
        statistic value for this set of parameters.

    See Also
    --------
    calc_flux, calc_sample_flux

    Notes
    -----
    This function is in the process of being documented and
    re-written, as it contains a lot of redundancy and the semantics
    of some of the options are unclear.

    The summary output displayed by this routine - giving the median
    and confidence ranges - is controlled by the standard Sherpa
    logging instance, and can be hidden by changing the logging to a
    level greater than "INFO" (e.g. with
    `sherpa.utils.logging.SherpaVerbosity`).

    """

    thawedpars = fit.model.thawedpars

    # Check the number of free parameters agrees with the samples argument,
    # noting that each row in samples is <flux> + <free pars> + <clip>. This is
    # a requirement for calling this routine, so is just handled as an
    # assert (as this is not intended for a user-level routine)
    #
    nthawed = len(thawedpars)
    npars = samples.shape[1] - 2
    assert nthawed == npars, (nthawed, npars)

    # Identify any row where a parameter lies outside the min/max range
    softmins = fit.model.thawedparmins
    softmaxs = fit.model.thawedparmaxes

    # We have to use columns 1 to n-1 of samples
    valid = numpy.ones(samples.shape[0], dtype=bool)
    for col, pmin, pmax in zip(samples.T[1:-1], softmins, softmaxs):
        valid &= (col >= pmin) & (col <= pmax)

    nrows = samples.shape[0]
    oflx = samples[:, 0]       # observed/absorbed flux
    if modelcomponent is None:
        iflx = oflx
    else:
        iflx = numpy.zeros(nrows)  # intrinsic/unabsorbed flux

    mystat = numpy.zeros((nrows, 1), dtype=samples.dtype)
    try:
        for nn in range(nrows):
            # Need to extract the subset that contains the parameters
            fit.model.thawedpars = samples[nn, 1:-1]
            if valid[nn] and modelcomponent is not None:
                iflx[nn] = calc_energy_flux(data, modelcomponent, lo=lo, hi=hi)

            mystat[nn, 0] = fit.calc_stat()

    finally:
        fit.model.thawedpars = thawedpars

    # We could avoid calculating the iflx values if modelcomponent is None
    # but it would complicate the code for essentially no time saving.
    #
    hwidth = confidence / 2
    result = []
    for flx in [oflx, iflx]:
        result.append(numpy.percentile(flx[valid],
                                       [50, 50 + hwidth, 50 - hwidth]))

    for lbl, arg in zip(['original model', 'model component'], result):
        med, usig, lsig = arg
        msg = '{} flux = {:g}, + {:g}, - {:g}'.format(lbl, med, usig - med, med - lsig)
        info(msg)

    samples = numpy.concatenate((samples, mystat), axis=1)
    result.append(samples)

    return result
