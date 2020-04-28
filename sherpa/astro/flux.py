#
#  Copyright (C) 2009, 2015, 2016, 2019, 2020
#       Smithsonian Astrophysical Observatory
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

import numpy
import numpy.random
import logging
from sherpa.astro.utils import calc_energy_flux
from sherpa.utils import parallel_map
from sherpa.utils.err import ArgumentErr, FitErr, ModelErr
from sherpa.sim import NormalParameterSampleFromScaleMatrix, \
    NormalParameterSampleFromScaleVector

__all__ = ['calc_flux', 'sample_flux', 'calc_sample_flux']


class CalcFluxWorker():
    """Internal class for use by calc_flux."""

    def __init__(self, method, data, src, lo, hi):
        self.method = method
        self.data = data
        self.src = src
        self.lo = lo
        self.hi = hi

    def __call__(self, sample):
        self.src.thawedpars = sample
        flux = self.method(self.data, self.src, self.lo, self.hi)
        return numpy.asarray([flux] + list(sample))


def calc_flux(data, src, samples, method=calc_energy_flux,
              lo=None, hi=None, numcores=None):
    """Calculate model fluxes from a sample of parameter values.

    Given a set of parameter values, calculate the model flux for
    each set.

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
        in the model then the array must have a size of num by n,
        where num is the number of fluxes to calculate.
    method : function, optional
        How to calculate the flux: assumed to be one of calc_energy_flux
        or calc_photon_flux
    lo : number or None, optional
        The lower edge of the dataspace range for the flux calculation.
        If None then the lower edge of the data grid is used.
    hi : number or None, optional
        The upper edge of the dataspace range for the flux calculation.
        If None then the upper edge of the data grid is used.
    numcores : int or None, optonal
        Should the analysis be split across multiple CPU cores?
        When set to None all available cores are used.

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

    Notes
    -----
    The ordering of the samples array matches that of the free
    parameters in the src parameter. That is::

        [p.fullname for p in src.pars if not p.frozen]

    """

    old_vals = src.thawedpars
    npar = len(old_vals)
    if npar != samples.shape[1]:
        raise ModelErr('numthawed', npar, samples.shape[1])

    worker = CalcFluxWorker(method, data, src, lo, hi)
    try:
        fluxes = parallel_map(worker, samples, numcores)
    finally:
        src.thawedpars = old_vals

    return numpy.asarray(fluxes)


def _sample_flux_get_samples_with_scales(fit, src, correlated, scales, num):
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
        the covariance matrix. The scales parameter can either represent
        the free parameters in fit.model or src; that is have size
        mfree or sfree (1D) or mfree by mfree or sfree by sfree (2D).
    num : int
        Tne number of samples to return. This must be 1 or greater.

    Returns
    -------
    samples : 2D NumPy array
        The dimensions are num by sfree. The ordering of the parameter
        values in each row matches that of src.

    Raises
    ------
    ArgumentErr
        If the scales argument contains invalid (e.g. None or IEEE
        non-finite values) values, or is the wrong shape.
    ModelErr
        If the scales argument has the wrong size (that is, it
        does not represent either sfree or mfree parameter values).

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

    # At this point either 1D or 2D square array. We require
    # either npar or mpar values. This unfortunately doesn't
    # work with the expected semantics of the numthawed error
    # type (which expects an integer), so use npar as the
    # error message for now.
    #
    if scales.shape[0] not in [npar, mpar]:
        raise ModelErr('numthawed', npar, scales.shape[0])

    # The scales argument to the sampler.get_sample method has to
    # match the size of the number of free parameters in the
    # fit.model expression.
    #
    # Possible conditions (npar > mpar is impossible):
    # a) npar == mpar
    # b) npar < mpar, scales.shape[0] = npar
    # c) npar < mpar, scales.shape[0] = mpar
    #
    # case a is easy, as nothing has to be done.
    #
    # case b is handled by temporarily freezing the "extra" free
    # parameters in the fit.model expression.
    #
    # case c is handled by post-processing the samples to remove
    # the "extra" parameters.
    #
    # At this point we requre that src be a subset of fit.model
    #
    frozen_pars = []
    subset = []
    if npar < mpar:
        if scales.shape[0] == npar:
            # Use a set rather than list to potentially improve
            # the 'in' check below, although the number of elements
            # is not assumed to be large.
            #
            spars = {p for p in src.pars if not p.frozen}
            for p in fit.model.pars:
                if p.frozen or p in spars:
                    continue

                frozen_pars.append(p)

        else:
            mpars = [p for p in fit.model.pars if not p.frozen]
            spars = [p for p in src.pars if not p.frozen]
            mpars = dict(zip(mpars, range(mpar)))
            subset = [mpars[p] for p in spars]

    if correlated:
        sampler = NormalParameterSampleFromScaleMatrix()
    else:
        sampler = NormalParameterSampleFromScaleVector()

    try:
        for p in frozen_pars:
            p.frozen = True

        samples = sampler.get_sample(fit, scales, num=num)
    finally:
        for p in frozen_pars:
            p.frozen = False

    if subset != []:
        samples = samples[:, subset]

    return samples


def _sample_flux_get_samples(fit, src, correlated, num):
    """Return the parameter samples, using fit to define the scales.

    The covariance method is used to estimate the errors for the
    parameter sampling.

    Parameters
    ----------
    fit : sherpa.fit.Fit instance
        The fit instance. The fit.model expression is assumed to include
        any necessary response information.
    src : sherpa.models.ArithmeticModel instance
        The model for which the flux is being calculated. This must be
        a subset of the fit.model expression, and should not include the
        response information. There must be at least one thawed parameter
        in this model. The number of free parameters in src is sfree.
    correlated : bool
        Are the parameters assumed to be correlated or not? If correlated
        is True then scales must be 2D.
    num : int
        Tne number of samples to return. This must be 1 or greater.

    Returns
    -------
    samples : 2D NumPy array
        The dimensions are num by sfree. The ordering of the parameter
        values in each row matches that of src.

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

    # Do we need to remove the free parameters that are in
    # fit.model but not in src from the samples?
    #
    if npar < mpar:
        mpars = [p for p in fit.model.pars if not p.frozen]
        spars = [p for p in src.pars if not p.frozen]
        mpars = dict(zip(mpars, range(mpar)))
        subset = [mpars[p] for p in spars]

        samples = samples[:, subset]

    return samples


def decompose(mdl):
    """Return the individual instances in the model expression.

    Parameters
    ----------
    mdl : sherpa.models.model.ArithmeticModel instance
        The model expression. This can be a single component
        or a combination.

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
    for cpt in mdl:
        if hasattr(cpt, 'parts'):
            continue

        out.add(cpt)

    return out


def sample_flux(fit, data, src,
                method=calc_energy_flux, correlated=False,
                num=1, lo=None, hi=None, numcores=None, samples=None):
    """Calculate model fluxes from a sample of parameter values.

    Draw parameter values from a normal distribution and then calculate
    the model flux for each set of parameter values. The values are
    drawn from normal distributions, and the distributions can either
    be independent or have correlations between the parameters.

    Parameters
    ----------
    fit : sherpa.fit.Fit instance
        The fit object. The src parameter is assumed to be a subset of
        the fit.model expression (to allow for calculating the flux of
        a model component), and the fit.data object matches the data
        object. The fit.model argument is expected to include instrumental
        models (for PHA data sets).
    data : sherpa.data.Data subclass
        The data object to use.
    src : sherpa.models.Arithmetic instance
        The source model (without instrument response for PHA data) that
        is used for calculating the flux.
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
    numcores : int or None, optonal
        Should the analysis be split across multiple CPU cores?
        When set to None all available cores are used.
    samples : 1D or 2D array, optional
        What are the errors on the parameters? If set to None then
        the covariance method is used to estimate the parameter errors.
        If given and correlated is True then samples must be a
        2D array, and contain the covariance matrix for the free
        parameters in src. If correlated is False then samples
        can either be sent the covariance matrix or a 1D array
        of the error values (i.e. the sigma of the normal distribution).
        If there are n free parameters then the 1D array has to have
        n elements and the 2D array n by n elements.

    Returns
    -------
    vals : 2D NumPy array
        The shape of samples is (num, nfree + 1), where nfree is the
        number of free parameters in src. Each row contains one
        iteration, and the columns are the calculated flux,
        followed by the free parameters.

    See Also
    --------
    calc_flux

    Notes
    -----
    The ordering of the samples array, and the columns in the output,
    matches that of the free parameters in the src expression. That is::

        [p.fullname for p in src.pars if not p.frozen]

    If src is a subset of the full source expression then samples,
    when not None, can either match the number of free parameters in
    the full source expression (that given by fit.model), or the src
    parameter.
    """

    if num <= 0:
        raise ArgumentErr('bad', 'num', 'must be a positive integer')

    # Ensure we have free parameters
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
        samples = _sample_flux_get_samples(fit, src, correlated, num)
    else:
        samples = _sample_flux_get_samples_with_scales(fit, src, correlated, scales, num)

    # Ensure that samples falls within the hard limits by
    # clipping the values.
    #
    hardmins = src.thawedparhardmins
    hardmaxs = src.thawedparhardmaxes
    for pvals, pmin, pmax in zip(samples.T, hardmins, hardmaxs):
        # do the clipping in place
        numpy.clip(pvals, pmin, pmax, out=pvals)

    return calc_flux(data, src, samples, method, lo, hi, numcores)


def calc_sample_flux(id, lo, hi, session, fit, data, samples, modelcomponent,
                     confidence):

    def simulated_pars_within_ranges(mysamples, mysoftmins, mysoftmaxs):
        num = len(mysoftmins)
        for ii in range(num):
            ii1 = ii + 1
            tmp = (mysamples[:, ii1] > mysoftmins[ii])
            selpmasym = mysamples[tmp]
            tmp = (selpmasym[:, ii1] < mysoftmaxs[ii])
            mysamples = selpmasym[tmp]
        return mysamples

    def print_sample_result(title, arg):

        print('%s = %g, + %g, - %g' % (title, arg[0], arg[1] - arg[0],
                                       arg[0] - arg[2]))
    #
    # For later restoration
    #
    orig_model = session.get_model(id)
    orig_model_vals = fit.model.thawedpars

    orig_source = session.get_source(id)
    orig_source_vals = orig_source.thawedpars

    try:

        softmins = fit.model.thawedparmins
        softmaxs = fit.model.thawedparmaxes
        mysim = simulated_pars_within_ranges(samples, softmins, softmaxs)

        size = len(mysim[:, 0])
        oflx = numpy.zeros(size)  # observed/absorbed flux
        iflx = numpy.zeros(size)  # intrinsic/unabsorbed flux
        thawedpars = [par for par in fit.model.pars if not par.frozen]

        logger = logging.getLogger("sherpa")
        orig_log_level = logger.level

        mystat = []
        for nn in range(size):
            logger.setLevel(logging.ERROR)
            session.set_source(id, orig_source)
            logger.setLevel(orig_log_level)
            oflx[nn] = mysim[nn, 0]
            for ii in range(len(thawedpars)):
                val = mysim[nn, ii + 1]
                session.set_par(thawedpars[ii].fullname, val)
            session.set_source(id, modelcomponent)
            iflx[nn] = session.calc_energy_flux(lo, hi)
            #####################################
            session.set_full_model(id, orig_model)
            mystat.append(session.calc_stat())
            #####################################
        oflxiflx = [oflx, iflx]

        myconfidence = (1.0 - confidence / 100.0) / 2.0
        result = []

        for x in oflxiflx:
            sf = numpy.sort(x)
            median = numpy.median(sf)
            upconfidence_index = int((1.0 - myconfidence) * size - 1)
            loconfidence_index = int(myconfidence * size - 1)
            upconfidence = sf[upconfidence_index]
            loconfidence = sf[loconfidence_index]
            result.append(numpy.array([median, upconfidence, loconfidence]))

        print_sample_result('original model flux', result[0])
        print_sample_result('model component flux', result[1])

        sampletmp = numpy.zeros((samples.shape[0], 1), dtype=samples.dtype)
        samples = numpy.concatenate((samples, sampletmp), axis=1)

        for index in range(size):
            samples[index][-1] = mystat[index]

        #samples = numpy.delete( samples, (size), axis=0 )
        result.append(samples)

        return result

    finally:

        session.set_full_model(id, orig_model)
        fit.model.thawedpars = orig_model_vals

        logger.setLevel(logging.ERROR)
        session.set_source(id, orig_source)
        logger.setLevel(orig_log_level)
