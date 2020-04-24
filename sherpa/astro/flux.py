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
from sherpa.utils.err import ArgumentErr
from sherpa.sim import NormalParameterSampleFromScaleMatrix, \
    NormalParameterSampleFromScaleVector

__all__ = ['calc_flux', 'sample_flux', 'calc_sample_flux']


class CalcFluxWorker():
    def __init__(self, fit, method, data, src, lo, hi):
        self.fit = fit
        self.method = method
        self.data = data
        self.src = src
        self.lo = lo
        self.hi = hi

    def __call__(self, sample):
        self.fit.model.thawedpars = sample
        flux = self.method(self.data, self.src, self.lo, self.hi)
        return [flux] + list(sample)


def calc_flux(fit, data, src, samples, method=calc_energy_flux,
              lo=None, hi=None, numcores=None):
    """Calculate model fluxes from a sample of parameter values.

    Given a set of parameter values, calculate the model flux for
    each set.

    Parameters
    ----------
    fit : sherpa.fit.Fit instance
        The fit object
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

    def evaluate(sample):
        fit.model.thawedpars = sample
        flux = method(data, src, lo, hi)
        return [flux] + list(sample)

    old_model_vals = fit.model.thawedpars
    try:
        fluxes = parallel_map(CalcFluxWorker(fit, method, data, src, lo, hi), samples, numcores)
    finally:
        fit.model.thawedpars = old_model_vals

    return numpy.asarray(fluxes)


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
        The fit object
    data : sherpa.data.Data subclass
        The data object to use.
    src : sherpa.models.Arithmetic instance
        The source model (without instrument response for PHA data)
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
        parameters in the fit. If correlated is False then samples
        can either be sent the covariance matrix or a 1D array
        of the error values (i.e. the sigma of the normal distribution).
        If there are n free parameters then the 1D array has to have
        n elements and the 2D array n by n elements.

    Returns
    -------
    vals : 2D NumPy array
        The shape of samples is (num, nfree + 1), where nfree is the
        number of free parameters in the model. Each row represents
        an iteration, and the columns are the calculated flux,
        followed by the free parameters.

    See Also
    --------
    calc_flux

    Notes
    -----
    The ordering of the samples array, and the columns in the output,
    matches that of the free parameters in the src parameter. That is::

        [p.fullname for p in src.pars if not p.frozen]

    """

    if num <= 0:
        raise ArgumentErr('bad', 'num', 'must be a positive integer')

    # What is the number of free parameters in the model expression?
    npar = len(src.thawedpars)

    #
    # The following function should be modified to take advantage of numpy
    #
    def within_limits(mysamples, mymins, mymaxs):
        num_par = mysamples.shape[1]
        for row in mysamples:
            for index in range(num_par):
                if row[index] < mymins[index]:
                    row[index] = mymins[index]
                if row[index] > mymaxs[index]:
                    row[index] = mymaxs[index]
        return mysamples

    sampler = NormalParameterSampleFromScaleVector()
    if correlated:
        sampler = NormalParameterSampleFromScaleMatrix()

    # Rename samples to scales as it is confusing here (as it does
    # not refer to the samples drawn from the distribution but the
    # widths of the parameters used to create the samples).
    #
    scales = samples
    if scales is not None:
        scales = numpy.asarray(scales)

        # A None value will cause scales to have a dtype of object,
        # which is not supported by isfinite, so check for this
        # first.
        #
        if None in scales:
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

        if correlated:
            if scales.ndim != 2:
                raise ArgumentErr('bad', 'scales',
                                  'when correlated=True, scales must be 2D')
        else:
            if scales.ndim == 2:
                # convert from covariance matrix
                scales = numpy.sqrt(scales.diagonal())
            elif scales.ndim != 1:
                raise ArgumentErr('bad', 'scales',
                                  'when correlated=False, scales must be 1D or 2D')

        # at this point either 1D or 2D square array
        #
        if scales.shape[0] != npar:
            raise ArgumentErr('bad', 'scales',
                              'does not match number of free parameters')

    samples = sampler.get_sample(fit, scales, num=num)

    hardmins = fit.model._get_thawed_par_hardmins()
    hardmaxs = fit.model._get_thawed_par_hardmaxes()
    samples = within_limits(samples, hardmins, hardmaxs)

    return calc_flux(fit, data, src, samples, method, lo, hi, numcores)


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

        softmins = fit.model._get_thawed_par_mins()
        softmaxs = fit.model._get_thawed_par_maxes()
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
