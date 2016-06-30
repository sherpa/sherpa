from __future__ import print_function
#
#  Copyright (C) 2009, 2015  Smithsonian Astrophysical Observatory
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
import numpy.random
import logging
from sherpa.astro.utils import calc_energy_flux
from sherpa.utils import parallel_map
from sherpa.sim import NormalParameterSampleFromScaleMatrix, \
    NormalParameterSampleFromScaleVector

__all__ = ['calc_flux', 'sample_flux', 'calc_sample_flux']


def calc_flux(fit, data, src, samples, method=calc_energy_flux,
              lo=None, hi=None, numcores=None):

    def evaluate(sample):
        fit.model.thawedpars = sample
        flux = method(data, src, lo, hi)
        return [flux] + list(sample)

    old_model_vals = fit.model.thawedpars
    try:
        fluxes = parallel_map(evaluate, samples, numcores)
    finally:
        fit.model.thawedpars = old_model_vals

    return numpy.asarray(fluxes)


def sample_flux(fit, data, src, method=calc_energy_flux, correlated=False,
                num=1, lo=None, hi=None, numcores=None, samples=None):

    #
    # The following function should be modified to take advantage of numpy
    #
    def within_limits(mysamples, mymins, mymaxs):
        num_par = mysamples.shape[1]
        for row in mysamples:
            for index in xrange(num_par):
                if row[index] < mymins[index]:
                    row[index] = mymins[index]
                if row[index] > mymaxs[index]:
                    row[index] = mymaxs[index]
        return mysamples

    sampler = NormalParameterSampleFromScaleVector()
    if correlated:
        sampler = NormalParameterSampleFromScaleMatrix()

    # If user entered a covariance matrix but wants to run with
    # correlated as false then extract the diagonal elements.
    if correlated == False and samples is not None:
        if numpy.ndarray == type(samples):
            samples = samples.diagonal(0)

    samples = sampler.get_sample(fit, samples, num=num)

    hardmins = fit.model._get_thawed_par_hardmins()
    hardmaxs = fit.model._get_thawed_par_hardmaxes()
    samples = within_limits(samples, hardmins, hardmaxs)

    return calc_flux(fit, data, src, samples, method, lo, hi, numcores)


def calc_sample_flux(id, lo, hi, session, fit, data, samples, modelcomponent,
                     confidence):

    def simulated_pars_within_ranges(mysamples, mysoftmins, mysoftmaxs):
        num = len(mysoftmins)
        for ii in xrange(num):
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
        for nn in xrange(size):
            logger.setLevel(logging.ERROR)
            session.set_source(id, orig_source)
            logger.setLevel(orig_log_level)
            oflx[nn] = mysim[nn, 0]
            for ii in xrange(len(thawedpars)):
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
            upconfidence = sf[(1.0 - myconfidence) * size - 1]
            loconfidence = sf[myconfidence * size - 1]
            result.append(numpy.array([median, upconfidence, loconfidence]))

        print_sample_result('original model flux', result[0])
        print_sample_result('model component flux', result[1])

        sampletmp = numpy.zeros((samples.shape[0], 1), dtype=samples.dtype)
        samples = numpy.concatenate((samples, sampletmp), axis=1)

        for index in xrange(size):
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
