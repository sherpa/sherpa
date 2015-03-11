# 
#  Copyright (C) 2011  Smithsonian Astrophysical Observatory
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
from itertools import izip

from sherpa.estmethods import Covariance, Confidence
from sherpa.utils.err import EstErr
from sherpa.utils import parallel_map, NoNewAttributesAfterInit


import logging
warning = logging.getLogger("sherpa").warning


__all__ = ['multivariate_t', 'multivariate_cauchy',
           'normal_sample', 'uniform_sample', 't_sample',
           'ParameterScaleVector','ParameterScaleMatrix',
           'UniformParameterSampleFromScaleVector',
           'NormalParameterSampleFromScaleVector',
           'NormalParameterSampleFromScaleMatrix',
           'StudentTParameterSampleFromScaleMatrix',
           'NormalSampleFromScaleMatrix', 'NormalSampleFromScaleVector',
           'UniformSampleFromScaleVector', 'StudentTSampleFromScaleMatrix',
           ]

def multivariate_t(mean, cov, df, size=None):
    """
    multivariate_t(mean, cov, df[, size])

    Draw random deviates from a multivariate Student's T distribution Such a
    distribution is specified by its mean covariance matrix, and degrees of
    freedom.  These parameters are analogous to the mean (average or "center"),
    variance (standard deviation, or "width," squared), and the degrees of
    freedom of the one-dimensional t distribution.

    Parameters
    ----------
    mean : 1-D array_like, length N
        Mean of the N-dimensional distribution
    cov : 2-D array_like, shape (N, N)
        Covariate matrix of the distribution.  Must be symmetric and
        positive semi-definite for "physically meaningful" results.
    df : int
        Degrees of freedom of the distribution
    size : tuple of ints, optional
        Given a shape of, for example, ``(m,n,k)``, ``m*n*k`` samples are
        generated, and packed in an `m`-by-`n`-by-`k` arrangement.  Because
        each sample is `N`-dimensional, the output shape is ``(m,n,k,N)``.
        If no shape is specified, a single (`N`-D) sample is returned.

    Returns
    -------
    out : ndarray
        The drawn samples, of shape *size*, if that was provided.  If not,
        the shape is ``(N,)``.

        In other words, each entry ``out[i,j,...,:]`` is an N-dimensional
        value drawn from the distribution.

    Is this right?  This needs to be checked!  A reference to the literature
    the better
    """
    df = float(df)
    mean = numpy.asarray(mean)
    normal = numpy.random.multivariate_normal(numpy.zeros_like(mean), cov, size)
    #x = numpy.sqrt(numpy.random.chisquare(df)/df)
    #numpy.divide(normal, x, normal)
    x = numpy.sqrt(numpy.random.chisquare(df,size)/df)
    numpy.divide(normal, x[numpy.newaxis].T, normal)
    numpy.add(mean, normal, normal)
    x = normal
    return x


def multivariate_cauchy(mean, cov, size=None):
    """
    This needs to be checked too! A reference to the literature the better
    """
    return multivariate_t(mean, cov, 1, size=None)



class ParameterScale(NoNewAttributesAfterInit):

    sigma = 1

    def get_scales(self, fit, myscale=None):
        raise NotImplementedError


class ParameterScaleVector(ParameterScale):

    def get_scales(self, fit, myscales=None):

        scales = []
        thawedpars = [par for par in fit.model.pars if not par.frozen]

        if None == myscales:

            oldestmethod = fit.estmethod

            covar = Covariance()
            covar.config['sigma'] = self.sigma
            fit.estmethod = Covariance()

            try:
                r = fit.est_errors()
            finally:
                fit.estmethod = oldestmethod

            for par, val, lo, hi in izip(thawedpars, r.parvals, r.parmins, r.parmaxes):
                scale = None
                if lo is not None and hi is not None:
                    scale = numpy.abs(lo)
                else:
                    warning("Covariance failed for '%s', trying Confidence..." %
                            par.fullname)

                    conf = Confidence()
                    conf.config['sigma'] = self.sigma
                    fit.estmethod = conf
                    try:
                        t = fit.est_errors(parlist = (par,))
                        if t.parmins[0] is not None and t.parmaxes[0] is not None:
                            scale = numpy.abs(t.parmins[0])

                        else:

                            if t.parmins[0] is None and t.parmaxes[0] is not None:
                                scale = numpy.abs(t.parmaxes[0])

                            else:

                                warning('1 sigma bounds for parameter ' +
                                        par.fullname +
                                        ' could not be found, using soft limit minimum')
                                if 0.0 == numpy.abs(par.min):
                                    scale = 1.0e-16
                                else:
                                    scale = numpy.abs(par.min)

                    finally:
                        fit.estmethod = oldestmethod
                scales.append(scale)

        else:
            if not numpy.iterable( myscales ):
                raise TypeError( "scales option must be iterable of length %d " % len( thawedpars ) )
            scales = map( abs, myscales )
        scales = numpy.asarray(scales).transpose()
        return scales


class ParameterScaleMatrix(ParameterScale):


    def get_scales(self, fit, myscales=None):

        def get_size_of_covar( pars ):
            thawedpars = [par for par in pars if not par.frozen]
            npar = len( thawedpars )
            msg = 'scales must be a numpy array of size (%d,%d)' % (npar, npar)
            return npar, msg

        if myscales == None:
            oldestmethod = fit.estmethod
            fit.estmethod = Covariance()

            try:
                r = fit.est_errors()
            finally:
                fit.estmethod = oldestmethod

            cov = r.extra_output

        else:
            if isinstance( myscales, (numpy.ndarray) ):
                npar, msg = get_size_of_covar( fit.model.pars )
                if ( npar, npar ) == myscales.shape:
                    cov = myscales
                else:
                    raise EstErr( msg )
            else:
                npar, msg = get_size_of_covar( fit.model.pars )
                raise EstErr( msg )
            cov = myscales

        if cov is None:
            raise EstErr('nocov')

        cov = numpy.asarray(cov)

        # Investigate spectral decomposition to avoid requirement that the cov be
        # semi-positive definite.  Nevermind, NumPy already uses SVD to generate
        # deviates from a multivariate normal.  An alternative is to use Cholesky
        # decomposition, but it assumes that the matrix is semi-positive definite.
        if numpy.min(numpy.linalg.eigvalsh(cov)) <= 0:
            raise TypeError("The covariance matrix is not positive definite")

        return cov


class ParameterSampleFromScaleVector(NoNewAttributesAfterInit):

    def __init__(self):
        self.scale = ParameterScaleVector()
        NoNewAttributesAfterInit.__init__(self)
    

    def get_sample(self):
        raise NotImplementedError



class ParameterSampleFromScaleMatrix(NoNewAttributesAfterInit):

    def __init__(self):
        self.scale = ParameterScaleMatrix()
        NoNewAttributesAfterInit.__init__(self)
    

    def get_sample(self):
        raise NotImplementedError


class UniformParameterSampleFromScaleVector(ParameterSampleFromScaleVector):

    def get_sample(self, fit, factor=4, num=1):
        vals = numpy.array(fit.model.thawedpars)
        scales = self.scale.get_scales(fit)
        samples = [numpy.random.uniform(val - factor*abs(scale),
                                        val + factor*abs(scale),
                                        int(num)) for val, scale in izip(vals, scales)]
        return numpy.asarray(samples).T


class NormalParameterSampleFromScaleVector(ParameterSampleFromScaleVector):

    def get_sample(self, fit, myscales=None, num=1):
        vals = numpy.array(fit.model.thawedpars)
        scales = self.scale.get_scales(fit, myscales)
        samples = [numpy.random.normal(val, scale, int(num)) for val, scale in izip(vals, scales)]
        return numpy.asarray(samples).T


class NormalParameterSampleFromScaleMatrix(ParameterSampleFromScaleMatrix):

    def get_sample(self, fit, mycov=None, num=1):
        vals = numpy.array(fit.model.thawedpars)
        cov = self.scale.get_scales(fit,mycov)
        return numpy.random.multivariate_normal(vals, cov, int(num))


class StudentTParameterSampleFromScaleMatrix(ParameterSampleFromScaleMatrix):

    def get_sample(self, fit, dof, num=1):
        vals = numpy.array(fit.model.thawedpars)
        cov = self.scale.get_scales(fit)
        return multivariate_t(vals, cov, dof, int(num))


def _sample_stat(fit, samples, numcores=None):

    oldvals = fit.model.thawedpars
 
    def evaluate(sample):
        fit.model.thawedpars = sample
        return fit.calc_stat()

    stats = None
    try:
        fit.model.startup()
        stats = numpy.asarray(parallel_map(evaluate, samples, numcores))
    finally:
        fit.model.teardown()
        fit.model.thawedpars = oldvals

    return numpy.concatenate([stats[:,numpy.newaxis], samples], axis=1)


class NormalSampleFromScaleMatrix(NormalParameterSampleFromScaleMatrix):

    def get_sample(self, fit, num=1, numcores=None):
        samples = NormalParameterSampleFromScaleMatrix.get_sample(self, fit, num=num)
        return _sample_stat(fit, samples, numcores)


class NormalSampleFromScaleVector(NormalParameterSampleFromScaleVector):

    def get_sample(self, fit, num=1, numcores=None):
        samples = NormalParameterSampleFromScaleVector.get_sample(self, fit, num=num)
        return _sample_stat(fit, samples, numcores)


class UniformSampleFromScaleVector(UniformParameterSampleFromScaleVector):

    def get_sample(self, fit, num=1, factor=4, numcores=None):
        samples = UniformParameterSampleFromScaleVector.get_sample(self, fit,
                                                                   factor=factor, num=num)
        return _sample_stat(fit, samples, numcores)


class StudentTSampleFromScaleMatrix(StudentTParameterSampleFromScaleMatrix):

    def get_sample(self, fit, num=1, dof=2, numcores=None):
        samples = StudentTParameterSampleFromScaleMatrix.get_sample(self, fit, dof, num)
        return _sample_stat(fit, samples, numcores)


def normal_sample(fit, num=1, sigma=1, correlate=True, numcores=None):
    """
    Calculate `num` samples of the current thawed parameters from a Normal
    distribution with a spread of `sigma`.

    `correlate` True uses a multi-variate normal on all
                False uses a uni-variate normal for each
    """
    sampler = NormalSampleFromScaleVector()
    sampler.scale.sigma = sigma
    if correlate:
        sampler = NormalSampleFromScaleMatrix()
    return sampler.get_sample(fit, num, numcores)


def uniform_sample(fit, num=1, factor=4, numcores=None):
    """
    Calculate `num` samples of the current thawed parameters from a Uniform
    distribution.
    """
    sampler = UniformSampleFromScaleVector()
    sampler.scale.sigma = 1
    return sampler.get_sample(fit, num, factor, numcores)


def t_sample(fit, num=1, dof=2, numcores=None):
    """
    Calculate `num` samples of the current thawed parameters from a Student's T
    distribution with degrees of freedom `dof`.
    """
    sampler = StudentTSampleFromScaleMatrix()
    return sampler.get_sample(fit, num, dof, numcores)
