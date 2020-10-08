#
#  Copyright (C) 2011, 2015, 2016, 2019, 2020  Smithsonian Astrophysical Observatory
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

import logging

import numpy
import numpy.random

from sherpa.estmethods import Covariance, Confidence
from sherpa.utils.err import EstErr
from sherpa.utils import parallel_map, NoNewAttributesAfterInit

warning = logging.getLogger("sherpa").warning


__all__ = ('multivariate_t', 'multivariate_cauchy',
           'normal_sample', 'uniform_sample', 't_sample',
           'ParameterScaleVector', 'ParameterScaleMatrix',
           'UniformParameterSampleFromScaleVector',
           'NormalParameterSampleFromScaleVector',
           'NormalParameterSampleFromScaleMatrix',
           'StudentTParameterSampleFromScaleMatrix',
           'NormalSampleFromScaleMatrix', 'NormalSampleFromScaleVector',
           'UniformSampleFromScaleVector', 'StudentTSampleFromScaleMatrix',
           )


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
    normal = numpy.random.multivariate_normal(
        numpy.zeros_like(mean), cov, size)
    # x = numpy.sqrt(numpy.random.chisquare(df)/df)
    # numpy.divide(normal, x, normal)
    x = numpy.sqrt(numpy.random.chisquare(df, size) / df)
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
    """Create the scaling used to generate parameters.

    The scaling generally refers to an error value (defaulting
    to one sigma) for each parameter.

    """

    # The sigma value to use
    sigma = 1

    def get_scales(self, fit, myscales=None):
        """Return the samples.

        Parameters
        ----------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to
            generate the samples, along with any possible error
            analysis.
        myscales : numpy array or None, optional
            The scales to use. If None then they are
            calculated from the fit.

        Returns
        -------
        scales : numpy array
            The scales array (npar elements, matching the free
            parameters in fit). It may be multi-dimensional.

        """
        raise NotImplementedError


class ParameterScaleVector(ParameterScale):
    """Uncorrelated errors for the parameters.

    """

    def get_scales(self, fit, myscales=None):
        """Return the samples.

        Parameters
        ----------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to
            generate the samples, along with any possible error
            analysis.
        myscales : numpy array or None, optional
            The scales to use. If None then they are
            calculated from the fit.

        Returns
        -------
        scales : numpy array
            One-dimensional scales array (npar elements, matching the
            free parameters in fit). The values are the sigma errors
            for the parameters (or the input values if given).

        """

        scales = []
        thawedpars = [par for par in fit.model.pars if not par.frozen]

        if myscales is None:

            oldestmethod = fit.estmethod

            covar = Covariance()
            covar.config['sigma'] = self.sigma
            fit.estmethod = Covariance()

            try:
                r = fit.est_errors()
            finally:
                fit.estmethod = oldestmethod

            for par, val, lo, hi in zip(thawedpars, r.parvals, r.parmins, r.parmaxes):
                scale = None
                if lo is not None and hi is not None:
                    scale = numpy.abs(lo)
                else:
                    wmsg = "Covariance failed for '{}',".format(par.fullname) + \
                           " trying Confidence..."
                    warning(wmsg)

                    conf = Confidence()
                    conf.config['sigma'] = self.sigma
                    fit.estmethod = conf
                    try:
                        t = fit.est_errors(parlist=(par,))
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
            if not numpy.iterable(myscales):
                emsg = "scales option must be iterable of " + \
                       "length {}".format(len(thawedpars))
                raise TypeError(emsg)
            scales = list(map(abs, myscales))
        scales = numpy.asarray(scales).transpose()
        return scales


class ParameterScaleMatrix(ParameterScale):
    """Correlated errors for the parameters.

    """

    def get_scales(self, fit, myscales=None):
        """Return the samples.

        Parameters
        ----------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to
            generate the samples, along with any possible error
            analysis.
        myscales : numpy array or None, optional
            The scales to use. If None then they are
            calculated from the fit.

        Returns
        -------
        scales : numpy array
            Two-dimensional scales array (npar by npar elements,
            matching the free parameters in fit). The values are the
            covariance matrix for the parameters (or the input values if
            given).

        """

        if myscales is None:
            oldestmethod = fit.estmethod
            fit.estmethod = Covariance()

            try:
                r = fit.est_errors()
            finally:
                fit.estmethod = oldestmethod

            cov = r.extra_output

        else:

            thawedpars = [par for par in fit.model.pars if not par.frozen]
            npar = len(thawedpars)
            msg = 'scales must be a numpy array of size ({0},{0})'.format(npar)

            if not isinstance(myscales, numpy.ndarray):
                raise EstErr(msg)

            if (npar, npar) != myscales.shape:
                raise EstErr(msg)

            cov = myscales

        if cov is None:
            raise EstErr('nocov')

        cov = numpy.asarray(cov)

        # Investigate spectral decomposition to avoid requirement that the cov be
        # semi-positive definite.  Nevermind, NumPy already uses SVD to generate
        # deviates from a multivariate normal.  An alternative is to use Cholesky
        # decomposition, but it assumes that the matrix is semi-positive
        # definite.
        if numpy.min(numpy.linalg.eigvalsh(cov)) <= 0:
            raise TypeError("The covariance matrix is not positive definite")

        return cov


class ParameterSample(NoNewAttributesAfterInit):
    """Create a set of parameter samples.

    """

    def get_sample(self, fit, num=1):
        """Return the samples.

        Parameters
        ----------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to generate
            the samples, along with any possible error analysis.
        num : int, optional
            The number of samples to return.

        Returns
        -------
        samples : 2D numpy array
            The array is num by npar size, where npar is the number of
            free parameters in the fit argument.

        """
        raise NotImplementedError

    def clip(self, fit, samples, clip='none'):
        """Clip the samples if out of bounds.

        Parameter
        --------
        fit : sherpa.fit.Fit instance
            Contains the thawed parameters used to generate the
            samples.
        samples : 2D numpy array
            The samples array, stored as a n by npar matrix. This
            array is changed in place.
        clip : {'none', 'hard', 'soft'} optional
            How should the values be clipped? The default ('none') has no
            clipping. The other methods restrict the values to lie within
            the hard or soft limits of the parameters.

        Returns
        -------
        clipped : 1D numpy array
            The clipped samples (may be unchanged) and a 1D boolean
            array indicating whether any sample in a row was clipped.

        """

        niter = samples.shape[0]
        clipped = numpy.zeros(niter, dtype=numpy.bool)
        if clip == 'none':
            return clipped

        # Values are clipped to lie within mins/maxs (inclusive)
        #
        if clip == 'hard':
            mins = fit.model.thawedparhardmins
            maxs = fit.model.thawedparhardmaxes
        elif clip == 'soft':
            mins = fit.model.thawedparmins
            maxs = fit.model.thawedparmaxes
        else:
            raise ValueError('invalid clip argument: sent {}'.format(clip))

        for pvals, pmin, pmax in zip(samples.T, mins, maxs):
            porig = pvals.copy()

            # do the clipping in place
            numpy.clip(pvals, pmin, pmax, out=pvals)

            # update the clipped array (which is True if a
            # value on the row has been clipped).
            #
            clipped |= (pvals != porig)

        return clipped


class ParameterSampleFromScaleVector(ParameterSample):
    """Samples drawn from uncorrelated parameters.
    """

    def __init__(self):
        self.scale = ParameterScaleVector()
        ParameterSample.__init__(self)


class ParameterSampleFromScaleMatrix(ParameterSample):
    """Samples drawn from correlated parameters.
    """

    def __init__(self):
        self.scale = ParameterScaleMatrix()
        ParameterSample.__init__(self)


class UniformParameterSampleFromScaleVector(ParameterSampleFromScaleVector):
    """Use a uniform distribution to sample parameters.

    The parameters are drawn from a uniform distribution which is set
    to `factor` times the parameter error (the lower bound is included
    but the upper bound is not).
    """

    def get_sample(self, fit, factor=4, num=1):
        """Return the parameter samples.

        Parameter
        ---------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to generate
            the samples, along with any possible error analysis.
        factor : number, optional
            The half-width of the uniform distribution is factor times
            the one-sigma error.
        num : int, optional
            The number of samples to return.

        Returns
        -------
        samples : 2D numpy array
            The array is num by npar size, where npar is the number of
            free parameters in the fit argument.

        """
        vals = numpy.array(fit.model.thawedpars)
        scales = self.scale.get_scales(fit)
        samples = [numpy.random.uniform(val - factor * abs(scale),
                                        val + factor * abs(scale),
                                        int(num)) for val, scale in zip(vals, scales)]
        return numpy.asarray(samples).T


class NormalParameterSampleFromScaleVector(ParameterSampleFromScaleVector):
    """Use a normal distribution to sample parameters (uncorrelated),

    The parameters are drawn from a normal distribution based on the
    parameter errors, and do not include any correlations between the
    parameters. The errors will be generated from the fit object or
    specified directly.

    """

    def get_sample(self, fit, myscales=None, num=1):
        """Return the parameter samples.

        Parameter
        ---------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to generate
            the samples, along with any possible error analysis.
        myscales : 1D numpy array or None, optional
            The error values (one sigma values) to use. If None then
            it is calculated from the fit object.
        num : int, optional
            The number of samples to return.

        Returns
        -------
        samples : 2D numpy array
            The array is num by npar size, where npar is the number of
            free parameters in the fit argument.

        """
        vals = numpy.array(fit.model.thawedpars)
        scales = self.scale.get_scales(fit, myscales)
        samples = [numpy.random.normal(
            val, scale, int(num)) for val, scale in zip(vals, scales)]
        return numpy.asarray(samples).T


class NormalParameterSampleFromScaleMatrix(ParameterSampleFromScaleMatrix):
    """Use a normal distribution to sample parameters (correlated),

    The parameters are drawn from a normal distribution based on the
    parameter errors, and include the correlations between the
    parameters. The errors will be generated from the fit object or
    specified directly as a covariance matrix.

    """

    def get_sample(self, fit, mycov=None, num=1):
        """Return the parameter samples.

        Parameter
        ---------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to generate
            the samples, along with any possible error analysis.
        mycov : 2D numpy array or None, optional
            The covariance matrix to use. If None then it is
            calculated from the fit object.
        num : int, optional
            The number of samples to return.

        Returns
        -------
        samples : 2D numpy array
            The array is num by npar size, where npar is the number of
            free parameters in the fit argument.

        """
        vals = numpy.array(fit.model.thawedpars)
        cov = self.scale.get_scales(fit, mycov)
        return numpy.random.multivariate_normal(vals, cov, int(num))


class StudentTParameterSampleFromScaleMatrix(ParameterSampleFromScaleMatrix):
    """Use a student's t-distribution to sample parameters (correlated),

    The parameters are drawn from a normal distribution based on the
    parameter errors, and include the correlations between the
    parameters. The errors will be generated from the fit object or
    specified directly as a covariance matrix.

    """

    def get_sample(self, fit, dof, num=1):
        """Return the parameter samples.

        Parameter
        ---------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to generate
            the samples, along with any possible error analysis.
        dof : int
            The degrees of freedom of the distribution.
        num : int, optional
            The number of samples to return.

        Returns
        -------
        samples : 2D numpy array
            The array is num by npar size, where npar is the number of
            free parameters in the fit argument.

        """
        vals = numpy.array(fit.model.thawedpars)
        cov = self.scale.get_scales(fit)
        return multivariate_t(vals, cov, dof, int(num))


class Evaluate():
    """
    Callable class for _sample_stat multiprocessing call
    This class used to be a nested function, which can't be pickled and results in
    python3 failing to execute the code.

    Note that this does not guarantee to reset the model
    parameters after being run.
    """
    def __init__(self, fit):
        self.fit = fit

    def __call__(self, sample):
        self.fit.model.thawedpars = sample
        return self.fit.calc_stat()


def _sample_stat(fit, samples, numcores=None, cache=True):
    """Calculate the statistic for each set of samples.

    Parameter
    ----------
    fit : sherpa.fit.Fit instance
        This defines the thawed parameters that are used to generate
        the samples, along with any possible error analysis.
    samples : 2D numpy array
        The samples array, stored as a npar by niter matrix.
    numcores : int or None, optional
        Should the calculation be done on multiple CPUs?  The default
        (None) is to rely on the parallel.numcores setting of the
        configuration file.
    cache : bool, optional
        Should the model cache be used?

    Returns
    -------
    vals : 2D numpy array
        A copy of the samples input with an extra row added to its
        start, giving the statistic value for that row.

    """

    oldvals = fit.model.thawedpars

    try:
        fit.model.startup(cache=cache)
        stats = numpy.asarray(parallel_map(Evaluate(fit), samples, numcores))
    finally:
        fit.model.teardown()
        fit.model.thawedpars = oldvals

    return numpy.concatenate([stats[:, numpy.newaxis], samples], axis=1)


class NormalSampleFromScaleMatrix(NormalParameterSampleFromScaleMatrix):
    """Use a normal distribution to sample statistic and parameters (correlated),

    The parameters are drawn from a normal distribution based on the
    parameter errors, and include the correlations between the
    parameters. The errors will be generated from the fit object or
    specified directly as a covariance matrix.

    """

    def get_sample(self, fit, num=1, numcores=None):
        """Return the statistic and parameter samples.

        Parameter
        ---------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to generate
            the samples, along with any possible error analysis.
        num : int, optional
            The number of samples to return.
        numcores : int or None, optional
            Should the calculation be done on multiple CPUs?
            The default (None) is to rely on the parallel.numcores
            setting of the configuration file.

        Returns
        -------
        samples : 2D numpy array
            The array is num by (npar + 1) size, where npar is the
            number of free parameters in the fit argument. The first
            element in each row is the statistic value, and the
            remaining are the parameter values.

        """

        samples = NormalParameterSampleFromScaleMatrix.get_sample(
            self, fit, num=num)
        return _sample_stat(fit, samples, numcores)


class NormalSampleFromScaleVector(NormalParameterSampleFromScaleVector):
    """Use a normal distribution to sample statistic and parameters (uncorrelated),

    The parameters are drawn from a normal distribution based on the
    parameter errors, and do not include any correlations between the
    parameters. The errors will be generated from the fit object or
    specified directly as a covariance matrix.

    """

    def get_sample(self, fit, num=1, numcores=None):
        """Return the statistic and parameter samples.

        Parameter
        ---------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to generate
            the samples, along with any possible error analysis.
        num : int, optional
            The number of samples to return.
        numcores : int or None, optional
            Should the calculation be done on multiple CPUs?
            The default (None) is to rely on the parallel.numcores
            setting of the configuration file.

        Returns
        -------
        samples : 2D numpy array
            The array is num by (npar + 1) size, where npar is the
            number of free parameters in the fit argument. The first
            element in each row is the statistic value, and the
            remaining are the parameter values.

        """
        samples = NormalParameterSampleFromScaleVector.get_sample(
            self, fit, num=num)
        return _sample_stat(fit, samples, numcores)


class UniformSampleFromScaleVector(UniformParameterSampleFromScaleVector):
    """Use a uniform distribution to sample statistic and parameters.

    The parameters are drawn from a uniform distribution which is set
    to `factor` times the parameter error (the lower bound is included
    but the upper bound is not).
    """

    def get_sample(self, fit, num=1, factor=4, numcores=None):
        """Return the statistic and parameter samples.

        Parameter
        ---------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to generate
            the samples, along with any possible error analysis.
        num : int, optional
            The number of samples to return.
        factor : number, optional
            The half-width of the uniform distribution is factor times
            the one-sigma error.
        numcores : int or None, optional
            Should the calculation be done on multiple CPUs?
            The default (None) is to rely on the parallel.numcores
            setting of the configuration file.

        Returns
        -------
        samples : 2D numpy array
            The array is num by (npar + 1) size, where npar is the
            number of free parameters in the fit argument. The first
            element in each row is the statistic value, and the
            remaining are the parameter values.

        """
        samples = UniformParameterSampleFromScaleVector.get_sample(self, fit,
                                                                   factor=factor, num=num)
        return _sample_stat(fit, samples, numcores)


class StudentTSampleFromScaleMatrix(StudentTParameterSampleFromScaleMatrix):
    """Use a student's t-distribution to sample statistic and parameters (correlated),

    The parameters are drawn from a normal distribution based on the
    parameter errors, and include the correlations between the
    parameters. The errors will be generated from the fit object or
    specified directly as a covariance matrix.

    """

    def get_sample(self, fit, num=1, dof=2, numcores=None):
        """Return the statistic and parameter samples.

        Parameter
        ---------
        fit : sherpa.fit.Fit instance
            This defines the thawed parameters that are used to generate
            the samples, along with any possible error analysis.
        num : int, optional
            The number of samples to return.
        dof : int
            The degrees of freedom of the distribution.
        numcores : int or None, optional
            Should the calculation be done on multiple CPUs?
            The default (None) is to rely on the parallel.numcores
            setting of the configuration file.

        Returns
        -------
        samples : 2D numpy array
            The array is num by (npar + 1) size, where npar is the
            number of free parameters in the fit argument. The first
            element in each row is the statistic value, and the
            remaining are the parameter values.

        """
        samples = StudentTParameterSampleFromScaleMatrix.get_sample(
            self, fit, dof, num)
        return _sample_stat(fit, samples, numcores)


def normal_sample(fit, num=1, sigma=1, correlate=True, numcores=None):
    """Sample the fit statistic by taking the parameter values
    from a normal distribution.

    For each iteration (sample), change the thawed parameters by
    drawing values from a uni- or multi-variate normal (Gaussian)
    distribution, and calculate the fit statistic.

    Parameters
    ----------
    fit :
       The fit results.
    num : int, optional
       The number of samples to use (default is `1`).
    sigma : number, optional
       The width of the normal distribution (the default
       is `1`).
    correlate : bool, optional
       Should a multi-variate normal be used, with parameters
       set by the covariance matrix (`True`) or should a
       uni-variate normal be used (`False`)?
    numcores : optional
       The number of CPU cores to use. The default is to use all
       the cores on the machine.

    Returns
    -------
    samples :
       A NumPy array table with the first column representing the
       statistic and later columns the parameters used.

    See Also
    --------
    t_sample : Sample from the Student's t-distribution.
    uniform_sample : Sample from a uniform distribution.

    Notes
    -----
    All thawed model parameters are sampled from the Gaussian
    distribution, where the mean is set as the best-fit parameter
    value and the variance is determined by the diagonal elements
    of the covariance matrix. The multi-variate Gaussian is
    assumed by default for correlated parameters, using the
    off-diagonal elements of the covariance matrix.

    """
    sampler = NormalSampleFromScaleVector()
    sampler.scale.sigma = sigma
    if correlate:
        sampler = NormalSampleFromScaleMatrix()
    return sampler.get_sample(fit, num, numcores)


def uniform_sample(fit, num=1, factor=4, numcores=None):
    """Sample the fit statistic by taking the parameter values
    from an uniform distribution.

    For each iteration (sample), change the thawed parameters by
    drawing values from a uniform distribution, and calculate the
    fit statistic.

    Parameters
    ----------
    fit :
       The fit results.
    num : int, optional
       The number of samples to use (default is `1`).
    factor : number, optional
       Multiplier to expand the scale parameter (default is `4`).
    numcores : optional
       The number of CPU cores to use. The default is to use all
       the cores on the machine.

    Returns
    -------
    samples :
       A NumPy array table with the first column representing the
       statistic and later columns the parameters used.

    See Also
    --------
    normal_sample : Sample from a normal distribution.
    t_sample : Sample from the Student's t-distribution.

    """
    sampler = UniformSampleFromScaleVector()
    sampler.scale.sigma = 1
    return sampler.get_sample(fit, num, factor, numcores)


def t_sample(fit, num=1, dof=2, numcores=None):
    """Sample the fit statistic by taking the parameter values from
    a Student's t-distribution.

    For each iteration (sample), change the thawed parameters
    by drawing values from a Student's t-distribution, and
    calculate the fit statistic.

    Parameters
    ----------
    fit :
       The fit results.
    num : int, optional
       The number of samples to use (default is `1`).
    dof : optional
       The number of degrees of freedom to use (the default
       is to use the number from the current fit).
    numcores : optional
       The number of CPU cores to use. The default is to use all
       the cores on the machine.

    Returns
    -------
    samples :
       A NumPy array table with the first column representing the
       statistic and later columns the parameters used.

    See Also
    --------
    normal_sample : Sample from the normal distribution.
    uniform_sample : Sample from a uniform distribution.

    """
    sampler = StudentTSampleFromScaleMatrix()
    return sampler.get_sample(fit, num, dof, numcores)
