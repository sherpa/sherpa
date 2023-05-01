#
#  Copyright (C) 2011, 2016, 2018, 2020, 2021, 2023, 2024
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

from collections import namedtuple
from io import StringIO

import numpy as np

import pytest

from sherpa.data import Data1D
from sherpa.models import Gauss1D, PowLaw1D
from sherpa.fit import Fit
from sherpa.stats import Cash, Chi2DataVar, CStat
from sherpa.optmethods import NelderMead, LevMar
from sherpa.estmethods import Covariance
from sherpa import sim
from sherpa.utils.err import EstErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.parallel import _multi, _ncpus


_max = np.finfo(np.float32).max
_tiny = np.finfo(np.float32).tiny
_eps = np.finfo(np.float32).eps


_fit_results_bench = {
    'rstat': 89.29503933428586,
    'qval': 0.0,
    'succeeded': 1,
    'numpoints': 100,
    'dof': 95,
    'nfev': 93,
    'statval': 8483.0287367571564,
    'parnames': ['p1.gamma', 'p1.ampl', 'g1.fwhm',
                 'g1.pos', 'g1.ampl'],
    'parvals': np.array(
        [1.0701938169914813,
         9.1826254677279469,
         2.5862083052721028,
         2.601619746022207,
         47.262657692418749])
    }

_x = np.arange(0.1, 10.1, 0.1)
_y = np.array(
    [114, 47, 35, 30, 40, 27, 30, 26, 24, 20, 26, 35,
     29, 28, 34, 36, 43, 39, 33, 47, 44, 46, 53, 56,
     52, 53, 49, 57, 49, 36, 33, 42, 49, 45, 42, 32,
     31, 34, 18, 24, 25, 11, 17, 17, 11,  9,  8,  5,
     4, 10,  3,  4,  6,  3,  0,  2,  4,  4,  0,  1,
     2,  0,  3,  3,  0,  2,  1,  2,  3,  0,  1,  0,
     1,  0,  0,  1,  3,  3,  0,  2,  0,  0,  1,  2,
     0,  1,  0,  1,  1,  0,  1,  1,  1,  1,  1,  1,
     1,  0,  1,  0]
    )
_err = np.ones(100) * 0.4


@pytest.fixture
def setup():
    data = Data1D('fake', _x, _y, _err)

    g1 = Gauss1D('g1')
    g1.fwhm.set(1.0, _tiny, _max, frozen=False)
    g1.pos.set(1.0, -_max, _max, frozen=False)
    g1.ampl.set(1.0, -_max, _max, frozen=False)
    p1 = PowLaw1D('p1')
    p1.gamma.set(1.0, -10, 10, frozen=False)
    p1.ampl.set(1.0, 0.0, _max, frozen=False)
    p1.ref.set(1.0, -_max, _max, frozen=True)
    model = p1 + g1

    method = LevMar()
    method.config['maxfev'] = 10000
    method.config['ftol'] = float(_eps)
    method.config['epsfcn'] = float(_eps)
    method.config['gtol'] = float(_eps)
    method.config['xtol'] = float(_eps)
    method.config['factor'] = float(100)

    # The optimizer here does not use a RNG which means we do not have
    # to worry about setting it here (to match the test data which was
    # written using the Legacy NumPy random API).
    #
    fit = Fit(data, model, Chi2DataVar(), method, Covariance())
    results = fit.fit()

    for key in ["succeeded", "numpoints", "nfev"]:
        assert _fit_results_bench[key] == int(getattr(results, key))

    for key in ["rstat", "qval", "statval", "dof"]:
        # used rel and abs tol of 1e-7 with numpy allclose
        assert float(getattr(results, key)) == pytest.approx(_fit_results_bench[key])

    for key in ["parvals"]:
        # used rel and abs tol of 1e-4 with numpy allclose
        assert getattr(results, key) == pytest.approx(_fit_results_bench[key])

    fields = ['data', 'model', 'method', 'fit', 'results',
              'covresults', 'dof', 'mu', 'num', 'rng']
    out = namedtuple('Results', fields)

    out.data = data
    out.model = model
    out.method = method
    out.fit = fit
    out.results = results
    out.covresults = fit.est_errors()
    out.dof = results.dof
    out.mu = np.array(results.parvals)
    out.cov = np.array(out.covresults.extra_output)
    out.num = 10
    out.rng = np.random.RandomState(23)
    return out


# These are regression tests, and check that the code behaves the
# same.  They will need updating if the code or random generator
# changes.
#
# There are also some tests that check we return a subset of the
# information (e.g. all-but-the-first column).
#
EXPECTED_T = np.asarray(
    [[8.48554745e+03, 1.06405130, 9.27921234, 2.58664931, 2.59985055, 4.72032018e+01],
     [8.49039767e+03, 1.07407403, 9.05126945, 2.59679740, 2.60386838, 4.73836300e+01],
     [8.48487406e+03, 1.07593411, 9.09514664, 2.59122066, 2.60017483, 4.72314632e+01],
     [8.49161845e+03, 1.08176457, 8.99021849, 2.58955534, 2.60555871, 4.73149762e+01],
     [8.48500890e+03, 1.07345136, 9.10233367, 2.59110813, 2.59826504, 4.71895128e+01],
     [8.48768940e+03, 1.07931966, 9.04095922, 2.58360678, 2.59945504, 4.73628249e+01],
     [8.48965960e+03, 1.06198138, 9.36079509, 2.57625718, 2.59810196, 4.71039215e+01],
     [8.49112060e+03, 1.06906602, 9.21644535, 2.57147393, 2.60053414, 4.75622912e+01],
     [8.48414894e+03, 1.07225846, 9.13799921, 2.58982250, 2.59983239, 4.71818523e+01],
     [8.48467565e+03, 1.07387648, 9.12069098, 2.59440151, 2.60188040, 4.71716468e+01]])

EXPECTED_UNIFORM = np.asarray(
    [[8.92261355e+03, 1.07123512, 8.57093063, 2.61429088, 2.61114484, 4.73497877e+01],
     [9.82725325e+03, 1.09710009, 9.65477406, 2.59719719, 2.60793168, 4.75320785e+01],
     [9.26055735e+03, 1.08617397, 9.65589944, 2.55228604, 2.59058740, 4.74161153e+01],
     [8.89791148e+03, 1.05709447, 8.93723905, 2.54271093, 2.60441234, 4.73194803e+01],
     [8.57046575e+03, 1.05340130, 9.29276194, 2.62471517, 2.59089248, 4.69458349e+01],
     [9.30936360e+03, 1.08140402, 9.77082772, 2.55498773, 2.58937941, 4.70072160e+01],
     [8.50747533e+03, 1.05015625, 9.60690127, 2.57938039, 2.61121593, 4.72686196e+01],
     [9.07680066e+03, 1.06371906, 8.64790722, 2.57283955, 2.60058787, 4.71609855e+01],
     [8.50393051e+03, 1.07730034, 8.93027393, 2.61841172, 2.59182705, 4.72459725e+01],
     [8.65618143e+03, 1.06489218, 8.92190134, 2.57999046, 2.60300976, 4.77114146e+01]])

EXPECTED_NORMAL = np.asarray(
    [[8.48932906e+03, 1.07521274, 9.12922751, 2.58401110, 2.60674560, 4.72643098e+01],
     [8.51403222e+03, 1.07038805, 9.28561349, 2.59758807, 2.60327743, 4.73778582e+01],
     [8.49257783e+03, 1.06434242, 9.23215258, 2.59206666, 2.60425951, 4.73144299e+01],
     [8.51807534e+03, 1.07733205, 9.19189905, 2.59504891, 2.60874881, 4.73133218e+01],
     [8.48944871e+03, 1.07547372, 9.02226669, 2.58883295, 2.59963583, 4.72538788e+01],
     [8.60701037e+03, 1.06228469, 9.02741627, 2.57583942, 2.59915910, 4.72570650e+01],
     [8.49815729e+03, 1.06742811, 9.25051210, 2.58472492, 2.59254295, 4.73907749e+01],
     [8.50229908e+03, 1.06163471, 9.35611275, 2.60000339, 2.60215199, 4.72510835e+01],
     [8.75384863e+03, 1.06024498, 8.90014878, 2.58809846, 2.60381520, 4.72108582e+01],
     [8.51787451e+03, 1.08352963, 9.03835426, 2.57289213, 2.60153812, 4.71947163e+01]])

EXPECTED_NORMAL2 = np.asarray(
    [[8.48548565e+03, 1.06412678, 9.27802554, 2.58664389, 2.59987229, 4.72039324e+01],
     [8.49055139e+03, 1.07411404, 9.04991502, 2.59690659, 2.60389156, 4.73848773e+01],
     [8.48477491e+03, 1.07577639, 9.09755020, 2.59108294, 2.60021453, 4.72323203e+01],
     [8.48969939e+03, 1.08038114, 9.01322326, 2.58915516, 2.60508776, 4.73087209e+01],
     [8.48517760e+03, 1.07358732, 9.09898261, 2.59131263, 2.59812503, 4.71864601e+01],
     [8.48702560e+03, 1.07863842, 9.05153453, 2.58380099, 2.59961664, 4.73553475e+01],
     [8.48985527e+03, 1.06186041, 9.36341964, 2.57611060, 2.59805014, 4.71015833e+01],
     [8.49014958e+03, 1.06913589, 9.21435006, 2.57238678, 2.60060140, 4.75437276e+01],
     [8.48425884e+03, 1.07235749, 9.13585873, 2.58999586, 2.59974666, 4.71779765e+01],
     [8.48460130e+03, 1.07379213, 9.12210957, 2.59421385, 2.60187443, 4.71737314e+01]])


def test_student_t(setup):
    out = sim.multivariate_t(setup.mu, setup.cov, setup.dof, setup.num,
                             rng=setup.rng)

    assert out == pytest.approx(EXPECTED_T[:, 1:])


def test_cauchy(setup):
    out = sim.multivariate_cauchy(setup.mu, setup.cov, setup.num,
                                  rng=setup.rng)

    expected = [1.04452632, 9.58622941, 2.58805112, 2.59422687, 47.01421166]

    assert out == pytest.approx(np.asarray(expected))


def test_parameter_scale_vector_checks_scale_iterable(setup):
    """Error check"""

    with pytest.raises(TypeError,
                       match="^scales option must be iterable of length 5$"):
        sim.ParameterScaleVector().get_scales(setup.fit, myscales=3)


@pytest.mark.parametrize("scales", [3, np.asarray([3])])
def test_parameter_scale_matrix_checks_scale_iterable(scales, setup):
    """Error check"""

    with pytest.raises(EstErr,
                       match=r"^scales must be a numpy array of size \(5,5\)$"):
        sim.ParameterScaleMatrix().get_scales(setup.fit, myscales=scales)


def test_parameter_scale_matrix_checks_scale_positive_definite(setup):
    """Error check"""

    scales = np.arange(-10, 15).reshape(5, 5)
    with pytest.raises(TypeError,
                       match="^The covariance matrix is not positive definite$"):
        sim.ParameterScaleMatrix().get_scales(setup.fit, myscales=scales)


def test_parameter_scale_vector(setup):
    ps = sim.ParameterScaleVector()
    out = ps.get_scales(setup.fit)

    expected = [0.00752475, 0.15368132, 0.01088586, 0.00362169, 0.12308473]

    assert out == pytest.approx(np.asarray(expected))


def test_parameter_scale_matrix(setup):
    ps = sim.ParameterScaleMatrix()
    out = ps.get_scales(setup.fit)

    expected = [[ 5.66219290e-05, -1.13204057e-03,  5.74775798e-05, -1.41279846e-05,   3.89927955e-04],
                [-1.13204057e-03,  2.36179487e-02, -1.24335374e-03,  3.07465479e-04,  -7.64855819e-03],
                [ 5.74775798e-05, -1.24335374e-03,  1.18501872e-04, -1.78172521e-05,  -8.27934918e-05],
                [-1.41279846e-05,  3.07465479e-04, -1.78172521e-05,  1.31166529e-05,  -8.41924293e-05],
                [ 3.89927955e-04, -7.64855819e-03, -8.27934918e-05, -8.41924293e-05,   1.51498514e-02]]

    assert out == pytest.approx(np.asarray(expected))


def test_parameter_sample_checks_clip_argument(setup):
    """Error check"""

    obj = sim.NormalParameterSampleFromScaleVector()
    with pytest.raises(ValueError,
                       match="^invalid clip argument: sent max$"):
        obj.clip(setup.fit, np.arange(20).reshape(4, 5),
                 clip="max")


def test_uniform_parameter_sample(setup):
    ps = sim.UniformParameterSampleFromScaleVector()
    out = ps.get_sample(setup.fit, num=setup.num, rng=setup.rng)

    assert out == pytest.approx(EXPECTED_UNIFORM[:, 1:])


def test_normal_parameter_sample_vector(setup):
    ps = sim.NormalParameterSampleFromScaleVector()
    out = ps.get_sample(setup.fit, num=setup.num, rng=setup.rng)

    assert out == pytest.approx(EXPECTED_NORMAL[:, 1:])


def test_normal_parameter_sample_matrix(setup):
    ps = sim.NormalParameterSampleFromScaleMatrix()
    out = ps.get_sample(setup.fit, num=setup.num, rng=setup.rng)

    assert out == pytest.approx(EXPECTED_NORMAL2[:, 1:])


def test_t_parameter_sample_matrix(setup):
    ps = sim.StudentTParameterSampleFromScaleMatrix()
    out = ps.get_sample(setup.fit, dof=setup.dof, num=setup.num,
                        rng=setup.rng)

    assert out == pytest.approx(EXPECTED_T[:, 1:])


def test_uniform_sample(setup):
    ps = sim.UniformSampleFromScaleVector()
    out = ps.get_sample(setup.fit, num=setup.num, rng=setup.rng)

    assert out == pytest.approx(EXPECTED_UNIFORM)


# This used to have the same name as test_uniform_sample,
# so it has been renamed.
#
def test_uniform_sample2(setup):
    out = sim.uniform_sample(setup.fit, num=setup.num, rng=setup.rng)

    assert out == pytest.approx(EXPECTED_UNIFORM)


def test_normal_sample_vector(setup):
    ps = sim.NormalSampleFromScaleVector()
    out = ps.get_sample(setup.fit, num=setup.num, rng=setup.rng)

    assert out == pytest.approx(EXPECTED_NORMAL)


def test_normal_sample_matrix(setup):
    ps = sim.NormalSampleFromScaleMatrix()
    out = ps.get_sample(setup.fit, num=setup.num, rng=setup.rng)

    assert out == pytest.approx(EXPECTED_NORMAL2)


def test_t_sample_matrix(setup):
    ps = sim.StudentTSampleFromScaleMatrix()
    out = ps.get_sample(setup.fit, num=setup.num, dof=setup.dof,
                        rng=setup.rng)

    assert out == pytest.approx(EXPECTED_T)


def test_normal_sample(setup):
    out = sim.normal_sample(setup.fit, num=setup.num,
                            correlate=False, rng=setup.rng)

    assert out == pytest.approx(EXPECTED_NORMAL)


def test_normal_sample_correlated(setup):
    out = sim.normal_sample(setup.fit, num=setup.num,
                            correlate=True, rng=setup.rng)

    assert out == pytest.approx(EXPECTED_NORMAL2)


def test_t_sample(setup):
    out = sim.t_sample(setup.fit, setup.num, setup.dof, rng=setup.rng)

    assert out == pytest.approx(EXPECTED_T)


RATIOS_ONE = np.asarray([2.43733721, 3.52628288, 4.18396336, 1.73659659, 3.70862308, 2.19009678,
                         0.99262327, 2.47817863, 4.29140984, 8.66538795, 0.75845443, 2.19029536,
                         1.08725896, 1.58262032, 2.91657742, 0.70781436, 0.79954232, 2.5850919,
                         2.0543363,  0.4747516,  8.69094441, 2.35362447, 2.23331886, 4.20676696,
                         3.56214367])

RATIOS_TWO = np.asarray([1.74752217e+00, 4.15056013,  3.96667032e+00, 1.45119325e+00,
                         1.99113224e+00, 4.60875346,  2.61179034e+00, 2.52248040e-01,
                         2.59828000e+00, 3.45398280,  5.44372792e+00, 5.62797204e+00,
                         6.80978133e+00, 1.06134553, -1.42108547e-12, 2.07549322e+00,
                         1.14993770e+00, 5.88951206,  2.61912975e+00, 1.31001158e+00,
                         4.63278539e-01, 2.21026052,  3.12162630e+00, 1.83602835e+00,
                         6.28386219e+00])


def test_lrt(setup):
    """There is a limited check of the results.

    These are regression tests, so will need to be updated if the code
    or random generator changes.

    It looks like each iteration can generate different results
    depending on the platform (linux vs macOS) due to numerical
    differences (assumed) - e.g. the fwhm may be 0.03 in one and 0.06
    in the other, which can then lead to very different amplitudes,
    even though the resulting fit is not very different.  Then,
    because the values are not reset after each fit - see issue #1746
    - we quickly get rather different ratio values.  At least, that is
    the current hypothesis. So this test is limited to a small number
    of iterations, where this difference does not seem to be a
    problem. To make it easier to write (since the values change as
    niter is changed, which makes it a loop to try and identify a
    sensible number of iterations), we run for 25 iterations but only
    check the first 3 numbers, as after that the ratios begin to
    differ.

    """

    results = sim.LikelihoodRatioTest.run(setup.fit, setup.fit.model.lhs,
                                          setup.fit.model, niter=25, numcores=1,
                                          rng=setup.rng)

    assert results.null == pytest.approx(1376.9504116259966)
    assert results.alt == pytest.approx(98.57840260844087)
    assert results.lr == pytest.approx(1278.3720090175557)
    assert results.ppp == pytest.approx(0.0)

    assert results.samples.shape == (25, 2)
    assert results.stats.shape == (25, 2)

    # TODO: do we still need to restrict the elements being checked?
    assert results.ratios[:3] == pytest.approx(RATIOS_ONE[:3])


def test_lrt_multicore(setup):
    """The multi-core version of test_lrt.

    Testing this depends on how many cores were used, as that
    determines how many separate generators were used.
    """

    results = sim.LikelihoodRatioTest.run(setup.fit, setup.fit.model.lhs,
                                          setup.fit.model, niter=25, numcores=2,
                                          rng=setup.rng)

    assert results.null == pytest.approx(1376.9504116259966)
    assert results.alt == pytest.approx(98.57840260844087)
    assert results.lr == pytest.approx(1278.3720090175557)

    # This value is from the RNG so the check may need to be relaxed.
    assert results.ppp == pytest.approx(0.0)

    assert results.samples.shape == (25, 2)
    assert results.stats.shape == (25, 2)

    # The result depends on whether multi-processing is enabled and
    # if there's more than one core available, which we can check with
    # the _ncpus setting (although it appears that we still need to
    # check whether multi-processing is enabled or not).
    #
    if _multi and _ncpus > 1:
        expected = RATIOS_TWO
    else:
        expected = RATIOS_ONE

    # TODO: do we still need to restrict the elements being checked?
    assert results.ratios[:3] == pytest.approx(expected[:3])


def test_mh(setup, caplog):

    setup.fit.method = NelderMead()
    setup.fit.stat = Cash()

    # As part of investigating #2063 we now check the fit output
    out = StringIO()
    setup.fit.fit(outfile=out)

    results = setup.fit.est_errors()
    cov = results.extra_output

    # Note when the covariance changes; this is more just as a check
    # hence not a full check.
    #
    diag = np.asarray([1.68857268e-03, 6.48507899e-01,
                       9.22029508e-03, 1.88340111e-03,
                       3.22767792e+00])
    assert np.diag(cov) == pytest.approx(diag)

    # Check outfile; expect
    #     nfev statistic pl.gamma pl.ampl g1.fwhm g1.pos g1.ampl
    #
    row0 = np.asarray([0, -8975.692, 1.070194, 9.182625, 2.586208, 2.601620, 47.26266])
    row19 = np.asarray([19, -8975.953, 1.069763, 9.219122,  2.570933, 2.586584,  47.36504])
    row20 = row19.copy()
    row20[0] = 20

    out.seek(0)
    fitvals = np.loadtxt(out)
    assert fitvals.shape == (21, 7)
    assert fitvals[0] == pytest.approx(row0)
    assert fitvals[19] == pytest.approx(row19)
    assert fitvals[20] == pytest.approx(row20)

    mcmc = sim.MCMC()
    for par in setup.fit.model.pars:
        mcmc.set_prior(par, sim.flat)
        prior = mcmc.get_prior(par)
        assert prior.__name__ == 'flat'

    mcmc.set_sampler('MH')

    opt = mcmc.get_sampler_opt('defaultprior')
    mcmc.set_sampler_opt('defaultprior', opt)
    # mcmc.set_sampler_opt('verbose', True)

    with SherpaVerbosity("ERROR"):
        stats, accept, params = mcmc.get_draws(setup.fit, cov,
                                               niter=1e2, rng=setup.rng)

    assert len(caplog.records) == 0

    assert len(stats) == 101
    assert len(accept) == 101
    assert params.shape == (5, 101)

    assert stats.mean() == pytest.approx(-8970.945533198466)
    assert accept.sum() == 81
    assert accept.min() == 0
    assert accept.max() == 1

    means = np.asarray([1.06497461, 9.2188531, 2.58191838, 2.57882817, 47.5158447])
    assert params.mean(axis=1) == pytest.approx(means)


def test_metropolisMH(setup, caplog):

    setup.fit.method = NelderMead()
    setup.fit.stat = CStat()

    # As part of investigating #2063 we now check the fit output
    out = StringIO()
    setup.fit.fit(outfile=out)

    results = setup.fit.est_errors()
    cov = results.extra_output

    # Note when the covariance changes; this is more just as a check
    # hence not a full check.
    #
    diag = np.asarray([1.68857268e-03, 6.48507899e-01,
                       9.22029508e-03, 1.88340111e-03,
                       3.22767792e+00])
    assert np.diag(cov) == pytest.approx(diag)

    # Check outfile; expect
    #     nfev statistic pl.gamma pl.ampl g1.fwhm g1.pos g1.ampl
    #
    row0 = np.asarray([0, 98.83959, 1.070194, 9.182625, 2.586208, 2.601620, 47.26266])
    row19 = np.asarray([19, 98.5784, 1.069763, 9.219122, 2.570933, 2.586584, 47.36504])
    row20 = row19.copy()
    row20[0] = 20

    out.seek(0)
    fitvals = np.loadtxt(out)
    assert fitvals.shape == (21, 7)
    assert fitvals[0] == pytest.approx(row0)
    assert fitvals[19] == pytest.approx(row19)
    assert fitvals[20] == pytest.approx(row20)

    mcmc = sim.MCMC()
    mcmc.set_sampler('MetropolisMH')
    # mcmc.set_sampler_opt('verbose', True)

    with SherpaVerbosity("ERROR"):
        stats, accept, params = mcmc.get_draws(setup.fit, cov,
                                               niter=1e2, rng=setup.rng)

    assert len(caplog.records) == 0

    assert len(stats) == 101
    assert len(accept) == 101
    assert params.shape == (5, 101)

    assert stats.mean() == pytest.approx(103.81260679391825)
    assert accept.sum() == 50
    assert accept.min() == 0
    assert accept.max() == 1

    means = np.asarray([1.06278015, 9.21533855, 2.5736483, 2.5853907, 47.27058904])
    assert params.mean(axis=1) == pytest.approx(means)
