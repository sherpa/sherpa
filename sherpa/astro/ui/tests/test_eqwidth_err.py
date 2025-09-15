#
#  Copyright (C) 2018, 2021, 2023, 2025
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

import numpy as np

import pytest

from sherpa.astro import ui
from sherpa.utils.err import IOErr, SessionErr
from sherpa.utils.testing import requires_data, requires_fits, requires_xspec


@requires_data
@requires_fits
@requires_xspec
def test_eqwith_err(make_data_path, restore_xspec_settings, clean_astro_ui):

    def check(a0, a1, a2):
        assert a0 == pytest.approx(0.16443033244310976, rel=1e-3)
        assert a1 == pytest.approx(0.09205564216156815, rel=1e-3)
        assert a2 == pytest.approx(0.23933118287470895, rel=1e-3)

    ui.set_method('neldermead')
    ui.set_stat('cstat')
    ui.set_xsabund('angr')
    ui.set_xsxsect('bcmc')

    ui.load_data(make_data_path('12845.pi'))
    ui.notice(0.5, 7)

    ui.set_model("xsphabs.gal*xszphabs.zabs*(powlaw1d.p1+xszgauss.g1)")
    ui.set_par(gal.nh, 0.08)
    ui.freeze(gal)

    ui.set_par(zabs.redshift, 0.518)
    ui.set_par(g1.redshift, 0.518)
    ui.set_par(g1.Sigma, 0.01)
    ui.freeze(g1.Sigma)
    ui.set_par(g1.LineE, min=6.0, max=7.0)

    ui.fit()

    ui.set_rng(np.random.RandomState(12345))
    result = ui.eqwidth(p1, p1 + g1, error=True, niter=100)
    check(result[0], result[1], result[2])
    params = result[3]

    ui.set_rng(np.random.RandomState(12345))
    result = ui.eqwidth(p1, p1 + g1, error=True, params=params, niter=100)
    check(result[0], result[1], result[2])

    parvals = ui.get_fit_results().parvals
    assert parvals[0] == pytest.approx(0.6111340686157877, rel=1.0e-3)
    assert parvals[1] == pytest.approx(1.6409785803466297, rel=1.0e-3)
    assert parvals[2] == pytest.approx(8.960926761312153e-05, rel=1.0e-3)
    assert parvals[3] == pytest.approx(6.620017726014523, rel=1.0e-3)
    assert parvals[4] == pytest.approx(1.9279114810359657e-06, rel=1.0e-3)


@requires_data
@requires_fits
@requires_xspec
def test_eqwith_err1(make_data_path, restore_xspec_settings, clean_astro_ui):

    def check1(e0, e1, e2):
        assert e0 == pytest.approx(0.028335201547206704, rel=1.0e-3)
        assert e1 == pytest.approx(-0.00744118799274448756, rel=1.0e-3)
        assert e2 == pytest.approx(0.0706249544851336, rel=1.0e-3)

    ui.set_xsabund('angr')
    ui.set_xsxsect('bcmc')

    ui.load_pha(make_data_path('3c273.pi'))
    ui.notice(0.5, 7.0)
    ui.set_stat("chi2datavar")
    ui.set_method("simplex")
    ui.set_model('powlaw1d.p1+gauss1d.g1')
    g1.fwhm = 0.1
    g1.pos = 2.0
    ui.freeze(g1.pos, g1.fwhm)
    ui.fit()

    ui.set_rng(np.random.RandomState(2345))
    e = ui.eqwidth(p1, p1 + g1, error=True, niter=100)
    check1(e[0], e[1], e[2])
    params = e[3]

    ui.set_rng(np.random.RandomState(2345))
    e = ui.eqwidth(p1, p1 + g1, error=True, params=params, niter=100)
    check1(e[0], e[1], e[2])

    parvals = ui.get_fit_results().parvals
    assert parvals[0] == pytest.approx(1.9055272902160334, rel=1.0e-3)
    assert parvals[1] == pytest.approx(0.00017387966749772638, rel=1.0e-3)
    assert parvals[2] == pytest.approx(1.279415076070516e-05, rel=1.0e-3)


def test_eqwidth_err_needs_fit(clean_astro_ui):
    """We get an error if fit has not been called"""

    ui.load_arrays(2, [1, 5], [1, 12], ui.Data1D)

    cmdl = ui.const1d.cmdl
    pmdl = ui.polynom1d.pmdl
    cmdl.c0 = 1.5
    pmdl.c0 = -5
    pmdl.c1 = 2
    ui.set_source(2, cmdl + pmdl)

    with pytest.raises(SessionErr) as exc:
        ui.eqwidth(cmdl, cmdl + pmdl, id=2, error=True)

    assert str(exc.value) == 'no fit has been performed'


@pytest.mark.parametrize('arg', ['params', 'covar_matrix'])
def test_eqwidth_err_arg_is_numpy(arg, clean_astro_ui):
    """Ensure argument is a NumPy argument."""

    ui.load_arrays('bob', [1, 5], [1, 12], ui.Data1D)

    cmdl = ui.const1d.cmdl
    pmdl = ui.polynom1d.pmdl
    cmdl.c0 = 1.5
    pmdl.c0 = -5
    pmdl.c1 = 2
    ui.set_source('bob', cmdl + pmdl)

    ui.fit('bob')

    arglist = {'id': 'bob', 'error': True, 'niter': 100, arg: 2.3}

    with pytest.raises(IOErr) as exc:
        ui.eqwidth(cmdl, cmdl + pmdl, **arglist)

    assert str(exc.value) == f'{arg} must be of type numpy.ndarray'


@pytest.mark.parametrize('arg', ['params', 'covar_matrix'])
def test_eqwidth_err_arg_is_2d(arg, clean_astro_ui):
    """Ensure argument is a 2D array."""

    ui.load_arrays('bob', [1, 5], [1, 12], ui.Data1D)

    cmdl = ui.const1d.cmdl
    pmdl = ui.polynom1d.pmdl
    cmdl.c0 = 1.5
    pmdl.c0 = -5
    pmdl.c1 = 2
    ui.set_source('bob', cmdl + pmdl)

    ui.fit('bob')

    arglist = {'id': 'bob', 'error': True, 'niter': 100,
               arg: np.arange(3)}

    with pytest.raises(IOErr) as exc:
        ui.eqwidth(cmdl, cmdl + pmdl, **arglist)

    assert str(exc.value) == f'{arg} must be 2d numpy.ndarray'


@pytest.mark.parametrize('arg', ['params', 'covar_matrix'])
def test_eqwidth_err_arg_size(arg, clean_astro_ui):
    """Ensure argument is the correct size."""

    ui.load_arrays('bob', [1, 5], [1, 12], ui.Data1D)

    cmdl = ui.const1d.cmdl
    pmdl = ui.polynom1d.pmdl
    cmdl.c0 = 1.5
    pmdl.c0 = -5
    pmdl.c1 = 2
    ui.set_source('bob', cmdl + pmdl)

    ui.fit('bob')

    arglist = {'id': 'bob', 'error': True, 'niter': 100,
               arg: np.arange(12).reshape(3, 4)}

    with pytest.raises(IOErr) as exc:
        ui.eqwidth(cmdl, cmdl + pmdl, **arglist)

    assert str(exc.value) == f'{arg} must be of dimension (2, x)'


def test_eqwidth_err_arg_square(clean_astro_ui):
    """Ensure argument is square.

    It looks like this is only checked for with covar_matrix argument
    """

    ui.load_arrays('bob', [1, 5], [1, 12], ui.Data1D)

    cmdl = ui.const1d.cmdl
    pmdl = ui.polynom1d.pmdl
    cmdl.c0 = 1.5
    pmdl.c0 = -5
    pmdl.c1 = 2
    ui.set_source('bob', cmdl + pmdl)

    ui.fit('bob')

    arglist = {'id': 'bob', 'error': True, 'niter': 100,
               'covar_matrix': np.arange(12).reshape(2, 6)}

    with pytest.raises(IOErr) as exc:
        ui.eqwidth(cmdl, cmdl + pmdl, **arglist)

    assert str(exc.value) == 'covar_matrix must be of dimension (2, 2)'


def setup_multi_id():
    """Setup the data for the multi-id case."""

    # We could have data with different responses but the simplest
    # case is to have the same response.
    #
    nchan = 200
    ui.dataspace1d(1, nchan, id="a", dstype=ui.DataPHA)
    ui.dataspace1d(1, nchan, id="b", dstype=ui.DataPHA)

    egrid = np.linspace(0.1, 2, nchan + 1)
    elo = egrid[:-1]
    ehi = egrid[1:]
    perfect_rmf = ui.create_rmf(elo, ehi)

    # Change the ARF just to get a different "response"
    arf_a = ui.create_arf(elo, ehi, np.full(nchan, 0.5))
    arf_b = ui.create_arf(elo, ehi, np.full(nchan, 2))

    ui.set_rmf("a", perfect_rmf)
    ui.set_rmf("b", perfect_rmf)

    ui.set_arf("a", arf_a)
    ui.set_arf("b", arf_b)

    # The "truth".
    #
    pl_t = ui.create_model_component("powlaw1d", "pl_t")
    gline_t = ui.create_model_component("gauss1d", "gline_t")

    # Ensure the counts are large enough for ~ gaussian statistics.
    # The line is roughly x=0.85 to 1.15.
    #
    pl_t.ampl = 5000
    pl_t.gamma = 1
    gline_t.pos = 1
    gline_t.fwhm = 0.1
    gline_t.ampl = 2000

    truth = pl_t - gline_t
    ui.set_source("a", truth)
    ui.set_source("b", truth)

    # Fake up the data
    #
    ui.set_rng(np.random.RandomState(382134))
    ui.fake_pha(id="a")
    ui.fake_pha(id="b")

    # To fit.
    pl = ui.create_model_component("powlaw1d", "pl")
    gline = ui.create_model_component("gauss1d", "gline")
    to_fit = pl - gline
    pl.ampl = 5000
    gline.pos = 1
    gline.pos.freeze()
    gline.fwhm = 0.1
    gline.fwhm.freeze()
    gline.ampl = 2000

    ui.set_source("a", to_fit)
    ui.set_source("b", to_fit)


def test_eqwidth_multi_id_chisq(clean_astro_ui):
    """Regression test for handling multiple ids (chi-square).

    This is only an issue when error=True. At the moment this is
    treated as a regression test as it is unclear what we want the
    results to be.

    """

    setup_multi_id()

    ui.set_stat("chi2datavar")
    ui.fit()
    fres = ui.get_fit_results()
    assert fres.datasets == ("a", "b")

    # Run eqwidth with no errors: the result should be the same
    # no matter the choice of id/otherids.
    #
    val1 = ui.eqwidth(pl, pl + gline, id="a", otherids=())
    val2 = ui.eqwidth(pl, pl + gline, id="b", otherids=())
    val3 = ui.eqwidth(pl, pl + gline, id="a", otherids=("b", ))
    val4 = ui.eqwidth(pl, pl + gline, id="b", otherids=("a", ))

    # These should be identical.
    assert val2 == val1
    assert val3 == val1
    assert val4 == val1

    # A regression to check if the calculation changes.
    #
    assert val1 == pytest.approx(0.040235587264827975)

    # Check the error results. We are using gaussian statistics
    # so this will not use get_draws but will sample the values
    # using the covariance matrix.
    #
    # Fix the random state before each run. Since the same model is
    # used for both datasets the random distributions should be the
    # same, so that the difference covariance matrices should be
    # apparent.
    #
    ui.set_rng(np.random.RandomState(287))
    resa = ui.eqwidth(pl, pl + gline, id="a", otherids=(), error=True, niter=100)
    ui.set_rng(np.random.RandomState(287))
    resb = ui.eqwidth(pl, pl + gline, id="b", otherids=(), error=True, niter=100)
    ui.set_rng(np.random.RandomState(287))
    resab = ui.eqwidth(pl, pl + gline, id="a", otherids=("b", ), error=True, niter=100)
    ui.set_rng(np.random.RandomState(287))
    resba = ui.eqwidth(pl, pl + gline, id="b", otherids=("a"), error=True, niter=100)

    # This is a "simple" check to see if the results are as
    # expected. For now all that we check against is the different
    # median values and assume that this means the correct data has
    # been used in creating the parameter values used for each
    # iteration.
    #
    # Ideally the distributions for the multi-id case would be
    # "tighter" than the single case, which is hard to check but can
    # be approximated by a "smaller" sigma limits range, and the
    # "a + b" and "b + a" would give the same results. However,
    # for now the otherids argument is ignored (issue #2379).
    #
    assert resab[4] == pytest.approx(resa[4])
    assert resba[4] == pytest.approx(resb[4])

    assert resa[0] == pytest.approx(0.040947827224568634)
    assert resb[0] == pytest.approx(0.04060370992021427)

    delta_a = resa[2] - resa[1]
    delta_ab = resab[2] - resab[1]
    assert delta_ab == delta_a

    delta_b = resb[2] - resb[1]
    delta_ba = resba[2] - resba[1]
    assert delta_ba == delta_b


def test_eqwidth_multi_id_poisson(clean_astro_ui):
    """Regression test for handling multiple ids (poisson).

    This is only an issue when error=True. At the moment this is
    treated as a regression test as it is unclear what we want the
    results to be.

    """

    setup_multi_id()

    ui.set_stat("cstat")
    ui.fit()
    fres = ui.get_fit_results()
    assert fres.datasets == ("a", "b")

    # Run eqwidth with no errors: the result should be the same
    # no matter the choice of id/otherids.
    #
    val1 = ui.eqwidth(pl, pl + gline, id="a", otherids=())
    val2 = ui.eqwidth(pl, pl + gline, id="b", otherids=())
    val3 = ui.eqwidth(pl, pl + gline, id="a", otherids=("b", ))
    val4 = ui.eqwidth(pl, pl + gline, id="b", otherids=("a", ))

    # These should be identical.
    assert val2 == val1
    assert val3 == val1
    assert val4 == val1

    # A regression to check if the calculation changes.
    #
    assert val1 == pytest.approx(0.0409343934694133)

    # Check the error results. This is calculated via get_draws,
    # which handles the multiple-id case.
    #
    # Fix the random state before each run, although we expect the
    # results to be different.
    #
    ui.set_rng(np.random.RandomState(287))
    resa = ui.eqwidth(pl, pl + gline, id="a", otherids=(), error=True, niter=100)
    ui.set_rng(np.random.RandomState(287))
    resb = ui.eqwidth(pl, pl + gline, id="b", otherids=(), error=True, niter=100)
    ui.set_rng(np.random.RandomState(287))
    resab = ui.eqwidth(pl, pl + gline, id="a", otherids=("b", ), error=True, niter=100)
    ui.set_rng(np.random.RandomState(287))
    resba = ui.eqwidth(pl, pl + gline, id="b", otherids=("a"), error=True, niter=100)

    # This is a "simple" check to see if the results are as
    # expected. For now all that we check against is the different
    # median values and assume that this means the correct data has
    # been used in creating the parameter values used for each
    # iteration.
    #
    # Ideally the distributions for the multi-id case would be
    # "tighter" than the single case, which is hard to check but can
    # be approximated by a "smaller" sigma limits range, and the
    # "a + b" and "b + a" would give the same results (the latter does
    # not hold at this time, so is it a bug in the code or our
    # understanding of the code?).
    #
    assert resa[0] == pytest.approx(0.04259715629705425)
    assert resb[0] == pytest.approx(0.043209134605777154)
    assert resab[0] == pytest.approx(0.04307866629966477)
    assert resba[0] == pytest.approx(0.04104380634099836)

    delta_a = resa[2] - resa[1]
    delta_ab = resab[2] - resab[1]
    assert delta_ab < delta_a

    delta_b = resb[2] - resb[1]
    delta_ba = resba[2] - resba[1]
    assert delta_ba < delta_b
