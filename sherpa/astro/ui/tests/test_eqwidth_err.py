#
#  Copyright (C) 2018, 2021, 2023
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
