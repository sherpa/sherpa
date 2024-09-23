#
#  Copyright (C) 2007, 2015 - 2021, 2023, 2024
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

import warnings

from numpy import sqrt

from pytest import approx
import pytest

from sherpa.astro import ui
from sherpa.astro.data import DataPHA

from sherpa.utils.parallel import ncpus
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_xspec, requires_group, requires_region, requires_wcs

try:
    from sherpa.astro import xspec
    has_xspec = True
except ImportError:
    has_xspec = False


try:
    from sherpa.astro.io import backend
    is_crates_io = backend.name == "crates"
except ImportError:
    is_crates_io = False


# This value should match the [ogip] minimum_energy setting in the
# Sherpa configuration file.
#
EMIN = 1.0e-10


@pytest.fixture
def fix_xspec(clean_astro_ui, hide_logging):
    """As of XSPEC 12.11.1 it's useful to fix abundance/cross sections
    rather than rely on reading from the user's directory

    This requires XSPEC support.
    """

    # As of XSPEC 12.10.1, it is safest to explicitly set
    # the state to a known value. For now just pick the
    # "easy" settings.
    #
    abund = xspec.get_xsabund()
    xsect = xspec.get_xsxsect()

    xspec.set_xsabund('angr')
    xspec.set_xsxsect('bcmc')

    yield

    xspec.set_xsabund(abund)
    xspec.set_xsxsect(xsect)


def check_thread(run_thread, thread, parallel, cmpfunc, parnames):
    """Run the thread test and call the comparison function"""

    if parallel:
        thread += '_ncpus'

    tlocals = run_thread(thread)
    fit_results = ui.get_fit_results()
    covarerr = sqrt(fit_results.extra_output['covar'].diagonal())

    # what parameters do we want to test
    #
    parvals = [tlocals[p] for p in parnames]

    # compare the results
    cmpfunc(fit_results, *parvals, covarerr)

    if not parallel or ncpus < 2:
        assert fit_results.extra_output['num_parallel_map'] == 0
    else:
        assert fit_results.extra_output['num_parallel_map'] > 0


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize('parallel', [False, True])
def test_pha_intro(parallel, run_thread, fix_xspec):

    # Currently there's no test of covarerr
    def cmp_thread(fit_result, p1, covarerr):
        assert fit_result.statval == approx(37.9079, rel=1e-4)
        assert fit_result.rstat == approx(0.902569, rel=1e-4)
        assert fit_result.qval == approx(0.651155, rel=1e-4)
        assert fit_result.nfev == 22
        assert fit_result.numpoints == 44
        assert fit_result.dof == 42
        p1.gamma.val == approx(2.15852, rel=1e-4)
        p1.ampl.val == approx(0.00022484, rel=1e-4)

        assert ui.calc_photon_flux() == approx(0.000469964, rel=1e-4)
        assert ui.calc_energy_flux() == approx(9.614847e-13, rel=1e-4)
        assert ui.calc_data_sum() == approx(706.85714092, rel=1e-4)
        assert ui.calc_model_sum() == approx(638.45693377, rel=1e-4)
        assert ui.calc_source_sum() == approx(0.046996409, rel=1e-4)

        calc = ui.eqwidth(p1, ui.get_source())
        assert calc == approx(-0.57731725, rel=1e-4)

        calc = ui.calc_kcorr([1, 1.2, 1.4, 1.6, 1.8, 2], 0.5, 2)
        expected = [0.93132747, 0.9352768, 0.94085917,
                    0.94738472, 0.95415463, 0.96121113]
        assert calc == approx(expected, rel=1e-4)

    check_thread(run_thread, 'pha_intro', parallel, cmp_thread, ['p1'])


@requires_data
@requires_fits
def test_pha_read(run_thread):
    run_thread('pha_read')
    assert type(ui.get_data()) == DataPHA


@requires_data
@requires_fits
@pytest.mark.parametrize('parallel', [False, True])
def test_basic(parallel, run_thread, clean_astro_ui):

    def cmp_thread(fit_results, m1, m2, covarerr):
        assert fit_results.nfev == 9
        assert fit_results.numpoints == 11
        assert fit_results.dof == 9

        assert covarerr[0] == approx(0.0192539, rel=1e-4)
        assert covarerr[1] == approx(0.00392255, rel=1e-4)
        assert fit_results.statval == approx(151.827, rel=1e-4)
        assert fit_results.rstat == approx(16.8697, rel=1e-4)
        assert fit_results.qval == approx(3.68798e-28, rel=1e-4)
        assert m1.c0.val == approx(1.49843, rel=1e-4)
        assert m1.c1.val == approx(0.1447, rel=1e-4)
        assert m1.c2.val == approx(0.0322936, rel=1e-4)
        assert m1.c3.val == approx(-0.00277729, rel=1e-4)
        assert m2.c0.val == approx(1.75548, rel=1e-4)
        assert m2.c1.val == approx(0.198455, rel=1e-4)

    # In data1.dat for this test, there is a comment with one
    # word at the beginning -- deliberately would break when reading
    # with DM ASCII kernel, but passes because we have Sherpa code
    # to bypass that.
    #
    check_thread(run_thread, 'basic', parallel, cmp_thread, ['m1', 'm2'])


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize('parallel', [False, True])
def test_simultaneous(parallel, run_thread, fix_xspec):

    def cmp_thread(fit_results, abs1, pl1, pl2, covarerr):
        assert fit_results.numpoints == 18
        assert fit_results.dof == 14

        assert covarerr[0] == approx(0.397769, rel=1e-3)
        assert covarerr[1] == approx(0.486058, rel=1e-3)
        assert covarerr[2] == approx(1.48213e-05, rel=1e-3)
        assert covarerr[3] == approx(1.54245e-05, rel=1e-3)
        assert fit_results.statval == approx(7.4429, rel=1e-4)
        assert fit_results.rstat == approx(0.531636, rel=1e-4)
        assert fit_results.qval == approx(0.916288, rel=1e-4)
        assert abs1.nh.val == approx(0.898162, rel=1e-2)
        assert pl1.gamma.val == approx(1.645, rel=1e-4)
        assert pl1.ampl.val == approx(2.28323e-05, 1e-3)
        assert pl2.ampl.val == approx(2.44585e-05, 1e-3)

    check_thread(run_thread, 'simultaneous', parallel, cmp_thread,
                 ['abs1', 'pl1', 'pl2'])


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize('parallel', [False, True])
def test_sourceandbg(parallel, run_thread, fix_xspec):

    def cmp_thread(fit_results, a1, b1, b2, covarerr):
        assert fit_results.numpoints == 1330
        assert fit_results.dof == 1325

        assert covarerr[0] == approx(0.012097, rel=1.07e-3)
        assert covarerr[1] == approx(0, rel=1e-3)
        assert covarerr[2] == approx(0.000280678, rel=1e-3)
        assert covarerr[3] == approx(0.00990783, rel=1e-3)
        assert covarerr[4] == approx(2.25746e-07, rel=1e-3)
        assert fit_results.statval == approx(947.5, rel=1e-4)
        assert fit_results.rstat == approx(0.715094, rel=1e-4)
        assert fit_results.qval == approx(1, rel=1e-4)
        assert a1.nh.val == approx(0.0342266, rel=1e-2)
        assert b1.kt.val == approx(20, rel=1e-2)
        assert b1.norm.val == approx(0.00953809, 1e-2)
        assert b2.kt.val == approx(0.563109, 1e-2)
        assert b2.norm.val == approx(1.16118e-05, 1e-2)

    check_thread(run_thread, 'sourceandbg', parallel, cmp_thread,
                 ['a1', 'b1', 'b2'])


@requires_data
@requires_fits
@requires_region
@requires_wcs
def test_spatial(run_thread, clean_astro_ui):
    tlocals = run_thread('spatial')
    g1 = tlocals['g1']
    g2 = tlocals['g2']

    results = ui.get_fit_results()
    assert results.nfev == 397
    assert results.numpoints == 4881
    assert results.dof == 4877

    assert results.statval == approx(-59229.749441, 1e-4)
    assert g1.fwhm.val == approx(61.5615, rel=1e-4)
    assert g1.xpos.val == approx(4070.45, rel=1e-4)
    assert g1.ypos.val == approx(4251.35, rel=1e-4)
    assert g1.ampl.val == approx(22.1269, rel=1e-4)
    assert g2.fwhm.val == approx(6.20409, rel=1e-4)
    assert g2.xpos.val == approx(4070.78, rel=1e-4)
    assert g2.ypos.val == approx(4249.33, rel=1e-4)
    assert g2.ampl.val == approx(226.563, rel=1e-4)


def cmp_radpro(fit_results, src, covarerr):
    """Check the radpro/radpro_dm tests"""

    assert fit_results.nfev == 92
    assert fit_results.numpoints == 38
    assert fit_results.dof == 35

    assert covarerr[0] == approx(9.37345, rel=1e-4)
    assert covarerr[1] == approx(0.512596, rel=1e-4)
    assert covarerr[2] == approx(0.0691102, rel=1e-4)
    assert fit_results.statval == approx(217.450, rel=1e-4)
    assert fit_results.rstat == approx(6.21287, rel=1e-4)
    assert fit_results.qval == approx(0.0, rel=1e-4)
    assert src.r0.val == approx(125.829, rel=1e-4)
    assert src.beta.val == approx(4.1633, rel=1e-4)
    assert src.xpos.val == approx(0.0, rel=1e-4)
    assert src.ampl.val == approx(4.42821, rel=1e-4)


@requires_data
@requires_fits
@pytest.mark.parametrize('parallel', [False, True])
def test_radpro(parallel, run_thread, clean_astro_ui):

    check_thread(run_thread, 'radpro', parallel, cmp_radpro, ['src'])


# This is test_radpro but it uses crates syntax to load the
# data, hence we have a very-slightly-different thread. We have
# not bothered to add a "_ncpus" version of the thread to test
# out the parallel handling.
#
@requires_data
@pytest.mark.skipif(not is_crates_io, reason='Test requires crates')
def test_radpro_dm(run_thread, clean_astro_ui):

    check_thread(run_thread, 'radpro_dm', False, cmp_radpro, ['src'])


@requires_data
@requires_fits
@requires_region
def test_psf2d(run_thread, clean_astro_ui):
    tlocals = run_thread('psf')
    g1 = tlocals['g1']

    results = ui.get_fit_results()
    assert results.nfev == 296
    assert results.numpoints == 4899
    assert results.dof == 4895

    assert results.statval == approx(4066.78, rel=1e-4)
    assert g1.fwhm.val == approx(2.80117, rel=1e-2)
    assert g1.ypos.val == approx(77.2271, rel=1e-2)
    assert g1.xpos.val == approx(88.661, rel=1e-2)
    assert g1.ampl.val == approx(166.649, rel=1e-2)


@requires_data
@requires_fits
@requires_region
def test_fpsf2d(run_thread, clean_astro_ui):
    tlocals = run_thread('fpsf')

    fres = ui.get_fit_results()
    assert fres.nfev == 581
    assert fres.numpoints == 4899
    assert fres.dof == 4895

    assert fres.statval == approx(-4053.6635, rel=1e-4)

    # assert tlocals['b1'].xlow.val == approx(-4.70832, rel=1e-4)
    # assert tlocals['b1'].xhi.val == approx(164.687, rel=1e-4)
    # assert tlocals['b1'].ylow.val == approx(0.83626, rel=1e-4)
    # assert tlocals['b1'].yhi.val == approx(142.603, rel=1e-4)
    # assert tlocals['b1'].ampl.val == approx(0.956766, rel=1e-4)

    g1mdl = tlocals['g1']
    assert g1mdl.fwhm.val == approx(6.420237, rel=1e-4)
    assert g1mdl.xpos.val == approx(88.940712, rel=1e-4)
    assert g1mdl.ypos.val == approx(76.577265, rel=1e-4)
    assert g1mdl.ampl.val == approx(36344.48324, rel=1e-4)


@requires_data
@requires_fits
def test_radpro_psf(run_thread, clean_astro_ui):
    tlocals = run_thread('radpro_psf')

    results = ui.get_fit_results()
    assert results.nfev == 48
    assert results.numpoints == 38
    assert results.dof == 35
    assert results.statval == approx(200.949, rel=1e-4)

    src = tlocals['src']
    assert src.r0.val == approx(83.0997, rel=1e-4)
    assert src.beta.val == approx(2.97737, rel=1e-4)
    assert src.ampl.val == approx(5.27604, rel=1e-4)


@requires_data
@requires_fits
@pytest.mark.parametrize('parallel', [False, True])
def test_linepro(parallel, run_thread, clean_astro_ui):

    def cmp_thread(fit_results, b1, covarerr):
        assert fit_results.nfev == 17
        assert fit_results.numpoints == 75
        assert fit_results.dof == 72

        assert covarerr[0] == approx(0.176282, rel=1e-4)
        assert covarerr[1] == approx(0.0019578, rel=1e-4)
        assert covarerr[2] == approx(0.495889, rel=1e-4)
        assert fit_results.statval == approx(203.34, rel=1e-4)
        assert b1.r0.val == approx(4.25557, rel=1e-4)
        assert b1.beta.val == approx(0.492232, rel=1e-4)
        assert b1.ampl.val == approx(11.8129, rel=1e-4)

    check_thread(run_thread, 'linepro', parallel, cmp_thread, ['b1'])


@requires_data
@requires_fits
@pytest.mark.parametrize('parallel', [False, True])
def test_kernel(parallel, run_thread, clean_astro_ui):

    def cmp_thread(fit_results, b1, covarerr):
        assert fit_results.nfev == 21
        assert fit_results.numpoints == 75
        assert fit_results.dof == 72

        assert covarerr[0] == approx(0.210895, rel=1e-4)
        assert covarerr[1] == approx(0.00154839, rel=1e-4)
        assert covarerr[2] == approx(0.0223859, rel=1e-4)
        assert fit_results.statval == approx(98.5793, rel=1e-4)
        assert b1.r0.val == approx(19.2278, rel=1e-4)
        assert b1.beta.val == approx(0.555464, rel=1e-4)
        assert b1.ampl.val == approx(1.93706, rel=1e-4)

    check_thread(run_thread, 'kernel', parallel, cmp_thread, ['b1'])


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize('parallel', [False, True])
def test_spectrum(parallel, run_thread, fix_xspec):

    def cmp_thread(fres, abs2, mek1, mek2, covarerr):
        assert fres.numpoints == 446
        assert fres.dof == 441

        assert covarerr[0] == approx(0.00148391, rel=1e-3)
        assert covarerr[1] == approx(0.0011518, rel=1e-3)
        assert covarerr[2] == approx(0.00377755, rel=1e-3)
        assert covarerr[3] == approx(0.00370543, rel=1e-3)
        assert covarerr[4] == approx(0.0016608, rel=1e-3)
        assert fres.statval == approx(0.0496819, rel=1e-4)
        assert abs2.nh.val == approx(1.1015, rel=1e-4)
        assert mek1.kt.val == approx(0.841025, rel=1e-4)
        assert mek1.norm.val == approx(0.699761, rel=1e-4)
        assert mek2.kt.val == approx(2.35845, rel=1e-4)
        assert mek2.norm.val == approx(1.03724, rel=1e-4)

    check_thread(run_thread, 'spectrum', parallel, cmp_thread,
                 ['abs2', 'mek1', 'mek2'])


@requires_data
@requires_fits
def test_histo(run_thread, clean_astro_ui):
    tlocals = run_thread('histo')
    g1 = tlocals['g1']

    results = ui.get_fit_results()
    assert results.nfev == 3062
    assert results.numpoints == 50
    assert results.dof == 47
    assert results.statval == approx(14.7264, rel=1e-4)
    assert g1.fwhm.val == approx(0.0232473, rel=1e-4)
    assert g1.pos.val == approx(1.26713, rel=1e-4)
    assert g1.ampl.val == approx(40.4503, rel=1e-4)


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize('parallel', [False, True])
def test_xmm(parallel, run_thread, fix_xspec):

    def cmp_thread(fres, intrin, phard, covarerr):
        assert fres.nfev == 95
        assert fres.numpoints == 161
        assert fres.dof == 158

        assert covarerr[0] == approx(0.9582243108294659, rel=1e-3)
        assert covarerr[1] == approx(0.1429841749646365, rel=1e-3)
        assert covarerr[2] == approx(0.00039469440047366997, rel=1e-3)
        assert fres.statval == approx(117.59572388261346, rel=1e-4)

        assert intrin.nh.val == approx(11.1155, 1e-2)
        assert phard.phoindex.val == approx(1.498, 1e-2)
        assert phard.norm.val == approx(0.0014224, 1e-2)

    check_thread(run_thread, 'xmm', parallel, cmp_thread,
                 ['intrin', 'phard'])


@requires_data
@requires_fits
@pytest.mark.parametrize('parallel', [False, True])
def test_grouped_ciao4_5(parallel, run_thread, clean_astro_ui):

    def cmp_thread(fres, aa, covarerr):
        assert fres.numpoints == 46

        assert fres.statval == approx(19.0311, rel=1e-4)

        assert covarerr[0] == approx(0.104215, rel=1e-4)
        assert covarerr[1] == approx(2.52228e-05, rel=1e-4)
        assert aa.gamma.val == approx(1.83531, rel=1e-4)
        assert aa.ampl.val == approx(0.000302545, rel=1e-4)

    check_thread(run_thread, 'grouped_ciao4.5', parallel, cmp_thread,
                 ['aa'])


@requires_data
@requires_fits
@requires_xspec
def test_proj(run_thread, fix_xspec):
    tlocals = run_thread('proj')
    g = tlocals['g']
    xs1 = tlocals['xs1']
    p = tlocals['p']

    fit_res1 = tlocals['fit_res1']
    covar_res1 = tlocals['covar_res1']
    proj_res1 = tlocals['proj_res1']

    fit_res2 = tlocals['fit_res2']
    covar_res2 = tlocals['covar_res2']
    proj_res2 = tlocals['proj_res2']

    # Fit 1
    assert fit_res1.statval == approx(64.3803, 1e-2)
    assert g.fwhm.val == approx(32.3536, rel=1e-2)
    assert g.pos.val == approx(807.863, rel=1e-2)
    assert g.ampl.val == approx(117.826, rel=1e-2)

    # Covar 1
    assert covar_res1.parmins[0] == approx(-0.49357, 1e-2)
    assert covar_res1.parmins[1] == approx(-0.264056, 1e-2)
    assert covar_res1.parmins[2] == approx(-2.58857, 1e-2)
    assert covar_res1.parmaxes[0] == approx(0.49357, 1e-2)
    assert covar_res1.parmaxes[1] == approx(0.264056, 1e-2)
    assert covar_res1.parmaxes[2] == approx(2.58857, 1e-2)

    # Projection 1
    assert proj_res1.parmins[0] == approx(-0.492815, 1e-2)
    assert proj_res1.parmins[1] == approx(-0.263945, 1e-2)
    assert proj_res1.parmins[2] == approx(-2.57957, 1e-2)
    assert proj_res1.parmaxes[0] == approx(0.49421, 1e-2)
    assert proj_res1.parmaxes[1] == approx(0.26389, 1e-2)
    assert proj_res1.parmaxes[2] == approx(2.59695, 1e-2)

    # Fit 2
    assert fit_res2.statval == approx(108.048, 1e-2)
    assert xs1.nh.val == approx(0.0532586, rel=1e-2)
    assert p.gamma.val == approx(1.4816, rel=1e-2)
    assert p.ampl.val == approx(0.302343, rel=1e-2)

    # Covar 2
    assert covar_res2.parmins[0] == approx(-0.0379975, 5e-1)
    assert covar_res2.parmins[1] == approx(-0.0877478, 5e-1)
    assert covar_res2.parmins[2] == approx(-0.0831945, 5e-1)
    assert covar_res2.parmaxes[0] == approx(0.0379975, 5e-1)
    assert covar_res2.parmaxes[1] == approx(0.0877478, 5e-1)
    assert covar_res2.parmaxes[2] == approx(0.0831945, 5e-1)

    # Projection 2
    assert proj_res2.parmins[0] == approx(-0.0385636, 1e-2)
    assert proj_res2.parmins[1] == approx(-0.0891261, 1e-2)
    assert proj_res2.parmins[2] == approx(-0.0737413, 1e-2)
    assert proj_res2.parmaxes[0] == approx(0.0388651, 1e-2)
    assert proj_res2.parmaxes[1] == approx(0.0896566, 1e-2)
    assert proj_res2.parmaxes[2] == approx(0.0981627, 1e-2)


@requires_data
@requires_fits
@requires_xspec
def test_counts(run_thread, fix_xspec):
    tlocals = run_thread('counts')

    assert tlocals['counts_data1'] == approx(52701.0, rel=1e-4)
    assert tlocals['counts_data2'] == approx(25032.0, rel=1e-4)
    assert tlocals['counts_model1'] == approx(73226263.55355, 1e-4)
    assert tlocals['eflux1'] == approx(1.495482e-08, 1e-4)
    assert tlocals['counts_model2'] == approx(46082543.21529, 1e-4)
    assert tlocals['eflux2'] == approx(1.39662954483e-08, 1e-3)
    assert tlocals['pflux1'] == approx(1.6178938637, 1e-2)


@requires_data
@requires_fits
@requires_xspec
def test_stats_all(run_thread, fix_xspec):
    tlocals = run_thread('stats_all')

    assert tlocals['stat_lsqr'] == approx(213746.236464, 1e-4)
    assert tlocals['stat_chi2m'] == approx(1232.0330242, 1e-4)
    assert tlocals['stat_cash'] == approx(-1163566.90688, 1e-4)
    assert tlocals['stat_chi2c'] == approx(1411.60090961, 1e-4)
    assert tlocals['stat_chi2g'] == approx(972.388468358, 1e-4)
    assert tlocals['stat_chi2d'] == approx(1204.69363458, 1e-4)
    assert tlocals['stat_chi2x'] == approx(1204.69363458, 1e-4)
    assert tlocals['stat_cstat'] == approx(1210.56896183, 1e-4)


@requires_data
@requires_fits
@requires_region
@requires_wcs
def test_lev3fft(run_thread, clean_astro_ui):
    tlocals = run_thread('lev3fft', scriptname='bar.py')

    assert tlocals['src'].fwhm.val == approx(1.48914, rel=1e-5)
    assert tlocals['src'].xpos.val == approx(3142.84, rel=1e-5)
    assert tlocals['src'].ypos.val == approx(4519.72, rel=1e-5)
    assert tlocals['src'].ampl.val == approx(9.11268, rel=1e-5)
    assert tlocals['bkg'].c0.val == approx(0.0118703, rel=1e-5)

    fres = ui.get_fit_results()
    assert fres.istatval == approx(6483.84, rel=1e-5)
    assert fres.statval == approx(558.987, rel=1e-4)
    assert fres.numpoints == 3307
    assert fres.dof == 3302


@requires_data
@requires_fits
@requires_xspec
@requires_group
def test_setfullmodel(run_thread, fix_xspec):
    # TODO: add some tests
    run_thread('setfullmodel')


@requires_data
@requires_fits
def test_bug13537(run_thread, clean_astro_ui):
    run_thread('bug13537')


@requires_data
@requires_fits
@requires_xspec
def test_xmm2(run_thread, fix_xspec):

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        run_thread('xmm2')

    # NOTE: if this test is run on its own it can generate three warnings,
    # with the first being a RuntimeWarnnig about numpy.ndarray size
    # changed. We filter this out as it's not at all clear what is going
    # on and we have filters in conftest to remove similar warnings
    #
    ws = [w for w in ws if not (w.category == RuntimeWarning and
                                str(w.message).startswith('numpy.ndarray size changed, may indicate binary incompatibility.'))]
    assert len(ws) == 2
    cats = set([w.category for w in ws])
    assert cats == set([UserWarning])

    # The order of reading the ARF and RMF is not guaranteed,
    # so do not force it here when testing the two warning
    # messages.
    #
    arffile = 'MNLup_2138_0670580101_EMOS1_S001_spec.arf'
    rmffile = 'MNLup_2138_0670580101_EMOS1_S001_spec.rmf'
    emsg_arf = "The minimum ENERG_LO in the ARF " + \
               "'{}' ".format(arffile) + \
               "was 0 and has been replaced by {}".format(EMIN)
    emsg_rmf = "The minimum ENERG_LO in the RMF " + \
               "'{}' ".format(rmffile) + \
               "was 0 and has been replaced by {}".format(EMIN)

    emsgs = set([emsg_arf, emsg_rmf])
    wmsgs = set([str(w.message) for w in ws])
    assert wmsgs == emsgs

    assert ui.get_data().channel[0] == approx(0.0)
    rmf = ui.get_rmf()
    arf = ui.get_arf()
    assert rmf.detchans == 800
    assert len(rmf.energ_lo) == 2400
    assert len(rmf.energ_hi) == 2400
    assert len(rmf.n_grp) == 2400
    assert len(rmf.f_chan) == 2394
    assert len(rmf.n_chan) == 2394
    assert len(rmf.matrix) == 1281216
    assert rmf.offset == 0
    assert len(rmf.e_min) == 800
    assert len(rmf.e_max) == 800
    assert len(arf.energ_lo) == 2400
    assert len(arf.energ_hi) == 2400
    assert len(arf.specresp) == 2400

    etol = EMIN / 100.0
    assert rmf.energ_lo[0] == approx(EMIN, rel=etol)
    assert arf.energ_lo[0] == approx(EMIN, rel=etol)


@requires_data
@requires_fits
@requires_xspec
def test_missmatch_arf(make_data_path):
    ui.load_pha(1, make_data_path("source1.pi"))
    ui.load_bkg(1, make_data_path("back1.pi"))
    ui.load_arf(1, make_data_path("arf_1024.fits"))
    ui.load_rmf(1, make_data_path("rmf_1024.fits"))
    ui.set_method('levmar')
    ui.set_model(ui.powlaw1d.p1 * ui.xswabs.abs1)
    ui.set_par('p1.ampl', 0.0001)
    ui.set_stat('cash')
    ui.fit()
    parvals = ui.get_fit_results().parvals
    assert parvals[0] == approx(1.47969, rel=1.0e-3)
    assert parvals[1] == approx(0.0019491, rel=1.0e-3)
    assert parvals[2] == approx(2.35452, rel=1.0e-3)


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.usefixtures("clean_astro_ui")
def test_thread_pileup(run_thread):

    models = run_thread('pileup')
    fr = ui.get_fit_results()
    covarerr = sqrt(fr.extra_output['covar'].diagonal())

    assert fr.statval == approx(53.6112, rel=1e-4)
    assert fr.rstat == approx(1.44895, rel=1e-4)
    assert fr.qval == approx(0.0379417, rel=1e-4)
    assert fr.numpoints == 42
    assert fr.dof == 37

    jdp = models['jdp']
    assert jdp.alpha.val == approx(0.522593, rel=1e-1)
    assert jdp.f.val == approx(0.913458, rel=1e-2)

    abs1 = models['abs1']
    assert abs1.nh.val == approx(6.12101, rel=1e-2)

    power = models['power']
    assert power.gamma.val == approx(1.41887, rel=1e-2)
    assert power.ampl.val == approx(0.00199457, rel=1e-2)

    # Move covariance checks to the end of the parameters
    # to check whether the fit is the same, even if the
    # covariance is different.
    #
    # Note that these checks seem very sensitive to the configuration
    # of the system (e.g. XSPEC version), so DJB has commented them
    # out as it is not clear how much value they provide (e.g.
    # the alpha value is ~0.52 and has an error ~ 680 - 690, which
    # is not very useful). Is this a good test case?
    #
    # assert covarerr[0] == approx(684.056, rel=1e-4)
    # assert covarerr[1] == approx(191.055, rel=1e-3)
    assert covarerr[2] == approx(0.632061, rel=1e-3)
    assert covarerr[3] == approx(0.290159, rel=1e-3)
    # assert covarerr[4] == approx(1.62529, rel=1e-3)

    # Issue #294 was a problem with serializing the pileup model
    # after a fit in Python 3 (but not Python 2). Add some basic
    # validation that the conversion to a string works. For the
    # pileup model we expect the standard model layout - e.g.
    #
    #   jdp
    #   parameter headers
    #   ---- ----- ...
    #   jdp.alpha ...
    #   ...
    #   jdp.nterms ...
    #   <blank line>
    #   1: ...
    #   ...
    #   7: ...
    #   *** pileup fraction: value
    #
    lines = str(jdp).split('\n')
    assert len(lines) == 19
    assert lines[10].strip() == ''
    assert lines[11].startswith('   1: ')
    assert lines[17].startswith('   7: ')
    assert lines[18].startswith('   *** pileup fraction: ')


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize('parallel', [False, True])
def test_proj_bubble(parallel, clean_astro_ui, fix_xspec, run_thread):

    # How sensitive are the results to the change from bcmc to vern
    # made in XSPEC 12.10.1? It looks like the mekal best-fit
    # temperature can jump from ~17.9 to 18.6, so require bcmc
    # in this test (handled by fix_xspec).
    #
    # Note that the error on kT is large, so we can expect that
    # changes to the system could change these results. In particular,
    # the covariance errors on kt are < 1 but from other error
    # analysis they are > 10 or even unbound, so it is likely that the
    # covariance error results can change significantly.
    #
    # The fit results change in XSPEC 12.10.0 since the mekal model
    # was changed (FORTRAN to C++). A 1% difference is used for the
    # parameter ranges from covar and proj (matches the tolerance for
    # the fit results). Note that this tolerance has been relaced to
    # 10% for the kT errors, as there is a significant change seen
    # with different XSPEC versions for the covariance results.
    #
    # fit_results is unused
    def cmp_thread(fit_results, mek1, covarerr):

        assert covarerr[0] == approx(0, rel=1e-4)
        assert covarerr[1] == approx(8.74608e-07, rel=1e-3)
        assert mek1.kt.val == approx(17.8849, rel=1e-2)
        assert mek1.norm.val == approx(4.15418e-06, rel=1e-2)

    check_thread(run_thread, 'proj_bubble', parallel, cmp_thread,
                 ['mek1'])

    # Covar
    #
    # TODO: should this check that parmaxes is -1 * parmins instead?
    covar = ui.get_covar_results()
    assert covar.parmins[0] == approx(-0.653884, rel=0.1)
    assert covar.parmins[1] == approx(-8.94436e-07, rel=0.01)
    assert covar.parmaxes[0] == approx(0.653884, rel=0.1)
    assert covar.parmaxes[1] == approx(8.94436e-07, rel=0.01)

    # Proj -- Upper bound of kT can't be found
    #
    proj = ui.get_proj_results()
    assert proj.parmins[0] == approx(-12.048069, rel=0.01)
    assert proj.parmins[1] == approx(-9.510913e-07, rel=0.01)
    assert proj.parmaxes[0] is None
    assert proj.parmaxes[1] == approx(2.403640e-06, rel=0.01)

    # Conf
    #
    conf = ui.get_conf_results()
    assert conf.parmins[0] == approx(-12.1073, rel=0.01)
    assert conf.parmins[1] == approx(-9.5568e-07, rel=0.01)
    assert conf.parmaxes[0] is None
    assert conf.parmaxes[1] == approx(2.39937e-06, rel=0.01)
