#
#  Copyright (C) 2018, 2019, 2021  Smithsonian Astrophysical Observatory
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

import pytest

from sherpa.utils.testing import requires_data, requires_fits, requires_xspec

from sherpa.models import Polynom1D, SimulFitModel
from sherpa.models.basic import Gauss1D
from sherpa.data import Data1D, DataSimulFit
from sherpa.optmethods import NelderMead, LevMar
from sherpa.stats import Cash, LeastSq
from sherpa.fit import Fit
from sherpa.astro import ui


@requires_data
@requires_fits
@requires_xspec
def test_ARFModelPHA(make_data_path, clean_astro_ui):
    ui.load_pha(make_data_path("3c120_meg_1.pha"))

    # remove the RMF to ensure this is an ARF-only analysis
    # (which is what is needed to trigger the bug that lead to #699)
    ui.get_data().set_rmf(None)

    ui.group_counts(20)
    ui.notice(0.5, 6)
    ui.subtract()
    ui.set_model(ui.xsphabs.abs1 * (ui.xsapec.bubble + ui.powlaw1d.p1))
    ui.set_xsabund('angr')
    ui.set_xsxsect('vern')
    abs1.nh = 0.163
    abs1.nh.freeze()
    p1.ampl = 0.017
    p1.gamma = 1.9
    bubble.kt = 0.5
    bubble.norm = 4.2e-5
    tol = 1.0e-2
    ui.set_method_opt('ftol', tol)
    ui.fit()
    result = ui.get_fit_results()
    assert result.numpoints == 1212
    assert result.dof == 1208


fit_same_poly_bench = {
    'numpoints': 153,
    'dof': 151,
    'istatval': 306.0,
    'statval': -19417.55547762454,
    'parnames': ('polynom1d.c0', 'polynom1d.c1'),
    'parvals': numpy.array([1.6874869548492573, 0.4796635263358297])
    }

fit_diff_poly_bench = {
    'numpoints': 153,
    'dof': 147,
    'istatval': 306.0,
    'statval': -19418.614329430166,
    'parnames': ('polynom1d.c0', 'polynom1d.c1', 'polynom1d.c0',
                 'polynom1d.c1', 'polynom1d.c0', 'polynom1d.c1'),
    'parvals': numpy.array([1.755035876061136, 0.46847103779088173,
                            1.9361738167035432, 0.4791446119040638,
                            1.342017385707699, 0.49194814027685246])
    }

fit_g2g2_bench = {
    'numpoints': 200,
    'dof': 194,
    'istatval': 4.274899468449953,
    'statval': 0.504555457862299,
    'parnames': ('gauss1d.fwhm', 'gauss1d.pos', 'gauss1d.ampl',
                 'gauss1d.fwhm', 'gauss1d.pos', 'gauss1d.ampl'),
    'parvals': numpy.array([1.1421173550314392, 1.1963309074397634,
                            0.9677194890720789, 4.1688184531922765,
                            -1.3386595864803255, 1.0047067939881642])
    }


@pytest.fixture
def setUp3(hide_logging):

    x = numpy.linspace(1.0, 101., num=101)[0::2]
    y1 = [1., 5., 2., 4., 7., 11., 9., 8., 12., 18., 12., 11., 13., 12., 13.,
          13., 20., 23., 16., 20., 24., 17., 21., 26., 22., 24., 24., 21., 28.,
          28., 26., 25., 34., 26., 34., 33., 25., 38., 31., 43., 35., 42., 50.,
          41., 43., 47., 57., 53., 60., 46., 54.]
    y2 = [0., 7., 6., 3., 5., 5., 9., 11., 13., 8., 14., 13., 14., 18., 11.,
          15., 17., 26., 15., 19., 25., 30., 15., 29., 16., 25., 27., 29., 36.,
          41., 22., 27., 33., 32., 45., 37., 38., 38., 34., 52., 40., 41., 31.,
          47., 38., 52., 57., 33., 48., 53., 45.]
    y3 = [1., 2., 4., 2., 5., 8., 15., 10., 13., 10., 16., 10., 13., 12., 16.,
          17., 17., 20., 23., 16., 25., 22., 19., 31., 26., 24., 21., 29., 36.,
          30., 33., 30., 37., 27., 36., 32., 42., 44., 39., 30., 40., 33., 39.,
          49., 56., 47., 46., 35., 63., 40., 57.]
    d1 = Data1D('1', x, y1)
    d2 = Data1D('2', x, y2)
    d3 = Data1D('3', x, y3)

    return d1, d2, d3


@pytest.fixture
def setUp2(hide_logging, reset_seed):

    x = numpy.linspace(-5., 5., 100)
    g1, g2 = Gauss1D(), Gauss1D()
    g1.fwhm = 1.14
    g1.pos = 1.2
    g2.fwhm = 4.13
    g2.pos = -1.3

    numpy.random.seed(0)
    y1 = g1(x) + numpy.random.normal(0.0, 0.05, x.shape)
    y2 = g2(x) + numpy.random.normal(0.0, 0.05, x.shape)

    d4 = Data1D('4', x, y1)
    d5 = Data1D('5', x, y2)

    return d4, d5


def compare_results(expected, got, tol=1e-4):

    for key in ["numpoints", "dof"]:
        assert int(getattr(got, key)) == expected[key]

    for key in ["istatval", "statval"]:
        assert float(getattr(got, key)) == pytest.approx(expected[key], rel=tol)

    assert got.parvals == pytest.approx(expected['parvals'], rel=tol)


def test_diff_cache(setUp3):
    d1, d2, d3 = setUp3

    poly1 = Polynom1D()
    poly2 = Polynom1D()
    poly3 = Polynom1D()
    poly1.pars[1].thaw()
    poly2.pars[1].thaw()
    poly3.pars[1].thaw()
    sdata = DataSimulFit('d123', (d1, d2, d3))
    smodel = SimulFitModel('diff', (poly1, poly2, poly3))
    sfit = Fit(sdata, smodel, method=NelderMead(), stat=Cash())
    result = sfit.fit()
    compare_results(fit_diff_poly_bench, result)


def test_same_cache(setUp3):
    d1, d2, d3 = setUp3

    poly = Polynom1D()
    poly.pars[1].thaw()
    sdata = DataSimulFit('d1d2d3', (d1, d2, d3))
    smodel = SimulFitModel('same', (poly, poly, poly))
    sfit = Fit(sdata, smodel, method=NelderMead(), stat=Cash())
    result = sfit.fit()
    compare_results(fit_same_poly_bench, result)


def test_gauss_gauss(setUp2):
    d4, d5 = setUp2

    g1, g2 = Gauss1D(), Gauss1D()
    g1.fwhm = 1.3
    g1.pos = 1.5
    g2.fwhm = 4.
    g2.pos = -2.0
    sdata = DataSimulFit('d4d5', (d4, d5))
    smodel = SimulFitModel('g1g2', (g1, g2))
    sfit = Fit(sdata, smodel, method=LevMar(), stat=LeastSq())
    result = sfit.fit()
    compare_results(fit_g2g2_bench, result)


def test_cache_copy(clean_astro_ui):

    # fake up a PHA data set
    chans = numpy.arange(1, 11, dtype=numpy.int8)
    counts = numpy.ones(chans.size)

    # bin width is not 0.1 but something slightly different
    ebins = numpy.linspace(0.1, 1.2, num=chans.size + 1)
    elo = ebins[:-1]
    ehi = ebins[1:]

    dset = ui.DataPHA('test', chans, counts)

    # make sure ARF isn't 1 to make sure it's being applied
    arf = ui.create_arf(elo, ehi, specresp=0.7 * numpy.ones(chans.size))

    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    ui.set_data(1, dset)
    ui.set_arf(1, arf)
    ui.set_rmf(1, rmf)

    ui.set_source(ui.const1d.mdl)

    # again not 1
    mdl.c0 = 8

    # Copy the values from the plot structures, since get_xxx_plot
    # returns the same object so m1.y == m2.y will not note a difference.
    #
    d1y = ui.get_data_plot().y.copy()
    m1y = ui.get_model_plot().y.copy()
    s1y = ui.get_source_plot().y.copy()

    d2y = ui.get_data_plot().y.copy()
    m2y = ui.get_model_plot().y.copy()
    s2y = ui.get_source_plot().y.copy()

    assert d1y == pytest.approx(d2y)
    assert m1y == pytest.approx(m2y)
    assert s1y == pytest.approx(s2y)
