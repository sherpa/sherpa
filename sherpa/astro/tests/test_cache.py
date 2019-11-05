from __future__ import print_function

#
#  Copyright (C) 2018, 2019  Smithsonian Astrophysical Observatory
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

from sherpa.utils.testing import SherpaTestCase, requires_data, requires_fits,\
    requires_xspec

from sherpa.models import Polynom1D, SimulFitModel
from sherpa.models.basic import Gauss1D
from sherpa.data import Data1D, DataSimulFit
from sherpa.optmethods import NelderMead, LevMar
from sherpa.stats import Cash, LeastSq
from sherpa.fit import Fit

import logging
logger = logging.getLogger("sherpa")

@requires_data
@requires_fits
@requires_xspec
class test_ARFModelPHA(SherpaTestCase):
    _fit_using_ARFModelPHA = {
        'numpoints': 1212,
        'dof': 1208
    }
    
    def setup(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        
    def tearDown(self):
        if hasattr(self, "_old_logger_level"):
            logger.setLevel(self._old_logger_level)

    def test_ARFModelPHA(self):
        from sherpa.astro import ui
        ui.load_pha(self.make_path("3c120_meg_1.pha"))
        ui.group_counts(20)
        ui.notice(0.5, 6)
        ui.subtract()
        ui.set_model(ui.xsphabs.abs1 * (ui.xsapec.bubble + ui.powlaw1d.p1))
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
        assert result.numpoints == self._fit_using_ARFModelPHA['numpoints']
        assert result.dof == self._fit_using_ARFModelPHA['dof']

class test_cache(SherpaTestCase):

    _fit_same_poly_bench = {
        'numpoints': 153,
        'dof': 151,
        'istatval': 306.0,
        'statval': -19417.55547762454,
        'parnames': ('polynom1d.c0', 'polynom1d.c1'),
        'parvals': numpy.array([1.6874869548492573, 0.4796635263358297])
        }

    _fit_diff_poly_bench = {
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

    _fit_g2g2_bench = {
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

    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        x = numpy.linspace(1.0, 101., num=101)[0::2]
        y1 = [ 1., 5., 2., 4., 7.,11., 9., 8.,12.,18.,12.,11.,13.,12.,13.,
               13.,20.,23.,16.,20.,24.,17.,21.,26.,22.,24.,24.,21.,28.,
               28.,26.,25.,34.,26.,34.,33.,25.,38.,31.,43.,35.,42.,50.,
               41.,43.,47.,57.,53.,60.,46.,54.]
        y2 = [ 0., 7., 6., 3., 5., 5., 9.,11.,13., 8.,14.,13.,14.,18.,11.,
               15.,17.,26., 15.,19.,25.,30.,15.,29.,16.,25.,27.,29.,36.,
               41.,22.,27.,33.,32.,45.,37.,38.,38.,34.,52.,40.,41.,31.,
               47.,38.,52.,57.,33.,48.,53.,45.]
        y3 = [ 1., 2., 4., 2., 5., 8.,15.,10.,13.,10.,16.,10.,13.,12.,16.,
               17.,17.,20., 23.,16.,25.,22.,19.,31.,26.,24.,21.,29.,36.,
               30.,33.,30.,37.,27.,36.,32., 42.,44.,39.,30.,40.,33.,39.,
               49.,56.,47.,46.,35.,63.,40.,57.]
        self.d1 = Data1D('1', x, y1)
        self.d2 = Data1D('2', x, y2)
        self.d3 = Data1D('3', x, y3)

        x = numpy.linspace(-5., 5., 100)
        g1, g2 = Gauss1D(), Gauss1D()
        g1.fwhm = 1.14
        g1.pos = 1.2
        g2.fwhm = 4.13
        g2.pos = -1.3
        numpy.random.seed(0)
        y1 = g1(x) + numpy.random.normal(0.0, 0.05, x.shape)
        y2 = g2(x) + numpy.random.normal(0.0, 0.05, x.shape)
        self.d4 = Data1D('4', x, y1)
        self.d5 = Data1D('5', x, y2)

    def tearDown(self):
        if hasattr(self, "_old_logger_level"):
            logger.setLevel(self._old_logger_level)

    def compare_results(self, arg1, arg2, tol=1e-4):

        for key in ["numpoints", "dof"]:
            assert arg1[key] == int(getattr(arg2, key))

        for key in ["istatval", "statval"]:
            numpy.testing.assert_allclose(float(arg1[key]),
                                          float(getattr(arg2, key)), tol)

        for key in ["parvals"]:
            try:
                numpy.testing.assert_allclose(arg1[key],
                                              getattr(arg2, key), tol)
            except AssertionError:
                print('parvals bench: ', arg1[key])
                print('parvals fit:   ', getattr(arg2, key))
                print('results', arg2)
                raise

    def test_diff_cache(self):
        poly1 = Polynom1D()
        poly2 = Polynom1D()
        poly3 = Polynom1D()
        poly1.pars[1].thaw()
        poly2.pars[1].thaw()
        poly3.pars[1].thaw()
        sdata = DataSimulFit('d123', (self.d1, self.d2, self.d3))
        smodel = SimulFitModel('diff', (poly1, poly2, poly3))
        sfit = Fit(sdata, smodel, method=NelderMead(), stat=Cash())
        result = sfit.fit()
        self.compare_results(self._fit_diff_poly_bench, result)

    def test_same_cache(self):
        poly = Polynom1D()
        poly.pars[1].thaw()
        sdata = DataSimulFit('d1d2d3', (self.d1, self.d2, self.d3))
        smodel = SimulFitModel('same', (poly, poly, poly))
        sfit = Fit(sdata, smodel, method=NelderMead(), stat=Cash())
        result = sfit.fit()
        self.compare_results(self._fit_same_poly_bench, result)

    def test_gauss_gauss(self):
        g1, g2 = Gauss1D(), Gauss1D()
        g1.fwhm = 1.3
        g1.pos = 1.5
        g2.fwhm = 4.
        g2.pos = -2.0
        sdata = DataSimulFit('d4d5', (self.d4, self.d5))
        smodel = SimulFitModel('g1g2', (g1, g2))
        sfit = Fit(sdata, smodel, method=LevMar(), stat=LeastSq())
        result = sfit.fit()
        self.compare_results(self._fit_g2g2_bench, result)
