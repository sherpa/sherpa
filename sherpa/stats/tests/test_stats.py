from __future__ import print_function
#
#  Copyright (C) 2007, 2015, 2016  Smithsonian Astrophysical Observatory
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

import os.path

import numpy

from sherpa.utils import requires_data, requires_xspec, requires_fits
from sherpa.utils import SherpaTestCase
from sherpa.data import Data1D

from sherpa.models import PowLaw1D, SimulFitModel
from sherpa.data import DataSimulFit
from sherpa.fit import Fit
from sherpa.stats import Stat, Cash, LeastSq, UserStat, WStat, \
    CStat, LeastSq, Chi2Gehrels, Chi2ConstVar, Chi2ModVar, Chi2XspecVar, \
    Chi2DataVar
from sherpa.optmethods import LevMar, NelderMead
from sherpa.utils.err import StatErr
from sherpa.astro import ui
import logging
logger = logging.getLogger("sherpa")


class MySimulStat(UserStat):
    def __init__(self, name='mysimulstat'):
        UserStat.__init__(self, name)

    @staticmethod
    def mycal_staterror(data):
        return numpy.ones_like(data)

    @staticmethod
    def my_simulstat(data, model, staterror, *args, **kwargs):
        data_size = kwargs['extra_args']['data_size']
        data1 = data[:data_size[0]]
        data2 = data[data_size[0]:]
        model1 = model[:data_size[0]]
        model2 = model[data_size[0]:]
        staterror1 = staterror[:data_size[0]]
        staterror2 = staterror[data_size[0]:]
        mystat1 = Chi2DataVar()
        mystat2 = Chi2DataVar()
        stat1, fvec1 = mystat1.calc_stat(data1, model1, staterror1)
        stat2, fvec2 = mystat2.calc_stat(data2, model2, staterror2)
        fvec = numpy.power((data - model) / staterror, 2)
        stat = numpy.sum(fvec)
        # print stat1 + stat2 - stat
        return (stat, fvec)
        return (stat1 + stat2, numpy.append(fvec1, fvec2))

    calc_stat = my_simulstat
    calc_staterror = mycal_staterror


class MyCashWithBkg(UserStat):

    def __init__(self, name='mycash'):
        UserStat.__init__(self, name)

    @staticmethod
    def mycal_staterror(data):
        return None

    @staticmethod
    def cash_withbkg(data, model, staterror, *args, **kwargs):
        fvec = model - (data * numpy.log(model))
        weight = kwargs.get('weight')
        if weight is not None:
            fvec = fvec * weight
        return 2.0 * sum(fvec), fvec

    calc_stat = cash_withbkg
    calc_staterror = mycal_staterror


class MyChiWithBkg(UserStat):

    def __init__(self, name='mychi'):
        UserStat.__init__(self, name)

    @staticmethod
    def mycal_staterror(data):
        return numpy.ones_like(data)

    @staticmethod
    def chi_withbkg(data, model, staterror, *args, **kwargs):
        fvec = ((data - model) / staterror)**2
        stat = fvec.sum()
        return (stat, fvec)

    calc_stat = chi_withbkg
    calc_staterror = mycal_staterror


class MyCashNoBkg(UserStat):

    def __init__(self, name='mycash'):
        UserStat.__init__(self, name)

    @staticmethod
    def mycal_staterror(data):
        return None

    @staticmethod
    def cash_nobkg(data, model, staterror, *args, **kwargs):
        fvec = model - (data * numpy.log(model))
        weight = kwargs.get('weight')
        if weight is not None:
            fvec = fvec * weight
        return 2.0 * sum(fvec), fvec

    calc_stat = cash_nobkg
    calc_staterror = mycal_staterror


class MyChiNoBkg(UserStat):

    def __init__(self, name='mychi'):
        UserStat.__init__(self, name)

    @staticmethod
    def mycal_staterror(data):
        return numpy.ones_like(data)

    @staticmethod
    def chi_nobkg(data, model, staterror, *args, **kwargs):
        fvec = ((data - model) / staterror)**2
        stat = fvec.sum()
        return (stat, fvec)

    calc_stat = chi_nobkg
    calc_staterror = mycal_staterror


@requires_fits
@requires_xspec
@requires_data
class test_stats(SherpaTestCase):

    _fit_simul_datavarstat_results_bench = {
        'succeeded': 1,
        'numpoints': 18,
        'dof': 15,
        'istatval': 1218.11457171,
        'statval': 204.883073969,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [65647.539439194588, 2.1440354994101929, 13955.023665227312])
        }

    _fit_chi2xspecvar_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 3488.58436522,
        'statval': 1167.11621982,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [4717.3082876288863, 1.785895698907098, 39702.914274813716])
        }

    _fit_chi2modvar_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 98751.1141165,
        'statval': 951.052518517,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [6323.954237402093, 1.5898717247339578, 25049.100925267721])
        }

    _fit_chi2constvar_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 11078.2610904,
        'statval': 1658.30318937,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [5408.2022582955542, 1.2222311487955908, 4612.0008873213301])
        }

    _fit_chi2gehrels_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 2295.18738409,
        'statval': 590.888903039,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [5077.8010218337085, 1.592875823400443, 19067.111802328174])
        }

    _fit_leastsq_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 100275.650273,
        'statval': 15010.2465747,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [5408.4207593703177, 1.2222288063374338, 4611.9525063414185])
        }

    _fit_cstat_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 16859.677457,
        'statval': 1173.95573689,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [5886.0660236942495, 1.6556198746259132, 30098.968589487202])
        }

    _fit_mycash_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 796.401435754,
        'statval': -14889.3202844,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [5886.0660236942495, 1.6556198746259132, 30098.968589487202])
        }

    _fit_mycashnobkg_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 796.401435754,
        'statval': -14889.3202844,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [5886.0660236942495, 1.6556198746259132, 30098.968589487202])
        }

    _fit_mychi_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 100275.650273,
        'statval': 15082.4817361,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [65.215835020062741, 1.2149346471169165, 4454.4695930173866])
        }

    _fit_mycashnobkg_results_bench = {
        'succeeded': 1,
        'numpoints': 1024,
        'dof': 1021,
        'istatval': 2198.3631781,
        'statval': 1716.74869273,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [295.11120384933781, 0.69990055680397523, 20.998971817852862])
        }

    _fit_mychinobkg_results_bench = {
        'succeeded': 1,
        'numpoints': 1024,
        'dof': 1021,
        'istatval': 5664.2486547,
        'statval': 7928.05674899,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [346.51084808235697, 0.24721168701021015, 7.9993714921823997])
        }

    _fit_wstat_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 14000.5250801,
        'statval': 1156.57701288,
        'rstat': 2.6107833248,
        'parnames': ('abs1.nH', 'abs1.gamma', 'abs1.ampl'),
        'parvals': numpy.array(
            [5864.278543739505, 1.6569575154646112, 29868.225197035885])}

    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        from sherpa.astro.io import read_pha
        from sherpa.astro.xspec import XSphabs

        pha_fname = self.make_path("9774.pi")
        self.data = read_pha(pha_fname)
        self.data.notice(0.5, 7.0)

        bkg_fname = self.make_path("9774_bg.pi")
        self.bkg = read_pha(bkg_fname)

        abs1 = XSphabs('abs1')
        p1 = PowLaw1D('p1')
        self.model = abs1 + p1

        self.model_mult = abs1 * p1
        pi2278 = self.make_path("pi2278.fits")
        pi2286 = self.make_path("pi2286.fits")
        self.data_pi2278 = read_pha(pi2278)
        self.data_pi2286 = read_pha(pi2286)

    def tearDown(self):
        if hasattr(self, "_old_logger_level"):
            logger.setLevel(self._old_logger_level)

    def compare_results(self, arg1, arg2, tol1=1.0e-6, tol2=1.0e-3):

        for key in ["succeeded", "numpoints", "dof"]:
            assert arg1[key] == int(getattr(arg2, key))

        for key in ["istatval", "statval"]:
            assert numpy.allclose(float(arg1[key]), float(getattr(arg2, key)),
                                  tol1, tol1)

        for key in ["parvals"]:
            try:
                assert numpy.allclose(arg1[key], getattr(arg2, key),
                                      tol2, tol2)
            except AssertionError:
                print('parvals bench: ', arg1[key])
                print('parvals fit:   ', getattr(arg2, key))
                print('results', arg2)
                raise

    def test_chi2xspecvar_stat(self):
        fit = Fit(self.data, self.model, Chi2XspecVar(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_chi2xspecvar_results_bench, results)

    def test_chi2modvar_stat(self):
        fit = Fit(self.data, self.model, Chi2ModVar(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_chi2modvar_results_bench, results)

    def test_chi2constvar_stat(self):
        fit = Fit(self.data, self.model, Chi2ConstVar(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_chi2constvar_results_bench, results)

    def test_chi2gehrels_stat(self):
        fit = Fit(self.data, self.model, Chi2Gehrels(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_chi2gehrels_results_bench, results)

    def test_leastsq_stat(self):
        fit = Fit(self.data, self.model, LeastSq(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_leastsq_results_bench, results)

    def test_cstat_stat(self):
        fit = Fit(self.data, self.model, CStat(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_cstat_results_bench, results)

    def test_cash_stat(self):
        fit = Fit(self.data, self.model, Cash(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_mycash_results_bench, results)

    def test_mycash_data_and_model_have_bkg(self):
        fit = Fit(self.data, self.model, MyCashWithBkg(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_mycash_results_bench, results)

    def test_mychi_data_and_model_have_bkg(self):
        fit = Fit(self.data, self.model, MyChiWithBkg(), LevMar())
        results = fit.fit()
        self.compare_results(self._fit_mychi_results_bench, results)

    def test_mycash_data_and_model_donothave_bkg(self):
        data = self.bkg
        fit = Fit(data, self.model, MyCashNoBkg(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_mycashnobkg_results_bench, results)

    def test_mychi_data_and_model_donothave_bkg(self):
        data = self.bkg
        fit = Fit(data, self.model, MyChiNoBkg(), LevMar())
        results = fit.fit()
        self.compare_results(self._fit_mychinobkg_results_bench, results,
                             tol1=1.0e-5)

    def test_mycash_datahasbkg_modelhasnobkg(self):
        fit = Fit(self.data, self.model, MyCashNoBkg(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_mycash_results_bench, results)

    def test_mycash_nobkgdata_modelhasbkg(self):
        data = self.bkg
        fit = Fit(data, self.model, MyCashWithBkg(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_mycashnobkg_results_bench, results)

    def test_mychi_datahasbkg_modelhasnobkg(self):
        fit = Fit(self.data, self.model, MyChiNoBkg(), LevMar())
        results = fit.fit()
        self.compare_results(self._fit_mychi_results_bench, results)

    def test_mychi_nobkgdata_modelhasbkg(self):
        data = self.bkg
        fit = Fit(data, self.model, MyChiWithBkg(), LevMar())
        results = fit.fit()
        self.compare_results(self._fit_mychinobkg_results_bench, results,
                             tol1=1.0e-5)

    def test_wstat(self):
        fit = Fit(self.data, self.model, WStat(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_wstat_results_bench, results)

    # The following test passes if run by itself but fails when run with others
    # def test_wstat1(self):
    #     pha_fname = self.make_path("stats/acisf09122_000N001_r0013_pha3.fits")
    #     ui.load_pha(pha_fname)
    #     #ui.set_analysis('energy')
    #     ui.ignore(None, None)
    #     ui.notice(1.0, 1.6)
    #     src = ui.xsphabs.gal * ui.xspowerlaw.pl
    #     gal.nh = 0.1
    #     pl.phoindex = 0.7
    #     pl.norm = 1e-4
    #     ui.set_source(src)
    #     ui.set_stat('wstat')
    #     assert numpy.allclose(46.455049531, ui.calc_stat(), 1.e-7, 1.e-7)

    def test_wstat_error(self):
        data = self.bkg
        data.notice(0.5, 7.0)
        fit = Fit(data, self.model, WStat(), NelderMead())
        self.assertRaises(StatErr, fit.fit)

    def test_chi2datavar(self):
        num = 3
        xy = numpy.array(range(num))
        ui.load_arrays(1, xy, xy, Data1D)
        ui.set_stat('chi2datavar')
        err = ui.get_staterror()
        numpy.testing.assert_allclose(err, numpy.sqrt(xy),
                                      rtol=1e-7, atol=1e-7)

    def test_get_stat_info(self):
        fname_3c273 = self.make_path("3c273.pi")
        ui.load_pha(fname_3c273)
        src = ui.xspowerlaw.pl
        ui.set_source(src)
        ui.guess('pl')
        ui.set_stat('wstat')
        stat_info = ui.get_stat_info()[0]
        assert stat_info.dof == 44
        assert stat_info.numpoints == 46

    def test_mysimul_stat_fit(self):
        data1 = self.data_pi2278
        data2 = self.data_pi2286
        model1 = self.model_mult
        model2 = self.model_mult
        data = DataSimulFit(name='data1data2', datasets=[data1, data2])
        model = SimulFitModel(name='model1model2', parts=[model1, model2])
        fit = Fit(data=data, model=model, stat=MySimulStat(),
                  method=NelderMead())
        result = fit.fit()
        self.compare_results(self._fit_simul_datavarstat_results_bench,
                             result)

    def test_simul_stat_fit(self):
        data1 = self.data_pi2278
        data2 = self.data_pi2286
        model1 = self.model_mult
        model2 = self.model_mult
        data = DataSimulFit(name='data1data2', datasets=[data1, data2])
        model = SimulFitModel(name='model1model2', parts=[model1, model2])
        fit = Fit(data=data, model=model, stat=MyChiNoBkg(),
                  method=NelderMead())
        result = fit.fit()
        self.compare_results(self._fit_simul_datavarstat_results_bench,
                             result)

    # bug #147
    def test_wstat_calc_stat_info(self):
        ui.load_pha("stat", self.make_path("3c273.pi"))
        ui.set_source("stat", ui.powlaw1d.p1)
        ui.set_stat("wstat")
        ui.fit("stat")
        ui.get_stat_info()
