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

from sherpa.utils.test import SherpaTest, SherpaTestCase, \
    requires_data, requires_xspec, requires_fits

from sherpa.models import PowLaw1D
from sherpa.fit import Fit
from sherpa.stats import Cash, UserStat, WStat
from sherpa.optmethods import LevMar, NelderMead
from sherpa.utils.err import StatErr

import logging
logger = logging.getLogger("sherpa")


class MyCashWithBkg(UserStat):

    def __init__(self, name='mycash'):
        UserStat.__init__(self, name)

    @staticmethod
    def mycal_staterror(data):
        return None

    @staticmethod
    def cash_withbkg(data, model, staterror=None, syserror=None, weight=None,
                     bkg=None):
        fvec = model - (data * numpy.log(model))
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
    def chi_withbkg(data, model, staterror=None, syserror=None, weight=None,
                    bkg=None):
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
    def cash_nobkg(data, model, staterror=None, syserror=None,
                   weight=None):
        fvec = model - (data * numpy.log(model))
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
    def chi_nobkg(data, model, staterror=None, syserror=None, weight=None):
        fvec = ((data - model) / staterror)**2
        stat = fvec.sum()
        return (stat, fvec)

    calc_stat = chi_nobkg
    calc_staterror = mycal_staterror


@requires_fits
@requires_xspec
@requires_data
class test_stats(SherpaTestCase):

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
        from sherpa.astro.xspec import XSphabs
        from sherpa.astro.io import read_pha

        pha_fname = self.make_path("9774.pi")
        self.data = read_pha(pha_fname)
        self.data.notice(0.5, 7.0)

        bkg_fname = self.make_path("9774_bg.pi")
        self.bkg = read_pha(bkg_fname)

        abs1 = XSphabs('abs1')
        p1 = PowLaw1D('p1')
        self.model = abs1 + p1

    def tearDown(self):
        if hasattr(self, "_old_logger_level"):
            logger.setLevel(self._old_logger_level)

    def compare_results(self, arg1, arg2):

        for key in ["succeeded", "numpoints", "dof"]:
            assert arg1[key] == int(getattr(arg2, key))

        for key in ["istatval", "statval"]:
            assert numpy.allclose(float(arg1[key]), float(getattr(arg2, key)),
                                  1.e-6, 1.e-6)

        for key in ["parvals"]:
            try:
                assert numpy.allclose(arg1[key], getattr(arg2, key),
                                      1.e-3, 1.e-3)
            except AssertionError:
                print 'parvals bench: ', arg1[key]
                print 'parvals fit:   ', getattr(arg2, key)
                print 'results', arg2
                raise

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
        self.compare_results(self._fit_mychinobkg_results_bench, results)

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
        self.compare_results(self._fit_mychinobkg_results_bench, results)

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


def tstme(datadir=None):
    import sherpa.stats as stats
    SherpaTest(stats).test(datadir=datadir)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        datadir = sys.argv[1]
    else:
        datadir = '/data/scialg/testdata'
    if not os.path.exists(datadir):
        datadir = None
    tstme(datadir)
