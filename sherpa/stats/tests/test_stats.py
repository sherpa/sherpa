#
#  Copyright (C) 2007, 2015  Smithsonian Astrophysical Observatory
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
import unittest

from sherpa.utils import SherpaTest, SherpaTestCase, test_data_missing
from sherpa.utils import has_package_from_list, has_fits_support

from sherpa.models import PowLaw1D
from sherpa.fit import Fit
from sherpa.stats import Stat, Cash, LeastSq, UserStat
from sherpa.optmethods import LevMar, NelderMead


class MyCash(UserStat):

    def __init__(self, name='mycash'):
        UserStat.__init__(self, name)

    @staticmethod
    def mycal_staterror(data):
        return None

    @staticmethod
    def mycalc_stat(data, model, staterror=None, syserror=None, weight=None,
                    bkg=None):
        assert bkg is not None
        fvec = model - (data * numpy.log(model))
        if weight is not None:
            fvec = fvec * weight
        return 2.0 * sum(fvec), fvec

    calc_stat = mycalc_stat
    calc_staterror = mycal_staterror


class MyChi(UserStat):

    def __init__(self, name='mychi'):
        UserStat.__init__(self, name)

    @staticmethod
    def mycal_staterror(data):
        return numpy.ones_like(data)

    @staticmethod
    def mycalc_stat(data, model, staterror=None, syserror=None, weight=None,
                    bkg=None):
        assert bkg is not None
        fvec = ((data - model) / staterror)**2
        stat = fvec.sum()
        return (stat, fvec)

    calc_stat = mycalc_stat
    calc_staterror = mycal_staterror


class MyCashNoBkg(UserStat):

    def __init__(self, name='mycash'):
        UserStat.__init__(self, name)

    @staticmethod
    def mycal_staterror(data):
        return None

    @staticmethod
    def mycalc_stat(data, model, staterror=None, syserror=None,
                    weight=None):
        fvec = model - (data * numpy.log(model))
        if weight is not None:
            fvec = fvec * weight
        return 2.0 * sum(fvec), fvec

    calc_stat = mycalc_stat
    calc_staterror = mycal_staterror


class MyChiNoBkg(UserStat):

    def __init__(self, name='mychi'):
        UserStat.__init__(self, name)

    @staticmethod
    def mycal_staterror(data):
        return numpy.ones_like(data)

    @staticmethod
    def mycalc_stat(data, model, staterror=None, syserror=None, weight=None):
        fvec = ((data - model) / staterror)**2
        stat = fvec.sum()
        return (stat, fvec)

    calc_stat = mycalc_stat
    calc_staterror = mycal_staterror


@unittest.skipIf(not has_fits_support(),
                 'need pycrates, pyfits or astropy.io.fits')
@unittest.skipIf(not has_package_from_list('sherpa.astro.xspec'),
                 "required sherpa.astro.xspec module missing")
@unittest.skipIf(test_data_missing(), "required test data missing")
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

    def setUp(self):
        try:
            from sherpa.astro.xspec import XSphabs
            from sherpa.astro.io import read_pha
        except:
            return

        pha_fname = self.make_path("stats/9774.pi")
        self.data = read_pha(pha_fname)
        self.data.notice(0.5, 7.0)

        bkg_fname = self.make_path("stats/9774_bg.pi")
        self.bkg = read_pha(bkg_fname)

        abs1 = XSphabs('abs1')
        p1 = PowLaw1D('p1')
        self.model = abs1 + p1

    def compare_results(self, arg1, arg2):

        for key in ["succeeded", "numpoints", "dof"]:
            assert arg1[key] == int(getattr(arg2, key))

        for key in ["istatval", "statval"]:
            assert numpy.allclose(float(arg1[key]), float(getattr(arg2, key)),
                                  1.e-7, 1.e-7)

        for key in ["parvals"]:
            try:
                assert numpy.allclose(arg1[key], getattr(arg2, key),
                                      1.e-4, 1.e-4)
            except AssertionError:
                print 'parvals bench: ', arg1[key]
                print 'parvals fit:   ', getattr(arg2, key)
                print 'results', arg2
                raise

    def test_cash_stat(self):
        fit = Fit(self.data, self.model, Cash(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_mycash_results_bench, results)

    def test_mycash_stat(self):
        fit = Fit(self.data, self.model, MyCash(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_mycash_results_bench, results)

    def test_mychi_stat(self):
        fit = Fit(self.data, self.model, MyChi(), LevMar())
        results = fit.fit()
        self.compare_results(self._fit_mychi_results_bench, results)

    def test_cashnobkg(self):
        data = self.bkg
        fit = Fit(data, self.model, MyCashNoBkg(), NelderMead())
        results = fit.fit()
        self.compare_results(self._fit_mycashnobkg_results_bench, results)

    def test_chinobkg(self):
        data = self.bkg
        fit = Fit(data, self.model, MyChiNoBkg(), LevMar())
        results = fit.fit()
        self.compare_results(self._fit_mychinobkg_results_bench, results)


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
