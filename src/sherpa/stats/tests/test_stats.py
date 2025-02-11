#
#  Copyright (C) 2007, 2015, 2016, 2018, 2020 - 2024
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

import logging
from typing import Union

import numpy as np

import pytest

from sherpa.data import Data, Data1D, DataSimulFit
from sherpa.astro.data import DataPHA

from sherpa.astro import ui
from sherpa.astro.instrument import Response1D
from sherpa.data import DataSimulFit
from sherpa.fit import Fit
from sherpa.models import Model, Const1D, PowLaw1D, SimulFitModel
from sherpa.optmethods import LevMar, NelderMead
from sherpa.stats import Cash, CStat, WStat, LeastSq, UserStat, \
    Chi2Gehrels, Chi2ConstVar, Chi2ModVar, Chi2XspecVar, Chi2DataVar, \
    StatResults
from sherpa.utils.err import StatErr
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_group, requires_xspec

logger = logging.getLogger("sherpa")


# Note that these tests extend UserStat in a way not intended, in that
# they overload methods directly rather than taking advantage of
# set_statfnuc/set_errfunc. Perhaps they should just extend Stat
# directly instead?
#
class MySimulStat(UserStat):

    def __init__(self, name: str = 'mysimulstat') -> None:
        super().__init__(name=name)

    @staticmethod
    def calc_staterror(data: np.ndarray) -> np.ndarray:
        return np.ones_like(data)

    # This is based on the original code, but it's not 100% clear
    # why some of the values are being calculated.
    #
    def calc_stat(self,
                  data: Union[Data, DataSimulFit],
                  model: Model
                  ) -> StatResults:

        # This is not generic code, so we can assume we have a
        # simulfit object.
        #
        assert isinstance(data, DataSimulFit)

        tofit = data.to_fit(staterrfunc=self.calc_staterror)
        modeldata = data.eval_model_to_fit(model)

        fitdata = tofit[0]
        staterror = tofit[1]

        fvec = np.power((fitdata - modeldata) / staterror, 2)
        stat = np.sum(fvec)

        mstat = 0.0
        mfvec = []

        # It is not clear what separating the data sets does
        # here, based on the original my_simulsat code.
        #
        stats = [Chi2DataVar(), Chi2DataVar()]
        for mystat, dset, mexpr in zip(stats, data.datasets, model.parts):

            thisstat, thisvec = mystat.calc_stat(dset, mexpr)
            mstat += thisstat
            mfvec.append(thisvec)

        return (stat, fvec)


class MyCashWithBkg(UserStat):

    def __init__(self, name: str = 'mycashwithbkg') -> None:
        super().__init__(name=name)

    @staticmethod
    def calc_staterror(data: np.ndarray) -> np.ndarray:
        return np.zeros_like(data)

    def calc_stat(self,
                  data: Union[Data, DataSimulFit],
                  model: Model
                  ) -> StatResults:

        tofit = data.to_fit(staterrfunc=None)
        modeldata = data.eval_model_to_fit(model)

        fitdata = tofit[0]
        fvec = modeldata - (fitdata * np.log(modeldata))
        return 2.0 * fvec.sum(), fvec


class MyChiWithBkg(UserStat):

    def __init__(self, name: str = 'mychiwithbkg') -> None:
        super().__init__(name=name)

    @staticmethod
    def calc_staterror(data: np.ndarray) -> np.ndarray:
        return np.ones_like(data)

    def calc_stat(self,
                  data: Union[Data, DataSimulFit],
                  model: Model
                  ) -> StatResults:

        tofit = data.to_fit(staterrfunc=self.calc_staterror)
        modeldata = data.eval_model_to_fit(model)

        fitdata = tofit[0]
        staterror = tofit[1]

        fvec = ((fitdata - modeldata) / staterror)**2
        return fvec.sum(), fvec


class MyCashNoBkg(UserStat):

    def __init__(self, name: str = 'mycashnobkg') -> None:
        super().__init__(name=name)

    @staticmethod
    def calc_staterror(data: np.ndarray) -> np.ndarray:
        return np.zeros_like(data)

    def calc_stat(self,
                  data: Union[Data, DataSimulFit],
                  model: Model
                  ) -> StatResults:

        tofit = data.to_fit(staterrfunc=None)
        modeldata = data.eval_model_to_fit(model)

        fitdata = tofit[0]
        fvec = modeldata - (fitdata * np.log(modeldata))

        return 2.0 * fvec.sum(), fvec



class MyChiNoBkg(UserStat):

    def __init__(self, name: str = 'mychinobkg') -> None:
        super().__init__(name=name)

    @staticmethod
    def calc_staterror(data: np.ndarray) -> np.ndarray:
        return np.ones_like(data)

    def calc_stat(self,
                  data: Union[Data, DataSimulFit],
                  model: Model
                  ) -> StatResults:

        tofit = data.to_fit(staterrfunc=self.calc_staterror)
        modeldata = data.eval_model_to_fit(model)

        fitdata = tofit[0]
        staterror = tofit[1]

        fvec = ((fitdata - modeldata) / staterror)**2
        return fvec.sum(), fvec


@pytest.fixture
def reset_xspec():

    from sherpa.astro import xspec

    # Ensure we have a known set of XSPEC settings.
    # At present this is just the abundance and cross-section,
    # since the cosmology settings do not affect any of the
    # models used here.
    #
    abund = xspec.get_xsabund()
    xsect = xspec.get_xsxsect()

    xspec.set_xsabund('angr')
    xspec.set_xsxsect('bcmc')

    yield

    xspec.set_xsabund(abund)
    xspec.set_xsxsect(xsect)


@pytest.fixture
def setup(make_data_path):

    from sherpa.astro.io import read_pha
    from sherpa.astro import xspec

    infile = make_data_path("9774.pi")
    data = read_pha(infile)
    data.notice(0.3, 7.0)

    # Change the exposure time to make the fitted amplitude
    # > 1
    #
    data.exposure = 1

    # Use the wabs model because it is unlikely to change
    # (as scientifically it is no-longer useful). The problem
    # with using something like the phabs model is that
    # changes to that model in XSPEC could change the results
    # here.
    #
    # We fit the log of the nH since this makes the numbers
    # a bit closer to O(1), and so checking should be easier.
    #
    abs1 = xspec.XSwabs('abs1')
    p1 = PowLaw1D('p1')
    factor = Const1D('factor')
    factor.integrate = False
    model = abs1 * p1 + 0 * factor

    factor.c0 = 0
    abs1.nh = 10**factor.c0

    # Ensure the nh limits are honoured by factor (upper limit only).
    # If you don't do this then the fit can fail because a value
    # outside the abs1.nh but within factor.c0 can be picked.
    #
    factor.c0.max = np.log10(abs1.nh.max)

    rsp = Response1D(data)
    return {'data': data, 'model': rsp(model)}


@pytest.fixture
def setup_group(setup):
    data = setup['data']
    data.group_counts(20, tabStops=~data.mask)
    return setup


@pytest.fixture
def setup_bkg(make_data_path):

    from sherpa.astro.io import read_pha, read_arf, read_rmf
    from sherpa.astro import xspec

    infile = make_data_path("9774_bg.pi")
    bkg = read_pha(infile)
    bkg.exposure = 1

    arf = read_arf(make_data_path('9774.arf'))
    rmf = read_rmf(make_data_path('9774.rmf'))
    bkg.set_arf(arf)
    bkg.set_rmf(rmf)

    bkg.set_analysis('energy')
    bkg.notice(0.5, 7.0)

    # We stay with a linear scale for the absorption model
    # here as the values don't seem to go below 0.1.
    #
    abs1 = xspec.XSwabs('abs1')
    p1 = PowLaw1D('p1')
    model = abs1 * p1

    p1.ampl = 1e-4

    rsp = Response1D(bkg)
    return {'bkg': bkg, 'model': rsp(model)}


@pytest.fixture
def setup_bkg_group(setup_bkg):
    data = setup_bkg['bkg']
    data.group_counts(10, tabStops=~data.mask)
    return setup_bkg


@pytest.fixture
def setup_two(make_data_path):

    from sherpa.astro.io import read_pha, read_arf, read_rmf
    from sherpa.astro import xspec

    abs1 = xspec.XSphabs('abs1')
    p1 = PowLaw1D('p1')

    model = abs1 * p1 * 1e-4

    pi2278 = make_data_path("pi2278.fits")
    pi2286 = make_data_path("pi2286.fits")
    data_pi2278 = read_pha(pi2278)
    data_pi2286 = read_pha(pi2286)

    data_pi2278.set_rmf(read_rmf(make_data_path('rmf2278.fits')))
    data_pi2278.set_arf(read_arf(make_data_path('arf2278.fits')))

    data_pi2286.set_rmf(read_rmf(make_data_path('rmf2286.fits')))
    data_pi2286.set_arf(read_arf(make_data_path('arf2286.fits')))

    rsp_pi2278 = Response1D(data_pi2278)
    rsp_pi2286 = Response1D(data_pi2286)
    return {'data_pi2278': data_pi2278, 'data_pi2286': data_pi2286,
            'model_pi2278': rsp_pi2278(model),
            'model_pi2286': rsp_pi2286(model)}


def compare_results(expected, got, tol=1e-6):

    for key in ["succeeded", "numpoints", "dof"]:
        assert expected[key] == int(getattr(got, key))

    for key in ["istatval", "statval"]:
        val = float(getattr(got, key))
        assert val == pytest.approx(expected[key], rel=tol)

    # Convert the "got" values to a NumPy array to make sure we can
    # display the output nicely if there's a failure (since using a
    # tuple, which got.parvals is, can lead to less-than-useful
    # diagnostic errors when this test fails).
    #
    assert np.asarray(got.parvals) == pytest.approx(expected['parvals'],
                                                    rel=tol)


@requires_fits
@requires_group
@requires_xspec
@requires_data
def test_chi2xspecvar_stat(hide_logging, reset_xspec, setup_group):
    fit = Fit(setup_group['data'], setup_group['model'], Chi2XspecVar(), NelderMead())
    results = fit.fit()

    _fit_chi2xspecvar_results_bench = {
        'succeeded': 1,
        'numpoints': 143,
        'dof': 140,
        'istatval': 3386.3567560855163,
        'statval': 145.28136361814595,
        'parvals': np.array(
            [1.7780969509842395, 5.359332855712486, -1.9150879306302135]
            )
    }

    compare_results(_fit_chi2xspecvar_results_bench, results, tol=2e-5)


@requires_fits
@requires_group
@requires_xspec
@requires_data
def test_chi2modvar_stat(hide_logging, reset_xspec, setup_group):
    fit = Fit(setup_group['data'], setup_group['model'], Chi2ModVar(), NelderMead())
    results = fit.fit()

    _fit_chi2modvar_results_bench = {
        'succeeded': 1,
        'numpoints': 143,
        'dof': 140,
        'istatval': 87612.31084053952,
        'statval': 151.57751844595734,
        'parvals': np.array(
            [1.7348189936236909, 5.530307076686988, -2.0533493383033363]
            )
    }

    compare_results(_fit_chi2modvar_results_bench, results)


@requires_fits
@requires_group
@requires_xspec
@requires_data
def test_chi2constvar_stat(hide_logging, reset_xspec, setup_group):
    fit = Fit(setup_group['data'], setup_group['model'], Chi2ConstVar(), LevMar())
    results = fit.fit()

    _fit_chi2constvar_results_bench = {
        'succeeded': 1,
        'numpoints': 143,
        'dof': 140,
        'istatval': 3903.1647954751857,
        'statval': 140.1384389790626,
        'parvals': np.array(
            [1.8081424949164122, 5.461611041944607, -1.9077365276482876]
        )
    }

    compare_results(_fit_chi2constvar_results_bench, results, tol=2e-4)


@requires_fits
@requires_group
@requires_xspec
@requires_data
def test_chi2gehrels_stat(hide_logging, reset_xspec, setup_group):
    fit = Fit(setup_group['data'], setup_group['model'], Chi2Gehrels(), NelderMead())
    results = fit.fit()

    _fit_chi2gehrels_results_bench = {
        'succeeded': 1,
        'numpoints': 143,
        'dof': 140,
        'istatval': 2410.6724538404405,
        'statval': 100.22120218142433,
        'parvals': np.array(
            [1.7831198951601808, 5.376678702350025, -1.9152513828120412]
            )
    }

    # The tolerance used to be 1e-5 but for macOS ARM it has jumped to
    # 8e-5. Should the test values be re-evaluated for our supported
    # platforms?
    #
    compare_results(_fit_chi2gehrels_results_bench, results, tol=8e-5)


@requires_fits
@requires_group
@requires_xspec
@requires_data
def test_leastsq_stat(hide_logging, reset_xspec, setup_group):
    fit = Fit(setup_group['data'], setup_group['model'], LeastSq(), LevMar())
    results = fit.fit()

    _fit_leastsq_results_bench = {
        'succeeded': 1,
        'numpoints': 143,
        'dof': 140,
        'istatval': 117067.64900554597,
        'statval': 4203.173180288109,
        'parvals': np.array(
            [1.808142494916457, 5.461611041944977, -1.907736527635154]
        )
    }

    compare_results(_fit_leastsq_results_bench, results, tol=2e-4)


@requires_fits
@requires_xspec
@requires_data
def test_cstat_stat(hide_logging, reset_xspec, setup):
    fit = Fit(setup['data'], setup['model'], CStat(), NelderMead())
    results = fit.fit()

    _fit_cstat_results_bench = {
        'succeeded': 1,
        'numpoints': 460,
        'dof': 457,
        'istatval': 21647.62293983995,
        'statval': 472.6585691450068,
        'parvals': np.array(
            [1.75021021282262, 5.474614304244775, -1.9985761873334102]
        )
    }

    compare_results(_fit_cstat_results_bench, results, tol=1e-4)


@requires_fits
@requires_xspec
@requires_data
@pytest.mark.parametrize('stat', [Cash, MyCashWithBkg, MyCashNoBkg])
def test_cash_stat(stat, hide_logging, reset_xspec, setup):
    fit = Fit(setup['data'], setup['model'], stat(), NelderMead())
    results = fit.fit()

    _fit_mycash_results_bench = {
        'succeeded': 1,
        'numpoints': 460,
        'dof': 457,
        'istatval': 4594.094357753832,
        'statval': -16580.870007925594,
        'parvals': np.array(
            [1.7502119905435731, 5.474756852726033, -1.998450403564366]
        )
    }

    compare_results(_fit_mycash_results_bench, results, tol=1e-4)


@requires_fits
@requires_group
@requires_xspec
@requires_data
@pytest.mark.parametrize('stat', [MyChiWithBkg, MyChiNoBkg])
def test_mychi_data(stat, hide_logging, reset_xspec, setup_group):
    fit = Fit(setup_group['data'], setup_group['model'], stat(), LevMar())
    results = fit.fit()

    _fit_mychi_results_bench = {
        'succeeded': 1,
        'numpoints': 143,
        'dof': 140,
        'istatval': 117067.64900554594,
        'statval': 4211.349359724583,
        'parvals': np.array(
            [1.8177747886737923, 5.448440759203273, -1.8728780046411722]
        )
    }

    compare_results(_fit_mychi_results_bench, results, tol=2e-4)


@requires_fits
@requires_xspec
@requires_data
@pytest.mark.parametrize('stat', [MyCashNoBkg, MyCashWithBkg])
def test_mycash(stat, hide_logging, reset_xspec, setup_bkg):
    fit = Fit(setup_bkg['bkg'], setup_bkg['model'], stat(), NelderMead())
    results = fit.fit()

    _fit_mycashnobkg_results_bench = {
        'succeeded': 1,
        'numpoints': 446,
        'dof': 443,
        'istatval': 16800.126444958027,
        'statval': -995.672570978315,
        'parvals': np.array(
            # note: 'nh, gamma, norm' not 'gamma, norm, lognh'
            [0.10852102397785003, 2.163584294153897, 1.7059851092468952]
        )
    }

    compare_results(_fit_mycashnobkg_results_bench, results, tol=3e-4)


@requires_fits
@requires_group
@requires_xspec
@requires_data
@pytest.mark.parametrize('stat', [MyChiWithBkg, MyChiNoBkg])
def test_mychi_bkg(stat, hide_logging, reset_xspec, setup_bkg_group):
    fit = Fit(setup_bkg_group['bkg'], setup_bkg_group['model'], stat(), LevMar())
    results = fit.fit()

    _fit_mychinobkg_results_bench = {
        'succeeded': 1,
        'numpoints': 70,
        'dof': 67,
        'istatval': 12368.806484278228,
        'statval': 799.9399745311307,
        'parvals': np.array(
            [0.1904013138796835, 2.497496167887353, 2.111511871780941]
        )
    }

    compare_results(_fit_mychinobkg_results_bench, results, tol=1e-4)


@requires_fits
@requires_xspec
@requires_data
def test_wstat(hide_logging, reset_xspec, setup):

    fit = Fit(setup['data'], setup['model'], WStat(), LevMar())
    results = fit.fit()

    _fit_wstat_results_bench = {
        'succeeded': 1,
        'numpoints': 460,
        'dof': 457,
        'istatval': 21647.48285025895,
        'statval': 472.6585709918982,
        'parvals': np.array(
            [1.750204250228727, 5.47466040324842, -1.9983562007031974]
        )
    }

    compare_results(_fit_wstat_results_bench, results, tol=2e-4)


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
#     assert np.allclose(46.455049531, ui.calc_stat(), 1.e-7, 1.e-7)


@requires_fits
@requires_xspec
@requires_data
def test_wstat_error(hide_logging, reset_xspec, setup_bkg):
    fit = Fit(setup_bkg['bkg'], setup_bkg['model'], WStat(), NelderMead())

    with pytest.raises(StatErr):
        fit.fit()


def test_chi2datavar(hide_logging):
    num = 3
    xy = np.array(range(num))
    ui.load_arrays(1, xy, xy, Data1D)
    ui.set_stat('chi2datavar')
    err = ui.get_staterror()
    assert err == pytest.approx(np.sqrt(xy))


@requires_fits
@requires_xspec
@requires_data
def test_get_stat_info(hide_logging, make_data_path):
    fname_3c273 = make_data_path("3c273.pi")
    ui.load_pha(fname_3c273)
    src = ui.xspowerlaw.pl
    ui.set_source(src)
    ui.guess('pl')
    ui.set_stat('wstat')
    stat_info = ui.get_stat_info()[0]
    assert stat_info.dof == 44
    assert stat_info.numpoints == 46


@requires_fits
@requires_xspec
@requires_data
@pytest.mark.parametrize('stat', [MySimulStat, MyChiNoBkg])
def test_simul_stat_fit(stat, hide_logging, reset_xspec, setup_two):
    data1 = setup_two['data_pi2278']
    data2 = setup_two['data_pi2286']
    model1 = setup_two['model_pi2278']
    model2 = setup_two['model_pi2286']
    data = DataSimulFit(name='data1data2', datasets=[data1, data2])
    model = SimulFitModel(name='model1model2', parts=[model1, model2])
    fit = Fit(data=data, model=model, stat=stat(),
              method=NelderMead())
    result = fit.fit()

    _fit_simul_datavarstat_results_bench = {
        'succeeded': 1,
        'numpoints': 18,
        'dof': 15,
        'istatval': 56609.70689926489,
        'statval': 126.1509268988255,
        'parvals': np.array(
            [0.8417576197443695, 1.6496933246579941, 0.2383939869443424]
        )
    }

    # The tolerance used to be 1e-6 but for macOS ARM it has jumped to
    # 2.1e-4. Should the test values be re-evaluated for our supported
    # platforms?
    #
    compare_results(_fit_simul_datavarstat_results_bench,
                    result, tol=2.1e-4)


@requires_data
@requires_fits
def test_wstat_calc_stat_info(hide_logging, make_data_path, clean_astro_ui):
    "bug #147"
    ui.load_pha("stat", make_data_path("3c273.pi"))
    ui.set_source("stat", ui.powlaw1d.p1)
    ui.set_stat("wstat")
    ui.fit("stat")
    ui.get_stat_info()


@pytest.mark.parametrize("sexp,bexp,sscal,bscal,yexp,dyexp",
                         [(1, 1, 1, 1,
                           [0, -1, -3, 1, 3, 0, -2, 2, 0],
                           [1, 1, 1.73205078, 1, 1.73205078, 1.41421354, 2, 2, 2.44948983]),
                          (1, 10, 1, 1,
                           [0, -0.1, -0.3, 1, 3, 0.9, 0.7, 2.9, 2.7],
                           [0.1, 0.1, 0.173205078, 1, 1.73205078, 1.0049876, 1.01488912, 1.73493516, 1.74068952]),
                          (1, 0.1, 1.0, 1.0,
                           [0, -10, -30, 1, 3, -9, -29, -7, -27],
                           [1, 10, 17.320507, 1, 1.73205078, 10.0498753, 17.3493519, 10.1488914, 17.4068947]),
                          (100, 10, 0.1, 0.1,
                           [0, -10, -30, 1, 3, -9, -29, -7, -27],
                           [1, 10, 17.3205078, 1, 1.73205081, 10.0498758, 17.3493519, 10.1488918, 17.4068958]),
                          (100, 10, 0.1, 0.2,
                           [0, -5, -15, 1, 3, -4, -14, -2, -12],
                           [1, 5, 8.66025388, 1, 1.73205081, 5.09901941, 8.7177977, 5.29150255, 8.83176103]),
                          (100, 50, 0.1, 0.4,
                           [0, -0.5, -1.5, 1, 3, 0.5, -0.5, 2.5, 1.5],
                           [0.5, 0.5, 0.866025407, 1, 1.73205081, 1.11803403, 1.32287564, 1.80277564, 1.93649177]),
                          (100, 200, 0.1, 0.2,
                           [0, -0.25, -0.75, 1, 3, 0.75, 0.25, 2.75, 2.25],
                           [0.25, 0.25, 0.433012703, 1, 1.73205081, 1.03077637, 1.08972471, 1.75, 1.78535711])
                          ])
def test_xspecvar_zero_handling(sexp, bexp, sscal, bscal, yexp, dyexp):
    """How does XSPEC variance handle 0 in source and/or background?

    The values were calculated using XSPEC 12.10.1m (HEASOFT 6.26.1)
    using the following commands to create the file foo.dat which
    contains (after three 'header' lines) the data 'x 0.5 y dy'

        data foo.fits
        iplot data
        wplot foo.dat
        quit

    where foo.fits is a fake PHA file set up to have the channel/count
    values used below (a CSC-style PHA file was used so that source
    and background were in the same file but a separate bgnd PHA file
    could also have been used). Remember that XSPEC plots rates and
    not counts, so when the source exposure is not 1s, you need to
    account for this, as the y and dy values checked in this routine
    are in counts.

    """

    stat = Chi2XspecVar()
    chans = np.arange(1, 10, dtype=np.int16)
    scnts = np.asarray([0, 0, 0, 1, 3, 1, 1, 3, 3], dtype=np.int16)
    bcnts = np.asarray([0, 1, 3, 0, 0, 1, 3, 1, 3], dtype=np.int16)

    s = DataPHA('src', chans, scnts, exposure=sexp, backscal=sscal)
    b = DataPHA('bkg', chans, bcnts, exposure=bexp, backscal=bscal)
    s.set_background(b)
    s.subtract()

    y, dy, other = s.to_fit(staterrfunc=stat.calc_staterror)
    assert other is None
    assert y == pytest.approx(yexp)
    assert dy == pytest.approx(dyexp)


def test_xspecvar_zero_handling_variable():
    """How does XSPEC variance handle variable BACKSCAL?

    test_xspecvar_zero_handling checks most things, but not a
    variable scaling (that is, some combination of AREASCAL and
    BACKSCAL is not constant for each bin). In this case we handle
    a single case just to check.
    """

    stat = Chi2XspecVar()
    chans = np.arange(1, 5, dtype=np.int16)
    cnts = np.zeros(4, dtype=np.int16)

    sscal = np.asarray([0.1, 0.2, 2, 5])
    bscal = np.asarray([0.1, 0.1, 4, 20])
    s = DataPHA('src', chans, cnts, exposure=100, backscal=sscal)
    b = DataPHA('bkg', chans, cnts, exposure=125, backscal=bscal)
    s.set_background(b)
    s.subtract()

    yexp = np.asarray([0, 0, 0, 0])
    dyexp = np.asarray([0.8, 1.0, 0.4, 0.2])

    y, dy, other = s.to_fit(staterrfunc=stat.calc_staterror)
    assert other is None
    assert y == pytest.approx(yexp)
    assert dy == pytest.approx(dyexp)
