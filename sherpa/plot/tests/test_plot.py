#
#  Copyright (C) 2007, 2015, 2018, 2019, 2020, 2021
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

import numpy
import pytest

import sherpa.all as sherpa
from sherpa.ui.utils import Session as BaseSession
from sherpa.astro.ui.utils import Session as AstroSession
from sherpa.models import basic
from sherpa import plot as sherpaplot
from sherpa.data import Data1D, Data1DInt
from sherpa.utils.testing import requires_data, requires_plotting, requires_pylab


_datax = numpy.array(
    [0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
     11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.,  20.,  21.,
     22.,  23.,  24.,  25.,  26.,  27.,  28.,  29.,  30.,  31.,  32.,
     33.,  34.,  35.,  36.,  37.,  38.,  39.,  40.,  41.,  42.,  43.,
     44.,  45.,  46.,  47.,  48.,  49.,  50.,  51.,  52.,  53.,  54.,
     55.,  56.,  57.,  58.,  59.,  60.,  61.,  62.,  63.,  64.,  65.,
     66.,  67.,  68.,  69.,  70.,  71.,  72.,  73.,  74.,  75.,  76.,
     77.,  78.,  79.,  80.,  81.,  82.,  83.,  84.,  85.,  86.,  87.,
     88.,  89.,  90.,  91.,  92.,  93.,  94.,  95.,  96.,  97.,  98.,
     99.])


_datay = numpy.array(
    [0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     0.,   1.,   0.,   1.,   0.,   0.,   1.,   1.,   1.,   1.,   0.,
     3.,   5.,   7.,   6.,   7.,   5.,   5.,   5.,  15.,  18.,  15.,
     16.,   7.,  23.,  13.,  16.,  21.,  15.,  24.,  22.,  15.,  15.,
     8.,  10.,  14.,  10.,  10.,  11.,   0.,   6.,   6.,   8.,   8.,
     2.,   3.,   1.,   1.,   2.,   1.,   0.,   0.,   0.,   1.,   0.,
     2.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     0.])


# TODO: many tests in this class do not perform any assertions
# True, but I am glad they were there despite not performing any assertions
# as they made me spot some regressions I wouldn't have otherwise spotted (OL)
#

SetupPlot = namedtuple('SetupPlot',
                       ['data', 'g1', 'f'])
SetupContour = namedtuple('SetupContour',
                          ['data', 'g1', 'f', 'levels'])
SetupConfidence = namedtuple('SetupConfidence',
                             ['data', 'g1', 'f', 'ip', 'iu', 'rp', 'ru'])


@pytest.fixture
def setup_plot():

    data = sherpa.Data1D('testdata', _datax, _datay)
    g1 = sherpa.Gauss1D('g1')
    f = sherpa.Fit(data, g1)

    return SetupPlot(data, g1, f)


@pytest.fixture
def setup_contour(make_data_path):

    data = sherpa.read_data(make_data_path('gauss2d.dat'),
                            ncols=3, sep='\t', dstype=sherpa.Data2D)
    g1 = sherpa.Gauss2D('g1')
    g1.ellip.freeze()
    g1.theta.freeze()
    f = sherpa.Fit(data, g1)
    levels = numpy.array([0.5, 2, 5, 10, 20])

    return SetupContour(data, g1, f, levels)


@pytest.fixture
def setup_confidence():
    """Note that this performs a fit, so can be slow"""

    data = sherpa.Data1D('testdata', _datax, _datay)
    g1 = sherpa.Gauss1D('g1')
    f = sherpa.Fit(data, g1)
    f.fit()
    ip = sherpa.IntervalProjection()
    iu = sherpa.IntervalUncertainty()
    rp = sherpa.RegionProjection()
    ru = sherpa.RegionUncertainty()

    return SetupConfidence(data, g1, f, ip, iu, rp, ru)


# These tests don't need a backend as there is no "plot" call,
# and so do not need a @requires_plotting decorator.
#
def test_dataplot(setup_plot):
    dp = sherpa.DataPlot()
    dp.prepare(setup_plot.data, setup_plot.f.stat)
    # dp.plot()


def test_modelplot(setup_plot):
    mp = sherpa.ModelPlot()
    mp.prepare(setup_plot.data, setup_plot.g1, setup_plot.f.stat)
    # mp.plot()


def test_residplot(setup_plot):
    rp = sherpa.ResidPlot()
    rp.prepare(setup_plot.data, setup_plot.g1, setup_plot.f.stat)
    # rp.plot()


def test_delchiplot(setup_plot):
    dp = sherpa.DelchiPlot()
    dp.prepare(setup_plot.data, setup_plot.g1, setup_plot.f.stat)
    # dp.plot()


def test_chisqrplot(setup_plot):
    cs = sherpa.ChisqrPlot()
    cs.prepare(setup_plot.data, setup_plot.g1, setup_plot.f.stat)
    # cs.plot()


def test_ratioplot(setup_plot):
    tp = sherpa.RatioPlot()
    tp.prepare(setup_plot.data, setup_plot.g1, setup_plot.f.stat)
    # tp.plot()


@requires_plotting
@pytest.mark.parametrize("plottype", [sherpa.DelchiPlot,
                                      sherpa.RatioPlot,
                                      sherpa.ResidPlot])
def test_ignore_ylog_prefs_direct(setup_plot, plottype):
    """Do we ignore the ylog preference setting?"""
    tp = plottype()
    tp.plot_prefs['ylog'] = True
    tp.prepare(setup_plot.data, setup_plot.g1, setup_plot.f.stat)
    tp.plot()
    assert not tp.plot_prefs['ylog']


@requires_plotting
@pytest.mark.parametrize("plottype", [sherpa.DelchiPlot,
                                      sherpa.RatioPlot,
                                      sherpa.ResidPlot])
def test_ignore_ylog_prefs(setup_plot, plottype):
    """Do we ignore the ylog preference setting?

    Use the get_prefs accessor.
    """
    tp = plottype()
    prefs = tp.get_prefs()
    prefs['ylog'] = True
    tp.prepare(setup_plot.data, setup_plot.g1, setup_plot.f.stat)
    tp.plot()
    assert not tp.get_prefs()['ylog']


def test_fitplot(setup_plot):
    dp = sherpa.DataPlot()
    dp.prepare(setup_plot.data, setup_plot.f.stat)

    mp = sherpa.ModelPlot()
    mp.prepare(setup_plot.data, setup_plot.g1, setup_plot.f.stat)

    fp = sherpa.FitPlot()
    fp.prepare(dp, mp)
    # fp.plot()


def test_splitplot(setup_plot):
    dp = sherpa.DataPlot()
    dp.prepare(setup_plot.data, setup_plot.f.stat)

    mp = sherpa.ModelPlot()
    mp.prepare(setup_plot.data, setup_plot.g1, setup_plot.f.stat)

    rp = sherpa.ResidPlot()
    rp.prepare(setup_plot.data, setup_plot.g1, setup_plot.f.stat)

    fp = sherpa.FitPlot()
    fp.prepare(dp, mp)

    # TODO: without calling addplot this is hust the same as test_splitcontour
    sp = sherpa.SplitPlot(2, 2)
    # sp.addplot(dp)
    # sp.addplot(mp)
    # sp.addplot(fp)
    # sp.addplot(rp)


@requires_data
def test_datacontour(setup_contour):
    dc = sherpa.DataContour()
    dc.prepare(setup_contour.data)
    dc.levels = setup_contour.levels
    # dc.contour()


@requires_data
def test_modelcontour(setup_contour):
    mc = sherpa.ModelContour()
    mc.prepare(setup_contour.data, setup_contour.g1, setup_contour.f.stat)
    mc.levels = setup_contour.levels
    # mc.contour()


@requires_data
def test_residcontour(setup_contour):
    rc = sherpa.ResidContour()
    rc.prepare(setup_contour.data, setup_contour.g1, setup_contour.f.stat)
    rc.levels = setup_contour.levels
    # rc.contour()


@requires_data
def test_ratiocontour(setup_contour):
    tc = sherpa.RatioContour()
    tc.prepare(setup_contour.data, setup_contour.g1, setup_contour.f.stat)
    tc.levels = setup_contour.levels
    # tc.contour()


@requires_data
def test_fitcontour(setup_contour):
    dc = sherpa.DataContour()
    dc.prepare(setup_contour.data)

    mc = sherpa.ModelContour()
    mc.prepare(setup_contour.data, setup_contour.g1, setup_contour.f.stat)

    fc = sherpa.FitContour()
    fc.prepare(dc, mc)
    # fc.contour()


@requires_data
def test_splitcontour(setup_contour):
    dc = sherpa.DataContour()
    dc.levels = setup_contour.levels
    dc.prepare(setup_contour.data)

    mc = sherpa.ModelContour()
    mc.levels = setup_contour.levels
    mc.prepare(setup_contour.data, setup_contour.g1, setup_contour.f.stat)

    fc = sherpa.FitContour()
    fc.prepare(dc, mc)

    rc = sherpa.ResidContour()
    rc.prepare(setup_contour.data, setup_contour.g1, setup_contour.f.stat)
    rc.levels = setup_contour.levels

    # TODO: without calling addcontour this is hust the same as test_splitplot
    sp = sherpa.SplitPlot(2, 2)
    # sp.addcontour(dc)
    # sp.addcontour(mc)
    # sp.addcontour(fc)
    # sp.addcontour(rc)


def test_interval_projection(setup_confidence):
    _ipx = numpy.array(
        [15.60720526,  15.92784424,  16.24848322,  16.56912221,
         16.88976119,  17.21040017,  17.53103916,  17.85167814,
         18.17231712,  18.49295611,  18.81359509,  19.13423407,
         19.45487306,  19.77551204,  20.09615102,  20.41679001,
         20.73742899,  21.05806798,  21.37870696,  21.69934594])
    _ipy = numpy.array(
        [40.09661435,  39.18194614,  38.37963283,  37.68723716,
         37.10218543,  36.62179549,  36.2432939,  35.96383078,
         35.78049297,  35.69031588,  35.69029438,  35.77739294,
         35.94855493,  36.20071145,  36.53078937,  36.9357187,
         37.41243938,  37.95790737,  38.56910025,  39.24302218])

    setup_confidence.ip.fac = 2
    setup_confidence.ip.calc(setup_confidence.f,
                             setup_confidence.g1.fwhm)

    assert setup_confidence.ip.x == pytest.approx(_ipx, abs=1e-4)
    assert setup_confidence.ip.y == pytest.approx(_ipy, abs=1e-4)
    # setup_confidence.ip.plot()


def test_interval_uncertainty(setup_confidence):
    _iux = numpy.array(
        [15.60720526,  15.92784424,  16.24848322,  16.56912221,
         16.88976119,  17.21040017,  17.53103916,  17.85167814,
         18.17231712,  18.49295611,  18.81359509,  19.13423407,
         19.45487306,  19.77551204,  20.09615102,  20.41679001,
         20.73742899,  21.05806798,  21.37870696,  21.69934594])

    _iuy = numpy.array(
        [42.2582845,  40.96299483,  39.80577195,  38.78821363,
         37.91185112,  37.17815809,  36.58855785,  36.14442877,
         35.84710827,  35.69789528,  35.69805141,  35.84880111,
         36.15133085,  36.60678759,  37.21627669,  37.98085933,
         38.90154977,  39.97931233,  41.21505839,  42.60964342])

    setup_confidence.iu.fac = 2
    setup_confidence.iu.calc(setup_confidence.f, setup_confidence.g1.fwhm)
    assert setup_confidence.iu.x == pytest.approx(_iux, abs=1e-4)
    assert setup_confidence.iu.y == pytest.approx(_iuy, abs=1e-4)
    # setup_confidence.iu.plot()


def test_region_projection(setup_confidence):
    _rpx0 = numpy.array(
        [11.03809974,  12.73036104,  14.42262235,  16.11488365, 17.80714495,
         19.49940625,  21.19166755,  22.88392885,  24.57619016, 26.26845146,
         11.03809974,  12.73036104,  14.42262235,  16.11488365, 17.80714495,
         19.49940625,  21.19166755,  22.88392885,  24.57619016, 26.26845146,
         11.03809974,  12.73036104,  14.42262235,  16.11488365, 17.80714495,
         19.49940625,  21.19166755,  22.88392885,  24.57619016, 26.26845146,
         11.03809974,  12.73036104,  14.42262235,  16.11488365, 17.80714495,
         19.49940625,  21.19166755,  22.88392885,  24.57619016, 26.26845146,
         11.03809974,  12.73036104,  14.42262235,  16.11488365, 17.80714495,
         19.49940625,  21.19166755,  22.88392885,  24.57619016, 26.26845146,
         11.03809974,  12.73036104,  14.42262235,  16.11488365, 17.80714495,
         19.49940625,  21.19166755,  22.88392885,  24.57619016, 26.26845146,
         11.03809974,  12.73036104,  14.42262235,  16.11488365, 17.80714495,
         19.49940625,  21.19166755,  22.88392885,  24.57619016, 26.26845146,
         11.03809974,  12.73036104,  14.42262235,  16.11488365, 17.80714495,
         19.49940625,  21.19166755,  22.88392885,  24.57619016, 26.26845146,
         11.03809974,  12.73036104,  14.42262235,  16.11488365, 17.80714495,
         19.49940625,  21.19166755,  22.88392885,  24.57619016, 26.26845146,
         11.03809974,  12.73036104,  14.42262235,  16.11488365, 17.80714495,
         19.49940625,  21.19166755,  22.88392885,  24.57619016, 26.26845146])

    _rpx1 = numpy.array(
        [8.75749218,   8.75749218,   8.75749218,   8.75749218, 8.75749218,
         8.75749218,   8.75749218,   8.75749218,   8.75749218, 8.75749218,
         10.54556708,  10.54556708,  10.54556708,  10.54556708, 10.54556708,
         10.54556708,  10.54556708,  10.54556708,  10.54556708, 10.54556708,
         12.33364199,  12.33364199,  12.33364199,  12.33364199, 12.33364199,
         12.33364199,  12.33364199,  12.33364199,  12.33364199, 12.33364199,
         14.12171689,  14.12171689,  14.12171689,  14.12171689, 14.12171689,
         14.12171689,  14.12171689,  14.12171689,  14.12171689, 14.12171689,
         15.9097918,  15.9097918,  15.9097918,  15.9097918, 15.9097918,
         15.9097918,  15.9097918,  15.9097918,  15.9097918, 15.9097918,
         17.6978667,  17.6978667,  17.6978667,  17.6978667, 17.6978667,
         17.6978667,  17.6978667,  17.6978667,  17.6978667, 17.6978667,
         19.48594161,  19.48594161,  19.48594161,  19.48594161, 19.48594161,
         19.48594161,  19.48594161,  19.48594161,  19.48594161, 19.48594161,
         21.27401651,  21.27401651,  21.27401651,  21.27401651, 21.27401651,
         21.27401651,  21.27401651,  21.27401651,  21.27401651, 21.27401651,
         23.06209142,  23.06209142,  23.06209142,  23.06209142, 23.06209142,
         23.06209142,  23.06209142,  23.06209142,  23.06209142, 23.06209142,
         24.85016632,  24.85016632,  24.85016632,  24.85016632, 24.85016632,
         24.85016632,  24.85016632,  24.85016632,  24.85016632, 24.85016632])

    _rpy = numpy.array(
        [121.18227727,  109.21389788,   98.48751045,   89.15211989,
         81.31437819,   75.0489422,   70.40478381,   67.40903363,
         66.06931267,   66.37505662,  107.09555406,   93.85962337,
         82.2173972,   72.37069659,   64.46288283,   58.59635882,
         54.8420303,   53.24411902,   53.82255037,   56.57387709,
         95.08597997,   80.98176715,   68.85606392,   58.97098329,
         51.51170797,   46.61282359,   44.37139326,   44.85269567,
         48.09257852,   54.09778178,   85.15342677,   70.58027832,
         58.4035077,   48.95296874,   42.46079117,   39.0982177,
         38.99272007,   42.23460729,   48.87925026,   58.9466444,
         77.29778983,   62.65512325,   50.85972221,   42.31664458,
         37.31008905,   36.05245858,   38.70590422,   45.38974383,
         56.18246641,   71.12038263,   71.51901508,   57.20627498,
         46.22471072,   39.06200464,   36.0595698,   37.4754866,
         43.51086985,   54.31802734,   70.00215676,   90.6189382,
         67.81704451,   54.23371422,   44.49847032,   39.18904432,
         38.70920964,   43.36725746,   53.40756108,   69.01940099,
         90.33827068,  117.44226953,   66.19183538,   53.73742613,
         45.68100001,   42.6977601,   45.25899046,   53.72773755,
         68.39593576,   89.49382219,  117.19077045,  151.59034586,
         66.64335816,   55.71739907,   49.77229917,   49.58814921,
         55.70889825,   68.55690091,   88.47596148,  115.74125836,
         150.55962734,  193.06314386,   69.17158712,   60.17362372,
         56.77236574,   59.86020953,   70.05892198,   87.8547272,
         113.64761294,  147.76168412,  190.44481909,  241.86064553])

    setup_confidence.rp.fac = 5
    setup_confidence.rp.calc(setup_confidence.f,
                             setup_confidence.g1.fwhm,
                             setup_confidence.g1.ampl)
    assert setup_confidence.rp.x0 == pytest.approx(_rpx0, abs=1e-4)
    assert setup_confidence.rp.x1 == pytest.approx(_rpx1, abs=1e-4)
    assert setup_confidence.rp.y == pytest.approx(_rpy, abs=1e-4)
    # setup_confidence.rp.contour()


def test_region_uncertainty(setup_confidence):
    _rux0 = numpy.array(
        [12.56113491,  13.91494395,  15.268753,  16.62256204, 17.97637108,
         19.33018012,  20.68398916,  22.0377982,  23.39160724, 24.74541629,
         12.56113491,  13.91494395,  15.268753,  16.62256204, 17.97637108,
         19.33018012,  20.68398916,  22.0377982,  23.39160724, 24.74541629,
         12.56113491,  13.91494395,  15.268753,  16.62256204, 17.97637108,
         19.33018012,  20.68398916,  22.0377982,  23.39160724, 24.74541629,
         12.56113491,  13.91494395,  15.268753,  16.62256204, 17.97637108,
         19.33018012,  20.68398916,  22.0377982,  23.39160724, 24.74541629,
         12.56113491,  13.91494395,  15.268753,  16.62256204, 17.97637108,
         19.33018012,  20.68398916,  22.0377982,  23.39160724, 24.74541629,
         12.56113491,  13.91494395,  15.268753,  16.62256204, 17.97637108,
         19.33018012,  20.68398916,  22.0377982,  23.39160724, 24.74541629,
         12.56113491,  13.91494395,  15.268753,  16.62256204, 17.97637108,
         19.33018012,  20.68398916,  22.0377982,  23.39160724, 24.74541629,
         12.56113491,  13.91494395,  15.268753,  16.62256204, 17.97637108,
         19.33018012,  20.68398916,  22.0377982,  23.39160724, 24.74541629,
         12.56113491,  13.91494395,  15.268753,  16.62256204, 17.97637108,
         19.33018012,  20.68398916,  22.0377982,  23.39160724, 24.74541629,
         12.56113491,  13.91494395,  15.268753,  16.62256204, 17.97637108,
         19.33018012,  20.68398916,  22.0377982,  23.39160724, 24.74541629])

    _rux1 = numpy.array(
        [10.36675959,  10.36675959,  10.36675959,  10.36675959, 10.36675959,
         10.36675959,  10.36675959,  10.36675959,  10.36675959, 10.36675959,
         11.79721952,  11.79721952,  11.79721952,  11.79721952, 11.79721952,
         11.79721952,  11.79721952,  11.79721952,  11.79721952, 11.79721952,
         13.22767944,  13.22767944,  13.22767944,  13.22767944, 13.22767944,
         13.22767944,  13.22767944,  13.22767944,  13.22767944, 13.22767944,
         14.65813936,  14.65813936,  14.65813936,  14.65813936, 14.65813936,
         14.65813936,  14.65813936,  14.65813936,  14.65813936, 14.65813936,
         16.08859929,  16.08859929,  16.08859929,  16.08859929, 16.08859929,
         16.08859929,  16.08859929,  16.08859929,  16.08859929, 16.08859929,
         17.51905921,  17.51905921,  17.51905921,  17.51905921, 17.51905921,
         17.51905921,  17.51905921,  17.51905921,  17.51905921, 17.51905921,
         18.94951914,  18.94951914,  18.94951914,  18.94951914, 18.94951914,
         18.94951914,  18.94951914,  18.94951914,  18.94951914, 18.94951914,
         20.37997906,  20.37997906,  20.37997906,  20.37997906, 20.37997906,
         20.37997906,  20.37997906,  20.37997906,  20.37997906, 20.37997906,
         21.81043898,  21.81043898,  21.81043898,  21.81043898, 21.81043898,
         21.81043898,  21.81043898,  21.81043898,  21.81043898, 21.81043898,
         23.24089891,  23.24089891,  23.24089891,  23.24089891, 23.24089891,
         23.24089891,  23.24089891,  23.24089891,  23.24089891, 23.24089891])

    _ruy = numpy.array(
        [96.58117835,   87.02838072,   78.58324208,   71.31800953,
         65.28984102,   60.54431967,   57.11719555,   55.03459229,
         54.31255884,   54.95668786,   85.96169801,   75.9815896,
         67.33096852,   60.09842943,   54.35398817,   50.15439725,
         47.54566158,   46.56317725,   47.23082527,   49.56009821,
         76.90300053,   66.7116158,   58.0882768,   51.13947262,
         45.94928247,   42.58681159,   41.10961032,   41.56372169,
         43.98222109,   48.38374461,   69.40508592,   59.21845932,
         50.85516692,   44.4411391,   40.07572392,   37.84156267,
         37.80904175,   40.0362256,   44.56674629,   51.42762706,
         63.46795418,   53.50212017,   45.63163889,   40.00342886,
         36.73331252,   35.91865051,   37.6439559,   41.980689,
         48.98440089,   58.69174556,   59.09160531,   49.56259833,
         42.41769271,   37.82634191,   35.92204826,   36.8180751,
         40.61435275,   47.39711188,   57.23518487,   70.1761001,
         56.27603931,   47.39989381,   41.21332836,   37.90987826,
         37.64193114,   40.53983645,   46.7202323,   56.28549424,
         69.31909824,   85.88069069,   55.02125618,   47.01400662,
         42.01854587,   40.25403789,   41.89296118,   47.08393454,
         55.96159455,   68.64583609,   85.236141,  105.80551733,
         55.32725591,   48.40493674,   44.83334521,   44.85882081,
         48.67513836,   56.45036939,   68.33843951,   84.47813741,
         104.98631314,  129.95058002,   57.19403852,   51.57268419,
         49.6577264,   51.72422701,   57.98846268,   68.63914099,
         83.85076718,  103.78239821,  128.56961467,  158.31587876])

    setup_confidence.ru.fac = 4
    setup_confidence.ru.calc(setup_confidence.f,
                             setup_confidence.g1.fwhm,
                             setup_confidence.g1.ampl)
    assert setup_confidence.ru.x0 == pytest.approx(_rux0, abs=1e-4)
    assert setup_confidence.ru.x1 == pytest.approx(_rux1, abs=1e-4)
    assert setup_confidence.ru.y == pytest.approx(_ruy, abs=1e-4)
    # setup_confidence.ru.contour()


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_source_component_arbitrary_grid(session):
    ui = session()

    # depending on how the test is used we need to have
    # registered a model with the Session.
    #
    ui._add_model_types(basic)

    model = ui.create_model_component('const1d', 'c')
    model.c0 = 10

    def tst(x, y, re_x, yy):
        ui.load_arrays(1, x, y)
        regrid_model = model.regrid(re_x)

        # could set the model to regrid to show can send in any source model
        ui.set_source(regrid_model)

        ui.plot_source_component(regrid_model)
        assert (ui._compsrcplot.x == x).all()
        assert ui._compsrcplot.y == pytest.approx(yy)

    x = [1, 2, 3]
    y = [1, 2, 3]
    re_x = [10, 20, 30]
    tst(x, y, re_x, [0, 0, 0])

    x = numpy.linspace(1, 10, 10)
    y = x
    re_x = numpy.linspace(5, 15, 15)
    tst(x, y, re_x, [0, 0, 0, 0, 10, 10, 10, 10, 10, 10])

    re_x = numpy.linspace(1, 5, 15)
    tst(x, y, re_x, [10, 10, 10, 10, 10,  0,  0,  0,  0,  0])

    re_x = numpy.linspace(3, 5, 15)
    tst(x, y, re_x, [0,  0, 10, 10, 10,  0,  0,  0,  0,  0])


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_model_arbitrary_grid_integrated(session):
    ui = session()
    model = basic.Const1D('c')
    model.c0 = 10

    def tst(x, y, re_x, yy):
        ui.load_arrays(1, x[0], x[1], y, Data1DInt)
        regrid_model = model.regrid(*re_x)
        ui.set_model(regrid_model)
        ui.plot_model()

        # should this use ui.get_model_plot()?
        assert ui._modelhistplot.xlo == pytest.approx(x[0])
        assert ui._modelhistplot.xhi == pytest.approx(x[1])
        assert ui._modelhistplot.y == pytest.approx(yy)

    tmp = numpy.arange(1, 5, 1)
    x = tmp[:-1], tmp[1:]
    y = x[0]
    tmp = numpy.arange(10, 50, 10)
    re_x = tmp[:-1], tmp[1:]
    tst(x, y, re_x, [0, 0, 0])

    tmp = numpy.arange(1, 20, 1)
    x = tmp[:-1], tmp[1:]
    y = x[0]
    tmp = numpy.arange(1, 20, 0.5)
    re_x = tmp[:-1], tmp[1:]
    tst(x, y, re_x, len(y) * [10])

    tmp = numpy.arange(1, 20, 1)
    x = tmp[:-1], tmp[1:]
    y = x[0]
    tmp = numpy.arange(10, 20, 0.5)
    re_x = tmp[:-1], tmp[1:]
    n = int(len(y) / 2)
    yy = numpy.append(n * [0.], n * [10.])
    tst(x, y, re_x, yy)


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_source_component_arbitrary_grid_int(session):
    ui = session()

    x = numpy.array([1, 2, 3]), numpy.array([2, 3, 4])
    y = [1.5, 2.5, 3.5]
    re_x = numpy.array([10, 20, 30]), numpy.array([20, 30, 40])

    ui.load_arrays(1, x[0], x[1], y, Data1DInt)
    model = basic.Const1D('c')
    model.c0 = 10

    regrid_model = model.regrid(*re_x)
    ui.plot_source_component(regrid_model)

    # should this use ui.get_source_component_plot()?
    assert ui._compsrchistplot.xlo == pytest.approx(x[0])
    assert ui._compsrchistplot.xhi == pytest.approx(x[1])
    assert ui._compsrchistplot.y == pytest.approx([0.0, 0.0, 0.0])


def test_numpy_histogram_density_vs_normed():
    from sherpa.astro import ui

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3])
    ui.set_source('const1d.c')
    c = ui.get_model_component('c')
    ui.fit()
    res = ui.eqwidth(c, c+c, error=True)
    ui.plot_pdf(res[4])
    plot = ui.get_pdf_plot()
    expected_x = numpy.linspace(2.5, 3.5, 13)
    expected_xlo, expected_xhi = expected_x[:-1], expected_x[1:]
    expected_y = [0, 0, 0, 0, 0, 0, 12, 0, 0, 0, 0, 0]
    assert plot.y == pytest.approx(expected_y)
    assert plot.xlo == pytest.approx(expected_xlo)
    assert plot.xhi == pytest.approx(expected_xhi)


def test_errors_with_no_stat():
    """Check we get no errors when stat is None"""

    d = Data1D('x', numpy.asarray([2, 4, 10]), numpy.asarray([2, 4, 0]))
    dp = sherpaplot.DataPlot()
    dp.prepare(d, stat=None)
    assert dp.yerr is None


def test_histogram_returns_x():
    """We support x accessor for histogram plots."""

    xlo = [10, 20, 40, 80]
    xhi = [20, 40, 60, 90]
    y = [1, 2, 3, 4]
    d = Data1DInt('xx', xlo, xhi, y)

    dp = sherpaplot.DataHistogramPlot()
    dp.prepare(d)

    xlo = numpy.asarray(xlo)
    xhi = numpy.asarray(xhi)

    assert dp.xlo == pytest.approx(xlo)
    assert dp.xhi == pytest.approx(xhi)
    assert dp.x == pytest.approx(0.5 * (xlo + xhi))

    # check x is not in the str output
    #
    for line in str(dp).split('\n'):
        toks = line.split()
        assert len(toks) > 2
        assert toks[0] != 'x'


def test_histogram_can_not_set_x():
    """We cannot change the x accessor"""

    xlo = [10, 20, 40, 80]
    xhi = [20, 40, 60, 90]
    y = [1, 2, 3, 4]
    d = Data1DInt('xx', xlo, xhi, y)

    dp = sherpaplot.DataHistogramPlot()
    dp.prepare(d)

    with pytest.raises(AttributeError):
        dp.x = xlo


def test_histogram_empty_x():
    """x accessor is None if there's no data"""

    dp = sherpaplot.DataHistogramPlot()
    assert dp.x is None


@requires_plotting
@pytest.mark.parametrize("pclass", [sherpaplot.DataPlot,
                                    sherpaplot.DataHistogramPlot,
                                    sherpaplot.DataContour,
                                    sherpaplot.Confidence2D])
def test_get_prefs_xlog(pclass):
    """Does get_prefs return a dictionary?

    I don't think this requires a specific backend.
    """

    prefs = pclass().get_prefs()
    assert not prefs['xlog']


@requires_pylab
@pytest.mark.parametrize("pclass,setting,value",
                         [(sherpaplot.Point, 'symbol', None),
                          (sherpaplot.SplitPlot, 'wspace', 0.3)])
def test_get_prefs_setting(pclass, setting, value):
    """Does get_prefs return a dictionary?

    The preferences depend on the specific backend.
    """

    prefs = pclass().get_prefs()
    assert prefs[setting] == value
