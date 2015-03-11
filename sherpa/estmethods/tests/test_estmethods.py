# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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
from sherpa.estmethods import *
from sherpa.utils import SherpaTestCase

# Test data arrays -- together this makes a line best fit with a
# 1D Gaussian function.
x = numpy.array(
    [ 758,  759,  760,  761,  762,  763,  764,  765,  766,  767,  768,  769,
      770,  771,  772,  773,  774,  775,  776,  777,  778,  779,  780,  781,
      782,  783,  784,  785,  786,  787,  788,  789,  790,  791,  792,  793,
      794,  795,  796,  797,  798,  799,  800,  801,  802,  803,  804,  805,
      806,  807,  808,  809,  810,  811,  812,  813,  814,  815,  816,  817,
      818,  819,  820,  821,  822,  823,  824,  825,  826,  827,  828,  829,
      830,  831,  832,  833,  834,  835,  836,  837,  838,  839,  840,  841,
      842,  843,  844,  845,  846,  847,  848,  849,  850,  851,  852,  853,
      854,  855,  856,  857,  858,  859,  860,  861,  862,  863,  864,  865,
      866,  867,  868,  869,  870,  871,  872,  873,  874,  875,  876,  877,
      878,  879,  880,  881])

y = numpy.array(
    [ 1,    0,    1,    0,    1,    2,    0,    3,    0,    1,    0,    5,
      0,    5,    3,    5,    6,    5,   11,   14,   11,   13,   12,   21,
      15,   24,   20,   29,   32,   43,   47,   49,   50,   64,   60,   72,
      61,   73,   83,   98,   99,  100,   94,   92,  121,  107,  126,  107,
      112,  123,  114,  126,  113,   86,  111,  126,   95,  119,   93,  119,
      93,   89,   75,   80,   71,   68,   59,   54,   61,   37,   21,   33,
      37,   32,   31,   22,   19,   25,   14,   13,   12,   10,    7,   10,
      5,    4,    8,    1,    5,    2,    1,    3,    1,    5,    0,    3,
      1,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0,    0,
      0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
      0,    0,    0,    1])

# These parameter values are known to be the best-fit parameter values
# of a Gaussian to the x- and y-arrays above.
fittedpars = numpy.array([32.3536, 807.863, 117.826])
limit_parnums = numpy.array([0,1,2])
maxpars = numpy.array([200, 900, 200])
minpars = numpy.array([1, 0, 0])
hardmaxpars = numpy.array([1.0e+120, 1.0e+120, 1.0e+120])
hardminpars = numpy.array([1.0e-120, -1.0e+120, -1.0e+120])

gfactor = 4.0 * 0.6931471805599453094172321214581766

def freeze_par(pars, parmins, parmaxes, i):
    return (pars, parmins, parmaxes)

def thaw_par(i):
    pass

def report_progress(i, lower, upper):
    pass

def get_par_name( ii ):
    pass

# Here's the 1D Gaussian function to use to generate predicted
# data, which will be compared to the y-array above
def gauss_func(p):
    return p[2] * numpy.exp(-1.0 * gfactor *
                            (x - p[1]) * (x - p[1]) / p[0] / p[0])

# Compare the y-array above to values calculated with gauss_func,
# and calculate chi-squared, using Gehrel's approximation for
# estimating errors.
def stat(p):
    errors = 1.0 + numpy.sqrt(y + 0.75)
    fvec = (y - gauss_func(p) ) / errors
    return ((fvec * fvec).sum(),)

# Easiest "fit" function for unit tests is actually just to
# return current parameter values.  In this test, that will
# have the effect of making projection act like the old
# uncertainty method, but that is OK for the unit test.
def fitter(scb, pars, parmins, parmaxs):
    return (1, pars, scb(pars)[0])

class test_estmethods(SherpaTestCase):

    def test_covar_failures(self):
        self.assertRaises(TypeError, Covariance().compute,
                          None, fitter, fittedpars,
                          minpars, maxpars,
                          hardminpars, hardmaxpars,
                          limit_parnums, freeze_par, thaw_par, report_progress)
        
        self.assertRaises(RuntimeError, Covariance().compute,
                          stat, fitter, None, minpars, maxpars,
                          hardminpars, hardmaxpars, limit_parnums, freeze_par,
                          thaw_par, report_progress, get_par_name)
        
        self.assertRaises(RuntimeError, Covariance().compute,
                          stat, fitter, numpy.array([1,2]),
                          minpars, maxpars, hardminpars, hardmaxpars,
                          limit_parnums, freeze_par, thaw_par,
                          report_progress, get_par_name)

        self.assertRaises(RuntimeError, Covariance().compute,
                          stat, fitter, fittedpars, numpy.array([1,2]),
                          maxpars, hardminpars, hardmaxpars, limit_parnums,
                          freeze_par, thaw_par, report_progress, get_par_name)
        
    def test_projection_failures(self):
        self.assertRaises(TypeError, Projection().compute,
                          stat, None, fittedpars, minpars, maxpars,
                          hardminpars, hardmaxpars, limit_parnums, freeze_par,
                          thaw_par, report_progress, get_par_name)

        self.assertRaises(RuntimeError, Projection().compute,
                          stat, fitter, None, minpars, maxpars,
                          hardminpars, hardmaxpars, limit_parnums,
                          freeze_par, thaw_par, report_progress, get_par_name)

        self.assertRaises(RuntimeError, Projection().compute,
                          stat, fitter, numpy.array([1,2]),
                          minpars, maxpars, hardminpars, hardmaxpars,
                          limit_parnums, freeze_par, thaw_par,
                          report_progress, get_par_name)

        self.assertRaises(RuntimeError, Projection().compute,
                          stat, fitter, fittedpars, numpy.array([1,2]),
                          maxpars, hardminpars, hardmaxpars, limit_parnums,
                          freeze_par, thaw_par, report_progress, get_par_name)
        
    def test_covar(self):
        standard = numpy.array([[ 0.4935702,  0.06857833, numpy.nan],
                                [ 0.06857833, 0.26405554, numpy.nan],
                                [ numpy.nan,  numpy.nan,  2.58857314]])
        results = Covariance().compute(stat, None, fittedpars,
                                       minpars, maxpars,
                                       hardminpars, hardmaxpars,
                                       limit_parnums, freeze_par, thaw_par,
                                       report_progress, get_par_name)
        self.assertEqualWithinTol(standard.diagonal(),
                                  #results[2].diagonal(), 1e-4)
                                  results[1], 1e-4)

    def test_projection(self):
        standard_elo = numpy.array([-0.39973743, -0.26390339, -2.08784716])
        standard_ehi = numpy.array([ 0.39580942,  0.26363223,  2.08789851])
        results = Projection().compute(stat, fitter, fittedpars,
                                       minpars, maxpars,
                                       hardminpars, hardmaxpars,
                                       limit_parnums, freeze_par, thaw_par,
                                       report_progress, get_par_name)
        self.assertEqualWithinTol(standard_elo,results[0], 1e-4)
        self.assertEqualWithinTol(standard_ehi,results[1], 1e-4)
