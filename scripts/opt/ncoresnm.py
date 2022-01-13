#!/usr/bin/env python

#
#  Copyright (C) 2019, 2020, 2021, 2022
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

"""
This used to be in sherpa.optmethods.ncoresnm

"""

from optparse import OptionParser

from sherpa.optmethods import ncoresnm

from testing import tst_opt, tst_unc_opt


parser = OptionParser()
parser.add_option("-N", "--num", dest="num", default=10,
                  type=int, help="set num")
parser.add_option("-s", "--single", dest="single", default=True,
                  action="store_false", help="run test using 1 core")
parser.add_option("-u", "--unc_opt", dest="unc_opt", default=True,
                  action="store_false", help="do not run tst_unc_opt")
parser.add_option("-o", "--opt", dest="global_func", default=True,
                  action="store_false", help="do not run tst_opt")
(options, args) = parser.parse_args()
npar = options.num

if npar % 2 != 0:
    raise ValueError("-N option must be an even number")


if options.single:
    # midnight = Midnight()
    algorithms = [ncoresnm.ncoresNelderMead()]
else:
    algorithms = [ncoresnm.NelderMead0(), ncoresnm.NelderMead1(), ncoresnm.NelderMead2(),
                  ncoresnm.NelderMead3(), ncoresnm.NelderMead4(), ncoresnm.NelderMead5(),
                  ncoresnm.NelderMead6(), ncoresnm.NelderMead7()]

if options.unc_opt:
    tst_unc_opt(algorithms, npar)

if options.global_func:
    tst_opt(algorithms, npar)
