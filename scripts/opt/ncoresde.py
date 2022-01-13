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
This used to be in sherpa.optmethods.ncoresde

"""

import argparse

from sherpa.optmethods import ncoresde

from testing import tst_opt, tst_unc_opt


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--difevo', action="store_true",
                    default=False, help='run simple difevo', dest="difevo")
parser.add_argument("-u", "--unc_opt", dest="unc_opt", default=True,
                    action="store_false", help="do not run tst_unc_opt")
parser.add_argument("-o", "--opt", dest="global_func", default=True,
                    action="store_false", help="do not run tst_opt")
parser.add_argument('-N', action="store", dest="num", default=4, type=int)

options = parser.parse_args()
# print('options =', options)
npar = options.num
if npar % 2 != 0:
    raise ValueError("-N option must be an even number")

if options.difevo:
    algo = [ncoresde.DifEvo()]
else:
    algo = [ncoresde.ncoresDifEvo()]

if options.unc_opt:
    tst_unc_opt(algo, npar)

if options.global_func:
    tst_opt(algo, npar)
