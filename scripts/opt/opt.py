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

from optparse import OptionParser

import numpy as np

from sherpa.optmethods import opt

from testing import Rosenbrock


parser = OptionParser()
parser.add_option("-N", "--npar", dest="npar", default=10,
                  type=int, help="set npar")
(options, args) = parser.parse_args()
npar = options.npar

x0 = np.array(npar * [-1.2, 1.0])
xmin = npar * [-1000, -1000]
xmax = npar * [1000, 1000]
factor = 10
seed = 234
simp = opt.SimplexNoStep(Rosenbrock, len(x0) + 1, x0, xmin, xmax, None,
                     seed, factor)
print('simp =\n', simp.simplex)

simp = opt.SimplexStep(Rosenbrock, len(x0) + 2, x0, xmin, xmax, x0 + 1.2,
                   seed, factor)
print('simp =\n', simp.simplex)

simp = opt.SimplexRandom(Rosenbrock, len(x0) + 5, x0, xmin, xmax, x0 + 1.2,
                     seed, factor)
print('simp =\n', simp.simplex)
