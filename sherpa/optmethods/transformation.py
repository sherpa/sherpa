#
#  Copyright (C) 2022
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

import numpy

from sherpa.models.parameter import hugeval

__all__ = ('check_transformation', 'NoTransformation', 'Transformation')

class NoTransformation:
    '''An invariant transformation.  Model evaluations are done in the
    ext(ernal) coordinate, but within the optimization algorthm calcuations
    are done in the int(ernal) coordinate system. So final output from the
    optimization must re-transform from internal to external coordinate.
    However, since this transformation is invariant the internal and external
    coordinates are identical.
    '''

    def __init__(self):
        self.lo = None
        self.hi = None
        return

    def calc_covar(self, m, n, xint, info, fjac):
        covar = None
        if info > 0:
            fjac = numpy.reshape(numpy.ravel(fjac, order='F'),
                                 (m, n), order='F')
            if m != n:
                covar = fjac[:n, :n]
            else:
                covar = fjac
        return covar

    def chain_rule(self, x):
        return numpy.ones_like(x)

    def ext2int(self, x):
        return x

    def int2ext(self, x):
        return x

    def set_bounds(self, lo, hi):
        self.lo = lo
        self.hi = hi


class Transformation(NoTransformation):
    ''' See the doc of the NoTransformation base class.
    This class enables an implementation of an unbounded optimization algorithm
    into a bounded one. It is capable of transforming each free parameter using
    the appropiate transformation depending on its boundary conditions of each
    free parameter.
    It is recommended that the class should only be used if a free parameter,
    after a fit, is sufficiently close to upper and/or lower limits. The reason
    for this is some Sherpa model parameters have large, upper and/or lower
    limits and one could lose resolution when transforming. A parameter
    with a low, equal to -hugeval, or high bound,  equal to hugeval,
    is considered to have only a low or high bound only.
    Moreover, there is a performance penalty for transforming
    between the external and internal parameter coordinates since a
    transformation has to be done for each parameter for each iteration
    of the fit.
    '''
    def __init__(self):
        NoTransformation.__init__(self)
        self.bounds = None
        self.int2ext_fcns = None
        self.ext2int_fcns = None
        self.chain_fcns = None
        return

    def init(self):
        '''For each parameter, decide which transformation to use.
        A model has multiple parameters since parameters might be
        bounded, while others are unbounded. There are four different
        scenarios:
            0) A parameter is unbounded
            1) A parameter only has a high bound
            2) A parameter only has a low bound
            3) A parameter has both low and high bound
        '''
        bounds = []
        int2ext = []
        ext2int = []
        chain = []
        # cause it's faster and more accurate than cmp with max hugeval
        myhugeval = 1.0e-1 * hugeval
        for a, b in zip(self.lo, self.hi):
            if a <= - myhugeval and b >= myhugeval:
                # 0) A parameter is unbounded
                bounds.append((None, None))
                int2ext.append(self.return_arg)
                ext2int.append(self.return_arg)
                chain.append(self.chain_rule_1)
            elif a <= - myhugeval:
                # 1) A parameter only has a high bound
                bounds.append((None, b))
                int2ext.append(self.int2ext_hi)
                ext2int.append(self.ext2int_hi)
                chain.append(self.chain_rule_hi)
            elif b >= myhugeval:
                # 2) A parameter only has a low bound
                bounds.append((a, None))
                int2ext.append(self.int2ext_lo)
                ext2int.append(self.ext2int_lo)
                chain.append(self.chain_rule_lo)
            else:
                # 3) A parameter has both low and high  bound
                bounds.append((a, b))
                int2ext.append(self.int2ext_both)
                ext2int.append(self.ext2int_both)
                chain.append(self.chain_rule_both)
        return bounds, int2ext, ext2int, chain

    def calc_covar(self, m, n, xint, info, fjac):
        '''apply the chain rule to the internal covariance matrix

        From https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm:
            ...The Jacobian matrix as defined above is not (in general)
            a square matrix, but a rectangular matrix of size m x n,
            where n is the number of parameter (size of the vector B).
            The matrix multiplication (J.T J) yields the required
            n x n square matarix...
        '''
        covar = super().calc_covar(m, n, xint, info, fjac)
        chain = self.chain_rule(xint)
        mychain = numpy.atleast_2d(chain)
        grad = numpy.dot(mychain.T, mychain)
        covar *= grad
        return covar

    def call_fcns(self, x, fcns):
        '''The common function (used by: chain_rule, ext2int and int2ext)
        to iterate through the appropiate transformation or chain_rule
        for each parameter
        '''
        result = numpy.empty_like(x)
        for ii, (fcn, arg) in enumerate(zip(fcns, x)):
            result[ii] = fcn(ii, arg)
        return result

    def chain_rule(self, x):
        return self.call_fcns(x, self.chain_fcns)

    def chain_rule_both(self, ii, x):
        '''The chain rule for a parameter that has both limits'''
        return 0.5 * numpy.cos(x) * (self.hi[ii] - self.lo[ii])

    def chain_rule_hi(self, ii, x):
        '''The chain rule for a parameter that only has an upper bound'''
        return 2.0 * (self.hi[ii] - x)

    def chain_rule_lo(self, ii, x):
        '''The chain rule for a parameter that only has an lower bound'''
        return 2.0 * (x - self.lo[ii])

    def chain_rule_1(self, ii, x):
        '''The chain rule for an unbounded parameter'''
        return 1.0

    def ext2int(self, x):
        return self.call_fcns(x, self.ext2int_fcns)

    def ext2int_both(self, ii, x):
        '''The inverse transformation for the internal to external parameter
        with upper and lower limits
        '''
        x_a = x - self.lo[ii]
        b_a = self.hi[ii] - self.lo[ii]
        arg = 2.0 * x_a / b_a  - 1.0
        return numpy.arcsin(arg)

    def ext2int_lo(self, ii, x):
        '''The inverse transformation for the internal to external parameter
        with lower limit
        '''
        return numpy.sqrt(x -self.lo[ii]) + self.lo[ii]

    def ext2int_hi(self, ii, x):
        '''The inverse transformation for the internal to external parameter
        with upper limit
        '''
        return self.hi[ii] - numpy.sqrt(self.hi[ii] - x)

    def int2ext_both(self, ii, x):
        '''The transformation from external to internal coordinate for a
        parameter with upper and lower limits via the following transformation:

        P    = lo + (hi - lo) * (sin(P   ) + 1) / 2.0
         ext                          int
        '''
        b_a = self.hi[ii] - self.lo[ii]
        sin_p1 = numpy.sin(x) + 1.0
        return self.lo[ii] + 0.5 *  b_a * sin_p1

    def int2ext_hi(self, ii, x):
        '''The transformation from external to internal coordinate for a
        parameter with only an upper limit via the following transformation:

        P    = hi - (hi - P   )^2
         ext               int
        '''
        return self.hi[ii] - (self.hi[ii] - x)**2

    def int2ext_lo(self, ii, x):
        '''The transformation from external to internal coordinate for a
        parameter with only a lower limit via the following transformation:

        P    = lo - (P    - lo)^2
         ext          int
        '''
        return self.lo[ii] + (x - self.lo[ii])**2

    def int2ext(self, x):
        return self.call_fcns(x, self.int2ext_fcns)

    def return_arg(self, ii, x):
        '''An unbounded parameter and therefore a no-op'''
        return x

    def set_bounds(self, lo, hi):
        super().set_bounds(lo, hi)
        self.bounds, self.int2ext_fcns, self.ext2int_fcns, \
            self.chain_fcns = self.init()
        return

def check_transformation(arg, xmin, xmax):
    tfmt = None
    if arg is False:
        tfmt = NoTransformation()
    elif arg is True or isinstance(arg, Transformation):
        tfmt = Transformation()
    elif inspect.isclass(arg) and issubclass(arg, Transformation):
        tfmt = arg
    else:
        raise TypeError("arg must inherits from class Transformation")
    tfmt.set_bounds(xmin, xmax)
    return tfmt
