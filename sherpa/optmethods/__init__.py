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

import logging
import numpy
from sherpa.utils import NoNewAttributesAfterInit, \
     get_keyword_names, get_keyword_defaults, print_fields
from sherpa.utils.err import FitErr
from sherpa.optmethods.optfcts import *

warning = logging.getLogger(__name__).warning


__all__ = ('GridSearch', 'OptMethod', 'LevMar', 'MonCar', 'NelderMead')


class OptMethod(NoNewAttributesAfterInit):

    def __init__(self, name, optfunc):        
        self.name = name
	self._optfunc = optfunc
        self.config = self.default_config
        NoNewAttributesAfterInit.__init__(self)

    def __getattr__(self, name):
        if name in self.__dict__.get('config', ()):
            return self.config[name]
        raise AttributeError("'%s' object has no attribute '%s'" %
                             (type(self).__name__, name))

    def __setattr__(self, name, val):
        if name in self.__dict__.get('config', ()):
            self.config[name] = val
        else:
            NoNewAttributesAfterInit.__setattr__(self, name, val)

    def __repr__(self):
        return ("<%s optimization method instance '%s'>" %
                (type(self).__name__, self.name))

    # Need to support users who have pickled sessions < CIAO 4.2
    def __setstate__(self, state):
        new_config = get_keyword_defaults(state.get('_optfunc'))
        old_config = state.get('config', {})

        # remove old kw args from opt method dict
        for key in old_config.keys():
            if key not in new_config:
                old_config.pop(key)

        # add new kw args with defaults
        for key, val in new_config.items():
            if key not in old_config:
                old_config[key] = val

        self.__dict__.update(state)

    def __str__(self):
        names = ['name']
        names.extend(get_keyword_names(self._optfunc))
        #names.remove('full_output')
        # Add the method's name to printed output
        # Don't add to self.config b/c name isn't a
        # fit function config setting
        add_name_config = {}
        add_name_config['name'] = self.name
        add_name_config.update(self.config)
        return print_fields(names, add_name_config)

    def _get_default_config(self):
        args = get_keyword_defaults(self._optfunc)
	return args
    default_config = property(_get_default_config)

    def fit(self, statfunc, pars, parmins, parmaxes, statargs=(),
            statkwargs={}):

        def cb(pars):
            return statfunc(pars, *statargs, **statkwargs)

	output = self._optfunc(cb, pars, parmins, parmaxes, **self.config)

        success = output[0]
        msg = output[3]
        if not success:
            warning('fit failed: %s' % msg)

        # Ensure that the best-fit parameters are in an array.  (If there's
        # only one, it might be returned as a bare float.)
        output = list(output)
        output[1] = numpy.asarray(output[1]).ravel()
        output = tuple(output)

	return output

class GridSearch(OptMethod):
    """A simple iterative method to support the template model interface,
    the method can be used for non-template model but it is very ineffecient
    for this purpose."""
    
    def __init__(self, name='gridsearch'):
	OptMethod.__init__(self, name, grid_search)

class LevMar(OptMethod):
    """
  LMDIF.

                                                                 Page 1

               Documentation for MINPACK subroutine LMDIF

                        Double precision version

                      Argonne National Laboratory

         Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

                               March 1980


 1. Purpose.

       The purpose of LMDIF is to minimize the sum of the squares of M
       nonlinear functions in N variables by a modification of the
       Levenberg-Marquardt algorithm.  The user must provide a subrou-
       tine which calculates the functions.  The Jacobian is then cal-
       culated by a forward-difference approximation.


 2. Subroutine and type statements.

       SUBROUTINE LMDIF(FCN,M,N,X,FVEC,FTOL,XTOL,GTOL,MAXFEV,EPSFCN,
      *                 DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC,
      *                 IPVT,QTF,WA1,WA2,WA3,WA4)
       INTEGER M,N,MAXFEV,MODE,NPRINT,INFO,NFEV,LDFJAC
       INTEGER IPVT(N)
       DOUBLE PRECISION FTOL,XTOL,GTOL,EPSFCN,FACTOR
       DOUBLE PRECISION X(N),FVEC(M),DIAG(N),FJAC(LDFJAC,N),QTF(N),
      *                 WA1(N),WA2(N),WA3(N),WA4(M)
       EXTERNAL FCN


 3. Parameters.

       Parameters designated as input parameters must be specified on
       entry to LMDIF and are not changed on exit, while parameters
       designated as output parameters need not be specified on entry
       and are set to appropriate values on exit from LMDIF.

       FCN is the name of the user-supplied subroutine which calculates
         the functions.  FCN must be declared in an EXTERNAL statement
         in the user calling program, and should be written as follows.

         SUBROUTINE FCN(M,N,X,FVEC,IFLAG)
         INTEGER M,N,IFLAG
         DOUBLE PRECISION X(N),FVEC(M)
         ----------
         CALCULATE THE FUNCTIONS AT X AND
         RETURN THIS VECTOR IN FVEC.
         ----------
         RETURN
         END


                                                                 Page 2

         The value of IFLAG should not be changed by FCN unless the
         user wants to terminate execution of LMDIF.  In this case set
         IFLAG to a negative integer.

       M is a positive integer input variable set to the number of
         functions.

       N is a positive integer input variable set to the number of
         variables.  N must not exceed M.

       X is an array of length N.  On input X must contain an initial
         estimate of the solution vector.  On output X contains the
         final estimate of the solution vector.

       FVEC is an output array of length M which contains the functions
         evaluated at the output X.

       FTOL is a nonnegative input variable.  Termination occurs when
         both the actual and predicted relative reductions in the sum
         of squares are at most FTOL.  Therefore, FTOL measures the
         relative error desired in the sum of squares.  Section 4 con-
         tains more details about FTOL.

       XTOL is a nonnegative input variable.  Termination occurs when
         the relative error between two consecutive iterates is at most
         XTOL.  Therefore, XTOL measures the relative error desired in
         the approximate solution.  Section 4 contains more details
         about XTOL.

       GTOL is a nonnegative input variable.  Termination occurs when
         the cosine of the angle between FVEC and any column of the
         Jacobian is at most GTOL in absolute value.  Therefore, GTOL
         measures the orthogonality desired between the function vector
         and the columns of the Jacobian.  Section 4 contains more
         details about GTOL.

       MAXFEV is a positive integer input variable.  Termination occurs
         when the number of calls to FCN is at least MAXFEV by the end
         of an iteration.

       EPSFCN is an input variable used in determining a suitable step
         for the forward-difference approximation.  This approximation
         assumes that the relative errors in the functions are of the
         order of EPSFCN.  If EPSFCN is less than the machine preci-
         sion, it is assumed that the relative errors in the functions
         are of the order of the machine precision.

       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
         internally set.  If MODE = 2, DIAG must contain positive
         entries that serve as multiplicative scale factors for the
         variables.

       MODE is an integer input variable.  If MODE = 1, the variables
         will be scaled internally.  If MODE = 2, the scaling is


                                                                 Page 3

         specified by the input DIAG.  Other values of MODE are equiva-
         lent to MODE = 1.

       FACTOR is a positive input variable used in determining the ini-
         tial step bound.  This bound is set to the product of FACTOR
         and the Euclidean norm of DIAG*X if nonzero, or else to FACTOR
         itself.  In most cases FACTOR should lie in the interval
         (.1,100.).  100. is a generally recommended value.

       NPRINT is an integer input variable that enables controlled
         printing of iterates if it is positive.  In this case, FCN is
         called with IFLAG = 0 at the beginning of the first iteration
         and every NPRINT iterations thereafter and immediately prior
         to return, with X and FVEC available for printing.  If NPRINT
         is not positive, no special calls of FCN with IFLAG = 0 are
         made.

       INFO is an integer output variable.  If the user has terminated
         execution, INFO is set to the (negative) value of IFLAG.  See
         description of FCN.  Otherwise, INFO is set as follows.

         INFO = 0  Improper input parameters.

         INFO = 1  Both actual and predicted relative reductions in the
                   sum of squares are at most FTOL.

         INFO = 2  Relative error between two consecutive iterates is
                   at most XTOL.

         INFO = 3  Conditions for INFO = 1 and INFO = 2 both hold.

         INFO = 4  The cosine of the angle between FVEC and any column
                   of the Jacobian is at most GTOL in absolute value.

         INFO = 5  Number of calls to FCN has reached or exceeded
                   MAXFEV.

         INFO = 6  FTOL is too small.  No further reduction in the sum
                   of squares is possible.

         INFO = 7  XTOL is too small.  No further improvement in the
                   approximate solution X is possible.

         INFO = 8  GTOL is too small.  FVEC is orthogonal to the
                   columns of the Jacobian to machine precision.

         Sections 4 and 5 contain more details about INFO.

       NFEV is an integer output variable set to the number of calls to
         FCN.

       FJAC is an output M by N array.  The upper N by N submatrix of
         FJAC contains an upper triangular matrix R with diagonal ele-
         ments of nonincreasing magnitude such that


                                                                 Page 4

                T     T           T
               P *(JAC *JAC)*P = R *R,

         where P is a permutation matrix and JAC is the final calcu-
         lated J"""


    def __init__(self, name='levmar'):
	OptMethod.__init__(self, name, lmdif)


class MonCar(OptMethod):

    def __init__(self, name='moncar'):
	OptMethod.__init__(self, name, montecarlo)


# Call Sherpa's Nelder-Mead implementation
class NelderMead(OptMethod):

    def __init__(self, name='simplex'):
	OptMethod.__init__(self, name, neldermead)


###############################################################################

## import sherpa.optmethods.myoptfcts
## from sherpa.optmethods.fminpowell import *
## from sherpa.optmethods.nmpfit import *

## from sherpa.optmethods.odrpack import odrpack
## from sherpa.optmethods.stogo import stogo
## from sherpa.optmethods.chokkan import chokkanlbfgs
## from sherpa.optmethods.odr import odrf77

## def myall( targ, arg ):
##     fubar = list( targ )
##     fubar.append( arg )
##     return tuple( fubar )

## __all__ = myall( __all__, 'Bobyqa' )
## __all__ = myall( __all__, 'Chokkan' )
## __all__ = myall( __all__, 'cppLevMar' )
## __all__ = myall( __all__, 'Dif_Evo' )
## __all__ = myall( __all__, 'MarLev' )
## __all__ = myall( __all__, 'MyMinim' )
## __all__ = myall( __all__, 'Nelder_Mead' )
## __all__ = myall( __all__, 'NMPFIT' )        
## __all__ = myall( __all__, 'Newuoa' )
## __all__ = myall( __all__, 'Odr' )
## __all__ = myall( __all__, 'OdrPack' )
## __all__ = myall( __all__, 'PortChi' )
## __all__ = myall( __all__, 'PortFct' )    
## __all__ = myall( __all__, 'ScipyPowell' )
## __all__ = myall( __all__, 'StoGo' )

## class Bobyqa(OptMethod):
##     def __init__(self, name='bobyqa'):
##         OptMethod.__init__(self, name, myoptfcts.bobyqa)

## class Chokkan(OptMethod):
##     def __init__(self, name='chokkan'):
##         OptMethod.__init__(self, name, chokkanlbfgs)

## class cppLevMar(OptMethod):

##    def __init__(self, name='clevmar'):
## 	OptMethod.__init__(self, name, optfcts.lmdif_cpp)
    
## class Dif_Evo(OptMethod):
##     def __init__(self, name='dif_evo'):
##         OptMethod.__init__(self, name, myoptfcts.dif_evo)
                 
## class MarLev(OptMethod):
##     def __init__(self, name='marlev'):
##         OptMethod.__init__(self, name, myoptfcts.marquadt_levenberg)

## class MyMinim(OptMethod):

##     def __init__(self, name='simplex'):
## 	OptMethod.__init__(self, name, minim)

## class Nelder_Mead(OptMethod):
##     def __init__(self, name='nelder_mead'):
##         OptMethod.__init__(self, name, myoptfcts.nelder_mead)

## class Newuoa(OptMethod):
##     def __init__(self, name='newuoa'):
##         OptMethod.__init__(self, name, myoptfcts.newuoa)

## class NMPFIT(OptMethod):
##     def __init__(self, name='pytools_nmpfit'):
##         OptMethod.__init__(self, name, nmpfit.pytools_nmpfit)

## class OdrPack(OptMethod):
##     def __init__(self, name='odrpack'):
##         OptMethod.__init__(self, name, odrpack)

## class Odr(OptMethod):
##     def __init__(self, name='odr'):
##         OptMethod.__init__(self, name, odrf77)

## class PortChi(OptMethod):
##     def __init__(self, name='dn2fb'):
##         OptMethod.__init__(self, name, myoptfcts.dn2fb)

## class PortFct(OptMethod):
##     def __init__(self, name='dmnfb'):
##         OptMethod.__init__(self, name, myoptfcts.dmnfb)
    
## class ScipyPowell(OptMethod):
##     def __init__(self, name='scipypowell'):
##         OptMethod.__init__(self, name, my_fmin_powell)

## class StoGo(OptMethod):
##     def __init__(self, name='stogo'):
## 	OptMethod.__init__(self, name, stogo)
        
###############################################################################
