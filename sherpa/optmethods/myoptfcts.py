import numpy
import random

#from sherpa.optmethods import _port
#from sherpa.optmethods import _powell
import _odrpack

from sherpa.optmethods import OptMethod
#from sherpa.astro.ui import ui

from sherpa.utils import func_counter, Knuth_close, parallel_map
from sherpa.optmethods.optfcts import _check_args, _get_saofit_msg, \
     _outside_limits, _my_is_nan, EPSILON, \
     neldermead, difevo_nm, FUNC_MAX, sao_fcmp, _narrow_limits

def select_result( results, npar, debug ):
    mynfev = 0
    myfval = FUNC_MAX
    myx = numpy.zeros( npar )
    for result in results:
        mynfev += int( result[ -1 ] )
        if result[ -2 ] < myfval:
            myx = result[ 0:npar ]
            myfval = result[ -2 ]

    if False != debug:
        print 'select_result: f%s=%e in %d nfevs' % \
              ( myx, myfval, mynfev )
    return myx, myfval, mynfev
        
def func_bounds( func, low, high, myinf=numpy.inf ):
    """In order to keep the corrent number of function evaluations:
    func_counter should be called before func_bounds. For example,
    the following code
      x0 = [ -1.2, 1.0 ]
      xmin = [ -10.0, -10.0 ]
      xmax = [ 10.0, 10.0 ]
      nfev, rosen = func_counter( Rosenbrock )
      rosenbrock = func_bounds( rosen, xmin, xmax )
      print rosenbrock( [-15.0, 1.0] ), nfev[0]
    should output:
    inf 0
    """
    def func_bounds_wrapper( x, *args ):
        #print 'low = ', low, '\thigh = ', high
        if _outside_limits( x, low, high ):
            return myinf
        return func( x, *args )
    return func_bounds_wrapper

#
################################### ODRPACK ###################################
#
#
# Odrpack
#
def odrpack( fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0 ):

    def stat_cb1( pars ):
        return fcn( pars )[ 1 ]

    x, xmin, xmax = _check_args(x0, xmin, xmax)
    
    n = numpy.asanyarray(stat_cb1(x)).size

    if maxfev is None:
        maxfev = 512 * len( x )

    orig_fcn = stat_cb1
    def stat_cb1( x_new ):
        fvec = orig_fcn( x_new )
        return fvec
    
    info, nfev, fval, covarerr = _odrpack.odrpack( stat_cb1, n, x, xmin, xmax,
                                                   verbose, ftol )

    ierr = 4
    if 1 <= info and info <= 3:
        ierr = 0

    status, msg = _get_saofit_msg( maxfev, ierr )
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev, 'covarerr': covarerr})

    return rv

#
################################### ODRPACK ###################################
#

#
#################################### PORT ####################################
#
def get_port_err_msg( ierr ):
    key = {
        3: (True, '***** PARAMETER-CONVERGENCE *****'),
        4: (True, '***** RELATIVE FUNCTION CONVERGENCE *****'),
        5: (True, '***** PARAMETER- AND RELATIVE FUNCTION CONVERGENCE *****'),
        6: (True, '***** ABSOLUTE FUNCTION CONVERGENCE *****'),
        7: (True, '***** SINGULAR CONVERGENCE *****'),
        8: (False, '***** FALSE CONVERGENCE *****'),
        9: (False, '***** FUNCTION EVALUATION LIMIT *****'),
        10: (False, '***** ITERATION LIMIT *****') }

    status, msg = key.get( ierr, (False, 'unknown status flag (%d)' % ierr))
    return status, msg

def dmnfb( fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0 ):
    
    x, xmin, xmax = _check_args( x0, xmin, xmax )

    if maxfev is None:
        maxfev = 512 * len( x0 )

    bounds = numpy.ndarray( (2,len(x0)), order='Fortran' )
    bounds[ 0, : ] = xmin.copy( )
    bounds[ 1, : ] = xmax.copy( )
    nfev, fval, ierr, dpptri, covarerr = _port.mydmnfb( x, bounds, ftol, ftol,
                                                        ftol, maxfev,
                                                        verbose, fcn )
    status, msg = get_port_err_msg( ierr )
    rv = ( status, x, fval )
    rv += (msg, {'info':status, 'nfev': nfev, 'covarerr': covarerr})
    
    return rv


def dn2fb( fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0 ):
    
    def stat_cb1( pars ):
        return fcn( pars )[ 1 ]

    x, xmin, xmax = _check_args( x0, xmin, xmax )

    m = numpy.asanyarray(stat_cb1(x)).size
    npar = len( x0 )
    if maxfev is None:
        maxfev = 512 * npar

    bounds = numpy.ndarray( ( 2, npar ), order='Fortran' )
    bounds[ 0, : ] = xmin.copy( )
    bounds[ 1, : ] = xmax.copy( )
    nfev, fval, ierr = _port.mydn2fb( m, x, bounds, ftol, ftol, ftol, maxfev,
                                      verbose, stat_cb1 )
    status, msg = get_port_err_msg( ierr )    
    rv = ( status, x, fval )
    rv += (msg, {'info':status, 'nfev': nfev,})
    
    return rv

#
#################################### PORT ####################################
#

#
#################################### MarLev ###################################
#

def marquadt_levenberg( fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None,
                         verbose=0 ):

    debug = False
    
    x, xmin, xmax = _check_args( x0, xmin, xmax )

    npar = len( x0 )
    if maxfev is None:
        maxfev = 512 * npar

    odrpack_result = odrpack( fcn, x, xmin, xmax, ftol=ftol, maxfev=maxfev,
                              verbose=verbose )
    x = numpy.asarray( odrpack_result[ 1 ], numpy.float_ )
    fval = odrpack_result[2]
    nfev = odrpack_result[4].get( 'nfev' )
    if debug:
        print 'odrpack: f%s = %e in %d nfevs' % ( x, fval, nfev )

    dn2fb_result = dn2fb( fcn, x, xmin, xmax, ftol=ftol, maxfev=maxfev-nfev, 
                          verbose=verbose )
    x = numpy.asarray( dn2fb_result[ 1 ], numpy.float_ )
    fval = dn2fb_result[2]
    nfev += dn2fb_result[4].get( 'nfev' )
    if debug:
        print 'dn2fb: f%s = %e in %d nfevs' % ( x, fval, dn2fb_result[4].get( 'nfev' ) )

    odrpack_result = odrpack( fcn, x, xmin, xmax, ftol=ftol,
                              maxfev=maxfev-nfev, verbose=verbose )
    x = numpy.asarray( odrpack_result[ 1 ], numpy.float_ )
    info = odrpack_result[4].get( 'info' )
    fval = odrpack_result[2]
    covarerr = odrpack_result[4].get( 'covarerr' )
    nfev += odrpack_result[4].get( 'nfev' )
    if debug:
        print 'odrpack: f%s = %e in %d nfevs' % ( x, fval, odrpack_result[4].get( 'nfev' ) )

    ierr = 4
    if 1 <= info and info <= 3:
        ierr = 0

    status, msg = _get_saofit_msg( maxfev, ierr )
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev, 'covarerr': covarerr})

    return rv

    
#
#################################### MarLev ###################################
#

#
################################### POWELL ###################################
#
## http://www.mat.univie.ac.at/~neum/glopt/contrib/bobyqa.txt
## From: "M.J.D. Powell" <mjdp@cam.ac.uk>
## Date: Wed, 7 Jan 2009 07:33:47 -0500
## Subject: Announcement of release of BOBYQA

## This Fortran software seeks the least value of a function of several variables
## without requiring any derivatives of the objective function. It was developed
## from my NEWUOA package for this calculation in the unconstrained case. The
## main new feature of BOBYQA, however, is that it allows lower and upper bounds
## on each variable. The name BOBYQA denotes Bound Optimization BY Quadratic
## Approximation.

## Please send an e-mail to me at mjdp@cam.ac.uk if you would like to receive a
## free copy of the Fortran software. As far as I know BOBYQA is the most
## powerful package available at present for minimizing functions of hundreds of
## variables without derivatives subject to simple bound constraints. There are
## no restrictions on its use. I would be delighted if it becomes valuable to
## much research and many applications.

## January 7th, 2009                          Michael J.D. Powell.
def bobyqa( fcn, x0, xmin, xmax, npt=None, rhobeg=None, rhoend=None,
            ftol=EPSILON, maxfev=None, verbose=-1 ):
    """    
C
C     This subroutine seeks the least value of a function of many variables,
C     by applying a trust region method that forms quadratic models by
C     interpolation. There is usually some freedom in the interpolation
C     conditions, which is taken up by minimizing the Frobenius norm of
C     the change to the second derivative of the model, beginning with the
C     zero matrix. The values of the variables are constrained by upper and
C     lower bounds. The arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in
C       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
C       recommended.
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
C       bounds, respectively, on X(I). The construction of quadratic models
C       requires XL(I) to be strictly less than XU(I) for each I. Further,
C       the contribution to a model from changes to the I-th variable is
C       damaged severely by rounding errors if XU(I)-XL(I) is too small.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND no greater than
C       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
C       expected change to a variable, while RHOEND should indicate the
C       accuracy that is required in the final values of the variables. An
C       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
C       is less than 2*RHOBEG.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
C
C     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
C     F to the value of the objective function for the current values of the
C     variables X(1),X(2),...,X(N), which are generated automatically in a
C     way that satisfies the bounds given in XL and XU.
C
C     Return if the value of NPT is unacceptable.
C
"""
    x, xmin, xmax = _check_args( x0, xmin, xmax )

    if rhobeg is None:
        rhobeg = 0.995 * min( min( xmax - xmin ) / 2.0, 10.0 )
    if rhoend is None:
        rhoend = ftol
    if maxfev is None:
        maxfev = 1024 * len(x)
    n = len( x )
    if npt is None:
        npt = 2 * n + 1
    if npt < n + 2 or npt > n * ( n + 3 ) / 2 + 1:
        npt = n + 2

##     orig_fcn = fcn
##     def fcn(x_new):
##         if _outside_limits(x_new, xmin, xmax):
##             return 1.0e128
##         return orig_fcn(x_new)

    mynfev, fcn_counter = func_counter( fcn )
    myfcn = func_bounds( fcn_counter, xmin, xmax, myinf=1.0e128 )

    ierr = 0
    if len(x) >= 2:
        nfev, fval = _powell.bobyqa( npt, x, xmin, xmax, rhobeg, rhoend,
                                     verbose, maxfev, myfcn )
    else:
        fval = fcn(x)
        nfev = 1
        ierr = 1

    if nfev >= maxfev:
        ierr = 3

    status, msg = _get_saofit_msg( maxfev, ierr )
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev})
    return rv

def newuoa( fcn, x0, xmin, xmax, npt=None, rhobeg=None, rhoend=None,
            ftol=EPSILON, maxfev=None, verbose=-1 ):
    """    
C
C     This subroutine seeks the least value of a function of many variables,
C     by a trust region method that forms quadratic models by interpolation.
C     There can be some freedom in the interpolation conditions, which is
C     taken up by minimizing the Frobenius norm of the change to the second
C     derivative of the quadratic model, beginning with a zero matrix. The
C     arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in the
C       interval [N+2,(N+1)(N+2)/2].
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C       RHOBEG should be about one tenth of the greatest expected change to a
C       variable, and RHOEND should indicate the accuracy that is required in
C       the final values of the variables.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
C
C     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
C     the value of the objective function for the variables X(1),X(2),...,X(N).
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
"""
    x, xmin, xmax = _check_args( x0, xmin, xmax )

    if rhobeg is None:
        rhobeg = 10.0
    if rhoend is None:
        rhoend = ftol
    if maxfev is None:
        maxfev = 1024 * len(x)
    n = len( x )
    if npt is None:
        npt = ( n * ( n + 5 ) + 7 ) / 4

    if npt < n + 2 or npt > n * ( n + 3 ) / 2 + 1:
        npt = 2*n + 1
        
##     orig_fcn = fcn
##     def fcn(x_new):
##         if _outside_limits(x_new, xmin, xmax):
##             return 1.0e128
##         return orig_fcn(x_new)

    mynfev, fcn_counter = func_counter( fcn )
    myfcn = func_bounds( fcn_counter, xmin, xmax, myinf=1.0e128 )
    
    ierr = 0
    if len(x) >= 2:
        nfev, fval = _powell.newuoa( npt, x, rhobeg, rhoend, verbose, maxfev,
                                     myfcn )
    else:
        fval = fcn(x)
        nfev = 1
        ierr = 1

    if nfev >= maxfev:
        ierr = 3

    status, msg = _get_saofit_msg( maxfev, ierr )
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev})
    return rv

def powell( fcn, x0, xmin, xmax, npt=None, rhobeg=None, rhoend=None,
            algorithm='bobyqa', ftol=EPSILON, maxfev=None, verbose=-1,):
    
    if 'newuoa' == algorithm:
        return newuoa( fcn, x0, xmax, ftol=ftol, maxfev=maxfev, npt=None,
                       rhobeg=rhobeg, rhoend=rhoend, verbose=verbose )
    else:
        return bobyqa( fcn, x0, xmax, ftol=ftol, maxfev=maxfev, npt=None,
                       rhobeg=rhobeg, rhoend=rhoend, verbose=verbose )        
#
################################### POWELL ###################################
#

#
################################## PMonCar ##################################
#
def mc_montecarlo( fcn, x0, xmin, xmax, seed=74815, population_size=None,
                   xprob=0.9, weighting_factor=0.8, ftol=EPSILON,
                   maxfev=None, verbose=0 ):

    def myopt( fcn, xxx, ftol, maxfev, seed, pop, xprob, strategy,
               weight, limit_factor=4.0, debug=False ):

        def random_start( xmin, xmax ):
            xx = []
            for ii in range( len(xmin ) ):
                xx.append( random.uniform( xmin[ ii ], xmax[ ii ] ) )
            return numpy.asarray( xx )

        def multicore_difevo_nm( myx, myxmin, myxmax, mymaxfev, mystrategy,
                                 debug ):

            def mc_difevo_nm( args ):

                npar = int( args[ 0 ] )
                mystrategy = int( args[ 1 ] )
                mymaxfev = int( args[ 2 ] )

                myx = args[ 3: 3 + npar ]
                myxmin = args[ 3 + npar: 3 + 2 * npar ]
                myxmax = args[ -npar : ]

##                 print 'mc_difevo_nm: args = ', args
##                 print 'mc_difevo_nm: npar = ', npar
##                 print 'mc_difevo_nm: mystrategy = ', mystrategy
##                 print 'mc_difevo_nm: mymaxfev = ', mymaxfev
##                 print 'mc_difevo_nm: x = ', myx
##                 print 'mc_difevo_nm: xmin = ', myxmin
##                 print 'mc_difevo_nm: xmax = ', myxmax
                 
                verbose = 0
                seed = 2005
                xprob=0.9
                weighting_factor=0.8
                result = difevo_nm( fcn, myx, myxmin, myxmax, ftol=ftol,
                                    maxfev=mymaxfev, verbose=verbose,
                                    seed=seed, population_size=None,
                                    xprob=xprob, strategy=mystrategy,
                                    weighting_factor=weighting_factor )
                mypar = numpy.asarray( result[ 1 ], numpy.float_ )
                myfval = result[2]
                mynfev = result[4].get( 'nfev' )
                return numpy.concatenate( (mypar, [myfval,mynfev]) )


            npar = len( myx )

            arg1 = numpy.concatenate( ([npar,mystrategy,mymaxfev],
                                       myx, myxmin, myxmax) )
            tmpx = random_start( myxmin, myxmax )
            arg2 = numpy.concatenate( ([npar,(mystrategy+5)%10,mymaxfev],
                                       tmpx, myxmin, myxmax) )
            args = [ arg1, arg2 ]
            results = parallel_map( mc_difevo_nm, args )

            return select_result( results, npar, debug )


        x = xxx[ 0 ]
        xmin = xxx[ 1 ]
        xmax = xxx[ 2 ]
        maxfev_per_iter = 512 * x.size


        ############################# NelderMead #############################
        mymaxfev = min( maxfev_per_iter, maxfev )
        result = neldermead( fcn, x, xmin, xmax, maxfev=mymaxfev, ftol=ftol,
                             finalsimplex=9 )
        x = numpy.asarray( result[ 1 ], numpy.float_ )
        nfval = result[2]
        nfev = result[4].get( 'nfev' )
        if verbose or False != debug:
            if debug:
                print 'neldermead: ',
            print 'f%s=%.14e in %d nfevs\n' % ( x, nfval, nfev )
        ############################# NelderMead #############################

        ############################## nmDifEvo ##############################
        myxmin, myxmax = _narrow_limits( limit_factor, [x,xmin,xmax],
                                        debug=debug )
        mymaxfev = min( maxfev_per_iter, maxfev - nfev )
        result = multicore_difevo_nm( x, myxmin, myxmax, mymaxfev, strategy,
                                      debug )
        x = result[ 0 ]
        ofval = result[ 1 ]
        nfev += result[ 2 ]
        if verbose or False != debug:
            if debug:
                print 'neldermead: ',
            print 'f%s=%.14e in %d nfevs\n' % ( x, ofval, nfev )
        ############################## nmDifEvo ##############################
            
        while nfev < maxfev:

            strategy = ( strategy + 1 ) % 10
            limit_factor *= 2
            myxmin, myxmax = _narrow_limits( limit_factor, [x,xmin,xmax],
                                             debug=debug )

            ############################ nmDifEvo #############################
            result = multicore_difevo_nm( x, myxmin, myxmax, mymaxfev,
                                          strategy, debug )
            x = result[ 0 ]
            nfval = result[ 1 ]
            nfev += result[ 2 ]
            if verbose or False != debug:
                print 'difevo_nm: f%s=%.14e in %d nfevs\n' % ( x, ofval,nfev)
            ############################ nmDifEvo #############################


            if verbose or False != debug:
                print 'ofval=%.14e\tnfval=%.14e\n' % (ofval, nfval)

            if sao_fcmp( ofval, nfval, ftol ) <= 0:
                return x, nfval, nfev
            ofval = nfval
            
        return x, nfval, nfev

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = max( 0.1, xprob )
    xprob = min( xprob, 1.0 )

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = max( 0.1, weighting_factor )
    weighting_factor = min( weighting_factor, 1.0 )

    strategy=0
    # make sure strategy is within [0,9]
    strategy = max( 0, strategy )
    strategy = min( strategy, 9 )

    if seed is None:
        seed = randint(0, 2147483648L) # pow(2,31) == 2147483648L
    if population_size is None:
        population_size = max( population_size, 16 * x.size )

    if maxfev is None:
        maxfev = 8192 * population_size


    x, fval, nfev = myopt( fcn, [x, xmin, xmax], numpy.sqrt(ftol), maxfev,
                           seed, population_size, xprob, strategy,
                           weighting_factor, limit_factor=2.0, debug=False )

    if nfev < maxfev:
        result = neldermead( fcn, x, xmin, xmax,
                             maxfev=min( 512*len(x), maxfev - nfev ),
                             ftol=ftol, finalsimplex=9 )

        x = numpy.asarray( result[ 1 ], numpy.float_ )
        fval = result[2]
        nfev += result[4].get( 'nfev' )

    ierr = 0
    if nfev >= maxfev:
        ierr = 3
    status, msg = _get_saofit_msg( maxfev, ierr )

    rv = (status, x, fval)
    rv += (msg, {'info': status, 'nfev': nfev})
    return rv

#
################################## PMonCar ##################################
#

class OptErr( Exception ):
    """Base class for exceptions in this module."""
    pass

class InputErr( OptErr ):
    """Exception raised for errors in the input."""
    def __init__( self, x, low, high ):
        self.x = x
        self.low = low
        self.high = high
    def __str__(self):
        msg = 'The input arrays: low (%s), x (%s) and high (%s) sizes must match' % ( self.low, self.x, self.high )
        return msg

class MaxfevErr( OptErr ):
    """Exception raised when the maximum number of function evaluation is exceeded"""
    def __init__( self, maxfev ):
        self.maxfev = maxfev
    def __str__( self ):
        msg = 'number of function evaluations has exceeded maxfev=%d' % self.maxfev
        return msg
    
class OutOfBoundErr( OptErr ):
    """Exception raised for errors in the input (low <= x <= high)"""
    def __init__( self, x, low, high ):
        self.x = x
        self.low = low
        self.high = high
    def __str__( self ):
        msg = 'The following condition must be valid: low (%s) <= x (%s) <= high (%s)' % ( self.low, self.x, self.high )
        return msg

class MyOpt( object ):
    
    def __init__( self, fcn ):
        self.nfev, self.fcn_counter = func_counter( fcn )

    def check_args( self, x, xmin, xmax ):
        x = numpy.array( x, numpy.float_ )
        xmin = numpy.asarray( xmin, numpy.float_ )
        xmax = numpy.asarray( xmax, numpy.float_ )
        if ( x.shape != xmin.shape ) or ( x.shape != xmax.shape ):
            raise InputError( x, xmin, xmax )
        if _outside_limits( x, xmin, xmax ):
            raise OutOfBoundErr( x, xmin, xmax )
        return x, xmin, xmax

class Polytope( MyOpt ):

    debug = True
    
    def __init__( self, fcn, x, xmin, xmax, step, initsimplex ):
        MyOpt.__init__( self, fcn )
        x, self.xmin, self.xmax = self.check_args( x, xmin, xmax )
        self.fcn = func_bounds( self.fcn_counter, xmin, xmax )
        self.polytope = self.init_simplex( x, step, initsimplex )
        
    def __getitem__( self, key ):
        return self.polytope[ key ]
    
    def __setitem__( self, key, item ):
        self.polytope[ key ] = item

    def __str__( self ):
        return self.polytope.__str__( )
    
    def calc_centroid( self, badindex ):
        return numpy.mean( self.polytope[ :badindex, : ], 0 )

    def check_convergence( self, tolerance, method ):
        def are_func_vals_close_enough( tolerance ):
            smallest_fct_val = self.polytope[ 0, -1]
            largest_fct_val = self.polytope[ -1, -1 ]
            return Knuth_close( smallest_fct_val, largest_fct_val,
                                tolerance )

        def is_fct_stddev_small_enough( tolerance ):
            return numpy.std( self.get_func_vals( ) ) < tolerance

        def is_max_length_small_enough( tol ):
            """
                       
               max  || x  - x  || <= tol max( 1.0, || x || )
                        i    0                         0
                        
            where 1 <= i <= n"""

            max_xi_x0 = -1.0
            x0 = self.polytope[ 0, :-1 ]

            npar_plus_1 = len( x0 ) + 1
            for ii in xrange( 1, npar_plus_1 ):
                xi = self.polytope[ ii, :-1 ]
                xi_x0 = xi - x0
                max_xi_x0 = max( max_xi_x0, numpy.dot( xi_x0, xi_x0 ) )
            return max_xi_x0 <= tol * max( 1.0, numpy.dot( x0, x0 ) )

        if 0 == method:
            if is_max_length_small_enough( tolerance ):
                return True
        elif 2 == method:
            if False == is_max_length_small_enough( tolerance ):
                return False
            stddev = is_fct_stddev_small_enough( tolerance )
            fctval = are_func_vals_close_enough( tolerance )
            return stddev and fctval
        else:
            if False == is_max_length_small_enough( tolerance ):
                return False
            stddev = is_fct_stddev_small_enough( tolerance )
            fctval = are_func_vals_close_enough( tolerance )
            return stddev or fctval
        return False

        num = 2.0 * abs( self.polytope[ 0, -1 ] - self.polytope[ -1, -1 ] )
        denom = abs( self.polytope[ 0, -1 ] ) + abs( self.polytope[ -1, -1 ] ) + 1.0
        if ( num / denom > tolerance ):
            return False
        func_vals = self.get_func_vals( )
        if numpy.std( func_vals ) > tolerance:
            return False
        return True

    # 4.
    # return of true ==> perform step 5 (shrink) otherwise do not shrink
    def contract_in_out( self, centroid, reflection_pt, rho_gamma,
                         contraction_coef, badindex, maxfev, verbose ):

        if self.polytope[ badindex - 1, -1 ] <= reflection_pt[ -1 ] and \
               reflection_pt[ -1 ] < self.polytope[ badindex, -1 ]:

            # Perform an outside contraction
	    #
	    # 4a. Outside.  If f  <= f  < f    (i.e., x  is stricly better then
	    #                   n     r    n+1         r
	    # x    ), perform an outside contraction: calculate
	    #  n+1
	    #      _                 _                       _
	    # x  = x + gamma * ( x - x ) = ( 1 + rho gamma ) x  - 
	    #  c                  r
            #                                                   rho gamma x
	    #                                                               n+1
	    #                                                             (2.6)
	    #
	    # and evaluate f  = f( x  ).
	    #               c       c
            outside_contraction_pt = self.move_vertex( centroid, rho_gamma,
                                                       badindex )

            # If f  <= f  , accept x  
            #     c     r           c
            if outside_contraction_pt[ -1 ] <= reflection_pt[ -1 ]:
                #self.polytope[ badindex ] = outside_contraction_pt[ : ]
                self.replace_vertex( badindex, outside_contraction_pt, maxfev )
                if verbose:
                    print '\taccept outside contraction point'
                    # and terminate the iteration
                return False
            else:
                # otherwise, go to step 5 (perform a shrink).
                return True
            
        elif reflection_pt[ -1 ] >= self.polytope[ badindex, -1 ]:

            # Perform an inside contraction
	    #
	    # 4b. Inside. If f  >= f   , perform an inside contraction: calculate
	    #                 r     n+1
	    #       _           _                         _
	    # x   = x - gamma ( x - x   ) = ( 1 - gamma ) x + gamma x     (2.7)
	    #  cc                    n+1                             n+1
	    #
	    # and evaluate f   = f( x   ).
	    #               cc       cc
	    #
            inside_contraction_pt = self.move_vertex( centroid,
                                                      - contraction_coef,
                                                      badindex )
                    
            if inside_contraction_pt[ -1 ] < self.polytope[ badindex, -1 ]:
                #
                # If f   < f   , accept x  
                #     cc    n+1          cc
                #
                self.replace_vertex( badindex, inside_contraction_pt, maxfev )
                if verbose:
                    print '\taccept inside contraction point'
                return False
            else:
                # otherwise, go to step 5 (perform a shrink).
                return True

        else:
            print 'something is wrong with contract_in_out'
            return True

    def get_func_vals( self ):
        return self.polytope[ :, -1 ]

    def order( self ):
        func_vals = self.get_func_vals( )
        # return a list of indices goes from largest to smallest
        index_from_largest_to_smallest = numpy.argsort( func_vals )
        # re-arrange so self.polytope goes from smallest to largest fct val
        self.polytope = numpy.take( self.polytope,
                                   index_from_largest_to_smallest, 0 )

    def init_simplex( self, x, step, initsimplex ):
        npar = len( x )
        npar_plus_1 = npar + 1
        simplex = numpy.zeros( (npar_plus_1,npar_plus_1), dtype=numpy.float_ )
        if step is None:
            # step = 1.2 * numpy.ones(x.shape, numpy.float_, numpy.isfortran(x))
            step = map( lambda fubar: 1.2 * fubar + 1.2, x )
        for ii in xrange( npar_plus_1 ):
            simplex[ ii , :-1 ] = numpy.array( x, copy=True )
        if 0 == initsimplex:
            for ii in xrange( npar ):
                simplex[ ii + 1 ,  ii ] += step[ ii ]
        else:
            delta = 0.1
            for ii in xrange( npar ):
                if 0.0 == simplex[ ii + 1 ,  ii ]:
                    simplex[ ii + 1 ,  ii ] = delta
                else:
                    simplex[ ii + 1 ,  ii ] *= ( 1.0 + delta )
        for ii in xrange( npar + 1 ):
            simplex[ ii , -1 ] = self.fcn( simplex[ ii , :-1 ] )
        return simplex

    def move_vertex( self, centroid, coef, badindex ):
        vertex = ( 1.0 + coef ) * centroid - coef * self.polytope[ badindex ]
        vertex[ -1 ] = self.fcn( vertex[ :-1 ] )
        return vertex

    def replace_vertex( self, badindex, newvertex, maxfev ):
        self.polytope[ badindex ] = newvertex[ : ]
        
    def shrink( self, shrink_coef, npar ):
        npars_plus_1 = npar + 1
        for ii in xrange( 1, npars_plus_1 ):
            self.polytope[ ii ] = self.polytope[ 0 ] + shrink_coef * ( self.polytope[ ii ] - self.polytope[ 0 ] )
            self.polytope[ ii, -1 ] = self.fcn( self.polytope[ ii, :-1 ] )


class DirectSearch( MyOpt ):

    def __init__( self, fcn ):
        MyOpt.__init__( self, fcn )
        self.expansion_coef   = 2.0          # chi
        self.contraction_coef = 0.5          # gamma
        self.reflection_coef  = 1.0          # rho
        self.shrink_coef      = 0.5          # sigma

class MyNelderMead( DirectSearch ):

    def __init__( self, fcn ):
        DirectSearch.__init__( self, fcn )
        
    def __call__( self, x, xmin, xmax, maxfev, tol, step, initsimplex,
                  finalsimplex, verbose, algorithm ):

        def optimize( polytope, badindex ):

            rho_chi = self.reflection_coef * self.expansion_coef
            rho_gamma = self.reflection_coef * self.contraction_coef
            while polytope.nfev[ 0 ] < maxfev:

                polytope.order( )
                if verbose:
                    print 'f%s=%e' % (polytope[ 0, :-1 ],polytope[ 0, -1 ])

                centroid = polytope.calc_centroid( badindex )

                
##                 print polytope
##                 print polytope.check_convergence( tol, finalsimplex ), 
                
                if polytope.check_convergence( tol, finalsimplex ):
                    break;

                reflection_pt = polytope.move_vertex( centroid,
                                                      self.reflection_coef,
                                                      badindex )

                if polytope[ 0, -1 ] <= reflection_pt[ -1 ] and \
                       reflection_pt[ -1 ] < polytope[ badindex - 1, -1 ]:
                    #
                    # If f  <= f  < f  ,
                    #     1     r    n
                    polytope.replace_vertex( badindex, reflection_pt, maxfev )
                    if verbose:
                        print '\taccept reflection point'

                elif reflection_pt[ -1 ] < polytope[ 0, -1 ]:
                    # 3. Expand. If f  < f  ,
                    #                r    1  
                    #
                    # calculate the expansion point x  :
                    #                                e
                    #      _               _                      _
                    # x  = x + chi * ( x - x ) =  ( 1 + rho chi ) x  - rho chi x      (2.5)
                    #  e                r                                       n+1
                    #
                    # and evaluate f  = f( x )
                    #               e       e
                    #
                    expansion_pt = polytope.move_vertex( centroid, rho_chi,
                                                         badindex )

                    if expansion_pt[ -1 ] < reflection_pt[ -1 ]:
                        #
                        # If f  < f  , accept x  and terminate the iteration;
                        #     e    r           e
                        polytope.replace_vertex( badindex, expansion_pt,
                                                 maxfev )
                        if verbose:
                            print '\taccept expansion point'
                    else:
                        #
                        # otherwise, (if f  >= f  ), accept x  and terminate the iteration.
                        #                 e     r            r
                        polytope.replace_vertex( badindex, reflection_pt,
                                                 maxfev )
                        if verbose:
                            print '\taccept reflection point'

                else: 

                    #
                    # 4. Contract. If f  >= f  perform a contraction between 
                    #                  r     n
                    # 
                    # centroid and the better of x    and x .
                    #                             n+1      r
                    #

                    if polytope.contract_in_out( centroid, reflection_pt,
                                                 rho_gamma,
                                                 self.contraction_coef,
                                                 badindex, maxfev, verbose ):
                        polytope.shrink( self.shrink_coef, len(x) )

            best_vertex = polytope[ 0 ]
            best_par = best_vertex[ :-1 ]
            best_val = best_vertex[ -1 ]
            return best_par, best_val, polytope.nfev[ 0 ]


        bad_index = -1
        myfcn = func_bounds( self.fcn_counter, xmin, xmax )
        polytope = algorithm( myfcn, x, xmin, xmax, step, initsimplex )
        return optimize( polytope, bad_index )


class PowellPolytope( Polytope ):

    def __init__( self, fcn, x, xmin, xmax, step, initsimplex ):
        Polytope.__init__( self, fcn, x, xmin, xmax, step, initsimplex )

    def replace_vertex( self, badindex, newvertex, maxfev ):
        #print newvertex[ -1 ], '\t', self.polytope[ 0, -1 ]
        if newvertex[ -1 ] < self.polytope[ 0, -1 ] and len( newvertex ) > 2 and sao_fcmp( newvertex[ -1 ], self.polytope[ 0, -1 ], 1.0e-2 ) < 0:
            result = newuoa( self.fcn, newvertex[ :-1 ], self.xmin, self.xmax,
                             maxfev = maxfev - self.nfev[ 0 ] )
            if self.debug:
                print '%e\tvs\t%e: %e @ %d nfevs' % ( newvertex[ -1 ], 
                                                      self.polytope[ 0, -1 ],
                                                      result[ 2 ],
                                                      result[ 4 ].get( 'nfev' ) )
                print sao_fcmp( newvertex[ -1 ], self.polytope[ 0, -1 ], 1.0e-1 ), Knuth_close( newvertex[ -1 ], self.polytope[ 0, -1 ], 1.0e-1 )
                                
            newvertex[ :-1 ] = numpy.asarray( result[ 1 ], numpy.float_ )
            newvertex[ -1 ] = result[ 2 ]
        self.polytope[ badindex ] = newvertex[ : ]

def nelder_mead( fcn, x0, xmin, xmax, step=None, initsimplex=0, finalsimplex=2,
                 algorithm=Polytope, maxfev=None, ftol=EPSILON,
                 verbose=0 ):

    def stat_cb0( pars ):
        return fcn( pars )[ 0 ]

    print 'nelder_mead: fcn(x0) = ', fcn( x0 )
    
    x0, xmin, xmax = _check_args(x0, xmin, xmax)

    nm = MyNelderMead( stat_cb0 )

    if maxfev is None:
        maxfev = 1024 * len( x0 )

    def my_nelder_mead( myx, mymaxfev, mystep, myfinalsimplex,
                        debug=True ):

        ofval = FUNC_MAX
        while nm.nfev[ 0 ] < mymaxfev:
            x, fval, nfev = nm( myx, xmin, xmax, mymaxfev, ftol, mystep,
                                initsimplex, finalsimplex, verbose, algorithm )
        
            if debug:
                print 'nelder_mead: f%s=%e at %d' % ( x, fval, nm.nfev[ 0 ] )
            if sao_fcmp( ofval, fval, ftol ) <= 0:
                return x, fval, nm.nfev[ 0 ]
            ofval = fval

        return x, fval, nfev
            
    myfinalsimplex = 2

    if step is None:
        step = map( lambda fubar: 1.2 * fubar + 1.2, x0 )

    x, fval, nfev = my_nelder_mead( x0, maxfev, step, myfinalsimplex )
    
    ierr = 0
    if nfev >= maxfev:
        ierr = 3

    status, msg = _get_saofit_msg( maxfev, ierr )
    rv = (status, x, fval)
    rv += (msg, {'info': status, 'nfev': nfev})
    return rv
    
    
    
#
############################ Differential Evolution ###########################
#

class DifferentialEvolution( MyOpt ):

    def __init__( self, fcn, seed ):
        MyOpt.__init__( self, fcn )
        random.seed( seed )

    def init_population( self, x0, xmin, xmax, population_size ):
        npar = len( xmin )
        population = numpy.zeros( (population_size, npar + 1),
                                  dtype=numpy.float_ )
        for i in xrange( population_size ):
            for j in xrange( npar ):
                population[ i, j ] = random.uniform( xmin[ j ], xmax[ j ] )
            population[ i, -1 ] = FUNC_MAX
        return population

    def __call__( self, x0, xmin, xmax, strategy, scale, xprob, maxfev, tol,
                  population_size, verbose ):

        npar = len( xmin )

        model_par = numpy.zeros( ( npar + 1, ), dtype=numpy.float_ )
        trial_solution = numpy.zeros( ( npar + 1, ), dtype=numpy.float_ )
        population = self.init_population( x0, xmin, xmax, population_size )
        x0, xmin, xmax = self.check_args( x0, xmin, xmax )
            
        fcn = func_bounds( self.fcn_counter, xmin, xmax )

        model_par[ :-1 ] = x0[ : ]
        model_par[ -1 ] = fcn( model_par[ :-1 ] )

        if 0 == strategy:
            strategy_fcn = self.Best1Exp
        elif 1== strategy:
            strategy_fcn = self.Rand1Exp
        elif 2 == strategy:
            strategy_fcn = self.RandToBest1Exp
        elif 3 == strategy:
            strategy_fcn = self.Best2Exp
        elif 4 == strategy:
            strategy_fcn = self.Rand2Exp
        elif 5 == strategy:
            strategy_fcn = self.Best1Bin
        elif 6 == strategy:
            strategy_fcn = self.Rand1Bin
        elif 7 == strategy:
            strategy_fcn = self.RandToBest1Bin
        elif 8 == strategy:
            strategy_fcn = self.Best2Bin
        elif 9 == strategy:
            strategy_fcn = self.Rand2Bin
        else:
            strategy_fcn = self.Best1Exp

        while self.nfev[ 0 ] < maxfev:

            for pop in xrange( population_size ):

                if numpy.std( population[ :, -1 ] ) < tol or \
                       self.nfev[ 0 ] >= maxfev:
                    return model_par[ :-1 ], model_par[ -1 ], self.nfev[ 0 ]
                
                trial_solution = strategy_fcn( pop, population, model_par,
                                               population[ pop ].copy(),
                                               scale, xprob )
                
                trial_solution[ -1 ] = fcn( trial_solution[ :-1 ] )

                if trial_solution[ -1 ] < population[ pop, -1 ]:
                    population[ pop ] = trial_solution.copy( )

                    if trial_solution[ -1 ] < model_par[ -1 ]:

                        mystep = map( lambda fubar: 1.2 * fubar + 1.2,
                                      trial_solution[ :-1 ] )
                        result = nelder_mead( fcn, trial_solution[ :-1 ], xmin,
                                              xmax, step=mystep,
                                              maxfev=maxfev - self.nfev[ 0 ] )
                        model_par[ :-1 ] = numpy.asarray( result[ 1 ],
                                                          numpy.float_ )
                        model_par[ -1 ] = result[ 2 ]
                        if verbose:
                            print 'f(', model_par[ :-1 ], ')=', \
                                  model_par[ -1 ], '\t', self.nfev[ 0 ]

        return model_par[ :-1 ], model_par[ -1 ], self.nfev[ 0 ]

    def select_candidate( self, population_size, pop, num ):
        return random.sample( range( pop ) +
                              range( pop + 1, population_size ), num )
        
    def Best1Exp( self, pop, population, model_par, trial_solution,
                  scale, xprob ):

        npar = population.shape[ 1 ] - 1
        n = random.randint( 0, npar - 1 )
        r1, r2 = self.select_candidate( population.shape[ 0 ], pop, 2 )

        while True:
            for i in xrange( npar ):
                if random.random( ) >= xprob:
                    return trial_solution
                trial_solution[ n ] = model_par[ n ] \
                                      + scale * ( population[ r1, n ] - \
                                                  population[ r2, n ] )
                n = (n + 1) % npar

        return trial_solution            


    def Rand1Exp( self, pop, population, model_par, trial_solution,
                  scale, xprob ):

        npar = population.shape[ 1 ] - 1
        n = random.randint( 0, npar - 1 )
        r1, r2, r3 = self.select_candidate( population.shape[ 0 ], pop,
                                            3 )

        while True:
            for i in xrange( npar ):
                if random.random( ) >= xprob:
                    return trial_solution
                trial_solution[ n ] = population[ r1, n ] + \
                                      scale * ( population[ r2, n ] - \
                                                population[ r3, n ] )
                n = (n + 1) % npar
            
        return trial_solution            


    def RandToBest1Exp( self, pop, population, model_par,
                        trial_solution, scale, xprob ):

        npar = population.shape[ 1 ] - 1
        n = random.randint( 0, npar - 1 )
        r1, r2 = self.select_candidate( population.shape[ 0 ], pop, 2 )

        while True:
            for i in xrange( npar ):
                if random.random( ) >= xprob:
                    return trial_solution
                trial_solution[ n ] += scale * \
                                       ( model_par[ n ] - \
                                         trial_solution[n] ) + \
                                         scale * ( population[ r1, n ] - \
                                                   population[ r2, n ] )
                n = (n + 1) % npar

        return trial_solution            


    def Best2Exp( self, pop, population, model_par, trial_solution,
                  scale, xprob ):

        npar = population.shape[ 1 ] - 1
        n = random.randint( 0, npar - 1 )
        r1, r2, r3, r4 = self.select_candidate( population.shape[ 0 ],
                                                pop, 4 )

        while True:
            for i in xrange( npar ):
                if random.random( ) >= xprob:
                    return trial_solution
                trial_solution[ n ] = model_par[ n ] + \
                                      scale * ( population[ r1, n ] + \
                                                population[ r2, n ] - \
                                                population[ r3, n ] - \
                                                population[ r4, n ] )
                n = (n + 1) % npar

        return trial_solution            


    def Rand2Exp( self, pop, population, model_par, trial_solution,
                  scale, xprob ):

        npar = population.shape[ 1 ] - 1
        n = random.randint( 0, npar - 1 )
        r1, r2, r3, r4, r5 = self.select_candidate( population.shape[ 0 ],
                                                    pop, 5 )

        while True:
            for i in xrange( npar ):
                if random.random( ) >= xprob:
                    return trial_solution
                trial_solution[ n ] = population[ r1, n ] + \
                                      scale * ( population[ r2, n ] + \
                                                population[ r3, n ] - \
                                                population[ r4, n ] - \
                                                population[ r5, n ] )
                n = (n + 1) % npar
            
        return trial_solution            


    def Best1Bin( self, pop, population, model_par, trial_solution,
                  scale, xprob ):

        npar = population.shape[ 1 ] - 1
        n = random.randint( 0, npar - 1 )
        r1, r2 = self.select_candidate( population.shape[ 0 ], pop, 2 )

        while True:
            for i in xrange( npar ):
                if random.random( ) >= xprob:
                    return trial_solution
                trial_solution[ n ] = model_par[ n ] + \
                                      scale * ( population[ r1, n ] - \
                                                population[ r2, n ] )
                n = (n + 1) % npar
            
        return trial_solution


    def Rand1Bin( self, pop, population, model_par, trial_solution,
                  scale, xprob ):

        npar = population.shape[ 1 ] - 1
        n = random.randint( 0, npar - 1 )
        r1, r2, r3 = self.select_candidate( population.shape[ 0 ], pop,
                                            3 )

        while True:
            for i in xrange( npar ):
                if random.random( ) >= xprob:
                    return trial_solution
                trial_solution[ n ] = population[ r1, n ] + \
                                      scale * ( population[ r2, n ] -\
                                                population[ r3, n ] )
                n = (n + 1) % npar
            
        return trial_solution            


    def RandToBest1Bin( self, pop, population, model_par,
                        trial_solution, scale, xprob ):

        npar = population.shape[ 1 ] - 1
        n = random.randint( 0, npar - 1 )
        r1, r2 = self.select_candidate( population.shape[ 0 ], pop, 2 )

        while True:
            for i in xrange( npar ):
                if random.random( ) >= xprob:
                    return trial_solution
                trial_solution[ n ] += scale * \
                                       ( model_par[ n ] -
                                         trial_solution[n]) + \
                                         scale * ( population[ r1, n ] - \
                                                   population[ r2, n ] )
                n = (n + 1) % npar
            
        return trial_solution            


    def Best2Bin( self, pop, population, model_par, trial_solution,
                  scale, xprob ):

        npar = population.shape[ 1 ] - 1
        n = random.randint( 0, npar - 1 )
        r1, r2, r3, r4 = self.select_candidate( population.shape[ 0 ],
                                                pop, 4 )

        while True:
            for i in xrange( npar ):
                if random.random( ) >= xprob:
                    return trial_solution
                trial_solution[ n ] = model_par[ n ] + \
                                      scale * ( population[ r1, n ] + \
                                                population[ r2, n ] - \
                                                population[ r3, n ] - \
                                                population[ r4, n ] )
                n = (n + 1) % npar
            
        return trial_solution            


    def Rand2Bin( self, pop, population, model_par, trial_solution,
                  scale, xprob ):

        npar = population.shape[ 1 ] - 1
        n = random.randint( 0, npar - 1 )
        r1, r2, r3, r4, r5 = self.select_candidate( population.shape[ 0 ],
                                                    pop, 5 )

        while True:
            for i in xrange( npar ):
                if random.random( ) >= xprob:
                    return trial_solution
                trial_solution[ n ] = population[ r1, n ] + \
                                      scale * ( population[ r2, n ] + \
                                                population[ r3, n ] - \
                                                population[ r4, n ] - \
                                                population[ r5, n ] )
                n = (n + 1) % npar
            
        return trial_solution


def dif_evo( fcn, x0, xmin, xmax, population_size=None, scale=0.5, strategy=0,
             xprob=0.9, seed=351, ftol=EPSILON, maxfev=None, verbose=False ):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    def stat_cb0( pars ):
        return fcn( pars )[ 0 ]

    de = DifferentialEvolution( stat_cb0, seed )

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = max( 0.1, xprob )
    xprob = min( xprob, 1.0 )

    # make sure that scale is within [0.1,1.0]
    scale = max( 0.1, scale )
    scale = min( scale, 1.0 )

    if maxfev is None:
        maxfev = 8192 * x.size
        
    if population_size is None:
        population_size = 12 * x.size

    try:
        ierr = 0

        maxfev_per_iter = 512 * x.size

        def my_dif_evo( myx, myxmin, myxmax, mystrategy, mymaxfev,
                        limit_factor=4.0, debug=True ):

            ofval = FUNC_MAX
            de_nfev = 0
            while de_nfev < mymaxfev:
                
                myxmin, myxmax = _narrow_limits( limit_factor,
                                                 [ myx, myxmin, myxmax ],
                                                 debug = False )

                print 'myxmin = ', myxmin
                print 'myxmax = ', myxmax

                myx, fval, tmp_nfev = de( myx, myxmin, myxmax, strategy, scale,
                                          xprob, maxfev_per_iter, ftol,
                                          population_size, verbose )
                de_nfev += tmp_nfev

                if debug:
                    print 'de: f%s=%.14e at %d nfevs' % ( myx, fval, de_nfev )

                if sao_fcmp( ofval, fval, ftol ) <= 0:
                    return myx, fval, de_nfev
                ofval = fval
                limit_factor *= 2

            return myx, fval, de_nfev                

        mystep = map( lambda fubar: 1.2 * fubar + 1.2, x0 )
        result = neldermead( fcn, x0, xmin, xmax, maxfev=maxfev,
                             ftol=ftol, finalsimplex=2, step=mystep )
        myx = numpy.asarray( result[ 1 ], numpy.float_ )
        nm_nfev = result[4].get( 'nfev' )

        xpar, fval, de_nfev = my_dif_evo( myx, xmin, xmax, strategy,
                                          maxfev - nm_nfev, limit_factor=8.0 )

        result = neldermead( fcn, xpar, xmin, xmax,
                             maxfev=maxfev-de_nfev-nm_nfev,
                             ftol=ftol, finalsimplex=2, step=mystep )
        xpar = numpy.asarray( result[ 1 ], numpy.float_ )
        nm_nfev += result[4].get( 'nfev' )

        
    except InputErr:
        ierr = 1
    except OutOfBoundErr:
        ierr = 2

    status, msg = _get_saofit_msg( maxfev, ierr )
    rv = ( status, xpar, fval )
    rv += (msg, {'info': ierr, 'nfev': de_nfev + nm_nfev})

    print rv
    return rv
    
    
#
############################ Differential Evolution ###########################
#


###############################################################################
## class Bobyqa(OptMethod):
##     def __init__(self, name='bobyqa'):
##         OptMethod.__init__(self, name, bobyqa)

## class Dif_Evo(OptMethod):
##     def __init__(self, name='dif_evo'):
##         OptMethod.__init__(self, name, dif_evo)

## class MarLev(OptMethod):
##     def __init__(self, name='marlev'):
##         OptMethod.__init__(self, name, marquadt_levenberg, True)

## class mcMonCar(OptMethod):
##     def __init__(self, name='mcMonCar'):
##         OptMethod.__init__(self, name, mc_montecarlo)

## class Nelder_Mead(OptMethod):
##     def __init__(self, name='nelder_mead'):
##         OptMethod.__init__(self, name, nelder_mead)

## class Newuoa(OptMethod):
##     def __init__(self, name='newuoa'):
##         OptMethod.__init__(self, name, newuoa)

## class OdrPack(OptMethod):
##     def __init__(self, name='odrpack'):
##         OptMethod.__init__(self, name, odrpack, True)

## class PortChi(OptMethod):
##     def __init__(self, name='dn2fb'):
##         OptMethod.__init__(self, name, dn2fb, True)

## class PortFct(OptMethod):
##     def __init__(self, name='dmnfb'):
##         OptMethod.__init__(self, name, dmnfb)
    
## ui._session._methods['bobyqa'] = Bobyqa( )
## ui._session._methods['dif_evo'] = Dif_Evo( )
## ui._session._methods['portchi'] = PortChi( )
## ui._session._methods['portfct'] = PortFct( )
## ui._session._methods['mcmoncar'] = mcMonCar( )
## ui._session._methods['newuoa'] = Newuoa( )
## ui._session._methods['odrpack'] = OdrPack( )
###############################################################################

if '__main__' == __name__:

    def print_result( result, name ):
        x = result[1]
        f = result[2]
        nfev = result[4].get( 'nfev' )
        print '%s%s = %f in %d nfevs' % ( name, x, f, nfev )

    def tst_odr( ):
        def Rosenbrock_fvec(x):
            x = numpy.asarray(x)
            f1 = 1.0-x[0]
            f2 = 10.0*(x[1]-x[0]*x[0])
            f3 = 1.0-x[2]
            f4 = 10.0*(x[3]-x[2]*x[2])
            return [f1,f2,f3,f4]

        def Brown_badly_scaled_fvec(x):
            fvec0 = x[ 0 ] - 1.0e6;
            fvec1 = x[ 1 ] - 2.0e-6;
            fvec2 = x[ 0 ] * x[ 1 ] - 2.0;
            return [fvec0, fvec1, fvec2]

        x0 = [-1.2, 1.0, -1.2, 1.0]
        xmin = [-10.0, -10.0, -10.0, -10.0]
        xmax = [10.0, 10.0, 10.0, 10.0]
        print_result( odrpack( Rosenbrock_fvec, x0, xmin, xmax ), 'odrpack' )

        x0 = [ 1.0, 1.0 ]
        xmin = [-1.0e8, -1.0e8]
        xmax = [1.0e8, 1.0e8]
        print_result( odrpack( Brown_badly_scaled_fvec, x0, xmin, xmax ),
                      'odrpack' )


    def Rosenbrock(x):  # The Rosenbrock function
        #print 'hi'
        x = numpy.asarray(x)
        val = numpy.sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0,axis=0)
        #print 'f%s=%f' % ( x, val )
        return val

    def FreudensteinRoth( x ):
        npar = len( x )
        sum = 0.0
        for ii in xrange( 0, npar, 2 ):
            fx0 = - 13.0 + x[ ii ] + x[ ii + 1 ] * (( 5.0 - x[ ii + 1 ] ) \
                                                    * x[ ii + 1 ] - 2.0 )
            fx1 = - 29.0 + x[ ii ] + x[ ii + 1 ] * (( x[ ii + 1 ] + 1.0 ) \
                                                    * x[ ii + 1 ] - 14.0 )
            sum += fx0*fx0 + fx1*fx1
        return sum

    def tst_nm( npar ):
        x0 = npar * [ -1.2, 1.0 ]
        xmin = npar * [ -10, -10. ]
        xmax = npar * [ 10, 10 ]
        maxfev = None
        ftol = EPSILON
        verbose = 0
        step = 4.0 
        initsimplex = 1
        finalsimplex = 0
        #algo = PowellPolytope
        algo = Polytope
        print nelder_mead( Rosenbrock, x0, xmin, xmax, step=step,
                           initsimplex=initsimplex, finalsimplex=finalsimplex,
                           algorithm=algo, maxfev=maxfev, ftol=ftol,
                           verbose=verbose )

    def tst_multicore_montecarlo( npar ):
        x0 = npar * [ -1.2, 1.0 ]
        xmin = npar * [ -10, -10. ]
        xmax = npar * [ 10, 10 ]
        print mc_montecarlo( Rosenbrock, x0, xmin, xmax, verbose=0 )
        
    def tst( algorithm, npar=4 ):
        x0 = npar * [ -1.2, 1.0 ]
        xmin = npar * [ -10, -10. ]
        xmax = npar * [ 10, 10 ]
        result = algorithm( Rosenbrock, x0, xmin, xmax )
        print_result( result, algorithm.__name__ )

        x0 = npar * [ 0.5, -2.0 ]
        xmin = npar * [ -20, -20. ]
        xmax = npar * [ 20., 20. ]
        result = algorithm( FreudensteinRoth, x0, xmin, xmax )
        print_result( result, algorithm.__name__ )

    def tst_de( npar ):
        import time
        x0 = npar * [ -1.2, 1.0 ]
        xmin = npar * [ -10, -10. ]
        xmax = npar * [ 10, 10 ]
        maxfev = 8192 * (npar * 2)
        tolerance = 1.0e-7
        scale_factor = 0.7
        xprob = 0.5
        population_size = npar * 32
        tol = 1.0e-7
        verbose = False

        for strategy in range( 1 ):
            start_time = time.time()
            print dif_evo( Rosenbrock, x0, xmin, xmax, strategy=strategy );
            print 'strategy = %d\t%f secs' % ( strategy, time.time() - start_time )
        
    
    import sys
    npar = 8
    if ( 2 == len(sys.argv) ):
        npar = int( sys.argv[1] )
##     tst( bobyqa, npar )
##     tst( newuoa, npar )
##     #tst( dmnfb, npar )
##     tst_odr( )
##     tst_nm( npar )
##     tst_multicore_montecarlo( npar )
    tst_de( npar )
