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
import random
import sys
from itertools import izip


import _minpack
import _minim
import _saoopt

from sherpa.utils import parallel_map
from sherpa.utils._utils import sao_fcmp

#
# Use FLT_EPSILON as default tolerance
#
EPSILON = numpy.float_(numpy.finfo(numpy.float32).eps)

#
# Maximum callback function value, used to indicate that the optimizer
# has exceeded parameter boundaries.  All the optimizers expect double
# precision arguments, so we use numpy.float_ instead of SherpaFloat.
#
# As of numpy 0.9.5, 'max' is a 0-D array, which seems like a bug.
#
FUNC_MAX = numpy.float_(numpy.finfo(numpy.float_).max)

def _check_args(x0, xmin, xmax):
    x = numpy.array(x0, numpy.float_)  # Make a copy
    xmin = numpy.asarray(xmin, numpy.float_)
    xmax = numpy.asarray(xmax, numpy.float_)

    if (x.shape != xmin.shape) or (x.shape != xmax.shape):
        raise TypeError('input array sizes do not match')

    _move_within_limits(x, xmin, xmax)
    
    return x, xmin, xmax

def _get_saofit_msg( maxfev, ierr ):
    key = {
        0: (True, 'successful termination'),
        1: (False, 'improper input parameters'),
        2: (False, 'initial parameter value is out of bounds'),
        3: (False,
            ('number of function evaluations has exceeded maxfev=%d' %
             maxfev))
        }
    return key.get( ierr, (False, 'unknown status flag (%d)' % ierr))

def _move_within_limits(x, xmin, xmax):
    below = numpy.flatnonzero(x < xmin)
    if below.size > 0:
        x[below] = xmin[below]

    above = numpy.flatnonzero(x > xmax)
    if above.size > 0:
        x[above] = xmax[above]

def _my_is_nan( x ):
    fubar = filter( lambda xx: xx != xx or xx is numpy.nan or numpy.isnan(xx) and numpy.isfinite(xx), x )
    if len( fubar ) > 0:
        return True
    else:
        return False

def _narrow_limits( myrange, xxx, debug ):

    def double_check_limits( myx, myxmin, myxmax ):
        for my_l,my_x,my_h in izip( myxmin, myx, myxmax ):
            if my_x < my_l:
                print 'x = ', my_x, ' is < lower limit = ', my_l
            if my_x > my_h:
                print 'x = ', my_x, ' is > upper limit = ', my_h

    def raise_min_limit( range, xmin, x, debug=False ):
        myxmin = numpy.asarray( map( lambda xx: xx - range * numpy.abs(xx),
                                     x ), numpy.float_ )
        if False != debug:
            print
            print 'raise_min_limit: myxmin=%s' % myxmin    
            print 'raise_min_limit: x=%s' % x
        below = numpy.flatnonzero(myxmin < xmin)
        if below.size > 0:
            myxmin[below] = xmin[below]
        if False != debug:
            print 'raise_min_limit: myxmin=%s' % myxmin    
            print 'raise_min_limit: x=%s' % x
            print
        return myxmin

    def lower_max_limit( range, x, xmax, debug=False ):
        myxmax = numpy.asarray( map( lambda xx: xx + range * numpy.abs(xx),
                                     x ), numpy.float_ )
        if False != debug:
            print
            print 'lower_max_limit: x=%s' % x
            print 'lower_max_limit: myxmax=%s' % myxmax    
        above = numpy.flatnonzero(myxmax > xmax)
        if above.size > 0:
            myxmax[above] = xmax[above]
        if False != debug:
            print 'lower_max_limit: x=%s' % x
            print 'lower_max_limit: myxmax=%s' % myxmax    
            print
        return myxmax

    x = xxx[ 0 ]
    xmin = xxx[ 1 ]
    xmax = xxx[ 2 ]

    if False != debug:
        print 'narrow_limits: xmin=%s' % xmin
        print 'narrow_limits: x=%s' % x        
        print 'narrow_limits: xmax=%s' % xmax
    myxmin = raise_min_limit( myrange, xmin, x, debug=False )
    myxmax = lower_max_limit( myrange, x, xmax, debug=False )

    if False != debug:
        print 'range = %d' % myrange
        print 'narrow_limits: myxmin=%s' % myxmin
        print 'narrow_limits: x=%s' % x        
        print 'narrow_limits: myxmax=%s\n' % myxmax

    double_check_limits( x, myxmin, myxmax )
            
    return myxmin, myxmax


def _outside_limits(x, xmin, xmax):
    return (numpy.any(x < xmin) or numpy.any(x > xmax))

def _same_par(a,b):
    b = numpy.array(b, numpy.float_)
    same = numpy.flatnonzero( a < b )
    if same.size == 0:
        return 1
    return 0

def _set_limits(x, xmin, xmax):
    below = numpy.nonzero(x < xmin)
    if below.size > 0:
        return 1

    above = numpy.nonzero(x > xmax)
    if above.size > 0:
        return 1

    return 0


__all__ = ('difevo', 'difevo_lm', 'difevo_nm', 'grid_search', 'lmdif', 'minim', 'montecarlo', 'neldermead')

def difevo(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0,
           seed=2005815, population_size=None, xprob=0.9,
           weighting_factor=0.8):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = max( 0.1, xprob )
    xprob = min( xprob, 1.0 )

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = max( 0.1, weighting_factor )
    weighting_factor = min( weighting_factor, 1.0 )

    if population_size is None:
        population_size = 16 * x.size

    if maxfev is None:
        maxfev = 1024 * x.size

    de = _saoopt.difevo( verbose, maxfev, seed, population_size, ftol, xprob,
                         weighting_factor, xmin, xmax, x, fcn )
    fval = de[ 1 ]
    nfev = de[ 2 ]
    ierr = de[ 3 ]

    if verbose:
        print 'difevo: f%s=%e in %d nfev' % ( x, fval, nfev )
    
    status, msg = _get_saofit_msg( maxfev, ierr )
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev})
        
    return rv

def difevo_lm(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0,
              seed=2005815, population_size=None, xprob=0.9,
              weighting_factor=0.8):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = max( 0.1, xprob )
    xprob = min( xprob, 1.0 )

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = max( 0.1, weighting_factor )
    weighting_factor = min( weighting_factor, 1.0 )

    if population_size is None:
        population_size = 16 * x.size

    if maxfev is None:
        maxfev = 1024 * x.size

    de = _saoopt.lm_difevo( verbose, maxfev, seed, population_size, ftol,
                            xprob, weighting_factor, xmin, xmax,
                            x, fcn, numpy.asanyarray(fcn(x)).size )
    fval = de[ 1 ]
    nfev = de[ 2 ]
    ierr = de[ 3 ]

    status, msg = _get_saofit_msg( maxfev, ierr )
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev})
        
    return rv

def difevo_nm(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0,
               seed=741985, population_size=None, xprob=0.9,
               weighting_factor=0.8):

    def stat_cb0( pars ):
        return fcn( pars )[ 0 ]

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = max( 0.1, xprob )
    xprob = min( xprob, 1.0 )

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = max( 0.1, weighting_factor )
    weighting_factor = min( weighting_factor, 1.0 )

    if population_size is None:
        population_size = max( population_size, 16 * x.size )

    if maxfev is None:
        maxfev = 1024 * population_size

    de = _saoopt.nm_difevo( verbose, maxfev, seed, population_size,
                            ftol, xprob, weighting_factor, xmin, xmax,
                            x, stat_cb0 )
    fval = de[ 1 ]
    nfev = de[ 2 ]
    ierr = de[ 3 ]

    if verbose:
        print 'difevo_nm: f%s=%e in %d nfev' % ( x, fval, nfev )
    
    status, msg = _get_saofit_msg( maxfev, ierr )
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev})
        
    return rv

def grid_search( fcn, x0, xmin, xmax, num=16, sequence=None, numcores=1,
                 maxfev=None, ftol=EPSILON, method=None, verbose=0 ):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    npar = len( x )

    def func( pars ):
        aaa = fcn( pars )[ 0 ]
        if verbose:
            print 'f%s=%g' % ( pars, aaa )
        return aaa

    def make_sequence( ranges, N ):
        list_ranges = list( ranges )
        for ii in range( npar ):
            list_ranges[ ii ] = tuple( list_ranges[ ii ] ) + ( complex( N ), )
            list_ranges[ ii ] = slice( *list_ranges[ ii ] )
        grid = numpy.mgrid[ list_ranges ]
        mynfev = pow( N, npar )
        grid = map( numpy.ravel, grid )
        sequence = []
        for index in xrange( mynfev ):
            tmp = []
            for xx in xrange( npar ):
                tmp.append( grid[ xx ][ index ] )
            sequence.append( tmp )
        return sequence

    def eval_stat_func( xxx ):
        return numpy.append( func( xxx ), xxx )

    if sequence is None:
        ranges = []
        for index in xrange( npar ):
            ranges.append( [ xmin[ index ], xmax[ index ] ] )
        sequence = make_sequence( ranges, num )
    else:
        if not numpy.iterable( sequence ):
            raise TypeError( "sequence option must be iterable" )
        else:
            for seq in sequence:
                if npar != len( seq ):
                    msg = "%s must be of length %d" % ( seq, npar )
                    raise TypeError( msg )
                    

    answer = eval_stat_func( x )
    sequence_results = parallel_map( eval_stat_func, sequence, numcores )
    for xresult in sequence_results[ 1: ]:
        if xresult[ 0 ] < answer[ 0 ]:
            answer = xresult

    fval = answer[ 0 ]
    x = answer[ 1: ]
    nfev = len( sequence_results ) + 1
    ierr = 0
    status, msg = _get_saofit_msg( ierr, ierr )
    rv = ( status, x, fval )
    rv += (msg, {'info': ierr, 'nfev': nfev })
    
    if ( 'NelderMead' == method or 'neldermead' == method or \
         'Neldermead' == method or 'nelderMead' == method ):
        #re.search( '^[Nn]elder[Mm]ead', method ):
        nm_result = neldermead( fcn, x, xmin, xmax, ftol=ftol, maxfev=maxfev,
                                verbose=verbose )
        tmp_nm_result = list(nm_result)
        tmp_nm_result_4 = tmp_nm_result[4]
        tmp_nm_result_4['nfev'] += nfev
        rv = tuple( tmp_nm_result )

    if ( 'LevMar' == method or 'levmar' == method or \
         'Levmar' == method or 'levMar' == method ):
        #re.search( '^[Ll]ev[Mm]ar', method ):
        levmar_result = lmdif( fcn, x, xmin, xmax, ftol=ftol, xtol=ftol,
                               gtol=ftol, maxfev=maxfev, verbose=verbose )
        tmp_levmar_result = list(levmar_result)
        tmp_levmar_result_4 = tmp_levmar_result[4]
        tmp_levmar_result_4['nfev'] += nfev
        rv = tuple( tmp_levmar_result )
    

    return rv

#
# Levenberg-Marquardt
#
def lmdif(fcn, x0, xmin, xmax, ftol=EPSILON, xtol=EPSILON, gtol=EPSILON,
          maxfev=None, epsfcn=EPSILON, factor=100.0, verbose=0):

    def par_at_boundary( low, val, high, tol ):
        for par_min, par_val, par_max in izip( low, val, high ):
            if sao_fcmp( par_val, par_min, tol ) == 0:
                return True
            if sao_fcmp( par_val, par_max, tol ) == 0:
                return True
        return False

    x, xmin, xmax = _check_args(x0, xmin, xmax)
    
    if maxfev is None:
       maxfev = 256 * len(x)

    def stat_cb0( pars ):
        return fcn( pars )[ 0 ]
    def stat_cb1( pars ):
        return fcn( pars )[ 1 ]

    m = numpy.asanyarray(stat_cb1(x)).size
    
    orig_fcn = stat_cb1
    error = []

    def stat_cb1(x_new, iflag):
        fvec = None
        try:
##             if _outside_limits(x_new, xmin, xmax) or _my_is_nan(x_new):
##                 fvec = numpy.empty((m,), x_new.dtype, numpy.isfortran(x_new))
##                 fvec[:] = numpy.sqrt(FUNC_MAX / m)
##             else:
            fvec = orig_fcn(x_new)
        except:
            error.append(sys.exc_info()[1])
            fvec = numpy.zeros((m,), x_new.dtype, numpy.isfortran(x_new))
            iflag = -1

        return fvec, iflag

    info, nfev, fval, covarerr = _minpack.mylmdif(stat_cb1, m, x, ftol, xtol, gtol, maxfev, epsfcn, factor, verbose, xmin, xmax)

    if par_at_boundary( xmin, x, xmax, xtol ):
        nm_result = neldermead( fcn, x, xmin, xmax, ftol=numpy.sqrt(ftol), maxfev=maxfev-nfev, finalsimplex=2, iquad=0, verbose=0 )
        nfev += nm_result[ 4 ][ 'nfev' ]
        x = nm_result[ 1 ]
        fval = nm_result[ 2 ]
##        if nm_result[ 2 ] < fval:
##            x = nm_result[ 1 ]
##            fval = nm_result[ 2 ]
##            info, mynfev, fval, covarerr = _minpack.mylmdif(stat_cb1, m, x, ftol, xtol, gtol, maxfev-nfev, epsfcn, factor, verbose, xmin, xmax)
 ##           nfev += mynfev
        
    if error:
        raise error.pop()

    key = {
        0: (False, 'improper input parameters'),
        1: (True,
            ('both actual and predicted relative reductions in the sum ' +
             'of squares are at most ftol=%g') % ftol),
        2: (True,
            ('relative error between two consecutive iterates is at ' +
             'most xtol=%g') % xtol),
        4: (True,
            ('the cosine of the angle between fvec and any column of ' +
             'the jacobian is at most gtol=%g in absolute value') % gtol),
        5: (False,
            ('number of calls to function has reached or exceeded '+
             'maxfev=%d') % maxfev),
        6: (False,
            ('ftol=%g is too small; no further reduction in the sum of ' +
             'squares is possible') % ftol),
        7: (False,
            ('xtol=%g is too small; no further improvement in the ' +
             'approximate solution is possible') % xtol),
        8: (False,
            ('gtol=%g is too small; fvec is orthogonal to the columns ' +
             'of the jacobian to machine precision') % gtol),
        }
    key[3] = (True, key[1][1] + ' and ' + key[2][1])
    status, msg = key.get(info, (False, 'unknown status flag (%d)' % info))

    if 0 == info:
        info = 1
    elif info >= 1 or info <= 4:
        info = 0
    else:
        info = 3
    status, msg = _get_saofit_msg( maxfev, info )
      
    rv = (status, x, fval)
    print_covar_err = False
    if print_covar_err:
        rv += (msg, {'info': info, 'nfev': nfev, 'covarerr': covarerr})
    else:
        rv += (msg, {'info': info, 'nfev': nfev})

    return rv


#
# Nelder Mead 
#
def minim(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, step=None,
          nloop=1, iquad=1, simp=None, verbose=-1):

    def stat_cb0( pars ):
        return fcn( pars )[ 0 ]

    x, xmin, xmax = _check_args(x0, xmin, xmax)
    
    if step is None:
        step = 0.4*numpy.ones(x.shape, numpy.float_, numpy.isfortran(x))
    if simp is None:
        simp = 1.0e-2 * ftol
    if maxfev is None:
        maxfev = 512 * len(x)

    orig_fcn = stat_cb0
    def stat_cb0(x_new):
        if _outside_limits(x_new, xmin, xmax) or _my_is_nan(x_new):        
            return FUNC_MAX
        return orig_fcn(x_new)
    
    fval, var, ifault, neval = _minim.minim(x, step, maxfev, verbose, ftol,
                                            nloop, iquad, simp, stat_cb0,
                                            xmin, xmax)
    key = {
        0: (True, 'successful termination'),
        1: (False,
            'number of function evaluations has exceeded maxfev=%d' % maxfev),
        2: (False, 'information matrix is not +ve semi-definite'),
        3: (False, 'number of parameters is less than 1'),
        4: (False, 'nloop=%d is less than 1' % nloop)
        }
    status, msg = key.get(ifault, (False, 'unknown status flag (%d)' % ifault))

    rv = (status, x, fval)
    print_covar_err = False
    if print_covar_err and 0 == ifault and numpy.any( var > 0.0 ):
        var = map( lambda x: x * numpy.sqrt(2.0), var )
        rv += (msg, {'covarerr': numpy.sqrt(var), 'info': ifault,
                     'nfev': neval})
    else:
        rv += (msg, {'info': ifault, 'nfev': neval})
    return rv


#
# Monte Carlo
#
def montecarlo(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0,
               seed=74815, population_size=None, xprob=0.9, 
               weighting_factor=0.8):

    def stat_cb0( pars ):
        return fcn( pars )[ 0 ]

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = max( 0.1, xprob )
    xprob = min( xprob, 1.0 )

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = max( 0.1, weighting_factor )
    weighting_factor = min( weighting_factor, 1.0 )

    random.seed( seed  )
    if seed is None:
        seed = random.randint(0, 2147483648L) # pow(2,31) == 2147483648L
    if population_size is None:
        population_size = max( population_size, 12 * x.size )

    if maxfev is None:
        maxfev = 8192 * population_size

    def myopt( myfcn, xxx, ftol, maxfev, seed, pop, xprob,
               weight, factor=4.0, debug=False ):

        x = xxx[ 0 ]
        xmin = xxx[ 1 ]
        xmax = xxx[ 2 ]
        maxfev_per_iter = 512 * x.size

        def random_start( xmin, xmax ):
            xx = []
            for ii in range( len(xmin ) ):
                xx.append( random.uniform( xmin[ ii ], xmax[ ii ] ) )
            return numpy.asarray( xx )

        ############################# NelderMead #############################
        mymaxfev = min( maxfev_per_iter, maxfev )
        if all( x == 0.0 ):
            mystep = map( lambda fubar: 1.2 + fubar, x )
        else: 
            mystep = map( lambda fubar: 1.2 * fubar, x )
        result = neldermead( myfcn, x, xmin, xmax, maxfev=mymaxfev, ftol=ftol,
                             finalsimplex=9, step=mystep )
        x = numpy.asarray( result[ 1 ], numpy.float_ )
        nfval = result[2]
        nfev = result[4].get( 'nfev' )
        if verbose or False != debug:
            print 'f_nm%s=%.14e in %d nfev' % ( x, nfval, nfev )
        ############################# NelderMead #############################

        ############################## nmDifEvo ##############################
        xmin, xmax = _narrow_limits( 4 * factor, [x,xmin,xmax], debug=False )
        mymaxfev = min( maxfev_per_iter, maxfev - nfev )
        result = difevo_nm( myfcn, x, xmin, xmax, ftol=ftol, maxfev=mymaxfev,
                            seed=seed, population_size=pop, xprob=xprob,
                            weighting_factor=weight )
        nfev += result[4].get( 'nfev' )
        x = numpy.asarray( result[1], numpy.float_ )
        nfval = result[2]
        if verbose or False != debug:
            print 'f_de_nm%s=%.14e in %d nfev' % ( x, result[2],
                                                   result[4].get('nfev'))
        ############################## nmDifEvo ##############################
            
        ofval = FUNC_MAX        
        while nfev < maxfev:

            xmin, xmax = _narrow_limits( factor, [x,xmin,xmax], debug=False )

            ############################ nmDifEvo #############################
            y = random_start( xmin, xmax )
            mymaxfev = min( maxfev_per_iter, maxfev - nfev )
            result = difevo_nm( myfcn, y, xmin, xmax, ftol=ftol,
                                maxfev=mymaxfev, seed=seed,
                                population_size=pop, xprob=xprob,
                                weighting_factor=weight )
            nfev += result[4].get( 'nfev' )
            if result[2] < nfval:
                nfval = result[2]
                x = numpy.asarray( result[1], numpy.float_ )
            if verbose or False != debug:
                print 'f_de_nm%s=%.14e in %d nfev' % \
                      ( x, result[2], result[4].get('nfev'))
            ############################ nmDifEvo #############################


            if False != debug:
                print 'ofval=%.14e\tnfval=%.14e\n' % (ofval, nfval)

            if sao_fcmp( ofval, nfval, ftol ) <= 0:
                return x, nfval, nfev
            ofval = nfval
            factor *= 2
            
        return x, nfval, nfev

    x, fval, nfev = myopt( fcn, [x, xmin, xmax], numpy.sqrt(ftol), maxfev,
                           seed, population_size, xprob, weighting_factor,
                           factor=2.0, debug=False )

    covarerr = None
    if nfev < maxfev:
        if all( x == 0.0 ):
            mystep = map( lambda fubar: 1.2 + fubar, x )
        else: 
            mystep = map( lambda fubar: 1.2 * fubar, x )
        result = neldermead( fcn, x, xmin, xmax,
                             maxfev=min( 512*len(x), maxfev - nfev ),
                             ftol=ftol, finalsimplex=9, step=mystep )
        covarerr = result[4].get( 'covarerr' )

        x = numpy.asarray( result[ 1 ], numpy.float_ )
        fval = result[2]
        nfev += result[4].get( 'nfev' )

    ierr = 0
    if nfev >= maxfev:
        ierr = 3
    status, msg = _get_saofit_msg( maxfev, ierr )

    rv = (status, x, fval)
    print_covar_err = False
    if print_covar_err and None != covarerr:
        rv += (msg, {'covarerr': covarerr, 'info': status, 'nfev': nfev})
    else:
        rv += (msg, {'info': status, 'nfev': nfev})
    return rv


#
# Nelder Mead 
#
def neldermead( fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None,
                initsimplex=0, finalsimplex=9, step=None, iquad=1,
                verbose=0 ):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    if step is None or ( numpy.iterable(step) and len(step) != len(x) ):
        step = 1.2*numpy.ones(x.shape, numpy.float_, numpy.isfortran(x))
    elif numpy.isscalar(step):
        step = step*numpy.ones(x.shape, numpy.float_, numpy.isfortran(x))

    def stat_cb0( pars ):
        return fcn( pars )[ 0 ]

    #
    # A safeguard just in case the initial simplex is outside the bounds
    #
    orig_fcn = stat_cb0
    def stat_cb0(x_new):
        if _outside_limits(x_new, xmin, xmax) or _my_is_nan(x_new):
            return FUNC_MAX
        return orig_fcn(x_new)

    debug = False             # for internal use only
    
    if numpy.isscalar(finalsimplex) and 0 == numpy.iterable(finalsimplex):
        finalsimplex = int( finalsimplex )
        if 0 == finalsimplex:
            finalsimplex = [ 1 ]
        elif 1 == finalsimplex:
            finalsimplex = [ 2 ]
        elif 2 == finalsimplex:
            finalsimplex = [ 0, 0 ]
        elif 3 == finalsimplex:            
            finalsimplex = [ 0, 1 ]
        elif 4 == finalsimplex:
            finalsimplex = [ 0, 1, 0 ]
        elif 5 == finalsimplex:
            finalsimplex = [ 0, 2, 0 ]
        elif 6 == finalsimplex:
            finalsimplex = [ 1, 1, 0 ]
        elif 7 == finalsimplex:
            finalsimplex = [ 2, 1, 0 ]
        elif 8 == finalsimplex:
            finalsimplex = [ 1, 2, 0 ]
        elif 9 == finalsimplex:
            finalsimplex = [ 0, 1, 1 ]            
        elif 10 == finalsimplex:
            finalsimplex = [ 0, 2, 1 ]            
        elif 11 == finalsimplex:
            finalsimplex = [ 1, 1, 1 ]                        
        elif 12 == finalsimplex:
            finalsimplex = [ 1, 2, 1 ]            
        elif 13 == finalsimplex:
            finalsimplex = [ 2, 1, 1 ]            
        else:
            finalsimplex = [ 2, 2, 2 ]
    elif ( False == numpy.isscalar(finalsimplex) and
           1 == numpy.iterable(finalsimplex) ):
        pass
    else:
        finalsimplex = [ 2, 2, 2 ]
        
    finalsimplex = numpy.asarray(finalsimplex, numpy.int_)

    if maxfev is None:
        maxfev = 1024 * len( x )

    if debug:
        print 'opfcts.py neldermead() finalsimplex=%s\tisscalar=%s\titerable=%d' % (finalsimplex,numpy.isscalar(finalsimplex), numpy.iterable(finalsimplex))

    def simplex( verbose, maxfev, init, final, tol, step, xmin, xmax, x,
                 myfcn, debug, ofval=FUNC_MAX ):

        tmpfinal = final[:]
        if len( final ) >= 3:
            tmpfinal = final[0:-1] # get rid of the last entry in the list

        xx,ff,nf,er = _saoopt.neldermead( verbose, maxfev, init, tmpfinal, tol,
                                          step, xmin, xmax, x, myfcn )

        if debug:
            print 'finalsimplex=%s, nfev=%d:\tf%s=%.20e' % (tmpfinal,nf,xx,ff)

        if len( final ) >= 3 and ff < 0.995 * ofval and nf < maxfev:
            myfinal = [final[-1]]
            x,fval,nfev,err = simplex( verbose, maxfev-nf, init, myfinal, tol,
                                       step, xmin, xmax, x, myfcn, debug,
                                       ofval=ff )
            return x,fval,nfev+nf,err
        else:
            return xx,ff,nf,er


    x, fval, nfev, ier = simplex( verbose, maxfev, initsimplex, finalsimplex,
                                  ftol, step, xmin, xmax, x, stat_cb0, debug )
    if debug:
        print 'f%s=%f in %d nfev' % ( x, fval, nfev )
    
    info=1
    covarerr=None
    if len( finalsimplex ) >= 3 and 0 != iquad:
        nelmea = minim( fcn, x, xmin, xmax, ftol=10.0*ftol, maxfev=maxfev-nfev-12, iquad=1 )
        nelmea_x = numpy.asarray( nelmea[1], numpy.float_ )
        nelmea_nfev = nelmea[4].get( 'nfev' )
        info = nelmea[4].get( 'info' )
        covarerr = nelmea[4].get( 'covarerr' )
        nfev += nelmea_nfev
        minim_fval = nelmea[ 2 ]
        if minim_fval < fval:
            x = nelmea_x
            fval = minim_fval
        if debug:
            print 'minim: f%s=%e %d nfev, info=%d' % (x,fval,nelmea_nfev,info)
            
    if nfev >= maxfev:
        ier = 3
    key = {
        0: (True, 'Optimization terminated successfully'),
        1: (False, 'improper input parameters'),
        2: (False, 'improper values for x, xmin or xmax'),
        3: (False,
            'number of function evaluations has exceeded %d' % maxfev)
        }
    status, msg = key.get( ier,
                           (False, 'unknown status flag (%d)' % ier) )

    rv = (status, x, fval)
    print_covar_err = False
    if print_covar_err and None != covarerr:
        rv += (msg, {'covarerr': covarerr, 'info': status, 'nfev': nfev})
    else:
        rv += (msg, {'info': status, 'nfev': nfev})
    return rv


def lmdif_cpp(fcn, x0, xmin, xmax, ftol=EPSILON, xtol=EPSILON, gtol=EPSILON,
              maxfev=None, epsfcn=EPSILON, factor=100.0, verbose=0):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    if maxfev is None:
        maxfev = 256 * len(x)

    def stat_cb0( pars ):
        return fcn( pars )[ 0 ]
    def stat_cb1( pars ):
        return fcn( pars )[ 1 ]

    m = numpy.asanyarray(stat_cb1(x)).size

    orig_fcn = stat_cb1
    error = []

    def stat_cb1(x_new):
        return orig_fcn(x_new)

    x, fval, nfev, info, covarerr = _saoopt.cpp_lmdif( stat_cb1, m, x, ftol, xtol, gtol, maxfev, epsfcn, factor, verbose, xmin, xmax)
    
    if error:
        raise error.pop()

    key = {
        0: (False, 'improper input parameters'),
        1: (True,
            ('both actual and predicted relative reductions in the sum ' +
             'of squares are at most ftol=%g') % ftol),
        2: (True,
            ('relative error between two consecutive iterates is at ' +
             'most xtol=%g') % xtol),
        4: (True,
            ('the cosine of the angle between fvec and any column of ' +
             'the jacobian is at most gtol=%g in absolute value') % gtol),
        5: (False,
            ('number of calls to function has reached or exceeded '+
             'maxfev=%d') % maxfev),
        6: (False,
            ('ftol=%g is too small; no further reduction in the sum of ' +
             'squares is possible') % ftol),
        7: (False,
            ('xtol=%g is too small; no further improvement in the ' +
             'approximate solution is possible') % xtol),
        8: (False,
            ('gtol=%g is too small; fvec is orthogonal to the columns ' +
             'of the jacobian to machine precision') % gtol),
        }
    key[3] = (True, key[1][1] + ' and ' + key[2][1])
    status, msg = key.get(info, (False, 'unknown status flag (%d)' % info))

    rv = (status, x, fval)
    rv += (msg, {'info': info, 'nfev': nfev, 'covarerr': covarerr })
    return rv
