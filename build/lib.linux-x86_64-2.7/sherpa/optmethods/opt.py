from copy import copy
import numpy
import time

import multiprocessing
from sherpa.utils._utils import sao_fcmp
from sherpa.utils import Knuth_close, func_counter, is_iterable, is_in, \
     divide_run_parallel
from sherpa.optmethods.optfcts import _check_args, EPSILON, _get_saofit_msg, \
     _outside_limits, _narrow_limit, neldermead
from sherpa.optmethods.optfcts import neldermead as sherpa_neldermead

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

def func_bounds( func, low, high ):
    def func_bounds_wrapper( x, *args ):
        if _outside_limits( x, low, high ):
            return numpy.inf
        return func( x, *args )
    return func_bounds_wrapper

class Opt( object ):
    
    def __init__( self, fcn ):
        self.nfev, self.fcn = func_counter( fcn )
        
    def check_args( self, x, xmin, xmax ):
        x = numpy.array( x, numpy.float_ )          # Make a copy
        xmin = numpy.asarray( xmin, numpy.float_ )  # Make a copy
        xmax = numpy.asarray( xmax, numpy.float_  ) # Make a copy
        if ( x.shape != xmin.shape ) or ( x.shape != xmax.shape ):
            raise InputError( x, xmin, xmax )
        if _outside_limits( x, xmin, xmax ):
            raise OutOfBoundErr( x, xmin, xmax )
        return x, xmin, xmax
    
class Polytope( Opt ):
    
    def __init__( self, fcn, x, xmin, xmax, step, initsimplex ):
        Opt.__init__( self, fcn )
        x, self.xmin, self.xmax = self.check_args( x, xmin, xmax )
        self.fcn = func_bounds( self.fcn, xmin, xmax )
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
    
    def contract_in_out( self, centroid, reflection_pt, rho_gamma,
                         contraction_coef, badindex, verbose ):

        if self.polytope[ badindex - 1, -1 ] <= reflection_pt[ -1 ] and \
               reflection_pt[ -1 ] < self.polytope[ badindex, -1 ]:
            
            # Perform outside contraction
            outside_contraction_pt = self.move_vertex( centroid, rho_gamma,
                                                       badindex )

            if outside_contraction_pt[ -1 ] <= reflection_pt[ -1 ]:
                #self.polytope[ badindex ] = outside_contraction_pt[ : ]
                self.replace_vertex( badindex, outside_contraction_pt )
                if verbose:
                    print '\taccept outside contraction point'
                return False
            else:
                return True
            
        elif reflection_pt[ -1 ] >= self.polytope[ badindex, -1 ]:

            # Perform an inside contraction
            inside_contraction_pt = self.move_vertex( centroid,
                                                      - contraction_coef,
                                                      badindex )
                    
            if inside_contraction_pt[ -1 ] < self.polytope[ badindex, -1 ]:
                #self.polytope[ badindex ] = inside_contraction_pt[ : ]
                self.replace_vertex( badindex, inside_contraction_pt )
                if verbose:
                    print '\taccept inside contraction point'
                return False
            else:
                return True

        else:
            print 'something is wrong with contract_in_out'
            return True
        
    def get_func_vals( self ):
        return self.polytope[ :, -1 ]

    def order_simplex( self ):
        func_vals = self.get_func_vals( )
        # return a list of indices goes from largest to smallest
        index_from_largest_to_smallest = numpy.argsort( func_vals )
        # re-arrange so self.polytope goes from smallest to largest fct val
        self.polytope = numpy.take( self.polytope,
                                   index_from_largest_to_smallest, 0 )

    def init_simplex( self, x, step, initsimplex ):
        npar = len(x)
        npar_plus_1 = npar + 1
        simplex = numpy.zeros( (npar_plus_1,npar_plus_1), dtype=numpy.float_ )
        if step is None:
            step = 0.4*numpy.ones(x.shape, numpy.float_, numpy.isfortran(x))
        for ii in xrange( npar_plus_1 ):
            simplex[ ii , :-1 ] = numpy.array( x, copy=True )
        if 0 == initsimplex:
            for ii in xrange( npar ):
                simplex[ ii + 1 ,  ii ] += step[ ii ]
        else:
            delta = 0.05
            for ii in xrange( nparx ):
                if 0.0 == simplex[ ii + 1 ,  ii ]:
                    simplex[ ii + 1 ,  ii ] = delta
                else:
                    simplex[ ii + 1 ,  ii ] *= (1.0+delta)
        for ii in xrange( npar + 1 ):
            simplex[ ii , -1 ] = self.fcn( simplex[ ii , :-1 ] )
        return simplex

    def move_vertex( self, centroid, coef, badindex ):
        vertex = ( 1.0 + coef ) * centroid - coef * self.polytope[ badindex ]
        vertex[ -1 ] = self.fcn( vertex[ :-1 ] )
        return vertex

    def replace_vertex( self, badindex, newvertex ):
        if newvertex[ -1 ] < 0.95 * self.polytope[0,-1]:
            print 'polytope[0]=', self.polytope[0,-1], '\tpolytope[', badindex, ']=', self.polytope[badindex,-1]
            #result = dmnfb( self.fcn, newvertex[ :-1 ], self.xmin, self.xmax )
            #result = newuoa( self.fcn, newvertex[ :-1 ], self.xmin, self.xmax )
            newvertex[ :-1 ] = numpy.asarray( result[ 1 ], numpy.float_ )
            newvertex[ -1 ] = result[ 2 ]
            print 'f(%s)=%f' % (newvertex[:-1],newvertex[-1])
        self.polytope[ badindex ] = newvertex[ : ]
        
    def shrink( self, shrink_coef, npar ):
        npars_plus_1 = npar + 1
        for ii in xrange( 1, npars_plus_1 ):
            self.polytope[ ii ] = self.polytope[ 0 ] + shrink_coef * ( self.polytope[ ii ] - self.polytope[ 0 ] )
            self.polytope[ ii, -1 ] = self.fcn( self.polytope[ ii, :-1 ] )

class classicSimplex( Polytope ):

    def __init__( self, fcn, npar, aflatarray, xmin, xmax ):
        self.fcn = fcn
        self.polytope = aflatarray
        self.polytope.shape = ( npar + 1, npar + 1 )
        self.xmin = xmin
        self.xmax = xmax

class DirectSearch( Opt ):

    def __init__( self, fcn ):
        Opt.__init__( self, fcn )
        self.expansion_coef   = 2.0          # chi
        self.contraction_coef = 0.5          # gamma
        self.reflection_coef  = 1.0          # rho
        self.shrink_coef      = 0.5          # sigma

class classicNelderMead( DirectSearch ):

    def __init__( self, fcn ):
        DirectSearch.__init__( self, fcn )
        
    def __call__( self, x, xmin, xmax, maxfev, tol, step, initsimplex,
                  finalsimplex, multicore, verbose ):

        try :

            def optimize( bad_index, fcn, npar, aflatarray ):

                rho_chi = self.reflection_coef * self.expansion_coef
                rho_gamma = self.reflection_coef * self.contraction_coef

                result = []
                for badindex in bad_index:

                    begin_nfev = self.nfev[0]

                    simplex = classicSimplex( fcn, npar, aflatarray,
                                              xmin, xmax )

                    shrinkme = False

                    centroid = simplex.calc_centroid( badindex )

                    reflection_pt = simplex.move_vertex( centroid,
                                                         self.reflection_coef,
                                                         badindex )

                    if simplex[ 0, -1 ] <= reflection_pt[ -1 ] and \
                           reflection_pt[ -1 ] < simplex[ badindex - 1, -1 ]:
                        #simplex[ badindex ] = reflection_pt[ : ]
                        simplex.replace_vertex( badindex, reflection_pt )
                        if verbose:
                            print '\taccept reflection point'

                    elif reflection_pt[ -1 ] < simplex[ 0, -1 ]:

                        # calculate the expansion point
                        expansion_pt = simplex.move_vertex( centroid, rho_chi,
                                                            badindex )

                        if expansion_pt[ -1 ] < reflection_pt[ -1 ]:
                            #simplex[ badindex ] = expansion_pt[ : ]
                            simplex.replace_vertex( badindex, expansion_pt )
                            if verbose:
                                print '\taccept expansion point'
                        else:
                            #simplex[ badindex ] = reflection_pt[ : ]
                            simplex.replace_vertex( badindex, reflection_pt )
                            if verbose:
                                print '\taccept reflection point'

                    else: 

                        shrinkme = simplex.contract_in_out( centroid,
                                                            reflection_pt,
                                                            rho_gamma,
                                                            self.contraction_coef,
                                                            badindex, verbose )
                    if shrinkme is False:
                        shrinkme = 0
                    else:
                        skrinkme = 1
                    result.append( numpy.append( simplex.polytope.ravel( ),
                                                 [self.nfev[0]-begin_nfev,
                                                  shrinkme, badindex] ) )
                        
                return numpy.asarray( result )

            x, xmin, xmax = self.check_args( x, xmin, xmax )

            polytope = Polytope( self.fcn, x, xmin, xmax, step, initsimplex )

            npar = len( x )
            bad_vertices = [ npar ]
            mynfev = 0

            while ( self.nfev[ 0 ] < maxfev ):

                polytope.order_simplex( )

                if verbose:
                    print 'f%s=%f' % (polytope[ 0, :-1 ],polytope[ 0, -1 ])

                if polytope.check_convergence( tol, finalsimplex ):
                    break

                if multicore:
                    badindex = multiprocessing.cpu_count()
                    if badindex == npar:
                        badindex = 1
                    bad_vertices = range( npar, npar - badindex, -1 )

                    results = divide_run_parallel( optimize, bad_vertices,
                                                   polytope.fcn, npar,
                                                   polytope.polytope.ravel( ) )
                else:
                    results = optimize( bad_vertices, polytope.fcn, npar,
                                        polytope.polytope.ravel( ) )

                datasize = npar * ( npar + 2 ) + 4
                nrow = results.size / datasize
                ncol = results.size / nrow
                results = results.reshape( nrow, ncol )

                num_not_modified = 0
                for result in results:
                    mynfev += int( result[-3] )
                    shrinkme = int( result[-2] )
                    index = int( result[-1] )
                    aflatarray = result[:-3]
                    if shrinkme:
                        num_not_modified += 1
                    else:
                        aflatarray.shape = ( npar + 1, npar + 1 )
                        #print index, polytope[ index, -1 ], aflatarray[ index, -1 ]
                        polytope[ index ] = numpy.array( aflatarray[ index ],
                                                         copy=True )

                if num_not_modified == len( bad_vertices ):
                    polytope.shrink( self.shrink_coef, npar )
                    if verbose:
                        print '\tshrink'

            ierr = 0
            if mynfev >= maxfev:
                ierr = 3

            vertex = polytope[ 0 ]
            bestpar = vertex[ :-1 ]
            bestval = vertex[ -1 ]
            return ierr, bestpar, bestval, mynfev
        except InputErr:
            vertex = polytope.get_vertex( 0 )
            bestpar = vertex[ :-1 ]
            bestval = vertex[ -1 ]
            return 1, bestpar, bestval, mynfev
        except OutOfBoundErr:
            return 2, x, numpy.nan, mynfev


def get_result( arg, maxfev ):
    ierr = arg[ 0 ]
    status, msg = _get_saofit_msg( maxfev, ierr )
    xopt = arg[ 1 ]
    fopt = arg[ 2 ]
    nfev = arg[ 3 ]
    rv = (status, xopt, fopt)
    rv += (msg, {'info': ierr, 'nfev': nfev})
    return rv


def optneldermead( afcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None,
                   initsimplex=0, finalsimplex=None, step=None,
                   multicore=False, verbose=None ):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    nfev, fcn = func_counter( afcn )
    
    if step is None or ( numpy.iterable(step) and len(step) != len(x) ):
        step = 1.2*numpy.ones(x.shape, numpy.float_, numpy.isfortran(x))
    elif numpy.isscalar(step):
        step = step*numpy.ones(x.shape, numpy.float_, numpy.isfortran(x))

    if maxfev is None:
        maxfev = 1024 * len( x )

    if finalsimplex is None:
        finalsimplex = [ 0, 2 ]

    if False == is_iterable( finalsimplex ):
        finalsimplex = [ finalsimplex ]
    
    # For internal use only:
    debug = True

    def myloop( xxx, mymaxfev, mystep ):

        nms = classicNelderMead( fcn )
        mynfev = 0
        for myfinalsimplex in finalsimplex:
            starttime = time.time( )            
            result = nms( xxx, xmin, xmax, mymaxfev - mynfev, ftol, mystep,
                          initsimplex, myfinalsimplex, multicore, verbose )
            mynfev += result[ 3 ]            
            if 0 != result[ 0 ]:
                return result[0], result[1], result[2], mynfev
            if debug:
                print 'neldermead::myloop: result = ', result, ' took ', time.time() - starttime, ' secs'
        if multicore:
            return result[0], result[1], result[2], mynfev + nfev[0]
        else:
            return result[0], result[1], result[2], mynfev

    def myneldermead( xxx, mymaxfev, mystep,
                      fold=numpy.float_(numpy.finfo(numpy.float_).max) ):
        result = myloop( xxx, mymaxfev, mystep )
        mynfev = result[ 3 ]
        fnew = result[ 2 ]
        if fnew < fold * 0.995 and 0 == Knuth_close( fnew, fold, ftol ) and \
               0 == Knuth_close( fnew, fold, ftol ) and \
               0 == result[0] and mynfev < mymaxfev:
            if debug:
                print 'neldermead: fnew = ', fnew, '\tfold = ', fold, \
                      '\tx = ', result[ 1 ]
            result = myneldermead( result[ 1 ], mymaxfev - mynfev, step,
                                   fold=fnew )
            mynfev += result[ 3 ]
        return result[0], result[1], result[2], mynfev

    #result = myneldermead( x, maxfev, step )
    result = myloop( x, maxfev, step )
    return get_result( result, maxfev )

import random
class DifEvo( Polytope ):

    def __init__( self, fcn ):
        Opt.__init__( self, fcn )
        self.trial_solution = None

    def check_convergence( self, polytope, tolerance ):

        def is_max_length_small_enough( polytope, tol ):
            """
                       
               max  || x  - x  || <= tol max( 1.0, || x || )
                        i    0                         0
                        
            where 1 <= i <= n"""

            max_xi_x0 = -1.0
            x0 = polytope[ 0, :-1 ]

            npar_plus_1 = len( x0 ) + 1
            for ii in xrange( 1, npar_plus_1 ):
                xi = polytope[ ii, :-1 ]
                xi_x0 = xi - x0
                max_xi_x0 = max( max_xi_x0, numpy.dot( xi_x0, xi_x0 ) )
            return max_xi_x0 <= tol * max( 1.0, numpy.dot( x0, x0 ) )

        if is_max_length_small_enough( polytope, tolerance ):
            return True
        func_vals = polytope[ :, -1 ]
        return numpy.std( func_vals ) < tolerance


    def init_population( self, x, xmin, xmax, population_size, seed ):
        npar = len(x)
        population = numpy.zeros( (population_size, npar + 1),
                                  dtype=numpy.float_ )
        random.seed( seed )
        for ii in xrange( population_size - 1 ):
            for jj in xrange( npar ):
                population[ ii, jj ] = random.uniform( xmin[ jj ],
                                                         xmax[ jj ] )
        population[ population_size - 1, :-1 ] = x[ : ]
            
        for ii in xrange( population_size ):
            population[ ii, -1 ] = self.fcn( population[ ii, :-1 ] )
            #print 'f%s=%f' % (population[ ii, :-1 ], population[ ii, -1 ])
        return population, population[ population_size - 1, -1 ]

    def select_samples( self, candidate, population_size, n ):

        mysample = random.sample( xrange( population_size ), n )
        #print 'select_samples: mysample = ', mysample
        if is_in( candidate, mysample ):
            return self.select_samples( candidate, population_size, n )
        else:
            return mysample

    def best1exp( self, candidate, population, model_par, xprob, scale_factor ):

        population_size = population.shape[ 0 ]
        npar = population.shape[ 1 ] - 1
        r = self.select_samples( candidate, population_size, 2 )
        n = random.randint( 0, npar - 1 )
            
        for ii in xrange( npar ):
            if random.random() >= xprob:
                return
            self.trial_solution[ n ] = model_par[ n ] + \
                                       scale_factor * ( population[ r[0], n ] - \
                                                        population[ r[1], n ] )
            n = ( n + 1 ) % npar
                                       
    def __call__( self, x, xmin, xmax, maxnfev, tolerance, population_size,
                  xprob, scale_factor, seed=81547, multicore=False,
                  verbose=None ):

        debug = False
        use_local_opt = True
        
        x, self.xmin, self.xmax = self.check_args( x, xmin, xmax )
        self.fcn = func_bounds( self.fcn, xmin, xmax )
        polytope, fstat = self.init_population( x, xmin, xmax, population_size, seed )
        strategy_function = self.best1exp
        
        npar = len( x )
        npar_plus_1 = npar + 1

        self.trial_solution = numpy.zeros( npar_plus_1, dtype=numpy.float_ )
        fct_vals = numpy.zeros( npar, dtype=numpy.float_ )

        model_par = x[ : ]
        while self.nfev[ 0 ] < maxnfev:
            for candidate in xrange( population_size ):

                self.trial_solution[ : ] = polytope[ candidate ][ : ]
                strategy_function( candidate, polytope, model_par,
                                   xprob, scale_factor )
                self.trial_solution[ -1 ] = self.fcn( self.trial_solution[ :-1 ] )
                if self.trial_solution[ -1 ] < polytope[ candidate, -1 ]:
                    polytope[ candidate ][ : ] = self.trial_solution[ : ]

                    if self.trial_solution[ npar ] < fstat:

                        if use_local_opt:
                            if debug:
                                print 'before nm: f%s=%f' % \
                                      ( self.trial_solution[ :-1 ], \
                                        self.trial_solution[ -1 ] )

                            result = sherpa_neldermead( self.fcn,
                                                        self.trial_solution[ :-1 ],
                                                        xmin, xmax,
                                                        ftol=tolerance,
                                                        maxfev=maxnfev - self.nfev[0],
                                                        initsimplex=0,
                                                        finalsimplex=[0,1],
                                                        step=None )
                            model_par[ : ] = result[ 1 ]
                            fstat = result[ 2 ]

                            if debug:
                                print 'after nm: f%s=%f' % \
                                      ( result[ 1 ], result[ 2 ])
                                print

                        else:
                            model_par[ : ] = self.trial_solution[ :-1 ]
                            fstat = self.trial_solution[ -1 ]
                            
                        if self.check_convergence( polytope, tolerance ) :
                            return 0, model_par, fstat, self.nfev[ 0 ]

                if self.nfev[ 0 ] >= maxnfev:
                    3, model_par, fstat, self.nfev[ 0 ]
                    
        return 3, model_par, fstat, self.nfev[ 0 ]


            
def difevo( fcn, x0, xmin, xmax, maxfev=None, ftol=EPSILON,xprob=0.9, 
            weighting_factor=0.8, population_size=None, multicore=True,
            verbose=0 ):

    x0, xmin, xmax = _check_args(x0, xmin, xmax)
    
    npar = len( x0 )
    if maxfev is None:
        maxfev = 4096 * npar * 16
    if population_size is None:
        population_size = 16 * npar

    starttime = time.time()
    nm_result = sherpa_neldermead( fcn, x0, xmin, xmax, ftol=ftol,
                                   maxfev=maxfev, initsimplex=0,
                                   finalsimplex=9, step=None, verbose=0)

    print nm_result
    print '%f secs' % ( time.time( ) - starttime )

    if nm_result[0]:
        
        dif_evo = DifEvo( fcn )

        maxfev_per_iter = min( maxfev - nm_result[4].get( 'nfev' ),
                               512 * npar )
        factor = 4.0
        myx = nm_result[ 1 ]
        myxmin, myxmax = _narrow_limit( 4 * factor, [myx,xmin, xmax],
                                        debug=True )
        starttime = time.time()
        de_result = dif_evo( myx, myxmin, myxmax, maxfev_per_iter, ftol,
                             population_size, xprob, weighting_factor )
        print de_result
        print '%f secs' % ( time.time( ) - starttime )
        
        de_result = list( de_result )
        de_result[ -1 ] += nm_result[4].get( 'nfev' )
        return get_result( de_result, maxfev )
    
    else:
        
        return nm_result

def FreudensteinRoth( x ):
    npar = len( x )
    sum = 0.0
    for ii in xrange( 0, npar, 2 ):
        fx0 = - 13.0 + x[ ii ] + x[ ii + 1 ] * (( 5.0 - x[ ii + 1 ] ) * x[ ii + 1 ] - 2.0 )
        fx1 = - 29.0 + x[ ii ] + x[ ii + 1 ] * (( x[ ii + 1 ] + 1.0 ) * x[ ii + 1 ] - 14.0 )
        sum += fx0*fx0 + fx1*fx1

    return sum

def rosen(x):  # The Rosenbrock function
    #print 'hi'
    x = numpy.asarray(x)
    val = numpy.sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0,axis=0)
    #print 'f%s=%f' % ( x, val )
    return val

def tst_prob0( multicore ):
    def prob0( x ):
        return 2.0 * x[ 0 ] * x[ 0 ] + x[ 1 ] * x[ 1 ] + 3.0 * x[ 2 ] * x[ 2 ]

    x0 = ( 2., -2., 1.0 )
    xmin = ( -10., -10, -10 )
    xmax = ( 10., 10, 10 )
    maxfev = 4096
    verbose = False
    ftol = numpy.float_(numpy.finfo(numpy.float32).eps)
    step = None
    initsimplex = 0
    finalsimplex = 0
    starttime = time.time( )
    print neldermead( prob0, x0, xmin, xmax, ftol, maxfev, initsimplex,
                      finalsimplex, step=None, verbose=False,
                      multicore=multicore )
    print '%f secs' % ( time.time( ) - starttime )
    solution = ( 0.0, 3.0e-6, 1.0e-6 )
    print prob0( solution )
    
def tst_prob1( multicore ):
    def prob1( x ):
        return 5.0 * x[ 0 ] * x[ 0 ] + 2.0 * x[ 1 ] * x[ 1 ] + \
               2.0 * x[ 2 ] * x[ 2 ] + 2.0 * x[ 0 ] * x[ 1 ] + \
               2.0 * x[ 1 ] * x[ 2 ] - 2.0 * x[ 0 ] * x[ 2 ] - 6.0 * x[ 2 ]

    x0 = ( 0.0, 0.0, 0.0 )
    xmin = ( -10., -10, -10 )
    xmax = ( 10., 10, 10 )
    maxfev = 4096
    verbose = False
    ftol = numpy.float_(numpy.finfo(numpy.float32).eps)
    step = None
    initsimplex = 0
    finalsimplex = 0
    starttime = time.time( )
    print neldermead( prob1, x0, xmin, xmax, ftol, maxfev, initsimplex,
                      finalsimplex, step=None, verbose=False,
                      multicore=multicore )
    print '%f secs' % ( time.time() - starttime )
    solution = ( 1.0, -2.0, 3.0 )
    print prob1( solution )

def tst_prob2( multicore ):
    def prob2( x ):
        return 2.0 * x[ 0 ] * x[ 0 ] + x[ 1 ] * x[ 1 ] - x[ 0 ] * x[ 1 ]

    x0 = ( 2.0, 2.0 )
    xmin = ( -10., -10 )
    xmax = ( 10., 10 )
    maxfev = 4096
    verbose = False
    ftol = numpy.float_(numpy.finfo(numpy.float32).eps)
    step = None
    initsimplex = 0
    finalsimplex = 0
    starttime = time.time( )
    print neldermead( prob2, x0, xmin, xmax, ftol, maxfev, initsimplex,
                      finalsimplex, step=None, verbose=False,
                      multicore=multicore )
    print '%f secs' % ( time.time() - starttime )
    solution = ( 0.0, 0.0 )
    print prob2( solution )
    
def tst_rosen( npar, multicore ):
    x0 = npar * [ -1.2, 1.0 ]
    xmin = npar * [ -10, -10. ]
    xmax = npar * [ 10, 10 ]
    maxfev = npar * 4096 * 16
    verbose = False
    ftol = numpy.float_(numpy.finfo(numpy.float32).eps)
    step = None
    initsimplex = 0
    finalsimplex = 0
    starttime = time.time( )
    print neldermead( rosen, x0, xmin, xmax, ftol, maxfev, initsimplex,
                      finalsimplex, step=step, verbose=verbose,
                      multicore=multicore )
    print '%f secs' % ( time.time( ) - starttime )

def tst_de_rosen( npar ):
    x0 = npar * [ -1.2, 1.0 ]
    xmin = npar * [ -10, -10. ]
    xmax = npar * [ 10, 10 ]
    maxfev = npar * 4096 * 16
    verbose = False
    ftol = numpy.float_(numpy.finfo(numpy.float32).eps)
    starttime = time.time( )
    print difevo( rosen, x0, xmin, xmax, maxfev, ftol )
    print '%f secs' % ( time.time( ) - starttime )
    
def tst_de_freudensteinroth( npar ):
    x0 = npar * [ 0.5, -2.0 ]
    xmin = npar * [ -20, -20. ]
    xmax = npar * [ 20., 20. ]
    maxfev = npar * 4096 * 16
    verbose = False
    ftol = numpy.float_(numpy.finfo(numpy.float32).eps)
    starttime = time.time( )
    print difevo( FreudensteinRoth, x0, xmin, xmax, maxfev, ftol )
    print '%f secs' % ( time.time( ) - starttime )

if __name__ == '__main__':
    import sys
    npar = 2
    if ( 2 == len(sys.argv) ):
        npar = int( sys.argv[1] )
    multicore = False
##     tst_prob0( multicore )        
##     tst_prob1( multicore )
##     tst_prob2( multicore )    
##    tst_rosen( npar, multicore )
    #tst_de_rosen( npar )
    tst_de_freudensteinroth( npar )
