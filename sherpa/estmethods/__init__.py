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
_ = numpy.seterr(invalid='ignore')

from sherpa.utils import NoNewAttributesAfterInit, print_fields, Knuth_close, is_iterable, list_to_open_interval, mysgn, quad_coef, apache_muller, bisection, demuller, zeroin, OutOfBoundErr, func_counter, _multi, _ncpus

import logging
import sherpa.estmethods._est_funcs
from itertools import izip

try:
    import multiprocessing
except:
    pass


__all__ = ('EstNewMin', 'Covariance', 'Confidence',
           'Projection', 'est_success', 'est_failure', 'est_hardmin',
           'est_hardmax', 'est_hardminmax', 'est_newmin', 'est_maxiter',
           'est_hitnan')

est_success       = 0
est_failure       = 1
est_hardmin       = 2
est_hardmax       = 3
est_hardminmax    = 4
est_newmin        = 5
est_maxiter       = 6
est_hitnan        = 7

# For every method listed here, we have the same goal:  derive confidence
# limits for thawed parameters.  Thawed parameters are allowed to vary
# during a fit; when a model has been fit to data, then the current
# parameter values are presumably the best-fit values.  Think of the
# best-fit values as being at the lowest point in a valley in parameter
# space--any step away, in any direction, means a worse fit (i.e., the
# value of the fit statistic is greater).  Confidence limits tell you
# how well constrained those best fit values are; i.e., are we in a deep,
# narrow valley in parameter space?  If so, we can be confident the limits
# are small.  But if the valley is shallow and broad, then the confidence
# limits will also be very broad.
#
# Every method is passed the same information:
#  the current values of all thawed parameters;
#  the soft limits of all thawed parameters;
#  the hard limits of all thawed parameters;
#  the list of parameters for which we are actually want confidence
#    limits (this can be a subset of all thawed parameters)
#  a reference to the statistic function;
#  a reference to the fitting function.


#class EstMethodError(SherpaError):
#    "Reached an error while computing parameter confidence limits"
#    pass
#
#class EstHardMin(EstMethodError):
#    "Reached a parameter hard minimum"
#    pass
#
#class EstHardMax(EstMethodError):
#    "Reached a parameter hard maximum"
#    pass
#
class EstNewMin(Exception):
    "Reached a new minimum fit statistic"
    pass
#
#class EstMaxIter(EstMethodError):
#    "Reached maxmimum iterations in scaling function"
#    pass
#
#class EstNaN(EstMethodError):
#    "Reached a NaN during computation"
#    pass

class EstMethod(NoNewAttributesAfterInit):

    # defined pre-instantiation for pickling
    config = {'sigma' : 1,
              'eps' : 0.01,
              'maxiters' : 200,
              'soft_limits' : False}

    def __init__(self, name, estfunc):
        self._estfunc = estfunc
        self.name = name

        # config should be defined pre-instantiation for pickling
        # however, for some unknown reason membership in self.__dict__
        # requires declaration in __init__()
        self.config = self.config.copy()

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
        return ("<%s error-estimation method instance '%s'>" %
                (type(self).__name__, self.name))

    def __str__(self):
        # Put name first always 
        keylist = self.config.keys()
        keylist = ['name'] + keylist
        full_config = {'name' : self.name}
        full_config.update(self.config)
        
        return print_fields(keylist, full_config)

    def __setstate__(self, state):
        self.__dict__.update(state)

        # obtain config values from object class
        self.__dict__['config'] = getattr(self.__class__(),'config',{})

        # update new config dict with user defined from old
        self.__dict__['config'].update(state.get('config',{}))


    def compute(self, statfunc, fitfunc, pars,
                parmins, parmaxes, parhardmins,
                parhardmaxes, limit_parnums, freeze_par, thaw_par,
                report_progress, get_par_name,
                statargs=(), statkwargs={}):

        def stat_cb(pars):
            return statfunc(pars)[0]
        
        def fit_cb(scb, pars, parmins, parmaxes, i):
            # parameter i is a no-op usually
            return fitfunc(scb, pars, parmins, parmaxes)[2]

        # remin means reminimize -- *generally* not done (only
        # proj needs to reminimize).
        #
        # covar still needs to pass a reminimize flag to
        # get_one_sided_interval; a value less than zero is
        # interpreted as "never reminimize", so pass that value here.

        remin = -1.0
        tol = -1.0
        return self._estfunc(pars, parmins, parmaxes, parhardmins,
                             parhardmaxes, self.sigma, self.eps,
                             tol,
                             self.maxiters, remin, limit_parnums,
                             stat_cb, fit_cb, report_progress)


class Covariance(EstMethod):

    def __init__(self, name='covariance'):
        EstMethod.__init__(self, name, covariance)

    
class Confidence(EstMethod):

    # defined pre-instantiation for pickling
    _added_config = {'remin': 0.01,
                     'fast': False,
                     'parallel': True,
                     'numcores' : _ncpus,
                     'maxfits' : 5,
                     'max_rstat' : 3,
                     'tol' : 0.2,
                     'verbose' : False,
                     'openinterval': False }

    def __init__(self, name='confidence'):
        EstMethod.__init__(self, name, confidence)

        # Update EstMethod.config dict with Confidence specifics
        self.config.update(self._added_config)

    def compute(self, statfunc, fitfunc, pars,
                parmins, parmaxes, parhardmins,
                parhardmaxes, limit_parnums, freeze_par, thaw_par,
                report_progress, get_par_name,
                statargs=(), statkwargs={}):

        def stat_cb(pars):
            return statfunc(pars)[0]
        
        def fit_cb(pars, parmins, parmaxes, i):
            # freeze model parameter i
            (current_pars,
             current_parmins,
             current_parmaxes) = freeze_par(pars, parmins, parmaxes, i)
            
            fit_pars = fitfunc(statfunc, current_pars,
                               current_parmins,
                               current_parmaxes)[1]
            # If stat is not chi-squared, and fit method is
            # lmdif, need to recalculate stat at end, just
            # like in sherpa/sherpa/fit.py:fit()
            stat = statfunc(fit_pars)[0]
            #stat = fitfunc(scb, pars, parmins, parmaxes)[2]
            # thaw model parameter i
            thaw_par(i)
            return stat

        #
        # convert stat call back to have the same signature as fit call back
        #
        def stat_cb_extra_args( fcn ):
            def stat_cb_wrapper( x, *args ):
                return fcn( x )
            return stat_cb_wrapper

        statcb = stat_cb_extra_args( stat_cb )
        if 1 == len( pars ):
            fitcb = statcb
        else:
            fitcb = fit_cb

        return self._estfunc(pars, parmins, parmaxes, parhardmins,
                             parhardmaxes, self.sigma, self.eps,
                             self.tol, self.maxiters, self.remin,
                             self.verbose, limit_parnums,
                             statcb, fitcb, report_progress, get_par_name,
                             self.parallel, self.numcores, self.openinterval)

class Projection(EstMethod):

    # defined pre-instantiation for pickling
    _added_config = {'remin': 0.01,
                     'fast': False,
                     'parallel':True,
                     'numcores' : _ncpus,
                     'maxfits' : 5,
                     'max_rstat' : 3,
                     'tol' : 0.2}

    def __init__(self, name='projection'):
        EstMethod.__init__(self, name, projection)

        # Update EstMethod.config dict with Projection specifics
        self.config.update(self._added_config)

    def compute(self, statfunc, fitfunc, pars,
                parmins, parmaxes, parhardmins,
                parhardmaxes, limit_parnums, freeze_par, thaw_par,
                report_progress, get_par_name,
                statargs=(), statkwargs={}):

        def stat_cb(pars):
            return statfunc(pars)[0]
        
        def fit_cb(pars, parmins, parmaxes, i):
            # freeze model parameter i
            (current_pars,
             current_parmins,
             current_parmaxes) = freeze_par(pars, parmins, parmaxes, i)
            fit_pars = fitfunc(statfunc, current_pars,
                               current_parmins,
                               current_parmaxes)[1]
            # If stat is not chi-squared, and fit method is
            # lmdif, need to recalculate stat at end, just
            # like in sherpa/sherpa/fit.py:fit()
            stat = statfunc(fit_pars)[0]
            #stat = fitfunc(scb, pars, parmins, parmaxes)[2]
            # thaw model parameter i
            thaw_par(i)
            return stat

        return self._estfunc(pars, parmins, parmaxes, parhardmins,
                             parhardmaxes, self.sigma, self.eps,
                             self.tol,
                             self.maxiters, self.remin, limit_parnums,
                             stat_cb, fit_cb, report_progress, get_par_name,
                             self.parallel, self.numcores)

def covariance(pars, parmins, parmaxes, parhardmins, parhardmaxes, sigma, eps,
               tol, maxiters, remin, limit_parnums, stat_cb, fit_cb,
               report_progress):
    # Do nothing with tol
    # Do nothing with report_progress (generally fast enough we don't
    # need to report back per-parameter progress)
    
    # Even though we only want limits on certain parameters, we have to
    # compute the matrix for *all* thawed parameters.  So we will do that,
    # and then pick the parameters of interest out of the result.

    try:
        info = _est_funcs.info_matrix(pars, parmins, parmaxes, parhardmins,
                                      parhardmaxes, sigma, eps, maxiters,
                                      remin, stat_cb)
    except EstNewMin:
        # catch the EstNewMin exception and attach the modified
        # parameter values to the exception obj.  These modified
        # parvals determine the new lower statistic.
        raise EstNewMin(pars)
    except:
        raise

    # Invert matrix, take its square root and multiply by sigma to get
    # parameter uncertainties; parameter uncertainties are the
    # diagonal elements of the matrix.

    # Use simpler matrix inversion function from numpy.  If that
    # doesn't work, assume it's an ill-conditioned or singular matrix,
    # and call pinv from numpy -- pinv will call the SVD function to
    # invert the matrix.  But call pinv *only* when inv is shown not
    # to work in a particular case -- use inv by default!

    # The reason for this is that pinv can give back very strange
    # results, when you don't *need* to use pinv, and it *also*
    # happens that the ratio between smallest and largest diagonal
    # elements approaches the machine precision for the data type.
    # The result is that errors that come from the largest diagonal
    # element are ludicrously small; you can't have a parameter value
    # of order 1.0, and an error of order 10^-30, for example.  The
    # simpler inv function for inverting matrices does not appear to
    # have the same issue.

    invfunc = numpy.linalg.inv
    inv_info = None
    
    try:
        inv_info = invfunc(info)

    except numpy.linalg.linalg.LinAlgError:
        # catch the SVD exception and exit gracefully
        inv_info = numpy.zeros_like(info)
        inv_info[:] = numpy.nan

    except:
        # Compatibility with pre-0.9.8 numpy
        if hasattr(numpy.linalg, 'pinv'):
            invfunc = numpy.linalg.pinv
        else:
            invfunc = numpy.linalg.generalized_inverse

        try:
            inv_info = invfunc(info)
        except numpy.linalg.linalg.LinAlgError:
            # catch the SVD exception and exit gracefully
            inv_info = numpy.zeros_like(info)
            inv_info[:] = numpy.nan
    
    diag = (sigma * numpy.sqrt(inv_info)).diagonal()

    # limit_parnums lists the indices of the array pars, that
    # correspond to the parameters of interest.  We will pick out
    # the diagonal elements corresponding to entries in limits_parnums,
    # and return only those bounds to the user.
    upper_bounds = []
    lower_bounds = []
    error_flags = []
    for num in limit_parnums:
        eflag = est_success
        ubound = diag[num]
        lbound = -diag[num]
        if (pars[num] + ubound < parhardmaxes[num]):
            pass
        else:
            ubound = numpy.nan
            eflag = est_hardmax
        if (pars[num] + lbound > parhardmins[num]):
            pass
        else:
            lbound = numpy.nan
            if (eflag == est_hardmax):
                eflag = est_hardminmax
            else:
                eflag = est_hardmin
        upper_bounds.append(ubound)
        lower_bounds.append(lbound)
        error_flags.append(eflag)
        
    return (numpy.array(lower_bounds), numpy.array(upper_bounds),
            numpy.array(error_flags), 0, inv_info)


def projection(pars, parmins, parmaxes, parhardmins, parhardmaxes, sigma, eps,
               tol, maxiters, remin, limit_parnums, stat_cb, fit_cb,
               report_progress, get_par_name, do_parallel, numcores):
    i = 0                                 # Iterate through parameters
                                          #  to be searched on
    numsearched = len(limit_parnums)      # Number of parameters to be
                                          #  searched on (*not* number
                                          #  of thawed parameters, just
                                          #  number we are searching on
                                          #  (i.e., len(limit_parnums))
    lower_limits = numpy.array([])        # Lower limits for parameters
                                          #  searched on 
    upper_limits = numpy.array([])        # Upper limits for parameters
                                          #  searched on
    eflags = numpy.array([], numpy.int)   # Fail status after search for
                                          # each parameter
    nfits = 0                             # Total number of fits

    # _est_funcs.projection can be called on any subset of the thawed
    # parameters.  So we made a change here to call _est_funcs.projection
    # once per parameter number listed in limit_parnums, instead of
    # calling _est_funcs.projection once, with the original limit_parnums
    # array.  This way, we can report pack progress after the confidence
    # limit search is completed for each parameter, without descending
    # into the C++ code.
    #
    # It does mean we have to take apart the tuple returned by each call
    # to _est_funcs.projection; take the data we've pulled out, and
    # upon exiting the while loop, constructing a new tuple to return.
    # SMD 03/17/2009

    # Keep references to numpy.append, _est_funcs.projection, because
    # we call these functions every time through the loop.
    append = numpy.append
    proj_func = _est_funcs.projection

    def func(i, singleparnum, lock=None):
        try:
            singlebounds = proj_func(pars, parmins, parmaxes,
                                     parhardmins, parhardmaxes,
                                     sigma, eps, tol, maxiters,
                                     remin, [singleparnum], stat_cb,
                                     fit_cb)
        except EstNewMin:
            # catch the EstNewMin exception and attach the modified
            # parameter values to the exception obj.  These modified
            # parvals determine the new lower statistic.
            raise EstNewMin(pars)
        except:
            raise

        if lock is not None:
            lock.acquire()

        report_progress(singleparnum, singlebounds[0], singlebounds[1])

        if lock is not None:
            lock.release()

        return (singlebounds[0][0], singlebounds[1][0], singlebounds[2][0],
                singlebounds[3], None)

    if len(limit_parnums) < 2 or not _multi or numcores < 2:
        do_parallel = False

    if not do_parallel:
        append = numpy.append
        lower_limits = numpy.array([])
        upper_limits = numpy.array([])
        eflags = numpy.array([], numpy.int)
        nfits = 0
        for i in range(len(limit_parnums)):
            singlebounds = func(i, limit_parnums[i])
            lower_limits = append(lower_limits, singlebounds[0])
            upper_limits = append(upper_limits, singlebounds[1])
            eflags = append(eflags, singlebounds[2])
            nfits = nfits + singlebounds[3]
        return (lower_limits, upper_limits, eflags, nfits, None)

    return parallel_est(func, limit_parnums, pars, numcores)

#################################confidence###################################

class ConfArgs( object ):
    """The class ConfArgs is responsible for the arguments to the fit
    call back function."""

    def __init__( self, xpars, smin, smax, hmin, hmax, target_stat ):
        self.ith_par = 0
        self.xpars = numpy.array( xpars, copy=True )
        self.slimit = ( numpy.array( smin, copy=True ),
                        numpy.array( smax, copy=True ) )
        self.hlimit = ( numpy.array( hmin, copy=True ),
                        numpy.array( hmax, copy=True ) )        
        self.target_stat = target_stat

    def __call__( self ):
        return ( self.ith_par, self.xpars, self.slimit, self.hlimit,
                 self.target_stat )

    def __str__( self ):
        a2s = numpy.array2string
        msg = ''
        msg += '# smin = ' + a2s(self.slimit[0],precision=6) + '\n'
        msg += '# smax = ' + a2s(self.slimit[1],precision=6) + '\n'
        msg += '# hmin = ' + a2s(self.hlimit[0],precision=6) + '\n'
        msg += '# hmax = ' + a2s(self.hlimit[1],precision=6) + '\n#\n'
        msg += '# Note: for the intermediate steps, the notation:\n'
        msg += '         par.name -/+: f( x ) = stat\n'
        msg += '# ==> `stat` is the statistic when parameter `par.name` is frozen at `x`\n'
        msg += '# while searching for the `lower/upper` confidence level, repectively.\n#'
        return msg

    def __rep__( self ):
        return ( "<%s ConfArgs method instance'%s'>" %
                 (type(self).__name__, self.name) )

    def get_par( self ):
        """return the current (worked on) par"""
        return self.xpars[ self.ith_par ]

    def get_hlimit( self, dir ):
        """ return the current (worked on) hard limit"""        
        return self.hlimit[ dir ][ self.ith_par ]

    def get_slimit( self, dir ):
        """ return the current (worked on) soft limit"""
        return self.slimit[ dir ][ self.ith_par ]    


class ConfBlog( object ):

    def __init__( self, blogger, prefix, verbose, lock, debug=False ):
        self.blogger = blogger
        self.prefix = prefix
        self.verbose = verbose
        self.lock = lock
        self.debug = debug

    def __str__( self ):
        return 'ConfBlog::__str__( )'

    def __rep__( self ):
        return 'ConfBlog::__rep__( )'


class ConfBracket( object ):
    """The class ConfBracket is reponsible for bracketing the root within
    the interval (a,b) where f(a)*f(b) < 0.0"""
    
    neg_pos = ( -1, 1 )

    class Limit( object ):
        def __init__( self, limit ):
            self.limit = limit

    class LowerLimit( Limit ):
        def __init__( self, limit ):
            ConfBracket.Limit.__init__( self, limit )
        def __str__( self ):
            str = 'LowerLimit: limit=%e' % self.limit            
            return str
        def is_beyond_limit( self, x ):
            if x < self.limit:
                return True
            else:
                return False

    class UpperLimit( Limit ):
        def __init__( self, limit ):
            ConfBracket.Limit.__init__( self, limit )
        def __str__( self ):
            str = 'UpperLimit: limit=%e' % self.limit
            return str
        def is_beyond_limit( self, x ):
            if x > self.limit:
                return True
            else:
                return False
        
    def __init__( self, myargs, trial_points ):
        self.myargs = myargs
        self.trial_points = trial_points
        self.fcn = None
        
    def __repr__( self ):
        return ("<%s Bracket error-estimation method instance '%s'>" %
                (type(self).__name__, self.name))

    def __call__( self, dir, iter, step_size, open_interval, maxiters, tol,
                  bloginfo ):

        #
        # Either 1) a root has been found (ConfRootZero), 2) an interval
        # where the root has been confined (ConfRootBracket) or 3) No
        # possible chance for a root (ConfRootNone), ie by trying points
        # upto/beyond the hard limit and no chance for a root has been found.
        #
        find = trace_fcn( self.find, bloginfo )
        return find( dir, iter, step_size, open_interval, maxiters, tol,
                     bloginfo )

    def find( self, dir, iter, step_size, open_interval, maxiters, tol,
              bloginfo, base=2.0 ):

        assert self.fcn != None, 'callback func has not been set'
            
        hlimit = [ ConfBracket.LowerLimit( self.myargs.get_hlimit( dir ) ),
                   ConfBracket.UpperLimit( self.myargs.get_hlimit( dir ) ) ]
        slimit = [ ConfBracket.LowerLimit( self.myargs.get_slimit( dir ) ),
                   ConfBracket.UpperLimit( self.myargs.get_slimit( dir ) ) ]

        xxx = self.trial_points[ 0 ]
        fff = self.trial_points[ 1 ]
        
        conf_step = ConfStep( xxx, fff )

        mymaxiters = maxiters
        if mymaxiters > 16:
            mymaxiters = 16

        plateau = 0
        max_plateau_iter = 5

        try:
            
            for iter in range( mymaxiters ):

                if 0 == iter:
                    x = conf_step.covar( dir, iter, step_size, base )
                elif 1 == iter:
                    x = conf_step.secant( dir, iter, step_size, base )
                    #x = conf_step.covar( dir, iter, step_size, base )
                else:
                    x = conf_step.quad( dir, iter, step_size, base, bloginfo )

                if x is None or numpy.isnan( x ):
                    return ConfRootNone( )

                # Make sure x is not beyond the **hard** limit
                if hlimit[ dir ].is_beyond_limit( x ):
                    
                    x = hlimit[ dir ].limit
                    f = self.fcn( x, self.myargs( ) )
                    #print 'find(): beyond hard limit: f(%.14e)=%.14e' % (x,f)
                    if abs( f ) <= tol:
                        return ConfRootZero( x )
                    if f >= 0.0:
                        return ConfRootBracket( self.fcn, self.trial_points,
                                                open_interval )
                    else:
                        return ConfRootNone( )

                elif slimit[ dir ].is_beyond_limit( x ):
                    
                    f = self.fcn( x, self.myargs( ) )
                    #print 'find(): beyond soft limit: f(%.14e)=%.14e' % (x,f)
                    if abs( f ) <= tol:
                        return ConfRootZero( x )
                    if f >= 0.0:
                        return ConfRootBracket( self.fcn, self.trial_points,
                                                open_interval )
                    elif f < fff[ -2 ]:
                        # if the fit beyond the soft limit is a better fit
                        # then the confidence for the parameter does not exist
                        return ConfRootNone( )

                else:

                    f = self.fcn( x, self.myargs( ) )
                    #print 'find(): f(%.14e)=%.14e' % (x,f)            
                    if abs( f ) <= tol:
                        return ConfRootZero( x )
                    elif f >= 0.0:
                        return ConfRootBracket( self.fcn, self.trial_points,
                                                open_interval )


                if Knuth_close( fff[-2], fff[-1], 1.0e-6 ):
                    plateau += 1
                    if plateau > max_plateau_iter:
                        #print 'find( %d ): plateau = %d', (iter,plateau)
                        return ConfRootNone( None )
                #else:
                #    if plateau > 0:
                #        plateau -= 1
                    
            return ConfRootNone( None )

        except OutOfBoundErr:
            return ConfRootNone( )



class ConfRootNone( object ):
    """The base class for the root of the confidence interval"""
    
    def __init__( self, root=None ):
        """If self.root == None, then 
        1) points up to the hard limits were tried and it was not possible
        to bracketed the solution.
        2) a parameter beyond the soft limit has been tried and the new stat
        was found to be **less** then the initial minimum"""
        self.root = root

    def __call__( self, tol, bloginfo ):
        return self.root

    def __str__( self ):
        return 'No possible root exist'

class ConfRootBracket( ConfRootNone ):
    """The class contains the bracket where the confidence root has
    been bracketed, ie where f(a)*f(b) < 0"""
    
    def __init__( self, fcn, trial_points, open_interval ):
        ConfRootNone.__init__( self, None )
        self.fcn = fcn
        self.trial_points = trial_points
        self.open_interval = open_interval

    def __call__( self, tol, bloginfo ):

        def warn_user_about_open_interval( list ):
            
            if bloginfo.lock is not None:
                bloginfo.lock.acquire()

            if 0 == bloginfo.verbose:
                prefix = '%s ' % bloginfo.prefix.lstrip()
            else:
                prefix = '%s ' % bloginfo.prefix

            interval = list_to_open_interval( list )
            bloginfo.blogger.info( prefix +
                                   'WARNING: The confidence level lies within '
                                   + interval )

            if bloginfo.lock is not None:
                bloginfo.lock.release()

            return
        
        xxx = self.trial_points[ 0 ]
        fff = self.trial_points[ 1 ]

        if mysgn( fff[ -2 ] ) == mysgn( fff[ -1 ] ):
            self.root = None
            return None

        myzeroin = trace_fcn( zeroin, bloginfo )
        answer = myzeroin( self.fcn, xxx[ -2 ], xxx[ -1 ], fa=fff[ -2 ],
                           fb=fff[ -1 ], maxfev=32, tol=tol )
        if abs( answer[0][1] ) > tol:
            xafa = answer[ 1 ][ 0 ]
            xa = xafa[ 0 ]
            fa = xafa[ 1 ]
            xbfb = answer[ 1 ][ 1 ]
            xb = xbfb[ 0 ]
            fb = xbfb[ 1 ]
            if mysgn( fa ) != mysgn( fb ):
                if False == self.open_interval:
                    warn_user_about_open_interval( [ xa, xb ] )
                    return ( xa + xb ) / 2.0
                else:
                    if xa < xb:
                        return ( xa, xb )
                    else:
                        return ( xb, xa )
            else:
                return None
        self.root = answer[ 0 ][ 0 ]
        return self.root
    
    def __str__( self ):
        str = 'root is within the interval ( f(%e)=%e, f(%e)=%e )' \
              % ( self.trial_points[ 0 ][ -2 ], self.trial_points[ 1 ][ -2 ],
                  self.trial_points[ 0 ][ -1 ], self.trial_points[ 1 ][ -1 ], )
        return str
    
class ConfRootZero( ConfRootNone ):
    """The class with the root/zero of the confidence interval"""
    
    def __init__( self, root ):
        ConfRootNone.__init__( self, root )

    def __str__( self ):
        str = 'root = %e' % self.root
        return str


class ConfStep( object ):

    def __init__( self, xtrial, ftrial ):
        self.xtrial = xtrial
        self.ftrial = ftrial

    def covar( self, dir, iter, stepsize, base ):
        return self.xtrial[ -1 ] + \
               ConfBracket.neg_pos[ dir ] * pow( base, iter ) * stepsize
        #return self.xtrial[ 0 ] + \
        #       ConfBracket.neg_pos[ dir ] * pow( base, iter ) * stepsize

    def Halley( self, coeffs, x, maxfev=8, tol=1.0e-3 ):

        for nfev in range( maxfev ):

            ax = coeffs[ 0 ] * x
            fval = ( ax + coeffs[ 1 ] ) * x + coeffs[ 2 ]
            if abs( fval ) <= tol:
                return [x, fval]

            fdif = 2.0 * ax + coeffs[ 1 ]
            fdif2 = 2.0
        
            numer = 2.0 * fval * fdif
            denom = 2.0 * fdif * fdif - fval * fdif2
            x -= numer / denom
            nfev += 1
            
        return [x, fval]

    def is_same_dir( self, dir, current_pos, proposed_pos ):
        delta = proposed_pos - current_pos
        return mysgn( delta ) == ConfBracket.neg_pos[ dir ]
    
    def quad( self, dir, iter, step_size, base, bloginfo ):

        coeffs = quad_coef( self.xtrial[ -3: ], self.ftrial[ -3: ] )
        delta = ConfBracket.neg_pos[ dir ]
        delta *= abs( self.xtrial[ -1 ] - self.xtrial[ -2 ] )
        lastx = self.xtrial[ -1 ]

        mroot = demuller( numpy.poly1d( coeffs ), lastx + delta,
                          lastx + 2 * delta, lastx + 3 * delta,
                          tol=1.0e-2 )
        
        xroot = mroot[ 0 ][ 0 ]
        if xroot is None or numpy.isnan( xroot ):
            return self.covar( dir, iter, step_size, base )
        
        try:
            Halley = trace_fcn( self.Halley, bloginfo )
            [xroot, froot] = Halley( coeffs, xroot, tol=1.0e-3 )
        except ZeroDivisionError:
            xroot = None

        if (None != xroot and False == numpy.isnan( xroot )) and \
               self.is_same_dir( dir, self.xtrial[ -1 ], xroot ):
            return xroot
        else:
            return self.covar( dir, iter, step_size, base )


    def secant( self, dir, iter, step_size, base ):

        xb = self.xtrial[ -2 ]
        fb = self.ftrial[ -2 ]
        xa = self.xtrial[ -1 ]
        fa = self.ftrial[ -1 ]

        if abs( fb ) > abs( fa ) or 0.0 == fa:
            return self.covar( dir, iter, step_size, base )

        s = fb / fa
        p = ( xa - xb ) * s
        if 1.0 == s:
            return self.covar( dir, iter, step_size, base )
        q = 1.0 - s
        x = xb - p / q

        if self.is_same_dir( dir, xa, x ):
            return x
        else:
            return self.covar( dir, iter, step_size, base )            


def trace_fcn( fcn, bloginfo ):

    if False == bloginfo.debug:
        return fcn

    from itertools import chain
    def echo( *args, **kwargs ):
        '''compact but more details then debugger'''
        name = fcn.__name__
        str = '%s%s(%s)' % (bloginfo.prefix,fcn.__name__, ", ".join(map(repr, chain(args, kwargs.values()))))
        bloginfo.blogger.info( str )
        return fcn( *args, **kwargs )
    
    def debugger( *args, **kwargs ):
        str = '%s%s( ' % ( bloginfo.prefix, fcn.__name__ )
        if len( args ) > 1:
            str += args[ 0 ].__str__( )
        for arg in args[ 1: ]:
            str = '%s, %s' % ( str, arg )
        for key in kwargs.iterkeys( ):
            value = kwargs[ key ]
            str = '%s, %s=%s' % ( str, key, value )
        val = fcn( *args, **kwargs )
        str += ' )  %s ' % val
        bloginfo.blogger.info( str )
        return val
    
    return debugger

def confidence(pars, parmins, parmaxes, parhardmins, parhardmaxes, sigma, eps,
               tol, maxiters, remin, verbose, limit_parnums, stat_cb,
               fit_cb, report_progress, get_par_name, do_parallel, numcores,
               open_interval):

    def get_prefix( index, name, minus_plus ):
        '''To print the prefix/indent when verbose is on'''
        prefix = [[],[]]
        blank = 3 * index * ' '
        for dir in range( 2 ):
            prefix[ dir ] = blank + name + ' ' + minus_plus[ dir ] + ':'
        return prefix

    def get_delta_root( arg, dir, par_at_min ):

        my_neg_pos = ConfBracket.neg_pos[ dir ]

        if is_iterable( arg ):
            return arg
            #return map( lambda x: my_neg_pos * abs( x - par_at_min ), arg )
        elif None != arg:
            arg -= par_at_min
            return my_neg_pos * abs( arg )
        else:
            return arg
        
    def get_step_size( error_scales, upper_scales, index, par ):
        
        if 0 != error_scales[ index ]:
            # if covar error is NaN then set it to fraction of the par value.
            ith_covar_err = 0.0625 * abs( par )
        else:
            ith_covar_err = abs( upper_scales[ index ] )
        if 0.0 == ith_covar_err:
            # just in case covar and/or par is 0
            ith_covar_err = 1.0e-6
            
        return ith_covar_err

    def monitor_func( fcn, history ):
        def myfunc( x, *args ):
            fval = fcn( x, *args )
            history[ 0 ].append( x )
            history[ 1 ].append( fval )
            return fval
        return myfunc

    def print_status( myblog, verbose, prefix, answer, lock ):

        if lock is not None:
            lock.acquire()
        
        if 0 == verbose:
            msg = '%s\t' % prefix.lstrip()
        else:
            msg = '%s\t' % prefix

        if is_iterable( answer ):
            msg += list_to_open_interval( answer )
        elif answer is None:
            msg += '-----'
        else:
            msg += '%g' % answer
        myblog( msg )

        if lock is not None:
            lock.release()

    #
    # Work in the translated coordinate. Hence the 'errors/confidence'
    # are the zeros/roots in the translated coordinate system.
    #
    def translated_fit_cb( fcn, myargs ):
        def translated_fit_cb_wrapper( x, *args ):
            hlimit = myargs.hlimit
            slimit = myargs.slimit
            hmin = hlimit[ 0 ]
            hmax = hlimit[ 1 ]
            xpars = myargs.xpars
            ith_par = myargs.ith_par
            # The parameter must be within the hard limits
            if x < hmin[ ith_par ] or x > hmax[ ith_par ]:
                raise OutOfBoundErr
            smin = slimit[ 0 ]
            smax = slimit[ 1 ]
            orig_ith_xpar = xpars[ ith_par ]
            xpars[ ith_par ] = x
            translated_stat = fcn( xpars, smin, smax, ith_par ) - myargs.target_stat
            xpars[ ith_par ] = orig_ith_xpar
            return translated_stat
        return translated_fit_cb_wrapper

    def verbose_fitcb( fcn, bloginfo ):
        if 0 == bloginfo.verbose:
            return fcn
        def verbose_fcn( x, *args ):
            fval = fcn( x, *args )
            str = '%s f( %e ) =' % ( bloginfo.prefix, x )
            if fval is None:
                str = '%s None' % str                
            else:
                str = '%s %e' % ( str, fval )
            bloginfo.blogger.info( str )
            return fval
        return verbose_fcn

    sherpablog = logging.getLogger( 'sherpa' ) # where to print progress report

    # Get minimum fit statistic, and calculate target statistic value
    orig_min_stat = stat_cb(pars)
    delta_stat = sigma * sigma
    target_stat = orig_min_stat + delta_stat

    lower_scales = None
    upper_scales = None
    error_scales = None
    nfits = 0
    results = None

    try:
        (lower_scales, upper_scales, error_scales, nfits,
         results) = covariance( pars, parmins, parmaxes, parhardmins,
                                parhardmaxes, 1.0, eps, tol, maxiters,
                                remin, limit_parnums, stat_cb,
                                fit_cb, report_progress )
    except EstNewMin, e:
        raise e
    except:
        error_scales = numpy.array( len( pars ) * [ est_hardminmax ] )

    debug = False                                 # for internal use only

    myargs = ConfArgs( pars, parmins, parmaxes, parhardmins, parhardmaxes,
                       target_stat )

    if 0 != verbose:
        msg = '#\n# f' + numpy.array2string(numpy.asarray(pars),precision=6)
        msg += ' = %e\n' % orig_min_stat
        msg += '# sigma = %e\n' % sigma
        msg += '# target_stat = %e\n' % target_stat
        msg += '# tol = %e\n' % eps
        msg += '%s' % myargs
        sherpablog.info( msg )
        
    dict = { }

    def func( counter, singleparnum, lock=None ):

        # nfev contains the number of times it was fitted
        nfev, counter_cb = func_counter( fit_cb )

        #
        # These are the bounds to be returned by this method
        #
        conf_int = [ [], [] ]
        error_flags = []

        #
        # If the user has requested a specific parameter to be 
        # calculated then 'ith_par' represents the index of the 
        # free parameter to deal with.
        #
        myargs.ith_par = singleparnum

        fitcb = translated_fit_cb( counter_cb, myargs )
        
        par_name = get_par_name( myargs.ith_par )

        ith_covar_err = get_step_size( error_scales, upper_scales, counter,
                                       pars[ myargs.ith_par ] )

        trial_points = [ [  ], [  ] ]
        fitcb = monitor_func( fitcb, trial_points )

        bracket = ConfBracket( myargs, trial_points )

        # the parameter name is set, may as well get the prefix
        prefix = get_prefix( counter, par_name, ['-', '+' ] )
        
        
        myfitcb = [ verbose_fitcb( fitcb,
                                   ConfBlog(sherpablog,prefix[0],verbose,lock) ),
                    verbose_fitcb( fitcb,
                                   ConfBlog(sherpablog,prefix[1],verbose,lock) ) ]
        
        for dir in range( 2 ):

            #
            # trial_points stores the history of the points for the
            # parameter which has been evaluated in order to locate
            # the root. Note the first point is 'given' since the info
            # of the minimum is crucial to the search.
            #
            bracket.trial_points[0].append( pars[ myargs.ith_par ] )
            bracket.trial_points[1].append( - delta_stat )

            myblog = ConfBlog( sherpablog, prefix[ dir ], verbose, lock,
                               debug )

            # have to set the callback func otherwise disaster.
            bracket.fcn = myfitcb[ dir ]
            root = bracket( dir, iter, ith_covar_err, open_interval, maxiters,
                            eps, myblog )

            myzero = root( eps, myblog )

            delta_zero = get_delta_root( myzero, dir, pars[ myargs.ith_par ] )
                
            conf_int[ dir ].append( delta_zero )

            status_prefix = get_prefix( counter, par_name, ['lower bound',
                                                            'upper bound' ] )
            print_status( myblog.blogger.info, verbose, status_prefix[ dir ],
                          delta_zero, lock )

        error_flags.append( est_success )

        #
        # include the minimum point to seperate the -/+ interval
        #
        dict[ par_name ] = trial_points

        return ( conf_int[ 0 ][0], conf_int[ 1 ][0], error_flags[0],
                 nfev[0], None )

    if len(limit_parnums) < 2 or not _multi or numcores < 2:
        do_parallel = False

    if not do_parallel:
        lower_limits = []
        upper_limits = []
        eflags = []
        nfits = 0
        for i in range(len(limit_parnums)):
            lower_limit, upper_limit, flags, nfit, extra = func(i, limit_parnums[i])
            lower_limits.append(lower_limit)
            upper_limits.append(upper_limit)
            eflags.append(flags)
            nfits += nfit
        return (lower_limits, upper_limits, eflags, nfits, None)

    return parallel_est(func, limit_parnums, pars, numcores)

#################################confidence###################################

def parallel_est(estfunc, limit_parnums, pars, numcores=_ncpus):

    tasks = []

    def worker(out_q, err_q, parids, parnums, parvals, lock):
        results = []
        for parid, singleparnum in izip(parids, parnums):
            try:
                result = estfunc(parid, singleparnum, lock)
                results.append( (parid, result) )
            except EstNewMin:
                # catch the EstNewMin exception and include the exception
                # class and the modified parameter values to the error queue.
                # These modified parvals determine the new lower statistic.
                # The exception class will be instaniated re-raised with the
                # parameter values attached.  C++ Python exceptions are not
                # picklable for use in the queue.
                err_q.put( EstNewMin(parvals) )
                return
            except Exception, e:
                #err_q.put( e.__class__() )
                err_q.put(e)
                return

        out_q.put(results)

    # The multiprocessing manager provides references to process-safe
    # shared objects like Queue and Lock
    manager = multiprocessing.Manager()
    out_q = manager.Queue()
    err_q = manager.Queue()
    lock  = manager.Lock()

    size = len(limit_parnums)
    parids = numpy.arange(size)

    # if len(limit_parnums) is less than numcores, only use length number of 
    # processes
    if size < numcores:
        numcores = size  

    # group limit_parnums into numcores-worth of chunks
    limit_parnums = numpy.array_split(limit_parnums, numcores)
    parids = numpy.array_split(parids, numcores)

    tasks = [multiprocessing.Process(target=worker,
                       args=(out_q, err_q, parid, parnum, pars, lock))
             for parid, parnum in izip(parids, limit_parnums)]

    return run_tasks(tasks, out_q, err_q, size)


def run_tasks(tasks, out_q, err_q, size):

    die = (lambda tasks : [task.terminate() for task in tasks
                           if task.exitcode is None])

    try:
        for task in tasks:
            task.start()

        for task in tasks:
            task.join()

    except KeyboardInterrupt, e:
        # kill all slave processes on ctrl-C
        die(tasks)
        raise e

    if not err_q.empty():
        die(tasks)
        raise err_q.get()

    lower_limits = size*[None]
    upper_limits = size*[None]
    eflags = size*[None]
    nfits = 0

    while not out_q.empty():
        for parid, singlebounds in out_q.get():
            # Have to guarantee that the tuple returned by projection
            # is always (array, array, array, int) for this to work.
            lower_limits[parid] = singlebounds[0]
            upper_limits[parid] = singlebounds[1]
            eflags[parid] = singlebounds[2]
            nfits += singlebounds[3]

    return (lower_limits, upper_limits, eflags, nfits, None)
