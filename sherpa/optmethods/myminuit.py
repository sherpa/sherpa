#
# https://github.com/scikit-hep/iminuit
#
import numpy
from iminuit import Minuit
from sherpa.optmethods.optfcts import EPSILON, _check_args, _get_saofit_msg
from sherpa.models.parameter import hugeval

__all__ = ('myminuit',)


def myminuit(fcn, x0, xmin, xmax, tol=EPSILON, maxfev=None, leastsqr=True,
             migrad=True, verbose=0):
    """if migrad == False then use Minuit simplex optimization"""

    def stat_cb0(pars):
        return fcn(pars)[0]

    x0, xmin, xmax = _check_args(x0, xmin, xmax)

    bounds = []
    x_tmp = None
    y_tmp = None
    myhugeval = hugeval / 10.0
    for xxx, yyy in zip(xmin, xmax):
        x_tmp = xxx
        # +/- inf to preven iminuit from transforming
        if xxx < - myhugeval:
            x_tmp = - numpy.inf
        y_tmp = yyy
        if yyy > myhugeval:
            y_tmp = numpy.inf
        bounds.append((x_tmp, y_tmp))

    midnight =  Minuit(stat_cb0, numpy.asarray(x0))
    midnight.limits = bounds
    midnight.errors = tol
    midnight.print_level = verbose
    if leastsqr:
        midnight.errordef = Minuit.LEAST_SQUARES
    else:
        midnight.errordef = Minuit.LIKELIHOOD

    if migrad:
        if maxfev is None:
            maxfev = 128 * len(x0)
        midnight.migrad(ncall=maxfev)
    else:
        if maxfev is None:
            maxfev = 512 * len(x0)
        midnight.simplex(ncall=maxfev)

    nfev = midnight.nfcn
    fval = midnight.fval
    # x = numpy.asarray([ val for val in midnight.values ])
    x = midnight.values
    
    if midnight.valid:
        status, msg = _get_saofit_msg(maxfev, 0)
        if midnight.accurate:
            pass
        else:
            msg += ", but uncertainties are unrealiable."
    else:
        fmin = midnight.fmin
        status, msg = _get_saofit_msg(maxfev, 4)
        if fmin.has_reached_call_limit:
            status, msg = _get_saofit_msg(maxfev, 3)
        if fmin.is_above_max_edm:
            msg += " Estimated distance to minimum too large."

    npar = len(x)
    if midnight.covariance is not None:
        covar = midnight.covariance
        rv = (status, x, fval, msg, {'nfev': nfev,
                                     'covar': covar, 'info': 1})
    else:
        rv = (status, x, fval, msg, {'nfev': nfev, 'info': 1})

    return rv
