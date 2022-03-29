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
    myhugeval = hugeval / 10.0
    for xxx, yyy in zip(xmin, xmax):
        x_tmp = xxx
        # +/- inf to prevent iminuit from transforming
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

    if midnight.valid:
        status, msg = _get_saofit_msg(maxfev, 0)
        if not midnight.accurate:
            msg += ", but uncertainties are unrealiable."

    else:
        fmin = midnight.fmin
        status, msg = _get_saofit_msg(maxfev, 4)
        if fmin.has_reached_call_limit:
            status, msg = _get_saofit_msg(maxfev, 3)
        if fmin.is_above_max_edm:
            msg += " Estimated distance to minimum too large."

    rvopt = {'nfev': midnight.nfcn, 'info': 1}
    rv = (status, midnight.values, midnight.fval, msg, rvopt)

    if midnight.covariance is not None:
        rvopt['covar'] = midnight.covariance

    return rv
