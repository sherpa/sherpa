# 
#  Copyright (C) 2008  Smithsonian Astrophysical Observatory
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
warning = logging.getLogger(__name__).warning
import numpy
from _utils import *
from _pileup import *
from itertools import izip
from sherpa.utils import SherpaFloat, get_position, filter_bins
from sherpa.utils.err import IOErr, DataErr

__all__ = ['arf_fold', 'rmf_fold', 'do_group', 'apply_pileup',
           'eqwidth', 'calc_photon_flux', 'calc_energy_flux',
           'calc_data_sum', 'calc_model_sum', 'shrink_effarea',
           'calc_data_sum2d','calc_model_sum2d', 'filter_resp',
           'calc_source_sum', 'compile_energy_grid',
           'expand_grouped_mask','resp_init', 'is_in', 'get_xspec_position']

try:
    from _region import *
    __all__.append('region_mask')
    __all__.append('Region')

except ImportError:
    warning('failed to import sherpa.astro.utils._region; Region routines ' +
            'will not be available')

_hc = 12.39841874  # nist.gov in [keV-Angstrom]

def get_xspec_position(y, x, xhi=None):
    if xhi is not None:
        if x[0] > x[-1] and xhi[0] > xhi[-1]:
            lo = _hc / xhi
            hi = _hc / x
            x = lo; xhi = hi
    else:
        if x[0] > x[-1]:
            x = _hc / x
    return get_position(y, x, xhi)


def compile_energy_grid(arglist):
    elo = numpy.unique(numpy.concatenate([indep[0] for indep in arglist]))
    ehi = numpy.unique(numpy.concatenate([indep[1] for indep in arglist]))

    in_elo = numpy.setdiff1d(elo,ehi)
    in_ehi = numpy.setdiff1d(ehi,elo)
    if len(in_elo) > 1:
        ehi = numpy.concatenate((ehi, in_elo[-(len(in_elo)-1):]))
    if len(in_ehi) > 1:
        elo = numpy.concatenate((elo, in_ehi[0:len(in_ehi)-1]))

    # FIXME since numpy.unique calls sort() underneath, may not need to
    # sort again here...
    elo.sort()
    ehi.sort()

    # determine index intervals using binary search in large src model
    # for n_th ARF energy bounds and populate a table to use later for
    # integration over said intervals
    htable = [(numpy.searchsorted(elo, arg[0]).astype(numpy.int32),
               numpy.searchsorted(ehi, arg[1]).astype(numpy.int32))
              for arg in arglist]

    return [elo, ehi, htable]

def bounds_check(lo, hi):
    if lo is not None and hi is not None and lo > hi:
        raise IOErr('boundscheck', lo, hi)

    if lo is not None and hi is None:
        hi = lo
        
    if hi is not None and lo is None:
        lo = hi

    return (lo, hi)


_charge_e = 1.60217653e-09 #elementary charge [ergs] 1 keV, nist.gov

def _flux( data, lo, hi, src, eflux=False, srcflux=False):
    lo, hi = bounds_check(lo, hi)

    axislist = None
    if hasattr(data, '_get_indep'):
        axislist = data._get_indep(filter=False)
    else:
        axislist = data.get_indep(filter=False)

    y = src(*axislist)

    if srcflux and len(axislist) > 1:
        y /= numpy.asarray(axislist[1]-axislist[0])

    dim = numpy.asarray(axislist).squeeze().ndim
    if eflux:
        # for energy flux, the sum of grid below must be in keV.
        energ = []
        for axis in axislist:
            grid = axis
            if hasattr(data, 'units') and data.units == 'wavelength':
                grid = data._hc/grid
            energ.append(grid)

        if dim == 1:
            y = numpy.asarray(0.5 * y * energ[0], SherpaFloat)
        elif dim == 2:
            y = numpy.asarray(0.5 * y * (energ[0]+energ[1]),
                            SherpaFloat)
        else:
            raise IOErr('>axes', "2")

    mask = filter_bins( (lo,), (hi,), (axislist[0],) )

    val = y.sum()
    if mask is not None:
        flux = y[mask]
        # flux density at a single bin -> divide by bin width.
        if dim == 2 and len(flux) == 1:
            flux /= numpy.abs(axislist[1][mask]-axislist[0][mask])
        val = flux.sum()

    if eflux:
        val *= _charge_e

    return val


def _counts( data, lo, hi, func, *args):
    lo, hi = bounds_check(lo, hi)
    old_filter = data.filter
    old_mask = data.mask
    old_quality_filter = getattr(data, 'quality_filter', None)
    try:
        data.notice()  # save and clear filter
        data.filter=None
        # filter_rsp = getattr(data, 'notice_response', None)
        # if filter_rsp is not None:
        #     filter_rsp(False)

        data.notice(lo, hi)
        counts = func(*args).sum()
        data.notice()
    finally:
        data.filter = old_filter
        data.mask = old_mask
        if old_quality_filter is not None:
            data.quality_filter = old_quality_filter

    return counts

def _counts2d(data, reg, func, *args):
    old_filter = data.filter
    old_mask = data.mask
    coord = getattr(data, 'coord', None)
    old_region = getattr(data, '_region', None)
    try:
        data.notice2d()  # save and clear filter

        data.notice2d(reg)
        counts = func(*args).sum()
        data.notice2d()
    finally:
        data.filter = old_filter
        data.mask = old_mask

        if coord is not None:
            data.set_coord(coord)

        if old_region is not None:
            data._region = old_region

    return counts

def calc_energy_flux( data, src, lo=None, hi=None):
    return _flux(data, lo, hi, src, eflux=True)

def calc_photon_flux( data, src, lo=None, hi=None):
    return _flux(data, lo, hi, src)

def calc_source_sum( data, src, lo=None, hi=None):
    return _flux(data, lo, hi, src, srcflux=True)

#def calc_source_sum2d( data, src, reg=None):
#    return _counts2d(data, reg, data.eval_model_to_fit, src)

def calc_data_sum(data, lo=None, hi=None):
    return _counts( data, lo, hi, data.apply_filter, data.get_dep() )

def calc_data_sum2d(data, reg=None):
    return _counts2d(data, reg, data.apply_filter, data.get_dep() )

def calc_model_sum(data, model, lo=None, hi=None):
    return _counts(data, lo, hi, data.eval_model_to_fit, model)

def calc_model_sum2d(data, model, reg=None):
    return _counts2d(data, reg, data.eval_model_to_fit, model)

def eqwidth(data, model, combo, lo=None, hi=None):

    lo, hi = bounds_check(lo, hi)

    my = None
    cy = None
    xlo = None
    xhi = None
    num = None
    eqw = 0.0
    if hasattr(data, 'get_response'):
        xlo, xhi = data._get_indep(filter=False)
        my = model(xlo,xhi)
        cy = combo(xlo,xhi)
        num = len(xlo)
    else:
        my = data.eval_model_to_fit( model )
        cy = data.eval_model_to_fit( combo )
        xlo = data.get_indep(filter=True)[0]
        num = len(xlo)

    mask = filter_bins( (lo,), (hi,), (xlo,) )
    if mask is not None:
        my = my[mask]
        cy = cy[mask]
        xlo = xlo[mask]
        num = len(xlo)

    for ebin, val in enumerate(xlo):
        if ebin < (num-1):
            eave = numpy.abs(xlo[ebin+1] - xlo[ebin])
        else:
            eave = numpy.abs(xlo[ebin-1] - xlo[ebin])
        if my[ebin] != 0.0:
            eqw += eave*(cy[ebin]-my[ebin])/my[ebin]

    return eqw


def calc_kcorr(data, model, z, obslo, obshi, restlo=None, resthi=None):

    if restlo is None:
        restlo = obslo
    if resthi is None:
        resthi = obshi

    if numpy.isscalar(z):
        z = numpy.array([z], dtype=float)
    else:
        z = numpy.asarray(z)

    if( 0 != sum(z[z<0]) ):
        raise IOErr('z<=0')

    if( obslo <= 0 or restlo <=0 or obshi <= obslo or resthi <= restlo ):
        raise IOErr('erange')

    if hasattr(data, 'get_response'):
        arf, rmf = data.get_response()
        elo = data.bin_lo
        ehi = data.bin_hi
        if arf is not None:
            elo = arf.energ_lo
            ehi = arf.energ_hi
        elif rmf is not None:
            elo = rmf.energ_lo
            ehi = rmf.energ_hi
    else:
        elo,ehi = data.get_indep()

    if elo is None or ehi is None:
        raise DataErr('noenergybins', data.name)

    emin = elo[0]
    emax = ehi[-1]

    if( restlo < emin or resthi > emax ):
        raise IOErr('energoverlap', emin, emax, 'rest-frame',
                    restlo, resthi, '')

    if( obslo*(1.0+z.min()) < emin ):
        raise IOErr('energoverlap', emin, emax, 'observed-frame',
                    restlo, resthi, "at a redshift of %f" % z.min())

    if( obshi*(1.0+z.max()) > emax ):
        raise IOErr('energoverlap', emin, emax, 'rest-frame',
                    restlo, resthi, "at a redshift of %f" % z.min())

    zplus1 = z + 1.0
    flux_rest = _flux( data, restlo, resthi, model, eflux=True)
    obs = numpy.asarray([_flux(data, obslo*zz, obshi*zz, model, eflux=True)
                         for zz in zplus1], dtype=float)
    kcorr = flux_rest / obs

    if len(kcorr) == 1:
        return kcorr[0]

    return kcorr
