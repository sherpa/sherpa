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

### Ahelp ingest: 2015-05-05 DJB
def calc_energy_flux( data, src, lo=None, hi=None):
    """Integrate the source model over a pass band.

    Calculate the integral of E * S(E) over a pass band, where E is
    the energy of the bin and S(E) the spectral model evaluated for
    that bin.

    Parameters
    ----------
    data :
       The data object to use.
    src :
       The source expression: this should not include any instrument
       responses.
    lo : number, optional
       The minimum limit of the band. Use `None`, the default,
       to use the low value of the data set.
    hi : number, optional
       The maximum limit of the band, which must be larger than
       `lo`. Use `None`, the default, to use the upper value of
       the data set.

    Returns
    -------
    flux :
       The flux from the source model integrated over the given
       band. For X-Spec style models the units will be erg/cm^2/s. If
       `hi` is `None` but `lo` is set then the flux density is
       returned at that point: erg/cm^2/s/keV or erg/cm^2/s/Angstrom
       depending on the analysis setting.

    See Also
    --------
    calc_data_sum : Sum up the data values over a pass band.
    calc_model_sum : Sum up the fitted model over a pass band.
    calc_source_sum: Sum up the source model over a pass band.
    calc_photon_flux : Integrate the source model over a pass band.

    Notes
    -----
    The units of `lo` and `hi` are determined by the analysis
    setting for the data set (e.g. `data.get_analysis`).

    Any existing filter on the data set is ignored by this function.

    The flux is calculated from the given source model, so if it
    includes an absorbing component then the result will represent
    the absorbed flux. The absorbing component can be removed, or
    set to absorb no photons, to get the un-absorbed flux.

    The units of the answer depend on the model components used in
    the source expression and the axis or axes of the data set.
    It is unlikely to give sensible results for 2D data sets.

    Examples
    --------

    Calculate the energy flux over the ranges 0.5 to 2 and 0.5 to
    7 keV:

    >>> data.set_analysis('energy')
    >>> calc_energy_flux(data, smodel, 0.5, 2)
    5.7224906878061796e-10
    >>> calc_energy_flux(data, smodel, 0.5, 7)
    1.3758131915063825e-09

    Calculate the energy flux density at 0.5 keV:

    >>> calc_energy_flux(data, smodel, 0.5)
    5.2573786652855304e-10

    """
    return _flux(data, lo, hi, src, eflux=True)

### Ahelp ingest: 2015-05-05 DJB
def calc_photon_flux( data, src, lo=None, hi=None):
    """Integrate the source model over a pass band.

    Calculate the integral of S(E) over a pass band, where S(E) is the
    spectral model evaluated for each bin.

    Parameters
    ----------
    data :
       The data object to use.
    src :
       The source expression: this should not include any instrument
       responses.
    lo : number, optional
       The minimum limit of the band. Use `None`, the default,
       to use the low value of the data set.
    hi : number, optional
       The maximum limit of the band, which must be larger than
       `lo`. Use `None`, the default, to use the upper value of
       the data set.

    Returns
    -------
    flux :
       The flux from the source model integrated over the given
       band. For X-Spec style models the units will be
       photon/cm^2/s. If `hi` is `None` but `lo` is set then the flux
       density is returned at that point: photon/cm^2/s/keV or
       photon/cm^2/s/Angstrom depending on the analysis setting.

    See Also
    --------
    calc_data_sum : Sum up the data values over a pass band.
    calc_model_sum : Sum up the fitted model over a pass band.
    calc_energy_flux : Integrate the source model over a pass band.
    calc_source_sum: Sum up the source model over a pass band.

    Notes
    -----
    The units of `lo` and `hi` are determined by the analysis
    setting for the data set (e.g. `data.get_analysis`).

    Any existing filter on the data set is ignored by this function.

    The flux is calculated from the given source model, so if it
    includes an absorbing component then the result will represent
    the absorbed flux. The absorbing component can be removed, or
    set to absorb no photons, to get the un-absorbed flux.

    The units of the answer depend on the model components used in
    the source expression and the axis or axes of the data set.
    It is unlikely to give sensible results for 2D data sets.

    Examples
    --------

    Calculate the photon flux over the ranges 0.5 to 2 and 0.5 to
    7 keV, and compared to the energy fluxes for the same bands:

    >>> data.set_analysis('energy')
    >>> calc_photon_flux(data, smodel, 0.5, 2)
    0.35190275
    >>> calc_photon_flux(data, smodel, 0.5, 7)
    0.49050927
    >>> calc_energy_flux(data, smodel, 0.5, 2)
    5.7224906878061796e-10
    >>> calc_energy_flux(data, smodel, 0.5, 7)
    1.3758131915063825e-09

    Calculate the photon flux density at 0.5 keV:

    >>> calc_photon_flux(data, smodel, 0.5)
    0.64978176

    """
    return _flux(data, lo, hi, src)

### Ahelp ingest: 2015-05-05 DJB
### DOC-TODO: compare to calc_photon_flux ?
def calc_source_sum( data, src, lo=None, hi=None):
    """Sum up the source model over a pass band.

    Sum up S(E) over a pass band, where S(E) is the spectral model
    evaluated for each bin.

    Parameters
    ----------
    data :
       The data object to use.
    src :
       The source expression.
    lo : number, optional
       The minimum limit of the band. Use `None`, the default, to use
       the low value of the data set.
    hi : number, optional
       The maximum limit of the band, which must be larger than
       `lo`. Use `None`, the default, to use the upper value of the
       data set.

    Returns
    -------
    signal : number
       The source model summed up over the given band. This does
       *not* include the bin width when using histogram-style
       ('integrated' data spaces), such as used with X-Spec
       emission - also known as additive - models.

    See Also
    --------
    calc_data_sum : Sum up the observed counts over a pass band.
    calc_model_sum : Sum up the fitted model over a pass band.
    calc_energy_flux : Integrate the source model over a pass band.
    calc_photon_flux : Integrate the source model over a pass band.

    Notes
    -----
    The units of `lo` and `hi` are determined by the analysis
    setting for the data set (e.g. `data.get_analysis`).

    Any existing filter on the data set - e.g. as created by
    `ignore` or `notice` - is ignored by this function.

    The units of the answer depend on the model components used in
    the source expression and the axis or axes of the data set.
    It is unlikely to give sensible results for 2D data sets.

    Examples
    --------

    Sum up the model over the data range 0.5 to 2:

    >>> calc_source_sum(data, smodel, 0.5, 2)
    139.12819041922018

    """
    return _flux(data, lo, hi, src, srcflux=True)

#def calc_source_sum2d( data, src, reg=None):
#    return _counts2d(data, reg, data.eval_model_to_fit, src)

### Ahelp ingest: 2015-05-05 DJB
def calc_data_sum(data, lo=None, hi=None):
    """Sum up the data values over a pass band.

    Parameters
    ----------
    data :
       The data object to use.
    lo : number, optional
       The minimum limit of the band. Use `None`, the default, to use
       the low value of the data set.
    hi : number, optional
       The maximum limit of the band, which must be larger than
       `lo`. Use `None`, the default, to use the upper value of the
       data set.

    Returns
    -------
    dsum : number
       The sum of the data values that lie within the given limits.
       If `hi` is `None` but `lo` is set then the data value of the
       bin containing the `lo` value are returned.  If a background
       estimate has been subtracted from the data set then the
       calculation will use the background-subtracted values.

    See Also
    --------
    calc_data_sum2d : Sum up the data values of a 2D data set.
    calc_model_sum : Sum up the fitted model over a pass band.
    calc_energy_flux : Integrate the source model over a pass band.
    calc_photon_flux : Integrate the source model over a pass band.
    calc_source_sum: Sum up the source model over a pass band.
    set_model : Set the source model expression for a data set.

    Notes
    -----
    The units of `lo` and `hi` are determined by the analysis
    setting for the data set (e.g. `data.get_analysis`).

    Any existing filter on the data set - e.g. as created by
    `ignore` or `notice` - is ignored by this function.

    If a grouping scheme has been applied to the data set that it will
    be used. This can change the results, since the first and last
    bins of the selected range may extend outside the requested range.

    Examples
    --------

    Calculate the number of counts over the ranges 0.5 to 2 and 0.5 to
    7 keV, first using the observed signal and then, for the 0.5 to 2
    keV band - the background-subtraced estimate:

    >>> set_analysis('energy')
    >>> calc_data_sum(data, 0.5, 2)
    745.0
    >>> calc_data_sum(data, 0.5, 7)
    60.0
    >>> data.subtract()
    >>> calc_data_sum(data, 0.5, 2)
    730.9179738207356

    """
    return _counts( data, lo, hi, data.apply_filter, data.get_dep() )

def calc_data_sum2d(data, reg=None):
    """Sum up the data values of a 2D data set.
    """
    return _counts2d(data, reg, data.apply_filter, data.get_dep() )

### Ahelp ingest: 2015-05-05 DJB
### DOC-TODO: How does this differe to calc_source_sum, since
###           model here includes the instrumental responses
###           so do we need both? Probably to do with how grids
###           are evaluated.
def calc_model_sum(data, model, lo=None, hi=None):
    """Sum up the fitted model over a pass band.

    Sum up S(E) over a pass band, where S(E) is the spectral model
    evaluated for each bin.

    Parameters
    ----------
    data :
       The data object to use.
    src :
       The source expression, which should include any instrumental
       responses.
    lo : number, optional
       The minimum limit of the band. Use `None`, the default, to use
       the low value of the data set.
    hi : number, optional
       The maximum limit of the band, which must be larger than
       `lo`. Use `None`, the default, to use the upper value of the
       data set.

    Returns
    -------
    signal : number
       The source model summed up over the given band. This does
       *not* include the bin width when using histogram-style
       ('integrated' data spaces), such as used with X-Spec
       emission - also known as additive - models.

    See Also
    --------
    calc_data_sum : Sum up the observed counts over a pass band.
    calc_energy_flux : Integrate the source model over a pass band.
    calc_photon_flux : Integrate the source model over a pass band.
    calc_source_sum: Sum up the source model over a pass band.

    Notes
    -----
    The units of `lo` and `hi` are determined by the analysis
    setting for the data set (e.g. `data.get_analysis`).

    Any existing filter on the data set - e.g. as created by
    `ignore` or `notice` - is ignored by this function.

    The units of the answer depend on the model components used in
    the source expression and the axis or axes of the data set.
    It is unlikely to give sensible results for 2D data sets.

    """
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


### Ahelp ingest: 2015-05-04 DJB
def calc_kcorr(data, model, z, obslo, obshi, restlo=None, resthi=None):
    """Calculate the K correction for a model.

    The K correction ([1]_, [2]_, [3]_, [4]_) is the numeric
    factor applied to measured energy fluxes to convert values in
    an observed energy band to that they are in a rest-frame
    energy band (that is, correct for the change in spectral shape
    between the rest-frame and observed-frame bands). This is
    often used when converting a flux into a luminosity.

    Parameters
    ----------
    data :
       The data object to use.
    model :
       The source expression: this should not include any instrument
       responses.
    z : number or array, >= 0
       The redshift, or redshifts, of the source.
    obslo : number
       The minimum energy of the observed band.
    obshi : number
       The maximum energy of the observed band, which must
       be larger than `obslo`.
    restlo : number or `None`
       The minimum energy of the rest-frame band. If `None` then
       use `obslo`.
    restlo : number or `None`
       The maximum energy of the rest-frame band. It must be
       larger than `restlo`. If `None` then use `obshi`.

    Returns
    -------
    kz : number or array of numbers

    Notes
    -----
    This is only defined when the analysis is in 'energy' units.

    If the model contains a redshift parameter then it should
    be set to `0`, rather than the source redshift.

    If the source model is at zero redshift, the observed energy
    band is olo to ohi, and the rest frame band is rlo to rhi
    (which need not match the observed band), then the K
    correction at a redshift z can be calculated as::

      frest = calc_energy_flux(data, model, rlo, rhi)
      fobs  = calc_energy_flux(data, model, olo*(1+z), ohi*(1+z))
      kz    = frest / fobs

    The energy ranges used - rlo to rhi and olo*(1+z) to ohi*(1+z)
    - should be fully covered by the data grid, otherwise the flux
    calculation will be truncated at the grid boundaries, leading
    to incorrect results.

    References
    ----------

    .. [1] "The K correction", Hogg, D.W., et al.
           http://arxiv.org/abs/astro-ph/0210394

    .. [2] Appendix B of Jones et al. 1998, ApJ, vol 495,
           p. 100-114.
           http://adsabs.harvard.edu/abs/1998ApJ...495..100J

    .. [3] "K and evolutionary corrections from UV to IR",
           Poggianti, B.M., A&AS, 1997, vol 122, p. 399-407.
           http://adsabs.harvard.edu/abs/1997A%26AS..122..399P

    .. [4] "Galactic evolution and cosmology - Probing the
           cosmological deceleration parameter", Yoshii, Y. &
           Takahara, F., ApJ, 1988, vol 326, p. 1-18.
           http://adsabs.harvard.edu/abs/1988ApJ...326....1Y

    """

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
