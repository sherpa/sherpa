#
#  Copyright (C) 2008, 2016, 2020, 2021, 2022
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

import logging

import numpy

from sherpa.utils import get_position, filter_bins
from sherpa.utils.err import IOErr, DataErr

from ._utils import arf_fold, do_group, expand_grouped_mask, \
    filter_resp, is_in, resp_init, rmf_fold, shrink_effarea
from ._pileup import apply_pileup
from sherpa.astro import hc, charge_e


__all__ = ['arf_fold', 'rmf_fold', 'do_group', 'apply_pileup',
           'eqwidth', 'calc_photon_flux', 'calc_energy_flux',
           'calc_data_sum', 'calc_model_sum', 'shrink_effarea',
           'calc_data_sum2d', 'calc_model_sum2d', 'filter_resp',
           'calc_source_sum', 'compile_energy_grid',
           'calc_kcorr',
           'expand_grouped_mask', 'resp_init', 'is_in',
           'get_xspec_position']


warning = logging.getLogger(__name__).warning


def reshape_2d_arrays(x0, x1):
    new_x0, new_x1 = numpy.meshgrid(x0, x1)
    return new_x0.ravel(), new_x1.ravel()


def get_xspec_position(y, x, xhi=None):
    if xhi is not None:
        if x[0] > x[-1] and xhi[0] > xhi[-1]:
            lo = hc / xhi
            hi = hc / x
            x, xhi = lo, hi
    else:
        if x[0] > x[-1]:
            x = hc / x
    return get_position(y, x, xhi)


def compile_energy_grid(arglist):
    '''Combine several grids (energy, channel, wavelength) into one

    This function combines several grids into one, such that the model
    does not have to be evaluated several times, just because the same energy
    point occurs in every grid.

    This function works under the assumption that evaluating the model
    is expensive and it is is worthwhile to do some bookkeeping to be
    able to derive the model for all intervals needed with the minimum
    number of model points.

    Parameters
    ----------
    arglist : list of tuples
        The list contains the input energy grids. Each tuple is a pair
        of arrays specifying the lower and upper limits for the bins in
        each input grid.

    Returns
    -------
    out : list
        The elements of the list are:
        - elo : array of lower bin values
        - ehi : array of upper bin values
        - htable : list of tuples, in the same format as the input.
             The entries in ``htable`` are indices into ``elo`` and ``ehi``
             that return the original input arrays.

    Examples
    --------

    >>> compile_energy_grid([([1,3,5], [3, 5, 7]), ([0, 1, 2], [1, 2, 3])])
    [array([0, 1, 2, 3, 5]),
     array([1, 2, 3, 5, 7]),
     [(array([1, 3, 4], dtype=int32), array([2, 3, 4], dtype=int32)),
      (array([0, 1, 2], dtype=int32), array([0, 1, 2], dtype=int32))]]

    '''
    elo = numpy.unique(numpy.concatenate([indep[0] for indep in arglist]))
    ehi = numpy.unique(numpy.concatenate([indep[1] for indep in arglist]))

    in_elo = numpy.setdiff1d(elo, ehi)
    in_ehi = numpy.setdiff1d(ehi, elo)
    if len(in_elo) > 1:
        ehi = numpy.concatenate((ehi, in_elo[1:]))
        ehi.sort()
    if len(in_ehi) > 1:
        elo = numpy.concatenate((elo, in_ehi[:- 1]))
        elo.sort()

    # determine index intervals using binary search in large src model
    # for n_th ARF energy bounds and populate a table to use later for
    # integration over said intervals
    htable = [(numpy.searchsorted(elo, arg[0]).astype(numpy.int32),
               numpy.searchsorted(ehi, arg[1]).astype(numpy.int32))
              for arg in arglist]

    return [elo, ehi, htable]


def bounds_check(lo, hi):
    """Ensure that the limits of a filter make sense.

    Note that if only one of them is given we always return
    that value first.
    """

    # TODO: this error message is not ideal, as it assumes the
    # units are energy when the values could be wavelengths,
    # channels, or for a non PHA dataset. See issue #1443
    #
    if lo is not None and hi is not None and lo > hi:
        raise IOErr('boundscheck', lo, hi)

    if lo is not None and hi is None:
        hi = lo

    if hi is not None and lo is None:
        lo = hi

    return (lo, hi)


def range_overlap_1dint(axislist, lo, hi):
    """Return 'overlap' fraction of a 1D int. grid for a range.

    Parameters
    ----------
    axislist : sequence of array, array
        The low and high edges of the bins (in that order). The
        two arrays must have the same size, in either
        ascending or descending order (the same for both), and
        have no gaps. The bins are inclusive for the low edge
        and exclusive for the upper edge.
    lo, hi : None or number
        The range bounds; either both are None or both are set
        with hi >= lo.

    Returns
    -------
    scale : None or numpy array
        An array indicating how much overlap there is between
        the lo-hi grid and each bin (on a scale of 0 to 1).
        There is one wrinkle: if lo = hi then the scale is
        set to 1 for that bin. If there is no match then None
        is returned.

    Notes
    -----
    The bins are assumed to be inclusive for the lower edge and
    exclusive for the higher edge. This means that if a flux
    density is requested (lo=hi) and it falls on an edge then the
    bin returned is the one starting at the requested value.

    The consequences are that if the grid starts at gmin and ends
    at gmax (low edge of first bin and high edge of last bin) then
    lo=hi=gmin will return [1, 0, ...] and lo=hi=gmax returns None.

    As a special case, when lo < hi (so a range is requested) and
    hi=gmin, None is returned (rather than an array of 0's).

    Comparison to the bin edges is done using the standard NumPy
    comparison orders (since it is done by the numpy.digitize
    routine). Given that RMF response grids (which are what this
    is going to be used with) have "interesting" values, this
    can cause surprising issues - for example, the test 3c273.pi
    dataset has an RMF that starts at 0.10000000149011612 rather
    than 0.1 keV, so a flux density calculated at 0.1 keV will
    return 0. An alternative aprproach is to compare to bin edges
    using sao_fcmp, but that would complicate the code.

    """

    # TODO: is compile_energy_grid relevant to this at all?

    # A number of asserts are used to ensure invariants are actually
    # true. They could be removed once this has shown to be
    # working.
    #
    axislo = axislist[0]
    axishi = axislist[1]
    nbins = axislo.size
    assert nbins == axishi.size
    assert axishi[0] > axislo[0]

    if lo is None and hi is None:
        return numpy.ones(axislo.size)

    assert lo is not None
    assert hi is not None
    density = lo == hi

    # To simplify the following, we require that the grid be
    # in ascending order.
    #
    ascending = axislo[1] > axislo[0]
    if ascending:
        edges = numpy.append(axislo, axishi[-1])
    else:
        edges = numpy.append(axishi[0], axislo)[::-1]

    axmin = edges[0]
    axmax = edges[-1]
    assert axmin < axmax, (axmin, axmax)

    # Some of these special cases could be checked before
    # creating edges, but leave here for now.
    #
    if lo >= axmax or hi < axmin:
        return None

    if lo <= axmin and hi >= axmax:
        return numpy.ones(nbins)

    # special case handling of hi == axmin but lo != hi
    #
    if hi == axmin and lo < hi:
        return None

    # At this point there should be no cases where lo=hi
    # and lo < axmin or >= axmax. This means we can replace
    # lo and hi by axmin and axmax (if they exceed the limits),
    # as that makes code below easier.
    #
    if density:
        assert lo >= axmin and lo < axmax

    lo = max(lo, axmin)
    hi = min(hi, axmax)

    # Find the bins corresponding to the start and end
    # points. See the digitize documentation for the
    # meaning of the return values.
    #
    bins = numpy.digitize([lo, hi], edges, right=False)
    blo = bins[0]
    bhi = bins[1]

    # since lo has been set to a minimum of the lower edge,
    # blo should always be >= 1. The upper limit can be
    # > nbins when it is equal to the upper limit.
    #
    assert blo > 0, blo
    # assert bhi <= nbins, (bhi, nbins)

    scale = numpy.zeros(nbins)

    ilo = blo - 1
    if density:
        assert blo == bhi
        scale[ilo] = 1
        return scale if ascending else scale[::-1]

    if blo == bhi:
        # a single bin
        scale[ilo] = (hi - lo) / (edges[blo] - edges[ilo])
        return scale if ascending else scale[::-1]

    # Fully included
    ihi = bhi - 1
    scale[blo:ihi] = 1.0

    # Low edge (may be fully included)
    assert lo >= edges[ilo], (lo, edges[ilo])
    assert lo <= edges[blo], (lo, edges[blo])
    scale[ilo] = (edges[blo] - lo) / (edges[blo] - edges[ilo])

    # High edge (may be fully included)
    if bhi <= nbins:
        assert hi >= edges[ihi], (hi, edges[ihi])
        assert hi <= edges[bhi], (hi, edges[bhi])
        scale[ihi] = (hi - edges[ihi]) / (edges[bhi] - edges[ihi])

    return scale if ascending else scale[::-1]


def _flux(data, lo, hi, src, eflux=False, srcflux=False):
    lo, hi = bounds_check(lo, hi)

    try:
        method = data._get_indep
    except AttributeError:
        method = data.get_indep

    axislist = method(filter=False)
    dim = numpy.asarray(axislist).squeeze().ndim
    if dim > 2:
        raise IOErr('>axes', "2")

    # assume this should not happen, so we do not have to worry
    # about a nice error message
    assert dim > 0

    # To make things simpler, evaluate on the full grid
    y = src(*axislist)

    if srcflux and dim == 2:
        y /= numpy.asarray(axislist[1] - axislist[0])

    if eflux:
        # for energy flux, the sum of grid below must be in keV.
        #
        energ = []
        convert = hasattr(data, 'units') and data.units == 'wavelength'

        for axis in axislist:
            grid = axis
            if convert:
                grid = hc / grid
            energ.append(grid)

        if dim == 2:
            ecorr = 0.5 * (energ[0] + energ[1])
        else:
            # why multiply by 0.5?
            ecorr = 0.5 * energ[0]

        y *= ecorr

    # What bins do we use for the calculation? Linear interpolation
    # is used for bin edges (for integrated data sets)
    #
    if dim == 1:
        mask = filter_bins((lo,), (hi,), (axislist[0],))
        assert mask is not None

        # no bin found
        if numpy.all(~mask):
            return 0.0

        # convert boolean to numbers
        scale = 1.0 * mask

    else:
        scale = range_overlap_1dint(axislist, lo, hi)
        if scale is None:
            return 0.0

        assert scale.max() > 0

    # Originally a flux density was calculated if both lo and hi
    # fell in the same bin, but this has been changed so that
    # we only calculate a density if the lo and hi values are the
    # same (which is set by bounds_check when a density is requested).
    #
    if lo is not None and dim == 2 and lo == hi:
        assert scale.sum() == 1, 'programmer error: sum={}'.format(scale.sum())
        y /= numpy.abs(axislist[1] - axislist[0])

    flux = (scale * y).sum()
    if eflux:
        flux *= charge_e

    return flux


def _counts(data, lo, hi, func, *args):
    """Sum up the data for the range lo to hi.

    The routine starts by saving the current filter, which is a
    combination of the mask field and - if set, which it is only for
    PHA data with ignore_bad called - the quality_filter field.

    The existing noticed range of the data set is ignored,
    and instead the range lo-hi is used (the meaning is slightly
    different to the notice method when one of the two fields
    is None but the other is set; see bounds_check).

    At this point the func method is called, with the remaining
    arguments. It is assumed that this is a method of the
    data object, but there is no requirement. The return value
    is summed up and this value will be returned by the routine.

    It is then time to restore the original filter, which is done in a
    finally block so that any errors calling func or the notice calls
    do not change the data object.

    The use of the name counts is a bit of a mis-nomer as this could
    be used on non-PHA data, but the user only sees this via
    calc_data_sum/calc_model_sum.

    """
    lo, hi = bounds_check(lo, hi)
    old_mask = data.mask
    old_quality_filter = getattr(data, 'quality_filter', None)
    try:
        data.notice()  # clear filter
        data.notice(lo, hi)
        counts = func(*args).sum()
    finally:
        # restore the old filter
        data.mask = old_mask
        if old_quality_filter is not None:
            data.quality_filter = old_quality_filter

    return counts


def _counts2d(data, reg, func, *args):
    """Sum up the data for the given region.

    This is the same as _counts but for a data object that supports
    the notice2d method (so currently DataIMG and DataIMGInt). It has
    the same structure as _counts except that the fields used to save
    and restore the state of the data object match the DataIMG
    interface (i.e. mask and optional _region rather than mask and an
    optional quality_filter)

    The use of the name counts is a bit of a mis-nomer as this could
    be used on non-PHA data, but the user only sees this via
    calc_data_sum2d/calc_model_sum2d.

    """
    old_mask = data.mask
    old_region = getattr(data, '_region', None)
    try:
        data.notice2d()  # save and clear filter
        data.notice2d(reg)
        counts = func(*args).sum()
    finally:
        # restore the old filter
        data.mask = old_mask
        if old_region is not None:
            data._region = old_region

    return counts


def calc_energy_flux(data, src, lo=None, hi=None):
    """Integrate the source model over a pass band.

    Calculate the integral of E * S(E) over a pass band, where E is
    the energy of the bin and S(E) the spectral model evaluated for
    that bin.

    Parameters
    ----------
    data
       The data object to use.
    src
       The source expression: this should not include any instrument
       responses.
    lo, hi : number, optional
       If both are None or both are set then calculate the flux
       over the given band. If only one is set then calculate
       the flux density at that point. The units for `lo` and `hi`
       are given by the current analysis setting of the `data`
       parameter.

    Returns
    -------
    flux
       The flux or flux density of the source model. For X-Spec
       models the flux units will be erg/cm^2/s and the flux
       density is either erg/cm^2/s/keV or erg/cm^2/s/Angstrom,
       depending on the analysis setting.

    See Also
    --------
    calc_data_sum : Sum up the data values over a pass band.
    calc_model_sum : Sum up the fitted model over a pass band.
    calc_source_sum : Sum up the source model over a pass band.
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


def calc_photon_flux(data, src, lo=None, hi=None):
    """Integrate the source model over a pass band.

    Calculate the integral of S(E) over a pass band, where S(E) is the
    spectral model evaluated for each bin.

    Parameters
    ----------
    data
       The data object to use.
    src
       The source expression: this should not include any instrument
       responses.
    lo, hi : number, optional
       If both are None or both are set then calculate the flux
       over the given band. If only one is set then calculate
       the flux density at that point. The units for `lo` and `hi`
       are given by the current analysis setting of the `data`
       parameter.

    Returns
    -------
    flux
       The flux or flux density of the source model. For X-Spec
       models the flux units will be photon/cm^2/s and the flux
       density is either photon/cm^2/s/keV or
       photon/cm^2/s/Angstrom, depending on the analysis setting.

    See Also
    --------
    calc_data_sum : Sum up the data values over a pass band.
    calc_model_sum : Sum up the fitted model over a pass band.
    calc_energy_flux : Integrate the source model over a pass band.
    calc_source_sum : Sum up the source model over a pass band.

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


# ## DOC-TODO: compare to calc_photon_flux ?
def calc_source_sum(data, src, lo=None, hi=None):
    """Sum up the source model over a pass band.

    Sum up S(E) over a pass band, where S(E) is the spectral model
    evaluated for each bin.

    Parameters
    ----------
    data
       The data object to use.
    src
       The source expression. This must not include the instrumental
       responses.
    lo, hi : number, optional
       If both are None or both are set then sum up over the given
       band. If only one is set then use the model value in the
       selected bin. The units for `lo` and `hi` are given by the
       current analysis setting of the `data` object.

    Returns
    -------
    signal : number
       The source model summed up over the given band or for
       a single bin.

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


# def calc_source_sum2d( data, src, reg=None):
#     return _counts2d(data, reg, data.eval_model_to_fit, src)

def calc_data_sum(data, lo=None, hi=None):
    """Sum up the data values over a pass band.

    Parameters
    ----------
    data
       The data object to use.
    lo, hi : number, optional
       If both are None or both are set then use the full dataset.
       If only one is set then use the data count for that bin.
       The units for `lo` and `hi` are given by the current analysis
       setting of the `data` parameter.

    Returns
    -------
    dsum : number
       The sum of the data values that lie within the given limits.
       If a background estimate has been subtracted from the data set
       then the calculation will use the background-subtracted values.

    See Also
    --------
    calc_data_sum2d : Sum up the data values of a 2D data set.
    calc_model_sum : Sum up the fitted model over a pass band.
    calc_energy_flux : Integrate the source model over a pass band.
    calc_photon_flux : Integrate the source model over a pass band.
    calc_source_sum : Sum up the source model over a pass band.
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
    return _counts(data, lo, hi, data.apply_filter, data.get_dep())


def calc_data_sum2d(data, reg=None):
    """Sum up the data values of a 2D data set.

    Parameters
    ----------
    data : sherpa.astro.data.DataIMG instance
       The data object to use.
    reg : str, optional
       The spatial filter to use. The default, ``None``, is to use the
       whole data set.

    Returns
    -------
    dsum : number
       The sum of the data values that lie within the given region.

    See Also
    --------
    calc_data_sum : Sum up the data values of a data set.
    calc_model_sum2d : Sum up the fitted model for a 2D data set.

    Notes
    -----
    The coordinate system of the region filter is determined by the
    coordinate setting for the data set (e.g. `data.coord`).

    Any existing filter on the data set - e.g. as created by
    `ignore2d` or `notice2d` - is ignored by this function.

    """
    return _counts2d(data, reg, data.apply_filter, data.get_dep())


# ## DOC-TODO: better comparison of calc_source_sum and calc_model_sum
# ##           needed (e.g. integration or results in PHA case?)
def calc_model_sum(data, model, lo=None, hi=None):
    """Sum up the fitted model over a pass band.

    Sum up S(E) over a pass band, where S(E) is the spectral model
    evaluated for each bin.

    Parameters
    ----------
    data
       The data object to use.
    model
       The source expression, which should include the instrumental
       responses.
    lo, hi : number, optional
       If both are None or both are set then sum up over the given
       band. If only one is set then use the model value in the
       selected bin. The units for `lo` and `hi` are given by the
       current analysis setting of the `data` object.

    Returns
    -------
    signal : number
       The source model summed up over the given band or for
       a single bin.

    See Also
    --------
    calc_data_sum : Sum up the observed counts over a pass band.
    calc_energy_flux : Integrate the source model over a pass band.
    calc_photon_flux : Integrate the source model over a pass band.
    calc_source_sum : Sum up the source model over a pass band.

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


# ## DOC-TODO: clean up whether the calc_model_* versions should or
# ##           should not contain the instrument response/PSF components.
# ##           Note: there is no calc_source_sum2d in this module, so
# ##           this needs looking at to see if the text is correct
def calc_model_sum2d(data, model, reg=None):
    """Sum up the fitted model for a 2D data set.

    Parameters
    ----------
    data : sherpa.astro.data.DataIMG instance
       The data object to use.
    model
       The source expression, which should not include the PSF model.
    reg : str, optional
       The spatial filter to use. The default, ``None``, is to use the
       whole data set.

    Returns
    -------
    dsum : number
       The sum of the model values, including any PSF convolution,
       that lie within the given region.

    See Also
    --------
    calc_data_sum2d : Sum up the data values of a 2D data set.
    calc_model_sum : Sum up the fitted model over a pass band.

    Notes
    -----
    The coordinate system of the region filter is determined by the
    coordinate setting for the data set (e.g. `data.coord`).

    Any existing filter on the data set - e.g. as created by
    `ignore2d` or `notice2d` - is ignored by this function.

    """
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
        my = model(xlo, xhi)
        cy = combo(xlo, xhi)
        num = len(xlo)
    else:
        my = data.eval_model_to_fit(model)
        cy = data.eval_model_to_fit(combo)
        xlo = data.get_indep(filter=True)[0]
        num = len(xlo)

    # TODO: should this follow _flux and handle the case when
    #       we have xlo, xhi differently?
    #
    mask = filter_bins((lo,), (hi,), (xlo,))
    if mask is not None:
        my = my[mask]
        cy = cy[mask]
        xlo = xlo[mask]
        num = len(xlo)

    for ebin, val in enumerate(xlo):
        if ebin < (num - 1):
            eave = numpy.abs(xlo[ebin + 1] - xlo[ebin])
        else:
            eave = numpy.abs(xlo[ebin - 1] - xlo[ebin])
        if my[ebin] != 0.0:
            eqw += eave * (cy[ebin] - my[ebin]) / my[ebin]

    return eqw


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
    data
       The data object to use.
    model
       The source expression: this should not include any instrument
       responses.
    z : number or array, >= 0
       The redshift, or redshifts, of the source.
    obslo : number
       The minimum energy of the observed band.
    obshi : number
       The maximum energy of the observed band, which must
       be larger than `obslo`.
    restlo : number or ``None``
       The minimum energy of the rest-frame band. If ``None`` then
       use `obslo`.
    restlo : number or ``None``
       The maximum energy of the rest-frame band. It must be
       larger than `restlo`. If ``None`` then use `obshi`.

    Returns
    -------
    kz : number or array of numbers

    Notes
    -----
    This is only defined when the analysis is in 'energy' units.

    If the model contains a redshift parameter then it should
    be set to ``0``, rather than the source redshift.

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

    if 0 != sum(z[z < 0]):
        raise IOErr('z<=0')

    if obslo <= 0 or restlo <= 0 or obshi <= obslo or resthi <= restlo:
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
        elo, ehi = data.get_indep()

    if elo is None or ehi is None:
        raise DataErr('noenergybins', data.name)

    emin = elo[0]
    emax = ehi[-1]

    if restlo < emin or resthi > emax:
        raise IOErr('energoverlap', emin, emax, 'rest-frame',
                    restlo, resthi, '')

    if obslo * (1.0 + z.min()) < emin:
        raise IOErr('energoverlap', emin, emax, 'observed-frame',
                    restlo, resthi, "at a redshift of %f" % z.min())

    if obshi * (1.0 + z.max()) > emax:
        raise IOErr('energoverlap', emin, emax, 'rest-frame',
                    restlo, resthi, "at a redshift of %f" % z.min())

    zplus1 = z + 1.0
    flux_rest = _flux(data, restlo, resthi, model, eflux=True)
    obs = numpy.asarray([_flux(data, obslo * zz, obshi * zz, model, eflux=True)
                         for zz in zplus1], dtype=float)
    kcorr = flux_rest / obs

    if len(kcorr) == 1:
        return kcorr[0]

    return kcorr
