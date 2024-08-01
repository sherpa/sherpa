#
#  Copyright (C) 2007, 2015, 2016, 2018 - 2024
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

"""
Routines related to estimating initial parameter values and limits.

.. versionadded:: 4.17.0
   In prior versions these routines were provided by sherpa.utils.

"""

from typing import Any, Optional, Literal, TypedDict, cast, overload

import numpy as np

__all__ = ("_guess_ampl_scale", "get_midpoint", "get_peak",
           "get_valley", "get_fwhm", "guess_fwhm",
           "param_apply_limits", "get_amplitude_position",
           "guess_amplitude", "guess_amplitude_at_ref",
           "guess_amplitude2d", "guess_reference",
           "get_position", "guess_position", "guess_bounds",
           "guess_reference")


# This is liable to be changed at any point (e.g. it could be changed
# to a dataclass).
#
class ValueAndRange(TypedDict):
    """Represent a value and range as a dict."""
    val: float
    min: Optional[float]
    max: Optional[float]


_guess_ampl_scale: float = 1.e+3
"""The scaling applied to a value to create its range.

The minimum and maximum values for a range are calculated by
dividing and multiplying the value by ``_guess_ampl_scale``.
"""


def get_midpoint(a: np.ndarray) -> float:
    """Estimate the middle of the data.

    Parameters
    ----------
    a : array_like
       The data points.

    Returns
    -------
    ans : scalar
       The middle of the data.

    See Also
    --------
    get_peak, guess_bounds

    Notes
    -----
    The estimate is based on the range of the data, and not
    the distribution of the points.
    """

    # return np.abs(a.max() - a.min())/2. + a.min()
    return np.abs(a.max() + a.min()) / 2.0


# TODO: the float return values below are "aspirational"
#
#   It would be nice to link the type of x to the return value but
#   that appears to be tricky to do (thanks to the fact we can store
#   objects as well as numeric values).
#
# TODO: xhi is often unused
#
def get_peak(y: np.ndarray,
             x: np.ndarray,
             xhi: Optional[np.ndarray] = None) -> float:
    """Estimate the peak position of the data.

    Parameters
    ----------
    y, x : array_like
       The data points.
    xhi : None or array_like, optional
       If given then the x array is taken to be the low-edge
       of each bin.

    Returns
    -------
    ans : scalar
       The X position of the peak.

    See Also
    --------
    get_fwhm, get_midpoint, get_valley

    Notes
    -----
    If there are multiple peaks of the same height then
    the location of the first peak is used.
    """

    return x[y.argmax()]


def get_valley(y: np.ndarray,
               x: np.ndarray,
               xhi: Optional[np.ndarray] = None) -> float:
    """Estimate the position of the minimum of the data.

    Parameters
    ----------
    y, x : array_like
       The data points.
    xhi : None or array_like, optional
       If given then the x array is taken to be the low-edge
       of each bin.

    Returns
    -------
    ans : scalar
       The X position of the minimum.

    See Also
    --------
    get_fwhm, get_peak

    Notes
    -----
    If there are multiple minima with the same value then
    the location of the first minimum is used.
    """
    return x[y.argmin()]


def get_fwhm(y: np.ndarray,
             x: np.ndarray,
             xhi: Optional[np.ndarray] = None) -> float:
    """Estimate the width of the data.

    This is only valid for positive data values (``y``).

    Parameters
    ----------
    y, x : array_like
       The data points. The x array must be in ascending order.
    xhi : None or array_like, optional
       If given then the x array is taken to be the low-edge
       of each bin. This is unused.

    Returns
    -------
    ans : scalar
       An estimate of the full-width half-maximum of the peak. If the
       data is negative, or no edge is found then half the X range is
       returned.

    See Also
    --------
    get_peak, get_valley, guess_fwhm

    Notes
    -----
    If there are multiple peaks of the same height then the first peak
    is used.

    The approach is to find the maximum position and then extend out
    to the first bins which fall below half the height. The difference
    of the two points is used. If only one side falls below the value
    then twice this separation is used. If the half-height is not
    reached then the value is set to be half the width of the
    x array. In all cases the upper-edge of the x arrays is ignored,
    if given.

    """

    # Pick half the width of the X array, purely as a guess.
    # The x array is required to be ordered, so we can just
    # take the first and last points.
    #
    guess_fwhm = (x[-1] - x[0]) / 2

    y_argmax = y.argmax()
    if y[y_argmax] <= 0:
        return guess_fwhm

    half_max_val = y[y_argmax] / 2.0
    x_max = x[y_argmax]

    # Where do the values fall below the half-height? The assumption
    # is that the arrays are not so large that evaluating the whole
    # array, rather than just looping out from the maximum location,
    # is not an expensive operation.
    #
    flags = (y - half_max_val) < 0

    # Find the distances from these points to the
    # maximum location.
    #
    dist = x[flags] - x_max

    # We want the maximum value of the negative distances,
    # and the minimum value of the positive distances.
    # There's no guarantee either exist.
    #
    try:
        ldist = -1 * max(dist[dist < 0])
    except ValueError:
        ldist = None

    try:
        rdist = min(dist[dist > 0])
    except ValueError:
        rdist = None

    # If we have both HWHM values then sum them and use that,
    # otherwise if we have one then double it.
    #
    if ldist is not None and rdist is not None:
        return ldist + rdist

    if ldist is not None:
        return 2 * ldist

    if rdist is not None:
        return 2 * rdist

    # No value, so use the guess.
    #
    return guess_fwhm


def guess_fwhm(y: np.ndarray,
               x: np.ndarray,
               xhi: Optional[np.ndarray] = None,
               scale: float = 1000) -> ValueAndRange:
    """Estimate the value and valid range for the FWHM of the data.

    Parameters
    ----------
    y, x : array_like
       The data points.
    xhi : None or array_like, optional
       If given then the x array is taken to be the low-edge
       of each bin.
    scale : number, optional
       The scaling factor applied to the value to calculate the
       minimum (divide) and maximum (multiply) value for the
       FWHM range.

    Returns
    -------
    ans : dict
       The keys are 'val', 'min', and 'max', which give the
       full-width half-maximum and its range.

    See Also
    --------
    get_fwhm

    Notes
    -----
    If there are multiple peaks of the same height then
    the first peak is used.
    """

    fwhm = get_fwhm(y, x, xhi)
    return {'val': fwhm, 'min': fwhm / scale, 'max': fwhm * scale}


def param_apply_limits(param_limits: ValueAndRange,
                       par,  # sherpa.models.parameter.Parameter
                       limits: bool = True,
                       values: bool = True) -> None:
    """Apply the given limits to a parameter.

    This is primarily used by the ``guess`` routine of a model
    to set one or more of its parameters to a given value or
    range.

    Parameters
    ----------
    param_limits : dict
    par : sherpa.models.parameter.Parameter instance
       If the parameter is frozen then nothing is changed.
    limits : bool, optional
       The parameter limits are not changed when ``values`` is
       ``True`` and ``limits`` is ``False``. In all other cases
       the limits are changed.
    values : bool, optional
       When ``True`` the parameter value is changed and the
       original value is stored (for use by the parameter's
       ``reset`` method).

    Examples
    --------

    Create an initial guess for the ``mdl.fwhm`` parameter,
    changing both the value and the soft limits, based on the
    ``x`` and ``y`` arrays.

    >>> vals = guess_fwhm(y, x)
    >>> param_apply_limits(vals, mdl.fwhm)

    Change the soft limits for the ``xpos`` and ``ypos`` parameters
    of the ``src`` model:

    >>> pos = guess_position(y, x0, x1)
    >>> param_apply_limits(pos[0], src.xpos, limits=True, values=False)
    >>> param_apply_limits(pos[1], src.ypos, limits=True, values=False)

    """
    # only guess thawed parameters!
    if par.frozen:
        return

    if limits and values:
        default_val = par.val
        par.set(param_limits['val'], param_limits['min'], param_limits['max'],
                default_min=par.min, default_max=par.max)

        # set original value as default outside the property interface
        par._default_val = default_val

        # set guessed flag to enable user-defined limit reset
        par._guessed = True

    elif values:
        default_val = par.val
        par.set(param_limits['val'])

        # set original value as default outside the property interface
        par._default_val = default_val

    else:
        par.set(min=param_limits['min'], max=param_limits['max'],
                default_min=par.min, default_max=par.max)

        # set guessed flag to enable user-defined limit reset
        par._guessed = True


# TODO: map the return types to the input type
def get_amplitude_position(arr: np.ndarray,
                           mean: bool = False
                           ) -> tuple[float, float, float, Any]:
    """
    Guess model amplitude and position of array
    returns (xval, xmin, xmax, xpos)

    """

    # TODO: better typing for xpos
    xpos: Any = 0
    xmin: float = 0
    xmax: float = 0
    xval: float = 0

    amax = arr.max()
    amin = arr.min()
    if((amax > 0.0 and amin >= 0.0) or
       (amax > 0.0 and amin < 0.0 and (abs(amin) <= amax))):
        xpos = arr.argmax()
        if mean:
            # TODO: should this be "xpos, "?
            xpos = np.where(arr == amax)

        xmax = amax * _guess_ampl_scale
        xmin = amax / _guess_ampl_scale
        xval = amax

    elif((amax > 0.0 and amin < 0.0 and abs(amin) > amax) or
         (amax == 0.0 and amin < 0.0) or (amax < 0.0)):
        xpos = arr.argmin()
        if mean:
            # TODO: should this be "xpos, "?
            xpos = np.where(arr == amin)

        xmax = amin / _guess_ampl_scale
        xmin = amin * _guess_ampl_scale
        xval = amin
    elif (amax == 0.0 and amin == 0.0):
        xpos = arr.argmax()
        if mean:
            # TODO: should this be "xpos, "?
            xpos = np.where(arr == amax)

        xmax = 100.0 / _guess_ampl_scale
        xmin = 0.0
        xval = 0.0

    return (xval, xmin, xmax, xpos)


def guess_amplitude(y: np.ndarray,
                    x: np.ndarray,
                    xhi: Optional[np.ndarray] = None
                    ) -> ValueAndRange:
    """
    Guess model parameter amplitude (val, min, max)

    """

    val, ymin, ymax, pos = get_amplitude_position(y)

    amin, amax = None, None
    if ymin != 0.0 or ymax != 0.0:
        amin = ymin
        amax = ymax

    if xhi is not None:
        # pyright needs the cast otherwise it thinks binsize is a
        # ndarray, presumably because pos is not well typed.
        #
        binsize = cast(float, np.abs(xhi[pos] - x[pos]))
        if amin is not None:
            amin /= binsize
        if amax is not None:
            amax /= binsize
        val /= binsize

    return {'val': val, 'min': amin, 'max': amax}


def guess_amplitude_at_ref(r: float,
                           y: np.ndarray,
                           x: np.ndarray,
                           xhi: Optional[np.ndarray] = None
                           ) -> ValueAndRange:
    """
    Guess model parameter amplitude (val, min, max)

    only valid for beta1d, bpl1d, powlaw1d

    """

    t = 1.0
    if x[1] > x[0] and r < x[0]:
        t = np.abs(y[0] + y[1]) / 2.0
    elif x[1] > x[0] and r > x[-1]:
        t = np.abs(y[-1] + y[-2]) / 2.0
    elif x[1] < x[0] and r > x[0]:
        t = np.abs(y[0] + y[1]) / 2.0
    elif x[1] < x[0] and r < x[-1]:
        t = np.abs(y[-1] + y[-2]) / 2.0
    else:
        for i in range(len(x) - 1):
            if ((r >= x[i] and r < x[i + 1]) or (r >= x[i + 1] and r < x[i])):
                t = np.abs(y[i] + y[i + 1]) / 2.0
                break

    if t == 0.0:
        totband = 0.0
        dv = 0.0
        i = 1
        for j in range(len(x) - 1):
            dv = x[i] - x[i - 1]
            t += y[i] * dv
            totband += dv
            i += 1
        t /= totband

    return {'val': t,
            'min': t / _guess_ampl_scale, 'max': t * _guess_ampl_scale}


@overload
def guess_amplitude2d(y: np.ndarray,
                      x0lo: np.ndarray,
                      x1lo: np.ndarray,
                      x0hi: Literal[None],
                      x1hi: Literal[None]
                      ) -> ValueAndRange:
    ...

@overload
def guess_amplitude2d(y: np.ndarray,
                      x0lo: np.ndarray,
                      x1lo: np.ndarray,
                      x0hi: np.ndarray,
                      x1hi: np.ndarray
                      ) -> ValueAndRange:
    ...

def guess_amplitude2d(y, x0lo, x1lo, x0hi=None, x1hi=None):
    """
    Guess 2D model parameter amplitude (val, min, max)

    """

    limits = guess_amplitude(y, x0lo)

    # if (x0hi is not None and x1hi is not None):
    #     binsize = np.abs((x0hi[0]-x0lo[0])*(x1hi[0]-x1lo[0]))
    #     if limits['min'] is not None:
    #         limits['min'] /= binsize
    #     if limits['max'] is not None:
    #         limits['max'] /= binsize
    #     limits['val'] /= binsize

    return limits


def guess_reference(pmin: float,
                    pmax: float,
                    x: np.ndarray,
                    xhi: Optional[np.ndarray] = None
                    ) -> ValueAndRange:
    """
    Guess model parameter reference (val, min, max)

    """

    xmin = x.min()
    xmax = x.max()

    if xmin >= 1:
        pmin = 1
    if xmax <= 1:
        pmax = 1

    val = 0.0
    if xmin < 1.0 and xmax > 1.0:
        val = 1.0
    else:
        refval = np.floor((xmin + xmax) / 2.0)
        if refval < pmin or refval > pmax:
            refval = (xmin + xmax) / 2.0
        val = refval

    return {'val': val, 'min': None, 'max': None}


def get_position(y: np.ndarray,
                 x: np.ndarray,
                 xhi: Optional[np.ndarray] = None
                 ) -> ValueAndRange:
    """
    Get 1D model parameter positions pos (val, min, max)

    """
    pos = get_amplitude_position(y, mean=True)
    xpos = pos[3]

    val = cast(float, np.mean(x[xpos]))
    xmin = x.min()
    if xhi is None:
        xmax = x.max()
    else:
        xmax = xhi.max()

    return {'val': val, 'min': xmin, 'max': xmax}


@overload
def guess_position(y: np.ndarray,
                   x0lo: np.ndarray,
                   x1lo: np.ndarray,
                   x0hi: Literal[None],
                   x1hi: Literal[None]
                   ) -> tuple[ValueAndRange, ValueAndRange]:
    ...

@overload
def guess_position(y: np.ndarray,
                   x0lo: np.ndarray,
                   x1lo: np.ndarray,
                   x0hi: np.ndarray,
                   x1hi: np.ndarray
                   ) -> tuple[ValueAndRange, ValueAndRange]:
    ...

def guess_position(y, x0lo, x1lo, x0hi=None, x1hi=None):
    """
    Guess 2D model parameter positions xpos, ypos ({val0, min0, max0},
                                                   {val1, min1, max1})

    """

    # pos = int(y.argmax())
    # return the average of location of brightest pixels
    pos = np.where(y == y.max())

    def mkr(x: np.ndarray, xmax: np.ndarray) -> ValueAndRange:
        return {'val': cast(float, np.mean(x[pos])),
                'min': x.min(),
                'max': xmax.max()
                }

    x0, x1 = x0lo, x1lo
    if x0hi is None or x1hi is None:
        r1max = x0
        r2max = x1
    else:
        r1max = x0hi
        r2max = x1hi

    return (mkr(x0, r1max), mkr(x1, r2max))


@overload
def guess_bounds(x: np.ndarray,
                 xhi: Literal[True]
                 ) -> tuple[ValueAndRange, ValueAndRange]:
    ...

@overload
def guess_bounds(x: np.ndarray,
                 xhi: Literal[False]
                 ) -> ValueAndRange:
    ...

def guess_bounds(x, xhi=True):
    """Guess the bounds of a parameter from the independent axis.

    Parameters
    ----------
    x : array_like
       The axis values.
    xhi : bool, optional
       When ``True``, the return value is two dictionaries,
       with values set to 1/3 and 2/3 along the axis, otherwise
       a single dictionary, with a value set to the mid-point
       is returned.

    Returns
    -------
    ans : dict or (dict, dict)
       When ``xhi`` is True then two dictionaries are returned,
       otherwise one. The keys are 'val', 'min', and 'max'.

    See Also
    --------
    get_midpoint

    """

    xmin = x.min()
    xmax = x.max()
    if not xhi:
        lo = xmin + (xmax - xmin) / 2.0
        out: ValueAndRange = {'val': lo, 'min': xmin, 'max': xmax}
        return out

    lo = xmin + (xmax - xmin) / 3.0
    hi = xmin + 2.0 * (xmax - xmin) / 3.0
    out1: ValueAndRange = {'val': lo, 'min': xmin, 'max': xmax}
    out2: ValueAndRange = {'val': hi, 'min': xmin, 'max': xmax}
    return (out1, out2)


@overload
def guess_radius(x0lo: np.ndarray,
                 x1lo: np.ndarray,
                 x0hi: Literal[None],
                 x1hi: Literal[None]
                 ) -> ValueAndRange:
    ...

@overload
def guess_radius(x0lo: np.ndarray,
                 x1lo: np.ndarray,
                 x0hi: np.ndarray,
                 x1hi: np.ndarray
                 ) -> ValueAndRange:
    ...

def guess_radius(x0lo, x1lo, x0hi=None, x1hi=None):
    """Guess the radius parameter of a 2D model.

    Parameters
    ----------
    x0lo, x1lo : array_like
       The independent axes of the grid, in 1D form.
    x0hi, x1hi : array_like or None, optional
       The upper bounds of each pixel.

    Returns
    -------
    ans : dict
       The keys are 'val', 'min', and 'max', which give the
       radius and its range, in units of the input grid (``x0lo``
       and ``x1lo``).

    Notes
    -----
    Currently only ``x0lo`` is used, and it is assumed to be arranged
    so that this axis varies fastest (that is ``x0lo[1] > x0lo[0]``)
    as well as representing square pixels of the same size.

    The scaling factor for the minimum and maximum values is taken
    from the module's `_guess_ampl_scale` variable.

    """
    # TODO: the following was the original code, but
    #   a) x1 isn't used
    #   b) there's no difference between the two branches
    # So, x0hi/x1hi are currently unused.
    #
    # if x0hi is None and x1hi is None:
    #     x0, x1 = x0lo, x1lo
    # else:
    #     x0, x1 = x0lo, x1lo
    x0 = x0lo

    delta = np.apply_along_axis(np.diff, 0, x0)[0]
    rad = np.abs(10 * delta)

    out: ValueAndRange = {'val': rad,
                          'min': rad / _guess_ampl_scale,
                          'max': rad * _guess_ampl_scale}
    return out
