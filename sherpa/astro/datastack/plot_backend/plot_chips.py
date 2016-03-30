#
#  Copyright (C) 2010, 2014, 2015  Smithsonian Astrophysical Observatory
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
Plotting routines for the data stack module provided by ChIPS.

"""

import pychips
from sherpa.plot import chips_backend

name = "chips_backend"


def initialize_backend():
    """Ensure that the plotting backend is initialized."""

    # FIXME We need to make sure ChIPS is initialized.
    # As a hack, we just call begin() and end() on the chips backend
    # To make sure the backend is initialized before we start
    # creating windows.
    try:
        chips_backend.begin()
    finally:
        chips_backend.end()


def initialize_plot(dataset, ids):
    """Create the plot window or figure for the given dataset.

    Parameters
    ----------
    dataset : str or int
       The dataset.
    ids : array_like
       The identifier array from the DataStack object.

    See Also
    --------
    select_plot

    """
    try:
        pychips.add_window(['id', dataset['id']])
    except RuntimeError:
        pass
    pychips.current_window(str(dataset['id']))


def select_plot(dataset, ids):
    """Select the plot window or figure for the given dataset.

    The plot for this dataset is assumed to have been created.

    Parameters
    ----------
    dataset : str or int
       The dataset.
    ids : array_like
       The identifier array from the DataStack object.

    See Also
    --------
    initialize_plot

    """
    pychips.set_current_window(dataset['id'])


def save_plot(*args):
    """Save the current plot."""
    pychips.print_window(*args)

# The original code provided matplotlib-only access to these
# functions. Now also provide ChIPS access too.

plot_savefig = pychips.print_window
plot_xlabel = pychips.set_plot_xlabel
plot_ylabel = pychips.set_plot_ylabel
plot_title = pychips.set_plot_title


# It does not quite match the matplotlib API:
#   - no support for the argument being a pair
#   - does not return the new limits
#   - active data range after calling this may not match
#     ChIPS expectations, but okay for here
def plot_xlim(xmin=None, xmax=None):
    """Set the limits of the X axis.
    If no parameters are set, then return the current limits.

    Parameters
    ----------
    xmin
       If set, change the minimum value
    xmax
       If set, change the maximum value
    """
    if xmin is None and xmax is None:
        return pychips.get_plot_range()[:2]

    if xmin is None:
        xmin = pychips.get_plot_range()[0]

    if xmax is None:
        xmax = pychips.get_plot_range()[1]

    pychips.limits(pychips.X_AXIS, xmin=xmin, xmax=xmax)
    return pychips.get_plot_range()[:2]


def plot_ylim(ymin=None, ymax=None):
    """Set the limits of the Y axis.
    If no parameters are set, then return the current limits.

    Parameters
    ----------
    ymin
       If set, change the minimum value
    ymax
       If set, change the maximum value
    """
    return_function = lambda: pychips.get_plot_range()[2:]

    if ymin is None and ymax is None:
        return return_function()

    if ymin is None:
        ymin = pychips.get_plot_range()[2]

    if ymax is None:
        ymax = pychips.get_plot_range()[3]

    pychips.limits(pychips.Y_AXIS, ymin=ymin, ymax=ymax)
    return return_function()


def plot_set_xscale(scale):
    """Change the scaling of the X axis.

    Parameters
    ----------
    scale : 'log', 'lin', 'logarithmic', 'linear'
    """
    lscale = scale.lower()
    if lscale in ['log', 'logarithmic']:
        pychips.log_scale(pychips.X_AXIS)
    elif lscale in ['lin', 'linear']:
        pychips.lin_scale(pychips.X_AXIS)
    else:
        raise ValueError("Unsupported axis scale: {0}".format(scale))


def plot_set_yscale(scale):
    """Change the scaling of the Y axis.

    Parameters
    ----------
    scale : 'log', 'lin', 'logarithmic', 'linear'
    """
    lscale = scale.lower()
    if lscale in ['log', 'logarithmic']:
        pychips.log_scale(pychips.Y_AXIS)
    elif lscale in ['lin', 'linear']:
        pychips.lin_scale(pychips.Y_AXIS)
    else:
        raise ValueError("Unsupported axis scale: {0}".format(scale))
