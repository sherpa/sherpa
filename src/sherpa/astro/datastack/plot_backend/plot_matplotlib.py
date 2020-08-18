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
Plotting routines for the data stack module provided by matplotlib.

"""

import matplotlib.pyplot as plt

name = "pylab_backend"


def initialize_backend():
    """Ensure that the plotting backend is initialized.
    """
    pass


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
    plt.figure(ids.index(dataset['id']) + 1)


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
    plt.figure(ids.index(dataset['id']) + 1)


def save_plot(*args, **kwargs):
    """Save the current plot."""
    plt.savefig(*args, **kwargs)


# How is this different from the _print_window/savefig methods
# of the DataStack class?
plot_savefig = plt.savefig
plot_xlabel = plt.xlabel
plot_ylabel = plt.ylabel
plot_title = plt.title
plot_xlim = plt.xlim
plot_ylim = plt.ylim
plot_set_xscale = plt.xscale
plot_set_yscale = plt.yscale
