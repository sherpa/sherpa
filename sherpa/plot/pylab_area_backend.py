#
#  Copyright (C) 2022 - 2024
#  MIT
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
'''This module contains a modified Pylab backend.

The backend is fully functional and can be used in real applications, but
its main purpose is to serve as an example of inheritance in backends.

For that purpose, it is convenient to split it out into its own file, so that
specific line numbers from this file can be directly included in the documentation.
'''
# Specific line numbers are included in docs/plots/backends.rst
# When this file is edited, you need to manually check that the line
# numbers given there are still correct.

from sherpa.plot.pylab_backend import PylabBackend
from sherpa.plot.backends import translate_args

__all__ = ('PylabErrorArea', )


class PylabErrorArea(PylabBackend):
    '''A Matplotlib backend displaying data uncertainties as shaded regions

    This class changes the behavior of the plotting in a way that is not
    possible with just setting parameters alone: For 1D data with y-errors
    the error range is not shown with error bars as in
    `sherpa.plot.backends.PylabBackend`, but instead with shaded regions.

    This class does not display x errors, even if they are given.
    '''

    @property
    def name(self):
        '''An easy-to-read string name for a plotting backend.'''
        return self.__class__.__name__

    @translate_args
    def plot(self, x, y, *,
             yerr=None, xerr=None, title=None,
             xlabel=None, ylabel=None,
             overplot=False, clearwindow=True,
             xerrorbars=False,
             yerrorbars=False,
             xlog=False,
             ylog=False,
             linestyle='solid',
             drawstyle='default',
             color=None,
             ecolor=None,
             marker='None',
             markerfacecolor=None,
             markersize=None,
             alpha=None,
             label=None,
             linewidth=None,
             capsize=None,
             # derives from PylabBackend, which has "barsabove", so plot might be called with
             # that parameter from e.g. histogram. So, we allow that parameter, but ignore it.
             barsabove=False,
             ):

        if x is None or y is None:
            return None

        axes = self.setup_axes(overplot, clearwindow)

        # Set up the axes
        if not overplot:
            self.setup_plot(axes, title, xlabel, ylabel, xlog=xlog, ylog=ylog)

        if yerrorbars and yerr is not None:
            # Are the errors symmetric or not? At the moment this is
            # decided by checking if yerr is a tuple or not.
            #
            if type(yerr) == tuple:
                ylo = y - yerr[0]
                yhi = y + yerr[1]
            else:
                ylo = y - yerr
                yhi = y + yerr

            axes.fill_between(x, ylo, yhi, alpha=0.2, linewidth=0,
                              color=ecolor)

        return axes.plot(x, y,
                         color=color,
                         alpha=alpha,
                         linestyle=linestyle,
                         linewidth=linewidth,
                         drawstyle=drawstyle,
                         marker=marker,
                         markersize=markersize,
                         markerfacecolor=markerfacecolor,
                         label=label)
