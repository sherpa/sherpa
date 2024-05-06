#
#  Copyright (C) 2010, 2015, 2017, 2019 - 2024
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
import io
import logging
from collections import ChainMap
from functools import partialmethod

import numpy

from matplotlib import pyplot as plt

from sherpa.utils.err import ArgumentErr, NotImplementedErr
from sherpa.utils import formatting
from sherpa.plot.utils import histogram_line
from sherpa.plot.backends import BasicBackend
from sherpa.plot.backends import kwargs_doc as orig_kwargs_doc
from sherpa.plot.backend_utils import (translate_args,
                                       add_kwargs_to_doc)

__all__ = ['PylabBackend']

logger = logging.getLogger(__name__)

# Most of the kwargs docstrings are still true, but some can be updated
# because matplotlib does support a wider range than just the
# set of Sherpa backend-independent options.
updated_kwarg_docs = {
    'color': ['str or tuple', 'Any matplotlib color'],
    'linecolor': ['str or tuple', 'Any matplotlib color'],
    'marker': ['str',
               '''"None" (as a string, no marker shown), "" (empty string, no marker shown),
or any matplotlib marker, e.g. any of `os+v` or others
(see matplotlib documentation).'''],
    'linestyle': ['str',
                  '''``'noline'``,
``'None'`` (as string, same as ``'noline'``),
``'solid'``, ``'dot'``, ``'dash'``, ``'dashdot'``, ``'-'`` (solid
line), ``':'`` (dotted), ``'--'`` (dashed), ``'-.'`` (dashdot),
``''`` (empty string, no line shown), `None` (default - usually
solid line) or any other matplotlib linestyle.'''],
    'drawstyle': ['str', 'matplolib drawstyle'],
}

kwargs_doc = ChainMap(updated_kwarg_docs, orig_kwargs_doc)


class PylabBackend(BasicBackend):
    """sherpa plotting backend using matplotlib"""

    name = 'pylab'
    '''Alternative string names for plotting backend.

    This differs from the class name, because external code such
    as ciao_contrib scripts might depend on this name being present.
    '''

    translate_dict = {'linestyle': {
        'noline': ' ',
        'solid': '-',
        'dot': ':',
        'dash': '--',
        'dashdot': '-.',
        'None': ' ', },
        'label': {None: '_nolegend_'},
    }
    '''Dict of keyword arguments that need to be translated for this backend.

    The keys in this dict are keyword arguments(e.g. ``'markerfacecolor'``)
    and the values are one of the following:
    - A dict where the keys are the backend-independent values and the values
    are the values expressed for this backend. For values not listed in the
    dict, no translation is done.
    - A callable. The callable translation function is called with any argument
    given and allows the backend arbitrary translations. It should, at the
    very least, accept all backend independent values for this parameter
    without error.

    Example:

       >>> translate_dict = {'markerfacecolor': {'k': (0., 0., 0.)},
       ...                   'alpha': lambda a: 256 * a}

    This translates the color 'k' to tuple of RGB values and alpha values
    to a number between 0 and 256.
    '''

    def colorlist(self, n):
        '''Generate the list of n colors for use in multi-line plots.

        Generally, the color will be ordered in some way and do not repeat or
        do so only after a large number of colors.
        Different backends might generate different lists of colors.

        Parameters
        ----------
        n : int
            Number of colors requested

        Returns
        -------
        colors : list
            list of color specifiers
        '''
        return plt.cm.inferno(numpy.linspace(0, 1, n))

    def clear_window(self):
        """Clear default pyplot figure."""
        plt.clf()

    def __exit__(self, exec_type, value, traceback):
        plt.draw()
        return False

    def setup_axes(self, overplot, clearwindow):
        """Return the axes object, creating it if necessary.

        Parameters
        ----------
        overplot : bool
        clearwindow : bool

        Returns
        -------
        axis
            The matplotlib axes object.
        """

        # Do we need to clear the window?
        if not overplot and clearwindow:
            self.clear_window()

        return plt.gca()

    @add_kwargs_to_doc(kwargs_doc)
    def setup_plot(self, axes, title=None, xlabel=None, ylabel=None,
                   xlog=False, ylog=False):
        """Basic plot setup.

        Parameters
        ----------
        axes
            The plot axes (output of setup_axes).
        {kwargs}

        """
        axes.set_xscale('log' if xlog else 'linear')
        axes.set_yscale('log' if ylog else 'linear')
        axes.set_title(title)
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)

    def find_zorder(self, axes):
        """Try to come up with a good zorder value

        Parameters
        ----------
        axes
            The plot axes

        Returns
        -------
        zorder : float or None
            The estimated zorder value.

        Notes
        -----
        The following is from https://github.com/sherpa/sherpa/issues/662
        which is circa matplotlib 3. The issue is how to ensure that
        a plot with multiple data values (e.g. a fit plot with
        data+error bars and model, and one with multiple fits),
        are drawn "sensibly" so that the user can make out the model.
        For this we really want the zorder of plot items to be determined
        by the order they are added. However, this appears to be complicated
        by the fact that the errorbar command creates multiple objects
        (for the point, error bars, and I believe caps), and not just
        a Line2D object. Once all the plot items are concatenated and
        sorted to ensure a consistent ordering - as discussed in
        https://github.com/matplotlib/matplotlib/issues/1622#issuecomment-142469866 -
        we lose a "nice" order in Sherpa (for our cases it does not
        seem to help that the error bar adds in a zorder Line2D of 2.1
        as compared to 2, which means that the point for the error bar
        is drawn on top of the error bar, but also above other lines
        added to the plot.

        One option would be to evaluate the current plot to find out
        what the minimum zorder is, and explicitly set larger
        values. Note that the default zindex values appear to be 2 and
        2.1 from plot/errorbar. This should work for the case of
        plot something
        overplot something
        but may fall apart as soon as users add their own features
        to the visualization, but that may be acceptable.
        """

        # Is this a sufficient set of objects?
        #
        objects = axes.lines + axes.collections

        # could do a list comprehension, but want to catch
        # AttributeErrors
        zs = []
        for o in objects:
            try:
                zo = o.get_zorder()
            except AttributeError:
                continue

            zs.append(zo)

        if len(zs) > 0:
            # Add a small offset for this set of data
            return max(zs) + 0.1

        return None

    @add_kwargs_to_doc(kwargs_doc)
    @translate_args
    def histo(self, xlo, xhi, y, *,
              yerr=None, title=None,
              xlabel=None, ylabel=None,
              overplot=False, clearwindow=True,
              xerrorbars=False,
              yerrorbars=False,
              ecolor=None,
              capsize=None,
              barsabove=False,
              xlog=False,
              ylog=False,
              linestyle='solid',
              drawstyle='default',
              color=None,
              alpha=None,
              marker='None',
              markerfacecolor=None,
              markersize=None,
              label=None,
              linewidth=None,
              linecolor=None,
              ):
        """Draw histogram data.

        The histogram is drawn as horizontal lines connecting the
        start and end points of each bin, with vertical lines connecting
        consecutive bins. Non-consecutive bins are drawn with a
        (NaN, NaN) between them so no line is drawn connecting them.

        Points are drawn at the middle of the bin, along with any
        error values.

        Note that the linecolor is not used, and is only included
        to support old code that may have set this option (use `color` instead).

        Parameters
        ----------
        x0 : array-like or scalar number
            lower bin boundary values
        x1 : array-like or scalar number
            upper bin boundary values
        y : array-like or scalar number
            y values, same dimension as `x0`.
        {kwargs}
        """
        if linecolor is not None:
            logger.warning(
                "The linecolor attribute ({}) is unused.".format(linecolor))
        x, y2 = histogram_line(xlo, xhi, y)

        # Note: this handles clearing the plot if needed.
        #
        objs = self.plot(x, y2, yerr=None, xerr=None,
                         title=title, xlabel=xlabel, ylabel=ylabel,
                         overplot=overplot, clearwindow=clearwindow,
                         xerrorbars=False, yerrorbars=False,
                         label=label,
                         ecolor=ecolor, capsize=capsize, barsabove=barsabove,
                         xlog=xlog, ylog=ylog,
                         linestyle=linestyle, linewidth=linewidth,
                         drawstyle=drawstyle,
                         color=color, marker=None, alpha=alpha,
                         )

        # Draw points and error bars at the mid-point of the bins.
        # Use the color used for the data plot: should it
        # also be used for marker[face]color and ecolor?
        #
        # NOTE: this ignores any xerr parameter.
        #
        xmid = 0.5 * (xlo + xhi)
        xerr = (xhi - xlo) / 2 if xerrorbars else None
        yerr = yerr if yerrorbars else None

        try:
            color = objs[0].get_color()
        except AttributeError:
            pass

        axes = plt.gca()
        zorder = self.find_zorder(axes)

        # Do not draw a line connecting the points
        #
        # Unlike plot, using errorbar for both cases.
        #

        axes.errorbar(xmid, y, yerr, xerr,
                      color=color,
                      alpha=alpha,
                      linestyle='',
                      marker=marker,
                      markersize=markersize,
                      markerfacecolor=markerfacecolor,
                      ecolor=ecolor,
                      capsize=capsize,
                      barsabove=barsabove,
                      linewidth=linewidth,
                      # We do not want two labels for the same dataset
                      # TODO: Do we want errorbars in the label?
                      # Currently, we gives a line with no marker
                      # label=label,
                      zorder=zorder)
        handles, _ = axes.get_legend_handles_labels()
        if len(handles) > 0:
            axes.legend()

    # Note that this is an internal method that is not wrapped in
    # @ translate_args
    # This is called from methods like plot, histo, etc.
    # so the arguments passed to this method are already translated and we do
    # not want to apply the translation a second time.
    def _set_line(self, line, linecolor=None, linestyle=None, linewidth=None):
        """Apply the line attributes, if set.

        Parameters
        ----------
        line : matplotlib.lines.Line2D
            The line to change
        linecolor, linestyle, linewidth : optional
            The attribute value or None.
        """
        def setf(label, newval):
            if newval is None:
                return

            func = getattr(line, f'set_{label}')
            func(newval)

        setf('color', linecolor)
        setf('linestyle', linestyle)
        setf('linewidth', linewidth)

    # There is no support for alpha in the Plot.vline class
    @add_kwargs_to_doc(kwargs_doc)
    @translate_args
    def vline(self, x, *,
              ymin=0, ymax=1,
              linecolor=None,
              linestyle=None,
              linewidth=None,
              overplot=False, clearwindow=True):
        """Draw a vertical line

        Parameters
        ----------
        x : float
            x position of the vertical line in data units
        {kwargs}
        """
        axes = self.setup_axes(overplot, clearwindow)

        line = axes.axvline(x, ymin, ymax)
        self._set_line(line, linecolor=linecolor, linestyle=linestyle,
                       linewidth=linewidth)

    # There is no support for alpha in the Plot.hline class

    @add_kwargs_to_doc(kwargs_doc)
    @translate_args
    def hline(self, y, *,
              xmin=0, xmax=1,
              linecolor=None,
              linestyle=None,
              linewidth=None,
              overplot=False, clearwindow=True):
        """Draw a horizontal line

        Parameters
        ----------
        y : float
            x position of the vertical line in data units
        {kwargs}
        """
        axes = self.setup_axes(overplot, clearwindow)

        line = axes.axhline(y, xmin, xmax)
        self._set_line(line, linecolor=linecolor, linestyle=linestyle,
                       linewidth=linewidth)

    @add_kwargs_to_doc(kwargs_doc)
    @translate_args
    def plot(self, x, y, *,
             yerr=None, xerr=None, title=None,
             xlabel=None, ylabel=None,
             overplot=False, clearwindow=True,
             xerrorbars=False,
             yerrorbars=False,
             ecolor=None,
             capsize=None,
             barsabove=False,
             xlog=False,
             ylog=False,
             linestyle='solid',
             drawstyle='default',
             color=None,
             marker='None',
             markerfacecolor=None,
             markersize=None,
             alpha=None,
             label=None,
             linewidth=None,
             linecolor=None,
             xaxis=None, ratioline=None,
             ):
        """Draw x, y data.

        This method combines a number of different ways to draw x/y data:
            - a line connecting the points
            - scatter plot of symbols
            - errorbars

        All three of them can be used together (symbols with errorbars connected
        by a line), but it is also possible to use only one or two of them. By
        default, a line is shown (``linestyle='solid'``), but marker and error
        bars are not (``marker='None'`` and ``xerrorbars=False`` as well as
        ``yerrorbars=False``).

        Note that the linecolor is not used, and is only included
        to support old code that may have set this option (use `color` instead).

        Parameters
        ----------
        x : array-like or scalar number
            x values
        y : array-like or scalar number
            y values, same dimension as `x`.
        {kwargs}
        ratioline, xaxis : None
            These parameters are deprecated and not used any longer.
        """
        if xaxis is not None:
            logger.warning('Keyword "xaxis" is deprecated and has no effect. xaxis are always drawn for delta_xxx plots.')

        if ratioline is not None:
            logger.warning('Keyword "ratioline" is deprecated and has no effect. Ratio lines are always drawn for ratio plots.')

        axes = self.setup_axes(overplot, clearwindow)

        if linecolor is not None:
            logger.warning(
                f"The linecolor attribute ({linecolor}) is unused.")

        # Set up the axes
        if not overplot:
            self.setup_plot(axes, title, xlabel, ylabel, xlog=xlog, ylog=ylog)

        # See the discussion in find_zorder
        zorder = None
        if overplot:
            zorder = self.find_zorder(axes)

        # Rely on color-cycling to work for both the "no errorbar" and
        # "errorbar" case.
        #
        # TODO: do we really need this, or can we just use errorbar
        #       for both cases?
        #
        if xerrorbars or yerrorbars:
            if markerfacecolor is None:
                markerfacecolor = color
            xerr = xerr if xerrorbars else None
            yerr = yerr if yerrorbars else None
            obj = axes.errorbar(x, y, yerr, xerr,
                                 label=label,
                                 color=color,
                                 linestyle=linestyle,
                                 linewidth=linewidth,
                                 marker=marker,
                                 markersize=markersize,
                                 markerfacecolor=markerfacecolor,
                                 alpha=alpha,
                                 ecolor=ecolor,
                                 capsize=capsize,
                                 barsabove=barsabove,
                                 zorder=zorder)

        else:
            obj = axes.plot(x, y,
                            color=color,
                            alpha=alpha,
                            linestyle=linestyle,
                            linewidth=linewidth,
                            label=label,
                            drawstyle=drawstyle,
                            marker=marker,
                            markersize=markersize,
                            markerfacecolor=markerfacecolor,
                            zorder=zorder)
        handles, _ = axes.get_legend_handles_labels()
        if len(handles) > 0:
            axes.legend()
        return obj


    @add_kwargs_to_doc(kwargs_doc)
    @translate_args
    def contour(self, x0, x1, y, *,
                levels=None, title=None,
                xlabel=None, ylabel=None,
                overcontour=False, clearwindow=True,
                xlog=False,
                ylog=False,
                alpha=None,
                linewidths=None,
                linestyles='solid',
                colors=None,
                label=None,
                ):
        """Draw 2D contour data.

        Parameters
        ----------
        x0 : array-like
            independent axis in the first dimenation (on regular grid, flattened)
        x1 : array-like
            independent axis in the second dimenation (on regular grid, flattened)
        y : array-like
            dependent axis (i.e. image values) (on regular grid, flattened)
        {kwargs}
        """
        x0 = numpy.unique(x0)
        x1 = numpy.unique(x1)
        y = numpy.asarray(y)

        if x0.size * x1.size != y.size:
            raise NotImplementedErr('contourgrids')

        y = y.reshape(x1.size, x0.size)

        axes = self.setup_axes(overcontour, clearwindow)

        # Set up the axes
        if not overcontour:
            self.setup_plot(axes, title, xlabel, ylabel, xlog=xlog, ylog=ylog)

        if levels is None:
            axes.contour(x0, x1, y, alpha=alpha,
                         colors=colors,
                         linewidths=linewidths,
                         linestyles=linestyles)
        else:
            axes.contour(x0, x1, y, levels, alpha=alpha,
                         colors=colors,
                         linewidths=linewidths,
                         linestyles=linestyles)
        handles, _ = axes.get_legend_handles_labels()
        if len(handles) > 0:
            axes.legend()

    @add_kwargs_to_doc(kwargs_doc)
    @translate_args
    def image(self, x0, x1, y, *,
               aspect='auto',
               title=None, xlabel=None, ylabel=None,
               clearwindow=True,
               overplot=False,
               **kwargs):
        """Draw 2D image data.

        .. warning::
           This function is a non-functional dummy. The documentation is provided
           as a template only.

        Parameters
        ----------
        x0 : array-like
            independent axis in the first dimension
        x1 : array-like
            independent axis in the second dimension
        y : array-like, with shape (len(x0), len(x1))
            dependent axis (i.e. image values) in 2D
            with shape (len(x0), len(x1))
        {kwargs}
        """

        axes = self.setup_axes(overplot, clearwindow)

        # Set up the axes
        if not overplot:
            self.setup_plot(axes, title, xlabel, ylabel)

        extent = (x0[0], x0[-1], x1[0], x1[-1])
        im = axes.imshow(y, origin='lower', extent=extent, aspect=aspect)

        # TODO: This should be optional and have parameters.
        # but, for now, it's only used in _repr_html_ for DataIMG,
        # so can hardcode.
        plt.colorbar(im, ax=axes)


    def set_subplot(self, row, col, nrows, ncols, clearaxes=True,
                    left=None,
                    right=None,
                    bottom=None,
                    top=None,
                    wspace=0.3,
                    hspace=0.4):
        """Select a plot space in a grid of plots or create new grid

        This method adds a new subplot in a grid of plots.

        Parameters
        ----------
        row, col : int
            index (starting at 0) of a subplot in a grid of plots
        nrows, ncols : int
            Number of rows and column in the plot grid
        clearaxes : bool
            If True, clear entire plotting area before adding the new
            subplot.
        """
        plt.subplots_adjust(left=left, right=right, bottom=bottom, top=top,
                            wspace=wspace, hspace=hspace)

        # historically there have been issues with numpy integers here
        nrows = int(nrows)
        ncols = int(ncols)
        num = int(row) * ncols + int(col) + 1

        plt.subplot(nrows, ncols, num)
        if clearaxes:
            plt.cla()

    def set_jointplot(self, row, col, nrows, ncols, create=True,
                      top=0, ratio=2):
        """Move to the plot, creating them if necessary.

        Parameters
        ----------
        row : int
            The row number, starting from 0.
        col : int
            The column number, starting from 0.
        nrows : int
            The number of rows.
        ncols : int
            The number of columns.
        create : bool, optional
            If True then create the plots
        top : int
            The row that is set to the ratio height, numbered from 0.
        ratio : float
            The ratio of the height of row number top to the other
            rows.
        """

        if (row < 0) or (row >= nrows):
            emsg = f'Invalid row number {row} (nrows={nrows})'
            raise ArgumentErr(emsg)

        if (col < 0) or (col >= ncols):
            emsg = f'Invalid column number {col} (ncols={ncols})'
            raise ArgumentErr(emsg)

        plotnum = row * ncols + col
        if create:
            # We only change the vertical spacing
            ratios = [1] * nrows
            ratios[top] = ratio
            gs = {'height_ratios': ratios}
            fig, axes = plt.subplots(nrows, ncols, sharex=True, num=1,
                                     gridspec_kw=gs)
            fig.subplots_adjust(hspace=0.05)

            # Change all but the bottom row. By setting sharex
            # this is probably unneeded as get_xticklabels is empty.
            #
            if axes.ndim == 2:
                axes = axes[:-1, :].flatten()
            else:
                axes = axes[:-1]

            plt.setp([a.get_xticklabels() for a in axes[:-1]],
                     visible=False)

        try:
            ax = plt.gcf().axes[plotnum]
        except IndexError:
            emsg = f"Unable to find a plot with row={row} col={col}"
            raise ArgumentErr(emsg) from None

        plt.sca(ax)

    def set_title(self, title: str) -> None:
        """Change the display title.

        Parameters
        ----------
        title : str
           The title text to use.

        """

        plt.title(title)

    # HTML representation as SVG plots

    def as_svg(self, func):
        """Create HTML representation of a plot

        The output is a SVG representation of the data, as a HTML
        svg element, or a single div pointing out that the
        plot object has not been prepared.

        Parameters
        ----------
        func : function
            The function, which takes no arguments, which will create the
            plot. It creates and returns the Figure.

        Returns
        -------
        plot : str or None
            The HTML, or None if there was an error (e.g. prepare not
            called).

        """

        svg = io.StringIO()
        try:
            fig = func()
            fig.savefig(svg, format='svg')
            plt.close(fig)

        except Exception as e:
            logger.debug("Unable to create SVG plot: %s", str(e))
            return None

        # strip out the leading text so this can be used in-line
        svg = svg.getvalue()
        idx = svg.find('<svg ')
        if idx == -1:
            logger.debug("SVG output does not contain '<svg ': %s", svg)
            return None

        return svg[idx:]

    def as_html_plot_or_contour(self, data, summary=None, func="plot"):
        """Create HTML representation of a contour

        The output is a SVG representation of the data, as a HTML
        svg element.

        Parameters
        ----------
        data : Contour instance
            The contour object to display. It has already had its prepare
            method called.
        summary : str or None, optional
            The summary of the detail. If not set then the data type
            name is used.
        func : str, optional
            The function to call on the data object to make the plot.

        Returns
        -------
        plot : str or None
            The HTML, or None if there was an error (e.g. prepare not
            called).

        """

        def plotfunc():
            fig = plt.figure()
            getattr(data, func)()
            return fig

        svg = self.as_svg(plotfunc)
        if svg is None:
            return None

        if summary is None:
            summary = type(data).__name__

        ls = [formatting.html_svg(svg, summary)]
        return formatting.html_from_sections(data, ls)

    as_html_plot = partialmethod(as_html_plot_or_contour, func="plot")
    as_html_contour = partialmethod(as_html_plot_or_contour, func="contour")
    as_html_image = as_html_plot
    as_html_histogram = as_html_plot
    as_html_pdf = as_html_plot
    as_html_cdf = as_html_plot
    as_html_lr = as_html_plot
    as_html_data = as_html_plot
    as_html_datacontour = as_html_contour
    as_html_model = as_html_plot
    as_html_modelcontour = as_html_contour
    as_html_fit = as_html_plot
    as_html_fitcontour = as_html_contour
    as_html_contour1d = as_html_plot
    as_html_contour2d = as_html_contour

    # Needed for datastack plotting wrapper
    def initialize_plot(self, dataset, ids):
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

    def select_plot(self, dataset, ids):
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
