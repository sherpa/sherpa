#
#  Copyright (C) 2023, 2024
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
'''
Plotting using bokeh.
'''

import logging
from collections import ChainMap
from functools import partialmethod
from itertools import cycle

import numpy as np

import bokeh
from bokeh.plotting import figure, show
from bokeh.models.scales import LinearScale, LogScale
from bokeh.models.annotations import Span
from bokeh.models import Whisker, Scatter
from bokeh.models import ColumnDataSource
from bokeh import embed
from bokeh.layouts import gridplot
from bokeh import palettes

from sherpa.utils.err import ArgumentErr, NotImplementedErr
from sherpa.utils import formatting
from sherpa.plot.utils import histogram_line

from sherpa.plot.backends import BasicBackend
from sherpa.plot.backends import kwargs_doc as orig_kwargs_doc
from sherpa.plot.backend_utils import (translate_args,
                                       add_kwargs_to_doc)

__all__ = ['BokehBackend', ]

logger = logging.getLogger(__name__)

# Most of the kwargs docstrings are still true, but some can be updated
# because bokeh does support a wider range than just the
# set of Sherpa backend-independent options.
updated_kwarg_docs = {
    'color': ['str or tuple', 'Any bokeh color'],
    'linecolor': ['string or tuple', 'Any bokeh color'],
    'marker': ['str',
               '''"None" (as a string, no marker shown), "" (empty string, no marker shown),
or any bokeh marker (see bokeh documentation).'''],
    'linestyle': ['str',
                  '''``'noline'``,
``'None'`` (as string, same as ``'noline'``),
``'solid'``, ``'dot'``, ``'dash'``, ``'dotdash'``, ``'-'`` (solid
line), ``':'`` (dotted), ``'--'`` (dashed), ``'-.'`` (dot-dashed),
``''`` (empty string, no line shown), `None` (default - usually
solid line) or any other bokeh linestyle.'''],
    'drawstyle': ['str', 'bokeh drawstyle'],
}

kwargs_doc = ChainMap(updated_kwarg_docs, orig_kwargs_doc)


class BokehBackend(BasicBackend):
    """Sherpa plotting backend for the bokeh package.

    :term:`bokeh` is a plotting library that generates a json representation of
    a plot, which is then rendered in a browser (either in a jupyter notebook or
    in a new tab for traditional python scripts). For the relatively simple cases
    supporting standard sherpa plotting commands, the json representation can be
    saved and displayed in a browser without the need for a running python kernel.
    Thus, plots like this can be embedded into webpages or as "interactive figures"
    in e.g. AAS journals.

    Notes
    -----
    The sherpa plotting API is procedural and most of the commands act on the
    "current plot" or "current figure". This fits in very well with the `matplotlib.pyplot`
    model, but requires some extra work for object-oriented plotting packages.

    :term:`bokeh` is such an object-oriented package, which, by itself,
    does not keep a plotting state. Usually, bokeh commands act
    on an axis object that is passed in as a parameter.
    Sherpa, on the other hand, does not keep track of those objects,
    it expects the plotting package to know what the "current" plot is.

    We solve this problem with attributes in the BokehBackend class:

    - `current_fig`: the current figure
    - `current_axis`: the current axis. This is a reference to one of the
      panels in the current figure.

    We follow a similar approach to default to cycling through colors
    like matplotlib does. The bokeh package does not have a similar
    mechanism, so we wrap `bokeh.plotting.figure` to add an additional
    attribute `_color_cycle` that is an `itertools.cycle` object
    cycling through the colors in the palette defined in the
    `palette` attribute of the BokehBackend class.

    """

    translate_colors = {
        "r": "red",
        "g": "green",
        "b": "blue",
        "k": "black",
        "w": "white",
        "c": "cyan",
        "y": "yellow",
        "m": "magenta",
    }

    translate_dict = {
        "color": translate_colors,
        "ecolor": translate_colors,
        "markerfacecolor": translate_colors,
        "linecolor": translate_colors,
        "linestyle": {
            "None": "noline",
            None: "solid",
            "solid": "solid",
            "dot": "dotted",
            "dash": "dashed",
            "longdash": "dashed",
            "-": "solid",
            ":": "dotted",
            "--": "dashed",
            "-.": "dotdash",
        },
        "linewidth": {
            None: 1,
        },
        "drawstyle": {
            "steps": "before",
            "steps-pre": "before",
            "steps-mid": "center",
            "steps-post": "after",
        },
        "marker": {
            #'None': 'circle',
            None: "",
            ".": "circle",
            "o": "circle",
            "+": "cross",
            "s": "square",
        },
        "markersize": {
            None: 5,
        },
        "alpha": {
            None: 1.0,
        },
        "capsize": {
            None: 0,
        },
        "linewidths": {None: 1},
        "linestyles": {None: "solid"},
        "colors": {None: "black"},
        "aspect": {"auto": None, "equal": 1.0},
    }
    '''Dict of keyword arguments that need to be translated for this backend.

    The keys in this dict are keyword arguments (e.g. ``'markerfacecolor'``)
    and the values are one of the following:
    - A dict where the keys are the backend-independent values and the values
    are the values expressed for this backend. For values not listed in the
    dict, no translation is done.
    - A callable. The callable translation function is called with any argument
    given and allows the backend arbitrary translations. It should, at the
    very least, accept all backend independent values for this parameter
    without error.

    In particular, sherpa uses None to mean "the default", while in bokeh
    None often means "don't show this", so None values are translated to
    a sensible default value (e.g. line thickness of 1)
    '''

    palette = palettes.Category10[10]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.clear_window()

    def __exit__(self, exec_type, value, traceback):
        '''Called from the UI after an interactive plot is done.'''
        return False

    def get_latex_for_string(self, txt):
        """Convert LaTeX formula

        Parameters
        ----------
        txt : str
            The text component in LaTeX form (e.g. r'\alpha^2'). It
            should not contain any non-LaTeX content.

        Returns
        -------
        latex : str
            The text modified as appropriate for a backend so that the LaTeX
            will be displayed properly.

        """
        return f"$${txt}$$"

    def _figure(self):
        fig = figure()
        fig._color_cycle = cycle(self.palette)
        return fig

    def clear_window(self):
        fig = self._figure()
        self.current_fig = gridplot([fig], ncols=1)
        self.current_axis = fig

    def setup_axes(self, overplot, clearwindow):
        """Return the axes object, creating it if necessary.

        Parameters
        ----------
        overplot : bool
        clearwindow : bool

        Returns
        -------
        axis
            The bokeh axes object.
        """
        if not overplot and clearwindow:
            self.clear_window()

        return self.current_axis

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
        axes.x_scale = LogScale() if xlog else LinearScale()
        axes.y_scale = LogScale() if ylog else LinearScale()
        # To-do: DO I need the "if" or does bokeh accept to set to None?
        if title:
            axes.title.text = title
        if xlabel:
            axes.xaxis.axis_label = xlabel
        if ylabel:
            axes.yaxis.axis_label = ylabel

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
              xlog=False,
              ylog=False,
              linestyle='None',
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
            logger.warning("The linecolor attribute (%s) is unused.", linecolor)

        x, y2 = histogram_line(xlo, xhi, y)
        # Note: this handles clearing the plot if needed.
        #
        objs = self.plot(x, y2, yerr=None, xerr=None,
                title=title, xlabel=xlabel, ylabel=ylabel,
                overplot=overplot, clearwindow=clearwindow,
                xerrorbars=False, yerrorbars=False,
                ecolor=ecolor, capsize=capsize,
                xlog=xlog, ylog=ylog,
                linestyle=linestyle,
                linewidth=linewidth,
                drawstyle=drawstyle,
                color=color, marker=None, alpha=alpha,
                xaxis=False, ratioline=False)

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
        except (AttributeError, IndexError):
            pass

        ebars = self.plot(xmid, y,
                          # Plot error y bars if present
                          xerrorbars=False,
                          yerrorbars=True,
                          yerr=yerr, xerr=xerr,
                          overplot=True,
                          clearwindow=False,
                          color=color,
                          alpha=alpha,
                          linestyle='noline',
                          marker=marker,
                          markersize=markersize,
                          markerfacecolor=markerfacecolor,
                          ecolor=ecolor,
                          capsize=capsize,
                         )

    # Note that this is an internal method that is not wrapped in
    # @ translate_args
    # This is called from methods like plot, histo, etc.
    # so the arguments passed to this method are already translated and we do
    # not want to apply the translation a second time.
    def _set_line(self, line, linecolor=None, linestyle=None, linewidth=None):
        """Apply the line attributes, if set.

        Parameters
        ----------
        line :
            The line to change
        linecolor, linestyle, linewidth : optional
            The attribute value or None.
        """
        if linecolor is not None:
            line.glyph.line_color = linecolor

        if linestyle is not None:
            line.glyph.line_dash = linestyle

        if linewidth is not None:
            line.glyph.line_width = linewidth

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

        line = Span(location=x, dimension='height',
                    line_width=linewidth,
                    line_color=linecolor,
                    line_dash=linestyle)
        axes.add_layout(line)

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

        line = Span(location=y, dimension='width',
                    line_width=linewidth,
                    line_color=linecolor,
                    line_dash=linestyle)
        axes.add_layout(line)

    @add_kwargs_to_doc(kwargs_doc)
    @translate_args
    def plot(self, x, y, *,
             yerr=None, xerr=None, title=None,
             xlabel=None, ylabel=None,
             overplot=False, clearwindow=True,
             xerrorbars=False,
             yerrorbars=False,
             ecolor=None,
             capsize=0,
             xlog=False,
             ylog=False,
             linestyle='solid',
             drawstyle='default',
             color=None,
             marker='None',
             markerfacecolor=None,
             markersize=1,
             alpha=1,
             label=None,
             linewidth=1,
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
        axes = self.setup_axes(overplot, clearwindow)

        # Set up the axes
        if not overplot:
            self.setup_plot(axes, title, xlabel, ylabel, xlog=xlog, ylog=ylog)


        if color is None:
            color = next(axes._color_cycle)
        if linecolor is not None:
            logger.warning("The linecolor attribute, set to {}, is unused.".format(linecolor))
        if markerfacecolor is None:
            markerfacecolor = color
        if ecolor is None:
            ecolor = color

        objs = []

        kwargs = {}
        if label is not None:
            kwargs['legend_label'] = label
        if linestyle != 'noline':
            if drawstyle == 'default':
                objs.append(axes.line(x, y,
                                      line_dash=linestyle,
                                      line_color=color,
                                      line_width=linewidth,
                                      alpha=alpha,
                                      **kwargs
                                      ))
            else:
                objs.append(axes.steps(x, y,
                                mode=drawstyle,
                                line_dash=linestyle,
                                line_color=color,
                                alpha=alpha,
                                **kwargs
                                ))
        if marker not in ('None', ''):
            source = ColumnDataSource({'x': np.atleast_1d(x),
                                       'y': np.atleast_1d(y)})
            glyph = Scatter(marker=marker,
                            line_alpha=alpha,
                            fill_alpha=alpha,
                            hatch_alpha=alpha,
                            line_color=markerfacecolor,
                            fill_color=markerfacecolor,
                            hatch_color=markerfacecolor,
                            size=markersize,
                            *kwargs
                            )
            axes.add_glyph(source, glyph)
            objs.append(glyph)

        if xerrorbars and xerr is not None:
            source = ColumnDataSource({'x': np.atleast_1d(x),
                                       'y': np.atleast_1d(y),
                                       'lower': np.atleast_1d(x - xerr),
                                       'upper': np.atleast_1d(x + xerr),
                                       })
            xwhisker = Whisker(
                source=source,
                base='y', lower='lower', upper='upper',
                               dimension='width',
                               line_alpha=alpha,
                               line_color=ecolor,
                               *kwargs
                               )
            xwhisker.upper_head.size=capsize
            xwhisker.lower_head.size=capsize
            axes.add_layout(xwhisker)
            objs.append(xwhisker)

        if yerrorbars and yerr is not None:
            # Are the errors symmetric or not? At the moment this is
            # decided by checking if yerr is a tuple or not.
            #
            if type(yerr) == tuple:
                ylo = np.atleast_1d(y - yerr[0])
                yhi = np.atleast_1d(y + yerr[1])
            else:
                ylo = np.atleast_1d(y - yerr)
                yhi = np.atleast_1d(y + yerr)

            source = ColumnDataSource({'x': np.atleast_1d(x),
                                       'y': np.atleast_1d(y),
                                       'lower': ylo,
                                       'upper': yhi,
                                       })
            ywhisker = Whisker(
                source=source,
                base='x', lower='lower', upper='upper',
                               dimension='height',
                               line_alpha=alpha, line_color=ecolor,
                               **kwargs
                               )
            ywhisker.upper_head.size=capsize
            ywhisker.lower_head.size=capsize
            axes.add_layout(ywhisker)
            objs.append(ywhisker)

        if xaxis:
            axes.axhline(y=0, xmin=0, xmax=1, color='k', linewidth=1)

        if ratioline:
            axes.axhline(y=1, xmin=0, xmax=1, color='k', linewidth=1)

        return objs

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
        x0 = np.unique(x0)
        x1 = np.unique(x1)
        y = np.asarray(y)

        if x0.size * x1.size != y.size:
            raise NotImplementedErr('contourgrids')

        y = y.reshape(x1.size, x0.size)

        axes = self.setup_axes(overcontour, clearwindow)

        # Set up the axes
        if not overcontour:
            self.setup_plot(axes, title, xlabel, ylabel, xlog=xlog, ylog=ylog)

        if levels is None:
            levels = 7
        if isinstance(levels, int):
            # Matplotlib has defaults for levels and sherpa can pass in "None" to mean
            # "use the default". Bokeh always requires to specify the levels, so
            # we have to define our own defaults here.
            # We try to pick something sensible, but not too complicated.
            levels = np.linspace(np.min(y), np.max(y), levels)[1: -1]


        axes.contour(x0, x1, y, levels, line_alpha=alpha,
                     line_color=colors, line_width=linewidths)

    @add_kwargs_to_doc(kwargs_doc)
    @translate_args
    def image(self, x0, x1, y, *,
               aspect=None,
               title=None, xlabel=None, ylabel=None,
               clearwindow=True,
               overplot=False,
               **kwargs):
        """Draw 2D image data.

        .. warning::
           This function is experimental and is currently only used in
           the rish display code.

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
        print(x0, x1, y, aspect, kwargs)
        axes = self.setup_axes(overplot, clearwindow)

        # Set up the axes
        if not overplot:
            self.setup_plot(axes, title, xlabel, ylabel)

        axes.image(image=[y], x=x0[0], y=x1[0], dw=x0[-1] - x0[0], dh=x1[-1] - x1[0],
                   level='image', palette="Sunset11")
        axes.aspect_ratio = aspect


    def _index_axis(self, row, col):
        '''Find index number of subplot in a grid of plots

        self.current_fig is a `bokeh.models.plots.GridPlot`
        object. It has a children attribute, which is a list of tuples, where
        each tuple is a (figure, row, col) tuple. This method returns the
        index of the tuple in the list that matches the given row and col or,
        if no such tuple exists, None.
        '''
        for i, fig in enumerate(self.current_fig.children):
            if fig[1] == row and fig[2] == col:
                return i
        return None

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
        index = self._index_axis(row, col)
        if index is not None and not clearaxes:
            # We found the right subplot, now make it the current one
            self.current_axis = self.current_fig.children[index]
            return

        # We found the subplot, but we want to replace it with a new one
        if index is not None and clearaxes:
            self.current_fig.children.pop(index)

        # At this stage, the subplot does not exist, so we make a new one
        # and add it to the grid.
        newf = self._figure()
        self.current_fig.children.append((newf, row, col))
        self.current_axis = newf



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
        # This function is written to work with an incomplete grid of plots,
        # so there is a whole lot of "if plotnum is not None" in here.
        # If we were to require that all grid are complete (i.e. a 3x2 grid)
        # always has 6 plots in a fixed order, the implementation would be
        # simpler.
        if create:
            figs = [self._figure() for n in range(nrows * ncols)]
            self.current_fig = gridplot(figs, ncols=ncols)
            self.current_axis = figs[row * ncols + col]

        # Space is dynamically allocated and lots of things can change the
        # height, e.g. a user may have changed bokeh defaults.
        # It's easy to change the height of a plot, but we first need to find
        # what height we consider "normal".
        # If we simply were to multiply the height of the top row by the ratio,
        # then the total height of the plot gets larger every time this method
        # is called, which, for ncols > 1, might be several times.
        # We pick the first plot in the grid that's not in row "top" and use
        # that as the "normal" height. That might fail for more complex grids,
        # but I declare that outside the scope of this function.
        height = None
        for fig in self.current_fig.children:
            if fig[1] != top:
                height = fig[0].height
                break
        # `height` is None if we only have one row, so we don't need to
        # change anything.
        if height is not None:
            for i in range(ncols):
                plotnum = self._index_axis(top, i)
                # We might have an incomplete grid of plots, and then
                # plotnum might be None.
                if plotnum is not None:
                    self.current_fig.children[plotnum][0].height = height * ratio

        # Axis sharing
        for c in range(ncols):
            basenum = None
            for r in range(nrows):
                index = self._index_axis(0, c)
                if basenum is None:
                    basenum = index
                else:
                    self.current_fig.children[index][
                        0
                    ].x_range = self.current_fig.children[basenum][0].x_range

        plotnum = self._index_axis(row, col)
        if plotnum is None:
            raise ArgumentErr(f"No plot at row {row} and column {col}." +
                              "Use create=True to create an entirely new, " +
                              "complete grid of plots.") from None
        self.current_axis = self.current_fig.children[plotnum][0]

    def set_title(self, title: str) -> None:
        """Change the display title.

        Parameters
        ----------
        title : str
           The title text to use.

        """

        # It is not at all obvious how to do this since is it the
        # first component or all components or ?
        #
        plotnum = self._index_axis(0, 0)
        if plotnum is None:
            # Should this error out?
            return

        self.current_fig.children[plotnum][0].title.text = title

    # HTML representation

    def as_html(self, func):
        """Create HTML representation of a plot

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
        try:
            func()
        except Exception as e:
            logger.debug("Unable to create bokeh plot: %s", str(e))
            return None

        image = embed.file_html(self.current_fig)

        # See https://docs.bokeh.org/en/latest/docs/reference/embed.html#bokeh.embed.components
        # for a list of other js files that might be needed.
        # Also, do we want to use this offline instead of CDN?
        return f'''
        <script src="https://cdn.bokeh.org/bokeh/release/bokeh-{bokeh.__version__}.min.js"></script>
        ''' + image


    def as_html_plot_or_contour(self, data, summary=None, plotting_func="plot"):
        """Create HTML representation of a plot

        Parameters
        ----------
        data : Plot instance
            The plot object to display. It has already had its prepare
            method called.
        summary : str or None, optional
            The summary of the detail. If not set then the data type
            name is used.
        plotting_func : str, optional
            The function to call on the data object to make the plot.

        Returns
        -------
        plot : str or None
            The HTML, or None if there was an error (e.g. prepare not
            called).
        """
        image = self.as_html(getattr(data, plotting_func))

        if image is None:
            return None

        if summary is None:
            summary = type(data).__name__

        ls = [formatting.html_svg(image, summary)]
        return formatting.html_from_sections(data, ls)

    as_html_plot = partialmethod(as_html_plot_or_contour, plotting_func="plot")
    as_html_contour = partialmethod(as_html_plot_or_contour, plotting_func="contour")
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
