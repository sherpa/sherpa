#
#  Copyright (C) 2010, 2015, 2017, 2019, 2020, 2021
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

import numpy

from matplotlib import pyplot as plt

from sherpa.utils import get_keyword_defaults
from sherpa.utils.err import ArgumentErr, NotImplementedErr
from sherpa.utils import formatting


__all__ = ('clear_window', 'point', 'plot', 'histo', 'contour', 'set_subplot',
           'init', 'get_split_plot_defaults', 'get_plot_defaults', 'begin',
           'end', 'get_data_plot_defaults', 'get_model_plot_defaults',
           'exceptions', 'get_fit_plot_defaults', 'get_resid_plot_defaults',
           'get_ratio_plot_defaults', 'get_contour_defaults',
           'get_data_contour_defaults', 'get_model_contour_defaults',
           'get_fit_contour_defaults', 'get_resid_contour_defaults',
           'get_ratio_contour_defaults', 'get_confid_plot_defaults',
           'get_confid_contour_defaults', 'set_window_redraw', 'set_jointplot',
           'get_model_histo_defaults', 'get_histo_defaults',
           'get_component_plot_defaults', 'get_component_histo_defaults',
           'vline', 'hline', 'get_scatter_plot_defaults',
           'get_cdf_plot_defaults', 'get_latex_for_string', 'name')


# Name of the backend
name = 'pylab'

logger = logging.getLogger(__name__)


def init():
    pass


def begin():
    pass


def end():
    set_window_redraw(True)
    if plt.isinteractive():
        plt.draw()


def exceptions():
    pass


_errorbar_defaults = get_keyword_defaults(plt.errorbar)


def clear_window():
    plt.clf()


def set_window_redraw(redraw):
    if redraw:
        plt.draw()


def setup_axes(overplot, clearwindow):
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
        clear_window()

    return plt.gca()


def setup_plot(axes, title, xlabel, ylabel, xlog=False, ylog=False):
    """Basic plot setup.

    Parameters
    ----------
    axes
        The plot axes (output of setup_axes).
    title, xlabel, ylabel : str or None
        The plot, x-axis, and y-axis titles. They are skipped if
        the empty string or None.
    xlog , ylog : bool
        Should the scale be logarithmic (True) or linear (False)?

    """

    axes.set_xscale('log' if xlog else 'linear')
    axes.set_yscale('log' if ylog else 'linear')

    if title:
        axes.set_title(title)
    if xlabel:
        axes.set_xlabel(xlabel)
    if ylabel:
        axes.set_ylabel(ylabel)


def find_zorder(axes):
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


def point(x, y, overplot=True, clearwindow=False,
          symbol=None, alpha=None,
          color=None):

    axes = setup_axes(overplot, clearwindow)

    if color is None:
        style = '{}'.format(symbol)
    else:
        style = '{}{}'.format(color, symbol)

    axes.plot(numpy.array([x]), numpy.array([y]), style, alpha=alpha)


def histo(xlo, xhi, y, yerr=None, title=None, xlabel=None, ylabel=None,
          overplot=False, clearwindow=True,
          xerrorbars=False,
          yerrorbars=False,
          ecolor=_errorbar_defaults['ecolor'],
          capsize=_errorbar_defaults['capsize'],
          barsabove=_errorbar_defaults['barsabove'],
          xlog=False,
          ylog=False,
          linestyle='solid',
          drawstyle='default',
          color=None,
          alpha=None,
          marker='None',
          markerfacecolor=None,
          markersize=None,
          linecolor=None):
    """Draw histogram data.

    The histogram is drawn as horizontal lines connecting the
    start and end points of each bin, with vertical lines connecting
    consecutive bins. Non-consecutive bins are drawn with a
    (Nan,NaN) between them so no line is drawn connecting them.

    Points are drawn at the middle of the bin, along with any
    error values.

    Note that the linecolor is not used, and is only included
    to support old code that may have set this option.

    Notes
    -----
    It is not clear if the same zorder problems that hit
    plot could be relevant here.

    """

    if linecolor is not None:
        logger.warning("The linecolor attribute ({}) is unused.".format(linecolor))

    # Draw the data as a histogram, manually creating the lines
    # from the low to high edge of each bin. An alternative
    # would be to create RectanglePatches, one for each bin,
    # but I don't want each bin to go down to 0. I do not find
    # the existing drawstyle options to be sufficient.
    #
    # See https://stackoverflow.com/questions/5347065/interweaving-two-numpy-arrays#5347492
    # for interleaving arrays.
    #
    def intersperse(a, b):
        out = numpy.empty((a.size + b.size, ), dtype=a.dtype)
        out[0::2] = a
        out[1::2] = b
        return out

    x = intersperse(xlo, xhi)
    y2 = intersperse(y, y)

    # We need to identify non-consecutive bins, as a nan needs
    # to be added between the segments.
    #
    idxs, = numpy.where(xhi[:-1] != xlo[1:])
    if idxs.size > 0:
        idxs = 2 * (idxs + 1)
        nans = [numpy.nan] * idxs.size

        # ensure the arrays are floats so we can add nan values
        #
        x = numpy.insert(x.astype(numpy.float64), idxs, nans)
        y2 = numpy.insert(y2.astype(numpy.float64), idxs, nans)

    # Note: this handles clearing the plot if needed.
    #
    objs = plot(x, y2, yerr=None, xerr=None,
                title=title, xlabel=xlabel, ylabel=ylabel,
                overplot=overplot, clearwindow=clearwindow,
                xerrorbars=False, yerrorbars=False,
                ecolor=ecolor, capsize=capsize, barsabove=barsabove,
                xlog=xlog, ylog=ylog,
                linestyle=linestyle,
                drawstyle=drawstyle,
                color=color, marker=None, alpha=alpha,
                xaxis=False, ratioline=False)

    # Draw points and error bars at the mid-point of the bins.
    # Use the color used for the data plot: should it
    # also be used for marker[face]color and ecolor?
    #
    xmid = 0.5 * (xlo + xhi)
    xerr = (xhi - xlo) / 2 if xerrorbars else None
    yerr = yerr if yerrorbars else None

    try:
        color = objs[0].get_color()
    except AttributeError:
        pass

    axes = plt.gca()
    zorder = find_zorder(axes)

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
                  zorder=zorder)


_linestyle_map = {
    'noline': ' ',
    'solid': '-',
    'dot': ':',
    'dash': '--',
    'dotdash': '-.',
}


def _set_line(line, linecolor=None, linestyle=None, linewidth=None):
    """Apply the line attributes, if set.

    Parameters
    ----------
    line : matplotlib.lines.Line2D
        The line to change
    linecolor, linestyle, linewidth : optional
        The attribute value or None.

    """

    def set(label, newval):
        func = getattr(line, 'set_' + label)
        func(newval)

    if linecolor is not None:
        set('color', linecolor)

    if linestyle is not None:
        set('linestyle', _linestyle_map[linestyle])

    if linewidth is not None:
        set('linewidth', linewidth)


# There is no support for alpha in the Plot.vline class
#
def vline(x, ymin=0, ymax=1,
          linecolor=None,
          linestyle=None,
          linewidth=None,
          overplot=False, clearwindow=True):

    axes = setup_axes(overplot, clearwindow)

    line = axes.axvline(x, ymin, ymax)
    _set_line(line, linecolor=linecolor, linestyle=linestyle,
              linewidth=linewidth)


# There is no support for alpha in the Plot.hline class
#
def hline(y, xmin=0, xmax=1,
          linecolor=None,
          linestyle=None,
          linewidth=None,
          overplot=False, clearwindow=True):

    axes = setup_axes(overplot, clearwindow)

    line = axes.axhline(y, xmin, xmax)
    _set_line(line, linecolor=linecolor, linestyle=linestyle,
              linewidth=linewidth)


def plot(x, y, yerr=None, xerr=None, title=None, xlabel=None, ylabel=None,
         overplot=False, clearwindow=True,
         xerrorbars=False,
         yerrorbars=False,
         ecolor=_errorbar_defaults['ecolor'],
         capsize=_errorbar_defaults['capsize'],
         barsabove=_errorbar_defaults['barsabove'],
         xlog=False,
         ylog=False,
         linestyle='solid',
         drawstyle='default',
         color=None,
         marker='None',
         markerfacecolor=None,
         markersize=None,
         alpha=None,
         xaxis=False,
         ratioline=False,
         linecolor=None):
    """Draw x,y data.

    Note that the linecolor is not used, and is only included
    to support old code that may have set this option.
    """

    if linecolor is not None:
        logger.warning("The linecolor attribute, set to {}, is unused.".format(linecolor))

    axes = setup_axes(overplot, clearwindow)

    # Set up the axes
    if not overplot:
        setup_plot(axes, title, xlabel, ylabel, xlog=xlog, ylog=ylog)

    # See the discussion in find_zorder
    zorder = None
    if overplot:
        zorder = find_zorder(axes)

    # Rely on color-cycling to work for both the "no errorbar" and
    # "errorbar" case.
    #
    # TODO: do we really need this, or can we just use errorbar
    #       for both cases?
    #
    if xerrorbars or yerrorbars:
        if markerfacecolor is None:
            markerfacecolor = color
        if xerr is not None:
            xerr = xerr / 2.
        xerr = xerr if xerrorbars else None
        yerr = yerr if yerrorbars else None
        objs = axes.errorbar(x, y, yerr, xerr,
                             color=color,
                             linestyle=linestyle,
                             marker=marker,
                             markersize=markersize,
                             markerfacecolor=markerfacecolor,
                             alpha=alpha,
                             ecolor=ecolor,
                             capsize=capsize,
                             barsabove=barsabove,
                             zorder=zorder)

    else:
        objs = axes.plot(x, y,
                         color=color,
                         alpha=alpha,
                         linestyle=linestyle,
                         drawstyle=drawstyle,
                         marker=marker,
                         markersize=markersize,
                         markerfacecolor=markerfacecolor,
                         zorder=zorder)

    """
    for var in ('linestyle', 'color', 'marker', 'markerfacecolor',
                'markersize'):
        val = locals()[var]
        if val is not None:
            print(' -- set_{}({})'.format(var, val))
            getattr(line, 'set_' + var)(val)
    """

    # Should the color for these lines be taken from the current axes?
    #
    # Using black (color='k') and the default line width (of 1) in
    # matplotlib 2.0 produces a slightly-discrepant plot, since the
    # axis (and tick labels) are drawn with a thickness of 0.8.
    #
    #
    try:
        lw = axes.spines['left'].get_linewidth()
    except (AttributeError, KeyError):
        lw = 1.0

    if xaxis:
        axes.axhline(y=0, xmin=0, xmax=1, color='k', linewidth=lw)

    if ratioline:
        axes.axhline(y=1, xmin=0, xmax=1, color='k', linewidth=lw)

    return objs


def contour(x0, x1, y, levels=None, title=None, xlabel=None, ylabel=None,
            overcontour=False, clearwindow=True,
            xlog=False,
            ylog=False,
            alpha=None,
            linewidths=None,
            colors=None):

    x0 = numpy.unique(x0)
    x1 = numpy.unique(x1)
    y = numpy.asarray(y)

    if x0.size * x1.size != y.size:
        raise NotImplementedErr('contourgrids')

    y = y.reshape(x1.size, x0.size)

    axes = setup_axes(overcontour, clearwindow)

    # Set up the axes
    if not overcontour:
        setup_plot(axes, title, xlabel, ylabel, xlog=xlog, ylog=ylog)

    if levels is None:
        axes.contour(x0, x1, y, alpha=alpha,
                     colors=colors, linewidths=linewidths)
    else:
        axes.contour(x0, x1, y, levels, alpha=alpha,
                     colors=colors, linewidths=linewidths)


def set_subplot(row, col, nrows, ncols, clearaxes=True,
                left=None,
                right=None,
                bottom=None,
                top=None,
                wspace=0.3,
                hspace=0.4):

    plt.subplots_adjust(left=left, right=right, bottom=bottom, top=top,
                        wspace=wspace, hspace=hspace)

    # historically there have been issues with numpy integers here
    nrows = int(nrows)
    ncols = int(ncols)
    num = int(row) * ncols + int(col) + 1

    plt.subplot(nrows, ncols, num)
    if clearaxes:
        plt.cla()


def set_jointplot(row, col, nrows, ncols, create=True,
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
        emsg = 'Invalid row number {} (nrows={})'.format(row, nrows)
        raise ArgumentErr(emsg)

    if (col < 0) or (col >= ncols):
        emsg = 'Invalid column number {} (ncols={})'.format(col, ncols)
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
        emsg = "Unable to find a plot with row={} col={}".format(row, col)
        raise ArgumentErr(emsg) from None

    plt.sca(ax)


def get_split_plot_defaults():
    return get_keyword_defaults(set_subplot, 1)


def get_plot_defaults():
    return get_keyword_defaults(plot, 7)


def get_point_defaults():
    return get_keyword_defaults(point, 2)


def get_histo_defaults():
    return get_keyword_defaults(histo, 6)


def get_confid_point_defaults():
    d = get_point_defaults()
    d['symbol'] = '+'
    return d


def get_data_plot_defaults():
    d = get_plot_defaults()

    d['yerrorbars'] = True
    d['linestyle'] = 'None'
    d['marker'] = '.'
    return d


def get_model_histo_defaults():
    d = get_histo_defaults()
    return d


def get_model_plot_defaults():
    d = get_plot_defaults()

    d['linestyle'] = '-'
    d['marker'] = 'None'

    return d


def get_fit_plot_defaults():
    return {}


def get_resid_plot_defaults():
    d = get_data_plot_defaults()
    d['xerrorbars'] = True
    d['capsize'] = 0
    # d['marker'] = '_'
    d['xaxis'] = True
    return d


def get_ratio_plot_defaults():
    d = get_data_plot_defaults()
    d['xerrorbars'] = True
    d['capsize'] = 0
    # d['marker'] = '_'
    d['ratioline'] = True
    return d


def get_confid_plot_defaults():
    d = get_plot_defaults()

    d['linestyle'] = '-'
    d['marker'] = 'None'
    return d


def get_contour_defaults():
    return get_keyword_defaults(contour, 6)


get_data_contour_defaults = get_contour_defaults
get_model_contour_defaults = get_contour_defaults


def get_fit_contour_defaults():
    return {}


get_confid_contour_defaults = get_data_contour_defaults
get_resid_contour_defaults = get_data_contour_defaults
get_ratio_contour_defaults = get_data_contour_defaults
get_component_plot_defaults = get_model_plot_defaults
get_component_histo_defaults = get_model_histo_defaults


def get_cdf_plot_defaults():
    d = get_model_plot_defaults()
    return d


def get_scatter_plot_defaults():
    d = get_data_plot_defaults()
    return d


def get_latex_for_string(txt):
    """Convert to LaTeX form for the matplotlib back end.

    Parameters
    ----------
    txt : str
        The text component in LaTeX form (e.g. r'\alpha^2'). It
        should not contain any non-LaTeX content.

    Returns
    -------
    latex : str
        The input text surrounded by $. Note that there's no
        attempt to protect any $ characters in txt.

    """

    return "${}$".format(txt)


# HTML representation as SVG plots
#

def as_svg(func):
    """Create HTML representation of a plot

    The output is a SVG representation of the data, as a HTML
    svg element, or a single div pointing out that the
    plot object has not been prepared,

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
        logger.debug("Unable to create SVG plot: {}".format(e))
        return None

    # strip out the leading text so this can be used in-line
    svg = svg.getvalue()
    idx = svg.find('<svg ')
    if idx == -1:
        logger.debug("SVG output does not contain '<svg ': {}".format(svg))
        return None

    return svg[idx:]


def as_html_plot(data, summary=None):
    """Create HTML representation of a plot

    The output is a SVG representation of the data, as a HTML
    svg element.

    Parameters
    ----------
    data : Plot instance
        The plot object to display. It has already had its prepare
        method called.
    summary : str or None, optional
        The summary of the detail. If not set then the data type
        name is used.

    Returns
    -------
    plot : str or None
        The HTML, or None if there was an error (e.g. prepare not
        called).

    """

    def plotfunc():
        fig = plt.figure()
        data.plot()
        return fig

    svg = as_svg(plotfunc)
    if svg is None:
        return None

    if summary is None:
        summary = type(data).__name__

    ls = [formatting.html_svg(svg, summary)]
    return formatting.html_from_sections(data, ls)


def as_html_contour(data, summary=None):
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

    Returns
    -------
    plot : str or None
        The HTML, or None if there was an error (e.g. prepare not
        called).

    """

    def plotfunc():
        fig = plt.figure()
        data.contour()
        return fig

    svg = as_svg(plotfunc)
    if svg is None:
        return None

    if summary is None:
        summary = type(data).__name__

    ls = [formatting.html_svg(svg, summary)]
    return formatting.html_from_sections(data, ls)


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
