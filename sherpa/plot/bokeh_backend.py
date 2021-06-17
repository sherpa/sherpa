#
#  Copyright (C) 2021
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

The sherpa plotting API is procedural and most of the commands act on the
"current plot" or "current figure". This fits in very well with the matplotlib.pyplot model, but requires some extra work for object-oriented plotting packages.

bokeh is such an object-oriented package, which, by itself, does not keep a plotting state. Usually, bokeh commands act on an axis object that is passed in as a parameter. Sherpa, on the other hand, does not keep track of those objects, it expects the plotting package to know what the "current" plot is. In this module, we solve this problem using a few module-level variables, which hold the "currently active plot".
'''

import io
import logging

import numpy

from bokeh.plotting import figure, show
from bokeh.models.scales import LinearScale, LogScale
from bokeh.models.annotations import Span
from bokeh.models import Whisker

from sherpa.utils import get_keyword_defaults
from sherpa.utils.err import ArgumentErr, NotImplementedErr
from sherpa.utils import formatting
from sherpa.plot.utils import histogram_line


__all__ = ('clear_window', 'point', 'plot', 'histo', 'contour', 'set_subplot',
           'get_split_plot_defaults', 'get_plot_defaults', 'begin',
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
name = 'bokeh'

logger = logging.getLogger(__name__)

# These will be used a global variables.
# Yuck. 
current_fig = []


def begin():
    pass


def end():
    show(current_fig)


def exceptions():
    pass


_errorbar_defaults = get_keyword_defaults(plt.errorbar)


def clear_window():
    global current_fig
    fig = figure()
    current_fig = {'nrows': 1, 'ncols': 1,
                   'current_axis': fig,
                   'all_axes': [fig]}


# Now, that the function is defined, we call it to make sure we have a window to start
clear_window()


def set_window_redraw(redraw):
    pass


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
    if not overplot and clearwindow:
        clear_window()

    return current_fig['current_axis']


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
    axes.x_scale = LogScale() if xlog else LinearScale()
    axes.y_scale = LogScale() if ylog else LinearScale()

    if title:
        axes.title.text = title
    if xlabel:
        axes.xaxis.label = xlabel
    if ylabel:
        axes.yaxis.label = ylabel


def point(x, y, overplot=True, clearwindow=False,
          symbol=None, alpha=None,
          color=None):

    # bokeh marker=None means "no marker", but in sherpa is means
    # "default maker"
    if symbol is None:
        symbol = 'circle'

    axis = setup_axes(overplot, clearwindow)
    axis.scatter(x, y, alpha=alpha, color=color, marker=symbol)


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

    x, y2 = histogram_line(xlo, xhi, y)
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

    axes = current_fig[-1]
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
    'noline': None,
    'solid': 'solid',
    'dot': 'dotted',
    'dash': 'dashed',
    'dotdash': 'dotdash',
}

_drawstyle_map = {
    'steps': 'before',
    'steps-pre': 'before',
    'steps-mid': 'center',
    'steps-post': 'after'
}


def _set_line(line, linecolor=None, linestyle=None, linewidth=None):
    """Apply the line attributes, if set.

    Parameters
    ----------
    line :
        The line to change
    linecolor, linestyle, linewidth : optional
        The attribute value or None.

    """

    def set(label, newval):
        func = getattr(line, 'set_' + label)
        func(newval)

    if linecolor is not None:
        line.glyph.line_color = linecolor

    if linestyle is not None:
        line.glyph.line_dash = _linestyle_map[linestyle]

    if linewidth is not None:
        line.glyph.line_width = linewidth


# Could add support for alpha
# xmin/xmax don't do anything. Do I need them for compatibility?
#
def vline(x, ymin=0, ymax=1,
          linecolor=None,
          linestyle=None,
          linewidth=None,
          overplot=False, clearwindow=True):

    axes = setup_axes(overplot, clearwindow)

    line = Span(location=x, dimension='height',
                line_width=linewidth,
                line_color=linecolor,
                line_dash=_linestyle_map[linestyle])
    axes.add_layout(line)


def hline(y, xmin=0, xmax=1,
          linecolor=None,
          linestyle=None,
          linewidth=None,
          overplot=False, clearwindow=True):

    axes = setup_axes(overplot, clearwindow)

    line = Span(location=y, dimension='width',
                line_width=linewidth,
                line_color=linecolor,
                line_dash=_linestyle_map[linestyle])
    axes.add_layout(line)


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

    objs = []
    # Rely on color-cycling to work for both the "no errorbar" and
    # "errorbar" case.
    # TODO: Check how color cycling in bokeh works
    if linestyle != 'noline':
        if drawstyle == 'default':
            objs.append(axes.line(x, y,
                                  line_dash=_linestyle_map[linestyle],
                                  line_color=color, alpha=alpha))
        else:
            objs(axes.steps(x, y,
                            mode=_drawstyle_map[drawstyle],
                            line_dash=_linestyle_map[linestyle],
                            line_color=color, alpha=alpha))
    if marker != 'None':
        objs.append(point(x, y, symbol=marker, alpha=alpha, color=color))

    if xerrorbars:
        if xerr is not None:
            xerr = xerr / 2.
        xwhisker = Whisker(base=y, lower=x - xerr, upper=x + xerr,
                           dimension='width',
                           alpha=alpha, line_color=ecolor)
        axes.add_layout(xwhisker)
        objs.append(xwhisker)

    if yerrorbars:
        ywhisker = Whisker(base=x, lower=y - yerr, upper=y + yerr,
                           dimension='width',
                           alpha=alpha, line_color=ecolor)
        axes.add_layout(ywhisker)
        objs.append(ywhisker)

    if xaxis:
        axes.axhline(y=0, xmin=0, xmax=1, color='k', linewidth=1)

    if ratioline:
        axes.axhline(y=1, xmin=0, xmax=1, color='k', linewidth=1)

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

    global current_fig
    if create:

        current_fig = {'nrows': nrows, 'ncols': ncols,
                       'all_axes': [figure() for n in range(nrows * ncols)]
                       }


        # We only change the vertical spacing
        ratios = [1] * nrows
        ratios[top] = ratio
        # TO-Do: implement ratios

        # Axis sharing
        for r in range(1, nrows):
            for c in range(ncols):
                current_fig['all_axes'][r * ncols + c].x_range = \
                    current_fig['all_axes'][c].x_range

    current_fig['current_axis'] = current_fig['all_axes'][plotnum]


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
