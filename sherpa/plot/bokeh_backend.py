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


color-cycling
'''

import io
import logging

import numpy

from bokeh.plotting import figure, show
from bokeh.models.scales import LinearScale, LogScale
from bokeh.models.annotations import Span
from bokeh.models import Whisker, Line, Scatter
from bokeh.models import ColumnDataSource
from bokeh.models.arrow_heads import TeeHead
# the next 3 are needed to update and dispaly notebook plots
# is output_notebook is set
from bokeh.io.state import curstate
from bokeh.layouts import gridplot
from bokeh.io.notebook import push_notebook

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
# Not really a way to pass them in, so make module-level variables
width = 250
height = 250

# Translations from sherpa/matplotlib parameter values to
# bokeh.
# For now, I'm collecting here, but will package into better
# form later and write as better spec.

# This is not a 1:1. Matplotlib has marker with no equivalent in bokeh
# (e.g. '_') and the other way around (e.g. 'triangle_pin').
# So try to translate those that do match. Importantly, translate a minimum
# set, defined to be:
# "so.+"
# Importantly: if not translated, just pass through, so that
# users can set backend-dependent options
_marker_map = {
    '*': 'star',
    'o': 'circle',
    'D': 'diamond',
    '.': 'dot',
    'h': 'hex',
    'v': 'inverted_triangle',
    '+': 'plus',
    's': 'square',
    '^': 'triangle',
    'x': 'x',
    '1': 'y',
    }

# Note: In matplotlib, None means "default". In bokeh, default is done by not
# passing in that specific keyword, so really, should just strip from
# kwargs and not translante.

def tr_marker(marker):
    '''
    This principle can be turned into a more general function that
    also applies to _linestyle map etc.
    And is actually applicable beyond bokeh.
    '''
    if marker in _marker_map:
        return _marker_map[marker]
    else:
        # Could check here that marker in bokeh.core.enums.MarkerType
        # but can also leave that validation to bokeh itself.
        return marker

def tr_markersize(size):
    if size is None:
        return 10
    else:
        return size


def tr_capsize(size):
    '''Defaults to no head'''
    if size is None:
        return 0
    else:
        return size
    
# pass though what's not in this list    
_linestyle_map = {
    'None': None,
    'noline': None,
    'solid': 'solid',  # Sherpa name happens to match
    'dot': 'dotted',
    'dash': 'dashed',
    'dotdash': 'dotdash',  # Sherpa name happens to match
    # Simple matplotlib specifiers
    '-': 'solid',
    ':': 'dotted',
    '--': 'dashed',
    '-.': 'dotdash',
    '': None,
}

def tr_linestyle(ls):
    if ls in _linestyle_map:
        return _linestyle_map[ls]
    return ls

_drawstyle_map = {
    'steps': 'before',
    'steps-pre': 'before',
    'steps-mid': 'center',
    'steps-post': 'after'
}

def tr_alpha(alpha):
    '''In matplotlib alpha=None means "default". 
    That is almost always 1 and bokeh has no special value
    for the default, so we set it to 1.
    '''
    if alpha is None:
        return 1.
    else:
        return alpha

def tr_color(color):
    if color is None:
        return '#1f77b4'
    # translante for common matplotlib shotcuts rgbymck ?
    return color


def begin():
    pass


def end():
    if len(current_fig['all_axes']) == 1:
        show(current_fig['all_axes'][0])
    else:
        raise NotImplementedError


def exceptions():
    pass


def clear_window():
    global current_fig
    fig = figure()
    current_fig = {'nrows': 1, 'ncols': 1,
                   'current_axis': fig,
                   'all_axes': [[fig]],
                   'handle': None}

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
        axes.xaxis.axis_label = xlabel
    if ylabel:
        axes.yaxis.axis_label = ylabel

def notebook():
    if curstate().notebook:
        handle = current_fig.get('handle', None)
        if handle is None:
            current_fig['handle'] = show(gridplot(current_fig['all_axes'],
                                                  width=width, height=height),
                                         notebook_handle=True)
        else:
            push_notebook()

    
        
def point(x, y, overplot=True, clearwindow=False,
          symbol=None, alpha=None,
          color=None, **kwargs):

    # bokeh marker=None means "no marker", but in sherpa is means
    # "default maker"
    if symbol is None:
        symbol = 'circle'

    axes = setup_axes(overplot, clearwindow)
    source = ColumnDataSource({'x': x, 'y': y})
    for a, val in [('alpha', tr_alpha(alpha)),
                    ('color', tr_color(color))]:
        for b in ['fill', 'hatch', 'line']:
            if b + '_' + a not in kwargs:
                   kwargs[b + '_' + a] = val
                
    glyph = Scatter(marker=tr_marker(symbol),
                    **kwargs)
    axes.add_glyph(source, glyph)
    notebook()
    return glyph


def histo(xlo, xhi, y, yerr=None, title=None, xlabel=None, ylabel=None,
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
    if linestyle is None:
        linestyle = 'solid'
        
    x, y2 = histogram_line(xlo, xhi, y)
    # Note: this handles clearing the plot if needed.
    # Note: No translation needed here. That is done inside plot.
    #import pbd
    #pdb.set_trace()
    obj1 = plot(x, y2, yerr=None, xerr=None,
                title=title, xlabel=xlabel, ylabel=ylabel,
                overplot=overplot, clearwindow=clearwindow,
                xerrorbars=False, yerrorbars=False,
                xlog=xlog, ylog=ylog,
                linestyle=linestyle,
                drawstyle=drawstyle,
                color=color, marker='None', alpha=alpha,
                xaxis=False, ratioline=False)

    # Draw points and error bars at the mid-point of the bins.
    # Use the color used for the data plot: should it
    # also be used for marker[face]color and ecolor?
    #
    xmid = 0.5 * (xlo + xhi)
    xerr = (xhi - xlo) if xerrorbars else None
    yerr = yerr if yerrorbars else None

    #if marker is None:
    #    marker = 'o'
    obj2 = plot(xmid, y, xerr=xerr, yerr=yerr,
                overplot=True, clearwindow=False,
                xerrorbars=xerrorbars, yerrorbars=yerrorbars,
                linestyle='None',
                #color=obj1[0].line_color,
                marker=marker, alpha=tr_alpha(alpha),
                ecolor=ecolor, capsize=capsize,
                xaxis=False, ratioline=False)
    
    notebook()


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
        line.glyph.line_dash = tr_linestyle(linestyle)

    if linewidth is not None:
        line.glyph.line_width = linewidth


# Could add support for alpha
# xmin/xmax don't do anything. Do I need them for compatibility?
def vline(x, ymin=0, ymax=1,
          linecolor=None,
          linestyle=None,
          linewidth=None,
          overplot=False, clearwindow=True):

    axes = setup_axes(overplot, clearwindow)

    line = Span(location=x, dimension='height',
                line_width=linewidth,
                line_color=tr_color(linecolor),
                line_dash=tr_linestyle(linestyle),
    )
    axes.add_layout(line)
    notebook()


def hline(y, xmin=0, xmax=1,
          linecolor=None,
          linestyle=None,
          linewidth=None,
          overplot=False, clearwindow=True):

    axes = setup_axes(overplot, clearwindow)

    line = Span(location=y, dimension='width',
                line_width=linewidth,
                line_color=tr_color(linecolor),
                line_dash=tr_linestyle(linestyle),
                )
    axes.add_layout(line)
    notebook()


def plot(x, y, yerr=None, xerr=None, title=None, xlabel=None, ylabel=None,
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
         xaxis=False,
         ratioline=False,
         linecolor=None):
    """Draw x,y data.

    Note that the linecolor is not used, and is only included
    to support old code that may have set this option.
    """

    # If plotting backend where classes, which head to use could be a class
    # attribute
    arrowhead = TeeHead(size=tr_capsize(capsize),
                        # Not doing linestyle here. Let's just draw those as lines.
                        # But we could...
                        line_color=tr_color(color),
                        line_alpha=tr_alpha(alpha),
                        )

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
    source = ColumnDataSource({'x': x, 'y': y})
    
    if linestyle not in ['None', 'noline']:
        if drawstyle == 'default':
            line = Line(x='x', y='y',
                        line_dash=tr_linestyle(linestyle),
                        line_color=tr_color(color),
                        line_alpha=tr_alpha(alpha))
            axes.add_glyph(source, line)
            objs.append(line)
        else:
            # steps does not return line, so following will fail
            objs(axes.steps(x, y,
                            mode=_drawstyle_map[drawstyle],
                            line_dash=tr_linestyle(linestyle),
                            line_color=tr_color(color), alpha=tr_alpha(alpha)))
    if marker not in ('None', None):
         objs.append(point(x, y, symbol=tr_marker(marker),
                          alpha=tr_alpha(alpha), color=tr_color(color),
                          size=tr_markersize(markersize),
                          fill_color=tr_color(markerfacecolor)
        ))

    if xerrorbars and xerr is not None:
        source.data['lower_x'] = x - xerr
        source.data['upper_x'] = x + xerr
        xwhisker = Whisker(base='y', lower='lower_x', upper='upper_x',
                           source=source,
                           dimension='width',
                           line_alpha=tr_alpha(alpha),
                           line_color=tr_color(ecolor),
                           lower_head=arrowhead,
                           upper_head=arrowhead,
        )
        axes.add_layout(xwhisker)
        objs.append(xwhisker)

    if yerrorbars and yerr is not None:
        source.data['lower_y'] = y - yerr
        source.data['upper_y'] = y + yerr
        ywhisker = Whisker(base='x', lower='lower_y', upper='upper_y',
                           source=source,
                           dimension='height',
                           line_alpha=tr_alpha(alpha),
                           line_color=tr_color(ecolor),
                           lower_head=arrowhead,
                           upper_head=arrowhead,
        )
        axes.add_layout(ywhisker)        
        objs.append(ywhisker)

    if xaxis:
        axes.axhline(y=0, xmin=0, xmax=1, color='black', linewidth=1)

    if ratioline:
        axes.axhline(y=1, xmin=0, xmax=1, color='black', linewidth=1)

    notebook()
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
    notebook()


def set_subplot(row, col, nrows, ncols, clearaxes=True,
                left=None,
                right=None,
                bottom=None,
                top=None,
                wspace=0.3,
                hspace=0.4):
    # lefft right, ... seems to be matplotlib specific.
    # Can we remove that from the function signature?

    global current_fig
    
    if (nrows != current_fig['nrows']) or (ncols != current_fig['ncols']):
        clear_window()
        current_fig = {'nrows': nrows, 'ncols': ncols,
                       'current_axis': None,
                       # The following does not work, because it creates
                       # multiple references to the same None
                       #'all_axes': [None * ncols] * nrows
                       'all_axes': [[None for i in range(ncols)] for j in range(nrows)]}
    if clearaxes or current_fig['all_axes'][row][col] is None:
        fig = figure()
        current_fig['all_axes'][row][col] = fig
    current_fig['current_axis'] = current_fig['all_axes'][row][col]


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

    d['linestyle'] = 'solid'
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

    d['linestyle'] = 'solid'
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
    # As of bokeh 2.4 latex needs to be entire line
    # not in the middle of the line.
    # Not worth working around it on our side,
    # because bokeh is already working on that and will likely be released
    # before this sherpa code is used.
    return "$${}$$".format(txt)


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
