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
'''Interface between Sherpa plotting objects and :term:`bokeh`

This module provides an interface between the Sherpa plot objects in
`sherpa.plot` and the :term:`bokeh` plotting package. The bokeh plotting
package formats data into fully interactive javascript, that can be displayed
in a web-browser. All plots allow interactive panning, zooming, scrolling and
the current state of thep lot can be exported to static file, e.g. a png.
Other interactive tools might be available depending on the plot type.

Currently active plot
=====================

The sherpa plotting API is procedural and most of the commands act on
the "current plot" or "current figure". This fits in very well with
the matplotlib.pyplot model, but requires some extra work for
object-oriented plotting packages.

:term:`bokeh` is such an object-oriented package, which, by itself,
does not keep a plotting state. Usually, bokeh commands act on an axis
object that is passed in as a parameter. Sherpa, on the other hand,
does not keep track of those objects, it expects the plotting package
to know what the "current" plot is. In this module, we solve this
problem using a module-level variable
`sherpa.plot.bokeh_backend.current_fig`, which holds the "currently
active plot". As a user, you can access this variable to obtain the
bokeh objects to futher modify its appearance or to output it in whole
or in parts to a file. In this example, we first prepare a Shera data
object and a Sherpa plot object (assuming that the bokeh backend is active)::

    >>> import numpy as np
    >>> from bokeh.plotting import show
    >>> from sherpa.data import Data1DInt
    >>> from sherpa.plot import DataPlot

    >>> edges = np.asarray([-10, -5, 5, 12, 17, 20, 30, 56, 60])
    >>> y = np.asarray([28, 62, 17, 4, 2, 4, 125, 55])
    >>> d = Data1DInt('example histogram', edges[:-1], edges[1:], y)

    >>> dplot = DataPlot()
    >>> dplot.prepare(d)
    >>> dplot.plot()

    >>> from sherpa.plot import backend
    >>> p = backend.current_fig['current_axes']
    >>> p.title.text = 'Fruit sizes'
    >>> p.title.text_font_size = "25px"
    >>> show(p)


Color-cycling
=============
Unlike :term:`matplotlib`, :term:`bokeh` does not cycle colors automatically
when new lines or symbols (collectively called "glyph" in bokeh) are added.
This means that data and model appear in the same color by default, but it
also makes distinguishing different models opr datasets hard, unless explicit
colors are set:

    >>> dplot.plot(color='red')

Use in the notebook
===================
For use in the `Jupyter notebook <https://jupyter.org/>`_ the following
command enables interactive output in every cell::

    >>> from bokeh.io import output_notebook
    >>> output_notebook()

See the `bokeh documentation <https://docs.bokeh.org/en/latest/docs/first_steps/first_steps_7.html>`_ 
for other output options.
'''
from contextlib import contextmanager
import io
import logging

import numpy as np

from bokeh.plotting import figure, show
from bokeh.models.scales import LinearScale, LogScale
from bokeh.models.annotations import Span
from bokeh.models import Whisker, Line, Scatter
from bokeh.models import ColumnDataSource
from bokeh.models.arrow_heads import TeeHead
# the next 3 are needed to update and display notebook plots
# if output_notebook is set
from bokeh.io.state import curstate
from bokeh.layouts import gridplot
from bokeh.io.notebook import push_notebook
from bokeh.embed import components

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
        return 20
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
    None: 'solid',
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
    # If we are in an active notebook, then the plot is already
    # displayed by the notebook specific function
    if not curstate().notebook:
        show(gridplot(current_fig['all_axes']))


def exceptions():
    pass


def set_subplot(row, col, nrows, ncols, clearaxes=True,
                left=None,
                right=None,
                bottom=None,
                top=None,
                wspace=0.3,
                hspace=0.4):
    # left right, ... seems to be matplotlib specific.
    # Can we remove that from the function signature?

    global current_fig
    
    if (nrows != current_fig['nrows']) or (ncols != current_fig['ncols']):
        current_fig = {'nrows': nrows, 'ncols': ncols,
                       # The following does not work, because it creates
                       # multiple references to the same None
                       #'all_axes': [None * ncols] * nrows
                       'all_axes': [[None for i in range(ncols)] for j in range(nrows)]}

    if clearaxes or current_fig['all_axes'][row][col] is None:
        fig = figure(width=width, height=height)
        current_fig['all_axes'][row][col] = fig
        # Make layout as late as possible on first display (in the notebook
        # function) but if it exisits already, update with new added plot.
        # Directly replacing children is not the best way, but bokeh
        # does not have a clean_figure() function and this seems to work
        # for the limited layouts that sherpa uses (square grids).
        
        # Also, it would be more efficient to replace the exiting child in
        # (row, col) instead of appending even if that (row, col) combination
        # might exisit already,
        # but this is a rare use case in sherpa and efficiency is not our
        # main concern here.
        if 'layout' in current_fig:
            current_fig['layout'].children[1].children.append((fig, row, col))
        
    current_fig['current_axes'] = current_fig['all_axes'][row][col]



def clear_window():
    global current_fig
    current_fig = {'nrows': 0, 'ncols': 0}
    set_subplot(0, 0, 1, 1)

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

    return current_fig['current_axes']


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
    global current_fig
    if curstate().notebook:
        if 'layout' not in current_fig:
            current_fig['layout'] = gridplot(current_fig['all_axes'])
        handle = current_fig.get('handle', None)
        if handle is None:
            current_fig['handle'] = show(current_fig['layout'],
                                         notebook_handle=True)
        else:
            push_notebook(handle=handle)

    
        
def point(x, y, overplot=True, clearwindow=False,
          symbol=None, alpha=None,
          color=None, **kwargs):

    # bokeh marker=None means "no marker", but in sherpa is means
    # "default maker"
    if symbol is None:
        symbol = 'circle'

    axes = setup_axes(overplot, clearwindow)
    # 2D will fail later, but we need at least one 1D, since
    # bokeh does not handle scalars here
    source = ColumnDataSource({'x': np.atleast_1d(x),
                               'y': np.atleast_1d(y)})
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
    xerr = (xhi - xlo) / 2 if xerrorbars else None
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
            # steps does return None)
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
        objs.append(hline(overplot=True, y=0, xmin=0, xmax=1,
                          linecolor='black', linewidth=1))

    if ratioline:
        objs.append(hline(overplot=True, y=1, xmin=0, xmax=1,
                          linecolor='black', linewidth=1))
    # It is possible to go through this method without plotting anything.
    # That happens when a histogram with "linestyle='None'" etc. and no
    # errors is plotted. In that case, updating the notebook display
    # may raise a warning, so don't do that.
    if objs:
        notebook()
    return objs


def contour(x0, x1, y, levels=None, title=None, xlabel=None, ylabel=None,
            overcontour=False, clearwindow=True,
            xlog=False,
            ylog=False,
            alpha=None,
            linewidths=None,
            colors=None):

    '''Contour plots are not available natively in bokeh yet, see
    https://discourse.bokeh.org/t/contouring-in-bokeh/8517

    However, we can address many use cases by plotting an image instead.'''

    if overcontour:
        raise NotImplementedError('Bokeh does not support contours, so instead levels are displayed as images. Try plotting the images first and then overlay other plot elements on top of that.')
    x0 = np.unique(x0)
    x1 = np.unique(x1)
    y = np.asarray(y)

    if x0.size * x1.size != y.size:
        raise NotImplementedErr('contourgrids')

    image = y.reshape(x1.size, x0.size)

    axes = setup_axes(overcontour, clearwindow)

    # Set up the axes
    if not overcontour:
        setup_plot(axes, title, xlabel, ylabel, xlog=xlog, ylog=ylog)

    p = current_fig['current_axes']
    p.x_range.range_padding = p.y_range.range_padding = 0

    # Need images on regular grid
    dx0 = np.diff(x0)
    dx1 = np.diff(x1)
    if np.allclose(dx0, dx0[0]) and (np.allclose(dx1, dx1[0])):
        x = x0[0]
        y = x1[0]
        dw = x0[-1] - x0[0]
        dh = x1[-1] - x1[0]
    else:
        raise NotImplementedError('Interpolation of non-regular grids is not implemented yet.')
    
    # must give a vector of image data for image parameter
    p.image(image=[image], x=x, y=y, dw=dw, dh=dh,
            palette="Spectral11", level="image", alpha=tr_alpha(alpha))

    notebook()




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

    set_subplot(row, col, nrows, ncols, clearaxes=create)
    if create:
        for plot in current_fig['all_axes'][top]:
            plot.height = height * ratio

        # Axis sharing
        for r in range(1, nrows):
            for c in range(ncols):
                if current_fig['all_axes'][r][c] is None:
                    current_fig['all_axes'][r][c] = figure(height=height, width=width)
                current_fig['all_axes'][r][c].x_range = \
                    current_fig['all_axes'][top][c].x_range


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




@contextmanager
def supress_notebook_plotting():
    '''Suppress notebook plotting independent of what the user requested

    We do not know if the user set a notebook output. However, we do know
    that in some situations we do not want an output in the notebook cell
    or a saved plot to appear, for example when we render an object just
    for the _repr_html_ display in the notebook.

    This contex manager temporarily resets the notebook output to be 
    silent.
    '''
    # Code to acquire resource, e.g.:
    notebook_state = curstate().notebook
    file_state = curstate().file
    curstate()._notebook = False
    curstate()._file = None
    try:
        yield None
    finally:
        # Code to release resource, e.g.:
        curstate()._notebook = notebook_state
        curstate()._file = file_state


def as_html_plot(data, summary=None):
    """Create HTML representation of a plot

    The output is a HTML representation of the data, as a HTML
    using a <script> element that renders into a <div> element.

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
    with supress_notebook_plotting():
        data.plot()
    
    if 'layout' not in current_fig:
        current_fig['layout'] = gridplot(current_fig['all_axes'])

    script, div = components(current_fig['layout'])

    if summary is None:
        summary = type(data).__name__

    ls = [formatting.html_svg(div, summary)]
    return formatting.html_from_sections(data, [script] + ls)



def as_html_contour(data, summary=None):
    """Create HTML representation of a contour

    The output is a HTML representation of the data, as a HTML
    using a <script> element that renders into a <div> element.

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
    with supress_notebook_plotting():
        data.contour()
    
    if 'layout' not in current_fig:
        current_fig['layout'] = gridplot(current_fig['all_axes'])

    script, div = components(current_fig['layout'])

    if summary is None:
        summary = type(data).__name__

    ls = [formatting.html_svg(div, summary)]
    return formatting.html_from_sections(data, [script] + ls)



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
