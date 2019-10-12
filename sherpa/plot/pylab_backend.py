#
#  Copyright (C) 2010, 2015, 2017, 2019  Smithsonian Astrophysical Observatory
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

# Although this is labelled pylab mode, use the pyplot interface
# (for the functionlity used here, they are the same).
#
import numpy
from matplotlib import pyplot as plt

from sherpa.utils import get_keyword_defaults
from sherpa.utils.err import NotImplementedErr


__all__ = ('clear_window','point','plot','histo','contour','set_subplot','init',
           'get_split_plot_defaults', 'get_plot_defaults', 'begin', 'end',
           'get_data_plot_defaults', 'get_model_plot_defaults', 'exceptions',
           'get_fit_plot_defaults', 'get_resid_plot_defaults',
           'get_ratio_plot_defaults', 'get_contour_defaults',
           'get_data_contour_defaults', 'get_model_contour_defaults',
           'get_fit_contour_defaults', 'get_resid_contour_defaults',
           'get_ratio_contour_defaults','get_confid_plot_defaults',
           'get_confid_contour_defaults', 'set_window_redraw', 'set_jointplot',
           'get_model_histo_defaults', 'get_histo_defaults',
           'get_component_plot_defaults','get_component_histo_defaults',
           'vline', 'hline', 'get_scatter_plot_defaults', 'get_cdf_plot_defaults',
           'get_latex_for_string', 'name')

name = 'pylab'


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


def point(x, y, overplot=True, clearwindow=False,
          symbol=None,
          color=None):

    axes = setup_axes(overplot, clearwindow)

    if color is None:
        style = '{}'.format(symbol)
    else:
        style = '{}{}'.format(color, symbol)

    axes.plot(numpy.array([x]), numpy.array([y]), style)


def histo(xlo, xhi, y, yerr=None, title=None, xlabel=None, ylabel=None,
          overplot=False, clearwindow=True,
          yerrorbars=False,
          ecolor=_errorbar_defaults['ecolor'],
          capsize=_errorbar_defaults['capsize'],
          barsabove=_errorbar_defaults['barsabove'],
          xlog=False,
          ylog=False,
          linestyle='solid',
          linecolor=None,
          drawstyle='steps-mid',
          color=None,
          marker='None',
          markerfacecolor=None,
          markersize=None):

    xmid = 0.5 * (xlo + xhi)
    plot(xmid, y, yerr=yerr, xerr=None,
         title=title, xlabel=xlabel, ylabel=ylabel,
         overplot=overplot, clearwindow=clearwindow,
         xerrorbars=False, yerrorbars=yerrorbars,
         ecolor=ecolor, capsize=capsize, barsabove=barsabove,
         xlog=xlog, ylog=ylog,
         linestyle=linestyle, linecolor=linecolor,
         drawstyle=drawstyle,
         color=color, marker=marker,
         markerfacecolor=markerfacecolor, markersize=markersize,
         xaxis=False, ratioline=False)


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


def vline(x, ymin=0, ymax=1,
          linecolor=None,
          linestyle=None,
          linewidth=None,
          overplot=False, clearwindow=True):

    axes = setup_axes(overplot, clearwindow)

    line = axes.axvline(x, ymin, ymax)
    _set_line(line, linecolor=linecolor, linestyle=linestyle,
              linewidth=linewidth)


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
         linecolor=None,
         drawstyle='steps',
         color=None,
         marker='None',
         markerfacecolor=None,
         markersize=None,
         xaxis=False,
         ratioline=False):

    axes = setup_axes(overplot, clearwindow)

    # Set up the axes
    if not overplot:
        setup_plot(axes, title, xlabel, ylabel, xlog=xlog, ylog=ylog)

    # The following is from https://github.com/sherpa/sherpa/issues/662
    # which is circa matplotlib 3. The issue is how to ensure that
    # a plot with multiple data values (e.g. a fit plot with
    # data+error bars and model, and one with multiple fits),
    # are drawn "sensibly" so that the user can make out the model.
    # For this we really want the zorder of plot items to be determined
    # by the order they are added. However, this appears to be complicated
    # by the fact that the errorbar command creates multiple objects
    # (for the point, error bars, and I believe caps), and not just
    # a Line2D object. Once all the plot items are concatenated and
    # sorted to ensure a consistent ordering - as discussed in
    # https://github.com/matplotlib/matplotlib/issues/1622#issuecomment-142469866 -
    # we lose a "nice" order in Sherpa (for our cases it does not
    # seem to help that the error bar adds in a zorder Line2D of 2.1
    # as compared to 2, which means that the point for the error bar
    # is drawn on top of the error bar, but also above other lines
    # added to the plot.
    #

    # One option would be to evaluate the current plot to find out
    # what the minimum zorder is, and explicitly set larger
    # values. Note that the default zindex values appear to be 2 and
    # 2.1 from plot/errorbar. This should work for the case of
    #   plot something
    #   overplot something
    # but may fall appart as soon as users add their own features
    # to the visualization, but that may be acceptable.
    #
    zorder = None
    if overplot:
        # Is this a sufficient set of objects?
        #
        objects = axes.lines + axes.collections

        # could do a list comprehension, but want to catch
        # AtributeErrors
        zs = []
        for o in objects:
            try:
                zo = o.get_zorder()
            except AttributeError:
                continue

            zs.append(zo)

        if len(zs) > 0:
            # Add a small offset for this set of data
            zorder = max(zs) + 0.1

    # Rely on color-cycling to work for both the "no errorbar" and
    # "errorbar" case.
    #
    if xerrorbars or yerrorbars:
        if markerfacecolor is None:
            markerfacecolor = color
        if xerr is not None:
            xerr = xerr / 2.
        xerr = xerr if xerrorbars else None
        yerr = yerr if yerrorbars else None
        axes.errorbar(x, y, yerr, xerr,
                      color=color,
                      linestyle=linestyle,
                      drawstyle=drawstyle,
                      marker=marker,
                      markersize=markersize,
                      markerfacecolor=markerfacecolor,
                      ecolor=ecolor,
                      capsize=capsize,
                      barsabove=barsabove,
                      zorder=zorder)

    else:
        axes.plot(x, y,
                  color=color,
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


def contour(x0, x1, y, levels=None, title=None, xlabel=None, ylabel=None,
            overcontour=False, clearwindow=True,
            xlog=False,
            ylog=False,
            linewidths=None,
            colors=None):

    x0 = numpy.unique(x0)
    x1 = numpy.unique(x1)
    y  = numpy.asarray(y)

    if x0.size * x1.size != y.size:
        raise NotImplementedErr('contourgrids')

    y = y.reshape(x1.size, x0.size)

    axes = setup_axes(overcontour, clearwindow)

    # Set up the axes
    if not overcontour:
        setup_plot(axes, title, xlabel, ylabel, xlog=xlog, ylog=ylog)

    if levels is None:
        axes.contour(x0, x1, y, colors=colors, linewidths=linewidths)
    else:
        axes.contour(x0, x1, y, levels, colors=colors,
                     linewidths=linewidths)


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


def set_jointplot(row, col, nrows, ncols, clearaxes=True,
                  top=1,
                  ratio=2):

    if clearaxes:
        # need to set axes[1] as current axes.
        axes = plt.gca()
        ax2 = axes.figure.axes[-1]
        plt.sca(ax2)

    else:
        # follow the chips backend and set plot number "top" (numbering
        # from 1) as the plot with a height of ratio, and the rest
        # with a value of 1. Note that the chips backend created
        # nrows * ncols plots, but here we only create nrows plots.
        #
        ratios = [1] * nrows
        ratios[top - 1] = ratio
        gs = {'height_ratios': ratios}
        f, axarr = plt.subplots(nrows, sharex=True, num=1,
                                gridspec_kw=gs)
        f.subplots_adjust(hspace=0.05)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]],
                 visible=False)

        # need to set axes[0] as current axes.
        plt.sca(axarr[0])


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
    d['linecolor'] = 'red'
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
