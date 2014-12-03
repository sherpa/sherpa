#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
from itertools import izip
import numpy
import pychips as chips
from sherpa.utils import get_keyword_defaults
from sherpa import get_config
from ConfigParser import ConfigParser, NoSectionError

config = ConfigParser()
config.read(get_config())

__all__ = ('clear_window', 'plot', 'histo', 'contour', 'point', 'set_subplot',
           'get_split_plot_defaults','get_confid_point_defaults',
           'get_plot_defaults','get_point_defaults', 'begin', 'end',
           'get_data_plot_defaults', 'get_model_plot_defaults',
           'get_fit_plot_defaults', 'get_resid_plot_defaults',
           'get_ratio_plot_defaults', 'get_contour_defaults', 'exceptions',
           'get_data_contour_defaults', 'get_model_contour_defaults',
           'get_fit_contour_defaults', 'get_resid_contour_defaults',
           'get_ratio_contour_defaults','get_confid_plot_defaults',
           'get_confid_contour_defaults', 'set_window_redraw', 'set_jointplot',
           'get_histo_defaults', 'get_model_histo_defaults',
           'get_component_plot_defaults', 'get_component_histo_defaults',
           'vline', 'hline', 'get_cdf_plot_defaults', 'get_scatter_plot_defaults')

_initialized = False # Set this True first time begin() is called

def _chips_wrap(func, *args, **kwargs):
    func(*args, **kwargs)

def _set_window_redraw(redraw):
    if chips.info_current() is not None:
        if chips.info_current().find('Window') != -1:
            chips.advanced.set_window_redraw(redraw)


def _clear_window():
    if chips.info_current() is not None:
        if chips.info_current().find('Frame') != -1:
            chips.erase()
        else:
            chips.add_frame()
    
    if chips.info() is None:
        chips.add_window()


def _point(x, y, overplot=True, clearwindow=False,
          style=chips.chips_plus,
          color=None,
          size=None,
          angle=None,
          fill=None):

    if (not overplot) and clearwindow:
        _clear_window()
    
    chips.add_point(x, y)

    for var in ('style','color', 'size', 'angle', 'fill'):
        val = locals()[var]
        if val is not None:
            if 'color' in var:
                val = _check_hex_color(val)
            getattr(chips.advanced, 'set_point_' + var)(val)


_attr_map = {
    'linecolor' : 'color',
    'linestyle' : 'style',
    'linewidth' : 'thickness',
    }

_linestyle_map = {

    'noline'  : chips.chips_noline,
    'solid'   : chips.chips_solid,
    'dot'     : chips.chips_shortdash,
    'dash'    : chips.chips_longdash,
    'dotdash' : chips.chips_dotlongdash,
    }


def _check_hex_color(val):
    if type(val) in (str, numpy.string_) and val.startswith('0x'):
        val = str(val).replace('0x','').rjust(6,'0')
    return val

def _vline(x, ymin=0, ymax=1,
          linecolor=None,
          linestyle=None,
          linewidth=None,
          overplot=False, clearwindow=True):

    if (not overplot) and clearwindow:
        _clear_window()

    chips.add_vline(x)

    for var in ('linecolor', 'linestyle', 'linewidth'):
        val = locals()[var]
        if val is not None:
            if 'color' in var:
                val = _check_hex_color(val)
            elif 'style' in var:
                val = _linestyle_map[val]
            getattr(chips.advanced, 'set_line_' + _attr_map[var])(val)


def _hline(y, xmin=0, xmax=1,
          linecolor=None,
          linestyle=None,
          linewidth=None,
          overplot=False, clearwindow=True):

    if (not overplot) and clearwindow:
        _clear_window()

    chips.add_hline(y)

    for var in ('linecolor', 'linestyle', 'linewidth'):
        val = locals()[var]
        if val is not None:
            if 'color' in var:
                val = _check_hex_color(val)
            elif 'style' in var:
                val = _linestyle_map[val]
            getattr(chips.advanced, 'set_line_' + _attr_map[var])(val)




def _plot(x, y, yerr=None, xerr=None, title=None, xlabel=None, ylabel=None,
         overplot=False, clearwindow=True,
         xerrorbars=False,
         yerrorbars=False,
         errstyle=None,
         errcolor=None,
         errthickness=None,
         xlog=False,
         ylog=False,
         linestyle=chips.chips_solid,
         linecolor=None,
         symbolstyle=chips.chips_none,
         symbolcolor=None,
         symbolsize=None,
         symbolfill=True,
         linethickness=None,
         xaxis=False,
         ratioline=False):


    if (not overplot) and clearwindow:
        _clear_window()

    if yerrorbars and (yerr is not None) and xerrorbars and (xerr is not None):
        xerr = xerr / 2.
        chips.add_curve(x, y, (yerr, yerr, xerr, xerr) )
    elif yerrorbars and (yerr is not None):
        chips.add_curve(x, y, yerr)
    else:
        chips.add_curve(x, y)

    for var in ('errstyle', 'errcolor', 'errthickness', 'linestyle',
                'linecolor', 'symbolstyle', 'symbolcolor', 'symbolsize',
                'symbolfill', 'linethickness'):
        val = locals()[var]
        if val is not None:
            if 'color' in var:
                val = _check_hex_color(val)
            getattr(chips.advanced, 'set_curve_' + var)(val)

    if not overplot:
        for log_axis, axis_id in izip((xlog, ylog),
                                      (chips.X_AXIS, chips.Y_AXIS)):
            if log_axis:
                chips.log_scale(axis_id)
            else:
                chips.linear_scale(axis_id)

        if title:
            ttl = title.replace('_', '\\_')
            chips.set_plot_title(ttl)
        if xlabel:
            xlbl = xlabel.replace('_', '\\_')
            chips.set_plot_xlabel(xlbl)
        if ylabel:
            ylbl = ylabel.replace('_', '\\_')
            chips.set_plot_ylabel(ylbl)

    if xaxis:
        chips.add_hline(0);

    if ratioline:
        chips.add_hline(1);

    #chips.limits(chips.X_AXIS, 'AUTO', 'AUTO')


def _histogram(xlo, xhi, y, yerr=None, title=None, xlabel=None, ylabel=None,
               overplot=False, clearwindow=True,
               yerrorbars=False,
               errstyle=None,
               errcolor=None,
               errthickness=None,
               fillcolor=None,
               fillopacity=None,
               fillstyle=None,
               xlog=False,
               ylog=False,
               linestyle=chips.chips_solid,
               linecolor=None,
               linethickness=None,
               symbolangle=None,
               symbolcolor=None,
               symbolfill=None,
               symbolsize=None,
               symbolstyle=chips.chips_none):

    if (not overplot) and clearwindow:
        _clear_window()

    if yerrorbars and yerr is not None:
        chips.add_histogram(xlo, xhi, y, yerr)
    else:
        chips.add_histogram(xlo, xhi, y)

    for var in ('errstyle', 'errcolor', 'errthickness',
                'fillcolor', 'fillopacity', 'fillstyle',
                'linestyle', 'linecolor', 'linethickness',
                'symbolangle', 'symbolcolor', 'symbolfill', 'symbolsize',
                'symbolstyle'):
        val = locals()[var]
        if val is not None:
            if 'color' in var:
                val = _check_hex_color(val)
            getattr(chips.advanced, 'set_histogram_' + var)(val)

    if not overplot:
        for log_axis, axis_id in izip((xlog, ylog),
                                      (chips.X_AXIS, chips.Y_AXIS)):
            if log_axis:
                chips.log_scale(axis_id)
            else:
                chips.linear_scale(axis_id)

        if title:
            ttl = title.replace('_', '\\_')
            chips.set_plot_title(ttl)
        if xlabel:
            xlbl = xlabel.replace('_', '\\_')
            chips.set_plot_xlabel(xlbl)
        if ylabel:
            ylbl = ylabel.replace('_', '\\_')
            chips.set_plot_ylabel(ylbl)

    #chips.limits(chips.X_AXIS, 'AUTO', 'AUTO')


def _contour(x0, x1, y, levels=None, title=None, xlabel=None, ylabel=None,
             overcontour=False, clearwindow=True,
             xlog=False,
             ylog=False,
             style=None,
             color=None,
             thickness=None,
             axis_pad=0.05):

    if (not overcontour) and clearwindow:
        _clear_window()

    # Catch NANs before sending to ChIPS
    bad = list(numpy.where(numpy.isnan(y)==True)).pop(0)
    bad_vals = numpy.array(y[bad])
    y[bad] = 0.0

    if levels is None:
        chips.add_contour(x0, x1, y)
    else:
        levels = numpy.asarray(levels, numpy.float_)
        chips.add_contour(x0, x1, y, levels)

    y[bad] = bad_vals

    for var in ('style', 'color', 'thickness'):
        val = locals()[var]
        if val is not None:
            if 'color' in var:
                val = _check_hex_color(val)
            getattr(chips.advanced, 'set_contour_' + var)(val)

    chips.advanced.set_axis_pad(axis_pad)

    chips.set_data_aspect_ratio()
    chips.limits(chips.X_AXIS, x0.min(), x0.max())
    chips.limits(chips.Y_AXIS, x1.min(), x1.max())

    if not overcontour:
        for log_axis, axis_id in izip((xlog, ylog),
                                      (chips.X_AXIS, chips.Y_AXIS)):
            if log_axis:
                chips.log_scale(axis_id)
            else:
                chips.linear_scale(axis_id)

        if title:
            ttl = title.replace('_', '\\_')
            chips.set_plot_title(ttl)
        if xlabel:
            xlbl = xlabel.replace('_', '\\_')
            chips.set_plot_xlabel(xlbl)
        if ylabel:
            ylbl = ylabel.replace('_', '\\_')
            chips.set_plot_ylabel(ylbl)


def _set_subplot(row, col, nrows, ncols, clearaxes=True,
                xgap=0.18,
                ygap=0.18):
    
    chips.add_plot()
    chips.grid_objects(ncols, nrows, xgap, ygap)


def _set_jointplot(row, col, nrows, ncols, clearaxes=True,
                  top=1,
                  ratio=2):
    
    # FIXME: misuse of kwarg clearaxes
    if not clearaxes:
        chips.strip_chart(nrows*ncols)
        chips.adjust_grid_yrelsize(top,ratio)
    else:
        chips.set_current_plot('plot2')

def init():
    # This function now a no-op; structure kept in case
    # we want init to do something in future.
    pass

def begin():
    global _initialized

    chips.lock()    
    chips.advanced.open_undo_buffer()
    if _initialized is False:
        try:
            overrides = config.items('chips')
            for item in overrides:
                chips.set_preference(item[0], item[1])
            # OL: No apparent reason to call add_window() here.
            # ChIPS is smart enough to open a window if none are available,
            # plus this code only gets executed if the user has a [chips] section in sherpa.rc
            # which is the exception rather than the rule.
            # chips.add_window() # Have Sherpa talk to its own
                               # chips window
        except NoSectionError:
            chips.unlock()
        except:
            chips.unlock()
            raise
        _initialized = True

def end():
    # Don't need to call redraw here ourselves, the
    # ChIPS undo buffer does what we need.
    chips.advanced.close_undo_buffer()
    chips.unlock()

def exceptions():
    chips.advanced.discard_undo_buffer()
    chips.erase()
    chips.unlock()

def clear_window(*args, **kwargs):
    _chips_wrap( _clear_window, *args, **kwargs)


def set_window_redraw(*args, **kwargs):
    _chips_wrap( _set_window_redraw, *args, **kwargs)


def point(*args, **kwargs):
    _chips_wrap(_point, *args, **kwargs)

def vline(*args, **kwargs):
    _chips_wrap(_vline, *args, **kwargs)

def hline(*args, **kwargs):
    _chips_wrap(_hline, *args, **kwargs)

def plot(*args, **kwargs):
    _chips_wrap(_plot, *args, **kwargs)


def histo(*args, **kwargs):
    _chips_wrap(_histogram, *args, **kwargs)


def contour(*args, **kwargs):
    _chips_wrap(_contour, *args, **kwargs)


def set_subplot(*args, **kwargs):
    _chips_wrap(_set_subplot, *args, **kwargs)


def set_jointplot(*args, **kwargs):
    _chips_wrap(_set_jointplot, *args, **kwargs)


def get_split_plot_defaults():
    return get_keyword_defaults(_set_subplot, 3)


def get_plot_defaults():
    return get_keyword_defaults(_plot, 7)


def get_point_defaults():
    return get_keyword_defaults(_point, 2)


def get_histo_defaults():
    return get_keyword_defaults(_histogram, 6)


def get_confid_point_defaults():
    d = get_point_defaults()
    d['style'] = chips.chips_plus
    d['size'] = 7
    return d


def get_data_plot_defaults():
    d = get_plot_defaults()

    d['yerrorbars'] = True
    d['errstyle'] = 'line'
    d['linestyle'] = chips.chips_noline
    d['symbolstyle'] = chips.chips_circle
    d['symbolsize'] = 3
    d['symbolfill'] = False

    return d


def get_model_plot_defaults():
    d = get_plot_defaults()

    d['linestyle'] = chips.chips_solid
    d['linethickness'] = 3
    d['linecolor'] = 'red'
    d['symbolstyle'] = chips.chips_none

    return d


def get_confid_plot_defaults():
    d = get_plot_defaults()

    d['linestyle'] = chips.chips_solid
    d['symbolstyle'] = chips.chips_none
    d['linethickness'] = 3
    d['linecolor'] = 'skyblue'
    return d


def get_fit_plot_defaults():
    return {}


def get_resid_plot_defaults():
    d = get_data_plot_defaults()
    d['xerrorbars'] = True
    d['errstyle'] = 'line'
    d['symbolstyle'] = chips.chips_diamond
    d['xaxis'] = True
    d['symbolsize'] = 3
    return d


def get_ratio_plot_defaults():
    d = get_data_plot_defaults()
    d['xerrorbars'] = True
    d['errstyle'] = 'line'
    d['symbolstyle'] = chips.chips_diamond
    d['ratioline'] = True
    d['symbolsize'] = 3
    return d


def get_contour_defaults():
    return get_keyword_defaults(contour, 6)

get_data_contour_defaults = get_contour_defaults

def get_model_contour_defaults():
    d = get_contour_defaults()
    d['style'] = None
    d['color'] = 'red'
    d['thickness'] = 3
    return d

def get_confid_contour_defaults():
    d = get_contour_defaults()
    d['style'] = None
    d['color'] = 'skyblue'
    d['thickness'] = 3
    d['axis_pad'] = 0.0
    return d


def get_fit_contour_defaults():
    return {}


get_resid_contour_defaults = get_data_contour_defaults
get_ratio_contour_defaults = get_data_contour_defaults

def get_model_histo_defaults():
    d = get_histo_defaults()
#    d['linestyle'] = chips.chips_solid
    d['linethickness'] = 2
    d['linecolor'] = 'red'

    return d

def get_component_plot_defaults():
    d = get_model_plot_defaults()
    d['linecolor'] = 'orange'
    return d

def get_component_histo_defaults():
    d = get_model_histo_defaults()
    d['linecolor'] = 'orange'
    return d

def get_cdf_plot_defaults():
    d = get_model_plot_defaults()
    d['linecolor'] = 'red'
    return d

def get_scatter_plot_defaults():
    d = get_data_plot_defaults()
    d['symbolsize'] = 1
    d['symbolfill'] = True
    d['yerrorbars'] = False
    return d
