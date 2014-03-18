#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
from sherpa.utils import get_keyword_defaults

__all__ = ('clear_window', 'plot', 'contour', 'point', 'set_subplot',
           'get_split_plot_defaults','get_confid_point_defaults',
           'get_plot_defaults','get_point_defaults', 'begin', 'end', 'init',
           'get_data_plot_defaults', 'get_model_plot_defaults',
           'get_fit_plot_defaults', 'get_resid_plot_defaults',
           'get_ratio_plot_defaults', 'get_contour_defaults',
           'get_data_contour_defaults', 'get_model_contour_defaults',
           'get_fit_contour_defaults', 'get_resid_contour_defaults',
           'get_ratio_contour_defaults','get_confid_plot_defaults',
           'get_confid_contour_defaults', 'set_window_redraw', 'set_jointplot',
           'get_histo_defaults', 'get_model_histo_defaults',
           'get_component_plot_defaults', 'get_component_histo_defaults',
           'vline', 'hline', 'get_scatter_plot_defaults', 'get_cdf_plot_defaults')

def point(*args, **kwargs):
    pass

clear_window = point
set_window_redraw = point
end = point
begin = point
init = point
exceptions = point

plot = point
histo = point
contour = point
set_subplot = point
set_jointplot = point

vline = point
hline = point

def get_split_plot_defaults():
    return get_keyword_defaults(set_subplot, 3)


def get_plot_defaults():
    return get_keyword_defaults(plot, 7)


def get_point_defaults():
    return get_keyword_defaults(point, 2)

def get_contour_defaults():
    return get_keyword_defaults(contour, 6)

def get_histo_defaults():
    return get_keyword_defaults(histo, 6)

def get_dummy_defaults():
    return {}

get_data_plot_defaults = get_dummy_defaults
get_model_plot_defaults = get_dummy_defaults
get_fit_plot_defaults = get_dummy_defaults
get_resid_plot_defaults = get_dummy_defaults
get_ratio_plot_defaults = get_dummy_defaults

get_data_contour_defaults = get_dummy_defaults
get_model_contour_defaults = get_dummy_defaults
get_fit_contour_defaults = get_dummy_defaults
get_resid_contour_defaults = get_dummy_defaults
get_ratio_contour_defaults = get_dummy_defaults

get_confid_point_defaults = get_dummy_defaults
get_confid_plot_defaults = get_dummy_defaults
get_confid_contour_defaults = get_dummy_defaults
get_model_histo_defaults = get_dummy_defaults
get_component_plot_defaults = get_dummy_defaults
get_component_histo_defaults = get_dummy_defaults
get_scatter_plot_defaults = get_dummy_defaults
get_cdf_plot_defaults = get_dummy_defaults
