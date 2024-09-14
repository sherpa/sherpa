#
#  Copyright (C) 2009, 2015, 2016, 2018 - 2024
#  Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

"""A visualization interface to Sherpa.

Classes provide access to common plotting tasks, which is done by the
plotting backend defined in the ``options.plot_pkg`` setting of the
Sherpa configuration file. Note that plot objects can be created and
used even when only the `sherpa.plot.backends.BasicBackend` is
available, it is just that no graphical display will be created.

Which backend is used?
----------------------

When this module is first imported, Sherpa tries to import the
backends installed with Sherpa in the order listed in the
``options.plot_pkg`` setting from the ``sherpa.rc`` startup file.
The first module that imports successfully is set as the active
backend. The following command prints the name of the backend:

   >>> from sherpa import plot
   >>> print(plot.backend.name)
   pylab

Change the backend
------------------

After the initial import, the backend can be changed by loading one of
the plotting backends shipped with sherpa (or any other module that
provides the same interface):

  >>> from sherpa.plot.pylab_backend import PylabBackend
  >>> plot.backend = PylabBackend()

Creating a plot
---------------

A plot backend can act as a context manager to apply a specific backend to just one plot
without globally changing the backend for the rest of the session::

  >>> from sherpa.plot.pylab_backend import PylabBackend
  >>> with PylabBackend():
  ...     # Now call the plot/overplot or contour/overcontour methods
  ...     obj.plot()   # doctest: +SKIP

This handles setting up the backend, handles any error handling,
and then ends the session.

"""

from __future__ import annotations

from configparser import ConfigParser
import contextlib
import copy
import logging
import importlib
from typing import Any, Literal, Optional, Sequence, Union, \
    overload

import numpy as np

from sherpa import get_config
from sherpa.data import Data, Data1D, Data1DInt, Data2D
from sherpa.estmethods import Covariance
from sherpa.models.model import Model
from sherpa.optmethods import LevMar, NelderMead
from sherpa.plot.backends import BaseBackend, BasicBackend, PLOT_BACKENDS
from sherpa.stats import Stat, Likelihood, LeastSq, Chi2XspecVar
from sherpa.utils import NoNewAttributesAfterInit, erf, \
    bool_cast, parallel_map, dataspace1d, histogram1d, get_error_estimates, \
    is_iterable_not_str
from sherpa.utils.err import ArgumentTypeErr, ConfidenceErr, \
    IdentifierErr, PlotErr, StatErr
from sherpa.utils.numeric_types import SherpaFloat

# PLOT_BACKENDS only contains backends in modules that are imported successfully
# but modules are not discovered by itself. Entrypoints would solve this problem
# but the current implementation does not have this capability.
# See docstring of sherpa.plot.backends.MetaBaseBackend for details.
#
for name in ["pylab", "pylab_area", "bokeh"]:
    try:
        importlib.import_module(f"sherpa.plot.{name}_backend")
    except ImportError:
        pass

config = ConfigParser()
config.read(get_config())

warning = logging.getLogger(__name__).warning

# TODO: why is this module globally changing the invalid mode of NumPy?
_ = np.seterr(invalid='ignore')

basicbackend = BasicBackend()
'''This backend has all the defaults that are backend-independent.
For all the plot classes in this module (e.g. FitPlot, JointPlot, ...)
the default settings are set in code like `plot_prefs = xxx.get_plot_defaults()`
where xxx is some backend.
'''

# This may be over-ridden below, depending on the options.plot_pkg setting.
#
backend: BaseBackend = basicbackend
'''Currently active backend module for plotting.'''

# In the current design, this code is executed when this
# module imported and then the values are shared between all instances of a class.
# That allows one to change those default "globally" (i.e. for all instances of the class)
# but it means that the default values are based on the backend that is active when
# this module is imported. That depends on the sherpa.rc setting and thus `plot_prefs`
# might contain e.g. matplotlib specific defaults that are not applicable when other
# backends are used.
# Thus, we currently initialize them with the BasicBackend that has only the
# backend-independent (i.e. those that works for any backend) defaults set.


plot_opt_str = config.get('options', 'plot_pkg', fallback='BasicBackend')
plot_opt = [o.strip() for o in plot_opt_str.split()]

for plottry in plot_opt:
    if plottry in PLOT_BACKENDS:
        backend = PLOT_BACKENDS[plottry]()
        break
    warning("Plotting backend '%s' not found or dependencies missing. Trying next option.", plottry)
else:
    # The backend will default to BasicBackend from above
    warning(f"Tried the following backends listed in {get_config()}: \n" +
            f"{plot_opt}\n" +
            "None of these imported correctly, so using the 'BasicBackend'.\n" +
            f"List of backends that have loaded and would be available: {[k for k in PLOT_BACKENDS]}")

__all__ = ('Plot', 'Contour', 'Point', 'Histogram',
           'MultiPlot',
           'HistogramPlot', 'DataHistogramPlot',
           'ModelHistogramPlot', 'SourceHistogramPlot',
           'PDFPlot', 'CDFPlot', 'LRHistogram',
           'SplitPlot', 'JointPlot',
           'DataPlot', 'TracePlot', 'ScatterPlot',
           'PSFKernelPlot',
           'DataContour', 'PSFKernelContour',
           'ModelPlot', 'ComponentModelPlot',
           'ComponentModelHistogramPlot', 'ComponentTemplateModelPlot',
           'SourcePlot', 'ComponentSourcePlot',
           'ComponentSourceHistogramPlot', 'ComponentTemplateSourcePlot',
           'PSFPlot',
           'ModelContour', 'PSFContour', 'SourceContour',
           'FitPlot',
           'FitContour',
           'DelchiPlot', 'DelchiHistogramPlot',
           'ChisqrPlot', 'ChisqrHistogramPlot',
           'ResidPlot', 'ResidHistogramPlot',
           'RatioPlot', 'RatioHistogramPlot',
           'ResidContour',
           'RatioContour',
           'Confidence1D', 'Confidence2D',
           'IntervalProjection', 'IntervalUncertainty',
           'RegionProjection', 'RegionUncertainty',
           'TemporaryPlottingBackend',
           'set_backend',
           )

_stats_noerr = ('cash', 'cstat', 'leastsq', 'wstat')
"""Statistics that do not make use of uncertainties."""


def set_backend(new_backend: Union[str, BaseBackend, type[BaseBackend]]) -> None:
    '''Set the Sherpa plotting backend.

    Plotting backends are registered in Sherpa with a string name.
    See the examples below for how to get a list of available backends.

    Parameters
    ----------
    new_backend : string, class, or instance
        Set a sherpa plotting backend. The backend can be passed in as an
        object, or as a convenience, as a simple string

    Example
    -------
    Set the backend to use Pylab from matplotlib for plotting. This is
    probably what most users need:

        >>> from sherpa.plot import set_backend
        >>> set_backend('pylab')

    Get a list of registered backends:

        >>> from sherpa.plot import PLOT_BACKENDS
        >>> PLOT_BACKENDS
        {'BaseBackend': <class 'sherpa.plot.backends.BaseBackend'>, ...}

    This list shows the names and the class for each backend.
    Details for each backend can be found in the Sherpa documentation or by
    inspecting the backend classes using normal Python functionality:

        >>> from sherpa.plot.backends import BasicBackend
        >>> help(BasicBackend)
        Help on class BasicBackend in module sherpa.plot.backends:
        <BLANKLINE>
        class BasicBackend(BaseBackend)
        |  A dummy backend for plotting.
        ...
    '''
    global backend

    if isinstance(new_backend, str):
        try:
            backend = PLOT_BACKENDS[new_backend]()
        except KeyError:
            raise IdentifierErr('noplotbackend', new_backend,
                                list(PLOT_BACKENDS.keys())) from None

        return

    # It's a class and that class is a subclass of BaseBackend
    if isinstance(new_backend, type) and issubclass(new_backend, BaseBackend):
        backend = new_backend()
        return

    if isinstance(new_backend, BaseBackend):
        backend = new_backend
        return

    raise ArgumentTypeErr('tempplotbackend', new_backend)


class TemporaryPlottingBackend(contextlib.AbstractContextManager):
    '''Set the Sherpa plotting backend as a context, e.g. for a single plot

    This changes the logging level globally for all modules in sherpa.

    Parameters
    ----------
    new_backend : string, class, or instance
        Set a sherpa plotting backend. The backend can be passed in as an
        instance of a plotting backend class. For simplicity, the user can
        also pass in a string naming a loaded backend class or the class
        itself; calling this context manager will then create an instance.

    Example
    -------

    >>> from sherpa.plot import TemporaryPlottingBackend, DataPlot
    >>> from sherpa.data import Data1D
    >>> with TemporaryPlottingBackend('pylab'):
    ...     x1 = [100, 200, 600, 1200]
    ...     y1 = [2000, 2100, 1400, 3050]
    ...     d1 = Data1D('oned', x1, y1)
    ...     plot1 = DataPlot()
    ...     plot1.prepare(d1)
    ...     plot1.plot()

    '''

    def __init__(self, new_backend: BaseBackend) -> None:
        self.backend = new_backend

    def __enter__(self) -> TemporaryPlottingBackend:
        global backend
        self.old = backend
        set_backend(self.backend)
        return self

    # As exc_type/val/tb are not used here we type them as Any, and the
    # return value is marked as an explicit False rather than bool since
    # it avoids warnings from mypy (as __exit__ is special).
    #
    def __exit__(self,
                 exc_type: Any, exc_val: Any, exc_tb: Any) -> Literal[False]:
        global backend
        backend = self.old
        return False


def _make_title(title: str, name: Optional[str] = '') -> str:
    """Return the plot title to use.

    Parameters
    ----------
    title : str
        The main title to use
    name : str or None
        The identifier for the dataset.

    Returns
    -------
    title : str
        If the name is empty or None then use title,
        otherwise use title + ' for ' + name.

    """

    if name in [None, '']:
        return title

    return f"{title} for {name}"


def _errorbar_warning(stat: Stat) -> str:
    """The warning message to display when error bars are being "faked".

    Parameters
    ----------
    stat : sherpa.stats.Stat instance
        The name attribute is used in the error message.

    Returns
    -------
    msg : str
        The warning message
    """

    return "The displayed errorbars have been supplied with the " + \
        "data or calculated using chi2xspecvar; the errors are not " + \
        f"used in fits with {stat.name}"


def calculate_errors(data: Data,
                     stat: Stat,
                     yerrorbars: bool = True
                     ) -> Optional[np.ndarray]:
    """Calculate errors from the statistics object."""

    if stat is None:
        return None

    # Do we warn about adding in error values? The logic here isn't
    # quite the same as the other times _errorbar_warning is used.
    # This suggests that the code really should be re-worked to
    # go through a common code path.
    #
    # Do we require that self.yerr is always set, even if the
    # value is not used in a plot? I am going to assume not
    # (since the existing code will not change self.yerr if
    # the stat.name is not in _stats_noerr but an exception is
    # raised by get_yerr).
    #
    # This assumes that the error bars calculated by data.to_plot are
    # over-ridden once a statistic is given. However, how does this
    # work if the statistic is a Chi2 variant - e.g. Chi2DataVar -
    # but the user has given explicit errors (i.e. they are not to
    # be calculated by the "DataVar" part but used as is). Does
    # data.get_yerr handle this for us, or are invalid errors
    # used here? It appears that the correct answers are being
    # returned, but should we only call data.get_yerr if yerr
    # is None/empty/whatever is returned by to_plot? This also
    # holds for the Resid/RatioPlot classes.
    #
    # Note that we should probably return a value for yerr if
    # we can, even if 'yerrorbars' is set to False, so that
    # downstream users can make use of the value even if the
    # plot doesn't. This is similar to how labels or xerr
    # attributes are created even if they don't get used by the
    # plot.
    #

    msg = _errorbar_warning(stat)
    if stat.name in _stats_noerr:
        if yerrorbars:
            warning(msg)

        return data.get_yerr(True, Chi2XspecVar.calc_staterror)

    try:
        return data.get_yerr(True, stat.calc_staterror)
    except Exception:
        # TODO: can we report a useful error here?
        #
        # It is possible that this is actually an unrelated
        # error: it's unclear what error class is expected to
        # be thrown, but likely ValueError as this is raised
        # by Chi2DataVar when sent values < 0
        # (note that over time the behavior has changed from
        #  <= 0 to < 0, but this error message has not been
        # changed).
        #
        if yerrorbars:
            warning("%s\nzeros or negative values found", msg)

        return None


@overload
def arr2str(x: None) -> None:
    ...

@overload
def arr2str(x: ArrayType) -> str:
    ...

def arr2str(x):
    """Convert an array to a string for __str__ calls

    Parameters
    ----------
    x : sequence or None

    Returns
    -------
    arrstr : str

    """

    if x is None:
        return x

    return np.array2string(np.asarray(x), separator=',',
                           precision=4, suppress_small=False)


# This is used in a similar manner to the way that sherpa.data._fields
# is used to control the string output of the Data objects. We extend
# that version to allow a decision of whether to display as an array -
# via arr2str - or as is. Given how our classes are arranged it is not
# worth trying to share infrastructure between the plot classes and
# the data classes to handle the similar __str__ handling.
#
def display_fields(obj,  # hard to provide an accurate type
                   fields: Sequence[str]) -> str:
    """Create the __str__ output for the object.

    Parameters
    ----------
    obj
       The object to convert to a string. It is expected to be
       derived from one of the plot classes.
    fields
       The fields to diplay. A field name ending in ! means display
       directly, otherwise the value is passed through arr2str.

    Returns
    -------
    strval
       The string representation (one field per line)

    Notes
    -----

    The default length for the labels is 6 but we check the labels so
    that a larger value is used if need be (except for any xxx_prefs
    label, just to better match the original output). This is wasted
    logic, in that it could be calculated at object creation, but it
    is not that much work to do here.

    """

    out = []
    size = 6
    for lbl in fields:
        if lbl.endswith('!'):
            lbl = lbl[:-1]
            scalar = True
        else:
            scalar = False

        val = getattr(obj, lbl)
        if not scalar:
            val = arr2str(val)

        out.append((lbl, val))
        if not lbl.endswith("_prefs"):
            size = max(size, len(lbl))

    fmt = f"{{:{size}s}} = {{}}"
    return "\n".join(fmt.format(*vals) for vals in out)


class Plot(NoNewAttributesAfterInit):
    "Base class for line plots"

    plot_prefs = basicbackend.get_plot_defaults()
    "The preferences for the plot."

    def __init__(self):
        """
        Initialize a Plot object.  All 1D line plot
        instances utilize Plot, which provides a generic
        interface to a backend.

        Once an instance of Plot is initialized no new
        attributes of the class can be made. (To eliminate
        the accidental creation of erroneous attributes)
        """
        self.plot_prefs = self.plot_prefs.copy()
        super().__init__()

    @staticmethod
    def vline(x, ymin=0, ymax=1,
              linecolor=None, linestyle=None, linewidth=None,
              overplot=False, clearwindow=True):
        "Draw a line at constant x, extending over the plot."
        backend.vline(x, ymin=ymin, ymax=ymax, linecolor=linecolor,
                      linestyle=linestyle, linewidth=linewidth,
                      overplot=overplot, clearwindow=clearwindow)

    @staticmethod
    def hline(y, xmin=0, xmax=1,
              linecolor=None, linestyle=None, linewidth=None,
              overplot=False, clearwindow=True):
        "Draw a line at constant y, extending over the plot."
        backend.hline(y, xmin=xmin, xmax=xmax, linecolor=linecolor,
                      linestyle=linestyle, linewidth=linewidth,
                      overplot=overplot, clearwindow=clearwindow)

    def _merge_settings(self, kwargs):
        """Return the plot preferences merged with user settings."""
        return {**self.plot_prefs, **kwargs}

    def plot(self, x, y, yerr=None, xerr=None, title=None, xlabel=None,
             ylabel=None, overplot=False, clearwindow=True,
             **kwargs):
        """Plot the data.

        Parameters
        ----------
        x, y
           The data values to plot. They should have the same length.
        yerr, xerr : optional
           The symmetric errors to apply to the y or x values, or `None`
           for no errors along the axis. If given, they should match the
           length of the data.
        title, xlabel, ylabel : optional, string
           The optional plot title and axis labels. These are ignored
           if overplot is set to `True`.
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           plot_prefs dictionary.

        See Also
        --------
        overplot

        """

        opts = self._merge_settings(kwargs)
        backend.plot(x, y, yerr=yerr, xerr=xerr,
                     title=title, xlabel=xlabel, ylabel=ylabel,
                     overplot=overplot, clearwindow=clearwindow,
                     **opts)

    def overplot(self, *args, **kwargs):
        """Add the data to an existing plot.

        This is the same as calling the plot method with
        overplot set to True.

        See Also
        --------
        plot

        """
        kwargs['overplot'] = True
        self.plot(*args, **kwargs)


class Contour(NoNewAttributesAfterInit):
    "Base class for contour plots"

    contour_prefs = basicbackend.get_contour_defaults()
    "The preferences for the plot."

    def __init__(self):
        """
        Initialize a Contour object.  All 2D contour plot
        instances utilize Contour, which provides a generic
        interface to a backend.

        Once an instance of Contour is initialized no new
        attributes of the class can be made. (To eliminate
        the accidental creation of erroneous attributes)
        """
        self.contour_prefs = self.contour_prefs.copy()
        super().__init__()

    def _merge_settings(self, kwargs):
        """Return the plot preferences merged with user settings."""
        return {**self.contour_prefs, **kwargs}

    def contour(self, x0, x1, y, title=None, xlabel=None,
                ylabel=None, overcontour=False, clearwindow=True,
                **kwargs):
        opts = self._merge_settings(kwargs)
        backend.contour(x0, x1, y, title=title,
                        xlabel=xlabel, ylabel=ylabel,
                        overcontour=overcontour,
                        clearwindow=clearwindow, **opts)

    def overcontour(self, *args, **kwargs):
        kwargs['overcontour'] = True
        self.contour(*args, **kwargs)


class Point(NoNewAttributesAfterInit):
    "Base class for point plots"

    point_prefs = basicbackend.get_point_defaults()
    "The preferences for the plot."

    def __init__(self):
        """
        Initialize a Point object.  All 1D point plot
        instances utilize Point, which provides a generic
        interface to a backend.

        Once an instance of Point is initialized no new
        attributes of the class can be made. (To eliminate
        the accidental creation of erroneous attributes)
        """
        self.point_prefs = self.point_prefs.copy()
        super().__init__()

    def _merge_settings(self, kwargs):
        """Return the plot preferences merged with user settings."""
        return {**self.point_prefs, **kwargs}

    def point(self, x, y, overplot=True, clearwindow=False, **kwargs):
        """Draw a point at the given location.

        Parameters
        ----------
        x, y
           The coordinates of the plot.
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           point_prefs dictionary.

        """
        opts = self._merge_settings(kwargs)
        backend.plot(x, y,
                     overplot=overplot, clearwindow=clearwindow,
                     **opts)


class Image(NoNewAttributesAfterInit):
    """Base class for image plots

    .. warning::
        This class is experimental and subject to change in the future.
        Currently, is is only used within _repr_html_ methods.
    """

    image_prefs = basicbackend.get_image_defaults()
    "The preferences for the plot."

    def __init__(self):
        """
        Initialize a Image object.

        Once an instance of Image is initialized no new
        attributes of the class can be made. (To eliminate
        the accidental creation of erroneous attributes.)
        """
        self.image_prefs = self.image_prefs.copy()
        super().__init__()

    def _merge_settings(self, kwargs):
        """Return the plot preferences merged with user settings."""
        return {**self.image_prefs, **kwargs}

    def plot(self, x0, x1, y, overplot=True, clearwindow=False, **kwargs):
        """Draw a point at the given location.

        Parameters
        ----------
        x0, x1 : array-like
            Value for the image axes.
        y : array-like, 2D
           The coordinates of the plot.
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           image_prefs dictionary.

        """
        opts = self._merge_settings(kwargs)
        backend.image(x0, x1, y,
                     overplot=overplot, clearwindow=clearwindow,
                     **opts)

class Histogram(NoNewAttributesAfterInit):
    "Base class for histogram plots"

    histo_prefs = basicbackend.get_histo_defaults()
    "The preferences for the plot."

    def __init__(self):
        """
        Initialize a Histogram object.  All 1D histogram plot
        instances utilize Histogram, which provides a generic
        interface to a backend.

        Once an instance of Histogram is initialized no new
        attributes of the class can be made. (To eliminate
        the accidental creation of erroneous attributes)
        """
        self.histo_prefs = self.histo_prefs.copy()
        super().__init__()

    def _merge_settings(self, kwargs):
        """Return the plot preferences merged with user settings."""
        return {**self.histo_prefs, **kwargs}

    def plot(self, xlo, xhi, y, yerr=None, title=None, xlabel=None,
             ylabel=None, overplot=False, clearwindow=True, **kwargs):
        """Plot the data.

        Parameters
        ----------
        xlo, xhi, y
           The data values to plot (the bin edges along the X axis and
           the bin values for the Y axis). They should have the same length.
        yerr : optional
           The symmetric errors to apply to the y values, or `None`
           for no errors along the axis. If given, they should match the
           length of the data.
        title, xlabel, ylabel : optional, string
           The optional plot title and axis labels. These are ignored
           if overplot is set to `True`.
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           plot_prefs dictionary.

        See Also
        --------
        overplot

        """

        opts = self._merge_settings(kwargs)
        backend.histo(xlo, xhi, y, yerr=yerr, title=title,
                      xlabel=xlabel, ylabel=ylabel, overplot=overplot,
                      clearwindow=clearwindow, **opts)

    def overplot(self, *args, **kwargs):
        kwargs['overplot'] = True
        self.plot(*args, **kwargs)


class HistogramPlot(Histogram):
    """Base class for histogram-style plots with a prepare method."""

    _fields: list[str] = ["xlo", "xhi", "y", "xlabel!", "ylabel!",
                          "title!", "histo_prefs!"]
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    def __init__(self):
        self.xlo = None
        self.xhi = None
        self.y = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        super().__init__()

    def __str__(self) -> str:
        return display_fields(self, self._fields)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the histogram plot."""
        return backend.as_html_histogram(self)

    @property
    def x(self):
        """Return (xlo + xhi) / 2

        This is intended to make it easier to swap between
        plot- and histogram-style plots by providing access
        to an X value.
        """

        if self.xlo is None or self.xhi is None:
            return None

        # As we do not (yet) require NumPy arrays, enforce it.
        xlo = np.asarray(self.xlo)
        xhi = np.asarray(self.xhi)
        return (xlo + xhi) / 2

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        """Plot the data.

        This will plot the data sent to the prepare method.

        Parameters
        ----------
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           plot_prefs dictionary.

        See Also
        --------
        prepare, overplot

        """

        super().plot(self.xlo, self.xhi, self.y, title=self.title,
                     xlabel=self.xlabel, ylabel=self.ylabel,
                     overplot=overplot, clearwindow=clearwindow,
                     **kwargs)

    @staticmethod
    def vline(x, ymin=0, ymax=1,
              linecolor=None, linestyle=None, linewidth=None,
              overplot=False, clearwindow=True):
        "Draw a line at constant x, extending over the plot."
        backend.vline(x, ymin=ymin, ymax=ymax, linecolor=linecolor,
                      linestyle=linestyle, linewidth=linewidth,
                      overplot=overplot, clearwindow=clearwindow)

    @staticmethod
    def hline(y, xmin=0, xmax=1,
              linecolor=None, linestyle=None, linewidth=None,
              overplot=False, clearwindow=True):
        "Draw a line at constant y, extending over the plot."
        backend.hline(y, xmin=xmin, xmax=xmax, linecolor=linecolor,
                      linestyle=linestyle, linewidth=linewidth,
                      overplot=overplot, clearwindow=clearwindow)


# I think we want slightly different histogram preferences
# than most (mark by point rather than line).
#
def get_data_hist_prefs():
    """Copy the data preferences to the histogram class"""

    hprefs = basicbackend.get_model_histo_defaults()
    dprefs = basicbackend.get_data_plot_defaults()
    for k, v in dprefs.items():
        if k in hprefs:
            hprefs[k] = v

    return hprefs


class DataHistogramPlot(HistogramPlot):
    """Create 1D histogram plots of data values."""

    histo_prefs = get_data_hist_prefs()
    "The preferences for the plot."

    def __init__(self):
        self.yerr = None
        super().__init__()

    @property
    def xerr(self):
        """Return abs(xhi - xlow) / 2

        The plotting backends actually calculate the xerr values
        explicitly rather than use this field, so it is provided as a
        property in case users need it.

        .. versionchanged:: 4.16.1
           This field is now a property.

        """

        if self.xlo is None or self.xhi is None:
            return None

        # As we do not (yet) require NumPy arrays, enforce it.
        xlo = np.asarray(self.xlo)
        xhi = np.asarray(self.xhi)
        return np.abs(xhi - xlo) / 2

    def prepare(self, data, stat=None):
        """Create the data to plot

        Parameters
        ----------
        data
            The Sherpa data object to display (it is assumed to be
            one dimensional and represent binned data).
        stat : optional
            The Sherpa statistics object to use to add Y error bars
            if the data has none.

        See Also
        --------
        plot
        """

        # Need a better way of accessing the binning of the data.
        # Maybe to_plot should return the lo/hi edges as a pair
        # here.
        #
        plot = data.to_plot()
        (_, self.y, self.yerr, _, self.xlabel, self.ylabel) = plot

        self.xlo, self.xhi = data.get_indep(True)

        if stat is not None:
            yerrorbars = self.histo_prefs.get('yerrorbars', True)
            self.yerr = calculate_errors(data, stat, yerrorbars)

        self.title = data.name

    # We have
    #
    # - the base class Histogram.plot can accept a yerr argument
    # - the superclass (HistogramPlot.plot) does not know
    #   anything about yerr
    #
    # so the superclass is over-ridden here to basically call
    # the base class.
    #
    def plot(self, overplot=False, clearwindow=True, **kwargs):
        """Plot the data.

        This will plot the data sent to the prepare method.

        Parameters
        ----------
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           plot_prefs dictionary.

        See Also
        --------
        prepare, overplot

        """

        Histogram.plot(self, self.xlo, self.xhi, self.y,
                       yerr=self.yerr, title=self.title,
                       xlabel=self.xlabel, ylabel=self.ylabel,
                       overplot=overplot, clearwindow=clearwindow, **kwargs)


class ModelHistogramPlot(HistogramPlot):
    """Display a model as a histogram."""

    # TODO: should this include
    # histo_prefs = shplot.basicbackend.get_model_histo_defaults()

    def __init__(self):
        super().__init__()
        self.title = 'Model'

    def prepare(self,
                data: Data1DInt,
                model: Model,
                stat: Optional[Stat] = None) -> None:
        """Create the plot.

        The stat parameter is ignored.
        """

        plot = data.to_plot(yfunc=model)
        (_, y, _, _, self.xlabel, self.ylabel) = plot

        # taken from get_x from Data1DInt
        indep = data.get_evaluation_indep(filter=True, model=model)

        self.xlo = indep[0]
        self.xhi = indep[1]
        self.y = y[1]


class SourceHistogramPlot(ModelHistogramPlot):
    """Display a source model as a histogram."""

    def __init__(self):
        super().__init__()
        self.title = 'Source'


class PDFPlot(HistogramPlot):
    """Display the probability density of an array.

    See Also
    --------
    CDFPlot
    """

    _fields: list[str] = ["points"] + HistogramPlot._fields
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    def __init__(self):
        self.points = None
        super().__init__()

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the PDF plot."""
        return backend.as_html_pdf(self)

    def prepare(self, points, bins=12, normed=True, xlabel="x", name="x"):
        """Create the data to plot.

        Parameters
        ----------
        points
            The values to display.
        bins : int, optional
            The number of bins to use.
        normed : bool, optional
            See the description of the density parameter in the
            numpy.histogram routine.
        xlabel : optional, string
            The label for the X axis.
        name : optional, string
            Used to create the plot title.

        See Also
        --------
        plot
        """

        self.points = points
        self.y, xx = np.histogram(points, bins=bins, density=normed)
        self.xlo = xx[:-1]
        self.xhi = xx[1:]
        self.ylabel = "probability density"
        self.xlabel = xlabel
        self.title = f"PDF: {name}"


class CDFPlot(Plot):
    """Display the cumulative distribution of an array.

    The cumulative distribution of the data is drawn along with
    vertical lines to indicate the median, 15.87%, and 84.13%
    percentiles.

    See Also
    --------
    PDFPlot

    Examples
    --------

    Show the cumulative distribution of 1000 randomly-distributed
    points from a Gumbel distribution:

    >>> import numpy as np
    >>> rng = np.random.default_rng()  # requires NumPy 1.17 or later
    >>> pts = rng.gumbel(loc=100, scale=20, size=1000)
    >>> plot = CDFPlot()
    >>> plot.prepare(pts)
    >>> plot.plot()

    """

    _fields: list[str] = ["points", "x", "y", "median!", "lower!",
                          "upper!", "xlabel!", "ylabel!", "title!",
                          "plot_prefs!"]
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    median_defaults = {"linestyle": 'dash', "linecolor": 'orange',
                       "linewidth": 1.5}
    """The options used to draw the median line."""

    lower_defaults = {"linestyle": 'dash', "linecolor": 'blue',
                      "linewidth": 1.5}
    """The options used to draw the 15.87% line."""

    upper_defaults = {"linestyle": 'dash', "linecolor": 'blue',
                      "linewidth": 1.5}
    """The options used to draw the 84.13% line."""

    plot_prefs = basicbackend.get_cdf_plot_defaults()
    """The plot options (the CDF and axes)."""

    def __init__(self):
        self.x = None
        self.y = None
        self.points = None
        self.median = None
        self.lower = None
        self.upper = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        super().__init__()

    def __str__(self) -> str:
        return display_fields(self, self._fields)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the CDF plot."""
        return backend.as_html_cdf(self)

    def prepare(self, points, xlabel="x", name="x"):
        """Create the data to plot.

        Parameters
        ----------
        points
            The values to display (they do not have to be sorted).
        xlabel : optional, string
            The label for the X axis.
        name : optional, string
            Used to create the Y-axis label and plot title.

        See Also
        --------
        plot
        """

        self.points = np.asarray(points)
        self.x = np.sort(points)
        (self.median, self.lower,
         self.upper) = get_error_estimates(self.x, True)
        xsize = len(self.x)
        self.y = (np.arange(xsize) + 1.0) / xsize
        self.xlabel = xlabel
        self.ylabel = f"p(<={xlabel})"
        self.title = f"CDF: {name}"

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        """Plot the data.

        This will plot the data sent to the prepare method.

        Parameters
        ----------
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           plot_prefs dictionary.

        See Also
        --------
        prepare, overplot

        """

        super().plot(self.x, self.y, title=self.title,
                     xlabel=self.xlabel, ylabel=self.ylabel,
                     overplot=overplot, clearwindow=clearwindow,
                     **kwargs)

        # Note: the user arguments are not applied to the vertical lines
        #
        super().vline(self.median, overplot=True, clearwindow=False,
                      **self.median_defaults)
        super().vline(self.lower, overplot=True, clearwindow=False,
                      **self.lower_defaults)
        super().vline(self.upper, overplot=True, clearwindow=False,
                      **self.upper_defaults)


class LRHistogram(HistogramPlot):
    "Derived class for creating 1D likelihood ratio distribution plots"


    _fields: list[str] = ["ratios", "lr!"] + HistogramPlot._fields
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    def __init__(self):
        self.ratios = None
        self.lr = None
        self.ppp = None
        super().__init__()

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the LRHistogram plot."""
        return backend.as_html_lr(self)

    def prepare(self, ratios, bins, niter, lr, ppp):
        "Create the data to plot"

        self.ppp = float(ppp)
        self.lr = float(lr)
        y = np.asarray(ratios)
        self.ratios = y
        self.xlo, self.xhi = dataspace1d(y.min(), y.max(),
                                         numbins=bins + 1)[:2]
        y = histogram1d(y, self.xlo, self.xhi)
        self.y = y / float(niter)
        self.title = "Likelihood Ratio Distribution"
        self.xlabel = "Likelihood Ratio"
        self.ylabel = "Frequency"

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        """Plot the data.

        This will plot the data sent to the prepare method.

        Parameters
        ----------
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           plot_prefs dictionary.

        See Also
        --------
        prepare, overplot

        """

        super().plot(overplot=overplot, clearwindow=clearwindow,
                     **kwargs)

        if self.lr is None:
            return

        # Note: the user arguments are not applied to the vertical line
        #
        if self.lr <= self.xhi.max() and self.lr >= self.xlo.min():
            super().vline(self.lr, linecolor="orange",
                          linestyle="solid", linewidth=1.5,
                          overplot=True, clearwindow=False)


# The handling of the overplot/contour keywords is not ideal. To
# address issue #2128 the UI code ended up accessing internal details
# of this class. Ideally the API would be cleaned up but it's not
# currently clear what the requirements are.
#
class SplitPlot(Plot, Contour):
    """Create multiple plots.

    Attributes
    ----------
    rows : int
       Number of rows of plots. The default is 2.
    cols : int
       Number of columns of plots. The default is 1.
    """

    _fields: list[str] = ["rows!", "cols!", "plot_prefs!"]
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    plot_prefs = basicbackend.get_split_plot_defaults()
    "The preferences for the plot."

    def __init__(self, rows=2, cols=1):
        self.reset(rows, cols)
        super().__init__()

    def __str__(self) -> str:
        return display_fields(self, self._fields)

    # Should this use the preferences to set the rows and cols?
    #
    def reset(self, rows=2, cols=1):
        "Prepare for a new set of plots or contours."
        self.rows = rows
        self.cols = cols
        self._reset_used()
        self._current = (0, 0)
        self._cleared_window = False

    def _clear_window(self):
        if not self._cleared_window:
            backend.clear_window()
            self._cleared_window = True

    def _reset_used(self):
        self._used = np.zeros((self.rows, self.cols), np.bool_)

    def _next_subplot(self):
        row, col = np.where(~self._used)
        if row.size != 0:
            row, col = row[0], col[0]
        else:
            self._reset_used()
            row, col = 0, 0
        return row, col

    def addplot(self, plot, *args, **kwargs):
        "Add the plot to the next space."
        row, col = self._next_subplot()
        self.plot(row, col, plot, *args, **kwargs)

    def addcontour(self, plot, *args, **kwargs):
        "Add the contour plot to the next space."
        row, col = self._next_subplot()
        self.contour(row, col, plot, *args, **kwargs)

    def plot(self, row, col, plot, *args, **kwargs):
        "Add the plot in the given space."
        self._clear_window()
        clearaxes = ((not kwargs.get('overplot', False)) and
                     kwargs.get('clearwindow', True))

        backend.set_subplot(row, col, self.rows, self.cols, clearaxes,
                            **self.plot_prefs)

        kwargs['clearwindow'] = False
        plot.plot(*args, **kwargs)

        self._used[row, col] = True
        self._current = (row, col)

    def contour(self, row, col, plot, *args, **kwargs):
        "Add the contour in the given space."
        self._clear_window()
        clearaxes = ((not kwargs.get('overcontour', False)) and
                     kwargs.get('clearwindow', True))

        backend.set_subplot(row, col, self.rows, self.cols, clearaxes,
                            **self.plot_prefs)

        kwargs['clearwindow'] = False
        plot.contour(*args, **kwargs)

        self._used[row, col] = True
        self._current = (row, col)

    def overlayplot(self, plot, *args, **kwargs):
        "Add the plot to the current space without destroying the contents."
        self.overplot(self._current[0], self._current[1], plot, *args,
                      **kwargs)

    def overlaycontour(self, plot, *args, **kwargs):
        "Add the contour to the current space without destroying the contents."
        self.overcontour(self._current[0], self._current[1], plot, *args,
                         **kwargs)

    # FIXME: work on overplot issue
    def overplot(self, row, col, plot, *args, **kwargs):
        "Add a plot to the given space without destroying the contents."
        kwargs['overplot'] = True
        self.plot(row, col, plot, *args, **kwargs)

    def overcontour(self, row, col, plot, *args, **kwargs):
        "Add a contour plot to the given space without destroying the contents."
        kwargs['overcontour'] = True
        self.contour(row, col, plot, *args, **kwargs)


# TODO: move logic from sherpa.ui.utils.Session._plot_jointplot
#       regarding the X axis here
#
class JointPlot(SplitPlot):
    """Multiple plots that share a common axis

    This supports two plots, where the top plot is twice as
    tall as the bottom plot.
    """

    def plottop(self, plot, *args, overplot=False, clearwindow=True,
                **kwargs):

        create = clearwindow and not overplot
        if create:
            self._clear_window()

        backend.set_jointplot(0, 0, self.rows, self.cols,
                              create=create)

        # The code used to check if the plot was an instance of
        # FitPlot, which has been updated to check for the presence
        # of attributes instead.
        #
        if hasattr(plot, 'xlabel'):
            plot.xlabel = ''
        elif hasattr(plot, 'dataplot') and hasattr(plot, 'modelplot'):
            dplot = plot.dataplot
            mplot = plot.modelplot
            if hasattr(dplot, 'xlabel'):
                dplot.xlabel = ''
            if hasattr(mplot, 'xlabel'):
                mplot.xlabel = ''

        kwargs['clearwindow'] = False
        plot.plot(*args, overplot=overplot, **kwargs)

    def plotbot(self, plot, *args, overplot=False, **kwargs):

        backend.set_jointplot(1, 0, self.rows, self.cols,
                              create=False)

        # FIXME: terrible hack to remove title from bottom
        plot.title = ''
        kwargs['clearwindow'] = False
        plot.plot(*args, overplot=overplot, **kwargs)


class DataPlot(Plot):
    """Create 1D plots of data values.

    Attributes
    ----------
    x : array_like
       The X value for each point (the independent variable).
    y : array_like
       The Y value for each point (the dependent variable).
    xerr : array_like
       The half-width of each X "bin", if set.
    yerr : array_like
       The error on the Y value, if set.
    xlabel, ylabel, title : str
       Plot labels.

    Examples
    --------

    Plot up an example dataset. The default appearance is to draw
    a symbol at each point, but no line connecting the points (the
    actual choice depends on the plot backend):

    >>> from sherpa.data import Data1D
    >>> from sherpa.plot import DataPlot
    >>> data = Data1D('a dataset', [10, 20, 25], [2, -7, 4])
    >>> dplot = DataPlot()
    >>> dplot.prepare(data)
    >>> dplot.plot()

    The plot attributes can be changed to adjust the appearance of
    the plot, and the data re-drawn. The following also shows how
    the plot preferences can be over-ridden to turn off the points
    and draw a dotted line connecting the points (this assumes that
    the Matplotlib backend is in use):

    >>> dplot.xlabel = 'length'
    >>> dplot.ylabel = 'data values'
    >>> dplot.plot(marker=' ', linestyle='dashed')

    """

    _fields: list[str] = ["x", "y", "yerr", "xerr", "xlabel!",
                          "ylabel!", "title!", "plot_prefs!"]
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    plot_prefs = basicbackend.get_data_plot_defaults()
    "The preferences for the plot."

    def __init__(self):
        self.x = None
        self.y = None
        self.yerr = None
        self.xerr = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        super().__init__()

    def __str__(self) -> str:
        return display_fields(self, self._fields)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the data plot."""
        return backend.as_html_data(self)

    def prepare(self, data, stat=None):
        """Create the data to plot

        Parameters
        ----------
        data
            The Sherpa data object to display (it is assumed to be
            one dimensional).
        stat : optional
            The Sherpa statistics object to use to add Y error bars
            if the data has none.

        See Also
        --------
        plot
        """

        (self.x, self.y, self.yerr, self.xerr, self.xlabel,
         self.ylabel) = data.to_plot()

        if stat is not None:
            yerrorbars = self.plot_prefs.get('yerrorbars', True)
            self.yerr = calculate_errors(data, stat, yerrorbars)

        self.title = data.name

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        """Plot the data.

        This will plot the data sent to the prepare method.

        Parameters
        ----------
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           plot_prefs dictionary.

        See Also
        --------
        prepare, overplot

        """
        super().plot(self.x, self.y, yerr=self.yerr, xerr=self.xerr,
                     title=self.title, xlabel=self.xlabel,
                     ylabel=self.ylabel, overplot=overplot,
                     clearwindow=clearwindow, **kwargs)


class TracePlot(DataPlot):

    plot_prefs = basicbackend.get_model_plot_defaults()

    def prepare(self, points, xlabel="x", name="x"):
        """The data to plot.

        Parameters
        ----------
        points
            The trace to plot (one dimensional). The X axis is
            set to consecutive integers, starting at 0.
        xlabel : optional
            This value is ignored.
        name : optional, string
            Used to create the Y-axis label and plot title.

        See Also
        --------
        plot
        """
        self.x = np.arange(len(points), dtype=SherpaFloat)
        self.y = points
        self.xlabel = "iteration"
        self.ylabel = name
        self.title = f"Trace: {name}"


class ScatterPlot(DataPlot):

    plot_prefs = basicbackend.get_scatter_plot_defaults()

    def prepare(self, x, y, xlabel="x", ylabel="y", name="(x,y)"):
        """The data to plot.

        Parameters
        ----------
        x, y
            The points to plot. The number of points in the two
            sequences should be the same.
        xlabel, ylabel, name : optional, str
            The axis labels and name for the plot title.

        See Also
        --------
        plot
        """
        self.x = np.asarray(x, dtype=SherpaFloat)
        self.y = np.asarray(y, dtype=SherpaFloat)
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = f"Scatter: {name}"


class PSFKernelPlot(DataPlot):
    "Derived class for creating 1D PSF kernel data plots"

    def prepare(self, psf, data=None, stat=None):
        """Create the data to plot

        Parameters
        ----------
        psf
            The Sherpa PSF to display
        data : optional
            The Sherpa data object to pass to the PSF.
        stat : optional
            The Sherpa statistics object to use to add Y error bars
            if the data has none.

        See Also
        --------
        plot
        """
        psfdata = psf.get_kernel(data)
        super().prepare(data=psfdata, stat=stat)
        # self.ylabel = 'PSF value'
        # self.xlabel = 'PSF Kernel size'
        self.title = 'PSF Kernel'


class DataContour(Contour):
    """Create contours of 2D data.

    Attributes
    ----------
    x0, x1 : array_like
       The coordinates of each point (the independent variables), as
       one-dimensional arrays.
    y : array_like
       The Y value for each point (the dependent variable), as a
       one-dimensional array.
    levels : array_like or `None`
       The values at which to draw contours.
    xlabel, ylabel, title : str
       Plot labels.

    """

    _fields: list[str] = ["x0", "x1", "y", "xlabel!", "ylabel!",
                          "title!", "levels!", "contour_prefs!"]
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    contour_prefs = basicbackend.get_data_contour_defaults()
    "The preferences for the plot."

    def __init__(self):
        self.x0 = None
        self.x1 = None
        self.y = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        self.levels = None
        super().__init__()

    def __str__(self) -> str:
        return display_fields(self, self._fields)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the contour plot."""
        return backend.as_html_datacontour(self)

    def prepare(self, data, stat=None):
        (self.x0, self.x1, self.y, self.xlabel,
         self.ylabel) = data.to_contour()
        self.title = data.name

    def contour(self, overcontour=False, clearwindow=True, **kwargs):
        super().contour(self.x0, self.x1, self.y, levels=self.levels,
                        title=self.title, xlabel=self.xlabel,
                        ylabel=self.ylabel, overcontour=overcontour,
                        clearwindow=clearwindow, **kwargs)


class PSFKernelContour(DataContour):
    "Derived class for creating 2D PSF Kernel contours"

    def prepare(self, psf, data=None, stat=None):
        psfdata = psf.get_kernel(data)
        super().prepare(data=psfdata)
        # self.xlabel = 'PSF Kernel size x0'
        # self.ylabel = 'PSF Kernel size x1'
        self.title = 'PSF Kernel'


class ModelPlot(Plot):
    """Create 1D plots of model values.

    Attributes
    ----------
    x : array_like
       The X value for each point (the independent variable).
    y : array_like
       The Y value for each point (the dependent variable).
    xerr : array_like
       The half-width of each X "bin", if set.
    yerr : array_like
       The error on the Y value, if set.
    xlabel, ylabel, title : str
       Plot labels.

    Examples
    --------

    Plot up an example dataset. The default appearance is to draw
    a line between each point:

    >>> from sherpa.data import Data1D
    >>> from sherpa.models.basic import StepLo1D
    >>> from sherpa.plot import ModelPlot
    >>> data = Data1D('a dataset', [10, 20, 25], [2, -7, 4])
    >>> model = StepLo1D()
    >>> model.xcut = 19
    >>> mplot = ModelPlot()
    >>> mplot.prepare(data, model)
    >>> mplot.plot()

    The plot attributes can be changed to adjust the appearance of
    the plot, and the data re-drawn. The following also shows how
    the plot preferences can be over-ridden to draw the model as
    squares with no line connecting them (this assumes that
    the Matplotlib backend is in use):

    >>> mplot.xlabel = 'length'
    >>> mplot.ylabel = 'model values'
    >>> mplot.plot(marker='s', linestyle='none')

    """

    _fields: list[str] = ["x", "y", "yerr", "xerr", "xlabel!",
                          "ylabel!", "title!", "plot_prefs!"]
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    plot_prefs = basicbackend.get_model_plot_defaults()
    "The preferences for the plot."

    def __init__(self):
        self.x = None
        self.y = None
        self.yerr = None
        self.xerr = None
        self.xlabel = None
        self.ylabel = None
        self.title = 'Model'
        super().__init__()

    def __str__(self) -> str:
        return display_fields(self, self._fields)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the model plot."""
        return backend.as_html_model(self)

    def prepare(self,
                data: Data1D,
                model: Model,
                stat: Optional[Stat] = None) -> None:
        """Create the data to plot

        Parameters
        ----------
        data
            The Sherpa data object to display (it is assumed to be
            one dimensional). This defines the grid over which the
            model is displayed.
        model
            The Sherpa model expression to evaluate and display.
        stat : optional
            This parameter is unused.

        See Also
        --------
        plot
        """
        (self.x, self.y, self.yerr, self.xerr,
         self.xlabel, self.ylabel) = data.to_plot(yfunc=model)
        self.y = self.y[1]

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        """Plot the data.

        This will plot the data sent to the prepare method.

        Parameters
        ----------
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           plot_prefs dictionary.

        See Also
        --------
        prepare, overplot

        """
        # This does not display yerr or xerr.
        super().plot(self.x, self.y, title=self.title,
                     xlabel=self.xlabel, ylabel=self.ylabel,
                     overplot=overplot, clearwindow=clearwindow,
                     **kwargs)


class ComponentModelPlot(ModelPlot):
    """Component model plots."""

    plot_prefs = basicbackend.get_component_plot_defaults()
    "The preferences for the plot."

    # TODO: ComponentSourcePlot calls data.to_component_plot but this
    # just uses data.to_plot. Does it matter?
    #
    def prepare(self, data, model, stat=None):
        super().prepare(data=data, model=model, stat=stat)
        self.title = f'Model component: {model.name}'


class ComponentModelHistogramPlot(ModelHistogramPlot):
    """Component model plots for histogram data.

    ..versionchanged:: 4.16.1
      The `histo_prefs` attribute is now properly set (plot_prefs is
      also no-longer set).

    """

    histo_prefs = basicbackend.get_component_histo_defaults()
    "The preferences for the plot."

    # TODO: ComponentSourceHistogramPlot calls data.to_component_plot
    # but this just uses data.to_plot. Does it matter?
    #
    def prepare(self, data, model, stat=None):
        super().prepare(data=data, model=model, stat=stat)
        self.title = f'Model component: {model.name}'


class ComponentTemplateModelPlot(ComponentModelPlot):

    def prepare(self, data, model, stat=None):
        self.x = model.get_x()
        self.y = model.get_y()
        self.xlabel = data.get_xlabel()
        self.ylabel = data.get_ylabel()
        self.title = f'Model component: {model.name}'


class SourcePlot(ModelPlot):
    """Create 1D plots of unconvolved model values.

    Attributes
    ----------
    x : array_like
       The X value for each point (the independent variable).
    y : array_like
       The Y value for each point (the model value).
    xerr : array_like
       The half-width of each X "bin", if set.
    yerr : array_like
       The error on the Y value, if set.
    xlabel, ylabel, title : str
       Plot labels.

    """

    def __init__(self):
        super().__init__()
        self.title = 'Source'


class ComponentSourcePlot(SourcePlot):
    """Component source plots."""

    plot_prefs = basicbackend.get_component_plot_defaults()

    def prepare(self, data, model, stat=None):
        (self.x, self.y, self.yerr, self.xerr,
         self.xlabel, self.ylabel) = data.to_component_plot(yfunc=model)
        self.y = self.y[1]
        self.title = f'Source model component: {model.name}'


class ComponentSourceHistogramPlot(SourceHistogramPlot):
    """Component source plots for histogram data.

    ..versionchanged:: 4.16.1
      The `histo_prefs` attribute is now properly set (plot_prefs is
      also no-longer set).

    """

    histo_prefs = basicbackend.get_component_histo_defaults()

    def prepare(self, data, model, stat=None):

        plot = data.to_component_plot(yfunc=model)
        (_, y, _, _, self.xlabel, self.ylabel) = plot

        # taken from get_x from Data1DInt
        indep = data.get_evaluation_indep(filter=True, model=model)

        self.xlo = indep[0]
        self.xhi = indep[1]
        self.y = y[1]

        self.title = f'Source model component: {model.name}'


class ComponentTemplateSourcePlot(ComponentSourcePlot):

    def prepare(self, data, model, stat=None):
        self.x = model.get_x()
        self.y = model.get_y()

        if np.iterable(data.mask):
            x = data.to_plot()[0]
            mask = (self.x > x.min()) & (self.x <= x.max())
            self.x = self.x[mask]
            self.y = self.y[mask]

        self.xlabel = data.get_xlabel()
        self.ylabel = data.get_ylabel()
        self.title = f'Source model component: {model.name}'


class PSFPlot(DataPlot):
    "Derived class for creating 1D PSF kernel data plots"

    def prepare(self, psf, data=None, stat=None):
        """Create the data to plot

        Parameters
        ----------
        psf
            The Sherpa PSF to display
        data : optional
            The Sherpa data object to pass to the PSF.
        stat : optional
            The Sherpa statistics object to use to add Y error bars
            if the data has none.

        See Also
        --------
        plot
        """
        psfdata = psf.get_kernel(data, False)
        super().prepare(data=psfdata, stat=stat)
        self.title = psf.kernel.name


class ModelContour(Contour):
    "Derived class for creating 2D model contours"

    _fields: list[str] = ["x0", "x1", "y", "xlabel!", "ylabel!",
                          "title!", "levels!", "contour_prefs!"]
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    contour_prefs = basicbackend.get_model_contour_defaults()
    "The preferences for the plot."

    def __init__(self):
        self.x0 = None
        self.x1 = None
        self.y = None
        self.xlabel = None
        self.ylabel = None
        self.title = 'Model'
        self.levels = None
        super().__init__()

    def __str__(self) -> str:
        return display_fields(self, self._fields)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the model contour plot."""
        return backend.as_html_modelcontour(self)

    def prepare(self, data, model, stat=None):
        """Prepare the data.

        .. versionchanged:: 4.16.1
           The stat argument is now unused and will likely be removed
           in a future release.

        """

        (self.x0, self.x1, self.y, self.xlabel,
         self.ylabel) = data.to_contour(yfunc=model)
        self.y = self.y[1]

    def contour(self, overcontour=False, clearwindow=True, **kwargs):
        super().contour(self.x0, self.x1, self.y, levels=self.levels,
                        title=self.title, xlabel=self.xlabel,
                        ylabel=self.ylabel, overcontour=overcontour,
                        clearwindow=clearwindow, **kwargs)


class PSFContour(DataContour):
    "Derived class for creating 2D PSF contours"

    def prepare(self, psf, data=None, stat=None):
        psfdata = psf.get_kernel(data, False)
        super().prepare(data=psfdata)
        self.title = psf.kernel.name


class SourceContour(ModelContour):
    "Derived class for creating 2D model contours"

    def __init__(self):
        super().__init__()
        self.title = 'Source'


class FitPlot(Plot):
    """Combine data and model plots for 1D data.

    Attributes
    ----------
    dataplot
       The Sherpa plot object used to display the data.
    modelplot
       The Sherpa plot object used to display the model.

    Examples
    --------

    If dplot and mplot are data and model plots respectively,
    such as those created in the examples for `DataPlot` and `ModelPlot`,
    then a combined plot can be created with:

    >>> from sherpa.data import Data1D
    >>> from sherpa.models.basic import StepLo1D
    >>> from sherpa.plot import DataPlot, ModelPlot, FitPlot
    >>> data = Data1D('a dataset', [10, 20, 25], [2, -7, 4])
    >>> dplot = DataPlot()
    >>> fplot = FitPlot()
    >>> model = StepLo1D()
    >>> model.xcut = 19
    >>> mplot = ModelPlot()
    >>> dplot.prepare(data)
    >>> mplot.prepare(data, model)
    >>> fplot.prepare(dplot, mplot)
    >>> fplot.plot()


    Keyword arguments can be given in the `plot` call, and these
    are passed through to both the data and model plots (in the
    following example the Matplotlib backend is assumed to be in use):

    >>> fplot.plot(color='k')

    """

    plot_prefs = basicbackend.get_fit_plot_defaults()
    """The preferences for the plot. Note that the display for the
    data and model plots are controlled by the preferences for
    the dataplot and modelplot objects, so this is currently
    unused.
    """

    def __init__(self):
        self.dataplot = None
        self.modelplot = None
        super().__init__()

    def __str__(self) -> str:
        data_title = None
        if self.dataplot is not None:
            data_title = self.dataplot.title

        model_title = None
        if self.modelplot is not None:
            model_title = self.modelplot.title

        return f"""dataplot   = {data_title}
{self.dataplot}

modelplot  = {model_title}
{self.modelplot}"""

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the fit plot."""
        return backend.as_html_fit(self)

    def prepare(self, dataplot, modelplot):
        """Create the data to plot

        Parameters
        ----------
        dataplot
            The Sherpa plot object representing the data (normally a
            DataPlot instance).
        modelplot
            The Sherpa plot object representing the model (normally a
            ModelPlot instance).

        See Also
        --------
        plot
        """
        self.dataplot = dataplot
        self.modelplot = modelplot

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        """Plot the data.

        This will plot the data sent to the prepare method.

        Parameters
        ----------
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           plot_prefs dictionary. They are sent to both the dataplot
           and modelplot instances.

        See Also
        --------
        prepare, overplot

        """
        if self.dataplot is None or self.modelplot is None:
            raise PlotErr("nodataormodel")

        # Note: the user preferences are sent to *both* the data and
        #       model plot. Is this a good idea?
        #
        self.dataplot.plot(overplot=overplot, clearwindow=clearwindow,
                           **kwargs)
        self.modelplot.overplot(**kwargs)


class FitContour(Contour):
    "Derived class for creating 2D combination data and model contours"

    contour_prefs = basicbackend.get_fit_contour_defaults()
    "The preferences for the plot."

    def __init__(self):
        self.datacontour = None
        self.modelcontour = None
        super().__init__()

    def __str__(self) -> str:
        data_title = None
        if self.datacontour is not None:
            data_title = self.datacontour.title

        model_title = None
        if self.modelcontour is not None:
            model_title = self.modelcontour.title

        return f"""datacontour = {data_title}
{self.datacontour}

modelcontour = {model_title}
{self.modelcontour}"""

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the fit contour plot."""
        return backend.as_html_fitcontour(self)

    def prepare(self, datacontour, modelcontour):
        self.datacontour = datacontour
        self.modelcontour = modelcontour

    def contour(self, overcontour=False, clearwindow=True, **kwargs):
        # Note: the user arguments are applied to both plots
        self.datacontour.contour(overcontour=overcontour,
                                 clearwindow=clearwindow, **kwargs)
        self.modelcontour.overcontour(**kwargs)



class BaseResidualPlot(ModelPlot):
    """Residuals of model + data.

    Subclasses need to implement `_calc_y` and `_title`.

    .. versionadded:: 4.16.1

    Notes
    -----
    The ylog setting is ignored, whether given as a preference or
    a keyword argument, so the Y axis is always drawn with a
    linear scale.

    Errors are supported on the Y axis.
    """

    residual_axis: float = 0
    """At what point on the y axis is the residual line drawn."""

    residual_color: str = 'k'
    """The color of the residual line."""

    residual_lw: float = 0.8
    """The line width of the residual line."""

    def _calc_y(self,
                data: Data1D,
                stat: Stat,
                ylist: tuple[np.ndarray, np.ndarray],
                staterr: Optional[np.ndarray]) -> None:
        """Define the self.y and self.yerr fields"""
        raise NotImplementedError()

    def _change_ylabel(self) -> None:
        """Change the self.ylabel field, if needed."""
        pass

    def _title(self, data: Data1D) -> None:
        """Set the self.title field"""
        raise NotImplementedError()

    def prepare(self,  # type: ignore[override]
                data: Data1D,
                model: Model,
                # Note that the superclass has stat as optional
                stat: Stat) -> None:

        plot = data.to_plot(yfunc=model)
        (self.x, y, staterr, self.xerr, self.xlabel, self.ylabel) = plot

        self._calc_y(data, stat, y, staterr)
        self._change_ylabel()
        self._title(data)

    def plot(self,  # type: ignore[override]
             overplot: bool = False,
             clearwindow: bool = True,
             **kwargs) -> None:

        x = self.x
        y = self.y

        # The y axis is only drawn with a linear scale.
        #
        self.plot_prefs['ylog'] = False
        kwargs.pop('ylog', True)

        # The superclass is not called since the y errors are
        # potentially included in the plot.
        #
        Plot.plot(self, x, y, yerr=self.yerr, xerr=self.xerr,
                  title=self.title, xlabel=self.xlabel,
                  ylabel=self.ylabel, overplot=overplot,
                  clearwindow=clearwindow, **kwargs)
        super().hline(y=self.residual_axis, xmin=0, xmax=1,
                      linecolor=self.residual_color,
                      linewidth=self.residual_lw,
                      overplot=True)


class BaseResidualHistogramPlot(ModelHistogramPlot):
    """Residuals of model + data for histogram data.

    Subclasses need to implement `_calc_y` and `_title`.

    .. versionadded:: 4.16.1

    Notes
    -----
    The ylog setting is ignored, whether given as a preference or
    a keyword argument, so the Y axis is always drawn with a
    linear scale.

    Errors are supported on the Y axis.
    """

    residual_axis: float = 0
    """At what point on the y axis is the residual line drawn."""

    residual_color: str = 'k'
    """The color of the residual line."""

    residual_lw: float = 0.8
    """The line width of the residual line."""

    def __init__(self) -> None:
        self.yerr: Optional[np.ndarray]
        self.yerr = None
        super().__init__()

    def _calc_x(self, data: Data1DInt, model: Model) -> None:
        """Define the xlo and xhi fields"""

        # taken from get_x from Data1DInt
        indep = data.get_evaluation_indep(filter=True, model=model)
        self.xlo = indep[0]
        self.xhi = indep[1]

    def _calc_y(self,
                data: Data1DInt,
                stat: Stat,
                ylist: tuple[np.ndarray, np.ndarray],
                staterr: Optional[np.ndarray]) -> None:
        """Define the self.y and self.yerr fields"""
        raise NotImplementedError()

    def _change_ylabel(self) -> None:
        """Change the self.ylabel field, if needed."""
        pass

    def _title(self, data: Data1DInt) -> None:
        """Set the self.title field"""
        raise NotImplementedError()

    def prepare(self,  # type: ignore[override]
                data: Data1DInt,
                model: Model,
                # Note that the superclass has stat as optional
                stat: Stat) -> None:

        plot = data.to_plot(model)
        (_, y, staterr, _, self.xlabel, self.ylabel) = plot

        self._calc_x(data, model)
        self._calc_y(data, stat, y, staterr)
        self._change_ylabel()
        self._title(data)

    def plot(self,  # type: ignore[override]
             overplot: bool = False,
             clearwindow: bool = True,
             **kwargs) -> None:

        xlo = self.xlo
        xhi = self.xhi
        y = self.y

        # The y axis is only drawn with a linear scale.
        #
        self.histo_prefs['ylog'] = False
        kwargs.pop('ylog', True)

        # The superclass is not called since the y errors are
        # potentially included in the plot.
        #
        Histogram.plot(self, xlo, xhi, y, yerr=self.yerr,
                       title=self.title, xlabel=self.xlabel,
                       ylabel=self.ylabel, overplot=overplot,
                       clearwindow=clearwindow, **kwargs)
        super().hline(y=self.residual_axis, xmin=0, xmax=1,
                      linecolor=self.residual_color,
                      linewidth=self.residual_lw,
                      overplot=True)


class BaseResidualContour(ModelContour):
    """Residuals of model + data for contour data.

    Subclasses need to implement `_calc_y` and `_title`.

    .. versionadded:: 4.16.1

    """

    def _calc_y(self,
                ylist: tuple[np.ndarray, np.ndarray]) -> None:
        """Define the self.y field"""
        raise NotImplementedError()

    def _title(self, data: Data2D) -> None:
        """Set the self.title field"""
        raise NotImplementedError()

    def prepare(self,
                data: Data2D,
                model: Model,
                stat = None) -> None:

        (self.x0, self.x1, ys, self.xlabel,
         self.ylabel) = data.to_contour(yfunc=model)

        self._calc_y(ys)
        self._title(data)


class DelchiPlot(BaseResidualPlot):
    """Create plots of the delta-chi value per point.

    The value of (data-model)/error is plotted for each point.

    Attributes
    ----------
    x : array_like
       The X value for each point.
    y : array_like
       The Y value for each point: (data-model)/error
    xerr : array_like
       The half-width of each X "bin", if set.
    yerr : array_like
       The error on the Y value (each element is `1`).
    xlabel, ylabel, title : str
       Plot labels.

    Notes
    -----
    The ylog setting is ignored, whether given as a preference or
    a keyword argument, so the Y axis is always drawn with a
    linear scale.
    """

    plot_prefs = basicbackend.get_resid_plot_defaults()
    "The preferences for the plot."

    def _calc_y(self,
                data: Data1D,
                stat: Stat,
                ylist: tuple[np.ndarray, np.ndarray],
                staterr: Optional[np.ndarray]) -> None:
        """Define the self.y and self.yerr fields"""

        if staterr is None:
            if stat.name in _stats_noerr:
                raise StatErr('badstat', "DelchiPlot", stat.name)
            staterr = data.get_yerr(True, stat.calc_staterror)

        self.y = (ylist[0] - ylist[1]) / staterr
        self.yerr = np.ones_like(self.y)

    def _change_ylabel(self) -> None:
        self.ylabel = 'Sigma'

    def _title(self, data: Data1D) -> None:
        self.title = _make_title('Sigma Residuals', data.name)


class DelchiHistogramPlot(BaseResidualHistogramPlot):
    """Create plots of the delta-chi value per point."""

    histo_prefs = basicbackend.get_resid_histo_defaults()
    "The preferences for the delchi plot."

    def _calc_y(self,
                data: Data1DInt,
                stat: Stat,
                ylist: tuple[np.ndarray, np.ndarray],
                staterr: Optional[np.ndarray]) -> None:

        if staterr is None:
            if stat.name in _stats_noerr:
                raise StatErr('badstat', "DelchiHistogramPlot", stat.name)

            staterr = data.get_yerr(True, stat.calc_staterror)

        self.y = (ylist[0] - ylist[1]) / staterr
        self.yerr = staterr / staterr

    def _change_ylabel(self) -> None:
        self.ylabel = 'Sigma'

    def _title(self, data: Data1DInt) -> None:
        self.title = _make_title('Sigma Residuals', data.name)


class ChisqrPlot(ModelPlot):
    """Create plots of the chi-square value per point.

    The value of ((data-model)/error)^2 is plotted for each point.

    Attributes
    ----------
    x : array_like
       The X value for each point.
    y : array_like
       The Y value for each point: ((data-model)/error)^2
    xerr : array_like
       The half-width of each X "bin", if set.
    yerr : array_like
       The error on the Y value. Will be `None` here.
    xlabel, ylabel, title : str
       Plot labels.
    """

    plot_prefs = basicbackend.get_model_plot_defaults()
    "The preferences for the plot."

    def _calc_chisqr(self, ylist, staterr):
        dy = ylist[0] - ylist[1]
        return dy * dy / (staterr * staterr)

    def prepare(self, data, model, stat):
        (self.x, y, _, self.xerr,
         self.xlabel, self.ylabel) = data.to_plot(model)

        if stat.name in _stats_noerr:
            raise StatErr('badstat', "ChisqrPlot", stat.name)
        staterr = data.get_yerr(True, stat.calc_staterror)

        self.y = self._calc_chisqr(y, staterr)
        # TODO: should this
        #   self.yerr = staterr * staterr
        # but we currently do not display errors
        self.ylabel = backend.get_latex_for_string(r'\chi^2')
        self.title = _make_title(
            backend.get_latex_for_string(r'\chi^2'), data.name)


class ChisqrHistogramPlot(ModelHistogramPlot):
    """Create plot of the ChiSq residuals per point."""

    histo_prefs = basicbackend.get_resid_histo_defaults()
    "The preferences for the residual plot."

    def _calc_x(self, data: Data1DInt, model: Model) -> None:
        """Define the xlo and xhi fields"""

        # taken from get_x from Data1DInt
        indep = data.get_evaluation_indep(filter=True, model=model)
        self.xlo = indep[0]
        self.xhi = indep[1]

    def _calc_chisqr(self,
                     ylist: tuple[np.ndarray, np.ndarray],
                     staterr: np.ndarray) -> np.ndarray:
        dy = ylist[0] - ylist[1]
        return dy * dy / (staterr * staterr)

    def prepare(self,  # type: ignore[override]
                data: Data1DInt,
                model: Model,
                stat: Stat) -> None:
        plot = data.to_plot(model)
        (_, y, _, _, self.xlabel, self.ylabel) = plot

        self._calc_x(data, model)

        # if staterr is None:
        if stat.name in _stats_noerr:
            raise StatErr('badstat', "ChisqrHistogramPlot", stat.name)

        staterr = data.get_yerr(True, stat.calc_staterror)
        self.y = self._calc_chisqr(y, staterr)

        self.ylabel = backend.get_latex_for_string(r'\chi^2')
        self.title = _make_title(
            backend.get_latex_for_string(r'\chi^2'), data.name)


class ResidPlot(BaseResidualPlot):
    """Create plots of the residuals (data - model) per point.

    The value of (data-model) is plotted for each point.

    Attributes
    ----------
    x : array_like
       The X value for each point.
    y : array_like
       The Y value for each point: data - model.
    xerr : array_like
       The half-width of each X "bin", if set.
    yerr : array_like
       The error on the Y value, if set.
    xlabel, ylabel, title : str
       Plot labels.

    Notes
    -----
    The ylog setting is ignored, whether given as a preference or
    a keyword argument, so the Y axis is always drawn with a
    linear scale.
    """

    plot_prefs = basicbackend.get_resid_plot_defaults()
    "The preferences for the plot."

    def _calc_y(self,
                data: Data1D,
                stat: Stat,
                ylist: tuple[np.ndarray, np.ndarray],
                staterr: Optional[np.ndarray]) -> None:
        """Define the self.y and self.yerr fields"""

        self.y = ylist[0] - ylist[1]

        if stat.name in _stats_noerr:
            calcfn = Chi2XspecVar.calc_staterror
            if self.plot_prefs.get('yerrorbars', True):
                warning(_errorbar_warning(stat))
        else:
            calcfn = stat.calc_staterror

        self.yerr = data.get_yerr(True, calcfn)

    def _change_ylabel(self) -> None:
        self.ylabel = 'Data - Model'

    def _title(self, data: Data1D) -> None:
        self.title = _make_title('Residuals', data.name)


class ResidHistogramPlot(BaseResidualHistogramPlot):
    """Create plots of the residual value per point."""

    histo_prefs = basicbackend.get_resid_histo_defaults()
    "The preferences for the residual plot."

    def _calc_y(self,
                data: Data1DInt,
                stat: Stat,
                ylist: tuple[np.ndarray, np.ndarray],
                staterr: Optional[np.ndarray]) -> None:

        self.y = ylist[0] - ylist[1]

        # See the discussion in DataPlot.prepare
        try:
            yerrorbars = self.histo_prefs['yerrorbars']
        except KeyError:
            yerrorbars = True

        if stat.name in _stats_noerr:
            calcyerr = Chi2XspecVar.calc_staterror
            if yerrorbars:
                warning(_errorbar_warning(stat))
        else:
            calcyerr = stat.calc_staterror

        self.yerr = data.get_yerr(True, calcyerr)

    def _change_ylabel(self) -> None:
        self.ylabel = 'Data - Model'

    def _title(self, data: Data1DInt) -> None:
        self.title = _make_title('Residuals', data.name)


class ResidContour(BaseResidualContour):
    "Derived class for creating 2D residual contours (data-model)"

    contour_prefs = basicbackend.get_resid_contour_defaults()
    "The preferences for the plot."

    def _calc_y(self,
                ylist: tuple[np.ndarray, np.ndarray]) -> None:
        self.y = ylist[0] - ylist[1]

    def _title(self, data):
        self.title = _make_title('Residuals', data.name)


class RatioPlot(BaseResidualPlot):
    """Create plots of the ratio of data to model per point.

    The value of data / model is plotted for each point.

    Attributes
    ----------
    x : array_like
       The X value for each point.
    y : array_like
       The Y value for each point: data / model.
    xerr : array_like
       The half-width of each X "bin", if set.
    yerr : array_like
       The error on the Y value, if set.
    xlabel, ylabel, title : str
       Plot labels.

    Notes
    -----
    The ylog setting is ignored, whether given as a preference or
    a keyword argument, so the Y axis is always drawn with a
    linear scale.
    """

    plot_prefs = basicbackend.get_ratio_plot_defaults()
    "The preferences for the plot."

    residual_axis = 1
    """At what point on the y axis is the residual line drawn."""

    def _calc_y(self,
                data: Data1D,
                stat: Stat,
                ylist: tuple[np.ndarray, np.ndarray],
                staterr: Optional[np.ndarray]) -> None:
        """Define the self.y and self.yerr fields"""

        dvals = np.array(ylist[0])
        mvals = np.asarray(ylist[1])
        bad = np.where(mvals == 0.0)
        dvals[bad] = 0.0
        mvals[bad] = 1.0
        self.y = dvals / mvals

        if stat.name in _stats_noerr:
            calcfn = Chi2XspecVar.calc_staterror
            if self.plot_prefs.get('yerrorbars', True):
                warning(_errorbar_warning(stat))
        else:
            calcfn = stat.calc_staterror

        self.yerr = data.get_yerr(True, calcfn)
        self.yerr = self.yerr / ylist[1]

    def _change_ylabel(self) -> None:
        self.ylabel = 'Data / Model'

    def _title(self, data: Data1D) -> None:
        self.title = _make_title('Ratio of Data to Model', data.name)


class RatioHistogramPlot(BaseResidualHistogramPlot):
    """Create plots of the ratio value per point."""

    # Turn on yerrorbars by default. What is the best way to do this?
    #
    histo_prefs = basicbackend.get_resid_histo_defaults() | {"yerrorbars": True}
    "The preferences for the ratio plot."

    residual_axis = 1

    def _calc_y(self,
                data: Data1DInt,
                stat: Stat,
                ylist: tuple[np.ndarray, np.ndarray],
                staterr: Optional[np.ndarray]) -> None:

        # should not need np.asarray here but leave in for safety
        dvals = np.asarray(ylist[0])
        mvals = np.asarray(ylist[1])
        bad = np.where(mvals == 0.0)
        dvals[bad] = 0.0
        mvals[bad] = 1.0
        self.y = dvals / mvals

        # See the discussion in DataPlot.prepare
        try:
            yerrorbars = self.histo_prefs['yerrorbars']
        except KeyError:
            yerrorbars = True

        if stat.name in _stats_noerr:
            calcyerr = Chi2XspecVar.calc_staterror
            if yerrorbars:
                warning(_errorbar_warning(stat))
        else:
            calcyerr = stat.calc_staterror

        self.yerr = data.get_yerr(True, calcyerr) / ylist[1]

    def _change_ylabel(self) -> None:
        self.ylabel = 'Data / Model'

    def _title(self, data: Data1DInt) -> None:
        self.title = _make_title('Ratio of Data to Model', data.name)


class RatioContour(BaseResidualContour):
    "Derived class for creating 2D ratio contours (data divided by model)"

    contour_prefs = basicbackend.get_ratio_contour_defaults()
    "The preferences for the plot."

    def _calc_y(self,
                ylist: tuple[np.ndarray, np.ndarray]) -> None:

        data = np.array(ylist[0])
        model = np.asarray(ylist[1])
        bad = np.where(model == 0.0)
        data[bad] = 0.0
        model[bad] = 1.0
        self.y = data / model

    def _title(self, data: Data2D) -> None:
        self.title = _make_title('Ratio of Data to Model', data.name)


def calc_par_range(minval: float,
                   maxval: float,
                   nloop: int,
                   delv: Optional[float] = None,
                   log: Optional[bool] = False) -> np.ndarray:
    """Calculate the parameter range to use.

    This assumes that the arguments have already been checked for
    validity.

    Parameters
    ----------
    minval, maxval : number
       The requested parameter range, with maxval > minval.
    nloop : int
       The number of values to create. This is not used if delv is
       set.
    delv : float or None, optional
       If set this is the spacing to use (and over-rides the nloop
       parameter). It must be positive and, when `log` is set, the
       spacing is 10**delv.
    log : bool, optional
       Should the spacing be linear or logarithmic?

    Returns
    -------
    pars : ndarray
       The parameter values.

    Raises
    ------
    ConfidenceErr
       If either minval or maxval are not positive definite when the
       `log` option is set.

    """

    if delv is not None:
        # This assumes that delv > eps, but this should be safe.
        # The float(eps) term is for mypy.
        #
        eps = np.finfo(np.float32).eps
        maxval = maxval + float(eps)

    if log:
        if minval <= 0.0 or maxval <= 0.0:
            raise ConfidenceErr('badarg', 'Log scale',
                                'on positive boundaries')

        minval = np.log10(minval)
        maxval = np.log10(maxval)

    if delv is None:
        x = np.linspace(minval, maxval, nloop)
    else:
        x = np.arange(minval, maxval, delv)

    if log:
        x = 10**x

    return x


class Confidence1D(DataPlot):
    """The base class for 1D confidence plots.

    .. versionchanged:: 4.16.1
       Handling of log-scaled axes and use of the delv argument has
       been improved, and the string output now includes the parameter
       value (if available).

    """



    _fields: list[str] = ["x", "y", "min!", "max!", "nloop!", "delv!",
                          "fac!", "log!", "parval!"]
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    plot_prefs = basicbackend.get_confid_plot_defaults()
    "The preferences for the plot."

    # This is only used for an error message
    conf_type = "unknown"
    "The type of confidence analysis."

    def __init__(self):
        self.min = None
        self.max = None
        self.nloop = 20
        self.delv = None
        self.fac = 1
        self.log = False
        self.parval = None
        self.stat = None
        self.numcores = None
        super().__init__()

    def __setstate__(self, state):
        self.__dict__.update(state)

        if 'stat' not in state:
            self.__dict__['stat'] = None

        if 'parval' not in state:
            self.__dict__['parval'] = None

        if 'numcores' not in state:
            self.__dict__['numcores'] = None

    def __str__(self) -> str:
        return display_fields(self, self._fields)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the confidence 1D plot."""
        return backend.as_html_contour1d(self)

    def prepare(self, min=None, max=None, nloop=20,
                delv=None, fac=1, log=False, numcores=None):
        """Set the data to plot.

        This defines the range over which the statistic will be
        calculated, but does not perform the evaluation.

        Parameters
        ----------
        min, max : number or None, optional
            The minimum and maximum parameter value to used. If either
            is not set then the range is calculated using the fac
            parameter.
        nloop : int, optional
            The number of points at which to evaluate the
            statistic. It must be greater than 1. This is used when
            delv is set to None.
        delv : number or None, optional
            The spacing of the parameter grid. This takes precedence
            over nloop.
        fac : number, optional
            Used when either min or max are not set. The parameter
            range in this case is taken to be fac times the separation
            of the covariance limits for the parameter (unless
            explicitly given).
        log : bool, optional
            Should the parameter be evaluated on a
            logarithmically-spaced grid rather than a linearly-spaced
            one?
        numcores : int or None, optional
            Should the parameter evaluation use multiple CPU cores if
            available?

        See Also
        --------
        calc

        """

        self.min = min
        self.max = max
        self.nloop = nloop
        self.delv = delv
        self.fac = fac
        self.log = log
        self.numcores = numcores

    def _interval_init(self, fit, par):
        """Calculate the grid to use for the parameter.

        Parameters
        ----------
        fit : sherpa.fit.Fit instance
            The current fit.
        par : sherpa.models.parameter.Parameter instance
            The parameter to analyze.

        Returns
        -------
        x : ndarray
            The parameter values to use (also set to self.x)

        """

        # Validate the values first (as much as we can).
        #
        if self.min is not None and not np.isscalar(self.min):
            raise ConfidenceErr('badarg', 'Parameter limits', 'scalars')

        if self.max is not None and not np.isscalar(self.max):
            raise ConfidenceErr('badarg', 'Parameter limits', 'scalars')

        if self.nloop <= 1:
            raise ConfidenceErr('badarg', 'Nloop parameter', '> 1')

        if self.delv is not None and self.delv <= 0:
            raise ConfidenceErr('badarg', 'delv parameter', '> 0')

        self.stat = fit.calc_stat()
        self.parval = par.val
        self.xlabel = par.fullname
        self.ylabel = 'Statistic Value'

        if self.min is None or self.max is None:
            oldestmethod = fit.estmethod
            fit.estmethod = Covariance()
            r = fit.est_errors()
            index = list(r.parnames).index(par.fullname)
            fit.estmethod = oldestmethod

            if self.min is None:
                self.min = par.min
                minval = r.parmins[index]
                if minval is not None and not np.isnan(minval):
                    self.min = par.val + minval

            if self.max is None:
                self.max = par.max
                maxval = r.parmaxes[index]
                if maxval is not None and not np.isnan(maxval):
                    self.max = par.val + maxval

            v = (self.max + self.min) / 2.
            dv = np.fabs(v - self.min)
            self.min = v - self.fac * dv
            self.max = v + self.fac * dv

        if self.min < par.min:
            self.min = par.min
        if self.max > par.max:
            self.max = par.max

        # check user limits for errors
        if self.min >= self.max:
            raise ConfidenceErr('badlimits')

        self.x = calc_par_range(self.min, self.max, self.nloop,
                                delv=self.delv, log=self.log)
        return self.x

    def calc(self, fit, par):
        """Evaluate the statistic for the parameter range.

        This requires prepare to have been called, and must be
        called before plot is called.

        Parameters
        ----------
        fit
            The Sherpa fit instance to use (defines the statistic
            and optimiser to use).
        par
            The parameter to iterate over.

        See Also
        --------
        plot, prepare

        Notes
        -----
        This method is assumed to be over-ridden in derived classes,
        where it will perform the statistic calculations needed to
        create the visualization. This version should be called from
        these classes as it validates the fit and par arguments.

        """

        if isinstance(fit.stat, LeastSq):
            raise ConfidenceErr('badargconf', fit.stat.name)

        # Check the parameter
        #  - is thawed
        #  - is part of the model expression
        #
        if par.frozen:
            raise ConfidenceErr('frozen', par.fullname,
                                f'interval {self.conf_type}')

        if par not in fit.model.pars:
            raise ConfidenceErr('thawed', par.fullname, fit.model.name)

        # Make sure xerr and yerr are cleared
        self.xerr = None
        self.yerr = None

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        if self.log:
            self.plot_prefs['xlog'] = True

        super().plot(overplot=overplot, clearwindow=clearwindow,
                     **kwargs)

        # Note: the user arguments are not applied to the lines
        #
        if self.stat is not None:
            super().hline(self.stat, linecolor="green",
                          linestyle="dash", linewidth=1.5,
                          overplot=True, clearwindow=False)

        if self.parval is not None:
            super().vline(self.parval, linecolor="orange",
                          linestyle="dash", linewidth=1.5,
                          overplot=True, clearwindow=False)

        if self.log:
            self.plot_prefs['xlog'] = False


class Confidence2D(DataContour, Point):
    """The base class for 2D confidence contours.

    .. versionchanged:: 4.16.1
       Handling of log-scaled axes and use of the delv argument has
       been improved.

    """

    _fields: list[str] = ["x0", "x1", "y", "min!", "max!", "nloop!",
                          "fac!", "delv!", "log!", "sigma!", "parval0!",
                          "parval1!", "levels!"]
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    contour_prefs = basicbackend.get_confid_contour_defaults()
    point_prefs = basicbackend.get_confid_point_defaults()

    # This is only used for an error message
    conf_type = "unknown"
    "The type of confidence analysis."

    def __init__(self):
        self.min = None
        self.max = None
        self.nloop = (10, 10)
        self.fac = 4
        self.delv = None
        self.log = (False, False)
        self.sigma = (1, 2, 3)
        self.parval0 = None
        self.parval1 = None
        self.stat = None
        self.numcores = None
        super().__init__()

    def __setstate__(self, state):
        self.__dict__.update(state)

        if 'stat' not in state:
            self.__dict__['stat'] = None

        if 'numcores' not in state:
            self.__dict__['numcores'] = None

    def __str__(self) -> str:
        return display_fields(self, self._fields)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the confidence 2D plot."""
        return backend.as_html_contour2d(self)

    def prepare(self, min=None, max=None, nloop=(10, 10),
                delv=None, fac=4, log=(False, False),
                sigma=(1, 2, 3), levels=None, numcores=None):
        """Set the data to plot.

        This defines the ranges over which the statistic will be
        calculated, but does not perform the evaluation.

        Parameters
        ----------
        min, max : sequence of number or None, optional
            The minimum and maximum parameter values to used. If set
            then they must contain two elements, and if not then the
            range is calculated using the fac parameter.
        nloop : sequence of int, optional
            The number of points at which to evaluate the statistic,
            where each value must be greater than 1. This is used when
            delv is set to None.
        delv : sequence of number or None, optional
            The spacing of the parameter grids, and if set it must
            contain two values each greater than 0. This takes
            precedence over nloop.
        fac : number, optional
            Used when either min or max are not set. The parameter
            range in this case is taken to be fac times the separation
            of the covariance limits for the parameter (unless
            explicitly given).
        log : sequence of bool, optional
            Should each parameter be evaluated on a
            logarithmically-spaced grid rather than a linearly-spaced
            one?
        sigma : sequence of number, optional
            The sigma values at which to draw contours. This is only
            used if levels is set to None.
        levels : sequence of number or None, optional
            The levels at which the contours are drawn. This over-rides
            the sigma setting.
        numcores : int or None, optional
            Should the parameter evaluation use multiple CPU cores if
            available?

        See Also
        --------
        calc

        """

        self.min = min
        self.max = max
        self.nloop = nloop
        self.delv = delv
        self.fac = fac
        self.log = log
        self.sigma = sigma
        self.levels = levels
        self.parval0 = None
        self.parval1 = None
        self.numcores = numcores

    def _region_init(self, fit, par0, par1):
        """Calculate the grid to use for the parameters.

        Parameters
        ----------
        fit : sherpa.fit.Fit instance
            The current fit.
        par0, par1 : sherpa.models.parameter.Parameter instance
            The parameters to analyze.

        Returns
        -------
        x : 2D ndarray
            The parameter values to use. Unlike the Confidence1D case
            this is not just self.x0 or self.x1, but is a 2D array
            where each row represent each pair of parameters, so the
            first column is par0 and the second column is par1.

        """

        def check2(value, pname, opt=True):
            """Check value is None (when opt set) or has size of 2"""

            if opt and value is None:
                return None

            if np.isscalar(value):
                raise ConfidenceErr('badarg', pname, 'a list')

            if len(value) != 2:
                raise ConfidenceErr('badarg', pname, 'a list of size 2')

            # Ensure we return an ndarray
            return np.asarray(value)

        # Issue #1093 points out that if min or max is a tuple we can
        # have a problem, as below the code can assign to an element
        # of one of them. So ensure we have a list not a tuple. The
        # exact setting of min/max (whether a tuple, list, or ndarray)
        # makes the following code a bit messy.
        #
        # See also #1967 which discusses what type we should accept
        # for fields like min and max.
        #
        if self.min is not None:
            try:
                self.min = list(self.min)
            except TypeError:
                raise ConfidenceErr(
                    'badarg', 'Parameter limits', 'a list') from None

        if self.max is not None:
            try:
                self.max = list(self.max)
            except TypeError:
                raise ConfidenceErr(
                    'badarg', 'Parameter limits', 'a list') from None

        # We ignore the return value for min/max
        check2(self.min, "Parameter limits")
        check2(self.max, "Parameter limits")
        nloop = check2(self.nloop, "Nloop parameter", opt=False)
        delv = check2(self.delv, "delv parameter")

        if np.any(nloop <= 1):
            raise ConfidenceErr('badarg', 'Nloop parameter',
                                'a list with elements > 1')

        if delv is not None and np.any(delv <= 0):
            raise ConfidenceErr('badarg', 'delv parameter',
                                'a list with elements > 0')

        self.stat = fit.calc_stat()
        self.xlabel = par0.fullname
        self.ylabel = par1.fullname
        self.parval0 = par0.val
        self.parval1 = par1.val

        if self.levels is None:
            stat = self.stat
            if self.sigma is None or np.isscalar(self.sigma):
                raise ConfidenceErr('needlist', 'sigma bounds')

            sigma = np.asarray(self.sigma)
            lvls = stat - (2. * np.log(1. - erf(sigma / np.sqrt(2.))))
            self.levels = np.asarray(lvls, dtype=SherpaFloat)

        if self.min is None or self.max is None:
            oldestmethod = fit.estmethod
            fit.estmethod = Covariance()

            r = fit.est_errors()
            fit.estmethod = oldestmethod

            index0 = list(r.parnames).index(par0.fullname)
            index1 = list(r.parnames).index(par1.fullname)

            if self.min is None:
                self.min = np.array([par0.min, par1.min])
                min0 = r.parmins[index0]
                min1 = r.parmins[index1]

                if min0 is not None and not np.isnan(min0):
                    self.min[0] = par0.val + min0

                if min1 is not None and not np.isnan(min1):
                    self.min[1] = par1.val + min1

            if self.max is None:
                self.max = np.array([par0.max, par1.max])
                max0 = r.parmaxes[index0]
                max1 = r.parmaxes[index1]

                if max0 is not None and not np.isnan(max0):
                    self.max[0] = par0.val + max0

                if max1 is not None and not np.isnan(max1):
                    self.max[1] = par1.val + max1

            # This assumes that self.min/max are ndarray
            v = (self.max + self.min) / 2.
            dv = np.fabs(v - self.min)
            self.min = v - self.fac * dv
            self.max = v + self.fac * dv

        hmin = np.array([par0.min, par1.min])
        hmax = np.array([par0.max, par1.max])

        for i in [0, 1]:
            # check user limits for errors
            if self.min[i] >= self.max[i]:
                raise ConfidenceErr('badlimits')

            if self.min[i] < hmin[i]:
                self.min[i] = hmin[i]
            if self.max[i] > hmax[i]:
                self.max[i] = hmax[i]

        if delv is None:
            # Just make the following code easier to write
            delv = [None, None]

        x0 = calc_par_range(self.min[0], self.max[0],
                            self.nloop[0], delv=delv[0],
                            log=self.log[0])
        x1 = calc_par_range(self.min[1], self.max[1],
                            self.nloop[1], delv=delv[1],
                            log=self.log[1])

        # Does it matter whether this is x0,x1 or x1,x0? Not for the
        # calculation of the parameter values, but the contour
        # plotting may care (but at least self.x0 and self.x1 should
        # match the ordering of self.y).
        #
        x0g, x1g = np.meshgrid(x0, x1)
        self.x0 = x0g.flatten()
        self.x1 = x1g.flatten()

        return np.array([self.x0, self.x1]).T

    def calc(self, fit, par0, par1):
        """Evaluate the statistic for the parameter range.

        This requires prepare to have been called, and must be
        called before contour is called.

        Parameters
        ----------
        fit
            The Sherpa fit instance to use (defines the statistic
            and optimiser to use).
        par0, par1
            The parameters to iterate over.

        See Also
        --------
        contour, prepare

        Notes
        -----
        This method is assumed to be over-ridden in derived classes,
        where it will perform the statistic calculations needed to
        create the visualization. This version should be called from
        these classes as it validates the fit and par arguments.

        """

        if isinstance(fit.stat, LeastSq):
            raise ConfidenceErr('badargconf', fit.stat.name)

        # Check that each parameter
        #  - is thawed
        #  - is part of the model expression
        #
        for par in [par0, par1]:
            if par.frozen:
                raise ConfidenceErr('frozen', par.fullname,
                                    f'region {self.conf_type}')

            if par not in fit.model.pars:
                raise ConfidenceErr('thawed', par.fullname, fit.model.name)

    # TODO: should this be overcontour rather than overplot?
    def contour(self, overplot=False, clearwindow=True, **kwargs):

        if self.log[0]:
            self.contour_prefs['xlog'] = True
        if self.log[1]:
            self.contour_prefs['ylog'] = True

        # Work around possible naming issue with overplot/overcontour.
        # Set overcontour if either overcontour or overplot are set.
        # The overcontour argument must be removed from kwargs if
        # sent, hence it us used first here (if second lazy evaluation
        # could mean the pop is never executed).
        #
        overcontour = kwargs.pop('overcontour', False) or overplot

        super().contour(overcontour=overcontour,
                        clearwindow=clearwindow, **kwargs)

        # Note: the user arguments are not applied to the point
        #
        point = Point()
        point.point_prefs = self.point_prefs
        point.point(self.parval0, self.parval1,
                    overplot=True, clearwindow=False)

        if self.log[0]:
            self.contour_prefs['xlog'] = False
        if self.log[1]:
            self.contour_prefs['ylog'] = False


class IntervalProjectionWorker:
    """Used to evaluate the model by IntervalProjection.

    .. versionchanged:: 4.16.1
       The calling convention was changed.

    See Also
    --------
    IntervalUncertaintyWorker

    """

    def __init__(self, par, fit, otherpars):
        self.par = par
        self.fit = fit
        self.otherpars = otherpars

    def __call__(self, val):
        self.par.val = val

        # It there are other parameters then we need to fit to get the
        # best-fit statistic.
        #
        if self.otherpars:
            r = self.fit.fit()
            return r.statval

        # If this was the only free parameter we can just calculate
        # the statistic value.
        #
        return self.fit.calc_stat()


class IntervalProjection(Confidence1D):
    """The Interval-Projection method.

    Evaluate the parameter value on a grid of points, allowing the
    other thawed parameters to be fit.

    .. versionchanged:: 4.16.1
       Handling of log-scaled axes has been improved and the string
       output now includes the parameter value (if available).

    See Also
    --------
    IntervalUncertainty

    """

    conf_type = "projection"

    def __init__(self):
        self.fast = True
        super().__init__()

    def prepare(self, fast=True, min=None, max=None, nloop=20,
                delv=None, fac=1, log=False, numcores=None):
        self.fast = fast
        super().prepare(min, max, nloop, delv, fac, log, numcores)

    def calc(self, fit, par, methoddict=None, cache=True):
        self.title = 'Interval-Projection'
        super().calc(fit=fit, par=par)

        # If "fast" option enabled, set fitting method to
        # lmdif if stat is chi-squared,
        # else set to neldermead.

        # If current method is not LM or NM, warn it is not a good
        # method for estimating parameter limits.
        #
        if not isinstance(fit.method, (NelderMead, LevMar)):
            warning("%s is inappropriate for confidence limit "
                    "estimation", fit.method.name)

        oldfitmethod = fit.method
        if bool_cast(self.fast) is True and methoddict is not None:
            if isinstance(fit.stat, Likelihood):
                if not isinstance(fit.method, NelderMead):
                    fit.method = methoddict['neldermead']
                    warning("Setting optimization to %s for interval "
                            "projection plot", fit.method.name)

            elif not isinstance(fit.method, LevMar):
                fit.method = methoddict['levmar']
                warning("Setting optimization to %s for interval "
                        "projection plot", fit.method.name)

        xvals = self._interval_init(fit, par)
        oldpars = fit.model.thawedpars
        par.freeze()

        # We know that par is thawed, so we can check to see whether a fit
        # is needed by looking for other parameters.
        #
        otherpars = len(oldpars) > 1

        # Store these before we enter the try block as they are used
        # in the finally block.
        #
        startup = fit.model.startup
        teardown = fit.model.teardown

        try:
            fit.model.startup(cache)

            # these calls are unnecessary for every fit
            fit.model.startup = return_none
            fit.model.teardown = return_none

            worker = IntervalProjectionWorker(par, fit, otherpars)
            res = parallel_map(worker, xvals, self.numcores)
            self.y = np.asarray(res)

        finally:
            # Set back data that we changed
            par.thaw()

            fit.model.startup = startup
            fit.model.teardown = teardown

            fit.model.teardown()
            fit.model.thawedpars = oldpars
            fit.method = oldfitmethod


class IntervalUncertaintyWorker:
    """Used to evaluate the model by IntervalUncertainty.

    .. versionchanged:: 4.16.1
       The calling convention was changed.

    See Also
    --------
    IntervalProjectionWorker

    """

    def __init__(self, par, fit):
        self.par = par
        self.fit = fit

    def __call__(self, val):
        self.par.val = val
        return self.fit.calc_stat()


class IntervalUncertainty(Confidence1D):
    """The Interval-Projection method.

    Evaluate the parameter value on a grid of points, where the other
    thawed parameters are *not* changed from their current values.

    .. versionchanged:: 4.16.1
       Handling of log-scaled axes has been improved and the string
       output now includes the parameter value (if available).

    See Also
    --------
    IntervalProjection

    """

    conf_type = "uncertainty"

    def calc(self, fit, par, methoddict=None, cache=True):
        self.title = 'Interval-Uncertainty'
        super().calc(fit=fit, par=par)

        thawed = [p for p in fit.model.pars if not p.frozen]
        oldpars = fit.model.thawedpars
        xvals = self._interval_init(fit, par)
        for p in thawed:
            p.freeze()

        try:
            fit.model.startup(cache)

            worker = IntervalUncertaintyWorker(par, fit)
            res = parallel_map(worker, xvals, self.numcores)
            self.y = np.asarray(res)

        finally:
            # Set back data that we changed
            for p in thawed:
                p.thaw()

            fit.model.teardown()
            fit.model.thawedpars = oldpars


class RegionProjectionWorker:
    """Used to evaluate the model by RegionProjection.

    .. versionchanged:: 4.16.1
       The calling convention was changed.

    See Also
    --------
    RegionUncertaintyWorker

    """

    def __init__(self, par0, par1, fit, otherpars):
        self.par0 = par0
        self.par1 = par1
        self.fit = fit
        self.otherpars = otherpars

    def __call__(self, pars):
        (self.par0.val, self.par1.val) = pars

        # It there are other parameters then we need to fit to get the
        # best-fit statistic.
        #
        if self.otherpars:
            r = self.fit.fit()
            return r.statval

        # If these were the only two free parameters we can just
        # calculate the statistic value.
        #
        return self.fit.calc_stat()


def return_none(cache=None):
    """dummy implementation of callback for multiprocessing"""
    return None


class RegionProjection(Confidence2D):
    """The Region-Projection method.

    Evaluate the statistic on a grid of points for two parameters,
    where the other thawed parameters are fit for each location.

    .. versionchanged:: 4.16.1
       Support for logarithmically-spaced grids has been improved.

    See Also
    --------
    RegionUncertainty

    """

    conf_type = "projection"

    def __init__(self):
        self.fast = True
        super().__init__()

    def prepare(self, fast=True, min=None, max=None, nloop=(10, 10),
                delv=None, fac=4, log=(False, False),
                sigma=(1, 2, 3), levels=None, numcores=None):
        self.fast = fast
        super().prepare(min, max, nloop, delv, fac, log, sigma,
                        levels=levels, numcores=numcores)

    def calc(self, fit, par0, par1, methoddict=None, cache=True):
        self.title = 'Region-Projection'
        super().calc(fit=fit, par0=par0, par1=par1)

        # If "fast" option enabled, set fitting method to
        # lmdif if stat is chi-squared,
        # else set to neldermead

        # If current method is not LM or NM, warn it is not a good
        # method for estimating parameter limits.
        #
        if not isinstance(fit.method, (NelderMead, LevMar)):
            warning("%s is inappropriate for confidence limit "
                    "estimation", fit.method.name)

        oldfitmethod = fit.method
        if bool_cast(self.fast) is True and methoddict is not None:
            if isinstance(fit.stat, Likelihood):
                if not isinstance(fit.method, NelderMead):
                    fit.method = methoddict['neldermead']
                    warning("Setting optimization to %s for region "
                            "projection plot", fit.method.name)

            elif not isinstance(fit.method, LevMar):
                fit.method = methoddict['levmar']
                warning("Setting optimization to %s for region "
                        "projection plot", fit.method.name)

        oldpars = fit.model.thawedpars
        otherpars = len(oldpars) > 2

        # Store these before we enter the try block as they are used
        # in the finally block.
        #
        startup = fit.model.startup
        teardown = fit.model.teardown

        try:
            fit.model.startup(cache)

            # these calls are unnecessary for every fit
            fit.model.startup = return_none
            fit.model.teardown = return_none

            grid = self._region_init(fit, par0, par1)

            par0.freeze()
            par1.freeze()

            worker = RegionProjectionWorker(par0, par1, fit, otherpars)
            results = parallel_map(worker, grid, self.numcores)
            self.y = np.asarray(results)

        finally:
            # Set back data after we changed it
            par0.thaw()
            par1.thaw()

            fit.model.startup = startup
            fit.model.teardown = teardown

            fit.model.teardown()
            fit.model.thawedpars = oldpars
            fit.method = oldfitmethod


class RegionUncertaintyWorker:
    """Used to evaluate the model by RegionUncertainty.

    .. versionchanged:: 4.16.1
       The calling convention was changed.

    See Also
    --------
    RegionProjectionWorker

    """

    def __init__(self, par0, par1, fit):
        self.par0 = par0
        self.par1 = par1
        self.fit = fit

    def __call__(self, pars):
        (self.par0.val, self.par1.val) = pars
        return self.fit.calc_stat()


class RegionUncertainty(Confidence2D):
    """The Region-Projection method.

    Evaluate the statistic on a grid of points for two parameters,
    where the other thawed parameters are *not* changed from their
    current values.

    .. versionchanged:: 4.16.1
       Support for logarithmically-spaced grids has been improved.

    See Also
    --------
    RegionProjection

    """

    conf_type = "uncertainty"

    def calc(self, fit, par0, par1, methoddict=None, cache=True):
        self.title = 'Region-Uncertainty'
        super().calc(fit=fit, par0=par0, par1=par1)

        thawed = [i for i in fit.model.pars if not i.frozen]
        oldpars = fit.model.thawedpars

        try:
            fit.model.startup(cache)

            grid = self._region_init(fit, par0, par1)

            for p in thawed:
                p.freeze()

            worker = RegionUncertaintyWorker(par0, par1, fit)
            result = parallel_map(worker, grid, self.numcores)
            self.y = np.asarray(result)

        finally:
            # Set back data after we changed it
            for p in thawed:
                p.thaw()

            fit.model.teardown()
            fit.model.thawedpars = oldpars


def get_per_plot_kwargs(nplots: int,
                        kwargs
                        ) -> list[dict[str, Any]]:
    """Separate the per-plot plot arguments.

    .. versionadded:: 4.17.0

    Parameters
    ----------
    nplots
       The number of plots.
    kwargs
       The user-specified keyword arguments. Any values which are
       sequences (not a string) need to contain nplots elements.

    Returns
    -------
    kwstores
       The keyword arguments for each plot.

    """

    # Allow kwargs to be specified per-plot. This is done by
    # checking if any value is an iterable (and not a string) and
    # extracting a single value per plot.
    #
    # Need these {} to be separate which means we can not just say
    # `kwstore = [{}] * nplots`.
    #
    ## kwstore: list[dict[str, Any]]
    kwstore = [{} for _ in range(nplots)]
    for key, val in kwargs.items():
        if is_iterable_not_str(val):
            nval = len(val)
            if nval != nplots:
                raise ValueError(f"keyword '{key}': expected "
                                 f"{nplots} elements but found "
                                 f"{nval}")

            for store, v in zip(kwstore, val):
                store[key] = v

        else:
            for store in kwstore:
                store[key] = val

    return kwstore


class MultiPlot:
    """Combine multiple line-style plots.

    Allow multiple line-style plots - so those plot classes
    that use the `plot` method to display - to be drawn in
    the same area. Each plot is added with the add method.

    .. versionadded:: 4.16.1

    """

    __slots__ = ("plots", "title")

    def __init__(self) -> None:
        self.plots: list[Union[Plot, HistogramPlot]] = []
        self.title = ""

    # The typing here says Plot but we actually want the sub-classes
    # like DataPlot (i.e.  those that use the prepare method to set up
    # the data to plot so that the plot method requires no data
    # arguments) rather than the actual Plot class.
    #
    def add(self, plot: Union[Plot, HistogramPlot]) -> None:
        """Add the plot to the list of data to plot.

        A copy of the plot object is stored, rather than the
        input argument. The `title` attribute can be set or
        changed.

        Parameters
        ----------
        plot : instance
           The plot or histogram object to add. It must have
           a `plot` method.

        """

        # A copy is stored since often the same underlying object is
        # returned by the UI routines. It also allows the code to
        # change the self.plots array and not worry about changing the
        # calling code.
        #
        self.plots.append(copy.deepcopy(plot))

    def plot(self,
             overplot: bool = False,
             clearwindow: bool = True,
             **kwargs) -> None:
        """Plot the data.

        .. versionchanged:: 4.17.0
           The keyword arguments can now be set per plot by sending in
           a sequence of values.

        Parameters
        ----------
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's plot_prefs
           dictionary. If the value is a scalar then the value is sent
           to each plot, but if it's a sequence (not a string) then it
           must match the number of plots and the values are used
           sequentially.

        """

        kwstore = get_per_plot_kwargs(len(self.plots),
                                      kwargs)
        for plot, store in zip(self.plots, kwstore):
            plot.plot(overplot=overplot, clearwindow=clearwindow,
                      **store)

            # To decide whether we draw a title or not we do rely on
            # the clearwindow setting, since it is not guaranteed to
            # be set (e.g. if this is being displayed within a
            # SplitPlot).  So the overplot setting is used to decide
            # whether to add in the title or not.
            #
            if not overplot:
                backend.set_title(self.title)

            overplot = True
            clearwindow = False
