#
#  Copyright (C) 2009, 2015, 2016, 2018, 2019, 2020, 2021, 2022, 2023
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
Sherpa configuration file. Note that plot objects can be created
and used even when there is no available plot backend. It is just
that no graphical display will be created.

Which backend is used?
----------------------

When this module is first imported, Sherpa tries to import the
backends installed with Sherpa in the order listed in the
``options.plot_pkg`` setting from the ``sherpa.rc`` startup file.
The first module that imports successfully is set as the active
backend. The following command prints the name and the location
on disk of that module::

   >>> from sherpa import plot
   >>> print(plot.backend)

Change the backend
------------------

After the initial import, the backend can be changed by loading one of
the plotting backends shipped with sherpa (or any other module that
provides the same interface):

  >>> import sherpa.plot.pylab_backend
  >>> plot.backend = sherpa.plot.pylab_backend

Creating a plot
---------------

The simplest approach is to use the `sherpa_plot` context manager::

    with sherpa_plot():
        # Now call the plot/overplot or contour/overcontour methods
        obj.plot()

This handles setting up the backend, handles any error handling,
and then ends the session. It is equivalent to::

    from sherpa import plot
    try:
        plot.backend.begin()
        # Now call the plot/overplot or contour/overcontour methods
        obj.plot()
    except:
        plot.backend.exceptions()
        raise
    else:
        plot.backend.end()

"""

from configparser import ConfigParser
from contextlib import contextmanager
import logging
import importlib

import numpy

from sherpa.utils import NoNewAttributesAfterInit, erf, SherpaFloat, \
    bool_cast, parallel_map, dataspace1d, histogram1d, get_error_estimates
from sherpa.utils.err import PlotErr, StatErr, ConfidenceErr
from sherpa.estmethods import Covariance
from sherpa.optmethods import LevMar, NelderMead
from sherpa.stats import Likelihood, LeastSq, Chi2XspecVar
from sherpa import get_config
config = ConfigParser()
config.read(get_config())

lgr = logging.getLogger(__name__)
warning = lgr.warning

# TODO: why is this module globally changing the invalid mode of NumPy?
_ = numpy.seterr(invalid='ignore')

plot_opt = config.get('options', 'plot_pkg', fallback='dummy')
plot_opt = [o.strip().lower() + '_backend' for o in plot_opt.split()]

backend = None
'''Currently active backend module for plotting.'''

for plottry in plot_opt:
    try:
        backend = importlib.import_module('.' + plottry,
                                          package='sherpa.plot')
        break
    except ImportError:
        pass
else:
    # None of the options in the rc file work, e.g. because it's an old file
    # that does not have dummy listed
    import sherpa.plot.dummy_backend as backend

__all__ = ('Plot', 'Contour', 'Point', 'Histogram',
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
           'DelchiPlot', 'ChisqrPlot', 'ResidPlot',
           'ResidContour',
           'RatioPlot',
           'RatioContour',
           'Confidence1D', 'Confidence2D',
           'IntervalProjection', 'IntervalUncertainty',
           'RegionProjection', 'RegionUncertainty',
           'sherpa_plot'
           )

_stats_noerr = ('cash', 'cstat', 'leastsq', 'wstat')


def _make_title(title, name=''):
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
    else:
        return "{} for {}".format(title, name)


def _errorbar_warning(stat):
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
        "used in fits with {}".format(stat.name)


def calculate_errors(data, stat, yerrorbars=True):
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
            warning(msg + "\nzeros or negative values found")

        return None


@contextmanager
def sherpa_plot():
    """Context manager to handle Sherpa plot actions with a backend.

    Handles setting up the backend and then handle any errors and
    clean up.

    .. versionadded:: 4.15.1

    Examples
    --------

    Plot the data in the dplot object:

    >>> with sherpa_plot():
    ...     dplot.plot()

    """

    try:
        backend.begin()
        yield

    except:
        backend.exceptions()
        raise
    else:
        backend.end()


class Plot(NoNewAttributesAfterInit):
    "Base class for line plots"
    plot_prefs = backend.get_plot_defaults()

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
        NoNewAttributesAfterInit.__init__(self)

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
    contour_prefs = backend.get_contour_defaults()

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
        NoNewAttributesAfterInit.__init__(self)

    def _merge_settings(self, kwargs):
        """Return the plot preferences merged with user settings."""
        return {**self.contour_prefs, **kwargs}

    def contour(self, x0, x1, y, levels=None, title=None, xlabel=None,
                ylabel=None, overcontour=False, clearwindow=True,
                **kwargs):
        opts = self._merge_settings(kwargs)
        backend.contour(x0, x1, y, levels, title, xlabel, ylabel, overcontour,
                        clearwindow, **opts)

    def overcontour(self, *args, **kwargs):
        kwargs['overcontour'] = True
        self.contour(*args, **kwargs)


class Point(NoNewAttributesAfterInit):
    "Base class for point plots"
    point_prefs = backend.get_point_defaults()

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
        NoNewAttributesAfterInit.__init__(self)

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
        backend.point(x, y,
                      overplot=overplot, clearwindow=clearwindow,
                      **opts)


class Histogram(NoNewAttributesAfterInit):
    "Base class for histogram plots"
    histo_prefs = backend.get_histo_defaults()

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
        NoNewAttributesAfterInit.__init__(self)

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

    def __init__(self):
        self.xlo = None
        self.xhi = None
        self.y = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        Histogram.__init__(self)

    def __str__(self):
        xlo = self.xlo
        if self.xlo is not None:
            xlo = numpy.array2string(numpy.asarray(self.xlo), separator=',',
                                     precision=4, suppress_small=False)

        xhi = self.xhi
        if self.xhi is not None:
            xhi = numpy.array2string(numpy.asarray(self.xhi), separator=',',
                                     precision=4, suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(numpy.asarray(self.y), separator=',',
                                   precision=4, suppress_small=False)

        return (('xlo    = %s\n' +
                 'xhi    = %s\n' +
                 'y      = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'histo_prefs = %s') %
                (xlo,
                 xhi,
                 y,
                 self.xlabel,
                 self.ylabel,
                 self.title,
                 self.histo_prefs))

    def _repr_html_(self):
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
        xlo = numpy.asarray(self.xlo)
        xhi = numpy.asarray(self.xhi)
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

        Histogram.plot(self, self.xlo, self.xhi, self.y, title=self.title,
                       xlabel=self.xlabel, ylabel=self.ylabel,
                       overplot=overplot, clearwindow=clearwindow, **kwargs)


# I think we want slightly different histogram preferences
# than most (mark by point rather than line).
#
def get_data_hist_prefs():
    """Copy the data preferences to the histogram class"""

    hprefs = backend.get_model_histo_defaults()
    dprefs = backend.get_data_plot_defaults()
    for k, v in dprefs.items():
        if k not in hprefs:
            continue

        hprefs[k] = v

    return hprefs


class DataHistogramPlot(HistogramPlot):
    """Create 1D histogram plots of data values."""

    histo_prefs = get_data_hist_prefs()

    def __init__(self):
        self.xerr = None
        self.yerr = None
        super().__init__()

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
        (_, self.y, self.yerr, self.xerr, self.xlabel,
         self.ylabel) = data.to_plot()

        self.xlo, self.xhi = data.get_indep(True)

        if stat is not None:
            yerrorbars = self.histo_prefs.get('yerrorbars', True)
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

        Histogram.plot(self, self.xlo, self.xhi, self.y,
                       yerr=self.yerr, title=self.title,
                       xlabel=self.xlabel, ylabel=self.ylabel,
                       overplot=overplot, clearwindow=clearwindow, **kwargs)


class ModelHistogramPlot(HistogramPlot):
    """Display a model as a histogram."""

    def __init__(self):
        super().__init__()
        self.title = 'Model'

    def prepare(self, data, model, stat=None):
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
        assert self.y.size == self.xlo.size


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

    def __init__(self):
        self.points = None
        HistogramPlot.__init__(self)

    def __str__(self):
        points = self.points
        if self.points is not None:
            points = numpy.array2string(numpy.asarray(self.points),
                                        separator=',', precision=4,
                                        suppress_small=False)

        return ('points = %s\n' % (points) + HistogramPlot.__str__(self))

    def _repr_html_(self):
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
        self.y, xx = numpy.histogram(points, bins=bins, density=normed)
        self.xlo = xx[:-1]
        self.xhi = xx[1:]
        self.ylabel = "probability density"
        self.xlabel = xlabel
        self.title = "PDF: {}".format(name)


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

    >>> rng = np.random.default_rng()  # requires NumPy 1.17 or later
    >>> pts = rng.gumbel(loc=100, scale=20, size=1000)
    >>> plot = CDFPlot()
    >>> plot.prepare(pts)
    >>> plot.plot()

    """

    median_defaults = dict(linestyle='dash', linecolor='orange',
                           linewidth=1.5)
    """The options used to draw the median line."""

    lower_defaults = dict(linestyle='dash', linecolor='blue',
                          linewidth=1.5)
    """The options used to draw the 15.87% line."""

    upper_defaults = dict(linestyle='dash', linecolor='blue',
                          linewidth=1.5)
    """The options used to draw the 84.13% line."""

    plot_prefs = backend.get_cdf_plot_defaults()
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
        Plot.__init__(self)

    def __str__(self):
        x = self.x
        if self.x is not None:
            x = numpy.array2string(self.x, separator=',', precision=4,
                                   suppress_small=False)
        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4,
                                   suppress_small=False)

        points = self.points
        if self.points is not None:
            points = numpy.array2string(self.points, separator=',',
                                        precision=4, suppress_small=False)

        return f"""points = {points}
x      = {x}
y      = {y}
median = {self.median}
lower  = {self.lower}
upper  = {self.upper}
xlabel = {self.xlabel}
ylabel = {self.ylabel}
title  = {self.title}
plot_prefs = {self.plot_prefs}"""

    def _repr_html_(self):
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

        self.points = numpy.asarray(points)
        self.x = numpy.sort(points)
        (self.median, self.lower,
         self.upper) = get_error_estimates(self.x, True)
        xsize = len(self.x)
        self.y = (numpy.arange(xsize) + 1.0) / xsize
        self.xlabel = xlabel
        self.ylabel = "p(<={})".format(xlabel)
        self.title = "CDF: {}".format(name)

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

        Plot.plot(self, self.x, self.y, title=self.title,
                  xlabel=self.xlabel, ylabel=self.ylabel,
                  overplot=overplot, clearwindow=clearwindow, **kwargs)

        # Note: the user arguments are not applied to the vertical lines
        #
        Plot.vline(self.median, overplot=True, clearwindow=False,
                   **self.median_defaults)
        Plot.vline(self.lower, overplot=True, clearwindow=False,
                   **self.lower_defaults)
        Plot.vline(self.upper, overplot=True, clearwindow=False,
                   **self.upper_defaults)


class LRHistogram(HistogramPlot):
    "Derived class for creating 1D likelihood ratio distribution plots"
    def __init__(self):
        self.ratios = None
        self.lr = None
        self.ppp = None
        HistogramPlot.__init__(self)

    def __str__(self):
        ratios = self.ratios
        if self.ratios is not None:
            ratios = numpy.array2string(numpy.asarray(self.ratios),
                                        separator=',', precision=4,
                                        suppress_small=False)

        return '\n'.join(['ratios = %s' % ratios,
                          'lr = %s' % str(self.lr),
                          HistogramPlot.__str__(self)])

    def _repr_html_(self):
        """Return a HTML (string) representation of the LRHistogram plot."""
        return backend.as_html_lr(self)

    def prepare(self, ratios, bins, niter, lr, ppp):
        "Create the data to plot"

        self.ppp = float(ppp)
        self.lr = float(lr)
        y = numpy.asarray(ratios)
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

        Histogram.plot(self, self.xlo, self.xhi, self.y, title=self.title,
                       xlabel=self.xlabel, ylabel=self.ylabel,
                       overplot=overplot, clearwindow=clearwindow, **kwargs)

        if self.lr is None:
            return

        # Note: the user arguments are not applied to the vertical line
        #
        if self.lr <= self.xhi.max() and self.lr >= self.xlo.min():
            Plot.vline(self.lr, linecolor="orange", linestyle="solid",
                       linewidth=1.5, overplot=True, clearwindow=False)


class SplitPlot(Plot, Contour):
    """Create multiple plots.

    Attributes
    ----------
    rows : int
       Number of rows of plots. The default is 2.
    cols : int
       Number of columns of plots. The default is 1.
    plot_prefs : dict
       The preferences for the plots. This depends on the plot
       backend.

    """
    plot_prefs = backend.get_split_plot_defaults()

    def __init__(self, rows=2, cols=1):
        self.reset(rows, cols)
        Plot.__init__(self)
        Contour.__init__(self)

    def __str__(self):
        return (('rows   = %s\n' +
                 'cols   = %s\n' +
                 'plot_prefs = %s') %
                (self.rows,
                 self.cols,
                 self.plot_prefs))

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
        self._used = numpy.zeros((self.rows, self.cols), numpy.bool_)

    def _next_subplot(self):
        row, col = numpy.where(~self._used)
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

    def __init__(self):
        SplitPlot.__init__(self)

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
    plot_prefs : dict
       The preferences for the plot.
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

    plot_prefs = backend.get_data_plot_defaults()

    def __init__(self):
        self.x = None
        self.y = None
        self.yerr = None
        self.xerr = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        Plot.__init__(self)

    def __str__(self):
        x = self.x
        if self.x is not None:
            x = numpy.array2string(self.x, separator=',', precision=4,
                                   suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4,
                                   suppress_small=False)

        yerr = self.yerr
        if self.yerr is not None:
            yerr = numpy.array2string(self.yerr, separator=',', precision=4,
                                      suppress_small=False)

        xerr = self.xerr
        if self.xerr is not None:
            xerr = numpy.array2string(self.xerr, separator=',', precision=4,
                                      suppress_small=False)

        return (('x      = %s\n' +
                 'y      = %s\n' +
                 'yerr   = %s\n' +
                 'xerr   = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'plot_prefs = %s') %
                (x,
                 y,
                 yerr,
                 xerr,
                 self.xlabel,
                 self.ylabel,
                 self.title,
                 self.plot_prefs))

    def _repr_html_(self):
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
        Plot.plot(self, self.x, self.y, yerr=self.yerr, xerr=self.xerr,
                  title=self.title, xlabel=self.xlabel, ylabel=self.ylabel,
                  overplot=overplot, clearwindow=clearwindow, **kwargs)


class TracePlot(DataPlot):

    plot_prefs = backend.get_model_plot_defaults()

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
        self.x = numpy.arange(len(points), dtype=SherpaFloat)
        self.y = points
        self.xlabel = "iteration"
        self.ylabel = name
        self.title = "Trace: {}".format(name)


class ScatterPlot(DataPlot):

    plot_prefs = backend.get_scatter_plot_defaults()

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
        self.x = numpy.asarray(x, dtype=SherpaFloat)
        self.y = numpy.asarray(y, dtype=SherpaFloat)
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = "Scatter: {}".format(name)


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
        DataPlot.prepare(self, psfdata, stat)
        # self.ylabel = 'PSF value'
        # self.xlabel = 'PSF Kernel size'
        self.title = 'PSF Kernel'


class DataContour(Contour):
    """Create contours of 2D data.

    Attributes
    ----------
    contour_prefs : dict
       The preferences for the plot.
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

    contour_prefs = backend.get_data_contour_defaults()

    def __init__(self):
        self.x0 = None
        self.x1 = None
        self.y = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        self.levels = None
        Contour.__init__(self)

    def __str__(self):
        x0 = self.x0
        if self.x0 is not None:
            x0 = numpy.array2string(self.x0, separator=',', precision=4,
                                    suppress_small=False)

        x1 = self.x1
        if self.x1 is not None:
            x1 = numpy.array2string(self.x1, separator=',', precision=4,
                                    suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4,
                                   suppress_small=False)

        return (('x0     = %s\n' +
                 'x1     = %s\n' +
                 'y      = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'levels = %s\n' +
                 'contour_prefs = %s') %
                (x0,
                 x1,
                 y,
                 self.xlabel,
                 self.ylabel,
                 self.title,
                 self.levels,
                 self.contour_prefs))

    def _repr_html_(self):
        """Return a HTML (string) representation of the contour plot."""
        return backend.as_html_datacontour(self)

    def prepare(self, data, stat=None):
        (self.x0, self.x1, self.y, self.xlabel,
         self.ylabel) = data.to_contour()
        self.title = data.name

    def contour(self, overcontour=False, clearwindow=True, **kwargs):
        Contour.contour(self, self.x0, self.x1, self.y,
                        self.levels, self.title, self.xlabel,
                        self.ylabel, overcontour=overcontour,
                        clearwindow=clearwindow, **kwargs)


class PSFKernelContour(DataContour):
    "Derived class for creating 2D PSF Kernel contours"

    def prepare(self, psf, data=None, stat=None):
        psfdata = psf.get_kernel(data)
        DataContour.prepare(self, psfdata)
        # self.xlabel = 'PSF Kernel size x0'
        # self.ylabel = 'PSF Kernel size x1'
        self.title = 'PSF Kernel'


class ModelPlot(Plot):
    """Create 1D plots of model values.

    Attributes
    ----------
    plot_prefs : dict
       The preferences for the plot.
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
    >>> from sherpa.models.basic import StepLot1D
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

    plot_prefs = backend.get_model_plot_defaults()

    def __init__(self):
        self.x = None
        self.y = None
        self.yerr = None
        self.xerr = None
        self.xlabel = None
        self.ylabel = None
        self.title = 'Model'
        Plot.__init__(self)

    def __str__(self):
        x = self.x
        if self.x is not None:
            x = numpy.array2string(self.x, separator=',', precision=4,
                                   suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4,
                                   suppress_small=False)

        yerr = self.yerr
        if self.yerr is not None:
            yerr = numpy.array2string(self.yerr, separator=',', precision=4,
                                      suppress_small=False)

        xerr = self.xerr
        if self.xerr is not None:
            xerr = numpy.array2string(self.xerr, separator=',', precision=4,
                                      suppress_small=False)

        return (('x      = %s\n' +
                 'y      = %s\n' +
                 'yerr   = %s\n' +
                 'xerr   = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'plot_prefs = %s') %
                (x,
                 y,
                 yerr,
                 xerr,
                 self.xlabel,
                 self.ylabel,
                 self.title,
                 self.plot_prefs))

    def _repr_html_(self):
        """Return a HTML (string) representation of the model plot."""
        return backend.as_html_model(self)

    def prepare(self, data, model, stat=None):
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
        Plot.plot(self, self.x, self.y, title=self.title, xlabel=self.xlabel,
                  ylabel=self.ylabel, overplot=overplot,
                  clearwindow=clearwindow, **kwargs)


class ComponentModelPlot(ModelPlot):

    plot_prefs = backend.get_component_plot_defaults()

    def prepare(self, data, model, stat=None):
        ModelPlot.prepare(self, data, model, stat)
        self.title = 'Model component: %s' % model.name


class ComponentModelHistogramPlot(ModelHistogramPlot):

    # Is this the correct setting?
    plot_prefs = backend.get_component_plot_defaults()

    def prepare(self, data, model, stat=None):
        super().prepare(data, model, stat)
        self.title = 'Model component: {}'.format(model.name)


class ComponentTemplateModelPlot(ComponentModelPlot):

    def prepare(self, data, model, stat=None):
        self.x = model.get_x()
        self.y = model.get_y()
        self.xlabel = data.get_xlabel()
        self.ylabel = data.get_ylabel()
        self.title = 'Model component: {}'.format(model.name)


class SourcePlot(ModelPlot):
    """Create 1D plots of unconvolved model values.

    Attributes
    ----------
    plot_prefs : dict
       The preferences for the plot.
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
        ModelPlot.__init__(self)
        self.title = 'Source'


class ComponentSourcePlot(SourcePlot):

    plot_prefs = backend.get_component_plot_defaults()

    def prepare(self, data, model, stat=None):
        (self.x, self.y, self.yerr, self.xerr,
         self.xlabel, self.ylabel) = data.to_component_plot(yfunc=model)
        self.y = self.y[1]
        self.title = 'Source model component: {}'.format(model.name)


class ComponentSourceHistogramPlot(SourceHistogramPlot):

    # Is this the correct setting?
    plot_prefs = backend.get_component_plot_defaults()

    def prepare(self, data, model, stat=None):

        plot = data.to_component_plot(yfunc=model)
        (_, y, _, _, self.xlabel, self.ylabel) = plot

        # taken from get_x from Data1DInt
        indep = data.get_evaluation_indep(filter=True, model=model)

        self.xlo = indep[0]
        self.xhi = indep[1]
        self.y = y[1]
        assert self.y.size == self.xlo.size

        self.title = 'Source model component: {}'.format(model.name)


class ComponentTemplateSourcePlot(ComponentSourcePlot):

    def prepare(self, data, model, stat=None):
        self.x = model.get_x()
        self.y = model.get_y()

        if numpy.iterable(data.mask):
            x = data.to_plot()[0]
            mask = (self.x > x.min()) & (self.x <= x.max())
            self.x = self.x[mask]
            self.y = self.y[mask]

        self.xlabel = data.get_xlabel()
        self.ylabel = data.get_ylabel()
        self.title = 'Source model component: {}'.format(model.name)


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
        DataPlot.prepare(self, psfdata, stat)
        self.title = psf.kernel.name


class ModelContour(Contour):
    "Derived class for creating 2D model contours"
    contour_prefs = backend.get_model_contour_defaults()

    def __init__(self):
        self.x0 = None
        self.x1 = None
        self.y = None
        self.xlabel = None
        self.ylabel = None
        self.title = 'Model'
        self.levels = None
        Contour.__init__(self)

    def __str__(self):
        x0 = self.x0
        if self.x0 is not None:
            x0 = numpy.array2string(self.x0, separator=',', precision=4,
                                    suppress_small=False)

        x1 = self.x1
        if self.x1 is not None:
            x1 = numpy.array2string(self.x1, separator=',', precision=4,
                                    suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4,
                                   suppress_small=False)

        return (('x0     = %s\n' +
                 'x1     = %s\n' +
                 'y      = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'levels = %s\n' +
                 'contour_prefs = %s') %
                (x0,
                 x1,
                 y,
                 self.xlabel,
                 self.ylabel,
                 self.title,
                 self.levels,
                 self.contour_prefs))

    def _repr_html_(self):
        """Return a HTML (string) representation of the model contour plot."""
        return backend.as_html_modelcontour(self)

    def prepare(self, data, model, stat):
        (self.x0, self.x1, self.y, self.xlabel,
         self.ylabel) = data.to_contour(yfunc=model)
        self.y = self.y[1]

    def contour(self, overcontour=False, clearwindow=True, **kwargs):
        Contour.contour(self, self.x0, self.x1, self.y, levels=self.levels,
                        title=self.title, xlabel=self.xlabel,
                        ylabel=self.ylabel, overcontour=overcontour,
                        clearwindow=clearwindow, **kwargs)


class PSFContour(DataContour):
    "Derived class for creating 2D PSF contours"

    def prepare(self, psf, data=None, stat=None):
        psfdata = psf.get_kernel(data, False)
        DataContour.prepare(self, psfdata)
        self.title = psf.kernel.name


class SourceContour(ModelContour):
    "Derived class for creating 2D model contours"
    def __init__(self):
        ModelContour.__init__(self)
        self.title = 'Source'


class FitPlot(Plot):
    """Combine data and model plots for 1D data.

    Attributes
    ----------
    plot_prefs : dict
       The preferences for the plot. Note that the display for the
       data and model plots are controlled by the preferences for
       the dataplot and modelplot objects, so this is currently
       unused.
    dataplot
       The Sherpa plot object used to display the data.
    modelplot
       The Sherpa plot object used to display the model.

    Examples
    --------

    If dplot and mplot are data and model plots respectively,
    such as those created in the examples for `DataPlot` and `ModelPlot`,
    then a combined plot can be created with:

    >>> from sherpa.plot import FitPlot
    >>> fplot = FitPlot()
    >>> fplot.prepare(dplot, mplot)
    >>> fplot.plot()

    Keyword arguments can be given in the `plot` call, and these
    are passed through to both the data and model plots (in the
    following example the Matplotlib backend is assumed to be in use):

    >>> fplot.plot(color='k')

    """

    plot_prefs = backend.get_fit_plot_defaults()

    def __init__(self):
        self.dataplot = None
        self.modelplot = None
        Plot.__init__(self)

    def __str__(self):
        data_title = None
        if self.dataplot is not None:
            data_title = self.dataplot.title

        model_title = None
        if self.modelplot is not None:
            model_title = self.modelplot.title

        return (('dataplot   = %s\n%s\n\nmodelplot  = %s\n%s') %
                (data_title,
                 self.dataplot,
                 model_title,
                 self.modelplot))

    def _repr_html_(self):
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
    contour_prefs = backend.get_fit_contour_defaults()

    def __init__(self):
        self.datacontour = None
        self.modelcontour = None
        Contour.__init__(self)

    def __str__(self):
        data_title = None
        if self.datacontour is not None:
            data_title = self.datacontour.title

        model_title = None
        if self.modelcontour is not None:
            model_title = self.modelcontour.title
        return (('datacontour = %s\n%s\n\nmodelcontour = %s\n%s') %
                (data_title,
                 self.datacontour,
                 model_title,
                 self.modelcontour))

    def _repr_html_(self):
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


class DelchiPlot(ModelPlot):
    """Create plots of the delta-chi value per point.

    The value of (data-model)/error is plotted for each point.

    Attributes
    ----------
    plot_prefs : dict
       The preferences for the plot.
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

    plot_prefs = backend.get_resid_plot_defaults()

    def _calc_delchi(self, ylist, staterr):
        return (ylist[0] - ylist[1]) / staterr

    def prepare(self, data, model, stat):
        (self.x, y, staterr, self.xerr,
         self.xlabel, self.ylabel) = data.to_plot(model)

        if staterr is None:
            if stat.name in _stats_noerr:
                raise StatErr('badstat', "DelchiPlot", stat.name)
            staterr = data.get_yerr(True, stat.calc_staterror)

        self.y = self._calc_delchi(y, staterr)
        self.yerr = staterr / staterr
        self.ylabel = 'Sigma'
        self.title = _make_title('Sigma Residuals', data.name)

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        self.plot_prefs['ylog'] = False
        kwargs.pop('ylog', True)
        Plot.plot(self, self.x, self.y, yerr=self.yerr, xerr=self.xerr,
                  title=self.title, xlabel=self.xlabel, ylabel=self.ylabel,
                  overplot=overplot, clearwindow=clearwindow, **kwargs)


class ChisqrPlot(ModelPlot):
    """Create plots of the chi-square value per point.

    The value of ((data-model)/error)^2 is plotted for each point.

    Attributes
    ----------
    plot_prefs : dict
       The preferences for the plot.
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
    plot_prefs = backend.get_model_plot_defaults()

    def _calc_chisqr(self, ylist, staterr):
        dy = ylist[0] - ylist[1]
        return dy * dy / (staterr * staterr)

    def prepare(self, data, model, stat):
        (self.x, y, staterr, self.xerr,
         self.xlabel, self.ylabel) = data.to_plot(model)

        # if staterr is None:
        if stat.name in _stats_noerr:
            raise StatErr('badstat', "ChisqrPlot", stat.name)
        staterr = data.get_yerr(True, stat.calc_staterror)

        self.y = self._calc_chisqr(y, staterr)
        self.ylabel = backend.get_latex_for_string(r'\chi^2')
        self.title = _make_title(backend.get_latex_for_string(r'\chi^2'), data.name)

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        Plot.plot(self, self.x, self.y, title=self.title,
                  xlabel=self.xlabel, ylabel=self.ylabel,
                  overplot=overplot, clearwindow=clearwindow, **kwargs)


class ResidPlot(ModelPlot):
    """Create plots of the residuals (data - model) per point.

    The value of (data-model) is plotted for each point.

    Attributes
    ----------
    plot_prefs : dict
       The preferences for the plot.
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

    plot_prefs = backend.get_resid_plot_defaults()

    def _calc_resid(self, ylist):
        return ylist[0] - ylist[1]

    def prepare(self, data, model, stat):
        (self.x, y, self.yerr, self.xerr,
         self.xlabel, self.ylabel) = data.to_plot(model)

        self.y = self._calc_resid(y)

        # See the discussion in DataPlot.prepare
        try:
            yerrorbars = self.plot_prefs['yerrorbars']
        except KeyError:
            yerrorbars = True

        if stat.name in _stats_noerr:
            self.yerr = data.get_yerr(True, Chi2XspecVar.calc_staterror)
            if yerrorbars:
                warning(_errorbar_warning(stat))
        else:
            self.yerr = data.get_yerr(True, stat.calc_staterror)

        # Some data sets (e.g. DataPHA, which shows the units) have a y
        # label that could (should?) be displayed (or added to the label).
        # To avoid a change in behavior, the label is only changed if
        # the "generic" Y axis label is used. To be reviewed.
        #
        if self.ylabel == 'y':
            self.ylabel = 'Data - Model'

        self.title = _make_title('Residuals', data.name)

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        self.plot_prefs['ylog'] = False
        kwargs.pop('ylog', True)
        Plot.plot(self, self.x, self.y, yerr=self.yerr, xerr=self.xerr,
                  title=self.title, xlabel=self.xlabel, ylabel=self.ylabel,
                  overplot=overplot, clearwindow=clearwindow, **kwargs)


class ResidContour(ModelContour):
    "Derived class for creating 2D residual contours (data-model)"
    contour_prefs = backend.get_resid_contour_defaults()

    def _calc_resid(self, ylist):
        return ylist[0] - ylist[1]

    def prepare(self, data, model, stat):
        (self.x0, self.x1, self.y, self.xlabel,
         self.ylabel) = data.to_contour(yfunc=model)

        self.y = self._calc_resid(self.y)
        self.title = _make_title('Residuals', data.name)

    def contour(self, overcontour=False, clearwindow=True, **kwargs):
        Contour.contour(self, self.x0, self.x1, self.y, levels=self.levels,
                        title=self.title,
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        overcontour=overcontour, clearwindow=clearwindow,
                        **kwargs)


class RatioPlot(ModelPlot):
    """Create plots of the ratio of data to model per point.

    The value of data / model is plotted for each point.

    Attributes
    ----------
    plot_prefs : dict
       The preferences for the plot.
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

    plot_prefs = backend.get_ratio_plot_defaults()

    def _calc_ratio(self, ylist):
        data = numpy.array(ylist[0])
        model = numpy.asarray(ylist[1])
        bad = numpy.where(model == 0.0)
        data[bad] = 0.0
        model[bad] = 1.0
        return (data / model)

    def prepare(self, data, model, stat):
        (self.x, y, self.yerr, self.xerr,
         self.xlabel, self.ylabel) = data.to_plot(model)

        self.y = self._calc_ratio(y)

        # See the discussion in DataPlot.prepare
        try:
            yerrorbars = self.plot_prefs['yerrorbars']
        except KeyError:
            yerrorbars = True

        if stat.name in _stats_noerr:
            self.yerr = data.get_yerr(True, Chi2XspecVar.calc_staterror)
            self.yerr = self.yerr / y[1]
            if yerrorbars:
                warning(_errorbar_warning(stat))
        else:
            staterr = data.get_yerr(True, stat.calc_staterror)
            self.yerr = staterr / y[1]

        self.ylabel = 'Data / Model'
        self.title = _make_title('Ratio of Data to Model', data.name)

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        self.plot_prefs['ylog'] = False
        kwargs.pop('ylog', True)
        Plot.plot(self, self.x, self.y, yerr=self.yerr, xerr=self.xerr,
                  title=self.title, xlabel=self.xlabel, ylabel=self.ylabel,
                  overplot=overplot, clearwindow=clearwindow, **kwargs)


class RatioContour(ModelContour):
    "Derived class for creating 2D ratio contours (data divided by model)"
    contour_prefs = backend.get_ratio_contour_defaults()

    def _calc_ratio(self, ylist):
        data = numpy.array(ylist[0])
        model = numpy.asarray(ylist[1])
        bad = numpy.where(model == 0.0)
        data[bad] = 0.0
        model[bad] = 1.0
        return (data / model)

    def prepare(self, data, model, stat):
        (self.x0, self.x1, self.y, self.xlabel,
         self.ylabel) = data.to_contour(yfunc=model)

        self.y = self._calc_ratio(self.y)
        self.title = _make_title('Ratio of Data to Model', data.name)

    def contour(self, overcontour=False, clearwindow=True, **kwargs):
        Contour.contour(self, self.x0, self.x1, self.y, levels=self.levels,
                        title=self.title, xlabel=self.xlabel,
                        ylabel=self.ylabel,
                        overcontour=overcontour, clearwindow=clearwindow,
                        **kwargs)


class Confidence1D(DataPlot):

    plot_prefs = backend.get_confid_plot_defaults()

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
        DataPlot.__init__(self)

    def __setstate__(self, state):
        self.__dict__.update(state)

        if 'stat' not in state:
            self.__dict__['stat'] = None

        if 'parval' not in state:
            self.__dict__['parval'] = None

        if 'numcores' not in state:
            self.__dict__['numcores'] = None

    def __str__(self):
        x = self.x
        if self.x is not None:
            x = numpy.array2string(self.x, separator=',', precision=4,
                                   suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4,
                                   suppress_small=False)

        return (('x     = %s\n' +
                 'y     = %s\n' +
                 'min   = %s\n' +
                 'max   = %s\n' +
                 'nloop = %s\n' +
                 'delv  = %s\n' +
                 'fac   = %s\n' +
                 'log   = %s') %
                (x,
                 y,
                 self.min,
                 self.max,
                 self.nloop,
                 self.delv,
                 self.fac,
                 self.log))

    def _repr_html_(self):
        """Return a HTML (string) representation of the confidence 1D plot."""
        return backend.as_html_contour1d(self)

    def prepare(self, min=None, max=None, nloop=20,
                delv=None, fac=1, log=False, numcores=None):
        """Set the data to plot.

        This defines the range over which the statistic will
        be calculated, but does not perform the evaluation.

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
                min = r.parmins[index]
                if min is not None and not numpy.isnan(min):
                    self.min = par.val + min

            if self.max is None:
                self.max = par.max
                max = r.parmaxes[index]
                if max is not None and not numpy.isnan(max):
                    self.max = par.val + max

            v = (self.max + self.min) / 2.
            dv = numpy.fabs(v - self.min)
            self.min = v - self.fac * dv
            self.max = v + self.fac * dv

        if not numpy.isscalar(self.min) or not numpy.isscalar(self.max):
            raise ConfidenceErr('badarg', 'Parameter limits', 'scalars')

        # check user limits for errors
        if self.min >= self.max:
            raise ConfidenceErr('badlimits')

        if self.nloop <= 1:
            raise ConfidenceErr('badarg', 'Nloop parameter', '> 1')

        if self.min < par.min:
            self.min = par.min
        if self.max > par.max:
            self.max = par.max

        if self.delv is None:
            self.x = numpy.linspace(self.min, self.max, self.nloop)
        else:
            eps = numpy.finfo(numpy.float32).eps
            self.x = numpy.arange(self.min, self.max + self.delv - eps,
                                  self.delv)

        x = self.x
        if self.log:
            if self.max <= 0.0 or self.min <= 0.0:
                raise ConfidenceErr('badarg', 'Log scale',
                                    'on positive boundaries')
            self.max = numpy.log10(self.max)
            self.min = numpy.log10(self.min)

            x = numpy.linspace(self.min, self.max, len(x))

        return x

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
        """

        if type(fit.stat) in (LeastSq,):
            raise ConfidenceErr('badargconf', fit.stat.name)

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        if self.log:
            self.plot_prefs['xlog'] = True

        Plot.plot(self, self.x, self.y, title=self.title, xlabel=self.xlabel,
                  ylabel=self.ylabel, overplot=overplot,
                  clearwindow=clearwindow, **kwargs)

        # Note: the user arguments are not applied to the lines
        #
        if self.stat is not None:
            Plot.hline(self.stat, linecolor="green", linestyle="dash",
                       linewidth=1.5, overplot=True, clearwindow=False)

        if self.parval is not None:
            Plot.vline(self.parval, linecolor="orange", linestyle="dash",
                       linewidth=1.5, overplot=True, clearwindow=False)

        if self.log:
            self.plot_prefs['xlog'] = False


class Confidence2D(DataContour, Point):

    contour_prefs = backend.get_confid_contour_defaults()
    point_prefs = backend.get_confid_point_defaults()

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
        DataContour.__init__(self)

    def __setstate__(self, state):
        self.__dict__.update(state)

        if 'stat' not in state:
            self.__dict__['stat'] = None

        if 'numcores' not in state:
            self.__dict__['numcores'] = None

    def __str__(self):
        x0 = self.x0
        if self.x0 is not None:
            x0 = numpy.array2string(self.x0, separator=',', precision=4,
                                    suppress_small=False)

        x1 = self.x1
        if self.x1 is not None:
            x1 = numpy.array2string(self.x1, separator=',', precision=4,
                                    suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4,
                                   suppress_small=False)

        return (('x0      = %s\n' +
                 'x1      = %s\n' +
                 'y       = %s\n' +
                 'min     = %s\n' +
                 'max     = %s\n' +
                 'nloop   = %s\n' +
                 'fac     = %s\n' +
                 'delv    = %s\n' +
                 'log     = %s\n' +
                 'sigma   = %s\n' +
                 'parval0 = %s\n' +
                 'parval1 = %s\n' +
                 'levels  = %s') %
                (x0,
                 x1,
                 y,
                 self.min,
                 self.max,
                 self.nloop,
                 self.fac,
                 self.delv,
                 self.log,
                 self.sigma,
                 self.parval0,
                 self.parval1,
                 self.levels))

    def _repr_html_(self):
        """Return a HTML (string) representation of the confidence 2D plot."""
        return backend.as_html_contour2d(self)

    def prepare(self, min=None, max=None, nloop=(10, 10),
                delv=None, fac=4, log=(False, False),
                sigma=(1, 2, 3), levels=None, numcores=None):
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

        # Issue #1093 points out that if min or max is a tuple
        # we can have a problem, as below the code can assign
        # to an element of one of them. So ensure we have a list
        # not a tuple.
        #
        if self.min is not None:
            try:
                self.min = list(self.min)
            except TypeError:
                raise ConfidenceErr('badarg', 'Parameter limits', 'a list') from None

        if self.max is not None:
            try:
                self.max = list(self.max)
            except TypeError:
                raise ConfidenceErr('badarg', 'Parameter limits', 'a list') from None

        self.stat = fit.calc_stat()
        self.xlabel = par0.fullname
        self.ylabel = par1.fullname
        self.parval0 = par0.val
        self.parval1 = par1.val

        if self.levels is None:
            stat = self.stat
            if self.sigma is None or numpy.isscalar(self.sigma):
                raise ConfidenceErr('needlist', 'sigma bounds')
            thelevels = numpy.zeros(len(self.sigma), SherpaFloat)
            for i in range(len(self.sigma)):
                thelevels[i] = stat - (2. * numpy.log(1. - erf(
                    self.sigma[i] / numpy.sqrt(2.))))
            self.levels = thelevels

        if self.min is None or self.max is None:
            oldestmethod = fit.estmethod
            fit.estmethod = Covariance()

            r = fit.est_errors()
            fit.estmethod = oldestmethod

            index0 = list(r.parnames).index(par0.fullname)
            index1 = list(r.parnames).index(par1.fullname)

            if self.min is None:
                self.min = numpy.array([par0.min, par1.min])
                min0 = r.parmins[index0]
                min1 = r.parmins[index1]

                if min0 is not None and not numpy.isnan(min0):
                    self.min[0] = par0.val + min0

                if min1 is not None and not numpy.isnan(min1):
                    self.min[1] = par1.val + min1

            if self.max is None:
                self.max = numpy.array([par0.max, par1.max])
                max0 = r.parmaxes[index0]
                max1 = r.parmaxes[index1]

                if max0 is not None and not numpy.isnan(max0):
                    self.max[0] = par0.val + max0

                if max1 is not None and not numpy.isnan(max1):
                    self.max[1] = par1.val + max1

            for i in [0, 1]:
                v = (self.max[i] + self.min[i]) / 2.
                dv = numpy.fabs(v - self.min[i])
                self.min[i] = v - self.fac * dv
                self.max[i] = v + self.fac * dv

        hmin = numpy.array([par0.min, par1.min])
        hmax = numpy.array([par0.max, par1.max])

        for i in [0, 1]:
            # check user limits for errors
            if numpy.isscalar(self.min) or numpy.isscalar(self.max):
                raise ConfidenceErr('badarg', 'Parameter limits', 'a list')

            if self.min[i] >= self.max[i]:
                raise ConfidenceErr('badlimits')

            if numpy.isscalar(self.nloop) or self.nloop[i] <= 1:
                raise ConfidenceErr('badarg', 'Nloop parameter',
                                    'a list with elements > 1')

            if self.min[i] < hmin[i]:
                self.min[i] = hmin[i]
            if self.max[i] > hmax[i]:
                self.max[i] = hmax[i]

        if self.delv is None:
            self.x0 = numpy.linspace(self.min[0], self.max[0], self.nloop[0])
            self.x1 = numpy.linspace(self.min[1], self.max[1], self.nloop[1])

        else:
            eps = numpy.finfo(numpy.float32).eps
            self.x0 = numpy.arange(self.min[0],
                                   self.max[0] + self.delv[0] - eps,
                                   self.delv[0])
            self.x1 = numpy.arange(self.min[1],
                                   self.max[1] + self.delv[1] - eps,
                                   self.delv[1])

        # x = numpy.array([self.x0, self.x1])
        x = [self.x0, self.x1]

        self.x0, self.x1 = numpy.meshgrid(self.x0, self.x1)
        self.x0 = self.x0.ravel()
        self.x1 = self.x1.ravel()

        for i in [0, 1]:
            if self.log[i]:
                if self.max[i] <= 0.0 or self.min[i] <= 0.0:
                    raise ConfidenceErr('badarg', 'Log scale',
                                        'on positive boundaries')
                self.max[i] = numpy.log10(self.max[i])
                self.min[i] = numpy.log10(self.min[i])
                x[i] = numpy.linspace(self.min[i], self.max[i], len(x[i]))

        x0, x1 = numpy.meshgrid(x[0], x[1])
        return numpy.array([x0.ravel(), x1.ravel()]).T

    def calc(self, fit, par0, par1):
        if type(fit.stat) in (LeastSq,):
            raise ConfidenceErr('badargconf', fit.stat.name)

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

        Contour.contour(self, self.x0, self.x1, self.y, levels=self.levels,
                        title=self.title, xlabel=self.xlabel,
                        ylabel=self.ylabel, overcontour=overcontour,
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


class IntervalProjectionWorker():
    def __init__(self, log, par, thawed, fit):
        self.log = log
        self.par = par
        self.thawed = thawed
        self.fit = fit

    def __call__(self, val):
        if self.log:
            val = numpy.power(10, val)
        self.par.val = val
        if len(self.thawed) > 1:
            r = self.fit.fit()
            return r.statval
        return self.fit.calc_stat()


class IntervalProjection(Confidence1D):

    def __init__(self):
        self.fast = True
        Confidence1D.__init__(self)

    def prepare(self, fast=True, min=None, max=None, nloop=20,
                delv=None, fac=1, log=False, numcores=None):
        self.fast = fast
        Confidence1D.prepare(self, min, max, nloop, delv, fac, log, numcores)

    def calc(self, fit, par, methoddict=None, cache=True):
        self.title = 'Interval-Projection'

        Confidence1D.calc(self, fit, par)

        if par.frozen:
            raise ConfidenceErr('frozen', par.fullname, 'interval projection')

        thawed = [i for i in fit.model.pars if not i.frozen]

        if par not in thawed:
            raise ConfidenceErr('thawed', par.fullname, fit.model.name)

        # If "fast" option enabled, set fitting method to
        # lmdif if stat is chi-squared,
        # else set to neldermead.

        # If current method is not LM or NM, warn it is not a good
        # method for estimating parameter limits.
        if type(fit.method) not in (NelderMead, LevMar):
            warning(fit.method.name + " is inappropriate for confidence " +
                    "limit estimation")

        oldfitmethod = fit.method
        if (bool_cast(self.fast) is True and methoddict is not None):
            if (isinstance(fit.stat, Likelihood)):
                if (type(fit.method) is not NelderMead):
                    fit.method = methoddict['neldermead']
                    warning("Setting optimization to " + fit.method.name +
                            " for interval projection plot")
            else:
                if (type(fit.method) is not LevMar):
                    fit.method = methoddict['levmar']
                    warning("Setting optimization to " + fit.method.name +
                            " for interval projection plot")

        xvals = self._interval_init(fit, par)
        oldpars = fit.model.thawedpars
        par.freeze()

        try:
            fit.model.startup(cache)

            # store the class methods for startup and teardown
            # these calls are unnecessary for every fit
            startup = fit.model.startup
            fit.model.startup = return_none
            teardown = fit.model.teardown
            fit.model.teardown = return_none

            self.y = numpy.asarray(parallel_map(IntervalProjectionWorker(self.log, par, thawed, fit),
                                                xvals,
                                                self.numcores)
                                   )

        finally:
            # Set back data that we changed
            par.thaw()

            fit.model.startup = startup
            fit.model.teardown = teardown

            fit.model.teardown()
            fit.model.thawedpars = oldpars
            fit.method = oldfitmethod


class IntervalUncertaintyWorker():
    def __init__(self, log, par, fit):
        self.log = log
        self.par = par
        self.fit = fit

    def __call__(self, val):
        if self.log:
            val = numpy.power(10, val)
        self.par.val = val
        return self.fit.calc_stat()


class IntervalUncertainty(Confidence1D):

    def calc(self, fit, par, methoddict=None, cache=True):
        self.title = 'Interval-Uncertainty'

        Confidence1D.calc(self, fit, par)
        if par.frozen:
            raise ConfidenceErr('frozen', par.fullname, 'interval uncertainty')

        thawed = [i for i in fit.model.pars if not i.frozen]

        if par not in thawed:
            raise ConfidenceErr('thawed', par.fullname, fit.model.name)

        oldpars = fit.model.thawedpars

        xvals = self._interval_init(fit, par)

        for i in thawed:
            i.freeze()

        try:
            fit.model.startup(cache)
            self.y = numpy.asarray(parallel_map(IntervalUncertaintyWorker(self.log, par, fit),
                                                xvals,
                                                self.numcores)
                                   )

        finally:
            # Set back data that we changed
            for i in thawed:
                i.thaw()
            fit.model.teardown()
            fit.model.thawedpars = oldpars


class RegionProjectionWorker():
    def __init__(self, log, par0, par1, thawed, fit):
        self.log = log
        self.par0 = par0
        self.par1 = par1
        self.thawed = thawed
        self.fit = fit

    def __call__(self, pars):
        for ii in [0, 1]:
            if self.log[ii]:
                pars[ii] = numpy.power(10, pars[ii])
        (self.par0.val, self.par1.val) = pars
        if len(self.thawed) > 2:
            r = self.fit.fit()
            return r.statval
        return self.fit.calc_stat()


def return_none(cache=None):
    """
    dummy implementation of callback for multiprocessing
    """
    return None


class RegionProjection(Confidence2D):

    def __init__(self):
        self.fast = True
        Confidence2D.__init__(self)

    def prepare(self, fast=True, min=None, max=None, nloop=(10, 10),
                delv=None, fac=4, log=(False, False),
                sigma=(1, 2, 3), levels=None, numcores=None):
        self.fast = fast
        Confidence2D.prepare(self, min, max, nloop, delv, fac, log,
                             sigma, levels, numcores)

    def calc(self, fit, par0, par1, methoddict=None, cache=True):
        self.title = 'Region-Projection'

        Confidence2D.calc(self, fit, par0, par1)
        if par0.frozen:
            raise ConfidenceErr('frozen', par0.fullname, 'region projection')
        if par1.frozen:
            raise ConfidenceErr('frozen', par1.fullname, 'region projection')

        thawed = [i for i in fit.model.pars if not i.frozen]

        if par0 not in thawed:
            raise ConfidenceErr('thawed', par0.fullname, fit.model.name)
        if par1 not in thawed:
            raise ConfidenceErr('thawed', par1.fullname, fit.model.name)

        # If "fast" option enabled, set fitting method to
        # lmdif if stat is chi-squared,
        # else set to neldermead

        # If current method is not LM or NM, warn it is not a good
        # method for estimating parameter limits.
        if type(fit.method) not in (NelderMead, LevMar):
            warning(fit.method.name + " is inappropriate for confidence " +
                    "limit estimation")

        oldfitmethod = fit.method
        if (bool_cast(self.fast) is True and methoddict is not None):
            if (isinstance(fit.stat, Likelihood)):
                if (type(fit.method) is not NelderMead):
                    fit.method = methoddict['neldermead']
                    warning("Setting optimization to " + fit.method.name +
                            " for region projection plot")
            else:
                if (type(fit.method) is not LevMar):
                    fit.method = methoddict['levmar']
                    warning("Setting optimization to " + fit.method.name +
                            " for region projection plot")

        oldpars = fit.model.thawedpars

        try:
            fit.model.startup(cache)

            # store the class methods for startup and teardown
            # these calls are unnecessary for every fit
            startup = fit.model.startup
            fit.model.startup = return_none
            teardown = fit.model.teardown
            fit.model.teardown = return_none

            grid = self._region_init(fit, par0, par1)

            par0.freeze()
            par1.freeze()

            self.y = numpy.asarray(parallel_map(RegionProjectionWorker(self.log, par0, par1, thawed, fit),
                                                grid,
                                                self.numcores)
                                   )

        finally:
            # Set back data after we changed it
            par0.thaw()
            par1.thaw()

            fit.model.startup = startup
            fit.model.teardown = teardown

            fit.model.teardown()
            fit.model.thawedpars = oldpars
            fit.method = oldfitmethod


class RegionUncertaintyWorker():
    def __init__(self, log, par0, par1, fit):
        self.log = log
        self.par0 = par0
        self.par1 = par1
        self.fit = fit

    def __call__(self, pars):
        for ii in [0, 1]:
            if self.log[ii]:
                pars[ii] = numpy.power(10, pars[ii])
        (self.par0.val, self.par1.val) = pars
        return self.fit.calc_stat()


class RegionUncertainty(Confidence2D):

    def calc(self, fit, par0, par1, methoddict=None, cache=True):
        self.title = 'Region-Uncertainty'

        Confidence2D.calc(self, fit, par0, par1)
        if par0.frozen:
            raise ConfidenceErr('frozen', par0.fullname, 'region uncertainty')
        if par1.frozen:
            raise ConfidenceErr('frozen', par1.fullname, 'region uncertainty')

        thawed = [i for i in fit.model.pars if not i.frozen]

        if par0 not in thawed:
            raise ConfidenceErr('thawed', par0.fullname, fit.model.name)
        if par1 not in thawed:
            raise ConfidenceErr('thawed', par1.fullname, fit.model.name)

        oldpars = fit.model.thawedpars

        try:
            fit.model.startup(cache)

            grid = self._region_init(fit, par0, par1)

            for i in thawed:
                i.freeze()

            self.y = numpy.asarray(parallel_map(RegionUncertaintyWorker(self.log, par0, par1, fit),
                                                grid,
                                                self.numcores)
                                   )

        finally:
            # Set back data after we changed it
            for i in thawed:
                i.thaw()
            fit.model.teardown()
            fit.model.thawedpars = oldpars
