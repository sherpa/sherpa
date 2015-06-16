# 
#  Copyright (C) 2009, 2015  Smithsonian Astrophysical Observatory
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

"""
A visualization interface to Sherpa
"""


from __future__ import division
import os
import logging
warning = logging.getLogger(__name__).warning

import numpy
_ = numpy.seterr(invalid='ignore')

from sherpa.utils import NoNewAttributesAfterInit, erf, SherpaFloat, \
    bool_cast, parallel_map, dataspace1d, histogram1d, get_error_estimates
from sherpa.utils.err import PlotErr, StatErr, ConfidenceErr
from sherpa.estmethods import Covariance
from sherpa.optmethods import LevMar, NelderMead
from sherpa.stats import Likelihood, LeastSq, Chi2XspecVar
from sherpa import get_config
from ConfigParser import ConfigParser

config = ConfigParser()
config.read(get_config())

plot_opt = config.get('options','plot_pkg')
plot_opt = str(plot_opt).strip().lower() + '_backend'
if plot_opt == 'matplotlib_backend':
    plot_opt = 'pylab_backend'

if plot_opt == 'none_backend':
    plot_opt = 'dummy_backend'

try:
    backend = __import__(plot_opt, globals(), locals(), [])
except:
    # if the user inputs a malformed backend or it is not found,
    # give a useful warning and fall back on dummy_backend of noops
    warning('failed to import sherpa.plot.%s;' % plot_opt +
            ' plotting routines will not be available')
    import dummy_backend as backend

backend.init()



__all__ = ('Plot', 'Contour', 'Point', 'SplitPlot', 'JointPlot',
           'DataPlot', 'DataContour', 'DelchiPlot', 'ComponentModelPlot',
           'ModelPlot', 'ModelContour', 'FitPlot', 'FitContour',
           'ResidPlot', 'ResidContour','RatioPlot', 'RatioContour',
           'IntervalProjection', 'IntervalUncertainty', 'ChisqrPlot',
           'RegionProjection', 'RegionUncertainty', 'ComponentSourcePlot',
           'PSFPlot','PSFContour','begin', 'end', 'exceptions', 'backend',
           'SourcePlot', 'SourceContour', 'Histogram')

_stats_noerr = ('cash', 'cstat', 'leastsq')

begin = backend.begin
end = backend.end
exceptions = backend.exceptions

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


    def plot(self, x, y, yerr=None, xerr=None, title=None, xlabel=None,
             ylabel=None, overplot=False, clearwindow=True):
        backend.plot(x, y, yerr, xerr, title, xlabel, ylabel, overplot,
                     clearwindow, **self.plot_prefs)

    def overplot(self, *args, **kwargs):
        "Add the data to an existing plot."
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
        
    def contour(self, x0, x1, y, levels=None, title=None, xlabel=None,
                ylabel=None, overcontour=False, clearwindow=True):
        backend.contour(x0, x1, y, levels, title, xlabel, ylabel, overcontour,
                        clearwindow, **self.contour_prefs)

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
        
    def point(self, x, y, overplot=True, clearwindow=False):
        backend.point(x, y, overplot, clearwindow, **self.point_prefs)


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

    def plot(self, xlo, xhi, y, yerr=None, title=None, xlabel=None,
             ylabel=None, overplot=False, clearwindow=True):
        backend.histo(xlo, xhi, y, yerr, title, xlabel, ylabel, overplot,
                     clearwindow, **self.histo_prefs)

    def overplot(self, *args, **kwargs):
        kwargs['overplot'] = True
        self.plot(*args, **kwargs)


class HistogramPlot(Histogram):

    def __init__(self):
        self.xlo = None
        self.xhi = None
        self.y  = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        Histogram.__init__(self)

    def __str__(self):
        xlo = self.xlo
        if self.xlo is not None:
            xlo = numpy.array2string(numpy.asarray(self.xlo), separator=',', precision=4, suppress_small=False)

        xhi = self.xhi
        if self.xhi is not None:
            xhi = numpy.array2string(numpy.asarray(self.xhi), separator=',', precision=4, suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(numpy.asarray(self.y), separator=',', precision=4, suppress_small=False)
        
        return (('xlo    = %s\n' +
                 'xhi    = %s\n' +
                 'y      = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'histo_prefs = %s') %
                ( xlo,
                  xhi,
                  y,
                  self.xlabel,
                  self.ylabel,
                  self.title,
                  self.histo_prefs))

    def plot(self, overplot=False, clearwindow=True):
        Histogram.plot(self, self.xlo, self.xhi, self.y, title=self.title,
                       xlabel=self.xlabel, ylabel=self.ylabel,
                       overplot=overplot, clearwindow=clearwindow)


class PDFPlot(HistogramPlot):

    def __init__(self):
        self.points = None
        HistogramPlot.__init__(self)


    def __str__(self):
        points = self.points
        if self.points is not None:
            points = numpy.array2string(numpy.asarray(self.points), separator=',', precision=4, suppress_small=False)

        return ('points = %s\n'%(points) +
                HistogramPlot.__str__(self))

    def prepare(self, points, bins=12, normed=True, xlabel="x", name="x"):
        self.points = points
        self.y, xx = numpy.histogram(points, bins=bins, normed=normed)
        self.xlo = xx[:-1]
        self.xhi = xx[1:]
        self.ylabel = "probability density"
        self.xlabel = xlabel
        self.title  = "PDF: %s"%name



class CDFPlot(Plot):
    
    median_defaults = dict(linestyle='dash', linecolor='orange',
                           linewidth=1.5)
    lower_defaults = dict(linestyle='dash', linecolor='blue',
                          linewidth=1.5)
    upper_defaults = dict(linestyle='dash', linecolor='blue',
                          linewidth=1.5)

    plot_prefs = backend.get_cdf_plot_defaults()
    
    def __init__(self):
        self.x  = None
        self.y  = None
        self.points = None
        self.median = None
        self.lower  = None
        self.upper  = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        Plot.__init__(self)


    def __str__(self):
        x = self.x
        if self.x is not None:
            x = numpy.array2string(self.x, separator=',', precision=4, suppress_small=False)
        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)

        points = self.points
        if self.points is not None:
            points = numpy.array2string(self.points, separator=',', precision=4, suppress_small=False)

        return (('points = %s\n' +
                 'x      = %s\n' +
                 'y      = %s\n' +
                 'median = %s\n' +
                 'lower  = %s\n' +
                 'upper  = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'plot_prefs = %s') %
                ( x,
                  y,
                  points,
                  self.median,
                  self.lower,
                  self.upper,
                  self.xlabel,
                  self.ylabel,
                  self.title,
                  self.plot_prefs))


    def prepare(self, points, xlabel="x", name="x"):
        self.points = points
        self.x = numpy.sort(points)
        (self.median, self.lower,
         self.upper) = get_error_estimates(self.x, True)
        xsize = len(self.x)
        self.y = (numpy.arange(xsize) + 1.0) / xsize
        self.xlabel = xlabel
        self.ylabel = "p(<=%s)"%(xlabel)
        self.title  = "CDF: %s"%(name)


    def plot(self, overplot=False, clearwindow=True):
        Plot.plot(self, self.x, self.y, title=self.title,
                  xlabel=self.xlabel, ylabel=self.ylabel,
                  overplot=overplot, clearwindow=clearwindow)
        Plot.vline(self.median, overplot=True, clearwindow=False,
                   **self.median_defaults)
        Plot.vline(self.lower, overplot=True, clearwindow=False, 
                   **self.lower_defaults)
        Plot.vline(self.upper, overplot=True, clearwindow=False,
                   **self.upper_defaults)


class LRHistogram(HistogramPlot):
    "Derived class for creating 1D likelihood ratio distribution plots"
    def __init__(self):
        self.ratios=None
        self.lr=None
        self.ppp=None
        HistogramPlot.__init__(self)


    def __str__(self):
        ratios = self.ratios
        if self.ratios is not None:
            ratios = numpy.array2string(numpy.asarray(self.ratios), separator=',', precision=4, suppress_small=False)

        return '\n'.join(['ratios = %s' % ratios,
                          'lr = %s' % str(self.lr),
                          HistogramPlot.__str__(self)])


    def prepare(self, ratios, bins, niter, lr, ppp):
        self.ppp = float(ppp)
        self.lr = float(lr)
        y = numpy.asarray(ratios)
        self.ratios = y
        self.xlo, self.xhi = dataspace1d(y.min(), y.max(), numbins=bins+1)[:2]
        y = histogram1d(y, self.xlo, self.xhi)
	self.y = y/float(niter)
        self.title = "Likelihood Ratio Distribution"
        self.xlabel = "Likelihood Ratio"
        self.ylabel = "Frequency"


    def plot(self, overplot=False, clearwindow=True):
        Histogram.plot(self, self.xlo, self.xhi, self.y, title=self.title,
                       xlabel=self.xlabel, ylabel=self.ylabel,
                       overplot=overplot, clearwindow=clearwindow)

        if self.lr is not None:
            if self.lr <= self.xhi.max() and self.lr >= self.xlo.min():
                Plot.vline(self.lr, linecolor="orange", linestyle="solid",
                           linewidth=1.5, overplot=True, clearwindow=False)


class SplitPlot(Plot,Contour):
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
                ( self.rows,
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
        row, col = numpy.where(self._used == False)
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

class JointPlot(SplitPlot):
    def __init__(self):
        SplitPlot.__init__(self)

    def _clear_window(self, row, col, erase=True):
        if not self._cleared_window:
            if erase:
                backend.clear_window()
            # FIXME: misuse of kwarg clearaxes for chips backend
            backend.set_jointplot(row, col, self.rows, self.cols, False)
            self._cleared_window = True
        else:
            # FIXME: misuse of kwarg clearaxes for chips backend
            backend.set_jointplot(row, col, self.rows, self.cols)
    
    def plottop(self, plot, *args, **kwargs):
        clearaxes = kwargs.get('clearwindow', True)
        self._clear_window(0,0,clearaxes)
        
        # FIXME: should not know about FitPlot, terrible hack to remove label
        if isinstance(plot,FitPlot):
            plot.dataplot.xlabel=''
            plot.modelplot.xlabel=''
        else:
            plot.xlabel=''

        kwargs['clearwindow'] = False
        plot.plot(*args, **kwargs)

    def plotbot(self, plot, *args, **kwargs):
        self._clear_window(1,0)
          
        # FIXME: terrible hack to remove title from bottom
        plot.title=''
        kwargs['clearwindow'] = False
        plot.plot(*args, **kwargs)


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

    """

    plot_prefs = backend.get_data_plot_defaults()
    
    def __init__(self):
        self.x  = None
        self.y  = None
        self.yerr = None
        self.xerr = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        Plot.__init__(self)

    def __str__(self):
        x = self.x
        if self.x is not None:
            x = numpy.array2string(self.x, separator=',', precision=4, suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)

        yerr = self.yerr
        if self.yerr is not None:
            yerr = numpy.array2string(self.yerr, separator=',', precision=4, suppress_small=False)

        xerr = self.xerr
        if self.xerr is not None:
            xerr = numpy.array2string(self.xerr, separator=',', precision=4, suppress_small=False)
            
        return (('x      = %s\n' +
                 'y      = %s\n' +
                 'yerr   = %s\n' +
                 'xerr   = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'plot_prefs = %s') %
                ( x,
                  y,
                  yerr,
                  xerr,
                  self.xlabel,
                  self.ylabel,
                  self.title,
                  self.plot_prefs))
        
    def prepare(self, data, stat=None):
        (self.x, self.y, self.yerr, self.xerr, self.xlabel,
         self.ylabel) = data.to_plot()

        #if self.yerr is None and stat is not None:
        if stat is not None:
            msg = ("The displayed errorbars have been supplied with the data or calculated using chi2xspecvar; the errors are not used in fits with %s" % stat.name)
            if stat.name in _stats_noerr:
                self.yerr = data.get_yerr(True, Chi2XspecVar.calc_staterror)
                warning(msg)
            else:
                try:
                    self.yerr = data.get_yerr(True, stat.calc_staterror)
                except:
                    warning(msg + "\nzeros or negative values found")

        self.title = data.name

    def plot(self, overplot=False, clearwindow=True):
        Plot.plot(self, self.x, self.y, self.yerr, self.xerr, self.title,
                  self.xlabel, self.ylabel, overplot, clearwindow)


class TracePlot(DataPlot):

    plot_prefs = backend.get_model_plot_defaults()
    
    def prepare(self, points, xlabel="x", name="x"):
        self.x = numpy.arange(len(points), dtype=SherpaFloat)
        self.y = points
        self.xlabel = "iteration"
        self.ylabel = name
        self.title  = "Trace: %s"%(name)


class ScatterPlot(DataPlot):

    plot_prefs = backend.get_scatter_plot_defaults()

    def prepare(self, x, y, xlabel="x", ylabel="y", name="(x,y)"):
        self.x = numpy.asarray(x, dtype=SherpaFloat)
        self.y = numpy.asarray(y, dtype=SherpaFloat)
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title  = "Scatter: %s"%(name)


class PSFKernelPlot(DataPlot):
    "Derived class for creating 1D PSF kernel data plots"

    def prepare(self, psf, data=None, stat=None):
        psfdata = psf.get_kernel(data)
        DataPlot.prepare(self, psfdata, stat)
        #self.ylabel = 'PSF value'
        #self.xlabel = 'PSF Kernel size'
        self.title  = 'PSF Kernel'


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
        self.y  = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        self.levels = None
        Contour.__init__(self)

    def __str__(self):
        x0 = self.x0
        if self.x0 is not None:
            x0 = numpy.array2string(self.x0, separator=',', precision=4, suppress_small=False)

        x1 = self.x1
        if self.x1 is not None:
            x1 = numpy.array2string(self.x1, separator=',', precision=4, suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)

        return (('x0     = %s\n' +
                 'x1     = %s\n' +
                 'y      = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'levels = %s\n' +
                 'contour_prefs = %s') %
                ( x0,
                  x1,
                  y,
                  self.xlabel,
                  self.ylabel,
                  self.title,
                  self.levels,
                  self.contour_prefs))
        
    def prepare(self, data, stat=None):
        (self.x0, self.x1, self.y, self.xlabel,
         self.ylabel) = data.to_contour()
        self.title = data.name
           
    def contour(self, overcontour=False, clearwindow=True):
        Contour.contour(self, self.x0, self.x1, self.y,
                        self.levels,self.title, self.xlabel,
                        self.ylabel, overcontour, clearwindow)

class PSFKernelContour(DataContour):
    "Derived class for creating 2D PSF Kernel contours"

    def prepare(self, psf, data=None, stat=None):
        psfdata = psf.get_kernel(data)
        DataContour.prepare(self, psfdata)
        #self.xlabel = 'PSF Kernel size x0'
        #self.ylabel = 'PSF Kernel size x1'
        self.title  = 'PSF Kernel'


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

    """

    plot_prefs = backend.get_model_plot_defaults()
    
    def __init__(self):
        self.x  = None
        self.y  = None
        self.yerr = None
        self.xerr = None
        self.xlabel = None
        self.ylabel = None
        self.title = 'Model'
        Plot.__init__(self)

    def __str__(self):
        x = self.x
        if self.x is not None:
            x = numpy.array2string(self.x, separator=',', precision=4, suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)

        yerr = self.yerr
        if self.yerr is not None:
            yerr = numpy.array2string(self.yerr, separator=',', precision=4, suppress_small=False)

        xerr = self.xerr
        if self.xerr is not None:
            xerr = numpy.array2string(self.xerr, separator=',', precision=4, suppress_small=False)
        
        return (('x      = %s\n' +
                 'y      = %s\n' +
                 'yerr   = %s\n' +
                 'xerr   = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'plot_prefs = %s') %
                ( x,
                  y,
                  yerr,
                  xerr,
                  self.xlabel,
                  self.ylabel,
                  self.title,
                  self.plot_prefs))

    def prepare(self, data, model, stat=None):
        (self.x, self.y, self.yerr, self.xerr,
         self.xlabel, self.ylabel) = data.to_plot(yfunc=model)
        self.y = self.y[1]
            
    def plot(self, overplot=False, clearwindow=True):
        Plot.plot(self, self.x, self.y, title=self.title, xlabel=self.xlabel,
                  ylabel=self.ylabel, overplot=overplot, clearwindow=clearwindow)


class ComponentModelPlot(ModelPlot):

    plot_prefs = backend.get_component_plot_defaults()

    def prepare(self, data, model, stat=None):
        ModelPlot.prepare(self, data, model, stat)
        self.title = 'Model component: %s' % model.name


class ComponentTemplateModelPlot(ComponentModelPlot):

    def prepare(self, data, model, stat=None):
        self.x = model.get_x()
        self.y = model.get_y()
        self.xlabel = data.get_xlabel()
        self.ylabel = data.get_ylabel()
        self.title = 'Model component: %s' % model.name

class SourcePlot(ModelPlot):
    """Create 1D plots of unconcolved model values.

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
        SourcePlot.prepare(self, data, model, stat)
        self.title = 'Source model component: %s' % model.name


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
        self.title = 'Source model component: %s' % model.name




class PSFPlot(DataPlot):
    "Derived class for creating 1D PSF kernel data plots"

    def prepare(self, psf, data=None, stat=None):
        psfdata = psf.get_kernel(data, False)
        DataPlot.prepare(self, psfdata, stat)
        self.title = psf.kernel.name


class ModelContour(Contour):
    "Derived class for creating 2D model contours"
    contour_prefs = backend.get_model_contour_defaults()
    
    def __init__(self):
        self.x0 = None
        self.x1 = None
        self.y  = None
        self.xlabel = None
        self.ylabel = None
        self.title = 'Model'
        self.levels = None
        Contour.__init__(self)

    def __str__(self):
        x0 = self.x0
        if self.x0 is not None:
            x0 = numpy.array2string(self.x0, separator=',', precision=4, suppress_small=False)
            
        x1 = self.x1
        if self.x1 is not None:
            x1 = numpy.array2string(self.x1, separator=',', precision=4, suppress_small=False)
            
        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)
            
        return (('x0     = %s\n' +
                 'x1     = %s\n' +
                 'y      = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'levels = %s\n' +
                 'contour_prefs = %s') %
                ( x0,
                  x1,
                  y,
                  self.xlabel,
                  self.ylabel,
                  self.title,
                  self.levels,
                  self.contour_prefs))

    def prepare(self, data, model, stat):
        (self.x0, self.x1, self.y, self.xlabel,
         self.ylabel) = data.to_contour(yfunc=model)
        self.y = self.y[1]

    def contour(self, overcontour=False, clearwindow=True):
        Contour.contour(self, self.x0, self.x1, self.y, levels=self.levels,
                        title=self.title, xlabel=self.xlabel, ylabel=self.ylabel,
                        overcontour=overcontour, clearwindow=clearwindow)


class PSFContour(DataContour):
    "Derived class for creating 2D PSF contours"

    def prepare(self, psf, data=None, stat=None):
        psfdata = psf.get_kernel(data, False)
        DataContour.prepare(self, psfdata)
        self.title  = psf.kernel.name


class SourceContour(ModelContour):
    "Derived class for creating 2D model contours"
    def __init__(self):
        ModelContour.__init__(self)
        self.title = 'Source'


class FitPlot(Plot):
    "Derived class for creating 1D combination data and model plots"
    plot_prefs = backend.get_fit_plot_defaults()
    
    def __init__(self):
        self.dataplot = None
        self.modelplot = None
        Plot.__init__(self)

    def __str__(self):
        data_title=None
        if self.dataplot is not None:
            data_title = self.dataplot.title
        
        model_title=None
        if self.modelplot is not None:
            model_title = self.modelplot.title
        
        return (('dataplot   = %s'+'\n%s\n\n'+'modelplot  = %s'+'\n%s') %
                (data_title,
                 self.dataplot,
                 model_title,
                 self.modelplot))

    def prepare(self, dataplot, modelplot):
        self.dataplot = dataplot
        self.modelplot = modelplot

    def plot(self, overplot=False, clearwindow=True):
        if self.dataplot is None or self.modelplot is None:
            raise PlotErr("nodataormodel")
        self.dataplot.plot(overplot, clearwindow)
        self.modelplot.overplot()


class FitContour(Contour):
    "Derived class for creating 2D combination data and model contours"
    contour_prefs = backend.get_fit_contour_defaults()

    def __init__(self):
        self.datacontour = None
        self.modelcontour = None
        Contour.__init__(self)

    def __str__(self):
        data_title=None
        if self.datacontour is not None:
            data_title = self.datacontour.title
        
        model_title=None
        if self.modelcontour is not None:
            model_title = self.modelcontour.title
        return (('datacontour = %s'+'\n%s\n\n'+'modelcontour = %s'+'\n%s') %
                (data_title,
                 self.datacontour,
                 model_title,
                 self.modelcontour))

    def prepare(self, datacontour, modelcontour):
        self.datacontour = datacontour
        self.modelcontour = modelcontour
        
    def contour(self, overcontour=False, clearwindow=True):
        self.datacontour.contour( overcontour, clearwindow)
        self.modelcontour.overcontour()


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

    """

    plot_prefs = backend.get_resid_plot_defaults()

    def _calc_delchi(self, ylist, staterr):
        return (ylist[0] - ylist[1])/staterr
    
    def prepare(self, data, model, stat):
        (self.x, y, staterr, self.xerr,
         self.xlabel, self.ylabel) = data.to_plot(model)
        
        if staterr is None:
            if stat.name in _stats_noerr:
                raise StatErr('badstat', "DelchiPlot", stat.name)
            staterr = data.get_yerr(True,stat.calc_staterror)

        self.y = self._calc_delchi(y, staterr)
        self.yerr = staterr/staterr
        self.ylabel = 'Sigma'
        self.title = _make_title('Sigma Residuals', data.name)

    def plot(self, overplot=False, clearwindow=True):
        Plot.plot(self, self.x, self.y, self.yerr, self.xerr, self.title,
                  self.xlabel, self.ylabel, overplot, clearwindow)

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
        return (ylist[0] - ylist[1])*(ylist[0] - ylist[1]) / (staterr*staterr)

    def prepare(self, data, model, stat):
        (self.x, y, staterr, self.xerr,
         self.xlabel, self.ylabel) = data.to_plot(model)

        #if staterr is None:
        if stat.name in _stats_noerr:
            raise StatErr('badstat', "ChisqrPlot", stat.name)
        staterr = data.get_yerr(True,stat.calc_staterror)

        self.y = self._calc_chisqr(y, staterr)
        # NOTE: these are not displayed as LaTeX by the matplotlib backend
        self.ylabel = '\chi^2'
        self.title = _make_title('\chi^2', data.name)

    def plot(self, overplot=False, clearwindow=True):
        Plot.plot(self, self.x, self.y, title=self.title, xlabel=self.xlabel,
                  ylabel=self.ylabel, overplot=overplot, clearwindow=clearwindow)


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
       The Y value for each point: data-model.
    xerr : array_like
       The half-width of each X "bin", if set.
    yerr : array_like
       The error on the Y value, if set.
    xlabel, ylabel, title : str
       Plot labels.

    """

    plot_prefs = backend.get_resid_plot_defaults()

    def _calc_resid(self, ylist):
        return ylist[0] - ylist[1]
    
    def prepare(self, data, model, stat):
        (self.x, y, self.yerr, self.xerr,
         self.xlabel, self.ylabel) = data.to_plot(model)

        self.y = self._calc_resid(y)
#        if self.yerr is None:
        if stat.name in _stats_noerr:
            self.yerr = data.get_yerr(True, Chi2XspecVar.calc_staterror)
            warning("The displayed errorbars have been supplied with the data or calculated using chi2xspecvar; the errors are not used in fits with %s" % stat.name)
        else:
            self.yerr = data.get_yerr(True,stat.calc_staterror)

        self.title = _make_title('Residuals', data.name)

    def plot(self, overplot=False, clearwindow=True):
        Plot.plot(self, self.x, self.y, self.yerr, self.xerr, self.title,
                  self.xlabel, self.ylabel, overplot, clearwindow)

    
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

    def contour(self, overcontour=False, clearwindow=True):
        Contour.contour(self, self.x0, self.x1, self.y, levels=self.levels,
                        title=self.title, xlabel=self.xlabel, ylabel=self.ylabel,
                        overcontour=overcontour, clearwindow=clearwindow)    


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
       The Y value for each point: data-model.
    xerr : array_like
       The half-width of each X "bin", if set.
    yerr : array_like
       The error on the Y value, if set.
    xlabel, ylabel, title : str
       Plot labels.

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
        #if self.yerr is None:
        if stat.name in _stats_noerr:
            self.yerr = data.get_yerr(True, Chi2XspecVar.calc_staterror)
            self.yerr = self.yerr/y[1]
            warning("The displayed errorbars have been supplied with the data or calculated using chi2xspecvar; the errors are not used in fits with %s" % stat.name)
        else:
            staterr = data.get_yerr(True,stat.calc_staterror)
            self.yerr = staterr/y[1]

        self.ylabel = 'Data / Model'
        self.title = _make_title('Ratio of Data to Model', data.name)

    def plot(self, overplot=False, clearwindow=True):
        Plot.plot(self, self.x, self.y, self.yerr, self.xerr, self.title,
                  self.xlabel, self.ylabel, overplot, clearwindow)

    
class RatioContour(ModelContour):
    "Derived class for creating 2D ratio contours (data:model)"
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

    def contour(self, overcontour=False, clearwindow=True):
        Contour.contour(self, self.x0, self.x1, self.y, levels=self.levels,
                        title=self.title, xlabel=self.xlabel,
                        ylabel=self.ylabel,
                        overcontour=overcontour, clearwindow=clearwindow)



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

        if not state.has_key('stat'):
            self.__dict__['stat'] = None

        if not state.has_key('parval'):
            self.__dict__['parval'] = None

        if not state.has_key('numcores'):
            self.__dict__['numcores'] = None

    def __str__(self):
        x = self.x
        if self.x is not None:
            x = numpy.array2string(self.x, separator=',', precision=4, suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)

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


    def prepare(self, min=None, max=None, nloop=20,
                delv=None, fac=1, log=False, numcores=None):
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

        if( not numpy.isscalar( self.min ) or
            not numpy.isscalar( self.max ) ):
            raise ConfidenceErr('badarg','Parameter limits', 'scalars')

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
            self.x = numpy.linspace(self.min,self.max,self.nloop)
        else:
            eps = numpy.finfo(numpy.float32).eps
            self.x = numpy.arange(self.min,self.max+self.delv-eps,self.delv)

        x = self.x
        if self.log:
            if self.max <= 0.0 or self.min <= 0.0:
                raise ConfidenceErr('badarg', 'Log scale',
                                    'on positive boundaries')
            self.max = numpy.log10(self.max)
            self.min = numpy.log10(self.min)

            x = numpy.linspace(self.min,self.max,len(x))

        return x


    def calc(self, fit, par):
        if type(fit.stat) in (LeastSq,):
            raise ConfidenceErr('badargconf', fit.stat.name)

    def plot(self, overplot=False, clearwindow=True):
        if self.log:
            self.plot_prefs['xlog']=True

        Plot.plot(self, self.x, self.y, title=self.title, xlabel=self.xlabel,
                  ylabel=self.ylabel, overplot=overplot,
                  clearwindow=clearwindow)

        if self.stat is not None:
            Plot.hline(self.stat, linecolor="green", linestyle="dash",
                       linewidth=1.5,overplot=True, clearwindow=False)

        if self.parval is not None:
            Plot.vline(self.parval, linecolor="orange", linestyle="dash",
                       linewidth=1.5,overplot=True, clearwindow=False)

        if self.log:
            self.plot_prefs['xlog']=False


class Confidence2D(DataContour, Point):

    contour_prefs = backend.get_confid_contour_defaults()
    point_prefs = backend.get_confid_point_defaults()

    def __init__(self):
        self.min = None
        self.max = None
        self.nloop = (10,10)
        self.fac = 4
        self.delv = None
        self.log = (False,False)
        self.sigma = (1,2,3)
        self.parval0 = None
        self.parval1 = None
        self.stat = None
        self.numcores = None
        DataContour.__init__(self)

    def __setstate__(self, state):
        self.__dict__.update(state)

        if not state.has_key('stat'):
            self.__dict__['stat'] = None

        if not state.has_key('numcores'):
            self.__dict__['numcores'] = None


    def __str__(self):
        x0 = self.x0
        if self.x0 is not None:
            x0 = numpy.array2string(self.x0, separator=',', precision=4, suppress_small=False)

        x1 = self.x1
        if self.x1 is not None:
            x1 = numpy.array2string(self.x1, separator=',', precision=4, suppress_small=False)

        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)

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


    def prepare(self, min=None, max=None, nloop=(10,10),
                delv=None, fac=4, log=(False,False),
                sigma=(1,2,3), levels=None, numcores=None):
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
                    self.sigma[i]/numpy.sqrt(2.))))
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

            # check user limits for errors
            if( numpy.isscalar( self.min ) or
                numpy.isscalar( self.max ) ):
                raise ConfidenceErr('badarg','Parameter limits', 'a list')

            for i in [0,1]:
                v = (self.max[i] + self.min[i]) / 2.
                dv = numpy.fabs(v - self.min[i])
                self.min[i] = v - self.fac * dv
                self.max[i] = v + self.fac * dv

        hmin = numpy.array([ par0.min, par1.min ])
        hmax = numpy.array([ par0.max, par1.max ])

        for i in [0,1]:
            # check user limits for errors
            if( numpy.isscalar( self.min ) or
                numpy.isscalar( self.max ) ):
                raise ConfidenceErr('badarg','Parameter limits', 'a list')

            if self.min[i] >= self.max[i]:
                raise ConfidenceErr('badlimits')

            if numpy.isscalar( self.nloop ) or self.nloop[i] <= 1:
                raise ConfidenceErr('badarg', 'Nloop parameter','a list with' +
                                      ' elements > 1')

            if self.min[i] < hmin[i]:
                self.min[i] = hmin[i]
            if self.max[i] > hmax[i]:
                self.max[i] = hmax[i]


        if self.delv is None:
            self.x0 = numpy.linspace(self.min[0],self.max[0],self.nloop[0])
            self.x1 = numpy.linspace(self.min[1],self.max[1],self.nloop[1])

        else:
            eps = numpy.finfo(numpy.float32).eps
            self.x0 = numpy.arange(self.min[0],self.max[0]+self.delv[0]-eps,
                                   self.delv[0])
            self.x1 = numpy.arange(self.min[1],self.max[1]+self.delv[1]-eps,
                                   self.delv[1])

        #x = numpy.array([self.x0, self.x1])
        x = [self.x0, self.x1]

        self.x0, self.x1 = numpy.meshgrid(self.x0, self.x1)
        self.x0 = self.x0.ravel()
        self.x1 = self.x1.ravel()

        for i in [0,1]:
            if self.log[i]:
                if self.max[i] <= 0.0 or self.min[i] <= 0.0:
                    raise ConfidenceErr('badarg', 'Log scale',
                                        'on positive boundaries')
                self.max[i] = numpy.log10(self.max[i])
                self.min[i] = numpy.log10(self.min[i])
                x[i] = numpy.linspace(self.min[i],self.max[i],len(x[i]))

        x0, x1 = numpy.meshgrid(x[0], x[1])
        return numpy.array([x0.ravel(), x1.ravel()]).T


    def calc(self, fit, par0, par1):
        if type(fit.stat) in (LeastSq,):
            raise ConfidenceErr('badargconf', fit.stat.name)


    def contour(self, overplot=False, clearwindow=True):

        if self.log[0]:
            self.contour_prefs['xlog']=True
        if self.log[1]:
            self.contour_prefs['ylog']=True

        Contour.contour(self, self.x0, self.x1, self.y, levels=self.levels,
                        title=self.title, xlabel=self.xlabel,
                        ylabel=self.ylabel, overcontour=overplot,
                        clearwindow=clearwindow)
        Point.point(self, self.parval0, self.parval1,
                    overplot=True, clearwindow=False)

        if self.log[0]:
            self.contour_prefs['xlog']=False
        if self.log[1]:
            self.contour_prefs['ylog']=False


class IntervalProjection(Confidence1D):

    def __init__(self):
        self.fast = True
        Confidence1D.__init__(self)


    def prepare(self, fast=True, min=None, max=None, nloop=20,
                delv=None, fac=1, log=False, numcores=None):
        self.fast = fast
        Confidence1D.prepare(self, min, max, nloop, delv, fac, log, numcores)


    def calc(self, fit, par, methoddict=None):
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
                    warning("Setting optimization to " + fit.method.name
                            + " for interval projection plot")
            else:
                if (type(fit.method) is not LevMar):
                    fit.method = methoddict['levmar']
                    warning("Setting optimization to " + fit.method.name
                            + " for interval projection plot")

        xvals = self._interval_init(fit, par)
        oldpars = fit.model.thawedpars
        par.freeze()

        def eval_proj(val):
            if self.log:
                val = numpy.power(10, val)
            par.val = val
            if len(thawed) > 1:
                r = fit.fit()
                return r.statval
            return fit.calc_stat()

        try:
            fit.model.startup()

            # store the class methods for startup and teardown
            # these calls are unnecessary for every fit
            startup = fit.model.startup
            fit.model.startup = lambda : None
            teardown = fit.model.teardown
            fit.model.teardown = lambda : None

            self.y = numpy.asarray(parallel_map(eval_proj, xvals,
                                                self.numcores))

        finally:
            # Set back data that we changed
            par.thaw()

            fit.model.startup = startup
            fit.model.teardown = teardown

            fit.model.teardown()
            fit.model.thawedpars = oldpars
            fit.method = oldfitmethod


class IntervalUncertainty(Confidence1D):

    def calc(self, fit, par, methoddict=None):
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

        def eval_uncert(val):
            if self.log:
                val = numpy.power(10, val)
            par.val = val
            return fit.calc_stat()

        try:
            fit.model.startup()
            self.y = numpy.asarray(parallel_map(eval_uncert, xvals,
                                                self.numcores))

        finally:
            # Set back data that we changed
            for i in thawed:
                i.thaw()
            fit.model.teardown()
            fit.model.thawedpars = oldpars


class RegionProjection(Confidence2D):

    def __init__(self):
        self.fast = True
        Confidence2D.__init__(self)


    def prepare(self, fast=True, min=None, max=None, nloop=(10,10),
                delv=None, fac=4, log=(False,False),
                sigma=(1,2,3), levels=None, numcores=None):
        self.fast = fast
        Confidence2D.prepare(self, min, max, nloop, delv, fac, log,
                             sigma, levels, numcores)


    def calc(self, fit, par0, par1, methoddict=None):
        self.title='Region-Projection'

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
                    warning("Setting optimization to " + fit.method.name
                            + " for region projection plot")
            else:
                if (type(fit.method) is not LevMar):
                    fit.method = methoddict['levmar']
                    warning("Setting optimization to " + fit.method.name
                            + " for region projection plot")


        def eval_proj(pars):
            for ii in [0,1]:
                if self.log[ii]:
                    pars[ii] = numpy.power(10, pars[ii])
            (par0.val, par1.val) = pars
            if len(thawed) > 2:
                r = fit.fit()
                return r.statval
            return fit.calc_stat()

        oldpars = fit.model.thawedpars

        try:
            fit.model.startup()
            
            # store the class methods for startup and teardown
            # these calls are unnecessary for every fit
            startup = fit.model.startup
            fit.model.startup = lambda : None
            teardown = fit.model.teardown
            fit.model.teardown = lambda : None

            grid = self._region_init(fit, par0, par1)

            par0.freeze()
            par1.freeze()

            self.y = numpy.asarray(parallel_map(eval_proj, grid,
                                                self.numcores))

        finally:
            # Set back data after we changed it
            par0.thaw()
            par1.thaw()

            fit.model.startup = startup
            fit.model.teardown = teardown

            fit.model.teardown()
            fit.model.thawedpars = oldpars
            fit.method = oldfitmethod


class RegionUncertainty(Confidence2D):

    def calc(self, fit, par0, par1, methoddict=None):
        self.title='Region-Uncertainty'

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

        def eval_uncert(pars):
            for ii in [0,1]:
                if self.log[ii]:
                    pars[ii] = numpy.power(10, pars[ii])
            (par0.val, par1.val) = pars
            return fit.calc_stat()

        oldpars = fit.model.thawedpars

        try:
            fit.model.startup()

            grid = self._region_init(fit, par0, par1)

            for i in thawed:
                i.freeze()

            self.y = numpy.asarray(parallel_map(eval_uncert, grid,
                                                self.numcores))

        finally:
            # Set back data after we changed it
            for i in thawed:
                i.thaw()
            fit.model.teardown()
            fit.model.thawedpars = oldpars
