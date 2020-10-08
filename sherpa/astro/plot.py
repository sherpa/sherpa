#
#  Copyright (C) 2010, 2015, 2016, 2019, 2020  Smithsonian Astrophysical Observatory
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
Classes for plotting, analysis of astronomical data sets
"""

from sherpa.astro.data import DataPHA
from sherpa.plot import DataPlot, ModelPlot, FitPlot, DelchiPlot, ResidPlot, \
    RatioPlot, ChisqrPlot, HistogramPlot, backend, Histogram
from sherpa.plot import ComponentSourcePlot as _ComponentSourcePlot
import sherpa.plot
from sherpa.astro.utils import bounds_check
from sherpa.utils.err import PlotErr, IOErr
from sherpa.utils import parse_expr, dataspace1d, histogram1d, filter_bins
from numpy import iterable, array2string, asarray
import logging

warning = logging.getLogger(__name__).warning

__all__ = ('DataPHAPlot', 'SourcePlot', 'ComponentModelPlot',
           'ComponentSourcePlot', 'ARFPlot', 'BkgDataPlot',
           'BkgFitPlot', 'BkgSourcePlot', 'BkgDelchiPlot', 'BkgResidPlot',
           'BkgRatioPlot', 'BkgChisqrPlot',
           'OrderPlot', 'ModelHistogram', 'BkgModelHistogram',
           'FluxHistogram', 'EnergyFluxHistogram', 'PhotonFluxHistogram')


def to_latex(txt):
    """Add any backend-specific markup to indicate LaTeX.

    Parameters
    ----------
    txt : str
       The LaTeX input (the contents are going to be interpreted
       as LaTeX so do not send in text that is not to be converted).

    Returns
    -------
    out : str
       The LaTeX text including the markup to indicate to the
       plotting backend to display as LaTeX.
    """

    return backend.get_latex_for_string(txt)


class DataPHAPlot(sherpa.plot.DataHistogramPlot):
    """Plot a PHA dataset."""

    histo_prefs = sherpa.plot.get_data_hist_prefs()

    def prepare(self, data, stat=None):

        # Need a better way of accessing the binning of the data.
        # Maybe to_plot should return the lo/hi edges as a pair
        # here.
        #
        (_, self.y, self.yerr, self.xerr, self.xlabel,
         self.ylabel) = data.to_plot()

        if stat is not None:
            yerrorbars = self.histo_prefs.get('yerrorbars', True)
            self.yerr = sherpa.plot.calculate_errors(data, stat, yerrorbars)

        self.title = data.name

        # Get the X axis data.
        #
        if data.units != 'channel':
            elo, ehi = data._get_ebins(group=False)
        else:
            elo, ehi = (data.channel, data.channel + 1.)

        self.xlo = data.apply_filter(elo, data._min)
        self.xhi = data.apply_filter(ehi, data._max)
        if data.units == 'wavelength':
            self.xlo = data._hc / self.xlo
            self.xhi = data._hc / self.xhi


class ModelPHAHistogram(HistogramPlot):
    """Plot a model for a PHA dataset.

    The filtering and grouping from the PHA datset
    are used to create the bins for the model.
    """

    histo_prefs = backend.get_model_histo_defaults()

    def __init__(self):
        HistogramPlot.__init__(self)
        self.title = 'Model'

    def prepare(self, data, model, stat=None):

        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        (_, self.y, _, _,
         self.xlabel, self.ylabel) = data.to_plot(yfunc=model)
        self.y = self.y[1]

        if data.units != 'channel':
            elo, ehi = data._get_ebins(group=False)
        else:
            elo, ehi = (data.channel, data.channel + 1.)

        self.xlo = data.apply_filter(elo, data._min)
        self.xhi = data.apply_filter(ehi, data._max)
        if data.units == 'wavelength':
            self.xlo = data._hc / self.xlo
            self.xhi = data._hc / self.xhi


class ModelHistogram(ModelPHAHistogram):
    """Plot a model for a PHA dataset with no grouping.

    The model is drawn at the native resolution of
    the instrument response, or ungrouped channels.
    The model is filtered to the minimum and maximum
    of the PHA dataset but any ignored ranges within
    this are ignored.
    """

    def prepare(self, data, model, stat=None):

        old_filter = parse_expr(data.get_filter())
        old_group = data.grouped
        new_filter = parse_expr(data.get_filter(group=False))
        try:
            if old_group:
                data.ungroup()
                for interval in new_filter:
                    data.notice(*interval)

            super().prepare(data, model, stat=stat)

        finally:
            if old_group:
                data.ignore()
                data.group()
                for interval in old_filter:
                    data.notice(*interval)


class SourcePlot(HistogramPlot):
    """Create 1D plots of unconcolved model values.

    Attributes
    ----------
    histo_prefs : dict
       The preferences for the plot.
    xlo, xhi : array_like
       The lower and upper edges for each bin (the independent variable).
    y : array_like
       The Y value for each point (the model value).
    xlabel, ylabel, title : str
       Plot labels.

    """

    histo_prefs = backend.get_model_histo_defaults()

    def __init__(self):
        self.units = None
        self.mask  = None
        HistogramPlot.__init__(self)
        self.title = 'Source'

    def prepare(self, data, src, lo=None, hi=None):
        # Note: src is source model before folding
        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        lo, hi = bounds_check(lo, hi)

        self.units = data.units
        if self.units == "channel":
            warning("Channel space is unappropriate for the PHA unfolded" +
                    " source model,\nusing energy.")
            self.units = "energy"

        self.xlabel = data.get_xlabel()
        self.title  = 'Source Model of %s' % data.name
        self.xlo, self.xhi = data._get_indep(filter=False)
        self.mask = filter_bins((lo,), (hi,), (self.xlo,))
        self.y = src(self.xlo, self.xhi)
        prefix_quant = 'E'
        quant = 'keV'

        if self.units == "wavelength":
            # No other labels use the LaTeX forms for lambda and
            # Angstrom, so use the text version here too.
            # prefix_quant = to_latex('\\lambda')
            # quant = to_latex('\\AA')
            prefix_quant = 'lambda'
            quant = 'Angstrom'
            (self.xlo, self.xhi) = (self.xhi, self.xlo)

        xmid = abs(self.xhi - self.xlo)

        sqr = to_latex('^2')

        self.xlabel = '%s (%s)' % (self.units.capitalize(), quant)
        self.ylabel = '%s  Photons/sec/cm' + sqr + '%s'

        if data.plot_fac == 0:
            self.y /= xmid
            self.ylabel = self.ylabel % ('f(%s)' % prefix_quant,
                                         '/%s ' % quant)

        elif data.plot_fac == 1:
            self.ylabel = self.ylabel % ('%s f(%s)' % (prefix_quant,
                                                       prefix_quant), '')

        elif data.plot_fac == 2:
            self.y *= xmid
            self.ylabel = self.ylabel % ('%s%s f(%s)' % (prefix_quant, sqr,
                                                         prefix_quant),
                                         ' %s ' % quant)
        else:
            raise PlotErr('plotfac', 'Source', data.plot_fac)

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        xlo = self.xlo
        xhi = self.xhi
        y = self.y

        if self.mask is not None:
            xlo = self.xlo[self.mask]
            xhi = self.xhi[self.mask]
            y = self.y[self.mask]

        Histogram.plot(self, xlo, xhi, y, title=self.title,
                       xlabel=self.xlabel, ylabel=self.ylabel,
                       overplot=overplot, clearwindow=clearwindow,
                       **kwargs)


class ComponentModelPlot(_ComponentSourcePlot, ModelHistogram):

    histo_prefs = backend.get_component_histo_defaults()

    def __init__(self):
        ModelHistogram.__init__(self)

    def __str__(self):
        return ModelHistogram.__str__(self)

    def prepare(self, data, model, stat=None):
        ModelHistogram.prepare(self, data, model, stat)
        self.title = 'Model component: %s' % model.name

    def _merge_settings(self, kwargs):
        return sherpa.plot.merge_settings(self.histo_prefs, kwargs)

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        ModelHistogram.plot(self, overplot=overplot,
                            clearwindow=clearwindow, **kwargs)


class ComponentSourcePlot(_ComponentSourcePlot, SourcePlot):

    histo_prefs = backend.get_component_histo_defaults()

    def __init__(self):
        SourcePlot.__init__(self)

    def __str__(self):
        return SourcePlot.__str__(self)

    def prepare(self, data, model, stat=None):
        SourcePlot.prepare(self, data, model)
        self.title = 'Source model component: %s' % model.name

    def _merge_settings(self, kwargs):
        return sherpa.plot.merge_settings(self.histo_prefs, kwargs)

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        SourcePlot.plot(self, overplot=overplot,
                        clearwindow=clearwindow, **kwargs)


class ARFPlot(HistogramPlot):
    """Create plots of the ancillary response file (ARF).

    Attributes
    ----------
    histo_prefs : dict
       The preferences for the plot.
    xlo, xhi : array_like
       The lower and upper edges of each bin.
    y : array_like
       The effective area (ARF) value for the bin.
    xlabel, ylabel, title : str
       Plot labels.

    """

    histo_prefs = backend.get_model_histo_defaults()

    def prepare(self, arf, data=None):
        """Fill the fields given the ARF.

        Parameters
        ----------
        arf :
           The ARF to plot
        data : DataPHA instance, optional
           The `units` attribute of this object is used
           to determine whether the X axis should be
           in Angstrom instead of KeV (the default).

        """
        self.xlo = arf.energ_lo
        self.xhi = arf.energ_hi
        self.y = arf.specresp

        self.title = arf.name
        self.xlabel = arf.get_xlabel()
        self.ylabel = arf.get_ylabel()

        if data is not None:
            if not isinstance(data, DataPHA):
                raise PlotErr('notpha', data.name)
            if data.units == "wavelength":
                self.xlabel = 'Wavelength (Angstrom)'
                self.xlo = data._hc / self.xlo
                self.xhi = data._hc / self.xhi


class BkgDataPlot(DataPHAPlot):
    "Derived class for creating plots of background counts"

    # Is this derived class worth it?
    pass


# was BkgModelPlot; this is not a good class name
class BkgModelPHAHistogram(ModelPHAHistogram):
    """Plot a background model for a PHA dataset.

    The filtering and grouping from the background of
    the PHA datset are used to create the bins for the model.
    """

    def __init__(self):
        ModelPHAHistogram.__init__(self)
        self.title = 'Background Model Contribution'


class BkgModelHistogram(ModelHistogram):
    """Plot a background model for a PHA dataset with no grouping.

    The model is drawn at the native resolution of
    the instrument response, or ungrouped channels.
    The model is filtered to the minimum and maximum
    of the background to the PHA dataset but any ignored
    ranges within this are ignored.
    """

    def __init__(self):
        ModelPHAHistogram.__init__(self)
        self.title = 'Background Model Contribution'


class BkgFitPlot(FitPlot):
    "Derived class for creating plots of background counts with fitted model"
    def __init__(self):
        FitPlot.__init__(self)


class BkgDelchiPlot(DelchiPlot):
    "Derived class for creating background plots of 1D delchi chi ((data-model)/error)"
    def __init__(self):
        DelchiPlot.__init__(self)


class BkgResidPlot(ResidPlot):
    "Derived class for creating background plots of 1D residual (data-model)"
    def __init__(self):
        ResidPlot.__init__(self)

    def prepare(self, data, model, stat):
        ResidPlot.prepare(self, data, model, stat)
        self.title = 'Residuals of %s - Bkg Model' % data.name


class BkgRatioPlot(RatioPlot):
    "Derived class for creating background plots of 1D ratio (data:model)"
    def __init__(self):
        RatioPlot.__init__(self)

    def prepare(self, data, model, stat):
        RatioPlot.prepare(self, data, model, stat)
        self.title = 'Ratio of %s : Bkg Model' % data.name


class BkgChisqrPlot(ChisqrPlot):
    "Derived class for creating background plots of 1D chi**2 ((data-model)/error)**2"
    def __init__(self):
        ChisqrPlot.__init__(self)


class BkgSourcePlot(SourcePlot):
    "Derived class for plotting the background unconvolved source model"
    def __init__(self):
        SourcePlot.__init__(self)


class OrderPlot(ModelHistogram):
    """
    Derived class for creating plots of the convolved source model using
    selected multiple responses
    """

    def __init__(self):
        self.orders = None
        self.colors = None
        self.use_default_colors = True
        ModelHistogram.__init__(self)

    def prepare(self, data, model, orders=None, colors=None):
        self.orders = data.response_ids

        if orders is not None:
            if iterable(orders):
                self.orders = list(orders)
            else:
                self.orders = [orders]

        if colors is not None:
            self.use_default_colors = False
            if iterable(colors):
                self.colors = list(colors)
            else:
                self.colors = [colors]
        else:
            self.colors = []
            top_color = '0xffffff'
            bot_color = '0x0000bf'
            num = len(self.orders)
            jump = (int(top_color, 16) - int(bot_color, 16)) // (num + 1)
            for order in self.orders:
                self.colors.append(top_color)
                top_color = hex(int(top_color, 16) - jump)

        if not self.use_default_colors and len(colors) != len(orders):
            raise PlotErr('ordercolors', len(orders), len(colors))

        old_filter = parse_expr(data.get_filter())
        old_group = data.grouped

        try:
            if old_group:
                data.ungroup()
                for interval in old_filter:
                    data.notice(*interval)

            self.xlo = []
            self.xhi = []
            self.y = []
            (xlo, y, yerr, xerr,
             self.xlabel, self.ylabel) = data.to_plot(model)
            y = y[1]
            if data.units != 'channel':
                elo, ehi = data._get_ebins(group=False)
                xlo = data.apply_filter(elo, data._min)
                xhi = data.apply_filter(ehi, data._max)
                if data.units == 'wavelength':
                    xlo = data._hc / xlo
                    xhi = data._hc / xhi
            else:
                xhi = xlo + 1.

            for order in self.orders:
                self.xlo.append(xlo)
                self.xhi.append(xhi)
                if len(data.response_ids) > 2:
                    if order < 1 or order > len(model.rhs.orders):
                        raise PlotErr('notorder', order)
                    y = data.apply_filter(model.rhs.orders[order - 1])
                    y = data._fix_y_units(y, True)
                    if data.exposure:
                        y = data.exposure * y
                self.y.append(y)

        finally:
            if old_group:
                data.ignore()
                data.group()
                for interval in old_filter:
                    data.notice(*interval)

        self.title = 'Model Orders %s' % str(self.orders)

        if len(self.xlo) != len(self.y):
            raise PlotErr("orderarrfail")

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        default_color = self.histo_prefs['linecolor']
        count = 0
        for xlo, xhi, y, color in \
                zip(self.xlo, self.xhi, self.y, self.colors):
            if count != 0:
                overplot = True
                self.histo_prefs['linecolor'] = color

            # Note: the user settings are sent to each plot
            Histogram.plot(self, xlo, xhi, y, title=self.title,
                           xlabel=self.xlabel, ylabel=self.ylabel,
                           overplot=overplot, clearwindow=clearwindow,
                           **kwargs)
            count += 1

        self.histo_prefs['linecolor'] = default_color


# TODO: we should probably derive from a histogram plot that
#       has less PHA heritage
#
class FluxHistogram(ModelHistogram):
    "Derived class for creating 1D flux distribution plots"

    def __init__(self):
        self.modelvals = None
        self.clipped = None
        self.flux = None
        ModelHistogram.__init__(self)

    def __str__(self):
        vals = self.modelvals
        if self.modelvals is not None:
            vals = array2string(asarray(self.modelvals), separator=',',
                                precision=4, suppress_small=False)

        clip = self.clipped
        if self.clipped is not None:
            # Could convert to boolean, but it is surprising for
            # anyone trying to access the clipped field
            clip = array2string(asarray(self.clipped), separator=',',
                                precision=4, suppress_small=False)

        flux = self.flux
        if self.flux is not None:
            flux = array2string(asarray(self.flux), separator=',',
                                precision=4, suppress_small=False)

        return '\n'.join(['modelvals = {}'.format(vals),
                          'clipped = {}'.format(clip),
                          'flux = {}'.format(flux),
                          ModelHistogram.__str__(self)])

    def prepare(self, fluxes, bins):
        """Define the histogram plot.

        Parameter
        ---------
        fluxes : numpy array
            The data, stored in a niter by (npar + 2) matrix, where
            each row is an iteration, the first column is the flux for
            that row, the next npar columns are the parameter values,
            and the last column indicates whether the row was clipped
            (1) or not (0).
        bins : int
            The number of bins to split the flux data into.

        """

        fluxes = asarray(fluxes)
        y = fluxes[:, 0]
        self.flux = y
        self.modelvals = fluxes[:, 1:-1]
        self.clipped = fluxes[:, -1]
        self.xlo, self.xhi = dataspace1d(y.min(), y.max(),
                                         numbins=bins + 1)[:2]
        y = histogram1d(y, self.xlo, self.xhi)
        self.y = y / float(y.max())


class EnergyFluxHistogram(FluxHistogram):
    "Derived class for creating 1D energy flux distribution plots"

    def __init__(self):
        FluxHistogram.__init__(self)
        self.title = "Energy flux distribution"
        self.xlabel = "Energy flux (ergs cm{} sec{})".format(
            to_latex('^{-2}'), to_latex('^{-1}'))
        self.ylabel = "Frequency"


class PhotonFluxHistogram(FluxHistogram):
    "Derived class for creating 1D photon flux distribution plots"

    def __init__(self):
        FluxHistogram.__init__(self)
        self.title = "Photon flux distribution"
        self.xlabel = "Photon flux (Photons cm{} sec{})".format(
            to_latex('^{-2}'), to_latex('^{-1}'))
        self.ylabel = "Frequency"
