#
#  Copyright (C) 2010, 2015, 2016, 2019 - 2024
#  Smithsonian Astrophysical Observatory
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

import logging

import numpy as np

from sherpa.astro import hc
from sherpa.astro.data import DataARF, DataPHA, DataRMF
from sherpa.astro.instrument import ARF1D, RMF1D
from sherpa.astro.utils import bounds_check
from sherpa.data import Data1DInt
from sherpa.models.basic import Delta1D
from sherpa.models.model import Model
from sherpa import plot as shplot
from sherpa.utils import parse_expr, dataspace1d, histogram1d, filter_bins, \
    sao_fcmp
from sherpa.utils.err import IOErr, PlotErr

warning = logging.getLogger(__name__).warning

__all__ = ('DataPHAPlot', 'ModelPHAHistogram', 'ModelHistogram',
           'SourcePlot', 'ComponentModelPlot', 'ComponentSourcePlot',
           'RatioPHAPlot', 'ResidPHAPlot', 'DelchiPHAPlot', 'ChisqrPHAPlot',
           'ARFPlot', 'RMFPlot',
           'BkgDataPlot', 'BkgModelPHAHistogram', 'BkgModelHistogram',
           'BkgFitPlot', 'BkgDelchiPlot', 'BkgResidPlot', 'BkgRatioPlot',
           'BkgChisqrPlot', 'BkgSourcePlot',
           'OrderPlot',
           'FluxHistogram', 'EnergyFluxHistogram', 'PhotonFluxHistogram',
           'DataIMGPlot',
           )


# Identify "close-enough" bin edges when plotting histograms
_tol = np.finfo(np.float32).eps


def _check_hist_bins(xlo: np.ndarray,
                     xhi: np.ndarray
                     ) -> tuple[np.ndarray, np.ndarray]:
    """Ensure lo/hi edges that are "close" are merged.

    Ensure that "close-enough" bin edges use the same value.  We do
    this for all bins, even those that are identical, as it's
    easier. The tolerance is taken to be the float32 "eps" setting, as
    this seems to work for the (limited) data sets I've seen. This is
    to fix issue #977.

    Parameters
    ----------
    xlo, xhi : array
        Lower and upper bin boundaries. Typically, ``xlo`` will contain the
        lower boundary and ``xhi`` the upper boundary, but this function can
        deal with situations where that is reversed. Both arrays have to be
        monotonically increasing or decreasing.

    Returns
    -------
    xlo, xhi : array
        xlo and xhi with values that were very close (within numerical
        tolerance) before changed such that they now match exactly.

    Notes
    -----
    Note that this holds even when plotting wavelength values, who
    have xlo/xhi in decreasing order, since the lo/hi values still
    hold.
    """
    if len(xlo) != len(xhi):
        # Not a Sherpa specific error, because this is more for developers.
        raise ValueError('Input arrays must have same length.')
    # Nothing to compare if input arrays are empty.
    if len(xlo) == 0:
        return xlo, xhi

    # Technically idx should be 0 or 1, with no -1 values. We
    # do not enforce this. What we do is to take all bins that
    # appear similar (sao_fcmp==0) and set the xlo[i+1] bin
    # to the xhi[i] value.
    #
    # Deal with xhi <-> xlo switches. Those can occor when converting
    # from energy to wavelength.
    # Deal with reversed order. Can happen when converting from energy
    # to wavelength, or if input PHA is not ordered in increasing energy.
    # But is both are happening at the same time, need to switch twice, which
    # is a no-op. So, we get to use the elusive Python XOR operator.
    if (xlo[0] > xhi[0]) ^ (xhi[0] > xhi[-1]):
        xlo, xhi = xhi, xlo

    # The input arguments may be read-only, so explicitly copy them.
    # It also makes sure we don't change any existing data.
    #
    xlo = xlo.copy()
    xhi = xhi.copy()

    equal = sao_fcmp(xlo[1:], xhi[:-1], _tol)
    idx, = np.where(equal == 0)
    xlo[idx + 1] = xhi[idx]

    return xlo, xhi


def calc_x(data: DataPHA) -> tuple[np.ndarray, np.ndarray]:
    """Calculate the X axis values

    Parameters
    ----------
    data : DataPHA
       The data object.

    Returns
    -------
    xlo, xhi : tuple of ndarray
       The low and high edges of each bin.

    """

    # Get the X axis data.
    #
    if data.units != 'channel':
        elo, ehi = data._get_ebins(group=False)
    else:
        elo, ehi = (data.channel, data.channel + 1.)

    xlo = data.apply_filter(elo, data._min)
    xhi = data.apply_filter(ehi, data._max)
    if data.units == 'wavelength':
        # Should this swap xlo and xhi here?
        xlo = hc / xlo
        xhi = hc / xhi

    return _check_hist_bins(xlo, xhi)


class DataPHAPlot(shplot.DataHistogramPlot):
    """Plot a PHA dataset."""

    histo_prefs = shplot.get_data_hist_prefs()
    "The preferences for the plot."

    def prepare(self, data, stat=None):

        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        # Need a better way of accessing the binning of the data.
        # Maybe to_plot should return the lo/hi edges as a pair
        # here.
        #
        plot = data.to_plot()
        (_, self.y, self.yerr, _, self.xlabel, self.ylabel) = plot

        if stat is not None:
            yerrorbars = self.histo_prefs.get('yerrorbars', True)
            self.yerr = shplot.calculate_errors(data, stat, yerrorbars)

        self.title = data.name

        self.xlo, self.xhi = calc_x(data)


class RatioPHAPlot(shplot.RatioHistogramPlot):
    """Plot ratio for a PHA dataset.

    .. versionadded:: 4.16.1

    """

    def _calc_x(self, data: Data1DInt, model: Model) -> None:
        """Define the xlo and xhi fields"""

        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        self.xlo, self.xhi = calc_x(data)


class ResidPHAPlot(shplot.ResidHistogramPlot):
    """Plot residuals for a PHA dataset.

    .. versionadded:: 4.16.1

    """

    def _calc_x(self, data: Data1DInt, model: Model) -> None:
        """Define the xlo and xhi fields"""

        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        self.xlo, self.xhi = calc_x(data)

    def _change_ylabel(self) -> None:
        # The original code had the y label displaying units rather
        # than 'Data - Model', which is what the super-class sets. So
        # we override the parent behaviour.
        pass


class DelchiPHAPlot(shplot.DelchiHistogramPlot):
    """Plot delchi residuals for a PHA dataset.

    .. versionadded:: 4.16.1

    """

    def _calc_x(self, data: Data1DInt, model: Model) -> None:
        """Define the xlo and xhi fields"""

        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        self.xlo, self.xhi = calc_x(data)


class ChisqrPHAPlot(shplot.ChisqrHistogramPlot):
    """Plot residuals for a PHA dataset.

    .. versionadded:: 4.16.1

    """

    def _calc_x(self, data: Data1DInt, model: Model) -> None:
        """Define the xlo and xhi fields"""

        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        self.xlo, self.xhi = calc_x(data)


class ModelPHAHistogram(shplot.HistogramPlot):
    """Plot a model for a PHA dataset.

    The filtering and grouping from the PHA dataset
    are used to create the bins for the model.
    """

    histo_prefs = shplot.basicbackend.get_model_histo_defaults()

    def __init__(self):
        super().__init__()
        self.title = 'Model'

    def prepare(self, data, model, stat=None):

        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        (_, self.y, _, _,
         self.xlabel, self.ylabel) = data.to_plot(yfunc=model)
        self.y = self.y[1]

        self.xlo, self.xhi = calc_x(data)


class ModelHistogram(ModelPHAHistogram):
    """Plot a model for a PHA dataset with no grouping.

    The model is drawn at the native resolution of the instrument
    response, or ungrouped channels.

    """

    def prepare(self, data, model, stat=None):

        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        # We could fit this into a single try/finally group but
        # it makes it harder to see what is going on so split
        # it out into separate plot types:
        #   - ungrouped
        #   - all data has been masked out (we let prepare
        #     throw any errors or plot no data)
        #   - grouped but no filter (or the filter doesn't
        #     remove any points).
        #   - grouped and filtered
        #
        if not data.grouped or data.mask is False or \
           (data.mask is not True and not data.mask.any()):
            super().prepare(data, model, stat=stat)
            return

        # At this point mask can be True or an array.
        #
        if data.mask is True or data.mask.all():
            try:
                data.ungroup()
                super().prepare(data, model, stat=stat)

            finally:
                data.group()

            return

        # We need to convert the filter expression from grouped
        # to ungrouped, apply it, create the plot, then restore
        # the old expression. Note that we can just copy over
        # the original mask value to restore the old filter.
        #
        old_mask = data.mask.copy()
        new_filter = parse_expr(data.get_filter(group=False))
        try:
            data.ungroup()
            for interval in new_filter:
                data.notice(*interval)

            super().prepare(data, model, stat=stat)

        finally:
            data.group()
            data.mask = old_mask


class SourcePlot(shplot.SourceHistogramPlot):
    """Create PHA plots of unconvolved model values.

    .. versionchanged:: 4.16.1
       The parent class is now SourceHistogramPlot rather than
       HistogramPlot.

    Attributes
    ----------
    xlo, xhi : array_like
       The lower and upper edges for each bin (the independent variable).
    y : array_like
       The Y value for each point (the model value).
    xlabel, ylabel, title : str
       Plot labels.

    """

    histo_prefs = shplot.basicbackend.get_model_histo_defaults()

    def __init__(self):
        self.units = None
        self.mask = None
        super().__init__()

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

        # Note that there is some effort required to match other PHA plots,
        # e.g. model plots, for:
        #  - plot counts vs rates
        #  - what factor to scale the Y axis by
        #
        self.xlabel = data.get_xlabel()
        self.title = f'Source Model of {data.name}'
        self.xlo, self.xhi = data._get_indep(filter=False)

        # Why do we not apply the mask at the end of prepare?
        #
        self.mask = filter_bins((lo,), (hi,), (self.xlo,))

        # The source model is assumed to not contain an instrument model,
        # and so it evaluates the expected number of photons/cm^2/s in
        # each bin.
        #
        self.y = src(self.xlo, self.xhi)
        prefix_quant = 'E'
        quant = 'keV'

        if self.units == "wavelength":
            # No other labels use the LaTeX forms for lambda and
            # Angstrom, so use the text version here too.
            # prefix_quant = shplot.backend.get_latex_for_string('\\lambda')
            # quant = shplot.backend.get_latex_for_string('\\AA')
            prefix_quant = 'lambda'
            quant = 'Angstrom'
            (self.xlo, self.xhi) = (self.xhi, self.xlo)

        # The y values are in photon/cm^2/s and we want
        #
        # - to divide by the bin width to get per keV/A
        # - to convert to counts when (not data.rate)
        # - to handle the factor value
        #
        # This is similar to DataPHA._fix_y_units but does not have to
        # worry about grouping data, and so is simpler.
        #
        self.y /= abs(self.xhi - self.xlo)

        if not data.rate and data.exposure:
            self.y *= data.exposure
            tlabel = ""
        else:
            tlabel = "/sec"

        # We treat plot_fac < 0 as 0.
        #
        if data.plot_fac <= 0:
            pre = f'f({prefix_quant})'
            post = f'/{quant}'

        elif data.plot_fac == 1:
            pre = f'{prefix_quant} f({prefix_quant})'
            post = ''

        elif data.plot_fac > 1:
            pterm = shplot.backend.get_latex_for_string(f'^{data.plot_fac}')
            pre = f'{prefix_quant}{pterm} f({prefix_quant})'
            post = f' {quant}'
            if data.plot_fac > 2:
                pterm = shplot.backend.get_latex_for_string(f'^{data.plot_fac - 1}')
                post += pterm

        scale = (self.xhi + self.xlo) / 2
        self.y *= scale ** data.plot_fac

        sqr = shplot.backend.get_latex_for_string('^2')
        self.xlabel = f'{self.units.capitalize()} ({quant})'
        self.ylabel = f'{pre}  Photons{tlabel}/cm{sqr}{post}'

        # Should self.mask be applied to self.xlo/hi/y here?
        # If so, the plot method could be removed.

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        xlo = self.xlo
        xhi = self.xhi
        y = self.y

        if self.mask is not None:
            xlo = self.xlo[self.mask]
            xhi = self.xhi[self.mask]
            y = self.y[self.mask]

        # We could temporarily over-write self.xlo, self.xhi, self.y
        # which would mean this could call super().plot(), or set up
        # the data in prepare and avoid this method completely.
        #
        shplot.Histogram.plot(self, xlo, xhi, y, title=self.title,
                              xlabel=self.xlabel, ylabel=self.ylabel,
                              overplot=overplot, clearwindow=clearwindow,
                              **kwargs)


# We do not derive from shplot.ComponentModelHistogramPlot to avoid
# confusion over what behavior is wanted.
#
class ComponentModelPlot(ModelHistogram):
    """The component model plot for DataPHA data.

    .. versionchanged:: 4.16.1

       The class no-longer derives from `sherpa.plot.ComponentSourcePlot`.
    """

    histo_prefs = shplot.basicbackend.get_component_histo_defaults()

    def prepare(self, data, model, stat=None):
        super().prepare(data=data, model=model, stat=stat)
        self.title = f'Model component: {model.name}'


# We do not derive from shplot.ComponentSourceHistogramPlot to avoid
# confusion over what behavior is wanted.
#
class ComponentSourcePlot(SourcePlot):
    """The component source plot for DataPHA data.

    .. versionchanged:: 4.16.1

       The class no-longer derives from `sherpa.plot.ComponentSourcePlot`.
    """

    histo_prefs = shplot.basicbackend.get_component_histo_defaults()

    def prepare(self, data, model, stat=None):
        super().prepare(data=data, src=model)
        self.title = f'Source model component: {model.name}'


class ARFPlot(shplot.HistogramPlot):
    """Create plots of the ancillary response file (ARF).

    Attributes
    ----------
    xlo, xhi : array_like
       The lower and upper edges of each bin.
    y : array_like
       The effective area (ARF) value for the bin.
    xlabel, ylabel, title : str
       Plot labels.

    """

    histo_prefs = shplot.basicbackend.get_model_histo_defaults()

    def prepare(self, arf, data=None):
        """Fill the fields given the ARF.

        Parameters
        ----------
        arf :
           The ARF to plot
        data : DataPHA instance, optional
           The `units` attribute of this object is used to determine
           whether the X axis should be in Angstrom instead of KeV
           (the default). If the units setting is "channel" then the
           X axis is shown using keV.

        """

        if not isinstance(arf, (DataARF, ARF1D)):
            raise IOErr(f"data set '{arf.name}' does not contain an ARF")

        self.xlo = arf.energ_lo
        self.xhi = arf.energ_hi
        self.y = arf.specresp

        self.title = arf.name
        self.xlabel = arf.get_xlabel()
        self.ylabel = arf.get_ylabel()

        if data is None:
            return

        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        # There is no sensible X axis when data.units is channel, so
        # leave as an energy. The alternative is to error out.
        #
        if data.units == "wavelength":
            self.xlabel = 'Wavelength (Angstrom)'
            self.xlo = hc / arf.energ_hi
            self.xhi = hc / arf.energ_lo


class RMFPlot(shplot.HistogramPlot):
    """Create plots of the ancillary response file (RMF).

    A full RMF is a matrix that is hard to visualize.
    Here, we select a few specific energies and show
    the response function for those energies as histograms.

    .. versionchanged:: 4.16.1
       The plot will now display with the analysis setting of the data
       argument sent to `prepare`. Previously it would always use
       energies for the X axis. The `energies` setting can be used to
       control what energies are chosen for the plot.

    """

   # Because this derived from NoNewAttributesAfterInit we need to
   # make the attributes here so they can be used in prepare

    xlo = None
    "array_like: The  lower edges of each bin."

    xhi = None
    "array_like: The  upper edges of each bin."

    y  = None
    "array_like: The response for a specific energy or channel."

    energies = None
    """The energies at which to draw the response (in keV).

    If set to None then `n_lines` energies will be selected to span
    the energy range of the response.
    """

    xlabel = ''
    "Label for X axis"

    ylabel = ''
    "Label for Y axis"

    title = ""
    "Title of plot"

    rmf_plot_prefs = shplot.basicbackend.get_rmf_plot_defaults()
    'Plot preferences'

    # TODO: Make that a plot preference
    # How many monochromatic lines to use
    n_lines = 5
    "The number of lines to draw (only used when `energies` is `None`)."

    labels = None
    "List of strings: The labels for each line."

    def _merge_settings(self, kwargs):
        return {**self.rmf_plot_prefs, **kwargs}

    def prepare(self, rmf, data=None):
        """Fill the fields given the RMF.

        The `n_lines` and `energies` fields can be set to control
        what lines are shown.

        .. versionchanged:: 4.16.1
           The `energies` field, when set, is used to select the
           energies to display (if left as `None` then `n_lines`
           energies are chosen automatically).

        Parameters
        ----------
        RMF :
           The RMF to plot
        data : DataPHA instance, optional
           The `units` attribute of this object is used to determine
           whether the X axis should be in Angstrom or channels,
           instead of KeV (the default).

        """

        if not isinstance(rmf, (DataRMF, RMF1D)):
            raise IOErr(f"data set '{rmf.name}' does not contain a RMF")

        # X access: unlike the ARF case it is possible to display in
        # channel units here.
        #
        if data is not None:
            if not isinstance(data, DataPHA):
                raise IOErr('notpha', data.name)

            units = data.units
        elif rmf.e_min is not None:
            units = "energy"
        else:
            units = "channel"

        # Assume that if the units are not channel then the RMF
        # contains the needed data.
        #
        if units == "energy":
            self.xlo = rmf.e_min
            self.xhi = rmf.e_max
            self.xlabel = "Energy (keV)"

        elif units == "wavelength":
            self.xlo = hc / rmf.e_max
            self.xhi = hc / rmf.e_min
            self.xlabel = "Wavelength (Angstrom)"

        else:
            self.xlo = np.arange(rmf.offset, rmf.detchans + rmf.offset)
            self.xhi = self.xlo + 1
            self.xlabel = "Channel"

        # For now let's just create log-spaced energies. There is no way
        # to select an intelligent range, so we can use something that
        # goes across most of the energ_lo/hi range. This is okay for
        # a Chandra ACIS response. The spacing calculation could use
        # wavelength limits for units=wavelength, but is it worth it?
        #
        elo, ehi = rmf.energ_lo, rmf.energ_hi
        l1 = np.log10(elo[0])
        l2 = np.log10(ehi[-1])

        # NOTE: we do not actually record the chosen energies other
        # than in self.labels, which is a lossy transform.
        #
        if self.energies is None:
            if self.n_lines < 1:
                raise ValueError("n_lines must be >= 1")

            dl = (l2 - l1) / (self.n_lines + 1)
            lines = l1 + dl * np.arange(1, self.n_lines + 1)
            energies = np.power(10, lines)

        else:
            # Minimal checks. We do not re-order this list as the user
            # may want the ordering for a particular reason.  This is
            # liable to change.
            #
            energies = [energy for energy in self.energies
                        if elo[0] <= energy < ehi[-1]]
            if len(energies) == 0:
                raise ValueError("energies must be "
                                 f">= {elo[0]} and < {ehi[-1]} keV")

        mdl = Delta1D()

        y = []
        self.labels = []
        for energy in energies:
            mdl.pos = energy
            y.append(rmf.apply_rmf(mdl(elo, ehi)))
            if units == "wavelength":
                self.labels.append(f'{hc / energy:.2g} Angstrom')
            else:
                # Note: use energy labels for channel space
                self.labels.append(f'{energy:.2g} keV')

        # __str__ and similar functions in the superclass
        # expect this to be an array, not a list
        self.y = np.stack(y)
        self.title = rmf.name

    def plot(self, overplot=False, clearwindow=True,
             **kwargs):
        """Plot the data.

        This will plot the data sent to the prepare method.

        Parameters
        ----------
        overplot : bool, optional
           If `True` then add the data to an existing plot, otherwise
           create a new plot.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plotcs)?
        **kwargs
           These values are passed on to the plot backend, and must
           match the names of the keys of the object's
           plot_prefs dictionary.

        See Also
        --------
        prepare, overplot

        """

        y_array = self.y
        labels = self.labels

        # Override the self.y array with each "energy".
        #
        for n, label in enumerate(labels):
            self.y = y_array[n, :]
            super().plot(
                overplot=overplot if n == 0 else True,
                clearwindow=clearwindow if n == 0 else False,
                label=label,
                **kwargs)
        self.y = y_array


class BkgDataPlot(DataPHAPlot):
    "Derived class for creating plots of background counts"

    # Is this derived class worth it?
    pass


# was BkgModelPlot; this is not a good class name
class BkgModelPHAHistogram(ModelPHAHistogram):
    """Plot a background model for a PHA dataset.

    The filtering and grouping from the background of
    the PHA dataset are used to create the bins for the model.
    """

    def __init__(self):
        super().__init__()
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
        super().__init__()
        self.title = 'Background Model Contribution'


class BkgFitPlot(shplot.FitPlot):
    "Derived class for creating plots of background counts with fitted model"
    pass


class BkgDelchiPlot(DelchiPHAPlot):
    """Derived class for creating background plots of PHA delchi chi ((data-model)/error).

    .. versionchanged:: 4.16.1
       The parent class is now DelchiPHAPlot rather than
       DelchiPlot.

    """

    # leave the title as the parent, which is
    # 'Sigma Residuals for <name>'.
    #
    pass


class BkgResidPlot(ResidPHAPlot):
    """Derived class for creating background plots of PHA residual (data-model).

    .. versionchanged:: 4.16.1
       The parent class is now ResidPHAPlot rather than
       ResidPlot.

    """

    def _title(self, data: Data1DInt) -> None:
        self.title = f'Residuals of {data.name} - Bkg Model'


class BkgRatioPlot(RatioPHAPlot):
    """Derived class for creating background plots of PHA ratio (data:model).

    .. versionchanged:: 4.16.1
       The parent class is now RatioPHAPlot rather than
       RatioPlot.

    """

    def _title(self, data: Data1DInt) -> None:
        self.title = f'Ratio of {data.name} : Bkg Model'


class BkgChisqrPlot(ChisqrPHAPlot):
    """Derived class for creating background plots of chi**2 ((data-model)/error)**2.

    .. versionchanged:: 4.16.1
       The parent class is now ChisqrPHAPlot rather than
       ChisqrPlot.

    """

    pass


class BkgSourcePlot(SourcePlot):
    "Derived class for plotting the background unconvolved source model"
    pass


class OrderPlot(ModelHistogram):
    """
    Derived class for creating plots of the convolved source model using
    selected multiple responses
    """

    def __init__(self):
        self.orders = None
        self.colors = None
        self.use_default_colors = True
        super().__init__()

    # Note: this does not accept a stat parameter.
    def prepare(self, data, model, orders=None, colors=None):

        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        self.orders = data.response_ids

        if orders is not None:
            if np.iterable(orders):
                self.orders = list(orders)
            else:
                self.orders = [orders]

        if colors is not None:
            self.use_default_colors = False
            if np.iterable(colors):
                self.colors = list(colors)
            else:
                self.colors = [colors]
        else:
            self.colors = shplot.backend.colorlist(len(self.orders))

        if not self.use_default_colors and len(self.colors) != len(self.orders):
            raise PlotErr('ordercolors', len(self.orders), len(self.colors))

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
            (xlo, y, _, _,
             self.xlabel, self.ylabel) = data.to_plot(model)
            y = y[1]

            # TODO: should this use calc_x? The logic isn't quite the
            # same but that may be a logical error in the following.
            #
            if data.units != 'channel':
                elo, ehi = data._get_ebins(group=False)
                xlo = data.apply_filter(elo, data._min)
                xhi = data.apply_filter(ehi, data._max)
                if data.units == 'wavelength':
                    xlo = hc / xlo
                    xhi = hc / xhi
            else:
                xhi = xlo + 1.

            for order in self.orders:
                self.xlo.append(xlo)
                self.xhi.append(xhi)
                # QUS: why check that response_ids > 2 and not 1 here?
                #
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

        self.title = f'Model Orders {self.orders}'

        if len(self.xlo) != len(self.y):
            raise PlotErr("orderarrfail")

    def plot(self, overplot=False, clearwindow=True, **kwargs):
        default_color = self.histo_prefs['color']
        count = 0
        for xlo, xhi, y, color in \
                zip(self.xlo, self.xhi, self.y, self.colors):
            if count != 0:
                overplot = True
                self.histo_prefs['color'] = color

            # Note: the user settings are sent to each plot
            shplot.Histogram.plot(self, xlo, xhi, y, title=self.title,
                                  xlabel=self.xlabel, ylabel=self.ylabel,
                                  overplot=overplot, clearwindow=clearwindow,
                                  **kwargs)
            count += 1

        self.histo_prefs['color'] = default_color


# TODO: we should probably derive from a histogram plot that
#       has less PHA heritage
#
class FluxHistogram(ModelHistogram):
    "Derived class for creating 1D flux distribution plots"

    _fields: list[str] = ["modelvals", "clipped", "flux"] + \
        ModelHistogram._fields
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    def __init__(self):
        self.modelvals = None
        self.clipped = None
        self.flux = None
        super().__init__()

    def prepare(self, fluxes, bins):
        """Define the histogram plot.

        Parameters
        ----------
        fluxes : numpy array
            The data, stored in a niter by (npar + 2) matrix, where
            each row is an iteration, the first column is the flux for
            that row, the next npar columns are the parameter values,
            and the last column indicates whether the row was clipped
            (1) or not (0).
        bins : int
            The number of bins to split the flux data into.

        """

        fluxes = np.asarray(fluxes)
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
        super().__init__()
        self.title = "Energy flux distribution"
        self.xlabel = f"Energy flux (ergs cm{shplot.backend.get_latex_for_string('^{-2}')} sec{shplot.backend.get_latex_for_string('^{-1}')})"
        self.ylabel = "Frequency"


class PhotonFluxHistogram(FluxHistogram):
    "Derived class for creating 1D photon flux distribution plots"

    def __init__(self):
        super().__init__()
        self.title = "Photon flux distribution"
        self.xlabel = f"Photon flux (Photons cm{shplot.backend.get_latex_for_string('^{-2}')} sec{shplot.backend.get_latex_for_string('^{-1}')})"
        self.ylabel = "Frequency"


class DataIMGPlot(shplot.Image):
    """Class for DataIMG plots.

    .. warning::
        This class is experimental and subject to change in the future.
        Currently, is is only used within _repr_html_ methods.
    """
    xlabel = 'x'
    ylabel = 'y'
    title = ''
    aspect = 'auto'

    # These are filled with array in the prepare stage, but we need
    # to define them here with some value because of NoNewAttributesAfterInit
    y = None
    x0 = None
    x1 = None

    def prepare(self, img):
        # Apply filter and coordinate system
        #
        self.y = img.get_img()

        # extent is left, right, bottom, top and describes the
        # outer-edge of the pixels.
        #
        ny, nx = img.shape
        coord = img.coord
        if coord in ['physical', 'world']:
            x0, y0 = img._logical_to_physical(0.5, 0.5)
            x1, y1 = img._logical_to_physical(nx + 0.5, ny + 0.5)
            lbl = 'physical'
            cdelt = img.sky.cdelt
            self.aspect = 'equal' if cdelt[1] == cdelt[0] else 'auto'

        else:
            x0, x1, y0, y1 = 0.5, nx + 0.5, 0.5, ny + 0.5
            self.aspect = 'equal'
            lbl = 'logical'

        self.x0 = np.linspace(x0, x1, nx, endpoint=True)
        self.x1 = np.linspace(y0, y1, ny, endpoint=True)

        # What is the filtered dataset?
        #
        if img.get_filter_expr() != '':
            x0, x1 = img.get_indep(filter=True)

            x0min, x0max = np.min(x0), np.max(x0)
            x1min, x1max = np.min(x1), np.max(x1)

            # Should add in half cdelt to pad these, but
            # it looks like it isn't necessary.
            # TODO: The old code (before converting it to a class)
            # used set_xlim and set_ylim
            # on the filtered limits, but we don't have that in the
            # general in the backend. We should have an approach
            # that is consistent with the line plots.
            # For now, define here, but don't use.
            filtered = (x0min, x1min, x0max, x1max)

        else:
            filtered = None

        # Other plot classes use something like img.get_xlabel()
        # but we dn't have that here, so for now this is what it is.
        self.xlabel = f'X ({lbl})'
        self.ylabel = f'Y ({lbl})'
        self.title = img.name

    def plot(self, overplot=False, clearwindow=True, **kwargs):

        super().plot(self.x0, self.x1, self.y, title=self.title,
                     xlabel=self.xlabel, ylabel=self.ylabel,
                     aspect=self.aspect, overplot=overplot,
                     clearwindow=clearwindow, **kwargs)
