#
#  Copyright (C) 2021, 2023 MIT
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
'''Simulate PHA datasets

The ``fake_pha`` routine is used to create simulated
`sherpa.astro.data.DataPHA` data objects.
'''

from warnings import warn

from sherpa.utils.err import ArgumentTypeErr, DataErr
from sherpa.utils.random import poisson_noise


__all__ = ('fake_pha', )


def fake_pha(data, model,
             is_source=None, pileup_model=None,
             add_bkgs=None, bkg_models=None, id=None,
             method=None, rng=None, include_bkg_data=False):
    """Simulate a PHA data set from a model.

    This function replaces the counts in a PHA dataset with simulated
    counts drawn from a model with Poisson noise. For the simulations,
    all the details already set up in the PHA dataset will be used,
    including the exposure time, one or more ARFs and RMFs, area and
    background scalings, grouping, and data quality arrys.

    .. versionchanged:: 4.16.1
       This routine has seen significant changes, and the is_source,
       pileup_model, add_bkgs, bkg_models, and id arguments are
       no-longer used. The include_bkg_data argument was added.

    .. versionchanged:: 4.16.0
       The method and rng parameters were added.

    Parameters
    ----------
    data : sherpa.astro.data.DataPHA
        The dataset (may be a background dataset).
    model : sherpa.models.model.ArithmeticModel
        The model that will be used for simulations. It must contain
        any background components, appropriately scaled, and include
        the relevant response.
    include_bkg_data : bool, optional
        Should the counts in the background datasets be included when
        calculating the predicted signal? As background datasets are
        often noisy it is generally better to include a model of the
        background in the model argument.
    method : callable or None
        The routine used to simulate the data. If None (the default)
        then sherpa.utils.random.poisson_noise is used, otherwise the
        function must accept a single ndarray and an optional rng
        argument, returning a ndarray of the same shape as the input.
    rng : numpy.random.Generator, numpy.random.RandomState, or None, optional
        Determines how random numbers are created. If set to None then
        the routines from `numpy.random` are used, and so can be
        controlled by calling `numpy.random.seed`.

    Notes
    -----

    The model must include the relevant response, and if a background
    model is required then it must be included. The model is evaluated
    to create the predicted counts. If the PHA dataset contains any
    backgrounds and include_bkg_data is set then the counts are added,
    after scaling, to the predicted counts. The predicted counts are
    then passed to the simulate routine to create the simulated data.

    If simulated backgrounds are required then they should be
    simulated with this routine and then added to the faked source
    dataset.

    This routine was significantly changed in Sherpa 4.16.1. To update
    old code, note that the following arguments are no-longer used:
    is_source, pileup_model, add_bkgs, bkg_models and id. The model
    argument must now include any background model terms, and the
    necessary scaling term for the conversion from the background to
    source aperture (a combination of the add_bkgs and
    bkg_models="model" arguments), and also include the response
    necessary to convert to channels and counts (a combination of the
    is_source and pileup_model arguments). Setting the bkg_models
    argument to a PHA dataset is now handled by adding the background
    to the data first, and then setting include_bkg_data to True.

    Examples
    --------

    Simulate the data from an absorbed APEC model, where the response
    is manually created (a "perfect" RMF and the ARF is flat, with a
    break at 3.2 keV), and then use the Sherpa plot objects to display
    the simulated data and model. Note that the exposure time must be
    set (if the model normalization is set). First the imports:

    >>> import numpy as np
    >>> from sherpa.astro.data import DataPHA
    >>> from sherpa.astro.instrument import Response1D, create_arf, create_delta_rmf
    >>> from sherpa.astro.xspec import XSphabs, XSapec

    Create the model to simulate (the response is added later):

    >>> gal = XSphabs()
    >>> clus = XSapec()
    >>> model = gal * clus
    >>> gal.nh = 0.12
    >>> clus.kt = 4.5
    >>> clus.abundanc = 0.3
    >>> clus.redshift = 0.23
    >>> clus.norm = 6.3e-4

    Create the ARF, RMF, and PHA datasets. For this simulation the ARF
    is flat, with a break at 3.2 keV, and the response is "perfect"
    (that is, there is no broadening of the model):

    >>> chans = np.arange(1, 101, dtype=np.int16)
    >>> egrid = np.linspace(1, 5, 101)
    >>> elo, ehi = egrid[:-1], egrid[1:]
    >>> yarf = 100 * (elo < 3.2) + 50 * (elo >= 3.2)
    >>> arf = create_arf(elo, ehi, yarf)
    >>> rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    >>> pha = DataPHA("faked", chans, chans * 0)
    >>> pha.set_response(arf=arf, rmf=rmf)
    >>> pha.exposure = 50000

    The response - a combination of the ARF and RMF, along with the
    exposure time - is applied to the model:

    >>> resp = Response1D(pha)
    >>> full_model = resp(model)

    The counts array of the PHA can now be simulated:

    >>> rng = numpy.random.default_rng()
    >>> fake_pha(pha, full_model, rng=rng)

    The following overlplots the model - used to simulate the data -
    on the simulated data:

    >>> pha.set_analysis("energy")
    >>> from sherpa.astro.plot import DataPHAPlot, ModelPHAHistogram
    >>> dplot = DataPHAPlot()
    >>> dplot.prepare(pha)
    >>> dplot.plot()
    >>> mplot = ModelPHAHistogram()
    >>> mplot.prepare(pha, full_model)
    >>> mplot.overplot()

    """

    # Warn users if they are using the old options.
    #
    if is_source is not None:
        warn("is_source is no-longer used, the model "
             "must contain a response",
             category=DeprecationWarning)

    if pileup_model is not None:
        warn("pileup_model is no-longer used, the model "
             "must contain any needed pileup model",
             category=DeprecationWarning)

    if add_bkgs is not None:
        warn("add_bkgs is no-longer used, as the model "
             "either contains a backgorund term or the "
             "backgound datasets are used",
             category=DeprecationWarning)

    if bkg_models is not None:
        warn("bkg_models is no-longer used, as the model "
             "must contain any background component",
             category=DeprecationWarning)

    if id is not None:
        warn("id is no-longer used",
             category=DeprecationWarning)

    # The assumption is that the response matches the PHA, that is the
    # number of channels matches. There is currently no explicit check
    # of this.
    #
    if len(data.response_ids) == 0:
        raise DataErr('normffake', data.name)

    if method is None:
        method = poisson_noise
    elif not callable(method):
        raise ArgumentTypeErr("badarg", "method", "a callable")

    # Evaluate the model. This assumes the model term contains a
    # response and any scaled background components.
    #
    model_prediction = data.eval_model(model)

    if include_bkg_data:
        # If there are background components then we sum up the data,
        # after scaling to match the source data set, and use that as a
        # source term. Note that get_background_scale will correct for
        # exposure time (as units=counts), scaling factors, and the number
        # of background components.
        #
        for bkg_id in data.background_ids:
            cts = data.get_background(bkg_id).counts
            scale = data.get_background_scale(bkg_id, units="counts")
            model_prediction += scale * cts

    data.counts = method(model_prediction, rng=rng)
