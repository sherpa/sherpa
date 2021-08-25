#
#  Copyright (C) 2021 MIT
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
import numpy as np
import sherpa
from sherpa.utils.err import DataErr
from sherpa.astro.background import get_response_for_pha

__all__ = ('fake_pha', )


def fake_pha(data, model,
             is_source=True, pileup_model=None,
             add_bkgs=False, bkg_models={}, id=None):
    """Simulate a PHA data set from a model.

    This function replaces the counts in a PHA dataset with simulated counts
    drawn from a model with Poisson noise. For the simulations, all the details
    already set up in the PHA dataset will be used, including the exposure
    time, one or more ARFs and RMFs, area and background scalings,
    grouping, and data quality arrys.

    Including a background component is optional; if requested, the background
    will be a Poisson draw from the average of all backgrounds that have been
    set for the input `sherpa.astro.data.DataPHA`. For each background
    component, the method can use either the PHA distribution in that
    background component or a model that is evaluated using the response
    set for that background component.
    The later case avoids adding extra noise from a Poisson draw from a
    distribution that might already have very few counts in the first place.

    The backgrounds itself are not changed by this function. To simulate
    backgrounds as well as the source spectrum, call this function on the
    source PHA dataset and the background PHA dataset(s) independently.

    Parameters
    ----------
    data : sherpa.astro.data.DataPHA
        The dataset (may be a background dataset).
    model : sherpa.models.model.ArithmeticModel instance
        The model that will be used for simulations.
    is_source : bool
        ``True`` means that the ``model`` does not contain response or
        background components and that these need to be added based on the
        ARF, RMF, and backgrounds set up for the data set. If ``False``, then
        the ``model`` contains any components to describe the instrument
        already.
    pileup_model : None or `sherpa.astro.models.JDPileup` instance
        Pileup Model for the source spectrum
    add_bkgs : bool
        If ``True`` backgrounds are added to the source counts.
    bkg_srcs : dict
        Keys in the dictionary need to be the background ids in the dataset
        ``data``, and the values are the corresponsing source models. For all
        background datasets that are listed in this dictionary, the
        background counts will be simulated based on the model, appropriately
        scaled (for area etc.) and added to the source. The same ``is_source``
        setting as for the source model applies here, i.e. if the source model
        already contains the ARF and the RMF, then the background models should
        do so, too. This setting has no effect if ``add_bkgs=False``.

        For all background ids not listed in this dictionary, the
        counts will be drawn from the PHA data of the background data set.
    id : str
        String with id number if called from UI layer. This is only used for
        certain error messages.

    Examples
    --------
        Estimate the signal from a 5000 second observation using the
        ARF and RMF from "src.arf" and "src.rmf" respectively:

        >>> set_source(1, xsphabs.gal * xsapec.clus)
        >>> gal.nh = 0.12
        >>> clus.kt, clus.abundanc = 4.5, 0.3
        >>> clus.redshift = 0.187
        >>> clus.norm = 1.2e-3
        >>> fake_pha(1, 'src.arf', 'src.rmf', 5000)

        Simulate a 1 mega second observation for the data and model
        from the default data set. The simulated data will include an
        estimated background component based on scaling the existing
        background observations for the source. The simulated data
        set, which has the same grouping as the default set, for
        easier comparison, is created with the 'sim' label and then
        written out to the file 'sim.pi':

        >>> arf = get_arf()
        >>> rmf = get_rmf()
        >>> bkg = get_bkg()
        >>> bscal = get_backscal()
        >>> grp = get_grouping()
        >>> qual = get_quality()
        >>> texp = 1e6
        >>> set_source('sim', get_source())
        >>> fake_pha('sim', arf, rmf, texp, backscal=bscal, bkg=bkg,
        ...          grouping=grp, quality=qual, grouped=True)
        >>> save_pha('sim', 'sim.pi')

    """
    if len(data.response_ids) == 0:
        raise DataErr('normffake', data.name)

    if is_source:
        model = get_response_for_pha(data, model, bkg_srcs={},
                                     pileup_model=pileup_model, id=id)

    # Get one RMF. Hopefully all of them have the same number of
    # channels, but that sanity check should really be done elsewhere.
    rmf0 = data.get_rmf(data.response_ids[0])
    data.channel = np.arange(1, rmf0.detchans + 1)

    # Calculate the source model, and take a Poisson draw based on
    # the source model.  That becomes the simulated data.
    data.counts = sherpa.utils.poisson_noise(data.eval_model(model))

    # Add in background counts:
    #  -- Scale each background properly given data's
    #     exposure time, BACKSCAL and AREASCAL
    #  -- Take average of scaled backgrounds
    #  -- Take a Poisson draw based on the average scaled background
    #  -- Add that to the simulated data counts
    #
    # Adding background counts is OPTIONAL, only done if user sets
    # "bkg" argument to fake_pha.  The reason is that the user could
    # well set a "source" model that does include a background
    # component.  In that case users should have the option to simulate
    # WITHOUT background counts being added in.
    #
    if add_bkgs:
        nbkg = len(data.background_ids)
        b = 0
        for bkg_id in data.background_ids:
            # we do (probably) want to filter and group the scale array
            scale = data.get_background_scale(bkg_id)
            if bkg_id in bkg_models:
                bkg_pha = data.get_background(bkg_id)
                bkg_model = bkg_models[bkg_id]
                if is_source:
                    bkg_model = get_response_for_pha(bkg_pha, bkg_model, id=id)

                # Exposure in background could differ from exposure in
                # source. But need to set here so eval_model works
                # correctly.
                # (At least I think that's how it works.)
                orig_exposure = bkg_pha.exposure
                bkg_pha.exposure = data.exposure
                # No Poisson here because we make a Poisson draw
                # later using the average of all backgrounds
                cts = bkg_pha.eval_model(bkg_model)
                bkg_pha.exposure = orig_exposure
            else:
                cts = data.get_background(bkg_id).counts
            b += scale * cts

        if nbkg > 0:
            b = b / nbkg
            b_poisson = sherpa.utils.poisson_noise(b)
            data.counts = data.counts + b_poisson
