#
#  Copyright (C) 2007, 2020  Smithsonian Astrophysical Observatory
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

"""Support multiple background components for a PHA dataset.

"""

from collections import defaultdict
from functools import reduce
import logging
import operator

import numpy as np

from sherpa.astro.instrument import PileupResponse1D
from sherpa.utils.err import ModelErr


__all__ = ('add_response', )


warning = logging.getLogger(__name__).warning


def add_response(session, id, data, model):
    """Create the response model describing the source and model.

    Include any background components and apply the response
    model for the dataset.

    Parameters
    ----------
    session : sherpa.astro.ui.utils.Session instance
    id : int ot str
        The identifier for the dataset.
    data : sherpa.astro.data.DataPHA instance
        The dataset (may be a background dataset).
    model : sherpa.models.model.ArithmeticModel instance
        The model (without response or background components)
        to match to data.

    Returns
    -------
    fullmodel : sherpa.models.model.ArithmeticModel
        The model including the necessary response models and
        background components.

    """

    resp = session._get_response(id, data)
    if data.subtracted:
        return resp(model)

    bkg_srcs = session._background_sources.get(id, {})
    if len(bkg_srcs) == 0:
        return resp(model)

    # What are the scaling factors for each background component?
    # This gives the value you multply the background model to
    # correct it to the source aperture, exposure time, and
    # area scaling.
    #
    scales = data._get_background_scales()

    # Check we have a background model for each background.
    #
    for key in scales.keys():
        if key not in bkg_srcs:
            raise ModelErr('nobkg', key, id)

    # Group by the scale factor: numpy arrays do not hash,
    # so we need a way to handle them. For now I am
    # going to assume that each vector is unique (i.e.
    # it is not worth combining components).
    # later on th
    # an issue here, in terms of memory or space?
    # I am going to assume it isn't for now, otherwise we could
    # go with something like a hash of the value as a key.
    #
    flattened = defaultdict(list)
    vectors = []
    for key, scale in scales.items():
        if np.isscalar(scale):
            flattened[scale].append(key)
        else:
            vectors.append((key, scale))

    for scale in sorted(flattened.keys(), reverse=True):
        ms = [bkg_srcs[k] for k in flattened[scale]]
        combined = reduce(operator.add, ms)
        model += scale * combined

    model = resp(model)
    if len(vectors) == 0:
        return model

    # Warn if a pileup model is being used. The error message here
    # is terrible.
    #
    # Should this be a Python Warning rather than a logged message?
    #
    if isinstance(resp, PileupResponse1D):
        bkg_ids = ",".join([str(k) for k, _ in vectors])
        wmsg = "model results for dataset {} ".format(id) + \
                "likely wrong: use of pileup model and scaling " + \
                "of bkg_id={}".format(bkg_ids)

        # warnings.warn(wmsg)
        warning(wmsg)

    for key, scale in vectors:
        # NOTE: with a NumPy array we need to
        # say mdl * array and not the other way
        # around, otherwise the mdl will get
        # broadcast to each element of array.
        #
        model += resp(bkg_srcs[key]) * scale

    return model
