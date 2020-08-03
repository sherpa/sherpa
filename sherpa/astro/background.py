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
from sherpa.models.parameter import Parameter
from sherpa.models.model import ArithmeticModel, modelCacher1d
from sherpa.utils.err import ModelErr


__all__ = ('ScaleArray', 'add_response')


warning = logging.getLogger(__name__).warning


class ScaleArray(ArithmeticModel):
    """Represent a scale array as a Sherpa model.

    The model always return the saved data (multiplied by the
    amplitude), whatever grid it is asked to evaluate. This is
    designed for handling the scaling factor to convert a background
    model to match the source aperture when the value is an array, not
    a scalar.

    Attributes
    ----------
    ampl
        The amplitude of the model.
    y
        The array of values.

    Examples
    --------

    Create a component and set the y attribute to contain the ``vals``
    array. The model always evaluates to the y array (multiplied by
    the ``ampl`` parameter, which defaults to 1):

    >>> sc = ScaleArray('sc')
    >>> sc.y = vals
    >>> all(sc([1, 2, 4, 5]) == vals)
    >>> all(sc([0.2, 0.3, 0.4], [0.3, 0.4, 0.5]) == vals)

    """

    def __init__(self, name='scalearray'):
        self._y = None
        self.ampl = Parameter(name, 'ampl', 1, frozen=True)
        ArithmeticModel.__init__(self, name, (self.ampl,))

    @property
    def y(self):
        """The data to return"""
        return self._y

    @y.setter
    def y(self, value):
        self._y = value

    @modelCacher1d
    def calc(self, p, x0, x1=None, *args, **kwargs):
        """The model data is always returned.

        The input grids are ignored, and the scaled
        model data is returned.
        """

        if self._y is None:
            raise ModelErr("Model data (y) is not set.")

        return p[0] * self._y


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

    # QUS: if this gets used to generate the response for the
    #      background then how does it pick up the correct response
    #      (ie when fit_bkg is used). Or does that get generated
    #      by a different code path?

    resp = session._get_response(id, data)
    if data.subtracted:
        return resp(model)

    bkg_srcs = session._background_sources.get(id, {})
    if len(bkg_srcs) == 0:
        return resp(model)

    # At this point we have background one or more background
    # components that need to be added to the overall model.
    # If the scale factors are all scalars then we can return
    #
    #   resp(model + sum (scale_i * bgnd_i))        [1]
    #
    # but if any are arrays then we have apply the scale factor
    # after applying the response, that is
    #
    #   resp(model) + sum(scale_i * resp(bgnd_i))   [2]
    #
    # This is because the scale values are in channel space,
    # and not the instrument response (i.e. the values used inside
    # the resp call).
    #
    # Note that if resp is not a linear response - for instance,
    # it is a pileup model - then we will not get the correct
    # answer if there's an array value for the scale factor
    # (i.e. equation [2] above). A warning message is created in this
    # case, but it is not treated as an error.
    #
    # For multiple background datasets we can have different models -
    # that is, one for each dataset - but it is expected that the
    # same model is used for all backgrounds (i.e. this is what
    # we 'optimise' for).

    # Identify the scalar and vector scale values for each
    # background dataset, and combine using the model as a key.
    #
    scales_scalar = defaultdict(list)
    scales_vector = defaultdict(list)
    for bkg_id in data.background_ids:
        try:
            bmdl = bkg_srcs[bkg_id]
        except KeyError:
            raise ModelErr('nobkg', bkg_id, id)

        scale = data.get_background_scale(bkg_id, group=False)

        if np.isscalar(scale):
            store = scales_scalar
        else:
            store = scales_vector

        store[bmdl].append(scale)

    # Combine the scalar terms, grouping by the model.
    #
    for mdl, scales in scales_scalar.items():
        scale = sum(scales)
        model += scale * mdl

    # Apply the instrument response.
    #
    model = resp(model)

    if len(scales_vector) == 0:
        return model

    # Warn if a pileup model is being used. The error message here
    # is terrible.
    #
    # Should this be a Python Warning rather than a logged message?
    #
    if isinstance(resp, PileupResponse1D):
        wmsg = "model results for dataset {} ".format(id) + \
                "likely wrong: use of pileup model and array scaling " + \
                "for the background"

        # warnings.warn(wmsg)
        warning(wmsg)

    # Combine the vector terms, grouping by the model. A trick here
    # is that,to make the string version of the model be readable,
    # we add a model to contain the scale values. This model is
    # similar to a tablemodel but is different enough that we use
    # the ScaleArray class.
    #
    # The choice of namespace for the array (i.e. how to keep
    # it from colliding with a user component whilst still fitting
    # into the Schema scheme) is an open question.
    #
    # Another question is how to remove the model once it's not
    # needed.
    #
    nvectors = len(scales_vector)
    for i, (mdl, scales) in enumerate(scales_vector.items(), 1):
        # special case the single-value case
        if nvectors == 1:
            name = 'scale{}'.format(id)
        else:
            name = 'scale{}_{}'.format(id, i)

        tbl = ScaleArray(name)
        tbl.y = sum(scales)

        # register the component before using it in an expression,
        # so that the str method will create 'scalearray.name'
        # rather than just 'name'.
        #
        session._add_model_component(tbl)

        model += tbl * resp(mdl)

    return model
