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

import numpy

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
    scales = {}
    for bkg_id in data.background_ids:
        if bkg_id not in bkg_srcs:
            raise ModelErr('nobkg', bkg_id, id)

        scales[bkg_id] = data.get_background_scale(bkg_id, group=False)

    # Group by the scale factor: numpy arrays do not hash,
    # so we need a way to handle them. For now I am
    # going to assume that each vector is unique (i.e.
    # it is not worth combining components).
    #
    # TODO: combine with the above loop.
    #
    flattened = defaultdict(list)
    vectors = []
    for key, scale in scales.items():
        if numpy.isscalar(scale):
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

        # Use a tablemodel to contain the array values:
        # this makes the string output nicer, but we
        # need to deal with the namespace for these
        # models, and cleaning up after ourselves. We actually
        # use a specialized model class, rather than a tablemodel,
        # since it has to just return the full data (in
        # the same way the RSP models do).
        #
        tbl = ScaleArray('scale{}_{}'.format(id, key))
        tbl.y = scale

        # Add the model to the global symbol list, so
        # users can interrogate it. The issues are:
        #  - name clash: it's okay to overwrite models created
        #    here, but what happens if a user has created a
        #    model with this name?
        #  - how do we clean up after the data goes out of
        #    scope (or the background scaling changes)
        #
        session._add_model_component(tbl)

        model += tbl * resp(bkg_srcs[key])

    return model
