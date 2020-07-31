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

from sherpa.models.parameter import Parameter
from sherpa.models.model import ArithmeticModel, modelCacher1d

from sherpa.utils.err import ModelErr

__all__ = ('ScaleArray', )


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
