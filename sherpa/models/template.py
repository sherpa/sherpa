#
#  Copyright (C) 2011, 2016, 2019 - 2021, 2023, 2024
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

# pylint: disable=invalid-name

"""Create models that interpolate between table models.

A `sherpa.models.basic.TableModel` allows an array of values to be
used as a model, where only the normalization is varied. The
`TemplateModel` class lets users take a number of table models,
associate them with one or more parameters, and the
`InterpolatingTemplateModel` class allows for interpolation between
these models, and hence parameter values. Users are expected to use
the `create_template_model` routine to create the template model.

The aim is to allow one or more parameters to have values associated
with numeric data, and then to vary the model based on the parameter
values. For example, consider the case of two parameters - Temperature
and Height - where there's data for a set of these parameters:

=========== ====== ===================
Temperature Height Data (x,y pairs)
=========== ====== ===================
1           10     (1,10), (2,30), ...
2           10     (1,20), (2,40), ...
2           15     (1,12), (3,46), ...
...         ...    ...
10          50     ...
=========== ====== ===================

These data values can be converted to `sherpa.models.basic.TableModel`
instances and then associated with the parameters with
`create_template_model`. The model can then be evaluated to allow the
Temperature and Height values to be fit to observed data. For users of
:term:`XSPEC`, this provides functionality similar to `XSPEC table
models
<https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html>`_.

Examples
========

The templates used here represent a single parameter, called `alpha`,
where the values are for the independent axis values of
[50, 100, 150, 250]:

===== =============
alpha values
===== =============
5     1, 2, 3, 4
17    3, 7, 4, 8
25    10, 15, 0, 12
===== =============

>>> import numpy as np
>>> from sherpa.models.template import create_template_model
>>> from sherpa.models.basic import TableModel
>>> parnames = ["alpha"]
>>> parvals = np.asarray([[5], [17], [25]])
>>> m1 = TableModel("m1")
>>> m2 = TableModel("m2")
>>> m3 = TableModel("m3")
>>> x = [50, 100, 150, 250]
>>> m1.load(x, [1, 2, 3, 4])
>>> m2.load(x, [3, 7, 4, 8])
>>> m3.load(x, [10, 15, 0, 12])
>>> tmps = [m1, m2, m3]
>>> mdl = create_template_model("model", parnames, parvals, tmps)

The model defaults to the first parameter value:

>>> print(mdl)
model
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   model.alpha  thawed            5            5           25

>>> mdl(x)
array([1., 2., 3., 4.])

If the parameter is changed then the templates are weighted to create
the output:

>>> mdl.alpha = 16
>>> mdl(x)
array([3.7, 7.8, 3.6, 8.4])

The model can also be interpolated onto a different X axis:

>>> mdl([70, 120, 200])
array([5.34, 6.12, 6.  ])

"""

from typing import Callable
import operator

import numpy as np

from sherpa.utils.err import ModelErr

from .parameter import Parameter
from .model import ArithmeticModel, modelCacher1d
from .basic import TableModel

__all__ = ('TemplateModel', 'InterpolatingTemplateModel',
           'KNNInterpolator', 'Template', 'add_interpolator',
           'create_template_model', 'reset_interpolators')


# This is reset by reset_interpolators below.
interpolators: dict[str, tuple[Callable, dict]] = {}


def create_template_model(modelname, names, parvals, templates,
                          template_interpolator_name='default'):
    """Create a TemplateModel model class.

    Parameters
    ----------
    modelname : str
        The name of the template model.
    names : sequence of str
        The parameters names.
    parvals : ndarray
        The parameter values, organised as a 2D ndarray of shape
        (nelem, npars), where nelem is the number of parameter values
        and npars the number of parameters (which must match the names
        parameter).
    templates : sequence of TableModel instances
        The model for each set of parameters (each row of parvals).
        It must match the first dimension of parvals.
    template_interpolator_name : str or None, optional
        The interpolator name. If None then there is no interpolation
        between templates.

    Returns
    -------
    model : TemplateModel instance
        The template model (may be a subclass of TemplateModel if
        temlpate_interpolater_name is not None).

    """

    if parvals.ndim != 2:
        raise ValueError(f"parvals must be 2D, sent {parvals.ndim}D")

    # It is not obvious why we have the data transposed, but leave as
    # is for now.
    #
    if parvals.shape[1] != len(names):
        raise ValueError(f"number of parvals and names do not match: {parvals.shape[1]} vs {len(names)}")

    # Create a list of parameters from input
    pars = []
    for name, parval in zip(names, parvals.T):
        minimum = min(parval)
        maximum = max(parval)
        # Initial parameter value is always first parameter value listed
        par = Parameter(modelname, name, parval[0],
                        min=minimum, max=maximum,
                        hard_min=minimum, hard_max=maximum)
        pars.append(par)

    # Create the templates table from input
    tm = TemplateModel(modelname, pars, parvals, templates)
    if template_interpolator_name is None:
        return tm

    try:
        func, args = interpolators[template_interpolator_name]
    except KeyError:
        raise ModelErr(f"Unknown template_interpolator_name={template_interpolator_name}") from None

    # Do not add keys to the stored dictionary
    args = args.copy()
    args['template_model'] = tm
    args['name'] = modelname
    return func(**args)


# Should this provide direct access to query to better match
# TemplateModel.  Or should we not export the parvals attribute, as we
# don't want to copy that interface?
#
class InterpolatingTemplateModel(ArithmeticModel):
    """Allow parameter interpolation for a TemplateModel.

    This class should not be used directly.

    Parameters
    ----------
    name : str
        The name of the model
    template_model : TemplateModel instance
        The templates to use.

    See Also
    --------
    TemplateModel, KNNInterpolator

    """

    def __init__(self, name, template_model):
        self.template_model = template_model
        self.parvals = template_model.parvals

        # The parameters are copied over from parvals, so we need to
        # add them to the object, as the parent class does not do this.
        #
        for par in template_model.pars:
            self.__dict__[par.name] = par

        ArithmeticModel.__init__(self, name, template_model.pars)

    def fold(self, data):
        """Ensure the templates match the data filtering.

        Parameters
        ----------
        data : sherpa.data.Data instance

        """
        self.template_model.fold(data)

    @modelCacher1d
    def calc(self, p, x0, x1=None, *args, **kwargs):
        interpolated_template = self.interpolate(p, x0)
        return interpolated_template(x0, x1, *args, **kwargs)

    def interpolate(self, point, x_out):
        """Calculate the model combination for the given parameter values.

        If the given point matches a template than that is used, otherwise
        weight the templates in some manner.

        Parameters
        ----------
        point : sequence of float
            The parameter values to use. This must match the size and
            ordering of the parameters of the model.
        x_out : sequence of float
            The x values to use for the model evaluation.

        Returns
        -------
        model : TableModel instance
            The model to use.

        """

        raise NotImplementedError


class KNNInterpolator(InterpolatingTemplateModel):
    """Use K-nearest neighbors interpolation for parameters.

    Select the template to use based on the K-nearest neighbors
    algorithm
    https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm

    Parameters
    ----------
    name : str
        The name of the model
    template_model : TemplateModel instance
        The templates to use.
    k : int or None, optional
        The k to use.
    order : int, optional
        The order to use.

    See Also
    --------
    TemplateModel

    """
    def __init__(self, name, template_model, k=None, order=2):
        if k is None:
            self.k = 2 * template_model.parvals.shape[1]
        else:
            self.k = k

        self.order = order
        InterpolatingTemplateModel.__init__(self, name, template_model)

    def _calc_distances(self, point):
        """What are the distances for the given set of parameters?"""
        distances = []
        for i, t_point in enumerate(self.template_model.parvals):
            dist = np.linalg.norm(point - t_point, self.order)
            distances.append((i, dist))

        return sorted(distances, key=operator.itemgetter(1))

    def interpolate(self, point, x_out):
        distances = self._calc_distances(point)
        if distances[0][1] == 0:
            return self.template_model.templates[distances[0][0]]

        k_distances = distances[:self.k]
        weights = [(idx, 1/np.array(distance)) for idx, distance in k_distances]
        y_out = np.zeros(len(x_out))
        for idx, weight in weights:
            y_out += self.template_model.templates[idx].calc((weight,), x_out)

        sum_weights = sum(1 / weight for (idx, weight) in k_distances)
        y_out /= sum_weights
        tm = TableModel('interpolated')
        tm.load(x_out, y_out)
        return tm


class Template(KNNInterpolator):
    """The Template class.

    See Also
    --------
    KNNInterpolator

    """

    pass


class TemplateModel(ArithmeticModel):
    """Combine TableModel instances.

    Parameters
    ----------
    name : str
        The name of the model
    pars : sequence of sherpa.models.parameter.Parameter
        The parameters.
    parvals : sequence of sequence of float
        The parameter values for each template.
    templates : sequence of TableModel
        The template for each element of parvals.

    See Also
    --------
    TemplateModel, KNNInterpolator

    """

    def __init__(self, name='templatemodel', pars=(), parvals=None,
                 templates=None):

        self.parvals = parvals if parvals is not None else []
        self.templates = templates if templates is not None else []

        if len(self.parvals) != len(self.templates):
            raise ModelErr("Number of parameter values and templates do not match")

        if len(self.parvals) == 0:
            raise ModelErr("parvals is empty or not set")

        self.index = {}

        # The parameters are copied over from parvals, so we need to
        # add them to the object, as the parent class does not do this.
        #
        for par in pars:
            self.__dict__[par.name] = par

        npars = len(pars)
        for ii, parval in enumerate(self.parvals):
            ngot = len(parval)
            if ngot != npars:
                raise ModelErr(f"parvals[{ii}] has length {ngot}, expected {npars}")

            self.index[tuple(parval)] = templates[ii]

        ArithmeticModel.__init__(self, name, pars)
        self.is_discrete = True

    def fold(self, data):
        """Ensure the templates match the data filtering.

        Parameters
        ----------
        data : sherpa.data.Data instance

        """
        for template in self.templates:
            template.fold(data)

    def get_x(self):
        """Return the x (independent) axis for the selected template.

        Returns
        -------
        x : sequence of float

        """
        p = tuple(par.val for par in self.pars)
        template = self.query(p)
        return template.get_x()

    def get_y(self):
        """Return the y (dependent) axis for the selected template.

        Returns
        -------
        y : sequence of float

        """
        p = tuple(par.val for par in self.pars)
        template = self.query(p)
        return template.get_y()

    def query(self, p):
        """Return the template for a given set of parameters.

        Parameters
        ----------
        p : sequence of float
            The parameter values to use. This must match the size and
            ordering of the parameters of the model.

        Returns
        -------
        model : TableModel instance

        Raises
        ------
        ModelErr
            There is no template that matches the given set of parameters.

        """
        try:
            return self.index[tuple(p)]
        except KeyError:
            raise ModelErr("Interpolation of template parameters was disabled for this model, but parameter values not in the template library have been requested. Please use gridsearch method and make sure the sequence option is consistent with the template library") from None

    @modelCacher1d
    def calc(self, p, x0, x1=None, *args, **kwargs):
        table_model = self.query(p)

        # return interpolated the spectrum according to the input grid
        # (x0, [x1])
        return table_model(x0, x1, *args, **kwargs)


def reset_interpolators() -> None:
    """Reset the list of interpolators to the default.

    If the list does not exist then recreate it.

    See Also
    --------
    add_interpolator

    """

    d = globals()
    d["interpolators"] = {
        'default': (Template, {'k': 2, 'order': 2})
    }


def add_interpolator(name: str, interpolator: Callable, **kwargs) -> None:
    """Add the interpolator to the list.

    Parameters
    ----------
    name : str
    interpolator
       An interpolator class.
    **kwargs
       The arguments for the interpolator.

    See Also
    --------
    reset_interpolators

    """

    interpolators[name] = (interpolator, kwargs)


reset_interpolators()
