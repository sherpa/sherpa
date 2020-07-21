#
# Copyright (C) 2015, 2016  Smithsonian Astrophysical Observatory
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
from sherpa.utils import public
from sherpa.utils.logging import config_logger
import sherpa
import re

from . import plot_backend as backend

logger = config_logger(__name__)

ID_STR = '__ID'

try:
    import stk
except:
    logger.warning("could not import stk library. CIAO stack files and syntax will be disabled")


@public
def model_wrapper(func):
    def wrapfunc(self, model):
        """Run a model-setting function for each of the datasets."""
        datasets = self.filter_datasets()
        try:
            # if model is passed as a string
            model = eval(model, globals(), sherpa.astro.ui.__dict__)
        except TypeError:
            pass
        except Exception as exc:
            raise type(exc)('Error converting model "{0}" to a sherpa '
                            'model object: {1}'.format(model, exc))
        for dataset in datasets:
            id_ = dataset['id']
            logger.info('Setting stack model using {0}() for id={1}'.format(
                func.__name__, id_))
            new_model, model_comps = create_stack_model(model, id_)
            func(id_, new_model)
            dataset['model_comps'].update(model_comps)
        return None

    wrapfunc.__name__ = func.__name__
    wrapfunc.__doc__ = func.__doc__
    return wrapfunc


@public
def create_stack_model(model, id_):
    """Create any model components for the given data set.

    Parameters
    ----------
    model
       The model expression to use (a string or a Sherpa model
       expression). The names of the model instances will be
       used to create the dataset-specific names, replacing any
       string matching the template with the data set identifier.
    id_ : string or integer
       The identifier for the data set.

    Returns
    -------
    instance, comps
       The instance for the model expression and a dict containing
       any model components.

    See Also
    --------
    set_template_id

    """
    return _create_stack_model(model, id_)


def _create_stack_model(model, model_id, model_components=None):
    model_comps = model_components or {}

    if hasattr(model, 'parts'):
        # Recursively descend through model and create new parts
        # (as needed) corresponding to the stacked model components.
        newparts = []
        for part in model.parts:
            newpart, new_comps = _create_stack_model(
                part, model_id, model_comps)
            model_comps.update(new_comps)
            newparts.append(newpart)

        if hasattr(model, 'op'):
            model = model.op(*newparts)
        elif hasattr(model, 'rmf') and hasattr(model, 'arf'):
            model = sherpa.astro.instrument.RSPModelPHA(rmf=model.rmf,
                                                        model=newparts[0],
                                                        arf=model.arf,
                                                        pha=model.pha)
        elif hasattr(model, 'rmf'):
            model = sherpa.astro.instrument.RMFModelPHA(rmf=model.rmf,
                                                        model=newparts[0],
                                                        pha=model.pha)
        elif hasattr(model, 'arf'):
            model = sherpa.astro.instrument.ARFModelPHA(model=newparts[0],
                                                        arf=model.arf,
                                                        pha=model.pha)
        else:
            raise ValueError(
                "Unexpected composite model {0} (not operator, ARF or RMF)".format(repr(model)))
    else:
        if hasattr(model, 'val'):
            model = model.val
        else:
            try:
                model_type, model_name_ID = model.name.split('.')
            except ValueError:
                raise ValueError(
                    'Model name "{0}" must be in format <model_type>.<name>'.format(model.name))

            model_name = re.sub(ID_STR, str(model_id), model_name_ID)
            if ID_STR in model_name_ID:
                try:
                    model = getattr(getattr(sherpa.astro.ui, model_type),
                                    model_name)
                except AttributeError:
                    # Must be a user model, so use add_model to put a
                    # modelwrapper function into namespace
                    sherpa.astro.ui.add_model(model.__class__)
                    model = sherpa.ui.create_model_component(
                        model_type, model_name)

            model_name_no_ID = re.sub(ID_STR, "", model_name_ID)
            model_comps[model_name_no_ID] = dict(model=model,
                                                 model_name=model_name)

    return model, model_comps


@public
def load_wrapper(load_func):
    """Override a native Sherpa data loading function."""

    def _load(self, *args, **kwargs):

        if len(args) == 1:
            id, arg = None, args[0]
            args = []
        if len(args) > 1:
            args = args[1:]

        if id is not None:
            if self._default_instance:
                self._load_func(load_func, id, arg, *args, **kwargs)
                return
            else:
                raise AttributeError(load_error_msg(id))

        # File Stacks. If the file argument is a stack file, expand
        # the file and call this function for each file in the stack.
        try:
            files = stk.build(arg)
            for file in files:
                self._load_func(load_func, file, *args, **kwargs)
        except:
            self._load_func(load_func, arg, *args, **kwargs)

    _load.__name__ = load_func.__name__
    _load.__doc__ = load_func.__doc__
    return _load


@public
def simple_wrapper(func):
    def wrapfunc(self, *args, **kwargs):
        """Apply an arbitrary Sherpa function to each of the datasets.

        Returns
        -------
        out : list
           The return value of each call, which may be ``None``.
        """
        datasets = self.filter_datasets()
        logger.info('Running {0} with args={1} and kwargs={2} for ids={3}'.format(
            func.__name__, args, kwargs, [x['id'] for x in datasets]))
        return [func(x['id'], *args, **kwargs) for x in datasets]

    wrapfunc.__name__ = func.__name__
    wrapfunc.__doc__ = func.__doc__
    return wrapfunc


@public
def fit_wrapper(func):
    def _fit(self, *args, **kwargs):
        """Fit or error analysis for all the datasets in the stack.

        Parameters
        ----------
        args
           Any additional arguments to be passed to the call.
        kwargs
           Any keyword arguments to be passed to the call.
        """
        ids = tuple(x['id'] for x in self.filter_datasets())
        func(*(ids + args), **kwargs)

    _fit.__name__ = func.__name__
    return _fit


@public
def plot_wrapper(func):
    def plot(self, *args, **kwargs):
        """Call the Sherpa plot ``func`` for each dataset.

        Parameters
        ----------
        func
           The Sherpa plot function
        args
           The arguments to the function.
        kwargs
           The keyword arguments to the function.
        """
        backend.initialize_backend()
        for dataset in self.filter_datasets():
            backend.initialize_plot(dataset, self.ids)
            func(dataset['id'], *args, **kwargs)

    plot.__name__ = func.__name__
    return plot


@public
def set_template_id(newid):
    """Change the template used for model components.

    When source models are created for a stack of data sets,
    the identifier for the stack can be added to the model
    name by including the template id value in its name.
    The default value is ``__ID``, but it can be changed by
    this routine.

    Parameters
    ----------
    newid : string
       The new template value.

    Examples
    --------

    This creates one ``const1d`` component - shared between all
    data sets in the datastack - and a separate ``polynom1d``
    instance for each data set. The name of these instances is
    ``poly`` and ends in the data set identifier.

    >>> set_source([], 'const1d.const * polynom1d.poly__ID')

    The polynomial-model components will now be called
    ``poly_`` and then the data set identifier.

    >>> set_template_id('ID')
    >>> set_source([], 'const1d.const * polynom1d.poly_ID')

    """
    global ID_STR
    ID_STR = newid


@public
def load_error_msg(id_):
    """The error message when an invalid id is given."""
    return "When called from a datastack instance, an ID cannot be " + \
        " provided to a load function ({0})".format(id_)
