# Copyright (c) 2010,2014,2015 Smithsonian Astrophysical Observatory
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Smithsonian Astrophysical Observatory nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""Support multiple data sets when using Sherpa commands.

The ``sherpa.astro.datastack`` module supports manipulating a stack
of related datasets and simultaneously fitting a model to them. It
provides stack-enabled (i.e. vectorized) versions of the key Sherpa
commands used to load data, set source models, get and set parameters,
fit, and plot. It is intended for use when simultaneously fitting a
common source model to multiple data sets.

Example: PHA data
-----------------

In the following example, ``src.lis`` is an ASCII text file
with the name of a PHA file to load on each line::

    from sherpa.astro import datastack
    from sherpa.astro import ui
    datastack.load_pha("@src.lis")

At this point the PHA files are loaded into data sets ``1`` to ``n``,
where ``n`` is the number of lines in the ``src.lis``. Any ancillary
files - such as background, ARF, and RMF - will be loaded in as if
the files were loaded separately.

The loaded data sets can be shown using::

    datastack.show_stack()

The module uses the special identifier ``[]`` to indicate all
members of a stack, so::

    datastack.set_source([], ui.xsphabs.gal * ui.xspowerlaw.pl)

will set each file to have the *same* model (in this case an
absorbed power law). Adding the suffix ``__ID`` to a component
name will create a separate component for each data set; so::

    src = ui.const1d.c__ID * ui.xsphabs.gal * ui.xspowerlaw.pl
    datastack.set_source([], src)
    ui.freeze(pl.norm)

will have a common absorbing component (``gal``) and power law
model (``pl``), but each data set has a separate constant term
labelled ``c`` followed by the data set identifier (e.g.
``c1`` and ``c2``). Since the normalization of the power law
component has been frozen the constant term represents the
normalization of each component (i.e. the model shape is
assumed constant, but its amplitude is not). These expressions
can be viewed using the command::

    ui.show_source()

The ``integrate`` flag of the constant model component should
be turned off (so that the component acts as a scalar term rather
than including the bin-width). The ``datastack`` module does not
provide a simple way to do this, so the setting of each component
has to be changed individually::

    for did in datastack.get_stack_ids():
        mdl = ui.get_model_component('c{}'.format(did))
        mdl.integrate = False

The ``datastack`` module provides versions of the ``sherpa.astro.ui``
module which accept ``[]``, so::

    datastack.subtract([])

will subtract the background from each data set. Some commands are
the same - so either of the following will filter the data::

    ui.notice(0.5, 7)
    datastack.notice(0.5, 7)

The data and model for each data set can be viewed with::

    datastack.plot_fit([])

and the model fit to all the data sets as normal::

    ui.fit()

Loading data
------------

Multiple inputs (referred to as a stack here) can be specified
either from an ASCII file - one file name per line - by placing
the ``@`` character before the file name, or as a comma-separated
list of names:

- ``load_data``
- ``load_ascii``
- ``load_pha``
- ``load_bkg``

The ``load_arrays`` function is slightly different, in that it
accepts a list of array arguments, one for each dataset.

Examples include::

    datastack.load_data("@srcs.lis")
    datastack.load_pha("obs1.pha,obs2.pha,obs3.pha")
    datastack.load_arrays([[x1, y1], [x2, y2]])

Identifying a stack
-------------------

When reading in a stack of data, the individual data sets
are numbered sequentially. These identifiers can be used to
select individual data sets using functions from the
``sherpa.astro.ui`` module. The functions from the ``datastack``
module work with a datastack identifier, which can be:

- ``[]``
- an iterable sequence of data set identifiers
- a datastack instance reference
- a subset of a datastack instance

So::

    datastack.plot_data([])
    datastack.plot_data([1,3])

plots all the data sets, and then just the first and third entries.
The following repeats the example, using a ``DataStack`` object::

    ds = datastack.DataStack()
    ds.load_data('@src.lis')
    datastack.plot_data(ds)
    datastack.plot_data(ds[1,3])

Note that when accessing a subset of a DataStack object - e.g.
``ds[1,3]`` - the numbers match the data set identifiers (and so
start at ``1``, not ``0``).

Setting a model for a stack
---------------------------

The functions

- ``set_source``
- ``set_model``
- ``set_bkg_model``
- ``set_full_model``
- ``set_bkg_full_model``

can be used to set the source expression for a datastack. This
expression can include models with components that are shared between
data sets and models which have a component per data set. The
later case are created by using the identifier ``__ID`` in the
name of the component. The following call will fit the sum of a
polynomial and gaussian model to the data, with the same parameters
used for each data set (the model components are called ``bgnd``
and ``src`` respectively)::

    datastack.set_source([], ui.polynom1d.bgnd + ui.gauss1d.src)

whereas::

    datastack.set_source([], ui.polynom1d.bgnd__ID + ui.gauss1d.src)

fits a single gaussian model (``src``) to all data sets, but allows the
polynomial to vary between datasets (with names ``bgnd1``, ``bgnd2``,
...).

Utility functions for data stacks
---------------------------------

The following functions are provided:

- ``get_stack_ids``, which returns the data set identifiers that are
  included in the data stack,
- ``show_stack``, which prints the data sets available in the data stack,
  together with some basic metadata,
- ``clear_models``, to remove all the model expressions from a data stack,
- and ``clear_stack``, to remove all the data sets that are part of the
  data stack.

Information about data sets which match a particular query are provided
by the ``query_by_header_keyword``, ``query_by_obsid``, and ``query``
functions.

"""

import sys
import types
import re
import ConfigParser
import numpy
import sherpa
import sherpa.astro.ui as ui
import logging
import stk


# Configure logging for the module
def _config_logger(name, level, stream):
    logger = logging.getLogger(name)
    for hdlr in logger.handlers:
        logger.removeHandler(hdlr)
    logger.setLevel(level)
    logger.propagate = False

    fmt = logging.Formatter('Datastack: %(message)s', None)
    hdlr = logging.StreamHandler(stream)
    # hdlr.setLevel(level)
    hdlr.setFormatter(fmt)
    logger.addHandler(hdlr)

    return logger

logger = _config_logger(__name__, level=logging.WARNING, stream=sys.stdout)

def set_stack_verbosity(level):
    logger.setLevel(level)

def set_stack_verbose(verbose=True):
    """Should stack functions print informational messages?

    Parameters
    ----------
    verbose : bool, opt
       If ``True`` then display messages at the logging level
       of ``INFO`` or above. If ``False`` only display
       ``WARNING`` messages or above.
    """
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

# Get plot package
_cp = ConfigParser.ConfigParser()
_cp.read(sherpa.get_config())
_plot_pkg =  _cp.get('options', 'plot_pkg')
if _plot_pkg == 'pylab':
    import matplotlib.pyplot as plt
elif _plot_pkg == 'chips':
    import pychips
    from sherpa.plot import chips_backend
else:
    raise ValueError('Unknown plot package {0}'.format(_plot_pkg))

# Global list of dataset ids in use
_all_dataset_ids = {}


id_str = '__ID'


def set_template_id(newid):
    global id_str
    id_str = newid

def create_stack_model(model, id_):
    model_comps = {}
    def _get_new_model(model, level=0):
        if hasattr(model, 'parts'):
            # Recursively descend through model and create new parts (as needed)
            # corresponding to the stacked model components.
            newparts = []
            for part in model.parts:
                newparts.append(_get_new_model(part, level+1))

            if hasattr(model, 'op'):
                return model.op(*newparts)
            elif isinstance(model, sherpa.astro.instrument.RSPModelPHA):
                return sherpa.astro.instrument.RSPModelPHA(rmf=model.rmf, model=newparts[0],
                                                        arf=model.arf, pha=model.pha)
            elif isinstance(model, sherpa.astro.instrument.RMFModelPHA):
                return sherpa.astro.instrument.RSPModelPHA(rmf=model.rmf, model=newparts[0],
                                                        arf=model.arf, pha=model.pha)
            elif isinstance(model, sherpa.astro.instrument.ARFModelPHA):
                return sherpa.astro.instrument.ARFModelPHA(rmf=model.rmf, model=newparts[0],
                                                        arf=model.arf, pha=model.pha)
            else:
                raise ValueError("Unexpected composite model {0} (not operator, ARF or RMF)".format(repr(model)))
        else:
            if isinstance(model, sherpa.models.model.ArithmeticConstantModel):
                return model.val

            try:
                model_type, model_name_ID = model.name.split('.')
            except ValueError:
                raise ValueError('Model name "{0}" must be in format <model_type>.<name>'.format(model.name))

            model_name = re.sub(id_str, str(id_), model_name_ID)
            if id_str in model_name_ID:
                try:
                    model = getattr(getattr(sherpa.astro.ui, model_type), model_name)
                except AttributeError:
                    # Must be a user model, so use add_model to put a modelwrapper function into namespace
                    sherpa.astro.ui.add_model(type(model))
                    model = eval('{0}.{1}'.format(model_type, model_name))

            model_name_no_ID = re.sub(id_str, "", model_name_ID)
            model_comps[model_name_no_ID] = dict(model=model,
                                                 model_name=model_name)

            return model

    return _get_new_model(model), model_comps

class DataStack(object):
    """
    Manipulate a stack of data in Sherpa.
    """
    def __init__(self):
        self.getitem_ids = None
        self.datasets = []
        self.dataset_ids = {}           # Access datasets by ID

    def __getitem__(self, item):
        """Overload datastack getitem ds[item(s)] to set self.filter_ids to a tuple
        corresponding to the specified items.
        """
        try:
            ids = (item + '',)
        except TypeError:
            try:
                ids = tuple(item)
            except TypeError:
                ids = (item, )
        self.getitem_ids = ids
        return self

    def __del__(self):
        for dataid in self.dataset_ids:
            try:
                del _all_dataset_ids[dataid]
            except:
                pass

    def clear_models(self):
        """Clear all model components in the stack.

        See Also
        --------
        clear_stack, clean, show_stack
        """
        for dataset in self.datasets:
            dataset['model_comps'] = {}

    def clear_stack(self):
        """Clear all datasets in the stack.

        See Also
        --------
        clear_models, clean, show_stack
        """
        for dataid in self.dataset_ids:
            del _all_dataset_ids[dataid]
        self.__init__()

    def show_stack(self):
        """Display basic information about the contents of the data stack.

        A brief summary of each data set in the stack is printed to the
        standard output channel. The information displayed depends on the
        type of each data set.
        """
        for dataset in self.filter_datasets():
            obsid = "N/A"
            time = "N/A"
            if hasattr(dataset['data'], 'header') and "OBS_ID" in dataset['data'].header.keys():
                obsid = dataset['data'].header['OBS_ID']
            if hasattr(dataset['data'], 'header') and "MJD_OBS" in dataset['data'].header.keys():
                time = dataset['data'].header['MJD_OBS']
            print('{0}: {1} {2}: {3} {4}: {5}'.format(dataset['id'], dataset['data'].name, 'OBS_ID', obsid, "MJD_OBS", time))

    def get_stack_ids(self):
        """Return the data set identifiers in the stack.

        Returns
        -------
        ids : array of int or str
           The data set identifiers.

        """
        return self.ids

    @property
    def ids(self):
        """List of ids corresponding to stack datasets.

        Returns
        -------
        ids : array of int or str
           The data set identifiers.

        """
        return [x['id'] for x in self.datasets]

    def _get_dataid(self):
        if self.getitem_ids:
            dataid = self.getitem_ids[0]
            self.getitem_ids = None
        else:
            dataid = 1
            while dataid in _all_dataset_ids:
                dataid += 1

        if dataid in self.dataset_ids:
            raise ValueError('Data ID = {0} is already in the DataStack'.format(dataid))

        return dataid

    def _add_dataset(self, dataid):
        dataset = dict(id=dataid, args=[], model_comps={}, data=ui.get_data(dataid))
        _all_dataset_ids[dataid] = dataset
        self.dataset_ids[dataid] = dataset
        self.datasets.append(dataset)

    def _load_func(self, func, *args, **kwargs):
        dataid = self._get_dataid()

        logger.info('Loading dataset id %s' % dataid)
        func(dataid, *args, **kwargs)

        self._add_dataset(dataid)

    # NOTE: this doc string is over-ridden by the contents of the
    #       sherpa.astro.ui version
    def load_arrays(self, *args, **kwargs):
        """Load multiple data arrays.

        This extends ``sherpa.astro.ui.load_arrays`` to load multiple
        data sets with one call. The array arguments are supplied as an
        array of arrays, one for each data set. The other arguments are
        as supported by ``sherpa.astro.ui.load_arrays``.

        Examples
        --------

        Create three data sets, using the array pairs ``x1, y1``,
        ``x2, y2``, and ``x3, y3``.

        >>> load_arrays([[x1, y1], [x2, y2], [x3, y3]])

        """

        if len(args) == 0:
            raise AttributeError("load_arrays takes at least one argument (got none).")

        if not hasattr(args[0], '__iter__'):
            id, args = args[0], args[1:]
        else:
            id = None

        if id is not None:
            if self is DATASTACK:
                ui.load_arrays(id, *args, **kwargs)
                return
            else:
                raise AttributeError("When called from a datastack instance, an ID cannot be provided to a load function ("+id+")")

        # Array Stack.
        for arrays in args[0]:
            dataid = self._get_dataid()
            logger.info('Loading dataset id %s' % dataid)
            ui.load_arrays(dataid, *arrays)
            self._add_dataset(dataid)


    def load_pha(self, id, arg=None, use_errors=False):
        if arg is None:
            id, arg = arg, id

        if id is not None:
            if self is DATASTACK:
                ui.load_pha(id, arg, use_errors)
                return
            else:
                raise AttributeError("When called from a datastack instance, an ID cannot be provided to a load function ("+id+")")

        # File Stacks. If the file argument is a stack file, expand the file and call this function for each file
        #   in the stack.
        try:
            files = stk.build(arg)
            for file in files:
                self._load_func(ui.load_pha, file, use_errors)
        except:
            self._load_func(ui.load_pha, arg, use_errors)


    def _load_func_factory(load_func):
        """Override a native Sherpa data loading function."""
        def _load(self, *args, **kwargs):

            if len(args)==1:
                id, arg = None, args[0]
                args=[]
            if len(args)>1:
                args = args[1:]

            if id is not None:
                if self is DATASTACK:
                    self._load_func(load_func, id, arg, *args, **kwargs)
                    return
                else:
                    raise AttributeError("When called from a datastack instance, an ID cannot be provided to a load function ("+id+")")

            # File Stacks. If the file argument is a stack file, expand the file and call this function for each file
            #   in the stack.
            try:
                files = stk.build(arg)
                for file in files:
                    self._load_func(load_func, file, *args, **kwargs)
            except:
                self._load_func(load_func, arg, *args, **kwargs)

        # def _load(self, *args, **kwargs):
        #     """Load a dataset and add to the datasets for stacked analysis.
        #     """
        #     dataid = self._get_dataid()
        #
        #     logger.info('Loading dataset id %s' % dataid)
        #     out = load_func(dataid, *args, **kwargs)
        #
        #     #dataset = dict(id=dataid, args=args, model_comps={}, data=ui.get_data(dataid))
        #     #dataset.update(kwargs)  # no sherpa load func 'args' keyword so no conflict
        #     self._add_dataset(dataid)
        #
        #     return out

        _load.__name__ = load_func.__name__
        _load.__doc__ = load_func.__doc__
        return _load

#    load_arrays = _load_func_factory(ui.load_arrays)
    load_ascii = _load_func_factory(ui.load_ascii)
    load_data = _load_func_factory(ui.load_data)
#    load_image = _load_func_factory(ui.load_image)
#    load_pha = _load_func_factory(ui.load_pha)
    load_bkg = _load_func_factory(ui.load_bkg)

    def _set_model_factory(func):
        def wrapfunc(self, model):
            """Run a model-setting function for each of the datasets."""
            datasets = self.filter_datasets()
            try:
                # if model is passed as a string
                model = eval(model, globals(), ui.__dict__)
            except TypeError:
                pass
            except Exception, exc:
                raise type(exc)('Error converting model "{0}" '
                                 'to a sherpa model object: {1}'.format(model, exc))
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

    set_source = _set_model_factory(ui.set_source)
    set_model = _set_model_factory(ui.set_model)
    set_bkg_model = _set_model_factory(ui.set_bkg_model)
    set_full_model = _set_model_factory(ui.set_full_model)
    set_bkg_full_model = _set_model_factory(ui.set_bkg_full_model)

    def filter_datasets(self):
        """Return filtered list of datasets as specified in the __getitem__
        argument (via self.getitem_ids which gets set in __getitem__).
        """
        if self.getitem_ids is None:
            return self.datasets

        filter_ids = self.getitem_ids
        self.getitem_ids = None

        try:
            return [self.dataset_ids[x] for x in filter_ids]
        except KeyError:
            raise ValueError('IDs = {0} not contained in dataset IDs = {1}'.format(filter_ids, self.ids))

    def _sherpa_cmd_factory(func):
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

    subtract = _sherpa_cmd_factory(ui.subtract)
    notice = _sherpa_cmd_factory(ui.notice_id)
    ignore = _sherpa_cmd_factory(ui.ignore_id)
    get_arf = _sherpa_cmd_factory(ui.get_arf)
    get_rmf = _sherpa_cmd_factory(ui.get_rmf)
    get_response = _sherpa_cmd_factory(ui.get_response)
    get_bkg_arf = _sherpa_cmd_factory(ui.get_bkg_arf)
    get_bkg_rmf = _sherpa_cmd_factory(ui.get_bkg_rmf)
    get_bkg = _sherpa_cmd_factory(ui.get_bkg)
    get_source = _sherpa_cmd_factory(ui.get_source)
    get_model = _sherpa_cmd_factory(ui.get_model)
    get_bkg_model = _sherpa_cmd_factory(ui.get_bkg_model)
    try:
        get_bkg_scale = _sherpa_cmd_factory(ui.get_bkg_scale)
    except AttributeError:
        pass  # not available for CIAO < 4.3
    group_adapt = _sherpa_cmd_factory(ui.group_adapt)
    group_adapt_snr = _sherpa_cmd_factory(ui.group_adapt_snr)
    group_bins = _sherpa_cmd_factory(ui.group_bins)
    group_counts = _sherpa_cmd_factory(ui.group_counts)
    group_snr = _sherpa_cmd_factory(ui.group_snr)
    group_width = _sherpa_cmd_factory(ui.group_width)
    load_arf = _sherpa_cmd_factory(ui.load_arf)
    load_rmf = _sherpa_cmd_factory(ui.load_rmf)
    load_bkg_arf = _sherpa_cmd_factory(ui.load_bkg_arf)
    load_bkg_rmf = _sherpa_cmd_factory(ui.load_bkg_rmf)
    load_filter = _sherpa_cmd_factory(ui.load_filter)
    load_grouping = _sherpa_cmd_factory(ui.load_grouping)

    def _sherpa_par(self, func, par, msg, *args, **kwargs):
        """Apply ``func(*args)`` to all model component or model component parameters named ``mcpar``.

        See thaw(), freeze(), set_par() and get_par() for examples.

        Parameters
        ----------
        func
           The Sherpa function to call. It should take a full parameter name
           specification and optional args, e.g. ``set_par`` since it can
           be called as ``set_par('mekal_7.kt', 2.0)``.
        par : str
           The name of the parameter of model component.
        msg : str
           The format string to inficate action.
        *args
           Any remaining arguments to send to the function.
        *args
           Any keyword arguments to send to the function.

        Returns
        -------
        out : list
           A list of the return value from each function call, which can
           be ``None``.
        """

        vals = par.split('.')
        name = vals[0]
        parname = (vals[1] if len(vals) > 1 else None)
        if len(vals) > 2:
            raise ValueError('Invalid parameter name specification "%s"' % par)

        retvals = []                       # return values
        processed = set()
        for dataset in self.filter_datasets():
            model_comps = dataset['model_comps']
            if name in model_comps:
                model_name = model_comps[name]['model_name']
                fullparname = '{0}.{1}'.format(model_name, parname) if parname else model_name
                if fullparname not in processed:
                    if msg is not None:
                        logger.info(msg % fullparname)
                    retvals.append(func(fullparname, *args, **kwargs))
                    processed.add(fullparname)

        return retvals

    # NOTE: this doc string is over-ridden by the contents of the
    #       sherpa.astro.ui version
    def thaw(self, *pars):
        """Apply thaw command to specified parameters for each dataset.

        Parameters
        ----------
        pars
           The names of the parameters or model components to thaw.
        """
        for par in pars:
            self._sherpa_par(ui.thaw, par, 'Thawing %s')

    # NOTE: this doc string is over-ridden by the contents of the
    #       sherpa.astro.ui version
    def freeze(self, *pars):
        """Apply freeze command to specified parameters for each dataset.

        Parameters
        ----------
        pars
           The names of the parameters or model components to freeze.
        """
        for par in pars:
            self._sherpa_par(ui.freeze, par, 'Freezing %s')

    def _get_parname_attr_pars(self, par, msg):
        parts = par.split('.')
        if len(parts) == 1:
            raise ValueError('par="%s" must be in the form "name.par" or "name.par.attr"' % par)
        parname = '.'.join(parts[:-1])
        attr = parts[-1]
        return parname, attr, self._sherpa_par(eval, parname, msg % attr)

    # NOTE: this doc string is over-ridden by the contents of the
    #       sherpa.astro.ui version
    def get_par(self, par):
        """Get parameter object for each dataset.

        Parameters
        ----------
        par : str
           The name of the parameters, incuding the model component name.
        """
        parname, attr, pars = self._get_parname_attr_pars(par, 'Getting %%s.%s')
        return numpy.array([getattr(x, attr) for x in pars])

    def set_par(self, par, val):
        """Set parameter attribute value for each dataset.

        Parameters
        ----------
        par : str
           The name of the parameter, including the attribute name.
        val
           The parameter attribute value.

        """
        parname, attr, pars = self._get_parname_attr_pars(par, 'Setting %%%%s.%%s = %s' % val)
        for x in pars:
            setattr(x, attr, val)

    def link(self, par):
        datasets = self.filter_datasets()
        name, parname = par.split('.')

        fullparname0 = '{0}.{1}'.format(datasets[0]['model_comps'][name]['model_name'], parname)
        for dataset in datasets[1:]:
            fullparname = '{0}.{1}'.format(dataset['model_comps'][name]['model_name'], parname)
            if fullparname != fullparname0:
                logger.info('Linking {0} => {1}'.format(fullparname, fullparname0))
                ui.link(fullparname, fullparname0)

    def unlink(self, par):
        self._sherpa_par(ui.unlink, par, 'Unlinking %s')

    def _sherpa_fit_func(func):
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
        _fit.__doc__ = func.__doc__

        return _fit

    fit_bkg = _sherpa_fit_func(ui.fit_bkg)
    fit = _sherpa_fit_func(ui.fit)
    conf = _sherpa_fit_func(ui.conf)


    def _print_window(self, *args, **kwargs):
        """Save figure for each dataset.

        Parameters
        ----------
        args
           The parameters sent to the plotting back end to create a hardcopy
           version of the plot. The first argument - if given - is taken to
           be the filename, and any ``#`` character in it will be replaced
           by the data set identifier (so that multiple files will be
           created).
        kwargs
           Any keyword arguments to pass to the plotting back end.
        """
        orig_args = args
        for dataset in self.filter_datasets():
            args = orig_args
            if len(args) > 0:
                filename = re.sub(r'#', str(dataset['id']), args[0])
                args = tuple([filename]) + args[1:]

            if _plot_pkg == 'chips':
                pychips.set_current_window(dataset['id'])
                func = pychips.print_window
            elif _plot_pkg == 'pylab':
                plt.figure(self.ids.index(dataset['id']) + 1)
                func = plt.savefig
            else:
                raise ValueError('Unknown plot package')

            func(*args, **kwargs)

    savefig = _print_window              # matplotlib alias for print_window

    def _sherpa_plot_func(func):
        def _sherpa_plot(self, *args, **kwargs):
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
            # FIXME We need to make sure ChIPS is initialized.
            # As a hack, we just call begin() and end() on the chips backend
            # To make sure the backend is initialized before we start creating windows.
            if _plot_pkg == 'chips':
                try:
                    chips_backend.begin()
                finally:
                    chips_backend.end()
            for dataset in self.filter_datasets():
                if _plot_pkg == 'chips':
                    try:
                        pychips.add_window(['id', dataset['id']])
                    except RuntimeError:
                        pass  # already exists
                    # window_id = pychips.ChipsId()
                    # window_id.window = dataset['id']
                    pychips.current_window(str(dataset['id']))
                elif _plot_pkg == 'pylab':
                    plt.figure(self.ids.index(dataset['id']) + 1)
                else:
                    raise ValueError('Unknown plot package')

                func(dataset['id'], *args, **kwargs)
        return _sherpa_plot

    # log_scale = _sherpa_plot_func(pychips.log_scale)
    # linear_scale = _sherpa_plot_func(pychips.linear_scale)

    plot_arf = _sherpa_plot_func(ui.plot_arf)
    plot_bkg_fit = _sherpa_plot_func(ui.plot_bkg_fit)
    plot_bkg_ratio = _sherpa_plot_func(ui.plot_bkg_ratio)
    plot_chisqr = _sherpa_plot_func(ui.plot_chisqr)
    plot_fit_delchi = _sherpa_plot_func(ui.plot_fit_delchi)
    plot_psf = _sherpa_plot_func(ui.plot_psf)
    plot_bkg = _sherpa_plot_func(ui.plot_bkg)
    plot_bkg_fit_delchi = _sherpa_plot_func(ui.plot_bkg_fit_delchi)
    plot_bkg_resid = _sherpa_plot_func(ui.plot_bkg_resid)
    plot_data = _sherpa_plot_func(ui.plot_data)
    plot_fit_resid = _sherpa_plot_func(ui.plot_fit_resid)
    plot_ratio = _sherpa_plot_func(ui.plot_ratio)
    plot_bkg_chisqr = _sherpa_plot_func(ui.plot_bkg_chisqr)
    plot_bkg_fit_resid = _sherpa_plot_func(ui.plot_bkg_fit_resid)
    plot_bkg_source = _sherpa_plot_func(ui.plot_bkg_source)
    plot_delchi = _sherpa_plot_func(ui.plot_delchi)
    plot_model = _sherpa_plot_func(ui.plot_model)
    plot_resid = _sherpa_plot_func(ui.plot_resid)
    plot_bkg_delchi = _sherpa_plot_func(ui.plot_bkg_delchi)
    plot_bkg_model = _sherpa_plot_func(ui.plot_bkg_model)
    plot_bkg_source = _sherpa_plot_func(ui.plot_bkg_source)
    plot_fit = _sherpa_plot_func(ui.plot_fit)
    plot_order = _sherpa_plot_func(ui.plot_order)
    plot_source = _sherpa_plot_func(ui.plot_source)

    def _matplotlib_func(func, axis_cmd=False):
        def _matplotlib(self, *args, **kwargs):
            """Call the matplotlib plot ``func`` for each dataset.

            Parameters
            ----------
            func
               The matplotlib plot function
            args
               The arguments to the function.
            kwargs
               The keyword arguments to the function.
            """
            orig_args = args
            for dataset in self.filter_datasets():
                args = orig_args
                if _plot_pkg != 'pylab':
                    raise ValueError('Plot package must be pylab')

                if len(args) > 0:
                    try:
                        arg0 = re.sub('#', str(dataset['id']), args[0])
                        args = tuple([arg0]) + args[1:]
                    except TypeError:
                        pass

                plt.figure(self.ids.index(dataset['id']) + 1)
                if axis_cmd:
                    ax = plt.gca()
                    getattr(ax, func)(*args, **kwargs)
                else:
                    func(*args, **kwargs)

        return _matplotlib

    if _plot_pkg == 'pylab':
        plot_savefig = _matplotlib_func(plt.savefig)
        plot_xlabel = _matplotlib_func(plt.xlabel)
        plot_ylabel = _matplotlib_func(plt.ylabel)
        plot_title = _matplotlib_func(plt.title)
        plot_xlim = _matplotlib_func(plt.xlim)
        plot_ylim = _matplotlib_func(plt.ylim)
        plot_set_xscale = _matplotlib_func('set_xscale', axis_cmd=True)
        plot_set_yscale = _matplotlib_func('set_yscale', axis_cmd=True)

    def query(self, func):
        """Return the data sets identified by a function.

        Parameters
        ----------
        func
           A function which accepts a Sherpa Data object and
           returns ``True`` if the data set matches.

        Returns
        -------
        matches : list of int or str
           A list of the data set identifiers that match. This
           list may be empty.

        See Also
        --------
        query_by_header_keyword, query_by_obsid

        Examples
        --------

        This query function selects those data sets which have
        been grouped, or had their background subtracted, or both.

        >>> myquery = lambda d: d.subtracted or d.grouped
        >>> query(myquery)
        []

        """

        output = []
        for dataset in self.filter_datasets():
            id = dataset['id']
            if func(ui.get_data(id)):
                output.append(id)
        return output

    def query_by_header_keyword(self, keyword, value):
        """Return the data sets which contain the given keyword value.

        Parameters
        ----------
        keyword : str
           The keyword to search for (it is looked for in the ``header``
           attribute of each data set, using a case-sensitive search).
        value
           The value of the keyword.

        Returns
        -------
        matches : list of int or str
           A list of the data set identifiers that match the query. This
           list may be empty.

        See Also
        --------
        query, query_by_obsid

        Examples
        --------

        >>> query_by_header_keyword('TELESCOP', 'CHANDRA')
        [1, 2]

        """

        def func(dataset):
            if hasattr(dataset, 'header'):
                if keyword in dataset.header.keys() and dataset.header[keyword] == str(value):
                    return True
            return False
        return self.query(func)

    def query_by_obsid(self, value):
        """Return the data sets which match the OBS_ID value.

        Parameters
        ----------
        value : int or str
           The value to match against the ``OBS_ID`` keyword in
           the header of each data set.

        Returns
        -------
        matches : list of int or str
           A list of the data set identifiers that match. This
           list may be empty.

        See Also
        --------
        query, query_by_header_keyword

        Examples
        --------

        >>> query_by_obsid(7867)
        [3]

        """

        return self.query_by_header_keyword('OBS_ID', value)


_always_wrapped = ('load_pha', 'load_arrays', 'load_ascii', 'load_data', 'load_bkg')

# Use this and subsequent loop to wrap every function in sherpa.astro.ui with a datastack version
def _sherpa_ui_wrap(func):
    def wrap(*args, **kwargs):
        wrapfunc = func
        if args:
            if isinstance(args[0], DataStack):
                datastack, args = args[0], args[1:]
            # If the first argument is a list and it's either empty or made of non-iterables, then it's a datastack definition.
            # If the list contains iterable it must be arrays for load_arrays.
            elif isinstance(args[0], list) and not (len(args[0])>0 and hasattr(args[0][0],'__iter__')):
                datastack, args = (DATASTACK[args[0]] if args[0] else DATASTACK), args[1:]
            else:
                if func.__name__ in _always_wrapped:
                # some (all?) load_* functions must always be wrapped for file stack syntax check
                # and for ensuring dataset id consistency.
                    datastack = DATASTACK
                else:
                    return func(*args, **kwargs) # No stack specifier so use native sherpa func

            try:
                wrapfunc = getattr(datastack, func.__name__)
            except AttributeError:
                raise AttributeError(
                    '{0} is not a stack-enabled function.'.format(func.__name__))

        return wrapfunc(*args, **kwargs)

    wrap.__name__ = func.__name__
    wrap.__doc__ = func.__doc__
    return wrap

def _datastack_wrap(func):
    def wrap(*args, **kwargs):
        if not args:
            args = ([],) + args

        if isinstance(args[0], DataStack):
            datastack, args = args[0], args[1:]
        elif isinstance(args[0], list) and not (len(args[0])>0 and hasattr(args[0][0],'__iter__')):
            datastack, args = (DATASTACK[args[0]] if args[0] else DATASTACK), args[1:]
        else:
            datastack = DATASTACK
            # raise TypeError('First argument to {0} must be a list or datastack')

        return getattr(datastack, func.__name__)(*args, **kwargs)

    wrap.__name__ = func.__name__
    wrap.__doc__ = func.__doc__
    return wrap

DATASTACK = DataStack()  # Default datastack

# Wrap all sherpa UI funcs and a few DataStack methods for command-line interface.
_module = sys.modules[__name__]
for attr in dir(ui):
    func = getattr(ui, attr)
    if type(func) == types.FunctionType:
        setattr(_module, attr, _sherpa_ui_wrap(func))

for funcname in ['clear_stack', 'show_stack', 'get_stack_ids', 'query', 'query_by_header_keyword', 'query_by_obsid']:
    setattr(_module, funcname, _datastack_wrap(getattr(DataStack, funcname)))

def clean():
    """Remove the models and data from the data stack and Sherpa.

    This function clears out the models and data set up in the data
    stack and in the Sherpa session.

    See Also
    --------
    clear_models, clear_stack
    sherpa.astro.ui.clean
    """
    DATASTACK.clear_models()
    DATASTACK.clear_stack()
    ui.clean()
    logger.warning("clean() will invalidate any existing DataStack instances by removing all the datasets from the " +
                   "Sherpa session")
