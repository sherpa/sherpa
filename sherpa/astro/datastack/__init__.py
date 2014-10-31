"""
Manipulate a stack of data in Sherpa.

The methods in the DataStack class provide a way to automatically apply
familiar Sherpa commands such as `set_par`_ or `freeze`_ or `plot_fit`_
to a stack of datasets.  This simplifies simultaneous fitting of
multiple datasets.

:Copyright: Smithsonian Astrophysical Observatory (2014)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
:Author: Omar Laurino (olaurino@head.cfa.harvard.edu)
"""
## Copyright (c) 2010, Smithsonian Astrophysical Observatory
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * Neither the name of the Smithsonian Astrophysical Observatory nor the
##       names of its contributors may be used to endorse or promote products
##       derived from this software without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS  
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
     """Configure whether stack functions print informational messages.
     :param verbose: print messages if True (default=True)
     :returns: None
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
        :returns: None
        """
        for dataset in self.datasets:
            dataset['model_comps'] = {}

    def clear_stack(self):
        """Clear all datasets in the stack.
        :returns: None
        """
        for dataid in self.dataset_ids:
            del _all_dataset_ids[dataid]
        self.__init__()

    def show_stack(self):
        """Show the data id and file name (where meaningful) for selected
        datasets in stack.
        :returns: None
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
        """Get the ids for all datasets in stack
        :returns: list of ids
        """
        return self.ids

    @property
    def ids(self):
        """List of ids corresponding to stack datasets"""
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

    def load_arrays(self, *args, **kwargs):
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
            """
            Run a model-setting function for each of the datasets.  
            :rtype: None
            """
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
            """
            Apply an arbitrary Sherpa function to each of the datasets.  
            :rtype: List of results
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

        :param func: Sherpa function that takes a full parameter name specification and
                     optional args, e.g. set_par() used as set_par('mekal_7.kt', 2.0)
        :param par: Param name or model compoent name ('mekal.kt' or 'mekal')
        :param msg: Format string to indicate action.
        :param *args: Optional function arguments

        :rtype: numpy array of function return values ordered by shell
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

    def thaw(self, *pars):
        """Apply thaw command to specified parameters for each dataset.

        :param pars: parameter specifiers in format <model_type>.<par_name>
        :rtype: None
        """
        for par in pars:
            self._sherpa_par(ui.thaw, par, 'Thawing %s')

    def freeze(self, *pars):
        """Apply freeze command to specified parameters for each dataset.

        :param pars: parameter specifiers in format <model_type>.<par_name>
        :rtype: None
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

    def get_par(self, par):
        """Get parameter attribute value for each dataset.

        :param par: parameter specifier in format <model_type>.<par_name>
        :rtype: None
        """
        parname, attr, pars = self._get_parname_attr_pars(par, 'Getting %%s.%s')
        return numpy.array([getattr(x, attr) for x in pars])

    def set_par(self, par, val):
        """Set parameter attribute value for each dataset.

        :param par: parameter spec: <model_type>.<attr> or <model_type>.<par>.<attr>
        :param val: parameter value
        :rtype: None
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
            """Fit simultaneously all the datasets in the stack using the current
            source models. 

            :args: additional args that get passed to the sherpa fit() routine
            :kwargs: additional keyword args that get passed to the sherpa fit() routine
            :rtype: None
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

        Here is the text for index.rst:

        In the ``print_window`` command if a filename is supplied (for saving to a
        set of files) it should have a ``#`` character which will be replaced by the
        dataset ``id`` in each case.

        :param args: list arguments to pass to print_window
        :param kwargs: named (keyword) arguments to pass to print_window
        :rtype: None
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
            """Call Sherpa plot ``func`` for each dataset.

            :param func: Sherpa plot function
            :param args: plot function list arguments
            :param kwargs: plot function named (keyword) arguments
            :rtype: None
            """
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
            """Call matplotlib plot ``func`` for each dataset.

            :param func: Sherpa plot function
            :param args: plot function list arguments
            :param kwargs: plot function named (keyword) arguments
            :rtype: None
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
        output = []
        for dataset in self.filter_datasets():
            id = dataset['id']
            if func(ui.get_data(id)):
                output.append(id)
        return output

    def query_by_header_keyword(self, keyword, value):
        def func(dataset):
            if hasattr(dataset, 'header'):
                if keyword in dataset.header.keys() and dataset.header[keyword] == str(value):
                    return True
            return False
        return self.query(func)

    def query_by_obsid(self, value):
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
    DATASTACK.clear_models()
    DATASTACK.clear_stack()
    ui.clean()
    logger.warning("clean() will invalidate any existing DataStack instances by removing all the datasets from the " +
                   "Sherpa session")