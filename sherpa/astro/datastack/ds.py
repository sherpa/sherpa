from __future__ import absolute_import
from __future__ import print_function
#
# Copyright (C) 2015, 2016, 2019  Smithsonian Astrophysical Observatory
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

import numpy
from sherpa.utils.logging import config_logger
from sherpa.utils.err import IOErr
from sherpa.astro import ui
from .utils import load_error_msg, load_wrapper, model_wrapper, \
    simple_wrapper, fit_wrapper, plot_wrapper
from . import plot_backend as backend

logger = config_logger(__name__)

try:
    import stk
except ImportError:
    logger.warning("could not import stk library. CIAO stack files and syntax will be disabled")

# Global list of dataset ids in use
_all_dataset_ids = {}


class DataStack(object):

    """Manipulate a stack of data in Sherpa.
    """

    _default_instantiated = False

    def __init__(self):
        if not DataStack._default_instantiated\
                or (hasattr(self, "_default_instance") and self._default_instance):
            DataStack._default_instantiated = True
            self._default_instance = True
        else:
            self._default_instance = False
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
            except KeyError:
                pass

    @property
    def ids(self):
        """List of ids corresponding to stack datasets.

        Returns
        -------
        ids : array of int or str
           The identifiers of all the datasets in the stack.

        """
        return [x['id'] for x in self.datasets]

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

    # QUS: should this use the logger instead of print?
    def show_stack(self):
        """Display basic information about the contents of the data stack.

        A brief summary of each data set in the stack is printed to the
        standard output channel. The information displayed depends on the
        type of each data set.
        """

        # The observation-date (MJD) format keyword was MJD_OBS but
        # has changed (CIAO 4.12) for Chandra data to MJD-OBS, so
        # need to look for both. The idea is to use the new key
        # and then fall back to the old one.
        #
        for dataset in self.filter_datasets():
            obsid = "N/A"
            time = "N/A"
            timekey = None
            if hasattr(dataset['data'], 'header'):
                try:
                    obsid = dataset['data'].header['OBS_ID']
                except KeyError:
                    pass

                for key in ['MJD-OBS', 'MJD_OBS']:
                    try:
                        time = dataset['data'].header[key]
                        timekey = key
                        break
                    except KeyError:
                        pass

            # This is informational, so do not necessarily need
            # to report the actual key, but for now try to do so.
            #
            if timekey is None:
                timekey = 'MJD-OBS'

            print('{0}: {1} {2}: {3} {4}: {5}'.format(dataset['id'],
                                                      dataset['data'].name,
                                                      'OBS_ID', obsid,
                                                      timekey, time))

    # QUS: should this return a copy of the list?
    def get_stack_ids(self):
        """Return the data set identifiers in the stack.

        Returns
        -------
        ids : array of int or str
           The data set identifiers.
        """
        return self.ids

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
            raise AttributeError(
                "load_arrays takes at least one argument (got none).")

        if not hasattr(args[0], '__iter__'):
            id, args = args[0], args[1:]
        else:
            id = None

        if id is not None:
            if self._default_instance:
                ui.load_arrays(id, *args, **kwargs)
                return
            else:
                raise AttributeError(load_error_msg(id))

        # Array Stack.
        for arrays in args[0]:
            dataid = self._get_dataid()
            logger.info('Loading dataset id %s' % dataid)
            ui.load_arrays(dataid, *arrays)
            self._add_dataset(dataid)

    # DOC-TODO This docstring can probably be expanded
    def load_pha(self, id, arg=None, use_errors=False):
        """Load multiple data arrays.

        This extends ``sherpa.astro.ui.load_arrays`` to load multiple
        data sets with one call.

        The usual ``filename`` argument can be a stack file with multiple
        data files defined in it. In this case, the load function will be
        called as many times as datasets are included in the stack file.
        """
        if arg is None:
            id, arg = arg, id

        if id is not None:
            if self._default_instance:
                ui.load_pha(id, arg, use_errors)
                return
            else:
                raise AttributeError(load_error_msg(id))

        # File Stacks. If the file argument is a stack file, expand the
        # file and call this function for each file in the stack.
        try:
            for infile in stk.build(arg):
                self._load_func(ui.load_pha, infile, use_errors)
        except (NameError, OSError, IOErr):
            self._load_func(ui.load_pha, arg, use_errors)

    def thaw(self, *pars):
        """Apply the thaw command to specified parameters for each dataset.

        Parameters
        ----------
        pars
           The names of the parameters or model components to thaw.
        """
        for par in pars:
            self._sherpa_par(ui.thaw, par, 'Thawing %s')

    def freeze(self, *pars):
        """Apply the freeze command to specified parameters for each dataset.

        Parameters
        ----------
        pars
           The names of the parameters or model components to freeze.
        """
        for par in pars:
            self._sherpa_par(ui.freeze, par, 'Freezing %s')

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
            raise ValueError(
                'IDs = {0} not contained in dataset IDs = {1}'.format(filter_ids, self.ids))

    def get_par(self, par):
        """Get parameter object for each dataset.

        Parameters
        ----------
        par : str
           The name of the parameters, incuding the model component name.
        """
        parname, attr, pars = self._get_parname_attr_pars(par,
                                                          'Getting %%s.%s')
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
        parname, attr, pars = self._get_parname_attr_pars(par,
                                                          'Setting %%%%s.%%s = %s' % val)
        for x in pars:
            setattr(x, attr, val)

    # QUS: can there ever be 1 or 0 datasets?
    def link(self, par):
        """Link the parameter across the data stacks.

        Link all occurences of the model parameter to the value in the
        first dataset in the stack.

        Parameters
        ----------
        par : string
           Parameter name, in the format "<model_type>.<par_name>".
        """
        datasets = self.filter_datasets()
        name, parname = par.split('.')

        def _name(dset):
            return dset['model_comps'][name]['model_name']

        fullparname0 = '{0}.{1}'.format(_name(datasets[0]), parname)
        for dataset in datasets[1:]:
            fullparname = '{0}.{1}'.format(_name(dataset), parname)
            if fullparname != fullparname0:
                logger.info('Linking {0} => {1}'.format(fullparname,
                                                        fullparname0))
                ui.link(fullparname, fullparname0)

    def unlink(self, par):
        """Unlink the parameter across the data stacks.

        Parameters
        ----------
        par : string
           Parameter name, in the format "<model_type>.<par_name>".
        """
        self._sherpa_par(ui.unlink, par, 'Unlinking %s')

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
                str_value = str(value)
                try:  # Python 3
                    header_value = str(dataset.header[keyword], "utf-8")
                except TypeError:  # Python 2
                    header_value = dataset.header[keyword]
                if keyword in dataset.header.keys() and \
                   header_value == str_value:
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

    load_ascii = load_wrapper(ui.load_ascii)
    load_data = load_wrapper(ui.load_data)
    load_bkg = load_wrapper(ui.load_bkg)
    set_source = model_wrapper(ui.set_source)
    set_model = model_wrapper(ui.set_model)
    set_bkg_model = model_wrapper(ui.set_bkg_model)
    set_full_model = model_wrapper(ui.set_full_model)
    set_bkg_full_model = model_wrapper(ui.set_bkg_full_model)
    subtract = simple_wrapper(ui.subtract)
    notice = simple_wrapper(ui.notice_id)
    ignore = simple_wrapper(ui.ignore_id)
    get_arf = simple_wrapper(ui.get_arf)
    get_rmf = simple_wrapper(ui.get_rmf)
    get_response = simple_wrapper(ui.get_response)
    get_bkg_arf = simple_wrapper(ui.get_bkg_arf)
    get_bkg_rmf = simple_wrapper(ui.get_bkg_rmf)
    get_bkg = simple_wrapper(ui.get_bkg)
    get_source = simple_wrapper(ui.get_source)
    get_model = simple_wrapper(ui.get_model)
    get_bkg_model = simple_wrapper(ui.get_bkg_model)
    get_bkg_scale = simple_wrapper(ui.get_bkg_scale)
    group_adapt = simple_wrapper(ui.group_adapt)
    group_adapt_snr = simple_wrapper(ui.group_adapt_snr)
    group_bins = simple_wrapper(ui.group_bins)
    group_counts = simple_wrapper(ui.group_counts)
    group_snr = simple_wrapper(ui.group_snr)
    group_width = simple_wrapper(ui.group_width)
    load_arf = simple_wrapper(ui.load_arf)
    load_rmf = simple_wrapper(ui.load_rmf)
    load_bkg_arf = simple_wrapper(ui.load_bkg_arf)
    load_bkg_rmf = simple_wrapper(ui.load_bkg_rmf)
    load_filter = simple_wrapper(ui.load_filter)
    load_grouping = simple_wrapper(ui.load_grouping)
    fit_bkg = fit_wrapper(ui.fit_bkg)
    fit = fit_wrapper(ui.fit)
    conf = fit_wrapper(ui.conf)
    plot_arf = plot_wrapper(ui.plot_arf)
    plot_bkg_fit = plot_wrapper(ui.plot_bkg_fit)
    plot_bkg_ratio = plot_wrapper(ui.plot_bkg_ratio)
    plot_chisqr = plot_wrapper(ui.plot_chisqr)
    plot_fit_delchi = plot_wrapper(ui.plot_fit_delchi)
    plot_psf = plot_wrapper(ui.plot_psf)
    plot_bkg = plot_wrapper(ui.plot_bkg)
    plot_bkg_fit_delchi = plot_wrapper(ui.plot_bkg_fit_delchi)
    plot_bkg_resid = plot_wrapper(ui.plot_bkg_resid)
    plot_data = plot_wrapper(ui.plot_data)
    plot_fit_resid = plot_wrapper(ui.plot_fit_resid)
    plot_ratio = plot_wrapper(ui.plot_ratio)
    plot_bkg_chisqr = plot_wrapper(ui.plot_bkg_chisqr)
    plot_bkg_fit_resid = plot_wrapper(ui.plot_bkg_fit_resid)
    plot_bkg_source = plot_wrapper(ui.plot_bkg_source)
    plot_delchi = plot_wrapper(ui.plot_delchi)
    plot_model = plot_wrapper(ui.plot_model)
    plot_resid = plot_wrapper(ui.plot_resid)
    plot_bkg_delchi = plot_wrapper(ui.plot_bkg_delchi)
    plot_bkg_model = plot_wrapper(ui.plot_bkg_model)
    plot_bkg_source = plot_wrapper(ui.plot_bkg_source)
    plot_fit = plot_wrapper(ui.plot_fit)
    plot_order = plot_wrapper(ui.plot_order)
    plot_source = plot_wrapper(ui.plot_source)
    plot_savefig = backend.plot_savefig
    plot_xlabel = backend.plot_xlabel
    plot_ylabel = backend.plot_ylabel
    plot_title = backend.plot_title
    plot_xlim = backend.plot_xlim
    plot_ylim = backend.plot_ylim
    plot_set_xscale = backend.plot_set_xscale
    plot_set_yscale = backend.plot_set_yscale

    def _get_dataid(self):
        if self.getitem_ids:
            dataid = self.getitem_ids[0]
            self.getitem_ids = None
        else:
            dataid = 1
            while dataid in _all_dataset_ids:
                dataid += 1

        if dataid in self.dataset_ids:
            raise ValueError(
                'Data ID = {0} is already in the DataStack'.format(dataid))

        return dataid

    def _add_dataset(self, dataid):
        dataset = dict(id=dataid, args=[], model_comps={},
                       data=ui.get_data(dataid))
        _all_dataset_ids[dataid] = dataset
        self.dataset_ids[dataid] = dataset
        self.datasets.append(dataset)

    def _load_func(self, func, *args, **kwargs):
        dataid = self._get_dataid()

        logger.info('Loading dataset id {0}'.format(dataid))
        func(dataid, *args, **kwargs)

        self._add_dataset(dataid)

    def _sherpa_par(self, func, par, msg, *args, **kwargs):
        """Apply a function to one or more  model parameters.

        The given function is called - with the given ``args``
        and ``kwargs`` to one or more parameters (the parameter
        name is the first argument of the function). See thaw(),
        freeze(), set_par() and get_par() for examples.

        Parameters
        ----------
        func
           Sherpa function that takes a full parameter name specification
           and optional args, e.g. set_par
        par : string
           Parameter or model component name (e.g. 'mekal.kt' or 'mekal').
        msg
           Format string to indicate action (it should include one "%s"
           which will get replaced by the parameter name in use).
        args
           Arguments to pass to ``func`` after the parameter name.
        kwargs
           Keyword arguments for ``func``.

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

        retvals = []
        processed = set()
        for dataset in self.filter_datasets():
            model_comps = dataset['model_comps']
            if name in model_comps:
                model_name = model_comps[name]['model_name']
                if parname:
                    fullparname = '{0}.{1}'.format(model_name, parname)
                else:
                    fullparname = model_name
                if fullparname not in processed:
                    if msg is not None:
                        logger.info(msg % fullparname)
                    retvals.append(func(fullparname, *args, **kwargs))
                    processed.add(fullparname)

        return retvals

    def _get_parname_attr_pars(self, par, msg):
        parts = par.split('.')
        if len(parts) == 1:
            raise ValueError(
                'par="%s" must be in the form "name.par" or "name.par.attr"' % par)
        parname = '.'.join(parts[:-1])
        attr = parts[-1]
        return parname, attr, self._sherpa_par(eval, parname, msg % attr)
