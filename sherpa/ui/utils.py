# 
#  Copyright (C) 2010  Smithsonian Astrophysical Observatory
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

import copy
import copy_reg
import cPickle as pickle
import inspect
from itertools import izip
import logging
import sys
import os
import re
import numpy
import sherpa.all
from sherpa.utils import SherpaFloat, NoNewAttributesAfterInit, export_method
from sherpa.utils.err import *

info = logging.getLogger(__name__).info
warning = logging.getLogger(__name__).warning

from sherpa import get_config
from ConfigParser import ConfigParser

import readline, inspect

config = ConfigParser()
config.read(get_config())
# Suppress printing of traceback in high-level UI.
# If needed for debugging, set sys.tracebacklimit to a
# nonzero value in that session--traceback is for debugging
# but is irritating to ordinary users.
sys.tracebacklimit = int(config.get('verbosity','level'))
numpy.set_printoptions(threshold=int(config.get('verbosity','arraylength')))


_builtin_symbols_ = sys.modules["__builtin__"].__dict__.keys()

__all__ = ('ModelWrapper', 'Session')



###############################################################################
#
# Errors and argument checking
#
###############################################################################


def _argument_type_error(argname, argdesc):
    raise ArgumentTypeErr('badarg', argname, argdesc)


def _check_type(arg, argtype, argname, argdesc, nottype=None):
    if ((not isinstance(arg, argtype)) or
        ((nottype is not None) and isinstance(arg, nottype))):
        _argument_type_error(argname, argdesc)


def _is_integer(val):
    return isinstance(val, (int, numpy.integer))


def _fix_array(arg, argname, ndims=1):
    try:
        arg = numpy.asarray(arg, SherpaFloat)
        if len(arg.shape) != ndims:
            raise TypeError
    except TypeError:
        _argument_type_error(argname, 'a %d-D array' % ndims)

    return arg


def _is_subclass(t1, t2):
    return isinstance(t1, type) and issubclass(t1, t2) and (t1 is not t2)


def _send_to_pager(all, filename=None, clobber=False):
    pager = None
    clobber=sherpa.utils.bool_cast(clobber)
    try:
        if filename is None:
            if (os.environ.has_key('PAGER') == True):
                pager = os.popen(os.environ['PAGER'], 'w')
            else:
                if (os.access('/bin/more', os.X_OK) == 1):
                    pager = os.popen('/bin/more', 'w')
                elif (os.access('/usr/bin/more', os.X_OK) == 1):
                    pager = os.popen('/usr/bin/more', 'w')
                else:
                    raise IOErr('nopager')

        # check if filename is StringIO obj
        elif hasattr(filename, 'write'):
            pager = filename
            print >> pager, all
            return

        else:
            _check_type(filename, basestring, 'filename', 'a string')
            if os.path.isfile(filename) and not clobber:
                raise IOErr('filefound', filename)
            pager = file(filename, 'w')
        print >> pager, all
    except:
        if (pager is not None):
            pager.close()
        raise
    else:
        if (pager is not None):
            pager.close()


###############################################################################
#
# Pickling support
#
###############################################################################



def construct_ufunc(modname, funcname):
    module = sys.modules.get(modname)
    if module is None:
        module = __import__(modname)
    return getattr(module, funcname)


def reduce_ufunc(func):
    modname  = getattr(func, '__module__', 'numpy')
    funcname = func.__name__
    if func is not getattr(sys.modules[modname], funcname, None):
        raise ValueError("module '%s' does not contain ufunc '%s'" %
                         (modname, funcname))
    return (construct_ufunc, (modname, funcname))


copy_reg.constructor(construct_ufunc)
copy_reg.pickle(numpy.ufunc, reduce_ufunc)



###############################################################################
#
# ModelWrapper
#
###############################################################################


class ModelWrapper(NoNewAttributesAfterInit):

    def __init__(self, session, modeltype, args=(), kwargs={}):
        self._session = session
        self.modeltype = modeltype
        self.args = args
        self.kwargs = kwargs
        NoNewAttributesAfterInit.__init__(self)

    def __call__(self, name):
        _check_type(name, basestring, 'name', 'a string')

        m = self._session._get_model_component(name)
        if (m is not None) and isinstance(m, self.modeltype):
            return m

        m = self.modeltype(name, *self.args, **self.kwargs)
        self._session._add_model_component(m)
        return m

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError("'%s\' object has no attribute '%s'" %
                                 (type(self).__name__, name))
        return self(name)

    def __repr__(self):
        return '<%s model type>' % self.modeltype.__name__

    def __str__(self):
        if self.modeltype.__doc__ is not None:
            return self.modeltype.__doc__
        return self.__repr__()


def _assign_obj_to_main(name, obj):
    sys.modules["__main__"].__dict__[name] = obj
    sys.modules["__builtin__"].__dict__[name] = obj


def _assign_model_to_main(name, model):
    # Ask sys what the __main__ module is; packages such
    # as IPython can add their own __main__ module.
    model.name = '%s.%s' % (type(model).__name__.lower(),name)
    _assign_obj_to_main(name, model)


###############################################################################
#
# Session
#
###############################################################################
class loggable(object):
	def __init__(self, with_id=False, with_keyword=False, with_name=None):
		self.with_id = with_id
		self.with_keyword = with_keyword
		self.with_name = with_name

	def __call__(self, func):
		if self.with_name:
			name = self.with_name
		else:
			name = func.__name__
		def log_decorator(*args, **kwargs):
			ret = func(*args, **kwargs)
			session = args[0]
			line = readline.get_history_item(readline.get_current_history_length())
			if self.with_id:
				the_args = inspect.getcallargs(func, *args, **kwargs)
				id = the_args['id']
				if self.with_keyword:    					model = the_args[self.with_keyword]    					if model is None:    						id = None  
				id = session._fix_id(id)  			    				if id is not None: # otherwise don't do anything and let normal error handling take action					if not session._calls_tracker.has_key(id):						session._calls_tracker[id] = dict()    					session._calls_tracker[id][name] = line
			else:    				session._calls_tracker[name] = line
            		return ret		log_decorator._original = func # this is needed because __init__.py will recreate the methods, see that file for info (look up 'decorator')        	return log_decorator


class Session(NoNewAttributesAfterInit):


    ###########################################################################
    # Standard methods
    ###########################################################################


    def __init__(self):
        self.clean()
        self._model_types = {}
        self._model_globals = numpy.__dict__.copy()
	self._calls_tracker = dict()
        NoNewAttributesAfterInit.__init__(self)
	global _session
	_session = self

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_model_globals']
        return state

    def __setstate__(self, state):
        self._model_globals = numpy.__dict__.copy()

        self._model_globals.update(state['_model_types'])

        if not state.has_key('_sources'):
            self.__dict__['_sources'] = state.pop('_models')

        self.__dict__.update(state)

    ###########################################################################
    # High-level utilities
    ###########################################################################


    def _export_names(self, gdict):
        allnames = []

        for name in dir(self):
            if ((not name.startswith('_')) or
                (name is '_sherpa_version') or
                (name is '_sherpa_version_string')):
                gdict[name] = export_method(getattr(self, name),
                                            modname=gdict.get('__name__'))
                allnames.append(name)

        gdict.update(self._model_types)
        allnames.extend(self._model_types.keys())

        return allnames

    ### Ahelp ingest: 2015-04-27 DJB
    def clean(self):
        """Clear out the current Sherpa session.

        The `clean` function removes all data sets and model
        assignments, and restores the default settings for the
        optimisation and fit statistic.

        See Also
        --------
        save : Save the current Sherpa session to a file.
        restore : Load in a Sherpa session from a file.
        sherpa.astro.ui.utils.save_all : Save the Sherpa session as an ASCII file.

        Examples
        --------

        >>> clean()

        """
        # Make version, version string names and formats identical
        # between Python and S-Lang interfaces.
        self._sherpa_version = sherpa.__version__
        self._sherpa_version_string = sherpa.__versionstr__
        
        self._default_id = 1
        self._paramprompt=False

        self._methods = {}
        self._itermethods = {'none': {'name':'none'},
                             'primini': {'name':'primini',
                                         'maxiters':10,
                                         'tol':1.0e-3},
                             'sigmarej': {'name':'sigmarej',
                                          'maxiters':5,
                                          'hrej':3,
                                          'lrej':3,
                                          'grow':0}}

        self._stats = {}
        self._estmethods = {}

        modules   = (sherpa.optmethods, sherpa.stats, sherpa.estmethods)
        basetypes = (sherpa.optmethods.OptMethod, sherpa.stats.Stat,
                     sherpa.estmethods.EstMethod)
        objdicts  = (self._methods, self._stats, self._estmethods)

        for mod, base, odict in izip(modules, basetypes, objdicts):
            for name in mod.__all__:
                cls = getattr(mod, name)
                if _is_subclass(cls, base):
                    odict[name.lower()] = cls()

        self._current_method = self._methods['levmar']
        self._current_itermethod = self._itermethods['none']
        self._current_stat   = self._stats['chi2gehrels']
        # Add simplex as alias to neldermead
        self._methods['simplex'] = self._methods['neldermead']

        self._data = {}
        self._psf = {}
        self._tbl_models = []
        self._psf_models = []

        self._model_autoassign_func = _assign_model_to_main
        self._model_components = {}
        self._models = {}
        self._sources = {}

        self._fit_results = None
        self._pvalue_results = None

        self._covariance_results = None
        self._confidence_results = None        
        self._projection_results = None

        self._pyblocxs = sherpa.sim.MCMC()

        self._splitplot = sherpa.plot.SplitPlot()
        self._jointplot = sherpa.plot.JointPlot()
        self._dataplot  = sherpa.plot.DataPlot()
        self._modelplot = sherpa.plot.ModelPlot()

        self._compmdlplot = sherpa.plot.ComponentModelPlot()
        self._compsrcplot = sherpa.plot.ComponentSourcePlot()
        #self._comptmplmdlplot = sherpa.plot.ComponentTemplateModelPlot()
        self._comptmplsrcplot = sherpa.plot.ComponentTemplateSourcePlot()

        self._sourceplot= sherpa.plot.SourcePlot()
        self._fitplot   = sherpa.plot.FitPlot()
        self._residplot = sherpa.plot.ResidPlot()
        self._delchiplot= sherpa.plot.DelchiPlot()
        self._chisqrplot= sherpa.plot.ChisqrPlot()
        self._ratioplot = sherpa.plot.RatioPlot()
        self._psfplot = sherpa.plot.PSFPlot()
        self._kernelplot = sherpa.plot.PSFKernelPlot()
        self._lrplot = sherpa.plot.LRHistogram()
        self._pdfplot = sherpa.plot.PDFPlot()
        self._cdfplot = sherpa.plot.CDFPlot()
        self._traceplot = sherpa.plot.TracePlot()
        self._scatterplot = sherpa.plot.ScatterPlot()

        self._datacontour  = sherpa.plot.DataContour()
        self._modelcontour = sherpa.plot.ModelContour()
        self._sourcecontour= sherpa.plot.SourceContour()
        self._fitcontour   = sherpa.plot.FitContour()
        self._residcontour = sherpa.plot.ResidContour()
        self._ratiocontour = sherpa.plot.RatioContour()
        self._psfcontour   = sherpa.plot.PSFContour()
        self._kernelcontour = sherpa.plot.PSFKernelContour()

        self._intproj = sherpa.plot.IntervalProjection()
	self._intunc  = sherpa.plot.IntervalUncertainty()
	self._regproj = sherpa.plot.RegionProjection()
	self._regunc  = sherpa.plot.RegionUncertainty()

        self._plot_types = {
            'data': self._dataplot,
            'model': self._modelplot,
            'source': self._sourceplot,
            'fit': self._fitplot,
            'resid': self._residplot,
            'ratio': self._ratioplot,
            'delchi': self._delchiplot,
            'chisqr': self._chisqrplot,
            'psf': self._psfplot,
            'kernel': self._kernelplot,
            'compsource' : self._compsrcplot,
            'compmodel' : self._compmdlplot
            }

        self._contour_types = {
            'data': self._datacontour,
            'model': self._modelcontour,
            'source': self._sourcecontour,
            'fit': self._fitcontour,
            'resid': self._residcontour,
            'ratio': self._ratiocontour,
            'psf' : self._psfcontour,
            'kernel': self._kernelcontour
            }
        
        self._dataimage = sherpa.image.DataImage()
        self._modelimage = sherpa.image.ModelImage()
        self._sourceimage = sherpa.image.SourceImage()
        self._ratioimage = sherpa.image.RatioImage()
        self._residimage = sherpa.image.ResidImage()
        self._psfimage  = sherpa.image.PSFImage()
        self._kernelimage  = sherpa.image.PSFKernelImage()
        self._mdlcompimage = sherpa.image.ComponentModelImage()
        self._srccompimage = sherpa.image.ComponentSourceImage()


    ### Ahelp ingest: 2015-04-27 DJB
    def save(self, filename='sherpa.save', clobber=False):
        """Save the current Sherpa session to a file.

        Parameters
        ----------
        filename : str, optional
           The name of the file to write the results to. The default
           is `sherpa.save`.
        clobber : bool, optional
           This flag controls whether an existing file can be
           overwritten (`True`) or if it raises an exception (`False`,
           the default setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        clean : Clear all stored session data.
        restore : Load in a Sherpa session from a file.
        sherpa.astro.ui.utils.save_all : Save the Sherpa session as an ASCII file.

        Notes
        -----
        The current Sherpa session is saved using the Python `pickle`
        module. The output is a binary file, which may not be portable
        between versions of Sherpa, but is platform independent, and
        contains all the data. This means that files created by `save`
        can be sent to collaborators to share results.

        Examples
        --------

        Save the current session to the file `sherpa.save`.

        >>> save()

        Save the current session to the file `bestfit.sherpa`,
        overwriting any existing version of the file.

        >>> save('bestfit.sherpa', clobber=True)

        """
        
        _check_type(filename, basestring, 'filename', 'a string')
        clobber=sherpa.utils.bool_cast(clobber)
        
        if os.path.isfile(filename) and not clobber:
            raise sherpa.utils.err.IOErr("filefound", filename)

        fout = file(filename, 'wb')
        try:
            pickle.dump(self, fout, 2)  # Use newer binary protocol
        finally:
            fout.close()

    ### Ahelp ingest: 2015-04-27 DJB
    def restore(self, filename='sherpa.save'):
        """Load in a Sherpa session from a file.

        Parameters
        ----------
        filename : str, optional
           The name of the file to read the results from. The default
           is `sherpa.save`.

        Raises
        ------
        IOError
           If `filename` does not exist.

        See Also
        --------
        clean : Clear all stored session data.
        save : Save the current Sherpa session to a file.

        Notes
        -----

        The input to `restore` must have been created with the `save`
        command. This is a binary file, which may not be portable
        between versions of Sherpa, but is platform independent. A
        warning message may be created if a file saved by an older
        (or newer) version of Sherpa is loaded. An example of such
        a message is

        WARNING: Could not determine whether the model is discrete.
        This probably means that you have restored a session saved with a previous version of Sherpa.
        Falling back to assuming that the model is continuous.

        Examples
        --------

        Load in the Sherpa session from `sherpa.save`.

        >>> restore()

        Load in the session from the given file:

        >>> restore('/data/m31/setup.sherpa')

        """
        _check_type(filename, basestring, 'filename', 'a string')

        fin = file(filename, 'rb')
        try:
            obj = pickle.load(fin)
        finally:
            fin.close()

        if not isinstance(obj, Session):
            raise ArgumentErr('nosession', filename)

        # Update optmethods, stats, and estmethods
        # obj.__dict__ should not clobber new classes!
        dicts = [self._methods,self._stats,self._estmethods]
        names = ['_methods','_stats','_estmethods']
        for name, dic in izip(names, dicts):
            # update current session with user definitions
            dic.update(obj.__dict__[name])
            # remove old items from pickle
            obj.__dict__.pop(name)

        # update current session with pickle
        self.__dict__.update(obj.__dict__)

        if self._model_autoassign_func is not None:
            for name, cmpt in self._model_components.items():
                self._model_autoassign_func(name, cmpt)


    def _get_show_data(self, id=None):
        data_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            data_str += 'Data Set: %s\n' % id
            data_str += self.get_data(id).__str__() + '\n\n'
        return data_str

    def _get_show_filter(self, id=None):
        filt_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            filt_str += 'Data Set Filter: %s\n' % id
            filt_str += self.get_data(id).get_filter_expr() + '\n\n'
        return filt_str

    def _get_show_model(self, id=None):
        model_str = ''
        ids = self.list_data_ids()
        mdl_ids = self.list_model_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            if id in mdl_ids:
                model_str += 'Model: %s\n' % id
                model_str += self._get_model(id).__str__() + '\n\n'
        return model_str

    def _get_show_source(self, id=None):
        model_str = ''
        ids = self.list_data_ids()
        src_ids = self._sources.keys()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            if id in src_ids:
                model_str += 'Model: %s\n' % id
                model_str += self._get_source(id).__str__() + '\n\n'
        return model_str


    def _get_show_kernel(self, id=None):
        kernel_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            if id in self._psf.keys():
                kernel_str += 'PSF Kernel: %s\n' % id
                # Show the PSF parameters
                kernel_str += self.get_psf(id).__str__() + '\n\n'
        return kernel_str


    def _get_show_psf(self, id=None):
        psf_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            if id in self._psf.keys():
                psf_str += 'PSF Model: %s\n' % id
                # Show the PSF dataset or PSF model
                psf_str += self.get_psf(id).kernel.__str__() + '\n\n'
        return psf_str


    def _get_show_method(self):
        return ('Optimization Method: %s\n%s\n' %
                (type(self._current_method).__name__,
                 self._current_method.__str__()))

    def _get_show_stat(self):
        return ('Statistic: %s\n%s\n' %
                (type(self._current_stat).__name__,
                 self._current_stat.__str__()))

    def _get_show_fit(self):
        fit_str = ''
        if self._fit_results != None:
            fit_str += self._get_show_method()
            fit_str += '\n'
            fit_str += self._get_show_stat()
            fit_str += '\n'
            fit_str += 'Fit:' 
            fit_str += self.get_fit_results().format() + '\n\n'
        return fit_str

    def _get_show_conf(self):
        conf_str = ''
        if self._confidence_results != None:
            conf_str += 'Confidence:'
            conf_str += self.get_conf_results().format() + '\n\n'
        return conf_str

    def _get_show_proj(self):
        proj_str = ''
        if self._projection_results != None:
            proj_str += 'Projection:'
            proj_str += self.get_proj_results().format() + '\n\n'
        return proj_str

    def _get_show_covar(self):
        covar_str = ''
        if self._covariance_results != None:
            covar_str += 'Covariance:'
            covar_str += self.get_covar_results().format() + '\n\n'
        return covar_str

    ### Ahelp ingest: 2015-04-27 DJB
    def show_stat(self, outfile=None, clobber=False):
        """Display the current fit statistic.

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        calc_stat : Calculate the fit statistic for a data set.
        calc_stat_info : Display the statistic values for the current models.
        get_stat : Return a fit-statistic method.
        show_all : Report the current state of the Sherpa session.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        Examples
        --------

        >>> set_stat('cash')
        >>> show_stat()
        Statistic: Cash
        Maximum likelihood function

        """
        all = ''
        all += self._get_show_stat()
        _send_to_pager(all, outfile, clobber)


    ### Ahelp ingest: 2015-04-27 DJB
    def show_method(self, outfile=None, clobber=False):
        """Display the current optimization method and options.

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        get_method : Return an optimization method.
        get_method_opt : Return one or all options of the current optimization method.
        show_all : Report the current state of the Sherpa session.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        Examples
        --------

        >>> set_method('levmar')
        >>> show_method()
        Optimization Method: LevMar
        name    = levmar
        ftol    = 1.19209289551e-07
        xtol    = 1.19209289551e-07
        gtol    = 1.19209289551e-07
        maxfev  = x
        epsfcn  = 1.19209289551e-07
        factor  = 100.0
        verbose = 0

        """
        all = ''
        all += self._get_show_method()
        _send_to_pager(all, outfile, clobber)


    ### Ahelp ingest: 2015-05-02 DJB
    def show_fit(self, outfile=None, clobber=False):
        """Summarize the fit results.

        Display the results of the last call to `fit`, including:
        optimization method, statistic, and details of the fit (it
        does not reflect any changes made after the fit, such as to
        the model expression or fit parameters).

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        fit : Fit one or more data sets.
        get_fit_results : Return the results of the last fit.
        list_data_ids : List the identifiers for the loaded data sets.
        list_model_ids : List of all the data sets with a source expression.
        show_all : Report the current state of the Sherpa session.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_fit()
        _send_to_pager(all, outfile, clobber)


    ### Ahelp ingest: 2015-05-02 DJB
    def show_data(self, id=None, outfile=None, clobber=False):
        """Summarize the available data sets.

        Display information on the data sets that have been
        loaded. The details depend on the type of the data set
        (e.g. 1D, image, PHA files).

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then all data sets are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        list_data_ids : List the identifiers for the loaded data sets.
        show_all : Report the current state of the Sherpa session.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_data(id)
        _send_to_pager(all, outfile, clobber)


    ### Ahelp ingest: 2015-05-02 DJB
    def show_filter(self, id=None, outfile=None, clobber=False):
        """Show any filters applied to a data set.

        Display any filters that have been applied to the independent
        axis or axes of the data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then all data sets are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        ignore : Exclude data from the fit.
        sherpa.astro.utils.ignore2d : Exclude a spatial region from an image.
        list_data_ids : List the identifiers for the loaded data sets.
        notice : Include data in the fit.
        sherpa.astro.utils.notice2d : Include a spatial region of an image.
        show_all : Report the current state of the Sherpa session.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_filter(id)
        _send_to_pager(all, outfile, clobber)


    ### Ahelp ingest: 2015-05-01 DJB
    def show_model(self, id=None, outfile=None, clobber=False):
        """Display the model expression used to fit a data set.

        This displays the model used to fit the data set, that is
        that is, the expression set by `set_model` or `set_source`
        combined with any instrumental responses, together with the
        parameter values of the model. The `show_source` function
        displays just the source expression, without the instrumental
        components (if any).

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then all source expressions are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        list_model_ids : List of all the data sets with a source expression.
        set_model : Set the source model expression for a data set.
        show_all : Report the current state of the Sherpa session.
        show_source : Display the source model expression for a data set.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        #all += self._get_show_kernel(id)
        all += self._get_show_psf(id)
        all += self._get_show_model(id)
        _send_to_pager(all, outfile, clobber)


    ### Ahelp ingest: 2015-05-01 DJB
    def show_source(self, id=None, outfile=None, clobber=False):
        """Display the source model expression for a data set.

        This displays the source model for a data set, that is, the
        expression set by `set_model` or `set_source`, as well as the
        parameter values for the model. The `show_model` function
        displays the model that is fit to the data; that is, it
        includes any instrument responses.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then all source expressions are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        list_model_ids : List of all the data sets with a source expression.
        set_model : Set the source model expression for a data set.
        show_all : Report the current state of the Sherpa session.
        show_model : Display the model expression used to fit a data set.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_source(id)
        _send_to_pager(all, outfile, clobber)

    ### Ahelp ingest: 2015-05-02 DJB
    ### DOC-TODO: how and where to describe the PSF/kernel difference
    ###           as the Notes section below is inadequate
    def show_kernel(self, id=None, outfile=None, clobber=False):
        """Display any kernel applied to a data set.

        The kernel represents the subset of the PSF model that is used
        to fit the data. The `show_psf` function shows the un-filtered
        version.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then all data sets are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        image_kernel : Plot the 2D kernel applied to a data set.
        list_data_ids : List the identifiers for the loaded data sets.
        load_psf : Create a PSF model.
        plot_kernel : Plot the 1D kernel applied to a data set.
        set_psf : Add a PSF model to a data set.
        show_all : Report the current state of the Sherpa session.
        show_psf : Display any PSF model applied to a data set.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        The point spread function (PSF) is defined by the full
        (unfiltered) PSF image or model expression evaluated over the
        full range of the dataset; both types of PSFs are established
        with `load_psf`.  The kernel is the subsection of the PSF
        image or model which is used to convolve the data: this is
        changed using `set_psf`.  While the kernel and PSF might be
        congruent, defining a smaller kernel helps speed the
        convolution process by restricting the number of points within
        the PSF that must be evaluated.

        """
        all = ''
        all += self._get_show_kernel(id)
        _send_to_pager(all, outfile, clobber)

    ### Ahelp ingest: 2015-05-02 DJB
    ### DOC-TODO: how and where to describe the PSF/kernel difference
    ###           as the Notes section below is inadequate
    def show_psf(self, id=None, outfile=None, clobber=False):
        """Display any PSF model applied to a data set.

        The PSF model represents the full model or data set that is
        applied to the source expression. The `show_kernel` function
        shows the filtered version.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then all data sets are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        image_psf : View the 2D PSF model applied to a data set.
        list_data_ids : List the identifiers for the loaded data sets.
        load_psf : Create a PSF model.
        plot_psf : Plot the 1D PSF model applied to a data set.
        set_psf : Add a PSF model to a data set.
        show_all : Report the current state of the Sherpa session.
        show_kernel : Display any kernel applied to a data set.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        The point spread function (PSF) is defined by the full
        (unfiltered) PSF image or model expression evaluated over the
        full range of the dataset; both types of PSFs are established
        with `load_psf`.  The kernel is the subsection of the PSF
        image or model which is used to convolve the data: this is
        changed using `set_psf`.  While the kernel and PSF might be
        congruent, defining a smaller kernel helps speed the
        convolution process by restricting the number of points within
        the PSF that must be evaluated.

        """
        all = ''
        #all += self._get_show_kernel(id)
        all += self._get_show_psf(id)
        _send_to_pager(all, outfile, clobber)

    ### Ahelp ingest: 2015-05-01 DJB
    def show_conf(self, outfile=None, clobber=False):
        """Display the results of the last conf evaluation.

        The output includes the best-fit model parameter values,
        associated confidence limits, choice of statistic, and details
        on the best fit location.

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        conf : Estimate confidence intervals using the confidence method.
        show_all : Report the current state of the Sherpa session.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_conf()
        _send_to_pager(all, outfile, clobber)

    ### Ahelp ingest: 2015-05-01 DJB
    def show_proj(self, outfile=None, clobber=False):
        """Display the results of the last proj evaluation.

        The output includes the best-fit model parameter values,
        associated confidence limits, choice of statistic, and details
        on the best fit location.

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        proj : Estimate confidence intervals using the projection method.
        show_all : Report the current state of the Sherpa session.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_proj()
        _send_to_pager(all, outfile, clobber)

    ### Ahelp ingest: 2015-05-01 DJB
    def show_covar(self, outfile=None, clobber=False):
        """Display the results of the last covar evaluation.

        The output includes the best-fit model parameter values,
        associated confidence limits, choice of statistic, and details
        on the best fit location.

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        covar : Estimate confidence intervals using the covariance method.
        show_all : Report the current state of the Sherpa session.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_covar()
        _send_to_pager(all, outfile, clobber)


    ### Ahelp ingest: 2015-05-02 DJB
    def show_all(self, id=None, outfile=None, clobber=False):
        """Report the current state of the Sherpa session.

        Display information about one or all of the data sets that
        have been loaded into the Sherpa session. The information
        shown includes that provided by the other `show_xxx` routines,
        and depends on the type of data that is loaded.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then all data sets are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        clean : Clear all stored session data.
        list_data_ids : List the identifiers for the loaded data sets.
        save : Save the current Sherpa session to a file.
        sherpa.astro.ui.utils.save_all : Save the Sherpa session as an ASCII file.
        sherpa.astro.ui.show_bkg
        sherpa.astro.ui.show_bkg_model
        sherpa.astro.ui.show_bkg_source
        show_conf
        show_covar
        show_data
        show_filter
        show_fit
        show_kernel
        show_method
        show_model
        show_proj
        show_psf
        show_source
        show_stat

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_data(id)
        all += self._get_show_model(id)
        all += self._get_show_fit()
        all += self._get_show_conf()
        all += self._get_show_proj()
        all += self._get_show_covar()
        _send_to_pager(all, outfile, clobber)

    ### Ahelp ingest: 2015-05-04 DJB
    def get_functions(self):
        """Return the functions provided by Sherpa.

        Returns
        -------
        functions : list of str

        See Also
        --------
        list_functions : Display the functions provided by Sherpa.

        """
        funcs = []
        for func in dir(self):
            if not func.startswith('_') and callable(getattr(self,func)):
                funcs.append(func)
        return funcs
        
    ### Ahelp ingest: 2015-04-27 DJB
    def list_functions(self, outfile=None, clobber=False):
        """Display the functions provided by Sherpa.

        Unlike the other `list_` commands, this does not
        return an array. Instead it acts like the `show_`
        family of commands. 

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        get_functions : Return the functions provided by Sherpa.
        show_all : Report the current state of the Sherpa session.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        funcs_list = self.get_functions()
        funcs = ''
        for func in funcs_list:
            funcs += '%s\n' % func

        _send_to_pager(funcs, outfile, clobber)


    ###########################################################################
    # IDs and general data management
    ###########################################################################


    @staticmethod
    def _valid_id(id):
        return (_is_integer(id) or isinstance(id, basestring))

    def _fix_id(self, id):
        if id is None:
            return self._default_id
        if not self._valid_id(id):
            raise ArgumentTypeErr('intstr')
        if id in self._plot_types.keys() or id in self._contour_types.keys():
            raise IdentifierErr('badid', id)
        return id

    def _get_item(self, id, itemdict, itemdesc, errdesc):
        id = self._fix_id(id)
        item = itemdict.get(id)
        if item is None:
            raise IdentifierErr('getitem', itemdesc, id, errdesc)
        return item

    def _set_item(self, id, item, itemdict, itemtype, itemname, itemdesc):
        id = self._fix_id(id)
        _check_type(item, itemtype, itemname, itemdesc)
        itemdict[id] = item

    ### Ahelp ingest: 2015-04-28 DJB
    def get_default_id(self):
        """Return the default data set identifier.

        The Sherpa data id ties data, model, fit, and plotting
        information into a data set easily referenced by id. The
        default identifier, used by many commands, is returned by this
        command and can be changed by `set_default_id`.

        Returns
        -------
        id : int or str
           The default data set identifier used by certain Sherpa
           functions when an identifier is not given, or set to
           `None`.

        See Also
        --------
        list_data_ids : List the identifiers for the loaded data sets.
        set_default_id : Set the default data set identifier.

        Notes
        -----
        The default Sherpa data set identifier is the integer
        `1`.

        """
        return self._default_id

    ### Ahelp ingest: 2015-04-28 DJB
    def set_default_id(self, id):
        """Set the default data set identifier.

        The Sherpa data id ties data, model, fit, and plotting
        information into a data set easily referenced by id. The
        default identifier, used by many commands, is changed by this
        command. The current setting can be found by using
        `get_default_id`.

        Parameters
        ----------
        id : int or str
           The default data set identifier to be used by certain
           Sherpa functions when an identifier is not given, or set to
           `None`.

        See Also
        --------
        get_default_id : Return the default data set identifier.
        list_data_ids : List the identifiers for the loaded data sets.

        Notes
        -----
        The default Sherpa data set identifier is the integer
        `1`.

        Examples
        --------

        After the following, many commands, such as `set_source`, will
        use `src` as the default data set identifier:

        >>> set_default_id('src')

        Restore the default data set identifier.

        >>> set_default_id(1)

        """
        self._default_id = self._fix_id(id)


    ###########################################################################
    # Optimization methods
    ###########################################################################


    ### Ahelp ingest: 2015-04-27 DJB
    def list_methods(self):
        """List the optimization methods.

        Returns
        -------
        methods : list of str
           A list of the names that can be used with
           `set_method`.

        See Also
        --------
        get_method_name : Return the name of the current optimization method.
        set_method : Set the optimization method.

        Examples
        --------

        >>> list_methods()
        ['gridsearch', 'levmar', 'moncar', 'neldermead', 'simplex']

        """
        keys = self._methods.keys()[:]
        keys.sort()
        return keys

    def _get_method_by_name(self, name):
        meth = self._methods.get(name.lower())
        if meth is None:
            raise ArgumentErr('badmethod', name)
        return meth

    ### Ahelp ingest: 2015-04-27 DJB
    def get_method(self, name=None):
        """Return an optimization method.

        Parameters
        ----------
        name : str, optional
           If not given, the current method is returned, otherwise
           it should be one of the names returned by the
           `list_methods` function.

        Returns
        -------
        method : object
           An object representing the optimization method.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `name` argument is not recognized.

        See Also
        --------
        get_method_opt : Get the options for the current optimization method.
        list_methods : List the supported optimization methods.
        set_method : Change the optimization method.
        set_method_opt : Change an option of the current optimization method.

        Examples
        --------

        The fields of the object returned by `get_method` can be
        used to view or change the method options.

        >>> method = ui.get_method()
        >>> print(method.name)
        levmar
        >>> print(method)
        name    = levmar
        ftol    = 1.19209289551e-07
        xtol    = 1.19209289551e-07
        gtol    = 1.19209289551e-07
        maxfev  = None
        epsfcn  = 1.19209289551e-07
        factor  = 100.0
        verbose = 0
        >>> method.maxfev = 5000

        """
        if name is None:
            return self._current_method
        _check_type(name, basestring, 'name', 'a string')
        return self._get_method_by_name(name)

    ### Ahelp ingest: 2015-04-27 DJB
    ### DOC-TODO: is this guaranteed to be the same as get_method().name
    ###           or get_method().name.lower() and, if so, shouldn't this be
    ###           how it is coded?
    def get_method_name(self):
        """Return the name of current Sherpa optimization method.

        Returns
        -------
        name : str
           The name of the current optimization method, in lower
           case. This may not match the value sent to `set_method`
           because some methods can be set by multiple names.

        See Also
        --------
        get_method : Return an optimization method.
        get_method_opt : Get the options for the current optimization method.

        Examples
        --------

        >>> get_method_name()
        'levmar'

        The 'neldermead' method can also be referred to as 'simplex':

        >>> set_method('simplex')
        >>> get_method_name()
        'neldermead'

        """
        return type(self.get_method()).__name__.lower()

    ### Ahelp ingest: 2015-04-25 DJB
    ### DOC-TODO: remove the list of supported methods once the
    ### relevant documenation has been updated.
    def set_method(self, meth):
        """Set the optimization method.

        The primary task of Sherpa is to fit a model M(p) to a set of
        observed data, where the vector p denotes the model
        parameters. An optimization method is one that is used to
        determine the vector of model parameter values, p0, for which
        the chosen fit statistic is minimized.

        Parameters
        ----------
        meth : str
           The name of the method (case is not important). The
           `list_methods` function returns the list of supported values.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `meth` argument is not recognized.

        See Also
        --------
        get_method_name : Return the name of the current optimization method.
        list_methods : List the supported optimization methods.
        set_stat : Set the fit statistic.

        Notes
        -----
        The available methods include:

        `levmar`
           The Levenberg-Marquardt method is an interface to the
           MINPACK subroutine lmdif to find the local minimum of
           nonlinear least squares functions of several variables by a
           modification of the Levenberg-Marquardt algorithm [1]_.

        `moncar`
           The implementation of the moncar method is based on [2]_.

        `neldermead`
           The implementation of the Nelder Mead Simplex direct search
           is based on [3]_.

        `simplex`
           This is another name for `neldermead`.

        References
        ----------

        .. [1] J.J. More, "The Levenberg Marquardt algorithm:
           implementation and theory," in Lecture Notes in Mathematics
           630: Numerical Analysis, G.A. Watson (Ed.), Springer-Verlag:
           Berlin, 1978, pp.105-116. 

        .. [2] Storn, R. and Price, K. "Differential Evolution: A
           Simple and Efficient Adaptive Scheme for Global Optimization
           over Continuous Spaces." J. Global Optimization 11, 341-359,
           1997. 

        .. [3] Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
           Paul E. Wright "Convergence Properties of the Nelder-Mead
           Simplex Algorithm in Low Dimensions", SIAM Journal on
           Optimization,Vol. 9, No. 1 (1998), pages 112-147.

        Examples
        --------

        >>> set_method('neldermead')

        """
        if isinstance(meth, basestring):
            meth = self._get_method_by_name(meth)
        else:
            _check_type(meth, sherpa.optmethods.OptMethod, 'meth',
                        'a method name or object')
        self._current_method = meth

    def _check_method_opt(self, optname):
        _check_type(optname, basestring, 'optname', 'a string')
        if optname not in self._current_method.config:
            raise ArgumentErr('badopt', optname, self.get_method_name())

    ### Ahelp ingest: 2015-04-27 DJB
    def get_method_opt(self, optname=None):
        """Return one or all of the options for the current optimization
        method.

        This is a helper function since the optimization options can also
        be read directly using the object returned by `get_method`.

        Parameters
        ----------
        optname : str, optional
           If not given, a dictionary of all the options are returned.
           When given, the individual value is returned.
            
        Returns
        -------
        value : dictionary or value
           
        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `optname` argument is not recognized.

        See Also
        --------
        get_method : Return an optimization method.
        set_method : Change the optimization method.
        set_method_opt : Change an option of the current optimization method.

        Examples
        --------

        >>> get_method_opt('maxfev') == None
        True

        >>> mopts = get_method_opt()
        >>> mopts['maxfev'] == None
        True

        """
        if optname is None:
            return self._current_method.config
        self._check_method_opt(optname)
        return self._current_method.config[optname]

    ### Ahelp ingest: 2015-04-27 DJB
    def set_method_opt(self, optname, val):
        """Set an option for the current optimization method.

        This is a helper function since the optimization options can also
        be set directly using the object returned by `get_method`.

        Parameters
        ----------
        optname : str
           The name of the option to set. The `get_method`
           and `get_method_opt` routines can be used to find
           out valid values for this argument.
        val :
           The new value for the option.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `optname` argument is not recognized.

        See Also
        --------
        get_method : Return an optimization method.
        get_method_opt : Return one or all options of the current optimization method.
        set_method : Change the optimization method.

        Examples
        --------

        Change the maxfev parameter for the current optimizer
        to 2000.

        >>> set_method_opt('maxfev', 2000)

        """
        self._check_method_opt(optname)
        self._current_method.config[optname] = val

    #### Iterative Fitting Methods for CIAO 4.3 testing

    ### Ahelp ingest: 2015-05-04 DJB
    def get_iter_method_name(self):
        """Return the name of the iterative fitting scheme.

        Returns
        -------
        name : str
           The name of the iterative fitting scheme set by
           `set_iter_method`.

        See Also
        --------
        list_iter_methods : List the iterative fitting schemes.
        set_iter_method : Set the iterative-fitting scheme used in the fit.

        """
        return self._current_itermethod['name']

    ### Ahelp ingest: 2015-05-04 DJB
    def get_iter_method_opt(self, optname=None):
        """Return one or all options for the iterative-fitting scheme.

        Parameters
        ----------
        optname : str, optional
           If not given, a dictionary of all the options are returned.
           When given, the individual value is returned.

        Returns
        -------
        value : dictionary or value
           The dictionary is empty when no iterative scheme is being
           used.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `optname` argument is not recognized.

        See Also
        --------
        get_iter_method_name : Return the name of the iterative fitting scheme.
        set_iter_method_opt : Set an option for the iterative-fitting scheme.
        set_iter_method : Set the iterative-fitting scheme used in the fit.

        """
        itermethod_opts = dict(self._current_itermethod)
        del itermethod_opts['name']
        
        if optname is None:
            return itermethod_opts

        _check_type(optname, basestring, 'optname', 'a string')
        if optname not in itermethod_opts:
            raise ArgumentErr('badopt', optname, self._current_itermethod['name'])
        return itermethod_opts[optname]

    ### Ahelp ingest: 2015-05-04 DJB
    def list_iter_methods(self):
        """List the iterative fitting schemes.

        Returns
        -------
        schemes : list of str
           A list of the names that can be used with
           `set_iter_method`.

        See Also
        --------
        get_iter_method_name : Return the name of the iterative fitting scheme.
        set_iter_method : Set the iterative-fitting scheme used in the fit.

        Examples
        --------

        >>> list_iter_methods()
        ['none', 'primini', 'sigmarej']

        """
        keys = self._itermethods.keys()[:]
        keys.sort()
        return keys

    ### Ahelp ingest: 2015-05-04 DJB
    ### DOC-TODO: this information is also in sherpa/fit.py
    ### DOC-TODO: this raises a ValueError rather than a Sherpa error class
    def set_iter_method(self, meth):
        """Set the iterative-fitting scheme used in the fit.

        Control whether an iterative scheme should be applied to
        the fit.

        Parameters
        ----------
        meth : { 'none', 'primini', 'sigmarej' }
           The name of the scheme used during the fit; 'none' means no
           scheme is used. It is only valid to change the scheme
           when a chi-square statistic is in use.

        Raises
        ------
        TypeError
           When the `meth` argument is not recognized.

        See Also
        --------
        fit : Fit a model to one or more data sets.
        get_iter_method_name : Return the name of the iterative fitting scheme.
        get_iter_method_opt : Return one or all options for the iterative-fitting scheme.
        list_iter_methods : List the iterative fitting schemes.
        set_iter_method_opt : Set an option for the iterative-fitting scheme.
        set_stat : Set the statistical method.

        Notes
        -----
        The parameters of each scheme are described in
        `set_iter_method_opt`.

        The `primini` scheme is used for re-calculating statistical
        errors, using the best-fit model parameters from the
        *previous* fit, until the fit can no longer be improved.

        This is a chi-square statistic where the variance is computed
        from model amplitudes derived in the previous iteration of the
        fit. This 'Iterative Weighting' ([1]_) attempts to remove
        biased estimates of model parameters which is inherent in
        chi-square2 statistics ([2]_).

        The variance in bin i is estimated to be:

        sigma^2_i^j = S(i, t_s^(j-1)) + (A_s/A_b)^2 B_off(i, t_b^(j-1))

        where j is the number of iterations that have been carried out
        in the fitting process, B_off is the background model
        amplitude in bin i of the off-source region, and t_s^(j-1) and
        t_b^(j-1) are the set of source and background model parameter
        values derived during the iteration previous to the current
        one. The variances are set to an array of ones on the first
        iteration.

        In addition to reducing parameter estimate bias, this
        statistic can be used even when the number of counts in each
        bin is small (< 5), although the user should proceed with
        caution.

        The `sigmarej` scheme is based on the IRAF `sfit` function
        [3]_, where after a fit data points are excluded if the value
        of `(data-model)/error)` exceeds a threshold, and the data
        re-fit. This removal of data points continues until the fit
        has converged. The error removal can be asymmetric, since
        there are separate parameters for the lower and upper limits.

        References
        ----------

        .. [1] "Multiparameter linear least-squares fitting to Poisson
               data one count at a time", Wheaton et al. 1995, ApJ 438,
               322
               http://adsabs.harvard.edu/abs/1995ApJ...438..322W

        .. [2] "Bias-Free Parameter Estimation with Few Counts, by
               Iterative Chi-Squared Minimization", Kearns, Primini, &
               Alexander, 1995, ADASS IV, 331
               http://adsabs.harvard.edu/abs/1995ASPC...77..331K

        .. [3] http://iraf.net/irafhelp.php?val=sfit

        """
        if isinstance(meth, basestring):
            if (self._itermethods.has_key(meth) == True):
                self._current_itermethod = self._itermethods[meth]
            else:
                raise TypeError(meth + ' is not an iterative fitting method')
        else:
            _argument_type_error(meth, 'a string')

    ### Ahelp ingest: 2015-05-04 DJB
    def set_iter_method_opt(self, optname, val):
        """Set an option for the iterative-fitting scheme.

        Parameters
        ----------
        optname : str
           The name of the option to set. The
           `get_iter_method_opt_method` routine can be used to find
           out valid values for this argument.
        val :
           The new value for the option.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `optname` argument is not recognized.

        See Also
        --------
        get_iter_method_name : Return the name of the iterative fitting scheme.
        get_iter_method_opt : Return one or all options for the iterative-fitting scheme.
        list_iter_methods : List the iterative fitting schemes.
        set_iter_method : Set the iterative-fitting scheme used in the fit.

        Notes
        -----
        The supported fields for the `primini` scheme are:

        `maxiters`
           The maximum number of iterations to perform.

        `tol`
           The iteration stops when the change in the best-fit
           statistic varies by less than this value.

        The supported fields for the `sigmarej` scheme are:

        `grow`
           The number of points adjacent to a rejected point that
           should also be removed. A value of `0` means that only the
           discrepant point is removed whereas a value of `1` means
           that the two adjacent points (one lower and one higher)
           will also be removed.

        `hrej`
           The rejection criterion in units of sigma, for data
           points above the model (`hrej` is >= 0).

        `lrej`
           The rejection criterion in units of sigma, for data
           points below the model (`lrej` is >= 0).

        `maxiters`
           The maximum number of iterations to perform. If this
           value is `0` then the fit will run until it has
           converged.

        Examples
        --------

        Reject any points that are more than 5 sigma away from the
        best fit model and re-fit.

        >>> set_iter_method('sigmarej')
        >>> set_iter_method_opt('lrej', 5)
        >>> set_iter_method_opt('hrej', 5)
        >>> fit()

        """
        _check_type(optname, basestring, 'optname', 'a string')
        if (optname not in self._current_itermethod or
            optname == 'name'):
            raise ArgumentErr('badopt', optname, self._current_itermethod['name'])
        self._current_itermethod[optname] = val

    ###########################################################################
    # Statistics
    ###########################################################################


    ### Ahelp ingest: 2015-04-27 DJB
    def list_stats(self):
        """List the fit statistics.

        Returns
        -------
        stat : list of str
           A list of the names that can be used with
           `set_stat`.

        See Also
        --------
        get_stat_name : Return the name of the current statistical method.
        set_stat : Set the statistical method.

        Examples
        --------

        >>> list_stats()
        ['cash',
         'chi2',
         'chi2constvar',
         'chi2datavar',
         'chi2gehrels',
         'chi2modvar',
         'chi2xspecvar',
         'cstat',
         'leastsq']

        """
        keys = self._stats.keys()[:]
        keys.sort()
        return keys

    def _get_stat_by_name(self, name):
        stat = self._stats.get(name.lower())
        if stat is None:
            raise ArgumentErr('badstat', name)
        return stat

    ### Ahelp ingest: 2015-04-27 DJB
    def get_stat(self, name=None):
        """Return a fit statisic.

        Parameters
        ----------
        name : str, optional
           If not given, the current fit statistic is returned, otherwise
           it should be one of the names returned by the
           `list_stats` function.

        Returns
        -------
        stat : object
           An object representing the fit statistic.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `name` argument is not recognized.

        See Also
        --------
        get_stat_name : Return the name of the current fit statistic.
        list_stats : List the fit statistics.
        set_stat : Change the fit statistic.

        Examples
        --------

        >>> stat = ui.stat()
        >>> stat
        Chi Squared with Gehrels variance
        >>> stat.name
        'chi2gehrels'

        """
        if name is None:
            return self._current_stat
        _check_type(name, basestring, 'name', 'a string')
        return self._get_stat_by_name(name)

    ### Ahelp ingest: 2015-04-27 DJB
    def get_stat_name(self):
        """Return the name of the current fit statistic.

        Returns
        -------
        name : str
           The name of the current fit statistic method, in lower
           case.

        See Also
        --------
        get_stat : Return a fit statistic.
        set_stat : Set the fit statistic.

        Examples
        --------

        >>> get_method_name()
        'chi2gehrels'

        >>> set_stat('cash')
        >>> get_stat_name()
        'cash'

        """
        return type(self.get_stat()).__name__.lower()

    ### Ahelp ingest: 2015-04-24 DJB
    ### DOC-TODO: remove the list of supported methods once the
    ### relevant documenation has been updated.
    def set_stat(self, stat):
        """Set the statistical method.

        Changes the method used to evaluate the fit statistic, that is
        the numerical measure that determines how closely the model
        represents the data.

        Parameters
        ----------
        stat : str
           The name of the statistic (case is not important). The
           `list_stats` function returns the list of supported values.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `stat` argument is not recognized.

        See Also
        --------
        calc_stat : Calculate the statistic value for a dataset.
        get_stat_name : Return the current statistic method.
        list_stats : List the supported fit statistics.

        Notes
        -----
        The available statistics include:

        `cash`
           A maximum likelihood function [1]_.

        `chi2`
           \chi^2 statistic using the supplied error values.

        `chi2constvar`
           \chi^2 with constant variance computed from the counts
           data.

        `chi2datavar`
           \chi^2 with data variance.

        `chi2gehrels`
           \chi^2 with gehrels method [2]_. This is the default method.

        `chi2modvar`
           \chi^2 with model amplitude variance.

        `chi2xspecvar`
           \chi^2 with data variance (XSPEC-style,
           variance = 1.0 if data less than or equal to 0.0).

        `cstat`
           A maximum likelihood function
           (the XSPEC implementation of the Cash function) [3]_.

        `leastsq`
           The least-squares statisic (the error is not used in
           this statistic).

        References
        ----------

        .. [1] Cash, W. "Parameter estimation in astronomy through
        application of the likelihood ratio", ApJ, vol 228,
        p. 939-947 (1979).
        http://adsabs.harvard.edu/abs/1979ApJ...228..939C

        .. [2] Gehrels, N. "Confidence limits for small numbers of
        events in astrophysical data", ApJ, vol 303,
        p. 336-346 (1986).
        http://adsabs.harvard.edu/abs/1986ApJ...303..336G

        .. [3] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html

        Examples
        --------

        >>> set_stat('cash')

        """
        if isinstance(stat, basestring):
            stat = self._get_stat_by_name(stat)
        else:
            _check_type(stat, sherpa.stats.Stat, 'stat',
                        'a statistic name or object')

        self._current_stat = stat


    ###########################################################################
    # Data sets
    ###########################################################################


    ### Ahelp ingest: 2015-04-29 DJB
    def list_data_ids(self):
        """List the identifiers for the loaded data sets.

        Returns
        -------
        ids : list of int or str
           A list of the data set identifiers that have been created
           by commands like `load_data` and `load_arrays`.

        See Also
        --------
        delete_data : Delete a data set by identifier.
        load_arrays : Create a data set from arrays of data.
        load_data : Create a data set from a file.

        Examples
        --------

        In this case only one data set has been loaded:

        >>> list_data_ids()
        [1]

        Two data sets have been loaded, using string identifiers:

        >>> list_data_ids()
        ['nucleus', 'jet']

        """
        keys = self._data.keys()[:]
        keys.sort()
        return keys

    ### Ahelp ingest: 2015-05-01 DJB
    def get_data(self, id=None):
        """Return the data set by identifier.

        The object returned by the call can be used to query and
        change properties of the data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.

        Returns
        -------
        data :
           An instance of a sherpa.Data.Data-derived class.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If no model expression has been set for the data set
           (with `set_model` or `set_source`).

        See Also
        --------
        copy_data : Copy a data set to a new identifier.
        delete_data : Delete a data set by identifier.
        load_data : Create a data set from a file.
        set_data : Set a data set.

        Examples
        --------

        >>> d = get_data()

        >>> dimg = get_data('img')

        >>> load_arrays('tst', [10, 20, 30], [5.4, 2.3, 9.8])
        >>> print(get_data('tst'))
        name      =
        x         = Int64[3]
        y         = Float64[3]
        staterror = None
        syserror  = None

        """
        return self._get_item(id, self._data, 'data set', 'has not been set')

    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: terrible synopsis
    def set_data(self, id, data=None):
        """Set a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        data : instance of a sherpa.Data.Data-derived class
           The new contents of the data set. This can be
           copied from an existing data set or loaded
           in from a file (e.g. `unpack_data`).

        See Also
        --------
        copy_data : Copy a data set to a new identifier.
        delete_data : Delete a data set by identifier.
        get_data : Return the data set by identifier.
        load_data : Create a data set from a file.
        unpack_data : Read in a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `data` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `data` parameters,
        respectively.

        Examples
        --------

        >>> d1 = get_data(2)
        >>> set_data(d1)

        Copy the background data from the default data set into
        a new data set identified as 'bkg':

        >>> set_data('bkg', get_bkg())

        """
        if data is None:
            id, data = data, id
        self._set_item(id, data, self._data, sherpa.data.Data, 'data',
                       'a data set')


    def _read_error(self, filename, *args, **kwargs):
        err = sherpa.io.get_ascii_data(filename, *args, **kwargs)[1].pop()
        return err


    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-06 DJB
    ### DOC-NOTE: is ncols really 2 here? Does it make sense?
    def load_staterror(self, id, filename=None, ncols=2, *args, **kwargs):
        """Load the statistical errors from a file.

        Read in a column or image from a file and use the values
        as the statistical errors for a data set. This over rides
        the errors calculated by any statistic, such as
        `chi2gehrels` or `chi2datavar`.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to read in. Supported formats depends
           on the I/O library in use (Crates or AstroPy) and the
           type of data set (e.g. 1D or 2D).
        ncols : int, optional
           The number of columns to read in (the first `ncols` columns
           in the file).
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.

        See Also
        --------
        get_staterror : Return the statistical error on the dependent axis of a data set.
        load_syserror : Load the systematic errors from a file.
        set_staterror : Set the statistical errors on the dependent axis of a data set.
        set_stat : Set the statistical method.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Read in the first column from 'tbl.dat':

        >>> load_staterror('tbl.dat')

        Use the column labelled 'col3'

        >>> load_staterror('tbl.dat', colkeys=['col3'])

        When using the Crates I/O library, the file name can include
        CIAO Data Model syntax, such as column selection:

        >>> load_staterror('tbl.dat[cols col3]')

        Read in the first column from the file 'errors.fits' as the
        statistical errors for the 'core' data set:

        >>> load_staterror('core', 'errors.fits')

        """
        if filename is None:
            id, filename = filename, id

        self.set_staterror(id, self._read_error(filename, ncols=ncols,
                                                *args, **kwargs))


    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-06 DJB
    ### DOC-NOTE: is ncols really 2 here? Does it make sense?
    def load_syserror(self, id, filename=None, ncols=2, *args, **kwargs):
        """Load the systematic errors from a file.

        Read in a column or image from a file and use the values
        as the systematic errors for a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to read in. Supported formats depends
           on the I/O library in use (Crates or AstroPy) and the
           type of data set (e.g. 1D or 2D).
        ncols : int, optional
           The number of columns to read in (the first `ncols` columns
           in the file).
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.

        See Also
        --------
        get_syserror : Return the systematic error on the dependent axis of a data set.
        load_staterror : Load the statistical errors from a file.
        set_syserror : Set the systematic errors on the dependent axis of a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Read in the first column from 'tbl.dat':

        >>> load_syserror('tbl.dat')

        Use the column labelled 'col3'

        >>> load_syserror('tbl.dat', colkeys=['col3'])

        When using the Crates I/O library, the file name can include
        CIAO Data Model syntax, such as column selection:

        >>> load_syserror('tbl.dat[cols col3]')

        Read in the first column from the file 'errors.fits' as the
        systematic errors for the 'core' data set:

        >>> load_syserror('core', 'errors.fits')

        """
        if filename is None:
            id, filename = filename, id

        self.set_syserror(id, self._read_error(filename, ncols=ncols,
                                               *args, **kwargs))


    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-06 DJB
    ### DOC-TODO: does ncols make sense here? (have removed for now)
    ### DOC-TODO: labelling as AstroPy; i.e. assuming conversion
    ###           from PyFITS lands soon.
    def load_filter(self, id, filename=None, ignore=False, ncols=2,
                    *args, **kwargs):
        """Load the filter array from a file and add to a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file that contains the filter
           information. This file can be a FITS table or an ASCII
           file. Selection of the relevant column depends on the I/O
           library in use (Crates or AstroPy).
        ignore : bool, optional
           If `False` (the default) then include bins with a non-zero
           filter value, otherwise exclude these bins.
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.

        See Also
        --------
        get_filter : Return the filter expression for a data set.
        ignore : Exclude data from the fit.
        notice : Include data in the fit.
        save_filter : Save the filter array to a file.
        set_filter : Set the filter array of a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Read in the first column of the file and apply it to the
        default data set:

        >>> load_filter('filt.dat')

        Select the `FILTER` column of the file:

        >>> load_filter(2, 'filt.dat', colkeys=['FILTER'])

        When using Crates as the I/O library, the above can
        also be written as

        >>> load_filter(2, 'filt.dat[cols filter]')

        """
        if filename is None:
            id, filename = filename, id

        self.set_filter(id,
            self._read_error(filename, ncols=ncols, *args, **kwargs),
                        ignore=ignore)


    ### Ahelp ingest: 2015-05-06 DJB
    def set_filter(self, id, val=None, ignore=False):
        """Set the filter array of a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        val : array
           The array of filter values (`0` or `1`). The size should
           match the array returned by `get_dep`.
        ignore : bool, optional
           If `False` (the default) then include bins with a non-zero
           filter value, otherwise exclude these bins.

        See Also
        --------
        get_dep : Return the dependent axis of a data set.
        get_filter : Return the filter expression for a data set.
        ignore : Exclude data from the fit.
        load_filter : Load the filter array from a file and add to a data set.
        notice : Include data in the fit.
        save_filter : Save the filter array to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `val` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `val` parameters,
        respectively.

        Examples
        --------

        Ignore those bins with a value less 20.

        >>> d = get_dep()
        >>> f = d >= 20
        >>> set_filter(f)

        """
        if val is None:
            val, id = id, val
        filter = numpy.asarray(val, dtype=numpy.bool_)
        d = self.get_data(id)
        if numpy.iterable(d.mask):
            if len(d.mask) == len(filter):
                if not ignore:
                    d.mask |= filter
                else:
                    d.mask &= ~filter
            else:
                raise sherpa.utils.err.DataErr('mismatch',
                                               len(d.mask), len(filter))
        else:
            if len(d.get_y(False)) == len(filter):
                if not ignore:
                    d.mask = filter
                else:
                    d.mask = ~filter
            else:
                raise sherpa.utils.err.DataErr('mismatch',
                                               len(d.get_y(False)), len(filter))


    # also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-11 DJB
    def set_dep(self, id, val=None):
        """Set the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        val : array
           The array of values for the dependent axis.

        See Also
        --------
        dataspace1d : Create the independent axis for a 1D data set.
        dataspace2d : Create the independent axis for a 2D data set.
        get_dep : Return the dependent axis of a data set.
        load_arrays : Create a data set from array values.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `val` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `val` parameters,
        respectively.

        Examples
        --------

        Create a 1D data set with values at (0,4), (2,10), (4,12),
        (6,8), (8,2), and (10,12):

        >>> dataspace1d(0, 10, 2, dstype=Data1D)
        >>> set_dep([4, 10, 12, 8, 2, 12])

        Set the values for the data set 'src':

        >>> set_dep('src', y1)

        """
        if val is None:
            val, id = id, val
        d = self.get_data(id)
        dep=None
        if numpy.iterable(val):
            dep = numpy.asarray(val, SherpaFloat)
        else:
            val = SherpaFloat(val)
            dep = numpy.array([val]*len(d.get_indep()[0]))
        d.y = dep


    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def set_staterror(self, id, val=None, fractional=False):
        """Set the statistical errors on the dependent axis of a data set.

        These values over-ride the errors calculated by any statistic,
        such as `chi2gehrels` or `chi2datavar`.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        val : array or scalar
           The systematic error.
        fractional : bool, optional
           If `False` (the default value), then the `val` parameter is
           the absolute value, otherwise the `val` parameter
           represents the fractional error, so the absolute value is
           calculated as `get_dep() * val` (and `val` must be
           a scalar).

        See Also
        --------
        load_staterror : Set the statistical errors on the dependent axis of a data set.
        load_syserror : Set the systematic errors on the dependent axis of a data set.
        set_syserror : Set the systematic errors on the dependent axis of a data set.
        get_error : Return the errors on the dependent axis of a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `val` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `val` parameters,
        respectively.

        Examples
        --------

        Set the statistical error for the default data set to the value
        in `dys` (a scalar or an array):

        >>> set_staterror(dys)

        Set the statistical error on the `core` data set to be 5% of
        the data values:

        >>> set_staterror('core', 0.05, fractional=True)

        """
        if val is None:
            val, id = id, val
        err=None
        d = self.get_data(id)
        fractional=sherpa.utils.bool_cast(fractional)
        if numpy.iterable(val):
            err = numpy.asarray(val, SherpaFloat)
        elif val is not None:
            val = SherpaFloat(val)
            if fractional:
                err = val*d.get_dep()
            else:
                err = numpy.array([val]*len(d.get_dep()))
        d.staterror = err


    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def set_syserror(self, id, val=None, fractional=False):
        """Set the systematic errors on the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        val : array or scalar
           The systematic error.
        fractional : bool, optional
           If `False` (the default value), then the `val` parameter is
           the absolute value, otherwise the `val` parameter
           represents the fractional error, so the absolute value is
           calculated as `get_dep() * val` (and `val` must be
           a scalar).

        See Also
        --------
        load_staterror : Set the statistical errors on the dependent axis of a data set.
        load_syserror : Set the systematic errors on the dependent axis of a data set.
        set_staterror : Set the statistical errors on the dependent axis of a data set.
        get_error : Return the errors on the dependent axis of a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `val` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `val` parameters,
        respectively.

        Examples
        --------

        Set the systematic error for the default data set to the value
        in `dys` (a scalar or an array):

        >>> set_syserror(dys)

        Set the systematic error on the `core` data set to be 5% of
        the data values:

        >>> set_syserror('core', 0.05, fractional=True)

        """
        if val is None:
            val, id = id, val
        err=None
        d = self.get_data(id)
        fractional=sherpa.utils.bool_cast(fractional)        
        if numpy.iterable(val):
            err = numpy.asarray(val, SherpaFloat)
        elif val is not None:
            val = SherpaFloat(val)
            if fractional:
                err = val*d.get_dep()
            else:
                err = numpy.array([val]*len(d.get_dep()))
        d.syserror = err

    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def get_staterror(self, id=None, filter=False):
        """Return the statistical error on the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.

        Returns
        -------
        axis : array
           The statistical error for each data point. This may be
           estimated from the data (e.g. with the `chi2gehrels`
           statistic) or have been set explicitly (`set_staterror`).

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        get_syserror : Return the systematic errors on the dependent axis of a data set.
        list_data_ids : List the identifiers for the loaded data sets.
        set_staterror : Set the statistical errors on the dependent axis of a data set.

        Examples
        --------

        If not explicitly given, the statistical errors on a data set
        may be calculated from the data values (the independent axis),
        depending on the chosen statistic:

        >>> load_arrays(1, [10,15,19], [4,5,9])
        >>> set_stat('chi2datavar')
        >>> get_staterror()
        array([ 2.        ,  2.23606798,  3.        ])
        >>> set_stat('chi2gehrels')
        >>> get_staterror()
        array([ 3.17944947,  3.39791576,  4.122499  ])

        If the statistical errors are set - either when the data set
        is created or with a call to `set_errors` - then these values
        will be used, no matter the statistic:

        >>> load_arrays(1, [10,15,19], [4,5,9], [2,3,5])
        >>> set_stat('chi2datavar')
        >>> get_staterror()
        array([2, 3, 5])
        >>> set_stat('chi2gehrels')
        >>> get_staterror()
        array([2, 3, 5])

        """
        return self.get_data(id).get_staterror(filter,
                                               self.get_stat().calc_staterror)


    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def get_syserror(self, id=None, filter=False):
        """Return the systematic error on the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.

        Returns
        -------
        axis : array
           The systematic error for each data point.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set has no systematic errors.
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        get_staterror : Return the statistical errors on the dependent axis of a data set.
        list_data_ids : List the identifiers for the loaded data sets.
        set_syserror : Set the systematic errors on the dependent axis of a data set.

        """
        d = self.get_data(id)
        id = self._fix_id(id)
        err = d.get_syserror(filter)
        if err is None or not numpy.iterable(err):
            raise sherpa.utils.err.DataErr('nosyserr', id)
        return err


    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def get_error(self, id=None, filter=False):
        """Return the errors on the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.

        Returns
        -------
        axis : array
           The error for each data point, formed by adding the
           statistical and systematic errors in quadrature.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        get_staterror : Return the statistical errors on the dependent axis of a data set.
        get_syserror : Return the systematic errors on the dependent axis of a data set.
        list_data_ids : List the identifiers for the loaded data sets.

        """
        return self.get_data(id).get_error(filter,
                                           self.get_stat().calc_staterror)

    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-05 DJB
    ### DOC-NOTE: shouldn't this expose a filter parameter?
    def get_indep(self, id=None):
        """Return the independent axes of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        axis : tuple of arrays
           The independent axis values. These are the values at which
           the model is evaluated during fitting.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_dep : Return the dependent axis of a data set.
        list_data_ids : List the identifiers for the loaded data sets.

        Examples
        --------

        For a one-dimensional data set, the X values are returned:

        >>> load_arrays(1, [10,15,19], [4,5,9])
        >>> get_indep()
        (array([10, 15, 19]),)

        For a 2D data set the X0 and X1 values are returned:

        >>> load_arrays(2, [10,15,12,19], [12,14,10,17], [4,5,9,-2], Data2D)
        >>> get_indep(2)
        (array([10, 15, 12, 19]), array([12, 14, 10, 17]))

        """
        return self.get_data(id).get_indep()


    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def get_dep(self, id=None, filter=False):
        """Return the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.

        Returns
        -------
        axis : array
           The dependent axis values. The model estimate is compared
           to these values during fitting.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        list_data_ids : List the identifiers for the loaded data sets.

        Examples
        --------

        >>> load_arrays(1, [10,15,19], [4,5,9])
        >>> ignore_id(1, 17, None)
        >>> get_dep()
        array([4, 5, 9])
        >>> get_dep(filter=True)
        array([4, 5])

        >>> load_arrays(2, [10,15,12,19], [12,14,10,17], [4,5,9,-2], Data2D)
        >>> get_dep(2)
        array([ 4,  5,  9, -2])

        """
        return self.get_data(id).get_y(filter)


    ### Ahelp ingest: 2015-05-01 DJB
    def get_dims(self, id=None, filter=False):
        """Return the dimensions of the data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        filter : bool, optional
           If `True` then apply any filter to the data set before
           returning the dimensions. The default is `False`.

        Returns
        -------
        dims : a tuple of int

        See Also
        --------
        ignore : Exclude data from the fit.
        sherpa.astro.utils.ignore2d : Exclude a spatial region from an image.
        notice : Include data in the fit.
        sherpa.astro.utils.notice2d : Include a spatial region of an image.

        """
        return self.get_data(id).get_dims(filter)

    ### Ahelp ingest: 2015-05-06 DJB
    ### DOC-NOTE: should there be a version in sherpa.astro.utils with a bkg_id
    ###           parameter?
    def get_filter(self, id=None):
        """Return the filter expression for a data set.

        This returns the filter expression, created by one or more
        calls to `ignore` and `notice`, for the data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        filter : str
           The empty string or a string expression representing the
           filter used. For PHA data dets the units are controlled by
           the analysis setting for the data set.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not exist.

        See Also
        --------
        ignore : Exclude data from the fit.
        load_filter : Load the filter array from a file and add to a data set.
        notice : Include data in the fit.
        save_filter : Save the filter array to a file.
        show_filter : Show any filters applied to a data set.
        set_filter : Set the filter array of a data set.

        Examples
        --------

        The default filter is the full dataset, given in the
        format `lowval:hival` (both are includive limits):

        >>> load_arrays(1, [10,15,20,25], [5,7,4,2])
        >>> get_filter()
        '10.0000:25.0000'

        The `notice` call restricts the data to the range between
        14 and 30. The resulting filter is the combination of this
        range and the data:

        >>> notice(14,30)
        >>> get_filter()
        '15.0000:25.0000'

        Ignoring the point at `x=20` means that only the points at
        `x=15` and `x=25` remain, so a comma-separated list is used:

        >>> ignore(19,22)
        >>> get_filter()
        '15.0000,25.0000'

        The filter equivalent to the per-bin array of filter values:

        >>> set_filter([1,1,0,1])
        >>> get_filter()
        '10.0000,15.0000,25.0000'

        """
        return self.get_data(id).get_filter()


    ### Ahelp ingest: 2015-05-01 DJB
    def copy_data(self, fromid, toid):
        """Copy a data set to a new identifier.

        This is a "deep" copy, so that once it has been
        made, changes to one of the data sets will not
        be made to the other data set.

        Parameters
        ----------
        fromid : int or str
           The input data set.
        toid : int or str
           The output data set.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If there is no data set with a `fromid` identifier.

        Examples
        --------

        >>> copy_data(1, 2)

        >>> copy_data(2, "orig")

        """
        data = self.get_data(fromid)
        data = copy.deepcopy(data)
        self.set_data(toid, data)

    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: this does not delete the source expression;
    ###           is this intended or a bug?
    def delete_data(self, id=None):
        """Delete a data set by identifier.

        The data set, and any associated structures - such as the ARF
        and RMF for PHA data sets - are removed.

        Parameters
        ----------
        id : int or str, optional
           The data set to delete. If not given then the default
           identifier is used, as returned by `get_default_id`.

        See Also
        --------
        clean : Clear all stored session data.
        copy_data : Copy a data set to a new identifier.
        delete_model : Delete the model expression from a data set.
        get_default_id : Return the default data set identifier.
        list_data_ids : List the identifiers for the loaded data sets.

        Notes
        -----
        The source expression is not removed by this function.

        Examples
        --------

        Delete the data from the default data set:

        >>> delete_data()

        Delete the data set identified as 'src':

        >>> delete_data('src')

        """
        id = self._fix_id(id)
        self._data.pop(id, None)


    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-11 DJB
    def dataspace1d(self, start, stop, step=1, numbins=None, 
                    id=None, dstype=sherpa.data.Data1DInt):
        """Create the independent axis for a 1D data set.

        Create an "empty" one-dimensional data set by defining the
        grid on which the points are defined (the independent axis).
        The values are set to 0.

        Parameters
        ----------
        start : number
           The minimum value of the axis.
        stop : number
           The maximum value of the axis.
        step : number, optional
           The separation between each grid point. This is not used if
           `numbins` is set.
        numbins : int, optional
           The number of grid points. This over-rides the `step`
           setting.
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        dstype : data class to use, optional
           What type of data is to be used. Supported values include
           `Data1DInt` (the default) and `Data1D`.

        See Also
        --------
        dataspace2d : Create the independent axis for a 2D data set.
        get_dep : Return the dependent axis of a data set.
        get_indep : Return the independent axes of a data set.
        set_dep : Set the dependent axis of a data set.

        Notes
        -----
        The meaning of the `stop` parameter depends on whether it is a
        binned or unbinned data set (as set by the `dstype`
        parameter).

        Examples
        --------

        Create a binned data set, starting at 1 and with a
        bin-width of 1.

        >>> dataspace1d(1, 5, 1)
        >>> print(get_indep())
        (array([ 1.,  2.,  3.,  4.]), array([ 2.,  3.,  4.,  5.]))

        This time for an un-binned data set:

        >>> dataspace1d(1, 5, 1, dstype=Data1D)
        >>> print(get_indep())
        (array([ 1.,  2.,  3.,  4.,  5.]),)

        Specify the number of bins rather than the grid spacing:

        >>> dataspace1d(1, 5, numbins=5, id=2)
        >>> (xlo, xhi) = get_indep(2)
        >>> xlo
        array([ 1. ,  1.8,  2.6,  3.4,  4.2])
        >>> xhi
        array([ 1.8,  2.6,  3.4,  4.2,  5. ])

        >>> dataspace1d(1, 5, numbins=5, id=3, dstype=Data1D)
        >>> (x, ) = get_indep(3)
        >>> x
        array([ 1.,  2.,  3.,  4.,  5.])

        """
        # support non-integrated grids with inclusive boundaries
        if dstype in (sherpa.data.Data1D,):
            stop += step

        xlo,xhi,y = sherpa.utils.dataspace1d(start, stop, step=step,
                                             numbins=numbins)
        args = [xlo,xhi,y]

        if dstype is not sherpa.data.Data1DInt:
            args = [xlo, y]

        self.set_data(id, dstype('dataspace1d', *args))


    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-11 DJB
    def dataspace2d(self, dims, id=None, dstype=sherpa.data.Data2D):
        """Create the independent axis for a 2D data set.

        Create an "empty" two-dimensional data set by defining the
        grid on which the points are defined (the independent axis).
        The values are set to 0.

        Parameters
        ----------
        dims : sequence of 2 number
           The dimensions of the grid in `(width,height)` order.
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        dstype : data class to use, optional
           What type of data is to be used. Supported values include
           `Data2D` (the default) and `Data2DInt`.

        See Also
        --------
        dataspace1d : Create the independent axis for a 1D data set.
        get_dep : Return the dependent axis of a data set.
        get_indep : Return the independent axes of a data set.
        set_dep : Set the dependent axis of a data set.

        Examples
        --------

        Create a 200 pixel by 150 pixel grid (number of columns by
        number of rows) and display it (each pixel has a value of 0):

        >>> dataspace2d([200,150])
        >>> image_data()

        Create a data space called "fakeimg":

        >>> dataspace2d([nx,ny], id="fakeimg")

        """
        x0, x1, y, shape = sherpa.utils.dataspace2d(dims)

        dataset = None
        if issubclass(dstype, sherpa.data.Data2DInt):
            dataset = dstype('dataspace2d', x0-0.5, x1-0.5, x0+0.5, x1+0.5,
                             y, shape)
        else:
            dataset = dstype('dataspace2d', x0, x1, y, shape)

        self.set_data(id, dataset)


    def fake(self, id=None, method=sherpa.utils.poisson_noise):
        """
        fake

        SYNOPSIS
           Fake source data using supplied noise function by data id

        SYNTAX

        Arguments:
           id      -  data id with corresponding model to evaluate and store
                      faked data.
                      default = default data id

           method  -  noise function
                      default = poisson_noise

        Returns:
           None

        DESCRIPTION
           Evalutes the source model by data id, adds noise defined by 'method'
           and stores the resulting array in the data set.  The data set can
           created from file or populated using dataspace1d or dataspace2d.

        EXAMPLES
           Blank data set with faked data.

              dataspace1d(0.1,10,0.1,dstype=Data1D)
              set_model(gauss1d.g1)
              fake()

        SEE ALSO
           dataspace1d, dataspace2d, set_model
        """
        data = self.get_data(id)
        model = self.get_model(id)
        self.set_dep(id, method(data.eval_model(model)))


    @staticmethod
    def _read_data(readfunc, filename, *args, **kwargs):
        _check_type(filename, basestring, 'filename', 'a string')
        return readfunc(filename, *args, **kwargs)

    # DOC-NOTE: also in sherpa.astro.utils
    ### DOC-TODO: What data types are supported here?
    ### Ahelp ingest: 2015-05-12 DJB
    def unpack_arrays(self, *args):
        """Create a sherpa data object from arrays of data.

        The object returned by `unpack_arrays` can be used in a
        `set_data` call.

        Parameters
        ----------
        a1, .., aN : array_like
           Arrays of data. The order, and number, is determined by
           the `dstype` parameter, and listed in the `load_arrays`
           routine.
        dstype :
           The data set type. The default is `Data1D` and values
           include: `Data1D`, `Data1DInt`, `Data2D`, and `Data2DInt`.

        Returns
        -------
        data
           The data set object matching the requested `dstype`.

        See Also
        --------
        get_data : Return the data set by identifier.
        load_arrays : Create a data set from array values.
        set_data : Set a data set.
        unpack_data : Create a sherpa data object from a file.

        Examples
        --------

        Create a 1D (unbinned) data set from the values in
        the x and y arrays. Use the returned object to create
        a data set labelled "oned":

        >>> x = [1,3,7,12]
        >>> y = [2.3,3.2,-5.4,12.1]
        >>> dat = unpack_arrays(x, y)
        >>> set_data("oned", dat)

        Include statistical errors on the data:

        >>> edat = unpack_arrays(x, y, dy)

        Create a "binned" 1D data set, giving the low,
        and high edges of the independent axis (xlo
        and xhi respectively) and the dependent values
        for this grid (y):

        >>> hdat = unpack_arrays(xlo, xhi, y, Data1DInt)

        """
        return sherpa.io.read_arrays(*args)


    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-12 DJB
    def unpack_data(self, filename, ncols=2, colkeys=None,
                    dstype=sherpa.data.Data1D, sep=' ', comment='#', require_floats=True):
        """Create a sherpa data object from a file.

        The object returned by `unpack_data` can be used in a
        `set_data` call.

        Parameters
        ----------
        filename : str
           The name of the file to read in. Supported formats depends
           on the I/O library in use (Crates or AstroPy) and the
           type of data set (e.g. 1D or 2D).
        ncols : int, optional
           The number of columns to read in (the first `ncols` columns
           in the file).
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        dstype : data class to use, optional
           What type of data is to be used. Supported values include
           `Data1D` (the default), `Data1DInt`, `Data2D`, and
           `Data2DInt`.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        require_floats : bool, optional
           If `True` (the default), non-numeric data values will
           raise a `ValueError`.

        Returns
        -------
        data
           The data set object.

        Raises
        ------
        ValueError
           If a column value can not be converted into a numeric value
           and the `require_floats` parameter is True.

        See Also
        --------
        get_data : Return the data set by identifier.
        load_arrays : Create a data set from array values.
        set_data : Set a data set.
        unpack_arrays : Create a sherpa data object from arrays of data.

        Examples
        --------

        Create a data object from the first two columns of the file
        "src.dat" and use it to create a Sherpa data set called "src":

        >>> dat = unpack_data('src.dat')
        >>> set_data('src', dat)

        Read in the first three columns - the independent axis (x),
        the dependent variable (y), and the error on y:

        >>> dat = unpack_data('src.dat', ncols=3)

        """
        return self._read_data(sherpa.io.read_data, filename, ncols, colkeys,
                               sep, dstype, comment, require_floats)


    # DOC-NOTE: also in sherpa.astro.utils
    def load_data(self, id, filename=None, ncols=2, colkeys=None,
                  dstype=sherpa.data.Data1D, sep=' ', comment='#', require_floats=True):
        """
        load_data

        SYNOPSIS
           Load a data text file by id

        SYNTAX

        Arguments:
           id         - data id

           filename   - filename with path

        Keyword Arguments:
           ncols      - number of columns
                        default = 2

           colkeys    - list of column names
                        default = None

           dstype     - dataset type desired
                        default = Data1D

           sep        - separation character between columns
                        default = ' '

           comment    - comment character
                        default = '#'

        Returns:
           None

        DESCRIPTION
           Load tabular data from column-based text file
           into a Sherpa dataset by data id.

        SEE ALSO
           load_arrays, unpack_data, unpack_arrays
        """
        if filename is None:
            id, filename = filename, id
        self.set_data(id, self.unpack_data(filename, ncols=ncols,
                                           colkeys=colkeys, dstype=dstype,
                                           sep=sep, comment=comment, require_floats=require_floats))

    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: rework the Data type notes section (also needed by unpack_arrays)
    ##@loggable(with_id = True)
    def load_arrays(self, id, *args):
        """Create a data set from array values.

        Parameters
        ----------
        id : int or str
           The identifier for the data set to use.
        *args :
           Two or more arrays, followed by the type of data set to
           create.

        See Also
        --------
        copy_data : Copy a data set to a new identifier.
        delete_data : Delete a data set by identifier.
        get_data : Return the data set by identifier.
        load_data : Create a data set from a file.
        set_data : Set a data set.
        unpack_arrays : Create a sherpa data object from arrays of data.

        Notes
        -----
        The data type identifier, which defaults to `Data1D`,
        determines the number, and order, of the required inputs.

        `Data1D`
           required fields: x, y
           optional fields: statistical error, systematic error

        `Data1DInt`
           required fields: xlo, xhi, y
           optional fields: statistical error, systematic error

        `Data2D`
           required fields: x0, x1, y
           optional fields: shape, statistical error, systematic error
           The `shape` argument should be a tuple giving the
           size of the data (ny,nx).

        `Data2DInt`
           required fields: x0lo, x1lo, x0hi, x1hi, y
           optional fields: shape, statistical error, systematic error
           The `shape` argument should be a tuple giving the
           size of the data (ny,nx).

        Examples
        --------

        Create a 1D data set with three points:

        >>> load_arrays(1, [10, 12, 15], [4.2, 12.1, 8.4])

        Create a 1D data set, with the identifier 'prof', from the
        arrays `x` (independent axis), `y` (dependent axis), and `dy`
        (statistical error on the dependent axis):

        >>> load_arrays('prof', x, y, dy)

        Explicitly define the type of the data set:

        >>> load_arrays('prof', x, y, dy, Data1D)

        Data set 1 is a histogram, where the bins cover the range
        1-3, 3-5, and 5-7 with values 4, 5, and 9 respectively.

        >>> load_arrays(1, [1,3,5], [3,5,7], [4,5,9], Data1DInt)

        """
        self.set_data(id, self.unpack_arrays(*args))


    def _save_type(self, objtype, id, filename, **kwargs):
        if filename is None:
            id, filename = filename, id
        d = self.get_data(id)

        args = None
        fields = None

        if type(d) in (sherpa.data.Data2D, sherpa.data.Data2DInt):
            backup = d.y
            if objtype == 'delchi':
                raise AttributeErr('badfunc', "save_delchi()", "images")

            imgtype = getattr(self,'get_' + objtype + '_image', None)
            if imgtype is None:
                raise AttributeErr('attributeerr', 'session', 'get_%s_image()' % objtype)
            obj = imgtype(id)
            args = [obj.y]
            fields = [str(objtype).upper()]

        else:
            plottype = getattr(self,'get_' + objtype + '_plot', None)
            if plottype is None:
                raise AttributeErr('attributeerr', 'session', 'get_%s_plot()' % objtype)
            obj = plottype(id)
            args = [obj.x, obj.y]
            fields = ["X", str(objtype).upper()]

        sherpa.io.write_arrays(filename, args, fields, **kwargs)


    # DOC-NOTE: also in sherpa.astro.utils with a different interface
    ### Ahelp ingest: 2015-05-12 DJB
    def save_arrays(self, filename, args, fields=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'):
        """Write a list of arrays to a file.

        Parameters
        ----------
        filename : str
           The name of the file to write the array to.
        args : array of arrays
           The arrays to write out.
        fields : array of str
           The column names (should match the size of `args`).
        clobber : bool, optional
           If `filename` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        linebreak : str, optional
           Indicate a new line. The default is '\n'.
        format : str, optional
           The format used to write out the numeric values. The
           default is '%g%.

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        save_data : Save the data to a file.
        save_image :
        save_table :

        Examples
        --------

        Write the x and y columns from the default data set to the
        file 'src.dat':

        >>> x = get_indep()
        >>> y = get_dep()
        >>> save_arrays('src.dat', [x,y])

        Use the column names "r" and "surbri" for the columns:

        >>> save_arrays('prof.txt', [x,y], fields=["r", "surbri"],
                        clobber=True)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        _check_type(filename, basestring, 'filename', 'a string')
        sherpa.io.write_arrays(filename, args, fields, sep, comment, clobber,
                               linebreak, format)

    def save_source(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'):
        """
        save_source

        SYNOPSIS
           Write the unconvolved source model to file

        SYNTAX

        Arguments:
           id         - data id
                        default = default data id

           filename   - filename with path

           clobber    - clobber the existing output file
                        default = False

           sep        - separation character between columns
                        default = ' '

           comment    - comment character
                        default = '#'

           linebreak  - line break character between rows
                        default = '\n'

           format     - array element format string
                        default = '%g'

        Returns:
           None

        DESCRIPTION
           Write the unconvolved source model to file.  NOTE that the source 
           model array written to file respects the filter.

        EXAMPLE

           save_source("source.dat")

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_model,
           save_resid, save_delchi
        """
        clobber=sherpa.utils.bool_cast(clobber)
        _check_type(filename, basestring, 'filename', 'a string')
        self._save_type('source', id, filename, clobber=clobber, sep=sep,
                        comment=comment, linebreak=linebreak, format=format)

    def save_model(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'):
        """
        save_model

        SYNOPSIS
           Write the convolved source model to file

        SYNTAX

        Arguments:
           id         - data id
                        default = default data id

           filename   - filename with path

           clobber    - clobber the existing output file
                        default = False

           sep        - separation character between columns
                        default = ' '

           comment    - comment character
                        default = '#'

           linebreak  - line break character between rows
                        default = '\n'

           format     - array element format string
                        default = '%g'

        Returns:
           None

        DESCRIPTION
           Write the convolved source model to file.  NOTE that the source 
           model array written to file respects the filter.

        EXAMPLE

           save_model("model.dat")

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_resid, save_delchi
        """
        clobber=sherpa.utils.bool_cast(clobber)
        _check_type(filename, basestring, 'filename', 'a string')
        self._save_type('model', id, filename, clobber=clobber, sep=sep,
                        comment=comment, linebreak=linebreak, format=format)

    # DOC-NOTE: also in sherpa.utils with a different interface
    ### Ahelp ingest: 2015-05-12 DJB
    def save_resid(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'):
        """Save the residuals (data-model) to a file.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to.
        clobber : bool, optional
           This flag controls whether an existing file can be
           overwritten (`True`) or if it raises an exception (`False`,
           the default setting).
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        linebreak : str, optional
           Indicate a new line. The default is '\n'.
        format : str, optional
           The format used to write out the numeric values. The
           default is '%g%.

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        save_chisqr :
        save_data : Save the data to a file.
        save_delchi : Save the ratio of residuals (data-model) to error to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns `X` and `RESID`. The
        residuals array respects any filter setting for the data set.

        Examples
        --------

        Write the residuals to the file "resid.dat":

        >>> save_resid('resid.dat')

        Write the residuals from the data set 'jet' to the
        file "resid.dat":

        >>> save_resid('jet', "resid.dat", clobber=True)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        _check_type(filename, basestring, 'filename', 'a string')
        self._save_type('resid', id, filename, clobber=clobber, sep=sep,
                        comment=comment, linebreak=linebreak, format=format)

    # DOC-NOTE: also in sherpa.utils with a different interface
    ### Ahelp ingest: 2015-05-12 DJB
    def save_delchi(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'):
        """Save the ratio of residuals (data-model) to error to a file.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to.
        clobber : bool, optional
           This flag controls whether an existing file can be
           overwritten (`True`) or if it raises an exception (`False`,
           the default setting).
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        linebreak : str, optional
           Indicate a new line. The default is '\n'.
        format : str, optional
           The format used to write out the numeric values. The
           default is '%g%.

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        save_chisqr :
        save_data : Save the data to a file.
        save_delchi : Save the ratio of residuals (data-model) to error to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns `X` and `DELCHI`. The
        residuals array respects any filter setting for the data set.

        Examples
        --------

        Write the residuals to the file "delchi.dat":

        >>> save_delchi('delchi.dat')

        Write the residuals from the data set 'jet' to the
        file "delchi.dat":

        >>> save_resid('jet', "delchi.dat", clobber=True)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        _check_type(filename, basestring, 'filename', 'a string')
        self._save_type('delchi', id, filename, clobber=clobber, sep=sep,
                        comment=comment, linebreak=linebreak, format=format)


    ### DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-12 DJB
    def save_data(self, id, filename=None, fields=None, sep=' ', comment='#',
                  clobber=False, linebreak='\n', format='%g'):
        """Save the data to a file.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to. The data is
           written out as an ASCII file.
        fields : array of str, optional
           The attributes of the data set to write out. If `None`,
           write out all the columns.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).
        linebreak : str, optional
           Indicate a new line. The default is '\n'.
        format : str, optional
           The format used to write out the numeric values. The
           default is '%g%.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If there is no matching data set.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        save_arrays : 
        save_delchi : Save the ratio of residuals (data-model) to error to a file.
        save_error : Save the errors to a file.
        save_filter : Save the filter array to a file.
        save_resid :
        save_staterror : Save the statistical errors to a file.
        save_syserror : Save the statistical errors to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Write the default data set out to the ASCII file 'src.dat':

        >>> save_data('src.dat')

        Only write out the x, y, and staterror columns for data set
        'rprof' to the file 'prof.out', over-writing it if it already
        exists:

        >>> save_data('rprof', 'prof.out', clobber=True,
                      fields=['x', 'y', 'staterror'])

        """
        clobber=sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        sherpa.io.write_data(filename, self.get_data(id), fields, sep,
                             comment, clobber, linebreak, format)


    # DOC-NOTE: also in sherpa.astro.utils with a different interface
    ### Ahelp ingest: 2015-05-06 DJB
    def save_filter(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'):
        """Save the filter array to a file.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        linebreak : str, optional
           Indicate a new line. The default is '\n'.
        format : str, optional
           The format used to write out the numeric values. The
           default is '%g%.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set has not been filtered.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        load_filter : Load the filter array from a file and add to a data set.
        save_data : Save the data to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns `X` and `FILTER`.

        Examples
        --------

        Write the filter from the default data set as an ASCII file:

        >>> save_filter('filt.dat')

        """
        clobber=sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        d = self.get_data(id)
        id = self._fix_id(id)
        if d.mask is False:
            raise sherpa.utils.err.DataErr('notmask')
        if not numpy.iterable(d.mask):
            raise sherpa.utils.err.DataErr('nomask', id)
        x = d.get_indep(filter=False)[0]
        mask = numpy.asarray(d.mask, numpy.int)
        self.save_arrays(filename, [x, mask], ['X', 'FILTER'],
                         clobber, sep, comment, linebreak, format)


    # DOC-NOTE: also in sherpa.astro.utils with a different interface
    ### Ahelp ingest: 2015-05-06 DJB
    def save_staterror(self, id, filename=None, clobber=False, sep=' ',
                       comment='#', linebreak='\n', format='%g'):
        """Save the statistical errors to a file.

        If the statistical errors have not been set explicitly, then
        the values calculated by the statistic - such as `chi2gehrels`
        or `chi2datavar` - will be used.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        linebreak : str, optional
           Indicate a new line. The default is '\n'.
        format : str, optional
           The format used to write out the numeric values. The
           default is '%g%.

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        load_staterror : Load the statistical errors from a file.
        save_error : Save the errors to a file.
        save_syserror : Save the systematic errors to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns `X` and `STAT_ERR`.

        Examples
        --------

        Write out the statistical errors from the default data set to the
        file 'errs.dat'.

        >>> save_staterror('errs.dat')

        Over-write the file it it already exists, and take the data
        from the data set "jet":

        >>> save_staterror('jet', 'err.out', clobber=True)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        x = self.get_data(id).get_indep(filter=False)[0]
        err = self.get_staterror(id, filter=False)
        self.save_arrays(filename, [x, err], ['X', 'STAT_ERR'],
                         clobber, sep, comment, linebreak, format)


    # DOC-NOTE: also in sherpa.astro.utils with a different interface
    ### Ahelp ingest: 2015-05-06 DJB
    def save_syserror(self, id, filename=None, clobber=False, sep=' ',
                       comment='#', linebreak='\n', format='%g'):
        """Save the statistical errors to a file.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        linebreak : str, optional
           Indicate a new line. The default is '\n'.
        format : str, optional
           The format used to write out the numeric values. The
           default is '%g%.

        Raises
        ------
        sherpa.utils.err.IOErr
           If the data set does not contain any systematic errors.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        load_syserror : Load the systematic errors from a file.
        save_error : Save the errors to a file.
        save_staterror : Save the statistical errors to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns `X` and `SYS_ERR`.

        Examples
        --------

        Write out the systematic errors from the default data set to the
        file 'errs.dat'.

        >>> save_syserror('errs.dat')

        Over-write the file it it already exists, and take the data
        from the data set "jet":

        >>> save_syserror('jet', 'err.out', clobber=True)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        x = self.get_data(id).get_indep(filter=False)[0]
        err = self.get_syserror(id, filter=False)
        self.save_arrays(filename, [x, err], ['X', 'SYS_ERR'],
                         clobber, sep, comment, linebreak, format)

    # DOC-NOTE: also in sherpa.astro.utils with a different interface
    ### Ahelp ingest: 2015-05-06 DJB
    def save_error(self, id, filename=None, clobber=False, sep=' ',
                       comment='#', linebreak='\n', format='%g'):
        """Save the errors to a file.

        The total errors for a data set are the quadrature combination
        of the statistical and systematic errors. The systematic
        errors can be 0. If the statistical errors have not been set
        explicitly, then the values calculated by the statistic - such
        as `chi2gehrels` or `chi2datavar` - will be used.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        linebreak : str, optional
           Indicate a new line. The default is '\n'.
        format : str, optional
           The format used to write out the numeric values. The
           default is '%g%.

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        load_staterror : Load the statistical errors from a file.
        load_syserror : Load the systematic errors from a file.
        save_data : Save the data to a file.
        save_staterror : Save the statistical errors to a file.
        save_syserror : Save the systematic errors to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns `X` and `ERR`.

        Examples
        --------

        Write out the errors from the default data set to the file
        'errs.dat'.

        >>> save_error('errs.dat')

        Over-write the file it it already exists, and take the data
        from the data set "jet":

        >>> save_error('jet', 'err.out', clobber=True)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        x = self.get_data(id).get_indep(filter=False)[0]
        err = self.get_error(id, filter=False)
        self.save_arrays(filename, [x, err], ['X', 'ERR'],
                         clobber, sep, comment, linebreak, format)


    def _notice_expr(self, expr=None, **kwargs):
        ids = self.list_data_ids()
        for vals in sherpa.utils.parse_expr(expr):
            self.notice_id(ids, *vals, **kwargs)


    def _notice_expr_id(self, ids, expr=None, **kwargs):
        for vals in sherpa.utils.parse_expr(expr):
            self.notice_id(ids, *vals, **kwargs)


    ### Ahelp ingest: 2015-05-06 DJB
    def notice(self, lo=None, hi=None, **kwargs):
        """Include data in the fit.

        Select one or more ranges of data to include by filtering on
        the independent axis value. The filter is applied to all data
        sets.

        Parameters
        ----------
        lo : number or str, optional
           The lower bound of the filter (when a number) or a string
           expression listing ranges in the form `a:b`, with multiple
           ranges allowed, where the ranges are separated by a
           `,`. The term `:b` means include everything up to, and
           including `b`, and `a:` means include everything that is
           higher than, or equal to, `a`.
        hi : number, optional
           The upper bound of the filter when `lo` is not a string.

        See Also
        --------
        notice_id : Include data for a data set.
        sherpa.astro.utils.notice2d : Include a spatial region in an image.
        ignore : Exclude data from the fit.
        show_filter : Show any filters applied to a data set.

        Notes
        -----
        The order of `ignore` and `notice` calls is important, and the
        results are a union, rather than intersection, of the
        combination.

        If `notice` is called on an un-filtered data set, then the
        ranges outside the noticed range are excluded: it can be
        thought of as if `ignore` had been used to remove all data
        points. If `notice` is called after a filter has been applied
        then the filter is applied to the existing data.

        For binned data sets, the bin is included if the noticed
        range falls anywhere within the bin.

        The units used depend on the `analysis` setting of the data
        set, if appropriate.

        To filter a 2D data set by a shape use `notice2d`.

        Examples
        --------

        Since the `notice` call is applied to an un-filtered
        data set, the filter choses only those points that lie
        within the range 12 <= X <= 18.

        >>> load_arrays(1, [10,15,20,30], [5,10,7,13])
        >>> notice(12, 28)
        >>> get_dep(filter=True)
        array([10,  7])

        As no limits are given, the whole data set is included:

        >>> notice()
        >>> get_dep(filter=True)
        array([ 5, 10,  7, 13])

        The `ignore` call excludes the first two points, but the
        `notice` call adds back in the second point:

        >>> ignore(None, 17)
        >>> notice(12, 16)
        >>> get_dep(filter=True)
        array([10,  7, 13])

        Only include data points in the range 8<=X<=12 and 18<=X=22:

        >>> ignore()
        >>> notice("8:12, 18:22")
        >>> get_dep(filter=True)
        array([5, 7])

        """
        if lo is not None and type(lo) in (str,numpy.string_):
            return self._notice_expr(lo, **kwargs)
        if self._data.items() == []:
            raise IdentifierErr("nodatasets")
        for d in self._data.values():
            d.notice(lo, hi, **kwargs)


    ### Ahelp ingest: 2015-05-06 DJB
    def ignore(self, lo=None, hi=None, **kwargs):
        """Exclude data from the fit.

        Select one or more ranges of data to exclude by filtering on
        the independent axis value. The filter is applied to all data
        sets.

        Parameters
        ----------
        lo : number or str, optional
           The lower bound of the filter (when a number) or a string
           expression listing ranges in the form `a:b`, with multiple
           ranges allowed, where the ranges are separated by a
           `,`. The term `:b` means exclude everything up to, and
           including `b`, and `a:` means exclude everything that is
           higher than, or equal to, `a`.
        hi : number, optional
           The upper bound of the filter when `lo` is not a string.

        See Also
        --------
        ignore_id : Exclude data from the fit for a data set.
        sherpa.astro.utils.ignore2d : Exclude a spatial region from an image.
        notice : Include data in the fit.
        show_filter : Show any filters applied to a data set.

        Notes
        -----
        The order of `ignore` and `notice` calls is important, and the
        results are a union, rather than intersection, of the
        combination.

        For binned data sets, the bin is excluded if the ignored
        range falls anywhere within the bin.

        The units used depend on the `analysis` setting of the data
        set, if appropriate.

        To filter a 2D data set by a shape use `ignore2d`.

        Examples
        --------

        Ignore all data points with an X value (the independent axis)
        between 12 and 18. For this one-dimensional data set, this
        means that the second bin is ignored:

        >>> load_arrays(1, [10,15,20,30], [5,10,7,13])
        >>> ignore(12, 18)
        >>> get_dep(filter=True)
        array([ 5,  7, 13])

        Filtering X values that are 25 or larger means that the last
        point is also ignored:

        >>> ignore(25, None)
        >>> get_dep(filter=True)
        array([ 5,  7])

        The `notice` call removes the previous filter, and then a
        multi-range filter is applied to exclude values between 8 and
        12 and 18 and 22:

        >>> notice()
        >>> ignore("8:12,18:22")
        >>> get_dep(filter=True)
        array([10, 13])

        """
        kwargs['ignore'] = True
        if lo is not None and type(lo) in (str,numpy.string_):
            return self._notice_expr(lo, **kwargs)
        self.notice(lo, hi, **kwargs)

    ### Ahelp ingest: 2015-05-06 DJB
    def notice_id(self, ids, lo=None, hi=None, **kwargs):
        """Include data from the fit for a data set.

        Select one or more ranges of data to include by filtering on
        the independent axis value. The filter is applied to the
        given data set, or data sets.

        Parameters
        ----------
        ids : int or str, or array of int or str
           The data set, or sets, to use.
        lo : number or str, optional
           The lower bound of the filter (when a number) or a string
           expression listing ranges in the form `a:b`, with multiple
           ranges allowed, where the ranges are separated by a
           `,`. The term `:b` means include everything up to, and
           including `b`, and `a:` means inlude everything that is
           higher than, or equal to, `a`.
        hi : number, optional
           The upper bound of the filter when `lo` is not a string.
        bkg_id : int or str, optional
           The filter will be applied to the associated background
           component of the data set if `bkg_id` is set. Only PHA
           data sets support this option; if not given, then the
           filter is applied to all background components as well
           as the source data.

        See Also
        --------
        ignore_id : Exclude data from the fit for a data set.
        sherpa.astro.utils.ignore2d : Exclude a spatial region from an image.
        notice : Include data in the fit.
        show_filter : Show any filters applied to a data set.

        Notes
        -----
        The order of `ignore` and `notice` calls is important.

        The units used depend on the `analysis` setting of the data
        set, if appropriate.

        To filter a 2D data set by a shape use `ignore2d`.

        Examples
        --------

        Include all data points with an X value (the independent axis)
        between 12 and 18 for data set 1:

        >>> notice_id(1, 12, 18)

        Include the range 0.5 to 7, for data sets 1,
        2, and 3:

        >>> notice_id([1,2,3], 0.5, 7)

        Apply the filter 0.5 to 2 and 2.2 to 7 to the data sets "core"
        and "jet":

        >>> notice_id(["core","jet"], "0.5:2, 2.2:7")

        """
        if self._valid_id(ids):
            ids = (ids,)
        else:
            try:
                ids = tuple(ids)
            except TypeError:
                _argument_type_error('ids',
                                     'an identifier or list of identifiers')

        if lo is not None and type(lo) in (str,numpy.string_):
            return self._notice_expr_id(ids, lo, **kwargs)
        for i in ids:
            self.get_data(i).notice(lo, hi, **kwargs)


    ### Ahelp ingest: 2015-05-06 DJB
    def ignore_id(self, ids, lo=None, hi=None, **kwargs):
        """Exclude data from the fit for a data set.

        Select one or more ranges of data to exclude by filtering on
        the independent axis value. The filter is applied to the given
        data set, or sets.

        Parameters
        ----------
        ids : int or str, or array of int or str
           The data set, or sets, to use.
        lo : number or str, optional
           The lower bound of the filter (when a number) or a string
           expression listing ranges in the form `a:b`, with multiple
           ranges allowed, where the ranges are separated by a
           `,`. The term `:b` means exclude everything up to, and
           including `b`, and `a:` means exclude everything that is
           higher than, or equal to, `a`.
        hi : number, optional
           The upper bound of the filter when `lo` is not a string.
        bkg_id : int or str, optional
           The filter will be applied to the associated background
           component of the data set if `bkg_id` is set. Only PHA
           data sets support this option; if not given, then the
           filter is applied to all background components as well
           as the source data.

        See Also
        --------
        ignore : Exclude data from the fit.
        sherpa.astro.utils.ignore2d : Exclude a spatial region from an image.
        notice_id : Include data from the fit for a data set.
        show_filter : Show any filters applied to a data set.

        Notes
        -----
        The order of `ignore` and `notice` calls is important.

        The units used depend on the `analysis` setting of the data
        set, if appropriate.

        To filter a 2D data set by a shape use `ignore2d`.

        Examples
        --------

        Ignore all data points with an X value (the independent axis)
        between 12 and 18 for data set 1:

        >>> ignore_id(1, 12, 18)

        Ignore the range up to 0.5 and 7 and above, for data sets 1,
        2, and 3:

        >>> ignore_id([1,2,3], None, 0.5)
        >>> ignore_id([1,2,3], 7, None)

        Apply the same filter as the previous example, but to
        data sets "core" and "jet":

        >>> ignore_id(["core","jet"], ":0.5,7:")

        """
        kwargs['ignore'] = True
        if lo is not None and type(lo) in (str,numpy.string_):
            return self._notice_expr_id(ids, lo, **kwargs)
        self.notice_id(ids, lo, hi, **kwargs)


    ###########################################################################
    # Models
    ###########################################################################

    ### Ahelp ingest: 2015-05-02 DJB
    def paramprompt(self, val=False):
        """Should the user be asked for the parameter values when creating a model?

        When `val` is `True`, calls to `set_model` will cause the user
        to be prompted for each parameter in the expression.  The
        prompt includes the parameter name and default value, in `[]`:
        the valid responses are

        - return  which accepts the default
        - value   which changes the parameter value
        - value, min  which changes the value and the minimum value
        - value, min, max  which changes the value, minimum, and
          maximum values

        The `value`, `min`, and `max` components are optional, so
        ",-5" will use the default parameter value and set its minimum
        to -5, while "2,,10" will change the parameter value to 2 and
        its maximum to 10, but leave the minimum at its default. If
        any value is invalid then the parameter is re-prompted.

        Parameters
        ----------
        val : bool, optional
           If `True`, the user will be prompted to enter each
           parameter value, including support for changing the minimum
           and maximum values, when a model component is created. The
           default is `False`.

        See Also
        --------
        set_model : Set the source model expression for a data set.
        set_par : Set the value, limits, or behavior of a model parameter.
        show_model : Display the model expression used to fit a data set.

        Notes
        -----
        Setting this to `True` only makes sense in an interactive
        environment.  It is designed to be similar to the parameter
        prompting provided by X-Spec [1]_.

        References
        ----------

        .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/

        Examples
        --------

        In the following, the default parameter settings are accepted
        for the `pl.gamma` parameter, the starting values for the
        `pl.ref` and `gline.pos` values are changed, the starting
        value and ranges of both the `pl.ampl` and `gline.ampl`
        parameters are set, and the `gline.fwhm` parameter is set to
        100, with its maximum changed to 10000.

        >>> paramprompt(True)
        >>> set_source(powlaw1d.pl + gauss1d.gline)
        pl.gamma parameter value [1]
        pl.ref parameter value [1] 4500
        pl.ampl parameter value [1] 1.0e-5,1.0e-8,0.01
        gline.fwhm parameter value [10] 100,,10000
        gline.pos parameter value [0] 4900
        gline.ampl parameter value [1] 1.0e-3,1.0e-7,1

        """
        self._paramprompt = sherpa.utils.bool_cast(val)


    def _add_model_types(self, module,
                         baselist=(sherpa.models.ArithmeticModel,)):
        if not isinstance(baselist, tuple):
            baselist = (baselist,)

        for name in module.__all__:
            cls = getattr(module, name)

            for base in baselist:
                if _is_subclass(cls, base):
                    break
            else:
                continue

            name = name.lower()
            self._model_types[name] = ModelWrapper(self, cls)
            self._model_globals.update(self._model_types)


    def add_model(self, modelclass, args=(), kwargs={}):
        """
        add_model

        SYNOPSIS
           Add a user-defined model class as a Sherpa model

        SYNTAX

        Arguments:
           modelclass     - User-defined model class

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           list_models
        """
        name = modelclass.__name__.lower()

        if not _is_subclass(modelclass, sherpa.models.ArithmeticModel):
            raise TypeError("model class '%s' is not a derived class" % name + 
                            " from sherpa.models.ArithmeticModel")

        self._model_types[name] = ModelWrapper(self, modelclass, args, kwargs)
        self._model_globals.update(self._model_types)
        _assign_obj_to_main(name, self._model_types[name])


    #
    # Model components
    #

    ### Ahelp ingest: 2015-05-04 DJB
    def get_model_autoassign_func(self):
        """Return the method used to create model component identifiers.

        Provides access to the function which is used by
        `create_model_component` and when creating model components
        directly to add an identifier in the current Python namespace.

        Returns
        -------
        func :
           The model function set by `set_model_autoassign_func`.

        See Also
        --------
        create_model_component : Create a model component.
        set_model : Set the source model expression for a data set.
        set_model_autoassign_func : Set the method used to create model component identifiers.

        """
        return self._model_autoassign_func

    ### Ahelp ingest: 2015-05-04 DJB
    ### DOC-TODO: what does func=None mean? If you try None then it
    ###           fails with AttributeError: 'Session' object attribute
    ###           '_model_autoassign_func' cannot be replaced with a non-callable attribute
    def set_model_autoassign_func(self, func=None):
        """Set the method used to create model component identifiers.

        When a model component is created, the default behavior is to
        add the component to the default Python namespace. This is
        controlled by a function which can be changed with this
        routine.

        Parameters
        ----------
        func : function reference
           The function to use: this should accept two arguments, a
           string (component name), and the model instance.

        See Also
        --------
        create_model_component : Create a model component.
        get_model_autoassign_func : Return the method used to create model component identifiers
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The default assignment function first renames a model
        component to include the model type and user-defined
        identifier. It then updates the '__main__' module's dictionary
        with the model identifier as the key and the model instance as
        the value. Similarly, it updates the '__builtin__' module's
        dictionary just like '__main__' for compatibility with
        IPython.

        """
        if (func is not None) and (not callable(func)):
            _argument_type_error('func', 'a function or other callable object')
        self._model_autoassign_func = func

    ### Ahelp ingest: 2015-04-27 DJB
    def list_models(self, show="all"):
        """List the available model types.

        Parameters
        ----------
        show : { 'all', '1d', '2d', 'xspec' }, optional
           What type of model should be returned. The default is
           'all'. An unrecognized value is treated as 'all'.

        Returns
        -------
        models : list of str

        See Also
        --------
        create_model_components : Create a model component.
        list_model_components : List the current model components.

        Examples
        --------

        >>> models = list_models()
        >>> models[0:5]
        ['absorptionedge',
         'absorptiongaussian',
         'absorptionlorentz',
         'absorptionvoigt',
         'accretiondisk']

        """
        keys = self._model_types.keys()[:]
        keys.sort()

        show = show.strip().lower()
        for item in ["arfmodel", "convolutionmodel", "multiresponsesummodel", "pileuprmfmodel",
                     "jdpileup", "rmfmodel", "tablemodel", "usermodel", "psfmodel"]:
            if item in keys:
                keys.remove(item)

        if show.startswith("xspec"):
            return filter(lambda x: x.startswith('xs'), keys)
        elif show.startswith("1d"):
            return filter(lambda x: (not x.startswith('xs')) and (not x.endswith("2d")) , keys)
        elif show.startswith("2d"):
            return filter(lambda x: x.endswith('2d'), keys)

        return keys

    ### Ahelp ingest: 2015-04-29 DJB
    def list_model_components(self):
        """List the names of all the model components.

        Models are created either directly - by using the form
        `mname.mid`, where `mname` is the name of the model, such as
        `gauss1d`, and `mid` is the name of the component - or with
        the `create_model_component` function, which accepts `mname`
        and `mid` as separate arguments. This function returns all the
        `mid` values that have been created.

        Returns
        -------
        ids : list of str
           The identifiers for all the model components that have been
           created. They do not need to be associated with a source
           expression (i.e. they do not need to have been included in
           a call to `set_model`).

        See Also
        --------
        create_model_component : Create a model component.
        delete_model_component : Delete a model component.
        list_models : List the available model types.
        list_model_ids : List of all the data sets with a source expression.
        set_model : Set the source model expression for a data set.

        Examples
        --------

        The `gal` and `pl` models are created - as versions
        of the `xsphabs` and `powlaw1d` model types - which means
        that the list of model components returned as `mids` will
        contain both strings.

        >>> set_model(xsphabs.gal * powlaw1d.pl)
        >>> mids = list_model_components()
        >>> 'gal' in mids
        True
        >>> 'pl' in mids
        True

        The model component does not need to be included as
        part of a source expression for it to be included in
        the output of this function:

        >>> create_model_component('gauss2d', 'gsrc')
        >>> 'gsrc' in list_model_components()
        True

        """
        keys = self._model_components.keys()[:]
        keys.sort()
        return keys

    def _add_model_component(self, cmpt):
        # If model component name is a model type name
        # or session function name, don't create it, raise
        # warning
        if (self._model_types.has_key(cmpt.name) is True):
            modeltype = cmpt.name
            del cmpt
            raise IdentifierErr('badidmodel', modeltype)

        if (cmpt.name.lower() in _builtin_symbols_):
            modeltype = cmpt.name
            del cmpt
            raise IdentifierErr('badidnative', modeltype)

        self._model_components[cmpt.name] = cmpt
        if self._model_autoassign_func is not None:
            self._model_autoassign_func(cmpt.name, cmpt)

    def _get_model_component(self, name, require=False):
        cmpt = self._model_components.get(name)
        require=sherpa.utils.bool_cast(require)
        if require and (cmpt is None):
           raise IdentifierErr('nomodelcmpt', name)
        return cmpt


    def _get_user_stat(self, statname):
        userstat = statname
        if not isinstance(statname, sherpa.stats.Stat):
            if (type(statname) is not str):
                _argument_type_error("stat name", "an instance or a string")
            userstat = self._get_model_component(statname, True)

        return userstat


    ### Ahelp ingest: 2015-05-04 DJB
    ### DOC-TODO: can send in a model variable, but this is just the
    ###           identity function, so not worth documenting
    def get_model_component(self, name):
        """Returns a model component given its name.

        Parameters
        ----------
        name : str
           The name of the model component.

        Returns
        -------
        component : a sherpa.models.model.Model instance
           The model component object.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If there is no model component with the given `name`.

        See Also
        --------
        create_model_component : Create a model component.
        get_model : Return the model expression for a data set.
        get_source : Return the source model expression for a data set.
        list_model_components : List the names of all the model components.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The model instances are named as modeltype.username, and it is
        the `username` component that is used here to access the
        instance.

        Examples
        --------

        When a model component is created, a variable is created that
        contains the model instance. The instance can also be returned
        with `get_model_component`, as shown here:

        >>> create_model_component('gauss1d', 'gline')
        >>> gmodel = get_model_component('gline')
        >>> gmodel.name
        'gauss1d.gline'
        >>> gmodel.pars
        (<Parameter 'fwhm' of model 'gline'>,
         <Parameter 'pos' of model 'gline'>,
         <Parameter 'ampl' of model 'gline'>)

        """   
        # If user mistakenly passes an actual model reference,
        # just return the reference
        if isinstance(name, sherpa.models.Model):
            return name

        _check_type(name, basestring, 'name', 'a string')
        return self._get_model_component(name, require=True)


    ### Ahelp ingest: 2015-04-29 DJB
    def create_model_component(self, typename=None, name=None):
        """Create a model component.

        Model components created by this function are set to their
        default values. Components can also be created directly using
        the syntax `typename.name`, such as in calls to `set_model`,
        when using the default model auto assignment setting (see
        `set_model_autoassign_func`).

        Parameters
        ----------
        typename : str
           The name of the model. This should match an entry from the
           return value of `list_models`, and defines the type of
           model.
        name : str
           The name used to refer to this instance, or component,
           of the model. A Python variable will be created with this
           name that can be used to inspect and change the model
           parameters, as well as use it in model expressions.

        See Also
        --------
        delete_model_component : Delete a model component.
        get_model_component : Returns a model component given its name.
        list_models : List the available model types.
        list_model_components : List the names of all the model components.
        set_model : Set the source model expression for a data set.
        set_model_autoassign_func : Set the method used to create model component identifiers.

        Notes
        -----
        This function can over-write an existing component. If the
        over-written component is part of a source expression - as set
        by `set_model` - then the model evaluation will still use the
        old model definition (and be able to change the fit
        parameters), but direct access to its parameters is not
        possible since the name now refers to the new component (this
        is true using direct access, such as `mname.parname`, or with
        `set_par`).

        Examples
        --------

        Create an instance of the `powlaw1d` model called `pl`,
        and then freeze its `gamma` parameter to 2.6.

        >>> create_model_component("powlaw1d", "pl")
        >>> pl.gamma = 2.6
        >>> freeze(pl.gamma)

        """

        # If user mistakenly passes an actual model reference,
        # just return (i.e., create_model_component(const1d.c1)
        # is redundant, so just return)
        if isinstance(typename, sherpa.models.Model) and name is None:
            return
        
        _check_type(typename, basestring, 'typename', 'a string')
        _check_type(name, basestring, 'name', 'a string')

        if typename is None:
            raise ArgumentErr('notype')
        if name is None:
            raise ArgumentErr('noname')
        typename = typename.lower()
        cls = self._model_types.get(typename)
        if cls is None:
            raise ArgumentErr('badtype', typename)

        self._model_components[name] = cls(name)
        #self._add_model_component(cls(name))

    ### Ahelp ingest: 2015-04-28 DJB
    def reset(self, model=None, id=None):
        """Reset the model parameters to their default settings.

        The `reset` function restores the parameter values to the
        default value set by `guess` or to the user-defined default.
        If the user set initial model values or soft limits -
        e.g. either with `set_par` or by using parameter prompting
        via `paramprompt` - then `reset` will restore these values
        and limits even after `guess` or `fit` has been called.

        Parameters
        ----------
        model : optional
           The model component or expression to reset. The default
           is to use all source expressions.
        id : int or string, optional
           The data set to use. The default is to use all
           data sets with a source expression.

        See Also
        --------
        fit : Fit one or more data sets.
        guess : Set model parameters to values matching the data. 
        paramprompt : Control how parameter values are set.
        set_par : Set a parameter value.

        Examples
        --------

        The following examples assume that the source model has been
        set using:

        >>> set_source(powlaw1d.pl * xsphabs.gal)

        Fit the model and then reset the values of both components
        (`pl` and `gal`):

        >>> fit()
        >>> reset()

        Reset just the parameters of the `pl` model component:

        >>> reset(pl)

        Reset all the components of the source expression for data
        set 2.

        >>> reset(get_source(2))

        """
        ids = [id]
        if id is None:
            ids = self._sources.keys()

        if model is None:
            for id in ids:
                model = self._get_source(id)
                model.reset()
        elif model is not None:
            model.reset()


    ### Ahelp ingest: 2015-04-29 DJB
    def delete_model_component(self, name):
        """Delete a model component.

        Parameters
        ----------
        name : str
           The name used to refer to this instance, or component, of
           the model. The corresponding Python variable will be
           deleted by this function.

        See Also
        --------
        create_model_component : Create a model component.
        delete_model : Delete the model expression for a data set.
        list_models : List the available model types.
        list_model_components : List the names of all the model components.
        set_model : Set the source model expression for a data set.
        set_model_autoassign_func :

        Notes
        -----
        It is an error to try to delete a component that is part of a
        model expression - i.e. included as part of an expression in a
        `set_model` or `set_source` call. In such a situation, use the
        `delete_model` function to remove the source expression before
        calling `delete_model_component`.

        Examples
        --------

        If a model instance called `pl` has been created - e.g. by
        `create_model_component('powlaw1d', 'pl') - then the
        following will remove it:

        >>> delete_model_component('pl')

        """
        _check_type(name, basestring, 'name', 'a string')
        mod = self._model_components.pop(name, None)
        if mod is not None:
            for key in self.list_model_ids():
                if mod in self._models[key] or mod in self._sources[key]:
                    warning("the model component '%s' is found in model %s" %
                            (mod.name, str(key) + " and cannot be deleted" ))
                    # restore the model component in use and return
                    self._model_components[name]=mod
                    return

            del sys.modules["__main__"].__dict__[name]
            del sys.modules["__builtin__"].__dict__[name]
        else:
            raise IdentifierErr('nomodelcmpt', name)
        
    # Back-compatibility
    #create_model = create_model_component

    #
    # Source models
    #

    def _eval_model_expression(self, expr, typestr='model'):
        try:
            return eval(expr, self._model_globals, self._model_components)
        except:
            raise ArgumentErr('badexpr', typestr, sys.exc_info()[1])

    ### Ahelp ingest: 2015-04-29 DJB
    def list_model_ids(self):
        """List of all the data sets with a source expression.

        Returns
        -------
        ids : list of int or str
           The identifiers for all the data sets which have a source
           expression set by `set_model` or `set_source`.

        See Also
        --------
        list_data_ids : List the identifiers for the loaded data sets.
        list_model_components : List the names of all the model components.
        set_model : Set the source model expression for a data set.

        """
        keys = self._models.keys()[:]
        keys.extend(self._sources.keys()[:])
        keys = list(set(keys))
        keys.sort()
        return keys

    # Return full model for fitting, plotting, etc.  Expects a corresponding
    # data set to be available.
    def _get_model_status(self, id=None):
        id = self._fix_id(id)
        src = self._sources.get(id)
        mdl = self._models.get(id)

        if src is None and mdl is None:
            raise IdentifierErr('getitem', 'model', id, 'has not been set')

        model = mdl
        is_source = False

        if mdl is None and src is not None:
            is_source = True
            model = src

        return (model, is_source)


    def _add_convolution_models(self, id, data, model, is_source):

        [tbl.fold(data) for tbl in self._tbl_models if tbl in model]

        if is_source:
            model = self._add_psf(id, data, model)
        else:
            [psf.fold(data) for psf in self._psf_models if psf in model]

        return model


    def _get_model(self, id=None):
        data = self.get_data(id)
        model, is_source = self._get_model_status(id)
        return self._add_convolution_models(id, data, model, is_source)


    def _get_source(self, id=None):
        id = self._fix_id(id)
        mdl = self._models.get(id, None)
        if mdl is not None:
            raise IdentifierErr("Convolved model\n'%s'\n is set for dataset %s. You should use get_model instead." % 
                    (mdl.name, str(id)))
        return self._get_item(id, self._sources, 'source',
                              'has not been set, consider using set_source()' +
                              ' or set_model()')

    ### Ahelp ingest: 2015-05-08 DJB
    def get_source(self, id=None):
        """Return the source model expression for a data set.

        This returns the model expression created by `set_model` or
        `set_source`. It does not include any instrument response.

        Parameters
        ----------
        id : int or str, optional
           The data set containing the source expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        model : a sherpa.models.Model object
           This can contain multiple model components. Changing
           attributes of this model changes the model used by the data
           set.

        See Also
        --------
        delete_model : Delete the model expression from a data set.
        get_model : Return the model expression for a data set.
        get_model_pars : Return the names of the parameters of a model.
        get_model_type : Describe a model expression.
        list_model_ids : List of all the data sets with a source expression.
        sherpa.astro.utils.set_bkg_model : Set the background model expression for a data set.
        set_model : Set the source model expression for a data set.
        set_full_model : Define the convolved model expression for a data set.
        show_model : Display the source model expression for a data set.

        Examples
        --------

        Return the source expression for the default data set:

        >>> src = get_source()
        >>> len(src.pars)
        5

        Set the source expression for data set 'obs2' to be equal to
        the model of data set 'obs1' multiplied by a scalar value:

        >>> set_source('obs2', const1d.norm * get_source('obs1'))

        """
        return self._get_source(id)

    ### Ahelp ingest: 2015-05-08 DJB
    def get_model(self, id=None):
        """Return the model expression for a data set.

        This returns the model expression for a data set, including
        any instrument response (e.g. PSF or ARF and RMF) whether
        created automatically or explicitly, with `set_full_model`.

        Parameters
        ----------
        id : int or str, optional
           The data set containing the source expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        model :
           This can contain multiple model components and any
           instrument response. Changing attributes of this model
           changes the model used by the data set.

        See Also
        --------
        delete_model : Delete the model expression from a data set.
        get_model_pars : Return the names of the parameters of a model.
        get_model_type : Describe a model expression.
        get_source : Return the source model expression for a data set.
        list_model_ids : List of all the data sets with a source expression.
        sherpa.astro.utils.set_bkg_model : Set the background model expression for a data set.
        set_model : Set the source model expression for a data set.
        set_full_model : Define the convolved model expression for a data set.
        show_model : Display the source model expression for a data set.

        Examples
        --------

        Return the model fitted to the default data set:

        >>> mdl = get_model()
        >>> len(mdl.pars)
        5

        """
        return self._get_model(id)


    def _runparamprompt(self, pars):

        # paramprompt
        if self._paramprompt:
            for par in pars:
                while True:
                    input = raw_input("%s parameter value [%g] " %
                                      (par.fullname, par.val))
                    if input != "":
                        val = min = max = None
                        count = input.count(',')

                        if count == 0:
                            try:
                                val = float(input)
                            except Exception, e:
                                info("Please provide a float value; " + str(e))
                                continue

                        elif count == 1:
                            try:
                                str_val,str_min = input.split(',')
                                if str_val != "":
                                    val = float(str_val)
                                if str_min != "":
                                    min = float(str_min)
                            except Exception, e:
                                info("Please provide a float value; " + str(e))
                                continue

                        elif count == 2:
                            try:
                                str_val,str_min,str_max = input.split(',')
                                
                                if str_val != "":
                                    val = float(str_val)
                                if str_min != "":
                                    min = float(str_min)
                                if str_max != "":
                                    max = float(str_max)
                            except Exception, e:
                                info("Please provide a float value; " + str(e))
                                continue
                        else:
                            info("Error: Please provide a comma separated list"+
                                 " of floats; e.g. val,min,max")
                            continue

                        try:
                            self.set_par(par,val,min,max)
                            break
                        except Exception, e:
                            info(str(e))
                            continue
                    else:
                        break

    # DOC-NOTE: also in sherpa.astro.utils
    ##@loggable(with_id=True, with_keyword='model')
    def set_full_model(self, id, model=None):
        """
        set_full_model

        SYNOPSIS
           Set a convolved Sherpa model by model id

        SYNTAX

        Arguments:
           id         - model id
                        default = default model id

           model      - Sherpa model

        Returns:
           None

        DESCRIPTION
           Add a Sherpa model to the list of current active Sherpa model
           by model id.

        SEE ALSO
           list_model_ids, get_model, delete_model, get_model_type,
           get_model_pars        
        """
        if model is None:
            id, model = model, id
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        self._set_item(id, model, self._models, sherpa.models.Model, 'model',
                       'a model object or model expression string')
        self._runparamprompt(model.pars)

    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: the .cache value appears to default to 5
    ##@loggable(with_id=True, with_keyword="model")
    def set_model(self, id, model=None):
        """Set the source model expression for a data set.

        The function is available as both `set_model` and
        `set_source`. The model fit to the data can be further
        modified by instrument responses which can be set explicitly -
        e.g. by `set_psf` - or be defined automatically by the type of
        data being used (e.g. the ARF and RMF of a PHA data set). The
        `set_full_model` command can be used to explicitly include the
        instrument response if necessary.

        Parameters
        ----------
        id : int or str, optional
           The data set containing the source expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        model : str or sherpa.models.Model object
           This defines the model used to fit the data. It can be a
           Python expression or a string version of it.

        See Also
        --------
        delete_model : Delete the model expression from a data set.
        fit : Fit one or more data sets.
        freeze : Fix model parameters so they are not changed by a fit.
        get_source : Return the source model expression for a data set.
        integrate1d : Integrate 1D source expressions.
        sherpa.astro.utils.set_bkg_model : Set the background model expression for a data set.
        set_full_model : Define the convolved model expression for a data set.
        show_model : Display the source model expression for a data set.
        thaw : Allow model parameters to be varied during a fit.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        PHA data sets will automatically apply the instrumental
        response (ARF and RMF) to the source expression. For some
        cases this is not useful - for example, when different
        responses should be applied to different model components - in
        which case `set_full_model` should be used instead.

        Model caching is available via the model `cache` attribute. A
        non-zero value for this attribute means that the results of
        evaluating the model will be cached if all the parameters are
        frozen, which may lead to a reduction in the time taken to
        evaluate a fit. A zero value turns off the cacheing.  The
        default setting for X-Spec and 1D analytic models is that
        `cache` is `5`, but `0` for the 2D analytic models.

        The `integrate1d` model can be used to apply a numerical
        integration to an arbitrary model expression.

        Examples
        --------

        Create an instance of the `powlaw1d` model type, called `pl`,
        and use it as the model for the default data set.

        >>> set_model(polynom1d.pl)

        Create a model for the default dataset which is the `xsphabs`
        model multiplied by the sum of an `xsapec` and `powlaw1d`
        models (the model components are identified by the labels
        `gal`, `clus`, and `pl`).

        >>> set_model(xsphabs.gal * (xsapec.clus + powlaw1d.pl))

        Repeat the previous example, using a string to define the
        model expression:

        >>> set_model('xsphabs.gal * (xsapec.clus + powlaw1d.pl)')

        Use the same model component (`src`, a `gauss2d` model)
        for the two data sets ('src1' and 'src2').

        >>> set_model('src1',  gauss2d.src + const2d.bgnd1)
        >>> set_model('src2', src + const2d.bgnd2)

        Share an expression - in this case three gaussian lines -
        between three data sets. The normalization of this line
        complex is allowed to vary in data sets 2 and 3 (the `norm2`
        and `norm3` components of the `const1d` model), and each data
        set has a separate `polynom1d` component (`bgnd1`, `bgnd2`,
        and `bgnd3`). The `c1` parameters of the `polynom1d` model
        components are thawed and then linked together (to reduce the
        number of free parameters):

        >>> lines = gauss1d.l1 + gauss1d.l2 + gauss1d.l3
        >>> set_model(1, lines + polynom1d.bgnd1)
        >>> set_model(2, lines * const1d.norm2 + polynom1d.bgnd2)
        >>> set_model(3, lines * const1d.norm3 + polynom1d.bgnd3)
        >>> thaw(bgnd1.c1, bgnd2.c1, bgnd3.c1)
        >>> link(bgnd2.c2, bgnd1.c1)
        >>> link(bgnd3.c3, bgnd1.c1)

        For this expression, the `gal` component is frozen, so it is
        not varied in the fit. The `cache` attribute is set to a
        non-zero value to ensure that it is cached during a fit (this
        is actually the default value for this model so it not
        normally needed).

        >>> set_model(xsphabs.gal * (xsapec.clus + powlaw1d.pl))
        >>> gal.nh = 0.0971
        >>> freeze(gal)
        >>> gal.cache = 1

        """
        if model is None:
            id, model = model, id
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        self._set_item(id, model, self._sources, sherpa.models.Model,
                       'source model',
                       'a model object or model expression string')

        self._runparamprompt(model.pars)

        # Delete any previous model set with set_full_model()
        id = self._fix_id(id)
        mdl = self._models.pop(id, None)
        if mdl is not None:
            warning("Clearing convolved model\n'%s'\nfor dataset %s" % 
                    (mdl.name, str(id)))


    set_source = set_model


    ### Ahelp ingest: 2015-04-29 DJB
    def delete_model(self, id=None):
        """Delete the model expression for a data set.

        This removes the model expression, created by `set_model`,
        for a data set. It does not delete the components of the
        expression.

        Parameters
        ----------
        id : int or str, optional
           The data set containing the source expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.

        See Also
        --------
        clean : Clear all stored session data.
        delete_data : Delete a data set by identifier.
        get_default_id : Return the default data set identifier.
        set_model : Set the source model expression for a data set.
        show_model : Display the source model expression for a data set.

        Examples
        --------

        Remove the model expression for the default data set:

        >>> delete_model()

        Remove the model expression for the data set with the
        identifier called 'src':

        >>> delete_model('src')

        """
        id = self._fix_id(id)
        self._models.pop(id, None)
        self._sources.pop(id, None)


    def _check_model(self, model):
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)
        _check_type(model, sherpa.models.Model, 'model',
                    'a model object or model expression string')
        return model

    ### Ahelp ingest: 2015-05-08 DJB
    def get_model_type(self, model):
        """Describe a model expression.

        Parameters
        ----------
        model : str or a sherpa.models.model.Model object

        Returns
        -------
        type : str
           The name of the model expression.

        See Also
        --------
        create_model_component : Create a model component.
        get_model : Return the model expression for a data set.
        get_model_pars : Return the names of the parameters of a model.
        get_source : Return the source model expression for a data set.

        Examples
        --------

        >>> create_model_component("powlaw1d", "pl")
        >>> get_model_type("pl")
        'powlaw1d'

        For expressions containing more than one component, the
        result is likely to be 'binaryopmodel'

        >>> get_model_type(const1d.norm * (polynom1d.poly + gauss1d.gline))
        'binaryopmodel'

        For sources with some form of an instrument model - such as a
        PSF convolution for an image or a PHA file with response
        information from the ARF and RMF - the response can depend on
        whether the expression contains this extra information or not:

        >>> get_model_type(get_source('spec'))
        'binaryopmodel'
        >>> get_model_type(get_model('spec'))
        'rspmodelpha'

        """
        model = self._check_model(model)
        return type(model).__name__.lower()

    ### Ahelp ingest: 2015-05-08 DJB
    def get_model_pars(self, model):
        """Return the names of the parameters of a model.

        Parameters
        ----------
        model : str or a sherpa.models.model.Model object

        Returns
        -------
        names : list of str
           The names of the parameters in the model expression.  These
           names do not include the name of the parent component.

        See Also
        --------
        create_model_component : Create a model component.
        get_model : Return the model expression for a data set.
        get_model_type : Describe a model expression.
        get_source : Return the source model expression for a data set.

        Examples
        --------

        >>> set_source(gauss2d.src + const2d.bgnd)
        >>> get_model_pars(get_source())
        ['fwhm', 'xpos', 'ypos', 'ellip', 'theta', 'ampl', 'c0']

        """
        model = self._check_model(model)
        return [p.name for p in model.pars]

    ### Ahelp ingest: 2015-05-01 DJB
    def get_num_par(self, id=None):
        """Return the number of parameters in a model expression.

        Parameters
        ----------
        id : int or str, optional
           The data set containing the model expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        npar : int

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If no model expression has been set for the data set
           (with `set_model` or `set_source`).

        See Also
        --------
        get_num_par_frozen : Return the number of frozen parameters.
        get_num_par_thawed : Return the number of thawed parameters.
        set_model : Set the source model expression for a data set.

        """
        return len(self._get_source(id).pars)
    
    ### Ahelp ingest: 2015-05-01 DJB
    def get_num_par_thawed(self, id=None):
        """Return the number of thawed parameters in a model expression.

        Parameters
        ----------
        id : int or str, optional
           The data set containing the model expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        npar : int

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If no model expression has been set for the data set
           (with `set_model` or `set_source`).

        See Also
        --------
        get_num_par : Return the number of parameters.
        get_num_par_frozen : Return the number of frozen parameters.
        set_model : Set the source model expression for a data set.

        """
        return len(self._get_source(id).thawedpars)

    ### Ahelp ingest: 2015-05-01 DJB
    def get_num_par_frozen(self, id=None):
        """Return the number of frozen parameters in a model expression.

        Parameters
        ----------
        id : int or str, optional
           The data set containing the model expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        npar : int

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If no model expression has been set for the data set
           (with `set_model` or `set_source`).

        See Also
        --------
        get_num_par : Return the number of parameters.
        get_num_par_thawed : Return the number of thawed parameters.
        set_model : Set the source model expression for a data set.

        """
        model = self._get_source(id)
        return len(model.pars)-len(model.thawedpars)

    #
    # "Special" models (user and table models)
    #

    def _read_user_model(self, filename, ncols=2, colkeys=None,
                         dstype=sherpa.data.Data1D, sep=' ', comment='#'):
        x = None
        y = None
        try:
            data = self.unpack_data(filename, ncols, colkeys,
                                    dstype, sep, comment)
            y = data.get_y()
            x = data.get_x()
        # we have to check for the case of a *single* column in ascii file
        # extract the single array from the read and bypass the dataset
        except TypeError:
            y = sherpa.io.get_ascii_data(filename, ncols=1, colkeys=colkeys,
                                           sep=sep, dstype=dstype,
                                           comment=comment)[1].pop()
        except:
            raise
        return (x,y)

    def load_template_model(self, modelname, templatefile, dstype=sherpa.data.Data1D,
                            sep=' ', comment='#', method=sherpa.utils.linear_interp, template_interpolator_name='default'):

        if sherpa.utils.is_binary_file(templatefile):
            raise sherpa.utils.err.IOErr('notascii', templatefile)

        names, cols = sherpa.io.read_file_data(templatefile,
                                               sep=sep, comment=comment, require_floats=False)

        if len(names) > len(cols):
            raise sherpa.utils.err.IOErr('toomanycols')

        names = [name.strip().lower() for name in names]

        filenames = None
        modelflags = None
        parnames = names[:]
        parvals = []
        for name, col in izip(names, cols):
            # Find the column with the filenames, remove it from the set of
            # parameter names
            if name.startswith('file'):
                filenames = col
                parnames.remove(name)
                continue
            # Find the column with the modelflags, remove it from the set of
            # parameter names
            if name.startswith('modelflag'):
                modelflags = col
                parnames.remove(name)
                continue
            parvals.append(numpy.array(col, dtype=SherpaFloat))

        parvals = numpy.asarray(parvals).T

        if len(parvals) == 0:
            raise sherpa.utils.err.IOErr('noparamcols', templatefile)

        if filenames is None:
            raise sherpa.utils.err.IOErr('reqcol', 'filename', templatefile)

        if modelflags is None:
            raise sherpa.utils.err.IOErr('reqcol', 'modelflag', templatefile)

        templates = []
        for filename in filenames:
            tnames, tcols = sherpa.io.read_file_data(filename,
                                                     sep=sep, comment=comment)
            if len(tcols) == 1:
                raise sherpa.utils.err.IOErr('onecolneedtwo', filename)
            elif len(tcols) == 2:
                tm = sherpa.models.TableModel(filename)
                tm.method = method  # interpolation method
                tm.load(*tcols)
                tm.ampl.freeze()
                templates.append(tm)
            else:
                raise sherpa.utils.err.IOErr('toomanycols', str(tnames), '2')

        assert( len(templates) == parvals.shape[0] )

        templatemodel = sherpa.models.create_template_model(modelname, parnames, parvals, 
                                                            templates, template_interpolator_name)
        self._tbl_models.append(templatemodel)
        self._add_model_component(templatemodel)


    ##@loggable()
    def load_template_interpolator(self, name, interpolator_class, **kwargs):
	sherpa.models.template.interpolators[name] = (interpolator_class, kwargs)


    def load_table_model(self, modelname, filename, ncols=2, colkeys=None,
                         dstype=sherpa.data.Data1D, sep=' ', comment='#',
                         method=sherpa.utils.linear_interp):
        """
        load_table_model

        SYNOPSIS
           Load a table model from file into a Sherpa session

        SYNTAX

        Arguments:
           modelname  - model label

           filename   - file from which table model data are read

           ncols      - number of columns to read from
                        default = 2

           colkeys    - column keys
                        default = None

           dstype     - Sherpa data class to contain table model data
                        default = sherpa.data.Data1D
        
           sep        - separator character
                        default = ' '

           comment    - comment character
                        default = '#'

           method     - interpolation method
                        default = linear {neville, linear}

        Returns:
           None

        DESCRIPTION
           Load data from a file, and put it in a new model.  This
           model can be used in fitting, just as models that containing
           functions can be used.
        
        SEE ALSO
           set_model, load_user_model, add_user_pars
        """
        tablemodel = sherpa.models.TableModel(modelname)
        # interpolation method
        tablemodel.method = method 
        tablemodel.filename = filename
        x, y = self._read_user_model(filename, ncols, colkeys,
                                     dstype, sep, comment)
        tablemodel.load(x,y)
        self._tbl_models.append(tablemodel)
        self._add_model_component(tablemodel)

    ##@loggable()
    def load_user_model(self, func, modelname, filename=None, ncols=2,
                        colkeys=None, dstype=sherpa.data.Data1D,
                        sep=' ', comment='#'):
        """
        load_user_model

        SYNOPSIS
           Load a user model from file into a Sherpa session

        SYNTAX

        Arguments:
           func       - reference to a user model function
           
           modelname  - model label

           filename   - file from which optional data are read
                        default = None

           ncols      - number of columns to read from
                        default = 2

           colkeys    - column keys
                        default = None

           dstype     - Sherpa data class to contain table model data
                        default = sherpa.data.Data1D
        
           sep        - separator character
                        default = ' '
        
           comment    - comment character
                        default = '#'

        Returns:
           None

        DESCRIPTION
           Take a function written by the user, and assign to a new
           user model class.  Instances of the new class can be created,
           and used as models during fits--just as ordinary Sherpa
           models can.  Optionally, data from a file can be attached to
           the model, and used in an arbitrary way by the user model
           function; but data from file is not required, the user model
           can be just a function.  After a user model is created,
           parameters need to be added with the add_user_pars function.
        
        SEE ALSO
           set_model, load_table_model, add_user_pars
        """
        usermodel = sherpa.models.UserModel(modelname)
        usermodel.calc = func
        usermodel._file = filename
        if (filename is not None):
            x, usermodel._y = self._read_user_model(filename, ncols, colkeys,
                                                    dstype, sep, comment)
        self._add_model_component(usermodel)

    ##@loggable()
    def add_user_pars(self, modelname, parnames,
                      parvals = None, parmins = None, parmaxs = None,
                      parunits = None, parfrozen = None):
        """
        add_user_pars

        SYNOPSIS
           Add parameters to a user model

        SYNTAX

        Arguments:
           modelname  - model label

           parnames   - list of parameter names

           parvals    - list of parameter values
                        default = None

           parmins    - list of parameter minimum values
                        default = None

           parmaxs    - list of parameter maxinum values
                        default = None

           parunits   - list of parameter units
                        default = None

           parfrozen  - list of frozen parameter names
                        default = None

        Returns:
           None

        DESCRIPTION
           Add a new set of parameters to a user model that can be 
           used in Sherpa model definitions.
        
        SEE ALSO
           set_model, load_user_model
        """
        pars = []
        vals=None
        if parvals is not None:
            vals = list(parvals)[:]
        mins=None
        if parmins is not None:
            mins = list(parmins)[:]
        maxs=None
        if parmaxs is not None:
            maxs = list(parmaxs)[:]
        units=None
        if parunits is not None:
            units = list(parunits)[:]
        frozen=None
        if parfrozen is not None:
            frozen = list(parfrozen)[:]
        
        for name in parnames:
            par = sherpa.models.Parameter(modelname, name, 0.0)
            if parvals is not None:
                par.val = vals.pop(0)
            if parmins is not None:
                par.min = mins.pop(0)
            if parmaxs is not None:
                par.max = maxs.pop(0)
            if parunits is not None:
                par.units = units.pop(0)
            if parfrozen is not None:
                par.frozen = frozen.pop(0)
            pars.append(par)

        if (type(modelname) is not str):
            _argument_type_error("model name", "a string")
        usermodel = self._get_model_component(modelname)
        # If not a user model, exit
        if (usermodel is None or
            type(usermodel) is not sherpa.models.UserModel):
            _argument_type_error(modelname, "a user model")

        # Create a new user model with the desired parameters,
        # and copy over calc, file and y from the old usermodel
        # (these are the only attributes we care about preserving)
        newusermodel = sherpa.models.UserModel(modelname, pars)
        newusermodel.calc = usermodel.calc
        newusermodel._file = usermodel._file
        newusermodel._y = usermodel._y

        # Remove old usermodel from collection of known model
        # instances, and add new user model to that collection
        self._model_components.pop(modelname, None)
        self._add_model_component(newusermodel)


    def load_user_stat(self, statname, calc_stat_func, calc_err_func=None,
                       priors={}):
        """
        load_user_stat

        SYNOPSIS
           Load a user defined statistic into a Sherpa session

        SYNTAX

        Arguments:
           statname        - reference to a user statistic

           calc_stat_func  - function to calculate and return the statistic
                             value and the statistic contribution per bin as
                             a tuple.

           calc_err_func   - function to calculate the statistical error.
                             default = None

           priors          - dictionary of model parameters and 
                             hyper-parameters for priors.
                             default = {}

        Returns:
           None

        DESCRIPTION
           Load a user-defined statistic from user_calc_stat and 
           user_calc_err functions with identifier statname.  Optionally,
           include a dictionary of hyper-parameters and references to 
           source model parameters for priors statistics.  The symbols
           supplied in the dictionary are available in the first argument
           of the user_calc_stat function signature for use at every model
           evaluation.

        SEE ALSO
           set_set, calc_stat, calc_chisqr
        """
        userstat = sherpa.stats.UserStat(calc_stat_func, 
                                         calc_err_func, statname)
        if priors:
            pars = [(key, priors.pop(key)) for key in priors.keys()
                    if isinstance(priors[key], sherpa.models.Parameter)]
            pars = dict(pars)
            userstat = sherpa.logposterior.Prior(calc_stat_func, priors, pars)

        _assign_obj_to_main(statname, userstat)


    # Back-compatibility
    #set_source = set_model

    #
    # Conv
    #

    ##@loggable()
    def load_conv(self, modelname, filename_or_model, *args, **kwargs):
        """
        load_conv

        SYNOPSIS
           load a file-based or model-based 1D kernel into a 1D convolution model

        SYNTAX

        Arguments:
           modelname - name of convolution kernel

           filename_or_model - filename with path for file-based kernel
                               a Sherpa model for a model-based kernel

           args      - additional arguments when reading a file kernel

           kwargs    - additional keyword arguments when reading a file
                       kernel

        Returns:
           None

        DESCRIPTION
           Create a convolution model object with identifier 'modelname' and 
           initializes the convolution kernel to be either a Sherpa dataset
           loaded from file or a Sherpa model.

           NOTE: load_conv() is intended for 1D convolution only.  It uses the
                 midpoint as the origin.

        SEE ALSO
           set_psf, get_psf, delete_psf
        """        
        kernel = filename_or_model
        if isinstance(filename_or_model, basestring):
            try:
                kernel = self._eval_model_expression(filename_or_model)
            except:
                try:
                    kernel = self.unpack_data(filename_or_model,
                                              *args, **kwargs)
                except:
                    raise

        conv = sherpa.instrument.ConvolutionKernel(kernel,modelname)
        self._add_model_component(conv)


    #
    # PSF1
    #

    ### Ahelp ingest: 2015-05-11 DJB
    ### DOC-TODO: what args/kwargs are supported?
    def load_psf(self, modelname, filename_or_model, *args, **kwargs):
        """Create a PSF model.

        Create a PSF model representing either an array of data, read
        from a file, or a model component (such as a gaussian). The
        `set_psf` function is used to associate this model with a data
        set.

        Parameters
        ----------
        modelname : str
           The identifier for this PSF model.
        filename_or_model : str or model instance
           The form of the PSF. This can be a file name, which will
           be read in using the chosen Sherpa I/O library, or a
           model component. When reading in a model, additional
           arguments can be used to specify the data format (see
           `load_table` and `load_image` for more information).

        See Also
        --------
        delete_psf : Delete the PSF model for a data set.
        load_image :
        load_table :
        set_psf : Add a PSF model to a data set.

        Examples
        --------

        Create a PSF model using a 2D gaussian:

        >>> load_psf('psf1', gauss2d.gpsf)
        >>> set_psf('psf1')
        >>> gpsf.fwhm = 4.2
        >>> gpsf.ellip = 0.2
        >>> gpsf.theta = 30 * np.pi / 180
        >>> image_psf()

        Create a PSF model from the data in the file
        `line_profile.fits` and apply it to the data set called
        `bgnd`:

        >>> load_psf('pmodel', 'line_profile.fits')
        >>> set_psf('bgnd', 'pmodel')

        """        
        kernel = filename_or_model
        if isinstance(filename_or_model, basestring):
            try:
                kernel = self._eval_model_expression(filename_or_model)
            except:
                try:
                    kernel = self.unpack_data(filename_or_model,
                                              *args, **kwargs)
                except:
                    raise

        psf = sherpa.instrument.PSFModel(modelname, kernel)
        self._add_model_component(psf)
        self._psf_models.append(psf)

    ### Ahelp ingest: 2015-05-11 DJB
    ### DOC-TODO: am I correct about the multiple use warning?
    ##@loggable(with_id=True, with_keyword='psf')
    def set_psf(self, id, psf=None):
        """Add a PSF model to a data set.

        After this call, the model that is fit to the data (as set by
        `set_model`) will be convolved by the given PSF model. The
        term "psf" is used in functions to refer to the data sent to
        this function whereas the term "kernel" refers to the data
        that is used in the actual convolution (this can be
        re-normalized and a sub-set of the PSF data).

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the
           default identifier is used, as returned by `get_default_id`.
        psf : str or sherpa.instrument.PSFModel instance
           The PSF model created by `load_psf`.

        See Also
        --------
        delete_psf : Delete the PSF model for a data set.
        get_psf : Return the PSF model defined for a data set.
        image_psf : Display the 2D PSF model for a data set in the image viewer.
        load_psf : Create a PSF model.
        plot_psf : Plot the 1D PSF model applied to a data set.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `psf` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `psf` parameters,
        respectively.

        A PSF component should only be applied to a single data set.
        This is not enforced by the system, and incorrect results
        can occur if this condition is not true.

        The point spread function (PSF) is defined by the full
        (unfiltered) PSF image loaded into Sherpa or the PSF model
        expression evaluated over the full range of the dataset; both
        types of PSFs are established with the `load_psf` command.
        The kernel is the subsection of the PSF image or model which
        is used to convolve the data. This subsection is created from
        the PSF when the size and center of the kernel are defined by
        the command `set_psf`. While the kernel and PSF might be
        congruent, defining a smaller kernel helps speed the
        convolution process by restricting the number of points within
        the PSF that must be evaluated.

        In a 1-D PSF model, a radial profile or 1-D model array is
        used to convolve (fold) the given source model using the Fast
        Fourier Transform (FFT) technique. In a 2-D PSF model, an
        image or 2-D model array is used.

        The parameters of a PSF model include:

        kernel
           The data used for the convolution (file name or model
           instance).

        size
           The number of pixels used in the convolution (this can be a
           subset of the full PSF). This is a scalar (1D) or a
           sequence (2D, width then height) value.

        center
           The center of the kernel. This is a scalar (1D) or a
           sequence (2D, width then height) value. The kernel centroid
           must always be at the center of the extracted sub-image,
           otherwise, systematic shifts will occur in the best-fit
           positions.

        radial
           Set to `1` to use a symmetric array. The default is `0` to
           reduce edge effects.

        norm
           Should the kernel be normalized so that it sums to 1?  This
           summation is done over the full data set (not the subset
           defined by the `size` parameter). The default is `1` (yes).

        Examples
        --------

        Use the data in the file `line_profile.fits` as the PSF for
        the default data set:

        >>> load_psf('psf1', 'line_profile.fits')
        >>> set_psf(psf1)

        Use the same PSF for different data sets:

        >>> load_psf('p1', 'psf.img')
        >>> load_psf('p2', 'psf.img')
        >>> set_psf(1, 'p1')
        >>> set_psf(2, 'p2')

        Restrict the convolution to a sub-set of the PSF data and
        compare the two:

        >>> set_psf(psf1)
        >>> psf1.size = (41,41)
        >>> image_psf()
        >>> image_kernel(newframe=True, tile=True)

        """
        if psf is None:
            id, psf, = psf, id
        id = self._fix_id(id)

        if isinstance(psf, basestring):
            psf = self._eval_model_expression(psf)

        self._set_item(id, psf, self._psf, sherpa.instrument.PSFModel, 'psf',
                       'a PSF model object or PSF model expression string')

        # fold the PSF with data and model if available, if not pass
        try:
            data = self.get_data(id)
            psf = self._psf.get(id, None)
            if psf is not None:
                psf.fold(data)

        except IdentifierErr:
            pass

        # If the PSF is Sherpa model and it implements the center
        # attribute, then populate the model position parameters
        # using the PSF center if the user has not already done so.
        # Note: PSFKernel only.
        if (psf.kernel is not None and callable(psf.kernel) and 
            psf.model is not None and isinstance(psf.model, sherpa.instrument.PSFKernel)):

            psf_center = psf.center
            if numpy.isscalar(psf_center):
                psf_center = [psf_center]
            try:
                center = psf.kernel.get_center()
                if (numpy.asarray(center) == 0.0).all():
                    psf.kernel.set_center(*psf_center, values=True)
            except NotImplementedError:
                pass


    ### Ahelp ingest: 2015-05-11 DJB
    def get_psf(self, id=None):
        """Return the PSF model defined for a data set.

        Return the parameter settings for the PSF model assigned to
        the data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        psf : sherpa.instrument.PSFModel instance

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If no PSF model has been set for the data set.

        See Also
        --------
        delete_psf : Delete the PSF model for a data set.
        image_psf : Display the 2D PSF model for a data set in the image viewer.
        load_psf : Create a PSF model.
        plot_psf : Plot the 1D PSF model applied to a data set.
        set_psf : Add a PSF model to a data set.

        Examples
        --------

        Change the size and center of the PSF for th default data set:

        >>> psf = get_psf()
        >>> psf.size = (21,21)
        >>> psf.center = (10,10)

        """
        return self._get_item(id, self._psf, 'psf model', 'has not been set')


    ### Ahelp ingest: 2015-05-11 DJB
    def delete_psf(self, id=None):
        """Delete the PSF model for a data set.

        Remove the PSF convolution applied to a source model.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the
           default identifier is used, as returned by `get_default_id`.

        See Also
        --------
        load_psf : Create a PSF model.
        set_psf : Add a PSF model to a data set.
        get_psf : Return the PSF model defined for a data set.

        Examples
        --------

        >>> delete_psf()

        >>> delete_psf('core')

        """
        id = self._fix_id(id)
        self._psf.pop(id, None)


    def _add_psf(self, id, data, model):
        id = self._fix_id(id)
        psf = self._psf.get(id, None)

        if psf is not None:
            model = psf(model)
            psf.fold(data)
        return model


    #
    # Parameters
    #

    def _check_par(self, par, argname='par'):
        if isinstance(par, basestring):
            par = self._eval_model_expression(par, 'parameter')
        _check_type(par, sherpa.models.Parameter, argname,
                    'a parameter object or parameter expression string')
        return par

    def get_par(self, par):
        """
        get_par

        SYNOPSIS
           Return a model parameter

        SYNTAX

        Arguments:
           par       - model parameter

        Returns:
           Sherpa model parameter

        DESCRIPTION
           Return a Sherpa model parameter given a parameter string.

        SEE ALSO
           set_par
        """
        return self._check_par(par)

    def set_par(self, par, val=None, min=None, max=None, frozen=None):
        """Set the value, limits, or behavior of a model parameter.

        set_par

        SYNOPSIS
           Set initial values for a model parameter

        SYNTAX

        Arguments:
           par       - model parameter

           val       - initial parameter value
                       default = None

           min       - minimum limit
                       default = None

           max       - maximum limit
                       default = None

           frozen    - is the parameter frozen?
                       default = None

        Returns:
           None

        DESCRIPTION
           Set initial values for parameter fields which include initial
           value, minimum limit, maximum limit, and whether it should be
           frozen during a fit.

        SEE ALSO
           get_par
        """
        self._check_par(par).set(val, min, max, frozen)

    def _freeze_thaw_par_or_model(self, par, action):
        if isinstance(par, basestring):
            par = self._eval_model_expression(par, 'parameter or model')

        _check_type(par, (sherpa.models.Parameter, sherpa.models.Model), 'par',
                    'a parameter or model object or expression string')

        if isinstance(par, sherpa.models.Parameter):
            getattr(par, action)()
        else:
            for p in par.pars:
                if not p.alwaysfrozen:
                    getattr(p, action)()

    ### Ahelp ingest: 2015-04-28 DJB
    ### DOC-TODO: is this the best way to document the arguments?
    def freeze(self, *args):
        """Fix model parameters so they are not changed by a fit.

        If called with no arguments, then all parameters
        of models in source expressions are frozen. The
        arguments can be parameters or models (in which case
        all parameters of the model are frozen).

        See Also
        --------
        fit : Fit one or more data sets.
        link : Link a parameter value to an associated value.
        thaw : Allow model parameters to be varied during a fit.
        unlink : Unlink a parameter value.

        Examples
        --------

        Fix the FWHM parameter of the line model (in this case a
        `gauss1d` model) so that it will not be varied in the fit.

        >>> set_source(const1d.bgnd + gauss1d.line)
        >>> line.fwhm = 2.1
        >>> freeze(line.fwhm)
        >>> fit()

        Freeze all parameters of the line model and then re-fit:

        >>> freeze(line)
        >>> fit()

        Freeze the nh parameter of the gal model and the abund
        parameter of the src model:

        >>> freeze(gal.nh, src.abund)

        """
        par = list(args)
        for p in par:
            self._freeze_thaw_par_or_model(p, 'freeze')

    ### Ahelp ingest: 2015-04-28 DJB
    ### DOC-TODO: is this the best way to document the arguments?
    def thaw(self, *args):
        """Allow model parameters to be varied during a fit.

        If called with no arguments, then all parameters
        of models in source expressions are thawed. The
        arguments can be parameters or models (in which case
        all parameters of the model are thawed).

        See Also
        --------
        fit : Fit one or more data sets.
        freeze : Fix model parameters so they are not changed by a fit.
        link : Link a parameter value to an associated value.
        unlink : Unlink a parameter value.

        Examples
        --------

        Ensure that the FWHM parameter of the line model (in this case a
        `gauss1d` model) will be varied in any fit.

        >>> set_source(const1d.bgnd + gauss1d.line)
        >>> thaw(line.fwhm)
        >>> fit()

        Thaw all parameters of the line model and then re-fit:

        >>> thaw(line)
        >>> fit()

        Thaw the nh parameter of the gal model and the abund
        parameter of the src model:

        >>> thaw(gal.nh, src.abund)

        """
        par = list(args)
        for p in par:
            self._freeze_thaw_par_or_model(p, 'thaw')

    ### Ahelp ingest: 2015-04-28 DJB
    def link(self, par, val):
        """Link a parameter to a value.

        A parameter can be linked to another parameter value, or
        function of that value, rather than be an independent value.
        As the linked-to values change, the parameter value will
        change.

        Parameters
        ----------
        par :
           The parameter to link.
        val :
           The value - wihch can be a numeric value or a function
           of other model parameters, to set `par` to.

        See Also
        --------
        freeze : Fix model parameters so they are not changed by a fit.
        set_par : Set a parameter value.
        thaw : Allow model parameters to be varied during a fit.
        unlink : Unlink a parameter value.

        Notes
        -----
        The `link` attribute of the parameter is set to match the
        mathematical expression used for `val`.

        For a parameter value to be varied during a fit, it must be
        part of one of the source expressions involved in the fit.
        So, in the following, the `src1.xpos` parameter will not
        be varied because the `src2` model - from which it takes
        its value - is not included in the source expression of any
        of the data sets being fit.

        >>> set_source(1, gauss1d.src1)
        >>> gauss1d.src2
        >>> link(src1.xpos, src2.xpos)
        >>> fit(1)

        One way to work around this is to include the model but
        with zero signal: for example

        >>> set_source(1, gauss1d.src1 + 0 * gauss1d.src2)

        Examples
        --------

        The `fwhm` parameter of the `g2` model is set to be the same as
        the `fwhm` parameter of the `g1` model.

        >>> link(g2.fwhm, g1.fwhm)

        Fix the `pos` parameter of `g2` to be 2.3 more than the `pos`
        parameter of the `g1` model.

        >>> gauss1d.g1
        >>> gauss1d.g2
        >>> g1.pos = 12.2
        >>> link(g2.pos, g1.pos + 2.3)
        >>> g2.pos.val
        14.5
        >>> g1.pos = 12.1
        >>> g2.pos.val
        14.399999999999999

        """
        par = self._check_par(par)
        if isinstance(val, basestring):
            val = self._eval_model_expression(val, 'parameter link')
        par.link = val

    ### Ahelp ingest: 2015-04-28 DJB
    def unlink(self, par):
        """Unlink a parameter value.

        Remove any parameter link - created by `link` - for the
        parameter. The parameter value is reset to the value it
        had before `link` was called.

        Parameters
        ----------
        par :
           The parameter to unlink. If the parameter is not linked
           then nothing happens.

        See Also
        --------
        freeze : Fix model parameters so they are not changed by a fit.
        link : Link a parameter to a value.
        set_par : Set a parameter value.
        thaw : Allow model parameters to be varied during a fit.

        Examples
        --------

        >>> unlink(bgnd.ampl)

        """
        par = self._check_par(par)
        par.unlink()


    ###########################################################################
    # Fitting
    ###########################################################################


    def _add_extra_data_and_models(self, ids, datasets, models):
        pass


    def _get_fit_ids(self, id, otherids):
        # If id==None, assume it means, simultaneous fit
        # to all data sets that have models assigned to
        # them.  Otherwise, assume that id and otherids
        # as input list the data sets to be fit.
        if id is None:
            all_ids = self.list_data_ids()
            if (len(all_ids) > 0):
                id = all_ids[0]
                del all_ids[0]
                otherids = tuple(all_ids)
            if id is None:
                id = self._fix_id(id)

        ids = [id]
        for nextid in otherids:
            nextid = self._fix_id(nextid)
            if nextid not in ids:
                ids.append(nextid)

        return ids


    def _get_fit_obj(self, datasets, models, estmethod):
        datasets   = tuple(datasets)
        models     = tuple(models)
        if len(datasets) == 1:
            d = datasets[0]
            m = models[0]
        else:
            d = sherpa.data.DataSimulFit('simulfit data', datasets)
            m = sherpa.models.SimulFitModel('simulfit model', models)
	if not self._current_method.name == 'gridsearch':
		if m.is_discrete:
		   raise ModelErr("You are trying to fit a model which has a discrete template model component with a continuous optimization method. Since CIAO4.6 this is not possible anymore. Please use gridsearch as the optimization method and make sure that the 'sequence' option is correctly set, or enable interpolation for the templates you are loading (which is the default behavior).")


        f = sherpa.fit.Fit(d, m, self._current_stat, self._current_method,
                           estmethod, self._current_itermethod)

        return f


    def _prepare_fit(self, id, otherids=()):

        # prep data ids for fitting
        ids = self._get_fit_ids(id, otherids)

        # Gather up lists of data objects and models to fit
        # to them.  Add to lists *only* if there actually is
        # a model to fit.  E.g., if data sets 1 and 2 exist,
        # but only data set 1 has a model, then "fit all" is
        # understood to mean "fit 1".  If both data sets have
        # models, then "fit all" means "fit 1 and 2 together".
        datasets = []
        models = []
        fit_to_ids = []
        for i in ids:
            ds = self.get_data(i)
            mod = None
            if self._models.has_key(i) or self._sources.has_key(i):
                mod = self._get_model(i)

            # The issue with putting a try/catch here is that if an exception
            # is thrown folding a model, it will be swallowed up and the user
            # will be confused why the model is not fit.
            # ex. PSF folding where parameter values are particular
            #
            #
            #try:
            #    mod = self._get_model(i)
            #except:
            #    mod = None
            if mod is not None:
                datasets.append(ds)
                models.append(mod)
                fit_to_ids.append(i)

        # If no data sets have models assigned to them, stop now.
        if len(models) < 1:
            raise IdentifierErr("nomodels")
        
        return fit_to_ids, datasets, models


    def _get_fit(self, id, otherids=(), estmethod=None):
        
        fit_to_ids, datasets, models = self._prepare_fit(id, otherids)

        self._add_extra_data_and_models(fit_to_ids, datasets, models)

        f = self._get_fit_obj(datasets, models, estmethod)

        fit_to_ids = tuple(fit_to_ids)

        return fit_to_ids, f


    def _get_stat_info(self):
	
        ids, datasets, models = self._prepare_fit(None)

        self._add_extra_data_and_models(ids, datasets, models)

        output = []
        if len(datasets) > 1:
            for id, d, m in izip(ids, datasets, models):
                f = sherpa.fit.Fit(d, m, self._current_stat)

                statinfo = f.calc_stat_info()
                statinfo.name = 'Dataset %s' % (str(id))
                statinfo.ids = (id,)
		
                output.append(statinfo)

        f = self._get_fit_obj(datasets, models, None)
        statinfo = f.calc_stat_info()
        if len(ids) == 1:
            statinfo.name = 'Dataset %s' % str(ids)
        else:
            statinfo.name = 'Datasets %s' % str(ids).strip("()")
        statinfo.ids = ids
        output.append(statinfo)

        return output


    ### Ahelp ingest: 2015-05-03 DJB
    def calc_stat_info(self):
        """Display the statistic values for the current models.

        Displays the statistic value for each data set, and the
        combined fit, using the current set of models, parameters, and
        ranges. The output is printed to stdout, and so is intended
        for use in interactive analysis. The `get_stat_info` function
        returns the values as a list of objects.

        See Also
        --------
        calc_stat : Calculate the fit statistic for a data set.
        get_stat_info : Return the statistic values for the current models.

        Notes
        -----
        If a fit to a particular data set has not been made, or values
        - such as parameter settings, the noticed data range, or
        choice of statistic - have been changed since the last fit,
        then the results for that data set may not be meaningful and
        will therefore bias the results for the simultaneous results.

        """
        output = self._get_stat_info()
        output = [statinfo.format() for statinfo in output]

        if len(output) > 1:
            info('\n\n'.join(output))
        else:
            info(output[0])


    ### Ahelp ingest: 2015-05-03 DJB
    def get_stat_info(self):
        """Return the statistic values for the current models.

        Calculate the statistic value for each data set, and the
        combined fit, using the current set of models, parameters, and
        ranges.

        Returns
        -------
        stats : array of sherpa.fit.StatInfoResults
           The values for each data set. If there are multiple model
           expressions then the last element will be the value for the
           combined data sets.

        See Also
        --------
        calc_stat : Calculate the fit statistic for a data set.
        calc_stat_info : Display the statistic values for the current models.
        get_fit_results : Return the results of the last fit.
        list_data_ids : List the identifiers for the loaded data sets.
        list_model_ids : List of all the data sets with a source expression.

        Notes
        -----
        If a fit to a particular data set has not been made, or values
        - such as parameter settings, the noticed data range, or
        choice of statistic - have been changed since the last fit,
        then the results for that data set may not be meaningful and
        will therefore bias the results for the simultaneous results.

        The fields of the object include:

        `name`
           The name of the data set, or sets, as a string.

        `ids`
           A sequence of the data set ids (it may be a tuple or
           array) included in the results.

        `bkg_ids`
           A sequence of the background data set ids (it may be a
           tuple or array) included in the results, if any.

        `statname`
           The name of the statistic function (as used in `set_stat`).

        `statval`
           The statistic value.

        `numpoints`
           The number of bins used in the fits.

        `dof`
           The number of degrees of freedom in the fit (the number of
           bins minus the number of free parameters).

        `qval`
           The Q-value (probability) that one would observe the
           reduced statistic value, or a larger value, if the assumed
           model is true and the current model parameters are the
           true parameter values. This will be `None` if the value
           can not be calculated with the current statistic (e.g.
           the Cash statistic).

        `rstat`
           The reduced statistic value (the `statval` field divided by
           `dof`). This is not calculated for all statistics.

        Examples
        --------

        >>> res = get_stat_info()
        >>> res[0].statval
        498.21750663761935
        >>> res[0].dof
        439

        """
        return self._get_stat_info()


    ### Ahelp ingest: 2015-05-03 DJB
    def get_fit_results(self):
        """Return the results of the last fit.

        Returns
        -------
        stats : sherpa.fit.FitResults instance
           The results of the last fit. It does not reflect any
           changes made to the model parameter, or settings, since the
           last fit.

        See Also
        --------
        calc_stat : Calculate the fit statistic for a data set.
        calc_stat_info : Display the statistic values for the current models.
        fit : Fit a model to one or more data sets.
        get_stat_info : Return the statistic values for the current models.
        set_iter_method : Set the iterative-fitting scheme used in the fit.

        Notes
        -----
        The fields of the object include:

        `datasets`
           A sequence of the data set ids included in the results.

        `itermethodname`
           What iterated-fit scheme was used, if any (as set by
           `set_iter_method`).

        `statname`
           The name of the statistic function (as used in `set_stat`).

        `succeeded`
           Was the fit successful (did it converge)?

        `parnames`
           A tuple of the parameter names that were varied in the fit
           (the thawed parameters in the model expression).

        `parvals`
           A tuple of the parameter values, in the same order as
           `parnames`.

        `statval`
           The statistic value after the fit.

        `istatval`
           The statistic value at the start of the fit.

        `dstatval`
           The change in the statistic value (`istatval - statval`).

        `numpoints`
           The number of bins used in the fits.

        `dof`
           The number of degrees of freedom in the fit (the number of
           bins minus the number of free parameters).

        `qval`
           The Q-value (probability) that one would observe the
           reduced statistic value, or a larger value, if the assumed
           model is true and the current model parameters are the
           true parameter values. This will be `None` if the value
           can not be calculated with the current statistic (e.g.
           the Cash statistic).

        `rstat`
           The reduced statistic value (the `statval` field divided by
           `dof`). This is not calculated for all statistics.

        `message`
           A message about the results of the fit (e.g. if the fit was
           unable to converge). The format and contents depend on the
           optimisation method.

        `nfev`
           The number of model evaluations made during the fit.

        Examples
        --------

        >>> res = get_fit_results()
        >>> res.statval
        498.21750663761935
        >>> res.dof
        439
        >>> res.parnames
        ('pl.gamma', 'pl.ampl', 'gline.fwhm', 'gline.pos', 'gline.ampl')
        >>> res.parvals
        (-0.20659543380329071, 0.00030398852609788524, 100.0, 4900.0, 0.001)

        """
        if self._fit_results == None:
            raise SessionErr('nofit', 'fit')
        else:
            return self._fit_results

    ### Ahelp ingest: 2015-04-28 DJB
    def guess(self, id=None, model=None, limits=True, values=True):
        """Estimate the parameter values and ranges given the loaded data.

        The `guess` function can change the parameter values and
        limits to match the loaded data. This is generally limited to
        changing the amplitude and position parameters (sometimes just
        the values and sometimes just the limits). The parameters that
        are changed depend on the type of model.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        model :
           Change the parameters of this model component. If `None`,
           then the source expression is assumed to consist of a single
           component, and that component is used.
        limits : bool
           Should the parameter limits be changed? The default is `True`.
        values : bool
           Should the parameter values be changed? The default is `True`.

        See Also
        --------
        get_default_id : Return the default data set identifier.
        reset : Reset the model parameters to their default settings.
        set_par : Set a parameter value.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        The `guess` function can reduce the time required to fit a
        data set by moving the parameters closer to a realistic
        solution. It can also be useful because it can set bounds on
        the parameter values based on the data: for instance, many
        two-dimensional models will limit their `xpos` and `ypos`
        values to lie within the data area. This can be done manually,
        but `guess` simplifies this, at least for those parameters
        that are supported. Instrument models - such as an ARF and
        RMF - should be set up *before* calling `guess`.

        Examples
        --------

        Since the source expression contains only one component,
        `guess` can be called with no arguments:

        >>> set_source(polynom1d.poly)
        >>> guess()

        In this case, `guess` is called on each component separately.

        >>> set_source(gauss1d.line + powlaw1d.cont)
        >>> guess(line)
        >>> guess(cont)

        In this example, the values of the `src` model component are
        guessed from the "src" data set, whereas the `bgnd` component
        is guessed from the "bgnd" data set.

        >>> set_source("src", gauss2d.src + const2d.bgnd)
        >>> set_source("bgnd", bgnd)
        >>> guess("src", src)
        >>> guess("bgnd", bgnd)

        """
        if model is None:
            id, model = model, id

        kwargs = { 'limits': limits, 'values': values }

        if model is not None:
            model = self._check_model(model)
            try:

                model.guess(*self.get_data(id).to_guess(), **kwargs)
            except NotImplementedError:
                info('WARNING: No guess found for %s' % model.name)
            return

        ids, f = self._get_fit(self._fix_id(id))
        try:
            f.guess(**kwargs)
        except NotImplementedError:
            #raise NotImplementedError('No guess found for model %s' %
            info('WARNING: No guess found for %s' %
                 self._get_model(id).name)


    ### Ahelp ingest: 2015-05-02 DJB
    def calc_stat(self, id=None, *otherids):
        """Calculate the fit statistic for a data set.

        Evaluate the model for one or more data sets, compare it to
        the data using the current statistic, and return the value.
        No fitting is done, as the current model parameter, and any
        filters, are used.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        otherids : int or str, optional
           Include multiple data sets in the calculation.

        Returns
        -------
        stat : number
           The current statistic value.

        See Also
        --------
        calc_chisqr : Calculate the per-bin chi-squared statistic.
        calc_stat_info : Display the statistic values for the current models.
        set_stat : Set the statistical method.

        Examples
        --------

        Calculate the statistic for the model and data in the default
        data set:

        >>> stat = calc_stat()

        Use the data sets labelled "core" and "jet":

        >>> stat = calc_stat("core", "jet")

        Calculate the statistic value using two different statistics:

        >>> set_stat('cash')
        >>> s1 = calc_stat()
        >>> set_stat('cstat')
        >>> s2 = calc_stat()

        """
        ids, f = self._get_fit(id, otherids)
        return f.calc_stat()

    ### Ahelp ingest: 2015-05-03 DJB
    def calc_chisqr(self, id=None, *otherids):
        """Calculate the per-bin chi-squared statistic.

        Evaluate the model for one or more data sets, compare it to
        the data using the current statistic, and return the value for
        each bin.  No fitting is done, as the current model parameter,
        and any filters, are used.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        otherids : int or str, optional
           Include multiple data sets in the calculation.

        Returns
        -------
        chisq : array or `None`
           The chi-square value for each bin of the data, using the
           current statistic (as set by `set_stat`).  A value of
           `None` is returned if the statistic is not a chi-square
           distribution.

        See Also
        --------
        calc_stat : Calculate the fit statistic for a data set.
        calc_stat_info : Display the statistic values for the current models.
        set_stat : Set the statistical method.

        """
        ids, f = self._get_fit(id, otherids)
        return f.calc_chisqr()

    def fit(self, id=None, *otherids, **kwargs):
        """Fit a model to one or more data sets.

        fit

        SYNOPSIS
           Perform fitting process using current optimization method and 
           current fit statistic.

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - List of other Sherpa data ids

           outfile   - filename and path of parameter value output vs. number
                       of function evaluations
                       default = None

           clobber   - boolean whether to clobber outfile
                       default = False

        Returns:
           Formatted fit results output 

        DESCRIPTION
           Initiate optimization of model parameter values by id(s).

        SEE ALSO
           get_fit_results, conf, proj, covar, set_iter_method,
           set_stat, set_method

        """
        ids, f = self._get_fit(id, otherids)
        res = f.fit(**kwargs)
        res.datasets = ids
        self._fit_results = res
        info(res.format())

    # Back-compatibility
    simulfit = fit


    #
    ## Simulation functions
    #

    def _run_pvalue(self, null_model, alt_model, conv_model=None,
                 id=1, otherids=(), num=500, bins=25, numcores=None):
        ids, fit = self._get_fit(id, otherids)
        
        pvalue = sherpa.sim.LikelihoodRatioTest.run
        results = pvalue(fit, null_model, alt_model, conv_model,
                      niter=num,
                      stat=self._current_stat,
                      method=self._current_method,
                      numcores=numcores)

        info(results.format())
        self._pvalue_results = results


    def get_pvalue_results(self):
        """
        get_pvalue_results

        SYNOPSIS
           Access the simulation results of the likelihood ratio test.

        SYNTAX

        Arguments:
           None

        Returns:
           Likelihood Ratio Test results 

        DESCRIPTION
           Access the simulation results of the likelihood ratio test.

        SEE ALSO
           plot_pvalue, get_pvalue_plot
        """
        return self._pvalue_results


    def plot_pvalue(self, null_model, alt_model, conv_model=None,
                    id=1, otherids=(), num=500, bins=25, numcores=None,
                    replot=False, overplot=False, clearwindow=True):
        """
        plot_pvalue

        SYNOPSIS
           Plot a histogram of likelihood ratios comparing fits of the null
           model to fits of the alternative model to faked data using poisson
           noise.  Computes the likelihood ratio on the real data and the 
           p-value.

        SYNTAX

        Arguments:

           null_model  - model representing the null hypothesis

           alt_model   - alternative model to compare to null

           conv_model  - convolution model to include for fitting.
                         default = None

           id          - Sherpa data id
                         default = default data id

           otherids    - List of other Sherpa data ids
                         default = ()

           num         - Number of iterations to run
                         default = 500

           bins        - Number of bins for the histogram
                         default = 25

           numcores    - Number of cpus to use during simulation
                         default = number of detected cpus

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False


        Returns:
           Likelihood Ratio Test results 

        DESCRIPTION
           Access the simulation results of the likelihood ratio test.

        SEE ALSO
           get_pvalue_results, get_pvalue_plot
        """
        if not sherpa.utils.bool_cast(replot) or self._pvalue_results is None:
            self._run_pvalue(null_model, alt_model, conv_model,
                          id, otherids, num, bins, numcores)

        results = self._pvalue_results
        lrplot = self._lrplot
        if not sherpa.utils.bool_cast(replot):
            lrplot.prepare(results.ratios, bins,
                           num, results.lr, results.ppp)

        try:
            sherpa.plot.begin()
            lrplot.plot(overplot=overplot, clearwindow=clearwindow)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    def get_pvalue_plot(self, null_model=None, alt_model=None, conv_model=None,
                     id=1, otherids=(), num=500, bins=25, numcores=None,
                     recalc=False):
        """
        get_pvalue_plot

        SYNOPSIS
           Acess the histogram plot of the likelihood ratios comparing fits 
           of the null model to fits of the alternative model to faked data
           using poisson noise.  Access the likelihood ratio on the real data
           and the p-value.

        SYNTAX

        Arguments:

           null_model  - model representing the null hypothesis

           alt_model   - alternative model to compare to null

           conv_model  - convolution model to include for fitting.
                         default = None

           id          - Sherpa data id
                         default = default data id

           otherids    - List of other Sherpa data ids
                         default = ()

           num         - Number of iterations to run
                         default = 500

           bins        - Number of bins for the histogram
                         default = 25

           numcores    - Number of cpus to use during simulation
                         default = number of detected cpus

           recalc      - Recalculate the likelihood ratio test simulation
                         default = False

        Returns:
           LRHistogram object

        DESCRIPTION
           Access the histogram plot of the likelihood ratio test.

        SEE ALSO
           plot_pvalue, get_pvalue_results

        """
        lrplot = self._lrplot
        if not recalc:
            return lrplot

        if null_model is None:
            raise TypeError("null model cannot be None")

        if alt_model is None:
            raise TypeError("alternative model cannot be None")

        self._run_pvalue(null_model, alt_model, conv_model,
                      id, otherids, num, bins, numcores)
        results = self._pvalue_results
        lrplot.prepare(results.ratios, bins,
                       num, results.lr, results.ppp)
        self._lrplot = lrplot
        return lrplot


    #
    ## Sampling functions
    #
    ### DOC-TODO: This copies from sherpa/sim/sample.py, so can
    ###           functools/... be used to share docstrings?
    ###           Unfortunately not quite a direct copy, so hard
    ###           to see how to do

    ### Ahelp ingest: 2015-04-30 DJB
    def normal_sample(self, num=1, sigma=1, correlate=True,
                      id=None, otherids=(), numcores=None):
        """Sample the fit statistic by taking the parameter values
        from a normal distribution.

        For each iteration (sample), change the thawed parameters by
        drawing values from a uni- or multi-variate normal (Gaussian)
        distribution, and calculate the fit statistic.

        Parameters
        ----------
        num : int, optional
           The number of samples to use (default is `1`).
        sigma : number, optional
           The width of the normal distribution (the default
           is `1`).
        correlate : bool, optional
           Should a multi-variate normal be used, with parameters
           set by the covariance matrix (`True`) or should a
           uni-variate normal be used (`False`)?
        id : int or str, optional
           The data set containing the model expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        otherids : sequence of int or str, optional
           For when multiple source expressions are being used.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        samples :
           A NumPy array table with the first column representing the
           statistic and later columns the parameters used.

        See Also
        --------
        fit : Fit a model to one or more data sets.
        set_model : Set the source model expression for a data set.
        set_stat : Set the statistical method.
        t_sample : Sample from the Student's t-distribution.
        uniform_sample : Sample from a uniform distribution.

        Notes
        -----
        All thawed model parameters are sampled from the Gaussian
        distribution, where the mean is set as the best-fit parameter
        value and the variance is determined by the diagonal elements
        of the covariance matrix. The multi-variate Gaussian is
        assumed by default for correlated parameters, using the
        off-diagonal elements of the covariance matrix.

        Examples
        --------

        The model fit to the default data set has three free
        parameters. The median value of the statistic calculated by
        `normal_sample` is returned:

        >>> ans = normal_sample(num=10000)
        >>> ans.shape
        (1000, 4)
        >>> np.median(ans[:,0])
        119.82959326927781

        """
        ids, fit = self._get_fit(id, otherids)
        return sherpa.sim.normal_sample(fit, num, sigma, correlate, numcores)


    ### Ahelp ingest: 2015-04-30 DJB
    ### DOC-TODO: improve the description of factor parameter
    def uniform_sample(self, num=1, factor=4,
                       id=None, otherids=(), numcores=None):
        """Sample the fit statistic by taking the parameter values
        from an uniform distribution.

        For each iteration (sample), change the thawed parameters by
        drawing values from a uniform distribution, and calculate the
        fit statistic.

        Parameters
        ----------
        num : int, optional
           The number of samples to use (default is `1`).
        factor : number, optional
           Multiplier to expand the scale parameter (default is `4`).
        id : int or str, optional
           The data set containing the model expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        otherids : sequence of int or str, optional
           For when multiple source expressions are being used.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        samples :
           A NumPy array table with the first column representing the
           statistic and later columns the parameters used.

        See Also
        --------
        fit : Fit a model to one or more data sets.
        normal_sample : Sample from a normal distribution.
        set_model : Set the source model expression for a data set.
        set_stat : Set the statistical method.
        t_sample : Sample from the Student's t-distribution.

        Examples
        --------

        The model fit to the default data set has three free
        parameters. The median value of the statistic calculated by
        `uniform_sample` is returned:

        >>> ans = uniform_sample(num=10000)
        >>> ans.shape
        (1000, 4)
        >>> np.median(ans[:,0])
        284.66534775948134

        """
        ids, fit = self._get_fit(id, otherids)
        return sherpa.sim.uniform_sample(fit, num, factor, numcores)


    ### Ahelp ingest: 2015-04-30 DJB
    def t_sample(self, num=1, dof=None, id=None, otherids=(), numcores=None):
        """Sample the fit statistic by taking the parameter values from
        a Student's t-distribution.

        For each iteration (sample), change the thawed parameters
        by drawing values from a Student's t-distribution, and
        calculate the fit statistic.

        Parameters
        ----------
        num : int, optional
           The number of samples to use (default is `1`).
        dof : optional
           The number of degrees of freedom to use (the default
           is to use the number from the current fit).
        id : int or str, optional
           The data set containing the model expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        otherids : sequence of int or str, optional
           For when multiple source expressions are being used.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        samples :
           A NumPy array table with the first column representing the
           statistic and later columns the parameters used.

        See Also
        --------
        fit : Fit a model to one or more data sets.
        normal_sample : Sample from the normal distribution.
        set_model : Set the source model expression for a data set.
        set_stat : Set the statistical method.
        uniform_sample : Sample from a uniform distribution.

        Examples
        --------

        The model fit to the default data set has three free
        parameters. The median value of the statistic calculated by
        `t_sample` is returned:

        >>> ans = t_sample(num=10000)
        >>> ans.shape
        (1000, 4)
        >>> np.median(ans[:,0])
        119.9764357725326

        """
        ids, fit = self._get_fit(id, otherids)
        if dof is None:
            dof = (len(fit.data.eval_model_to_fit(fit.model)) -
                   len(fit.model.thawedpars))
        return sherpa.sim.t_sample(fit, num, dof, numcores)


    ###########################################################################
    # Error estimation
    ###########################################################################


    ### Ahelp ingest: 2015-05-08 DJB
    ### DOC-TODO: how best to document the settings?
    ### DOC-TODO: have I got soft_limits described correctly?
    def get_covar(self):
        """Return the covariance estimation object.

        Returns
        -------
        covar : object

        See Also
        --------
        covar : Estimate confidence intervals using the covariance method.
        get_covar_opt : Return one or all of the options for the covariance method.
        set_covar_opt : Set an option of the covar estimation object.

        Notes
        -----
        The attributes of the covariance object include:

        `eps`
           The precision of the calculated limits. The default is
           `0.01`.

        `maxiters`
           The maximum number of iterations allowed before stopping
           for that parameter. The default is `200`.

        `sigma`
           What is the error limit being calculated. The default is
           `1`.

        `soft_limits`
           Should the search be restricted to the soft limits of the
           parameters (`True`), or can parameter values go out all the
           way to the hard limits if necessary (`False`).  The default
           is `False`

        Examples
        --------

        >>> print(get_covar())
        name        = covariance
        sigma       = 1
        maxiters    = 200
        soft_limits = False
        eps         = 0.01

        Change the `sigma` field of the `covar` method to 1.9.

        >>> cv = get_covar()
        >>> cv.sigma = 1.6

        """
        return self._estmethods['covariance']

    ### Ahelp ingest: 2015-04-27 DJB
    ### DOC-TODO: how best to document the settings?
    ### DOC-TODO: have I got soft_limits described correctly?
    ### DOC-TODO: when verbose=True how is extra output displayed?
    ###           stdout, stderr, sherpa logging?
    def get_conf(self):
        """Return the confidence-interval estimation object.

        Returns
        -------
        conf : object
           
        See Also
        --------
        conf : Estimate confidence intervals using the confidence method.
        get_conf_opt : Return one or all of the options for the confidence interval method.
        set_conf_opt : Set an option of the conf estimation object.

        Notes
        -----
        The attributes of the confidence-interval object include:

        `eps`
           The precision of the calculated limits. The default is
           `0.01`.
        
        `fast`
           If `True` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is `False`.

        `max_rstat`
           If the reduced chi square is larger than this value, do not
           use (only used with chi-square statistics). The default is
           `3`.

        `maxfits`
           The maximum number of re-fits allowed (that is, when the
           `remin` filter is met). The default is `5`.

        `maxiters`
           The maximum number of iterations allowed when bracketing
           limits, before stopping for that parameter. The default is
           `200`.

        `numcores`
           The number of computer cores to use when evaluating results
           in parallel. This is only used if `parallel` is `True`.
           The default is to use all cores.

        `openinterval`
           How the `conf` method should cope with intervals that do
           not converge (that is, when the `maxiters` limit has been
           reached). The default is `False`.

        `parallel`
           If there is more than one free parameter then the results
           can be evaluated in parallel, to reduce the time required.
           The default is `True`.

        `remin`
           The minimum difference in statistic value for a new fit
           location to be considered better than the current best fit
           (which starts out as the starting location of the fit at
           the time `conf` is called). The default is `0.01`.

        `sigma`
           What is the error limit being calculated. The default is
           `1`.

        `soft_limits`
           Should the search be restricted to the soft limits of the
           parameters (`True`), or can parameter values go out all the
           way to the hard limits if necessary (`False`).  The default
           is `False`

        `tol`
           The tolerance for the fit. The default is `0.2`.

        `verbose`
           Should extra information be displayed during fitting?
           The default is `False`.

        Examples
        --------

        >>> print(get_conf())
        name         = confidence
        numcores     = 8
        verbose      = False
        openinterval = False
        max_rstat    = 3
        maxiters     = 200
        soft_limits  = False
        eps          = 0.01
        fast         = False
        maxfits      = 5
        remin        = 0.01
        tol          = 0.2
        sigma        = 1
        parallel     = True

        Change the `remin` field of the `conf` method to 0.05.

        >>> cf = get_conf()
        >>> cf.remin = 0.05

        """
        return self._estmethods['confidence']

    ### Ahelp ingest: 2015-05-07 DJB
    def get_proj(self):
        """Return the confidence-interval estimation object.

        .. note:: The `conf` function should be used instead of `proj`.

        Returns
        -------
        proj : object

        See Also
        --------
        conf : Estimate confidence intervals for fit parameters.
        get_proj_opt : Return one or all of the options for the confidence interval method.
        proj : Estimate confidence intervals for fit parameters.
        set_proj_opt : Set an option of the proj estimation object.

        Notes
        -----
        The attributes of the object include:

        `eps`
           The precision of the calculated limits. The default is
           `0.01`.

        `fast`
           If `True` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is `False`.

        `max_rstat`
           If the reduced chi square is larger than this value, do not
           use (only used with chi-square statistics). The default is
           `3`.

        `maxfits`
           The maximum number of re-fits allowed (that is, when the
           `remin` filter is met). The default is `5`.

        `maxiters`
           The maximum number of iterations allowed when bracketing
           limits, before stopping for that parameter. The default is
           `200`.

        `numcores`
           The number of computer cores to use when evaluating results
           in parallel. This is only used if `parallel` is `True`.
           The default is to use all cores.

        `parallel`
           If there is more than one free parameter then the results
           can be evaluated in parallel, to reduce the time required.
           The default is `True`.

        `remin`
           The minimum difference in statistic value for a new fit
           location to be considered better than the current best fit
           (which starts out as the starting location of the fit at
           the time `proj` is called). The default is `0.01`.

        `sigma`
           What is the error limit being calculated. The default is
           `1`.

        `soft_limits`
           Should the search be restricted to the soft limits of the
           parameters (`True`), or can parameter values go out all the
           way to the hard limits if necessary (`False`).  The default
           is `False`

        `tol`
           The tolerance for the fit. The default is `0.2`.

        Examples
        --------

        >>> print(get_proj())
        name        = projection
        numcores    = 8
        max_rstat   = 3
        maxiters    = 200
        soft_limits = False
        eps         = 0.01
        fast        = False
        maxfits     = 5
        remin       = 0.01
        tol         = 0.2
        sigma       = 1
        parallel    = True

        """
        return self._estmethods['projection']

    # New "wrappers" to access estimation methods without calling
    # get_proj(), etc.

    def _check_estmethod_opt(self, estmethod, optname):
        _check_type(optname, basestring, 'optname', 'a string')
        if optname not in estmethod.config:
            raise ArgumentErr('badopt', optname, estmethod.name)

    def _get_estmethod_opt(self, methodname, optname=None):
        meth = self._estmethods.get(methodname.lower())
        if meth is None:
            raise ArgumentErr('badconf', methodname)
        if optname is None:
            return meth.config
        self._check_estmethod_opt(meth,optname)
        return meth.config[optname]

    def _set_estmethod_opt(self, methodname, optname, val):
        meth = self._estmethods.get(methodname.lower())
        if meth is None:
            raise ArgumentErr('badconf', methodname)
        self._check_estmethod_opt(meth, optname)
        meth.config[optname] = val
        
    ### Ahelp ingest: 2015-05-08 DJB
    def get_covar_opt(self, name=None):
        """Return one or all of the options for the covariance
        method.

        This is a helper function since the options can also
        be read directly using the object returned by `get_covar`.

        Parameters
        ----------
        name : str, optional
           If not given, a dictionary of all the options are returned.
           When given, the individual value is returned.

        Returns
        -------
        value : dictionary or value

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `name` argument is not recognized.

        See Also
        --------
        covar : Estimate confidence intervals using the covariance method.
        get_covar : Return the covariance estimation object.
        set_covar_opt : Set an option of the covar estimation object.

        Examples
        --------

        >>> get_covar_opt('sigma')
        1

        >>> copts = get_covar_opt()
        >>> copts['sigma']
        1

        """
        return self._get_estmethod_opt('covariance', name)

    ### Ahelp ingest: 2015-04-27 DJB
    def get_conf_opt(self, name=None):
        """Return one or all of the options for the confidence interval
        method.

        This is a helper function since the options can also
        be read directly using the object returned by `get_conf`.

        Parameters
        ----------
        name : str, optional
           If not given, a dictionary of all the options are returned.
           When given, the individual value is returned.
            
        Returns
        -------
        value : dictionary or value
           
        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `name` argument is not recognized.

        See Also
        --------
        conf : Estimate confidence intervals using the confidence method.
        get_conf : Return the confidence-interval estimation object.
        set_conf_opt : Set an option of the conf estimation object.

        Examples
        --------

        >>> get_conf_opt('verbose')
        False

        >>> copts = get_conf_opt()
        >>> copts['verbose']
        False

        """
        return self._get_estmethod_opt('confidence', name)
    
    ### Ahelp ingest: 2015-05-07 DJB
    def get_proj_opt(self, name=None):
        """Return one or all of the options for the confidence interval
        method.

        .. note:: The `conf` function should be used instead of `proj`.

        This is a helper function since the options can also
        be read directly using the object returned by `get_proj`.

        Parameters
        ----------
        name : str, optional
           If not given, a dictionary of all the options are returned.
           When given, the individual value is returned.

        Returns
        -------
        value : dictionary or value

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `name` argument is not recognized.

        See Also
        --------
        conf : Estimate confidence intervals for fit parameters.
        proj : Estimate confidence intervals for fit parameters.
        get_proj : Return the confidence-interval estimation object.
        set_proj_opt : Set an option of the proj estimation object.

        Examples
        --------

        >>> get_proj_opt('sigma')
        1

        >>> popts = get_proj_opt()
        >>> popts['sigma']
        1

        """
        return self._get_estmethod_opt('projection', name)
    
    ### Ahelp ingest: 2015-05-08 DJB
    def set_covar_opt(self, name, val):
        """Set an option for the covariance method.

        This is a helper function since the options can also
        be set directly using the object returned by `get_covar`.

        Parameters
        ----------
        name : str
           The name of the option to set. The `get_covar`
           routine can be used to find out valid values for
           this argument.
        val :
           The new value for the option.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `name` argument is not recognized.

        See Also
        --------
        covar : Estimate confidence intervals using the covariance method.
        get_covar : Return the covar estimation object.
        get_covar_opt : Return one or all options of the covar estimation object.

        Examples
        --------

        >>> set_covar_opt('sigma', 1.6)

        """
        self._set_estmethod_opt('covariance', name, val)

    ### Ahelp ingest: 2015-04-27 DJB
    def set_conf_opt(self, name, val):
        """Set an option for the confidence interval method.

        This is a helper function since the options can also
        be set directly using the object returned by `get_conf`.

        Parameters
        ----------
        name : str
           The name of the option to set. The `get_conf`
           routine can be used to find out valid values for
           this argument.
        val :
           The new value for the option.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `name` argument is not recognized.

        See Also
        --------
        conf : Estimate confidence intervals using the confidence method.
        get_conf : Return the conf estimation object.
        get_conf_opt : Return one or all options of the conf estimation object.

        Examples
        --------

        >>> set_conf_opt('parallel', False)

        """
        self._set_estmethod_opt('confidence', name, val)

    ### Ahelp ingest: 2015-05-08 DJB
    def set_proj_opt(self, name, val):
        """Set an option for the projection method.

        .. note:: The `conf` function should be used instead of `proj`.

        This is a helper function since the options can also
        be set directly using the object returned by `get_proj`.

        Parameters
        ----------
        name : str
           The name of the option to set. The `get_proj`
           routine can be used to find out valid values for
           this argument.
        val :
           The new value for the option.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `name` argument is not recognized.

        See Also
        --------
        conf : Estimate confidence intervals using the confidence method.
        proj : Estimate confidence intervals using the projection method.
        get_proj : Return the proj estimation object.
        get_proj_opt : Return one or all options of the proj estimation object.

        Examples
        --------

        >>> set_proj_opt('parallel', False)

        """
        self._set_estmethod_opt('projection', name, val)

    ### Ahelp ingest: 2015-05-08 DJB
    def get_covar_results(self):
        """Return the results of the last `covar` run.

        Returns
        -------
        results : sherpa.fit.ErrorEstResults object

        Raises
        ------
        sherpa.utils.err.SessionErr
           If no `covar` call has been made.

        See Also
        --------
        get_covar_opt : Return one or all of the options for the covariance method.
        set_covar_opt : Set an option of the covar estimation object.

        Notes
        -----
        The fields of the object include:

        `datasets`
           A tuple of the data sets used in the analysis.

        `methodname`
           This will be `covariance`.

        `iterfitname`
           The name of the iterated-fit method used, if any.

        `fitname`
           The name of the optimization method used.

        `statname`
           The name of the fit statistic used.

        `sigma`
           The `sigma` value used to calculate the confidence
           intervals.

        `percent`
           The percentage of the signal contained within the
           confidence intervals (calculated from the `sigma`
           value assuming a normal distribution).

        `parnames`
           A tuple of the parameter names included in the analysis.

        `parvals`
           A tuple of the best-fit parameter values, in the same
           order as `parnames`.

        `parmins`
           A tuple of the lower error bounds, in the same
           order as `parnames`.

        `parmaxes`
           A tuple of the upper error bounds, in the same
           order as `parnames`.

        `nfits`

        There is also an `extra_output` field which is used to return
        the covariance matrix.

        Examples
        --------

        >>> res = get_covar_results()
        >>> print(res)
        datasets    = (1,)
        methodname  = covariance
        iterfitname = none
        fitname     = levmar
        statname    = chi2gehrels
        sigma       = 1
        percent     = 68.2689492137
        parnames    = ('bgnd.c0',)
        parvals     = (10.228675427602724,)
        parmins     = (-2.4896739438296795,)
        parmaxes    = (2.4896739438296795,)
        nfits       = 0

        In this case, of a single parameter, the covariance
        matrix is just the variance of the parameter:

        >>> copt.extra_output
        array([[ 6.19847635]])

        """
        if self._covariance_results == None:
            raise SessionErr('noaction', "covariance")
        else:
            return self._covariance_results

    ### Ahelp ingest: 2015-04-27 DJB
    ### DOC_TODO: what is the best description for nfits?
    def get_conf_results(self):
        """Return the results of the last `conf` run.

        Returns
        -------
        results : sherpa.fit.ErrorEstResults object
           
        Raises
        ------
        sherpa.utils.err.SessionErr
           If no `conf` call has been made.

        See Also
        --------
        get_conf_opt : Return one or all of the options for the confidence interval method.
        set_conf_opt : Set an option of the conf estimation object.

        Notes
        -----
        The fields of the object include:

        `datasets`
           A tuple of the data sets used in the analysis.

        `methodname`
           This will be `confidence`.

        `iterfitname`
           The name of the iterated-fit method used, if any.

        `fitname`
           The name of the optimization method used.

        `statname`
           The name of the fit statistic used.

        `sigma`
           The `sigma` value used to calculate the confidence
           intervals.

        `percent`
           The percentage of the signal contained within the
           confidence intervals (calculated from the `sigma`
           value assuming a normal distribution).

        `parnames`
           A tuple of the parameter names included in the analysis.

        `parvals`
           A tuple of the best-fit parameter values, in the same
           order as `parnames`.

        `parmins`
           A tuple of the lower error bounds, in the same
           order as `parnames`.

        `parmaxes`
           A tuple of the upper error bounds, in the same
           order as `parnames`.

        `nfits`


        Examples
        --------

        >>> res = get_conf_results()
        >>> print(res)
        datasets    = (1,)
        methodname  = confidence
        iterfitname = none
        fitname     = levmar
        statname    = chi2gehrels
        sigma       = 1
        percent     = 68.2689492137
        parnames    = ('p1.gamma', 'p1.ampl')
        parvals     = (2.1585155113403327, 0.00022484014787994827)
        parmins     = (-0.082785567348122591, -1.4825550342799376e-05)
        parmaxes    = (0.083410634144100104, 1.4825550342799376e-05)
        nfits       = 13

        The following converts the above into a dictionary where the
        keys are the parameter names and the values are the tuple
        (best-fit value, lower-limit, upper-limit):
        
        >>> pvals1 = zip(res.parvals, res.parmins, res.parmaxes)
        >>> pvals2 = [(v, v+l, v+h) for (v,l,h) in pvals1]
        >>> dres = dict(zip(res.parnames, pvals2))
        >>> dres['p1.gamma']
        (2.1585155113403327, 2.07572994399221, 2.241926145484433)

        """
        if self._confidence_results == None:
            raise SessionErr('noaction', "confidence")
        else:
            return self._confidence_results

    ### Ahelp ingest: 2015-05-07 DJB
    def get_proj_results(self):
        """Return the results of the last `proj` run.

        .. note:: The `conf` function should be used instead of `proj`.

        Returns
        -------
        results : sherpa.fit.ErrorEstResults object

        Raises
        ------
        sherpa.utils.err.SessionErr
           If no `proj` call has been made.

        See Also
        --------
        conf : Estimate confidence intervals for fit parameters.
        proj : Estimate confidence intervals for fit parameters.
        get_proj_opt : Return one or all of the options for the projection method.
        set_proj_opt : Set an option of the proj estimation object.

        Notes
        -----
        The fields of the object include:

        `datasets`
           A tuple of the data sets used in the analysis.

        `methodname`
           This will be `projection`.

        `iterfitname`
           The name of the iterated-fit method used, if any.

        `fitname`
           The name of the optimization method used.

        `statname`
           The name of the fit statistic used.

        `sigma`
           The `sigma` value used to calculate the confidence
           intervals.

        `percent`
           The percentage of the signal contained within the
           confidence intervals (calculated from the `sigma`
           value assuming a normal distribution).

        `parnames`
           A tuple of the parameter names included in the analysis.

        `parvals`
           A tuple of the best-fit parameter values, in the same
           order as `parnames`.

        `parmins`
           A tuple of the lower error bounds, in the same
           order as `parnames`.

        `parmaxes`
           A tuple of the upper error bounds, in the same
           order as `parnames`.

        `nfits`

        Examples
        --------

        >>> res = get_proj_results()
        >>> print(res)
        datasets    = ('src',)
        methodname  = projection
        iterfitname = none
        fitname     = levmar
        statname    = chi2gehrels
        sigma       = 1
        percent     = 68.2689492137
        parnames    = ('bgnd.c0',)
        parvals     = (9.1958148476800918,)
        parmins     = (-2.0765029551804268,)
        parmaxes    = (2.0765029551935186,)
        nfits       = 0

        """
        if self._projection_results == None:
            raise SessionErr('noaction', "projection")
        else:
            return self._projection_results

    def _est_errors(self, args, methodname):
        # Any argument that is a model parameter should be detected
        # and added to the list of parameters for which we want limits.
        # Else, the argument is an integer or string denoting a
        # data set ID.
        id = None
        parlist = []
        otherids = ()
        for arg in args:
            if type(arg) is sherpa.models.Parameter:
                if arg.frozen is False:
                    parlist.append(arg)
                else:
                    raise sherpa.utils.err.ParameterErr('frozen', arg.fullname)
            else:
                if id is None:
                    id = arg
                else:
                    otherids = otherids + (arg,)
        if (len(parlist) == 0):
            parlist = None
        ids, f = self._get_fit(id, otherids, self._estmethods[methodname])
        res = f.est_errors(self._methods, parlist)
        res.datasets = ids
        info(res.format())
        return res

    ### Ahelp ingest: 2015-05-08 DJB
    ### DOC-TODO: include screen output of covar() ?
    def covar(self, *args):
        """Estimate the confidence intervals for parameters using the
        covariance method.

        The `covar` command computes confidence interval bounds for
        the specified model parameters in the dataset, using the
        covariance matrix of the statistic. The `get_covar` and
        `set_covar_opt` commands can be used to configure the error
        analysis; an example being changing the 'sigma' field to `1.6`
        (i.e. 90%) from its default value of `1`.  The output from the
        routine is displayed on screen, and the `get_covar_results`
        routine can be used to retrieve the results.

        Parameters
        ----------
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        parameters : optional
           The default is to calculate the confidence limits on all
           thawed parameters of the model, or models, for all the
           data sets. The evaluation can be restricted by listing
           the parameters to use. Note that each parameter should be
           given as a separate argument, rather than as a list.
           For example `covar(g1.ampl, g1.sigma)`.

        See Also
        --------
        covar : Estimate the confidence intervals using the confidence method.
        get_covar : Return the covariance estimation object.
        get_covar_results : Return the results of the last `covar` run.
        int_proj : Plot the statisitic value as a single parameter is varied.
        int_unc : Plot the statisitic value as a single parameter is varied.
        reg_proj : Plot the statistic value as two parameters are varied.
        reg_unc : Plot the statistic value as two parameters are varied.
        set_covar_opt : Set an option of the `covar` estimation object.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with multiple `ids` or `parameters` values, the
        order is unimportant, since any argument that is not defined
        as a model parameter is assumed to be a data id.

        The `covar` command is different to `conf`, in that in that
        all other thawed parameters are fixed, rather than being
        allowed to float to new best-fit values.  While `conf` is more
        general (e.g. allowing the user to examine the parameter space
        away from the best-fit point), it is in the strictest sense no
        more accurate than `covar` for determining confidence
        intervals.

        An estimated confidence interval is accurate if and only if:

        1. the chi^2 or logL surface in parameter space is
           approximately shaped like a multi-dimensional paraboloid,
           and

        2. the best-fit point is sufficiently far from parameter space
           boundaries.

        One may determine if these conditions hold, for example, by
        plotting the fit statistic as a function of each parameter's
        values (the curve should approximate a parabola) and by
        examining contour plots of the fit statistics made by varying
        the values of two parameters at a time (the contours should be
        elliptical, and parameter space boundaries should be no closer
        than approximately 3 sigma from the best-fit point). The
        `int_proj` and `reg_proj` commands may be used for this.

        If either of the conditions given above does not hold, then
        the output from `covar` may be meaningless except to give an
        idea of the scale of the confidence intervals. To accurately
        determine the confidence intervals, one would have to
        reparameterize the model, use Monte Carlo simulations, or
        Bayesian methods.

        As `covar` estimates intervals for each parameter
        independently, the relationship between sigma and the change
        in statistic value delta_S can be particularly simple: sigma =
        the square root of delta_S for statistics sampled from the
        chi-square distribution and for the Cash statistic, and is
        approximately equal to the square root of (2 * delta_S) for
        fits based on the general log-likelihood. The default setting
        is to calculate the one-sigma interval, which can be changed
        with the `sigma` option to `set_covar_opt` or `get_covar`.

        Examples
        --------

        Evaluate confidence intervals for all thawed parameters in all
        data sets with an associated source model. The results are
        then stored in the variable `res`.

        >>> covar()
        >>> res = get_covar_results()

        Only evaluate the parametes associated with data set 2.

        >>> covar(2)

        Only evaluate the intervals for the `pos.xpos` and `pos.ypos`
        parameters:

        >>> covar(pos.xpos, pos.ypos)

        Change the limits to be 1.6 sigma (90%) rather than the default
        1 sigma.

        >>> get_covar().sigma = 1.6
        >>> covar()

        Only evaluate the `clus.kt` parameter for the data sets with
        identifiers "obs1", "obs5", and "obs6". This will still use
        the 1.6 sigma setting from the previous run.

        >>> covar("obs1", ["obs5","obs6"], clus.kt)

        """
        self._covariance_results = self._est_errors(args, 'covariance')
    
    ### Ahelp ingest: 2015-04-27 DJB
    ### DOC-TODO: include screen output of conf() ?
    def conf(self, *args):
        """Estimate the confidence intervals for parameters using the
        confidence method.

        The `conf` command computes confidence interval bounds for the
        specified model parameters in the dataset.  A given
        parameter's value is varied along a grid of values while the
        values of all the other thawed parameters are allowed to float
        to new best-fit values. The `get_conf` and `set_conf_opt`
        commands can be used to configure the error analysis; an
        example being changing the 'sigma' field to `1.6` (i.e. 90%)
        from its default value of `1`.  The output from the routine is
        displayed on screen, and the `get_conf_results` routine can be
        used to retrieve the results.

        Parameters
        ----------
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        parameters : optional
           The default is to calculate the confidence limits on all
           thawed parameters of the model, or models, for all the
           data sets. The evaluation can be restricted by listing
           the parameters to use. Note that each parameter should be
           given as a separate argument, rather than as a list.
           For example `conf(g1.ampl, g1.sigma)`.

        See Also
        --------
        covar : Estimate the confidence intervals using the covariance method.
        get_conf : Return the confidence-interval estimation object.
        get_conf_results : Return the results of the last `conf` run.
        int_proj : Plot the statisitic value as a single parameter is varied.
        int_unc : Plot the statisitic value as a single parameter is varied.
        reg_proj : Plot the statistic value as two parameters are varied.
        reg_unc : Plot the statistic value as two parameters are varied.
        set_conf_opt : Set an option of the `conf` estimation object.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with multiple `ids` or `parameters` values, the
        order is unimportant, since any argument that is not defined
        as a model parameter is assumed to be a data id.

        The `conf` command is different to `covar`, in that in that
        all other thawed parameters are allowed to float to new
        best-fit values, instead of being fixed to the initial
        best-fit values as they are in `covar`.  While `conf` is more
        general (e.g. allowing the user to examine the parameter space
        away from the best-fit point), it is in the strictest sense no
        more accurate than `covar` for determining confidence
        intervals.

        An estimated confidence interval is accurate if and only if:

        1. the chi^2 or logL surface in parameter space is
           approximately shaped like a multi-dimensional paraboloid,
           and

        2. the best-fit point is sufficiently far from parameter space
           boundaries.

        One may determine if these conditions hold, for example, by
        plotting the fit statistic as a function of each parameter's
        values (the curve should approximate a parabola) and by
        examining contour plots of the fit statistics made by varying
        the values of two parameters at a time (the contours should be
        elliptical, and parameter space boundaries should be no closer
        than approximately 3 sigma from the best-fit point). The
        `int_proj` and `reg_proj` commands may be used for this.

        If either of the conditions given above does not hold, then
        the output from `conf` may be meaningless except to give an
        idea of the scale of the confidence intervals. To accurately
        determine the confidence intervals, one would have to
        reparameterize the model, use Monte Carlo simulations, or
        Bayesian methods.

        As the calculation can be computer intensive, the default
        behavior is to use all available CPU cores to speed up the
        analysis. This can be changed be varying the `numcores` option
        - or setting `parallel` to `False` - either with
        `set_conf_opt` or `get_conf`.

        As `conf` estimates intervals for each parameter
        independently, the relationship between sigma and the change
        in statistic value delta_S can be particularly simple: sigma =
        the square root of delta_S for statistics sampled from the
        chi-square distribution and for the Cash statistic, and is
        approximately equal to the square root of (2 * delta_S) for
        fits based on the general log-likelihood. The default setting
        is to calculate the one-sigma interval, which can be changed
        with the `sigma` option to `set_conf_opt` or `get_conf`.

        The limit calculated by `conf` is basically a 1-dimensional
        root in the translated coordinate system (translated by the
        value of the statistic at the minimum plus sigma^2).  The
        Taylor series expansion of the multi-dimensional function at
        the minimum is:

        f(x + dx) ~ f(x) + grad( f(x) )^T dx + dx^T Hessian( f(x) ) dx + ...

        where x is understood to be the n-dimensional vector
        representing the free parameters to be fitted and the
        super-script 'T' is the transpose of the row-vector. At or
        near the minimum, the gradient of the function is zero or
        negligible, respectively. So the leading term of the expansion
        is quadratic.  The best root finding algorithm for a curve
        which is approximately parabolic is Muller's method [1]_.
        Muller's method is a generalization of the secant method [2]_:
        the secant method is an iterative root finding method that
        approximates the function by a straight line through two
        points, whereas Muller's method is an iterative root finding
        method that approxmiates the function by a quadratic
        polynomial through three points.

        Three data points are the minimum input to Muller's root
        finding method. The first point to be submitted to the
        Muller's root finding method is the point at the minimum. To
        strategically choose the other two data points, the confidence
        function uses the output from covariance as the second data
        point. To generate the third data points for the input to
        Muller's root finding method, the secant root finding method
        is used since it only requires two data points to generate the
        next best approximation of the root.

        However, there are cases where `conf` cannot locate the root
        even though the root is bracketed within an interval (perhaps
        due to the bad resolution of the data). In such cases, when
        the option `openinterval` is set to `False` (which is the
        default), the routine will print a warning message about not
        able to find the root within the set tolerance and the
        function will return the average of the open interval which
        brackets the root. If `openinterval` is set to `True` then
        `conf` will print the minimal open interval which brackets the
        root (not to be confused with the lower and upper bound of the
        confidence interval). The most accurate thing to do is to
        return an open interval where the root is localized/bracketed
        rather then the average of the open interval (since the
        average of the interval is not a root within the specified
        tolerance).

        References
        ----------

        .. [1] Muller, David E., "A Method for Solving Algebraic
               Equations Using an Automatic Computer," MTAC, 10
               (1956), 208-215.

        .. [2] Numerical Recipes in Fortran, 2nd edition, 1986, Press
               et al., p. 347

        Examples
        --------

        Evaluate confidence intervals for all thawed parameters in all
        data sets with an associated source model. The results are
        then stored in the variable `res`.

        >>> conf()
        >>> res = get_conf_results()

        Only evaluate the parametes associated with data set 2.

        >>> conf(2)

        Only evaluate the intervals for the `pos.xpos` and `pos.ypos`
        parameters:

        >>> conf(pos.xpos, pos.ypos)

        Change the limits to be 1.6 sigma (90%) rather than the default
        1 sigma.

        >>> get_conf().sigma = 1.6
        >>> conf()

        Only evaluate the `clus.kt` parameter for the data sets with
        identifiers "obs1", "obs5", and "obs6". This will still use
        the 1.6 sigma setting from the previous run.

        >>> conf("obs1", ["obs5","obs6"], clus.kt)

        """
        self._confidence_results = self._est_errors(args, 'confidence')

    ### Ahelp ingest: 2015-05-06 DJB
    ### DOC-TODO: add a deprecation note?
    def proj(self, *args):
        """Estimate the confidence intervals for parameters using the
        projection method.

        .. note:: The `conf` function should be used instead of `proj`.

        The `proj` command computes confidence interval bounds for the
        specified model parameters in the dataset. A given parameter's
        value is varied along a grid of values while the values of all
        the other thawed parameters are allowed to float to new
        best-fit values. The `get_proj` and `set_proj_opt`
        commands can be used to configure the error analysis; an
        example being changing the 'sigma' field to `1.6` (i.e. 90%)
        from its default value of `1`.  The output from the routine is
        displayed on screen, and the `get_proj_results` routine can be
        used to retrieve the results.

        Parameters
        ----------
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        parameters : optional
           The default is to calculate the confidence limits on all
           thawed parameters of the model, or models, for all the
           data sets. The evaluation can be restricted by listing
           the parameters to use. Note that each parameter should be
           given as a separate argument, rather than as a list.
           For example `proj(g1.ampl, g1.sigma)`.

        See Also
        --------
        conf : Estimate the confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        get_proj : Return the confidence-interval estimation object.
        get_proj_results : Return the results of the last `proj` run.
        int_proj : Plot the statisitic value as a single parameter is varied.
        reg_proj : Plot the statistic value as two parameters are varied.
        set_proj_opt : Set an option of the `proj` estimation object.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with multiple `ids` or `parameters` values, the
        order is unimportant, since any argument that is not defined
        as a model parameter is assumed to be a data id.

        The `proj` command is different to `covar`, in that all other
        thawed parameters are allowed to float to new best-fit values,
        instead of being fixed to the initial best-fit values.  While
        `proj` is more general (e.g. allowing the user to examine the
        parameter space away from the best-fit point), it is in the
        strictest sense no more accurate than `covar` for
        determining confidence intervals.

        An estimated confidence interval is accurate if and only if:

        1. the chi^2 or logL surface in parameter space is
           approximately shaped like a multi-dimensional paraboloid,
           and

        2. the best-fit point is sufficiently far from parameter space
           boundaries.

        One may determine if these conditions hold, for example, by
        plotting the fit statistic as a function of each parameter's
        values (the curve should approximate a parabola) and by
        examining contour plots of the fit statistics made by varying
        the values of two parameters at a time (the contours should be
        elliptical, and parameter space boundaries should be no closer
        than approximately 3 sigma from the best-fit point). The
        `int_proj` and `reg_proj` commands may be used for this.

        If either of the conditions given above does not hold, then
        the output from `proj` may be meaningless except to give an
        idea of the scale of the confidence intervals. To accurately
        determine the confidence intervals, one would have to
        reparameterize the model, use Monte Carlo simulations, or
        Bayesian methods.

        As the calculation can be computer intensive, the default
        behavior is to use all available CPU cores to speed up the
        analysis. This can be changed be varying the `numcores` option
        - or setting `parallel` to `False` - either with
        `set_proj_opt` or `get_proj`.

        As `proj` estimates intervals for each parameter
        independently, the relationship between sigma and the change
        in statistic value delta_S can be particularly simple: sigma =
        the square root of delta_S for statistics sampled from the
        chi-square distribution and for the Cash statistic, and is
        approximately equal to the square root of (2 * delta_S) for
        fits based on the general log-likelihood. The default setting
        is to calculate the one-sigma interval, which can be changed
        with the `sigma` option to `set_proj_opt` or `get_proj`.

        """
        self._projection_results = self._est_errors(args, 'projection')

    # Aliases
    get_covariance_results = get_covar_results
    get_confidence_results = get_conf_results
    get_projection_results = get_proj_results
    covariance = covar
    confidence = conf
    projection = proj


    ###########################################################################
    # PyBLoCXS routines for Markov Chain Monte Carlo
    #
    # DOC-TODO: should this use functools.wraps or something similar,
    #           to avoid copying the docs?
    # DOC-TODO: integrate the existing pyblocks python documentation - e.g.
    #           http://hea-www.harvard.edu/AstroStat/pyBLoCXS/#high-level-user-interface-functions
    ###########################################################################


    ### Ahelp ingest: 2015-04-30 DJB
    def set_sampler_opt(self, opt, value):
        """Set an option for the current pyBLoCXS sampler.

        Parameters
        ----------
        opt : str
           The option to change. Use `get_sampler` to view the
           available options for the current sampler.
        value :
           The value for the option.

        See Also
        --------
        get_sampler : Return the current pyBLoCXS sampler options.
        set_prior: Set the prior function to use with a parameter.
        set_sampler : Set the pyBLoCXS sampler.

        Notes
        -----
        The options depend on the sampler [1]_. The options include:

        `defaultprior`
           Set to `False` when the default prior (flat, between the
           parameter's soft limits) should not be used. Use
           `set_prior` to set the form of the prior for each
           parameter.

        `inv`
           A bool, or array of bools, to indicate which parameter is
           on the inverse scale.

        `log`
           A bool, or array of bools, to indicate which parameter is
           on the logarithm (natural log) scale.

        `original`
           A bool, or array of bools, to indicate which parameter is
           on the original scale.

        `p_M`
           The proportion of jumps generatd by the Metropolis
           jumping rule.

        `priorshape`
           An array of bools indicating which parameters have a
           user-defined prior functions set with `set_prior`.

        `scale`
           Multiply the output of `covar` by this factor and
           use the result as the scale of the t-distribution.

        References
        ----------

        .. [1] http://hea-www.harvard.edu/AstroStat/pyBLoCXS/#high-level-user-interface-functions

        Examples
        --------

        >>> set_sampler_opt('scale', 3)

        """
        self._pyblocxs.set_sampler_opt(opt, value)

    ### Ahelp ingest: 2015-04-30 DJB
    def get_sampler_opt(self, opt):
        """Return an option of the current pyBLoCXS sampler.

        Returns
        -------
        opt : str
           The name of the option. The fields depend on the current
           sampler.

        See Also
        --------
        get_sampler : Return the current pyBLoCXS sampler options.
        set_sampler_opt : Set an option for the current pyBLoCXS sampler.

        Examples
        --------

        >>> get_sampler_opt('log')
        False

        """
        return self._pyblocxs.get_sampler_opt(opt)

    ### Ahelp ingest: 2015-04-30 DJB
    def get_sampler_name(self):
        """Return the name of the current pyBLoCXS sampler.

        Returns
        -------
        name : str

        See Also
        --------
        get_sampler : Return the current pyBLoCXS sampler options.
        set_sampler : Set the pyBLoCXS sampler.

        Examples
        --------

        >>> get_sampler_name()
        'MetropolisMH'

        """
        return self._pyblocxs.get_sampler_name()
        
    ### Ahelp ingest: 2015-04-30 DJB
    def set_sampler(self, sampler):
        """Set the pyBLoCXS sampler.

        The sampler determines the type of jumping rule to
        be used when running the MCMC analysis.

        Parameters
        ----------
        sampler : str or sherpa.sim.Sampler instance
           When a string, the name of the sampler to use (case
           insensitive). The supported options are given by the
           `list_samplers` function.

        See Also
        --------
        get_draws : Run the pyBLoCXS MCMC algorithm.
        list_samplers : List the pyBLoCXS samplers.
        set_sampler : Set the pyBLoCXS sampler.
        set_sampler_opt : Set an option for the current pyBLoCXS sampler.

        Notes
        -----
        The jumping rules are [1]_:

        MH
           The 'MH' option refers to the Metropolis-Hastings rule,
           which always jumps from the best-fit location.

        MetropolisMH
           This is the Metropolis with Metropolis-Hastings algorithm,
           that jumps from the best-fit with probability 'p_M',
           otherwise it jumps from the last accepted jump. The
           value of `p_M` can be changed using `set_sampler_opt`.

        PragBayes
           This is used when the effective area calibration
           uncertainty is to be included in the calculation. At each
           nominal MCMC iteration, a new calibration product is
           generated, and a series of N (the `nsubiters` option) MCMC
           sub-iteration steps are carried out, choosing between
           Metropolis and Metropolis-Hastings types of samplers with
           probability `p_M`.  Only the last of these sub-iterations
           are kept in the chain.  The `nsubiters` and `p_M` values
           can be changed using `set_sampler_opt`.

        FullBayes
           Another sampler for use when including uncertainties due
           to the effective area.

        References
        ----------

        .. [1] http://hea-www.harvard.edu/AstroStat/pyBLoCXS/#high-level-user-interface-functions

        Examples
        --------

        >>> set_sampler('metropolismh')

        """
        self._pyblocxs.set_sampler(sampler)

    ### Ahelp ingest: 2015-04-30 DJB
    def get_sampler(self):
        """Return the current pyBLoCXS sampler options.

        Returns
        -------
        options : dict
           A copy of the  options for the chosen sampler.  Use
           `set_sampler_opt` to change these values. The fields depend
           on the current sampler.

        See Also
        --------
        get_sampler_name : Return the name of the current pyBLoCXS sampler.
        get_sampler_opt : Return an option of the current pyBLoCXS sampler.
        set_sampler : Set the pyBLoCXS sampler.
        set_sampler_opt : Set an option for the current pyBLoCXS sampler.

        """
        return self._pyblocxs.get_sampler()

    ### Ahelp ingest: 2015-04-30 DJB
    ### DOC-TODO: should set_sampler_opt be mentioned here?
    def set_prior(self, par, prior):
        """Set the prior function to use with a parameter.

        The pyBLoCXS Markov Chain Monte Carlo (MCMC) algorithm [1]_
        supports Bayesian Low-Count X-ray Spectral analysis. By
        default, a flat prior is used for each parameter in the fit,
        varying between its soft minimum and maximum values.  The
        `set_prior` function is used to change the form of the prior
        for a parameter.

        Parameters
        ----------
        par : sherpa.models.parameter.Parameter instance
           A parameter of a model instance.
        prior : function or sherpa.models.model.Model instance
           The function to use for a prior. It must accept a
           single argument and return a value of the same size
           as the input.

        See Also
        --------
        get_draws : Run the pyBLoCXS MCMC algorithm.
        get_prior : Set the prior function to use with a parameter.
        set_sampler : Set the pyBLoCXS sampler.

        References
        ----------

        .. [1] http://hea-www.harvard.edu/AstroStat/pyBLoCXS/#high-level-user-interface-functions

        Examples
        --------

        Set the prior for the `kT` parameter of the `therm` instance
        to be a gaussian, centered on 1.7 keV and with a FWHM of 0.35
        keV:

        >>> create_model_component('xsapec', 'therm')
        >>> create_model_component('gauss1d', 'p_temp')
        >>> p_temp.pos = 1.7
        >>> p.temo_fwhm = 0.35
        >>> set_prior(therm.kT, p_temp)

        Create a function (`lognorm`) and use it as the prior the the
        `nH` parameter of the `abs1` instance:

        >>> create_model_component('xsphabs', 'abs1')
        >>> def lognorm(x):
           # center on 10^20 cm^2 with a sigma of 0.5
           sigma = 0.5
           x0 = 20
           # nH is in units of 10^-22 so convert
           dx = np.log10(x) + 22 - x0
           norm = sigma / np.sqrt(2 * np.pi)
           return norm * np.exp(-0.5*dx*dx/(sigma*sigma))

        >>> set_prior(abs1.nH, lognorm)

        """
        self._pyblocxs.set_prior(par, prior)

    ### Ahelp ingest: 2015-04-30 DJB
    def get_prior(self, par):
        """Return the prior function for a parameter.

        Parameters
        ----------
        par : sherpa.models.parameter.Parameter
           A parameter of a model instance.

        Returns
        -------
        prior :
           The function or parameter instance set by
           a previous call to `set_prior`.

        Raises
        ------
        ValueError
           If a prior has not been set for the parameter.

        See Also
        --------
        set_prior : Set the prior function to use with a parameter.

        Examples
        --------

        >>> pfunc = get_prior(bgnd.c0)

        """
        return self._pyblocxs.get_prior(par)

    ### Ahelp ingest: 2015-04-30 DJB
    ### DOC-TODO: include examples once this returns something useful
    def list_priors(self):
        """Return the priors set for model parameters, if any.

        Returns
        -------
        priors : string
           A string representation of the dictionary mapping between
           parameters (keys) and priot functions (values).

        See Also
        --------
        get_prior : Return the prior function for a parameter.
        set_prior : Set the prior function to use with a parameter.

        """
        return self._pyblocxs.list_priors()

    ### Ahelp ingest: 2015-04-30 DJB
    def list_samplers(self):
        """List the pyBLoCXS samplers.

        Returns
        -------
        samplers : list of str
           A list of the names (in lower case) that can be used with
           `set_sampler`.

        See Also
        --------
        get_sampler_name : Return the name of the current pyBLoCXS sampler.
        set_sampler : Set the pyBLoCXS sampler.

        Examples
        --------

        >>> list_samplers()
        ['metropolismh', 'fullbayes', 'mh', 'pragbayes']

        """
        return self._pyblocxs.list_samplers()

    def get_draws(self, id=None, otherids=(), niter=1000):

        ids, fit = self._get_fit(id, otherids)

        # Allow the user to jump from a user defined point in parameter space?
        # Meaning let the user set up parameter space without fitting first.

        #fit_results = self.get_fit_results()
        #if fit_results is None:
        #    raise TypeError("Fit has not been run")

        covar_results = self.get_covar_results()
        if covar_results is None:
            raise TypeError("Covariance has not been calculated")

        covar_matrix = covar_results.extra_output

        stats, accept, params = self._pyblocxs.get_draws(fit, covar_matrix, niter=niter)
        return (stats, accept, params)


    ###########################################################################
    # Basic plotting
    ###########################################################################
    #
    # Plot object access
    #

    def get_split_plot(self):
        """
        get_split_plot

        SYNOPSIS
           Return a Sherpa split plot

        SYNTAX

        Arguments:
           None

        Returns:
           Sherpa SplitPlot plot

        DESCRIPTION
           The Sherpa split plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              rows         - number of rows of plots
                             default = 2

              cols         - number of columns of plots
                             default = 1

              plot_prefs   - dictionary of plotting preferences
                  None

           Functions:

              addplot(self, plot, *args, **kwargs)
                 add a plot to the series in the split plot panel

              addcontour(self, plot, *args, **kwargs)
                 add a contour plot to the series in the split plot panel

              plot(self, row, col, plot, *args, **kwargs)
                 send the split plot panel to the visualizer

              contour(self, row, col, plot, *args, **kwargs)
                 send the split plot panel to the visualizer

              overlayplot(self, plot, *args, **kwargs)
                 plot over current plot

              overlaycontour(self, plot, *args, **kwargs)
                 plot contour over current contour plot

              overplot(self, row, col, plot, *args, **kwargs)
                 plot over current plot at specific coordinates

              overcontour(self, row, col, plot, *args, **kwargs)
                plot contour over current contour plot at specific coordinates

        SEE ALSO
           plot, plot_fit_resid, plot_fit_delchi
        """
        return self._splitplot

    ### Ahelp ingest: 2015-04-29 DJB
    def get_data_plot(self, id=None):
        """Return the data used by plot_data.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        data : a sherpa.plot.DataPlot instance
           An object representing the data used to create the plot by
           `plot_data`. The relationship between the returned values
           and the values in the data set depend on the data type. For
           example PHA data are plotted in units controlled by
           `set_analysis`, but are stored as channels and counts, and
           may have been grouped and the background estimate removed.

        See Also
        --------
        get_data_plot_prefs : Return the preferences for plot_data.
        get_default_id : Return the default data set identifier.
        plot_data : Plot the data values.

        """
        self._prepare_plotobj(id, self._dataplot)
        return self._dataplot

    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: discussion of preferences needs better handling
    ###           of how it interacts with the chosen plot backend.
    def get_data_plot_prefs(self):
        """Return the preferences for plot_data.

        Returns
        -------
        prefs : dict
           Changing the values of this dictionary will change any new
           data plots. This dictionary will be empty if no plot
           backend is available.

        See Also
        --------
        plot_data : Plot the data values.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Notes
        -----
        The meaning of the fields depend on the chosen plot backend.
        A value of `None` means to use the default value for that
        attribute, unless indicated otherwise. These preferences
        are used by the following commands: `plot_data`, `plot_bkg`,
        `plot_ratio`, and the "fit" variants, such as `plot_fit`,
        `plot_fit_resid`, and `plot_bkg_fit`.

        `errcolor`
           The color to draw error bars. The default is `None`.

        `errstyle`
           How to draw errors. The default is `line`.

        `errthickness`
           What thickness of line to draw error bars. The default is
           `None`.

        `linecolor`
           What color to use for the line connecting the data points.
           The default is `None`.

        `linestyle`
           How should the line connecting the data points be drawn.
           The default is `0`, which means no line is drawn.

        `linethickness`
           What thickness should be used to draw the line connecting
           the data points. The default is `None`.

        `ratioline`
           Should a horizontal line be drawn at y=1?  The default is
           `False`.

        `symbolcolor`
           What color to draw the symbol representing the data points.
           The default is `None`.

        `symbolfill`
           Should the symbol be drawn filled? The default is `False`.

        `symbolsize`
           What size is the symbol drawn. The default is `3`.

        `symbolstyle`
           What style is used for the symbols. The default is `4`
           which means `circle` for the ChIPS back end.

        `xaxis`
           The default is `False`

        `xerrorbars`
           Should error bars be drawn for the X axis. The default is
           `False`.

        `xlog`
           Should the X axis be drawn with a logarithmic scale? The
           default is `False`. This field can also be changed with the
           `set_xlog` and `set_xlinear` functions.

        `yerrorbars`
           Should error bars be drawn for the Y axis. The default is
           `True`.

        `ylog`
           Should the Y axis be drawn with a logarithmic scale? The
           default is `False`. This field can also be changed with the
           `set_ylog` and `set_ylinear` functions.

        Examples
        --------

        After these commands, any data plot will use a green symbol
        and not display Y error bars.

        >>> prefs = get_data_plot_prefs()
        >>> prefs['symbolcolor'] = 'green'
        >>> prefs['yerrorbars'] = False

        """
        return self._dataplot.plot_prefs

    # also in sherpa.astro.utils (copies this docstring)
    ### Ahelp ingest: 2015-05-12 DJB
    def get_model_plot(self, id=None):
        """Return the data used by plot_model.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        data :
           An object representing the data used to create the plot by
           `plot_model`. The return value depends on the data
           set (e.g. 1D binned or un-binned).

        See Also
        --------
        get_model_plot_prefs : Return the preferences for plot_model.
        plot_model : Plot the model for a data set.

        Examples
        --------

        >>> mplot = get_model_plot()

        """
        self._prepare_plotobj(id, self._modelplot)
        return self._modelplot

    # also in sherpa.astro.utils (does not copy this docstring)
    ### Ahelp ingest: 2015-05-12 DJB
    def get_source_plot(self, id=None):
        """Return the data used by plot_source.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        data :
           An object representing the data used to create the plot by
           `plot_source`. The return value depends on the data
           set (e.g. 1D binned or un-binned).

        See Also
        --------
        get_model_plot : Return the data used by plot_model.
        plot_model : Plot the model for a data set.
        plot_source : Plot the source expression for a data set.

        Examples
        --------

        >>> splot = get_source_plot()

        >>> splot1 = get_source_plot(id='jet')
        >>> splot2 = get_source_plot(id='core')

        """
        self._prepare_plotobj(id, self._sourceplot)
        return self._sourceplot


    # sherpa.astro.utils version copies this docstring
    def get_model_component_plot(self, id, model=None):
        """Return the data used by plot_model_component.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to use (the name, if a string).

        Returns
        -------
        data :
           An object representing the data used to create the plot by
           `plot_model_component`. The return value depends on the
           data set (e.g. 1D binned or un-binned).

        See Also
        --------
        get_model_plot : Return the data used by plot_model.
        plot_model : Plot the model for a data set.
        plot_model_component : Plot a component of the model for a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        Examples
        --------

        Return the plot data for the `pl` component used in the
        default data set:

        >>> cplot = get_model_component(pl)

        Return the full source model (`fplot`) and then for the
        components `gal * pl` and `gal * gline`, for the data set
        'jet':

        >>> fmodel = xsphabs.gal * (powlaw1d.pl + gauss1d.gline)
        >>> set_source('jet', fmodel)
        >>> fit('jet')
        >>> fplot = get_model('jet')
        >>> plot1 = get_model_component('jet', pl*gal)
        >>> plot2 = get_model_component('jet', gline*gal)

        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        self._prepare_plotobj(id, self._compmdlplot, model=model)
        return self._compmdlplot


    # sherpa.astro.utils version copies this docstring
    def get_source_component_plot(self, id, model=None):
        """Return the data used by plot_source_component.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to use (the name, if a string).

        Returns
        -------
        data :
           An object representing the data used to create the plot by
           `plot_source_component`. The return value depends on the
           data set (e.g. 1D binned or un-binned).

        See Also
        --------
        get_source_plot : Return the data used by plot_source.
        plot_source : Plot the source expression for a data set.
        plot_source_component : Plot a component of the source expression for a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        Examples
        --------

        Return the plot data for the `pl` component used in the
        default data set:

        >>> cplot = get_source_component(pl)

        Return the full source model (`fplot`) and then for the
        components `gal * pl` and `gal * gline`, for the data set
        'jet':

        >>> fmodel = xsphabs.gal * (powlaw1d.pl + gauss1d.gline)
        >>> set_source('jet', fmodel)
        >>> fit('jet')
        >>> fplot = get_source('jet')
        >>> plot1 = get_source_component('jet', pl*gal)
        >>> plot2 = get_source_component('jet', gline*gal)

        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        if isinstance(model, sherpa.models.TemplateModel):
            self._prepare_plotobj(id, self._comptmplsrcplot, model=model)
            return self._comptmplsrcplot

        self._prepare_plotobj(id, self._compsrcplot, model=model)
        return self._compsrcplot


    ### Ahelp ingest: 2015-05-12 DJB
    def get_model_plot_prefs(self):
        """Return the preferences for plot_model.

        Returns
        -------
        prefs : dict
           Changing the values of this dictionary will change any new
           model plots. This dictionary will be empty if no plot
           backend is available.

        See Also
        --------
        plot_model : Plot the model for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Notes
        -----
        The meaning of the fields depend on the chosen plot backend.
        A value of `None` means to use the default value for that
        attribute, unless indicated otherwise. These preferences are
        used by the following commands: `plot_model`, `plot_ratio`,
        `plot_bkg_model`, and the "fit" variants, such as `plot_fit`,
        `plot_fit_resid`, and `plot_bkg_fit`.

        `errcolor`
           The color to draw error bars. The default is `None`.

        `errstyle`
           How to draw errors. The default is `None`.

        `errthickness`
           What thickness of line to draw error bars. The default is
           `None`.

        `linecolor`
           What color to use for the line connecting the data points.
           The default is `red`.

        `linestyle`
           How should the line connecting the data points be drawn.
           The default is `1`, which means a solid line is drawn.

        `linethickness`
           What thickness should be used to draw the line connecting
           the data points. The default is `3`.

        `ratioline`
           Should a horizontal line be drawn at y=1?  The default is
           `False`.

        `symbolcolor`
           What color to draw the symbol representing the data points.
           The default is `None`.

        `symbolfill`
           Should the symbol be drawn filled? The default is `True`.

        `symbolsize`
           What size is the symbol drawn. The default is `None`.

        `symbolstyle`
           What style is used for the symbols. The default is `0`,
           which means no symbol is used.

        `xaxis`
           The default is `False`

        `xerrorbars`
           Should error bars be drawn for the X axis. The default is
           `False`.

        `xlog`
           Should the X axis be drawn with a logarithmic scale? The
           default is `False`. This field can also be changed with the
           `set_xlog` and `set_xlinear` functions.

        `yerrorbars`
           Should error bars be drawn for the Y axis. The default is
           `False`.

        `ylog`
           Should the Y axis be drawn with a logarithmic scale? The
           default is `False`. This field can also be changed with the
           `set_ylog` and `set_ylinear` functions.

        Examples
        --------

        After these commands, any model plot will use a green line
        to display the model:

        >>> prefs = get_model_plot_prefs()
        >>> prefs['linecolor'] = 'green'

        """
        return self._modelplot.plot_prefs

    def get_fit_plot(self, id=None):
        """
        get_fit_plot

        SYNOPSIS
           Return a Sherpa fit plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa FitPlot plot

        DESCRIPTION
           The Sherpa fit plot object holds a reference to a data plot and
           model plot instance.

           Attributes:
              dataplot

              modelplot

        SEE ALSO
           plot_fit
        """
        self._prepare_plotobj(id, self._fitplot)
        return self._fitplot

    ### Ahelp ingest: 2015-05-11 DJB
    def get_resid_plot(self, id=None):
        """Return the data used by plot_resid.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        resid_data : a sherpa.plot.ResidPlot instance

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_chisqr_plot : Return the data used by plot_chisqr.
        get_delchi_plot : Return the data used by plot_delchi.
        get_ratio_plot : Return the data used by plot_ratio.
        plot_resid : Plot the residuals (data - model) for a data set.

        Examples
        --------

        Return the residual data for the default data set:

        >>> rplot = get_resid_plot()
        >>> np.min(rplot.y)
        -2.9102595936209896
        >>> np.max(rplot.y)
        4.0897404063790104

        """
        self._prepare_plotobj(id, self._residplot)
        return self._residplot

    ### Ahelp ingest: 2015-05-11 DJB
    def get_delchi_plot(self,id=None):
        """Return the data used by plot_delchi.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        resid_data : a sherpa.plot.DelchiPlot instance

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_chisqr_plot : Return the data used by plot_chisqr.
        get_ratio_plot : Return the data used by plot_ratio.
        get_resid_plot : Return the data used by plot_resid.
        plot_delchi : Plot the ratio of residuals to error for a data set.

        Examples
        --------

        Return the residual data, measured in units of the error, for
        the default data set:

        >>> rplot = get_delchi_plot()
        >>> np.min(rplot.y)
        -2.85648373819671875
        >>> np.max(rplot.y)
        2.89477053577520982

        """
        self._prepare_plotobj(id, self._delchiplot)
        return self._delchiplot

    ### Ahelp ingest: 2015-05-11 DJB
    def get_chisqr_plot(self,id=None):
        """Return the data used by plot_chisqr.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        resid_data : a sherpa.plot.ChisqrPlot instance

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_delchi_plot : Return the data used by plot_delchi.
        get_ratio_plot : Return the data used by plot_ratio.
        get_resid_plot : Return the data used by plot_resid.
        plot_chisqr : Plot the chi-squared value for each point in a data set.

        Examples
        --------

        Return the residual data, measured as chi square, for the
        default data set:

        >>> rplot = get_chisqr_plot()
        >>> np.min(rplot.y)
        0.0005140622701128954
        >>> np.max(rplot.y)
        8.379696454792295

        """
        self._prepare_plotobj(id, self._chisqrplot)
        return self._chisqrplot
    
    ### Ahelp ingest: 2015-05-11 DJB
    def get_ratio_plot(self, id=None):
        """Return the data used by plot_ratio.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        resid_data : a sherpa.plot.RatioPlot instance

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_chisqr_plot : Return the data used by plot_chisqr.
        get_delchi_plot : Return the data used by plot_delchi.
        get_resid_plot : Return the data used by plot_resid.
        plot_ratio : Plot the ratio of data to model for a data set.

        Examples
        --------

        Return the ratio of the data to the model for the default data
        set:

        >>> rplot = get_ratio_plot()
        >>> np.min(rplot.y)
        0.6320905073750186
        >>> np.max(rplot.y)
        1.5170172177000447

        """
        self._prepare_plotobj(id, self._ratioplot)
        return self._ratioplot
    
    ### Ahelp ingest: 2015-05-12 DJB
    def get_data_contour(self, id=None):
        """Return the data used by contour_data.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        resid_data : a sherpa.plot.DataContour instance
           The `y` attribute contains the residual values and the `x0`
           and `x1` arrays the corresponsing coordinate values, as
           one-dimensional arrays.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set is not 2D.
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_data_image : Return the data used by image_data.
        contour_data : Contour the values of an image data set.
        image_data : Display a data set in the image viewer.

        Examples
        --------

        Return the data for the default data set:

        >>> dinfo = get_data_contour()

        """
        self._prepare_plotobj(id, self._datacontour)
        return self._datacontour

    ### Ahelp ingest: 2015-05-12 DJB
    def get_data_contour_prefs(self):
        """Return the preferences for contour_data.

        Returns
        -------
        prefs : dict
           Changing the values of this dictionary will change any new
           contour plots. The default is an empty dictionary.

        See Also
        --------
        contour_data : Contour the values of an image data set.

        Notes
        -----
        The meaning of the fields depend on the chosen plot backend.
        A value of `None` (or not set) means to use the default value
        for that attribute, unless indicated otherwise.

        `color`
           The color to draw the contours. The default is `None`.

        `style`
           How to draw the contours. The default is `None`.

        `thickness`
           What thickness of line to draw the contours. The default is
           `None`.

        `xlog`
           Should the X axis be drawn with a logarithmic scale? The
           default is `False`.

        `ylog`
           Should the Y axis be drawn with a logarithmic scale? The
           default is `False`.

        Examples
        --------

        Change the contours to be drawn in 'green':

        >>> contour_data()
        >>> prefs = get_data_contour_prefs()
        >>> prefs['color'] = 'green'
        >>> contour_data()

        """
        return self._datacontour.contour_prefs

    ### Ahelp ingest: 2015-05-12 DJB
    def get_model_contour(self, id=None):
        """Return the data used by contour_model.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        resid_data : a sherpa.plot.ModelContour instance
           The `y` attribute contains the model values and the `x0`
           and `x1` arrays the corresponsing coordinate values, as
           one-dimensional arrays.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set is not 2D.
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_model_image : Return the data used by image_model.
        contour_model : Contour the values of the model, including any PSF.
        image_model : Display the model for a data set in the image viewer.

        Examples
        --------

        Return the model pixel values for the default data set:

        >>> minfo = get_model_contour()

        """
        self._prepare_plotobj(id, self._modelcontour)
        return self._modelcontour

    ### Ahelp ingest: 2015-05-12 DJB
    def get_source_contour(self, id=None):
        """Return the data used by contour_source.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        resid_data : a sherpa.plot.SourceContour instance
           The `y` attribute contains the model values and the `x0`
           and `x1` arrays the corresponsing coordinate values, as
           one-dimensional arrays.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set is not 2D.
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_source_image : Return the data used by image_source.
        contour_source : Contour the values of the model, without any PSF.
        image_source : Display the source expression for a data set in the image viewer.

        Examples
        --------

        Return the source model pixel values for the default data set:

        >>> sinfo = get_source_contour()

        """
        self._prepare_plotobj(id, self._sourcecontour)
        return self._sourcecontour

    ### Ahelp ingest: 2015-05-12 DJB
    def get_model_contour_prefs(self):
        """Return the preferences for contour_model.

        Returns
        -------
        prefs : dict
           Changing the values of this dictionary will change any new
           contour plots.

        See Also
        --------
        contour_model : Contour the values of the model, including any PSF.

        Notes
        -----
        The meaning of the fields depend on the chosen plot backend.
        A value of `None` (or not set) means to use the default value
        for that attribute, unless indicated otherwise.

        `color`
           The color to draw the contours. The default is `red`.

        `style`
           How to draw the contours. The default is `None`.

        `thickness`
           What thickness of line to draw the contours. The default is
           `3`.

        `xlog`
           Should the X axis be drawn with a logarithmic scale? The
           default is `False`.

        `ylog`
           Should the Y axis be drawn with a logarithmic scale? The
           default is `False`.

        Examples
        --------

        Change the contours for the model to be drawn in 'orange':

        >>> prefs = get_model_contour_prefs()
        >>> prefs['color'] = 'orange'
        >>> contour_data()
        >>> contour_model(overcontour=True)

        """
        return self._modelcontour.contour_prefs

    def get_fit_contour(self, id=None):
        """
        get_fit_contour

        SYNOPSIS
           Return a Sherpa fit contour

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa FitContour plot

        DESCRIPTION
           The Sherpa fit contour object holds a reference to a data contour
           and model contour instance.

           Attributes:
              datacontour

              modelcontour

        SEE ALSO
           contour_fit
        """
        self._prepare_plotobj(id, self._fitcontour)
        return self._fitcontour

    ### Ahelp ingest: 2015-05-11 DJB
    def get_resid_contour(self, id=None):
        """Return the data used by contour_resid.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        resid_data : a sherpa.plot.ResidContour instance
           The `y` attribute contains the residual values and the `x0`
           and `x1` arrays the corresponsing coordinate values, as
           one-dimensional arrays.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set is not 2D.
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_ratio_contour : Return the data used by contour_ratio.
        get_resid_image : Return the data used by image_resid.
        contour_resid : Contour the residuals of the fit.
        image_resid : Display the residuals (data - model) for a data set in the image viewer.

        Examples
        --------

        Return the residual data for the default data set:

        >>> rinfo = get_resid_contour()

        """
        self._prepare_plotobj(id, self._residcontour)
        return self._residcontour

    ### Ahelp ingest: 2015-05-11 DJB
    def get_ratio_contour(self, id=None):
        """Return the data used by contour_ratio.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        ratio_data : a sherpa.plot.RatioContour instance
           The `y` attribute contains the ratio values and the `x0`
           and `x1` arrays the corresponsing coordinate values, as
           one-dimensional arrays.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set is not 2D.
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_ratio_image : Return the data used by image_ratio.
        get_resid_contour : Return the data used by contour_resid.
        contour_ratio : Contour the ratio of data to model.
        image_ratio : Display the ratio (data/model) for a data set in the image viewer.

        Examples
        --------

        Return the ratio data for the default data set:

        >>> rinfo = get_ratio_contour()

        """
        self._prepare_plotobj(id, self._ratiocontour)
        return self._ratiocontour

    ### Ahelp ingest: 2015-05-11 DJB
    def get_psf_contour(self, id=None):
        """Return the data used by contour_psf.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        psf_data : a sherpa.plot.PSFContour instance

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set is not 2D.
        sherpa.utils.err.IdentifierErr
           If a PSF model has not been created for the data set.

        See Also
        --------
        get_kernel_contour : Return the data used by contour_kernel.
        contour_kernel : Contour the kernel applied to the model of an image data set.
        contour_psf : Contour the PSF applied to the model of an image data set.

        Examples
        --------

        Return the contour data for the PSF for the default data set:

        >>> cplot = get_psf_contour()

        """
        self._prepare_plotobj(id, self._psfcontour)
        return self._psfcontour


    ### Ahelp ingest: 2015-05-11 DJB
    def get_kernel_contour(self, id=None):
        """Return the data used by contour_kernel.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        psf_data : a sherpa.plot.PSFKernelContour instance

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set is not 2D.
        sherpa.utils.err.IdentifierErr
           If a PSF model has not been created for the data set.

        See Also
        --------
        get_psf_contour : Return the data used by contour_psf.
        contour_kernel : Contour the kernel applied to the model of an image data set.
        contour_psf : Contour the PSF applied to the model of an image data set.

        Examples
        --------

        Return the contour data for the kernel for the default data
        set:

        >>> kplot = get_kernel_contour()

        """
        self._prepare_plotobj(id, self._kernelcontour)
        return self._kernelcontour


    ### Ahelp ingest: 2015-05-11 DJB
    def get_psf_plot(self, id=None):
        """Return the data used by plot_psf.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        psf_plot : a sherpa.plot.PSFPLot instance

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If a PSF model has not been created for the data set.

        See Also
        --------
        get_kernel_plot : Return the data used by plot_kernel.
        plot_kernel : Plot the 1D kernel applied to a data set.
        plot_psf : Plot the 1D PSF model applied to a data set.

        Examples
        --------

        Return the plot data and then create a plot with it:

        >>> pplot = get_psf_plot()
        >>> pplot.plot()

        """
        self._prepare_plotobj(id, self._psfplot)
        return self._psfplot


    ### Ahelp ingest: 2015-05-11 DJB
    def get_kernel_plot(self, id=None):
        """Return the data used by plot_kernel.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        kernel_plot : a sherpa.plot.PSFKernelPLot instance

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If a PSF model has not been created for the data set.

        See Also
        --------
        get_psf_plot : Return the data used by plot_psf.
        plot_kernel : Plot the 1D kernel applied to a data set.
        plot_psf : Plot the 1D PSF model applied to a data set.

        Examples
        --------

        Return the plot data and then create a plot with it:

        >>> kplot = get_kernel_plot()
        >>> kplot.plot()

        """
        self._prepare_plotobj(id, self._kernelplot)
        return self._kernelplot

    #
    # Line plots
    #

    def _prepare_plotobj(self, id, plotobj, model=None):
        id = self._fix_id(id)
        if isinstance(plotobj, sherpa.plot.FitPlot):
            plotobj.prepare(self._prepare_plotobj(id, self._dataplot),
                            self._prepare_plotobj(id, self._modelplot))
        elif isinstance(plotobj, sherpa.plot.FitContour):
            plotobj.prepare(self._prepare_plotobj(id, self._datacontour),
                            self._prepare_plotobj(id, self._modelcontour))
        else:
            if(isinstance(plotobj, sherpa.plot.PSFPlot) or
               isinstance(plotobj, sherpa.plot.PSFContour) or
               isinstance(plotobj, sherpa.plot.PSFKernelPlot) or 
               isinstance(plotobj, sherpa.plot.PSFKernelContour)):
                plotobj.prepare(self.get_psf(id), self.get_data(id))
            elif(isinstance(plotobj, sherpa.plot.DataPlot) or
                 isinstance(plotobj, sherpa.plot.DataContour)):
                plotobj.prepare(self.get_data(id), self.get_stat())
            elif(isinstance(plotobj, sherpa.plot.ComponentModelPlot) or 
                 isinstance(plotobj, sherpa.plot.ComponentSourcePlot)):
                plotobj.prepare(self.get_data(id), model, self.get_stat())
            elif(isinstance(plotobj, sherpa.plot.SourcePlot) or
                 isinstance(plotobj, sherpa.plot.SourceContour)):
                plotobj.prepare(self.get_data(id), self._get_source(id),
                                self.get_stat())
            else:
                # Using _get_fit becomes very complicated using simulfit
                # models and datasets
                #
                #ids, f = self._get_fit(id)
                plotobj.prepare(self.get_data(id), self.get_model(id),
                                self.get_stat())

        return plotobj


    def _multi_plot(self, args, plotmeth='plot'):
        if len(args) == 0:
            raise ArgumentTypeErr('plotargs')

        plots = []
        type = getattr(self, '_%s_types' % plotmeth)
        args = list(args)
        while args:
            plottype = args.pop(0)
            _check_type(plottype, basestring, 'plottype', 'a string')
            plottype = plottype.lower()

            if plottype not in type:
                raise ArgumentErr('badplottype', plottype)

            if args and (args[0] not in type):
                id = args.pop(0)
            else:
                id = None

            plots.append((id, type[plottype]))

        if len(plots) == 1:
            plotmeth = getattr(self, '_' + plotmeth)
            plotmeth(*plots[0])
            return

        nrows = 2
        ncols = int((len(plots) + 1) / 2.0)
        sp = self._splitplot
        sp.reset(nrows, ncols)
        plotmeth = getattr(sp, 'add' + plotmeth)

        plot_objs=[copy.deepcopy(self._prepare_plotobj(*plot))
                   for plot in plots]

        try:
            sherpa.plot.begin()
            while plot_objs:
                plotmeth(plot_objs.pop(0))
#            while plots:
#                plotmeth(self._prepare_plotobj(*plots.pop(0)))

        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()

    def _plot(self, id, plotobj, *args,  **kwargs):
        #if len(args) > 0:
        #    raise SherpaError("cannot create plot for multiple data sets")
        obj = plotobj

        if not sherpa.utils.bool_cast(kwargs.pop('replot',False)):
            obj = self._prepare_plotobj(id, plotobj, *args)
        try:
            sherpa.plot.begin()
            obj.plot(**kwargs)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()
    
    def _overplot(self, id, plotobj, *args, **kwargs):
        #if len(args) > 0:
        #    raise SherpaError("cannot overplot for multiple data sets")
        obj = plotobj

        if not sherpa.utils.bool_cast(kwargs.pop('replot',False)):
            obj = self._prepare_plotobj(id, plotobj, *args)
        try:
            sherpa.plot.begin()
            obj.overplot(**kwargs)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    def _set_plot_item(self, plottype, item, value):
        keys = self._plot_types.keys()[:]

        if plottype.strip().lower() != "all":
            if plottype not in keys:
                raise sherpa.utils.err.PlotErr('wrongtype', plottype, str(keys))
            keys = [plottype]

        for key in keys:
            plot = self._plot_types[key]
            if _is_subclass(plot.__class__, sherpa.plot.Plot):
                plot.plot_prefs[item] = value
            elif _is_subclass(plot.__class__, sherpa.plot.Histogram):
                plot.histo_prefs[item] = value


    ### Ahelp ingest: 2015-04-28 DJB
    def set_xlog(self, plottype="all"):
        """New plots will display a logarithmically-scaled X axis.

        This setting only affects plots created after the call to
        `set_xlog`.

        Parameters
        ----------
        plottype : optional
           The type of plot that is to use a log-scaled X axis. The
           options are the same as accepted by `plot`, together with
           the 'all' option (which is the default setting).

        See Also
        --------
        plot : Create one or more plot types.
        set_xlinear : New plots will display a linear X axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Use a logarithmic scale for the X axis of `data` plots:

        >>> set_xlog('data')
        >>> plot('data', 'arf')

        All plots use a logarithmic scale for the X axis.

        >>> set_xlog()
        >>> plot_fit()

        """
        self._set_plot_item(plottype, 'xlog', True)


    ### Ahelp ingest: 2015-04-28 DJB
    def set_ylog(self, plottype="all"):
        """New plots will display a logarithmically-scaled Y axis.

        This setting only affects plots created after the call to
        `set_ylog`.

        Parameters
        ----------
        plottype : optional
           The type of plot that is to use a log-scaled X axis. The
           options are the same as accepted by `plot`, together with
           the 'all' option (which is the default setting).

        See Also
        --------
        plot : Create one or more plot types.
        set_xlog : New plots will display a logarithmically-scaled x axis.
        set_ylinear : New plots will display a linear Y axis.

        Examples
        --------

        Use a logarithmic scale for the Y axis of `data` plots:

        >>> set_ylog('data')
        >>> plot('data', 'arf')

        All plots use a logarithmic scale for the Y axis.

        >>> set_ylog()
        >>> plot_fit()

        """
        self._set_plot_item(plottype, 'ylog', True)


    ### Ahelp ingest: 2015-04-28 DJB
    def set_xlinear(self, plottype="all"):
        """New plots will display a linear X axis.

        This setting only affects plots created after the call to
        `set_xlinear`.

        Parameters
        ----------
        plottype : optional
           The type of plot that is to use a log-scaled X axis. The
           options are the same as accepted by `plot`, together with
           the 'all' option (which is the default setting).

        See Also
        --------
        plot : Create one or more plot types.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.

        Examples
        --------

        Use a linear X axis for 'data' plots:

        >>> set_xlinear('data')
        >>> plot('data', 'arf')

        All plots use a linear scale for the X axis.

        >>> set_xlinear()
        >>> plot_fit()

        """
        self._set_plot_item(plottype, 'xlog', False)


    ### Ahelp ingest: 2015-04-28 DJB
    def set_ylinear(self, plottype="all"):
        """New plots will display a linear Y axis.

        This setting only affects plots created after the call to
        `set_ylinear`.

        Parameters
        ----------
        plottype : optional
           The type of plot that is to use a log-scaled X axis. The
           options are the same as accepted by `plot`, together with
           the 'all' option (which is the default setting).

        See Also
        --------
        plot : Create one or more plot types.
        set_xlinear : New plots will display a linear X axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Use a linear Y axis for 'data' plots:

        >>> set_ylinear('data')
        >>> plot('data', 'arf')

        All plots use a linear scale for the Y axis.

        >>> set_ylinear()
        >>> plot_fit()

        """
        self._set_plot_item(plottype, 'ylog', False)


    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: how to describe optional plot types
    ### DOC-TODO: should we add plot_order
    ### DOC-TODO: how to list information/examples about the backends?
    ###           have some introductory text, but prob. need a link
    ###           to more information
    def plot(self, *args):
        """Create one or more plot types.

        The plot function creates one or more plots, depending on the
        arguments it is sent: a plot type, followed by an optional
        data set identifier, and this can be repeated. If no data set
        identifier is given for a plot type, the default identifier -
        as returned by `get_default_id` - is used.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           The data set does not support the requested plot type.

        See Also
        ---------
        get_default_id : Return the default data set identifier.
        sherpa.astro.ui.set_analysis : Set the units used when fitting and displaying spectral data.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Notes
        -----
        The supported plot types depend on the data set type, and
        include the following list. There are also individual
        functions, with `plot_` prepended to the plot type, such as
        `plot_data` (the `bkg` variants use a prefix of
        `plot_bkg_`). There are also several multiple-plot commands
        (e.g. `plot_fit_resid`).

        `arf`
           The ARF for the data set (only for `DataPHA` dataset).

        `bkg`
           The background.

        `bkgchisqr`
           The chi-squared statistic calculated for each bin when
           fitting the background.

        `bkgdelchi`
           The residuals for each bin, calculated as (data-model)
           divided by the error, for the background.

        `bkgfit`
           The data (as points) and the convolved model (as a line),
           for the background data set.

        `bkgmodel`
           The convolved background model.

        `bkgratio`
           The residuals for each bin, calculated as data/model,
           for the background data set.

        `bkgresid`
           The residuals for each bin, calculated as (data-model),
           for the background data set.

        `bkgsource`
           The un-convolved background model.

        `chisqr`
           The chi-squared statistic calculated for each bin.

        `data`
           The data (which may be background subtracted).

        `delchi`
           The residuals for each bin, calculated as (data-model)
           divided by the error.

        `fit`
           The data (as points) and the convolved model (as a line).
 
        `kernel`
           The PSF kernel associated with the data set.

        `model`
           The convolved model.

        `psf`
           The unfiltered PSF kernel associated with the data set.

        `ratio`
           The residuals for each bin, calculated as data/model.

        `resid`
           The residuals for each bin, calculated as (data-model).

        `source`
           The un-convolved model.

        The plots can be specialized for a particular data type,
        such as the `set_analysis` command controlling the units
        used for PHA data sets.

        See the documentation for the individual routines for
        information on how to configure the plots.

        The plot capabilities depend on what plotting backend, if any,
        is installed. If there is none available, a warning message
        will be displayed when `sherpa.ui` or `sherpa.astro.ui` is
        imported, and the `plot` set of commands will not create any
        plots. The choice of back end is made by changing the
        `options.plot_pkg` setting in the Sherpa configuration file.

        Examples
        --------

        Plot the data for the default data set. This is the same as
        `plot_data`.

        >>> plot("data")

        Plot the data for data set 2.

        >>> plot("data", 2)

        Plot the data and ARF for the default data set, in two
        seaparate plots.

        >>> plot("data", "arf")

        Plot the fit (data and model) for data sets 1 and 2, in two
        separate plots.

        >>> plot("fit", 1, "fit", 2)

        Plot the fit (data and model) for data sets "fit" and "jet",
        in two separate plots.

        >>> plot("fit", "nucleus", "fit", "jet")

        """
        self._multi_plot(args)        

    ### Ahelp ingest: 2015-04-29 DJB
    def plot_data(self, id=None, **kwargs):
        """Plot the data values.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_data`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        See Also
        --------
        get_data_plot : Return the data used by plot_data.
        get_data_plot_prefs : Return the preferences for plot_data.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        sherpa.astro.ui.set_analysis : Set the units used when fitting and displaying spectral data.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Plot the data from the default data set:

        >>> plot_data()

        Plot the data from data set 1:

        >>> plot_data(1)

        Plot the data from data set labelled "jet" and then overplot
        the "core" data set. The `set_xlog` command is used to select
        a logarithmic scale for the X axis.

        >>> set_xlog("data")
        >>> plot_data("jet")
        >>> plot_data("core", overplot=True)

        """
        self._plot(id, self._dataplot, **kwargs)
    
    # DOC-NOTE: also in sherpa.astro.utils
    ### Ahelp ingest: 2015-05-11 DJB
    def plot_model(self, id=None, **kwargs):
        """Plot the model for a data set.

        This function plots the model for a data set, which includes
        any instrument response (e.g. a convolution created by
        `set_psf`).

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_model`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        See Also
        --------
        get_model_plot : Return the data used by plot_model.
        get_model_plot_prefs : Return the preferences for plot_model.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_model_component : Plot a component of the model for a data set.
        plot_source : Plot the source expression for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Plot the convolved source model for the default data set:

        >>> plot_model()

        Overplot the model for data set 2 on data set 1:

        >>> plot_model(1)
        >>> plot_model(2, overplot=True)

        Create the equivalent of `plot_fit('jet')`:

        >>> plot_data('jet')
        >>> plot_model('jet', overplot=True)

        """
        self._plot(id, self._modelplot, **kwargs)

    # DOC-NOTE: also in sherpa.astro.utils, for now copies this text
    #           but does the astro version support a bkg_id parameter?
    ### Ahelp ingest: 2015-05-11 DJB
    def plot_source_component(self, id, model=None, **kwargs):
        """Plot a component of the source expression for a data set.

        This function evaluates and plots a component of the model
        expression for a data set, without any instrument response.
        Use `plot_model_component` to include any response.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to display (the name, if a string).
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_source_component`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        See Also
        --------
        get_source_component_plot : Return the data used by plot_source_component.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_model_component : Plot a component of the model for a data set.
        plot_source : Plot the source expression for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        Examples
        --------

        Overplot the `pl` component of the source expression for
        the default data set:

        >>> plot_source()
        >>> plot_source_component(pl, overplot=True)

        """

        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        plotobj = self._compsrcplot
        if isinstance(model, sherpa.models.TemplateModel):
            plotobj = self._comptmplsrcplot

        self._plot(id, plotobj, model, **kwargs)


    # DOC-NOTE: also in sherpa.astro.utils, for now copies this text
    #           but does the astro version support a bkg_id parameter?
    ### Ahelp ingest: 2015-05-11 DJB
    def plot_model_component(self, id, model=None, **kwargs):
        """Plot a component of the model for a data set.

        This function evaluates and plots a component of the model
        expression for a data set, including any instrument response.
        Use `plot_source_component` to display without any response.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to display (the name, if a string).
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_model_component`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        See Also
        --------
        get_model_component_plot : Return the data used by plot_model_component.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_source_component : Plot a component of the source expression for a data set.
        plot_model : Plot the model for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        Examples
        --------

        Overplot the `pl` component of the model expression for
        the default data set:

        >>> plot_model()
        >>> plot_model_component(pl, overplot=True)

        Display the results for the 'jet' data set (data and model),
        and then overplot the `pl` component evaluated for the 'jet'
        and 'core' data sets:

        >>> plot_fit('jet')
        >>> plot_model_component('jet', pl, overplot=True)
        >>> plot_model_component('core', pl, overplot=True)

        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        is_source = self._get_model_status(id)[1]
        model = self._add_convolution_models(id, self.get_data(id),
                                             model, is_source)

        self._plot(id, self._compmdlplot, model, **kwargs)


    # DOC-NOTE: also in sherpa.astro.utils, but with extra lo/hi arguments
    ### Ahelp ingest: 2015-05-11 DJB
    def plot_source(self, id=None, **kwargs):
        """Plot the source expression for a data set.

        This function plots the source model for a data set. This does
        not include any instrument response (e.g. a convolution
        created by `set_psf`).

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_source`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        See Also
        --------
        get_source_plot : Return the data used by plot_source.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_model : Plot the model for a data set.
        plot_source_component : Plot a component of the source expression for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Plot the unconvolved source model for the default data set:

        >>> plot_source()

        Overplot the source model for data set 2 on data set 1:

        >>> plot_source(1)
        >>> plot_source(2, overplot=True)

        """
        id = self._fix_id(id)
        mdl = self._models.get(id, None)
        if mdl is not None:
            raise IdentifierErr("Convolved model\n'%s'\n is set for dataset %s. You should use plot_model instead." % 
                    (mdl.name, str(id)))
        self._plot(id, self._sourceplot, **kwargs)

    ### Ahelp ingest: 2015-05-11 DJB
    def plot_fit(self, id=None, **kwargs):
        """Plot the fit results (data, model) for a data set.

        This function creates a plot containing the data and the model
        (including any instrument response) for a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_fit`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_fit_plot : Return the data used by plot_fit.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_fit_delchi : Plot the fit results, and the residuals, for a data set.
        plot_fit_resid : Plot the fit results, and the residuals, for a data set.
        plot_data : Plot the data values.
        plot_model : Plot the model for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Plot the fit results for the default data set:

        >>> plot_fit()

        Overplot the 'core' results on those from the 'jet' data set,
        using a logarithmic scale for the X axis:

        >>> set_xlog()
        >>> plot_fit('jet')
        >>> plot_fit('core', overplot=True)

        """
        self._plot(id, self._fitplot, **kwargs)

    ### Ahelp ingest: 2015-05-11 DJB
    def plot_resid(self, id=None, **kwargs):
        """Plot the residuals (data - model) for a data set.

        This function displays the residuals (data - model) for a data
        set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_resid`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_resid_plot : Return the data used by plot_resid.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_chisqr : Plot the chi-squared value for each point in a data set.
        plot_delchi : Plot the ratio of residuals to error for a data set.
        plot_ratio : Plot the ratio of data to model for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Plot the residuals for the default data set:

        >>> plot_resid()

        Overplot the residuals from the 'core' data set on those
        from the 'jet' dataset:

        >>> plot_resid('jet')
        >>> plot_resid('core', overplot=True)

        Add the residuals to the plot of the data, for the default
        data set:

        >>> plot_data()
        >>> plot_resid(overplot=True)

        """
        self._plot(id, self._residplot, **kwargs)

    ### Ahelp ingest: 2015-05-11 DJB
    def plot_chisqr(self, id=None, **kwargs):
        """Plot the chi-squared value for each point in a data set.

        This function displays the square of the residuals (data -
        model) divided by the error, for a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_chisqr`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_chisqr_plot : Return the data used by plot_chisqr.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_delchi : Plot the ratio of residuals to error for a data set.
        plot_ratio : Plot the ratio of data to model for a data set.
        plot_resid : Plot the residuals (data - model) for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Plot the chi-quare values for each point in the default data
        set:

        >>> plot_chisqr()

        Overplot the values from the 'core' data set on those
        from the 'jet' dataset:

        >>> plot_chisqr('jet')
        >>> plot_chisqr('core', overplot=True)

        """
        self._plot(id, self._chisqrplot, **kwargs)

    ### Ahelp ingest: 2015-05-11 DJB
    def plot_delchi(self, id=None, **kwargs):
        """Plot the ratio of residuals to error for a data set.

        This function displays the residuals (data - model) divided by
        the error, for a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_delchi`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_delchi_plot : Return the data used by plot_delchi.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_chisqr : Plot the chi-squared value for each point in a data set.
        plot_ratio : Plot the ratio of data to model for a data set.
        plot_resid : Plot the residuals (data - model) for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Plot the residuals for the default data set, divided by
        the error value for each bin:

        >>> plot_delchi()

        Overplot the values from the 'core' data set on those
        from the 'jet' dataset:

        >>> plot_delchi('jet')
        >>> plot_delchi('core', overplot=True)

        """
        self._plot(id, self._delchiplot, **kwargs)
        
    ### Ahelp ingest: 2015-05-11 DJB
    def plot_ratio(self, id=None, **kwargs):
        """Plot the ratio of data to model for a data set.

        This function displays the ratio data / model for a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_ratio`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_ratio_plot : Return the data used by plot_ratio.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_chisqr : Plot the chi-squared value for each point in a data set.
        plot_delchi : Plot the ratio of residuals to error for a data set.
        plot_resid : Plot the residuals (data - model) for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Plot the ratio of data to model for the default data set:

        >>> plot_ratio()

        Overplot the ratios from the 'core' data set on those from the
        'jet' dataset:

        >>> plot_ratio('jet')
        >>> plot_ratio('core', overplot=True)

        """
        self._plot(id, self._ratioplot, **kwargs)

    ### Ahelp ingest: 2015-05-11 DJB
    def plot_psf(self, id=None, **kwargs):
        """Plot the 1D PSF model applied to a data set.

        The `plot_kernel` function shows the data used to convolve
        the model.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_psf`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If a PSF model has not been created for the data set.

        See Also
        --------
        get_psf_plot : Return the data used by plot_psf.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_kernel : Plot the 1D kernel applied to a data set.
        set_psf : Add a PSF model to a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Create a model (a step function) that is convolved by
        a gaussian, and display the PSF:

        >>> dataspace1d(1, 10, step=1, dstype=Data1D)
        >>> set_model(steplo1d.stp)
        >>> stp.xcut = 4.4
        >>> load_psf('psf1', gauss1d.gline)
        >>> set_psf('psf1')
        >>> gline.fwhm = 1.2
        >>> plot_psf()

        """
        self._plot(id, self._psfplot, **kwargs)


    ### Ahelp ingest: 2015-05-11 DJB
    def plot_kernel(self, id=None, **kwargs):
        """Plot the 1D kernel applied to a data set.

        The `plot_psf` function shows the full PSF, from which the
        kernel is derived.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_kernel`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If a PSF model has not been created for the data set.

        See Also
        --------
        get_kernel_plot : Return the data used by plot_kernel.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_psf : Plot the 1D PSF model applied to a data set.
        set_psf : Add a PSF model to a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Create a model (a step function) that is convolved by
        a gaussian, and display the kernel overplotted on the
        PSF:

        >>> dataspace1d(1, 10, step=1, dstype=Data1D)
        >>> set_model(steplo1d.stp)
        >>> stp.xcut = 4.4
        >>> load_psf('psf1', gauss1d.gline)
        >>> set_psf('psf1')
        >>> gline.fwhm = 1.2
        >>> plot_psf()
        >>> plot_kernel(overplot=True)

        """
        self._plot(id, self._kernelplot, **kwargs)


    ### Ahelp ingest: 2015-05-11 DJB
    def plot_fit_resid(self, id=None, replot=False, overplot=False,
                       clearwindow=True):
        """Plot the fit results, and the residuals, for a data set.

        This creates two plots - the first from `plot_fit` and the
        second from `plot_resid` - for a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_fit_resid`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.
        clearwindow : bool, optional
           When using ChIPS for plotting, should the existing frame
           be cleared before creating the plot?

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_fit_plot : Return the data used by plot_fit.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_fit : Plot the fit results for a data set.
        plot_fit_delchi : Plot the fit results, and the residuals, for a data set.
        plot_data : Plot the data values.
        plot_model : Plot the model for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Plot the results for the default data set:

        >>> plot_fit_resid()

        Overplot the 'core' results on those from the 'jet' data set,
        using a logarithmic scale for the X axis:

        >>> set_xlog()
        >>> plot_fit_resid('jet')
        >>> plot_fit_resid('core', overplot=True)

        """
        self._jointplot.reset()
        fp = self._fitplot
        rp = self._residplot
        if not sherpa.utils.bool_cast(replot):
            fp = self._prepare_plotobj(id, fp)
            rp = self._prepare_plotobj(id, rp)
        try:
            sherpa.plot.begin()            
            self._jointplot.plottop(fp, overplot=overplot,
                                    clearwindow=clearwindow)

            oldval = rp.plot_prefs['xlog']
            if ((self._dataplot.plot_prefs.has_key('xlog') and
                 self._dataplot.plot_prefs['xlog']) or 
                (self._modelplot.plot_prefs.has_key('xlog') and
                 self._modelplot.plot_prefs['xlog'])):
                rp.plot_prefs['xlog']=True

            self._jointplot.plotbot(rp, overplot=overplot)

            rp.plot_prefs['xlog'] = oldval
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()

    ### Ahelp ingest: 2015-05-11 DJB
    def plot_fit_delchi(self, id=None, replot=False, overplot=False,
                        clearwindow=True):
        """Plot the fit results, and the residuals, for a data set.

        This creates two plots - the first from `plot_fit` and the
        second from `plot_delchi` - for a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_fit_delchi`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.
        clearwindow : bool, optional
           When using ChIPS for plotting, should the existing frame
           be cleared before creating the plot?

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_fit_plot : Return the data used by plot_fit.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_fit : Plot the fit results for a data set.
        plot_fit_resid : Plot the fit results, and the residuals, for a data set.
        plot_data : Plot the data values.
        plot_model : Plot the model for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Examples
        --------

        Plot the results for the default data set:

        >>> plot_fit_delchi()

        Overplot the 'core' results on those from the 'jet' data set,
        using a logarithmic scale for the X axis:

        >>> set_xlog()
        >>> plot_fit_delchi('jet')
        >>> plot_fit_delchi('core', overplot=True)

        """
        self._jointplot.reset()
        fp = self._fitplot
        dp = self._delchiplot
        if not sherpa.utils.bool_cast(replot):
            fp = self._prepare_plotobj(id, fp)
            dp = self._prepare_plotobj(id, dp)
        try:
            sherpa.plot.begin()           
            self._jointplot.plottop(fp, overplot=overplot,
                                    clearwindow=clearwindow)

            oldval = dp.plot_prefs['xlog']
            if ((self._dataplot.plot_prefs.has_key('xlog') and
                 self._dataplot.plot_prefs['xlog']) or 
                (self._modelplot.plot_prefs.has_key('xlog') and
                 self._modelplot.plot_prefs['xlog'])):
                dp.plot_prefs['xlog']=True

            self._jointplot.plotbot(dp, overplot=overplot)

            dp.plot_prefs['xlog'] = oldval
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    #
    ## Statistical plotting routines
    #

    ### Ahelp ingest: 2015-04-30 DJB
    def plot_pdf(self, points, name="x", xlabel="x", bins=12, normed=True, 
                 replot=False, overplot=False, clearwindow=True ):
        """Plot the probability density function of an array of values.

        Create and plot the probability density function (PDF) of
        the input array.

        Parameters
        ----------
        points : array
           The values used to create the probability density function.
        name : str, optional
           The label to use as part of the plot title.
        xlabel : str, optional
           The label for the X axis
        bins : int, optional
           The number of bins to use to create the PDF.
        normed : bool, optional
           Should the PDF be normalized (the default is `True`).
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_pdf`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.
        clearwindow : bool, optional
           When using ChIPS for plotting, should the existing frame
           be cleared before creating the plot?

        See Also
        --------
        get_draws : Run the pyBLoCXS MCMC algorithm.
        get_pdf_plot : Return the data used to plot the last PDF.
        plot_cdf : Plot the cumulative density function of an array.
        plot_scatter : Create a scatter plot.

        Examples
        --------

        >>> mu, sigma, n = 100, 15, 500
        >>> x = np.random.normal(loc=mu, scale=sigma, size=n)
        >>> plot_pdf(x, bins=25)

        >>> plot_pdf(x, normed=False, xlabel="mu", name="Simulations")

        """
        if not sherpa.utils.bool_cast(replot):
            self._pdfplot.prepare(points, bins, normed, xlabel, name)

        try:
            sherpa.plot.begin()
            self._pdfplot.plot(overplot, clearwindow)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    def get_pdf_plot(self):
        """Return the data used to plot the last PDF.

        Returns
        -------
        plot : sherpa.plot.PDFPlot instance
           An object containing the data used by the last call to
           `plot_pdf`. The fields will be `None` if the function
           has not been called.

        See Also
        --------
        plot_pdf : Plot the probability density function of an array.

        """
        return self._pdfplot


    ### Ahelp ingest: 2015-04-30 DJB
    def plot_cdf(self, points, name="x", xlabel="x", 
                 replot=False, overplot=False, clearwindow=True ):
        """Plot the cumulative density function of an array of values.

        Create and plot the cumulative density function (CDF) of
        the input array. Median and upper- and lower- quartiles
        are marked on the plot.

        Parameters
        ----------
        points : array
           The values used to create the cumulative density function.
        name : str, optional
           The label to use as part of the plot title.
        xlabel : str, optional
           The label for the X and part of the Y axes.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_pdf`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.
        clearwindow : bool, optional
           When using ChIPS for plotting, should the existing frame
           be cleared before creating the plot?

        See Also
        --------
        get_cdf_plot : Return the data used to plot the last CDF.
        get_draws : Run the pyBLoCXS MCMC algorithm.
        plot_pdf : Plot the probability density function of an array.
        plot_scatter : Create a scatter plot.

        Examples
        --------

        >>> mu, sigma, n = 100, 15, 500
        >>> x = np.random.normal(loc=mu, scale=sigma, size=n)
        >>> plot_cdf(x)

        >>> plot_cdf(x, xlabel="x pos", name="Simulations")

        """

        if not sherpa.utils.bool_cast(replot):
            self._cdfplot.prepare(points, xlabel, name)

        try:
            sherpa.plot.begin()
            self._cdfplot.plot(overplot, clearwindow)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    def get_cdf_plot(self):
        """Return the data used to plot the last CDF.

        Returns
        -------
        plot : sherpa.plot.CDFPlot instance
           An object containing the data used by the last call to
           `plot_cdf`. The fields will be `None` if the function
           has not been called.

        See Also
        --------
        plot_cdf : Plot the cumulative density function of an array.

        """
        return self._cdfplot


    ### Ahelp ingest: 2015-04-30 DJB
    ### DOC-TODO: what does xlabel do?
    ### DOC-TODO: is clearwindow a ChIPS-only setting?
    def plot_trace(self, points, name="x", xlabel="x", 
                   replot=False, overplot=False, clearwindow=True ):
        """Create a trace plot of row number versus value.

        Dispay a plot of the `points` array values (Y axis) versus row
        number (X axis). This can be useful to view how a value
        changes, such as the value of a parameter returned by
        `get_draws`.

        Parameters
        ----------
        points : array
           The values to plot on the Y axis.
        name : str, optional
           The label to use on the Y axis and as part of the plot
           title.
        xlabel : str, optional
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_trace`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.
        clearwindow : bool, optional
           When using ChIPS for plotting, should the existing frame
           be cleared before creating the plot?

        See Also
        --------
        get_draws : Run the pyBLoCXS MCMC algorithm.
        get_trace_plot : Return the data used to plot the last trace.
        plot_cdf : Plot the cumulative density function of an array.
        plot_pdf : Plot the probability density function of an array.
        plot_scatter : Create a scatter plot.

        Examples
        --------

        Plot the trace of the 500 elements in the `x` array:

        >>> mu, sigma = 100, 15
        >>> x = mu + sigma * np.random.randn(500)
        >>> plot_trace(x)

        Use "ampl" as the Y axis label:

        >>> plot_trace(ampl, name='ampl')

        """

        if not sherpa.utils.bool_cast(replot):
            self._traceplot.prepare(points, xlabel, name)

        try:
            sherpa.plot.begin()
            self._traceplot.plot(overplot, clearwindow)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    def get_trace_plot(self):
        """Return the data used to plot the last trace.

        Returns
        -------
        plot : sherpa.plot.TracePlot instance
           An object containing the data used by the last call to
           `plot_trace`. The fields will be `None` if the function
           has not been called.

        See Also
        --------
        plot_trace : Create a trace plot of row number versus value.

        """
        return self._traceplot


    ### Ahelp ingest: 2015-04-30 DJB
    def plot_scatter(self, x, y, name="(x,y)", xlabel="x", ylabel="y",
                   replot=False, overplot=False, clearwindow=True ):
        """Create a scatter plot.

        Parameters
        ----------
        x : array
           The values to plot on the X axis.
        y : array
           The values to plot on the Y axis. This must match the size
           of the `x` array.
        name : str, optional
           The plot title.
        xlabel : str, optional
           The label for the X axis.
        ylabel : str, optional
           The label for the Y axis.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_scatter`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.
        clearwindow : bool, optional
           When using ChIPS for plotting, should the existing frame
           be cleared before creating the plot?

        See Also
        --------
        get_scatter_plot : Return the data used to plot the last scatter plot.
        plot_trace : Create a trace plot of row number versus value.

        Examples
        --------

        Plot the X and Y points:

        >>> mu, sigma, n = 100, 15, 500
        >>> x = mu + sigma * np.random.randn(n)
        >>> y = mu + sigma * np.random.randn(n)
        >>> plot_scatter(x, y)

        Change the axis labels and the plot title:

        >>> plot_scatter(nh, kt, xlabel='nH', ylabel='kT', name='Simulations')

        """
        if not sherpa.utils.bool_cast(replot):
            self._scatterplot.prepare(x, y, xlabel, ylabel, name)

        try:
            sherpa.plot.begin()
            self._scatterplot.plot(overplot, clearwindow)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    def get_scatter_plot(self):
        """Return the data used to plot the last scatter plot.

        Returns
        -------
        plot : sherpa.plot.ScatterPlot instance
           An object containing the data used by the last call to
           `plot_scatter`. The fields will be `None` if the function
           has not been called.

        See Also
        --------
        plot_scatter : Create a scatter plot.

        """
        return self._scatterplot



    #
    # Contours
    #

    def _contour(self, id, plotobj, **kwargs):
        #if len(args) > 0:
        #raise SherpaError("cannot create contour plot for multiple data sets")
        obj = plotobj

        if not sherpa.utils.bool_cast(kwargs.pop('replot',False)):
            obj = self._prepare_plotobj(id, plotobj)
        try:
            sherpa.plot.begin()
            obj.contour(**kwargs)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()

    def _overcontour(self, id, plotobj, **kwargs):
        #if len(args) > 0:
        #raise SherpaError("cannot overplot contours for multiple data sets")
        obj = plotobj

        if not sherpa.utils.bool_cast(kwargs.pop('replot',False)):
            obj = self._prepare_plotobj(id, plotobj)
        try:
            sherpa.plot.begin()
            obj.overcontour(**kwargs)
        except:
            sherpa.plot.exceptions()
            raise
        else:            
            sherpa.plot.end()

    ### Ahelp ingest: 2015-05-07 DJB
    ### DOC-TODO: how to describe optional plot types
    ### DOC-TODO: how to list information/examples about the backends?
    ###           have some introductory text, but prob. need a link
    ###           to more information
    def contour(self, *args):
        """Create a contour plot for an image data set.

        Create one or more contour plots, depending on the arguments
        it is set: a plot type, followed by an optional data set
        identifier, and this can be repeated. If no data set
        identifier is given for a plot type, the default identifier -
        as returned by `get_default_id` - is used. This is for
        2D data sets.

        Raises
        ------
        sherpa.utils.err.DataErr
           The data set does not support the requested plot type.

        See Also
        ---------
        contour_data : Contour the values of an image data set.
        contour_fit : Contour the fit to a data set.
        contour_fit_resid :
        contour_kernel : Contour the kernel applied to the model of an image data set.
        contour_model : Contour the values of the model, including any PSF.
        contour_psf : Contour the PSF applied to the model of an image data set.
        contour_ratio : Contour the ratio of data to model.
        contour_resid : Contour the residuals of the fit.
        contour_source : Contour the values of the model, without any PSF.
        get_default_id : Return the default data set identifier.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.

        Notes
        -----
        The supported plot types depend on the data set type, and
        include the following list. There are also individual
        functions, with `contour_` prepended to the plot type, such as
        `contour_data` and the `contour_fit_resid` variant.

        `data`
           The data.

        `fit`
           Contours of the data and the source model.

        `fit_resid`
           Two plots: the first is the contours of the data and the
           source model and the second is the residuals.

        `kernel`
           The kernel.

        `model`
           The source model including any PSF convolution set by
           `set_psf`.

        `psf`
           The PSF.

        `ratio`
           Contours of the ratio image, formed by dividing the data by
           the model.

        `resid`
           Contours of the residual image, formed by subtracting the
           model from the data.

        `source`
           The source model (without any PSF convolution set by
           `set_psf`).

        Examples
        --------

        >>> contour('data')

        >>> contour('data', 1, 'data', 2)

        >>> contour('data', 'model')

        """
        self._multi_plot(args, 'contour')

    ### Ahelp ingest: 2015-05-07 DJB
    def contour_data(self, id=None, **kwargs):
        """Contour the values of an image data set.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `contour_data`. The default is `False`.
        overcontour : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new contour plot. The default is `False`.

        See Also
        --------
        get_data_contour : Return the data used by contour_data.
        get_data_contour_prefs : Return the preferences for contour_data.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.

        Examples
        --------

        Plot the data from the default data set:

        >>> contour_data()

        Contour the data and then overplot the data from the second
        data set:

        >>> contour_data()
        >>> contour_data(2, overcontour=True)

        """
        self._contour(id, self._datacontour, **kwargs)
        
    ### Ahelp ingest: 2015-05-07 DJB
    def contour_model(self, id=None, **kwargs):
        """Contour the values of the model, including any PSF.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the model. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `contour_model`. The default is `False`.
        overcontour : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new contour plot. The default is `False`.

        See Also
        --------
        get_model_contour : Return the data used by contour_model.
        get_model_contour_prefs : Return the preferences for contour_model.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.
        set_psf : Add a PSF model to a data set.

        Examples
        --------

        Plot the model from the default data set:

        >>> contour_model()

        Compare the model without and with the PSF component,
        for the "img" data set:

        >>> contour_source("img")
        >>> contour_model("img", overcontour=True)

        """
        self._contour(id, self._modelcontour, **kwargs)

    ### Ahelp ingest: 2015-05-07 DJB
    def contour_source(self, id=None, **kwargs):
        """Contour the values of the model, without any PSF.

        The preferences are the same as `contour_model`.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the model. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `contour_source`. The default is `False`.
        overcontour : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new contour plot. The default is `False`.

        See Also
        --------
        get_source_contour : Return the data used by contour_source.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.
        set_psf : Add a PSF model to a data set.

        Examples
        --------

        Plot the model from the default data set:

        >>> contour_source()

        Compare the model without and with the PSF component,
        for the "img" data set:

        >>> contour_model("img")
        >>> contour_source("img", overcontour=True)

        """
        self._contour(id, self._sourcecontour, **kwargs)

    ### Ahelp ingest: 2015-05-08 DJB
    def contour_fit(self, id=None, **kwargs):
        """Contour the fit to a data set.

        Overplot the model - including any PSF - on the data. The
        preferences are the same as `contour_data` and
        `contour_model`.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data and model. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `contour_fit`. The default is `False`.
        overcontour : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new contour plot. The default is `False`.

        See Also
        --------
        get_fit_contour : Return the data used by contour_fit.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.

        Examples
        --------

        Plot the fit for the default data set:

        >>> contour_fit()

        Overplot the fit to data set 's2' on that of the default data
        set:

        >>> contour_fit()
        >>> contour_fit('s2', overcontour=True)

        """
        self._contour(id, self._fitcontour, **kwargs)

    ### Ahelp ingest: 2015-05-08 DJB
    def contour_resid(self, id=None, **kwargs):
        """Contour the residuals of the fit.

        The residuals are formed by subtracting the current model -
        including any PSF - from the data.  The preferences are the
        same as `contour_data`.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data and model. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `contour_resid`. The default is `False`.
        overcontour : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new contour plot. The default is `False`.

        See Also
        --------
        get_resid_contour : Return the data used by contour_resid.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.

        Examples
        --------

        Plot the residuals from the default data set:

        >>> contour_resid()

        Overplot the residuals on the model:

        >>> contour_model('img')
        >>> contour_resid('img', overcontour=True)

        """
        self._contour(id, self._residcontour, **kwargs)
    
    ### Ahelp ingest: 2015-05-08 DJB
    def contour_ratio(self, id=None, **kwargs):
        """Contour the ratio of data to model.

        The ratio image is formed by dividing the data by the current
        model, including any PSF. The preferences are the same as
        `contour_data`.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data and model. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `contour_ratio`. The default is `False`.
        overcontour : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new contour plot. The default is `False`.

        See Also
        --------
        get_ratio_contour : Return the data used by contour_ratio.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.

        Examples
        --------

        Plot the ratio from the default data set:

        >>> contour_ratio()

        Overplot the ratio on the residuals:

        >>> contour_resid('img')
        >>> contour_ratio('img', overcontour=True)

        """
        self._contour(id, self._ratiocontour, **kwargs)

    ### Ahelp ingest: 2015-05-08 DJB
    def contour_psf(self, id=None, **kwargs):
        """Contour the PSF applied to the model of an image data set.

        If the data set has no PSF applied to it, the model is
        displayed.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the model. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `contour_psf`. The default is `False`.
        overcontour : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new contour plot. The default is `False`.

        See Also
        --------
        get_psf_contour : Return the data used by contour_psf.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        contour_kernel : Contour the kernel applied to the model of an image data set.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.
        set_psf : Add a PSF model to a data set.

        """
        self._contour(id, self._psfcontour, **kwargs)


    ### Ahelp ingest: 2015-05-08 DJB
    def contour_kernel(self, id=None, **kwargs):
        """Contour the kernel applied to the model of an image data set.

        If the data set has no PSF applied to it, the model is
        displayed.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the model. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `contour_kernel`. The default is `False`.
        overcontour : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new contour plot. The default is `False`.

        See Also
        --------
        get_psf_contour : Return the data used by contour_psf.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        contour_psf : Contour the PSF applied to the model of an image data set.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.
        set_psf : Add a PSF model to a data set.

        """
        self._contour(id, self._kernelcontour, **kwargs)

    ### Ahelp ingest: 2015-05-08 DJB
    def contour_fit_resid(self, id=None, replot=False, overcontour=False):
        """Contour the fit and the residuals to a data set.

        Overplot the model - including any PSF - on the data. In a
        separate plot contour the residuals. The preferences are the
        same as `contour_data` and `contour_model`.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data and model. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `contour_fit_resid`. The default is `False`.
        overcontour : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new contour plot. The default is `False`.

        See Also
        --------
        get_fit_contour : Return the data used by contour_fit.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        contour_fit : Contour the fit to a data set.
        contour_resid : Contour the residuals of the fit.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.

        Examples
        --------

        Plot the fit and residuals for the default data set:

        >>> contour_fit_resid()

        """
        self._splitplot.reset()
        fc = self._fitcontour
        rc = self._residcontour
        if not sherpa.utils.bool_cast(replot):
            fc = self._prepare_plotobj(id, fc)
            rc = self._prepare_plotobj(id, rc)
        try:
            sherpa.plot.begin()
            self._splitplot.addcontour(fc, overcontour=overcontour)
            self._splitplot.addcontour(rc, overcontour=overcontour)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    ###########################################################################
    # Projection and uncertainty plots
    ###########################################################################

    ### Ahelp ingest: 2015-05-07 DJB
    ### DOC-NOTE: I am not convinced that this code is working when recalc=True
    ### DOC-NOTE: needs to support the fast option of int_proj
    def get_int_proj(self, par=None, id=None, otherids=None, recalc=False,
                     min=None, max=None, nloop=20, delv=None, fac=1, 
                     log=False, numcores=None):
        """Return the interval-projection object.

        This returns (and optionally calculates) the data used to
        display the `int_proj` plot.

        Parameters
        ----------
        par
           The parameter to plot.
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        recalc : bool, optional
           The default value (`False`) means that the results from the
           last call to `int_proj` (or `get_int_proj`) are returned,
           ignoring the other parameter values. Otherwise, the
           statistic curve is re-calculated, but not plotted.
        min : number, optional
           The minimum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : number, optional
           The maximum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : int, optional
           The number of steps to use. This is used when `delv` is set
           to `None`.
        delv : number, optional
           The step size for the parameter. Setting this over-rides
           the `nloop` parameter. The default is `None`.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : bool, optional
           Should the step size be logarithmically spaced? The
           default (`False`) is to use a linear grid.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        iproj : sherpa.plot.IntervalProjection instance
           The fields of this object can be used to re-create the plot
           created by `int_proj`.

        See Also
        --------
        conf : Estimate the confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        int_proj : Calculate and plot the fit statistic versus fit parameter value.
        int_unc : Calculate and plot the fit statistic versus fit parameter value.
        reg_proj : Plot the statistic value as two parameters are varied.

        Examples
        --------

        Return the results of the `int_proj` run:

        >>> int_proj(src.xpos)
        >>> iproj = get_int_proj()
        >>> min(iproj.y)
        119.55942437129544

        Create the data without creating a plot:

        >>> iproj = get_int_proj(pl.gamma, recalc=True)

        Control how the data is created

        >>> iproj = get_int_proj(pl.gamma, id="src", min=12, max=14,
                                 nloop=51, recalc=True)

        """
        if sherpa.utils.bool_cast(recalc):
            par = self._check_par(par)
            if otherids is None:
                otherids = ()
            ids, fit = self._get_fit(id, otherids)
            self._intproj.prepare(min, max, nloop, delv, fac,
                                  sherpa.utils.bool_cast(log), numcores)
            self._intproj.calc(fit,par,self._methods)
        return self._intproj

    ### Ahelp ingest: 2015-05-07 DJB
    ### DOC-NOTE: Check that this works (since get_int_proj may not) when recalc=True
    def get_int_unc(self, par=None, id=None, otherids=None, recalc=False,
                    min=None, max=None, nloop=20, delv=None, fac=1, log=False,
                    numcores=None):
        """Return the interval-uncertainty object.

        This returns (and optionally calculates) the data used to
        display the `int_unc` plot.

        Parameters
        ----------
        par
           The parameter to plot.
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        recalc : bool, optional
           The default value (`False`) means that the results from the
           last call to `int_proj` (or `get_int_proj`) are returned,
           ignoring the other parameter values. Otherwise, the
           statistic curve is re-calculated, but not plotted.
        min : number, optional
           The minimum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : number, optional
           The maximum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : int, optional
           The number of steps to use. This is used when `delv` is set
           to `None`.
        delv : number, optional
           The step size for the parameter. Setting this over-rides
           the `nloop` parameter. The default is `None`.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : bool, optional
           Should the step size be logarithmically spaced? The
           default (`False`) is to use a linear grid.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        iunc : sherpa.plot.IntervalUncertainty instance
           The fields of this object can be used to re-create the plot
           created by `int_unc`.

        See Also
        --------
        conf : Estimate the confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        int_proj : Calculate and plot the fit statistic versus fit parameter value.
        int_unc : Calculate and plot the fit statistic versus fit parameter value.
        reg_proj : Plot the statistic value as two parameters are varied.

        Examples
        --------

        Return the results of the `int_unc` run:

        >>> int_unc(src.xpos)
        >>> iunc = get_int_unc()
        >>> min(iunc.y)
        119.55942437129544

        Create the data without creating a plot:

        >>> iunc = get_int_unc(pl.gamma, recalc=True)

        Control how the data is created

        >>> iunc = get_int_unc(pl.gamma, id="src", min=12, max=14,
                               nloop=51, recalc=True)

        """
        if sherpa.utils.bool_cast(recalc):
            par = self._check_par(par)
            if otherids is None:
                otherids = ()
            ids, fit = self._get_fit(id, otherids)
            self._intunc.prepare(min, max, nloop, delv, fac,
                                 sherpa.utils.bool_cast(log), numcores)
            self._intunc.calc(fit,par)
        return self._intunc

    ### Ahelp ingest: 2015-05-08 DJB
    def get_reg_proj(self, par0=None, par1=None, id=None, otherids=None,
                     recalc=False, fast=True, min=None, max=None, 
                     nloop=(10,10),delv=None, fac=4, log=(False,False),
                     sigma=(1,2,3), levels=None, numcores=None):
        """Return the region-projection object.

        This returns (and optionally calculates) the data used to
        display the `reg_proj` contour plot.

        Parameters
        ----------
        par0, par1
           The parameters to plot on the X and Y axes, respectively.
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        recalc : bool, optional
           The default value (`False`) means that the results from the
           last call to `reg_proj` (or `get_reg_proj`) are returned,
           ignoring the other parameter values. Otherwise, the
           statistic curve is re-calculated, but not plotted.
        fast : bool, optional
           If `True` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is `False`.
        min : pair of numbers, optional
           The minimum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : pair of number, optional
           The maximum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : pair of int, optional
           The number of steps to use. This is used when `delv` is set
           to `None`.
        delv : pair of number, optional
           The step size for the parameter. Setting this over-rides
           the `nloop` parameter. The default is `None`.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : pair of bool, optional
           Should the step size be logarithmically spaced? The
           default (`False`) is to use a linear grid.
        sigma : sequence of number, optional
           The levels at which to draw the contours. The units are the
           change in significance relative to the starting value,
           in units of sigma.
        levels : sequence of number, optional
           The numeric values at which to draw the contours. This
           over-rides the `sigma` parameter, if set (the default is
           `None`).
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        rproj : sherpa.plot.RegionProjection instance
           The fields of this object can be used to re-create the plot
           created by `reg_proj`.

        See Also
        --------
        conf : Estimate the confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        int_proj : Calculate and plot the fit statistic versus fit parameter value.
        int_unc : Calculate and plot the fit statistic versus fit parameter value.
        reg_proj : Plot the statistic value as two parameters are varied.
        reg_unc : Plot the statistic value as two parameters are varied.

        Examples
        --------

        Return the results for the `reg_proj` run:

        >>> reg_proj(src.xpos, src.ypos)
        >>> rproj = get_reg_proj()

        Create the data without creating a plot:

        >>> rproj = get_reg_proj(pl.gamma, gal.nh, recalc=True)

        Control how the data is created:

        >>> rproj = get_reg_proj(pl.gamma, gal.nh, id="src",
                                 min=(0.5,0.01), max=(2.5,1),
                                 nloop=(51,51), log=(False,True),
                                 recalc=True)

        """
        if sherpa.utils.bool_cast(recalc):
            par0 = self._check_par(par0, 'par0')
            par1 = self._check_par(par1, 'par1')
            if otherids is None:
                otherids = ()
            ids, fit = self._get_fit(id, otherids)
            self._regproj.prepare(fast, min, max, nloop, delv, fac,
                                  sherpa.utils.bool_cast(log), sigma, levels,
                                  numcores)
            self._regproj.calc(fit,par0,par1,self._methods)
        return self._regproj
    
    ### Ahelp ingest: 2015-05-08 DJB
    def get_reg_unc(self, par0=None, par1=None, id=None, otherids=None,
                    recalc=False, min=None, max=None, nloop=(10,10), delv=None,
                    fac=4, log=(False,False), sigma=(1,2,3), levels=None,
                    numcores=None):
        """Return the region-uncertainty object.

        This returns (and optionally calculates) the data used to
        display the `reg_unc` contour plot.

        Parameters
        ----------
        par0, par1
           The parameters to plot on the X and Y axes, respectively.
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        recalc : bool, optional
           The default value (`False`) means that the results from the
           last call to `reg_unc` (or `get_reg_unc`) are returned,
           ignoring the other parameter values. Otherwise, the
           statistic curve is re-calculated, but not plotted.
        fast : bool, optional
           If `True` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is `False`.
        min : pair of numbers, optional
           The minimum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : pair of number, optional
           The maximum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : pair of int, optional
           The number of steps to use. This is used when `delv` is set
           to `None`.
        delv : pair of number, optional
           The step size for the parameter. Setting this over-rides
           the `nloop` parameter. The default is `None`.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : pair of bool, optional
           Should the step size be logarithmically spaced? The
           default (`False`) is to use a linear grid.
        sigma : sequence of number, optional
           The levels at which to draw the contours. The units are the
           change in significance relative to the starting value,
           in units of sigma.
        levels : sequence of number, optional
           The numeric values at which to draw the contours. This
           over-rides the `sigma` parameter, if set (the default is
           `None`).
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        rproj : sherpa.plot.RegionUncertainty instance
           The fields of this object can be used to re-create the plot
           created by `reg_unc`.

        See Also
        --------
        conf : Estimate the confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        int_proj : Calculate and plot the fit statistic versus fit parameter value.
        int_unc : Calculate and plot the fit statistic versus fit parameter value.
        reg_proj : Plot the statistic value as two parameters are varied.
        reg_unc : Plot the statistic value as two parameters are varied.

        Examples
        --------

        Return the results for the `reg_unc` run:

        >>> reg_unc(src.xpos, src.ypos)
        >>> runc = get_reg_unc()

        Create the data without creating a plot:

        >>> runc = get_reg_unc(pl.gamma, gal.nh, recalc=True)

        Control how the data is created:

        >>> runc = get_reg_unc(pl.gamma, gal.nh, id="src",
                               min=(0.5,0.01), max=(2.5,1),
                               nloop=(51,51), log=(False,True),
                               recalc=True)

        """
        if sherpa.utils.bool_cast(recalc):
            par0 = self._check_par(par0, 'par0')
            par1 = self._check_par(par1, 'par1')
            
            if otherids is None:
                otherids = ()
            ids, fit = self._get_fit(id, otherids)
            self._regunc.prepare(min, max, nloop, delv, fac,
                                 sherpa.utils.bool_cast(log), sigma, levels,
                                 numcores)
            self._regunc.calc(fit,par0,par1)
        return self._regunc

    def _int_plot(self, plotobj, par, **kwargs):
        prepare_dict = sherpa.utils.get_keyword_defaults(plotobj.prepare)
        plot_dict = sherpa.utils.get_keyword_defaults(plotobj.plot)
        for key in kwargs.keys():
            if prepare_dict.has_key(key):
                prepare_dict[key] = kwargs[key]
            if plot_dict.has_key(key):
                plot_dict[key] = kwargs[key]
                
        if sherpa.utils.bool_cast(kwargs['replot']):
            self._plot(id, plotobj, replot=True, **plot_dict)
            return
        
        par = self._check_par(par)
        if kwargs['otherids'] is None:
            kwargs['otherids'] = ()
        ids, fit = self._get_fit(kwargs['id'], kwargs['otherids'])
        prepare_dict['log'] = sherpa.utils.bool_cast(prepare_dict['log'])
        plotobj.prepare(**prepare_dict)
        plotobj.calc(fit,par,self._methods)
        # replot but have calculated already, differing interfaces to prepare
        self._plot(id, plotobj, replot=True, **plot_dict)

    
    ### Ahelp ingest: 2015-05-07 DJB
    ### DOC-NOTE: I am not convinced I have fac described correctly
    ### DOC-NOTE: same synopsis as int_unc
    def int_proj(self, par, id=None, otherids=None, replot=False, fast=True,
                 min=None, max=None, nloop=20, delv=None, fac=1, log=False,
                 numcores=None, overplot=False):
        """Calculate and plot the fit statistic versus fit parameter value.

        Create a confidence plot of the fit statistic as a function of
        parameter value. Dashed lines are added to indicate the
        current statistic value and the parameter value at this
        point. The parameter value is varied over a grid of points and
        the free parameters re-fit. It is expected that this is run
        after a successful fit, so that the parameter values are at
        the best-fit location.

        Parameters
        ----------
        par
           The parameter to plot.
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `int_proj`. The default is `False`.
        fast : bool, optional
           If `True` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is `False`.
        min : number, optional
           The minimum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : number, optional
           The maximum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : int, optional
           The number of steps to use. This is used when `delv` is set
           to `None`.
        delv : number, optional
           The step size for the parameter. Setting this over-rides
           the `nloop` parameter. The default is `None`.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : bool, optional
           Should the step size be logarithmically spaced? The
           default (`False`) is to use a linear grid.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        See Also
        --------
        conf : Estimate the confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        get_int_proj : Return the interval-projection object.
        int_unc : Calculate and plot the fit statistic versus fit parameter value.
        reg_proj : Plot the statistic value as two parameters are varied.

        Notes
        -----
        The difference to `int_unc` is that at each step, a fit is
        made to the remaining thawed parameters in the source
        model. This makes the result a more-accurate rendering of the
        projected shape of the hypersurface formed by the statistic,
        but the run-time is longer than, the results of `int_unc`,
        which does not vary any other parameter. If there are no free
        parameters in the source expression, other than the parameter
        being plotted, then the results will be the same.

        Examples
        --------

        Vary the `gamma` parameter of the `p1` model component for
        all data sets with a source expression.

        >>> int_proj(p1.gamma)

        Use only the data in data set 1:

        >>> int_proj(p1.gamma, id=1)

        Use two data sets ('obs1' and 'obs2'):

        >>> int_proj(clus.kt, id='obs1', otherids=['obs2'])

        Vary the `bgnd.c0` parameter between 1e-4 and 2e-4,
        using 41 points:

        >>> int_proj(bgnd.c0, min=1e-4, max=2e-4, step=41)

        This time define the step size, rather than the number of
        steps to use:

        >>> int_proj(bgnd.c0, min=1e-4, max=2e-4, delv=2e-6)

        Overplot the `int_proj` results for the parameter on top of
        the `int_unc` values:

        >>> int_unc(mdl.xpos)
        >>> int_proj(mdl.xpos, overplot=True)

        """
        self._int_plot(self._intproj, par, id=id, otherids=otherids,
                       replot=replot, fast=fast, min=min, max=max, nloop=nloop,
                       delv=delv, fac=fac, log=log, numcores=numcores, 
                       overplot=overplot)

    ### Ahelp ingest: 2015-05-07 DJB
    ### DOC-NOTE: I am not convinced I have fac described correctly
    ### DOC-NOTE: same synopsis as int_proj
    def int_unc(self, par, id=None, otherids=None, replot=False, min=None,
                 max=None, nloop=20, delv=None, fac=1, log=False,
                 numcores=None, overplot=False):
        """Calculate and plot the fit statistic versus fit parameter value.

        Create a confidence plot of the fit statistic as a function of
        parameter value. Dashed lines are added to indicate the
        current statistic value and the parameter value at this
        point. The parameter value is varied over a grid of points and
        the statistic evaluated while holding the other parameters
        fixed. It is expected that this is run after a successful fit,
        so that the parameter values are at the best-fit location.

        Parameters
        ----------
        par
           The parameter to plot.
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `int_proj`. The default is `False`.
        min : number, optional
           The minimum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : number, optional
           The maximum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : int, optional
           The number of steps to use. This is used when `delv` is set
           to `None`.
        delv : number, optional
           The step size for the parameter. Setting this over-rides
           the `nloop` parameter. The default is `None`.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : bool, optional
           Should the step size be logarithmically spaced? The
           default (`False`) is to use a linear grid.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        See Also
        --------
        conf : Estimate the confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        get_int_unc : Return the interval-uncertainty object.
        int_proj : Calculate and plot the fit statistic versus fit parameter value.
        reg_unc : Plot the statistic value as two parameters are varied.

        Notes
        -----
        The difference to `int_proj` is that at each step *only* the
        single parameter value is varied while all other parameters
        remain at their starting value. This makes the result a
        less-accurate rendering of the projected shape of the
        hypersurface formed by the statistic, but the run-time is
        likely shorter than, the results of `int_proj`, which fits the
        model to the remaining thawed parameters at each step. If
        there are no free parameters in the source expression, other
        than the parameter being plotted, then the results will be the
        same.

        Examples
        --------

        Vary the `gamma` parameter of the `p1` model component for
        all data sets with a source expression.

        >>> int_unc(p1.gamma)

        Use only the data in data set 1:

        >>> int_unc(p1.gamma, id=1)

        Use two data sets ('obs1' and 'obs2'):

        >>> int_unc(clus.kt, id='obs1', otherids=['obs2'])

        Vary the `bgnd.c0` parameter between 1e-4 and 2e-4,
        using 41 points:

        >>> int_unc(bgnd.c0, min=1e-4, max=2e-4, step=41)

        This time define the step size, rather than the number of
        steps to use:

        >>> int_unc(bgnd.c0, min=1e-4, max=2e-4, delv=2e-6)

        Overplot the `int_unc` results for the parameter on top of
        the `int_proj` values:

        >>> int_proj(mdl.xpos)
        >>> int_unc(mdl.xpos, overplot=True)

        """
        self._int_plot(self._intunc, par, id=id, otherids=otherids,
                       replot=replot, min=min, max=max, nloop=nloop,
                       delv=delv, fac=fac, log=log, numcores=numcores, 
                       overplot=overplot)

    def _reg_plot(self, plotobj, par0, par1, **kwargs):
        prepare_dict = sherpa.utils.get_keyword_defaults(plotobj.prepare)
        cont_dict = sherpa.utils.get_keyword_defaults(plotobj.contour)
        for key in kwargs.keys():
            if prepare_dict.has_key(key):
                prepare_dict[key] = kwargs[key]
            if cont_dict.has_key(key):
                cont_dict[key] = kwargs[key]
                
        if sherpa.utils.bool_cast(kwargs['replot']):
            self._contour(id, plotobj, replot=True, **cont_dict)
            return

        par0 = self._check_par(par0, 'par0')
        par1 = self._check_par(par1, 'par1')
        
        if kwargs['otherids'] is None:
            kwargs['otherids'] = ()
        ids, fit = self._get_fit(kwargs['id'], kwargs['otherids'])
        prepare_dict['log'] = sherpa.utils.bool_cast(prepare_dict['log'])
        plotobj.prepare(**prepare_dict)
        plotobj.calc(fit,par0,par1,self._methods)
        # replot but have calculated already, differing interfaces to prepare
        self._contour(id, plotobj, replot=True, **cont_dict)


    ### Ahelp ingest: 2015-05-07 DJB
    ### DOC-TODO: how is sigma converted into delta_stat
    def reg_proj(self, par0, par1, id=None, otherids=None, replot=False,
                 fast=True, min=None, max=None, nloop=(10,10), delv=None, fac=4,
                 log=(False,False), sigma=(1,2,3), levels=None, numcores=None, 
                 overplot=False):
        """Plot the statistic value as two parameters are varied.

        Create a confidence plot of the fit statistic as a function of
        parameter value. Dashed lines are added to indicate the
        current statistic value and the parameter value at this
        point. The parameter value is varied over a grid of points and
        the free parameters re-fit. It is expected that this is run
        after a successful fit, so that the parameter values are at
        the best-fit location.

        Parameters
        ----------
        par0, par1
           The parameters to plot on the X and Y axes, respectively.
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `int_proj`. The default is `False`.
        fast : bool, optional
           If `True` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is `False`.
        min : pair of numbers, optional
           The minimum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : pair of number, optional
           The maximum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : pair of int, optional
           The number of steps to use. This is used when `delv` is set
           to `None`.
        delv : pair of number, optional
           The step size for the parameter. Setting this over-rides
           the `nloop` parameter. The default is `None`.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : pair of bool, optional
           Should the step size be logarithmically spaced? The
           default (`False`) is to use a linear grid.
        sigma : sequence of number, optional
           The levels at which to draw the contours. The units are the
           change in significance relative to the starting value,
           in units of sigma.
        levels : sequence of number, optional
           The numeric values at which to draw the contours. This
           over-rides the `sigma` parameter, if set (the default is
           `None`).
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        See Also
        --------
        conf : Estimate the confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        get_reg_proj : Return the interval-projection object.
        int_proj : Calculate and plot the fit statistic versus fit parameter value.
        reg_unc : Plot the statistic value as two parameters are varied.

        Notes
        -----
        The difference to `reg_unc` is that at each step, a fit is
        made to the remaining thawed parameters in the source
        model. This makes the result a more-accurate rendering of the
        projected shape of the hypersurface formed by the statistic,
        but the run-time is longer than, the results of `reg_unc`,
        which does not vary any other parameter. If there are no
        free parameters in the model, other than the parameters
        being plotted, then the results will be the same.

        Examples
        --------

        Vary the `xpos` and `ypos` parameters of the `gsrc` model
        component for all data sets with a source expression.

        >>> reg_proj(gsrc.xpos, gsrc.ypos)

        Use only the data in data set 1:

        >>> reg_proj(gsrc.xpos, gsrc.ypos, id=1)

        Only display the one- and three-sigma contours:

        >>> reg_proj(gsrc.xpos, gsrc.ypos, sigma=(1,3))

        Display contours at values of 5, 10, and 20 more than the
        statistic value of the source model for data set 1:

        >>> s0 = calc_stat(id=1)
        >>> lvls = s0 + np.asarray([5, 10, 20])
        >>> reg_proj(gsrc.xpos, gsrc.ypos, levels=lvls, id=1)

        Increase the limits of the plot and the number of steps along
        each axis:

        >>> reg_proj(gsrc.xpos, gsrc.ypos, id=1, fac=6, nloop=(41,41))

        Compare the `ampl` parameters of the `g` and `b` model
        components, for data sets 'core' and 'jet', over the given
        ranges:

        >>> reg_proj(g.ampl, b.ampl, min=(0,1e-4), max=(0.2,5e-4),
                     nloop=(51,51), id='core', otherids=['jet'])

        """
        self._reg_plot(self._regproj, par0, par1, id=id, otherids=otherids,
                       replot=replot, fast=fast, min=min, max=max, nloop=nloop,
                       delv=delv, fac=fac, log=log, sigma=sigma, levels=levels,
                       numcores=numcores, overplot=overplot)

    ### Ahelp ingest: 2015-05-07 DJB
    ### DOC-TODO: how is sigma converted into delta_stat
    def reg_unc(self, par0, par1, id=None, otherids=None, replot=False,
                min=None, max=None, nloop=(10,10), delv=None, fac=4,
                log=(False,False), sigma=(1,2,3), levels=None, numcores=None, 
                overplot=False):
        """Plot the statistic value as two parameters are varied.

        Create a confidence plot of the fit statistic as a function of
        parameter value. Dashed lines are added to indicate the
        current statistic value and the parameter value at this
        point. The parameter value is varied over a grid of points and
        the statistic evaluated while holding the other parameters
        fixed. It is expected that this is run after a successful fit,
        so that the parameter values are at the best-fit location.

        Parameters
        ----------
        par0, par1
           The parameters to plot on the X and Y axes, respectively.
        id : str or int, optional
        otherids : list of str or int, optional
           The `id` and `otherids` arguments determine which data set
           or data sets are used. If not given, all data sets which
           have a defined source model are used.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `int_proj`. The default is `False`.
        min : pair of numbers, optional
           The minimum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : pair of number, optional
           The maximum parameter value for the calcutation. The
           default value of `None` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : pair of int, optional
           The number of steps to use. This is used when `delv` is set
           to `None`.
        delv : pair of number, optional
           The step size for the parameter. Setting this over-rides
           the `nloop` parameter. The default is `None`.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : pair of bool, optional
           Should the step size be logarithmically spaced? The
           default (`False`) is to use a linear grid.
        sigma : sequence of number, optional
           The levels at which to draw the contours. The units are the
           change in significance relative to the starting value,
           in units of sigma.
        levels : sequence of number, optional
           The numeric values at which to draw the contours. This
           over-rides the `sigma` parameter, if set (the default is
           `None`).
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        See Also
        --------
        conf : Estimate the confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        get_reg_unc : Return the interval-uncertainty object.
        int_unc : Calculate and plot the fit statistic versus fit parameter value.
        reg_proj : Plot the statistic value as two parameters are varied.

        Notes
        -----

        The difference to `reg_proj` is that at each step *only* the
        pair of parameters are varied, while all the other parameters
        remain at their starting value. This makes the result a
        less-accurate rendering of the projected shape of the
        hypersurface formed by the statistic, but the run-time is
        likely shorter than, the results of `reg_proj`, which fits the
        model to the remaining thawed parameters at each step. If
        there are no free parameters in the model, other than the
        parameters being plotted, then the results will be the same.

        Examples
        --------

        Vary the `xpos` and `ypos` parameters of the `gsrc` model
        component for all data sets with a source expression.

        >>> reg_unc(gsrc.xpos, gsrc.ypos)

        Use only the data in data set 1:

        >>> reg_unc(gsrc.xpos, gsrc.ypos, id=1)

        Only display the one- and three-sigma contours:

        >>> reg_unc(gsrc.xpos, gsrc.ypos, sigma=(1,3))

        Display contours at values of 5, 10, and 20 more than the
        statistic value of the source model for data set 1:

        >>> s0 = calc_stat(id=1)
        >>> lvls = s0 + np.asarray([5, 10, 20])
        >>> reg_unc(gsrc.xpos, gsrc.ypos, levels=lvls, id=1)

        Increase the limits of the plot and the number of steps along
        each axis:

        >>> reg_unc(gsrc.xpos, gsrc.ypos, id=1, fac=6, nloop=(41,41))

        Compare the `ampl` parameters of the `g` and `b` model
        components, for data sets 'core' and 'jet', over the given
        ranges:

        >>> reg_unc(g.ampl, b.ampl, min=(0,1e-4), max=(0.2,5e-4),
                    nloop=(51,51), id='core', otherids=['jet'])

        Overplot the results on the `reg_proj` plot:

        >>> reg_proj(s1.c0, s2.xpos)
        >>> reg_unc(s1.c0, s2.xpos, overplot=True)

        """
        self._reg_plot(self._regunc, par0, par1, id=id, otherids=otherids,
                       replot=replot, min=min, max=max, nloop=nloop, delv=delv,
                       fac=fac, log=log, sigma=sigma, levels=levels, 
                       numcores=numcores, overplot=overplot)


    # Aliases
    #interval_projection = int_proj
    #interval_uncertainty = int_unc
    #region_projection = reg_proj
    #region_uncertainty = reg_unc

    ###########################################################################
    # Basic imaging
    ###########################################################################


    def _prepare_imageobj(self, id, imageobj, model=None):
        if( isinstance(imageobj, sherpa.image.ComponentModelImage) or
            isinstance(imageobj, sherpa.image.ComponentSourceImage)):
            imageobj.prepare_image(self.get_data(id), model)
        elif (isinstance(imageobj, sherpa.image.PSFImage) or
              isinstance(imageobj, sherpa.image.PSFKernelImage)):
            imageobj.prepare_image(self.get_psf(id), self.get_data(id))
        elif isinstance(imageobj, sherpa.image.DataImage):
            imageobj.prepare_image(self.get_data(id))
        elif isinstance(imageobj, sherpa.image.SourceImage):
            imageobj.prepare_image(self.get_data(id), self.get_source(id))
        else:
            imageobj.prepare_image(self.get_data(id), self.get_model(id))
        return imageobj

    #
    # Image object access
    #

    ### Ahelp ingest: 2015-05-12 DJB
    def get_data_image(self, id=None):
        """Return the data used by image_data.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        data_img : a sherpa.image.DataImage instance
           The `y` attribute contains the ratio values as a 2D NumPy
           array.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set is not 2D.
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        contour_data : Contour the values of an image data set.
        image_data : Display a data set in the image viewer.

        Examples
        --------

        Return the image data for the default data set:

        >>> dinfo = get_data_image()
        >>> dinfo.y.shape
        (150, 175)

        """
        self._prepare_imageobj(id, self._dataimage)
        return self._dataimage
    
    ### Ahelp ingest: 2015-05-12 DJB
    def get_model_image(self, id=None):
        """Return the data used by image_model.

        Evaluate the source expression for the image pixels -
        including any PSF convolution defined by `set_psf` - and
        return the results.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        src_img : a sherpa.image.ModelImage instance
           The `y` attribute contains the source model values as a 2D
           NumPy array.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_source_image : Return the data used by image_source.
        contour_model : Contour the values of the model, including any PSF.
        image_model : Display the model for a data set in the image viewer.
        set_psf : Add a PSF model to a data set.

        Examples
        --------

        Calculate the residuals (data - model) for the default
        data set:

        >>> minfo = get_model_image()
        >>> dinfo = get_data_image()
        >>> resid = dinfo.y - minfo.y

        """
        self._prepare_imageobj(id, self._modelimage)
        return self._modelimage

    ### Ahelp ingest: 2015-05-12 DJB
    ### DOC-TODO: it looks like get_source_image doesn't raise DataErr with
    ###           a non-2D data set
    def get_source_image(self, id=None):
        """Return the data used by image_source.

        Evaluate the source expression for the image pixels - without
        any PSF convolution - and return the results.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        src_img : a sherpa.image.SourceImage instance
           The `y` attribute contains the source model values as a 2D
           NumPy array.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_model_image : Return the data used by image_model.
        contour_source : Contour the values of the model, without any PSF.
        image_source : Display the source expression for a data set in the image viewer.

        Examples
        --------

        Return the model data for the default data set:

        >>> sinfo = get_source_image()
        >>> sinfo.y.shape
        (150, 175)

        """
        self._prepare_imageobj(id, self._sourceimage)
        return self._sourceimage


    ### Ahelp ingest: 2015-05-12 DJB
    def get_model_component_image(self, id, model=None):
        """Return the data used by image_model_component.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to display (the name, if a string).

        Returns
        -------
        cpt_img : a sherpa.image.ComponentModelImage instance
           The `y` attribute contains the component model values as a
           2D NumPy array.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_source_component_image : Return the data used by image_source_component.
        get_model_image : Return the data used by image_model.
        image_model : Display the model for a data set in the image viewer.
        image_model_component : Display a component of the model in the image viewer.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        Examples
        --------

        Return the gsrc component values for the default data set:

        >>> minfo = get_model_component_image(gsrc)

        Get the `bgnd` model pixel values for data set 2:

        >>> minfo = get_model_component_image(2, bgnd)

        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        self._prepare_imageobj(id, self._mdlcompimage, model=model)
        return self._mdlcompimage

    ### Ahelp ingest: 2015-05-12 DJB
    def get_source_component_image(self, id, model=None):
        """Return the data used by image_source_component.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to display (the name, if a string).

        Returns
        -------
        cpt_img : a sherpa.image.ComponentSourceImage instance
           The `y` attribute contains the component model values as a
           2D NumPy array.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_model_component_image : Return the data used by image_model_component.
        get_source_image : Return the data used by image_source.
        image_source : Display the source expression for a data set in the image viewer.
        image_source_component : Display a component of the source expression in the image viewer.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        Examples
        --------

        Return the gsrc component values for the default data set:

        >>> sinfo = get_source_component_image(gsrc)

        Get the `bgnd` model pixel values for data set 2:

        >>> sinfo = get_source_component_image(2, bgnd)

        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        self._prepare_imageobj(id, self._srccompimage, model=model)
        return self._srccompimage

    ### Ahelp ingest: 2015-05-11 DJB
    def get_ratio_image(self, id=None):
        """Return the data used by image_ratio.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        ratio_img : a sherpa.image.RatioImage instance
           The `y` attribute contains the ratio values as a 2D NumPy
           array.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set is not 2D.
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_resid_image : Return the data used by image_resid.
        contour_ratio : Contour the ratio of data to model.
        image_ratio : Display the ratio (data/model) for a data set in the image viewer.

        Examples
        --------

        Return the ratio data for the default data set:

        >>> rinfo = get_ratio_image()

        """
        self._prepare_imageobj(id, self._ratioimage)
        return self._ratioimage
    
    ### Ahelp ingest: 2015-05-11 DJB
    def get_resid_image(self, id=None):
        """Return the data used by image_resid.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        resid_img : a sherpa.image.ResidImage instance
           The `y` attribute contains the residual values as a 2D
           NumPy array.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set is not 2D.
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_ratio_image : Return the data used by image_ratio.
        contour_resid : Contour the residuals of the fit.
        image_resid : Display the residuals (data - model) for a data set in the image viewer.

        Examples
        --------

        Return the residual data for the default data set:

        >>> rinfo = get_resid_image()

        """
        self._prepare_imageobj(id, self._residimage)
        return self._residimage

    ### Ahelp ingest: 2015-05-11 DJB
    def get_psf_image(self, id=None):
        """Return the data used by image_psf.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        psf_data : a sherpa.image.PSFImage instance

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If a PSF model has not been created for the data set.

        See Also
        --------
        get_kernel_image : Return the data used by image_kernel.
        image_kernel : Display the 2D kernel for a data set in the image viewer.
        image_psf : Display the 2D PSF model for a data set in the image viewer.

        Examples
        --------

        Return the image data for the PSF for the default data set:

        >>> iplot = get_psf_image()
        >>> iplot.y.shape
        (175, 200)

        """
        self._prepare_imageobj(id, self._psfimage)
        return self._psfimage


    ### Ahelp ingest: 2015-05-11 DJB
    def get_kernel_image(self, id=None):
        """Return the data used by image_kernel.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        psf_data : a sherpa.image.PSFKernelImage instance

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If a PSF model has not been created for the data set.

        See Also
        --------
        get_psf_image : Return the data used by image_psf.
        image_kernel : Display the 2D kernel for a data set in the image viewer.
        image_psf : Display the 2D PSF model for a data set in the image viewer.

        Examples
        --------

        Return the image data for the kernel for the default data set:

        >>> lplot = get_kernel_image()
        >>> iplot.y.shape
        (51, 51)

        """
        self._prepare_imageobj(id, self._kernelimage)
        return self._kernelimage

    #
    # Images
    #

    def _image(self, id, imageobj, shape, newframe, tile, model=None):
        self._prepare_imageobj(id, imageobj, model).image(shape, newframe, tile)

    ### Ahelp ingest: 2015-05-09 DJB
    def image_data(self, id=None, newframe=False, tile=False):
        """Display a data set in the image viewer.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If `False`, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If `False`, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_data_image :
        image_close : Close the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_open : Open the image viewer.
        image_source : Display the model for a data set in the image viewer.

        Notes
        -----
        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        Display the data in default data set.

        >>> image_data()

        Display data set 2 in a new frame so that the data in the
        current frame is not destroyed. The new data will be displayed
        in a single frame (i.e. the only data shown by the viewer).

        >>> image_data(2, newframe=True)

        Display data sets 'i1' and 'i2' side by side:

        >>> image_data('i1')
        >>> image_data('i2', newframe=True, tile=True)

        """
        self._image(id, self._dataimage, None,
                    newframe, tile)

    ### Ahelp ingest: 2015-05-09 DJB
    def image_model(self, id=None, newframe=False, tile=False):
        """Display the model for a data set in the image viewer.

        This function evaluates and displays the model expression for
        a data set, including any instrument response (e.g. PSF or ARF
        and RMF) whether created automatically or with
        `set_full_model`.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If `False`, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If `False`, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_model_image :
        image_close : Close the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_model_component : Display a component of the model in the image viewer.
        image_open : Open the image viewer.
        image_source : Display the model for a data set in the image viewer.
        image_source_component : Display a component of the source expression in the image viewer.

        Notes
        -----
        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        Display the model for the default data set.

        >>> image_model()

        Display the model for data set 2 in a new frame so that the
        data in the current frame is not destroyed. The new data will
        be displayed in a single frame (i.e. the only data shown by
        the viewer).

        >>> image_model(2, newframe=True)

        Display the models for data sets 'i1' and 'i2' side by side:

        >>> image_model('i1')
        >>> image_model('i2', newframe=True, tile=True)

        """
        self._image(id, self._modelimage, None,
                    newframe, tile)


    ### Ahelp ingest: 2015-05-11 DJB
    def image_source_component(self, id, model=None, newframe=False,
                               tile=False):
        """Display a component of the source expression in the image viewer.

        This function evaluates and displays a component of the model
        expression for a data set, without any instrument response.
        Use `image_model_component` to include any response.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to display (the name, if a string).
        newframe : bool, optional
           Create a new frame for the data? If `False`, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If `False`, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_source_component_image : Return the data used by image_source_component.
        image_close : Close the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_model : Display the model for a data set in the image viewer.
        image_model_component : Display a component of the model in the image viewer.
        image_open : Open the image viewer.
        image_source : Display the source expression for a data set in the image viewer.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        Display the full source model and then just the `gsrc`
        component for the default data set:

        >>> image_source()
        >>> image_source_component(gsrc)

        Display the 'clus' and 'bgnd' components of the model for the
        'img' data set side by side:

        >>> image_source_component('img', 'clus')
        >>> image_source_component('img', 'bgnd', newframe=True,
                                   tile=True)

        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        self._image(id, self._srccompimage, None, newframe, tile, model=model)


    ### Ahelp ingest: 2015-05-11 DJB
    def image_model_component(self, id, model=None, newframe=False, tile=False):
        """Display a component of the model in the image viewer.

        This function evaluates and displays a component of the model
        expression for a data set, including any instrument response.
        Use `image_source_component` to exclude the response.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to display (the name, if a string).
        newframe : bool, optional
           Create a new frame for the data? If `False`, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If `False`, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_model_component_image :
        image_close : Close the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_model : Display the model for a data set in the image viewer.
        image_open : Open the image viewer.
        image_source : Display the source expression for a data set in the image viewer.
        image_source_component : Display a component of the source expression in the image viewer.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        Display the full source model and then just the `gsrc`
        component for the default data set:

        >>> image_model()
        >>> image_model_component(gsrc)

        Display the 'clus' component of the model for the 'img' data
        set side by side without the with any instrument response
        (such as convolution with a PSF model):

        >>> image_source_component('img', 'clus')
        >>> image_model_component('img', 'clus', newframe=True,
                                  tile=True)

        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        is_source = self._get_model_status(id)[1]
        model = self._add_convolution_models(id, self.get_data(id),
                                             model, is_source)

        self._image(id, self._mdlcompimage, None, newframe, tile, model=model)


    ### Ahelp ingest: 2015-05-09 DJB
    def image_source(self, id=None, newframe=False, tile=False):
        """Display the source expression for a data set in the image viewer.

        This function evaluates and displays the model expression for
        a data set, without any instrument response.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If `False`, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If `False`, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_source_image : Return the data used by image_source.
        image_close : Close the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_model : Display the model for a data set in the image viewer.
        image_model_component : Display a component of the model in the image viewer.
        image_open : Open the image viewer.
        image_source_component : Display a component of the source expression in the image viewer.

        Notes
        -----
        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        Display the source model for the default data set.

        >>> image_source()

        Display the source model for data set 2 in a new frame so that
        the data in the current frame is not destroyed. The new data
        will be displayed in a single frame (i.e. the only data shown
        by the viewer).

        >>> image_source(2, newframe=True)

        Display the source models for data sets 'i1' and 'i2' side by
        side:

        >>> image_source('i1')
        >>> image_source('i2', newframe=True, tile=True)

        """
        self._image(id, self._sourceimage, None, newframe, tile)

    ### Ahelp ingest: 2015-05-09 DJB
    ### DOC-TODO: does newframe make sense here?
    def image_fit(self, id=None, newframe=True, tile=True, deleteframes=True):
        """Display the data, model, and residuals for a data set in the image viewer.

        This function displays the data, model (including any
        instrument response), and the residuals (data - model), for a
        data set.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If `False`, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If `False`, the default, then
           only a single frame is displayed.
        deleteframes : bool, optional
           Should existing frames be deleted? The default is `True`.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        image_close : Close the image viewer.
        image_data : Display a data set in the image viewer.
        image_model : Display the model for a data set in the image viewer.
        image_open : Open the image viewer.
        image_resid : Display the residuals (data - model) for a data set in the image viewer.

        Notes
        -----
        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        Display the fit results - that is, the data, model, and
        residuals - for the default data set.

        >>> image_fit()

        Do not tile the frames (the three frames are loaded, but only
        the last displayed, the residuals), and then change the frame
        being displayed to the second one (the model).

        >>> image_fit('img', tile=False)
        >>> image_xpaset('frame 2')

        """
        self._prepare_imageobj(id, self._dataimage)
        self._prepare_imageobj(id, self._modelimage)
        self._prepare_imageobj(id, self._residimage)
        deleteframes=sherpa.utils.bool_cast(deleteframes)
        if deleteframes is True:
            sherpa.image.Image.open()
            sherpa.image.Image.delete_frames()
        self._dataimage.image(None, False, tile)
        self._modelimage.image(None, newframe, tile)
        self._residimage.image(None, newframe, tile)

    ### Ahelp ingest: 2015-05-09 DJB
    def image_resid(self, id=None, newframe=False, tile=False):
        """Display the residuals (data - model) for a data set in the image viewer.

        This function displays the residuals (data - model) for a data
        set.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If `False`, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If `False`, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_resid_image : Return the data used by image_resid.
        image_close : Close the image viewer.
        image_data : Display a data set in the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_model : Display the model for a data set in the image viewer.
        image_open : Open the image viewer.
        image_ratio : Display the ratio (data/model) for a data set in the image viewer.
        image_source : Display the model for a data set in the image viewer.

        Notes
        -----
        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        Display the residuals for the default data set.

        >>> image_resid()

        Display the residuals for data set 2 in a new frame so that
        the data in the current frame is not destroyed. The new data
        will be displayed in a single frame (i.e. the only data shown
        by the viewer).

        >>> image_resid(2, newframe=True)

        Display the residuals for data sets 'i1' and 'i2' side by
        side:

        >>> image_resid('i1')
        >>> image_resid('i2', newframe=True, tile=True)

        """
        self._image(id, self._residimage, None,
                    newframe, tile)

    ### Ahelp ingest: 2015-05-09 DJB
    def image_ratio(self, id=None, newframe=False, tile=False):
        """Display the ratio (data/model) for a data set in the image viewer.

        This function displays the ratio data/model for a data
        set.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If `False`, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If `False`, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_ratio_image : Return the data used by image_ratio.
        image_close : Close the image viewer.
        image_data : Display a data set in the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_model : Display the model for a data set in the image viewer.
        image_open : Open the image viewer.
        image_resid : Display the residuals (data - model) for a data set in the image viewer.
        image_source : Display the model for a data set in the image viewer.

        Notes
        -----
        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        Display the ratio (data/model) for the default data set.

        >>> image_ratio()

        """
        self._image(id, self._ratioimage, None,
                    newframe, tile)

    ### Ahelp ingest: 2015-05-11 DJB
    ### DOC-TODO: what gets displayed when there is no PSF?
    def image_psf(self, id=None, newframe=False, tile=False):
        """Display the 2D PSF model for a data set in the image viewer.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If `False`, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If `False`, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_psf_image :
        image_close : Close the image viewer.
        image_data : Display a data set in the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_model : Display the model for a data set in the image viewer.
        image_open : Open the image viewer.
        image_source : Display the model for a data set in the image viewer.
        plot_psf : Plot the 1D PSF model applied to a data set.

        Notes
        -----
        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        >>> image_psf()

        >>> image_psf(2)

        """
        self._image(id, self._psfimage, None, newframe, tile)


    ### Ahelp ingest: 2015-05-11 DJB
    ### DOC-TODO: what gets displayed when there is no PSF?
    ### DOC-TODO: where to point to for PSF/kernel discussion/description
    ###           (as it appears in a number of places)?
    def image_kernel(self, id=None, newframe=False, tile=False):
        """Display the 2D kernel for a data set in the image viewer.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If `False`, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If `False`, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_kernel_image :
        image_close : Close the image viewer.
        image_data : Display a data set in the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_model : Display the model for a data set in the image viewer.
        image_open : Open the image viewer.
        image_source : Display the model for a data set in the image viewer.
        plot_kernel : Plot the 1D kernel applied to a data set.

        Notes
        -----
        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        >>> image_kernel()

        >>> image_kernel(2)

        """
        self._image(id, self._kernelimage, None, newframe, tile)

    # Manage these functions (open, close, delete frames, regions, XPA)
    # through unbound functions of the Image class--always talking to
    # the same instance of the Image backend, so OK for now

    ### Ahelp ingest: 2015-05-08 DJB
    def image_deleteframes(self):
        """Delete all the frames open in the image viewer.

        Delete all the frames - in other words, images - being
        displayed in the image viewer (e.g. as created by
        `image_data` or `image_fit`).

        See Also
        --------
        image_close : Close the image viewer.
        image_getregion : Return the region defined in the image viewer.
        image_open : Create the image viewer.
        image_setregion : Set the region to display in the image viewer.
        image_xpaget : Return the result of an XPA call to the image viewer.
        image_xpaset : Send an XPA command to the image viewer.

        Examples
        --------

        >>> image_deleteframes()

        """
        sherpa.image.Image.delete_frames()
        
    ### Ahelp ingest: 2015-05-08 DJB
    def image_open(self):
        """Start the image viewer.

        The image viewer will be started, if found. Calling this
        function when the viewer has already been started will not
        cause a second viewer to be started. The image viewer
        will be started automatically by any of the commands
        like `image_data`.

        See Also
        --------
        image_close : Close the image viewer.
        image_deleteframes : Delete all the frames open in the image viewer.
        image_getregion : Return the region defined in the image viewer.
        image_setregion : Set the region to display in the image viewer.
        image_xpaget : Return the result of an XPA call to the image viewer.
        image_xpaset : Send an XPA command to the image viewer.

        Notes
        -----
        Image visualization is optional, and provided by the
        DS9 application [1]_.

        References
        ----------

        .. [1] http://ds9.si.edu/site/Home.html

        Examples
        --------

        >>> image_open()

        """
        sherpa.image.Image.open()

    ### Ahelp ingest: 2015-05-08 DJB
    def image_close(self):
        """Close the image viewer.

        Close the image viewer created by a previous call to one
        of the `image_xxx` functions.

        See Also
        --------
        image_deleteframes : Delete all the frames open in the image viewer.
        image_getregion : Return the region defined in the image viewer.
        image_open : Start the image viewer.
        image_setregion : Set the region to display in the image viewer.
        image_xpaget : Return the result of an XPA call to the image viewer.
        image_xpaset : Send an XPA command to the image viewer.

        Examples
        --------

        >>> image_close()

        """
        sherpa.image.Image.close()

    ### Ahelp ingest: 2015-05-11 DJB
    ### DOC-TODO: what is the "default" coordinate system
    def image_getregion(self, coord=''):
        """Return the region defined in the image viewer.

        The regions defined in the current frame are returned.

        Parameters
        ----------
        coord : str, optional
           The coordinate system to use.

        Returns
        -------
        region : str
           The region, or regions, or the empty string.

        Raises
        ------
        sherpa.utils.err.DS9Err
           Invalid coordinate system.

        See Also
        --------
        image_setregion : Set the region to display in the image viewer.
        image_xpaget : Return the result of an XPA call to the image viewer.
        image_xpaset : Send an XPA command to the image viewer.

        Examples
        --------

        >>> image_getregion()
        'circle(123,128,12.377649);-box(130,121,14,14,329.93142);'

        >>> image_getregion('physical')
        'circle(3920.5,4080.5,396.08476);-rotbox(4144.5,3856.5,448,448,329.93142);'

        """
        return sherpa.image.Image.get_region(coord)

    ### Ahelp ingest: 2015-05-11 DJB
    ### DOC-TODO: what is the "default" coordinate system
    def image_setregion(self, reg, coord=''):
        """Set the region to display in the image viewer.

        Parameters
        ----------
        reg : str
           The region to display.
        coord : str, optional
           The coordinate system to use.

        Raises
        ------
        sherpa.utils.err.DS9Err
           Invalid coordinate system.

        See Also
        --------
        image_getregion : Return the region defined in the image viewer.
        image_xpaget : Return the result of an XPA call to the image viewer.
        image_xpaset : Send an XPA command to the image viewer.

        Examples
        --------

        Add a circle, in the physical coordinate system, to the data
        from the default data set:

        >>> image_data()
        >>> image_setregion('circle(4234.53,3245.29,46.74)', 'physical')

        Copy the region from the current frame, create a new frame
        displaying the residuals from data set 'img', and then display
        the region on it:

        >>> r = image_getregion()
        >>> image_resid('img', newframe=True)
        >>> image_setregion(r)

        """
        sherpa.image.Image.set_region(reg, coord)

    ### Ahelp ingest: 2015-05-09 DJB
    ### DOC-TODO: check the ds9 link when it is working
    ### DOC-TODO: is there a link of ds9 commands we can link to?
    def image_xpaget(self, arg):
        """Return the result of an XPA call to the image viewer.

        Send a query to the image viewer.

        Parameters
        ----------
        arg : str
           A command to send to the image viewer via XPA.

        Returns
        -------
        returnval : str

        Raises
        ------
        sherpa.utils.err.DS9Err
           The image viewer is not running.
        sherpa.utils.err.RuntimeErr
           If the command is not recognized.

        See Also
        --------
        image_close : Close the image viewer.
        image_getregion : Return the region defined in the image viewer.
        image_open : Create the image viewer.
        image_setregion : Set the region to display in the image viewer.
        image_xpaset : Send an XPA command to the image viewer.

        Notes
        -----

        This XPA access point [1]_ of the ds9 image viewer lets
        commands and queries to be sent to the viewer.

        References
        ----------

        .. [1] http://ds9.si.edu/ref/xpa.html

        Examples
        --------

        Return the current zoom setting of the active frame:

        >>> image_xpaget('zoom')
        '1\n'

        """
        return sherpa.image.Image.xpaget(arg)

    ### Ahelp ingest: 2015-05-09 DJB
    ### DOC-TODO: check the ds9 link when it is working
    ### DOC-TODO: is there a link of ds9 commands we can link to?
    def image_xpaset(self, arg, data=None):
        """Return the result of an XPA call to the image viewer.

        Send a command to the image viewer.

        Parameters
        ----------
        arg : str
           A command to send to the image viewer via XPA.
        data : optional
           The data for the command.

        Raises
        ------
        sherpa.utils.err.DS9Err
           The image viewer is not running.
        sherpa.utils.err.RuntimeErr
           If the command is not recognized or could not be completed.

        See Also
        --------
        image_close : Close the image viewer.
        image_getregion : Return the region defined in the image viewer.
        image_open : Create the image viewer.
        image_setregion : Set the region to display in the image viewer.
        image_xpaset : Send an XPA command to the image viewer.

        Notes
        -----

        This XPA access point [1]_ of the ds9 image viewer lets
        commands and queries to be sent to the viewer.

        References
        ----------

        .. [1] http://ds9.si.edu/ref/xpa.html

        Examples
        --------

        Change the zoom setting of the active frame:

        >>> image_xpaset('zoom 4')

        Overlay the coordinate grid on the current frame:

        >>> image_xpaset('grid yes')

        Add the region file `src.reg` to the display:

        >>> image_xpaset('regions src.reg')

        Create a png version of the image being displayed:

        >>> image_xpaset('saveimage png /tmp/img.png')

        """
        return sherpa.image.Image.xpaset(arg, data)

#    def log_model_call(self, id, model, func_name):
#	id_ = id
#	if model is None:
#		id_ = None
#	id_ = self._fix_id(id_)
#	if not self._calls_tracker.has_key(id_):
#		self._calls_tracker[id_] = dict()
#	line = readline.get_history_item(readline.get_current_history_length())
#	self._calls_tracker[id_][func_name] = line

