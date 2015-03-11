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

    def clean(self):
        """
        clean

        SYNOPSIS
           Clear all stored session data

        SYNTAX

        Arguments:
           None

        Returns:
           None

        DESCRIPTION
           Clear all stored session data from internal structure.

        SEE ALSO
           save, restore
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


    def save(self, filename='sherpa.save', clobber=False):
        """
        save

        SYNOPSIS
           Save the current Sherpa session to a pickled file

        SYNTAX

        Arguments:
           filename   - name of save file
                        default = 'sherpa.save'           

           clobber    - clobber the file if it exists
                        default = False

        Returns:
           None

        DESCRIPTION
           Save the current Sherpa session information to a pickled file
           to be restored for future use.

        SEE ALSO
           restore, clean
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

    def restore(self, filename='sherpa.save'):
        """
        restore

        SYNOPSIS
           Restore a previous Sherpa session from a pickled file

        SYNTAX

        Arguments:
           filename   - name of saved file
                        default = 'sherpa.save'           

        Returns:
           None

        DESCRIPTION
           Restore previous Sherpa session information from a pickled file
           for continued use.

        SEE ALSO
           save, clean
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

    def show_stat(self, outfile=None, clobber=False):
        """
        show_stat

        SYNOPSIS
           Show the current Sherpa statistic

        SYNTAX

        Arguments:
           outfile   - filename to capture the output
                      default = None

           clobber  - overwrite outfile if exists
                      default = False

        Returns:
           None

        DESCRIPTION
           Show the current Sherpa fit statistic.

           Examples:
              show_stat()

              show_stat("sherpa.stat", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_model,
           show_conf, show_proj, show_data, show_covar
        """
        all = ''
        all += self._get_show_stat()
        _send_to_pager(all, outfile, clobber)


    def show_method(self, outfile=None, clobber=False):
        """
        show_method

        SYNOPSIS
           Show the current Sherpa method

        SYNTAX

        Arguments:
           outfile   - filename to capture the output
                      default = None

           clobber  - overwrite outfile if exists
                      default = False

        Returns:
           None

        DESCRIPTION
           Show the current Sherpa optimization method.

           Examples:
              show_method()

              show_method("sherpa.method", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_model,
           show_conf, show_proj, show_data, show_covar
        """
        all = ''
        all += self._get_show_method()
        _send_to_pager(all, outfile, clobber)


    def show_fit(self, outfile=None, clobber=False):
        """
        show_fit

        SYNOPSIS
           Show results from the last Sherpa fit performed

        SYNTAX

        Arguments:
           outfile  - filename to capture the output
                     default = None

           clobber - overwrite outfile if exists
                     default = False

        Returns:
           None

        DESCRIPTION
           Show results from the last Sherpa fit performed.  Fit results
           can be written to file.

           Examples:
              show_fit()

              show_fit("sherpa.fit", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_model, show_conf,
           show_proj, show_data, show_covar
        """
        all = ''
        all += self._get_show_fit()
        _send_to_pager(all, outfile, clobber)


    def show_data(self, id=None, outfile=None, clobber=False):
        """
        show_data

        SYNOPSIS
           Show Sherpa datasets

        SYNTAX

        Arguments:
           id      - data id
                     default = All available data ids

           outfile  - filename to capture the output
                     default = None

           clobber - overwrite outfile if exists
                     default = False

        Returns:
           None

        DESCRIPTION
           Show all current Sherpa datasets or by Sherpa data id.

           Examples:
              show_data()

              show_data(1)

              show_data(2, "sherpa.data", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_model, show_conf,
           show_proj, show_fit, show_covar
        """
        all = ''
        all += self._get_show_data(id)
        _send_to_pager(all, outfile, clobber)


    def show_filter(self, id=None, outfile=None, clobber=False):
        """
        show_filter

        SYNOPSIS
           Show filters on Sherpa datasets

        SYNTAX

        Arguments:
           id      - data id
                     default = All available data ids

           outfile  - filename to capture the output
                     default = None

           clobber - overwrite outfile if exists
                     default = False

        Returns:
           None

        DESCRIPTION
           Show all current filters on Sherpa datasets or by Sherpa data id.

           Examples:
              show_filter()

              show_filter(1)

              show_filter(2, "sherpa.filter", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_model,
           show_conf, show_proj, show_fit, show_covar, show_data
        """
        all = ''
        all += self._get_show_filter(id)
        _send_to_pager(all, outfile, clobber)


    def show_model(self, id=None, outfile=None, clobber=False):
        """
        show_model

        SYNOPSIS
           Show Sherpa models

        SYNTAX

        Arguments:
           id      - data id
                     default = All available data ids

           outfile  - filename to capture the output
                     default = None

           clobber - overwrite outfile if exists
                     default = False

        Returns:
           None

        DESCRIPTION
           Show all current Sherpa models or by Sherpa data id.

           Examples:
              show_model()

              show_model(1)

              show_model(2, "sherpa.model", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_data,
           show_conf, show_proj, show_fit, show_covar
        """
        all = ''
        #all += self._get_show_kernel(id)
        all += self._get_show_psf(id)
        all += self._get_show_model(id)
        _send_to_pager(all, outfile, clobber)


    def show_source(self, id=None, outfile=None, clobber=False):
        """
        show_source

        SYNOPSIS
           Show Sherpa sources

        SYNTAX

        Arguments:
           id      - data id
                     default = All available data ids

           outfile  - filename to capture the output
                     default = None

           clobber - overwrite outfile if exists
                     default = False

        Returns:
           None

        DESCRIPTION
           Show all current Sherpa source models or by Sherpa data id.

           Examples:
              show_source()

              show_source(1)

              show_source(2, "sherpa.source", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_data,
           show_conf, show_proj, show_fit, show_covar
        """
        all = ''
        all += self._get_show_source(id)
        _send_to_pager(all, outfile, clobber)

    def show_kernel(self, id=None, outfile=None, clobber=False):
        """
        show_kernel

        SYNOPSIS
           Show Sherpa PSF kernels

        SYNTAX

        Arguments:
           id      - data id
                     default = All available data ids

           outfile  - filename to capture the output
                     default = None

           clobber - overwrite outfile if exists
                     default = False

        Returns:
           None

        DESCRIPTION
           Show all current Sherpa PSF kernel models or by Sherpa data id.

           Examples:
              show_kernel()

              show_kernel(1)

              show_kernel(2, "sherpa.kernel", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_data,
           show_conf, show_proj, show_fit, show_covar, show_model, show_source
        """
        all = ''
        all += self._get_show_kernel(id)
        _send_to_pager(all, outfile, clobber)

    def show_psf(self, id=None, outfile=None, clobber=False):
        """
        show_psf

        SYNOPSIS
           Show Sherpa PSF model with PSF kernel

        SYNTAX

        Arguments:
           id      - data id
                     default = All available data ids

           outfile  - filename to capture the output
                     default = None

           clobber - overwrite outfile if exists
                     default = False

        Returns:
           None

        DESCRIPTION
           Show all current Sherpa PSF models with PSF kernels or by Sherpa data id.

           Examples:
              show_psf()

              show_psf(1)

              show_psf(2, "sherpa.psf", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_data,
           show_conf, show_proj, show_fit, show_covar
        """
        all = ''
        #all += self._get_show_kernel(id)
        all += self._get_show_psf(id)
        _send_to_pager(all, outfile, clobber)

    def show_conf(self, outfile=None, clobber=False):
        """
        show_conf

        SYNOPSIS
           Show results from last time confidence was run

        SYNTAX

        Arguments:
           outfile  - filename to capture the output
                     default = None

           clobber - overwrite outfile if exists
                     default = False

        Returns:
           None

        DESCRIPTION
           Shows results from the last time confidence was run to
           determine parameter confidence limits.

           Examples:
              show_conf()

              show_conf("sherpa.conf", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_data, show_model,
           show_fit, show_covar
        """
        all = ''
        all += self._get_show_conf()
        _send_to_pager(all, outfile, clobber)

    def show_proj(self, outfile=None, clobber=False):
        """
        show_proj

        SYNOPSIS
           Show results from last time projection was run

        SYNTAX

        Arguments:
           outfile  - filename to capture the output
                     default = None

           clobber - overwrite outfile if exists
                     default = False

        Returns:
           None

        DESCRIPTION
           Shows results from the last time projection was run to
           determine parameter confidence limits.

           Examples:
              show_proj()

              show_proj("sherpa.proj", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_data, show_model,
           show_fit, show_covar
        """
        all = ''
        all += self._get_show_proj()
        _send_to_pager(all, outfile, clobber)

    def show_covar(self, outfile=None, clobber=False):
        """
        show_covar

        SYNOPSIS
           Show results from the last time covariance was run to
           determine parameter confidence limits

        SYNTAX

        Arguments:
           outfile  - filename to capture the output
                     default = None

           clobber - overwrite outfile if exists
                     default = False

        Returns:
           None

        DESCRIPTION
           Shows results from the last time covariance was run to
           determine parameter confidence limits

           Examples:
              show_covar()

              show_covar("sherpa.covar", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_all, show_data, show_model,
           show_fit, show_conf, show_proj
        """
        all = ''
        all += self._get_show_covar()
        _send_to_pager(all, outfile, clobber)


    def show_all(self, id=None, outfile=None, clobber=False):
        """
        show_all

        SYNOPSIS
           Show current state of Sherpa fitting session

        SYNTAX

        Arguments:
           filename   - name of saved file
                        default = None

        Returns:
           None

        DESCRIPTION
           Show current state of Sherpa fitting session including Opt Method,
           Statistic, and associated Data Sets and Models by Sherpa data id.
           If no data id is given then all available data ids will be used.

           Examples:
              show_all()

              show_all(1)

              show_all(1, "sherpa.session", True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, list_functions, show_data, show_model, show_fit,
           show_conf, show_proj, show_covar
        """
        all = ''
        all += self._get_show_data(id)
        all += self._get_show_model(id)
        all += self._get_show_fit()
        all += self._get_show_conf()
        all += self._get_show_proj()
        all += self._get_show_covar()
        _send_to_pager(all, outfile, clobber)

    def get_functions(self):
        """
        get_functions
        
        SYNOPSIS
           Return all available Sherpa functions in a list
        
        SYNTAX
        
        Arguments:
           None
        
        Returns:
           List of all Sherpa function names

        DESCRIPTION
           Returns a list containing names of all functions
           defined in the high-level Sherpa UI.

        SEE ALSO
            list_functions
        """
        funcs = []
        for func in dir(self):
            if not func.startswith('_') and callable(getattr(self,func)):
                funcs.append(func)
        return funcs
        
    def list_functions(self, outfile=None, clobber=False):
        """
        list_functions

        SYNOPSIS
           List all available Sherpa functions

        SYNTAX

        Arguments:
           outfile    - name of saved file
                        default = None

           clobber    - overwrite outfile if exists
                        default = False

        Returns:
           None

        DESCRIPTION
           List all available Sherpa functions from currrent Sherpa namespace.

              list_functions()

              list_functions("funcs.txt")

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           save, clean, show_all
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

    def get_default_id(self):
        """
        get_default_id

        SYNOPSIS
           Return the default Sherpa data id

        SYNTAX

        Arguments:
           None

        Returns:
           Sherpa data id

        DESCRIPTION
           The Sherpa data id ties data, model, fit, and plotting information
           into a dataset easily referenced by id.  The id can be a user
           defined string or integer.

        SEE ALSO
           set_default_id, list_data_ids
        """
        return self._default_id

    def set_default_id(self, id):
        """
        set_default_id

        SYNOPSIS
           Set the default Sherpa data id

        SYNTAX

        Arguments:
           id         - new Sherpa data id

        Returns:
           None

        DESCRIPTION
           The Sherpa data id ties data, model, fit, and plotting information
           into a dataset easily referenced by id.  The id can be a user
           defined string or integer.

        SEE ALSO
           get_default_id, list_data_ids
        """
        self._default_id = self._fix_id(id)


    ###########################################################################
    # Optimization methods
    ###########################################################################


    def list_methods(self):
        """
        list_methods

        SYNOPSIS
           List the available optimization methods in Sherpa

        SYNTAX

        Arguments:
           None

        Returns:
           list of methods

        DESCRIPTION

        SEE ALSO
           get_method, get_method_name, set_method, get_method_opt,
           set_method_opt
        """
        keys = self._methods.keys()[:]
        keys.sort()
        return keys

    def _get_method_by_name(self, name):
        meth = self._methods.get(name.lower())
        if meth is None:
            raise ArgumentErr('badmethod', name)
        return meth

    def get_method(self, name=None):
        """
        get_method

        SYNOPSIS
           Return the Sherpa optimization method by name

        SYNTAX

        Arguments:
           name       - name of opt method
                        default = default method

        Returns:
           Sherpa opt method

        DESCRIPTION

        SEE ALSO
           list_methods, get_method_name, set_method, get_method_opt,
           set_method_opt
        """
        if name is None:
            return self._current_method
        _check_type(name, basestring, 'name', 'a string')
        return self._get_method_by_name(name)

    def get_method_name(self):
        """
        get_method_name

        SYNOPSIS
           Return the name of current Sherpa optimization method

        SYNTAX

        Arguments:
           None

        Returns:
           name

        DESCRIPTION

        SEE ALSO
           list_methods, get_method, set_method, get_method_opt,
           set_method_opt
        """
        return type(self.get_method()).__name__.lower()

    def set_method(self, meth):
        """
        set_method

        SYNOPSIS
           Set the Sherpa optimization method by name

        SYNTAX

        Arguments:
           name       - name of opt method

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           list_methods, get_method, get_method_name, get_method_opt,
           set_method_opt
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

    def get_method_opt(self, optname=None):
        """
        get_method_opt

        SYNOPSIS
           Return a Sherpa optimization method option by name

        SYNTAX

        Arguments:
           name       - opt method option name

        Returns:
           opt method option value

        DESCRIPTION

        SEE ALSO
           list_methods, get_method, get_method_name, set_method,
           set_method_opt
        """
        if optname is None:
            return self._current_method.config
        self._check_method_opt(optname)
        return self._current_method.config[optname]

    def set_method_opt(self, optname, val):
        """
        set_method_opt

        SYNOPSIS
           Set a Sherpa optimization method option by name

        SYNTAX

        Arguments:
           name       - opt method option name
           val        - opt method option value

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           list_methods, get_method, get_method_name, set_method,
           get_method_opt
        """
        self._check_method_opt(optname)
        self._current_method.config[optname] = val

    #### Iterative Fitting Methods for CIAO 4.3 testing
    def get_iter_method_name(self):
        """
        get_iter_method_name

        SYNOPSIS
           Return the name of the current Sherpa iterative fitting method

        SYNTAX

        Arguments:

        Returns:
           Name of Sherpa iterative fitting method
           
        DESCRIPTION

        SEE ALSO
           list_iter_methods, set_iter_method, get_iter_method_opt,
           set_iter_method_opt
        """
        return self._current_itermethod['name']

    def get_iter_method_opt(self, optname=None):
        """
        get_iter_method_opt

        SYNOPSIS
           Return a Sherpa iterative fitting method option by name

        SYNTAX

        Arguments:
           name       - option name

        Returns:
           iterative fitting method option value
        
        DESCRIPTION

        SEE ALSO
           list_iter_methods, get_iter_method_name, set_iter_method,
           set_iter_method_opt
        """
        itermethod_opts = dict(self._current_itermethod)
        del itermethod_opts['name']
        
        if optname is None:
            return itermethod_opts

        _check_type(optname, basestring, 'optname', 'a string')
        if optname not in itermethod_opts:
            raise ArgumentErr('badopt', optname, self._current_itermethod['name'])
        return itermethod_opts[optname]

    def list_iter_methods(self):
        """
        list_iter_methods

        SYNOPSIS
           List the available iterative fitting methods in Sherpa

        SYNTAX

        Arguments:
           None

        Returns:
           list of iterative fitting methods

        DESCRIPTION

        SEE ALSO
           get_iter_method_name, set_iter_method, get_iter_method_opt,
           set_iter_method_opt
        """
        keys = self._itermethods.keys()[:]
        keys.sort()
        return keys

    def set_iter_method(self, meth):
        """
        set_iter_method

        SYNOPSIS
           Set the Sherpa iterative fitting method by name

        SYNTAX

        Arguments:
           name       - name of iterative method ['none'|'primini'|'sigmarej']

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           list_iter_methods, get_iter_method_name, get_iter_method_opt,
           set_iter_method_opt
        """
        if isinstance(meth, basestring):
            if (self._itermethods.has_key(meth) == True):
                self._current_itermethod = self._itermethods[meth]
            else:
                raise TypeError(meth + ' is not an iterative fitting method')
        else:
            _argument_type_error(meth, 'a string')

    def set_iter_method_opt(self, optname, val):
        """
        set_iter_method_opt

        SYNOPSIS
           Set a Sherpa iterative fitting method option by name

        SYNTAX

        Arguments:
           name       - option name
           val        - option value

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           list_methods, get_method, get_method_name, set_method,
           get_method_opt
        """
        _check_type(optname, basestring, 'optname', 'a string')
        if (optname not in self._current_itermethod or
            optname == 'name'):
            raise ArgumentErr('badopt', optname, self._current_itermethod['name'])
        self._current_itermethod[optname] = val

    ###########################################################################
    # Statistics
    ###########################################################################


    def list_stats(self):
        """
        list_stats

        SYNOPSIS
           List statistics available in Sherpa

        SYNTAX

        Arguments:
           None

        Returns:
           list of available stats

        DESCRIPTION
           Available statistics include

           * 'chi2constvar'  \chi^2 with constant variance computed
                             from the counts data.

           * 'chi2modvar'    \chi^2 with model amplitude variance.

           * 'chi2gehrels'   \chi^2 with gehrels method (Sherpa default).

           * 'chi2datavar'   \chi^2 with data variance.

           * 'chi2xspecvar'  \chi^2 with data variance (XSPEC-style,
                             variance = 1.0 if data less than or equal to 0.0).

           * 'cstat'         CStat - A maximum likelihood function
                             (XSPEC implementation of Cash).

           * 'cash'          Cash  - A maximum likelihood function.

        SEE ALSO
           get_stat, get_stat_name, set_stat
        """
        keys = self._stats.keys()[:]
        keys.sort()
        return keys

    def _get_stat_by_name(self, name):
        stat = self._stats.get(name.lower())
        if stat is None:
            raise ArgumentErr('badstat', name)
        return stat

    def get_stat(self, name=None):
        """
        get_stat

        SYNOPSIS
           Return a Sherpa statistic by name

        SYNTAX

        Arguments:
           name       - name of statistic

        Returns:
           statistic

        DESCRIPTION
           Available statistics include

           * 'chi2constvar'  \chi^2 with constant variance computed
                             from the counts data.

           * 'chi2modvar'    \chi^2 with model amplitude variance.

           * 'chi2gehrels'   \chi^2 with gehrels method (Sherpa default).

           * 'chi2datavar'   \chi^2 with data variance.

           * 'chi2xspecvar'  \chi^2 with data variance (XSPEC-style,
                             variance = 1.0 if data less than or equal to 0.0).

           * 'cstat'         CStat - A maximum likelihood function
                             (XSPEC implementation of Cash).

           * 'cash'          Cash  - A maximum likelihood function.

        SEE ALSO
           list_stats, get_stat_name, set_stat
        """
        if name is None:
            return self._current_stat
        _check_type(name, basestring, 'name', 'a string')
        return self._get_stat_by_name(name)

    def get_stat_name(self):
        """
        get_stat_name

        SYNOPSIS
           Return the current Sherpa statistic by name

        SYNTAX

        Arguments:
           None

        Returns:
           statistic name

        DESCRIPTION
           Available statistics include

           * 'chi2constvar'  \chi^2 with constant variance computed
                             from the counts data.

           * 'chi2modvar'    \chi^2 with model amplitude variance.

           * 'chi2gehrels'   \chi^2 with gehrels method (Sherpa default).

           * 'chi2datavar'   \chi^2 with data variance.

           * 'chi2xspecvar'  \chi^2 with data variance (XSPEC-style,
                             variance = 1.0 if data less than or equal to 0.0).

           * 'cstat'         CStat - A maximum likelihood function
                             (XSPEC implementation of Cash).

           * 'cash'          Cash  - A maximum likelihood function.

        SEE ALSO
           list_stats, get_stat, set_stat
        """
        return type(self.get_stat()).__name__.lower()

    def set_stat(self, stat):
        """
        set_stat

        SYNOPSIS
           Set the current Sherpa statistic by name

        SYNTAX

        Arguments:
           name       - name of statistic

        Returns:
           None

        DESCRIPTION
           Available statistics include

           * 'chi2constvar'  \chi^2 with constant variance computed
                             from the counts data.

           * 'chi2modvar'    \chi^2 with model amplitude variance.

           * 'chi2gehrels'   \chi^2 with gehrels method (Sherpa default).

           * 'chi2datavar'   \chi^2 with data variance.

           * 'chi2xspecvar'  \chi^2 with data variance (XSPEC-style,
                             variance = 1.0 if data less than or equal to 0.0).

           * 'cstat'         CStat - A maximum likelihood function
                             (XSPEC implementation of Cash).

           * 'cash'          Cash  - A maximum likelihood function.

        SEE ALSO
           list_stats, get_stat, get_stat_name
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


    def list_data_ids(self):
        """
        list_data_ids

        SYNOPSIS
           List the available Sherpa data ids

        SYNTAX

        Arguments:
           None
          
        Returns:
           list of data ids

        DESCRIPTION
           The Sherpa data id ties data, model, fit, and plotting information
           into a dataset easily referenced by id.  The id can be a user
           defined string or integer.

        SEE ALSO
           set_default_id, get_default_id
        """
        keys = self._data.keys()[:]
        keys.sort()
        return keys

    def get_data(self, id=None):
        """
        get_data

        SYNOPSIS
           Return a Sherpa dataset by data id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

        Returns:
           dataset

        DESCRIPTION
           The Sherpa data id ties data, model, fit, and plotting information
           into a dataset easily referenced by id.  The id can be a user
           defined string or integer.

        SEE ALSO
           list_data_ids, set_data, copy_data, delete_data,
           read_data, load_data
        """
        return self._get_item(id, self._data, 'data set', 'has not been set')

    def set_data(self, id, data=None):
        """
        set_data

        SYNOPSIS
           Set a dataset by Sherpa data id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           data       - data object

        Returns:
           None

        DESCRIPTION
           The Sherpa data id ties data, model, fit, and plotting information
           into a dataset easily referenced by id.  The id can be a user
           defined string or integer.

        SEE ALSO
           list_data_ids, get_data, copy_data, delete_data,
           read_data, load_data
        """
        if data is None:
            id, data = data, id
        self._set_item(id, data, self._data, sherpa.data.Data, 'data',
                       'a data set')


    def _read_error(self, filename, *args, **kwargs):
        err = sherpa.io.get_ascii_data(filename, *args, **kwargs)[1].pop()
        return err


    def load_staterror(self, id, filename=None, ncols=2, *args, **kwargs):
        """
        load_staterror

        SYNOPSIS
           Load the statistical errors for a dataset from file

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filename   - filename with path

           ncols      - number of columns to read from
                        default = 2

           colkeys    - column keys
                        default = None

           sep        - separator character
                        default = ' '

           comment    - comment character
                        default = '#'

        Returns:
           None

        DESCRIPTION
           Load the statistical error for a dataset from file by data id.  
           Users can specify the column name by using the colkeys argument to 
           set the statistical error.

        EXAMPLE
           load_staterror("data.dat", colkeys=["STAT_ERR"])

        SEE ALSO
           load_syserror, set_staterror, set_syserror
        """
        if filename is None:
            id, filename = filename, id

        self.set_staterror(id, self._read_error(filename, ncols=ncols,
                                                *args, **kwargs))


    def load_syserror(self, id, filename=None, ncols=2, *args, **kwargs):
        """
        load_syserror

        SYNOPSIS
           Load the systematic errors for a dataset from file

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filename   - filename with path

           ncols      - number of columns to read from
                        default = 2

           colkeys    - column keys
                        default = None

           sep        - separator character
                        default = ' '

           comment    - comment character
                        default = '#'

        Returns:
           None

        DESCRIPTION
           Load the systematic error for a dataset from file by data id.  Users 
           can specify the column name by using the colkeys argument to set the
           systematic error.

        EXAMPLE
           load_syserror("data.dat", colkeys=["SYS_ERR"])

        SEE ALSO
           load_staterror, set_staterror, set_syserror
        """
        if filename is None:
            id, filename = filename, id

        self.set_syserror(id, self._read_error(filename, ncols=ncols,
                                               *args, **kwargs))


    def load_filter(self, id, filename=None, ignore=False, ncols=2,
                    *args, **kwargs):
        """
        load_filter

        SYNOPSIS
           Load the dataset filter from file

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filename   - filename with path

           ignore     - non-zero values ignore instead of notice
                        default = False

           ncols      - number of columns to read from
                        default = 2

           colkeys    - column keys
                        default = None

           sep        - separator character
                        default = ' '

           comment    - comment character
                        default = '#'

        Returns:
           None

        DESCRIPTION
           Load the filter for a dataset from file by data id.

        EXAMPLE
           load_filter("data.dat", colkeys=["FILTER"])

        SEE ALSO
            set_filter
        """
        if filename is None:
            id, filename = filename, id

        self.set_filter(id,
            self._read_error(filename, ncols=ncols, *args, **kwargs),
                        ignore=ignore)


    def set_filter(self, id, val=None, ignore=False):
        """
        set_filter

        SYNOPSIS
           Set the dataset filter by data id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           val        - array of 0s or 1s

           ignore     - non-zero values ignore instead of notice
                        default = False

        Returns:
           None

        DESCRIPTION
           Set the filter of a dataset by data id.  

        EXAMPLE
           set_filter([0, 1, 1, ...])

        SEE ALSO
           load_filter
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


    def set_dep(self, id, val=None):
        """
        set_dep

        SYNOPSIS
           Set the dependent variable of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           val        - dependent variable array or scalar

        Returns:
           None

        DESCRIPTION
           Set the dependent variable of a data set by data id.

        EXAMPLE
           set_dep([1,2,3,...])

           set_dep(1,1)

        SEE ALSO
           get_dep, get_indep
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


    def set_staterror(self, id, val=None, fractional=False):
        """
        set_staterror

        SYNOPSIS
           Set the statistical errors of a dataset by data id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           val        - array or scalar error values

           fractional - use fractional portion of dependent array as error,
                        val must be a scalar value
                        default = False

        Returns:
           None

        DESCRIPTION
           Set the statistical error of a dataset by data id.  Users can specify
           the entire error as an array or as a single value to be repeated for
           every bin.  Also, setting the fractional argument will use the single
           value as the fractional portion of the dependent array as the error.

        EXAMPLE
           set_staterror([0.040, 0.039, 0.041, ...])

           set_staterror(2, 0.04)

           set_staterror(0.05, fractional=True)

        SEE ALSO
           set_syserror
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


    def set_syserror(self, id, val=None, fractional=False):
        """
        set_syserror

        SYNOPSIS
           Set the systematic errors of a dataset by data id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           val        - array or scalar error values

           fractional - use fractional portion of dependent array as error,
                        val must be a scalar value
                        default = False

        Returns:
           None

        DESCRIPTION
           Set the systematic error of a dataset by data id.  Users can specify
           the entire error as an array or as a single value to be repeated for
           every bin.  Also, setting the fractional argument will use the single
           value as the fractional portion of the dependent array as the error.

        EXAMPLE
           set_syserror([0.040, 0.039, 0.041, ...])

           set_syserror(2, 0.04)

           set_syserror(0.05, fractional=True)

        SEE ALSO
           set_syserror
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


    def get_staterror(self, id=None, filter=False):
        """
        get_staterror

        SYNOPSIS
           Get the statistical errors of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filter     - apply filter
                        default = False

        Returns:
           Statistical error array

        DESCRIPTION
           Get the statistical error of a dataset by data id.

        EXAMPLE
           get_staterror()

           get_staterror(1)

           get_staterror(1, True)

        SEE ALSO
           set_syserror, get_syserror
        """
        return self.get_data(id).get_staterror(filter,
                                               self.get_stat().calc_staterror)


    def get_syserror(self, id=None, filter=False):
        """
        get_syserror

        SYNOPSIS
           Get the systematic errors of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filter     - apply filter
                        default = False

        Returns:
           Systematic error array

        DESCRIPTION
           Get the systematic error of a dataset by data id.

        EXAMPLE
           get_syserror()

           get_syserror(1)

           get_syserror(1, True)

        SEE ALSO
           set_syserror, get_syserror
        """
        d = self.get_data(id)
        id = self._fix_id(id)
        err = d.get_syserror(filter)
        if err is None or not numpy.iterable(err):
            raise sherpa.utils.err.DataErr('nosyserr', id)
        return err


    def get_error(self, id=None, filter=False):
        """
        get_error

        SYNOPSIS
           Get the total errors of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filter     - apply filter
                        default = False

        Returns:
           Total error array

        DESCRIPTION
           Get the total error (statistical + systematic in quadrature) of a
           dataset by data id.

        EXAMPLE
           get_error()

           get_error(1)

           get_error(1, True)

        SEE ALSO
           set_syserror, get_syserror, set_staterror, set_syserror
        """
        return self.get_data(id).get_error(filter,
                                           self.get_stat().calc_staterror)

    def get_indep(self, id=None):
        """
        get_indep

        SYNOPSIS
           Get the independent grid of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

        Returns:
           Array of the independent variable

        DESCRIPTION
           Get the data set independend grid by data id.

        EXAMPLE
           get_indep()

           get_indep(1)

        SEE ALSO

        """
        return self.get_data(id).get_indep()


    def get_dep(self, id=None, filter=False):
        """
        get_dep

        SYNOPSIS
           Get the dependent variable of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filter     - apply filter
                        default = False

        Returns:
           Array of the dependent variable

        DESCRIPTION
           Get the dependent variable array of data set by data id.

        EXAMPLE
           get_dep()

           get_dep(1)

        SEE ALSO

        """
        return self.get_data(id).get_y(filter)


    def get_dims(self, id=None, filter=False):
        """
        get_dims

        SYNOPSIS
           Get the dimensionality of a dataset by id

        SYNTAX

        Arguments:
           id        - session data id
                       default = default data id

           filter    - apply filter
                       default = False

        Returns:
           List of dimensional lengths

        DESCRIPTION
           Get the dimensionality (dim0, dim1, ...) of data set by data id.

        EXAMPLE
           get_dims()

           get_dims(1, True)

        SEE ALSO

        """
        return self.get_data(id).get_dims(filter)


    def get_filter(self, id=None):
        """
        get_filter

        SYNOPSIS
           Get the filter of a dataset by id

        SYNTAX

        Arguments:
           id        - session data id
                       default = default data id

        Returns:
           filter string

        DESCRIPTION
           Get the filter expression of data set by data id.

        EXAMPLE
           get_filter()

           get_filter(1)

        SEE ALSO

        """
        return self.get_data(id).get_filter()


    def copy_data(self, fromid, toid):
        """
        copy_data

        SYNOPSIS
           Copy a dataset by data id to a new data id (deep copy)

        SYNTAX

        Arguments:
           fromid     - source data id
           toid       - destination data id

        Returns:
           None

        DESCRIPTION
           The Sherpa data id ties data, model, fit, and plotting information
           into a dataset easily referenced by id.  The id can be a user
           defined string or integer.

        SEE ALSO
           list_data_ids, get_data, set_data, delete_data,
           read_data, load_data
        """
        data = self.get_data(fromid)
        data = copy.deepcopy(data)
        self.set_data(toid, data)

    def delete_data(self, id=None):
        """
        delete_data

        SYNOPSIS
           Delete a dataset by data id

        SYNTAX

        Arguments:
           id         - Sherpa data id
                        default = default data id

        Returns:
           None

        DESCRIPTION
           The Sherpa data id ties data, model, fit, and plotting information
           into a dataset easily referenced by id.  The id can be a user
           defined string or integer.

        SEE ALSO
           list_data_ids, get_data, set_data, copy_data,
           read_data, load_data
        """
        id = self._fix_id(id)
        self._data.pop(id, None)


    def dataspace1d(self, start, stop, step=1, numbins=None, 
                    id=None, dstype=sherpa.data.Data1DInt):
        """
        dataspace1d

        SYNOPSIS
           Populates a blank 1D Sherpa data set by data id

        SYNTAX

        Arguments:
           start   -  lower bound of grid

           stop    -  upper bound of grid

           step    -  bin width size
                      default is 1

           numbins -  number of bins desired
                      default is None

           id      -  Sherpa data id
                      defaut is default data id

           dstype  -  Type of data set to use
                      default is Data1DInt

        Returns:
           None

        DESCRIPTION
           Populates a blank 1D Sherpa data set with the specified grid
           by Sherpa data id.
           
           Specifying a dataspace using step size:
           if numbins is None (default) -> numpy.arange(start,stop,step)

           Specifying a dataspace by indicating the number of bins:
           if numbins is not None -> numpy.linspace(start, stop, numbins)

        EXAMPLES
           Blank integrated data set
           
              dataspace1d(0.1,10,0.1)

           Blank non-integrated data set

              dataspace1d(0.1,10,0.1,dstype=Data1D)

        SEE ALSO
           dataspace2d
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


    def dataspace2d(self, dims, id=None, dstype=sherpa.data.Data2D):
        """
        dataspace2d

        SYNOPSIS
           Populates a blank 2D Sherpa image data set by data id

        SYNTAX

        Arguments:
           dims    -  array of image dimensions, i.e. [width,height]

           id      -  Sherpa data id
                      defaut is default data id

           dstype  -  Type of data set to use
                      default is Data2D

        Returns:
           None

        DESCRIPTION
           Populates a blank 2D Sherpa image data set by Sherpa data id.

        SEE ALSO
           dataspace1d
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

    def unpack_arrays(self, *args):
        """
        unpack_arrays
        
        SYNOPSIS
           Read NumPy arrays into a dataset

        SYNTAX

        Arguments:
           array0     - first NumPy array

           ...

           arrayN     - last NumPy array

           dstype     - dataset type desired
                        default = Data1D

        Returns:
           Sherpa dataset

        DESCRIPTION
           Read NumPy arrays into a Sherpa dataset.

        SEE ALSO
           unpack_data, load_data, load_arrays
        """
        return sherpa.io.read_arrays(*args)


    def unpack_data(self, filename, ncols=2, colkeys=None,
                    dstype=sherpa.data.Data1D, sep=' ', comment='#', require_floats=True):
        """
        unpack_data

        SYNOPSIS
           Read a data text file into a dataset

        SYNTAX

        Arguments:
           filename   - name of text file

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
           Sherpa dataset

        DESCRIPTION
           Read tabular data from column based text file into a
           Sherpa dataset given a filename.

        SEE ALSO
           load_data, unpack_arrays, load_arrays
        """
        return self._read_data(sherpa.io.read_data, filename, ncols, colkeys,
                               sep, dstype, comment, require_floats)


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

    ##@loggable(with_id = True)
    def load_arrays(self, id, *args):
        """
        load_arrays
        
        SYNOPSIS
           Load NumPy arrays into a dataset

        SYNTAX

        Arguments:
           id         - data id
        
           array0     - first NumPy array

           ...

           arrayN     - last NumPy array

           dstype     - dataset type desired
                        default = Data1D

        Returns:
           None

        DESCRIPTION
           Load NumPy arrays into a Sherpa dataset by data id.

        SEE ALSO
           load_data, unpack_arrays
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


    def save_arrays(self, filename, args, fields=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'):
        """
        save_arrays

        SYNOPSIS
           Write a list of arrays to file as columns

        SYNTAX

        Arguments:
           filename   - filename with path

           args       - list of arrays that correspond to columns

           fields     - list of column headings
                        default = None

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
           Write a list of arrays to file as columns.

        EXAMPLE

           save_arrays("foo.dat", [a,b,c], fields=['a','b','c'])

        SEE ALSO
           save_image, save_data, save_table, save_source, save_model,
           save_resid, save_delchi
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

    def save_resid(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'):
        """
        save_resid

        SYNOPSIS
           Write the data - model residuals to file

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
           Write the data - model residuals to file. NOTE that the residuals
           array written to file respects the filter.

        EXAMPLE

           save_resid("resid.dat")

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_delchi
        """
        clobber=sherpa.utils.bool_cast(clobber)
        _check_type(filename, basestring, 'filename', 'a string')
        self._save_type('resid', id, filename, clobber=clobber, sep=sep,
                        comment=comment, linebreak=linebreak, format=format)

    def save_delchi(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'):
        """
        save_delchi

        SYNOPSIS
           Write the delta chi squared residuals to file

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
           Write the delta chi squared residuals to file.  NOTE that the 
           delta chi squared residuals array written to file respects the
           filter.

        EXAMPLE

           save_delchi("delchi.dat")

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_resid
        """
        clobber=sherpa.utils.bool_cast(clobber)
        _check_type(filename, basestring, 'filename', 'a string')
        self._save_type('delchi', id, filename, clobber=clobber, sep=sep,
                        comment=comment, linebreak=linebreak, format=format)


    def save_data(self, id, filename=None, fields=None, sep=' ', comment='#',
                  clobber=False, linebreak='\n', format='%g'):
        """
        save_data

        SYNOPSIS
           Write tabular data by id

        SYNTAX

        Arguments:
           id         - dataset ID
                        default = default data id
           filename   - filename with path

        Keyword Arguments:

           fields     - dataset attribute names
                        default = None

           sep        - separation character
                        default = ' '

           comment    - comment character
                        default = '#'

           clobber    - clobber output file
                        default = False

           linebreak  - new line character
                        default = '\n'

           format     - format strings for array element
                        default = '%g'

        Returns:
           None

        DESCRIPTION
           Write tabular data to a column-based text file from a
           Sherpa dataset by id.

        SEE ALSO
           save_pha, save_arf, save_rmf, save_table, save_image
        """
        clobber=sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        sherpa.io.write_data(filename, self.get_data(id), fields, sep,
                             comment, clobber, linebreak, format)


    def save_filter(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'):
        """
        save_filter

        SYNOPSIS
           Write the data set filter to file

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
           Write the data set filter to file.

        EXAMPLE

           save_filter("filter.dat")

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_delchi
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


    def save_staterror(self, id, filename=None, clobber=False, sep=' ',
                       comment='#', linebreak='\n', format='%g'):
        """
        save_staterror

        SYNOPSIS
           Write the data set statistical errors to file

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
           Write the data set statistical errors to file by id.

        EXAMPLE

           save_staterror("staterror.dat")

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_delchi
        """
        clobber=sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        x = self.get_data(id).get_indep(filter=False)[0]
        err = self.get_staterror(id, filter=False)
        self.save_arrays(filename, [x, err], ['X', 'STAT_ERR'],
                         clobber, sep, comment, linebreak, format)


    def save_syserror(self, id, filename=None, clobber=False, sep=' ',
                       comment='#', linebreak='\n', format='%g'):
        """
        save_syserror

        SYNOPSIS
           Write the data set systematic errors to file

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
           Write the data set systematic errors to file by id.

        EXAMPLE

           save_syserror("syserror.dat")

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_delchi
        """
        clobber=sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        x = self.get_data(id).get_indep(filter=False)[0]
        err = self.get_syserror(id, filter=False)
        self.save_arrays(filename, [x, err], ['X', 'SYS_ERR'],
                         clobber, sep, comment, linebreak, format)

    def save_error(self, id, filename=None, clobber=False, sep=' ',
                       comment='#', linebreak='\n', format='%g'):
        """
        save_error

        SYNOPSIS
           Write the total errors of a data set to file

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
           Write the total errors (statistical + systematic in quadrature) of a
           data set to file by id.

        EXAMPLE

           save_error("error.dat")

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_delchi, save_syserror, save_staterror
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


    def notice(self, lo=None, hi=None, **kwargs):
        """
        notice

        SYNOPSIS
           Exclusive 1D notice on interval(s) for all available
           Sherpa data ids

        SYNTAX

        Arguments:

           lo    -   lower bound OR interval expression string
                     default = None

           hi    -   upper bound
                     default = None

        Returns:
           None

        DESCRIPTION

           notice()

           notice(0.5, 7.0)

           notice("0.5:7.0, 8.0:10.0")

           notice(":0.5, 7.0:")

        SEE ALSO
           notice_id, ignore, ignore_id
        """
        if lo is not None and type(lo) in (str,numpy.string_):
            return self._notice_expr(lo, **kwargs)
        if self._data.items() == []:
            raise IdentifierErr("nodatasets")
        for d in self._data.values():
            d.notice(lo, hi, **kwargs)


    def ignore(self, lo=None, hi=None, **kwargs):
        """
        ignore

        SYNOPSIS
           Exclusive 1D ignore on interval(s) for all available 
           Sherpa data ids

        SYNTAX

        Arguments:

           lo    -   lower bound OR interval expression string
                     default = None

           hi    -   upper bound
                     default = None

        Returns:
           None

        DESCRIPTION

           ignore()

           ignore(0.5, 7.0)

           ignore(":0.5, 7.0:")

           ignore(":0.5, 7.0:")

        SEE ALSO
           notice_id, notice, ignore_id
        """
        kwargs['ignore'] = True
        if lo is not None and type(lo) in (str,numpy.string_):
            return self._notice_expr(lo, **kwargs)
        self.notice(lo, hi, **kwargs)

    def notice_id(self, ids, lo=None, hi=None, **kwargs):
        """
        notice_id

        SYNOPSIS
           Exclusive 1D notice on interval(s) for specific Sherpa data id(s)

        SYNTAX

        Arguments:

           ids   -  list of specific Sherpa data ids

           lo    -  lower bound OR interval expression string
                    default = None

           hi    -  upper bound
                    default = None

        Returns:
           None

        DESCRIPTION

           notice_id(1)

           notice_id(1, 0.5, 7.0)

           notice_id(2, "0.5:7.0, 8.0:10.0")

           notice_id([2,3], ":0.5, 7.0:")

        SEE ALSO
           notice, ignore, ignore_id
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


    def ignore_id(self, ids, lo=None, hi=None, **kwargs):
        """
        ignore_id

        SYNOPSIS
           Exclusive 1D ignore on interval(s) for specific Sherpa data id(s)

        SYNTAX

        Arguments:

           ids   -  list of specific Sherpa data ids

           lo    -  lower bound OR interval expression string
                    default = None

           hi    -  upper bound
                    default = None

        Returns:
           None

        DESCRIPTION

           ignore_id(1)

           ignore_id(1, 0.5, 7.0)

           ignore_id(2, ":0.5, 7.0:")

           ignore_id([2,3], ":0.5, 7.0:")

        SEE ALSO
           notice, ignore, notice_id
        """
        kwargs['ignore'] = True
        if lo is not None and type(lo) in (str,numpy.string_):
            return self._notice_expr_id(ids, lo, **kwargs)
        self.notice_id(ids, lo, hi, **kwargs)


    ###########################################################################
    # Models
    ###########################################################################

    def paramprompt(self, val=False):
        """
        paramprompt

        SYNOPSIS
           Prompts the user for initial, minimum, and maximum parameter values

        SYNTAX

        Arguments:

           val    - boolean value indicating prompt behavior
                    default = False

        Returns:
           None

        DESCRIPTION
           Prompts for the initial, minimum, and maximum parameter values.  User provides
           the values in a respectively delimited with commas.  The minimum and maximum
           values are optional.

        EXAMPLE
            abs1.nh parameter value 0.07

            abs1.nh parameter value  0.07, 1, 10

        SEE ALSO
           set_model
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

    def get_model_autoassign_func(self):
        """
        get_model_autoassign_func

        SYNOPSIS
           Return a function pointer to the model assignment function

        SYNTAX

        Arguments:
           None

        Returns:
           function ptr

        DESCRIPTION

        SEE ALSO
           set_model_autoassign_func
        """
        return self._model_autoassign_func

    def set_model_autoassign_func(self, func=None):
        """
        set_model_autoassign_func

        SYNOPSIS
           Sets the model assignment function to a user defined function
           pointer

        SYNTAX

        Arguments:
           func       - model assignment function pointer
                        default = None

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           get_model_autoassign_func
        """
        if (func is not None) and (not callable(func)):
            _argument_type_error('func', 'a function or other callable object')
        self._model_autoassign_func = func

    def list_models(self, show="all"):
        """
        list_models

        SYNOPSIS
           List the available Sherpa model types

        SYNTAX

        Arguments:
           show  -  filter list by keywords "all", "xspec", "1d", and "2d"
                    default = "all"

        Returns:
           list of available model types

        DESCRIPTION

        SEE ALSO
           list_model_components, create_model_component,
           delete_model_component
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

    def list_model_components(self):
        """
        list_model_components

        SYNOPSIS
           List the Sherpa model components of active models

        SYNTAX

        Arguments:
           None

        Returns:
           list of available model components

        DESCRIPTION

        SEE ALSO
           list_models, create_model_component,
           delete_model_component
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


    def get_model_component(self, name):
        """
        get_model_component

        SYNOPSIS
           Access a Sherpa model component by name

        SYNTAX

        Arguments:
           name       - component label as a string

        Returns:
           Sherpa model component

        DESCRIPTION

        SEE ALSO
           list_models, list_model_components,
           delete_model_component, create_model_component
        """   
        # If user mistakenly passes an actual model reference,
        # just return the reference
        if isinstance(name, sherpa.models.Model):
            return name

        _check_type(name, basestring, 'name', 'a string')
        return self._get_model_component(name, require=True)


    def create_model_component(self, typename=None, name=None):
        """
        create_model_component

        SYNOPSIS
           Create a new Sherpa model component

        SYNTAX

        Arguments:
           typename   - name of component type
           name       - component label

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           list_models, list_model_components,
           delete_model_component
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

    def reset(self, model=None, id=None):
        """
        reset

        SYNOPSIS
           Resets model parameter values to defaults

        SYNTAX

        Arguments:
           name       - model instance

        Returns:
           None

        DESCRIPTION
           Resets model parameter values to defaults or user defined defaults.

           Example 1:
               reset( get_model() )

           Example 2:
               powlaw1d.p1
               beta1d.b1
               reset( p1*b1 )

           Example 3:
               beta1d.b1
               reset( b1 )

        SEE ALSO
           list_models, list_model_components,
           delete_model_component
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


    def delete_model_component(self, name):
        """
        delete_model_component

        SYNOPSIS
           Delete a Sherpa model component from active models

        SYNTAX

        Arguments:
           name       - component label

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           list_models, list_model_components,
           create_model_component        
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

    def list_model_ids(self):
        """
        list_model_ids

        SYNOPSIS
           List the Sherpa session model ids

        SYNTAX

        Arguments:
           None

        Returns:
           list of model ids

        DESCRIPTION
           List all the current active Sherpa session model ids.

        SEE ALSO
           get_model, set_model, delete_model, get_model_type, get_model_pars
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

    def get_source(self, id=None):
        """
        get_source

        SYNOPSIS
           Return a Sherpa model by model id

        SYNTAX

        Arguments:
           id         - model id
                        default = default model id

        Returns:
           Sherpa model

        DESCRIPTION
           Retrieve a Sherpa model by model id

        SEE ALSO
           list_model_ids, set_model, delete_model, get_model_type,
           get_model_pars
        """
        return self._get_source(id)

    def get_model(self, id=None):
        """
        get_model

        SYNOPSIS
           Return a Sherpa model by model id

        SYNTAX

        Arguments:
           id         - model id
                        default = default model id

        Returns:
           Sherpa model

        DESCRIPTION
           Retrieve the full convolved Sherpa model by model id

        SEE ALSO
           list_model_ids, set_model, delete_model, get_model_type,
           get_model_pars, get_model
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

    ##@loggable(with_id=True, with_keyword="model")
    def set_model(self, id, model=None):
        """
        set_model

        SYNOPSIS
           Set a Sherpa source model by model id

        SYNTAX

        Arguments:
           id         - model id
                        default = default model id

           model      - Sherpa model

        Returns:
           None

        DESCRIPTION
           Add a Sherpa source model to the list of current active 
           Sherpa source models by model id.

        SEE ALSO
           list_model_ids, get_model, delete_model, get_model_type,
           get_model_pars        
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


    def delete_model(self, id=None):
        """
        delete_model

        SYNOPSIS
           Delete a Sherpa model by model id

        SYNTAX

        Arguments:
           id         - id of model
                        default = default model id

        Returns:
           None

        DESCRIPTION
           Delete a Sherpa model from the list of currently active models by
           model id.

        SEE ALSO
           list_model_ids, set_model, get_model, get_model_type, get_model_pars
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

    def get_model_type(self, model):
        """
        get_model_type

        SYNOPSIS
           Return a Sherpa model type by model or model expression string

        SYNTAX

        Arguments:
           model      - model variable

        Returns:
           type of model

        DESCRIPTION
           Get the Sherpa model type by the model variable or model expression
           string.

           Example 1

               foo = gauss1d.foo
               get_model_type( foo )
           'gauss1d'
           
           Example 2
           
               get_model_type( gauss1d.foo * const1d.bar )
           'binaryopmodel'

           Example 3 ( from astro package )

               arf = get_arf()
               rmf = get_rmf()
               src = xsphabs.abs1 + powlaw1d.p1
               foo = (src.apply( arf.apply_arf )).apply( rmf.apply_rmf )
               get_model_type(foo)
           'nestedmodel'

        SEE ALSO
           list_model_ids, set_model, get_model, delete_model, get_model_pars
        """
        model = self._check_model(model)
        return type(model).__name__.lower()

    def get_model_pars(self, model):
        """
        get_model_pars

        SYNOPSIS
           Return a list of Sherpa model parameters by model or model 
           expression

        SYNTAX

        Arguments:
           model      - label of model object

        Returns:
           list of model parameters

        DESCRIPTION
           Get a list of Sherpa model parameters by model variable or model
           expression string.

           Example 1

               get_model_pars( gauss1d.foo.apply( psf1d.pp ) )
           ['fwhm', 'pos', 'ampl', 'xsize', 'xoff']

        SEE ALSO
           list_model_ids, set_model, get_model, delete_model, get_model_type
        """
        model = self._check_model(model)
        return [p.name for p in model.pars]

    def get_num_par(self, id=None):
        """
        get_num_par

        SYNOPSIS
           Return the number of parameters in a Sherpa model

        SYNTAX

        Arguments:
           id         - id of model
                        default = default model id

        Returns:
           Number of model parameters

        DESCRIPTION
           Returns the number of parameters in the model regardless
           of combination
        
        SEE ALSO
           get_num_par_thawed, get_num_par_frozen
        """
        return len(self._get_source(id).pars)
    
    def get_num_par_thawed(self, id=None):
        """
        get_num_par_thawed

        SYNOPSIS
           Return the number of thawed parameters in a Sherpa model

        SYNTAX

        Arguments:
           id         - id of model
                        default = default model id

        Returns:
           Number of thawed model parameters

        DESCRIPTION
           Returns the number of thawed parameters in the model regardless
           of combination
        
        SEE ALSO
           get_num_par, get_num_par_frozen
        """
        return len(self._get_source(id).thawedpars)

    def get_num_par_frozen(self, id=None):
        """
        get_num_par_frozen

        SYNOPSIS
           Return the number of frozen parameters in a Sherpa model

        SYNTAX

        Arguments:
           id         - id of model
                        default = default model id

        Returns:
           Number of frozen model parameters

        DESCRIPTION
           Returns the number of frozen parameters in the model regardless
           of combination
        
        SEE ALSO
           get_num_par, get_num_par_thawed
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
    # PSF
    #

    
    def load_psf(self, modelname, filename_or_model, *args, **kwargs):
        """
        load_psf

        SYNOPSIS
           load a file-based or model-based kernel into a PSF model

        SYNTAX

        Arguments:
           modelname - name of PSF

           filename_or_model - filename with path for file-based kernel
                               a Sherpa model for a model-based kernel

           args      - additional arguments when reading a file kernel

           kwargs    - additional keyword arguments when reading a file
                       kernel

        Returns:
           None

        DESCRIPTION
           Create a PSF model object with identifier 'modelname' and 
           initializes the PSF kernel to be either a Sherpa dataset
           loaded from file or a Sherpa model.

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

        psf = sherpa.instrument.PSFModel(modelname, kernel)
        self._add_model_component(psf)
        self._psf_models.append(psf)

    ##@loggable(with_id=True, with_keyword='psf')
    def set_psf(self, id, psf=None):
        """
        set_psf

        SYNOPSIS
           Add a PSF model by Sherpa data id

        SYNTAX

        Arguments:
           id   - Sherpa data id
                  default = default data id

           psf  - Sherpa PSF model

        Returns:
           None

        DESCRIPTION
           Add a PSF model to the instrument list by Sherpa data id.

        SEE ALSO
           get_psf, load_psf, delete_psf
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


    def get_psf(self, id=None):
        """
        get_psf

        SYNOPSIS
           Return a PSF model by Sherpa data id

        SYNTAX

        Arguments:
           id   - Sherpa data id
                  default = default data id

        Returns:
           None

        DESCRIPTION
           Return a PSF model from the instrument list by Sherpa data id.

        SEE ALSO
           set_psf, load_psf, delete_psf
        """
        return self._get_item(id, self._psf, 'psf model', 'has not been set')


    def delete_psf(self, id=None):
        """
        delete_psf

        SYNOPSIS
           Delete the specified PSF model by Sherpa data id.

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           None

        DESCRIPTION
           Delete the PSF model from the instrument list by Sherpa data id.

        SEE ALSO
           set_psf, get_psf, load_psf
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
        """
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

    def freeze(self, *args):
        """
        freeze

        SYNOPSIS
           Freeze a list of parameters

        SYNTAX

        Arguments:
           args      - list of Sherpa parameters

        Returns:
           None

        DESCRIPTION
           Freeze a list of Sherpa parameters.  Frozen parameters will not
           vary during a fit.

        SEE ALSO
           thaw, link, unlink
        """
        par = list(args)
        for p in par:
            self._freeze_thaw_par_or_model(p, 'freeze')

    def thaw(self, *args):
        """
        thaw

        SYNOPSIS
           Thaw a list of parameters

        SYNTAX

        Arguments:
           args      - list of Sherpa parameters

        Returns:
           None

        DESCRIPTION
           Thaw a list of Sherpa parameters.  Thawed parameters will vary
           during a fit.

        SEE ALSO
           freeze, link, unlink
        """
        par = list(args)
        for p in par:
            self._freeze_thaw_par_or_model(p, 'thaw')

    def link(self, par, val):
        """
        link

        SYNOPSIS
           Link a parameter with an associated value

        SYNTAX

        Arguments:
           par       - Sherpa parameter to be linked

           val       - value

        Returns:
           None

        DESCRIPTION
           Link a Sherpa parameter with a value.  Linked parameters will
           assume this value during the fit.

           Note: linked values may be variable or static.

        SEE ALSO
           freeze, thaw, unlink
        """
        par = self._check_par(par)
        if isinstance(val, basestring):
            val = self._eval_model_expression(val, 'parameter link')
        par.link = val

    def unlink(self, par):
        """
        unlink

        SYNOPSIS
           Unlink a parameter

        SYNTAX

        Arguments:
           par       - Sherpa parameter to be unlinked

        Returns:
           None

        DESCRIPTION
           Unlink a Sherpa parameter. See link for more information.

        SEE ALSO
           freeze, thaw, link
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


    def calc_stat_info(self):
        """
        calc_stat_info

        SYNOPSIS
           Shows calculated information of the goodness-of-fit

        SYNTAX

        Arguments:
           None

        Returns:
           None

        DESCRIPTION
           Shows calculated results of the current goodness-of-fit by data id.

           Example 1:

              calc_stat_info()
              Dataset               = Sherpa data id
              Statistic             = Fit statistic
              Fit statistic value   = Fit statistic value
              Data points           = Number of data points
              Degrees of freedom    = (number of points - number of thawed
                                       parameters)
              Probability [Q-value] = Null hypothesis probability
              Reduced statistic     = Reduced statistic value (statval/dof)

        SEE ALSO
           get_stat_info
        """
        output = self._get_stat_info()
        output = [statinfo.format() for statinfo in output]

        if len(output) > 1:
            info('\n\n'.join(output))
        else:
            info(output[0])


    def get_stat_info(self):
        """
        calc_stat_info

        SYNOPSIS
           Access calculated information of the goodness-of-fit

        SYNTAX

        Arguments:
           None

        Returns:
           List of StatInfoResults objects

        DESCRIPTION
           Return a list of calculated results of the current goodness-of-fit
           by data id.

           Example 1:

              print get_stat_info()[0]
              name      = Sherpa dataset label
              ids       = Sherpa dataset ids
              bkg_ids   = Sherpa background dataset ids
              statname  = Fit statistic
              statval   = Fit statistic value
              numpoints = Number of data points
              dof       = (numpoints - number of thawed parameters)
              qval      = Null hypothesis probability
              rstat     = Reduced statistic value (statval/dof)

        SEE ALSO
           calc_stat_info
        """
        return self._get_stat_info()


    def get_fit_results(self):
        """
        get_fit_results

        SYNOPSIS
           Return results from the last fit performed

        SYNTAX

        Arguments:
           None

        Returns:
           Sherpa fit_results object

        DESCRIPTION
           Printing out a Sherpa fit_results object displays results from a
           fit in table form.

           Example 1:
           
              print get_fit_results()
              succeeded = boolean of fit success
              parnames  = list of thawed parameter names
              parvals   = list of thawed parameter values
              statval   = statistic value
              numpoints = number of points on grid
              dof       = degrees of freedom
              qval      = probability                Note: N/A for Cash,CStat
              rstat     = reduced statistic value    Note: N/A for Cash,CStat
              message   = message from optimization method
              nfev      = number of function evalutions

           Example 2:

              print get_fit_results().format()
              Statistic value = statval at function evaluation nfev
              Data points = numpoints
              Degrees of freedom = dof
              Probability [Q-value] = qval
              Reduced statistic  = rstat

        SEE ALSO
           freeze, thaw, link
        """
        if self._fit_results == None:
            raise SessionErr('nofit', 'fit')
        else:
            return self._fit_results

    def guess(self, id=None, model=None, limits=True, values=True):
        """
        guess

        SYNOPSIS
           Guess the parameters values of model from the data

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           model     - Sherpa model to perform guess
                       default = None

        Returns:
           None

        DESCRIPTION
           Guess works only for single, non-composite models

           Examples:
              set_model(powlaw1d.p1)
              guess()

              set_model("src", powlaw1d.p1)
              guess("src")

              set_model(powlaw1d.p1*gauss1d.g1)
              guess(g1)
              guess(p1)

        SEE ALSO
           set_model, fit
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


    def calc_stat(self, id=None, *otherids):
        """
        calc_stat

        SYNOPSIS
           Return the statistic value

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - List of other Sherpa data ids

        Returns:
           Statistic value

        DESCRIPTION
           
        SEE ALSO
           calc_chisqr, get_stat, set_stat
        """
        ids, f = self._get_fit(id, otherids)
        return f.calc_stat()

    def calc_chisqr(self, id=None, *otherids):
        """
        calc_chisqr

        SYNOPSIS
           Return the chi squared statistic contribution by bin

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - List of other Sherpa data ids

        Returns:
           Statistic array

        DESCRIPTION
           NOTE:  Only available for Chi Squared statistics
           
        SEE ALSO
           calc_stat, get_stat, set_stat
        """
        ids, f = self._get_fit(id, otherids)
        return f.calc_chisqr()

    def fit(self, id=None, *otherids, **kwargs):
        """
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
           get_fit_results, conf, proj, covar
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

    def normal_sample(self, num=1, sigma=1, correlate=True,
                      id=None, otherids=(), numcores=None):
        """
        normal_sample

        SYNOPSIS
           Samples the current thawed parameter from a uni-variate or
           multi-variate normal distribution and computes the current
           statistic value on the set.

        SYNTAX

        Arguments:

           num         - Number of samples to calculate
                         default = 1

           sigma       - Spread of the normal distribution
                         default = 1

           correlate   - If True, sample from multi-variate normal using 
                         covariance; if False, sample from uni-variate normal.
                         default = True

           id          - Sherpa data id
                         default = default data id

           otherids    - List of other Sherpa data ids
                         default = ()

           numcores    - Number of cpus to use to calculate the statistic
                         default = number of detected cpus

        Returns:
           A NumPy array table with the first column representing the
           statistic and subsequent columns as the parameters.

        DESCRIPTION
           

        SEE ALSO
           uniform_sample, t_sample

        """
        ids, fit = self._get_fit(id, otherids)
        return sherpa.sim.normal_sample(fit, num, sigma, correlate, numcores)


    def uniform_sample(self, num=1, factor=4,
                       id=None, otherids=(), numcores=None):
        """
        uniform_sample

        SYNOPSIS
           Samples the current thawed parameter from a uniform distribution
           and computes the current statistic value on the set.

        SYNTAX

        Arguments:

           num         - Number of samples to calculate
                         default = 1

           factor      - Multiplier to expand the parameter scale
                         default = 4

           id          - Sherpa data id
                         default = default data id

           otherids    - List of other Sherpa data ids
                         default = ()

           numcores    - Number of cpus to use to calculate the statistic
                         default = number of detected cpus

        Returns:
           A NumPy array table with the first column representing the
           statistic and subsequent columns as the parameters.

        DESCRIPTION
           

        SEE ALSO
           normal_sample, t_sample

        """
        ids, fit = self._get_fit(id, otherids)
        return sherpa.sim.uniform_sample(fit, num, factor, numcores)


    def t_sample(self, num=1, dof=None, id=None, otherids=(), numcores=None):
        """
        t_sample

        SYNOPSIS
           Samples the current thawed parameter from a Student's t 
           distribution and computes the current statistic value on the set.

        SYNTAX

        Arguments:

           num         - Number of samples to calculate
                         default = 1

           dof         - Degrees of freedom, will calculate from current fit
                         by default.
                         default = None

           id          - Sherpa data id
                         default = default data id

           otherids    - List of other Sherpa data ids
                         default = ()

           numcores    - Number of cpus to use to calculate the statistic
                         default = number of detected cpus

        Returns:
           A NumPy array table with the first column representing the
           statistic and subsequent columns as the parameters.

        DESCRIPTION
           

        SEE ALSO
           normal_sample, uniform_sample
        """
        ids, fit = self._get_fit(id, otherids)
        if dof is None:
            dof = (len(fit.data.eval_model_to_fit(fit.model)) -
                   len(fit.model.thawedpars))
        return sherpa.sim.t_sample(fit, num, dof, numcores)


    ###########################################################################
    # Error estimation
    ###########################################################################


    def get_covar(self):
        """
        get_covar

        SYNOPSIS
           Access current covar estimation method object

        SYNTAX

        Arguments:
           None

        Returns:
           Current covar estimation method object

        DESCRIPTION
           Estimation method objects include the following attributes:

           * sigma                      - default = 1

           * eps                        - default = 0.01

           * maxiters                   - default = 200

           * remin                      - default = 0.01

           * soft_limits                - default = False

        SEE ALSO
           conf, proj, covar, get_covar_results, get_proj_results, get_proj,
           get_conf_results, get_conf
        """
        return self._estmethods['covariance']

    def get_conf(self):
        """
        get_conf

        SYNOPSIS
           Access current conf estimation method object

        SYNTAX

        Arguments:
           None

        Returns:
           Current conf estimation method object

        DESCRIPTION
           Estimation method objects include the following attributes:

           * sigma                      - default = 1

           * eps                        - default = 0.01

           * maxiters                   - default = 200

           * remin                      - default = 0.01

           * maxfits                    - default = 5
                      
           * max_rstat                  - default = 3

           * soft_limits                - default = False

           * fast                       - default = True

           * tol                        - default = 0.2

           * verbose                    - default = False           

        SEE ALSO
           conf, proj, covar, get_covar_results, get_conf_results,
           get_proj_results, get_covar
        """
        return self._estmethods['confidence']

    def get_proj(self):
        """
        get_proj

        SYNOPSIS
           Access current proj estimation method object

        SYNTAX

        Arguments:
           None

        Returns:
           Current proj estimation method object

        DESCRIPTION
           Estimation method objects include the following attributes:

           * sigma                      - default = 1

           * eps                        - default = 0.01

           * maxiters                   - default = 200

           * remin                      - default = 0.01

           * maxfits                    - default = 5
                      
           * max_rstat                  - default = 3

           * soft_limits                - default = False

           * fast                       - default = True

           * tol                        - default = 0.2

        SEE ALSO
           proj, covar, get_covar_results, get_proj_results, get_covar
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
        
    def get_covar_opt(self, name=None):
        """
        get_covar_opt

        SYNOPSIS
           Return a covariance option by name

        SYNTAX

        Arguments:
           name       - covariance option name

        Returns:
           covariance option value

        DESCRIPTION
           If given no argument, returns dictionary of all options
           that govern how covariance is run.  If given the name
           of an option, returns the value of that option.

        SEE ALSO
           covar, set_covar_opt
        """
        return self._get_estmethod_opt('covariance', name)

    def get_conf_opt(self, name=None):
        """
        get_conf_opt

        SYNOPSIS
           Return a confidence option by name

        SYNTAX

        Arguments:
           name       - confidence option name

        Returns:
           confidence option value

        DESCRIPTION
           If given no argument, returns dictionary of all options
           that govern how confidence is run.  If given the name
           of an option, returns the value of that option.

        SEE ALSO
           conf, set_conf_opt
        """
        return self._get_estmethod_opt('confidence', name)
    
    def get_proj_opt(self, name=None):
        """
        get_proj_opt

        SYNOPSIS
           Return a projection option by name

        SYNTAX

        Arguments:
           name       - projection option name

        Returns:
           projection option value

        DESCRIPTION
           If given no argument, returns dictionary of all options
           that govern how projection is run.  If given the name
           of an option, returns the value of that option.

        SEE ALSO
           proj, set_proj_opt
        """
        return self._get_estmethod_opt('projection', name)
    
    def set_covar_opt(self, name, val):
        """
        set_covar_opt
        
        SYNOPSIS
           Set a covariance option by name
        
        SYNTAX
        
        Arguments:
           name       - covariance option name
           val        - covariance option value
        
        Returns:
           None
        
        DESCRIPTION
           For the named covariance option, set that option to the new
           value.
        
        SEE ALSO
           covar, get_covar_opt
        """
        self._set_estmethod_opt('covariance', name, val)

    def set_conf_opt(self, name, val):
        """
        set_conf_opt
        
        SYNOPSIS
           Set a confidence option by name
        
        SYNTAX
        
        Arguments:
           name       - confidence option name
           val        - confidence option value
        
        Returns:
           None
        
        DESCRIPTION
           For the named confidence option, set that option to the new
           value.
        
        SEE ALSO
           conf, get_conf_opt
        """
        self._set_estmethod_opt('confidence', name, val)

    def set_proj_opt(self, name, val):
        """
        set_proj_opt
        
        SYNOPSIS
           Set a projection option by name
        
        SYNTAX
        
        Arguments:
           name       - projection option name
           val        - projection option value
        
        Returns:
           None
        
        DESCRIPTION
           For the named projection option, set that option to the new
           value.
        
        SEE ALSO
           proj, get_proj_opt
        """
        self._set_estmethod_opt('projection', name, val)

    # Functions to get confidence limit results after last run
    # of named confidence limit function.
    def get_covar_results(self):
        """
        get_covar_results

        SYNOPSIS
           Access covariance estimation results object

        SYNTAX

        Arguments:
           None

        Returns:
           Covar estimation results object

        DESCRIPTION
           Access results from the last time covar was run.  The results
           include the following attributes:

           * datasets                        - Data sets in fit
           
           * methodname                      - Estimation method name

           * fitname                         - Fitting method name

           * statname                        - Statistic name

           * sigma                           - Change in statistic

           * parnames                        - Model parameter names

           * parvals                         - Model parameter fit values

           * parmins                         - Model parameter minimum values

           * parmaxes                        - Model parameter maximum values

           * warnings                        - Warning messages

        SEE ALSO
           conf, proj, covar, get_proj, get_conf_results,
           get_proj_results, get_covar
        """
        if self._covariance_results == None:
            raise SessionErr('noaction', "covariance")
        else:
            return self._covariance_results

    def get_conf_results(self):
        """
        get_conf_results

        SYNOPSIS
           Access confidence estimation results object

        SYNTAX

        Arguments:
           None

        Returns:
           Conf estimation results object

        DESCRIPTION
           Access results from the last time confidence was run.  The results
           include the following attributes:

           * datasets                        - Data sets in fit
           
           * methodname                      - Estimation method name

           * fitname                         - Fitting method name

           * statname                        - Statistic name

           * sigma                           - Change in statistic

           * parnames                        - Model parameter names

           * parvals                         - Model parameter fit values

           * parmins                         - Model parameter minimum values

           * parmaxes                        - Model parameter maximum values

           * warnings                        - Warning messages

        SEE ALSO
           conf, proj, covar, get_conf, get_proj, get_covar_results, get_covar
        """
        if self._confidence_results == None:
            raise SessionErr('noaction', "confidence")
        else:
            return self._confidence_results

    def get_proj_results(self):
        """
        get_proj_results

        SYNOPSIS
           Access projection estimation results object

        SYNTAX

        Arguments:
           None

        Returns:
           Proj estimation results object

        DESCRIPTION
           Access results from the last time projection was run.  The results
           include the following attributes:

           * datasets                        - Data sets in fit
           
           * methodname                      - Estimation method name

           * fitname                         - Fitting method name

           * statname                        - Statistic name

           * sigma                           - Change in statistic

           * parnames                        - Model parameter names

           * parvals                         - Model parameter fit values

           * parmins                         - Model parameter minimum values

           * parmaxes                        - Model parameter maximum values

           * warnings                        - Warning messages

        SEE ALSO
           proj, covar, get_proj, get_covar_results, get_covar
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

    def covar(self, *args):
        """
        covar

        SYNOPSIS
           Estimate confidence intervals for selected thawed parameters

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - Other Sherpa data ids

           pars      - Sherpa model parameters

        Returns:
           Formatted estimation results output

        DESCRIPTION

        SEE ALSO
           proj, get_proj, get_covar, get_covar_results,
           get_proj_results
        """
        self._covariance_results = self._est_errors(args, 'covariance')
    
    def conf(self, *args):
        """
        conf

        SYNOPSIS
           Estimate confidence intervals for selected thawed parameters

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - Other Sherpa data ids

           pars      - Sherpa model parameters

        Returns:
           Formatted estimation results output

        DESCRIPTION
           Confidence interval bounds are determined for each selected
           parameter in turn. A given parameter's value is varied along a grid
           of values while the values of all the other nominally thawed
           parameters are allowed to float to new best-fit values (compare to
           covar, where the values of all the other nominally thawed parameters
           remain fixed to their best-fit values). This method of estimating
           confidence interval bounds gives truly accurate results only in
           special cases:

           An estimated confidence interval is accurate if and only if:

           1. the chi^2 or logL surface in parameter space is approximately
              shaped like a multi-dimensional paraboloid, and
           2. the best-fit point is sufficiently far from parameter space
              boundaries.

           One may determine if these conditions hold, for example, by plotting
           the fit statistic as a function of each parameter's values (the
           curve should approximate a parabola) and by examining contour plots
           of the fit statistics made by varying the values of two parameters
           at a time (the contours should be elliptical, and parameter space
           boundaries should be no closer than approximately 3 sigma from the
           best-fit point).

           If no arguments are given this function, it is assumed that the
           data id is the default data id, and that limits should be computed
           for all thawed parameters.

           If arguments are given, each argument is examined to see if it
           is a Sherpa model parameter.  If model parameters are given, it
           is assumed that limits for only those parameters should be
           computed.  Any argument that is not a model parameter is assumed
           to be a data id.

        SEE ALSO
           covar, get_conf, get_proj, get_covar, get_covar_results,
           get_conf_results, get_proj_results
        """
        self._confidence_results = self._est_errors(args, 'confidence')

    def proj(self, *args):
        """
        proj

        SYNOPSIS
           Estimate confidence intervals for selected thawed parameters

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - Other Sherpa data ids

           pars      - Sherpa model parameters

        Returns:
           Formatted estimation results output

        DESCRIPTION
           Confidence interval bounds are determined for each selected
           parameter in turn. A given parameter's value is varied along a grid
           of values while the values of all the other nominally thawed
           parameters are allowed to float to new best-fit values (compare to
           covar, where the values of all the other nominally thawed parameters
           remain fixed to their best-fit values). This method of estimating
           confidence interval bounds gives truly accurate results only in
           special cases:

           An estimated confidence interval is accurate if and only if:

           1. the chi^2 or logL surface in parameter space is approximately
              shaped like a multi-dimensional paraboloid, and
           2. the best-fit point is sufficiently far from parameter space
              boundaries.

           One may determine if these conditions hold, for example, by plotting
           the fit statistic as a function of each parameter's values (the
           curve should approximate a parabola) and by examining contour plots
           of the fit statistics made by varying the values of two parameters
           at a time (the contours should be elliptical, and parameter space
           boundaries should be no closer than approximately 3 sigma from the
           best-fit point).

           If no arguments are given this function, it is assumed that the
           data id is the default data id, and that limits should be computed
           for all thawed parameters.

           If arguments are given, each argument is examined to see if it
           is a Sherpa model parameter.  If model parameters are given, it
           is assumed that limits for only those parameters should be
           computed.  Any argument that is not a model parameter is assumed
           to be a data id.

        SEE ALSO
           covar, get_proj, get_covar, get_covar_results,
           get_proj_results
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
    ###########################################################################


    def set_sampler_opt(self, opt, value):
        self._pyblocxs.set_sampler_opt(opt, value)

    def get_sampler_opt(self, opt):
        return self._pyblocxs.get_sampler_opt(opt)

    def get_sampler_name(self):
        return self._pyblocxs.get_sampler_name()
        
    def set_sampler(self, sampler):
        self._pyblocxs.set_sampler(sampler)

    def get_sampler(self):
        return self._pyblocxs.get_sampler()

    def set_prior(self, par, prior):
        self._pyblocxs.set_prior(par, prior)

    def get_prior(self, par):
        return self._pyblocxs.get_prior(par)

    def list_priors(self):
        return self._pyblocxs.list_priors()

    def list_samplers(self):
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

    def get_data_plot(self, id=None):
        """
        get_data_plot

        SYNOPSIS
           Return a Sherpa data plot

        SYNTAX

        Arguments:
           id        - Sherpa data id

        Returns:
           Sherpa DataPlot object

        DESCRIPTION
           The Sherpa data plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

           Functions:

              prepare()
                 populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_data
        """
        self._prepare_plotobj(id, self._dataplot)
        return self._dataplot

    def get_data_plot_prefs(self):
        """
        get_data_plot_prefs

        SYNOPSIS
           Return data plot preferences

        SYNTAX

        Arguments:
           None

        Returns:
           Dictionary of data plot preferences

        DESCRIPTION
              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - 'line'
                 errthickness   - None
                 linecolor      - None
                 linestyle      - 0
                 linethickness  - None
                 ratioline      - N/A
                 symbolcolor    - None
                 symbolfill     - False
                 symbolsize     - 3
                 symbolstyle    - 4
                 xaxis          - N/A
                 xerrorbars     - False
                 xlog           - False
                 yerrorbars     - True
                 ylog           - False

           Examples:

               get_data_plot_prefs()
           {'errstyle': 'line', 'symbolfill': False, 'symbolstyle': 4,
            'linestyle': 0, 'symbolsize': 3, 'yerrorbars': True}

               get_data_plot_prefs()['xlog']=True

               get_data_plot_prefs()
           {'errstyle': 'line', 'symbolfill': False, 'symbolstyle': 4,
            'xlog': True, 'linestyle': 0, 'symbolsize': 3, 'yerrorbars': True}


        SEE ALSO
           plot_data, get_data_plot
        """
        return self._dataplot.plot_prefs

    def get_model_plot(self, id=None):
        """
        get_model_plot

        SYNOPSIS
           Return a Sherpa model plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa ModelPlot object

        DESCRIPTION
           The Sherpa model plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

           Functions:

              prepare()
                 calculate the model and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           plot_model
        """
        self._prepare_plotobj(id, self._modelplot)
        return self._modelplot

    def get_source_plot(self, id=None):
        """
        get_source_plot

        SYNOPSIS
           Return a Sherpa source plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa SourcePlot object

        DESCRIPTION
           The Sherpa source plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

           Functions:

              prepare()
                 calculate the source and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           plot_source
        """
        self._prepare_plotobj(id, self._sourceplot)
        return self._sourceplot


    def get_model_component_plot(self, id, model=None):
        """
        get_model_component_plot

        SYNOPSIS
           Return a Sherpa model component plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           model     - Sherpa model component

        Returns:
           Sherpa ComponentModelPlot object

        DESCRIPTION
           The Sherpa model component plot object holds references to 
           various plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

           Functions:

              prepare()
                 calculate the source and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           plot_model_component, plot_source_component
        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        self._prepare_plotobj(id, self._compmdlplot, model=model)
        return self._compmdlplot


    def get_source_component_plot(self, id, model=None):
        """
        get_source_component_plot

        SYNOPSIS
           Return a Sherpa source model component plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           model     - Sherpa source model component

        Returns:
           Sherpa ComponentSourcePlot object

        DESCRIPTION
           The Sherpa source model component plot object holds references to 
           various plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

           Functions:

              prepare()
                 calculate the source and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           plot_model_component, plot_source_component
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


    def get_model_plot_prefs(self):
        """
        get_model_plot_prefs

        SYNOPSIS
           Return model plot preferences

        SYNTAX

        Arguments:
           None

        Returns:
           Dictionary of model plot preferences

        DESCRIPTION
              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - None
                 errthickness   - None
                 linecolor      - 'red'
                 linestyle      - 1
                 linethickness  - 3
                 ratioline      - N/A
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - None
                 symbolstyle    - 0
                 xaxis          - N/A
                 xerrorbars     - False
                 xlog           - False
                 yerrorbars     - False
                 ylog           - False

           Examples:

               get_model_plot_prefs()
           {'symbolstyle': 0, 'linethickness': 3, 'linestyle': 1,
            'linecolor': 'red'}

               get_model_plot_prefs()['linecolor']='yellow'

               get_model_plot_prefs()
           {'symbolstyle': 0, 'linethickness': 3, 'linestyle': 1,
            'linecolor': 'yellow'}


        SEE ALSO
           plot_model, get_model_plot
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

    def get_resid_plot(self, id=None):
        """
        get_resid_plot

        SYNOPSIS
           Return a Sherpa resid plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa ResidPlot object

        DESCRIPTION
           The Sherpa resid plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - 'line'
                 errthickness   - None
                 linecolor      - None
                 linestyle      - 0
                 linethickness  - None
                 ratioline      - False
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - 3
                 symbolstyle    - 2
                 xaxis          - True
                 xerrorbars     - True
                 xlog           - False
                 yerrorbars     - True
                 ylog           - False

           Functions:

              prepare()
                 calculate the residuals and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_resid
        """
        self._prepare_plotobj(id, self._residplot)
        return self._residplot

    def get_delchi_plot(self,id=None):
        """
        get_delchi_plot

        SYNOPSIS
           Return a Sherpa delta chi plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa DelchiPlot object

        DESCRIPTION
           The Sherpa delta chi plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - 'line'
                 errthickness   - None
                 linecolor      - None
                 linestyle      - 0
                 linethickness  - None
                 ratioline      - False
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - 3
                 symbolstyle    - 2
                 xaxis          - True
                 xerrorbars     - True
                 xlog           - False
                 yerrorbars     - True
                 ylog           - False

           Functions:

              prepare()
                 calculate the delta chi and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_delchi
        """
        self._prepare_plotobj(id, self._delchiplot)
        return self._delchiplot

    def get_chisqr_plot(self,id=None):
        """
        get_chisqr_plot

        SYNOPSIS
           Return a Sherpa chi square plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa ChisqrPlot object

        DESCRIPTION
           The Sherpa chi square plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - None
                 linecolor      - 'red'
                 errthickness   - None
                 linestyle      - 1
                 linethickness  - 3
                 ratioline      - N/A
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - None
                 symbolstyle    - 0
                 xaxis          - N/A
                 xerrorbars     - False
                 xlog           - False
                 yerrorbars     - False
                 ylog           - False

           Functions:

              prepare()
                 calculate the chi square and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_chisqr
        """
        self._prepare_plotobj(id, self._chisqrplot)
        return self._chisqrplot
    
    def get_ratio_plot(self, id=None):
        """
        get_ratio_plot

        SYNOPSIS
           Return a Sherpa ratio plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa RatioPlot object

        DESCRIPTION
           The Sherpa ratio plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - 'line'
                 errthickness   - None
                 linecolor      - None
                 linestyle      - 0
                 linethickness  - None
                 ratioline      - True
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - 3
                 symbolstyle    - 2
                 xaxis          - False
                 xerrorbars     - True
                 xlog           - False
                 yerrorbars     - True
                 ylog           - False                 

           Functions:

              prepare()
                 calculate the ratios and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_ratio
        """
        self._prepare_plotobj(id, self._ratioplot)
        return self._ratioplot
    
    def get_data_contour(self, id=None):
        """
        get_data_contour

        SYNOPSIS
           Return a Sherpa data contour

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa DataContour object

        DESCRIPTION
           The Sherpa data contour object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x0           - independent variable array

              x1           - independent variable array

              y            - dependent variable array

              levels       - list of contour slices 

           Functions:

              prepare()
                 populate the data arrays

              contour( overcontour=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           contour_data
        """
        self._prepare_plotobj(id, self._datacontour)
        return self._datacontour

    def get_data_contour_prefs(self):
        """
        get_data_contour_prefs

        SYNOPSIS
           Return data contour preferences

        SYNTAX

        Arguments:
           None

        Returns:
           Dictionary of data contour preferences

        DESCRIPTION
              contour_prefs   - dictionary of plotting preferences

                 color      - None
                 thickness  - None
                 style      - None
                 xlog       - False
                 ylog       - False

           Examples:

               get_data_contour_prefs()
           {}

               get_data_contour_prefs()['color']='blue'

               get_data_contour_prefs()
           {'color': 'blue'}


        SEE ALSO
           contour_data, get_data_contour
        """
        return self._datacontour.contour_prefs

    def get_model_contour(self, id=None):
        """
        get_model_contour

        SYNOPSIS
           Return a Sherpa model contour

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa ModelContour object

        DESCRIPTION
           The Sherpa model contour object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x0           - independent variable array

              x1           - independent variable array

              y            - dependent variable array

              levels       - list of contour slices 

           Functions:

              prepare()
                 populate the data arrays

              contour( overcontour=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           contour_model
        """
        self._prepare_plotobj(id, self._modelcontour)
        return self._modelcontour

    def get_source_contour(self, id=None):
        """
        get_source_contour

        SYNOPSIS
           Return a Sherpa source contour

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa SourceContour object

        DESCRIPTION
           The Sherpa source contour object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x0           - independent variable array

              x1           - independent variable array

              y            - dependent variable array

              levels       - list of contour slices 

           Functions:

              prepare()
                 populate the data arrays

              contour( overcontour=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           contour_source
        """
        self._prepare_plotobj(id, self._sourcecontour)
        return self._sourcecontour

    def get_model_contour_prefs(self):
        """
        get_model_contour_prefs

        SYNOPSIS
           Return model contour preferences

        SYNTAX

        Arguments:
           None

        Returns:
           Dictionary of model contour preferences

        DESCRIPTION
              contour_prefs   - dictionary of plotting preferences

                 color      - 'red'
                 thickness  - None
                 style      - None
                 xlog       - False
                 ylog       - False

           Examples:

               get_model_contour_prefs()
           {'color': 'red', 'style': None, 'thickness': 3}

               get_model_contour_prefs()['xlog']=True

               get_model_contour_prefs()['ylog']=True

               get_model_contour_prefs()
           {'color': 'red', 'style': None, 'ylog': True, 'xlog': True,
            'thickness': 3}


        SEE ALSO
           contour_model, get_model_contour
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

    def get_resid_contour(self, id=None):
        """
        get_resid_contour

        SYNOPSIS
           Return a Sherpa resid contour

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa ResidContour object

        DESCRIPTION
           The Sherpa resid contour object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x0           - independent variable array

              x1           - independent variable array

              y            - dependent variable array

              levels       - list of contour slices 

              contour_prefs   - dictionary of plotting preferences

                 color      - None
                 thickness  - None
                 style      - None
                 xlog       - False
                 ylog       - False

           Functions:

              prepare()
                 populate the data arrays

              contour( overcontour=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           contour_resid
        """
        self._prepare_plotobj(id, self._residcontour)
        return self._residcontour

    def get_ratio_contour(self, id=None):
        """
        get_ratio_contour

        SYNOPSIS
           Return a Sherpa ratio contour

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa RatioContour object

        DESCRIPTION
           The Sherpa ratio contour object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x0           - independent variable array

              x1           - independent variable array

              y            - dependent variable array

              levels       - list of contour slices 

              contour_prefs   - dictionary of plotting preferences

                 color      - None
                 thickness  - None
                 style      - None
                 xlog       - False
                 ylog       - False

           Functions:

              prepare()
                 populate the data arrays

              contour( overcontour=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           contour_ratio
        """
        self._prepare_plotobj(id, self._ratiocontour)
        return self._ratiocontour

    def get_psf_contour(self, id=None):
        """
        get_psf_contour

        SYNOPSIS
           Return a Sherpa PSF contour

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa PSFContour object

        DESCRIPTION
           The Sherpa PSF contour object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x0           - independent variable array

              x1           - independent variable array

              y            - dependent variable array

              levels       - list of contour slices 

              contour_prefs   - dictionary of plotting preferences

                 color      - None or 'red'
                 thickness  - None
                 style      - None
                 xlog       - False
                 ylog       - False

           Functions:

              prepare()
                 populate the data arrays

              contour( overcontour=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           contour_psf
        """
        self._prepare_plotobj(id, self._psfcontour)
        return self._psfcontour


    def get_kernel_contour(self, id=None):
        """
        get_kernel_contour

        SYNOPSIS
           Return a Sherpa PSF kernel contour

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa PSFKernelContour object

        DESCRIPTION
           The Sherpa PSF kernel contour object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x0           - independent variable array

              x1           - independent variable array

              y            - dependent variable array

              levels       - list of contour slices 

              contour_prefs   - dictionary of plotting preferences

                 color      - None or 'red'
                 thickness  - None
                 style      - None
                 xlog       - False
                 ylog       - False

           Functions:

              prepare()
                 populate the data arrays

              contour( overcontour=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           contour_kernel
        """
        self._prepare_plotobj(id, self._kernelcontour)
        return self._kernelcontour


    def get_psf_plot(self, id=None):
        """
        get_psf_plot

        SYNOPSIS
           Return a Sherpa PSF plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa PSFPlot object

        DESCRIPTION
           The Sherpa PSF plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - 'line' or None
                 errthickness   - None
                 linecolor      - None or 'red'
                 linestyle      - 0 or 1
                 linethickness  - None or 3
                 ratioline      - N/A
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - 3 or None
                 symbolstyle    - 4 or 0
                 xaxis          - N/A
                 xerrorbars     - False
                 xlog           - False
                 yerrorbars     - True or False
                 ylog           - False

           Functions:

              prepare()
                 populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_psf, plot_data, plot_model
        """
        self._prepare_plotobj(id, self._psfplot)
        return self._psfplot


    def get_kernel_plot(self, id=None):
        """
        get_kernel_plot

        SYNOPSIS
           Return a Sherpa PSF kernel plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa PSFKernelPlot object

        DESCRIPTION
           The Sherpa PSF kernel plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - 'line' or None
                 errthickness   - None
                 linecolor      - None or 'red'
                 linestyle      - 0 or 1
                 linethickness  - None or 3
                 ratioline      - N/A
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - 3 or None
                 symbolstyle    - 4 or 0
                 xaxis          - N/A
                 xerrorbars     - False
                 xlog           - False
                 yerrorbars     - True or False
                 ylog           - False

           Functions:

              prepare()
                 populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_kernel, plot_data, plot_model
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


    def set_xlog(self, plottype="all"):
        self._set_plot_item(plottype, 'xlog', True)


    def set_ylog(self, plottype="all"):
        self._set_plot_item(plottype, 'ylog', True)


    def set_xlinear(self, plottype="all"):
        self._set_plot_item(plottype, 'xlog', False)


    def set_ylinear(self, plottype="all"):
        self._set_plot_item(plottype, 'ylog', False)


    def plot(self, *args):
        """
        plot

        SYNOPSIS
           Send a combination plot to the visualizer

        SYNTAX

        Arguments:
           plot0       - string of first plot type

           id0         - Sherpa data id
                         default = default data id

           ...

           plotn       - string of nth plot type

           idn         - Sherpa data id
                         default = default data id

        Returns:
           None

        DESCRIPTION
           Visualize multiple plots by Sherpa data ids.

           Applicable types include: 'data', 'model', 'fit', 'resid',
                                     'ratio', 'delchi', 'chisqr','psf'

           Example 1:

               plot('data', 'model')

           Example 2: using ids

               plot('data', 1, 'model', 1)

        SEE ALSO
           plot_data, plot_model, plot_fit, plot_resid, plot_ratio,
           plot_delchi, plot_chisqr
        """
        self._multi_plot(args)        

    def plot_data(self, id=None, **kwargs):
        """
        plot_data

        SYNOPSIS
           Send a data plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a dataset by Sherpa data id.

        SEE ALSO
           get_data_plot, plot_model, plot_fit, plot_fit_resid, plot_fit_delchi
        """
        self._plot(id, self._dataplot, **kwargs)
    
    def plot_model(self, id=None, **kwargs):
        """
        plot_model

        SYNOPSIS
           Send a model plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a dataset model by Sherpa data id.

        SEE ALSO
           get_model_plot, plot_data, plot_fit, plot_fit_resid, plot_fit_delchi
        """
        self._plot(id, self._modelplot, **kwargs)


    def plot_source_component(self, id, model=None, **kwargs):
        """
        plot_source_component

        SYNOPSIS
           Send a source model component plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           model       - Sherpa model component or expression string

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a source model component by Sherpa data id.

        SEE ALSO
           get_model_component_plot, plot_model_component, plot_data,
           plot_fit, plot_fit_resid, plot_fit_delchi
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


    def plot_model_component(self, id, model=None, **kwargs):
        """
        plot_model_component

        SYNOPSIS
           Send a model component plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           model       - Sherpa model component or expression string

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a model component by Sherpa data id.

        SEE ALSO
           get_model_component_plot, plot_source_component, plot_data,
           plot_fit, plot_fit_resid, plot_fit_delchi
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


    def plot_source(self, id=None, **kwargs):
        """
        plot_source

        SYNOPSIS
           Send a source plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a dataset source by Sherpa data id.

        SEE ALSO
           get_source_plot, plot_data, plot_fit, plot_fit_resid, plot_fit_delchi
        """
        id = self._fix_id(id)
        mdl = self._models.get(id, None)
        if mdl is not None:
            raise IdentifierErr("Convolved model\n'%s'\n is set for dataset %s. You should use plot_model instead." % 
                    (mdl.name, str(id)))
        self._plot(id, self._sourceplot, **kwargs)

    def plot_fit(self, id=None, **kwargs):
        """
        plot_fit

        SYNOPSIS
           Send a fit plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a dataset and dataset model by Sherpa data id.

        SEE ALSO
           get_fit_plot, plot_model, plot_data, plot_fit_resid, plot_fit_delchi
        """
        self._plot(id, self._fitplot, **kwargs)

    def plot_resid(self, id=None, **kwargs):
        """
        plot_resid

        SYNOPSIS
           Send a residuals plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the residuals (dataset minus dataset model) by Sherpa data
           id.

        SEE ALSO
           get_resid_plot, plot_ratio, plot_delchi, plot_fit_resid,
           plot_fit_delchi, plot_chisqr, plot_fit, plot_data, plot_model
        """
        self._plot(id, self._residplot, **kwargs)

    def plot_chisqr(self, id=None, **kwargs):
        """
        plot_chisqr

        SYNOPSIS
           Send a chi^2 plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the chi^2 (residuals divided by dataset uncertainties, the
           quantity squared) by Sherpa data id.

        SEE ALSO
           get_chisqr_plot, plot_resid, plot_ratio, plot_fit_resid,
           plot_fit_delchi, plot_delchi, plot_fit, plot_data, plot_model
        """
        self._plot(id, self._chisqrplot, **kwargs)

    def plot_delchi(self, id=None, **kwargs):
        """
        plot_delchi

        SYNOPSIS
           Send a delta chi plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the delta chi (residuals divided by dataset uncertainties)
           by Sherpa data id.

        SEE ALSO
           get_delchi_plot, plot_resid, plot_ratio, plot_fit_resid,
           plot_fit_delchi, plot_chisqr, plot_fit, plot_data, plot_model
        """
        self._plot(id, self._delchiplot, **kwargs)
        
    def plot_ratio(self, id=None, **kwargs):
        """
        plot_ratio

        SYNOPSIS
           Send a ratio plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the ratio (dataset divided by dataset model) by Sherpa
           data id.

        SEE ALSO
           get_ratio_plot, plot_resid, plot_delchi, plot_fit_resid,
           plot_fit_delchi, plot_chisqr, plot_fit, plot_data, plot_model
        """
        self._plot(id, self._ratioplot, **kwargs)

    def plot_psf(self, id=None, **kwargs):
        """
        plot_psf

        SYNOPSIS
           Send a PSF plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the PSF dataset or PSF model by Sherpa data id.

        SEE ALSO
           get_psf_plot, plot_data
        """
        self._plot(id, self._psfplot, **kwargs)


    def plot_kernel(self, id=None, **kwargs):
        """
        plot_kernel

        SYNOPSIS
           Send a PSF kernel plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the extracted sub-kernel dataset or extracted kernel model 
           by Sherpa data id.

        SEE ALSO
           get_psf_plot, plot_data
        """
        self._plot(id, self._kernelplot, **kwargs)


    def plot_fit_resid(self, id=None, replot=False, overplot=False,
                       clearwindow=True):
        """
        plot_fit_resid

        SYNOPSIS
           Send fit and residuals plots to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the fit plot and residuals plot in a joint plot
           window by Sherpa data id.

        SEE ALSO
           plot_resid, plot_delchi, plot_ratio, plot_chisqr, plot_fit,
           plot_data, plot_model, plot_fit_delchi
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

    def plot_fit_delchi(self, id=None, replot=False, overplot=False,
                        clearwindow=True):
        """
        plot_fit_delchi

        SYNOPSIS
           Send fit and delta chi plots to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the fit plot and delta chi plot in a joint plot
           window by Sherpa data id.

        SEE ALSO
           plot_resid, plot_delchi, plot_ratio, plot_chisqr, plot_fit,
           plot_data, plot_model, plot_fit_resid
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

    def plot_pdf(self, points, name="x", xlabel="x", bins=12, normed=True, 
                 replot=False, overplot=False, clearwindow=True ):

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
        return self._pdfplot


    def plot_cdf(self, points, name="x", xlabel="x", 
                 replot=False, overplot=False, clearwindow=True ):

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
        return self._cdfplot


    def plot_trace(self, points, name="x", xlabel="x", 
                   replot=False, overplot=False, clearwindow=True ):

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
        return self._traceplot


    def plot_scatter(self, x, y, name="(x,y)", xlabel="x", ylabel="y",
                   replot=False, overplot=False, clearwindow=True ):

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

    def contour(self, *args):
        """
        contour

        SYNOPSIS
           Send a combination contour plot to the visualizer

        SYNTAX

        Arguments:
           contour0    - string of first plot type

           id0         - Sherpa data id
                         default = default data id

           ...

           contourn    - string of nth plot type

           idn         - Sherpa data id
                         default = default data id

        Returns:
           None

        DESCRIPTION
           Visualize multiple contour plots by Sherpa data ids.

           Applicable types include: 'data', 'model', 'fit', 'resid',
                                     'ratio', 'psf'

           Example 1:

               contour('data', 'model')

           Example 2: using ids

               contour('data', 1, 'model', 1)

        SEE ALSO
           contour_fit, contour_data, contour_model, contour_resid,
           contour_ratio, contour_fit_resid
        """
        self._multi_plot(args, 'contour')

    def contour_data(self, id=None, **kwargs):
        """
        contour_data

        SYNOPSIS
           Send a data contour plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overcontour - Contour data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a dataset by Sherpa data id.

        SEE ALSO
           get_data_contour, contour_model, contour_fit, contour_fit_resid
        """
        self._contour(id, self._datacontour, **kwargs)
        
    def contour_model(self, id=None, **kwargs):
        """
        contour_model

        SYNOPSIS
           Send a model contour plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overcontour - Contour data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a dataset model by Sherpa data id.

        SEE ALSO
           get_model_contour, contour_data, contour_fit, contour_fit_resid
        """
        self._contour(id, self._modelcontour, **kwargs)

    def contour_source(self, id=None, **kwargs):
        """
        contour_source

        SYNOPSIS
           Send a source contour plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overcontour - Contour data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a dataset source by Sherpa data id.

        SEE ALSO
           get_source_contour, contour_data, contour_fit, contour_fit_resid
        """
        self._contour(id, self._sourcecontour, **kwargs)

    def contour_fit(self, id=None, **kwargs):
        """
        contour_fit

        SYNOPSIS
           Send a fit contour plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overcontour - Contour data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a dataset and dataset model by Sherpa data id.

        SEE ALSO
           get_fit_contour, contour_model, contour_data, contour_fit_resid
        """
        self._contour(id, self._fitcontour, **kwargs)

    def contour_resid(self, id=None, **kwargs):
        """
        contour_resid

        SYNOPSIS
           Send a residuals contour plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overcontour - Contour data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the residuals (dataset minus dataset model) by Sherpa data
           id.

        SEE ALSO
           get_resid_contour, contour_ratio, contour_fit_resid, contour_fit,
           contour_data, contour_model
        """
        self._contour(id, self._residcontour, **kwargs)
    
    def contour_ratio(self, id=None, **kwargs):
        """
        contour_ratio

        SYNOPSIS
           Send a ratio plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overcontour - Contour data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the ratio (dataset divided by dataset model) by Sherpa
           data id.

        SEE ALSO
           get_ratio_contour, contour_resid, contour_fit_resid, contour_fit,
           contour_data, contour_model
        """
        self._contour(id, self._ratiocontour, **kwargs)

    def contour_psf(self, id=None, **kwargs):
        """
        contour_psf

        SYNOPSIS
           Send a PSF contour to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overcontour - Contour data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the PSF image or PSF model contour by Sherpa data id.

        SEE ALSO
           get_ratio_contour, contour_resid, contour_fit_resid, contour_fit,
           contour_data, contour_model
        """
        self._contour(id, self._psfcontour, **kwargs)


    def contour_kernel(self, id=None, **kwargs):
        """
        contour_kernel

        SYNOPSIS
           Send a PSF kernel contour to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overcontour - Contour data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the PSF sub-kernel contour by Sherpa data id.

        SEE ALSO
           get_ratio_contour, contour_resid, contour_fit_resid, contour_fit,
           contour_data, contour_model
        """
        self._contour(id, self._kernelcontour, **kwargs)


    def contour_fit_resid(self, id=None, replot=False, overcontour=False):
        """
        contour_fit_resid

        SYNOPSIS
           Send fit and residual contours plot to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           replot      - Send cached data arrays to visualizer
                         default = False

           overcontour - Contour data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the fit contour and residuals contour in a joint contour
           window by Sherpa data id.

        SEE ALSO
           contour_resid, contour_ratio, contour_fit, contour_data,
           contour_model
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

    def get_int_proj(self, par=None, id=None, otherids=None, recalc=False,
                     min=None, max=None, nloop=20, delv=None, fac=1, 
                     log=False, numcores=None):
        """
        get_int_proj

        SYNOPSIS
           Return a confidence plot of fit statistic vs. a thawed parameter
           value.  At each step a fit is performed to obtain a new statistic
           if other thawed parameter(s) exist in the source model, otherwise,
           calculate the statistic (see get_int_unc).

        SYNTAX

        Arguments:
           par       - source model parameter
                       default = None

           id        - Sherpa data id
                       default = default data id

           otherids  - list of ids required for simultaneous fit
                       default=None

           recalc    - calculate confidence data
                       default=False

           min       - minimum bound
                       default=None

           max       - maximum bound
                       default=None

           nloop     - bin size, used in calculating stepsize
                       default=20

           delv      - stepsize, calculated by default
                       default=None

           fac       - factor used to expand or condense interval,
                       default=1

           log       - boolean to use log space for interval
                       default=False

           numcores  - specify the number of cores for parallel processing.
                       All available cores are used by default.
                       default=None

        Returns:
           int_proj object

        DESCRIPTION

           Example: for users who do not want to create plots

               print get_int_proj( par, recalc=True )

        SEE ALSO
           int_unc, reg_proj, reg_unc, get_int_inc
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

    def get_int_unc(self, par=None, id=None, otherids=None, recalc=False,
                    min=None, max=None, nloop=20, delv=None, fac=1, log=False,
                    numcores=None):
        """
        get_int_unc

        SYNOPSIS
           Return a confidence plot of fit statistic vs. parameter value.  At
           each step calculate the statistic with the other parameter(s) frozen
           at best fit values.

        SYNTAX

        Arguments:
           par       - source model parameter
                       default = None

           id        - Sherpa data id
                       default = default data id

           otherids  - list of ids required for simultaneous fit
                       default=None

           recalc    - calculate confidence data
                       default=False

           min       - minimum bound
                       default=None

           max       - maximum bound
                       default=None

           nloop     - bin size, used in calculating stepsize
                       default=20

           delv      - stepsize, calculated by default
                       default=None

           fac       - factor used to expand or condense interval,
                       default=1

           log       - boolean to use log space for interval
                       default=False

           numcores  - specify the number of cores for parallel processing.
                       All available cores are used by default.
                       default=None

        Returns:
           int_unc object

        DESCRIPTION

           Example: for users who do not want to create plots

               print get_int_unc( par, recalc=True )

        SEE ALSO
           int_proj, reg_proj, reg_unc, get_int_proj
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

    def get_reg_proj(self, par0=None, par1=None, id=None, otherids=None,
                     recalc=False, fast=True, min=None, max=None, 
                     nloop=(10,10),delv=None, fac=4, log=(False,False),
                     sigma=(1,2,3), levels=None, numcores=None):
        """
        get_reg_proj

        SYNOPSIS
           Return a confidence contour of fit statistic vs. two thawed
           parameter values.  At each step a fit is performed to obtain a new
           statistic if other thawed parameter(s) exist in the source model,
           otherwise, calculate the statistic (see get_reg_unc).

        SYNTAX

        Arguments:
           par0      - first source model parameter
                     - default = None

           par1      - second source model parameter
                       default = None

           id        - Sherpa data id
                       default = default data id

           otherids  - list of ids required for simultaneous fit
                       default=None

           recalc    - calculate confidence data
                       default=False

           fast      - change opt method to levmar for Chi2 statistics
                       default=True

           min       - list of minimums [min par0, min par1]
                       default=None

           max       - list of maximums [max par0, max par1]
                       default=None

           nloop     - list of bin sizes, used in calculating stepsize for each
                       dimension
                       default=(10,10)

           delv      - list of stepsizes, calculated by default
                       default=None

           fac       - factor used to expand or condense interval,
                       default=4

           log       - list of booleans to use log space for interval
                       default=(False,False)

           sigma     - list of sigmas used to calculate the confidence levels
                       (slices)
                       default=(1,2,3)

           levels    - confidence level values
                       default=None

           numcores  - specify the number of cores for parallel processing.
                       All available cores are used by default.
                       default=None

        Returns:
           reg_proj object

        DESCRIPTION

           Example: for users who do not want to create contours:

                print get_reg_proj( par0, par1, recalc=True )

        SEE ALSO
           int_unc, int_proj, reg_unc, get_reg_unc
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
    
    def get_reg_unc(self, par0=None, par1=None, id=None, otherids=None,
                    recalc=False, min=None, max=None, nloop=(10,10), delv=None,
                    fac=4, log=(False,False), sigma=(1,2,3), levels=None,
                    numcores=None):
        """
        get_reg_unc

        SYNOPSIS
           Return a confidence contour of fit statistic vs. two thawed
           parameter values.  At each step calculate the statistic with the
           other parameter(s) frozen at best fit values.

        SYNTAX

        Arguments:
           par0      - first source model parameter
                       default = None

           par1      - second source model parameter
                       default = None

           id        - Sherpa data id
                       default = default data id

           otherids  - list of ids required for simultaneous fit
                       default=None

           recalc    - calculate confidence data
                       default=False

           min       - list of minimums [min par0, min par1]
                       default=None

           max       - list of maximums [max par0, max par1]
                       default=None

           nloop     - list of bin sizes, used in calculating stepsize for each
                       dimension
                       default=(10,10)

           delv      - list of stepsizes, calculated by default
                       default=None

           fac       - factor used to expand or condense interval,
                       default=4

           log       - list of booleans to use log space for interval
                       default=(False,False)

           sigma     - list of sigmas used to calculate the confidence levels
                       (slices)
                       default=(1,2,3)

           levels    - confidence level values
                       default=None

           numcores  - specify the number of cores for parallel processing.
                       All available cores are used by default.
                       default=None

        Returns:
           reg_unc object 

        DESCRIPTION

           Example: for users who do not want to create contours:

              print get_reg_unc( par0, par1, recalc=True )

        SEE ALSO
           int_unc, int_proj, reg_proj, get_reg_proj
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

    
    def int_proj(self, par, id=None, otherids=None, replot=False, fast=True,
                 min=None, max=None, nloop=20, delv=None, fac=1, log=False,
                 numcores=None, overplot=False):
        """
        int_proj

        SYNOPSIS
           Create a confidence plot of fit statistic vs. a thawed parameter
           value.  At each step a fit is performed to obtain a new statistic
           if other thawed parameter(s) exist in the source model, otherwise,
           calculate the statistic (see int_unc).

        SYNTAX

        Arguments:
           par       - source model parameter

        Keyword arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - list of ids required for simultaneous fit
                       default=None

           replot    - replot the previously calculated data in cache
                       default=False

           fast      - change opt method to levmar for Chi2 statistics
                       default=True

           min       - minimum bound
                       default=None

           max       - maximum bound
                       default=None

           nloop     - bin size, used in calculating stepsize
                       default=20

           delv      - stepsize, calculated by default
                       default=None

           fac       - factor used to expand or condense interval,
                       default=1

           log       - boolean to use log space for interval
                       default=False

           numcores  - specify the number of cores for parallel processing.
                       All available cores are used by default.
                       default=None

           overplot  - plot over existing plot
                       default=False

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           int_unc, reg_proj, reg_unc
        """
        self._int_plot(self._intproj, par, id=id, otherids=otherids,
                       replot=replot, fast=fast, min=min, max=max, nloop=nloop,
                       delv=delv, fac=fac, log=log, numcores=numcores, 
                       overplot=overplot)

    def int_unc(self, par, id=None, otherids=None, replot=False, min=None,
                 max=None, nloop=20, delv=None, fac=1, log=False,
                 numcores=None, overplot=False):
        """
        int_unc

        SYNOPSIS
           Create a confidence plot of fit statistic vs. parameter value.  At
           each step calculate the statistic with the other parameter(s) frozen
           at best fit values.

        SYNTAX

        Arguments:
           par       - source model parameter

        Keyword arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - list of ids required for simultaneous fit
                       default=None

           replot    - replot the previously calculated data in cache
                       default=False

           min       - minimum bound
                       default=None

           max       - maximum bound
                       default=None

           nloop     - bin size, used in calculating stepsize
                       default=20

           delv      - stepsize, calculated by default
                       default=None

           fac       - factor used to expand or condense interval,
                       default=1

           log       - boolean to use log space for interval
                       default=False

           numcores  - specify the number of cores for parallel processing.
                       All available cores are used by default.
                       default=None

           overplot  - plot over existing plot
                       default=False

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           int_proj, reg_proj, reg_unc
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


    def reg_proj(self, par0, par1, id=None, otherids=None, replot=False,
                 fast=True, min=None, max=None, nloop=(10,10), delv=None, fac=4,
                 log=(False,False), sigma=(1,2,3), levels=None, numcores=None, 
                 overplot=False):
        """
        reg_proj

        SYNOPSIS
           Create a confidence contour of fit statistic vs. two thawed
           parameter values.  At each step a fit is performed to obtain a new
           statistic if other thawed parameter(s) exist in the source model,
           otherwise, calculate the statistic (see reg_unc).

        SYNTAX

        Arguments:
           par0      - first source model parameter
           par1      - second source model parameter

        Keyword arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - list of ids required for simultaneous fit
                       default=None

           replot    - replot the previously calculated data in cache
                       default=False

           fast      - change opt method to levmar for Chi2 statistics
                       default=True

           min       - list of minimums [min par0, min par1]
                       default=None

           max       - list of maximums [max par0, max par1]
                       default=None

           nloop     - list of bin sizes, used in calculating stepsize for each
                       dimension
                       default=(10,10)

           delv      - list of stepsizes, calculated by default
                       default=None

           fac       - factor used to expand or condense interval,
                       default=4

           log       - list of booleans to use log space for interval
                       default=(False,False)

           sigma     - list of sigmas used to calculate the confidence levels
                       (slices)
                       default=(1,2,3)

           levels    - confidence level values
                       default=None

           numcores  - specify the number of cores for parallel processing.
                       All available cores are used by default.
                       default=None

           overplot  - plot over existing plot
                       default=False

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           int_unc, int_proj, reg_unc
        """
        self._reg_plot(self._regproj, par0, par1, id=id, otherids=otherids,
                       replot=replot, fast=fast, min=min, max=max, nloop=nloop,
                       delv=delv, fac=fac, log=log, sigma=sigma, levels=levels,
                       numcores=numcores, overplot=overplot)

    def reg_unc(self, par0, par1, id=None, otherids=None, replot=False,
                min=None, max=None, nloop=(10,10), delv=None, fac=4,
                log=(False,False), sigma=(1,2,3), levels=None, numcores=None, 
                overplot=False):
        """
        reg_unc

        SYNOPSIS
           Create a confidence contour of fit statistic vs. two thawed
           parameter values.  At each step calculate the statistic with the
           other parameter(s) frozen at best fit values.

        SYNTAX

        Arguments:
           par0      - first source model parameter
           par1      - second source model parameter

        Keyword arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - list of ids required for simultaneous fit
                       default=None

           replot    - replot the previously calculated data in cache
                       default=False

           min       - list of minimums [min par0, min par1]
                       default=None

           max       - list of maximums [max par0, max par1]
                       default=None

           nloop     - list of bin sizes, used in calculating stepsize for each
                       dimension
                       default=(10,10)

           delv      - list of stepsizes, calculated by default
                       default=None

           fac       - factor used to expand or condense interval,
                       default=4

           log       - list of booleans to use log space for interval
                       default=(False,False)

           sigma     - list of sigmas used to calculate the confidence levels
                       (slices)
                       default=(1,2,3)

           levels    - confidence level values
                       default=None

           numcores  - specify the number of cores for parallel processing.
                       All available cores are used by default.
                       default=None

           overplot  - plot over existing plot
                       default=False

        Returns:
           None

        DESCRIPTION

        SEE ALSO
           int_unc, int_proj, reg_proj
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

    def get_data_image(self, id=None):
        """
        get_data_image

        SYNOPSIS
           Return a Sherpa data image obj

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa DataImage object

        DESCRIPTION
           The data image object holds the reference to the image array.

           Attributes:
              y            - image array

        SEE ALSO
           image_data
        """
        self._prepare_imageobj(id, self._dataimage)
        return self._dataimage
    
    def get_model_image(self, id=None):
        """
        get_model_image

        SYNOPSIS
           Return a Sherpa model image obj

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa ModelImage object

        DESCRIPTION
           The model image object holds the reference to the image array.

           Attributes:
              y            - image array

        SEE ALSO
           image_model
        """
        self._prepare_imageobj(id, self._modelimage)
        return self._modelimage

    def get_source_image(self, id=None):
        """
        get_source_image

        SYNOPSIS
           Return a Sherpa source image obj

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa SourceImage object

        DESCRIPTION
           The source image object holds the reference to the image array.

           Attributes:
              y            - image array

        SEE ALSO
           image_source
        """
        self._prepare_imageobj(id, self._sourceimage)
        return self._sourceimage


    def get_model_component_image(self, id, model=None):
        """
        get_model_component_image

        SYNOPSIS
           Return a Sherpa model image obj

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa ModelImage object

        DESCRIPTION
           The model image object holds the reference to the image array.

           Attributes:
              y            - image array

        SEE ALSO
           image_model
        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        self._prepare_imageobj(id, self._mdlcompimage, model=model)
        return self._mdlcompimage

    def get_source_component_image(self, id, model=None):
        """
        get_source_component_image

        SYNOPSIS
           Return a Sherpa model component image obj

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa ModelImage object

        DESCRIPTION
           The model component image object holds the reference to the 
           image array.

           Attributes:
              y            - image array

        SEE ALSO
           image_model
        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        self._prepare_imageobj(id, self._srccompimage, model=model)
        return self._srccompimage


    def get_ratio_image(self, id=None):
        """
        get_ratio_image

        SYNOPSIS
           Return a Sherpa ratio image obj

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa RatioImage object

        DESCRIPTION
           The ratio image object holds the reference to the image array.

           Attributes:
              y            - image array

        SEE ALSO
           image_ratio
        """
        self._prepare_imageobj(id, self._ratioimage)
        return self._ratioimage
    
    def get_resid_image(self, id=None):
        """
        get_resid_image

        SYNOPSIS
           Return a Sherpa resid image obj

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa ResidImage object

        DESCRIPTION
           The resid image object holds the reference to the image array.

           Attributes:
              y            - image array

        SEE ALSO
           image_resid
        """
        self._prepare_imageobj(id, self._residimage)
        return self._residimage

    def get_psf_image(self, id=None):
        """
        get_psf_image

        SYNOPSIS
           Return a Sherpa PSF image obj

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa PSFImage object

        DESCRIPTION
           The PSF image object holds the reference to the image array.

           Attributes:
              y            - image array

        SEE ALSO
           image_psf
        """
        self._prepare_imageobj(id, self._psfimage)
        return self._psfimage


    def get_kernel_image(self, id=None):
        """
        get_kernel_image

        SYNOPSIS
           Return a Sherpa PSF kernel image obj

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           Sherpa PSFKernelImage object

        DESCRIPTION
           The PSF kernel image object holds the reference to the image array.

           Attributes:
              y            - image array

        SEE ALSO
           image_kernel
        """
        self._prepare_imageobj(id, self._kernelimage)
        return self._kernelimage

    #
    # Images
    #

    def _image(self, id, imageobj, shape, newframe, tile, model=None):
        self._prepare_imageobj(id, imageobj, model).image(shape, newframe, tile)
        
    def image_data(self, id=None, newframe=False, tile=False):
        """
        image_data

        SYNOPSIS
           Send a data image to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           newframe    - Add a new frame
                         default = False

           tile        - Tile image frame
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize an image dataset by Sherpa data id.

        SEE ALSO
           get_data_image, image_model, image_fit, image_resid,
           image_ratio, image_fit_resid, image_psf
        """
        self._image(id, self._dataimage, None,
                    newframe, tile)
    
    def image_model(self, id=None, newframe=False, tile=False):
        """
        image_model

        SYNOPSIS
           Send a model image to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           newframe    - Add a new frame
                         default = False

           tile        - Tile image frame
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a model image by Sherpa data id.

        SEE ALSO
           get_model_image, image_data, image_fit, image_resid,
           image_ratio, image_fit_resid, image_psf, image_source
        """
        self._image(id, self._modelimage, None,
                    newframe, tile)


    def image_source_component(self, id, model=None, newframe=False,
                               tile=False):
        """
        image_source_component

        SYNOPSIS
           Send a source model component image to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           source       - Sherpa source component to image

           newframe    - Add a new frame
                         default = False

           tile        - Tile image frame
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a source model component image by Sherpa data id.

        SEE ALSO
           get_model_image, image_data, image_fit, image_resid,
           image_ratio, image_fit_resid, image_psf, image_source
        """
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        self._image(id, self._srccompimage, None, newframe, tile, model=model)


    def image_model_component(self, id, model=None, newframe=False, tile=False):
        """
        image_model_component

        SYNOPSIS
           Send a model component image to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           model       - Sherpa model component to image

           newframe    - Add a new frame
                         default = False

           tile        - Tile image frame
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a model component image by Sherpa data id.

        SEE ALSO
           get_model_image, image_data, image_fit, image_resid,
           image_ratio, image_fit_resid, image_psf, image_source
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


    def image_source(self, id=None, newframe=False, tile=False):
        """
        image_source

        SYNOPSIS
           Send a source image to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           newframe    - Add a new frame
                         default = False

           tile        - Tile image frame
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a source image by Sherpa data id.

        SEE ALSO
           get_source_image, image_data, image_fit, image_resid,
           image_ratio, image_fit_resid, image_psf, image_model
        """
        self._image(id, self._sourceimage, None, newframe, tile)

    def image_fit(self, id=None, newframe=True, tile=True, deleteframes=True):
        """
        image_fit

        SYNOPSIS
           Send a fit image to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           newframe    - Add a new frame
                         default = False

           tile        - Tile image frame
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a fit image by Sherpa data id.

        SEE ALSO
           get_fit_image, image_data, image_model, image_resid, image_ratio,
           image_fit_resid, image_psf
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

    def image_resid(self, id=None, newframe=False, tile=False):
        """
        image_resid

        SYNOPSIS
           Send a resid image to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           newframe    - Add a new frame
                         default = False

           tile        - Tile image frame
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a resid image by Sherpa data id.

        SEE ALSO
           get_resid_image, image_data, image_model, image_fit, image_ratio,
           image_fit_resid, image_psf
        """
        self._image(id, self._residimage, None,
                    newframe, tile)

    def image_ratio(self, id=None, newframe=False, tile=False):
        """
        image_ratio

        SYNOPSIS
           Send a ratio image to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           newframe    - Add a new frame
                         default = False

           tile        - Tile image frame
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a ratio image by Sherpa data id.

        SEE ALSO
           get_ratio_image, image_data, image_model, image_fit, image_resid,
           image_fit_resid, image_psf
        """
        self._image(id, self._ratioimage, None,
                    newframe, tile)

    def image_psf(self, id=None, newframe=False, tile=False):
        """
        image_psf

        SYNOPSIS
           Send a PSF image to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           newframe    - Add a new frame
                         default = False

           tile        - Tile image frame
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a PSF image by Sherpa data id.

        SEE ALSO
           get_psf_image, image_data, image_model,
           image_fit, image_resid, image_ratio, image_fit_resid
        """
        self._image(id, self._psfimage, None, newframe, tile)


    def image_kernel(self, id=None, newframe=False, tile=False):
        """
        image_kernel

        SYNOPSIS
           Send a PSF kernel image to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           newframe    - Add a new frame
                         default = False

           tile        - Tile image frame
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a PSF sub-kernel image by Sherpa data id.

        SEE ALSO
           get_kernel_image, image_data, image_model,
           image_fit, image_resid, image_ratio, image_fit_resid
        """
        self._image(id, self._kernelimage, None, newframe, tile)

    # Manage these functions (open, close, delete frames, regions, XPA)
    # through unbound functions of the Image class--always talking to
    # the same instance of the Image backend, so OK for now
    def image_deleteframes(self):
        """
        image_deleteframes

        SYNOPSIS
           Delete frames in image visualizer

        SYNTAX

        Arguments:
           None

        Returns:
           None

        DESCRIPTION
           Delete all frames in the image visualizer window.

        SEE ALSO
           image_open, image_close, image_getregion, image_setregion,
           image_xpaget, image_xpaset
        """
        sherpa.image.Image.delete_frames()
        
    def image_open(self):
        """
        image_open

        SYNOPSIS
           Open a image visualizer window

        SYNTAX

        Arguments:
           None

        Returns:
           None

        DESCRIPTION
           Open a window for image visualization.

        SEE ALSO
           image_deleteframes, image_close, image_getregion, image_setregion,
           image_xpaget, image_xpaset
        """
        sherpa.image.Image.open()

    def image_close(self):
        """
        image_close

        SYNOPSIS
           Close a image visualizer window

        SYNTAX

        Arguments:
           None

        Returns:
           None

        DESCRIPTION
           Close a image visualizer window.

        SEE ALSO
           image_open, image_deleteframes, image_getregion, image_setregion,
           image_xpaget, image_xpaset
        """
        sherpa.image.Image.close()

    def image_getregion(self, coord=''):
        """
        image_getregion

        SYNOPSIS
           Return a visualizer image region

        SYNTAX

        Arguments:
           None

        Returns:
           visualizer image region

        DESCRIPTION
           Return a visualizer image region.

        SEE ALSO
           image_open, image_close, image_deleteframes, image_setregion,
           image_xpaget, image_xpaset
        """
        return sherpa.image.Image.get_region(coord)

    def image_setregion(self, reg, coord=''):
        """
        image_setregion

        SYNOPSIS
           Set a visualizer image region

        SYNTAX

        Arguments:
           reg      - image visualizer region

        Returns:
           None

        DESCRIPTION
           Set a visualizer image region.

        SEE ALSO
           image_open, image_close, image_getregion, image_deleteframes,
           image_xpaget, image_xpaset
        """
        sherpa.image.Image.set_region(reg, coord)

    def image_xpaget(self, arg):
        """
        image_xpaget

        SYNOPSIS
           Return an XPA data stream

        SYNTAX

        Arguments:
           arg      - XPA agrument

        Returns:
           XPA data stream

        DESCRIPTION
           Retrieve an XPA data stream from image visualizer.

        SEE ALSO
           image_open, image_close, image_getregion, image_setregion,
           image_deleteframes, image_xpaset
        """
        return sherpa.image.Image.xpaget(arg)

    def image_xpaset(self, arg, data=None):
        """
        image_deleteframes

        SYNOPSIS
           Send an XPA data stream

        SYNTAX

        Arguments:
           arg      - XPA agrument

           data     - data to be sent as stdin to the
                      XPA argument (None by default)
        Returns:
           None

        DESCRIPTION
           Send an XPA data stream to image visualizer.

        SEE ALSO
           image_open, image_close, image_getregion, image_setregion,
           image_xpaget, image_deleteframes
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

