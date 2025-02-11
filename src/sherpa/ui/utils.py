#
#  Copyright (C) 2010, 2015 - 2024
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

from __future__ import annotations

from configparser import ConfigParser
import copy
import copyreg as copy_reg
from dataclasses import dataclass
import importlib
import logging
import os
import pickle
import sys
from typing import Any, Callable, Literal, Optional, Sequence, \
    TypeVar, Union, overload

import numpy as np

from sherpa import get_config
import sherpa.all
from sherpa.data import Data, DataSimulFit
from sherpa.estmethods import EstMethod
from sherpa.fit import Fit, FitResults
from sherpa.models.basic import TableModel
from sherpa.models.model import Model, SimulFitModel
from sherpa.models.template import add_interpolator, create_template_model, \
    reset_interpolators
from sherpa.optmethods import OptMethod
from sherpa.plot import Plot, MultiPlot, set_backend, get_per_plot_kwargs
from sherpa.stats import Stat, UserStat
from sherpa.utils import NoNewAttributesAfterInit, is_subclass, \
    export_method, send_to_pager
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, \
    DataErr, IdentifierErr, IOErr, ModelErr, ParameterErr, PlotErr, \
    SessionErr
from sherpa.utils.numeric_types import SherpaFloat
from sherpa.utils.random import RandomType
from sherpa.utils.types import IdType

info = logging.getLogger(__name__).info
warning = logging.getLogger(__name__).warning

config = ConfigParser()
config.read(get_config())

np.set_printoptions(threshold=int(config.get('verbosity',
                                             'arraylength',
                                             fallback=1000000)))

__all__ = ('ModelWrapper', 'Session')

BUILTINS = sys.modules["builtins"]
_builtin_symbols_ = tuple(BUILTINS.__dict__.keys())


ModelType = Union[Model, str]
T = TypeVar("T")


###############################################################################
#
# Errors and argument checking
#
###############################################################################


def _check_type(arg, argtype, argname: str, argdesc: str) -> None:
    if isinstance(arg, argtype):
        return

    raise ArgumentTypeErr('badarg', argname, argdesc)


def _check_str_type(arg: str, argname: str) -> None:
    """Ensure that arg (with name argname) is a string"""
    if _is_str(arg):
        return

    raise ArgumentTypeErr('badarg', argname, "a string")


# With Python 3.10 these can return TypeGuard annotations.
#
def _is_integer(val):
    return isinstance(val, (int, np.integer))


def _is_str(val):
    return isinstance(val, (str, ))


def get_plot_prefs(plotobj):
    """Return the preferences for the plot object.

    The current preference design has the attribute name be different
    depending if this is a "histogram" style compared to a "line"
    style, so wrap this logic up. If there is no matching preference
    then raise an AttributeError. If the object contains both line and
    histogram preferences then the histogram preferences are returned.

    """

    try:
        return plotobj.histo_prefs
    except AttributeError:
        try:
            return plotobj.plot_prefs
        except AttributeError:
            raise AttributeError("plot object has no preferences") from None


def _get_filter(data):
    """Report the filter for report_filter_change.

    This takes care of calling get_filter with the correct
    arguments (needed in case anyone wants to use the string
    in a call to notice/ignore), and it handles the case of
    the case of all-data-being-filtered leading to an error
    for certain data types.

    Parameters
    ----------
    data : sherpa.data.Data1D instance

    Returns
    -------
    msg : str or None
        The filter expression, with "all-data-filtered" indicated by
        the empty string, or None when there was an error (which
        shouldn't happen, but we don't want to fail a notice or ignore
        call when reporting the change).

    """

    # See #1430 for a discussion on why we have different behavior
    # when all the data is excluded: it should be unified.
    #
    # We follow the approach used for the DataIMG case be relying on
    # the mask attribute to determine what the output should
    # be. Unfortunately this is only helpful for the "all data
    # excluded" case, as we do not have the equivalent of "Field()" to
    # indicate "use all the data", and we can not remove the DataErr
    # exception check from get_filter, just in case.
    #
    if data.mask is False:
        return ""

    try:
        return data.get_filter(delim=':', format='%g')
    except DataErr as exc:
        # It may be possible that this happens - e.g. that data.mask
        # is a ndarray of all False values - but it is not 100%
        # obvious.  This is left in just in case.
        #
        if str(exc) == "mask excludes all data":
            return ""

        # It's not clear if this can happen, but assume it can.
        #
        return None

    except TypeError:
        # This can happen when the dataset is Data2D, as get_filter
        # does not accept kwargs, so we skip over it (as the handling
        # of filters for Data2D is not clear).
        #
        return None

    except IndexError:
        # This is a known failure case with ignore_bad handling of
        # DataPHA. Is this still an issue? It is safer to leave in.
        #
        return None


def report_filter_change(idstr: str,
                         ofilter: Optional[str],
                         nfilter: Optional[str],
                         xlabel: Optional[str] = None
                         ) -> None:
    """Report the filter change for ignore/filter.

    Parameters
    ----------
    idstr : str
       The dataset identifier (expected to be "dataset <idval>" but
       may be "dataset <idval>: background <bkg_id>" for DataPHA).
    ofilter, nfilter : str or None
       The filter string before and after filtering (the output
       of _get_filter).
    xlabel : str or None
       The units of the filter (if set).

    Notes
    -----

    A filter expression of "" (the empty string) is taken to mean all
    data has been removed, and converted to something more readable
    here.

    If either ofilter or nfilter is None then we assume that some
    error has happened. The reporting of this depends on what has gone
    wrong. That is, are both None (so broken before and after), is
    only nfilter None (so now broken), or is only ofilter None (so it
    has somehow been "fixed")?

    """

    ostr = f"{idstr}: "

    if ofilter is None and nfilter is not None:
        # Let the user know the filter has changed but do not create a
        # warning message. Chose "<broken>" as the original filter
        # state since this is not going to match nfilter, and
        # indicates that the user should be careful.
        #
        ofilter = "<broken>"

    if ofilter is None or nfilter is None:
        emsg = f"{ostr}1D filter has failed"
        logging.getLogger(__name__).error(emsg)
        return

    # Make it easy to handle labels being optional
    if xlabel is None:
        label = ""
    else:
        label = f" {xlabel}"

    # We have to be careful because we can get a filter expression
    # of '' when there's no data. The current handling of this case
    # is not robust (sometimes it errors out, sometimes it is an
    # empty string).
    #
    if ofilter == "":
        if nfilter == "":
            ostr += "no data (unchanged)"

        else:
            ostr += f"no data -> {nfilter}{label}"

    elif nfilter == "":
        ostr += f"{ofilter}{label} -> no data"

    else:
        ostr += f"{ofilter}"
        if ofilter == nfilter:
            ostr += f"{label} (unchanged)"
        else:
            ostr += f" -> {nfilter}{label}"

    info(ostr)


def notice_data_range(get_data: Callable[[IdType], Data],
                      ids: Sequence[IdType],
                      lo, hi,
                      kwargs) -> None:
    """Filter each dataset and report the change in filter.

    Parameters
    ----------
    get_data : callable
        The routine to use to return the data object given a dataset
        identifier.
    ids : sequence of int or str
        The identifiers to process, in order. It is required that
        get_data can be called on each item in this sequence.
    lo, hi : number or None
        The arguments for `notice`.
    kwargs
        The extra arguments to pass to the Data `notice`, and
        the "bkg_id" identifier if the data to be filtered is a
        background PHA component instead.

    Notes
    -----
    This could be part of the Session class, but we don't need more
    methods there, so see how well things work with it out here. If it
    were a class method we could move the handling of the bkg_id
    argument to the astro Session class.

    """

    # The bkg_id argument is astro-specific. The need
    # to use it to access the correct data object is
    # problematic.
    #
    bkg_id = kwargs.pop("bkg_id", None)

    for idval in ids:
        idstr = f"dataset {idval}"
        data = get_data(idval)
        if bkg_id is not None:
            data = data.get_background(bkg_id)
            idstr += f": background {bkg_id}"

        ofilter = _get_filter(data)
        data.notice(lo, hi, **kwargs)
        nfilter = _get_filter(data)

        try:
            xlabel = data.get_xlabel()
        except AttributeError:
            # Data2D case
            xlabel = None

        report_filter_change(idstr, ofilter, nfilter, xlabel)


@overload
def calc_multiplot_size(rows: int,
                        cols: Optional[int],
                        nplots: int
                        ) -> tuple[int, int]:
    ...

@overload
def calc_multiplot_size(rows: Optional[int],
                        cols: int,
                        nplots: int
                        ) -> tuple[int, int]:
    ...

def calc_multiplot_size(rows, cols, nplots):
    """How many rows and columns to use for multi plot/contours?

    Parameters
    ----------
    rows, cols
       The requested sizes, if set (one of them must be set).
    nplots
       The number of plots.

    Returns
    -------
    (nrows, ncols)

    """

    def get(n: int) -> int:
        return int(np.ceil(nplots / n))

    # Since at least one of rows or cols must be set.
    #
    nrows = rows if rows is not None else get(cols)
    ncols = cols if cols is not None else get(rows)

    if nrows * ncols < nplots:
        # Go for the nearest (upper) square value for the number
        # of columns, and calculate the number of rows to match
        # the data.
        #
        ncols = int(np.ceil(np.sqrt(nplots)))
        nrows = get(ncols)

    return nrows, ncols


###############################################################################
#
# Pickling support
#
###############################################################################


def construct_ufunc(modname, funcname):
    module = sys.modules.get(modname)
    if module is None:
        module = importlib.import_module(modname)
    return getattr(module, funcname)


def reduce_ufunc(func):
    modname = getattr(func, '__module__', 'numpy')
    funcname = func.__name__
    if func is not getattr(sys.modules[modname], funcname, None):
        raise ValueError("module '%s' does not contain ufunc '%s'" %
                         (modname, funcname))
    return (construct_ufunc, (modname, funcname))


copy_reg.constructor(construct_ufunc)
copy_reg.pickle(np.ufunc, reduce_ufunc)


###############################################################################
#
# I/O routines
#
###############################################################################


def read_template_model(modelname, templatefile,
                        sep=' ', comment='#',
                        method=sherpa.utils.linear_interp,
                        template_interpolator_name='default'):
    """Read in a set of templates and create a template model.

    Parameters
    ----------
    modelname : str
       The identifier for this table model.
    templatefile : str
       The name of the file to read in. This file lists the template
       data files.
    sep : str, optional
       The separator character. The default is ``' '``.
    comment : str, optional
       The comment character. The default is ``'#'``.
    method : func
       The interpolation method to use to map the input data onto the
       coordinate grid of the data set. Linear, nearest-neighbor, and
       polynomial schemes are provided in the sherpa.utils module.
    template_interpolator_name : str
       The method used to interpolate within the set of templates.
       The default is ``default``. A value of ``None`` turns off the
       interpolation; in this case the grid-search optimiser must be
       used to fit the data.

    Returns
    -------
    model
        The template model.

    """

    if sherpa.utils.is_binary_file(templatefile):
        raise IOErr("notascii", templatefile)

    names, cols = sherpa.io.read_file_data(templatefile, sep=sep,
                                           comment=comment,
                                           require_floats=False)

    ncols = len(cols)
    nnames = len(names)
    if ncols != nnames:
        raise IOErr("wrongnumcols", ncols, nnames)

    names = [name.strip().lower() for name in names]

    filenames = None
    modelflags = None
    parnames = names[:]
    parvals = []
    for name, col in zip(names, cols):
        # Find the column with the filenames, remove it from the set of
        # parameter names
        if name.startswith("file"):
            filenames = col
            parnames.remove(name)
            continue

        # Find the column with the modelflags, remove it from the set of
        # parameter names
        if name.startswith("modelflag"):
            modelflags = col
            parnames.remove(name)
            continue

        parvals.append(np.array(col, dtype=SherpaFloat))

    parvals = np.asarray(parvals).T

    if len(parvals) == 0:
        raise IOErr("noparamcols", templatefile)

    if filenames is None:
        raise IOErr("reqcol", "filename", templatefile)

    if modelflags is None:
        raise IOErr("reqcol", "modelflag", templatefile)

    templates = []
    for filename in filenames:
        tnames, tcols = sherpa.io.read_file_data(filename, sep=sep, comment=comment)

        ntcols = len(tcols)
        if ntcols == 1:
            raise IOErr("onecolneedtwo", filename)
        if ntcols != 2:
            raise IOErr("wrongnumcols", 2, ntcols)

        tm = TableModel(filename)
        tm.method = method  # interpolation method
        tm.load(*tcols)
        tm.ampl.freeze()
        templates.append(tm)

    assert len(templates) == parvals.shape[0]

    return create_template_model(modelname, parnames, parvals,
                                 templates,
                                 template_interpolator_name=template_interpolator_name)


###############################################################################
#
# ModelWrapper
#
###############################################################################


class ModelWrapper(NoNewAttributesAfterInit):
    """Wrap up a model class so we can create instances easily.

    Creates a model instance with a given name - you can say mdl.mname
    or mdl("mname") - and either syntax stores the model instance in
    the session object with the name "mname", over-writing any
    previous component with this name.  If mdl1 and mdl2 have been
    created by ModelWrapper then a user can say mdl1.n1 + mdl2.n2 to
    create model instances named n1 and n2 and then return the
    expression which represents their sum.

    """

    # How do we type modeltype so that we can label __call__ as
    # returning an instance of it, and still identify as deriving from
    # Model?
    #
    def __init__(self,
                 session: Session,
                 modeltype,
                 args=(), kwargs=None) -> None:
        # This is an internal class so do not bother with
        # sherpa.utils.err exceptions.
        #
        if not isinstance(session, Session):
            raise ValueError(f"session={session} is not a Session instance")

        if not is_subclass(modeltype, Model):
            raise ValueError(f"modeltype={modeltype} is not a Model class")

        self._session = session
        self.modeltype = modeltype
        self.args = args
        self.kwargs = kwargs if kwargs else {}

        # Edit the docstring of the new object to hide the
        # ModelWrapper text and replace it with a reference to the
        # wrapped class.
        #
        mname = modeltype.__name__.lower()
        prefix = "an" if mname[0] in "aeiou" else "a"

        # Grab the first line of the model description. It might be
        # helpful to also list the parameters but this is not easy to
        # access without actually creating an instance of the model.
        #
        if modeltype.__doc__ is None:
            mdesc = ""
        else:
            # Could try to ensure this is a sentence - e.g. ends in a
            # "." - but rely on the documentation of the model classes
            # to enforce that.
            #
            mdesc = modeltype.__doc__.split("\n")[0]
            mdesc += "\n\n    "

        self.__doc__ = f"""Create {prefix} {mname} model instance.

    {mdesc}Instances can be created either as an attribute of {mname},
    as long as the attribute does not begin with an underscore,
    or by calling {mname} directly.

    Examples
    --------

    The model, here called mdl, is returned but it's also stored in
    the session and can be returned with `get_model_component`:

    >>> m1 = {mname}.mdl

    If the model has already been created with the same name then
    the old version will be returned, rather than creating a new
    instance:

    >>> m1 = {mname}.mdl
    >>> m2 = {mname}("mdl")
    >>> m1 == m2
    True

"""

        NoNewAttributesAfterInit.__init__(self)

    # TODO: can we say that this is the same class as sent to
    # the constructor? A TypeVar could be used, but how do we
    # constrain it to Model?
    #
    def __call__(self, name):
        _check_str_type(name, "name")

        m = self._session._get_model_component(name)
        if isinstance(m, self.modeltype):
            return m

        m = self.modeltype(name, *self.args, **self.kwargs)
        self._session._add_model_component(m)
        return m

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

        return self(name)

    def __repr__(self) -> str:
        return f'<{self.modeltype.__name__} model type>'

    def __str__(self) -> str:
        if self.modeltype.__doc__ is not None:
            # Use the documentation from wrapped model if available,
            # rather than a normal "repr" call.
            #
            return self.modeltype.__doc__

        return self.__repr__()


def _assign_obj_to_main(name: str, obj) -> None:
    """Create a "global" symbol.

    Parameters
    ----------
    name : str
        The name to use.
    obj : object
        The value for the symbol.

    """
    sys.modules["__main__"].__dict__[name] = obj
    BUILTINS.__dict__[name] = obj


def _assign_model_to_main(name: str, model: Model) -> None:
    """Ensure the model is added to the "global" symbol table.

    Parameters
    ----------
    name : str
        The name to use.
    model : sherpa.models.model.Model instance
        The model. The name field will be changed to match the
        "<instance>.<name>" scheme.

    """
    # Ask sys what the __main__ module is; packages such
    # as IPython can add their own __main__ module.
    model.name = '%s.%s' % (type(model).__name__.lower(), name)
    _assign_obj_to_main(name, model)


def _remove_obj_from_main(name: str) -> None:
    """Remove the symbol."""

    # Is this sufficient or overkill?
    #
    for base in [BUILTINS, sys.modules["__main__"]]:
        try:
            del base.__dict__[name]
        except KeyError:
            pass


def set_dep(data, val) -> None:
    """Set the dependent axis.

    Parameters
    ----------
    data : sherpa.data.Data instance
    val : scalar or array_like
        If a scalar then use the same value for all elements.

    """

    if np.iterable(val):
        dep = np.asarray(val, SherpaFloat)
    else:
        val = SherpaFloat(val)
        dep = np.array([val] * data.size)

    data.dep = dep


def set_error(data, field, val, fractional=False) -> None:
    """Set the error field.

    Parameters
    ----------
    data : sherpa.data.Data instance
    field : {"staterror", "syserror"}
    val : None, a scalar, or array_like
        The statistical error.
    fractional : bool, optional
        If `val` is a scalar then set errors to
        `val * data.y`.

    Notes
    -----
    The `val` field, when an array, must match the data.size
    field, even for grouped PHA data.

    """

    if val is None:
        err = None

    elif np.iterable(val):
        err = np.asarray(val, SherpaFloat)

    else:
        val = SherpaFloat(val)
        if sherpa.utils.bool_cast(fractional):
            err = val * data.get_dep()
        else:
            err = np.array([val] * data.size)

    setattr(data, field, err)


def set_filter(data, val, ignore=False) -> None:
    """Set the filter field.

    Parameters
    ----------
    data : sherpa.data.Data instance
    val : array_like
        The filter array (bool-like).
    ignore : bool, optional
        If set the filter should be inverted.

    Notes
    -----
    The `val` field should match the mask array, that is the
    "effective" size of the data. This is only different from
    data.size when the data is grouped (and at least one group
    contains multiple channels).

    """

    val = np.asarray(val, dtype=np.bool_)
    nval = len(val)

    # Note that we do not use data.size as the check, as that fails
    # for grouped PHA data. Is there some way to have an "effective
    # size" value?
    #
    nexp = len(data.get_y(False))
    if np.iterable(data.mask):
        if nexp != nval:
            raise sherpa.utils.err.DataErr('mismatchn', 'data', 'filter',
                                           nexp, nval)

        if not ignore:
            data.mask |= val
        else:
            data.mask &= ~val

        return

    # The mask attribute checks the length, so we do not have to here.
    #
    if not ignore:
        data.mask = val
    else:
        data.mask = ~val


@dataclass
class FitStore:
    """Store per-dataset information for a fit"""

    idval : IdType
    data : Data
    model : Model


def get_components_helper(getfunc: Callable[..., Plot],
                          model: Model,
                          idval: IdType) -> MultiPlot:
    """Handle get_source/model_components_plot.

    Iterate through each term in the source model.
    """

    out = MultiPlot()
    cpts = sherpa.models.model.model_deconstruct(model)
    for cpt in cpts:
        plotobj = getfunc(id=idval, model=cpt, recalc=True)
        out.add(plotobj)

    out.title = "Component plot"
    return out


###############################################################################
#
# Session
#
###############################################################################
class Session(NoNewAttributesAfterInit):

    ###########################################################################
    # Standard methods
    ###########################################################################

    def __init__(self) -> None:
        self.clean()
        self._model_types: dict[str, ModelWrapper] = {}
        self._model_globals = np.__dict__.copy()
        NoNewAttributesAfterInit.__init__(self)
        global _session
        _session = self

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_model_globals']
        return state

    def __setstate__(self, state):
        self._model_globals = np.__dict__.copy()

        self._model_globals.update(state['_model_types'])

        if '_sources' not in state:
            self.__dict__['_sources'] = state.pop('_models')

        self.__dict__.update(state)

    ###########################################################################
    # High-level utilities
    ###########################################################################

    def _export_names(self, gdict):
        allnames = []

        for name in dir(self):
            if ((not name.startswith('_')) or
                    (name == '_sherpa_version') or
                    (name == '_sherpa_version_string')):
                gdict[name] = export_method(getattr(self, name),
                                            modname=gdict.get('__name__'))
                allnames.append(name)

        gdict.update(self._model_types)
        allnames.extend(self._model_types.keys())

        return allnames

    def clean(self) -> None:
        """Clear out the current Sherpa session.

        The `clean` function removes all data sets and model
        assignments, and restores the default settings for the
        optimisation and fit statistic.

        .. versionchanged:: 4.15.0
           The model names are now removed from the global symbol
           table.

        See Also
        --------
        save : Save the current Sherpa session to a file.
        restore : Load in a Sherpa session from a file.
        sherpa.astro.ui.save_all : Save the Sherpa session as an ASCII file.

        Examples
        --------

        >>> clean()

        After the call to `clean`, the `line` and `bgnd` variables
        will be removed, so accessing them would cause a NameError.

        >>> set_source(gauss1d.line + const1d.bgnd)
        >>> bgnd.c0.min = 0
        >>> print(line)
        >>> clean()

        """

        # Best-guess attempt at removing all the symbols we have
        # created (e.g. from creating model instances). We only do
        # this when the _assign_model_to_main routine is in place:
        # this is not ideal, as it could have been changed after some
        # symbols were created, but assume that if the user is going
        # to do this then they can deal with the consequences.
        #
        # It looks like _model_components stores all the various models,
        # such as table and psf models.
        #
        try:
            if self._model_autoassign_func == _assign_model_to_main:
                for name in self._model_components.keys():
                    _remove_obj_from_main(name)

        except AttributeError:
            pass

        self._sherpa_version = sherpa.__version__
        self._sherpa_version_string = sherpa.__version__

        self._default_id: IdType = 1
        self._paramprompt = False

        self._methods: dict[str, OptMethod] = {}
        self._itermethods = {'none': {'name': 'none'},
                             'sigmarej': {'name': 'sigmarej',
                                          'maxiters': 5,
                                          'hrej': 3,
                                          'lrej': 3,
                                          'grow': 0}}

        self._stats: dict[str, Stat] = {}
        self._estmethods: dict[str, EstMethod] = {}

        modules = (sherpa.optmethods, sherpa.stats, sherpa.estmethods)
        basetypes = (OptMethod, Stat, EstMethod)
        objdicts = (self._methods, self._stats, self._estmethods)

        for mod, base, odict in zip(modules, basetypes, objdicts):
            for name in mod.__all__:
                cls = getattr(mod, name)
                if is_subclass(cls, base):
                    odict[name.lower()] = cls()

        # Note: levmar does not support the rng option so this
        # can be done before set_rng is called.
        #
        self._current_method = self._methods['levmar']
        self._current_itermethod = self._itermethods['none']
        self._current_stat = self._stats['chi2gehrels']
        # Add simplex as alias to neldermead
        self._methods['simplex'] = self._methods['neldermead']

        reset_interpolators()

        # Should some of these dictionaries have more-restrictive
        # types for the value? The current hierarchy doesn't help, as
        # the "base" class is often a "virtual" class which ends up
        # missing some of the methods that the actual instances rely
        # on, or that it has different parameters than its
        # sub-classes (e.g. the plot classes), or that maybe a
        # more-specialized class than Model would be helpful
        # (e.g. _tbl_models).
        #
        self._data: dict[IdType, Data] = {}
        self._psf: dict[IdType, Model] = {}
        self._tbl_models: list[Model] = []
        self._psf_models: list[Model] = []

        self._model_autoassign_func = _assign_model_to_main
        self._model_components: dict[str, Model] = {}
        self._models: dict[IdType, Model] = {}
        self._sources: dict[IdType, Model] = {}

        self._fit_results = None
        self._pvalue_results = None

        self._covariance_results = None
        self._confidence_results = None
        self._projection_results = None

        self._pyblocxs = sherpa.sim.MCMC()

        self._splitplot = sherpa.plot.SplitPlot()
        self._jointplot = sherpa.plot.JointPlot()

        # self._comptmplmdlplot = sherpa.plot.ComponentTemplateModelPlot()
        self._comptmplsrcplot = sherpa.plot.ComponentTemplateSourcePlot()

        self._lrplot = sherpa.plot.LRHistogram()
        self._pdfplot = sherpa.plot.PDFPlot()
        self._cdfplot = sherpa.plot.CDFPlot()
        self._traceplot = sherpa.plot.TracePlot()
        self._scatterplot = sherpa.plot.ScatterPlot()

        self._intproj = sherpa.plot.IntervalProjection()
        self._intunc = sherpa.plot.IntervalUncertainty()
        self._regproj = sherpa.plot.RegionProjection()
        self._regunc = sherpa.plot.RegionUncertainty()

        self._set_plot_types()
        self._set_contour_types()
        self._set_image_types()

        # Reset the generator.
        #
        self.set_rng(None)

    def _set_plot_types(self) -> None:
        """Set up the plot types."""

        # The keys of _plot_types are used to define:
        # a) the mapping from get_<key>_plot to the plot objects
        #    (that is, there must be a matching get_<key>_plot method)
        # b) valid arguments for the plot() call
        # c) the arguments that set_xlog/... accept
        # d) a set of names that can not be used as a dataset identifier
        #    (because of point b); see also _contour_types
        #
        # Note that not all plot_xxx commands use this structure.
        #
        # Unlike the contour case, we have different plot classes to
        # handle different types of plot. There are either plot types
        # with only one plot class, such as "fit", or those that
        # support a version that depends on the data type, such as
        # "data". A list is used in both cases, with the first element
        # being the "generic" case and, for those that support it, the
        # second being the Data1DInt plot class.
        #
        # There is an argument to be made to have the "singleton"
        # plots, such as "fit", storing the object directly, and not
        # as a single-element list, but it was felt that having a
        # consistent access pattern was cleaner.
        #
        self._plot_types: dict[str, list[Any]] = {
            'data': [sherpa.plot.DataPlot(), sherpa.plot.DataHistogramPlot()],
            'model': [sherpa.plot.ModelPlot(), sherpa.plot.ModelHistogramPlot()],
            'model_component': [sherpa.plot.ComponentModelPlot(), sherpa.plot.ComponentModelHistogramPlot()],
            'model_components': [sherpa.plot.MultiPlot()],
            'source': [sherpa.plot.SourcePlot(), sherpa.plot.SourceHistogramPlot()],
            'source_component': [sherpa.plot.ComponentSourcePlot(), sherpa.plot.ComponentSourceHistogramPlot()],
            'source_components': [sherpa.plot.MultiPlot()],
            'fit': [sherpa.plot.FitPlot()],
            'resid': [sherpa.plot.ResidPlot(), sherpa.plot.ResidHistogramPlot()],
            'ratio': [sherpa.plot.RatioPlot(), sherpa.plot.RatioHistogramPlot()],
            'delchi': [sherpa.plot.DelchiPlot(), sherpa.plot.DelchiHistogramPlot()],
            'chisqr': [sherpa.plot.ChisqrPlot(), sherpa.plot.ChisqrHistogramPlot()],
            'psf': [sherpa.plot.PSFPlot()],
            'kernel': [sherpa.plot.PSFKernelPlot()]
        }

        # Set up aliases so that calls to set_xlog/.. will still
        # succeed, but the user will be told to use the new
        # name. These keys are also used, along with _plot_type_names,
        # to determine the set of forbidden identifiers.
        #
        self._plot_types_alias: dict[str, str] = {
            "compsource": "source_component",
            "compmodel": "model_component"
        }

    def _set_contour_types(self) -> None:
        """Set up the contour types."""

        # This is used by the get_<key>_contour calls to access the
        # relevant contour class. The keys define the labels that can be
        # used in calls to contour(), and are also used to determine
        # the set of forbidden identifiers.
        #
        self._contour_types: dict[str, Any] = {
            "data": sherpa.plot.DataContour(),
            "model": sherpa.plot.ModelContour(),
            "source": sherpa.plot.SourceContour(),
            "fit": sherpa.plot.FitContour(),
            "resid": sherpa.plot.ResidContour(),
            "ratio": sherpa.plot.RatioContour(),
            "psf": sherpa.plot.PSFContour(),
            "kernel": sherpa.plot.PSFKernelContour()
        }

    def _set_image_types(self) -> None:
        """Set up the image types."""

        # This is used by the get_<key>_image calls to access the
        # relevant image class. The keys are not included in any check
        # of valid identifiers. unlike the plot and contour cases, as
        # there is
        #
        # - no image() call that acts like plot() or contour();
        # - and no set_xlog/... like call to change the image displays.
        #
        self._image_types: dict[str, Any] = {
            'data': sherpa.image.DataImage(),
            'model': sherpa.image.ModelImage(),
            'source': sherpa.image.SourceImage(),
            'ratio': sherpa.image.RatioImage(),
            'resid': sherpa.image.ResidImage(),
            'psf': sherpa.image.PSFImage(),
            'kernel': sherpa.image.PSFKernelImage(),
            'model_component': sherpa.image.ComponentModelImage(),
            'source_component': sherpa.image.ComponentSourceImage()
        }

    def get_rng(self) -> Optional[RandomType]:
        """Return the RNG generator in use.

        The return can be None, which means that the routines in
        `numpy.random` are used, and thus are affected by calls to
        `numpy.random.seed`, otherwise the supplied generator is used
        to create random numbers. See
        https://numpy.org/doc/stable/reference/random/legacy.html for
        more information.

        .. versionadded:: 4.16.0

        See Also
        --------
        set_rng

        """

        return self._rng

    def set_rng(self, rng: Optional[RandomType]) -> None:
        """Set the RNG generator.

        .. versionadded:: 4.16.0
           This replaces the seed argument for certain routines and
           the need to explicitly call `numpy.random.seed` in others.

        Parameters
        ----------
        rng : numpy.random.Generator, numpy.random.RandomState, or None
           Determines how random numbers are created. If set to None
           then the routines in `numpy.random` are used, and so can be
           controlled by calling `numpy.random.seed`.

        See Also
        --------
        get_rng

        """

        if rng is not None and not isinstance(rng,
                                              (np.random.Generator,
                                               np.random.RandomState)):
            # Do not include RandomState in the error message as it is
            # really meant for testing/old code.
            raise ArgumentTypeErr("badarg", "rng", "a Generator or None")

        self._rng = rng

    def set_plot_backend(self, backend) -> None:
        """Change the plot backend.

        This will reset any plot structures, such as that returned by
        get_data_plot.

        .. versionadded:: 4.16.0

        Parameters
        ----------
        backend : str
            The name of the plot backend.

        """

        set_backend(backend)

        # Re-create all the plot objects
        self._set_plot_types()
        self._set_contour_types()

    def save(self, filename='sherpa.save', clobber=False) -> None:
        """Save the current Sherpa session to a file.

        Parameters
        ----------
        filename : str, optional
           The name of the file to write the results to. The default
           is 'sherpa.save'.
        clobber : bool, optional
           This flag controls whether an existing file can be
           overwritten (``True``) or if it raises an exception (``False``,
           the default setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

        See Also
        --------
        clean : Clear all stored session data.
        restore : Load in a Sherpa session from a file.
        sherpa.astro.ui.save_all : Save the Sherpa session as an ASCII file.

        Notes
        -----
        The current Sherpa session is saved using the Python `pickle`
        module. The output is a binary file, which may not be portable
        between versions of Sherpa, but is platform independent, and
        contains all the data. This means that files created by `save`
        can be sent to collaborators to share results.

        Examples
        --------

        Save the current session to the file 'sherpa.save'.

        >>> save()

        Save the current session to the file 'bestfit.sherpa',
        overwriting any existing version of the file.

        >>> save('bestfit.sherpa', clobber=True)

        """

        _check_str_type(filename, "filename")
        clobber = sherpa.utils.bool_cast(clobber)

        if os.path.isfile(filename) and not clobber:
            raise IOErr("filefound", filename)

        fout = open(filename, 'wb')
        try:
            # Use the default version rather than fix a version
            pickle.dump(self, fout)
        finally:
            fout.close()

    def restore(self, filename='sherpa.save') -> None:
        """Load in a Sherpa session from a file.

        .. warning::
             Security risk: The imported functions and objects
             could contain arbitrary Python code and be malicious.
             Never use this function on untrusted input.

        Parameters
        ----------
        filename : str, optional
           The name of the file to read the results from. The default
           is 'sherpa.save'.

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
        a message is::

          WARNING: Could not determine whether the model is discrete.
          This probably means that you have restored a session saved with a previous version of Sherpa.
          Falling back to assuming that the model is continuous.

        Examples
        --------

        Load in the Sherpa session from 'sherpa.save'.

        >>> restore()

        Load in the session from the given file:

        >>> restore('/data/m31/setup.sherpa')

        """
        _check_str_type(filename, "filename")

        fin = open(filename, 'rb')
        try:
            obj = pickle.load(fin)
        finally:
            fin.close()

        if not isinstance(obj, Session):
            raise ArgumentErr('nosession', filename)

        # Update optmethods, stats, and estmethods
        # obj.__dict__ should not clobber new classes!
        dicts = [self._methods, self._stats, self._estmethods]
        names = ['_methods', '_stats', '_estmethods']
        for name, dic in zip(names, dicts):
            # update current session with user definitions
            dic.update(obj.__dict__[name])
            # remove old items from pickle
            obj.__dict__.pop(name)

        # update current session with pickle
        self.__dict__.update(obj.__dict__)

        # TODO: is it possible for _model_autoassign_func to be None,
        # since NoNewAttributesAfterInit stops it from being set?
        #
        if self._model_autoassign_func is not None:
            for name, cmpt in self._model_components.items():
                self._model_autoassign_func(name, cmpt)

    def _get_show_data(self, id: Optional[IdType] = None) -> str:
        data_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            data_str += 'Data Set: %s\n' % id
            data_str += str(self.get_data(id)) + '\n\n'
        return data_str

    def _get_show_filter(self, id: Optional[IdType] = None) -> str:
        filt_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            filt_str += 'Data Set Filter: %s\n' % id
            filt_str += self.get_data(id).get_filter_expr() + '\n\n'
        return filt_str

    def _get_show_model(self, id: Optional[IdType] = None) -> str:
        model_str = ''
        ids = self.list_data_ids()
        mdl_ids = self.list_model_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            if id in mdl_ids:
                model_str += 'Model: %s\n' % id
                model_str += str(self.get_model(id)) + '\n\n'
        return model_str

    def _get_show_source(self, id: Optional[IdType] = None) -> str:
        model_str = ''
        ids = self.list_data_ids()
        src_ids = self._sources.keys()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            if id in src_ids:
                model_str += 'Model: %s\n' % id
                model_str += str(self.get_source(id)) + '\n\n'
        return model_str

    def _get_show_kernel(self, id: Optional[IdType] = None) -> str:
        kernel_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            if id in self._psf.keys():
                kernel_str += 'PSF Kernel: %s\n' % id
                # Show the PSF parameters
                kernel_str += str(self.get_psf(id)) + '\n\n'
        return kernel_str

    def _get_show_psf(self, id: Optional[IdType] = None) -> str:
        psf_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            if id in self._psf.keys():
                psf_str += 'PSF Model: %s\n' % id
                # Show the PSF dataset or PSF model
                psf_str += str(self.get_psf(id).kernel) + '\n\n'
        return psf_str

    def _get_show_method(self) -> str:
        return ('Optimization Method: %s\n%s\n' %
                (type(self._current_method).__name__,
                 str(self._current_method)))

    def _get_show_stat(self) -> str:
        return ('Statistic: %s\n%s\n' %
                (type(self._current_stat).__name__,
                 str(self._current_stat)))

    def _get_show_fit(self) -> str:
        if self._fit_results is None:
            return ''

        fit_str = self._get_show_method()
        fit_str += '\n'
        fit_str += self._get_show_stat()
        fit_str += '\n'
        fit_str += 'Fit:'
        fit_str += self.get_fit_results().format() + '\n\n'
        return fit_str

    def _get_show_conf(self) -> str:
        if self._confidence_results is None:
            return ''

        conf_str = 'Confidence:'
        conf_str += self.get_conf_results().format() + '\n\n'
        return conf_str

    def _get_show_proj(self) -> str:
        if self._projection_results is None:
            return ''

        proj_str = 'Projection:'
        proj_str += self.get_proj_results().format() + '\n\n'
        return proj_str

    def _get_show_covar(self) -> str:
        if self._covariance_results is None:
            return ''

        covar_str = 'Covariance:'
        covar_str += self.get_covar_results().format() + '\n\n'
        return covar_str

    def show_stat(self, outfile=None, clobber=False) -> None:
        """Display the current fit statistic.

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        calc_stat : Calculate the fit statistic for a data set.
        calc_stat_info : Display the statistic values for the current models.
        get_stat : Return a fit-statistic method.
        show_all : Report the current state of the Sherpa session.

        Examples
        --------

        >>> set_stat('cash')
        >>> show_stat()
        Statistic: Cash
        Maximum likelihood function

        """
        txt = self._get_show_stat()
        send_to_pager(txt, outfile, clobber)

    def show_method(self, outfile=None, clobber=False) -> None:
        """Display the current optimization method and options.

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        get_method : Return an optimization method.
        get_method_opt : Return one or all options of the current optimization method.
        show_all : Report the current state of the Sherpa session.

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
        txt = self._get_show_method()
        send_to_pager(txt, outfile, clobber)

    def show_fit(self, outfile=None, clobber=False) -> None:
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
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        fit : Fit one or more data sets.
        get_fit_results : Return the results of the last fit.
        list_data_ids : List the identifiers for the loaded data sets.
        list_model_ids : List of all the data sets with a source expression.
        show_all : Report the current state of the Sherpa session.

        """
        txt = self._get_show_fit()
        send_to_pager(txt, outfile, clobber)

    def show_data(self,
                  id: Optional[IdType] = None,
                  outfile=None, clobber=False) -> None:
        """Summarize the available data sets.

        Display information on the data sets that have been
        loaded. The details depend on the type of the data set
        (e.g. 1D, image, PHA files).

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then all data sets are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        list_data_ids : List the identifiers for the loaded data sets.
        show_all : Report the current state of the Sherpa session.

        """
        txt = self._get_show_data(id)
        send_to_pager(txt, outfile, clobber)

    def show_filter(self,
                    id: Optional[IdType] = None,
                    outfile=None, clobber=False) -> None:
        """Show any filters applied to a data set.

        Display any filters that have been applied to the independent
        axis or axes of the data set.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then all data sets are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        ignore : Exclude data from the fit.
        sherpa.astro.ui.ignore2d : Exclude a spatial region from an image.
        list_data_ids : List the identifiers for the loaded data sets.
        notice : Include data in the fit.
        sherpa.astro.ui.notice2d : Include a spatial region of an image.
        show_all : Report the current state of the Sherpa session.

        """
        txt = self._get_show_filter(id)
        send_to_pager(txt, outfile, clobber)

    def show_model(self,
                   id: Optional[IdType] = None,
                   outfile=None, clobber=False) -> None:
        """Display the model expression used to fit a data set.

        This displays the model used to fit the data set, that is,
        the expression set by `set_model` or `set_source`
        combined with any instrumental responses, together with the
        parameter values of the model. The `show_source` function
        displays just the source expression, without the instrumental
        components (if any).

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then all source expressions are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        list_model_ids : List of all the data sets with a source expression.
        set_model : Set the source model expression for a data set.
        show_all : Report the current state of the Sherpa session.
        show_source : Display the source model expression for a data set.

        """
        txt = self._get_show_psf(id)
        txt += self._get_show_model(id)
        send_to_pager(txt, outfile, clobber)

    def show_source(self,
                    id: Optional[IdType] = None,
                    outfile=None, clobber=False) -> None:
        """Display the source model expression for a data set.

        This displays the source model for a data set, that is, the
        expression set by `set_model` or `set_source`, as well as the
        parameter values for the model. The `show_model` function
        displays the model that is fit to the data; that is, it
        includes any instrument responses.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then all source expressions are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        list_model_ids : List of all the data sets with a source expression.
        set_model : Set the source model expression for a data set.
        show_all : Report the current state of the Sherpa session.
        show_model : Display the model expression used to fit a data set.

        """
        txt = self._get_show_source(id)
        send_to_pager(txt, outfile, clobber)

    # DOC-TODO: how and where to describe the PSF/kernel difference
    # as the Notes section below is inadequate
    def show_kernel(self,
                    id: Optional[IdType] = None,
                    outfile=None, clobber=False) -> None:
        """Display any kernel applied to a data set.

        The kernel represents the subset of the PSF model that is used
        to fit the data. The `show_psf` function shows the un-filtered
        version.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then all data sets are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

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
        txt = self._get_show_kernel(id)
        send_to_pager(txt, outfile, clobber)

    # DOC-TODO: how and where to describe the PSF/kernel difference
    # as the Notes section below is inadequate
    def show_psf(self,
                 id: Optional[IdType] = None,
                 outfile=None, clobber=False) -> None:
        """Display any PSF model applied to a data set.

        The PSF model represents the full model or data set that is
        applied to the source expression. The `show_kernel` function
        shows the filtered version.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then all data sets are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

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
        txt = self._get_show_psf(id)
        send_to_pager(txt, outfile, clobber)

    def show_conf(self, outfile=None, clobber=False) -> None:
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
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
        show_all : Report the current state of the Sherpa session.

        """
        txt = self._get_show_conf()
        send_to_pager(txt, outfile, clobber)

    def show_proj(self, outfile=None, clobber=False) -> None:
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
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        proj : Estimate parameter confidence intervals using the projection method.
        show_all : Report the current state of the Sherpa session.

        """
        txt = self._get_show_proj()
        send_to_pager(txt, outfile, clobber)

    def show_covar(self, outfile=None, clobber=False) -> None:
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
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        covar : Estimate parameter confidence intervals using the covariance method.
        show_all : Report the current state of the Sherpa session.

        """
        txt = self._get_show_covar()
        send_to_pager(txt, outfile, clobber)

    def show_all(self,
                 id: Optional[IdType] = None,
                 outfile=None, clobber=False) -> None:
        """Report the current state of the Sherpa session.

        Display information about one or all of the data sets that
        have been loaded into the Sherpa session. The information
        shown includes that provided by the other ``show_xxx`` routines,
        and depends on the type of data that is loaded.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then all data sets are
           displayed.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        clean : Clear all stored session data.
        list_data_ids : List the identifiers for the loaded data sets.
        save : Save the current Sherpa session to a file.
        sherpa.astro.ui.save_all : Save the Sherpa session as an ASCII file.
        sherpa.astro.ui.show_bkg : Show the details of the PHA background data sets.
        sherpa.astro.ui.show_bkg_model : Display the background model expression for a data set.
        sherpa.astro.ui.show_bkg_source : Display the background model expression for a data set.
        show_conf : Display the results of the last conf evaluation.
        show_covar : Display the results of the last covar evaluation.
        show_data : Summarize the available data sets.
        show_filter : Show any filters applied to a data set.
        show_fit : Summarize the fit results.
        show_kernel : Display any kernel applied to a data set.
        show_method : Display the current optimization method and options.
        show_model : Display the model expression used to fit a data set.
        show_proj : Display the results of the last proj evaluation.
        show_psf : Display any PSF model applied to a data set.
        show_source : Display the source model expression for a data set.
        show_stat : Display the current fit statistic.

        """
        txt = self._get_show_data(id)
        txt += self._get_show_model(id)
        txt += self._get_show_fit()
        txt += self._get_show_conf()
        txt += self._get_show_proj()
        txt += self._get_show_covar()
        send_to_pager(txt, outfile, clobber)

    def get_functions(self) -> list[str]:
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
            if not func.startswith('_') and callable(getattr(self, func)):
                funcs.append(func)
        return funcs

    def list_functions(self, outfile=None, clobber=False) -> None:
        """Display the functions provided by Sherpa.

        Unlike the other ``list_xxx`` commands, this does not
        return an array. Instead it acts like the ``show_xxx``
        family of commands.

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is ``False``.

        See Also
        --------
        get_functions : Return the functions provided by Sherpa.
        show_all : Report the current state of the Sherpa session.

        """
        funcs_list = self.get_functions()
        funcs = ''
        for func in funcs_list:
            funcs += '%s\n' % func

        send_to_pager(funcs, outfile, clobber)

    ###########################################################################
    # IDs and general data management
    ###########################################################################

    @staticmethod
    def _valid_id(id: Any) -> bool:  # TODO: mark as TypeGuard[IdType] instead
        """Is the identifier valid for Sherpa?

        This does not treat None as a valid identifier.
        """
        return (_is_integer(id) or _is_str(id))

    def _fix_id(self, id: Optional[IdType]) -> IdType:
        """Validate the dataset id.

        The identifier can be any string or integer except for the
        plot and contour types that are supported by the `plot` and
        `contour` methods.

        Parameters
        ----------
        id : int or str or None
            The dataset identifier. If set to None then the default
            identifier, returned by `get_default_id`, is used.

        Returns
        -------
        id : int or str
            The identifier to use (it will only differ from the input
            parameter was set to None).

        Raises
        ------
        sherpa.utils.err.ArgumentTypeErr
            If the identifier was not a string or an integer.
        sherpa.utils.err.IdentifierErr
            If the identifier was invalid.

        See Also
        --------
        get_default_id, set_default_id

        Notes
        -----
        Since there is currently no way to set the default background
        id of the DataPHA class (e.g. in unpack_pha) we do not use the
        _default_id setting here.

        """

        if id is None:
            return self._default_id

        if not self._valid_id(id):
            raise ArgumentTypeErr("intstr")

        if _is_integer(id):
            return id

        if self._check_plottype(id) or self._check_contourtype(id):
            raise IdentifierErr("badid", id)

        return id

    def _get_plottype(self, plottype: str) -> str:
        """Return the name to refer to a given plot type.

        This supports aliases for the plot name. If an alias is used
        then a warning message is logged, telling the user the name
        they should use instead.

        Parameters
        ----------
        plottype : str
            The requested plot type, such as "data".

        Returns
        -------
        answer : str
            The actual plot type (which may match the input).

        Raises
        ------
        PlotErr
            The plottype argument is not valid.

        """

        if plottype in self._plot_types:
            return plottype

        try:
            answer = self._plot_types_alias[plottype]
        except KeyError:
            allowed = list(self._plot_types)
            raise PlotErr("wrongtype", plottype, str(allowed)) from None

        # This could be a warnings.warn(..., DeprecationWarning)
        # message but then most users would not see the message, so we
        # use the logger interface instead. It does mean that the
        # message will be repeated each time it is used.
        #
        warning("The argument '%s' is deprecated and "
                "'%s' should be used instead", plottype, answer)
        return answer

    def _check_plottype(self, plottype: IdType) -> bool:
        """Is this a valid plot type (including aliases)?"""

        return plottype in self._plot_types or \
            plottype in self._plot_types_alias

    def _get_contourtype(self, plottype: str):
        """Return the name to refer to a given contour type."""

        if plottype in self._contour_types:
            return plottype

        allowed = list(self._contour_types)
        raise PlotErr("wrongtype", plottype, str(allowed))

    def _check_contourtype(self, plottype: IdType) -> bool:
        """Is this a valid contour type?"""

        return plottype in self._contour_types

    # The assumption is that itemdict has a single type for its
    # values.
    #
    def _get_item(self,
                  id: Optional[IdType],
                  itemdict: dict[IdType, T],
                  itemdesc: str,
                  errdesc: str
                  ) -> T:
        id = self._fix_id(id)
        item = itemdict.get(id)
        if item is None:
            raise IdentifierErr('getitem', itemdesc, id, errdesc)
        return item

    def _set_item(self,
                  id: Optional[IdType],
                  item: T,
                  itemdict: dict[IdType, T],
                  itemtype,
                  itemname: str,
                  itemdesc: str
                  ) -> None:
        id = self._fix_id(id)
        _check_type(item, itemtype, itemname, itemdesc)
        itemdict[id] = item

    def get_default_id(self) -> IdType:
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
           ``None``.

        See Also
        --------
        list_data_ids : List the identifiers for the loaded data sets.
        set_default_id : Set the default data set identifier.

        Notes
        -----
        The default Sherpa data set identifier is the integer 1.

        Examples
        --------

        Display the default identifier:

        >>> print(get_default_id())

        Store the default identifier and use it as an argument to
        call another Sherpa routine:

        >>> defid = get_default_id()
        >>> load_arrays(defid, x, y)

        """
        return self._default_id

    def set_default_id(self, id: IdType) -> None:
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
           ``None``.

        See Also
        --------
        get_default_id : Return the default data set identifier.
        list_data_ids : List the identifiers for the loaded data sets.

        Notes
        -----
        The default Sherpa data set identifier is the integer 1.

        Examples
        --------

        After the following, many commands, such as `set_source`, will
        use 'src' as the default data set identifier:

        >>> set_default_id('src')

        Restore the default data set identifier.

        >>> set_default_id(1)

        """
        self._default_id = self._fix_id(id)

    ###########################################################################
    # Optimization methods
    ###########################################################################

    def list_methods(self) -> list[str]:
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
        keys = list(self._methods.keys())
        keys.sort()
        return keys

    def _get_method_by_name(self, name: str) -> OptMethod:
        meth = self._methods.get(name.lower())
        if meth is None:
            raise ArgumentErr('badmethod', name)
        return meth

    def get_method(self,
                   name: Optional[str] = None
                   ) -> OptMethod:
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
           If the ``name`` argument is not recognized.

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

        _check_str_type(name, "name")
        return self._get_method_by_name(name)

    # DOC-TODO: is this guaranteed to be the same as get_method().name
    # or get_method().name.lower() and, if so, shouldn't this be
    # how it is coded?
    def get_method_name(self) -> str:
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

    # DOC-TODO: remove the list of supported methods once the
    # relevant documentation has been updated.
    #
    def set_method(self,
                   meth: Union[OptMethod, str]
                   ) -> None:
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
           If the ``meth`` argument is not recognized.

        See Also
        --------
        get_method_name : Return the name of the current optimization method.
        list_methods : List the supported optimization methods.
        set_stat : Set the fit statistic.

        Notes
        -----
        The available methods include:

        ``levmar``
           The Levenberg-Marquardt method is an interface to the
           MINPACK subroutine lmdif to find the local minimum of
           nonlinear least squares functions of several variables by a
           modification of the Levenberg-Marquardt algorithm [1].

        ``moncar``
           The implementation of the moncar method is based on [2].

        ``neldermead``
           The implementation of the Nelder Mead Simplex direct search
           is based on [3].

        ``simplex``
           This is another name for ``neldermead``.

        References
        ----------

        1. J.J. More, "The Levenberg Marquardt algorithm:
           implementation and theory," in Lecture Notes in Mathematics
           630: Numerical Analysis, G.A. Watson (Ed.), Springer-Verlag:
           Berlin, 1978, pp.105-116.

        2. Storn, R. and Price, K. "Differential Evolution: A
           Simple and Efficient Adaptive Scheme for Global Optimization
           over Continuous Spaces." J. Global Optimization 11, 341-359,
           1997.

        3. Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
           Paul E. Wright "Convergence Properties of the Nelder-Mead
           Simplex Algorithm in Low Dimensions", SIAM Journal on
           Optimization,Vol. 9, No. 1 (1998), pages 112-147.

        Examples
        --------

        >>> set_method('neldermead')

        """
        if _is_str(meth):
            method = self._get_method_by_name(meth)
        else:
            _check_type(meth, OptMethod, 'meth',
                        'a method name or object')
            method = meth

        self._current_method = method

        # Do we need to set the RNG argument? This is not ideal, but
        # is a band-aid while we work out how to handle the RNG
        # handling in the UI layer. One option would only to update
        # this if the rng setting is None, but that would be
        # problematic to users who change the method (as it could
        # leave an "old" RNG state lying around), although it does
        # mean that if the user has explicitly set rng then it will be
        # over-written here. An alternative would be to remove the
        # "rng" setting here (i.e. hide it somehow), but then how does
        # this get passed to the method?
        #
        if "rng" in self._current_method.config:
            self._current_method.config["rng"] = self.get_rng()

    def _check_method_opt(self, optname: str) -> None:
        _check_str_type(optname, "optname")
        if optname not in self._current_method.config:
            raise ArgumentErr('badopt', optname, self.get_method_name())

    def get_method_opt(self,
                       optname: Optional[str] = None
                       ) -> Any:
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
           If the ``optname`` argument is not recognized.

        See Also
        --------
        get_method : Return an optimization method.
        set_method : Change the optimization method.
        set_method_opt : Change an option of the current optimization method.

        Examples
        --------

        >>> get_method_opt('maxfev') is None
        True

        >>> mopts = get_method_opt()
        >>> mopts['maxfev'] is None
        True

        """
        if optname is None:
            return self._current_method.config
        self._check_method_opt(optname)
        return self._current_method.config[optname]

    def set_method_opt(self, optname: str, val: Any) -> None:
        """Set an option for the current optimization method.

        This is a helper function since the optimization options can also
        be set directly using the object returned by `get_method`.

        Parameters
        ----------
        optname : str
           The name of the option to set. The `get_method`
           and `get_method_opt` routines can be used to find
           out valid values for this argument.
        val
           The new value for the option.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the ``optname`` argument is not recognized.

        See Also
        --------
        get_method : Return an optimization method.
        get_method_opt : Return one or all options of the current optimization method.
        set_method : Change the optimization method.

        Examples
        --------

        Change the ``maxfev`` parameter for the current optimizer
        to 2000.

        >>> set_method_opt('maxfev', 2000)

        """
        self._check_method_opt(optname)
        self._current_method.config[optname] = val

    def get_iter_method_name(self) -> str:
        """Return the name of the iterative fitting scheme.

        Returns
        -------
        name : {'none', 'sigmarej'}
           The name of the iterative fitting scheme set by
           `set_iter_method`.

        See Also
        --------
        list_iter_methods : List the iterative fitting schemes.
        set_iter_method : Set the iterative-fitting scheme used in the fit.

        Examples
        --------

        >>> print(get_iter_method_name())

        """
        return self._current_itermethod['name']

    def get_iter_method_opt(self,
                            optname: Optional[str] = None
                            ) -> Any:
        """Return one or all options for the iterative-fitting scheme.

        The options available for the iterative-fitting methods are
        described in `set_iter_method_opt`.

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
           If the ``optname`` argument is not recognized.

        See Also
        --------
        get_iter_method_name : Return the name of the iterative fitting scheme.
        set_iter_method_opt : Set an option for the iterative-fitting scheme.
        set_iter_method : Set the iterative-fitting scheme used in the fit.

        Examples
        --------

        Display the settings of the current iterative-fitting method:

        >>> print(get_iter_method_opt())

        Switch to the sigmarej scheme and find out the current settings:

        >>> set_iter_method('sigmarej')
        >>> opts = get_iter_method_opt()

        Return the 'maxiters' setting (if applicable):

        >>> get_iter_method_opt('maxiters')

        """
        itermethod_opts = dict(self._current_itermethod)
        del itermethod_opts['name']

        if optname is None:
            return itermethod_opts

        _check_str_type(optname, "optname")
        if optname not in itermethod_opts:
            raise ArgumentErr(
                'badopt', optname, self._current_itermethod['name'])
        return itermethod_opts[optname]

    def list_iter_methods(self) -> list[str]:
        """List the iterative fitting schemes.

        Returns
        -------
        schemes : list of str
           A list of the names that can be used with `set_iter_method`.

        See Also
        --------
        get_iter_method_name : Return the name of the iterative fitting scheme.
        set_iter_method : Set the iterative-fitting scheme used in the fit.

        Examples
        --------

        >>> list_iter_methods()
        ['none', 'sigmarej']

        """
        keys = list(self._itermethods.keys())
        keys.sort()
        return keys

    # DOC-TODO: this information is also in sherpa/fit.py
    # DOC-TODO: this raises a ValueError rather than a Sherpa error class
    def set_iter_method(self, meth: str) -> None:
        """Set the iterative-fitting scheme used in the fit.

        Control whether an iterative scheme should be applied to
        the fit.

        .. versionchanged:: 4.14.1
           The "primini" scheme has been removed from Sherpa.

        Parameters
        ----------
        meth : { 'none', 'sigmarej' }
           The name of the scheme used during the fit; 'none' means no
           scheme is used. It is only valid to change the scheme
           when a chi-square statistic is in use.

        Raises
        ------
        TypeError
           When the ``meth`` argument is not recognized.

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
        The parameters of the schemes are described in
        `set_iter_method_opt`.

        This is a chi-square statistic where the variance is computed
        from model amplitudes derived in the previous iteration of the
        fit. This 'Iterative Weighting' ([1]) attempts to remove
        biased estimates of model parameters.

        The variance in bin i is estimated to be::

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

        The ``sigmarej`` scheme is based on the
        `IRAF ``sfit`` function <https://iraf.readthedocs.io/en/latest/tasks/noao/imred/specred/sfit.html>`_,
        where after a fit data points are excluded if the value
        of ``(data-model)/error)`` exceeds a threshold, and the data
        re-fit. This removal of data points continues until the fit
        has converged. The error removal can be asymmetric, since
        there are separate parameters for the lower and upper limits.

        References
        ----------

        1. `"Multiparameter linear least-squares fitting to Poisson
           data one count at a time", Wheaton et al. 1995, ApJ 438, 322
           <https://adsabs.harvard.edu/abs/1995ApJ...438..322W>`_

        Examples
        --------

        Switch to the 'sigmarej' scheme for iterative fitting and
        change the low and high rejection limits to 4 and 3
        respectively:

        >>> set_iter_method('sigmarej')
        >>> set_iter_method_opt('lrej') = 4
        >>> set_iter_method_opt('hrej') = 3

        Remove any iterative-fitting method:

        >>> set_iter_method('none')

        """
        _check_str_type(meth, "meth")

        if meth not in self._itermethods:
            raise TypeError(f'{meth} is not an iterative fitting method')

        self._current_itermethod = self._itermethods[meth]

    def set_iter_method_opt(self, optname: str, val: Any) -> None:
        """Set an option for the iterative-fitting scheme.

        Parameters
        ----------
        optname : str
           The name of the option to set. The
           `get_iter_method_opt` routine can be used to find
           out valid values for this argument.
        val
           The new value for the option.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the ``optname`` argument is not recognized.

        See Also
        --------
        get_iter_method_name : Return the name of the iterative fitting scheme.
        get_iter_method_opt : Return one or all options for the iterative-fitting scheme.
        list_iter_methods : List the iterative fitting schemes.
        set_iter_method : Set the iterative-fitting scheme used in the fit.

        Notes
        -----
        The supported fields for the ``sigmarej`` scheme are:

        grow
           The number of points adjacent to a rejected point that
           should also be removed. A value of ``0`` means that only the
           discrepant point is removed whereas a value of ``1`` means
           that the two adjacent points (one lower and one higher)
           will also be removed.

        hrej
           The rejection criterion in units of sigma, for data
           points above the model (it must be >= 0).

        lrej
           The rejection criterion in units of sigma, for data
           points below the model (it must be >= 0).

        maxiters
           The maximum number of iterations to perform. If this
           value is ``0`` then the fit will run until it has
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
        _check_str_type(optname, "optname")
        if (optname not in self._current_itermethod or
                optname == 'name'):
            raise ArgumentErr(
                'badopt', optname, self._current_itermethod['name'])

        self._current_itermethod[optname] = val

    ###########################################################################
    # Statistics
    ###########################################################################

    def list_stats(self) -> list[str]:
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
         'leastsq',
         'wstat']

        """
        keys = list(self._stats.keys())
        keys.sort()
        return keys

    def _get_stat_by_name(self, name: str) -> Stat:
        stat = self._stats.get(name.lower())
        if stat is None:
            raise ArgumentErr('badstat', name)
        return stat

    def get_stat(self,
                 name: Optional[str] = None
                 ) -> Stat:
        """Return the fit statisic.

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
           If the ``name`` argument is not recognized.

        See Also
        --------
        get_stat_name : Return the name of the current fit statistic.
        list_stats : List the fit statistics.
        set_stat : Change the fit statistic.

        Examples
        --------

        Return the currently-selected statistic, display its name, and
        read the help documentation for it:

        >>> stat = get_stat()
        >>> stat.name
        'chi2gehrels'
        >>> help(stat)

        Read the help for the "wstat" statistic:

        >>> help(get_stat('wstat'))

        """
        if name is None:
            return self._current_stat

        _check_str_type(name, "name")
        return self._get_stat_by_name(name)

    def get_stat_name(self) -> str:
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

        >>> get_stat_name()
        'chi2gehrels'

        >>> set_stat('cash')
        >>> get_stat_name()
        'cash'

        """
        return type(self.get_stat()).__name__.lower()

    def set_stat(self,
                 stat: Union[str, Stat]
                 ) -> None:
        """Set the statistical method.

        Changes the method used to evaluate the fit statistic, that is
        the numerical measure that determines how closely the model
        represents the data.

        Parameters
        ----------
        stat : str or sherpa.stats.Stat instance
           When a string, the name of the statistic (case is not
           important): see `list_stats()` for supported
           values. Otherwise an instance of the statistic to use.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the ``stat`` argument is not recognized.

        See Also
        --------
        calc_stat : Calculate the statistic value for a dataset.
        get_stat_name : Return the current statistic method.
        list_stats : List the supported fit statistics.
        load_user_stat : Create a user-defined statistic.

        Notes
        -----
        The available statistics include:

        cash
           A maximum likelihood function [1].

        chi2
           Chi-squared statistic using the supplied error values.

        chi2constvar
           Chi-squared with constant variance computed from the counts
           data.

        chi2datavar
           Chi-squared with data variance. If the data has 0 counts then
           the error for that bin is 0.

        chi2gehrels
           Chi-squared with gehrels method [2]. This is the default method.

        chi2modvar
           Chi-squared with model amplitude variance.

        chi2xspecvar
           Chi-squared with data variance to match XSPEC. Errors from
           zero-count channels (source or background) are ignored if
           the other channel (background or source) contains counts,
           or replaced by a minimum value (when both source or
           background are empty). It should not be used when a model
           is fit to the background rather than the background is
           subtracted from the data.

        cstat
           A maximum likelihood function (the XSPEC implementation of
           the Cash function) [3]. This does *not* include support
           for including the background.

        wstat
           A maximum likelihood function which includes the background
           data as part of the fit (i.e. for when it is not being
           explicitly modelled) (the XSPEC implementation of the Cash
           function) [3].

        leastsq
           The least-squares statisic (the error is not used in this
           statistic).

        References
        ----------

        1. `Cash, W. "Parameter estimation in astronomy through
           application of the likelihood ratio", ApJ, vol 228,
           p. 939-947 (1979).
           <https://adsabs.harvard.edu/abs/1979ApJ...228..939C>`_

        2. `Gehrels, N. "Confidence limits for small numbers of
           events in astrophysical data", ApJ, vol 303, p. 336-346 (1986).
           <https://adsabs.harvard.edu/abs/1986ApJ...303..336G>`_

        3. https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html

        Examples
        --------

        >>> set_stat('cash')

        """
        if _is_str(stat):
            statobj = self._get_stat_by_name(stat)
        else:
            _check_type(stat, Stat, 'stat',
                        'a statistic name or object')
            statobj = stat

        self._current_stat = statobj

    ###########################################################################
    # Data sets
    ###########################################################################

    def list_data_ids(self) -> list[IdType]:
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
        keys = list(self._data.keys())
        keys.sort(key=str)  # always sort by string value.
        return keys

    def get_data(self, id: Optional[IdType] = None) -> Data:
        """Return the data set by identifier.

        The object returned by the call can be used to query and
        change properties of the data set.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.

        Returns
        -------
        instance : sherpa.data.Data
           The data instance.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           No data has been loaded for this data set.

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

    def _get_data(self, id: Optional[IdType]) -> Optional[Data]:
        """Return a data set or None.

        The same as get_data except that it returns None if the
        dataset does not exist.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.

        Returns
        -------
        instance : sherpa.data.Data or None

        """

        return self._data.get(self._fix_id(id))

    # DOC-TODO: terrible synopsis
    def set_data(self, id, data=None) -> None:

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
    # DOC-NOTE: is ncols really 2 here? Does it make sense?
    def load_staterror(self, id, filename=None, ncols=2,
                       *args, **kwargs) -> None:
        """Load the statistical errors from an ASCII file.

        Read in a column or image from a file and use the values
        as the statistical errors for a data set. This over rides
        the errors calculated by any statistic, such as
        ``chi2gehrels`` or ``chi2datavar``.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the ASCII file to read in.
        ncols : int, optional
           The number of columns to read in (the first ``ncols`` columns
           in the file).
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           ``None``.
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.

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

        See `unpack_data` for a description of the supported
        file format.

        Examples
        --------

        Read in the first column from 'tbl.dat':

        >>> load_staterror('tbl.dat')

        Use the column labelled 'col3'

        >>> load_staterror('tbl.dat', colkeys=['col3'])

        Read in the first column from the file 'errors.dat' as the
        statistical errors for the 'core' data set:

        >>> load_staterror('core', 'errors.dat')

        """
        if filename is None:
            id, filename = filename, id

        self.set_staterror(id, self._read_error(filename, ncols=ncols,
                                                *args, **kwargs))

    # DOC-NOTE: also in sherpa.astro.utils
    # DOC-NOTE: is ncols really 2 here? Does it make sense?
    def load_syserror(self, id, filename=None, ncols=2,
                      *args, **kwargs) -> None:
        """Load the systematic errors from an ASCII file.

        Read in a column or image from a file and use the values
        as the systematic errors for a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the ASCII file to read in.
        ncols : int, optional
           The number of columns to read in (the first ``ncols`` columns
           in the file).
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           ``None``.
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.

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

        See `unpack_data` for a description of the supported
        file format.

        Examples
        --------

        Read in the first column from 'tbl.dat':

        >>> load_syserror('tbl.dat')

        Use the column labelled 'col3'

        >>> load_syserror('tbl.dat', colkeys=['col3'])

        Read in the first column from the file 'errors.dat' as the
        systematic errors for the 'core' data set:

        >>> load_syserror('core', 'errors.dat')

        """
        if filename is None:
            id, filename = filename, id

        self.set_syserror(id, self._read_error(filename, ncols=ncols,
                                               *args, **kwargs))

    # DOC-NOTE: also in sherpa.astro.utils
    # DOC-TODO: does ncols make sense here? (have removed for now)
    def load_filter(self, id, filename=None, ignore=False, ncols=2,
                    *args, **kwargs) -> None:
        """Load the filter array from an ASCII file and add to a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the ASCII file that contains the filter
           information.
        ignore : bool, optional
           If ``False`` (the default) then include bins with a non-zero
           filter value, otherwise exclude these bins.
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           ``None``.
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.

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

        See `unpack_data` for a description of the supported
        file format.

        Examples
        --------

        Read in the first column of the file and apply it to the
        default data set:

        >>> load_filter('filt.dat')

        Select the FILTER column of the file:

        >>> load_filter(2, 'filt.dat', colkeys=['FILTER'])

        """
        if filename is None:
            id, filename = filename, id

        self.set_filter(id,
                        self._read_error(
                            filename, ncols=ncols, *args, **kwargs),
                        ignore=ignore)

    def set_filter(self, id, val=None, ignore=False) -> None:
        """Set the filter array of a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        val : array
           The array of filter values (``0`` or ``1``). The size should
           match the array returned by `get_dep`.
        ignore : bool, optional
           If ``False`` (the default) then include bins with a non-zero
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
        >>> set_filter(d >= 20)

        """
        if val is None:
            val, id = id, val

        d = self.get_data(id)
        set_filter(d, val, ignore=ignore)

    # also in sherpa.astro.utils
    def set_dep(self, id, val=None) -> None:
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

        set_dep(d, val)

    # DOC-NOTE: also in sherpa.utils
    def set_staterror(self, id, val=None, fractional=False) -> None:
        """Set the statistical errors on the dependent axis of a data set.

        These values override the errors calculated by any statistic,
        such as ``chi2gehrels`` or ``chi2datavar``.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        val : array or scalar
           The systematic error.
        fractional : bool, optional
           If ``False`` (the default value), then the `val` parameter is
           the absolute value, otherwise the `val` parameter
           represents the fractional error, so the absolute value is
           calculated as ``get_dep() * val`` (and `val` must be
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
        in ``dys`` (a scalar or an array):

        >>> set_staterror(dys)

        Set the statistical error on the 'core' data set to be 5% of
        the data values:

        >>> set_staterror('core', 0.05, fractional=True)

        """
        if val is None:
            val, id = id, val

        d = self.get_data(id)
        set_error(d, "staterror", val, fractional=fractional)

    # DOC-NOTE: also in sherpa.astro.utils
    def set_syserror(self, id, val=None, fractional=False) -> None:
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
           If ``False`` (the default value), then the `val` parameter is
           the absolute value, otherwise the `val` parameter
           represents the fractional error, so the absolute value is
           calculated as ``get_dep() * val`` (and `val` must be
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
        in ``dys`` (a scalar or an array):

        >>> set_syserror(dys)

        Set the systematic error on the 'core' data set to be 5% of
        the data values:

        >>> set_syserror('core', 0.05, fractional=True)

        """
        if val is None:
            val, id = id, val
        err = None

        d = self.get_data(id)
        set_error(d, "syserror", val, fractional=fractional)

    # DOC-NOTE: also in sherpa.astro.utils
    def get_staterror(self,
                      id: Optional[IdType] = None,
                      filter=False):
        """Return the statistical error on the dependent axis of a data set.

        The function returns the statistical errors on the values
        (dependenent axis) of a data set. These may have been set
        explicitly - either when the data set was created or with
        a call to `set_staterror` - or as defined by the chosen
        fit statistic (such as "chi2gehrels").

        Parameters
        ----------
        id : int, str, or None, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is ``False``.

        Returns
        -------
        staterrors : array
           The statistical error for each data point. This may be
           estimated from the data (e.g. with the ``chi2gehrels``
           statistic) or have been set explicitly (`set_staterror`).
           The size of this array depends on the `filter` argument.

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

        Notes
        -----
        The default behavior is to not apply any filter defined on the
        independent axes to the results, so that the return value is for
        all points (or bins) in the data set. Set the `filter` argument
        to `True` to apply this filter.

        Examples
        --------

        If not explicitly given, the statistical errors on a data set
        may be calculated from the data values (the independent axis),
        depending on the chosen statistic:

        >>> load_arrays(1, [10, 15, 19], [4, 5, 9])
        >>> set_stat('chi2datavar')
        >>> get_staterror()
        array([ 2.        ,  2.23606798,  3.        ])
        >>> set_stat('chi2gehrels')
        >>> get_staterror()
        array([ 3.17944947,  3.39791576,  4.122499  ])

        If the statistical errors are set - either when the data set
        is created or with a call to `set_staterror` - then these values
        will be used, no matter the statistic:

        >>> load_arrays(1, [10, 15, 19], [4, 5, 9], [2, 3, 5])
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
    def get_syserror(self,
                     id: Optional[IdType] = None,
                     filter=False):
        """Return the systematic error on the dependent axis of a data set.

        The function returns the systematic errors on the values
        (dependenent axis) of a data set. It is an error if called
        on a data set with no systematic errors (which are set
        with `set_syserror`).

        Parameters
        ----------
        id : int, str, or None, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is ``False``.

        Returns
        -------
        syserrors : array
           The systematic error for each data point. The size of this
           array depends on the `filter` argument.

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

        Notes
        -----
        The default behavior is to not apply any filter defined on the
        independent axes to the results, so that the return value is for
        all points (or bins) in the data set. Set the `filter` argument
        to `True` to apply this filter.

        Examples
        --------

        Return the systematic error for the default data set:

        >>> yerr = get_syserror()

        Return an array that has been filtered to match the data:

        >>> yerr = get_syserror(filter=True)

        Return the filtered errors for data set "core":

        >>> yerr = get_syserror("core", filter=True)

        """
        d = self.get_data(id)
        id = self._fix_id(id)
        err = d.get_syserror(filter)
        if err is None:
            raise sherpa.utils.err.DataErr('nosyserr', id)

        return err

    # DOC-NOTE: also in sherpa.astro.utils
    def get_error(self,
                  id: Optional[IdType] = None,
                  filter=False):
        """Return the errors on the dependent axis of a data set.

        The function returns the total errors (a quadrature addition
        of the statistical and systematic errors) on the values
        (dependent axis) of a data set. The individual components
        can be retrieved with the `get_staterror` and `get_syserror`
        functions.

        Parameters
        ----------
        id : int, str, or None, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is ``False``.

        Returns
        -------
        errors : array
           The error for each data point, formed by adding the
           statistical and systematic errors in quadrature.
           The size of this array depends on the `filter` argument.

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

        Notes
        -----
        The default behavior is to not apply any filter defined on the
        independent axes to the results, so that the return value is for
        all points (or bins) in the data set. Set the `filter` argument
        to `True` to apply this filter.

        Examples
        --------

        Return the error values for the default data set, ignoring any
        filter applied to it:

        >>> err = get_error()

        Ensure that the return values are for the selected (filtered)
        points in the default data set (the return array may be smaller
        than in the previous example):

        >>> err = get_error(filter=True)

        Find the errors for the "core" data set:

        >>> err = get_error('core', filter=True)

        """
        return self.get_data(id).get_error(filter,
                                           self.get_stat().calc_staterror)

    # DOC-NOTE: also in sherpa.astro.utils
    # DOC-NOTE: shouldn't this expose a filter parameter?
    def get_indep(self,
                  id: Optional[IdType] = None):
        """Return the independent axes of a data set.

        This function returns the coordinates of each point, or pixel,
        in the data set.

        Parameters
        ----------
        id : int, str, or None, optional
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

        >>> load_arrays(1, [10, 15, 19], [4, 5, 9])
        >>> get_indep()
        (array([10, 15, 19]),)

        For a 2D data set the X0 and X1 values are returned:

        >>> x0 = [10, 15, 12, 19]
        >>> x1 = [12, 14, 10, 17]
        >>> y = [4, 5, 9, -2]
        >>> load_arrays(2, x0, x1, y, Data2D)
        >>> get_indep(2)
        (array([10, 15, 12, 19]), array([12, 14, 10, 17]))

        """
        return self.get_data(id).get_indep()

    # DOC-NOTE: also in sherpa.astro.utils
    def get_dep(self,
                id: Optional[IdType] = None,
                filter=False):
        """Return the dependent axis of a data set.

        This function returns the data values (the dependent axis)
        for each point or pixel in the data set.

        Parameters
        ----------
        id : int, str, or None, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is ``False``.

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

        >>> load_arrays(1, [10, 15, 19], [4, 5, 9])
        >>> get_dep()
        array([4, 5, 9])

        >>> x0 = [10, 15, 12, 19]
        >>> x1 = [12, 14, 10, 17]
        >>> y = [4, 5, 9, -2]
        >>> load_arrays(2, x0, x1, y, Data2D)
        >>> get_dep(2)
        array([ 4,  5,  9, -2])

        If the ``filter`` flag is set then the return will be limited to
        the data that is used in the fit:

        >>> load_arrays(1, [10, 15, 19], [4, 5, 9])
        >>> ignore_id(1, 17, None)
        >>> get_dep()
        array([4, 5, 9])
        >>> get_dep(filter=True)
        array([4, 5])

        """
        return self.get_data(id).get_y(filter)

    def get_dims(self,
                 id: Optional[IdType] = None,
                 filter=False):
        """Return the dimensions of the data set.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        filter : bool, optional
           If ``True`` then apply any filter to the data set before
           returning the dimensions. The default is ``False``.

        Returns
        -------
        dims : a tuple of int

        See Also
        --------
        ignore : Exclude data from the fit.
        sherpa.astro.ui.ignore2d : Exclude a spatial region from an image.
        notice : Include data in the fit.
        sherpa.astro.ui.notice2d : Include a spatial region of an image.

        Examples
        --------

        Display the dimensions for the default data set:

        >>> print(get_dims())

        Find the number of bins in dataset 'a2543' without and with
        any filters applied to it:

        >>> nall = get_dims('a2543')
        >>> nfilt = get_dims('a2543', filter=True)

        """
        return self.get_data(id).get_dims(filter)

    # DOC-NOTE: should there be a version in sherpa.astro.utils with a bkg_id
    # parameter?
    def get_filter(self,
                   id: Optional[IdType] = None,
                   format: Optional[str] = None,
                   delim: Optional[str] = None
                   ) -> str:
        """Return the filter expression for a data set.

        This returns the filter expression, created by one or more
        calls to `ignore` and `notice`, for the data set.

        .. versionchanged:: 4.17.0
           The format and delim arguments can now be set.

        .. versionchanged:: 4.14.0
           The filter expressions have been tweaked for Data1DInt and
           PHA data sets (when using energy or wavelength units) and
           now describe the full range of the bins, rather than the
           mid-points.

        Parameters
        ----------
        id : int, str, or None, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        format : str or None, optional
           If set, use this rather than the default format value for
           the dataset.
        delim : str or None, optional
           If set, use this rather than the default delim value for
           the dataset.

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

        The default filter is the full dataset, given in the format
        ``lowval:hival`` (for a `Data1D` dataset like this these are
        inclusive limits):

        >>> load_arrays(1, [10, 15, 20, 25], [5, 7, 4, 2])
        >>> get_filter()
        '10.0000:25.0000'

        Change the formatting of the output:

        >>> get_filter(format="%d", delim="-")
        "10-25"

        The `notice` call restricts the data to the range between
        14 and 30. The resulting filter is the combination of this
        range and the data:

        >>> notice(14, 30)
        >>> get_filter()
        '15.0000:25.0000'

        Ignoring the point at ``x=20`` means that only the points at
        ``x=15`` and ``x=25`` remain, so a comma-separated list is used:

        >>> ignore(19, 22)
        >>> get_filter()
        '15.0000,25.0000'

        The filter equivalent to the per-bin array of filter values:

        >>> set_filter([1, 1, 0, 1])
        >>> get_filter()
        '10.0000:15.0000,25.0000'

        For an integrated data set (Data1DInt and DataPHA with energy
        or wavelength units)

        >>> load_arrays(1, [10, 15, 20, 25], [15, 20, 23, 30], [5, 7, 4, 2], Data1DInt)
        >>> get_filter()
        '10.0000:30.0000'

        For integrated datasets the limits are now inclusive only for
        the lower limit, but in this the end-point ends within a bin
        so is is included:

        >>> notice(17, 28)
        >>> get_filter()
        '15.0000:30.0000'

        There is no data in the range 23 to 24 so the ignore doesn't
        change anything:

        >>> ignore(23, 24)
        >>> get_filter()
        '15.0000:30.0000'

        However it does match the range 22 to 23 and so changes
        the filter:

        >>> ignore(22, 23)
        >>> get_filter()
        '15.0000:20.0000,25:000:30.0000'

        Return the filter for data set 3:

        >>> get_filter(3)

        """
        kwargs = {}
        if format is not None:
            kwargs["format"] = format
        if delim is not None:
            kwargs["delim"] = delim

        return self.get_data(id).get_filter(**kwargs)

    def copy_data(self, fromid: IdType, toid: IdType) -> None:
        """Copy a data set, creating a new identifier.

        After copying the data set, any changes made to the
        original data set (that is, the `fromid` identifier)
        will not be reflected in the new (the `toid`
        identifier) data set.

        Parameters
        ----------
        fromid : int or str
           The input data set.
        toid : int or str
           The output data set.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If there is no data set with a ``fromid`` identifier.

        See Also
        --------
        delete_data, load_data, set_data

        Examples
        --------

        >>> copy_data(1, 2)

        Rename the data set with identifier 2 to "orig", and
        then delete the old data set:

        >>> copy_data(2, "orig")
        >>> delete_data(2)

        """
        data = self.get_data(fromid)
        data = copy.deepcopy(data)
        self.set_data(toid, data)

    # DOC-TODO: this does not delete the source expression;
    # is this intended or a bug?
    def delete_data(self, id: Optional[IdType] = None) -> None:
        """Delete a data set by identifier.

        The data set, and any associated structures - such as the ARF
        and RMF for PHA data sets - are removed.

        Parameters
        ----------
        id : int, str, or None, optional
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
        idval = self._fix_id(id)
        self._data.pop(idval, None)

    # DOC-NOTE: also in sherpa.astro.utils
    def dataspace1d(self, start, stop, step=1, numbins=None,
                    id: Optional[IdType] = None,
                    dstype=sherpa.data.Data1DInt
                    ) -> None:
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
           ``numbins`` is set.
        numbins : int, optional
           The number of grid points. This overrides the ``step``
           setting.
        id : int, str, or None, optional
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
        The meaning of the ``stop`` parameter depends on whether it is a
        binned or unbinned data set (as set by the ``dstype``
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

        xlo, xhi, y = sherpa.utils.dataspace1d(start, stop, step=step,
                                               numbins=numbins)
        args = [xlo, xhi, y]

        if dstype is not sherpa.data.Data1DInt:
            args = [xlo, y]

        self.set_data(id, dstype('dataspace1d', *args))

    # DOC-NOTE: also in sherpa.astro.utils
    def dataspace2d(self, dims,
                    id: Optional[IdType] = None,
                    dstype=sherpa.data.Data2D
                    ) -> None:
        """Create the independent axis for a 2D data set.

        Create an "empty" two-dimensional data set by defining the
        grid on which the points are defined (the independent axis).
        The values are set to 0.

        Parameters
        ----------
        dims : sequence of 2 number
           The dimensions of the grid in ``(width,height)`` order.
        id : int, str, or None, optional
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

        >>> dataspace2d([200, 150])
        >>> image_data()

        Create a data space called "fakeimg":

        >>> dataspace2d([nx,ny], id="fakeimg")

        """
        x0, x1, y, shape = sherpa.utils.dataspace2d(dims)

        dataset = None
        if issubclass(dstype, sherpa.data.Data2DInt):
            dataset = dstype('dataspace2d', x0 - 0.5, x1 - 0.5, x0 + 0.5, x1 + 0.5,
                             y, shape)
        else:
            dataset = dstype('dataspace2d', x0, x1, y, shape)

        self.set_data(id, dataset)

    def fake(self,
             id: Optional[IdType] = None,
             method=sherpa.utils.poisson_noise
             ) -> None:
        """Simulate a data set.

        Take a data set, evaluate the model for each bin, and then use
        this value to create a data value from each bin. The default
        behavior is to use a Poisson distribution, with the model
        value as the expectation value of the distribution.

        Parameters
        ----------
        id : int, str, or None, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        method : func
           The function used to create a random realisation of
           a data set.

        See Also
        --------
        dataspace1d : Create the independent axis for a 1D data set.
        dataspace2d : Create the independent axis for a 2D data set.
        get_dep : Return the dependent axis of a data set.
        load_arrays : Create a data set from array values.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The function for the ``method`` argument accepts a single
        argument, the data values, and should return an array of the
        same shape as the input, with the data values to use.

        The function can be called on any data set, it does not need
        to have been created with `dataspace1d` or `dataspace2d`.

        Specific data set types may have their own, specialized,
        version of this function.

        Examples
        --------
        Create a random realisation of the model - a constant plus
        gaussian line - for the range x=-5 to 5.

        >>> dataspace1d(-5, 5, 0.5, dstype=Data1D)
        >>> set_source(gauss1d.gline + const1d.bgnd)
        >>> bgnd.c0 = 2
        >>> gline.fwhm = 4
        >>> gline.ampl = 5
        >>> gline.pos = 1
        >>> fake()
        >>> plot_data()
        >>> plot_model(overplot=True)

        For a 2D data set, display the simulated data, model, and
        residuals:

        >>> dataspace2d([150, 80], id='fakeimg')
        >>> set_source('fakeimg', beta2d.src + polynom2d.bg)
        >>> src.xpos, src.ypos = 75, 40
        >>> src.r0, src.alpha = 15, 2.3
        >>> src.ellip, src.theta = 0.4, 1.32
        >>> src.ampl = 100
        >>> bg.c, bg.cx1, bg.cy1 = 3, 0.4, 0.3
        >>> fake('fakeimg')
        >>> image_fit('fakeimg')

        """
        data = self.get_data(id)
        model = self.get_model(id)
        ndep = method(data.eval_model(model), rng=self.get_rng())
        self.set_dep(id, ndep)

    @staticmethod
    def _read_data(readfunc, filename, *args, **kwargs):
        _check_str_type(filename, "filename")
        return readfunc(filename, *args, **kwargs)

    # DOC-NOTE: also in sherpa.astro.utils
    # DOC-TODO: What data types are supported here?
    def unpack_arrays(self, *args):
        """Create a sherpa data object from arrays of data.

        The object returned by `unpack_arrays` can be used in a
        `set_data` call.

        Parameters
        ----------
        args : array_like
           Arrays of data. The order, and number, is determined by
           the `dstype` parameter, and listed in the `load_arrays`
           routine.
        dstype
           The data set type. The default is `Data1D` and values
           include: `Data1D`, `Data1DInt`, `Data2D`, and `Data2DInt`.
           It is expected to be derived from `sherpa.data.BaseData`.

        Returns
        -------
        instance
           The data set object matching the requested ``dstype``.

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

        >>> x = [1, 3, 7, 12]
        >>> y = [2.3, 3.2, -5.4, 12.1]
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

    def unpack_data(self, filename, ncols=2, colkeys=None,
                    dstype=sherpa.data.Data1D, sep=' ', comment='#',
                    require_floats=True):
        """Create a sherpa data object from an ASCII file.

        This function is used to read in columns from an ASCII
        file and convert them to a Sherpa data object.

        Parameters
        ----------
        filename : str
           The name of the ASCII file to read in.
        ncols : int, optional
           The number of columns to read in (the first `ncols` columns
           in the file).
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           ``None``.
        dstype : data class to use, optional
           What type of data is to be used. Supported values include
           `Data1D` (the default), `Data1DInt`, `Data2D`, and
           `Data2DInt`.
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        require_floats : bool, optional
           If `True` (the default), non-numeric data values will
           raise a `ValueError`.

        Returns
        -------
        instance
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
        load_data : Load a data set from a file.
        set_data : Set a data set.
        unpack_arrays : Create a sherpa data object from arrays of data.

        Notes
        -----
        The file reading is performed by `sherpa.io.get_ascii_data`,
        which reads in each line from the file, strips out any unsupported
        characters (replacing them by the `sep` argument), skips
        empty lines, and then identifies whether it is a comment or data
        line.

        The list of unsupported characters are: tab, new line,
        carriage return, comma, semi-colon, colon, space, and "|".

        The last comment line before the data is used to define the
        column names, splitting the line by the `sep` argument.
        If there are no comment lines then the columns are named
        starting at `col1`, `col2`, increasing up to the number of
        columns.

        Data lines are separated into columns - splitting by the
        `sep` comment - and then converted to NumPy arrays.
        If the `require_floats` argument is ``True`` then the
        column will be converted to the `sherpa.utils.SherpaFloat`
        type, with an error raised if this fails.

        An error is raised if the number of columns per row
        is not constant.

        If the `colkeys` argument is used then a case-sensitive
        match is used to determine what columns to return.

        Examples
        --------

        Create a data object from the first two columns of the file
        "src.dat" and use it to create a Sherpa data set called "src":

        >>> dat = unpack_data('src.dat')
        >>> set_data('src', dat)

        Read in the first three columns - the independent axis (x),
        the dependent variable (y), and the error on y:

        >>> dat = unpack_data('src.dat', ncols=3)

        Read in the X and Y columns from the file. The last line
        before the data must contain the column names:

        >>> dat = unpack_data('src.dat', colkeys=['X', 'Y'])

        Read in a histogram:

        >>> cols = ['XLO', 'XHI', 'Y']
        >>> idat = unpack_data('hist.dat', colkeys=cols,
        ...                    dstype=ui.Data1DInt)

        """
        return self._read_data(sherpa.io.read_data, filename, ncols, colkeys,
                               sep, dstype, comment, require_floats)

    # DOC-NOTE: also in sherpa.astro.utils
    def load_data(self, id, filename=None, ncols=2, colkeys=None,
                  dstype=sherpa.data.Data1D, sep=' ', comment='#',
                  require_floats=True) -> None:
        """Load a data set from an ASCII file.

        Parameters
        ----------
        id : int or str
           The identifier for the data set to use.
        filename : str
           The name of the ASCII file to read in.
        ncols : int, optional
           The number of columns to read in (the first ``ncols`` columns
           in the file).
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           ``None``.
        dstype : data class to use, optional
           What type of data is to be used. Supported values include
           `Data1D` (the default), `Data1DInt`, `Data2D`, and
           `Data2DInt`.
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        require_floats : bool, optional
           If ``True`` (the default), non-numeric data values will
           raise a `ValueError`.

        Raises
        ------
        ValueError
           If a column value can not be converted into a numeric value
           and the ``require_floats`` parameter is True.

        See Also
        --------
        get_data : Return the data set by identifier.
        load_arrays : Create a data set from array values.
        unpack_arrays : Create a sherpa data object from arrays of data.
        unpack_data : Create a sherpa data object from a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        See `unpack_data` for a description of the supported
        file format.

        Examples
        --------

        >>> load_data('tbl.dat')

        >>> load_data('hist.dat', dstype=Data1DInt)

        >>> cols = ['rmid', 'sur_bri', 'sur_bri_err']
        >>> load_data(2, 'profile.dat', colkeys=cols)

        """
        if filename is None:
            id, filename = filename, id
        self.set_data(id, self.unpack_data(filename, ncols=ncols,
                                           colkeys=colkeys, dstype=dstype,
                                           sep=sep, comment=comment,
                                           require_floats=require_floats))

    # DOC-NOTE: also in sherpa.astro.utils
    # DOC-TODO: rework the Data type notes section (also needed by unpack_arrays)
    # @loggable(with_id = True)
    def load_arrays(self,
                    id: IdType,
                    *args) -> None:
        """Create a data set from array values.

        Parameters
        ----------
        id : int or str
           The identifier for the data set to use.
        *args
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

        +------------+-----------------+--------------------+
        | Identifier | Required Fields |   Optional Fields  |
        +============+=================+====================+
        | Data1D     | x, y            | statistical error, |
        |            |                 | systematic error   |
        +------------+-----------------+--------------------+
        | Data1DInt  | xlo, xhi, y     | statistical error, |
        |            |                 | systematic error   |
        +------------+-----------------+--------------------+
        | Data2D     | x0, x1, y       | shape,             |
        |            |                 | statistical error, |
        |            |                 | systematic error   |
        +------------+-----------------+--------------------+
        | Data2DInt  | x0lo, x1lo,     | shape,             |
        |            | x0hi, x1hi, y   | statistical error, |
        |            |                 | systematic error   |
        +------------+-----------------+--------------------+

        The ``shape`` argument should be a tuple giving the size of
        the data ``(ny,nx)``.

        Examples
        --------

        Create a 1D data set with three points:

        >>> load_arrays(1, [10, 12, 15], [4.2, 12.1, 8.4])

        Create a 1D data set, with the identifier 'prof', from the
        arrays ``x`` (independent axis), ``y`` (dependent axis), and
        ``dy`` (statistical error on the dependent axis):

        >>> load_arrays('prof', x, y, dy)

        Explicitly define the type of the data set:

        >>> load_arrays('prof', x, y, dy, Data1D)

        Data set 1 is a histogram, where the bins cover the range
        1-3, 3-5, and 5-7 with values 4, 5, and 9 respectively.

        >>> load_arrays(1, [1, 3, 5], [3, 5, 7], [4, 5, 9], Data1DInt)

        """
        self.set_data(id, self.unpack_arrays(*args))

    def _save_type(self, objtype: str,
                   id, filename, **kwargs) -> None:
        if filename is None:
            id, filename = filename, id

        _check_str_type(filename, "filename")

        d = self.get_data(id)

        args = None
        fields = None

        if type(d) in (sherpa.data.Data2D, sherpa.data.Data2DInt):
            if objtype == 'delchi':
                # raise AttributeErr('badfunc', "save_delchi()", "images")
                raise AttributeError("save_delchi() can not be used " +
                                     "with 2D datasets")

            funcname = f"get_{objtype}_image"
            imgtype = getattr(self, funcname)

            obj = imgtype(id)
            args = [obj.y]
            fields = [str(objtype).upper()]

        else:
            funcname = f"get_{objtype}_plot"
            plottype = getattr(self, funcname)

            obj = plottype(id)
            args = [obj.x, obj.y]
            fields = ["X", str(objtype).upper()]

        sherpa.io.write_arrays(filename, args, fields, **kwargs)

    # DOC-NOTE: also in sherpa.astro.utils with a different interface
    def save_arrays(self, filename, args, fields=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'
                    ) -> None:
        """Write a list of arrays to an ASCII file.

        Parameters
        ----------
        filename : str
           The name of the file to write the array to.
        args : array of arrays
           The arrays to write out.
        fields : array of str
           The column names (should match the size of `args`).
        clobber : bool, optional
           If `filename` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        linebreak : str, optional
           Indicate a new line. The default is ``'\\n'``.
        format : str, optional
           The format used to write out the numeric values. The
           default is ``'%g%'``.

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

        See Also
        --------
        save_data : Save the data to a file.

        Examples
        --------

        Write the x and y columns from the default data set to the
        file 'src.dat':

        >>> x = get_indep()
        >>> y = get_dep()
        >>> save_arrays('src.dat', [x, y])

        Use the column names "r" and "surbri" for the columns:

        >>> save_arrays('prof.txt', [x, y], fields=["r", "surbri"],
        ...             clobber=True)

        """
        clobber = sherpa.utils.bool_cast(clobber)
        _check_str_type(filename, "filename")
        sherpa.io.write_arrays(filename, args, fields, sep, comment, clobber,
                               linebreak, format)

    # DOC-NOTE: also in sherpa.utils with a different interface
    def save_source(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'
                    ) -> None:
        """Save the model values to a file.

        The model is evaluated on the grid of the data set, but does
        not include any instrument response (such as a PSF).

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
           overwritten (``True``) or if it raises an exception (``False``,
           the default setting).
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        linebreak : str, optional
           Indicate a new line. The default is ``'\\n'``.
        format : str, optional
           The format used to write out the numeric values. The
           default is ``'%g%'``.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If no model has been set for this data set.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

        See Also
        --------
        save_data : Save the data to a file.
        save_model : Save the model values to a file.
        set_full_model : Define the convolved model expression for a data set.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns ``X`` and ``SOURCE`` (for 1D
        data). The residuals array respects any filter setting for the
        data set.

        Examples
        --------

        Write the model to the file "model.dat":

        >>> save_source('model.dat')

        Write the model from the data set 'jet' to the
        file "jet.mdl":

        >>> save_source('jet', "jet.mdl", clobber=True)

        """
        clobber = sherpa.utils.bool_cast(clobber)
        self._save_type('source', id, filename, clobber=clobber, sep=sep,
                        comment=comment, linebreak=linebreak, format=format)

    # DOC-NOTE: also in sherpa.utils with a different interface
    def save_model(self, id, filename=None, clobber=False, sep=' ',
                   comment='#', linebreak='\n', format='%g'
                   ) -> None:
        """Save the model values to a file.

        The model is evaluated on the grid of the data set, including
        any instrument response (such as a PSF).

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to.
        clobber : bool, optional
           If `filename` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        linebreak : str, optional
           Indicate a new line. The default is ``'\\n'``.
        format : str, optional
           The format used to write out the numeric values. The
           default is ``'%g%'``.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If no model has been set for this data set.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

        See Also
        --------
        save_data : Save the data to a file.
        save_source : Save the model values to a file.
        set_model : Set the source model expression for a data set.
        set_full_model : Define the convolved model expression for a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns ``X`` and ``MODEL`` (for 1D
        data). The residuals array respects any filter setting for the
        data set.

        Examples
        --------

        Write the model to the file "model.dat":

        >>> save_model('model.dat')

        Write the model from the data set 'jet' to the
        file "jet.mdl":

        >>> save_model('jet', "jet.mdl", clobber=True)

        """
        clobber = sherpa.utils.bool_cast(clobber)
        self._save_type('model', id, filename, clobber=clobber, sep=sep,
                        comment=comment, linebreak=linebreak, format=format)

    # DOC-NOTE: also in sherpa.utils with a different interface
    def save_resid(self, id, filename=None, clobber=False, sep=' ',
                   comment='#', linebreak='\n', format='%g'
                   ) -> None:
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
           If `filename` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        linebreak : str, optional
           Indicate a new line. The default is ``'\\n'``.
        format : str, optional
           The format used to write out the numeric values. The
           default is ``'%g%'``.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If no model has been set for this data set.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

        See Also
        --------
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

        The output file contains the columns ``X`` and ``RESID``. The
        residuals array respects any filter setting for the data set.

        Examples
        --------

        Write the residuals to the file "resid.dat":

        >>> save_resid('resid.dat')

        Write the residuals from the data set 'jet' to the
        file "resid.dat":

        >>> save_resid('jet', "resid.dat", clobber=True)

        """
        clobber = sherpa.utils.bool_cast(clobber)
        self._save_type('resid', id, filename, clobber=clobber, sep=sep,
                        comment=comment, linebreak=linebreak, format=format)

    # DOC-NOTE: also in sherpa.utils with a different interface
    def save_delchi(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'
                    ) -> None:
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
           If `filename` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        linebreak : str, optional
           Indicate a new line. The default is ``'\\n'``.
        format : str, optional
           The format used to write out the numeric values. The
           default is ``'%g%'``.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If no model has been set for this data set.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

        See Also
        --------
        save_data : Save the data to a file.
        save_resid : Save the residuals (data-model) to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns ``X`` and ``DELCHI``. The
        residuals array respects any filter setting for the data set.

        Examples
        --------

        Write the residuals to the file "delchi.dat":

        >>> save_delchi('delchi.dat')

        Write the residuals from the data set 'jet' to the
        file "delchi.dat":

        >>> save_resid('jet', "delchi.dat", clobber=True)

        """
        clobber = sherpa.utils.bool_cast(clobber)
        self._save_type('delchi', id, filename, clobber=clobber, sep=sep,
                        comment=comment, linebreak=linebreak, format=format)

    # DOC-NOTE: also in sherpa.astro.utils
    def save_data(self, id, filename=None, fields=None, sep=' ', comment='#',
                  clobber=False, linebreak='\n', format='%g'
                  ) -> None:
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
           The attributes of the data set to write out. If ``None``,
           write out all the columns.
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        clobber : bool, optional
           If `filename` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).
        linebreak : str, optional
           Indicate a new line. The default is ``'\\n'``.
        format : str, optional
           The format used to write out the numeric values. The
           default is ``'%g%'``.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If there is no matching data set.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

        See Also
        --------
        save_arrays : Write a list of arrays to a file.
        save_delchi : Save the ratio of residuals (data-model) to error to a file.
        save_error : Save the errors to a file.
        save_filter : Save the filter array to a file.
        save_resid : Save the residuals (data-model) to a file.
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
        ...           fields=['x', 'y', 'staterror'])

        """
        clobber = sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id

        _check_str_type(filename, "filename")
        sherpa.io.write_data(filename, self.get_data(id), fields, sep,
                             comment, clobber, linebreak, format)

    # DOC-NOTE: also in sherpa.astro.utils with a different interface
    def save_filter(self, id, filename=None, clobber=False, sep=' ',
                    comment='#', linebreak='\n', format='%g'
                    ) -> None:
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
           If `filename` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        linebreak : str, optional
           Indicate a new line. The default is ``'\\n'``.
        format : str, optional
           The format used to write out the numeric values. The
           default is ``'%g%'``.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set has not been filtered.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

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

        The output file contains the columns ``X`` and ``FILTER``.

        Examples
        --------

        Write the filter from the default data set as an ASCII file:

        >>> save_filter('filt.dat')

        """
        clobber = sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id

        _check_str_type(filename, "filename")
        d = self.get_data(id)
        id = self._fix_id(id)

        if d.mask is False:
            raise sherpa.utils.err.DataErr('notmask')
        if not np.iterable(d.mask):
            raise sherpa.utils.err.DataErr('nomask', id)

        x = d.get_indep(filter=False)[0]
        mask = np.asarray(d.mask, int)
        self.save_arrays(filename, [x, mask], fields=['X', 'FILTER'],
                         clobber=clobber, sep=sep, comment=comment,
                         linebreak=linebreak, format=format)

    # DOC-NOTE: also in sherpa.astro.utils with a different interface
    def save_staterror(self, id, filename=None, clobber=False, sep=' ',
                       comment='#', linebreak='\n', format='%g'
                       ) -> None:
        """Save the statistical errors to a file.

        If the statistical errors have not been set explicitly, then
        the values calculated by the statistic - such as ``chi2gehrels``
        or ``chi2datavar`` - will be used.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to.
        clobber : bool, optional
           If `filename` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        linebreak : str, optional
           Indicate a new line. The default is ``'\\n'``.
        format : str, optional
           The format used to write out the numeric values. The
           default is ``'%g%'``.

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

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

        The output file contains the columns ``X`` and ``STAT_ERR``.

        Examples
        --------

        Write out the statistical errors from the default data set to the
        file 'errs.dat'.

        >>> save_staterror('errs.dat')

        Over-write the file it it already exists, and take the data
        from the data set "jet":

        >>> save_staterror('jet', 'err.out', clobber=True)

        """
        clobber = sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id

        _check_str_type(filename, "filename")
        x = self.get_data(id).get_indep(filter=False)[0]
        err = self.get_staterror(id, filter=False)
        self.save_arrays(filename, [x, err], fields=['X', 'STAT_ERR'],
                         clobber=clobber, sep=sep, comment=comment,
                         linebreak=linebreak, format=format)

    # DOC-NOTE: also in sherpa.astro.utils with a different interface
    def save_syserror(self, id, filename=None, clobber=False, sep=' ',
                      comment='#', linebreak='\n', format='%g'
                      ) -> None:
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
           If `filename` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        linebreak : str, optional
           Indicate a new line. The default is ``'\\n'``.
        format : str, optional
           The format used to write out the numeric values. The
           default is ``'%g%'``.

        Raises
        ------
        sherpa.utils.err.IOErr
           If the data set does not contain any systematic errors.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

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

        The output file contains the columns ``X`` and ``SYS_ERR``.

        Examples
        --------

        Write out the systematic errors from the default data set to the
        file 'errs.dat'.

        >>> save_syserror('errs.dat')

        Over-write the file it it already exists, and take the data
        from the data set "jet":

        >>> save_syserror('jet', 'err.out', clobber=True)

        """
        clobber = sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id

        _check_str_type(filename, "filename")
        x = self.get_data(id).get_indep(filter=False)[0]
        err = self.get_syserror(id, filter=False)
        self.save_arrays(filename, [x, err], fields=['X', 'SYS_ERR'],
                         clobber=clobber, sep=sep, comment=comment,
                         linebreak=linebreak, format=format)

    # DOC-NOTE: also in sherpa.astro.utils with a different interface
    def save_error(self, id, filename=None, clobber=False, sep=' ',
                   comment='#', linebreak='\n', format='%g'
                   ) -> None:
        """Save the errors to a file.

        The total errors for a data set are the quadrature combination
        of the statistical and systematic errors. The systematic
        errors can be 0. If the statistical errors have not been set
        explicitly, then the values calculated by the statistic - such
        as ``chi2gehrels`` or ``chi2datavar`` - will be used.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to.
        clobber : bool, optional
           If `filename` is not ``None``, then this flag controls
           whether an existing file can be overwritten (``True``)
           or if it raises an exception (``False``, the default
           setting).
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        linebreak : str, optional
           Indicate a new line. The default is ``'\\n'``.
        format : str, optional
           The format used to write out the numeric values. The
           default is ``'%g%'``.

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is ``False``.

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

        The output file contains the columns ``X`` and ``ERR``.

        Examples
        --------

        Write out the errors from the default data set to the file
        'errs.dat'.

        >>> save_error('errs.dat')

        Over-write the file it it already exists, and take the data
        from the data set "jet":

        >>> save_error('jet', 'err.out', clobber=True)

        """
        clobber = sherpa.utils.bool_cast(clobber)
        if filename is None:
            id, filename = filename, id

        _check_str_type(filename, "filename")
        x = self.get_data(id).get_indep(filter=False)[0]
        err = self.get_error(id, filter=False)
        self.save_arrays(filename, [x, err], fields=['X', 'ERR'],
                         clobber=clobber, sep=sep, comment=comment,
                         linebreak=linebreak, format=format)

    def _notice_expr(self,
                     expr: Optional[str] = None,
                     **kwargs) -> None:
        ids = self.list_data_ids()
        for vals in sherpa.utils.parse_expr(expr):
            self.notice_id(ids, *vals, **kwargs)

    def _notice_expr_id(self,
                        ids: Union[IdType, Sequence[IdType]],
                        expr: Optional[str] = None,
                        **kwargs) -> None:
        for vals in sherpa.utils.parse_expr(expr):
            self.notice_id(ids, *vals, **kwargs)

    # DOC-NOTE: inclusion of bkg_id is technically wrong, as it
    # should only be in the sherpa.astro.ui version, but it is not
    # worth creating a copy of the routine just for this.
    #
    def notice(self, lo=None, hi=None, **kwargs) -> None:
        """Include data in the fit.

        Select one or more ranges of data to include by filtering on
        the independent axis value. The filter is applied to all data
        sets.

        .. versionchanged:: 4.15.0
           The change in the filter is now reported for each dataset.

        .. versionchanged:: 4.14.0
           Integrated data sets - so Data1DInt and DataPHA when using
           energy or wavelengths - now ensure that the `hi` argument
           is exclusive and better handling of the `lo` argument when it
           matches a bin edge. This can result in the same filter selecting
           a smaller number of bins than in earlier versions of Sherpa.

        Parameters
        ----------
        lo : number or str, optional
           The lower bound of the filter (when a number) or a string
           expression listing ranges in the form ``a:b``, with
           multiple ranges allowed, where the ranges are separated by
           a ``,``. The term ``:b`` means include everything up to
           ``b`` (an exclusive limit for integrated datasets), and
           ``a:`` means include everything that is higher than, or
           equal to, ``a``.
        hi : number, optional
           The upper bound of the filter when ``lo`` is not a string.
        bkg_id : int or str, optional
           The filter will be applied to the associated background
           component of the data set if ``bkg_id`` is set. Only PHA
           data sets support this option; if not given, then the
           filter is applied to all background components as well
           as the source data.

        See Also
        --------
        notice_id : Include data for a data set.
        sherpa.astro.ui.notice2d : Include a spatial region in an image.
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

        For binned data sets, the bin is included if the noticed range
        falls anywhere within the bin, but excluding the ``hi`` value
        (except for PHA data sets when using ``channel`` units).

        The units used depend on the ``analysis`` setting of the data
        set, if appropriate.

        To filter a 2D data set by a shape use `notice2d`.

        The report of the change in the filter expression can be
        controlled with the `SherpaVerbosity` context manager, as
        shown in the examples below.

        Examples
        --------

        Since the `notice` call is applied to an un-filtered
        data set, the filter chooses only those points that lie
        within the range 12 <= X <= 18.

        >>> load_arrays(1, [10, 15, 20, 30], [5, 10, 7, 13])
        >>> notice(12, 28)
        dataset 1: 10:30 -> 15:20 x
        >>> get_dep(filter=True)
        array([10,  7])

        As no limits are given, the whole data set is included:

        >>> notice()
        dataset 1: 15:20 -> 10:30 x
        >>> get_dep(filter=True)
        array([ 5, 10,  7, 13])

        The `ignore` call excludes the first two points, but the
        `notice` call adds back in the second point:

        >>> ignore(hi=17)
        dataset 1: 10:30 -> 20:30 x
        >>> notice(12, 16)
        dataset 1: 20:30 -> 15:30 x
        >>> get_dep(filter=True)
        array([10,  7, 13])

        Only include data points in the range 8<=X<=12 and 18<=X=22:

        >>> ignore()
        dataset 1: 15:30 x -> no data
        >>> notice("8:12, 18:22")
        dataset 1: no data -> 10 x
        dataset 1: 10 -> 10,20 x
        >>> get_dep(filter=True)
        array([5, 7])

        The messages from `notice` and `ignore` use the standard
        Sherpa logging infrastructure, and so can be ignored by using
        `SherpaVerbosity`:

        >>> from sherpa.utils.logging import SherpaVerbosity
        >>> with SherpaVerbosity("WARN"):
        ...     notice()
        ...

        """
        if len(self._data) == 0:
            raise IdentifierErr("nodatasets")

        # TODO: do we still expect to get bytes here?
        if lo is not None and isinstance(lo, (str, np.bytes_)):
            return self._notice_expr(lo, **kwargs)

        # Jump through the data sets in "order".
        #
        notice_data_range(self.get_data, self.list_data_ids(), lo, hi, kwargs)

    # DOC-NOTE: inclusion of bkg_id is technically wrong, as it
    # should only be in the sherpa.astro.ui version, but it is not
    # worth creating a copy of the routine just for this.
    #
    def ignore(self, lo=None, hi=None, **kwargs) -> None:
        """Exclude data from the fit.

        Select one or more ranges of data to exclude by filtering on
        the independent axis value. The filter is applied to all data
        sets.

        .. versionchanged:: 4.15.0
           The change in the filter is now reported for each dataset.

        .. versionchanged:: 4.14.0
           Integrated data sets - so Data1DInt and DataPHA when using
           energy or wavelengths - now ensure that the `hi` argument
           is exclusive and better handling of the `lo` argument when it
           matches a bin edge. This can result in the same filter selecting
           a smaller number of bins than in earlier versions of Sherpa.

        Parameters
        ----------
        lo : number or str, optional
           The lower bound of the filter (when a number) or a string
           expression listing ranges in the form ``a:b``, with
           multiple ranges allowed, where the ranges are separated by
           a ``,``. The term ``:b`` means exclude everything up to
           ``b`` (an exclusive limit for integrated datasets), and
           ``a:`` means exclude everything that is higher than, or
           equal to, ``a``.
        hi : number, optional
           The upper bound of the filter when ``lo`` is not a string.
        bkg_id : int or str, optional
           The filter will be applied to the associated background
           component of the data set if ``bkg_id`` is set. Only PHA
           data sets support this option; if not given, then the
           filter is applied to all background components as well
           as the source data.

        See Also
        --------
        ignore_id : Exclude data from the fit for a data set.
        sherpa.astro.ui.ignore2d : Exclude a spatial region from an image.
        notice : Include data in the fit.
        show_filter : Show any filters applied to a data set.

        Notes
        -----
        The order of `ignore` and `notice` calls is important, and the
        results are a union, rather than intersection, of the
        combination.

        For binned data sets, the bin is excluded if the ignored
        range falls anywhere within the bin.

        The units used depend on the ``analysis`` setting of the data
        set, if appropriate.

        To filter a 2D data set by a shape use `ignore2d`.

        The report of the change in the filter expression can be
        controlled with the `SherpaVerbosity` context manager, as
        shown in the examples below.

        Examples
        --------

        Ignore all data points with an X value (the independent axis)
        between 12 and 18. For this one-dimensional data set, this
        means that the second bin is ignored:

        >>> load_arrays(1, [10, 15, 20, 30], [5, 10, 7, 13])
        >>> ignore(12, 18)
        dataset 1: 10:30 -> 10,20:30 x
        >>> get_dep(filter=True)
        array([ 5,  7, 13])

        Filtering X values that are 25 or larger means that the last
        point is also ignored:

        >>> ignore(lo=25)
        dataset 1: 10,20:30 -> 10,20 x
        >>> get_dep(filter=True)
        array([ 5,  7])

        The `notice` call removes the previous filter, and then a
        multi-range filter is applied to exclude values between 8 and
        12 and 18 and 22:

        >>> notice()
        dataset 1: 10,20 -> 10:30 x
        >>> ignore("8:12, 18:22")
        dataset 1: 10:30 -> 15:30 x
        dataset 1: 15:30 -> 15,30 x
        >>> get_dep(filter=True)
        array([10, 13])

        The `SherpaVerbosity` context manager can be used to hide the
        screen output:

        >>> from sherpa.utils.logging import SherpaVerbosity
        >>> with SherpaVerbosity("WARN"):
        ...     ignore(hi=12)
        ...

        """
        kwargs['ignore'] = True
        self.notice(lo, hi, **kwargs)

    # DOC-NOTE: inclusion of bkg_id is technically wrong, as it
    # should only be in the sherpa.astro.ui version, but it is not
    # worth creating a copy of the routine just for this.
    #
    def notice_id(self,
                  ids: Union[IdType, Sequence[IdType]],
                  lo=None, hi=None, **kwargs
                  ) -> None:
        """Include data from the fit for a data set.

        Select one or more ranges of data to include by filtering on
        the independent axis value. The filter is applied to the
        given data set, or data sets.

        .. versionchanged:: 4.15.0
           The change in the filter is now reported for the dataset.

        .. versionchanged:: 4.14.0
           Integrated data sets - so Data1DInt and DataPHA when using
           energy or wavelengths - now ensure that the `hi` argument
           is exclusive and better handling of the `lo` argument when it
           matches a bin edge. This can result in the same filter selecting
           a smaller number of bins than in earlier versions of Sherpa.

        Parameters
        ----------
        ids : int or str, or array of int or str
           The data set, or sets, to use.
        lo : number or str, optional
           The lower bound of the filter (when a number) or a string
           expression listing ranges in the form ``a:b``, with
           multiple ranges allowed, where the ranges are separated by
           a ``,``. The term ``:b`` means include everything up to
           ``b`` (an exclusive limit for integrated datasets), and
           ``a:`` means include everything that is higher than, or
           equal to, ``a``.
        hi : number, optional
           The upper bound of the filter when ``lo`` is not a string.
        bkg_id : int or str, optional
           The filter will be applied to the associated background
           component of the data set if ``bkg_id`` is set. Only PHA
           data sets support this option; if not given, then the
           filter is applied to all background components as well
           as the source data.

        See Also
        --------
        ignore_id : Exclude data from the fit for a data set.
        sherpa.astro.ui.ignore2d : Exclude a spatial region from an image.
        notice : Include data in the fit.
        show_filter : Show any filters applied to a data set.

        Notes
        -----
        The order of `ignore` and `notice` calls is important.

        The units used depend on the ``analysis`` setting of the data
        set, if appropriate.

        To filter a 2D data set by a shape use `ignore2d`.

        The report of the change in the filter expression can be
        controlled with the `SherpaVerbosity` context manager, as
        shown in the examples below.

        Examples
        --------

        Include all data points with an X value (the independent axis)
        between 12 and 18 for data set 1:

        >>> notice_id(1, 12, 18)
        dataset 1: 10:30 -> 15:20 x

        Include the range 0.5 to 7, for data sets 1, 2, and 3 (the
        screen output will depend on the existing data and filters
        applied to them):

        >>> notice_id([1, 2, 3], 0.5, 7)
        dataset 1: 0.00146:14.9504 -> 0.4818:9.0374 Energy (keV)
        dataset 2: 0.00146:14.9504 -> 0.4964:13.6072 Energy (keV)
        dataset 3: 0.00146:14.9504 -> 0.4234:9.3878 Energy (keV)

        Apply the filter 0.5 to 2 and 2.2 to 7 to the data sets "core"
        and "jet", and hide the screen output:

        >>> from sherpa.utils.logging import SherpaVerbosity
        >>> with SherpaVerbsity("WARN"):
        ...     notice_id(["core", "jet"], "0.5:2, 2.2:7")
        ...

        """
        if ids is None:
            raise ArgumentTypeErr('badarg', 'ids',
                                  'an identifier or list of identifiers')

        if self._valid_id(ids):
            idvals = (ids,)
        else:
            try:
                idvals = tuple(ids)
            except TypeError:
                raise ArgumentTypeErr('badarg', 'ids',
                                      'an identifier or list of identifiers') from None

        # TODO: do we still expect to get bytes here?
        if lo is not None and isinstance(lo, (str, np.bytes_)):
            return self._notice_expr_id(idvals, lo, **kwargs)

        # Unlike notice() we do not sort the id list as this
        # was set by the user.
        #
        notice_data_range(self.get_data, idvals, lo, hi, kwargs)

    # DOC-NOTE: inclusion of bkg_id is technically wrong, as it
    # should only be in the sherpa.astro.ui version, but it is not
    # worth creating a copy of the routine just for this.
    #
    def ignore_id(self,
                  ids: Union[IdType, Sequence[IdType]],
                  lo=None, hi=None, **kwargs
                  ) -> None:
        """Exclude data from the fit for a data set.

        Select one or more ranges of data to exclude by filtering on
        the independent axis value. The filter is applied to the given
        data set, or sets.

        .. versionchanged:: 4.15.0
           The change in the filter is now reported for the dataset.

        .. versionchanged:: 4.14.0
           Integrated data sets - so Data1DInt and DataPHA when using
           energy or wavelengths - now ensure that the `hi` argument
           is exclusive and better handling of the `lo` argument when it
           matches a bin edge. This can result in the same filter selecting
           a smaller number of bins than in earlier versions of Sherpa.

        Parameters
        ----------
        ids : int or str, or array of int or str
           The data set, or sets, to use.
        lo : number or str, optional
           The lower bound of the filter (when a number) or a string
           expression listing ranges in the form ``a:b``, with
           multiple ranges allowed, where the ranges are separated by
           a ``,``. The term ``:b`` means exclude everything up to
           ``b`` (an exclusive limit for integrated datasets), and
           ``a:`` means exclude everything that is higher than, or
           equal to, ``a``.
        hi : number, optional
           The upper bound of the filter when ``lo`` is not a string.
        bkg_id : int or str, optional
           The filter will be applied to the associated background
           component of the data set if ``bkg_id`` is set. Only PHA
           data sets support this option; if not given, then the
           filter is applied to all background components as well
           as the source data.

        See Also
        --------
        ignore : Exclude data from the fit.
        sherpa.astro.ui.ignore2d : Exclude a spatial region from an image.
        notice_id : Include data from the fit for a data set.
        show_filter : Show any filters applied to a data set.

        Notes
        -----
        The order of `ignore` and `notice` calls is important.

        The units used depend on the ``analysis`` setting of the data
        set, if appropriate.

        To filter a 2D data set by a shape use `ignore2d`.

        Examples
        --------

        Ignore all data points with an X value (the independent axis)
        between 12 and 18 for data set 1:

        >>> ignore_id(1, 12, 18)
        dataset 1: 10:30 -> 10,20:30 x

        Ignore the range up to 0.5 and 7 and above, for data sets 1,
        2, and 3:

        >>> ignore_id([1, 2, 3], hi=0.5)
        dataset 1: 0.00146:14.9504 -> 0.584:14.9504 Energy (keV)
        dataset 2: 0.00146:14.9504 -> 0.6424:14.9504 Energy (keV)
        dataset 3: 0.00146:14.9504 -> 0.511:14.9504 Energy (keV)
        >>> ignore_id([1, 2, 3], lo=7)
        dataset 1: 0.584:14.9504 -> 0.584:4.4384 Energy (keV)
        dataset 2: 0.6424:14.9504 -> 0.6424:5.1392 Energy (keV)
        dataset 3: 0.511:14.9504 -> 0.511:4.526 Energy (keV)

        Apply the same filter as the previous example, but to data
        sets "core" and "jet", and hide the screen output:

        >>> from sherpa.utils.logging import SherpaVerbosity
        >>> with SherpaVerbsity("WARN"):
        ...     ignore_id(["core", "jet"], ":0.5,7:")
        ...

        """
        kwargs['ignore'] = True
        self.notice_id(ids, lo, hi, **kwargs)

    ###########################################################################
    # Models
    ###########################################################################

    def paramprompt(self, val=False) -> None:
        """Should the user be asked for the parameter values when creating a model?

        When `val` is ``True``, calls to `set_model` will cause the user
        to be prompted for each parameter in the expression.  The
        prompt includes the parameter name and default value, in ``[]``:
        the valid responses are

        - return  which accepts the default
        - value   which changes the parameter value
        - value, min  which changes the value and the minimum value
        - value, min, max  which changes the value, minimum, and
          maximum values

        The ``value``, ``min``, and ``max`` components are optional, so
        ",-5" will use the default parameter value and set its minimum
        to -5, while "2,,10" will change the parameter value to 2 and
        its maximum to 10, but leave the minimum at its default. If
        any value is invalid then the parameter is re-prompted.

        Parameters
        ----------
        val : bool, optional
           If ``True``, the user will be prompted to enter each
           parameter value, including support for changing the minimum
           and maximum values, when a model component is created. The
           default is ``False``.

        See Also
        --------
        set_model : Set the source model expression for a data set.
        set_par : Set the value, limits, or behavior of a model parameter.
        show_model : Display the model expression used to fit a data set.

        Notes
        -----
        Setting this to ``True`` only makes sense in an interactive
        environment.  It is designed to be similar to the parameter
        prompting provided by `XSPEC <https://heasarc.gsfc.nasa.gov/xanadu/xspec/>`_.

        Examples
        --------

        In the following, the default parameter settings are accepted
        for the ``pl.gamma`` parameter, the starting values for the
        ``pl.ref`` and ``gline.pos`` values are changed, the starting
        value and ranges of both the ``pl.ampl`` and ``gline.ampl``
        parameters are set, and the ``gline.fwhm`` parameter is set to
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
                         baselist=(sherpa.models.ArithmeticModel,)
                         ) -> None:
        if not isinstance(baselist, tuple):
            baselist = (baselist,)

        for name in module.__all__:
            cls = getattr(module, name)

            for base in baselist:
                if is_subclass(cls, base):
                    break
            else:
                continue

            name = name.lower()
            self._model_types[name] = ModelWrapper(self, cls)
            self._model_globals.update(self._model_types)

    def add_model(self, modelclass, args=(), kwargs={}) -> None:
        """Create a user-defined model class.

        Create a model from a class. The name of the class can then be
        used to create model components - e.g.  with `set_model` or
        `create_model_component` - as with any existing Sherpa model.

        Parameters
        ----------
        modelclass
           A class derived from `sherpa.models.model.ArithmeticModel`. This
           class defines the functional form and the parameters of the
           model.
        args
           Arguments for the class constructor.
        kwargs
           Keyword arguments for the class constructor.

        See Also
        --------
        create_model_component : Create a model component.
        list_models : List the available model types.
        load_table_model : Load tabular data and use it as a model component.
        load_user_model : Create a user-defined model.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The `load_user_model` function is designed to make it easy to
        add a model, but the interface is not the same as the existing
        models (such as having to call both `load_user_model` and
        `add_user_pars` for each new instance).  The `add_model`
        function is used to add a model as a Python class, which is
        more work to set up, but then acts the same way as the
        existing models.

        Examples
        --------

        The following example creates a model type called "mygauss1d"
        which will behave exactly the same as the existing "gauss1d"
        model.  Normally the class used with `add_model` would add new
        functionality.

        >>> from sherpa.models import Gauss1D
        >>> class MyGauss1D(Gauss1D):
        ...     pass
        ...
        >>> add_model(MyGauss1D)
        >>> set_source(mygauss1d.g1 + mygauss1d.g2)

        """
        name = modelclass.__name__.lower()

        if not is_subclass(modelclass, sherpa.models.ArithmeticModel):
            raise TypeError(f"model class '{name}' is not a derived class "
                            "from sherpa.models.ArithmeticModel")

        self._model_types[name] = ModelWrapper(self, modelclass, args, kwargs)
        self._model_globals.update(self._model_types)
        _assign_obj_to_main(name, self._model_types[name])

    #
    # Model components
    #

    def get_model_autoassign_func(self
                                  ) -> Callable[[str, Model], None]:
        """Return the method used to create model component identifiers.

        Provides access to the function which is used by
        `create_model_component` and when creating model components
        directly to add an identifier in the current Python namespace.

        Returns
        -------
        func
           The model function set by `set_model_autoassign_func`.

        See Also
        --------
        create_model_component : Create a model component.
        set_model : Set the source model expression for a data set.
        set_model_autoassign_func : Set the method used to create model component identifiers.

        """
        return self._model_autoassign_func

    # DOC-TODO: what does func=None mean? If you try None then it
    # fails with AttributeError: 'Session' object attribute
    # '_model_autoassign_func' cannot be replaced with a non-callable attribute
    def set_model_autoassign_func(self,
                                  func: Optional[Callable[[str, Model], None]] = None
                                  ) -> None:
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
            raise ArgumentTypeErr('badarg', 'func',
                                  'a function or other callable object')

        self._model_autoassign_func = func

    def list_models(self, show: str = "all") -> list[str]:
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

        >>> list_models('2d')
        ['beta2d',
         'box2d',
         'const2d',
         'delta2d',
         'devaucouleurs2d',
         'disk2d',
         'gauss2d',
         'lorentz2d',
         'normgauss2d',
         'polynom2d',
         'scale2d',
         'sersic2d',
         'shell2d',
         'sigmagauss2d']

        """
        keys = list(self._model_types.keys())
        keys.sort()

        show = show.strip().lower()
        for item in ["arfmodel", "convolutionmodel", "multiresponsesummodel", "pileuprmfmodel",
                     "jdpileup", "rmfmodel", "tablemodel", "usermodel", "psfmodel"]:
            if item in keys:
                keys.remove(item)

        select = None
        if show.startswith("xspec"):
            select = lambda x: x.startswith('xs')
        elif show.startswith("1d"):
            select = lambda x: (not x.startswith('xs')) and (not x.endswith("2d"))
        elif show.startswith("2d"):
            select = lambda x: x.endswith('2d')

        if select is None:
            return keys

        return list(filter(select, keys))

    def list_model_components(self) -> list[str]:
        """List the names of all the model components.

        Models are created either directly - by using the form
        ``mname.mid``, where ``mname`` is the name of the model, such as
        ``gauss1d``, and ``mid`` is the name of the component - or with
        the `create_model_component` function, which accepts ``mname``
        and ``mid`` as separate arguments. This function returns all the
        ``mid`` values that have been created.

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

        The ``gal`` and ``pl`` model components are created - as versions
        of the ``xsphabs`` and ``powlaw1d`` model types - which means
        that the list of model components returned as ``mids`` will
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
        keys = list(self._model_components.keys())
        keys.sort()
        return keys

    def _add_model_component(self, cmpt: Model) -> None:
        # If model component name is a model type name
        # or session function name, don't create it, raise
        # warning
        if cmpt.name in self._model_types:
            modeltype = cmpt.name
            del cmpt
            raise IdentifierErr('badidmodel', modeltype)

        if cmpt.name.lower() in _builtin_symbols_:
            modeltype = cmpt.name
            del cmpt
            raise IdentifierErr('badidnative', modeltype)

        self._model_components[cmpt.name] = cmpt
        if self._model_autoassign_func is not None:
            self._model_autoassign_func(cmpt.name, cmpt)

    @overload
    def _get_model_component(self,
                             name: str,
                             require: Literal[False]
                             ) -> Optional[Model]:
        ...

    @overload
    def _get_model_component(self,
                             name: str,
                             require: Literal[True]
                             ) -> Model:
        ...

    def _get_model_component(self, name, require=False):
        """Access the model component by name.

        Parameters
        ----------
        name : str
            The name of the model component.
        require : bool, optional
            If `True` then the name must exist.

        Returns
        -------
        component : sherpa.models.model.Model instance or None
            `None` is returned if the name does not exist
            and the require argument is `False`.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
            If the name does not exist and require is `True`.

        """
        cmpt = self._model_components.get(name)
        require = sherpa.utils.bool_cast(require)
        if require and (cmpt is None):
            raise IdentifierErr('nomodelcmpt', name)
        return cmpt

    # DOC-TODO: can send in a model variable, but this is just the
    # identity function, so not worth documenting
    def get_model_component(self, name: str) -> Model:
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
           If there is no model component with the given ``name``.

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
        the ``username`` component that is used here to access the
        instance.

        Examples
        --------

        When a model component is created, a variable is created that
        contains the model instance. The instance can also be returned
        with `get_model_component`, which can then be queried or used
        to change the model settings:

        >>> create_model_component('gauss1d', 'gline')
        >>> gmodel = get_model_component('gline')
        >>> gmodel.name
        'gauss1d.gline'
        >>> print([p.name for p in gmodel.pars])
        ['fwhm', 'pos', 'ampl']
        >>> gmodel.fwhm.val = 12.2
        >>> gmodel.fwhm.freeze()

        """
        # If user mistakenly passes an actual model reference,
        # just return the reference
        if isinstance(name, Model):
            return name

        _check_str_type(name, "name")
        return self._get_model_component(name, require=True)

    def create_model_component(self, typename=None, name=None):
        """Create a model component.

        Model components created by this function are set to their
        default values. Components can also be created directly using
        the syntax ``typename.name``, such as in calls to `set_model`
        and `set_source` (unless you have called `set_model_autoassign_func`
        to change the default model auto-assignment setting).

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

        Returns
        -------
        model : the sherpa.models.Model object created

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
        is true using direct access, such as ``mname.parname``, or with
        `set_par`).

        Examples
        --------

        Create an instance of the ``powlaw1d`` model called ``pl``,
        and then freeze its ``gamma`` parameter to 2.6.

        >>> create_model_component("powlaw1d", "pl")
        >>> pl.gamma = 2.6
        >>> freeze(pl.gamma)

        Create a blackbody model called bb, check that it is
        recognized as a component, and display its parameters:

        >>> create_model_component("bbody", "bb")
        >>> list_model_components()
        >>> print(bb)
        >>> print(bb.ampl)

        """

        # If user mistakenly passes an actual model reference,
        # just return (i.e., create_model_component(const1d.c1)
        # is redundant, so just return)
        if isinstance(typename, Model) and name is None:
            return typename

        _check_str_type(typename, "typename")
        _check_str_type(name, "name")

        typename = typename.lower()
        cls = self._model_types.get(typename)
        if cls is None:
            raise ArgumentErr('badtype', typename)

        model = cls(name)
        self._model_components[name] = model
        return model

    def reset(self, model=None,
              id: Optional[IdType] = None
              ) -> None:
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
        id : int, str, or None, optional
           The data set to use. The default is to use all
           data sets with a source expression.

        See Also
        --------
        fit : Fit one or more data sets.
        guess : Set model parameters to values matching the data.
        paramprompt : Control how parameter values are set.
        set_par : Set the value, limits, or behavior of a model parameter.

        Examples
        --------

        The following examples assume that the source model has been
        set using:

        >>> set_source(powlaw1d.pl * xsphabs.gal)

        Fit the model and then reset the values of both components
        (``pl`` and ``gal``):

        >>> fit()
        >>> reset()

        Reset just the parameters of the ``pl`` model component:

        >>> reset(pl)

        Reset all the components of the source expression for data
        set 2.

        >>> reset(get_source(2))

        """
        if model is not None:
            model.reset()
            return

        if id is None:
            ids = list(self._sources.keys())
        else:
            ids = [id]

        for id in ids:
            self.get_source(id).reset()

    def delete_model_component(self, name: str) -> None:
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
        set_model_autoassign_func : Set the method used to create model component identifiers.

        Notes
        -----
        It is an error to try to delete a component that is part of a
        model expression - i.e. included as part of an expression in a
        `set_model` or `set_source` call. In such a situation, use the
        `delete_model` function to remove the source expression before
        calling `delete_model_component`.

        Examples
        --------

        If a model instance called ``pl`` has been created - e.g. by
        ``create_model_component('powlaw1d', 'pl')`` - then the
        following will remove it:

        >>> delete_model_component('pl')

        """
        _check_str_type(name, "name")
        mod = self._model_components.pop(name, None)
        if mod is None:
            raise IdentifierErr('nomodelcmpt', name)

        # If the component is part of a model expression we
        # warn the user but make no change.
        #
        for key in self.list_model_ids():

            # We can not guarantee that key exists in both
            # _models and _sources - in fact, it shouldn't,
            # so we need to be a bit careful.
            #
            has_model = key in self._models and mod in self._models[key]
            has_source = key in self._sources and mod in self._sources[key]

            if has_model or has_source:
                warning(f"the model component '{mod.name}' is found in model {key}" +
                        " and cannot be deleted")
                # restore the model component in use and return
                self._model_components[name] = mod
                return

        _remove_obj_from_main(name)

    # Back-compatibility
    # create_model = create_model_component

    #
    # Source models
    #

    def _eval_model_expression(self, expr, typestr='model'):
        try:
            return eval(expr, self._model_globals, self._model_components)
        except Exception as exc:
            raise ArgumentErr('badexpr', typestr, sys.exc_info()[1]) from exc

    def list_model_ids(self) -> list[IdType]:
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
        list_psf_ids : List of all the data sets with a PSF.
        set_model : Set the source model expression for a data set.

        """
        keys = list(self._models.keys())
        keys.extend(list(self._sources.keys()))
        keys = list(set(keys))
        return sorted(keys, key=str)

    # Return full model for fitting, plotting, etc.  Expects a corresponding
    # data set to be available.
    def _get_model_status(self,
                          id: Optional[IdType] = None):
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

    def _add_convolution_models(self,
                                id: Optional[IdType],
                                data, model, is_source):
        """Add in "hidden" components to the model expression.

        This handles PSF and table models (ensuring that the
        model is folded on the dataset and adding the response
        if necessary).
        """

        def fold(mdls):
            for mdl in mdls:
                if mdl in model:
                    mdl.fold(data)

        fold(self._tbl_models)

        if is_source:
            model = self._add_psf(id, data, model)
        else:
            fold(self._psf_models)

        return model

    def get_source(self,
                   id: Optional[IdType] = None
                   ) -> Model:
        """Return the source model expression for a data set.

        This returns the model expression created by `set_model` or
        `set_source`. It does not include any instrument response.

        Parameters
        ----------
        id : int, str, or None, optional
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
        sherpa.astro.ui.set_bkg_model : Set the background model expression for a data set.
        set_model : Set the source model expression for a data set.
        set_full_model : Define the convolved model expression for a data set.
        show_model : Display the source model expression for a data set.

        Examples
        --------

        Return the source expression for the default data set,
        display it, and then find the number of parameters in it:

        >>> src = get_source()
        >>> print(src)
        >>> len(src.pars)
        5

        Set the source expression for data set 'obs2' to be equal to
        the model of data set 'obs1' multiplied by a scalar value:

        >>> set_source('obs2', const1d.norm * get_source('obs1'))

        """

        idval = self._fix_id(id)
        mdl = self._models.get(idval, None)
        if mdl is not None:
            raise IdentifierErr("Convolved model\n"
                                f"'{mdl.name}'\n is set for "
                                f"dataset {idval}. You should use "
                                "get_model instead.")
        return self._get_item(idval, self._sources, 'source',
                              'has not been set, consider using '
                              'set_source() or set_model()')

    def get_model(self,
                  id: Optional[IdType] = None
                  ) -> Model:
        """Return the model expression for a data set.

        This returns the model expression for a data set, including
        any instrument response (e.g. PSF or ARF and RMF) whether
        created automatically or explicitly, with `set_full_model`.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set containing the source expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        instance
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
        sherpa.astro.ui.set_bkg_model : Set the background model expression for a data set.
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

        data = self.get_data(id)
        model, is_source = self._get_model_status(id)
        return self._add_convolution_models(id, data, model, is_source)

    def _runparamprompt(self, pars):

        if not self._paramprompt:
            return

        for par in pars:
            while True:
                inval = input(f"{par.fullname} parameter value [{par.val}] ")
                if inval == "":
                    break

                val = None
                minval = None
                maxval = None
                tokens = [t.strip() for t in inval.split(',')]
                ntokens = len(tokens)

                if ntokens > 3:
                    info("Error: Please provide a comma-separated "
                         "list of floats; e.g. val,min,max")
                    continue

                if tokens[0] != '':
                    try:
                        val = float(tokens[0])
                    except Exception as e:
                        info("Please provide a float value; %s", e)
                        continue

                if ntokens > 1 and tokens[1] != '':
                    try:
                        minval = float(tokens[1])
                    except Exception as e:
                        info("Please provide a float value; %s", e)
                        continue

                if ntokens > 2 and tokens[2] != '':
                    try:
                        maxval = float(tokens[2])
                    except Exception as e:
                        info("Please provide a float value; %s", e)
                        continue

                try:
                    self.set_par(par, val, minval, maxval)
                    break
                except Exception as e:
                    info(str(e))
                    continue

    # DOC-NOTE: also in sherpa.astro.utils
    # DOC-TODO: what examples/info should be talked about here?
    # (e.g. no PHA/ARF/RMF)
    # @loggable(with_id=True, with_keyword='model')
    def set_full_model(self, id, model=None):
        """Define the convolved model expression for a data set.

        The model expression created by `set_model` can be modified by
        "instrumental effects", such as a PSF set by `set_psf`.  The
        `set_full_model` function is for when this is not sufficient,
        and full control is needed. An example of when this would be
        if different PSF models should be applied to different source
        components.

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
        fit : Fit one or more data sets.
        set_psf : Add a PSF model to a data set.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        Some functions - such as `plot_source` - may not work for
        model expressions created by `set_full_model`.

        Examples
        --------

        Apply different PSFs to different components, as well as an
        unconvolved component:

        >>> load_psf("psf1", "psf1.dat")
        >>> load_psf("psf2", "psf2.dat")
        >>> smodel = psf1(gauss2d.src1) + psf2(beta2d.src2) + const2d.bgnd
        >>> set_full_model("src", smodel)

        """
        if model is None:
            id, model = model, id

        model = self._check_model(model)
        idval = self._fix_id(id)
        self._models[idval] = model
        self._runparamprompt(model.pars)

    # DOC-TODO: the .cache value appears to default to 5
    # @loggable(with_id=True, with_keyword="model")
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
        sherpa.astro.ui.set_bkg_model : Set the background model expression for a data set.
        set_full_model : Define the convolved model expression for a data set.
        show_model : Display the source model expression for a data set.
        set_par : Set the value, limits, or behavior of a model parameter.
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

        Model caching is available via the model ``cache`` attribute. A
        non-zero value for this attribute means that the results of
        evaluating the model will be cached if all the parameters are
        frozen, which may lead to a reduction in the time taken to
        evaluate a fit. A zero value turns off the caching.  The
        default setting for X-Spec and 1D analytic models is that
        ``cache`` is ``5``, but ``0`` for the 2D analytic models.

        The `integrate1d` model can be used to apply a numerical
        integration to an arbitrary model expression.

        Examples
        --------

        Create an instance of the `powlaw1d` model type, called ``pl``,
        and use it as the model for the default data set.

        >>> set_model(polynom1d.pl)

        Create a model for the default dataset which is the `xsphabs`
        model multiplied by the sum of an `xsapec` and `powlaw1d`
        models (the model components are identified by the labels
        ``gal``, ``clus``, and ``pl``).

        >>> set_model(xsphabs.gal * (xsapec.clus + powlaw1d.pl))

        Repeat the previous example, using a string to define the
        model expression:

        >>> set_model('xsphabs.gal * (xsapec.clus + powlaw1d.pl)')

        Use the same model component (``src``, a `gauss2d` model)
        for the two data sets ('src1' and 'src2').

        >>> set_model('src1',  gauss2d.src + const2d.bgnd1)
        >>> set_model('src2', src + const2d.bgnd2)

        Share an expression - in this case three gaussian lines -
        between three data sets. The normalization of this line
        complex is allowed to vary in data sets 2 and 3 (the ``norm2``
        and ``norm3`` components of the `const1d` model), and each data
        set has a separate `polynom1d` component (``bgnd1``, ``bgnd2``,
        and ``bgnd3``). The ``c1`` parameters of the `polynom1d` model
        components are thawed and then linked together (to reduce the
        number of free parameters):

        >>> lines = gauss1d.l1 + gauss1d.l2 + gauss1d.l3
        >>> set_model(1, lines + polynom1d.bgnd1)
        >>> set_model(2, lines * const1d.norm2 + polynom1d.bgnd2)
        >>> set_model(3, lines * const1d.norm3 + polynom1d.bgnd3)
        >>> thaw(bgnd1.c1, bgnd2.c1, bgnd3.c1)
        >>> link(bgnd2.c2, bgnd1.c1)
        >>> link(bgnd3.c3, bgnd1.c1)

        For this expression, the ``gal`` component is frozen, so it is
        not varied in the fit. The ``cache`` attribute is set to a
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

        model = self._check_model(model)
        idval = self._fix_id(id)

        # Fit(data, model) does the dimensionality validation. Note
        # that we assume that any convolution-style model (e.g. a PSF
        # or a PHA response) does not change the dimensionality check.
        #
        data = self._data.get(idval)
        if data is not None:
            Fit(data, model)

        self._sources[idval] = model
        self._runparamprompt(model.pars)

        # Delete any previous model set with set_full_model()
        mdl = self._models.pop(idval, None)
        if mdl is not None:
            warning("Clearing convolved model\n'%s'\nfor dataset %s",
                    mdl.name, idval)

    set_source = set_model

    def delete_model(self,
                     id: Optional[IdType] = None
                     ) -> None:
        """Delete the model expression for a data set.

        This removes the model expression, created by `set_model`,
        for a data set. It does not delete the components of the
        expression.

        Parameters
        ----------
        id : int, str, or None, optional
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
        idval = self._fix_id(id)
        self._models.pop(idval, None)
        self._sources.pop(idval, None)

    def _check_model(self, model: ModelType) -> Model:
        """Return a Model instance.

        Converts a string to the model instance if necessary.

        Parameters
        ----------
        model : `str` or a `sherpa.models.model.Model` object

        Returns
        -------
        model : `sherpa.models.model.Model` instance

        Raises
        ------
        ArgumentTypeError
            The argument is not a string or a model
        ArgumentErr
            The argument is a string but does not reference a model

        """
        if _is_str(model):
            mdl = self._eval_model_expression(model)
        else:
            mdl = model
        _check_type(mdl, Model, 'model',
                    'a model object or model expression string')
        return mdl

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

    def get_model_pars(self, model: ModelType) -> list[str]:
        """Return the names of the parameters of a model.

        .. versionadded:: 4.16.1
           Linked components are now included without having to
           include them in the model expression.

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

        >>> mdl = gauss2d.src + const2d.bgnd
        >>> get_model_pars(mdl)
        ['fwhm', 'xpos', 'ypos', 'ellip', 'theta', 'ampl', 'c0']

        The return value is unchanged if the two models are linked
        together, such as setting the background amplitude to 1% of
        the source amplitude:

        >>> bgnd.c0 = 0.01 * src.ampl
        >>> get_model_pars(mdl)
        ['fwhm', 'xpos', 'ypos', 'ellip', 'theta', 'ampl', 'c0']

        If a model expression contains linked parameters that are not
        part of the model expression then they will also be included
        (in this case both the const2d and scale1d parameters are
        named 'c0', hence the duplication):

        >>> scale1d.sep
        <Scale1D model instance 'scale1d.sep'>
        >>> src.ypos = src.xpos + sep.c0
        >>> get_model_pars(mdl)
        ['fwhm', 'xpos', 'ypos', 'ellip', 'theta', 'ampl', 'c0', 'c0']

        """
        mdl = self._check_model(model)
        # TODO: we could add the component name
        names = [p.name for p in mdl.pars]
        names.extend(p.name for p in mdl.lpars)
        return names

    def get_num_par(self, id: Optional[IdType] = None) -> int:
        """Return the number of parameters in a model expression.

        The `get_num_par` function returns the number of parameters,
        both frozen and thawed, in the model assigned to a data set.

        .. versionadded:: 4.16.1
           Linked components are now included without having to
           include them in the model expression.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set containing the model expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        npar : int
           The number of parameters in the model expression. This
           sums up all the parameters of the components in the
           expression, and includes both frozen and thawed components.

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

        Examples
        --------

        Return the total number of parameters for the default data set:

        >>> print(get_num_par())

        Find the number of parameters for the model associated with
        the data set called "jet":

        >>> njet = get_num_par('jet')

        """
        mdl = self.get_source(id)
        return len(mdl.pars) + len(mdl.lpars)

    def get_num_par_thawed(self, id: Optional[IdType] = None) -> int:
        """Return the number of thawed parameters in a model expression.

        The `get_num_par_thawed` function returns the number of
        thawed parameters in the model assigned to a data set.

        .. versionadded:: 4.16.1
           Linked components are now included without having to
           include them in the model expression.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set containing the model expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        npar : int
           The number of parameters in the model expression. This
           sums up all the thawed parameters of the components in the
           expression.

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

        Examples
        --------

        Return the number of thawed parameters for the default data set:

        >>> print(get_num_par_thawed())

        Find the number of thawed parameters for the model associated
        with the data set called "jet":

        >>> njet = get_num_par_thawed('jet')

        """
        mdl = self.get_source(id)
        return len(mdl.thawedpars)

    def get_num_par_frozen(self, id: Optional[IdType] = None) -> int:
        """Return the number of frozen parameters in a model expression.

        The `get_num_par_frozen` function returns the number of
        frozen parameters in the model assigned to a data set.

        .. versionadded:: 4.16.1
           Linked components are now included without having to
           include them in the model expression.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set containing the model expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.

        Returns
        -------
        npar : int
           The number of parameters in the model expression. This
           sums up all the frozen parameters of the components in the
           expression.

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

        Examples
        --------

        Return the number of frozen parameters for the default data set:

        >>> print(get_num_par_frozen())

        Find the number of frozen parameters for the model associated
        with the data set called "jet":

        >>> njet = get_num_par_frozen('jet')

        """
        return self.get_num_par(id) - self.get_num_par_thawed(id)

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
            x = data.get_x()
            y = data.get_y()

        except TypeError:
            # we have to check for the case of a *single* column in the file
            # extract the single array from the read and bypass the dataset
            y = sherpa.io.get_ascii_data(filename, ncols=1, colkeys=colkeys,
                                         sep=sep, dstype=dstype,
                                         comment=comment)[1].pop()

        return (x, y)

    # DOC-TODO: I am not sure I have the data format correct.
    # DOC-TODO: description of template interpolation needs a lot of work.
    def load_template_model(self, modelname, templatefile, dstype=sherpa.data.Data1D,
                            sep=' ', comment='#', method=sherpa.utils.linear_interp, template_interpolator_name='default'):
        """Load a set of templates and use it as a model component.

        A template model can be considered to be an extension
        of the table model supported by `load_table_model`. In
        the template case, a set of models (the "templates")
        are read in and then compared to the data, with the
        best-fit being used to return a set of parameters.

        Parameters
        ----------
        modelname : str
           The identifier for this table model.
        templatefile : str
           The name of the file to read in. This file lists the
           template data files.
        dstype : data class to use, optional
           What type of data is to be used. This is currently unused.
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.
        method : func
           The interpolation method to use to map the input data onto
           the coordinate grid of the data set. Linear,
           nearest-neighbor, and polynomial schemes are provided in
           the sherpa.utils module.
        template_interpolator_name : str
           The method used to interpolate within the set of templates.
           The default is ``default``. A value of ``None`` turns off the
           interpolation; in this case the grid-search optimiser
           must be used to fit the data.

        See Also
        --------
        load_conv : Load a 1D convolution model.
        load_psf : Create a PSF model
        load_table_model : Load tabular data and use it as a model component.
        load_template_interpolator : Set the template interpolation scheme.
        set_model : Set the source model expression for a data set.
        set_full_model : Define the convolved model expression for a data set.

        Notes
        -----
        Examples of interpolation schemes provided by `sherpa.utils`
        are: `linear_interp`, `nearest_interp`, and `neville`.

        The template index file is the argument to
        `load_template_model`, and is used to list the data files. It
        is an ASCII file with one line per template, and each line
        containing the model parameters (numeric values), followed by
        the MODELFLAG column and then the file name for the data file
        (its name must begin with FILE). The MODELFLAG column is used
        to indicate whether a file should be used or not; a value of
        1 means that the file should be used, and a value of 0
        that the line should be ignored. The parameter names are set
        by the column names.

        The data file - the last column of the template index file -
        is read in and the first two columns used to set up the x and
        y values (`Data1D`) or xlo, xhi, and y values (`Data1DInt`).
        These files must be in ASCII format.

        The ``method`` parameter determines how the template data values
        are interpolated onto the source data grid.

        The ``template_interpolator_name`` parameter determines how the
        dependent axis (Y) values are interpolated when the parameter
        values are varied. This interpolation can be turned off by
        using a value of ``None``, in which case the grid-search
        optimiser *must* be used. See `load_template_interpolator` for
        how to create a valid interpolator. The "default" interpolator
        uses `sherpa.models.KNNInterpolator` with k=2 and order=2.

        Examples
        --------

        Load in the templates from the file "index.tmpl" as the model
        component "kerr", and set it as the source model for the
        default data set. The optimisation method is switched to
        use a grid search for the parameters of this model.

        >>> load_template_model("kerr", "index.tmpl")
        >>> set_source(kerr)
        >>> set_method('gridsearch')
        >>> set_method_opt('sequence', kerr.parvals)
        >>> fit()

        Fit a constant plus the templates, using the neville scheme
        for integrating the template onto the data grid. The
        Monte-Carlo based optimiser is used.

        >>> load_template_model('tbl', 'table.lis',
        ...                     sherpa.utils.neville)
        >>> set_source(tbl + const1d.bgnd)
        >>> set_method('moncar')

        """

        # Note: dstype is not used
        templatemodel = read_template_model(modelname, templatefile,
                                            sep=sep, comment=comment,
                                            method=method,
                                            template_interpolator_name=template_interpolator_name)

        self._tbl_models.append(templatemodel)
        self._add_model_component(templatemodel)

    # DOC-TODO: description of template interpolation needs a lot of work.
    # @loggable()
    def load_template_interpolator(self, name, interpolator_class, **kwargs):
        """Set the template interpolation scheme.

        Parameters
        ----------
        name : str
        interpolator_class
           An interpolator class.
        **kwargs
           The arguments for the interpolator.

        See Also
        --------
        load_template_model : Load a set of templates and use it as a model component.

        Examples
        --------

        Create an interpolator name that can be used as the
        ``template_interpolator_name`` argument to
        `load_template_model`.

        >>> from sherpa.models import KNNInterpolator
        >>> load_template_interpolator('myint', KNNInterpolator, k=4, order=3)

        """
        add_interpolator(name, interpolator_class, **kwargs)

    def load_table_model(self, modelname, filename, ncols=2, colkeys=None,
                         dstype=sherpa.data.Data1D, sep=' ', comment='#',
                         method=sherpa.utils.linear_interp):
        """Load ASCII tabular data and use it as a model component.

        A table model is defined on a grid of points which is
        interpolated onto the independent axis of the data set. The
        model has a single parameter, ``ampl``, which is used to scale
        the data, and it can be fixed or allowed to vary during a fit.

        Parameters
        ----------
        modelname : str
           The identifier for this table model.
        filename : str
           The name of the ASCII file to read in.
        ncols : int, optional
           The number of columns to read in (the first ``ncols`` columns
           in the file). It should be 1 or 2.
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           ``None``, which uses the first ``ncols`` columns in the file.
           The default column names are col followed by the column number,
           so ``col1`` for the first column.
        dstype : data class to use, optional
           What type of data is to be used. Supported values include
           `Data1D` (the default) and `Data1DInt`.
        sep : str, optional
           The separator character for columns. The default is ``' '``.
        comment : str, optional
           Lines starting with this string are ignored. The default
           is ``'#'``.
        method : func
           The interpolation method to use to map the input data onto
           the coordinate grid of the data set. Linear,
           nearest-neighbor, and polynomial schemes are provided in
           the sherpa.utils module.

        See Also
        --------
        load_conv : Load a 1D convolution model.
        load_psf : Create a PSF model
        load_template_model : Load a set of templates and use it as a model component.
        set_model : Set the source model expression for a data set.
        set_full_model : Define the convolved model expression for a data set.

        Notes
        -----
        Examples of interpolation schemes provided by `sherpa.utils`
        are: `linear_interp`, `nearest_interp`, `neville`, and
        `neville2d`.

        See `unpack_data` for a description of the supported
        file format.

        When reading in two columns, the data will be re-ordered
        so that the first column read in (the independent axis)
        is numerically increasing.

        If ``ncols=1``, only the model values (dependent axis) are
        read in. In this case, the data set to which the model
        is applied - via ``set_source`` - must have the same number
        of data points as the model.

        When used with an integrated data set (for example,
        `Data1DInt`), then the first column of the table - the
        independent axis - should be the left-edge of the bin,
        and the second column is the integrated value for that
        bin.

        Examples
        --------

        Load in the data from filt.dat and use it to multiply
        the source model (a power law and a gaussian). Allow
        the amplitude for the table model to vary between 1
        and 1e6, starting at 1e3.

        >>> load_table_model('filt', 'filt.dat')
        >>> set_source(filt * (powlaw1d.pl + gauss1d.gline))
        >>> set_par(filt.ampl, 1e3, min=1, max=1e6)

        """

        x, y = self._read_user_model(filename, ncols, colkeys,
                                     dstype, sep, comment)

        tablemodel = TableModel(modelname)
        tablemodel.method = method
        tablemodel.filename = filename
        tablemodel.load(x, y)
        self._tbl_models.append(tablemodel)
        self._add_model_component(tablemodel)

    # also in sherpa.astro.utils
    # DOC-TODO: how is the _y value used if set
    # @loggable()
    def load_user_model(self, func, modelname, filename=None, ncols=2,
                        colkeys=None, dstype=sherpa.data.Data1D,
                        sep=' ', comment='#'):
        """Create a user-defined model.

        Assign a name to a function; this name can then be used as any
        other name of a model component, either in a source expression
        - such as with `set_model` - or to change a parameter
        value. The `add_user_pars` function should be called after
        `load_user_model` to set up the parameter names and
        defaults.

        Parameters
        ----------
        func : func
           The function that evaluates the model.
        modelname : str
           The name to use to refer to the model component.
        filename : str, optional
           Set this to include data from this file in the model. The
           file should contain two columns, and the second column is
           stored in the ``_y`` attribute of the model.
        ncols : int, optional
           The number of columns to read in (the first ``ncols`` columns
           in the file).
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           ``None``.
        dstype : data class to use, optional
           What type of data is to be used. Supported values include
           `Data1D` (the default), `Data1DInt`, `Data2D`, and
           `Data2DInt`.
        sep : str, optional
           The separator character. The default is ``' '``.
        comment : str, optional
           The comment character. The default is ``'#'``.

        See Also
        --------
        add_model : Create a user-defined model class.
        add_user_pars : Add parameter information to a user model.
        load_table_model : Load tabular data and use it as a model component.
        load_template_model : Load a set of templates and use it as a model component.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The `load_user_model` function is designed to make it easy to
        add a model, but the interface is not the same as the existing
        models (such as having to call both `load_user_model` and
        `add_user_pars` for each new instance).  The `add_model`
        function is used to add a model as a Python class, which is
        more work to set up, but then acts the same way as the
        existing models.

        The function used for the model depends on the dimensions of
        the data. For a 1D model, the signature is::

           def func1d(pars, x, xhi=None):

        where, if xhi is not None, then the dataset is binned and the
        x argument is the low edge of each bin. The pars argument is
        the parameter array - the names, defaults, and limits can be
        set with `add_user_pars` - and should not be changed.  The
        return value is an array the same size as x.

        For 2D models, the signature is::

           def func2d(pars, x0, x1, x0hi=None, x1hi=None):

        There is no way using this interface to indicate that the
        model is for 1D or 2D data.

        The format for the input file, and how to control what columns
        to read, are described in the help for the `unpack_data` function.

        Examples
        --------

        Create a two-parameter model of the form "y = mx + c",
        where the intercept is the first parameter and the slope the
        second, set the parameter names and default values, then
        use it in a source expression::

            >>> def func1d(pars, x, xhi=None):
            ...     if xhi is not None:
            ...         x = (x + xhi) / 2
            ...     return x * pars[1] + pars[0]
            ...
            >>> load_user_model(func1d, "myfunc")
            >>> add_user_pars(myfunc, ["c", "m"], [0, 1])
            >>> set_source(myfunc + gauss1d.gline)

        """
        usermodel = sherpa.models.UserModel(modelname)
        usermodel.calc = func
        usermodel._file = filename
        if filename is not None:
            _, usermodel._y = self._read_user_model(filename, ncols, colkeys,
                                                    dstype, sep, comment)
        self._add_model_component(usermodel)

    # @loggable()
    def add_user_pars(self, modelname, parnames,
                      parvals=None, parmins=None, parmaxs=None,
                      parunits=None, parfrozen=None
                      ) -> None:
        """Add parameter information to a user model.

        Parameters
        ----------
        modelname : str
           The name of the user model (created by `load_user_model`).
        parnames : array of str
           The names of the parameters. The order of all the parameter
           arrays must match that expected by the model function (the
           first argument to `load_user_model`).
        parvals : array of number, optional
           The default values of the parameters. If not given each
           parameter is set to 0.
        parmins : array of number, optional
           The minimum values of the parameters (hard limit). The
           default value is -3.40282e+38.
        parmaxs : array of number, optional
           The maximum values of the parameters (hard limit). The
           default value is 3.40282e+38.
        parunits : array of str, optional
           The units of the parameters. This is only used in screen
           output (i.e. is informational in nature).
        parfrozen : array of bool, optional
           Should each parameter be frozen. The default is that
           all parameters are thawed.

        See Also
        --------
        add_model : Create a user-defined model class.
        load_user_model : Create a user-defined model.
        set_par : Set the value, limits, or behavior of a model parameter.

        Notes
        -----

        The parameters must be specified in the order that the function
        expects. That is, if the function has two parameters,
        pars[0]='slope' and pars[1]='y_intercept', then the call to
        add_user_pars must use the order ["slope", "y_intercept"].

        Examples
        --------

        Create a user model for the function `profile` called
        "myprof", which has two parameters called "core" and "ampl",
        both of which will start with a value of 0.

        >>> load_user_model(profile, "myprof")
        >>> add_user_pars("myprof", ["core", "ampl"])

        Set the starting values, minimum values, and whether or not
        the parameter is frozen by default, for the "prof" model:

        >>> pnames = ["core", "ampl", "intflag"]
        >>> pvals = [10, 200, 1]
        >>> pmins = [0.01, 0, 0]
        >>> pfreeze = [False, False, True]
        >>> add_user_pars("prof", pnames, pvals,
        ...               parmins=pmins, parfrozen=pfreeze)

        """

        _check_str_type(modelname, "model name")

        usermodel = self._get_model_component(modelname)
        if (usermodel is None or
            not isinstance(usermodel, sherpa.models.UserModel)):
            raise ArgumentTypeErr('badarg', modelname, "a user model")

        pars = []
        vals = None
        if parvals is not None:
            vals = list(parvals)
        mins = None
        if parmins is not None:
            mins = list(parmins)
        maxs = None
        if parmaxs is not None:
            maxs = list(parmaxs)
        units = None
        if parunits is not None:
            units = list(parunits)
        frozen = None
        if parfrozen is not None:
            frozen = list(parfrozen)

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

    # DOC-TODO: Improve priors documentation
    def load_user_stat(self, statname, calc_stat_func, calc_err_func=None,
                       priors={}
                       ) -> None:
        """Create a user-defined statistic.

        The choice of statistics - that is, the numeric value that is
        minimised during the fit - can be extended by providing a
        function to calculate a numeric value given the data. The
        statistic is given a name and then can be used just like any
        of the pre-defined statistics.

        Parameters
        ----------
        statname : str
           The name to use for the new statistic when calling
           `set_stat`.
        calc_stat_func : func
           The function that calculates the statistic.
        calc_err_func : func, optional
           How to calculate the statistical error on a data point.
        priors : dict
           A dictionary of hyper-parameters for the priors.

        See Also
        --------
        calc_stat : Calculate the fit statistic for a data set.
        set_stat : Set the statistical method.

        Notes
        -----
        The ``calc_stat_func`` should have the following signature::

           def func(data, model, staterror=None, syserrr=None, weight=None)

        where data is the array of dependent values, model the array
        of the predicted values, staterror and syserror are arrays of
        statistical and systematic errors respectively (if valid), and
        weight an array of weights. The return value is the pair
        (stat_value, stat_per_bin), where stat_value is a scalar and
        stat_per_bin is an array the same length as data.

        The ``calc_err_func`` should have the following signature::

           def func(data)

        and returns an array the same length as data.

        Examples
        --------

        Define a chi-square statistic with the label "qstat":

            >>> def qstat(d, m, staterr=None, syserr=None, w=None):
            ...     if staterr is None:
            ...         staterr = 1
            ...     c = ((d-m) / staterr)
            ...     return ((c*c).sum(), c)
            ...
            >>> load_user_stat("qstat", qstat)
            >>> set_stat("qstat")

        """
        userstat = UserStat(calc_stat_func, calc_err_func, statname)
        if priors:
            assert False
            # TODO: clean this up
            pars = [(key, priors.pop(key)) for key in priors.keys()
                    if isinstance(priors[key], sherpa.models.Parameter)]
            kwargs = dict(pars)
            userstat = sherpa.logposterior.Prior(calc_stat_func, priors, kwargs)

        _assign_obj_to_main(statname, userstat)

    # Back-compatibility
    # set_source = set_model

    #
    # Conv
    #

    # DOC-NOTE: why isn't the "flux" of the convolved model ~
    # that of the unconvolved model?
    # DOC-NOTE: better description of conv vs psf
    # @loggable()
    def load_conv(self, modelname, filename_or_model, *args, **kwargs):
        """Load a 1D convolution model.

        The convolution model can be defined either by a data set,
        read from a file, or an analytic model, using a Sherpa model
        instance. A source model can be convolved with this model
        by including ``modelname`` in the `set_model` call, using the
        form::

           modelname(modelexpr)

        Parameters
        ----------
        modelname : str
           The identifier for this PSF model.
        filename_or_model : str or model instance
           This can be the name of an ASCII file or a Sherpa model component.
        args
           Arguments for `unpack_data` if `filename_or_model` is a file.
        kwargs
           Keyword arguments for `unpack_data` if `filename_or_model` is
           a file.

        See Also
        --------
        delete_psf : Delete the PSF model for a data set.
        load_psf : Create a PSF model.
        load_table_model : Load tabular data and use it as a model component.
        set_full_model : Define the convolved model expression for a data set.
        set_model : Set the source model expression for a data set.
        set_psf : Add a PSF model to a data set.

        Examples
        --------

        Create a 1D data set, assign a box model - which is flat
        between the xlow and xhi values and zero elsewhere - and then
        display the model values. Then add in a convolution component
        by a gaussian and overplot the resulting source model with two
        different widths.

        >>> dataspace1d(-10, 10, 0.5, id='tst', dstype=Data1D)
        >>> set_source('tst', box1d.bmdl)
        >>> bmdl.xlow = -2
        >>> bmdl.xhi = 3
        >>> plot_source('tst')
        >>> load_conv('conv', normgauss1d.gconv)
        >>> gconv.fwhm = 2
        >>> set_source('tst', conv(bmdl))
        >>> plot_source('tst', overplot=True)
        >>> gconv.fwhm = 5
        >>> plot_source('tst', overplot=True)

        Create a convolution component called "cmodel" which uses the
        data in the file "conv.dat", which should have two columns
        (the X and Y values).

        >>> load_conv('cmodel', 'conv.dat')

        """
        kernel = filename_or_model
        if _is_str(filename_or_model):
            try:
                kernel = self._eval_model_expression(filename_or_model)
            except Exception:
                kernel = self.unpack_data(filename_or_model,
                                          *args, **kwargs)

        conv = sherpa.instrument.ConvolutionKernel(kernel, modelname)
        self._add_model_component(conv)

    #
    # PSF1
    #

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
           This can be the name of an ASCII file or a Sherpa model component.
        args
           Arguments for `unpack_data` if `filename_or_model` is a file.
        kwargs
           Keyword arguments for `unpack_data` if `filename_or_model` is
           a file.

        See Also
        --------
        delete_psf : Delete the PSF model for a data set.
        load_conv : Load a 1D convolution model.
        load_table_model : Load tabular data and use it as a model component.
        set_full_model : Define the convolved model expression for a data set.
        set_model : Set the source model expression for a data set.
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

        Create a PSF model from the data in the ASCII file
        'line_profile.dat' and apply it to the data set called
        'bgnd':

        >>> load_psf('pmodel', 'line_profile.dat')
        >>> set_psf('bgnd', 'pmodel')

        """
        kernel = filename_or_model
        if _is_str(filename_or_model):
            try:
                kernel = self._eval_model_expression(filename_or_model)
            except Exception:
                kernel = self.unpack_data(filename_or_model,
                                          *args, **kwargs)

        psf = sherpa.instrument.PSFModel(modelname, kernel)
        self._add_model_component(psf)
        self._psf_models.append(psf)

    # DOC-TODO: am I correct about the multiple use warning?
    # @loggable(with_id=True, with_keyword='psf')
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
        psf : str or `sherpa.instrument.PSFModel` instance
           The PSF model created by `load_psf`.

        See Also
        --------
        delete_psf : Delete the PSF model for a data set.
        get_psf : Return the PSF model defined for a data set.
        image_psf : Display the 2D PSF model for a data set in the image viewer.
        load_psf : Create a PSF model.
        plot_psf : Plot the 1D PSF model applied to a data set.
        set_full_model : Define the convolved model expression for a data set.
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
           Set to ``1`` to use a symmetric array. The default is ``0`` to
           reduce edge effects.

        norm
           Should the kernel be normalized so that it sums to 1?  This
           summation is done over the full data set (not the subset
           defined by the ``size`` parameter). The default is ``1`` (yes).

        Examples
        --------

        Use the data in the ASCII file 'line_profile.dat' as the PSF for
        the default data set:

        >>> load_psf('psf1', 'line_profile.dat')
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

        if _is_str(psf):
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
            if np.isscalar(psf_center):
                psf_center = [psf_center]
            try:
                center = psf.kernel.get_center()
                if (np.asarray(center) == 0.0).all():
                    psf.kernel.set_center(*psf_center, values=True)
            except NotImplementedError:
                pass

    def get_psf(self, id: Optional[IdType] = None):
        """Return the PSF model defined for a data set.

        Return the parameter settings for the PSF model assigned to
        the data set.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        psf : a `sherpa.instrument.PSFModel` instance

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If no PSF model has been set for the data set.

        See Also
        --------
        delete_psf : Delete the PSF model for a data set.
        image_psf : Display the 2D PSF model for a data set in the image viewer.
        list_psf_ids : List of all the data sets with a PSF.
        load_psf : Create a PSF model.
        plot_psf : Plot the 1D PSF model applied to a data set.
        set_psf : Add a PSF model to a data set.

        Examples
        --------

        Change the size and center of the PSF for the default data set:

        >>> psf = get_psf()
        >>> psf.size = (21, 21)
        >>> psf.center = (10, 10)

        """
        return self._get_item(id, self._psf, 'psf model', 'has not been set')

    def delete_psf(self, id: Optional[IdType] = None) -> None:
        """Delete the PSF model for a data set.

        Remove the PSF convolution applied to a source model.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the
           default identifier is used, as returned by `get_default_id`.

        See Also
        --------
        list_psf_ids : List of all the data sets with a PSF.
        load_psf : Create a PSF model.
        set_psf : Add a PSF model to a data set.
        get_psf : Return the PSF model defined for a data set.

        Examples
        --------

        >>> delete_psf()

        >>> delete_psf('core')

        """
        idval = self._fix_id(id)
        self._psf.pop(idval, None)

    def list_psf_ids(self) -> list[IdType]:
        """List of all the data sets with a PSF.

        .. versionadded:: 4.12.2

        Returns
        -------
        ids : list of int or str
           The identifiers for all the data sets which have a PSF
           model set by `set_psf`.

        See Also
        --------
        list_data_ids : List the identifiers for the loaded data sets.
        list_model_ids : List of all the data sets with a source expression.
        set_psf : Add a PSF model to a data set.

        """
        keys = list(self._psf.keys())
        return sorted(keys, key=str)

    def _add_psf(self, id: Optional[IdType], data, model):
        idval = self._fix_id(id)
        psf = self._psf.get(idval, None)

        if psf is not None:
            model = psf(model)
            psf.fold(data)
        return model

    #
    # Parameters
    #

    def _check_par(self, par, argname='par'):
        if _is_str(par):
            par = self._eval_model_expression(par, 'parameter')
        _check_type(par, sherpa.models.Parameter, argname,
                    'a parameter object or parameter expression string')
        return par

    # DOC-NOTE: I have not documented that par can be an actual parameter
    # since in that case get_par(x) === x, so seems somewhat pointless!
    def get_par(self, par):
        """Return a parameter of a model component.

        Parameters
        ----------
        par : str
           The name of the parameter, using the format
           "componentname.parametername".

        Returns
        -------
        par : a `sherpa.models.parameter.Parameter` instance
           The parameter values - e.g. current value, limits, and
           whether it is frozen - can be changed using this
           object.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the `par` argument is invalid: the model component does
           not exist or the given model has no parameter with that
           name.

        See Also
        --------
        set_par : Set the value, limits, or behavior of a model parameter.

        Examples
        --------

        Return the "c0" parameter of the "bgnd" model component
        and change it to be frozen:

        >>> p = get_par('bgnd.c0')
        >>> p.frozen = True

        """
        return self._check_par(par)

    # DOC-NOTE: I have not documented that par can be an actual parameter
    # since you can just change the values directly then (although
    # may have to car about order of operations)
    def set_par(self, par, val=None, min=None, max=None, frozen=None):
        """Set the value, limits, or behavior of a model parameter.

        Parameters
        ----------
        par : str
           The name of the parameter, using the format
           "componentname.parametername".
        val : number, optional
           Change the current value of the parameter.
        min : number, optional
           Change the minimum value of the parameter (the soft limit).
        max : number, optional
           Change the maximum value of the parameter (the soft limit).
        frozen : bool, optional
           Freeze (``True``) or thaw (``False``) the parameter.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the ``par`` argument is invalid: the model component does
           not exist or the given model has no parameter with that
           name.

        See Also
        --------
        freeze : Fix model parameters so they are not changed by a fit.
        get_par : Return a parameter of a model component.
        link : Link a parameter value to an associated value.
        thaw : Allow model parameters to be varied during a fit.
        unlink : Unlink a parameter value.

        Notes
        -----
        The parameter object can be used to change these values
        directly, by setting the attribute with the same
        name as the argument - so that::

           set_par('emis.flag', val=2, frozen=True)

        is the same as::

           emis.flag.val = 2
           emis.flag.frozen = True

        Examples
        --------

        Change the parameter value to 23.

        >>> set_par('bgnd.c0', 23)

        Restrict the line.ampl parameter to be between 1e-4
        and 10 and to have a value of 0.1.

        >>> set_par('line.ampl', 0.1, min=1e-4, max=10)

        """
        self._check_par(par).set(val, min, max, frozen)

    def freeze(self, *args):
        """Fix model parameters so they are not changed by a fit.

        The arguments can be parameters or models, in which case all
        parameters of the model are frozen. If no arguments are
        given then nothing is changed.

        Parameters
        ----------
        args : sequence of str or Parameter or Model
            The parameters or models to freeze.

        See Also
        --------
        fit : Fit one or more data sets.
        link : Link a parameter value to an associated value.
        set_par : Set the value, limits, or behavior of a model parameter.
        thaw : Allow model parameters to be varied during a fit.
        unlink : Unlink a parameter value.

        Notes
        -----

        The `thaw` function can be used to reverse this setting,
        so that parameters can be varied in a fit.

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
        for par in list(args):
            if _is_str(par):
                par = self._eval_model_expression(par, 'parameter or model')

            try:
                par.freeze()
            except AttributeError as exc:
                raise ArgumentTypeErr('badarg', 'par',
                                      'a parameter or model object or expression string') from exc

    def thaw(self, *args):
        """Allow model parameters to be varied during a fit.

        The arguments can be parameters or models, in which case all
        parameters of the model are thawed. If no arguments are
        given then nothing is changed.

        Parameters
        ----------
        args : sequence of str or Parameter or Model
            The parameters or models to thaw.

        See Also
        --------
        fit : Fit one or more data sets.
        freeze : Fix model parameters so they are not changed by a fit.
        link : Link a parameter value to an associated value.
        set_par : Set the value, limits, or behavior of a model parameter.
        unlink : Unlink a parameter value.

        Notes
        -----

        The `freeze` function can be used to reverse this setting,
        so that parameters are "frozen" and so remain constant during
        a fit.

        Certain parameters may be marked as "always frozen", in which
        case using the parameter in a call to `thaw` will raise an
        error. If the model is sent to `thaw` then the "always frozen"
        parameter will be skipped.

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
        for par in list(args):
            if _is_str(par):
                par = self._eval_model_expression(par, 'parameter or model')

            try:
                par.thaw()
            except AttributeError as exc:
                raise ArgumentTypeErr('badarg', 'par',
                                      'a parameter or model object or expression string') from exc

    def link(self, par, val):
        """Link a parameter to a value.

        A parameter can be linked to another parameter value, or
        function of that value, rather than be an independent value.
        As the linked-to values change, the parameter value will
        change.

        .. versionchanged:: 4.16.1
           Source models no longer have to contain the linked parameter.

        Parameters
        ----------
        par : str or Parameter
           The parameter to link.
        val
           The value - which can be a numeric value or a function
           of other model parameters, to set `par` to.

        See Also
        --------
        freeze : Fix model parameters so they are not changed by a fit.
        set_par : Set the value, limits, or behavior of a model parameter.
        thaw : Allow model parameters to be varied during a fit.
        unlink : Unlink a parameter value.

        Notes
        -----
        The ``link`` attribute of the parameter is set to match the
        mathematical expression used for `val`.

        In the following, the fit will vary the ``pos`` parameter
        even though the ``src2`` component is not part of the source
        expression (this behavior changed in the 4.16.1 release):

        >>> set_source(1, gauss1d.src1)
        >>> gauss1d.src2
        >>> link(src1.pos, src2.pos)
        >>> fit(1)

        The ``lpars`` attribute of a source model will include these
        "extra" parameters:

        >>> get_source(1).lpars
        (<Parameter 'pos' of model 'src2'>,)

        Examples
        --------

        The ``fwhm`` parameter of the ``g2`` model is set to be the
        same as the ``fwhm`` parameter of the ``g1`` model.

        >>> link(g2.fwhm, g1.fwhm)

        Fix the ``pos`` parameter of ``g2`` to be 2.3 more than the
        ``pos`` parameter of the ``g1`` model.

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
        if _is_str(val):
            val = self._eval_model_expression(val, 'parameter link')
        par.link = val

    def unlink(self, par):
        """Unlink a parameter value.

        Remove any parameter link - created by `link` - for the
        parameter. The parameter value is reset to the value it
        had before `link` was called.

        Parameters
        ----------
        par : str or Parameter
           The parameter to unlink. If the parameter is not linked
           then nothing happens.

        See Also
        --------
        freeze : Fix model parameters so they are not changed by a fit.
        link : Link a parameter to a value.
        set_par : Set the value, limits, or behavior of a model parameter.
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

    def _get_fit_ids(self,
                     id: Optional[IdType],
                     otherids: Optional[Sequence[IdType]] = None
                     ) -> list[IdType]:
        """Return the identifiers that will be used for a fit.

        This routine does not ensure that the dataset actually exists.

        Parameters
        ----------
        id: int or str or None
            If None then this will select all data sets.
        otherids: sequence of int or str or None, or None
            When id is not None, the other identifiers to use.

        Returns
        -------
        idvals : list of int or str
            The identifiers that may contain data. It will not
            be empty.

        """

        if id is None:
            ids = self.list_data_ids()
            if len(ids) == 0:
                raise IdentifierErr("getitem", "data set", self._default_id, "has not been set")

            return ids

        ids = [id]
        if otherids is None:
            return ids

        for idval in otherids:
            idval = self._fix_id(idval)
            if idval not in ids:
                ids.append(idval)

        return ids

    def _get_fit_obj(self,
                     store: Sequence[FitStore],
                     estmethod, numcores=1
                     ) -> tuple[tuple[IdType, ...], Fit]:
        """Create the fit object given the data and models.

        Parameters
        ----------
        store : list of FitStore
            The per-dataset data and model information.
        estmethod : `sherpa.estmethods.EstMethod` or None
            Passed to the Fit object.
        numcores : int, optional
            The number of CPU cores to use (this is used when
            evaluating the models for multiple data sets).

        Returns
        -------
        ids, fit : tuple, `sherpa.fit.Fit` instance
            The datasets used and the fit object.

        Notes
        -----

        The data needed is passed around as a list of tuples as it
        makes it easy to pass in extra information in derived classes.
        It could be updated to pass around an object that exposes
        given fields (e.g. idval, dataset, and model) but at the
        moment it is not needed.

        """

        if not self._current_method.name == 'gridsearch':
            for s in store:
                if s.model.is_discrete:
                    raise ModelErr(
                        "You are trying to fit a model which has a discrete template model component with a continuous optimization method. Since CIAO4.6 this is not possible anymore. Please use gridsearch as the optimization method and make sure that the 'sequence' option is correctly set, or enable interpolation for the templates you are loading (which is the default behavior).")

        # Data and DataSimulFit do not have a common base class.
        #
        d: Union[Data, DataSimulFit]
        if len(store) == 1:
            d = store[0].data
            m = store[0].model
        else:
            datasets = [s.data for s in store]
            d = DataSimulFit('simulfit data', datasets, numcores)

            models = [s.model for s in store]
            m = SimulFitModel('simulfit model', models)

        # Ensure the id value is not repeated, but keep the order (so can not
        # just convert to a set and back again).
        #
        idvals = []
        for s in store:
            idval = s.idval
            if idval not in idvals:
                idvals.append(idval)

        return tuple(idvals), Fit(d, m, self._current_stat, self._current_method,
                                  estmethod, self._current_itermethod)

    def _prepare_fit(self,
                     id: Optional[IdType],
                     otherids: Sequence[IdType] = ()
                     ) -> list[FitStore]:
        """Ensure we have all the requested ids, datasets, and models.

        This checks whether the dataset is loaded and has an
        associated source or model. If datasets are explicitly listed
        then they must contain both data and a model, but when id is
        None then those data sets which have no model are skipped.

        Parameters
        ----------
        id: int or str or None
            If `None` then this fits all data simultaneously.
        otherids: sequence of int or str or None, or None
            When id is not None, the other identifiers to use.

        Returns
        -------
        store : list of FitStore

        Raises
        ------
        IdentifierErr
            If there are no datasets with an associated model.

        """

        ids = self._get_fit_ids(id, otherids)

        # If an id is given then it must have data but does not have to
        # have a model (to keep with existing behavior).
        #
        # At this point ids is not empty.
        #
        out = []
        for idval in ids:
            data = self.get_data(idval)
            try:
                model = self.get_model(idval)
            except IdentifierErr:
                continue

            out.append(FitStore(idval, data, model))

        # Ensure we have something to fit.
        #
        if len(out) == 0:
            raise IdentifierErr("nomodels")

        return out

    def _get_fit(self,
                 id: Optional[IdType],
                 otherids: Sequence[IdType] = (),
                 estmethod=None, numcores=1
                 ) -> tuple[tuple[IdType, ...], Fit]:
        """Create the fit object for the given identifiers.

        Given the identifiers (the id and otherids arguments), find
        the data and models and return a Fit object.

        Parameters
        ----------
        id : int or str or None
            The identifier to fit. A value of None means all available
            datasets with models.
        otherids : sequence of int or str
            Additional identifiers to fit. Ignored when id is None.
        estmethod : `sherpa.estmethods.EstMethod` or None
            Passed to the Fit object.
        numcores : int, optional
            The number of CPU cores to use (this is used when
            evaluating the models for multiple data sets).

        Returns
        -------
        ids, fit : tuple, `sherpa.fit.Fit` instance
            The datasets used (it may not include all the values from
            id and otherids as those datasets without associated
            models will be skipped) and the fit object.

        Raises
        ------
        IdentifierErr
            If there are no datasets with an associated model.

        """

        store = self._prepare_fit(id, otherids)
        return self._get_fit_obj(store, estmethod, numcores)

    def _get_stat_info(self):
        """Return the stat info structures.

        For each identifier with a dataset and model, calculate the
        current statistics (stored in a sherpa.fit.StatInfoResults
        object), and then - when there are multiple such identifiers -
        a combined result.

        Returns
        -------
        stats : list of `sherpa.Fit.StatInfoResults`

        Raises
        ------
        IdentifierErr
            If there are no datasets with an associated model.

        """

        store = self._prepare_fit(None)

        output = []

        # Report the per-dataset statistics before the combined data
        # (only relevant when there's more than one entry).
        #
        if len(store) > 1:
            for s in store:
                f = Fit(s.data, s.model, self._current_stat)

                statinfo = f.calc_stat_info()
                statinfo.name = f'Dataset {s.idval}'
                statinfo.ids = (s.idval, )

                output.append(statinfo)

        idvals, f = self._get_fit_obj(store, estmethod=None)
        statinfo = f.calc_stat_info()
        statinfo.ids = list(idvals)  # TODO: list or tuple?
        if len(store) == 1:
            statinfo.name = f'Dataset {statinfo.ids}'  # TODO: do we want to use ids[0]?
        else:
            statinfo.name = f'Datasets {statinfo.ids}'

        output.append(statinfo)
        return output

    def calc_stat_info(self):
        """Display the statistic values for the current models.

        Displays the statistic value for each data set, and the
        combined fit, using the current set of models, parameters, and
        ranges. The output is printed to stdout, and so is intended
        for use in interactive analysis. The `get_stat_info`
        function returns the same information but as an array of
        Python structures.

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

        The information returned by `calc_stat_info` includes:

        Dataset
           The dataset identifier (or identifiers).

        Statistic
           The name of the statistic used to calculate the results.

        Fit statistic value
           The current fit statistic value.

        Data points
           The number of bins used in the fit.

        Degrees of freedom
            The number of bins minus the number of thawed parameters.

        Some fields are only returned for a subset of statistics:

        Probability (Q-value)
            A measure of the probability that one would observe the
            reduced statistic value, or a larger value, if the assumed
            model is true and the best-fit model parameters are the true
            parameter values.

        Reduced statistic
            The fit statistic value divided by the number of degrees of
            freedom.

        Examples
        --------

        >>> calc_stat_info()

        """
        output = self._get_stat_info()
        output = [statinfo.format() for statinfo in output]

        if len(output) > 1:
            info('\n\n'.join(output))
        else:
            info(output[0])

    def get_stat_info(self):
        """Return the statistic values for the current models.

        Calculate the statistic value for each data set, and the
        combined fit, using the current set of models, parameters, and
        ranges.

        Returns
        -------
        stats : array of `sherpa.fit.StatInfoResults`
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

        The return value of `get_stat_info` differs to `get_fit_results`
        since it includes values for each data set, individually, rather
        than just the combined results.

        The fields of the object include:

        name
           The name of the data set, or sets, as a string.

        ids
           A sequence of the data set ids (it may be a tuple or
           array) included in the results.

        bkg_ids
           A sequence of the background data set ids (it may be a
           tuple or array) included in the results, if any.

        statname
           The name of the statistic function (as used in `set_stat`).

        statval
           The statistic value.

        numpoints
           The number of bins used in the fits.

        dof
           The number of degrees of freedom in the fit (the number of
           bins minus the number of free parameters).

        qval
           The Q-value (probability) that one would observe the
           reduced statistic value, or a larger value, if the assumed
           model is true and the current model parameters are the
           true parameter values. This will be ``None`` if the value
           can not be calculated with the current statistic (e.g.
           the Cash statistic).

        rstat
           The reduced statistic value (the ``statval`` field divided
           by ``dof``). This is not calculated for all statistics.

        Examples
        --------

        >>> res = get_stat_info()
        >>> res[0].statval
        498.21750663761935
        >>> res[0].dof
        439

        """
        return self._get_stat_info()

    def get_fit_results(self) -> FitResults:
        """Return the results of the last fit.

        This function returns the results from the most-recent fit.
        The returned value includes information on the parameter values
        and fit statistic.

        Returns
        -------
        stats : a `sherpa.fit.FitResults` instance
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

        datasets
           A sequence of the data set ids included in the results.

        itermethodname
           What iterated-fit scheme was used, if any (as set by
           `set_iter_method`).

        statname
           The name of the statistic function (as used in `set_stat`).

        succeeded
           Was the fit successful (did it converge)?

        parnames
           A tuple of the parameter names that were varied in the fit
           (the thawed parameters in the model expression).

        parvals
           A tuple of the parameter values, in the same order as
           ``parnames``.

        statval
           The statistic value after the fit.

        istatval
           The statistic value at the start of the fit.

        dstatval
           The change in the statistic value (``istatval - statval``).

        numpoints
           The number of bins used in the fits.

        dof
           The number of degrees of freedom in the fit (the number of
           bins minus the number of free parameters).

        qval
           The Q-value (probability) that one would observe the
           reduced statistic value, or a larger value, if the assumed
           model is true and the current model parameters are the
           true parameter values. This will be ``None`` if the value
           can not be calculated with the current statistic (e.g.
           the Cash statistic).

        rstat
           The reduced statistic value (the ``statval`` field divided
           by ``dof``). This is not calculated for all statistics.

        message
           A message about the results of the fit (e.g. if the fit was
           unable to converge). The format and contents depend on the
           optimisation method.

        nfev
           The number of model evaluations made during the fit.

        Examples
        --------

        Display the fit results:

        >>> print(get_fit_results())

        Inspect the fit results:

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
        if self._fit_results is None:
            raise SessionErr('nofit', 'fit')

        return self._fit_results

    def guess(self, id=None, model=None, limits=True, values=True):
        """Estimate the parameter values and ranges given the loaded data.

        The guess function can change the parameter values and
        limits to match the loaded data. This is generally limited to
        changing the amplitude and position parameters (sometimes just
        the values and sometimes just the limits). The parameters that
        are changed depend on the type of model.

        .. versionchanged:: 4.17.0
           The guess routine will now work with composite models and
           those which include an instrumental response, such as an
           ARF.  It only works on individual models, so the values,
           and limits, guessed are only approximate.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        model
           Change the parameters of this model component. If ``None``,
           then the source expression is assumed to consist of a single
           component, and that component is used.
        limits : bool
           Should the parameter limits be changed? The default is
           ``True``.
        values : bool
           Should the parameter values be changed? The default is
           ``True``.

        See Also
        --------
        get_default_id : Return the default data set identifier.
        reset : Reset the model parameters to their default settings.
        set_par : Set the value, limits, or behavior of a model parameter.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        The guess function can reduce the time required to fit a
        data set by moving the parameters closer to a realistic
        solution. It can also be useful because it can set bounds on
        the parameter values based on the data: for instance, many
        two-dimensional models will limit their ``xpos`` and ``ypos``
        values to lie within the data area. This can be done manually,
        but `guess` simplifies this, at least for those parameters
        that are supported. Instrument models - such as an ARF and
        RMF - should be set up *before* calling guess.

        Examples
        --------

        Since the source expression contains only one component,
        guess can be called with no arguments:

        >>> set_source(polynom1d.poly)
        >>> guess()

        Both components - that is, gal and pl - will be passed to the
        guess routine, but not all models can have their parameters
        guessed, and there is no attempt to recognize that the models
        are combined (in this case by multiplication):

        >>> set_source(xsphabs.gal * powlaw1d.pl)
        >>> guess()

        In this case, guess is called on each component separately.

        >>> set_source(gauss1d.line + powlaw1d.cont)
        >>> guess(line)
        >>> guess(cont)

        In this example, the values of the ``src`` model component are
        guessed from the "src" data set, whereas the ``bgnd`` component
        is guessed from the "bgnd" data set.

        >>> set_source("src", gauss2d.src + const2d.bgnd)
        >>> set_source("bgnd", bgnd)
        >>> guess("src", src)
        >>> guess("bgnd", bgnd)

        The above could also have been written as:

        >>> guess("src", src)
        >>> guess("bgnd")

        Set the source model for the default dataset. Guess is run to
        determine the values of the model component "p1" and the limits
        of the model component "g1":

        >>> set_source(powlaw1d.p1 + gauss1d.g1)
        >>> guess(p1, limits=False)
        >>> guess(g1, values=False)

        """

        # We can have
        #    id=None, model=None
        #    id=X, model=None
        #       -> is this id=X or model=X
        #    id=X, model=Y
        #
        if model is None:
            # Try to find out if this is a known identifier
            if id not in self.list_data_ids():
                id, model = model, id

        idval = self._fix_id(id)
        kwargs = {'limits': limits, 'values': values}

        if model is not None:
            model = self._check_model(model)
            try:
                model.guess(*self.get_data(idval).to_guess(), **kwargs)
            except NotImplementedError:
                warning('No guess found for %s', model.name)
            return

        ids, f = self._get_fit(idval)
        try:
            f.guess(**kwargs)
        except NotImplementedError:
            warning('No guess found for %s', self.get_model(idval).name)

    def calc_stat(self,
                  id: Optional[IdType] = None,
                  *otherids: IdType):
        """Calculate the fit statistic for a data set.

        Evaluate the model for one or more data sets, compare it to
        the data using the current statistic, and return the value.
        No fitting is done, as the current model parameter, and any
        filters, are used.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        *otherids : int or str, optional
           Other data sets to use in the calculation.

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

        Find the statistic for data set 3:

        >>> stat = calc_stat(3)

        When fitting to multiple data sets, you can get the contribution
        to the total fit statistic from only one data set, or from
        several by listing the datasets explicitly. The following finds
        the contribution from the data sets labelled "core" and "jet":

        >>> stat = calc_stat("core", "jet")

        Calculate the statistic value using two different statistics:

        >>> set_stat('cash')
        >>> s1 = calc_stat()
        >>> set_stat('cstat')
        >>> s2 = calc_stat()

        """
        ids, f = self._get_fit(id, otherids)
        return f.calc_stat()

    def calc_chisqr(self,
                    id: Optional[IdType] = None,
                    *otherids: IdType):
        """Calculate the per-bin chi-squared statistic.

        Evaluate the model for one or more data sets, compare it to the
        data using the current statistic, and return an array of
        chi-squared values for each bin. No fitting is done, as the
        current model parameter, and any filters, are used.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        *otherids : int or str, optional
           Other data sets to use in the calculation.

        Returns
        -------
        chisq : array or ``None``
           The chi-square value for each bin of the data, using the
           current statistic (as set by `set_stat`).  A value of
           ``None`` is returned if the statistic is not a chi-square
           distribution.

        See Also
        --------
        calc_stat : Calculate the fit statistic for a data set.
        calc_stat_info : Display the statistic values for the current models.
        set_stat : Set the statistical method.

        Notes
        -----

        The output array length equals the sum of the arrays lengths
        of the requested data sets.

        Examples
        --------

        When called with no arguments, the return value is the chi-squared
        statistic for each bin in the data sets which have a defined model.

        >>> calc_chisqr()

        Supplying a specific data set ID to calc_chisqr - such as "1" or
        "src" - will return the chi-squared statistic array for only that
        data set.

        >>> calc_chisqr(1)
        >>> calc_chisqr("src")

        Restrict the calculation to just datasets 1 and 3:

        >>> calc_chisqr(1, 3)

        """
        ids, f = self._get_fit(id, otherids)
        return f.calc_chisqr()

    # also in sherpa.astro.utils
    def fit(self,
            id: Optional[IdType] = None,
            *otherids: IdType,
            **kwargs
            ) -> None:
        """Fit a model to one or more data sets.

        Use forward fitting to find the best-fit model to one or more
        data sets, given the chosen statistic and optimization
        method. The fit proceeds until the results converge or the
        number of iterations exceeds the maximum value (these values
        can be changed with `set_method_opt`). An iterative scheme can
        be added using `set_iter_method` to try and improve the
        fit. The final fit results are displayed to the screen and can
        be retrieved with `get_fit_results`.

        .. versionchanged:: 4.17.0
           The outfile parameter can now be sent a Path object or a
           file handle instead of a string.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are fit simultaneously.
        *otherids : int or str, optional
           Other data sets to use in the calculation.
        outfile : str, Path, IO object, or None, optional
           If set, then the fit results will be written to a file with
           this name. The file contains the per-iteration fit results.
        clobber : bool, optional
           This flag controls whether an existing file can be
           overwritten (``True``) or if it raises an exception
           (``False``, the default setting). This is only used if
           `outfile` is set to a string or Path object.

        Raises
        ------
        sherpa.utils.err.FitErr
           If `filename` already exists and `clobber` is ``False``.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
        contour_fit : Contour the fit to a data set.
        covar : Estimate the confidence intervals using the confidence method.
        freeze : Fix model parameters so they are not changed by a fit.
        get_fit_results : Return the results of the last fit.
        plot_fit : Plot the fit results (data, model) for a data set.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        set_stat : Set the statistical method.
        set_method : Change the optimization method.
        set_method_opt : Change an option of the current optimization method.
        set_full_model : Define the convolved model expression for a data set.
        set_iter_method : Set the iterative-fitting scheme used in the fit.
        set_model : Set the model expression for a data set.
        show_fit : Summarize the fit results.
        thaw : Allow model parameters to be varied during a fit.

        Notes
        -----
        If outfile is sent a file handle then it is not closed by this
        routine.

        Examples
        --------

        Simultaneously fit all data sets with models and then
        store the results in the variable fres:

        >>> fit()
        >>> fres = get_fit_results()

        Fit just the data set 'img':

        >>> fit('img')

        Simultaneously fit data sets 1, 2, and 3:

        >>> fit(1, 2, 3)

        Fit data set 'jet' and write the fit results to the text file
        'jet.fit', over-writing it if it already exists:

        >>> fit('jet', outfile='jet.fit', clobber=True)

        Store the per-iteration values in a StringIO object and
        extract the data into the variable txt (this avoids the need
        to create a file):

        >>> from io import StringIO
        >>> out = StringIO()
        >>> fit(outfile=out)
        >>> txt = out.getvalue()

        """
        ids, f = self._get_fit(id, otherids)
        res = f.fit(**kwargs)
        res.datasets = ids
        self._fit_results = res
        info(res.format())

    # Back-compatibility
    # DOC-NOTE: can this be noted as deprecated now?
    simulfit = fit

    #
    # Simulation functions
    #

    def get_pvalue_results(self):
        """Return the data calculated by the last plot_pvalue call.

        The `get_pvalue_results` function returns the likelihood ratio test
        results computed by the `plot_pvalue` command, which compares fits
        of the null model to fits of the alternative model using faked data
        with Poisson noise. The likelihood ratio based on the observed data is
        returned, along with the p-value, used to reject or accept the null
        model.

        .. versionchanged:: 4.17.0
           The "wstat" statistic can now be used with this routine.

        .. versionchanged:: 4.15.1
           The parnames and parvals attributes have been added. They
           are intended to debug problem cases and so are not
           displayed by default.

        Returns
        -------
        plot : None or a `sherpa.sim.simulate.LikelihoodRatioResults` instance
           If `plot_pvalue` or `get_pvalue_plot` have been called then
           the return value is a `sherpa.sim.simulate.LikelihoodRatioResults`
           instance, otherwise `None` is returned.

        See Also
        --------
        plot_value : Compute and plot a histogram of likelihood ratios by simulating data.
        get_pvalue_plot : Return the data used by plot_pvalue.

        Notes
        -----
        The fields of the returned (`LikelihoodRatioResults`) object are:

        ratios
           The calculated likelihood ratio for each iteration.

        stats
           The calculated fit statistics for each iteration, stored
           as the null model and then the alt model in a nsim by 2
           array.

        samples
           The parameter samples array for each simulation, stored in
           a nsim by npar array.

        lr
           The likelihood ratio of the observed data for the null and
           alternate models.

        ppp
           The p value of the observed data for the null and alternate
           models.

        null
           The fit statistic of the null model on the observed data.

        alt
           The fit statistic of the alternate model on the observed data.

        parnames
           The names of the fitted parameters in the alternate
           model. This will be larger than the number of parameters
           returned in the samples field.

        parvals
           The best-fit parameter values to the alternate model for
           each simulation, stored as a nsim by len(parnames) array.

        Examples
        --------

        Return the results of the last pvalue analysis and display the
        results - first using the `format` method, which provides a
        summary of the data, and then a look at the individual fields
        in the returned object. The last call displays the contents
        of one of the fields (`ppp`).

        >>> res = get_pvalue_results()
        >>> print(res.format())
        >>> print(res)
        >>> print(res.ppp)

        Display the ratio values to check they look sensible (such as
        not dropping to a long range of 0's, although this can also
        suggest the alternate model is not preferred to the null
        model):

        >>> plot_trace(res.ratios, name="ratios")

        Look at the cumulative distribution of the ratios:

        >>> plot_cdf(res.ratios, name="ratios")

        The parvals field shows the fitted parameter values for the
        alternate model at each iteration:

        >>> plot_trace(res.parvals[:, 0], name=res.parnames[0])
        >>> plot_trace(res.parvals[:, 1], name=res.parnames[1])

        """
        return self._pvalue_results

    # DOC-TODO: improve discussion of how the simulations are done.
    def plot_pvalue(self, null_model, alt_model, conv_model=None,
                    id: IdType = 1,
                    otherids: Sequence[IdType] = (),
                    num=500, bins=25, numcores=None,
                    replot=False, overplot=False, clearwindow=True,
                    **kwargs):
        """Compute and plot a histogram of likelihood ratios by simulating data.

        Compare the likelihood of the null model to an alternative model
        by running a number of simulations to calibrate the likelihood ratio test statistics.
        The distribution of the simulated likelihood ratios is plotted and compared to the likelihoods
        of the two models fit to the observed data. The fit
        statistic must be set to a likelihood-based method, such
        as "cash" or "cstat". Screen output is created as well as the plot; these values
        can be retrieved with `get_pvalue_results`.

        The algorithm is based on the description in Sec.5.2 in "Statistics,
        Handle with Care: Detecting Multiple Model Components with the Likelihood Ratio Test"
        by Protassov et al., 2002, The Astrophysical Journal, 571, 545; <doi:10.1086/339856>

        .. versionchanged:: 4.17.0
           The "wstat" statistic can now be used with this routine.

        Parameters
        ----------
        null_model
           The model expression for the null hypothesis.
        alt_model
           The model expression for the alternative hypothesis.
        conv_model : optional
           An expression used to modify the model so that it can be
           compared to the data (e.g. a PSF or PHA response).
        id : int or str, optional
           The data set that provides the data. The default is 1.
        otherids : sequence of int or str, optional
           Other data sets to use in the calculation.
        num : int, optional
           The number of simulations to run. The default is 500.
        bins : int, optional
           The number of bins to use to create the histogram. The
           default is 25.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_pvalue`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        Raises
        ------
        TypeError
           An invalid statistic.

        See Also
        --------
        get_pvalue_plot : Return the data used by plot_pvalue.
        get_pvalue_results : Return the data calculated by the last plot_pvalue call.

        Notes
        -----
        Each simulation involves creating a data set using the observed
        data simulated with Poisson noise.

        For the likelihood ratio test to be valid, the following
        conditions must hold:

           1. The null model is nested within the alternative model.

           2. The extra parameters of the alternative model have
              Gaussian (normal) distributions that are not truncated
              by the boundaries of the parameter spaces.

        Examples
        --------

        Use the likelihood ratio to see if the data in data set 1 has
        a statistically-significant gaussian component:

        >>> create_model_component('powlaw1d', 'pl')
        >>> create_model_component('gauss1d', 'gline')
        >>> plot_pvalue(pl, pl + gline)

        Use 1000 simulations and use the data from data sets
        'core', 'jet1', and 'jet2':

        >>> mdl1 = pl
        >>> mdl2 = pl + gline
        >>> plot_pvalue(mdl1, mdl2, id='core', otherids=('jet1', 'jet2'),
        ...             num=1000)

        Apply a convolution to the models before fitting:

        >>> rsp = get_psf()
        >>> plot_pvalue(mdl1, mdl2, conv_model=rsp)

        """

        lrplot = self.get_pvalue_plot(null_model=null_model, alt_model=alt_model,
                                      conv_model=conv_model, id=id, otherids=otherids,
                                      num=num, bins=bins, numcores=numcores,
                                      recalc=not replot)
        self._plot(lrplot, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def get_pvalue_plot(self, null_model=None, alt_model=None, conv_model=None,
                        id: IdType = 1,
                        otherids: Sequence[IdType] = (),
                        num=500, bins=25, numcores=None,
                        recalc=False):
        """Return the data used by plot_pvalue.

        Access the data arrays and preferences defining the histogram plot
        produced by the `plot_pvalue` function, a histogram of the likelihood
        ratios comparing fits of the null model to fits of the alternative
        model using faked data with Poisson noise. Data returned includes the
        likelihood ratio computed using the observed data, and the p-value,
        used to reject or accept the null model.

        .. versionchanged:: 4.17.0
           The "wstat" statistic can now be used with this routine.

        Parameters
        ----------
        null_model
           The model expression for the null hypothesis.
        alt_model
           The model expression for the alternative hypothesis.
        conv_model : optional
           An expression used to modify the model so that it can be
           compared to the data (e.g. a PSF or PHA response).
        id : int or str, optional
           The data set that provides the data. The default is 1.
        otherids : sequence of int or str, optional
           Other data sets to use in the calculation.
        num : int, optional
           The number of simulations to run. The default is 500.
        bins : int, optional
           The number of bins to use to create the histogram. The
           default is 25.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.
        recalc : bool, optional
           The default value (``False``) means that the results from the
           last call to `plot_pvalue` or `get_pvalue_plot` are
           returned. If ``True``, the values are re-calculated.

        Returns
        -------
        plot : a `sherpa.plot.LRHistogram` instance

        See Also
        --------
        get_pvalue_results, plot_pvalue, set_rng

        Notes
        -----
        The `set_rng` routine is used to control how the random
        numbers are generated.

        Examples
        --------

        Return the values from the last call to `plot_pvalue`:

        >>> pvals = get_pvalue_plot()
        >>> pvals.ppp
        0.472

        Run 500 simulations for the two models and print the
        results:

        >>> pvals = get_pvalue_plot(mdl1, mdl2, recalc=True, num=500)
        >>> print(pvals)

        """

        lrplot = self._lrplot
        if not recalc:
            return lrplot

        if null_model is None:
            raise TypeError("null model cannot be None")

        if alt_model is None:
            raise TypeError("alternative model cannot be None")

        ids, fit = self._get_fit(id, otherids)

        pvalue = sherpa.sim.LikelihoodRatioTest.run
        results = pvalue(fit, null_model, alt_model,
                         conv_mdl=conv_model,
                         stat=self._current_stat,
                         method=self._current_method,
                         niter=num,
                         numcores=numcores,
                         rng=self.get_rng())

        info(results.format())
        self._pvalue_results = results

        lrplot.prepare(ratios=results.ratios, bins=bins,
                       niter=num, lr=results.lr, ppp=results.ppp)
        return lrplot

    #
    # Sampling functions
    #
    # DOC-TODO: This copies from sherpa/sim/sample.py, so can
    # docstrings be shared (whether directly or via functools)?
    # Unfortunately not quite a direct copy, so hard
    # to see how to do

    def normal_sample(self, num=1, sigma=1, correlate=True,
                      id: Optional[IdType] = None,
                      otherids: Sequence[IdType] = (),
                      numcores=None):
        """Sample the fit statistic by taking the parameter values
        from a normal distribution.

        For each iteration (sample), change the thawed parameters by
        drawing values from a uni- or multi-variate normal (Gaussian)
        distribution, and calculate the fit statistic.

        Parameters
        ----------
        num : int, optional
           The number of samples to use (default is 1).
        sigma : number, optional
           The width of the normal distribution (the default
           is 1).
        correlate : bool, optional
           Should a multi-variate normal be used, with parameters
           set by the covariance matrix (``True``) or should a
           uni-variate normal be used (``False``)?
        id : int, str, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : sequence of int or str, optional
           Other data sets to use in the calculation.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        samples
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

    # DOC-TODO: improve the description of factor parameter
    def uniform_sample(self, num=1, factor=4,
                       id: Optional[IdType] = None,
                       otherids: Sequence[IdType] = (),
                       numcores=None):
        """Sample the fit statistic by taking the parameter values
        from an uniform distribution.

        For each iteration (sample), change the thawed parameters by
        drawing values from a uniform distribution, and calculate the
        fit statistic.

        Parameters
        ----------
        num : int, optional
           The number of samples to use (default is 1).
        factor : number, optional
           Multiplier to expand the scale parameter (default is 4).
        id : int, str, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : sequence of int or str, optional
           Other data sets to use in the calculation.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        samples
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

    def t_sample(self, num=1, dof=None,
                 id: Optional[IdType] = None,
                 otherids: Sequence[IdType] = (),
                 numcores=None):
        """Sample the fit statistic by taking the parameter values from
        a Student's t-distribution.

        For each iteration (sample), change the thawed parameters
        by drawing values from a Student's t-distribution, and
        calculate the fit statistic.

        Parameters
        ----------
        num : int, optional
           The number of samples to use (default is 1).
        dof : optional
           The number of degrees of freedom to use (the default
           is to use the number from the current fit).
        id : int, str, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : sequence of int or str, optional
           Other data sets to use in the calculation.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        samples
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

    # DOC-TODO: how best to document the settings?
    # DOC-TODO: have I got soft_limits described correctly?
    def get_covar(self):
        """Return the covariance estimation object.

        Returns
        -------
        covar : object

        See Also
        --------
        covar : Estimate parameter confidence intervals using the covariance method.
        get_covar_opt : Return one or all of the options for the covariance method.
        set_covar_opt : Set an option of the covar estimation object.

        Notes
        -----
        The attributes of the covariance object include:

        ``eps``
           The precision of the calculated limits. The default is
           0.01.

        ``maxiters``
           The maximum number of iterations allowed before stopping
           for that parameter. The default is 200.

        ``sigma``
           What is the error limit being calculated. The default is
           1.

        ``soft_limits``
           Should the search be restricted to the soft limits of the
           parameters (``True``), or can parameter values go out all the
           way to the hard limits if necessary (``False``).  The default
           is ``False``

        Examples
        --------

        >>> print(get_covar())
        name        = covariance
        sigma       = 1
        maxiters    = 200
        soft_limits = False
        eps         = 0.01

        Change the ``sigma`` field to 1.9.

        >>> cv = get_covar()
        >>> cv.sigma = 1.6

        """
        return self._estmethods['covariance']

    # DOC-TODO: how best to document the settings?
    # DOC-TODO: have I got soft_limits described correctly?
    # DOC-TODO: when verbose=True how is extra output displayed?
    # stdout, stderr, sherpa logging?
    def get_conf(self):
        """Return the confidence-interval estimation object.

        Returns
        -------
        conf : object

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
        get_conf_opt : Return one or all of the options for the confidence interval method.
        set_conf_opt : Set an option of the conf estimation object.

        Notes
        -----
        The attributes of the confidence-interval object include:

        ``eps``
           The precision of the calculated limits. The default is
           0.01.

        ``fast``
           If ``True`` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is ``False``.

        ``max_rstat``
           If the reduced chi square is larger than this value, do not
           use (only used with chi-square statistics). The default is
           3.

        ``maxfits``
           The maximum number of re-fits allowed (that is, when the
           ``remin`` filter is met). The default is 5.

        ``maxiters``
           The maximum number of iterations allowed when bracketing
           limits, before stopping for that parameter. The default is
           200.

        ``numcores``
           The number of computer cores to use when evaluating results
           in parallel. This is only used if ``parallel`` is ``True``.
           The default is to use all cores.

        ``openinterval``
           How the `conf` method should cope with intervals that do
           not converge (that is, when the ``maxiters`` limit has been
           reached). The default is ``False``.

        ``parallel``
           If there is more than one free parameter then the results
           can be evaluated in parallel, to reduce the time required.
           The default is ``True``.

        ``remin``
           The minimum difference in statistic value for a new fit
           location to be considered better than the current best fit
           (which starts out as the starting location of the fit at
           the time `conf` is called). The default is 0.01.

        ``sigma``
           What is the error limit being calculated. The default is
           1.

        ``soft_limits``
           Should the search be restricted to the soft limits of the
           parameters (``True``), or can parameter values go out all the
           way to the hard limits if necessary (``False``).  The default
           is ``False``

        ``tol``
           The tolerance for the fit. The default is 0.2.

        ``verbose``
           Should extra information be displayed during fitting?
           The default is ``False``.

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

        Change the ``remin`` field to 0.05.

        >>> cf = get_conf()
        >>> cf.remin = 0.05

        """
        return self._estmethods['confidence']

    def get_proj(self):
        """Return the confidence-interval estimation object.

        .. note:: The `conf` function should be used instead of `proj`.

        Returns
        -------
        proj : object

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
        get_proj_opt : Return one or all of the options for the confidence interval method.
        proj : Estimate confidence intervals for fit parameters.
        set_proj_opt : Set an option of the proj estimation object.

        Notes
        -----
        The attributes of the object include:

        ``eps``
           The precision of the calculated limits. The default is
           0.01.

        ``fast``
           If ``True`` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is ``False``.

        ``max_rstat``
           If the reduced chi square is larger than this value, do not
           use (only used with chi-square statistics). The default is
           3.

        ``maxfits``
           The maximum number of re-fits allowed (that is, when the
           ``remin`` filter is met). The default is 5.

        ``maxiters``
           The maximum number of iterations allowed when bracketing
           limits, before stopping for that parameter. The default is
           200.

        ``numcores``
           The number of computer cores to use when evaluating results
           in parallel. This is only used if ``parallel`` is ``True``.
           The default is to use all cores.

        ``parallel``
           If there is more than one free parameter then the results
           can be evaluated in parallel, to reduce the time required.
           The default is ``True``.

        ``remin``
           The minimum difference in statistic value for a new fit
           location to be considered better than the current best fit
           (which starts out as the starting location of the fit at
           the time `proj` is called). The default is 0.01.

        ``sigma``
           What is the error limit being calculated. The default is
           1.

        ``soft_limits``
           Should the search be restricted to the soft limits of the
           parameters (``True``), or can parameter values go out all the
           way to the hard limits if necessary (``False``).  The default
           is ``False``

        ``tol``
           The tolerance for the fit. The default is 0.2.

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
        _check_str_type(optname, "optname")
        if optname not in estmethod.config:
            raise ArgumentErr('badopt', optname, estmethod.name)

    def _get_estmethod_opt(self, methodname, optname=None):
        meth = self._estmethods.get(methodname.lower())
        if meth is None:
            raise ArgumentErr('badconf', methodname)
        if optname is None:
            return meth.config
        self._check_estmethod_opt(meth, optname)
        return meth.config[optname]

    def _set_estmethod_opt(self, methodname, optname, val):
        meth = self._estmethods.get(methodname.lower())
        if meth is None:
            raise ArgumentErr('badconf', methodname)
        self._check_estmethod_opt(meth, optname)
        meth.config[optname] = val

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
           If the ``name`` argument is not recognized.

        See Also
        --------
        covar : Estimate parameter confidence intervals using the covariance method.
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
           If the ``name`` argument is not recognized.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
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
           If the ``name`` argument is not recognized.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
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
        val
           The new value for the option.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the ``name`` argument is not recognized.

        See Also
        --------
        covar : Estimate parameter confidence intervals using the covariance method.
        get_covar : Return the covar estimation object.
        get_covar_opt : Return one or all options of the covar estimation object.

        Examples
        --------

        >>> set_covar_opt('sigma', 1.6)

        """
        self._set_estmethod_opt('covariance', name, val)

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
        val
           The new value for the option.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the ``name`` argument is not recognized.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
        get_conf : Return the conf estimation object.
        get_conf_opt : Return one or all options of the conf estimation object.

        Examples
        --------

        >>> set_conf_opt('parallel', False)

        """
        self._set_estmethod_opt('confidence', name, val)

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
        val
           The new value for the option.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the ``name`` argument is not recognized.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
        proj : Estimate parameter confidence intervals using the projection method.
        get_proj : Return the proj estimation object.
        get_proj_opt : Return one or all options of the proj estimation object.

        Examples
        --------

        >>> set_proj_opt('parallel', False)

        """
        self._set_estmethod_opt('projection', name, val)

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

        ``datasets``
           A tuple of the data sets used in the analysis.

        ``methodname``
           This will be 'covariance'.

        ``iterfitname``
           The name of the iterated-fit method used, if any.

        ``fitname``
           The name of the optimization method used.

        ``statname``
           The name of the fit statistic used.

        ``sigma``
           The sigma value used to calculate the confidence
           intervals.

        ``percent``
           The percentage of the signal contained within the
           confidence intervals (calculated from the ``sigma``
           value assuming a normal distribution).

        ``parnames``
           A tuple of the parameter names included in the analysis.

        ``parvals``
           A tuple of the best-fit parameter values, in the same
           order as ``parnames``.

        ``parmins``
           A tuple of the lower error bounds, in the same
           order as ``parnames``.

        ``parmaxes``
           A tuple of the upper error bounds, in the same
           order as ``parnames``.

        ``nfits``
           The number of model evaluations.

        There is also an ``extra_output`` field which is used to return
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

        >>> res.extra_output
        array([[ 6.19847635]])

        """
        if self._covariance_results is None:
            raise SessionErr('noaction', "covariance")

        return self._covariance_results

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

        ``datasets``
           A tuple of the data sets used in the analysis.

        ``methodname``
           This will be 'confidence'.

        ``iterfitname``
           The name of the iterated-fit method used, if any.

        ``fitname``
           The name of the optimization method used.

        ``statname``
           The name of the fit statistic used.

        ``sigma``
           The sigma value used to calculate the confidence
           intervals.

        ``percent``
           The percentage of the signal contained within the
           confidence intervals (calculated from the ``sigma``
           value assuming a normal distribution).

        ``parnames``
           A tuple of the parameter names included in the analysis.

        ``parvals``
           A tuple of the best-fit parameter values, in the same
           order as ``parnames``.

        ``parmins``
           A tuple of the lower error bounds, in the same
           order as ``parnames``.

        ``parmaxes``
           A tuple of the upper error bounds, in the same
           order as ``parnames``.

        ``nfits``
           The number of model evaluations.

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
        >>> pvals2 = [(v, v+l, v+h) for (v, l, h) in pvals1]
        >>> dres = dict(zip(res.parnames, pvals2))
        >>> dres['p1.gamma']
        (2.1585155113403327, 2.07572994399221, 2.241926145484433)

        """
        if self._confidence_results is None:
            raise SessionErr('noaction', "confidence")

        return self._confidence_results

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
        conf : Estimate parameter confidence intervals using the confidence method.
        proj : Estimate confidence intervals for fit parameters.
        get_proj_opt : Return one or all of the options for the projection method.
        set_proj_opt : Set an option of the proj estimation object.

        Notes
        -----
        The fields of the object include:

        ``datasets``
           A tuple of the data sets used in the analysis.

        ``methodname``
           This will be 'projection'.

        ``iterfitname``
           The name of the iterated-fit method used, if any.

        ``fitname``
           The name of the optimization method used.

        ``statname``
           The name of the fit statistic used.

        ``sigma``
           The sigma value used to calculate the confidence
           intervals.

        ``percent``
           The percentage of the signal contained within the
           confidence intervals (calculated from the ``sigma``
           value assuming a normal distribution).

        ``parnames``
           A tuple of the parameter names included in the analysis.

        ``parvals``
           A tuple of the best-fit parameter values, in the same
           order as ``parnames``.

        ``parmins``
           A tuple of the lower error bounds, in the same
           order as ``parnames``.

        ``parmaxes``
           A tuple of the upper error bounds, in the same
           order as ``parnames``.

        ``nfits``
           The number of model evaluations.

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
        if self._projection_results is None:
            raise SessionErr('noaction', "projection")

        return self._projection_results

    def _est_errors(self, args, methodname):
        """Evaluate the errors for the given arguments.

        The formatted output of the estimation is logged at the
        info level.

        Parameters
        ----------
        args : sequence of sherpa.models.Parameter, sherpa.models.Model, int, or str
            The dataset (when an integer or str) to evaluate, the
            model parameter to apply the error estimate on (must be
            thawed), or a model from which all the thawed parameters
            are used.
        methodname : str
            One of the valid error estimates (sub-classes of
            sherpa.estmethods.EstMethod).

        Returns
        -------
        result : sherpa.fit.ErrorEstResults instance

        """

        id = None
        parlist = []
        otherids = []
        for arg in args:
            if isinstance(arg, sherpa.models.Parameter):
                if arg.frozen:
                    raise ParameterErr('frozen', arg.fullname)

                parlist.append(arg)
                continue

            if isinstance(arg, Model):
                norig = len(parlist)
                for par in arg.pars:
                    if not par.frozen:
                        parlist.append(par)

                # If there were no free parameters then error out.
                # Should this be a ParameterErr or a ModelErr? Neither
                # have an existing label for this case. Pick ParameterErr
                # to match the case when a single parameter is frozen.
                #
                if len(parlist) == norig:
                    emsg = f"Model '{arg.name}' has no thawed parameters"
                    raise ParameterErr(emsg)

                continue

            if id is None:
                id = arg
            else:
                otherids.append(arg)

        if len(parlist) == 0:
            parlist = None

        ids, f = self._get_fit(id, otherids, self._estmethods[methodname])
        res = f.est_errors(self._methods, parlist)
        res.datasets = ids
        info(res.format())
        return res

    # DOC-TODO: include screen output of covar() ?
    def covar(self, *args):
        """Estimate parameter confidence intervals using the covariance method.

        The `covar` command computes confidence interval bounds for
        the specified model parameters in the dataset, using the
        covariance matrix of the statistic. The `get_covar` and
        `set_covar_opt` commands can be used to configure the error
        analysis; an example being changing the ``sigma`` field to 1.6
        (i.e. 90%) from its default value of 1.  The output from the
        routine is displayed on screen, and the `get_covar_results`
        routine can be used to retrieve the results.

        Parameters
        ----------
        id : int or str, optional
           The data set, or sets, that provides the data. If not given
           then all data sets with an associated model are used
           simultaneously.
        parameter : sherpa.models.parameter.Parameter, optional
           The default is to calculate the confidence limits on all
           thawed parameters of the model, or models, for all the data
           sets. The evaluation can be restricted by listing the
           parameters to use. Note that each parameter should be given
           as a separate argument, rather than as a list.  For example
           ``covar(g1.ampl, g1.sigma)``.
        model : sherpa.models.model.Model, optional
           Select all the thawed parameters in the model.

        See Also
        --------
        covar : Estimate the confidence intervals using the confidence method.
        get_covar : Return the covariance estimation object.
        get_covar_results : Return the results of the last `covar` run.
        int_proj : Plot the statistic value as a single parameter is varied.
        int_unc : Plot the statistic value as a single parameter is varied.
        reg_proj : Plot the statistic value as two parameters are varied.
        reg_unc : Plot the statistic value as two parameters are varied.
        set_covar_opt : Set an option of the `covar` estimation object.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with multiple ``ids`` or ``parameters`` values, the
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
        with the ``sigma`` option to `set_covar_opt` or `get_covar`.

        Examples
        --------

        Evaluate confidence intervals for all thawed parameters in all
        data sets with an associated source model. The results are
        then stored in the variable ``res``.

        >>> covar()
        >>> res = get_covar_results()

        Only evaluate the parameters associated with data set 2.

        >>> covar(2)

        Only evaluate the intervals for the ``pos.xpos`` and ``pos.ypos``
        parameters:

        >>> covar(pos.xpos, pos.ypos)

        Change the limits to be 1.6 sigma (90%) rather than the default
        1 sigma.

        >>> set_covar_ope('sigma', 1.6)
        >>> covar()

        Only evaluate the ``clus.kt`` parameter for the data sets with
        identifiers "obs1", "obs5", and "obs6". This will still use
        the 1.6 sigma setting from the previous run.

        >>> covar("obs1", "obs5", "obs6", clus.kt)

        Estimate the errors for all the thawed parameters from the
        ``line`` model and the ``clus.kt`` parameter for datasets 1,
        3, and 4:

        >>> covar(1, 3, 4, line, clus.kt)

        """
        self._covariance_results = self._est_errors(args, 'covariance')

    # DOC-TODO: include screen output of conf() ?
    def conf(self, *args):
        """Estimate parameter confidence intervals using the confidence method.

        The `conf` command computes confidence interval bounds for the
        specified model parameters in the dataset.  A given
        parameter's value is varied along a grid of values while the
        values of all the other thawed parameters are allowed to float
        to new best-fit values. The `get_conf` and `set_conf_opt`
        commands can be used to configure the error analysis; an
        example being changing the 'sigma' field to ``1.6`` (i.e. 90%)
        from its default value of ``1``.  The output from the routine is
        displayed on screen, and the `get_conf_results` routine can be
        used to retrieve the results.

        Parameters
        ----------
        id : int or str, optional
           The data set, or sets, that provides the data. If not given
           then all data sets with an associated model are used
           simultaneously.
        parameter : sherpa.models.parameter.Parameter, optional
           The default is to calculate the confidence limits on all
           thawed parameters of the model, or models, for all the data
           sets. The evaluation can be restricted by listing the
           parameters to use. Note that each parameter should be given
           as a separate argument, rather than as a list.  For example
           ``conf(g1.ampl, g1.sigma)``.
        model : sherpa.models.model.Model, optional
           Select all the thawed parameters in the model.

        See Also
        --------
        covar : Estimate the confidence intervals using the covariance method.
        get_conf : Return the confidence-interval estimation object.
        get_conf_results : Return the results of the last `conf` run.
        int_proj : Plot the statistic value as a single parameter is varied.
        int_unc : Plot the statistic value as a single parameter is varied.
        reg_proj : Plot the statistic value as two parameters are varied.
        reg_unc : Plot the statistic value as two parameters are varied.
        set_conf_opt : Set an option of the `conf` estimation object.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with multiple ``ids`` or ``parameters`` values, the
        order is unimportant, since any argument that is not defined
        as a model parameter is assumed to be a data id.

        The `conf` function is different to `covar`, in that in that
        all other thawed parameters are allowed to float to new
        best-fit values, instead of being fixed to the initial
        best-fit values as they are in `covar`.  While `conf` is more
        general (e.g. allowing the user to examine the parameter space
        away from the best-fit point), it is in the strictest sense no
        more accurate than `covar` for determining confidence
        intervals.

        The `conf` function is a replacement for the `proj` function,
        which uses a different algorithm to estimate parameter
        confidence limits.

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
        analysis. This can be changed be varying the ``numcores`` option
        - or setting ``parallel`` to ``False`` - either with
        `set_conf_opt` or `get_conf`.

        As `conf` estimates intervals for each parameter
        independently, the relationship between sigma and the change
        in statistic value delta_S can be particularly simple: sigma =
        the square root of delta_S for statistics sampled from the
        chi-square distribution and for the Cash statistic, and is
        approximately equal to the square root of (2 * delta_S) for
        fits based on the general log-likelihood. The default setting
        is to calculate the one-sigma interval, which can be changed
        with the ``sigma`` option to `set_conf_opt` or `get_conf`.

        The limit calculated by `conf` is basically a 1-dimensional
        root in the translated coordinate system (translated by the
        value of the statistic at the minimum plus sigma^2).  The
        Taylor series expansion of the multi-dimensional function at
        the minimum is::

           f(x + dx) ~ f(x) + grad( f(x) )^T dx + dx^T Hessian( f(x) ) dx + ...

        where x is understood to be the n-dimensional vector
        representing the free parameters to be fitted and the
        super-script 'T' is the transpose of the row-vector. At or
        near the minimum, the gradient of the function is zero or
        negligible, respectively. So the leading term of the expansion
        is quadratic.  The best root finding algorithm for a curve
        which is approximately parabolic is Muller's method [1].
        Muller's method is a generalization of the secant method [2]:
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
        the option ``openinterval`` is set to ``False`` (which is the
        default), the routine will print a warning message about not
        able to find the root within the set tolerance and the
        function will return the average of the open interval which
        brackets the root. If ``openinterval`` is set to ``True`` then
        `conf` will print the minimal open interval which brackets the
        root (not to be confused with the lower and upper bound of the
        confidence interval). The most accurate thing to do is to
        return an open interval where the root is localized/bracketed
        rather then the average of the open interval (since the
        average of the interval is not a root within the specified
        tolerance).

        References
        ----------

        1. Muller, David E., "A Method for Solving Algebraic
           Equations Using an Automatic Computer," MTAC, 10
           (1956), 208-215.

        2. Numerical Recipes in Fortran, 2nd edition, 1986, Press
           et al., p. 347

        Examples
        --------

        Evaluate confidence intervals for all thawed parameters in all
        data sets with an associated source model. The results are
        then stored in the variable ``res``.

        >>> conf()
        >>> res = get_conf_results()

        Only evaluate the parameters associated with data set 2:

        >>> conf(2)

        Only evaluate the intervals for the ``pos.xpos`` and ``pos.ypos``
        parameters:

        >>> conf(pos.xpos, pos.ypos)

        Change the limits to be 1.6 sigma (90%) rather than the default
        1 sigma.

        >>> set_conf_opt('sigma', 1.6)
        >>> conf()

        Only evaluate the ``clus.kt`` parameter for the data sets with
        identifiers "obs1", "obs5", and "obs6". This will still use
        the 1.6 sigma setting from the previous run.

        >>> conf("obs1", "obs5", "obs6", clus.kt)

        Only use two cores when evaluating the errors for the parameters
        used in the model for data set 3:

        >>> set_conf_opt('numcores', 2)
        >>> conf(3)

        Estimate the errors for all the thawed parameters from the
        ``line`` model and the ``clus.kt`` parameter for datasets 1,
        3, and 4:

        >>> conf(1, 3, 4, line, clus.kt)

        """
        self._confidence_results = self._est_errors(args, 'confidence')

    # DOC-TODO: add a deprecation note?
    def proj(self, *args):
        """Estimate parameter confidence intervals using the projection method.

        .. note:: The `conf` function should be used instead of `proj`.

        The `proj` command computes confidence interval bounds for the
        specified model parameters in the dataset. A given parameter's
        value is varied along a grid of values while the values of all
        the other thawed parameters are allowed to float to new
        best-fit values. The `get_proj` and `set_proj_opt`
        commands can be used to configure the error analysis; an
        example being changing the 'sigma' field to ``1.6`` (i.e. 90%)
        from its default value of ``1``.  The output from the routine is
        displayed on screen, and the `get_proj_results` routine can be
        used to retrieve the results.

        Parameters
        ----------
        id : int or str, optional
           The data set, or sets, that provides the data. If not given
           then all data sets with an associated model are used
           simultaneously.
        parameter : sherpa.models.parameter.Parameter, optional
           The default is to calculate the confidence limits on all
           thawed parameters of the model, or models, for all the data
           sets. The evaluation can be restricted by listing the
           parameters to use. Note that each parameter should be given
           as a separate argument, rather than as a list.  For example
           ``proj(g1.ampl, g1.sigma)``.
        model : sherpa.models.model.Model, optional
           Select all the thawed parameters in the model.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        get_proj : Return the confidence-interval estimation object.
        get_proj_results : Return the results of the last `proj` run.
        int_proj : Plot the statistic value as a single parameter is varied.
        reg_proj : Plot the statistic value as two parameters are varied.
        set_proj_opt : Set an option of the `proj` estimation object.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with multiple ``ids`` or ``parameters`` values, the
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
        analysis. This can be changed be varying the ``numcores`` option
        - or setting ``parallel`` to ``False`` - either with
        `set_proj_opt` or `get_proj`.

        As `proj` estimates intervals for each parameter
        independently, the relationship between sigma and the change
        in statistic value delta_S can be particularly simple: sigma =
        the square root of delta_S for statistics sampled from the
        chi-square distribution and for the Cash statistic, and is
        approximately equal to the square root of (2 * delta_S) for
        fits based on the general log-likelihood. The default setting
        is to calculate the one-sigma interval, which can be changed
        with the ``sigma`` option to `set_proj_opt` or `get_proj`.

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
    ###########################################################################

    def set_sampler_opt(self, opt, value):
        """Set an option for the current MCMC sampler.

        Parameters
        ----------
        opt : str
           The option to change. Use `get_sampler` to view the
           available options for the current sampler.
        value
           The value for the option.

        See Also
        --------
        get_sampler : Return the current MCMC sampler options.
        set_prior: Set the prior function to use with a parameter.
        set_sampler : Set the MCMC sampler.

        Notes
        -----
        The options depend on the sampler. The options include:

        defaultprior
           Set to ``False`` when the default prior (flat, between the
           parameter's soft limits) should not be used. Use
           `set_prior` to set the form of the prior for each
           parameter.

        inv
           A bool, or array of bools, to indicate which parameter is
           on the inverse scale.

        log
           A bool, or array of bools, to indicate which parameter is
           on the logarithm (natural log) scale.

        original
           A bool, or array of bools, to indicate which parameter is
           on the original scale.

        p_M
           The proportion of jumps generated by the Metropolis
           jumping rule.

        priorshape
           An array of bools indicating which parameters have a
           user-defined prior functions set with `set_prior`.

        scale
           Multiply the output of `covar` by this factor and
           use the result as the scale of the t-distribution.

        Examples
        --------

        >>> set_sampler_opt('scale', 3)

        """
        self._pyblocxs.set_sampler_opt(opt, value)

    def get_sampler_opt(self, opt):
        """Return an option of the current MCMC sampler.

        Returns
        -------
        opt : str
           The name of the option. The fields depend on the current
           sampler.

        See Also
        --------
        get_sampler : Return the current MCMC sampler options.
        set_sampler_opt : Set an option for the current MCMC sampler.

        Examples
        --------

        >>> get_sampler_opt('log')
        False

        """
        return self._pyblocxs.get_sampler_opt(opt)

    def get_sampler_name(self):
        """Return the name of the current MCMC sampler.

        Returns
        -------
        name : str

        See Also
        --------
        get_sampler : Return the current MCMC sampler options.
        set_sampler : Set the MCMC sampler.

        Examples
        --------

        >>> get_sampler_name()
        'MetropolisMH'

        """
        return self._pyblocxs.get_sampler_name()

    def set_sampler(self, sampler):
        """Set the MCMC sampler.

        The sampler determines the type of jumping rule to
        be used when running the MCMC analysis.

        Parameters
        ----------
        sampler : str or `sherpa.sim.Sampler` instance
           When a string, the name of the sampler to use (case
           insensitive). The supported options are given by the
           `list_samplers` function.

        See Also
        --------
        get_draws : Run the pyBLoCXS MCMC algorithm.
        list_samplers : List the MCMC samplers.
        set_sampler : Set the MCMC sampler.
        set_sampler_opt : Set an option for the current MCMC sampler.

        Notes
        -----
        The jumping rules are:

        MH
           The Metropolis-Hastings rule, which jumps from the best-fit
           location, even if the previous iteration had moved away from
           it.

        MetropolisMH
           This is the Metropolis with Metropolis-Hastings algorithm,
           that jumps from the best-fit with probability ``p_M``,
           otherwise it jumps from the last accepted jump. The
           value of ``p_M`` can be changed using `set_sampler_opt`.

        PragBayes
           This is used when the effective area calibration
           uncertainty is to be included in the calculation. At each
           nominal MCMC iteration, a new calibration product is
           generated, and a series of N (the ``nsubiters`` option) MCMC
           sub-iteration steps are carried out, choosing between
           Metropolis and Metropolis-Hastings types of samplers with
           probability ``p_M``.  Only the last of these sub-iterations
           are kept in the chain.  The ``nsubiters`` and ``p_M`` values
           can be changed using `set_sampler_opt`.

        FullBayes
           Another sampler for use when including uncertainties due
           to the effective area.

        Examples
        --------

        >>> set_sampler('metropolismh')

        """
        self._pyblocxs.set_sampler(sampler)

    def get_sampler(self):
        """Return the current MCMC sampler options.

        Returns the options for the current pyBLoCXS MCMC sampling
        method (jumping rules).

        Returns
        -------
        options : dict
           A copy of the  options for the chosen sampler.  Use
           `set_sampler_opt` to change these values. The fields depend
           on the current sampler.

        See Also
        --------
        get_sampler_name : Return the name of the current MCMC sampler.
        get_sampler_opt : Return an option of the current MCMC sampler.
        set_sampler : Set the MCMC sampler.
        set_sampler_opt : Set an option for the current MCMC sampler.

        Examples
        --------

        >>> print(get_sampler())

        """
        return self._pyblocxs.get_sampler()

    # DOC-TODO: should set_sampler_opt be mentioned here?
    def set_prior(self, par, prior):
        """Set the prior function to use with a parameter.

        The default prior used by `get_draws` for each parameter
        is flat, varying between the soft minimum and maximum
        values of the parameter (as given by the `min` and
        `max` attributes of the parameter object). The `set_prior`
        function is used to change the form of the prior for a
        parameter, and `get_prior` returns the current prior for
        a parameter.

        Parameters
        ----------
        par : a `sherpa.models.parameter.Parameter` instance
           A parameter of a model instance.
        prior : function or sherpa.models.model.Model instance
           The function to use for a prior. It must accept a
           single argument and return a value of the same size
           as the input.

        See Also
        --------
        get_draws : Run the pyBLoCXS MCMC algorithm.
        get_prior : Return the prior function for a parameter (MCMC).
        set_sampler : Set the MCMC sampler.

        Examples
        --------

        Set the prior for the ``kT`` parameter of the ``therm`` component
        to be a gaussian, centered on 1.7 keV and with a FWHM of 0.35
        keV:

        >>> create_model_component('xsapec', 'therm')
        >>> create_model_component('gauss1d', 'p_temp')
        >>> p_temp.pos = 1.7
        >>> p_temp.fwhm = 0.35
        >>> set_prior(therm.kT, p_temp)

        Create a function (``lognorm``) and use it as the prior of the
        ``nH`` parameter of the ``abs1`` model component::

            >>> def lognorm(x):
            ...     nh = 20
            ...     sigma = 0.5  # use a sigma of 0.5
            ...     # nH is in units of 10^-22 so convert
            ...     dx = np.log10(x) + 22 - nh
            ...     norm = sigma / np.sqrt(2 * np.pi)
            ...     return norm * np.exp(-0.5 * dx * dx / (sigma * sigma))
            ...
            >>> create_model_component('xsphabs', 'abs1')
            >>> set_prior(abs1.nH, lognorm)

        """

        # NOTE: the second piece of code is indented in the example
        #       above because otherwise sphinx seems to think that the
        #       colon at the end of the "def lognorm" line ends the
        #       code block.

        self._pyblocxs.set_prior(par, prior)

    def get_prior(self, par):
        """Return the prior function for a parameter (MCMC).

        The default behavior of the pyBLoCXS MCMC sampler (run by the
        `get_draws` function) is to use a flat prior for each parameter.
        The `get_prior` routine finds the current prior assigned to
        a parameter, and `set_prior` is used to change it.

        Parameters
        ----------
        par : a `sherpa.models.parameter.Parameter` instance
           A parameter of a model instance.

        Returns
        -------
        prior
           The parameter prior set by a previous call to `set_prior`.
           This may be a function or model instance.

        Raises
        ------
        ValueError
           If a prior has not been set for the parameter.

        See Also
        --------
        set_prior : Set the prior function to use with a parameter.

        Examples
        --------

        >>> prior = get_prior(bgnd.c0)
        >>> print(prior)

        """
        return self._pyblocxs.get_prior(par)

    def list_priors(self):
        """Return the priors set for model parameters, if any.

        Returns
        -------
        priors : dict
           The dictionary of mappings between parameters (keys)
           and prior functions (values) created by `set_prior`.

        See Also
        --------
        get_prior : Return the prior function for a parameter (MCMC).
        set_prior : Set the prior function to use with a parameter.

        Examples
        --------

        In this example a prior on the ``PhoIndex`` parameter of the
        ``pl`` instance has been set to be a gaussian:

        >>> list_priors()
        {'pl.PhoIndex': <Gauss1D model instance 'gauss1d.gline'>}

        """
        return self._pyblocxs.list_priors()

    def list_samplers(self):
        """List the MCMC samplers.

        Returns
        -------
        samplers : list of str
           A list of the names (in lower case) that can be used with
           `set_sampler`.

        See Also
        --------
        get_sampler_name : Return the name of the current MCMC sampler.

        Examples
        --------

        >>> list_samplers()
        ['metropolismh', 'fullbayes', 'mh', 'pragbayes']

        """
        return self._pyblocxs.list_samplers()

    # DOC-TODO: add pointers on what to do with the return values
    def get_draws(self,
                  id: Optional[IdType] = None,
                  otherids: Sequence[IdType] = (),
                  niter=1000, covar_matrix=None):
        """Run the pyBLoCXS MCMC algorithm.

        The function runs a Markov Chain Monte Carlo (MCMC) algorithm
        designed to carry out Bayesian Low-Count X-ray Spectral
        (BLoCXS) analysis. It explores the model parameter space at
        the suspected statistic minimum (i.e.  after using `fit`). The return
        values include the statistic value, parameter values, and an
        acceptance flag indicating whether the row represents a jump from the
        current location or not. For more information see the
        `sherpa.sim` module and the reference given below.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : sequence of int or str, optional
           Other data sets to use in the calculation.
        niter : int, optional
           The number of draws to use. The default is ``1000``.
        covar_matrix : 2D array, optional
           The covariance matrix to use. If ``None`` then the
           result from `get_covar_results().extra_output` is used.

        Returns
        -------
        stats, accept, params
           The results of the MCMC chain. The stats and accept arrays
           contain ``niter+1`` elements, with the first row being the
           starting values. The params array has ``(nparams, niter+1)``
           elements, where nparams is the number of free parameters in
           the model expression, and the first column contains the
           values that the chain starts at. The accept array contains
           boolean values, indicating whether the jump, or step, was
           accepted (``True``), so the parameter values and statistic
           change, or it wasn't, in which case there is no change to
           the previous row. The `sherpa.utils.get_error_estimates`
           routine can be used to calculate the credible one-sigma
           interval from the params array.

        See Also
        --------
        covar, fit, get_sampler, plot_cdf, plot_pdf, plot_scatter,
        plot_trace, set_prior, set_rng, set_sampler

        Notes
        -----
        The chain is run using fit information associated with the
        specified data set, or sets, the currently set sampler
        (`set_sampler`) and parameter priors (`set_prior`), for a
        specified number of iterations. The model should have been fit
        to find the best-fit parameters, and `covar` called, before
        running `get_draws`. The results from `get_draws` is used to
        estimate the parameter distributions.

        The `set_rng` routine is used to control how the random
        numbers are generated.

        References
        ----------

        `"Analysis of Energy Spectra with Low Photon Counts via
        Bayesian Posterior Simulation", van Dyk, D.A., Connors, A.,
        Kashyap, V.L., & Siemiginowska, A.  2001, Ap.J., 548, 224
        <https://adsabs.harvard.edu/abs/2001ApJ...548..224V>`_

        Examples
        --------
        Fit a source and then run a chain to investigate the parameter
        distributions.  The distribution of the stats values created
        by the chain is then displayed, using `plot_trace`, and the
        parameter distributions for the first two thawed parameters
        are displayed. The first one as a cumulative distribution
        using `plot_cdf` and the second one as a probability
        distribution using `plot_pdf`. Finally the acceptance fraction
        (number of draws where the chain moved) is displayed.
        Note that in a full analysis session a burn-in period would
        normally be removed from the chain before using the
        results.

        >>> fit()
        >>> covar()
        >>> stats, accept, params = get_draws(1, niter=1e4)
        >>> plot_trace(stats, name='stat')
        >>> names = [p.fullname for p in get_source().pars if not p.frozen]
        >>> plot_cdf(params[0,:], name=names[0], xlabel=names[0])
        >>> plot_pdf(params[1,:], name=names[1], xlabel=names[1])
        >>> accept[:-1].sum() * 1.0 / len(accept - 1)
        0.4287

        The following runs the chain on multiple data sets, with
        identifiers 'core', 'jet1', and 'jet2':

        >>> stats, accept, params = get_draws('core', ['jet1', 'jet2'], niter=1e4)

        """

        ids, fit = self._get_fit(id, otherids)

        # Allow the user to jump from a user defined point in parameter space?
        # Meaning let the user set up parameter space without fitting first.

        # fit_results = self.get_fit_results()
        # if fit_results is None:
        #    raise TypeError("Fit has not been run")

        if covar_matrix is None:
            covar_results = self.get_covar_results()
            if covar_results is None:
                raise TypeError("Covariance has not been calculated")

            covar_matrix = covar_results.extra_output

        return self._pyblocxs.get_draws(fit, covar_matrix,
                                        niter=niter, rng=self.get_rng())

    ###########################################################################
    # Basic plotting
    ###########################################################################
    #
    # Plot object access
    #

    def get_split_plot(self):
        """Return the plot attributes for displays with multiple plots.

        Returns
        -------
        splot : a `sherpa.plot.SplitPlot` instance

        See Also
        --------
        contour, plot

        Examples
        --------

        Change the layout of the plot and contour commands to display
        three vertical plots:

        >>> sp = get_split_plot()
        >>> sp.rows = 3
        >>> sp.cols = 1

        """
        return self._splitplot

    def get_data_plot(self,
                      id: Optional[IdType] = None,
                      recalc=True):
        """Return the data used by plot_data.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_data` (or `get_data_plot`) are returned, otherwise
           the data is re-generated.

        Returns
        -------
        data : a `sherpa.plot.DataPlot` instance
           An object representing the data used to create the plot by
           `plot_data`. The relationship between the returned values
           and the values in the data set depend on the data type. For
           example PHA data are plotted in units controlled by
           `sherpa.astro.ui.set_analysis`, but are stored as channels
           and counts, and may have been grouped and the background
           estimate removed.

        See Also
        --------
        get_data_plot_prefs : Return the preferences for plot_data.
        get_default_id : Return the default data set identifier.
        plot_data : Plot the data values.

        """

        # Allow an answer to be returned if recalc is False and no
        # data has been loaded. However, it's not obvious what the
        # answer should be if recalc=False and the dataset has
        # changed type since get_data_plot was last called.
        #
        if recalc:
            data = self.get_data(id)
        else:
            data = self._get_data(id)

        # This uses the implicit conversion of bool to 0 or 1.
        #
        idx = isinstance(data, sherpa.data.Data1DInt)
        plotobj = self._plot_types["data"][idx]
        if recalc:
            plotobj.prepare(data, self.get_stat())

        return plotobj

    # DOC-TODO: discussion of preferences needs better handling
    # of how it interacts with the chosen plot backend.
    #
    def get_plot_prefs(self,
                       plottype: str,
                       id: Optional[IdType] = None,
                       **kwargs):
        """Return the preferences for the given plot type.

        .. versionadded:: 4.16.0

        Parameters
        ----------
        plottype : str
           The type of plt, such as "data", "model", or "resid". The
           "fit" argument is not supported.
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        prefs : dict
           Changing the values of this dictionary will change any new
           data plots. This dictionary will be empty if no plot
           backend is available.

        See Also
        --------
        get_data_plot_prefs, get_model_plot_prefs,
        set_xlinear, set_xlog, set_ylinear, set_ylog

        Notes
        -----
        The meaning of the fields depend on the chosen plot backend.
        A value of ``None`` means to use the default value for that
        attribute, or not to use that setting.

        The "fit" argument can not be used, even though there is a
        get_fit_plot call. Either use the "data" or "model" arguments
        to access the desired plot type, or use get_fit_plot() and
        access the dataplot and modelplot attributes directly.

        Examples
        --------

        After these commands, any data plot will use a green symbol
        and not display Y error bars.

        >>> prefs = get_plot_prefs("data")
        >>> prefs['color'] = 'green'
        >>> prefs['yerrorbars'] = False

        """

        # By accepting kwargs we can send in the bkg_id argument from
        # the astro version without needing to override this call.
        #
        # This relies on each key in _plot_types having a
        # corresponding get_<key>_plot() call.
        #
        ptype = self._get_plottype(plottype)

        # To also catch the astro layer we do not check exactly, just
        # that fit is an argument. Technically we could adjust the
        # error message to include 'bkg'/'bkg_model' but this is
        # excessive.
        #
        if ptype.find("fit") > -1:
            raise ArgumentErr(f"Use 'data' or 'model' instead of '{plottype}'")

        get = getattr(self, f"get_{ptype}_plot")
        plotobj = get(id, recalc=False, **kwargs)
        return get_plot_prefs(plotobj)

    def get_data_plot_prefs(self,
                            id: Optional[IdType] = None):
        """Return the preferences for plot_data.

        The plot preferences may depend on the data set,
        so it is now an optional argument.

        .. versionchanged:: 4.12.2
           The id argument has been given.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        prefs : dict
           Changing the values of this dictionary will change any new
           data plots. This dictionary will be empty if no plot
           backend is available.

        See Also
        --------
        get_plot_prefs, plot_data, set_xlinear, set_xlog, set_ylinear,
        set_ylog

        Notes
        -----
        The meaning of the fields depend on the chosen plot backend.
        A value of ``None`` means to use the default value for that
        attribute, unless indicated otherwise. These preferences
        are used by the following commands: `plot_data`, `plot_bkg`,
        `plot_ratio`, and the "fit" variants, such as `plot_fit`,
        `plot_fit_resid`, and `plot_bkg_fit`.

        The following preferences are recognized by the matplotlib
        backend:

        ``alpha``
           The transparency value used to draw the line or symbol,
           where 0 is fully transparent and 1 is fully opaque.

        ``barsabove``
           The barsabove argument for the matplotlib errorbar function.

        ``capsize``
           The capsize argument for the matplotlib errorbar function.

        ``color``
           The color to use (will be over-ridden by more-specific
           options below). The default is ``None``.

        ``ecolor``
           The color to draw error bars. The default is ``None``.

        ``linestyle``
           How should the line connecting the data points be drawn.
           The default is 'None', which means no line is drawn.

        ``marker``
           What style is used for the symbols. The default is ``'.'``
           which indicates a point.

        ``markerfacecolor``
           What color to draw the symbol representing the data points.
           The default is ``None``.

        ``markersize``
           What size is the symbol drawn. The default is ``None``.

        ``xerrorbars``
           Should error bars be drawn for the X axis. The default is
           ``False``.

        ``xlog``
           Should the X axis be drawn with a logarithmic scale? The
           default is ``False``. This field can also be changed with the
           `set_xlog` and `set_xlinear` functions.

        ``yerrorbars``
           Should error bars be drawn for the Y axis. The default is
           ``True``.

        ``ylog``
           Should the Y axis be drawn with a logarithmic scale? The
           default is ``False``. This field can also be changed with the
           `set_ylog` and `set_ylinear` functions.

        Examples
        --------

        After these commands, any data plot will use a green symbol
        and not display Y error bars.

        >>> prefs = get_data_plot_prefs()
        >>> prefs['color'] = 'green'
        >>> prefs['yerrorbars'] = False

        """

        plotobj = self.get_data_plot(id, recalc=False)
        return get_plot_prefs(plotobj)

    # also in sherpa.astro.utils (copies this docstring)
    def get_model_plot(self,
                       id: Optional[IdType] = None,
                       recalc=True):
        """Return the data used to create the model plot.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_model` (or `get_model_plot`) are returned, otherwise
           the data is re-generated.

        Returns
        -------
        instance
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
        >>> print(mplot)

        """

        if recalc:
            data = self.get_data(id)
        else:
            data = self._get_data(id)

        idx = isinstance(data, sherpa.data.Data1DInt)
        plotobj = self._plot_types["model"][idx]
        if recalc:
            plotobj.prepare(data, self.get_model(id), self.get_stat())

        return plotobj

    # also in sherpa.astro.utils (does not copy this docstring)
    def get_source_plot(self,
                        id: Optional[IdType] = None,
                        recalc=True):
        """Return the data used to create the source plot.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_source` (or `get_source_plot`) are returned, otherwise
           the data is re-generated.

        Returns
        -------
        instance
           An object representing the data used to create the plot by
           `plot_source`. The return value depends on the data
           set (e.g. 1D binned or un-binned).

        See Also
        --------
        get_model_plot : Return the data used to create the model plot.
        plot_model : Plot the model for a data set.
        plot_source : Plot the source expression for a data set.

        Examples
        --------

        Retrieve the source plot information for the default data
        set and then display it:

        >>> splot = get_source_plot()
        >>> print(splot)

        Return the plot data for data set 2, and then use it to create
        a plot:

        >>> s2 = get_source_plot(2)
        >>> s2.plot()

        Display the two source plots for the 'jet' and 'core' datasets
        on the same plot:

        >>> splot1 = get_source_plot(id='jet')
        >>> splot2 = get_source_plot(id='core')
        >>> splot1.plot()
        >>> splot2.overplot()

        """

        idval = self._fix_id(id)
        mdl = self._models.get(idval, None)
        if mdl is not None:
            raise IdentifierErr(f"Convolved model\n'{mdl.name}'\n"
                                f" is set for dataset {idval}. "
                                "You should use get_model_plot instead.")

        if recalc:
            data = self.get_data(idval)
        else:
            data = self._get_data(idval)

        idx = isinstance(data, sherpa.data.Data1DInt)
        plotobj = self._plot_types["source"][idx]
        if recalc:
            plotobj.prepare(data, self.get_source(idval), self.get_stat())

        return plotobj

    def get_model_component_plot(self, id, model=None, recalc=True):
        """Return the data used to create the model-component plot.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to use (the name, if a string).
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_model_component` (or `get_model_component_plot`)
           are returned, otherwise the data is re-generated.

        Returns
        -------
        instance
           An object representing the data used to create the plot by
           `plot_model_component`. The return value depends on the
           data set (e.g. 1D binned or un-binned).

        See Also
        --------
        get_model_plot : Return the data used to create the model plot.
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

        Return the plot data for the ``pl`` component used in the
        default data set:

        >>> cplot = get_model_component_plot(pl)

        Return the full source model (``fplot``) and then for the
        components ``gal * pl`` and ``gal * gline``, for the data set
        'jet':

        >>> fmodel = xsphabs.gal * (powlaw1d.pl + gauss1d.gline)
        >>> set_source('jet', fmodel)
        >>> fit('jet')
        >>> fplot = get_model_plot('jet')
        >>> plot1 = get_model_component_plot('jet', pl*gal)
        >>> plot2 = get_model_component_plot('jet', gline*gal)

        """
        if model is None:
            id, model = model, id
        model = self._check_model(model)

        if recalc:
            data = self.get_data(id)
        else:
            data = self._get_data(id)

        idx = isinstance(data, sherpa.data.Data1DInt)
        plotobj = self._plot_types["model_component"][idx]
        if recalc:
            plotobj.prepare(data, model, self.get_stat())

        return plotobj

    def get_model_components_plot(self,
                                  id: Optional[IdType] = None
                                  ) -> MultiPlot:
        """Return the data used by plot_model_components.

        .. versionadded:: 4.16.1

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        plot : MultiPlot
           A plot object containing the individual plot objects.

        See Also
        --------
        get_model_plot, plot_model, plot_model_component,
        plot_model_components

        Notes
        -----

        Unlike get_model_component this routine does not accept
        either model or recalc arguments.

        Examples
        --------

        Return the plot data for the individual components used in the
        default data set. In this case it will report plots for
        ``gal * pl`` and ``gal * gline``:

        >>> set_source(xsphabs.gal * (powlaw1d.pl + gauss1d.gline))
        >>> cplots = get_model_components_plot()

        """

        idval = self._fix_id(id)
        model = self.get_source(idval)
        return get_components_helper(self.get_model_component_plot,
                                     model=model, idval=idval)

    # sherpa.astro.utils version copies this docstring
    def get_source_component_plot(self, id, model=None, recalc=True):
        """Return the data used by plot_source_component.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to use (the name, if a string).
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_source_component` (or `get_source_component_plot`)
           are returned, otherwise the data is re-generated.

        Returns
        -------
        instance
           An object representing the data used to create the plot by
           `plot_source_component`. The return value depends on the
           data set (e.g. 1D binned or un-binned).

        See Also
        --------
        get_source_plot : Return the data used to create the source plot.
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

        Return the plot data for the ``pl`` component used in the
        default data set:

        >>> cplot = get_source_component_plot(pl)

        Return the full source model (``fplot``) and then for the
        components ``gal * pl`` and ``gal * gline``, for the data set
        'jet':

        >>> fmodel = xsphabs.gal * (powlaw1d.pl + gauss1d.gline)
        >>> set_source('jet', fmodel)
        >>> fit('jet')
        >>> fplot = get_source('jet')
        >>> plot1 = get_source_component_plot('jet', pl*gal)
        >>> plot2 = get_source_component_plot('jet', gline*gal)

        """
        if model is None:
            id, model = model, id
        model = self._check_model(model)

        if recalc:
            data = self.get_data(id)
        else:
            data = self._get_data(id)

        if isinstance(model, sherpa.models.TemplateModel):
            plotobj = self._comptmplsrcplot
        else:
            idx = isinstance(data, sherpa.data.Data1DInt)
            plotobj = self._plot_types["source_component"][idx]

        if recalc:
            plotobj.prepare(data, model, self.get_stat())

        return plotobj

    def get_source_components_plot(self,
                                   id: Optional[IdType] = None
                                   ) -> MultiPlot:
        """Return the data used by plot_source_components.

        .. versionadded:: 4.16.1

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        plot : MultiPlot
           A plot object containing the individual plot objects.

        See Also
        --------
        get_source_plot, plot_source, plot_source_component,
        plot_source_components

        Notes
        -----

        Unlike get_source_component this routine does not accept
        either model or recalc arguments.

        Examples
        --------

        Return the plot data for the individual components used in the
        default data set. In this case it will report plots for
        ``gal * pl`` and ``gal * gline``:

        >>> set_source(xsphabs.gal * (powlaw1d.pl + gauss1d.gline))
        >>> cplots = get_source_components_plot()

        """

        idval = self._fix_id(id)
        model = self.get_source(idval)
        return get_components_helper(self.get_source_component_plot,
                                     model=model, idval=idval)

    def get_model_plot_prefs(self, id: Optional[IdType] = None):
        """Return the preferences for plot_model.

        The plot preferences may depend on the data set,
        so it is now an optional argument.

        .. versionchanged:: 4.12.2
           The id argument has been given.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        prefs : dict
           Changing the values of this dictionary will change any new
           model plots. This dictionary will be empty if no plot
           backend is available.

        See Also
        --------
        get_plot_prefs, plot_model, set_xlinearm set_xlog, set_ylinear,
        set_ylog

        Notes
        -----
        The meaning of the fields depend on the chosen plot backend.
        A value of ``None`` means to use the default value for that
        attribute, unless indicated otherwise. These preferences are
        used by the following commands: `plot_model`, `plot_ratio`,
        `plot_bkg_model`, and the "fit" variants, such as `plot_fit`,
        `plot_fit_resid`, and `plot_bkg_fit`.

        The preferences recognized by the matplotlib backend are the
        same as for `get_data_plot_prefs`.

        Examples
        --------

        After these commands, any model plot will use a green line
        to display the model:

        >>> prefs = get_model_plot_prefs()
        >>> prefs['color'] = 'green'

        """

        plotobj = self.get_model_plot(id, recalc=False)
        return get_plot_prefs(plotobj)

    def get_fit_plot(self,
                     id: Optional[IdType] = None,
                     recalc=True):
        """Return the data used to create the fit plot.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_fit` (or `get_fit_plot`) are returned, otherwise
           the data is re-generated.

        Returns
        -------
        data : a `sherpa.plot.FitPlot` instance
           An object representing the data used to create the plot by
           `plot_fit`. It contains the data from `get_data_plot`
           and `get_model_plot` in the ``dataplot`` and ``modelplot``
           attributes.

        See Also
        --------
        get_data_plot_prefs : Return the preferences for plot_data.
        get_model_plot_prefs : Return the preferences for plot_model.
        get_default_id : Return the default data set identifier.
        plot_data : Plot the data values.
        plot_model : Plot the model for a data set.

        Examples
        --------

        Create the data needed to create the "fit plot" for the default
        data set and display it:

        >>> fplot = get_fit_plot()
        >>> print(fplot)

        Return the plot data for data set 2, and then use it to create
        a plot:

        >>> f2 = get_fit_plot(2)
        >>> f2.plot()

        The fit plot consists of a combination of a data plot and a
        model plot, which are captured in the `dataplot` and `modelplot`
        attributes of the return value. These can be used to display
        the plots individually, such as:

        >>> f2.dataplot.plot()
        >>> f2.modelplot.plot()

        or, to combine the two:

        >>> f2.dataplot.plot()
        >>> f2.modelplot.overplot()

        """

        plotobj = self._plot_types["fit"][0]
        if recalc:
            dataobj = self.get_data_plot(id, recalc=recalc)
            modelobj = self.get_model_plot(id, recalc=recalc)
            plotobj.prepare(dataobj, modelobj)

        return plotobj

    def get_resid_plot(self,
                       id: Optional[IdType] = None,
                       recalc=True):
        """Return the data used by plot_resid.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_resid` (or `get_resid_plot`) are returned, otherwise
           the data is re-generated.

        Returns
        -------
        resid_data : a `sherpa.plot.ResidPlot` instance

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

        Display the contents of the residuals plot for data set 2:

        >>> print(get_resid_plot(2))

        Overplot the residuals plot from the 'core' data set on the 'jet'
        data set:

        >>> r1 = get_resid_plot('jet')
        >>> r2 = get_resid_plot('core')
        >>> r1.plot()
        >>> r2.overplot()

        """

        # Allow an answer to be returned if recalc is False and no
        # data has been loaded. However, it's not obvious what the
        # answer should be if recalc=False and the dataset has
        # changed type since get_delchi_plot was last called.
        #
        if recalc:
            data = self.get_data(id)
        else:
            data = self._get_data(id)

        # This uses the implicit conversion of bool to 0 or 1.
        #
        idx = isinstance(data, sherpa.data.Data1DInt)
        plotobj = self._plot_types["resid"][idx]
        if recalc:
            plotobj.prepare(data, self.get_model(id), self.get_stat())

        return plotobj

    def get_delchi_plot(self,
                        id: Optional[IdType] = None,
                        recalc=True):
        """Return the data used by plot_delchi.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_delchi` (or `get_delchi_plot`) are returned, otherwise
           the data is re-generated.

        Returns
        -------
        resid_data : a `sherpa.plot.DelchiPlot` instance

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

        Display the contents of the residuals plot for data set 2:

        >>> print(get_delchi_plot(2))

        Overplot the residuals plot from the 'core' data set on the 'jet'
        data set:

        >>> r1 = get_delchi_plot('jet')
        >>> r2 = get_delchi_plot('core')
        >>> r1.plot()
        >>> r2.overplot()

        """

        # Allow an answer to be returned if recalc is False and no
        # data has been loaded. However, it's not obvious what the
        # answer should be if recalc=False and the dataset has
        # changed type since get_delchi_plot was last called.
        #
        if recalc:
            data = self.get_data(id)
        else:
            data = self._get_data(id)

        # This uses the implicit conversion of bool to 0 or 1.
        #
        idx = isinstance(data, sherpa.data.Data1DInt)
        plotobj = self._plot_types["delchi"][idx]
        if recalc:
            plotobj.prepare(data, self.get_model(id), self.get_stat())

        return plotobj

    def get_chisqr_plot(self,
                        id: Optional[IdType] = None,
                        recalc=True):
        """Return the data used by plot_chisqr.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_chisqr` (or `get_chisqr_plot`) are returned, otherwise
           the data is re-generated.

        Returns
        -------
        resid_data : a `sherpa.plot.ChisqrPlot` instance

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

        Display the contents of the residuals plot for data set 2:

        >>> print(get_chisqr_plot(2))

        Overplot the residuals plot from the 'core' data set on the 'jet'
        data set:

        >>> r1 = get_chisqr_plot('jet')
        >>> r2 = get_chisqr_plot('core')
        >>> r1.plot()
        >>> r2.overplot()

        """

        # Allow an answer to be returned if recalc is False and no
        # data has been loaded. However, it's not obvious what the
        # answer should be if recalc=False and the dataset has
        # changed type since get_delchi_plot was last called.
        #
        if recalc:
            data = self.get_data(id)
        else:
            data = self._get_data(id)

        # This uses the implicit conversion of bool to 0 or 1.
        #
        idx = isinstance(data, sherpa.data.Data1DInt)
        plotobj = self._plot_types["chisqr"][idx]
        if recalc:
            plotobj.prepare(data, self.get_model(id), self.get_stat())

        return plotobj

    def get_ratio_plot(self,
                       id: Optional[IdType] = None,
                       recalc=True):
        """Return the data used by plot_ratio.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_ratio` (or `get_ratio_plot`) are returned, otherwise
           the data is re-generated.

        Returns
        -------
        ratio_data : a `sherpa.plot.RatioPlot` instance

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

        Display the contents of the ratio plot for data set 2:

        >>> print(get_ratio_plot(2))

        Overplot the ratio plot from the 'core' data set on the 'jet'
        data set:

        >>> r1 = get_ratio_plot('jet')
        >>> r2 = get_ratio_plot('core')
        >>> r1.plot()
        >>> r2.overplot()

        """

        # Allow an answer to be returned if recalc is False and no
        # data has been loaded. However, it's not obvious what the
        # answer should be if recalc=False and the dataset has
        # changed type since get_delchi_plot was last called.
        #
        if recalc:
            data = self.get_data(id)
        else:
            data = self._get_data(id)

        # This uses the implicit conversion of bool to 0 or 1.
        #
        idx = isinstance(data, sherpa.data.Data1DInt)
        plotobj = self._plot_types["ratio"][idx]
        if recalc:
            plotobj.prepare(data, self.get_model(id), self.get_stat())

        return plotobj

    def get_data_contour(self,
                         id: Optional[IdType] = None,
                         recalc=True):
        """Return the data used by contour_data.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `contour_data` (or `get_data_contour`) are returned,
           otherwise the data is re-generated.

        Returns
        -------
        resid_data : a `sherpa.plot.DataContour` instance
           The ``y`` attribute contains the residual values and the ``x0``
           and ``x1`` arrays contain the corresponding coordinate values, as
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

        plotobj = self._contour_types["data"]
        if recalc:
            plotobj.prepare(self.get_data(id), self.get_stat())

        return plotobj

    def get_contour_prefs(self,
                          contourtype: str,
                          id: Optional[IdType] = None):
        """Return the preferences for the given contour type.

        .. versionadded:: 4.16.0

        Parameters
        ----------
        contourtype : str
           The contour type, such as "data", "model", or "resid". The
           "fit" argument is not supported.
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.

        Returns
        -------
        prefs : dict
           Changing the values of this dictionary will change any new
           contour plots. This dictionary will be empty if no plot
           backend is available.

        See Also
        --------
        get_data_contour_prefs, get_model_contour_prefs

        Notes
        -----
        The meaning of the fields depend on the chosen plot backend.
        A value of ``None`` means to use the default value for that
        attribute, or not to use that setting.

        The "fit" argument can not be used, even though there is a
        get_fit_contour call. Either use the "data" or "model" arguments
        to access the desired plot type, or use get_fit_contour() and
        access the datacontour and modelcontour attributes directly.

        Examples
        --------

        Change the contours of the model to be drawn partly opaque (with
        the matplotlib backend):

        >>> prefs = get_contour_prefs("data")
        >>> prefs['alpha'] = 0.7
        >>> contour_data()
        >>> contour_model(overcontour=True)

        """

        # This is an identity check but will throw an error if the
        # argument is invalid. It requires that each contour type has
        # a get_<contourtype>_contour() call.
        #
        ctype = self._get_contourtype(contourtype)
        if ctype == "fit":
            raise ArgumentErr("Use 'data' or 'model' instead of 'fit'")

        get = getattr(self, f"get_{ctype}_contour")
        return get(id, recalc=False).contour_prefs

    def get_data_contour_prefs(self):
        """Return the preferences for contour_data.

        Returns
        -------
        prefs : dict
           Changing the values of this dictionary will change any new
           contour plots. The default is an empty dictionary.

        See Also
        --------
        contour_data, get_contour_prefs

        Notes
        -----
        The meaning of the fields depend on the chosen plot backend.
        A value of ``None`` (or not set) means to use the default value
        for that attribute, unless indicated otherwise.

        ``alpha``
           The transparency value used to draw the contours,
           where 0 is fully transparent and 1 is fully opaque.

        ``colors``
           The colors to draw the contours.

        ``linewidths``
           What thickness of line to draw the contours.

        ``xlog``
           Should the X axis be drawn with a logarithmic scale? The
           default is ``False``.

        ``ylog``
           Should the Y axis be drawn with a logarithmic scale? The
           default is ``False``.

        Examples
        --------

        Change the contours to be drawn in 'green':

        >>> contour_data()
        >>> prefs = get_data_contour_prefs()
        >>> prefs['color'] = 'green'
        >>> contour_data()

        """
        return self.get_data_contour(id, recalc=False).contour_prefs

    def get_model_contour(self,
                          id: Optional[IdType] = None,
                          recalc=True):
        """Return the data used by contour_model.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `contour_model` (or `get_model_contour`) are returned,
           otherwise the data is re-generated.

        Returns
        -------
        model_data : a `sherpa.plot.ModelContour` instance
           The ``y`` attribute contains the model values and the ``x0``
           and ``x1`` arrays contain the corresponding coordinate values, as
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

        plotobj = self._contour_types["model"]
        if recalc:
            plotobj.prepare(self.get_data(id), self.get_model(id), self.get_stat())

        return plotobj

    def get_source_contour(self,
                           id: Optional[IdType] = None,
                           recalc=True):
        """Return the data used by contour_source.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `contour_source` (or `get_source_contour`) are returned,
           otherwise the data is re-generated.

        Returns
        -------
        source_data : a `sherpa.plot.SourceContour` instance
           The ``y`` attribute contains the model values and the ``x0``
           and ``x1`` arrays contain the corresponding coordinate values, as
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

        plotobj = self._contour_types["source"]
        if recalc:
            plotobj.prepare(self.get_data(id), self.get_source(id), self.get_stat())

        return plotobj

    def get_model_contour_prefs(self):
        """Return the preferences for contour_model.

        Returns
        -------
        prefs : dict
           Changing the values of this dictionary will change any new
           contour plots.

        See Also
        --------
        contour_model, get_contour_prefs

        Notes
        -----
        The meaning of the fields depend on the chosen plot backend.
        A value of ``None`` (or not set) means to use the default value
        for that attribute, unless indicated otherwise.

        ``alpha``
           The transparency value used to draw the contours,
           where 0 is fully transparent and 1 is fully opaque.

        ``colors``
           The colors to draw the contours.

        ``linewidths``
           What thickness of line to draw the contours.

        ``xlog``
           Should the X axis be drawn with a logarithmic scale? The
           default is ``False``.

        ``ylog``
           Should the Y axis be drawn with a logarithmic scale? The
           default is ``False``.

        Examples
        --------

        Change the contours for the model to be drawn in 'orange':

        >>> prefs = get_model_contour_prefs()
        >>> prefs['color'] = 'orange'
        >>> contour_data()
        >>> contour_model(overcontour=True)

        """
        return self.get_model_contour(id, recalc=False).contour_prefs

    def get_fit_contour(self,
                        id: Optional[IdType] = None,
                        recalc=True):
        """Return the data used by contour_fit.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `contour_fit` (or `get_fit_contour`) are returned,
           otherwise the data is re-generated.

        Returns
        -------
        fit_data : a `sherpa.plot.FitContour` instance
           An object representing the data used to create the plot by
           `contour_fit`. It contains the data from `get_data_contour`
           and `get_model_contour` in the ``datacontour`` and
           ``modelcontour`` attributes.

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
        get_model_image : Return the data used by image_model.
        contour_data : Contour the values of an image data set.
        contour_model : Contour the values of the model, including any PSF.
        image_data : Display a data set in the image viewer.
        image_model : Display the model for a data set in the image viewer.

        Examples
        --------

        Return the contour data for the default data set:

        >>> finfo = get_fit_contour()

        """

        plotobj = self._contour_types["fit"]
        if recalc:
            dataobj = self.get_data_contour(id, recalc=recalc)
            modelobj = self.get_model_contour(id, recalc=recalc)
            plotobj.prepare(dataobj, modelobj)

        return plotobj

    def get_resid_contour(self,
                          id: Optional[IdType] = None,
                          recalc=True):
        """Return the data used by contour_resid.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `contour_resid` (or `get_resid_contour`) are returned,
           otherwise the data is re-generated.

        Returns
        -------
        resid_data : a `sherpa.plot.ResidContour` instance
           The ``y`` attribute contains the residual values and the ``x0``
           and ``x1`` arrays contain the corresponding coordinate values, as
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

        plotobj = self._contour_types["resid"]
        if recalc:
            plotobj.prepare(self.get_data(id), self.get_model(id), self.get_stat())

        return plotobj

    def get_ratio_contour(self,
                          id: Optional[IdType] = None,
                          recalc=True):
        """Return the data used by contour_ratio.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `contour_ratio` (or `get_ratio_contour`) are returned,
           otherwise the data is re-generated.

        Returns
        -------
        ratio_data : a `sherpa.plot.RatioContour` instance
           The ``y`` attribute contains the ratio values and the ``x0``
           and ``x1`` arrays contain the corresponding coordinate values, as
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

        plotobj = self._contour_types["ratio"]
        if recalc:
            plotobj.prepare(self.get_data(id), self.get_model(id), self.get_stat())

        return plotobj

    def get_psf_contour(self,
                        id: Optional[IdType] = None,
                        recalc=True):
        """Return the data used by contour_psf.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `contour_psf` (or `get_psf_contour`) are returned,
           otherwise the data is re-generated.

        Returns
        -------
        psf_data : a `sherpa.plot.PSFContour` instance

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

        plotobj = self._contour_types["psf"]
        if recalc:
            plotobj.prepare(self.get_psf(id), self.get_data(id))

        return plotobj

    def get_kernel_contour(self,
                           id: Optional[IdType] = None,
                           recalc=True):
        """Return the data used by contour_kernel.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `contour_kernel` (or `get_kernel_contour`) are returned,
           otherwise the data is re-generated.

        Returns
        -------
        psf_data : a `sherpa.plot.PSFKernelContour` instance

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

        plotobj = self._contour_types["kernel"]
        if recalc:
            plotobj.prepare(self.get_psf(id), self.get_data(id))

        return plotobj

    def get_psf_plot(self,
                     id: Optional[IdType] = None,
                     recalc=True):
        """Return the data used by plot_psf.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_psf` (or `get_psf_plot`) are returned, otherwise
           the data is re-generated.

        Returns
        -------
        psf_plot : a `sherpa.plot.PSFPlot` instance

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

        plotobj = self._plot_types["psf"][0]
        if recalc:
            plotobj.prepare(self.get_psf(id), self.get_data(id))

        return plotobj

    def get_kernel_plot(self,
                        id: Optional[IdType] = None,
                        recalc=True):
        """Return the data used by plot_kernel.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        recalc : bool, optional
           If ``False`` then the results from the last call to
           `plot_kernel` (or `get_kernel_plot`) are returned, otherwise
           the data is re-generated.

        Returns
        -------
        kernel_plot : a `sherpa.plot.PSFKernelPlot` instance

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

        plotobj = self._plot_types["kernel"][0]
        if recalc:
            plotobj.prepare(self.get_psf(id), self.get_data(id))

        return plotobj

    #
    # Line plots
    #
    def _multi_plot(self,
                    args: Sequence[IdType],
                    plotmeth: Literal["plot", "contour"] = "plot",
                    rows: Optional[int] = None,
                    cols: Optional[int] = None,
                    **kwargs) -> None:
        """Handle the plot() or contour() call.

        The arguments are split up into groups - a "command" followed
        by optional arguments - and then each group is used to create
        a plot/contour.

        Parameters
        ----------
        args
            The arguments to the call. It can not be empty.
        plotmeth
            The call.
        rows, cols
            The number of rows or columns (to override split plot).
        kwargs
            The keyword arguments to apply to each plot or contour.

        Notes
        -----
        Each group consists of the name of the plot/contour, followed
        by optional arguments. So

            self._multi_plot(["data", "data", "org"], plotmeth="plot")

        has two "data" style plots, the first with the default
        identifier and the second with the identifier "org". It is the
        need to scan through the arguments that leads to the ban on
        identifiers matching the plot/contour "command names".

        Each plot/contour command matches a key in the
        _plot/contour_type_names dictionary, and the value is used to
        get the underlying "plot object" via get_<value>_plot or
        get_<value>_contour.

        Multiple arguments can be sent, but they rely on positional
        ordering only, so

            args = ["model_component", mdl, "model_component", 2, mdl]
            plotmeth = "plot"

        only works if

            get_model_component_plot(mdl)
            get_model_component_plot(2, mdl)

        are meaningful calls.

        The number of rows and columns depends on the rows and cols
        arguments:

        - if both are set then they are used
        - if either are set then the other is calculated
        - if both are unset then we use the split_plot preferences

        There is then a check, so that they get recalculated if there
        are too many plots. In this case an approach is used to

        - try to go with a square display
        - fill in columns first

        The handling of overplot/contour has been improved - see issue
        #2128 - but this involves using internal features of the
        SplitPlot class. This could be improved, once we have more
        experience of how this all works.

        """

        # This is an internal routine so it is not expected to be
        # called incorrectly, but make sure we catch any such problem.
        #
        if len(args) == 0:
            raise ArgumentTypeErr("plotargs")

        if plotmeth not in ["plot", "contour"]:
            raise ArgumentErr(f"Unsupported plotmeth={plotmeth}")

        if rows is not None and (rows < 1 or not _is_integer(rows)):
            raise ArgumentErr(f"rows must be a positive integer, not {rows}")

        if cols is not None and (cols < 1 or not _is_integer(cols)):
            raise ArgumentErr(f"cols must be a positive integer, not {cols}")

        get = getattr(self, f"_get_{plotmeth}type")
        check = getattr(self, f"_check_{plotmeth}type")

        plots = []
        largs = list(args)
        while largs:
            plottype = largs.pop(0)
            _check_str_type(plottype, "plottype")
            plottype = get(plottype.lower())

            # Collect the arguments for the get_<>_plot/contour
            # call. Loop through until we hit a supported
            # plot/contour type.
            #
            getargs = []
            while largs:
                if check(largs[0]):
                    break

                getargs.append(largs.pop(0))

            funcname = f"get_{plottype}_{plotmeth}"
            getfunc = getattr(self, funcname)

            # Need to make sure we have a copy of each plot
            # object to support plots like
            #    plot("data", 1, "data", 2)
            #
            plots.append(copy.deepcopy(getfunc(*getargs)))

        nplots = len(plots)

        # What is the overplot/contour handling? To address #2128 we
        # do not need to worry about the clearwindow option.  Note
        # that the value (assumed to be a scalar) is not removed from
        # kwargs.
        #
        over = kwargs.get(f"over{plotmeth}", False)
        # clearwindow = kwargs.get("clearwindow", True)

        # Deconstruct the keyword arguments.
        #
        kwstore = get_per_plot_kwargs(nplots, kwargs)

        # Store the original values in case they are overwritten.
        #
        sp = self._splitplot
        nrows_orig = sp.rows
        ncols_orig = sp.cols

        plotfunc = getattr(sp, f"add{plotmeth}")

        # We only care about the size when creating the display, that
        # is, if overplot/contour is set then we assume the size is
        # known. There is a check to make sure the user isn't trying
        # to overplot too many plots, but we allow less plots to be
        # overplot than the original.
        #
        # Should this care about clearwindow too?
        #
        if over:
            # Dig inside the SplitPlot class to find this out.
            #
            nrows, ncols = sp._used.shape
            nplots_used = sp._used.sum()
            if nplots_used < nplots:
                raise ArgumentErr(f"Too many elements for over{plotmeth}")

        elif rows is None and cols is None:
            # Special case the single-plot call
            if nplots == 1:
                nrows = 1
                ncols = 1
            else:
                nrows, ncols = calc_multiplot_size(nrows_orig,
                                                   ncols_orig, nplots)

        else:
            nrows, ncols = calc_multiplot_size(rows, cols, nplots)

        plotfunc = getattr(sp, f"add{plotmeth}")
        sp.reset(rows=nrows, cols=ncols)

        # This is a hack as the SplitPlot API does not provide a way
        # to do this. We want to avoid calling the _clear_window
        # method when calling addplot/contour.
        #
        if over:
            sp._cleared_window = True

        with sherpa.plot.backend:
            for plot, store in zip(plots, kwstore):
                plotfunc(plot, **store)

        # Restore the settings. Is this needed?
        #
        sp.rows = nrows_orig
        sp.cols = ncols_orig

    def _plot(self, plotobj, **kwargs) -> None:
        """Display a plot object

        Parameters
        ----------
        plotobj
            The plot to display
        **kwargs
            The arguments to the plot method of plotobj.

        See Also
        --------
        _contour

        """

        with sherpa.plot.backend:
            plotobj.plot(**kwargs)

    def _set_plot_item(self, plottype: str,
                       item, value) -> None:
        """Change a plot setting.

        This will change all related plot classes; that is setting
        plottype to "data" will change all the associated plot classes
        for "data".

        Parameters
        ----------
        plottype : str
            A valid plot (currently taken from _plot_types) or "all".
        item : str
            The plot setting to change. If there is no such item then
            nothing happens.
        value
            The value to set the item to. There is no validation of
            this value.

        Raises
        ------
        sherpa.utils.err.PlotErr
           An invalid plottype value.

        """

        _check_str_type(plottype, "plottype")
        allowed = list(self._plot_types)

        plottype = plottype.strip().lower()
        if plottype == "all":
            keys = allowed
        else:
            keys = [self._get_plottype(plottype)]

        for key in keys:
            for plot in self._plot_types[key]:
                try:
                    get_plot_prefs(plot)[item] = value
                except AttributeError:
                    # skip it if no preference setting
                    pass

    def set_xlog(self, plottype: str = "all") -> None:
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

        Use a logarithmic scale for the X axis of 'data' plots:

        >>> set_xlog('data')
        >>> plot('data', 'arf')

        All plots use a logarithmic scale for the X axis.

        >>> set_xlog()
        >>> plot_fit()

        """
        self._set_plot_item(plottype, 'xlog', True)

    def set_ylog(self, plottype: str = "all") -> None:
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

        Use a logarithmic scale for the Y axis of 'data' plots:

        >>> set_ylog('data')
        >>> plot('data', 'arf')

        All plots use a logarithmic scale for the Y axis.

        >>> set_ylog()
        >>> plot_fit()

        """
        self._set_plot_item(plottype, 'ylog', True)

    def set_xlinear(self, plottype: str = "all") -> None:
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

    def set_ylinear(self, plottype: str = "all") -> None:
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

    # This documentation is really for the sherpa.astro.ui version, as
    # some of this is not relevant for the sherpa.ui code, but leave
    # like this.
    #
    def plot(self,
             *args,
             # At present export_method does not support this
             # rows: Optional[int] = None,
             # cols: Optional[int] = None,
             **kwargs) -> None:
        """Create one or more plot types.

        The plot function creates one or more plots, depending on the
        arguments it is sent: a plot type, followed by optional
        identifiers, and this can be repeated. If no data set
        identifier is given for a plot type, the default identifier -
        as returned by `get_default_id` - is used.

        .. versionchanged:: 4.17.0
           The keyword arguments can now be set per plot by using a
           sequence of values. The layout can be changed with the
           rows and cols arguments and the automatic calculation
           no longer forces two rows. Handling of the overplot flag
           has been improved.

        .. versionchanged:: 4.15.0
           A number of labels, such as "bkgfit", are marked as
           deprecated and using them will cause a warning message to
           be displayed, indicating the new label to use.

        .. versionchanged:: 4.12.2
           Keyword arguments, such as alpha and ylog, can be sent to
           each plot.

        Parameters
        ----------
        args
           The plot names and identifiers.
        rows, cols
           The number of rows and columns (if set).
        kwargs
           The plot arguments applied to each plot.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           The label is invalid.

        See Also
        ---------
        get_default_id, get_split_plot, set_xlinear, set_xlog,
        set_ylinear, set_ylog

        Notes
        -----
        The supported plot types depend on the data set type, and
        include the following list. There are also individual
        functions, with ``plot_`` prepended to the plot type, such as
        `plot_data`. There are also several multiple-plot commands,
        such as `plot_fit_ratio`, `plot_fit_resid`, and
        `plot_fit_delchi`.

        ``arf``
           The ARF for the data set (only for `DataPHA` data sets).

        ``bkg``
           The background.

        ``bkg_chisqr``
           The chi-squared statistic calculated for each bin when
           fitting the background.

        ``bkg_delchi``
           The residuals for each bin, calculated as (data-model)
           divided by the error, for the background.

        ``bkg_fit``
           The data (as points) and the convolved model (as a line),
           for the background data set.

        ``bkg_model``
           The convolved background model.

        ``bkg_ratio``
           The residuals for each bin, calculated as data/model,
           for the background data set.

        ``bkg_resid``
           The residuals for each bin, calculated as (data-model),
           for the background data set.

        ``bkg_source``
           The un-convolved background model.

        ``chisqr``
           The chi-squared statistic calculated for each bin.

        ``data``
           The data (which may be background subtracted).

        ``delchi``
           The residuals for each bin, calculated as (data-model)
           divided by the error.

        ``fit``
           The data (as points) and the convolved model (as a line).

        ``kernel``
           The PSF kernel associated with the data set.

        ``model``
           The convolved model.

        ``model_component``
           Part of the full model expression (convolved).

        ``model_components``
           Parts of the full model expression (convolved).

        ``order``
          Plot the model for a selected response

        ``psf``
           The unfiltered PSF kernel associated with the data set.

        ``ratio``
           The residuals for each bin, calculated as data/model.

        ``resid``
           The residuals for each bin, calculated as (data-model).

        ``source``
           The un-convolved model.

        ``source_component``
           Part of the full model expression (un-convolved).

        ``source_components``
           Parts of the full model expression (un-convolved).

        The plots can be specialized for a particular data type,
        such as the `set_analysis` command controlling the units
        used for PHA data sets.

        Given a plot name, such as "data", the remaining arguments up
        to the next plot name match those from the corresponding
        plot_xxx call (in this case plot_data), ignoring the replot,
        overplot, and clearwindow arguments. So the call

        >>> plot("data", "bkg", 1, "up", ylog=True)

        can be thought of as combining the plots created by calling
        plot_data(ylog=True) and plot_bkg(1, "up", ylog=True).

        The plot capabilities depend on what plotting backend, if any,
        is installed. If there is none available, a warning message
        will be displayed when `sherpa.ui` or `sherpa.astro.ui` is
        imported, and the `plot` set of commands will not create any
        plots. The choice of back end is made by changing the
        ``options.plot_pkg`` setting in the Sherpa configuration file.

        The keyword arguments are sent to each plot (so care must be
        taken to ensure they are valid for all plots).

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

        Draw the data and model plots both with a log-scale for the
        y axis:

        >>> plot("data", "model", ylog=True)

        Plot the background data components "up" and "down" for
        dataset 1:

        >>> plot("bkg", 1, "up", "bkg", 1, "down")

        Draw both data and model for the default dataset in black, but
        with partial opacity:

        >>> plot("data", "model", color="black", alpha=0.5)

        Draw the two plots in black but with different opacities:

        >>> plot("data", "model", color="black", alpha=[1, 0.5])

        Label each plot (the output depends on the backend and the
        plot options):

        >>> plot("data", "model", label=["data", "model"])

        Draw the two plots with different y-axis scalings:

        >>> plot("data", 2, "model", 2, ylog=[False, True])

        Change the layout to a single column of plots:

        >>> plot("data", "data", 2, cols=1)

        Use a two-column by three-row display (although in this case
        only one of the rows or cols arguments needed to be given):

        >>> plot("data", "data", 2, "model", "model", 2,
        ...      "resid", "resid", 2, rows=3, cols=2)

        Create a display for three plots, vertically aligned,
        but only display plots in the first two:

        >>> plot("data", "model", cols=1, rows=3)

        Draw the data and residuals for the default dataset and then
        overplot those from dataset 2:

        >>> plot("data", "resid", cols=1, color="black")
        >>> plot("data", 2, "resid", 2, overplot=True, color="black", alpha=0.5)

        """

        rows = kwargs.pop("rows", None)
        cols = kwargs.pop("cols", None)
        self._multi_plot(args, rows=rows, cols=cols, **kwargs)

    def plot_data(self,
                  id: Optional[IdType] = None,
                  replot=False, overplot=False,
                  clearwindow=True, **kwargs) -> None:
        r"""Plot the data values.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_data`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

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

        Notes
        -----
        The additional arguments supported by `plot_data` are the same
        as the keywords of the dictionary returned by `get_data_plot_prefs`.

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

        The following example requires that the Matplotlib backend
        is selected, and uses a Matplotlib function to create a
        subplot (in this case one filling the bottom half of the
        plot area) and then calls `plot_data` with the `clearwindow`
        argument set to `False` to use this subplot. If the
        `clearwindow` argument had not been used then the plot area
        would have been cleared and the plot would have filled the
        area.

        >>> plt.subplot(2, 1, 2)
        >>> plot_data(clearwindow=False)

        Additional arguments can be given that are passed to the
        plot backend: the supported arguments match the keywords
        of the dictionary returned by `get_data_plot_prefs`. Examples
        include (for the Matplotlib backend): adding a "cap" to
        the error bars:

        >>> plot_data(capsize=4)

        changing the symbol to a square:

        >>> plot_data(marker='s')

        using a dotted line to connect the points:

        >>> plot_data(linestyle='dotted')

        and plotting multiple data sets on the same plot, using
        a log scale for the Y axis, setting the alpha transparency
        for each plot, and explicitly setting the
        colors of the last two datasets:

        >>> plot_data(ylog=True, alpha=0.7)
        >>> plot_data(2, overplot=True, alpha=0.7, color='brown')
        >>> plot_data(3, overplot=True, alpha=0.7, color='purple')

        Set the labels used for the X and Y axes for the data. In this
        example the matplotlib backend is used and so the LaTeX
        support is used to display an Angstrom symbol as part of the X
        axis label. Note that the labels will be retained for other
        plots, including other plot types such as plot_model() or
        plot_fit_resid().

        >>> d = get_data()
        >>> d.set_xlabel(r"x axis [$\AA$]")
        >>> d.set_ylabel("y axis")
        >>> plot_data()

        """

        plotobj = self.get_data_plot(id, recalc=not replot)
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    # DOC-NOTE: also in sherpa.astro.utils
    #  - we include a description of the DataPHA handling here
    #    even though its only relevant to sherpa.astro.ui
    #
    def plot_model(self,
                   id: Optional[IdType] = None,
                   replot=False, overplot=False,
                   clearwindow=True, **kwargs) -> None:
        """Plot the model for a data set.

        This function plots the model for a data set, which includes
        any instrument response (e.g. a convolution created by
        `set_psf`).

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_model`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        See Also
        --------
        get_model_plot : Return the data used to create the model plot.
        get_model_plot_prefs : Return the preferences for plot_model.
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
        The additional arguments supported by `plot_model` are the same
        as the keywords of the dictionary returned by `get_model_plot_prefs`.

        For PHA data sets the model plot created by `plot_model` differs to
        the model plot created by `plot_fit`: the fit version uses the grouping
        of the data set whereas the `plot_model` version shows the ungrouped
        data (that is, it uses the instrumental grid). The filters used are the
        same in both cases.

        Examples
        --------

        Plot the convolved source model for the default data set:

        >>> plot_model()

        Overplot the model for data set 2 on data set 1:

        >>> plot_model(1)
        >>> plot_model(2, overplot=True)

        Create the equivalent of ``plot_fit('jet')``:

        >>> plot_data('jet')
        >>> plot_model('jet', overplot=True)

        Additional arguments can be given that are passed to the
        plot backend: the supported arguments match the keywords
        of the dictionary returned by `get_model_plot_prefs`. The
        following plots the model using a log scale for both axes,
        and then overplots the model from data set 2 using a dashed
        line and slightly transparent:

        >>> plot_model(xlog=True, ylog=True)
        >>> plot_model(2, overplot=True, alpha=0.7, linestyle='dashed')

        """

        plotobj = self.get_model_plot(id, recalc=not replot)
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    # DOC-NOTE: also in sherpa.astro.utils, for now copies this text
    #           but does the astro version support a bkg_id parameter?
    def plot_source_component(self, id, model=None, replot=False,
                              overplot=False, clearwindow=True,
                              **kwargs) -> None:
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
           Set to ``True`` to use the values calculated by the last
           call to `plot_source_component`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        See Also
        --------
        get_source_component_plot : Return the data used by plot_source_component.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_model_component : Plot a component of the model for a data set.
        plot_source_components : Plot all the components of a source.
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

        The additional keyword arguments match the keywords of the
        dictionary returned by get_model_plot_prefs.

        Examples
        --------

        Overplot the ``pl`` component of the source expression for
        the default data set:

        >>> plot_source()
        >>> plot_source_component(pl, overplot=True)

        """

        if model is None:
            id, model = model, id

        plotobj = self.get_source_component_plot(id, model,
                                                 recalc=not replot)

        # Note: if replot=True then the value of model is ignored,
        #       which is probably surprising to users
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def plot_source_components(self,
                               id: Optional[IdType] = None,
                               overplot=False,
                               clearwindow=True,
                               **kwargs) -> None:
        """Plot all the components of a source.

        Display the individual components of a source expression.

        .. versionchanged:: 4.17.0
           The keyword arguments can now be set per plot by using a
           sequence of values.

        .. versionadded:: 4.16.1

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot,
           otherwise create a new plot. The default is ``False``. This
           is only used for the first component.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating
           this new plot (e.g. for multi-panel plots)? This is only
           used for the first component.

        See Also
        --------
        get_source_component_plots, plot_model_components, plot_source_component

        Notes
        -----
        The additional keyword arguments match the keywords of the
        dictionary returned by get_model_plot_prefs.

        Examples
        --------

        Display two plots, the first for the `gal * pl` component and
        the second for `gal * line`:

        >>> set_source(xsphabs.gal * (powlaw1d.pl + xsgaussian.line))
        >>> plot_source_components(alpha=0.6)

        Plot the combined source and then overplot the two components
        in black, partly opaque, and using dotted and dashed line
        styles:

        >>> plot_source(label="combined")
        >>> plot_source_components(overplot=True, color="black",
        ...                        linestyle=["dotted", "dashed"],
        ...                        label=["model 1", "model 2"],
        ...                        alpha=0.5)

        """

        plots = self.get_source_components_plot(id)
        self._plot(plots, overplot=overplot,
                   clearwindow=clearwindow, **kwargs)

    def plot_model_component(self, id, model=None, replot=False,
                             overplot=False, clearwindow=True,
                             **kwargs) -> None:
        """Plot a component of the model for a data set.

        This function evaluates and plots a component of the model
        expression for a data set, including any instrument response.
        Use `plot_source_component` to display without any response.
        For PHA data, the response model is automatically added by the
        routine unless the model contains a response.

        Parameters
        ----------
        id : int or str, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        model : str or sherpa.models.model.Model instance
           The component to display (the name, if a string).
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_model_component`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        See Also
        --------
        get_model_component_plot : Return the data used to create the model-component plot.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_model_components : Plot all the components of a model.
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

        The additional keyword arguments match the keywords of the
        dictionary returned by get_model_plot_prefs.

        Examples
        --------

        Overplot the ``pl`` component of the model expression for
        the default data set:

        >>> plot_model()
        >>> plot_model_component(pl, overplot=True)

        Display the results for the 'jet' data set (data and model),
        and then overplot the ``pl`` component evaluated for the 'jet'
        and 'core' data sets:

        >>> plot_fit('jet')
        >>> plot_model_component('jet', pl, overplot=True)
        >>> plot_model_component('core', pl, overplot=True)

        For PHA data sets the response is automatically added, but it
        can also be explicitly included, which will create the same
        plot:

        >>> plot_model_component(pl)
        >>> rsp = get_response()
        >>> plot_model_component(rsp(pl))

        """

        if model is None:
            id, model = model, id

        plotobj = self.get_model_component_plot(id, model,
                                                recalc=not replot)

        # Note: if replot=True then the value of model is ignored,
        #       which is probably surprising to users
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def plot_model_components(self,
                              id: Optional[IdType] = None,
                              overplot=False,
                              clearwindow=True,
                              **kwargs) -> None:
        """Plot all the components of a model.

        Display the individual model components of a source expression.

        .. versionchanged:: 4.17.0
           The keyword arguments can now be set per plot by using a
           sequence of values.

        .. versionadded:: 4.16.1

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot,
           otherwise create a new plot. The default is ``False``. This
           is only used for the first component.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating
           this new plot (e.g. for multi-panel plots)? This is only
           used for the first component.

        See Also
        --------
        get_model_components_plot, plot_model_component, plot_source_components

        Notes
        -----
        The additional keyword arguments match the keywords of the
        dictionary returned by get_model_plot_prefs.

        Examples
        --------

        Display two plots, the first for the `gal * pl` component and
        the second for `gal * line`:

        >>> set_source(xsphabs.gal * (powlaw1d.pl + xsgaussian.line))
        >>> plot_model_components(alpha=0.6)

        Plot the combined model and then overplot the two components
        in black, partly opaque, and using dotted and dashed line
        styles:

        >>> plot_model(label="combined")
        >>> plot_model_components(overplot=True, color="black",
        ...                       linestyle=["dotted", "dashed"],
        ...                       label=["model 1", "model 2"],
        ...                       alpha=0.5)

        """

        plots = self.get_model_components_plot(id)
        self._plot(plots, overplot=overplot,
                   clearwindow=clearwindow, **kwargs)

    # DOC-NOTE: also in sherpa.astro.utils, but with extra lo/hi arguments
    def plot_source(self,
                    id: Optional[IdType] = None,
                    replot=False,
                    overplot=False, clearwindow=True,
                    **kwargs) -> None:
        """Plot the source expression for a data set.

        This function plots the source model for a data set. This does
        not include any instrument response (e.g. a convolution
        created by `set_psf`).

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_source`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        See Also
        --------
        get_source_plot : Return the data used to create the source plot.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_model : Plot the model for a data set.
        plot_source_component : Plot a component of the source expression for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Notes
        -----
        The additional arguments supported by `plot_source` are the same
        as the keywords of the dictionary returned by `get_model_plot_prefs`.

        Examples
        --------

        Plot the unconvolved source model for the default data set:

        >>> plot_source()

        Overplot the source model for data set 2 on data set 1:

        >>> plot_source(1)
        >>> plot_source(2, overplot=True)

        Additional arguments can be given that are passed to the
        plot backend: the supported arguments match the keywords
        of the dictionary returned by `get_model_plot_prefs`. The
        following plots the source using a log scale for both axes,
        and then overplots the source from data set "jet" using a dashed
        line:

        >>> plot_source(xlog=True, ylog=True)
        >>> plot_source('jet', overplot=True, linestyle='dashed')

        """

        # Add a warning here for plot_source, otherwise the error is
        # about get_source_plot. Perhaps we should just use a single
        # error message?
        #
        idval = self._fix_id(id)
        mdl = self._models.get(idval, None)
        if mdl is not None:
            raise IdentifierErr(f"Convolved model\n'{mdl.name}'\n"
                                f" is set for dataset {idval}. "
                                "You should use plot_model instead.")

        plotobj = self.get_source_plot(idval, recalc=not replot)
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def plot_fit(self,
                 id: Optional[IdType] = None,
                 replot=False, overplot=False,
                 clearwindow=True,
                 **kwargs) -> None:
        """Plot the fit results (data, model) for a data set.

        This function creates a plot containing the data and the model
        (including any instrument response) for a data set.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_fit`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_fit_plot : Return the data used to create the fit plot.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_fit_delchi : Plot the fit results, and the residuals, for a data set.
        plot_fit_ratio : Plot the fit results, and the ratio of data to model, for a data set.
        plot_fit_resid : Plot the fit results, and the residuals, for a data set.
        plot_data : Plot the data values.
        plot_model : Plot the model for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Notes
        -----
        The additional arguments supported by `plot_fit` are the same
        as the keywords of the dictionary returned by `get_data_plot_prefs`.

        Examples
        --------

        Plot the fit results for the default data set:

        >>> plot_fit()

        Overplot the 'core' results on those from the 'jet' data set,
        using a logarithmic scale for the X axis:

        >>> set_xlog()
        >>> plot_fit('jet')
        >>> plot_fit('core', overplot=True)

        Keyword arguments can be given to override the plot preferences;
        for example the following sets the y axis to a log scale, but
        only for this plot:

        >>> plot_fit(ylog=True)

        The color can be changed for both the data and model using
        (note that the keyword name and supported values depends on
        the plot backend; this example assumes that Matplotlib is being
        used):

        >>> plot_fit(color='orange')

        Draw the fits for two datasets, setting the second one partially
        transparent (this assumes Matplotlib is used):

        >>> plot_fit(1)
        >>> plot_fit(2, alpha=0.7, overplot=True)

        """

        plotobj = self.get_fit_plot(id, recalc=not replot)
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def plot_resid(self,
                   id: Optional[IdType] = None,
                   replot=False, overplot=False,
                   clearwindow=True,
                   **kwargs) -> None:
        """Plot the residuals (data - model) for a data set.

        This function displays the residuals (data - model) for a data
        set.

        .. versionchanged:: 4.12.0
           The Y axis is now always drawn using a linear scale.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_resid`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

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

        Notes
        -----
        The additional arguments supported by `plot_resid` are the same
        as the keywords of the dictionary returned by `get_data_plot_prefs`.

        The ylog setting is ignored, and the Y axis is drawn using a
        linear scale.

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

        Additional arguments can be given that are passed to the
        plot backend: the supported arguments match the keywords
        of the dictionary returned by `get_data_plot_prefs`. The
        following sets the cap length for the ends of the error bars:

        >>> plot_resid(capsize=5)

        """

        plotobj = self.get_resid_plot(id, recalc=not replot)
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def plot_chisqr(self,
                    id: Optional[IdType] = None,
                    replot=False, overplot=False,
                    clearwindow=True,
                    **kwargs) -> None:
        """Plot the chi-squared value for each point in a data set.

        This function displays the square of the residuals (data -
        model) divided by the error, for a data set.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_chisqr`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

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

        plotobj = self.get_chisqr_plot(id, recalc=not replot)
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def plot_delchi(self,
                    id: Optional[IdType] = None,
                    replot=False, overplot=False,
                    clearwindow=True,
                    **kwargs) -> None:
        """Plot the ratio of residuals to error for a data set.

        This function displays the residuals (data - model) divided by
        the error, for a data set.

        .. versionchanged:: 4.12.0
           The Y axis is now always drawn using a linear scale.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_delchi`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

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

        Notes
        -----
        The additional arguments supported by `plot_delchi` are the same
        as the keywords of the dictionary returned by `get_data_plot_prefs`.

        The ylog setting is ignored, and the Y axis is drawn using a
        linear scale.

        Examples
        --------

        Plot the residuals for the default data set, divided by
        the error value for each bin:

        >>> plot_delchi()

        Overplot the values from the 'core' data set on those
        from the 'jet' dataset:

        >>> plot_delchi('jet')
        >>> plot_delchi('core', overplot=True)

        Additional arguments can be given that are passed to the
        plot backend: the supported arguments match the keywords
        of the dictionary returned by `get_data_plot_prefs`. The
        following sets error bars to be orange and the marker to
        be a circle (larger than the default one):

        >>> plot_delchi(ecolor='orange', marker='o')

        """

        plotobj = self.get_delchi_plot(id, recalc=not replot)
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def plot_ratio(self,
                   id: Optional[IdType] = None,
                   replot=False, overplot=False,
                   clearwindow=True,
                   **kwargs) -> None:
        """Plot the ratio of data to model for a data set.

        This function displays the ratio data / model for a data set.

        .. versionchanged:: 4.12.0
           The Y axis is now always drawn using a linear scale.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_ratio`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

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

        Notes
        -----
        The additional arguments supported by `plot_ratio` are the same
        as the keywords of the dictionary returned by `get_data_plot_prefs`.

        The ylog setting is ignored, and the Y axis is drawn using a
        linear scale.

        Examples
        --------

        Plot the ratio of data to model for the default data set:

        >>> plot_ratio()

        Overplot the ratios from the 'core' data set on those from the
        'jet' dataset:

        >>> plot_ratio('jet')
        >>> plot_ratio('core', overplot=True)

        Additional arguments can be given that are passed to the
        plot backend: the supported arguments match the keywords
        of the dictionary returned by `get_data_plot_prefs`. The
        following sets the X axis to a log scale and draws a solid
        line between the points:

        >>> plot_ratio(xlog=True, linestyle='solid')

        """

        plotobj = self.get_ratio_plot(id, recalc=not replot)
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def plot_psf(self,
                 id: Optional[IdType] = None,
                 replot=False, overplot=False,
                 clearwindow=True,
                 **kwargs) -> None:
        """Plot the 1D PSF model applied to a data set.

        The `plot_kernel` function shows the data used to convolve
        the model.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_psf`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

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

        plotobj = self.get_psf_plot(id, recalc=not replot)
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def plot_kernel(self,
                    id: Optional[IdType] = None,
                    replot=False, overplot=False,
                    clearwindow=True,
                    **kwargs) -> None:
        """Plot the 1D kernel applied to a data set.

        The `plot_psf` function shows the full PSF, from which the
        kernel is derived.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_kernel`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

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

        plotobj = self.get_kernel_plot(id, recalc=not replot)
        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def _jointplot2(self, plot1, plot2,
                    overplot=False, clearwindow=True,
                    **kwargs) -> None:
        """Create a joint plot, vertically aligned, fit data on the top.

        Parameters
        ----------
        plot1 : sherpa.plot.Plot instance
           The plot to appear in the top panel.
        plot2 : sherpa.plot.Plot instance
           The plot to appear in the bottom panel.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        """

        self._jointplot.reset()

        with sherpa.plot.backend:

            # Note: the user preferences are set to both plots
            #
            self._jointplot.plottop(plot1, overplot=overplot,
                                    clearwindow=clearwindow, **kwargs)

            # The two plots are intended to have the same scaling
            # on the X axis (log or linear), and the approach is
            # to use log if either of the components of the top
            # plot (i.e. it is assumed that this is a fit plot)
            # use a log scale.
            #
            # This is complicated by the need to access the individual
            # plot preferences of plot1, and then check for different
            # types of plot objects.
            #
            p2prefs = get_plot_prefs(plot2)
            oldval = p2prefs['xlog']
            dprefs = get_plot_prefs(plot1.dataplot)
            mprefs = get_plot_prefs(plot1.modelplot)

            if dprefs['xlog'] or mprefs['xlog']:
                p2prefs['xlog'] = True

            self._jointplot.plotbot(plot2, overplot=overplot, **kwargs)

            p2prefs['xlog'] = oldval

    def plot_fit_resid(self,
                       id: Optional[IdType] = None,
                       replot=False, overplot=False,
                       clearwindow=True,
                       **kwargs) -> None:
        """Plot the fit results, and the residuals, for a data set.

        This creates two plots - the first from `plot_fit` and the
        second from `plot_resid` - for a data set.

        .. versionchanged:: 4.12.2
           The ``overplot`` option now works.

        .. versionchanged:: 4.12.0
           The Y axis of the residual plot is now always drawn using a
           linear scale.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the previous values. The default is
           ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_fit_plot : Return the data used to create the fit plot.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_fit : Plot the fit results for a data set.
        plot_fit_delchi : Plot the fit results, and the residuals, for a data set.
        plot_fit_ratio : Plot the fit results, and the ratio of data to model, for a data set.
        plot_data : Plot the data values.
        plot_model : Plot the model for a data set.
        plot_resid : Plot the residuals (data - model) for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Notes
        -----
        The additional arguments supported by `plot_fit_resid` are the same
        as the keywords of the dictionary returned by `get_data_plot_prefs`,
        and are applied to both plots.

        For the residual plot, the ylog setting is ignored, and the Y axis
        is drawn using a linear scale.

        Examples
        --------

        Plot the results for the default data set:

        >>> plot_fit_resid()

        Overplot the 'core' results on those from the 'jet' data set,
        using a logarithmic scale for the X axis:

        >>> set_xlog()
        >>> plot_fit_resid('jet')
        >>> plot_fit_resid('core', overplot=True)

        Additional arguments can be given that are passed to the
        plot backend: the supported arguments match the keywords
        of the dictionary returned by `get_data_plot_prefs`. The
        following sets the data in both plots to be drawn in a
        blue color, have caps on the error bars, but to only
        draw the y error bars:

        >>> plot_fit_resid(capsize=4, color='skyblue', xerrorbars=False)

        """

        plot1obj = self.get_fit_plot(id, recalc=not replot)
        plot2obj = self.get_resid_plot(id, recalc=not replot)
        self._jointplot2(plot1obj, plot2obj,
                         overplot=overplot, clearwindow=clearwindow,
                         **kwargs)

    def plot_fit_ratio(self,
                       id: Optional[IdType] = None,
                       replot=False, overplot=False,
                       clearwindow=True,
                       **kwargs) -> None:
        """Plot the fit results, and the ratio of data to model, for a data set.

        This creates two plots - the first from `plot_fit` and the
        second from `plot_ratio` - for a data set.

        .. versionchanged:: 4.12.2
           The ``overplot`` option now works.

        .. versionadded:: 4.12.0

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_fit_ratio`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_fit_plot : Return the data used to create the fit plot.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_fit : Plot the fit results for a data set.
        plot_fit_resid : Plot the fit results, and the residuals, for a data set.
        plot_fit_delchi : Plot the fit results, and the residuals, for a data set.
        plot_data : Plot the data values.
        plot_model : Plot the model for a data set.
        plot_ratio : Plot the ratio of data to model for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Notes
        -----
        The additional arguments supported by `plot_fit_ratio` are the same
        as the keywords of the dictionary returned by `get_data_plot_prefs`,
        and are applied to both plots.

        For the ratio plot, the ylog setting is ignored, and the Y axis
        is drawn using a linear scale.

        Examples
        --------

        Plot the results for the default data set:

        >>> plot_fit_ratio()

        Overplot the 'core' results on those from the 'jet' data set,
        using a logarithmic scale for the X axis:

        >>> set_xlog()
        >>> plot_fit_ratio('jet')
        >>> plot_fit_ratio('core', overplot=True)

        Additional arguments can be given that are passed to the
        plot backend: the supported arguments match the keywords
        of the dictionary returned by `get_data_plot_prefs`. The
        following sets the plots to use square symbols (this includes
        the model as well as data in the top plot) and turns off
        any line between plots, when using the Matplotlib backend:

        >>> plot_fit_ratio(marker='s', linestyle='none')

        """

        plot1obj = self.get_fit_plot(id, recalc=not replot)
        plot2obj = self.get_ratio_plot(id, recalc=not replot)
        self._jointplot2(plot1obj, plot2obj,
                         overplot=overplot, clearwindow=clearwindow,
                         **kwargs)

    def plot_fit_delchi(self,
                        id: Optional[IdType] = None,
                        replot=False, overplot=False,
                        clearwindow=True,
                        **kwargs) -> None:
        """Plot the fit results, and the residuals, for a data set.

        This creates two plots - the first from `plot_fit` and the
        second from `plot_delchi` - for a data set.

        .. versionchanged:: 4.12.2
           The ``overplot`` option now works.

        .. versionchanged:: 4.12.0
           The Y axis of the delchi plot is now always drawn using a
           linear scale.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_fit_delchi`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_fit_plot : Return the data used to create the fit plot.
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.
        plot_fit : Plot the fit results for a data set.
        plot_fit_ratio : Plot the fit results, and the ratio of data to model, for a data set.
        plot_fit_resid : Plot the fit results, and the residuals, for a data set.
        plot_data : Plot the data values.
        plot_delchi : Plot the ratio of residuals to error for a data set.
        plot_model : Plot the model for a data set.
        set_xlinear : New plots will display a linear X axis.
        set_xlog : New plots will display a logarithmically-scaled X axis.
        set_ylinear : New plots will display a linear Y axis.
        set_ylog : New plots will display a logarithmically-scaled Y axis.

        Notes
        -----
        The additional arguments supported by `plot_fit_delchi` are the same
        as the keywords of the dictionary returned by `get_data_plot_prefs`,
        and are applied to both plots.

        For the delchi plot, the ylog setting is ignored, and the Y axis
        is drawn using a linear scale.

        Examples
        --------

        Plot the results for the default data set:

        >>> plot_fit_delchi()

        Overplot the 'core' results on those from the 'jet' data set,
        using a logarithmic scale for the X axis:

        >>> set_xlog()
        >>> plot_fit_delchi('jet')
        >>> plot_fit_delchi('core', overplot=True)

        Additional arguments can be given that are passed to the
        plot backend: the supported arguments match the keywords
        of the dictionary returned by `get_data_plot_prefs`. The
        following sets the error bars to be drawn in gray when
        using the Matplotlib backend:

        >>> plot_fit_delchi(ecolor='gray')

        """

        plot1obj = self.get_fit_plot(id, recalc=not replot)
        plot2obj = self.get_delchi_plot(id, recalc=not replot)
        self._jointplot2(plot1obj, plot2obj,
                         overplot=overplot, clearwindow=clearwindow,
                         **kwargs)

    #
    # Statistical plotting routines
    #

    def plot_pdf(self, points, name="x", xlabel="x", bins=12, normed=True,
                 replot=False, overplot=False, clearwindow=True,
                 **kwargs) -> None:
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
           Should the PDF be normalized (the default is ``True``).
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_pdf`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

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

        plotobj = self.get_pdf_plot()
        if not replot:
            plotobj.prepare(points, bins, normed, xlabel, name)

        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def get_pdf_plot(self):
        """Return the data used to plot the last PDF.

        Returns
        -------
        plot : a `sherpa.plot.PDFPlot` instance
           An object containing the data used by the last call to
           `plot_pdf`. The fields will be ``None`` if the function
           has not been called.

        See Also
        --------
        plot_pdf : Plot the probability density function of an array.

        """
        return self._pdfplot

    def plot_cdf(self, points, name="x", xlabel="x",
                 replot=False, overplot=False, clearwindow=True,
                 **kwargs) -> None:
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
           Set to ``True`` to use the values calculated by the last
           call to `plot_cdf`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

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

        plotobj = self.get_cdf_plot()
        if not replot:
            plotobj.prepare(points, xlabel, name)

        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def get_cdf_plot(self):
        """Return the data used to plot the last CDF.

        Returns
        -------
        plot : a `sherpa.plot.CDFPlot` instance
           An object containing the data used by the last call to
           `plot_cdf`. The fields will be ``None`` if the function
           has not been called.

        See Also
        --------
        plot_cdf : Plot the cumulative density function of an array.

        """
        return self._cdfplot

    # DOC-TODO: what does xlabel do?
    def plot_trace(self, points, name="x", xlabel="x",
                   replot=False, overplot=False, clearwindow=True,
                   **kwargs) -> None:
        """Create a trace plot of row number versus value.

        Display a plot of the ``points`` array values (Y axis) versus row
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
           Set to ``True`` to use the values calculated by the last
           call to `plot_trace`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

        See Also
        --------
        get_draws : Run the pyBLoCXS MCMC algorithm.
        get_trace_plot : Return the data used to plot the last trace.
        plot_cdf : Plot the cumulative density function of an array.
        plot_pdf : Plot the probability density function of an array.
        plot_scatter : Create a scatter plot.

        Examples
        --------

        Plot the trace of the 500 elements in the ``x`` array:

        >>> mu, sigma = 100, 15
        >>> x = mu + sigma * np.random.randn(500)
        >>> plot_trace(x)

        Use "ampl" as the Y axis label:

        >>> plot_trace(ampl, name='ampl')

        """

        plotobj = self.get_trace_plot()
        if not replot:
            plotobj.prepare(points, xlabel, name)

        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def get_trace_plot(self):
        """Return the data used to plot the last trace.

        Returns
        -------
        plot : a `sherpa.plot.TracePlot` instance
           An object containing the data used by the last call to
           `plot_trace`. The fields will be ``None`` if the function
           has not been called.

        See Also
        --------
        plot_trace : Create a trace plot of row number versus value.

        """
        return self._traceplot

    def plot_scatter(self, x, y, name="(x,y)", xlabel="x", ylabel="y",
                     replot=False, overplot=False, clearwindow=True,
                     **kwargs) -> None:
        """Create a scatter plot.

        Parameters
        ----------
        x : array
           The values to plot on the X axis.
        y : array
           The values to plot on the Y axis. This must match the size
           of the ``x`` array.
        name : str, optional
           The plot title.
        xlabel : str, optional
           The label for the X axis.
        ylabel : str, optional
           The label for the Y axis.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `plot_scatter`. The default is ``False``.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.
        clearwindow : bool, optional
           Should the existing plot area be cleared before creating this
           new plot (e.g. for multi-panel plots)?

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

        plotobj = self.get_scatter_plot()
        if not replot:
            plotobj.prepare(x, y, xlabel, ylabel, name)

        self._plot(plotobj, overplot=overplot, clearwindow=clearwindow,
                   **kwargs)

    def get_scatter_plot(self):
        """Return the data used to plot the last scatter plot.

        Returns
        -------
        plot : a `sherpa.plot.ScatterPlot` instance
           An object containing the data used by the last call to
           `plot_scatter`. The fields will be ``None`` if the function
           has not been called.

        See Also
        --------
        plot_scatter : Create a scatter plot.

        """
        return self._scatterplot

    #
    # Contours
    #

    def _contour(self, plotobj, overcontour=False, **kwargs) -> None:
        """Display a plot object

        Parameters
        ----------
        plotobj
            The contour plot to display
        **kwargs
            The arguments to the contour method of plotobj.

        See Also
        --------
        _plot

        """

        with sherpa.plot.backend:
            plotobj.contour(overcontour=overcontour, **kwargs)

    def contour(self,
                *args,
                # At present export_method does not support this
                # rows: Optional[int] = None,
                # cols: Optional[int] = None,
                **kwargs) -> None:
        """Create a contour plot for an image data set.

        Create one or more contour plots, depending on the arguments
        it is set: a plot type, followed by an optional data set
        identifier, and this can be repeated. If no data set
        identifier is given for a plot type, the default identifier -
        as returned by `get_default_id` - is used. This is for
        2D data sets.

        .. versionchanged:: 4.17.0
           The keyword arguments can now be set per plot by using a
           sequence of values. The layout can be changed with the
           rows and cols arguments and the automatic calculation
           no longer forces two rows. Handling of the overcontour
           flag has been improved.

        .. versionchanged:: 4.12.2
           Keyword arguments, such as alpha, can be sent to
           each plot.

        Parameters
        ----------
        args
           The contour-plot names and identifiers.
        rows, cols
           The number of rows and columns (if set).
        kwargs
           The plot arguments applied to each contour plot.

        Raises
        ------
        sherpa.utils.err.DataErr
           The data set does not support the requested plot type.

        See Also
        ---------
        contour_data, contour_fit, contour_fit_resid, contour_kernel,
        contour_model, contour_psf, contour_ratio, contour_resid,
        contour_source, get_default_id, get_split_plot

        Notes
        -----
        The supported plot types depend on the data set type, and
        include the following list. There are also individual
        functions, with ``contour_`` prepended to the plot type, such as
        `contour_data` and the `contour_fit_resid` variant:

        ``data``
           The data.

        ``fit``
           Contours of the data and the source model.

        ``fit_resid``
           Two plots: the first is the contours of the data and the
           source model and the second is the residuals.

        ``kernel``
           The kernel.

        ``model``
           The source model including any PSF convolution set by
           `set_psf`.

        ``psf``
           The PSF.

        ``ratio``
           Contours of the ratio image, formed by dividing the data by
           the model.

        ``resid``
           Contours of the residual image, formed by subtracting the
           model from the data.

        ``source``
           The source model (without any PSF convolution set by
           `set_psf`).

        The keyword arguments are sent to each plot (so care must be
        taken to ensure they are valid for all plots).

        Examples
        --------

        >>> contour('data')

        >>> contour('data', 1, 'data', 2)

        >>> contour('data', 'model')

        >>> contour('data', 'model', 'fit', 'resid')

        >>> contour('data', 'model', alpha=0.7)

        Use a single column rather than single row to display the
        contour plots:

        >>> contour('data', 'model', cols=1)

        """
        rows = kwargs.pop("rows", None)
        cols = kwargs.pop("cols", None)
        self._multi_plot(args, plotmeth='contour', rows=rows,
                         cols=cols, **kwargs)

    def contour_data(self,
                     id: Optional[IdType] = None,
                     replot=False, overcontour=False,
                     **kwargs) -> None:
        """Contour the values of an image data set.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `contour_data`. The default is ``False``.
        overcontour : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new contour plot. The default is ``False``.

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

        plotobj = self.get_data_contour(id, recalc=not replot)
        self._contour(plotobj, overcontour=overcontour, **kwargs)

    def contour_model(self,
                      id: Optional[IdType] = None,
                      replot=False, overcontour=False,
                      **kwargs) -> None:
        """Create a contour plot of the model.

        Displays a contour plot of the values of the model,
        evaluated on the data, including any PSF kernel convolution
        (if set).

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the model. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `contour_model`. The default is ``False``.
        overcontour : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new contour plot. The default is ``False``.

        See Also
        --------
        get_model_contour : Return the data used by contour_model.
        get_model_contour_prefs : Return the preferences for contour_model.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        contour_source : Create a contour plot of the unconvolved spatial model.
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

        plotobj = self.get_model_contour(id, recalc=not replot)
        self._contour(plotobj, overcontour=overcontour, **kwargs)

    def contour_source(self,
                       id: Optional[IdType] = None,
                       replot=False, overcontour=False,
                       **kwargs) -> None:
        """Create a contour plot of the unconvolved spatial model.

        Displays a contour plot of the values of the model,
        evaluated on the data, without any PSF kernel convolution
        applied. The preferences are the same as `contour_model`.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the model. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `contour_source`. The default is ``False``.
        overcontour : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new contour plot. The default is ``False``.

        See Also
        --------
        get_source_contour : Return the data used by contour_source.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        contour_model : Create a contour plot of the model.
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

        plotobj = self.get_source_contour(id, recalc=not replot)
        self._contour(plotobj, overcontour=overcontour, **kwargs)

    def contour_fit(self,
                    id: Optional[IdType] = None,
                    replot=False, overcontour=False,
                    **kwargs) -> None:
        """Contour the fit to a data set.

        Overplot the model - including any PSF - on the data. The
        preferences are the same as `contour_data` and
        `contour_model`.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data and model. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `contour_fit`. The default is ``False``.
        overcontour : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new contour plot. The default is ``False``.

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

        plotobj = self.get_fit_contour(id, recalc=not replot)
        self._contour(plotobj, overcontour=overcontour, **kwargs)

    def contour_resid(self,
                      id: Optional[IdType] = None,
                      replot=False, overcontour=False,
                      **kwargs) -> None:
        """Contour the residuals of the fit.

        The residuals are formed by subtracting the current model -
        including any PSF - from the data.  The preferences are the
        same as `contour_data`.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data and model. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `contour_resid`. The default is ``False``.
        overcontour : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new contour plot. The default is ``False``.

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

        plotobj = self.get_resid_contour(id, recalc=not replot)
        self._contour(plotobj, overcontour=overcontour, **kwargs)

    def contour_ratio(self,
                      id: Optional[IdType] = None,
                      replot=False, overcontour=False,
                      **kwargs) -> None:
        """Contour the ratio of data to model.

        The ratio image is formed by dividing the data by the current
        model, including any PSF. The preferences are the same as
        `contour_data`.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data and model. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `contour_ratio`. The default is ``False``.
        overcontour : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new contour plot. The default is ``False``.

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

        plotobj = self.get_ratio_contour(id, recalc=not replot)
        self._contour(plotobj, overcontour=overcontour, **kwargs)

    def contour_psf(self,
                    id: Optional[IdType] = None,
                    replot=False, overcontour=False,
                    **kwargs) -> None:
        """Contour the PSF applied to the model of an image data set.

        If the data set has no PSF applied to it, the model is
        displayed.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the model. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `contour_psf`. The default is ``False``.
        overcontour : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new contour plot. The default is ``False``.

        See Also
        --------
        get_psf_contour : Return the data used by contour_psf.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        contour_kernel : Contour the kernel applied to the model of an image data set.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.
        set_psf : Add a PSF model to a data set.

        """

        plotobj = self.get_psf_contour(id, recalc=not replot)
        self._contour(plotobj, overcontour=overcontour, **kwargs)

    def contour_kernel(self,
                       id: Optional[IdType] = None,
                       replot=False, overcontour=False,
                       **kwargs) -> None:
        """Contour the kernel applied to the model of an image data set.

        If the data set has no PSF applied to it, the model is
        displayed.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the model. If not given then the
           default identifier is used, as returned by `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `contour_kernel`. The default is ``False``.
        overcontour : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new contour plot. The default is ``False``.

        See Also
        --------
        get_psf_contour : Return the data used by contour_psf.
        get_default_id : Return the default data set identifier.
        contour : Create one or more plot types.
        contour_psf : Contour the PSF applied to the model of an image data set.
        sherpa.astro.ui.set_coord : Set the coordinate system to use for image analysis.
        set_psf : Add a PSF model to a data set.

        """

        plotobj = self.get_kernel_contour(id, recalc=not replot)
        self._contour(plotobj, overcontour=overcontour, **kwargs)

    def contour_fit_resid(self,
                          id: Optional[IdType] = None,
                          replot=False, overcontour=False,
                          **kwargs) -> None:
        """Contour the fit and the residuals to a data set.

        Overplot the model - including any PSF - on the data. In a
        separate plot contour the residuals. The preferences are the
        same as `contour_data` and `contour_model`.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set that provides the data and model. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `contour_fit_resid`. The default is ``False``.
        overcontour : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new contour plot. The default is ``False``.

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

        plot1obj = self.get_fit_contour(id, recalc=not replot)
        plot2obj = self.get_resid_contour(id, recalc=not replot)

        # This does not use _jointplot because the X axis is not
        # obviously shared between the two plots.
        #
        self._splitplot.reset(rows=2, cols=1)
        with sherpa.plot.backend:

            # Note: the user settings are applied to both contours
            self._splitplot.addcontour(plot1obj, overcontour=overcontour, **kwargs)
            self._splitplot.addcontour(plot2obj, overcontour=overcontour, **kwargs)


    ###########################################################################
    # Projection and uncertainty plots
    ###########################################################################

    # DOC-NOTE: I am not convinced that this code is working when recalc=True
    # DOC-NOTE: needs to support the fast option of int_proj
    def get_int_proj(self, par=None,
                     id: Optional[IdType] = None,
                     otherids: Optional[Sequence[IdType]] = None,
                     recalc=False,
                     fast=True, min=None, max=None, nloop=20, delv=None, fac=1,
                     log=False, numcores=None):
        """Return the interval-projection object.

        This returns (and optionally calculates) the data used to
        display the `int_proj` plot. Note that if the the `recalc`
        parameter is `False` (the default value) then all other parameters
        are ignored and the results of the last `int_proj` call are
        returned.

        .. versionchanged:: 4.16.1
           The log parameter can now be set to `True`.

        Parameters
        ----------
        par
           The parameter to plot. This argument is only used if `recalc` is
           set to `True`.
        id : str, int, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : list of str or int, or None, optional
           Other data sets to use in the calculation.
        recalc : bool, optional
           The default value (``False``) means that the results from the
           last call to `int_proj` (or `get_int_proj`) are returned,
           ignoring *all other* parameter values. Otherwise, the
           statistic curve is re-calculated, but not plotted.
        fast : bool, optional
           If ``True`` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is ``False``.
        min : number, optional
           The minimum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : number, optional
           The maximum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : int, optional
           The number of steps to use. This is used when `delv` is set
           to ``None``.
        delv : number, optional
           The step size for the parameter. Setting this overrides
           the `nloop` parameter. The default is ``None``.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : bool, optional
           Should the step size be logarithmically spaced? The
           default (``False``) is to use a linear grid.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        iproj : a `sherpa.plot.IntervalProjection` instance
           The fields of this object can be used to re-create the plot
           created by `int_proj`.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
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

        Since the `recalc` parameter has not been changed to `True`, the
        following will return the results for the last call to `int_proj`,
        which may not have been for the src.ypos parameter:

        >>> iproj = get_int_proj(src.ypos)

        Create the data without creating a plot:

        >>> iproj = get_int_proj(pl.gamma, recalc=True)

        Specify the range and step size for the parameter, in this
        case varying linearly between 12 and 14 with 51 values:

        >>> iproj = get_int_proj(src.r0, id="src", min=12, max=14,
        ...                      nloop=51, recalc=True)

        """

        plotobj = self._intproj
        if recalc:
            par = self._check_par(par)
            if otherids is None:
                otherids = ()
            ids, fit = self._get_fit(id, otherids)
            plotobj.prepare(min=min, max=max, nloop=nloop, delv=delv,
                            fac=fac, log=log, numcores=numcores)
            plotobj.calc(fit, par, self._methods)

        return plotobj

    # DOC-NOTE: Check that this works (since get_int_proj may not) when
    # recalc=True
    def get_int_unc(self, par=None,
                    id: Optional[IdType] = None,
                    otherids: Optional[Sequence[IdType]] = None,
                    recalc=False,
                    min=None, max=None, nloop=20, delv=None, fac=1, log=False,
                    numcores=None):
        """Return the interval-uncertainty object.

        This returns (and optionally calculates) the data used to
        display the `int_unc` plot. Note that if the the `recalc`
        parameter is `False` (the default value) then all other parameters
        are ignored and the results of the last `int_unc` call are
        returned.

        .. versionchanged:: 4.16.1
           The log parameter can now be set to `True`.

        Parameters
        ----------
        par
           The parameter to plot. This argument is only used if `recalc` is
           set to `True`.
        id : str, int, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : list of str or int, or None, optional
           Other data sets to use in the calculation.
        recalc : bool, optional
           The default value (``False``) means that the results from the
           last call to `int_proj` (or `get_int_proj`) are returned,
           ignoring *all other* parameter values. Otherwise, the
           statistic curve is re-calculated, but not plotted.
        min : number, optional
           The minimum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : number, optional
           The maximum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : int, optional
           The number of steps to use. This is used when `delv` is set
           to ``None``.
        delv : number, optional
           The step size for the parameter. Setting this overrides
           the `nloop` parameter. The default is ``None``.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : bool, optional
           Should the step size be logarithmically spaced? The
           default (``False``) is to use a linear grid.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        iunc : a `sherpa.plot.IntervalUncertainty` instance
           The fields of this object can be used to re-create the plot
           created by `int_unc`.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
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

        Since the `recalc` parameter has not been changed to `True`, the
        following will return the results for the last call to `int_unc`,
        which may not have been for the src.ypos parameter:

        >>> iunc = get_int_unc(src.ypos)

        Create the data without creating a plot:

        >>> iunc = get_int_unc(pl.gamma, recalc=True)

        Specify the range and step size for the parameter, in this
        case varying linearly between 12 and 14 with 51 values:

        >>> iunc = get_int_unc(src.r0, id="src", min=12, max=14,
        ...                    nloop=51, recalc=True)

        """

        plotobj = self._intunc
        if recalc:
            par = self._check_par(par)
            if otherids is None:
                otherids = ()
            ids, fit = self._get_fit(id, otherids)
            plotobj.prepare(min=min, max=max, nloop=nloop, delv=delv,
                            fac=fac, log=sherpa.utils.bool_cast(log),
                            numcores=numcores)
            plotobj.calc(fit, par)

        return plotobj

    def get_reg_proj(self, par0=None, par1=None,
                     id: Optional[IdType] = None,
                     otherids: Optional[Sequence[IdType]] = None,
                     recalc=False, fast=True, min=None, max=None,
                     nloop=(10, 10), delv=None, fac=4, log=(False, False),
                     sigma=(1, 2, 3), levels=None, numcores=None):
        """Return the region-projection object.

        This returns (and optionally calculates) the data used to
        display the `reg_proj` contour plot. Note that if the the `recalc`
        parameter is `False` (the default value) then all other parameters
        are ignored and the results of the last `reg_proj` call are
        returned.

        .. versionchanged:: 4.16.1
           The log parameter can now be set to `True` for one or both
           parameters.

        Parameters
        ----------
        par0, par1
           The parameters to plot on the X and Y axes, respectively.
           These arguments are only used if recalc is set to `True`.
        id : str, int, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : list of str or int, or None, optional
           Other data sets to use in the calculation.
        recalc : bool, optional
           The default value (``False``) means that the results from the
           last call to `reg_proj` (or `get_reg_proj`) are returned,
           ignoring *all other* parameter values. Otherwise, the
           statistic curve is re-calculated, but not plotted.
        fast : bool, optional
           If ``True`` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is ``False``.
        min : pair of numbers, optional
           The minimum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : pair of number, optional
           The maximum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : pair of int, optional
           The number of steps to use. This is used when `delv` is set
           to ``None``.
        delv : pair of number, optional
           The step size for the parameter. Setting this overrides
           the `nloop` parameter. The default is ``None``.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : pair of bool, optional
           Should the step size be logarithmically spaced? The
           default (``False``) is to use a linear grid.
        sigma : sequence of number, optional
           The levels at which to draw the contours. The units are the
           change in significance relative to the starting value,
           in units of sigma.
        levels : sequence of number, optional
           The numeric values at which to draw the contours. This
           overrides the `sigma` parameter, if set (the default is
           ``None``).
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        rproj : a `sherpa.plot.RegionProjection` instance
           The fields of this object can be used to re-create the plot
           created by `reg_proj`.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        int_proj : Calculate and plot the fit statistic versus fit parameter value.
        int_unc : Calculate and plot the fit statistic versus fit parameter value.
        reg_proj : Plot the statistic value as two parameters are varied.
        reg_unc : Plot the statistic value as two parameters are varied.

        Examples
        --------

        Return the results for the `reg_proj` run for the `xpos` and `ypos`
        parameters of the `src` component, for the default data set:

        >>> reg_proj(src.xpos, src.ypos)
        >>> rproj = get_reg_proj()

        Since the `recalc` parameter has not been changed to `True`, the
        following will return the results for the last call to `reg_proj`,
        which may not have been for the r0 and alpha parameters:

        >>> rprog = get_reg_proj(src.r0, src.alpha)

        Create the data without creating a plot:

        >>> rproj = get_reg_proj(pl.gamma, gal.nh, recalc=True)

        Specify the range and step size for both the parameters,
        in this case pl.gamma should vary between 0.5 and 2.5, with
        gal.nh between 0.01 and 1, both with 51 values and the
        nH range done over a log scale:

        >>> rproj = get_reg_proj(pl.gamma, gal.nh, id="src",
        ...                      min=(0.5, 0.01), max=(2.5, 1),
        ...                      nloop=(51, 51), log=(False, True),
        ...                      recalc=True)

        """

        plotobj = self._regproj
        if recalc:
            par0 = self._check_par(par0, 'par0')
            par1 = self._check_par(par1, 'par1')
            if otherids is None:
                otherids = ()
            ids, fit = self._get_fit(id, otherids)
            plotobj.prepare(fast=fast, min=min, max=max,
                            nloop=nloop, delv=delv, fac=fac,
                            log=log, sigma=sigma, levels=levels,
                            numcores=numcores)
            plotobj.calc(fit, par0, par1, self._methods)

        return plotobj

    def get_reg_unc(self, par0=None, par1=None,
                    id: Optional[IdType] = None,
                    otherids: Optional[Sequence[IdType]] = None,
                    recalc=False, min=None, max=None, nloop=(10, 10),
                    delv=None, fac=4, log=(False, False), sigma=(1, 2, 3),
                    levels=None, numcores=None):
        """Return the region-uncertainty object.

        This returns (and optionally calculates) the data used to
        display the `reg_unc` contour plot.  Note that if the the `recalc`
        parameter is `False` (the default value) then all other parameters
        are ignored and the results of the last `reg_unc` call are
        returned.

        .. versionchanged:: 4.16.1
           The log parameter can now be set to `True` for one or both
           parameters.

        Parameters
        ----------
        par0, par1
           The parameters to plot on the X and Y axes, respectively.
           These arguments are only used if `recalc` is set to `True`.
        id : str, int, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : list of str or int, or None, optional
           Other data sets to use in the calculation.
        recalc : bool, optional
           The default value (``False``) means that the results from the
           last call to `reg_unc` (or `get_reg_unc`) are returned,
           ignoring *all other* parameter values. Otherwise, the
           statistic curve is re-calculated, but not plotted.
        fast : bool, optional
           If ``True`` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is ``False``.
        min : pair of numbers, optional
           The minimum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : pair of number, optional
           The maximum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : pair of int, optional
           The number of steps to use. This is used when `delv` is set
           to ``None``.
        delv : pair of number, optional
           The step size for the parameter. Setting this overrides
           the `nloop` parameter. The default is ``None``.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : pair of bool, optional
           Should the step size be logarithmically spaced? The
           default (``False``) is to use a linear grid.
        sigma : sequence of number, optional
           The levels at which to draw the contours. The units are the
           change in significance relative to the starting value,
           in units of sigma.
        levels : sequence of number, optional
           The numeric values at which to draw the contours. This
           overrides the `sigma` parameter, if set (the default is
           ``None``).
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.

        Returns
        -------
        rproj : a `sherpa.plot.RegionUncertainty` instance
           The fields of this object can be used to re-create the plot
           created by `reg_unc`.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
        covar : Estimate the confidence intervals using the covariance method.
        int_proj : Calculate and plot the fit statistic versus fit parameter value.
        int_unc : Calculate and plot the fit statistic versus fit parameter value.
        reg_proj : Plot the statistic value as two parameters are varied.
        reg_unc : Plot the statistic value as two parameters are varied.

        Examples
        --------

        Return the results for the `reg_unc` run for the `xpos` and `ypos`
        parameters of the `src` component, for the default data set:

        >>> reg_unc(src.xpos, src.ypos)
        >>> runc = get_reg_unc()

        Since the `recalc` parameter has not been changed to `True`, the
        following will return the results for the last call to `reg_unc`,
        which may not have been for the r0 and alpha parameters:

        >>> runc = get_reg_unc(src.r0, src.alpha)

        Create the data without creating a plot:

        >>> runc = get_reg_unc(pl.gamma, gal.nh, recalc=True)

        Specify the range and step size for both the parameters,
        in this case pl.gamma should vary between 0.5 and 2.5, with
        gal.nh between 0.01 and 1, both with 51 values and the
        nH range done over a log scale:

        >>> runc = get_reg_unc(pl.gamma, gal.nh, id="src",
        ...                    min=(0.5, 0.01), max=(2.5, 1),
        ...                    nloop=(51, 51), log=(False, True),
        ...                    recalc=True)

        """

        plotobj = self._regunc
        if recalc:
            par0 = self._check_par(par0, 'par0')
            par1 = self._check_par(par1, 'par1')

            if otherids is None:
                otherids = ()
            ids, fit = self._get_fit(id, otherids)
            plotobj.prepare(min=min, max=max,
                            nloop=nloop, delv=delv, fac=fac,
                            log=log, sigma=sigma, levels=levels,
                            numcores=numcores)
            plotobj.calc(fit, par0, par1)

        return plotobj

    # DOC-NOTE: I am not convinced I have fac described correctly
    # DOC-NOTE: same synopsis as int_unc
    def int_proj(self, par,
                 id: Optional[IdType] = None,
                 otherids: Optional[Sequence[IdType]] = None,
                 replot=False, fast=True,
                 min=None, max=None, nloop=20, delv=None, fac=1, log=False,
                 numcores=None,
                 overplot=False) -> None:
        """Calculate and plot the fit statistic versus fit parameter value.

        Create a confidence plot of the fit statistic as a function of
        parameter value. Dashed lines are added to indicate the
        current statistic value and the parameter value at this
        point. The parameter value is varied over a grid of points and
        the free parameters re-fit. It is expected that this is run
        after a successful fit, so that the parameter values identify
        the best-fit location.

        .. versionchanged:: 4.16.1
           The log parameter can now be set to `True`.

        Parameters
        ----------
        par
           The parameter to plot.
        id : int, str, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : sequence of int or str, or None, optional
           Other data sets to use in the calculation.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `int_proj`. The default is ``False``.
        fast : bool, optional
           If ``True`` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is ``False``.
        min : number, optional
           The minimum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : number, optional
           The maximum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : int, optional
           The number of steps to use. This is used when `delv` is set
           to ``None``.
        delv : number, optional
           The step size for the parameter. Setting this overrides
           the `nloop` parameter. The default is ``None``.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : bool, optional
           Should the step size be logarithmically spaced? The
           default (``False``) is to use a linear grid.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
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

        Vary the ``gamma`` parameter of the ``p1`` model component for
        all data sets with a source expression.

        >>> int_proj(p1.gamma)

        Use only the data in data set 1:

        >>> int_proj(p1.gamma, id=1)

        Use two data sets ('obs1' and 'obs2'):

        >>> int_proj(clus.kt, id='obs1', otherids=['obs2'])

        Vary the ``bgnd.c0`` parameter between 1e-4 and 2e-4,
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

        plotobj = self.get_int_proj(par=par, id=id, otherids=otherids,
                                    recalc=not replot, fast=fast,
                                    min=min, max=max, nloop=nloop,
                                    delv=delv, fac=fac, log=log, numcores=numcores)
        self._plot(plotobj, overplot=overplot)

    # DOC-NOTE: I am not convinced I have fac described correctly
    # DOC-NOTE: same synopsis as int_proj
    def int_unc(self, par,
                id: Optional[IdType] = None,
                otherids: Optional[Sequence[IdType]] = None,
                replot=False, min=None,
                max=None, nloop=20, delv=None, fac=1, log=False,
                numcores=None,
                overplot=False) -> None:
        """Calculate and plot the fit statistic versus fit parameter value.

        Create a confidence plot of the fit statistic as a function of
        parameter value. Dashed lines are added to indicate the
        current statistic value and the parameter value at this
        point. The parameter value is varied over a grid of points and
        the statistic evaluated while holding the other parameters
        fixed. It is expected that this is run after a successful fit,
        so that the parameter values identify the best-fit location.

        .. versionchanged:: 4.16.1
           The log parameter can now be set to `True`.

        Parameters
        ----------
        par
           The parameter to plot.
        id : int, str, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : sequence of int or str, or None, optional
           Other data sets to use in the calculation.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `int_proj`. The default is ``False``.
        min : number, optional
           The minimum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : number, optional
           The maximum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : int, optional
           The number of steps to use. This is used when `delv` is set
           to ``None``.
        delv : number, optional
           The step size for the parameter. Setting this overrides
           the `nloop` parameter. The default is ``None``.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : bool, optional
           Should the step size be logarithmically spaced? The
           default (``False``) is to use a linear grid.
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
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

        Vary the ``gamma`` parameter of the ``p1`` model component for
        all data sets with a source expression.

        >>> int_unc(p1.gamma)

        Use only the data in data set 1:

        >>> int_unc(p1.gamma, id=1)

        Use two data sets ('obs1' and 'obs2'):

        >>> int_unc(clus.kt, id='obs1', otherids=['obs2'])

        Vary the ``bgnd.c0`` parameter between 1e-4 and 2e-4,
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

        plotobj = self.get_int_unc(par=par, id=id, otherids=otherids,
                                   recalc=not replot, min=min, max=max,
                                   nloop=nloop, delv=delv, fac=fac,
                                   log=log, numcores=numcores)
        self._plot(plotobj, overplot=overplot)

    # DOC-TODO: how is sigma converted into delta_stat
    def reg_proj(self, par0, par1,
                 id: Optional[IdType] = None,
                 otherids: Optional[Sequence[IdType]] = None,
                 replot=False,
                 fast=True, min=None, max=None, nloop=(10, 10), delv=None,
                 fac=4, log=(False, False), sigma=(1, 2, 3), levels=None,
                 numcores=None,
                 overplot=False) -> None:
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
        id : int, str, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : sequence of int or str, or None, optional
           Other data sets to use in the calculation.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `int_proj`. The default is ``False``.
        fast : bool, optional
           If ``True`` then the fit optimization used may be changed from
           the current setting (only for the error analysis) to use
           a faster optimization method. The default is ``False``.
        min : pair of numbers, optional
           The minimum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : pair of number, optional
           The maximum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : pair of int, optional
           The number of steps to use. This is used when `delv` is set
           to ``None``.
        delv : pair of number, optional
           The step size for the parameter. Setting this overrides
           the `nloop` parameter. The default is ``None``.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : pair of bool, optional
           Should the step size be logarithmically spaced? The
           default (``False``) is to use a linear grid.
        sigma : sequence of number, optional
           The levels at which to draw the contours. The units are the
           change in significance relative to the starting value,
           in units of sigma.
        levels : sequence of number, optional
           The numeric values at which to draw the contours. This
           overrides the `sigma` parameter, if set (the default is
           ``None``).
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
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

        Vary the ``xpos`` and ``ypos`` parameters of the ``gsrc`` model
        component for all data sets with a source expression.

        >>> reg_proj(gsrc.xpos, gsrc.ypos)

        Use only the data in data set 1:

        >>> reg_proj(gsrc.xpos, gsrc.ypos, id=1)

        Only display the one- and three-sigma contours:

        >>> reg_proj(gsrc.xpos, gsrc.ypos, sigma=(1, 3))

        Display contours at values of 5, 10, and 20 more than the
        statistic value of the source model for data set 1:

        >>> s0 = calc_stat(id=1)
        >>> lvls = s0 + np.asarray([5, 10, 20])
        >>> reg_proj(gsrc.xpos, gsrc.ypos, levels=lvls, id=1)

        Increase the limits of the plot and the number of steps along
        each axis:

        >>> reg_proj(gsrc.xpos, gsrc.ypos, id=1, fac=6, nloop=(41, 41))

        Compare the ``ampl`` parameters of the ``g`` and ``b`` model
        components, for data sets 'core' and 'jet', over the given
        ranges:

        >>> reg_proj(g.ampl, b.ampl, min=(0, 1e-4), max=(0.2, 5e-4),
        ...          nloop=(51, 51), id='core', otherids=['jet'])

        """

        plotobj = self.get_reg_proj(par0, par1, id=id, otherids=otherids,
                                    recalc=not replot, fast=fast,
                                    min=min, max=max, nloop=nloop,
                                    delv=delv, fac=fac, log=log,
                                    sigma=sigma, levels=levels,
                                    numcores=numcores)
        self._contour(plotobj, overcontour=overplot)

    # DOC-TODO: how is sigma converted into delta_stat
    def reg_unc(self, par0, par1,
                id: Optional[IdType] = None,
                otherids: Optional[Sequence[IdType]] = None,
                replot=False,
                min=None, max=None, nloop=(10, 10), delv=None, fac=4,
                log=(False, False), sigma=(1, 2, 3), levels=None,
                numcores=None,
                overplot=False) -> None:
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
        id : int, str, or None, optional
           The data set that provides the data. If not given then
           all data sets with an associated model are used simultaneously.
        otherids : sequence of int or str, or None, optional
           Other data sets to use in the calculation.
        replot : bool, optional
           Set to ``True`` to use the values calculated by the last
           call to `int_proj`. The default is ``False``.
        min : pair of numbers, optional
           The minimum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        max : pair of number, optional
           The maximum parameter value for the calculation. The
           default value of ``None`` means that the limit is calculated
           from the covariance, using the `fac` value.
        nloop : pair of int, optional
           The number of steps to use. This is used when `delv` is set
           to ``None``.
        delv : pair of number, optional
           The step size for the parameter. Setting this overrides
           the `nloop` parameter. The default is ``None``.
        fac : number, optional
           When `min` or `max` is not given, multiply the covariance
           of the parameter by this value to calculate the limit
           (which is then added or subtracted to the parameter value,
           as required).
        log : pair of bool, optional
           Should the step size be logarithmically spaced? The
           default (``False``) is to use a linear grid.
        sigma : sequence of number, optional
           The levels at which to draw the contours. The units are the
           change in significance relative to the starting value,
           in units of sigma.
        levels : sequence of number, optional
           The numeric values at which to draw the contours. This
           overrides the `sigma` parameter, if set (the default is
           ``None``).
        numcores : optional
           The number of CPU cores to use. The default is to use all
           the cores on the machine.
        overplot : bool, optional
           If ``True`` then add the data to an existing plot, otherwise
           create a new plot. The default is ``False``.

        See Also
        --------
        conf : Estimate parameter confidence intervals using the confidence method.
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

        Vary the ``xpos`` and ``ypos`` parameters of the ``gsrc`` model
        component for all data sets with a source expression.

        >>> reg_unc(gsrc.xpos, gsrc.ypos)

        Use only the data in data set 1:

        >>> reg_unc(gsrc.xpos, gsrc.ypos, id=1)

        Only display the one- and three-sigma contours:

        >>> reg_unc(gsrc.xpos, gsrc.ypos, sigma=(1, 3))

        Display contours at values of 5, 10, and 20 more than the
        statistic value of the source model for data set 1:

        >>> s0 = calc_stat(id=1)
        >>> lvls = s0 + np.asarray([5, 10, 20])
        >>> reg_unc(gsrc.xpos, gsrc.ypos, levels=lvls, id=1)

        Increase the limits of the plot and the number of steps along
        each axis:

        >>> reg_unc(gsrc.xpos, gsrc.ypos, id=1, fac=6, nloop=(41, 41))

        Compare the ``ampl`` parameters of the ``g`` and ``b`` model
        components, for data sets 'core' and 'jet', over the given
        ranges:

        >>> reg_unc(g.ampl, b.ampl, min=(0, 1e-4), max=(0.2, 5e-4),
        ...         nloop=(51, 51), id='core', otherids=['jet'])

        Overplot the results on the `reg_proj` plot:

        >>> reg_proj(s1.c0, s2.xpos)
        >>> reg_unc(s1.c0, s2.xpos, overplot=True)

        """

        plotobj = self.get_reg_unc(par0, par1, id=id, otherids=otherids,
                                   recalc=not replot, min=min, max=max, nloop=nloop,
                                   delv=delv, fac=fac, log=log,
                                   sigma=sigma, levels=levels,
                                   numcores=numcores)
        self._contour(plotobj, overcontour=overplot)

    # Aliases
    # interval_projection = int_proj
    # interval_uncertainty = int_unc
    # region_projection = reg_proj
    # region_uncertainty = reg_unc

    ###########################################################################
    # Basic imaging
    ###########################################################################

    #
    # Image object access
    #

    def get_data_image(self,
                       id: Optional[IdType] = None):
        """Return the data used by image_data.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        data_img : a `sherpa.image.DataImage` instance
           The ``y`` attribute contains the ratio values as a 2D NumPy
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
        imageobj = self._image_types['data']
        data = self.get_data(id)
        imageobj.prepare_image(data)
        return imageobj

    def get_model_image(self,
                        id: Optional[IdType] = None):
        """Return the data used by image_model.

        Evaluate the source expression for the image pixels -
        including any PSF convolution defined by `set_psf` - and
        return the results.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        src_img : a `sherpa.image.ModelImage` instance
           The ``y`` attribute contains the source model values as a 2D
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
        imageobj = self._image_types['model']
        data = self.get_data(id)
        model = self.get_model(id)
        imageobj.prepare_image(data, model)
        return imageobj


    # DOC-TODO: it looks like get_source_image doesn't raise DataErr with
    # a non-2D data set
    def get_source_image(self,
                         id: Optional[IdType] = None):
        """Return the data used by image_source.

        Evaluate the source expression for the image pixels - without
        any PSF convolution - and return the results.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        src_img : a `sherpa.image.SourceImage` instance
           The ``y`` attribute contains the source model values as a 2D
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
        imageobj = self._image_types['source']
        data = self.get_data(id)
        source = self.get_source(id)
        imageobj.prepare_image(data, source)
        return imageobj

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
        cpt_img : a `sherpa.image.ComponentModelImage` instance
           The ``y`` attribute contains the component model values as a
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
        model = self._check_model(model)

        # Ensure the convolution models get applied
        is_source = self._get_model_status(id)[1]
        model = self._add_convolution_models(id, self.get_data(id),
                                             model, is_source)

        imageobj = self._image_types['model_component']
        data = self.get_data(id)
        imageobj.prepare_image(data, model)
        return imageobj

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
        cpt_img : a `sherpa.image.ComponentSourceImage` instance
           The ``y`` attribute contains the component model values as a
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

        Get the 'bgnd' model pixel values for data set 2:

        >>> sinfo = get_source_component_image(2, bgnd)

        """
        if model is None:
            id, model = model, id
        model = self._check_model(model)

        imageobj = self._image_types['source_component']
        data = self.get_data(id)
        imageobj.prepare_image(data, model)
        return imageobj

    def get_ratio_image(self,
                        id: Optional[IdType] = None):
        """Return the data used by image_ratio.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        ratio_img : a `sherpa.image.RatioImage` instance
           The ``y`` attribute contains the ratio values as a 2D NumPy
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
        imageobj = self._image_types['ratio']
        data = self.get_data(id)
        model = self.get_model(id)
        imageobj.prepare_image(data, model)
        return imageobj

    def get_resid_image(self,
                        id: Optional[IdType] = None):
        """Return the data used by image_resid.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        resid_img : a `sherpa.image.ResidImage` instance
           The ``y`` attribute contains the residual values as a 2D
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
        imageobj = self._image_types['resid']
        data = self.get_data(id)
        model = self.get_model(id)
        imageobj.prepare_image(data, model)
        return imageobj

    def get_psf_image(self,
                      id: Optional[IdType] = None):
        """Return the data used by image_psf.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        psf_data : a `sherpa.image.PSFImage` instance

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
        imageobj = self._image_types['psf']
        psf = self.get_psf(id)
        data = self.get_data(id)
        imageobj.prepare_image(psf, data)
        return imageobj

    def get_kernel_image(self,
                         id: Optional[IdType] = None):
        """Return the data used by image_kernel.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default identifier is
           used, as returned by `get_default_id`.

        Returns
        -------
        psf_data : a `sherpa.image.PSFKernelImage` instance

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
        imageobj = self._image_types['kernel']
        psf = self.get_psf(id)
        data = self.get_data(id)
        imageobj.prepare_image(psf, data)
        return imageobj

    #
    # Images
    #

    def image_data(self,
                   id: Optional[IdType] = None,
                   newframe=False,
                   tile=False) -> None:
        """Display a data set in the image viewer.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If ``False``, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If ``False``, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_data_image : Return the data used by image_data.
        image_close : Close the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_open : Open the image viewer.
        image_source : Display the model for a data set in the image viewer.

        Notes
        -----
        Image visualization is optional, and provided by the
        `DS9 application <https://ds9.si.edu/>`_.

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
        imageobj = self.get_data_image(id)
        imageobj.image(newframe=newframe, tile=tile)

    def image_model(self,
                    id: Optional[IdType] = None,
                    newframe=False,
                    tile=False) -> None:
        """Display the model for a data set in the image viewer.

        This function evaluates and displays the model expression for
        a data set, including any instrument response (e.g. PSF or ARF
        and RMF) whether created automatically or with
        `set_full_model`.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If ``False``, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If ``False``, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_model_image : Return the data used by image_model.
        image_close : Close the image viewer.
        image_fit : Display the data, model, and residuals for a data set in the image viewer.
        image_model_component : Display a component of the model in the image viewer.
        image_open : Open the image viewer.
        image_source : Display the model for a data set in the image viewer.
        image_source_component : Display a component of the source expression in the image viewer.

        Notes
        -----
        Image visualization is optional, and provided by the
        `DS9 application <https://ds9.si.edu/>`_.

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
        imageobj = self.get_model_image(id)
        imageobj.image(newframe=newframe, tile=tile)

    def image_source_component(self, id, model=None, newframe=False,
                               tile=False) -> None:
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
           Create a new frame for the data? If ``False``, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If ``False``, the default, then
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
        `DS9 application <https://ds9.si.edu/>`_.

        Examples
        --------

        Display the full source model and then just the 'gsrc'
        component for the default data set:

        >>> image_source()
        >>> image_source_component(gsrc)

        Display the 'clus' and 'bgnd' components of the model for the
        'img' data set side by side:

        >>> image_source_component('img', 'clus')
        >>> image_source_component('img', 'bgnd', newframe=True,
        ...                        tile=True)

        """
        imageobj = self.get_source_component_image(id, model)
        imageobj.image(newframe=newframe, tile=tile)

    def image_model_component(self, id, model=None, newframe=False,
                              tile=False) -> None:
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
           Create a new frame for the data? If ``False``, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If ``False``, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist or a source expression has
           not been set.

        See Also
        --------
        get_model_component_image : Return the data used by image_model_component.
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
        `DS9 application <https://ds9.si.edu/>`_.

        Examples
        --------

        Display the full source model and then just the 'gsrc'
        component for the default data set:

        >>> image_model()
        >>> image_model_component(gsrc)

        Display the 'clus' component of the model for the 'img' data
        set side by side without the with any instrument response
        (such as convolution with a PSF model):

        >>> image_source_component('img', 'clus')
        >>> image_model_component('img', 'clus', newframe=True,
        ...                       tile=True)

        """
        imageobj = self.get_model_component_image(id, model)
        imageobj.image(newframe=newframe, tile=tile)

    def image_source(self,
                     id: Optional[IdType] = None,
                     newframe=False,
                     tile=False) -> None:
        """Display the source expression for a data set in the image viewer.

        This function evaluates and displays the model expression for
        a data set, without any instrument response.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If ``False``, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If ``False``, the default, then
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
        `DS9 application <https://ds9.si.edu/>`_.

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
        imageobj = self.get_source_image(id)
        imageobj.image(newframe=newframe, tile=tile)

    # DOC-TODO: does newframe make sense here?
    def image_fit(self,
                  id: Optional[IdType] = None,
                  newframe=True, tile=True,
                  deleteframes=True) -> None:
        """Display the data, model, and residuals for a data set in the image viewer.

        This function displays the data, model (including any
        instrument response), and the residuals (data - model), for a
        data set.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If ``False``, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If ``False``, the default, then
           only a single frame is displayed.
        deleteframes : bool, optional
           Should existing frames be deleted? The default is ``True``.

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
        `DS9 application <https://ds9.si.edu/>`_.

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
        data = self.get_data_image(id)
        model = self.get_model_image(id)
        resid = self.get_resid_image(id)
        deleteframes = sherpa.utils.bool_cast(deleteframes)
        if deleteframes is True:
            sherpa.image.Image.open()
            sherpa.image.Image.delete_frames()
        data.image(None, False, tile)
        model.image(None, newframe, tile)
        resid.image(None, newframe, tile)

    def image_resid(self,
                    id: Optional[IdType] = None,
                    newframe=False,
                    tile=False) -> None:
        """Display the residuals (data - model) for a data set in the image viewer.

        This function displays the residuals (data - model) for a data
        set.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If ``False``, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If ``False``, the default, then
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
        `DS9 application <https://ds9.si.edu/>`_.

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
        imageobj = self.get_resid_image(id)
        imageobj.image(newframe=newframe, tile=tile)

    def image_ratio(self,
                    id: Optional[IdType] = None,
                    newframe=False,
                    tile=False) -> None:
        """Display the ratio (data/model) for a data set in the image viewer.

        This function displays the ratio data/model for a data
        set.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If ``False``, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If ``False``, the default, then
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
        `DS9 application <https://ds9.si.edu/>`_.

        Examples
        --------

        Display the ratio (data/model) for the default data set.

        >>> image_ratio()

        """
        imageobj = self.get_ratio_image(id)
        imageobj.image(newframe=newframe, tile=tile)

    # DOC-TODO: what gets displayed when there is no PSF?
    def image_psf(self,
                  id: Optional[IdType] = None,
                  newframe=False,
                  tile=False) -> None:
        """Display the 2D PSF model for a data set in the image viewer.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If ``False``, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If ``False``, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_psf_image : Return the data used by image_psf.
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
        `DS9 application <https://ds9.si.edu/>`_.

        Examples
        --------

        >>> image_psf()

        >>> image_psf(2)

        """
        imageobj = self.get_psf_image(id)
        imageobj.image(newframe=newframe, tile=tile)

    # DOC-TODO: what gets displayed when there is no PSF?
    # DOC-TODO: where to point to for PSF/kernel discussion/description
    # (as it appears in a number of places)?
    def image_kernel(self,
                     id: Optional[IdType] = None,
                     newframe=False,
                     tile=False) -> None:
        """Display the 2D kernel for a data set in the image viewer.

        The image viewer is automatically started if it is not
        already open.

        Parameters
        ----------
        id : int, str, or None, optional
           The data set. If not given then the default
           identifier is used, as returned by `get_default_id`.
        newframe : bool, optional
           Create a new frame for the data? If ``False``, the default,
           then the data will be displayed in the current frame.
        tile : bool, optional
           Should the frames be tiles? If ``False``, the default, then
           only a single frame is displayed.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_kernel_image : Return the data used by image_kernel.
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
        `DS9 application <https://ds9.si.edu/>`_.

        Examples
        --------

        >>> image_kernel()

        >>> image_kernel(2)

        """
        imageobj = self.get_kernel_image(id)
        imageobj.image(newframe=newframe, tile=tile)

    # Manage these functions (open, close, delete frames, regions, XPA)
    # through unbound functions of the Image class--always talking to
    # the same instance of the Image backend, so OK for now

    def image_deleteframes(self) -> None:
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

    def image_open(self) -> None:
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
        `DS9 application <https://ds9.si.edu/>`_.

        Examples
        --------

        >>> image_open()

        """
        sherpa.image.Image.open()

    def image_close(self) -> None:
        """Close the image viewer.

        Close the image viewer created by a previous call to one
        of the ``image_xxx`` functions.

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

    # DOC-TODO: what is the "default" coordinate system
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

    # DOC-TODO: what is the "default" coordinate system
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

    # DOC-TODO: check the ds9 link when it is working
    # DOC-TODO: is there a link of ds9 commands we can link to?
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
        The XPA access point of the ds9 image viewer lets commands and
        queries to be sent to the viewer.

        Examples
        --------

        Return the current zoom setting of the active frame:

        >>> image_xpaget('zoom')
        '1\\n'

        """
        return sherpa.image.Image.xpaget(arg)

    # DOC-TODO: check the ds9 link when it is working
    # DOC-TODO: is there a link of ds9 commands we can link to?
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
        The XPA access point of the ds9 image viewer lets commands and
        queries to be sent to the viewer.


        Examples
        --------

        Change the zoom setting of the active frame:

        >>> image_xpaset('zoom 4')

        Overlay the coordinate grid on the current frame:

        >>> image_xpaset('grid yes')

        Add the region file 'src.reg' to the display:

        >>> image_xpaset('regions src.reg')

        Create a png version of the image being displayed:

        >>> image_xpaset('saveimage png /tmp/img.png')

        """
        return sherpa.image.Image.xpaset(arg, data)
