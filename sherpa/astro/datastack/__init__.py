# Copyright (c) 2010, 2014, 2015, 2020, 2021 Smithsonian Astrophysical Observatory
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

Acknowledgements
----------------

The Datastack code was developed by Thomas Aldcroft and originally
provided as an external package for Sherpa, available at
http://cxc.harvard.edu/contrib/datastack/. The code was added to
Sherpa in version 4.7.1.

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

import logging
import sys
import types

from sherpa.astro import ui
from sherpa.utils.logging import config_logger
from sherpa.utils import public
from sherpa.astro.datastack.ds import DataStack

from .utils import set_template_id

logger = config_logger(__name__)

__all__ = ['set_template_id', 'DataStack']


@public
def set_stack_verbosity(level):
    """Change the logging level.

    Informational messages from the datastack module are displayed using
    the Python logging system. This routine determines the severity
    of the messages that are displayed.

    Parameters
    ----------
    level
       The logging level to use (e.g. `logging.INFO` or
       `logging.WARNING`).

    See Also
    --------
    set_stack_verbose

    Examples
    --------
    >>> set_stack_verbosity(logging.ERROR)

    """
    logger.setLevel(level)


@public
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


_always_wrapped = ('load_pha', 'load_arrays', 'load_ascii', 'load_data',
                   'load_bkg')


# Use this and subsequent loop to wrap every function in sherpa.astro.ui
# with a datastack version
def _sherpa_ui_wrap(func):
    def wrap(*args, **kwargs):
        wrapfunc = func
        if args:
            if isinstance(args[0], DataStack):
                datastack, args = args[0], args[1:]

            # If the first argument is a list and it's either empty or
            # made of non-iterables, then it's a datastack definition.
            # If the list contains iterable it must be arrays for load_arrays.
            elif isinstance(args[0], list) and \
                    not (len(args[0]) > 0 and hasattr(args[0][0], '__iter__')):
                datastack = DATASTACK[args[0]] if args[0] else DATASTACK
                args = args[1:]
            else:
                if func.__name__ in _always_wrapped:
                    # some (all?) load_* functions must always be wrapped
                    # for file stack syntax check and for ensuring dataset
                    # id consistency.
                    datastack = DATASTACK
                else:
                    # No stack specifier so use native sherpa func
                    return func(*args, **kwargs)

            try:
                wrapfunc = getattr(datastack, func.__name__)
            except AttributeError:
                raise AttributeError(
                    '{0} is not a stack-enabled function.'.format(func.__name__))

        return wrapfunc(*args, **kwargs)

    wrap.__name__ = func.__name__
    if not hasattr(DataStack, func.__name__):
        doc = func.__doc__
    else:
        doc = getattr(DataStack, func.__name__).__doc__
    wrap.__doc__ = doc
    return wrap


def _datastack_wrap(func):
    def wrap(*args, **kwargs):
        if not args:
            args = ([],) + args

        if isinstance(args[0], DataStack):
            datastack, args = args[0], args[1:]
        elif isinstance(args[0], list) and \
                not (len(args[0]) > 0 and hasattr(args[0][0], '__iter__')):
            datastack = DATASTACK[args[0]] if args[0] else DATASTACK
            args = args[1:]
        else:
            datastack = DATASTACK

        return getattr(datastack, func.__name__)(*args, **kwargs)

    wrap.__name__ = func.__name__
    if not hasattr(DataStack, func.__name__):
        doc = func.__doc__
    else:
        doc = getattr(DataStack, func.__name__).__doc__
    wrap.__doc__ = doc
    return wrap


# The default datastack
DATASTACK = DataStack()

# Wrap all sherpa UI funcs and a few DataStack methods for the
# command-line interface.
_module = sys.modules[__name__]
for attr in dir(ui):
    func = getattr(ui, attr)
    if isinstance(func, types.FunctionType):
        setattr(_module, attr, public(_sherpa_ui_wrap(func)))

for funcname in ['clear_stack', 'show_stack', 'get_stack_ids',
                 'query', 'query_by_header_keyword', 'query_by_obsid']:
    setattr(_module, funcname, public(
        _datastack_wrap(getattr(DataStack, funcname))))


@public
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
    logger.warning("clean() will invalidate any existing DataStack " +
                   "instances by removing all the datasets from the " +
                   "Sherpa session")
