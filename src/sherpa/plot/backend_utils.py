#
#  Copyright (C) 2022, 2023
#  MIT
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
'''This module has utilities for backends.

They are not useful on their own, but can help to define new plotting backends.
'''

import functools
from inspect import signature
import re

__all__ = ('translate_args', 'add_kwargs_to_doc', 'get_keyword_defaults')

kwargs_indent = re.compile(r"\n(?P<spaces>\s*){kwargs}",flags=re.M)

def find_indent(doc):
    match = kwargs_indent.findall(doc)
    return len(match[0])


def translate_args(func):
    '''A decorator to translate function arguments.

    When a method decorated with this decorator is called,
    the decorator inspects the arguments that are passed into the
    function. For each argument that is found as a key in the object's
    ``translate_args`` dictionary, the value of that argument is
    translated. The items in that dictionary can be functions or
    dictionaries.

    The purpose of this decorator is to support backends that use different
    syntax for the same option, e.g. one backend might call a color
    "red", while another uses the tuple (1, 0, 0) to describe the same color.

    Example
    -------
    In this example, the input 'r' or 'b' will be translated into an rgb tuple
    before the ``plot`` function is called. Other values (e.g. ``color=(0, 0, 0)``)
    will be passed through unchanged so that the user can also make use of any other
    color specification that the backend allows.

    >>> from sherpa.plot import backend_utils
    >>> class Plotter:
    ...     translate_args = {'color': {'r': (1,0,0), 'b': (0,0,1)}}
    ...
    ...     @backend_utils.translate_args
    ...     def plot(color=None):
    ...         print('RGB color tuple is: ', color)
    '''
    @functools.wraps(func)
    def inner(self, *args, **kwargs):
        for kw, val in kwargs.items():
            if kw in self.translate_dict:
                transl = self.translate_dict[kw]

                if callable(transl):
                    # It's a function
                    kwargs[kw] = transl(val)
                else:
                    # It should be a dict
                    # Let's check if val is one of those that need to
                    # be translated
                    if val in transl:
                        kwargs[kw] = transl[val]
        return func(self, *args, **kwargs)

    return inner


def add_kwargs_to_doc(param_doc):
    '''Add documentation for keyword parameters

    The plotting functions for each backend take a large number of keyword
    arguments that are repeated over several methods. This decorator can
    be used to fill in the description for those keyword arguments into
    the method docstring to reduce repetition, keep the docstrings readable
    in the code files, and ensure consistency between different functions.

    The decorator uses string formatting and fills all keyword parameters
    where `{kwargs}` appears in the docstring. This allows all other parts
    of the docstring to be written normally.

    A typical use case is to define a dictionary of all parameters for a
    backend, e.g. describe all values that `'color'` can take and then pass
    that same dictionary to the decorator in all methods.

    The appropriate number of white space is inserted to generate the
    numpydoc format, but long lines are not wrapped. Instead, the line
    wrapping is preserved to maintain special markup (e.g. lists). If lines
    are very long, they should contain newlines in the `param_doc`.

    Parameters
    ----------
    param_doc : dict
        Keys in the dictionary are parameters names and values are tuples.
        The first element of the tuple is the parameters type, the second one
        the description.

    Examples
    --------

    >>> param_doc = {'c' : ['int', 'thickness of line'],
    ...              'title' : ['string', 'Title of figure (only use if `overplot=False`)'],
    ...              'color': ['string or number',
    ...                  'any matplotlib color with a really long text attached to it that will not fit in one line of text in the docstring']}
    >>> @add_kwargs_to_doc(param_doc)
    ... def test_func2(a, *, title=None, color='None'):
    ...     """Func that does nothing
    ...
    ...     Parameters
    ...     ----------
    ...     a : int
    ...         Our stuff
    ...     {kwargs}
    ...     """
    ...     pass
    >>> help(test_func2)
    Help on function test_func2 in module sherpa.plot.backend_utils:
    <BLANKLINE>
    test_func2(a, *, title=None, color='None')
        Func that does nothing
    <BLANKLINE>
        Parameters
        ----------
        a : int
            Our stuff
        title : string, default=None
            Title of figure (only use if `overplot=False`)
        color : string or number, default=None
            any matplotlib color with a really long text attached to it that will not fit in one line of text in the docstring

    '''
    def set_docstring(obj):
        sig = signature(obj)
        indent = find_indent(obj.__doc__)
        out = []

        for p, par in sig.parameters.items():
            if par.kind == par.KEYWORD_ONLY:
                pdoc = param_doc.get(p, ['', ''])
                out.append(' ' * indent + f'{p} : {pdoc[0]}, default={sig.parameters[p].default}')
            if par.kind == par.VAR_KEYWORD:
                pdoc = param_doc.get(p, ['', 'All other keyword parameters are passed to the plotting library.'])
                out.append(' ' * indent + f'{p} : dict, optional')
            if par.kind in (par.KEYWORD_ONLY, par.VAR_KEYWORD):
                for line in pdoc[1].split('\n'):
                    out.append(' ' * (indent + 4) + line)

        out = '\n'.join(out)
        # out[indent:] is needed because `   {kwargs}` already has spaces in front of it.
        # We don't want to add those back a second time.
        obj.__doc__ = obj.__doc__.format(kwargs=out[indent:])
        return obj
    return set_docstring


def get_keyword_defaults(func, ignore_args=['title', 'xlabel', 'ylabel',
                                            'overplot', 'overcontour',
                                            'clearwindow', 'clearaxes',
                                            'xerr', 'yerr']):
    '''Get default values for keyword arguments

    This method differs from `sherpa.utils.get_keyword_defaults`, which inspects all
    arguments, while this function looks only at keyword arguments. Also,
    in `sherpa.utils.get_keyword_defaults` arguments can be skipped by order,
    while here they are skipped by name (using the `ignore_args`) parameter.
    Thus, this function is better suited to the plotting backends, which use
    several keyword-only arguments.

    Parameters
    ----------
    func : callable
        function or method to inspect
    ignore_args : list
        Any keyword arguments with names listed here will be ignored

    Returns
    -------
    default_values : dict
        Dictionary with argument names and default values

    See also
    --------
    sherpa.utils.get_keyword_defaults
    '''
    default_values = {}
    sig = signature(func)

    for param in sig.parameters.values():
        if (param.default is not param.empty and
                param.name not in ignore_args):
            default_values[param.name] = param.default
    return default_values
