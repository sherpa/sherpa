#
#  Copyright (C) 2022
#      MIT
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
from sherpa.plot.backend_utils import (translate_args,
                                       add_kwargs_to_doc,
                                       get_keyword_defaults)

class ExampleClass():
    def __init__(self, translate_dict={}):
        self.translate_dict = translate_dict

    @translate_args
    def func(self, a, b, c='qwer'):
        return a, b, c


def test_translate_dict():
    """Test that a translatedict filled with dicts works"""
    t = ExampleClass({'a': {'val1': 'val2'},
                      'c': {None: 5, 5: 6}})
    # No translation of value not in dict
    assert (3, 4, 2.) == t.func(3, 4, c=2.)
    # translate value for just one of the values
    assert ('val1', 4, 2.) == t.func('val1', 4, 2.)
    # translate all
    assert ('val2', 7, 6) == t.func(a='val1', b=7, c=5)
    # translate does only translate keyword arguments
    assert ('val1', 7, 6) == t.func('val1', 7, c=5)


def test_translate_func():
    """Test that a translatedict filled with functions works"""
    t = ExampleClass({'a': {'val1': 'val2'},
                      'c': lambda x: x * 2})
    assert ('a', 'b', 16) == t.func('a', 'b', c=8)


def test_keyword_defaults():
    """Check that we get a dictionary of the default values defined in a function"""
    def func(a, b=5, c=None):
        pass
    assert get_keyword_defaults(func) == {'b': 5, 'c': None}


def test_keyword_defaults_method():
    """Repeat the test above for a method"""
    class MyClass():
        def func(self, a, b=5, c=None):
            pass

        def keyword(self):
            return get_keyword_defaults(self.func)

    assert get_keyword_defaults(MyClass.func) == {'b': 5, 'c': None}
    myinstance = MyClass()
    assert MyClass.keyword(MyClass) == {'b': 5, 'c': None}


param_doc = {'c' : ['int', 'thickness of line'],
             'title' : ['string', '''Title of figure (only use if `overplot=False`)'''],
             'color': ['string or number', 'any matplotlib color with a really long text attached to it that will not fit in the one line of text in the docstring']
             }

class A():
    @add_kwargs_to_doc(param_doc)
    def func(a, *, title=None, color='None', **kwargs):
        '''Method that does nothing

        more text here

        Parameters
        ----------
        a : int
            Our stuff
        {kwargs}

        Returns
        -------
        something
        '''
        pass


def test_modify_doctring():
    '''Check that kwarg descriptions are properly inserted into the docstring.'''
    a = A()
    expected = '''Method that does nothing

        more text here

        Parameters
        ----------
        a : int
            Our stuff
        title : string, default=None
            Title of figure (only use if `overplot=False`)
        color : string or number, default=None
            any matplotlib color with a really long text attached to it that will not fit in the one line of text in the docstring
        kwargs : dict, optional
            All other keyword parameters are passed to the plotting library.

        Returns
        -------
        something
        '''
    assert a.func.__doc__ == expected
 