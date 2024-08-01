#
#  Copyright (C) 2013, 2016, 2018, 2021  Smithsonian Astrophysical Observatory
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

from sherpa.utils.err import SherpaErr
from sherpa.utils.err import ModelErr


class OldSherpaErr(Exception):
    "Old class for all Sherpa exceptions"

    def __init__(self, dict, key, *args):

        if key in dict:
            errmsg = dict[key] % args
        else:
            errmsg = "unknown key '%s'" % key
        Exception.__init__(self, errmsg)


dict = {'simple': 'simple message', 'arg': 'argument: %s'}


def test1():
    """verify that a correct call of the new constructor has the same result of the old one"""
    err = SherpaErr(dict, 'simple')
    old_err = OldSherpaErr(dict, 'simple')
    assert str(err) == str(old_err)


def test2():
    """same as before, but with string placeholders"""
    err = SherpaErr(dict, 'arg', 'foo')
    assert 'argument: foo' == str(err)


def test3():
    """verify that a call without a key results in a generic message being produced"""
    err = SherpaErr(dict)
    assert 'Generic Error' == str(err)


def test4():
    """verify the user's expected behavior, i.e. a string is provided as error message"""
    err = SherpaErr(dict, 'My Error')
    assert 'My Error' == str(err)


def test5():
    """verify the user provided example, which exercises a derived class"""
    err = ModelErr("Unable to frobnicate model %s" % 'modelname')
    assert 'Unable to frobnicate model modelname' == str(err)
