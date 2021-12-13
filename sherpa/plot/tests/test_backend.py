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
'''This module contains tests for the plotting backend infrastructure.

Those tests don't create a plot, but check that backends can be selected (or not),
that backends perform checks on inputs etc.
'''
import logging
from tkinter import Place
import pytest

from sherpa.plot.backends import (IndepOnlyBackend, BaseBackend, PLOT_BACKENDS,
                                  backend_indep_colors)
from sherpa.plot import set_backend
from sherpa import plot


def test_IndepOnlyBackend_raises_for_values():
    '''The IndepOnlyBackend is supposed to raise exceptions when
    called with parameters that are not in the
    backend-independent list.'''
    back = IndepOnlyBackend()
    with pytest.raises(ValueError, match='but got xxx'):
        back.plot(1, 2, linestyle='xxx')

    with pytest.raises(TypeError,
        match='contour got keyword argument color, which is not part of the named keyword arguments'):
        back.contour(1, 2, 3, color='yyy')


def test_IndepOnlyBackend_raises_for_arguments():
    '''The IndepOnlyBackend is supposed to raise exceptions when
    called with keywords that are not in the
    backend-independent list.'''
    back = IndepOnlyBackend()
    with pytest.raises(TypeError,
                       match='plot got keyword argument notthis'):
        back.plot(1, 2, notthis=5)


def test_warning_backend_unkown(caplog):
    '''Check that a warning appears if a backend is not registered'''
    with caplog.at_level(logging.WARNING):
        set_backend('qweraef')
    assert 'qweraef is not a known plotting backend' in caplog.text


def test_backend_is_set():
    '''check that set_backend sets backend instances'''
    # We don't want this test to mess up the backends for other tests.
    old_backend = plot.backend
    try:
        # Don't know what old_backend is, since that depends on how the tests
        # are run, but if we set_backend twice, we know for sure it's beem
        # changed.
        set_backend('BaseBackend')
        assert isinstance(plot.backend, BaseBackend)
        set_backend('IndepOnlyBackend')
        assert isinstance(plot.backend, IndepOnlyBackend)
    finally:
        plot.backend = old_backend


def test_backend_registry(caplog):
    '''Check that a warning is issued for dublicate plotting backend names'''

    class StupidName(BaseBackend):
        pass
    assert 'StupidName' in PLOT_BACKENDS.keys()

    class StupidName2(BaseBackend):
        name = 'othername'
    assert 'othername' in PLOT_BACKENDS.keys()

    with caplog.at_level(logging.WARNING):
        class OtherStupidName(BaseBackend):
            name = "StupidName"
    print(PLOT_BACKENDS)
    assert 'StupidName is already a registered name' in caplog.text
    assert PLOT_BACKENDS['StupidName'] == StupidName
    del PLOT_BACKENDS['StupidName']
    del PLOT_BACKENDS['othername']


def test_colorlist():
    '''Check that we can generate a list of n colors'''
    back = BaseBackend()
    clist = back.colorlist(25)
    assert len(clist) == 25
    assert all([c in backend_indep_colors for c in clist])
