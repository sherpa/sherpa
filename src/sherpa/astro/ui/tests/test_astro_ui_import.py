#
#  Copyright (C) 2022
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

"""While frowned-upon in polite Python society, many Sherpa users
will be using the

    from sherpa.astro.ui import *

form, either drectly or indirectly (when using the sherpa application
in CIAO). So this adds some basic tests of this behavior, as there are
a few "interesting" things to review.

"""

import sys

import pytest

from sherpa.astro.ui import *
from sherpa.astro.ui.utils import Session
import sherpa.models.basic
from sherpa.ui.utils import ModelWrapper


# Select models from sherpa.models.basic and sherpa.astro.models
#
MODEL_NAMES = ["const1d", "gauss2d", "lorentz1d"]


def test_secret_session_is_created():
    """This is a regression test.

    The sherpa.ui module does not export _session but astro does.
    """

    assert isinstance(_session, Session)


@pytest.mark.parametrize("model", MODEL_NAMES)
def test_model_is_known(model):
    """Check we have loaded some models"""

    assert model in list_models()


@pytest.mark.parametrize("model", MODEL_NAMES)
def test_model_is_defined(model):
    """Check we have loaded some models"""

    assert model in globals()

    sym = globals()[model]
    assert isinstance(sym, ModelWrapper)


@pytest.mark.parametrize("as_string", [True, False])
def test_model_identifiers_set_globally(as_string):
    """Check we create a global symbol for the models.

    See also the same test in
      sherpa/astro/ui/tests/test_astro_session.py
      sherpa/astro/ui/tests/test_astro_ui_unit.py

    """

    # The "global" symbol table depends on what has been run before. We
    # could try and make sure that we are "clean", but this makes checking
    # what this test is doing hard to do, so we remove the symbols just
    # in case.
    #
    for name in ["mdl1", "mdl2"]:
        try:
            del sys.modules["__main__"].__dict__[name]
        except KeyError:
            pass

    dataspace1d(1, 10, 1)

    for store in [globals(), locals(), sys.modules["__main__"].__dict__]:
        assert "mdl1" not in store
        assert "mdl2" not in store

    if as_string:
        set_source("const1d.mdl1 + gauss1d.mdl2")
    else:
        set_source(const1d.mdl1 + gauss1d.mdl2)

    for store in [globals(), locals()]:
        assert "mdl1" not in store
        assert "mdl2" not in store

    assert "mdl1" in sys.modules["__main__"].__dict__
    assert "mdl2" in sys.modules["__main__"].__dict__

    assert isinstance(sys.modules["__main__"].__dict__["mdl1"],
                      sherpa.models.basic.Const1D)
    assert isinstance(sys.modules["__main__"].__dict__["mdl2"],
                      sherpa.models.basic.Gauss1D)

    clean()


def test_delete_model_removes_global_identifier():
    """Check we create a global symbol for the models."""

    create_model_component("const1d", "mdl1")
    create_model_component("gauss1d", "mdl2")

    assert "mdl1" in sys.modules["__main__"].__dict__
    assert "mdl2" in sys.modules["__main__"].__dict__

    delete_model_component("mdl1")
    delete_model_component("mdl2")

    assert "mdl1" not in sys.modules["__main__"].__dict__
    assert "mdl2" not in sys.modules["__main__"].__dict__


def test_what_happens_if_the_same_identifier_is_reused():
    """Regression test.

    It's not obvious what we want to do here.
    """

    # This could error out or set the "first" name as the winner, or
    # the "last", or something else. Test the current behavior.
    #
    combined = gauss1d.bob + const1d.bob
    assert isinstance(bob, sherpa.models.basic.Const1D)

    delete_model_component("bob")
