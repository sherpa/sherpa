#
#  Copyright (C) 2016, 2018, 2020 - 2024
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

# pylint: disable=missing-function-docstring
# pylint: disable=invalid-name

from io import StringIO
import logging
import os
import re
import sys

import numpy

import pytest

from sherpa.astro.data import DataPHA
from sherpa.astro.instrument import create_arf, create_delta_rmf
from sherpa.astro import io
from sherpa.astro.ui.utils import Session as AstroSession
from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt
from sherpa.instrument import ConvolutionKernel
from sherpa.io import get_ascii_data
from sherpa.models import ArithmeticModel, Const1D
import sherpa.models.basic
from sherpa.ui.utils import Session
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, DataErr, \
    IdentifierErr, IOErr, ModelErr, PlotErr, SessionErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_data, requires_fits, requires_group


@pytest.fixture
def skip_if_no_io(request):
    """For tests with a session fixture that need to be skipped if
    the I/O backend is not present ONLY for AstroSession.
    """

    if request.getfixturevalue("session") == Session:
        return

    # If the I/O backend is not the dummy backend then we assume that
    # we can return.
    #
    if io.backend.name != "dummy":
        return

    pytest.skip(reason="FITS backend required")



@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("as_string", [True, False])
def test_model_identifiers_set_globally(session, as_string):
    """Check we create a global symbol for the models.

    See also the same test in
      sherpa/astro/ui/tests/test_astro_ui_import.py
      sherpa/astro/ui/tests/test_astro_ui_unit.py

    There shouldn't be anything do in the astro class that changes the
    behavior here, but run tests on both just in case.

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

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.dataspace1d(1, 10, 1)

    for store in [globals(), locals(), sys.modules["__main__"].__dict__]:
        assert "mdl1" not in store
        assert "mdl2" not in store

    if as_string:
        s.set_source("const1d.mdl1 + gauss1d.mdl2")
    else:
        # Unlike 'from sherpa[.astro].ui import *', the models are not
        # added to the global symbol table, so this can not use
        #
        #     s.set_source(const1d.mdl1 + gauss1d.mdl2)
        #
        # so the create_model_component call is used instead.
        #
        s.create_model_component("const1d", "mdl1")
        s.create_model_component("gauss1d", "mdl2")

    for store in [globals(), locals()]:
        assert "mdl1" not in store
        assert "mdl2" not in store

    assert "mdl1" in sys.modules["__main__"].__dict__
    assert "mdl2" in sys.modules["__main__"].__dict__

    assert isinstance(sys.modules["__main__"].__dict__["mdl1"],
                      sherpa.models.basic.Const1D)
    assert isinstance(sys.modules["__main__"].__dict__["mdl2"],
                      sherpa.models.basic.Gauss1D)

    s.clean()


def test_astro_model_identifiers_set_globally():
    """Check we create a global symbol for the background/pileup models."""

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

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.dataspace1d(1, 10, 1, dstype=DataPHA)

    # Technically gauss1d is not a pileup model but it is accepted.
    #
    s.set_bkg_source("const1d.mdl1")
    s.set_pileup_model("gauss1d.mdl2")

    for store in [globals(), locals()]:
        assert "mdl1" not in store
        assert "mdl2" not in store

    assert "mdl1" in sys.modules["__main__"].__dict__
    assert "mdl2" in sys.modules["__main__"].__dict__

    assert isinstance(sys.modules["__main__"].__dict__["mdl1"],
                      sherpa.models.basic.Const1D)
    assert isinstance(sys.modules["__main__"].__dict__["mdl2"],
                      sherpa.models.basic.Gauss1D)

    s.clean()


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_does_clean_remove_model_identifiers(session):
    """Check whether clean() will remove model identifiers.

    There shouldn't be anything do in the astro class that changes the
    behavior here, but run tests on both just in case.

    """

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.dataspace1d(1, 10, 1)

    s.create_model_component("const1d", "mdl1")
    s.create_model_component("gauss1d", "mdl2")

    s.clean()

    for store in [globals(), locals(), sys.modules["__main__"].__dict__]:
        assert "mdl1" not in store
        assert "mdl2" not in store


def test_astro_does_clean_remove_model_identifiers():
    """Check a few astro-related symbols."""

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.dataspace1d(1, 10, 1, dstype=DataPHA)

    # Technically gauss1d is not a pileup model but it is accepted.
    #
    s.set_bkg_source("const1d.mdl1")
    s.set_pileup_model("gauss1d.mdl2")

    s.clean()

    for store in [globals(), locals(), sys.modules["__main__"].__dict__]:
        assert "mdl1" not in store
        assert "mdl2" not in store


@pytest.mark.parametrize("label", ["arf", "rmf"])
def test_get_response_no_arf(label):
    """Just check we error out"""

    s = AstroSession()
    s.load_arrays(1, [1, 2, 3, 4, 5], [10, 11, 12, 13, 14], DataPHA)

    msg = f"^{label.upper()} data set 1 in PHA data set 1 has not been set$"
    with pytest.raises(IdentifierErr,
                       match=msg):
        getattr(s, f"get_{label}")()


# bug #303
def test_show_bkg_model():
    session = AstroSession()
    session.load_arrays(1, [1, 2], [1, 2])
    session.show_bkg_model()
    session.show_bkg_model('xx')
    session.show_bkg_source()
    session.show_bkg_source('xx')


# bug #303
@requires_data
@requires_fits
def test_show_bkg_model_with_bkg(make_data_path):
    session = AstroSession()
    session.load_data('foo', make_data_path('3c273.pi'))
    session.show_bkg_model()
    session.show_bkg_model('foo')


def check_nth_caplog(caplog, lname, lvl, msg, pos=-1):
    """Check the given element contains the expected data."""

    lname_got, lvl_got, msg_got = caplog.record_tuples[pos]
    assert lname_got == lname
    assert lvl_got == lvl
    assert msg_got == msg


def clc_filter(caplog, msg, astro=False, pos=-1):
    """Special case for the ignore/notice filter check"""

    reporter = "sherpa.astro.ui.utils" if astro else "sherpa.ui.utils"
    check_nth_caplog(caplog, reporter, logging.INFO, msg, pos=pos)


# Fix 476 - this should be in sherpa/ui/tests/test_session.py
@requires_group
def test_zero_division_calc_stat(caplog):
    ui = AstroSession()
    x = numpy.arange(100)
    y = numpy.zeros(100)
    ui.load_arrays(1, x, y, DataPHA)

    assert len(caplog.record_tuples) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.group_counts(1, 100)

    assert ui.get_data().grouped

    assert len(caplog.record_tuples) == 1
    clc_filter(caplog, "dataset 1: 0:99 Channel (unchanged)")

    ui.set_full_model(1, Const1D("const"))

    # in principle I wouldn't need to call calc_stat_info(), I could just
    # use _get_stat_info to reproduce the issue, However, _get_stat_info is not a public
    # method, so I want to double check that calc_stat_info does not throw an exception.
    # So, first we try to run calc_stat_info and make sure there are no exceptions.
    # Then, since calc_stat_info only logs something and doesn't return anything, we use
    # a white box approach to get the result from _get_stat_info.
    ui.calc_stat_info()
    assert ui._get_stat_info()[0].rstat is numpy.nan


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_iter_method_name_default(session):
    s = session()
    assert s.get_iter_method_name() == "none"


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("opt", ["none", "sigmarej"])
def test_set_iter_method_valid(session, opt):
    s = session()
    s.set_iter_method(opt)
    assert s.get_iter_method_name() == opt


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_set_iter_method_unknown_string(session):
    s = session()
    with pytest.raises(TypeError,
                       match="^not a method is not an iterative fitting method$"):
        s.set_iter_method("not a method")


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_set_iter_method_not_a_string(session):
    s = session()
    with pytest.raises(ArgumentTypeErr,
                       match="^'meth' must be a string$"):
        s.set_iter_method(23)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_iter_method_opt_default(session):
    s = session()
    assert s.get_iter_method_opt() == {}


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_iter_method_opt_sigmarej(session):
    s = session()
    s.set_iter_method("sigmarej")
    out = s.get_iter_method_opt()

    keys = set(out.keys())
    assert keys == set(["grow", "lrej", "hrej", "maxiters"])


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_iter_method_opt_sigmarej_lrej(session):
    s = session()
    s.set_iter_method("sigmarej")
    assert s.get_iter_method_opt("lrej") == pytest.approx(3)


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("opt,key", [("none", "lrej"), ("sigmarej", "fast")])
def test_get_iter_method_opt_unknown(session, opt, key):
    s = session()
    s.set_iter_method(opt)
    with pytest.raises(ArgumentErr,
                       match=f"^'{key}' is not a valid option for method {opt}$"):
        s.get_iter_method_opt(key)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_set_iter_method_opt_sigmarej_lrej(session):
    s = session()
    s.set_iter_method("sigmarej")
    s.set_iter_method_opt("lrej", 5)
    assert s.get_iter_method_opt("lrej") == pytest.approx(5)


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("setting", ['chisqr', 'compmodel', 'compsource', 'data',
                                     "model_component", "source_component",
                                     'delchi', 'fit', 'kernel', 'model',
                                     'psf', 'ratio', 'resid', 'source'])
def test_id_checks_session(session, setting):
    """Do some common identifiers fail?"""

    s = session()
    with pytest.raises(IdentifierErr,
                       match=f"identifier '{setting}' is a reserved word"):
        s.load_arrays(setting, [1, 2], [1, 2])


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("setting", ['cdf', 'energy', 'lr', 'photon', 'pdf', 'scatter', 'trace'])
def test_id_checks_session_unexpected(session, setting):
    """These identifiers are allowed. Should they be?"""

    s = session()
    s.load_arrays(setting, [1, 2], [1, 2])
    d = s.get_data(setting)
    assert isinstance(d, Data1D)


@pytest.mark.parametrize("session,success",
                         [(Session, True), (AstroSession, False)])
@pytest.mark.parametrize("setting", ['arf', 'bkg', 'bkgchisqr', 'bkgdelchi', 'bkgfit',
                                     'bkgmodel', 'bkgratio', 'bkgresid', 'bkgsource',
                                     "bkg_model", "bkg_source", "bkg_resid", "bkg_ratio",
                                     "bkg_delchi", "bkg_chisqr", "bkg_fit",
                                     'order',
                                     "astrocompsource", "astrocompmodel", "astrodata",
                                     "astrosource", "astromodel"])
def test_id_checks_astro_session(session, success, setting):
    """Do some common identifiers fail for astro but not default?"""

    s = session()
    if success:
        s.load_arrays(setting, [1, 2], [1, 2])
        d = s.get_data(setting)
        assert isinstance(d, Data1D)
    else:
        with pytest.raises(IdentifierErr,
                           match=f"identifier '{setting}' is a reserved word"):
            s.load_arrays(setting, [1, 2], [1, 2])


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("plottype", ["source_component", "compsource",
                                      "model_component", "compmodel"])
def test_plot_component(session, plottype):
    """Can we call plot with a "component" call.

    This is a regression test to see if compsource/model works
    with the plot call. As a check we include the "full" names
    (e.g. source_component).
    """

    s = session()
    s._add_model_types(sherpa.models.basic)
    s.load_arrays(1, [1, 2, 3], [5, 2, 3])

    mdl = s.create_model_component("const1d", "mdl")
    s.set_source(mdl)

    # All we do is check we can call the routine. We do not check it
    # has done anything sensible, but we do check you can call it
    # with and without a dataset identifier.
    #
    s.plot(plottype, 1, mdl)
    s.plot(plottype, mdl)


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("plottype", ["source_component", "compsource",
                                      "model_component", "compmodel"])
def test_plot_component_fails(session, plottype):
    """Can we call plot with a "component" call and get it to error out

    This is a regression test to see both whether we can use
    the compsource/model label, and to check how it errors out
    when no argument is given.
    """

    s = session()
    s._add_model_types(sherpa.models.basic)
    s.load_arrays(1, [1, 2, 3], [5, 2, 3])

    mdl = s.create_model_component("const1d", "mdl")
    s.set_source(mdl)

    # This is a low-level error that we don't catch and convert (at
    # least at the moment), so just check the current behavior.
    #
    with pytest.raises(TypeError,
                       match=r"_plot\(\) missing 1 required positional argument: 'id'"):
        s.plot(plottype)


@pytest.mark.parametrize("label", ["chisqr", "delchi", "fit", "model", "ratio",
                                   "resid", "source"])
def test_astro_plot_bkgxxx(label):
    """A regression test of plot("bkg<label>")

    This is astro-specific and is just a check to see if we can
    make the call, not to check what the actual plot looks like.

    See also test_astro_plot_bkg_xxx which has been separated out to
    allow different behavior over time.

    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    data = DataPHA("data", [1, 2, 3], [5, 2, 3])
    bkg = DataPHA("bkg", [1, 2, 3], [2, 1, 2])

    egrid = numpy.asarray([0.1, 0.2, 0.3, 0.4])
    arf = create_arf(egrid[:-1], egrid[1:])

    s.set_data(data)
    s.set_bkg(bkg)

    s.set_arf(arf, bkg_id=1)

    mdl = s.create_model_component("const1d", "mdl")
    s.set_bkg_source(mdl)

    s.plot(f"bkg{label}")


@pytest.mark.parametrize("label", ["chisqr", "delchi", "fit", "model", "ratio",
                                   "resid", "source"])
def test_astro_plot_bkg_xxx(label):
    """A regression test of plot("bkg_<label>")

    This is astro-specific and is just a check to see if we can
    make the call, not to check what the actual plot looks like.

    See also test_astro_plot_bkgxxx which has been separated out to
    allow different behavior over time.
    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    data = DataPHA("data", [1, 2, 3], [5, 2, 3])
    bkg = DataPHA("bkg", [1, 2, 3], [2, 1, 2])

    egrid = numpy.asarray([0.1, 0.2, 0.3, 0.4])
    arf = create_arf(egrid[:-1], egrid[1:])

    s.set_data(data)
    s.set_bkg(bkg)

    s.set_arf(arf, bkg_id=1)

    mdl = s.create_model_component("const1d", "mdl")
    s.set_bkg_source(mdl)

    s.plot(f"bkg_{label}")


def save_ascii_file(s, kwargs, idval, outfile, savefunc, syserr=False):
    """create data and save a file based on it"""

    s._add_model_types(sherpa.models.basic)

    cpt = s.create_model_component("const1d", "cpt")
    cpt.integrate = True
    cpt.c0 = 2

    # Select bin edges that lead to integer for the middle values,
    # along with the model per-bin values (2 * width).
    #
    s.load_arrays(idval, [1, 3, 5], [3, 5, 9], [3, 4, 5], Data1DInt)
    dy = [1, 2, 6]

    errfunc = s.set_syserror if syserr else s.set_staterror

    if idval is None:
        errfunc(dy)
        s.set_source(cpt)
        savefunc(str(outfile), **kwargs)
    else:
        errfunc(idval, dy)
        s.set_staterror(idval, dy)
        s.set_source(idval, cpt)
        savefunc(idval, str(outfile), **kwargs)


def get_data(path, **kwargs):
    """wrapper routine to handle path conversion and dropping last argument."""

    ans = get_ascii_data(str(path), **kwargs)
    return ans[0], ans[1]


@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_model_ascii(session, kwargs, idval, tmp_path, skip_if_no_io):
    """Can we use save_model?"""

    s = session()
    outfile = tmp_path / "model.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_model)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "MODEL"]
    assert len(data) == 2
    assert data[0] == pytest.approx([2, 4, 7])
    assert data[1] == pytest.approx([4, 4, 8])


@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_source_ascii(session, kwargs, idval, tmp_path, skip_if_no_io):
    """Can we use save_source?"""

    s = session()
    outfile = tmp_path / "source.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_source)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "SOURCE"]
    assert len(data) == 2
    assert data[0] == pytest.approx([2, 4, 7])
    assert data[1] == pytest.approx([4, 4, 8])


@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_resid_ascii(session, kwargs, idval, tmp_path, skip_if_no_io):
    """Can we use save_resid?"""

    s = session()
    outfile = tmp_path / "resid.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_resid)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "RESID"]
    assert len(data) == 2
    assert data[0] == pytest.approx([2, 4, 7])
    assert data[1] == pytest.approx([-1, 0, -3])


@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_delchi_ascii(session, kwargs, idval, tmp_path, skip_if_no_io):
    """Can we use save_delchi?"""

    s = session()
    outfile = tmp_path / "resid.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_delchi)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "DELCHI"]
    assert len(data) == 2
    assert data[0] == pytest.approx([2, 4, 7])
    assert data[1] == pytest.approx([-1, 0, -0.5])


@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_data_ascii(session, kwargs, idval, tmp_path, skip_if_no_io):
    """Can we use save_data?"""

    s = session()
    outfile = tmp_path / "data.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_data)

    names, data = get_data(outfile, ncols=4)
    assert names == ["XLO", "XHI", "Y", "STATERROR"]
    assert len(data) == 4
    assert data[0] == pytest.approx([1, 3, 5])
    assert data[1] == pytest.approx([3, 5, 9])
    assert data[2] == pytest.approx([3, 4, 5])
    assert data[3] == pytest.approx([1, 2, 6])


@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_staterror_ascii(session, kwargs, idval, tmp_path, skip_if_no_io):
    """Can we use save_staterror?

    The output of this looks wrong (just showing the XLO values not
    XMID or XLO,XHI) but treat as a regression test. See #1526.

    """

    s = session()
    outfile = tmp_path / "errs.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_staterror)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "STAT_ERR"]
    assert len(data) == 2
    assert data[0] == pytest.approx([1, 3, 5])
    assert data[1] == pytest.approx([1, 2, 6])


@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_syserror_ascii(session, kwargs, idval, tmp_path, skip_if_no_io):
    """Can we use save_syserror?

    The output of this looks wrong (just showing the XLO values not
    XMID or XLO,XHI) but treat as a regression test. See #1526.

    """

    s = session()
    outfile = tmp_path / "errs.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_syserror, syserr=True)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "SYS_ERR"]
    assert len(data) == 2
    assert data[0] == pytest.approx([1, 3, 5])
    assert data[1] == pytest.approx([1, 2, 6])


@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_syserror_no_data_ascii(session, kwargs, idval, tmp_path, skip_if_no_io):
    """Check save_syserror needs a systematic-error column."""

    s = session()
    outfile = tmp_path / "errs.dat"

    mid = 1 if idval is None else idval
    with pytest.raises(DataErr,
                       match=f"^data set '{mid}' does not specify systematic errors$"):
        save_ascii_file(s, kwargs, idval, outfile, s.save_syserror)


@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_error_ascii(session, kwargs, idval, tmp_path, skip_if_no_io):
    """Can we use save_error?

    This does not bother testing the various combinations of
    statistical and systematic errors as this is assumed tested
    elsewhere.

    The output of this looks wrong (just showing the XLO values not
    XMID or XLO,XHI) but treat as a regression test.

    """

    s = session()
    outfile = tmp_path / "errs.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_error)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "ERR"]
    assert len(data) == 2
    assert data[0] == pytest.approx([1, 3, 5])
    assert data[1] == pytest.approx([1, 2, 6])


@pytest.mark.parametrize("session,flag",
                         [(Session, False),
                          (AstroSession, True)])
def test_show_data(session, flag):
    """Is show_data doing anything sensible?"""

    s = session()

    s.load_arrays(1, [-300, -100], [-150, -50], [-200, -100], Data1DInt)
    s.load_arrays(2, [1, 2, 10, 20], [3, 5, 12, 20])

    out = StringIO()
    s.show_data(outfile=out)

    toks = out.getvalue().split("\n")
    assert toks[0] == "Data Set: 1"
    idx = 1
    if flag:
        assert toks[idx] == "Filter: -300.0000--50.0000 x"
        idx += 1

    assert toks[idx] == "name      = "
    idx += 1
    assert toks[idx] == "xlo       = Int64[2]"
    idx += 1
    assert toks[idx] == "xhi       = Int64[2]"
    idx += 1
    assert toks[idx] == "y         = Int64[2]"
    idx += 1
    assert toks[idx] == "staterror = None"
    idx += 1
    assert toks[idx] == "syserror  = None"
    idx += 1
    assert toks[idx] == ""
    idx += 1
    assert toks[idx] == "Data Set: 2"
    idx += 1

    if flag:
        assert toks[idx] == "Filter: 1.0000-20.0000 x"
        idx += 1

    assert toks[idx] == "name      = "
    idx += 1
    assert toks[idx] == "x         = Int64[4]"
    idx += 1
    assert toks[idx] == "y         = Int64[4]"
    idx += 1
    assert toks[idx] == "staterror = None"
    idx += 1
    assert toks[idx] == "syserror  = None"
    idx += 1
    assert toks[idx] == ""
    idx += 1
    assert toks[idx] == ""
    idx += 1
    assert toks[idx] == ""

    n = 19 if flag else 17
    assert len(toks) == n


def test_show_data_datapha_no_bkg_no_response():
    """Is show_data doing anything sensible with PHA data (no background or response)"""

    s = AstroSession()

    s.set_default_id(2)

    s.load_arrays(1, [1, 2, 3, 4, 5], [10, 11, 12, 13, 14], DataPHA)
    s.load_arrays(2, [1, 2, 3, 4, 5], [20, 21, 22, 23, 24], DataPHA)

    s.set_exposure(200)
    s.set_exposure(1, 100)

    s.ignore_id(1, lo=2, hi=3)

    out = StringIO()
    s.show_data(outfile=out)

    toks = out.getvalue().split("\n")
    assert toks[0] == "Data Set: 1"
    assert toks[1] == "Filter: 1,4-5 Channel"
    assert toks[2] == "Noticed Channels: 1,4-5"
    assert toks[3] == "name           = "
    assert toks[4] == "channel        = Int64[5]"
    assert toks[5] == "counts         = Int64[5]"
    assert toks[6] == "staterror      = None"
    assert toks[7] == "syserror       = None"
    assert toks[8] == "bin_lo         = None"
    assert toks[9] == "bin_hi         = None"
    assert toks[10] == "grouping       = None"
    assert toks[11] == "quality        = None"
    assert toks[12] == "exposure       = 100.0"
    assert toks[13] == "backscal       = None"
    assert toks[14] == "areascal       = None"
    assert toks[15] == "grouped        = False"
    assert toks[16] == "subtracted     = False"
    assert toks[17] == "units          = channel"
    assert toks[18] == "rate           = True"
    assert toks[19] == "plot_fac       = 0"
    assert toks[20] == "response_ids   = []"
    assert toks[21] == "background_ids = []"
    assert toks[22] == ""
    assert toks[23] == "Data Set: 2"
    assert toks[24] == "Filter: 1-5 Channel"
    assert toks[25] == "Noticed Channels: 1-5"
    assert toks[26] == "name           = "
    assert toks[27] == "channel        = Int64[5]"
    assert toks[28] == "counts         = Int64[5]"
    assert toks[29] == "staterror      = None"
    assert toks[30] == "syserror       = None"
    assert toks[31] == "bin_lo         = None"
    assert toks[32] == "bin_hi         = None"
    assert toks[33] == "grouping       = None"
    assert toks[34] == "quality        = None"
    assert toks[35] == "exposure       = 200.0"
    assert toks[36] == "backscal       = None"
    assert toks[37] == "areascal       = None"
    assert toks[38] == "grouped        = False"
    assert toks[39] == "subtracted     = False"
    assert toks[40] == "units          = channel"
    assert toks[41] == "rate           = True"
    assert toks[42] == "plot_fac       = 0"
    assert toks[43] == "response_ids   = []"
    assert toks[44] == "background_ids = []"
    assert toks[45] == ""
    assert toks[46] == ""
    assert toks[47] == ""

    assert len(toks) == 48


def test_show_data_datapha_bkg_no_response():
    """Is show_data doing anything sensible with PHA data (muptiple backgrounds, no response)"""

    s = AstroSession()

    chans = numpy.arange(1, 6, dtype=int)
    counts = numpy.asarray([10, 20, 15, 12, 10], dtype=int)
    data = DataPHA("src", chans, counts)
    bkg1 = DataPHA("down", chans, counts)
    bkg2 = DataPHA("up", chans, counts)

    s.set_data(data)
    s.set_bkg(1, bkg1)
    s.set_bkg(1, bkg2, bkg_id=2)

    s.set_exposure(400)
    s.set_exposure(200, bkg_id=1)
    s.set_exposure(100, bkg_id=2)

    # filter one of the background components
    s.ignore_id(lo=2, hi=3, ids=1, bkg_id=1)

    out = StringIO()
    s.show_data(outfile=out)

    toks = out.getvalue().split("\n")
    assert toks[0] == "Data Set: 1"
    assert toks[1] == "Filter: 1-5 Channel"
    assert toks[2] == "Bkg Scale 1: 1"  # TODO: what are these meant to be?
    assert toks[3] == "Bkg Scale 2: 2"  # TODO: what are these meant to be?
    assert toks[4] == "Noticed Channels: 1-5"
    assert toks[5] == "name           = src"
    assert toks[6] == "channel        = Int64[5]"
    assert toks[7] == "counts         = Int64[5]"
    assert toks[8] == "staterror      = None"
    assert toks[9] == "syserror       = None"
    assert toks[10] == "bin_lo         = None"
    assert toks[11] == "bin_hi         = None"
    assert toks[12] == "grouping       = None"
    assert toks[13] == "quality        = None"
    assert toks[14] == "exposure       = 400.0"
    assert toks[15] == "backscal       = None"
    assert toks[16] == "areascal       = None"
    assert toks[17] == "grouped        = False"
    assert toks[18] == "subtracted     = False"
    assert toks[19] == "units          = channel"
    assert toks[20] == "rate           = True"
    assert toks[21] == "plot_fac       = 0"
    assert toks[22] == "response_ids   = []"
    assert toks[23] == "background_ids = [1, 2]"
    assert toks[24] == ""
    assert toks[25] == "Background Data Set: 1:1"
    assert toks[26] == "Filter: 1,4-5 Channel"
    assert toks[27] == "Noticed Channels: 1,4-5"
    assert toks[28] == "name           = down"
    assert toks[29] == "channel        = Int64[5]"
    assert toks[30] == "counts         = Int64[5]"
    assert toks[31] == "staterror      = None"
    assert toks[32] == "syserror       = None"
    assert toks[33] == "bin_lo         = None"
    assert toks[34] == "bin_hi         = None"
    assert toks[35] == "grouping       = None"
    assert toks[36] == "quality        = None"
    assert toks[37] == "exposure       = 200.0"
    assert toks[38] == "backscal       = None"
    assert toks[39] == "areascal       = None"
    assert toks[40] == "grouped        = False"
    assert toks[41] == "subtracted     = False"
    assert toks[42] == "units          = channel"
    assert toks[43] == "rate           = True"
    assert toks[44] == "plot_fac       = 0"
    assert toks[45] == "response_ids   = []"
    assert toks[46] == "background_ids = []"
    assert toks[47] == ""
    assert toks[48] == "Background Data Set: 1:2"
    assert toks[49] == "Filter: 1-5 Channel"
    assert toks[50] == "Noticed Channels: 1-5"
    assert toks[51] == "name           = up"
    assert toks[52] == "channel        = Int64[5]"
    assert toks[53] == "counts         = Int64[5]"
    assert toks[54] == "staterror      = None"
    assert toks[55] == "syserror       = None"
    assert toks[56] == "bin_lo         = None"
    assert toks[57] == "bin_hi         = None"
    assert toks[58] == "grouping       = None"
    assert toks[59] == "quality        = None"
    assert toks[60] == "exposure       = 100.0"
    assert toks[61] == "backscal       = None"
    assert toks[62] == "areascal       = None"
    assert toks[63] == "grouped        = False"
    assert toks[64] == "subtracted     = False"
    assert toks[65] == "units          = channel"
    assert toks[66] == "rate           = True"
    assert toks[67] == "plot_fac       = 0"
    assert toks[68] == "response_ids   = []"
    assert toks[69] == "background_ids = []"
    assert toks[70] == ""
    assert toks[71] == ""
    assert toks[72] == ""
    assert len(toks) == 73


def test_show_bkg_datapha_no_response():
    """Is show_bkg doing anything sensible with PHA data (single background, no response)"""

    s = AstroSession()

    chans = numpy.arange(1, 6, dtype=int)
    counts = numpy.asarray([10, 20, 15, 12, 10], dtype=int)
    data = DataPHA("src", chans, counts)
    bkg = DataPHA("down", chans, counts)

    s.set_data(data)
    s.set_bkg(1, bkg)

    s.set_exposure(400)
    s.set_exposure(200, bkg_id=1)

    # filter one of the background components
    s.ignore_id(lo=2, hi=3, ids=1, bkg_id=1)

    out = StringIO()
    s.show_bkg(outfile=out)

    toks = out.getvalue().split("\n")
    assert toks[0] == "Background Data Set: 1:1"
    assert toks[1] == "Filter: 1,4-5 Channel"
    assert toks[2] == "Noticed Channels: 1,4-5"
    assert toks[3] == "name           = down"
    assert toks[4] == "channel        = Int64[5]"
    assert toks[5] == "counts         = Int64[5]"
    assert toks[6] == "staterror      = None"
    assert toks[7] == "syserror       = None"
    assert toks[8] == "bin_lo         = None"
    assert toks[9] == "bin_hi         = None"
    assert toks[10] == "grouping       = None"
    assert toks[11] == "quality        = None"
    assert toks[12] == "exposure       = 200.0"
    assert toks[13] == "backscal       = None"
    assert toks[14] == "areascal       = None"
    assert toks[15] == "grouped        = False"
    assert toks[16] == "subtracted     = False"
    assert toks[17] == "units          = channel"
    assert toks[18] == "rate           = True"
    assert toks[19] == "plot_fac       = 0"
    assert toks[20] == "response_ids   = []"
    assert toks[21] == "background_ids = []"
    assert toks[22] == ""
    assert toks[23] == ""
    assert toks[24] == ""

    assert len(toks) == 25


def test_show_data_datapha_bkg():
    """Is show_data doing anything sensible with PHA data (background and responses)"""

    s = AstroSession()

    chans = numpy.arange(1, 6, dtype=int)
    counts = numpy.asarray([10, 20, 15, 12, 10], dtype=int)
    data = DataPHA("src", chans, counts)
    bkg = DataPHA("bkg", chans, counts)

    # Pick a variety of bin edges. Just pick a RMF-only example.
    #
    edges = numpy.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5])
    src_rmf = create_delta_rmf(edges[:-1], edges[1:], name="srmf",
                               e_min=edges[:-1], e_max=edges[1:])
    bkg_rmf = create_delta_rmf(edges[:-1], edges[1:], name="brmf",
                               e_min=edges[:-1], e_max=edges[1:])

    s.set_data(data)
    s.set_bkg(1, bkg)

    s.set_exposure(400)
    s.set_exposure(200, bkg_id=1)

    s.set_rmf(src_rmf)
    s.set_rmf(bkg_rmf, bkg_id=1)

    out = StringIO()
    s.show_data(outfile=out)

    toks = out.getvalue().split("\n")
    assert toks[0] == "Data Set: 1"
    assert toks[1] == "Filter: 0.1000-1.5000 Energy (keV)"
    assert toks[2] == "Bkg Scale: 2"
    assert toks[3] == "Noticed Channels: 1-5"
    assert toks[4] == "name           = src"
    assert toks[5] == "channel        = Int64[5]"
    assert toks[6] == "counts         = Int64[5]"
    assert toks[7] == "staterror      = None"
    assert toks[8] == "syserror       = None"
    assert toks[9] == "bin_lo         = None"
    assert toks[10] == "bin_hi         = None"
    assert toks[11] == "grouping       = None"
    assert toks[12] == "quality        = None"
    assert toks[13] == "exposure       = 400.0"
    assert toks[14] == "backscal       = None"
    assert toks[15] == "areascal       = None"
    assert toks[16] == "grouped        = False"
    assert toks[17] == "subtracted     = False"
    assert toks[18] == "units          = energy"
    assert toks[19] == "rate           = True"
    assert toks[20] == "plot_fac       = 0"
    assert toks[21] == "response_ids   = [1]"
    assert toks[22] == "background_ids = [1]"
    assert toks[23] == ""
    assert toks[24] == "RMF Data Set: 1:1"
    assert toks[25] == "name     = srmf"
    assert toks[26] == "energ_lo = Float64[5]"
    assert toks[27] == "energ_hi = Float64[5]"
    assert toks[28] == "n_grp    = Int16[5]"
    assert toks[29] == "f_chan   = Int16[5]"
    assert toks[30] == "n_chan   = Int16[5]"
    assert toks[31] == "matrix   = Float32[5]"
    assert toks[32] == "e_min    = Float64[5]"
    assert toks[33] == "e_max    = Float64[5]"
    assert toks[34] == "detchans = 5"
    assert toks[35] == "offset   = 1"
    assert toks[36] == "ethresh  = None"
    assert toks[37] == ""
    assert toks[38] == "Background Data Set: 1:1"
    assert toks[39] == "Filter: 0.1000-1.5000 Energy (keV)"
    assert toks[40] == "Noticed Channels: 1-5"
    assert toks[41] == "name           = bkg"
    assert toks[42] == "channel        = Int64[5]"
    assert toks[43] == "counts         = Int64[5]"
    assert toks[44] == "staterror      = None"
    assert toks[45] == "syserror       = None"
    assert toks[46] == "bin_lo         = None"
    assert toks[47] == "bin_hi         = None"
    assert toks[48] == "grouping       = None"
    assert toks[49] == "quality        = None"
    assert toks[50] == "exposure       = 200.0"
    assert toks[51] == "backscal       = None"
    assert toks[52] == "areascal       = None"
    assert toks[53] == "grouped        = False"
    assert toks[54] == "subtracted     = False"
    assert toks[55] == "units          = energy"
    assert toks[56] == "rate           = True"
    assert toks[57] == "plot_fac       = 0"
    assert toks[58] == "response_ids   = [1]"
    assert toks[59] == "background_ids = []"
    assert toks[60] == ""
    assert toks[61] == "Background RMF Data Set: 1:1"
    assert toks[62] == "name     = brmf"
    assert toks[63] == "energ_lo = Float64[5]"
    assert toks[64] == "energ_hi = Float64[5]"
    assert toks[65] == "n_grp    = Int16[5]"
    assert toks[66] == "f_chan   = Int16[5]"
    assert toks[67] == "n_chan   = Int16[5]"
    assert toks[68] == "matrix   = Float32[5]"
    assert toks[69] == "e_min    = Float64[5]"
    assert toks[70] == "e_max    = Float64[5]"
    assert toks[71] == "detchans = 5"
    assert toks[72] == "offset   = 1"
    assert toks[73] == "ethresh  = None"
    assert toks[74] == ""
    assert toks[75] == ""
    assert toks[76] == ""

    assert len(toks) == 77


def test_show_bkg_source_output():
    """Very basic checks"""

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)
    s._add_model_types(sherpa.astro.models)

    chans = numpy.arange(1, 6, dtype=int)
    counts = numpy.asarray([10, 20, 15, 12, 10], dtype=int)
    data = DataPHA("src", chans, counts)
    bkg = DataPHA("bkg", chans, counts)

    s.set_data(data)
    s.set_bkg(1, bkg)

    # Pick a variety of bin edges. Just pick a RMF-only example.
    #
    edges = numpy.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5])
    src_rmf = create_delta_rmf(edges[:-1], edges[1:], name="srmf",
                               e_min=edges[:-1], e_max=edges[1:])
    bkg_rmf = create_delta_rmf(edges[:-1], edges[1:], name="brmf",
                               e_min=edges[:-1], e_max=edges[1:])

    s.set_rmf(src_rmf)
    s.set_rmf(bkg_rmf, bkg_id=1)

    s.set_exposure(400)
    s.set_exposure(200, bkg_id=1)

    other = s.create_model_component("lorentz1d", "other")
    s.set_bkg_source(other)

    out = StringIO()
    s.show_bkg_source(outfile=out)

    toks = out.getvalue().split("\n")
    assert toks[0] == "Background Source: 1:1"
    assert toks[1] == "lorentz1d.other"
    assert toks[2] == "   Param        Type          Value          Min          Max      Units"
    assert toks[3] == "   -----        ----          -----          ---          ---      -----"
    assert toks[4] == "   other.fwhm   thawed           10            0  3.40282e+38           "
    assert toks[5] == "   other.pos    thawed            1 -3.40282e+38  3.40282e+38           "
    assert toks[6] == "   other.ampl   thawed            1 -3.40282e+38  3.40282e+38           "
    assert toks[7] == ""
    assert toks[8] == ""
    assert toks[9] == ""

    assert len(toks) == 10

    out = StringIO()
    s.show_bkg_model(outfile=out)

    toks = out.getvalue().split("\n")
    assert toks[0] == "Background Model: 1:1"
    assert toks[1] == "apply_rmf(200.0 * lorentz1d.other)"
    assert toks[2] == "   Param        Type          Value          Min          Max      Units"
    assert toks[3] == "   -----        ----          -----          ---          ---      -----"
    assert toks[4] == "   other.fwhm   thawed           10            0  3.40282e+38           "
    assert toks[5] == "   other.pos    thawed            1 -3.40282e+38  3.40282e+38           "
    assert toks[6] == "   other.ampl   thawed            1 -3.40282e+38  3.40282e+38           "
    assert toks[7] == ""
    assert toks[8] == ""
    assert toks[9] == ""

    assert len(toks) == 10

    out = StringIO()
    s.show_bkg_model(bkg_id=1, outfile=out)

    toks = out.getvalue().split("\n")
    assert toks[0] == "Background Model: 1:1"
    assert toks[1] == "apply_rmf(200.0 * lorentz1d.other)"
    assert toks[2] == "   Param        Type          Value          Min          Max      Units"
    assert toks[3] == "   -----        ----          -----          ---          ---      -----"
    assert toks[4] == "   other.fwhm   thawed           10            0  3.40282e+38           "
    assert toks[5] == "   other.pos    thawed            1 -3.40282e+38  3.40282e+38           "
    assert toks[6] == "   other.ampl   thawed            1 -3.40282e+38  3.40282e+38           "
    assert toks[7] == ""
    assert toks[8] == ""
    assert toks[9] == ""

    assert len(toks) == 10

    out = StringIO()
    s.show_bkg_source(bkg_id=1, outfile=out)

    toks = out.getvalue().split("\n")
    assert toks[0] == "Background Source: 1:1"
    assert toks[1] == "lorentz1d.other"
    assert toks[2] == "   Param        Type          Value          Min          Max      Units"
    assert toks[3] == "   -----        ----          -----          ---          ---      -----"
    assert toks[4] == "   other.fwhm   thawed           10            0  3.40282e+38           "
    assert toks[5] == "   other.pos    thawed            1 -3.40282e+38  3.40282e+38           "
    assert toks[6] == "   other.ampl   thawed            1 -3.40282e+38  3.40282e+38           "
    assert toks[7] == ""
    assert toks[8] == ""
    assert toks[9] == ""

    assert len(toks) == 10


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("idval", [None, "foo"])
def test_show_filter(idval, session):
    """Is show_filter doing anything sensible?"""

    s = session()
    s.set_default_id("foo")

    s.load_arrays("foo", [1, 3, 5], [3, 4, 8], [1, 2, 3], Data1DInt)
    s.ignore(2, 4)

    out = StringIO()
    s.show_filter(id=idval, outfile=out)

    toks = out.getvalue().split("\n")

    assert toks[0] == "Data Set Filter: foo"
    assert toks[1] == "5.0000-8.0000 x"
    assert toks[2] == ""
    assert toks[3] == ""
    assert toks[4] == ""

    assert len(toks) == 5


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("idval", [None, 1])
@pytest.mark.parametrize("label", ["source", "model"])
def test_show_model(idval, label, session):
    """Is show_model/source doing anything sensible?

    In this case the output is the same for source and model.

    """

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(1, [1, 3, 5], [1, 2, 3], Data1D)
    s.load_arrays("not-interested", [200, 300], [200, 300])

    s.create_model_component("gauss1d", "g1")
    s1 = s.create_model_component("scale1d", "s1")
    s.set_source(s1)

    show = getattr(s, f"show_{label}")
    out = StringIO()
    show(id=idval, outfile=out)

    toks = out.getvalue().split("\n")

    assert toks[0] == "Model: 1"
    assert toks[1] == "scale1d.s1"

    # Assume the remaining three lines are the model, so just do a
    # minimal check.
    #
    assert toks[2].strip().startswith("Param ")
    assert toks[3].strip().startswith("----- ")
    assert toks[4].strip().startswith("s1.c0 ")

    assert toks[5] == ""
    assert toks[6] == ""
    assert toks[7] == ""

    assert len(toks) == 8


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_show_method(session):
    """Is show_method doing anything sensible?"""

    s = session()

    # Select the method so we are not affected if the default changes.
    #
    s.set_method("levmar")

    out = StringIO()
    s.show_method(outfile=out)

    toks = out.getvalue().split("\n")
    assert toks[0] == "Optimization Method: LevMar"
    assert toks[1] == "name     = levmar"
    assert toks[2] == "ftol     = 1.1920928955078125e-07"
    assert toks[3] == "xtol     = 1.1920928955078125e-07"
    assert toks[4] == "gtol     = 1.1920928955078125e-07"
    assert toks[5] == "maxfev   = None"
    assert toks[6] == "epsfcn   = 1.1920928955078125e-07"
    assert toks[7] == "factor   = 100.0"
    assert toks[8] == "numcores = 1"
    assert toks[9] == "verbose  = 0"
    assert toks[10] == ""
    assert toks[11] == ""

    assert len(toks) == 12


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_show_stat(session):
    """Is show_stat doing anything sensible?"""

    s = session()

    # Select the statistic so we are not affected if the default changes.
    #
    s.set_stat("leastsq")

    out = StringIO()
    s.show_stat(outfile=out)

    # The statistic visualization is rather large and will change
    # whenever the docstring changes, so just check the start of the
    # output.
    #
    toks = out.getvalue().split("\n")
    assert toks[0] == "Statistic: LeastSq"
    assert toks[1] == "Least Squared Statistic."
    assert toks[2] == ""


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_show_fit(session):
    """Cover a number of internal routines to check show_fit works."""

    s = session()
    s._add_model_types(sherpa.models.basic)

    # Have a model for datasets 1 and 2 but only fit to 2
    s.load_arrays(1, [-300, -100], [-200, -100])
    s.load_arrays(2, [1, 2, 10, 20], [3, 5, 12, 20])

    mdl1 = s.create_model_component("const1d", "mdl1")
    mdl2 = s.create_model_component("gauss1d", "mdl2")
    mdl3 = s.create_model_component("polynom1d", "mdl3")

    s.set_source(1, mdl1 + mdl2)
    s.set_source(2, mdl1 + mdl3)
    mdl3.c0 = 0
    mdl3.c0.freeze()
    mdl3.c1.thaw()

    s.set_stat("leastsq")
    s.set_method("simplex")

    s.fit(2)
    assert s.calc_stat(2) == pytest.approx(0.9473684210526329)
    assert mdl1.c0.val == pytest.approx(2.804511278195664)
    assert mdl3.c1.val == pytest.approx(0.8721804511277418)

    out = StringIO()
    s.show_fit(outfile=out)

    # The following is somethiing that will need to be changed
    # whenever an underlying object has it's string representation
    # changed, but it should not be hard to do. The main issue is
    # going to be numeric precision.
    #
    toks = out.getvalue().split("\n")
    assert toks[0] == "Optimization Method: NelderMead"
    assert toks[1] == "name         = simplex"
    assert toks[2] == "ftol         = 1.1920928955078125e-07"
    assert toks[3] == "maxfev       = None"
    assert toks[4] == "initsimplex  = 0"
    assert toks[5] == "finalsimplex = 9"
    assert toks[6] == "step         = None"
    assert toks[7] == "iquad        = 1"
    assert toks[8] == "verbose      = 0"
    assert toks[9] == "reflect      = True"
    assert toks[10] == ""
    assert toks[11] == "Statistic: LeastSq"
    assert toks[12] == "Least Squared Statistic."
    assert toks[13] == ""
    assert toks[14] == "    The least-square statistic is equivalent to a chi-square"
    assert toks[15] == "    statistic where the error on each point - sigma(i) - is 1."
    assert toks[16] == ""
    assert toks[17] == "    "
    assert toks[18] == ""
    assert toks[19] == "Fit:Dataset               = 2"
    assert toks[20] == "Method                = neldermead"
    assert toks[21] == "Statistic             = leastsq"
    assert toks[22] == "Initial fit statistic = 502"
    assert toks[23] == "Final fit statistic   = 0.947368 at function evaluation 233"
    assert toks[24] == "Data points           = 4"
    assert toks[25] == "Degrees of freedom    = 2"
    assert toks[26] == "Change in statistic   = 501.053"
    assert toks[27] == "   mdl1.c0        2.80451     "
    assert toks[28] == "   mdl3.c1        0.87218     "
    assert toks[29] == ""
    assert toks[30] == ""
    assert toks[31] == ""

    assert len(toks) == 32


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("label", ["proj", "unc"])
def test_get_int_xxx_recalc_false(session, label):
    """Check we call the recalc=False path, even with no data loaded."""

    s = session()
    getfunc= getattr(s, f"get_int_{label}")
    plotobj = getfunc(recalc=False)

    name = "Uncertainty" if label == "unc" else "Projection"
    expected = getattr(sherpa.plot, f"Interval{name}")
    assert isinstance(plotobj, expected)

    # check it's empty (only need a single field check)
    assert plotobj.x is None


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("label", ["proj", "unc"])
def test_get_reg_xxx_recalc_false(session, label):
    """Check we call the recalc=False path, even with no data loaded."""

    s = session()
    getfunc= getattr(s, f"get_reg_{label}")
    plotobj = getfunc(recalc=False)

    name = "Uncertainty" if label == "unc" else "Projection"
    expected = getattr(sherpa.plot, f"Region{name}")
    assert isinstance(plotobj, expected)

    # check it's empty (only need a single field check)
    assert plotobj.x0 is None


def setup_template_model(session, make_data_path, interp="default"):
    """Create the template model.

    The interp argument is passed to load_template_model as the
    template_interpolator_name field: we currently don't have a lot of
    documentation on what this is, but there is a code path which
    allows it to be set to None, which currently triggers different
    code paths, so it has been added to test this.

    The logic is test_309() from
    sherpa/models/tests/test_template_unit.py as this appears to be
    the only test of template models with the UI we have.

    """

    ynorm = 1e9  # make the values nicer

    # Need to load the data from the same directory as the index
    basedir = os.getcwd()
    os.chdir(make_data_path(""))
    try:
        session.load_template_model("bbtemp", "bb_index.dat",
                                    template_interpolator_name=interp)
    finally:
        os.chdir(basedir)

    bbtemp = session.get_model_component("bbtemp")
    return bbtemp, ynorm


@requires_data
@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("label", ["model", "source"])
def test_get_xxx_component_plot_with_templates_recalc_true(session, label, make_data_path, skip_if_no_io):
    """What is the intended behavior for template models?

    This is a regression test as it is not obvious what the intended
    logic is. It was added to ensure a code path was tested, but
    this code path is highly-dependent on how the underlying code
    is written.

    """

    # This sets up a template model for the default id, so try
    # calling with this model but a different dataset.
    #
    s = session()
    dname = make_data_path('load_template_with_interpolation-bb_data.dat')
    s.load_data(dname)

    bbtemp, ynorm = setup_template_model(s, make_data_path)
    s.set_source(ynorm * bbtemp)

    with pytest.raises(IdentifierErr,
                       match="data set 2 has not been set"):
        getattr(s, f"get_{label}_component_plot")(2, bbtemp, recalc=True)


@requires_data
@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("label", ["model", "source"])
def test_get_xxx_component_plot_with_templates_recalc_true_no_interp(session, label, make_data_path, skip_if_no_io):
    """What is the intended behavior for template models?

    Try to catch all code paths.

    """

    # This sets up a template model for the default id, so try
    # calling with this model but a different dataset.
    #
    s = session()
    dname = make_data_path('load_template_with_interpolation-bb_data.dat')
    s.load_data(dname)

    bbtemp, ynorm = setup_template_model(s, make_data_path, interp=None)
    s.set_source(ynorm * bbtemp)

    with pytest.raises(IdentifierErr,
                       match="data set 2 has not been set"):
        getattr(s, f"get_{label}_component_plot")(2, bbtemp, recalc=True)


@requires_data
@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("label", ["model", "source"])
def test_get_xxx_component_plot_with_templates_recalc_false(session, label, make_data_path, skip_if_no_io):
    """What is the intended behavior for template models?

    This is a regression test as it is not obvious what the intended
    logic is. It was added to ensure a code path was tested, but
    this code path is highly-dependent on how the underlying code
    is written.

    """

    # This sets up a template model for the default id, so try
    # calling with this model but a different dataset.
    #
    s = session()
    dname = make_data_path('load_template_with_interpolation-bb_data.dat')
    s.load_data(dname)

    bbtemp, ynorm = setup_template_model(s, make_data_path)
    s.set_source(ynorm * bbtemp)

    # The plot results should be empty
    #
    mplot = getattr(s, f"get_{label}_plot")(2, recalc=False)
    cplot = getattr(s, f"get_{label}_component_plot")(2, bbtemp, recalc=False)

    mclass = getattr(sherpa.plot, f"{label.capitalize()}Plot")
    assert isinstance(mplot, mclass)

    cclass = getattr(sherpa.plot, f"Component{label.capitalize()}Plot")
    assert isinstance(cplot, cclass)

    # Check it's not a sub-class.
    #
    subclass = getattr(sherpa.plot, f"ComponentTemplate{label.capitalize()}Plot")
    assert not isinstance(cplot, subclass)

    assert mplot.x is None
    assert cplot.x is None
    assert mplot.title == label.capitalize()
    assert cplot.title == label.capitalize()


@requires_data
@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("label", ["model", "source"])
def test_get_xxx_component_plot_with_templates_recalc_false_no_interp(session, label, make_data_path, skip_if_no_io):
    """What is the intended behavior for template models?

    This is a regression test.

    """

    # This sets up a template model for the default id, so try
    # calling with this model but a different dataset.
    #
    s = session()
    dname = make_data_path('load_template_with_interpolation-bb_data.dat')
    s.load_data(dname)

    bbtemp, ynorm = setup_template_model(s, make_data_path, interp=None)
    s.set_source(ynorm * bbtemp)

    # The plot results should be empty
    #
    mplot = getattr(s, f"get_{label}_plot")(2, recalc=False)
    cplot = getattr(s, f"get_{label}_component_plot")(2, bbtemp, recalc=False)

    mclass = getattr(sherpa.plot, f"{label.capitalize()}Plot")
    assert isinstance(mplot, mclass)

    cclass = getattr(sherpa.plot, f"Component{label.capitalize()}Plot")
    assert isinstance(cplot, cclass)

    # The logic here is not-at-all obvious, but just test the current
    # behavior.
    #
    subclass = getattr(sherpa.plot, f"ComponentTemplate{label.capitalize()}Plot")
    if label == "source":
        assert isinstance(cplot, subclass)
    else:
        assert not isinstance(cplot, subclass)

    assert mplot.x is None
    assert cplot.x is None
    assert mplot.title == label.capitalize()
    assert cplot.title == label.capitalize()


@requires_data
@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("label", ["model", "source"])
def test_get_xxx_component_plot_with_templates_data1d(session, label, make_data_path, skip_if_no_io):
    """What is the intended behavior for template models?  Data1D

    As this is not a PHA data set the source and model values should
    be the same (although plot labels will be different).

    This is a regression test as it is not obvious what the intended
    logic is.

    """

    s = session()
    dname = make_data_path('load_template_with_interpolation-bb_data.dat')
    s.load_data(dname)
    bbtemp, ynorm = setup_template_model(s, make_data_path)
    s.set_source(bbtemp * ynorm)

    # Because the source model is a composite, the y values of the
    # two should differ by ynorm.
    #
    mplot = getattr(s, f"get_{label}_plot")()
    cplot = getattr(s, f"get_{label}_component_plot")(bbtemp)

    mclass = getattr(sherpa.plot, f"{label.capitalize()}Plot")
    assert isinstance(mplot, mclass)

    cclass = getattr(sherpa.plot, f"Component{label.capitalize()}Plot")
    assert isinstance(cplot, cclass)

    # Check it's not a sub-class.
    #
    subclass = getattr(sherpa.plot, f"ComponentTemplate{label.capitalize()}Plot")
    assert not isinstance(cplot, subclass)

    # Make sure we have a valid set of y values, even if we do
    # not check the actual values.
    #
    assert mplot.y == pytest.approx(cplot.y * ynorm)
    assert numpy.all(mplot.y > 0)


@requires_data
@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("label", ["model", "source"])
def test_get_xxx_component_plot_with_templates_data1d_no_interp(session, label, make_data_path, skip_if_no_io):
    """What is the intended behavior for template models?  Data1D, no interpolator

    As this is not a PHA data set the source and model values should
    be the same (although plot labels will be different).

    This is a regression test as it is not obvious what the intended
    logic is.

    """

    s = session()
    dname = make_data_path('load_template_with_interpolation-bb_data.dat')
    s.load_data(dname)
    bbtemp, ynorm = setup_template_model(s, make_data_path, interp=None)
    s.set_source(bbtemp * ynorm)

    # Because the source model is a composite, the y values of the
    # two should differ by ynorm.
    #
    mplot = getattr(s, f"get_{label}_plot")()
    cplot = getattr(s, f"get_{label}_component_plot")(bbtemp)

    mclass = getattr(sherpa.plot, f"{label.capitalize()}Plot")
    assert isinstance(mplot, mclass)

    cclass = getattr(sherpa.plot, f"Component{label.capitalize()}Plot")
    assert isinstance(cplot, cclass)

    # The logic here is not-at-all obvious, but just test the current
    # behavior.
    #
    subclass = getattr(sherpa.plot, f"ComponentTemplate{label.capitalize()}Plot")
    if label == "source":
        assert isinstance(cplot, subclass)
    else:
        assert not isinstance(cplot, subclass)

    # Make sure we have a valid set of y values, even if we do
    # not check the actual values.
    #
    assert mplot.y == pytest.approx(cplot.y * ynorm)
    assert numpy.all(mplot.y > 0)


@requires_data
@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("label", ["model", "source"])
def test_get_xxx_component_plot_with_templates_data1dint(session, label, make_data_path):
    """What is the intended behavior for template models?  Data1DInt"""

    s = session()

    # The actual data isn't too relevant, but pick similar x range
    # to load_template_with_interpolation-bb_data.dat
    #
    edges = numpy.arange(1, 3, 0.1) * 1e14
    s.load_arrays(1, edges[:-1], edges[1:], edges[1:] * 0, Data1DInt)

    bbtemp, ynorm = setup_template_model(s, make_data_path)
    s.set_source(bbtemp * ynorm)

    # Because the source model is a composite, the y values of the
    # two should differ by ynorm.
    #
    mplot = getattr(s, f"get_{label}_plot")()
    cplot = getattr(s, f"get_{label}_component_plot")(bbtemp)

    mclass = getattr(sherpa.plot, f"{label.capitalize()}HistogramPlot")
    assert isinstance(mplot, mclass)

    cclass = getattr(sherpa.plot, f"Component{label.capitalize()}HistogramPlot")
    assert isinstance(cplot, cclass)

    # Check it's not a sub-class.
    #
    subclass = getattr(sherpa.plot, f"ComponentTemplate{label.capitalize()}Plot")
    assert not isinstance(cplot, subclass)

    # Make sure we have a valid set of y values, even if we do
    # not check the actual values.
    #
    assert mplot.y == pytest.approx(cplot.y * ynorm)
    assert numpy.all(mplot.y > 0)


@requires_data
@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("label", ["model", "source"])
def test_get_xxx_component_plot_with_templates_data1dint_no_interp(session, label, make_data_path):
    """What is the intended behavior for template models?  Data1DInt, no interpolator"""

    s = session()

    # The actual data isn't too relevant, but pick similar x range
    # to load_template_with_interpolation-bb_data.dat
    #
    edges = numpy.arange(1, 3, 0.1) * 1e14
    s.load_arrays(1, edges[:-1], edges[1:], edges[1:] * 0, Data1DInt)

    bbtemp, ynorm = setup_template_model(s, make_data_path, interp=None)
    s.set_source(bbtemp * ynorm)

    # Because the source model is a composite, the y values of the
    # two should differ by ynorm.
    #
    mplot = getattr(s, f"get_{label}_plot")()
    cplot = getattr(s, f"get_{label}_component_plot")(bbtemp)

    mclass = getattr(sherpa.plot, f"{label.capitalize()}HistogramPlot")
    assert isinstance(mplot, mclass)

    # The logic here is not-at-all obvious, but just test the current
    # behavior.
    #
    # NOTE THE DIFFERENCE TO Data1DInt with interpolator when label=spirce
    #
    if label == "source":
        assert isinstance(cplot, sherpa.plot.ComponentTemplateSourcePlot)

        # Unlike the normal Data1DInt case we can not compare the two
        # plots, so just add some regression tests for now.
        #
        assert len(mplot.y) == 19
        assert len(cplot.y) == 100

        assert mplot.y[0] == pytest.approx(cplot.y[0] * ynorm)
        assert numpy.all(mplot.y > 0)
        assert numpy.all(cplot.y > 0)

    else:

        assert isinstance(cplot, sherpa.plot.ComponentModelHistogramPlot)

        # Check it's not a sub-class.
        #
        assert not isinstance(cplot, sherpa.plot.ComponentTemplateModelPlot)

        # Make sure we have a valid set of y values, even if we do
        # not check the actual values.
        #
        assert mplot.y == pytest.approx(cplot.y * ynorm)
        assert numpy.all(mplot.y > 0)


@requires_data
@pytest.mark.parametrize("interp", ["default", None])
def test_get_source_component_plot_with_templates_datapha(interp, make_data_path):
    """What is the intended behavior for template models?  Data1PHA, with and without interpolator

    It would be nice to be able to check both source and model
    commands with a single test, but the behavior is not guaranteed to
    match.  This is a regression test as it's not clear what the
    behavior should be.

    """

    s = AstroSession()

    # The actual data isn't too relevant, but pick similar x range
    # to load_template_with_interpolation-bb_data.dat
    #
    chans = numpy.arange(1, 100, dtype=numpy.int16)
    s.load_arrays(1, chans, chans * 0, DataPHA)

    edges = numpy.linspace(1, 3, num=len(chans) + 1) * 1e4
    arf = create_arf(edges[:-1], edges[1:])
    s.set_arf(arf)

    bbtemp, ynorm = setup_template_model(s, make_data_path, interp=interp)
    s.set_source(bbtemp * ynorm)

    # Because the source model is a composite, the y values of the
    # two should differ by ynorm.
    #
    mplot = s.get_source_plot()
    cplot = s.get_source_component_plot(bbtemp)

    assert isinstance(mplot, sherpa.astro.plot.SourcePlot)
    assert isinstance(cplot, sherpa.astro.plot.ComponentSourcePlot)

    # Make sure we have a valid set of y values, even if we do
    # not check the actual values.
    #
    assert mplot.y == pytest.approx(cplot.y * ynorm)
    assert numpy.all(mplot.y > 0)


@requires_data
@pytest.mark.parametrize("interp", ["default", None])
def test_get_model_component_plot_with_templates_datapha(interp, make_data_path):
    """What is the intended behavior for template models?  Data1PHA, with and without interpolator

    It would be nice to be able to check both source and model
    commands with a single test, but the behavior is not guaranteed to
    match.  This is a regression test as it's not clear what the
    behavior should be.

    """

    s = AstroSession()

    # The actual data isn't too relevant, but pick similar x range
    # to load_template_with_interpolation-bb_data.dat
    #
    chans = numpy.arange(1, 100, dtype=numpy.int16)
    s.load_arrays(1, chans, chans * 0, DataPHA)

    edges = numpy.linspace(1, 3, num=len(chans) + 1) * 1e4
    arf = create_arf(edges[:-1], edges[1:])
    s.set_arf(arf)

    bbtemp, ynorm = setup_template_model(s, make_data_path, interp=interp)
    s.set_source(bbtemp * ynorm)

    # Because the source model is a composite, the y values of the
    # two should differ by ynorm.
    #
    mplot = s.get_model_plot()
    cplot = s.get_model_component_plot(bbtemp)

    assert isinstance(mplot, sherpa.astro.plot.ModelHistogram)
    assert isinstance(cplot, sherpa.astro.plot.ComponentModelPlot)

    # Make sure we have a valid set of y values, even if we do
    # not check the actual values.
    #
    assert mplot.y == pytest.approx(cplot.y * ynorm)
    assert numpy.all(mplot.y > 0)


@requires_data
@pytest.mark.parametrize("interp", ["default", None])
def test_compare_get_model_component_plot_with_templates(interp, make_data_path):
    """Check what happens comparing DataPHA with Data1DInt and Data1D

    This is a regression test as it is not clear what the intended
    behavior is.

    """

    NBINS = 100
    SPECRESP = 0.8
    BIN_WIDTH = 200

    s = AstroSession()

    # Create a dataset with the same "x" axis. This is not technically
    # possible, but knowing how model evaluation works (the first
    # array is normally used for 1D non-integrated) we can try
    # and get the same values used.
    #
    chans = numpy.arange(1, NBINS + 1, dtype=numpy.int16)
    s.load_arrays("pha", chans, chans * 0, DataPHA)

    # Give the ARF a non-unit response so we can see if it has been
    # applied (or, can infer it has been applied since the model
    # component will apply it).
    #
    edges = numpy.arange(10000, 30001, BIN_WIDTH)
    arf = create_arf(edges[:-1], edges[1:], numpy.ones(NBINS) * SPECRESP)
    s.set_arf("pha", arf)

    ones = numpy.ones(NBINS)
    s.load_arrays("int", edges[:-1], edges[1:], ones, Data1DInt)
    s.load_arrays("1d", edges[:-1], ones, Data1D)

    bbtemp, ynorm = setup_template_model(s, make_data_path, interp=interp)
    s.set_source("pha", bbtemp * ynorm)
    s.set_source("int", bbtemp * ynorm)
    s.set_source("1d", bbtemp * ynorm)

    mplot_pha = s.get_model_plot("pha")
    cplot_pha = s.get_model_component_plot("pha", bbtemp)

    mplot_int = s.get_model_plot("int")
    cplot_int = s.get_model_component_plot("int", bbtemp)

    mplot_1d = s.get_model_plot("1d")
    cplot_1d = s.get_model_component_plot("1d", bbtemp)

    # In case the pathway is different from using the default dataset
    # identifier, check some basic things about the response.
    #
    assert numpy.all(cplot_pha.y > 0)
    assert numpy.all(cplot_int.y > 0)
    assert numpy.all(cplot_1d.y > 0)

    assert cplot_int.xlo == pytest.approx(cplot_pha.xlo)
    assert cplot_1d.x == pytest.approx(cplot_pha.xlo)

    # It looks like the model is not integrated across the bin for
    # the 1DInt case.
    #
    assert cplot_1d.y == pytest.approx(cplot_int.y)

    # The DataPHA result, when compared to Data1D is:
    #
    # a) multiplied by SPECRESP
    # b) divided by BIN_WIDTH
    #
    r = cplot_int.y / cplot_pha.y
    assert r == pytest.approx(ones * BIN_WIDTH / SPECRESP)


@requires_data
def test_get_source_component_plot_with_templates_datapha_no_response(make_data_path):
    """What is the behavior with no response?

    It would be nice to be able to check both source and model
    commands with a single test, but the behavior is not guaranteed to
    match.  This is a regression test as it's not clear what the
    behavior should be.

    """

    s = AstroSession()

    chans = numpy.arange(1, 100, dtype=numpy.int16)
    s.load_arrays(1, chans, chans * 0, DataPHA)

    bbtemp, ynorm = setup_template_model(s, make_data_path)
    s.set_source(bbtemp * ynorm)

    with pytest.raises(DataErr,
                       match="Response does not specify energy bins"):
        s.get_source_plot()

    with pytest.raises(DataErr,
                       match="Response does not specify energy bins"):
        s.get_source_component_plot(bbtemp)


@requires_data
def test_get_model_component_plot_with_templates_datapha_no_response(make_data_path):
    """What is the behavior with no response?

    It would be nice to be able to check both source and model
    commands with a single test, but the behavior is not guaranteed to
    match.  This is a regression test as it's not clear what the
    behavior should be.

    """

    s = AstroSession()

    chans = numpy.arange(1, 100, dtype=numpy.int16)
    s.load_arrays(1, chans, chans * 0, DataPHA)

    bbtemp, ynorm = setup_template_model(s, make_data_path)
    s.set_source(bbtemp * ynorm)

    with pytest.raises(DataErr,
                       match="^No instrument response found for dataset $"):
        s.get_model_plot()

    cplot = s.get_model_component_plot(bbtemp)
    assert cplot.title == "Model component: template.bbtemp"
    assert numpy.all(cplot.y > 0)


def check_stat_info_basic(sinfo, name, ids, numpoints, statval):
    """Check the "basic" result"""

    assert sinfo.name == name
    assert sinfo.ids == ids
    assert sinfo.bkg_ids is None  # TODO: why does Session class have this?
    assert sinfo.numpoints == numpoints
    assert sinfo.dof == (numpoints - 1)
    assert sinfo.qval is None
    assert sinfo.rstat is None

    assert sinfo.statname == "leastsq"
    assert sinfo.statval == pytest.approx(statval)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_stat_info_basic_one(session):
    """Check get_stat_info with one dataset"""

    s = session()
    s._add_model_types(sherpa.models.basic)
    s.set_stat("leastsq")

    data = Data1D("example", [1, 2, 5], [3, 7, 6])
    s.set_data(data)

    cpt = s.create_model_component("scale1d", "scale")
    cpt.c0 = 5
    s.set_source(cpt)

    sinfo = s.get_stat_info()
    assert len(sinfo) == 1

    check_stat_info_basic(sinfo[0], "Dataset [1]", [1], 3, 9)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_stat_info_basic_two(session):
    """Check get_stat_info with two datasets"""

    s = session()
    s._add_model_types(sherpa.models.basic)
    s.set_stat("leastsq")

    data1 = Data1D("example1", [1, 2, 5], [3, 7, 6])
    data2 = Data1D("example2", [1, 2, 5], [4, 6, 5])
    s.set_data(1, data1)
    s.set_data(2, data2)

    cpt = s.create_model_component("scale1d", "scale")
    cpt.c0 = 5
    s.set_source(1, cpt)
    s.set_source(2, cpt)

    sinfo = s.get_stat_info()
    assert len(sinfo) == 3

    check_stat_info_basic(sinfo[0], "Dataset 1", (1, ), 3, 9)
    check_stat_info_basic(sinfo[1], "Dataset 2", (2, ), 3, 2)
    check_stat_info_basic(sinfo[2], "Datasets [1, 2]", [1, 2], 6, 11)


def test_get_stat_info_astro_one(caplog):
    """Check get_stat_info with one dataset and a background"""

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)
    s.set_stat("leastsq")

    data = DataPHA("example", [1, 2, 3], [3, 7, 6])
    bkg = DataPHA("background", [1, 2, 3], [1, 1, 2])

    egrid = numpy.asarray([0.1, 0.2, 0.4, 0.8])
    arf = create_arf(egrid[:-1], egrid[1:])
    data.set_arf(arf)
    bkg.set_arf(arf)
    data.set_background(bkg)

    s.set_data(data)
    s.subtract()

    cpt = s.create_model_component("scale1d", "scale")
    cpt.c0 = 5
    s.set_source(cpt)

    assert len(caplog.record_tuples) == 0
    sinfo = s.get_stat_info()
    assert len(caplog.record_tuples) == 0
    assert len(sinfo) == 1

    check_stat_info_basic(sinfo[0], "Dataset [1]", [1], 3, 11)


def test_get_stat_info_astro_two(caplog):
    """Check get_stat_info with two datasets and a background"""

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)
    s.set_stat("leastsq")

    data1 = DataPHA("example1", [1, 2, 3], [3, 7, 6])
    data2 = DataPHA("example2", [1, 2, 3], [4, 6, 5])

    bkg2 = DataPHA("background2", [1, 2, 3], [1, 1, 2])

    egrid = numpy.asarray([0.1, 0.2, 0.4, 0.8])
    arf = create_arf(egrid[:-1], egrid[1:])
    data1.set_arf(arf)
    data2.set_arf(arf)

    bkg2.set_arf(arf)
    data2.set_background(bkg2)

    s.set_data(1, data1)
    s.set_data(2, data2)

    s.subtract(2)

    cpt = s.create_model_component("scale1d", "scale")
    cpt.c0 = 5
    s.set_source(1, cpt)
    s.set_source(2, cpt)

    assert len(caplog.record_tuples) == 0
    sinfo = s.get_stat_info()
    assert len(caplog.record_tuples) == 0
    assert len(sinfo) == 3

    check_stat_info_basic(sinfo[0], "Dataset 1", (1, ), 3, 9)
    check_stat_info_basic(sinfo[1], "Dataset 2", (2, ), 3, 8)
    check_stat_info_basic(sinfo[2], "Datasets [1, 2]", [1, 2], 6, 17)


class DummyClass:
    "A dummy class"
    pass


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_add_model_errors_out(session):
    """The model class needs to be derived from ArithmeticModel."""

    msg = "^model class 'dummyclass' is not a derived class from " + \
        "sherpa.models.ArithmeticModel$"

    s = session()
    with pytest.raises(TypeError, match=msg):
        s.add_model(DummyClass)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_source_with_convolved_model(session):
    """Check we get an error.

    The obvous way to check this would be with a PHA dataset
    but then we can not test the ui.Session, so here we
    use set_full_model to simulate this (the fact that
    this is not really a case that needs set_full_model is
    not relevant here).

    It turns out that with Session we do not need to set a dataset but
    we do with AstroSession.

    """

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(1, [1, 2], [2, 4])

    gmdl = s.create_model_component("gauss1d", "gmdl")
    s.set_full_model(gmdl)

    msg = "^Convolved model\n'gauss1d.gmdl'\n is set for " + \
        "dataset 1. You should use get_model instead.$"

    with pytest.raises(IdentifierErr, match=msg):
        s.get_source()


def test_set_full_model_pha_warning_no_response(caplog):
    """Check we do not get the warning.

    This requires a PHA dataset with no response.
    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(1, [1, 2], [2, 4], DataPHA)

    gmdl = s.create_model_component("gauss1d", "gmdl")

    with caplog.at_level(logging.INFO, logger='sherpa'):
        s.set_full_model(gmdl)

    assert len(caplog.record_tuples) == 0


def test_set_full_model_pha_warning_no_response_bkg(caplog):
    """Check we do not get the warning.

    This requires a PHA dataset with no response.
    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(1, [1, 2], [2, 4], DataPHA)
    s.set_bkg(1, DataPHA("b", [1, 2], [1, 2]))

    gmdl = s.create_model_component("gauss1d", "gmdl")

    with caplog.at_level(logging.INFO, logger='sherpa'):
        s.set_bkg_full_model(gmdl)

    assert len(caplog.record_tuples) == 0


def test_set_full_model_pha_warning_response(caplog):
    """Check we do get the warning.

    This requires a PHA dataset with a response.
    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(1, [1, 2], [2, 4], DataPHA)

    egrid = numpy.asarray([0.1, 0.2, 0.3])
    s.set_arf(create_arf(egrid[:-1], egrid[1:], [10, 20]))

    gmdl = s.create_model_component("gauss1d", "gmdl")
    assert len(caplog.record_tuples) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        s.set_full_model(gmdl)

    assert len(caplog.record_tuples) == 1
    loc, lvl, msg = caplog.record_tuples[0]
    assert loc == "sherpa.astro.ui.utils"
    assert lvl == logging.WARNING
    assert msg == "PHA source model 'gauss1d.gmdl' \n" + \
        "does not have an associated instrument model; consider using \n" + \
        "set_source() instead of set_full_model() to include associated \n" + \
        "instrument automatically"


def test_set_full_model_pha_warning_response_bkg():
    """Check we do get the warning. Actually, this is an error.

    This requires a PHA dataset with a response.
    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(1, [1, 2], [2, 4], DataPHA)
    s.set_bkg(1, DataPHA("b", [1, 2], [1, 2]))

    egrid = numpy.asarray([0.1, 0.2, 0.3])
    s.set_arf(create_arf(egrid[:-1], egrid[1:], [10, 20]))
    s.set_arf(create_arf(egrid[:-1], egrid[1:], [5, 10]), bkg_id=1)

    gmdl = s.create_model_component("gauss1d", "gmdl")

    msg = "PHA background source model 'gauss1d.gmdl' \n" + \
          " does not have an associated instrument model; consider using\n" + \
          " set_bkg_source() instead of set_bkg_model() to include associated\n" + \
          " instrument automatically"
    with pytest.raises(TypeError, match=re.escape(msg)):
        s.set_bkg_full_model(gmdl)


def test_bkg_model_warns_after_full(caplog):
    """Check we get a warning"""

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    # It looks like you can't use load_arrays to set the background.
    #
    s.load_arrays(1, [1, 2, 3], [10, 11, 12], DataPHA)

    bkg = DataPHA("bg", [1, 2, 3], [2, 1, 3])
    s.set_bkg(bkg)

    # It does not matter what model is being used.
    #
    cpt1 = s.create_model_component("gauss1d", "g1")
    cpt2 = s.create_model_component("polynom1d", "p1")

    assert len(caplog.record_tuples) == 0
    s.set_bkg_full_model(cpt1)
    assert len(caplog.record_tuples) == 0
    s.set_bkg_model(cpt2)
    assert len(caplog.record_tuples) == 1

    loc, lvl, msg = caplog.record_tuples[0]
    assert loc == "sherpa.astro.ui.utils"
    assert lvl == logging.WARNING
    assert msg == "Clearing background convolved model\n'gauss1d.g1'\nfor dataset 1 background 1"


# Note: the Session message can be either fit() or Fit.fit(), so
# fortunately it is easy to check (it is not obvious what causes
# the difference, but it has been seen on a CI run).
#
@pytest.mark.parametrize("session,msg",
                         [(Session, r"fit\(\) got an unexpected keyword argument 'unknown_argument'"),
                          (AstroSession, "unknown keyword argument: 'unknown_argument'")])
def test_fit_checks_kwarg(session, msg):
    """Check what happens if fit is sent an unknown argument.

    This is just a regression test so we know if anything ever changes.

    """

    s = session()
    s._add_model_types(sherpa.models.basic)
    s.set_stat("leastsq")

    data = Data1D("example", [1, 2, 5], [3, 7, 6])
    s.set_data(data)

    cpt = s.create_model_component("scale1d", "scale")
    cpt.c0 = 5
    s.set_source(cpt)

    with pytest.raises(TypeError, match=msg):
        s.fit(unknown_argument=True)


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("alias,original",
                         [("compsource", "source_component"),
                          ("compmodel", "model_component")])
def test_plot_alias_warning(session, alias, original, caplog):
    """Check we get a deprecated warning from using a plot alias.

    Support for aliases is intended to be short-term, but ensure
    they are tested. This is only relevant for the set_xlog/...
    family of commands.
    """

    s = session()
    assert len(caplog.record_tuples) == 0
    s.set_xlog(alias)
    assert len(caplog.record_tuples) == 1

    loc, lvl, msg = caplog.record_tuples[0]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.WARNING
    assert msg == f"The argument '{alias}' is deprecated and '{original}' should be used instead"


@pytest.mark.parametrize("session,success", [(Session, False), (AstroSession, True)])
@pytest.mark.parametrize("key", ["model", "fit", "source", "ratio", "resid", "delchi", "chisqr"])
def test_astro_plot_alias_warning(session, success, key, caplog):
    """Check we get a deprecated warning from using a plot alias.

    Support for aliases is intended to be short-term, but ensure
    they are tested. This is only relevant for the set_xlog/...
    family of commands.
    """

    alias = f"bkg{key}"
    original = f"bkg_{key}"

    s = session()
    assert len(caplog.record_tuples) == 0

    if success:
        s.set_xlog(alias)
        assert len(caplog.record_tuples) == 1

        loc, lvl, msg = caplog.record_tuples[0]
        assert loc == "sherpa.ui.utils"
        assert lvl == logging.WARNING
        assert msg == f"The argument '{alias}' is deprecated and '{original}' should be used instead"

    else:
        # Check this errors out (i.e. is not an alias).
        #
        with pytest.raises(PlotErr, match=rf"^Plot type '{alias}' not found in \[.*\]$"):
            s.set_xlog(alias)


def test_notice_warning(caplog):
    """Check we get a warning from notice

    This requires PHA files with different analysis settings
    """

    s = AstroSession()
    s.load_arrays(1, [1, 2, 3], [1, 2, 3], DataPHA)
    s.load_arrays(2, [1, 2, 3], [1, 2, 3], DataPHA)

    egrid = numpy.asarray([0.5, 0.7, 1.0, 2.0])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    s.set_rmf(2, rmf)

    assert len(caplog.record_tuples) == 0
    s.set_analysis(2, "energy")
    assert len(caplog.record_tuples) == 1
    clc_filter(caplog, "dataset 2: 0.5:2 Energy (keV)", astro=True)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        s.notice(lo=2)

    assert len(caplog.record_tuples) == 4
    loc, lvl, msg = caplog.record_tuples[1]
    assert loc == "sherpa.astro.ui.utils"
    assert lvl == logging.WARNING
    assert msg == "not all PHA datasets have equal analysis quantities: channel, energy"

    clc_filter(caplog, "dataset 1: 1:3 -> 2:3 Channel", pos=2)
    clc_filter(caplog, "dataset 2: 0.5:2 Energy (keV) -> no data", pos=3)


def test_ignore_warning(caplog):
    """Check we get a warning from ignore (was a test of #1641)

    This used to be more-different to test_ignore_warning (before
    #1641 was fixed), but keep it a separate test.
    """

    s = AstroSession()
    s.load_arrays(1, [1, 2, 3], [1, 2, 3], DataPHA)
    s.load_arrays(2, [1, 2, 3], [1, 2, 3], DataPHA)

    egrid = numpy.asarray([0.5, 0.7, 1.0, 2.0])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    s.set_rmf(2, rmf)

    assert len(caplog.record_tuples) == 0
    s.set_analysis(2, "energy")
    assert len(caplog.record_tuples) == 1
    clc_filter(caplog, "dataset 2: 0.5:2 Energy (keV)", astro=True)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        s.ignore(lo=2)

    assert len(caplog.record_tuples) == 4
    loc, lvl, msg = caplog.record_tuples[1]
    assert loc == "sherpa.astro.ui.utils"
    assert lvl == logging.WARNING
    assert msg == "not all PHA datasets have equal analysis quantities: channel, energy"

    clc_filter(caplog, "dataset 1: 1:3 -> 1 Channel", pos=2)
    clc_filter(caplog, "dataset 2: 0.5:2 Energy (keV) (unchanged)", pos=3)



def test_set_analysis_messages(caplog):
    """What happens with multiple datasets.

    Note that set_analysis fail if a dataset is not a DataPHA, or
    there's no response, so that reduces the number of things to
    check.

    """

    s = AstroSession()
    s.load_arrays(2, [1, 2, 3], [1, 2, 3], DataPHA)
    s.load_arrays("foo", [1, 2, 3], [1, 2, 3], DataPHA)

    egrid = numpy.asarray([0.5, 0.7, 1.0, 2.0])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    s.set_rmf(2, rmf)
    s.set_rmf("foo", rmf)

    assert s.get_data(2).units == "energy"  # just check
    assert s.get_data("foo").units == "energy"  # just check

    assert len(caplog.record_tuples) == 0
    with SherpaVerbosity("INFO"):
        s.set_analysis("foo", "channel")

    assert len(caplog.record_tuples) == 1
    clc_filter(caplog, "dataset foo: 1:3 Channel", astro=True)

    with SherpaVerbosity("INFO"):
        s.set_analysis("wave")

    assert len(caplog.record_tuples) == 3
    clc_filter(caplog, "dataset 2: 6.19921:24.7968 Wavelength (Angstrom)", astro=True, pos=1)
    clc_filter(caplog, "dataset foo: 6.19921:24.7968 Wavelength (Angstrom)", astro=True, pos=2)

    # check the "no data" handling
    #
    s.ignore_id("foo")
    assert len(caplog.record_tuples) == 4
    clc_filter(caplog, "dataset foo: 6.19921:24.7968 Wavelength (Angstrom) -> no data", pos=3)

    with SherpaVerbosity("INFO"):
        s.set_analysis("wave")

    assert len(caplog.record_tuples) == 6
    clc_filter(caplog, "dataset 2: 6.19921:24.7968 Wavelength (Angstrom)", astro=True, pos=4)
    clc_filter(caplog, "dataset foo: no data", astro=True, pos=5)


def test_set_analysis_no_data():
    """What happens with no data"""

    s = AstroSession()
    with pytest.raises(IdentifierErr,
                       match="^No data sets found$"):
        s.set_analysis("energy")


@pytest.mark.parametrize("session,emsg",
                         [(Session, r"save_delchi\(\) can not be used with 2D datasets"),
                          (AstroSession, r"save_delchi\(\) does not apply for images")])
@pytest.mark.parametrize("idval", [None, "bob"])
def test_save_delchi_data2d(session, emsg, idval, tmp_path):
    """Check it errors out

    Technically there's no reason why we can't have errors for a 2D
    dataset but we currently error out.

    """

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.dataspace2d([4, 3], dstype=Data2D, id=idval)

    mdl = s.create_model_component("const2d", "mdl")
    s.set_source(idval, mdl)

    out = tmp_path / "not-created.dat"
    with pytest.raises(AttributeError, match=emsg):
        if idval is None:
            s.save_delchi(str(out))
        else:
            s.save_delchi(idval, str(out))


def check_save_ascii2d(session, expected, out, savefunc, idval, kwargs, check_str):
    """Checking the output is tricky.

    If crates is in use with the AstroSession backend then the call
    will fail as we don't allow ASCII images to be written out. Why,
    I don't know. So catch this so we can see if it ever changes.

    The assumption is that if pycrates is available it is being used
    in the test. This is not 100% correct in general, but should be
    good enough for the test.

    """

    if session == AstroSession and io.backend.name == "crates":
        if idval is None:
            with pytest.raises(IOErr,
                               match="writing images in ASCII is not supported"):
                savefunc(str(out), **kwargs)
        else:
            with pytest.raises(IOErr,
                               match="writing images in ASCII is not supported"):
                savefunc(idval, str(out), **kwargs)

        return

    if idval is None:
        savefunc(str(out), **kwargs)
    else:
        savefunc(idval, str(out), **kwargs)

    check_str(out.read_text(), expected)


@pytest.mark.parametrize("session,kwargs,expected",
                         [(Session, {"comment": "!! "}, ["!! SOURCE", "7 11", ""]),
                          (AstroSession, {"ascii": True}, ["7", "11", ""])])
@pytest.mark.parametrize("idval", [None, "bob"])
def test_save_source_ascii_data2d(session, kwargs, expected, idval, tmp_path, skip_if_no_io, check_str):
    """Basic check it works

    The output depends on the backend which is a bit distressing, and
    it's just a bit different enough to make it hard to check both
    source and model with the same test.

    """

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.dataspace2d([2, 1], dstype=Data2D, id=idval)

    mdl = s.create_model_component("polynom2d", "mdl")
    mdl.cy1 = 2
    mdl.cx1 = 4
    s.set_source(idval, mdl)

    out = tmp_path / "created.dat"
    check_save_ascii2d(session, expected, out, s.save_source, idval, kwargs, check_str)


@pytest.mark.parametrize("session,kwargs,expected",
                         [(Session, {"comment": ""}, ["MODEL", "7 11", ""]),
                          (AstroSession, {"ascii": True}, ["7", "11", ""])])
@pytest.mark.parametrize("idval", [None, "bob"])
def test_save_model_ascii_data2d(session, kwargs, expected, idval, tmp_path, skip_if_no_io, check_str):
    """Basic check it works

    The output depends on the backend which is a bit distressing, and
    it's just a bit different enough to make it hard to check both
    source and model with the same test.

    """

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.dataspace2d([2, 1], dstype=Data2D, id=idval)

    mdl = s.create_model_component("polynom2d", "mdl")
    mdl.cy1 = 2
    mdl.cx1 = 4
    s.set_source(idval, mdl)

    out = tmp_path / "created.dat"
    check_save_ascii2d(session, expected, out, s.save_model, idval, kwargs, check_str)


@pytest.mark.parametrize("session,kwargs,expected",
                         [(Session, {}, ["#RESID", "-7 -11", ""]),
                          (AstroSession, {"ascii": True}, ["-7", "-11", ""])])
@pytest.mark.parametrize("idval", [None, "bob"])
def test_save_resid_ascii_data2d(session, kwargs, expected, idval, tmp_path, skip_if_no_io, check_str):
    """Basic check it works

    The output depends on the backend which is a bit distressing.
    """

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.dataspace2d([2, 1], dstype=Data2D, id=idval)

    mdl = s.create_model_component("polynom2d", "mdl")
    mdl.cy1 = 2
    mdl.cx1 = 4
    s.set_source(idval, mdl)

    out = tmp_path / "created.dat"
    check_save_ascii2d(session, expected, out, s.save_resid, idval, kwargs, check_str)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_load_template_model_error_mismatched_header(session, tmp_path):
    """Check error handling: col names and col values differ"""

    mfile = tmp_path / "model.dat"
    mfile.write_text("# XPAR YPAR MODELFLAG FILENAME\n" +
                     f"1 2 foo.dat\n")

    s = session()
    with pytest.raises(IOErr,
                       match="Expected 3 columns but found 4"):
        s.load_template_model("tmp", str(mfile))


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_load_template_model_error_no_pars(session, tmp_path):
    """Check error handling: no parameter column"""

    mfile = tmp_path / "model.dat"
    mfile.write_text("# MODELFLAG FILENAME\n" +
                     f"1 foo.dat\n")

    s = session()
    with pytest.raises(IOErr,
                       match="No parameter columns found in .*/model.dat"):
        s.load_template_model("tmp", str(mfile))


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_load_template_model_error_no_modelfile(session, tmp_path):
    """Check error handling: no modelfile"""

    mfile = tmp_path / "model.dat"
    mfile.write_text("# XPAR YPAR\n" +
                     f"1 2\n1 3\n")

    s = session()
    with pytest.raises(IOErr,
                       match="Required column 'filename' not found in .*/model.dat"):
        s.load_template_model("tmp", str(mfile))


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_load_template_model_error_no_modelflag(session, tmp_path):
    """Check error handling: no modelfile"""

    mfile = tmp_path / "model.dat"
    mfile.write_text("# XPAR YPAR FILENAME\n" +
                     f"1 2 a.dat\n1 3 b.dat\n")

    s = session()
    with pytest.raises(IOErr,
                       match="Required column 'modelflag' not found in .*/model.dat"):
        s.load_template_model("tmp", str(mfile))


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_load_template_model_error_1_col(session, tmp_path):
    """Check error handling: template file has < 2 columns"""

    # Note: we need the actual path in the template file
    foo12 = tmp_path / "foo12.dat"
    foo12.write_text("1\n2\n")

    mfile = tmp_path / "model.dat"
    mfile.write_text("# XPAR YPAR MODELFLAG FILENAME\n" +
                     f"1 2 1 {foo12}\n")

    s = session()
    with pytest.raises(IOErr,
                       match="Only found 1 column in .*.dat, need at least 2"):
        s.load_template_model("tmp", str(mfile))


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_load_template_model_error_multi_col(session, tmp_path):
    """Check error handling: template file has > 2 columns"""

    # Note: we need the actual path in the template file
    foo12 = tmp_path / "foo12.dat"
    foo12.write_text("1 2 3\n2 4 5\n")

    mfile = tmp_path / "model.dat"
    mfile.write_text("# XPAR YPAR MODELFLAG FILENAME\n" +
                     f"1 2 1 {foo12}\n")

    s = session()
    with pytest.raises(IOErr,
                       match="Expected 2 columns but found 3"):
        s.load_template_model("tmp", str(mfile))


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_method_unknown(session):
    """Check we error out"""

    s = session()
    with pytest.raises(ArgumentErr,
                       match="'bob' is not a valid method"):
        s.get_method("bob")


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_method_opt_unknown(session):
    """Check we error out"""

    s = session()
    # make sure we use a known method
    s.set_method("LEVMAR")
    with pytest.raises(ArgumentErr,
                       match="'bob' is not a valid option for method levmar"):
        s.get_method_opt("bob")


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_method_opt_known(session):
    """Check we get a result.

    This is a regression test so it will need updating if we update
    the levmar settings.

    """

    s = session()
    # make sure we use a known method
    s.set_method("LEVMAR")
    assert s.get_method_opt("factor") == pytest.approx(100.0)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_set_iter_method_opt_unknown(session):
    """Check we error out"""

    s = session()
    with pytest.raises(ArgumentErr,
                       match="'bob' is not a valid option for method none"):
        s.set_iter_method_opt("bob", True)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_stat_unknown(session):
    """Check we error out"""

    s = session()
    with pytest.raises(ArgumentErr,
                       match="'bob' is not a valid statistic"):
        s.get_stat("bob")


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_indep(session):
    """Somehow we were not calling the Session version"""

    s = session()
    s.load_arrays(1, [2, 3, 4], [9, 8, 7])
    got = s.get_indep()
    assert len(got) == 1
    assert isinstance(got, tuple)
    assert got[0] == pytest.approx([2, 3, 4])


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_dims(session):
    """Somehow we were not calling the Session version"""

    s = session()
    s.dataspace2d([2, 3], dstype=Data2DInt)
    got = s.get_dims()
    assert isinstance(got, tuple)
    assert got == (2, 3)


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("etype", ["stat", "sys"])
@pytest.mark.parametrize("idval", [None, "bob"])
def test_load_xxxerror(session, etype, idval, tmp_path, skip_if_no_io):
    """Check we are called

    Mainly to ensure we have coverage of the Session class.
    """

    statfile = tmp_path / "stats.dat"
    statfile.write_text("1.1\n2.2\n0.0\n")

    s = session()
    s.load_arrays(idval, [1, 2, 3], [4, 5, 6])

    load = getattr(s, f"load_{etype}error")
    get = getattr(s, f"get_{etype}error")

    if idval is None:
        load(str(statfile))
    else:
        load(idval, str(statfile))

    assert get(idval) == pytest.approx([1.1, 2.2, 0.0])


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("idval", [None, "bob"])
def test_fake_fixed(session, idval):
    """Basic check of fake"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(idval, [2, 5, 9], [12, 13, 14])
    mdl = s.create_model_component("polynom1d", "mdl")
    s.set_source(idval, mdl)
    mdl.c0 = 2
    mdl.c1 = 1

    # We use a custom method which does not add noise
    #
    def custom(x, rng=None):
        return x

    assert s.get_dep(idval) == pytest.approx([12, 13, 14])

    if idval is None:
        s.fake(method=custom)
    else:
        s.fake(idval, method=custom)

    assert s.get_dep(idval) == pytest.approx([4, 7, 11])


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("idval", [None, "bob"])
def test_fake_random(session, idval):
    """Basic check of fake"""

    s = session()
    s._add_model_types(sherpa.models.basic)
    s.set_rng(numpy.random.RandomState(735))

    s.load_arrays(idval, [2, 5, 9], [12, 13, 14])
    mdl = s.create_model_component("polynom1d", "mdl")
    s.set_source(idval, mdl)
    mdl.c0 = 2
    mdl.c1 = 1

    assert s.get_dep(idval) == pytest.approx([12, 13, 14])

    if idval is None:
        s.fake()
    else:
        s.fake(idval)

    # values will need changing if the random seed <-> generator changes
    #
    assert s.get_dep(idval) == pytest.approx([9, 4, 11])


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_save_filter_excluded_everything(session, tmp_path):
    """Check we error out"""

    s = session()
    s.load_arrays(1, [10, 15, 20], [1, 2, 3])

    # Hide the messages from ignore
    with SherpaVerbosity('WARN'):
        s.ignore()

    ofile = tmp_path / 'not-created.dat'
    with pytest.raises(DataErr,
                       match="mask excludes all data"):
        s.save_filter(str(ofile))


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_save_filter_noticed_everything(session, tmp_path):
    """Check we error out"""

    s = session()
    s.set_default_id("bob")
    s.load_arrays("bob", [10, 15, 20], [1, 2, 3])

    ofile = tmp_path / 'not-created.dat'
    with pytest.raises(DataErr,
                       match="data set 'bob' has no filter"):
        s.save_filter(str(ofile))


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_add_model_types_simple(session):
    """coverage check"""

    s = session()

    assert s.list_models() == []
    s._add_model_types(sherpa.models.basic)

    # This will need updating if models are added to basic
    assert len(s.list_models()) == 31


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_add_model_types_scalar(session):
    """coverage check"""

    s = session()
    s._add_model_types(sherpa.models.basic,
                       baselist=ArithmeticModel)

    # This will need updating if models are added to basic
    assert len(s.list_models()) == 31


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_model_component_sent_str(session):
    """Just to match test_get_model_component_sent_model"""

    s = session()
    s._add_model_types(sherpa.models.basic)
    mdl = s.create_model_component("gauss1d", "g1")

    assert s.get_model_component("g1") == mdl


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_model_component_sent_model(session):
    """coverage check"""

    s = session()
    s._add_model_types(sherpa.models.basic)
    mdl = s.create_model_component("gauss1d", "g1")

    assert s.get_model_component(mdl) == mdl


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_create_model_component_sent_model(session):
    """coverage check"""

    s = session()
    s._add_model_types(sherpa.models.basic)
    mdl = s.create_model_component("gauss1d", "g1")

    assert s.create_model_component(mdl) == mdl


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_create_model_component_sent_unknown(session):
    """coverage check"""

    s = session()
    s._add_model_types(sherpa.models.basic)
    with pytest.raises(ArgumentErr,
                       match="'notamodel1d' is not a valid model type"):
        s.create_model_component("notamodel1d", "g1")


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_reset_sent_model(session):
    """Check model behavior"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    mdl = s.create_model_component("const1d", "foo")
    mdl.c0 = 2

    # Manually change the model so that reset can change it
    #
    mdl.thawedpars = [5]

    assert mdl.c0.val == pytest.approx(5)

    s.reset(mdl)

    assert mdl.c0.val == pytest.approx(2)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_reset_sent_id(session):
    """Check behavior when sent an id"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    mdl1 = s.create_model_component("const1d", "mdl1")
    mdl1.c0 = 2

    mdl2 = s.create_model_component("gauss1d", "mdl2")
    mdl2.pos = 10

    s.set_source(1, mdl1)
    s.set_source(2, mdl2)

    # Manually change the model so that reset can change it
    #
    mdl1.thawedpars = [5]
    mdl2.thawedpars = [5, 5, 5]

    s.reset(id=2)

    assert mdl1.c0.val == pytest.approx(5)
    assert mdl2.pos.val == pytest.approx(10)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_reset_sent_nothing(session):
    """Check behavior when sent nothing"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    mdl1 = s.create_model_component("const1d", "mdl1")
    mdl1.c0 = 2

    mdl2 = s.create_model_component("gauss1d", "mdl2")
    mdl2.pos = 10

    s.set_source(1, mdl1)
    s.set_source(2, mdl2)

    # Manually change the model so that reset can change it
    #
    mdl1.thawedpars = [5]
    mdl2.thawedpars = [5, 5, 5]

    s.reset()

    assert mdl1.c0.val == pytest.approx(2)
    assert mdl2.pos.val == pytest.approx(10)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_model_type_basic(session):
    """Check behavior"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    mdl1 = s.create_model_component("const1d", "c1")
    mdl2 = s.create_model_component("gauss1d", "g1")

    assert s.get_model_type(mdl1) == "const1d"
    assert s.get_model_type(mdl2) == "gauss1d"

    assert s.get_model_type(-mdl1) == "unaryopmodel"
    assert s.get_model_type(mdl1 + mdl2) == "binaryopmodel"


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_model_pars_str(session):
    """Check behavior"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    mdl = s.create_model_component("gauss1d", "g1")
    assert s.get_model_pars("g1") == ["fwhm", "pos", "ampl"]


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_model_pars_model(session):
    """Check behavior"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    mdl = s.create_model_component("gauss1d", "g1")
    assert s.get_model_pars(mdl) == ["fwhm", "pos", "ampl"]


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_num_par(session):
    """Check behavior"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    mdl1 = s.create_model_component("const1d", "mdl1")
    mdl2 = s.create_model_component("gauss1d", "mdl2")

    s.set_source(1, mdl1)
    s.set_source(2, mdl2)

    assert s.get_num_par() == 1
    assert s.get_num_par(2) == 3


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_model_component_checks_model_name(session):
    """Check corner case"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    with pytest.raises(IdentifierErr,
                       match="'gauss1d' is a model type, not allowed as a model name"):
        s.create_model_component("gauss1d", "gauss1d")


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_model_component_checks_builtin_name(session):
    """Check corner case"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    with pytest.raises(IdentifierErr,
                       match="'exit' is reserved for the native Python function"):
        s.create_model_component("gauss1d", "exit")


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_conf_opt_none(session):
    """Check corner case

    We could also check get_covar_opt/get_proj_opt but not worth it.
    """

    s = session()
    out = s.get_conf_opt()
    assert isinstance(out, dict)
    assert len(out) > 0

    # This will need to be updated if the setting changes
    assert not out["fast"]


def check_text_output(path, header, coldata):
    """The output format for the various backends is subtly different.

    path is the pathlib object representing the file
    """

    # Use the same logic as test_astro_ui_unit.py.
    #
    from sherpa.astro import io
    if io.backend.name == "crates":
        expected = f"#TEXT/SIMPLE\n# {header}\n"
    elif io.backend.name == "pyfits":
        expected = f"# {header}\n"
    else:
        assert False, f"Unknown I/O backend: {io.backend.name}"

    expected += coldata
    assert path.read_text() == expected


@requires_fits
@pytest.mark.parametrize("idval", [None, 1])
def test_save_filter_astro(idval, tmp_path):
    """Compare to test_session::test_save_filter_simple"""

    s = AstroSession()
    s.load_arrays(idval, [1, 2, 3], [12, 15, 2])
    s.get_data(idval).mask = [True, False, True]

    outfile = tmp_path / "filter.dat"
    if idval is None:
        s.save_filter(str(outfile))
    else:
        s.save_filter(idval, str(outfile))

    check_text_output(outfile, "X FILTER", "1 1\n2 0\n3 1\n")


@requires_fits
@pytest.mark.parametrize("idval", [None, 1])
def test_save_grouping(idval, tmp_path):
    """Check corner case"""

    s = AstroSession()
    s.load_arrays(idval, [1, 2, 3], [12, 15, 2], DataPHA)
    s.set_grouping(idval, [1, -1, 1])

    outfile = tmp_path / "grp.dat"
    if idval is None:
        s.save_grouping(str(outfile))
    else:
        s.save_grouping(idval, str(outfile))

    check_text_output(outfile, "CHANNEL GROUPS", "1 1\n2 -1\n3 1\n")


@requires_fits
@pytest.mark.parametrize("idval", [None, 1])
def test_save_quality(idval, tmp_path):
    """Check corner case"""

    s = AstroSession()
    s.load_arrays(idval, [1, 2, 3], [12, 15, 2], DataPHA)
    s.set_quality(idval, [0, 2, 0])

    outfile = tmp_path / "qual.dat"
    if idval is None:
        s.save_quality(str(outfile))
    else:
        s.save_quality(idval, str(outfile))

    check_text_output(outfile, "CHANNEL QUALITY", "1 0\n2 2\n3 0\n")


@requires_fits
@pytest.mark.parametrize("idval", [None, 1])
def test_save_table(idval, tmp_path):
    """Check corner case"""

    s = AstroSession()
    s.load_arrays(idval, [1, 2, 3], [12, 15, 2], DataPHA)
    s.set_quality(idval, [0, 2, 0])

    outfile = tmp_path / "table.dat"
    if idval is None:
        s.save_table(str(outfile), ascii=True)
    else:
        s.save_table(idval, str(outfile), ascii=True)

    check_text_output(outfile, "CHANNEL COUNTS QUALITY", "1 12 0\n2 15 2\n3 2 0\n")


def test_save_all_not_in_session():
    """This is to note that if we ever change this then test_save_all_empty_stdout can be parametrized"""

    s = Session()
    assert not hasattr(s, "save_all")


def test_save_all_empty_stdout(caplog, capfd):
    """Check corner case.

    Probably a good idea to check we get screen output.

    For some reason Session does not have save_all.

    """

    s = AstroSession()
    s.save_all()
    assert len(caplog.record_tuples) == 0

    out, err = capfd.readouterr()
    assert err == ""

    # do not check the full output
    assert out.startswith("import numpy\nfrom sherpa.astro.ui import *\n\n")


def test_calc_data_sum_of_bkg():
    """Check we recognize the background"""

    s = AstroSession()
    s.load_arrays(1, [1, 2, 3], [4, 5, 6], DataPHA)
    s.set_bkg(DataPHA("x", [1, 2, 3], [9, 2, 1]))

    assert s.calc_data_sum(2, 3, bkg_id=1) == 3


def test_check_get_source_and_model_with_background():
    """Check the background is not included in get_source but is for get_model

    See issue #1686
    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)
    s.load_arrays(1, [1, 2, 3], [4, 5, 6], DataPHA)
    s.set_bkg(DataPHA("x", [1, 2, 3], [9, 2, 1]))

    s.set_exposure(100)
    s.set_exposure(200, bkg_id=1)

    s.set_backscal(2)
    s.set_backscal(8, bkg_id=1)

    egrid = numpy.asarray([0.1, 0.2, 0.4, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)

    # Need a response for get_model.
    #
    s.set_rmf(rmf)
    s.set_rmf(rmf, bkg_id=1)

    g1 = s.create_model_component("gauss1d", "g1")
    g1.ampl = 10000

    g2 = s.create_model_component("polynom1d", "g2")
    g2.c1 = 2

    s.set_source(g1)
    s.set_bkg_source(g2)

    assert s.get_source().name == "gauss1d.g1"
    assert s.get_bkg_source().name == "polynom1d.g2"

    # Note: DJB has repeatedly got concerned about how to scale the
    # background to match the source region. As noted in #888 the
    # source expressions (g1 and g2 in this case) evaluate to rates,
    # so there is no need to account for the difference in the
    # exposure times for the source and background datasets. However,
    # there *is* a need to account for the aperture differences, since
    # the models are not normalized by area. Hence the inclusion of
    # the backscal ratio (i.e. 2 / 8) in the model expression.
    #
    assert s.get_model().name == "apply_rmf(100.0 * (gauss1d.g1 + 0.25 * polynom1d.g2))"
    assert s.get_bkg_model().name == "apply_rmf(200.0 * polynom1d.g2)"


def test_calc_model_sum_of_bkg():
    """Check we recognize the background"""

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)
    s.load_arrays(1, [1, 2, 3], [4, 5, 6], DataPHA)
    s.set_bkg(DataPHA("x", [1, 2, 3], [9, 2, 1]))

    egrid = numpy.asarray([0.1, 0.2, 0.4, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)

    s.set_rmf(rmf)
    s.set_rmf(rmf, bkg_id=1)

    s.set_analysis("channel")

    # Pick the source model that is very different to the background
    # model.
    #
    g1 = s.create_model_component("gauss1d", "g1")
    g1.ampl = 10000

    g2 = s.create_model_component("polynom1d", "g2")
    g2.c1 = 2

    s.set_source(g1)
    s.set_bkg_source(g2)

    # The model is y = 1 + 2 * x, so for channels 2:3 it evaluates to
    # [5, 7]. We then have to filter in the oddness that is PHA files,
    # which in this case is the fact we have a RMF file applied. So
    # this is a regression test as I currently can not work out why we
    # get this value (note that changing egrid changes this result, so
    # the bin widths are important).
    #
    assert s.calc_model_sum(2, 3, bkg_id=1) == pytest.approx(1.2)

    # safety check
    assert s.calc_model_sum(2, 3) == pytest.approx(5954.871098716832)


def test_calc_source_sum_of_bkg():
    """Check we recognize the background"""

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)
    s.load_arrays(1, [1, 2, 3], [4, 5, 6], DataPHA)
    s.set_bkg(DataPHA("x", [1, 2, 3], [9, 2, 1]))

    egrid = numpy.asarray([0.1, 0.2, 0.4, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)

    s.set_rmf(rmf)
    s.set_rmf(rmf, bkg_id=1)

    s.set_analysis("channel")

    # Pick the source model that is very different to the background
    # model.
    #
    g1 = s.create_model_component("gauss1d", "g1")
    g1.ampl = 10000

    g2 = s.create_model_component("polynom1d", "g2")
    g2.c1 = 2

    s.set_source(g1)
    s.set_bkg_source(g2)

    # This is even more confusing than test_calc_model_sum_of_bkg.
    # Why is this returning 0?
    #
    assert s.calc_source_sum(2, 3, bkg_id=1) == pytest.approx(0)

    # Also what is going on here?
    assert s.calc_source_sum(2, 3) == pytest.approx(0)


def test_fit_bkg_no_models():
    """Corner case"""

    s = AstroSession()
    s.load_arrays(1, [1, 2, 3], [1, 2, 3], DataPHA)
    with pytest.raises(IdentifierErr,
                       match="^model stack is empty$"):
        s.fit_bkg()


def test_bkg_fit_data_with_no_model():
    """Coverage check

    This was aimed at ensuring a particular error condition was checked,
    but it turns out to be hard to trigger, so we just leave this as is
    as it is an edge case.

    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(1, [1, 2, 3], [1, 2, 3], DataPHA)
    s.set_bkg(DataPHA("ex", [1, 2, 3], [2, 0, 1]))
    s.set_bkg(DataPHA("ex", [1, 2, 3], [2, 0, 1]), bkg_id=2)

    egrid = numpy.asarray([0.1, 0.2, 0.4, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)

    s.set_rmf(rmf)
    s.set_rmf(rmf, bkg_id=1)

    g1 = s.create_model_component("gauss1d", "g1")
    s.set_source(g1)

    g2 = s.create_model_component("gauss1d", "g2")
    s.set_bkg_source(g2)

    with pytest.raises(ModelErr,
                       match="background model 2 for data set 1 has not been set"):
        s.get_stat_info()


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_do_not_fold_all_models(session, tmp_path, skip_if_no_io):
    """Corner case of _add_convolution_models

    We add a check that the fold model has been called, but it is
    rather "ugly" as there's no "nice" way to find this out, so we
    call the model with filtered data and check whether it errors out
    when it should.

    """

    # In #1662 when the test was added the comment character was '!' and
    # not '#', but that causes issues with the pycrates backend, so it
    # was changed as it is not an important part of the test.
    #
    tblfile = tmp_path / "tbl.dat"
    tblfile.write_text("#Y\n10\n5\n2\n1\n")

    s = session()
    s.load_table_model("tbl1", str(tblfile), comment="#", ncols=1)
    s.load_table_model("tbl2", str(tblfile), comment="#", ncols=1)

    tbl1 = s.get_model_component("tbl1")
    tbl2 = s.get_model_component("tbl2")

    def does_error(mdl):
        with pytest.raises(ModelErr, match="Mismatch between table model and data, "):
            mdl([1, 2, 3])

    does_error(tbl1)
    does_error(tbl2)

    s.load_arrays(1, [1, 2, 3, 4], [9, 3, 3, 2])
    s.load_arrays("x", [1, 2, 3, 4], [9, 3, 3, 2])
    s.set_source("x", "tbl2")

    # Filter the data so we know the data has been folded.
    # This only works if the table model had no X axis.
    #
    s.ignore_id("x", lo=4)

    does_error(tbl1)
    does_error(tbl2)

    # This call is the main thing being tested; the does_error/tbl2
    # checks are done as extras.
    #
    check = s.get_model("x")
    assert check == tbl2

    does_error(tbl1)
    assert tbl2([1, 2, 3]) == pytest.approx([10, 5, 2])


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_delete_psf_not_there(session):
    """delete_psf does nothing if no PSF exist"""

    s = session()
    assert s.list_psf_ids() == []
    s.delete_psf()
    assert s.list_psf_ids() == []


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_delete_psf_there(session):
    """delete_psf removes the PSF"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.dataspace1d(1, 10)
    s.load_psf("bob", "gauss1d.mdl")
    s.set_psf("bob")

    s.dataspace1d(1, 10, id=2)
    s.load_psf("fred", "gauss1d.mdl")
    s.set_psf(2, "fred")

    assert s.list_psf_ids() == [1, 2]
    s.delete_psf(2)
    assert s.list_psf_ids() == [1]


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_fit_ignores_repeated_otherid(session, caplog):

    s = session()
    s._add_model_types(sherpa.models.basic)

    # Hmm, using the same id but as a string is a bit dangerous,
    # so let's see what happens.
    #
    s.load_arrays(1, [1, 2], [2, 5])
    s.load_arrays("1", [10, 20], [7, 8])

    s.set_source(1, "const1d.con")
    s.set_source("1", "con")

    s.set_stat("leastsq")
    with SherpaVerbosity("WARN"):
        s.fit(1, "1", 1, "1")

    res = s.get_fit_results()
    assert res.datasets == (1, "1")
    assert res.numpoints == 4
    assert res.dof == 3
    assert res.parnames == ("con.c0", )
    assert len(res.parvals) == 1
    assert res.parvals[0] == pytest.approx(5.5)

    # We can also check a corner case of calc_stat_info, for
    # handling multiple datasets.
    #
    assert len(caplog.record_tuples) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        s.calc_stat_info()

    assert len(caplog.record_tuples) == 1
    loc, lvl, msg = caplog.record_tuples[0]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO

    toks = msg.split("\n\n")
    assert len(toks) == 3

    # This is a case where using the same id as an integer and string
    # makes it hard to understand the results!
    #
    assert toks[0].startswith("Dataset               = 1\n")
    assert "Fit statistic value   = 12.5\n" in toks[0]
    assert toks[0].endswith("Degrees of freedom    = 1")

    assert toks[1].startswith("Dataset               = 1\n")
    assert "Fit statistic value   = 8.5\n" in toks[1]
    assert toks[1].endswith("Degrees of freedom    = 1")

    assert toks[2].startswith("Datasets              = [1, '1']\n")
    assert "Fit statistic value   = 21\n" in toks[2]
    assert toks[2].endswith("Degrees of freedom    = 3")

    # We can also check out calc_chisqr
    #
    d1 = [12.25, 0.25]
    d2 = [2.25, 6.25]
    assert s.calc_chisqr() == pytest.approx(d1 + d2)  # array concatenation not value addition
    assert s.calc_chisqr(1) == pytest.approx(d1)
    assert s.calc_chisqr("1") == pytest.approx(d2)
    assert s.calc_chisqr(1, "1") == pytest.approx(d1 + d2)
    assert s.calc_chisqr("1", 1, 1, "1") == pytest.approx(d2 + d1)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_draws_wants_errors(session):

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(1, [1, 2], [2, 5])
    s.set_source(1, "const1d.foo")
    s.fit()

    with pytest.raises(SessionErr,
                       match="covariance has not been performed"):
        s.get_draws(niter=1)


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("method", ["proj", "unc"])
def test_get_int_xxx_otherids(method, session):
    """Corner case.

    This does not need a plotting backend.
    """

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(1, [1, 2, 3, 4], [2, 5, 7, 12])
    s.load_arrays(2, [10, 20, 35, 40], [7, 8, 6, 12])

    s.set_staterror(1, [1, 1, 1, 1])
    s.set_staterror(2, [1, 1, 1, 1])

    s.set_source(1, "const1d.con")
    s.set_source(2, "con")

    with SherpaVerbosity("WARN"):
        s.fit(1, 2)

    # A minimal check of the returned values.
    #
    get = getattr(s, f"get_int_{method}")
    plot = get(id=2, par="con.c0", otherids=[1], recalc=True,
               min=5, max=8, nloop=5)

    assert plot.min == pytest.approx(5)
    assert plot.max == pytest.approx(8)
    assert plot.nloop == 5
    assert plot.x == pytest.approx([5, 5.75, 6.5, 7.25, 8])
    assert plot.y == pytest.approx([125, 101, 86, 80, 83])


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("method", ["proj", "unc"])
def test_get_reg_xxx_otherids(method, session):
    """Corner case.

    This does not need a plotting backend.
    """

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.load_arrays(1, [1, 2, 3, 4], [2, 5, 7, 12])
    s.load_arrays(2, [10, 20, 35, 40], [7, 8, 6, 12])

    s.set_staterror(1, [1, 1, 1, 1])
    s.set_staterror(2, [1, 1, 1, 1])

    s.set_source(1, "polynom1d.mdl")
    s.set_source(2, "mdl")

    mdl = s.get_model_component("mdl")
    mdl.c1.thaw()

    with SherpaVerbosity("WARN"):
        s.fit(1, 2)

    # A minimal check of the returned values.
    #
    get = getattr(s, f"get_reg_{method}")
    plot = get("mdl.c0", "mdl.c1", id=2, otherids=[1], recalc=True,
               min=[4, 0], max=[8, 0.2], nloop=[5, 5])

    assert plot.parval0 == pytest.approx(mdl.c0.val)
    assert plot.parval1 == pytest.approx(mdl.c1.val)
    assert plot.y.min() == pytest.approx(66.35)
    assert plot.y.max() == pytest.approx(184.8)


@pytest.mark.parametrize("etype", ["stat", "sys"])
def test_set_xxxerror_bkg(etype):

    s = AstroSession()
    s.load_arrays(1, [1, 2, 3], [2, 3, 2], DataPHA)
    s.set_bkg(DataPHA("bkg", [1, 2, 3], [1, 0, 2]))

    setfunc = getattr(s, f"set_{etype}error")
    setfunc([2, 4, 3], bkg_id=1)

    getfunc = getattr(s, f"get_{etype}error")
    assert getfunc(bkg_id=1) == pytest.approx([2, 4, 3])


@pytest.mark.parametrize("idval", [None,
                                   pytest.param(1, marks=pytest.mark.xfail),
                                   pytest.param("bob", marks=pytest.mark.xfail)])
@pytest.mark.parametrize("field,startval",
                         [("exposure", 23.0),
                          ("backscal", 0.23),
                          ("areascal", 0.23)])
def test_set_pha_value_none(idval, field, startval):
    """Corner case.

    Unfortunately, the API makes it hard to set the value to None:
    it only works for the default id. The errors when the idval is set
    differ (if an integer then the id gets set as the value, if a
    string there's a string conversion error).

    """

    s = AstroSession()
    s.load_arrays(idval, [1, 2, 3], [2, 4, 3], DataPHA)

    setfunc = getattr(s, f"set_{field}")
    getfunc = getattr(s, f"get_{field}")

    setfunc(idval, startval)
    assert getfunc(idval) == pytest.approx(startval)

    setfunc(idval, None)
    assert getfunc(idval) is None


@requires_fits
def test_pack_pha():
    """Very basic check: does it work?"""

    s = AstroSession()
    s.load_arrays(1, [1, 2, 3], [4, 5, 2], DataPHA)

    # We do the absolute-minimum as a check
    assert s.pack_pha() is not None


@requires_fits
def test_pack_table():
    """Very basic check: does it work?"""

    s = AstroSession()
    s.load_arrays(1, [1, 2, 3], [4, 5, 2], DataPHA)

    # We do the absolute-minimum as a check
    assert s.pack_table() is not None


@requires_fits
def test_pack_image():
    """Very basic check: does it work?"""

    s = AstroSession()
    s.dataspace2d([2, 3])

    # We do the absolute-minimum as a check
    assert s.pack_image() is not None


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_show_kernel(session, idval):

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.dataspace1d(1, 10, id=idval)
    s.load_psf("bob", "gauss1d.mdl")
    if idval is None:
        s.set_psf("bob")
    else:
        s.set_psf(idval, "bob")

    out = StringIO()
    if idval is None:
        s.show_kernel(outfile=out)
    else:
        s.show_kernel(idval, outfile=out)

    idstr = "1" if idval is None else idval

    toks = out.getvalue().split("\n")
    assert toks[0] == f"PSF Kernel: {idstr}"
    assert toks[1] == "psfmodel.bob"
    assert len(toks) == 12


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_show_kernel_multi(session):

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.dataspace1d(1, 10, id=1)
    s.dataspace1d(1, 10, id="twenty")
    s.load_psf("bob", "gauss1d.mdl")
    s.load_psf("fred", "const1d.mdl2")

    bob = s.get_model_component("bob")
    fred = s.get_model_component("fred")
    s.set_psf(bob)
    s.set_psf("twenty", fred)

    out = StringIO()
    s.show_kernel(outfile=out)

    toks = out.getvalue().split("\n\n")
    assert toks[0].startswith("PSF Kernel: 1\npsfmodel.bob\n")
    assert toks[1].startswith("PSF Kernel: twenty\npsfmodel.fred\n")
    assert toks[2] == "\n"
    assert len(toks) == 3


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_load_conv_model_instance(session):
    """Corner case"""

    s = session()
    s._add_model_types(sherpa.models.basic)

    ngl = s.create_model_component("normgauss1d", "ngl")
    s.list_model_components() == ["ngl"]

    s.load_conv("bobby", ngl)
    s.list_model_components() == ["ngl", "bobby"]

    got = s.get_model_component("bobby")
    assert isinstance(got, ConvolutionKernel)
    assert got.name == "convolutionkernel.bobby"
    assert got.kernel == ngl


def check_fit_results(s, datasets, parnames, parvals, istatval, statval, numpoints, dof):
    """Check the fit results"""

    fres = s.get_fit_results()
    assert fres.datasets == datasets
    assert fres.succeeded
    assert fres.parnames == parnames
    assert fres.parvals == pytest.approx(parvals)
    assert fres.istatval == pytest.approx(istatval)
    assert fres.statval == pytest.approx(statval)
    assert fres.numpoints == numpoints
    assert fres.dof == dof


@pytest.mark.parametrize("id1", [1, 2])
@pytest.mark.parametrize("session", [Session, AstroSession])
def test_fit_data1d_mix_data(id1, session):
    """Check what happens with a mix of datasets with/without data.

    This was suggested by #1721 - I have not done a rigorous check to
    see if we have something like this already.

    """

    s = session()
    s._add_model_types(sherpa.models.basic)

    s.set_stat("leastsq")
    s.set_method("simplex")

    s.load_arrays(1, [1, 2], [1, 2])
    s.load_arrays(2, [1, 2], [1, 2])
    s.load_arrays(3, [10, 20], [0.7, 2.1])

    mdl = s.create_model_component("scale1d", "mdl")
    s.set_source(id1, mdl)
    s.set_source(3, mdl)

    s.fit()
    check_fit_results(s, (id1, 3), ("mdl.c0", ),
                      [1.45], 2.3, 1.49, 4, 3)


@pytest.mark.parametrize("id1", [1, 2])
def test_fit_datapha_mix_data(id1):
    """Check what happens with a mix of datasets with/without data.

    See also test_fit_data1d_mix_data

    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.set_stat("leastsq")
    s.set_method("simplex")

    s.load_arrays(1, [1, 2], [2, 4], DataPHA)
    s.load_arrays(2, [1, 2], [2, 4], DataPHA)
    s.load_arrays(3, [1, 2], [3, 2], DataPHA)

    # We allow ARF-only analysis, so test it. The fit seems to ignore
    # the bin widths in this case.
    #
    egrid = numpy.asarray([0.2, 0.5, 0.9])
    arf = create_arf(egrid[:-1], egrid[1:])
    s.set_arf(1, arf)
    s.set_arf(2, arf)
    s.set_arf(3, arf)

    mdl = s.create_model_component("scale1d", "mdl")
    s.set_source(id1, mdl)
    s.set_source(3, mdl)

    s.fit()
    check_fit_results(s, (id1, 3), ("mdl.c0", ),
                      [2.75], 15, 2.75, 4, 3)


@pytest.mark.parametrize("id1", [1, 2])
def test_fit_datapha_mix_data_with_bkg(id1):
    """Check what happens with a mix of datasets with/without data.

    See also test_fit_datapha_mix_data and test_fit_datapha_mix_data_only_bkg

    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.set_stat("leastsq")
    s.set_method("simplex")

    s.load_arrays(1, [1, 2], [2, 4], DataPHA)
    s.load_arrays(2, [1, 2], [2, 4], DataPHA)
    s.load_arrays(3, [1, 2], [3, 2], DataPHA)

    s.set_bkg(1, DataPHA("b1", [1, 2], [1, 0]))
    s.set_bkg(2, DataPHA("b2", [1, 2], [1, 0]))

    # We allow ARF-only analysis, so test it. The fit seems to ignore
    # the bin widths in this case.
    #
    egrid = numpy.asarray([0.2, 0.5, 0.9])
    arf = create_arf(egrid[:-1], egrid[1:])
    s.set_arf(1, arf)
    s.set_arf(2, arf)
    s.set_arf(3, arf)

    smdl = s.create_model_component("scale1d", "smdl")
    s.set_source(id1, smdl)
    s.set_source(3, smdl)

    bmdl = s.create_model_component("scale1d", "bmdl")
    s.set_bkg_source(id1, bmdl)

    s.fit()
    check_fit_results(s, (id1, 3), ("smdl.c0", "bmdl.c0"),
                      [2.5, 0.5], 10, 3, 6, 4)


@pytest.mark.parametrize("id1", [1, 2])
def test_fit_datapha_mix_data_only_bkg(id1):
    """Check what happens with a mix of datasets with/without data.

    See also test_fit_datapha_mix_data_with_bkg

    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.set_stat("leastsq")
    s.set_method("simplex")

    s.load_arrays(1, [1, 2], [2, 4], DataPHA)
    s.load_arrays(2, [1, 2], [2, 4], DataPHA)
    s.load_arrays(3, [1, 2], [3, 2], DataPHA)
    s.load_arrays(4, [1, 2], [10, 20], DataPHA)

    s.set_bkg(1, DataPHA("b1", [1, 2], [1, 0]))
    s.set_bkg(2, DataPHA("b2", [1, 2], [1, 0]))
    s.set_bkg(4, DataPHA("b4", [1, 2], [1, 0]))

    # We allow ARF-only analysis, so test it. The fit seems to ignore
    # the bin widths in this case.
    #
    egrid = numpy.asarray([0.2, 0.5, 0.9])
    arf = create_arf(egrid[:-1], egrid[1:])
    s.set_arf(1, arf)
    s.set_arf(2, arf)
    s.set_arf(3, arf)

    # Interesting that we need to do this. Is it a bug? Perhaps not
    # technically, but we no longer need to do this, but leave in to
    # compare to test_fit_datapha_mix_data_only_bkg_no_response.
    #
    s.set_arf(1, arf, bkg_id=1)
    s.set_arf(2, arf, bkg_id=1)

    smdl = s.create_model_component("scale1d", "smdl")
    s.set_source(id1, smdl)
    s.set_source(3, smdl)

    bmdl = s.create_model_component("scale1d", "bmdl")
    s.set_bkg_source(id1, bmdl)

    s.fit_bkg()
    check_fit_results(s, (id1, ), ("bmdl.c0", ),
                      [0.5], 1, 0.5, 2, 1)


@pytest.mark.parametrize("id1", [1, 2])
def test_fit_datapha_mix_data_only_bkg_no_response(id1):
    """Do we grab the response from the source data?

    See also test_fit_datapha_mix_data_only_bkg - this is essentially
    the same but does not set the background responses. Should this
    work or not (ie this is a regression test)?

    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.set_stat("leastsq")
    s.set_method("simplex")

    s.load_arrays(1, [1, 2], [2, 4], DataPHA)
    s.load_arrays(2, [1, 2], [2, 4], DataPHA)
    s.load_arrays(3, [1, 2], [3, 2], DataPHA)
    s.load_arrays(4, [1, 2], [10, 20], DataPHA)

    s.set_bkg(1, DataPHA("b1", [1, 2], [1, 0]))
    s.set_bkg(2, DataPHA("b2", [1, 2], [1, 0]))
    s.set_bkg(4, DataPHA("b4", [1, 2], [1, 0]))

    # We allow ARF-only analysis, so test it. The fit seems to ignore
    # the bin widths in this case.
    #
    egrid = numpy.asarray([0.2, 0.5, 0.9])
    arf = create_arf(egrid[:-1], egrid[1:])
    s.set_arf(1, arf)
    s.set_arf(2, arf)
    s.set_arf(3, arf)

    smdl = s.create_model_component("scale1d", "smdl")
    s.set_source(id1, smdl)
    s.set_source(3, smdl)

    bmdl = s.create_model_component("scale1d", "bmdl")
    s.set_bkg_source(id1, bmdl)

    s.fit_bkg()
    check_fit_results(s, (id1, ), ("bmdl.c0", ),
                      [0.5], 1, 0.5, 2, 1)


def test_fit_datapha_mix_data_with_bkg_sensible():
    """A more-standard set of backgrounds

    Each source has a single background.
    """

    s = AstroSession()
    s._add_model_types(sherpa.models.basic)

    s.set_stat("leastsq")
    s.set_method("simplex")

    s.load_arrays(1, [1, 2], [2, 4], DataPHA)
    s.load_arrays(2, [1, 2], [5, 3], DataPHA)
    s.load_arrays(3, [1, 2], [3, 2], DataPHA)

    s.set_bkg(1, DataPHA("b1", [1, 2], [1, 0]))
    s.set_bkg(2, DataPHA("b2", [1, 2], [0, 1]))
    s.set_bkg(3, DataPHA("b2", [1, 2], [1, 1]))

    # We allow ARF-only analysis, so test it. The fit seems to ignore
    # the bin widths in this case.
    #
    egrid = numpy.asarray([0.2, 0.5, 0.9])
    arf = create_arf(egrid[:-1], egrid[1:])
    s.set_arf(1, arf)
    s.set_arf(2, arf)
    s.set_arf(3, arf)

    smdl = s.create_model_component("scale1d", "smdl")
    s.set_source(1, smdl)
    s.set_source(2, smdl)
    s.set_source(3, smdl)

    bmdl = s.create_model_component("scale1d", "bmdl")
    s.set_bkg_source(1, bmdl)
    s.set_bkg_source(2, bmdl)
    s.set_bkg_source(3, bmdl)

    s.fit()
    check_fit_results(s, (1, 2, 3), ("smdl.c0", "bmdl.c0"),
                      [2.5, 2/3], 17, 8 + 1/6, 12, 10)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_fit_works(session):
    """Check the default output.

    This will need updating when the behavior changes.
    """

    s = session()
    assert s.get_rng() is None


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_set_fit_works(session):
    """Basic check."""

    s = session()

    rng = numpy.random.RandomState()
    s.set_rng(rng)
    got = s.get_rng()
    assert isinstance(got, numpy.random.RandomState)

    s.set_rng(None)
    got = s.get_rng()
    assert got is None


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_set_fit_checks_arg(session):
    """Basic check."""

    s = session()
    with pytest.raises(ArgumentTypeErr,
                       match="^'rng' must be a Generator or None$"):
        s.set_rng(1234)
