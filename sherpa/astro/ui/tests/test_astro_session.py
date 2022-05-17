#
#  Copyright (C) 2016, 2018, 2020, 2021, 2022
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

import numpy

import pytest

from sherpa.astro.data import DataPHA
from sherpa.astro.ui.utils import Session as AstroSession
from sherpa.data import Data1D, Data1DInt
from sherpa.io import get_ascii_data
from sherpa.models import Const1D
import sherpa.models.basic
from sherpa.ui.utils import Session
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, DataErr, IdentifierErr
from sherpa.utils.testing import requires_data, requires_fits, requires_group


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


# Fix 476 - this should be in sherpa/ui/tests/test_session.py
@requires_group
def test_zero_division_calc_stat():
    ui = AstroSession()
    x = numpy.arange(100)
    y = numpy.zeros(100)
    ui.load_arrays(1, x, y, DataPHA)
    ui.group_counts(1, 100)
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
    s = Session()
    assert s.get_iter_method_name() == "none"


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("opt", ["none", "sigmarej"])
def test_set_iter_method_valid(session, opt):
    s = Session()
    s.set_iter_method(opt)
    assert s.get_iter_method_name() == opt


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_set_iter_method_unknown_string(session):
    s = Session()
    with pytest.raises(TypeError,
                       match="^not a method is not an iterative fitting method$"):
        s.set_iter_method("not a method")


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_set_iter_method_not_a_string(session):
    s = Session()
    with pytest.raises(ArgumentTypeErr,
                       match="^'meth' must be a string$"):
        s.set_iter_method(23)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_iter_method_opt_default(session):
    s = Session()
    assert s.get_iter_method_opt() == {}


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_iter_method_opt_sigmarej(session):
    s = Session()
    s.set_iter_method("sigmarej")
    out = s.get_iter_method_opt()

    keys = set(out.keys())
    assert keys == set(["grow", "lrej", "hrej", "maxiters"])


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_get_iter_method_opt_sigmarej_lrej(session):
    s = Session()
    s.set_iter_method("sigmarej")
    assert s.get_iter_method_opt("lrej") == pytest.approx(3)


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("opt,key", [("none", "lrej"), ("sigmarej", "fast")])
def test_get_iter_method_opt_unknown(session, opt, key):
    s = Session()
    s.set_iter_method(opt)
    with pytest.raises(ArgumentErr,
                       match=f"^'{key}' is not a valid option for method {opt}$"):
        s.get_iter_method_opt(key)


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_set_iter_method_opt_sigmarej_lrej(session):
    s = Session()
    s.set_iter_method("sigmarej")
    s.set_iter_method_opt("lrej", 5)
    assert s.get_iter_method_opt("lrej") == pytest.approx(5)


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("setting", ['chisqr', 'compmodel', 'compsource', 'data',
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
                                     'order',
                                     # "energy", "photon",  these are currently both valid for astro
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


@requires_fits  # only for AstroSession, but not worth being clever
@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_model_ascii(session, kwargs, idval, tmp_path):
    """Can we use save_model?"""

    s = session()
    outfile = tmp_path / "model.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_model)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "MODEL"]
    assert len(data) == 2
    assert data[0] == pytest.approx([2, 4, 7])
    assert data[1] == pytest.approx([4, 4, 8])


@requires_fits  # only for AstroSession, but not worth being clever
@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_source_ascii(session, kwargs, idval, tmp_path):
    """Can we use save_source?"""

    s = session()
    outfile = tmp_path / "source.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_source)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "SOURCE"]
    assert len(data) == 2
    assert data[0] == pytest.approx([2, 4, 7])
    assert data[1] == pytest.approx([4, 4, 8])


@requires_fits  # only for AstroSession, but not worth being clever
@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_resid_ascii(session, kwargs, idval, tmp_path):
    """Can we use save_resid?"""

    s = session()
    outfile = tmp_path / "resid.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_resid)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "RESID"]
    assert len(data) == 2
    assert data[0] == pytest.approx([2, 4, 7])
    assert data[1] == pytest.approx([-1, 0, -3])


@requires_fits  # only for AstroSession, but not worth being clever
@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_delchi_ascii(session, kwargs, idval, tmp_path):
    """Can we use save_delchi?"""

    s = session()
    outfile = tmp_path / "resid.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_delchi)

    names, data = get_data(outfile, ncols=2)
    assert names == ["X", "DELCHI"]
    assert len(data) == 2
    assert data[0] == pytest.approx([2, 4, 7])
    assert data[1] == pytest.approx([-1, 0, -0.5])


@requires_fits  # only for AstroSession, but not worth being clever
@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_data_ascii(session, kwargs, idval, tmp_path):
    """Can we use save_data?"""

    s = session()
    outfile = tmp_path / "data.dat"
    save_ascii_file(s, kwargs, idval, outfile, s.save_data)

    print(outfile.read_text())

    names, data = get_data(outfile, ncols=4)
    assert names == ["XLO", "XHI", "Y", "STATERROR"]
    assert len(data) == 4
    assert data[0] == pytest.approx([1, 3, 5])
    assert data[1] == pytest.approx([3, 5, 9])
    assert data[2] == pytest.approx([3, 4, 5])
    assert data[3] == pytest.approx([1, 2, 6])


@requires_fits  # only for AstroSession, but not worth being clever
@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_staterror_ascii(session, kwargs, idval, tmp_path):
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


@requires_fits  # only for AstroSession, but not worth being clever
@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_syserror_ascii(session, kwargs, idval, tmp_path):
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


@requires_fits  # only for AstroSession, but not worth being clever
@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_syserror_no_data_ascii(session, kwargs, idval, tmp_path):
    """Check save_syserror needs a systematic-error column."""

    s = session()
    outfile = tmp_path / "errs.dat"

    mid = 1 if idval is None else idval
    with pytest.raises(DataErr,
                       match=f"^data set '{mid}' does not specify systematic errors$"):
        save_ascii_file(s, kwargs, idval, outfile, s.save_syserror)


@requires_fits  # only for AstroSession, but not worth being clever
@pytest.mark.parametrize("session,kwargs",
                         [(Session, {}), (AstroSession, {"ascii": True})])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_save_error_ascii(session, kwargs, idval, tmp_path):
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
