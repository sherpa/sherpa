#
#  Copyright (C) 2023
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

"""Can we write out XSPEC table models

This can be done without the need for XSPEC. However, in
order to check the models are sensible we need XSPEC to
read them in.

This is not a complete test of the functionality of table models, as
that is assumed to be handled by XSPEC. This is to check we can write
out the correct data and that it is handled correctly, but, for
instance, we assume that if something works for additive models then
it's likely to work for multiplicative models as well (we have other
table-model tests that provide some of those checks).

"""

import numpy as np

import pytest

from sherpa.astro.io import xstable
from sherpa.utils.err import IOErr
from sherpa.utils.testing import requires_fits, requires_xspec


@pytest.mark.parametrize("elo,ehi,msg",
                         [ ([], [],
                            "egrid_lo/hi can not be empty"),
                           ([1, 2], [1, 2, 3],
                            "egrid_lo/hi do not match: 2 vs 3"),
                           ([1, 2, 3], [1, 2],
                            "egrid_lo/hi do not match: 3 vs 2"),
                           ([1, 2], [],
                            "egrid_lo/hi do not match: 2 vs 0"),
                           ([], [1, 2],
                            "egrid_lo/hi do not match: 0 vs 2"),
                           ([1, 2, 4], [2, 3, 5],
                            "egrid_lo/hi are not consecutive"),
                           ([1, 3, 4], [2, 4, 5],
                            "egrid_lo/hi are not consecutive"),
                           ([1, 2, 3, 3], [2, 3, 4, 5],
                            "egrid_lo is not monotonically increasing"),
                           ([1, 2, 3, 4], [2, 2, 3, 4],
                            "egrid_hi is not monotonically increasing"),
                           ([4, 3, 2, 1], [5, 4, 3, 2],
                            "egrid_lo is not monotonically increasing")
                         ])
def test_check_valid_egrid(elo, ehi, msg):
    """Check we catch some error-grid errors"""

    p0 = xstable.Param("nop", 12, -1, 12, 12, values=[12])
    model = [1, 2, 3, 4, 5]
    with pytest.raises(ValueError, match=f"^{msg}$"):
        xstable.make_xstable_model("fake", elo, ehi,
                                   params=[p0], spectra=[model])


def test_check_single_spectrum_size():
    """Check we error out"""

    p0 = xstable.Param("nop", 12, -1, 12, 12, values=[12])
    model = [1, 2, 3, 4, 5]
    with pytest.raises(ValueError, match="^Spectrum should have 3 elements but has 5$"):
        xstable.make_xstable_model("fake", [1, 2, 3], [2, 3, 4],
                                   params=[p0], spectra=[model])



def test_check_intparam_is_not_empty():
    """I think we have to have at least one integrated parameter"""

    model = [1, 2, 3]
    with pytest.raises(ValueError, match="^params can not be empty$"):
        xstable.make_xstable_model("fake", [1, 2, 3], [2, 3, 4],
                                   params=[], spectra=[model])


@pytest.mark.parametrize("initial,delta,hardmin,softmin,softmax,hardmax,msg",
                         [(10, -1, 0, 0, 1, 1,
                           "initial=10 outside softmin=0 to softmax=1"),
                          (1, -1, 0.5, 0, 1, 1,
                           "softmin=0 < hardmin=0.5"),
                          (1, -1, 0, 0, 1, 0.5,
                           "softmax=1 > hardmax=0.5"),
                          (1, -1, 1, 1, 0, 0,
                           "softmin=1 > softmax=0"),
                         ])
def test_check_param_valid(initial, delta, hardmin, softmin, softmax, hardmax, msg):
    """Simple error checks"""

    with pytest.raises(ValueError, match=f"^Parameter nop {msg}$"):
        xstable.BaseParam("nop", initial, delta, hardmin, hardmax,
                          softmin=softmin, softmax=softmax)


# Note these are subtly different to test_check_param_valid since some
# of them have softmin/max set to None.
#
@pytest.mark.parametrize("initial,delta,hardmin,softmin,softmax,hardmax,values, msg",
                         [(10, -1, 0, None, None, 1, [0, 1],
                           "initial=10 outside softmin=0 to softmax=1"),
                          (1, -1, 0, None, None, 1, [-10, 1],
                           "hardmin=0 != first value=-10"),
                          (1, -1, 0, None, None, 1, [0, 10],
                           "hardmax=1 != last value=10"),
                          (1, -1, 0.5, 0, 1, 1, [0, 1],
                           "softmin=0 < hardmin=0.5"),
                          (1, -1, 0, 0, 1, 0.5, [0, 1],
                           "softmax=1 > hardmax=0.5"),
                          (1, -1, 1, 1, 0, 0, [0, 1],
                           "softmin=1 > softmax=0"),
                          (1, -1, 0, None, 1, 1, [],
                           "has no values"),
                          (1, -1, 0, 0, None, 1, [0, 1, 1],
                           "values are not monotonically increasing")
                         ])
def test_check_intparam_valid(initial, delta, hardmin, softmin, softmax, hardmax, values, msg):
    """Simple error checks"""

    with pytest.raises(ValueError, match=f"^Parameter nop {msg}$"):
        xstable.Param("nop", initial, delta, hardmin, hardmax,
                      softmin=softmin, softmax=softmax, values=values)


def mk(n):
    return xstable.Param(n, 0, 0.01, 0, 1, values=[0, 1])


def mka(n):
    return xstable.BaseParam(n, 0, 0.01, 0, 1)


@pytest.mark.parametrize("params,addparams",
                         [([mk("a"), mk("b"), mk("A")], None),
                          ([mk("a"), mk("b"), mk("x")], [mka("q"), mka("B")])
                         ])
def test_check_names_are_unique(params, addparams):
    """Ensure parameter names are unique"""

    elo = [1, 2, 3]
    ehi = [2, 3, 4]
    model = [1, 2, 3]
    with pytest.raises(ValueError, match="^Parameter names are not unique$"):
        xstable.make_xstable_model("fake", elo, ehi,
                                   params=params, addparams=addparams,
                                   spectra=[model] * 8)


def test_check_expected_spectra_size():
    """Expect 3 spectra but send in 2"""

    p0 = xstable.Param("nop", 1, -1, 0, 1, values=[0, 0.5, 1])
    model1 = [1, 2, 3]
    model2 = [2, 5, 6]
    with pytest.raises(ValueError, match="^Expected 3 spectra, found 2$"):
        xstable.make_xstable_model("fake", [1, 2, 3], [2, 3, 4],
                                   params=[p0], spectra=[model1, model2])


def test_check_additional_match():
    """If we send in addparams we need addspectra and vice versa

    This check is written with the knowledge the check is "symmetric",
    so that if we send in addparams and addspectra is empty there's
    also a failure.

    """

    p0 = xstable.Param("nop", 1, -1, 0, 1, values=[0, 1])
    model1 = [1, 2, 3]
    model2 = [2, 5, 6]
    with pytest.raises(ValueError, match="Mismatch between addparams and addspectra sizes: 0 1"):
        xstable.make_xstable_model("fake", [1, 2, 3], [2, 3, 4],
                                   params=[p0], spectra=[model1, model2],
                                   addspectra=[[model1, model2]])


def test_check_additional_match_number_of_spectra():
    """Check we send in the correct number of additional spectra."""

    p0 = xstable.Param("nop", 1, -1, 0, 1, values=[0, 1])
    ap0 = xstable.BaseParam("bob", 1, -1, 0, 1)
    model1 = [1, 2, 3]
    model2 = [2, 5, 6]
    with pytest.raises(ValueError, match="^Expected 2 spectra for additional parameter bob, found 1$"):
        xstable.make_xstable_model("fake", [1, 2, 3], [2, 3, 4],
                                   params=[p0], spectra=[model1, model2],
                                   addparams=[ap0],
                                   addspectra=[[model1]])


def test_check_additional_match_spectra_size():
    """Check the additional spectra nelem."""

    p0 = xstable.Param("nop", 1, -1, 0, 1, values=[0, 1])
    ap0 = xstable.BaseParam("bob", 1, -1, 0, 1)
    model1 = [1, 2, 3]
    model2 = [2, 5, 6]
    with pytest.raises(ValueError, match="^Spectrum for parameter bob should have 3 elements but has 2$"):
        xstable.make_xstable_model("fake", [1, 2, 3], [2, 3, 4],
                                   params=[p0], spectra=[model1, model2],
                                   addparams=[ap0],
                                   addspectra=[[model1, [2, 3]]])


def test_check_nxfp_match():
    """When NXFP is set we need more spectra"""

    p0 = xstable.Param("nop", 1, -1, 0, 1, values=[0, 1])
    ap0 = xstable.BaseParam("bob", 1, -1, 0, 1)
    model1 = [1, 2, 3]
    model2 = [2, 5, 6]
    with pytest.raises(ValueError, match="^Expected 4 spectra, found 2$"):
        xstable.make_xstable_model("fake", [1, 2, 3], [2, 3, 4],
                                   params=[p0],
                                   spectra=[model1, model2],
                                   xfxp=["foo: 1", "foo: 2"])


def make_single_model(elo, ehi, model, addmodel=True,
                      redshift=False, escale=False,
                      lolim=0, hilim=0):
    """Create a single-spectrum model"""

    p0 = xstable.Param("nop", 12, -1, 12, 12, values=[12])
    return xstable.make_xstable_model("fake", elo, ehi,
                                      params=[p0], spectra=[model],
                                      addmodel=addmodel,
                                      redshift=redshift,
                                      escale=escale,
                                      lolim=lolim, hilim=hilim)


def load_tmod(lbl, filename):
    """Read in an XSPEC table model.

    The test must use the requires_xspec fixture.
    """

    from sherpa.astro.xspec import XSTableModel, read_xstable_model

    mdl = read_xstable_model("foo", filename)
    assert isinstance(mdl, XSTableModel)
    return mdl


@requires_xspec
@requires_fits
def test_atable_single_spectrum(tmp_path):
    """Write out a model with a single spectrum.

    There is no parameter interpolation here!
    """

    outpath = tmp_path / "xs.mod"
    outfile = str(outpath)

    elo = [0.5, 0.7, 0.9]
    ehi = [0.7, 0.9, 1.0]
    model = [2, 4, 5]

    hdus = make_single_model(elo, ehi, model)
    xstable.write_xstable_model(outfile, hdus)

    mdl = load_tmod("foo", outfile)
    assert mdl.addmodel

    assert hasattr(mdl, "nop")
    assert not hasattr(mdl, "escale")
    assert not hasattr(mdl, "redshift")
    assert hasattr(mdl, "norm")
    assert len(mdl.pars) == 2

    assert mdl.nop.val == pytest.approx(12)
    assert mdl.nop.min == pytest.approx(12)
    assert mdl.nop.max == pytest.approx(12)
    assert mdl.nop.frozen

    # Evaluate the model with no interpolation
    #
    assert mdl(elo, ehi) == pytest.approx(model)

    # no apply normalization (the only parameter)
    #
    mdl.norm = 10
    assert mdl(elo, ehi) == pytest.approx([20, 40, 50])

    # Now check regridding and what happens outside the data range.
    #
    e2lo = [0.3, 0.4, 0.6, 0.8, 1.1]
    e2hi = [0.4, 0.6, 0.8, 1.1, 1.2]

    # This is what I expected, but XSPEC 12.13.1a does not return this
    # (I have asked HEASARC for confirmation but let's just treat this
    # as a regression test for now).
    #
    # expected = [0, 10, 10 + 20, 20 + 50, 0]
    expected = [0, 0, 10 + 20, 20 + 50, 0]
    assert mdl(e2lo, e2hi) == pytest.approx(expected)


@requires_xspec
@requires_fits
def test_mtable_single_spectrum(tmp_path):
    """Write out a model with a single spectrum.

    There is no parameter interpolation here!
    """

    outpath = tmp_path / "xs.mod"
    outfile = str(outpath)

    elo = [0.5, 0.7, 0.9]
    ehi = [0.7, 0.9, 1.0]
    model = [2, 4, 5]

    hdus = make_single_model(elo, ehi, model, addmodel=False)
    xstable.write_xstable_model(outfile, hdus)

    mdl = load_tmod("foo", outfile)
    assert not mdl.addmodel

    assert hasattr(mdl, "nop")
    assert not hasattr(mdl, "escale")
    assert not hasattr(mdl, "redshift")
    assert not hasattr(mdl, "norm")
    assert len(mdl.pars) == 1

    assert mdl.nop.val == pytest.approx(12)
    assert mdl.nop.min == pytest.approx(12)
    assert mdl.nop.max == pytest.approx(12)
    assert mdl.nop.frozen

    # Evaluate the model with no interpolation
    #
    assert mdl(elo, ehi) == pytest.approx(model)

    # Now check regridding and what happens outside the data range.
    #
    e2lo = [0.3, 0.4, 0.6, 0.8, 1.1]
    e2hi = [0.4, 0.6, 0.8, 1.1, 1.2]

    # This is a regression test based on the behavior of XSPEC 12.13.1a.
    #
    expected = [0, 0, 3, 14 / 3, 0]
    assert mdl(e2lo, e2hi) == pytest.approx(expected)


@requires_xspec
@requires_fits
def test_atable_interpolated(tmp_path):
    """Write out a model with multiple parameters.

    Two parameters:
         a   0, 5, 10          linear
         b   0.01, 0.1, 1, 10  logarithmic

    This is more a regression test than a check that the interpolation
    happens correctly (since we rely on XSPEC to do all of this for us
    anyway).

    """

    outpath = tmp_path / "xs.mod"
    outfile = str(outpath)

    egrid = np.arange(0.5, 1.6, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]
    emid = (elo + ehi) / 2

    p1 = xstable.Param("a", 1, 0.01, 0, 10, values=[0, 5, 10])
    p2 = xstable.Param("b", 0.2, -0.01, 0.01, 10,
                       values=[0.01, 0.1, 1, 10], loginterp=True)

    def mkspec(pb, pa):
        return (20 - pa) + 10**pb * emid

    pvals = [(0.01, 0), (0.1, 0), (1, 0), (10, 0),
             (0.01, 5), (0.1, 5), (1, 5), (10, 5),
             (0.01, 10), (0.1, 10), (1, 10), (10, 10)]
    spectra = [mkspec(*pv) for pv in pvals]

    hdus = xstable.make_xstable_model("fake", elo, ehi,
                                      params=[p1, p2], spectra=spectra)
    xstable.write_xstable_model(outfile, hdus)

    mdl = load_tmod("foo", outfile)
    assert mdl.addmodel

    assert hasattr(mdl, "a")
    assert hasattr(mdl, "b")
    assert not hasattr(mdl, "escale")
    assert not hasattr(mdl, "redshift")
    assert hasattr(mdl, "norm")
    assert len(mdl.pars) == 3

    assert mdl.a.val == pytest.approx(1)
    assert mdl.a.min == pytest.approx(0)
    assert mdl.a.max == pytest.approx(10)
    assert not mdl.a.frozen

    assert mdl.b.val == pytest.approx(0.2)
    assert mdl.b.min == pytest.approx(0.01)
    assert mdl.b.max == pytest.approx(10)
    assert mdl.b.frozen

    # Evaluate the model with no interpolation (ie 5,0.1 is one of the
    # cases we sent in).
    #
    mdl.a = 5
    mdl.b = 0.1
    expected = mkspec(0.1, 5)
    assert mdl(elo, ehi) == pytest.approx(expected)

    # Evaluate the model with interpolation. We pick a parameter
    # setting close enough we can use mkspec to generate the expected
    # result, but the required tolerance in the check is now large (as
    # the parameter values are not well chosen in this case). The
    # tolerance may need to be increased for other platforms (this was
    # generated with XSPEC 12.13.1a and Linux/x86).
    #
    mdl.a = 2
    mdl.b = 0.09
    expected = mkspec(0.09, 2)
    assert mdl(elo, ehi) == pytest.approx(expected, rel=2e-3)


@requires_xspec
@requires_fits
def test_atable_interpolated_with_additional(tmp_path):
    """Write out a model with multiple parameters and "additional" ones.

    Two parameters:
         a   0, 5, 10          linear
         b   0.01, 0.1, 1, 10  logarithmic

    then two additional parameters
         pa
         pb

    This is more a regressio test than a check that the interpolation
    happens correctly (since we rely on XSPEC to do all of this for
    us anyway).

    """

    outpath = tmp_path / "xs.mod"
    outfile = str(outpath)

    egrid = np.arange(0.5, 1.6, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]
    emid = (elo + ehi) / 2

    p1 = xstable.Param("a", 1, 0.01, 0, 10, values=[0, 5, 10])
    p2 = xstable.Param("b", 0.2, -0.01, 0.01, 10,
                       values=[0.01, 0.1, 1, 10], loginterp=True)

    ap1 = xstable.BaseParam("pa", 1, -0.01, 0, 10)
    ap2 = xstable.BaseParam("pb", 1, -0.01, 0, 10)

    def mkspec(pb, pa):
        return (20 - pa) + 10**pb * emid

    # I do not understand how the additional spectra work with the
    # parameter ranges, so just do something.
    #
    def mkaspec1(v):
        return v * emid

    def mkaspec2(v):
        return v * np.ones_like(emid)

    pvals = [(0.01, 0), (0.1, 0), (1, 0), (10, 0),
             (0.01, 5), (0.1, 5), (1, 5), (10, 5),
             (0.01, 10), (0.1, 10), (1, 10), (10, 10)]
    spectra = [mkspec(*pv) for pv in pvals]
    addspectra = [[mkaspec1(pv[0]) for pv in pvals],
                  [mkaspec2(pv[1]) for pv in pvals]]

    hdus = xstable.make_xstable_model("fake", elo, ehi,
                                      params=[p1, p2], spectra=spectra,
                                      addparams=[ap1, ap2],
                                      addspectra=addspectra)
    xstable.write_xstable_model(outfile, hdus)

    mdl = load_tmod("foo", outfile)
    assert mdl.addmodel

    assert hasattr(mdl, "a")
    assert hasattr(mdl, "b")
    assert hasattr(mdl, "pa")
    assert hasattr(mdl, "pb")
    assert not hasattr(mdl, "escale")
    assert not hasattr(mdl, "redshift")
    assert hasattr(mdl, "norm")
    assert len(mdl.pars) == 5

    assert mdl.a.val == pytest.approx(1)
    assert mdl.a.min == pytest.approx(0)
    assert mdl.a.max == pytest.approx(10)
    assert not mdl.a.frozen

    assert mdl.b.val == pytest.approx(0.2)
    assert mdl.b.min == pytest.approx(0.01)
    assert mdl.b.max == pytest.approx(10)
    assert mdl.b.frozen

    assert mdl.pa.val == pytest.approx(1)
    assert mdl.pa.min == pytest.approx(0)
    assert mdl.pa.max == pytest.approx(10)
    assert mdl.pa.frozen

    assert mdl.pb.val == pytest.approx(1)
    assert mdl.pb.min == pytest.approx(0)
    assert mdl.pb.max == pytest.approx(10)
    assert mdl.pb.frozen

    # A regression test.
    #
    mdl.a = 2
    mdl.b = 0.09
    mdl.pa = 2
    mdl.pb = 2
    expected = [22.791948, 22.93594, 23.079931, 23.223923,
                23.367912, 23.511902, 23.655893, 23.799883,
                23.943874, 24.087864]
    assert mdl(elo, ehi) == pytest.approx(expected)


# Used by test_addmodel_redshift/escale.
#
Z_POS = 69
Z_AMPL = 1.9814928


@requires_xspec
@requires_fits
def test_addmodel_resdhift(tmp_path, xsmodel):
    """Can we add a redshift parameter to a table model?

    Very-limited testing.
    """

    outpath = tmp_path / "xs.mod"
    outfile = str(outpath)

    egrid = np.arange(0.1, 2, 0.01)
    elo = egrid[:-1]
    ehi = egrid[1:]

    fake = xsmodel("gaussian")
    fake.linee = 1.6
    fake.norm = 50
    model = fake(elo, ehi)

    # check peak is where we expect
    assert np.argmax(model) == 150
    assert model[150] == pytest.approx(1.991392)

    hdus = make_single_model(elo, ehi, model, redshift=True)
    xstable.write_xstable_model(outfile, hdus)

    mdl = load_tmod("foo", outfile)
    assert mdl.addmodel

    assert hasattr(mdl, "nop")
    assert not hasattr(mdl, "escale")
    assert hasattr(mdl, "redshift")
    assert hasattr(mdl, "norm")
    assert len(mdl.pars) == 3

    assert mdl.nop.val == pytest.approx(12)
    assert mdl.nop.min == pytest.approx(12)
    assert mdl.nop.max == pytest.approx(12)
    assert mdl.nop.frozen

    # Evaluate the model with no interpolation
    #
    assert mdl(elo, ehi) == pytest.approx(model)

    # Check the redshift
    #
    mdl.redshift = 1

    got = mdl(elo, ehi)

    # The peak should have changed
    #
    assert np.argmax(got) == Z_POS
    assert got[Z_POS] == pytest.approx(Z_AMPL)


@requires_xspec
@requires_fits
def test_addmodel_escale(tmp_path, xsmodel):
    """Can we add a escale parameter to a table model?

    Very-limited testing.
    """

    outpath = tmp_path / "xs.mod"
    outfile = str(outpath)

    egrid = np.arange(0.1, 2, 0.01)
    elo = egrid[:-1]
    ehi = egrid[1:]

    fake = xsmodel("gaussian")
    fake.linee = 1.6
    fake.norm = 50
    model = fake(elo, ehi)

    # check peak is where we expect
    assert np.argmax(model) == 150
    assert model[150] == pytest.approx(1.991392)

    hdus = make_single_model(elo, ehi, model, escale=True)
    xstable.write_xstable_model(outfile, hdus)

    mdl = load_tmod("foo", outfile)
    assert mdl.addmodel

    assert hasattr(mdl, "nop")
    assert hasattr(mdl, "escale")
    assert not hasattr(mdl, "redshift")
    assert hasattr(mdl, "norm")
    assert len(mdl.pars) == 3

    assert mdl.nop.val == pytest.approx(12)
    assert mdl.nop.min == pytest.approx(12)
    assert mdl.nop.max == pytest.approx(12)
    assert mdl.nop.frozen

    # Evaluate the model with no interpolation
    #
    assert mdl(elo, ehi) == pytest.approx(model)

    # Check the shift (it is not clear to me what escale actually does
    # so this should be taken to be a regression test).
    #
    mdl.escale = 0.5

    got = mdl(elo, ehi)

    # The peak should have changed. We get the same shift as z=1
    # but the normalization changes (so we tie this test to the
    # redshift test).
    #
    assert np.argmax(got) == Z_POS
    assert got[Z_POS] == pytest.approx(2 * Z_AMPL)


@requires_xspec
@requires_fits
def test_atable_with_xfxp(tmp_path):
    """Write out a model with XFXP keywords.

    Note that Sherpa can not use these models. It is not entirely
    clear what the keys are meant to be.

    """

    outpath = tmp_path / "xs.mod"
    outfile = str(outpath)

    elo = [0.5, 0.7, 0.9]
    ehi = [0.7, 0.9, 1.0]
    model1 = np.asarray([2, 4, 5])
    model2 = np.asarray([5, 10, 40])

    p0 = xstable.Param("nop", 12, -1, 12, 20, values=[12, 20])
    hdus = xstable.make_xstable_model("fake", elo, ehi,
                                      params=[p0],
                                      spectra=[model1, model1 / 10,
                                               model2, model2 / 10],
                                      xfxp=["key: a", "key: b"])
    xstable.write_xstable_model(outfile, hdus)

    # We need a significant amount of work to support these models, so
    # we error out for now.
    #
    with pytest.raises(IOErr, match="^No support for NXFLTEXP=2 in "):
        load_tmod("foo", outfile)
