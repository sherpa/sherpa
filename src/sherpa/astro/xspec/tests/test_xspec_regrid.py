#
#  Copyright (C) 2020, 2023, 2024
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

"""Check we can regrid XSPEC models.

Basic testing is based on the XSPEC wabs and powerlaw models, since
they are "simple" and unlikely to change significantly. They are
multiplicative and additive models respectively.

These tests do not exercise the whole set of regrid options, since
it is expected that these should just work (from other tests). This
is to check sepcialized behavior with the XSPEC models.

"""

import numpy as np
import pytest

from sherpa.models.basic import Const1D, Gauss1D
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_xspec


# To make sure we track current behavior, we check that routines fail
# in the expected way rather than just marking them as xfail.  This
# way, when support is finally unlocked, we should find out from the
# tests (rather that just going from xfail to xpass which is easy to
# miss). It also lets us check why a test fails (in case the failure
# mode changes due to other parts of Sherpa).
#
IntegrateError = "'integrate' is an invalid keyword argument for this function"


@requires_xspec
@pytest.mark.parametrize("mname", ["wabs", "powerlaw"])
def test_regrid_name(mname, xsmodel):
    """What is the name of the regrid component"""

    ebase = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]

    mdl = xsmodel(mname, "base")
    regrid = mdl.regrid(ebase[:-1], ebase[1:])

    assert mdl.name == 'base'
    assert regrid.name == 'regrid1d(base)'


@requires_xspec
def test_regrid_name_combined():
    """What is the name of the regrid component"""

    from sherpa.astro import xspec

    ebase = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]

    c1 = xspec.XSwabs()
    c2 = xspec.XSpowerlaw()
    mdl = c1 * c2
    regrid = mdl.regrid(ebase[:-1], ebase[1:])

    assert mdl.name == 'wabs * powerlaw'
    assert regrid.name == 'regrid1d(wabs * powerlaw)'


@requires_xspec
@pytest.mark.parametrize("mname", ["wabs", "powerlaw"])
def test_regrid_does_not_require_bins(mname, xsmodel):
    """The regrid method does not require lo,hi bins but running it does"""

    ebase = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]

    mdl = xsmodel(mname, "base")
    regrid = mdl.regrid(ebase)

    emsg = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg):
        with pytest.raises(TypeError) as exc:
            regrid([0.4, 0.5, 0.6])

    # assert str(exc.value) == 'calc() requires pars,lo,hi arguments, sent 2 arguments'
    # assert str(exc.value).endswith('() takes no keyword arguments')
    assert str(exc.value) == IntegrateError


@requires_data
@requires_fits
@requires_xspec
def test_regrid_table_requires_bins(make_data_path):
    """Does the regrid method require lo,hi bins for XSPEC tables?"""

    from sherpa.astro import xspec

    path = make_data_path('testpcfabs.mod')
    tbl = xspec.read_xstable_model('bar', path)

    ebase = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]

    regrid = tbl.regrid(ebase)

    emsg = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg):
        y = regrid([0.52, 0.57, 0.62])

    # Answer calculated with XSPEC 12.12.1
    assert y == pytest.approx([0.50037217, 0.50017911, 0.50086778])


@requires_xspec
@pytest.mark.parametrize("mname", ["wabs", "powerlaw"])
def test_regrid_identity(mname, xsmodel):
    """Check regrid returns the same data when grids are equal"""

    elo = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85]
    ehi = [0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]

    mdl = xsmodel(mname, "base")
    regrid = mdl.regrid(elo, ehi)

    # scale y values to make them closer to unity
    y1 = 100 * mdl(elo, ehi)
    with pytest.raises(TypeError, match=IntegrateError):
        y2 = 100 * regrid(elo, ehi)

    # assert y2 == pytest.approx(y1)


@requires_xspec
def test_regrid_identity_combined():
    """Check regrid returns the same data when grids are equal"""

    from sherpa.astro import xspec

    elo = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85]
    ehi = [0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]

    c1 = xspec.XSwabs()
    c2 = xspec.XSpowerlaw()
    c1.nh = 0.05
    c2.norm = 100

    # The mdl.regrid call fails with a TypeError:
    # 'integrate' is an invalid keyword argument for this function
    #
    mdl = c1 * c2
    regrid = mdl.regrid(elo, ehi)

    # scale y values to make them closer to unity
    y1 = mdl(elo, ehi)
    with pytest.raises(TypeError, match=IntegrateError):
        y2 = regrid(elo, ehi)

    # assert y2 == pytest.approx(y1)


@requires_xspec
def test_additive():
    """Simple test of an additive model"""

    from sherpa.astro import xspec

    ebase = np.arange(1.1, 1.5, 0.01)
    egrid = np.arange(1.1, 1.5, 0.05)

    mdl = xspec.XSpowerlaw('add')
    regrid = mdl.regrid(ebase[:-1], ebase[1:])

    # will change regrid too
    mdl.norm = 100

    y1 = mdl(egrid[:-1], egrid[1:])
    with pytest.raises(TypeError, match=IntegrateError):
        y2 = regrid(egrid[:-1], egrid[1:])

    # assert y2 == pytest.approx(y1)


# ideally we would test all three cases
#   no overlap
#   overlaps low
#   overlaps high
#
overlap_none = np.arange(10, 15, 0.1)
overlap_low = np.arange(0.8, 1.3, 0.05)
overlap_high = np.arange(1.2, 1.8, 0.05)

ans_none = np.zeros(overlap_none.size - 1)
ans_low = np.asarray([0, 0, 0, 0, 0, 0,
                      4.44517626, 4.25596144, 4.08219945])
ans_high = np.asarray([4.08219945, 3.92207132, 3.7740328, 3.63676442,
                       3.50913198, 2.72125635,  # this is a partial bin
                       0, 0, 0, 0, 0, 0])


@requires_xspec
@pytest.mark.parametrize("egrid,yexp",
                         [(overlap_none, ans_none),
                          (overlap_low, ans_low),
                          (overlap_high, ans_high)
                         ])
def test_additive_overlap(egrid, yexp):
    """Simple test of an additive model.

    There is only partial overlap of the two grids.
    """

    from sherpa.astro import xspec

    ebase = np.arange(1.1, 1.5, 0.01)

    mdl = xspec.XSpowerlaw()
    regrid = mdl.regrid(ebase[:-1], ebase[1:])

    # will change regrid too
    mdl.norm = 100
    with pytest.raises(TypeError, match=IntegrateError):
        y2 = regrid(egrid[:-1], egrid[1:])

    # assert y2 == pytest.approx(yexp)


# To make it easy to compare the multiplicative tests we want
# the center of the bins to match. To do this I am just shifting
# the bins, so for an output of
#
#   0.1-0.2   entered at 0.15
#   0.2-0.3              0.25
#   0.3-0.4              0.35 ...
#
# we use a "regrid" that looks like
#   0.125-0.175  this is centered at 0.15
#   0.175-0.225
#   0.225-0.275  this is centered at 0.25
#   0.275-0.325
#

@requires_xspec
def test_multiplicative():
    """Simple test of a multiplicative model"""

    from sherpa.astro import xspec

    ebase = np.arange(0.475, 1.025, 0.05)
    egrid = np.arange(0.5, 1.0, 0.1)

    mdl = xspec.XSwabs('mul')
    regrid = mdl.regrid(ebase[:-1], ebase[1:])

    mdl.nh = 0.05

    y1 = mdl(egrid[:-1], egrid[1:])
    with pytest.raises(TypeError, match=IntegrateError):
        y2 = regrid(egrid[:-1], egrid[1:])

    # assert y2 == pytest.approx(y1, rel=0.04)


# The expected values have been calculated based on akima interpolation
# of the grided values, and differ from evaluating mdl on egrid (for
# those bins which overlap the ebase range).
#
overlap2_low = np.arange(0.4, 0.9, 0.1)
overlap2_high = np.arange(0.7, 1.3, 0.1)

ans2_low = np.asarray([0, 0.641675, 0.7399567, 0.7984249])
ans2_high = np.asarray([0.7984249, 0.84757835, 0.8698429, 0, 0, 0])


@requires_xspec
@pytest.mark.parametrize("egrid,yexp",
                         [(overlap_none, ans_none),
                          (overlap2_low, ans2_low),
                          (overlap2_high, ans2_high)
                         ])
def test_multiplicative_overlap(egrid, yexp):
    """Simple test of a multiplicative model

    There is only partial overlap of the two grids.
    """

    from sherpa.astro import xspec

    ebase = np.arange(0.475, 1.025, 0.05)

    mdl = xspec.XSwabs()
    regrid = mdl.regrid(ebase[:-1], ebase[1:])

    # will change regrid too
    mdl.nh = 0.05

    with pytest.raises(TypeError, match=IntegrateError):
        y2 = regrid(egrid[:-1], egrid[1:])

    # assert y2 == pytest.approx(yexp)


@requires_xspec
@pytest.mark.parametrize("name1,par1,val1,name2,par2,val2",
                         [('wabs', 'nh', 0.05,
                           'powerlaw', 'norm', 100),
                          ('powerlaw', 'norm', 100,
                           'wabs', 'nh', 0.05)])
def test_combined(name1, par1, val1, name2, par2, val2, xsmodel):
    """Simple test of an additive model * multiplicative model

    We want to check it works for both a * b and
    b * a (there is logic in the BinaryOp model
    that favors the left-hand model, so we want to check
    it is handled correctly here).
    """

    ebase = np.arange(0.475, 1.025, 0.05)
    egrid = np.arange(0.5, 1.0, 0.1)

    com1 = xsmodel(name1, "com1")
    com2 = xsmodel(name2, "com2")
    setattr(com1, par1, val1)
    setattr(com2, par2, val2)

    mdl = com1 * com2
    regrid = mdl.regrid(ebase[:-1], ebase[1:])

    y1 = mdl(egrid[:-1], egrid[1:])
    with pytest.raises(TypeError, match=IntegrateError):
        y2 = regrid(egrid[:-1], egrid[1:])

    # assert y2 == pytest.approx(y1, rel=0.04)


@requires_xspec
@pytest.mark.parametrize("name,par,val",
                         [('wabs', 'nh', 0.05),
                          ('powerlaw', 'norm', 100)])
def test_combined_arithmetic_left(name, par, val, xsmodel):
    """Simple test of an additive model * multiplicative model"""

    ebase = np.arange(0.475, 1.025, 0.05)
    egrid = np.arange(0.5, 1.0, 0.1)

    com = xsmodel(name, "com")
    setattr(com, par, val)

    mdl = 4 * com
    regrid = mdl.regrid(ebase[:-1], ebase[1:])

    yexp = 4 * com(egrid[:-1], egrid[1:])
    with pytest.raises(TypeError, match=IntegrateError):
        y1 = regrid(egrid[:-1], egrid[1:])

    # assert y1 == pytest.approx(yexp, rel=0.04)


@requires_xspec
@pytest.mark.parametrize("name,par,val",
                         [('wabs', 'nh', 0.05),
                          ('powerlaw', 'norm', 100)])
def test_combined_arithmetic_right(name, par, val, xsmodel):
    """Simple test of an additive model * multiplicative model"""

    ebase = np.arange(0.475, 1.025, 0.05)
    egrid = np.arange(0.5, 1.0, 0.1)

    com = xsmodel(name, "com")
    setattr(com, par, val)

    mdl = com * 4
    regrid = mdl.regrid(ebase[:-1], ebase[1:])

    yexp = 4 * com(egrid[:-1], egrid[1:])
    with pytest.raises(TypeError, match=IntegrateError):
        y1 = regrid(egrid[:-1], egrid[1:])

    # assert y1 == pytest.approx(yexp, rel=0.04)


@requires_xspec
def test_multi_combined_additive():
    """Can we handle a "deep" binop tree?

    What happens with (2 * mul) / (3 + mul)? Does it get
    treated as additive?
    """

    from sherpa.astro import xspec

    ebase = np.arange(0.475, 1.025, 0.05)
    egrid = np.arange(0.5, 1.0, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]

    m1 = xspec.XSgaussian('m1')
    m2 = xspec.XSgaussian('m2')

    mdl = (2 * m1) / (3 + m2)
    with pytest.raises(RecursionError):
        regrid = mdl.regrid(ebase[:-1], ebase[1:])

    # expected = mdl(elo, ehi)
    # with pytest.raises(TypeError, match=IntegrateError):
    #     got = regrid(elo, ehi)

    ## assert got == pytest.approx(expected)


@requires_xspec
def test_multi_combined_multiplicative():
    """Can we handle a "deep" binop tree?

    What happens with (2 * mul) / (3 + mul)? Does it get
    treated as multiplicative?
    """

    from sherpa.astro import xspec

    ebase = np.arange(0.475, 1.025, 0.05)
    egrid = np.arange(0.5, 1.0, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]

    m1 = xspec.XSconstant('m1')
    m2 = xspec.XSconstant('m2')

    mdl = (2 * m1) / (3 + m2)
    with pytest.raises(RecursionError):
        regrid = mdl.regrid(ebase[:-1], ebase[1:])

    # expected = mdl(elo, ehi)
    # with pytest.raises(TypeError, match=IntegrateError):
    #     got = regrid(elo, ehi)

    # # assert got == pytest.approx(expected)


@requires_xspec
@pytest.mark.parametrize('sherpa_first', [True, False])
def test_sherpa_mul_xspec_add(sherpa_first):
    """Check sherpa (multiplicative) * xspec (additive)"""

    from sherpa.astro import xspec

    ebase = np.arange(0.475, 1.025, 0.05)
    egrid = np.arange(0.5, 1.0, 0.1)
    eg1 = egrid[:-1]
    eg2 = egrid[1:]

    msherpa = Const1D('m1')
    msherpa.c0 = 2
    msherpa.integrate = False

    mxspec = xspec.XSpowerlaw('m2')
    mxspec.phoindex = 1.8

    if sherpa_first:
        comb = msherpa * mxspec
    else:
        comb = mxspec * msherpa

    rcomb = comb.regrid(ebase[:-1], ebase[1:])

    # Check the results make sense, but do not check the absolute values
    yexp = 2 * mxspec(eg1, eg2)

    # Note the tolerance is different to test_sherpa_add_xspec_mul
    assert comb(eg1, eg2) == pytest.approx(yexp)
    with pytest.raises(TypeError, match=IntegrateError):
        assert rcomb(eg1, eg2) == pytest.approx(yexp, rel=0.006)


@requires_xspec
@pytest.mark.parametrize('sherpa_first', [True, False])
def test_sherpa_add_xspec_mul(sherpa_first):
    """Check sherpa (additive) * xspec (multiplicative)"""

    from sherpa.astro import xspec

    ebase = np.arange(0.475, 1.025, 0.05)
    egrid = np.arange(0.5, 1.0, 0.1)
    eg1 = egrid[:-1]
    eg2 = egrid[1:]

    msherpa = Gauss1D('m1')
    msherpa.pos = 1.02
    msherpa.fwhm = 0.4
    msherpa.integrate = True

    mxspec = xspec.XSconstant('m2')
    mxspec.factor = 2

    if sherpa_first:
        comb = msherpa * mxspec
    else:
        comb = mxspec * msherpa

    rcomb = comb.regrid(ebase[:-1], ebase[1:])

    # Check the results make sense, but do not check the absolute values
    yexp = 2 * msherpa(eg1, eg2)

    # Note the tolerance is different to test_sherpa_mul_xspec_add
    assert comb(eg1, eg2) == pytest.approx(yexp)
    with pytest.raises(TypeError, match=IntegrateError):
        assert rcomb(eg1, eg2) == pytest.approx(yexp, rel=0.07)


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize('name,iflag,tol',
                         [('xspec-tablemodel-RCS.mod', True, 1e-6),
                          pytest.param('testpcfabs.mod', False, 1e-3, marks=pytest.mark.xfail)])  # XFAIL: iflag is wrong and values do not match
def test_regrid_table(name, iflag, tol, make_data_path):
    """Can we regrid a table model?

    We test out both additive and multiplicative models.
    """

    from sherpa.astro import xspec as xs

    infile = make_data_path(name)
    tbl = xs.read_xstable_model('mod', infile)

    assert tbl.ndim == 1
    assert tbl.integrate == iflag

    egrid = np.arange(0.5, 1.5, 0.1)
    eg1 = egrid[:-1]
    eg2 = egrid[1:]

    exp = tbl(eg1, eg2)
    assert (exp > 0).all()

    ebase = np.arange(0.45, 1.55, 0.05)
    rtbl = tbl.regrid(ebase[:-1], ebase[1:])
    y = rtbl(eg1, eg2)

    assert y == pytest.approx(exp, rel=tol)


@requires_xspec
@pytest.mark.parametrize("mname", ["wabs", "powerlaw"])
def test_regrid_convolution_requires_bins(mname, xsmodel):
    """Does the regrid method require lo,hi bins?

    This just checks the current behavior, which is about what is
    allowed and what warnings are created, not the returned values.

    """

    ebase = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]

    conv = xsmodel("cflux")
    mdl = xsmodel(mname, "base")
    cmdl = conv(mdl)

    rmdl = cmdl.regrid(ebase)
    with pytest.warns(FutureWarning) as record:
        rmdl([0.57, 0.62, 0.66, 0.72])

    assert len(record) == 3
    assert record[0].message.args[0] == 'calc() requires pars,lo,hi arguments, sent 2 arguments'
    assert record[1].message.args[0] == 'calc() requires pars,rhs,lo,hi arguments, sent 3 arguments'
    assert record[2].message.args[0] == 'calc() requires pars,lo,hi arguments, sent 2 arguments'


@requires_xspec
@pytest.mark.parametrize('model,tol',
                         [('powerlaw', 1e-3),
                          ('wabs', 0.72)])  # was tol=0.31 at one point
def test_regrid_convolution_single(model, tol, xsmodel):
    """regriding an additive model or multiplicative.
    """

    ebase = np.arange(0.4, 10, 0.1)
    egrid = np.arange(0.5, 10, 0.2)
    eg1 = egrid[:-1]
    eg2 = egrid[1:]

    conv = xsmodel("cflux")
    conv.emin = 1
    conv.emax = 9

    mdl = xsmodel(model)
    cmdl = conv(mdl)

    rmdl = cmdl.regrid(ebase[:-1], ebase[1:])

    exp = cmdl(eg1, eg2)
    assert (exp > 0).all()

    y = rmdl(eg1, eg2)
    assert y == pytest.approx(exp, rel=tol)


@requires_xspec
@pytest.mark.parametrize('name1,name2',
                         [('wabs', 'powerlaw'),
                          ('powerlaw', 'wabs')])
def test_regrid_convolution_combined(name1, name2, xsmodel):
    """regriding convolved(mdl1 * mdl2)

    The tolerance is very large
    """

    ebase = np.arange(0.4, 10, 0.1)
    egrid = np.arange(0.5, 10, 0.2)
    eg1 = egrid[:-1]
    eg2 = egrid[1:]

    conv = xsmodel("cflux")
    conv.emin = 1
    conv.emax = 9

    mdl1 = xsmodel(name1)
    mdl2 = xsmodel(name2)
    cmdl = conv(mdl1 * mdl2)

    rmdl = cmdl.regrid(ebase[:-1], ebase[1:])

    exp = cmdl(eg1, eg2)
    assert (exp > 0).all()

    y = rmdl(eg1, eg2)
    assert y == pytest.approx(exp, rel=0.35)
