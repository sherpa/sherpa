#
#  Copyright (C) 2020, 2021, 2022
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

"""
Add tests for sherpa.instrument that are not covered elsewhere.
"""

import logging

import numpy as np

import pytest

from sherpa.data import Data1D, Data2D
from sherpa.instrument import Kernel, PSFModel, RadialProfileKernel
from sherpa.models.basic import Box1D, Box2D, Const2D, Const1D, Gauss1D, Gauss2D, StepLo1D
from sherpa.utils.err import PSFErr


def make_none():
    """We can use this to create no kernel"""

    return None


def make_1d_model():
    """We want to create a new model each call.

    This could be made a fixture but this is easier.
    """

    box = Box1D('box1')
    box.xlow = 2
    box.xhi = 5
    return box


def make_2d_model():
    """We want to create a new model each call."""

    box = Box2D('box2')
    box.xlow = 2
    box.xhi = 5
    box.yhi = 4
    return box


def make_1d_data():
    """We want to create a new object each call."""

    k = np.asarray([2, 3, 5, 1, 0, 2, 1])
    x = np.arange(k.size)
    return Data1D('oned', x, k)


def make_2d_data():
    """We want to create a new object each call."""

    y, x = np.mgrid[1:3, 1:4]
    x = x.flatten()
    y = y.flatten()
    z = np.ones(x.size)
    return Data2D('twod', x, y, z)


def test_psf1d_show():
    """Test the __str__ method

    Loop through basic data and then add in the options
    (but not all possible combinations).
    """

    def check(x, n, ans):
        toks = x[n].split()
        assert toks == ans

    NAME = ['psfmodel']
    PARAMS = ['Param', 'Type', 'Value', 'Min', 'Max', 'Units']
    LINES = ['-----', '----', '-----', '---', '---', '-----']
    KERNEL = ['psfmodel.kernel', 'frozen', 'oned']
    SIZE = ['psfmodel.size', 'frozen', '7', '7', '7']
    CENTER = ['psfmodel.center', 'frozen', '3', '3', '3']
    RADIAL = ['psfmodel.radial', 'frozen', '0', '0', '1']
    NORM = ['psfmodel.norm', 'frozen', '1', '0', '1']

    m = PSFModel()

    # basic settings
    out = str(m).split('\n')

    assert len(out) == 5
    check(out, 0, NAME)
    check(out, 1, PARAMS)
    check(out, 2, LINES)
    check(out, 3, RADIAL)
    check(out, 4, NORM)

    m.kernel = make_1d_data()

    # before a fold you don't get the size and center parameters
    out = str(m).split('\n')

    assert len(out) == 6
    check(out, 0, NAME)
    check(out, 1, PARAMS)
    check(out, 2, LINES)
    check(out, 3, KERNEL)
    check(out, 4, RADIAL)
    check(out, 5, NORM)

    dfold = Data1D('fold', np.arange(10), np.zeros(10))
    m.fold(dfold)

    out = str(m).split('\n')

    assert len(out) == 8
    check(out, 0, NAME)
    check(out, 1, PARAMS)
    check(out, 2, LINES)
    check(out, 3, KERNEL)
    check(out, 4, SIZE)
    check(out, 5, CENTER)
    check(out, 6, RADIAL)
    check(out, 7, NORM)


def test_psf1d_model_show():
    """What happens when the kernel is a model?"""

    box1 = make_1d_model()
    m = PSFModel("pmodel1", box1)

    dfold = Data1D('fold', np.arange(10), np.zeros(10))
    m.fold(dfold)

    out = str(m).split("\n")
    assert len(out) == 8
    assert out[0] == "pmodel1"
    assert out[1] == "   Param        Type          Value          Min          Max      Units"
    assert out[2] == "   -----        ----          -----          ---          ---      -----"
    assert out[3] == "   pmodel1.kernel frozen         box1"
    assert out[4] == "   pmodel1.size frozen           10           10           10"
    assert out[5] == "   pmodel1.center frozen            5            5            5"
    assert out[6] == "   pmodel1.radial frozen            0            0            1           "
    assert out[7] == "   pmodel1.norm frozen            1            0            1           "


def test_psf2d_data_show():
    """What happens when the kernel is a Data2D instance?"""

    y, x = np.mgrid[1:4, 1:3]
    x = x.flatten()
    y = y.flatten()
    z = np.ones(x.size)
    data2 = Data2D('data2', x, y, z)

    m = PSFModel("pdata2", data2)

    yy, xx = np.mgrid[1:10, 1:9]
    xx = xx.flatten()
    yy = yy.flatten()
    zz = np.arange(xx.size)
    dfold = Data2D('fold', xx, yy, zz)
    m.fold(dfold)

    out = str(m).split("\n")
    assert len(out) == 8
    assert out[0] == "pdata2"
    assert out[1] == "   Param        Type          Value          Min          Max      Units"
    assert out[2] == "   -----        ----          -----          ---          ---      -----"
    assert out[3] == "   pdata2.kernel frozen        data2"
    assert out[4] == "   pdata2.size  frozen       (6, 6)       (6, 6)       (6, 6)"
    assert out[5] == "   pdata2.center frozen       (3, 3)       (3, 3)       (3, 3)"
    assert out[6] == "   pdata2.radial frozen            0            0            1           "
    assert out[7] == "   pdata2.norm  frozen            1            0            1           "


def test_psf2d_model_show():
    """What happens when the kernel is a model?"""

    m = PSFModel("pmodel2", make_2d_model())

    yy, xx = np.mgrid[1:10, 1:9]
    xx = xx.flatten()
    yy = yy.flatten()
    zz = np.arange(xx.size)
    dfold = Data2D('fold', xx, yy, zz)
    m.fold(dfold)

    print(m)
    out = str(m).split("\n")
    assert len(out) == 8
    assert out[0] == "pmodel2"
    assert out[1] == "   Param        Type          Value          Min          Max      Units"
    assert out[2] == "   -----        ----          -----          ---          ---      -----"
    assert out[3] == "   pmodel2.kernel frozen         box2"
    assert out[4] == "   pmodel2.size frozen     (72, 72)     (72, 72)     (72, 72)"
    assert out[5] == "   pmodel2.center frozen     (36, 36)     (36, 36)     (36, 36)"
    assert out[6] == "   pmodel2.radial frozen            0            0            1           "
    assert out[7] == "   pmodel2.norm frozen            1            0            1           "


def test_psf1d_empty_pars():
    """What does .pars mean for an empty PSFModel?"""

    m = PSFModel()
    assert m.pars == ()


def test_psf1d_pars():
    """What does .pars mean for a PSFModel?"""

    b = Box1D()
    m = PSFModel(kernel=b)
    assert m.pars == ()


def test_psf1d_convolved_pars():
    """What does .pars mean for a PSFModel applied to a model?"""

    b1 = Box1D('b1')
    m = PSFModel(kernel=b1)
    b2 = Box1D('b2')
    c = m(b2)

    b1.xlow = 1
    b1.xhi = 10
    b1.ampl = 0.2

    b2.xlow = 4
    b2.xhi = 8
    b2.ampl = 0.4

    bpars = b1.pars + b2.pars
    assert len(bpars) == 6

    cpars = c.pars
    assert len(cpars) == 6
    for bpar, cpar in zip(bpars, cpars):
        assert cpar == bpar


def test_psf1d_no_kernel():
    """Error out if there's no kernel"""

    m = PSFModel('bob')
    b = Box1D()
    with pytest.raises(PSFErr) as exc:
        m(b)

    assert "PSF kernel has not been set" == str(exc.value)


def test_psf1d_fold_no_kernel():
    """Error out if there's no kernel"""

    m = PSFModel('bob')
    dfold = Data1D('fold', np.arange(10), np.zeros(10))
    with pytest.raises(PSFErr) as exc:
        m.fold(dfold)

    assert "model 'bob' does not have an associated PSF function" == str(exc.value)


def test_psf1d_no_fold():
    """Error out if there's no kernel"""

    box = Box1D()
    psf = PSFModel('bob', box)

    cpt = Gauss1D()
    sm = psf(cpt)

    with pytest.raises(PSFErr) as exc:
        sm([1, 2, 3, 4, 5])

    assert "PSF model has not been folded" == str(exc.value)


def test_psf1d_kernel_data(caplog):
    """Access the kernel data: subkernel=True"""

    k = np.asarray([2, 5, 3])
    x = 5 + np.arange(k.size)
    d = Data1D('my-kernel', x, k)

    m = PSFModel(kernel=d)

    dfold = Data1D('fold', np.arange(10), np.zeros(10))

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ans = m.get_kernel(dfold)

    assert len(caplog.records) == 1
    r = caplog.record_tuples[0]
    assert r[0] == 'sherpa.instrument'
    assert r[1] == logging.INFO
    assert r[2] == "PSF frac: 1.0"

    assert isinstance(ans, Data1D)
    assert ans.name == 'kernel'
    # integers, so treat as exact
    assert (ans.x == np.arange(5, 8)).all()
    assert ans.y == pytest.approx(k / k.sum())
    assert ans.staterror is None
    assert ans.syserror is None


def test_psf1d_kernel_model(caplog):
    """Access the kernel data: subkernel=True"""

    k = Box1D()
    k.xlow = 1
    k.xhi = 3
    m = PSFModel(kernel=k)

    dfold = Data1D('fold', np.arange(10), np.zeros(10))

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ans = m.get_kernel(dfold)

    assert len(caplog.records) == 1
    r = caplog.record_tuples[0]
    assert r[0] == 'sherpa.instrument'
    assert r[1] == logging.INFO
    assert r[2] == "PSF frac: 1.0"

    assert isinstance(ans, Data1D)
    assert ans.name == 'kernel'

    # integers, so treat as exact
    assert (ans.x == np.arange(10)).all()

    # box1D between 1 and 3 inclusive
    y = np.asarray([0, 1, 1, 1, 0, 0, 0, 0, 0, 0])
    assert ans.y == pytest.approx(y / y.sum())

    assert ans.staterror is None
    assert ans.syserror is None


def test_psf2d_kernel_data(caplog):
    """Access the kernel data: no subkernel"""

    x1, x0 = np.mgrid[-5:-3, 10:13]
    s = x0.shape

    k = np.ones(s)
    k[1, 1] = 2

    x0 = x1.flatten()
    x1 = x1.flatten()
    k = k.flatten()
    d = Data2D('my-kernel', x0, x1, k, shape=s)

    m = PSFModel('x', kernel=d)

    d1, d0 = np.mgrid[5:15, -5:5]
    ds = d0.shape
    z = np.zeros(ds)
    dfold = Data2D('fold', d0.flatten(), d1.flatten(), z.flatten(), shape=ds)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ans = m.get_kernel(dfold)

    assert len(caplog.records) == 1
    r = caplog.record_tuples[0]
    assert r[0] == 'sherpa.instrument'
    assert r[1] == logging.INFO
    assert r[2] == "PSF frac: 1.0"

    assert isinstance(ans, Data2D)
    assert ans.name == 'kernel'
    assert (ans.shape == [2, 3]).all()

    # integers, so treat as exact
    assert (ans.x0 == x0).all()
    assert (ans.x1 == x1).all()

    assert ans.y == pytest.approx(k / k.sum())

    assert ans.staterror is None
    assert ans.syserror is None


@pytest.mark.xfail  # see #1428
def test_psf2d_kernel_model(caplog):
    """Access the kernel data: no subkernel"""

    k = Box2D()
    k.xlow = -1
    k.xhi = 1
    k.ylow = 10
    k.yhi = 12
    m = PSFModel('x', kernel=k)

    d1, d0 = np.mgrid[5:15, -5:4]
    ds = d0.shape

    d0 = d0.flatten()
    d1 = d1.flatten()
    z = np.zeros(ds).flatten()

    dfold = Data2D('fold', d0, d1, z, shape=ds)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ans = m.get_kernel(dfold)

    assert len(caplog.records) == 1
    r = caplog.record_tuples[0]
    assert r[0] == 'sherpa.instrument'
    assert r[1] == logging.INFO
    assert r[2] == "PSF frac: 1.0"

    assert isinstance(ans, Data2D)
    assert ans.name == 'kernel'
    assert (ans.shape == [10, 9]).all()

    # integers, so treat as exact
    assert (ans.x0 == d0).all()
    assert (ans.x1 == d1).all()

    # This is the old check, converted to what I assume is meant to be
    # the correct form - see #1428 - but ans.y is all zeros except for
    # a single 1 whilst k is not (it has 9 True values and the rest
    # False). I do not have time to identify the truth here.
    #
    k = (ans.x0 >= -1) & (ans.x0 <= 1) & (ans.x1 >= 10) & (ans.x1 <= 12)
    assert ans.y == pytest.approx(k / k.sum())

    assert ans.staterror is None
    assert ans.syserror is None


def test_psf1d_flat():
    """This is based on
    sherpa.models.tests.test_regrid_unit.test_regrid1_works_with_convolution_style
    but I wanted to make sure we have an explicit check of the underlying
    code.
    """

    cmdl = Const1D()
    cmdl.c0 = -500

    gsmooth = Gauss1D()
    gsmooth.fwhm = 3
    psf = PSFModel('psf', gsmooth)

    x = np.arange(5, 23, 3)
    d = Data1D('fake', x, x * 0)
    psf.fold(d)

    smoothed = psf(cmdl)
    y = smoothed(x)

    assert y == pytest.approx([-500] * 6)


def test_psf1d_step():
    """This is based on
    sherpa.models.tests.test_regrid_unit.test_regrid1_works_with_convolution_style
    but I wanted to make sure we have an explicit check of the underlying
    code.
    """

    smdl = StepLo1D()
    smdl.xcut = 12.5
    smdl.ampl = 10

    gsmooth = Gauss1D()
    gsmooth.fwhm = 3
    psf = PSFModel('psf', gsmooth)

    x = np.arange(5, 23, 3)
    d = Data1D('fake', x, x * 0)
    psf.fold(d)

    smoothed = psf(smdl)
    y = smoothed(x)

    assert y == pytest.approx([10.0, 10.0, 10.0, 0, 0, 0], abs=1e-4)


@pytest.mark.xfail  # see #1334
def test_psf1d_step_v2():
    """Trying to track down why we have seen different behavior in
    test_regrid_unit.py.
    """

    smdl = StepLo1D()
    smdl.xcut = 100
    smdl.ampl = 10

    gsmooth = Gauss1D()
    psf = PSFModel('psf', gsmooth)

    x = np.arange(0, 200, 0.5)
    d = Data1D('fake', x, x * 0)
    psf.fold(d)

    smoothed = psf(smdl)
    y = smoothed(x)

    # So the output is not easy to describe analytically, hence
    # we just check parts of it.
    #
    assert y[(x >= 19.5) & (x <= 100)] == pytest.approx([10] * 162, abs=1e-4)
    assert y[x >= 119] == pytest.approx([0] * 162, abs=1e-4)

    # check that the x <= 19 values are in ascending order
    y1 = y[x <= 19]
    assert (y1[1:] > y1[:-1]).all()


def test_psf1d_combined():
    """This is based on
    sherpa.models.tests.test_regrid_unit.test_regrid1_works_with_convolution_style
    but I wanted to make sure we have an explicit check of the underlying
    code.
    """

    smdl = StepLo1D()
    smdl.xcut = 12.5
    smdl.ampl = 10

    cmdl = Const1D()
    cmdl.c0 = -500

    imdl = smdl + cmdl

    gsmooth = Gauss1D()
    gsmooth.fwhm = 3
    psf = PSFModel('psf', gsmooth)

    x = np.arange(5, 23, 3)
    d = Data1D('fake', x, x * 0)
    psf.fold(d)

    smoothed = psf(imdl)
    y = smoothed(x)

    assert y == pytest.approx([-490, -490, -490, -500, -500, -500], rel=7e-3)


def test_psf1d_combined_v2():
    """See test_psf1d_step_v2"""

    smdl = StepLo1D()
    smdl.xcut = 100
    smdl.ampl = 10

    cmdl = Const1D()
    cmdl.c0 = -500

    imdl = smdl + cmdl

    gsmooth = Gauss1D()
    psf = PSFModel('psf', gsmooth)

    x = np.arange(0, 200, 0.5)
    d = Data1D('fake', x, x * 0)
    psf.fold(d)

    smoothed = psf(imdl)
    y = smoothed(x)

    # So the output is not easy to describe analytically, hence
    # we just check parts of it.
    #
    assert y[(x >= 19.5) & (x <= 100)] == pytest.approx([-490] * 162, abs=1e-4)
    assert y[x >= 119] == pytest.approx([-500] * 162, abs=1e-4)

    # check that the x <= 19 values are in ascending order
    y1 = y[x <= 19]
    assert (y1[1:] > y1[:-1]).all()


def test_psf1d_model_given_2d_dataset():
    """Do we error out or not?

    This is a regression test.
    """

    psf = PSFModel('psf', make_1d_model())
    d = make_2d_data()

    # Does this error out?
    psf.fold(d)

    smdl = StepLo1D()
    smdl.xcut = 100
    smdl.ampl = 10

    cmdl = Const1D()
    cmdl.c0 = -500

    imdl = smdl + cmdl

    # Or maybe this errors out?
    smoothed = psf(imdl)
    with pytest.raises(TypeError) as err:
        smoothed(np.arange(1, 7))

    # It is not at all obvious why we have 36 for the source dim.
    #
    assert str(err.value) == "input array sizes do not match dimensions, source size: 6 vs source dim: 36"


def test_psf1d_data_given_2d_dataset():
    """Do we error out or not?

    This is a regression test.
    """

    psf = PSFModel('psf', make_1d_data())
    d = make_2d_data()

    # Does this error out?
    psf.fold(d)

    smdl = StepLo1D()
    smdl.xcut = 100
    smdl.ampl = 10

    cmdl = Const1D()
    cmdl.c0 = -500

    imdl = smdl + cmdl

    # Or maybe this errors out?
    smoothed = psf(imdl)
    with pytest.raises(TypeError) as err:
        smoothed(np.arange(1, 7))

    assert str(err.value) == "input array sizes do not match, dims_src: 2 vs dims_kern: 1"


def test_psf2d_model_given_1d_dataset():
    """Do we error out or not?

    This is a regression test.
    """

    psf = PSFModel('psf', make_2d_model())
    d = make_1d_data()

    # Does this error out?
    psf.fold(d)

    smdl = Gauss2D()
    smdl.ampl = 1000

    cmdl = Const2D()
    cmdl.c0 = -500

    imdl = smdl + cmdl

    # Or maybe this errors out?
    smoothed = psf(imdl)
    with pytest.raises(TypeError) as err:
        smoothed(np.arange(1, 7))

    assert str(err.value) == "function missing required argument 'x1lo' (pos 3)"


def test_psf2d_data_given_1d_dataset():
    """Do we error out or not?

    This is a regression test.
    """

    psf = PSFModel('psf', make_2d_data())
    d = make_1d_data()

    # Does this error out?
    psf.fold(d)

    smdl = Gauss2D()
    smdl.ampl = 1000

    cmdl = Const2D()
    cmdl.c0 = -500

    imdl = smdl + cmdl

    # Or maybe this errors out?
    smoothed = psf(imdl)
    with pytest.raises(TypeError) as err:
        smoothed(np.arange(1, 7))

    assert str(err.value) == "function missing required argument 'x1lo' (pos 3)"


@pytest.mark.parametrize("kernel_func", [make_1d_model, make_1d_data])
@pytest.mark.parametrize("field", ["size", "origin", "center"])
@pytest.mark.parametrize("vals", [[2, 3], [2, 3, 4]])
def test_psf1d_set_invalid_field(kernel_func, field, vals):
    """What happens if we send in a 2D array?

    This is a regression test (since it currently does not error out).
    """

    kernel = kernel_func()
    m = PSFModel("p1", kernel)

    setattr(m, field, vals)
    assert getattr(m, field) == pytest.approx(vals)


@pytest.mark.parametrize("kernel_func", [make_2d_model, make_2d_data])
@pytest.mark.parametrize("field", ["size", "origin", "center"])
@pytest.mark.parametrize("vals", [4, [2, 3, 4]])
def test_psf2d_set_invalid_field(kernel_func, field, vals):
    """What happens if we send in a 2D array?

    This is a regression test (since it currently does not error out).
    """

    kernel = kernel_func()
    m = PSFModel("p2", kernel)

    setattr(m, field, vals)
    assert getattr(m, field) == pytest.approx(vals)


@pytest.mark.parametrize("kernel_func", [make_none, make_1d_model, make_1d_data, make_2d_model, make_2d_data])
def test_psf_set_model_to_model(kernel_func):
    """Can we change the model field? model

    This is a regression test.
    """

    m = PSFModel(kernel=kernel_func())
    with pytest.raises(AttributeError) as err:
        m.model = Box1D()

    assert str(err.value) == "'PSFModel' object attribute 'model' cannot be replaced with a callable attribute"


@pytest.mark.parametrize("kernel_func", [make_none, make_1d_model, make_1d_data, make_2d_model, make_2d_data])
def test_psf_set_model_to_bool(kernel_func):
    """Can we change the model field? boolean

    This is a regression test.
    """

    m = PSFModel(kernel=kernel_func())
    # This does not error out
    m.model = False


@pytest.mark.parametrize("dshape,kshape",
                         [pytest.param(None, None, marks=pytest.mark.xfail), (None, [1]), pytest.param([1], None, marks=pytest.mark.xfail),
                          ([1], [1, 2]), ([1, 2], [1]),
                          ([], [1]), ([1], []),
                          ([1, 3, 4], ([1, 2, 3]))])
def test_kernel_checks_arguments(dshape, kshape):
    """Do these error out?

    This is a regression test.
    """

    # XFAIL: when kshape is None a TypeError is raised from len(None)
    Kernel(dshape, kshape)


def test_radialprofile_is_1d():
    """Do these error out?

    This is a regression test.
    """

    RadialProfileKernel([10, 10], [4, 3])
