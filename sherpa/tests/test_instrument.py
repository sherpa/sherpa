#
#  Copyright (C) 2020, 2021, 2022, 2023
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
from sherpa.instrument import Kernel, PSFModel, RadialProfileKernel, \
    ConvolutionKernel, PSFKernel
from sherpa.models.basic import Box1D, Box2D, Const1D, Const2D, Gauss1D, \
    StepLo1D, TableModel
from sherpa.utils.err import PSFErr


def check_lines(out, lines):
    """Check out, when split on new line, matches lines.

    Parameters
    ----------
    out : str
        Text to check
    lines : list of str
        The expected output of out.split("\n")
    """

    olines = out.split("\n")
    for oline, expected in zip(olines, lines):
        assert oline == expected

    assert len(olines) == len(lines)


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

    check_lines(str(m),
                [ "pmodel1"
                , "   Param        Type          Value          Min          Max      Units"
                , "   -----        ----          -----          ---          ---      -----"
                , "   pmodel1.kernel frozen         box1"
                , "   pmodel1.size frozen           10           10           10"
                , "   pmodel1.center frozen            5            5            5"
                , "   pmodel1.radial frozen            0            0            1           "
                , "   pmodel1.norm frozen            1            0            1           "
                ])


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

    check_lines(str(m),
                [ "pdata2"
                , "   Param        Type          Value          Min          Max      Units"
                , "   -----        ----          -----          ---          ---      -----"
                , "   pdata2.kernel frozen        data2"
                , "   pdata2.size  frozen       (6, 6)       (6, 6)       (6, 6)"
                , "   pdata2.center frozen       (3, 3)       (3, 3)       (3, 3)"
                , "   pdata2.radial frozen            0            0            1           "
                , "   pdata2.norm  frozen            1            0            1           "
                ])


def test_psf2d_model_show():
    """What happens when the kernel is a model?"""

    m = PSFModel("pmodel2", make_2d_model())

    yy, xx = np.mgrid[1:10, 1:9]
    xx = xx.flatten()
    yy = yy.flatten()
    zz = np.arange(xx.size)
    dfold = Data2D('fold', xx, yy, zz)
    m.fold(dfold)

    check_lines(str(m),
                [ "pmodel2"
                , "   Param        Type          Value          Min          Max      Units"
                , "   -----        ----          -----          ---          ---      -----"
                , "   pmodel2.kernel frozen         box2"
                , "   pmodel2.size frozen     (72, 72)     (72, 72)     (72, 72)"
                , "   pmodel2.center frozen     (36, 36)     (36, 36)     (36, 36)"
                , "   pmodel2.radial frozen            0            0            1           "
                , "   pmodel2.norm frozen            1            0            1           "
                ])


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


def test_psfmodel_kernel_has_no_dimension():
    """It is expected that this will error out, but just where.

    It was intended to catch a different error condition but
    it is hard to trigger without writing very low-level test
    code, so just check what happens here.
    """

    x = np.arange(0.5, 10.5, 1)
    y = np.ones_like(x)
    data = Data1D("data-data", x, y)

    m = PSFModel(kernel=TableModel())
    with pytest.raises(PSFErr,
                       match="PSF model dimension must be <= 2"):
        m.get_kernel(data)


def test_psf1d_no_kernel():
    """Error out if there's no kernel"""

    m = PSFModel('bob')
    b = Box1D()
    with pytest.raises(PSFErr,
                       match="PSF kernel has not been set"):
        m(b)


def test_psf1d_fold_no_kernel():
    """Error out if there's no kernel"""

    m = PSFModel('bob')
    dfold = Data1D('fold', np.arange(10), np.zeros(10))
    with pytest.raises(PSFErr,
                       match="model 'bob' does not have an associated PSF function"):
        m.fold(dfold)


def test_psf1d_no_fold():
    """Error out if there's no kernel"""

    box = Box1D()
    psf = PSFModel('bob', box)

    cpt = Gauss1D()
    sm = psf(cpt)

    with pytest.raises(PSFErr,
                       match="PSF model has not been folded"):
        sm([1, 2, 3, 4, 5])


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
    """Check that we error out"""

    psf = PSFModel('psf', make_1d_model())
    d = make_2d_data()

    with pytest.raises(PSFErr,
                       match="kernel 'box1' and data 'twod' do not match: 1D vs 2D"):
        psf.fold(d)


def test_psf1d_data_given_2d_dataset():
    """Check that we error out"""

    psf = PSFModel('psf', make_1d_data())
    d = make_2d_data()

    with pytest.raises(PSFErr,
                       match="kernel 'oned' and data 'twod' do not match: 1D vs 2D"):
        psf.fold(d)


def test_psf2d_model_given_1d_dataset():
    """Check that we error out"""

    psf = PSFModel('psf', make_2d_model())
    d = make_1d_data()

    with pytest.raises(PSFErr,
                       match="kernel 'box2' and data 'oned' do not match: 2D vs 1D"):
        psf.fold(d)


def test_psf2d_data_given_1d_dataset():
    """Check that we error out"""

    psf = PSFModel('psf', make_2d_data())
    d = make_1d_data()

    with pytest.raises(PSFErr,
                       match="kernel 'twod' and data 'oned' do not match: 2D vs 1D"):
        psf.fold(d)


@pytest.mark.parametrize("kernel_func", [make_1d_model, make_1d_data])
@pytest.mark.parametrize("field", ["size", "origin", "center"])
@pytest.mark.parametrize("vals,ndim", [([2, 3], 2), ([2, 3, 4], 3)])
def test_psf1d_set_invalid_field(kernel_func, field, vals, ndim):
    """Check that we error out"""

    kernel = kernel_func()
    m = PSFModel("p1", kernel)

    with pytest.raises(PSFErr,
                       match=f"kernel 'p1' and data '{field}' do not match: 1D vs {ndim}D"):
        setattr(m, field, vals)


@pytest.mark.parametrize("kernel_func", [make_2d_model, make_2d_data])
@pytest.mark.parametrize("field", ["size", "origin", "center"])
@pytest.mark.parametrize("vals,ndim", [(4, 1), ([2, 3, 4], 3)])
def test_psf2d_set_invalid_field(kernel_func, field, vals, ndim):
    """Check that we error out"""

    kernel = kernel_func()
    m = PSFModel("p2", kernel)

    with pytest.raises(PSFErr,
                       match=f"kernel 'p2' and data '{field}' do not match: 2D vs {ndim}D"):
        setattr(m, field, vals)


@pytest.mark.parametrize("kernel_func", [make_none, make_1d_model, make_1d_data, make_2d_model, make_2d_data])
def test_psf_set_model_to_model(kernel_func):
    """Can we change the model field? model

    This is a regression test.
    """

    m = PSFModel(kernel=kernel_func())
    with pytest.raises(AttributeError,
                       match="'PSFModel' object attribute 'model' " +
                       "cannot be replaced with a callable attribute"):
        m.model = Box1D()


@pytest.mark.parametrize("kernel_func", [make_none, make_1d_model, make_1d_data, make_2d_model, make_2d_data])
def test_psf_set_model_to_bool(kernel_func):
    """Can we change the model field? boolean

    This is a regression test.
    """

    m = PSFModel(kernel=kernel_func())
    # This used to check the message that was raised, but it changes
    # with Python versions (3.9, 3.10, 3.11 all had different messages)
    # and so we now just check the error is raised.
    with pytest.raises(AttributeError):
        m.model = False


@pytest.mark.parametrize("dshape,kshape,exc,msg",
                         [(None, None, TypeError,
                           "dshape must be a sequence"),
                          (None, [1], TypeError,
                           "dshape must be a sequence"),
                          ([1], None, TypeError,
                           "kshape must be a sequence"),
                          ([], [], ValueError,
                           "0D kernel is not supported"),
                          ([1], [1, 2], ValueError,
                           "dshape and kshape must be the same size, not 1 and 2"),
                          ([1, 2], [1], ValueError,
                           "dshape and kshape must be the same size, not 2 and 1"),
                          ([], [1], ValueError,
                           "dshape and kshape must be the same size, not 0 and 1"),
                          ([1], [], ValueError,
                           "dshape and kshape must be the same size, not 1 and 0"),
                          ([1, 3, 4], ([1, 2, 3]), PSFErr,
                           "PSF model dimension must be <= 2")])
def test_kernel_checks_arguments(dshape, kshape, exc, msg):
    """Check that we error out"""

    with pytest.raises(exc, match=msg):
        Kernel(dshape, kshape)


def test_radialprofile_is_1d():
    """Check that we error out"""

    with pytest.raises(PSFErr,
                       match="^Radial profile requires 1D data, not 2D$"):
        RadialProfileKernel([10, 10], [4, 3])


def test_kernel_repr():

    assert repr(Kernel([5], [3])) == "<Kernel kernel instance>"


def test_convolutionkernel_repr():

    assert repr(ConvolutionKernel(None)) == "<ConvolutionKernel kernel instance>"


def test_kernel_show_1d():

    kern = Kernel([5], [3])
    check_lines(str(kern),
                [ "dshape   = [5]"
                , "kshape   = [3]"
                , "skshape  = None"
                , "norm     = False"
                , "origin   = [0.]"
                , "frozen   = True"
                , "center   = None"
                , "args     = []"
                , "kwargs   = {}"
                , "renorm_shape  = None"
                , "renorm   = None"
                , "do_pad   = False"
                , "pad_mask = None"
                , "frac     = None"
                ])


def test_convolutionkernel_show_2d():

    mdl = Const2D()
    mdl.c0.set(2, min=2, max=2)
    kern = ConvolutionKernel(mdl)
    check_lines(str(kern),
                [ "Convolution Kernel:"
                , "const2d"
                , "   Param        Type          Value          Min          Max      Units"
                , "   -----        ----          -----          ---          ---      -----"
                , "   const2d.c0   thawed            2            2            2           "
                ])


def test_psfkernel_show_1d():

    kern = PSFKernel([5], [3])
    check_lines(str(kern),
                [ "dshape   = [5]"
                , "kshape   = [3]"
                , "skshape  = None"
                , "norm     = True"
                , "origin   = None"
                , "frozen   = True"
                , "center   = None"
                , "args     = []"
                , "kwargs   = {}"
                , "renorm_shape  = None"
                , "renorm   = None"
                , "do_pad   = False"
                , "pad_mask = None"
                , "frac     = None"
                , "is_model = False"
                , "size     = None"
                , "lo       = None"
                , "hi       = None"
                , "width    = None"
                , "radial   = 0"
                ])


def test_radialprofilekernel_show_1d():

    kern = RadialProfileKernel([5], [3])
    check_lines(str(kern),
                [ "dshape   = [5]"
                , "kshape   = [3]"
                , "skshape  = None"
                , "norm     = True"
                , "origin   = None"
                , "frozen   = True"
                , "center   = None"
                , "args     = []"
                , "kwargs   = {}"
                , "renorm_shape  = None"
                , "renorm   = None"
                , "do_pad   = False"
                , "pad_mask = None"
                , "frac     = None"
                , "is_model = False"
                , "size     = None"
                , "lo       = None"
                , "hi       = None"
                , "width    = None"
                , "radial   = 1"
                , "radialsize = None"
                ])


def test_empty_psfmodel_has_no_kernel():

    assert PSFModel().kernel is None


def test_psfmodel_set_kernel():

    mdl = PSFModel(kernel=Box1D())
    assert isinstance(mdl.kernel, Box1D)


def test_psfmodel_can_not_set_scalar_kernel():
    """This is really a check of a non-callable kernel

    This is a regression test.
    """

    with pytest.raises(PSFErr,
                       match="model 'psfmodel' does not have an associated PSF function"):
        PSFModel().kernel = 23


def test_psfmodel_clearing_model_kernel_to_none():
    """The code suggests this should be possible, but ither
    parts of the system means it can't happen, at least like
    this.
    """

    mdl = PSFModel(kernel=Box1D())

    with pytest.raises(AttributeError,
                       match="'PSFModel' object attribute 'kernel' cannot be replaced with a non-callable attribute"):
        mdl.kernel = None


def test_psfmodel_clearing_data_kernel_to_none():
    """We can clear a data kernel, but what else happens?

    This is a regression test.
    """

    x = np.asarray([1, 2, 3])
    y = np.asarray([0.7, 0.2, 0.1])
    mdl = PSFModel(kernel=Data1D("foo", x, y))

    assert mdl.ndim == 1

    # DJB does not understand what these fields really mean,
    # but let's just check if they are changed or not.
    #
    mdl.center = 3
    mdl.origin = 2

    mdl.kernel = None
    assert mdl.kernel is None
    assert mdl.ndim is None
    assert mdl.center is None
    assert mdl.origin is None


def test_psfmodel_data_kernel_keep_dim():
    """What changes when the data dimensionality is the same?"""

    x = np.asarray([1, 2, 3])
    y = np.asarray([0.7, 0.2, 0.1])
    mdl = PSFModel(kernel=Data1D("foo", x, y))
    assert mdl.ndim == 1
    assert len(mdl.kernel) == 3

    # DJB does not understand what these fields really mean,
    # but let's just check if they are changed or not.
    #
    mdl.center = 3
    mdl.origin = 2

    x2 = np.asarray([300, 500, 1000, 2000])
    y2 = np.asarray([10, 14, 12, 2])
    mdl.kernel = Data1D("bar", x2, y2)
    assert isinstance(mdl.kernel, Data1D)
    assert mdl.ndim == 1
    assert len(mdl.kernel) == 4
    assert mdl.center == pytest.approx(3)
    assert mdl.origin == pytest.approx(2)


def test_psfmodel_data_kernel_change_dim():
    """What changes when the data dimensionality changes?"""

    x = np.asarray([1, 2, 3])
    y = np.asarray([0.7, 0.2, 0.1])
    mdl = PSFModel(kernel=Data1D("foo", x, y))
    assert len(mdl.kernel) == 3

    # DJB does not understand what these fields really mean,
    # but let's just check if they are changed or not.
    #
    mdl.center = 3
    mdl.origin = 2

    x0 = np.asarray([300, 500, 1000, 2000, 5000])
    x1 = np.asarray([10, 14, 12, 2, 7])
    mdl.kernel = Data2D("bar", x0, x1, np.ones(5))
    assert isinstance(mdl.kernel, Data2D)
    assert mdl.ndim == 2
    assert len(mdl.kernel) == 5
    assert mdl.center is None
    assert mdl.origin is None


def test_psfmodel_get_kernel_subkernel():
    """Check code has been tested

    It's not obvious what is going on here, so just check this
    particular code path has been run.
    """

    kernel = Data1D("input",
                    np.asarray([1, 2, 3]),
                    np.asarray([7, 2, 1]))

    mdl = PSFModel(kernel=kernel)

    x = np.asarray([1, 2, 3, 4, 5])
    y = np.asarray([10, 14, 12, 2, 7])
    response = mdl.get_kernel(Data1D("other", x, y))

    assert isinstance(response, Data1D)
    assert response.name == "kernel"
    assert response.x == pytest.approx([1, 2, 3])
    assert response.y == pytest.approx([0.7, 0.2, 0.1])


def test_psfmodel_get_kernel_no_subkernel():
    """Check code has been tested

    It's not obvious what is going on here, so just check this
    particular code path has been run.
    """

    kernel = Data1D("input",
                    np.asarray([1, 2, 3]),
                    np.asarray([7, 2, 1]))

    mdl = PSFModel(kernel=kernel)

    x = np.asarray([1, 2, 3, 4, 5])
    y = np.asarray([10, 14, 12, 2, 7])
    response = mdl.get_kernel(Data1D("other", x, y), subkernel=False)

    assert isinstance(response, Data1D)
    assert response.name == "kernel"
    assert response.x == pytest.approx([1, 2, 3])
    assert response.y == pytest.approx([7, 2, 1])
