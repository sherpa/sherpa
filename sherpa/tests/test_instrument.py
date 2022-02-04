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
from sherpa.instrument import PSFModel
from sherpa.models.basic import Box1D, Box2D, Const1D, Gauss1D, StepLo1D
from sherpa.utils.err import PSFErr


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
    KERNEL = ['psfmodel.kernel', 'frozen', 'my-kernel']
    SIZE = ['psfmodel.size', 'frozen', '3', '3', '3']
    CENTER = ['psfmodel.center', 'frozen', '1', '1', '1']
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

    k = np.asarray([2, 5, 3])
    x = np.arange(k.size)
    d = Data1D('my-kernel', x, k)
    m.kernel = d

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
