#
#  Copyright (C) 2021
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

"""Basic tests of sherpa.utils._psf"""

import numpy as np

import pytest

from sherpa.models.basic import Gauss1D, StepLo1D
from sherpa.utils import _psf


def test_exports():
    expected = ['extract_kernel', 'get_padsize', 'normalize',
                'pad_bounding_box', 'pad_data', 'set_origin',
                'tcdData', 'unpad_data']

    found = [n for n in dir(_psf) if not n.startswith('_')]
    assert found == expected


def test_create_tcdData():
    """Check we get an object with the expected methods"""

    expected = ['clear_kernel_fft', 'convolve']

    out = _psf.tcdData()
    found = [n for n in dir(out) if not n.startswith('_')]

    assert found == expected


def test_tcdData_convolve_error1():

    out = _psf.tcdData()
    data = [1, 2, 3]
    kernel = [1, 2, 2, 1]
    with pytest.raises(TypeError) as te:
        out.convolve(data, kernel, [3], [2, 2], [0])

    assert str(te.value) == 'input array sizes do not match, dims_src: 1 vs dims_kern: 2'


@pytest.mark.parametrize('center', [[0], 0])
def test_tcdData_convolve_error2(center):

    out = _psf.tcdData()
    data = [1, 2, 3, 4]
    kernel = [1, 2, 2, 1]
    with pytest.raises(TypeError) as te:
        out.convolve(data, kernel, [2, 2], [2, 2], center)

    assert str(te.value) == 'input array sizes do not match, dims_kern: 2 vs center: 1'


def test_tcdData_convolve_error3():

    out = _psf.tcdData()
    data = [1, 2, 3]
    kernel = [1, 2, 2, 1]
    with pytest.raises(TypeError) as te:
        out.convolve(data, kernel, [2, 2], [2, 2], [0, 0])

    assert str(te.value) == 'input array size do not match dimensions, source size: 3 vs source dim: 4'


def test_tcdData_convolve_error4():

    out = _psf.tcdData()
    data = [1, 2, 3, 2]
    kernel = [1, 2, 2]
    with pytest.raises(TypeError) as te:
        out.convolve(data, kernel, [2, 2], [2, 2], [0, 0])

    assert str(te.value) == 'input array size do not match dimensions, kernel size: 3 vs kernel dim: 4'


def test_tcdData_convolve_error_3d():

    out = _psf.tcdData()
    data = [1, 2, 3, 4, 5, 6, 7, 8]
    kernel = [1, 2, 2, 1]
    with pytest.raises(TypeError) as te:
        out.convolve(data, kernel, [2, 2, 2], [2, 2, 1], [0, 0, 0])

    assert str(te.value) == 'Padding dimension not supported'


@pytest.mark.parametrize('actual,expected',
                         [(-1, 2), (0, 2), (1, 2), (2, 2), (3, 3), (6, 6),
                          (7, 8), (32000, 32000), (32399, 32400),
                          (32400, 32400)])
def test_get_padsize(actual, expected):
    """These are taken from the hard-coded table, so will need updating
    wenever we update it in the code.
    """

    got = _psf.get_padsize(actual)
    assert got == expected


def test_get_padsize_too_large():
    """Pick one larger than the current maximum"""

    with pytest.raises(TypeError) as te:
        _psf.get_padsize(32401)

    assert str(te.value) == 'Padding dimension length 32401 not supported'


def test_pad_data_basic():
    """Not 100% sure what it's meant to do but check something."""

    out = _psf.pad_data([1, 2, 3, 4, 5, 6, 7], [7], [10])
    assert len(out) == 10
    assert isinstance(out, np.ndarray)
    assert out == pytest.approx([1, 2, 3, 4, 5, 6, 7, 0, 0, 0])


def test_unpad_data_basic():
    """Not 100% sure what it's meant to do but check something.

    Check we can "invert" pad_data
    """

    out = _psf.unpad_data([1, 2, 3, 4, 5, 6, 7, 0, 0, 0], [10], [7])
    assert len(out) == 7
    assert isinstance(out, np.ndarray)
    assert out == pytest.approx([1, 2, 3, 4, 5, 6, 7])


@pytest.mark.parametrize('func,i,j',
                         [(_psf.pad_data, 1, 2),
                          (_psf.unpad_data, 2, 1)])
def test_xpad_data_error1(func, i, j):

    with pytest.raises(TypeError) as te:
        func([1, 2, 3], [3], [1, 3])

    assert str(te.value) == f'input array sizes do not match, shape: {i} vs padshape: {j}'


@pytest.mark.parametrize('func,i',
                         [(_psf.pad_data, 0),
                          (_psf.unpad_data, 1)])
def test_xpad_data_error2(func, i):

    with pytest.raises(TypeError) as te:
        func([1, 2, 3], [3, 1], [1, 3])

    assert str(te.value) == f'pad size is smaller than data shape, padshape[{i}]: 1 < shape[{i}]: 3'


def test_pad_data_error3():

    with pytest.raises(TypeError) as te:
        _psf.pad_data([1, 2, 3], [2], [3])

    assert str(te.value) == 'input array size do not match dimensions, kernel size: 3 vs kernel dim: 2'


def test_unpad_data_error3():

    with pytest.raises(TypeError) as te:
        _psf.unpad_data([1, 2, 3], [4], [3])

    assert str(te.value) == 'input array size do not match dimensions, kernel size: 3 vs kernel dim: 4'


@pytest.mark.parametrize('func', [_psf.pad_data, _psf.unpad_data])
def test_xpad_data_error4(func):

    with pytest.raises(TypeError) as te:
        func([1, 2, 3, 4, 5, 6, 7, 8], [2, 2, 2], [2, 2, 2])

    if func.__name__ == 'pad_data':
        emsg = 'padding kernel failed - dimension unsupported'
    else:
        emsg = 'unpadding kernel failed-dimension unsupported'

    assert str(te.value) == emsg


@pytest.mark.parametrize('flags,expected',
                         [([True] * 4, [1, 2, 3, 4]),
                          ([True] * 10, [1, 2, 3, 4, 0, 0, 0, 0, 0, 0]),
                          ([True, True, False, False], [1, 2, 0, 0]),
                          ([True, False, False, True], [1, 0, 0, 2]),
                          ([True, False, False, True, False, False, False],
                           [1, 0, 0, 2, 0, 0, 0]),
                          ([True, False, False, True, False, True, False, True, False, False],
                           [1, 0, 0, 2, 0, 3, 0, 4, 0, 0])
                          ])
def test_pad_bounding_box_basic(flags, expected):
    """Not 100% sure what it's meant to do."""

    out = _psf.pad_bounding_box([1, 2, 3, 4], flags)
    assert out == pytest.approx(expected)


def test_pad_bounding_box_error():

    with pytest.raises(TypeError) as te:
        _psf.pad_bounding_box([1, 2, 3, 4, 5],
                              [4, 3, 0, 1])

    assert str(te.value) == 'kernel size: 5 is > than mask size: 4'


def test_convolve_simple_1d():
    """Very-low-level test"""

    data = np.asarray([0, 0, 0, 1, 2, 4, 3, 0, 0, 1, 0])
    kernel = np.asarray([2, 4, 1])

    tcd = _psf.tcdData()
    out = tcd.convolve(data, kernel, data.shape, kernel.shape, [0])

    # As the kernel is not normalized, we know the output is larger than
    # the input: data.sum() == 11, kernel.sum() == 7. This isn't
    # really needed given the check against np.convolve, but keep
    # in as a check.
    #
    assert out.sum() == pytest.approx(11 * 7)

    # the convolution doesn't quite match numpy.convolve
    expected = np.convolve(data, kernel, mode='same')
    expected = np.roll(expected, 1)
    assert out == pytest.approx(expected)

    # and now with the kernel internally stored
    out2 = tcd.convolve(data, kernel, data.shape, kernel.shape, [0])
    assert out2 == pytest.approx(out)


@pytest.mark.xfail  # see #1334
def test_convolve_combined_1d():
    """Try to replicate the logic of test_psf1d_step_v2
    from sherpa/tests/test_instrument.py
    """

    smdl = StepLo1D()
    smdl.xcut = 100
    smdl.ampl = 10

    gsmooth = Gauss1D()

    x = np.arange(0, 200, 0.5)

    data = smdl(x)
    kernel = gsmooth(x)
    kernel /= kernel.sum()

    tcd = _psf.tcdData()
    y = tcd.convolve(data, kernel, data.shape, kernel.shape, [0])

    # So the output is not easy to describe analytically, hence
    # we just check parts of it.
    #
    assert y[(x >= 19.5) & (x <= 100)] == pytest.approx([10] * 162, abs=1e-4)
    assert y[x >= 119] == pytest.approx([0] * 162, abs=1e-4)

    # check that the x <= 19 values are in ascending order
    y1 = y[x <= 19]
    assert (y1[1:] > y1[:-1]).all()

    # and now with the kernel internally stored
    y2 = tcd.convolve(data, kernel, data.shape, kernel.shape, [0])
    assert y2 == pytest.approx(y)
