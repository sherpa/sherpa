#
#  Copyright (C) 2026
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

from types import SimpleNamespace

import numpy as np

import pytest


# Allow the tests to be skipped if WCS support is not present.
#
_wcs = pytest.importorskip("sherpa.astro.utils._wcs")


@pytest.mark.parametrize("method", ["pix2world", "world2pix"])
@pytest.mark.parametrize("x0,x1",
                         [ ([1, 2], [1]),
                           ([1], [1, 2]),
                           ([], [2, 3]),
                           ([2, 3], []),
                           (None, [1, 2]),
                           ([1, 2], None),
                          ])
def test_check_array_size(method, x0, x1):
    """Do we error out if x0/x1 are not correct?"""

    # This test is written assuming the size checks are before other
    # checks, so things like worrying about all the ctypes can be ignored
    # here.
    #
    func = getattr(_wcs, method)
    with pytest.raises(TypeError,
                       match="^input pixel arrays x0,x1 not of equal length$"):
        func("LINEAR", x0, x1, [1, 1], [10, 10], [2, 2], 0, 2000, 2000)


@pytest.mark.parametrize("method", ["pix2world", "world2pix"])
@pytest.mark.parametrize("ncrpix,ncrval,ncdelt",
                         [ (2, 2, 1),
                           (2, 1, 2),
                           (1, 2, 2),
                           pytest.param(1, 1, 1, marks=pytest.mark.xfail),  # does not error out here as the arrays have the same length
                           [2, 2, 0],
                           (2, 0, 2),
                           (0, 2, 2),
                          ])
def test_check_wcs_size_mismatch(method, ncrpix, ncrval, ncdelt):
    """Do we error out if crpix/val/delt size are not correct?"""

    # This test is written assuming the size checks are before other
    # checks, so things like worrying about all the ctypes can be ignored
    # here.
    #
    func = getattr(_wcs, method)

    crpix = [1, 2][:ncrpix]
    crval = [10, 20][:ncrval]
    cdelt = [2, 1][:ncdelt]
    with pytest.raises(TypeError,
                       match="^input array sizes do not match, cdelt: "):
        func("LINEAR", [1, 2], [2, 1], crpix, crval, cdelt, 0, 2000, 2000)


def setup_linear():

    return SimpleNamespace(
        ctype = "LINEAR",
        x0 = np.asarray([2.5, 10, 1.5, -3]),
        x1 = np.asarray([3.0, 4, 1.5, 2]),
        w0 = np.asarray([118.0, 103, 120, 129]),
        w1 = np.asarray([0.5, 3.5, -4, -2.5]),
        crpix = [1.5, 1.5],
        crval = [120, -4],
        cdelt = [-2, 3],
        crota = 0,
        equinox = 2000,
        epoch = 2000,
    )


def setup_wcs():

    return SimpleNamespace(
        ctype = "WCS",
        x0 = np.asarray([2.5, 10, 1.5, -3]),
        x1 = np.asarray([3.0, 4, 1.5, 2]),
        w0 = np.asarray([119.79940227, 118.29479111, 120., 120.90228831]),
        w1 = np.asarray([-4.44996358, -4.74786331, -4., -4.14948687]),
        crpix = [1.5, 1.5],
        crval = [120, -4],
        cdelt = [-0.2, -0.3],
        crota = 0,
        equinox = 2000,
        epoch = 2000,
    )


def setup_tanp():

    return SimpleNamespace(
        ctype = "TAN-P",
        x0 = np.asarray([128.5, 130, 120, 128.5]),
        x1 = np.asarray([128.5, 120, 40, 128.5]),
        w0 = np.asarray([270.0, 280.0079798, 264.51384357, 270.0]),
        w1 = np.asarray([0.0, 0.03774769, 0.3888148, 0.0]),
        # looks like we want crval to be 0 and cdelt to be positive/equal
        crpix = [128.5, 128.5],
        crval = [0, 0],
        cdelt = [0.00437333, 0.00437333],
        crota = 0,
        equinox = 2000,
        epoch = 2000,
    )


@pytest.mark.parametrize("setup", [setup_linear, setup_wcs, setup_tanp])
def test_pix2world_simple(setup):
    """Simple check for pix2world"""

    opts = setup()
    args = [opts.crpix, opts.crval, opts.cdelt,
            opts.crota, opts.equinox, opts.epoch]

    x0 = opts.x0.copy()
    x1 = opts.x1.copy()

    w0, w1 = _wcs.pix2world(opts.ctype, x0, x1, *args)
    assert w0 == pytest.approx(opts.w0)
    assert w1 == pytest.approx(opts.w1)

    # Are the original arguments changed? This is a regression test
    # and so we can change the behaviour.
    assert x0 == pytest.approx(w0)
    assert x1 == pytest.approx(w1)


# XFAIL: world2pix isn't set up for TANP
@pytest.mark.parametrize("setup", [setup_linear, setup_wcs, pytest.param(setup_tanp, marks=pytest.mark.xfail)])
def test_world2pix_simple(setup):
    """Simple check for world2pix"""

    opts = setup()
    args = [opts.crpix, opts.crval, opts.cdelt,
            opts.crota, opts.equinox, opts.epoch]

    w0 = opts.w0.copy()
    w1 = opts.w1.copy()

    x0, x1 = _wcs.world2pix(opts.ctype, w0, w1, *args)
    assert x0 == pytest.approx(opts.x0)
    assert x1 == pytest.approx(opts.x1)

    # Are the original arguments changed? This is a regression test
    # and so we can change the behaviour.
    assert w0 == pytest.approx(x0)
    assert w1 == pytest.approx(x1)


@pytest.mark.parametrize("setup", [setup_linear, setup_wcs, setup_tanp])
def test_pix2world_argument_change_when_integer(setup):
    """Corner-case handling

    There is no reason for this other than it is the
    current behaviour, so add a regression test.
    """

    opts = setup()
    args = [opts.crpix, opts.crval, opts.cdelt,
            opts.crota, opts.equinox, opts.epoch]

    x0 = np.asarray(opts.x0, dtype=int)
    x1 = np.asarray(opts.x1, dtype=int)
    x0_ = x0.copy()
    x1_ = x1.copy()

    w0, w1 = _wcs.pix2world(opts.ctype, x0, x1, *args)
    assert x0 == pytest.approx(x0_)
    assert x1 == pytest.approx(x1_)


@pytest.mark.parametrize("setup", [setup_linear, setup_wcs, setup_tanp])
def test_world2pix_argument_change_when_integer(setup):
    """Corner-case handling

    There is no reason for this other than it is the
    current behaviour, so add a regression test.
    """

    opts = setup()
    args = [opts.crpix, opts.crval, opts.cdelt,
            opts.crota, opts.equinox, opts.epoch]

    w0 = np.asarray(opts.w0, dtype=int)
    w1 = np.asarray(opts.w1, dtype=int)
    w0_ = w0.copy()
    w1_ = w1.copy()

    x0, x1 = _wcs.world2pix(opts.ctype, w0, w1, *args)
    assert w0 == pytest.approx(w0_)
    assert w1 == pytest.approx(w1_)
