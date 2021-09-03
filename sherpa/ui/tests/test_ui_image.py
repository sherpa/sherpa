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

"""Basic tests of the image functionality in sherpa.ui. The hope is that we
can run some tests even without DS9/XPA available. Although tests are
run with both the base and Astro session classes, there's no tests of
Astro-specific functionality (e.g. DataIMG).

Notes
-----
A number of tests send data to an external DS9 process (which will be
started if needed) and then check the values displayed by the process.
These tests are highly-likely to fail if the tests are run in
parallel, such as with

    pytest -n auto

as there is no support in our DS9 interface to either create a new DS9
instance per worker, or to "lock" a DS9 so that the tests can be run
in serial.

"""

import logging

import numpy as np

import pytest

from sherpa.ui.utils import Session as BaseSession
from sherpa.astro.ui.utils import Session as AstroSession

from sherpa.data import Data2D
from sherpa.models import basic

from sherpa.stats import Chi2Gehrels
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, IdentifierErr
from sherpa.utils.testing import requires_ds9


def example_data():
    """Create an example data set."""

    x1, x0 = np.mgrid[-4:5, 6:15]

    shp = x1.shape

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / (1 + np.sqrt((x0 - 11)**2 + (x1 + 2)**2))

    return Data2D('example', x0, x1, y, shape=shp)


def example_model():
    """Create an example data set."""

    gmdl = basic.Gauss2D('gmdl')
    gmdl.xpos = 9
    gmdl.ypos = -1
    gmdl.fwhm = 3
    gmdl.ampl = 100

    cmdl = basic.Const2D()
    cmdl.c0 = 2

    return gmdl + cmdl


def example_psf():
    """Create an example PSF."""

    # Note that the PSF coordinates have to match the data.
    #
    psf = basic.Box2D()
    psf.xlow = 8
    psf.xhi = 12
    psf.ylow = 0
    psf.yhi = 3
    return psf


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("shape", [None, (9, 9)])
def test_load_arrays2d(session, shape):
    """Does load_arrays work with 2D data?"""

    d = example_data()

    s = session()
    if shape is None:
        s.load_arrays(1, d.x0, d.x1, d.y, Data2D)
    else:
        s.load_arrays(1, d.x0, d.x1, d.y, shape, Data2D)

    got = s.get_data()
    assert isinstance(got, Data2D)

    # use pytest.approx even though expect equality here as it's
    # easier to read the output if there's an error rather than
    # np.all(a == b)
    #
    assert got.x0 == pytest.approx(d.x0)
    assert got.x1 == pytest.approx(d.x1)
    assert got.y == pytest.approx(d.y)
    if shape is None:
        assert got.shape is None
    else:
        assert got.shape == pytest.approx(shape)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_data_image(session):
    from sherpa.image import DataImage

    s = session()
    d = example_data()
    s.set_data(d)

    obj = s.get_data_image()
    assert isinstance(obj, DataImage)
    assert obj.name == 'Data'
    assert obj.eqpos is None
    assert obj.sky is None

    y = d.y.reshape((9, 9))
    assert obj.y == pytest.approx(y)

    # sanity check
    assert y[2, 5] == pytest.approx(100.0)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_model_image(session):
    from sherpa.image import ModelImage

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    obj = s.get_model_image()
    assert isinstance(obj, ModelImage)
    assert obj.name == 'Model'
    assert obj.eqpos is None
    assert obj.sky is None

    y = d.eval_model(m).reshape((9, 9))
    assert obj.y == pytest.approx(y)

    # sanity check
    assert y[3, 3] == pytest.approx(102.0)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_model_component_image(session):
    from sherpa.image import ComponentModelImage

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    # rely on m being (gmdl + cmdl) and we want the
    # gaussian component.
    #
    obj = s.get_model_component_image(m.parts[0])
    assert isinstance(obj, ComponentModelImage)
    assert obj.name == 'Model_component'
    assert obj.eqpos is None
    assert obj.sky is None

    y = d.eval_model(m.parts[0]).reshape((9, 9))
    assert obj.y == pytest.approx(y)

    # sanity check
    assert y[3, 3] == pytest.approx(100.0)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_model_component_image_with_convolution(session):
    """What happens if there's a convolution component in play?"""

    from sherpa.image import ComponentModelImage

    s = session()
    d = example_data()
    m = example_model()
    p = example_psf()

    s.set_data(d)
    s.set_source(m)
    s.load_psf('psf', p)
    s.set_psf('psf')

    # rely on m being (gmdl + cmdl) and we want the
    # gaussian component.
    #
    obj = s.get_model_component_image(m.parts[0])
    assert isinstance(obj, ComponentModelImage)
    assert obj.name == 'Model_component'
    assert obj.eqpos is None
    assert obj.sky is None

    # Note that the convolution component is applied
    # to the data.
    psf = s.get_model_component('psf')
    y = d.eval_model(psf(m.parts[0])).reshape((9, 9))
    assert obj.y == pytest.approx(y)

    # sanity check
    assert y[3, 3] == pytest.approx(42.25302)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_source_image(session):
    """Note that here there's no difference of source and model"""
    from sherpa.image import SourceImage

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    obj = s.get_source_image()
    assert isinstance(obj, SourceImage)
    assert obj.name == 'Source'
    assert obj.eqpos is None
    assert obj.sky is None

    y = d.eval_model(m).reshape((9, 9))
    assert obj.y == pytest.approx(y)

    # sanity check
    assert y[3, 3] == pytest.approx(102.0)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_source_component_image(session):
    """Note that here there's no difference of source and model"""
    from sherpa.image import ComponentSourceImage

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    obj = s.get_source_component_image(m.parts[0])
    assert isinstance(obj, ComponentSourceImage)
    assert obj.name == 'Source_component'
    assert obj.eqpos is None
    assert obj.sky is None

    y = d.eval_model(m.parts[0]).reshape((9, 9))
    assert obj.y == pytest.approx(y)

    # sanity check
    assert y[3, 3] == pytest.approx(100.0)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_source_component_image_with_convolution(session):
    """There is a difference thanks to the PSF"""
    from sherpa.image import ComponentSourceImage

    s = session()
    d = example_data()
    m = example_model()
    p = example_psf()
    s.set_data(d)
    s.set_source(m)
    s.load_psf('psf', p)
    s.set_psf(psf)

    obj = s.get_source_component_image(m.parts[0])
    assert isinstance(obj, ComponentSourceImage)
    assert obj.name == 'Source_component'
    assert obj.eqpos is None
    assert obj.sky is None

    # Note that the PSF convolution is not applied here
    y = d.eval_model(m.parts[0]).reshape((9, 9))
    assert obj.y == pytest.approx(y)

    # sanity check
    assert y[3, 3] == pytest.approx(100.0)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_resid_image(session):
    from sherpa.image import ResidImage

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    obj = s.get_resid_image()
    assert isinstance(obj, ResidImage)
    assert obj.name == 'Residual'
    assert obj.eqpos is None
    assert obj.sky is None

    y = (d.y - d.eval_model(m)).reshape((9, 9))
    assert obj.y == pytest.approx(y)

    # sanity check the two peaks: data and model
    assert y[2, 5] == pytest.approx(76.5689)
    assert y[3, 3] == pytest.approx(-71.0983)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_ratio_image(session):
    from sherpa.image import RatioImage

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    obj = s.get_ratio_image()
    assert isinstance(obj, RatioImage)
    assert obj.name == 'Ratio'
    assert obj.eqpos is None
    assert obj.sky is None

    y = (d.y / d.eval_model(m)).reshape((9, 9))
    assert obj.y == pytest.approx(y)

    # sanity check the two peaks: data and model
    assert y[2, 5] == pytest.approx(4.26783)
    assert y[3, 3] == pytest.approx(0.302958)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_psf_image(session):
    from sherpa.image import PSFImage

    s = session()
    d = example_data()
    m = example_model()
    p = example_psf()
    s.set_data(d)
    s.set_source(m)
    s.load_psf('pmdl', p)
    s.set_psf('pmdl')

    obj = s.get_psf_image()
    assert isinstance(obj, PSFImage)
    assert obj.name == 'box2d'
    assert obj.eqpos is None
    assert obj.sky is None

    y = np.zeros((9, 9))
    y[5:7, 3:6] = 1
    assert obj.y == pytest.approx(y)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_kernel_image(session):
    from sherpa.image import PSFKernelImage

    s = session()
    d = example_data()
    m = example_model()
    p = example_psf()
    s.set_data(d)
    s.set_source(m)
    s.load_psf('pmdl', p)
    s.set_psf('pmdl')

    obj = s.get_kernel_image()
    assert isinstance(obj, PSFKernelImage)
    assert obj.name == 'PSF_Kernel'
    assert obj.eqpos is None
    assert obj.sky is None

    # Unlike the PSF, this is normalized to 1 by default
    y = np.zeros((9, 9))
    y[5:7, 3:6] = 1 / 6
    assert obj.y == pytest.approx(y)


def check_xpa(backend, command, expected):
    """Do we get the expected response?"""

    ans = backend.xpaget(command)
    assert ans == expected


def check_xpa_data(backend):
    grid = '30.9017\n100\n50\n50\n33.3333\n41.4214\n'
    check_xpa(backend, 'data image 6 3 3 2 yes', grid)


def check_xpa_model(backend):
    grid = '2.72334\n23.4311\n6.59292\n31.1632\n2.53156\n8.25\n'
    check_xpa(backend, 'data image 6 3 3 2 yes', grid)


def check_xpa_model_component(backend):
    # This is check_xpa_model values - 2 (and a tweak because of
    # the number of dp in one bin)
    grid = '0.72334\n21.4311\n4.59292\n29.1632\n0.531559\n6.25\n'
    check_xpa(backend, 'data image 6 3 3 2 yes', grid)


def check_xpa_model_component_with_convolution(backend):
    # This is check_xpa_model values - 2 and then convolved by
    # the PSF (which is purposefully strange)
    grid = '1.20076\n6.4275\n2.13285\n18.631\n0.414249\n6.18236\n'
    check_xpa(backend, 'data image 6 3 3 2 yes', grid)


def check_xpa_resid(backend):
    grid = '28.1784\n76.5689\n43.4071\n18.8368\n30.8018\n33.1714\n'
    check_xpa(backend, 'data image 6 3 3 2 yes', grid)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_data(session):
    from sherpa.image import backend

    s = session()
    d = example_data()
    s.set_data(d)

    s.image_data()

    check_xpa(backend, 'tile', 'no\n')
    check_xpa_data(backend)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_model(session):
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    s.image_model()

    check_xpa(backend, 'tile', 'no\n')
    check_xpa_model(backend)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_model_component(session):
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    s.image_model_component(m.parts[0])

    check_xpa(backend, 'tile', 'no\n')
    check_xpa_model_component(backend)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_model_component_with_convolution(session):
    """The convolution component changes the data."""
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    p = example_psf()
    s.set_data(d)
    s.set_source(m)
    s.load_psf('psf', p)
    s.set_psf(psf)

    s.image_model_component(m.parts[0])

    check_xpa(backend, 'tile', 'no\n')
    check_xpa_model_component_with_convolution(backend)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_source(session):
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    s.image_source()

    check_xpa(backend, 'tile', 'no\n')

    # In this case the values are the same as test_image_model
    check_xpa_model(backend)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_source_component(session):
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    s.image_source_component(m.parts[0])

    check_xpa(backend, 'tile', 'no\n')

    # In this case the values are the same as test_image_model_component
    check_xpa_model_component(backend)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_source_component_with_convolution(session):
    """Unlike model, this does not change the component."""
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    p = example_psf()
    s.set_data(d)
    s.set_source(m)
    s.load_psf('psf', p)
    s.set_psf(psf)

    s.image_source_component(m.parts[0])

    check_xpa(backend, 'tile', 'no\n')
    check_xpa_model_component(backend)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_resid(session):
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    s.image_resid()

    check_xpa(backend, 'tile', 'no\n')
    check_xpa_resid(backend)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_ratio(session):
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    s.image_ratio()

    check_xpa(backend, 'tile', 'no\n')

    grid = '11.347\n4.26783\n7.58389\n1.60446\n13.1671\n5.02077\n'
    check_xpa(backend, 'data image 6 3 3 2 yes', grid)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_fit(session):
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    s.set_data(d)
    s.set_source(m)

    s.image_fit()

    # Technically the frame values could differ depending on what
    # previous checks have been run.
    #
    check_xpa(backend, 'tile', 'yes\n')
    check_xpa(backend, 'frame', '3\n')
    check_xpa(backend, 'frame active', '1 2 3 \n')

    # Check frame 1 - data
    #       frame 2 - model
    #       frame 3 - resif
    backend.xpaset('frame 1')
    check_xpa_data(backend)
    backend.xpaset('frame 2')
    check_xpa_model(backend)
    backend.xpaset('frame 3')
    check_xpa_resid(backend)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_psf(session):
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    p = example_psf()

    s.set_data(d)
    s.set_source(m)
    s.load_psf('pmdl', p)
    s.set_psf('pmdl')
    s.image_psf()

    check_xpa(backend, 'tile', 'no\n')

    # Check a different region than other tests
    grid = '0\n1\n0\n1\n1\n0\n0\n0\n1\n1\n0\n0\n1\n0\n0\n0\n'
    check_xpa(backend, 'data image 4 5 4 4 yes', grid)


@requires_ds9
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_image_kernel(session):
    from sherpa.image import backend

    s = session()
    d = example_data()
    m = example_model()
    p = example_psf()

    s.set_data(d)
    s.set_source(m)
    s.load_psf('pmdl', p)
    s.set_psf('pmdl')
    s.image_kernel()

    check_xpa(backend, 'tile', 'no\n')

    # Check the same reagion as test_image_psf
    grid = '0\n0.166667\n0\n0.166667\n0.166667\n0\n0\n0\n0.166667\n0.166667\n0\n0\n0.166667\n0\n0\n0\n'
    check_xpa(backend, 'data image 4 5 4 4 yes', grid)
