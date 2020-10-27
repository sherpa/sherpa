#
#  Copyright (C) 2019, 2020  Smithsonian Astrophysical Observatory
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
These tests exercise an issue with Sherpa summarized in
https://github.com/sherpa/sherpa/issues/43

A very simple test case which would capture this is the following.

Assume you have a 2D source with a simple, circular, Gaussian emission. Assume you have a PSF defined at a
higher resolution than the images to fit will have. And assume that this PSF is also a simple, circular
Gaussian with different parameters than the source.

We can create a simulated image IMG by convolving the model M, evaluated at the PSF's resolution,
with the PSF, and then rebin it at the desired image resolution to obtain a new image IMG_lo. As the
convolution is between two Gaussians, the resulting image will also be the image of a Gaussian,
whose values we can predict, and in particular the variance of IMG will be the sum in quadrature of the
variances of M and PSF. For simplicity, we might need to assume an exact ratio in pixel sizes.

When the current Sherpa is run, it will apply the PSF to the image as if they were defined at the same
resolution. Sherpa will overestimate the size of the PSF and thus underestimate the size of the source.
We can analytically derive the (incorrect) fit results Sherpa would calculate now, while
the correct results we should get after the improvements are already known by us setting up the image.

The tests in this module implement the following scenario.

We will keep updating the file with more tests for error handling, edge and corner cases, etc.

We might have a different file for integration tests with more realistic data, as we will need
to exercise actual images and PSFs with WCS headers.
"""

from math import sqrt

import attr
import numpy as np

from pytest import approx, fixture, mark

from sherpa.astro.ui.utils import Session
from sherpa.astro.data import DataIMG
from sherpa.astro.instrument import PSFModel
from sherpa.instrument import PSFSpace2D
from sherpa.models import SigmaGauss2D
from sherpa.models.regrid import EvaluationSpace2D


DATA_PIXEL_SIZE = 2


@attr.s
class FixtureConfiguration():
    image_size = attr.ib()
    psf_size = attr.ib()
    source_amplitude = attr.ib()
    source_sigma = attr.ib()
    psf_sigma = attr.ib()
    resolution_ratio = attr.ib()
    psf_amplitude = attr.ib(default=1)
    image_resolution = attr.ib(default=1)

    def __attrs_post_init__(self):
        self.source_position = self.image_size / 2
        self.psf_position = self.psf_size / 2


@attr.s
class FixtureData():
    image = attr.ib()
    psf = attr.ib()
    psf_model = attr.ib()
    configuration = attr.ib()

    def __attrs_post_init__(self):
        self.data_space = EvaluationSpace2D(*self.image.get_indep(filter=False))


def generate_psf_space_configurations():
    return (FixtureConfiguration(image_size=500, psf_size=100, source_amplitude=100,
                                 source_sigma=50, psf_sigma=5, resolution_ratio=1.5),
            FixtureConfiguration(image_size=300, psf_size=150, source_amplitude=100,
                                 source_sigma=1, psf_sigma=15, resolution_ratio=1.5),
            FixtureConfiguration(image_size=500, psf_size=100, source_amplitude=100,
                                 source_sigma=50, psf_sigma=5, resolution_ratio=2.5),
            FixtureConfiguration(image_size=300, psf_size=150, source_amplitude=100,
                                 source_sigma=1, psf_sigma=15, resolution_ratio=2.5),
            FixtureConfiguration(image_size=500, psf_size=100, source_amplitude=100,
                                 source_sigma=50, psf_sigma=5, resolution_ratio=0.5),
            FixtureConfiguration(image_size=300, psf_size=150, source_amplitude=100,
                                 source_sigma=1, psf_sigma=15, resolution_ratio=0.5),
            FixtureConfiguration(image_size=500, psf_size=100, source_amplitude=100,
                                 source_sigma=50, psf_sigma=5, resolution_ratio=0.7),
            FixtureConfiguration(image_size=300, psf_size=150, source_amplitude=100,
                                 source_sigma=1, psf_sigma=15, resolution_ratio=0.7),
            FixtureConfiguration(image_size=500, psf_size=100, source_amplitude=100,
                                 source_sigma=50, psf_sigma=5, resolution_ratio=1),
            FixtureConfiguration(image_size=300, psf_size=150, source_amplitude=100,
                                 source_sigma=1, psf_sigma=15, resolution_ratio=2),
            FixtureConfiguration(image_size=500, psf_size=100, source_amplitude=100,
                                 source_sigma=50, psf_sigma=5, resolution_ratio=3),
            FixtureConfiguration(image_size=300, psf_size=150, source_amplitude=100,
                                 source_sigma=1, psf_sigma=15, resolution_ratio=4),
            FixtureConfiguration(image_size=500, psf_size=100, source_amplitude=100,
                                 source_sigma=50, psf_sigma=5, resolution_ratio=5),
            FixtureConfiguration(image_size=300, psf_size=150, source_amplitude=100,
                                 source_sigma=1, psf_sigma=15, resolution_ratio=6),
            FixtureConfiguration(image_size=500, psf_size=100, source_amplitude=100,
                                 source_sigma=50, psf_sigma=5, resolution_ratio=7),
            FixtureConfiguration(image_size=300, psf_size=150, source_amplitude=100,
                                 source_sigma=1, psf_sigma=15, resolution_ratio=8),
            )


def generate_rebinning_configurations():
    return (FixtureConfiguration(image_size=128, psf_size=64, source_amplitude=100,
                                 source_sigma=10, psf_sigma=5, resolution_ratio=2),
            FixtureConfiguration(image_size=64, psf_size=128, source_amplitude=100,
                                 source_sigma=10, psf_sigma=5, resolution_ratio=2),
            FixtureConfiguration(image_size=256, psf_size=128, source_amplitude=100,
                                 source_sigma=5, psf_sigma=5, resolution_ratio=3),
            FixtureConfiguration(image_size=128, psf_size=64, source_amplitude=100,
                                 source_sigma=10, psf_sigma=5, resolution_ratio=3),
            FixtureConfiguration(image_size=64, psf_size=512, source_amplitude=100,
                                 source_sigma=10, psf_sigma=20, resolution_ratio=4),
            FixtureConfiguration(image_size=128, psf_size=512, source_amplitude=100,
                                 source_sigma=5, psf_sigma=20, resolution_ratio=4),
            FixtureConfiguration(image_size=128, psf_size=64, source_amplitude=100,
                                 source_sigma=10, psf_sigma=5, resolution_ratio=1),
            FixtureConfiguration(image_size=64, psf_size=128, source_amplitude=100,
                                 source_sigma=10, psf_sigma=5, resolution_ratio=1),
            FixtureConfiguration(image_size=128, psf_size=64, source_amplitude=100,
                                 source_sigma=5, psf_sigma=5, resolution_ratio=1),
            FixtureConfiguration(image_size=128, psf_size=64, source_amplitude=100,
                                 source_sigma=10, psf_sigma=5, resolution_ratio=1),
            FixtureConfiguration(image_size=64, psf_size=512, source_amplitude=100,
                                 source_sigma=10, psf_sigma=5, resolution_ratio=1),
            FixtureConfiguration(image_size=128, psf_size=512, source_amplitude=100,
                                 source_sigma=5, psf_sigma=5, resolution_ratio=1),
            )


@mark.parametrize("psf_fixture", generate_rebinning_configurations(), indirect=True)
def test_psf_resolution_bug(psf_fixture):
    session, source, expected_sigma = psf_fixture
    session.fit()

    assert source.sigma_a.val == expected_sigma
    assert source.sigma_b.val == expected_sigma


@mark.parametrize("configuration", generate_psf_space_configurations())
def test_psf_space(configuration):

    fixture_data = make_images(configuration)

    psf_space = PSFSpace2D(fixture_data.data_space, fixture_data.psf_model, data_pixel_size=fixture_data.image.sky.cdelt)

    # Test that the PSF space has the same boundaries of the data space, but a number of
    # bins equal to the bins the data image would have if it had the same pixel size as the PSF.
    # Also, we want this space to be larger, if not the same, than the data space, to avoid introducing
    # more boundary effects in the convolution.
    n_bins = configuration.image_size
    assert psf_space.start == (0, 0)
    assert psf_space.x_axis.size == psf_space.y_axis.size == n_bins * configuration.resolution_ratio
    assert psf_space.end == (approx(n_bins - 1 / configuration.resolution_ratio),
                             approx(n_bins - 1 / configuration.resolution_ratio))


def test_rebin_int_no_int():
    """
    Test that the correct rebin_* function is called depending on whether the pixel size is an integer or close to
    an integer with a tolerance that can be changed by the user.
    """
    from sherpa.models.regrid import rebin_2d
    from unittest import mock

    rebin_int = mock.MagicMock()
    rebin_no_int = mock.MagicMock()

    to_space = mock.MagicMock()
    y = mock.MagicMock()

    with mock.patch('sherpa.models.regrid.rebin_int', rebin_int):
        # The pixel ratio is 2, perfect integer, rebin_int should be called
        from_space = mock.MagicMock(data_2_psf_pixel_size_ratio=(0.5, 0.5))
        rebin_2d(y, from_space, to_space)
        assert rebin_int.called
        rebin_int.reset_mock()

    with mock.patch('sherpa.models.regrid.rebin_int', rebin_int):
        # Not a perfect integer, but close to an integer within the default tolerance
        from_space = mock.MagicMock(data_2_psf_pixel_size_ratio=(0.333, 0.333))
        rebin_2d(y, from_space, to_space)
        assert rebin_int.called
        rebin_int.reset_mock()

    with mock.patch('sherpa.models.regrid.rebin_no_int', rebin_no_int):
        # Same case as above, but I am changing the tolerance so the ratio is not equal to an integer within
        # the new tolerance.
        from_space = mock.MagicMock(data_2_psf_pixel_size_ratio=(0.333, 0.333))
        from sherpa.models import regrid
        regrid.PIXEL_RATIO_THRESHOLD = 0.000001
        rebin_2d(y, from_space, to_space)
        assert rebin_no_int.called


def symmetric_gaussian_image(amplitude, sigma, position, n_bins):
    model = SigmaGauss2D()
    model.ampl = amplitude
    model.sigma_a = sigma
    model.sigma_b = sigma
    model.xpos = position
    model.ypos = position

    arrays = np.arange(0, n_bins), np.arange(0, n_bins)
    x_array, y_array = np.meshgrid(*arrays)
    x_array, y_array = x_array.flatten(), y_array.flatten()

    return model(x_array, y_array).flatten(), x_array, y_array


def make_images(configuration):

    psf, psf_model = make_psf(configuration)
    data_image = make_image(configuration)

    return FixtureData(data_image, psf, psf_model, configuration)


def make_image(configuration):
    source_position = configuration.source_position

    # The convolution of two gaussians is a gaussian with a stddev which is the sum
    # of those convolved, so we model the image as a Gaussian itself.
    image_sigma = sqrt(configuration.source_sigma ** 2 + configuration.psf_sigma ** 2)
    image_amplitude = configuration.source_amplitude * configuration.psf_amplitude
    image, image_x, image_y = symmetric_gaussian_image(amplitude=image_amplitude, sigma=image_sigma,
                                                       position=source_position, n_bins=configuration.image_size)

    data_image = DataIMG("image", image_x, image_y, image,
                         shape=(configuration.image_size, configuration.image_size),
                         sky=WcsStub([DATA_PIXEL_SIZE, DATA_PIXEL_SIZE])
                         )

    return data_image


def make_psf(configuration):
    # The psf parameters in terms of the input parameters. To simulate the different pixel size
    # we set the sigma of the psf to be a multiple (according to the ratio input)
    # of the actual sigma in data pixel units.
    # This should correspond, when ratio>1 to the case when the PSF's resolution is bigger
    # than the image. If the ratio is 1 than they have the same pixel size.
    psf_sigma = configuration.psf_sigma * configuration.resolution_ratio
    psf_amplitude = configuration.psf_amplitude
    psf_position = configuration.psf_position

    # We model the PSF as a gaussian as well. Note we are using the psf_sigma variable,
    # that is the one which is scaled through the ration. In other terms, we are working in
    # units of Data Pixels to simulate the conditions of the bug, when the ratio != 1.
    psf, psf_x, psf_y = symmetric_gaussian_image(amplitude=psf_amplitude, sigma=psf_sigma,
                                                 position=psf_position, n_bins=configuration.psf_size)
    # Normalize PSF
    norm_psf = psf / psf.sum()

    cdelt = DATA_PIXEL_SIZE / configuration.resolution_ratio
    psf_wcs = WcsStub([cdelt, cdelt])

    # Create a Sherpa PSF model object using the psf arrays
    sherpa_kernel = DataIMG('kernel_data',
                            psf_x, psf_y, norm_psf,
                            shape=(configuration.psf_size, configuration.psf_size),
                            sky=psf_wcs
                            )

    psf_model = PSFModel('psf_model', kernel=sherpa_kernel)
    psf_model.norm = 1
    psf_model.origin = (psf_position + 1, psf_position + 1)

    return psf, psf_model


@fixture
def psf_fixture(request):
    configuration = request.param

    fixture_data = make_images(configuration)

    ui = Session()

    ui.set_data(1, fixture_data.image)

    exact_expected_sigma = configuration.source_sigma
    approx_expected_sigma = approx(exact_expected_sigma, rel=3e-2)

    # Set the source model as a 2D Gaussian, and set the PSF in Sherpa
    source_position = configuration.source_position
    sherpa_source = SigmaGauss2D('source')
    sherpa_source.ampl = configuration.source_amplitude
    sherpa_source.sigma_a = exact_expected_sigma
    sherpa_source.sigma_b = exact_expected_sigma
    sherpa_source.xpos = source_position
    sherpa_source.ypos = source_position
    ui.set_source(sherpa_source)
    ui.set_psf(fixture_data.psf_model)

    return ui, sherpa_source, approx_expected_sigma


@attr.s
class WcsStub():
    cdelt = attr.ib()
