# Copyright 2018 Smithsonian Astrophysical Observatory
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
# following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
# disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
# products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
import pytest

from pytest import approx

from sherpa.models import Const2D
from sherpa.astro.ui.utils import Session
from sherpa.utils.testing import requires_data, requires_fits


@pytest.fixture
def setup(make_data_path):
    image = make_data_path('sim_0.5bin_2g+c.fits')
    psf_bin1 = make_data_path('psf_0.0_00_bin1.img')
    psf_bin05 = make_data_path('psf_0.0_00_bin0.5.img')

    return Session(), image, psf_bin05, psf_bin1


@requires_fits
@requires_data
def test_psf_rebin_warning(setup):
    ui, image, _, psf_bin1 = setup

    ui.load_image(image)
    ui.load_psf('psf', psf_bin1)

    with pytest.raises(AttributeError):
        ui.set_psf('psf')


@requires_fits
@requires_data
def test_psf_rebin_no_warning(setup):
    ui, image, psf_bin05, _ = setup

    ui.load_image(image)
    ui.load_psf('psf', psf_bin05)
    ui.set_psf('psf')


# https://github.com/gammapy/gammapy/issues/1905
@requires_fits
@requires_data
def test_psf_sky_none(setup):
    ui, image, psf_bin05, _ = setup

    ui.load_image(image)
    ui.load_psf('psf', psf_bin05)

    ui.set_psf('psf')

    ui.set_model(Const2D("m"))

    ui.get_psf().kernel.sky = None

    # In https://github.com/gammapy/gammapy/issues/1905 this would fail because of a call to kernel.sky.cdelt
    with pytest.warns(UserWarning):
        ui.fit()


@requires_fits
@requires_data
def test_calc_source_model_sum2d(setup):
    ui, image, psf_bin05, _ = setup

    ui.load_image(image)
    ui.get_data().sky.cdelt = (1, 1)

    ui.load_psf('psf', psf_bin05)
    ui.set_psf('psf')

    ui.set_source(Const2D())
    ui.fit()
    assert ui.calc_source_sum2d() == approx(ui.calc_model_sum2d(), rel=0.03)
