#
# Copyright (C) 2020, 2021  Smithsonian Astrophysical Observatory
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

"""Very-basic tests of the HTML representation of objects.

"""

import numpy as np

import pytest

from sherpa.astro import data
from sherpa import plot
from sherpa.astro.instrument import create_delta_rmf
from sherpa.utils.testing import requires_data, requires_fits


TEST_HEADER = {}
TEST_HEADER['OBJECT'] = 'The best object'
TEST_HEADER['TELESCOP'] = 'Chandra'
TEST_HEADER['HDUCLAS3'] = 'test for RMF'
TEST_HEADER['HDUCLAS4'] = 'test for PHA'
TEST_HEADER['EXPOSURE'] = 123
TEST_HEADER['LO_THRES'] = 0.001


def check(r, summary, name, label, nmeta):
    """Very limited checks of the structure"""
    assert r is not None

    if plot.backend.name == 'pylab':
        assert '<summary>{} Plot</summary>'.format(summary) in r
        assert '<svg ' in r

    else:
        assert '<summary>{} Data ('.format(summary) in r
        assert '<div class="dataname">{}</div>'.format(label) in r
        assert '<svg ' not in r

    assert '<summary>Summary (' in r
    if nmeta == 0:
        assert '<summary>Metadata' not in r
    else:
        assert '<summary>Metadata ({})'.format(nmeta) in r

    assert '<div class="dataval">{}</div>'.format(name) in r


@pytest.mark.parametrize('header', [None, {}, TEST_HEADER])
def test_arf(header, override_plot_backend):
    ebins = np.arange(0.1, 1, 0.1)

    # fake wavelength conversion
    wbins = 10 / ebins

    y = np.linspace(10, 80, ebins.size - 1)
    d = data.DataARF('arf 1', ebins[:-1], ebins[1:], y,
                     bin_lo=wbins[1:], bin_hi=wbins[:-1],
                     exposure=100.1, header=header)
    r = d._repr_html_()

    nmeta = 0 if (header is None) or (header == {}) else 2
    check(r, 'ARF', 'arf 1', 'SPECRESP', nmeta=nmeta)

    assert '<div class="dataname">Energy range</div>' in r
    assert '<div class="dataname">Wavelength range</div>' in r
    assert '<div class="dataval">0.1 - 0.9 keV, bin size 0.1 keV</div>' in r
    assert '<div class="dataval">12.5 - 50 &#8491;, bin size 1.38889 - 50 &#8491;</div>' in r


@pytest.mark.parametrize('header', [None, {}, TEST_HEADER])
def test_rmf(header, override_plot_backend):
    ebins = np.arange(0.1, 1, 0.1)
    elo, ehi = ebins[:-1], ebins[1:]
    d = create_delta_rmf(elo, ehi, header=header)
    r = d._repr_html_()

    nmeta = 0 if (header is None) or (header == {}) else 4
    check(r, 'RMF', 'delta-rmf', 'MATRIX', nmeta=nmeta)



@pytest.mark.parametrize('header', [None, {}, TEST_HEADER])
def test_pha(header, override_plot_backend):
    d = data.DataPHA('x x',
                     np.arange(1, 5, dtype=np.int16),
                     np.arange(1, 5, dtype=np.int16),
                     header=header)
    r = d._repr_html_()
    if header is None:
        nmeta = 5
    elif header == {}:
        nmeta = 0
    else:
        nmeta = 4
    check(r, 'PHA', 'x x', 'COUNTS', nmeta=nmeta)


@requires_data
@requires_fits
@pytest.mark.parametrize("subtract", [False, True])
@pytest.mark.parametrize("group", [False, True])
def test_pha_real(subtract, group, make_data_path, override_plot_backend):
    """This is grouped and has a background"""
    from sherpa.astro.io import read_pha
    d = read_pha(make_data_path('3c273.pi'))
    d.name = '3c273.pi'

    if not group:
        d.ungroup()

    d.notice(0.5, 7)
    if subtract:
        d.subtract()

    r = d._repr_html_()

    check(r, 'PHA', '3c273.pi', 'COUNTS', nmeta=11)

    assert '<div class="dataval">2000-01-10T06:47:15</div>' in r
    assert '<div class="dataname">Background</div>' in r
    assert '<div class="dataname">Grouping</div>' in r

    if subtract:
        assert '<div class="dataval">Subtracted</div>' in r
    else:
        assert '<div class="dataval">Not subtracted</div>' in r

    if group:
        assert '<div class="dataval">Applied (46 groups)</div>' in r
        assert '<div class="dataval">0.4672-9.8696 Energy (keV) with 42 groups</div>' in r
    else:
        assert '<div class="dataval">Not applied</div>' in r
        assert '<div class="dataval">0.4964-7.0080 Energy (keV) with 446 channels</div>' in r


@requires_data
@requires_fits
def test_arf_real(make_data_path, override_plot_backend):
    """Read in an ARF"""
    from sherpa.astro.io import read_arf
    d = read_arf(make_data_path('9774.arf'))
    d.name = '9774.arf'

    r = d._repr_html_()

    check(r, 'ARF', '9774.arf', 'Exposure', nmeta=6)

    assert '<div class="dataname">Identifier</div><div class="dataval">9774.arf</div>' in r
    assert '<div class="dataname">Exposure</div><div class="dataval">75141.2 s</div>' in r
    assert '<div class="dataname">Number of bins</div><div class="dataval">1078</div>' in r
    assert '<div class="dataname">Energy range</div><div class="dataval">0.22 - 11 keV, bin size 0.0100002 keV</div>' in r
    assert '<div class="dataname">Area range</div><div class="dataval">0.533159 - 681.944 cm<sup>2</sup></div>' in r

    assert '<div class="dataname">Object</div><div class="dataval">3C 186</div>' in r
    assert '<div class="dataname">Observation date</div><div class="dataval">2007-12-06T05:30:00</div>' in r


@requires_data
@requires_fits
def test_rmf_real(make_data_path, override_plot_backend):
    """Read in a RMF"""
    from sherpa.astro.io import read_rmf
    d = read_rmf(make_data_path('9774.rmf'))
    d.name = '9774.rmf'

    r = d._repr_html_()

    check(r, 'RMF', '9774.rmf', 'DETCHANS', nmeta=6)

    assert '<div class="dataname">Identifier</div><div class="dataval">9774.rmf</div>' in r
    assert '<div class="dataname">Number of channels</div><div class="dataval">1024</div>' in r
    assert '<div class="dataname">Number of energies</div><div class="dataval">1078</div>' in r
    assert '<div class="dataname">Energy range</div><div class="dataval">0.22 - 11 keV, bin size 0.0100002 keV</div>' in r
    assert '<div class="dataname">Channel range</div><div class="dataval">1 - 1024</div>' in r

    assert '<div class="dataname">The channel type</div><div class="dataval">PI</div>' in r
    assert '<div class="dataname">Matrix contents</div><div class="dataval">REDIST</div>' in r


@pytest.mark.parametrize('header', [None, {}, TEST_HEADER])
def test_img(header,
             override_plot_backend, old_numpy_printing):
    y, x = np.mgrid[1:4, 2:4]
    z = np.arange(x.size)
    d = data.DataIMG('x x', x.flatten(), y.flatten(), z,
                     shape=x.shape, header=header)
    r = d._repr_html_()

    # structure doesn't quite match the other cases
    assert r is not None

    if plot.backend.name == 'pylab':
        assert '<summary>DataIMG Plot</summary>' in r
        assert '<svg ' in r

    else:
        assert '<summary>DataIMG Data (' in r
        assert '<div class="dataname">X0</div>' in r
        assert '<svg ' not in r

    if (header is None) or (header == {}):
        assert '<summary>Metadata' not in r
    else:
        assert '<summary>Metadata (3)</summary>' in r


@requires_data
@requires_fits
@pytest.mark.parametrize('coord', ['logical', 'physical'])
def test_img_real(coord, make_data_path,
                  override_plot_backend, old_numpy_printing):
    """Use an image from a file (easy to set up)"""
    from sherpa.astro.io import read_image
    d = read_image(make_data_path('acisf07999_000N001_r0035_regevt3_srcimg.fits'))
    d.name = 'regevt3'
    d.set_coord(coord)

    r = d._repr_html_()

    assert r is not None

    if plot.backend.name == 'pylab':
        assert '<summary>DataIMG Plot</summary>' in r

    assert '<summary>Coordinates: physical (3)</summary>' in r
    assert '<summary>Coordinates: world (6)</summary>' in r
    assert '<summary>Metadata (6)</summary>' in r

    assert '<div class="dataval">[ 0.5  0.5]</div>' in r
    assert '<div class="dataval">[ 4096.5  4096.5]</div>' in r
    assert '<div class="dataval">[-0.000137  0.000137]</div>' in r
    assert '<div class="dataval">destreak - CAT3.0.2</div>' in r


# does this needs a decorator if region support is not available?
@requires_data
@requires_fits
@pytest.mark.parametrize('coord', ['logical', 'physical'])
def test_img_real_filtered(coord, make_data_path,
                           override_plot_backend, old_numpy_printing):
    """Filter the image.

    There is no significant check that the filtering has worked
    """
    from sherpa.astro.io import read_image
    d = read_image(make_data_path('acisf07999_000N001_r0035_regevt3_srcimg.fits'))
    d.name = 'regevt3'

    d.set_coord('logical')
    d.notice2d("circle(90, 190, 20)")
    d.set_coord(coord)

    r = d._repr_html_()

    assert r is not None

    if plot.backend.name == 'pylab':
        assert '<summary>DataIMG Plot</summary>' in r

    assert '<summary>Coordinates: physical (3)</summary>' in r
    assert '<summary>Coordinates: world (6)</summary>' in r
    assert '<summary>Metadata (6)</summary>' in r

    assert '<div class="dataval">[ 0.5  0.5]</div>' in r
    assert '<div class="dataval">[ 4096.5  4096.5]</div>' in r
    assert '<div class="dataval">[-0.000137  0.000137]</div>' in r
    assert '<div class="dataval">destreak - CAT3.0.2</div>' in r
