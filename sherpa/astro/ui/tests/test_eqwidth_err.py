#
#  Copyright (C) 2018  Smithsonian Astrophysical Observatory
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

from sherpa.utils import requires_data, requires_fits, requires_xspec
import numpy

from pytest import approx

from sherpa.astro.ui import *


@requires_data
@requires_fits
@requires_xspec
def test_eqwith_err(make_data_path):
    set_method('neldermead')
    set_stat('cstat')
    load_data(make_data_path('12845.pi'))
    notice(0.5,7)

    set_model("xsphabs.gal*xszphabs.zabs*(powlaw1d.p1+xszgauss.g1)")
    set_par(gal.nh,0.08)
    freeze(gal)

    set_par(zabs.redshift, 0.518)
    set_par(g1.redshift, 0.518)
    set_par(g1.Sigma, 0.01)
    freeze(g1.Sigma)
    set_par(g1.LineE,min=6.0,max=7.0)
    fit()
    numpy.random.seed(12345)
    result = eqwidth(p1,p1+g1, error=True, niter=100)
    assert result[0] == approx(0.16443033244310976, rel=1e-3)
    assert result[1] == approx(0.09205564216156815, rel=1e-3)
    assert result[2] == approx(0.23933118287470895, rel=1e-3)
