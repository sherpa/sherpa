#
#  Copyright (C) 2019, 2020, 2021  Smithsonian Astrophysical Observatory
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

"""This is labelled ui but really tests astro/ui

"""

import pytest

from sherpa.utils.testing import requires_data, \
    requires_xspec, requires_fits

from sherpa.astro import ui


@pytest.mark.xfail
@requires_xspec
@requires_fits
@requires_data
def test_plot_pvalue(make_data_path, clean_astro_ui, hide_logging):
    """Check plot_pvalue with PHA data."""

    fname = make_data_path("qso.pi")
    ui.load_pha(fname)

    ui.set_stat('cstat')
    ui.set_method("neldermead")

    ui.group_counts(10)
    ui.notice(0.3, 8)

    ui.set_model("xsphabs.abs1*(xspowerlaw.p1 +gauss1d.g1)")

    # move the fit close to the best fit to save a small amount
    # of time.
    abs1.nh = 0.05
    p1.phoindex = 1.28
    p1.norm = 2e-4
    g1.ampl = 1.8e-5

    g1.pos = 3.
    ui.freeze(g1.pos)
    g1.fwhm = 0.1
    ui.freeze(g1.fwhm)

    # Could we reduce the number of bins to save evaluation time?
    # We do want a non-default num value when checking the shapes
    # of the output attributes.
    #
    ui.fit()
    ui.plot_pvalue(p1, p1 + g1, num=100, bins=20)

    tmp = ui.get_pvalue_results()

    assert tmp.null == pytest.approx(210.34566845619273)
    assert tmp.alt == pytest.approx(207.66618095925094)
    assert tmp.lr == pytest.approx(2.679487496941789)

    # Have we returned the correct info?
    #
    # Is it worth checking the stored data (aka how randomised is this
    # output)?
    #
    assert tmp.samples.shape == (100, 2)
    assert tmp.stats.shape == (100, 2)
    assert tmp.ratios.shape == (100, )

    # Check the plot
    #
    tmp = ui.get_pvalue_plot()

    assert tmp.lr == pytest.approx(2.679487496941789)

    assert tmp.xlabel == 'Likelihood Ratio'
    assert tmp.ylabel == 'Frequency'
    assert tmp.title == 'Likelihood Ratio Distribution'

    # It would be nice to check the values here
    #
    assert tmp.ratios.shape == (100, )
    assert tmp.xlo.shape == (21, )
    assert tmp.xhi.shape == (21, )
    assert tmp.y.shape == (21, )
