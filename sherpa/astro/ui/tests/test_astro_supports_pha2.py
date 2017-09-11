#
#  Copyright (C) 2016  Smithsonian Astrophysical Observatory
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
Check that we support PHA2 data sets (this is a format used by Chandra
to represent the multiple orders - and potentially - arms of a
grating dataset).

At present most of the specialized handling of PHA2 files is found in
the sherpa.astro.ui module, which is why the tests are placed here.

"""

import pytest

from sherpa.utils import requires_data, requires_fits

from sherpa.astro import ui
from sherpa.astro.data import DataPHA


def validate_pha(d, bkg=True):
    """Check the data is a PHA with expected structure.

    Parameters
    ----------
    d
        The data object to check
    bkg : bool, optional
        If True then the data set is expected to have two
        background components.
    """

    exptime = 77716.294300038993
    assert isinstance(d, DataPHA)

    # There are counts, and we have the expected number
    # of channels (the channels are stored as floats
    # even though they are integers)
    #
    assert d.channel.min() == pytest.approx(1.0)
    assert d.channel.max() == pytest.approx(8192.0)
    assert d.channel.sum() == pytest.approx(8193.0 * 8192.0 / 2.0)

    assert len(d.bin_lo) == 8192
    assert len(d.bin_hi) == 8192

    # check the bins are increasing (so bin_lo[i] < bin_hi[i])
    # but in descending order
    assert (d.bin_lo < d.bin_hi).all()
    assert (d.bin_lo[:-1] > d.bin_lo[1:]).all()

    # are they consecutive?
    assert (d.bin_lo[:-1] == d.bin_hi[1:]).all()

    assert d.counts.sum() > 0
    assert (d.counts >= 0).all()

    assert d.exposure == pytest.approx(exptime)

    assert d.staterror is None
    assert d.syserror is None

    if bkg:
        assert d.backscal == pytest.approx(1.0)
        assert d.areascal == pytest.approx(1.0)
    else:
        # The BACKSKUP and BACKSKDN keywords are set to this
        assert d.backscal == pytest.approx(4.01882840)

        # Why is areascal set to None?
        assert d.areascal is None

    assert not(d.grouped)
    assert not(d.subtracted)

    assert d.units == 'channel'
    assert d.rate
    assert d.plot_fac == 0

    assert d.response_ids == []
    if bkg:
        assert d.background_ids == [1, 2]
    else:
        assert d.background_ids == []


@requires_data
@requires_fits
def test_load_pha2(make_data_path):
    """Basic test that a pha2 file can be read in."""

    basename = '3c120_pha2'

    ui.clean()
    orig_ids = ui.list_data_ids()
    assert orig_ids == []

    # The file is stored gzip-encoded
    infile = make_data_path(basename)
    ui.load_pha(infile)

    pha_ids = ui.list_data_ids()
    assert len(pha_ids) == 12

    # list_data_ids doesn't guarantee an order
    # Do an explicit check, rather than via a set (testing
    # all at once) to make it easier to see what is missing
    # (if any)
    #
    for i in range(1, 13):
        assert i in pha_ids

    for i in range(1, 13):
        d = ui.get_data(i)
        validate_pha(d, bkg=True)

        # There is no indication of what "part" this data set
        # represents in the file name
        #
        assert d.name == infile

        b = ui.get_bkg(i, bkg_id=1)
        validate_pha(b, bkg=False)
        assert b.name == infile

        b = ui.get_bkg(i, bkg_id=2)
        validate_pha(b, bkg=False)
        assert b.name == infile

    ui.clean()


# TODO: how to test what messages are logged when load_pha is
#       called with a pha2 file?


@requires_data
@requires_fits
def test_load_pha2_compare_meg_order1(make_data_path):
    """Do we read in the MEG +/-1 orders?"""

    # The MEG -1 order is dataset 9
    # The MEG +1 order is dataset 10
    #
    pha2file = make_data_path('3c120_pha2')
    meg_p1file = make_data_path('3c120_meg_1.pha')
    meg_m1file = make_data_path('3c120_meg_-1.pha')

    ui.clean()
    ui.load_pha('meg_p1', meg_p1file)
    ui.load_pha('meg_m1', meg_m1file)

    orig_ids = set(ui.list_data_ids())
    assert 'meg_p1' in orig_ids
    assert 'meg_m1' in orig_ids

    ui.load_pha(pha2file)

    for n, lbl in zip([9, 10], ["-1", "1"]):
        h = '3c120_meg_{}'.format(lbl)
        ui.load_arf(n, make_data_path(h + '.arf'))
        ui.load_rmf(n, make_data_path(h + '.rmf'))

    # check that loading the pha2 file doesn't overwrite existing
    # data
    new_ids = set(ui.list_data_ids())

    for i in range(1, 13):
        orig_ids.add(i)

    assert orig_ids == new_ids

    # Check that the same model gives the same statistic
    # value; this should check that the data and response are
    # read in, that grouping and filtering work, and that
    # model evaluation is the same, without having to
    # check these steps individually.
    #
    # The model is not meant to be physically meaningful,
    # just one that reasonably represents the data and
    # can be evaluated without requiring XSPEC.
    #
    pmdl = ui.create_model_component('powlaw1d', 'pmdl')
    pmdl.gamma = 0.318
    pmdl.ampl = 2.52e-3

    ncts = 20
    for i in [9, 10, "meg_m1", "meg_p1"]:
        ui.set_analysis(i, 'wave')
        ui.group_counts(i, ncts)
        ui.notice_id(i, 2, 12)
        ui.set_source(i, pmdl)

    ui.set_stat('chi2datavar')
    s9 = ui.calc_stat(9)
    s10 = ui.calc_stat(10)

    sm1 = ui.calc_stat('meg_m1')
    sp1 = ui.calc_stat('meg_p1')

    # Since these should be the same, we use an equality test
    # rather than approximation. At least until it becomes
    # a problem.
    #
    assert s9 == sm1
    assert s10 == sp1

    # The values were calculated using CIAO 4.9, Linux64, with
    # Python 3.5.
    #
    assert s9 == pytest.approx(1005.4378559390879)
    assert s10 == pytest.approx(1119.980439489647)
