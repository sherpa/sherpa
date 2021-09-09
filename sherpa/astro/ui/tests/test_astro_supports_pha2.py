#
#  Copyright (C) 2017, 2018, 2020, 2021
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
Check that we support PHA2 data sets (this is a format used by Chandra
to represent the multiple orders - and potentially - arms of a
grating dataset).

At present most of the specialized handling of PHA2 files is found in
the sherpa.astro.ui module, which is why the tests are placed here.

It is also a place to test the handling of multiple datasets, in
particular backgrounds, since even the PHA1 versions of the grating
data have two background components.
"""

import logging
import pytest

import numpy as np

from sherpa.utils.testing import requires_data, requires_fits
from sherpa.utils.err import IdentifierErr

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


one_to_12 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
hundred_to_11 = [99 + i for i in one_to_12]
xids = ["x{}".format(i) for i in one_to_12]


@requires_data
@requires_fits
@pytest.mark.parametrize("loader", [ui.load_pha, ui.load_data])
@pytest.mark.parametrize("id0,ids",
                         [(None, one_to_12),
                          (1, one_to_12),
                          (100, hundred_to_11),
                          ("x", xids)])
def test_load_pha2(loader, id0, ids, make_data_path, caplog, clean_astro_ui):
    """Basic test that a pha2 file can be read in."""

    basename = '3c120_pha2'

    orig_ids = ui.list_data_ids()
    assert orig_ids == []

    # The file is stored gzip-encoded
    infile = make_data_path(basename)
    if id0 is None:
        loader(infile)
    else:
        loader(id0, infile)

    pha_ids = ui.list_data_ids()
    assert len(pha_ids) == 12

    # list_data_ids doesn't guarantee an order
    # Do an explicit check, rather than via a set (testing
    # all at once) to make it easier to see what is missing
    # (if any)
    #
    for i in ids:
        assert i in pha_ids

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

    # Test Log messages
    msg_one = "systematic errors were not found in file '{}'".format(infile)

    # Editors can remove trailing spaces from lines, so split into
    # separate lines so the space after the file name is included.
    # Perhaps this space should be removed from the warning message?
    #
    msg_two = "statistical errors were found in file '{}' \n".format(infile) + \
              "but not used; to use them, re-read with use_errors=True"

    msg_three = "read background_up into a dataset from file {}".format(infile)
    msg_four = "read background_down into a dataset from file {}".format(infile)

    msg_five = "Multiple data sets have been input: " + \
               "{}-{}".format(ids[0], ids[11])

    assert caplog.record_tuples == [
        ('sherpa.astro.io', logging.WARNING, msg_one),
        ('sherpa.astro.io', logging.INFO, msg_two),
        ('sherpa.astro.io', logging.INFO, msg_three),
        ('sherpa.astro.io', logging.INFO, msg_four),
        ('sherpa.astro.ui.utils', logging.INFO, msg_five),
    ]


@requires_data
@requires_fits
def test_load_pha2_compare_meg_order1(make_data_path, clean_astro_ui):
    """Do we read in the MEG +/-1 orders?"""

    # The MEG -1 order is dataset 9
    # The MEG +1 order is dataset 10
    #
    pha2file = make_data_path('3c120_pha2')
    meg_p1file = make_data_path('3c120_meg_1.pha')
    meg_m1file = make_data_path('3c120_meg_-1.pha')

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


# See #397
#
@requires_data
@requires_fits
def test_list_bkg_ids(make_data_path, clean_astro_ui):
    """Does list_bkg_ids return a list"""

    basename = '3c120_pha2'

    infile = make_data_path(basename)
    ui.load_pha(infile)

    def validate(bids):
        """The documentation doesn't promise that the ids are
        listed in order, so do not assume so in the test.
        """
        assert len(bids) == 2
        assert set(bids) == set([1, 2])
        assert isinstance(bids, list)

    bids = ui.list_bkg_ids()
    validate(bids)

    for i in range(1, 13):
        bids = ui.list_bkg_ids(i)
        validate(bids)


# See #397
#
@requires_data
@requires_fits
def test_list_response_ids_pha2(make_data_path, clean_astro_ui):
    """Does list_response_ids return a list when input is pha2"""

    basename = '3c120_pha2'
    fakeid = 3

    infile = make_data_path(basename)
    ui.load_pha(infile)

    def validate(rids):
        "No response is read in with PHA2"
        assert rids == []

    def validate_err(i, excinfo):
        emsg = "background data set {} ".format(fakeid) + \
               "in PHA data set {} has not been set".format(i)
        assert str(excinfo.value) == emsg

    rids = ui.list_response_ids()
    validate(rids)

    with pytest.raises(IdentifierErr) as excinfo:
        ui.list_response_ids(bkg_id=fakeid)

    validate_err(1, excinfo)

    for i in range(1, 13):
        rids = ui.list_response_ids(i)
        validate(rids)

        for bkgid in [1, 2]:
            rids = ui.list_response_ids(i, bkg_id=bkgid)
            validate(rids)

        with pytest.raises(IdentifierErr) as excinfo:
            ui.list_response_ids(i, bkg_id=fakeid)

        validate_err(i, excinfo)


# See #397
#
@requires_data
@requires_fits
def test_list_response_ids_pha1(make_data_path, clean_astro_ui):
    """Does list_response_ids return a list when input is pha1"""

    basename = '3c120_heg_1.pha'
    fakeid = 3

    infile = make_data_path(basename)
    ui.load_pha(infile)

    def validate(rids):
        assert len(rids) == 1
        assert 1 in rids
        assert isinstance(rids, list)

    def validate_err(i, excinfo):
        emsg = "background data set {} ".format(fakeid) + \
               "in PHA data set {} has not been set".format(i)
        assert str(excinfo.value) == emsg

    rids = ui.list_response_ids()
    validate(rids)

    for bid in [1, 2]:
        rids = ui.list_response_ids(bkg_id=bid)
        validate(rids)

    with pytest.raises(IdentifierErr) as excinfo:
        ui.list_response_ids(bkg_id=fakeid)

    validate_err(1, excinfo)

    ui.load_pha('heg1', infile)
    rids = ui.list_response_ids('heg1')
    validate(rids)

    for bid in [1, 2]:
        rids = ui.list_response_ids('heg1', bkg_id=bid)
        validate(rids)

    with pytest.raises(IdentifierErr) as excinfo:
        ui.list_response_ids('heg1', bkg_id=fakeid)

    validate_err('heg1', excinfo)


@requires_data
@requires_fits
def test_746(make_data_path, clean_astro_ui):
    """Test https://github.com/sherpa/sherpa/issues/746

    Something in #444 (reverted in #759) caused:

      - the fit to fail (niter=2)
      - the line amplitude not to change significantly
      - the statistic reported by the fit to be different to
        that returned by calc_stat

    Something with how the cache code handles analysis=wave
    appears to be the problem. This test takes the line data
    from 746 and adds it into the existing PHA2 file we have
    (3c120) to replicate the problem. Fortunately this
    wavelength range contains essentially no counts, so we
    can just add in the counts to make a fake line and check
    we can fit it.
    """

    ui.load_pha(make_data_path('3c120_pha2.gz'))
    ui.load_arf(10, make_data_path('3c120_meg_1.arf.gz'))
    ui.load_rmf(10, make_data_path('3c120_meg_1.rmf.gz'))

    # Add in the line
    d10 = ui.get_data(10)
    idx = np.arange(4068, 4075, dtype=int)
    d10.counts[idx] = [1, 1, 1, 2, 3, 1, 1]

    # group the data
    ui.group_width(id=10, num=2)

    ui.set_analysis('wave')
    ui.notice(21.4, 21.7)

    # internal check that getting the expected data
    expected = np.zeros(30)
    expected[[8, 9, 10, 11]] = [2, 3, 4, 1]
    expected[25] = 1  # this count is from the 3c120 data
    d = d10.get_dep(filter=True)
    assert d == pytest.approx(expected)

    line = ui.create_model_component('gauss1d', 'name')

    # treat the line as a delta function
    line.fwhm.val = 0.00001
    line.fwhm.frozen = True
    line.pos = 21.6
    line.pos.frozen = True

    line.ampl = 2

    ui.set_source(10, line)

    # the original fit used levmar, so use that here
    # (it looks like it also fails with simplex)
    ui.set_method('levmar')
    ui.set_stat('cstat')

    sinit = ui.calc_stat(10)
    ui.fit(10)
    fr = ui.get_fit_results()
    sfinal = ui.calc_stat(10)

    # Did the
    #   - fit improve
    #   - take more than two iterations
    #   - report the same statistic values as calc_stat
    #   - change the fit parameter
    #
    assert sfinal < sinit

    assert fr.succeeded
    assert fr.nfev > 2

    assert fr.istatval == sinit
    assert fr.statval == sfinal

    assert line.ampl.val > 100
    assert len(fr.parvals) == 1
    assert fr.parvals[0] == line.ampl.val

    # some simple checks to throw in because we can
    assert fr.parnames == ('name.ampl', )
    assert fr.datasets == (10,)
    assert fr.statname == 'cstat'
    assert fr.methodname == 'levmar'
    assert fr.itermethodname == 'none'
    assert fr.numpoints == 30
    assert fr.dof == 29

    # Now add in some "absolute" checks to act as regression tests.
    # If these fail then it doesn't necessarily mean something bad
    # has happened.
    #
    assert fr.nfev == 15
    assert sinit == pytest.approx(82.72457294394245)
    assert sfinal == pytest.approx(15.395474325257306)
    assert line.ampl.val == pytest.approx(113.97801045821409)
