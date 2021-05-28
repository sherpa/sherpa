#
#  Copyright (C) 2017, 2018, 2021  Smithsonian Astrophysical Observatory
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

# Additional tests of filtering.

import logging

import pytest

import numpy as np

from sherpa.astro import ui
from sherpa.utils.testing import requires_data, requires_fits


def setup_model(make_data_path):
    """Set up a model that is reasonably close to the data.

    Returns the expected statistic values for various filters.
    """

    infile = make_data_path('q1127_src1_grp30.pi')

    ui.clean()
    ui.load_pha(infile)
    ui.subtract()

    ui.set_stat('chi2datavar')
    ui.set_source(ui.powlaw1d.pl)

    pl = ui.get_model_component('pl')
    pl.ampl = 5.28e-4
    pl.gamma = 1.04

    # These statistic values were created using CIAO 4.9 on a
    # Ubuntu machine. The quality=2 values are for high energies
    # (above ~ 10 keV or so), and so a filter of 0.5-8.0 keV should
    # give the same answer with or without ignore_bad.
    #
    return {'all': 2716.7086246284807,
            'bad': 2716.682482792285,
            '0.5-8.0': 1127.7165108405597}


@requires_fits
@requires_data
def test_filter_basic(make_data_path):
    """Test out issue 361 without ignore_bad calls.

    This should not trigger the bug behind issue 361.
    """

    stats = setup_model(make_data_path)

    s1 = ui.calc_stat()
    assert s1 == pytest.approx(stats['all'])

    ui.ignore_bad()
    s2 = ui.calc_stat()
    assert s2 == pytest.approx(stats['bad'])

    ui.notice(None, None)
    s3 = ui.calc_stat()
    assert s3 == pytest.approx(stats['all'])

    ui.notice(0.5, 8.0)
    s4 = ui.calc_stat()
    assert s4 == pytest.approx(stats['0.5-8.0'])


@requires_fits
@requires_data
def test_filter_notice_bad_361(make_data_path, caplog):
    """Test out issue 361: notice then ignore bad.

    Since the data is grouped, the ignore_bad call is expected to
    drop the filter expression, with a warning message.
    """

    stats = setup_model(make_data_path)

    # We don't care about these warnings, so I could just
    # store the number and check that we get an extra one
    # below, but use this as a canary to check when the
    # system changes.
    #
    assert len(caplog.records) == 5

    ui.notice(0.5, 8.0)
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.ignore_bad()

    assert len(caplog.records) == 6

    lname, lvl, msg = caplog.record_tuples[5]
    assert lname == 'sherpa.astro.data'
    assert lvl == logging.WARNING
    assert msg == 'filtering grouped data with quality ' + \
        'flags, previous filters deleted'

    s1 = ui.calc_stat()
    assert s1 == pytest.approx(stats['bad'])


@requires_fits
@requires_data
def test_filter_bad_notice_361(make_data_path):
    """Test out issue 361: ignore bad then notice.

    This is expected to fail with NumPy version 1.13 and should cause
    a DeprecationWarning with earlier versions (I do not know when
    the warning was added).
    """

    stats = setup_model(make_data_path)

    ui.ignore_bad()
    ui.notice(0.5, 8.0)
    s1 = ui.calc_stat()
    assert s1 == pytest.approx(stats['0.5-8.0'])


@requires_fits
@requires_data
def test_filter_bad_ungrouped(make_data_path, clean_astro_ui):
    """Check behavior when the data is ungrouped.

    This is a test of the current behavior, to check that
    values still hold. It may be necessary to change this
    test if we change the quality handling.
    """

    infile = make_data_path('q1127_src1_grp30.pi')
    ui.load_pha(infile)
    pha = ui.get_data()
    assert pha.quality_filter is None
    assert pha.mask is True

    assert ui.get_dep().shape == (439, )
    ui.ungroup()
    assert ui.get_dep().shape == (1024, )
    assert pha.quality_filter is None
    assert pha.mask is True

    ui.ignore_bad()
    assert ui.get_dep().shape == (1024, )
    assert pha.quality_filter is None

    expected = np.ones(1024, dtype=bool)
    expected[996:1025] = False
    assert pha.mask == pytest.approx(expected)

    # At this point we've changed the mask array so Sherpa thinks
    # we've applied a filter, so a notice is not going to change
    # anything. See issue #1169
    #
    ui.notice(0.5, 7)
    assert pha.mask == pytest.approx(expected)

    # We need to ignore to change the mask.
    #
    ui.ignore(None, 0.5)
    ui.ignore(7, None)
    expected[0:35] = False
    expected[479:1025] = False
    assert pha.mask == pytest.approx(expected)


@requires_fits
@requires_data
def test_filter_bad_grouped(make_data_path, clean_astro_ui):
    """Check behavior when the data is grouped.

    This is a test of the current behavior, to check that
    values still hold. It may be necessary to change this
    test if we change the quality handling.
    """

    infile = make_data_path('q1127_src1_grp30.pi')
    ui.load_pha(infile)
    pha = ui.get_data()
    assert pha.quality_filter is None
    assert pha.mask is True

    assert ui.get_dep().shape == (439, )
    assert pha.quality_filter is None
    assert pha.mask is True

    # The last group is marked as quality=2 and so calling
    # ignore_bad means we lose that group.
    #
    ui.ignore_bad()
    assert ui.get_dep().shape == (438, )
    assert pha.mask is True

    expected = np.ones(1024, dtype=bool)
    expected[996:1025] = False
    assert pha.quality_filter == pytest.approx(expected)

    # What happens when we filter the data? Unlike #1169
    # we do change the noticed range.
    #
    ui.notice(0.5, 7)
    assert pha.quality_filter == pytest.approx(expected)

    # The mask has been filtered to remove the bad channels
    # (this is grouped data)
    expected = np.ones(438, dtype=bool)
    expected[0:15] = False
    expected[410:438] = False
    assert pha.mask == pytest.approx(expected)

    expected = np.ones(996, dtype=bool)
    expected[0:34] = False
    expected[481:996] = False
    assert pha.get_mask() == pytest.approx(expected)
