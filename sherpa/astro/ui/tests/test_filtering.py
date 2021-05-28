#
#  Copyright (C) 2017, 2018  Smithsonian Astrophysical Observatory
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


# Check that a logging message is created: this is based on the
# code in https://stackoverflow.com/a/1049375 and
# https://stackoverflow.com/a/20553331
#
# Alternatively, we could require textfixtures. Can this be
# simplified once we drop Python versions < 3.4?
#
# Could this cause problems when running the whole test suite?
#
class MockLoggingHandler(logging.Handler):
    """Mock logging handler to check for expected logs.

    Messages are available from an instance's ``messages`` dict,
    in order, indexed by a lowercase log level string (e.g., 'debug',
    'info', etc.).
    """

    def __init__(self, *args, **kwargs):
        self.messages = {'debug': [], 'info': [], 'warning': [], 'error': [],
                         'critical': []}
        super(MockLoggingHandler, self).__init__(*args, **kwargs)

    def emit(self, record):
        "Store a message from ``record`` in the instance's ``messages`` dict."
        self.acquire()
        try:
            self.messages[record.levelname.lower()].append(record.getMessage())
        finally:
            self.release()

    def reset(self):
        self.acquire()
        try:
            for message_list in self.messages.values():
                message_list.clear()
        finally:
            self.release()


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
def test_filter_notice_bad_361(make_data_path):
    """Test out issue 361: notice then ignore bad.

    Since the data is grouped, the ignore_bad call is expected to
    drop the filter expression, with a warning message.
    """

    logger = logging.getLogger('sherpa')
    hdlr = MockLoggingHandler(level='WARNING')
    logger.addHandler(hdlr)

    stats = setup_model(make_data_path)

    ui.notice(0.5, 8.0)
    ui.ignore_bad()
    s1 = ui.calc_stat()
    assert s1 == pytest.approx(stats['bad'])

    msgs = hdlr.messages
    assert msgs['warning'] == ['filtering grouped data with quality ' +
                               'flags, previous filters deleted']
    for k in ['debug', 'info', 'error', 'critical']:
        assert msgs[k] == []


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
