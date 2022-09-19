#
#  Copyright (C) 2017, 2018, 2021, 2022
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

# Additional tests of filtering.

import logging

import pytest

import numpy as np

from sherpa.astro.data import DataIMG, DataPHA
from sherpa.astro import ui
from sherpa.astro.ui.utils import Session as AstroSession
from sherpa.data import Data1DInt, Data2D
from sherpa.ui.utils import Session
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

    assert len(caplog.records) == 8

    lname, lvl, msg = caplog.record_tuples[5]
    assert lname == 'sherpa.ui.utils'
    assert lvl == logging.INFO
    assert msg == 'dataset 1: 0.0073:14.9504 -> 0.4964:8.0592 Energy (keV)'

    lname, lvl, msg = caplog.record_tuples[6]
    assert lname == 'sherpa.astro.data'
    assert lvl == logging.WARNING
    assert msg == 'filtering grouped data with quality ' + \
        'flags, previous filters deleted'

    lname, lvl, msg = caplog.record_tuples[7]
    assert lname == 'sherpa.ui.utils'
    assert lvl == logging.INFO
    assert msg == 'dataset 1: 0.4964:8.0592 -> 0.0073:14.9504 Energy (keV)'

    s1 = ui.calc_stat()
    assert s1 == pytest.approx(stats['bad'])


@requires_fits
@requires_data
def test_filter_bad_notice_361(make_data_path):
    """Test out issue 361: ignore bad then notice.

    There have been changes in NumPy that this originally was written
    to catch (version 1.13). This has been left in, but it is
    currently a regression test as the handling of bad values means
    that the notice call fails when trying to write out the filter,
    hence the try/except handling there.

    """

    stats = setup_model(make_data_path)

    ui.ignore_bad()

    # With support for #1558 - writing out the filter as we apply them
    # - we now get this call failing, so catch this failure to allow
    # the test to continue. This should be fixed, but it's a
    # long-standing problem.
    #
    try:
        ui.notice(0.5, 8.0)
    except IndexError:
        pass

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
def test_filter_bad_grouped(make_data_path, clean_astro_ui, caplog):
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
    assert len(caplog.records) == 5
    ui.ignore_bad()
    assert len(caplog.records) == 6

    assert ui.get_dep().shape == (438, )
    assert pha.mask is True

    expected = np.ones(1024, dtype=bool)
    expected[996:1025] = False
    assert pha.quality_filter == pytest.approx(expected)

    # Do we really want this (added in #1562)? For now assume so.
    #
    lname, lvl, msg = caplog.record_tuples[5]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 0.0073:14.9504 Energy (keV) (unchanged)"

    # What happens when we filter the data? Unlike #1169
    # we do change the noticed range.
    #
    # This should create a filter message, but thanks to ignore_bad we
    # get an IndexError, which has been caught and that causes the
    # filter message to not be displayed.
    #
    ui.notice(0.5, 7)
    assert len(caplog.records) == 6

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


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("expr,result",
                         [("4:13", "5:13"),
                          (":11", "-100:10"),
                          ("5:", "5:100"),
                          (":", None),
                          (":-100", "-100"),
                          ("100:", "100"),
                          ("2", "2"),
                          ("11", ""),
                          (":-150", ""),
                          ("200:", "")
                          ])
def test_notice_string_data1d(session, expr, result, caplog):
    """Check we can call notice with a string argument: Data1D"""

    x = np.asarray([-100, 1, 2, 3, 5, 10, 12, 13, 100])
    y = np.ones_like(x)

    s = session()
    s.load_arrays(1, x, y)

    assert len(caplog.record_tuples) == 0
    s.notice(expr)
    assert len(caplog.record_tuples) == 1

    expected = "-100:100" if result is None else result
    assert s.get_data().get_filter(format='%d') == expected

    expected = "dataset 1: -100:100 "
    if result is None:
        expected += "x (unchanged)"
    elif result == "":
        expected += "x -> no data"
    else:
        expected += f"-> {result} x"

    loc, lvl, msg = caplog.record_tuples[0]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == expected


@pytest.mark.parametrize("session", [Session, AstroSession])
@pytest.mark.parametrize("expr,result",
                         [("4:13", "5:13"),
                          (":11", "-100:11"),
                          ("5:", "5:100"),
                          (":", None),
                          (":-100", ""),
                          ("100:", ""),
                          ("2", ""),
                          ("11", ""),
                          (":-150", ""),
                          ("200:", "")
                          ])
def test_notice_string_data1dint(session, expr, result, caplog):
    """Check we can call notice with a string argument: Data1DInt"""

    xlo = np.asarray([-100, 1, 2, 3, 5, 10, 12, 13, 99])
    xhi = np.asarray([-99, 2, 3, 4, 10, 11, 13, 90, 100])
    y = np.ones_like(xlo)

    s = session()
    s.load_arrays(1, xlo, xhi, y, Data1DInt)

    assert len(caplog.record_tuples) == 0
    s.notice(expr)
    assert len(caplog.record_tuples) == 1

    expected = "-100:100" if result is None else result
    assert s.get_data().get_filter(format='%d') == expected

    expected = "dataset 1: -100:100 "
    if result is None:
        expected += "x (unchanged)"
    elif result == "":
        expected += "x -> no data"
    else:
        expected += f"-> {result} x"

    loc, lvl, msg = caplog.record_tuples[0]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == expected


@pytest.mark.parametrize("expr,result",
                         [("4:13", "4:13"),
                          (":11", "1:11"),
                          ("5:", "5:19"),
                          (":", None),
                          (":-100", ""),
                          ("100:", ""),
                          ("2", "2"),
                          ("11", "11"),
                          (":-150", ""),
                          ("200:", "")
                          ])
def test_notice_string_datapha(expr, result, caplog):
    """Check we can call notice with a string argument: DataPHA"""

    x = np.arange(1, 20)
    y = np.ones_like(x)

    s = AstroSession()
    s.load_arrays(1, x, y, DataPHA)

    assert len(caplog.record_tuples) == 0
    s.notice(expr)
    assert len(caplog.record_tuples) == 1

    expected = "1:19" if result is None else result
    assert s.get_data().get_filter(format='%d') == expected

    expected = "dataset 1: 1:19 "
    if result is None:
        expected += "Channel (unchanged)"
    elif result == "":
        expected += "Channel -> no data"
    else:
        expected += f"-> {result} Channel"

    loc, lvl, msg = caplog.record_tuples[0]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == expected


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_notice_reporting_data1d(session, caplog):
    """Unit-style test of logging of notice/ignore: Data1D"""

    x = np.asarray([-100, 1, 2, 3, 5, 10, 12, 13, 100])
    y = np.ones_like(x)

    s = session()
    s.load_arrays(1, x, y)

    assert len(caplog.record_tuples) == 0
    s.notice()
    assert len(caplog.record_tuples) == 1

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: -100:100 x (unchanged)"

    s.notice(lo=2.5)
    assert len(caplog.record_tuples) == 2

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: -100:100 -> 3:100 x"

    s.ignore(6, 15)
    assert len(caplog.record_tuples) == 3

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 3:100 -> 3:5,100 x"

    s.notice_id(1, 12.5, 15)
    assert len(caplog.record_tuples) == 4

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 3:5,100 -> 3:5,13:100 x"

    s.ignore_id(1)
    assert len(caplog.record_tuples) == 5

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 3:5,13:100 x -> no data"

    s.ignore()
    assert len(caplog.record_tuples) == 6

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: no data (unchanged)"


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_notice_reporting_data1dint(session, caplog):
    """Unit-style test of logging of notice/ignore: Data1DInt"""

    xlo = np.asarray([-100, 1, 2, 3, 5, 10, 12, 13, 99])
    xhi = np.asarray([-99, 2, 3, 4, 10, 11, 13, 90, 100])
    y = np.ones_like(xlo)

    s = session()
    s.load_arrays(1, xlo, xhi, y, Data1DInt)

    assert len(caplog.record_tuples) == 0
    s.notice(lo=2.5)
    assert len(caplog.record_tuples) == 1

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: -100:100 -> 2:100 x"

    s.ignore(6, 15)
    assert len(caplog.record_tuples) == 2

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 2:100 -> 2:4,99:100 x"

    s.notice_id(1, 12.5, 15)
    assert len(caplog.record_tuples) == 3

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 2:4,99:100 -> 2:4,12:100 x"

    s.ignore_id(1)
    assert len(caplog.record_tuples) == 4

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 2:4,12:100 x -> no data"


def test_notice_reporting_datapha(caplog):
    """Unit-style test of logging of notice/ignore: DataPHA"""

    x = np.arange(1, 20)
    y = np.ones_like(x)

    s = AstroSession()
    s.load_arrays(1, x, y, DataPHA)

    assert len(caplog.record_tuples) == 0
    s.notice(lo=3)
    assert len(caplog.record_tuples) == 1

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 1:19 -> 3:19 Channel"

    s.ignore(6, 15)
    assert len(caplog.record_tuples) == 2

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 3:19 -> 3:5,16:19 Channel"

    s.notice_id(1, 13, 15)
    assert len(caplog.record_tuples) == 3

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 3:5,16:19 -> 3:5,13:19 Channel"

    s.ignore_id(1)
    assert len(caplog.record_tuples) == 4

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: 3:5,13:19 -> No noticed bins Channel"

    s.ignore()
    assert len(caplog.record_tuples) == 5

    loc, lvl, msg = caplog.record_tuples[-1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: No noticed bins Channel (unchanged)"


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_notice_reporting_data2d(session, caplog):
    """Unit-style test of logging of notice/ignore: Data2D

    At the moment it's not clear what the notice/ignore call is
    meant to do, as the ui layer doesn't really support the Data2D
    case.
    """

    x1, x0 = np.mgrid[10:20, 1:7]
    shape = x0.shape
    x0 = x0.flatten()
    x1 = x1.flatten()

    y = np.ones_like(x0)

    s = session()
    s.load_arrays(1, x0, x1, y, shape, Data2D)

    assert len(caplog.record_tuples) == 0
    assert s.get_filter() == ""

    # At the moment it's not obvious what the arguments mean, so just
    # act as a regression test. In particuar, we get no caplog output.
    #
    s.notice(lo=13)
    assert len(caplog.record_tuples) == 0
    assert s.get_filter() == ""

    s.ignore(6, 15)
    assert len(caplog.record_tuples) == 0
    assert s.get_filter() == ""

    s.notice_id(1, 13, 15)
    assert len(caplog.record_tuples) == 0
    assert s.get_filter() == ""

    s.ignore_id(1)
    assert len(caplog.record_tuples) == 0
    assert s.get_filter() == ""


def test_notice2d_reporting(caplog):
    """Check handling of notice2d/ignore2d/... reports"""

    # Remember, although we give a grid range here,
    # the filtering below is in logical coordinates.
    #
    x1, x0 = np.mgrid[10:20, 1:7]
    shape = x0.shape
    x0 = x0.flatten()
    x1 = x1.flatten()

    y = np.ones_like(x0)

    s = AstroSession()
    s.load_arrays(1, x0, x1, y, shape, DataIMG)

    assert len(caplog.record_tuples) == 0
    s.notice2d("circle(3, 15, 5)")
    assert len(caplog.record_tuples) == 1

    loc, lvl, msg = caplog.record_tuples[0]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: Field() -> Circle(3,15,5)"

    s.notice2d_id(1, "rect(4, 12, 7, 14)")
    assert len(caplog.record_tuples) == 2

    loc, lvl, msg = caplog.record_tuples[1]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: Circle(3,15,5) -> Circle(3,15,5)|Rectangle(4,12,7,14)"

    # I was trying to remove all the filters, but this doesn't seem to
    # do it. Maybe the region logic is not able to recognize that this
    # is now the null filter.
    #
    s.ignore2d("rect(-5, -5, 25, 25)")
    assert len(caplog.record_tuples) == 3

    loc, lvl, msg = caplog.record_tuples[2]
    assert loc == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "dataset 1: Circle(3,15,5)|Rectangle(4,12,7,14) -> Circle(3,15,5)&!Rectangle(-5,-5,25,25)|Rectangle(4,12,7,14)&!Rectangle(-5,-5,25,25)"
