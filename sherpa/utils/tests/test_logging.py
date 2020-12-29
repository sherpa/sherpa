#
#  Copyright (C) 2020  Smithsonian Astrophysical Observatory
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
import pytest
import logging
import numpy as np
import sherpa
from sherpa.utils.logging import config_logger, SherpaVerbosity
from sherpa.astro import ui


def test_logging_verbosity_contextmanager(caplog):
    # Previous tests might have changed this to a non-default
    # level already. So, we first make sure the root logger
    # is set right.
    sherpalogger = logging.getLogger('sherpa')
    sherpalogger.setLevel('WARNING')

    logger = logging.getLogger('sherpa.some_module')
    logger.warning('1: should be seen')
    assert len(caplog.records) == 1

    with SherpaVerbosity('ERROR'):
        logger.warning('2: Should not be seen')
    assert len(caplog.records) == 1

    logger.warning('3: should be seen')
    assert len(caplog.records) == 2
