from __future__ import absolute_import
#
# Copyright (C) 2015  Smithsonian Astrophysical Observatory
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
from StringIO import StringIO
import logging

from sherpa.utils import SherpaTestCase
from sherpa.utils.logging import config_logger


class LoggingUnitTest(SherpaTestCase):
    INFO_MESSAGE = "an info"
    WARNING_MESSAGE = "a warning"
    TEMPLATE = "logging_test: {}\n"
    LOGGER_NAME = "logging_test"
    ALT_TEMPLATE = "template"

    def _test_log(self, expected):
        self.assertEqual(expected, self.stream.getvalue())

    def setUp(self):
        self.stream = StringIO()

    def tearDown(self):
        self.stream.close()

    def test_warning(self):
        logger = config_logger(self.LOGGER_NAME, logging.WARNING, self.stream)

        # Test that info messages are not printed
        logger.info(self.INFO_MESSAGE)
        self.assertEqual("", self.stream.getvalue())

        # Test that warning messages are printed
        logger.warning(self.WARNING_MESSAGE)
        self.assertEqual(
            self.TEMPLATE.format(self.WARNING_MESSAGE), self.stream.getvalue())

    def test_info(self):
        # Now test that INFOs get through when setting the level to INFO
        logger = config_logger(self.LOGGER_NAME, logging.INFO, self.stream)

        # Test that warning messages are printed
        logger.warning(self.WARNING_MESSAGE)
        log = self.TEMPLATE.format(self.WARNING_MESSAGE)
        self.assertEqual(log, self.stream.getvalue())

        # Test that info messages are printed too
        logger.info(self.INFO_MESSAGE)
        log = "{}{}".format(log, self.TEMPLATE.format(self.INFO_MESSAGE))
        self.assertEqual(log, self.stream.getvalue())

    def test_template(self):
        logger = config_logger(
            self.LOGGER_NAME, template=self.ALT_TEMPLATE, stream=self.stream)
        logger.warning(self.WARNING_MESSAGE)
        self.assertEqual(self.ALT_TEMPLATE + "\n", self.stream.getvalue())
