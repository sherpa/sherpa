#
# Copyright (C) 2015-2020  Smithsonian Astrophysical Observatory
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
import logging
import sys
import contextlib


def config_logger(name, level=logging.WARNING, stream=sys.stdout, template=None):
    if template is None:
        template = '{}: %(message)s'.format(name)

    logger = logging.getLogger(name)
    for handler in logger.handlers:
        logger.removeHandler(handler)

    fmt = logging.Formatter(template, None)
    handler = logging.StreamHandler(stream)
    handler.setFormatter(fmt)

    logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False

    return logger


class SherpaVerbosity(contextlib.AbstractContextManager):
    '''Set the output logging level for sherpa as a context.

    This changes the logging level globally for all modules in sherpa.

    Parameters
    ----------
    loglevel : string or int
        New level for logging. Allowed strings are
        ``DEBUG``, ``INFO``, ``WARNING``, ``ERROR``, and ``CRITICAL``
    '''
    def __init__(self, level):
        self.level = level
        self.sherpalog = logging.getLogger('sherpa')

    def __enter__(self):
        self.old = self.sherpalog.level
        self.sherpalog.setLevel(self.level)

    def __exit__(self, *args):
        self.sherpalog.setLevel(self.old)
