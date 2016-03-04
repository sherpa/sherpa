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
import logging
import sys


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
