# 
#  Copyright (C) 2014  Smithsonian Astrophysical Observatory
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


from build import build
from clean import clean
from install import install
from sdist import sdist
from sherpa_config import sherpa_config
from xspec_config import xspec_config

commands = {
             'build': build,
             'clean' : clean,
             'install' : install,
             'sdist' : sdist,
             'sherpa_config' : sherpa_config,
             'xspec_config' : xspec_config,
            }
