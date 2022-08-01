#
#  Copyright (C) 2014, 2016, 2022
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


# It looks like setuptools sdist is currently incomplete, see
# https://github.com/numpy/numpy/pull/7131
# but the replacement test suggests using distutils,
# which has an add_defaults method.
#
# from setuptools.command.sdist import sdist as _sdist
# from numpy.distutils.command.sdist import sdist as _sdist
from distutils.command.sdist import sdist as _sdist
from .deps import clean_deps


class sdist(_sdist):

    def run(self):
        clean_deps()
        # There is no build_configure step for xspec_config
        self.get_finalized_command('sherpa_config', True).build_configure()
        super().run()
