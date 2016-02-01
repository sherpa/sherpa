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


from distutils.command.sdist import sdist as _sdist
from numpy.distutils.misc_util import get_data_files
from .deps import clean_deps

class sdist(_sdist):

    def add_defaults(self):
        _sdist.add_defaults(self)

        dist = self.distribution

        if dist.has_data_files():
            for data in dist.data_files:
                self.filelist.extend(get_data_files(data))

        if dist.has_headers():
            headers = []
            for h in dist.headers:
                if isinstance(h,str): headers.append(h)
                else: headers.append(h[1])
            self.filelist.extend(headers)

        return

    def run(self):
        clean_deps()
        self.get_finalized_command('sherpa_config', True).build_configure()
        _sdist.run(self)
