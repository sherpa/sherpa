from __future__ import print_function
#
#  Copyright (C) 2015  Smithsonian Astrophysical Observatory
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

import shutil
import os

try:
    from numpy.distutils.command.develop import develop as _develop

    class develop(_develop):

        def run(self):
            _develop.run(self)
            sherpa_config = self.get_finalized_command('sherpa_config', True)
            self.announce("install stk and group extensions locally")
            if not sherpa_config.disable_stk:
                shutil.copyfile(sherpa_config.stk_location, os.path.join(os.getcwd(), 'stk.so'))
            if not sherpa_config.disable_group:
                shutil.copyfile(sherpa_config.group_location, os.path.join(os.getcwd(), 'group.so'))

except ImportError:
    from distutils.core import Command

    class develop(Command):

        user_options = []

        def run(self):
            print("develop command is not available without setuptools")

        def initialize_options(self):
            pass

        def finalize_options(self):
            pass
