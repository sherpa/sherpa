#
#  Copyright (C) 2024
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


from setuptools.command.build import build as _build
from setuptools.command.build_ext import build_ext as _build_ext


class build(_build):
    """Ensure the configuration details and code are available/run.

    Needed for the build_ext code, to copy over stk/group libraries
    for editable installs (although it seems like it would make sense
    to include these config options somewhere now we don't have the
    setup.cfg command aliases).

    """

    sub_commands = [("sherpa_config", None),
                    ("xspec_config", None)] + _build.sub_commands


class build_ext(_build_ext):
    """This is a hack to copy over the stk/group libraries for
    editable installs. There must be a better way to do this.

    """

    def run(self):

        # Copy across the stk and group libraries to somewhere
        # "useful" for editable installs. This is really not a
        # great setup.
        #
        if self.editable_mode:
            sherpa_config = self.get_finalized_command('sherpa_config', True)
            if not sherpa_config.disable_group:
                self.copy_file(sherpa_config.group_location,
                               "group.so")
                self.announce("install group extension locally")

            if not sherpa_config.disable_stk:
                self.copy_file(sherpa_config.stk_location,
                               "stk.so")
                self.announce("install stk extension locally")

        super().run()
