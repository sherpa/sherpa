#
#  Copyright (C) 2014-2018, 2020-2025
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

from setuptools import Command
import re

from .extensions import build_ext

# What versions of XSPEC do we support? I am not sure what the
# naming of the XSPEC components are, but let's stick with
# major, minor, and micro. We drop the patch level - e.g.
# "c" in "12.12.0c" as that is not helpful to track here.
#
SUPPORTED_VERSIONS = [(12, 13, 0), (12, 13, 1),
                      (12, 14, 0), (12, 14, 1),
                      (12, 15, 0), (12, 15, 1)]


# We could use packaging.versions.Version here, but for our needs we
# can get away with a tuple of integers. That is, we do not need the
# full support for PEP-440.
#
MIN_VERSION = min(SUPPORTED_VERSIONS)
MAX_VERSION = max(SUPPORTED_VERSIONS)


def clean(xs):
    "Remove all '' entries from xs, returning the new list."
    return [x for x in xs if x != '']


def get_version(version):
    """Strip out any XSPEC patch level.

    So '12.12.0c' gets converted to '12.12.0', and then to (12, 12,
    0). This is helpful as then it makes version comparison easier, as
    we can rely on the standard tuple ordering.

    Parameters
    ----------
    version : str
        The XSPEC version string, of the form "12.12.0c", so it can
        include the XSPEC patch level.

    Returns
    -------
    (major, minor, micro) : tuple of int
        The XSPEC patchlevel is ignored.

    """

    # XSPEC versions do not match PEP 440, so strip out the trailing
    # text (which indicates the XSPEC patch level).
    #
    match = re.search(r'^(\d+)\.(\d+)\.(\d+)', version)
    if match is None:
        raise ValueError(f"Invalid XSPEC version string: {version}")

    return (int(match[1]), int(match[2]), int(match[3]))


class xspec_config(Command):
    description = "Configure XSPEC Models external module (optional) "
    user_options = [
                    ('with_xspec', None, "Whether sherpa must build the XSPEC module (default False)"),
                    ('xspec_version', None, "the XSPEC version (default 12.12.0)"),
                    ('xspec_lib_dirs', None, "Where the xspec libraries are located, if with_xspec is True"),
                    ('xspec_libraries', None, "Name of the libraries that should be linked for xspec"),
                    ('cfitsio_lib_dirs', None, "Where the cfitsio libraries are located, if with_xspec is True"),
                    ('cfitsio_libraries', None, "Name of the libraries that should be linked for cfitsio"),
                    ('ccfits_lib_dirs', None, "Where the CCfits libraries are located, if with_xspec is True"),
                    ('ccfits_libraries', None, "Name of the libraries that should be linked for CCfits"),
                    ('wcslib_lib_dirs', None, "Where the WCSLIB libraries are located, if with_xspec is True"),
                    ('wcslib_libraries', None, "Name of the libraries that should be linked for WCSLIB"),
                    ('gfortran_lib_dirs', None, "Where the gfortran libraries are located, if with_xspec is True"),
                    ('gfortran_libraries', None, "Name of the libraries that should be linked for gfortran"),
                    ]

    def initialize_options(self):
        self.with_xspec = False
        self.xspec_version = '12.14.1'
        self.xspec_include_dirs = ''
        self.xspec_lib_dirs = ''
        # This is set up for how CIAO builds XSPEC; other users may require more libraries
        self.xspec_libraries = 'XSFunctions XSUtil XS'
        self.cfitsio_include_dirs = ''
        self.cfitsio_lib_dirs = ''
        self.cfitsio_libraries = ''
        self.ccfits_include_dirs = ''
        self.ccfits_lib_dirs = ''
        self.ccfits_libraries = ''
        self.wcslib_include_dirs = ''
        self.wcslib_lib_dirs = ''
        self.wcslib_libraries = ''
        self.gfortran_include_dirs = ''
        self.gfortran_lib_dirs = ''
        self.gfortran_libraries = ''

    def finalize_options(self):
        pass

    def run(self):
        if not self.with_xspec:
            return

        macros = []
        if self.xspec_version:
            self.announce(f"Found XSPEC version: {self.xspec_version}", level=2)
            xspec_version = get_version(self.xspec_version)

            if xspec_version < MIN_VERSION:
                raise ValueError(f"XSPEC Version {xspec_version} is less than {MIN_VERSION}, which is the earliest supported version for Sherpa")

            for version in SUPPORTED_VERSIONS:
                if xspec_version >= version:
                    major, minor, micro = version
                    macros += [(f'XSPEC_{major}_{minor}_{micro}', None)]

            if xspec_version > MAX_VERSION:
                self.warn(f"XSPEC Version is greater than {MAX_VERSION}, which is the latest supported version for Sherpa")

        extension = build_ext(self, 'xspec', 'xspec', 'cfitsio', 'ccfits',
                              'wcslib', 'gfortran', define_macros=macros)
        self.distribution.ext_modules.append(extension)
