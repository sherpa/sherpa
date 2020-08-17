#
#  Copyright (C) 2014-2017, 2018, 2020
#       Smithsonian Astrophysical Observatory
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

import sys
import os

from distutils.version import LooseVersion
from distutils.cmd import Command
from .extensions import build_ext, build_lib_arrays

def clean(xs):
    "Remove all '' entries from xs, returning the new list."
    return [x for x in xs if x != '']


class xspec_config(Command):
    description = "Configure XSPEC Models external module (optional) "
    user_options = [
                    ('with-xspec', None, "Whether sherpa must build the XSPEC module (default False)"),
                    ('xspec-version', None, "the XSPEC version (default 12.9.0)"),
                    ('xspec-lib-dirs', None, "Where the xspec libraries are located, if with-xspec is True"),
                    ('xspec-libraries', None, "Name of the libraries that should be linked for xspec"),
                    ('cfitsio-lib-dirs', None, "Where the cfitsio libraries are located, if with-xspec is True"),
                    ('cfitsio-libraries', None, "Name of the libraries that should be linked for cfitsio"),
                    ('ccfits-lib-dirs', None, "Where the CCfits libraries are located, if with-xspec is True"),
                    ('ccfits-libraries', None, "Name of the libraries that should be linked for CCfits"),
                    ('wcslib-lib-dirs', None, "Where the WCSLIB libraries are located, if with-xspec is True"),
                    ('wcslib-libraries', None, "Name of the libraries that should be linked for WCSLIB"),
                    ('gfortran-lib-dirs', None, "Where the gfortran libraries are located, if with-xspec is True"),
                    ('gfortran-libraries', None, "Name of the libraries that should be linked for gfortran"),
                    ]

    def initialize_options(self):
        self.with_xspec = False
        self.xspec_version = '12.9.0'
        self.xspec_include_dirs = ''
        self.xspec_lib_dirs = ''
        self.xspec_libraries = 'XSFunctions XSModel XSUtil XS'
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
        package = 'sherpa.astro.xspec'
        dist_packages = self.distribution.packages
        dist_data = self.distribution.package_data

        if self.with_xspec:
            if package not in dist_packages:
                dist_packages.append(package)

            if package not in dist_data:
                dist_data[package] = ['tests/test_*.py']

            ld1, inc1, l1 = build_lib_arrays(self, 'xspec')
            ld2, inc2, l2 = build_lib_arrays(self, 'cfitsio')
            ld3, inc3, l3 = build_lib_arrays(self, 'ccfits')
            ld4, inc4, l4 = build_lib_arrays(self, 'wcslib')
            ld5, inc5, l5 = build_lib_arrays(self, 'gfortran')

            ld = clean(ld1 + ld2 + ld3 + ld4 + ld5)
            inc = clean(inc1 + inc2 + inc3 + inc4 + inc5)
            l = clean(l1 + l2 + l3 + l4 + l5)

            # I do not know if l2/l3 are guaranteed to be indexable
            # entries with at least one element in them.
            #
            try:
                cfitsio = l2[0]
            except IndexError:
                cfitsio = None

            try:
                ccfits = l3[0]
            except IndexError:
                ccfits = None

            xspec_raw_version = self.xspec_version

            macros = []

            if xspec_raw_version:
                self.announce("Found XSPEC version: {}".format(xspec_raw_version), 2)
                xspec_version = LooseVersion(xspec_raw_version)

                if xspec_version < LooseVersion("12.9.0"):
                    self.warn("XSPEC Version is less than 12.9.0, which is the minimal supported version for Sherpa")

                # I am not sure what the naming of the XSPEC components are,
                # but let's stick with major, minor, and patch.
                #
                for major, minor, patch in [(12, 9, 0), (12, 9, 1),
                                            (12, 10, 0), (12, 10, 1),
                                            (12, 11, 0), (12, 11, 1)]:
                    version = '{}.{}.{}'.format(major, minor, patch)
                    macro = 'XSPEC_{}_{}_{}'.format(major, minor, patch)
                    if xspec_version >= LooseVersion(version):
                        macros += [(macro, None)]

                # Since there are patches (e.g. 12.10.0c), look for the
                # "next highest version.
                if xspec_version >= LooseVersion("12.11.2"):
                    self.warn("XSPEC Version is greater than 12.11.1, which is the latest supported version for Sherpa")

            extension = build_ext('xspec', ld, inc, l, define_macros=macros)

            self.distribution.ext_modules.append(extension)

        else:
            if package in dist_packages:
                dist_packages.remove(package)
            if package in dist_data:
                del dist_data[package]
