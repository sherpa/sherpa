#
#  Copyright (C) 2014 - 2016, 2020, 2022, 2024
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
import os
import sys

from .extensions import build_ext
from .deps import build_deps

version = f'{sys.version_info[0]}.{sys.version_info[1]}'


class sherpa_config(Command):
    description = "Configure Sherpa build options. If in doubt, ignore this command and stick to defaults. See setup.cfg for more information."
    user_options = [
                    ('fftw', None, "Whether Sherpa should build the embedded fftw3 library, which is the default behavior: set to 'local' to make Sherpa link against existing libraries on the system.)"),
                    ('fftw-include-dirs', None, "Where the fftw3 headers are located, if fftw is 'local'"),
                    ('fftw-lib-dirs', None, "Where the fftw3 libraries are located, if fftw is 'local'"),
                    ('fftw-libraries', None, "Name of the libraries that should be linked as fftw3"),
                    ('region', None, "Whether Sherpa should build the embedded region library, which is the default behavior: set to 'local' to make Sherpa link against existing libraries on the system.)"),
                    ('region-include-dirs', None, "Where the region headers are located, if region is 'local'"),
                    ('region-lib-dirs', None, "Where the region libraries are located, if region is 'local'"),
                    ('region-libraries', None, "Name of the libraries that should be linked as region"),
                    ('region-use-cxc-parser', None, "If set to True than use the CXC Data Model library to parse region files"),
                    ('wcs', None, "Whether Sherpa should build the embedded wcs library, which is the default behavior: set to 'local' to make Sherpa link against existing libraries on the system.)"),
                    ('wcs-include-dirs', None, "Where the wcs subroutines headers are located"),
                    ('wcs-lib-dirs', None, "Where the wcs subroutines libraries are located"),
                    ('wcs-libraries', None, "Name of the libraries that should be linked as wcs"),
                    ('group-location', None, "Location of the group  python module"),
                    ('disable-group', None, "Disable the group module install"),
                    ('stk-location', None, "Location of the stack library python module"),
                    ('disable-stk', None, "Disable the stack library module install"),
                    ('install-dir', None, "Directory where external dependencies must be installed (--prefix)"),
                    ('configure', None, "Additional configure flags for the external dependencies"),
                    ('group-cflags', None, "Additional cflags for building the grouping library"),
                    ]

    def initialize_options(self):
        self.install_dir = os.path.join(os.getcwd(), 'build')
        self.fftw = None
        self.fftw_include_dirs = None
        self.fftw_lib_dirs = None
        self.fftw_libraries = 'fftw3'
        self.region = None
        self.region_include_dirs = None
        self.region_lib_dirs = None
        self.region_libraries = 'region'
        self.region_use_cxc_parser = False
        self.wcs = None
        self.wcs_include_dirs = None
        self.wcs_lib_dirs = None
        self.wcs_libraries = 'wcssubs'
        self.group_location = None
        self.disable_group = False
        self.configure = '--disable-maintainer-mode --enable-stuberrorlib --disable-shared --enable-shared=libgrp,stklib'
        self.group_cflags = None
        self.stk_location = None
        self.disable_stk = False

    def finalize_options(self):
        incdir = os.path.join(self.install_dir, 'include')
        libdir = os.path.join(self.install_dir, 'lib')
        pydir = os.path.join(libdir, f'python{version}', 'site-packages')

        if self.fftw_include_dirs is None:
            self.fftw_include_dirs = incdir

        if self.fftw_lib_dirs is None:
            self.fftw_lib_dirs = libdir

        if self.region_include_dirs is None:
            self.region_include_dirs = incdir

        if self.region_lib_dirs is None:
            self.region_lib_dirs = libdir

        # It is not clear that the other boolean options work
        # correctly when set to False/false.
        #
        if not isinstance(self.region_use_cxc_parser, bool):
            self.region_use_cxc_parser = self.region_use_cxc_parser.lower() == "true"

        if self.wcs_include_dirs is None:
            self.wcs_include_dirs = incdir

        if self.wcs_lib_dirs is None:
            self.wcs_lib_dirs = libdir

        if self.group_location is None:
            self.group_location = os.path.join(pydir, 'group.so')

        if self.stk_location is None:
            self.stk_location = os.path.join(pydir, 'stk.so')

    def build_configure(self):

        # can we find ascdm.h (if so we can support FITS region files as long
        # as the other settings are set up to include ascdm)?
        #
        region_macros = None
        if self.region_use_cxc_parser:
            region_macros = [("USE_CXCDM_PARSER", None)]

        psfext = build_ext(self, 'psf', 'fftw')
        wcsext = build_ext(self, 'wcs')
        regext = build_ext(self, 'region', define_macros=region_macros)
        self.distribution.ext_modules.append(psfext)
        self.distribution.ext_modules.append(wcsext)
        self.distribution.ext_modules.append(regext)

        configure = ['./configure', '--prefix=' + self.install_dir,
                     '--with-pic', '--enable-standalone']

        if self.group_cflags is not None:
            configure.append(f'GROUP_CFLAGS="{self.group_cflags}"')
        if self.configure != 'None':
            configure.extend(self.configure.split(' '))
        if self.fftw != 'local':
            configure.append('--enable-fftw')
        if self.region != 'local':
            configure.append('--enable-region')

        # Note that we hack in the "correct" location for the group
        # and stack libraries. This appears to only get used for
        # normal installs: editable installs do not seem to care about
        # the data_files setting.
        #
        libdir = os.path.join('lib',
                              f'python{version}',
                              'site-packages')
        dfiles = []
        if not self.disable_group:
            configure.append('--enable-group')
            dfiles.append((libdir, [self.group_location, ]))

        if not self.disable_stk:
            configure.append('--enable-stk')
            dfiles.append((libdir, [self.stk_location, ]))

        if dfiles:
            if self.distribution.data_files is None:
                self.distribution.data_files = dfiles
            else:
                self.distribution.data_diles.extend(dfiles)

        if self.wcs != 'local':
            configure.append('--enable-wcs')

        self.warn(f'built configure string {configure}')
        return configure

    def run(self):
        configure = self.build_configure()
        build_deps(configure)
