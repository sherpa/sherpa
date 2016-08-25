#
#  Copyright (C) 2014, 2015, 2016  Smithsonian Astrophysical Observatory
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


from numpy.distutils.core import Command
from .extensions import build_ext, build_lib_arrays

def clean(xs):
    "Remove all '' entries from xs, returning the new list."
    return [x for x in xs if x != '']

class xspec_config(Command):
    description = "Configure XSPEC Models external module (optional) "
    user_options = [
                    ('with-xspec', None, "Whether sherpa must build the XSPEC module (default False)"),
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
        self.xspec_include_dirs = ''
        self.xspec_lib_dirs = ''
        self.xspec_libraries = 'XSFunctions XSModel XSUtil XS'
        self.cfitsio_include_dirs = ''
        self.cfitsio_lib_dirs = ''
        self.cfitsio_libraries = 'cfitsio'
        self.ccfits_include_dirs = ''
        self.ccfits_lib_dirs = ''
        self.ccfits_libraries = 'CCfits'
        self.wcslib_include_dirs = ''
        self.wcslib_lib_dirs = ''
        self.wcslib_libraries = 'wcs'
        self.gfortran_include_dirs = ''
        self.gfortran_lib_dirs = ''
        self.gfortran_libraries = 'gfortran'

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

            self.distribution.ext_modules.append(build_ext('xspec', ld, inc, l))

        else:
            if package in dist_packages:
                dist_packages.remove(package)
            if package in dist_data:
                del dist_data[package]
