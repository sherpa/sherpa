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

import os

from numpy.distutils.core import Command
from .extensions import build_ext, build_lib_arrays

def clean(xs):
    "Remove all '' entries from xs, returning the new list."
    return [x for x in xs if x != '']

# Try to support "optional" models; that is, those models
# that are added (or deleted) in a version. At present
# (XSPEC 12.9.0) this is only the nlapec model, so special
# case it rather than come up with a generic method.
#
# The current approach to determining if the nlapec model is available
# is to do a quick query on the model.dat file. However, this is
# not specified in the build configuration - the best guess we have
# here are:
#  1) 'xspec-lib-dirs' + '../../spectral/manager/model.dat'
#  2) $HEADAS + '../spectral/manager/model.dat'
#
# For now, the former approach is taken. This is not as flexible,
# but requires less changes to the current build (since, at present,
# there is no explicit requirement that the HEADAS environment variable
# is set up for the build - although perhaps there should be).
#

def does_nlapec_model_exist(libdir):
    """Can we find nlapec model in model.dat?

    Parameters
    ----------
    libdir : str
        The location of the XSPEC libraries. The model.dat file is
        assumed to be located relative to this directory:
        ../spectral/manager/model.dat.

    Returns
    -------
    flag : bool
        True if model.dat can be found and nlapec is found in it,
        False otherwise.
    """

    fname = os.path.join(libdir, '../../spectral/manager/model.dat')
    try:
        cts = open(fname, 'r').read()
    except IOError:
        return False

    # Do not use a full parser for the model.dat file
    for l in cts.split("\n"):
        toks = l.split()
        # limited check that we have found the nlapec model line
        if len(toks) < 7 or toks[0] != 'nlapec' or toks[4] != 'C_nlapec':
            continue
        return True

    return False


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

            # Assume the first element in XSPEC library path
            # is the one to use.
            if does_nlapec_model_exist(ld1[0]):
                macros = [('XSPEC_HAS_NLAPEC_MODEL', None)]
            else:
                macros = None

            ext = build_ext('xspec', ld, inc, l, define_macros=macros)
            self.distribution.ext_modules.append(ext)

        else:
            if package in dist_packages:
                dist_packages.remove(package)
            if package in dist_data:
                del dist_data[package]
