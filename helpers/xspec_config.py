#
#  Copyright (C) 2014-2017, 2018  Smithsonian Astrophysical Observatory
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
from numpy.distutils.core import Command
from .extensions import build_ext, build_lib_arrays

def clean(xs):
    "Remove all '' entries from xs, returning the new list."
    return [x for x in xs if x != '']


def _get_xspec_version(libname, mode=0):
    """Call xs_getVersion provided by the library.

    Parameters
    ----------
    libname : str
        The full name to the library; passed to ctypes.CDLL
    mode : int, optional
        The mode parameter sent to ctypes.CDLL

    Returns
    -------
    ver : str or None
        Returns the version string of XSPEC or, if there is a problem,
        None.

    """

    import ctypes

    # This would likely need to be updated if we supported building outside
    # of Linux or OS-X.
    try:
        xslib = ctypes.CDLL(libname, mode=mode)
        version = b" "*256  # Is there a better way?
        xslib.xs_getVersion(version, len(version))
        return version.decode('ascii').rstrip(' \t\r\n\0')

    except:
        return None


def _find_xspec_version_1210(library_folders, library_extension,
                             cfitsio='cfitsio', ccfits='CCfits'):
    """Find the XSPEC version in XSPEC 12.10 (and hopefully later).

    Parameters
    ----------
    library_folders : list of str
        The directories to look in.
    library_extension : str
        The extenssion of the dylib (expected to be 'so' or 'dylib').
    cfitsio , ccfits : str or None, optional
        The name of the cfitsio and ccfits libraries (this is the library
        name without the "lib" prefix and the ".xxx" suffix). If either
        is set to None then no attempt to find the version is made.

    Returns
    -------
    version : str or None
        None is returned if the version could not be found.

    """

    import ctypes

    if cfitsio is None or ccfits is None:
        return None

    # Following the approach of
    # https://github.com/sherpa/sherpa/issues/436#issuecomment-371944527
    # which requires installing a lot of symbols into the global name space.
    #
    # This feels very fragile.
    #
    for name in [cfitsio, ccfits, 'XSUtil']:
        for folder in library_folders:
            libname = os.path.join(folder,
                                   "lib{}.{}".format(name, library_extension))
            try:
                ctypes.CDLL(libname, mode=ctypes.RTLD_GLOBAL)
                break
            except:
                pass

        else:
            # Unable to load this library; no point in continuing
            return None

    library_base_name = 'libXSFunctions'
    library_name = '{}.{}'.format(library_base_name, library_extension)

    for folder in library_folders:
        libname = os.path.join(folder, library_name)
        ver = _get_xspec_version(libname, mode=1)
        if ver is not None:
            return ver

    return None


def _find_xspec_version_129x(library_folders, library_extension):
    """Find the XSPEC version in XSPEC 12.9.x.

    Parameters
    ----------
    library_folders : list of str
        The directories to look in.
    library_extension : str
        The extenssion of the dylib (expected to be 'so' or 'dylib').

    Returns
    -------
    version : str or None
        None is returned if the version could not be found.

    """

    library_base_name = 'libXSUtil'
    library_name = '{}.{}'.format(library_base_name, library_extension)

    for folder in library_folders:
        libname = os.path.join(folder, library_name)
        ver = _get_xspec_version(libname)
        if ver is not None:
            return ver

    return None


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

            xspec_raw_version = self._find_xspec_version(ld,
                                                         cfitsio=cfitsio,
                                                         ccfits=ccfits
                                                         )

            macros = []

            if xspec_raw_version:
                self.announce("Found XSPEC version: {}".format(xspec_raw_version), 2)
                xspec_version = LooseVersion(xspec_raw_version)

                if xspec_version >= LooseVersion("12.9.0"):
                    macros += [('XSPEC_12_9_0', None)]
                else:
                    self.warn("XSPEC Version is less than 12.9.0, which is the minimal supported version for Sherpa")

                if xspec_version >= LooseVersion("12.9.1"):
                    macros += [('XSPEC_12_9_1', None)]

                if xspec_version >= LooseVersion("12.10.0"):
                    macros += [('XSPEC_12_10_0', None)]

                # Since there are patches (e.g. 12.10.0c), look for the
                # "next highest version.
                if xspec_version >= LooseVersion("12.10.1"):
                    self.warn("XSPEC Version is greater than 12.10.0, which is the latest supported version for Sherpa")

            extension = build_ext('xspec', ld, inc, l, define_macros=macros)

            self.distribution.ext_modules.append(extension)

        else:
            if package in dist_packages:
                dist_packages.remove(package)
            if package in dist_data:
                del dist_data[package]

    def _find_xspec_version(self, library_folders,
                            cfitsio=None,
                            ccfits=None):
        """Try and find the XSPEC version."""

        # Determining the full library name according to the platform.
        # I couldn't find a simpler, more portable way of doing this.
        if sys.platform == 'darwin':
            library_extension = 'dylib'
        else:  # assume linux
            library_extension = 'so'

        ver = _find_xspec_version_1210(library_folders, library_extension,
                                       cfitsio=cfitsio, ccfits=ccfits)
        if ver is not None:
            return ver

        return _find_xspec_version_129x(library_folders, library_extension)
