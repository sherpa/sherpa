from __future__ import print_function
#
#  Copyright (C) 2007,2014,2015,2016  Smithsonian Astrophysical Observatory
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

"""

Modeling and fitting package for scientific data analysis

Sherpa is a modeling and fitting package for scientific data analysis.
It includes general tools of interest to all users as well as
specialized components for particular disciplines (e.g. astronomy).

Note that the top level sherpa package does not import any
subpackages.  This allows the user to import a particular component
(e.g. `sherpa.optmethods`) without pulling in any others.  To import all
the standard subpackages, use ``import sherpa.all`` or
``from sherpa.all import *``.

"""


import logging
import os
import os.path
import sys


__all__ = ('get_include',)

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


class Formatter(logging.Formatter):
    def format(self, record):
        if record.levelno > logging.INFO:
            msg = '%s: %s' % (record.levelname, record.msg)
        else:
            msg = record.msg
        return msg

log = logging.getLogger('sherpa')
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(Formatter())
log.addHandler(handler)
log.setLevel(logging.INFO)

del Formatter, log, handler


def get_include():
    "Get the root path for installed Sherpa header files"

    return os.path.join(os.path.dirname(__file__), 'include')


def get_config():
    "Get the path for the installed Sherpa configuration file"

    filename = "sherpa-standalone.rc"

    home_dir = None
    config = None

    # If NOSHERPARC is set, read in system config file
    # ignore any user config file
    if (('NOSHERPARC' in os.environ) == True):
        return os.path.join(os.path.dirname(__file__), filename)

    # If SHERPARC is set, read in config file from there,
    # and ignore default location
    if (('SHERPARC' in os.environ) == True):
        config = os.environ.get('SHERPARC')
        if os.path.isfile(config):
            return config

    # SHERPARC was not set, so look for .sherpa.rc in default
    # location, which is user's home directory.
    home_dir = os.environ.get('HOME')
    config = os.path.join(home_dir, '.'+filename)

    if os.path.isfile(config):
        return config

    # If no user config file is set, fall back to system config file
    return os.path.join(os.path.dirname(__file__), filename)


def smoke(verbosity=0, require_failure=False, fits=None, xspec=False, ds9=False):
    """
    Run Sherpa's "smoke" test. The smoke test is a simple test that
        ensures the Sherpa installation is functioning. It is not a complete
        test suite, but it fails if obvious issues are found.

    Parameters
    ----------
    xspec : boolean
        Require xspec module when running tests. Tests requiring xspec may still run if the xspec module is present.
    fits : str
        Require a fits module with this name to be present before running the smoke test.
        This option makes sure that when the smoke test is run the required modules are present.
        Note that tests requiring fits may still run if any fits backend is available, and they might
        still fail on their own.
    require_failure : boolean
        For debugging purposes, the smoke test may be required to always fail. Defaults to False.
    verbosity : int
        The level of verbosity of this test

    Returns
    -------
    The method raises the SystemExit if errors are found during the smoke test
    """
    from sherpa.astro.utils import smoke
    smoke.run(verbosity=verbosity, require_failure=require_failure, fits=fits, xspec=xspec, ds9=ds9)


def _smoke_cli(verbosity=0, require_failure=False, fits=None, xspec=False, ds9=False):
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-v", "--verbosity", dest="verbosity",
                      help="verbosity level")
    parser.add_option("-0", "--require-failure", dest="require_failure", action="store_true",
                      help="require smoke test to fail (for debugging)")
    parser.add_option("-f", "--fits-module", dest="fits", action="store",
                      help="require a specific fits module to be present")
    parser.add_option("-x", "--require-xspec", dest="xspec", action="store_true",
                      help="require xspec module")
    parser.add_option("-d", "--require-ds9", dest="ds9", action="store_true",
                      help="require DS9")

    options, _ = parser.parse_args()

    xspec = options.xspec or xspec
    verbosity = options.verbosity or verbosity
    require_failure = options.require_failure or require_failure
    fits = options.fits or fits
    ds9 = options.ds9 or ds9

    smoke(verbosity=verbosity, require_failure=require_failure, fits=fits, xspec=xspec, ds9=ds9)


def _install_test_deps():
    def install(package_name):
        try:
            import pip
            pip.main(['install', package_name])
        except:
            print("""Cannot import pip or install packages with it.
            You need pytest, and possibly pytest-cov, in order to run the tests.
            If you downloaded the source code, please run 'pip install -r test_requirements.txt'
            from the source directory first.
            """)
            raise

    deps = ['pytest', 'mock']
    pytest_plugins =  ['pytest-catchlog',]

    installed_plugins = []

    for dep in deps:
        try:
            __import__(dep)
        except ImportError:
            install(dep)

    for plugin_name in pytest_plugins:
        module = plugin_name.replace("-", "_")
        try:
            __import__(module)
        except ImportError:
            install(plugin_name)
            installed_plugins.append(module)

    return installed_plugins


def clitest():
    plugins = _install_test_deps()
    import pytest
    import os
    sherpa_dir = os.path.dirname(__file__)
    errno = pytest.main([sherpa_dir, '-rs'], plugins=plugins)
    sys.exit(errno)
