#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007,2014)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
"""

Modeling and fitting package for scientific data analysis

Sherpa is a modeling and fitting package for scientific data analysis.
It includes general tools of interest to all users as well as
specialized components for particular disciplines (e.g. astronomy).

Note that the top level sherpa package does not import any
subpackages.  This allows the user to import a particular component
(e.g. sherpa.optmethods) without pulling in any others.  To import all
the standard subpackages, use 'import sherpa.all' or
'from sherpa.all import *'.

"""


import logging
import os
import os.path
import sys



__all__ = ('banner', 'get_include', 'test')

__versionstr__ = '4.7.1'
__version__ = 40701

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


_banner = """
-----------------------------------------------------
Welcome to Sherpa: CXC's Modeling and Fitting Package
-----------------------------------------------------
Version: Sherpa %s Monday, July 19, 2010
""" % __versionstr__


def _banner_fancy(file=sys.stdout):
    "Print a welcome message to the specified file object"
    print >> file, _banner

def banner(file=sys.stdout):
    "No-op function; override with fancier version if desired to print banner"
    pass

def get_include():
    "Get the root path for installed Sherpa header files"
    
    return os.path.join(os.path.dirname(__file__), 'include')

def get_config():
    "Get the path for the installed Sherpa configuration file"

    home_dir = None
    config = None

    # If NOSHERPARC is set, read in system config file
    # ignore any user config file
    if (os.environ.has_key('NOSHERPARC') == True):
        return os.path.join(os.path.dirname(__file__), 'sherpa.rc')
    
    # If SHERPARC is set, read in config file from there,
    # and ignore default location
    if (os.environ.has_key('SHERPARC') == True):
        config = os.environ.get('SHERPARC')
        if os.path.isfile(config):
            return config
        
    # SHERPARC was not set, so look for .sherpa.rc in default
    # location, which is user's home directory.
    home_dir = os.environ.get('HOME')
    config = os.path.join(home_dir, '.sherpa.rc')

    if os.path.isfile(config):
        return config

    # If no user config file is set, fall back to system config file
    return os.path.join(os.path.dirname(__file__), 'sherpa.rc')

    
def test(level=1, verbosity=1, datadir=None):
    """
    
    Run the Sherpa test suite, testing all available subpackages
    (including the discipline-specific ones)
    
    """
    # import sherpa.all
    # import sherpa.astro.all
    from sherpa.utils import SherpaTest
    SherpaTest().test(level, verbosity, datadir)

def clitest():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-l", "--level", dest="level",
                  help="Test Level", default=1)
    parser.add_option("-v", "--verbosity", dest="verbosity",
                  help="Test Verbosity", default=1)
    parser.add_option("-d", "--datadir", dest="datadir",
                  help="Test Level", default=None)
    (options, _) = parser.parse_args()
    test(options.level, options.verbosity, options.datadir)

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
