#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2014)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from build import build
from clean import clean
from install import install
from sdist import sdist
from sherpa_config import sherpa_config
from xspec_config import xspec_config

commands = {
             'build': build,
             'clean' : clean,
             'install' : install,
             'sdist' : sdist,
             'sherpa_config' : sherpa_config,
             'xspec_config' : xspec_config,
            }