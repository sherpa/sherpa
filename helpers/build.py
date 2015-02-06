#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2014)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from numpy.distutils.command.build import build as _build
from deps import build_deps

import os

class build(_build):
    def run(self):
        configure = self.get_finalized_command('sherpa_config', True).build_configure()
        self.get_finalized_command('xspec_config', True).run()
        if not os.path.exists('extern/built'):
            build_deps(configure)
        _build.run(self)