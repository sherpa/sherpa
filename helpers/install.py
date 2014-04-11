#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2014)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from numpy.distutils.command.install import install as _install
import os

from deps import build_deps

class install(_install):
    def run(self):
        if not os.path.exists('extern/built'):
            configure = self.get_finalized_command('sherpa_config', True).build_configure()
            self.get_finalized_command('xspec_config', True).run()
            build_deps(configure)
        _install.run(self)