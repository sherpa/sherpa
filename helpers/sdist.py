#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2014)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from distutils.command.sdist import sdist as _sdist
from numpy.distutils.misc_util import get_data_files
from deps import clean_deps

class sdist(_sdist):

    def add_defaults(self):
        _sdist.add_defaults(self)

        dist = self.distribution

        if dist.has_data_files():
            for data in dist.data_files:
                self.filelist.extend(get_data_files(data))

        if dist.has_headers():
            headers = []
            for h in dist.headers:
                if isinstance(h,str): headers.append(h)
                else: headers.append(h[1])
            self.filelist.extend(headers)

        return

    def run(self):
        clean_deps()
        self.get_finalized_command('sherpa_config', True).build_configure()
        _sdist.run(self)