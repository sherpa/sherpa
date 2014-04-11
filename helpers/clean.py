from distutils.command.clean import clean as _clean
from deps import clean_deps

class clean(_clean):
    def run(self):
        _clean.run(self)
        clean_deps()