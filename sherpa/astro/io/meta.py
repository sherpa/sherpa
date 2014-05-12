#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2009)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

import inspect
import sys
from sherpa.utils import NoNewAttributesAfterInit
from sherpa.utils.err import NotImplementedErr


__all__ = ('Meta',)


class Meta(NoNewAttributesAfterInit):

    def __init__(self):
        self.__header={}
        NoNewAttributesAfterInit.__init__(self)

    def __getitem__(self, name):
        return self.__header[name]

    def __setitem__(self, name, val):
        self.__header[name] = val

    def get(self, name, default=None):
        return self.__header.get(name, default)

    def pop(self, name):
        val = None
        try:
            val = self.__header.pop(name)
        except KeyError:
            return val

        return val

    def keys(self):
        keys = self.__header.keys()[:]
        keys.sort()
        return keys

    def has_key(self, key):
        return self.__header.has_key(key)

    def values(self):
        return [self.__header[key] for key in self.keys()]

    def copy(self):
        return self.__header.copy()

    def __repr__(self):
        return ("<%s instance>" % type(self).__name__)

    def __str__(self):
        s = ['']
        for name in self.keys():
            val = self.__header[name]

            if val is None or str(val).lower() == 'none':
                continue

            s.append(' %-13s = %s' % (name, val))

        return '\n'.join(s)
