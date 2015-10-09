# 
#  Copyright (C) 2009  Smithsonian Astrophysical Observatory
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
