from __future__ import print_function
#
#  Copyright (C) 2015, 2016  Smithsonian Astrophysical Observatory
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


def set_xspec_chatter():
    try:
        from sherpa.astro import xspec
        xspec.set_xschatter(50)
    except ImportError:
        # Well, it looks like xspec is not available. Passing.
        pass

try:
    from setuptools.command.test import test

    import sys


    class PyTest(test):
        user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

        def initialize_options(self):
            test.initialize_options(self)
            self.pytest_args = []

        def finalize_options(self):
            test.finalize_options(self)
            self.test_args = []
            self.test_suite = True
            if not self.pytest_args:
                self.pytest_args = 'sherpa'
            self.pytest_args = self.pytest_args.split(' ')

        def run_tests(self):
            set_xspec_chatter()
            # import here, cause outside the eggs aren't loaded
            import pytest
            errno = pytest.main(self.pytest_args)
            sys.exit(errno)

except ImportError:
    from distutils.cmd import Command

    class PyTest(Command):
        user_options = []

        def initialize_options(self):
            pass

        def finalize_options(self):
            pass

        def run(self):
            print("test command is not available without setuptools")
