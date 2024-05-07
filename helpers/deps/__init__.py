#
#  Copyright (C) 2014, 2016, 2022
#  Smithsonian Astrophysical Observatory
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


from contextlib import contextmanager
from multiprocessing import cpu_count
import os
from subprocess import call
import sys


@contextmanager
def in_dir(cwd):
    """Run the context in the given directory, which is assumed to exist."""

    owd = os.getcwd()
    os.chdir(cwd)
    try:
        yield

    finally:
        os.chdir(owd)


def make(command, opts=None, check=False):
    """Call make with the given command and possible options."""

    args = ["make"]
    if opts is not None:
        args.extend(opts)

    args.append(command)
    rval = call(args)
    if check and rval != 0:
        sys.exit(1)


def clean_deps():
    """Ensure the dependencies are cleaned"""

    with in_dir("extern"):
        make("uninstall")
        make("distclean")
        try:
            os.remove('built')
        except:
            pass


def build_deps(configure):
    if os.path.exists('extern/built'):
        return

    with in_dir("extern"):
        os.chmod(configure[0], 0o755)

        env = os.environ.copy()
        env['PYTHON'] = sys.executable
        out = call(configure, env=env)
        if out != 0:
            sys.exit(out)

        counter = cpu_count() + 1
        # cflags = '-fPIC'
        opts = [
            # "CFLAGS={cflags}",
            f"-j{counter}"
        ]
        make("install", opts, check=True)

        # Create the built file to indicate the dependencies have been
        # created.
        #
        open('built', 'w').close()
