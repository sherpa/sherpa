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
    print(f">> changing directory: {owd} to {cwd}"); sys.stdout.flush()
    try:
        yield

    finally:
        print(f"<< restoring directory: {owd}"); sys.stdout.flush()
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
    print(">> In configure step"); sys.stdout.flush()
    if os.path.exists('extern/built'):
        return

    print(f">> about to change to extern from {os.getcwd()}"); sys.stdout.flush()
    with in_dir("extern"):
        print(">>>"); sys.stdout.flush()
        os.system("ls -l")
        print("<<<"); sys.stdout.flush()
        print(f">> configure:    {configure}")
        print(f">> configure[0]: {configure[0]}")

        # If this succeeds then configure should exist ...
        os.chmod(configure[0], 0o755)

        print(">>>>>>"); sys.stdout.flush()
        os.system("ls -l configure")
        print("<<<<<<"); sys.stdout.flush()

        env = os.environ.copy()
        print(f">> python: {env.get('PYTHON')}}  -> {sys.executable}")
        env['PYTHON'] = sys.executable
        env['PWD'] = os.getcwd()  # is this needed
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
