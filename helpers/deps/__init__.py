# 
#  Copyright (C) 2014  Smithsonian Astrophysical Observatory
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


from subprocess import call
from multiprocessing import cpu_count
import os
import sys

def clean_deps():
    prefix = os.getcwd()
    os.chdir('extern')
    call(['make', 'uninstall'])
    call(['make', 'distclean'])
    try:
        os.remove('built')
    except:
        pass
    os.chdir(prefix)

def build_deps(configure):
    if not os.path.exists('extern/built'):
        prefix=os.getcwd()
        os.chdir('extern')
        os.chmod(configure[0], 0o755)
        env = os.environ.copy()
        env['PYTHON'] = sys.executable
        out = call(configure, env=env)
        if out != 0: exit(out)
    #    cflags = '-fPIC'
    #    out = call(['make', 'CFLAGS='+cflags,'-j'+str(cpu_count()+1), 'install'])
        out = call(['make', '-j'+str(cpu_count()+1), 'install'])
        if out != 0: exit(out)
        open('built', 'w').close()
        os.chdir(prefix)
