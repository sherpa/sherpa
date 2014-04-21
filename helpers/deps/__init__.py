#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2014)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from subprocess import call
from multiprocessing import cpu_count
import os

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
    prefix=os.getcwd()
    os.chdir('extern')
    os.chmod(configure[0], 0755)
    out = call(configure)
    if out != 0: exit(out)
#    cflags = '-fPIC'
#    out = call(['make', 'CFLAGS='+cflags,'-j'+str(cpu_count()+1), 'install'])
    out = call(['make', '-j'+str(cpu_count()+1), 'install'])
    if out != 0: exit(out)
    open('built', 'w').close()
    os.chdir(prefix)