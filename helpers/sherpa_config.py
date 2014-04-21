#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2014)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from extensions import build_ext, build_lib_arrays
from deps import build_deps

from numpy.distutils.core import Command
import os

class sherpa_config(Command):
        description = "Configure Sherpa build options. If in doubt, ignore this command and stick to defaults. See setup.cfg for more information."
        user_options = [
                        ('fftw', None, "Whether Sherpa should build the embedded fftw3 library, which is the default behavior: set to 'local' to make Sherpa link against existing libraries on the system.)"),
                        ('fftw-include-dirs', None, "Where the fftw3 headers are located, if fftw is 'local'"),
                        ('fftw-lib-dirs', None, "Where the fftw3 libraries are located, if fftw is 'local'"),
                        ('fftw-libraries', None, "Name of the libraries that should be linked as fftw3"),
                        ('region', None, "Whether Sherpa should build the embedded region library, which is the default behavior: set to 'local' to make Sherpa link against existing libraries on the system.)"),
                        ('region-include-dirs', None, "Where the region headers are located, if region is 'local'"),
                        ('region-lib-dirs', None, "Where the region libraries are located, if region is 'local'"),
                        ('region-libraries', None, "Name of the libraries that should be linked as region"),
                        ('wcs-include-dirs', None, "Where the wcs subroutines headers are located"),
                        ('wcs-lib-dirs', None, "Where the wcs subroutines libraries are located"),
                        ('wcs-libraries', None, "Name of the libraries that should be linked as wcs"),
                        ('group-location', None, "Location of the group  python module"),
                        ('disable-group', None, "Disable the group module install"),
                        ('install-dir', None, "Directory where external dependencies must be installed (--prefix)"),
                        ('configure', None, "Additional configure flags for the external dependencies"),
                        ('group_cflags', None, "Additional cflags for building the grouping library"),
                        ]

        def initialize_options(self):
            self.install_dir=os.getcwd()+'/build'
            self.fftw=None
            self.fftw_include_dirs=None
            self.fftw_lib_dirs=None
            self.fftw_libraries='fftw3'
            self.region=None
            self.region_include_dirs=None
            self.region_lib_dirs=None
            self.region_libraries='region'
            self.wcs_include_dirs=None
            self.wcs_lib_dirs=None
            self.wcs_libraries='wcs'
            self.group_location=None
            self.disable_group=False
            self.configure='--enable-stuberrorlib --disable-shared --enable-shared=libgrp'
            self.group_cflags=None

        def finalize_options(self):
            if self.fftw_include_dirs is None:
                self.fftw_include_dirs=self.install_dir+'/include'

            if self.fftw_lib_dirs is None:
                self.fftw_lib_dirs=self.install_dir+'/lib'

            if self.region_include_dirs is None:
                self.region_include_dirs=self.install_dir+'/include'

            if self.region_lib_dirs is None:
                self.region_lib_dirs=self.install_dir+'/lib'

            if self.wcs_include_dirs is None:
                self.wcs_include_dirs=self.install_dir+'/include'

            if self.wcs_lib_dirs is None:
                self.wcs_lib_dirs=self.install_dir+'/lib'

            if self.group_location is None:
                self.group_location=self.install_dir+'/lib/python2.7/site-packages/group.so'

        def build_configure(self):
            configure = ['./configure', '--prefix='+self.install_dir, '--with-pic']
            if self.group_cflags is not None:
                configure.append('GROUP_CFLAGS="'+self.group_cflags+'"')
            if self.configure != 'None':
                configure.extend(self.configure.split(' '))
            if self.fftw != 'local':
                configure.append('--enable-fftw')
            self.distribution.ext_modules.append(build_ext('psf', *build_lib_arrays(self, 'fftw')))
            self.distribution.ext_modules.append(build_ext('wcs', *build_lib_arrays(self, 'wcs')))
            ld1, inc1, l1 = build_lib_arrays(self, 'wcs')
            if self.region != 'local':
                configure.append('--enable-region')
            ld2, inc2, l2 = build_lib_arrays(self, 'region')
            ld, inc, l = (ld1+ld2, inc1+inc2, l1+l2)
            self.distribution.ext_modules.append(build_ext('region', ld, inc, l))

#            if self.disable_group:
#                to_remove = None
#                for mod in self.distribution.ext_modules:
#                    if mod.name == 'group':
#                        to_remove = mod
#                if to_remove is not None:
#                    self.distribution.ext_modules.remove(to_remove)
#                    self.warn('Removing Group module from list of extension modules')

            if not self.disable_group:
                self.distribution.data_files.append(('', [self.group_location,]))

            self.warn('built configure string' + str(configure))

            return configure

        def run(self):
            configure = self.build_configure()
            build_deps(configure)


