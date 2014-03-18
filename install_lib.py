#    ARCS Software 1.0
#
#    COPYRIGHT AND PERMISSION NOTICE
#    Copyright (c) 2006 California Institute of Technology.
#    All rights reserved.
#
#    Permission is hereby granted, free of charge, to any person obtaining a 
#    copy of this software and associated documentation files (the 
#    "Software"), to deal in the Software without restriction, including 
#    without limitation the rights to use, copy, modify, merge, publish, 
#    distribute, and/or sell copies of the Software, and to permit persons 
#    to whom the Software is furnished to do so, provided that the above 
#    copyright notice(s) and this permission notice appear in all copies of 
#    the Software and that both the above copyright notice(s) and this 
#    permission notice appear in supporting documentation.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
#    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
#    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT 
#    OF THIRD PARTY RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
#    HOLDERS INCLUDED IN THIS NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL 
#    INDIRECT OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING 
#    FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, 
#    NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION 
#    WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
#    Except as contained in this notice, the name of a copyright holder 
#    shall not be used in advertising or otherwise to promote the sale, use 
#    or other dealings in this Software without prior written authorization 
#    of the copyright holder.
#
#    All source code included in this distribution is covered by this notice,
##   unless specifically stated otherwise within each file. See each file within
#    each release for specific copyright holders.
#
#    ARCS is the name of an instrument under construction with U.S. DOE 
#    funding at Oak Ridge National Laboratory.
#


#import install
from distutils.command.install_lib import install_lib as _base

class install_lib(_base):

    # added a user option for c library installation directory.
    # distutils is supposed to only deal with python modules
    # and python extensions. It does not consider c libraries.
    # For a project like ARCS or DANSE, c libraries should
    # be built and put into a separate directory.
    # c library is not in the distutils SCHEME. We have to
    # overload install_lib to do this for us. This is
    # really a crappy approach :(
    user_options = _base.user_options
    user_options.append(
        ('clib-install-dir=', None, "directory to install c libraries to")
        )

    def run(self):
        _base.run(self)

        ## print out hints to set up environment variables
        ## this might not be necessary any more because
        ## in install.py we will save a sh script with
        ## env vars.
        print "* c libraries were installed into %s" % self.clib_install_dir
        #print "* python modules were installed into %s" % self.install_dir
        import os, sys
        if os.name == "nt":
            dll_ld_path = "PATH"
            pass
        elif sys.platform[:6] == "darwin":
            dll_ld_path = "DYLD_LIBRARY_PATH"
            pass
        else:
            dll_ld_path = "LD_LIBRARY_PATH"
            pass
        print "* You may need to add %s to environment variable %s."%(self.clib_install_dir, dll_ld_path)
        #print ">>> You may also need to change some environment variable so that compilers can find the right path to your c libraries."
        #print ">>> You may need to add %s to environment variable PYTHONPATH." % self.install_dir
        print
        return
        
    
    def initialize_options (self):
        #add an option
        #it is necessary to set any new option to None
        self.clib_install_dir = None
        _base.initialize_options(self)
        return


    def finalize_options (self):
        _base.finalize_options(self)
        
        self.set_undefined_options('install',
                                   ('install_base','clib_install_dir'))
        if self.clib_install_dir != None:
            import os
            if os.name == "nt":
                 #for windows installation, it is better to put dlls to
                 #a different directory than python Lib directory
                self.clib_install_dir = os.path.join( 
                    self.clib_install_dir, 'libs' )
            else:
                self.clib_install_dir = os.path.join( 
                    self.clib_install_dir, 'lib' )

        return

    def install(self):
        outfiles = _base.install(self)
        #distutils is used to install python modules, not c libraries
        #so it does not have a command to install c libs.
        #for a hack, we put installation of c libs here.
        if self.distribution.has_c_libraries():
            build_clib = self.get_finalized_command('build_clib')
            self.install_clibs( build_clib.libraries )
        return outfiles

    def install_clibs(self, libs):
        for (lib_name, build_info) in libs:
            self.install_clib(lib_name)
            pass
        return

    def install_clib( self, lib_name):
        #this is the part that deals with bulding c libraries
        build_clib = self.get_finalized_command('build_clib')
        #first it determines whether the library is a library
        #with headers only. if that is the case, no build is necessary
        if build_clib.isLibraryOnlyWithHeaders( lib_name ): return

        #darwin is different from other unix flavors
        if _isDarwin(): lib_type = "dylib"
        else: lib_type = "shared"
        filename = build_clib.compiler.library_filename( lib_name, lib_type = lib_type )
        self._install_clib( filename )
        #windows is always special!
        #install static lib for windows
        import os
        if os.name == "nt":
            filename = build_clib.compiler.library_filename( lib_name, lib_type = "static")
            self._install_clib( filename)
            pass
        return

    def _install_clib( self, lib_filename):
        build_clib = self.get_finalized_command('build_clib')
        from os.path import join
        dest = self.clib_install_dir
        self.mkpath( dest )
        self.copy_file( join( build_clib.build_clib, lib_filename ), dest )
        return

        
def _isDarwin():
    import sys
    return sys.platform[:6] == "darwin"
