# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                   Jiao Lin
#                      California Institute of Technology
#              (C) 2005 All Rights Reserved  All Rights Reserved
#
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
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


import _addons
from numpy.distutils.command.build_clib import build_clib as _base
from distutils.command.build_clib import ListType, TupleType, log


class build_clib( _base ):

    def __init__(self, *args, **kwds):
        _base.__init__(self, *args, **kwds)
        # this is needed by _addons._copy_header_tree(s)
        self.outfiles = []
        # some "library" only has headers in it. we need to put a flag on
        # such libraries
        self._headersOnly = {}
        return


    def finalize_options (self):
        # disutils by default set the building directory of clib to build_temp
        # because clib is thought to be temp we want it to be a little more
        # than temp
        if self.build_clib is None:
            # copy the distutils way to create directory name
            from distutils.util import get_platform
            import sys
            plat_specifier = ".%s-%s" % (get_platform(), sys.version[0:3])
            import os
            # the build tree base path can be obtained from command object
            # "build", which will always be called when doing
            #   setup.py build or setup.py install
            build = self.distribution.get_command_obj( "build" )
            #builiding tree for c library is now under clib.<platform name>
            self.build_clib = os.path.join(build.build_base,
                                           'clib' + plat_specifier)
            pass
        # call the original class method
        _base.finalize_options(self)
        return


    def isLibraryOnlyWithHeaders( self, libname ):
        # determine if the given library is only with header files (no sources)
        return self._headersOnly.get( libname )


    def build_libraries (self, libraries):
        # overload this function to change the way how libraries are
        # built

        # build
        for (lib_name, build_info) in libraries:
            self._build_library(lib_name, build_info)
            pass
        return


    def get_library_names (self):
        #overload this function so that it can deal with
        #header-only libraries correctly
        if not self.libraries:
            return None

        lib_names = []
        for (lib_name, build_info) in self.libraries:
            #only "real" libraries are returned
            #libaries without sources are ignored
            if not self._headersOnly.get(lib_name):
                lib_names.append(lib_name)
        return lib_names


    def _build_library (self, lib_name, build_info):
        #build one library
        #get sources
        sources = build_info.get('sources')
        #make sure "sources" is a list
        sources = self._check_source_list( lib_name, sources)

        # added building directory to lib dirs so that later-built libraries
        # can use libraries built earlier
        build_info['libdirs'].append(self.build_clib)

        #build
        if len(sources): # with sources
            log.info("building '%s' library", lib_name)
            #compile sources 
            objects = self._build_objs( lib_name, sources, build_info)
            #link
            self._link_objs( lib_name, objects, build_info)
            #
            self._headersOnly[lib_name] = False
        else: #no sources, only headers
            self._headersOnly[lib_name] = True
            pass
        
        #copy headers to the build tree
        #so that binding codes can find them
        self._copy_header_trees(
            self.distribution.headers, self.build_clib)
        return


    def _check_source_list(self, lib_name, sources):
        #this is the standard distutils way to check things
        #I just copy and paste 
        if sources is None or type(sources) not in (ListType, TupleType):
            raise DistutilsSetupError, \
                  ("in 'libraries' option (library '%s'), " +
                   "'sources' must be present and must be " +
                   "a list of source filenames") % lib_name
        return list(sources)


    def _build_objs(self, lib_name, sources, build_info):
        # First, compile the source code to object files in the library
        # directory.  (This should probably change to putting object
        # files in a temporary build directory.)
        macros = build_info.get('macros')
        include_dirs = build_info.get('include_dirs')
        #add the build temp tree to include so that headers can be
        #found
        try:
            include_dirs.append(self.build_clib)
        except:
            include_dirs = [self.build_clib]
        #now compile everything
        objects = self.compiler.compile(sources,
                                        output_dir=self.build_temp,
                                        macros=macros,
                                        include_dirs=include_dirs,
                                        debug=self.debug)
        return objects


    def _link_objs( self, lib_name, objects, build_info):
        #extra linking arguments
        extra_link_args = build_info.get('extra_link_args')
        if extra_link_args is None: extra_link_args = []

        #
        self._check_linker_for_darwin()

        #build shared library
        #self.compiler.set_executable( compiler_cxx = 'c++' )
        if _isDarwin():
            #darwin platform needs special treatment
            self._link_darwin_dylib(
                objects, lib_name,
                output_dir=self.build_clib,
                debug=self.debug,
                #build_temp=self.build_temp,
                extra_preargs=extra_link_args,
                target_lang = 'c++',
                libraries = build_info['libs'],
                library_dirs = build_info['libdirs'] 
                )
        else:
            self.compiler.link_shared_lib(
                objects, lib_name,
                output_dir=self.build_clib,
                debug=self.debug,
                #build_temp=self.build_temp,
                extra_preargs=extra_link_args,
                target_lang = 'c++',
                libraries = build_info['libs'],
                library_dirs = build_info['libdirs'] 
                )

        #build static library just in case of windows system
        import os
        if os.name == "nt" :
            # Now "link" the object files together into a static library.
            # (On Unix at least, this isn't really linking -- it just
            # builds an archive.  Whatever.)
            self.compiler.create_static_lib(objects, lib_name,
                                            output_dir=self.build_clib,
                                            debug=self.debug)

        return

    def _link_darwin_dylib(self, objects, output_libname,
                           output_dir=None,
                           libraries=None,
                           library_dirs=None,
                           runtime_library_dirs=None,
                           export_symbols=None,
                           debug=0,
                           extra_preargs=None,
                           extra_postargs=None,
                           build_temp=None,
                           target_lang=None):
        """link dynamic loading library for darwin platform. this is necessary
        because different from most unix system, .so and .dylib is different
        in Mac OS X."""
        compiler = self.compiler
        compiler.link(
            compiler.SHARED_LIBRARY, objects,
            compiler.library_filename(output_libname, lib_type='dylib'),
            output_dir,
            libraries, library_dirs, runtime_library_dirs,
            export_symbols, debug,
            extra_preargs, extra_postargs, build_temp, target_lang)


    def _check_linker_for_darwin(self):
        if _isDarwin():
            #for darwin platform, dylib should be built with -dynamiclib
            #and without -bundle
            self._saved_linker_so = self.compiler.linker_so[:]
            linker_so = self.compiler.linker_so
            if '-bundle' in linker_so:
                del ( linker_so[ linker_so.index('-bundle') ] )
                pass
            _append( linker_so, '-dynamiclib')
            _append( linker_so, '-single_module')
            return
        return

            
def _isDarwin():
    import sys
    return sys.platform[:6] == "darwin"


def _append( aList, element ):
    if not (element in aList): aList.append(element)



# version
__id__ = "$Id: build_clib.py 176 2006-10-15 16:35:31Z linjiao $"

# End of file
