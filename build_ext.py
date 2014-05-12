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


from distutils.command.build_ext import build_ext as _base, ListType, TupleType, customize_compiler


class build_ext(_base):

    def run(self):
        # If we were asked to build any C/C++ libraries, make sure that the
        # directory where we put them is in the library search path for
        # include
        if self.distribution.has_c_libraries():
            build_clib = self.get_finalized_command('build_clib')
            self.include_dirs.append(build_clib.build_clib)
            self.library_dirs.append(build_clib.build_clib)
            pass

        #_base.run(self)
        #return


        from distutils.ccompiler import new_compiler

        # 'self.extensions', as supplied by setup.py, is a list of
        # Extension instances.  See the documentation for Extension (in
        # distutils.extension) for details.
        #
        # For backwards compatibility with Distutils 0.8.2 and earlier, we
        # also allow the 'extensions' list to be a list of tuples:
        #    (ext_name, build_info)
        # where build_info is a dictionary containing everything that
        # Extension instances do except the name, with a few things being
        # differently named.  We convert these 2-tuples to Extension
        # instances as needed.

        if not self.extensions:
            return


        
        # If we were asked to build any C/C++ libraries, make sure that the
        # directory where we put them is in the library search path for
        # linking extensions.

        # it is not always good to blindly add all c libraries to link list
        # so we should comment out following codes
##         if self.distribution.has_c_libraries():
##             build_clib = self.get_finalized_command('build_clib')
##             self.libraries.extend(build_clib.get_library_names() or [])
##             self.library_dirs.append(build_clib.build_clib)


        # Setup the CCompiler object that we'll use to do all the
        # compiling and linking
        self.compiler = new_compiler(compiler=self.compiler,
                                     verbose=self.verbose,
                                     dry_run=self.dry_run,
                                     force=self.force)
        customize_compiler(self.compiler)

        # And make sure that any compile/link-related options (which might
        # come from the command-line or from the setup script) are set in
        # that CCompiler object -- that way, they automatically apply to
        # all compiling and linking done here.
        if self.include_dirs is not None:
            self.compiler.set_include_dirs(self.include_dirs)
        if self.define is not None:
            # 'define' option is a list of (name,value) tuples
            for (name,value) in self.define:
                self.compiler.define_macro(name, value)
        if self.undef is not None:
            for macro in self.undef:
                self.compiler.undefine_macro(macro)
        if self.libraries is not None:
            self.compiler.set_libraries(self.libraries)
        if self.library_dirs is not None:
            self.compiler.set_library_dirs(self.library_dirs)
        if self.rpath is not None:
            self.compiler.set_runtime_library_dirs(self.rpath)
        if self.link_objects is not None:
            self.compiler.set_link_objects(self.link_objects)

        # Now actually compile and link everything.
        self.build_extensions()

        return


# version
__id__ = "$Id: build_ext.py 146 2005-08-09 17:03:59Z linjiao $"

# End of file
