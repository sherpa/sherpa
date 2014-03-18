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



def _copy_header_trees( self, headerTrees, dest ):
    '''
    headerTrees must be a list of 3-tuples of  (name, sourcedir, exts):
    
      name: header package name
      sourcedir: header package source directory
      exts: header extensions
      
    dest: destination root directory
    '''
    from os.path import join
    for packageName, sourceDir, exts in headerTrees:
        _copy_header_tree( self, sourceDir, exts, join(dest, packageName) )
        pass
    return



def _copy_header_tree( self, sourceDir, exts, destDir ):
    '''copy a header tree
    '''
    from os.path import join
    packageDest = destDir

    for header in _headers(sourceDir, exts):
        remainder = _getPathRemainder( header, sourceDir )
        remainderDir = _getDirectory( remainder )
        destdir = join( packageDest, remainderDir )
        #print header, destdir
        self.mkpath(destdir)
        (out, _) = self.copy_file(header, destdir)
        self.outfiles.append(out)
        pass
    pass


def _getPathRemainder( path, directory ):
    """remove directory from path
     _getPathRemainder( 'a/b/c/1.cc', 'a/b' ) --> 'c/1.cc'
    """
    if not path.startswith( directory ):
        raise 'I am confused. Is %s in directory %s?'%(path, directory)
    path = path.replace( directory, '' )
    from os.path import sep
    #remove leading '/' if necessary
    if path.startswith(sep): path = path[1:]
    return path

def _getDirectory( path2file ):
    """remove the file from the path
     _getDirectory( 'a/b/1.cc' ) --> 'a/b'
    """
    from os.path import split
    return split(path2file)[0]


#extend the member function of Command
#now every command class has these two member functions
from distutils.cmd import Command
Command._copy_header_trees = _copy_header_trees
Command._copy_header_tree = _copy_header_tree



def _headers( path, exts = ['.h', '*.icc', '*.hh'] ):
    #return all header files under path
    from _fileLister import _recursiveListByExts
    return _recursiveListByExts( path, exts )


# version
__id__ = "$Id: _addons.py 176 2006-10-15 16:35:31Z linjiao $"

# End of file 
