#
#  Copyright (C) 2006-2010, 2016-2021, 2025-2026
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

"""Interface for viewing images with the ds9 image viewer.

Loosely based on XPA, by Andrew Williams, with the original code by
ROwen 2004-2005 and then from the Sherpa team from 2006. This code has
been simplified to only support the features that Sherpa needs.

"""

from collections.abc import Mapping
import os
import sys
import time
from typing import Any
import warnings
import subprocess

import numpy as np

from sherpa.utils.err import RuntimeErr, TypeErr

__all__ = ["setup", "xpaget", "xpaset", "DS9Win"]


def _findUnixApp(appName: str) -> str:
    """Search PATH to find first directory that has the application.
    Return the path if found.
    Raise RuntimeError if not found.
    """
    appPath = ''
    for path in os.environ['PATH'].split(':'):
        if os.access(path + '/' + appName, os.X_OK):
            appPath = path
            break

    if appPath == '' or not appPath.startswith("/"):
        raise RuntimeErr('notonpath', appName)

    return appPath


def _findDS9AndXPA() -> tuple[str, str]:
    """Locate ds9 and xpa, and add to PATH if not already there.

    Returns:
    - ds9Dir        directory containing ds9 executable
    - xpaDir        directory containing xpaget and (presumably)
                            the other xpa executables

    Raise RuntimeError if ds9 or xpa are not found.
    """
    ds9Dir = _findUnixApp("ds9")
    xpaDir = _findUnixApp("xpaget")

    return (ds9Dir, xpaDir)


def setup(doRaise: bool = True,
          debug: bool = False
          ) -> str | None:
    """Search for xpa and ds9 and set globals accordingly.
    Return None if all is well, else return an error string.
    The return value is also saved in global variable _SetupError.

    Sets globals:
    - _SetupError        same value as returned
    - _Popen                subprocess.Popen, if ds9 and xpa found,
                                    else a variant that searches for ds9 and xpa
                                    first and then runs subprocess.Popen if found
                                    else raises an exception
                                    This permits the user to install ds9 and xpa
                                    and use this module without reloading it
    """
    global _SetupError, _Popen, _ex
    _SetupError = None
    try:
        ds9Dir, xpaDir = _findDS9AndXPA()
        if debug:
            print(f"ds9Dir={repr(ds9Dir)}\npaDir={repr(xpaDir)}")
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e:
        _ex = e
        _SetupError = f"DS9Win unusable: {e}"
        ds9Dir = xpaDir = None

    if _SetupError:
        class _Popen(subprocess.Popen):
            def __init__(self, *args, **kargs):
                setup(doRaise=True)
                super().__init__(*args, **kargs)

        if doRaise:
            raise RuntimeErr('badwin', _ex)

    else:
        _Popen = subprocess.Popen

    return _SetupError


errStr = setup(doRaise=True, debug=False)
if errStr:
    warnings.warn(errStr)

_ArrayKeys = ("dim", "dims", "xdim", "ydim", "zdim", "bitpix", "skip", "arch")
_DefTemplate = "sherpa"

_OpenCheckInterval = 0.2  # seconds
_MaxOpenTime = 60.0  # seconds


def xpaget(cmd: str,  # Do not try to type the "this can be a list" version
           template: str = _DefTemplate,
           doRaise: bool = True
           ) -> str:
    """Executes a simple xpaget command:
            xpaget <template> <cmd>
    returning the reply.

    Inputs:
    - cmd                command to execute; may be a string or a list
    - template        xpa template; can be the ds9 window title
                            (as specified in the -title command-line option)
                            host:port, etc.
    - doRaise        if True, raise RuntimeError if there is a communications error,
                            else issue a UserWarning warning

    Raises RuntimeError or issues a warning (depending on doRaise)
    if anything is written to stderr.
    """
    # Would be better to make a sequence rather than have to quote arguments
    fullCmd = f'xpaget {template} "{cmd}"'

    with _Popen(args=fullCmd,
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE) as p:
        try:
            p.stdin.close()
            errMsg = p.stderr.read()
            if errMsg:
                errMsgStr = errMsg.decode()
                if doRaise:
                    raise RuntimeErr('cmdfail', fullCmd, errMsgStr)

                fullErrMsg = f"{repr(fullCmd)} failed: {errMsgStr}"
                warnings.warn(fullErrMsg)

            return_value = p.stdout.read()
            return return_value.decode()

        finally:
            p.stdout.close()
            p.stderr.close()


def xpaset(cmd: str,
           data: str | bytes | None = None,
           dataFunc=None,
           template: str = _DefTemplate,
           doRaise: bool = True
           ) -> None:
    """Executes a simple xpaset command:
            xpaset -p <template> <cmd>
    or else feeds data to:
            xpaset <template> <cmd>

    The command must not return any output for normal completion.

    Inputs:
    - cmd                command to execute
    - data                data to write to xpaset's stdin; ignored if dataFunc specified.
                            If data[-1] is not \n then a final \n is appended.
    - dataFunc        a function that takes one argument, a file-like object,
                            and writes data to that file. If specified, data is ignored.
                            Warning: if a final \n is needed, dataFunc must supply it.
    - template        xpa template; can be the ds9 window title
                            (as specified in the -title command-line option)
                            host:port, etc.
    - doRaise        if True, raise RuntimeError if there is a communications error,
                            else issue a UserWarning warning

    Raises RuntimeError or issues a warning (depending on doRaise)
    if anything is written to stdout or stderr.
    """
    # Would be better to make a sequence rather than have to quote arguments
    if data or dataFunc:
        fullCmd = f'xpaset {template} "{cmd}"'
    else:
        fullCmd = f'xpaset -p {template} "{cmd}"'

    with _Popen(args=fullCmd,
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT) as p:
        try:
            try:
                data = bytearray(data, "UTF-8")
            except Exception:
                pass

            if data:
                p.stdin.write(data)
                if data[-1] != b'\n':
                    p.stdin.write(b'\n')
            p.stdin.close()
            reply = p.stdout.read()
            if reply:
                errMsgStr = reply.strip().decode()
                if doRaise:
                    raise RuntimeErr('cmdfail', fullCmd, errMsgStr)

                fullErrMsg = f"{repr(fullCmd)} failed: {errMsgStr}"
                warnings.warn(fullErrMsg)

        finally:
            p.stdin.close()  # redundant
            p.stdout.close()


def _computeCnvDict():
    """Compute array type conversion dict.
    Each item is: unsupported type: type to which to convert.

    ds9 supports UInt8, Int16, Int32, Float32 and Float64.
    """

    cnvDict = {
        np.int8: np.int16,
        np.uint16: np.int32,
        np.uint32: np.float32,  # ds9 can't handle 64 bit integer data
        np.int64: np.float64,
    }

    # TODO: should this check for 'uint64' since 'uint64=' is not a
    #       valid attribute name
    if hasattr(np, "uint64="):
        cnvDict[np.uint64] = np.float64

    return cnvDict


_CnvDict = _computeCnvDict()
_FloatTypes = (np.float32, np.float64)
_ComplexTypes = (np.complex64, np.complex128)


def _formatOptions(kargs: Mapping[str, Any]) -> str:
    """Returns a string: "key1=val1,key2=val2,..."
    (where keyx and valx are string representations)
    """
    arglist = [f"{k}={str(v)}" for k,v in kargs.items()]
    return ','.join(arglist)


def _splitDict(inDict,
               keys
               ):
    """Splits a dictionary into two parts:
    - outDict contains any keys listed in "keys";
      this is returned by the function
    - inDict has those keys removed (this is the dictionary passed in;
      it is modified by this call)
    """
    outDict = {}
    for key in keys:
        if key in inDict:
            outDict[key] = inDict.pop(key)
    return outDict


class DS9Win:
    """An object that talks to a particular window on ds9

    Inputs:
    - template:        window name (see ds9 docs for talking to a remote ds9)
    - doOpen: open ds9 using the desired template, if not already open;
                    MacOS X warning: opening ds9 requires ds9 to be on your PATH;
                    this may not be true by default;
                    see the module documentation above for workarounds.
    - doRaise        if True, raise RuntimeError if there is a communications error,
                    else issue a UserWarning warning.
                    Note: doOpen always raises RuntimeError on failure!
    """
    def __init__(self,
                 template: str = _DefTemplate,
                 doOpen: bool = True,
                 doRaise: bool = True
                 ) -> None:
        self.template = str(template)
        self.doRaise = bool(doRaise)
        self.alreadyOpen = self.isOpen()
        if doOpen:
            self.doOpen()

    def doOpen(self) -> None:
        """Open the ds9 window (if necessary).

        Raise OSError or RuntimeError on failure, even if doRaise is False.
        """
        if self.isOpen():
            return

        # We want to fork ds9. This is possible with os.fork, but
        # it doesn't work on Windows. At present Sherpa does not
        # run on Windows, so it is not a serious problem, but it is
        # not clear if it is an acceptable, or sensible, option.
        #
        p = _Popen(
            args=('ds9', '-title', self.template, '-port', "0"),
            cwd=None,
            close_fds=True, stdin=None, stdout=None, stderr=None
        )

        startTime = time.time()
        while True:
            time.sleep(_OpenCheckInterval)
            if self.isOpen():
                # Trick to stop a ResourceWarning warning to be created when
                # running sherpa/tests/test_image.py
                #
                # Adapted from https://hg.python.org/cpython/rev/72946937536e
                p.returncode = 0
                return
            if time.time() - startTime > _MaxOpenTime:
                raise RuntimeErr('nowin', self.template)

    def isOpen(self) -> bool:
        """Return True if this ds9 window is open
        and available for communication, False otherwise.
        """
        try:
            xpaget('mode', template=self.template, doRaise=True)
            return True
        except RuntimeErr:
            return False

    def showArray(self,
                  arr,
                  **kargs) -> None:
        """Display a 2-d or 3-d grayscale integer numarray arrays.
        3-d images are displayed as data cubes, meaning one can
        view a single z at a time or play through them as a movie,
        that sort of thing.

        Inputs:
        - arr: a numarray array; must be 2-d or 3-d:
                2-d arrays have index order (y, x)
                3-d arrays are loaded as a data cube index order (z, y, x)
        kargs: see Extra Keyword Arguments in the module doc string for information.
        Keywords that specify array info (see doc for showBinFile for the list)
        are ignored, because array info is determined from the array itself.

        Data types:
        - UInt8, Int16, Int32 and floating point types sent unmodified.
        - All other integer types are converted before transmission.
        - Complex types are rejected.

        Raises ValueError if arr's elements are not some kind of integer.
        Raises RuntimeError if ds9 is not running or returns an error message.
        """
        if not hasattr(arr, "dtype") or not hasattr(arr, "astype"):
            arr = np.array(arr)

        if np.iscomplexobj(arr):
            raise TypeErr('nocomplex')

        ndim = arr.ndim
        if ndim not in (2, 3):
            raise RuntimeErr('only2d3d')

        dimNames = ["z", "y", "x"][3 - ndim:]

        # if necessary, convert array type
        cnvType = _CnvDict.get(arr.dtype.type)
        if cnvType:
            # print "converting array from %s to %s" % (arr.type(), cnvType)
            arr = arr.astype(cnvType)

        # determine byte order of array
        # First check if array endianness is not native--if
        # not, use the nonnative endianness
        # If the byteorder is native, then use the system
        # endianness
        if arr.dtype.byteorder == '>':
            isBigendian = True
        elif arr.dtype.byteorder == '<':
            isBigendian = False
        else:
            if sys.byteorder == 'big':
                isBigendian = True
            else:
                isBigendian = False

        # compute bits/pix; ds9 uses negative values for floating values
        bitsPerPix = arr.itemsize * 8

        # if np.issubclass_(arr.dtype.type, float):
        if arr.dtype.type in _FloatTypes:
            # array is float; use negative value
            bitsPerPix = -bitsPerPix

        # remove array info keywords from kargs; we compute all that
        _splitDict(kargs, _ArrayKeys)

        # generate array info keywords; note that numarray
        # 2-d images are in order [y, x]
        # 3-d images are in order [z, y, x]
        arryDict = {}
        for axis, size in zip(dimNames, arr.shape):
            arryDict[f"{axis}dim"] = size

        arryDict["bitpix"] = bitsPerPix
        if isBigendian:
            arryDict["arch"] = 'bigendian'
        else:
            arryDict["arch"] = 'littleendian'

        self.xpaset(
            cmd=f'array [{_formatOptions(arryDict)}]',
            data=arr.tobytes(),
        )

        for keyValue in kargs.items():
            self.xpaset(cmd=' '.join(keyValue))

    def xpaget(self,
               cmd: str
               ) -> str:
        """Execute a simple xpaget command and return the reply.

        Inputs:
        - cmd                command to execute

        Raises RuntimeError if anything is written to stderr.
        """
        return xpaget(
            cmd=cmd,
            template=self.template,
            doRaise=self.doRaise,
        )

    def xpaset(self,
               cmd: str,
               data: str | bytes | None = None,
               dataFunc=None
               ) -> None:
        """Executes a simple xpaset command.

        The command must not return any output for normal completion.

        Inputs:
        - cmd                command to execute
        - data                data to write to xpaset's stdin; ignored if dataFunc specified
        - dataFunc        a function that takes one argument, a file-like object,
                                and writes data to that file. If specified, data is ignored.

        Raises RuntimeError if anything is written to stdout or stderr.
        """
        xpaset(
            cmd=cmd,
            data=data,
            dataFunc=dataFunc,
            template=self.template,
            doRaise=self.doRaise,
        )


if __name__ == "__main__":
    errStr = setup(doRaise=True, debug=True)
    if errStr:
        print(errStr)
    else:
        ds9Win = DS9Win("Test")
