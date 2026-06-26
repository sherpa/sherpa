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

.. versionchanged:: 4.19.0

   Functionality not used by Sherpa has been removed, including the
   setup process (supporting multiple options), the showFITSFile
   method, and removing the unused dataFunc argument to `xpaset`.  The
   code assumes that DS9 8.7 or later is in use (if an older version
   is used then some image displays may not be correct).  XPA
   communication now defaults to the "local" method unless the
   XPA_METHOD environment variable is set.

"""

from collections.abc import Mapping
import os
import shlex
import sys
import time
from typing import Any
import subprocess

import numpy as np

from sherpa.utils.err import RuntimeErr, TypeErr

__all__ = ["xpaget", "xpaset", "DS9Win"]


def _findUnixApp(appName: str) -> None:
    """Search PATH to find first directory that has the application.

    The call will raise a RuntimeErr if appName can not be found.

    Parameters
    ----------
    appName
       The application name

    """
    appPath = ''
    for path in os.environ['PATH'].split(':'):
        if os.access(path + '/' + appName, os.X_OK):
            appPath = path
            break

    if appPath == '' or not appPath.startswith("/"):
        raise RuntimeErr('notonpath', appName)


# If ds9 and the xpa tools are accessible (only xpaget is checked for)
# then things are fine. If not, error out.
#
try:
    _findUnixApp("ds9")
    _findUnixApp("xpaget")
except RuntimeErr as e:
    raise RuntimeErr('badwin', e) from e


_DefTemplate: str = "sherpa"

_OpenCheckInterval: float = 0.2  # seconds
_MaxOpenTime: float = 60.0  # seconds


def xpaget(cmd: str,
           template: str = _DefTemplate,
           method: str | None = None
           ) -> str:
    """Executes a simple xpaget command, returning the reply.

    Parameters
    ----------
    cmd
       The XPA command.
    template
       The target of the XPA call. It can be the ds9 window title,
       a string giving "host:port", or other supported forms.
    method
       The communication method (optional).

    Returns
    -------
    response
       The respose from DS9 to the query.

    """

    fullCmd = ['xpaget']
    if method is not None:
        fullCmd.extend(['-m', method])

    fullCmd.extend([template, cmd])

    with subprocess.Popen(args=fullCmd,
                          shell=False,
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE) as p:
        try:
            p.stdin.close()
            errMsg = p.stderr.read()
            if errMsg:
                errMsgStr = errMsg.decode()
                cmdStr = shlex.join(fullCmd)
                raise RuntimeErr('cmdfail', cmdStr, errMsgStr)

            return_value = p.stdout.read()
            return return_value.decode()

        finally:
            p.stdout.close()
            p.stderr.close()


def xpaset(cmd: str,
           data: str | bytes | None = None,
           template: str = _DefTemplate,
           method: str | None = None
           ) -> None:
    """Executes a single xpaset command.

    Parameters
    ----------
    cmd
       The XPA command.
    data
       Extra data to send via stdout (a trailing new-line character is
       added if needed).
    template
       The target of the XPA call. It can be the ds9 window title,
       a string giving "host:port", or other supported forms.
    method
       The communication method (optional).

    """

    fullCmd = ['xpaset']
    if method is not None:
        fullCmd.extend(['-m', method])

    if not data:
        fullCmd.append('-p')

    fullCmd.extend([template, cmd])
    with subprocess.Popen(args=fullCmd,
                          shell=False,
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT) as p:
        try:
            try:
                data = bytearray(data, "UTF-8")
            except TypeError:
                pass

            if data:
                p.stdin.write(data)
                if data[-1] != b'\n':
                    p.stdin.write(b'\n')
            p.stdin.close()
            reply = p.stdout.read()
            if reply:
                errMsgStr = reply.strip().decode()
                cmdStr = shlex.join(fullCmd)
                raise RuntimeErr('cmdfail', cmdStr, errMsgStr)


        finally:
            p.stdin.close()  # redundant
            p.stdout.close()


# What data types need to be converted before sending to DS9?
# DS9 8.7 supports
#
#     Integer types     Floating-point types
#     -------------     --------------------
#     unsigned 8-bit          32 bit
#              16-bit         64-bit
#     unsigned 16-bit
#              32-bit
#              64-bit
#
# There is no attempt to make this conversion match the capabilities
# of the DS9 being used.
#
_CnvDict = {
    np.int8: np.int16,
    np.uint16: np.int32,
    np.uint32: np.int64,
}

if hasattr(np, "uint64"):
    _CnvDict[np.uint64] = np.float64

_FloatTypes = (np.float32, np.float64)


def _formatOptions(kargs: Mapping[str, Any]) -> str:
    """Returns a string: "key1=val1,key2=val2,..."
    (where keyx and valx are string representations)
    """
    arglist = [f"{k}={str(v)}" for k,v in kargs.items()]
    return ','.join(arglist)


class DS9Win:
    """An object that talks to a particular window on ds9

    .. versionchanged:: 4.19.0
       The XPA communication method now defaults to "local" unless the
       XPA_METHOD environment variable is set when the class is
       created.

    Parameters
    ----------
    template
       window name (see ds9 docs for talking to a remote ds9)
    doOpen
       open ds9 using the desired template, if not already open; MacOS
       X warning: opening ds9 requires ds9 to be on your PATH; this
       may not be true by default; see the module documentation above
       for workarounds.

    """
    def __init__(self,
                 template: str = _DefTemplate,
                 doOpen: bool = True
                 ) -> None:
        self.template = str(template)

        # What communication method to use? CIAO defaults to using the
        # "local" method, so follow this (as there have been problems
        # on macOS with the default method of "inet"). However, this
        # is only done if the XPA_METHOD environment variable is not
        # set (note there is no check whether the variable is set to
        # anything sensible). The aim is to use the same method for
        # each communication (as this should not change once DS9 has
        # been started). The xpa_method field could be set to the
        # os.environ variable, but for now leave as is.
        #
        self.xpa_method = None if "XPA_METHOD" in os.environ else "local"

        self.alreadyOpen = self.isOpen()
        if doOpen:
            self.doOpen()

    def doOpen(self) -> None:
        """Open the ds9 window (if necessary).

        .. versionchanged:: 4.19.0
           The communication method is set to "local" unless the
           XPA_METHOD environment variable is set.

        """
        if self.isOpen():
            return

        args = ['ds9', '-title', self.template, '-port', "0"]
        if self.xpa_method is not None:
            args.extend(["-xpa", self.xpa_method])

        with subprocess.Popen(
                args=args,
                shell=False,
                cwd=None,
                close_fds=True,
                stdin=None,
                stdout=None,
                stderr=None) as p:

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
        """Is the DS9 window open and responding to XPA queries?"""
        try:
            _ = self.xpaget('mode')
            return True
        except RuntimeErr:
            return False

    def showArray(self, arr) -> None:
        """Display a 2-d or 3-d grayscale integer numarray arrays.

        3-d images are displayed as data cubes, meaning one can
        view a single z at a time or play through them as a movie,
        that sort of thing.

        .. versionchanged:: 4.19.0
           The datatype conversion used assumes that DS9 8.7 or later
           is in use.

        Parameters
        ----------
        arr
           An array (expected to be NumPy but need not be) that must
           be 2-d (axis order is y,x) or 3-d (axis order is z,y,x).

        Notes
        -----
        Complex data types cause an invalid. Other data types may be
        modified before sending to DS9.

        """

        if not hasattr(arr, "dtype") or not hasattr(arr, "astype"):
            arr = np.array(arr)

        if np.iscomplexobj(arr):
            raise TypeErr('nocomplex')

        ndim = arr.ndim
        if ndim not in (2, 3):
            raise RuntimeErr('only2d3d')

        # if necessary, convert array type
        cnvType = _CnvDict.get(arr.dtype.type)
        if cnvType:
            arr = arr.astype(cnvType)

        # determine byte order of array
        # First check if array endianness is not native--if
        # not, use the nonnative endianness
        # If the byteorder is native, then use the system
        # endianness
        if arr.dtype.byteorder == '>':
            arch = "bigendian"
        elif arr.dtype.byteorder == '<':
            arch = "littleendian"
        elif sys.byteorder == 'big':
            arch = "bigendian"
        else:
            arch = "littleendian"

        # compute bits/pix; ds9 uses negative values for floating values
        bitsPerPix = arr.itemsize * 8

        if arr.dtype.type in _FloatTypes:
            # array is float; use negative value
            bitsPerPix = -bitsPerPix

        # generate array info keywords; note that numarray
        # 2-d images are in order [y, x]
        # 3-d images are in order [z, y, x]
        arryDict = {}
        dimNames = ["z", "y", "x"][3 - ndim:]
        for axis, size in zip(dimNames, arr.shape):
            arryDict[f"{axis}dim"] = size

        arryDict["bitpix"] = bitsPerPix
        arryDict["arch"] = arch

        self.xpaset(
            cmd=f'array [{_formatOptions(arryDict)}]',
            data=arr.tobytes(),
        )

    def xpaget(self,
               cmd: str
               ) -> str:
        """Execute a simple xpaget command and return the reply.

        Parameters
        ----------
        cmd
           The XPA command.

        Returns
        -------
        response
           The respose from DS9 to the query.

        """
        return xpaget(
            cmd=cmd,
            template=self.template,
            method=self.xpa_method
        )

    def xpaset(self,
               cmd: str,
               data: str | bytes | None = None
               ) -> None:
        """Executes a simple xpaset command.

        Parameters
        ----------
        cmd
           The XPA command.
        data
           Extra data to send via stdout (a trailing new-line
           character is added if needed).

        """
        xpaset(
            cmd=cmd,
            data=data,
            template=self.template,
            method=self.xpa_method
        )
