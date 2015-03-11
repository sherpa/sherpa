# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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


import os
import sys
import IPython.ipapi
import IPython.iplib

try:
   if __ciao_ahelp_context__:
      pass
except:
   __ciao_ahelp_context__ = None
   
def ciao_exception_handler(self, type, value, traceback):
    if str(value).find('chips ') == 0:
       print >> sys.stderr, '%s' % (value,)
    else:
       print >> sys.stderr, '%s: %s' % (type.__name__, value)
IPython.ipapi.get().set_custom_exc((Exception,), ciao_exception_handler)

# If readline.clear_history() undefined, add a no-op
# IPython is supposed to do so itself, doesn't as of 0.8.0.
if hasattr(IPython.iplib.readline, "clear_history") is False:
   def clear_history():
      pass
   IPython.iplib.readline.clear_history = clear_history

def ahelp (*args):
   import os
   cmd = ''
   for elem in args:
      cmd = cmd + ' ' + elem
   if cmd.find ('-x') == -1:
      os.system ("ahelp -x '/py[.]*/' "+cmd)
   else:
      os.system ("ahelp "+cmd)


del IPython, ciao_exception_handler

def script(filename="sherpa.log", clobber=False):
   """
   script

   SYNOPSIS
      Save Sherpa session commands to file

   SYNTAX

   Arguments:
      filename  - script filename
                  default = 'sherpa.log'
      clobber   - clobber file flag
                  default = False

   Returns:
      None

   DESCRIPTION
      Save Sherpa commands from current session to a script file.

   SEE ALSO
      save
   """
   hist = _ip.IP.input_hist

   if os.path.isfile(filename) and not clobber:
      raise IOError("file '%s' exists and clobber is not set" % filename)

   fp = file(filename,"w")
   fp.writelines(hist)
   fp.close()
