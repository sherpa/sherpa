#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

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
