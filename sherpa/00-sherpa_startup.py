# -*- Mode: Shell-Script -*-  Not really, but shows comments correctly
# 
#  Copyright (C) 2012  Smithsonian Astrophysical Observatory
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

#

import os
import sys
import copy

# Take any instance of '' out of sys.path, since that instructs import
# statements to first look in the current working directory for modules
# to load.  So, scripts named "sherpa.py", "group.py", and so on, will
# override system modules of the same name, that we always need to load.
# This method is more robust because it takes out '' from anywhere in
# sys.path

# Save the old sys.path and restore it after the imports of CIAO modules
# are complete, so that the user can later import user modules from the
# current working directory.  SMD 10/26/12

old_sys_path = copy.deepcopy(sys.path)
while ('' in sys.path):
    sys.path.remove('')

from pychips import *
from pycrates import *
from pychips.hlui import *
from pycrates.hlui import *
from ipython_cxc import *
from ahelp_interface import *
from sherpa.astro.ui import *

import ipython_cxc
import ahelp_interface

def _initialize_sherpa_app_():
  sherpa_version = 0
  chips_version = 0
  chips_version_str = "0"
  crates_version = 0
  crates_version_str = '0.0.0'
  sherpa_path = ''

  try:
    from pycrates import __version__ as crates_version
  except:
    pass

  try:
    from pychips import __version__ as chips_version
  except:
    pass

  try:
    from pycrates import __versionstr__ as crates_version_str
  except:
    pass

  try:
    from pychips import __versionstr__ as chips_version_str
  except:
    pass

  try:
    from sherpa import __version__ as sherpa_version
  except:
    pass

  try:
    from sherpa import __file__ as sherpa_path
    sherpa_path = os.path.dirname(sherpa_path)
  except:
    pass
  site_path = sherpa_path.replace('/sherpa','',-1)

  if sherpa_version > 0 and crates_version < 40501:
     print "Warning: Importing CRATES version {0}; CRATES version 4.5.1 is now available, consider upgrading".format (crates_version_str)

  if sherpa_version > 0 and chips_version < 40501:
     print "Warning: Importing ChIPS version {0}; ChIPS version 4.5.1 is now available, consider upgrading".format (chips_version_str)



_initialize_sherpa_app_()

ahelp_interface.__ciao_ahelp_context__ = 'py.sherpa'
ipython_cxc.init ("sherpa")

set_preference_autoload(True)

sys.path = old_sys_path

# Replace native help text, for Sherpa HLUI functions, with pointer to
# ahelp text for these functions.
for func in get_functions(): 
  try:
    eval(func).__doc__="Type ahelp("+func+") for more details."
  except:
    pass
