#
#  Copyright (C) 2015  Smithsonian Astrophysical Observatory
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

"""Test the functions in sherpa.astro.ui.serialize, also testing the
corresponding functions in sherpa.astro.ui.utils.
"""

import re
import StringIO
import tempfile
import unittest

from sherpa.utils import SherpaTest, SherpaTestCase, has_package_from_list, \
    test_data_missing
from sherpa.astro import ui
# from sherpa.astro.ui import serialize

import logging
logger = logging.getLogger('sherpa')

has_xspec = has_package_from_list("sherpa.astro.xspec")

# The tests can either check that the output ASCII is identical
# to a canonical form, or try to execute the saved file and
# check that the restored version matches the input version
# for selected metrics (e.g. data values, source expression,
# parameter values). There's merit in having both (since the
# latter is also a check that the saved form catches all
# relevant information).
#
# It is likely that the first case - comparing against a canonical
# form - is going to be limited by numerical accuracy. If the
# tests against a canonical form are to be used then it probably
# makes sense to store the output in external files, and read in
# the content, rather than write them directly to this file (although
# then have to deal with multiple versions depending on the modules
# that are available).

# A representation of the default Sherpa state
_canonical_empty = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets



######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("verbose", 0)
set_method_opt("factor", 100.0)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("xtol", 1.19209289551e-07)
set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("ftol", 1.19209289551e-07)


######### Set Model Components and Parameters



######### Set Source, Pileup and Background Models


"""

# Change a few settings for statistic/method
_canonical_empty_stats = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets



######### Set Statistic

set_stat("leastsq")


######### Set Fitting Method

set_method("neldermead")

set_method_opt("iquad", 1)
set_method_opt("initsimplex", 0)
set_method_opt("verbose", 1)
set_method_opt("step", None)
set_method_opt("finalsimplex", 9)
set_method_opt("maxfev", 5000)
set_method_opt("ftol", 1.19209289551e-07)


######### Set Model Components and Parameters



######### Set Source, Pileup and Background Models


"""

# use @@ as a marker to indicate that the location of the data directory
# needs to be inserted into the text before testing.
#
_canonical_pha_basic = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha(1, "@@/threads/pha_intro/3c273.pi")

######### Set Image Coordinates

if get_data(1).grouping is not None and not get_data(1).grouped:
    ######### Group Data
    group(1)

######### Data Spectral Responses

load_arf(1, "@@/threads/pha_intro/3c273.arf", resp_id=1)
load_rmf(1, "@@/threads/pha_intro/3c273.rmf", resp_id=1)

######### Load Background Data Sets

load_bkg(1, "@@/threads/pha_intro/3c273_bg.pi", bkg_id=1)
if get_bkg(1, 1).grouping is not None and not get_bkg(1, 1).grouped:
    ######### Group Background
    group(1, 1)

######### Background Spectral Responses

load_arf(1, "@@/threads/pha_intro/3c273.arf", resp_id=1, bkg_id=1)
load_rmf(1, "@@/threads/pha_intro/3c273.rmf", resp_id=1, bkg_id=1)

######### Set Energy or Wave Units

set_analysis(1, 'energy', "rate", 0)
if not get_data(1).subtracted:
    ######### Subtract Background Data
    subtract(1)

######### Filter Data

notice_id(1, "0.518300011754:8.219800233841")


######### Set Statistic

set_stat("chi2datavar")


######### Set Fitting Method

set_method("levmar")

set_method_opt("verbose", 0)
set_method_opt("factor", 100.0)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("xtol", 1.19209289551e-07)
set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("ftol", 1.19209289551e-07)


######### Set Model Components and Parameters

create_model_component("xsapec", "src")
src.integrate = True

src.redshift.default_val = 0.0
src.redshift.default_min = -0.999
src.redshift.default_max = 10.0
src.redshift.val     = 0.0
src.redshift.min     = -0.999
src.redshift.max     = 10.0
src.redshift.units   = ""
src.redshift.frozen  = True

src.Abundanc.default_val = 1.0
src.Abundanc.default_min = 0.0
src.Abundanc.default_max = 5.0
src.Abundanc.val     = 1.0
src.Abundanc.min     = 0.0
src.Abundanc.max     = 5.0
src.Abundanc.units   = ""
src.Abundanc.frozen  = True

src.kT.default_val = 1.0
src.kT.default_min = 0.0080000000000000002
src.kT.default_max = 64.0
src.kT.val     = 1.0
src.kT.min     = 0.0080000000000000002
src.kT.max     = 64.0
src.kT.units   = "keV"
src.kT.frozen  = False

src.norm.default_val = 1.0
src.norm.default_min = 0.0
src.norm.default_max = 9.9999999999999998e+23
src.norm.val     = 1.0
src.norm.min     = 0.0
src.norm.max     = 9.9999999999999998e+23
src.norm.units   = ""
src.norm.frozen  = False

create_model_component("powlaw1d", "pl")
pl.integrate = True

pl.ampl.default_val = 1.0
pl.ampl.default_min = 0.0
pl.ampl.default_max = 3.4028234663852886e+38
pl.ampl.val     = 1.0
pl.ampl.min     = 0.0
pl.ampl.max     = 3.4028234663852886e+38
pl.ampl.units   = ""
pl.ampl.frozen  = False

pl.ref.default_val = 1.0
pl.ref.default_min = -3.4028234663852886e+38
pl.ref.default_max = 3.4028234663852886e+38
pl.ref.val     = 1.0
pl.ref.min     = -3.4028234663852886e+38
pl.ref.max     = 3.4028234663852886e+38
pl.ref.units   = ""
pl.ref.frozen  = True

pl.gamma.default_val = 1.0
pl.gamma.default_min = -10.0
pl.gamma.default_max = 10.0
pl.gamma.val     = 1.0
pl.gamma.min     = -10.0
pl.gamma.max     = 10.0
pl.gamma.units   = ""
pl.gamma.frozen  = False

create_model_component("xsphabs", "gal")
gal.integrate = True

gal.nH.default_val = 1.0
gal.nH.default_min = 0.0
gal.nH.default_max = 100000.0
gal.nH.val     = 1.0
gal.nH.min     = 0.0
gal.nH.max     = 100000.0
gal.nH.units   = "10^22 atoms / cm^2"
gal.nH.frozen  = False



######### Set Source, Pileup and Background Models

set_source(1, (xsphabs.gal * (powlaw1d.pl + xsapec.src)))




"""

_canonical_pha_grouped = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha("grp", "@@/threads/pha_intro/3c273.pi")

######### Set Image Coordinates


######### Data grouping flags

set_grouping("grp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16))

######### Data quality flags

set_quality("grp", val=numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], numpy.int16))
if get_data("grp").grouping is not None and not get_data("grp").grouped:
    ######### Group Data
    group("grp")

######### Data Spectral Responses

load_arf("grp", "@@/threads/pha_intro/3c273.arf", resp_id=1)
load_rmf("grp", "@@/threads/pha_intro/3c273.rmf", resp_id=1)

######### Load Background Data Sets

load_bkg("grp", "@@/threads/pha_intro/3c273_bg.pi", bkg_id=1)

######### Background grouping flags

set_grouping("grp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)

######### Background quality flags

set_quality("grp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)
if get_bkg("grp", 1).grouping is not None and not get_bkg("grp", 1).grouped:
    ######### Group Background
    group("grp", 1)

######### Background Spectral Responses

load_arf("grp", "@@/threads/pha_intro/3c273.arf", resp_id=1, bkg_id=1)
load_rmf("grp", "@@/threads/pha_intro/3c273.rmf", resp_id=1, bkg_id=1)

######### Set Energy or Wave Units

set_analysis("grp", 'energy', "rate", 0)
if not get_data("grp").subtracted:
    ######### Subtract Background Data
    subtract("grp")

######### Filter Data

notice_id("grp", "0.489099994302:6.131999969482")


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("verbose", 0)
set_method_opt("factor", 100.0)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("xtol", 1.19209289551e-07)
set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("ftol", 1.19209289551e-07)


######### Set Model Components and Parameters

create_model_component("powlaw1d", "gpl")
gpl.integrate = True

gpl.ampl.default_val = 1.0
gpl.ampl.default_min = 0.0
gpl.ampl.default_max = 3.4028234663852886e+38
gpl.ampl.val     = 1.0
gpl.ampl.min     = 0.0
gpl.ampl.max     = 3.4028234663852886e+38
gpl.ampl.units   = ""
gpl.ampl.frozen  = False

gpl.ref.default_val = 1.0
gpl.ref.default_min = -3.4028234663852886e+38
gpl.ref.default_max = 3.4028234663852886e+38
gpl.ref.val     = 1.0
gpl.ref.min     = -3.4028234663852886e+38
gpl.ref.max     = 3.4028234663852886e+38
gpl.ref.units   = ""
gpl.ref.frozen  = True

gpl.gamma.default_val = 1.0
gpl.gamma.default_min = -10.0
gpl.gamma.default_max = 10.0
gpl.gamma.val     = 1.0
gpl.gamma.min     = -10.0
gpl.gamma.max     = 5.0
gpl.gamma.units   = ""
gpl.gamma.frozen  = False

create_model_component("xsphabs", "ggal")
ggal.integrate = True

ggal.nH.default_val = 2.0
ggal.nH.default_min = 0.0
ggal.nH.default_max = 100000.0
ggal.nH.val     = 2.0
ggal.nH.min     = 0.0
ggal.nH.max     = 100000.0
ggal.nH.units   = "10^22 atoms / cm^2"
ggal.nH.frozen  = True



######### Set Source, Pileup and Background Models

set_source("grp", (xsphabs.ggal * powlaw1d.gpl))




"""


# The serialization of the user model is not ideal, but check that
# we return something useful.
#
# This is also a good test of serializing data sets that was
# created by load_arrays (currently not very well).
#
_canonical_usermodel = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_arrays(3,
            [1.0, 12.2, 2.0, 14.0],
            [4, 8, 12, 4],
            Data1D)

######### Set Image Coordinates


######### Data Spectral Responses


######### Load Background Data Sets


######### Set Energy or Wave Units


######### Filter Data

notice_id(3, "1.0000:14.0000")


######### Set Statistic

set_stat("cash")


######### Set Fitting Method

set_method("neldermead")

set_method_opt("iquad", 1)
set_method_opt("initsimplex", 0)
set_method_opt("verbose", 0)
set_method_opt("step", None)
set_method_opt("finalsimplex", 9)
set_method_opt("maxfev", None)
set_method_opt("ftol", 1.19209289551e-07)


######### Set Model Components and Parameters

create_model_component("sin", "sin_model")
sin_model.integrate = True

sin_model.ampl.default_val = 1.0
sin_model.ampl.default_min = 1.0000000000000001e-05
sin_model.ampl.default_max = 3.4028234663852886e+38
sin_model.ampl.val     = 1.0
sin_model.ampl.min     = 1.0000000000000001e-05
sin_model.ampl.max     = 3.4028234663852886e+38
sin_model.ampl.units   = ""
sin_model.ampl.frozen  = False

sin_model.period.default_val = 1.0
sin_model.period.default_min = 1e-10
sin_model.period.default_max = 10.0
sin_model.period.val     = 1.0
sin_model.period.min     = 1e-10
sin_model.period.max     = 10.0
sin_model.period.units   = ""
sin_model.period.frozen  = False

sin_model.offset.default_val = 0.0
sin_model.offset.default_min = 0.0
sin_model.offset.default_max = 3.4028234663852886e+38
sin_model.offset.val     = 0.0
sin_model.offset.min     = 0.0
sin_model.offset.max     = 3.4028234663852886e+38
sin_model.offset.units   = ""
sin_model.offset.frozen  = False

print("Found user model 'mymodel'; please check it is saved correctly.")
def mymodel_func(pars, x, xhi=None):
    return pars[0] + pars[1] * x

load_user_model(mymodel_func, "mymodel")
add_user_model("mymodel",
               parnames=['c', 'm'],
               parvals=[2.0, 0.5],
               parmins=[-10.0, 0.0],
               parmaxs=[10.0, 5.5],
               parunits=['m', ''],
               parfrozen=[False, True]
               )

mymodel.integrate = True

mymodel.c.default_val = 2.0
mymodel.c.default_min = -3.4028234663852886e+38
mymodel.c.default_max = 3.4028234663852886e+38
mymodel.c.val     = 2.0
mymodel.c.min     = -10.0
mymodel.c.max     = 10.0
mymodel.c.units   = "m"
mymodel.c.frozen  = False

mymodel.m.default_val = 0.5
mymodel.m.default_min = -3.4028234663852886e+38
mymodel.m.default_max = 3.4028234663852886e+38
mymodel.m.val     = 0.5
mymodel.m.min     = 0.0
mymodel.m.max     = 5.5
mymodel.m.units   = ""
mymodel.m.frozen  = True



######### Set Source, Pileup and Background Models

set_full_model(3, (sin.sin_model + usermodel.mymodel))


"""

if has_xspec:
    _canonical_xspec_extra = """######### XSPEC Module Settings

set_xschatter(0)
set_xsabund("angr")
set_xscosmo(70, 0, 0.73)
set_xsxsect("bcmc")
"""

    _canonical_empty += _canonical_xspec_extra
    _canonical_empty_stats += _canonical_xspec_extra
    _canonical_pha_basic += _canonical_xspec_extra
    _canonical_pha_grouped += _canonical_xspec_extra
    _canonical_usermodel += _canonical_xspec_extra


def _dump_lines(cts):
    """Dump outcts, an array of strings, to stdout, with line numbering"""
    print("***\n")
    for i, l in enumerate(cts):
        print("{:02d} {}".format(i, l))
    print("***\n")


class test_ui(SherpaTestCase):

    def setUp(self):
        ui.clean()  # I thought the test harness did this anyway
        self._old_logging_level = logger.level
        logger.setLevel(logging.ERROR)
        if has_xspec:
            # try and avoid messages from setting XSPEC values, but
            # it does not seem to work
            self._xspec_chatter = ui.get_xschatter()
            ui.set_xschatter(0)

    def tearDown(self):
        logger.setLevel(self._old_logging_level)
        if has_xspec:
            ui.set_xschatter(self._xspec_chatter)

        ui.clean()

    def _add_datadir_path(self, output):
        """Replace any @@ characters by the value of self.datadir,
        making sure that the replacement text does not end in a /."""

        dname = self.datadir
        if dname.endswith('/'):
            dname = dname[:-1]

        return re.sub('@@', dname, output, count=0)

    def _compare_lines(self, expected, got):
        """Check that each line in got matches that in
        expected. This is to provide a more-readable
        set of error messages (may prefer a diff-style
        analysis).
        """

        elines = expected.split('\n')
        glines = got.split('\n')

        # _dump_lines(elines)
        # _dump_lines(glines)

        for e, g in zip(elines, glines):
            self.assertEqual(e, g)

        # Do the line length after checking for the file
        # contents as it is easier to see what the difference
        # is this way around, since a difference in the
        # number of lines is often not very informative.
        self.assertEqual(len(elines), len(glines))

    def _compare(self, expected):
        """Run save_all and check the output (saved to a
        StringIO object) to the string value expected.
        """
        output = StringIO.StringIO()
        ui.save_all(outfh=output)
        output = output.getvalue()
        self._compare_lines(expected, output)

    def _restore(self):
        """Run save_all then call clean and try to restore
        the Sherpa state from the saved file. Will raise
        a test failure if there was an error when
        executing the save file.
        """

        output = StringIO.StringIO()
        ui.save_all(outfh=output)
        output = output.getvalue()
        ui.clean()
        try:
            exec(output)
            success = True
            e = "no exception"
        except Exception as e:
            success = False

        self.assertTrue(success, msg="exception={}".format(e))

    def _setup_pha_basic(self):
        """Load up a PHA file and make "simple" changes to the
        Sherpa state. Returns the name of the file that is
        loaded and the canonical output.
        """

        ui.clean()
        fname = self.make_path('threads', 'pha_intro', '3c273.pi')
        ui.load_pha(1, fname)
        ui.subtract()
        ui.set_stat('chi2datavar')
        ui.notice(0.5, 7)
        ui.set_source(ui.xsphabs.gal * (ui.powlaw1d.pl +
                                        ui.xsapec.src))
        return fname, self._add_datadir_path(_canonical_pha_basic)

    def _setup_pha_grouped(self):
        """Add in grouping and a few different choices.

        Returns the name of the file that is
        loaded, the new grouping and quality arrays,
        and the canonical output.
        """

        ui.clean()
        fname = self.make_path('threads', 'pha_intro', '3c273.pi')
        ui.load_pha('grp', fname)
        channels = ui.get_data('grp').channel

        exclude = (channels < 20) | (channels > 800)
        qual = exclude * 1

        ui.subtract('grp')
        ui.group_counts('grp', 10, tabStops=exclude)
        ui.set_quality('grp', exclude)

        grp = ui.get_data('grp').grouping

        ui.set_stat('chi2gehrels')
        ui.notice_id('grp', 0.5, 6)
        ui.set_source('grp', ui.xsphabs.ggal * ui.powlaw1d.gpl)
        gpl.gamma.max = 5
        ui.set_par('ggal.nh', val=2.0, frozen=True)

        return fname, (grp, qual), \
            self._add_datadir_path(_canonical_pha_grouped)

    def _setup_usermodel(self):
        """Try a user model.
        """

        ui.clean()
        # Note: array is not sorted on purpose, and float/int
        # values.
        ui.load_arrays(3, [1, 12.2, 2, 14], [4, 8, 12, 4])

        def mymodel_func(pars, x, xhi=None):
            return pars[0] + pars[1] * x

        ui.load_user_model(mymodel_func, "mymodel")
        ui.add_user_pars("mymodel",
                         parnames=["c", "m"],
                         parvals=[2, 0.5],
                         parmins=[-10, 0],
                         parmaxs=[10, 5.5],
                         parunits=["m", ""],
                         parfrozen=[False, True])

        ui.set_source(3, ui.sin.sin_model + mymodel)

        ui.set_stat('cash')
        ui.set_method('simplex')

    def test_restore_empty(self):
        "Can the empty state be evaluated?"

        ui.clean()

        # At present the only check is that the file can be
        # loaded.
        self._restore()

    def test_canonical_empty(self):
        self._compare(_canonical_empty)

    def test_canonical_empty_outfile(self):
        tfile = tempfile.NamedTemporaryFile(suffix='.sherpa')
        ui.save_all(tfile.name, clobber=True)
        with open(tfile.name, 'r') as fh:
            output = fh.read()
        self._compare_lines(_canonical_empty, output)

    def test_canonical_empty_stats(self):

        ui.set_stat('leastsq')

        ui.set_method('simplex')
        ui.set_method_opt('maxfev', 5000)
        ui.set_method_opt('verbose', 1)

        self._compare(_canonical_empty_stats)

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_canonical_pha_basic(self):

        _, canonical = self._setup_pha_basic()
        self._compare(canonical)

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_restore_pha_basic(self):
        "Can the state be evaluated?"

        fname, _ = self._setup_pha_basic()
        statval = ui.calc_stat()

        self._restore()

        self.assertEqual([1], ui.list_data_ids())
        self.assertEqual(fname, ui.get_data(1).name)
        self.assertTrue(ui.get_data().subtracted,
                        msg='Data should be subtracted')

        src_expr = ui.get_source()
        self.assertEqual(src_expr.name,
                         '(xsphabs.gal * (powlaw1d.pl + xsapec.src))')
        self.assertEqual(gal.name, 'xsphabs.gal')
        self.assertEqual(pl.name, 'powlaw1d.pl')
        self.assertEqual(src.name, 'xsapec.src')

        self.assertEqual(ui.calc_stat(), statval)

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_canonical_pha_grouped(self):

        _, _, canonical = self._setup_pha_grouped()
        self._compare(canonical)

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_restore_pha_grouped(self):
        "Can the state be evaluated?"

        fname, (grp,qual), _ = self._setup_pha_grouped()
        statval = ui.calc_stat('grp')

        self._restore()

        self.assertEqual(['grp'], ui.list_data_ids())
        self.assertEqual(fname, ui.get_data('grp').name)
        self.assertTrue(ui.get_data('grp').subtracted,
                        msg='Data should be subtracted')

        # TODO: add in a test of grouping and quality arrays

        src_expr = ui.get_source('grp')
        self.assertEqual(src_expr.name,
                         '(xsphabs.ggal * powlaw1d.gpl)')
        self.assertTrue(ggal.nh.frozen, msg="is ggal.nh frozen?")
        self.assertEqual(ggal.nh.val, 2.0)
        self.assertEqual(gpl.gamma.max, 5.0)

        self.assertEqual(ui.calc_stat('grp'), statval)

    # Since the code can NOT be restored, we currently do not
    # try to load the script
    def test_canonical_usermodel(self):

        self._setup_usermodel()
        self._compare(_canonical_usermodel)

if __name__ == '__main__':

    import sys
    if len(sys.argv) > 1:
        datadir = sys.argv[1]
    else:
        datadir = None

    SherpaTest(ui).test(datadir=datadir)
