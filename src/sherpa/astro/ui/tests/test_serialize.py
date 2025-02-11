#
#  Copyright (C) 2015, 2016, 2018 - 2021, 2023, 2024
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

"""Test the functions in sherpa.astro.ui.serialize, also testing the
corresponding functions in sherpa.astro.ui.utils.
"""

from io import StringIO
import re
import warnings

import numpy
from numpy.testing import assert_array_equal

import pytest

from sherpa.astro.models import JDPileup
from sherpa.astro import ui
from sherpa.astro.io.wcs import WCS

from sherpa.models.basic import TableModel

from sherpa.utils.err import ArgumentErr, DataErr, \
    IdentifierErr, IOErr, StatErr
from sherpa.utils.testing import get_datadir, requires_data, \
    requires_xspec, has_package_from_list, requires_fits, \
    requires_group, requires_region, requires_wcs


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

# Note that the XSPEC defaults changed in XSPEC 12.10.1 (xsxsect is
# now vern rather than bcmc, but it can depend on what files a
# user has - e.g. ~/.xspec/Xspec.init will over-ride it), so there
# are explicit statements to set the values before the tests so
# that we have a "known" state.
#

# A representation of the default Sherpa state
_canonical_empty = """import numpy
from sherpa.astro.ui import *


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

# Change a few settings for statistic/method
_canonical_empty_stats = """import numpy
from sherpa.astro.ui import *


######### Set Statistic

set_stat("leastsq")


######### Set Fitting Method

set_method("neldermead")

set_method_opt("finalsimplex", 9)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("initsimplex", 0)
set_method_opt("iquad", 1)
set_method_opt("maxfev", 5000)
set_method_opt("reflect", True)
set_method_opt("step", None)
set_method_opt("verbose", 1)

"""

_canonical_empty_iterstat = """import numpy
from sherpa.astro.ui import *


######### Set Statistic

set_stat("leastsq")


######### Set Fitting Method

set_method("neldermead")

set_method_opt("finalsimplex", 9)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("initsimplex", 0)
set_method_opt("iquad", 1)
set_method_opt("maxfev", 5000)
set_method_opt("reflect", True)
set_method_opt("step", None)
set_method_opt("verbose", 1)


######### Set Iterative Fitting Method

set_iter_method("sigmarej")

set_iter_method_opt("grow", 1)
set_iter_method_opt("hrej", 3)
set_iter_method_opt("lrej", 3)
set_iter_method_opt("maxiters", 5)

"""

# use @@ as a marker to indicate that the location of the data directory
# needs to be inserted into the text before testing.
#
_canonical_pha_basic = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha(1, "@@/3c273.pi")
group(1)

######### Data Spectral Responses

load_arf(1, "@@/3c273.arf", resp_id=1)
load_rmf(1, "@@/3c273.rmf", resp_id=1)

######### Load Background Data Sets

load_bkg(1, "@@/3c273_bg.pi", bkg_id=1)
group(1, bkg_id=1)

######### Background Spectral Responses

load_arf(1, "@@/3c273.arf", resp_id=1, bkg_id=1)
load_rmf(1, "@@/3c273.rmf", resp_id=1, bkg_id=1)

######### Set Energy or Wave Units

set_analysis(1, quantity="energy", type="rate", factor=0)
subtract(1)

######### Filter Data

notice_id(1, "0.467200011015:9.869600296021")


######### Set Statistic

set_stat("chi2datavar")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)


######### Set Model Components and Parameters

create_model_component("xsapec", "src")
src.integrate = True

src.kT.default_val = 1.0
src.kT.default_min = 0.0080000000000000002
src.kT.default_max = 64.0
src.kT.val     = 1.0
src.kT.min     = 0.0080000000000000002
src.kT.max     = 64.0
src.kT.units   = "keV"
src.kT.frozen  = False

src.Abundanc.default_val = 1.0
src.Abundanc.default_min = 0.0
src.Abundanc.default_max = 5.0
src.Abundanc.val     = 1.0
src.Abundanc.min     = 0.0
src.Abundanc.max     = 5.0
src.Abundanc.units   = ""
src.Abundanc.frozen  = True

src.Redshift.default_val = 0.0
src.Redshift.default_min = -0.999
src.Redshift.default_max = 10.0
src.Redshift.val     = 0.0
src.Redshift.min     = -0.999
src.Redshift.max     = 10.0
src.Redshift.units   = ""
src.Redshift.frozen  = True

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

pl.gamma.default_val = 1.0
pl.gamma.default_min = -10.0
pl.gamma.default_max = 10.0
pl.gamma.val     = 1.0
pl.gamma.min     = -10.0
pl.gamma.max     = 10.0
pl.gamma.units   = ""
pl.gamma.frozen  = False

pl.ref.default_val = 1.0
pl.ref.default_min = -3.4028234663852886e+38
pl.ref.default_max = 3.4028234663852886e+38
pl.ref.val     = 1.0
pl.ref.min     = -3.4028234663852886e+38
pl.ref.max     = 3.4028234663852886e+38
pl.ref.units   = ""
pl.ref.frozen  = True

pl.ampl.default_val = 1.0
pl.ampl.default_min = 0.0
pl.ampl.default_max = 3.4028234663852886e+38
pl.ampl.val     = 1.0
pl.ampl.min     = 0.0
pl.ampl.max     = 3.4028234663852886e+38
pl.ampl.units   = ""
pl.ampl.frozen  = False

create_model_component("xsphabs", "gal")
gal.integrate = True

gal.nH.default_val = 1.0
gal.nH.default_min = 0.0
gal.nH.default_max = 1000000.0
gal.nH.val     = 1.0
gal.nH.min     = 0.0
gal.nH.max     = 1000000.0
gal.nH.units   = "10^22 atoms / cm^2"
gal.nH.frozen  = False



######### Set Source, Pileup and Background Models

set_source(1, xsphabs.gal * (powlaw1d.pl + xsapec.src))

"""

_canonical_pha_grouped = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha("grp", "@@/3c273.pi")

######### Data grouping flags

set_grouping("grp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16))

######### Data quality flags

set_quality("grp", val=numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], numpy.int16))
group("grp")

######### Data Spectral Responses

load_arf("grp", "@@/3c273.arf", resp_id=1)
load_rmf("grp", "@@/3c273.rmf", resp_id=1)

######### Load Background Data Sets

load_bkg("grp", "@@/3c273_bg.pi", bkg_id=1)

######### Background grouping flags

set_grouping("grp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)

######### Background quality flags

set_quality("grp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)
group("grp", bkg_id=1)

######### Background Spectral Responses

load_arf("grp", "@@/3c273.arf", resp_id=1, bkg_id=1)
load_rmf("grp", "@@/3c273.rmf", resp_id=1, bkg_id=1)

######### Set Energy or Wave Units

set_analysis("grp", quantity="energy", type="rate", factor=0)
subtract("grp")

######### Filter Data

notice_id("grp", "0.467200011015:6.365600109100")


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)


######### Set Model Components and Parameters

create_model_component("powlaw1d", "gpl")
gpl.integrate = True

gpl.gamma.default_val = 1.0
gpl.gamma.default_min = -10.0
gpl.gamma.default_max = 10.0
gpl.gamma.val     = 1.0
gpl.gamma.min     = -10.0
gpl.gamma.max     = 5.0
gpl.gamma.units   = ""
gpl.gamma.frozen  = False

gpl.ref.default_val = 1.0
gpl.ref.default_min = -3.4028234663852886e+38
gpl.ref.default_max = 3.4028234663852886e+38
gpl.ref.val     = 1.0
gpl.ref.min     = -3.4028234663852886e+38
gpl.ref.max     = 3.4028234663852886e+38
gpl.ref.units   = ""
gpl.ref.frozen  = True

gpl.ampl.default_val = 1.0
gpl.ampl.default_min = 0.0
gpl.ampl.default_max = 3.4028234663852886e+38
gpl.ampl.val     = 1.0
gpl.ampl.min     = 0.0
gpl.ampl.max     = 3.4028234663852886e+38
gpl.ampl.units   = ""
gpl.ampl.frozen  = False

create_model_component("xsphabs", "ggal")
ggal.integrate = True

ggal.nH.default_val = 2.0
ggal.nH.default_min = 0.0
ggal.nH.default_max = 1000000.0
ggal.nH.val     = 2.0
ggal.nH.min     = 0.0
ggal.nH.max     = 1000000.0
ggal.nH.units   = "10^22 atoms / cm^2"
ggal.nH.frozen  = True



######### Set Source, Pileup and Background Models

set_source("grp", xsphabs.ggal * powlaw1d.gpl)

"""

_canonical_pha_back = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha("bgrp", "@@/3c273.pi")
group("bgrp")

######### Data Spectral Responses

load_arf("bgrp", "@@/3c273.arf", resp_id=1)
load_rmf("bgrp", "@@/3c273.rmf", resp_id=1)

######### Load Background Data Sets

load_bkg("bgrp", "@@/3c273_bg.pi", bkg_id=1)

######### Background grouping flags

set_grouping("bgrp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)

######### Background quality flags

set_quality("bgrp", val=numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], numpy.int16), bkg_id=1)
group("bgrp", bkg_id=1)

######### Background Spectral Responses

load_arf("bgrp", "@@/3c273.arf", resp_id=1, bkg_id=1)
load_rmf("bgrp", "@@/3c273.rmf", resp_id=1, bkg_id=1)

######### Set Energy or Wave Units

set_analysis("bgrp", quantity="energy", type="rate", factor=0)

######### Filter Data

notice_id("bgrp", "0.467200011015:6.570000171661")
notice_id("bgrp", bkg_id=1)
notice_id("bgrp", "1.605999946594:8.760000228882", bkg_id=1)


######### Set Statistic

set_stat("chi2xspecvar")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)


######### Set Model Components and Parameters

create_model_component("powlaw1d", "gpl")
gpl.integrate = True

gpl.gamma.default_val = 1.0
gpl.gamma.default_min = -10.0
gpl.gamma.default_max = 10.0
gpl.gamma.val     = 1.0
gpl.gamma.min     = -5.0
gpl.gamma.max     = 10.0
gpl.gamma.units   = ""
gpl.gamma.frozen  = False

gpl.ref.default_val = 1.0
gpl.ref.default_min = -3.4028234663852886e+38
gpl.ref.default_max = 3.4028234663852886e+38
gpl.ref.val     = 1.0
gpl.ref.min     = -3.4028234663852886e+38
gpl.ref.max     = 3.4028234663852886e+38
gpl.ref.units   = ""
gpl.ref.frozen  = True

gpl.ampl.default_val = 1.0
gpl.ampl.default_min = 0.0
gpl.ampl.default_max = 3.4028234663852886e+38
gpl.ampl.val     = 1.0
gpl.ampl.min     = 0.0
gpl.ampl.max     = 3.4028234663852886e+38
gpl.ampl.units   = ""
gpl.ampl.frozen  = False

create_model_component("xsphabs", "ggal")
ggal.integrate = True

ggal.nH.default_val = 2.0
ggal.nH.default_min = 0.0
ggal.nH.default_max = 1000000.0
ggal.nH.val     = 2.0
ggal.nH.min     = 0.0
ggal.nH.max     = 1000000.0
ggal.nH.units   = "10^22 atoms / cm^2"
ggal.nH.frozen  = True

create_model_component("steplo1d", "bstep")
bstep.integrate = True

bstep.xcut.default_val = 0.0
bstep.xcut.default_min = -3.4028234663852886e+38
bstep.xcut.default_max = 3.4028234663852886e+38
bstep.xcut.val     = 0.0
bstep.xcut.min     = -3.4028234663852886e+38
bstep.xcut.max     = 3.4028234663852886e+38
bstep.xcut.units   = ""
bstep.xcut.frozen  = False

bstep.ampl.default_val = 1.0
bstep.ampl.default_min = 0.0
bstep.ampl.default_max = 3.4028234663852886e+38
bstep.ampl.val     = 1.0
bstep.ampl.min     = 0.0
bstep.ampl.max     = 3.4028234663852886e+38
bstep.ampl.units   = ""
bstep.ampl.frozen  = False

create_model_component("polynom1d", "bpoly")
bpoly.integrate = True

bpoly.c0.default_val = 1.0
bpoly.c0.default_min = -3.4028234663852886e+38
bpoly.c0.default_max = 3.4028234663852886e+38
bpoly.c0.val     = 1.0
bpoly.c0.min     = -3.4028234663852886e+38
bpoly.c0.max     = 3.4028234663852886e+38
bpoly.c0.units   = ""
bpoly.c0.frozen  = True

bpoly.c1.default_val = 0.0
bpoly.c1.default_min = -3.4028234663852886e+38
bpoly.c1.default_max = 3.4028234663852886e+38
bpoly.c1.val     = 0.0
bpoly.c1.min     = -3.4028234663852886e+38
bpoly.c1.max     = 3.4028234663852886e+38
bpoly.c1.units   = ""
bpoly.c1.frozen  = True

bpoly.c2.default_val = 0.0
bpoly.c2.default_min = -3.4028234663852886e+38
bpoly.c2.default_max = 3.4028234663852886e+38
bpoly.c2.val     = 0.0
bpoly.c2.min     = -3.4028234663852886e+38
bpoly.c2.max     = 3.4028234663852886e+38
bpoly.c2.units   = ""
bpoly.c2.frozen  = True

bpoly.c3.default_val = 0.0
bpoly.c3.default_min = -3.4028234663852886e+38
bpoly.c3.default_max = 3.4028234663852886e+38
bpoly.c3.val     = 0.0
bpoly.c3.min     = -3.4028234663852886e+38
bpoly.c3.max     = 3.4028234663852886e+38
bpoly.c3.units   = ""
bpoly.c3.frozen  = True

bpoly.c4.default_val = 0.0
bpoly.c4.default_min = -3.4028234663852886e+38
bpoly.c4.default_max = 3.4028234663852886e+38
bpoly.c4.val     = 0.0
bpoly.c4.min     = -3.4028234663852886e+38
bpoly.c4.max     = 3.4028234663852886e+38
bpoly.c4.units   = ""
bpoly.c4.frozen  = True

bpoly.c5.default_val = 0.0
bpoly.c5.default_min = -3.4028234663852886e+38
bpoly.c5.default_max = 3.4028234663852886e+38
bpoly.c5.val     = 0.0
bpoly.c5.min     = -3.4028234663852886e+38
bpoly.c5.max     = 3.4028234663852886e+38
bpoly.c5.units   = ""
bpoly.c5.frozen  = True

bpoly.c6.default_val = 0.0
bpoly.c6.default_min = -3.4028234663852886e+38
bpoly.c6.default_max = 3.4028234663852886e+38
bpoly.c6.val     = 0.0
bpoly.c6.min     = -3.4028234663852886e+38
bpoly.c6.max     = 3.4028234663852886e+38
bpoly.c6.units   = ""
bpoly.c6.frozen  = True

bpoly.c7.default_val = 0.0
bpoly.c7.default_min = -3.4028234663852886e+38
bpoly.c7.default_max = 3.4028234663852886e+38
bpoly.c7.val     = 0.0
bpoly.c7.min     = -3.4028234663852886e+38
bpoly.c7.max     = 3.4028234663852886e+38
bpoly.c7.units   = ""
bpoly.c7.frozen  = True

bpoly.c8.default_val = 0.0
bpoly.c8.default_min = -3.4028234663852886e+38
bpoly.c8.default_max = 3.4028234663852886e+38
bpoly.c8.val     = 0.0
bpoly.c8.min     = -3.4028234663852886e+38
bpoly.c8.max     = 3.4028234663852886e+38
bpoly.c8.units   = ""
bpoly.c8.frozen  = True

bpoly.offset.default_val = 0.0
bpoly.offset.default_min = -3.4028234663852886e+38
bpoly.offset.default_max = 3.4028234663852886e+38
bpoly.offset.val     = 0.0
bpoly.offset.min     = -3.4028234663852886e+38
bpoly.offset.max     = 3.4028234663852886e+38
bpoly.offset.units   = ""
bpoly.offset.frozen  = True



######### Set Source, Pileup and Background Models

set_source("bgrp", xsphabs.ggal * powlaw1d.gpl)

set_bkg_source("bgrp", steplo1d.bstep + polynom1d.bpoly, bkg_id=1)


######### XSPEC Module Settings

set_xschatter(0)
set_xsabund("lodd")
set_xscosmo(72, 0.02, 0.71)
set_xsxsect("vern")
"""

_canonical_pha_no_response = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha("x", "@@/source1.pi")

######### Set Energy or Wave Units

set_analysis("x", quantity="channel", type="counts", factor=0)

######### Filter Data

notice_id("x", "40:299,320:400")


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

_canonical_pha_full_model = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha(1, "@@/3c273.pi")
group(1)

######### Data Spectral Responses

load_arf(1, "@@/3c273.arf", resp_id=1)
load_rmf(1, "@@/3c273.rmf", resp_id=1)

######### Load Background Data Sets

load_bkg(1, "@@/3c273_bg.pi", bkg_id=1)
group(1, bkg_id=1)

######### Background Spectral Responses

load_arf(1, "@@/3c273.arf", resp_id=1, bkg_id=1)
load_rmf(1, "@@/3c273.rmf", resp_id=1, bkg_id=1)

######### Set Energy or Wave Units

set_analysis(1, quantity="energy", type="rate", factor=0)

######### Filter Data

notice_id(1, "0.992799997330:6.570000171661")


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)


######### Set Model Components and Parameters

create_model_component("powlaw1d", "pl")
pl.integrate = True

pl.gamma.default_val = 1.7
pl.gamma.default_min = -10.0
pl.gamma.default_max = 10.0
pl.gamma.val     = 1.7
pl.gamma.min     = -10.0
pl.gamma.max     = 10.0
pl.gamma.units   = ""
pl.gamma.frozen  = False

pl.ref.default_val = 1.0
pl.ref.default_min = -3.4028234663852886e+38
pl.ref.default_max = 3.4028234663852886e+38
pl.ref.val     = 1.0
pl.ref.min     = -3.4028234663852886e+38
pl.ref.max     = 3.4028234663852886e+38
pl.ref.units   = ""
pl.ref.frozen  = True

pl.ampl.default_val = 0.0001
pl.ampl.default_min = 0.0
pl.ampl.default_max = 3.4028234663852886e+38
pl.ampl.val     = 0.0001
pl.ampl.min     = 0.0
pl.ampl.max     = 3.4028234663852886e+38
pl.ampl.units   = ""
pl.ampl.frozen  = False

create_model_component("polynom1d", "con")
con.integrate = True

con.c0.default_val = 0.001
con.c0.default_min = -3.4028234663852886e+38
con.c0.default_max = 3.4028234663852886e+38
con.c0.val     = 0.001
con.c0.min     = -3.4028234663852886e+38
con.c0.max     = 3.4028234663852886e+38
con.c0.units   = ""
con.c0.frozen  = False

con.c1.default_val = 0.002
con.c1.default_min = -3.4028234663852886e+38
con.c1.default_max = 3.4028234663852886e+38
con.c1.val     = 0.002
con.c1.min     = -3.4028234663852886e+38
con.c1.max     = 3.4028234663852886e+38
con.c1.units   = ""
con.c1.frozen  = True

con.c2.default_val = 0.0
con.c2.default_min = -3.4028234663852886e+38
con.c2.default_max = 3.4028234663852886e+38
con.c2.val     = 0.0
con.c2.min     = -3.4028234663852886e+38
con.c2.max     = 3.4028234663852886e+38
con.c2.units   = ""
con.c2.frozen  = True

con.c3.default_val = 0.0
con.c3.default_min = -3.4028234663852886e+38
con.c3.default_max = 3.4028234663852886e+38
con.c3.val     = 0.0
con.c3.min     = -3.4028234663852886e+38
con.c3.max     = 3.4028234663852886e+38
con.c3.units   = ""
con.c3.frozen  = True

con.c4.default_val = 0.0
con.c4.default_min = -3.4028234663852886e+38
con.c4.default_max = 3.4028234663852886e+38
con.c4.val     = 0.0
con.c4.min     = -3.4028234663852886e+38
con.c4.max     = 3.4028234663852886e+38
con.c4.units   = ""
con.c4.frozen  = True

con.c5.default_val = 0.0
con.c5.default_min = -3.4028234663852886e+38
con.c5.default_max = 3.4028234663852886e+38
con.c5.val     = 0.0
con.c5.min     = -3.4028234663852886e+38
con.c5.max     = 3.4028234663852886e+38
con.c5.units   = ""
con.c5.frozen  = True

con.c6.default_val = 0.0
con.c6.default_min = -3.4028234663852886e+38
con.c6.default_max = 3.4028234663852886e+38
con.c6.val     = 0.0
con.c6.min     = -3.4028234663852886e+38
con.c6.max     = 3.4028234663852886e+38
con.c6.units   = ""
con.c6.frozen  = True

con.c7.default_val = 0.0
con.c7.default_min = -3.4028234663852886e+38
con.c7.default_max = 3.4028234663852886e+38
con.c7.val     = 0.0
con.c7.min     = -3.4028234663852886e+38
con.c7.max     = 3.4028234663852886e+38
con.c7.units   = ""
con.c7.frozen  = True

con.c8.default_val = 0.0
con.c8.default_min = -3.4028234663852886e+38
con.c8.default_max = 3.4028234663852886e+38
con.c8.val     = 0.0
con.c8.min     = -3.4028234663852886e+38
con.c8.max     = 3.4028234663852886e+38
con.c8.units   = ""
con.c8.frozen  = True

con.offset.default_val = 0.0
con.offset.default_min = -3.4028234663852886e+38
con.offset.default_max = 3.4028234663852886e+38
con.offset.val     = 0.0
con.offset.min     = -3.4028234663852886e+38
con.offset.max     = 3.4028234663852886e+38
con.offset.units   = ""
con.offset.frozen  = True



######### Set Source, Pileup and Background Models

set_full_model(1, (apply_rmf(apply_arf((38564.6089269 * powlaw1d.pl))) + polynom1d.con))

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


######### Set Statistic

set_stat("cash")


######### Set Fitting Method

set_method("neldermead")

set_method_opt("finalsimplex", 9)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("initsimplex", 0)
set_method_opt("iquad", 1)
set_method_opt("maxfev", None)
set_method_opt("reflect", True)
set_method_opt("step", None)
set_method_opt("verbose", 0)


######### Set Model Components and Parameters

create_model_component("sin", "sin_model")
sin_model.integrate = True

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

sin_model.ampl.default_val = 1.0
sin_model.ampl.default_min = 1.0000000000000001e-05
sin_model.ampl.default_max = 3.4028234663852886e+38
sin_model.ampl.val     = 1.0
sin_model.ampl.min     = 1.0000000000000001e-05
sin_model.ampl.max     = 3.4028234663852886e+38
sin_model.ampl.units   = ""
sin_model.ampl.frozen  = False

print("Found user model 'mymodel'; please check it is saved correctly.")
def mymodel_func(pars, x, xhi=None):
    return pars[0] + pars[1] * x

load_user_model(mymodel_func, "mymodel")
add_user_pars("mymodel",
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

set_source(3, sin.sin_model + usermodel.mymodel)

"""

# An image file with no filter and no model
#
_canonical_img_no_filter_no_model = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_image(1, "@@/img.fits")

######### Set Image Coordinates

set_coord(1, 'logical')


######### Set Statistic

set_stat("cstat")


######### Set Fitting Method

set_method("neldermead")

set_method_opt("finalsimplex", 9)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("initsimplex", 0)
set_method_opt("iquad", 1)
set_method_opt("maxfev", None)
set_method_opt("reflect", True)
set_method_opt("step", None)
set_method_opt("verbose", 0)

"""

# An image file with filter and model
#
_canonical_img_filter_model = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_image(1, "@@/img.fits")

######### Set Image Coordinates

set_coord(1, 'logical')

######### Filter Data

notice2d_id(1, "Circle(50,50,30)&!RotBox(30,30,10,5,45)|Ellipse(40,75,30,20,320)")


######### Set Statistic

set_stat("cstat")


######### Set Fitting Method

set_method("neldermead")

set_method_opt("finalsimplex", 9)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("initsimplex", 0)
set_method_opt("iquad", 1)
set_method_opt("maxfev", None)
set_method_opt("reflect", True)
set_method_opt("step", None)
set_method_opt("verbose", 0)


######### Set Model Components and Parameters

create_model_component("gauss2d", "gmdl")
gmdl.integrate = True

gmdl.fwhm.default_val = 10.0
gmdl.fwhm.default_min = 1.1754943508222875e-38
gmdl.fwhm.default_max = 3.4028234663852886e+38
gmdl.fwhm.val     = 10.0
gmdl.fwhm.min     = 1.1754943508222875e-38
gmdl.fwhm.max     = 3.4028234663852886e+38
gmdl.fwhm.units   = ""
gmdl.fwhm.frozen  = False

gmdl.xpos.default_val = 51.0
gmdl.xpos.default_min = -3.4028234663852886e+38
gmdl.xpos.default_max = 3.4028234663852886e+38
gmdl.xpos.val     = 51.0
gmdl.xpos.min     = -3.4028234663852886e+38
gmdl.xpos.max     = 3.4028234663852886e+38
gmdl.xpos.units   = ""
gmdl.xpos.frozen  = False

gmdl.ypos.default_val = 49.0
gmdl.ypos.default_min = -3.4028234663852886e+38
gmdl.ypos.default_max = 3.4028234663852886e+38
gmdl.ypos.val     = 49.0
gmdl.ypos.min     = -3.4028234663852886e+38
gmdl.ypos.max     = 3.4028234663852886e+38
gmdl.ypos.units   = ""
gmdl.ypos.frozen  = False

gmdl.ellip.default_val = 0.80000000000000004
gmdl.ellip.default_min = 0.0
gmdl.ellip.default_max = 0.999
gmdl.ellip.val     = 0.80000000000000004
gmdl.ellip.min     = 0.0
gmdl.ellip.max     = 0.999
gmdl.ellip.units   = ""
gmdl.ellip.frozen  = True

gmdl.theta.default_val = 1.2
gmdl.theta.default_min = -6.2831853071795862
gmdl.theta.default_max = 6.2831853071795862
gmdl.theta.val     = 1.2
gmdl.theta.min     = -6.2831853071795862
gmdl.theta.max     = 6.2831853071795862
gmdl.theta.units   = "radians"
gmdl.theta.frozen  = True

gmdl.ampl.default_val = 10.0
gmdl.ampl.default_min = -3.4028234663852886e+38
gmdl.ampl.default_max = 3.4028234663852886e+38
gmdl.ampl.val     = 10.0
gmdl.ampl.min     = -3.4028234663852886e+38
gmdl.ampl.max     = 3.4028234663852886e+38
gmdl.ampl.units   = ""
gmdl.ampl.frozen  = False

create_model_component("scale2d", "bmdl")
bmdl.integrate = False

bmdl.c0.default_val = 2.0
bmdl.c0.default_min = -3.4028234663852886e+38
bmdl.c0.default_max = 3.4028234663852886e+38
bmdl.c0.val     = 2.0
bmdl.c0.min     = -3.4028234663852886e+38
bmdl.c0.max     = 3.4028234663852886e+38
bmdl.c0.units   = ""
bmdl.c0.frozen  = False



######### Set Source, Pileup and Background Models

set_source(1, gauss2d.gmdl + scale2d.bmdl)

"""

# An image file with no filter, model and PSF
# See issue #1873
#
# This is currently wrong, since it does not include the "set the psf"
# code, but it's unclear what this should be, so just describe the
# existing output (the test that uses this fails because the
# reconstruction is not correct).
#
_canonical_img_no_filter_model_psf = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_image(1, "@@/img.fits")

######### Set Image Coordinates

set_coord(1, 'logical')


######### Set Statistic

set_stat("cstat")


######### Set Fitting Method

set_method("neldermead")

set_method_opt("finalsimplex", 9)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("initsimplex", 0)
set_method_opt("iquad", 1)
set_method_opt("maxfev", None)
set_method_opt("reflect", True)
set_method_opt("step", None)
set_method_opt("verbose", 0)


######### Set Model Components and Parameters

load_psf("p0", "@@/psf_0.0_00_bin1.img")
p0.size = (26, 26)
p0.center = (13, 13)

create_model_component("gauss2d", "g1")
g1.integrate = True

g1.fwhm.default_val = 5.0
g1.fwhm.default_min = 1.1754943508222875e-38
g1.fwhm.default_max = 3.4028234663852886e+38
g1.fwhm.val     = 5.0
g1.fwhm.min     = 1.1754943508222875e-38
g1.fwhm.max     = 3.4028234663852886e+38
g1.fwhm.units   = ""
g1.fwhm.frozen  = False

g1.xpos.default_val = 49.0
g1.xpos.default_min = -3.4028234663852886e+38
g1.xpos.default_max = 3.4028234663852886e+38
g1.xpos.val     = 49.0
g1.xpos.min     = -3.4028234663852886e+38
g1.xpos.max     = 3.4028234663852886e+38
g1.xpos.units   = ""
g1.xpos.frozen  = False

g1.ypos.default_val = 52.0
g1.ypos.default_min = -3.4028234663852886e+38
g1.ypos.default_max = 3.4028234663852886e+38
g1.ypos.val     = 52.0
g1.ypos.min     = -3.4028234663852886e+38
g1.ypos.max     = 3.4028234663852886e+38
g1.ypos.units   = ""
g1.ypos.frozen  = False

g1.ellip.default_val = 0.40000000000000002
g1.ellip.default_min = 0.0
g1.ellip.default_max = 0.999
g1.ellip.val     = 0.40000000000000002
g1.ellip.min     = 0.0
g1.ellip.max     = 0.999
g1.ellip.units   = ""
g1.ellip.frozen  = True

g1.theta.default_val = 1.2
g1.theta.default_min = -6.2831853071795862
g1.theta.default_max = 6.2831853071795862
g1.theta.val     = 1.2
g1.theta.min     = -6.2831853071795862
g1.theta.max     = 6.2831853071795862
g1.theta.units   = "radians"
g1.theta.frozen  = True

g1.ampl.default_val = 100.0
g1.ampl.default_min = -3.4028234663852886e+38
g1.ampl.default_max = 3.4028234663852886e+38
g1.ampl.val     = 100.0
g1.ampl.min     = -3.4028234663852886e+38
g1.ampl.max     = 3.4028234663852886e+38
g1.ampl.units   = ""
g1.ampl.frozen  = False



######### Associate PSF models with the datasets

set_psf(1, p0)

######### Set Source, Pileup and Background Models

set_source(1, gauss2d.g1)

"""

_canonical_table_model = """import numpy
from sherpa.astro.ui import *


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)


######### Set Model Components and Parameters

load_table_model("tbl", "@@/test_rmfimg.fits")
tbl.integrate = True

tbl.ampl.default_val = 10.0
tbl.ampl.default_min = -3.4028234663852886e+38
tbl.ampl.default_max = 3.4028234663852886e+38
tbl.ampl.val     = 10.0
tbl.ampl.min     = 0.0
tbl.ampl.max     = 20.0
tbl.ampl.units   = ""
tbl.ampl.frozen  = True


"""

_canonical_xstable_model = """import numpy
from sherpa.astro.ui import *


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)


######### Set Model Components and Parameters

load_xstable_model("tbl", "@@/testpcfabs.mod")
tbl.integrate = True

tbl.nh.default_val = 2.0
tbl.nh.default_min = 0.0
tbl.nh.default_max = 1000.0
tbl.nh.val     = 2.0
tbl.nh.min     = 0.0
tbl.nh.max     = 1000.0
tbl.nh.units   = ""
tbl.nh.frozen  = False

tbl.fract.default_val = 0.20000000000000001
tbl.fract.default_min = 0.0
tbl.fract.default_max = 1.0
tbl.fract.val     = 0.20000000000000001
tbl.fract.min     = 0.0
tbl.fract.max     = 1.0
tbl.fract.units   = ""
tbl.fract.frozen  = True


"""

_canonical_pileup_model = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_arrays(1,
            [1, 2, 3],
            [2, 5, 2],
            DataPHA)

######### Set Energy or Wave Units

set_analysis(1, quantity="channel", type="rate", factor=0)
load_arrays(2,
            [1, 2, 3],
            [3, 4, 1],
            DataPHA)

######### Set Energy or Wave Units

set_analysis(2, quantity="channel", type="rate", factor=0)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)


######### Set Model Components and Parameters

create_model_component("jdpileup", "pmod")
pmod.integrate = True

pmod.alpha.default_val = 0.5
pmod.alpha.default_min = 0.0
pmod.alpha.default_max = 1.0
pmod.alpha.val     = 0.5
pmod.alpha.min     = 0.0
pmod.alpha.max     = 1.0
pmod.alpha.units   = ""
pmod.alpha.frozen  = False

pmod.g0.default_val = 1.0
pmod.g0.default_min = 1.1754943508222875e-38
pmod.g0.default_max = 1.0
pmod.g0.val     = 1.0
pmod.g0.min     = 1.1754943508222875e-38
pmod.g0.max     = 1.0
pmod.g0.units   = ""
pmod.g0.frozen  = True

pmod.f.default_val = 0.94999999999999996
pmod.f.default_min = 0.90000000000000002
pmod.f.default_max = 1.0
pmod.f.val     = 0.94999999999999996
pmod.f.min     = 0.90000000000000002
pmod.f.max     = 1.0
pmod.f.units   = ""
pmod.f.frozen  = False

pmod.n.default_val = 1.0
pmod.n.default_min = 1.1754943508222875e-38
pmod.n.default_max = 100.0
pmod.n.val     = 1.0
pmod.n.min     = 1.1754943508222875e-38
pmod.n.max     = 100.0
pmod.n.units   = ""
pmod.n.frozen  = True

pmod.ftime.default_val = 3.2410000000000001
pmod.ftime.default_min = 1.1754943508222875e-38
pmod.ftime.default_max = 5.0
pmod.ftime.val     = 3.2410000000000001
pmod.ftime.min     = 1.1754943508222875e-38
pmod.ftime.max     = 5.0
pmod.ftime.units   = "sec"
pmod.ftime.frozen  = True

pmod.fracexp.default_val = 0.98699999999999999
pmod.fracexp.default_min = 0.0
pmod.fracexp.default_max = 1.0
pmod.fracexp.val     = 0.98699999999999999
pmod.fracexp.min     = 0.0
pmod.fracexp.max     = 1.0
pmod.fracexp.units   = ""
pmod.fracexp.frozen  = True

pmod.nterms.default_val = 30.0
pmod.nterms.default_min = 1.0
pmod.nterms.default_max = 100.0
pmod.nterms.val     = 30.0
pmod.nterms.min     = 1.0
pmod.nterms.max     = 100.0
pmod.nterms.units   = ""
pmod.nterms.frozen  = True



######### Set Source, Pileup and Background Models

set_pileup_model(1, jdpileup.pmod)
set_pileup_model(2, jdpileup.pmod)
"""

_canonical_dataspace1d_int = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_arrays(1,
            [1.0, 3.0, 5.0, 7.0, 9.0],
            [3.0, 5.0, 7.0, 9.0, 11.0],
            [2.0, 5.0, 6.0, 0.0, 2.0],
            Data1DInt)

######### Filter Data

notice_id(1, "1.0000:7.0000,9.0000:11.0000")


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

_canonical_dataspace2d_img = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_arrays(1,
            [1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0],
            [1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0],
            [1.0, 2.0, 1.0, 0.0, 2.0, 3.0, 5.0, 2.0, 1.0, 0.0, 1.0, 1.0],
            (3, 4),
            DataIMG)

######### Set Image Coordinates

set_coord(1, 'logical')

######### Filter Data

notice2d_id(1, "Field()&!Box(1.4,2.5,1.2,1.4)")


######### Set Statistic

set_stat("cash")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

_canonical_load_arrays_simple = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_arrays("f",
            [-50, -20],
            [-20000.0, 300000.0],
            Data1D)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

_canonical_load_arrays_pha = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_arrays(1,
            [1, 2, 3, 4, 5],
            [12, 2, 1, 0, 1],
            DataPHA)
set_exposure(1, 100)
set_backscal(1, 0.002)
set_areascal(1, 0.001)

######### Data grouping flags

set_grouping(1, val=numpy.array([1, 1, -1, 0, 1], numpy.int16))

######### Data quality flags

set_quality(1, val=numpy.array([0, 0, 0, 0, 2], numpy.int16))
group(1)

######### Set Energy or Wave Units

set_analysis(1, quantity="channel", type="rate", factor=0)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

_canonical_load_arrays_data2d = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_arrays(1,
            [1, 1, 2],
            [1, 2, 2],
            [3, 4, 5],
            Data2D)
set_staterror(1, [0.1, 0.1, 0.2])


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

# Set an XSPEC parameter beyond the default hard limits: min
_canonical_xspec_hard_limit_min = """import numpy
from sherpa.astro.ui import *


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1)
set_method_opt("factor", 1)
set_method_opt("ftol", 1)
set_method_opt("gtol", 1)
set_method_opt("maxfev", 1)
set_method_opt("numcores", 1)
set_method_opt("verbose", 1)
set_method_opt("xtol", 1)


######### Set Model Components and Parameters

create_model_component("xspowerlaw", "mdl")
mdl.integrate = True

mdl.PhoIndex.hard_min    = -5.0
mdl.PhoIndex.default_val = -5.0
mdl.PhoIndex.default_min = -3.0
mdl.PhoIndex.default_max = 10.0
mdl.PhoIndex.val     = -5.0
mdl.PhoIndex.min     = -5.0
mdl.PhoIndex.max     = 10.0
mdl.PhoIndex.units   = ""
mdl.PhoIndex.frozen  = True

mdl.norm.default_val = 1.0
mdl.norm.default_min = 0.0
mdl.norm.default_max = 9.9999999999999998e+23
mdl.norm.val     = 1.0
mdl.norm.min     = 0.0
mdl.norm.max     = 100.0
mdl.norm.units   = ""
mdl.norm.frozen  = False


"""

# Set an XSPEC parameter beyond the default hard limits: max
_canonical_xspec_hard_limit_max = """import numpy
from sherpa.astro.ui import *


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1)
set_method_opt("factor", 1)
set_method_opt("ftol", 1)
set_method_opt("gtol", 1)
set_method_opt("maxfev", 1)
set_method_opt("numcores", 1)
set_method_opt("verbose", 1)
set_method_opt("xtol", 1)


######### Set Model Components and Parameters

create_model_component("xspowerlaw", "mdl")
mdl.integrate = True

mdl.PhoIndex.hard_max    = 15.0
mdl.PhoIndex.default_val = 15.0
mdl.PhoIndex.default_min = -3.0
mdl.PhoIndex.default_max = 10.0
mdl.PhoIndex.val     = 15.0
mdl.PhoIndex.min     = -3.0
mdl.PhoIndex.max     = 15.0
mdl.PhoIndex.units   = ""
mdl.PhoIndex.frozen  = True

mdl.norm.default_val = 1.0
mdl.norm.default_min = 0.0
mdl.norm.default_max = 9.9999999999999998e+23
mdl.norm.val     = 1.0
mdl.norm.min     = 0.0
mdl.norm.max     = 100.0
mdl.norm.units   = ""
mdl.norm.frozen  = False


"""

_canonical_link_par = """import numpy
from sherpa.astro.ui import *


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)


######### Set Model Components and Parameters

create_model_component("const1d", "sep")
sep.integrate = True

sep.c0.default_val = 100.0
sep.c0.default_min = -3.4028234663852886e+38
sep.c0.default_max = 3.4028234663852886e+38
sep.c0.val     = 100.0
sep.c0.min     = 100.0
sep.c0.max     = 100.0
sep.c0.units   = ""
sep.c0.frozen  = False

create_model_component("const1d", "m2")
m2.integrate = True

m2.c0.default_val = 300.0
m2.c0.default_min = -3.4028234663852886e+38
m2.c0.default_max = 3.4028234663852886e+38
m2.c0.val     = 300.0
m2.c0.min     = 10.0
m2.c0.max     = 500.0
m2.c0.units   = ""
m2.c0.frozen  = True

create_model_component("const1d", "m1")
m1.integrate = True

m1.c0.default_val = 200.0
m1.c0.default_min = -3.4028234663852886e+38
m1.c0.default_max = 3.4028234663852886e+38
m1.c0.val     = 200.0
m1.c0.min     = 0.0
m1.c0.max     = 1000.0
m1.c0.units   = ""
m1.c0.frozen  = False


link(m2.c0, m1.c0 + sep.c0)

"""

_canonical_load_data = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_data(1, "@@/data1.dat", ncols=3)


######### Set Statistic

set_stat("chi2datavar")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

_canonical_load_data_basic = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_data(1, "@@/double.dat")

######### Filter Data

notice_id(1, "20.1000:24.9000,27.1000:89.9000")


######### Set Statistic

set_stat("chi2")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

_canonical_pha_multiple_backgrounds = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha(1, "@@/3c120_meg_1.pha")

######### Data Spectral Responses

load_arf(1, "@@/3c120_meg_1.arf", resp_id=1)
load_rmf(1, "@@/3c120_meg_1.rmf", resp_id=1)

######### Load Background Data Sets

load_bkg(1, "@@/3c120_meg_1.pha", bkg_id=1)

######### Background Spectral Responses

load_arf(1, "@@/3c120_meg_1.arf", resp_id=1, bkg_id=1)
load_rmf(1, "@@/3c120_meg_1.rmf", resp_id=1, bkg_id=1)
load_bkg(1, "@@/3c120_meg_1.pha", bkg_id=2)

######### Background Spectral Responses

load_arf(1, "@@/3c120_meg_1.arf", resp_id=1, bkg_id=2)
load_rmf(1, "@@/3c120_meg_1.rmf", resp_id=1, bkg_id=2)

######### Set Energy or Wave Units

set_analysis(1, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(1, "2.000000000000:12.000000000000")
notice_id(1, bkg_id=1)
notice_id(1, "2.000000000000:12.000000000000", bkg_id=1)
notice_id(1, bkg_id=2)
notice_id(1, "2.000000000000:12.000000000000", bkg_id=2)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

_canonical_pha2 = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha(1, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(1, "@@/3c120_pha2", bkg_id=1)
load_bkg(1, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(1, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(1, "2.000000000000:4.000000000000")
notice_id(1, bkg_id=1)
notice_id(1, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(1, bkg_id=2)
notice_id(1, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(10, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(10, "@@/3c120_pha2", bkg_id=1)
load_bkg(10, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(10, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(10, "2.000000000000:4.000000000000")
notice_id(10, bkg_id=1)
notice_id(10, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(10, bkg_id=2)
notice_id(10, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(11, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(11, "@@/3c120_pha2", bkg_id=1)
load_bkg(11, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(11, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(11, "2.000000000000:4.000000000000")
notice_id(11, bkg_id=1)
notice_id(11, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(11, bkg_id=2)
notice_id(11, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(12, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(12, "@@/3c120_pha2", bkg_id=1)
load_bkg(12, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(12, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(12, "2.000000000000:4.000000000000")
notice_id(12, bkg_id=1)
notice_id(12, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(12, bkg_id=2)
notice_id(12, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(2, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(2, "@@/3c120_pha2", bkg_id=1)
load_bkg(2, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(2, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(2, "2.000000000000:4.000000000000")
notice_id(2, bkg_id=1)
notice_id(2, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(2, bkg_id=2)
notice_id(2, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(3, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(3, "@@/3c120_pha2", bkg_id=1)
load_bkg(3, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(3, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(3, "2.000000000000:4.000000000000")
notice_id(3, bkg_id=1)
notice_id(3, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(3, bkg_id=2)
notice_id(3, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(4, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(4, "@@/3c120_pha2", bkg_id=1)
load_bkg(4, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(4, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(4, "2.000000000000:4.000000000000")
notice_id(4, bkg_id=1)
notice_id(4, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(4, bkg_id=2)
notice_id(4, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(5, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(5, "@@/3c120_pha2", bkg_id=1)
load_bkg(5, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(5, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(5, "2.000000000000:4.000000000000")
notice_id(5, bkg_id=1)
notice_id(5, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(5, bkg_id=2)
notice_id(5, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(6, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(6, "@@/3c120_pha2", bkg_id=1)
load_bkg(6, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(6, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(6, "2.000000000000:4.000000000000")
notice_id(6, bkg_id=1)
notice_id(6, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(6, bkg_id=2)
notice_id(6, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(7, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(7, "@@/3c120_pha2", bkg_id=1)
load_bkg(7, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(7, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(7, "2.000000000000:4.000000000000")
notice_id(7, bkg_id=1)
notice_id(7, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(7, bkg_id=2)
notice_id(7, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(8, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(8, "@@/3c120_pha2", bkg_id=1)
load_bkg(8, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(8, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(8, "2.000000000000:4.000000000000")
notice_id(8, bkg_id=1)
notice_id(8, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(8, bkg_id=2)
notice_id(8, "2.000000000000:4.000000000000", bkg_id=2)
load_pha(9, "@@/3c120_pha2")

######### Load Background Data Sets

load_bkg(9, "@@/3c120_pha2", bkg_id=1)
load_bkg(9, "@@/3c120_pha2", bkg_id=2)

######### Set Energy or Wave Units

set_analysis(9, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(9, "2.000000000000:4.000000000000")
notice_id(9, bkg_id=1)
notice_id(9, "2.000000000000:4.000000000000", bkg_id=1)
notice_id(9, bkg_id=2)
notice_id(9, "2.000000000000:4.000000000000", bkg_id=2)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("gtol", 1.19209289551e-07)
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)

"""

_canonical_pha_csc = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha("csc", "@@/acisf01575_001N001_r0085_pha3.fits.gz")

######### Data grouping flags

set_grouping("csc", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, 1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16))

######### Data quality flags

set_quality("csc", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16))
group("csc")

######### Data Spectral Responses

load_arf("csc", "@@/acisf01575_001N001_r0085_arf3.fits", resp_id=1)
load_rmf("csc", "@@/acisf01575_001N001_r0085_rmf3.fits", resp_id=1)

######### Load Background Data Sets

load_bkg("csc", "@@/acisf01575_001N001_r0085_pha3.fits", bkg_id=1)

######### Background grouping flags

set_grouping("csc", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)

######### Background quality flags

set_quality("csc", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)
group("csc", bkg_id=1)

######### Background Spectral Responses

load_arf("csc", "@@/acisf01575_001N001_r0085_arf3.fits", resp_id=1, bkg_id=1)
load_rmf("csc", "@@/acisf01575_001N001_r0085_rmf3.fits", resp_id=1, bkg_id=1)

######### Set Energy or Wave Units

set_analysis("csc", quantity="energy", type="rate", factor=0)

######### Filter Data

notice_id("csc", "0.496399998665:7.007999897003")
notice_id("csc", bkg_id=1)
notice_id("csc", "0.394199997187:8.000800132751", bkg_id=1)


######### Set Statistic

set_stat("cstat")


######### Set Fitting Method

set_method("neldermead")

set_method_opt("finalsimplex", 9)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("initsimplex", 0)
set_method_opt("iquad", 1)
set_method_opt("maxfev", None)
set_method_opt("reflect", True)
set_method_opt("step", None)
set_method_opt("verbose", 0)


######### Set Model Components and Parameters

create_model_component("powlaw1d", "spl")
spl.integrate = True

spl.gamma.default_val = 1.6399999999999999
spl.gamma.default_min = -10.0
spl.gamma.default_max = 10.0
spl.gamma.val     = 1.6399999999999999
spl.gamma.min     = -10.0
spl.gamma.max     = 10.0
spl.gamma.units   = ""
spl.gamma.frozen  = False

spl.ref.default_val = 1.0
spl.ref.default_min = -3.4028234663852886e+38
spl.ref.default_max = 3.4028234663852886e+38
spl.ref.val     = 1.0
spl.ref.min     = -3.4028234663852886e+38
spl.ref.max     = 3.4028234663852886e+38
spl.ref.units   = ""
spl.ref.frozen  = True

spl.ampl.default_val = 3.8999999999999999e-05
spl.ampl.default_min = 0.0
spl.ampl.default_max = 3.4028234663852886e+38
spl.ampl.val     = 3.8999999999999999e-05
spl.ampl.min     = 0.0
spl.ampl.max     = 3.4028234663852886e+38
spl.ampl.units   = ""
spl.ampl.frozen  = False

create_model_component("xsphabs", "gal")
gal.integrate = True

gal.nH.default_val = 0.23999999999999999
gal.nH.default_min = 0.0
gal.nH.default_max = 1000000.0
gal.nH.val     = 0.23999999999999999
gal.nH.min     = 0.0
gal.nH.max     = 1000000.0
gal.nH.units   = "10^22 atoms / cm^2"
gal.nH.frozen  = False

create_model_component("powlaw1d", "bpl")
bpl.integrate = True

bpl.gamma.default_val = 0.87
bpl.gamma.default_min = -10.0
bpl.gamma.default_max = 10.0
bpl.gamma.val     = 0.87
bpl.gamma.min     = -10.0
bpl.gamma.max     = 10.0
bpl.gamma.units   = ""
bpl.gamma.frozen  = False

bpl.ref.default_val = 1.0
bpl.ref.default_min = -3.4028234663852886e+38
bpl.ref.default_max = 3.4028234663852886e+38
bpl.ref.val     = 1.0
bpl.ref.min     = -3.4028234663852886e+38
bpl.ref.max     = 3.4028234663852886e+38
bpl.ref.units   = ""
bpl.ref.frozen  = True

bpl.ampl.default_val = 2.1399999999999998e-06
bpl.ampl.default_min = 0.0
bpl.ampl.default_max = 3.4028234663852886e+38
bpl.ampl.val     = 2.1399999999999998e-06
bpl.ampl.min     = 0.0
bpl.ampl.max     = 3.4028234663852886e+38
bpl.ampl.units   = ""
bpl.ampl.frozen  = False



######### Set Source, Pileup and Background Models

set_source("csc", xsphabs.gal * powlaw1d.spl)

set_bkg_source("csc", powlaw1d.bpl, bkg_id=1)

"""

if has_xspec:
    from sherpa.astro import xspec

    _canonical_extra = """
######### XSPEC Module Settings

set_xschatter(0)
set_xsabund("angr")
set_xscosmo(70, 0, 0.73)
set_xsxsect("bcmc")
"""

    _canonical_pha_basic += _canonical_extra
    _canonical_pha_grouped += _canonical_extra
    _canonical_xstable_model += _canonical_extra
    _canonical_xspec_hard_limit_min += _canonical_extra
    _canonical_xspec_hard_limit_max += _canonical_extra
    _canonical_pha_csc += _canonical_extra

    del _canonical_extra


@pytest.fixture(autouse=True)
def setup(hide_logging, old_numpy_printing, clean_astro_ui):
    """Setup for all the tests"""

    if has_xspec:
        old_xspec = xspec.get_xsstate()

        # ensure we have the same settings as the test cases
        # used (changes to XSPEC may have changed the defaults)
        xspec.set_xschatter(0)
        xspec.set_xsabund('angr')
        xspec.set_xsxsect('bcmc')
    else:
        old_xspec = None

    # run the test
    yield

    if old_xspec is not None:
        xspec.set_xsstate(old_xspec)


def add_datadir_path(output):
    """Replace any @@ characters by the value of self.datadir,
    making sure that the replacement text does not end in a /."""

    dname = get_datadir()
    if dname.endswith('/'):
        dname = dname[:-1]

    return re.sub('@@', dname, output, count=0)


def compileit(output):
    # Let it just throw an exception in case of failure.
    compile(output, "test.py", "exec")


def compare_lines(expected, got):
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
        assert e == g

    # Do the line length after checking for the file
    # contents as it is easier to see what the difference
    # is this way around, since a difference in the
    # number of lines is often not very informative.
    #
    assert len(elines) == len(glines)


def compare(expected):
    """Run save_all and check the output (saved to a
    StringIO object) to the string value expected.
    """
    output = StringIO()
    ui.save_all(output)
    output = output.getvalue()

    # check the output is a valid Python program.
    # this check does not guard against potential issues,
    # but ensures that the program can compile.
    #
    compileit(output)
    compare_lines(expected, output)


def restore():
    """Run save_all then call clean and try to restore
    the Sherpa state from the saved file. Will raise
    a test failure if there was an error when
    executing the save file.
    """

    output = StringIO()
    ui.save_all(output)
    output = output.getvalue()
    ui.clean()

    exec(output)


def setup_pha_basic(make_data_path):
    """Load up a PHA file and make "simple" changes to the
    Sherpa state. Returns the name of the file that is
    loaded and the canonical output.
    """

    fname = make_data_path('3c273.pi')
    ui.load_pha(1, fname)
    ui.subtract()
    ui.set_stat('chi2datavar')
    ui.notice(0.5, 7)
    ui.set_source(ui.xsphabs.gal * (ui.powlaw1d.pl +
                                    ui.xsapec.src))
    return fname, add_datadir_path(_canonical_pha_basic)


def setup_pha_grouped(make_data_path):
    """Add in grouping and a few different choices.

    Returns the name of the file that is
    loaded, the new grouping and quality arrays,
    and the canonical output.
    """

    fname = make_data_path('3c273.pi')
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
    ui.powlaw1d.gpl.gamma.max = 5
    ui.set_par('ggal.nh', val=2.0, frozen=True)

    return fname, (grp, qual), \
        add_datadir_path(_canonical_pha_grouped)


def setup_pha_back(make_data_path):
    """Fit the background, rather than subtract it.
    """

    fname = make_data_path('3c273.pi')
    ui.load_pha('bgrp', fname)

    # Note: do not group the source dataset

    bchannels = ui.get_bkg('bgrp').channel

    bexclude = (bchannels < 10) | (bchannels > 850)
    bqual = bexclude * 1

    ui.group_counts('bgrp', 10, tabStops=bexclude, bkg_id=1)
    ui.set_quality('bgrp', bexclude, bkg_id=1)

    bgrp = ui.get_bkg('bgrp').grouping

    ui.set_stat('chi2xspecvar')

    # This call sets the noticed range for both source and
    # background data sets.
    ui.notice_id('bgrp', 0.5, 6)

    # Remove the "source" filter
    ui.notice_id('bgrp', None, None, bkg_id=1)
    ui.notice_id('bgrp', 2, 7, bkg_id=1)

    ui.set_source('bgrp', ui.xsphabs.ggal * ui.powlaw1d.gpl)
    ui.set_bkg_source('bgrp', ui.steplo1d.bstep + ui.polynom1d.bpoly)

    ui.set_xsabund('lodd')
    ui.set_xsxsect('vern')
    ui.set_xscosmo(72, 0.02, 0.71)

    ui.powlaw1d.gpl.gamma.min = -5
    ui.freeze(ui.polynom1d.bpoly.c0)

    ui.set_par('ggal.nh', val=2.0, frozen=True)

    return fname, (bgrp, bqual), \
        add_datadir_path(_canonical_pha_back)


def setup_usermodel():
    """Try a user model.
    """

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

    mymodel = ui.get_model_component("mymodel")
    ui.set_source(3, ui.sin.sin_model + mymodel)

    ui.set_stat('cash')
    ui.set_method('simplex')


def test_compile_failure():
    with pytest.raises(Exception):
        compileit("foo bar")


def test_restore_empty():
    "Can the empty state be evaluated?"

    # At present the only check is that the file can be
    # loaded.
    restore()


def test_save_all_checks_clobber(tmp_path):
    """We should error out if the output file exists"""
    outfile = tmp_path / "save.sherpa"
    outfile.write_text("x")

    with pytest.raises(IOErr,
                       match="^file '.*' exists and clobber is not set$"):
        ui.save_all(str(outfile))


def test_save_all_does_clobber(tmp_path):
    """We should clobber the file if told too"""
    outfile = tmp_path / "save.sherpa"
    outfile.write_text("x")
    ui.save_all(str(outfile), clobber=True)

    # All we care about is that the original text has been replaced.
    #
    cts = outfile.read_text()
    assert cts.startswith("import numpy\n")


def test_canonical_empty():
    "Contents of empty state are as expected"
    compare(_canonical_empty)


def test_canonical_empty_outfile(tmp_path):
    "Can read in a save file"

    outfile = tmp_path / "save.sherpa"
    ui.save_all(str(outfile))
    output = outfile.read_text()
    compare_lines(_canonical_empty, output)


def test_canonical_empty_stats():
    "Change several settings but load no data"

    ui.set_stat('leastsq')

    ui.set_method('simplex')
    ui.set_method_opt('maxfev', 5000)
    ui.set_method_opt('verbose', 1)

    compare(_canonical_empty_stats)


def test_canonical_empty_iterstat():
    "Check iterated-fit setting"

    ui.set_stat('leastsq')

    ui.set_method('simplex')
    ui.set_method_opt('maxfev', 5000)
    ui.set_method_opt('verbose', 1)

    ui.set_iter_method('sigmarej')
    ui.set_iter_method_opt('grow', 1)

    compare(_canonical_empty_iterstat)


@requires_data
@requires_xspec
@requires_fits
def test_canonical_pha_basic(make_data_path):

    _, canonical = setup_pha_basic(make_data_path)
    compare(canonical)


@requires_data
@requires_xspec
@requires_fits
def test_restore_pha_basic(make_data_path):
    "Can the state be evaluated?"

    fname, _ = setup_pha_basic(make_data_path)
    statval = ui.calc_stat()

    restore()

    assert ui.list_data_ids() == [1]
    assert ui.get_data(1).name == fname
    assert ui.get_data().subtracted, 'Data should be subtracted'

    src_expr = ui.get_source()
    assert src_expr.name == 'xsphabs.gal * (powlaw1d.pl + xsapec.src)'
    assert ui.xsphabs.gal.name == 'xsphabs.gal'
    assert ui.powlaw1d.pl.name == 'powlaw1d.pl'
    assert ui.xsapec.src.name == 'xsapec.src'

    assert ui.calc_stat() == pytest.approx(statval)


@requires_data
@requires_xspec
@requires_fits
@requires_group
def test_canonical_pha_grouped(make_data_path):

    _, _, canonical = setup_pha_grouped(make_data_path)
    compare(canonical)


@requires_data
@requires_xspec
@requires_fits
@requires_group
def test_restore_pha_grouped(make_data_path):
    "Can the state be evaluated?"

    fname, (grp, qual), _ = setup_pha_grouped(make_data_path)
    statval = ui.calc_stat('grp')

    restore()

    assert ui.list_data_ids() == ['grp']
    assert ui.get_data('grp').name == fname
    assert ui.get_data('grp').subtracted, 'Data should be subtracted'

    g = ui.get_grouping('grp')
    q = ui.get_quality('grp')
    assert g.dtype == numpy.int16
    assert q.dtype == numpy.int16

    assert_array_equal(grp, g, err_msg='grouping column')
    assert_array_equal(qual, q, err_msg='grouping column')

    src_expr = ui.get_source('grp')
    assert src_expr.name == 'xsphabs.ggal * powlaw1d.gpl'
    assert ui.xsphabs.ggal.nh.frozen, "is ggal.nh frozen?"
    assert ui.xsphabs.ggal.nh.val == 2.0
    assert ui.powlaw1d.gpl.gamma.max == 5.0

    assert ui.calc_stat('grp') == pytest.approx(statval)


@requires_data
@requires_xspec
@requires_fits
@requires_group
def test_canonical_pha_back(make_data_path):

    _, _, canonical = setup_pha_back(make_data_path)
    compare(canonical)


@requires_data
@requires_xspec
@requires_fits
@requires_group
def test_restore_pha_back(make_data_path):
    "Can the state be evaluated?"

    fname, (bgrp, bqual), _ = setup_pha_back(make_data_path)
    statval = ui.calc_stat('bgrp')

    assert ui.get_data("bgrp").get_filter(format="%.2f") == "0.47:6.57"
    assert ui.get_bkg("bgrp").get_filter(format="%.2f") == "1.61:8.76"

    restore()

    assert ui.list_data_ids() == ['bgrp']
    assert ui.get_data('bgrp').name == fname
    assert not ui.get_data('bgrp').subtracted, 'Data should not be subtracted'
    assert not ui.get_bkg('bgrp').subtracted, 'Background should not be subtracted'

    # TODO: at present the source is grouped; is this "correct"?
    assert ui.get_data('bgrp').grouped, 'Data should be grouped'  # FIXME?
    assert ui.get_bkg('bgrp').grouped, 'Background should be grouped'

    # g = ui.get_grouping('bgrp')
    # q = ui.get_quality('bgrp')
    # The data types are '>i2' / int16
    # assert g.dtype == numpy.int16
    # assert q.dtype == numpy.int16

    # TODO set up correct grouping bins...
    # nchan = ui.get_data('bgrp').channel.size
    # assert_array_equal(g, numpy.ones(nchan), err_msg='src grouping')
    # assert_array_equal(q, numpy.zeros(nchan), err_msg='src quality')

    bg = ui.get_grouping('bgrp', bkg_id=1)
    bq = ui.get_quality('bgrp', bkg_id=1)
    assert bg.dtype == numpy.int16
    assert bq.dtype == numpy.int16

    assert_array_equal(bg, bgrp, err_msg='bgnd grouping')
    assert_array_equal(bq, bqual, err_msg='bgnd quality')

    assert ui.get_data("bgrp").get_filter(format="%.2f") == "0.47:6.57"
    assert ui.get_bkg("bgrp").get_filter(format="%.2f") == "1.61:8.76"

    src_expr = ui.get_source('bgrp')
    assert src_expr.name == 'xsphabs.ggal * powlaw1d.gpl'

    bg_expr = ui.get_bkg_source('bgrp')
    assert bg_expr.name == 'steplo1d.bstep + polynom1d.bpoly'

    assert ui.xsphabs.ggal.nh.frozen, "is ggal.nh frozen?"
    assert ui.polynom1d.bpoly.c0.frozen, "is bpoly.c0 frozen?"
    assert ui.xsphabs.ggal.nh.val == 2.0
    assert ui.powlaw1d.gpl.gamma.min == -5.0

    assert ui.get_xsabund() == 'lodd'
    assert ui.get_xsxsect() == 'vern'
    cosmo = ui.get_xscosmo()
    assert cosmo[0] == pytest.approx(72.0)
    assert cosmo[1] == pytest.approx(0.02)
    assert cosmo[2] == pytest.approx(0.71)

    assert ui.calc_stat('bgrp') == pytest.approx(statval)


@requires_data
@requires_fits
def test_pha_no_response(make_data_path):
    """Check when there's no ARF / RMF or grouping, and switch to counts"""

    infile = make_data_path("source1.pi")
    ui.load_pha("x", infile)
    ui.notice(40, 400)
    ui.ignore(300, 319)
    ui.set_analysis("x", quantity="channel", type="counts")

    assert ui.get_dep("x", filter=True).size == 341
    assert ui.get_dep("x", filter=True).sum() == 2904
    assert ui.get_dep("x").sum() == 3328

    compare(add_datadir_path(_canonical_pha_no_response))

    restore()

    assert ui.get_dep("x", filter=True).size == 341
    assert ui.get_dep("x", filter=True).sum() == 2904
    assert ui.get_dep("x").sum() == 3328


def test_canonical_usermodel():
    "Can we save a usermodel?"
    setup_usermodel()
    compare(_canonical_usermodel)


def test_restore_usermodel():
    "Can the reload a usermodel"

    setup_usermodel()
    statval = ui.calc_stat(3)
    restore()

    # TODO: For the moment the source expression is created, in
    # the serialized form, using set_full_model. This should
    # be changed so that get_source can be used below.
    #
    # src_expr = ui.get_source(3)
    src_expr = ui.get_model(3)
    assert src_expr.name == 'sin.sin_model + usermodel.mymodel'
    mymodel = ui.get_model_component("mymodel")
    assert mymodel.m.frozen, "is mymodel.m frozen?"
    assert mymodel.c.val == 2.0
    assert mymodel.c.units == "m"
    assert mymodel.m.max == 5.5
    assert mymodel.m.units == ""

    assert ui.calc_stat(3) == pytest.approx(statval)


@requires_data
@requires_fits
def test_restore_img_no_filter_no_model(make_data_path):
    """Check issue #437"""

    ui.load_image(make_data_path('img.fits'))
    ui.set_stat('cstat')
    ui.set_method('simplex')

    sorig = ui.calc_data_sum2d()

    # sanity check
    assert sorig == pytest.approx(5041.44)
    assert ui.get_filter() == ''

    compare(add_datadir_path(_canonical_img_no_filter_no_model))

    restore()
    snew = ui.calc_data_sum2d()
    assert snew == pytest.approx(sorig)
    assert ui.get_filter() == ''


@requires_data
@requires_fits
@requires_region
def test_restore_img_filter_model(make_data_path):
    """Simple image check"""

    ui.load_image(make_data_path('img.fits'))

    ui.set_stat('cstat')
    ui.set_method('simplex')

    ui.notice2d('circle(50, 50, 30)')
    ui.notice2d('ellipse(40,75, 30, 20, 320)')
    ui.ignore2d('rotbox(30, 30, 10, 5, 45)')

    gmdl = ui.create_model_component('gauss2d', 'gmdl')
    bmdl = ui.create_model_component('scale2d', 'bmdl')
    ui.set_source(gmdl + bmdl)

    gmdl.xpos = 51
    gmdl.ypos = 49
    gmdl.ellip = 0.8
    gmdl.theta = 1.2
    gmdl.ampl = 10
    bmdl.c0 = 2

    forig = ui.get_filter()
    sorig = ui.calc_data_sum2d(forig)
    corig = ui.calc_stat()

    # sanity check
    assert forig == 'Circle(50,50,30)&!RotBox(30,30,10,5,45)|Ellipse(40,75,30,20,320)'
    assert sorig == pytest.approx(1861.978)
    assert corig == pytest.approx(7122.753868262877)

    compare(add_datadir_path(_canonical_img_filter_model))

    restore()
    fnew = ui.get_filter()
    snew = ui.calc_data_sum2d(fnew)
    cnew = ui.calc_stat()
    assert fnew == forig
    assert snew == pytest.approx(sorig)
    assert cnew == pytest.approx(corig)


@requires_data
@requires_fits
def test_restore_img_no_filter_model_psf(make_data_path, recwarn):
    """Can we save a PSF setup? See issue #1873"""

    # This is not a great image to fit any data to!
    ui.load_image(make_data_path('img.fits'))
    ui.load_psf("p0", make_data_path("psf_0.0_00_bin1.img"))

    ui.set_stat('cstat')
    ui.set_method('simplex')

    ui.set_source(ui.gauss2d.g1)
    g1.fwhm = 5
    g1.xpos = 49
    g1.ypos = 52
    g1.ellip = 0.4
    g1.theta = 1.2
    g1.ampl = 100

    sdata = ui.get_source_image().y
    mdata = ui.get_model_image().y
    assert mdata == pytest.approx(sdata)

    # Check we get the warning message (so we can also check when the
    # script is loaded).
    #
    wmsg = "Data Image does not have a pixel size. Sherpa will assume the pixel size is the same as the PSF"
    with pytest.warns(UserWarning, match=wmsg):
        ui.set_psf(p0)

    # This is just to check that the data is different in the expected
    # manner, e.g. the source has been convolved by the PSF so the
    # peak is different, but the overall signal is essentially the
    # same (since the edge effects can be neglected here).
    #
    maxval = 77.30269712891258
    sdata = ui.get_source_image().y
    mdata = ui.get_model_image().y
    assert sdata.max() == pytest.approx(100)
    assert mdata.max() == pytest.approx(maxval)
    assert mdata.sum() == pytest.approx(sdata.sum())  # the tolerance may need to be tweaked

    assert ui.get_filter() == ''

    compare(add_datadir_path(_canonical_img_no_filter_model_psf))

    with pytest.warns(UserWarning, match=wmsg):
        restore()

    sdata = ui.get_source_image().y
    mdata = ui.get_model_image().y
    assert sdata.max() == pytest.approx(100)
    assert mdata.max() == pytest.approx(maxval)
    assert mdata.sum() == pytest.approx(sdata.sum())

    assert ui.get_filter() == ''

    # There appears to be odd behavior with pytest and whether the
    # warnings are "cleared" here or not (however it is done, I have
    # stuck with pytest.warns above but similar issues are seen if we
    # use warnings.catch_warnings), so let's just clear the warnings
    # so conftest.capture_all_warnings doesn't error out.
    #
    recwarn.clear()


@requires_data
@requires_fits
def test_restore_table_model(make_data_path):
    """Note: this only sets the table model"""

    ui.load_table_model("tbl", make_data_path('test_rmfimg.fits'))
    tbl.ampl.set(10, min=0, max=20, frozen=True)

    compare(add_datadir_path(_canonical_table_model))

    restore()

    assert isinstance(tbl, TableModel)
    assert tbl.name == "tablemodel.tbl"
    assert tbl.ampl.val == 10
    assert tbl.ampl.min == 0
    assert tbl.ampl.max == 20
    assert tbl.ampl.frozen

    assert tbl.get_y().size == 1024 * 900
    assert tbl.get_y().max() == pytest.approx(0.2316664457321167)


@requires_xspec
@requires_data
@requires_fits
def test_restore_xstable_model(make_data_path):
    """Note: this only sets the table model

    This is a regression test as it currently does not do the right
    thing.

    """

    ui.load_xstable_model("tbl", make_data_path('testpcfabs.mod'))
    tbl.nh = 2
    tbl.fract = 0.2
    ui.freeze(tbl.fract)

    # Comparing the string version is a simlple way to check the
    # parameters.
    #
    orig = str(tbl)

    compare(add_datadir_path(_canonical_xstable_model))

    restore()

    assert isinstance(tbl, xspec.XSTableModel)
    assert str(tbl) == orig


def test_restore_pileup_model():
    """Note: this is not a realistic pileup-model case.

    It is assumed that other tests check the other parts we'd have
    for an actual pileup case.
    """

    ui.load_arrays(1, [1, 2, 3], [2, 5, 2], ui.DataPHA)
    ui.load_arrays(2, [1, 2, 3], [3, 4, 1], ui.DataPHA)
    pmod = ui.create_model_component("jdpileup", "pmod")
    ui.set_pileup_model(pmod)
    ui.set_pileup_model(2, pmod)

    compare(_canonical_pileup_model)

    restore()

    mod1 = ui.get_pileup_model(1)
    mod2 = ui.get_pileup_model(2)
    assert mod1.name == "jdpileup.pmod"
    assert mod2.name == "jdpileup.pmod"
    assert isinstance(mod1, JDPileup)


def test_restore_dataspace1d_int():
    """Can we restore a dataspace1d case?"""

    ui.dataspace1d(1, 10, step=2)
    ui.set_dep([2, 5, 6, 0, 2])
    ui.ignore(7, 8)

    fstr = "1.0000:7.0000,9.0000:11.0000"
    expected = [2, 5, 6, 2]
    assert ui.get_filter() == fstr
    assert ui.get_dep(filter=True) == pytest.approx(expected)

    compare(_canonical_dataspace1d_int)

    restore()

    assert ui.get_filter() == fstr
    assert ui.get_dep(filter=True) == pytest.approx(expected)


@requires_region
def test_restore_dataspace2d_img():
    """Can we restore a dataspace2d case?"""

    ui.set_stat("cash")

    ui.dataspace2d([4, 3])
    ui.set_dep([1, 2, 1, 0, 2, 3, 5, 2, 1, 0, 1, 1])
    ui.ignore2d("box(1.4,2.5,1.2,1.4)")

    fstr = "Field()&!Box(1.4,2.5,1.2,1.4)"
    expected = [1, 2, 1, 0, 5, 2, 1, 1]
    assert ui.get_filter() == fstr
    assert ui.get_dep(filter=True) == pytest.approx(expected)

    compare(_canonical_dataspace2d_img)

    restore()

    assert ui.get_filter() == fstr
    assert ui.get_dep(filter=True) == pytest.approx(expected)


def test_restore_load_arrays_simple():
    """Can we re-create a load_arrays/Data1D case

    The test_restore_dataspace1d_int call has checked we
    can restore a Data1DInt case.
    """

    ui.load_arrays("f", [-50, -20], [-2e4, 3e5])

    compare(_canonical_load_arrays_simple)

    restore()

    assert ui.list_data_ids() == ["f"]
    assert len(ui.get_indep("f")) == 1
    assert ui.get_indep("f")[0] == pytest.approx([-50, -20])
    assert ui.get_dep("f") == pytest.approx([-2e4, 3e5])


@requires_group
def test_restore_load_arrays_pha():
    """Can we re-create a load_arrays/DataPHA case?"""

    dset = ui.DataPHA("ex", [1, 2, 3, 4, 5], [12, 2, 1, 0, 1])
    dset.exposure = 100
    dset.backscal = 0.002
    dset.areascal = 0.001

    ui.set_data(dset)
    ui.group_counts(3, tabStops=numpy.asarray([0, 0, 0, 1, 0]))

    compare(_canonical_load_arrays_pha)

    restore()

    pha = ui.get_data()
    assert isinstance(pha, ui.DataPHA)
    assert pha.exposure == pytest.approx(100)
    assert pha.backscal == pytest.approx(0.002)
    assert pha.areascal == pytest.approx(0.001)

    assert pha.grouping == pytest.approx([1, 1, -1, 0, 1])
    assert pha.quality == pytest.approx([0, 0, 0, 0, 2])

    # This is corrected by thebin width and the areascal values
    assert ui.get_dep() == pytest.approx([12000, 1500, 0, 1000])

    assert pha.counts == pytest.approx([12, 2, 1, 0, 1])


def test_restore_load_arrays_data2d():
    """Can we re-create a load_arrays/Data2D case

    Note that test_restore_dataspace2d_img is a basic DataIMG test.

    """

    dset = ui.Data2D("ex", [1, 1, 2], [1, 2, 2], [3, 4, 5], None, [0.1, 0.1, 0.2])
    ui.set_data(dset)

    assert ui.get_dep() == pytest.approx([3, 4, 5])
    assert ui.get_staterror() == pytest.approx([0.1, 0.1, 0.2])

    compare(_canonical_load_arrays_data2d)

    restore()

    assert isinstance(ui.get_data(), ui.Data2D)
    assert ui.get_data().shape is None
    assert ui.get_dep() == pytest.approx([3, 4, 5])
    assert ui.get_staterror() == pytest.approx([0.1, 0.1, 0.2])


@requires_xspec
def test_canonical_xspec_hard_limit_min():
    "Can we save an XSPEC model with the hard limit extended: min"

    # Reset the optimiser parameters to make them easy to check,
    # even if they are un-usable.
    #
    for key in ui.get_method_opt().keys():
        ui.set_method_opt(key, 1)

    ui.create_model_component('xspowerlaw', 'mdl')
    mdl.phoindex.set(val=-5, hard_min=-5, frozen=True)
    mdl.norm.max = 100

    compare(_canonical_xspec_hard_limit_min)


@requires_xspec
def test_canonical_xspec_hard_limit_max():
    "Can we save an XSPEC model with the hard limit extended: max"

    # Reset the optimiser parameters to make them easy to check,
    # even if they are un-usable.
    #
    for key in ui.get_method_opt().keys():
        ui.set_method_opt(key, 1)

    ui.create_model_component('xspowerlaw', 'mdl')
    mdl.phoindex.set(val=15, hard_max=15, frozen=True)
    mdl.norm.max = 100

    compare(_canonical_xspec_hard_limit_max)


@requires_xspec
def test_restore_xspec_hard_limit_min():
    "Can the reload a XSPEC model with changed hard limits: min"

    # TODO: Why do we need to capture the return value when we
    # did not need to in test_canonical_xspec_hard_limit_min?
    #
    mdl = ui.create_model_component('xspowerlaw', 'mdl')
    mdl.phoindex.set(val=-5, hard_min=-5, frozen=True)
    mdl.norm.max = 100
    mdl = None

    restore()

    mdl = ui.get_model_component('mdl')
    assert mdl.name == 'xspowerlaw.mdl'
    assert mdl.PhoIndex.val == pytest.approx(-5)
    assert mdl.PhoIndex.min == pytest.approx(-5)
    assert mdl.PhoIndex.hard_min == pytest.approx(-5)
    assert mdl.PhoIndex.frozen

    assert mdl.norm.max == pytest.approx(100)


@requires_xspec
def test_restore_xspec_hard_limit_max():
    "Can the reload a XSPEC model with changed hard limits: max"

    # TODO: Why do we need to capture the return value when we
    # did not need to in test_canonical_xspec_hard_limit_min?
    #
    mdl = ui.create_model_component('xspowerlaw', 'mdl')
    mdl.phoindex.set(val=15, hard_max=15)
    mdl.norm.max = 100
    mdl = None

    restore()

    mdl = ui.get_model_component('mdl')
    assert mdl.name == 'xspowerlaw.mdl'
    assert mdl.PhoIndex.val == pytest.approx(15)
    assert mdl.PhoIndex.max == pytest.approx(15)
    assert mdl.PhoIndex.hard_max == pytest.approx(15)
    assert not mdl.PhoIndex.frozen

    assert mdl.norm.max == pytest.approx(100)


def test_restore_linker_parameter_with_function(clean_astro_ui):
    """Check we can save and restore a numpy ufunc"""

    m1 = ui.create_model_component("gauss1d", "m1")
    m2 = ui.create_model_component("scale1d", "m2")

    delta = numpy.sqrt((m2.c0 - 23)**2)
    m1.fwhm = 2 * numpy.exp(delta / 14)

    m2.c0 = 30

    expected = 2 * numpy.exp(0.5)

    # safety check
    assert m1.fwhm.link is not None
    assert m1.fwhm.val == pytest.approx(expected)

    restore()

    got = ui.list_model_components()
    assert len(got) == 2
    assert "m1" in got
    assert "m2" in got

    mm1 = ui.get_model_component("m1")
    assert mm1.name == "gauss1d.m1"
    assert mm1.fwhm.link is not None
    assert mm1.fwhm.val == pytest.approx(expected)


def test_link_par():
    """Check we can set up a simple parameter link

    Similar to test_restore_linker_parameter_with_function but
    the link is simpler (not a numpy ufunc).
    """

    m1 = ui.create_model_component("const1d", "m1")
    m2 = ui.create_model_component("const1d", "m2")
    sep = ui.create_model_component("const1d", "sep")

    sep.c0.set(100, min=100, max=100)
    m1.c0.set(200, min=0, max=1000)
    m2.c0.set(min=10, max=500)
    ui.link(m2.c0, m1.c0 + sep.c0)

    assert m2.c0.val == 300
    assert m2.c0.min == 10
    assert m2.c0.max == 500
    assert m2.c0.link is not None
    assert m2.c0.link.name == "m1.c0 + sep.c0"

    compare(_canonical_link_par)

    restore()

    assert m2.c0.val == 300
    assert m2.c0.min == 10
    assert m2.c0.max == 500
    assert m2.c0.link is not None
    assert m2.c0.link.name == "m1.c0 + sep.c0"


@pytest.mark.xfail
@requires_data
@requires_fits
def test_pha_full_model(make_data_path):
    """A basic set_full_model test.

    Related to issue #1050
    """

    ui.load_pha(make_data_path("3c273.pi"))
    ui.notice(1, 6)

    pl = ui.create_model_component("powlaw1d", "pl")
    con = ui.create_model_component("polynom1d", "con")

    pl.gamma = 1.7
    pl.ampl = 1e-4
    con.c0 = 0.001
    con.c1 = 0.002

    rsp = ui.get_response()
    model = rsp(pl) + con

    ui.set_full_model(model)

    # Can not add lo/hi arguments as calc_model_sum fails.  This is a
    # regression test.
    #
    expected = 1501.0484798733592
    assert ui.calc_model_sum() == pytest.approx(expected)

    with pytest.raises(IdentifierErr,
                       match=". You should use get_model instead.$"):
        ui.get_source()

    assert ui.get_model().name == "apply_rmf(apply_arf(38564.6089269 * powlaw1d.pl)) + polynom1d.con"

    compare(add_datadir_path(_canonical_pha_full_model))

    # check we get the expected error
    #
    with pytest.raises(NameError,
                       match="^name 'apply_rmf' is not defined$"):
        restore()

    assert ui.calc_model_sum() == pytest.approx(expected)  # XFAIL because restore failed

    with pytest.raises(IdentifierErr,
                       match=". You should use get_model instead.$"):
        ui.get_source()

    assert ui.get_model().name == "apply_rmf(apply_arf(38564.6089269 * powlaw1d.pl)) + polynom1d.con"


@requires_data
@requires_fits
def test_load_data(make_data_path):
    """Check load_data path for Data1D case with staterror"""

    # Unfortunately we need to care about the backend here.
    #
    from sherpa.astro import io

    ui.set_stat("chi2datavar")

    infile = make_data_path("data1.dat")
    ui.load_data(infile, ncols=3)

    expected_x = numpy.arange(0.5, 11)
    expected_y = [1.6454, 1.7236, 1.9472, 2.2348, 2.6187, 2.8642,
                  3.1263, 3.2073, 3.2852, 3.3092, 3.4496]
    expected_dy = [0.04114] * 11

    x = ui.get_indep()
    assert len(x) == 1
    assert x[0] == pytest.approx(expected_x)
    assert ui.get_dep() == pytest.approx(expected_y)
    assert ui.get_staterror() == pytest.approx(expected_dy)
    with pytest.raises(DataErr,
                       match="data set '1' does not specify systematic errors"):
        ui.get_syserror()

    if io.backend.name == "crates":
        term = "@@/data1.dat"
        idx = _canonical_load_data.find(term)
        expected_output = _canonical_load_data[:idx] + term + \
            "[opt colnames=none]" + \
            _canonical_load_data[idx + len(term):]
    else:
        expected_output = _canonical_load_data

    expected_output = add_datadir_path(expected_output)
    compare(expected_output)

    restore()

    x = ui.get_indep()
    assert len(x) == 1
    assert x[0] == pytest.approx(expected_x)
    assert ui.get_dep() == pytest.approx(expected_y)

    # The code has been adjusted to also read in the staterror
    # column, so issue #1876 doesn't hold here, but this is really
    # a work around.
    #
    assert ui.get_staterror() == pytest.approx(expected_dy)

    with pytest.raises(DataErr,
                       match="data set '1' does not specify systematic errors"):
        ui.get_syserror()


@requires_data
@requires_fits
def test_load_data_basic(make_data_path):
    """Check load_data path for Data1D case with no errors"""

    ui.set_stat("chi2")

    infile = make_data_path("double.dat")
    ui.load_data(infile)
    ui.ignore(hi=20)
    ui.ignore(lo=90)
    ui.ignore(25, 27)

    expected = '20.1000:24.9000,27.1000:89.9000'
    ptp = 0.9973281046
    assert ui.get_filter() == expected

    x = ui.get_indep()
    assert len(x) == 1
    assert len(x[0]) == 1000
    assert x[0][0] == pytest.approx(0)
    assert x[0][-1] == pytest.approx(99.9)
    assert numpy.ptp(ui.get_dep()) == pytest.approx(ptp)

    with pytest.raises(StatErr,
                       match="^If you select chi2 as the statistic, all datasets must provide a staterror column$"):
        ui.get_staterror()

    with pytest.raises(DataErr,
                       match="data set '1' does not specify systematic errors"):
        ui.get_syserror()

    expected_output = add_datadir_path(_canonical_load_data_basic)
    compare(expected_output)

    restore()

    expected = '20.1000:24.9000,27.1000:89.9000'
    assert ui.get_filter() == expected

    x = ui.get_indep()
    assert len(x) == 1
    assert len(x[0]) == 1000
    assert x[0][0] == pytest.approx(0)
    assert x[0][-1] == pytest.approx(99.9)
    assert numpy.ptp(ui.get_dep()) == pytest.approx(ptp)

    with pytest.raises(StatErr,
                       match="^If you select chi2 as the statistic, all datasets must provide a staterror column$"):
        ui.get_staterror()

    with pytest.raises(DataErr,
                       match="data set '1' does not specify systematic errors"):
        ui.get_syserror()


@requires_data
@requires_fits
def test_restore_pha_multiple_backgrounds(make_data_path):
    """Can we restore a grating dataset with multiple backgrounds?

    See issue #320

    This is a regression test so we can see as soon as things have
    changed, rather than marking it as xfail.

    """

    # Note: not including .gz for the file name
    ui.load_pha(make_data_path("3c120_meg_1.pha"))

    ui.set_analysis("wave")
    ui.notice(2, 12)

    def check_data():
        nelem = 2000
        assert len(ui.get_dep()) == 8192
        assert len(ui.get_dep(filter=True)) == nelem
        assert len(ui.get_dep(filter=True, bkg_id=1)) == nelem
        assert len(ui.get_dep(filter=True, bkg_id=2)) == nelem
        assert ui.get_dep(filter=True).sum() == 31911
        assert ui.get_dep(filter=True, bkg_id=1).sum() == 782
        assert ui.get_dep(filter=True, bkg_id=2).sum() == 516

        assert ui.get_data().get_filter(format="%.2f") == "2.00:12.00"

    check_data()

    expected_output = add_datadir_path(_canonical_pha_multiple_backgrounds)
    compare(expected_output)

    restore()

    assert ui.get_dep(filter=True, bkg_id=1).sum() == 31911  # TODO fix this
    # check_data()


@requires_data
@requires_fits
def test_restore_pha2(make_data_path):
    """Can we restore a pha2 file?

    See issue #1882.

    This is a regression test so we can see as soon as things have
    changed, rather than marking it as xfail.

    """

    # Note: not including .gz for the file name
    ui.load_pha(make_data_path("3c120_pha2"))
    ui.set_analysis("wave")
    ui.notice(2, 4)

    def check_data():
        assert ui.list_data_ids() == pytest.approx([1, 10, 11, 12, 2, 3, 4, 5, 6, 7, 8, 9])

        svals = {1: 193, 2: 306, 3: 4305, 4: 4579, 5: 269, 6: 164,
                 7: 386, 8: 90, 9: 3933, 10: 4372, 11: 98, 12: 350}
        for idx, total in svals.items():
            assert ui.get_dep(idx, filter=True).sum() == total

        b1vals = {1: 64, 2: 57, 3: 133, 4: 96, 5: 48, 6: 55,
                  7: 65, 8: 38, 9: 106, 10: 86, 11: 16, 12: 36}
        for idx, total in b1vals.items():
            assert ui.get_dep(idx, bkg_id=1, filter=True).sum() == total

        b2vals = {1: 74, 2: 46, 3: 108, 4: 65, 5: 37, 6: 44,
                  7: 44, 8: 46, 9: 67, 10: 69, 11: 16, 12: 34}
        for idx, total in b2vals.items():
            assert ui.get_dep(idx, bkg_id=2, filter=True).sum() == total

    check_data()

    expected_output = add_datadir_path(_canonical_pha2)
    compare(expected_output)

    restore()

    assert len(ui.list_data_ids()) == 23   # TODO fix this
    # check_data()


@requires_data
@requires_group
@requires_fits
@requires_xspec
def test_restore_pha_csc(make_data_path):
    """Can we restore a CSC file.

    This just checks that results of reading in. It doesn't
    check the actual text of the saved file.
    """

    ui.set_stat("cstat")
    ui.set_method("simplex")
    ui.load_pha("csc", make_data_path("acisf01575_001N001_r0085_pha3.fits.gz"))

    # We can check out some specialized cases, where we have a
    # different filter for the source and background data. Even though
    # we are using cstat we still bin (to 3 counts per bin).
    #
    # Note that this data is ungrouped.
    #
    # We check the sum of the grouping array as a reasonable way to
    # check the grouping is correct, without having to check all 1024
    # elements.
    #
    expected_sgrp = -110
    expected_bgrp = -457

    ui.notice(0.5, 7)
    ui.group_counts("csc", 3, tabStops=~ui.get_data("csc").get_mask())
    assert ui.get_grouping("csc").sum() == expected_sgrp

    ui.notice_id("csc", 0.4, 8, bkg_id=1)
    ui.group_counts("csc", 3, tabStops=~ui.get_bkg("csc").get_mask(), bkg_id=1)
    assert ui.get_grouping("csc", bkg_id=1).sum() == expected_bgrp
    assert ui.get_grouping("csc").sum() == expected_sgrp  # check it hasn't changed

    sfilter = "0.50:7.01"
    bfilter = "0.39:8.00"
    assert ui.get_data("csc").get_filter(format="%.2f") == sfilter
    assert ui.get_bkg("csc").get_filter(format="%.2f") == bfilter

    ui.set_bkg_source("csc", ui.powlaw1d.bpl)
    ui.set_source("csc", ui.xsphabs.gal * ui.powlaw1d.spl)

    bpl.gamma = 0.87
    bpl.ampl = 2.14e-6
    spl.gamma = 1.64
    spl.ampl = 3.9e-5
    gal.nh = 0.24

    expected_stat = 185.3380144554129
    assert ui.calc_stat("csc") == pytest.approx(expected_stat)

    expected_output = add_datadir_path(_canonical_pha_csc)
    compare(expected_output)

    restore()

    assert ui.get_grouping("csc", bkg_id=1).sum() == expected_bgrp
    assert ui.get_grouping("csc").sum() == expected_sgrp

    assert ui.get_data("csc").get_filter(format="%.2f") == sfilter
    assert ui.get_bkg("csc").get_filter(format="%.2f") == bfilter
    assert ui.calc_stat("csc") == pytest.approx(expected_stat)


def test_filter1d_excluded1():
    """Check what happens if all data filtered out.

    This test the mask=False case.
    """

    ui.load_arrays(1, [1, 2], [10, 20])
    ui.ignore()
    assert ui.get_data().mask is False

    restore()

    with pytest.raises(DataErr, match="^mask excludes all data$"):
        ui.get_dep(filter=True)


def test_filter1d_excluded2():
    """Check what happens if all data filtered out.

    This test the mask=[False, False] case.
    """

    ui.load_arrays(1, [1, 2], [10, 20])
    ui.ignore(hi=1)
    ui.ignore(lo=1)
    assert ui.get_data().mask is not False

    restore()

    with pytest.raises(DataErr, match="^mask excludes all data$"):
        ui.get_dep(filter=True)


def test_filter2d_excluded1():
    """Check what happens if all data filtered out.

    This test the mask=False case.
    """

    ui.load_arrays(1,
                   [1, 1, 1, 2, 2, 2],
                   [1, 2, 1, 2, 1, 2],
                   [1, 2, 3, 4, 5, 6],
                   (2, 3), ui.DataIMG)
    ui.ignore2d()
    assert ui.get_data().mask is False

    restore()

    with pytest.raises(DataErr, match="^mask excludes all data$"):
        ui.get_dep(filter=True)

    assert ui.get_data().shape[0] == 2
    assert ui.get_data().shape[1] == 3
    assert len(ui.get_data().shape) == 2


@requires_region
def test_filter2d_excluded2():
    """Check what happens if all data filtered out.

    This test the mask=[False, ...] case.
    """

    ui.load_arrays(1,
                   [1, 1, 1, 2, 2, 2],
                   [1, 2, 1, 2, 1, 2],
                   [1, 2, 3, 4, 5, 6],
                   (2, 3), ui.DataIMG)
    ui.ignore2d("Field()")
    assert ui.get_data().mask is not False

    restore()

    with pytest.raises(DataErr, match="^mask excludes all data$"):
        ui.get_dep(filter=True)


@requires_wcs
def test_fake_image_wcs():
    """Can we restore an image with WCS?"""

    # nx=2, ny=3
    #
    x1, x0 = numpy.mgrid[1:4, 1:3]
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = numpy.arange(6) * 10 + 10

    sky = WCS("physical", "LINEAR", [100, 200], [1, 1], [10, 10])
    eqpos = WCS("world", "WCS", [30, 50], [100, 200], [-0.01, 0.01])
    data = ui.DataIMG("faked", x0, x1, y,
                      shape=(2, 3), sky=sky, eqpos=eqpos)

    ui.set_data(data)

    restore()

    assert ui.get_coord() == "logical"

    sky = ui.get_data().sky
    eqpos = ui.get_data().eqpos
    assert isinstance(sky, WCS)
    assert isinstance(eqpos, WCS)

    assert sky.name == "physical"
    assert eqpos.name == "world"
    assert sky.type == "LINEAR"
    assert eqpos.type == "WCS"

    assert sky.crval == pytest.approx([100, 200])
    assert eqpos.cdelt == pytest.approx([-0.01, 0.01])

    g0, g1 = ui.get_indep()
    assert g0 == pytest.approx(x0)
    assert g1 == pytest.approx(x1)


@requires_wcs
@requires_region
def test_fake_image_wcs_coord():
    """Can we restore an image with WCS when coord has been changed?

    We also add a spatial filter to check that we aren't messing
    up x0/y0.
    """

    # nx=3, ny=2
    #
    x1, x0 = numpy.mgrid[1:3, 1:4]
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = numpy.arange(6) * 10 + 10

    sky = WCS("physical", "LINEAR", [100, 200], [1, 1], [10, 10])
    eqpos = WCS("world", "WCS", [30, 50], [100, 200], [-0.1, 0.1])
    data = ui.DataIMG("faked", x0, x1, y,
                      shape=(2, 3), sky=sky, eqpos=eqpos)

    ui.set_data(data)
    ui.set_coord("physical")

    ui.ignore2d("circle(110, 210, 6)")

    def check():

        assert ui.get_coord() == "physical"
        assert ui.get_filter() == "Field()&!Circle(110,210,6)"
        assert ui.get_dep(filter=True) == pytest.approx([10, 20, 30, 40, 60])
        a0, a1 = ui.get_indep()
        assert a0 == pytest.approx([100, 110, 120] * 2)
        assert a1 == pytest.approx([200, 200, 200, 210, 210, 210])

    check()

    restore()

    check()


@pytest.mark.parametrize("defid", ["foo", 2])
def test_store_default_id(defid):
    "Do we record the default id when it's been changed?"

    x = [3, 45]
    y = [2, 19]

    assert ui.get_default_id() == 1
    ui.set_default_id(defid)
    ui.load_arrays(defid, x, y)

    restore()
    assert ui.get_default_id() == defid
    assert ui.list_data_ids() == [defid]
    with pytest.raises(IdentifierErr,
                       match="^data set 1 has not been set$"):
        assert ui.get_data(1)

    d = ui.get_data()
    assert d.name == ""
    assert d.x == pytest.approx(x)
    assert d.y == pytest.approx(y)
