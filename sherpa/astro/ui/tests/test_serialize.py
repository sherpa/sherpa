#
#  Copyright (C) 2015, 2016, 2018-2021, 2023-2025
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

import numpy as np
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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

"""

# Change a few settings for statistic/method
_canonical_empty_stats = """import numpy
from sherpa.astro.ui import *


######### Set Statistic

set_stat("leastsq")


######### Set Fitting Method

set_method("neldermead")

set_method_opt("finalsimplex", 9)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
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
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
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


######### Load Background Data Sets

group(1, bkg_id=1)

######### Background Spectral Responses


######### Set Energy or Wave Units

set_analysis(1, quantity="energy", type="rate", factor=0)
subtract(1)

######### Filter Data

notice_id(1, "0.467200011015:9.869600296021")


######### Set Statistic

set_stat("chi2datavar")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP


######### Set Model Components and Parameters

create_model_component("xsapec", "src")
src.integrate = True

src.kT.default_val = 1.0
src.kT.default_min = 0.0080000000000000002  # doctest: +FLOAT_CMP
src.kT.default_max = 64.0
src.kT.val     = 1.0
src.kT.min     = 0.0080000000000000002  # doctest: +FLOAT_CMP
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
src.norm.default_max = 9.9999999999999998e+23  # doctest: +FLOAT_CMP
src.norm.val     = 1.0
src.norm.min     = 0.0
src.norm.max     = 9.9999999999999998e+23  # doctest: +FLOAT_CMP
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

_canonical_pha_basic_load = """import numpy
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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

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


######### Load Background Data Sets


######### Background grouping flags

set_grouping("grp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)

######### Background quality flags

set_quality("grp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)
group("grp", bkg_id=1)

######### Background Spectral Responses


######### Set Energy or Wave Units

set_analysis("grp", quantity="energy", type="rate", factor=0)
subtract("grp")

######### Filter Data

notice_id("grp", "0.467200011015:6.365600109100")


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP


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


######### Load Background Data Sets


######### Background grouping flags

set_grouping("bgrp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)

######### Background quality flags

set_quality("bgrp", val=numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], numpy.int16), bkg_id=1)
group("bgrp", bkg_id=1)

######### Background Spectral Responses


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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP


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
set_xsxset("APECROOT", "3.0.9")
set_xsxset("NEIAPECROOT", "3.0.9")
set_xsxset("NEIVERS", "3.0.4")
"""

_canonical_pha_load_bkg = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_arrays("x",
            [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0, 129.0, 130.0, 131.0, 132.0, 133.0, 134.0, 135.0, 136.0, 137.0, 138.0, 139.0, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 161.0, 162.0, 163.0, 164.0, 165.0, 166.0, 167.0, 168.0, 169.0, 170.0, 171.0, 172.0, 173.0, 174.0, 175.0, 176.0, 177.0, 178.0, 179.0, 180.0, 181.0, 182.0, 183.0, 184.0, 185.0, 186.0, 187.0, 188.0, 189.0, 190.0, 191.0, 192.0, 193.0, 194.0, 195.0, 196.0, 197.0, 198.0, 199.0, 200.0, 201.0, 202.0, 203.0, 204.0, 205.0, 206.0, 207.0, 208.0, 209.0, 210.0, 211.0, 212.0, 213.0, 214.0, 215.0, 216.0, 217.0, 218.0, 219.0, 220.0, 221.0, 222.0, 223.0, 224.0, 225.0, 226.0, 227.0, 228.0, 229.0, 230.0, 231.0, 232.0, 233.0, 234.0, 235.0, 236.0, 237.0, 238.0, 239.0, 240.0, 241.0, 242.0, 243.0, 244.0, 245.0, 246.0, 247.0, 248.0, 249.0, 250.0, 251.0, 252.0, 253.0, 254.0, 255.0, 256.0, 257.0, 258.0, 259.0, 260.0, 261.0, 262.0, 263.0, 264.0, 265.0, 266.0, 267.0, 268.0, 269.0, 270.0, 271.0, 272.0, 273.0, 274.0, 275.0, 276.0, 277.0, 278.0, 279.0, 280.0, 281.0, 282.0, 283.0, 284.0, 285.0, 286.0, 287.0, 288.0, 289.0, 290.0, 291.0, 292.0, 293.0, 294.0, 295.0, 296.0, 297.0, 298.0, 299.0, 300.0, 301.0, 302.0, 303.0, 304.0, 305.0, 306.0, 307.0, 308.0, 309.0, 310.0, 311.0, 312.0, 313.0, 314.0, 315.0, 316.0, 317.0, 318.0, 319.0, 320.0, 321.0, 322.0, 323.0, 324.0, 325.0, 326.0, 327.0, 328.0, 329.0, 330.0, 331.0, 332.0, 333.0, 334.0, 335.0, 336.0, 337.0, 338.0, 339.0, 340.0, 341.0, 342.0, 343.0, 344.0, 345.0, 346.0, 347.0, 348.0, 349.0, 350.0, 351.0, 352.0, 353.0, 354.0, 355.0, 356.0, 357.0, 358.0, 359.0, 360.0, 361.0, 362.0, 363.0, 364.0, 365.0, 366.0, 367.0, 368.0, 369.0, 370.0, 371.0, 372.0, 373.0, 374.0, 375.0, 376.0, 377.0, 378.0, 379.0, 380.0, 381.0, 382.0, 383.0, 384.0, 385.0, 386.0, 387.0, 388.0, 389.0, 390.0, 391.0, 392.0, 393.0, 394.0, 395.0, 396.0, 397.0, 398.0, 399.0, 400.0, 401.0, 402.0, 403.0, 404.0, 405.0, 406.0, 407.0, 408.0, 409.0, 410.0, 411.0, 412.0, 413.0, 414.0, 415.0, 416.0, 417.0, 418.0, 419.0, 420.0, 421.0, 422.0, 423.0, 424.0, 425.0, 426.0, 427.0, 428.0, 429.0, 430.0, 431.0, 432.0, 433.0, 434.0, 435.0, 436.0, 437.0, 438.0, 439.0, 440.0, 441.0, 442.0, 443.0, 444.0, 445.0, 446.0, 447.0, 448.0, 449.0, 450.0, 451.0, 452.0, 453.0, 454.0, 455.0, 456.0, 457.0, 458.0, 459.0, 460.0, 461.0, 462.0, 463.0, 464.0, 465.0, 466.0, 467.0, 468.0, 469.0, 470.0, 471.0, 472.0, 473.0, 474.0, 475.0, 476.0, 477.0, 478.0, 479.0, 480.0, 481.0, 482.0, 483.0, 484.0, 485.0, 486.0, 487.0, 488.0, 489.0, 490.0, 491.0, 492.0, 493.0, 494.0, 495.0, 496.0, 497.0, 498.0, 499.0, 500.0, 501.0, 502.0, 503.0, 504.0, 505.0, 506.0, 507.0, 508.0, 509.0, 510.0, 511.0, 512.0, 513.0, 514.0, 515.0, 516.0, 517.0, 518.0, 519.0, 520.0, 521.0, 522.0, 523.0, 524.0, 525.0, 526.0, 527.0, 528.0, 529.0, 530.0, 531.0, 532.0, 533.0, 534.0, 535.0, 536.0, 537.0, 538.0, 539.0, 540.0, 541.0, 542.0, 543.0, 544.0, 545.0, 546.0, 547.0, 548.0, 549.0, 550.0, 551.0, 552.0, 553.0, 554.0, 555.0, 556.0, 557.0, 558.0, 559.0, 560.0, 561.0, 562.0, 563.0, 564.0, 565.0, 566.0, 567.0, 568.0, 569.0, 570.0, 571.0, 572.0, 573.0, 574.0, 575.0, 576.0, 577.0, 578.0, 579.0, 580.0, 581.0, 582.0, 583.0, 584.0, 585.0, 586.0, 587.0, 588.0, 589.0, 590.0, 591.0, 592.0, 593.0, 594.0, 595.0, 596.0, 597.0, 598.0, 599.0, 600.0, 601.0, 602.0, 603.0, 604.0, 605.0, 606.0, 607.0, 608.0, 609.0, 610.0, 611.0, 612.0, 613.0, 614.0, 615.0, 616.0, 617.0, 618.0, 619.0, 620.0, 621.0, 622.0, 623.0, 624.0, 625.0, 626.0, 627.0, 628.0, 629.0, 630.0, 631.0, 632.0, 633.0, 634.0, 635.0, 636.0, 637.0, 638.0, 639.0, 640.0, 641.0, 642.0, 643.0, 644.0, 645.0, 646.0, 647.0, 648.0, 649.0, 650.0, 651.0, 652.0, 653.0, 654.0, 655.0, 656.0, 657.0, 658.0, 659.0, 660.0, 661.0, 662.0, 663.0, 664.0, 665.0, 666.0, 667.0, 668.0, 669.0, 670.0, 671.0, 672.0, 673.0, 674.0, 675.0, 676.0, 677.0, 678.0, 679.0, 680.0, 681.0, 682.0, 683.0, 684.0, 685.0, 686.0, 687.0, 688.0, 689.0, 690.0, 691.0, 692.0, 693.0, 694.0, 695.0, 696.0, 697.0, 698.0, 699.0, 700.0, 701.0, 702.0, 703.0, 704.0, 705.0, 706.0, 707.0, 708.0, 709.0, 710.0, 711.0, 712.0, 713.0, 714.0, 715.0, 716.0, 717.0, 718.0, 719.0, 720.0, 721.0, 722.0, 723.0, 724.0, 725.0, 726.0, 727.0, 728.0, 729.0, 730.0, 731.0, 732.0, 733.0, 734.0, 735.0, 736.0, 737.0, 738.0, 739.0, 740.0, 741.0, 742.0, 743.0, 744.0, 745.0, 746.0, 747.0, 748.0, 749.0, 750.0, 751.0, 752.0, 753.0, 754.0, 755.0, 756.0, 757.0, 758.0, 759.0, 760.0, 761.0, 762.0, 763.0, 764.0, 765.0, 766.0, 767.0, 768.0, 769.0, 770.0, 771.0, 772.0, 773.0, 774.0, 775.0, 776.0, 777.0, 778.0, 779.0, 780.0, 781.0, 782.0, 783.0, 784.0, 785.0, 786.0, 787.0, 788.0, 789.0, 790.0, 791.0, 792.0, 793.0, 794.0, 795.0, 796.0, 797.0, 798.0, 799.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            DataPHA)

######### Load Background Data Sets

load_bkg("x", "@@/MNLup_2138_0670580101_EMOS1_S001_specbg.fits", bkg_id="foo", use_errors=True)

######### Background Spectral Responses

load_arf("x", "@@/MNLup_2138_0670580101_EMOS1_S001_spec.arf", resp_id=1, bkg_id="foo")
load_rmf("x", "@@/MNLup_2138_0670580101_EMOS1_S001_spec.rmf", resp_id=1, bkg_id="foo")

######### Set Energy or Wave Units

set_analysis("x", quantity="channel", type="rate", factor=0)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("gridsearch")

set_method_opt("ftol", 1.1920928955078125e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("method", None)
set_method_opt("num", 16)
set_method_opt("numcores", 1)
set_method_opt("sequence", None)
set_method_opt("verbose", 0)

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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP


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
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
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
sin_model.ampl.default_min = 1.0000000000000001e-05  # doctest: +FLOAT_CMP
sin_model.ampl.default_max = 3.4028234663852886e+38
sin_model.ampl.val     = 1.0
sin_model.ampl.min     = 1.0000000000000001e-05  # doctest: +FLOAT_CMP
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
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
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
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
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

gmdl.ellip.default_val = 0.80000000000000004  # doctest: +FLOAT_CMP
gmdl.ellip.default_min = 0.0
gmdl.ellip.default_max = 0.999
gmdl.ellip.val     = 0.80000000000000004  # doctest: +FLOAT_CMP
gmdl.ellip.min     = 0.0
gmdl.ellip.max     = 0.999
gmdl.ellip.units   = ""
gmdl.ellip.frozen  = True

gmdl.theta.default_val = 1.2
gmdl.theta.default_min = -6.2831853071795862  # doctest: +FLOAT_CMP
gmdl.theta.default_max = 6.2831853071795862  # doctest: +FLOAT_CMP
gmdl.theta.val     = 1.2
gmdl.theta.min     = -6.2831853071795862  # doctest: +FLOAT_CMP
gmdl.theta.max     = 6.2831853071795862  # doctest: +FLOAT_CMP
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
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
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

g1.ellip.default_val = 0.40000000000000002  # doctest: +FLOAT_CMP
g1.ellip.default_min = 0.0
g1.ellip.default_max = 0.999
g1.ellip.val     = 0.40000000000000002  # doctest: +FLOAT_CMP
g1.ellip.min     = 0.0
g1.ellip.max     = 0.999
g1.ellip.units   = ""
g1.ellip.frozen  = True

g1.theta.default_val = 1.2
g1.theta.default_min = -6.2831853071795862  # doctest: +FLOAT_CMP
g1.theta.default_max = 6.2831853071795862  # doctest: +FLOAT_CMP
g1.theta.val     = 1.2
g1.theta.min     = -6.2831853071795862  # doctest: +FLOAT_CMP
g1.theta.max     = 6.2831853071795862  # doctest: +FLOAT_CMP
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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP


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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP


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

tbl.fract.default_val = 0.20000000000000001  # doctest: +FLOAT_CMP
tbl.fract.default_min = 0.0
tbl.fract.default_max = 1.0
tbl.fract.val     = 0.20000000000000001  # doctest: +FLOAT_CMP
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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP


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

pmod.f.default_val = 0.94999999999999996  # doctest: +FLOAT_CMP
pmod.f.default_min = 0.90000000000000002  # doctest: +FLOAT_CMP
pmod.f.default_max = 1.0
pmod.f.val     = 0.94999999999999996  # doctest: +FLOAT_CMP
pmod.f.min     = 0.90000000000000002  # doctest: +FLOAT_CMP
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

pmod.ftime.default_val = 3.2410000000000001  # doctest: +FLOAT_CMP
pmod.ftime.default_min = 1.1754943508222875e-38
pmod.ftime.default_max = 5.0
pmod.ftime.val     = 3.2410000000000001  # doctest: +FLOAT_CMP
pmod.ftime.min     = 1.1754943508222875e-38
pmod.ftime.max     = 5.0
pmod.ftime.units   = "sec"
pmod.ftime.frozen  = True

pmod.fracexp.default_val = 0.98699999999999999  # doctest: +FLOAT_CMP
pmod.fracexp.default_min = 0.0
pmod.fracexp.default_max = 1.0
pmod.fracexp.val     = 0.98699999999999999  # doctest: +FLOAT_CMP
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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

"""

_canonical_load_arrays_pha_response = """import numpy
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

######### Data Spectral Responses

load_arf(1, "test-arf", resp_id=1)
load_rmf(1, "delta-rmf", resp_id=1)

######### Set Energy or Wave Units

set_analysis(1, quantity="energy", type="rate", factor=0)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

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
mdl.norm.default_max = 9.9999999999999998e+23  # doctest: +FLOAT_CMP
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
mdl.norm.default_max = 9.9999999999999998e+23  # doctest: +FLOAT_CMP
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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP


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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

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

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

"""

_canonical_pha_multiple_backgrounds = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha(1, "@@/3c120_meg_1.pha")

######### Data Spectral Responses


######### Load Background Data Sets


######### Background Spectral Responses


######### Background Spectral Responses


######### Set Energy or Wave Units

set_analysis(1, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(1, "2.000000050569:11.999999841898")
notice_id(1, bkg_id=1)
notice_id(1, "2.000000050569:11.999999841898", bkg_id=1)
notice_id(1, bkg_id=2)
notice_id(1, "2.000000050569:11.999999841898", bkg_id=2)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

"""

_canonical_pha2 = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

# Load PHA2 into: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
load_pha(1, "@@/3c120_pha2")

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(1, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(1, "3793:6192")
notice_id(1, bkg_id=1)
notice_id(1, "3793:6192", bkg_id=1)
notice_id(1, bkg_id=2)
notice_id(1, "3793:6192", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(10, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(10, "7593:7992")
notice_id(10, bkg_id=1)
notice_id(10, "7593:7992", bkg_id=1)
notice_id(10, bkg_id=2)
notice_id(10, "7593:7992", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(11, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(11, "6793:7592")
notice_id(11, bkg_id=1)
notice_id(11, "6793:7592", bkg_id=1)
notice_id(11, bkg_id=2)
notice_id(11, "6793:7592", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(12, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(12, "5993:7192")
notice_id(12, bkg_id=1)
notice_id(12, "5993:7192", bkg_id=1)
notice_id(12, bkg_id=2)
notice_id(12, "5993:7192", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(2, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(2, "5393:6992")
notice_id(2, bkg_id=1)
notice_id(2, "5393:6992", bkg_id=1)
notice_id(2, bkg_id=2)
notice_id(2, "5393:6992", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(3, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(3, "6993:7792")
notice_id(3, bkg_id=1)
notice_id(3, "6993:7792", bkg_id=1)
notice_id(3, bkg_id=2)
notice_id(3, "6993:7792", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(4, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(4, "6993:7792")
notice_id(4, bkg_id=1)
notice_id(4, "6993:7792", bkg_id=1)
notice_id(4, bkg_id=2)
notice_id(4, "6993:7792", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(5, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(5, "5393:6992")
notice_id(5, bkg_id=1)
notice_id(5, "5393:6992", bkg_id=1)
notice_id(5, bkg_id=2)
notice_id(5, "5393:6992", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(6, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(6, "3793:6192")
notice_id(6, bkg_id=1)
notice_id(6, "3793:6192", bkg_id=1)
notice_id(6, bkg_id=2)
notice_id(6, "3793:6192", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(7, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(7, "5993:7192")
notice_id(7, bkg_id=1)
notice_id(7, "5993:7192", bkg_id=1)
notice_id(7, bkg_id=2)
notice_id(7, "5993:7192", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(8, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(8, "6793:7592")
notice_id(8, bkg_id=1)
notice_id(8, "6793:7592", bkg_id=1)
notice_id(8, bkg_id=2)
notice_id(8, "6793:7592", bkg_id=2)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(9, quantity="channel", type="rate", factor=0)

######### Filter Data

notice_id(9, "7593:7992")
notice_id(9, bkg_id=1)
notice_id(9, "7593:7992", bkg_id=1)
notice_id(9, bkg_id=2)
notice_id(9, "7593:7992", bkg_id=2)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

"""

_canonical_pha2_delete = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

# Load PHA2 into: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
load_pha(1, "@@/3c120_pha2")
delete_data(1)
delete_data(2)
delete_data(5)
delete_data(6)
delete_data(7)
delete_data(8)
delete_data(11)
delete_data(12)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(10, quantity="channel", type="rate", factor=0)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(3, quantity="channel", type="rate", factor=0)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(4, quantity="channel", type="rate", factor=0)

######### Load Background Data Sets


######### Set Energy or Wave Units

set_analysis(9, quantity="channel", type="rate", factor=0)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.19209289551e-07)  # doctest: +FLOAT_CMP

"""

_canonical_pha_hetg = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha(1, "@@/3c120_heg_-1.pha")

######### Data Spectral Responses


######### Load Background Data Sets


######### Background Spectral Responses


######### Background Spectral Responses


######### Set Energy or Wave Units

set_analysis(1, quantity="wavelength", type="rate", factor=0)

######### Filter Data

notice_id(1, "5.000000222571:10.000000445141")
notice_id(1, bkg_id=1)
notice_id(1, "5.000000222571:10.000000445141", bkg_id=1)
notice_id(1, bkg_id=2)
notice_id(1, "5.000000222571:10.000000445141", bkg_id=2)


######### Set Statistic

set_stat("chi2gehrels")


######### Set Fitting Method

set_method("levmar")

set_method_opt("epsfcn", 1.1920928955078125e-07)  # doctest: +FLOAT_CMP
set_method_opt("factor", 100.0)
set_method_opt("ftol", 1.1920928955078125e-07)  # doctest: +FLOAT_CMP
set_method_opt("gtol", 1.1920928955078125e-07)  # doctest: +FLOAT_CMP
set_method_opt("maxfev", None)
set_method_opt("numcores", 1)
set_method_opt("verbose", 0)
set_method_opt("xtol", 1.1920928955078125e-07)  # doctest: +FLOAT_CMP

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


######### Load Background Data Sets


######### Background grouping flags

set_grouping("csc", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)

######### Background quality flags

set_quality("csc", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16), bkg_id=1)
group("csc", bkg_id=1)

######### Background Spectral Responses


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
set_method_opt("ftol", 1.19209289551e-07)  # doctest: +FLOAT_CMP
set_method_opt("initsimplex", 0)
set_method_opt("iquad", 1)
set_method_opt("maxfev", None)
set_method_opt("reflect", True)
set_method_opt("step", None)
set_method_opt("verbose", 0)


######### Set Model Components and Parameters

create_model_component("powlaw1d", "spl")
spl.integrate = True

spl.gamma.default_val = 1.6399999999999999  # doctest: +FLOAT_CMP
spl.gamma.default_min = -10.0
spl.gamma.default_max = 10.0
spl.gamma.val     = 1.6399999999999999  # doctest: +FLOAT_CMP
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

spl.ampl.default_val = 3.8999999999999999e-05  # doctest: +FLOAT_CMP
spl.ampl.default_min = 0.0
spl.ampl.default_max = 3.4028234663852886e+38
spl.ampl.val     = 3.8999999999999999e-05  # doctest: +FLOAT_CMP
spl.ampl.min     = 0.0
spl.ampl.max     = 3.4028234663852886e+38
spl.ampl.units   = ""
spl.ampl.frozen  = False

create_model_component("xsphabs", "gal")
gal.integrate = True

gal.nH.default_val = 0.23999999999999999  # doctest: +FLOAT_CMP
gal.nH.default_min = 0.0
gal.nH.default_max = 1000000.0
gal.nH.val     = 0.23999999999999999  # doctest: +FLOAT_CMP
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

bpl.ampl.default_val = 2.1399999999999998e-06  # doctest: +FLOAT_CMP
bpl.ampl.default_min = 0.0
bpl.ampl.default_max = 3.4028234663852886e+38
bpl.ampl.val     = 2.1399999999999998e-06  # doctest: +FLOAT_CMP
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
set_xsxset("APECROOT", "3.0.9")
set_xsxset("NEIAPECROOT", "3.0.9")
set_xsxset("NEIVERS", "3.0.4")
"""

    _canonical_pha_basic += _canonical_extra
    _canonical_pha_grouped += _canonical_extra
    _canonical_xstable_model += _canonical_extra
    _canonical_xspec_hard_limit_min += _canonical_extra
    _canonical_xspec_hard_limit_max += _canonical_extra
    _canonical_pha_csc += _canonical_extra

    del _canonical_extra


@pytest.fixture(autouse=True)
def setup(hide_logging, clean_astro_ui):
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


def compare(check_str, expected, **kwargs):
    """Run save_all and check the output (saved to a
    StringIO object) to the string value expected.

    Since check_str is sent in, lines in expected can end with '#
    doctest: +FLOAT_CMP' to allow for approximate numeric checks (but
    limited by the check_str fixture).

    """
    output = StringIO()
    ui.save_all(output, **kwargs)
    output = output.getvalue()

    # check the output is a valid Python program.
    # this check does not guard against potential issues,
    # but ensures that the program can compile.
    #
    compileit(output)
    check_str(output, expected.split("\n"))


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


def test_canonical_empty(check_str):
    "Contents of empty state are as expected"
    compare(check_str, _canonical_empty)


def test_canonical_empty_outfile(tmp_path, check_str):
    "Can read in a save file"

    outfile = tmp_path / "save.sherpa"
    ui.save_all(str(outfile))
    output = outfile.read_text()
    check_str(output, _canonical_empty.split("\n"))


def test_canonical_empty_stats(check_str):
    "Change several settings but load no data"

    ui.set_stat('leastsq')

    ui.set_method('simplex')
    ui.set_method_opt('maxfev', 5000)
    ui.set_method_opt('verbose', 1)

    compare(check_str, _canonical_empty_stats)


def test_canonical_empty_iterstat(check_str):
    "Check iterated-fit setting"

    ui.set_stat('leastsq')

    ui.set_method('simplex')
    ui.set_method_opt('maxfev', 5000)
    ui.set_method_opt('verbose', 1)

    ui.set_iter_method('sigmarej')
    ui.set_iter_method_opt('grow', 1)

    compare(check_str, _canonical_empty_iterstat)


@requires_data
@requires_xspec
@requires_fits
def test_canonical_pha_basic(make_data_path, check_str):

    _, canonical = setup_pha_basic(make_data_path)
    compare(check_str, canonical)


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
def test_canonical_pha_load(make_data_path, check_str):
    """Do we explicitly load the ancillary files?"""

    # This is setup_pha_basic but without the source model.
    #
    fname = make_data_path('3c273.pi')
    ui.load_pha(1, fname)
    ui.subtract()
    ui.set_stat('chi2datavar')
    ui.notice(0.5, 7)

    canonical = add_datadir_path(_canonical_pha_basic_load)
    compare(check_str, canonical, auto_load=False)


@requires_data
@requires_xspec
@requires_fits
def test_canonical_pha_basic_errors(make_data_path, check_str):
    """Do we acknowledge the use_errors setting in load_pha?"""

    fname = make_data_path('3c273.pi')
    ui.set_default_id("sos")
    ui.load_pha(fname, use_errors=True)
    ui.notice(0.5, 6)
    ui.subtract()

    # Approximate fit just so the statistic isn't too large.
    gal = ui.create_model_component("xsphabs", "gal")
    pl = ui.create_model_component("powlaw1d", "pl")
    ui.set_source(gal * pl)
    gal.nh = 0.04
    pl.gamma = 2.03
    pl.ampl = 1.96e-4

    # The assumption here is that the error values are different to
    # those from use_errors=False. This is not explicitly tested in
    # this set of tests.
    #
    statval = ui.calc_stat()
    evals = ui.get_staterror()

    restore()

    assert ui.get_default_id() == "sos"

    # Check the "error values" from the file are being used.
    #
    assert ui.calc_stat() == pytest.approx(statval)
    assert ui.get_staterror() == pytest.approx(evals)


@requires_data
@requires_xspec
@requires_fits
@requires_group
def test_canonical_pha_grouped(make_data_path, check_str):

    _, _, canonical = setup_pha_grouped(make_data_path)
    compare(check_str, canonical)


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
    assert g.dtype == np.int16
    assert q.dtype == np.int16

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
def test_canonical_pha_back(make_data_path, check_str):

    _, _, canonical = setup_pha_back(make_data_path)
    compare(check_str, canonical)


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
    # assert g.dtype == np.int16
    # assert q.dtype == np.int16

    # TODO set up correct grouping bins...
    # nchan = ui.get_data('bgrp').channel.size
    # assert_array_equal(g, np.ones(nchan), err_msg='src grouping')
    # assert_array_equal(q, np.zeros(nchan), err_msg='src quality')

    bg = ui.get_grouping('bgrp', bkg_id=1)
    bq = ui.get_quality('bgrp', bkg_id=1)
    assert bg.dtype == np.int16
    assert bq.dtype == np.int16

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
def test_canonical_pha_load_bkg(make_data_path, check_str):
    """Check an explicit load_bkg call.

    This is not a realistic situation, but does test the
    load_bkg call.

    """

    ui.set_method("gridsearch")

    ui.dataspace1d(0, 799, id="x", dstype=ui.DataPHA)

    # NOTE: use_errors setting is currently ignored; for now this is
    # treated as a regression test (i.e. _canonical_pha_load_bkg will
    # need to be updated once this is fixed). Issue #2383.
    #
    BASENAME = "MNLup_2138_0670580101_EMOS1_S001_specbg.fits"
    bgfile = make_data_path(BASENAME)
    ui.load_bkg("x", bgfile, bkg_id="foo", use_errors=True)

    # Load a response, but for the background only. These responses
    # cause warning messages to get displayed, so hide them.
    #
    RMFNAME = "MNLup_2138_0670580101_EMOS1_S001_spec.rmf"
    ARFNAME = "MNLup_2138_0670580101_EMOS1_S001_spec.arf"
    rmffile = make_data_path(RMFNAME)
    arffile = make_data_path(ARFNAME)
    with warnings.catch_warnings(record=True):
        ui.load_rmf("x", rmffile, bkg_id="foo")
        ui.load_arf("x", arffile, bkg_id="foo")

    # The name of the dataspace1d call is not saved, hence the
    # version-specific check.
    #
    def check_data(name: str) -> None:
        assert ui.list_data_ids() == ["x"]
        data = ui.get_data("x")
        assert data.name == name
        assert data.channel == pytest.approx(np.arange(0, 800))
        assert data.counts == pytest.approx([0] * 800)

        assert ui.list_bkg_ids("x") == ["foo"]
        bkg = ui.get_bkg("x", bkg_id="foo")
        assert bkg.name.endswith(f"/{BASENAME}")
        assert bkg.channel == pytest.approx(np.arange(0, 800))
        # compare to the sum of the counds calculated with dmstat
        # on the background file
        assert bkg.counts.sum() == pytest.approx(2880)
        assert bkg.counts.max() == pytest.approx(42)

        rmf = ui.get_rmf("x", bkg_id="foo")
        assert rmf.name.endswith(f"/{RMFNAME}")

        arf = ui.get_arf("x", bkg_id="foo")
        assert arf.name.endswith(f"/{ARFNAME}")

        # The source has no response
        #
        with pytest.raises(IdentifierErr):
            ui.get_rmf("x")

        with pytest.raises(IdentifierErr):
            ui.get_arf("x")

    check_data(name="dataspace1d")

    expected_output = add_datadir_path(_canonical_pha_load_bkg)
    compare(check_str, expected_output)

    # Need to avoid the warning messages from the ARF/RMF here too.
    #
    with warnings.catch_warnings(record=True):
        restore()

    check_data(name="")


@requires_data
@requires_fits
def test_pha_no_response(make_data_path, check_str):
    """Check when there's no ARF / RMF or grouping, and switch to counts"""

    infile = make_data_path("source1.pi")
    ui.load_pha("x", infile)
    ui.notice(40, 400)
    ui.ignore(300, 319)
    ui.set_analysis("x", quantity="channel", type="counts")

    assert ui.get_dep("x", filter=True).size == 341
    assert ui.get_dep("x", filter=True).sum() == 2904
    assert ui.get_dep("x").sum() == 3328

    compare(check_str, add_datadir_path(_canonical_pha_no_response))

    restore()

    assert ui.get_dep("x", filter=True).size == 341
    assert ui.get_dep("x", filter=True).sum() == 2904
    assert ui.get_dep("x").sum() == 3328


def test_canonical_usermodel(check_str):
    "Can we save a usermodel?"
    setup_usermodel()
    compare(check_str, _canonical_usermodel)


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
def test_restore_img_no_filter_no_model(make_data_path, check_str):
    """Check issue #437"""

    ui.load_image(make_data_path('img.fits'))
    ui.set_stat('cstat')
    ui.set_method('simplex')

    sorig = ui.calc_data_sum2d()

    # sanity check
    assert sorig == pytest.approx(5041.44)
    assert ui.get_filter() == ''

    compare(check_str, add_datadir_path(_canonical_img_no_filter_no_model))

    restore()
    snew = ui.calc_data_sum2d()
    assert snew == pytest.approx(sorig)
    assert ui.get_filter() == ''


@requires_data
@requires_fits
@requires_region
def test_restore_img_filter_model(make_data_path, check_str):
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

    compare(check_str, add_datadir_path(_canonical_img_filter_model))

    restore()
    fnew = ui.get_filter()
    snew = ui.calc_data_sum2d(fnew)
    cnew = ui.calc_stat()
    assert fnew == forig
    assert snew == pytest.approx(sorig)
    assert cnew == pytest.approx(corig)


@requires_data
@requires_fits
def test_restore_img_no_filter_model_psf(make_data_path, recwarn, check_str):
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

    compare(check_str, add_datadir_path(_canonical_img_no_filter_model_psf))

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
def test_restore_table_model(make_data_path, check_str):
    """Note: this only sets the table model"""

    ui.load_table_model("tbl", make_data_path('test_rmfimg.fits'))
    tbl.ampl.set(10, min=0, max=20, frozen=True)

    compare(check_str, add_datadir_path(_canonical_table_model))

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
def test_restore_xstable_model(make_data_path, check_str):
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

    compare(check_str, add_datadir_path(_canonical_xstable_model))

    restore()

    assert isinstance(tbl, xspec.XSTableModel)
    assert str(tbl) == orig


def test_restore_pileup_model(check_str):
    """Note: this is not a realistic pileup-model case.

    It is assumed that other tests check the other parts we'd have
    for an actual pileup case.
    """

    ui.load_arrays(1, [1, 2, 3], [2, 5, 2], ui.DataPHA)
    ui.load_arrays(2, [1, 2, 3], [3, 4, 1], ui.DataPHA)
    pmod = ui.create_model_component("jdpileup", "pmod")
    ui.set_pileup_model(pmod)
    ui.set_pileup_model(2, pmod)

    compare(check_str, _canonical_pileup_model)

    restore()

    mod1 = ui.get_pileup_model(1)
    mod2 = ui.get_pileup_model(2)
    assert mod1.name == "jdpileup.pmod"
    assert mod2.name == "jdpileup.pmod"
    assert isinstance(mod1, JDPileup)


def test_restore_dataspace1d_int(check_str):
    """Can we restore a dataspace1d case?"""

    ui.dataspace1d(1, 10, step=2)
    ui.set_dep([2, 5, 6, 0, 2])
    ui.ignore(7, 8)

    fstr = "1.0000:7.0000,9.0000:11.0000"
    expected = [2, 5, 6, 2]
    assert ui.get_filter() == fstr
    assert ui.get_dep(filter=True) == pytest.approx(expected)

    compare(check_str, _canonical_dataspace1d_int)

    restore()

    assert ui.get_filter() == fstr
    assert ui.get_dep(filter=True) == pytest.approx(expected)


@requires_region
def test_restore_dataspace2d_img(check_str):
    """Can we restore a dataspace2d case?"""

    ui.set_stat("cash")

    ui.dataspace2d([4, 3])
    ui.set_dep([1, 2, 1, 0, 2, 3, 5, 2, 1, 0, 1, 1])
    ui.ignore2d("box(1.4,2.5,1.2,1.4)")

    fstr = "Field()&!Box(1.4,2.5,1.2,1.4)"
    expected = [1, 2, 1, 0, 5, 2, 1, 1]
    assert ui.get_filter() == fstr
    assert ui.get_dep(filter=True) == pytest.approx(expected)

    compare(check_str, _canonical_dataspace2d_img)

    restore()

    assert ui.get_filter() == fstr
    assert ui.get_dep(filter=True) == pytest.approx(expected)


def test_restore_load_arrays_simple(check_str):
    """Can we re-create a load_arrays/Data1D case

    The test_restore_dataspace1d_int call has checked we
    can restore a Data1DInt case.
    """

    ui.load_arrays("f", [-50, -20], [-2e4, 3e5])

    compare(check_str, _canonical_load_arrays_simple)

    restore()

    assert ui.list_data_ids() == ["f"]
    assert len(ui.get_indep("f")) == 1
    assert ui.get_indep("f")[0] == pytest.approx([-50, -20])
    assert ui.get_dep("f") == pytest.approx([-2e4, 3e5])


@requires_group
def test_restore_load_arrays_pha(check_str):
    """Can we re-create a load_arrays/DataPHA case?"""

    dset = ui.DataPHA("ex", [1, 2, 3, 4, 5], [12, 2, 1, 0, 1])
    dset.exposure = 100
    dset.backscal = 0.002
    dset.areascal = 0.001

    ui.set_data(dset)
    ui.group_counts(3, tabStops=np.asarray([0, 0, 0, 1, 0]))

    compare(check_str, _canonical_load_arrays_pha)

    restore()

    pha = ui.get_data()
    assert isinstance(pha, ui.DataPHA)
    assert pha.exposure == pytest.approx(100)
    assert pha.backscal == pytest.approx(0.002)
    assert pha.areascal == pytest.approx(0.001)

    assert pha.grouping == pytest.approx([1, 1, -1, 0, 1])
    assert pha.quality == pytest.approx([0, 0, 0, 0, 2])

    # This is corrected by the bin width and the areascal values
    assert ui.get_dep() == pytest.approx([12000, 1500, 0, 1000])

    assert pha.counts == pytest.approx([12, 2, 1, 0, 1])


@requires_fits  # only needed because of incorrect serialization
@requires_group
def test_restore_load_arrays_pha_response(check_str):
    """Can we re-create a load_arrays/DataPHA case?

    This tests whether we can restore manually-created responses.  At
    the moment we can not, so this is a regression test.

    """

    dset = ui.DataPHA("ex", [1, 2, 3, 4, 5], [12, 2, 1, 0, 1])
    dset.exposure = 100
    dset.backscal = 0.002
    dset.areascal = 0.001

    ui.set_data(dset)
    ui.group_counts(3, tabStops=np.asarray([0, 0, 0, 1, 0]))

    # Create an ARF and RMF
    #
    egrid = np.asarray([0.1, 0.2, 0.4, 0.8, 1.2, 1.6])
    elo = egrid[:-1]
    ehi = egrid[1:]
    ui.set_arf(ui.create_arf(elo, ehi))
    ui.set_rmf(ui.create_rmf(elo, ehi))

    compare(check_str, _canonical_load_arrays_pha_response)

    # The error response likely depends on the backend, but it happens
    # because the ARF "file name" - in this case "test-arf" - does not
    # exist, hence the restoration fails.
    #
    with pytest.raises(IOErr):
        restore()

    # TODO: come up with tests once the state can be restored


def test_restore_load_arrays_data2d(check_str):
    """Can we re-create a load_arrays/Data2D case

    Note that test_restore_dataspace2d_img is a basic DataIMG test.

    """

    dset = ui.Data2D("ex", [1, 1, 2], [1, 2, 2], [3, 4, 5], None, [0.1, 0.1, 0.2])
    ui.set_data(dset)

    assert ui.get_dep() == pytest.approx([3, 4, 5])
    assert ui.get_staterror() == pytest.approx([0.1, 0.1, 0.2])

    compare(check_str, _canonical_load_arrays_data2d)

    restore()

    assert isinstance(ui.get_data(), ui.Data2D)
    assert ui.get_data().shape is None
    assert ui.get_dep() == pytest.approx([3, 4, 5])
    assert ui.get_staterror() == pytest.approx([0.1, 0.1, 0.2])


@requires_xspec
def test_canonical_xspec_hard_limit_min(check_str):
    "Can we save an XSPEC model with the hard limit extended: min"

    # Reset the optimiser parameters to make them easy to check,
    # even if they are un-usable.
    #
    for key in ui.get_method_opt().keys():
        ui.set_method_opt(key, 1)

    ui.create_model_component('xspowerlaw', 'mdl')
    mdl.phoindex.set(val=-5, hard_min=-5, frozen=True)
    mdl.norm.max = 100

    compare(check_str, _canonical_xspec_hard_limit_min)


@requires_xspec
def test_canonical_xspec_hard_limit_max(check_str):
    "Can we save an XSPEC model with the hard limit extended: max"

    # Reset the optimiser parameters to make them easy to check,
    # even if they are un-usable.
    #
    for key in ui.get_method_opt().keys():
        ui.set_method_opt(key, 1)

    ui.create_model_component('xspowerlaw', 'mdl')
    mdl.phoindex.set(val=15, hard_max=15, frozen=True)
    mdl.norm.max = 100

    compare(check_str, _canonical_xspec_hard_limit_max)


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

    delta = np.sqrt((m2.c0 - 23)**2)
    m1.fwhm = 2 * np.exp(delta / 14)

    m2.c0 = 30

    expected = 2 * np.exp(0.5)

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


def test_link_par(check_str):
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

    compare(check_str, _canonical_link_par)

    restore()

    assert m2.c0.val == 300
    assert m2.c0.min == 10
    assert m2.c0.max == 500
    assert m2.c0.link is not None
    assert m2.c0.link.name == "m1.c0 + sep.c0"


@pytest.mark.xfail
@requires_data
@requires_fits
def test_pha_full_model(make_data_path, check_str):
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

    compare(check_str, add_datadir_path(_canonical_pha_full_model))

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
def test_load_data(make_data_path, check_str):
    """Check load_data path for Data1D case with staterror"""

    # Unfortunately we need to care about the backend here.
    #
    from sherpa.astro import io

    ui.set_stat("chi2datavar")

    infile = make_data_path("data1.dat")
    ui.load_data(infile, ncols=3)

    expected_x = np.arange(0.5, 11)
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
    compare(check_str, expected_output)

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
def test_load_data_basic(make_data_path, check_str):
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
    assert np.ptp(ui.get_dep()) == pytest.approx(ptp)

    with pytest.raises(StatErr,
                       match="^If you select chi2 as the statistic, all datasets must provide a staterror column$"):
        ui.get_staterror()

    with pytest.raises(DataErr,
                       match="data set '1' does not specify systematic errors"):
        ui.get_syserror()

    expected_output = add_datadir_path(_canonical_load_data_basic)
    compare(check_str, expected_output)

    restore()

    expected = '20.1000:24.9000,27.1000:89.9000'
    assert ui.get_filter() == expected

    x = ui.get_indep()
    assert len(x) == 1
    assert len(x[0]) == 1000
    assert x[0][0] == pytest.approx(0)
    assert x[0][-1] == pytest.approx(99.9)
    assert np.ptp(ui.get_dep()) == pytest.approx(ptp)

    with pytest.raises(StatErr,
                       match="^If you select chi2 as the statistic, all datasets must provide a staterror column$"):
        ui.get_staterror()

    with pytest.raises(DataErr,
                       match="data set '1' does not specify systematic errors"):
        ui.get_syserror()


@requires_data
@requires_fits
def test_restore_pha_multiple_backgrounds(make_data_path, check_str):
    """Can we restore a grating dataset with multiple backgrounds?"""

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
    compare(check_str, expected_output)

    restore()

    check_data()


@requires_data
@requires_fits
def test_restore_pha2(make_data_path, check_str):
    """Can we restore a pha2 file?"""

    # Note: not including .gz for the file name
    ui.load_pha(make_data_path("3c120_pha2"))

    # Before #1564 we could say
    #
    #    ui.set_analysis("wave")
    #    ui.notice(2, 4)
    #
    for idval in [1, 6]:
        ui.notice_id(idval, 3793, 6192)

    for idval in [2, 5]:
        ui.notice_id(idval, 5393, 6992)

    for idval in [3, 4]:
        ui.notice_id(idval, 6993, 7792)

    for idval in [7, 12]:
        ui.notice_id(idval, 5993, 7192)

    for idval in [8, 11]:
        ui.notice_id(idval, 6793, 7592)

    for idval in [9, 10]:
        ui.notice_id(idval, 7593, 7992)

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
    compare(check_str, expected_output)

    restore()

    check_data()


@requires_data
@requires_fits
def test_restore_pha2_delete(make_data_path, check_str):
    """Can we restore a pha2 file after deleting files?

    This is a regression test so we can see as soon as things have
    changed, rather than marking it as xfail.

    """

    # Note: not including .gz for the file name
    ui.load_pha(make_data_path("3c120_pha2"))

    # Select just the first order (|TG_M| = 1) datasets.
    #
    for idval in [1, 2, 5, 6, 7, 8, 11, 12]:
        ui.delete_data(idval)

    def check_data():
        assert ui.list_data_ids() == pytest.approx([10, 3, 4, 9])

    check_data()

    expected_output = add_datadir_path(_canonical_pha2_delete)
    compare(check_str, expected_output)

    restore()

    check_data()


@requires_data
@requires_fits
def test_restore_pha_hetg(make_data_path, check_str):
    """Can we restore a HETG file?

    See issue #320.

    This is a regression test so we can see as soon as things have
    changed, rather than marking it as xfail.

    """

    # Note: not including .gz for the file name
    ui.load_pha(make_data_path("3c120_heg_-1.pha"))
    ui.set_analysis("wave")
    ui.notice(5, 10)

    # Make sure we have a copy of the count values.
    #
    c = ui.get_data().counts.copy()
    b1 = ui.get_bkg(bkg_id=1).counts.copy()
    b2 = ui.get_bkg(bkg_id=2).counts.copy()

    # Check they are different
    assert (c != b1).any()
    assert (c != b2).any()
    assert (b1 != b2).any()

    def check_data():
        assert ui.list_data_ids() == pytest.approx([1])
        data = ui.get_data(1)
        assert data.background_ids == pytest.approx([1, 2])
        assert data.response_ids == pytest.approx([1])

        assert ui.get_arf().name.endswith("/3c120_heg_-1.arf")
        assert ui.get_rmf().name.endswith("/3c120_heg_-1.rmf")

        # Although the backgrounds do not have a response, they get
        # automatically set to the source response.
        #
        for bkg_id in [1, 2]:
            bkg = ui.get_bkg(bkg_id=bkg_id)
            assert bkg.response_ids == [1], f"bkg_id={bkg_id}"

            assert ui.get_arf(bkg_id=bkg_id).name.endswith("/3c120_heg_-1.arf")
            assert ui.get_rmf(bkg_id=bkg_id).name.endswith("/3c120_heg_-1.rmf")

        assert ui.get_data().counts == pytest.approx(c)
        assert ui.get_bkg(bkg_id=1).counts == pytest.approx(b1)
        assert ui.get_bkg(bkg_id=2).counts == pytest.approx(b2)

    check_data()

    expected_output = add_datadir_path(_canonical_pha_hetg)
    compare(check_str, expected_output)

    restore()

    check_data()


@requires_data
@requires_group
@requires_fits
@requires_xspec
def test_restore_pha_csc(make_data_path, check_str):
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
    compare(check_str, expected_output)

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
    x1, x0 = np.mgrid[1:4, 1:3]
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.arange(6) * 10 + 10

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
    x1, x0 = np.mgrid[1:3, 1:4]
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.arange(6) * 10 + 10

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


@requires_data
@requires_fits
@requires_group
def test_copy_data_pha(make_data_path):
    """Can we copy a PHA dataset and track it?"""

    # This file does not have ANCR/BACK/RESPFILE settings.
    #
    ui.load_pha(make_data_path("source.pi"), use_errors=True)
    ui.load_rmf(make_data_path("source.rmf"))
    ui.load_arf(make_data_path("source.arf"))

    # Is the filter and grouping copied over?
    ui.notice(0.5, 7)
    ui.group_counts(20)

    # What can we easily check to see that we have the filtered and
    # grouped dataset?
    #
    filter_expr = ui.get_filter()
    dep = ui.get_dep(filter=True)
    evals = ui.get_staterror()

    ui.copy_data(1, 2)
    ui.delete_data(1)

    def check_data():
        assert ui.list_data_ids() == [2]
        assert ui.list_bkg_ids(2) == []
        assert ui.get_data(2).name.endswith("/source.pi")
        assert ui.get_arf(2).name.endswith("/source.arf")
        assert ui.get_rmf(2).name.endswith("/source.rmf")

        assert ui.get_filter(2) == filter_expr
        assert ui.get_dep(2, filter=True) == pytest.approx(dep)
        assert ui.get_staterror(2) == pytest.approx(evals)

    check_data()

    restore()

    check_data()
