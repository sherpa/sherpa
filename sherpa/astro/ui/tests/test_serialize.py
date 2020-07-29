#
#  Copyright (C) 2015, 2016, 2018, 2019, 2020  Smithsonian Astrophysical Observatory
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

# TODO add tests:
#    set_full_model
#    multiple sources
#    linked parameters
#    iter fit method
#    check the noticed range after restoring it
#    pha dataset; wavelength analysis
#    pha2 dataset
#    Chandra Source Catalog "pha3" dataset (src and bgnd in different
#      extensions); DougBurke has noted that astropy backend fails if
#      response files are gzipp-ed [may be acceptable behavior] but
#      want test to check basic behavior
#    psf model
#    table model
#

from io import StringIO
import logging
import re
import tempfile

import numpy
from numpy.testing import assert_array_equal

import pytest

from sherpa.utils.testing import SherpaTestCase, requires_data, \
    requires_xspec, has_package_from_list, requires_fits, requires_group
from sherpa.astro import ui
# from sherpa.astro.ui import serialize

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

# Note that the XSPEC defaults changed in XSPEC 12.10.1 (xsxsect is
# now vern rather than bcmc, but it can depend on what files a
# user has - e.g. ~/.xspec/Xspec.init will over-ride it), so there
# are explicit statements to set the values before the tests so
# that we have a "known" state.
#

# A representation of the default Sherpa state
_canonical_empty = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets



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

set_method_opt("finalsimplex", 9)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("initsimplex", 0)
set_method_opt("iquad", 1)
set_method_opt("maxfev", 5000)
set_method_opt("step", None)
set_method_opt("verbose", 1)


######### Set Model Components and Parameters



######### Set Source, Pileup and Background Models

"""

# use @@ as a marker to indicate that the location of the data directory
# needs to be inserted into the text before testing.
#
_canonical_pha_basic = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha(1, "@@/3c273.pi")

######### Set Image Coordinates

if get_data(1).grouping is not None and not get_data(1).grouped:
    ######### Group Data
    group(1)

######### Data Spectral Responses

load_arf(1, "@@/3c273.arf", resp_id=1)
load_rmf(1, "@@/3c273.rmf", resp_id=1)

######### Load Background Data Sets

load_bkg(1, "@@/3c273_bg.pi", bkg_id=1)
if get_bkg(1, 1).grouping is not None and not get_bkg(1, 1).grouped:
    ######### Group Background
    group(1, 1)

######### Background Spectral Responses

load_arf(1, "@@/3c273.arf", resp_id=1, bkg_id=1)
load_rmf(1, "@@/3c273.rmf", resp_id=1, bkg_id=1)

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

src.redshift.default_val = 0.0
src.redshift.default_min = -0.999
src.redshift.default_max = 10.0
src.redshift.val     = 0.0
src.redshift.min     = -0.999
src.redshift.max     = 10.0
src.redshift.units   = ""
src.redshift.frozen  = True

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

load_pha("grp", "@@/3c273.pi")

######### Set Image Coordinates


######### Data grouping flags

set_grouping("grp", val=numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], numpy.int16))

######### Data quality flags

set_quality("grp", val=numpy.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], numpy.int16))
if get_data("grp").grouping is not None and not get_data("grp").grouped:
    ######### Group Data
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
if get_bkg("grp", 1).grouping is not None and not get_bkg("grp", 1).grouped:
    ######### Group Background
    group("grp", 1)

######### Background Spectral Responses

load_arf("grp", "@@/3c273.arf", resp_id=1, bkg_id=1)
load_rmf("grp", "@@/3c273.rmf", resp_id=1, bkg_id=1)

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
ggal.nH.default_max = 100000.0
ggal.nH.val     = 2.0
ggal.nH.min     = 0.0
ggal.nH.max     = 100000.0
ggal.nH.units   = "10^22 atoms / cm^2"
ggal.nH.frozen  = True



######### Set Source, Pileup and Background Models

set_source("grp", (xsphabs.ggal * powlaw1d.gpl))



"""

_canonical_pha_back = """import numpy
from sherpa.astro.ui import *

######### Load Data Sets

load_pha("bgrp", "@@/3c273.pi")

######### Set Image Coordinates

if get_data("bgrp").grouping is not None and not get_data("bgrp").grouped:
    ######### Group Data
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
if get_bkg("bgrp", 1).grouping is not None and not get_bkg("bgrp", 1).grouped:
    ######### Group Background
    group("bgrp", 1)

######### Background Spectral Responses

load_arf("bgrp", "@@/3c273.arf", resp_id=1, bkg_id=1)
load_rmf("bgrp", "@@/3c273.rmf", resp_id=1, bkg_id=1)

######### Set Energy or Wave Units

set_analysis("bgrp", 'energy', "rate", 0)

######### Filter Data

notice_id("bgrp", "0.518300011754:6.234200000763")
notice_id("bgrp", None, None, bkg_id=1)
notice_id("bgrp", "1.890699982643:7.723400115967", bkg_id=1)


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
ggal.nH.default_max = 100000.0
ggal.nH.val     = 2.0
ggal.nH.min     = 0.0
ggal.nH.max     = 100000.0
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

set_source("bgrp", (xsphabs.ggal * powlaw1d.gpl))

set_bkg_source("bgrp", (steplo1d.bstep + polynom1d.bpoly), bkg_id=1)


######### XSPEC Module Settings

set_xschatter(0)
set_xsabund("lodd")
set_xscosmo(72, 0.02, 0.71)
set_xsxsect("vern")
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

set_method_opt("finalsimplex", 9)
set_method_opt("ftol", 1.19209289551e-07)
set_method_opt("initsimplex", 0)
set_method_opt("iquad", 1)
set_method_opt("maxfev", None)
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

set_source(3, (sin.sin_model + usermodel.mymodel))

"""

if has_xspec:
    _canonical_extra = """
######### XSPEC Module Settings

set_xschatter(0)
set_xsabund("angr")
set_xscosmo(70, 0, 0.73)
set_xsxsect("bcmc")
"""

else:
    _canonical_extra = ""

_canonical_empty += _canonical_extra
_canonical_empty_stats += _canonical_extra
_canonical_pha_basic += _canonical_extra
_canonical_pha_grouped += _canonical_extra
_canonical_usermodel += _canonical_extra


# This extends the clean_astro_ui logic but it isn't obvious
# to me how to chain them together, so it re-implements
# clean_astro_ui.
#
@pytest.fixture(autouse=True)
def setup(clean_astro_ui):

    # Setup for the test
    try:
        numpy.set_printoptions(legacy='1.13')
    except TypeError:  # numpy < 1.14
        pass

    old_logging_level = logger.level
    logger.setLevel(logging.ERROR)

    if has_xspec:
        from sherpa.astro import xspec
        old_xspec = xspec.get_xsstate()

        # ensure we have the same settings as the test cases
        # used (changes to XSPEC may have changed the defaults)
        xspec.set_xschatter(0)
        xspec.set_xsabund('angr')
        xspec.set_xsxsect('bcmc')
    else:
        old_xspec = None

    ui.clean()

    # run the test
    yield

    # Restore the Sherpa/related settings
    ui.clean()

    if old_xspec is not None:
        xspec.set_xsstate(old_xspec)

    try:
        numpy.set_printoptions(legacy=False)
    except TypeError:  # numpy < 1.14
        pass

    logger.setLevel(old_logging_level)


def add_datadir_path(output):
    """Replace any @@ characters by the value of self.datadir,
    making sure that the replacement text does not end in a /."""

    # it would be nice not to use SherpaTestCase
    dname = SherpaTestCase.datadir
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
    try:
        exec(output)
        success = True
        e = "no exception"
    except Exception as e:
        success = False

    assert success, "exception={}".format(e)


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


def test_canonical_empty():
    "Contents of empty state are as expected"
    compare(_canonical_empty)


def test_canonical_empty_outfile():
    "Can read in a save file"
    tfile = tempfile.NamedTemporaryFile(suffix='.sherpa')
    ui.save_all(tfile.name, clobber=True)
    with open(tfile.name, 'r') as fh:
        output = fh.read()

    compare_lines(_canonical_empty, output)


def test_canonical_empty_stats():
    "Change several settings but load no data"

    ui.set_stat('leastsq')

    ui.set_method('simplex')
    ui.set_method_opt('maxfev', 5000)
    ui.set_method_opt('verbose', 1)

    compare(_canonical_empty_stats)


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
    assert src_expr.name == '(xsphabs.gal * (powlaw1d.pl + xsapec.src))'
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
    assert src_expr.name == '(xsphabs.ggal * powlaw1d.gpl)'
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

    # TODO: check noticed range

    src_expr = ui.get_source('bgrp')
    assert src_expr.name == '(xsphabs.ggal * powlaw1d.gpl)'

    bg_expr = ui.get_bkg_source('bgrp')
    assert bg_expr.name == '(steplo1d.bstep + polynom1d.bpoly)'

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
    assert src_expr.name == '(sin.sin_model + usermodel.mymodel)'
    mymodel = ui.get_model_component("mymodel")
    assert mymodel.m.frozen,"is mymodel.m frozen?"
    assert mymodel.c.val == 2.0
    assert mymodel.c.units == "m"
    assert mymodel.m.max == 5.5
    assert mymodel.m.units == ""

    assert ui.calc_stat(3) == pytest.approx(statval)
