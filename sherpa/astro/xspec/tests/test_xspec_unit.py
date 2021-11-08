#
#  Copyright (C) 2016-2018, 2019, 2020, 2021
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

# Supplement test_xspec.py with py.test tests
#

import copy
import logging
import os
from tempfile import NamedTemporaryFile, gettempdir

import pytest

import numpy as np
from numpy.testing import assert_almost_equal

from sherpa.models.basic import Const1D
from sherpa.utils.testing import requires_data, requires_fits, requires_xspec
from sherpa.utils.err import ParameterErr


# It is hard to test many of the state routines, since it requires
# a full understanding of how they are implemented; the simplest
# way is to check those that change the state and ensure values
# pass through a round trip - i.e. that if you set it to a given
# value then request it, you get back the value you set.
#
# There is currently no (or limited) checks for invalid inputs.
#

# The following can depend on the XSPEC version; how much of a
# contract do we want to make with the initialization code to set
# to our values versus the default XSPEC settings?
#
# For XSPEC 12.10.1 and later, the default settings depend on
#  - the user's ~/.xspec/Xspec.init file
#  - the $HEADAS/../spectral/manaer/Xspec.init file
#  - in-built dfeaults
#
# This means that it is now hard to reliably check the default
# values. So, we now just chcek that the default values are
# one of the expected values.
#
# Note that a used could set up their own abundance table,
# in which case it is not obvious what to do.
#
DEFAULT_ABUND = ['angr', 'aspl', 'feld', 'aneb', 'grsa', 'wilm', 'lodd']
DEFAULT_XSECT = ['bcmc', 'obcm', 'vern']

# XSPEC presmably has its own default for these, but Sherpa explicitly
# sets values for these parameters.
#
DEFAULT_CHATTER = 0
DEFAULT_COSMO = (70.0, 0.0, 0.73)

# The XSET names are (currently) not validated, so pick a name that
# will not be used by an actual model so we can set it.
#
DEFAULT_XSET_NAME = 'SHERPA-TEST-DUMMY-NAME'
DEFAULT_XSET_VALUE = ''

# The number of elements in the abundance table
ELEMENT_NAMES = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
                 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
                 'Cu', 'Zn']
NELEM = len(ELEMENT_NAMES)


@pytest.mark.xfail
@requires_xspec
def test_chatter_default():
    """Check the expected default setting for chatter.

    Ideally this test would be run before any other
    tests of XSPEC are made (i.e. any XSPEC code is called).

    For some reason the chatter setting appears to default to 50,
    rather than 0. It is not at all clear why, so has been marked
    as xfail. DJB now thinks that this is due to the Sherpa test
    setup, which is designed to make the models be verbose in
    case there's an error.

    """

    from sherpa.astro import xspec

    oval = xspec.get_xschatter()
    assert oval == DEFAULT_CHATTER


@requires_xspec
def test_version():
    """Can we get at the XSPEC version?

    There is limited testing of the return value.
    """

    from sherpa.astro import xspec

    v = xspec.get_xsversion()
    assert isinstance(v, (str, ))
    assert len(v) > 0

    # Could check it's of the form a.b.c[optional] but leave that for
    # now.


@requires_xspec
def test_abund_default():
    """Check the expected default setting for the abundance.

    Ideally this test would be run before any other
    tests of XSPEC are made (i.e. any XSPEC code is called).
    """

    from sherpa.astro import xspec

    oval = xspec.get_xsabund()
    assert oval in DEFAULT_ABUND


@requires_xspec
def test_xset_default():
    """Check the expected default setting for the xset setting.

    Ideally this test would be run before any other
    tests of XSPEC are made (i.e. any XSPEC code is called).

    This is tricky since XSPEC does not return a value until
    one has been set (that is, if it returns '' then the
    setting is presumably taken to be "use the default value"
    by the model code).
    """

    from sherpa.astro import xspec

    # The test is case insensitive, but this test doesn't really
    # test this out (since it is expected to return '' whatever
    # the input name is).
    #
    name = DEFAULT_XSET_NAME.lower()
    oval = xspec.get_xsxset(name)
    assert oval == DEFAULT_XSET_VALUE


@requires_xspec
def test_xsect_default():
    """Check the expected default setting for the xsect setting.

    Ideally this test would be run before any other
    tests of XSPEC are made (i.e. any XSPEC code is called).
    """

    from sherpa.astro import xspec

    oval = xspec.get_xsxsect()
    assert oval in DEFAULT_XSECT


@requires_xspec
def test_manager_path_default():
    """Check the expected default setting for the manager path.

    Ideally this test would be run before any other
    tests of XSPEC are made (i.e. any XSPEC code is called).
    """

    # Is this always going to be correct?
    default_path = os.path.join(os.environ['HEADAS'],
                                '../spectral/manager')

    from sherpa.astro import xspec

    oval = xspec.get_xspath_manager()

    # Normalize the paths to remove double / - e.g. if HEADAS is
    # set to /a/b/c/ rather than /a/b/c.
    #
    assert os.path.normpath(oval) == os.path.normpath(default_path)


@requires_xspec
def test_model_path_default():
    """Check the expected default setting for the model data path.

    Ideally this test would be run before any other
    tests of XSPEC are made (i.e. any XSPEC code is called).
    """

    from sherpa.astro import xspec

    # Is this always going to be correct?
    #
    try:
        default_path = os.environ['XSPEC_MDATA_DIR']
    except KeyError:
        default_path = os.path.join(os.environ['HEADAS'],
                                    '../spectral/modelData/')

    oval = xspec.get_xspath_model()

    # Normalize the paths to remove double / - e.g. if HEADAS is
    # set to /a/b/c/ rather than /a/b/c.
    #
    assert os.path.normpath(oval) == os.path.normpath(default_path)


@requires_xspec
def test_cosmo_default():
    """Check the expected default setting for the cosmology settings.

    Ideally this test would be run before any other
    tests of XSPEC are made (i.e. any XSPEC code is called).
    """

    from sherpa.astro import xspec

    oval = xspec.get_xscosmo()

    # Since this is a tuple of numbers, check individually
    assert len(oval) == 3
    assert oval[0] == pytest.approx(DEFAULT_COSMO[0])
    assert oval[1] == pytest.approx(DEFAULT_COSMO[1])
    assert oval[2] == pytest.approx(DEFAULT_COSMO[2])


@requires_xspec
def test_abund_element():
    """Can we access the elemental settings?
    """

    from sherpa.astro import xspec

    oval = xspec.get_xsabund()
    try:
        xspec.set_xsabund('wilm')
        h = xspec.get_xsabund('H')
        he = xspec.get_xsabund('He')
        si = xspec.get_xsabund('Si')
        ar = xspec.get_xsabund('Ar')
        k = xspec.get_xsabund('K')
        fe = xspec.get_xsabund('Fe')

    finally:
        xspec.set_xsabund(oval)

    # These values were found from HEASOFT version 6.19
    # spectral/manager/abundances.dat
    # The values are given to two decimal places in this file.
    # It is not worth testing all settings, since we are not
    # testing the XSPEC implementation itself, just our use of it.
    #
    assert h == pytest.approx(1.0)
    assert he == pytest.approx(9.77e-2)
    assert si == pytest.approx(1.86e-05)
    assert ar == pytest.approx(2.57e-06)
    assert k == pytest.approx(0.0)
    assert fe == pytest.approx(2.69e-05)


def validate_xspec_setting(getfunc, setfunc, newval, altval):
    """Check we can change an XSPEC setting.

    Parameters
    ----------
    getfunc : function
        The XSPEC function to query the setting: it returns a
        value and has no arguments.
    setfunc : function
        The XSPEC function to change the setting: it takes a
        single argument and returns nothing.
    newval, altval
        The value to use (newval) and an alternative (altval) if
        the current setting is already at newval (this is perhaps
        a bit excessive but it avoids issues if other tests have
        changed things).
    """

    oval = getfunc()
    if oval == newval:
        nval = altval
    else:
        nval = newval

    try:
        setfunc(nval)
        xval = getfunc()
    finally:
        setfunc(oval)

    assert xval == nval

    # As a sanity check ensure we are back at the starting point
    assert getfunc() == oval


def validate_xspec_state_setting(key, newval, altval):
    """Check we can change an XSPEC setting via the state mechanism

    Parameters
    ----------
    key : string
        The name of the setting (e.g. 'abund' or 'xsect').
    newval, altval
        The value to use (newval) and an alternative (altval) if
        the current setting is already at newval (this is perhaps
        a bit excessive but it avoids issues if other tests have
        changed things).
    """

    from sherpa.astro import xspec

    ostate = xspec.get_xsstate()

    def getfunc():
        return xspec.get_xsstate()[key]

    def setfunc(val):
        nstate = ostate.copy()
        nstate[key] = val
        xspec.set_xsstate(nstate)

    validate_xspec_setting(getfunc, setfunc, newval, altval)

    assert xspec.get_xsstate() == ostate


@requires_xspec
def test_chatter_change():
    """Can we change the chatter setting."""

    from sherpa.astro import xspec

    validate_xspec_setting(xspec.get_xschatter,
                           xspec.set_xschatter,
                           25, 35)


@requires_xspec
def test_abund_change_string():
    """Can we change the abundance setting: string

    This only checks that we can use one of the hard-coded
    abundance names. It does not check the file I/O.
    """

    from sherpa.astro import xspec

    validate_xspec_setting(xspec.get_xsabund,
                           xspec.set_xsabund,
                           'grsa', 'wilm')


@requires_xspec
def test_abund_change_file():
    """Can we change the abundance setting: file

    This test hard-codes the number of elements expected in the
    file.
    """

    from sherpa.astro import xspec

    elems = {n: i * 0.1 for i, n in enumerate(ELEMENT_NAMES)}

    tfh = NamedTemporaryFile(mode='w', suffix='.xspec')
    for n in ELEMENT_NAMES:
        tfh.write("{}\n".format(elems[n]))

    tfh.flush()

    oval = xspec.get_xsabund()
    try:
        xspec.set_xsabund(tfh.name)

        abund = xspec.get_xsabund()
        out = {n: xspec.get_xsabund(n)
               for n in ELEMENT_NAMES}

    finally:
        xspec.set_xsabund(oval)

    assert abund == 'file'
    for n in ELEMENT_NAMES:
        assert out[n] == pytest.approx(elems[n])


@requires_xspec
def test_xset_change():
    """Can we change the xset setting.
    """

    from sherpa.astro import xspec

    def getfunc():
        return xspec.get_xsxset(DEFAULT_XSET_NAME)

    def setfunc(val):
        xspec.set_xsxset(DEFAULT_XSET_NAME.lower(), val)

    val1 = 'dummy value'
    val2 = 'a different setting'
    validate_xspec_setting(getfunc, setfunc, val1, val2)

    # A separate part of the XSET interface is that the settings
    # are recorded in the XSPEC state maintained by the xspec
    # module, so check that the stored value is included in this.
    #
    modelvals = xspec.get_xsstate()['modelstrings']
    assert DEFAULT_XSET_NAME in modelvals

    # Is it worth changing the code so we know which to check for?
    assert modelvals[DEFAULT_XSET_NAME] in [val1, val2]


@requires_xspec
def test_xsect_change():
    """Can we change the xsect setting: string

    This only checks that we can use one of the hard-coded
    abundance names. It does not check the file I/O.
    """

    from sherpa.astro import xspec

    validate_xspec_setting(xspec.get_xsxsect,
                           xspec.set_xsxsect,
                           'obcm', 'vern')


@requires_xspec
def test_cosmo_change():
    """Can we change the cosmology settings.
    """

    from sherpa.astro import xspec

    old_h0, old_q0, old_l0 = xspec.get_xscosmo()

    new_h0 = 51.0
    new_q0 = 0.2
    new_l0 = 0.48

    if old_h0 == pytest.approx(new_h0):
        new_h0 -= 10.0
    if old_q0 == pytest.approx(new_q0):
        new_q0 -= 0.05
    if old_l0 == pytest.approx(new_l0):
        new_l0 += 0.01

    try:
        xspec.set_xscosmo(new_h0, new_q0, new_l0)
        nval_h0, nval_q0, nval_l0 = xspec.get_xscosmo()
    finally:
        xspec.set_xscosmo(old_h0, old_q0, old_l0)

    assert nval_h0 == pytest.approx(new_h0)
    assert nval_q0 == pytest.approx(new_q0)
    assert nval_l0 == pytest.approx(new_l0)


@requires_xspec
def test_path_manager_change():
    """Can we change the manager-path setting?
    """

    from sherpa.astro import xspec

    validate_xspec_setting(xspec.get_xspath_manager,
                           xspec.set_xspath_manager,
                           '/dev/null',
                           gettempdir())


# Note that the XSPEC state is used in test_xspec.py, but only
# to save/restore the state after each test. There is no
# explicit test there of the functionality. The state tests here
# are very basic.
#

@requires_xspec
def test_get_xsstate_keys():
    """Check get_xsstate returns the expected keys.

    Checking the values here are hard, unless we save/restore
    the state in the requires_xspec decorator or essentially
    replicate the implementation of get_xsstate.
    """

    from sherpa.astro import xspec

    ostate = xspec.get_xsstate()
    assert isinstance(ostate, dict)

    for key in ["abund", "chatter", "cosmo", "xsect",
                "modelstrings", "paths"]:
        assert key in ostate


@requires_xspec
def test_set_xsstate_missing_key():
    """Check set_xsstate does nothing if required key is missing.

    """

    from sherpa.astro import xspec

    ostate = xspec.get_xsstate()

    for val in ostate.values():
        assert val is not None

    # paths is not a required key
    #
    req_keys = ["abund", "chatter", "cosmo", "xsect",
                "modelstrings"]

    fake = {'abund': ostate['abund'] + '_copy',
            'xsect': ostate['xsect'] + '_copy',
            'chatter': -10,
            'cosmo': (0.0, 0.0),  # two elements will cause a failure
            'modelstrings': {'foo': 2, 'bar': None},
            'paths': {'manager': '/dev/null'}}

    for key in req_keys:

        copy = fake.copy()
        del copy[key]
        xspec.set_xsstate(copy)

        nstate = xspec.get_xsstate()
        assert nstate == ostate


@requires_xspec
def test_set_xsstate_abund():
    """Check set_xsstate works for abundance.
    """

    validate_xspec_state_setting('abund', 'lodd', 'wilm')


@requires_xspec
def test_set_xsstate_xsect():
    """Check set_xsstate works for cross sections.
    """

    validate_xspec_state_setting('xsect', 'vern', 'obcm')


@requires_xspec
def test_set_xsstate_chatter():
    """Check set_xsstate works for chatter.
    """

    validate_xspec_state_setting('chatter', 5, 15)


@requires_xspec
def test_set_xsstate_xset():
    """Check set_xsstate works for an xset command.
    """

    from sherpa.astro import xspec

    ostate = xspec.get_xsstate()

    key = 'a-test-keyword'
    val = '/foo/bar/baz.pha'
    while key in ostate['modelstrings']:
        key += "a"

    ukey = key.upper()

    # There should be no value for this key (since it isn't
    # in modelstrings by construction).
    #
    assert key not in xspec.modelstrings
    assert xspec.get_xsxset(key) == ''

    nstate = copy.deepcopy(ostate)
    nstate['modelstrings'][key] = val
    xspec.set_xsstate(nstate)

    assert xspec.get_xsxset(key) == val
    assert ukey in xspec.modelstrings
    assert xspec.modelstrings[ukey] == val

    xspec.set_xsstate(ostate)

    # Unfortunately, due to there being no attempt at clearing out the
    # XSET settings (e.g. removing existing settings before restoring
    # the state), the following tests fail.
    #
    # TODO: the code should probably be updated to fix this
    #
    # assert xspec.get_xsxset(key) == ''
    # assert xspec.get_xsstate() == ostate

    xspec.set_xsxset(key, '')
    del xspec.modelstrings[ukey]
    assert xspec.get_xsstate() == ostate


@requires_xspec
def test_set_xsstate_path_manager():
    """Check set_xsstate works for the manager path
    """

    from sherpa.astro import xspec

    ostate = xspec.get_xsstate()
    opath = xspec.get_xspath_manager()

    spath = ostate['paths'].get('manager', None)

    # This is just an internal validation check
    if spath is not None:
        assert spath == opath

    if opath == 'b/a':
        npath = 'a/b'
    else:
        npath = 'b/a'

    nstate = copy.deepcopy(ostate)
    nstate['paths']['manager'] = npath
    xspec.set_xsstate(nstate)

    assert xspec.get_xspath_manager() == npath

    xspec.set_xsstate(ostate)

    # Similar to the state xset tests, using an empty
    # dictionary for paths does not clear out/reset the
    # manager path. In this case it's not obvious what
    # should be done (as there's no obvious default value
    # to use, unless we fall back to the FNINIT-created
    # value, which is not ideal since there's no guarantee
    # that we will notice any changes to that logic).
    #
    # This is an edge case.
    #
    # assert xspec.get_xspath_manager() == opath
    # assert xspec.get_xsstate() == ostate

    xspec.set_xspath_manager(opath)
    # should really clear out xspec.xspecpaths


@requires_data
@requires_fits
@requires_xspec
def test_read_xstable_model(make_data_path):
    """Limited test (only one file).

    Evaluation tests using this model are in
    sherpa.astro.xspec.tests.test_xspec.
    """

    from sherpa.astro import xspec

    path = make_data_path('xspec-tablemodel-RCS.mod')
    tbl = xspec.read_xstable_model('bar', path)

    assert tbl.name == 'bar'
    assert isinstance(tbl, xspec.XSTableModel)
    assert tbl.addmodel
    assert tbl.integrate

    assert len(tbl.pars) == 4
    assert tbl.pars[0].name == 'tau'
    assert tbl.pars[1].name == 'beta'
    assert tbl.pars[2].name == 't'
    assert tbl.pars[3].name == 'norm'

    assert_almost_equal(tbl.tau.val, 1)
    assert_almost_equal(tbl.tau.min, 1)
    assert_almost_equal(tbl.tau.max, 10)

    assert_almost_equal(tbl.beta.val, 0.1)
    assert_almost_equal(tbl.beta.min, 0.1)
    assert_almost_equal(tbl.beta.max, 0.5)

    assert_almost_equal(tbl.t.val, 0.1)
    assert_almost_equal(tbl.t.min, 0.1)
    assert_almost_equal(tbl.t.max, 1.3)

    assert_almost_equal(tbl.norm.val, 1)
    assert_almost_equal(tbl.norm.min, 0)
    assert_almost_equal(tbl.norm.max, 1e24)

    for p in tbl.pars:
        assert not(p.frozen)


@requires_xspec
@pytest.mark.parametrize("clsname", ["powerlaw", "wabs"])
def test_xspec_model_requires_bins(clsname):
    """Ensure you can not call with a single grid for the energies.

    You used to be able to do this (in Sherpa 4.13 and earlier).
    """

    from sherpa.astro import xspec

    mdl = getattr(xspec, 'XS{}'.format(clsname))()

    emsg = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg):
        mdl([0.1, 0.2, 0.3, 0.4])


@requires_xspec
@pytest.mark.parametrize("clsname", ["powerlaw", "wabs"])
def test_xspec_model_requires_bins_low_level(clsname):
    """Ensure you can not call with a single grid for the energies (calc).

    You used to be able to do this (in Sherpa 4.13 and earlier).
    """

    from sherpa.astro import xspec

    mdl = getattr(xspec, 'XS{}'.format(clsname))()

    emsg = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg):
        mdl.calc([p.val for p in mdl.pars], [0.1, 0.2, 0.3, 0.4])


@requires_xspec
@pytest.mark.parametrize("clsname", ["powerlaw", "wabs"])
def test_xspec_model_requires_bins_very_low_level(clsname):
    """Check we can use a single grid for direct access (_calc).

    This is the flip side to test-xspec_model_requires_bins_low_level
    """

    from sherpa.astro import xspec

    mdl = getattr(xspec, 'XS{}'.format(clsname))()

    # pick a range which does not evaluate to 0 for wabs
    egrid = np.arange(0.5, 1.0, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]

    y1 = mdl(elo, ehi)
    assert y1.size == elo.size

    y2 = mdl._calc([p.val for p in mdl.pars], egrid)
    assert y2.size == elo.size + 1

    # Scale by the median value so we have values ~ 1 for comparison
    ymed = np.median(y1)
    y1 /= ymed
    y2 /= ymed

    # This should be an exact match, but it's easier done with pytest.approx
    assert y2[:-1] == pytest.approx(y1)

    # Last element of _calc is 0
    assert y2[-1] == 0.0


@requires_fits
@requires_data
@requires_xspec
def test_xspec_tablemodel_requires_bin_edges(make_data_path, clean_astro_ui):
    """Check we can not call a table model with a single grid.

    This used to be supported in Sherpa 4.13 and before.
    """

    from sherpa.astro import xspec

    path = make_data_path('xspec-tablemodel-RCS.mod')
    tbl = xspec.read_xstable_model('bar', path)

    emsg = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg):
        tbl([0.1, 0.2, 0.3, 0.4])


@requires_fits
@requires_data
@requires_xspec
def test_xspec_tablemodel_requires_bin_edges_low_level(make_data_path, clean_astro_ui):
    """Check we can not call a table model with a single grid (calc).

    This used to be supported in Sherpa 4.13 and before.
    """

    from sherpa.astro import xspec

    path = make_data_path('xspec-tablemodel-RCS.mod')
    tbl = xspec.read_xstable_model('bar', path)

    emsg = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg):
        tbl.calc([p.val for p in tbl.pars], [0.1, 0.2, 0.3, 0.4])


@requires_xspec
def test_xspec_convolutionmodel_requires_bin_edges():
    """Check we can not call a convolution model with a single grid.

    This used to be supported in Sherpa 4.13 and before.
    """

    import sherpa.astro.xspec as xs

    m1 = xs.XSpowerlaw()
    m2 = xs.XScflux()
    mdl = m2(m1)

    emsg = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg):
        mdl([0.1, 0.2, 0.3, 0.4])


@requires_xspec
def test_xspec_convolutionmodel_requires_bin_edges_low_level():
    """Check we can not call a convolution model with a single grid (calc).

    This used to be supported in Sherpa 4.13 and before.
    """

    import sherpa.astro.xspec as xs

    m1 = xs.XSpowerlaw()
    m2 = xs.XScflux()
    mdl = m2(m1)

    emsg = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg):
        mdl.calc([p.val for p in mdl.pars], [0.1, 0.2, 0.3, 0.4])


@requires_data
@requires_fits
@requires_xspec
def test_evaluate_xspec_additive_model_beyond_grid(make_data_path):
    """Can we extend an additive table model beyond its grid?

    XSPEC 12.10.0 (if not manually patched) will crash if an
    XSPEC table model is evaluated beyond its grid (this is only
    an issue for programs like Sherpa that use XSPEC as a "library").
    """

    from sherpa.astro import xspec

    if xspec.get_xsversion().startswith('12.10.0'):
        pytest.skip('Test known to crash XSPEC 12.10.0')

    path = make_data_path('xspec-tablemodel-RCS.mod')
    tbl = xspec.read_xstable_model('bar', path)

    egrid = np.arange(0.1, 11, 0.01)
    elo = egrid[:-1]
    ehi = egrid[1:]
    y = tbl(elo, ehi)

    # Several simple regression tests.
    assert y[0] == pytest.approx(0.27216572)
    assert y.max() == pytest.approx(0.3047457)
    assert y.min() == 0.0

    # Is the following worth it?
    minval = 1.2102469e-11
    assert y[y > 0].min() == pytest.approx(minval)
    assert y[967] == pytest.approx(minval)

    zeros = np.where(y <= 0)
    assert (zeros[0] == np.arange(968, 1089)).all()


@requires_data
@requires_fits
@requires_xspec
def test_create_xspec_multiplicative_model(make_data_path):
    """Can we load multiplicative table models?
    """

    from sherpa.astro import xspec

    path = make_data_path('testpcfabs.mod')
    tbl = xspec.read_xstable_model('bar', path)

    assert tbl.name == 'bar'
    assert isinstance(tbl, xspec.XSTableModel)
    assert not tbl.addmodel
    assert tbl.integrate

    # Apparently we lose the case of the parameter names;
    # should investigate
    #
    assert len(tbl.pars) == 2
    assert tbl.pars[0].name == 'nh'
    assert tbl.pars[1].name == 'fract'

    assert tbl.nh.val == pytest.approx(1)
    assert tbl.nh.min == pytest.approx(0)
    assert tbl.nh.max == pytest.approx(1000)

    assert tbl.fract.val == pytest.approx(0.5)
    assert tbl.fract.min == pytest.approx(0)
    assert tbl.fract.max == pytest.approx(1)

    for p in tbl.pars:
        assert not(p.frozen)


@requires_data
@requires_fits
@requires_xspec
def test_evaluate_xspec_multiplicative_model(make_data_path):
    """Can we evaluate multiplicative table models?

    This is a limited test - in that it does not attempt to
    test the full set of grid inputs that we do with additive
    table models (and other XSPEC models) - as it is assumed that
    this logic hsa been tested.
    """

    from sherpa.astro import xspec

    if xspec.get_xsversion().startswith('12.10.0'):
        pytest.skip('Test known to crash XSPEC 12.10.0')

    path = make_data_path('testpcfabs.mod')
    tbl = xspec.read_xstable_model('bar', path)

    # This extends beyond the range of the model grid
    egrid = np.arange(0.1, 17, 1.0)
    elo = egrid[:-1]
    ehi = egrid[1:]

    # The expected values, evaluated with XSPEC 12.10.1b using
    # C++ code (i.e. not the Sherpa interface).
    #
    # It appears that the -1 is 1 in earlier versions.
    #
    yexp = np.asarray([0.511674,
                       0.730111,
                       0.898625,
                       0.95572,
                       0.977472,
                       0.987328,
                       0.992138,
                       0.990245,
                       0.992846,
                       0.994674,
                       0.995997,
                       0.996945,
                       0.997616,
                       0.998104,
                       0.998454,
                       -1])

    # Note, xspec 12.10.0 should not be seen here as explicitly
    # excluded above.
    xver = xspec.get_xsversion()
    if xver.startswith('12.9.'):
        yexp[-1] = 1.0

    y = tbl(elo, ehi)

    assert_almost_equal(y, yexp, decimal=6)


@requires_xspec
def test_ismabs_parameter_name_clashes():
    """Check the work around for the ismabs XSPEC 12.9.1 name clashes.

    The model.dat for ismabs has parameter names SiI and SII, which
    refer to different parameters (also SiII and SIII), but which
    Sherpa would treat as linked parameters. This test is provided to
    make sure that any documentation/code is updated if the chosen
    scheme to address this is updated (it is technically not needed, but is
    left in as a check that any future auto-generated XSPEC model
    handles these parameter names).

    Since this test does not evaluate the model, it should run with XSPEC
    12.9.0.
    """

    from sherpa.astro import xspec

    mdl = xspec.XSismabs()
    assert len(mdl.pars) == 31

    # List of expected names taken from XSPEC 12.9.1 model.dat file.
    #
    names = ["H", "HeII"]
    for el in ["C", "N", "O", "Ne", "Mg", "Si_", "S_", "Ar", "Ca"]:
        for i in ["I", "II", "III"]:
            names.append(el + i)
    names.extend(["Fe", "redshift"])
    assert len(names) == 31  # this tests the test, not the module!

    for par, expected in zip(mdl.pars, names):
        assert par.name == expected

        # Just check that there is no link between any of the parameters,
        # as would be the case if they were called SiI and SII (for example).
        assert par.link is None

    # It would be nice to be able to say the following, but at present
    # not sure how to enable this.
    #
    for name in ["SiI", "SII", "siii"]:
        with pytest.raises(AttributeError):
            getattr(mdl, name)


@requires_data
@requires_fits
@requires_xspec
def test_xstbl_link_parameter_evaluation(make_data_path):
    """See also sherpa/models/test_parameter::test_link_parameter_setting

    This is meant to replicate issue #742, where we want to ensure
    that parameter limits are respected, otherwise it is likely that
    the model evaluation will crash (since the table-model code will
    likely be indexing into unalocated memory).

    DJB has checked that this code causes a segfault on linux without
    a fix for #742.
    """

    from sherpa.astro import xspec

    path = make_data_path('xspec-tablemodel-RCS.mod')
    tbl = xspec.read_xstable_model('bar', path)

    # The tau parameter (first one) has a range of 1 to 10
    # - safety check that this still holds, so we know
    # that we are violating this limit when we set lmdl.c0
    # to 20
    #
    assert tbl.tau.min == pytest.approx(1)
    assert tbl.tau.max == pytest.approx(10)

    lmdl = Const1D()

    grid = np.arange(1, 6)
    glo = grid[:-1]
    ghi = grid[1:]

    tbl.tau = lmdl.c0
    lmdl.c0 = 2

    # just a safety check that we can change the parameter via
    # a link and run the model
    assert tbl.tau.val == pytest.approx(2)
    y2 = tbl(glo, ghi)
    assert (y2 > 0).all()

    # Test the fix for #742
    lmdl.c0 = 20
    emsg = 'parameter bar.tau has a maximum of 10'
    with pytest.raises(ParameterErr, match=emsg):
        tbl(glo, ghi)


@requires_xspec
@pytest.mark.parametrize("clsname", ["powerlaw", "wabs"])
def test_integrate_setting(clsname):
    """Can we change the integrate setting?

    It's not obvious what the integrate setting is meant to do for
    XSPEC models, so let's check what we can do with it.

    """

    from sherpa.astro import xspec

    egrid = np.arange(0.1, 1.0, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]

    mdl = getattr(xspec, 'XS{}'.format(clsname))()
    assert mdl.integrate

    y1 = mdl(elo, ehi)

    mdl.integrate = False
    assert not mdl.integrate

    y2 = mdl(elo, ehi)

    # Assume the integrate setting is ignored. To ensure we
    # can test small values use the log of the data, and
    # as we can have zero values (from xswabs) replace them
    # with a sentinel value
    #
    y1[y1 <= 0] = 1e-10
    y2[y2 <= 0] = 1e-10

    y1 = np.log10(y1)
    y2 = np.log10(y2)
    assert y2 == pytest.approx(y1)


@requires_xspec
@pytest.mark.parametrize("clsname", ["powerlaw", "wabs"])
def test_integrate_setting_con(clsname):
    """Can we change the integrate setting of convolution models.

    This is test_integrate_setting after wrapping the model by a
    convolution model. At present the convolution model does not have
    an integrate setting.

    Note that the convolved model doesn't make much sense physically -
    at least when it's xscflux(xswabs) - but we just care about the
    evaluation process here.

    """

    from sherpa.astro import xspec

    egrid = np.arange(0.1, 1.0, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]

    conv = xspec.XScflux()
    assert conv.integrate

    omdl = getattr(xspec, 'XS{}'.format(clsname))()
    assert omdl.integrate

    # Convolution models do not have an integrate setting
    mdl = conv(omdl)
    with pytest.raises(AttributeError):
        mdl.integrate

    y1 = mdl(elo, ehi)

    # as mdl does not have an integrate setting, just change
    # omdl
    omdl.integrate = False
    assert not omdl.integrate

    y2 = mdl(elo, ehi)

    # Assume the integrate setting is ignored. To ensure we
    # can test small values use the log of the data, and
    # as we can have zero values (from xswabs) replace them
    # with a sentinel value
    #
    y1[y1 <= 0] = 1e-10
    y2[y2 <= 0] = 1e-10

    y1 = np.log10(y1)
    y2 = np.log10(y2)
    assert y2 == pytest.approx(y1)


@requires_xspec
@pytest.mark.parametrize("val", [1.5, 7.5])
def test_xsparameter_within_limit(val):
    """What happens if we pass outside soft but within hard range?

    """

    from sherpa.astro.xspec import XSParameter

    p = XSParameter('temp', 'p', 4, min=2, max=7, hard_min=1, hard_max=8)
    assert not p.frozen

    p.val = val
    assert p.val == pytest.approx(val)
    assert not p.frozen


@requires_xspec
@pytest.mark.parametrize("val", [1, 8])
def test_xsparameter_at_limit(val):
    """What happens if we pass outside soft but within hard range?

    """

    from sherpa.astro.xspec import XSParameter

    p = XSParameter('temp', 'p', 4, min=2, max=7, hard_min=1, hard_max=8)

    p.val = val
    assert p.val == pytest.approx(val)
    assert not p.frozen


@requires_xspec
@pytest.mark.parametrize("base", [False, True])
@pytest.mark.parametrize("val", [0, 10])
def test_xsparameter_exceed_limit(base, val):
    """What happens if we set a value outside the hard range?

    """

    from sherpa.astro.xspec import XSBaseParameter, XSParameter

    cls = XSBaseParameter if base else XSParameter
    p = cls('temp', 'p', 4, min=2, max=7, hard_min=1, hard_max=8)
    assert not p.frozen

    with pytest.raises(ParameterErr):
        p.val = val

    assert p.val == pytest.approx(4)
    assert not p.frozen


@requires_xspec
@pytest.mark.parametrize("val", [0, 10])
def test_xsparameter_change_limit_larger(val, caplog):
    """Can we change the hard limits to increase the range?

    """

    from sherpa.astro.xspec import XSParameter

    p = XSParameter('temp', 'p', 4, min=2, max=7, hard_min=1, hard_max=8)

    p.hard_min = 0
    p.hard_max = 10
    p.val = val

    assert len(caplog.records) == 0

    assert p.val == pytest.approx(val)
    assert not p.frozen


@requires_xspec
@pytest.mark.parametrize("label,limit", [("min", 3), ("max", 6)])
def test_xsparameter_change_limit_smaller(label, limit, caplog):
    """Can we change the hard limits to decrease the range?

    The new range includes the current value.
    """

    from sherpa.astro.xspec import XSParameter

    p = XSParameter('temp', 'p', 5, min=2, max=7, hard_min=1, hard_max=8)

    assert len(caplog.records) == 0

    setattr(p, f'hard_{label}', limit)

    assert len(caplog.records) == 0

    assert p.val == pytest.approx(5)
    assert not p.frozen
    assert getattr(p, f'hard_{label}') == pytest.approx(limit)
    assert getattr(p, f'{label}') == pytest.approx(limit)


@requires_xspec
@pytest.mark.parametrize("label,limit", [("min", 6), ("max", 4)])
def test_xsparameter_change_limit_smaller_warning(label, limit, caplog):
    """Can we change the hard limits to decrease the range?

    We want to check what happens when the new limit range
    does not include the current value.
    """

    from sherpa.astro.xspec import XSParameter

    p = XSParameter('temp', 'p', 5, min=2, max=7, hard_min=1, hard_max=8)

    assert len(caplog.records) == 0

    setattr(p, f'hard_{label}', limit)

    assert len(caplog.records) == 1

    assert p.val == pytest.approx(limit)
    assert not p.frozen
    assert getattr(p, f'hard_{label}') == pytest.approx(limit)
    assert getattr(p, f'{label}') == pytest.approx(limit)

    if label == 'min':
        l1 = 'less'
        l2 = 'minimum'
    else:
        l1 = 'greater'
        l2 = 'maximum'

    expected = f"parameter temp.p {l1} than new {l2}; reset to {float(limit)}"

    module, level, msg = caplog.record_tuples[0]
    assert module == 'sherpa.astro.xspec'
    assert level == logging.WARNING
    assert msg == expected


@requires_xspec
@pytest.mark.parametrize("val", [0, 10])
def test_xsparameter_set(val, caplog):
    """Can we use set?"""

    from sherpa.astro.xspec import XSParameter

    p = XSParameter('temp', 'p', 4, min=2, max=7, hard_min=1, hard_max=8)

    p.set(hard_min=0, hard_max=10, val=val)

    assert p.val == pytest.approx(val)
    assert p.hard_min == pytest.approx(0)
    assert p.min == pytest.approx(0)
    assert p.hard_max == pytest.approx(10)
    assert p.max == pytest.approx(10)
    assert len(caplog.records) == 0


@requires_xspec
@pytest.mark.parametrize("label,limit", [("min", 6), ("max", 4)])
def test_xsparameter_set_change_val(label, limit, caplog):
    """Can we use set?"""

    from sherpa.astro.xspec import XSParameter

    p = XSParameter('temp', 'p', 5, min=2, max=7, hard_min=1, hard_max=8)

    kwargs = {f'hard_{label}': limit}
    p.set(**kwargs)

    assert len(caplog.records) == 1

    assert p.val == pytest.approx(limit)
    if label == 'min':
        assert p.hard_min == pytest.approx(limit)
        assert p.min == pytest.approx(limit)
        assert p.hard_max == pytest.approx(8)
        assert p.max == pytest.approx(8)
    else:
        assert p.hard_min == pytest.approx(1)
        assert p.min == pytest.approx(1)
        assert p.hard_max == pytest.approx(limit)
        assert p.max == pytest.approx(limit)

    if label == 'min':
        l1 = 'less'
        l2 = 'minimum'
    else:
        l1 = 'greater'
        l2 = 'maximum'

    expected = f"parameter temp.p {l1} than new {l2}; reset to {float(limit)}"

    module, level, msg = caplog.record_tuples[0]
    assert module == 'sherpa.astro.xspec'
    assert level == logging.WARNING
    assert msg == expected


@requires_xspec
@pytest.mark.parametrize("base", [True, False])
def test_xsparameter_limits(base):
    """Do we record the limits correctly?"""

    from sherpa.astro.xspec import XSBaseParameter, XSParameter

    cls = XSBaseParameter if base else XSParameter

    p = cls('temp', 'p', 4, min=2, max=7, hard_min=1, hard_max=8)

    # Note that there's no "API" for this, just these attributes
    assert p._xspec_soft_min == pytest.approx(2)
    assert p._xspec_soft_max == pytest.approx(7)

    assert p.min == pytest.approx(1)
    assert p.hard_min == pytest.approx(1)

    assert p.max == pytest.approx(8)
    assert p.hard_max == pytest.approx(8)
