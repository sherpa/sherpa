#
#  Copyright (C) 2016 - 2021, 2023, 2024
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
import re

import pytest

import numpy as np

from sherpa.models.basic import Const1D
from sherpa.utils.testing import requires_data, requires_fits, requires_xspec
from sherpa.utils.err import ArgumentErr, ParameterErr


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
# values. So, we now just check that the default values are
# one of the expected values.
#
# Note that a used could set up their own abundance table,
# in which case it is not obvious what to do.
#
DEFAULT_ABUND = ['angr', 'aspl', 'feld', 'aneb', 'grsa', 'wilm', 'lodd']
DEFAULT_XSECT = ['bcmc', 'obcm', 'vern']

# XSPEC defaults - we now set the chatter to 10 and the cosmology is
# fixed (although it may get set given some discussions DJB has had
# with the XSPEC developers circa XSPEC 12.12.0).
#
DEFAULT_CHATTER = 10
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


@requires_xspec
def test_chatter_default():
    """Check the expected default setting for chatter.

    Ideally this test would be run before any other
    tests of XSPEC are made (i.e. any XSPEC code is called).

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
def test_expected_elements():
    """If this fails then something is wrong!"""

    from sherpa.astro import xspec

    # At the moment this just checks I can repeat the information in
    # two places, but the aim is to change the interface to access
    # this information from XSPEC itself, at which point it becomes a
    # regression test.
    #
    xselems = xspec.get_xselements()
    for idx, elem in enumerate(ELEMENT_NAMES, 1):
        assert xselems[elem] == idx


@requires_xspec
def test_abund_angr_doc():
    """Check the doc string for angr.

    This is assumed to be constant.
    """

    from sherpa.astro import xspec

    doc = xspec.get_xsabund_doc("angr")
    assert doc == "Anders E. & Grevesse N. Geochimica et Cosmochimica Acta 53, 197 (1989)"


@requires_xspec
def test_abund_selected_doc():
    """Check the doc string for the selected table.

    This is assumed to be constant.
    """

    from sherpa.astro import xspec

    oval = xspec.get_xsabund()
    try:
        xspec.set_xsabund("lodd")
        doc = xspec.get_xsabund_doc()

    finally:
        xspec.set_xsabund(oval)

    assert doc == "Lodders, K. ApJ 591, 1220 (2003)"


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


def check_abundances(h, he, si, ar, k, fe):
    """Check wilm abundances.

    These values were found from HEASOFT version 6.19
    spectral/manager/abundances.dat

    The values are given to two decimal places in this file.
    It is not worth testing all settings, since we are not
    testing the XSPEC implementation itself, just our use of it.

    """

    assert h == pytest.approx(1.0)
    assert he == pytest.approx(9.77e-2)
    assert si == pytest.approx(1.86e-05)
    assert ar == pytest.approx(2.57e-06)
    assert k == pytest.approx(0.0)
    assert fe == pytest.approx(2.69e-05)


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

    check_abundances(h, he, si, ar, k, fe)


@requires_xspec
def test_abund_get_invalid_element(caplog):
    """Check what happens if sent the wrong element name"""

    from sherpa.astro import xspec

    # TODO: TypeError is not the best error type here.
    with pytest.raises(TypeError,
                       match="^could not find element 'O3'$"):
        xspec.get_xsabund("O3")

    assert len(caplog.records) == 0


@requires_xspec
def test_abund_set_invalid_name(caplog):
    """Check what happens if sent an unknown table

    It is unlikely that the name "foo-foo" will become valid.
    """

    from sherpa.astro import xspec

    with pytest.raises(ValueError,
                       match="^Cannot read file 'foo-foo'.  It may not exist or contains invalid data$"):
        xspec.set_xsabund("foo-foo")

    assert len(caplog.records) == 0


@requires_xspec
def test_xsect_set_invalid_name(caplog):
    """Check what happens if sent an unknown table

    It is unlikely that the name "foo-foo" will become valid.
    """

    from sherpa.astro import xspec

    with pytest.raises(ValueError,
                       match="^could not set XSPEC photoelectric cross-section to 'foo-foo'$"):
        xspec.set_xsxsect("foo-foo")

    assert len(caplog.records) == 0


@requires_xspec
def test_abund_get_dict():
    """Can we access the elemental settings?

    Make the same checks as test_abund_element
    """

    from sherpa.astro import xspec

    oval = xspec.get_xsabund()
    assert oval != 'wilm'
    try:
        xspec.set_xsabund('wilm')
        abunds = xspec.get_xsabundances()

    finally:
        xspec.set_xsabund(oval)

    check_abundances(abunds['H'],
                     abunds['He'],
                     abunds['Si'],
                     abunds['Ar'],
                     abunds['K'],
                     abunds['Fe'])


@requires_xspec
def test_abund_set_dict():
    """Can we set a dict of abundances?

    We only set a sub-set of elements, with the rest being set to 0.

    """

    from sherpa.astro import xspec

    oval = xspec.get_xsabund()

    # This is a pre-condition. It isn't really needed, but it is a better
    # test if it holds.
    #
    #
    assert oval != "file"

    abundances = {'He': 0.1, 'H': 1.0, 'O': 0.2, 'Fe': 0.3, 'Mn': 0.4}
    try:
        xspec.set_xsabundances(abundances)
        assert xspec.get_xsabund() == "file"

        got = []
        for elem in ELEMENT_NAMES:
            got.append(xspec.get_xsabund(elem))

    finally:
        xspec.set_xsabund(oval)

    assert got[0] == pytest.approx(1.0)
    assert got[1] == pytest.approx(0.1)
    assert got[7] == pytest.approx(0.2)
    assert got[25] == pytest.approx(0.3)
    assert got[24] == pytest.approx(0.4)

    for i in range(2, 7):
        assert got[i] == pytest.approx(0.0)

    for i in range(8, 24):
        assert got[i] == pytest.approx(0.0)

    for i in range(26, 30):
        assert got[i] == pytest.approx(0.0)

    assert xspec.get_xsabund() == oval


@requires_xspec
def test_abund_set_dict_invalid_key():
    """What happens if the key is not a number?"""

    from sherpa.astro import xspec

    ovals = xspec.get_xsabundances()

    abundances = {'He': 0.1, 'H': 1.0, 'O': 0.2, 'Fe': 0.3, 'Mn': 0.4,
                  'FooBar': 'value is ignored',
                  'Cr': 1.2}

    with pytest.raises(ArgumentErr,
                       match="^Invalid element name: 'FooBar'$"):
        xspec.set_xsabundances(abundances)

    # table should not have changed
    nvals = xspec.get_xsabundances()
    assert nvals == ovals


@requires_xspec
def test_abund_set_dict_invalid_value():
    """What happens if the value is not a number?"""

    from sherpa.astro import xspec

    ovals = xspec.get_xsabundances()

    # The abundances dict is processed in the order the keys are added,
    # so 'Cr' (invalid value) is processed by set_xsabundances before
    # 'FooBar' (invalid key).
    #
    abundances = {'He': 0.1, 'H': 1.0, 'O': 0.2, 'Fe': 0.3, 'Mn': 0.4,
                  'Cr': '1.2s', # number with text after it
                  'FooBar': 'value is ignored'}

    # This error comes from Sherpa (actually NumPy) and so the text
    # could change with NumPy version. Assume that to be unlikely
    # until we find it is a problem. This error path may change once
    # we can use the FunctionUtility interface.
    #
    with pytest.raises(ValueError,
                       match="^could not convert string to float: '1.2s'$"):
        xspec.set_xsabundances(abundances)

    # table should not have changed
    nvals = xspec.get_xsabundances()
    assert nvals == ovals


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
def test_abund_change_file(tmp_path):
    """Can we change the abundance setting: file

    This test hard-codes the number of elements expected in the
    file.
    """

    from sherpa.astro import xspec

    elems = {n: i * 0.1 for i, n in enumerate(ELEMENT_NAMES)}

    out = ""
    for n in ELEMENT_NAMES:
        out += f"{elems[n]}\n"

    tempfile = tmp_path / "abundances.xspec"
    tempfile.write_text(out)

    oval = xspec.get_xsabund()
    try:
        xspec.set_xsabund(str(tempfile))

        abund = xspec.get_xsabund()
        out = {n: xspec.get_xsabund(n)
               for n in ELEMENT_NAMES}
        out2 = xspec.get_xsabundances()

    finally:
        xspec.set_xsabund(oval)

    assert abund == 'file'
    for n in ELEMENT_NAMES:
        assert out[n] == pytest.approx(elems[n])

    assert out2 == out


@requires_xspec
def test_abund_change_file_subset(tmp_path):
    """What happens if send in too-few elements?"""

    from sherpa.astro import xspec

    elems = {n: i * 0.1 for i, n in enumerate(ELEMENT_NAMES)
             if i < 10}

    tmpname = tmp_path / "abunds.xspec"
    with open(tmpname, "w") as tfh:
        for v in elems.values():
            tfh.write(f"{v}\n")

    oval = xspec.get_xsabund()
    try:
        xspec.set_xsabund(str(tmpname))

        abund = xspec.get_xsabund()
        out = {n: xspec.get_xsabund(n)
               for n in ELEMENT_NAMES}
        out2 = xspec.get_xsabundances()

    finally:
        xspec.set_xsabund(oval)

    assert abund == 'file'
    for i, n in enumerate(ELEMENT_NAMES):
        if i < 10:
            assert out[n] == pytest.approx(elems[n])
        else:
            assert out[n] == pytest.approx(0)

    assert out2 == out


@requires_xspec
def test_abund_change_file_subset(tmp_path):
    """What happens if send in too-few elements?"""

    from sherpa.astro import xspec

    elems = {n: i * 0.1 for i, n in enumerate(ELEMENT_NAMES)
             if i < 10}

    tmpname = tmp_path / "abunds.xspec"
    with open(tmpname, "w") as tfh:
        for v in elems.values():
            tfh.write(f"{v}\n")

    oval = xspec.get_xsabund()
    try:
        xspec.set_xsabund(str(tmpname))

        abund = xspec.get_xsabund()
        out = {n: xspec.get_xsabund(n)
               for n in ELEMENT_NAMES}

    finally:
        xspec.set_xsabund(oval)

    assert abund == 'file'
    for i, n in enumerate(ELEMENT_NAMES):
        if i < 10:
            assert out[n] == pytest.approx(elems[n])
        else:
            assert out[n] == pytest.approx(0)


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
def test_path_manager_change(tmp_path):
    """Can we change the manager-path setting?
    """

    from sherpa.astro import xspec

    validate_xspec_setting(xspec.get_xspath_manager,
                           xspec.set_xspath_manager,
                           '/dev/null',
                           str(tmp_path))


# Note that the XSPEC state is used in test_xspec.py, but only
# to save/restore the state after each test. There is no
# explicit test there of the functionality. The state tests here
# are very basic.
#

@requires_xspec
@pytest.mark.parametrize("arg", [{}, True])
def test_set_xsstate_no_op(arg):
    """Just check we can sent in a "useless" argument.

    All we do is check the call can be made, we do not check whether
    anything has changed because of it.

    """

    from sherpa.astro import xspec

    xspec.set_xsstate(arg)


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

    def check(par, val, minval, maxval):
        assert par.val == pytest.approx(val)
        assert par.min == pytest.approx(minval)
        assert par.max == pytest.approx(maxval)

    check(tbl.tau, 1, 1, 10)
    check(tbl.beta, 0.1, 0.1, 0.5)
    check(tbl.t, 0.1, 0.1, 1.3)
    check(tbl.norm, 1, 0, 1e24)

    for p in tbl.pars:
        assert not(p.frozen)


@requires_xspec
@pytest.mark.parametrize("clsname", ["powerlaw", "wabs"])
def test_xspec_model_requires_bins(clsname, xsmodel):
    """Ensure you can not call with a single grid for the energies.

    You used to be able to do this (in Sherpa 4.13 and earlier).
    """

    mdl = xsmodel(clsname)

    emsg = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg):
        mdl([0.1, 0.2, 0.3, 0.4])


@requires_xspec
@pytest.mark.parametrize("clsname", ["powerlaw", "wabs"])
def test_xspec_model_requires_bins_low_level(clsname, xsmodel):
    """Ensure you can not call with a single grid for the energies (calc).

    You used to be able to do this (in Sherpa 4.13 and earlier).
    """

    mdl = xsmodel(clsname)

    emsg = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg):
        mdl.calc([p.val for p in mdl.pars], [0.1, 0.2, 0.3, 0.4])


@requires_xspec
@pytest.mark.parametrize("clsname", ["powerlaw", "wabs"])
def test_xspec_model_requires_bins_very_low_level(clsname, xsmodel):
    """Check we can use a single grid for direct access (_calc).

    This is the flip side to test-xspec_model_requires_bins_low_level
    """

    mdl = xsmodel(clsname)

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
def test_xspec_tablemodel_requires_bin_edges(make_data_path):
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
def test_xspec_tablemodel_requires_bin_edges_low_level(make_data_path):
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

    # We get warnings from m1 evaluated on the grid and then the
    # convolution model m2.
    #
    emsg1 = r'calc\(\) requires pars,rhs,lo,hi arguments, sent 3 arguments'
    emsg2 = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg1):
        with pytest.warns(FutureWarning, match=emsg2):
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

    # We get warnings from m1 evaluated on the grid and then the
    # convolution model m2.
    #
    emsg1 = r'calc\(\) requires pars,rhs,lo,hi arguments, sent 3 arguments'
    emsg2 = r'calc\(\) requires pars,lo,hi arguments, sent 2 arguments'
    with pytest.warns(FutureWarning, match=emsg1):
        with pytest.warns(FutureWarning, match=emsg2):
            mdl.calc([p.val for p in mdl.pars], [0.1, 0.2, 0.3, 0.4])


@requires_data
@requires_fits
@requires_xspec
def test_evaluate_xspec_additive_model_beyond_grid(make_data_path):
    """Can we extend an additive table model beyond its grid?"""

    from sherpa.astro import xspec

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
    this logic has been tested.
    """

    from sherpa.astro import xspec

    path = make_data_path('testpcfabs.mod')
    tbl = xspec.read_xstable_model('bar', path)

    # This extends beyond the range of the model grid
    egrid = np.arange(0.1, 17, 1.0)
    elo = egrid[:-1]
    ehi = egrid[1:]

    # The expected values, evaluated with XSPEC 12.10.1b using
    # C++ code (i.e. not the Sherpa interface).
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

    y = tbl(elo, ehi)

    # The final element is -1 prior to 12.13.1 but 1 in 12.13.1,
    # (probably because of an issue we reported), so let's drop that
    # element for the check.
    #
    assert y[:-1] == pytest.approx(yexp[:-1])

    # Check the final value. Should we provide a nice way to
    # extrct the Sherpa version?
    #
    from sherpa.astro.xspec import get_xsversion
    version = get_xsversion()
    match = re.search(r'^(\d+)\.(\d+)\.(\d+)', version)
    if match is None:
        raise ValueError(f"Invalid XSPEC version string: {version}")

    v = (int(match[1]), int(match[2]), int(match[3]))
    if v < (12, 13, 1):
        assert y[-1] == pytest.approx(-1.0)
    else:
        assert y[-1] == pytest.approx(1.0)


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

    As of XSPEC 12.13.0 (ish), the model.dat file now lists the
    underscore version for all the parameters (e.g. He_II, C_I,
    Ca_III), which changes this test somewhat. Note that as we don't
    autogenerate the interface we have to decide on a naming scheme
    (we don't go to the effort of making it depend on the XSPEC
    version), and so we have decided to go with the XSPEC 12.13.1 /
    HEASOFT 6.32 names. We have aliased the old names, so we can check
    they work too.

    """

    from sherpa.astro import xspec

    mdl = xspec.XSismabs()
    assert len(mdl.pars) == 31

    # List of expected names taken from XSPEC 12.13.1 model.dat file.
    #
    names = ["H", "He_II"]
    for el in ["C", "N", "O", "Ne", "Mg", "Si", "S", "Ar", "Ca"]:
        for i in ["I", "II", "III"]:
            names.append(f"{el}_{i}")
    names.extend(["Fe", "redshift"])
    assert len(names) == 31  # this tests the test, not the module!

    for par, expected in zip(mdl.pars, names):
        assert par.name == expected

        # Just check that there is no link between any of the parameters,
        # as would be the case if they were called SiI and SII (for example).
        assert par.link is None

    # Check that the aliases work and raise a deprecation warning.
    #
    # Note that we include Si and S in the list, just so that we can
    # keep the ordering with the names array, but we skip these as
    # they do not have an alias (hence the name begins with SKIP).
    #
    aliases = ["HeII"]
    for el in ["C", "N", "O", "Ne", "Mg", "SKIP-Si", "SKIP-S", "Ar", "Ca"]:
        for i in ["I", "II", "III"]:
            aliases.append(f"{el}{i}")

    for name, alias in zip(names[1:], aliases):
        if alias.startswith("SKIP"):
            continue

        assert alias != name  # safety check

        emsg = f"^Parameter name {alias} is deprecated for " + \
            f"model XSismabs, use {name} instead$"
        with pytest.warns(DeprecationWarning, match=emsg):
            par = getattr(mdl, alias)

        assert par.name == name

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
def test_integrate_setting(clsname, xsmodel):
    """Can we change the integrate setting?

    It's not obvious what the integrate setting is meant to do for
    XSPEC models, so let's check what we can do with it.

    """

    egrid = np.arange(0.1, 1.0, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]

    mdl = xsmodel(clsname)
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
def test_integrate_setting_con(clsname, xsmodel):
    """Can we change the integrate setting of convolution models.

    This is test_integrate_setting after wrapping the model by a
    convolution model. At present the convolution model does not have
    an integrate setting.

    Note that the convolved model doesn't make much sense physically -
    at least when it's xscflux(xswabs) - but we just care about the
    evaluation process here.

    """

    egrid = np.arange(0.1, 1.0, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]

    conv = xsmodel("cflux")
    assert conv.integrate

    omdl = xsmodel(clsname)
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
    """What happens if we pass outside soft but at hard limit?

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


@requires_xspec
@pytest.mark.parametrize("clsname", ["XSpowerlaw", "XSphabs"])
def test_model_can_send_spectrumnumber_indiv(clsname):
    """Check we can send spectrumNumber to individual models.

    All we do is check that the model can be called, and just for
    additive/multiplicative models that are not expected to use
    the spectrumNumber argument.
    """

    from sherpa.astro import xspec

    mdl = getattr(xspec, clsname)("tmp")
    egrid = np.arange(0.3, 0.4, 0.01)
    # This fails
    mdl(egrid[:-1], egrid[1:], spectrumNumber=2)


@requires_xspec
def test_model_can_send_spectrumnumber_indiv_con():
    """Check we can send spectrumNumber to individual models: convolution

    All we do is check that the model can be called.
    """

    from sherpa.astro import xspec

    base = xspec.XSpowerlaw("tmp")
    con = xspec.XScflux("tmp2")
    mdl = con(base)
    egrid = np.arange(0.3, 0.4, 0.01)

    # Does this send the argument to the wrapped model? There's no way
    # to know when calling the actual XSPEC models.
    #
    mdl(egrid[:-1], egrid[1:], spectrumNumber=2)


@requires_xspec
def test_model_can_send_spectrumnumber_combine():
    """Check the spectrumNumber is sent through.

    Mock XSPEC model classes are used to check that we send through
    the spectrumNumber when combining models. This uses the
    "old-style" of declaring the model (explicit setting of _calc).

    Due to the way we have to write these tests (with xspec not
    imported at the top-level) this is not a unit test as we test
    a number of features.
    """

    from sherpa.astro import xspec
    from sherpa.models import Parameter

    # We use the first parameter to identify what model is being
    # called (relying on the test to change the index value of the
    # model components).
    #
    args = []
    def test(cls, pars, lo, hi, spectrumNumber=5):
        args.append((spectrumNumber, pars[0]))
        return pars[1] * np.ones_like(lo)

    def testcon(cls, pars, fluxes, lo, hi, spectrumNumber=9):
        args.append((spectrumNumber, "con"))
        return fluxes + pars[1]

    class TestSpectrumNumber(xspec.XSAdditiveModel):
        _calc = test

        def __init__(self, name="test"):
            self.index = Parameter(name, "index", 1, 0, 2)
            self.norm = Parameter(name, "norm", 1, 0, 1)
            super().__init__(name, (self.index, self.norm))

    class TestConvSpectrumNumber(xspec.XSConvolutionKernel):
        _calc = testcon

        def __init__(self, name="test"):
            self.index = Parameter(name, "index", 1, 0, 2)
            self.con = Parameter(name, "con", 1, 0, 100)
            super().__init__(name, (self.index, self.con))

    m1 = TestSpectrumNumber("x1")
    m1.index = 1
    m1.norm = 0.2
    m2 = TestSpectrumNumber("2")
    m2.index = 2
    m2.norm = 0.4

    comb = m1 + m2

    # As we evaluate the models multiple times with the same
    # arguments we need to ensure the models are not cached.
    #
    m1._use_caching = False
    m2._use_caching = False
    # comb._use_caching = False  no cache for composite models

    egrid = np.arange(1, 4)
    elo = egrid[:-1]
    ehi = egrid[1:]

    assert len(args) == 0
    y1 = m1(elo, ehi)
    assert y1 == pytest.approx([0.2, 0.2])
    assert len(args) == 1
    print(args)
    assert args[0][0] == 5
    assert args[0][1] == pytest.approx(1)

    args.clear()
    y2 = m2(elo, ehi, spectrumNumber=1)
    assert y2 == pytest.approx([0.4, 0.4])
    assert len(args) == 1
    assert args[0][0] == 1
    assert args[0][1] == pytest.approx(2)

    args.clear()
    y = comb(elo, ehi, spectrumNumber=3)
    assert y == pytest.approx([0.6, 0.6])

    # The order of evaluation should be guaranteed (i.e we can test
    # that [0][1] is 1).
    #
    assert len(args) == 2
    assert args[0][0] == 3
    assert args[0][1] == pytest.approx(1)
    assert args[1][0] == 3
    assert args[1][1] == pytest.approx(2)

    # Try with a convolution model.
    #
    con = TestConvSpectrumNumber("x")
    con.con = 10
    cmdl = con(comb)

    args.clear()
    ycon = cmdl(elo, ehi, spectrumNumber=6)
    assert ycon == pytest.approx([10.6, 10.6])

    # NOTE: these are not the answers we want, but test them so we know
    #       when they change.
    #
    assert len(args) == 3
    assert args[0][0] == 5  # should be 6
    assert args[0][1] == pytest.approx(1)
    assert args[1][0] == 5  # should be 6
    assert args[1][1] == pytest.approx(2)
    assert args[2][0] == 9  # should be 6
    assert args[2][1] == "con"


@requires_xspec
def test_model_can_send_spectrumnumber_combine_non_xspec():
    """Check the spectrumNumber is sent through with non-XSPEC models.

    An extension of test_model_can_send_spectrumnumber_combine()
    """

    from sherpa.astro import xspec
    from sherpa.models import Parameter
    from sherpa.models.basic import Scale1D

    args = []
    def test(cls, pars, lo, hi, spectrumNumber=5):
        args.append((spectrumNumber, pars[0]))
        return pars[1] * np.ones_like(lo)

    class TestSpectrumNumber2(xspec.XSAdditiveModel):
        _calc = test

        def __init__(self, name="test"):
            self.index = Parameter(name, "index", 1, 0, 2)
            self.norm = Parameter(name, "norm", 1, 0, 1)
            super().__init__(name, (self.index, self.norm))

    m1 = TestSpectrumNumber2("m1")
    m1.norm = 0.5
    m1.index = 2

    m2 = Scale1D("m2")
    m2.c0 = 1.5

    m1._use_caching = False
    m2._use_caching = False

    # These should evaluate to the same thing.
    #
    comb12 = m1 + m2
    comb21 = m2 + m1

    elo = np.asarray([0.5, 1, 2])
    ehi = np.asarray([1, 2, 4])

    assert len(args) == 0
    y12 = comb12(elo, ehi)
    assert len(args) == 1
    assert args[0][0] == 5
    assert args[0][1] == pytest.approx(2)

    args.clear()

    y21 = comb21(elo, ehi, spectrumNumber=4)
    assert len(args) == 1
    assert args[0][0] == 4
    assert args[0][1] == pytest.approx(2)

    assert y21 == pytest.approx(y12)
    assert y12 == pytest.approx([2, 2, 2])


@requires_xspec
def test_xspec_model_kwarg_cached_not_needed():
    """Does the cache recognize when a kwarg has changed?

    This is for a XSPEC model that is known to not depend on the
    spectrumNumber/ifl argument (i.e. the definition in model.dat has
    a second 0 after the model type, in this case saying "add 0 0" or
    "add 0" as the 0 is inferred).

    At the moment Sherpa does not wrap any of the models (smaug,
    polcost, pollin, polpow, pileup) which have a 1 instead of a 0.

    Technically we could make the cache ignore the spectrumNumber
    value, since it does not change the model output, but it is tricky
    to do, so we do make each call be different. This is a
    pessimisation (i.e. it potentially slows down model evaluation)
    and could be changed in the future. It depends on how often the
    spectrumNumber argument will be added to model calls.

    """

    from sherpa.astro import xspec

    mdl = xspec.XSpowerlaw()

    # check the cache is working
    assert mdl._cache_ctr['hits'] == 0
    assert mdl._cache_ctr['misses'] == 0

    egrid = np.arange(0.2, 0.6, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]

    y1 = mdl(elo, ehi)
    assert mdl._cache_ctr['hits'] == 0
    assert mdl._cache_ctr['misses'] == 1

    y2 = mdl(elo, ehi, spectrumNumber=2)
    assert mdl._cache_ctr['hits'] == 0
    assert mdl._cache_ctr['misses'] == 2

    assert y2 == pytest.approx(y1)

    y3 = mdl(elo, ehi, spectrumNumber=1)
    assert mdl._cache_ctr['hits'] == 0
    assert mdl._cache_ctr['misses'] == 3

    assert y3 == pytest.approx(y1)


@requires_xspec
@requires_data
@requires_fits
@pytest.mark.parametrize("addmodel", [0, 1])
@pytest.mark.parametrize("redshift", [0, 1])
@pytest.mark.parametrize("escale", [0, 1])
def test_table_mod_negative_delta_1850(addmodel, redshift, escale, make_data_path):
    """Is delta<0 recognized as freezing a parameter?

    The tables in xspec_table_models are numbered in such a way that
    we can loop through them all, and hence check some other things at
    the same time (to avoid having to read things in multiple times).

    """

    from sherpa.astro import xspec

    name = f"smod{addmodel}{redshift}{escale}.tmod"
    infile = make_data_path(f"xspec_table_models/{name}")

    tbl = xspec.read_xstable_model('tbl', infile)

    # issue #1850; is the parameter recognized as frozen
    assert tbl.lscale.frozen

    parnames = [p.name for p in tbl.pars]
    assert parnames[0] == "lscale"

    # Ordering taken from XSPEC 12.13.1a (unreleased at the time of
    # writing of the test).
    #
    idx = 1
    if escale:
        assert parnames[idx] == "Escale"
        assert tbl.escale.frozen
        idx += 1

    if redshift:
        assert parnames[idx] == "redshift"
        assert tbl.redshift.frozen
        idx += 1

    if addmodel:
        assert parnames[idx] == "norm"
        assert not tbl.norm.frozen
        idx += 1

    assert len(parnames) == idx


@requires_xspec
@requires_data
@requires_fits
def test_table_mod_add(make_data_path):
    """Check the additive model is behaving as expected. No z/escale

    """

    from sherpa.astro import xspec

    name = "smod100.tmod"
    infile = make_data_path(f"xspec_table_models/{name}")

    tbl = xspec.read_xstable_model('tbl', infile)

    egrid = np.asarray([0.2, 0.3, 0.5, 0.7, 1.1, 1.6, 2.0, 2.2, 2.5])
    elo = egrid[:-1]
    ehi = egrid[1:]

    expected = [0, 0, 22.5, 6, 54, 56, 39.5, 0]
    assert tbl(elo, ehi) == pytest.approx(expected)

    tbl.lscale = 0
    expected = [0, 0, 15, 12, 52, 48, 15, 0]
    assert tbl(elo, ehi) == pytest.approx(expected)


# Values taken from XSPEC 12.13.1 with
#
#    dummyrsp 0.2 2.6 24 linear
#    mo atable{smod1xx.tmod}
#    newpar xxxxxx
#    iplot model
#    wdata
#
# Note that this gives photon/cm^2/s/keV and model evaluation
# gives photon/cm^2/s.
#
ADD_TABLE_BASIC = [0, 0, 0, 75, 150] + [15] * 4 + [100] * 4 + \
    [140] * 5 + [375, 20, 0, 0, 0, 0]
ADD_TABLE_Z1 = [0, 82.5, 15, 57.5, 100, 120, 140, 140, 197.5] + [0] * 15
ADD_TABLE_E2 = [0] * 8 + [37.5] * 2 + [75] * 2 + [7.5] * 8 + [50] * 4


@requires_xspec
@requires_data
@requires_fits
def test_table_mod_add_redshift(make_data_path):
    """Check the additive model is behaving as expected with redshift. No escale.

    """

    from sherpa.astro import xspec

    name = "smod110.tmod"
    infile = make_data_path(f"xspec_table_models/{name}")

    tbl = xspec.read_xstable_model('tbl', infile)

    assert tbl.redshift.val == pytest.approx(0)

    egrid = np.linspace(0.2, 2.6, 25)
    elo = egrid[:-1]
    ehi = egrid[1:]
    de = 0.1

    assert tbl(elo, ehi) / de == pytest.approx(ADD_TABLE_BASIC)

    tbl.redshift = 1
    assert tbl(elo, ehi) / de == pytest.approx(ADD_TABLE_Z1, rel=2e-6)


@requires_xspec
@requires_data
@requires_fits
def test_table_mod_add_escale(make_data_path):
    """Check the additive model is behaving as expected with escale. No z.

    """

    from sherpa.astro import xspec

    name = "smod101.tmod"
    infile = make_data_path(f"xspec_table_models/{name}")

    tbl = xspec.read_xstable_model('tbl', infile)

    assert tbl.escale.val == pytest.approx(1.0)

    egrid = np.linspace(0.2, 2.6, 25)
    elo = egrid[:-1]
    ehi = egrid[1:]
    de = 0.1

    assert tbl(elo, ehi) / de == pytest.approx(ADD_TABLE_BASIC)

    tbl.escale = 2
    assert tbl(elo, ehi) / de == pytest.approx(ADD_TABLE_E2, rel=2e-6)


@requires_xspec
@requires_data
@requires_fits
def test_table_mod_add_escale_redshift(make_data_path):
    """Check the additive model is behaving as expected with escale+redshift.

    Note, prior to XSPEC 12.13.1a XSPEC had the escale and redshift
    parameters the wrong way round.

    """

    from sherpa.astro import xspec

    name = "smod111.tmod"
    infile = make_data_path(f"xspec_table_models/{name}")

    tbl = xspec.read_xstable_model('tbl', infile)

    assert tbl.escale.val == pytest.approx(1)
    assert tbl.redshift.val == pytest.approx(0)

    egrid = np.linspace(0.2, 2.6, 25)
    elo = egrid[:-1]
    ehi = egrid[1:]
    de = 0.1

    assert tbl(elo, ehi) / de == pytest.approx(ADD_TABLE_BASIC)

    tbl.redshift = 1
    assert tbl(elo, ehi) / de == pytest.approx(ADD_TABLE_Z1, rel=2e-6)

    tbl.redshift = 0
    tbl.escale = 2
    assert tbl(elo, ehi) / de == pytest.approx(ADD_TABLE_E2, rel=2e-6)


@requires_xspec
@requires_data
@requires_fits
def test_table_mod_mul(make_data_path):
    """Check the multiplicative model is behaving as expected. No z/escale

    """

    from sherpa.astro import xspec

    name = "smod000.tmod"
    infile = make_data_path(f"xspec_table_models/{name}")

    tbl = xspec.read_xstable_model('tbl', infile)

    egrid = np.asarray([0.2, 0.3, 0.5, 0.7, 1.1, 1.6, 2.0, 2.2, 2.5])
    elo = egrid[:-1]
    ehi = egrid[1:]

    expected = [0, 0, 11.25, 6, 45, 70, 19.75, 5]
    assert tbl(elo, ehi) == pytest.approx(expected)

    tbl.lscale = 0
    expected = [0, 0, 7.5, 12, 40 + 10/3, 60, 7.5, 5]
    assert tbl(elo, ehi) == pytest.approx(expected)


# Values taken from XSPEC 12.13.1 with
#
#    dummyrsp 0.2 2.6 24 linear
#    mo mtable{smod0xx.tmod}
#    newpar xxxxx
#    iplot model
#    wdata
#
# As this is multiplicative there is no need to worry about
# the bin width.
#
MUL_TABLE_BASIC = [0, 0, 0, 7.5, 15] + [6] * 4 + [40] * 4 + \
    [70] * 5 + [37.5, 2, 5, 5, 5, 5]
MUL_TABLE_Z1 = [0, 13.2, 6, 23, 40, 50 + 10/3, 70, 70, 19.75] + [5] * 15
MUL_TABLE_E2 = [0] * 8 + [7.5] * 2 + [15] * 2 + [6] * 8 + [40] * 4


@requires_xspec
@requires_data
@requires_fits
def test_table_mod_mul_redshift(make_data_path):
    """Check the multiplicative model is behaving as expected with redshift. No escale.

    """

    from sherpa.astro import xspec

    name = "smod010.tmod"
    infile = make_data_path(f"xspec_table_models/{name}")

    tbl = xspec.read_xstable_model('tbl', infile)

    assert tbl.redshift.val == pytest.approx(0)

    egrid = np.linspace(0.2, 2.6, 25)
    elo = egrid[:-1]
    ehi = egrid[1:]

    assert tbl(elo, ehi) == pytest.approx(MUL_TABLE_BASIC)

    tbl.redshift = 1
    assert tbl(elo, ehi) == pytest.approx(MUL_TABLE_Z1, rel=2e-6)


@requires_xspec
@requires_data
@requires_fits
def test_table_mod_mul_escale(make_data_path):
    """Check the multiplicative model is behaving as expected with escale. No z.

    """

    from sherpa.astro import xspec

    name = "smod001.tmod"
    infile = make_data_path(f"xspec_table_models/{name}")

    tbl = xspec.read_xstable_model('tbl', infile)

    assert tbl.escale.val == pytest.approx(1.0)

    egrid = np.linspace(0.2, 2.6, 25)
    elo = egrid[:-1]
    ehi = egrid[1:]

    assert tbl(elo, ehi) == pytest.approx(MUL_TABLE_BASIC)

    tbl.escale = 2
    assert tbl(elo, ehi) == pytest.approx(MUL_TABLE_E2, rel=2e-6)


@requires_xspec
@requires_data
@requires_fits
def test_table_mod_mul_escale_redshift(make_data_path):
    """Check the multiplicative model is behaving as expected with escale+redshift.

    Note, prior to XSPEC 12.13.1a XSPEC had the escale and redshift
    parameters the wrong way round.

    """

    from sherpa.astro import xspec

    name = "smod011.tmod"
    infile = make_data_path(f"xspec_table_models/{name}")

    tbl = xspec.read_xstable_model('tbl', infile)

    assert tbl.escale.val == pytest.approx(1.0)
    assert tbl.redshift.val == pytest.approx(0)

    egrid = np.linspace(0.2, 2.6, 25)
    elo = egrid[:-1]
    ehi = egrid[1:]

    assert tbl(elo, ehi) == pytest.approx(MUL_TABLE_BASIC)

    tbl.redshift = 1
    assert tbl(elo, ehi) == pytest.approx(MUL_TABLE_Z1, rel=2e-6)

    tbl.redshift = 0
    tbl.escale = 2
    assert tbl(elo, ehi) == pytest.approx(MUL_TABLE_E2, rel=2e-6)
