#
#  Copyright (C) 2016, 2017  Smithsonian Astrophysical Observatory
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

from tempfile import NamedTemporaryFile

import pytest
import six

from numpy.testing import assert_almost_equal

from sherpa.utils import requires_data, requires_fits, requires_xspec


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
DEFAULT_ABUND = 'angr'
DEFAULT_XSECT = 'bcmc'

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
    as xfail.

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
    assert isinstance(v, six.string_types)
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
    assert oval == DEFAULT_ABUND


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
    assert oval == DEFAULT_XSECT


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
    assert getfunc() == oval


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
