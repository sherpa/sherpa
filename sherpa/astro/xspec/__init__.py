#
#  Copyright (C) 2010, 2015-2018, 2019, 2020
#                Smithsonian Astrophysical Observatory
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

"""Support for XSPEC models.

Sherpa supports versions 12.10.1, 12.10.0, 12.9.1, and 12.9.0 of XSPEC [1]_,
and can be built against the model library or the full application.  There is
no guarantee of support for older or newer versions of XSPEC.

To be able to use most routines from this module, the HEADAS environment
variable must be set. The `get_xsversion` function can be used to return the
XSPEC version - including patch level - the module is using::

   >>> from sherpa.astro import xspec
   >>> xspec.get_xsversion()
   '12.10.1b'

Initializing XSPEC
------------------

The XSPEC model library is initalized so that the cosmology parameters
are set to H_0=70, q_0=0.0, and lambda_0=0.73 (they can be changed with
`set_xscosmo`).

The other settings - for example for the abundance and cross-section
tables - follow the standard rules for XSPEC. For XSPEC versions prior
to 12.10.1, this means that the abundance table uses the ``angr``
setting and the cross sections the ``bcmc`` setting (see `set_xsabund`
and `set_xsxsect` for full details). As of XSPEC 12.10.1, the values
are now taken from the user's XSPEC configuration file - either
``~/.xspec/Xspec.init`` or ``$HEADAS/../spectral/manager/Xspec.init`` -
for these settings. The default value for the photo-ionization table
in this case is now ``vern`` rather than ``bcmc``.

References
----------

.. [1] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/index.html

"""

from __future__ import absolute_import

import string
from sherpa.models import Parameter, modelCacher1d, RegriddableModel1D
from sherpa.models.parameter import hugeval

from sherpa.utils import guess_amplitude, param_apply_limits, bool_cast
from sherpa.utils.err import ParameterErr
from sherpa.astro.utils import get_xspec_position

from .utils import ModelMeta, version_at_least, equal_or_greater_than
from . import _xspec


# Python wrappers around the exported functions from _xspec. This
# provides a more-accurate function signature to the user, makes
# the documentation easier to write, and makes it available even
# when the compiled code has not been compiled (e.g. for a Sphinx
# documentation run).
#
def get_xsabund(element=None):
    """Return the X-Spec abundance setting or elemental abundance.

    Parameters
    ----------
    element : str, optional
       When not given, the abundance table name is returned.
       If a string, then it must be an element name from:
       'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
       'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
       'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
       'Cu', 'Zn'. Case is important.

    Returns
    -------
    val : str or float
       When ``element`` is ``None``, the abundance table name is
       returned (see `set_xsabund`); the string 'file' is used
       when the abundances were read from a file. A numeric value
       is returned when an element name is given. This value is the
       elemental abundance relative to H.

    See Also
    --------
    set_xsabund

    Examples
    --------

    Return the current abundance setting, which in this case
    is 'angr', the default value for X-Spec:

    >>> get_xsabund()
    'angr'

    The `set_xsabund` function has been used to read in the
    abundances from a file, so the routine now returns the
    string 'file':

    >>> set_xsabund('abund.dat')
    >>> get_xsabund()
    'file'

    >>> get_xsabund('He')
    0.09769999980926514

    """

    if element is None:
        return _xspec.get_xsabund()
    else:
        return _xspec.get_xsabund(element)


def get_xschatter():
    """Return the chatter level used by X-Spec.

    Returns
    -------
    chatter : int
       The chatter setting used by the X-Spec routines.

    See Also
    --------
    set_xschatter

    Examples
    --------

    >>> get_xschatter()
    0
    """

    return _xspec.get_xschatter()


def get_xscosmo():
    """Return the X-Spec cosmology settings.

    Returns
    -------
    (h0, q0, l0)
       The Hubble constant, in km/s/Mpc, the deceleration parameter,
       and the cosmological constant.

    See Also
    --------
    set_xscosmo

    Examples
    --------

    >>> get_xscosmo()
    (70.0, 0.0, 0.7300000190734863)
    """

    return _xspec.get_xscosmo()


def get_xsversion():
    """Return the version of the X-Spec model library in use.

    Returns
    -------
    version : str
       The version of the X-Spec model library used by Sherpa [1]_.

    References
    ----------

    .. [1] http://heasarc.nasa.gov/docs/xanadu/xspec/

    Examples
    --------

    >>> get_xsversion()
    '12.9.1p'
    """

    return _xspec.get_xsversion()


def get_xsxsect():
    """Return the cross sections used by X-Spec models.

    Returns
    -------
    val : str
       The value of the photoelectric absorption setting: one of
       'bcmc', 'obcm', and 'vern'.

    See Also
    --------
    set_xsxsect

    Examples
    --------

    >>> get_xsxsect()
    'bcmc'
    """

    return _xspec.get_xsxsect()


def set_xsabund(abundance):
    """Set the elemental abundances used by X-Spec models.

    Set the abundance table used in the X-Spec plasma emission and
    photoelectric absorption models. It is equivalent to the X-Spec
    ``abund`` command [1]_.

    Parameters
    ----------
    abundance : str
       A file name, format described below, or one of the pre-defined
       names listed in the Notes section below.

    See Also
    --------
    get_xsabund, get_xsversion, set_xschatter

    Notes
    -----
    The pre-defined abundance tables are:
     - 'angr', from [2]_
     - 'aspl', from [3]_
     - 'feld', from [4]_, except for elements not listed which
       are given 'grsa' abundances
     - 'aneb', from [5]_
     - 'grsa', from [6]_
     - 'wilm', from [7]_, except for elements not listed which
       are given zero abundance
     - 'lodd', from [8]_

    The values for these tables are given at [1]_.

    Data files should be in ASCII format, containing a single
    numeric (floating-point) column of the abundance values,
    relative to Hydrogen.

    The screen output of this function is controlled by the
    X-Spec chatter setting (`set_xschatter`).

    References
    ----------
    .. [1] http://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html
           Note that this may refer to a newer version than the
           compiled version used by Sherpa; use `get_xsversion` to
           check.

    .. [2] Anders E. & Grevesse N. (1989, Geochimica et
           Cosmochimica Acta 53, 197)
           http://adsabs.harvard.edu/abs/1989GeCoA..53..197A

    .. [3] Asplund M., Grevesse N., Sauval A.J. & Scott P.
           (2009, ARAA, 47, 481)
           http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A

    .. [4] Feldman U.(1992, Physica Scripta 46, 202)
           http://adsabs.harvard.edu/abs/1992PhyS...46..202F

    .. [5] Anders E. & Ebihara (1982, Geochimica et Cosmochimica
           Acta 46, 2363)
           http://adsabs.harvard.edu/abs/1982GeCoA..46.2363A

    .. [6] Grevesse, N. & Sauval, A.J. (1998, Space Science
           Reviews 85, 161)
           http://adsabs.harvard.edu/abs/1998SSRv...85..161G

    .. [7] Wilms, Allen & McCray (2000, ApJ 542, 914)
           http://adsabs.harvard.edu/abs/2000ApJ...542..914W

    .. [8] Lodders, K (2003, ApJ 591, 1220)
           http://adsabs.harvard.edu/abs/2003ApJ...591.1220L

    Examples
    --------

    >>> set_xsabund('lodd')
     Solar Abundance Vector set to lodd:  Lodders, K. ApJ 591, 1220 (2003)

    >>> set_xsabund('abund.dat')
     Solar Abundance Vector set to file:  User defined abundance vector / no description specified

    """

    _xspec.set_xsabund(abundance)


def set_xschatter(level):
    """Set the chatter level used by X-Spec.

    Set the chatter setting used by the X-Spec routines for determining
    what information gets printed to the screen. It is equivalent to
    the X-Spec ``chatter`` command [1]_.

    Parameters
    ----------
    level : int
       The higher the value of ``level``, the more screen output will
       be created by X-Spec routines. A value of ``0`` hides
       most information while ``25`` will generate a lot of
       debug output.

    See Also
    --------
    get_xschatter, get_xsversion

    Notes
    -----
    The default chatter setting used by Sherpa is ``0``, which is lower
    than - so, creates less screen output - the default value used by
    X-Spec (``10``).

    There is no way to change the X-Spec "log chatter" setting.

    References
    ----------
    .. [1] http://heasarc.nasa.gov/xanadu/xspec/manual/XSchatter.html
           Note that this may refer to a newer version than the
           compiled version used by Sherpa; use `get_xsversion` to
           check.

    Examples
    --------

    Set the chatter level to the default used by X-Spec:

    >>> set_xschatter(10)
    """

    _xspec.set_xschatter(level)


def set_xscosmo(h0, q0, l0):
    """Set the cosmological parameters used by X-Spec models.

    Set the cosmological parameters (H_0, q_0, lambda_0) used
    by X-Spec. It is equivalent to the X-Spec ``cosmo``
    command [1]_. The default values are h0=70, q0=0, and l0=0.73

    Parameters
    ----------
    h0 : number
       The Hubble constant in km/s/Mpc.
    q0 : number
       The deceleration parameter.
    l0 : number
       The cosmological constant. If this is non-zero then the q0
       parameter is ignored and the Universe is assumed to be flat.

    See Also
    --------
    get_xscosmo, get_xsversion

    References
    ----------
    .. [1] http://heasarc.nasa.gov/xanadu/xspec/manual/XScosmo.html
           Note that this may refer to a newer version than the
           compiled version used by Sherpa; use `get_xsversion` to
           check.

    Examples
    --------

    >>> set_xscosmo(73, 0, 0.73)
    """

    _xspec.set_xscosmo(h0, q0, l0)


def set_xsxsect(name):
    """Set the cross sections used by X-Spec models.

    Set the X-Spec photoelectric absorption cross-sections
    setting, which changes the cross-sections used by all
    X-Spec absorption models *except* for `XSwabs`. It is
    equivalent to the X-Spec ``xsect`` command [1]_.

    Parameters
    ----------
    name : { 'bcmc', 'obcm', 'vern' }
       The options are: 'bcmc' from [2]_ with a new
       He cross-section based on [3]_; 'obcm' which is,
       the same as 'bcmc', but with the He cross-section
       from [2]_, or 'vern' [4]_.

    See Also
    --------
    get_xsversion, get_xsxsect, set_xschatter

    Notes
    -----
    The screen output of this function is controlled by the
    X-Spec chatter setting (`set_xschatter`).

    References
    ----------
    .. [1] http://heasarc.nasa.gov/xanadu/xspec/manual/XSxsect.html
           Note that this may refer to a newer version than the
           compiled version used by Sherpa; use `get_xsversion` to
           check.

    .. [2] Balucinska-Church & McCammon (1992; Ap.J.400, 699).
           http://adsabs.harvard.edu/abs/1992ApJ...400..699B

    .. [3] Yan, M., Sadeghpour, H. R., Dalgarno, A. 1998,
           Ap.J. 496, 1044
           http://adsabs.harvard.edu/abs/1998ApJ...496.1044Y

    .. [4] Verner et. al., 1996, Ap.J., 465, 487.
           http://adsabs.harvard.edu/abs/1996ApJ...465..487V

    Examples
    --------

    >>> set_xsxsect('vern')
     Cross Section Table set to vern:  Verner, Ferland, Korista, and Yakovlev 1996
    """

    _xspec.set_xsxsect(name)


# Wrap the XSET function in Python, so that we can keep a record of
# the strings the user sent as specific XSPEC model strings (if any) during
# the session.  Only store if setting was successful.
# See:
# http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSxset.html
modelstrings = {}

# Store any path changes
xspecpaths = {}


def get_xsxset(name):
    """Return the X-Spec model setting.

    Parameters
    ----------
    name : str
       The name of the setting (converted to upper case before being
       sent to X-Spec). There is no check that the name is valid.

    Returns
    -------
    val : str
       Returns the value set by a previous call to `set_xsxset` or the
       empty string, if the value has not been previously set.

    See Also
    --------
    set_xsxset

    Notes
    -----
    Due to the way X-Spec model settings work, `get_xsxset` will only
    return a value if it has previously been set with a call to
    `set_xsxset`. There is no way to retrive the default value of a
    setting.

    Examples
    --------

    >>> set_xsxset("POW_EMIN", "0.5")
    >>> get_xsxset("pow_emin")
    '0.5'

    """
    name = name.upper()
    return _xspec.get_xsxset(name)


def set_xsxset(name, value):
    """Set a X-Spec model setting.

    Set variables used by X-Spec models. It is equivalent to the
    X-Spec `xset` command [1]_, but only for setting the model
    database settings. See `set_xsabund`, `set_xscosmo`, and
    `set_xsxsect` for the other settings.

    Parameters
    ----------
    name : str
       The name of the setting. It is converted to upper case before
       being used. There is no check that the name is valid.
    value : str
       The new value of the setting. It must be given as a string.

    See Also
    --------
    get_xsxset, get_xsversion, set_xsabund, set_xschatter, set_xscosmo,
    set_xsxsect

    Notes
    -----
    The available settings are listed at [1]_. Not all the X-Spec
    model types are supported by Sherpa - for instance X-Spec "mixing
    models" - so changing some of these settings will make no
    difference. The X-Spec chatter setting can be increased with
    `set_xschatter` if it is not clear if a setting is being used.

    The model settings are stored so that they can be included in the
    output of `sherpa.astro.ui.utils.save_all`.

    References
    ----------

    .. [1] http://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html
           Note that this may refer to a newer version than the
           compiled version used by Sherpa; use `get_xsversion` to
           check.

    Examples
    --------

    >>> set_xsxset('NEIVERS', '2.0')

    >>> set_xsxset('NEIAPECROOT', '/data/spectral/modelData/APEC_nei_v11')

    >>> set_xsxset('POW_EMIN', '0.5')
    >>> set_xsxset('POW_EMAX', '2.0')

    """
    name = name.upper()
    _xspec.set_xsxset(name, value)
    if get_xsxset(name) != "":
        modelstrings[name] = get_xsxset(name)


def get_xspath_manager():
    """Return the path to the files describing the XSPEC models.

    Returns
    -------
    path : str
       The path to the manager directory containing the various model-data
       files used by XSPEC.

    See Also
    --------
    get_xspath_model, set_xspath_manager

    Examples
    --------

    >>> get_xspath_manager()
    '/usr/local/heasoft-6.22/x86_64-unknown-linux-gnu-libc2.24/../spectral/manager'
    """

    return _xspec.get_xspath_manager()


def get_xspath_model():
    """Return the path to the model data files.

    Returns
    -------
    path : str
       The path to the directory containing the files used by
       the XSPEC models.

    See Also
    --------
    get_xspath_manager

    Examples
    --------

    >>> get_xspath_model()
    '/usr/local/heasoft-6.22/x86_64-unknown-linux-gnu-libc2.24/../spectral/modelData'
    """

    return _xspec.get_xspath_model()


def set_xspath_manager(path):
    """Set the path to the files describing the XSPEC models.

    Parameters
    ----------
    path : str
        The new path.

    See Also
    --------
    get_xspath_manager : Return the path to the files describing the XSPEC models.

    Examples
    --------
    >>> set_xspath_manager('/data/xspec/spectral/manager')
    """

    _xspec.set_xspath_manager(path)
    spath = get_xspath_manager()
    if spath != path:
        raise IOError("Unable to set the XSPEC manager path " +
                      "to '{}'".format(path))

    xspecpaths['manager'] = path


# Provide XSPEC module state as a dictionary.  The "cosmo" state is
# a 3-tuple, and "modelstrings" is a dictionary of model strings
# applicable to certain models.  The abund and xsect settings are
# strings.  The chatter setting is an integer.  Please see the
# XSPEC manual concerning the following commands: abund, chatter,
# cosmo, xsect, and xset.
# http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/Control.html
# http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/Setting.html
#
# The path dictionary contains the manager path, which can be
# explicitly set. It could also contain the model path, but there
# is no XSPEC routine to change that; instead a user would set the
# XSPEC_MDATA_DIR environment variable before starting XSPEC.
# Should this path be included?
#
def get_xsstate():
    """Return the state of the XSPEC module.

    Returns
    -------
    state : dict
        The current settings for the XSPEC module, including but not
        limited to: the abundance and cross-section settings, parameters
        for the cosmological model, any XSET parameters that have been
        set, and changes to the paths used by the model library.

    See Also
    --------
    get_xsabund, get_xschatter, get_xscosmo, get_xsxsect, get_xsxset,
    set_xsstate
    """

    # Do not return the internal dictionary but a copy of it.
    return {"abund": get_xsabund(),
            "chatter": get_xschatter(),
            "cosmo": get_xscosmo(),
            "xsect": get_xsxsect(),
            "modelstrings": modelstrings.copy(),
            "paths": xspecpaths.copy()}


def set_xsstate(state):
    """Restore the state of the XSPEC module.

    Parameters
    ----------
    state : dict
        The current settings for the XSPEC module. This is expected to
        match the return value of ``get_xsstate``, and so uses the
        keys: 'abund', 'chatter', 'cosmo', 'xsect', 'modelstrings',
        and 'paths'.

    See Also
    --------
    get_xsstate, set_xsabund, set_xschatter, set_xscosmo, set_xsxsect,
    set_xsxset

    Notes
    -----
    The state of the XSPEC module will only be changed if all
    the required keys in the dictionary are present. All keys apart
    from 'paths' are required.
    """

    if type(state) == dict and \
       'abund' in state and \
       'chatter' in state and \
       'cosmo' in state and \
       'xsect' in state and \
       'modelstrings' in state:

        h0, q0, l0 = state["cosmo"]

        set_xsabund(state["abund"])
        set_xschatter(state["chatter"])
        set_xscosmo(h0, q0, l0)
        set_xsxsect(state["xsect"])
        for name in state["modelstrings"].keys():
            set_xsxset(name, state["modelstrings"][name])

        # This is optional to support re-loading state information
        # from a version of XSPEC which did not provide the path
        # information.
        #
        try:
            managerpath = state['paths']['manager']
        except KeyError:
            managerpath = None

        if managerpath is not None:
            set_xspath_manager(managerpath)


def read_xstable_model(modelname, filename):
    """Create a XSPEC table model.

    XSPEC additive (atable, [1]_) and multiplicative (mtable, [2]_)
    table models are supported.

    Parameters
    ----------
    modelname : str
       The identifier for this model component.
    filename : str
       The name of the FITS file containing the data, which should
       match the XSPEC table model definition [3]_.

    Returns
    -------
    tablemodel : XSTableModel instance

    References
    ----------

    .. [1] http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSmodelAtable.html

    .. [2] http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSmodelMtable.html

    .. [3] http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html

    Examples
    --------

    Load in the XSPEC table model from the file 'bbrefl_1xsolar.fits'
    and create a model labelled 'xmdl', which is then returned:

    >>> mdl = read_xstable_model('xmdl', 'bbrefl_1xsolar.fits')
    >>> print(mdl)

    """

    # TODO: how to avoid loading this if no backend is available
    import sherpa.astro.io
    read_tbl = sherpa.astro.io.backend.get_table_data
    read_hdr = sherpa.astro.io.backend.get_header_data

    blkname = 'PRIMARY'
    hdrkeys = ['HDUCLAS1', 'REDSHIFT', 'ADDMODEL']
    hdr = read_hdr(filename, blockname=blkname, hdrkeys=hdrkeys)

    addmodel = bool_cast(hdr[hdrkeys[2]])
    addredshift = bool_cast(hdr[hdrkeys[1]])

    # TODO: change Exception to something more useful
    if str(hdr[hdrkeys[0]]).upper() != 'XSPEC TABLE MODEL':
        raise Exception("Not an XSPEC table model")

    blkname = 'PARAMETERS'
    colkeys = ['NAME', 'INITIAL', 'DELTA', 'BOTTOM', 'TOP',
               'MINIMUM', 'MAXIMUM']
    hdrkeys = ['NINTPARM', 'NADDPARM']

    (colnames, cols,
     name, hdr) = read_tbl(filename, colkeys=colkeys, hdrkeys=hdrkeys,
                           blockname=blkname, fix_type=False)
    nint = int(hdr[hdrkeys[0]])
    return XSTableModel(filename, modelname, *cols,
                        nint=nint, addmodel=addmodel,
                        addredshift=addredshift)


# The model classes are added to __all__ at the end of the file
#
# Note that not all routines from _xspec are re-exported here.
#
__all__ = ('get_xschatter', 'get_xsabund', 'get_xscosmo', 'get_xsxsect',
           'set_xschatter', 'set_xsabund', 'set_xscosmo', 'set_xsxsect',
           'get_xsversion', 'set_xsxset', 'get_xsxset', 'set_xsstate',
           'get_xsstate')


def _f77_or_c_12100(name):
    """Models which changed from Fortran to C linkage in 12.10.0"""

    return "C_" + name if equal_or_greater_than("12.10.0") else name


class XSModel(RegriddableModel1D, metaclass=ModelMeta):
    """The base class for XSPEC models.

    It is expected that sub-classes are used to represent the
    five different types of XSPEC model (additive, multiplicative,
    convolution, pile up, mixing, and tables), although not all are
    currently supported in Sherpa.

    Notes
    -----
    The XSPEC models are evaluated on a one-dimensional, integrated,
    contiguous grid. When the ``calc`` method is called with both
    low and high bin values, the arrays are converted into a single
    array matching the XSPEC calling convention - that is elo_0,
    elo_1, ..., elo_n for n bins (so the last value is the upper edge
    of the last bin) - adding in any bins to account for a non-contiguous
    input. This array is used to evaluate the model, and then the
    return value is created by removing any extra bins that had to be
    added to account for non-contiguous input values.

    If used on an unbinned dataset, so only one array is sent to
    ``calc``, then the input values are taken to match the XSPEC
    calling convention - i.e. a contiguous grid where the last element
    represents the upper edge of the last bin. This means that for
    an input grid of ``n`` points, the returned array will contain
    ``n`` values, but the last element will be zero.
    """

    version_enabled = True

    @modelCacher1d
    def calc(self, *args, **kwargs):
        # Ensure output is finite (Keith Arnaud mentioned that XSPEC
        # does this as a check). This is done at this level (Python)
        # rather than in the C++ interface since:
        #  - it is easier
        #  - it allows for the user to find out what bins are bad,
        #    by directly calling the _calc function of a model
        out = self._calc(*args, **kwargs)

        # This check is being skipped in the 4.8.0 release as it
        # has had un-intended consequences. It should be re-evaluated
        # once we have had more experience with the issue.
        #
        # if not numpy.isfinite(out).all():
        #     # TODO: Should this be using a "Sherpa error class"?
        #     # I am not convinced that FloatingPointError is the best
        #     # type.
        #     msg = "model {} has created NaN or +/-Inf value".format(
        #           self.name)
        #     raise FloatingPointError(msg)

        return out


class XSTableModel(XSModel):
    """Interface to XSPEC table models.

    XSPEC supports loading in user-supplied data files for use
    as a table model [1]_. This class provides a low-level
    way to access this functionality. A simpler interface is provided
    by ``read_xstable_model`` and ``sherpa.astro.ui.load_xstable_model``.

    Parameters
    ----------
    filename : str
        The name of the FITS file containing the data for the XSPEC
        table model; the format is described in [2]_.
    name : str
        The name to use for the instance of the table model.
    parnames : sequence
        The parameter names. This corresponds to the "NAME" column from
        the "PARAMETER" block of the input file. Any invalid characters
        in each name will be replaced by the '_' character.
    initvals : sequence
        The initial values for each parameter. This corresponds to the
        "INITIAL" column from the "PARAMETER" block of the input file.
    delta : sequence
        The delta value for each parameter. This corresponds to the
        "DELTA" column from the "PARAMETER" block of the input file.
    mins, maxes, hardmins, hardmaxes : sequence
        The valid range of each parameter. These correspond to the
        "BOTTOM", "TOP", "MINIMUM", and "MAXIMUM" columns from the
        "PARAMETER" block of the input file.
    nint : int
        The first ``nint`` parameters are marked as thawed by default,
        the remaining default to frozen.
    addmodel : bool
        Is this an additive model (``True``) or multiplicative model
        (``False``)? It should be set to the value of the "ADDMODEL"
        keyword of the primary header of the input file.
    addredshift : bool
        If ``True`` then a redshift parameter is added to the parameters.
        It should be set to the value of the "REDSHIFT" keyword of the
        primary header of the input file.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html
    .. [2] http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html

    """

    def __init__(self, filename, name='xstbl', parnames=(),
                 initvals=(), delta=(), mins=(), maxes=(), hardmins=(),
                 hardmaxes=(), nint=0, addmodel=False, addredshift=False):

        # make translation table to turn reserved characters into '_'
        bad = string.punctuation + string.whitespace
        tbl = str.maketrans(bad, '_' * len(bad))

        pars = []
        for ii in range(len(parnames)):
            isfrozen = True
            if nint > 0:
                isfrozen = False

            # We don't explicitly demand this is a string, so support
            # byte strings (given that the values are read from the
            # FITS file this is a possibility).
            #
            try:
                parname = str(parnames[ii], 'utf-8')
            except TypeError:
                parname = parnames[ii]

            parname = parname.strip().lower().translate(tbl)
            par = Parameter(name, parname, initvals[ii],
                            mins[ii], maxes[ii],
                            hardmins[ii], hardmaxes[ii], frozen=isfrozen)
            self.__dict__[parname] = par
            pars.append(par)
            nint -= 1

        self.filename = filename
        self.addmodel = addmodel

        if addredshift:
            self.redshift = Parameter(name, 'redshift', 0., 0., 5.,
                                      0.0, hugeval, frozen=True)
            pars.append(self.redshift)

        if addmodel:
            self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0,
                                  hugeval)
            pars.append(self.norm)

        XSModel.__init__(self, name, pars)

    def fold(*args, **kwargs):
        pass

    @modelCacher1d
    def calc(self, p, *args, **kwargs):
        # The function used depends on XSPEC version and, prior
        # to XSPEC 12.10.1, the type of table.
        #
        # Note that this is lacking support for "exp" models.
        # It should also not be a run-time decision, since the
        # logic could be in __init__, but that can be changed
        # at a later date.
        #
        for param, value in zip(self.pars, p):
            if value < param._hard_min:
                raise ParameterErr('edge', self.name, 'minimum', param._hard_min)
            if value > param._hard_max:
                raise ParameterErr('edge', self.name, 'maximum', param._hard_max)

        if hasattr(_xspec, 'tabint'):
            tabtype = 'add' if self.addmodel else 'mul'
            return _xspec.tabint(p,
                                 filename=self.filename, tabtype=tabtype,
                                 *args, **kwargs)

        if self.addmodel:
            func = _xspec.xsatbl
        else:
            func = _xspec.xsmtbl

        return func(p, filename=self.filename, *args, **kwargs)


class XSAdditiveModel(XSModel):
    """The base class for XSPEC additive models.

    The XSPEC additive models are listed at [1]_.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/Additive.html

    """

    def guess(self, dep, *args, **kwargs):
        if hasattr(self, 'norm'):
            norm = guess_amplitude(dep, *args)
            param_apply_limits(norm, self.norm, **kwargs)


class XSMultiplicativeModel(XSModel):
    """The base class for XSPEC multiplicative models.

    The XSPEC multiplicative models are listed at [1]_.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/Multiplicative.html

    """

    pass


@version_at_least("12.10.1")
class XSagnsed(XSAdditiveModel):
    """The XSPEC agnsed model: AGN SED model

    The model is described at [1]_.

    Attributes
    ----------
    mass
        The black hole mass, in solar units.
    dist
        The comoving (proper) distance, in Mpc.
    logmdot
        log of mdot, where mdot = Mdot/Mdot_Edd and
        eta Mdot_Edd c^2 = L_Edd
    astar
        The black hole spin (dimensionless)
    cosi
        The cosine of the inclination angle i for the warm
        Comptonising component and the outer disc.
    kTe_hot
        The electron temperature for the hot Comptonisation component
        in keV. If negative then only the hot Comptonisation component
        is used.
    kTe_warm
        The electron temperature for the warm Comptonisation component
        in keV. If negative then only the warm Comptonisation component
        is used.
    Gamma_hot
        The spectral index of the hot Comptonisation component. If
        negative, the code will use the value calculated via eq.(2) of
        KD18 (see [1]_).
    Gamma_warm
        The spectral index of the warm Comptonisation component. If
        negative then only the outer disc component is used.
    R_hot
        The outer radius of the hot Comptonisation component in Rg.
    R_warm
        The outer radius of the warm Comptonisation component in Rg.
    logrout
        The log of the outer radius of the disc in units of Rg. If
        negative, the code will use the self gravity radius as calculated
        from Laor & Netzer (1989) (see [1]_).
    Htmax
        The upper limit of the scale height for the hot Comptonisation
        component in Rg. If smaller than R_hot, the hot Comptonisation
        region is a sphere of radius Htmax by keeping Ldiss_hot determined
        by R_hot via eq.(2) of KD18 (see [1]_).
    reprocess
        If this parameter is 0, reprocessing is not considered. If it
        is 1, reprocessing is included.
    redshift
        The redshift.
    norm
        The normalization of the model.

    See Also
    --------
    XSqsosed

    Notes
    -----
    This model is only available when used with XSPEC 12.10.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelAgnsed.html

    """

    __function__ = "agnsed"

    def __init__(self, name='agnsed'):
        self.mass = Parameter(name, 'mass', 1e7, 1.0, 1e10, 1.0, 1e10,
                              'solar', frozen=True)
        self.dist = Parameter(name, 'dist', 100, 0.01, 1e9, 0.01, 1e9,
                              'Mpc', frozen=True)
        self.logmdot = Parameter(name, 'logmdot', -1, -10, 2, -10, 2)
        self.astar = Parameter(name, 'astar', 0.0, -1, 0.998, -1, 0.998,
                               frozen=True)
        self.cosi = Parameter(name, 'cosi', 0.5, 0.05, 1.0, 0.05, 1.0,
                              frozen=True)
        # TODO: allow negative values
        self.kTe_hot = Parameter(name, 'kTe_hot', 100.0, 10, 300, 10, 300,
                                 'keV(_pl)', frozen=True)
        self.kTe_warm = Parameter(name, 'kTe_warm', 0.2, 0.1, 0.5, 0.1, 0.5,
                                  'keV(_sc)')
        self.Gamma_hot = Parameter(name, 'Gamma_hot', 1.7, 1.3, 3, 1.3, 3,
                                   '(_calc)')
        self.Gamma_warm = Parameter(name, 'Gamma_warm', 2.7, 2, 5, 2, 10,
                                    '(_disk)')

        self.R_hot = Parameter(name, 'R_hot', 10, 6, 500, 6, 500, 'Rg')
        self.R_warm = Parameter(name, 'R_warm', 20, 6, 500, 6, 500, 'Rg')

        self.logrout = Parameter(name, 'logrout', -1, -3, 7, -3, 7,
                                 '(_selfg)', frozen=True)

        self.Htmax = Parameter(name, 'Htmax', 100, 6, 200, 6, 200,
                               'Rg', frozen=True)

        self.reprocess = Parameter(name, 'reprocess', 1, 0, 1, 0, 1,
                                   '0off/1on', alwaysfrozen=True)

        self.redshift = Parameter(name, 'redshift', 0, 0, 1, 0, 1,
                                  frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.mass, self.dist, self.logmdot, self.astar, self.cosi,
                self.kTe_hot, self.kTe_warm, self.Gamma_hot, self.Gamma_warm,
                self.R_hot, self.R_warm, self.logrout, self.Htmax,
                self.reprocess, self.redshift, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


class XSapec(XSAdditiveModel):
    """The XSPEC apec model: APEC emission spectrum.

    The model is described at [1]_. The ``set_xsabund`` and ``get_xsabund``
    functions change and return the current settings for the relative
    abundances of the metals. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keywords "APECROOT", "APECTHERMAL", "APECVELOCITY",
    and "APEC_TRACE_ABUND".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function and the "APEC_TRACE_ABUND" xset
        keyword.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbapec, XSbvapec, XSbvvapec, XSnlapec, XSsnapec, XSvapec, XSvvapec

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelApec.html

    """
    __function__ = "C_apec" if equal_or_greater_than("12.9.1") else "xsaped"

    def __init__(self, name='apec'):
        self.kT = Parameter(name, 'kT', 1., 0.008, 64.0, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.redshift, self.norm))


class XSnlapec(XSAdditiveModel):
    """The XSPEC nlapec model: continuum-only APEC emission spectrum.

    The model is described at [1]_. The ``set_xsabund`` and ``get_xsabund``
    functions change and return the current settings for the relative
    abundances of the metals. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keywords "APECROOT", "APECTHERMAL", "APECVELOCITY",
    and "APEC_TRACE_ABUND".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function and the "APEC_TRACE_ABUND" xset
        keyword.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSapec, XSbapec, XSbvapec, XSbvvapec, XSvapec, XSvvapec

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNlapec.html

    """

    __function__ = "C_nlapec"

    def __init__(self, name='nlapec'):
        self.kT = Parameter(name, 'kT', 1., 0.008, 64.0, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.redshift, self.norm))


class XSbapec(XSAdditiveModel):
    """The XSPEC bapec model: velocity broadened APEC thermal plasma model.

    The model is described at [1]_. The ``set_xsabund`` and ``get_xsabund``
    functions change and return the current settings for the relative
    abundances of the metals. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keywords "APECROOT" and "APEC_TRACE_ABUND".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function and the "APEC_TRACE_ABUND" xset
        keyword.
    Redshift
        The redshift of the plasma.
    Velocity
        The gaussian sigma of the velocity broadening, in km/s.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSapec, XSbtapec, XSbvapec, XSbvvapec, XSvapec, XSvvapec

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBapec.html

    """

    __function__ = "C_bapec" if equal_or_greater_than("12.9.1") else "xsbape"

    def __init__(self, name='bapec'):
        self.kT = Parameter(name, 'kT', 1., 0.008, 64.0, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.Velocity = Parameter(name, 'Velocity', 0., 0., 1.e6, 0.0, hugeval, 'km/s', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.Redshift, self.Velocity, self.norm))


class XSbbody(XSAdditiveModel):
    """The XSPEC bbody model: blackbody spectrum.

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The temperature of the object, in keV.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbbodyrad, XSzbbody

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBbody.html

    """

    __function__ = "xsblbd"

    def __init__(self, name='bbody'):
        self.kT = Parameter(name, 'kT', 3.0, 1.e-2, 100., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval, 'L39 / (D10)**2')
        XSAdditiveModel.__init__(self, name, (self.kT, self.norm))


class XSbbodyrad(XSAdditiveModel):
    """The XSPEC bbodyrad model: blackbody spectrum, area normalized.

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The temperature of the object, in keV.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbbody, XSzbbody

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBbodyrad.html

    """

    __function__ = "xsbbrd"

    def __init__(self, name='bbodyrad'):
        self.kT = Parameter(name, 'kT', 3., 1e-3, 100, 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.norm))


# DOC-NOTE: the XSPEC documentation has a different parameter order to
#           the code.
class XSbexrav(XSAdditiveModel):
    """The XSPEC bexrav model: reflected e-folded broken power law, neutral medium.

    The model is described at [1]_.

    Attributes
    ----------
    Gamma1
        The power-law index of the first power-law component.
    breakE
        The break energy, in keV.
    Gamma2
        The power-law index of the second power-law component.
    foldE
        The e-folding energy (Ec) in keV. If zero there is no cut off.
    rel_refl
        The reflection scaling parameter (a value of 1 for an
        isotropic source above the disk).
    cosIncl
        The cosine of the inclination angle.
    abund
        The abundance of the elements heaver than He relative to their
        solar abundance, as set by the ``set_xsabund`` function.
    Fe_abund
        The iron abundance relative to the solar abundance, as set by
        the ``set_xsabund`` function.
    redshift
        The redshift of the source.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbexriv

    Notes
    -----
    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the BEXRAV_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBexrav.html

    """

    __function__ = "C_xsbexrav"

    def __init__(self, name='bexrav'):
        self.Gamma1 = Parameter(name, 'Gamma1', 2., -9., 9., -hugeval, hugeval)
        self.breakE = Parameter(name, 'breakE', 10., 0.1, 1000., 0.0, hugeval, 'keV')
        self.Gamma2 = Parameter(name, 'Gamma2', 2., -9., 9., -hugeval, hugeval)
        self.foldE = Parameter(name, 'foldE', 100., 1., 1.e6, 0.0, hugeval, 'keV')
        self.rel_refl = Parameter(name, 'rel_refl', 0., 0., 10., 0.0, hugeval)
        self.cosIncl = Parameter(name, 'cosIncl', 0.45, 0.05, 0.95, 0.0, hugeval, frozen=True)
        self.abund = Parameter(name, 'abund', 1., 0.0, 1.e6, 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.0, 1.e6, 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Gamma1, self.breakE, self.Gamma2, self.foldE, self.rel_refl, self.cosIncl, self.abund, self.Fe_abund, self.redshift, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.breakE, **kwargs)


class XSbexriv(XSAdditiveModel):
    """The XSPEC bexriv model: reflected e-folded broken power law, ionized medium.

    The model is described at [1]_.

    Attributes
    ----------
    Gamma1
        The power-law index of the first power-law component.
    breakE
        The break energy, in keV.
    Gamma2
        The power-law index of the second power-law component.
    foldE
        The e-folding energy (Ec) in keV. If zero there is no cut off.
    rel_refl
        The reflection scaling parameter (a value of 1 for an
        isotropic source above the disk).
    redshift
        The redshift of the source.
    abund
        The abundance of the elements heaver than He relative to their
        solar abundance, as set by the ``set_xsabund`` function.
    Fe_abund
        The iron abundance relative to the solar abundance, as set by
        the ``set_xsabund`` function.
    cosIncl
        The cosine of the inclination angle.
    T_disk
        The disk temperature in K.
    xi
        The disk ionization parameter: see [1]_ for an explanation.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbexrav

    Notes
    -----
    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the BEXRIV_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBexriv.html

    """

    __function__ = "C_xsbexriv"

    def __init__(self, name='bexriv'):
        self.Gamma1 = Parameter(name, 'Gamma1', 2., -9., 9., -hugeval, hugeval)
        self.breakE = Parameter(name, 'breakE', 10., 0.1, 1000., 0.0, hugeval, 'keV')
        self.Gamma2 = Parameter(name, 'Gamma2', 2., -9., 9., -hugeval, hugeval)
        self.foldE = Parameter(name, 'foldE', 100., 1., 1.e6, 0.0, hugeval, 'keV')
        self.rel_refl = Parameter(name, 'rel_refl', 0., 0., 1.e6, 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.abund = Parameter(name, 'abund', 1., 0.0, 1.e6, 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.0, 1.e6, 0.0, hugeval, frozen=True)
        self.cosIncl = Parameter(name, 'cosIncl', 0.45, 0.05, 0.95, 0.0, hugeval, frozen=True)
        self.T_disk = Parameter(name, 'T_disk', 3.e4, 1.e4, 1.e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 1., 0., 1.e3, 0.0, hugeval, 'erg cm/s')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Gamma1, self.breakE, self.Gamma2, self.foldE, self.rel_refl, self.redshift, self.abund, self.Fe_abund, self.cosIncl, self.T_disk, self.xi, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.breakE, **kwargs)


class XSbknpower(XSAdditiveModel):
    """The XSPEC bknpower model: broken power law.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndx1
        The power law photon index for energies less than BreakE.
    BreakE
        The break energy, in keV.
    PhoIndx2
        The power law photon index for energies greater than BreakE.
    norm
        The normalization of the model. See [1]_ for details, as its
        meaning depends on whether the "POW_EMIN" or "POW_EMAX"
        keywords have been set with ``set_xsxset``.

    See Also
    --------
    XSbkn2pow, XScutoffpl, XSpowerlaw, XSzbknpower, XSzcutoffpl, XSzpowerlw

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBknpower.html

    """

    __function__ = "C_brokenPowerLaw"

    def __init__(self, name='bknpower'):
        self.PhoIndx1 = Parameter(name, 'PhoIndx1', 1., -2., 9., -hugeval, hugeval)
        self.BreakE = Parameter(name, 'BreakE', 5., 1.e-2, 1.e6, 0.0, hugeval, 'keV')
        self.PhoIndx2 = Parameter(name, 'PhoIndx2', 2., -2., 9., -hugeval, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.PhoIndx1, self.BreakE, self.PhoIndx2, self.norm)

        XSAdditiveModel.__init__(self, name, pars)

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.BreakE, **kwargs)


class XSbkn2pow(XSAdditiveModel):
    """The XSPEC bkn2pow model: broken power law with two breaks.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndx1
        The power law photon index for energies less than BreakE1.
    BreakE1
        The first break energy, in keV.
    PhoIndx2
        The power law photon index for energies greater than BreakE1
        but less than BreakE2.
    BreakE2
        The second break energy, in keV.
    PhoIndx3
        The power law photon index for energies greater than BreakE2.
    norm
        The normalization of the model. See [1]_ for details, as its
        meaning depends on whether the "POW_EMIN" or "POW_EMAX"
        keywords have been set with ``set_xsxset``.

    See Also
    --------
    XSbknpower, XScutoffpl, XSpowerlaw, XSzcutoffpl, XSzpowerlw

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBkn2pow.html

    """

    __function__ = "C_broken2PowerLaw"

    def __init__(self, name='bkn2pow'):
        self.PhoIndx1 = Parameter(name, 'PhoIndx1', 1., -2., 9., -hugeval, hugeval)
        self.BreakE1 = Parameter(name, 'BreakE1', 5., 1.e-2, 1.e6, 0.0, hugeval, 'keV')
        self.PhoIndx2 = Parameter(name, 'PhoIndx2', 2., -2., 9., -hugeval, hugeval)
        self.BreakE2 = Parameter(name, 'BreakE2', 10., 1.e-2, 1.e6, 0.0, hugeval, 'keV')
        self.PhoIndx3 = Parameter(name, 'PhoIndx3', 3., -2., 9., -hugeval, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndx1, self.BreakE1, self.PhoIndx2, self.BreakE2, self.PhoIndx3, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.BreakE1, **kwargs)
        param_apply_limits(pos, self.BreakE2, **kwargs)


class XSbmc(XSAdditiveModel):
    """The XSPEC bmc model: Comptonization by relativistic matter.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``logA`` parameter has been renamed ``log_A`` to match the
       XSPEC definition. The name ``logA`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    kT
        The temperature of the thermal photon source in keV.
    alpha
        The energy spectral index.
    log_A
        The log of the A parameter: see [1]_ for more details.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBmc.html

    """

    __function__ = "xsbmc"

    def __init__(self, name='bmc'):
        self.kT = Parameter(name, 'kT', 1., 1.e-2, 100., 0.0, hugeval, 'keV')
        self.alpha = Parameter(name, 'alpha', 1., 1.e-2, 4.0, 0.0, hugeval)
        self.log_A = Parameter(name, 'log_A', 0.0, -6.0, 6.0, -hugeval, hugeval, aliases=["logA"])
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.alpha, self.log_A, self.norm))


class XSbremss(XSAdditiveModel):
    """The XSPEC bremss model: thermal bremsstrahlung.

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The plasma temperature in keV.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSvbremss, XSzbremss

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBremss.html

    """

    __function__ = "xsbrms"

    def __init__(self, name='bremss'):
        self.kT = Parameter(name, 'kT', 7.0, 1.e-4, 100., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.norm))


@version_at_least("12.10.0")
class XSbrnei(XSAdditiveModel):
    """The XSPEC brnei model: velocity-broadened non-equilibrium recombining collisional plasma.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    kT_init
        The initial temperature of the plasma, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    Tau
        The ionization timescale in units of s/cm^3.
    Redshift
        The redshift of the plasma.
    Velocity
        Velocity broadening in km/s.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbvrnei, XSbvvrnei, XSnei, XSgnei, XSrnei, XSvrnei, XSvvrnei

    Notes
    -----
    This model is only available when used with XSPEC 12.10.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBrnei.html

    """

    __function__ = "C_brnei"

    def __init__(self, name='brnei'):
        self.kT = Parameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = Parameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.Redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.Velocity = Parameter(name, 'Velocity', 0.0, 0.0, 1.0e6, 0.0, 1.0e6,
                                  units='km/s', frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.Abundanc, self.Tau, self.Redshift, self.Velocity, self.norm))


class XSbvapec(XSAdditiveModel):
    """The XSPEC bvapec model: velocity broadened APEC thermal plasma model.

    The model is described at [1]_, with ``XSbapec`` describing how
    it is implemented in Sherpa.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    He, C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element in solar units.
    Redshift
        The redshift of the plasma.
    Velocity
        The gaussian sigma of the velocity broadening, in km/s.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSapec, XSbapec, XSbvvapec, XSvapec, XSvvapec

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBapec.html

    """

    __function__ = "C_bvapec" if equal_or_greater_than("12.9.1") else "xsbvpe"

    def __init__(self, name='bvapec'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0, hugeval, 'keV')
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.Velocity = Parameter(name, 'Velocity', 0., 0., 1.e6, 0.0, hugeval, 'km/s', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Redshift, self.Velocity, self.norm))


@version_at_least("12.10.0")
class XSbvrnei(XSAdditiveModel):
    """The XSPEC bvrnei model: velocity-broadened non-equilibrium recombining collisional plasma.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    kT_init
        The initial temperature of the plasma, in keV.
    H, He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element in solar units.
    Tau
        The ionization timescale in units of s/cm^3.
    Redshift
        The redshift of the plasma.
    Velocity
        The gaussian sigma of the velocity broadening, in km/s.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbrnei, XSbvvrnei, XSnei, XSgnei, XSrnei, XSvrnei, XSvvrnei

    Notes
    -----
    This model is only available when used with XSPEC 12.10.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBrnei.html

    """

    __function__ = "C_bvrnei"

    def __init__(self, name='bvrnei'):
        self.kT = Parameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = Parameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = Parameter(name, 'H', 1., 0., 1., 0.0, 1.0, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -0.999, 10., frozen=True)
        self.Velocity = Parameter(name, 'Velocity', 0., 0., 1.e6, 0.0, 1.e6, 'km/s', frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.Redshift, self.Velocity, self.norm))


class XSc6mekl(XSAdditiveModel):
    """The XSPEC c6mekl model: differential emission measure using Chebyshev representations with multi-temperature mekal.

    The model is described at [1]_.

    Attributes
    ----------
    CPcoef1, CPcoef2, CPcoef3, CPcoef4, CPcoef5, CPcoef6
        Chebyshev polynomial coefficients.
    nH
        H density, in cm^-3.
    abundanc
        The abundance relative to Solar, as set by ``set_xsabund``.
    redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model.

    See Also
    --------
    XSc6pmekl, XSc6pvmkl, XSc6vmekl

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelC6mekl.html

    """

    __function__ = _f77_or_c_12100("c6mekl")

    def __init__(self, name='c6mekl'):
        self.CPcoef1 = Parameter(name, 'CPcoef1', 1.0, -1, 1, -hugeval, hugeval)
        self.CPcoef2 = Parameter(name, 'CPcoef2', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef3 = Parameter(name, 'CPcoef3', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef4 = Parameter(name, 'CPcoef4', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef5 = Parameter(name, 'CPcoef5', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef6 = Parameter(name, 'CPcoef6', 0.5, -1, 1, -hugeval, hugeval)
        self.nH = Parameter(name, 'nH', 1.0, 1.e-5, 1.e19, 0.0, hugeval, 'cm^-3', True)
        self.abundanc = Parameter(name, 'abundanc', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5,
                                              self.CPcoef6, self.nH, self.abundanc, self.redshift, self.switch,
                                              self.norm))


class XSc6pmekl(XSAdditiveModel):
    """The XSPEC c6pmekl model: differential emission measure using Chebyshev representations with multi-temperature mekal.

    The model is described at [1]_. It differs from ``XSc6mekl`` by
    by using the exponential of the 6th order Chebyshev polynomial.

    Attributes
    ----------
    CPcoef1, CPcoef2, CPcoef3, CPcoef4, CPcoef5, CPcoef6
        Chebyshev polynomial coefficients.
    nH
        H density, in cm^-3.
    abundanc
        The abundance relative to Solar, as set by ``set_xsabund``.
    redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model.

    See Also
    --------
    XSc6mekl, XSc6pvmkl, XSc6vmekl

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelC6mekl.html

    """

    __function__ = _f77_or_c_12100("c6pmekl")

    def __init__(self, name='c6pmekl'):
        self.CPcoef1 = Parameter(name, 'CPcoef1', 1.0, -1, 1, -hugeval, hugeval)
        self.CPcoef2 = Parameter(name, 'CPcoef2', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef3 = Parameter(name, 'CPcoef3', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef4 = Parameter(name, 'CPcoef4', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef5 = Parameter(name, 'CPcoef5', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef6 = Parameter(name, 'CPcoef6', 0.5, -1, 1, -hugeval, hugeval)
        self.nH = Parameter(name, 'nH', 1.0, 1.e-5, 1.e19, 0.0, hugeval, 'cm^-3', True)
        self.abundanc = Parameter(name, 'abundanc', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5,
                                              self.CPcoef6, self.nH, self.abundanc, self.redshift, self.switch,
                                              self.norm))


class XSc6pvmkl(XSAdditiveModel):
    """The XSPEC c6pvmkl model: differential emission measure using Chebyshev representations with multi-temperature mekal.

    The model is described at [1]_. It differs from ``XSc6vmekl`` by
    by using the exponential of the 6th order Chebyshev polynomial.

    Attributes
    ----------
    CPcoef1, CPcoef2, CPcoef3, CPcoef4, CPcoef5, CPcoef6
        Chebyshev polynomial coefficients.
    nH
        H density, in cm^-3.
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance relative to Solar.
    redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model.

    See Also
    --------
    XSc6mekl, XSc6pmekl, XSc6vmekl

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelC6mekl.html

    """

    __function__ = _f77_or_c_12100("c6pvmkl")

    def __init__(self, name='c6pvmkl'):
        self.CPcoef1 = Parameter(name, 'CPcoef1', 1.0, -1, 1, -hugeval, hugeval)
        self.CPcoef2 = Parameter(name, 'CPcoef2', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef3 = Parameter(name, 'CPcoef3', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef4 = Parameter(name, 'CPcoef4', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef5 = Parameter(name, 'CPcoef5', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef6 = Parameter(name, 'CPcoef6', 0.5, -1, 1, -hugeval, hugeval)
        self.nH = Parameter(name, 'nH', 1.0, 1.e-5, 1.e19, 0.0, hugeval, 'cm^-3', True)
        self.He = Parameter(name, 'He', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5,
                                              self.CPcoef6, self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na,
                                              self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni,
                                              self.redshift, self.switch, self.norm))


class XSc6vmekl(XSAdditiveModel):
    """The XSPEC c6vmekl model: differential emission measure using Chebyshev representations with multi-temperature mekal.

    The model is described at [1]_.

    Attributes
    ----------
    CPcoef1, CPcoef2, CPcoef3, CPcoef4, CPcoef5, CPcoef6
        Chebyshev polynomial coefficients.
    nH
        H density, in cm^-3.
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance relative to Solar.
    redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model.

    See Also
    --------
    XSc6mekl, XSc6pmekl, XSc6pvmkl

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelC6mekl.html

    """

    __function__ = _f77_or_c_12100("c6vmekl")

    def __init__(self, name='c6vmekl'):
        self.CPcoef1 = Parameter(name, 'CPcoef1', 1.0, -1, 1, -hugeval, hugeval)
        self.CPcoef2 = Parameter(name, 'CPcoef2', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef3 = Parameter(name, 'CPcoef3', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef4 = Parameter(name, 'CPcoef4', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef5 = Parameter(name, 'CPcoef5', 0.5, -1, 1, -hugeval, hugeval)
        self.CPcoef6 = Parameter(name, 'CPcoef6', 0.5, -1, 1, -hugeval, hugeval)
        self.nH = Parameter(name, 'nH', 1.0, 1.e-5, 1.e19, 0.0, hugeval, 'cm^-3', True)
        self.He = Parameter(name, 'He', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5,
                                              self.CPcoef6, self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na,
                                              self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni,
                                              self.redshift, self.switch, self.norm))


class XScemekl(XSAdditiveModel):
    """The XSPEC cemekl model: plasma emission, multi-temperature using mekal.

    The model is described at [1]_.

    Attributes
    ----------
    alpha
        The power-law index of the emissivity function.
    Tmax
        The maxmimum temperature, in keV.
    nH
        H density, in cm^-3.
    abundanc
        The abundance relative to Solar, as set by ``set_xsabund``.
    redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model.

    See Also
    --------
    XScevmkl

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCemekl.html

    """

    __function__ = "cemekl"

    def __init__(self, name='cemekl'):
        self.alpha = Parameter(name, 'alpha', 1.0, 0.01, 10, 0.0, hugeval, frozen=True)
        self.Tmax = Parameter(name, 'Tmax', 1.0, 0.01, 1.e2, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1.0, 1.e-5, 1.e19, 0.0, hugeval, 'cm^-3', True)
        self.abundanc = Parameter(name, 'abundanc', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 1, 0, 1, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.alpha, self.Tmax, self.nH, self.abundanc, self.redshift,
                                              self.switch, self.norm))


class XScevmkl(XSAdditiveModel):
    """The XSPEC cevmkl model: plasma emission, multi-temperature using mekal.

    The model is described at [1]_.

    Attributes
    ----------
    alpha
        The power-law index of the emissivity function.
    Tmax
        The maxmimum temperature, in keV.
    nH
        H density, in cm^-3.
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance relative to Solar.
    redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model.

    See Also
    --------
    XScemekl

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCemekl.html

    """

    __function__ = "C_cemVMekal"

    def __init__(self, name='cevmkl'):
        self.alpha = Parameter(name, 'alpha', 1.0, 0.01, 10, 0.0, hugeval, frozen=True)
        self.Tmax = Parameter(name, 'Tmax', 1.0, 0.01, 1.e2, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1.0, 1.e-5, 1.e19, 0.0, hugeval, 'cm^-3', True)
        self.He = Parameter(name, 'He', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 1, 0, 1, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.alpha, self.Tmax, self.nH, self.He, self.C, self.N, self.O, self.Ne,
                                              self.Na, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe,
                                              self.Ni, self.redshift, self.switch, self.norm))


class XScflow(XSAdditiveModel):
    """The XSPEC cflow model: cooling flow.

    The model is described at [1]_. The results of this model depend
    on the cosmology settings set with ``set_xscosmo``.

    Attributes
    ----------
    slope
        The power-law index of the emissivity function.
    lowT
        The minimum temperature, in keV.
    highT
        The maxmimum temperature, in keV.
    Abundanc
        The abundance relative to Solar, as set by ``set_xsabund``.
    redshift
        The redshift of the plasma.
    norm
        The mass accretion rate (solar mass per year).

    See Also
    --------
    XScevmkl, XSmkcflow

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCflow.html

    """

    __function__ = "C_xscflw"

    def __init__(self, name='cflow'):
        self.slope = Parameter(name, 'slope', 0., -5., 5., -hugeval, hugeval)
        self.lowT = Parameter(name, 'lowT', 0.1, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.highT = Parameter(name, 'highT', 4., 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', .1, 1.e-10, 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.slope, self.lowT, self.highT, self.Abundanc, self.redshift, self.norm))


class XScompbb(XSAdditiveModel):
    """The XSPEC compbb model: Comptonization, black body.

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The blackbody temperature, in keV.
    kTe
        The electron temperature of the hot plasma, in keV.
    tau
        The optical depth of the plasma.
    norm
        The model normalization: it is the same definition as
        used for the ``XSbbodyrad`` model.

    See Also
    --------
    XSbbodyrad

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCompbb.html

    """

    __function__ = "compbb"

    def __init__(self, name='compbb'):
        self.kT = Parameter(name, 'kT', 1.0, 1.e-2, 100., 0.0, hugeval, 'keV')
        self.kTe = Parameter(name, 'kTe', 50, 1., 200., 0.0, hugeval, 'keV', True)
        self.tau = Parameter(name, 'tau', 0.1, 0.0, 10., 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kTe, self.tau, self.norm))


class XScompLS(XSAdditiveModel):
    """The XSPEC compLS model: Comptonization, Lamb & Sanford.

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The blackbody temperature, in keV.
    tau
        The optical depth of the plasma.
    norm
        The model normalization.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCompls.html

    """

    __function__ = "compls"

    def __init__(self, name='compls'):
        self.kT = Parameter(name, 'kT', 2., .01, 10., 0.0, hugeval, 'keV')
        self.tau = Parameter(name, 'tau', 10, .001, 100., 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.tau, self.norm))


# DOC-NOTE: The XSPEC documentation for a number of parameters suggest
#           that they can be negative, but the model.dat limits suggest
#           otherwise (e.g. ktbb and tau_y).
#
class XScompPS(XSAdditiveModel):
    """The XSPEC compPS model: Comptonization, Poutanen & Svenson.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``tauy`` and ``HRcyl`` parameters have been renamed to ``tau_y``
       and ``HovR_cyl`` respectively to match the XSPEC definition. The names
       ``tauy`` and ``HRcyl`` can still be used to access the parameters, but
       they will be removed in a future release.

    Attributes
    ----------
    kTe
        The electron temperature in keV.
    EleIndex
        The electron power-law index.
    Gmin
        The minimum Lorentz factor gamma.
    GMax
        The maximum Lorentz factor gamma: see [1]_ for more details.
    kTbb
        The temperature of the soft photons, in keV.
    tau_y
        The vertical optical depth of the corona: see [1]_ for more
        details.
    geom
        The geometry to use; see [1]_ for more details.
    HovR_cyl
        The value of H/R, when a cylinder geometry is used
        (abs(geom) = 2).
    cosIncl
        The cosine of the inclination angle.
    cov_frac
        The covering fraction of the cold clouds (only used when
        abs(geom) < 4).
    rel_refl
        The amount of reflection (Omega / (2 pi))
    Fe_ab_re
        The iron abundance, in units of solar.
    Me_ab
        The abundance of heavy elements, in units of solar.
    xi
        The disk ionization parameter.
    Tdisk
        The disk temperature for reflection, in K.
    Betor10
        The reflection emissivity law: see [1]_ for more details.
    Rin
        The inner radius of the disk in Schwarzschild units.
    Rout
        The outer radius of the disk in Schwarzschild units.
    redshift
        The redshift of the source.
    norm
        The normalization of the model.

    Notes
    -----
    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the COMPPS_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCompps.html

    """

    __function__ = "C_xscompps"

    def __init__(self, name='compps'):
        self.kTe = Parameter(name, 'kTe', 100., 20., 1.e5, 0.0, hugeval, 'keV')
        self.EleIndex = Parameter(name, 'EleIndex', 2., 0.0, 5., 0.0, hugeval, frozen=True)
        self.Gmin = Parameter(name, 'Gmin', -1., -1., 10., -hugeval, hugeval, frozen=True)
        self.Gmax = Parameter(name, 'Gmax', 1.e3, 10., 1.e4, 0.0, hugeval, frozen=True)
        self.kTbb = Parameter(name, 'kTbb', 0.1, 0.001, 10., 0.0, hugeval, 'keV', True)
        self.tau_y = Parameter(name, 'tau_y', 1.0, 0.05, 3.0, 0.0, hugeval, aliases=["tauy"])
        self.geom = Parameter(name, 'geom', 0.0, -5.0, 4.0, -hugeval, hugeval, frozen=True)
        self.HovR_cyl = Parameter(name, 'HovR_cyl', 1.0, 0.5, 2.0, 0.0, hugeval, frozen=True, aliases=["HRcyl"])
        self.cosIncl = Parameter(name, 'cosIncl', 0.5, 0.05, 0.95, 0.0, hugeval, frozen=True)
        self.cov_frac = Parameter(name, 'cov_frac', 1.0, 0.0, 1.0, 0.0, hugeval, frozen=True)
        self.rel_refl = Parameter(name, 'rel_refl', 0., 0., 1.e4, 0.0, hugeval, frozen=True)
        self.Fe_ab_re = Parameter(name, 'Fe_ab_re', 1., 0.1, 10., 0.0, hugeval, frozen=True)
        self.Me_ab = Parameter(name, 'Me_ab', 1., 0.1, 10., 0.0, hugeval, frozen=True)
        self.xi = Parameter(name, 'xi', 0., 0., 1.e5, 0.0, hugeval, frozen=True)
        self.Tdisk = Parameter(name, 'Tdisk', 1.e6, 1.e4, 1.e6, 0.0, hugeval, 'K', True)
        self.Betor10 = Parameter(name, 'Betor10', -10., -10., 10., -hugeval, hugeval, frozen=True)
        self.Rin = Parameter(name, 'Rin', 10., 6.001, 1.e3, 0.0, hugeval, 'Rs', True)
        self.Rout = Parameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, hugeval, 'Rs', True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kTe, self.EleIndex, self.Gmin, self.Gmax, self.kTbb, self.tau_y, self.geom, self.HovR_cyl, self.cosIncl, self.cov_frac, self.rel_refl, self.Fe_ab_re, self.Me_ab, self.xi, self.Tdisk, self.Betor10, self.Rin, self.Rout, self.redshift, self.norm))


class XScompST(XSAdditiveModel):
    """The XSPEC compST model: Comptonization, Sunyaev & Titarchuk.

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The temperature, in keV.
    tau
        The optical depth of the plasma.
    norm
        The model normalization: see [1]_ for more details.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCompst.html

    """

    __function__ = "compst"

    def __init__(self, name='compst'):
        self.kT = Parameter(name, 'kT', 2., .01, 100., 0.0, hugeval, 'keV')
        self.tau = Parameter(name, 'tau', 10, .001, 100., 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.tau, self.norm))


# DOC-NOTE: the approx parameter is described in XSPEC as being allowed
#           to be negative, but the limits do not permit this.
#
class XScompTT(XSAdditiveModel):
    """The XSPEC compTT model: Comptonization, Titarchuk.

    The model is described at [1]_.

    Attributes
    ----------
    redshift
        The redshift of the source.
    T0
        The input soft photon (Wien) temperature in keV.
    kT
        The plasma temperature, in keV.
    taup
        The plasma optical depth.
    approx
        The geometry setting: a value less than or equal to 1 is a
        disk, above 1 it is a sphere. See [1]_ for more details on
        how this parameter affects the model. It must remain frozen.
    norm
        The normalization of the model.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelComptt.html

    """

    __function__ = "xstitg"

    def __init__(self, name='comptt'):
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.T0 = Parameter(name, 'T0', 0.1, .01, 100., 0.0, hugeval, 'keV')
        self.kT = Parameter(name, 'kT', 50., 2.0, 500., 0.0, hugeval, 'keV')
        self.taup = Parameter(name, 'taup', 1., .01, 100., 0.0, hugeval)
        self.approx = Parameter(name, 'approx', 1.0, 0.0, 5.0, 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.redshift, self.T0, self.kT, self.taup, self.approx, self.norm))


@version_at_least("12.10.1")
class XScph(XSAdditiveModel):
    """The XSPEC cph model: Cooling + heating model for cool core clusters

    The model is described at [1]_.

    Attributes
    ----------
    peakT
        The peak temperature, in keV.
    Abund
        The abundance, relative to Solar.
    Redshift
        The redshift.
    switch
        If 0 calculate, if 1 interpolate, if 2 use AtomDB data.
    norm
        The mass accretion rate, in solar mass per year.

    See Also
    --------
    XSmkcflow, XSvcph

    Notes
    -----
    This model is only available when used with XSPEC 12.10.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCph.html

    """

    __function__ = "C_cph"

    def __init__(self, name='cph'):
        self.peakT = Parameter(name, 'peakT', 2.2, 1e-1, 1e2, 1e-1, 1e2,
                               'keV')
        self.Abund = Parameter(name, 'Abund', 1, 0, 1000, 0, 1000,
                               frozen=True)
        # In XSPEC 12.10.1 the redshift value defaults to 0 but the
        # minimum is 1e-6, so switch to 0.1
        self.Redshift = Parameter(name, 'Redshift', 0.1, 1e-6, 50, 1e-6, 50,
                                  frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2,
                                alwaysfrozen=True)

        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.peakT, self.Abund, self.Redshift, self.switch, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


@version_at_least("12.10.1")
class XSvcph(XSAdditiveModel):
    """The XSPEC vcph model: Cooling + heating model for cool core clusters

    The model is described at [1]_.

    Attributes
    ----------
    peakT
        The peak temperature, in keV.
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element, with respect to Solar.
    Redshift
        The redshift.
    switch
        If 0 calculate, if 1 interpolate, if 2 use AtomDB data.
    norm
        The mass accretion rate, in solar mass per year.

    See Also
    --------
    XSmkcflow, XScph

    Notes
    -----
    This model is only available when used with XSPEC 12.10.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCph.html

    """

    __function__ = "C_vcph"

    def __init__(self, name='vcph'):
        self.peakT = Parameter(name, 'peakT', 2.2, 1e-1, 1e2, 1e-1, 1e2,
                               'keV')

        self.He = Parameter(name, 'He', 1.0, 0., 1000., 0.0, 1000.0,
                            frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 1000., 0.0, 1000.0,
                           frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 1000., 0.0, 1000.0,
                           frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 1000., 0.0, 1000.0,
                           frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 1000., 0.0, 1000.0,
                            frozen=True)
        self.Na = Parameter(name, 'Na', 1.0, 0., 1000., 0.0, 1000.0,
                            frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 1000., 0.0, 1000.0,
                            frozen=True)
        self.Al = Parameter(name, 'Al', 1.0, 0., 1000., 0.0, 1000.0,
                            frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 1000., 0.0, 1000.0,
                            frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 1000., 0.0, 1000.0,
                           frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 1000., 0.0, 1000.0,
                            frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 1000., 0.0, 1000.0,
                            frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 1000., 0.0, 1000.0,
                            frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 1000., 0.0, 1000.0,
                            frozen=True)

        # In XSPEC 12.10.1 the redshift value defaults to 0 but the
        # minimum is 1e-6, so switch to 0.1
        self.Redshift = Parameter(name, 'Redshift', 0.1, 1e-6, 50, 1e-6, 50,
                                  frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2,
                                alwaysfrozen=True)

        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.peakT,
                self.He, self.C, self.N, self.O, self.Ne, self.Na,
                self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca,
                self.Fe, self.Ni,
                self.Redshift, self.switch, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


class XScutoffpl(XSAdditiveModel):
    """The XSPEC cutoffpl model: power law, high energy exponential cutoff.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndex
        The power law photon index.
    HighECut
        The e-folding energy of the exponential rolloff, in keV.
    norm
        The normalization of the model. See [1]_ for details, as its
        meaning depends on whether the "POW_EMIN" or "POW_EMAX"
        keywords have been set with ``set_xsxset``.

    See Also
    --------
    XSbknpower, XSbkn2pow, XSpowerlaw, XSzcutoffpl, XSzpowerlw

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCutoffpl.html

    """

    __function__ = "C_cutoffPowerLaw"

    def __init__(self, name='cutoffpl'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 1., -2., 9., -hugeval, hugeval)
        self.HighECut = Parameter(name, 'HighECut', 15., 1., 500., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.HighECut, self.norm))


class XSdisk(XSAdditiveModel):
    """The XSPEC disk model: accretion disk, black body.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``NSmass`` parameter has been renamed to ``CenMass`` to match
       the XSPEC definition. The name ``NSmass`` can still be used to
       access the parameter, but this name will be removed in a future
       release.

    Attributes
    ----------
    accrate
        The accretion rate, in Eddington luminosities.
    CenMass
        The central mass, in solar mass units.
    Rinn
        The inner disk radius in gravitational units (three
        Schwarzschild radii).
    norm
        The normalization of the model: see [1]_ for details.

    See Also
    --------
    XSdiskbb, XSdiskm, XSdisko

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDisk.html

    """

    __function__ = "disk"

    def __init__(self, name='disk'):
        self.accrate = Parameter(name, 'accrate', 1., 1e-3, 9., 0.0, hugeval)
        self.CenMass = Parameter(name, 'CenMass', 1.4, .4, 10., 0.0, hugeval, units='Msun', frozen=True, aliases=["NSmass"])
        self.Rinn = Parameter(name, 'Rinn', 1.03, 1.01, 1.03, 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.accrate, self.CenMass, self.Rinn, self.norm))


class XSdiskir(XSAdditiveModel):
    """The XSPEC diskir model: Irradiated inner and outer disk.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``LcLd`` parameter has been renamed ``LcovrLd`` to match the
       XSPEC definition. The name ``LcLd`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    kT_disk
        The temperature of the innermost part of the unilluminated
        disk, in keV.
    Gamma
        The asymptotic power-law photon index.
    kT_e
        The electron temperature (high-energy rollover) in keV.
    LcovrLd
        The ratio of the luminosity in the Compton tail to that of
        the unilluminated disk.
    fin
        The fraction of luminosity in the Compton tail which is
        thermalized in the inner disk. Generally fix at 0.1 as
        appropriate for an albedo of 0.3 and solid angle of 0.3.
    rirr
        The radius of the Compton illuminated disk in terms of the
        inner disk radius.
    fout
        The fraction of bolometric flux which is thermalized in the
        outer disk.
    logrout
        The log (base 10) of the outer disk radius in terms of the
        inner disk radius.
    norm
        The model normalization: it is the same definition as
        used for the ``XSdiskbb`` model.

    See Also
    --------
    XSdiskbb

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDiskir.html

    """

    __function__ = "diskir"

    def __init__(self, name='diskir'):
        self.kT_disk = Parameter(name, 'kT_disk', 1.0, 0.01, 5., 0.0, hugeval, 'keV')
        self.Gamma = Parameter(name, 'Gamma', 1.7, 1.001, 5., 0.0, hugeval)
        self.kT_e = Parameter(name, 'kT_e', 100., 5., 1.e3, 0.0, hugeval, 'keV')
        self.LcovrLd = Parameter(name, 'LcovrLd', 0.1, 0., 10., 0.0, hugeval, aliases=["LcLd"])
        self.fin = Parameter(name, 'fin', 1.e-1, 0.0, 1., 0.0, hugeval, frozen=True)
        self.rirr = Parameter(name, 'rirr', 1.2, 1.0001, 10., 1.0001, hugeval)
        self.fout = Parameter(name, 'fout', 1.e-4, 0.0, 1.e-1, 0.0, hugeval)
        self.logrout = Parameter(name, 'logrout', 5.0, 3.0, 7.0, 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kT_disk, self.Gamma, self.kT_e, self.LcovrLd, self.fin, self.rirr, self.fout, self.logrout, self.norm))


class XSdiskbb(XSAdditiveModel):
    """The XSPEC diskbb model: accretion disk, multi-black body components.

    The model is described at [1]_.

    Attributes
    ----------
    Tin
        The temperature at the inner disk radius, in keV.
    norm
        The normalization of the model: see [1]_ for details.

    See Also
    --------
    XSdiskpbb, XSdiskpn

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDiskbb.html

    """

    __function__ = "xsdskb"

    def __init__(self, name='diskbb'):
        self.Tin = Parameter(name, 'Tin', 1., 0., 1000., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Tin, self.norm))


class XSdiskline(XSAdditiveModel):
    """The XSPEC diskline model: accretion disk line emission, relativistic.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``RinM`` and ``RoutM`` parameters have been renamed ``Rin_M``
       and ``Rout_M`` respectively to match the XSPEC definition. The names
       ``RinM`` and ``RoutM`` can still be used to access the parameters,
       but they will be removed in a future release.

    Attributes
    ----------
    LineE
        The line energy in keV.
    Betor10
        The power law dependence of emissivity: see [1]_ for more details.
    Rin_M
        The inner radius, in units of GM^2/c.
    Rout_M
        The outer radius, in units of GM^2/c.
    Incl
        The inclination, in degrees.
    norm
        The model normalization in photon/cm^2/s.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDiskline.html

    """

    __function__ = "C_diskline" if equal_or_greater_than("12.10.1") else "xsdili"

    def __init__(self, name='diskline'):
        self.LineE = Parameter(name, 'LineE', 6.7, 0., 100., 0.0, hugeval, 'keV')
        self.Betor10 = Parameter(name, 'Betor10', -2., -10., 20., -hugeval, hugeval, frozen=True)
        self.Rin_M = Parameter(name, 'Rin_M', 10., 6., 1000., 0.0, hugeval, frozen=True, aliases=["RinM"])
        self.Rout_M = Parameter(name, 'Rout_M', 1000., 0., 1000000., 0.0, hugeval, frozen=True, aliases=["RoutM"])
        self.Incl = Parameter(name, 'Incl', 30., 0., 90., 0.0, hugeval, 'deg')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.LineE, self.Betor10, self.Rin_M, self.Rout_M, self.Incl, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


class XSdiskm(XSAdditiveModel):
    """The XSPEC diskm model: accretion disk with gas pressure viscosity.

    The model is described at [1]_.

    Attributes
    ----------
    NSmass
        The accretion rate, in Eddington luminosities.
    NSmass
        The central mass, in solar mass units.
    Rinn
        The inner disk radius in gravitational units (three
        Schwarzschild radii).
    alpha
        The viscosity.
    norm
        The normalization of the model: see [1]_ for details.

    See Also
    --------
    XSdisk, XSdisko

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDiskm.html

    """

    __function__ = "diskm"

    def __init__(self, name='diskm'):
        self.accrate = Parameter(name, 'accrate', 1., 1e-3, 9., 0.0, hugeval)
        self.NSmass = Parameter(name, 'NSmass', 1.4, .4, 10., 0.0, hugeval, units='Msun', frozen=True)
        self.Rinn = Parameter(name, 'Rinn', 1.03, 1.01, 1.03, 0.0, hugeval, frozen=True)
        self.alpha = Parameter(name, 'alpha', 1., .01, 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.accrate, self.NSmass, self.Rinn, self.alpha, self.norm))


class XSdisko(XSAdditiveModel):
    """The XSPEC disko model: accretion disk, inner, radiation pressure viscosity.

    The model is described at [1]_.

    Attributes
    ----------
    NSmass
        The accretion rate, in Eddington luminosities.
    NSmass
        The central mass, in solar mass units.
    Rinn
        The inner disk radius in gravitational units (three
        Schwarzschild radii).
    alpha
        The viscosity.
    norm
        The normalization of the model: see [1]_ for details.

    See Also
    --------
    XSdisk, XSdiskm

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDisko.html

    """

    __function__ = "disko"

    def __init__(self, name='disko'):
        self.accrate = Parameter(name, 'accrate', 1., 1e-3, 9., 0.0, hugeval)
        self.NSmass = Parameter(name, 'NSmass', 1.4, .4, 10., 0.0, hugeval, units='Msun', frozen=True)
        self.Rinn = Parameter(name, 'Rinn', 1.03, 1.01, 1.03, 0.0, hugeval, frozen=True)
        self.alpha = Parameter(name, 'alpha', 1., .01, 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.accrate, self.NSmass, self.Rinn, self.alpha, self.norm))


class XSdiskpbb(XSAdditiveModel):
    """The XSPEC diskpbb model: accretion disk, power-law dependence for T(r).

    The model is described at [1]_.

    Attributes
    ----------
    Tin
        The temperature at the inner disk radius, in keV.
    p
        The exponent of the radial dependence of the disk temperature.
    norm
        The normalization of the model: see [1]_ for details.

    See Also
    --------
    XSdiskbb, XSdiskpn

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDiskpbb.html

    """

    __function__ = "diskpbb"

    def __init__(self, name='diskpbb'):
        self.Tin = Parameter(name, 'Tin', 1.0, 0.1, 10.0, 0.0, hugeval, 'keV')
        self.p = Parameter(name, 'p', 0.75, 0.5, 1.0, 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Tin, self.p, self.norm))


class XSdiskpn(XSAdditiveModel):
    """The XSPEC diskpn model: accretion disk, black hole, black body.

    The model is described at [1]_.

    Attributes
    ----------
    T_max
        The maximum temperature in the disk, in keV.
    R_in
        The inner disk radius, in units of Rg (GM/c^2).
    norm
        The normalization of the model: see [1]_ for details.

    See Also
    --------
    XSdiskbb, XSdiskpbb

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDiskpn.html

    """

    __function__ = "xsdiskpn"

    def __init__(self, name='diskpn'):
        self.T_max = Parameter(name, 'T_max', 1., 1e-3, 100, 0.0, hugeval, 'keV')
        self.R_in = Parameter(name, 'R_in', 6., 6., 1000., 0.0, hugeval, 'R_g')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.T_max, self.R_in, self.norm))


class XSequil(XSAdditiveModel):
    """The XSPEC equil model: collisional plasma, ionization equilibrium.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSnei, XSgnei, XSpshock, XSvequil

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelEquil.html

    """

    __function__ = "C_equil"

    def __init__(self, name='equil'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.redshift, self.norm))


class XSexpdec(XSAdditiveModel):
    """The XSPEC expdec model: exponential decay.

    The model is described at [1]_.

    Attributes
    ----------
    factor
        The exponential factor.
    norm
        The normalization of the model.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelExpdec.html

    """

    __function__ = "xsxpdec"

    def __init__(self, name='expdec'):
        self.factor = Parameter(name, 'factor', 1.0, 0., 100.0, 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.factor, self.norm))


class XSezdiskbb(XSAdditiveModel):
    """The XSPEC ezdiskbb model: multiple blackbody disk model with zero-torque inner boundary.

    The model is described at [1]_.

    Attributes
    ----------
    T_max
        The maximum temperature in the disk, in keV.
    norm
        The normalization of the model: see [1]_ for more details.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelEzdiskbb.html

    """

    __function__ = "ezdiskbb"

    def __init__(self, name='ezdiskbb'):
        self.T_max = Parameter(name, 'T_max', 1., 0.01, 100., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.T_max, self.norm))


class XSgaussian(XSAdditiveModel):
    """The XSPEC gaussian model: gaussian line profile.

    The model is described at [1]_.

    Attributes
    ----------
    LineE
        The line energy, in keV.
    Sigma
        The line width, in keV. A value of zero means a delta function.
    norm
        The flux in the line, in units of photon/cm^2/s.

    See Also
    --------
    XSagauss, XSlorentz, XSvoigt, XSzagauss, XSzgauss

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGaussian.html

    """

    __function__ = "C_gaussianLine" if equal_or_greater_than("12.9.1") else "xsgaul"

    def __init__(self, name='gaussian'):
        self.LineE = Parameter(name, 'LineE', 6.5, 0., 1.e6, 0.0, hugeval, 'keV')
        self.Sigma = Parameter(name, 'Sigma', 0.1, 0., 10., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Sigma, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


class XSgnei(XSAdditiveModel):
    """The XSPEC gnei model: collisional plasma, non-equilibrium, temperature evolution.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    .. note:: Deprecated in Sherpa 4.10.0

       The ``kT_ave`` parameter has been renamed ``meanKT`` to match the
       XSPEC definition. The name ``kT_ave`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    Tau
        The ionization timescale in units of s/cm^3.
    meankT
        The ionization timescale averaged plasma temperature in keV.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSnei, XSvgnei, XSvvgnei

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGnei.html

    """

    __function__ = "C_gnei"

    def __init__(self, name='gnei'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.meankT = Parameter(name, 'meankT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV', aliases=["kT_ave"])
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.Tau, self.meankT, self.redshift, self.norm))


class XSgrad(XSAdditiveModel):
    """The XSPEC grad model: accretion disk, Schwarzschild black hole.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``TclTef`` parameter has been renamed ``TclovTef`` to match the
       XSPEC definition. The name ``logA`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    D
        The distance to the source in kpc.
    i
        The disk inclination angle, in degrees. A face-on disk has
        i=0.
    Mass
        The mass of the central object, in solar masses.
    Mdot
        The mass accretion rate in units of 10^18 g/s.
    TclovTef
        The spectral hardening factor, Tcol/Teff. See [1]_ for more
        details.
    refflag
        A flag to control the relativistic effects: if positive then
        a relativistic calculation is used; if zero or negative then
        a Newtonian calculation is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model. It should be fixed to 1.

    See Also
    --------
    XSkerbb

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGrad.html

    """

    __function__ = "grad"

    def __init__(self, name='grad'):
        self.D = Parameter(name, 'D', 10.0, 0.0, 10000., 0.0, hugeval, 'kpc', True)
        self.i = Parameter(name, 'i', 0.0, 0.0, 90.0, 0.0, hugeval, 'deg', True)
        self.Mass = Parameter(name, 'Mass', 1.0, 0.0, 100.0, 0.0, hugeval, 'solar')
        self.Mdot = Parameter(name, 'Mdot', 1.0, 0.0, 100.0, 0.0, hugeval, '1e18')
        self.TclovTef = Parameter(name, 'TclovTef', 1.7, 1.0, 10.0, 0.0, hugeval, frozen=True, aliases=["TclTef"])
        self.refflag = Parameter(name, 'refflag', 1.0, -1.0, 1.0, -hugeval, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.D, self.i, self.Mass, self.Mdot, self.TclovTef, self.refflag, self.norm))


@version_at_least("12.10.0")
class XSgrbcomp(XSAdditiveModel):
    """The XSPEC grbcomp model: Comptonization for GRB prompt emission.

    The model is described at [1]_.

    Attributes
    ----------
    kTs
        Temperature of the seed blackbody spectrum in keV.
    gamma
        If set to 3 the seed soft spectrum is a blackbody, otherwise it
        approximates a modified blackbody.
    kTe
        Electron temperature of the subrelativistic outflow in keV.
    tau
        Radial optical depth of the subrelativistic outflow.
    beta
        Bulk outflow velocity of the thermal electrons.
    fbflag
        If set to 0 then only the first-order bulk Comptonization term is
        considered, otherwise if set to 1 then the second-order term
        is computed (see [1]_ for more details).
    log_A
        The geometrical covering factor which determines the relative
        weights of the seed and comptonized spectra to the total flux.
    z
        Redshift.
    a_boost
        The energy index of the Green's function with which the formerly
        comptonization spectrum is convolved.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSgrbm

    Notes
    -----
    This model is only available when used with XSPEC 12.10.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGrbcomp.html

    """

    __function__ = "xsgrbcomp"

    def __init__(self, name='grbcomp'):
        self.kTs = Parameter(name, 'kTs', 1.0, 0.0, 20., 0.0, 20.0, units='keV')
        self.gamma = Parameter(name, 'gamma', 3.0, 0.0, 10.0, 0.0, 10.0)
        self.kTe = Parameter(name, 'kTe', 100.0, 0.2, 2000., 0.2, 2000., units='keV')
        self.tau = Parameter(name, 'tau', 5.0, 0.0, 200., 0.0, 200.)
        self.beta = Parameter(name, 'beta', 0.2, 0.0, 1.0, 0.0, 1.0)
        self.fbflag = Parameter(name, 'fbflag', 0.0, 0.0, 1.0, 0.0, 1.0, frozen=True)
        self.log_A = Parameter(name, 'log_A', 5.0, -8., 8., -8., 8., frozen=True)
        self.z = Parameter(name, 'z', 0.0, 0.0, 10., 0.0, 10., frozen=True)
        self.a_boost = Parameter(name, 'a_boost', 5.0, 0., 30., 0., 30., frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kTs, self.gamma, self.kTe, self.tau, self.beta, self.fbflag, self.log_A, self.z, self.a_boost, self.norm))


class XSgrbm(XSAdditiveModel):
    """The XSPEC grbm model: gamma-ray burst continuum.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``temp`` parameter has been renamed ``tem`` to match the
       XSPEC definition. The name ``temp`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    alpha
        The first powerlaw index.
    beta
        The second powerlaw index.
    tem
        The characteristic energy, in keV.
    norm
        The normalization of the model.

    See Also
    --------
    XSgrbcomp

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGrbm.html

    """

    __function__ = "xsgrbm"

    def __init__(self, name='grbm'):
        self.alpha = Parameter(name, 'alpha', -1., -3., +2., -hugeval, hugeval)
        self.beta = Parameter(name, 'beta', -2., -5., +2., -hugeval, hugeval)
        self.tem = Parameter(name, 'tem', +300., +50., +1000., 0.0, hugeval, 'keV', aliases=["temp"])
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.alpha, self.beta, self.tem, self.norm))


class XSkerrbb(XSAdditiveModel):
    """The XSPEC kerrbb model: multi-temperature blackbody model for thin accretion disk around a Kerr black hole.

    The model is described at [1]_.

    Attributes
    ----------
    eta
        The ratio of the disk power produced by a torque at the disk
        inner boundary to the disk power arising from accretion. See
        [1]_ for more details.
    a
        The specific angular momentum of the black hole in units of the
        black hole mass M (when G=c=1). It should be in the range [0, 1).
    i
        The disk inclination angle, in degrees. A face-on disk has
        i=0. It must be less than or equal to 85 degrees.
    Mbh
        The mass of the black hole, in solar masses.
    Mdd
        The "effective" mass accretion rate in units of 10^18 g/s.
        See [1]_ for more details.
    Dbh
        The distance from the observer to the black hole, in units of kpc.
    hd
        The spectral hardening factor, Tcol/Teff. See [1]_ for more
        details.
    rflag
        A flag to switch on or off the effect of self irradiation:
        when greater than zero the self irradition is included,
        otherwise it is not. This parameter can not be thawed.
    lflag
        A flag to switch on or off the effect of limb darkening:
        when greater than zero the disk emission is assumed to be
        limb darkened, otherwise it is isotropic.
        This parameter can not be thawed.
    norm
        The normalization of the model. It should be fixed to 1
        if the inclination, mass, and distance are frozen.

    See Also
    --------
    XSgrad

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKerrbb.html

    """

    __function__ = "C_kerrbb"

    def __init__(self, name='kerrbb'):
        self.eta = Parameter(name, 'eta', 0., 0., 1.0, 0.0, hugeval, frozen=True)
        self.a = Parameter(name, 'a', 0., -1., 0.9999, -hugeval, hugeval)
        self.i = Parameter(name, 'i', 30., 0., 85., 0.0, hugeval, 'deg', True)
        self.Mbh = Parameter(name, 'Mbh', 1., 0., 100., 0.0, hugeval, 'Msun')
        self.Mdd = Parameter(name, 'Mdd', 1., 0., 1000., 0.0, hugeval, 'Mdd0')
        self.Dbh = Parameter(name, 'Dbh', 10., 0., 10000., 0.0, hugeval, 'kpc', True)
        self.hd = Parameter(name, 'hd', 1.7, 1., 10., 0.0, hugeval, frozen=True)
        self.rflag = Parameter(name, 'rflag', 1., -100., 100., -hugeval, hugeval, frozen=True)
        self.lflag = Parameter(name, 'lflag', 0., -100., 100., -hugeval, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.eta, self.a, self.i, self.Mbh, self.Mdd, self.Dbh, self.hd, self.rflag, self.lflag, self.norm))


class XSkerrd(XSAdditiveModel):
    """The XSPEC kerrd model: optically thick accretion disk around a Kerr black hole.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``TcolTeff`` parameter has been renamed ``TcoloTeff`` to match the
       XSPEC definition. The name ``TcolTeff`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    distance
        The distance, in units of kpc.
    TcoloTeff
        The spectral hardening factor, Tcol/Teff. See [1]_ for more
        details.
    M
        The mass of the central object, in solar masses.
    Mdot
        The mass accretion rate in units of 10^18 g/s.
    Incl
        The disk inclination angle, in degrees. A face-on disk has
        Incl=0.
    Rin
        The inner radius, in units of GM/c^2. The last stable orbit
        is 1.235.
    Rout
        The outer radius, in units of GM/c^2.
    norm
        The normalization of the model. It should be fixed to 1.

    See Also
    --------
    XSlaor

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKerrd.html

    """

    __function__ = "C_kerrd" if equal_or_greater_than("12.10.0") else "C_kerrdisk"

    def __init__(self, name='kerrd'):
        self.distance = Parameter(name, 'distance', 1., 0.01, 1000., 0.0, hugeval, 'kpc', True)
        self.TcoloTeff = Parameter(name, 'TcoloTeff', 1.5, 1.0, 2.0, 0.0, hugeval, frozen=True, aliases=["TcolTeff"])
        self.M = Parameter(name, 'M', 1.0, 0.1, 100., 0.0, hugeval, 'solar')
        self.Mdot = Parameter(name, 'Mdot', 1.0, 0.01, 100., 0.0, hugeval, '1e18')
        self.Incl = Parameter(name, 'Incl', 30., 0., 90., 0.0, hugeval, 'deg', True)
        self.Rin = Parameter(name, 'Rin', 1.235, 1.235, 100., 0.0, hugeval, 'Rg', True)
        self.Rout = Parameter(name, 'Rout', 1e5, 1e4, 1e8, 0.0, hugeval, 'Rg', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.distance, self.TcoloTeff, self.M, self.Mdot, self.Incl, self.Rin, self.Rout, self.norm))


class XSkerrdisk(XSAdditiveModel):
    """The XSPEC kerrdisk model: accretion disk line emission with BH spin as free parameter.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``r_brg``, ``Rinms``, and ``Routms`` parameters have been
       renamed to ``r_br_g``, ``Rin_ms``, and ``Rout_ms`` respectively
       to match the XSPEC definition. The names ``r_brg``, ``Rinms``,
       and ``Routms`` can still be used to access the parameters, but
       they will be removed in a future release.

    Attributes
    ----------
    lineE
        The rest-frame line energy, in keV.
    Index1
        The emissivity index for the inner disk.
    Index2
        The emissivity index for the outer disk.
    r_br_g
        The break radius separating the inner and outer portions of the
        disk, in gravitational radii.
    a
        The dimensionless black hole spin.
    Incl
        The disk inclination angle, in degrees. A face-on disk has
        Incl=0.
    Rin_ms
        The inner radius of the disk, in units of the radius of
        marginal stability.
    Rout_ms
        The outer radius of the disk, in units of the radius of
        marginal stability.
    z
        The redshift of the source.
    norm
        The flux in the line, in units of photon/cm^2/s.

    See Also
    --------
    XSdiskline, XSlaor

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKerrdisk.html

    """

    __function__ = "spin"

    def __init__(self, name='kerrdisk'):
        self.lineE = Parameter(name, 'lineE', 6.4, 0.1, 100., 0.0, hugeval, 'keV', frozen=True)
        self.Index1 = Parameter(name, 'Index1', 3., -10., 10., -hugeval, hugeval, frozen=True)
        self.Index2 = Parameter(name, 'Index2', 3., -10., 10., -hugeval, hugeval, frozen=True)
        self.r_br_g = Parameter(name, 'r_br_g', 6.0, 1.0, 400., 0.0, hugeval, frozen=True, aliases=["r_brg"])
        self.a = Parameter(name, 'a', 0.998, 0.01, 0.998, 0.01, hugeval)
        self.Incl = Parameter(name, 'Incl', 30., 0., 90., 0.0, hugeval, 'deg', True)
        self.Rin_ms = Parameter(name, 'Rin_ms', 1.0, 1.0, 400., 0.0, hugeval, frozen=True, aliases=["Rinms"])
        self.Rout_ms = Parameter(name, 'Rout_ms', 400., 1.0, 400., 0.0, hugeval, frozen=True, aliases=["Routms"])
        self.z = Parameter(name, 'z', 0., 0., 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.lineE, self.Index1, self.Index2, self.r_br_g, self.a, self.Incl, self.Rin_ms, self.Rout_ms, self.z, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.lineE, **kwargs)


@version_at_least("12.10.1")
class XSkyrline(XSAdditiveModel):
    """The XSPEC kyrline model: relativistic line from axisymmetric accretion disk

    The model is described at [1]_.

    Attributes
    ----------
    a
        Black Hole angular momentum, in units of GM/c.
    theta_o
        The observer inclination, where 0 is pole on. The units are
        degrees.
    rin
        The inner radius (units of GM/c^2).
    ms
        A flag that determines whether to integrate from rin (0)
        or to integrate emission from above the marginally-stable orbit
        only (1)
    rout
        The outer radius (units of GM/c^2).
    Erest
        The line energy in keV.
    alpha
        The accretion disk emissivity scales as r^(-alpha) for r < rb.
    beta
        The accretion disk emissivity scales as r^(-beta) for r > rb.
    rb
        The boundary radius between the inner and outer emissivity laws
        (units of GM/c^2).
    zshift
        The overall Doppler shift.
    limb
        0 means isotropic emission, 1 means Laor's limb darkening
        (1 + 0.26 \\mu), 2 means Haardt's limb brightening (ln(1 + 1/\\mu))
    norm
        The normalization.

    Notes
    -----
    This model is only available when used with XSPEC 12.10.1 or later.

    Early releases of XSPEC 12.10.1 refer to the first parameter as
    a/M. As this is not a valid Python name, the parameter has been
    renamed "a" to better match other XSPEC models (after consultation
    with Keith Arnaud).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKyrline.html

    """

    __function__ = "kyrline"

    def __init__(self, name='kyrline'):
        self.a = Parameter(name, 'a', 0.9982, 0, 1, 0, 1, 'GM/c')
        self.theta_o = Parameter(name, 'theta_o', 30, 0, 89, 0, 89, 'deg')
        self.rin = Parameter(name, 'rin', 1, 1, 1000, 1, 1000, 'GM/c^2',
                             frozen=True)
        self.ms = Parameter(name, 'ms', 1, 0, 1, 0, 1, alwaysfrozen=True)
        self.rout = Parameter(name, 'rout', 400, 1, 1000, 1, 1000, 'GM/c^2',
                              frozen=True)
        self.Erest = Parameter(name, 'Erest', 6.4, 1, 99, 1, 99, 'keV',
                               frozen=True)
        self.alpha = Parameter(name, 'alpha', 3, -20, 20, -20, 20, frozen=True)
        self.beta = Parameter(name, 'beta', 3, -20, 20, -20, 20, frozen=True)
        self.rb = Parameter(name, 'rb', 400, 1, 1000, 1, 1000, 'GM/c^2',
                            frozen=True)
        self.zshift = Parameter(name, 'zshift', 0, -0.999, 10, -0.999, 10,
                                frozen=True)
        self.limb = Parameter(name, 'limb', 1, 0, 2, 0, 2, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.a, self.theta_o, self.rin, self.ms, self.rout,
                self.Erest, self.alpha, self.beta, self.rb, self.zshift,
                self.limb, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


class XSlaor(XSAdditiveModel):
    """The XSPEC laor model: accretion disk, black hole emission line.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``RinG`` and ``RoutG`` parameters have been renamed to
       ``Rin_G`` and ``Rout_G`` respectively to match the XSPEC definition.
       The names ``RinG`` and ``RoutG`` can still be used to access the
       parameters, but they will be removed in a future release.

    Attributes
    ----------
    lineE
        The rest-frame line energy, in keV.
    Index
        The power law dependence of emissivity (scales as R^-Index).
    Rin_G
        The inner radius, in units of GM/c^2.
    Rout_G
        The outer radius, in units of GM/c^2.
    Incl
        The disk inclination angle, in degrees. A face-on disk has
        Incl=0.
    norm
        The flux in the line, in units of photon/cm^2/s.

    See Also
    --------
    XSlaor2

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLaor.html

    """

    __function__ = "C_laor" if equal_or_greater_than("12.10.1") else "C_xslaor"

    def __init__(self, name='laor'):
        self.lineE = Parameter(name, 'lineE', 6.4, 0., 100., 0.0, hugeval, 'keV')
        self.Index = Parameter(name, 'Index', 3., -10., 10., -hugeval, hugeval, frozen=True)
        self.Rin_G = Parameter(name, 'Rin_G', 1.235, 1.235, 400., 0.0, hugeval, frozen=True, aliases=["RinG"])
        self.Rout_G = Parameter(name, 'Rout_G', 400., 1.235, 400., 0.0, hugeval, frozen=True, aliases=["RoutG"])
        self.Incl = Parameter(name, 'Incl', 30., 0., 90., 0.0, hugeval, 'deg', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.lineE, self.Index, self.Rin_G, self.Rout_G, self.Incl, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.lineE, **kwargs)


class XSlaor2(XSAdditiveModel):
    """The XSPEC laor2 model: accretion disk with broken-power law emissivity profile, black hole emission line.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``RinG`` and ``RoutG`` parameters have been renamed to
       ``Rin_G`` and ``Rout_G`` respectively to match the XSPEC definition.
       The names ``RinG`` and ``RoutG`` can still be used to access the
       parameters, but they will be removed in a future release.

    Attributes
    ----------
    lineE
        The rest-frame line energy, in keV.
    Index
        The power law dependence of emissivity (scales as R^-Index).
    Rin_G
        The inner radius, in units of GM/c^2.
    Rout_G
        The outer radius, in units of GM/c^2.
    Incl
        The disk inclination angle, in degrees. A face-on disk has
        Incl=0.
    Rbreak
        The radius at which the emissivity power-law index changes.
    Index1
        The emissivity power-law index for r>Rbreak.
    norm
        The flux in the line, in units of photon/cm^2/s.

    See Also
    --------
    XSlaor

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLaor2.html

    """

    __function__ = "C_laor2"

    def __init__(self, name='laor2'):
        self.lineE = Parameter(name, 'lineE', 6.4, 0., 100., 0.0, hugeval, 'keV')
        self.Index = Parameter(name, 'Index', 3., -10., 10., -hugeval, hugeval, frozen=True)
        self.Rin_G = Parameter(name, 'Rin_G', 1.235, 1.235, 400., 0.0, hugeval, frozen=True, aliases=["RinG"])
        self.Rout_G = Parameter(name, 'Rout_G', 400., 1.235, 400., 0.0, hugeval, frozen=True, aliases=["RoutG"])
        self.Incl = Parameter(name, 'Incl', 30., 0., 90., 0.0, hugeval, 'deg', True)
        self.Rbreak = Parameter(name, 'Rbreak', 20., 1.235, 400., 0.0, hugeval, frozen=True)
        self.Index1 = Parameter(name, 'Index1', 3., -10., 10., -hugeval, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.lineE, self.Index, self.Rin_G, self.Rout_G, self.Incl, self.Rbreak, self.Index1, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.lineE, **kwargs)


class XSlorentz(XSAdditiveModel):
    """The XSPEC lorentz model: lorentz line profile.

    The model is described at [1]_.

    Attributes
    ----------
    LineE
        The line energy, in keV.
    Width
        The FWHM of the line, in keV.
    norm
        The flux in the line, in units of photon/cm^2/s.

    See Also
    --------
    XSgaussian, XSvoigt

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLorentz.html

    """

    __function__ = "C_lorentzianLine" if equal_or_greater_than("12.9.1") else "xslorz"

    def __init__(self, name='lorentz'):
        self.LineE = Parameter(name, 'LineE', 6.5, 0., 1.e6, 0.0, hugeval, 'keV')
        self.Width = Parameter(name, 'Width', 0.1, 0., 10., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Width, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


class XSmeka(XSAdditiveModel):
    """The XSPEC meka model: emission, hot diffuse gas (Mewe-Gronenschild).

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    nH
        H density, in cm^-3.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSvmeka

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelMeka.html

    """

    __function__ = "C_meka" if equal_or_greater_than("12.9.1") else "xsmeka"

    def __init__(self, name='meka'):
        self.kT = Parameter(name, 'kT', 1., 1.e-3, 1.e2, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1., 1.e-5, 1.e19, 0.0, hugeval, 'cm-3', True)
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.nH, self.Abundanc, self.redshift, self.norm))


class XSmekal(XSAdditiveModel):
    """The XSPEC mekal model: emission, hot diffuse gas (Mewe-Kaastra-Liedahl).

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    nH
        H density, in cm^-3.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSapec, XSvmekal

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelMekal.html

    """

    __function__ = "C_mekal" if equal_or_greater_than("12.9.1") else "xsmekl"

    def __init__(self, name='mekal'):
        self.kT = Parameter(name, 'kT', 1., 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1., 1.e-5, 1.e19, 0.0, hugeval, 'cm-3', True)
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kT, self.nH, self.Abundanc, self.redshift, self.switch, self.norm))


class XSmkcflow(XSAdditiveModel):
    """The XSPEC mkcflow model: cooling flow, mekal.

    The model is described at [1]_. The results of this model depend
    on the cosmology settings set with ``set_xscosmo``.

    Attributes
    ----------
    lowT
        The minimum temperature, in keV.
    highT
        The maxmimum temperature, in keV.
    Abundanc
        The abundance relative to Solar, as set by ``set_xsabund``.
    redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The mass accretion rate (solar mass per year).

    See Also
    --------
    XSapec, XScflow, XScevmkl, XScph, XSvcph, XSvmcflow

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelMkcflow.html

    """

    __function__ = "C_xsmkcf"

    def __init__(self, name='mkcflow'):
        self.lowT = Parameter(name, 'lowT', 0.1, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.highT = Parameter(name, 'highT', 4., 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.lowT, self.highT, self.Abundanc, self.redshift, self.switch,
                                              self.norm))


class XSnei(XSAdditiveModel):
    """The XSPEC nei model: collisional plasma, non-equilibrium, constant temperature.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    Tau
        The ionization timescale in units of s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSgnei, XSvnei, XSvvnei

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNei.html

    """

    __function__ = "C_nei"

    def __init__(self, name='nei'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.Tau, self.redshift, self.norm))


class XSrnei(XSAdditiveModel):
    """The XSPEC rnei model: non-equilibrium recombining collisional plasma.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    kT_init
        The initial temperature of the plasma, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    Tau
        The ionization timescale in units of s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSnei, XSgnei, XSvrnei, XSvvrnei

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRnei.html

    """

    __function__ = "C_rnei"

    def __init__(self, name='rnei'):
        self.kT = Parameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = Parameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.Abundanc, self.Tau, self.redshift, self.norm))


class XSvrnei(XSAdditiveModel):
    """The XSPEC vrnei model: non-equilibrium recombining collisional plasma.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    kT_init
        The initial temperature of the plasma, in keV.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element, with respect to Solar.
    Tau
        The ionization timescale in units of s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSrnei, XSvvrnei

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRnei.html

    """

    __function__ = "C_vrnei"

    def __init__(self, name='vrnei'):
        self.kT = Parameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = Parameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = Parameter(name, 'H', 1.0, 0., 1., 0.0, 10000.0, frozen=True)
        self.He = Parameter(name, 'He', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.redshift, self.norm))


class XSvvrnei(XSAdditiveModel):
    """The XSPEC vvrnei model: non-equilibrium recombining collisional plasma.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    kT_init
        The initial temperature of the plasma, in keV.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element, with respect to Solar.
    Tau
        The ionization timescale in units of s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSrnei, XSvrnei

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRnei.html

    """

    __function__ = "C_vvrnei"

    def __init__(self, name='vvrnei'):
        self.kT = Parameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = Parameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.redshift, self.norm))


class XSnpshock(XSAdditiveModel):
    """The XSPEC npshock model: shocked plasma, plane parallel, separate ion, electron temperatures.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT_a
        The mean shock temperature, in keV.
    kT_b
        The electron temperature immediately behind the shock
        front, in keV. See [1]_ for a discussion of the behavior
        of kT_a and kT_b.
    Abundanc
        The metal abundance, as defined by the
        ``set_xsabund`` function.
    Tau_l
        The lower limit on the ionization timescale in s/cm^3.
    Tau_u
        The upper limit on the ionization timescale in s/cm^3.
    redshift
        The redshift of the source.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSsedov, XSvnpshock, XSvvnpshock

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNpshock.html

    """

    __function__ = "C_npshock"

    def __init__(self, name='npshock'):
        self.kT_a = Parameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.kT_b = Parameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau_l = Parameter(name, 'Tau_l', 0.0, 0., 5.e13, 0.0, hugeval, 's/cm^3', True)
        self.Tau_u = Parameter(name, 'Tau_u', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.Abundanc, self.Tau_l, self.Tau_u, self.redshift, self.norm))


class XSnsa(XSAdditiveModel):
    """The XSPEC nsa model: neutron star atmosphere.

    The model is described at [1]_.

    Attributes
    ----------
    LogT_eff
        The log of Teff, the unredshifted effective temperature.
    M_ns
        The neutron star gravitational mass, in units of the solar mass.
    R_ns
        The neutron star radius, in km.
    MagField
        The neutron star magnetic field strength, in Gauss. It must be
        fixed at one of 0, 1e12, or 1e13 G.
    norm
        The normalization is 1/D^2, where D is the distance to the
        neutron star in pc.

    See Also
    --------
    XSnsagrav

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNsa.html

    """

    __function__ = "nsa"

    def __init__(self, name='nsa'):
        self.LogT_eff = Parameter(name, 'LogT_eff', 6.0, 5.0, 7.0, 0.0, hugeval, 'K')
        self.M_ns = Parameter(name, 'M_ns', 1.4, 0.5, 2.5, 0.0, hugeval, 'Msun')
        self.R_ns = Parameter(name, 'R_ns', 10.0, 5.0, 20., 0.0, hugeval, 'km')
        self.MagField = Parameter(name, 'MagField', 0.0, 0.0, 5.e13, 0.0, hugeval, 'G', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LogT_eff, self.M_ns, self.R_ns, self.MagField, self.norm))


class XSnsagrav(XSAdditiveModel):
    """The XSPEC nsagrav model: NS H atmosphere model for different g.

    The model is described at [1]_.

    Attributes
    ----------
    LogT_eff
        The log of Teff, the unredshifted effective temperature.
    NSmass
        The neutron star gravitational mass, in units of the solar mass.
    NSrad
        The "true" neutron star radius, in km.
    norm
        The normalization is 1/D^2, where D is the distance to the
        neutron star in pc.

    See Also
    --------
    XSnsa

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNsagrav.html

    """

    __function__ = "nsagrav"

    def __init__(self, name='nsagrav'):
        self.LogT_eff = Parameter(name, 'LogT_eff', 6.0, 5.5, 6.5, 0.0, hugeval, 'K')
        self.NSmass = Parameter(name, 'NSmass', 1.4, 0.3, 2.5, 0.0, hugeval, 'Msun')
        self.NSrad = Parameter(name, 'NSrad', 10.0, 6.0, 20., 0.0, hugeval, 'km')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LogT_eff, self.NSmass, self.NSrad, self.norm))


class XSnsatmos(XSAdditiveModel):
    """The XSPEC nsatmos model: NS Hydrogen Atmosphere model with electron conduction and self-irradiation.

    The model is described at [1]_.

    Attributes
    ----------
    LogT_eff
        The log of Teff, the unredshifted effective temperature.
    M_ns
        The neutron star gravitational mass, in units of the solar mass.
    R_ns
        The "true" neutron star radius, in km.
    dist
        The distance to the neutron star, in kpc.
    norm
        The fraction of the neutron star surface emitting.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNsatmos.html

    """

    __function__ = "nsatmos"

    def __init__(self, name='nsatmos'):
        self.LogT_eff = Parameter(name, 'LogT_eff', 6.0, 5.0, 6.5, 0.0, hugeval, 'K')
        self.M_ns = Parameter(name, 'M_ns', 1.4, 0.5, 3.0, 0.0, hugeval, 'Msun')
        self.R_ns = Parameter(name, 'R_ns', 10.0, 5.0, 30., 0.0, hugeval, 'km')
        self.dist = Parameter(name, 'dist', 10.0, 0.1, 100.0, 0.0, hugeval, 'kpc')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LogT_eff, self.M_ns, self.R_ns, self.dist, self.norm))


class XSnsmax(XSAdditiveModel):
    """The XSPEC nsmax model: Neutron Star Magnetic Atmosphere.

    The model is described at [1]_. It has been superceeded by
    ``XSnsmaxg``.

    Attributes
    ----------
    LogTeff
        The log of Teff, the unredshifted surface effective temperature.
    redshift
        This is 1 + zg, the gravitational redshift.
    specfile
        Which model to use: see [1]_ for more details.
    norm
        The normalization of the model: see [1]_ for more details.

    See Also
    --------
    XSnsmaxg

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNsmax.html

    """

    __function__ = _f77_or_c_12100("nsmax")

    def __init__(self, name='nsmax'):
        self.logTeff = Parameter(name, 'logTeff', 6.0, 5.5, 6.8, 0.0, hugeval, 'K')
        self.redshift = Parameter(name, 'redshift', 0.1, 1.0e-5, 1.5, 1.0e-5, 2.0)
        self.specfile = Parameter(name, 'specfile', 1200, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.logTeff, self.redshift, self.specfile, self.norm))


class XSnsmaxg(XSAdditiveModel):
    """The XSPEC nsmaxg model: neutron star with a magnetic atmosphere.

    The model is described at [1]_.

    Attributes
    ----------
    LogTeff
        The log of Teff, the unredshifted surface effective temperature.
    M_ns
        The neutron star gravitational mass, in units of the solar mass.
    R_ns
        The neutron star radius, in km.
    dist
        The distance to the neutron star, in kpc.
    specfile
        Which model to use: see [1]_ for more details.
    norm
        The normalization: see [1]_ for more details.

    See Also
    --------
    XSnsmax, XSnsx

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNsmaxg.html

    """

    __function__ = _f77_or_c_12100("nsmaxg")

    def __init__(self, name='nsmaxg'):
        self.logTeff = Parameter(name, 'logTeff', 6.0, 5.5, 6.9, 5.5, 6.9, units='K')
        self.M_ns = Parameter(name, 'M_ns', 1.4, 0.5, 3.0, 0.5, 3.0, units='Msun')
        self.R_ns = Parameter(name, 'R_ns', 10.0, 5.0, 30.0, 5.0, 30.0, units='km')
        self.dist = Parameter(name, 'dist', 1.0, 0.01, 100.0, 0.01, 100.0, units='kpc')
        self.specfile = Parameter(name, 'specfile', 1200, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.logTeff, self.M_ns, self.R_ns, self.dist, self.specfile, self.norm))


class XSnsx(XSAdditiveModel):
    """The XSPEC nsx model: neutron star with a non-magnetic atmosphere.

    The model is described at [1]_.

    Attributes
    ----------
    LogTeff
        The log of Teff, the unredshifted surface effective temperature.
    M_ns
        The neutron star gravitational mass, in units of the solar mass.
    R_ns
        The neutron star radius, in km.
    dist
        The distance to the neutron star, in kpc.
    specfile
        Which model to use: see [1]_ for more details.
    norm
        The normalization: see [1]_ for more details.

    See Also
    --------
    XSnsmaxg

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNsx.html

    """

    __function__ = _f77_or_c_12100("nsx")

    def __init__(self, name='nsx'):
        self.logTeff = Parameter(name, 'logTeff', 6.0, 5.5, 6.7, 5.5, 6.7, units='K')
        self.M_ns = Parameter(name, 'M_ns', 1.4, 0.5, 3.0, 0.5, 3.0, units='Msun')
        self.R_ns = Parameter(name, 'R_ns', 10.0, 5.0, 30.0, 5.0, 30.0, units='km')
        self.dist = Parameter(name, 'dist', 1.0, 0.01, 100.0, 0.01, 100.0, units='kpc')
        self.specfile = Parameter(name, 'specfile', 6, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.logTeff, self.M_ns, self.R_ns, self.dist, self.specfile, self.norm))


class XSnteea(XSAdditiveModel):
    """The XSPEC nteea model: non-thermal pair plasma.

    The model is described at [1]_.

    Attributes
    ----------
    l_nth
        The nonthermal electron compactness.
    l_bb
        The blackbody compactness.
    f_refl
        The scaling factor for reflection. This is 1 for an isotropic
        source above the disk.
    kT_bb
        The blackbody temperature in eV.
    g_max
        The maximum Lorentz factor.
    l_th
        The thermal compactness. Set to 0 for a pure nonthermal plasma.
    tau_p
        The Thomson optical depth of ionization electrons.
    G_inj
        The electron injection index (0 for monoenergetic injection).
    g_min
        The minimum Lorentz factor of the power law injection (not used
        for monoenergetic injection).
    g_0
        The minimum Lorentz factor for nonthermal reprocessing, in the
        range (1, g_min].
    radius
        The radius in cm (for Coulomb/bremsstrahlung only).
    pair_esc
        The pair escape rate in c.
    cosIncl
        The cosine of the inclination angle.
    Fe_abund
        The iron abundance relative to that set by the ``set_xsabund``
        function.
    redshift
        The redshift of the source.
    norm
        The normalization of the model: see [1]_ for more details.

    Notes
    -----
    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the NTEEA_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNteea.html

    """

    __function__ = "C_xsnteea"

    def __init__(self, name='nteea'):
        self.l_nth = Parameter(name, 'l_nth', 100., 0., 1.e4, 0.0, hugeval)
        self.l_bb = Parameter(name, 'l_bb', 100., 0., 1.e4, 0.0, hugeval)
        self.f_refl = Parameter(name, 'f_refl', 0., 0., 4., 0.0, hugeval)
        self.kT_bb = Parameter(name, 'kT_bb', 10., 1., 100., 0.0, hugeval, frozen=True)
        self.g_max = Parameter(name, 'g_max', 1.e3, 5., 1.e4, 0.0, hugeval, frozen=True)
        self.l_th = Parameter(name, 'l_th', 0., 0., 1.e4, 0.0, hugeval, frozen=True)
        self.tau_p = Parameter(name, 'tau_p', 0., 0., 10., 0.0, hugeval, frozen=True)
        self.G_inj = Parameter(name, 'G_inj', 0., 0., 5., 0.0, hugeval, frozen=True)
        self.g_min = Parameter(name, 'g_min', 1.3, 1., 1.e3, 0.0, hugeval, frozen=True)
        self.g_0 = Parameter(name, 'g_0', 1.3, 1., 5., 0.0, hugeval, frozen=True)
        self.radius = Parameter(name, 'radius', 1.e13, 1.e5, 1.e16, 0.0, hugeval, frozen=True)
        self.pair_esc = Parameter(name, 'pair_esc', 0., 0., 1., 0.0, hugeval, frozen=True)
        self.cosIncl = Parameter(name, 'cosIncl', 0.45, 0.05, 0.95, 0.0, hugeval)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.1, 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.l_nth, self.l_bb, self.f_refl, self.kT_bb, self.g_max, self.l_th, self.tau_p, self.G_inj, self.g_min, self.g_0, self.radius, self.pair_esc, self.cosIncl, self.Fe_abund, self.redshift, self.norm))


class XSnthComp(XSAdditiveModel):
    """The XSPEC nthComp model: Thermally comptonized continuum.

    The model is described at [1]_.

    Attributes
    ----------
    Gamma
        The asymptotic power-law photon index.
    kT_e
        The electron temperature (high-energy rollover) in keV.
    kT_bb
        The seed photon temperature (low-energy rollover) in keV.
    inp_type
        The seed photon source: when 0 use a blackbody, when 1 use
        a disk-blackbody.
    redshift
        The redshift.
    norm
        The normalization of the model: see [1]_ for more details.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNthcomp.html

    """

    __function__ = "C_nthcomp"

    def __init__(self, name='nthcomp'):
        self.Gamma = Parameter(name, 'Gamma', 1.7, 1.001, 5., 1.001, 10.)
        self.kT_e = Parameter(name, 'kT_e', 100., 5., 1.e3, 1., 1.e3, 'keV')
        self.kT_bb = Parameter(name, 'kT_bb', 0.1, 1.e-3, 10., 1.e-3, 10., 'keV', True)
        self.inp_type = Parameter(name, 'inp_type', 0., 0., 1., 0., 1., '0/1', True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, 10., frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Gamma, self.kT_e, self.kT_bb, self.inp_type, self.redshift, self.norm))


class XSpegpwrlw(XSAdditiveModel):
    """The XSPEC pegpwrlw model: power law, pegged normalization.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndex
        The power law photon index.
    eMin
        The lower peg for the normalization, in keV.
    eMax
        The upper peg for the normalization, in keV.
    norm
        The flux, in units of 10^-12 erg/cm^2/s, over the range
        eMin to eMax. If eMin=eMax then it is the flux in
        micro-Jy at eMin.

    See Also
    --------
    XSbknpower, XSbkn2pow, XScutoffpl, XSpowerlaw, XSzpowerlw

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPegpwrlw.html

    """

    __function__ = "xspegp"

    def __init__(self, name='pegpwrlw'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 1., -2., 9., -hugeval, hugeval)
        self.eMin = Parameter(name, 'eMin', 2., -100., 1.e10, -hugeval, hugeval, 'keV', True)
        self.eMax = Parameter(name, 'eMax', 10., -100., 1.e10, -hugeval, hugeval, 'keV', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.eMin, self.eMax, self.norm))


class XSpexrav(XSAdditiveModel):
    """The XSPEC pexrav model: reflected powerlaw, neutral medium.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndex
        The first power-law photon index.
    foldE
        The cut-off energy (E_c) in keV. Set to 0 for no cut off.
    rel_refl
        The reflection scaling parameter (a value between 0 and 1
        for an isotropic source above the disk, less than 0 for no
        reflected component).
    redshift
        The redshift of the source.
    abund
        The abundance of the elements heaver than He relative to their
        solar abundance, as set by the ``set_xsabund`` function.
    Fe_abund
        The iron abundance relative to the solar abundance, as set by
        the ``set_xsabund`` function.
    cosIncl
        The cosine of the inclination angle in degrees.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSpexriv, XSpexmon

    Notes
    -----
    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the PEXRAV_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPexrav.html

    """

    __function__ = "C_xspexrav"

    def __init__(self, name='pexrav'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 2., -9., 9., -hugeval, hugeval)
        self.foldE = Parameter(name, 'foldE', 100., 1., 1.e6, 0.0, hugeval, 'keV')
        self.rel_refl = Parameter(name, 'rel_refl', 0., 0., 1.e6, -1.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.abund = Parameter(name, 'abund', 1., 0.0, 1.e6, 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.0, 1.e6, 0.0, hugeval, frozen=True)
        self.cosIncl = Parameter(name, 'cosIncl', 0.45, 0.05, 0.95, 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.foldE, self.rel_refl, self.redshift, self.abund, self.Fe_abund, self.cosIncl, self.norm))


class XSpexriv(XSAdditiveModel):
    """The XSPEC pexriv model: reflected powerlaw, neutral medium.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndex
        The first power-law photon index.
    foldE
        The cut-off energy (E_c) in keV. Set to 0 for no cut off.
    rel_refl
        The reflection scaling parameter (a value between 0 and 1
        for an isotropic source above the disk, less than 0 for no
        reflected component).
    redshift
        The redshift of the source.
    abund
        The abundance of the elements heaver than He relative to their
        solar abundance, as set by the ``set_xsabund`` function.
    Fe_abund
        The iron abundance relative to the solar abundance, as set by
        the ``set_xsabund`` function.
    cosIncl
        The cosine of the inclination angle in degrees.
    T_disk
        The disk temperature in K.
    xi
        The disk ionization parameter: see [1]_ for more details.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSpexrav, XSpexmon

    Notes
    -----
    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the PEXRIV_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPexriv.html

    """

    __function__ = "C_xspexriv"

    def __init__(self, name='pexriv'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 2., -9., 9., -hugeval, hugeval)
        self.foldE = Parameter(name, 'foldE', 100., 1., 1.e6, 0.0, hugeval, 'keV')
        self.rel_refl = Parameter(name, 'rel_refl', 0., 0., 1.e6, -1.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.abund = Parameter(name, 'abund', 1., 0.0, 1.e6, 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.0, 1.e6, 0.0, hugeval, frozen=True)
        self.cosIncl = Parameter(name, 'cosIncl', 0.45, 0.05, 0.95, 0.0, hugeval, frozen=True)
        self.T_disk = Parameter(name, 'T_disk', 3.e4, 1.e4, 1.e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 1., 0., 1.e3, 0.0, hugeval, 'erg cm/s')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.foldE, self.rel_refl, self.redshift, self.abund, self.Fe_abund, self.cosIncl, self.T_disk, self.xi, self.norm))


class XSplcabs(XSAdditiveModel):
    """The XSPEC plcabs model: powerlaw observed through dense, cold matter.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The column density, in units of 10^22 cm^-2.
    nmax
        The maximum number of scatterings. This parameter can not be
        thawed.
    FeAbun
        The iron abundance.
    FeKedge
        The energy of the Fe K edge, in keV.
    PhoIndex
        The power law photon index.
    HighECut
        The high-energy cut-off threshold energy, in keV.
    foldE
        The high-energy cut-off e-folding energy, in keV.
    acrit
        The critical albedo for switching to elastic scattering. See
        [1]_ for more details.
    FAST
        If set to a value above 1, use a mean energy shift instead of
        integration. This parameter can not be thawed.
    redshift
        The redshift of the source.
    norm
        The normalization of the model.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPlcabs.html

    """

    __function__ = "xsp1tr"

    def __init__(self, name='plcabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.nmax = Parameter(name, 'nmax', 1, alwaysfrozen=True)
        self.FeAbun = Parameter(name, 'FeAbun', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.FeKedge = Parameter(name, 'FeKedge', 7.11, 7., 10., 0.0, hugeval, 'KeV', True)
        self.PhoIndex = Parameter(name, 'PhoIndex', 2., -2., 9., -hugeval, hugeval)
        self.HighECut = Parameter(name, 'HighECut', 95., 1., 100., 0.01, 200., 'keV', True)
        self.foldE = Parameter(name, 'foldE', 100., 1., 1.e6, 0.0, hugeval, frozen=True)
        self.acrit = Parameter(name, 'acrit', 1., 0.0, 1.0, 0.0, hugeval, frozen=True)
        self.FAST = Parameter(name, 'FAST', 0, alwaysfrozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.nH, self.nmax, self.FeAbun, self.FeKedge, self.PhoIndex,
                                              self.HighECut, self.foldE, self.acrit, self.FAST, self.redshift,
                                              self.norm))


class XSpowerlaw(XSAdditiveModel):
    """The XSPEC powerlaw model: power law photon spectrum.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndex
        The power law photon index.
    norm
        The normalization of the model. See [1]_ for details, as its
        meaning depends on whether the "POW_EMIN" or "POW_EMAX"
        keywords have been set with ``set_xsxset``.

    See Also
    --------
    XSbknpower, XSbkn2pow, XScutoffpl, XSpegpwrlw, XSzpowerlw

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPowerlaw.html

    """

    __function__ = "C_powerLaw"

    def __init__(self, name='powerlaw'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 1., -2., 9., -hugeval, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.norm))


class XSposm(XSAdditiveModel):
    """The XSPEC posm model: positronium continuum.

    The model is described at [1]_.

    Attributes
    ----------
    norm
        The normalization of the model: see [1]_ for details.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPosm.html

    """

    __function__ = "xsposm"

    def __init__(self, name='posm'):
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.norm,))


class XSpshock(XSAdditiveModel):
    """The XSPEC pshock model: plane-parallel shocked plasma, constant temperature.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    Tau_l
        The lower limit on the ionization timescale, in s/cm^3.
    Tau_u
        The upper limit on the ionization timescale, in s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSvpshock, XSvvpshock

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPshock.html

    """

    __function__ = "C_pshock"

    def __init__(self, name='pshock'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau_l = Parameter(name, 'Tau_l', 0.0, 0., 5.e13, 0.0, hugeval, 's/cm^3', True)
        self.Tau_u = Parameter(name, 'Tau_u', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.Tau_l, self.Tau_u, self.redshift, self.norm))


@version_at_least("12.10.1")
class XSqsosed(XSAdditiveModel):
    """The XSPEC qsosed model: AGN SED model

    The model is described at [1]_.

    Attributes
    ----------
    mass
        The black hole mass, in solar units.
    dist
        The comoving (proper) distance, in Mpc.
    logmdot
        log of mdot, where mdot = Mdot/Mdot_Edd and
        eta Mdot_Edd c^2 = L_Edd
    astar
        The black hole spin (dimensionless)
    cosi
        The cosine of the inclination angle i for the warm
        Comptonising component and the outer disc.
    redshift
        The redshift.
    norm
        The normalization of the model.

    See Also
    --------
    XSagnsed

    Notes
    -----
    This model is only available when used with XSPEC 12.10.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelAgnsed.html

    """

    __function__ = "qsosed"

    def __init__(self, name='qsosed'):
        self.mass = Parameter(name, 'mass', 1e7, 1e5, 1e10, 1e5, 1e10,
                              'solar', frozen=True)
        self.dist = Parameter(name, 'dist', 100, 0.01, 1e9, 0.01, 1e9,
                              'Mpc', frozen=True)
        self.logmdot = Parameter(name, 'logmdot', -1, -1.65, 0.39, -1.65, 0.39)
        self.astar = Parameter(name, 'astar', 0.0, -1, 0.998, -1, 0.998,
                               frozen=True)
        self.cosi = Parameter(name, 'cosi', 0.5, 0.05, 1.0, 0.05, 1.0,
                              frozen=True)

        self.redshift = Parameter(name, 'redshift', 0, 0, 5, 0, 5,
                                  frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.mass, self.dist, self.logmdot, self.astar, self.cosi,
                self.redshift, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


class XSraymond(XSAdditiveModel):
    """The XSPEC raymond model: emission, hot diffuse gas, Raymond-Smith.

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSvraymond

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRaymond.html

    """

    __function__ = "C_raysmith" if equal_or_greater_than("12.9.1") else "xsrays"

    def __init__(self, name='raymond'):
        self.kT = Parameter(name, 'kT', 1., 0.008, 64.0, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.redshift, self.norm))


class XSredge(XSAdditiveModel):
    """The XSPEC redge model: emission, recombination edge.

    The model is described at [1]_.

    Attributes
    ----------
    edge
        The threshold energy, in keV.
    kT
        The plasma temperature, in keV.
    norm
        The flux in the line, in units of photon/cm^2/s.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRedge.html

    """

    __function__ = "xredge"

    def __init__(self, name='redge'):
        self.edge = Parameter(name, 'edge', 1.4, 0.001, 100., 0.0, hugeval, 'keV')
        self.kT = Parameter(name, 'kT', 1., 0.001, 100., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.edge, self.kT, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.edge, **kwargs)


class XSrefsch(XSAdditiveModel):
    """The XSPEC refsch model: reflected power law from ionized accretion disk.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndex
        The power-law photon index.
    foldE
        The cut-off energy (E_c) in keV. Set to 0 for no cut off.
    rel_refl
        The reflection scaling parameter (a value between 0 and 1
        for an isotropic source above the disk, less than 0 for no
        direct component).
    redshift
        The redshift of the source.
    abund
        The abundance of the elements heaver than He relative to their
        solar abundance, as set by the ``set_xsabund`` function.
    Fe_abund
        The iron abundance relative to the solar abundance, as set by
        the ``set_xsabund`` function.
    Incl
        The inclination angle in degrees.
    T_disk
        The disk temperature in K.
    xi
        The disk ionization parameter: see [1]_ for more details.
    Betor10
        The power law dependence of emissivity: see [1]_ for more details.
    Rin
        The inner radius, in units of GM^2/c.
    Rout
        The outer radius, in units of GM^2/c.
    accuracy
        The internal model accuracy: points of spectrum per energy decade.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSdiskline, XSpexriv

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRefsch.html

    """

    __function__ = "xsrefsch"

    def __init__(self, name='refsch'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 2., -9., 9., -hugeval, hugeval)
        self.foldE = Parameter(name, 'foldE', 100., 1., 1.e6, 0.0, hugeval, 'keV')
        self.rel_refl = Parameter(name, 'rel_refl', 0., 0., 2., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.abund = Parameter(name, 'abund', 1., 0.5, 10., 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.1, 10., 0.0, hugeval, frozen=True)
        self.Incl = Parameter(name, 'Incl', 30., 19., 87., 0.0, hugeval, 'deg', True)
        self.T_disk = Parameter(name, 'T_disk', 3.e4, 1.e4, 1.e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 1., 0., 1.e3, 0.0, hugeval, 'erg cm/s')
        self.Betor10 = Parameter(name, 'Betor10', -2., -10., 20., -hugeval, hugeval, frozen=True)
        self.Rin = Parameter(name, 'Rin', 10., 6., 1000., 0.0, hugeval, 'R_g', True)
        self.Rout = Parameter(name, 'Rout', 1000., 0., 1000000., 0.0, hugeval, 'R_g', True)
        self.accuracy = Parameter(name, 'accuracy', 30., 30., 100000., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.foldE, self.rel_refl, self.redshift, self.abund, self.Fe_abund, self.Incl, self.T_disk, self.xi, self.Betor10, self.Rin, self.Rout, self.accuracy, self.norm))


class XSsedov(XSAdditiveModel):
    """The XSPEC sedov model: sedov model, separate ion/electron temperature.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT_a
        The mean shock temperature, in keV.
    kT_b
        The electron temperature immediately behind the shock
        front, in keV. See [1]_ for a discussion of the behavior
        of kT_a and kT_b.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    Tau
        The ionization timescale in units of s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSpshock, XSvsedov, XSvvsedov

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSedov.html

    """

    __function__ = "C_sedov"

    def __init__(self, name='sedov'):
        self.kT_a = Parameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.kT_b = Parameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.Abundanc, self.Tau, self.redshift, self.norm))


class XSsirf(XSAdditiveModel):
    """The XSPEC sirf model: self-irradiated funnel.

    The model is described at [1]_.

    Attributes
    ----------
    tin
        The inner temperature (at the inner, inside-the-funnel
        photosphere), in keV.
    rin
        The inner (inner, inside-the-fulle photosphere) radius in
        "spherisation radius" units (see [1]_ for details).
    rout
        The outer photosphere radius in "spherisation radius" units
    theta
        The half-opening angle of the cone, in degrees.
    incl
        The inclination angle of the funnel, in degrees. Affects mainly
        self-occultation and relativistic boost effects.
    valpha
        The velocity law exponent.
    gamma
        The adiabatic index. It affects the inner, hotter parts of the
        flow, therefore we set is to 4/3 by default.
    mdot
        The mass ejection rate in Eddington (critical) units.
    irrad
        The number of iterations for irradiation.
    norm
        The normalization of the model.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSirf.html

    """

    __function__ = "C_sirf"

    def __init__(self, name='sirf'):
        self.tin = Parameter(name,    'tin', 1., 0.01, 100., 0.01, 1000., 'keV')
        self.rin = Parameter(name,    'rin', 1.e-2, 1.e-5, 1.0, 1.e-6, 10, 'rsph')
        self.rout = Parameter(name,   'rout', 100., 0.1, 1.e8, 0.1, 1.e8, 'rsph')
        self.theta = Parameter(name,  'theta', 22.9, 1., 89., 0., 90., 'deg')
        self.incl = Parameter(name,   'incl', 0., -90., 90., -90., 90., 'deg', True)
        self.valpha = Parameter(name, 'valpha', -0.5, -1.0, 2., -1.5, 5., frozen=True)
        self.gamma = Parameter(name,  'gamma', 1.333, 0.5, 10., 0.5, 10., frozen=True)
        self.mdot = Parameter(name,   'mdot', 1000., 0.5, 1.e7, 0.5, 1.e7, frozen=True)
        self.irrad = Parameter(name,  'irrad', 2., 0., 10., 0., 20., frozen=True)
        self.norm = Parameter(name,   'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.tin, self.rin, self.rout, self.theta,
                                              self.incl, self.valpha, self.gamma, self.mdot,
                                              self.irrad, self.norm))


class XSsrcut(XSAdditiveModel):
    """The XSPEC srcut model: synchrotron spectrum, cutoff power law.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``breakfreq`` parameter has been renamed ``break_`` to match the
       XSPEC definition. The name ``breakfreq`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    alpha
        The radio spectral index.
    break_
        The break frequency: approximately the frequency at which the
        flux has dropped by a factor of 10 from a straight power law.
    norm
        The 1 Ghz flux in Jy.

    See Also
    --------
    XSsresc

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSrcut.html

    """

    __function__ = "srcut"

    def __init__(self, name='srcut'):
        self.alpha = Parameter(name, 'alpha', 0.5, 0.3, 0.8, 0.0, hugeval)
        self.break_ = Parameter(name, 'break_', 2.42E17, 1.E15, 1.E19, 0.0, hugeval, 'Hz', aliases=["breakfreq"])
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.alpha, self.break_, self.norm))


class XSsresc(XSAdditiveModel):
    """The XSPEC sresc model: synchrotron spectrum, cut off by particle escape.

    The model is described at [1]_.

    Attributes
    ----------
    alpha
        The radio spectral index.
    breakfreq
        The break frequency: approximately the frequency at which the
        flux has dropped by a factor of 6 from a straight power law.
        See [1]_ for more details.
    norm
        The 1 Ghz flux in Jy.

    See Also
    --------
    XSsrcut

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSresc.html

    """

    __function__ = "sresc"

    def __init__(self, name='sresc'):
        self.alpha = Parameter(name, 'alpha', 0.5, 0.3, 0.8, 0.0, hugeval)
        self.rolloff = Parameter(name, 'rolloff', 2.42E17, 1.E15, 1.E19, 0.0, hugeval, 'Hz')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.alpha, self.rolloff, self.norm))


@version_at_least("12.10.0")
class XSssa(XSAdditiveModel):
    """The XSPEC ssa model: Strangeon star atmosphere.

    The model is described at [1]_.

    Attributes
    ----------
    te
        The electron temperature, in keV.
    y
        A combination of ion density, temperature, and stellar radius.
        See [1]_ for a decription. The units are 10^42 keV/km/cm^6.
    norm
        This represents (R/d)^2, where d is the distance to the star in
        units of 10 kpc.

    Notes
    -----
    This model is only available when used with XSPEC 12.10.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSsa.html

    """

    __function__ = "ssa"

    def __init__(self, name='ssa'):
        self.te = Parameter(name, 'te', 0.1, 0.01, 0.5, 0.01, 0.5)
        self.y = Parameter(name, 'y', 0.7, 1e-4, 1e3, 1e-4, 1e3)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.te, self.y, self.norm))


class XSstep(XSAdditiveModel):
    """The XSPEC step model: step function convolved with gaussian.

    The model is described at [1]_.

    Attributes
    ----------
    Energy
        The start energy, in keV.
    Sigma
        The gaussian sigma, in keV.
    norm
        The step amplitude.

    See Also
    --------
    XSgaussian

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelStep.html

    """

    __function__ = "xsstep"

    def __init__(self, name='step'):
        self.Energy = Parameter(name, 'Energy', 6.5, 0., 100., 0.0, hugeval, 'keV')
        self.Sigma = Parameter(name, 'Sigma', 0.1, 0., 10., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Energy, self.Sigma, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.Energy, **kwargs)


class XSvapec(XSAdditiveModel):
    """The XSPEC vapec model: APEC emission spectrum.

    The model is described at [1]_, with ``XSapec`` describing how
    it is implemented in Sherpa.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    He, C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element in solar units.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSapec, XSbapec, XSbvapec, XSbvvapec, XSvvapec

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelApec.html

    """

    __function__ = "C_vapec" if equal_or_greater_than("12.9.1") else "xsvape"

    def __init__(self, name='vapec'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0, hugeval, 'keV')
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.norm))


class XSvbremss(XSAdditiveModel):
    """The XSPEC vbremss model: thermal bremsstrahlung.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``HeH`` parameter has been renamed ``HeovrH`` to match the
       XSPEC definition. The name ``HeH`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    kT
        The plasma temperature in keV.
    HeovrH
        The ratio n(He) / n(H).
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbremss, XSzbremss

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBremss.html

    """

    __function__ = "xsbrmv"

    def __init__(self, name='vbremss'):
        self.kT = Parameter(name, 'kT', 3.0, 1.e-2, 100., 0.0, hugeval, 'keV')
        self.HeovrH = Parameter(name, 'HeovrH', 1.0, 0., 100., 0.0, hugeval, aliases=["HeH"])
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kT, self.HeovrH, self.norm))


class XSvequil(XSAdditiveModel):
    """The XSPEC vequil model: collisional plasma, ionization equilibrium.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element in solar units.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelEquil.html

    """

    __function__ = "C_vequil"

    def __init__(self, name='vequil'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.He = Parameter(name, 'He', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.norm))


class XSvgnei(XSAdditiveModel):
    """The XSPEC vgnei model: collisional plasma, non-equilibrium, temperature evolution.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    .. note:: Deprecated in Sherpa 4.10.0

       The ``kT_ave`` parameter has been renamed ``meanKT`` to match the
       XSPEC definition. The name ``kT_ave`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element, with respect to Solar.
    Tau
        The ionization timescale in units of s/cm^3.
    meankT
        The ionization timescale averaged plasma temperature in keV.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSgnei, XSvnei, XSvvgnei

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGnei.html

    """

    __function__ = "C_vgnei"

    def __init__(self, name='vgnei'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.H = Parameter(name, 'H', 1.0, 0., 1., 0.0, hugeval, frozen=True)
        self.He = Parameter(name, 'He', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.meankT = Parameter(name, 'meankT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV', aliases=["kT_ave"])
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.meankT, self.redshift, self.norm))


class XSvvgnei(XSAdditiveModel):
    """The XSPEC vvgnei model: collisional plasma, non-equilibrium, temperature evolution.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    .. note:: Deprecated in Sherpa 4.10.0

       The ``kT_ave`` parameter has been renamed ``meanKT`` to match the
       XSPEC definition. The name ``kT_ave`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element, with respect to Solar.
    Tau
        The ionization timescale in units of s/cm^3.
    meankT
        The ionization timescale averaged plasma temperature in keV.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSgnei, XSvgnei, XSvvnei

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGnei.html

    """

    __function__ = "C_vvgnei"

    def __init__(self, name='vvgnei'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.meankT = Parameter(name, 'meankT', 1.0, 0.0808, 79.9, 0.0808, 79.9, 'keV', aliases=["kT_ave"])
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.meankT, self.redshift, self.norm))


class XSvmeka(XSAdditiveModel):
    """The XSPEC vmeka model: emission, hot diffuse gas (Mewe-Gronenschild).

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    nH
        H density, in cm^-3.
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance relative to Solar.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSmeka

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelMeka.html

    """

    __function__ = "C_vmeka" if equal_or_greater_than("12.9.1") else "xsvmek"

    def __init__(self, name='vmeka'):
        self.kT = Parameter(name, 'kT', 1., 1.e-3, 1.e2, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1., 1.e-5, 1.e19, 0.0, hugeval, 'cm-3', True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.norm))


class XSvmekal(XSAdditiveModel):
    """The XSPEC vmekal model: emission, hot diffuse gas (Mewe-Kaastra-Liedahl).

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    nH
        H density, in cm^-3.
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance relative to Solar.
    redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSapec, XSvmekal

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelMekal.html

    """

    __function__ = "C_vmekal" if equal_or_greater_than("12.9.1") else "xsvmkl"

    def __init__(self, name='vmekal'):
        self.kT = Parameter(name, 'kT', 1., 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1., 1.e-5, 1.e19, 0.0, hugeval, 'cm-3', True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kT, self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na,
                                              self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni,
                                              self.redshift, self.switch, self.norm))


class XSvmcflow(XSAdditiveModel):
    """The XSPEC vmcflow model: cooling flow, mekal.

    The model is described at [1]_. The results of this model depend
    on the cosmology settings set with ``set_xscosmo``.

    Attributes
    ----------
    lowT
        The minimum temperature, in keV.
    highT
        The maxmimum temperature, in keV.
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance relative to Solar.
    redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The mass accretion rate (solar mass per year).

    See Also
    --------
    XSapec, XScflow, XScevmkl, XSmkcflow

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelMkcflow.html

    """

    __function__ = "C_xsvmcf"

    def __init__(self, name='vmcflow'):
        self.lowT = Parameter(name, 'lowT', 0.1, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.highT = Parameter(name, 'highT', 4., 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.lowT, self.highT, self.He, self.C, self.N, self.O, self.Ne, self.Na,
                                              self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni,
                                              self.redshift, self.switch, self.norm))


class XSvnei(XSAdditiveModel):
    """The XSPEC vnei model: collisional plasma, non-equilibrium, constant temperature.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element, with respect to Solar.
    Tau
        The ionization timescale in units of s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSgnei, XSnei, XSvvnei

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNei.html

    """

    __function__ = "C_vnei"

    def __init__(self, name='vnei'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.H = Parameter(name, 'H', 1.0, 0., 1., 0.0, hugeval, frozen=True)
        self.He = Parameter(name, 'He', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.redshift, self.norm))


class XSvvnei(XSAdditiveModel):
    """The XSPEC vvnei model: collisional plasma, non-equilibrium, constant temperature.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element, with respect to Solar.
    Tau
        The ionization timescale in units of s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSgnei, XSnei, XSvnei

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNei.html

    """

    __function__ = "C_vvnei"

    def __init__(self, name='vvnei'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.redshift, self.norm))


class XSvnpshock(XSAdditiveModel):
    """The XSPEC vnpshock model: shocked plasma, plane parallel, separate ion, electron temperatures.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT_a
        The mean shock temperature, in keV.
    kT_b
        The electron temperature immediately behind the shock
        front, in keV. See [1]_ for a discussion of the behavior
        of kT_a and kT_b.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element, with respect to Solar.
    Tau_l
        The lower limit on the ionization timescale in s/cm^3.
    Tau_u
        The upper limit on the ionization timescale in s/cm^3.
    redshift
        The redshift of the source.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSnpshock, XSvvnpshock

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNpshock.html

    """

    __function__ = "C_vnpshock"

    def __init__(self, name='vnpshock'):
        self.kT_a = Parameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.kT_b = Parameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0, hugeval, 'keV')
        self.H = Parameter(name, 'H', 1.0, 0., 1., 0.0, hugeval, frozen=True)
        self.He = Parameter(name, 'He', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau_l = Parameter(name, 'Tau_l', 0.0, 0., 5.e13, 0.0, hugeval, 's/cm^3', True)
        self.Tau_u = Parameter(name, 'Tau_u', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau_l, self.Tau_u, self.redshift, self.norm))


class XSvvnpshock(XSAdditiveModel):
    """The XSPEC vvnpshock model: shocked plasma, plane parallel, separate ion, electron temperatures.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT_a
        The mean shock temperature, in keV.
    kT_b
        The electron temperature immediately behind the shock
        front, in keV. See [1]_ for a discussion of the behavior
        of kT_a and kT_b.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element, with respect to Solar.
    Tau_l
        The lower limit on the ionization timescale in s/cm^3.
    Tau_u
        The upper limit on the ionization timescale in s/cm^3.
    redshift
        The redshift of the source.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSnpshock, XSvnpshock

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNpshock.html

    """

    __function__ = "C_vvnpshock"

    def __init__(self, name='vvnpshock'):
        self.kT_a = Parameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_b = Parameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0100, 79.9, units='keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Tau_l = Parameter(name, 'Tau_l', 0.0, 0.0, 5.0e13, 0.0, 5.0e13, units='s/cm^3', frozen=True)
        self.Tau_u = Parameter(name, 'Tau_u', 1.0e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau_l, self.Tau_u, self.redshift, self.norm))


class XSvpshock(XSAdditiveModel):
    """The XSPEC vpshock model: plane-parallel shocked plasma, constant temperature.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element in solar units.
    Tau_l
        The lower limit on the ionization timescale, in s/cm^3.
    Tau_u
        The upper limit on the ionization timescale, in s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSpshock, XSvvpshock

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPshock.html

    """

    __function__ = "C_vpshock"

    def __init__(self, name='vpshock'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.H = Parameter(name, 'H', 1.0, 0., 1., 0.0, hugeval, frozen=True)
        self.He = Parameter(name, 'He', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau_l = Parameter(name, 'Tau_l', 0.0, 0., 5.e13, 0.0, hugeval, 's/cm^3', True)
        self.Tau_u = Parameter(name, 'Tau_u', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau_l, self.Tau_u, self.redshift, self.norm))


class XSvvpshock(XSAdditiveModel):
    """The XSPEC vvpshock model: plane-parallel shocked plasma, constant temperature.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element, with respect to Solar.
    Tau_l
        The lower limit on the ionization timescale, in s/cm^3.
    Tau_u
        The upper limit on the ionization timescale, in s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSequil, XSpshock, XSvpshock

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPshock.html

    """

    __function__ = "C_vvpshock"

    def __init__(self, name='vvpshock'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Tau_l = Parameter(name, 'Tau_l', 0.0, 0.0, 5.0e13, 0.0, 5.0e13, units='s/cm^3', frozen=True)
        self.Tau_u = Parameter(name, 'Tau_u', 1.0e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10.0, -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau_l, self.Tau_u, self.redshift, self.norm))


class XSvraymond(XSAdditiveModel):
    """The XSPEC vraymond model: emission, hot diffuse gas, Raymond-Smith.

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni
        The abundance relative to Solar.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSraymond

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRaymond.html

    """

    __function__ = "C_vraysmith" if equal_or_greater_than("12.9.1") else "xsvrys"

    def __init__(self, name='vraymond'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.He = Parameter(name, 'He', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.norm))


class XSvsedov(XSAdditiveModel):
    """The XSPEC vsedov model: sedov model, separate ion/electron temperature.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT_a
        The mean shock temperature, in keV.
    kT_b
        The electron temperature immediately behind the shock
        front, in keV. See [1]_ for a discussion of the behavior
        of kT_a and kT_b.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element, with respect to Solar.
    Tau
        The ionization timescale in units of s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSsedov, XSvvsedov

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSedov.html

    """

    __function__ = "C_vsedov"

    def __init__(self, name='vsedov'):
        self.kT_a = Parameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.kT_b = Parameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0, hugeval, 'keV')
        self.H = Parameter(name, 'H', 1.0, 0., 1., 0.0, hugeval, frozen=True)
        self.He = Parameter(name, 'He', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.redshift, self.norm))


class XSvvsedov(XSAdditiveModel):
    """The XSPEC vvsedov model: sedov model, separate ion/electron temperature.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT_a
        The mean shock temperature, in keV.
    kT_b
        The electron temperature immediately behind the shock
        front, in keV. See [1]_ for a discussion of the behavior
        of kT_a and kT_b.
    H
        The H abundance: it should be set to 0 to switch on and
        1 to switch off the free-free continuum.
    He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element, with respect to Solar.
    Tau
        The ionization timescale in units of s/cm^3.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSsedov, XSvsedov

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSedov.html

    """

    __function__ = "C_vvsedov"

    def __init__(self, name='vvsedov'):
        self.kT_a = Parameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0808, hugeval, 'keV')
        self.kT_b = Parameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.01, 79.9, 'keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, 10000.0, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.redshift, self.norm))


class XSzbbody(XSAdditiveModel):
    """The XSPEC zbbody model: blackbody spectrum.

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The temperature of the object, in keV.
    redshift
        The redshift of the object.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbbody, XSbbodyrad

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBbody.html

    """

    __function__ = "xszbod"

    def __init__(self, name='zbbody'):
        self.kT = Parameter(name, 'kT', 3.0, 1.e-2, 100., 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.redshift, self.norm))


@version_at_least("12.10.1")
class XSzbknpower(XSAdditiveModel):
    """The XSPEC bknpower model: broken power law.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndx1
        The power law photon index for energies less than BreakE.
    BreakE
        The break energy, in keV.
    PhoIndx2
        The power law photon index for energies greater than BreakE.
    Redshift
        The redshift.
    norm
        The normalization of the model. See [1]_ for details, as its
        meaning depends on whether the "POW_EMIN" or "POW_EMAX"
        keywords have been set with ``set_xsxset``.

    See Also
    --------
    XSbknpower

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBknpower.html

    """

    __function__ = "C_zBrokenPowerLaw"

    def __init__(self, name='zbknpower'):
        self.PhoIndx1 = Parameter(name, 'PhoIndx1', 1., -2., 9., -hugeval, hugeval)
        self.BreakE = Parameter(name, 'BreakE', 5., 1.e-2, 1.e6, 0.0, hugeval, 'keV')
        self.PhoIndx2 = Parameter(name, 'PhoIndx2', 2., -2., 9., -hugeval, hugeval)
        self.Redshift = Parameter(name, 'Redshift', 0, -0.999, 10, -0.999, 10,
                                  frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.PhoIndx1, self.BreakE, self.PhoIndx2, self.Redshift,
                self.norm)

        XSAdditiveModel.__init__(self, name, pars)

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.BreakE, **kwargs)


class XSzbremss(XSAdditiveModel):
    """The XSPEC zbremss model: thermal bremsstrahlung.

    The model is described at [1]_.

    Attributes
    ----------
    kT
        The plasma temperature in keV.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbremss, XSvbremss

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBremss.html

    """

    __function__ = "xszbrm"

    def __init__(self, name='zbremss'):
        self.kT = Parameter(name, 'kT', 7.0, 1.e-4, 100., 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.redshift, self.norm))


@version_at_least("12.10.0")
class XSzcutoffpl(XSAdditiveModel):
    """The XSPEC zcutoffpl model: power law, high energy exponential cutoff.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndex
        The power law photon index.
    HighECut
        The e-folding energy of the exponential rolloff, in keV.
    Redshift
        The redshift of the source.
    norm
        The normalization of the model. See [1]_ for details, as its
        meaning depends on whether the "POW_EMIN" or "POW_EMAX"
        keywords have been set with ``set_xsxset``.

    See Also
    --------
    XSbknpower, XSbkn2pow, XSzcutoffpl, XSpowerlaw, XSzpowerlw

    Notes
    -----
    This model is only available when used with XSPEC 12.10.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCutoffpl.html

    """

    __function__ = "C_zcutoffPowerLaw"

    def __init__(self, name='zcutoffpl'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 1., -2., 9., -hugeval, hugeval)
        self.HighECut = Parameter(name, 'HighECut', 15., 1., 500., 0.0, hugeval, 'keV')
        self.Redshift = Parameter(name, 'Redshift', 0.0, -0.999, 10., -0.999, 10., frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.HighECut, self.Redshift, self.norm))


class XSzgauss(XSAdditiveModel):
    """The XSPEC zgauss model: gaussian line profile.

    The model is described at [1]_.

    Attributes
    ----------
    LineE
        The line energy, in keV.
    Sigma
        The line width, in keV. A value of zero means a delta function.
    redshift
        The redshift of the line.
    norm
        The flux in the line, in units of photon/cm^2/s.

    See Also
    --------
    XSagauss, XSgaussian, XSzagauss

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGaussian.html

    """

    __function__ = "C_xszgau"

    def __init__(self, name='zgauss'):
        self.LineE = Parameter(name, 'LineE', 6.5, 0., 1.e6, 0.0, hugeval, 'keV')
        self.Sigma = Parameter(name, 'Sigma', 0.1, 0., 10., 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Sigma, self.redshift, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


@version_at_least("12.10.1")
class XSzlogpar(XSAdditiveModel):
    """The XSPEC zlogpar model: log-parabolic blazar model.

    The model is described at [1]_.

    Attributes
    ----------
    alpha
        The slope at the pivot energy.
    beta
        The curvature term.
    pivotE
        The pivot energy, in keV.
    Redshift
        The redshift.
    norm
        The normalization of the model: see [1]_ for more details.

    See Also
    --------
    XSlogpar

    Notes
    -----
    This model is only available when used with XSPEC 12.10.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLogpar.html

    """

    __function__ = "C_zLogpar"

    def __init__(self, name='zlogpar'):
        self.alpha = Parameter(name, 'alpha', 1.5, 0., 4., 0.0, hugeval)
        self.beta = Parameter(name, 'beta', 0.2, -4., 4., -hugeval, hugeval)
        self.pivotE = Parameter(name, 'pivotE', 1.0, units='keV',
                                alwaysfrozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0, -0.999, 10, -0.999, 10,
                                  frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.alpha, self.beta, self.pivotE, self.Redshift, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


class XSzpowerlw(XSAdditiveModel):
    """The XSPEC zpowerlw model: redshifted power law photon spectrum.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndex
        The power law photon index.
    redshift
        The redshift.
    norm
        The normalization of the model. See [1]_ for details, as its
        meaning depends on whether the "POW_EMIN" or "POW_EMAX"
        keywords have been set with ``set_xsxset``.

    See Also
    --------
    XSbknpower, XSbkn2pow, XScutoffpl, XSpowerlaw, XSzcutoffpl

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPowerlaw.html

    """

    __function__ = "C_zpowerLaw"

    def __init__(self, name='zpowerlw'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 1., -2., 9., -hugeval, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.redshift, self.norm))


class XSabsori(XSMultiplicativeModel):
    """The XSPEC absori model: ionized absorber.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndex
        The power-law photon index.
    nH
        The Hydrogen column, in units of 10^22 cm^-2.
    Temp_abs
        The absorber temperature, in K.
    xi
        The absorber ionization state. See [1]_ for more details.
    redshift
        The redshift of the absorber.
    Fe_abund
        The iron abundance with respect to solar, as set by the
        ``set_xsabund`` function.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelAbsori.html

    """

    __function__ = "C_xsabsori"

    def __init__(self, name='absori'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 2., 0., 4., 0.0, hugeval, frozen=True)
        self.nH = Parameter(name, 'nH', 1., 0., 100., 0.0, hugeval, '10^22 atoms / cm^2')
        self.Temp_abs = Parameter(name, 'Temp_abs', 3.e4, 1.e4, 1.e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 1., 0., 1.e3, 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0., 1.e6, 0.0, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.PhoIndex, self.nH, self.Temp_abs, self.xi, self.redshift, self.Fe_abund))


class XSacisabs(XSMultiplicativeModel):
    """The XSPEC acisabs model: Chandra ACIS q.e. decay.

    The model is described at [1]_.

    Attributes
    ----------
    Tdays
        Days between Chandra launch and ACIS observation.
    norm
    tauinf
        Slope of linear quantum efficiency decay.
    tefold
        Offset of linear quantum efficiency decay.
    nC
        Number of carbon atoms in hydrocarbon.
    nH
        Number of hydrogen atoms in hydrocarbon.
    nO
        Number of oxygen atoms in hydrocarbon.
    nN
        Number of nitrogen atoms in hydrocarbon.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelAcisabs.html

    """

    __function__ = _f77_or_c_12100("acisabs")

    def __init__(self, name='acisabs'):
        self.Tdays = Parameter(name, 'Tdays', 850., 0., 10000., 0.0, hugeval, 'days', True)
        self.norm = Parameter(name, 'norm', 0.00722, 0., 1., 0.0, hugeval, frozen=True)
        self.tauinf = Parameter(name, 'tauinf', 0.582, 0., 1., 0.0, hugeval, frozen=True)
        self.tefold = Parameter(name, 'tefold', 620., 1., 10000., 0.0, hugeval, 'days', True)
        self.nC = Parameter(name, 'nC', 10., 0., 50., 0.0, hugeval, frozen=True)
        self.nH = Parameter(name, 'nH', 20., 1., 50., 0.0, hugeval, frozen=True)
        self.nO = Parameter(name, 'nO', 2., 0., 50., 0.0, hugeval, frozen=True)
        self.nN = Parameter(name, 'nN', 1., 0., 50., 0.0, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.Tdays, self.norm, self.tauinf, self.tefold, self.nC, self.nH, self.nO, self.nN))


class XSconstant(XSMultiplicativeModel):
    """The XSPEC constant model: energy-independent factor.

    The model is described at [1]_.

    Attributes
    ----------
    factor
        The value of the model.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelConstant.html

    """

    __function__ = "xscnst"

    def __init__(self, name='constant'):
        self.factor = Parameter(name, 'factor', 1., 0.0, 1.e10, 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.factor,))


class XScabs(XSMultiplicativeModel):
    """The XSPEC cabs model: Optically-thin Compton scattering.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The Hydrogen column, in units of 10^22 cm^-2.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCabs.html

    """

    __function__ = "xscabs"

    def __init__(self, name='cabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        XSMultiplicativeModel.__init__(self, name, (self.nH,))


class XScyclabs(XSMultiplicativeModel):
    """The XSPEC cyclabs model: absorption line, cyclotron.

    The model is described at [1]_.

    Attributes
    ----------
    Depth0
        The depth of the fundamental.
    E0
        The cyclotron energy, in keV.
    Width0
        The width of the fundamental, in keV.
    Depth2
        The depth of the second harmonic.
    Width2
        The width of the second harmonic, in keV.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCyclabs.html

    """

    __function__ = "xscycl"

    def __init__(self, name='cyclabs'):
        self.Depth0 = Parameter(name, 'Depth0', 2.0, 0., 100., 0.0, hugeval)
        self.E0 = Parameter(name, 'E0', 30.0, 1.0, 100., 0.0, hugeval, 'keV')
        self.Width0 = Parameter(name, 'Width0', 10.0, 1.0, 100., 0.0, hugeval, 'keV', True)
        self.Depth2 = Parameter(name, 'Depth2', 0.0, 0., 100., 0.0, hugeval, frozen=True)
        self.Width2 = Parameter(name, 'Width2', 20.0, 1.0, 100., 0.0, hugeval, 'keV', True)
        XSMultiplicativeModel.__init__(self, name, (self.Depth0, self.E0, self.Width0, self.Depth2, self.Width2))


class XSdust(XSMultiplicativeModel):
    """The XSPEC dust model: dust scattering.

    The model is described at [1]_.

    Attributes
    ----------
    Frac
        The scattering fraction at 1 keV.
    Halosz
        The size of the halo at 1 keV in units of the detector beamsize.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDust.html

    """

    __function__ = "xsdust"

    def __init__(self, name='dust'):
        self.Frac = Parameter(name, 'Frac', 0.066, 0., 1., 0.0, hugeval, frozen=True)
        self.Halosz = Parameter(name, 'Halosz', 2., 0., 1.e5, 0.0, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.Frac, self.Halosz))


class XSedge(XSMultiplicativeModel):
    """The XSPEC edge model: absorption edge.

    The model is described at [1]_.

    Attributes
    ----------
    edgeE
        The threshold edge, in keV.
    MaxTau
        The absorption depth at the threshold.

    See Also
    --------
    XSsmedge, XSzedge

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelEdge.html

    """

    __function__ = "xsedge"

    def __init__(self, name='edge'):
        self.edgeE = Parameter(name, 'edgeE', 7.0, 0., 100., 0.0, hugeval, 'keV')
        self.MaxTau = Parameter(name, 'MaxTau', 1., 0., 5., 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.edgeE, self.MaxTau))


class XSexpabs(XSMultiplicativeModel):
    """The XSPEC expabs model: exponential roll-off at low E.

    The model is described at [1]_.

    Attributes
    ----------
    LowECut
        The e-folding energy for the absorption, in keV.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelExpabs.html

    """

    __function__ = "xsabsc"

    def __init__(self, name='expabs'):
        self.LowECut = Parameter(name, 'LowECut', 2., 0., 100., 0.0, hugeval, 'keV')
        XSMultiplicativeModel.__init__(self, name, (self.LowECut,))


class XSexpfac(XSMultiplicativeModel):
    """The XSPEC expfac model: exponential modification.

    The model is described at [1]_.

    Attributes
    ----------
    Ampl
        The amplitude of the effect.
    Factor
        The exponential factor.
    StartE
        The start energy of the modification, in keV.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelExpfac.html

    """

    __function__ = "xsexp"

    def __init__(self, name='expfac'):
        self.Ampl = Parameter(name, 'Ampl', 1., 0., 1.e5, 0.0, hugeval)
        self.Factor = Parameter(name, 'Factor', 1., 0., 1.e5, 0.0, hugeval)
        self.StartE = Parameter(name, 'StartE', 0.5, 0., 1.e5, 0.0, hugeval, 'keV', True)
        XSMultiplicativeModel.__init__(self, name, (self.Ampl, self.Factor, self.StartE))


class XSgabs(XSMultiplicativeModel):
    """The XSPEC gabs model: gaussian absorption line.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``Tau`` parameter has been renamed ``Strength`` to match the
       XSPEC definition. The name ``Tau`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    LineE
        The line energy, in keV.
    Sigma
        The line width (sigma), in keV.
    Strength
        The line depth. The optical depth at the line center is
        Strength / (sqrt(2 pi) * Sigma).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGabs.html

    """

    __function__ = "C_gaussianAbsorptionLine"

    def __init__(self, name='gabs'):
        self.LineE = Parameter(name, 'LineE', 1.0, 0., 1.e6, 0.0, hugeval, 'keV')
        self.Sigma = Parameter(name, 'Sigma', 0.01, 0., 10., 0.0, hugeval, 'keV')
        self.Strength = Parameter(name, 'Strength', 1.0, 0., 1.e6, 0.0, hugeval, aliases=["Tau"])

        XSMultiplicativeModel.__init__(self, name, (self.LineE, self.Sigma, self.Strength))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


class XShighecut(XSMultiplicativeModel):
    """The XSPEC highecut model: high-energy cutoff.

    The model is described at [1]_.

    Attributes
    ----------
    cutoffE
        The cut-off energy, in keV.
    foldE
        The e-folding energy, in keV.

    See Also
    --------
    XSzhighect

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelHighecut.html

    """

    __function__ = "xshecu"

    def __init__(self, name='highecut'):
        self.cutoffE = Parameter(name, 'cutoffE', 10., 1.e-2, 1.e6, 0.0, hugeval, 'keV')
        self.foldE = Parameter(name, 'foldE', 15., 1.e-2, 1.e6, 0.0, hugeval, 'keV')
        XSMultiplicativeModel.__init__(self, name, (self.cutoffE, self.foldE))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.cutoffE, **kwargs)


class XShrefl(XSMultiplicativeModel):
    """The XSPEC hrefl model: reflection model.

    The model is described at [1]_.

    Attributes
    ----------
    thetamin
        The minimum angle, in degrees, between source photons incident
        on the slab and the slab normal.
    thetamax
        The maximum angle, in degrees, between source photons incident
        on the slab and the slab normal.
    thetamobs
        The angle, in degrees, between the observer's line of sight and
        the slab normal.
    Feabun
        The iron abundance relative to solar.
    FeKedge
        The iron K edge energy, in keV.
    Escfrac
        The fraction of the direct flux seen by the observer: see [1]_
        for more details.
    covfac
        The normalization of the reflected continuum: see [1]_ for more
        details.
    redshift
        The redshift.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelHrefl.html

    """

    __function__ = "xshrfl"

    def __init__(self, name='hrefl'):
        self.thetamin = Parameter(name, 'thetamin', 0., 0.0, 90., 0.0, hugeval, frozen=True)
        self.thetamax = Parameter(name, 'thetamax', 90., 0.0, 90., 0.0, hugeval, frozen=True)
        self.thetaobs = Parameter(name, 'thetaobs', 60., 0.0, 90., 0.0, hugeval)
        self.Feabun = Parameter(name, 'Feabun', 1., 0.0, 100., 0.0, hugeval, frozen=True)
        self.FeKedge = Parameter(name, 'FeKedge', 7.11, 7.0, 10., 0.0, hugeval, 'keV', True)
        self.Escfrac = Parameter(name, 'Escfrac', 1.0, 0.0, 500., 0.0, hugeval)
        self.covfac = Parameter(name, 'covfac', 1.0, 0.0, 500., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.thetamin, self.thetamax, self.thetaobs, self.Feabun, self.FeKedge, self.Escfrac, self.covfac, self.redshift))


class XSnotch(XSMultiplicativeModel):
    """The XSPEC notch model: absorption line, notch.

    The model is described at [1]_.

    Attributes
    ----------
    LineE
        The line energy, in keV.
    Width
        The line width, in keV.
    CvrFract
        The covering fraction.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNotch.html

    """

    __function__ = "xsntch"

    def __init__(self, name='notch'):
        self.LineE = Parameter(name, 'LineE', 3.5, 0., 20., 0.0, hugeval, 'keV')
        self.Width = Parameter(name, 'Width', 1., 0., 20., 0.0, hugeval, 'keV')
        self.CvrFract = Parameter(name, 'CvrFract', 1., 0., 1., 0.0, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.LineE, self.Width, self.CvrFract))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


class XSpcfabs(XSMultiplicativeModel):
    """The XSPEC pcfabs model: partial covering fraction absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    CvrFract
        The covering fraction.

    See Also
    --------
    XSphabs, XSpwab, XSzpcfabs

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPcfabs.html

    """

    __function__ = "xsabsp"

    def __init__(self, name='pcfabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.CvrFract = Parameter(name, 'CvrFract', 0.5, 0.05, 0.95, 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.CvrFract))


class XSphabs(XSMultiplicativeModel):
    """The XSPEC phabs model: photoelectric absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).

    See Also
    --------
    XSvphabs, XSzphabs, XSzvphabs

    Notes
    -----
    The `set_xsxsect` function changes the cross sections used by this
    model. The `set_xsabund` function changes the relative abundances of
    the elements.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPhabs.html

    """

    __function__ = "xsphab"

    def __init__(self, name='phabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        XSMultiplicativeModel.__init__(self, name, (self.nH,))


class XSplabs(XSMultiplicativeModel):
    """The XSPEC plabs model: power law absorption.

    The model is described at [1]_.

    Attributes
    ----------
    index
        The power-law index.
    coef
        The normalization.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPlabs.html

    """

    __function__ = "xsplab"

    def __init__(self, name='plabs'):
        self.index = Parameter(name, 'index', 2.0, 0.0, 5., 0.0, hugeval)
        self.coef = Parameter(name, 'coef', 1.0, 0.0, 100., 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.index, self.coef))


class XSpwab(XSMultiplicativeModel):
    """The XSPEC pwab model: power-law distribution of neutral absorbers.

    The model is described at [1]_.

    Attributes
    ----------
    nHmin
        The minimum equivalent hydrogen column (in units of
        10^22 atoms/cm^2).
    nHmax
        The maximum equivalent hydrogen column (in units of
        10^22 atoms/cm^2).
    beta
        The power law index for the covering fraction.

    See Also
    --------
    XSpcfabs, XSwabs

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPwab.html

    """

    __function__ = "C_xspwab"

    def __init__(self, name='pwab'):
        self.nHmin = Parameter(name, 'nHmin', 1., 1.e-7, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.nHmax = Parameter(name, 'nHmax', 2., 1.e-7, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.beta = Parameter(name, 'beta', 1.0, -10., 10, -hugeval, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nHmin, self.nHmax, self.beta))


class XSredden(XSMultiplicativeModel):
    """The XSPEC redden model: interstellar extinction.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``EBV`` parameter has been renamed ``E_BmV`` to match the
       XSPEC definition. The name ``EBV`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    E_BmV
        The value of E(B-v) for the line of sight to the source.

    See Also
    --------
    XSzredden

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRedden.html

    """

    __function__ = "xscred"

    def __init__(self, name='redden'):
        self.E_BmV = Parameter(name, 'E_BmV', 0.05, 0., 10., 0.0, hugeval, aliases=["EBV"])

        XSMultiplicativeModel.__init__(self, name, (self.E_BmV,))


class XSsmedge(XSMultiplicativeModel):
    """The XSPEC smedge model: smeared edge.

    The model is described at [1]_.

    Attributes
    ----------
    edgeE
        The threshold edge, in keV.
    MaxTau
        The maximum absorption factor at the threshold.
    index
        The index for photo-electric cross-section.
    width
        The smearing width, in keV.

    See Also
    --------
    XSedge, XSzedge

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSmedge.html

    """

    __function__ = "xssmdg"

    def __init__(self, name='smedge'):
        self.edgeE = Parameter(name, 'edgeE', 7.0, 0.1, 100., 0.0, hugeval, 'keV')
        self.MaxTau = Parameter(name, 'MaxTau', 1., 0., 5., 0.0, hugeval)
        self.index = Parameter(name, 'index', -2.67, -10., 10., -hugeval, hugeval, frozen=True)
        self.width = Parameter(name, 'width', 10., 0.01, 100., 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.edgeE, self.MaxTau, self.index, self.width))


class XSspexpcut(XSMultiplicativeModel):
    """The XSPEC spexpcut model: super-exponential cutoff absorption.

    The model is described at [1]_.

    Attributes
    ----------
    Ecut
        The e-folding energy for the absorption, in keV.
    alpha
        The exponent index: see [1]_ for more details.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSpexpcut.html

    """

    __function__ = "C_superExpCutoff"

    def __init__(self, name='spexpcut'):
        self.Ecut = Parameter(name, 'Ecut', 10.0, 0.0, 1e6, 0.0, hugeval, 'keV')
        self.alpha = Parameter(name, 'alpha', 1.0, -5.0, 5.0, -hugeval, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.Ecut, self.alpha))


class XSspline(XSMultiplicativeModel):
    """The XSPEC spline model: spline modification.

    The model is described at [1]_.

    Attributes
    ----------
    Estart
        The start x value (energy), in keV.
    YStart
        The start y value.
    Yend
        The end y value.
    YPstart
        The start dy/dx value.
    YPEnd
        The end dy/dx value.
    Eend
        The end x value (energy), in keV.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSpline.html

    """

    __function__ = "xsspln"

    def __init__(self, name='spline'):
        self.Estart = Parameter(name, 'Estart', 0.1, 0., 100., 0.0, hugeval, 'keV')
        self.Ystart = Parameter(name, 'Ystart', 1., -1.e6, 1.e6, -hugeval, hugeval)
        self.Yend = Parameter(name, 'Yend', 1., -1.e6, 1.e6, -hugeval, hugeval)
        self.YPstart = Parameter(name, 'YPstart', 0., -1.e6, 1.e6, -hugeval, hugeval)
        self.YPend = Parameter(name, 'YPend', 0., -1.e6, 1.e6, -hugeval, hugeval)
        self.Eend = Parameter(name, 'Eend', 15., 0., 100., 0.0, hugeval, 'keV')
        XSMultiplicativeModel.__init__(self, name, (self.Estart, self.Ystart, self.Yend, self.YPstart, self.YPend, self.Eend))


class XSSSS_ice(XSMultiplicativeModel):
    """The XSPEC sss_ice model: Einstein SSS ice absorption.

    The model is described at [1]_.

    Attributes
    ----------
    clumps
        The ice thickness parameter.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSssice.html

    """

    __function__ = "xssssi"

    def __init__(self, name='sss_ice'):
        self.clumps = Parameter(name, 'clumps', 0.0, 0., 10., 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.clumps,))


class XSswind1(XSMultiplicativeModel):
    """The XSPEC swind1 model: absorption by partially ionized material with large velocity shear.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``logxi`` parameter has been renamed ``log_xi`` to match the
       XSPEC definition. The name ``logxi`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    column
        The column density, in units of 10^22 cm^2.
    log_xi
        The log of xi: see [1]_ for more details.
    sigma
        The gaussian sigma for velocity smearing (v/c).
    redshift
        The redshift of the source.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSwind1.html

    """

    __function__ = _f77_or_c_12100("swind1")

    def __init__(self, name='swind1'):
        self.column = Parameter(name, 'column', 6., 3., 50., 0.0, hugeval)
        self.log_xi = Parameter(name, 'log_xi', 2.5, 2.1, 4.1, 0.0, hugeval, aliases=["logxi"])
        self.sigma = Parameter(name, 'sigma', 0.1, 0., .5, 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)

        XSMultiplicativeModel.__init__(self, name, (self.column, self.log_xi, self.sigma, self.redshift))


class XSTBabs(XSMultiplicativeModel):
    """The XSPEC TBabs model: ISM grain absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).

    See Also
    --------
    XSTBfeo, XSTBgas, XSTBgrain, XSTBpcf, XSTBrel, XSTBvarabs, XSzTBabs

    Notes
    -----
    The `set_xsabund` function changes the relative abundances of
    the elements, in particular the "wilm" setting.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbabs"

    def __init__(self, name='tbabs'):
        self.nH = Parameter(name, 'nH', 1., 0., 1E5, 0.0, hugeval, '10^22 atoms / cm^2')
        XSMultiplicativeModel.__init__(self, name, (self.nH,))


class XSTBgrain(XSMultiplicativeModel):
    """The XSPEC TBgrain model: ISM grain absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    h2
        The equivalent molecular hydrogen column (in units of
         10^22 atoms/cm^2).
    rho
        The grain density, in g/cm^3.
    amin
        The minimum grain size, in micro-meters.
    amax
        The maximum grain size, in micro-meters.
    PL
        The power-law index of grain sizes.

    See Also
    --------
    XSTBabs, XSTBfeo, XSTBgas, XSTBpcf, XSTBrel, XSTBvarabs, XSzTBabs

    Notes
    -----
    The `set_xsabund` function changes the relative abundances of
    the elements, in particular the "wilm" setting.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbgrain"

    def __init__(self, name='tbgrain'):
        self.nH = Parameter(name, 'nH', 1., 0., 1E5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.h2 = Parameter(name, 'h2', 0.2, 0., 1., 0.0, hugeval, frozen=True)
        self.rho = Parameter(name, 'rho', 1., 0., 5., 0.0, hugeval, 'g/cm^3', True)
        self.amin = Parameter(name, 'amin', 0.025, 0., 0.25, 0.0, hugeval, 'mum', True)
        self.amax = Parameter(name, 'amax', 0.25, 0., 1., 0.0, hugeval, 'mum', True)
        self.PL = Parameter(name, 'PL', 3.5, 0., 5., 0.0, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.h2, self.rho, self.amin, self.amax, self.PL))


class XSTBvarabs(XSMultiplicativeModel):
    """The XSPEC TBvarabs model: ISM grain absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Cl, Ar, Ca, Cr, Fe, Co, Ni
        The abundance of the element in solar units.
    H2
        The equivalent molecular hydrogen column (in units of
         10^22 atoms/cm^2).
    rho
        The grain density, in g/cm^3.
    amin
        The minimum grain size, in micro-meters.
    amax
        The maximum grain size, in micro-meters.
    PL
        The power-law index of grain sizes.
    H_dep, He_dep, C_dep, N_dep, O_dep, Ne_dep, Na_dep, Mg_dep, Al_dep,
    Si_dep, S_dep, Cl_dep, Ar_dep, Ca_dep, Cr_dep, Fe_dep, Co_dep, Ni_dep
        The grain depletion fraction of the element.
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSTBabs, XSTBfeo, XSTBgas, XSTBgrain, XSTBpcf, XSTBrel, XSzTBabs

    Notes
    -----
    The `set_xsabund` function changes the relative abundances of
    the elements, in particular the "wilm" setting.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbvabs"

    def __init__(self, name='tbvarabs'):
        self.nH = Parameter(name, 'nH', 1., 0., 1E5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.He = Parameter(name, 'He', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.H2 = Parameter(name, 'H2', 0.2, 0., 1., 0.0, hugeval, frozen=True)
        self.rho = Parameter(name, 'rho', 1., 0., 5., 0.0, hugeval, 'g/cm^3', True)
        self.amin = Parameter(name, 'amin', 0.025, 0., 0.25, 0.0, hugeval, 'mum', True)
        self.amax = Parameter(name, 'amax', 0.25, 0., 1., 0.0, hugeval, 'mum', True)
        self.PL = Parameter(name, 'PL', 3.5, 0., 5., 0.0, hugeval, frozen=True)
        self.H_dep = Parameter(name, 'H_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.He_dep = Parameter(name, 'He_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.C_dep = Parameter(name, 'C_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.N_dep = Parameter(name, 'N_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.O_dep = Parameter(name, 'O_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Ne_dep = Parameter(name, 'Ne_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Na_dep = Parameter(name, 'Na_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Mg_dep = Parameter(name, 'Mg_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Al_dep = Parameter(name, 'Al_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Si_dep = Parameter(name, 'Si_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.S_dep = Parameter(name, 'S_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Cl_dep = Parameter(name, 'Cl_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Ar_dep = Parameter(name, 'Ar_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Ca_dep = Parameter(name, 'Ca_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Cr_dep = Parameter(name, 'Cr_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Fe_dep = Parameter(name, 'Fe_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Co_dep = Parameter(name, 'Co_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.Ni_dep = Parameter(name, 'Ni_dep', 1., 0., 1., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Cl, self.Ar, self.Ca, self.Cr, self.Fe, self.Co, self.Ni, self.H2, self.rho, self.amin, self.amax, self.PL, self.H_dep, self.He_dep, self.C_dep, self.N_dep, self.O_dep, self.Ne_dep, self.Na_dep, self.Mg_dep, self.Al_dep, self.Si_dep, self.S_dep, self.Cl_dep, self.Ar_dep, self.Ca_dep, self.Cr_dep, self.Fe_dep, self.Co_dep, self.Ni_dep, self.redshift))


class XSuvred(XSMultiplicativeModel):
    """The XSPEC uvred model: interstellar extinction, Seaton Law.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``EBV`` parameter has been renamed ``E_BmV`` to match the
       XSPEC definition. The name ``EBV`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    E_BmV
        The value of E(B-v) for the line of sight to the source.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelUvred.html

    """

    __function__ = "xsred"

    def __init__(self, name='uvred'):
        self.E_BmV = Parameter(name, 'E_BmV', 0.05, 0., 10., 0.0, hugeval, aliases=["EBV"])

        XSMultiplicativeModel.__init__(self, name, (self.E_BmV,))


class XSvarabs(XSMultiplicativeModel):
    """The XSPEC varabs model: photoelectric absorption.

    The model is described at [1]_.

    Attributes
    ----------
    H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Cl, Ar, Ca, Cr, Fe, Co, Ni
        See [1]_ for a description of the units.

    See Also
    --------
    XSvphabs, XSzvarabs

    Notes
    -----
    The `set_xsxsect` function changes the cross sections used by this
    model. The `set_xsabund` function changes the relative abundances of
    the elements.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelVarabs.html

    """

    __function__ = "xsabsv"

    def __init__(self, name='varabs'):
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, hugeval, 'sH22', True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, 'sHe22', True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, 'sC22', True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, 'sN22', True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, 'sO22', True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, 'sNe22', True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, 'sNa22', True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, 'sMg22', True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, 'sAl22', True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, 'sSi22', True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, 'sS22', True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, hugeval, 'sCl22', True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, 'sAr22', True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, 'sCa22', True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, hugeval, 'sCr22', True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, 'sFe22', True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, hugeval, 'sCo22', True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, 'sNi22', True)
        XSMultiplicativeModel.__init__(self, name, (self.H, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Cl, self.Ar, self.Ca, self.Cr, self.Fe, self.Co, self.Ni))


class XSvphabs(XSMultiplicativeModel):
    """The XSPEC vphabs model: photoelectric absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Cl, Ar, Ca, Cr, Fe, Co, Ni
        The abundance of the element in solar units.

    See Also
    --------
    XSphabs, XSvarabs, XSzphabs, XSzvphabs

    Notes
    -----
    The `set_xsxsect` function changes the cross sections used by this
    model. The `set_xsabund` function changes the relative abundances of
    the elements.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPhabs.html

    """

    __function__ = "xsvphb"

    def __init__(self, name='vphabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Cl, self.Ar, self.Ca, self.Cr, self.Fe, self.Co, self.Ni))


class XSwabs(XSMultiplicativeModel):
    """The XSPEC wabs model: photoelectric absorption, Wisconsin cross-sections.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).

    See Also
    --------
    XSzwabs

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelWabs.html

    """

    __function__ = "xsabsw"

    def __init__(self, name='wabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        XSMultiplicativeModel.__init__(self, name, (self.nH,))


class XSwndabs(XSMultiplicativeModel):
    """The XSPEC wndabs model: photo-electric absorption, warm absorber.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    WindowE
        The window energy, in keV.

    See Also
    --------
    XSzwndabs

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelWndabs.html

    """

    __function__ = "xswnab"

    def __init__(self, name='wndabs'):
        self.nH = Parameter(name, 'nH', 1., 0., 10., 0.0, hugeval, '10^22 atoms / cm^2')
        self.WindowE = Parameter(name, 'WindowE', 1., .05, 20., 0.0, hugeval, 'keV')
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.WindowE))


class XSxion(XSMultiplicativeModel):
    """The XSPEC xion model: reflected spectrum of photo-ionized accretion disk/ring.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``lxld`` parameter has been renamed ``lxovrld`` to match the
       XSPEC definition. The name ``lxld`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    height
        The height of the source above the disk (in Schwarzschild radii).
    lxovrld
        The ratio of the X-ray source luminosity to that of the disk.
    rate
        The accretion rate (in Eddington units).
    cosAng
        The cosine of the inclination angle (1 is face on).
    inner
        The inner radius of the disk (in Schwarzschild radii).
    outer
        The er radius of the disk (in Schwarzschild radii).
    index
        The photon index of the source.
    redshift
        The redshift of the absorber.
    Feabun
        The Fe abundance relative to Solar: see [1]_ for more details.
    E_cut
        The exponential high energy cut-off energy for the source.
    Ref_type
        See [1]_ for details.
    Ref_smear
        See [1]_ for details.
    Geometry
        See [1]_ for details.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelXion.html

    """

    __function__ = "xsxirf"

    def __init__(self, name='xion'):
        self.height = Parameter(name, 'height', 5., 0.0, 1.e2, 0.0, hugeval, 'r_s')
        self.lxovrld = Parameter(name, 'lxovrld', 0.3, 0.02, 100, 0.0, hugeval, aliases=["lxld"])
        self.rate = Parameter(name, 'rate', 0.05, 1.e-3, 1., 0.0, hugeval)
        self.cosAng = Parameter(name, 'cosAng', 0.9, 0., 1., 0.0, hugeval)
        self.inner = Parameter(name, 'inner', 3., 2., 1.e3, 0.0, hugeval, 'r_s')
        self.outer = Parameter(name, 'outer', 100., 2.1, 1.e5, 0.0, hugeval, 'r_s')
        self.index = Parameter(name, 'index', 2.0, 1.6, 2.2, 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.Feabun = Parameter(name, 'Feabun', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.E_cut = Parameter(name, 'E_cut', 150., 20., 300., 0.0, hugeval, 'keV')
        self.Ref_type = Parameter(name, 'Ref_type', 1., 1., 3., 0.0, hugeval, frozen=True)
        self.Rel_smear = Parameter(name, 'Rel_smear', 4., 1., 4., 0.0, hugeval, frozen=True)
        self.Geometry = Parameter(name, 'Geometry', 1., 1., 4., 0.0, hugeval, frozen=True)

        XSMultiplicativeModel.__init__(self, name, (self.height, self.lxovrld, self.rate, self.cosAng, self.inner, self.outer, self.index, self.redshift, self.Feabun, self.E_cut, self.Ref_type, self.Rel_smear, self.Geometry))


class XSzdust(XSMultiplicativeModel):
    """The XSPEC zdust model: extinction by dust grains.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``EBV`` parameter has been renamed ``E_BmV`` to match the
       XSPEC definition. The name ``EBV`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    method
        The model to use: 1 is Milky Way, 2 is LMC, 3 is SMC.
        This parameter can not be thawed.
    E_BmV
        The color excess, E(B-V).
    Rv
        The ratio of total to selective extinction.
    redshift
        The redshift of the absorber.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZdust.html

    """

    __function__ = "mszdst"

    def __init__(self, name='zdust'):
        self.method = Parameter(name, 'method', 1, 1, 3, 1, 3, alwaysfrozen=True)
        self.E_BmV = Parameter(name, 'E_BmV', 0.1, 0.0, 100., 0.0, hugeval, aliases=["EBV"])
        self.Rv = Parameter(name, 'Rv', 3.1, 0.0, 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0.0, 0.0, 20., 0.0, hugeval, 'z', True)

        XSMultiplicativeModel.__init__(self, name, (self.method, self.E_BmV, self.Rv, self.redshift))


class XSzedge(XSMultiplicativeModel):
    """The XSPEC zedge model: absorption edge.

    The model is described at [1]_.

    Attributes
    ----------
    edgeE
        The threshold edge, in keV.
    MaxTau
        The absorption depth at the threshold.
    redshift
        The redshift of the edge.

    See Also
    --------
    XSedge

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelEdge.html

    """

    __function__ = "xszedg"

    def __init__(self, name='zedge'):
        self.edgeE = Parameter(name, 'edgeE', 7.0, 0., 100., 0.0, hugeval, 'keV')
        self.MaxTau = Parameter(name, 'MaxTau', 1., 0., 5., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.edgeE, self.MaxTau, self.redshift))


class XSzhighect(XSMultiplicativeModel):
    """The XSPEC zhighect model: high-energy cutoff.

    The model is described at [1]_.

    Attributes
    ----------
    cutoffE
        The cut-off energy, in keV.
    foldE
        The e-folding energy, in keV.
    redshift
        The redshift.

    See Also
    --------
    XShighecut

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelHighecut.html

    """

    __function__ = "xszhcu"

    def __init__(self, name='zhighect'):
        self.cutoffE = Parameter(name, 'cutoffE', 10., 1.e-2, 100., 0.0, hugeval, 'keV')
        self.foldE = Parameter(name, 'foldE', 15., 1.e-2, 100., 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.cutoffE, self.foldE, self.redshift))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.cutoffE, **kwargs)


class XSzpcfabs(XSMultiplicativeModel):
    """The XSPEC zpcfabs model: partial covering fraction absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    CvrFract
        The covering fraction.
    redshift
        The redshift.

    See Also
    --------
    XSpcfabs, XSphabs

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPcfabs.html

    """

    __function__ = "xszabp"

    def __init__(self, name='zpcfabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.CvrFract = Parameter(name, 'CvrFract', 0.5, 0.05, 0.95, 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.CvrFract, self.redshift))


class XSzphabs(XSMultiplicativeModel):
    """The XSPEC zphabs model: photoelectric absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSphabs, XSvphabs, XSzvphabs

    Notes
    -----
    The `set_xsxsect` function changes the cross sections used by this
    model. The `set_xsabund` function changes the relative abundances of
    the elements.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPhabs.html

    """

    __function__ = "xszphb"

    def __init__(self, name='zphabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.redshift))


class XSzxipcf(XSMultiplicativeModel):
    """The XSPEC zxipcf model: partial covering absorption by partially ionized material.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``logxi`` parameter has been renamed ``log_xi`` to match the
       XSPEC definition. The name ``logxi`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    nH
        The column density, in units of 10^22 cm^2.
    log_xi
        The log of xi: see [1]_ for more details.
    CvrFract
        The covering fraction.
    redshift
        The redshift of the source.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZxipcf.html

    """

    __function__ = _f77_or_c_12100("zxipcf")

    def __init__(self, name='zxipcf'):
        self.Nh = Parameter(name, 'Nh', 10, 0.05, 500, 0.0, hugeval, '10^22 atoms / cm^2')
        self.log_xi = Parameter(name, 'log_xi', 3, -3, 6, -hugeval, hugeval, aliases=["logxi"])
        self.CvrFract = Parameter(name, 'CvrFract', 0.5, 0., 1., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)

        XSMultiplicativeModel.__init__(self, name, (self.Nh, self.log_xi, self.CvrFract, self.redshift))


class XSzredden(XSMultiplicativeModel):
    """The XSPEC zredden model: redshifted version of redden.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``EBV`` parameter has been renamed ``E_BmV`` to match the
       XSPEC definition. The name ``EBV`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    E_BmV
        The value of E(B-v) for the line of sight to the source.
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSredden

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZredden.html

    """

    __function__ = "xszcrd"

    def __init__(self, name='zredden'):
        self.E_BmV = Parameter(name, 'E_BmV', 0.05, 0., 10., 0.0, hugeval, aliases=["EBV"])
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)

        XSMultiplicativeModel.__init__(self, name, (self.E_BmV, self.redshift))


class XSzsmdust(XSMultiplicativeModel):
    """The XSPEC zsmdust model: extinction by dust grains in starburst galaxies.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``EBV`` parameter has been renamed ``E_BmV`` to match the
       XSPEC definition. The name ``EBV`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    E_BmV
        The value of E(B-v) for the line of sight to the source.
    ExtIndex
        The spectral index of the extinction curve.
    Rv
        The ratio of total to selective extinction.
    redshift
        The redshift of the absorber.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZsmdust.html

    """

    __function__ = "msldst"

    def __init__(self, name='zsmdust'):
        self.E_BmV = Parameter(name, 'E_BmV', 0.1, 0.0, 100., 0.0, hugeval, aliases=["EBV"])
        self.ExtIndex = Parameter(name, 'ExtIndex', 1.0, -10.0, 10., -hugeval, hugeval)
        self.Rv = Parameter(name, 'Rv', 3.1, 0.0, 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0.0, 0.0, 20., 0.0, hugeval, 'z', True)

        XSMultiplicativeModel.__init__(self, name, (self.E_BmV, self.ExtIndex, self.Rv, self.redshift))


class XSzTBabs(XSMultiplicativeModel):
    """The XSPEC zTBabs model: ISM grain absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSTBabs, XSTBfeo, XSTBgas, XSTBgrain, XSTBpcf, XSTBrel, XSTBvarabs

    Notes
    -----
    The `set_xsabund` function changes the relative abundances of
    the elements, in particular the "wilm" setting.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_ztbabs"

    def __init__(self, name='ztbabs'):
        self.nH = Parameter(name, 'nH', 1., 0., 1E5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.redshift))


class XSzvarabs(XSMultiplicativeModel):
    """The XSPEC zvarabs model: photoelectric absorption.

    The model is described at [1]_.

    Attributes
    ----------
    H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Cl, Ar, Ca, Cr, Fe, Co, Ni
        See [1]_ for a description of the units.
    redshift
        The redshift of the absorber

    See Also
    --------
    XSvarabs, XSzvphabs

    Notes
    -----
    The `set_xsxsect` function changes the cross sections used by this
    model. The `set_xsabund` function changes the relative abundances of
    the elements.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelVarabs.html

    """

    __function__ = "xszvab"

    def __init__(self, name='zvarabs'):
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, hugeval, 'sH22', True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, 'sHe22', True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, 'sC22', True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, 'sN22', True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, 'sO22', True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, 'sNe22', True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, 'sNa22', True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, 'sMg22', True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, 'sAl22', True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, 'sSi22', True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, 'sS22', True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, hugeval, 'sCl22', True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, 'sAr22', True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, 'sCa22', True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, hugeval, 'sCr22', True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, 'sFe22', True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, hugeval, 'sCo22', True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, 'sNi22', True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.H, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Cl, self.Ar, self.Ca, self.Cr, self.Fe, self.Co, self.Ni, self.redshift))


class XSzvfeabs(XSMultiplicativeModel):
    """The XSPEC zvfeabs model: photoelectric absorption with free Fe edge energy.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    metals
        The abundance relative to solar.
    FEabun
        The iron abundance relative to solar.
    FEKedge
        The Fe K edge energy, in keV.
    redshift
        The redshift of the absorber.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZvfeabs.html

    """

    __function__ = "xszvfe"

    def __init__(self, name='zvfeabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.metals = Parameter(name, 'metals', 1., 0.0, 100., 0.0, hugeval)
        self.FEabun = Parameter(name, 'FEabun', 1., 0.0, 100., 0.0, hugeval)
        self.FEKedge = Parameter(name, 'FEKedge', 7.11, 7.0, 9.5, 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.metals, self.FEabun, self.FEKedge, self.redshift))


class XSzvphabs(XSMultiplicativeModel):
    """The XSPEC zvphabs model: photoelectric absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Cl, Ar, Ca, Cr, Fe, Co, Ni
        The abundance of the element in solar units.
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSphabs, XSvphabs, XSzphabs

    Notes
    -----
    The `set_xsxsect` function changes the cross sections used by this
    model. The `set_xsabund` function changes the relative abundances of
    the elements.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPhabs.html

    """

    __function__ = "xszvph"

    def __init__(self, name='zvphabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Cl, self.Ar, self.Ca, self.Cr, self.Fe, self.Co, self.Ni, self.redshift))


class XSzwabs(XSMultiplicativeModel):
    """The XSPEC zwabs model: photoelectric absorption, Wisconsin cross-sections.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSwabs

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelWabs.html

    """

    __function__ = "xszabs"

    def __init__(self, name='zwabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.redshift))


class XSzwndabs(XSMultiplicativeModel):
    """The XSPEC zwndabs model: photo-electric absorption, warm absorber.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    WindowE
        The window energy, in keV.
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSwndabs

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelWndabs.html

    """

    __function__ = "xszwnb"

    def __init__(self, name='zwndabs'):
        self.nH = Parameter(name, 'nH', 1., 0., 10., 0.0, hugeval, '10^22 atoms / cm^2')
        self.WindowE = Parameter(name, 'WindowE', 1., .05, 20., 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.WindowE, self.redshift))


### New additive and multiplicative models, as of XSPEC 12.7


class XScplinear(XSAdditiveModel):
    """The XSPEC cplinear model: a non-physical piecewise-linear model for low count background spectra.

    The model is described at [1]_.

    Attributes
    ----------
    energy00, energy01, energy02, energy03, energy04, energy05,
    energy06, energy07, energy08, energy09
        Energy in keV.
    log_rate01, log_rate02, log_rate03, log_rate04, log_rate05,
    log_rate06, log_rate07, log_rate08, log_rate09
        Log of the rate.
    norm
        The normalization of the model.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCplinear.html

    """

    __function__ = "C_cplinear"

    def __init__(self, name='cplinear'):
        self.energy00 = Parameter(name, 'energy00', 0.5, min=0.0, max=10.0, units='keV', alwaysfrozen=True)
        self.energy01 = Parameter(name, 'energy01', 1., min=0.0, max=10.0, units='keV', alwaysfrozen=True)
        self.energy02 = Parameter(name, 'energy02', 1.5, min=0.0, max=10.0, units='keV', alwaysfrozen=True)
        self.energy03 = Parameter(name, 'energy03', 2., min=0.0, max=10.0, units='keV', alwaysfrozen=True)
        self.energy04 = Parameter(name, 'energy04', 3., min=0.0, max=10.0, units='keV', alwaysfrozen=True)
        self.energy05 = Parameter(name, 'energy05', 4., min=0.0, max=10.0, units='keV', alwaysfrozen=True)
        self.energy06 = Parameter(name, 'energy06', 5., min=0.0, max=10.0, units='keV', alwaysfrozen=True)
        self.energy07 = Parameter(name, 'energy07', 6., min=0.0, max=10.0, units='keV', alwaysfrozen=True)
        self.energy08 = Parameter(name, 'energy08', 7., min=0.0, max=10.0, units='keV', alwaysfrozen=True)
        self.energy09 = Parameter(name, 'energy09', 8., min=0.0, max=10.0, units='keV', alwaysfrozen=True)
        self.log_rate00 = Parameter(name, 'log_rate00', 0., -19.0, 19.0, -hugeval, hugeval, frozen=True)
        self.log_rate01 = Parameter(name, 'log_rate01', 1., -19.0, 19.0, -hugeval, hugeval, frozen=True)
        self.log_rate02 = Parameter(name, 'log_rate02', 0., -19.0, 19.0, -hugeval, hugeval, frozen=True)
        self.log_rate03 = Parameter(name, 'log_rate03', 1., -19.0, 19.0, -hugeval, hugeval, frozen=True)
        self.log_rate04 = Parameter(name, 'log_rate04', 0., -19.0, 19.0, -hugeval, hugeval, frozen=True)
        self.log_rate05 = Parameter(name, 'log_rate05', 1., -19.0, 19.0, -hugeval, hugeval, frozen=True)
        self.log_rate06 = Parameter(name, 'log_rate06', 0., -19.0, 19.0, -hugeval, hugeval, frozen=True)
        self.log_rate07 = Parameter(name, 'log_rate07', 1., -19.0, 19.0, -hugeval, hugeval, frozen=True)
        self.log_rate08 = Parameter(name, 'log_rate08', 0., -19.0, 19.0, -hugeval, hugeval, frozen=True)
        self.log_rate09 = Parameter(name, 'log_rate09', 1., -19.0, 19.0, -hugeval, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.energy00, self.energy01, self.energy02, self.energy03,
                                              self.energy04, self.energy05, self.energy06, self.energy07,
                                              self.energy08, self.energy09, self.log_rate00, self.log_rate01,
                                              self.log_rate02, self.log_rate03, self.log_rate04, self.log_rate05,
                                              self.log_rate06, self.log_rate07, self.log_rate08, self.log_rate09,
                                              self.norm))


class XSeqpair(XSAdditiveModel):
    """The XSPEC eqpair model: Paolo Coppi's hybrid (thermal/non-thermal) hot plasma emission models.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``l_hl_s``, ``l_ntl_h``, and ``AbHe`` parameters have been
       renamed to ``l_hovl_s``, ``l_ntol_h``, and ``Ab_met`` respectively
       to match the XSPEC definition. The names ``l_hl_s``, ``l_ntl_h``,
       and ``AbHe`` can still be used to access the parameters, but
       they will be removed in a future release.

    Attributes
    ----------
    l_hovl_s
        The ratio of the hard to soft compactness, l_h / l_s.
    l_bb
        The soft photon compactness.
    kT_bb
        The temperature of the blackbody if greater than 0.
        When less than zero then the absolute value is used as
        the T_max parameter of the ``XSdispkpn`` model. The units
        are in eV.
    l_ntol_h
        The fraction of power supplied to energetic particles which
        goes into accelerating non-thermal particles, l_nt / l_h.
    tau_p
        The Thomson scattering depth.
    radius
        The size of the scattering region in cm.
    g_min
        The minimum Lorentz factor of the pairs.
    g_max
        The maximum Lorentz factor of the pairs.
    G_inj
        If less than zero then the non-thermal spectrum is assumed
        mono-energetic at g_max, otherwise a power law is used from
        g_min to g_max.
    pairinj
        If zero then accelerated particles are electrons from thermal
        pool. If one then accelerated particles are electrons and
        positrons.
    cosIncl
        The cosine of the inclination angle of the reflecting material
        to the line of sight.
    Refl
        The fraction of the scattering region's emission intercepted
        by reflecting material.
    Fe_abund
        The iron abundance with respect to solar.
    Ab_met
        The abundance of the other metals with respect to solar.
    T_disk
        The temperature of the reflecting disk, in K.
    xi
        The ionization parameter of the reflector.
    Beta
        The power-law index with radius of disk reflection
        emissivity.
    Rin
        The inner radius of the reflecting material, in units of
        GM/c^2.
    Rout
        The outer radius of the reflecting material, in units of
        GM/c^2.
    redshift
        The redshift of the source.
    norm
        The normalization of the model: see [1]_ for more details.

    See Also
    --------
    XScompth, XSeqtherm

    Notes
    -----
    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the EQPAIR_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelEqpair.html

    """

    __function__ = "C_xseqpair"

    def __init__(self, name='eqpair'):
        self.l_hovl_s = Parameter(name, 'l_hovl_s', 1., 1e-6, 1.e6, 0.0, hugeval, aliases=["l_hl_s"])
        self.l_bb = Parameter(name, 'l_bb', 100., 0., 1.e4, 0.0, hugeval)
        self.kT_bb = Parameter(name, 'kT_bb', 200., 1., 4e5, 0.0, hugeval, 'eV', True)
        self.l_ntol_h = Parameter(name, 'l_ntol_h', 0.5, 0., 0.9999, 0.0, hugeval, aliases=["l_ntl_h"])
        self.tau_p = Parameter(name, 'tau_p', 0.1, 1e-4, 10., 0.0, hugeval, frozen=True)
        self.radius = Parameter(name, 'radius', 1.e7, 1.e5, 1.e16, 0.0, hugeval, 'cm', True)
        self.g_min = Parameter(name, 'g_min', 1.3, 1.2, 1.e3, 0.0, hugeval, frozen=True)
        self.g_max = Parameter(name, 'g_max', 1.e3, 5., 1.e4, 0.0, hugeval, frozen=True)
        self.G_inj = Parameter(name, 'G_inj', 2., 0., 5., 0.0, hugeval, frozen=True)
        self.pairinj = Parameter(name, 'pairinj', 0., 0., 1., 0.0, hugeval, frozen=True)
        self.cosIncl = Parameter(name, 'cosIncl', 0.50, 0.05, 0.95, 0.0, hugeval, frozen=True)
        self.Refl = Parameter(name, 'Refl', 1., 0., 2., 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.1, 10., 0.0, hugeval, frozen=True)
        self.Ab_met = Parameter(name, 'Ab_met', 1.0, 0.1, 10., 0.0, hugeval, frozen=True, aliases=["AbHe"])
        self.T_disk = Parameter(name, 'T_disk', 1.e6, 1e4, 1e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 0.0, 0.0, 1000.0, 0.0, hugeval)
        self.Beta = Parameter(name, 'Beta', -10., -10., 10., -hugeval, hugeval, frozen=True)
        self.Rin = Parameter(name, 'Rin', 10., 6.001, 1.e3, 0.0, hugeval, 'M', True)
        self.Rout = Parameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, hugeval, 'M', True)
        self.redshift = Parameter(name, 'redshift', 0., 0., 4., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.l_hovl_s, self.l_bb, self.kT_bb, self.l_ntol_h, self.tau_p, self.radius, self.g_min, self.g_max, self.G_inj, self.pairinj, self.cosIncl, self.Refl, self.Fe_abund, self.Ab_met, self.T_disk, self.xi, self.Beta, self.Rin, self.Rout, self.redshift, self.norm))


class XSeqtherm(XSAdditiveModel):
    """The XSPEC eqtherm model: Paolo Coppi's hybrid (thermal/non-thermal) hot plasma emission models.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``l_hl_s``, ``l_ntl_h``, and ``AbHe`` parameters have been
       renamed to ``l_hovl_s``, ``l_ntol_h``, and ``Ab_met`` respectively
       to match the XSPEC definition. The names ``l_hl_s``, ``l_ntl_h``,
       and ``AbHe`` can still be used to access the parameters, but
       they will be removed in a future release.

    Attributes
    ----------
    l_hovl_s
        The ratio of the hard to soft compactness, l_h / l_s.
    l_bb
        The soft photon compactness.
    kT_bb
        The temperature of the blackbody if greater than 0.
        When less than zero then the absolute value is used as
        the T_max parameter of the ``XSdispkpn`` model. The units
        are in eV.
    l_ntol_h
        The fraction of power supplied to energetic particles which
        goes into accelerating non-thermal particles, l_nt / l_h.
    tau_p
        The Thomson scattering depth.
    radius
        The size of the scattering region in cm.
    g_min
        The minimum Lorentz factor of the pairs.
    g_max
        The maximum Lorentz factor of the pairs.
    G_inj
        If less than zero then the non-thermal spectrum is assumed
        mono-energetic at g_max, otherwise a power law is used from
        g_min to g_max.
    pairinj
        If zero then accelerated particles are electrons from thermal
        pool. If one then accelerated particles are electrons and
        positrons.
    cosIncl
        The cosine of the inclination angle of the reflecting material
        to the line of sight.
    Refl
        The fraction of the scattering region's emission intercepted
        by reflecting material.
    Fe_abund
        The iron abundance with respect to solar.
    Ab_met
        The abundance of the other metals with respect to solar.
    T_disk
        The temperature of the reflecting disk, in K.
    xi
        The ionization parameter of the reflector.
    Beta
        The power-law index with radius of disk reflection
        emissivity.
    Rin
        The inner radius of the reflecting material, in units of
        GM/c^2.
    Rout
        The outer radius of the reflecting material, in units of
        GM/c^2.
    redshift
        The redshift of the source.
    norm
        The normalization of the model: see [1]_ for more details.

    See Also
    --------
    XScompth, XSeqpair

    Notes
    -----
    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the EQPAIR_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelEqpair.html

    """

    __function__ = "C_xseqth"

    def __init__(self, name='eqtherm'):
        self.l_hovl_s = Parameter(name, 'l_hovl_s', 1., 1e-6, 1.e6, 0.0, hugeval, aliases=["l_hl_s"])
        self.l_bb = Parameter(name, 'l_bb', 100., 0., 1.e4, 0.0, hugeval)
        self.kT_bb = Parameter(name, 'kT_bb', 200., 1., 4e5, 0.0, hugeval, 'eV', True)
        self.l_ntol_h = Parameter(name, 'l_ntol_h', 0.5, 0., 0.9999, 0.0, hugeval, aliases=["l_ntl_h"])
        self.tau_p = Parameter(name, 'tau_p', 0.1, 1e-4, 10., 0.0, hugeval, frozen=True)
        self.radius = Parameter(name, 'radius', 1.e7, 1.e5, 1.e16, 0.0, hugeval, 'cm', True)
        self.g_min = Parameter(name, 'g_min', 1.3, 1.2, 1.e3, 0.0, hugeval, frozen=True)
        self.g_max = Parameter(name, 'g_max', 1.e3, 5., 1.e4, 0.0, hugeval, frozen=True)
        self.G_inj = Parameter(name, 'G_inj', 2., 0., 5., 0.0, hugeval, frozen=True)
        self.pairinj = Parameter(name, 'pairinj', 0., 0., 1., 0.0, hugeval, frozen=True)
        self.cosIncl = Parameter(name, 'cosIncl', 0.50, 0.05, 0.95, 0.0, hugeval, frozen=True)
        self.Refl = Parameter(name, 'Refl', 1., 0., 2., 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.1, 10., 0.0, hugeval, frozen=True)
        self.Ab_met = Parameter(name, 'Ab_met', 1.0, 0.1, 10., 0.0, hugeval, frozen=True, aliases=["AbHe"])
        self.T_disk = Parameter(name, 'T_disk', 1.e6, 1e4, 1e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 0.0, 0.0, 1000.0, 0.0, hugeval)
        self.Beta = Parameter(name, 'Beta', -10., -10., 10., -hugeval, hugeval, frozen=True)
        self.Rin = Parameter(name, 'Rin', 10., 6.001, 1.e3, 0.0, hugeval, 'M', True)
        self.Rout = Parameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, hugeval, 'M', True)
        self.redshift = Parameter(name, 'redshift', 0., 0., 4., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.l_hovl_s, self.l_bb, self.kT_bb, self.l_ntol_h, self.tau_p, self.radius, self.g_min, self.g_max, self.G_inj, self.pairinj, self.cosIncl, self.Refl, self.Fe_abund, self.Ab_met, self.T_disk, self.xi, self.Beta, self.Rin, self.Rout, self.redshift, self.norm))


# It is not obvious from the XSPEC documentation what theta, showbb,
# and RefOn are.
#
class XScompth(XSAdditiveModel):
    """The XSPEC compth model: Paolo Coppi's hybrid (thermal/non-thermal) hot plasma emission models.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``AbHe`` parameter has been renamed ``Ab_met`` to match the
       XSPEC definition. The name ``AbHe`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    theta
    showbb
    kT_bb
        The temperature of the blackbody if greater than 0.
        When less than zero then the absolute value is used as
        the T_max parameter of the ``XSdispkpn`` model. The units
        are in eV.
    RefOn
    tau_p
        The Thomson scattering depth.
    radius
        The size of the scattering region in cm.
    g_min
        The minimum Lorentz factor of the pairs.
    g_max
        The maximum Lorentz factor of the pairs.
    G_inj
        If less than zero then the non-thermal spectrum is assumed
        mono-energetic at g_max, otherwise a power law is used from
        g_min to g_max.
    pairinj
        If zero then accelerated particles are electrons from thermal
        pool. If one then accelerated particles are electrons and
        positrons.
    cosIncl
        The cosine of the inclination angle of the reflecting material
        to the line of sight.
    Refl
        The fraction of the scattering region's emission intercepted
        by reflecting material.
    Fe_abund
        The iron abundance with respect to solar.
    Ab_met
        The abundance of the other metals with respect to solar.
    T_disk
        The temperature of the reflecting disk, in K.
    xi
        The ionization parameter of the reflector.
    Beta
        The power-law index with radius of disk reflection
        emissivity.
    Rin
        The inner radius of the reflecting material, in units of
        GM/c^2.
    Rout
        The outer radius of the reflecting material, in units of
        GM/c^2.
    redshift
        The redshift of the source.
    norm
        The normalization of the model: see [1]_ for more details.

    See Also
    --------
    XSeqpair, XSeqtherm

    Notes
    -----
    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the EQPAIR_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelEqpair.html

    """

    __function__ = "C_xscompth"

    def __init__(self, name='compth'):
        self.theta = Parameter(name, 'theta', 1., 1e-6, 1.e6, 0.0, hugeval, 'keV')
        self.showbb = Parameter(name, 'showbb', 1.0, 0., 1.e4, 0.0, hugeval, frozen=True)
        self.kT_bb = Parameter(name, 'kT_bb', 200., 1., 4e5, 0.0, hugeval, 'eV', True)
        self.RefOn = Parameter(name, 'RefOn', -1.0, -2.0, 2.0, -hugeval, hugeval, frozen=True)
        self.tau_p = Parameter(name, 'tau_p', 0.1, 1e-4, 10., 0.0, hugeval, frozen=True)
        self.radius = Parameter(name, 'radius', 1.e7, 1.e5, 1.e16, 0.0, hugeval, 'cm', True)
        self.g_min = Parameter(name, 'g_min', 1.3, 1.2, 1.e3, 0.0, hugeval, frozen=True)
        self.g_max = Parameter(name, 'g_max', 1.e3, 5., 1.e4, 0.0, hugeval, frozen=True)
        self.G_inj = Parameter(name, 'G_inj', 2., 0., 5., 0.0, hugeval, frozen=True)
        self.pairinj = Parameter(name, 'pairinj', 0., 0., 1., 0.0, hugeval, frozen=True)
        self.cosIncl = Parameter(name, 'cosIncl', 0.50, 0.05, 0.95, 0.0, hugeval, frozen=True)
        self.Refl = Parameter(name, 'Refl', 1., 0., 2., 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.1, 10., 0.0, hugeval, frozen=True)
        self.Ab_met = Parameter(name, 'Ab_met', 1.0, 0.1, 10., 0.0, hugeval, frozen=True, aliases=["AbHe"])
        self.T_disk = Parameter(name, 'T_disk', 1.e6, 1e4, 1e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 0.0, 0.0, 1000.0, 0.0, hugeval)
        self.Beta = Parameter(name, 'Beta', -10., -10., 10., -hugeval, hugeval, frozen=True)
        self.Rin = Parameter(name, 'Rin', 10., 6.001, 1.e3, 0.0, hugeval, 'M', True)
        self.Rout = Parameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, hugeval, 'M', True)
        self.redshift = Parameter(name, 'redshift', 0., 0., 4., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.theta, self.showbb, self.kT_bb, self.RefOn, self.tau_p, self.radius, self.g_min, self.g_max, self.G_inj, self.pairinj, self.cosIncl, self.Refl, self.Fe_abund, self.Ab_met, self.T_disk, self.xi, self.Beta, self.Rin, self.Rout, self.redshift, self.norm))


class XSbvvapec(XSAdditiveModel):
    """The XSPEC bvvapec model: velocity broadened APEC thermal plasma model.

    The model is described at [1]_, with ``XSbapec`` describing how
    it is implemented in Sherpa.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element in solar units.
    Redshift
        The redshift of the plasma.
    Velocity
        The gaussian sigma of the velocity broadening, in km/s.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSapec, XSbapec, XSbvapec, XSvapec, XSvvapec

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBapec.html

    """

    __function__ = "C_bvvapec" if equal_or_greater_than("12.9.1") else "xsbvvp"

    def __init__(self, name='bvvapec'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0, hugeval, 'keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -hugeval, hugeval, frozen=True)
        self.Velocity = Parameter(name, 'Velocity', 0., 0., 1.e6, 0.0, hugeval, 'km/s', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Redshift, self.Velocity, self.norm))


@version_at_least("12.10.0")
class XSbvvrnei(XSAdditiveModel):
    """The XSPEC bvvrnei model: velocity-broadened non-equilibrium recombining collisional plasma.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keyword "NEIVERS".

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    kT_init
        The initial temperature of the plasma, in keV.
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element in solar units.
    Tau
        The ionization timescale in units of s/cm^3.
    Redshift
        The redshift of the plasma.
    Velocity
        The gaussian sigma of the velocity broadening, in km/s.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbrnei, XSbvrnei, XSnei, XSgnei, XSrnei, XSvrnei, XSvvrnei

    Notes
    -----
    This model is only available when used with XSPEC 12.10.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBrnei.html

    """

    __function__ = "C_bvvrnei"

    def __init__(self, name='bvvrnei'):
        self.kT = Parameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = Parameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        maxval = 1000.0
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -0.999, 10., frozen=True)
        self.Velocity = Parameter(name, 'Velocity', 0., 0., 1.e6, 0.0, 1.e6,
                                  units='km/s', frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.Redshift, self.Velocity, self.norm))


class XSvvapec(XSAdditiveModel):
    """The XSPEC vvapec model: APEC emission spectrum.

    The model is described at [1]_, with ``XSapec`` describing how
    it is implemented in Sherpa.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element in solar units.
    Redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSapec, XSbapec, XSbvapec, XSbvvapec, XSvapec

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelApec.html

    """

    __function__ = "C_vvapec" if equal_or_greater_than("12.9.1") else "xsvvap"

    def __init__(self, name='vvapec'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0, hugeval, 'keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -hugeval, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Redshift, self.norm))


class XSzigm(XSMultiplicativeModel):
    """The XSPEC zigm model: UV/Optical attenuation by the intergalactic medium.

    The model is described at [1]_.

    Attributes
    ----------
    redshift
        The redshift of the absorber.
        This parameter can not be thawed.
    model
        The model to use: 0 is Madau, 1 is Meiksin.
        This parameter can not be thawed.
    lyman_limit
        Should photoelectric absorption be included (1), or not (0).
        This parameter can not be thawed.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZigm.html

    """

    __function__ = "zigm"

    def __init__(self, name='zigm'):
        self.redshift = Parameter(name, 'redshift', 0.0, alwaysfrozen=True)
        self.model = Parameter(name, 'model', 0, alwaysfrozen=True)
        self.lyman_limit = Parameter(name, 'lyman_limit', 1, alwaysfrozen=True)

        XSMultiplicativeModel.__init__(self, name, (self.redshift, self.model, self.lyman_limit))


## Here are the seven new additive models from XSPEC 12.7.1


class XSgadem(XSAdditiveModel):
    """The XSPEC gadem model: plasma emission, multi-temperature with gaussian distribution of emission measure.

    The model is described at [1]_. The ``set_xsabund`` and ``get_xsabund``
    functions change and return the current settings for the relative
    abundances of the metals. See the ``XSapec`` documentation for settings
    relevant to the APEC model (i.e. when ``switch=2``).

    Attributes
    ----------
    Tmean
        The mean temperature for the gaussian emission measure
        distribution, in keV.
    Tsigma
        The sigma of the temperature distribution for the gaussian
        emission measure, in keV.
    nH
        H density, in cm^-3.
    abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function.
    Redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model.

    See Also
    --------
    XSapec, XSvgadem

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGadem.html

    """

    __function__ = "C_gaussDem"

    def __init__(self, name='gadem'):
        self.Tmean = Parameter(name, 'Tmean', 4.0, 0.01, 10, 0.0, hugeval, 'keV', True)
        self.Tsigma = Parameter(name, 'Tsigma', 0.1, 0.01, 1.e2, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1.0, 1.e-5, 1.e19, 0.0, hugeval, 'cm^-3', True)
        self.abundanc = Parameter(name, 'abundanc', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -hugeval, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 2, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.Tmean, self.Tsigma, self.nH, self.abundanc, self.Redshift,
                                              self.switch, self.norm))


class XSvgadem(XSAdditiveModel):
    """The XSPEC vgadem model: plasma emission, multi-temperature with gaussian distribution of emission measure.

    The model is described at [1]_. The ``set_xsabund`` and ``get_xsabund``
    functions change and return the current settings for the relative
    abundances of the metals. See the ``XSapec`` documentation for settings
    relevant to the APEC model (i.e. when ``switch=2``).

    Attributes
    ----------
    Tmean
        The mean temperature for the gaussian emission measure
        distribution, in keV.
    Tsigma
        The sigma of the temperature distribution for the gaussian
        emission measure, in keV.
    nH
        H density, in cm^-3.
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element in solar units.
    Redshift
        The redshift of the plasma.
    switch
        If 0, the mekal code is run to evaluate the model; if 1
        then interpolation of the mekal data is used; if 2 then
        interpolation of APEC data is used. See [1]_ for more details.
        This parameter can not be thawed.
    norm
        The normalization of the model.

    See Also
    --------
    XSapec, XSgadem

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGadem.html

    """

    __function__ = "C_vgaussDem"

    def __init__(self, name='vgadem'):
        self.Tmean = Parameter(name, 'Tmean', 4.0, 0.01, 10, 0.0, hugeval, 'keV', True)
        self.Tsigma = Parameter(name, 'Tsigma', 0.1, 0.01, 1.e2, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1.0, 1.e-5, 1.e19, 0.0, hugeval, 'cm^-3', True)
        self.He = Parameter(name, 'He', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -hugeval, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 2, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.Tmean, self.Tsigma, self.nH, self.He, self.C, self.N, self.O,
                                              self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca,
                                              self.Fe, self.Ni, self.Redshift, self.switch, self.norm))


class XSeplogpar(XSAdditiveModel):
    """The XSPEC eplogpar model: log-parabolic blazar model with nu-Fnu normalization.

    The model is described at [1]_.

    Attributes
    ----------
    Ep
        The peak energy, in keV, of the nu-Fnu curve.
    beta
        The curvature term.
    norm
        The flux in nu-Fnu units at the energy Ep.

    See Also
    --------
    XSlogpar

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelEplogpar.html

    """

    __function__ = "eplogpar"

    def __init__(self, name='eplogpar'):
        self.Ep = Parameter(name, 'Ep', .1, 1.e-6, 1.e2, 0.0, hugeval, 'keV')
        self.beta = Parameter(name, 'beta', 0.2, -4., 4., -hugeval, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Ep, self.beta, self.norm))


class XSlogpar(XSAdditiveModel):
    """The XSPEC logpar model: log-parabolic blazar model.

    The model is described at [1]_.

    Attributes
    ----------
    alpha
        The slope at the pivot energy.
    beta
        The curvature term.
    pivotE
        The pivot energy, in keV.
    norm
        The normalization of the model: see [1]_ for more details.

    See Also
    --------
    XSeplogpar, XSzlogpar

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLogpar.html

    """

    __function__ = "C_logpar" if equal_or_greater_than("12.10.1") else "logpar"

    def __init__(self, name='logpar'):
        self.alpha = Parameter(name, 'alpha', 1.5, 0., 4., 0.0, hugeval)
        self.beta = Parameter(name, 'beta', 0.2, -4., 4., -hugeval, hugeval)
        self.pivotE = Parameter(name, 'pivotE', 1.0, units='keV',
                                alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.alpha, self.beta, self.pivotE, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


class XSoptxagn(XSAdditiveModel):
    """The XSPEC optxagn model: Colour temperature corrected disc and energetically coupled Comptonisation model for AGN.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``logLLEdd`` parameter has been renamed ``logLoLedd`` to match
       the XSPEC definition. The name ``logLLEdd`` can still be used to
       access the parameter, but this name will be removed in a future
       release.

    Attributes
    ----------
    mass
        The black hole mass in solar masses.
    dist
        The comoving (proper) distance in Mpc.
    logLoLEdd
        The Eddington ratio.
    astar
        The dimensionless black hole spin.
    rcor
        The coronal radius in Rg=GM/c^2. See [1]_ for more details.
    logrout
        The log of the outer radius of the disk in units of Rg. See
        [1]_ for more details.
    kT_e
        The electron temperature for the soft Comptonisation component
        (soft excess), in keV.
    tau
        The optical depth of the soft Comptonisation component. If this
        parameter is negative then only the soft Compton component is used.
    Gamma
        The spectral index of the hard Comptonisation component
        ('power law') which has temperature fixed to 100 keV.
    fpl
        The fraction of the power below rcor which is emitted in the hard
        comptonisation component. If this parameter is negative then only
        the hard Compton component is used.
    fcol
        The colour temperature correction to apply to the disc blackbody
        emission for radii below rcor with effective temperature > tscat.
    tscat
        The effective temperature limit, in K, used in the colour
        temperature correction.
    Redshift
        The redshift.
    norm
        The normalization of the model. It must be frozen.

    See Also
    --------
    XSoptxagnf

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelOptxagn.html

    """

    __function__ = "optxagn"

    def __init__(self, name='optxagn'):
        self.mass = Parameter(name, 'mass', 1e7, 1.0, 1.e9, 0.0, hugeval, 'solar', True)
        self.dist = Parameter(name, 'dist', 100, 0.01, 1.e9, 0.0, hugeval, 'Mpc', True)
        self.logLoLEdd = Parameter(name, 'logLoLEdd', -1., -10., 2, -hugeval, hugeval, aliases=["logLLEdd"])
        self.astar = Parameter(name, 'astar', 0.0, 0., 0.998, 0.0, hugeval, frozen=True)
        self.rcor = Parameter(name, 'rcor', 10.0, 1., 100., 0.0, hugeval, 'rg')
        self.logrout = Parameter(name, 'logrout', 5.0, 3.0, 7.0, 0.0, hugeval, frozen=True)
        self.kT_e = Parameter(name, 'kT_e', 0.2, 0.01, 10, 0.0, hugeval, 'keV')
        self.tau = Parameter(name, 'tau', 10., 0.1, 100, 0.0, hugeval)
        self.Gamma = Parameter(name, 'Gamma', 2.1, 0.5, 5., 0.0, hugeval)
        self.fpl = Parameter(name, 'fpl', 1.e-4, 0.0, 1.e-1, 0.0, hugeval)
        self.fcol = Parameter(name, 'fcol', 2.4, 1.0, 5, 0.0, hugeval, frozen=True)
        self.tscat = Parameter(name, 'tscat', 1.e5, 1e4, 1e5, 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0., 0., 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.mass, self.dist, self.logLoLEdd, self.astar, self.rcor, self.logrout, self.kT_e, self.tau, self.Gamma, self.fpl, self.fcol, self.tscat, self.Redshift, self.norm))


class XSoptxagnf(XSAdditiveModel):
    """The XSPEC optxagnf model: Colour temperature corrected disc and energetically coupled Comptonisation model for AGN.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``logLLEdd`` parameter has been renamed ``logLoLedd`` to match
       the XSPEC definition. The name ``logLLEdd`` can still be used to
       access the parameter, but this name will be removed in a future
       release.

    Attributes
    ----------
    mass
        The black hole mass in solar masses.
    dist
        The comoving (proper) distance in Mpc.
    logLoLEdd
        The Eddington ratio.
    astar
        The dimensionless black hole spin.
    rcor
        The coronal radius in Rg=GM/c^2. See [1]_ for more details.
    logrout
        The log of the outer radius of the disk in units of Rg. See
        [1]_ for more details.
    kT_e
        The electron temperature for the soft Comptonisation component
        (soft excess), in keV.
    tau
        The optical depth of the soft Comptonisation component. If this
        parameter is negative then only the soft Compton component is used.
    Gamma
        The spectral index of the hard Comptonisation component
        ('power law') which has temperature fixed to 100 keV.
    fpl
        The fraction of the power below rcor which is emitted in the hard
        comptonisation component. If this parameter is negative then only
        the hard Compton component is used.
    Redshift
        The redshift.
    norm
        The normalization of the model. It must be frozen.

    See Also
    --------
    XSoptxagn

    Notes
    -----
    The minimum allowed value for the Gamma parameter has been changed
    from 0.5 to 1.05 to match the XSPEC 12.10.0 model.dat file.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelOptxagn.html

    """

    __function__ = "optxagnf"

    def __init__(self, name='optxagnf'):
        self.mass = Parameter(name, 'mass', 1e7, 1.0, 1.e9, 0.0, hugeval, 'solar', True)
        self.dist = Parameter(name, 'dist', 100, 0.01, 1.e9, 0.0, hugeval, 'Mpc', True)
        self.logLoLEdd = Parameter(name, 'logLoLEdd', -1., -10., 2, -hugeval, hugeval, aliases=["logLLEdd"])
        self.astar = Parameter(name, 'astar', 0.0, 0., 0.998, 0.0, hugeval, frozen=True)
        self.rcor = Parameter(name, 'rcor', 10.0, 1., 100., 0.0, hugeval, 'rg')
        self.logrout = Parameter(name, 'logrout', 5.0, 3.0, 7.0, 0.0, hugeval, frozen=True)
        self.kT_e = Parameter(name, 'kT_e', 0.2, 0.01, 10, 0.0, hugeval, 'keV')
        self.tau = Parameter(name, 'tau', 10., 0.1, 100, 0.0, hugeval)
        self.Gamma = Parameter(name, 'Gamma', 2.1, 1.05, 5., 1.05, 10.0)
        self.fpl = Parameter(name, 'fpl', 1.e-4, 0.0, 1., 0.0, hugeval)
        self.Redshift = Parameter(name, 'Redshift', 0., 0., 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.mass, self.dist, self.logLoLEdd, self.astar, self.rcor, self.logrout, self.kT_e, self.tau, self.Gamma, self.fpl, self.Redshift, self.norm))


# DOC-NOTE: the parameter order in the XSPEC documentation is very different
#           to the model.dat file. Or, rel_refl is actually labelled as scale
#           and foldE as Ec.
#
class XSpexmon(XSAdditiveModel):
    """The XSPEC pexmon model: neutral Compton reflection with self-consistent Fe and Ni lines.

    The model is described at [1]_.

    Attributes
    ----------
    PhoIndex
        The power-law photon index.
    foldE
        The cut-off energy (E_c) in keV. Set to 0 for no cut off.
    rel_refl
        The reflection scaling parameter (a value of 1 for an
        isotropic source above the disk, less than 0 for no direct
        component).
    redshift
        The redshift of the source.
    abund
        The abundance of the elements heaver than He relative to their
        solar abundance, as set by the ``set_xsabund`` function.
    Fe_abund
        The iron abundance relative to the solar abundance, as set by
        the ``set_xsabund`` function.
    Incl
        The inclination angle in degrees.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSpexrav, XSpexriv

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPexmon.html

    """

    __function__ = "pexmon"

    def __init__(self, name='pexmon'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 2., 1.1, 2.5, 0.0, hugeval)
        self.foldE = Parameter(name, 'foldE', 1000., 1., 1.e6, 0.0, hugeval, 'keV', True)
        self.rel_refl = Parameter(name, 'rel_refl', -1, -1.e6, 1.e6, -hugeval, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., 0., 4., 0.0, hugeval, frozen=True)
        self.abund = Parameter(name, 'abund', 1., 0.0, 1.e6, 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.0, 100., 0.0, hugeval, frozen=True)
        self.Incl = Parameter(name, 'Incl', 60., 0., 85.0, 0.0, hugeval, 'deg')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.foldE, self.rel_refl, self.redshift, self.abund, self.Fe_abund, self.Incl, self.norm))


class XSagauss(XSAdditiveModel):
    """The XSPEC agauss model: gaussian line profile in wavelength space.

    The model is described at [1]_.

    Attributes
    ----------
    LineE
        The line wavelength, in Angstrom.
    Sigma
        The line width, in Angstrom. A value of zero means a delta function.
    norm
        The flux in the line, in units of photon/cm^2/s.

    See Also
    --------
    XSgaussian, XSzagauss, XSzgauss

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelAgauss.html

    """

    __function__ = "C_agauss"

    def __init__(self, name='agauss'):
        self.LineE = Parameter(name, 'LineE', 10.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.Sigma = Parameter(name, 'Sigma', 1.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Sigma, self.norm))


class XSzagauss(XSAdditiveModel):
    """The XSPEC zagauss model: gaussian line profile in wavelength space.

    The model is described at [1]_.

    Attributes
    ----------
    LineE
        The line wavelength, in Angstrom.
    Sigma
        The line width, in Angstrom. A value of zero means a delta function.
    Redshift
        The redshift of the line.
    norm
        The flux in the line, in units of photon/cm^2/s.

    See Also
    --------
    XSagauss, XSgaussian, XSzgauss

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelAgauss.html

    """

    __function__ = "C_zagauss"

    def __init__(self, name='zagauss'):
        self.LineE = Parameter(name, 'LineE', 10.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.Sigma = Parameter(name, 'Sigma', 1.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10.0, -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Sigma, self.Redshift, self.norm))


class XScompmag(XSAdditiveModel):
    """The XSPEC compmag model: Thermal and bulk Comptonization for cylindrical accretion onto the polar cap of a magnetized neutron star.

    The model is described at [1]_.

    Attributes
    ----------
    kTbb
        The seed blackbody temperature, in keV.
    kTe
        The electron temperature of the accretion column, in keV.
    tau
        The vertical optical depth of the accretion column, with
        electron cross-section equal to 10^-3 of the Thomson cross-section.
    eta
        The index of the velocity profile when the accretion velocity
        increases towards the neutron star (valid when betaflag is 1).
    beta0
        The terminal velocity of the accreting matter at the neutron star
        surface (valid when betaflag is 1).
    r0
        The radius of the accretion column in units of the neutron star
        Schwarzschild radius.
    A
        The albedo at the surface of the neutron star.
    betaflag
        A flag for setting the velocity profile of the accretion
        column, described at [1]_. It has values of 1 or 2 and can
        not be thawed.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCompmag.html

    """

    __function__ = "xscompmag"

    def __init__(self, name='compmag'):
        self.kTbb = Parameter(name, 'kTbb', 1.0, 0.2, 10.0, 0.2, 10.0, 'keV')
        self.kTe = Parameter(name, 'kTe', 5.0, 0.2, 2000.0, 0.2, 2000.0, 'keV')
        self.tau = Parameter(name, 'tau', 0.5, 0.0, 10.0, 0.0, 10.0)
        self.eta = Parameter(name, 'eta', 0.5, 0.01, 1.0, 0.01, 1.0)
        self.beta0 = Parameter(name, 'beta0', 0.57, 1.0e-4, 1.0, 1.0e-4, 1.0)
        self.r0 = Parameter(name, 'r0', 0.25, 1.0e-4, 100.0, 1.0e-4, 100.0)
        self.A = Parameter(name, 'A', 0.001, 0, 1, 0, 1, frozen=True)
        self.betaflag = Parameter(name, 'betaflag', 1, 0, 2, 0, 2, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kTbb, self.kTe, self.tau, self.eta, self.beta0, self.r0, self.A, self.betaflag, self.norm))


class XScomptb(XSAdditiveModel):
    """The XSPEC comptb model: Thermal and bulk Comptonization of a seed blackbody-like spectrum.

    The model is described at [1]_.

    Attributes
    ----------
    kTs
        The temperature of the seed photons, in keV.
    gamma
        The index of the seed photon spectrum.
    alpha
        The energy index of the Comptonization spectrum.
    delta
        The bulk parameter, efficiency of bulk over thermal
        Comptonization.
    kTe
        The temperature of the electrons, in keV.
    log_A
        The log of the illuminating factor parameter (A).
    norm
        The normalization of the seed photon spectrum: it is the
        same definition as used for the ``XSbbody`` model.

    See Also
    --------
    XSbbody

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelComptb.html

    """

    __function__ = "xscomptb"

    def __init__(self, name='comptb'):
        self.kTs = Parameter(name, 'kTs', 1.0, 0.1, 10.0, 0.1, 10.0, units='keV')
        self.gamma = Parameter(name, 'gamma', 3.0, 1.0, 10.0, 1.0, 10.0, frozen=True)
        self.alpha = Parameter(name, 'alpha', 2.0, 0.0, 400.0, 0.0, 400.0)
        self.delta = Parameter(name, 'delta', 20.0, 0.0, 200.0, 0.0, 200.0)
        self.kTe = Parameter(name, 'kTe', 5.0, 0.2, 2000.0, 0.2, 2000.0, units='keV')
        self.log_A = Parameter(name, 'log_A', 0.0, -8.0, 8.0, -8.0, 8.0)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kTs, self.gamma, self.alpha, self.delta, self.kTe, self.log_A, self.norm))


class XSheilin(XSMultiplicativeModel):
    """The XSPEC heilin model: Voigt absorption profiles for He I series.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``redshift`` parameter has been renamed ``z`` to match the
       XSPEC definition. The name ``redshift`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    nHeI
        The He I column density, in 10^22 atoms/cm^2.
    b
        The b value, in km/s.
    z
        The redshift of the absorber.

    See Also
    --------
    XSlyman

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelHeilin.html

    """

    __function__ = "xsphei"

    def __init__(self, name='heilin'):
        self.nHeI = Parameter(name, 'nHeI', 1.e-5, 0.0, 1.e6, 0.0, 1.0e6, '10^22 atoms / cm^2')
        self.b = Parameter(name, 'b', 10.0, 1.0, 1.0e5, 1.0, 1.0e6, units='km/s')
        self.z = Parameter(name, 'z', 0.0, -1.0e-3, 1.0e5, -1.0e-3, 1.0e5, aliases=["redshift"])

        # TODO: correct self.nHei to self.nHeI
        XSMultiplicativeModel.__init__(self, name, (self.nHei, self.b, self.z))


class XSlyman(XSMultiplicativeModel):
    """The XSPEC lyman model: Voigt absorption profiles for H I or He II Lyman series.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``nHeI`` and ``redshift`` parameters have been renamed to
       ``n`` and ``z`` respectively to match the XSPEC definition. The
       names ``nHeI`` and ``redshift`` can still be used to access the
       parameters, but they will be removed in a future release.

    Attributes
    ----------
    n
        The H I or He II column density, in 10^22 atoms/cm^2.
    b
        The b value, in km/s.
    z
        The redshift of the absorber.
    ZA
        The atomic number of the species being calculated.

    See Also
    --------
    XSheilin

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLyman.html

    """

    __function__ = "xslyman"

    def __init__(self, name='lyman'):
        self.n = Parameter(name, 'n', 1.e-5, 0.0, 1.0e6, 0.0, 1.0e6, '10^22 atoms / cm^2', aliases=["nHeI"])
        self.b = Parameter(name, 'b', 10.0, 1.0, 1.0e5, 1.0, 1.0e6, units='km/s')
        self.z = Parameter(name, 'z', 0.0, -1.0e-3, 1.0e5, -1.0e-3, 1.0e5, aliases=["redshift"])
        self.ZA = Parameter(name, 'ZA', 1.0, 1.0, 2.0, 1.0, 2.0)

        XSMultiplicativeModel.__init__(self, name, (self.n, self.b, self.z, self.ZA))


class XSzbabs(XSMultiplicativeModel):
    """The XSPEC zbabs model: EUV ISM attenuation.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``redshift`` parameter has been renamed ``z`` to match the
       XSPEC definition. The name ``redshift`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    nH
        The H column density, in 10^22 atoms/cm^2.
    nHeI
        The He I column density, in 10^22 atoms/cm^2.
    nHeI
        The He II column density, in 10^22 atoms/cm^2.
    z
        The redshift of the absorber.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZbabs.html

    """

    __function__ =  "xszbabs"

    def __init__(self, name='zbabs'):
        self.nH = Parameter(name, 'nH', 1.e-4, 0.0, 1.0e5, 0.0, 1.0e6, '10^22 atoms / cm^2')
        self.nHeI = Parameter(name, 'nHeI', 1.e-5, 0.0, 1.0e5, 0.0, 1.0e6, '10^22 atoms / cm^2')
        self.nHeII = Parameter(name, 'nHeII', 1.e-6, 0.0, 1.0e5, 0.0, 1.0e6, '10^22 atoms / cm^2')
        self.z = Parameter(name, 'z', 0.0, 0.0, 1.0e5, 0.0, 1.0e6, aliases=["redshift"])

        XSMultiplicativeModel.__init__(self, name, (self.nH, self.nHeI, self.nHeII, self.z))


@version_at_least("12.9.1")
class XSbtapec(XSAdditiveModel):
    """The XSPEC btapec model: velocity broadened APEC emission spectrum with separate continuum and line temperatures.

    The model is described at [1]_. The ``set_xsabund`` and ``get_xsabund``
    functions change and return the current settings for the relative
    abundances of the metals. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keywords "APECROOT" and "APEC_TRACE_ABUND".

    Attributes
    ----------
    kT
        Continuum temperature, in keV.
    kTi
        Line temperature, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function and the "APEC_TRACE_ABUND" xset
        keyword (which defaults to 1.0).
    Redshift
        The redshift of the plasma.
    Velocity
        The gaussian sigma of the velocity broadening, in km/s.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbapec, XSbvtapec, XSbvvtapec

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBtapec.html

    """

    __function__ = "C_btapec"

    def __init__(self, name='btapec'):
        self.kT = Parameter(name, 'kT', 1.0, 0.008, 64.0, 0.008, 64.0, 'keV')
        self.kTi = Parameter(name, 'kTi', 1.0, 0.008, 64.0, 0.008, 64.0, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0.0, 5.0, 0.0, 5.0, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10, frozen=True)
        self.Velocity = Parameter(name, 'Velocity', 0.0, 0.0, 1.0e6, 0.0, 1.0e6, 'km/s', frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.Abundanc, self.Redshift, self.Velocity, self.norm))


@version_at_least("12.9.1")
class XSbvtapec(XSAdditiveModel):
    """The XSPEC bvtapec model: velocity broadened APEC emission spectrum with separate continuum and line temperatures.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keywords "APECROOT" and "APEC_TRACE_ABUND".

    Attributes
    ----------
    kT
        Continuum temperature, in keV.
    kTi
        Line temperature, in keV.
    He, C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element in solar units. The trace element
        abundances are determined by the APEC_TRACE_ABUND xset
        keyword (the default is 1.0).
    Redshift
        The redshift of the plasma.
    Velocity
        The gaussian sigma of the velocity broadening, in km/s.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbapec, XSbtapec, XSbvvtapec

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBtapec.html

    """

    __function__ = "C_bvtapec"

    def __init__(self, name='bvtapec'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, 'keV')
        self.kTi = Parameter(name, 'kTi', 6.5, 0.0808, 68.447, 0.0808, 68.447, 'keV')
        self.He = Parameter(name, 'He', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10)
        self.Velocity = Parameter(name, 'Velocity', 0.0, 0.0, 1.0e6, 0.0, 1.0e6, 'km/s')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.He, self.C, self.N, self.O,
                                  self.Ne, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca,
                                  self.Fe, self.Ni, self.Redshift, self.Velocity, self.norm))


@version_at_least("12.9.1")
class XSbvvtapec(XSAdditiveModel):
    """The XSPEC bvvtapec model: velocity broadened APEC emission spectrum with separate continuum and line temperatures.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET "APECROOT"
    parameter.

    Attributes
    ----------
    kT
        Continuum temperature, in keV.
    kTi
        Line temperature, in keV.
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element in solar units.
    Redshift
        The redshift of the plasma.
    Velocity
        The gaussian sigma of the velocity broadening, in km/s.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbapec, XSbtapec

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBtapec.html

    """

    __function__ = "C_bvvtapec"

    def __init__(self, name='bvvtapec'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, 'keV')
        self.kTi = Parameter(name, 'kTi', 6.5, 0.0808, 68.447, 0.0808, 68.447, 'keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10)
        self.Velocity = Parameter(name, 'Velocity', 0.0, 0.0, 1.0e6, 0.0, 1.0e6, 'km/s')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Redshift, self.Velocity, self.norm))


@version_at_least("12.9.1")
class XStapec(XSAdditiveModel):
    """The XSPEC tapec model: APEC emission spectrum with separate continuum and line temperatures.

    The model is described at [1]_. The ``set_xsabund`` and ``get_xsabund``
    functions change and return the current settings for the relative
    abundances of the metals. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keywords "APECROOT", "APECTHERMAL", "APECVELOCITY",
    "APECNOLINES", and "APEC_TRACE_ABUND".

    Attributes
    ----------
    kT
        Continuum temperature, in keV.
    kTi
        Line temperature, in keV.
    Abundanc
        The metal abundance of the plasma, as defined by the
        ``set_xsabund`` function and the "APEC_TRACE_ABUND" xset
        keyword (which defaults to 1.0).
    Redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbtapec, XSvtapec, XSvvtapec

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTapec.html

    """

    __function__ = "C_tapec"

    def __init__(self, name='tapec'):
        self.kT = Parameter(name, 'kT', 1.0, 0.008, 64.0, 0.008, 64.0, 'keV')
        self.kTi = Parameter(name, 'kTi', 1.0, 0.008, 64.0, 0.008, 64.0, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0.0, 5.0, 0.0, 5.0)
        self.Redshift = Parameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.Abundanc, self.Redshift, self.norm))


@version_at_least("12.9.1")
class XSvtapec(XSAdditiveModel):
    """The XSPEC vtapec model: APEC emission spectrum with separate continuum and line temperatures.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keywords "APECROOT", "APECTHERMAL", "APECVELOCITY",
    "APECNOLINES", and "APEC_TRACE_ABUND".

    Attributes
    ----------
    kT
        Continuum temperature, in keV.
    kTi
        Line temperature, in keV.
    He, C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element in solar units. The trace element
        abundances are determined by the APEC_TRACE_ABUND xset
        keyword (the default is 1.0).
    Redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbvtapec, XStapec, XSvvtapec

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTapec.html

    """

    __function__ = "C_vtapec"

    def __init__(self, name='vtapec'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, 'keV')
        self.kTi = Parameter(name, 'kTi', 6.5, 0.0808, 68.447, 0.0808, 68.447, 'keV')
        self.He = Parameter(name, 'He', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.He, self.C, self.N, self.O,
                                  self.Ne, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca,
                                  self.Fe, self.Ni, self.Redshift, self.norm))


@version_at_least("12.9.1")
class XSvvtapec(XSAdditiveModel):
    """The XSPEC vvtapec model: APEC emission spectrum with separate continuum and line temperatures.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters, in
    particular the keywords "APECROOT", "APECTHERMAL", "APECVELOCITY",
    and "APECNOLINES".

    Attributes
    ----------
    kT
        Continuum temperature, in keV.
    kTi
        Line temperature, in keV.
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element in solar units.
    Redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSbvvtapec, XStapec, XSvtapec

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTapec.html

    """

    __function__ = "C_vvtapec"

    def __init__(self, name='vvtapec'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, 'keV')
        self.kTi = Parameter(name, 'kTi', 6.5, 0.0808, 68.447, 0.0808, 68.447, 'keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Redshift, self.norm))


@version_at_least("12.9.1")
class XScarbatm(XSAdditiveModel):
    """The XSPEC carbatm model: Nonmagnetic carbon atmosphere of a neutron star.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET "CARBATM"
    parameter.

    Attributes
    ----------
    T
        Effective temperature, in MK.
    NSmass
        Neutron star mass, in M_sol.
    NSrad
        Neutron star radius, in km.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XShatm

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCarbatm.html

    """

    __function__ = "C_carbatm"

    def __init__(self, name='carbatm'):
        self.T = Parameter(name, 'T', 2.0, 1.0, 4.0, 1.0, 4.0, 'MK')
        self.NSmass = Parameter(name, 'NSmass', 1.4, 0.6, 2.8, 0.6, 2.8,
                                'Msun')
        self.NSrad = Parameter(name, 'NSrad', 10.0, 6.0, 23.0, 6.0, 23.0, 'km')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.T, self.NSmass, self.NSrad, self.norm))


@version_at_least("12.9.1")
class XShatm(XSAdditiveModel):
    """The XSPEC hatm model: Nonmagnetic hydrogen atmosphere of a neutron star.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET "HATM"
    parameter.

    Attributes
    ----------
    T
        Effective temperature, in MK.
    NSmass
        Neutron star mass, in M_sol.
    NSrad
        Neutron star radius, in km.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XScarbatm

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelHatm.html

    """

    __function__ = "C_hatm"

    def __init__(self, name='hatm'):
        self.T = Parameter(name, 'T', 3.0, 0.5, 10.0, 0.5, 10.0, 'MK')
        self.NSmass = Parameter(name, 'NSmass', 1.4, 0.6, 2.8, 0.6, 2.8,
                                'Msun')
        self.NSrad = Parameter(name, 'NSrad', 10.0, 5.0, 23.0, 5.0, 23.0, 'km')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.T, self.NSmass, self.NSrad, self.norm))


@version_at_least("12.10.0")
class XSjet(XSAdditiveModel):
    """The XSPEC jet model: Leptonic relativistic jet model.

    The model is described at [1]_.

    Attributes
    ----------
    mass
        Black hole mass, in solar masses.
    Dco
        Comoving distance, as in XSoptxagnf.
    log_mdot
        log (accretion power / Eddington luminosity), in units of
        log L / logLEdd. See [1]_ for a full discussion of this parameter.
    thetaobs
        The inclination angle (deg) between the jet axis and line of sight.
    BulkG
        Bulk lorentz factor of the jet.
    phi
        The angular size scale - in radians - of the jet acceleration region
        as seen from the black hole. See [1]_ for more details.
    zdiss
        The vertical distance from the black hole of the jet dissipation
        region, in units of Rg=GM/c^2. See [1]_ for more details.
    B
        The magnetic field (in Gauss) in the jet.
    logPrel
        The log of the power injected in relativisitic particles, which is
        measured in ergs/s.
    gmin_inj
        The minimum lorentz factor of the injected electrons.
    gbreak
        The lorentz factor of the break in injected electron distribution.
    gmax
        The maximum lorentz factor.
    s1
        The injected index of the electron distribution below the break.
    s2
        The injected index of the electron distribution above the break.
    z
        The cosmological redshift corresponding to the Dco parameter.
    norm
        MUST BE FIXED AT UNITY as the jet spectrum normalisation is set
        by the relativisitic particle power.

    See Also
    --------
    XSoptxagnf

    Notes
    -----
    This model is only available when used with XSPEC 12.10.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelJet.html

    """

    __function__ = "jet"

    def __init__(self, name='jet'):
        self.mass = Parameter(name, 'mass', 1e9, 1., 1e10, 1., 1e10,
                              units='solar', frozen=True)
        self.Dco = Parameter(name, 'Dco', 3350.6, 1., 1e8, 1., 1e8,
                             units='Mpc', frozen=True)
        self.log_mdot = Parameter(name, 'log_mdot', -1., -5., 2., -5., 2.,
                                 units='logL/lEdd')
        self.thetaobs = Parameter(name, 'thetaobs', 3., 0., 90., 0., 90.,
                                  units='deg', frozen=True)
        self.BulkG = Parameter(name, 'BulkG', 13., 1., 100., 1., 100, frozen=True)
        self.phi = Parameter(name, 'phi', 0.1, 1e-2, 1e2, 1e-2, 1e2,
                             units='rad', frozen=True)
        self.zdiss = Parameter(name, 'zdiss', 1275., 10., 1e4, 10., 1e4,
                               units='Rg', frozen=True)
        self.B = Parameter(name, 'B', 2.6, 1e-2, 15., 1e-2, 15.,
                           units='Gau', frozen=True)
        self.logPrel = Parameter(name, 'logPrel', 43.3, 40., 48., 40., 48., frozen=True)
        self.gmin_inj = Parameter(name, 'gmin_inj', 1.0, 1., 1e3, 1., 1e3, frozen=True)
        self.gbreak = Parameter(name, 'gbreak', 300., 10., 1e4, 10., 1e4, frozen=True)
        self.gmax = Parameter(name, 'gmax', 3e3, 1e3, 1e6, 1e3, 1e6, frozen=True)
        self.s1 = Parameter(name, 's1', 1., -1., 1., -1., 1., frozen=True)
        self.s2 = Parameter(name, 's2', 2.7, 1., 5., 1., 5., frozen=True)
        self.z = Parameter(name, 'z', 0.0, 0., 10., 0., 10., frozen=True)
        # Note: alwaysfrozen is set for norm based on the documentation,
        # since there's no way to determine this from the model.dat file
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval,
                              alwaysfrozen=True)
        XSAdditiveModel.__init__(self, name,
                                 (self.mass, self.Dco, self.log_mdot, self.thetaobs,
                                  self.BulkG, self.phi, self.zdiss, self.B, self.logPrel,
                                  self.gmin_inj, self.gbreak, self.gmax, self.s1,
                                  self.s2, self.z, self.norm))


@version_at_least("12.9.1")
class XSismabs(XSMultiplicativeModel):
    """The XSPEC ismabs model: A high resolution ISM absorption model with variable columns for individual ions.

    The model is described at [1]_.

    Attributes
    ----------
    H
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    HeII, CI, CII, CIII, NI, NII, NIII, OI, OII, OIII, NeI, NeII, NeIII,
    MgI, MgII, MgIII, Si_I, Si_II, Si_III, S_I, S_II, S_III,
    ArI, ArII, ArIII, CaI, CaII, CaIII, Fe
        The column for the species (in units of 10^16 atoms/cm^2).
    redshift
        The redshift of the absorber.

    Notes
    -----
    As Sherpa parameter names are case insensitive the parameters
    for Silicon and Sulfur include an underscore character after the
    element name to avoid conflict: that is Si_I and S_I refer to
    the XSPEC SiI and SI parameters respectively.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelIsmabs.html

    """

    __function__ = "ismabs"

    def __init__(self, name='ismabs'):

        self.H = Parameter(name, 'H', 0.1, 0.0, 1e5, 0, 1e6, '10^22')
        self.HeII = Parameter(name, 'HeII', 0.0, 0.0, 1e5, 0, 1e6,
                              '10^16', frozen=True)
        self.CI = Parameter(name, 'CI', 33.1, 0.0, 1e5, 0, 1e6,
                            '10^16')
        self.CII = Parameter(name, 'CII', 0.0, 0.0, 1e5, 0, 1e6,
                             '10^16', frozen=True)
        self.CIII = Parameter(name, 'CIII', 0.0, 0.0, 1e5, 0, 1e6,
                              '10^16', frozen=True)
        self.NI = Parameter(name, 'NI', 8.32, 0.0, 1e5, 0, 1e6,
                            '10^16')
        self.NII = Parameter(name, 'NII', 0.0, 0.0, 1e5, 0, 1e6,
                             '10^16', frozen=True)
        self.NIII = Parameter(name, 'NIII', 0.0, 0.0, 1e5, 0, 1e6,
                              '10^16', frozen=True)
        self.OI = Parameter(name, 'OI', 67.6, 0.0, 1e5, 0, 1e6,
                            '10^16')
        self.OII = Parameter(name, 'OII', 0.0, 0.0, 1e5, 0, 1e6,
                             '10^16', frozen=True)
        self.OIII = Parameter(name, 'OIII', 0.0, 0.0, 1e5, 0, 1e6,
                              '10^16', frozen=True)
        self.NeI = Parameter(name, 'NeI', 12.0, 0.0, 1e5, 0, 1e6,
                             '10^16')
        self.NeII = Parameter(name, 'NeII', 0.0, 0.0, 1e5, 0, 1e6,
                              '10^16', frozen=True)
        self.NeIII = Parameter(name, 'NeIII', 0.0, 0.0, 1e5, 0, 1e6,
                               '10^16', frozen=True)
        self.MgI = Parameter(name, 'MgI', 3.8, 0.0, 1e5, 0, 1e6,
                             '10^16')
        self.MgII = Parameter(name, 'MgII', 0.0, 0.0, 1e5, 0, 1e6,
                              '10^16', frozen=True)
        self.MgIII = Parameter(name, 'MgIII', 0.0, 0.0, 1e5, 0, 1e6,
                               '10^16', frozen=True)
        # SiI and SI conflict, so add in underscores to differentiate.
        #
        self.Si_I = Parameter(name, 'Si_I', 3.35, 0.0, 1e5, 0, 1e6,
                              '10^16')
        self.Si_II = Parameter(name, 'Si_II', 0.0, 0.0, 1e5, 0, 1e6,
                               '10^16', frozen=True)
        self.Si_III = Parameter(name, 'Si_III', 0.0, 0.0, 1e5, 0, 1e6,
                                '10^16', frozen=True)
        self.S_I = Parameter(name, 'S_I', 2.14, 0.0, 1e5, 0, 1e6,
                             '10^16')
        self.S_II = Parameter(name, 'S_II', 0.0, 0.0, 1e5, 0, 1e6,
                              '10^16', frozen=True)
        self.S_III = Parameter(name, 'S_III', 0.0, 0.0, 1e5, 0, 1e6,
                               '10^16', frozen=True)
        self.ArI = Parameter(name, 'ArI', 0.25, 0.0, 1e5, 0, 1e6,
                             '10^16')
        self.ArII = Parameter(name, 'ArII', 0.0, 0.0, 1e5, 0, 1e6,
                              '10^16', frozen=True)
        self.ArIII = Parameter(name, 'ArIII', 0.0, 0.0, 1e5, 0, 1e6,
                               '10^16', frozen=True)
        self.CaI = Parameter(name, 'CaI', 0.22, 0.0, 1e5, 0, 1e6,
                             '10^16')
        self.CaII = Parameter(name, 'CaII', 0.0, 0.0, 1e5, 0, 1e6,
                              '10^16', frozen=True)
        self.CaIII = Parameter(name, 'CaIII', 0.0, 0.0, 1e5, 0, 1e6,
                               '10^16', frozen=True)
        self.Fe = Parameter(name, 'Fe', 3.16, 0.0, 1e5, 0, 1e6, '10^16')
        self.redshift = Parameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
                                  frozen=True)
        XSMultiplicativeModel.__init__(self, name,
                                       (self.H, self.HeII,
                                        self.CI,
                                        self.CII,
                                        self.CIII,
                                        self.NI,
                                        self.NII,
                                        self.NIII,
                                        self.OI,
                                        self.OII,
                                        self.OIII,
                                        self.NeI,
                                        self.NeII,
                                        self.NeIII,
                                        self.MgI,
                                        self.MgII,
                                        self.MgIII,
                                        self.Si_I,
                                        self.Si_II,
                                        self.Si_III,
                                        self.S_I,
                                        self.S_II,
                                        self.S_III,
                                        self.ArI,
                                        self.ArII,
                                        self.ArIII,
                                        self.CaI,
                                        self.CaII,
                                        self.CaIII,
                                        self.Fe, self.redshift))


@version_at_least("12.9.1")
class XSslimbh(XSAdditiveModel):
    """The XSPEC slimbh model: Stationary slim accretion disk.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET "SLIMBB_DEBUG",
    "SLIMBB_DIR", and "SLIMBB_TABLE" parameters.

    Attributes
    ----------
    M
        The mass of the black hole, in M_sol.
    a
        The black hole spin parameter.
    lumin
        The total disk luminosity in Eddington units.
    alpha
        alpha-viscosity
    inc
        The inclination, in degrees.
    D
        The distance to the source, in kpc.
    f_hard
        The hardening factor. A negative value means to use TLUSTY spectra.
    lflag
        A flag to control the effect of limb-darkening. If greater than
        zero then the disk emission is limb-darkened, otherwise it is assumed
        to be isotropic.
    vflag
        A flag to contorl the surface profile. If greater than zero,
        raytracing is done from the actual photosphere, otherwise the
        spectrum is raytraced from the equatorial plane ignoring the
        height profile of the disk.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSlimbh.html

    """

    __function__ = "slimbbmodel"

    def __init__(self, name='slimbh'):
        self.M = Parameter(name, 'M', 10.0, 0.0, 1000.0, 0.0, 1000.0, 'Msun',
                           frozen=True)
        self.a = Parameter(name, 'a', 0.0, 0.0, 0.999, 0.0, 0.999, 'GM/c')
        self.lumin = Parameter(name, 'lumin', 0.5, 0.05, 1.0, 0.05, 1.0,
                               'L_Edd')
        self.alpha = Parameter(name, 'alpha', 0.1, 0.005, 0.1, 0.001, 0.1,
                               frozen=True)
        self.inc = Parameter(name, 'inc', 60.0, 0.0, 85.0, 0.0, 85.0,
                             'deg', frozen=True)
        self.D = Parameter(name, 'D', 10.0, 0.0, 1e4, 0.0, 1e4, 'kpc',
                           frozen=True)
        self.f_hard = Parameter(name, 'f_hard', -1.0, -10.0, 10.0, -10.0, 10.0,
                                frozen=True)
        self.lflag = Parameter(name, 'lflag', 1.0, -1, 1, alwaysfrozen=True)
        self.vflag = Parameter(name, 'vflag', 1.0, -1, 1, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval,
                              frozen=True)
        XSAdditiveModel.__init__(self, name,
                                 (self.M, self.a, self.lumin, self.alpha,
                                  self.inc, self.D, self.f_hard,
                                  self.lflag, self.vflag, self.norm))


@version_at_least("12.9.1")
class XSsnapec(XSAdditiveModel):
    """The XSPEC snapec model: galaxy cluster spectrum using SN yields.

    The model is described at [1]_. The ``set_xsxset`` and ``get_xsxset``
    functions are used to set and query the XSPEC XSET parameters to control
    the APEC model.

    Attributes
    ----------
    kT
        The temperature of the plasma, in keV.
    N_SNe
        The number of SNe (in units of 10^9).
    R
        The percentage of SN1a.
    SNIModelIndex
        SNIa yield model: see [1]_ for more details.
    SNIIModelIndex
        SNIIa yield model: see [1]_ for more details.
    redshift
        The redshift of the plasma.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSapec

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSnapec.html

    """

    __function__ = "C_snapec"

    def __init__(self, name='snapec'):
        self.kT = Parameter(name, 'kT', 1., 0.008, 64.0, 0.008, 64.0, 'keV')
        self.N_SNe = Parameter(name, 'N_SNe', 1., 0.0, 1e20, 0, 1e20, '10^9')
        # QUS: if this is a percentage, why is the maximum not 100?
        self.R = Parameter(name, 'R', 1., 0.0, 1e20, 0, 1e20)
        self.SNIModelIndex = Parameter(name, 'SNIModelIndex', 1.,
                                       0.0, 125, 0, 125, alwaysfrozen=True)
        self.SNIIModelIndex = Parameter(name, 'SNIIModelIndex', 1.,
                                        0.0, 125, 0, 125, alwaysfrozen=True)
        self.redshift = Parameter(name, 'redshift', 0., 0.0, 10., 0.0, 10.0,
                                  frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.N_SNe, self.R,
                                  self.SNIModelIndex, self.SNIIModelIndex,
                                  self.redshift, self.norm))


@version_at_least("12.9.1")
class XSTBfeo(XSMultiplicativeModel):
    """The XSPEC TBfeo model: ISM grain absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    O
        Oxygen abundance relative to Solar.
    Fe
        Iron abundance relative to Solar.
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSTBabs, XSTBgas, XSTBgrain, XSTBpcf, XSTBrel, XSTBvarabs, XSzTBabs

    Notes
    -----
    The `set_xsabund` function changes the relative abundances of
    the elements, in particular the "wilm" setting.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbfeo"

    def __init__(self, name='tbfeo'):
        self.nH = Parameter(name, 'nH', 1., 0., 1.e5, 0.0, 1.0e6, '10^22')
        self.O = Parameter(name, 'O', 1., 0.0, 5.0, -1.0e38, 1.0e38)
        self.Fe = Parameter(name, 'Fe', 1., 0.0, 5.0, -1.0e38, 1.0e38)
        self.redshift = Parameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
                                  frozen=True)
        XSMultiplicativeModel.__init__(self, name,
                                       (self.nH, self.O, self.Fe,
                                        self.redshift))


@version_at_least("12.9.1")
class XSTBgas(XSMultiplicativeModel):
    """The XSPEC TBgas model: ISM grain absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSTBabs, XSTBfeo, XSTBgrain, XSTBpcf, XSTBrel, XSTBvarabs, XSzTBabs

    Notes
    -----
    The `set_xsabund` function changes the relative abundances of
    the elements, in particular the "wilm" setting.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbgas"

    def __init__(self, name='tbgas'):
        self.nH = Parameter(name, 'nH', 1., 0., 1.e5, 0.0, 1.0e6, '10^22')
        self.redshift = Parameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
                                  frozen=True)
        XSMultiplicativeModel.__init__(self, name,
                                       (self.nH, self.redshift))


@version_at_least("12.9.1")
class XSTBpcf(XSMultiplicativeModel):
    """The XSPEC TBpcf model: ISM grain absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    pcf
        Partial covering fraction of the absorber.
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSTBabs, XSTBfeo, XSTBgas, XSTBgrain, XSTBrel, XSTBvarabs, XSzTBabs

    Notes
    -----
    The `set_xsabund` function changes the relative abundances of
    the elements, in particular the "wilm" setting.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbpcf"

    def __init__(self, name='tbpcf'):
        self.nH = Parameter(name, 'nH', 1., 0., 1.e5, 0.0, 1.0e6, '10^22')
        self.pcf = Parameter(name, 'pcf', 0.5, 0, 1.0, 0, 1.0)
        self.redshift = Parameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
                                  frozen=True)
        XSMultiplicativeModel.__init__(self, name,
                                       (self.nH, self.pcf, self.redshift))


@version_at_least("12.9.1")
class XSTBrel(XSMultiplicativeModel):
    """The XSPEC TBrel model: ISM grain absorption.

    The model is described at [1]_.

    Attributes
    ----------
    nH
        The equivalent hydrogen column (in units of 10^22 atoms/cm^2).
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Cl, Ar, Ca, Cr, Fe, Co, Ni
        The abundance of the element in solar units.
    H2
        The equivalent molecular hydrogen column (in units of
         10^22 atoms/cm^2).
    rho
        The grain density, in g/cm^3.
    amin
        The minimum grain size, in micro-meters.
    amax
        The maximum grain size, in micro-meters.
    PL
        The power-law index of grain sizes.
    H_dep, He_dep, C_dep, N_dep, O_dep, Ne_dep, Na_dep, Mg_dep, Al_dep,
    Si_dep, S_dep, Cl_dep, Ar_dep, Ca_dep, Cr_dep, Fe_dep, Co_dep, Ni_dep
        The grain depletion fraction of the element.
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSTBabs, XSTBfeo, XSTBgas, XSTBgrain, XSTBpcf, XSTBvarabs, XSzTBabs

    Notes
    -----
    The `set_xsabund` function changes the relative abundances of
    the elements, in particular the "wilm" setting.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbrel"

    def __init__(self, name='tbrel'):
        self.nH = Parameter(name, 'nH', 1., -1e5, 1e5, -1e6, 1.0e6, '10^22')
        self.He = Parameter(name, 'He', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.H2 = Parameter(name, 'H2', 0.2, 0., 1., 0.0, 1.0, frozen=True)
        self.rho = Parameter(name, 'rho', 1., 0., 5., 0.0, 5.0, 'g/cm^3',
                             frozen=True)
        self.amin = Parameter(name, 'amin', 0.025, 0., 0.25, 0.0, 0.25,
                              'mum', frozen=True)
        self.amax = Parameter(name, 'amax', 0.25, 0., 1., 0.0, 1.0,
                              'mum', frozen=True)
        self.PL = Parameter(name, 'PL', 3.5, 0., 5., 0.0, 5.0, frozen=True)
        self.H_dep = Parameter(name, 'H_dep', 1., 0., 1., 0.0, 1.0,
                               frozen=True)
        self.He_dep = Parameter(name, 'He_dep', 1., 0., 1., 0.0, 1.0,
                                frozen=True)
        self.C_dep = Parameter(name, 'C_dep', 0.5, 0., 1., 0.0, 1.0,
                               frozen=True)
        self.N_dep = Parameter(name, 'N_dep', 1., 0., 1., 0.0, 1.0,
                               frozen=True)
        self.O_dep = Parameter(name, 'O_dep', 0.6, 0., 1., 0.0, 1.0,
                               frozen=True)
        self.Ne_dep = Parameter(name, 'Ne_dep', 1., 0., 1., 0.0, 1.0,
                                frozen=True)
        self.Na_dep = Parameter(name, 'Na_dep', 0.25, 0., 1., 0.0, 1.0,
                                frozen=True)
        self.Mg_dep = Parameter(name, 'Mg_dep', 0.2, 0., 1., 0.0, 1.0,
                                frozen=True)
        self.Al_dep = Parameter(name, 'Al_dep', 0.02, 0., 1., 0.0, 1.0,
                                frozen=True)
        self.Si_dep = Parameter(name, 'Si_dep', 0.1, 0., 1., 0.0, 1.0,
                                frozen=True)
        self.S_dep = Parameter(name, 'S_dep', 0.6, 0., 1., 0.0, 1.0,
                               frozen=True)
        self.Cl_dep = Parameter(name, 'Cl_dep', 0.5, 0., 1., 0.0, 1.0,
                                frozen=True)
        self.Ar_dep = Parameter(name, 'Ar_dep', 1., 0., 1., 0.0, 1.0,
                                frozen=True)
        self.Ca_dep = Parameter(name, 'Ca_dep', 0.003, 0., 1., 0.0, 1.0,
                                frozen=True)
        self.Cr_dep = Parameter(name, 'Cr_dep', 0.03, 0., 1., 0.0, 1.0,
                                frozen=True)
        self.Fe_dep = Parameter(name, 'Fe_dep', 0.3, 0., 1., 0.0, 1.0,
                                frozen=True)
        self.Co_dep = Parameter(name, 'Co_dep', 0.05, 0., 1., 0.0, 1.0,
                                frozen=True)
        self.Ni_dep = Parameter(name, 'Ni_dep', 0.04, 0., 1., 0.0, 1.0,
                                frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
                                  frozen=True)
        pars = (self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na,
                self.Mg, self.Al, self.Si, self.S, self.Cl, self.Ar, self.Ca,
                self.Cr, self.Fe, self.Co, self.Ni, self.H2, self.rho,
                self.amin, self.amax, self.PL, self.H_dep, self.He_dep,
                self.C_dep, self.N_dep, self.O_dep, self.Ne_dep, self.Na_dep,
                self.Mg_dep, self.Al_dep, self.Si_dep, self.S_dep, self.Cl_dep,
                self.Ar_dep, self.Ca_dep, self.Cr_dep, self.Fe_dep,
                self.Co_dep, self.Ni_dep, self.redshift)
        XSMultiplicativeModel.__init__(self, name, pars)


@version_at_least("12.9.1")
class XSvoigt(XSAdditiveModel):
    """The XSPEC voigt model: Voigt line profile.

    The model is described at [1]_.

    Attributes
    ----------
    LineE
        The line energy, in keV.
    Sigma
        The gaussian line width, in keV.
    Gamma
        The lorentzian FWHM line width, in keV.
    norm
        The flux in the line, in units of photon/cm^2/s.

    See Also
    --------
    XSgauss, XSlorentz

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelVoigt.html

    """

    __function__ = "C_voigtLine"

    def __init__(self, name='voigt'):
        self.LineE = Parameter(name, 'LineE', 6.5, 0., 1.e6, 0.0, hugeval, 'keV')
        self.Sigma = Parameter(name, 'Sigma', 0.01, 0., 10., 0.0, hugeval, 'keV')
        self.Gamma = Parameter(name, 'Gamma', 0.01, 0., 10., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.LineE, self.Sigma, self.Gamma,
                                  self.norm))


@version_at_least("12.9.1")
class XSxscat(XSMultiplicativeModel):
    """The XSPEC xscat model: dust scattering.

    The model is described at [1]_.

    Attributes
    ----------
    NH
        Interstellar hydrogen column density (in units of 10^22 atoms/cm^2).
    Xpos
        The relative position of the dust along the line of sight. A value
        of 0 means the dust is at the observer, and of 0.95 means that the
        dust is 5% from the source.
    Rext
        The radius of the circular extraction region, in arcsec.
    DustModel
        The dust model used: see [1]_ for more information.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelXscat.html

    """

    __function__ = "C_xscatmodel"

    def __init__(self, name='xscat'):
        self.NH = Parameter(name, 'NH', 1., 0., 1000.0, 0.0, 1000.0, '10^22')
        self.Xpos = Parameter(name, 'Xpos', 0.5, 0, 0.95, 0, 0.95)
        self.Rext = Parameter(name, 'Rext', 10.0, 0, 115.0, 0, 119.0, 'arcsec',
                              frozen=True)
        # The maxmimum number of models depends on the data file, so pick
        # a value that is unlikely to be exceeded (the max at the time
        # of writing was 3).
        self.DustModel = Parameter(name, 'DustModel', 1, 1, 100, 1, 100,
                                   alwaysfrozen=True)
        XSMultiplicativeModel.__init__(self, name,
                                       (self.NH, self.Xpos, self.Rext,
                                        self.DustModel))


# Add model classes to __all__
__all__ += tuple(n for n in globals() if n.startswith('XS'))
