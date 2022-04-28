#
#  Copyright (C) 2010, 2015-2018, 2019, 2020, 2021, 2022
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

"""Support for XSPEC models.

Sherpa supports versions 12.12.1, 12.12.0, 12.11.1, 12.11.0, 12.10.1, 12.10.0, 12.9.1, and 12.9.0
of XSPEC [1]_, and can be built against the model library or the full
application.  There is no guarantee of support for older or newer
versions of XSPEC.

To be able to use most routines from this module, the HEADAS environment
variable must be set. The `get_xsversion` function can be used to return the
XSPEC version - including patch level - the module is using::

   >>> from sherpa.astro import xspec
   >>> xspec.get_xsversion()
   '12.11.1'

Initializing XSPEC
------------------

The XSPEC model library is initalized so that the cosmology parameters
are set to H_0=70, q_0=0.0, and lambda_0=0.73 (they can be changed with
`set_xscosmo`).

The other settings - for example for the abundance and cross-section
tables - follow the standard rules for XSPEC. For XSPEC versions prior
to 12.10.1, this means that the abundance table uses the 'angr'
setting and the cross sections the 'bcmc' setting (see `set_xsabund`
and `set_xsxsect` for full details). As of XSPEC 12.10.1, the values
are now taken from the user's XSPEC configuration file - either
``~/.xspec/Xspec.init`` or ``$HEADAS/../spectral/manager/Xspec.init`` -
for these settings. The default value for the photo-ionization table
in this case is now 'vern' rather than 'bcmc'.

The default chatter setting - used by models to inform users of
issues - was set to 0 (which hid the messages) until Sherpa 4.14.0,
when it was changed to 10 (to match XSPEC).

Supported models
----------------

The additive [2]_, multiplicative [3]_, and convolution [4]_ models
from the XSPEC model library are supported, except for the `smaug`
model [5]_, since it requires use of information from the XFLT keywords
in the data file).

Parameter values
----------------

XSPEC parameters have soft and hard limits but they are different
to the Sherpa meaning:

- the XSPEC hard limit is more-like a Sherpa parameter with the
  same soft and hard limits;

- and it is possible to change the hard limits.

To support these, XSPEC models use the `XSBaseParameter` and
`XSParameter` classes, while the "norm" parameter added to additive
models remains a `sherpa.models.parameter.Parameter` instance.

References
----------

.. [1] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/index.html

.. [2] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/Additive.html

.. [3] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/Multiplicative.html

.. [4] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/Convolution.html

.. [5] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSmodelSmaug.html

"""


import logging
import string
import warnings

import numpy as np

from sherpa.models import ArithmeticModel, ArithmeticFunctionModel, \
    CompositeModel, Parameter, modelCacher1d, RegriddableModel1D
from sherpa.models.parameter import hugeval

from sherpa.utils import SherpaFloat, guess_amplitude, param_apply_limits, bool_cast
from sherpa.utils.err import ParameterErr
from sherpa.astro.utils import get_xspec_position

from .utils import ModelMeta, version_at_least, equal_or_greater_than
from . import _xspec


info = logging.getLogger(__name__).info
warning = logging.getLogger(__name__).warning

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
       When `element` is `None`, the abundance table name is
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
    10
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

    .. [1] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/

    Examples
    --------

    >>> get_xsversion()
    '12.11.0m'
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
     - 'lpgp', from [9]_ (photospheric, requires XSPEC 12.12.0 or later)
     - 'lpgs', from [9]_ (proto-solar, requires XSPEC 12.12.0 or later)

    The values for these tables are given at [1]_.

    Data files should be in ASCII format, containing a single
    numeric (floating-point) column of the abundance values,
    relative to Hydrogen.

    The screen output of this function is controlled by the
    X-Spec chatter setting (`set_xschatter`).

    References
    ----------
    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSabund.html
           Note that this may refer to a newer version than the
           compiled version used by Sherpa; use `get_xsversion` to
           check.

    .. [2] Anders E. & Grevesse N. (1989, Geochimica et
           Cosmochimica Acta 53, 197)
           https://adsabs.harvard.edu/abs/1989GeCoA..53..197A

    .. [3] Asplund M., Grevesse N., Sauval A.J. & Scott P.
           (2009, ARAA, 47, 481)
           https://adsabs.harvard.edu/abs/2009ARA%26A..47..481A

    .. [4] Feldman U.(1992, Physica Scripta 46, 202)
           https://adsabs.harvard.edu/abs/1992PhyS...46..202F

    .. [5] Anders E. & Ebihara (1982, Geochimica et Cosmochimica
           Acta 46, 2363)
           https://adsabs.harvard.edu/abs/1982GeCoA..46.2363A

    .. [6] Grevesse, N. & Sauval, A.J. (1998, Space Science
           Reviews 85, 161)
           https://adsabs.harvard.edu/abs/1998SSRv...85..161G

    .. [7] Wilms, Allen & McCray (2000, ApJ 542, 914)
           https://adsabs.harvard.edu/abs/2000ApJ...542..914W

    .. [8] Lodders, K (2003, ApJ 591, 1220)
           https://adsabs.harvard.edu/abs/2003ApJ...591.1220L

    .. [9] Lodders K., Palme H., Gail H.P., Landolt-Börnstein,
           New Series, vol VI/4B, pp 560–630 (2009)
           https://ui.adsabs.harvard.edu/abs/2009LanB...4B..712L/abstract

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

    .. versionchanged:: 4.14.0
       The default chatter setting has been bumped from 0 to 10 to
       match XSPEC. Users will see extra screen output the first
       time some XSPEC models are evaluated.

    Parameters
    ----------
    level : int
       The higher the value of ``level``, the more screen output will
       be created by X-Spec routines. A value of ``0`` hides most
       information while ``25`` will generate a lot of debug
       output. The starting value is ``10``.

    See Also
    --------
    get_xschatter, get_xsversion

    Notes
    -----
    There is no way to change the X-Spec "log chatter" setting.

    References
    ----------
    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSchatter.html
           Note that this may refer to a newer version than the
           compiled version used by Sherpa; use `get_xsversion` to
           check.

    Examples
    --------

    Set the chatter level to hide most, if not all, output from X-Spec
    models:

    >>> set_xschatter(0)

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
    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XScosmo.html
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
    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSxsect.html
           Note that this may refer to a newer version than the
           compiled version used by Sherpa; use `get_xsversion` to
           check.

    .. [2] Balucinska-Church & McCammon (1992; Ap.J.400, 699).
           https://adsabs.harvard.edu/abs/1992ApJ...400..699B

    .. [3] Yan, M., Sadeghpour, H. R., Dalgarno, A. 1998,
           Ap.J. 496, 1044
           https://adsabs.harvard.edu/abs/1998ApJ...496.1044Y

    .. [4] Verner et. al., 1996, Ap.J., 465, 487.
           https://adsabs.harvard.edu/abs/1996ApJ...465..487V

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
# https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSxset.html
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
    output of `sherpa.astro.ui.save_all`.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSabund.html
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
# https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/Control.html
# https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/Setting.html
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


def read_xstable_model(modelname, filename, etable=False):
    """Create a XSPEC table model.

    XSPEC additive (atable, [1]_), multiplicative (mtable, [2]_), and
    exponential (etable, [3]_) table models are supported.

    .. versionchanged:: 4.14.0
       The etable argument has been added to allow exponential table
       models to be used.

    Parameters
    ----------
    modelname : str
       The identifier for this model component.
    filename : str
       The name of the FITS file containing the data, which should
       match the XSPEC table model definition [4]_.
    etable : bool, optional
       Set if this is an etable (as there's no way to determine this
       from the file itself). Defaults to False.

    Returns
    -------
    tablemodel : XSTableModel instance

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSmodelAtable.html

    .. [2] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSmodelMtable.html

    .. [3] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSmodelMtable.html

    .. [4] https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html

    Examples
    --------

    Load in the XSPEC table model from the file 'bbrefl_1xsolar.fits'
    and create a model labelled 'xmdl', which is then returned:

    >>> mdl = read_xstable_model('xmdl', 'bbrefl_1xsolar.fits')
    >>> print(mdl)

    When the file is an etable model the etable parameter must be set:

    >>> emdl = read_xstable_model('xmdl', 'etable.fits', etable=True)
    >>> print(emdl)

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
                        addredshift=addredshift, etable=etable)


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


class XSBaseParameter(Parameter):
    """An XSPEC parameter.

    XSPEC has soft and hard parameter limits, which are the ones sent
    in as the `min`, `max`, `hard_min`, and `hard_max` parameters.
    However, Sherpa's soft limits are more-like the XSPEC hard limits,
    and it is possible in XSPEC to change a model's hard limits. This
    class therefore:

    - stores the input `min` and `max` values as the
      _xspec_soft_min and _xspec_soft_max attributes and the
      `hard_min` and `hard_max` values as _xspec_hard_min and
      _xspec_hard_max attributes;

    - sets the underlying `min` and `max` values to the XSPEC hard
      limits;

    - and sets the underlying `hard_min` and `hard_max` values to
      the XSPEC hard limits.

    Note that you can not change the hard limits; for that see
    `XSParameter`.

    See Also
    --------
    XSParameter

    Examples
    --------

    >>> p = XSBaseParameter('mod', 'p', 2, min=1, max=9, hard_min=0, hard_max=10)
    >>> p.min
    0.0
    >>> p.hard_min
    0.0
    >>> p.max
    10.0
    >>> p.hard_max
    10.0
    >>> p.val = 20
    sherpa.utils.err.ParameterErr: parameter mod.p has a maximum of 10

    """

    def __init__(self, modelname, name, val, min=-hugeval, max=hugeval,
                 hard_min=-hugeval, hard_max=hugeval, units='',
                 frozen=False, alwaysfrozen=False, hidden=False, aliases=None):

        self._xspec_soft_min = min
        self._xspec_soft_max = max
        self._xspec_hard_min = hard_min
        self._xspec_hard_max = hard_max
        super().__init__(modelname, name, val, min=hard_min, max=hard_max,
                         hard_min=hard_min, hard_max=hard_max, units=units,
                         frozen=frozen, alwaysfrozen=alwaysfrozen,
                         hidden=hidden, aliases=aliases)


class XSParameter(XSBaseParameter):
    """An XSPEC parameter where you exceed the hard limits.

    This parameter allows a user to change the hard limits (`hard_min`
    and `hard_max`). XSPEC allows the hard limits to be extended, and
    this is used in a few models to trigger different behavior (e.g.
    setting the value negative). The `hard_min_changed` and
    `hard_max_changed` methods can be used to determine if the limits
    have been changed.

    See Also
    --------
    XSBaseParameter

    Notes
    -----
    Some XSPEC parameter values are documented as changing behavior
    when the value is outside the XSPEC hard limits from the model.dat
    file (normally the hard minimum is 0 and setting the parameter
    negative changes the model in some way). XSPEC allows a user to
    change the hard limits so we do the same here.

    Setting the hard limit will also change the corresponding soft limit.

    Setting a value outside the original hard limits can cause models
    to fail or even crash the interpreter. These parameters should
    probably be frozen.

    Examples
    --------

    >>> p = XSParameter('mod', 'p', 2, min=1, max=9, hard_min=0, hard_max=10)
    >>> p.frozen
    False
    >>> p.min
    0.0
    >>> p.hard_min
    0.0
    >>> p.max
    10.0
    >>> p.hard_max
    10.0
    >>> p.val = 20
    sherpa.utils.err.ParameterErr: parameter mod.p has a maximum of 10
    >>> p.max = 30
    sherpa.utils.err.ParameterErr: parameter mod.p has a hard maximum of 10
    >>> p.hard_max = 30
    >>> p.max
    30.0
    >>> p.hard_max
    30.0
    >>> p.val = 20
    >>> p.frozen
    False

    """

    def _set_hard_min(self, val):

        # Ensure we are not selecting a value that doesn't make sense.
        #
        val = SherpaFloat(val)
        if val >= self.max:
            raise ParameterErr('edge', self.fullname,
                               'maximum', self.max)

        # If we have increased the minimum value so that the current
        # value is now excluded we need to update it.
        #
        if self.val < val:
            self.val = val
            warning(f'parameter {self.fullname} less than new minimum; reset to {self.val}')

        # Set both soft and hard limits. We use _min so we don't
        # trigger any validation of the value (e.g. when the new
        # minimum is larger than the old minimum).
        #
        self._min = val
        self._hard_min = val

    hard_min = property(Parameter._get_hard_min, _set_hard_min,
                        doc='The hard minimum of the parameter.\n\n' +
                        'Unlike normal parameters the `hard_min` value can be changed (and\n' +
                        'will also change the corresponding `min` value at the same time).\n' +
                        'This is needed to support the small-number of XSPEC models that\n' +
                        'use a value outside the default hard range as a way to control the\n' +
                        'model. Unfortunately some models can crash when using values like\n' +
                        'this so take care.\n\n' +
                        'See Also\n' +
                        '--------\n' +
                        'hard_max\n')

    def _set_hard_max(self, val):

        # Ensure we are not selecting a value that doesn't make sense.
        #
        val = SherpaFloat(val)
        if val <= self.min:
            raise ParameterErr('edge', self.fullname,
                               'minimum', self.min)

        # If we have decreased the maximum value so that the current
        # value is now excluded we need to update it.
        #
        if self.val > val:
            self.val = val
            warning(f'parameter {self.fullname} greater than new maximum; reset to {self.val}')

        # Set both soft and hard limits. We use _max so we don't
        # trigger any validation of the value (e.g. when the new
        # maximum is smaller than the old maximum).
        #
        self._max = val
        self._hard_max = val

    hard_max = property(Parameter._get_hard_max, _set_hard_max,
                        doc='The hard maximum of the parameter.\n\n' +
                        'Unlike normal parameters the `hard_max` value can be changed (and\n' +
                        'will also change the corresponding `max` value at the same time).\n' +
                        'This is needed to support the small-number of XSPEC models that\n' +
                        'use a value outside the default hard range as a way to control the\n' +
                        'model. Unfortunately some models can crash when using values like\n' +
                        'this so take care.\n\n' +
                        'See Also\n' +
                        '--------\n' +
                        'hard_min\n')

    def set(self, val=None, min=None, max=None, frozen=None,
            default_val=None, default_min=None, default_max=None,
            hard_min=None, hard_max=None):
        """Change a parameter setting.

        The hard limits can be changed, which will also change the
        matching soft limit. Note that XSPEC models can cause a crash
        if sent an un-supported value so use this feature carefully;
        it is likely that these parameters should also be frozen but
        this is not enforced.

        Parameters
        ----------
        val : number or None, optional
            The new parameter value.
        min, max : number or None, optional
            The new parameter range.
        frozen : bool or None, optional
            Should the frozen flag be set?
        default_val : number or None, optional
            The new default parameter value.
        default_min, default_max : number or None, optional
            The new default parameter limits.
        hard_min, hard_max : numer or None, optional
            Changing the hard limits will also change the matching
            soft limit (`min` or `max`).

        """

        if hard_min is not None:
            self.hard_min = hard_min

        if hard_max is not None:
            self.hard_max = hard_max

        super().set(val=val, min=min, max=max, frozen=frozen,
                    default_val=default_val,
                    default_min=default_min,
                    default_max=default_max)

    def hard_min_changed(self):
        """Has the hard limit (min) been changed from it's default value?"""
        return self._xspec_hard_min != self.hard_min

    def hard_max_changed(self):
        """Has the hard limit (max) been changed from it's default value?"""
        return self._xspec_hard_max != self.hard_max


class XSModel(RegriddableModel1D, metaclass=ModelMeta):
    """The base class for XSPEC models.

    It is expected that sub-classes are used to represent the
    five different types of XSPEC model (additive, multiplicative,
    convolution, pile up, mixing, and tables), although not all are
    currently supported in Sherpa.

    Notes
    -----
    The XSPEC models are evaluated on a one-dimensional, integrated,
    contiguous grid. When the `calc` method is called with both
    low and high bin values, the arrays are converted into a single
    array matching the XSPEC calling convention - that is elo_0,
    elo_1, ..., elo_n for n bins (so the last value is the upper edge
    of the last bin) - adding in any bins to account for a non-contiguous
    input. This array is used to evaluate the model, and then the
    return value is created by removing any extra bins that had to be
    added to account for non-contiguous input values.

    If used on an unbinned dataset, so only one array is sent to
    `calc`, then the input values are taken to match the XSPEC
    calling convention - i.e. a contiguous grid where the last element
    represents the upper edge of the last bin. This means that for
    an input grid of ``n`` points, the returned array will contain
    ``n`` values, but the last element will be zero.
    """

    version_enabled = True

    def __init__(self, name, pars):

        # Validate the parameters argument to check that the
        # names are unique (name insensitive) and members of
        # the class (i.e. there are attributes with this name).
        #
        # This could be done at a higher level of the class structure
        # but given the support for composite models it can't be
        # done at the Model layer. It probably makes sense to do this
        # at the ArithmeticModel layer but this is the place where we
        # make the most changes and we know these constraints are
        # correct here so do so here.
        #
        # These are development errors so use asserts rather than
        # raising an error.
        #
        seen = set()
        for par in pars:
            pname = par.name.lower()

            # unique names
            assert pname not in seen, (par.name, name)
            seen.add(pname)

            # there's actually an attribute with this name
            attr = getattr(self, par.name)
            assert attr.name == par.name, (par.name, name)

        super().__init__(name, pars)


    @modelCacher1d
    def calc(self, *args, **kwargs):
        """Calculate the model given the parameters and grid.

        Notes
        -----
        XSPEC models must always be evaluated with low and high bin
        edges. Although supported by the XSPEC model interface the
        ability to evaluate using an XSPEC-style grid (n+1 values for
        n bins which we pad with a 0), we do not allow this here since
        it complicates the handling of the regrid method.

        Keyword arguments are ignored.
        """

        nargs = len(args)
        if nargs != 3:
            emsg = f"calc() requires pars,lo,hi arguments, sent {nargs} arguments"
            warnings.warn(emsg, FutureWarning)
            # raise TypeError(emsg)

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
    by `read_xstable_model` and `sherpa.astro.ui.load_xstable_model`.

    .. versionchanged:: 4.14.0
       The etable argument has been added to allow exponential table
       models to be used.

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
        Is this an additive model (`True`) or multiplicative model
        (`False`)? It should be set to the value of the "ADDMODEL"
        keyword of the primary header of the input file. When False
        the etable keyword is used to distinguish between mtable and
        etable models.
    addredshift : bool
        If `True` then a redshift parameter is added to the parameters.
        It should be set to the value of the "REDSHIFT" keyword of the
        primary header of the input file.
    etable : bool
        When addmodel is False this defines whether the file is a
        mtable model (`False`, the default) or an etable model (`True`).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html
    .. [2] https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html

    """

    def __init__(self, filename, name='xstbl', parnames=(),
                 initvals=(), delta=(), mins=(), maxes=(), hardmins=(),
                 hardmaxes=(), nint=0,
                 addmodel=False, addredshift=False, etable=False):

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

            # For table models we know you can crash the system if you change
            # the hard limits, so we use XSBaseParameter rather than
            # XSParameter.
            #
            parname = parname.strip().lower().translate(tbl)
            par = XSBaseParameter(name, parname, initvals[ii],
                                  mins[ii], maxes[ii],
                                  hardmins[ii], hardmaxes[ii], frozen=isfrozen)
            self.__dict__[parname] = par
            pars.append(par)
            nint -= 1

        self.filename = filename
        self.addmodel = addmodel
        self.etable = etable

        if addredshift:
            self.redshift = XSBaseParameter(name, 'redshift', 0., 0., 5.,
                                            0.0, 5, frozen=True)
            pars.append(self.redshift)

        if addmodel:
            # Normalization parameters are not true XSPEC parameters and
            # so we do not need to use XSParameter.
            #
            self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0,
                                  hugeval)
            pars.append(self.norm)

        XSModel.__init__(self, name, pars)

    def fold(*args, **kwargs):
        pass

    @modelCacher1d
    def calc(self, p, *args, **kwargs):

        # Note the kwargs is ignored
        nargs = 1 + len(args)
        if nargs != 3:
            emsg = f"calc() requires pars,lo,hi arguments, sent {nargs} arguments"
            warnings.warn(emsg, FutureWarning)
            # raise TypeError(emsg)

        # The function used depends on XSPEC version and, prior
        # to XSPEC 12.10.1, the type of table.
        #
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
            tabtype = 'add' if self.addmodel else 'exp' if self.etable else 'mul'
            return _xspec.tabint(p, *args,
                                 filename=self.filename, tabtype=tabtype)

        # Technicaly this should be updated to add support for etable models
        # but this interface is no-longer supported so there's no way for
        # me to add the necessary code.
        #
        if self.addmodel:
            func = _xspec.xsatbl
        else:
            func = _xspec.xsmtbl

        return func(p, *args, filename=self.filename)


# TODO: we should add the norm parameter in the __init__ call, unless
# the object contains a norm attribute, to avoid having to repeat this
# for every additive model.
#
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

    .. [1] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/Multiplivative.html

    """

    pass


class XSConvolutionKernel(XSModel):
    """The base class for XSPEC convolution models.

    The XSPEC convolution models are listed at [1]_.

    .. versionadded:: 4.12.2

    Notes
    -----

    As these models are applied to the result of other models, this
    model class isn't called directly, but creates a wrapper instance
    (`XSConvolutionModel`) that is. This wrapping is done by applying
    the kernel to the model expression using normal function
    application, that is the following creates a model called ``mdl``
    which aplies the comvolution model (in this case `XScflux`) to the
    model expression (an absorbed APEC model):

    >>> cmdl = XScflux()
    >>> src = XSapec()
    >>> gal = XSphabs()
    >>> mdl = cmdl(gal  * src)

    These expressions can then be nested, so for instance you can
    apply a convolution model to only part of the model expression, as
    in:

    >>> mdl = gal * cmdl(src)

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/Convolution.html

    Examples
    --------

    Create an instance of the XSPEC cflux convolution model, and set
    its parameters:

    >>> from sherpa.astro import xspec
    >>> cmdl = xspec.XScflux()
    >>> cmdl.emin = 0.5
    >>> cmdl.emax = 7.0
    >>> cmdl.lg10flux = -12.1
    >>> print(cmdl)
    xscflux
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       xscflux.Emin frozen          0.5            0        1e+06        keV
       xscflux.Emax frozen           10            0        1e+06        keV
       xscflux.lg10Flux thawed          -12         -100          100        cgs

    The convolution models are not evaluated directly. Instead they
    are applied to a model expression which you specify as the
    argument to the convolved model. The following creates two
    different models: mdl1 applies the convolution to the
    full model expression (absorbed powerlaw), and mdl
    applies the convolution to the powerlaw component and then
    multiplies this by the absorption model.

    >>> gal = xspec.XSphabs()
    >>> pl = xspec.XSpowerlaw()
    >>> mdl1 = cmdl(gal * pl)
    >>> mdl2 = gal * cmdl(pl)

    For the XScflux use case it is important to freeze the normalization
    parameter of the emission model (i.e. additive component) to 1,
    to ensure the lg10Flux value is calculated correctly. That is:

    >>> pl.norm = 1
    >>> pl.norm.freeze()

    """

    def __repr__(self):
        return "<{} kernel instance '{}'>".format(type(self).__name__,
                                                  self.name)

    def __call__(self, model):
        return XSConvolutionModel(model, self)

    def calc(self, pars, rhs, *args, **kwargs):
        """Evaluate the convolved model.

        Note that this method is not cached by
        sherpa.models.modelCacher1d (may change in the future).

        Parameters
        ----------
        pars : sequence of numbers
            The parameters of the convolved model. The first npars
            parameters (where npars is the lenth of the objecs pars
            attribute) are applied to the convolution model, and the
            remaining are passed to the rhs model.
        rhs : sherpa.models.model.ArithmeticModel
            The model that is being convolved.
        *args
            The model grid. There should be two arrays (the low and
            high edges of the bin) to make sure the wrapped model is
            evaluated correctly.
        **kwargs
            At present all additional keyword arguments are dropped.

        """

        nargs = 2 + len(args)
        if nargs != 4:
            emsg = f"calc() requires pars,rhs,lo,hi arguments, sent {nargs} arguments"
            warnings.warn(emsg, FutureWarning)
            # raise TypeError(emsg)

        npars = len(self.pars)
        lpars = pars[:npars]
        rpars = pars[npars:]

        # We do not pass kwargs on to either rhs or self._calc;
        # this is needed because when used with regrid support the
        # kwargs may include keywords like 'integrate' but this
        # can then cause problems downstream (that is, they then
        # get sent to the xspec compiled routine, which doesn't
        # support them). This may be a problem with the XSPEC
        # interface (i.e. it should not pass on kwargs to the compiled
        # code), but for now stop it here.
        #
        # DJB is worried that if there is a nested set of grids
        # then this could be a problem, but has no experience or
        # tests to back this up.
        #
        # fluxes = np.asarray(rhs(rpars, *args, **kwargs))
        fluxes = np.asarray(rhs(rpars, *args))
        return self._calc(lpars, fluxes, *args)


# It makes sense to derive from XSModel to say "this is XSPEC",
# but does the extra machinery it provides cause a problem here?
#
class XSConvolutionModel(CompositeModel, XSModel):
    """Evaluate a model and pass it to an XSPEC convolution model.

    Calculate the convolved data - that is evaluate the model and then
    pass it to the wrapper model which applies the convolution model.

    .. versionadded:: 4.12.2

    Parameters
    ----------
    model : `sherpa.models.model.ArithmeticModel`
        The model whose results, when evaluated, are passed to
        the convolution model.
    wrapper : `sherpa.astro.xspec.XSConvolutionKernel`
        The XSPEC convolution model.

    Examples
    --------

    The following evaluates two models (creating the y1 and y2
    arrays), where y1 applies the `XScflux` convolution model to the
    combined absorption times powerlaw model, and y2 applies the
    convolution model to only the power-law model, and then multiples
    this by the absorption model. In the following mdl1 and mdl2
    are instances of XSConvolutionModel:

    >>> import numpy as np
    >>> from sherpa.astro import xspec
    >>> cmdl = xspec.XScflux()
    >>> gal = xspec.XSphabs()
    >>> pl = xspec.XSpowerlaw()
    >>> pl.norm.freeze()
    >>> mdl1 = cmdl(gal * pl)
    >>> mdl2 = gal * cmdl(pl)
    >>> cmdl.emin = 0.5
    >>> cmdl.emax = 7.0
    >>> cmdl.lg10flux = -12.1
    >>> egrid = np.arange(0.1, 10, 0.01)
    >>> elo, ehi = egrid[:-1], egrid[1:]
    >>> y1 = mdl1(elo, ehi)
    >>> y2 = mdl2(elo, ehi)

    Display the combined model:

    >>> print(mdl1)
    xscflux((phabs * powerlaw))
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       xscflux.Emin frozen          0.5            0        1e+06        keV
       xscflux.Emax frozen           10            0        1e+06        keV
       xscflux.lg10Flux thawed          -12         -100          100        cgs
       phabs.nH     thawed            1            0       100000 10^22 atoms / cm^2
       powerlaw.PhoIndex thawed            1           -2            9
       powerlaw.norm frozen            1            0        1e+24

    >>> print(mdl2)
    (phabs * xscflux(powerlaw))
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       phabs.nH     thawed            1            0       100000 10^22 atoms / cm^2
       xscflux.Emin frozen          0.5            0        1e+06        keV
       xscflux.Emax frozen           10            0        1e+06        keV
       xscflux.lg10Flux thawed          -12         -100          100        cgs
       powerlaw.PhoIndex thawed            1           -2            9
       powerlaw.norm frozen            1            0        1e+24

    """

    @staticmethod
    def wrapobj(obj):
        if isinstance(obj, ArithmeticModel):
            return obj
        else:
            return ArithmeticFunctionModel(obj)

    def __init__(self, model, wrapper):
        self.model = self.wrapobj(model)
        self.wrapper = wrapper
        CompositeModel.__init__(self,
                                "{}({})".format(self.wrapper.name,
                                                self.model.name),
                                (self.wrapper, self.model))

    # for now this is not cached
    def calc(self, p, *args, **kwargs):
        """Evaluate the convolved model on a grid.

        Parameters
        ----------
        p : sequence of numbers
            The parameters of the model, matching the ``pars``
            field. This will start with the convolution model
            parameters (if any) and then the model.
        *args
            The model grid. There should be two arrays (the low and
            high edges of the bin) to make sure the wrapped model is
            evaluated correctly.
        **kwargs
            Additional keyword arguments.

        """

        nargs = 1 + len(args)
        if nargs != 3:
            emsg = f"calc() requires pars,lo,hi arguments, sent {nargs} arguments"
            warnings.warn(emsg, FutureWarning)
            # raise TypeError(emsg)

        return self.wrapper.calc(p, self.model.calc,
                                 *args, **kwargs)


# Models from model.dat - try to follow the ordering of that file
#

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
        self.LineE = XSParameter(name, 'LineE', 10.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.Sigma = XSParameter(name, 'Sigma', 1.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Sigma, self.norm))


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
    XSagnslim, XSqsosed

    Notes
    -----
    This model is only available when used with XSPEC 12.10.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelAgnsed.html

    """

    __function__ = "agnsed"

    def __init__(self, name='agnsed'):
        self.mass = XSParameter(name, 'mass', 1e7, 1.0, 1e10, 1.0, 1e10,
                                'solar', frozen=True)
        self.dist = XSParameter(name, 'dist', 100, 0.01, 1e9, 0.01, 1e9,
                                'Mpc', frozen=True)
        self.logmdot = XSParameter(name, 'logmdot', -1, -10, 2, -10, 2)
        self.astar = XSParameter(name, 'astar', 0.0, -1, 0.998, -1, 0.998,
                                 frozen=True)
        self.cosi = XSParameter(name, 'cosi', 0.5, 0.05, 1.0, 0.05, 1.0,
                                frozen=True)
        # TODO: allow negative values
        self.kTe_hot = XSParameter(name, 'kTe_hot', 100.0, 10, 300, 10, 300,
                                   'keV(_pl)', frozen=True)
        self.kTe_warm = XSParameter(name, 'kTe_warm', 0.2, 0.1, 0.5, 0.1, 0.5,
                                    'keV(_sc)')
        self.Gamma_hot = XSParameter(name, 'Gamma_hot', 1.7, 1.3, 3, 1.3, 3,
                                     '(_calc)')
        self.Gamma_warm = XSParameter(name, 'Gamma_warm', 2.7, 2, 5, 2, 10,
                                      '(_disk)')

        self.R_hot = XSParameter(name, 'R_hot', 10, 6, 500, 6, 500, 'Rg')
        self.R_warm = XSParameter(name, 'R_warm', 20, 6, 500, 6, 500, 'Rg')

        self.logrout = XSParameter(name, 'logrout', -1, -3, 7, -3, 7,
                                   '(_selfg)', frozen=True)

        self.Htmax = XSParameter(name, 'Htmax', 10, 6, 10, 6, 10,
                                 'Rg', frozen=True)

        self.reprocess = XSParameter(name, 'reprocess', 1, 0, 1, 0, 1,
                                     '0off/1on', alwaysfrozen=True)

        self.redshift = XSParameter(name, 'redshift', 0, 0, 1, 0, 1,
                                    frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.mass, self.dist, self.logmdot, self.astar, self.cosi,
                self.kTe_hot, self.kTe_warm, self.Gamma_hot, self.Gamma_warm,
                self.R_hot, self.R_warm, self.logrout, self.Htmax,
                self.reprocess, self.redshift, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


@version_at_least("12.11.0")
class XSagnslim(XSAdditiveModel):
    """The XSPEC agnslim model: AGN super-Eddington accretion model

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The default logmdot parameter value has changed from -1 to 1 to
       match XSPEC.

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
        The spectral index of the hot Comptonisation component.
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
    rin
        The inner radius of the disc in Rg. If this parameter is -1
        (the default), the model will use the radius calculated from KD19.
        This must be greater than R_hot for mdot greater than 6 and greater
        than R_isco for mdot less than 6.
    redshift
        The redshift.
    norm
        The normalization of the model. This must be fixed at 1.

    See Also
    --------
    XSagnsed, XSqsosed

    Notes
    -----
    This model is only available when used with XSPEC 12.11.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelAgnslim.html

    """

    __function__ = "agnslim"

    def __init__(self, name='agnslim'):
        self.mass = XSParameter(name, 'mass', 1e7, 1.0, 1e10, 1.0, 1e10,
                                'solar', frozen=True)
        self.dist = XSParameter(name, 'dist', 100, 0.01, 1e9, 0.01, 1e9,
                                'Mpc', frozen=True)
        self.logmdot = XSParameter(name, 'logmdot', 1, -10, 3, -10, 3)
        self.astar = XSParameter(name, 'astar', 0.0, 0, 0.998, 0, 0.998,
                                 frozen=True)
        self.cosi = XSParameter(name, 'cosi', 0.5, 0.05, 1.0, 0.05, 1.0,
                                frozen=True)
        self.kTe_hot = XSParameter(name, 'kTe_hot', 100.0, 10, 300, 10, 300,
                                   'keV(-pl)', frozen=True)
        self.kTe_warm = XSParameter(name, 'kTe_warm', 0.2, 0.1, 0.5, 0.1, 0.5,
                                    'keV(-sc)')
        self.Gamma_hot = XSParameter(name, 'Gamma_hot', 2.4, 1.3, 3, 1.3, 3)
        self.Gamma_warm = XSParameter(name, 'Gamma_warm', 3.0, 2, 5, 2, 10,
                                      '(-disk)')

        self.R_hot = XSParameter(name, 'R_hot', 10, 2, 500, 2, 500, units='Rg')
        self.R_warm = XSParameter(name, 'R_warm', 20, 2, 500, 2, 500, units='Rg')

        self.logrout = XSParameter(name, 'logrout', -1, -3, 7, -3, 7,
                                   '(-selfg)', frozen=True)

        self.rin = XSParameter(name, 'rin', -1, -1, 100, -1, 100,
                               frozen=True)  # TODO: make alwaysfrozen?

        self.redshift = XSParameter(name, 'redshift', 0, 0, 5, 0, 5,
                                    frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.mass, self.dist, self.logmdot, self.astar, self.cosi,
                self.kTe_hot, self.kTe_warm, self.Gamma_hot, self.Gamma_warm,
                self.R_hot, self.R_warm, self.logrout, self.rin,
                self.redshift, self.norm)

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
    XSbapec, XSbvapec, XSbvvapec, XSnlapec, XSsnapec, XSvapec, XSvvapec, XSwdem

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelApec.html

    """
    __function__ = "C_apec" if equal_or_greater_than("12.9.1") else "xsaped"

    def __init__(self, name='apec'):
        self.kT = XSParameter(name, 'kT', 1., 0.008, 64.0, 0.008, 64.0, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1., 0., 5., 0.0, 5, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.kT = XSParameter(name, 'kT', 1., 0.008, 64.0, 0.008, 64, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1., 0., 5., 0.0, 5, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.Velocity = XSParameter(name, 'Velocity', 0., 0., 1.e6, 0.0, 1e6, units='km/s', frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.Redshift, self.Velocity, self.norm))


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
        self.kT = XSParameter(name, 'kT', 1.0, 0.008, 64.0, 0.008, 64.0, units='keV')
        self.kTi = XSParameter(name, 'kTi', 1.0, 0.008, 64.0, 0.008, 64.0, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1.0, 0.0, 5.0, 0.0, 5.0, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10, frozen=True)
        self.Velocity = XSParameter(name, 'Velocity', 0.0, 0.0, 1.0e6, 0.0, 1.0e6, 'km/s', frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.Abundanc, self.Redshift, self.Velocity, self.norm))


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
        self.kT = XSParameter(name, 'kT', 3.0, 1.e-2, 100., 1e-4, 200, units='keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
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
        self.kT = XSParameter(name, 'kT', 3., 1e-3, 100, 1e-4, 200, units='keV')
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
        The abundance of the elements heavier than He relative to their
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
    XSbexriv, XSreflect

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
        self.Gamma1 = XSParameter(name, 'Gamma1', 2., -9., 9., -10, 10)
        self.breakE = XSParameter(name, 'breakE', 10., 0.1, 1000., 0.1, 1000, units='keV')
        self.Gamma2 = XSParameter(name, 'Gamma2', 2., -9., 9., -10, 10)
        self.foldE = XSParameter(name, 'foldE', 100., 1., 1.e6, 1.0, 1e6, units='keV')
        self.rel_refl = XSParameter(name, 'rel_refl', 0., 0., 10., 0.0, 10)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.45, 0.05, 0.95, 0.05, 0.95, frozen=True)
        self.abund = XSParameter(name, 'abund', 1., 0.0, 1.e6, 0.0, 1e6, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0.0, 1.e6, 0.0, 1e6, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        The abundance of the elements heavier than He relative to their
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
    XSbexrav, XSireflect

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
        self.Gamma1 = XSParameter(name, 'Gamma1', 2., -9., 9., -10, 10)
        self.breakE = XSParameter(name, 'breakE', 10., 0.1, 1000., 0.1, 1000, units='keV')
        self.Gamma2 = XSParameter(name, 'Gamma2', 2., -9., 9., -10, 10)
        self.foldE = XSParameter(name, 'foldE', 100., 1., 1.e6, 1, 1e6, units='keV')
        self.rel_refl = XSParameter(name, 'rel_refl', 0., 0., 1.e6, 0.0, 1e6)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.abund = XSParameter(name, 'abund', 1., 0.0, 1.e6, 0.0, 1e6, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0.0, 1.e6, 0.0, 1e6, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.45, 0.05, 0.95, 0.05, 0.95, frozen=True)
        self.T_disk = XSParameter(name, 'T_disk', 3.e4, 1.e4, 1.e6, 1e4, 1e6, units='K', frozen=True)
        self.xi = XSParameter(name, 'xi', 1., 0., 1.e3, 0.0, 5e3, units='erg cm/s')
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
        self.PhoIndx1 = XSParameter(name, 'PhoIndx1', 1., -2., 9., -3, 10)
        self.BreakE = XSParameter(name, 'BreakE', 5., 1.e-2, 1.e6, 0.0, 1e6, units='keV')
        self.PhoIndx2 = XSParameter(name, 'PhoIndx2', 2., -2., 9., -3, 10)
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
        self.PhoIndx1 = XSParameter(name, 'PhoIndx1', 1., -2., 9., -3, 10)
        self.BreakE1 = XSParameter(name, 'BreakE1', 5., 1.e-2, 1.e6, 0.0, 1e6, units='keV')
        self.PhoIndx2 = XSParameter(name, 'PhoIndx2', 2., -2., 9., -3, 10)
        self.BreakE2 = XSParameter(name, 'BreakE2', 10., 1.e-2, 1.e6, 0.0, 1e6, units='keV')
        self.PhoIndx3 = XSParameter(name, 'PhoIndx3', 3., -2., 9., -3, 10)
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
        self.kT = XSParameter(name, 'kT', 1., 1.e-2, 100., 1e-4, 200, units='keV')
        self.alpha = XSParameter(name, 'alpha', 1., 1.e-2, 4.0, 1e-4, 6)
        self.log_A = XSParameter(name, 'log_A', 0.0, -6.0, 6.0, -8, 8, aliases=["logA"])
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
        self.kT = XSParameter(name, 'kT', 7.0, 1.e-4, 100., 1e-4, 200, units='keV')
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
        self.kT = XSParameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = XSParameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.Redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.Velocity = XSParameter(name, 'Velocity', 0.0, 0.0, 1.0e6, 0.0, 1.0e6,
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
        self.kT = XSParameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.Velocity = XSParameter(name, 'Velocity', 0., 0., 1.e6, 0.0, 1e6, units='km/s', frozen=True)
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
        self.kT = XSParameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = XSParameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1., 0.0, 1.0, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 10000., frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.Redshift = XSParameter(name, 'Redshift', 0., -0.999, 10., -0.999, 10., frozen=True)
        self.Velocity = XSParameter(name, 'Velocity', 0., 0., 1.e6, 0.0, 1.e6, 'km/s', frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.Redshift, self.Velocity, self.norm))


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
        self.kT = XSParameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.kTi = XSParameter(name, 'kTi', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10)
        self.Velocity = XSParameter(name, 'Velocity', 0.0, 0.0, 1.0e6, 0.0, 1.0e6, units='km/s')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.He, self.C, self.N, self.O,
                                  self.Ne, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca,
                                  self.Fe, self.Ni, self.Redshift, self.Velocity, self.norm))


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
        self.kT = XSParameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.Velocity = XSParameter(name, 'Velocity', 0., 0., 1.e6, 0.0, 1e6, units='km/s', frozen=True)
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
        self.kT = XSParameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = XSParameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        maxval = 1000.0
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, maxval, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.Redshift = XSParameter(name, 'Redshift', 0., -0.999, 10., -0.999, 10., frozen=True)
        self.Velocity = XSParameter(name, 'Velocity', 0., 0., 1.e6, 0.0, 1.e6,
                                    units='km/s', frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.Redshift, self.Velocity, self.norm))


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
        self.kT = XSParameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.kTi = XSParameter(name, 'kTi', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10)
        self.Velocity = XSParameter(name, 'Velocity', 0.0, 0.0, 1.0e6, 0.0, 1.0e6, units='km/s')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Redshift, self.Velocity, self.norm))


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
        self.CPcoef1 = XSParameter(name, 'CPcoef1', 1.0, -1, 1, -1, 1)
        self.CPcoef2 = XSParameter(name, 'CPcoef2', 0.5, -1, 1, -1, 1)
        self.CPcoef3 = XSParameter(name, 'CPcoef3', 0.5, -1, 1, -1, 1)
        self.CPcoef4 = XSParameter(name, 'CPcoef4', 0.5, -1, 1, -1, 1)
        self.CPcoef5 = XSParameter(name, 'CPcoef5', 0.5, -1, 1, -1, 1)
        self.CPcoef6 = XSParameter(name, 'CPcoef6', 0.5, -1, 1, -1, 1)
        self.nH = XSParameter(name, 'nH', 1.0, 1.e-5, 1.e19, 1e-6, 1e20, units='cm^-3', frozen=True)
        self.abundanc = XSParameter(name, 'abundanc', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 1, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5,
                                              self.CPcoef6, self.nH, self.abundanc, self.redshift, self.switch,
                                              self.norm))


class XSc6pmekl(XSAdditiveModel):
    """The XSPEC c6pmekl model: differential emission measure using Chebyshev representations with multi-temperature mekal.

    The model is described at [1]_. It differs from ``XSc6mekl`` by
    using the exponential of the 6th order Chebyshev polynomial.

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
        self.CPcoef1 = XSParameter(name, 'CPcoef1', 1.0, -1, 1, -1, 1)
        self.CPcoef2 = XSParameter(name, 'CPcoef2', 0.5, -1, 1, -1, 1)
        self.CPcoef3 = XSParameter(name, 'CPcoef3', 0.5, -1, 1, -1, 1)
        self.CPcoef4 = XSParameter(name, 'CPcoef4', 0.5, -1, 1, -1, 1)
        self.CPcoef5 = XSParameter(name, 'CPcoef5', 0.5, -1, 1, -1, 1)
        self.CPcoef6 = XSParameter(name, 'CPcoef6', 0.5, -1, 1, -1, 1)
        self.nH = XSParameter(name, 'nH', 1.0, 1.e-5, 1.e19, 1e-6, 1e20, units='cm^-3', frozen=True)
        self.abundanc = XSParameter(name, 'abundanc', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 1, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5,
                                              self.CPcoef6, self.nH, self.abundanc, self.redshift, self.switch,
                                              self.norm))


class XSc6pvmkl(XSAdditiveModel):
    """The XSPEC c6pvmkl model: differential emission measure using Chebyshev representations with multi-temperature mekal.

    The model is described at [1]_. It differs from ``XSc6vmekl`` by
    using the exponential of the 6th order Chebyshev polynomial.

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
        self.CPcoef1 = XSParameter(name, 'CPcoef1', 1.0, -1, 1, -1, 1)
        self.CPcoef2 = XSParameter(name, 'CPcoef2', 0.5, -1, 1, -1, 1)
        self.CPcoef3 = XSParameter(name, 'CPcoef3', 0.5, -1, 1, -1, 1)
        self.CPcoef4 = XSParameter(name, 'CPcoef4', 0.5, -1, 1, -1, 1)
        self.CPcoef5 = XSParameter(name, 'CPcoef5', 0.5, -1, 1, -1, 1)
        self.CPcoef6 = XSParameter(name, 'CPcoef6', 0.5, -1, 1, -1, 1)
        self.nH = XSParameter(name, 'nH', 1.0, 1.e-5, 1.e19, 1e-6, 1e20, units='cm^-3', frozen=True)
        self.He = XSParameter(name, 'He', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Na = XSParameter(name, 'Na', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Al = XSParameter(name, 'Al', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 1, alwaysfrozen=True)
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
        self.CPcoef1 = XSParameter(name, 'CPcoef1', 1.0, -1, 1, -1, 1)
        self.CPcoef2 = XSParameter(name, 'CPcoef2', 0.5, -1, 1, -1, 1)
        self.CPcoef3 = XSParameter(name, 'CPcoef3', 0.5, -1, 1, -1, 1)
        self.CPcoef4 = XSParameter(name, 'CPcoef4', 0.5, -1, 1, -1, 1)
        self.CPcoef5 = XSParameter(name, 'CPcoef5', 0.5, -1, 1, -1, 1)
        self.CPcoef6 = XSParameter(name, 'CPcoef6', 0.5, -1, 1, -1, 1)
        self.nH = XSParameter(name, 'nH', 1.0, 1.e-5, 1.e19, 1e-6, 1e20, units='cm^-3', frozen=True)
        self.He = XSParameter(name, 'He', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Na = XSParameter(name, 'Na', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Al = XSParameter(name, 'Al', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 1, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5,
                                              self.CPcoef6, self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na,
                                              self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni,
                                              self.redshift, self.switch, self.norm))


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
        self.T = XSParameter(name, 'T', 2.0, 1.0, 4.0, 1.0, 4.0, units='MK')
        self.NSmass = XSParameter(name, 'NSmass', 1.4, 0.6, 2.8, 0.6, 2.8,
                                  units='Msun')
        self.NSrad = XSParameter(name, 'NSrad', 10.0, 6.0, 23.0, 6.0, 23.0, units='km')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.T, self.NSmass, self.NSrad, self.norm))


class XScemekl(XSAdditiveModel):
    """The XSPEC cemekl model: plasma emission, multi-temperature using mekal.

    The model is described at [1]_.

    Attributes
    ----------
    alpha
        The power-law index of the emissivity function.
    Tmax
        The maximum temperature, in keV.
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
        self.alpha = XSParameter(name, 'alpha', 1.0, 0.01, 10, 0.01, 20, frozen=True)
        self.Tmax = XSParameter(name, 'Tmax', 1.0, 0.01, 1.e2, 0.01, 1e2, units='keV')
        self.nH = XSParameter(name, 'nH', 1.0, 1.e-5, 1.e19, 1e-6, 1e20, units='cm^-3', frozen=True)
        self.abundanc = XSParameter(name, 'abundanc', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 1, 0, 1, 0, 1, alwaysfrozen=True)
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
        The maximum temperature, in keV.
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
        self.alpha = XSParameter(name, 'alpha', 1.0, 0.01, 10, 0.01, 20, frozen=True)
        self.Tmax = XSParameter(name, 'Tmax', 1.0, 0.01, 1.e2, 0.01, 1e2, units='keV')
        self.nH = XSParameter(name, 'nH', 1.0, 1.e-5, 1.e19, 1e-6, 1e20, units='cm^-3', frozen=True)
        self.He = XSParameter(name, 'He', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Na = XSParameter(name, 'Na', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Al = XSParameter(name, 'Al', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 1, 0, 1, 0, 1, alwaysfrozen=True)
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
        The maximum temperature, in keV.
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
        self.slope = XSParameter(name, 'slope', 0., -5., 5., -5, 5)
        self.lowT = XSParameter(name, 'lowT', 0.1, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.highT = XSParameter(name, 'highT', 4., 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1., 0.0, 5., 0.0, 5)
        self.redshift = XSParameter(name, 'redshift', .1, 1.e-10, 10., 1e-10, 10, frozen=True)
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
        self.kT = XSParameter(name, 'kT', 1.0, 1.e-2, 100., 1e-4, 200, units='keV')
        self.kTe = XSParameter(name, 'kTe', 50, 1., 200., 1.0, 200, units='keV', frozen=True)
        self.tau = XSParameter(name, 'tau', 0.1, 0.0, 10., 0.0, 10)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kTe, self.tau, self.norm))


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
        self.kTbb = XSParameter(name, 'kTbb', 1.0, 0.2, 10.0, 0.2, 10.0, units='keV')
        self.kTe = XSParameter(name, 'kTe', 5.0, 0.2, 2000.0, 0.2, 2000.0, units='keV')
        self.tau = XSParameter(name, 'tau', 0.5, 0.0, 10.0, 0.0, 10.0)
        self.eta = XSParameter(name, 'eta', 0.5, 0.01, 1.0, 0.01, 1.0)
        self.beta0 = XSParameter(name, 'beta0', 0.57, 1.0e-4, 1.0, 1.0e-4, 1.0)
        self.r0 = XSParameter(name, 'r0', 0.25, 1.0e-4, 100.0, 1.0e-4, 100.0)
        self.A = XSParameter(name, 'A', 0.001, 0, 1, 0, 1, frozen=True)
        self.betaflag = XSParameter(name, 'betaflag', 1, 0, 2, 0, 2, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kTbb, self.kTe, self.tau, self.eta, self.beta0, self.r0, self.A, self.betaflag, self.norm))


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
        self.kT = XSParameter(name, 'kT', 2., .01, 10., 1e-3, 20, units='keV')
        self.tau = XSParameter(name, 'tau', 10, .001, 100., 1e-4, 200)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.tau, self.norm))


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
        self.kTe = XSParameter(name, 'kTe', 100., 20., 1.e5, 20.0, 1e5, units='keV')
        self.EleIndex = XSParameter(name, 'EleIndex', 2., 0.0, 5., 0.0, 5, frozen=True)
        self.Gmin = XSParameter(name, 'Gmin', -1., -1., 10., -1, 10, frozen=True)
        self.Gmax = XSParameter(name, 'Gmax', 1.e3, 10., 1.e4, 10, 1e4, frozen=True)
        self.kTbb = XSParameter(name, 'kTbb', 0.1, 0.001, 10., 0.001, 10, units='keV', frozen=True)
        self.tau_y = XSParameter(name, 'tau_y', 1.0, 0.05, 3.0, 0.05, 3, aliases=["tauy"])
        self.geom = XSParameter(name, 'geom', 0.0, -5.0, 4.0, -5, 4, frozen=True)
        self.HovR_cyl = XSParameter(name, 'HovR_cyl', 1.0, 0.5, 2.0, 0.5, 2, frozen=True, aliases=["HRcyl"])
        self.cosIncl = XSParameter(name, 'cosIncl', 0.5, 0.05, 0.95, 0.05, 0.95, frozen=True)
        self.cov_frac = XSParameter(name, 'cov_frac', 1.0, 0.0, 1.0, 0.0, 1, frozen=True)
        self.rel_refl = XSParameter(name, 'rel_refl', 0., 0., 1.e4, 0.0, 1e4, frozen=True)
        self.Fe_ab_re = XSParameter(name, 'Fe_ab_re', 1., 0.1, 10., 0.1, 10, frozen=True)
        self.Me_ab = XSParameter(name, 'Me_ab', 1., 0.1, 10., 0.1, 10, frozen=True)
        self.xi = XSParameter(name, 'xi', 0., 0., 1.e5, 0.0, 1e5, frozen=True)
        self.Tdisk = XSParameter(name, 'Tdisk', 1.e6, 1.e4, 1.e6, 1e4, 1e6, units='K', frozen=True)
        self.Betor10 = XSParameter(name, 'Betor10', -10., -10., 10., -10, 10, frozen=True)
        self.Rin = XSParameter(name, 'Rin', 10., 6.001, 1.e3, 6.001, 1e4, units='Rs', frozen=True)
        self.Rout = XSParameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, 1e6, units='Rs', frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.kT = XSParameter(name, 'kT', 2., .01, 100., 1e-3, 100, units='keV')
        self.tau = XSParameter(name, 'tau', 10, .001, 100., 1e-4, 200)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.tau, self.norm))


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
        self.kTs = XSParameter(name, 'kTs', 1.0, 0.1, 10.0, 0.1, 10.0, units='keV')
        self.gamma = XSParameter(name, 'gamma', 3.0, 1.0, 10.0, 1.0, 10.0, frozen=True)
        self.alpha = XSParameter(name, 'alpha', 2.0, 0.0, 400.0, 0.0, 400.0)
        self.delta = XSParameter(name, 'delta', 20.0, 0.0, 200.0, 0.0, 200.0)
        self.kTe = XSParameter(name, 'kTe', 5.0, 0.2, 2000.0, 0.2, 2000.0, units='keV')
        self.log_A = XSParameter(name, 'log_A', 0.0, -8.0, 8.0, -8.0, 8.0)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kTs, self.gamma, self.alpha, self.delta, self.kTe, self.log_A, self.norm))


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
        self.theta = XSParameter(name, 'theta', 1., 1e-6, 1.e6, 1e-6, 1e6, units='keV')
        self.showbb = XSParameter(name, 'showbb', 1.0, 0., 1.e4, 0.0, 1e4, frozen=True)
        self.kT_bb = XSParameter(name, 'kT_bb', 200., 1., 4e5, 1, 4e5, units='eV', frozen=True)
        self.RefOn = XSParameter(name, 'RefOn', -1.0, -2.0, 2.0, -2, 2, frozen=True)
        self.tau_p = XSParameter(name, 'tau_p', 0.1, 1e-4, 10., 1e-4, 10, frozen=True)
        self.radius = XSParameter(name, 'radius', 1.e7, 1.e5, 1.e16, 1e5, 1e16, units='cm', frozen=True)
        self.g_min = XSParameter(name, 'g_min', 1.3, 1.2, 1.e3, 1.2, 1e3, frozen=True)
        self.g_max = XSParameter(name, 'g_max', 1.e3, 5., 1.e4, 5, 1e4, frozen=True)
        self.G_inj = XSParameter(name, 'G_inj', 2., 0., 5., 0.0, 5, frozen=True)
        self.pairinj = XSParameter(name, 'pairinj', 0., 0., 1., 0.0, 1, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.50, 0.05, 0.95, 0.05, 0.95, frozen=True)
        self.Refl = XSParameter(name, 'Refl', 1., 0., 2., 0.0, 2, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0.1, 10., 0.1, 10, frozen=True)
        self.Ab_met = XSParameter(name, 'Ab_met', 1.0, 0.1, 10., 0.1, 10, frozen=True, aliases=["AbHe"])
        self.T_disk = XSParameter(name, 'T_disk', 1.e6, 1e4, 1e6, 1e4, 1e6, units='K', frozen=True)
        self.xi = XSParameter(name, 'xi', 0.0, 0.0, 1000.0, 0.0, 5000)
        self.Beta = XSParameter(name, 'Beta', -10., -10., 10., -10, 10, frozen=True)
        self.Rin = XSParameter(name, 'Rin', 10., 6.001, 1.e3, 6.001, 1e4, units='M', frozen=True)
        self.Rout = XSParameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, 1e6, units='M', frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., 0., 4., 0.0, 4, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.theta, self.showbb, self.kT_bb, self.RefOn, self.tau_p, self.radius, self.g_min, self.g_max, self.G_inj, self.pairinj, self.cosIncl, self.Refl, self.Fe_abund, self.Ab_met, self.T_disk, self.xi, self.Beta, self.Rin, self.Rout, self.redshift, self.norm))


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
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.T0 = XSParameter(name, 'T0', 0.1, .01, 100., 1e-3, 100, units='keV')
        self.kT = XSParameter(name, 'kT', 50., 2.0, 500., 2.0, 500, units='keV')
        self.taup = XSParameter(name, 'taup', 1., .01, 100., 0.01, 200)
        self.approx = XSParameter(name, 'approx', 1.0, 0.0, 5.0, 0.0, 200, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.redshift, self.T0, self.kT, self.taup, self.approx, self.norm))


@version_at_least("12.10.1")
class XScph(XSAdditiveModel):
    """The XSPEC cph model: Cooling + heating model for cool core clusters

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The default Redshift parameter value has changed from 0.1 to 0
       to match XSPEC.

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
        self.peakT = XSParameter(name, 'peakT', 2.2, 1e-1, 1e2, 1e-1, 1e2,
                                 units='keV')
        self.Abund = XSParameter(name, 'Abund', 1, 0, 1000, 0, 1000,
                                 frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, 0.0, 50, 0.0, 50,
                                    frozen=True)
        self.switch = XSParameter(name, 'switch', 1,
                                  alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.peakT, self.Abund, self.Redshift, self.switch, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


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
        self.energy00 = XSParameter(name, 'energy00', 0.5, 0.0, 10.0, 0, 10, units='keV', alwaysfrozen=True)
        self.energy01 = XSParameter(name, 'energy01', 1., 0.0, 10.0, 0, 10, units='keV', alwaysfrozen=True)
        self.energy02 = XSParameter(name, 'energy02', 1.5, 0.0, 10.0, 0, 10, units='keV', alwaysfrozen=True)
        self.energy03 = XSParameter(name, 'energy03', 2., 0.0, 10.0, 0, 10, units='keV', alwaysfrozen=True)
        self.energy04 = XSParameter(name, 'energy04', 3., 0.0, 10.0, 0, 10, units='keV', alwaysfrozen=True)
        self.energy05 = XSParameter(name, 'energy05', 4., 0.0, 10.0, 0, 10, units='keV', alwaysfrozen=True)
        self.energy06 = XSParameter(name, 'energy06', 5., 0.0, 10.0, 0, 10, units='keV', alwaysfrozen=True)
        self.energy07 = XSParameter(name, 'energy07', 6., 0.0, 10.0, 0, 10, units='keV', alwaysfrozen=True)
        self.energy08 = XSParameter(name, 'energy08', 7., 0.0, 10.0, 0, 10, units='keV', alwaysfrozen=True)
        self.energy09 = XSParameter(name, 'energy09', 8., 0.0, 10.0, 0, 10, units='keV', alwaysfrozen=True)
        self.log_rate00 = XSParameter(name, 'log_rate00', 0., -19.0, 19.0, -20, 20, frozen=True)
        self.log_rate01 = XSParameter(name, 'log_rate01', 1., -19.0, 19.0, -20, 20, frozen=True)
        self.log_rate02 = XSParameter(name, 'log_rate02', 0., -19.0, 19.0, -20, 20, frozen=True)
        self.log_rate03 = XSParameter(name, 'log_rate03', 1., -19.0, 19.0, -20, 20, frozen=True)
        self.log_rate04 = XSParameter(name, 'log_rate04', 0., -19.0, 19.0, -20, 20, frozen=True)
        self.log_rate05 = XSParameter(name, 'log_rate05', 1., -19.0, 19.0, -20, 20, frozen=True)
        self.log_rate06 = XSParameter(name, 'log_rate06', 0., -19.0, 19.0, -20, 20, frozen=True)
        self.log_rate07 = XSParameter(name, 'log_rate07', 1., -19.0, 19.0, -20, 20, frozen=True)
        self.log_rate08 = XSParameter(name, 'log_rate08', 0., -19.0, 19.0, -20, 20, frozen=True)
        self.log_rate09 = XSParameter(name, 'log_rate09', 1., -19.0, 19.0, -20, 20, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.energy00, self.energy01, self.energy02, self.energy03,
                                              self.energy04, self.energy05, self.energy06, self.energy07,
                                              self.energy08, self.energy09, self.log_rate00, self.log_rate01,
                                              self.log_rate02, self.log_rate03, self.log_rate04, self.log_rate05,
                                              self.log_rate06, self.log_rate07, self.log_rate08, self.log_rate09,
                                              self.norm))


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
        self.PhoIndex = XSParameter(name, 'PhoIndex', 1., -2., 9., -3, 10)
        self.HighECut = XSParameter(name, 'HighECut', 15., 1., 500., 0.01, 500, units='keV')
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
        self.accrate = XSParameter(name, 'accrate', 1., 1e-3, 9., 1e-4, 10)
        self.CenMass = XSParameter(name, 'CenMass', 1.4, .4, 10., 0.1, 20, units='Msun', frozen=True, aliases=["NSmass"])
        self.Rinn = XSParameter(name, 'Rinn', 1.03, 1.01, 1.03, 1.0, 1.04, frozen=True)
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
        self.kT_disk = XSParameter(name, 'kT_disk', 1.0, 0.01, 5., 0.01, 5, units='keV')
        self.Gamma = XSParameter(name, 'Gamma', 1.7, 1.001, 5., 1.001, 10)
        self.kT_e = XSParameter(name, 'kT_e', 100., 5., 1.e3, 1.0, 1000, units='keV')
        self.LcovrLd = XSParameter(name, 'LcovrLd', 0.1, 0., 10., 0.0, 10, aliases=["LcLd"])
        self.fin = XSParameter(name, 'fin', 1.e-1, 0.0, 1., 0.0, 1, frozen=True)
        self.rirr = XSParameter(name, 'rirr', 1.2, 1.0001, 10., 1.0001, 10)
        self.fout = XSParameter(name, 'fout', 1.e-4, 0.0, 1.e-1, 0.0, 0.1)
        self.logrout = XSParameter(name, 'logrout', 5.0, 3.0, 7.0, 3, 7)
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
        self.Tin = XSParameter(name, 'Tin', 1., 0., 1000., 0.0, 1000, units='keV')
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

    See Also
    --------
    XSrdblur

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDiskline.html

    """

    __function__ = "C_diskline" if equal_or_greater_than("12.10.1") else "xsdili"

    def __init__(self, name='diskline'):
        self.LineE = XSParameter(name, 'LineE', 6.7, 0., 100., 0.0, 100, units='keV')
        self.Betor10 = XSParameter(name, 'Betor10', -2., -10., 20., -10, 20, frozen=True)
        self.Rin_M = XSParameter(name, 'Rin_M', 10., 6., 1000., 6.0, 10000, frozen=True, aliases=["RinM"])
        self.Rout_M = XSParameter(name, 'Rout_M', 1000., 0., 1000000., 0.0, 10000000, frozen=True, aliases=["RoutM"])
        self.Incl = XSParameter(name, 'Incl', 30., 0., 90., 0.0, 90, units='deg')
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
    accrate
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
        self.accrate = XSParameter(name, 'accrate', 1., 1e-3, 9., 1e-4, 10)
        self.NSmass = XSParameter(name, 'NSmass', 1.4, .4, 10., 0.1, 20, units='Msun', frozen=True)
        self.Rinn = XSParameter(name, 'Rinn', 1.03, 1.01, 1.03, 1.0, 1.04, frozen=True)
        self.alpha = XSParameter(name, 'alpha', 1., .01, 10., 0.001, 20, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.accrate, self.NSmass, self.Rinn, self.alpha, self.norm))


class XSdisko(XSAdditiveModel):
    """The XSPEC disko model: accretion disk, inner, radiation pressure viscosity.

    The model is described at [1]_.

    Attributes
    ----------
    accrate
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
        self.accrate = XSParameter(name, 'accrate', 1., 1e-3, 9., 1e-4, 10)
        self.NSmass = XSParameter(name, 'NSmass', 1.4, .4, 10., 0.1, 20, units='Msun', frozen=True)
        self.Rinn = XSParameter(name, 'Rinn', 1.03, 1.01, 1.03, 1.0, 1.04, frozen=True)
        self.alpha = XSParameter(name, 'alpha', 1., .01, 10., 0.001, 20, frozen=True)
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
        self.Tin = XSParameter(name, 'Tin', 1.0, 0.1, 10.0, 0.1, 10, units='keV')
        self.p = XSParameter(name, 'p', 0.75, 0.5, 1.0, 0.5, 1)
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
        self.T_max = XSParameter(name, 'T_max', 1., 1e-3, 100, 1e-4, 200, units='keV')
        self.R_in = XSParameter(name, 'R_in', 6., 6., 1000., 6.0, 1000, units='R_g')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.T_max, self.R_in, self.norm))


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
        self.Ep = XSParameter(name, 'Ep', .1, 1.e-6, 1.e2, 1e-10, 1e4, units='keV')
        self.beta = XSParameter(name, 'beta', 0.2, -4., 4., -4, 4)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Ep, self.beta, self.norm))


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
        self.l_hovl_s = XSParameter(name, 'l_hovl_s', 1., 1e-6, 1.e6, 1e-6, 1e6, aliases=["l_hl_s"])
        self.l_bb = XSParameter(name, 'l_bb', 100., 0., 1.e4, 0.0, 1e4)
        self.kT_bb = XSParameter(name, 'kT_bb', 200., 1., 4e5, 1, 4e5, units='eV', frozen=True)
        self.l_ntol_h = XSParameter(name, 'l_ntol_h', 0.5, 0., 0.9999, 0.0, 0.9999, aliases=["l_ntl_h"])
        self.tau_p = XSParameter(name, 'tau_p', 0.1, 1e-4, 10., 1e-4, 10, frozen=True)
        self.radius = XSParameter(name, 'radius', 1.e7, 1.e5, 1.e16, 1e5, 1e16, units='cm', frozen=True)
        self.g_min = XSParameter(name, 'g_min', 1.3, 1.2, 1.e3, 1.2, 1e3, frozen=True)
        self.g_max = XSParameter(name, 'g_max', 1.e3, 5., 1.e4, 5, 1e4, frozen=True)
        self.G_inj = XSParameter(name, 'G_inj', 2., 0., 5., 0.0, 5, frozen=True)
        self.pairinj = XSParameter(name, 'pairinj', 0., 0., 1., 0.0, 1, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.50, 0.05, 0.95, 0.05, 0.95, frozen=True)
        self.Refl = XSParameter(name, 'Refl', 1., 0., 2., 0.0, 2, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0.1, 10., 0.1, 10, frozen=True)
        self.Ab_met = XSParameter(name, 'Ab_met', 1.0, 0.1, 10., 0.1, 10, frozen=True, aliases=["AbHe"])
        self.T_disk = XSParameter(name, 'T_disk', 1.e6, 1e4, 1e6, 1e4, 1e6, units='K', frozen=True)
        self.xi = XSParameter(name, 'xi', 0.0, 0.0, 1000.0, 0.0, 5000)
        self.Beta = XSParameter(name, 'Beta', -10., -10., 10., -10, 10, frozen=True)
        self.Rin = XSParameter(name, 'Rin', 10., 6.001, 1.e3, 6.001, 1e4, units='M', frozen=True)
        self.Rout = XSParameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, 1e6, units='M', frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., 0., 4., 0.0, 4, frozen=True)
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
        self.l_hovl_s = XSParameter(name, 'l_hovl_s', 1., 1e-6, 1.e6, 1e-6, 1e6, aliases=["l_hl_s"])
        self.l_bb = XSParameter(name, 'l_bb', 100., 0., 1.e4, 0.0, 1e4)
        self.kT_bb = XSParameter(name, 'kT_bb', 200., 1., 4e5, 1, 4e5, units='eV', frozen=True)
        self.l_ntol_h = XSParameter(name, 'l_ntol_h', 0.5, 0., 0.9999, 0.0, 0.9999, aliases=["l_ntl_h"])
        self.tau_p = XSParameter(name, 'tau_p', 0.1, 1e-4, 10., 1e-4, 10, frozen=True)
        self.radius = XSParameter(name, 'radius', 1.e7, 1.e5, 1.e16, 1e5, 1e16, units='cm', frozen=True)
        self.g_min = XSParameter(name, 'g_min', 1.3, 1.2, 1.e3, 1.2, 1e3, frozen=True)
        self.g_max = XSParameter(name, 'g_max', 1.e3, 5., 1.e4, 5, 1e4, frozen=True)
        self.G_inj = XSParameter(name, 'G_inj', 2., 0., 5., 0.0, 5, frozen=True)
        self.pairinj = XSParameter(name, 'pairinj', 0., 0., 1., 0.0, 1, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.50, 0.05, 0.95, 0.05, 0.95, frozen=True)
        self.Refl = XSParameter(name, 'Refl', 1., 0., 2., 0.0, 2, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0.1, 10., 0.1, 10, frozen=True)
        self.Ab_met = XSParameter(name, 'Ab_met', 1.0, 0.1, 10., 0.1, 10, frozen=True, aliases=["AbHe"])
        self.T_disk = XSParameter(name, 'T_disk', 1.e6, 1e4, 1e6, 1e4, 1e6, units='K', frozen=True)
        self.xi = XSParameter(name, 'xi', 0.0, 0.0, 1000.0, 0.0, 5000)
        self.Beta = XSParameter(name, 'Beta', -10., -10., 10., -10, 10, frozen=True)
        self.Rin = XSParameter(name, 'Rin', 10., 6.001, 1.e3, 6.001, 1e4, units='M', frozen=True)
        self.Rout = XSParameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, 1e6, units='M', frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., 0., 4., 0.0, 4, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.l_hovl_s, self.l_bb, self.kT_bb, self.l_ntol_h, self.tau_p, self.radius, self.g_min, self.g_max, self.G_inj, self.pairinj, self.cosIncl, self.Refl, self.Fe_abund, self.Ab_met, self.T_disk, self.xi, self.Beta, self.Rin, self.Rout, self.redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.factor = XSParameter(name, 'factor', 1.0, 0., 100.0, 0.0, 100)
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
        self.T_max = XSParameter(name, 'T_max', 1., 0.01, 100., 0.01, 100, units='keV')
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
        self.LineE = XSParameter(name, 'LineE', 6.5, 0., 1.e6, 0.0, 1e6, units='keV')
        self.Sigma = XSParameter(name, 'Sigma', 0.1, 0., 10., 0.0, 20, units='keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Sigma, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


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
        self.Tmean = XSParameter(name, 'Tmean', 4.0, 0.01, 10, 0.01, 20, units='keV', frozen=True)
        self.Tsigma = XSParameter(name, 'Tsigma', 0.1, 0.01, 1.e2, 0.01, 1e2, units='keV')
        self.nH = XSParameter(name, 'nH', 1.0, 1.e-5, 1.e19, 1e-6, 1e20, units='cm^-3', frozen=True)
        self.abundanc = XSParameter(name, 'abundanc', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.Tmean, self.Tsigma, self.nH, self.abundanc, self.Redshift,
                                              self.switch, self.norm))


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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 1e8, 5e13, units='s/cm^3')
        self.meankT = XSParameter(name, 'meankT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV', aliases=["kT_ave"])
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
    XSkerrbb

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGrad.html

    """

    __function__ = "grad"

    def __init__(self, name='grad'):
        self.D = XSParameter(name, 'D', 10.0, 0.0, 10000., 0.0, 10000, units='kpc', frozen=True)
        self.i = XSParameter(name, 'i', 0.0, 0.0, 90.0, 0.0, 90, units='deg', frozen=True)
        self.Mass = XSParameter(name, 'Mass', 1.0, 0.0, 100.0, 0.0, 100, units='solar')
        self.Mdot = XSParameter(name, 'Mdot', 1.0, 0.0, 100.0, 0.0, 100, units='1e18')
        self.TclovTef = XSParameter(name, 'TclovTef', 1.7, 1.0, 10.0, 1.0, 10, frozen=True, aliases=["TclTef"])
        self.refflag = XSParameter(name, 'refflag', 1.0, -1.0, 1.0, -1, 1, frozen=True)
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
    XSgrbjet, XSgrbm

    Notes
    -----
    This model is only available when used with XSPEC 12.10.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGrbcomp.html

    """

    __function__ = "xsgrbcomp"

    def __init__(self, name='grbcomp'):
        self.kTs = XSParameter(name, 'kTs', 1.0, 0.0, 20., 0.0, 20.0, units='keV')
        self.gamma = XSParameter(name, 'gamma', 3.0, 0.0, 10.0, 0.0, 10.0)
        self.kTe = XSParameter(name, 'kTe', 100.0, 0.2, 2000., 0.2, 2000., units='keV')
        self.tau = XSParameter(name, 'tau', 5.0, 0.0, 200., 0.0, 200.)
        self.beta = XSParameter(name, 'beta', 0.2, 0.0, 1.0, 0.0, 1.0)
        self.fbflag = XSParameter(name, 'fbflag', 0.0, 0.0, 1.0, 0.0, 1.0, frozen=True)
        self.log_A = XSParameter(name, 'log_A', 5.0, -8., 8., -8., 8., frozen=True)
        self.z = XSParameter(name, 'z', 0.0, 0.0, 10., 0.0, 10., frozen=True)
        self.a_boost = XSParameter(name, 'a_boost', 5.0, 0., 30., 0., 30., frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kTs, self.gamma, self.kTe, self.tau, self.beta, self.fbflag, self.log_A, self.z, self.a_boost, self.norm))


@version_at_least("12.12.0")
class XSgrbjet(XSAdditiveModel):
    """The XSPEC grbjet model: Two-phase Comptonization model of soft thermal seed photons for GRB prompt emission

    The model is described at [1]_.

    Attributes
    ----------
    thobs
        The observing viewing angle in degrees.
    thjet
        The jet half-opening angle in degrees.
    gamma
        The jet gamma Lorentz factor.
    r12
        The jet radius in 10^12 cm.
    p1
        The low-energy index of the coming frame broken powerlaw
        spectrum.
    p2
        The high-energy index of the coming frame broken powerlaw
        spectrum.
    E0
        The break energy in keV.
    delta
        The smoothness of the transition between the two powerlaws.
    index_pl
        The energy index of the comoving-frame cutoff powerlaw
        spectrum.
    ecut
        The cut-off energy in keV.
    ktbb
        The comoving frame blackbody temperature in keV.
    model
        The comoving frame emissivity law: 1 is broken powerlaw,
        2 is cutoff powerlaw, and 3 is blackbody.
    redshift
        The source redshift.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSgrbcomp, XSgrbm

    Notes
    -----
    This model is only available when used with XSPEC 12.12.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGrbjet.html

    """

    __function__ = "xsgrbjet"

    def __init__(self, name='grbjet'):
        self.thobs = XSParameter(name, 'thobs', 5.0, min=0.0, max=30.0, hard_min=0.0, hard_max=30.0, frozen=True)
        self.thjet = XSParameter(name, 'thjet', 10.0, min=2.0, max=20.0, hard_min=2.0, hard_max=20.0, frozen=True)
        self.gamma = XSParameter(name, 'gamma', 200.0, min=1.0, max=500.0, hard_min=1.0, hard_max=500.0)
        self.r12 = XSParameter(name, 'r12', 1.0, min=0.1, max=100.0, hard_min=0.1, hard_max=100.0, frozen=True)
        self.p1 = XSParameter(name, 'p1', 0.0, min=-2.0, max=1.0, hard_min=-2.0, hard_max=1.0)
        self.p2 = XSParameter(name, 'p2', 1.5, min=1.1, max=10.0, hard_min=1.1, hard_max=10.0)
        self.E0 = XSParameter(name, 'E0', 1.0, min=0.1, max=1000.0, hard_min=0.1, hard_max=1000.0, units='keV')
        self.delta = XSParameter(name, 'delta', 0.2, min=0.01, max=1.5, hard_min=0.01, hard_max=1.5, frozen=True)
        self.index_pl = XSParameter(name, 'index_pl', 0.8, min=0.0, max=1.5, hard_min=0.0, hard_max=1.5, frozen=True)
        self.ecut = XSParameter(name, 'ecut', 20.0, min=0.1, max=1000.0, hard_min=0.1, hard_max=1000.0, frozen=True, units='keV')
        self.ktbb = XSParameter(name, 'ktbb', 1.0, min=0.1, max=1000.0, hard_min=0.1, hard_max=1000.0, frozen=True, units='keV')
        self.model = XSParameter(name, 'model', 1, alwaysfrozen=True)
        self.redshift = XSParameter(name, 'redshift', 2.0, min=0.01, max=10.0, hard_min=0.001, hard_max=10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, min=0.0, max=1e+24, hard_min=0.0, hard_max=1e+24)
        XSAdditiveModel.__init__(self, name, (self.thobs, self.thjet, self.gamma, self.r12,
                                              self.p1, self.p2, self.E0, self.delta,
                                              self.index_pl, self.ecut, self.ktbb,
                                              self.model, self.redshift, self.norm))


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
    XSgrbcomp, XSgrbjet

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGrbm.html

    """

    __function__ = "xsgrbm"

    def __init__(self, name='grbm'):
        self.alpha = XSParameter(name, 'alpha', -1., -3., +2., -10, 5)
        self.beta = XSParameter(name, 'beta', -2., -5., +2., -10, 10)
        self.tem = XSParameter(name, 'tem', +300., +50., +1000., 10, 10000, units='keV', aliases=["temp"])
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.alpha, self.beta, self.tem, self.norm))


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
        self.T = XSParameter(name, 'T', 3.0, 0.5, 10.0, 0.5, 10.0, units='MK')
        self.NSmass = XSParameter(name, 'NSmass', 1.4, 0.6, 2.8, 0.6, 2.8,
                                  units='Msun')
        self.NSrad = XSParameter(name, 'NSrad', 10.0, 5.0, 23.0, 5.0, 23.0, units='km')
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
        self.mass = XSParameter(name, 'mass', 1e9, 1., 1e10, 1., 1e10,
                                units='solar', frozen=True)
        self.Dco = XSParameter(name, 'Dco', 3350.6, 1., 1e8, 1., 1e8,
                               units='Mpc', frozen=True)
        self.log_mdot = XSParameter(name, 'log_mdot', -1., -5., 2., -5., 2.,
                                    units='logL/LEdd')
        self.thetaobs = XSParameter(name, 'thetaobs', 3., 0., 90., 0., 90.,
                                    units='deg', frozen=True)
        self.BulkG = XSParameter(name, 'BulkG', 13., 1., 100., 1., 100, frozen=True)
        self.phi = XSParameter(name, 'phi', 0.1, 1e-2, 1e2, 1e-2, 1e2,
                               units='rad', frozen=True)
        self.zdiss = XSParameter(name, 'zdiss', 1275., 10., 1e4, 10., 1e4,
                                 units='Rg', frozen=True)
        self.B = XSParameter(name, 'B', 2.6, 1e-2, 15., 1e-2, 15.,
                             units='Gau', frozen=True)
        self.logPrel = XSParameter(name, 'logPrel', 43.3, 40., 48., 40., 48., frozen=True)
        self.gmin_inj = XSParameter(name, 'gmin_inj', 1.0, 1., 1e3, 1., 1e3, frozen=True)
        self.gbreak = XSParameter(name, 'gbreak', 300., 10., 1e4, 10., 1e4, frozen=True)
        self.gmax = XSParameter(name, 'gmax', 3e3, 1e3, 1e6, 1e3, 1e6, frozen=True)
        self.s1 = XSParameter(name, 's1', 1., -1., 1., -1., 1., frozen=True)
        self.s2 = XSParameter(name, 's2', 2.7, 1., 5., 1., 5., frozen=True)
        self.z = XSParameter(name, 'z', 0.0, 0., 10., 0., 10., frozen=True)
        # Note: alwaysfrozen is set for norm based on the documentation,
        # since there's no way to determine this from the model.dat file
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval,
                                alwaysfrozen=True)
        XSAdditiveModel.__init__(self, name,
                                 (self.mass, self.Dco, self.log_mdot, self.thetaobs,
                                  self.BulkG, self.phi, self.zdiss, self.B, self.logPrel,
                                  self.gmin_inj, self.gbreak, self.gmax, self.s1,
                                  self.s2, self.z, self.norm))


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
        black hole mass M (when G=c=1). It should be in the range [-1, 1).
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
    XSgrad, XSzkerrbb

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKerrbb.html

    """

    __function__ = "C_kerrbb"

    def __init__(self, name='kerrbb'):
        self.eta = XSParameter(name, 'eta', 0., 0., 1.0, 0.0, 1, frozen=True)
        self.a = XSParameter(name, 'a', 0., -1., 0.9999, -1, 0.9999)
        self.i = XSParameter(name, 'i', 30., 0., 85., 0.0, 85, units='deg', frozen=True)
        self.Mbh = XSParameter(name, 'Mbh', 1., 0., 100., 0.0, 100, units='Msun')
        self.Mdd = XSParameter(name, 'Mdd', 1., 0., 1000., 0.0, 1000, units='Mdd0')
        self.Dbh = XSParameter(name, 'Dbh', 10., 0., 10000., 0.0, 10000, units='kpc', frozen=True)
        self.hd = XSParameter(name, 'hd', 1.7, 1., 10., 1, 10, frozen=True)
        self.rflag = XSParameter(name, 'rflag', 1., alwaysfrozen=True)
        self.lflag = XSParameter(name, 'lflag', 0., alwaysfrozen=True)
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
        self.distance = XSParameter(name, 'distance', 1., 0.01, 1000., 0.01, 1000, units='kpc', frozen=True)
        self.TcoloTeff = XSParameter(name, 'TcoloTeff', 1.5, 1.0, 2.0, 1.0, 2, frozen=True, aliases=["TcolTeff"])
        self.M = XSParameter(name, 'M', 1.0, 0.1, 100., 0.1, 100, units='solar')
        self.Mdot = XSParameter(name, 'Mdot', 1.0, 0.01, 100., 0.01, 100, units='1e18')
        self.Incl = XSParameter(name, 'Incl', 30., 0., 90., 0.0, 90, units='deg', frozen=True)
        self.Rin = XSParameter(name, 'Rin', 1.235, 1.235, 100., 1.235, 100, units='Rg', frozen=True)
        self.Rout = XSParameter(name, 'Rout', 1e5, 1e4, 1e8, 1e4, 1e8, units='Rg', frozen=True)
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
    XSdiskline, XSkerrconv, XSlaor

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKerrdisk.html

    """

    __function__ = "spin"

    def __init__(self, name='kerrdisk'):
        self.lineE = XSParameter(name, 'lineE', 6.4, 0.1, 100., 0.1, 100, units='keV', frozen=True)
        self.Index1 = XSParameter(name, 'Index1', 3., -10., 10., -10, 10, frozen=True)
        self.Index2 = XSParameter(name, 'Index2', 3., -10., 10., -10, 10, frozen=True)
        self.r_br_g = XSParameter(name, 'r_br_g', 6.0, 1.0, 400., 1.0, 400, frozen=True, aliases=["r_brg"])
        self.a = XSParameter(name, 'a', 0.998, 0.01, 0.998, 0.01, 0.998)
        self.Incl = XSParameter(name, 'Incl', 30., 0., 90., 0.0, 90, units='deg', frozen=True)
        self.Rin_ms = XSParameter(name, 'Rin_ms', 1.0, 1.0, 400., 1.0, 400, frozen=True, aliases=["Rinms"])
        self.Rout_ms = XSParameter(name, 'Rout_ms', 400., 1.0, 400., 1.0, 400, frozen=True, aliases=["Routms"])
        self.z = XSParameter(name, 'z', 0., 0., 10., 0.0, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.lineE, self.Index1, self.Index2, self.r_br_g, self.a, self.Incl, self.Rin_ms, self.Rout_ms, self.z, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.lineE, **kwargs)


@version_at_least("12.10.1")
class XSkyconv(XSConvolutionKernel):
    """The XSPEC kyconv convolution model: convolution using a relativistic line from axisymmetric accretion disk

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    a
        a/M, the black-hole angular momentum in (GM/c).
    theta_o
        The observer inclination in degrees (0 is pole on).
    rin
        The inner radius, in GM/c^2.
    ms
        Flag: 0 means integrate from rin, 1 means integrate the emission
        from above the marginally-stable orbit only.
    rout
        The outer radius, in GM/c^2.
    alpha
        The accretion disk emissivity scales as r^-alpha for r < rb.
    beta
        The accretion disk emissivity scales as r^-beta for r > rb.
    rb
        The boundary radius between inner and outer emissivity laws,
        in units of GM/c^2.
    zshift
        The  overall Doppler shift.
    limb
        Flag: 0 means isotropic emission, 1 means Laor's limb darkening
        (1 + 0.26 \\mu), 2 means Haardt's limb brightening (log(1 + 1 / \\mu)).
    ne_loc
        The number of grid points in local energy (energy resolution of
        local flux, the grid is equidistant in logarithmic scale).
    normal
        Flag: 0 means normalize total flux to unity, > 0 means normalize
        to unity at the energy given by the parameter value, and
        < 0 means unnormalized.

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    This model is only available when used with XSPEC 12.10.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKyconv.html

    """

    __function__ = "kyconv"

    def __init__(self, name='xskyconv'):
        self.a = XSParameter(name, 'a', 0.9982, min=0.0, max=1.0,
                             hard_min=0.0, hard_max=1.0, frozen=False,
                             units='GM/c')
        self.theta_o = XSParameter(name, 'theta_o', 30.0, min=0.0, max=89.0,
                                   hard_min=0.0, hard_max=89.0, frozen=False,
                                   units='deg')
        self.rin = XSParameter(name, 'rin', 1.0, min=1.0, max=1000.0,
                               hard_min=1.0, hard_max=1000.0, frozen=True,
                               units='GM/c^2')
        self.ms = XSParameter(name, 'ms', 1.0, min=0.0, max=1.0,
                              hard_min=0.0, hard_max=1.0, frozen=True)
        self.rout = XSParameter(name, 'rout', 400.0, min=1.0, max=1000.0,
                                hard_min=1.0, hard_max=1000.0, frozen=True,
                                units='GM/c^2')
        self.alpha = XSParameter(name, 'alpha', 3.0, min=-20.0, max=20.0,
                                 hard_min=-20.0, hard_max=20.0, frozen=True)
        self.beta = XSParameter(name, 'beta', 3.0, min=-20.0, max=20.0,
                                hard_min=-20.0, hard_max=20.0, frozen=True)
        self.rb = XSParameter(name, 'rb', 400.0, min=1.0, max=1000.0,
                              hard_min=1.0, hard_max=1000.0, frozen=True,
                              units='GM/c^2')
        self.zshift = XSParameter(name, 'zshift', 0.0, min=-0.999, max=10.0,
                                  hard_min=-0.999, hard_max=10.0, frozen=True)
        self.limb = XSParameter(name, 'limb', 0.0, min=0.0, max=2.0,
                                hard_min=0.0, hard_max=2.0, frozen=True)
        self.ne_loc = XSParameter(name, 'ne_loc', 100.0, min=3.0, max=5000.0,
                                  hard_min=3.0, hard_max=5000.0, frozen=True)
        self.normal = XSParameter(name, 'normal', 1.0, min=-1.0, max=100.0,
                                  hard_min=-1.0, hard_max=100.0, frozen=True)

        pars = (self.a, self.theta_o, self.rin, self.ms, self.rout,
                self.alpha, self.beta, self.rb, self.zshift, self.limb,
                self.ne_loc, self.normal)
        XSConvolutionKernel.__init__(self, name, pars)


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
        self.a = XSParameter(name, 'a', 0.9982, 0, 1, 0, 1, units='GM/c')
        self.theta_o = XSParameter(name, 'theta_o', 30, 0, 89, 0, 89, units='deg')
        self.rin = XSParameter(name, 'rin', 1, 1, 1000, 1, 1000, units='GM/c^2',
                               frozen=True)
        self.ms = XSParameter(name, 'ms', 1, 0, 1, 0, 1, alwaysfrozen=True)
        self.rout = XSParameter(name, 'rout', 400, 1, 1000, 1, 1000, units='GM/c^2',
                                frozen=True)
        self.Erest = XSParameter(name, 'Erest', 6.4, 1, 99, 1, 99, units='keV',
                                 frozen=True)
        self.alpha = XSParameter(name, 'alpha', 3, -20, 20, -20, 20, frozen=True)
        self.beta = XSParameter(name, 'beta', 3, -20, 20, -20, 20, frozen=True)
        self.rb = XSParameter(name, 'rb', 400, 1, 1000, 1, 1000, units='GM/c^2',
                              frozen=True)
        self.zshift = XSParameter(name, 'zshift', 0, -0.999, 10, -0.999, 10,
                                  frozen=True)
        self.limb = XSParameter(name, 'limb', 1, 0, 2, 0, 2, frozen=True)
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
    XSkdblur, XSlaor2

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLaor.html

    """

    __function__ = "C_laor" if equal_or_greater_than("12.10.1") else "C_xslaor"

    def __init__(self, name='laor'):
        self.lineE = XSParameter(name, 'lineE', 6.4, 0., 100., 0.0, 100, units='keV')
        self.Index = XSParameter(name, 'Index', 3., -10., 10., -10, 10, frozen=True)
        self.Rin_G = XSParameter(name, 'Rin_G', 1.235, 1.235, 400., 1.235, 400, frozen=True, aliases=["RinG"])
        self.Rout_G = XSParameter(name, 'Rout_G', 400., 1.235, 400., 1.235, 400, frozen=True, aliases=["RoutG"])
        self.Incl = XSParameter(name, 'Incl', 30., 0., 90., 0.0, 90, units='deg', frozen=True)
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
    XSkdblur2, XSlaor

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLaor2.html

    """

    __function__ = "C_laor2"

    def __init__(self, name='laor2'):
        self.lineE = XSParameter(name, 'lineE', 6.4, 0., 100., 0.0, 100, units='keV')
        self.Index = XSParameter(name, 'Index', 3., -10., 10., -10, 10, frozen=True)
        self.Rin_G = XSParameter(name, 'Rin_G', 1.235, 1.235, 400., 1.235, 400, frozen=True, aliases=["RinG"])
        self.Rout_G = XSParameter(name, 'Rout_G', 400., 1.235, 400., 1.235, 400, frozen=True, aliases=["RoutG"])
        self.Incl = XSParameter(name, 'Incl', 30., 0., 90., 0.0, 90, units='deg', frozen=True)
        self.Rbreak = XSParameter(name, 'Rbreak', 20., 1.235, 400., 1.235, 400, frozen=True)
        self.Index1 = XSParameter(name, 'Index1', 3., -10., 10., -10, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.lineE, self.Index, self.Rin_G, self.Rout_G, self.Incl, self.Rbreak, self.Index1, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.lineE, **kwargs)


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
        self.alpha = XSParameter(name, 'alpha', 1.5, 0., 4., 0.0, 4)
        self.beta = XSParameter(name, 'beta', 0.2, -4., 4., -4, 4)
        self.pivotE = XSParameter(name, 'pivotE', 1.0, units='keV',
                                  alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.alpha, self.beta, self.pivotE, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


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
        self.LineE = XSParameter(name, 'LineE', 6.5, 0., 1.e6, 0.0, 1e6, units='keV')
        self.Width = XSParameter(name, 'Width', 0.1, 0., 10., 0.0, 20, units='keV')
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
        self.kT = XSParameter(name, 'kT', 1., 1.e-3, 1.e2, 1e-3, 1e2, units='keV')
        self.nH = XSParameter(name, 'nH', 1., 1.e-5, 1.e19, 1e-6, 1e20, units='cm-3', frozen=True)
        self.Abundanc = XSParameter(name, 'Abundanc', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.kT = XSParameter(name, 'kT', 1., 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.nH = XSParameter(name, 'nH', 1., 1.e-5, 1.e19, 1e-6, 1e20, units='cm-3', frozen=True)
        self.Abundanc = XSParameter(name, 'Abundanc', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 1, 0, 1, 0, 1, alwaysfrozen=True)
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
        The maximum temperature, in keV.
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
        self.lowT = XSParameter(name, 'lowT', 0.1, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.highT = XSParameter(name, 'highT', 4., 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1., 0., 5., 0.0, 5)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 1, 0, 1, 0, 1, alwaysfrozen=True)
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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 1e8, 5e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.Tau, self.redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 1., 0.008, 64.0, 0.008, 64, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1., 0., 5., 0.0, 5, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.redshift, self.norm))


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
        self.kT_a = XSParameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_b = XSParameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0100, 79.9, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau_l = XSParameter(name, 'Tau_l', 0.0, 0., 5.e13, 0.0, 5e13, units='s/cm^3', frozen=True)
        self.Tau_u = XSParameter(name, 'Tau_u', 1.e11, 1.e8, 5.e13, 1e8, 5e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.LogT_eff = XSParameter(name, 'LogT_eff', 6.0, 5.0, 7.0, 5, 7, units='K')
        self.M_ns = XSParameter(name, 'M_ns', 1.4, 0.5, 2.5, 0.5, 2.5, units='Msun')
        self.R_ns = XSParameter(name, 'R_ns', 10.0, 5.0, 20., 5.0, 20, units='km')
        self.MagField = XSParameter(name, 'MagField', 0.0, 0.0, 5.e13, 0.0, 5e13, units='G', frozen=True)
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
        self.LogT_eff = XSParameter(name, 'LogT_eff', 6.0, 5.5, 6.5, 5.5, 6.5, units='K')
        self.NSmass = XSParameter(name, 'NSmass', 1.4, 0.3, 2.5, 0.3, 2.5, units='Msun')
        self.NSrad = XSParameter(name, 'NSrad', 10.0, 6.0, 20., 6.0, 20, units='km')
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
        self.LogT_eff = XSParameter(name, 'LogT_eff', 6.0, 5.0, 6.5, 5.0, 6.5, units='K')
        self.M_ns = XSParameter(name, 'M_ns', 1.4, 0.5, 3.0, 0.5, 3, units='Msun')
        self.R_ns = XSParameter(name, 'R_ns', 10.0, 5.0, 30., 5.0, 30, units='km')
        self.dist = XSParameter(name, 'dist', 10.0, 0.1, 100.0, 0.1, 100, units='kpc')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LogT_eff, self.M_ns, self.R_ns, self.dist, self.norm))


class XSnsmax(XSAdditiveModel):
    """The XSPEC nsmax model: Neutron Star Magnetic Atmosphere.

    The model is described at [1]_. It has been superseded by
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
        self.logTeff = XSParameter(name, 'logTeff', 6.0, 5.5, 6.8, 5.5, 6.8, units='K')
        self.redshift = XSParameter(name, 'redshift', 0.1, 1.0e-5, 1.5, 1.0e-5, 2.0)
        self.specfile = XSParameter(name, 'specfile', 1200, alwaysfrozen=True)
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
        self.logTeff = XSParameter(name, 'logTeff', 6.0, 5.5, 6.9, 5.5, 6.9, units='K')
        self.M_ns = XSParameter(name, 'M_ns', 1.4, 0.5, 3.0, 0.5, 3.0, units='Msun')
        self.R_ns = XSParameter(name, 'R_ns', 10.0, 5.0, 30.0, 5.0, 30.0, units='km')
        self.dist = XSParameter(name, 'dist', 1.0, 0.01, 100.0, 0.01, 100.0, units='kpc')
        self.specfile = XSParameter(name, 'specfile', 1200, alwaysfrozen=True)
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
        self.logTeff = XSParameter(name, 'logTeff', 6.0, 5.5, 6.7, 5.5, 6.7, units='K')
        self.M_ns = XSParameter(name, 'M_ns', 1.4, 0.5, 3.0, 0.5, 3.0, units='Msun')
        self.R_ns = XSParameter(name, 'R_ns', 10.0, 5.0, 30.0, 5.0, 30.0, units='km')
        self.dist = XSParameter(name, 'dist', 1.0, 0.01, 100.0, 0.01, 100.0, units='kpc')
        self.specfile = XSParameter(name, 'specfile', 6, alwaysfrozen=True)
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
        self.l_nth = XSParameter(name, 'l_nth', 100., 0., 1.e4, 0.0, 1e4)
        self.l_bb = XSParameter(name, 'l_bb', 100., 0., 1.e4, 0.0, 1e4)
        self.f_refl = XSParameter(name, 'f_refl', 0., 0., 4., 0.0, 4)
        self.kT_bb = XSParameter(name, 'kT_bb', 10., 1., 100., 1.0, 100, frozen=True)
        self.g_max = XSParameter(name, 'g_max', 1.e3, 5., 1.e4, 5.0, 1e4, frozen=True)
        self.l_th = XSParameter(name, 'l_th', 0., 0., 1.e4, 0.0, 1e4, frozen=True)
        self.tau_p = XSParameter(name, 'tau_p', 0., 0., 10., 0.0, 10, frozen=True)
        self.G_inj = XSParameter(name, 'G_inj', 0., 0., 5., 0.0, 5, frozen=True)
        self.g_min = XSParameter(name, 'g_min', 1.3, 1., 1.e3, 1, 1e3, frozen=True)
        self.g_0 = XSParameter(name, 'g_0', 1.3, 1., 5., 1.0, 5, frozen=True)
        self.radius = XSParameter(name, 'radius', 1.e13, 1.e5, 1.e16, 1e5, 1e16, frozen=True)
        self.pair_esc = XSParameter(name, 'pair_esc', 0., 0., 1., 0.0, 1, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.45, 0.05, 0.95, 0.05, 0.95)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0.1, 10., 0.1, 10, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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

    See Also
    --------
    XSthcomp

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNthcomp.html

    """

    __function__ = "C_nthcomp"

    def __init__(self, name='nthcomp'):
        self.Gamma = XSParameter(name, 'Gamma', 1.7, 1.001, 5., 1.001, 10.)
        self.kT_e = XSParameter(name, 'kT_e', 100., 5., 1.e3, 1., 1.e3, units='keV')
        self.kT_bb = XSParameter(name, 'kT_bb', 0.1, 1.e-3, 10., 1.e-3, 10., units='keV', frozen=True)
        self.inp_type = XSParameter(name, 'inp_type', 0., 0., 1., 0., 1., units='0/1', frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10., frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Gamma, self.kT_e, self.kT_bb, self.inp_type, self.redshift, self.norm))


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
        self.mass = XSParameter(name, 'mass', 1e7, 1.0, 1.e9, 1.0, 1e9, units='solar', frozen=True)
        self.dist = XSParameter(name, 'dist', 100, 0.01, 1.e9, 0.01, 1e9, units='Mpc', frozen=True)
        self.logLoLEdd = XSParameter(name, 'logLoLEdd', -1., -10., 2, -10, 2, aliases=["logLLEdd"])
        self.astar = XSParameter(name, 'astar', 0.0, 0., 0.998, 0.0, 0.998, frozen=True)
        self.rcor = XSParameter(name, 'rcor', 10.0, 1., 100., 1.0, 100, units='rg')
        self.logrout = XSParameter(name, 'logrout', 5.0, 3.0, 7.0, 3.0, 7, frozen=True)
        self.kT_e = XSParameter(name, 'kT_e', 0.2, 0.01, 10, 0.01, 10, units='keV')
        self.tau = XSParameter(name, 'tau', 10., 0.1, 100, 0.1, 100)
        self.Gamma = XSParameter(name, 'Gamma', 2.1, 0.5, 5., 0.5, 10)
        self.fpl = XSParameter(name, 'fpl', 1.e-4, 0.0, 1.e-1, 0.0, 1e-1)
        self.fcol = XSParameter(name, 'fcol', 2.4, 1.0, 5, 1.0, 5, frozen=True)
        self.tscat = XSParameter(name, 'tscat', 1.e5, 1e4, 1e5, 1e4, 1e5, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0., 0., 10., 0.0, 10, frozen=True)
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
        self.mass = XSParameter(name, 'mass', 1e7, 1.0, 1.e9, 1.0, 1e9, units='solar', frozen=True)
        self.dist = XSParameter(name, 'dist', 100, 0.01, 1.e9, 0.01, 1e9, units='Mpc', frozen=True)
        self.logLoLEdd = XSParameter(name, 'logLoLEdd', -1., -10., 2, -10, 2, aliases=["logLLEdd"])
        self.astar = XSParameter(name, 'astar', 0.0, 0., 0.998, 0.0, 0.998, frozen=True)
        self.rcor = XSParameter(name, 'rcor', 10.0, 1., 100., 1.0, 100, units='rg')
        self.logrout = XSParameter(name, 'logrout', 5.0, 3.0, 7.0, 3.0, 7, frozen=True)
        self.kT_e = XSParameter(name, 'kT_e', 0.2, 0.01, 10, 0.01, 10, units='keV')
        self.tau = XSParameter(name, 'tau', 10., 0.1, 100, 0.1, 100)
        self.Gamma = XSParameter(name, 'Gamma', 2.1, 1.05, 5., 1.05, 10.0)
        self.fpl = XSParameter(name, 'fpl', 1.e-4, 0.0, 1., 0.0, 1)
        self.Redshift = XSParameter(name, 'Redshift', 0., 0., 10., 0.0, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.mass, self.dist, self.logLoLEdd, self.astar, self.rcor, self.logrout, self.kT_e, self.tau, self.Gamma, self.fpl, self.Redshift, self.norm))


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
        self.PhoIndex = XSParameter(name, 'PhoIndex', 1., -2., 9., -3, 10)
        self.eMin = XSParameter(name, 'eMin', 2., -100., 1.e10, -100, 1e10, units='keV', frozen=True)
        self.eMax = XSParameter(name, 'eMax', 10., -100., 1.e10, -100, 1e10, units='keV', frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.eMin, self.eMax, self.norm))


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
        The abundance of the elements heavier than He relative to their
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
        self.PhoIndex = XSParameter(name, 'PhoIndex', 2., 1.1, 2.5, 1.1, 2.5)
        self.foldE = XSParameter(name, 'foldE', 1000., 1., 1.e6, 1.0, 1e6, units='keV', frozen=True)
        self.rel_refl = XSParameter(name, 'rel_refl', -1, -1.e6, 1.e6, -1e6, 1e6, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., 0., 4., 0.0, 4, frozen=True)
        self.abund = XSParameter(name, 'abund', 1., 0.0, 1.e6, 0.0, 1e6, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0.0, 100., 0.0, 100, frozen=True)
        self.Incl = XSParameter(name, 'Incl', 60., 0., 85.0, 0.0, 85, units='deg')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.foldE, self.rel_refl, self.redshift, self.abund, self.Fe_abund, self.Incl, self.norm))


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
        The abundance of the elements heavier than He relative to their
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
    XSpexriv, XSpexmon, XSreflect

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
        self.PhoIndex = XSParameter(name, 'PhoIndex', 2., -9., 9., -10, 10)
        self.foldE = XSParameter(name, 'foldE', 100., 1., 1.e6, 1.0, 1e6, units='keV')
        self.rel_refl = XSParameter(name, 'rel_refl', 0., 0., 1.e6, 0, 1e6)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.abund = XSParameter(name, 'abund', 1., 0.0, 1.e6, 0.0, 1e6, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0.0, 1.e6, 0.0, 1e6, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.45, 0.05, 0.95, 0.05, 0.95, frozen=True)
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
        The abundance of the elements heavier than He relative to their
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
    XSireflect, XSpexrav, XSpexmon

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
        self.PhoIndex = XSParameter(name, 'PhoIndex', 2., -9., 9., -10, 10)
        self.foldE = XSParameter(name, 'foldE', 100., 1., 1.e6, 1.0, 1e6, units='keV')
        self.rel_refl = XSParameter(name, 'rel_refl', 0., 0., 1.e6, 0, 1e6)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.abund = XSParameter(name, 'abund', 1., 0.0, 1.e6, 0.0, 1e6, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0.0, 1.e6, 0.0, 1e6, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.45, 0.05, 0.95, 0.05, 0.95, frozen=True)
        self.T_disk = XSParameter(name, 'T_disk', 3.e4, 1.e4, 1.e6, 1e4, 1e6, units='K', frozen=True)
        self.xi = XSParameter(name, 'xi', 1., 0., 1.e3, 0.0, 5e3, units='erg cm/s')
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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.nmax = XSParameter(name, 'nmax', 1, alwaysfrozen=True)
        self.FeAbun = XSParameter(name, 'FeAbun', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.FeKedge = XSParameter(name, 'FeKedge', 7.11, 7., 10., 7.0, 10, units='KeV', frozen=True)
        self.PhoIndex = XSParameter(name, 'PhoIndex', 2., -2., 9., -3, 10)
        self.HighECut = XSParameter(name, 'HighECut', 95., 1., 100., 0.01, 200., units='keV', frozen=True)
        self.foldE = XSParameter(name, 'foldE', 100., 1., 1.e6, 1.0, 1e6, frozen=True)
        self.acrit = XSParameter(name, 'acrit', 1., 0.0, 1.0, 0.0, 1, frozen=True)
        self.FAST = XSParameter(name, 'FAST', 0, alwaysfrozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.PhoIndex = XSParameter(name, 'PhoIndex', 1., -2., 9., -3, 10)
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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau_l = XSParameter(name, 'Tau_l', 0.0, 0., 5.e13, 0.0, 5e13, units='s/cm^3', frozen=True)
        self.Tau_u = XSParameter(name, 'Tau_u', 1.e11, 1.e8, 5.e13, 1e8, 5e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.mass = XSParameter(name, 'mass', 1e7, 1e5, 1e10, 1e5, 1e10,
                                'solar', frozen=True)
        self.dist = XSParameter(name, 'dist', 100, 0.01, 1e9, 0.01, 1e9,
                                'Mpc', frozen=True)
        self.logmdot = XSParameter(name, 'logmdot', -1, -1.65, 0.39, -1.65, 0.39, units='Ledd')
        self.astar = XSParameter(name, 'astar', 0.0, -1, 0.998, -1, 0.998,
                                 frozen=True)
        self.cosi = XSParameter(name, 'cosi', 0.5, 0.05, 1.0, 0.05, 1.0,
                                frozen=True)

        self.redshift = XSParameter(name, 'redshift', 0, 0, 5, 0, 5,
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
        self.kT = XSParameter(name, 'kT', 1., 0.008, 64.0, 0.008, 64, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1., 0., 5., 0.0, 5, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.edge = XSParameter(name, 'edge', 1.4, 0.001, 100., 0.001, 100, units='keV')
        self.kT = XSParameter(name, 'kT', 1., 0.001, 100., 0.001, 100, units='keV')
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
        The abundance of the elements heavier than He relative to their
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
        self.PhoIndex = XSParameter(name, 'PhoIndex', 2., -9., 9., -10, 10)
        self.foldE = XSParameter(name, 'foldE', 100., 1., 1.e6, 1.0, 1e6, units='keV')
        self.rel_refl = XSParameter(name, 'rel_refl', 0., 0., 2., 0.0, 2)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.abund = XSParameter(name, 'abund', 1., 0.5, 10., 0.5, 10, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0.1, 10., 0.1, 10, frozen=True)
        self.Incl = XSParameter(name, 'Incl', 30., 19., 87., 19, 87, units='deg', frozen=True)
        self.T_disk = XSParameter(name, 'T_disk', 3.e4, 1.e4, 1.e6, 1e4, 1e6, units='K', frozen=True)
        self.xi = XSParameter(name, 'xi', 1., 0., 1.e3, 0.0, 5e3, units='ergcm/s')
        self.Betor10 = XSParameter(name, 'Betor10', -2., -10., 20., -10, 20, frozen=True)
        self.Rin = XSParameter(name, 'Rin', 10., 6., 1000., 6.0, 10000, units='R_g', frozen=True)
        self.Rout = XSParameter(name, 'Rout', 1000., 0., 1000000., 0.0, 10000000, units='R_g', frozen=True)
        self.accuracy = XSParameter(name, 'accuracy', 30., 30., 100000., 30.0, 100000, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.foldE, self.rel_refl, self.redshift, self.abund, self.Fe_abund, self.Incl, self.T_disk, self.xi, self.Betor10, self.Rin, self.Rout, self.accuracy, self.norm))


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
        self.kT = XSParameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = XSParameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.Abundanc, self.Tau, self.redshift, self.norm))


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
        self.kT_a = XSParameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_b = XSParameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0100, 79.9, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 1e8, 5e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.tin = XSParameter(name,    'tin', 1., 0.01, 100., 0.01, 1000., units='keV')
        self.rin = XSParameter(name,    'rin', 1.e-2, 1.e-5, 1.0, 1.e-6, 10, units='rsph')
        self.rout = XSParameter(name,   'rout', 100., 0.1, 1.e8, 0.1, 1.e8, units='rsph')
        self.theta = XSParameter(name,  'theta', 22.9, 1., 89., 0., 90., units='deg')
        self.incl = XSParameter(name,   'incl', 0., -90., 90., -90., 90., units='deg', frozen=True)
        self.valpha = XSParameter(name, 'valpha', -0.5, -1.0, 2., -1.5, 5., frozen=True)
        self.gamma = XSParameter(name,  'gamma', 1.333, 0.5, 10., 0.5, 10., frozen=True)
        self.mdot = XSParameter(name,   'mdot', 1000., 0.5, 1.e7, 0.5, 1.e7, frozen=True)
        self.irrad = XSParameter(name,  'irrad', 2., 0., 10., 0., 20., frozen=True)
        self.norm = Parameter(name,   'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.tin, self.rin, self.rout, self.theta,
                                              self.incl, self.valpha, self.gamma, self.mdot,
                                              self.irrad, self.norm))


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
        A flag to control the surface profile. If greater than zero,
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
        self.M = XSParameter(name, 'M', 10.0, 0.0, 1000.0, 0.0, 1000.0, units='Msun',
                             frozen=True)
        self.a = XSParameter(name, 'a', 0.0, 0.0, 0.999, 0.0, 0.999, units='GM/c')
        self.lumin = XSParameter(name, 'lumin', 0.5, 0.05, 1.0, 0.05, 1.0,
                                 units='L_Edd')
        self.alpha = XSParameter(name, 'alpha', 0.1, 0.005, 0.1, 0.005, 0.1,
                                 frozen=True)
        self.inc = XSParameter(name, 'inc', 60.0, 0.0, 85.0, 0.0, 85.0,
                               units='deg', frozen=True)
        self.D = XSParameter(name, 'D', 10.0, 0.0, 1e4, 0.0, 1e4, units='kpc',
                             frozen=True)
        self.f_hard = XSParameter(name, 'f_hard', -1.0, -10.0, 10.0, -10.0, 10.0,
                                  frozen=True)
        self.lflag = XSParameter(name, 'lflag', 1.0, alwaysfrozen=True)
        self.vflag = XSParameter(name, 'vflag', 1.0, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval,
                                frozen=True)
        XSAdditiveModel.__init__(self, name,
                                 (self.M, self.a, self.lumin, self.alpha,
                                  self.inc, self.D, self.f_hard,
                                  self.lflag, self.vflag, self.norm))


# NOTE: we do not support the smaug model yet


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
        The percentage of SNIa.
    SNIModelIndex
        SNIa yield model: see [1]_ for more details.
    SNIIModelIndex
        SNII yield model: see [1]_ for more details.
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
        self.kT = XSParameter(name, 'kT', 1., 0.008, 64.0, 0.008, 64.0, units='keV')
        self.N_SNe = XSParameter(name, 'N_SNe', 1., 0.0, 1e20, 0, 1e20, units='10^9')
        # QUS: if this is a percentage, why is the maximum not 100?
        self.R = XSParameter(name, 'R', 1., 0.0, 1e20, 0, 1e20)
        self.SNIModelIndex = XSParameter(name, 'SNIModelIndex', 1.,
                                         0.0, 125, 0, 125, alwaysfrozen=True)
        self.SNIIModelIndex = XSParameter(name, 'SNIIModelIndex', 1.,
                                          0.0, 125, 0, 125, alwaysfrozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., 0.0, 10., 0.0, 10.0,
                                    frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.N_SNe, self.R,
                                  self.SNIModelIndex, self.SNIIModelIndex,
                                  self.redshift, self.norm))


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
        self.alpha = XSParameter(name, 'alpha', 0.5, 0.3, 0.8, 1e-5, 1)
        self.break_ = XSParameter(name, 'break_', 2.42E17, 1.E15, 1.E19, 1e10, 1e25, 'Hz', aliases=["breakfreq"])
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
        The 1 GHz flux in Jy.

    See Also
    --------
    XSsrcut

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSresc.html

    """

    __function__ = "sresc"

    def __init__(self, name='sresc'):
        self.alpha = XSParameter(name, 'alpha', 0.5, 0.3, 0.8, 1e-5, 1)
        self.rolloff = XSParameter(name, 'rolloff', 2.42E17, 1.E15, 1.E19, 1e10, 1e25, units='Hz')
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
        self.te = XSParameter(name, 'te', 0.1, 0.01, 0.5, 0.01, 0.5)
        self.y = XSParameter(name, 'y', 0.7, 1e-4, 1e3, 1e-4, 1e3)
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
        self.Energy = XSParameter(name, 'Energy', 6.5, 0., 100., 0.0, 100, units='keV')
        self.Sigma = XSParameter(name, 'Sigma', 0.1, 0., 10., 0.0, 20, units='keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Energy, self.Sigma, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.Energy, **kwargs)


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
        self.kT = XSParameter(name, 'kT', 1.0, 0.008, 64.0, 0.008, 64.0, units='keV')
        self.kTi = XSParameter(name, 'kTi', 1.0, 0.008, 64.0, 0.008, 64.0, units='keV')
        self.Abundanc = XSParameter(name, 'Abundanc', 1.0, 0.0, 5.0, 0.0, 5.0)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.Abundanc, self.Redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.kT = XSParameter(name, 'kT', 3.0, 1.e-2, 100., 1e-4, 200, units='keV')
        self.HeovrH = XSParameter(name, 'HeovrH', 1.0, 0., 100., 0.0, 100, aliases=["HeH"])
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kT, self.HeovrH, self.norm))


@version_at_least("12.10.1")
class XSvcph(XSAdditiveModel):
    """The XSPEC vcph model: Cooling + heating model for cool core clusters

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The default Redshift parameter value has changed from 0.1 to 0
       to match XSPEC.

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
        self.peakT = XSParameter(name, 'peakT', 2.2, 1e-1, 1e2, 1e-1, 1e2,
                                 units='keV')

        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 1000.0,
                              frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 1000.0,
                             frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 1000.0,
                             frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 1000.0,
                             frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 1000.0,
                              frozen=True)
        self.Na = XSParameter(name, 'Na', 1.0, 0., 1000., 0.0, 1000.0,
                              frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 1000.0,
                              frozen=True)
        self.Al = XSParameter(name, 'Al', 1.0, 0., 1000., 0.0, 1000.0,
                              frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 1000.0,
                              frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 1000.0,
                             frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 1000.0,
                              frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 1000.0,
                              frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 1000.0,
                              frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 1000.0,
                              frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, 0.0, 50, 0.0, 50,
                                    frozen=True)
        self.switch = XSParameter(name, 'switch', 1,
                                  alwaysfrozen=True)

        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        pars = (self.peakT,
                self.He, self.C, self.N, self.O, self.Ne, self.Na,
                self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca,
                self.Fe, self.Ni,
                self.Redshift, self.switch, self.norm)

        XSAdditiveModel.__init__(self, name, pars)


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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.norm))


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
        self.Tmean = XSParameter(name, 'Tmean', 4.0, 0.01, 10, 0.01, 20, units='keV', frozen=True)
        self.Tsigma = XSParameter(name, 'Tsigma', 0.1, 0.01, 1.e2, 0.01, 1e2, units='keV')
        self.nH = XSParameter(name, 'nH', 1.0, 1.e-5, 1.e19, 1e-6, 1e20, units='cm^-3', frozen=True)
        self.He = XSParameter(name, 'He', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Na = XSParameter(name, 'Na', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Al = XSParameter(name, 'Al', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 10., 0.0, 10, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.Tmean, self.Tsigma, self.nH, self.He, self.C, self.N, self.O,
                                              self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca,
                                              self.Fe, self.Ni, self.Redshift, self.switch, self.norm))


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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1.0, 0., 1., 0.0, 1, frozen=True)
        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 1e8, 5e13, units='s/cm^3')
        self.meankT = XSParameter(name, 'meankT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV', aliases=["kT_ave"])
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.meankT, self.redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 1., 1.e-3, 1.e2, 1e-3, 100, units='keV')
        self.nH = XSParameter(name, 'nH', 1., 1.e-5, 1.e19, 1e-6, 1e20, units='cm-3', frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.kT = XSParameter(name, 'kT', 1., 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.nH = XSParameter(name, 'nH', 1., 1.e-5, 1.e19, 1e-6, 1e20, units='cm-3', frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 1, alwaysfrozen=True)
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
        The maximum temperature, in keV.
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
        self.lowT = XSParameter(name, 'lowT', 0.1, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.highT = XSParameter(name, 'highT', 4., 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., 0.0, 10., 0, 10, frozen=True)
        self.switch = XSParameter(name, 'switch', 1, alwaysfrozen=True)
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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1.0, 0., 1., 0.0, 1, frozen=True)
        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 1e8, 5e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.redshift, self.norm))


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
        self.kT_a = XSParameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_b = XSParameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0100, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1.0, 0., 1., 0.0, 1, frozen=True)
        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau_l = XSParameter(name, 'Tau_l', 0.0, 0., 5.e13, 0.0, 5e13, units='s/cm^3', frozen=True)
        self.Tau_u = XSParameter(name, 'Tau_u', 1.e11, 1.e8, 5.e13, 1e8, 5e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau_l, self.Tau_u, self.redshift, self.norm))


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

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelVoigt.html

    """

    __function__ = "C_voigtLine"

    def __init__(self, name='voigt'):
        self.LineE = XSParameter(name, 'LineE', 6.5, 0., 1.e6, 0.0, 1e6, units='keV')
        self.Sigma = XSParameter(name, 'Sigma', 0.01, 0., 10., 0.0, 20, units='keV')
        self.Gamma = XSParameter(name, 'Gamma', 0.01, 0., 10., 0.0, 20, units='keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.LineE, self.Sigma, self.Gamma,
                                  self.norm))


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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1.0, 0., 1., 0.0, 1, frozen=True)
        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau_l = XSParameter(name, 'Tau_l', 0.0, 0., 5.e13, 0.0, 5e13, units='s/cm^3', frozen=True)
        self.Tau_u = XSParameter(name, 'Tau_u', 1.e11, 1.e8, 5.e13, 1e8, 5e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau_l, self.Tau_u, self.redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 6.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = XSParameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1.0, 0., 1., 0.0, 1.0, frozen=True)
        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 10000.0, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.redshift, self.norm))


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
        self.kT_a = XSParameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_b = XSParameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0100, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1.0, 0., 1., 0.0, 1, frozen=True)
        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 1e8, 5e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.kTi = XSParameter(name, 'kTi', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.He = XSParameter(name, 'He', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, 0., 1000., 0.0, 1000, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.He, self.C, self.N, self.O,
                                  self.Ne, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca,
                                  self.Fe, self.Ni, self.Redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.meankT = XSParameter(name, 'meankT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV', aliases=["kT_ave"])
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.meankT, self.redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.redshift, self.norm))


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
        self.kT_a = XSParameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_b = XSParameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0100, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Tau_l = XSParameter(name, 'Tau_l', 0.0, 0.0, 5.0e13, 0.0, 5.0e13, units='s/cm^3', frozen=True)
        self.Tau_u = XSParameter(name, 'Tau_u', 1.0e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau_l, self.Tau_u, self.redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Tau_l = XSParameter(name, 'Tau_l', 0.0, 0.0, 5.0e13, 0.0, 5.0e13, units='s/cm^3', frozen=True)
        self.Tau_u = XSParameter(name, 'Tau_u', 1.0e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10.0, -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau_l, self.Tau_u, self.redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = XSParameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.redshift, self.norm))


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
        self.kT_a = XSParameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_b = XSParameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.01, 79.9, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, 1000.0, frozen=True)
        self.Tau = XSParameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.kTi = XSParameter(name, 'kTi', 6.5, 0.0808, 68.447, 0.0808, 68.447, units='keV')
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Li = XSParameter(name, 'Li', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Be = XSParameter(name, 'Be', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.B = XSParameter(name, 'B', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.F = XSParameter(name, 'F', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.P = XSParameter(name, 'P', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.K = XSParameter(name, 'K', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.V = XSParameter(name, 'V', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, -0.999, 10, -0.999, 10)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name,
                                 (self.kT, self.kTi, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Redshift, self.norm))


@version_at_least("12.12.0")
class XSvvwdem(XSAdditiveModel):
    """The XSPEC vvwdem model: plasma emission, multi-temperature with power-law distribution of emission measure.

    The model is described at [1]_.

    Attributes
    ----------
    Tmax
        The maximum temperature for power-law emission measure
        distribution.
    beta
        The ratio of minimum to maxmum temperature.
    inv_slope
        The inverse of the slope (labelled p in the XSPEC documentation).
    nH
        Fixed at 1 for most applications.
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
        The abundance of the element in solar units.
    Redshift
        The redshift of the plasma.
    switch
        What model to use: 0 calculates with MEKAL, 1 interpolates
        with MEKAL, and 2 interpoates with APEC.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSvwdem, XSwdem

    Notes
    -----
    This model is only available when used with XSPEC 12.12.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelWdem.html

    """
    __function__ = "C_vvwDem"

    def __init__(self, name='vvwdem'):
        self.Tmax = XSParameter(name, 'Tmax', 1.0, min=0.01, max=10.0, hard_min=0.01, hard_max=20.0, units='keV')
        self.beta = XSParameter(name, 'beta', 0.1, min=0.01, max=1.0, hard_min=0.01, hard_max=1.0)
        # can not use p for the name as it conflicts with P
        self.inv_slope = XSParameter(name, 'inv_slope', 0.25, min=-1.0, max=10.0, hard_min=-1.0, hard_max=10.0)
        self.nH = XSParameter(name, 'nH', 1.0, min=1e-05, max=1e+19, hard_min=1e-06, hard_max=1e+20, frozen=True, units='cm^-3')
        self.H = XSParameter(name, 'H', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.He = XSParameter(name, 'He', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Li = XSParameter(name, 'Li', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Be = XSParameter(name, 'Be', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.B = XSParameter(name, 'B', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.F = XSParameter(name, 'F', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Na = XSParameter(name, 'Na', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Al = XSParameter(name, 'Al', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.P = XSParameter(name, 'P', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.K = XSParameter(name, 'K', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Sc = XSParameter(name, 'Sc', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Ti = XSParameter(name, 'Ti', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.V = XSParameter(name, 'V', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Mn = XSParameter(name, 'Mn', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Co = XSParameter(name, 'Co', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Cu = XSParameter(name, 'Cu', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Zn = XSParameter(name, 'Zn', 1.0, min=0.0, max=1000.0, hard_min=0.0, hard_max=1000.0, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, min=-0.999, max=10.0, hard_min=-0.999, hard_max=10.0, frozen=True)
        self.switch = XSParameter(name, 'switch', 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, min=0.0, max=1e+24, hard_min=0.0, hard_max=1e+24)
        XSAdditiveModel.__init__(self, name, (self.Tmax, self.beta, self.inv_slope, self.nH,
                                              self.H, self.He, self.Li, self.Be, self.B,
                                              self.C, self.N, self.O, self.F, self.Ne,
                                              self.Na, self.Mg, self.Al, self.Si, self.P,
                                              self.S, self.Cl, self.Ar, self.K, self.Ca,
                                              self.Sc, self.Ti, self.V, self.Cr, self.Mn,
                                              self.Fe, self.Co, self.Ni, self.Cu, self.Zn,
                                              self.redshift, self.switch, self.norm))


@version_at_least("12.12.0")
class XSvwdem(XSAdditiveModel):
    """The XSPEC vwdem model: plasma emission, multi-temperature with power-law distribution of emission measure.

    The model is described at [1]_.

    Attributes
    ----------
    Tmax
        The maximum temperature for power-law emission measure
        distribution.
    beta
        The ratio of minimum to maxmum temperature.
    inv_slope
        The inverse of the slope (labelled p in the XSPEC documentation).
    nH
        Fixed at 1 for most applications.
    He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        The abundance of the element in solar units.
    Redshift
        The redshift of the plasma.
    switch
        What model to use: 0 calculates with MEKAL, 1 interpolates
        with MEKAL, and 2 interpoates with APEC.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSvvwdem, XSwdem

    Notes
    -----
    This model is only available when used with XSPEC 12.12.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelWdem.html

    """
    __function__ = "C_vwDem"

    def __init__(self, name='vwdem'):
        self.Tmax = XSParameter(name, 'Tmax', 1.0, min=0.01, max=10.0, hard_min=0.01, hard_max=20.0, units='keV')
        self.beta = XSParameter(name, 'beta', 0.1, min=0.01, max=1.0, hard_min=0.01, hard_max=1.0)
        # can not use p for the name as it conflicts with P for the XSvvwdem model
        self.inv_slope = XSParameter(name, 'inv_slope', 0.25, min=-1.0, max=10.0, hard_min=-1.0, hard_max=10.0)
        self.nH = XSParameter(name, 'nH', 1.0, min=1e-05, max=1e+19, hard_min=1e-06, hard_max=1e+20, frozen=True, units='cm^-3')
        self.He = XSParameter(name, 'He', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.C = XSParameter(name, 'C', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.N = XSParameter(name, 'N', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.O = XSParameter(name, 'O', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Na = XSParameter(name, 'Na', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Al = XSParameter(name, 'Al', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Si = XSParameter(name, 'Si', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.S = XSParameter(name, 'S', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, min=-0.999, max=10.0, hard_min=-0.999, hard_max=10.0, frozen=True)
        self.switch = XSParameter(name, 'switch', 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, min=0.0, max=1e+24, hard_min=0.0, hard_max=1e+24)
        XSAdditiveModel.__init__(self, name, (self.Tmax, self.beta, self.inv_slope, self.nH,
                                              self.He, self.C, self.N, self.O, self.Ne,
                                              self.Na, self.Mg, self.Al, self.Si, self.S,
                                              self.Ar, self.Ca, self.Fe, self.Ni,
                                              self.redshift, self.switch, self.norm))


@version_at_least("12.12.0")
class XSwdem(XSAdditiveModel):
    """The XSPEC wdem model: plasma emission, multi-temperature with power-law distribution of emission measure.

    The model is described at [1]_.

    Attributes
    ----------
    Tmax
        The maximum temperature for power-law emission measure
        distribution.
    beta
        The ratio of minimum to maxmum temperature.
    inv_slope
        The inverse of the slope (labelled p in the XSPEC documentation).
    nH
        Fixed at 1 for most applications.
    abundanc
        The abundance relative to solar.
    Redshift
        The redshift of the plasma.
    switch
        What model to use: 0 calculates with MEKAL, 1 interpolates
        with MEKAL, and 2 interpoates with APEC.
    norm
        The normalization of the model: see [1]_ for an explanation
        of the units.

    See Also
    --------
    XSapec, XSmekal, XSvwdem, XSvvdem

    Notes
    -----
    This model is only available when used with XSPEC 12.12.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelWdem.html

    """
    __function__ = "C_wDem"

    def __init__(self, name='wdem'):
        self.Tmax = XSParameter(name, 'Tmax', 1.0, min=0.01, max=10.0, hard_min=0.01, hard_max=20.0, units='keV')
        self.beta = XSParameter(name, 'beta', 0.1, min=0.01, max=1.0, hard_min=0.01, hard_max=1.0)
        # can not use p for the name as it conflicts with P for the XSvvwdem model
        self.inv_slope = XSParameter(name, 'inv_slope', 0.25, min=-1.0, max=10.0, hard_min=-1.0, hard_max=10.0)
        self.nH = XSParameter(name, 'nH', 1.0, min=1e-05, max=1e+19, hard_min=1e-06, hard_max=1e+20, frozen=True, units='cm^-3')
        self.abundanc = XSParameter(name, 'abundanc', 1.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, min=-0.999, max=10.0, hard_min=-0.999, hard_max=10.0, frozen=True)
        self.switch = XSParameter(name, 'switch', 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, min=0.0, max=1e+24, hard_min=0.0, hard_max=1e+24)
        XSAdditiveModel.__init__(self, name, (self.Tmax, self.beta, self.inv_slope, self.nH,
                                              self.abundanc, self.redshift, self.switch,
                                              self.norm))


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
        self.LineE = XSParameter(name, 'LineE', 10.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.Sigma = XSParameter(name, 'Sigma', 1.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.Redshift = XSParameter(name, 'Redshift', 0., -0.999, 10.0, -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Sigma, self.Redshift, self.norm))


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
        self.kT = XSParameter(name, 'kT', 3.0, 1.e-2, 100., 1e-4, 200, units='keV')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.redshift, self.norm))


@version_at_least("12.10.1")
class XSzbknpower(XSAdditiveModel):
    """The XSPEC zbknpower model: broken power law.

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
        self.PhoIndx1 = XSParameter(name, 'PhoIndx1', 1., -2., 9., -3, 10)
        self.BreakE = XSParameter(name, 'BreakE', 5., 1.e-2, 1.e6, 0.0, 1e6, units='keV')
        self.PhoIndx2 = XSParameter(name, 'PhoIndx2', 2., -2., 9., -3, 10)
        self.Redshift = XSParameter(name, 'Redshift', 0, -0.999, 10, -0.999, 10,
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
        self.kT = XSParameter(name, 'kT', 7.0, 1.e-4, 100., 1e-4, 200, units='keV')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.PhoIndex = XSParameter(name, 'PhoIndex', 1., -2., 9., -3, 10)
        self.HighECut = XSParameter(name, 'HighECut', 15., 1., 500., 0.01, 500, units='keV')
        self.Redshift = XSParameter(name, 'Redshift', 0.0, -0.999, 10., -0.999, 10., frozen=True)
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
        self.LineE = XSParameter(name, 'LineE', 6.5, 0., 1.e6, 0.0, 1e6, units='keV')
        self.Sigma = XSParameter(name, 'Sigma', 0.1, 0., 10., 0.0, 20, units='keV')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Sigma, self.redshift, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


@version_at_least("12.11.0")
class XSzkerrbb(XSAdditiveModel):
    """The XSPEC zkerrbb model: multi-temperature blackbody model for thin accretion disk around a Kerr black hole.

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The fcol parameter was incorrectly labelled as hd: both names
       can be used to access this parameter.

       The default a, Mbh, fcol, and lflag parameter values have
       changed from 0, 1, 1.7, and 0 to 0.5, 1e7, 2.0, and 1 to match
       XSPEC.

    Attributes
    ----------
    eta
        The ratio of the disk power produced by a torque at the disk
        inner boundary to the disk power arising from accretion. See
        [1]_ for more details.
    a
        The specific angular momentum of the black hole in units of the
        black hole mass M (when G=c=1). It should be in the range [-1, 1).
    i
        The disk inclination angle, in degrees. A face-on disk has
        i=0. It must be less than or equal to 85 degrees.
    Mbh
        The mass of the black hole, in solar masses.
    Mdd
        The "effective" mass accretion rate in units of M_solar/year.
        See [1]_ for more details.
    z
        The redshift of the black hole
    fcol
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
    XSkerrbb

    Notes
    -----
    This model is only available when used with XSPEC 12.11.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKerrbb.html

    """

    __function__ = "C_zkerrbb"

    def __init__(self, name='zkerrbb'):
        self.eta = XSParameter(name, 'eta', 0., 0., 1.0, 0.0, 1, frozen=True)
        self.a = XSParameter(name, 'a', 0.5, -0.99, 0.9999, -0.99, 0.9999)
        self.i = XSParameter(name, 'i', 30., 0., 85., 0.0, 85, units='deg', frozen=True)
        self.Mbh = XSParameter(name, 'Mbh', 1e7, 3, 1e10, 3.0, 1e10, units='M_sun')
        self.Mdd = XSParameter(name, 'Mdd', 1., 1e-4, 1e4, 1e-5, 1e5, units='M0yr')
        self.z = XSParameter(name, 'z', 0.01, 0., 10., 0.0, 10, frozen=True)
        self.fcol = XSParameter(name, 'fcol', 2.0, -100, 100, -100, 100, frozen=True,
                                # Parameter was mis-labelled until 4.14.0
                                aliases=['hd'])
        self.rflag = XSParameter(name, 'rflag', 1., alwaysfrozen=True)
        self.lflag = XSParameter(name, 'lflag', 1., alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.eta, self.a, self.i, self.Mbh, self.Mdd, self.z, self.fcol, self.rflag, self.lflag, self.norm))


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
        self.alpha = XSParameter(name, 'alpha', 1.5, 0., 4., 0.0, 4)
        self.beta = XSParameter(name, 'beta', 0.2, -4., 4., -4, 4)
        self.pivotE = XSParameter(name, 'pivotE', 1.0, units='keV',
                                  alwaysfrozen=True)
        self.Redshift = XSParameter(name, 'Redshift', 0, -0.999, 10, -0.999, 10,
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
        self.PhoIndex = XSParameter(name, 'PhoIndex', 1., -2., 9., -3, 10)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.PhoIndex = XSParameter(name, 'PhoIndex', 2., 0., 4., 0.0, 4, frozen=True)
        self.nH = XSParameter(name, 'nH', 1., 0., 100., 0.0, 100, units='10^22 atoms / cm^2')
        self.Temp_abs = XSParameter(name, 'Temp_abs', 3.e4, 1.e4, 1.e6, 1e4, 1e6, units='K', frozen=True)
        self.xi = XSParameter(name, 'xi', 1., 0., 1.e3, 0.0, 5000)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1., 0., 1.e6, 0.0, 1e6, frozen=True)
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
        self.Tdays = XSParameter(name, 'Tdays', 850., 0., 10000., 0.0, 10000, units='days', frozen=True)
        self.norm = Parameter(name, 'norm', 0.00722, 0., 1., 0.0, hugeval, frozen=True)
        self.tauinf = XSParameter(name, 'tauinf', 0.582, 0., 1., 0.0, 1, frozen=True)
        self.tefold = XSParameter(name, 'tefold', 620., 1., 10000., 1.0, 10000, units='days', frozen=True)
        self.nC = XSParameter(name, 'nC', 10., 0., 50., 0.0, 50, frozen=True)
        self.nH = XSParameter(name, 'nH', 20., 1., 50., 1.0, 50, frozen=True)
        self.nO = XSParameter(name, 'nO', 2., 0., 50., 0.0, 50, frozen=True)
        self.nN = XSParameter(name, 'nN', 1., 0., 50., 0.0, 50, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.Tdays, self.norm, self.tauinf, self.tefold, self.nC, self.nH, self.nO, self.nN))


class XSconstant(XSMultiplicativeModel):
    """The XSPEC constant model: energy-independent factor.

    The model is described at [1]_.

    Attributes
    ----------
    factor
        The value of the model.

    See Also
    --------
    XSlogconst, XSlog10con

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelConstant.html

    """

    __function__ = "xscnst"

    def __init__(self, name='constant'):
        self.factor = XSParameter(name, 'factor', 1., 0.0, 1.e10, 0.0, 1e10)
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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
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
        self.Depth0 = XSParameter(name, 'Depth0', 2.0, 0., 100., 0.0, 100)
        self.E0 = XSParameter(name, 'E0', 30.0, 1.0, 100., 1.0, 100, units='keV')
        self.Width0 = XSParameter(name, 'Width0', 10.0, 1.0, 100., 1.0, 100, units='keV', frozen=True)
        self.Depth2 = XSParameter(name, 'Depth2', 0.0, 0., 100., 0.0, 100, frozen=True)
        self.Width2 = XSParameter(name, 'Width2', 20.0, 1.0, 100., 1.0, 100, units='keV', frozen=True)
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
        self.Frac = XSParameter(name, 'Frac', 0.066, 0., 1., 0.0, 1, frozen=True)
        self.Halosz = XSParameter(name, 'Halosz', 2., 0., 1.e5, 0.0, 1e5, frozen=True)
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
        self.edgeE = XSParameter(name, 'edgeE', 7.0, 0., 100., 0.0, 100, units='keV')
        self.MaxTau = XSParameter(name, 'MaxTau', 1., 0., 5., 0.0, 10)
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
        self.LowECut = XSParameter(name, 'LowECut', 2., 0., 100., 0.0, 200, units='keV')
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
        self.Ampl = XSParameter(name, 'Ampl', 1., 0., 1.e5, 0.0, 1e6)
        self.Factor = XSParameter(name, 'Factor', 1., 0., 1.e5, 0.0, 1e6)
        self.StartE = XSParameter(name, 'StartE', 0.5, 0., 1.e5, 0.0, 1e6, units='keV', frozen=True)
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
        self.LineE = XSParameter(name, 'LineE', 1.0, 0., 1.e6, 0.0, 1e6, units='keV')
        self.Sigma = XSParameter(name, 'Sigma', 0.01, 0., 10., 0.0, 20, units='keV')
        self.Strength = XSParameter(name, 'Strength', 1.0, 0., 1.e6, 0.0, 1e6, aliases=["Tau"])

        XSMultiplicativeModel.__init__(self, name, (self.LineE, self.Sigma, self.Strength))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


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
        self.nHeI = XSParameter(name, 'nHeI', 1.e-5, 0.0, 1.e6, 0.0, 1.0e6, units='10^22 atoms / cm^2')
        self.b = XSParameter(name, 'b', 10.0, 1.0, 1.0e5, 1.0, 1.0e6, units='km/s')
        self.z = XSParameter(name, 'z', 0.0, -1.0e-3, 1.0e5, -1.0e-3, 1.0e5, aliases=["redshift"])

        XSMultiplicativeModel.__init__(self, name, (self.nHei, self.b, self.z))


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
        self.cutoffE = XSParameter(name, 'cutoffE', 10., 1.e-2, 1.e6, 1e-4, 1e6, units='keV')
        self.foldE = XSParameter(name, 'foldE', 15., 1.e-2, 1.e6, 1e-4, 1e6, units='keV')
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
        self.thetamin = XSParameter(name, 'thetamin', 0., 0.0, 90., 0.0, 90, frozen=True)
        self.thetamax = XSParameter(name, 'thetamax', 90., 0.0, 90., 0.0, 90, frozen=True)
        self.thetaobs = XSParameter(name, 'thetaobs', 60., 0.0, 90., 0.0, 90)
        self.Feabun = XSParameter(name, 'Feabun', 1., 0.0, 100., 0.0, 200, frozen=True)
        self.FeKedge = XSParameter(name, 'FeKedge', 7.11, 7.0, 10., 7.0, 10, units='keV', frozen=True)
        self.Escfrac = XSParameter(name, 'Escfrac', 1.0, 0.0, 500., 0.0, 1000)
        self.covfac = XSParameter(name, 'covfac', 1.0, 0.0, 500., 0.0, 1000)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.thetamin, self.thetamax, self.thetaobs, self.Feabun, self.FeKedge, self.Escfrac, self.covfac, self.redshift))


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

    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelIsmabs.html

    """

    __function__ = "ismabs"

    def __init__(self, name='ismabs'):

        self.H = XSParameter(name, 'H', 0.1, 0.0, 1e5, 0, 1e6, units='10^22')
        self.HeII = XSParameter(name, 'HeII', 0.0, 0.0, 1e5, 0, 1e6,
                                units='10^16', frozen=True)
        self.CI = XSParameter(name, 'CI', 33.1, 0.0, 1e5, 0, 1e6,
                              units='10^16')
        self.CII = XSParameter(name, 'CII', 0.0, 0.0, 1e5, 0, 1e6,
                               units='10^16', frozen=True)
        self.CIII = XSParameter(name, 'CIII', 0.0, 0.0, 1e5, 0, 1e6,
                                units='10^16', frozen=True)
        self.NI = XSParameter(name, 'NI', 8.32, 0.0, 1e5, 0, 1e6,
                              units='10^16')
        self.NII = XSParameter(name, 'NII', 0.0, 0.0, 1e5, 0, 1e6,
                               units='10^16', frozen=True)
        self.NIII = XSParameter(name, 'NIII', 0.0, 0.0, 1e5, 0, 1e6,
                                units='10^16', frozen=True)
        self.OI = XSParameter(name, 'OI', 67.6, 0.0, 1e5, 0, 1e6,
                              units='10^16')
        self.OII = XSParameter(name, 'OII', 0.0, 0.0, 1e5, 0, 1e6,
                               units='10^16', frozen=True)
        self.OIII = XSParameter(name, 'OIII', 0.0, 0.0, 1e5, 0, 1e6,
                                units='10^16', frozen=True)
        self.NeI = XSParameter(name, 'NeI', 12.0, 0.0, 1e5, 0, 1e6,
                               units='10^16')
        self.NeII = XSParameter(name, 'NeII', 0.0, 0.0, 1e5, 0, 1e6,
                                units='10^16', frozen=True)
        self.NeIII = XSParameter(name, 'NeIII', 0.0, 0.0, 1e5, 0, 1e6,
                                 units='10^16', frozen=True)
        self.MgI = XSParameter(name, 'MgI', 3.8, 0.0, 1e5, 0, 1e6,
                               units='10^16')
        self.MgII = XSParameter(name, 'MgII', 0.0, 0.0, 1e5, 0, 1e6,
                                units='10^16', frozen=True)
        self.MgIII = XSParameter(name, 'MgIII', 0.0, 0.0, 1e5, 0, 1e6,
                                 units='10^16', frozen=True)
        # SiI and SI conflict, so add in underscores to differentiate.
        #
        self.Si_I = XSParameter(name, 'Si_I', 3.35, 0.0, 1e5, 0, 1e6,
                                units='10^16')
        self.Si_II = XSParameter(name, 'Si_II', 0.0, 0.0, 1e5, 0, 1e6,
                                 units='10^16', frozen=True)
        self.Si_III = XSParameter(name, 'Si_III', 0.0, 0.0, 1e5, 0, 1e6,
                                  units='10^16', frozen=True)
        self.S_I = XSParameter(name, 'S_I', 2.14, 0.0, 1e5, 0, 1e6,
                               units='10^16')
        self.S_II = XSParameter(name, 'S_II', 0.0, 0.0, 1e5, 0, 1e6,
                                units='10^16', frozen=True)
        self.S_III = XSParameter(name, 'S_III', 0.0, 0.0, 1e5, 0, 1e6,
                                 units='10^16', frozen=True)
        self.ArI = XSParameter(name, 'ArI', 0.25, 0.0, 1e5, 0, 1e6,
                               units='10^16')
        self.ArII = XSParameter(name, 'ArII', 0.0, 0.0, 1e5, 0, 1e6,
                                units='10^16', frozen=True)
        self.ArIII = XSParameter(name, 'ArIII', 0.0, 0.0, 1e5, 0, 1e6,
                                 units='10^16', frozen=True)
        self.CaI = XSParameter(name, 'CaI', 0.22, 0.0, 1e5, 0, 1e6,
                               units='10^16')
        self.CaII = XSParameter(name, 'CaII', 0.0, 0.0, 1e5, 0, 1e6,
                                units='10^16', frozen=True)
        self.CaIII = XSParameter(name, 'CaIII', 0.0, 0.0, 1e5, 0, 1e6,
                                 units='10^16', frozen=True)
        self.Fe = XSParameter(name, 'Fe', 3.16, 0.0, 1e5, 0, 1e6, units='10^16')
        self.redshift = XSParameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
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


@version_at_least("12.11.0")
class XSismdust(XSMultiplicativeModel):
    """The XSPEC ismdust model: Extinction due to a power-law distribution of dust grains.

    The model is described at [1]_.

    Attributes
    ----------
    msil
        The dust mass column for silicate (in units of 10^-4).
    mgra
        The dust mass column for graphite (in units of 10^-4).
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSolivineabs

    Notes
    -----
    This model is only available when used with XSPEC 12.11.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelIsmdust.html

    """

    __function__ = "ismdust"

    def __init__(self, name='ismdust'):
        self.msil = XSParameter(name, 'msil', 1.0, 0.0, 1e4, 0, 1e5, units='10^-4')
        self.mgra = XSParameter(name, 'mgra', 1.0, 0.0, 1e4, 0, 1e5, units='10^-4')
        self.redshift = XSParameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
                                    frozen=True)
        XSMultiplicativeModel.__init__(self, name,
                                       (self.msil, self.mgra,
                                        self.redshift))


@version_at_least("12.11.0")
class XSlogconst(XSMultiplicativeModel):
    """The XSPEC logconst model: Constant in log units.

    The model is described at [1]_.

    Attributes
    ----------
    logfact
        The constant factor in natural log.

    See Also
    --------
    XSconstant, XSlog10con

    Notes
    -----
    This model is only available when used with XSPEC 12.11.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLogconst.html

    """

    __function__ = "C_logconst"

    def __init__(self, name='logconst'):
        self.logfact = XSParameter(name, 'logfact', 0.0, -20.0, 20, -20, 20)
        XSMultiplicativeModel.__init__(self, name, (self.logfact, ))


@version_at_least("12.11.0")
class XSlog10con(XSMultiplicativeModel):
    """The XSPEC log10con model: Constant in base 10 log units.

    The model is described at [1]_.

    Attributes
    ----------
    log10fac
        The constant factor in base 10 log.

    See Also
    --------
    XSconstant, XSlogconst

    Notes
    -----
    This model is only available when used with XSPEC 12.11.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLog10con.html

    """

    __function__ = "C_log10con"

    def __init__(self, name='log10con'):
        self.log10fac = XSParameter(name, 'log10fac', 0.0, -20.0, 20, -20, 20)
        XSMultiplicativeModel.__init__(self, name, (self.log10fac, ))


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
        self.n = XSParameter(name, 'n', 1.e-5, 0.0, 1.0e6, 0.0, 1.0e6, units='10^22 atoms / cm^2', aliases=["nHeI"])
        self.b = XSParameter(name, 'b', 10.0, 1.0, 1.0e5, 1.0, 1.0e6, units='km/s')
        self.z = XSParameter(name, 'z', 0.0, -1.0e-3, 1.0e5, -1.0e-3, 1.0e5, aliases=["redshift"])
        self.ZA = XSParameter(name, 'ZA', 1.0, 1.0, 2.0, 1.0, 2.0)

        XSMultiplicativeModel.__init__(self, name, (self.n, self.b, self.z, self.ZA))


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
        self.LineE = XSParameter(name, 'LineE', 3.5, 0., 20., 0.0, 20, units='keV')
        self.Width = XSParameter(name, 'Width', 1., 0., 20., 0.0, 20, units='keV')
        self.CvrFract = XSParameter(name, 'CvrFract', 1., 0., 1., 0.0, 1, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.LineE, self.Width, self.CvrFract))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


@version_at_least("12.11.0")
class XSolivineabs(XSMultiplicativeModel):
    """The XSPEC olivineabs model: Absorption due to olivine.

    The model is described at [1]_.

    Attributes
    ----------
    moliv
        The dust mass column for olivine grains (in units of 10^-4).
    redshift
        The redshift of the absorber.

    See Also
    --------
    XSismdust

    Notes
    -----
    This model is only available when used with XSPEC 12.11.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelOlivineabs.html

    """

    __function__ = "olivineabs"

    def __init__(self, name='olivineabs'):
        self.moliv = XSParameter(name, 'moliv', 1.0, 0.0, 1e4, 0, 1e5, units='10^-4')
        self.redshift = XSParameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
                                    frozen=True)
        XSMultiplicativeModel.__init__(self, name,
                                       (self.moliv, self.redshift))


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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.CvrFract = XSParameter(name, 'CvrFract', 0.5, 0.05, 0.95, 0.0, 1)
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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
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
        self.index = XSParameter(name, 'index', 2.0, 0.0, 5., 0.0, 5)
        self.coef = XSParameter(name, 'coef', 1.0, 0.0, 100., 0.0, 100)
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
    XSpcfabs, XSwabs, XSzxipab

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPwab.html

    """

    __function__ = "C_xspwab"

    def __init__(self, name='pwab'):
        self.nHmin = XSParameter(name, 'nHmin', 1., 1.e-7, 1.e5, 1e-7, 1e6, units='10^22 atoms / cm^2')
        self.nHmax = XSParameter(name, 'nHmax', 2., 1.e-7, 1.e5, 1e-7, 1e6, units='10^22 atoms / cm^2')
        self.beta = XSParameter(name, 'beta', 1.0, -10., 10, -10, 20, frozen=True)
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
        The value of E(B-V) for the line of sight to the source.

    See Also
    --------
    XSzredden

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRedden.html

    """

    __function__ = "xscred"

    def __init__(self, name='redden'):
        self.E_BmV = XSParameter(name, 'E_BmV', 0.05, 0., 10., 0.0, 10, aliases=["EBV"])

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
        self.edgeE = XSParameter(name, 'edgeE', 7.0, 0.1, 100., 0.1, 100, units='keV')
        self.MaxTau = XSParameter(name, 'MaxTau', 1., 0., 5., 0.0, 10)
        self.index = XSParameter(name, 'index', -2.67, -10., 10., -10, 10, frozen=True)
        self.width = XSParameter(name, 'width', 10., 0.01, 100., 0.01, 100)
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
        self.Ecut = XSParameter(name, 'Ecut', 10.0, 0.0, 1e6, 0.0, 1e6, units='keV')
        self.alpha = XSParameter(name, 'alpha', 1.0, -5.0, 5.0, -5, 5)
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
        self.Estart = XSParameter(name, 'Estart', 0.1, 0., 100., 0.0, 100, units='keV')
        self.Ystart = XSParameter(name, 'Ystart', 1., -1.e6, 1.e6, -1e6, 1e6)
        self.Yend = XSParameter(name, 'Yend', 1., -1.e6, 1.e6, -1e6, 1e6)
        self.YPstart = XSParameter(name, 'YPstart', 0., -1.e6, 1.e6, -1e6, 1e6)
        self.YPend = XSParameter(name, 'YPend', 0., -1.e6, 1.e6, -1e6, 1e6)
        self.Eend = XSParameter(name, 'Eend', 15., 0., 100., 0.0, 100, units='keV')
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
        self.clumps = XSParameter(name, 'clumps', 0.0, 0., 10., 0.0, 10)
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
        self.column = XSParameter(name, 'column', 6., 3., 50., 3.0, 50)
        self.log_xi = XSParameter(name, 'log_xi', 2.5, 2.1, 4.1, 2.1, 4.1, aliases=["logxi"])
        self.sigma = XSParameter(name, 'sigma', 0.1, 0., 0.5, 0.0, 0.5)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)

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
        self.nH = XSParameter(name, 'nH', 1., 0., 1E5, 0.0, 1e6, units='10^22 atoms / cm^2')
        XSMultiplicativeModel.__init__(self, name, (self.nH,))


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

    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbfeo"

    def __init__(self, name='tbfeo'):
        self.nH = XSParameter(name, 'nH', 1., 0., 1.e5, 0.0, 1.0e6, units='10^22')
        self.O = XSParameter(name, 'O', 1., 0.0, 5.0, -1.0e38, 1.0e38)
        self.Fe = XSParameter(name, 'Fe', 1., 0.0, 5.0, -1.0e38, 1.0e38)
        self.redshift = XSParameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
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

    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbgas"

    def __init__(self, name='tbgas'):
        self.nH = XSParameter(name, 'nH', 1., 0., 1.e5, 0.0, 1.0e6, units='10^22')
        self.redshift = XSParameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
                                    frozen=True)
        XSMultiplicativeModel.__init__(self, name,
                                       (self.nH, self.redshift))


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
        self.nH = XSParameter(name, 'nH', 1., 0., 1E5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.h2 = XSParameter(name, 'h2', 0.2, 0., 1., 0.0, 1, frozen=True)
        self.rho = XSParameter(name, 'rho', 1., 0., 5., 0.0, 5, units='g/cm^3', frozen=True)
        self.amin = XSParameter(name, 'amin', 0.025, 0., 0.25, 0.0, 0.25, units='mum', frozen=True)
        self.amax = XSParameter(name, 'amax', 0.25, 0., 1., 0.0, 1, units='mum', frozen=True)
        self.PL = XSParameter(name, 'PL', 3.5, 0., 5., 0.0, 5, frozen=True)
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
        self.nH = XSParameter(name, 'nH', 1., 0., 1E5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.He = XSParameter(name, 'He', 1., 0., 1., 0.0, 1, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1., 0.0, 1, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1., 0.0, 1, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1., 0.0, 1, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1., 0.0, 1, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1., 0.0, 1, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1., 0.0, 1, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1., 0.0, 1, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1., 0.0, 1, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1., 0.0, 1, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1., 0.0, 1, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1., 0.0, 1, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1., 0.0, 1, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1., 0.0, 1, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1., 0.0, 1, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1., 0.0, 1, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1., 0.0, 1, frozen=True)
        self.H2 = XSParameter(name, 'H2', 0.2, 0., 1., 0.0, 1, frozen=True)
        self.rho = XSParameter(name, 'rho', 1., 0., 5., 0.0, 5, units='g/cm^3', frozen=True)
        self.amin = XSParameter(name, 'amin', 0.025, 0., 0.25, 0.0, 0.25, units='mum', frozen=True)
        self.amax = XSParameter(name, 'amax', 0.25, 0., 1., 0.0, 1, units='mum', frozen=True)
        self.PL = XSParameter(name, 'PL', 3.5, 0., 5., 0.0, 5, frozen=True)
        self.H_dep = XSParameter(name, 'H_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.He_dep = XSParameter(name, 'He_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.C_dep = XSParameter(name, 'C_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.N_dep = XSParameter(name, 'N_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.O_dep = XSParameter(name, 'O_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Ne_dep = XSParameter(name, 'Ne_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Na_dep = XSParameter(name, 'Na_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Mg_dep = XSParameter(name, 'Mg_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Al_dep = XSParameter(name, 'Al_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Si_dep = XSParameter(name, 'Si_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.S_dep = XSParameter(name, 'S_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Cl_dep = XSParameter(name, 'Cl_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Ar_dep = XSParameter(name, 'Ar_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Ca_dep = XSParameter(name, 'Ca_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Cr_dep = XSParameter(name, 'Cr_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Fe_dep = XSParameter(name, 'Fe_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Co_dep = XSParameter(name, 'Co_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.Ni_dep = XSParameter(name, 'Ni_dep', 1., 0., 1., 0.0, 1, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Cl, self.Ar, self.Ca, self.Cr, self.Fe, self.Co, self.Ni, self.H2, self.rho, self.amin, self.amax, self.PL, self.H_dep, self.He_dep, self.C_dep, self.N_dep, self.O_dep, self.Ne_dep, self.Na_dep, self.Mg_dep, self.Al_dep, self.Si_dep, self.S_dep, self.Cl_dep, self.Ar_dep, self.Ca_dep, self.Cr_dep, self.Fe_dep, self.Co_dep, self.Ni_dep, self.redshift))


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

    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbpcf"

    def __init__(self, name='tbpcf'):
        self.nH = XSParameter(name, 'nH', 1., 0., 1.e5, 0.0, 1.0e6, units='10^22')
        self.pcf = XSParameter(name, 'pcf', 0.5, 0, 1.0, 0, 1.0)
        self.redshift = XSParameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
                                    frozen=True)
        XSMultiplicativeModel.__init__(self, name,
                                       (self.nH, self.pcf, self.redshift))


@version_at_least("12.9.1")
class XSTBrel(XSMultiplicativeModel):
    """The XSPEC TBrel model: ISM grain absorption.

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The default nH parameter value has changed from 1.0 to 0.0 to
       match XSPEC.

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

    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelTbabs.html

    """

    __function__ = "C_tbrel"

    def __init__(self, name='tbrel'):
        self.nH = XSParameter(name, 'nH', 0.0, -1e5, 1e5, -1e6, 1.0e6, units='10^22 atoms / cm^2')
        self.He = XSParameter(name, 'He', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 5., 0.0, 1e38, frozen=True)
        self.H2 = XSParameter(name, 'H2', 0.2, 0., 1., 0.0, 1.0, frozen=True)
        self.rho = XSParameter(name, 'rho', 1., 0., 5., 0.0, 5.0, units='g/cm^3',
                               frozen=True)
        self.amin = XSParameter(name, 'amin', 0.025, 0., 0.25, 0.0, 0.25,
                                units='mum', frozen=True)
        self.amax = XSParameter(name, 'amax', 0.25, 0., 1., 0.0, 1.0,
                                units='mum', frozen=True)
        self.PL = XSParameter(name, 'PL', 3.5, 0., 5., 0.0, 5.0, frozen=True)
        self.H_dep = XSParameter(name, 'H_dep', 1., 0., 1., 0.0, 1.0,
                                 frozen=True)
        self.He_dep = XSParameter(name, 'He_dep', 1., 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.C_dep = XSParameter(name, 'C_dep', 0.5, 0., 1., 0.0, 1.0,
                                 frozen=True)
        self.N_dep = XSParameter(name, 'N_dep', 1., 0., 1., 0.0, 1.0,
                                 frozen=True)
        self.O_dep = XSParameter(name, 'O_dep', 0.6, 0., 1., 0.0, 1.0,
                                 frozen=True)
        self.Ne_dep = XSParameter(name, 'Ne_dep', 1., 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.Na_dep = XSParameter(name, 'Na_dep', 0.25, 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.Mg_dep = XSParameter(name, 'Mg_dep', 0.2, 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.Al_dep = XSParameter(name, 'Al_dep', 0.02, 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.Si_dep = XSParameter(name, 'Si_dep', 0.1, 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.S_dep = XSParameter(name, 'S_dep', 0.6, 0., 1., 0.0, 1.0,
                                 frozen=True)
        self.Cl_dep = XSParameter(name, 'Cl_dep', 0.5, 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.Ar_dep = XSParameter(name, 'Ar_dep', 1., 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.Ca_dep = XSParameter(name, 'Ca_dep', 0.003, 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.Cr_dep = XSParameter(name, 'Cr_dep', 0.03, 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.Fe_dep = XSParameter(name, 'Fe_dep', 0.3, 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.Co_dep = XSParameter(name, 'Co_dep', 0.05, 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.Ni_dep = XSParameter(name, 'Ni_dep', 0.04, 0., 1., 0.0, 1.0,
                                  frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., 0.0, 10., -1.0, 10.0,
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
        The value of E(B-V) for the line of sight to the source.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelUvred.html

    """

    __function__ = "xsred"

    def __init__(self, name='uvred'):
        self.E_BmV = XSParameter(name, 'E_BmV', 0.05, 0., 10., 0.0, 10, aliases=["EBV"])

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
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 10000, units='sH22', frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 10000, units='sHe22', frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 10000, units='sC22', frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 10000, units='sN22', frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 10000, units='sO22', frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 10000, units='sNe22', frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 10000, units='sNa22', frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 10000, units='sMg22', frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 10000, units='sAl22', frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 10000, units='sSi22', frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 10000, units='sS22', frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 10000, units='sCl22', frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 10000, units='sAr22', frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 10000, units='sCa22', frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 10000, units='sCr22', frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 10000, units='sFe22', frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 10000, units='sCo22', frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 10000, units='sNi22', frozen=True)
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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
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
        self.nH = XSParameter(name, 'nH', 1., 0., 10., 0.0, 20, units='10^22 atoms / cm^2')
        self.WindowE = XSParameter(name, 'WindowE', 1., .05, 20., 0.03, 20, units='keV')
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
        The outer radius of the disk (in Schwarzschild radii).
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
        self.height = XSParameter(name, 'height', 5., 0.0, 1.e2, 0.0, 1e2, units='r_s')
        self.lxovrld = XSParameter(name, 'lxovrld', 0.3, 0.02, 100, 0.02, 100, aliases=["lxld"])
        self.rate = XSParameter(name, 'rate', 0.05, 1.e-3, 1., 1e-3, 1)
        self.cosAng = XSParameter(name, 'cosAng', 0.9, 0., 1., 0.0, 1)
        self.inner = XSParameter(name, 'inner', 3., 2., 1.e3, 2.0, 1000, units='r_s')
        self.outer = XSParameter(name, 'outer', 100., 2.1, 1.e5, 2.1, 1e5, units='r_s')
        self.index = XSParameter(name, 'index', 2.0, 1.6, 2.2, 1.6, 2.2)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        self.Feabun = XSParameter(name, 'Feabun', 1., 0., 5., 0.0, 5, frozen=True)
        self.E_cut = XSParameter(name, 'E_cut', 150., 20., 300., 20.0, 300, units='keV')
        self.Ref_type = XSParameter(name, 'Ref_type', 1., 1., 3., 1.0, 3, frozen=True)
        self.Rel_smear = XSParameter(name, 'Rel_smear', 4., 1., 4., 1.0, 4, frozen=True)
        self.Geometry = XSParameter(name, 'Geometry', 1., 1., 4., 1.0, 4, frozen=True)

        XSMultiplicativeModel.__init__(self, name, (self.height, self.lxovrld, self.rate, self.cosAng, self.inner, self.outer, self.index, self.redshift, self.Feabun, self.E_cut, self.Ref_type, self.Rel_smear, self.Geometry))


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

    Notes
    -----
    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelXscat.html

    """

    __function__ = "C_xscatmodel"

    def __init__(self, name='xscat'):
        self.NH = XSParameter(name, 'NH', 1., 0., 1000.0, 0.0, 1000.0, units='10^22')
        self.Xpos = XSParameter(name, 'Xpos', 0.5, 0, 0.99, 0, 0.999)
        self.Rext = XSParameter(name, 'Rext', 10.0, 0, 235.0, 0, 240.0, 'arcsec',
                                frozen=True)
        self.DustModel = XSParameter(name, 'DustModel', 1,
                                     alwaysfrozen=True)
        XSMultiplicativeModel.__init__(self, name,
                                       (self.NH, self.Xpos, self.Rext,
                                        self.DustModel))


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
    nHeII
        The He II column density, in 10^22 atoms/cm^2.
    z
        The redshift of the absorber.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZbabs.html

    """

    __function__ = "xszbabs"

    def __init__(self, name='zbabs'):
        self.nH = XSParameter(name, 'nH', 1.e-4, 0.0, 1.0e5, 0.0, 1.0e6, units='10^22 atoms / cm^2')
        self.nHeI = XSParameter(name, 'nHeI', 1.e-5, 0.0, 1.0e5, 0.0, 1.0e6, units='10^22 atoms / cm^2')
        self.nHeII = XSParameter(name, 'nHeII', 1.e-6, 0.0, 1.0e5, 0.0, 1.0e6, units='10^22 atoms / cm^2')
        self.z = XSParameter(name, 'z', 0.0, 0.0, 1.0e5, 0.0, 1.0e6, aliases=["redshift"])

        XSMultiplicativeModel.__init__(self, name, (self.nH, self.nHeI, self.nHeII, self.z))


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
        self.method = XSParameter(name, 'method', 1, 1, 3, 1, 3, alwaysfrozen=True)
        self.E_BmV = XSParameter(name, 'E_BmV', 0.1, 0.0, 100., 0.0, 100, aliases=["EBV"])
        self.Rv = XSParameter(name, 'Rv', 3.1, 0.0, 10., 0.0, 10, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0.0, 0.0, 20., 0.0, 20, frozen=True)

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
        self.edgeE = XSParameter(name, 'edgeE', 7.0, 0., 100., 0.0, 100, units='keV')
        self.MaxTau = XSParameter(name, 'MaxTau', 1., 0., 5., 0.0, 10)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.cutoffE = XSParameter(name, 'cutoffE', 10., 1.e-2, 100., 1e-4, 200, units='keV')
        self.foldE = XSParameter(name, 'foldE', 15., 1.e-2, 100., 1e-4, 200, units='keV')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.cutoffE, self.foldE, self.redshift))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.cutoffE, **kwargs)


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
        self.redshift = XSParameter(name, 'redshift', 0.0, alwaysfrozen=True)
        self.model = XSParameter(name, 'model', 0, alwaysfrozen=True)
        self.lyman_limit = XSParameter(name, 'lyman_limit', 1, alwaysfrozen=True)

        XSMultiplicativeModel.__init__(self, name, (self.redshift, self.model, self.lyman_limit))


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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.CvrFract = XSParameter(name, 'CvrFract', 0.5, 0.05, 0.95, 0.0, 1)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.redshift))


@version_at_least("12.12.0")
class XSzxipab(XSMultiplicativeModel):
    """The XSPEC zxipab model: power-law distribution of ionized absorbers.

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
    XSpwab

    Notes
    -----
    This model is only available when used with XSPEC 12.12.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZxipab.html

    """

    __function__ = "zxipab"

    def __init__(self, name='zxipab'):
        self.nHmin = XSParameter(name, 'nHmin', 0.01, min=1e-07, max=1000.0, hard_min=1e-07, hard_max=1000000.0, units='10^22')
        self.nHmax = XSParameter(name, 'nHmax', 10.0, min=1e-07, max=1000.0, hard_min=1e-07, hard_max=1000000.0, units='10^22')
        self.beta = XSParameter(name, 'beta', 0.0, min=-10.0, max=10.0, hard_min=-10.0, hard_max=10.0)
        self.log_xi = XSParameter(name, 'log_xi', 3.0, min=-3.0, max=6.0, hard_min=-3.0, hard_max=6.0)
        self.redshift = XSParameter(name, 'redshift', 0.0, min=0.0, max=10.0, hard_min=0.0, hard_max=10.0, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nHmin, self.nHmax, self.beta,
                                                    self.log_xi, self.redshift))


class XSzxipcf(XSMultiplicativeModel):
    """The XSPEC zxipcf model: partial covering absorption by partially ionized material.

    The model is described at [1]_.

    .. note:: Deprecated in Sherpa 4.10.0

       The ``logxi`` parameter has been renamed ``log_xi`` to match the
       XSPEC definition. The name ``logxi`` can still be used to access
       the parameter, but this name will be removed in a future release.

    Attributes
    ----------
    Nh
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
        self.Nh = XSParameter(name, 'Nh', 10, 0.05, 500, 0.05, 500, units='10^22 atoms / cm^2')
        self.log_xi = XSParameter(name, 'log_xi', 3, -3, 6, -3, 6, aliases=["logxi"])
        self.CvrFract = XSParameter(name, 'CvrFract', 0.5, 0., 1., 0.0, 1)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)

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
        The value of E(B-V) for the line of sight to the source.
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
        self.E_BmV = XSParameter(name, 'E_BmV', 0.05, 0., 10., 0.0, 10, aliases=["EBV"])
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)

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
        The value of E(B-V) for the line of sight to the source.
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
        self.E_BmV = XSParameter(name, 'E_BmV', 0.1, 0.0, 100., 0.0, 100, aliases=["EBV"])
        self.ExtIndex = XSParameter(name, 'ExtIndex', 1.0, -10.0, 10., -10, 10)
        self.Rv = XSParameter(name, 'Rv', 3.1, 0.0, 10., 0.0, 10, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0.0, 0.0, 20., 0.0, 20, units='z', frozen=True)

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
        self.nH = XSParameter(name, 'nH', 1., 0., 1E5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.H = XSParameter(name, 'H', 1., 0., 1000., 0.0, 10000, units='sH22', frozen=True)
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 10000, units='sHe22', frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 10000, units='sC22', frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 10000, units='sN22', frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 10000, units='sO22', frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 10000, units='sNe22', frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 10000, units='sNa22', frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 10000, units='sMg22', frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 10000, units='sAl22', frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 10000, units='sSi22', frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 10000, units='sS22', frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 10000, units='sCl22', frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 10000, units='sAr22', frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 10000, units='sCa22', frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 10000, units='sCr22', frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 10000, units='sFe22', frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 10000, units='sCo22', frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 10000, units='sNi22', frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.metals = XSParameter(name, 'metals', 1., 0.0, 100., 0.0, 100)
        self.FEabun = XSParameter(name, 'FEabun', 1., 0.0, 100., 0.0, 100)
        self.FEKedge = XSParameter(name, 'FEKedge', 7.11, 7.0, 9.5, 7.0, 9.5, units='keV')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.He = XSParameter(name, 'He', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.C = XSParameter(name, 'C', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.N = XSParameter(name, 'N', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.O = XSParameter(name, 'O', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ne = XSParameter(name, 'Ne', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Na = XSParameter(name, 'Na', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Mg = XSParameter(name, 'Mg', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Al = XSParameter(name, 'Al', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Si = XSParameter(name, 'Si', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.S = XSParameter(name, 'S', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cl = XSParameter(name, 'Cl', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ar = XSParameter(name, 'Ar', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ca = XSParameter(name, 'Ca', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Cr = XSParameter(name, 'Cr', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Fe = XSParameter(name, 'Fe', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Co = XSParameter(name, 'Co', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.Ni = XSParameter(name, 'Ni', 1., 0., 1000., 0.0, 1000, frozen=True)
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.nH = XSParameter(name, 'nH', 1., 0.0, 1.e5, 0.0, 1e6, units='10^22 atoms / cm^2')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
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
        self.nH = XSParameter(name, 'nH', 1., 0., 10., 0.0, 20, units='10^22 atoms / cm^2')
        self.WindowE = XSParameter(name, 'WindowE', 1., .05, 20., 0.03, 20, units='keV')
        self.redshift = XSParameter(name, 'redshift', 0., -0.999, 10., -0.999, 10, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.WindowE, self.redshift))


class XScflux(XSConvolutionKernel):
    """The XSPEC cflux convolution model: calculate flux

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Emin
        Minimum energy over which the flux is calculated.
    Emax
        Maximum energy over which the flux is calculated.
    lg10Flux
        log (base 10) of the flux in erg/cm^2/s

    See Also
    --------
    XSclumin, XScpflux

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    See [1]_ for the meaning and restrictions, in particular the
    necessity of freezing the amplitude, or normalization, of the
    emission component (or components) at 1.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCflux.html

    Examples
    --------

    With the following definitions:

    >>> cflux = XScflux()
    >>> absmdl = XSphabs()
    >>> plmdl = XSpowerlaw()
    >>> gmdl = XSgaussian()
    >>> srcmdl = plmdl + gmdl

    then the model can be applied in a number of ways, such as:

    >>> mdl1 = cflux(absmdl * srcmdl)
    >>> mdl2 = absmdl * cflux(srcmdl)
    >>> mdl3 = absmdl * (plmdl + cflux(gmdl))

    """

    _calc = _xspec.C_cflux

    def __init__(self, name='xscflux'):
        self.Emin = XSParameter(name, 'Emin', 0.5, min=0.0, max=1e6,
                                hard_min=0.0, hard_max=1e6, frozen=True,
                                units='keV')
        self.Emax = XSParameter(name, 'Emax', 10.0, min=0.0, max=1e6,
                                hard_min=0.0, hard_max=1e6, frozen=True,
                                units='keV')
        self.lg10Flux = XSParameter(name, 'lg10Flux', -12.0, min=-100.0,
                                    max=100.0, hard_min=-100.0, hard_max=100.0,
                                    frozen=False, units='cgs')
        XSConvolutionKernel.__init__(self, name, (self.Emin,
                                                  self.Emax,
                                                  self.lg10Flux
                                                  ))


@version_at_least("12.9.1")
class XSclumin(XSConvolutionKernel):
    """The XSPEC clumin convolution model: calculate luminosity

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The default lg10Lum parameter value has changed from -40 to 40
       to match XSPEC.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Emin
        Minimum energy over which the luminosity is calculated.
    Emax
        Maximum energy over which the luminosity is calculated.
    Redshift
        redshift of the source
    lg10Lum
        log (base 10) of the luminosity in erg/s

    See Also
    --------
    XScflux, XScpflux

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    See [1]_ for the meaning and restrictions, in particular the
    necessity of freezing the amplitude, or normalization, of the
    emission component (or components) at 1.

    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelClumin.html

    Examples
    --------

    With the following definitions:

    >>> clumin = XSclumin()
    >>> absmdl = XSphabs()
    >>> plmdl = XSpowerlaw()
    >>> gmdl = XSgaussian()
    >>> srcmdl = plmdl + gmdl

    then the model can be applied in a number of ways, such as:

    >>> mdl1 = clumin(absmdl * srcmdl)
    >>> mdl2 = absmdl * clumin(srcmdl)
    >>> mdl3 = absmdl * (plmdl + clumin(gmdl))

    """

    __function__ = "C_clumin"

    def __init__(self, name='xsclumin'):
        self.Emin = XSParameter(name, 'Emin', 0.5, min=0.0, max=1e6,
                                hard_min=0.0, hard_max=1e6, frozen=True,
                                units='keV')
        self.Emax = XSParameter(name, 'Emax', 10.0, min=0.0, max=1e6,
                                hard_min=0.0, hard_max=1e6, frozen=True,
                                units='keV')
        self.Redshift = XSParameter(name, 'Redshift', 0, min=-0.999, max=10,
                                    hard_min=-0.999, hard_max=10, frozen=True)
        self.lg10Lum = XSParameter(name, 'lg10Lum', 40.0, min=-100.0,
                                   max=100.0, hard_min=-100.0, hard_max=100.0,
                                   frozen=False, units='cgs')
        XSConvolutionKernel.__init__(self, name, (self.Emin,
                                                  self.Emax,
                                                  self.Redshift,
                                                  self.lg10Lum
                                                  ))


class XScpflux(XSConvolutionKernel):
    """The XSPEC cpflux convolution model: calculate photon flux

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Emin
        Minimum energy over which the photon flux is calculated.
    Emax
        Maximum energy over which the photon flux is calculated.
    Flux
        photon flux in photon/cm^2/s

    See Also
    --------
    XScflux, XSclumin

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    See [1]_ for the meaning and restrictions, in particular the
    necessity of freezing the amplitude, or normalization, of the
    emission component (or components) at 1.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelCpflux.html

    """

    _calc = _xspec.C_cpflux

    def __init__(self, name='xscpflux'):
        self.Emin = XSParameter(name, 'Emin', 0.5, min=0.0, max=1e6,
                                hard_min=0.0, hard_max=1e6, frozen=True,
                                units='keV')
        self.Emax = XSParameter(name, 'Emax', 10.0, min=0.0, max=1e6,
                                hard_min=0.0, hard_max=1e6, frozen=True,
                                units='keV')
        self.Flux = XSParameter(name, 'Flux', 1.0, min=0.0, max=1e10,
                                hard_min=0.0, hard_max=1e10,
                                frozen=False, units='')
        XSConvolutionKernel.__init__(self, name, (self.Emin,
                                                  self.Emax,
                                                  self.Flux
                                                  ))


class XSgsmooth(XSConvolutionKernel):
    """The XSPEC gsmooth convolution model: gaussian smoothing

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The XSPEC model parameter was renamed from Sig@6keV to Sig_6keV
       so this is now the name of the parameter. The old name (SigAt6keV)
       can still be used but will be removed in a future release.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Sig_6keV
        gaussian sigma at 6 keV
    Index
        power of energy for sigma variation

    See Also
    --------
    XSlsmooth

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelGsmooth.html

    """

    _calc = _xspec.C_gsmooth

    def __init__(self, name='xsgsmooth'):
        self.Sig_6keV = XSParameter(name, 'Sig_6keV', 1.0, min=0.0, max=10.0,
                                    hard_min=0.0, hard_max=20.0,
                                    frozen=False, units='keV', aliases=['SigAt6keV'])
        self.Index = XSParameter(name, 'Index', 0.0, min=-1.0, max=1.0,
                                 hard_min=-1.0, hard_max=1.0, frozen=True)
        XSConvolutionKernel.__init__(self, name, (self.Sig_6keV, self.Index))


class XSireflect(XSConvolutionKernel):
    """The XSPEC ireflect convolution model: reflection from ionized material

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    rel_refl
        The reflection scaling factor (1 for isotropic source above disk).
    Redshift
        The redshift of the source.
    abund
        The abundance of elements heavier than He relative to their solar
        abundances.
    Fe_abund
        The iron abundance relative to the above.
    cosIncl
        The cosine of the inclination angle.
    T_disk
        The disk temperature in K.
    xi
        The disk ionization parameter: see [1]_ for an explanation.

    See Also
    --------
    XSbexriv, XSpexriv

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the IREFLECT_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelIreflect.html

    """

    _calc = _xspec.C_ireflct

    def __init__(self, name='xsireflect'):
        self.rel_refl = XSParameter(name, 'rel_refl', 0.0, min=-1.0, max=1e6,
                                    hard_min=-1.0, hard_max=1e6, frozen=False)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, min=-0.999, max=10.0,
                                    hard_min=-0.999, hard_max=10.0, frozen=True)
        self.abund = XSParameter(name, 'abund', 1.0, min=0.0, max=1e6,
                                 hard_min=0.0, hard_max=1e6, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1.0, min=0.0, max=1e6,
                                    hard_min=0.0, hard_max=1e6, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.45, min=0.05, max=0.95,
                                   hard_min=0.05, hard_max=0.95, frozen=True)
        self.T_disk = XSParameter(name, 'T_disk', 3e4, min=1e4, max=1e6,
                                  hard_min=1e4, hard_max=1e6, frozen=True,
                                  units='K')
        self.xi = XSParameter(name, 'xi', 1.0, min=0.0, max=1e3, hard_min=0.0,
                              hard_max=5e3, frozen=True, units='erg cm/s')
        XSConvolutionKernel.__init__(self, name, (self.rel_refl,
                                                  self.Redshift,
                                                  self.abund,
                                                  self.Fe_abund,
                                                  self.cosIncl,
                                                  self.T_disk,
                                                  self.xi
                                                  ))


class XSkdblur(XSConvolutionKernel):
    """The XSPEC kdblur convolution model: convolve with the laor model

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Index
        The power law dependence of emissivity (scales as R^-Index).
    Rin_G
        The inner radius, in units of GM/c^2.
    Rout_G
        The outer radius, in units of GM/c^2.
    Incl
        The disk inclination angle, in degrees. A face-on disk has
        Incl=0.

    See Also
    --------
    XSkdblur2, XSlaor

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKdblur.html

    """

    _calc = _xspec.C_kdblur

    def __init__(self, name='xskdblur'):
        self.Index = XSParameter(name, 'Index', 3.0, min=-10.0, max=10.0,
                                 hard_min=-10.0, hard_max=10.0, frozen=True)
        self.Rin_G = XSParameter(name, 'Rin_G', 4.5, min=1.235, max=400.0,
                                 hard_min=1.235, hard_max=400.0, frozen=True)
        self.Rout_G = XSParameter(name, 'Rout_G', 100.0, min=1.235, max=400.0,
                                  hard_min=1.235, hard_max=400.0, frozen=True)
        self.Incl = XSParameter(name, 'Incl', 30.0, min=0.0, max=90.0,
                                hard_min=0.0, hard_max=90.0, frozen=False,
                                units='deg')
        XSConvolutionKernel.__init__(self, name, (self.Index,
                                                  self.Rin_G,
                                                  self.Rout_G,
                                                  self.Incl
                                                  ))


class XSkdblur2(XSConvolutionKernel):
    """The XSPEC kdblur2 convolution model: convolve with the laor2 model

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
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

    See Also
    --------
    XSkdblur, XSlaor2

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKdblur2.html

    """

    _calc = _xspec.C_kdblur2

    def __init__(self, name='xskdblur2'):
        self.Index = XSParameter(name, 'Index', 3.0, min=-10.0, max=10.0,
                                 hard_min=-10.0, hard_max=10.0, frozen=True)
        self.Rin_G = XSParameter(name, 'Rin_G', 4.5, min=1.235, max=400.0,
                                 hard_min=1.235, hard_max=400.0, frozen=True)
        self.Rout_G = XSParameter(name, 'Rout_G', 100.0, min=1.235, max=400.0,
                                  hard_min=1.235, hard_max=400.0, frozen=True)
        self.Incl = XSParameter(name, 'Incl', 30.0, min=0.0, max=90.0,
                                hard_min=0.0, hard_max=90.0, frozen=False,
                                units='deg')
        self.Rbreak = XSParameter(name, 'Rbreak', 20.0, min=1.235, max=400.0,
                                  hard_min=1.235, hard_max=400.0, frozen=True)
        self.Index1 = XSParameter(name, 'Index1', 3.0, min=-10.0, max=10.0,
                                  hard_min=-10.0, hard_max=10.0, frozen=True)
        XSConvolutionKernel.__init__(self, name, (self.Index,
                                                  self.Rin_G,
                                                  self.Rout_G,
                                                  self.Incl,
                                                  self.Rbreak,
                                                  self.Index1))


class XSkerrconv(XSConvolutionKernel):
    """The XSPEC kerrconv convolution model: accretion disk line shape with BH spin as free parameter

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The first two parameters have been renamed `Index1` and
       `Index2` to match the XSPEC definition (they had been called
       `Index` and `Index1`).

    .. versionadded:: 4.12.2

    Attributes
    ----------
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

    See Also
    --------
    XSdiskline, XSkerrdisk, XSlaor

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelKerrconv.html

    """

    _calc = _xspec.C_spinconv

    def __init__(self, name='xskerrconv'):
        self.Index1 = XSParameter(name, 'Index1', 3.0, min=-10.0, max=10.0,
                                 hard_min=-10.0, hard_max=10.0, frozen=True)
        self.Index2 = XSParameter(name, 'Index2', 3.0, min=-10.0, max=10.0,
                                  hard_min=-10.0, hard_max=10.0, frozen=True)
        self.r_br_g = XSParameter(name, 'r_br_g', 6.0, min=1.0, max=400.0,
                                  hard_min=1.0, hard_max=400.0, frozen=True)
        self.a = XSParameter(name, 'a', 0.998, min=0.0, max=0.998,
                             hard_min=0.0, hard_max=0.998, frozen=False)
        self.Incl = XSParameter(name, 'Incl', 30.0, min=0.0, max=90.0,
                                hard_min=0.0, hard_max=90.0, frozen=False,
                                units='deg')
        self.Rin_ms = XSParameter(name, 'Rin_ms', 1.0, min=1.0, max=400.0,
                                  hard_min=1.0, hard_max=400.0, frozen=True)
        self.Rout_ms = XSParameter(name, 'Rout_ms', 400.0, min=1.0, max=400.0,
                                   hard_min=1.0, hard_max=400.0, frozen=True)
        XSConvolutionKernel.__init__(self, name, (self.Index1,
                                                  self.Index2,
                                                  self.r_br_g,
                                                  self.a,
                                                  self.Incl,
                                                  self.Rin_ms,
                                                  self.Rout_ms
                                                  ))


class XSlsmooth(XSConvolutionKernel):
    """The XSPEC lsmooth convolution model: lorentzian smoothing

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The XSPEC model parameter was renamed from Sig@6keV to Sig_6keV
       so this is now the name of the parameter. The old name (SigAt6keV)
       can still be used but will be removed in a future release.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Sig_6keV
        lorentzian sigma at 6 keV
    Index
        power of energy for sigma variation

    See Also
    --------
    XSgsmooth

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelLsmooth.html

    """

    _calc = _xspec.C_lsmooth

    def __init__(self, name='xslsmooth'):
        self.Sig_6keV = XSParameter(name, 'Sig_6keV', 1.0, min=0.0, max=10.0,
                                    hard_min=0.0, hard_max=20.0, frozen=False,
                                    units='keV', aliases=['SigAt6keV'])
        self.Index = XSParameter(name, 'Index', 0.0, min=-1.0, max=1.0,
                                 hard_min=-1.0, hard_max=1.0, frozen=True)
        XSConvolutionKernel.__init__(self, name, (self.Sig_6keV, self.Index))


class XSpartcov(XSConvolutionKernel):
    """The XSPEC partcov convolution model: partial covering

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    CvrFract
        The covering fraction (0 to 1).

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelPartcov.html

    """
    _calc = _xspec.C_PartialCovering

    def __init__(self, name='xspartcov'):
        self.CvrFract = XSParameter(name, 'CvrFract', 0.5, min=0.05, max=0.95,
                                  hard_min=0.0, hard_max=1.0, frozen=False)
        XSConvolutionKernel.__init__(self, name, (self.CvrFract,))


class XSrdblur(XSConvolutionKernel):
    """The XSPEC rdblur convolution model: convolve with the diskline model shape

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Betor10
        The power law dependence of emissivity: see [1]_ for more details.
    Rin_M
        The inner radius, in units of GM^2/c.
    Rout_M
        The outer radius, in units of GM^2/c.
    Incl
        The inclination, in degrees.

    See Also
    --------
    XSdiskline

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRdblur.html

    """
    _calc = _xspec.C_rdblur

    def __init__(self, name='xsrdblur'):
        self.Betor10 = XSParameter(name, 'Betor10', -2.0, min=-10.0, max=20.0,
                                   hard_min=-10.0, hard_max=20.0, frozen=True)
        self.Rin_M = XSParameter(name, 'Rin_M', 10.0, min=6.0, max=1000.0,
                                 hard_min=6.0, hard_max=10000.0, frozen=True)
        self.Rout_M = XSParameter(name, 'Rout_M', 1000.0, min=0.0,
                                  max=1000000.0, hard_min=0.0,
                                  hard_max=10000000.0, frozen=True)
        self.Incl = XSParameter(name, 'Incl', 30.0, min=0.0, max=90.0,
                                hard_min=0.0, hard_max=90.0, frozen=False,
                                units='deg')
        XSConvolutionKernel.__init__(self, name, (self.Betor10,
                                                  self.Rin_M,
                                                  self.Rout_M,
                                                  self.Incl
                                                  ))


class XSreflect(XSConvolutionKernel):
    """The XSPEC reflect convolution model: reflection from neutral material

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    rel_refl
        The reflection scaling parameter (a value between 0 and 1
        for an isotropic source above the disk, less than 0 for no
        reflected component).
    Redshift
        The redshift of the source.
    abund
        The abundance of the elements heavier than He relative to their
        solar abundance, as set by the ``set_xsabund`` function.
    Fe_abund
        The iron abundance relative to the solar abundance, as set by
        the ``set_xsabund`` function.
    cosIncl
        The cosine of the inclination angle in degrees.

    See Also
    --------
    XSbexrav, XSpexrav

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the REFLECT_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelReflect.html

    """
    _calc = _xspec.C_reflct

    def __init__(self, name='xsreflect'):
        self.rel_refl = XSParameter(name, 'rel_refl', 0.0, min=-1.0, max=1e6,
                                    hard_min=-1.0, hard_max=1e6, frozen=False)
        self.Redshift = XSParameter(name, 'Redshift', 0.0, min=-0.999, max=10.0,
                                    hard_min=-0.999, hard_max=10.0, frozen=True)
        self.abund = XSParameter(name, 'abund', 1.0, min=0.0, max=1e6,
                                 hard_min=0.0, hard_max=1e6, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1.0, min=0.0, max=1e6,
                                    hard_min=0.0, hard_max=1e6, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.45, min=0.05, max=0.95,
                                   hard_min=0.05, hard_max=0.95, frozen=True)
        XSConvolutionKernel.__init__(self, name, (self.rel_refl,
                                                  self.Redshift,
                                                  self.abund,
                                                  self.Fe_abund,
                                                  self.cosIncl
                                                  ))


@version_at_least("12.9.1")
class XSrfxconv(XSConvolutionKernel):
    """The XSPEC rfxconv convolution model: angle-dependent reflection from an ionized disk

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    rel_refl
        The relative reflection normalization. If negative then only
        the reflected component is returned.
    redshift
        The redshift of the source.
    Fe_abund
        The iron abundance relative to the solar abundance, as set by
        the ``set_xsabund`` function. All other elements are assumed to
        have Solar abundances.
    cosIncl
        The cosine of the inclination angle in degrees.
    log_xi
        The ionization parameter used by the table models.

    See Also
    --------
    XSxilconv

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the RFXCONV_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRfxconv.html

    """

    __function__ = "C_rfxconv"

    def __init__(self, name='xsrfxconv'):
        self.rel_refl = XSParameter(name, 'rel_refl', -1.0, min=-1.0, max=1e6,
                                    hard_min=-1.0, hard_max=1e6)
        self.redshift = XSParameter(name, 'redshift', 0.0, min=0.0, max=4.0,
                                    hard_min=0.0, hard_max=4.0, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1.0, min=0.5, max=3,
                                    hard_min=0.5, hard_max=3, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.5, min=0.05, max=0.95,
                                   hard_min=0.05, hard_max=0.95, frozen=True)
        self.log_xi = XSParameter(name, 'log_xi', 1.0, min=1.0, max=6.0,
                                  hard_min=1.0, hard_max=6.0)
        XSConvolutionKernel.__init__(self, name, (self.rel_refl,
                                                  self.redshift,
                                                  self.Fe_abund,
                                                  self.cosIncl,
                                                  self.log_xi
                                                  ))


class XSrgsxsrc(XSConvolutionKernel):
    """The XSPEC rgsxsrc convolution model: convolve an RGS spectrum for extended emission.

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    order
        The order, which must be -1 to -3 inclusive.

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    The ``set_xsxset`` function must be used to set the RGS_XSOURCE_FILE
    value to point to a file as described in [1]_.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelRgxsrc.html

    """

    _calc = _xspec.rgsxsrc

    def __init__(self, name='xsrgsxsrc'):
        self.order = XSParameter(name, 'order', -1.0, min=-3.0, max=-1,
                                 hard_min=-3.0, hard_max=-1, frozen=True)
        XSConvolutionKernel.__init__(self, name, (self.order,))


class XSsimpl(XSConvolutionKernel):
    """The XSPEC simpl convolution model: comptonization of a seed spectrum.

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Gamma
        The photon power law index.
    FracSctr
        The scattered fraction (between 0 and 1).
    UpScOnly
        A flag to switch between up-scattering only (when greater than
        zero) or both up- and down-scattering (when zero or less).

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelSimpl.html

    """

    _calc = _xspec.C_simpl

    def __init__(self, name='xssimpl'):
        self.Gamma = XSParameter(name, 'Gamma', 2.3, min=1.1, max=4.0,
                                 hard_min=1.0, hard_max=5.0, frozen=False)
        self.FracSctr = XSParameter(name, 'FracSctr', 0.05, min=0.0, max=0.4,
                                    hard_min=0.0, hard_max=1.0, frozen=False)
        self.UpScOnly = XSParameter(name, 'UpScOnly', 1.0, min=0.0, max=100.0,
                                    hard_min=0.0, hard_max=100.0, frozen=True)
        XSConvolutionKernel.__init__(self, name, (self.Gamma,
                                                  self.FracSctr,
                                                  self.UpScOnly
                                                  ))


@version_at_least("12.11.0")
class XSthcomp(XSConvolutionKernel):
    """The XSPEC thcomp convolution model: Thermally comptonized continuum.

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    .. note:: Parameter renames in XSPEC 12.11.1

       In XSPEC 12.11.1 the parameters are now called ``Gamma_tau`` and
       ``cov_frac`` rather than ``gamma_tau`` and ``FracSctr``. The case of
       ``Gamma_tau`` does not matter for Sherpa, but the scattering-fraction
       parameter has been renamed (with an alias in place for users used to
       the XSPEC 12.11.0 names).

    Attributes
    ----------
    Gamma_tau
        The low-energy power-law photon index when positive and the Thomson
        optical depth (multiplied by -1) when negative.
    kT_e
        The electron temperature (high energy rollover)
    cov_frac
        The scattering fraction (between 0 and 1). If 1 then all of the
        seed photons will be Comptonized, and if 0 then only the original
        seed photons will be seen.
    z
        redshift

    See Also
    --------
    XSnthcomp

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    This model is only available when used with XSPEC 12.11.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelThcomp.html

    """

    __function__ = "thcompf"

    def __init__(self, name='xsthcomp'):
        # TODO: allow negative
        self.Gamma_tau = XSParameter(name, 'Gamma_tau', 1.7, min=1.001, max=5.0,
                                     hard_min=1.001, hard_max=10.0, frozen=False)
        self.kT_e = XSParameter(name, 'kT_e', 50.0, min=0.5, max=150.0,
                                hard_min=0.5, hard_max=150.0, units='keV',
                                frozen=False)
        # In XSPEC 12.11.0 the parameter was named FracSctr, and thanks to a typo
        # it was stored as the parameter FracStr.
        self.cov_frac = XSParameter(name, 'cov_frac', 1.0, min=0.0, max=1.0,
                                    hard_min=0.0, hard_max=1.0, frozen=False,
                                    aliases=["FracSctr", "FracStr"])
        self.z = XSParameter(name, 'z', 0.0, min=0.0, max=5.0,
                             hard_min=0.0, hard_max=5.0, frozen=True)

        XSConvolutionKernel.__init__(self, name, (self.Gamma_tau,
                                                  self.kT_e,
                                                  self.cov_frac,
                                                  self.z
                                                  ))


@version_at_least("12.9.1")
class XSvashift(XSConvolutionKernel):
    """The XSPEC vashift convolution model: velocity shift an additive model.

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The Velocity parameter was incorrectly labelled as Redshift:
       both names can be used to access this parameter.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Velocity
        The velocity, in km/s.

    See Also
    --------
    XSvmshift, XSzashift

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelVashift.html

    """

    __function__ = "C_vashift"

    def __init__(self, name='xsvashift'):
        self.Velocity = XSParameter(name, 'Velocity', 0.0,
                                    min=-1e4, max=1e4,
                                    hard_min=-1e4, hard_max=1e4,
                                    units='km/s', frozen=True,
                                    # Parameter was mis-labelled until 4.14.0
                                    aliases=['Redshift'])
        XSConvolutionKernel.__init__(self, name, (self.Velocity,))


@version_at_least("12.9.1")
class XSvmshift(XSConvolutionKernel):
    """The XSPEC vmshift convolution model: velocity shift a multiplicative model.

    The model is described at [1]_.

    .. versionchanged:: 4.14.0
       The Velocity parameter was incorrectly labelled as Redshift:
       both names can be used to access this parameter.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Velocity
        The velocity, in km/s.

    See Also
    --------
    XSvashift, XSzmshift

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelVmshift.html

    """

    __function__ = "C_vmshift"

    def __init__(self, name='xsvmshift'):
        self.Velocity = XSParameter(name, 'Velocity', 0.0,
                                    min=-1e4, max=1e4,
                                    hard_min=-1e4, hard_max=1e4,
                                    units='km/s', frozen=True,
                                    # Parameter was mis-labelled until 4.14.0
                                    aliases=['Redshift'])
        XSConvolutionKernel.__init__(self, name, (self.Velocity,))


@version_at_least("12.9.1")
class XSxilconv(XSConvolutionKernel):
    """The XSPEC xilconv convolution model: angle-dependent reflection from an ionized disk

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    rel_refl
        The relative reflection normalization. If negative then only
        the reflected component is returned.
    redshift
        The redshift of the source.
    Fe_abund
        The iron abundance relative to the solar abundance, as set by
        the ``set_xsabund`` function. All other elements are assumed to
        have Solar abundances.
    cosIncl
        The cosine of the inclination angle in degrees.
    log_xi
        The ionization parameter used by the table models.
    cutoff
        The exponential cut-off energy, in keV.

    See Also
    --------
    XSrfxconv

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    The precision of the numerical integration can be changed by using
    the ``set_xsxset`` function to set the value of the XILCONV_PRECISION
    keyword, which defines the fractional precision. The default is 0.01
    (1%).

    This model is only available when used with XSPEC 12.9.1 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelXilconv.html

    """

    __function__ = "C_xilconv"

    def __init__(self, name='xsxilconv'):
        self.rel_refl = XSParameter(name, 'rel_refl', -1.0, min=-1.0, max=1e6,
                                    hard_min=-1.0, hard_max=1e6)
        self.redshift = XSParameter(name, 'redshift', 0.0, min=0.0, max=4.0,
                                    hard_min=0.0, hard_max=4.0, frozen=True)
        self.Fe_abund = XSParameter(name, 'Fe_abund', 1.0, min=0.5, max=3.0,
                                    hard_min=0.5, hard_max=3.0, frozen=True)
        self.cosIncl = XSParameter(name, 'cosIncl', 0.5, min=0.05, max=0.95,
                                   hard_min=0.05, hard_max=0.95, frozen=True)
        self.log_xi = XSParameter(name, 'log_xi', 1.0, min=1.0, max=6,
                                  hard_min=1.0, hard_max=6)
        self.cutoff = XSParameter(name, 'cutoff', 300.0, min=20.0, max=300.0,
                                  hard_min=20.0, hard_max=300.0,
                                  units='keV', frozen=True)
        XSConvolutionKernel.__init__(self, name, (self.rel_refl,
                                                  self.redshift,
                                                  self.Fe_abund,
                                                  self.cosIncl,
                                                  self.log_xi,
                                                  self.cutoff
                                                  ))


class XSzashift(XSConvolutionKernel):
    """The XSPEC zashift convolution model: redshift an additive model.

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Redshift
        The redshift.

    See Also
    --------
    XSvashift, XSzmshift

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZashift.html

    """

    _calc = _xspec.C_zashift

    def __init__(self, name='xszashift'):
        self.Redshift = XSParameter(name, 'Redshift', 0.0, min=-0.999, max=10.0,
                                    hard_min=-0.999, hard_max=10, frozen=True)
        XSConvolutionKernel.__init__(self, name, (self.Redshift,))


class XSzmshift(XSConvolutionKernel):
    """The XSPEC zmshift convolution model: redshift a multiplicative model.

    The model is described at [1]_.

    .. versionadded:: 4.12.2

    Attributes
    ----------
    Redshift
        The redshift.

    See Also
    --------
    XSvmshift, XSzashift

    Notes
    -----
    Unlike XSPEC, the convolution model is applied directly to the model, or
    models, rather than using the multiplication symbol.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelZmshift.html

    """

    _calc = _xspec.C_zmshift

    def __init__(self, name='xszmshift'):
        self.Redshift = XSParameter(name, 'Redshift', 0.0, min=-0.999, max=10.0,
                                    hard_min=-0.999, hard_max=10, frozen=True)
        XSConvolutionKernel.__init__(self, name, (self.Redshift,))


@version_at_least("12.11.0")
class XSbwcycl(XSAdditiveModel):
    """The XSPEC bwcycl model: Becker-Wolff self-consistent cyclotron line model.

    The model is described at [1]_. Please review the restrictions
    on the model and parameter values at this reference before using
    the model.

    .. versionchanged:: 4.14.0
       The D parameter now displays the correct units (kpc) rather
       than km.

    Attributes
    ----------
    Radius
        The radius of the Neutron star, in km. Keep frozen.
    Mass
        The mass of the Neutron star, in solar units. Keep frozen.
    csi
        Parameter linked to the photon escape time (order of some
        unities).
    delta
        Ratio between bulk and thermal Comptonization importances.
    B
        The magnetic field in units of 10^12 G.
    Mdot
        The mass accretion rate, in 10^17 g/s.
    Te
        The electron temperature in units of keV.
    r0
        The column radius in m.
    D
        The source distance in kpc. Keep frozen.
    BBnorm
        The normalization of the blackbody seed photon component
        (fix it to zero at first).
    CYCnorm
        The normalization of the cyclotron emission seed photon
        component (fix it to one).
    FFnorm
        The normalization of the Bremsstrahlung emission seed photon
        component (fix it to one).
    norm
        The normalization of the model (fix it to one).

    Notes
    -----
    This model is only available when used with XSPEC 12.11.0 or later.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelBwcycl.html

    """

    __function__ = "beckerwolff"  # "c_beckerwolff"  do not have a direct interface to c_xxx

    def __init__(self, name='bwcycl'):
        self.Radius = XSParameter(name, 'Radius', 10, 5, 20, 5, 20, units='km', frozen=True)
        self.Mass = XSParameter(name, 'Mass', 1.4, 1, 3, 1, 3, units='Solar', frozen=True)
        self.csi = XSParameter(name, 'csi', 1.5, 0.01, 20, 0.01, 20)
        self.delta = XSParameter(name, 'delta', 1.8, 0.01, 20, 0.01, 20)
        self.B = XSParameter(name, 'B', 4, 0.01, 100, 0.01, 100, units='1e12G')
        self.Mdot = XSParameter(name, 'Mdot', 1, 1e-6, 1e6, 1e-6, 1e6, units='1e17g/s')
        self.Te = XSParameter(name, 'Te', 5, 0.1, 100, 0.1, 100, units='keV')
        self.r0 = XSParameter(name, 'r0', 44, 10, 1000, 10, 1000, units='m')
        self.D = XSParameter(name, 'D', 5, 1, 20, 1, 20, units='kpc', frozen=True)
        self.BBnorm = XSParameter(name, 'BBnorm', 0.0, 0, 100, 0, 100, frozen=True)
        self.CYCnorm = XSParameter(name, 'CYCnorm', 1.0, -1, 100, -1, 100, frozen=True)
        self.FFnorm = XSParameter(name, 'FFnorm', 1.0, -1, 100, -1, 100, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval, frozen=True)

        XSAdditiveModel.__init__(self, name, (self.Radius, self.Mass, self.csi, self.delta,
                                              self.B, self.Mdot, self.Te, self.r0, self.D,
                                              self.BBnorm, self.CYCnorm, self.FFnorm,
                                              self.norm))


# Add model classes to __all__
#
# Should this remove the "base" classes, such as
# XSModel and XSConvolutionModel?
#
__all__ += tuple(n for n in globals() if n.startswith('XS'))
