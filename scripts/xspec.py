#
#  Copyright (C) 2012-2016, 2020-2026
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

"""Parser for XSPEC model files.

The XSPEC library [1]_ uses ASCII files to define models [2]_, and it
can be useful to be able to read these files either to identify
changes to the Sherpa code to support a new XSPEC release [3]_ or for
writing a module for an XSPEC user model.

.. versionchanged:: 4.19.0
   The interface to a particular model is now named after the model
   name rather than the actual function name (e.g. the apec model is
   now called "apec" rather than something like "C_apec"). This makes
   it easier to identify which routine to call. The docstring for the
   model includes the name of the function and the number of
   parameters it requires. Several symbols related to XSPEC versions
   have been added to this module. There is now support for handling
   model definitions including the grad=xxx argument added in XSPEC
   13.0.0. Basic support for downloading data files for models has
   been added.

References
----------

.. [1] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/index.html

.. [2] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html

.. [3] https://sherpa.readthedocs.io/en/latest/developer/index.html#update-the-xspec-bindings

"""

from collections import Counter
from collections.abc import Callable, Sequence
from contextlib import suppress
from dataclasses import dataclass
import logging
from pathlib import Path
import re
import string
import time


__all__ = ("SUPPORTED_VERSIONS", "MIN_VERSION", "MAX_VERSION",
           "XSPECcode", "XSPECDataRow",
           "parse_xspec_model_description",
           "create_xspec_code", "get_version",
           "get_model_data_version_path",
           "read_model_data_version_file",
           "find_missing_model_data_files",
           "download_missing_model_data_files"
           )


debug = logging.getLogger(__name__).debug
info = logging.getLogger(__name__).info
warning = logging.getLogger(__name__).warning

# Represent the XSPEC version (without a patch level).
Version = tuple[int, int, int]

# I am not sure what the naming of the XSPEC components are, but let's
# stick with major, minor, and micro. We drop the patch level - e.g.
# "c" in "12.12.0c" as that is not helpful to track here.
#
SUPPORTED_VERSIONS: list[Version] = [
    (12, 13, 0), (12, 13, 1),
    (12, 14, 0), (12, 14, 1),
    (12, 15, 0), (12, 15, 1),
    (13, 0, 0)
]
"""What versions of XSPEC are supported by Sherpa?

Newer versions of XSPEC may be usable with Sherpa but there is no
guarantee of support (e.g. some models may not be usable). It is
unlikely that older versions will be usable.

"""

# We could use packaging.versions.Version here, but for our needs we
# can get away with a tuple of integers. That is, we do not need the
# full support for PEP-440.
#
MIN_VERSION = min(SUPPORTED_VERSIONS)
"""The minimum supported XSPEC version."""

MAX_VERSION = max(SUPPORTED_VERSIONS)
"""The maximum supported XSPEC version."""


def get_version(version: str) -> Version:
    """Strip out any XSPEC patch level.

    So '12.12.0c' gets converted to '12.12.0', and then to (12, 12,
    0). This is helpful as then it makes version comparison easier, as
    we can rely on the standard tuple ordering.

    Parameters
    ----------
    version : str
        The XSPEC version string, of the form "12.12.0c", so it can
        include the XSPEC patch level.

    Returns
    -------
    (major, minor, micro) : tuple of int
        The XSPEC patchlevel is ignored.

    """

    # Do not bother decoding the patch level. XSPEC seems to stop
    # parsing as soon as a non-numeric character is hit, which this
    # regexp also does.
    #
    matches = re.search(r'^(\d+)\.(\d+)\.(\d+)', version)
    if matches is None:
        raise ValueError(f"Invalid XSPEC version string: {version}")

    return (int(matches[1]), int(matches[2]), int(matches[3]))


@dataclass
class XSPECcode:
    """The code components needed to compile the XSPEC user model."""

    python: str
    """The Python code"""

    compiled: str
    """The C++ code"""


class ModelDefinition:
    """Represent the model definition from an XSPEC model file.

    Parameters
    ----------
    name : str
       The model name.
    clname : str
       The class name used to represent this model in Sherpa.
    funcname : str
       The name of the function from the model file (so it should
       include any prefix like C_).
    flags : sequence of int
       The flags value.
    elo : float
       The minimum energy supported by this model (unused).
    ehi : float
       The maximum energy supported by this model (unused).
    pars : sequence of ParameterDefinition
       Any parameter values. It is expected this is not empty.
    initString : str or None, optional
        The default string to send to the model.
    grad : str or None, optional
        The grad argument.

    See Also
    --------
    AddModelDefinition, MulModelDefinition, ConModelDefinition,
    MixModelDefinition, AcnModelDefinition, AmxModelDefinition

    Notes
    -----
    Do not instantiate this class directly.

    """

    modeltype: str
    language: str

    def __init__(self, name: str, clname: str, funcname: str,
                 flags: list[int], elo: float, ehi: float,
                 pars: list["ParameterDefinition"],
                 initString: str | None = None,
                 grad: str | None = None
                 ) -> None:
        assert self.modeltype is not None, \
            "ModelDefinition should not be directly created."
        self.name = name
        self.clname = clname
        self.funcname = funcname
        self.flags = flags
        self.elo = elo
        self.ehi = ehi
        self.pars = pars

        # This will probably need to be changed if mixing models
        # (mix or amx) are supported.
        #
        # The use of strings for the language is not ideal; really
        # should use some form of an enumeration.
        #
        # Note that FORTRAN function names are converted to lower case
        # as the name should be case insensitive but I found some
        # issues when mixed-case was used.
        #
        if self.funcname.startswith('F_'):
            self.language = 'Fortran - double precision'
            self.funcname = self.funcname[2:].lower()
        elif self.funcname.startswith('c_'):
            self.language = 'C style'
            self.funcname = self.funcname[2:]
        elif self.funcname.startswith('C_'):
            self.language = 'C++ style'
            self.funcname = self.funcname[2:]
        else:
            self.language = 'Fortran - single precision'
            self.funcname = self.funcname.lower()

        if initString is not None and self.language.startswith('F'):
            initString = None

        self.initString = initString
        self.grad = grad

    def __repr__(self) -> str:
        return f"<{self.modeltype}:{self.name}:{self.funcname}:{len(self.pars)} pars>"

    def __str__(self) -> str:
        pars = "\n".join([str(p) for p in self.pars])
        return f"{self.modeltype}.{self.name} " +  \
            f"function={self.funcname}\n{self.language}\n{pars}"


class AddModelDefinition(ModelDefinition):
    """XSPEC additive models.

    See [1]_ for examples.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Additive.html
    """

    modeltype = "Add"


class MulModelDefinition(ModelDefinition):
    """XSPEC multiplicative models.

    See [1]_ for examples.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Multiplicative.html
    """

    modeltype = "Mul"


class ConModelDefinition(ModelDefinition):
    """XSPEC convolution models.

    See [1]_ for examples.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Convolution.html
    """

    modeltype = "Con"


class MixModelDefinition(ModelDefinition):
    """XSPEC mixing models.

    See [1]_ for examples. These are currently unsupported in Sherpa.

    References
    ----------

    .. [1] https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Mixing.html
    """

    modeltype = "Mix"


class AcnModelDefinition(ModelDefinition):
    """XSPEC Acn model: pile-up models.

    These are currently unsupported in Sherpa.

    """

    modeltype = "Acn"


# Found in looking through
#   heasoft-6.16/Xspec/src/tools/initpackage/ModelMap.cxx
class AmxModelDefinition(ModelDefinition):
    """XSPEC Amx model: a combination of mixing and pile-up models.

    These are currently unsupported in Sherpa.

    """
    modeltype = "Amx: apparently a combination of mixing and pile-up models"


class ParameterDefinition:
    """Represent an XSPEC parameter.

    Parameters
    ----------
    name : str
        The parameter name.
    default : float
        The default value
    units : str or None, optional
        The unit field. There is no check this meets any standard.
    softmin, softmax, hardmin, hardmax : float or None
        The minimum and maximum values for the parameter (using the
        XSPEC definition of soft and hard, not Sherpa).
    delta : float or None, optional
        The delta parameter. At present this is only used to determine
        if the parameter is frozen by default (delta < 0).

    See Also
    --------
    BasicParametertDefinition, SwitchParameterDefinition,
    ScaleParameterDefinition

    Notes
    -----
    Do not instantiate this class directly.

    We are missing support for periodic parameters (that is parameters
    that end with a P) as it is unclear how to handle them in Sherpa.

    """

    paramtype: str

    def __init__(self,
                 name: str,
                 default: float | int,
                 units: str | None = None,
                 *,
                 softmin: float | None = None,
                 softmax: float | None = None,
                 hardmin: float | None = None,
                 hardmax: float | None = None,
                 delta: float | None = None
                 ) -> None:
        assert self.paramtype is not None, \
            'ParameterDefinition should not be directly created'

        self.name = name
        self.default = default
        self.units = units

        self.softmin = softmin
        self.softmax = softmax
        self.hardmin = hardmin
        self.hardmax = hardmax
        self.delta = None if delta is None else abs(delta)

    def __str__(self) -> str:
        return f"{self.name} = {self.default}"

    def param_string(self) -> str:
        out = f"XSParameter(name, '{self.name}', {self.default}"

        for (pval, pname) in [(self.softmin, "min"),
                              (self.softmax, "max"),
                              (self.hardmin, "hard_min"),
                              (self.hardmax, "hard_max")]:
            if pval is not None:
                out += f", {pname}={pval}"

        if self.units is not None:
            out += f", units='{self.units}'"

        out += ", alwaysfrozen=True)"
        return out


class SwitchParameterDefinition(ParameterDefinition):
    """A "switch" parameter.

    These are for parameter values that change how the model evaluates
    and are not changed during a fit.

    """

    paramtype = "Switch"


# Do we handle this type of parameter correctly?
#
class ScaleParameterDefinition(ParameterDefinition):
    """A "scale" parameter.
    """

    paramtype = "Scale"

    def __str__(self) -> str:
        out = super().__str__()
        if self.units is not None:
            out += f" units={self.units}"
        return out


class BasicParameterDefinition(ParameterDefinition):
    """A parameter.

    Most XSPEC parameters use this.

    .. versionchanged:: 4.18.0
       The norm parameter is no-longer included for additive models.

    """

    paramtype = "Basic"

    def __init__(self, name: str, default: float, units: str | None,
                 *,
                 softmin: float, softmax: float,
                 hardmin: float | None, hardmax: float | None,
                 delta: float) -> None:

        self.name = name

        self.units = units
        self.softmin = softmin
        self.softmax = softmax

        # What to do with hard limits?
        #
        if hardmin is None:
            raise ValueError(f"{name} - missing hardmin")
        if hardmax is None:
            raise ValueError(f"{name} - missing hardmax")

        self.hardmin = hardmin
        self.hardmax = hardmax

        if default < self.softmin:
            self.default = softmin
        elif default > self.softmax:
            self.default = softmax
        else:
            self.default = default

        if delta < 0.0:
            self.frozen = True
            self.delta = abs(delta)
        else:
            self.frozen = False
            self.delta = delta

    def __str__(self) -> str:
        out = f"{self.name} = {self.default} ({self.softmin} to {self.softmax})"
        if self.units is not None:
            out += f" units={self.units}"
        if self.frozen:
            out += " frozen"
        return out

    def param_string(self) -> str:

        out = f"XSParameter(name, '{self.name}', {self.default}, "
        out += f"min={self.softmin}, max={self.softmax}, "
        out += f"hard_min={self.hardmin}, hard_max={self.hardmax}"
        if self.frozen:
            out += ", frozen=True"
        if self.units is not None:
            out += f", units='{self.units}'"
        out += ")"
        return out


def read_model_definition(fh,
                          namefunc: Callable[[str], str]
                          ) -> ModelDefinition | None:
    """Parse the next model definition.

    The code attempts to handle the wide variety of model definitions
    found in both the XSPEC model.dat file and in user models but may
    error out in cases that are supported by XSPEC.

    .. versionchanged:: 4.19.0
       The code now supports XSPEC 13.0.0 "grad=xxx" arguments.

    .. versionchanged:: 4.18.0
       Additive models no-longer contain a norm parameter.

    Parameters
    ----------
    fh : file-like
        It should be set to the end of the last model parsed, or the
        start of the file (any leading empty lines are skipped).
    namefunc : callable
        The routine used to convert an XSPEC model name, such as
        "apec", into the Sherpa class name.

    Returns
    -------
    model : ModelDefinition or None
        A representation of the model or None if the end of the
        file has been reached.

    Notes
    -----
    XSPEC additive models do not contain a normalization parameter, so
    one is added for these cases.

    The model will fail if it contains periodic parameters (that is
    parameters that end with a P) as it is unclear how to handle them
    in Sherpa.

    """

    hdrline = ''
    while hdrline == '':
        hdrline = fh.readline()
        if hdrline == '':
            return None

        hdrline = hdrline.strip()

    # The header line, up to XSPEC 13.0.0, was
    #   modelname npars elo ehi funcname modeltype i1 [i2 [initString]]
    # There is now a (currently undocumented) "grad=xxx" option added
    # to the end of the line, and it is unclear what is treated as optional.
    #
    toks = hdrline.split()
    ntoks = len(toks)
    if ntoks < 7 or ntoks > 10:
        raise ValueError("Expected: modelname npars elo ehi funcname "
                         "modeltype i1 [i2 [initString [grad=xxx]]] "
                         f"but sent:\n{hdrline}")

    name = toks[0]
    clname = namefunc(name)
    npars = int(toks[1])
    if npars < 0:
        raise ValueError(f"Number of parameters is {npars}:\n{hdrline}")

    elo = float(toks[2])
    ehi = float(toks[3])
    funcname = toks[4]
    modeltype = toks[5]

    # It is not clear how exactly the "optional" parts are handled, so try
    # to be generic. If the first two arguments are numeric then assume
    # they are the flags.
    #
    flags = []
    grad = None
    initString = None

    if toks[6] in ["0", "1"]:
        flags.append(int(toks[6]))
    else:
        raise ValueError(f"Expected 0 or 1, found {toks[6]} in:\n{hdrline}")

    if ntoks > 7:
        if toks[7] in ["0", "1"]:
            flags.append(int(toks[7]))
        elif toks[7].startswith("grad="):
            grad = toks[7][5:]
        else:
            initString = toks[7]

    if ntoks > 8:
        if toks[8].startswith("grad="):
            if grad is not None:
                raise ValueError(f"multiple grad= arguments in:\n{hdrline}")

            grad = toks[8][5:]
        elif initString is None:
            initString = toks[8]
        else:
            raise ValueError(f"multiple initString arguments in:\n{hdrline}")

    if ntoks > 9:
        if toks[9].startswith("grad="):
            if grad is not None:
                raise ValueError(f"multiple grad= arguments in:\n{hdrline}")

            grad = toks[9][5:]
        elif initString is None:
            initString = toks[9]
        else:
            raise ValueError(f"multiple initString arguments in:\n{hdrline}")

    pars: list[ParameterDefinition] = []
    while len(pars) < npars:
        pline = fh.readline().strip()

        # When using StringIO we don't get an EOF error, instead it
        # returns the empty string.
        if pline == '':
            nmiss = npars - len(pars)
            raise ValueError(f'model={name} missing {nmiss} parameters')

        pars.append(process_parameter_definition(pline, model=name))

    # Need to define this type for mypy, so make it optional
    factory: type[ModelDefinition] | None = None

    if modeltype == "add":
        factory = AddModelDefinition

    elif modeltype == "mul":
        factory = MulModelDefinition

    elif modeltype == "con":
        factory = ConModelDefinition

    elif modeltype == "mix":
        factory = MixModelDefinition

    elif modeltype == "acn":
        factory = AcnModelDefinition

    elif modeltype == "amx":
        factory = AmxModelDefinition

    else:
        raise ValueError(f"Unexpected model type {modeltype} in:\n{hdrline}")

    # Safety check on the parameter names. We do not make this an
    # error because the user can change the Python parameter names
    # (which we have to do for the XSPEC ismabs model).
    #
    ctr = Counter([par.name.lower() for par in pars])
    for pname, count in ctr.items():
        if count == 1:
            continue

        warning(f"model={name} re-uses parameter name {pname}")

    return factory(name, clname, funcname, flags, elo, ehi, pars,
                   initString=initString, grad=grad)


def mpop(array: list[str]) -> float | None:
    """Pop first element from array (converting to float),
    returning None if empty.
    """

    try:
        return float(array.pop(0))
    except IndexError:
        return None


def pop(array: list[str]) -> float:
    """Pop first element from array (converting to float).

    Raises
    ------
    IndexError
        If there is no element to pop.
    """

    return float(array.pop(0))


def process_parameter_definition(pline: str, model: str) -> ParameterDefinition:
    """Process a parameter description.

    Parameters
    ----------
    pline : str
        The parameter definition
    model : str
        The name of the model to which the parameter definition
        belongs, and is only used in error messages.

    Returns
    -------
    param : ParameterDefinition

    Notes
    -----
    Parameter names are automatically converted to support Python
    attribute-name rules (XSPEC has, as of XSPEC 12.11 or so, got
    better about removing such characters but occasionally it is
    needed, and anything goes with user models).

    """

    if pline.endswith("P"):
        raise ValueError("Periodic parameters are unsupported; "
                         f"model={model}:\n{pline}\n")

    toks = pline.split()
    orig_parname = toks.pop(0)

    if orig_parname.startswith('<') and orig_parname.endswith('>'):
        name = orig_parname[1:-1] + "_ave"
    elif orig_parname.startswith(('$', '*')):
        name = orig_parname[1:]
    else:
        name = orig_parname

    name = name.replace('@', 'At')

    # replace foo(bar) with foo_bar
    # (do this before the following, otherwise have foo_bar_)
    #
    if name.endswith(')'):
        lpos = name.rfind('(')
        if lpos != -1:
            name = name[:lpos] + "_" + name[lpos + 1:-1]

    # Replace unsupported characters with '_'. I'd like
    # to use .translate(), but I am too lazy to see how
    # this works.
    valid_chars = string.ascii_letters + string.digits + '_'

    def cconv(c):
        return c if c in valid_chars else '_'

    name = "".join(map(cconv, name))

    if name in ["break", "lambda", "type"]:
        name += "_"

    if orig_parname.startswith('$'):
        # switch parameter
        # the XSPEC documentation say that switches only have 2
        # arguments but the model.dat from it's own model definitions
        # includes these cases:
        #
        # $switch    1     0       0     1      1       -1
        # $method   " "   1       1       1       3       3     -0.01
        # $model    " "     0
        #
        ntoks = len(toks)
        if ntoks == 1:
            idefault = int(toks[0])
            return SwitchParameterDefinition(name, idefault)

        if ntoks == 6:
            idefault = int(toks.pop(0))
            hardmin = float(toks.pop(0))
            softmin = float(toks.pop(0))
            softmax = float(toks.pop(0))
            hardmax = float(toks.pop(0))
            delta   = float(toks.pop(0))
            return SwitchParameterDefinition(name, idefault, None,
                                             softmin=softmin, softmax=softmax,
                                             hardmin=hardmin, hardmax=hardmax,
                                             delta=delta)

        if ntoks > 6:
            # ignore units for now
            delta   = float(toks.pop())
            hardmax = float(toks.pop())
            softmax = float(toks.pop())
            softmin = float(toks.pop())
            hardmin = float(toks.pop())
            idefault = int(toks.pop())
            return SwitchParameterDefinition(name, idefault, None,
                                             softmin=softmin, softmax=softmax,
                                             hardmin=hardmin, hardmax=hardmax,
                                             delta=delta)

        if toks[0].startswith('"'):
            # assume something like '$model " " val'
            # Technically the value should be an int but you can see '1.'
            # in the XSPEC model.dat (HEASARC 6.28)
            # default = int(toks.pop())
            val = toks.pop().removesuffix('.')
            idefault = int(val)
            return SwitchParameterDefinition(name, idefault)

        raise NotImplementedError(f"(switch) model={model} pline=\n{pline}")

    # Handle units
    units: str | None = None

    val = toks.pop(0)
    if val.startswith('"'):
        units = val[1:]
        if units.endswith('"'):
            units = units[:-1]

        else:
            flag = True
            unit_list = [units]
            while flag:
                try:
                    val = toks.pop(0)
                except IndexError as exc:
                    raise ValueError("Unable to parse units; model="
                                     f"{model}\n{pline}") from exc

                if val.endswith('"'):
                    val = val[:-1]
                    flag = False

                unit_list.append(val)

            units = ' '.join(unit_list).strip()

    else:
        units = val

    if units.strip() == '':
        units = None

    if orig_parname.startswith('*'):
        # scale parameter
        default = float(toks.pop(0))

        # Create new variables otherwise mypy doesn't like the fact
        # that these are maybe's.
        #
        s_hardmin = mpop(toks)
        s_softmin = mpop(toks)
        s_softmax = mpop(toks)
        s_hardmax = mpop(toks)
        s_delta   = mpop(toks)

        return ScaleParameterDefinition(name, default, units,
                                        softmin=s_softmin, softmax=s_softmax,
                                        hardmin=s_hardmin, hardmax=s_hardmax,
                                        delta=s_delta)

    if len(toks) != 6:
        raise ValueError(f"Expected 6 values after units; model={model}"
                         f"\n{pline}")

    default = pop(toks)
    hardmin = pop(toks)
    softmin = pop(toks)
    softmax = pop(toks)
    hardmax = pop(toks)
    delta = pop(toks)

    return BasicParameterDefinition(name, default, units,
                                    softmin=softmin, softmax=softmax,
                                    hardmin=hardmin, hardmax=hardmax,
                                    delta=delta)


def add_xs_prefix(inval: str) -> str:
    """Returns XS prepended to inval"""
    return f"XS{inval}"


def parse_xspec_model_description(modelfile,
                                  namefunc: Callable[[str], str] = add_xs_prefix
                                  ) -> list[ModelDefinition]:
    """Given an XSPEC model file - e.g. the lmodel.dat file -
    return information about the models it contains.

    Parameters
    ----------
    modelfile : str or os.PathLike or file-like
        The name of the model file (often called model.dat or
        lmodel.dat) or a file-like object containing the file
    namefunc : callable, optional
        The routine used to convert an XSPEC model name, such as
        "apec", into the Sherpa class name. The default function
        prepends 'XS' to the name.

    Returns
    -------
    models : list of ModelDefinition
        A representation of each model. This will include models that
        Sherpa does not support at this time (e.g. mixing models).

    Raises
    ------
    ValueError
        An invalid or unsupported parameter line, or an unrecognized
        model type, was found.

    """

    emsg = 'namefunc must be a callable which takes and returns a string'
    try:
        ans = namefunc('x')
    except TypeError:
        raise ValueError(emsg) from None

    if not isinstance(ans, str):
        raise ValueError(emsg)

    def process_fh(fh):
        out = []
        while True:
            # If there is a problem reading in a model definition then
            # we do not try to recover - e.g. by wrapping this in a
            # try/except block - since it is not clear how to skip over
            # the "invalid" model definition so that we can move to the
            # next model (well, some simple heuristics could be applied,
            # but leave off developing these until it turns out to be
            # a problem).
            #
            # A simple option would be to just stop parsing as soon as
            # there is a problem, but process any parsed model.
            #
            mdl = read_model_definition(fh, namefunc=namefunc)
            if mdl is None:
                break

            out.append(mdl)

        return out

    # Check if we have a StringIO instance
    #
    if hasattr(modelfile, 'read'):
        with modelfile as fh:
            out = process_fh(fh)
    else:
        with open(modelfile, "r") as fh:
            out = process_fh(fh)

    return out


def simple_wrap(modelname: str,
                mdl: ModelDefinition,
                internal: Version | None = None
                ) -> str:
    """Create the Python class wrapping this model.

    This creates the "starting point" for the user (it can be used
    without further work but the documentation will be poor).

    Parameters
    ----------
    modelname : str
        The XSPEC parent model class (without the leading 'XS').
    mdl : ModelDefinition
        The model.
    internal
        Is this for sherpa.astro.xspec?

    Returns
    -------
    mtext : str
        The model class.

    """

    t1 = ' ' * 4
    t2 = ' ' * 8
    out = "\n"

    if internal is None:
        label = f"XSPEC {modelname}: {mdl.name}"

    else:
        out += f'@version_at_least("{internal[0]}.{internal[1]}.{internal[2]}")\n'

        label = f"The XSPEC {mdl.name}"
        if isinstance(mdl, ConModelDefinition):
            label += " convolution"

        label += " model:  TBD"

    out += f"class {mdl.clname}(XS{modelname}):\n"
    out += f'{t1}"""{label}\n\n'

    if internal is not None:
        # NOTE: the Sherpa version could be guessed at, but better to
        # make it obvious that it needs replacing.
        #
        out += f'{t1}The model is described at [1]_.\n\n'
        out += f'{t1}.. versionadded:: ???\n'
        out += f'{t1}   This model requires XSPEC ' + \
            f'{internal[0]}.{internal[1]}.{internal[2]} or later.\n\n'

    # These are not parameters to __init__ but attributes, so label
    # them as such (using Parameters leads to subtly-different output
    # from sphinx).
    #
    out += f'{t1}Attributes\n'
    out += f'{t1}----------\n'
    for par in mdl.pars:
        out += f'{t1}{par.name}\n'

    # Add in the norm parameter for additive models.
    if isinstance(mdl, AddModelDefinition):
        out += f'{t1}norm\n'

    if internal is not None:
        # This may not be the correct URL, such as the redshift variant, but
        # but it should be close.
        #
        # This code could try to access the URL to note ones where the
        # name is not right.
        cname = mdl.name.capitalize()

        out += '\n'
        out += f'{t1}References\n'
        out += f'{t1}----------\n\n'
        out += f'{t1}.. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodel{cname}.html\n'

    out += f'\n{t1}"""\n\n'

    if internal is None:
        out += f"{t1}_module = _models\n\n"

    out += f"{t1}_xspec_name = '{mdl.name}'\n\n"

    out += f"{t1}def __init__(self, name='{mdl.name}'):\n"
    parnames = []
    for par in mdl.pars:
        # Skip norm if an additive model
        if par.name == "norm" and mdl.modeltype == "Add":
            continue

        out += f"{t2}self.{par.name} = {par.param_string()}\n"
        parnames.append(f"self.{par.name}")

    npars = len(parnames)
    if mdl.modeltype != "Add":
        assert npars > 0, f'Expected at least 1 parameter for {modelname} model'

    if npars == 0:
        pstr = "()"
    elif npars == 1:
        pstr = f"({parnames[0]},)"
    else:
        pstr = f"({', '.join(parnames)})"

    out += "\n"
    if mdl.modeltype == "Add":
        out += f"{t2}# norm parameter is automatically added by XSAdditiveModel\n"
    out += f"{t2}pars = {pstr}\n"
    out += f"{t2}XS{modelname}.__init__(self, name, pars)\n"

    nflags = len(mdl.flags)

    # If the model needs to be recalculated-per-spectrum turn off the
    # caching. This needs to be done after the parent class has been
    # initialized.
    #
    if nflags > 1 and mdl.flags[1] == 1:
        out += f"{t2}self.cache = 0\n"
        # Still warn the user that this is not tested.
        out += f"{t2}warnings.warn('support for models like xs{mdl.name.lower()} "
        out += "(recalculated per spectrum) is untested.')\n"

    # warn about untested models?
    #
    if nflags > 0 and mdl.flags[0] == 1:
        out += f"{t2}warnings.warn('support for models like xs{mdl.name.lower()} "
        out += "(variances are calculated by the model) is untested.')\n"

    out += "\n"
    return out


def additive_wrap(mdl: ModelDefinition,
                  internal: Version | None = None) -> str:
    """Return a string representing the Python code used to wrap
    up access to an Additive user model.
    """

    return simple_wrap('AdditiveModel', mdl, internal=internal)


def multiplicative_wrap(mdl: ModelDefinition,
                        internal: Version | None = None
                        ) -> str:
    """Return a string representing the Python code used to wrap
    up access to an Multiplicative user model.
    """

    return simple_wrap('MultiplicativeModel', mdl, internal=internal)


def convolution_wrap(mdl: ModelDefinition,
                     internal: Version | None = None
                     ) -> str:
    """Return a string representing the Python code used to wrap
    up access to a Convolution user model.
    """

    return simple_wrap('ConvolutionKernel', mdl, internal=internal)


def model_to_python(mdl: ModelDefinition,
                    internal: Version | None = None) -> str:
    """Return a string representing the Python code used to wrap
    up access to the given user model.

    Parameters
    ----------
    mdl : ModelDefinition
    internal : optional
       If set then this is for sherpa.astro.xspec.

    Returns
    -------
    mtext : str
        The model class definition.

    Raises
    ------
    ValueError
        The model is unsupported by Sherpa.

    """

    if mdl.modeltype == "Add":
        return additive_wrap(mdl, internal=internal)

    elif mdl.modeltype == "Mul":
        return multiplicative_wrap(mdl, internal=internal)

    elif mdl.modeltype == "Con":
        return convolution_wrap(mdl, internal=internal)

    else:
        raise ValueError(f"No wrapper for model={mdl.name} "
                         f"type={mdl.modeltype}")


def model_to_compiled(mdl: ModelDefinition) -> tuple[str, str]:
    """Return a string representing the C++ code needed to build the module.

    .. versionchanged:: 4.19.0
       The wrapcode has been updated to include the model name as well
       as the function name and number of parameters, as the library
       routines are now named to match the model name rather than the
       function name.

    Parameters
    ----------
    mdl : ModelDefinition

    Returns
    -------
    wrapcode, defcode : tuple of str, str
        The code needed to build the Python wrapper and any
        definition code needed (the latter can be the empty string).

    Raises
    ------
    ValueError
        The model is unsupported by Sherpa.

    """

    is_fortran = mdl.language.startswith('Fortran')

    # The wrapper code (the Python-accessible function to call this
    # model).
    #
    wrapcode = '  XSPECMODELFCT'
    if mdl.modeltype == "Con":
        wrapcode += '_CON'
        # only have to deal with F77 or not (may need to update)
        if mdl.language == 'Fortran - single precision':
            wrapcode += '_F77'

    elif mdl.modeltype == "Add":
        if not is_fortran:
            wrapcode += '_C'

    elif mdl.modeltype == "Mul":
        # Do we have any double-precision C/C++ models to worry about?
        if is_fortran:
            if mdl.language == 'Fortran - double precision':
                wrapcode += '_DBL'
        else:
            wrapcode += '_C'

    else:
        # This should have been raised by model_to_python
        raise ValueError("Unsupported model")

    funcname = mdl.funcname
    if mdl.language == 'C++ style':
        funcname = f'C_{funcname}'

    # Add in information about the parameters (the number and the
    # names as a single string).
    #
    pnames = ' '.join([p.name for p in mdl.pars])
    wrapcode += f'({mdl.name}, {funcname}, {len(mdl.pars)}, "{pnames}"),'

    # Do we need to define this model? Originally this was only
    # for FORTRAN routines but it may be worth just always
    # declaring it.
    #
    defcode = ''
    if is_fortran:
        defcode = '  xs'
        if mdl.language == 'Fortran - single precision':
            defcode += "f77"
        elif mdl.language == 'Fortran - double precision':
            defcode += "F77"
        else:
            raise RuntimeError(f"Internal error: {mdl.language}")

        defcode += f"Call {mdl.funcname}_;"

    elif mdl.language == "C++ style":
        # Fake up the C++ wrapper as this does not seem to be done for
        # us (not 100% sure about this but it seems to be necessary).
        #
        defcode = f'  XSCCall {mdl.funcname};\n'
        defcode += f'  void C_{mdl.funcname}'
        defcode += '(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr) {\n'
        defcode += f'    const size_t nPar = {len(mdl.pars)};\n'
        defcode += f'    cppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError, initStr, nPar, {mdl.funcname});\n'
        defcode += '  }'

    elif mdl.language == "C style":
        defcode = f"  xsccCall {mdl.funcname};"

    else:
        raise RuntimeError(f"Internal error: {mdl.language}")

    return wrapcode, defcode


def models_to_compiled(mdls: list[ModelDefinition],
                       name: str = "_models") -> str:
    """Return the C++ code needed to build the module.

    Parameters
    ----------
    mdls : list of ModelDefinition
    name : str, optional
        The name of the source / compiled model

    Returns
    -------
    mcode : str
        The wrapper code.

    Raises
    ------
    ValueError
        The model is unsupported by Sherpa.

    Notes
    -----
    Comments are added before each section to make it easier to
    identify (if post processing is needed). The sections are

        // Includes
        // Defines
        // Wrapper
        // Module

    """

    defcode_list = []
    wrapcode_list = []
    has_cxx = False
    for mdl in mdls:
        w, d = model_to_compiled(mdl)

        wrapcode_list.append(w)
        if d != '':
            defcode_list.append(d)

        has_cxx |= (mdl.language == "C++ style")

    defcode = '\n'.join(defcode_list)
    wrapcode = '\n'.join(wrapcode_list)

    def marker(label):
        # Ensure we have a consistent form for these markers
        return f"// {label}\n\n"

    # What includes are needed?
    #
    out = marker("Includes")

    # The Sherpa extension includes which will include Python.h
    # so must be done before other includes, such as iostream.
    #
    out += '#include "sherpa/astro/xspec_extension.hh"\n\n'

    out += '#include <iostream>\n\n'
    out += '#include <xsTypes.h>\n'
    out += '#include <XSFunctions/Utilities/funcType.h>\n\n'

    # The Sherpa/XSPEC interface uses a number of defines to control
    # behavior. These should not be needed for user models, but set
    # them up. Note that they depend on the available XSPEC library,
    # which means that this can only be run when XSPEC support is
    # present (and the output will depend on the XSPEC model library
    # in use).
    #
    from sherpa.astro import xspec
    versionstr = xspec.get_xsversion()
    xspec_version = get_version(versionstr)

    for version in SUPPORTED_VERSIONS:
        if xspec_version >= version:
            major, minor, micro = version
            out += f'#define XSPEC_{major}_{minor}_{micro}\n'

    out += '\n'
    out += marker("Defines")

    # Do we need to define cppModelWrapper? For XSPEC 12.12.1/12.13.0
    # we have to.
    #
    if has_cxx:
        out += 'void cppModelWrapper(const double* energy, int nFlux, const double* params,\n'
        out += '  int spectrumNumber, double* flux, double* fluxError, const char* initStr,\n'
        out += '  int nPar, void (*cppFunc)(const RealArray&, const RealArray&,\n'
        out += '  int, RealArray&, RealArray&, const string&));\n'
        out += '\n'

    out += 'extern "C" {\n'
    out += f'{defcode}\n'
    out += '}\n\n'

    out += marker("Wrapper")
    out += 'static PyMethodDef Wrappers[] = {\n'
    out += f'{wrapcode}\n'
    out += '  { NULL, NULL, 0, NULL }\n'
    out += '};\n\n'

    # Now the Python module
    #
    out += marker("Module")
    out += 'static struct PyModuleDef wrapper_module = {\n'
    out += '  PyModuleDef_HEAD_INIT,\n'
    out += f'  "{name}",\n'
    out += '  NULL,\n'
    out += '  -1,\n'
    out += '  Wrappers,\n'
    out += '};\n\n'

    out += f'PyMODINIT_FUNC PyInit_{name}(void) {{\n'
    out += '  import_array();\n'
    out += '  return PyModule_Create(&wrapper_module);\n'
    out += '}\n'

    return out


def create_xspec_code(models: list[ModelDefinition],
                      name: str = "_models",
                      internal: Version | None = None
                      ) -> XSPECcode:
    """Create the Python classes and C++ code for the models.

    Create the code fragments needed to build the XSPEC interface
    but they are not complete.

    Parameters
    ----------
    models : list of ModelDefiniton
    name : str, optional
        The name of the module.
    internal : optional
        If this is for sherpa.astro.xspec then what version of XSPEC
        is it for? If used it assumes the model is new to this version.

    Returns
    -------
    code : XSPECcode
        The code is accessible as the 'python' and 'compiled' fields.

    Notes
    -----

    We skip any model functions that are used in multiple models, as
    this was an error in the XSPEC 12.8.2 model.dat which caused the
    eplogpar to call the wrong function. This has been fixed but we
    add a check here just in case.

    """

    ctr = Counter([mdl.funcname for mdl in models])
    invalidnames = [n for n, c in ctr.items() if c > 1]
    if len(invalidnames) > 0:
        newmodels = []
        for mdl in models:
            if mdl.funcname not in invalidnames:
                newmodels.append(mdl)
                continue

            warning(f"Skipping model {mdl.name} as it calls " +
                    f"{mdl.funcname} which is used by " +
                    f"{ctr[mdl.funcname]} different models")

        models = newmodels
        del newmodels

    # Strip out unsupported models
    #
    mdls = []
    langs = set()
    requires_warnings = False
    for mdl in models:
        if mdl.modeltype in ['Mix', 'Acn']:
            warning(f"Skipping {mdl.name} as model type = {mdl.modeltype}")
            continue

        # The following check should never fire, but leave in
        if mdl.language not in ['Fortran - single precision',
                                'Fortran - double precision',  # un-tested
                                'C style', 'C++ style']:
            warning(f"Skipping {mdl.name} as language = {mdl.language}")
            continue

        nflags = len(mdl.flags)
        if nflags > 0:
            if mdl.flags[0] == 1:
                warning(f"{mdl.name} calculates model variances; this is untested/unsupported in Sherpa")
                requires_warnings = True

            if nflags > 1 and mdl.flags[1] == 1:
                warning(f"{mdl.name} needs to be re-calculated per spectrum; this is untested.")
                requires_warnings = True

        langs.add(mdl.language)
        mdls.append(mdl)

    nmdl = len(mdls)
    if nmdl == 0:
        raise ValueError("No supported models were found!")

    if requires_warnings:
        python = "import warnings\n"
    else:
        python = ""

    python += "\n\n".join([model_to_python(mdl, internal=internal) for mdl in mdls])
    compiled = models_to_compiled(mdls, name=name)
    return XSPECcode(python=python, compiled=compiled)


def find_supported_xspec_models(enabled: bool = True
                                ) -> list[str]:
    """The names of the supported XSPEC models.

    This routine requires that XSPEC support is enabled and the return
    value depends on the version of the XSPEC library in use.

    .. versionadded: 4.19.9

    Parameters
    ----------
    enabled : bool, optional
       Should the list be filtered to only those models that can be
       run with the provided XSPEC library?

    Returns
    -------
    models : list of str
       The list of the supported model names. Note that the return
       value matches the XSPEC model.dat (e.g. "TBabs" for the XSTBabs
       model).

    """

    from sherpa.astro import xspec as xs

    def is_proper_subclass(obj, cls):
        if obj in cls:
            return False
        return issubclass(obj, cls)

    valid_names: set[str] = set()
    for clname in dir(xs):
        if not clname.startswith("XS"):
            continue

        cls = getattr(xs, clname)
        if is_proper_subclass(cls, (xs.XSAdditiveModel,
                                    xs.XSMultiplicativeModel,
                                    xs.XSConvolutionKernel)):
            if enabled and not cls.version_enabled:
                continue

            valid_names.add(cls._xspec_name)

    return sorted(valid_names)


@dataclass
class XSPECDataRow:
    """A mapping from model to a data file it requires.

    This is read from the modelDataFiles.csv used by XSPEC 13.0.0
    and later, but split by model.

    .. versionadded: 4.19.9

    """

    filename: str
    """The name of the data file.

    The data file is expected to be found in the location returned
    by `sherpa.astro.xspec.get_xspath_model()`.
    """

    model: str | None
    """The model that uses the file, if any.

    This name matches that given in the file returned by
    `get_model_data_version_path()`, and not that used by Sherpa
    (e.g. "apec" rather than "XSapec" or "xsapec").
    """

    model_version: str | None
    """The version of the data file."""

    xspec_version_start: Version
    """The start of the version range (this drops any patch level)."""

    xspec_version_end: Version | None
    """The end of the version range (this drops any patch level), af any.

    A value of "Unknown" is converted to `None`.
    """

    release: str
    """The release number (expected to be "1" for current and "0" otherwise)."""

    checksum: str
    """The checksum for the file."""

    filesize: int
    """The file size."""


def get_model_data_version_path() -> Path | None:
    """Return the location of the XSPEC model data version file.

    This is only expected to return a value with XSPEC 13.0.0 or
    later.

    .. versionadded: 4.19.9

    Returns
    -------
    path : Path or None
       The location of the CSV file (if it exists) or None.

    See Also
    --------
    read_model_data_version_file

    """

    from sherpa.astro import xspec

    csvfile = Path(xspec.get_xspath_manager()) / "modelDataFiles.csv"
    debug("CSV location: %s", csvfile)
    if not csvfile.is_file():
        return None

    return csvfile


def read_model_data_version_file(latest: bool = True,
                                 version: str | None = None,
                                 models: Sequence[str] | None = None
                                 ) -> list[XSPECDataRow] | None:
    """Return the XSPEC model data file information, if present.

    This is only expected to return information with XSPEC 13.0.0
    and later.

    .. versionadded: 4.19.9

    Parameters
    ----------
    latest : bool, optional
       Should the search be restricted to the latest data files?
    version : str or None, optional
       What version filter should be applied? If None then all versions,
       otherwise limit to just this version (which is expected to be
       a value like '12.13.0' or '12.14.1e'). This parameter is only
       used if latest=False.
    models : sequence of str or None, optional
       If not None, then this is a list of model names to restrict the
       search to. Thes name use the values from the XSPEC `model.dat`
       file (converted to lower case) and not the Sherpa name.

    Returns
    -------
    rows : list of XSpecDataRow or None
       The data rows, split up by model, so this differs from how the
       data is stored in the CSV file, and it means that the same file
       can appear in multiple rows. Only models supported by Sherpa
       are included.

    See Also
    --------
    get_model_data_version_path

    Examples
    --------

    Find the rows that refer to the latest version of the data files:

    >>> latest_rows = read_model_data_version_file()

    What rows relate to the XSPEC apec model (note that this uses the
    name from the `model.dat` file):

    >>> apec_rows = read_model_data_version_file(models=["apec"])

    Manually select the apec rows:

    >>> latest_rows = read_model_data_version_file()
    >>> apec_rows = [r for r in latest_rows if r.model == "apec"]

    Find all rows:

    >>> all_rows = read_model_data_version_file(latest=False)

    Find the rows for XSPEC version 12.15.1:

    >>> rows = read_model_data_version_file(latest=False, version='12.15.1')

    """

    import csv

    csvfile = get_model_data_version_path()
    if csvfile is None:
        debug("no csv file found!")
        return None

    # Is a version check needed? This is allowed even if latest=True.
    vcheck : Version | None
    if version is None:
        vcheck = None
    else:
        vcheck = get_version(version)

    # What XSPEC model names do we need to care about? The check uses
    # the lower-case version of the model name from model.dat.
    #
    valid_names = {m.lower() for m in find_supported_xspec_models()}

    if models is None:
        req_models = None
    else:
        req_models = [m.lower() for m in models]

    def mk(rec: dict[str, str]) -> list[XSPECDataRow]:
        """Create zero, one, or more values from the row."""

        # If only the latest versions are being used then skip
        # old models.
        #
        if latest and rec["release"] != "1":
            debug(" - not the latest release")
            return []

        # Convert the input data to the arguments for XSPECDataRow
        #
        kwargs = {}

        # Explicit version clean up.
        # - convert Unknown to None
        # - drop trailing patch level
        #
        vstart = rec['xspec_version_start']
        v1 = get_version(vstart)
        if vcheck is not None and v1 > vcheck:
            debug("- too new %s", vstart)
            return []

        vend = rec['xspec_version_end']
        if vend in ["", "Unknown"]:
            v2 = None
        else:
            v2 = get_version(vend)

        if vcheck is not None and v2 is not None and v2 < vcheck:
            debug("- too old %s", vend)
            return []

        kwargs['xspec_version_start'] = v1
        kwargs['xspec_version_end'] = v2
        kwargs['filesize'] = int(rec['filesize'])

        # Copy over the remaining arguments.
        #
        for name in ["filename", "release", "checksum"]:
            kwargs[name] = rec[name]

        for name in ["model_version"]:
            kwargs[name] = rec[name] if rec[name] != '' else None

        # models can be
        # - empty
        # - a single name
        # - multiple models stored as "mdl1, ..., mdl"
        #
        rmodels = rec['models']
        if rmodels == '':
            # Should the unlabelled models be included or not?
            if req_models is not None:
                return []

            mnames = [None]
        else:
            all_mnames = [m.strip() for m in rmodels.split(",")]
            debug("models: %s", str(all_mnames))

            # Filter by those models that Sherpa supports.
            #
            mnames = [m for m in all_mnames if m.lower() in valid_names]
            nall = len(all_mnames)
            nmdl = len(mnames)
            if nmdl < nall:
                debug("- removed %d unsupported models", nall - nmdl)

            # Filter by requested user model
            if req_models is not None:
                mnames = [m for m in mnames if m.lower() in req_models]

        return [XSPECDataRow(model=m, **kwargs) for m in mnames]

    # Store those models that have missing data. The assumption is that
    # the CSV has columns for
    #    filename
    #    models
    #    model_version
    #    xspec_version_start
    #    xspec_version_end
    #    release
    #    checksum
    #    filesize
    #
    # What does a version like "12.14.0x" mean (since there was no
    # patch level x release for 12.14.0). It looks like ftgetmodeldata
    # just skips the patch-level value so we will do that too.
    #
    out = []
    with csvfile.open(mode='rt') as fh:
        reader = csv.reader(fh)
        header = None

        for row in reader:
            if header is None:
                header = row
                debug("header: %s", header)
                continue

            toks = dict(zip(header, row))
            infile = toks['filename']
            debug("file: %s", infile)

            # Convert to zero to n records.
            #
            out.extend(mk(toks))

    if len(out) == 0:
        # This is not expected to happen, but just in case
        debug("no data found in csv file!")
        return None

    return out


def find_missing_model_data_files(rows: Sequence[XSPECDataRow]
                                  ) -> list[XSPECDataRow] | None:
    """What XSPEC models are missing data?

    Filters the input list to remove those rows where the data
    file exists and has the correct size.

    .. versionadded:: 4.19.0

    Parameters
    ----------
    rows : sequence of XSPECDataRow
       The list of rows to check.

    Returns
    -------
    missing : list of XSPECDataRow or None
       The rows that represent missing data files, or None if no files
       are missing. Only the first appearance of a missing file is
       returned.

    See Also
    --------
    read_model_data_version_file

    Examples
    --------

    >>> latest_rows = read_model_data_version_file()
    >>> missing_rows = find_missing_model_data_files(latest_rows)

    """

    from sherpa.astro.xspec import get_xspath_model

    mpath = Path(get_xspath_model())
    out = []
    seen = set()
    for row in rows:

        # If we have seen this before (from another model) then
        # skip.
        #
        if row.filename in seen:
            continue

        # Check that the model has the expected size.
        modpath = mpath / row.filename
        with suppress(OSError):
            if modpath.stat().st_size == row.filesize:
                seen.add(row.filename)
                continue

        out.append(row)

    if len(out) == 0:
        return None

    return out


def download_missing_model_data_files(rows: Sequence[XSPECDataRow],
                                      remotedir: str | None = None
                                      ) -> list[XSPECDataRow] | None:
    """Download the data files for the XSPEC models.

    This skips any row where the files are already downloaded. This
    is similar to the `ftgetmodeldata
    <https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/help/ftgetmodeldata.html>`_
    tool from HEASARC.

    .. versionadded:: 4.19.0

    Parameters
    ----------
    rows : sequence of XSPECDataRow
       The list of rows to download.
    remotedir : str or None, optional
       The URL where the data files are stored. The default uses the
       HEASARC location.

    Returns
    -------
    failed : list of XSPECDataRow or None
       If not None then this represents those files that could not
       be downloaded.

    See Also
    --------
    find_missing_model_data_files, read_model_data_version_file

    Notes
    -----
    The data files are stored in the directory returned by
    `sherpa.astro.xspec.get_xspath_model()`. There is no support for
    partial downloads as existing files are over-written if they do
    not have the expected file size.

    Examples
    --------

    >>> latest_rows = read_model_data_version_file()
    >>> download_missing_model_data_files(latest_rows)
    ...

    """

    from http import client
    from shutil import copyfileobj
    from urllib.request import Request

    from sherpa._version import short_version
    from sherpa.astro.xspec import get_xspath_model

    # Exit early if nothing to do.
    if len(rows) == 0:
        return

    if remotedir is None:
        base_url = "https://heasarc.gsfc.nasa.gov/FTP/software/xspec/spectral/modelData"
    else:
        # Assume this is a URL
        base_url = remotedir

    # It is likely to be https, but allow for http.
    #
    match Request(base_url).type:
        case "http":
            Connection = client.HTTPConnection

        case "https":
            Connection = client.HTTPSConnection

        case _:
            raise ValueError(f"Unexpected remotedir URL: {base_url}")

    # Allow tracking of these requests
    #
    hdr = {"User-Agent": f"sherpa/{short_version} ftgetmodeldata equivalent"}

    mpath = Path(get_xspath_model())
    mpath_exists = mpath.is_dir()
    info("Data directory: %s", str(mpath))
    errs = []
    seen = set()
    for row in rows:

        # If we have seen this before (from another model) then
        # skip.
        #
        if row.filename in seen:
            continue

        # If the file exists and has the correct size then do nothing.
        #
        out = mpath / row.filename
        debug("Checking %s", str(out))
        with suppress(OSError):
            got = out.stat().st_size
            if got == row.filesize:
                seen.add(row.filename)
                continue

            # If the unlink fails then it will raise an OSError which
            # will be suppressed.
            #
            debug("- file size mis-match: %s vs %d", got, row.filesize)
            out.unlink()

        # Does the outtput directory exist? If this fails then bail out
        # immediately.
        #
        if not mpath_exists:
            debug("Creating output directory [%s]", str(mpath))
            mpath.mkdir()
            mpath_exists = True

        # Use http.client rather than urllib.request to allow for
        # streaming the response.
        #
        url = f"{base_url}/{row.filename}"
        info("Downloading %s %d", url, row.filesize)
        url_req = Request(url)

        stime = time.localtime()
        conn = Connection(url_req.host)
        conn.request("GET", url_req.selector, headers=hdr)
        resp = conn.getresponse()
        with out.open('wb') as outfh:
            debug("- copying %d bytes", row.filesize)
            copyfileobj(resp, outfh)

        etime = time.localtime()
        debug("- download time %.1f s", time.mktime(etime) - time.mktime(stime))

        # Check the download size.
        #
        with suppress(OSError):
            got = out.stat().st_size
            if got == row.filesize:
                continue

            debug("- file size mis-match: %d vs %d", got, row.filesize)

        errs.append(row)

    if len(errs) == 0:
        return None

    return errs
