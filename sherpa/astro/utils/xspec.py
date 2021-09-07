#
#  Copyright (C) 2012, 2013, 2014, 2015, 2016, 2020, 2021
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

References
----------

.. [1] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/index.html

.. [2] http://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html

.. [3] https://sherpa.readthedocs.io/en/latest/developer/index.html#update-the-xspec-bindings

"""

from collections import Counter, namedtuple
import logging
import string


__all__ = ("parse_xspec_model_description", "create_xspec_code")


warning = logging.getLogger(__name__).warning


class ModelDefinition():
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

    See Also
    --------
    AddModelDefinition, MulModelDefinition, ConModelDefinition,
    MixModelDefintion, AcnModelDefinition, AmxModelDefinition

    Notes
    -----
    Do not instantiate this class directly.

    """

    modeltype = None
    language = None

    def __init__(self, name, clname, funcname, flags, elo, ehi, pars,
                 initString=None):
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
        if self.funcname.startswith('F_'):
            self.language = 'Fortran - double precision'
            self.funcname = self.funcname[2:]
        elif self.funcname.startswith('c_'):
            self.language = 'C style'
            self.funcname = self.funcname[2:]
        elif self.funcname.startswith('C_'):
            self.language = 'C++ style'
            self.funcname = self.funcname[2:]
        else:
            self.language = 'Fortran - single precision'

        if initString is not None and self.language.startswith('F'):
            initString = None

        self.initString = initString

    def __repr__(self):
        return f"<{self.modeltype}:{self.name}:{self.funcname}:{len(self.pars)} pars>"

    def __str__(self):
        pars = "\n".join([str(p) for p in self.pars])
        return f"{self.modeltype}.{self.name} " +  \
            f"function={self.funcname}\n{self.language}\n{pars}"


class AddModelDefinition(ModelDefinition):
    """XSPEC additive models.

    See [1]_ for examples.

    References
    ----------

    .. [1] http://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Additive.html
    """

    modeltype = "Add"


class MulModelDefinition(ModelDefinition):
    """XSPEC multiplicative models.

    See [1]_ for examples.

    References
    ----------

    .. [1] http://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Multiplicative.html
    """

    modeltype = "Mul"


class ConModelDefinition(ModelDefinition):
    """XSPEC convolution models.

    See [1]_ for examples.

    References
    ----------

    .. [1] http://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Convolution.html
    """

    modeltype = "Con"


class MixModelDefinition(ModelDefinition):
    """XSPEC mixing models.

    See [1]_ for examples. These are currently unsupported in Sherpa.

    References
    ----------

    .. [1] http://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Mixing.html
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


class ParameterDefinition():
    """Represent an XSPEC parameter.

    Parameters
    ----------
    name : str
        The parameter name.
    default : float
        The default value
    units : str or None, optinal
        The unit field. There is no check this meets any standard.
    softmin, softmax, hardmin, hardmax : float or None
        The minimum and maximum values for the parameter (using the
        XSPEC definition of soft and hard, not Sherpa).
    delta : floar or None, optional
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

    paramtype = None

    def __init__(self, name, default, units=None,
                 softmin=None, softmax=None,
                 hardmin=None, hardmax=None, delta=None):
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

    def __str__(self):
        return f"{self.name} = {self.default}"

    def param_string(self):
        out = f"Parameter(name, '{self.name}', {self.default}"

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

    def __str__(self):
        out = super().__str__()
        if self.units is not None:
            out += " units={}".format(self.units)
        return out


class BasicParameterDefinition(ParameterDefinition):
    """A parameter.

    Most XSPEC parameters use this.

    """

    paramtype = "Basic"

    def __init__(self, name, default, units, softmin, softmax,
                 hardmin, hardmax, delta):

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

    def __str__(self):
        out = f"{self.name} = {self.default} ({self.softmin} to {self.softmax})"
        if self.units is not None:
            out += f" units={self.units}"
        if self.frozen:
            out += " frozen"
        return out

    def param_string(self):
        out = f"Parameter(name, '{self.name}', {self.default}, "
        out += f"min={self.softmin}, max={self.softmax}, "
        out += f"hard_min={self.hardmin}, hard_max={self.hardmax}"
        if self.frozen:
            out += ", frozen=True"
        if self.units is not None:
            out += f", units='{self.units}'"
        out += ")"
        return out


def read_model_definition(fh, namefunc):
    """Parse the next model definition.

    The code attempts to handle the wide variety of model definitions
    found in both the XSPEC model.dat file and in user models but may
    error out in cases that are supported by XSPEC.

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

    toks = hdrline.split()
    ntoks = len(toks)
    if ntoks < 7 or ntoks > 9:
        raise ValueError("Expected: modelname npars elo ehi funcname modeltype i1 [i2 [initString]] but sent:\n{}".format(hdrline))

    name = toks[0]
    clname = namefunc(name)
    npars = int(toks[1])
    if npars < 0:
        raise ValueError("Number of parameters is {}:\n{}".format(npars, hdrline))

    elo = float(toks[2])
    ehi = float(toks[3])
    funcname = toks[4]
    modeltype = toks[5]

    if ntoks == 9:
        initString = toks.pop()
    else:
        initString = None

    flags = [int(t) for t in toks[6:]]

    pars = []
    while len(pars) < npars:
        pline = fh.readline().strip()

        # When using StringIO we don't get an EOF error, instead it
        # returns the empty string.
        if pline == '':
            nmiss = npars - len(pars)
            raise ValueError(f'model={name} missing {nmiss} parameters')

        pars.append(process_parameter_definition(pline, model=name))

    if modeltype == "add":
        nstr = 'norm " " 1.0 0.0 0.0 1.0e24 1.0e24 0.1'
        pars.append(process_parameter_definition(nstr, model=name))
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
        raise ValueError("Unexpected model type {} in:\n{}".format(modeltype,
                                                                   hdrline))

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
                   initString=initString)


def mpop(array, defval=None):
    """Pop first element from array (converting to float),
    returning defval if empty.
    """

    try:
        return float(array.pop(0))
    except IndexError:
        return defval


def pop(array):
    """Pop first element from array (converting to float).

    Raises
    ------
    IndexError
        If there is no element to pop.
    """

    return float(array.pop(0))


def process_parameter_definition(pline, model):
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
        raise ValueError("Periodic parameters are unsupported; model={}:\n{}\n".format(model, pline))

    toks = pline.split()
    orig_parname = toks.pop(0)

    if orig_parname.startswith('<') and orig_parname.endswith('>'):
        name = orig_parname[1:-1] + "_ave"
    elif orig_parname.startswith('$') or orig_parname.startswith('*'):
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
            default = int(toks[0])
            return SwitchParameterDefinition(name, default)

        if ntoks == 6:
            default = int(toks.pop(0))
            hardmin = float(toks.pop(0))
            softmin = float(toks.pop(0))
            softmax = float(toks.pop(0))
            hardmax = float(toks.pop(0))
            delta   = float(toks.pop(0))
            return SwitchParameterDefinition(name, default, None,
                                             softmin, softmax,
                                             hardmin, hardmax, delta)

        if ntoks > 6:
            # ignore units for now
            delta   = float(toks.pop())
            hardmax = float(toks.pop())
            softmax = float(toks.pop())
            softmin = float(toks.pop())
            hardmin = float(toks.pop())
            default = int(toks.pop())
            return SwitchParameterDefinition(name, default, None,
                                             softmin, softmax,
                                             hardmin, hardmax, delta)

        if toks[0].startswith('"'):
            # assume something like '$model " " val'
            # Technically the value should be an int but you can see '1.'
            # in the XSPEC model.dat (HEASARC 6.28)
            # default = int(toks.pop())
            val = toks.pop()
            if val.endswith('.'):
                val = val[:-1]
            default = int(val)
            return SwitchParameterDefinition(name, default)

        raise NotImplementedError("(switch) model={} pline=\n{}".format(model, pline))

    # Handle units
    val = toks.pop(0)
    if val.startswith('"'):
        units = val[1:]
        if units.endswith('"'):
            units = units[:-1]

        else:
            flag = True
            units = [units]
            while flag:
                try:
                    val = toks.pop(0)
                except IndexError as exc:
                    raise ValueError("Unable to parse units; model={}\n{}".format(model, pline)) from exc

                if val.endswith('"'):
                    val = val[:-1]
                    flag = False

                units.append(val)

            units = ' '.join(units).strip()

    else:
        units = val

    if units.strip() == '':
        units = None

    if orig_parname.startswith('*'):
        # scale parameter
        default = float(toks.pop(0))

        hardmin = mpop(toks)
        softmin = mpop(toks)
        softmax = mpop(toks)
        hardmax = mpop(toks)
        delta   = mpop(toks)

        return ScaleParameterDefinition(name, default, units,
                                        softmin, softmax,
                                        hardmin, hardmax, delta)

    if len(toks) != 6:
        raise ValueError("Expected 6 values after units; model={}\n{}".format(model, pline))

    default = pop(toks)
    hardmin = pop(toks)
    softmin = pop(toks)
    softmax = pop(toks)
    hardmax = pop(toks)
    delta = pop(toks)

    return BasicParameterDefinition(name, default, units,
                                    softmin, softmax,
                                    hardmin, hardmax, delta)


def add_xs_prefix(inval):
    """Returns XS prepended to inval"""
    return f"XS{inval}"


def parse_xspec_model_description(modelfile, namefunc=add_xs_prefix):
    """Given an XSPEC model file - e.g. the lmodel.dat file -
    return information about the models it contains.

    Parameters
    ----------
    modelfile : str or os.PathLike or file-like
        The name of the model file (often called model.dat or
        lmodel.dat) or a file-like object containing the file
    namefunc : callable or None, optional
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
            # the "invalid" model definiton so that we can move to the
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


def simple_wrap(modelname, mdl):
    """Create the Python class wrapping this model.

    This creates the "starting point" for the user (it can be used
    without further work but the documentation will be poor).

    Parameters
    ----------
    modelname : str
        The XSPEC parent model class (without the leading 'XS').
    mdl : ModelDefinition
        The model.

    Returns
    -------
    mtext : str
        The model class.

    """

    t1 = ' ' * 4
    t2 = ' ' * 8
    out = f"\nclass {mdl.clname}(XS{modelname}):\n"
    out += f'{t1}"""XSPEC {modelname}: {mdl.name}\n\n'
    out += f'{t1}Parameters\n'
    out += f'{t1}----------\n'
    for par in mdl.pars:
        out += f'{t1}{par.name}\n'

    out += f'\n{t1}"""\n'

    if mdl.language == 'C++ style':
        funcname = f"C_{mdl.funcname}"
    else:
        funcname = mdl.funcname

    out += f"{t1}_calc = _models.{funcname}\n"

    out += "\n"
    out += f"{t1}def __init__(self, name='{mdl.name}'):\n"
    parnames = []
    for par in mdl.pars:
        out += f"{t2}self.{par.name} = {par.param_string()}\n"
        parnames.append(f"self.{par.name}")

    assert len(parnames) > 0, f'Expected at least 1 parameter for {modelname} model'
    if len(parnames) == 1:
        pstr = f"({parnames[0]},)"
    else:
        pstr = f"({','.join(parnames)})"

    # warn about untested models?
    nflags = len(mdl.flags)
    if nflags > 0:
        if mdl.flags[0] == 1:
            out += f"{t2}warnings.warn('support for models like {mdl.clname.lower()} "
            out += "(variances are calculated by the model) is untested.')\n"

        if nflags > 1 and mdl.flags[1] == 1:
            out += "{t2}warnings.warn('support for models like {mdl.clname.lower()} "
            out += "(recalculated per spectrum) is untested.')\n"

    out += f"{t2}XS{modelname}.__init__(self, name, {pstr})\n"
    out += "\n"
    return out


def additive_wrap(mdl):
    """Return a string representing the Python code used to wrap
    up access to an Additive user model.
    """

    return simple_wrap('AdditiveModel', mdl)


def multiplicative_wrap(mdl):
    """Return a string representing the Python code used to wrap
    up access to an Multiplicative user model.
    """

    return simple_wrap('MultiplicativeModel', mdl)


def convolution_wrap(mdl):
    """Return a string representing the Python code used to wrap
    up access to a Convolution user model.
    """

    return simple_wrap('ConvolutionKernel', mdl)


def model_to_python(mdl):
    """Return a string representing the Python code used to wrap
    up access to the given user model.

    Parameters
    ----------
    mdl : ModelDefinition

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
        return additive_wrap(mdl)

    elif mdl.modeltype == "Mul":
        return multiplicative_wrap(mdl)

    elif mdl.modeltype == "Con":
        return convolution_wrap(mdl)

    else:
        raise ValueError("No wrapper for model={} type={}".format(mdl.name, mdl.modeltype))


def model_to_compiled(mdl):
    """Return a string representing the C++ code needed to build the module.

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

        wrapcode += '_NORM'

    elif mdl.modeltype == "Mul":
        if not is_fortran:
            wrapcode += '_C'

    else:
        # This should have been raised by model_to_python
        raise ValueError("Unsupported model")

    funcname = mdl.funcname
    if mdl.language == 'C++ style':
        funcname = f'C_{funcname}'

    wrapcode += f'({funcname}, {len(mdl.pars)}),'

    # Do we need to define this model? We assume we can use the same logic
    # as for the XSPEC model library - we only need to define FORTRAN
    # code.
    #
    # Do we need to include C style?
    #
    defcode = ''
    if is_fortran:
        defcode = f'  void {mdl.funcname}'
        if mdl.language == 'Fortran - single precision':
            # Note add _ (assume this is needed)
            defcode += '_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);'
        elif mdl.language == 'Fortran - double precision':
            # Note add _ (assume this is needed)
            defcode += '_(double* ear, int* ne, double* param, int* ifl, double* photar, double* photer);'

    return wrapcode, defcode


def models_to_compiled(mdls):
    """Return the C++ code needed to build the module.

    Parameters
    ----------
    mdls : list of ModelDefinition

    Returns
    -------
    mcode : str
        The wrapper code.

    Raises
    ------
    ValueError
        The model is unsupported by Sherpa.

    """

    defcode = []
    wrapcode = []
    for mdl in mdls:
        w, d = model_to_compiled(mdl)

        wrapcode.append(w)
        if d != '':
            defcode.append(d)

    defcode = '\n'.join(defcode)
    wrapcode = '\n'.join(wrapcode)

    out = '// Defines\n\n'
    out += 'extern "C" {\n'
    out += f'{defcode}\n}}\n\n// Wrapper\n\n'
    out += 'static PyMethodDef Wrappers[] = {\n'
    out += f'{wrapcode}'
    out += '\n  { NULL, NULL, 0, NULL }\n}\n'
    return out


def create_xspec_code(models):
    """Create the Python classes and C++ code for the models.

    Create the code fragments needed to build the XSPEC interface
    but they are not complete.

    Parameters
    ----------
    models : list of ModelDefiniton

    Returns
    -------
    code : namedtuple
        The code stored in a namedtuple with fields 'python' and
        'compiled'.

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

            if nflags > 1 and mdl.flags[1] == 1:
                warning(f"{mdl.name} needs to be re-calculated per spectrum; this is untested.")

        langs.add(mdl.language)
        mdls.append(mdl)

    nmdl = len(mdls)
    if nmdl == 0:
        raise ValueError("No supported models were found!")

    out = namedtuple('XSPECcode', ['python', 'compiled'])

    out.python = "\n\n".join([model_to_python(mdl) for mdl in mdls])
    out.compiled = models_to_compiled(mdls)
    return out
