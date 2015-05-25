# 
#  Copyright (C) 2010  Smithsonian Astrophysical Observatory
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

import string
from sherpa.models import Parameter, ArithmeticModel, modelCacher1d
from sherpa.models.parameter import hugeval
import sherpa.astro.xspec._xspec
from sherpa.utils import guess_amplitude, param_apply_limits
from sherpa.astro.utils import get_xspec_position
from sherpa.astro.xspec._xspec import get_xschatter, get_xsabund, get_xscosmo, \
     get_xsxsect, set_xschatter, set_xsabund, set_xscosmo, set_xsxsect, \
     get_xsversion

# Wrap the XSET function in Python, so that we can keep a record of
# the strings the user sent as specific XSPEC model strings (if any) during
# the session.  Only store if setting was successful.
# See:
# http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSxset.html
modelstrings = {}

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
       empty string, if not set.

    See Also
    --------
    set_xsxset : Set a X-Spec model setting.

    Notes
    -----
    Due to the way X-Spec model settings work, `get_xsxset` will only
    return a value if it has previously been set with a call to
    `set_xsxset`. There is no way to retrive the default value of a
    setting.
    """
    name = name.upper()
    return sherpa.astro.xspec._xspec.get_xsxset(name)

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
    get_xsxset : Return the value of X-Spec model settings.
    get_xsversion : Return the version of X-Spec used by this module.
    set_xsabund : Set the X-Spec abundance table.
    set_xschatter : Set the X-Spec chatter level.
    set_xscosmo : Set the X-Spec cosmology settings.
    set_xsxsect : Set the X-Spec photoelectric absorption cross-sections.

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
    sherpa.astro.xspec._xspec.set_xsxset(name, value)
    if (get_xsxset(name) != ""):
        modelstrings[name] = get_xsxset(name)

# Provide XSPEC module state as a dictionary.  The "cosmo" state is
# a 3-tuple, and "modelstrings" is a dictionary of model strings
# applicable to certain models.  The abund and xsect settings are
# strings.  The chatter setting is an integer.  Please see the
# XSPEC manual concerning the following commands: abund, chatter,
# cosmo, xsect, and xset.
# http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/Control.html
# http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/Setting.html
def get_xsstate():
    return {"abund": get_xsabund(),
            "chatter": get_xschatter(),
            "cosmo": get_xscosmo(),
            "xsect": get_xsxsect(),
            "modelstrings": modelstrings}

def set_xsstate(state):
    if (type(state) == dict and
        state.has_key('abund') and
        state.has_key('chatter') and
        state.has_key('cosmo') and
        state.has_key('xsect') and
        state.has_key('modelstrings')):
        set_xsabund(state["abund"])
        set_xschatter(state["chatter"])
        set_xscosmo(state["cosmo"][0], state["cosmo"][1], state["cosmo"][2])
        set_xsxsect(state["xsect"])
        for name in state["modelstrings"].keys():
            set_xsxset(name, state["modelstrings"][name])

# The model classes are added to __all__ at the end of the file
__all__ = ('get_xschatter', 'get_xsabund', 'get_xscosmo', 'get_xsxsect',
           'set_xschatter', 'set_xsabund', 'set_xscosmo', 'set_xsxsect',
           'get_xsversion', 'set_xsxset', 'get_xsxset', 'set_xsstate',
           'get_xsstate')

class XSModel(ArithmeticModel):

    @modelCacher1d
    def calc(self, *args, **kwargs):
        return self._calc(*args, **kwargs)


class XSTableModel(XSModel):

    def __init__(self, filename, name='xstbl', parnames=(),
                 initvals=(), delta=(), mins=(), maxes=(), hardmins=(),
                 hardmaxes=(), nint=0, addmodel=False, addredshift=False):

        # make translation table to turn reserved characters into '_'
        bad = string.punctuation+string.whitespace
        tbl = string.maketrans(bad, '_'*len(bad))

        pars = []
        for ii in xrange(len(parnames)):
            isfrozen = True
            if nint > 0:
                isfrozen = False

            parname = parnames[ii].strip().lower().translate(tbl)
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
        func = _xspec.xsmtbl
        if self.addmodel:
            func = _xspec.xsatbl

        return func(p, filename=self.filename, *args, **kwargs)



class XSAdditiveModel(XSModel):

    def guess(self, dep, *args, **kwargs):
        if hasattr(self, 'norm'):
            norm = guess_amplitude(dep, *args)
            param_apply_limits(norm, self.norm, **kwargs)


class XSMultiplicativeModel(XSModel):
    pass


class XSapec(XSAdditiveModel):

    _calc =  _xspec.xsaped

    def __init__(self, name='apec'):
        self.kT = Parameter(name, 'kT', 1., 0.008, 64.0, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.redshift, self.norm))


class XSbapec(XSAdditiveModel):

    _calc =  _xspec.xsbape

    def __init__(self, name='bapec'):
        self.kT = Parameter(name, 'kT', 1., 0.008, 64.0, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.Velocity = Parameter(name, 'Velocity', 0., 0., 1.e6, 0.0, hugeval, 'km/s', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.Redshift, self.Velocity, self.norm))


class XSbbody(XSAdditiveModel):

    _calc =  _xspec.xsblbd

    def __init__(self, name='bbody'):
        self.kT = Parameter(name, 'kT', 3.0, 1.e-2, 100., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval, 'L39 / (D10)**2')
        XSAdditiveModel.__init__(self, name, (self.kT, self.norm))


class XSbbodyrad(XSAdditiveModel):

    _calc =  _xspec.xsbbrd

    def __init__(self, name='bbodyrad'):
        self.kT = Parameter(name, 'kT', 3., 1e-3, 100, 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.norm))


class XSbexrav(XSAdditiveModel):

    _calc =  _xspec.C_xsbexrav

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

    _calc =  _xspec.C_xsbexriv

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

    _calc =  _xspec.C_brokenPowerLaw

    def __init__(self, name='bknpower'):
        self.PhoIndx1 = Parameter(name, 'PhoIndx1', 1., -2., 9., -hugeval, hugeval)
        self.BreakE = Parameter(name, 'BreakE', 5., 1.e-2, 1.e6, 0.0, hugeval, 'keV')
        self.PhoIndx2 = Parameter(name, 'PhoIndx2', 2., -2., 9., -hugeval, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndx1, self.BreakE, self.PhoIndx2, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.BreakE, **kwargs)


class XSbkn2pow(XSAdditiveModel):

    _calc =  _xspec.C_broken2PowerLaw

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

    _calc =  _xspec.xsbmc

    def __init__(self, name='bmc'):
        self.kT = Parameter(name, 'kT', 1., 1.e-2, 100., 0.0, hugeval, 'keV')
        self.alpha = Parameter(name, 'alpha', 1., 1.e-2, 4.0, 0.0, hugeval)
        self.logA = Parameter(name, 'logA', 0.0, -6.0, 6.0, -hugeval, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.alpha, self.logA, self.norm))


class XSbremss(XSAdditiveModel):

    _calc =  _xspec.xsbrms

    def __init__(self, name='bremss'):
        self.kT = Parameter(name, 'kT', 7.0, 1.e-4, 100., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.norm))


class XSbvapec(XSAdditiveModel):

    _calc =  _xspec.xsbvpe

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


class XSc6mekl(XSAdditiveModel):

    _calc =  _xspec.c6mekl

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
        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5, self.CPcoef6, self.nH, self.abundanc, self.redshift, self.switch, self.norm))


class XSc6pmekl(XSAdditiveModel):

    _calc =  _xspec.c6pmekl

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
        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5, self.CPcoef6, self.nH, self.abundanc, self.redshift, self.switch, self.norm))


class XSc6pvmkl(XSAdditiveModel):

    _calc =  _xspec.c6pvmkl

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
        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5, self.CPcoef6, self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.switch, self.norm))


class XSc6vmekl(XSAdditiveModel):

    _calc =  _xspec.c6vmekl

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
        XSAdditiveModel.__init__(self, name, (self.CPcoef1, self.CPcoef2, self.CPcoef3, self.CPcoef4, self.CPcoef5, self.CPcoef6, self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.switch, self.norm))


class XScemekl(XSAdditiveModel):

    _calc =  _xspec.cemekl

    def __init__(self, name='cemekl'):
        self.alpha = Parameter(name, 'alpha', 1.0, 0.01, 10, 0.0, hugeval, frozen=True)
        self.Tmax = Parameter(name, 'Tmax', 1.0, 0.01, 1.e2, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1.0, 1.e-5, 1.e19, 0.0, hugeval, 'cm^-3', True)
        self.abundanc = Parameter(name, 'abundanc', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 1, 0, 1, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.alpha, self.Tmax, self.nH, self.abundanc, self.redshift, self.switch, self.norm))


class XScevmkl(XSAdditiveModel):

    _calc =  _xspec.C_cemVMekal

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
        XSAdditiveModel.__init__(self, name, (self.alpha, self.Tmax, self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.switch, self.norm))


class XScflow(XSAdditiveModel):

    _calc =  _xspec.C_xscflw

    def __init__(self, name='cflow'):
        self.slope = Parameter(name, 'slope', 0., -5., 5., -hugeval, hugeval)
        self.lowT = Parameter(name, 'lowT', 0.1, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.highT = Parameter(name, 'highT', 4., 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', .1, 1.e-10, 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.slope, self.lowT, self.highT, self.Abundanc, self.redshift, self.norm))


class XScompbb(XSAdditiveModel):

    _calc =  _xspec.compbb

    def __init__(self, name='compbb'):
        self.kT = Parameter(name, 'kT', 1.0, 1.e-2, 100., 0.0, hugeval, 'keV')
        self.kTe = Parameter(name, 'kTe', 50, 1., 200., 0.0, hugeval, 'keV', True)
        self.tau = Parameter(name, 'tau', 0.1, 0.0, 10., 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kTe, self.tau, self.norm))


class XScompLS(XSAdditiveModel):

    _calc =  _xspec.compls

    def __init__(self, name='compls'):
        self.kT = Parameter(name, 'kT', 2., .01, 10., 0.0, hugeval, 'keV')
        self.tau = Parameter(name, 'tau', 10, .001, 100., 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.tau, self.norm))


class XScompPS(XSAdditiveModel):

    _calc =  _xspec.C_xscompps

    def __init__(self, name='compps'):
        self.kTe = Parameter(name, 'kTe', 100., 20., 1.e5, 0.0, hugeval, 'keV')
        self.EleIndex = Parameter(name, 'EleIndex', 2., 0.0, 5., 0.0, hugeval, frozen=True)
        self.Gmin = Parameter(name, 'Gmin', -1., -1., 10., -hugeval, hugeval, frozen=True)
        self.Gmax = Parameter(name, 'Gmax', 1.e3, 10., 1.e4, 0.0, hugeval, frozen=True)
        self.kTbb = Parameter(name, 'kTbb', 0.1, 0.001, 10., 0.0, hugeval, 'keV', True)
        self.tauy = Parameter(name, 'tauy', 1.0, 0.05, 3.0, 0.0, hugeval)
        self.geom = Parameter(name, 'geom', 0.0, -5.0, 4.0, -hugeval, hugeval, frozen=True)
        self.HRcyl = Parameter(name, 'HRcyl', 1.0, 0.5, 2.0, 0.0, hugeval, frozen=True)
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
        XSAdditiveModel.__init__(self, name, (self.kTe, self.EleIndex, self.Gmin, self.Gmax, self.kTbb, self.tauy, self.geom, self.HRcyl, self.cosIncl, self.cov_frac, self.rel_refl, self.Fe_ab_re, self.Me_ab, self.xi, self.Tdisk, self.Betor10, self.Rin, self.Rout, self.redshift, self.norm))


class XScompST(XSAdditiveModel):

    _calc =  _xspec.compst

    def __init__(self, name='compst'):
        self.kT = Parameter(name, 'kT', 2., .01, 100., 0.0, hugeval, 'keV')
        self.tau = Parameter(name, 'tau', 10, .001, 100., 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.tau, self.norm))


class XScompTT(XSAdditiveModel):

    _calc =  _xspec.xstitg

    def __init__(self, name='comptt'):
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.T0 = Parameter(name, 'T0', 0.1, .01, 100., 0.0, hugeval, 'keV')
        self.kT = Parameter(name, 'kT', 50., 2.0, 500., 0.0, hugeval, 'keV')
        self.taup = Parameter(name, 'taup', 1., .01, 100., 0.0, hugeval)
        self.approx = Parameter(name, 'approx', 1.0, 0.0, 5.0, 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.redshift, self.T0, self.kT, self.taup, self.approx, self.norm))


class XScutoffpl(XSAdditiveModel):

    _calc =  _xspec.C_cutoffPowerLaw

    def __init__(self, name='cutoffpl'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 1., -2., 9., -hugeval, hugeval)
        self.HighECut = Parameter(name, 'HighECut', 15., 1., 500., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.HighECut, self.norm))


class XSdisk(XSAdditiveModel):

    _calc =  _xspec.disk

    def __init__(self, name='disk'):
        self.accrate = Parameter(name, 'accrate', 1., 1e-3, 9., 0.0, hugeval)
        self.NSmass = Parameter(name, 'NSmass', 1.4, .4, 10., 0.0, hugeval, units='Msun', frozen=True)
        self.Rinn = Parameter(name, 'Rinn', 1.03, 1.01, 1.03, 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.accrate, self.NSmass, self.Rinn, self.norm))


class XSdiskir(XSAdditiveModel):

    _calc =  _xspec.diskir

    def __init__(self, name='diskir'):
        self.kT_disk = Parameter(name, 'kT_disk', 1.0, 0.01, 5., 0.0, hugeval, 'keV')
        self.Gamma = Parameter(name, 'Gamma', 1.7, 1.001, 5., 0.0, hugeval)
        self.kT_e = Parameter(name, 'kT_e', 100., 5., 1.e3, 0.0, hugeval, 'keV')
        self.LcLd = Parameter(name, 'LcLd', 0.1, 0., 10., 0.0, hugeval)
        self.fin = Parameter(name, 'fin', 1.e-1, 0.0, 1., 0.0, hugeval, frozen=True)
        self.rirr = Parameter(name, 'rirr', 1.2, 1.0001, 10., 1.0001, hugeval)
        self.fout = Parameter(name, 'fout', 1.e-4, 0.0, 1.e-1, 0.0, hugeval)
        self.logrout = Parameter(name, 'logrout', 5.0, 3.0, 7.0, 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_disk, self.Gamma, self.kT_e, self.LcLd, self.fin, self.rirr, self.fout, self.logrout, self.norm))


class XSdiskbb(XSAdditiveModel):

    _calc =  _xspec.xsdskb

    def __init__(self, name='diskbb'):
        self.Tin = Parameter(name, 'Tin', 1., 0., 1000., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Tin, self.norm))


class XSdiskline(XSAdditiveModel):

    _calc =  _xspec.xsdili

    def __init__(self, name='diskline'):
        self.LineE = Parameter(name, 'LineE', 6.7, 0., 100., 0.0, hugeval, 'keV')
        self.Betor10 = Parameter(name, 'Betor10', -2., -10., 20., -hugeval, hugeval, frozen=True)
        self.RinM = Parameter(name, 'RinM', 10., 6., 1000., 0.0, hugeval, frozen=True)
        self.RoutM = Parameter(name, 'RoutM', 1000., 0., 1000000., 0.0, hugeval, frozen=True)
        self.Incl = Parameter(name, 'Incl', 30., 0., 90., 0.0, hugeval, 'deg')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Betor10, self.RinM, self.RoutM, self.Incl, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


class XSdiskm(XSAdditiveModel):

    _calc =  _xspec.diskm

    def __init__(self, name='diskm'):
        self.accrate = Parameter(name, 'accrate', 1., 1e-3, 9., 0.0, hugeval)
        self.NSmass = Parameter(name, 'NSmass', 1.4, .4, 10., 0.0, hugeval, units='Msun', frozen=True)
        self.Rinn = Parameter(name, 'Rinn', 1.03, 1.01, 1.03, 0.0, hugeval, frozen=True)
        self.alpha = Parameter(name, 'alpha', 1., .01, 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.accrate, self.NSmass, self.Rinn, self.alpha, self.norm))


class XSdisko(XSAdditiveModel):

    _calc =  _xspec.disko

    def __init__(self, name='disko'):
        self.accrate = Parameter(name, 'accrate', 1., 1e-3, 9., 0.0, hugeval)
        self.NSmass = Parameter(name, 'NSmass', 1.4, .4, 10., 0.0, hugeval, units='Msun', frozen=True)
        self.Rinn = Parameter(name, 'Rinn', 1.03, 1.01, 1.03, 0.0, hugeval, frozen=True)
        self.alpha = Parameter(name, 'alpha', 1., .01, 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.accrate, self.NSmass, self.Rinn, self.alpha, self.norm))


class XSdiskpbb(XSAdditiveModel):

    _calc =  _xspec.diskpbb

    def __init__(self, name='diskpbb'):
        self.Tin = Parameter(name, 'Tin', 1.0, 0.1, 10.0, 0.0, hugeval, 'keV')
        self.p = Parameter(name, 'p', 0.75, 0.5, 1.0, 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Tin, self.p, self.norm))


class XSdiskpn(XSAdditiveModel):

    _calc =  _xspec.xsdiskpn

    def __init__(self, name='diskpn'):
        self.T_max = Parameter(name, 'T_max', 1., 1e-3, 100, 0.0, hugeval, 'keV')
        self.R_in = Parameter(name, 'R_in', 6., 6., 1000., 0.0, hugeval, 'R_g')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.T_max, self.R_in, self.norm))


class XSequil(XSAdditiveModel):

    _calc =  _xspec.C_equil

    def __init__(self, name='equil'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.redshift, self.norm))


class XSexpdec(XSAdditiveModel):

    _calc =  _xspec.xsxpdec

    def __init__(self, name='expdec'):
        self.factor = Parameter(name, 'factor', 1.0, 0., 100.0, 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.factor, self.norm))


class XSezdiskbb(XSAdditiveModel):

    _calc =  _xspec.ezdiskbb

    def __init__(self, name='ezdiskbb'):
        self.T_max = Parameter(name, 'T_max', 1., 0.01, 100., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.T_max, self.norm))


class XSgaussian(XSAdditiveModel):

    _calc =  _xspec.xsgaul

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

    _calc =  _xspec.C_gnei

    def __init__(self, name='gnei'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.kT_ave = Parameter(name, 'kT_ave', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.Tau, self.kT_ave, self.redshift, self.norm))


class XSgrad(XSAdditiveModel):

    _calc =  _xspec.grad

    def __init__(self, name='grad'):
        self.D = Parameter(name, 'D', 10.0, 0.0, 10000., 0.0, hugeval, 'kpc', True)
        self.i = Parameter(name, 'i', 0.0, 0.0, 90.0, 0.0, hugeval, 'deg', True)
        self.Mass = Parameter(name, 'Mass', 1.0, 0.0, 100.0, 0.0, hugeval, 'solar')
        self.Mdot = Parameter(name, 'Mdot', 1.0, 0.0, 100.0, 0.0, hugeval, '1e18')
        self.TclTef = Parameter(name, 'TclTef', 1.7, 1.0, 10.0, 0.0, hugeval, frozen=True)
        self.refflag = Parameter(name, 'refflag', 1.0, -1.0, 1.0, -hugeval, hugeval, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.D, self.i, self.Mass, self.Mdot, self.TclTef, self.refflag, self.norm))


class XSgrbm(XSAdditiveModel):

    _calc =  _xspec.xsgrbm

    def __init__(self, name='grbm'):
        self.alpha = Parameter(name, 'alpha', -1., -3., +2., -hugeval, hugeval)
        self.beta = Parameter(name, 'beta', -2., -5., +2., -hugeval, hugeval)
        self.temp = Parameter(name, 'temp', +300., +50., +1000., 0.0, hugeval, 'keV')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.alpha, self.beta, self.temp, self.norm))


class XSkerrbb(XSAdditiveModel):

    _calc =  _xspec.C_kerrbb

    def __init__(self, name='kerrbb'):
        self.eta = Parameter(name, 'eta', 0., 0., 1.0, 0.0, hugeval, frozen=True)
        self.a = Parameter(name, 'a', 0., -1., 0.9999, -hugeval, hugeval)
        self.i = Parameter(name, 'i', 30., 0., 85., 0.0, hugeval, 'deg', True)
        self.Mbh = Parameter(name, 'Mbh', 1., 0., 100., 0.0, hugeval, 'Msun')
        self.Mdd = Parameter(name, 'Mdd', 1., 0., 1000., 0.0, hugeval, 'Mdd0')
        self.Dbh = Parameter(name, 'Dbh', 10., 0., 10000., 0.0, hugeval, 'kpc', True)
        self.hd = Parameter(name, 'hd', 1.7, 1., 10., 0.0, hugeval, frozen=True)
        self.rflag = Parameter(name, 'rflag', 1., -100., 100., -hugeval, hugeval, alwaysfrozen=True)
        self.lflag = Parameter(name, 'lflag', 0., -100., 100., -hugeval, hugeval, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.eta, self.a, self.i, self.Mbh, self.Mdd, self.Dbh, self.hd, self.rflag, self.lflag, self.norm))


class XSkerrd(XSAdditiveModel):

    _calc =  _xspec.C_kerrdisk

    def __init__(self, name='kerrd'):
        self.distance = Parameter(name, 'distance', 1., 0.01, 1000., 0.0, hugeval, 'kpc', True)
        self.TcolTeff = Parameter(name, 'TcolTeff', 1.5, 1.0, 2.0, 0.0, hugeval, frozen=True)
        self.M = Parameter(name, 'M', 1.0, 0.1, 100., 0.0, hugeval, 'solar')
        self.Mdot = Parameter(name, 'Mdot', 1.0, 0.01, 100., 0.0, hugeval, '1e18')
        self.Incl = Parameter(name, 'Incl', 30., 0., 90., 0.0, hugeval, 'deg', True)
        self.Rin = Parameter(name, 'Rin', 1.235, 1.235, 100., 0.0, hugeval, 'Rg', True)
        self.Rout = Parameter(name, 'Rout', 1e5, 1e4, 1e8, 0.0, hugeval, 'Rg', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.distance, self.TcolTeff, self.M, self.Mdot, self.Incl, self.Rin, self.Rout, self.norm))


class XSkerrdisk(XSAdditiveModel):

    _calc =  _xspec.spin

    def __init__(self, name='kerrdisk'):
        self.lineE = Parameter(name, 'lineE', 6.4, 0.1, 100., 0.0, hugeval, 'keV')
        self.Index1 = Parameter(name, 'Index1', 3., -10., 10., -hugeval, hugeval, frozen=True)
        self.Index2 = Parameter(name, 'Index2', 3., -10., 10., -hugeval, hugeval, frozen=True)
        self.r_brg = Parameter(name, 'r_brg', 6.0, 1.0, 400., 0.0, hugeval, frozen=True)
        self.a = Parameter(name, 'a', 0.998, 0.0, 0.998, 0.0, hugeval)
        self.Incl = Parameter(name, 'Incl', 30., 0., 90., 0.0, hugeval, 'deg', True)
        self.Rinms = Parameter(name, 'Rinms', 1.0, 1.0, 400., 0.0, hugeval, frozen=True)
        self.Routms = Parameter(name, 'Routms', 400., 1.0, 400., 0.0, hugeval, frozen=True)
        self.z = Parameter(name, 'z', 0., 0., 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.lineE, self.Index1, self.Index2, self.r_brg, self.a, self.Incl, self.Rinms, self.Routms, self.z, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.lineE, **kwargs)


class XSlaor(XSAdditiveModel):

    _calc =  _xspec.C_xslaor

    def __init__(self, name='laor'):
        self.lineE = Parameter(name, 'lineE', 6.4, 0., 100., 0.0, hugeval, 'keV')
        self.Index = Parameter(name, 'Index', 3., -10., 10., -hugeval, hugeval, frozen=True)
        self.RinG = Parameter(name, 'RinG', 1.235, 1.235, 400., 0.0, hugeval, frozen=True)
        self.RoutG = Parameter(name, 'RoutG', 400., 1.235, 400., 0.0, hugeval, frozen=True)
        self.Incl = Parameter(name, 'Incl', 30., 0., 90., 0.0, hugeval, 'deg', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.lineE, self.Index, self.RinG, self.RoutG, self.Incl, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.lineE, **kwargs)


class XSlaor2(XSAdditiveModel):

    _calc =  _xspec.C_laor2

    def __init__(self, name='laor2'):
        self.lineE = Parameter(name, 'lineE', 6.4, 0., 100., 0.0, hugeval, 'keV')
        self.Index = Parameter(name, 'Index', 3., -10., 10., -hugeval, hugeval, frozen=True)
        self.RinG = Parameter(name, 'RinG', 1.235, 1.235, 400., 0.0, hugeval, frozen=True)
        self.RoutG = Parameter(name, 'RoutG', 400., 1.235, 400., 0.0, hugeval, frozen=True)
        self.Incl = Parameter(name, 'Incl', 30., 0., 90., 0.0, hugeval, 'deg', True)
        self.Rbreak = Parameter(name, 'Rbreak', 20., 1.235, 400., 0.0, hugeval, frozen=True)
        self.Index1 = Parameter(name, 'Index1', 3., -10., 10., -hugeval, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.lineE, self.Index, self.RinG, self.RoutG, self.Incl, self.Rbreak, self.Index1, self.norm))

    def guess(self, dep, *args, **kwargs):
        XSAdditiveModel.guess(self, dep, *args, **kwargs)
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.lineE, **kwargs)


class XSlorentz(XSAdditiveModel):

    _calc =  _xspec.xslorz

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

    _calc =  _xspec.xsmeka

    def __init__(self, name='meka'):
        self.kT = Parameter(name, 'kT', 1., 1.e-3, 1.e2, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1., 1.e-5, 1.e19, 0.0, hugeval, 'cm-3', True)
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.nH, self.Abundanc, self.redshift, self.norm))


class XSmekal(XSAdditiveModel):

    _calc =  _xspec.xsmekl

    def __init__(self, name='mekal'):
        self.kT = Parameter(name, 'kT', 1., 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1., 1.e-5, 1.e19, 0.0, hugeval, 'cm-3', True)
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.nH, self.Abundanc, self.redshift, self.switch, self.norm))


class XSmkcflow(XSAdditiveModel):

    _calc =  _xspec.C_xsmkcf

    def __init__(self, name='mkcflow'):
        self.lowT = Parameter(name, 'lowT', 0.1, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.highT = Parameter(name, 'highT', 4., 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.lowT, self.highT, self.Abundanc, self.redshift, self.switch, self.norm))


class XSnei(XSAdditiveModel):

    _calc =  _xspec.C_nei

    def __init__(self, name='nei'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.Tau, self.redshift, self.norm))


class XSrnei(XSAdditiveModel):

    _calc =  _xspec.C_rnei

    def __init__(self, name='rnei'):
        self.kT = Parameter(name, 'kT', 0.5, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.kT_init = Parameter(name, 'kT_init', 1.0, 0.0808, 79.9, 0.0808, 79.9, units='keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, 10000, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.0e8, 5.0e13, 1.0e8, 5.0e13, units='s/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.kT_init, self.Abundanc, self.Tau, self.redshift, self.norm))


class XSvrnei(XSAdditiveModel):

    _calc =  _xspec.C_vrnei

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

    _calc =  _xspec.C_vvrnei

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

    _calc =  _xspec.C_npshock

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

    _calc =  _xspec.nsa

    def __init__(self, name='nsa'):
        self.LogT_eff = Parameter(name, 'LogT_eff', 6.0, 5.0, 7.0, 0.0, hugeval, 'K')
        self.M_ns = Parameter(name, 'M_ns', 1.4, 0.5, 2.5, 0.0, hugeval, 'Msun')
        self.R_ns = Parameter(name, 'R_ns', 10.0, 5.0, 20., 0.0, hugeval, 'km')
        self.MagField = Parameter(name, 'MagField', 0.0, 0.0, 5.e13, 0.0, hugeval, 'G', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LogT_eff, self.M_ns, self.R_ns, self.MagField, self.norm))


class XSnsagrav(XSAdditiveModel):

    _calc =  _xspec.nsagrav

    def __init__(self, name='nsagrav'):
        self.LogT_eff = Parameter(name, 'LogT_eff', 6.0, 5.5, 6.5, 0.0, hugeval, 'K')
        self.NSmass = Parameter(name, 'NSmass', 1.4, 0.3, 2.5, 0.0, hugeval, 'Msun')
        self.NSrad = Parameter(name, 'NSrad', 10.0, 6.0, 20., 0.0, hugeval, 'km')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LogT_eff, self.NSmass, self.NSrad, self.norm))


class XSnsatmos(XSAdditiveModel):

    _calc =  _xspec.nsatmos

    def __init__(self, name='nsatmos'):
        self.LogT_eff = Parameter(name, 'LogT_eff', 6.0, 5.0, 6.5, 0.0, hugeval, 'K')
        self.M_ns = Parameter(name, 'M_ns', 1.4, 0.5, 3.0, 0.0, hugeval, 'Msun')
        self.R_ns = Parameter(name, 'R_ns', 10.0, 5.0, 30., 0.0, hugeval, 'km')
        self.dist = Parameter(name, 'dist', 10.0, 0.1, 100.0, 0.0, hugeval, 'kpc')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LogT_eff, self.M_ns, self.R_ns, self.dist, self.norm))


class XSnsmax(XSAdditiveModel):

    _calc =  _xspec.nsmax

    def __init__(self, name='nsmax'):
        self.logTeff = Parameter(name, 'logTeff', 6.0, 5.5, 6.8, 0.0, hugeval, 'K')
        self.redshift = Parameter(name, 'redshift', 0.1, 0.0, 1.5, 0.0, 2.0)
        self.specfile = Parameter(name, 'specfile', 1200, 0, 1.0e6, 0, 1.0e6, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.logTeff, self.redshift, self.specfile, self.norm))


class XSnsmaxg(XSAdditiveModel):

    _calc =  _xspec.nsmaxg

    def __init__(self, name='nsmaxg'):
        self.logTeff = Parameter(name, 'logTeff', 6.0, 5.5, 6.9, 5.5, 6.9, units='K')
        self.M_ns = Parameter(name, 'M_ns', 1.4, 0.5, 3.0, 0.5, 3.0, units='Msun')
        self.R_ns = Parameter(name, 'R_ns', 10.0, 5.0, 30.0, 5.0, 30.0, units='km')
        self.dist = Parameter(name, 'dist', 1.0, 0.01, 100.0, 0.01, 100.0, units='kpc')
        self.specfile = Parameter(name, 'specfile', 1200, 0, 1.0e6, 0, 1.0e6, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.logTeff, self.M_ns, self.R_ns, self.dist, self.specfile, self.norm))


class XSnsx(XSAdditiveModel):

    _calc =  _xspec.nsx

    def __init__(self, name='nsx'):
        self.logTeff = Parameter(name, 'logTeff', 6.0, 5.5, 6.7, 5.5, 6.7, units='K')
        self.M_ns = Parameter(name, 'M_ns', 1.4, 0.5, 3.0, 0.5, 3.0, units='Msun')
        self.R_ns = Parameter(name, 'R_ns', 10.0, 5.0, 30.0, 5.0, 30.0, units='km')
        self.dist = Parameter(name, 'dist', 1.0, 0.01, 100.0, 0.01, 100.0, units='kpc')
        self.specfile = Parameter(name, 'specfile', 6, 0, 1.0e6, 0, 1.0e6, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.logTeff, self.M_ns, self.R_ns, self.dist, self.specfile, self.norm))

class XSnteea(XSAdditiveModel):

    _calc =  _xspec.C_xsnteea

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

    _calc =  _xspec.C_nthcomp

    def __init__(self, name='nthcomp'):
        self.Gamma = Parameter(name, 'Gamma', 1.7, 1.001, 5., 1.001, 10.)
        self.kT_e = Parameter(name, 'kT_e', 100., 5., 1.e3, 1., 1.e3, 'keV')
        self.kT_bb = Parameter(name, 'kT_bb', 0.1, 1.e-3, 10., 1.e-3, 10., 'keV', True)
        self.inp_type = Parameter(name, 'inp_type', 0., 0., 1., 0., 1., '0/1', True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, 10., frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Gamma, self.kT_e, self.kT_bb, self.inp_type, self.redshift, self.norm))


class XSpegpwrlw(XSAdditiveModel):

    _calc =  _xspec.xspegp

    def __init__(self, name='pegpwrlw'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 1., -2., 9., -hugeval, hugeval)
        self.eMin = Parameter(name, 'eMin', 2., -100., 1.e10, -hugeval, hugeval, 'keV', True)
        self.eMax = Parameter(name, 'eMax', 10., -100., 1.e10, -hugeval, hugeval, 'keV', True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.eMin, self.eMax, self.norm))


class XSpexrav(XSAdditiveModel):

    _calc =  _xspec.C_xspexrav

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

    _calc =  _xspec.C_xspexriv

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

    _calc =  _xspec.xsp1tr

    def __init__(self, name='plcabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.nmax = Parameter(name, 'nmax', 1, alwaysfrozen=True)
        self.FeAbun = Parameter(name, 'FeAbun', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.FeKedge = Parameter(name, 'FeKedge', 7.11, 7., 10., 0.0, hugeval, 'KeV', True)
        self.PhoIndex = Parameter(name, 'PhoIndex', 2., -2., 9., -hugeval, hugeval)
        self.HighECut = Parameter(name, 'HighECut', 25., 1., 50., 0.0, 90., 'keV', True)
        self.foldE = Parameter(name, 'foldE', 100., 1., 1.e6, 0.0, hugeval, frozen=True)
        self.acrit = Parameter(name, 'acrit', 1., 0.0, 1.0, 0.0, hugeval, frozen=True)
        self.FAST = Parameter(name, 'FAST', 0, alwaysfrozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.nH, self.nmax, self.FeAbun, self.FeKedge, self.PhoIndex, self.HighECut, self.foldE, self.acrit, self.FAST, self.redshift, self.norm))


class XSpowerlaw(XSAdditiveModel):

    _calc =  _xspec.C_powerLaw

    def __init__(self, name='powerlaw'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 1., -2., 9., -hugeval, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.norm))


class XSposm(XSAdditiveModel):

    _calc =  _xspec.xsposm

    def __init__(self, name='posm'):
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.norm,))


class XSpshock(XSAdditiveModel):

    _calc =  _xspec.C_pshock

    def __init__(self, name='pshock'):
        self.kT = Parameter(name, 'kT', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau_l = Parameter(name, 'Tau_l', 0.0, 0., 5.e13, 0.0, hugeval, 's/cm^3', True)
        self.Tau_u = Parameter(name, 'Tau_u', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.Tau_l, self.Tau_u, self.redshift, self.norm))


class XSraymond(XSAdditiveModel):

    _calc =  _xspec.xsrays

    def __init__(self, name='raymond'):
        self.kT = Parameter(name, 'kT', 1., 0.008, 64.0, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.Abundanc, self.redshift, self.norm))


class XSredge(XSAdditiveModel):

    _calc =  _xspec.xredge

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

    _calc =  _xspec.xsrefsch

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

    _calc =  _xspec.C_sedov

    def __init__(self, name='sedov'):
        self.kT_a = Parameter(name, 'kT_a', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.kT_b = Parameter(name, 'kT_b', 0.5, 0.0100, 79.9, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1.0, 0., 1000., 0.0, hugeval, frozen=True)
        self.Tau = Parameter(name, 'Tau', 1.e11, 1.e8, 5.e13, 0.0, hugeval, 's/cm^3')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT_a, self.kT_b, self.Abundanc, self.Tau, self.redshift, self.norm))

class XSsirf(XSAdditiveModel):

    _calc =  _xspec.C_sirf

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

    _calc =  _xspec.srcut

    def __init__(self, name='srcut'):
        self.alpha = Parameter(name, 'alpha', 0.5, 0.3, 0.8, 0.0, hugeval)
        self.breakfreq = Parameter(name, 'breakfreq', 2.42E17, 1.E15, 1.E19, 0.0, hugeval, 'Hz')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.alpha, self.breakfreq, self.norm))


class XSsresc(XSAdditiveModel):

    _calc =  _xspec.sresc

    def __init__(self, name='sresc'):
        self.alpha = Parameter(name, 'alpha', 0.5, 0.3, 0.8, 0.0, hugeval)
        self.rolloff = Parameter(name, 'rolloff', 2.42E17, 1.E15, 1.E19, 0.0, hugeval, 'Hz')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.alpha, self.rolloff, self.norm))


class XSstep(XSAdditiveModel):

    _calc =  _xspec.xsstep

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

    _calc =  _xspec.xsvape

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

    _calc =  _xspec.xsbrmv

    def __init__(self, name='vbremss'):
        self.kT = Parameter(name, 'kT', 3.0, 1.e-2, 100., 0.0, hugeval, 'keV')
        self.HeH = Parameter(name, 'HeH', 1.0, 0., 100., 0.0, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.HeH, self.norm))


class XSvequil(XSAdditiveModel):

    _calc =  _xspec.C_vequil

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

    _calc =  _xspec.C_vgnei

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
        self.kT_ave = Parameter(name, 'kT_ave', 1.0, 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Tau, self.kT_ave, self.redshift, self.norm))


class XSvvgnei(XSAdditiveModel):

    _calc =  _xspec.C_vvgnei

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
        self.kT_ave = Parameter(name, 'kT_ave', 1.0, 0.0808, 79.9, 0.0808, 79.9, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Tau, self.kT_ave, self.redshift, self.norm))

class XSvmeka(XSAdditiveModel):

    _calc =  _xspec.xsvmek

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

    _calc =  _xspec.xsvmkl

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
        XSAdditiveModel.__init__(self, name, (self.kT, self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.switch, self.norm))


class XSvmcflow(XSAdditiveModel):

    _calc =  _xspec.C_xsvmcf

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
        XSAdditiveModel.__init__(self, name, (self.lowT, self.highT, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.switch, self.norm))


class XSvnei(XSAdditiveModel):

    _calc =  _xspec.C_vnei

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

    _calc =  _xspec.C_vvnei

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

    _calc =  _xspec.C_vnpshock

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

    _calc =  _xspec.C_vvnpshock

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

    _calc =  _xspec.C_vpshock

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

    _calc =  _xspec.C_vvpshock

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

    _calc =  _xspec.xsvrys

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

    _calc =  _xspec.C_vsedov

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

    _calc =  _xspec.C_vvsedov

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

    _calc =  _xspec.xszbod

    def __init__(self, name='zbbody'):
        self.kT = Parameter(name, 'kT', 3.0, 1.e-2, 100., 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.redshift, self.norm))


class XSzbremss(XSAdditiveModel):

    _calc =  _xspec.xszbrm

    def __init__(self, name='zbremss'):
        self.kT = Parameter(name, 'kT', 7.0, 1.e-4, 100., 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kT, self.redshift, self.norm))


class XSzgauss(XSAdditiveModel):

    _calc =  _xspec.C_xszgau

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


class XSzpowerlw(XSAdditiveModel):

    _calc =  _xspec.C_zpowerLaw

    def __init__(self, name='zpowerlw'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 1., -2., 9., -hugeval, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.redshift, self.norm))


class XSabsori(XSMultiplicativeModel):

    _calc =  _xspec.C_xsabsori

    def __init__(self, name='absori'):
        self.PhoIndex = Parameter(name, 'PhoIndex', 2., 0., 4., 0.0, hugeval, frozen=True)
        self.nH = Parameter(name, 'nH', 1., 0., 100., 0.0, hugeval, '10^22 atoms / cm^2')
        self.Temp_abs = Parameter(name, 'Temp_abs', 3.e4, 1.e4, 1.e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 1., 0., 1.e3, 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0., 1.e6, 0.0, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.PhoIndex, self.nH, self.Temp_abs, self.xi, self.redshift, self.Fe_abund))


class XSacisabs(XSMultiplicativeModel):

    _calc =  _xspec.acisabs

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

    _calc =  _xspec.xscnst

    def __init__(self, name='constant'):
        self.factor = Parameter(name, 'factor', 1., 0.0, 1.e10, 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.factor,))


class XScabs(XSMultiplicativeModel):

    _calc =  _xspec.xscabs

    def __init__(self, name='cabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        XSMultiplicativeModel.__init__(self, name, (self.nH,))


class XScyclabs(XSMultiplicativeModel):

    _calc =  _xspec.xscycl

    def __init__(self, name='cyclabs'):
        self.Depth0 = Parameter(name, 'Depth0', 2.0, 0., 100., 0.0, hugeval)
        self.E0 = Parameter(name, 'E0', 30.0, 1.0, 100., 0.0, hugeval, 'keV')
        self.Width0 = Parameter(name, 'Width0', 10.0, 1.0, 100., 0.0, hugeval, 'keV', True)
        self.Depth2 = Parameter(name, 'Depth2', 0.0, 0., 100., 0.0, hugeval, frozen=True)
        self.Width2 = Parameter(name, 'Width2', 20.0, 1.0, 100., 0.0, hugeval, 'keV', True)
        XSMultiplicativeModel.__init__(self, name, (self.Depth0, self.E0, self.Width0, self.Depth2, self.Width2))


class XSdust(XSMultiplicativeModel):

    _calc =  _xspec.xsdust

    def __init__(self, name='dust'):
        self.Frac = Parameter(name, 'Frac', 0.066, 0., 1., 0.0, hugeval, frozen=True)
        self.Halosz = Parameter(name, 'Halosz', 2., 0., 1.e5, 0.0, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.Frac, self.Halosz))


class XSedge(XSMultiplicativeModel):

    _calc =  _xspec.xsedge

    def __init__(self, name='edge'):
        self.edgeE = Parameter(name, 'edgeE', 7.0, 0., 100., 0.0, hugeval, 'keV')
        self.MaxTau = Parameter(name, 'MaxTau', 1., 0., 5., 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.edgeE, self.MaxTau))


class XSexpabs(XSMultiplicativeModel):

    _calc =  _xspec.xsabsc

    def __init__(self, name='expabs'):
        self.LowECut = Parameter(name, 'LowECut', 2., 0., 100., 0.0, hugeval, 'keV')
        XSMultiplicativeModel.__init__(self, name, (self.LowECut,))


class XSexpfac(XSMultiplicativeModel):

    _calc =  _xspec.xsexp

    def __init__(self, name='expfac'):
        self.Ampl = Parameter(name, 'Ampl', 1., 0., 1.e5, 0.0, hugeval)
        self.Factor = Parameter(name, 'Factor', 1., 0., 1.e5, 0.0, hugeval)
        self.StartE = Parameter(name, 'StartE', 0.5, 0., 1.e5, 0.0, hugeval, 'keV', True)
        XSMultiplicativeModel.__init__(self, name, (self.Ampl, self.Factor, self.StartE))


class XSgabs(XSMultiplicativeModel):

    _calc =  _xspec.C_gaussianAbsorptionLine

    def __init__(self, name='gabs'):
        self.LineE = Parameter(name, 'LineE', 1.0, 0., 1.e6, 0.0, hugeval, 'keV')
        self.Sigma = Parameter(name, 'Sigma', 0.01, 0., 10., 0.0, hugeval, 'keV')
        self.Tau = Parameter(name, 'Tau', 1.0, 0., 1.e6, 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.LineE, self.Sigma, self.Tau))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


class XShighecut(XSMultiplicativeModel):

    _calc =  _xspec.xshecu

    def __init__(self, name='highecut'):
        self.cutoffE = Parameter(name, 'cutoffE', 10., 1.e-2, 1.e6, 0.0, hugeval, 'keV')
        self.foldE = Parameter(name, 'foldE', 15., 1.e-2, 1.e6, 0.0, hugeval, 'keV')
        XSMultiplicativeModel.__init__(self, name, (self.cutoffE, self.foldE))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.cutoffE, **kwargs)


class XShrefl(XSMultiplicativeModel):

    _calc =  _xspec.xshrfl

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

    _calc =  _xspec.xsntch

    def __init__(self, name='notch'):
        self.LineE = Parameter(name, 'LineE', 3.5, 0., 20., 0.0, hugeval, 'keV')
        self.Width = Parameter(name, 'Width', 1., 0., 20., 0.0, hugeval, 'keV')
        self.CvrFract = Parameter(name, 'CvrFract', 1., 0., 1., 0.0, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.LineE, self.Width, self.CvrFract))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.LineE, **kwargs)


class XSpcfabs(XSMultiplicativeModel):

    _calc =  _xspec.xsabsp

    def __init__(self, name='pcfabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.CvrFract = Parameter(name, 'CvrFract', 0.5, 0.05, 0.95, 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.CvrFract))


class XSphabs(XSMultiplicativeModel):

    _calc = _xspec.xsphab

    def __init__(self, name='phabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        XSMultiplicativeModel.__init__(self, name, (self.nH,))


class XSplabs(XSMultiplicativeModel):

    _calc =  _xspec.xsplab

    def __init__(self, name='plabs'):
        self.index = Parameter(name, 'index', 2.0, 0.0, 5., 0.0, hugeval)
        self.coef = Parameter(name, 'coef', 1.0, 0.0, 100., 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.index, self.coef))


class XSpwab(XSMultiplicativeModel):

    _calc =  _xspec.C_xspwab

    def __init__(self, name='pwab'):
        self.nHmin = Parameter(name, 'nHmin', 1., 1.e-7, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.nHmax = Parameter(name, 'nHmax', 2., 1.e-7, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.beta = Parameter(name, 'beta', 1.0, -10., 10, -hugeval, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nHmin, self.nHmax, self.beta))


class XSredden(XSMultiplicativeModel):

    _calc =  _xspec.xscred

    def __init__(self, name='redden'):
        self.EBV = Parameter(name, 'EBV', 0.05, 0., 10., 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.EBV,))


class XSsmedge(XSMultiplicativeModel):

    _calc =  _xspec.xssmdg

    def __init__(self, name='smedge'):
        self.edgeE = Parameter(name, 'edgeE', 7.0, 0.1, 100., 0.0, hugeval, 'keV')
        self.MaxTau = Parameter(name, 'MaxTau', 1., 0., 5., 0.0, hugeval)
        self.index = Parameter(name, 'index', -2.67, -10., 10., -hugeval, hugeval, frozen=True)
        self.width = Parameter(name, 'width', 10., 0.01, 100., 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.edgeE, self.MaxTau, self.index, self.width))


class XSspexpcut(XSMultiplicativeModel):

    _calc =  _xspec.C_superExpCutoff

    def __init__(self, name='spexpcut'):
        self.Ecut = Parameter(name, 'Ecut', 10.0, 0.0, 1e6, 0.0, hugeval, 'keV')
        self.alpha = Parameter(name, 'alpha', 1.0, -5.0, 5.0, -hugeval, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.Ecut, self.alpha))


class XSspline(XSMultiplicativeModel):

    _calc =  _xspec.xsspln

    def __init__(self, name='spline'):
        self.Estart = Parameter(name, 'Estart', 0.1, 0., 100., 0.0, hugeval, 'keV')
        self.Ystart = Parameter(name, 'Ystart', 1., -1.e6, 1.e6, -hugeval, hugeval)
        self.Yend = Parameter(name, 'Yend', 1., -1.e6, 1.e6, -hugeval, hugeval)
        self.YPstart = Parameter(name, 'YPstart', 0., -1.e6, 1.e6, -hugeval, hugeval)
        self.YPend = Parameter(name, 'YPend', 0., -1.e6, 1.e6, -hugeval, hugeval)
        self.Eend = Parameter(name, 'Eend', 15., 0., 100., 0.0, hugeval, 'keV')
        XSMultiplicativeModel.__init__(self, name, (self.Estart, self.Ystart, self.Yend, self.YPstart, self.YPend, self.Eend))


class XSSSS_ice(XSMultiplicativeModel):

    _calc =  _xspec.xssssi

    def __init__(self, name='sss_ice'):
        self.clumps = Parameter(name, 'clumps', 0.0, 0., 10., 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.clumps,))


class XSswind1(XSMultiplicativeModel):

    _calc =  _xspec.swind1

    def __init__(self, name='swind1'):
        self.column = Parameter(name, 'column', 6., 3., 50., 0.0, hugeval)
        self.logxi = Parameter(name, 'logxi', 2.5, 2.1, 4.1, 0.0, hugeval)
        self.sigma = Parameter(name, 'sigma', 0.1, 0., .5, 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.column, self.logxi, self.sigma, self.redshift))


class XSTBabs(XSMultiplicativeModel):

    _calc =  _xspec.C_tbabs

    def __init__(self, name='tbabs'):
        self.nH = Parameter(name, 'nH', 1., 0., 1E5, 0.0, hugeval, '10^22 atoms / cm^2')
        XSMultiplicativeModel.__init__(self, name, (self.nH,))


class XSTBgrain(XSMultiplicativeModel):

    _calc =  _xspec.C_tbgrain

    def __init__(self, name='tbgrain'):
        self.nH = Parameter(name, 'nH', 1., 0., 1E5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.h2 = Parameter(name, 'h2', 0.2, 0., 1., 0.0, hugeval, frozen=True)
        self.rho = Parameter(name, 'rho', 1., 0., 5., 0.0, hugeval, 'g/cm^3', True)
        self.amin = Parameter(name, 'amin', 0.025, 0., 0.25, 0.0, hugeval, 'mum', True)
        self.amax = Parameter(name, 'amax', 0.25, 0., 1., 0.0, hugeval, 'mum', True)
        self.PL = Parameter(name, 'PL', 3.5, 0., 5., 0.0, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.h2, self.rho, self.amin, self.amax, self.PL))


class XSTBvarabs(XSMultiplicativeModel):

    _calc =  _xspec.C_tbvabs

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

    _calc =  _xspec.xsred

    def __init__(self, name='uvred'):
        self.EBV = Parameter(name, 'EBV', 0.05, 0., 10., 0.0, hugeval)
        XSMultiplicativeModel.__init__(self, name, (self.EBV,))


class XSvarabs(XSMultiplicativeModel):

    _calc =  _xspec.xsabsv

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

    _calc =  _xspec.xsvphb

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

    _calc =  _xspec.xsabsw

    def __init__(self, name='wabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        XSMultiplicativeModel.__init__(self, name, (self.nH,))


class XSwndabs(XSMultiplicativeModel):

    _calc =  _xspec.xswnab

    def __init__(self, name='wndabs'):
        self.nH = Parameter(name, 'nH', 1., 0., 10., 0.0, hugeval, '10^22 atoms / cm^2')
        self.WindowE = Parameter(name, 'WindowE', 1., .05, 20., 0.0, hugeval, 'keV')
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.WindowE))


class XSxion(XSMultiplicativeModel):

    _calc =  _xspec.xsxirf

    def __init__(self, name='xion'):
        self.height = Parameter(name, 'height', 5., 0.0, 1.e2, 0.0, hugeval, 'r_s')
        self.lxld = Parameter(name, 'lxld', 0.3, 0.02, 100, 0.0, hugeval)
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
        XSMultiplicativeModel.__init__(self, name, (self.height, self.lxld, self.rate, self.cosAng, self.inner, self.outer, self.index, self.redshift, self.Feabun, self.E_cut, self.Ref_type, self.Rel_smear, self.Geometry))


class XSzdust(XSMultiplicativeModel):

    _calc =  _xspec.mszdst

    def __init__(self, name='zdust'):
        self.method = Parameter(name, 'method', 1, 1, 3, 1, 3, alwaysfrozen=True)
        self.EBV = Parameter(name, 'EBV', 0.1, 0.0, 100., 0.0, hugeval)
        self.Rv = Parameter(name, 'Rv', 3.1, 0.0, 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0.0, 0.0, 20., 0.0, hugeval, 'z', True)
        XSMultiplicativeModel.__init__(self, name, (self.method, self.EBV, self.Rv, self.redshift))


class XSzedge(XSMultiplicativeModel):

    _calc =  _xspec.xszedg

    def __init__(self, name='zedge'):
        self.edgeE = Parameter(name, 'edgeE', 7.0, 0., 100., 0.0, hugeval, 'keV')
        self.MaxTau = Parameter(name, 'MaxTau', 1., 0., 5., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.edgeE, self.MaxTau, self.redshift))


class XSzhighect(XSMultiplicativeModel):

    _calc =  _xspec.xszhcu

    def __init__(self, name='zhighect'):
        self.cutoffE = Parameter(name, 'cutoffE', 10., 1.e-2, 100., 0.0, hugeval, 'keV')
        self.foldE = Parameter(name, 'foldE', 15., 1.e-2, 100., 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.cutoffE, self.foldE, self.redshift))

    def guess(self, dep, *args, **kwargs):
        pos = get_xspec_position(dep, *args)
        param_apply_limits(pos, self.cutoffE, **kwargs)


class XSzpcfabs(XSMultiplicativeModel):

    _calc =  _xspec.xszabp

    def __init__(self, name='zpcfabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.CvrFract = Parameter(name, 'CvrFract', 0.5, 0.05, 0.95, 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.CvrFract, self.redshift))


class XSzphabs(XSMultiplicativeModel):

    _calc = _xspec.xszphb

    def __init__(self, name='zphabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.redshift))


class XSzxipcf(XSMultiplicativeModel):

    _calc =  _xspec.zxipcf

    def __init__(self, name='zxipcf'):
        self.Nh = Parameter(name, 'Nh', 10, 0.05, 500, 0.0, hugeval, '10^22 atoms / cm^2')
        self.logxi = Parameter(name, 'logxi', 3, -3, 6, -hugeval, hugeval)
        self.CvrFract = Parameter(name, 'CvrFract', 0.5, 0., 1., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.Nh, self.logxi, self.CvrFract, self.redshift))


class XSzredden(XSMultiplicativeModel):

    _calc =  _xspec.xszcrd

    def __init__(self, name='zredden'):
        self.EBV = Parameter(name, 'EBV', 0.05, 0., 10., 0.0, hugeval)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.EBV, self.redshift))


class XSzsmdust(XSMultiplicativeModel):

    _calc =  _xspec.msldst

    def __init__(self, name='zsmdust'):
        self.EBV = Parameter(name, 'EBV', 0.1, 0.0, 100., 0.0, hugeval)
        self.ExtIndex = Parameter(name, 'ExtIndex', 1.0, -10.0, 10., -hugeval, hugeval)
        self.Rv = Parameter(name, 'Rv', 3.1, 0.0, 10., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0.0, 0.0, 20., 0.0, hugeval, 'z', True)
        XSMultiplicativeModel.__init__(self, name, (self.EBV, self.ExtIndex, self.Rv, self.redshift))


class XSzTBabs(XSMultiplicativeModel):

    _calc =  _xspec.C_ztbabs

    def __init__(self, name='ztbabs'):
        self.nH = Parameter(name, 'nH', 1., 0., 1E5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.redshift))


class XSzvarabs(XSMultiplicativeModel):

    _calc =  _xspec.xszvab

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

    _calc =  _xspec.xszvfe

    def __init__(self, name='zvfeabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.metals = Parameter(name, 'metals', 1., 0.0, 100., 0.0, hugeval)
        self.FEabun = Parameter(name, 'FEabun', 1., 0.0, 100., 0.0, hugeval)
        self.FEKedge = Parameter(name, 'FEKedge', 7.11, 7.0, 9.5, 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.metals, self.FEabun, self.FEKedge, self.redshift))


class XSzvphabs(XSMultiplicativeModel):

    _calc =  _xspec.xszvph

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

    _calc =  _xspec.xszabs

    def __init__(self, name='zwabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.redshift))


class XSzwndabs(XSMultiplicativeModel):

    _calc =  _xspec.xszwnb

    def __init__(self, name='zwndabs'):
        self.nH = Parameter(name, 'nH', 1., 0., 10., 0.0, hugeval, '10^22 atoms / cm^2')
        self.WindowE = Parameter(name, 'WindowE', 1., .05, 20., 0.0, hugeval, 'keV')
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.WindowE, self.redshift))


### New additive and multiplicative models, as of XSPEC 12.7


class XScplinear(XSAdditiveModel):

    _calc = _xspec.C_cplinear

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
        XSAdditiveModel.__init__(self, name, (self.energy00, self.energy01, self.energy02, self.energy03, self.energy04, self.energy05, self.energy06, self.energy07, self.energy08, self.energy09, self.log_rate00, self.log_rate01, self.log_rate02, self.log_rate03, self.log_rate04, self.log_rate05, self.log_rate06, self.log_rate07, self.log_rate08, self.log_rate09, self.norm))


class XSeqpair(XSAdditiveModel):

    _calc = _xspec.C_xseqpair

    def __init__(self, name='eqpair'):
        self.l_hl_s = Parameter(name, 'l_hl_s', 1., 1e-6, 1.e6, 0.0, hugeval)
        self.l_bb = Parameter(name, 'l_bb', 100., 0., 1.e4, 0.0, hugeval)
        self.kT_bb = Parameter(name, 'kT_bb', 200., 1., 4e5, 0.0, hugeval, 'eV', True)
        self.l_ntl_h = Parameter(name, 'l_ntl_h', 0.5, 0., 0.9999, 0.0, hugeval)
        self.tau_p = Parameter(name, 'tau_p', 0.1, 1e-4, 10., 0.0, hugeval, frozen=True)
        self.radius = Parameter(name, 'radius', 1.e7, 1.e5, 1.e16, 0.0, hugeval, 'cm', True)
        self.g_min = Parameter(name, 'g_min', 1.3, 1.2, 1.e3, 0.0, hugeval, frozen=True)
        self.g_max = Parameter(name, 'g_max', 1.e3, 5., 1.e4, 0.0, hugeval, frozen=True)
        self.G_inj = Parameter(name, 'G_inj', 2., 0., 5., 0.0, hugeval, frozen=True)
        self.pairinj = Parameter(name, 'pairinj', 0., 0., 1., 0.0, hugeval, frozen=True)
        self.cosIncl = Parameter(name, 'cosIncl', 0.50, 0.05, 0.95, 0.0, hugeval, frozen=True)
        self.Refl = Parameter(name, 'Refl', 1., 0., 2., 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.1, 10., 0.0, hugeval, frozen=True)
        self.AbHe = Parameter(name, 'AbHe', 1.0, 0.1, 10., 0.0, hugeval, frozen=True)
        self.T_disk = Parameter(name, 'T_disk', 1.e6, 1e4, 1e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 0.0, 0.0, 1000.0, 0.0, hugeval)
        self.Beta = Parameter(name, 'Beta', -10., -10., 10., -hugeval, hugeval, frozen=True)
        self.Rin = Parameter(name, 'Rin', 10., 6.001, 1.e3, 0.0, hugeval, 'M', True)
        self.Rout = Parameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, hugeval, 'M', True)
        self.redshift = Parameter(name, 'redshift', 0., 0., 4., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.l_hl_s, self.l_bb, self.kT_bb, self.l_ntl_h, self.tau_p, self.radius, self.g_min, self.g_max, self.G_inj, self.pairinj, self.cosIncl, self.Refl, self.Fe_abund, self.AbHe, self.T_disk, self.xi, self.Beta, self.Rin, self.Rout, self.redshift, self.norm))


class XSeqtherm(XSAdditiveModel):

    _calc = _xspec.C_xseqth

    def __init__(self, name='eqtherm'):
        self.l_hl_s = Parameter(name, 'l_hl_s', 1., 1e-6, 1.e6, 0.0, hugeval)
        self.l_bb = Parameter(name, 'l_bb', 100., 0., 1.e4, 0.0, hugeval)
        self.kT_bb = Parameter(name, 'kT_bb', 200., 1., 4e5, 0.0, hugeval, 'eV', True)
        self.l_ntl_h = Parameter(name, 'l_ntl_h', 0.5, 0., 0.9999, 0.0, hugeval)
        self.tau_p = Parameter(name, 'tau_p', 0.1, 1e-4, 10., 0.0, hugeval, frozen=True)
        self.radius = Parameter(name, 'radius', 1.e7, 1.e5, 1.e16, 0.0, hugeval, 'cm', True)
        self.g_min = Parameter(name, 'g_min', 1.3, 1.2, 1.e3, 0.0, hugeval, frozen=True)
        self.g_max = Parameter(name, 'g_max', 1.e3, 5., 1.e4, 0.0, hugeval, frozen=True)
        self.G_inj = Parameter(name, 'G_inj', 2., 0., 5., 0.0, hugeval, frozen=True)
        self.pairinj = Parameter(name, 'pairinj', 0., 0., 1., 0.0, hugeval, frozen=True)
        self.cosIncl = Parameter(name, 'cosIncl', 0.50, 0.05, 0.95, 0.0, hugeval, frozen=True)
        self.Refl = Parameter(name, 'Refl', 1., 0., 2., 0.0, hugeval, frozen=True)
        self.Fe_abund = Parameter(name, 'Fe_abund', 1., 0.1, 10., 0.0, hugeval, frozen=True)
        self.AbHe = Parameter(name, 'AbHe', 1.0, 0.1, 10., 0.0, hugeval, frozen=True)
        self.T_disk = Parameter(name, 'T_disk', 1.e6, 1e4, 1e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 0.0, 0.0, 1000.0, 0.0, hugeval)
        self.Beta = Parameter(name, 'Beta', -10., -10., 10., -hugeval, hugeval, frozen=True)
        self.Rin = Parameter(name, 'Rin', 10., 6.001, 1.e3, 0.0, hugeval, 'M', True)
        self.Rout = Parameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, hugeval, 'M', True)
        self.redshift = Parameter(name, 'redshift', 0., 0., 4., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.l_hl_s, self.l_bb, self.kT_bb, self.l_ntl_h, self.tau_p, self.radius, self.g_min, self.g_max, self.G_inj, self.pairinj, self.cosIncl, self.Refl, self.Fe_abund, self.AbHe, self.T_disk, self.xi, self.Beta, self.Rin, self.Rout, self.redshift, self.norm))


class XScompth(XSAdditiveModel):

    _calc = _xspec.C_xscompth

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
        self.AbHe = Parameter(name, 'AbHe', 1.0, 0.1, 10., 0.0, hugeval, frozen=True)
        self.T_disk = Parameter(name, 'T_disk', 1.e6, 1e4, 1e6, 0.0, hugeval, 'K', True)
        self.xi = Parameter(name, 'xi', 0.0, 0.0, 1000.0, 0.0, hugeval)
        self.Beta = Parameter(name, 'Beta', -10., -10., 10., -hugeval, hugeval, frozen=True)
        self.Rin = Parameter(name, 'Rin', 10., 6.001, 1.e3, 0.0, hugeval, 'M', True)
        self.Rout = Parameter(name, 'Rout', 1.e3, 0., 1.e6, 0.0, hugeval, 'M', True)
        self.redshift = Parameter(name, 'redshift', 0., 0., 4., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.theta, self.showbb, self.kT_bb, self.RefOn, self.tau_p, self.radius, self.g_min, self.g_max, self.G_inj, self.pairinj, self.cosIncl, self.Refl, self.Fe_abund, self.AbHe, self.T_disk, self.xi, self.Beta, self.Rin, self.Rout, self.redshift, self.norm))


class XSbvvapec(XSAdditiveModel):

    _calc = _xspec.xsbvvp

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


class XSvvapec(XSAdditiveModel):

    _calc = _xspec.xsvvap

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

    _calc = _xspec.zigm

    def __init__(self, name='zigm'):
        self.redshift = Parameter(name, 'redshift', 0.0, alwaysfrozen=True)
        self.model = Parameter(name, 'model', 0, 0, 1, 0, 1, alwaysfrozen=True)
        self.lyman_limit = Parameter(name, 'lyman_limit', 1, 0, 1, 0, 1, alwaysfrozen=True)
        XSMultiplicativeModel.__init__(self, name, (self.redshift, self.model, self.lyman_limit))


## Here are the seven new additive models from XSPEC 12.7.1


class XSgadem(XSAdditiveModel):

    _calc = _xspec.C_gaussDem

    def __init__(self, name='gadem'):
        self.Tmean = Parameter(name, 'Tmean', 4.0, 0.01, 10, 0.0, hugeval, 'keV', True)
        self.Tsigma = Parameter(name, 'Tsigma', 0.1, 0.01, 1.e2, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1.0, 1.e-5, 1.e19, 0.0, hugeval, 'cm^-3', True)
        self.abundanc = Parameter(name, 'abundanc', 1.0, 0., 10., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -hugeval, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 2, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Tmean, self.Tsigma, self.nH, self.abundanc, self.Redshift, self.switch, self.norm))


class XSvgadem(XSAdditiveModel):

    _calc = _xspec.C_vgaussDem

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
        XSAdditiveModel.__init__(self, name, (self.Tmean, self.Tsigma, self.nH, self.He, self.C, self.N, self.O, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.Redshift, self.switch, self.norm))


class XSeplogpar(XSAdditiveModel):

    _calc = _xspec.eplogpar

    def __init__(self, name='eplogpar'):
        self.Ep = Parameter(name, 'Ep', .1, 1.e-6, 1.e2, 0.0, hugeval, 'keV')
        self.beta = Parameter(name, 'beta', 0.2, -4., 4., -hugeval, hugeval)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.Ep, self.beta, self.norm))


class XSlogpar(XSAdditiveModel):

    _calc = _xspec.logpar

    def __init__(self, name='logpar'):
        self.alpha = Parameter(name, 'alpha', 1.5, 0., 4., 0.0, hugeval)
        self.beta = Parameter(name, 'beta', 0.2, -4., 4., -hugeval, hugeval)
        self.pivotE = Parameter(name, 'pivotE', 1.0, units='keV', alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.alpha, self.beta, self.pivotE, self.norm))


class XSoptxagn(XSAdditiveModel):

    _calc = _xspec.optxagn

    def __init__(self, name='optxagn'):
        self.mass = Parameter(name, 'mass', 1e7, 1.0, 1.e9, 0.0, hugeval, 'solar', True)
        self.dist = Parameter(name, 'dist', 100, 0.01, 1.e9, 0.0, hugeval, 'Mpc', True)
        self.logLLEdd = Parameter(name, 'logLLEdd', -1., -10., 2, -hugeval, hugeval)
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
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval, frozen=True)
        XSAdditiveModel.__init__(self, name, (self.mass, self.dist, self.logLLEdd, self.astar, self.rcor, self.logrout, self.kT_e, self.tau, self.Gamma, self.fpl, self.fcol, self.tscat, self.Redshift, self.norm))


class XSoptxagnf(XSAdditiveModel):

    _calc = _xspec.optxagnf

    def __init__(self, name='optxagnf'):
        self.mass = Parameter(name, 'mass', 1e7, 1.0, 1.e9, 0.0, hugeval, 'solar', True)
        self.dist = Parameter(name, 'dist', 100, 0.01, 1.e9, 0.0, hugeval, 'Mpc', True)
        self.logLLEdd = Parameter(name, 'logLLEdd', -1., -10., 2, -hugeval, hugeval)
        self.astar = Parameter(name, 'astar', 0.0, 0., 0.998, 0.0, hugeval, frozen=True)
        self.rcor = Parameter(name, 'rcor', 10.0, 1., 100., 0.0, hugeval, 'rg')
        self.logrout = Parameter(name, 'logrout', 5.0, 3.0, 7.0, 0.0, hugeval, frozen=True)
        self.kT_e = Parameter(name, 'kT_e', 0.2, 0.01, 10, 0.0, hugeval, 'keV')
        self.tau = Parameter(name, 'tau', 10., 0.1, 100, 0.0, hugeval)
        self.Gamma = Parameter(name, 'Gamma', 2.1, 0.5, 5., 0.0, hugeval)
        self.fpl = Parameter(name, 'fpl', 1.e-4, 0.0, 1., 0.0, hugeval)
        self.Redshift = Parameter(name, 'Redshift', 0., 0., 10., 0.0, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval, frozen=True)
        XSAdditiveModel.__init__(self, name, (self.mass, self.dist, self.logLLEdd, self.astar, self.rcor, self.logrout, self.kT_e, self.tau, self.Gamma, self.fpl, self.Redshift, self.norm))


class XSpexmon(XSAdditiveModel):

    _calc = _xspec.pexmon

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

    _calc =  _xspec.C_agauss

    def __init__(self, name='agauss'):
        self.LineE = Parameter(name, 'LineE', 10.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.Sigma = Parameter(name, 'Sigma', 1.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Sigma, self.norm))

class XSzagauss(XSAdditiveModel):

    _calc =  _xspec.C_zagauss

    def __init__(self, name='zagauss'):
        self.LineE = Parameter(name, 'LineE', 10.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.Sigma = Parameter(name, 'Sigma', 1.0, 0.0, 1.0e6, 0.0, 1.0e6, units='A')
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10.0, -0.999, 10.0, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.LineE, self.Sigma, self.Redshift, self.norm))

class XScompmag(XSAdditiveModel):

    _calc =  _xspec.xscompmag

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

    _calc =  _xspec.xscomptb

    def __init__(self, name='comptb'):
        self.kTs = Parameter(name, 'kTs', 1.0, 0.1, 10.0, 0.1, 10.0, units='keV')
        self.gamma = Parameter(name, 'gamma', 3.0, 1.0, 10.0, 1.0, 10.0)
        self.alpha = Parameter(name, 'alpha', 2.0, 0.0, 400.0, 0.0, 400.0)
        self.delta = Parameter(name, 'delta', 20.0, 0.0, 200.0, 0.0, 200.0)
        self.kTe = Parameter(name, 'kTe', 5.0, 0.2, 2000.0, 0.2, 2000.0, units='keV')
        self.log_A = Parameter(name, 'log_A', 0.0, -8.0, 8.0, -8.0, 8.0)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        XSAdditiveModel.__init__(self, name, (self.kTs, self.gamma, self.alpha, self.delta, self.kTe, self.log_A, self.norm))


class XSheilin(XSMultiplicativeModel):

    _calc =  _xspec.xsphei

    def __init__(self, name='heilin'):
        self.nHeI = Parameter(name, 'nHeI', 1.e-5, 0.0, 1.e6, 0.0, 1.0e6, '10^22 atoms / cm^2')
        self.b = Parameter(name, 'b', 10.0, 1.0, 1.0e5, 1.0, 1.0e6, units='km/s')
        self.redshift = Parameter(name, 'redshift', 0.0, -1.0e-3, 1.0e5, -1.0e-3, 1.0e5)
        XSMultiplicativeModel.__init__(self, name, (self.nHei, self.b, self.redshift))

class XSlyman(XSMultiplicativeModel):

    _calc =  _xspec.xslyman

    def __init__(self, name='lyman'):
        self.nHeI = Parameter(name, 'nHeI', 1.e-5, 0.0, 1.0e6, 0.0, 1.0e6, '10^22 atoms / cm^2')
        self.b = Parameter(name, 'b', 10.0, 1.0, 1.0e5, 1.0, 1.0e6, units='km/s')
        self.redshift = Parameter(name, 'redshift', 0.0, -1.0e-3, 1.0e5, -1.0e-3, 1.0e5)
        self.ZA = Parameter(name, 'ZA', 1.0, 1.0, 2.0, 1.0, 2.0)
        XSMultiplicativeModel.__init__(self, name, (self.nHeI, self.b, self.redshift, self.ZA))

class XSzbabs(XSMultiplicativeModel):

    _calc =  _xspec.xszbabs

    def __init__(self, name='zbabs'):
        self.nH = Parameter(name, 'nH', 1.e-4, 0.0, 1.0e5, 0.0, 1.0e6, '10^22 atoms / cm^2')
        self.nHeI = Parameter(name, 'nHeI', 1.e-5, 0.0, 1.0e5, 0.0, 1.0e6, '10^22 atoms / cm^2')
        self.nHeII = Parameter(name, 'nHeII', 1.e-6, 0.0, 1.0e5, 0.0, 1.0e6, '10^22 atoms / cm^2')
        self.redshift = Parameter(name, 'redshift', 0.0, 0.0, 1.0e5, 0.0, 1.0e6)
        XSMultiplicativeModel.__init__(self, name, (self.nH, self.nHeI, self.nHeII, self.redshift))

# Add model classes to __all__
__all__ += tuple(n for n in globals() if n.startswith('XS'))
