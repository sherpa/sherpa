#
#  Copyright (C) 2023, 2024
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

"""Create XSPEC table models.

This module supports creating and writing out `XSPEC table models
<https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html>`_,
which can then be read in by the
`sherpa.astro.xspec.read_xstable_model` and
`sherpa.astro.ui.load_xstable_model` routines.

These routines should be considered experimental, as it is not
obvious that the interface is as usable as it could be.

.. versionchanged:: 4.17.0

   The FITS-like data types such as HeaderItem, Column, and
   TableBlock, are now defined in the sherpa.astro.io.types module.
   Some types have been changed as part of this work.

Example
-------

The following example creates an additive table model
that represents a gaussian model where the only parameters are
the line position and the normalization.

>>> import numpy as np
>>> from sherpa.astro.io import xstable
>>> from sherpa.astro.xspec import XSgaussian

We shall use the XSPEC gaussian model `sherpa.astro.xspec.XSgaussian`
to create the template models.

>>> mdl = XSgaussian()
>>> print(mdl)
gaussian
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   gaussian.LineE thawed          6.5            0        1e+06        keV
   gaussian.Sigma thawed          0.1            0           20        keV
   gaussian.norm thawed            1            0        1e+24

The models will be evaluated over the 0.1 to 2 keV range, with
a bin spacing of 0.01 keV (note that the last bin ends at 1.99 and
not 2.0 keV thanks to the behavior of `numpy.arange`).

>>> egrid = np.arange(0.1, 2, 0.01)
>>> elo = egrid[:-1]
>>> ehi = egrid[1:]
>>> emid = (elo + ehi) / 2

The table model will interpolate over the line position over
the range of 0 to 2.4 inclusive, with a spacing of 0.1 keV:

>>> linepos = np.arange(0, 2.5, 0.1)
>>> minval = linepos[0]
>>> maxval = linepos[-1]

The model is evaluated over the grid defined by ``elo`` and ``ehi``
for each element of ``linepos``, which is used to set the
`sherpa.astro.xspec.XSgaussian.LineE` parameter:

>>> models = []
>>> for lp in linepos:
...     mdl.linee = lp
...     ymodel = mdl(elo, ehi)
...     models.append(ymodel)
...

This model has a single interpolated parameter which we call
``"pos"``.  Note that the delta parameter - here set to 0.01 - is only
used by Sherpa to decide if the parameter is frozen (value is
negative) or not, but is used by the optimiser in XSPEC. The values
field is set to those used to generate the spectral models.

>>> param = xstable.Param("pos", 1, 0.01, minval, maxval,
...                       values=linepos)

With this we can create the necessary information to create the
table model (`make_xstable_model`) and then write it out
as a FITS file (`write_xstable_model`):

>>> hdus = xstable.make_xstable_model("gaussy", elo, ehi, params=[param],
...                                   spectra=models)
>>> xstable.write_xstable_model("example.mod", hdus, clobber=True)

This file can then be used in Sherpa with either
`sherpa.astro.ui.load_xstable_model` or
`sherpa.astro.xspec.read_xstable_model`. For example:

>>> from sherpa.astro import ui
>>> ui.load_xstable_model("mygauss", "example.mod")
>>> print(mygauss)
xstablemodel.mygauss
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   mygauss.pos  thawed            1            0          2.4
   mygauss.norm thawed            1            0        1e+24

Notes
-----

XSPEC models can be created without XSPEC support in Sherpa, but XSPEC
support is needed to read the files in.

For additive models it is assumed that the model values - that is,
each bin - have units of photon/cm^2/s. This is easy to accidentally
change - e.g.  in the example above if mdl.norm were changed to a
value other than 1 everything would work but Sherpa and XSPEC would
infer incorrect fluxes or luminosities.

"""

from dataclasses import dataclass, field
import itertools
import importlib.metadata
import time
from typing import Optional, Sequence, Union

import numpy as np

import sherpa

from .types import HeaderItem, Header, Column, TableBlock, BlockList


__all__ = ("BaseParam", "Param", "make_xstable_model",
           "write_xstable_model")


@dataclass
class BaseParam:
    """Represent a parameter.

    The absolute DELTA value is used by XSPEC but ignored by Sherpa.
    The sign of DELTA is used (by both Sherpa and XSPEC) to decide if
    the parameter should be frozen (DELTA is negative) or not by
    default.

    """

    name: str
    """The parameter name.

    This must be unique and ideally a valid attribute name in Python,
    although the latter is not enforced.
    """

    initial: float
    """The default value for the parameter."""

    delta: float
    """The delta value for fitting.

    If negative then the parameter defaults to being frozen. Sherpa
    does not use this value (other than checking the sign) but XSPEC
    does.
    """

    hardmin: float
    """The minimum value for the parameter."""

    hardmax: float
    """The maximum value for the parameter."""

    softmin: Optional[float] = None
    """The "soft" minimum. Defaults to hardmin if not given.."""

    softmax: Optional[float] = None
    """The "soft" maximum. Defaults to hardmax if not given."""

    def __post_init__(self):
        """Validate the settings"""

        # Do we need to set up softmin/max?
        #
        if self.softmax is None:
            self.softmax = self.hardmax

        if self.softmin is None:
            self.softmin = self.hardmin

        # There should probably be a check that the name is not
        # an invalid Python attribute name.
        #
        if self.name.strip() == "":
            raise ValueError("name can not be empty")

        if self.softmin < self.hardmin:
            raise ValueError(f"Parameter {self.name} softmin={self.softmin} < hardmin={self.hardmin}")

        if self.softmax > self.hardmax:
            raise ValueError(f"Parameter {self.name} softmax={self.softmax} > hardmax={self.hardmax}")

        # We allow min == max
        if self.softmin > self.softmax:
            raise ValueError(f"Parameter {self.name} softmin={self.softmin} > softmax={self.softmax}")

        if self.initial < self.softmin or self.initial > self.softmax:
            raise ValueError(f"Parameter {self.name} initial={self.initial} outside softmin={self.softmin} to softmax={self.softmax}")

        if self.delta == 0:
            raise ValueError("delta can not be 0")


@dataclass
class Param(BaseParam):
    """Represent an interpolated parameter."""

    values: list[float] = field(default_factory=list)
    """The parameter values used to create the model spectra.

    Is is required that they range from hardmin to hardmax and are in
    monotonically increasing order (so it can not be an empty list).
    """

    loginterp: bool = False
    """Are the values logarithmically interpolated (True) or linearly (False).

    Defaults to linear interpolation.
    """

    def __post_init__(self):
        """Validate the settings"""
        super().__post_init__()

        if len(self.values) == 0:
            raise ValueError(f"Parameter {self.name} has no values")

        if self.values[0] != self.hardmin:
            raise ValueError(f"Parameter {self.name} hardmin={self.hardmin} != first value={self.values[0]}")

        if self.values[-1] != self.hardmax:
            raise ValueError(f"Parameter {self.name} hardmax={self.hardmax} != last value={self.values[-1]}")

        if np.any(np.diff(self.values) <= 0):
            raise ValueError(f"Parameter {self.name} values are not monotonically increasing")


def key(name: str,
        value: Union[str, bool, int, float],
        desc: Optional[str] = None,
        unit: Optional[str] = None) -> HeaderItem:
    "Easily make a HeaderItem"
    return HeaderItem(name=name, value=value, desc=desc, unit=unit)


def col(name: str,
        values: np.ndarray,
        desc: Optional[str] = None,
        unit: Optional[str] = None) -> Column:
    "Easily make a column"
    return Column(name=name, values=values, desc=desc, unit=unit)


def xstable_primary(name: str,
                    unit: str,
                    addmodel: bool,
                    redshift: bool,
                    escale: bool,
                    lolim: float,
                    hilim: float,
                    xfxp: Optional[Sequence[str]] = None) -> Header:
    """The PRIMARY block for an xspec table model."""

    header = [
        key("HDUCLASS", "OGIP",
            "format conforms to OGIP standard"),
        key("HDUCLAS1", "XSPEC TABLE MODEL",
            "model spectra for XSPEC"),
        key("HDUVERS", "1.2.0",
            "version of format"),
        key("MODlNAME", name, "Model name"),
        key("MODLUNIT", unit, "Model units"),
        key("ADDMODEL", addmodel,
            "If T then this is an additive table model"),
        key("REDSHIFT", redshift,
            "If T then redshift will be a parameter"),
        key("ESCALE", escale,
            "If T then escale will be a parameter"),
        key("LOELIMIT", lolim,
            "Model value for energies below ENERGIES block"),
        key("HIELIMIT", hilim,
            "Model value for energies above ENERGIES block")
    ]

    if xfxp:
        header.append(key("NXFLTEXP", len(xfxp),
                          "Number of model spectra per grid point"))
        for idx, val in enumerate(xfxp, 1):
            header.append(key(f"XFXP{idx:04d}", val,
                              f"The expression for model {idx}"))

    # Let users know who created the file and when.
    #
    header.append(key("CREATOR",
                      f"sherpa {importlib.metadata.version('sherpa')}",
                      "Code that created the file"))
    header.append(key("DATE",
                      time.strftime("%Y-%m-%dT%H:%M:%S"),
                      "Date file created"))

    return Header(header)


def xstable_parameters(params: list[Param],
                       addparams: Optional[list[BaseParam]]) -> TableBlock:
    """The PARAMETERS block for an xspec table model."""

    nint = len(params)
    nadd = 0 if addparams is None else len(addparams)
    ntotal = nint + nadd
    header = [
        key("HDUCLASS", "OGIP",
            "format conforms to OGIP standard"),
        key("HDUCLAS1", "XSPEC TABLE MODEL",
            "model spectra for XSPEC"),
        key("HDUCLAS2", "PARAMETERS",
            "extension containing parameter info"),
        key("HDUVERS", "1.0.0", "version of format"),
        key("NINTPARM", nint,
            "Number of interpolation parameters"),
        key("NADDPARM", nadd,
            "Number of additional parameters")
    ]

    # Provide the parameters and then the additional parameters.
    #
    def get(fld):
        f1 = [getattr(p, fld) for p in params]
        if addparams is None:
            return f1

        return f1 + [getattr(p, fld) for p in addparams]

    methods = [1 if p.loginterp else 0 for p in params] + [0] * nadd
    nvalues = [len(p.values) for p in params] + [0] * nadd

    # Should the VALUE field use a variable-length field or not?
    #
    max_nvalues = max(nvalues)
    if nint > 1 and max_nvalues > 3:
        # The main array stores objects as each row can have different
        # lengths.
        #
        values = np.zeros(ntotal, dtype=object)
        for idx, p in enumerate(params):
            values[idx] = np.asarray(p.values, dtype=np.float32)

        # Empty arrays for the additional arrays
        for idx in range(nint, ntotal):
            values[idx] = np.asarray([], dtype=np.float32)

    else:
        values = np.zeros((ntotal, max_nvalues), dtype=np.float32)
        for idx, p in enumerate(params):
            values[idx, :len(p.values)] = p.values

    cols = [
        col("NAME", np.asarray(get("name"), dtype="U12"),
            "name of the parameter"),
        col("METHOD", np.asarray(methods, dtype=np.int32),
            "0 if linear, 1 if logarithmic"),
        col("INITIAL", np.asarray(get("initial"), dtype=np.float32),
            "Starting point"),
        col("DELTA", np.asarray(get("delta"), dtype=np.float32),
            "If negative frozen by default"),
        col("MINIMUM", np.asarray(get("hardmin"), dtype=np.float32),
            "hard lower limit"),
        col("BOTTOM", np.asarray(get("softmin"), dtype=np.float32),
            "soft lower limit"),
        col("TOP", np.asarray(get("softmax"), dtype=np.float32),
            "soft upper limit"),
        col("MAXIMUM", np.asarray(get("hardmax"), dtype=np.float32),
            "hard upper limit"),
        col("NUMBVALS", np.asarray(nvalues, dtype=np.int32),
            "number of tabulated parameter values"),
        col("VALUE", values, "tabulated parameter values")
    ]
    return TableBlock(name="PARAMETERS", header=Header(header),
                      columns=cols)


def xstable_energies(energ_lo: np.ndarray,
                     energ_hi: np.ndarray) -> TableBlock:
    """The ENERGIES block for an xspec table model."""

    header = [
        key("HDUCLASS", "OGIP",
            "format conforms to OGIP standard"),
        key("HDUCLAS1", "XSPEC TABLE MODEL",
            "model spectra for XSPEC"),
        key("HDUCLAS2", "ENERGIES",
            "extension containing energy bin info"),
        key("HDUVERS", "1.0.0", "version of format")
    ]
    cols = [
        col("ENERG_LO", np.asarray(energ_lo, dtype=np.float32),
            "Minimum energy of the bin", unit="keV"),
        col("ENERG_HI", np.asarray(energ_hi, dtype=np.float32),
            "Maximum energy of the bin", unit="keV"),
    ]
    return TableBlock(name="ENERGIES", header=Header(header),
                      columns=cols)


def xstable_spectra(paramvals: list[tuple[float, ...]],
                    spectra: list[np.ndarray],
                    addparam: Optional[list[BaseParam]],
                    addspectra: Optional[list[list[np.ndarray]]],
                    units: Optional[str]) -> TableBlock:
    """The SPECTRA block for an xspec table model."""

    header = [
        key("HDUCLASS", "OGIP",
            "format conforms to OGIP standard"),
        key("HDUCLAS1", "XSPEC TABLE MODEL",
            "model spectra for XSPEC"),
        key("HDUCLAS2", "MODEL SPECTRA",
            "extension containing model spectra"),
        key("HDUVERS", "1.0.0",
            "version of format")
    ]

    # The sizes have already been checked.
    #
    unit_opt = None if units == "" else units
    cols = [
        col("PARAMVAL", np.asarray(paramvals, dtype=np.float32),
            "Parameter values for spectrum"),
        col("INTPSPEC", np.asarray(spectra, dtype=np.float32),
            "Spectrum for interpolated parameters",
            unit=unit_opt)
    ]

    # Do we have any additional parameters?
    #
    if addparam is not None and addspectra is not None:
        zipped = zip(addparam, addspectra)
        for idx, (param, aspec) in enumerate(zipped, 1):
            cols.append(col(f"ADDSP{idx:03d}",
                            np.asarray(aspec, dtype=np.float32),
                            desc=f"Additional spectrum {idx:03d}: {param.name}",
                            unit=unit_opt))

    return TableBlock(name="SPECTRA", header=Header(header),
                      columns=cols)


# We could either send in the parameter values broken up by parameter,
# or as a list of those used to create spectra (either way, we need to
# deconstruct/reconstruct the other form).  For now we send them in
# per-parameter, but it could be changed if this approach is found to
# be more awkward for users.
#
def make_xstable_model(name: str,
                       egrid_lo: np.ndarray,
                       egrid_hi: np.ndarray,
                       params: list[Param],
                       spectra: list[np.ndarray],
                       addparams: Optional[list[BaseParam]] = None,
                       addspectra: Optional[list[list[np.ndarray]]] = None,
                       addmodel: bool = True,
                       redshift: bool = False,
                       escale: bool = False,
                       lolim: float = 0,
                       hilim: float = 0,
                       units: Optional[str] = None,
                       xfxp: Optional[Sequence[str]] = None
                       ) -> BlockList:
    """Create the blocks for a XSPEC table model.

    An XSPEC table model can be either additive or multiplicative,
    contain 1 or more interpolated parameters, 0 or more additional
    parameters, and support additional parameters (redshift and
    escale), as well as setting the value to use outside the tabulated
    energy range. The table defines the energy grid the models are
    defined on, and the models used (one per set of parameters unless
    the xfxp argument is set).

    Parameters
    ----------
    name : str
       The model name (stored in the table). Any spaces are removed.
    egrid_lo, egrid_hi : sequence of float
       The energy grid (in keV) used to define the spectra. It is
       required that len(egrid_lo) == len(egrid_hi), the arrays are in
       increasing order, that they are consecutive so that egrid_hi[i]
       == egrid_lo[i + 1], and that egrid_lo[0] > 0.
    params : sequence of Param
       The definitions of the parameters used to interpolate the
       models. It must have at least one element.
    spectra : sequence of sequence of float
       The spectra for each set of parameters in paramvals. It
       is a 2D array of shape(nrows, len(egrid_lo)), and each row
       matches the corresponding parameter grouping:
       paramvals[0].values[i], paramvals[1].values[j], ..
       where the first parameter loops the slowest. The number of rows
       is the multiplication of the number of parameter values.
    addparams : sequence of BaseParam or None, optional
       The definitions of the additional parameters; that is those
       that are not interpolated over.
    addspectra : sequence of sequence of sequence of float or None, optional
       The spectra for the additional parameters. It must is a 3D
       array of shape (nrows, len(addparams), len(egrid_lo)).
    addmodel : bool, optional
       Is this an additive model (True) or a multiplicative one (False).
       The default is True (additive).
    redshift : bool, optional
       Should the redshift parameter be added? The default is False.
    escale : bool, optional
       Should the escale parameter be added? The default is False.
    lolim, hilim : bool, optional
       The value to be used when energies are below or above the values
       in egrid_lo or egrid_hi respectively.
    units : str or None, optional
       The model units. For additive models the default is
       "photons/cm^2/s" and for multiplicative models it is "".
    xfxp : sequence of str or None
       If there are multiple spectra per parameter grid point, these
       values give the expression to use (the number of elements is
       the number of spectra per grid point). It is expected this is
       None or has a length more than 1.

    Returns
    -------
    hdus : BlockList
       The header and column data needed to create a FITS file.

    See Also
    --------
    write_xstable_model

    Notes
    -----

    This supports version 1.2.0 of the `XSPEC table model
    specification
    <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_92_009/ogip_92_009.html>`_,
    although this code should be considered experimental.

    It is not designed to be memory efficient, so for models that use
    a large number of spectra, other tools may make more sense.

    """

    nxfxp = 1 if xfxp is None else len(xfxp)

    # Not sure what the restrictions are, so just remove any spaces.
    #
    name_str = name.translate({32: None})

    if units is None:
        unit_str = "photons/cm^2/s" if addmodel else ""
    else:
        unit_str = units

    # Check the parameter ranges make sense.
    #
    if len(params) == 0:
        raise ValueError("params can not be empty")

    # Check the parameter names are unique. We do a case-insensitive
    # comparison, but we write out the user-specified name.
    #
    names = [p.name.lower() for p in params]
    if addparams is not None:
        names += [p.name.lower() for p in addparams]

    snames = set(names)
    if len(names) != len(snames):
        raise ValueError("Parameter names are not unique")

    # Check the additional parameters
    #
    nadd = 0 if addparams is None else len(addparams)
    naddspec = 0 if addspectra is None else len(addspectra)
    if nadd != naddspec:
        raise ValueError(f"Mismatch between addparams and addspectra sizes: {nadd} {naddspec}")

    # What are the parameter combinations? The last parameter
    # loops fastest. If XFXP keyword are given then we have to
    # replicate each set of parameters nxfxp times.
    #
    pvals_indiv = [p.values for p in params]
    if nxfxp > 1:
        pvals_indiv.append([0] * nxfxp)

    pvalues = list(itertools.product(*pvals_indiv))
    if nxfxp > 1:
        # remove the fake parameter used to duplicate each row
        pvalues = [pv[:-1] for pv in pvalues]

    # Does pvalues match the expected number of spectra?
    #
    npvs = len(pvalues)
    nspectra = len(spectra)
    if npvs != nspectra:
        raise ValueError(f"Expected {npvs} spectra, found {nspectra}")

    if addspectra is not None:
        for idx, aspec in enumerate(addspectra):
            naspec = len(aspec)
            if npvs != naspec:
                # We should use a type that combines
                # addparams/spectra but for now make mypy happy.
                #
                assert addparams is not None
                raise ValueError(f"Expected {npvs} spectra for additional parameter {addparams[idx].name}, found {naspec}")

    # Check the energy grid makes sense
    nlo = len(egrid_lo)
    nhi = len(egrid_hi)
    if nlo != nhi:
        raise ValueError(f"egrid_lo/hi do not match: {nlo} vs {nhi}")

    if nlo == 0:
        raise ValueError("egrid_lo/hi can not be empty")

    if np.any(np.diff(egrid_lo) <= 0):
        raise ValueError("egrid_lo is not monotonically increasing")

    if np.any(np.diff(egrid_hi) <= 0):
        raise ValueError("egrid_hi is not monotonically increasing")

    # Technically we should do a Knuth-like numerical comparison here,
    # but the assumption is that we created them so that they are
    # consecutive, so it's okay for this check (and to try and stop
    # yet-more-Xray-files-from-being-not-quite-valid(TM)).
    #
    if np.any(egrid_lo[1:] != egrid_hi[:-1]):
        raise ValueError("egrid_lo/hi are not consecutive")

    # Check we have the correct number of values in each spectrum.
    #
    for sp in spectra:
        nsp = len(sp)
        if nsp != nlo:
            raise ValueError(f"Spectrum should have {nlo} elements "
                             f"but has {nsp}")

    if addspectra is not None:
        for idx, aspecs in enumerate(addspectra):
            for asp in aspecs:
                nsp = len(asp)
                if nsp != nlo:
                    # We should use a type that combines
                    # addparams/spectra but for now make mypy happy.
                    #
                    assert addparams is not None
                    raise ValueError(f"Spectrum for parameter "
                                     f"{addparams[idx].name} should "
                                     f"have {nlo} elements but has {nsp}")

    # Create the separate blocks.
    #
    primary = xstable_primary(name_str, unit_str, bool(addmodel),
                              bool(redshift), bool(escale),
                              float(lolim), float(hilim),
                              xfxp=xfxp)

    out = [xstable_parameters(params, addparams),
           xstable_energies(egrid_lo, egrid_hi),
           xstable_spectra(pvalues, spectra, addparams, addspectra,
                           units=unit_str)
           ]
    return BlockList(blocks=list(out), header=primary)


def write_xstable_model(filename: str,
                        hdus: BlockList,
                        clobber: bool = False) -> None:
    """Write a XSPEC table model to disk.

    Parameters
    ----------
    filename : str
        The filename to create.
    hdus : BlockList
        The output of make_xstable_model.
    clobber : bool, optional
       If `True` then the output file will be over-written if it
       already exists, otherwise an error is raised.

    See Also
    --------
    make_xstable_model

    """

    from sherpa.astro import io
    io.backend.set_hdus(filename, hdus, clobber=clobber)
