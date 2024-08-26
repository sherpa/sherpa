#!/usr/bin/env python

#
#  Copyright (C) 2021, 2022, 2024
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

"""Usage:

  ./check_xspec_update.py infile

    --nohard
    --noswitch

Aim:

Given an XSPEC model.dat file report on differences to the supported
XSPEC models in Sherpa. This is to indicate what changes are needed to
update the XSPEC module to support a new XSPEC module (but it does not
handle things like conditional support for a model, or changes to how
the models are built). In other words, it simplifies the process of
supporting a new XSPEC release but does not do everything!

At the moment this requires a Sherpa build with a working XSPEC
module; this could be worked around but it's not felt to be worth it
at this time.

It could take the infile from $HEADAS/../spectral/manager/model.dat
if not set, but I want to be explicit about what file is being used,
and it's too easy to have the environment set for a different version
of XSPEC.

"""

import sys
from typing import Sequence

from sherpa.astro.utils.xspec import ModelDefinition, \
    parse_xspec_model_description
from sherpa.astro import xspec
from sherpa.models.parameter import hugeval


def compare_xspec_models(models: Sequence[ModelDefinition],
                         hard: bool = True,
                         switch: bool = True
                         ) -> None:
    """Check the sherpa.astro.xspec models to those in models.

    This is intended to see what updates are needed when updating the
    XSPEC version supported in Sherpa. The output points out
    differences but does not do a detailed comparison.

    Parameters
    ----------
    models : list of ModelDefinition
    hard : bool, optional
        Do we bother checking the hard limits (default True)?
    switch : bool, optional
        Do we bother checking the limits for the switch parameter
        (default True)?

    """

    # Check if we know about all the models in models
    #
    seen = set()
    for mdl in models:
        try:
            xscls = getattr(xspec, mdl.clname)
        except AttributeError:
            print(f"We do not support {mdl.name} ({mdl.modeltype}; {mdl.language}; {mdl.funcname})\n")
            continue

        seen.add(mdl.clname)
        xs = xscls()

        if mdl.modeltype == 'Add' and not isinstance(xs, xspec.XSAdditiveModel):
            print(f"Model {mdl.name} has switched to additive\n")
            continue

        if mdl.modeltype == 'Mul' and not isinstance(xs, xspec.XSMultiplicativeModel):
            print(f"Model {mdl.name} has switched to multiplicative\n")
            continue

        if mdl.modeltype == 'Con' and not isinstance(xs, xspec.XSConvolutionKernel):
            print(f"Model {mdl.name} has switched to convolution\n")
            continue

        if len(xs.pars) != len(mdl.pars):
            print(f"Model {mdl.name} parameters: {len(xs.pars)} -> {len(mdl.pars)}\n")
            continue

        # The previous changes are deemed serious enough that we need
        # to look at the whole model. Now we are concerned with
        # smaller changes (so can have multiple ones for a model).
        #
        reports = []

        # parse_xspec_user_model has removed the "language-type" flags
        # on the function name, so we need to restore them.
        #
        funcname = mdl.funcname
        if mdl.language == 'C++ style':
            funcname = f'C_{funcname}'
        elif mdl.language == "Fortran - double precision":
            funcname = f"F_{funcname}"

        if xs._calc.__name__ != funcname:
            reports.append(f"function name change: {xs._calc.__name__} to {funcname}")

        # VERY LIMITED CHECK OF PARAMETER VALUES
        #
        # We have two settings,the soft/hard limits (which are, as of
        # #1259, the same) and the original values we set, which are
        # stored in the _xpsec_xoft/hard_min/max attributes).
        #
        for i, (xpar, par) in enumerate(zip(xs.pars, mdl.pars), 1):
            if xpar.name != par.name:
                reports.append(f"par {i} name: {xpar.name} -> {par.name}")

            units = '' if par.units is None else par.units
            if xpar.units != units:
                reports.append(f"par {xpar.name} units: {xpar.units} -> {units}")

            if xpar.val != par.default:
                reports.append(f"par {xpar.name} default: {xpar.val} -> {par.default}")

            # We skip the Parameter instances which are the norm parameters
            #
            if isinstance(xpar, xspec.XSBaseParameter):
                # How do these values compare to the values from the old
                # model.dat file?
                #
                # First check the frozen report, as the remaining checks are
                # for min/max ranges. Note that we don't always have a frozen
                # attribute for a parameter.
                #
                try:
                    frozen = par.frozen
                except AttributeError:
                    frozen = True

                if xpar.frozen != frozen:
                    reports.append(f"par {xpar.name} frozen: {xpar.frozen} -> {frozen}")


                # The actual soft limits should match the model.dat hard limits
                # (and should also be reported in the _xspec_soft... checks).
                # We do special case the situation where the par limits are None
                # (i.e. undefined) as long as the xpar values are +/- hugeval
                #
                # For the moment this does not care about the switch flag
                #
                # Do we want to skip the min/max checks? For now we only
                # skip the switch settings if they are all (in the model.dat
                # file) not set.
                #
                if [par.softmin, par.hardmin, par.softmax, par.hardmax] == [None] * 4:

                    if not switch and xpar.name == "switch":
                        continue

                    for attr in ["min", "hard_min"]:
                        got = getattr(xpar, attr)
                        if got != -hugeval:
                            reports.append(f"par {xpar.name}.{attr} is not -hugeval but {got}")

                    for attr in ["max", "hard_max"]:
                        got = getattr(xpar, attr)
                        if got != hugeval:
                            reports.append(f"par {xpar.name}.{attr} is not hugeval but {got}")

                    continue

                if xpar.min != par.hardmin:
                    reports.append(f"par {xpar.name} min: {xpar.min} -> {par.hardmin}")

                if xpar.max != par.hardmax:
                    reports.append(f"par {xpar.name} max: {xpar.max} -> {par.hardmax}")

                if xpar._xspec_soft_min != par.softmin:
                    reports.append(f"par {xpar.name} softmin: {xpar._xspec_soft_min} -> {par.softmin}")

                if xpar._xspec_soft_max != par.softmax:
                    reports.append(f"par {xpar.name} softmax: {xpar._xspec_soft_max} -> {par.softmax}")

                if hard:
                    if xpar.hard_min != par.hardmin:
                        reports.append(f"par {xpar.name} hardmin: {xpar.hard_min} -> {par.hardmin}")

                    if xpar.hard_max != par.hardmax:
                        reports.append(f"par {xpar.name} hardmax: {xpar.hard_max} -> {par.hardmax}")

            elif xpar.name != 'norm':
                raise ValueError(f"Unexpected parameter for {mdl.clname}\n{xpar}")

        if len(reports) == 0:
            continue

        print(f"Model {mdl.name}")
        for i, report in enumerate(reports, 1):
            print(f"  [{i:2d}] {report}")

        print("")

    # Check to see if we have lost any models.
    #
    for name in dir(xspec):
        if not name.startswith('XS') or name.endswith('Model') or \
           name.endswith('Kernel') or name.endswith('Parameter'):
            continue

        try:
            xs = getattr(xspec, name)
        except TypeError:
            continue

        if name not in seen:
            print(f"Unable to find {name} in model.dat")


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    # The argument names don't fit the code (i.e. they'd be better
    # inverted) but this is better documentation for the user.
    #
    parser.add_argument("infile", type=str,
                        help="The spectral/manager/model.dat file")
    parser.add_argument("--nohard", action="store_true",
                        help="do not check hard ranges")
    parser.add_argument("--noswitch", action="store_true",
                        help="do not check limits of switch parameter")

    args = parser.parse_args()

    # This errors out in case of significant model-support issues
    # (e.g. a periodic model parameter) but can also just be an
    # issue with a poorly-specified interchange format (the model.dat
    # file) which may just need changes to this routine.
    #
    models = parse_xspec_model_description(args.infile)
    compare_xspec_models(models,
                         hard=not args.nohard,
                         switch=not args.noswitch)
