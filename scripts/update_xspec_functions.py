#!/usr/bin/env python

#
#  Copyright (C) 2024
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

  ./update_spec_functions.py infile

Aim:

Given an XSPEC model.dat file produce the output needed to be
included in sherpa/astro/xspec/src/_xspec.cc to fill up the
XSpecMethods array. The output is to stdout.

This code should be placed in _xspec.cc between the

  // Start model definitions

  // End model definitions

marks but please note that you will need to add back in the

    #ifdef XSPEC_x_y_x
    #else
    #endif

lines as there's no way to identify them (i.e. add them automatically)
from just a single model.dat file.

"""

from typing import Optional
import sys

from sherpa.astro.utils import xspec


def report_model_definitions(models: list[xspec.ModelDefinition]) -> None:
    """Print to stdout the code needed to bind to the models."""

    prev = None
    for model in models:
        prev = print_model(model, prev)


def print_model(model: xspec.ModelDefinition,
                prev: Optional[str]) -> Optional[str]:
    """Print out a model"""

    try:
        code = xspec.model_to_compiled(model)
    except ValueError:
        # Assume an unsupported model
        # sys.stderr.write(f"# Unsupported: {model.name} / {model.funcname}\n")
        return prev

    # Identify, with a comment, the per-spectrum models.
    #
    try:
        per_spectrum = model.flags[1] == 1
    except IndexError:
        per_spectrum = False

    # If we have changed "type" then print a new line
    if prev is not None and prev != model.modeltype:
        print("")

    # Add in the model name as a comment
    label = f"XS{model.name}"
    if per_spectrum:
        label += " [per spectrum]"

    print(f"{code[0]:50s} // {label}")
    return model.modeltype


if __name__ == "__main__":

    if len(sys.argv) != 2:
        sys.stderr.write(f"Usage: {sys.argv[0]} infile\n")
        sys.exit(1)

    # This errors out in case of significant model-support issues
    # (e.g. a periodic model parameter) but can also just be an
    # issue with a poorly-specified interchange format (the model.dat
    # file) which may just need changes to this routine.
    #
    mdefs = xspec.parse_xspec_model_description(sys.argv[1])
    report_model_definitions(mdefs)
