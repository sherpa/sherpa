#!/usr/bin/env python

#
#  Copyright (C) 2025
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

  ./update_xspec_docs.py infile

Aim:

Display updated contents for parts of the

    docs/model_classes/astro_xspec.rst

file. The input file is the model.dat file.

"""

import sys

from sherpa.astro.utils import xspec


def restrict_models(models: list[xspec.ModelDefinition]
                    ) -> list[xspec.ModelDefinition]:
    """What models do we support?"""

    out = []
    for model in models:
        if not isinstance(model, (xspec.AddModelDefinition,
                                  xspec.MulModelDefinition,
                                  xspec.ConModelDefinition)):
            continue

        try:
            if model.flags[1] == 1:
                # drop per-spectrum models
                continue

        except:
            pass

        out.append(model)

    assert len(out) > 0, "What happened to the models?"
    return out


def report_classes(models: list[xspec.ModelDefinition]) -> None:
    """List the classes"""

    def p1(x):
        print(f"   {x}")

    def p2(x):
        print(f"      {x}")


    p1(".. rubric:: Classes")
    print()
    p1(".. autosummary::")
    p2(":toctree: api")
    print()

    names = [m.name for m in models]

    for name in sorted(names):
        p2(f"XS{name}")

    print("")


def report_inheritance(models: list[xspec.ModelDefinition],
                       cls: xspec.ModelDefinition
                       ) -> None:
    """List the inheritance-diagram line for these models."""

    wanted = [m.name for m in models if isinstance(m, cls)]
    names = sorted(wanted)
    diagram = " ".join([f"XS{n}" for n in names])

    print("")
    print(f".. inheritance-diagram:: {diagram}")
    print("   :parts: 1")


def report_docs(models: list[xspec.ModelDefinition]) -> None:
    """Given the list of models, what should the docs say?"""

    wanted = restrict_models(models)
    report_classes(wanted)

    print("Class Inheritance Diagram")
    print("=========================")

    report_inheritance(wanted, xspec.AddModelDefinition)
    report_inheritance(wanted, xspec.MulModelDefinition)
    report_inheritance(wanted, xspec.ConModelDefinition)


if __name__ == "__main__":

    if len(sys.argv) != 2:
        sys.stderr.write(f"Usage: {sys.argv[0]} infile\n")
        sys.exit(1)

    # The script could instead process sherpa/astro/xspec/__init__.py
    # but it is better to be explicit here.
    #
    mdefs = xspec.parse_xspec_model_description(sys.argv[1])
    report_docs(mdefs)
