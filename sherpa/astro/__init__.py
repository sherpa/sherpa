#
#  Copyright (C) 2007, 2021  Smithsonian Astrophysical Observatory
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
'''Physical constants for astronomical spectra

This module contains some physical constants needed to perform conversions
that are common in the `sherpa.astro` module, e.g. from energy to wavelength.
'''

hc = 12.39841874  # nist.gov in [keV-Angstrom]
"""The product of Plank constant (h) and speed of light (c) using values from nist.gov in units appropriate for converting between keV and Angstroms."""

charge_e = 1.60217653e-09
'''elementary charge from nist.gov in units appropriate for converting between keV and erg.'''
