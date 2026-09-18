#
#  Copyright (C) 2022
#  MIT
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
from sherpa.models import model

class PolarizationDependence(model.RegriddableModel1D):
    """A one-dimensional trapezoid.

    The model parameters are:

    stokes
        0, 1, 2 for Stokes I, U, Q
    """
    def __init__(self, name, *args):
        self.stokes = model.Parameter(name, 'stokes', 0, alwaysfrozen=True, hidden=True)
        args = args[0] + (self.stokes, )
        super().__init__(name, args)

    def calc(self, frac, ang):
        # check here that frac is in ange 0..1
        if self.stokes.val == 0:
            # Stokes I
            return 1
        elif self.stokes.val == 1:
            # Stokes Q
            return frac * np.sin(np.deg2rad(2 * ang))
        elif self.stokes.val == 2:
            # Stokes U
            return frac * np.cos(np.deg2rad(2 * ang))
        else:
            raise ValueError('Set stokes to one of 0, 1, 2 for Stokes I, Q, U')


class PolConst(PolarizationDependence):
    '''A model with a constant polarization angle and fraction'''
    def __init__(self, name='polconst'):        
        self.polangle = model.Parameter(name, 'polangle', 10.)
        self.polfrac = model.Parameter(name, 'polfrac', .3, min=0, hard_min=0, max=1, hard_max=1)

        super().__init__(name,
                        (self.polangle, self.polfrac))

    def calc(self, pars, x, *args, **kwargs):
        """Evaluate the model"""
        (ang, frac, stokes) = pars
        return super().calc(frac, ang)


class PolLinearDep(PolarizationDependence):
    """Polarization where angle and fraction depend linearly on energy.
    
    The model parameters are:

    stokes
        0, 1, 2 for Stokes I, U, Q
    polangle
        Polarization angle extrapolated to 0 keV.
        (Note that I did nothing fancy for angles wrapping about 360 deg.)
    polangslope
        Slope of polarization angle with energy
    polfrac
        Polarization fraction extrapolated to 0 keV
    pofracslope
        Slope of the polarization fraction with energy
    """

    def __init__(self, name='pollineardep'):        
        self.polangle = model.Parameter(name, 'polangle', 0)
        self.polangslope = model.Parameter(name, 'polangslope', 10.)
        # demostrate some min/max settings. Since I chose polfrac to be extrapolated to 0 keV
        # the value could be less than 0 or more than 1, as long at it's between
        # 0 and 1 in the range 2..8 keV. 
        # So, for demonstration pupose only.
        self.polfrac = model.Parameter(name, 'polfrac', .3, min=0, hard_min=0)
        self.polfracslope = model.Parameter(name, 'polfracslope', .01)

        super().__init__(name,
                        (self.polangle, self.polangslope,
                         self.polfrac, self.polfracslope))

    def calc(self, pars, x, *args, **kwargs):
        """Evaluate the model"""

        # If given an integrated data set, use the center of the bin
        # Since we use linear relations, that's correct
        if len(args) == 1:
            x = (x + args[0]) / 2
        
        (ang, angslope, frac, fracslope, stokes) = pars
        
        return super().calc(frac + x * fracslope, ang + x * angslope)