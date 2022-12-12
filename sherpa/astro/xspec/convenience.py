#
#  Copyright (C) 2022
#  MIT
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

try:
    from . import _xspec
except ImportError as ie:
    # It would be nicer for the user to say "from None" rather than
    # "from ie" here but it may lose important information (in the
    # case there is some problem initializing the XSPEC library,
    # rather than it just not being built).
    #
    raise ImportError("XSPEC support is not enabled") from ie

import numpy as np
from sherpa.models import model

__all__ = ['MultiTvapec', 'VAPECgrid']


class MultiTvapec(model.RegriddableModel1D):
    def __init__(self, name='multiTvapec', n_components=1):
        '''
        Alternative implementation idea: 
            Make many VAPEC models, link parameters and keep a list of them in here.
            That way, each temperature component would have a dedicated VAPEC instance. 
            If there is some internal cashing, that might speed up computations - but not that much,
            since the abundances are already matched.
        '''
        # Why do I need this line? Should be done at the end in ArithmeticModel.__init__?
        self.name = name
        
        self.n_components = n_components
        self.vapec = xspec.XSvapec(name=f'xs')
        pars = []
        for i in range(n_components):

            p = model.Parameter(name, f'kT_{i + 1}',
                                 # +1 because can't set to 0 for i=0
                                val=min(i + 1, self.vapec.kT.hard_max),
                                min=self.vapec.kT.min,
                                hard_min=self.vapec.kT.hard_min,
                                max=self.vapec.kT.max,
                                hard_max=self.vapec.kT.hard_max)
            setattr(self, f'kT_{i + 1}', p)
            pars.append(p)
        for i in range(n_components):
            p = model.Parameter(name, f'norm_{i + 1}',
                                val=self.vapec.norm.val,
                                min=self.vapec.norm.min,
                                hard_min=self.vapec.norm.hard_min,
                                max=self.vapec.norm.max,
                                hard_max=self.vapec.norm.hard_max)
            setattr(self, f'norm_{i + 1}', p)
            pars.append(p)            
            
        for j in range(1, len(self.vapec.pars) - 1):
            par = self.vapec.pars[j]
            p = model.Parameter(name, par.name, val=par.val,
                    min=par.min, hard_min=par.hard_min,
                    max=par.max, hard_max=par.hard_max,
                               frozen=True)
            setattr(self, par.name, p)
            pars.append(p)
        model.ArithmeticModel.__init__(self, name, tuple(pars))
        
    def calc(self, pars, *args, **kwargs):
        components = []
        # Can I replace the following loop with sherpa.utils.parallel_map?
        # Or do I have to use different vapec objects to avoid some race condition?
        # Is this even the time critical step? Should check that...
        for i in range(self.n_components):
            components.append(self.vapec.calc([pars[i]] + \
                                              pars[2 * self.n_components:] + \
                                              [pars[self.n_components + i]],
                                             *args, **kwargs))
        return np.sum(components, axis=0)


class VAPECgrid(MultiTvapec):
    '''VAPEC grid at fixed temperatures'''
    def __init__(self, name='multiTvapec', temp_grid=[.5, 2., 5.]):
        self.temp_grid = temp_grid
        super().__init__(name=name, n_components=len(temp_grid))
        for i, t in enumerate(temp_grid):
            temp_par = getattr(self, f'kT_{i + 1}')
            temp_par.val = t
            temp_par.frozen = True