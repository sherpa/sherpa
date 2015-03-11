# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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

import numpy
from sherpa.models.parameter import Parameter, tinyval
from sherpa.models.model import ArithmeticModel, CompositeModel, modelCacher1d
from sherpa.astro.utils import apply_pileup
from sherpa.utils.err import ModelErr
from sherpa.utils import *
import sherpa.astro.models._modelfcts

__all__ = ('Atten', 'BBody', 'BBodyFreq', 'Beta1D', 'BPL1D', 'Dered', 'Edge',
           'LineBroad', 'Lorentz1D', 'NormBeta1D', 'Schechter',
           'Beta2D', 'DeVaucouleurs2D', 'HubbleReynolds', 'Lorentz2D',
           'JDPileup', 'MultiResponseSumModel', 'Sersic2D', 'Disk2D', 
           'Shell2D')


class Atten(ArithmeticModel):

    def __init__(self, name='atten'):
        self.hcol = Parameter(name, 'hcol', 1e+20, 1e+17, 1e+24, tinyval)
        self.heiRatio = Parameter(name, 'heiRatio', 0.1, 0, 1)
        self.heiiRatio = Parameter(name, 'heiiRatio', 0.01, 0, 1)
        ArithmeticModel.__init__(self, name,
                                 (self.hcol, self.heiRatio, self.heiiRatio))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        # atten should act like xsphabs, never integrate.
        kwargs['integrate']=False
        return _modelfcts.atten(*args, **kwargs)


class BBody(ArithmeticModel):

    def __init__(self, name='bbody'):
        self.space = Parameter(name, 'space', 0, 0, 1, 0, 1,
                               '0 - energy | 1 - wave', alwaysfrozen=True)
        self.kT = Parameter(name, 'kT', 1, 0.1, 1000, 0, 1e+10, 'keV')
        self.ampl = Parameter(name, 'ampl', 1, 1e-20, 1e+20, 1e-20, 1e+20)
        ArithmeticModel.__init__(self, name, (self.space, self.kT, self.ampl))


    def guess(self, dep, *args, **kwargs):
        if args[0][0] > args[0][-1]:
            self.space.val = 1
        Emax = get_peak(dep, *args)
        tMax = Emax/1.594
        kt = {'val':tMax,
              'min':tMax/_guess_ampl_scale,
              'max':tMax*_guess_ampl_scale}
        param_apply_limits(kt, self.kt, **kwargs)

        norm = guess_amplitude(dep, *args)
        modampl = norm['val']*(numpy.exp(Emax/tMax)-1)/numpy.square(Emax)
        mod = {'val': modampl,
               'min': modampl/_guess_ampl_scale,
               'max': modampl*_guess_ampl_scale}
        param_apply_limits(mod, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.bbody(*args, **kwargs)


class BBodyFreq(ArithmeticModel):

    def __init__(self, name='bbodyfreq'):
        self.T = Parameter(name, 'T', 1e+06, 1000, 1e+10, 1000, 1e+10)
        self.ampl = Parameter(name, 'ampl', 1, 0, hard_min=0)
        ArithmeticModel.__init__(self, name, (self.T, self.ampl))


    def guess(self, dep, *args, **kwargs):
        vmax = get_peak(dep, *args)
        tMax = vmax/5.88e+10
        t = {'val':tMax,
             'min':tMax/_guess_ampl_scale,
             'max':tMax*_guess_ampl_scale}
        norm = guess_amplitude(dep, *args)
        c_cm = 2.99792458e+10
        h_erg = 6.6260693e-27
        factor = numpy.exp(2.82)*numpy.square(c_cm)/h_erg/2.
        modampl = norm['val']*factor/numpy.power(vmax,3.)
        mod = {'val': modampl,
               'min': modampl/_guess_ampl_scale,
               'max': modampl*_guess_ampl_scale}
        param_apply_limits(mod, self.ampl, **kwargs)
        param_apply_limits(t, self.t, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.bbodyfreq(*args, **kwargs)


class Beta1D(ArithmeticModel):

    def __init__(self, name='beta1d'):
        self.r0 = Parameter(name, 'r0', 1, tinyval, hard_min=tinyval)
        self.beta = Parameter(name, 'beta', 1, 1e-05, 10, 1e-05, 10)
        self.xpos = Parameter(name, 'xpos', 0, 0, frozen=True)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.r0, self.beta, self.xpos, self.ampl))

    def get_center(self):
        return (self.xpos.val,)

    def set_center(self, xpos, *args, **kwargs):
        self.xpos.set(xpos)

    def guess(self, dep, *args, **kwargs):
        pos = get_position(dep, *args)
        param_apply_limits(pos, self.xpos, **kwargs)

        ref = guess_reference(self.r0.min, self.r0.max, *args)
        param_apply_limits(ref, self.r0, **kwargs)

        norm = guess_amplitude_at_ref(self.r0.val, dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)


    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.beta1d(*args, **kwargs)


class BPL1D(ArithmeticModel):

    def __init__(self, name='bpl1d'):
        self.gamma1 = Parameter(name, 'gamma1', 0, -10, 10)
        self.gamma2 = Parameter(name, 'gamma2', 0, -10, 10)
        self.eb = Parameter(name, 'eb', 100, 0, 1000, 0)
        self.ref = Parameter(name, 'ref', 1, frozen=True)
        self.ampl = Parameter(name, 'ampl', 1, 1e-20)
        ArithmeticModel.__init__(self, name,
                                 (self.gamma1, self.gamma2,
                                  self.eb, self.ref, self.ampl))


    def guess(self, dep, *args, **kwargs):
        pos = get_position(dep, *args)
        param_apply_limits(pos, self.eb, **kwargs)

        ref = guess_reference(self.ref.min, self.ref.max, *args)
        param_apply_limits(ref, self.ref, **kwargs)

        norm = guess_amplitude_at_ref(self.ref.val, dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.bpl1d(*args, **kwargs)


class Dered(ArithmeticModel):

    def __init__(self, name='dered'):
        self.rv = Parameter(name, 'rv', 10, 1e-10, 1000, tinyval)
        self.nhgal = Parameter(name, 'nhgal', 1e-07, 1e-07, 100000)
        ArithmeticModel.__init__(self, name, (self.rv, self.nhgal))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.dered(*args, **kwargs)


class Edge(ArithmeticModel):

    def __init__(self, name='edge'):
        self.space = Parameter(name, 'space', 0, 0, 1, 0, 1,
                               '0 - energy | 1 - wave')
        self.thresh = Parameter(name, 'thresh', 1, 0, hard_min=0)
        self.abs = Parameter(name, 'abs', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.space, self.thresh, self.abs))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.edge(*args, **kwargs)


class LineBroad(ArithmeticModel):

    def __init__(self, name='linebroad'):
        self.ampl = Parameter(name, 'ampl', 1, 0, hard_min=0)
        self.rest = Parameter(name, 'rest', 1000, tinyval, hard_min=tinyval)
        self.vsini = Parameter(name, 'vsini', tinyval, tinyval,
                               hard_min=tinyval)
        ArithmeticModel.__init__(self, name,
                                 (self.ampl, self.rest, self.vsini))


    def guess(self, dep, *args, **kwargs):
        ref = guess_reference(self.rest.min, self.rest.max, *args)
        param_apply_limits(ref, self.rest, **kwargs)

        norm = guess_amplitude_at_ref(self.rest.val, dep, *args)
        fwhm = get_fwhm(dep, *args)
        c_km = 2.99792458e+5
        if self.rest.val != 0:
            vsini = 2.*fwhm*c_km/numpy.sqrt(3.)/self.rest.val
            vs = {'val': vsini,
                  'min': vsini/_guess_ampl_scale,
                  'max': vsini*_guess_ampl_scale }
            param_apply_limits(vs, self.vsini, **kwargs)

        modampl = norm['val']*numpy.pi*self.vsini.val*self.rest.val/2./c_km
        mod = {'val': modampl,
               'min': modampl/_guess_ampl_scale,
               'max': modampl*_guess_ampl_scale}
        param_apply_limits(mod, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.linebroad(*args, **kwargs)


class Lorentz1D(ArithmeticModel):

    def __init__(self, name='lorentz1d'):
        self.fwhm = Parameter(name, 'fwhm', 10, 0, hard_min=0)
        self.pos = Parameter(name, 'pos', 1)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.fwhm, self.pos, self.ampl))

    def get_center(self):
        return (self.pos.val,)

    def set_center(self, pos, *args, **kwargs):
        self.pos.set(pos)

    def guess(self, dep, *args, **kwargs):
        pos = get_position(dep, *args)
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(pos, self.pos, **kwargs)
        param_apply_limits(fwhm, self.fwhm, **kwargs)

        norm = guess_amplitude(dep, *args)
        if fwhm != 10:
            aprime = norm['val']*self.fwhm.val*numpy.pi/2.
            ampl = {'val': aprime,
                    'min': aprime/_guess_ampl_scale,
                    'max': aprime*_guess_ampl_scale}
            param_apply_limits(ampl, self.ampl, **kwargs)
        else:
            param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.lorentz1d(*args, **kwargs)


class NormBeta1D(ArithmeticModel):

    def __init__(self, name='normbeta1d'):
        self.pos = Parameter(name, 'pos', 0)
        self.width = Parameter(name, 'width', 1, tinyval, hard_min=tinyval)
        self.index = Parameter(name, 'index', 2.5, 0.5, 1000, 0.5)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.pos, self.width, self.index, self.ampl))

    def get_center(self):
        return (self.pos.val,)

    def set_center(self, pos, *args, **kwargs):
        self.pos.set(pos)

    def guess(self, dep, *args, **kwargs):
        ampl = guess_amplitude(dep, *args)
        pos = get_position(dep, *args)
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(pos, self.pos, **kwargs)
        norm = (fwhm['val']*numpy.sqrt(numpy.pi)*
                numpy.exp(lgam(self.index.val-0.5)-lgam(self.index.val)))
        for key in ampl.keys():
            ampl[key] *= norm
        param_apply_limits(ampl, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.nbeta1d(*args, **kwargs)


class Schechter(ArithmeticModel):

    def __init__(self, name='schechter'):
        self.alpha = Parameter(name, 'alpha', 10)
        self.ref = Parameter(name, 'ref', 1)
        self.norm = Parameter(name, 'norm', 1)
        ArithmeticModel.__init__(self, name, (self.alpha, self.ref, self.norm))


    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.norm, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        if not self.integrate:
            raise ModelErr('alwaysint', self.name)
        return _modelfcts.schechter(*args, **kwargs)


class Beta2D(ArithmeticModel):

    def __init__(self, name='beta2d'):
        self.r0 = Parameter(name, 'r0', 10, tinyval, hard_min=tinyval)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999,
                               frozen=True)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians', True)
        self.ampl = Parameter(name, 'ampl', 1)
        self.alpha = Parameter(name, 'alpha', 1, -10, 10)
        ArithmeticModel.__init__(self, name,
                                 (self.r0, self.xpos, self.ypos, self.ellip,
                                  self.theta, self.ampl, self.alpha))
        self.cache = 0

    def get_center(self):
        return (self.xpos.val, self.ypos.val)

    def set_center(self, xpos, ypos, *args, **kwargs):
        self.xpos.set(xpos)
        self.ypos.set(ypos)

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        rad = guess_radius(*args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(rad, self.r0, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.beta2d(*args, **kwargs)


class DeVaucouleurs2D(ArithmeticModel):

    def __init__(self, name='devaucouleurs2d'):
        self.r0 = Parameter(name, 'r0', 10, 0, hard_min=0)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians')
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.r0, self.xpos, self.ypos, self.ellip,
                                  self.theta, self.ampl))
        self.cache = 0

    def get_center(self):
        return (self.xpos.val, self.ypos.val)

    def set_center(self, xpos, ypos, *args, **kwargs):
        self.xpos.set(xpos)
        self.ypos.set(ypos)

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        rad = guess_radius(*args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(rad, self.r0, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.devau(*args, **kwargs)


class HubbleReynolds(ArithmeticModel):

    def __init__(self, name='hubblereynolds'):
        self.r0 = Parameter(name, 'r0', 10, 0, hard_min=0)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians')
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.r0, self.xpos, self.ypos, self.ellip,
                                  self.theta, self.ampl))
        self.cache = 0

    def get_center(self):
        return (self.xpos.val, self.ypos.val)

    def set_center(self, xpos, ypos, *args, **kwargs):
        self.xpos.set(xpos)
        self.ypos.set(ypos)

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        rad = guess_radius(*args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(rad, self.r0, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.hr(*args, **kwargs)


class Lorentz2D(ArithmeticModel):

    def __init__(self, name='lorentz2d'):
        self.fwhm = Parameter(name, 'fwhm', 10, tinyval, hard_min=tinyval)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999,
                               frozen=True)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians',frozen=True)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.fwhm, self.xpos, self.ypos, self.ellip,
                                  self.theta, self.ampl))
        self.cache = 0

    def get_center(self):
        return (self.xpos.val, self.ypos.val)

    def set_center(self, xpos, ypos, *args, **kwargs):
        self.xpos.set(xpos)
        self.ypos.set(ypos)

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.lorentz2d(*args, **kwargs)


class JDPileup(ArithmeticModel):

    def __init__(self, name='jdpileup'):
        self.alpha = Parameter(name, 'alpha', 0.5, 0, 1, 0, 1)
        self.g0 = Parameter(name, 'g0', 1, tinyval, 1, tinyval, 1, frozen=True)
        self.f = Parameter(name, 'f', 0.95, 0.9, 1, 0, 1)
        self.n = Parameter(name, 'n', 1, tinyval, 100, tinyval, 2048,
                           alwaysfrozen=True)
        self.ftime = Parameter(name, 'ftime', 3.241, tinyval, 5, tinyval, 100,
                               'sec', alwaysfrozen=True)
        self.fracexp = Parameter(name, 'fracexp', 0.987, 0, 1, 0, 1,
                                 alwaysfrozen=True)
        self.nterms = Parameter(name, 'nterms', 30, 1, 100, 1, 100,
                                alwaysfrozen=True)
        self._results = None
        ArithmeticModel.__init__(self, name,
                                 (self.alpha, self.g0, self.f, self.n,
                                  self.ftime, self.fracexp, self.nterms))

    def __str__(self):
        s = ArithmeticModel.__str__(self)

        if self._results is not None:
            pileup_fractions, integral_ae, num_terms = self._results

            sum = 0.0
            if num_terms > 0:
                sum = pileup_fractions[1]

            sum_piled = pileup_fractions[2:num_terms+1].sum()
            sum += sum_piled

            pn = numpy.exp(-integral_ae)

            s += '\n\n'

            for i in xrange(1, num_terms+1):
                pn *= integral_ae / float(i)
                s += '   %d: %g  %g\n' % (i, pn, pileup_fractions[i]/sum)

            s += '   *** pileup fraction: %g' % (sum_piled/sum)

        return s

    def calc(self, p, arf_source, exposure_time, min_energy, max_energy,
             specresp, model):
        (alpha, g0, psf_frac, num_regions, frame_time, fracexpo,
         max_num_terms) = p
        
        out = apply_pileup(arf_source, exposure_time, int(max_num_terms),
                           min_energy, max_energy, specresp, fracexpo,
                           frame_time, alpha, g0, num_regions, psf_frac, model)
        self._results = out[1:]
        return out[0]


class MultiResponseSumModel(ArithmeticModel):
    pass

class Sersic2D(ArithmeticModel):

    def __init__(self, name='sersic2d'):
        self.r0 = Parameter(name, 'r0', 10, 0, hard_min=0)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians')
        self.ampl = Parameter(name, 'ampl', 1)
        self.n = Parameter(name,'n', 1, .1, 10, 0.01, 100, frozen=True )
        ArithmeticModel.__init__(self, name,
                                 (self.r0, self.xpos, self.ypos, self.ellip,
                                  self.theta, self.ampl, self.n))
        self.cache = 0

    def get_center(self):
        return (self.xpos.val, self.ypos.val)

    def set_center(self, xpos, ypos, *args, **kwargs):
        self.xpos.set(xpos)
        self.ypos.set(ypos)

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        rad = guess_radius(*args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(rad, self.r0, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.sersic(*args, **kwargs)

### disk2d and shell2d models 
### Contributed by Christoph Deil of the HESS project
### Added to CIAO 4.6 for Dec. 2013 release SMD

class Disk2D(ArithmeticModel):
    def __init__(self, name='disk2d'):
        self.xpos = Parameter(name, 'xpos', 0) # p[0]
        self.ypos = Parameter(name, 'ypos', 0) # p[1]
        self.ampl = Parameter(name, 'ampl', 1) # p[2]
        self.r0 = Parameter(name, 'r0', 1, 0) # p[3]
        ArithmeticModel.__init__(self, name, (self.xpos, self.ypos, self.ampl, self.r0))
 
    def calc(self, p, x, y, *args, **kwargs):
        # Compute radii
        r2 = (x - p[0]) ** 2 + (y - p[1]) ** 2
 
        # Return ampl when r2 <= r0 else return 0
        return numpy.select([r2 <= p[3] ** 2], [p[2]])

class Shell2D(ArithmeticModel):
    def __init__(self, name='shell2d'):
        self.xpos = Parameter(name, 'xpos', 0) # p[0]
        self.ypos = Parameter(name, 'ypos', 0) # p[1]
        self.ampl = Parameter(name, 'ampl', 1) # p[2]
        self.r0 = Parameter(name, 'r0', 1, 0) # p[3]
        self.width = Parameter(name, 'width', 0.1, 0)
        ArithmeticModel.__init__(self, name, (self.xpos, self.ypos, self.ampl, self.r0, self.width))
 
    def calc(self, p, x, y, *args, **kwargs):
        """Homogeneously emitting spherical shell,
        projected along the z-direction
        (this is not 100% correct for very large shells on the sky)."""
        (xpos, ypos, ampl, r_0, width) = p
        
        r2 = (x - xpos) * (x - xpos) + (y - ypos) * (y - ypos)
        r_out = r_0 + width
        r_in2, r_out2 = r_0 * r_0, r_out * r_out
        # r_in3, r_out3 = r_in * r_in2, r_out * r_out2
        # We only call abs() in sqrt() to avoid warning messages.
        sphere_out = numpy.sqrt(numpy.abs(r_out2 - r2))
        sphere_in = numpy.sqrt(numpy.abs(r_in2 - r2))
        # Note: for r > r_out 'numpy.select' fills automatically zeros!
        non_normalized = numpy.select([r2 <= r_in2, r2 <= r_out2],
                                      [sphere_out - sphere_in, sphere_out])
        # Computed with Mathematica:
        integral = 2 * numpy.pi / 3 * (r_out ** 3 - r_0 ** 3)
        # integral = 1
        return ampl * integral * non_normalized
