#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2010)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
import numpy
from parameter import Parameter, tinyval
from model import ArithmeticModel, modelCacher1d, CompositeModel, \
    ArithmeticFunctionModel
from sherpa.utils.err import ModelErr
from sherpa.utils import *
import _modelfcts

import logging
warning = logging.getLogger(__name__).warning


__all__ = ('Box1D', 'Const1D', 'Cos', 'Delta1D', 'Erf', 'Erfc', 'Exp', 'Exp10',
           'Gauss1D', 'Log', 'Log10', 'LogParabola', 'NormGauss1D', 'Poisson',
           'Polynom1D', 'PowLaw1D', 'Scale1D', 'Sin', 'Sqrt', 'StepHi1D',
           'StepLo1D', 'Tan', 'Box2D', 'Const2D', 'Delta2D', 'Gauss2D', 
           'SigmaGauss2D',
           'NormGauss2D', 'Polynom2D', 'Scale2D', 'UserModel', 'TableModel',
           'Integrate1D')

DBL_EPSILON = numpy.finfo(numpy.float).eps


class Box1D(ArithmeticModel):

    def __init__(self, name='box1d'):
        self.xlow = Parameter(name, 'xlow', 0)
        self.xhi = Parameter(name, 'xhi', 0)
        self.ampl = Parameter(name, 'ampl', 1, -1, 1)
        ArithmeticModel.__init__(self, name, (self.xlow, self.xhi, self.ampl))

    def guess(self, dep, *args, **kwargs):
        lo, hi = guess_bounds(args[0])
        norm = guess_amplitude(dep, *args)
        param_apply_limits(lo, self.xlow, **kwargs)
        param_apply_limits(hi, self.xhi, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.box1d(*args, **kwargs)


class Const1D(ArithmeticModel):

    def __init__(self, name='const1d'):
        self.c0 = Parameter(name, 'c0', 1)
        ArithmeticModel.__init__(self, name, (self.c0,))

    def guess(self, dep, *args, **kwargs):
        min = dep.min()
        max = dep.max()
        ylo = 0
        if numpy.abs(min - 0) > DBL_EPSILON:
            ylo = min/100.
            if min < 0: ylo = 0
        yhi = 0
        if numpy.abs(max - 0) > DBL_EPSILON:
            yhi = -10*max
            if max > 0: yhi = 100*max
        param_apply_limits({ 'val':(max+min)/2., 'min':ylo, 'max':yhi },
                           self.c0, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.const1d(*args, **kwargs)


class Cos(ArithmeticModel):

    def __init__(self, name='cos'):
        self.period = Parameter(name, 'period', 1, 1e-10, 10, tinyval)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.ampl = Parameter(name, 'ampl', 1, 1e-05, hard_min=0)
        ArithmeticModel.__init__(self, name,
                                 (self.period, self.offset, self.ampl))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.cos(*args, **kwargs)


class Delta1D(ArithmeticModel):

    def __init__(self, name='delta1d'):
        self.pos = Parameter(name, 'pos', 0)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.pos, self.ampl))

    def get_center(self):
        return (self.pos.val,)

    def set_center(self, pos, *args, **kwargs):
        self.pos.set(pos)

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        pos = get_position(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(pos, self.pos, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.delta1d(*args, **kwargs)


class Erf(ArithmeticModel):

    def __init__(self, name='erf'):
        self.ampl = Parameter(name, 'ampl', 1, 0)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.sigma = Parameter(name, 'sigma', 1, 1e-10, 10, tinyval)
        ArithmeticModel.__init__(self, name,
                                 (self.ampl, self.offset, self.sigma))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.erf(*args, **kwargs)


class Erfc(ArithmeticModel):

    def __init__(self, name='erfc'):
        self.ampl = Parameter(name, 'ampl', 1, 0)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.sigma = Parameter(name, 'sigma', 1, 1e-10, 10, tinyval)
        ArithmeticModel.__init__(self, name,
                                 (self.ampl, self.offset, self.sigma))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.erfc(*args, **kwargs)


class Exp(ArithmeticModel):

    def __init__(self, name='exp'):
        self.offset = Parameter(name, 'offset', 0)
        self.coeff = Parameter(name, 'coeff', -1)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.offset, self.coeff, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.exp(*args, **kwargs)


class Exp10(ArithmeticModel):

    def __init__(self, name='exp10'):
        self.offset = Parameter(name, 'offset', 0)
        self.coeff = Parameter(name, 'coeff', -1)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.offset, self.coeff, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.exp10(*args, **kwargs)


class Gauss1D(ArithmeticModel):

    def __init__(self, name='gauss1d'):
        self.fwhm = Parameter(name, 'fwhm', 10, tinyval, hard_min=tinyval)
        self.pos = Parameter(name, 'pos', 0)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos, self.ampl))

    def get_center(self):
        return (self.pos.val,)

    def set_center(self, pos, *args, **kwargs):
        self.pos.set(pos)

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        pos = get_position(dep, *args)
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(pos, self.pos, **kwargs)
        param_apply_limits(fwhm, self.fwhm, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.gauss1d(*args, **kwargs)


class Log(ArithmeticModel):

    def __init__(self, name='log'):
        self.offset = Parameter(name, 'offset', 0)
        self.coeff = Parameter(name, 'coeff', -1)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.offset, self.coeff, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.log(*args, **kwargs)


class Log10(ArithmeticModel):

    def __init__(self, name='log10'):
        self.offset = Parameter(name, 'offset', 0)
        self.coeff = Parameter(name, 'coeff', -1)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.offset, self.coeff, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.log10(*args, **kwargs)


class LogParabola(ArithmeticModel):

    def __init__(self, name='logparabola'):
        self.ref = Parameter(name, 'ref', 1, alwaysfrozen=True)
        self.c1 = Parameter(name, 'c1', 1)
        self.c2 = Parameter(name, 'c2', 1)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name, (self.ref,self.c1,
                                              self.c2,self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.logparabola(*args, **kwargs)


_gfactor = numpy.sqrt(numpy.pi/(4*numpy.log(2)))

class NormGauss1D(ArithmeticModel):

    def __init__(self, name='normgauss1d'):
        self.fwhm = Parameter(name, 'fwhm', 10, tinyval, hard_min=tinyval)
        self.pos = Parameter(name, 'pos', 0)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos, self.ampl))

    def get_center(self):
        return (self.pos.val,)

    def set_center(self, pos, *args, **kwargs):
        self.pos.set(pos)

    def guess(self, dep, *args, **kwargs):
        ampl = guess_amplitude(dep, *args)
        pos = get_position(dep, *args)
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(pos, self.pos, **kwargs)
        param_apply_limits(fwhm, self.fwhm, **kwargs)

        # Apply normalization factor to guessed amplitude
        norm = numpy.sqrt(numpy.pi/_gfactor)*self.fwhm.val
        for key in ampl.keys():
            if ampl[key] is not None:
                ampl[key] *= norm
        param_apply_limits(ampl, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.ngauss1d(*args, **kwargs)


class Poisson(ArithmeticModel):

    def __init__(self, name='poisson'):
        self.mean = Parameter(name, 'mean', 1, 1e-05, hard_min=tinyval)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.mean, self.ampl))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        pos = get_position(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(pos, self.mean, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.poisson(*args, **kwargs)


class Polynom1D(ArithmeticModel):

    def __init__(self, name='polynom1d'):
        pars = []
        
        for i in xrange(9):
            pars.append(Parameter(name, 'c%d' % i, 0, frozen=True))
        pars[0].val = 1
        pars[0].frozen = False
        for p in pars:
            setattr(self, p.name, p)

        self.offset = Parameter(name, 'offset', 0, frozen=True)
        pars.append(self.offset)

        ArithmeticModel.__init__(self, name, pars)


    def guess(self, dep, *args, **kwargs):
        xmin = args[0].min()
        xmax = args[0].max()
        ymin = dep.min()
        ymax = dep.max()
        dydx = (ymax-ymin)/(xmax-xmin)
        dydx2 = (ymax-ymin)/((xmax-xmin)*(xmax-xmin))

        xlo = 0
        if numpy.abs(xmin - 0) >= DBL_EPSILON:
            xlo = -xmin
            if xmin < 0:
                xlo = xmin

        xhi = 0
        if numpy.abs(xmax - 0) >= DBL_EPSILON:
            xhi = -xmax
            if xmax > 0:
                xhi = xmax

        ylo = 0
        if numpy.abs(ymin - 0) >= DBL_EPSILON:
            ylo = -ymin
            if ymin < 0:
                ylo = ymin

        yhi = 0
        if numpy.abs(ymax - 0) >= DBL_EPSILON:
            yhi = -ymax
            if ymax > 0:
                yhi = ymax

        c0 = {'val': (ymax+ymin)/2.0, 'min': ylo, 'max': yhi}
        c1 = {'val': 0.0, 'min': -100*dydx, 'max': 100*dydx}
        c2 = {'val': 0.0, 'min': -100*dydx2, 'max': 100*dydx2}
        c3 = {'val': 0.0, 'min': ylo, 'max': yhi}
        off = {'val': 0.0, 'min': xlo, 'max': xhi}

        param_apply_limits(c0, self.c0, **kwargs)
        param_apply_limits(c1, self.c1, **kwargs)
        param_apply_limits(c2, self.c2, **kwargs)
        param_apply_limits(c3, self.c3, **kwargs)
        param_apply_limits(c3, self.c4, **kwargs)
        param_apply_limits(c3, self.c5, **kwargs)
        param_apply_limits(c3, self.c6, **kwargs)
        param_apply_limits(c3, self.c7, **kwargs)
        param_apply_limits(c3, self.c8, **kwargs)
        param_apply_limits(off, self.offset, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.poly1d(*args, **kwargs)


class PowLaw1D(ArithmeticModel):

    def __init__(self, name='powlaw1d'):
        self.gamma = Parameter(name, 'gamma', 1, -10, 10)
        self.ref = Parameter(name, 'ref', 1, alwaysfrozen=True)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name, (self.gamma, self.ref, self.ampl))


    def guess(self, dep, *args, **kwargs):
        ref = guess_reference(self.ref.min, self.ref.max, *args)
        param_apply_limits(ref, self.ref, **kwargs)
        norm = guess_amplitude_at_ref(self.ref.val, dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

        
    @modelCacher1d
    def calc(self, pars, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        if kwargs['integrate']:
            # avoid numerical issues with C pow() function close to zero,
            # 0.0 +- ~1.e-14.  PowLaw1D integrated has multiple calls to 
            # pow(X, 1.0 - gamma).  So gamma values close to 1.0 +- 1.e-10
            # should be be 1.0 to avoid errors propagating in the calculated
            # model.
            if sao_fcmp(pars[0], 1.0, 1.e-10) == 0:
                # pars { gamma, ref, ampl }
                pars[0] = 1.0

        return _modelfcts.powlaw(pars, *args, **kwargs)


class Scale1D(Const1D):

    def __init__(self, name='scale1d'):
        Const1D.__init__(self, name)
        self.integrate=False


class Sin(ArithmeticModel):

    def __init__(self, name='sin'):
        self.period = Parameter(name, 'period', 1, 1e-10, 10, tinyval)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.ampl = Parameter(name, 'ampl', 1, 1e-05, hard_min=0)
        ArithmeticModel.__init__(self, name,
                                 (self.period, self.offset, self.ampl))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.sin(*args, **kwargs)


class Sqrt(ArithmeticModel):

    def __init__(self, name='sqrt'):
        self.offset = Parameter(name, 'offset', 0)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name, (self.offset, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.sqrt(*args, **kwargs)


class StepHi1D(ArithmeticModel):

    def __init__(self, name='stephi1d'):
        self.xcut = Parameter(name, 'xcut', 0)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name, (self.xcut, self.ampl))


    def guess(self, dep, *args, **kwargs):
        cut = guess_bounds(args[0], False)
        norm = guess_amplitude(dep, *args)
        param_apply_limits(cut, self.xcut, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.stephi1d(*args, **kwargs)


class StepLo1D(ArithmeticModel):

    def __init__(self, name='steplo1d'):
        self.xcut = Parameter(name, 'xcut', 0)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name, (self.xcut, self.ampl))


    def guess(self, dep, *args, **kwargs):
        cut = guess_bounds(args[0], False)
        norm = guess_amplitude(dep, *args)
        param_apply_limits(cut, self.xcut, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.steplo1d(*args, **kwargs)


class Tan(ArithmeticModel):

    def __init__(self, name='tan'):
        self.period = Parameter(name, 'period', 1, 1e-10, 10, tinyval)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.ampl = Parameter(name, 'ampl', 1, 1e-05, hard_min=0)
        ArithmeticModel.__init__(self, name,
                                 (self.period, self.offset, self.ampl))


    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.tan(*args, **kwargs)


class Box2D(ArithmeticModel):

    def __init__(self, name='box2d'):
        self.xlow = Parameter(name, 'xlow', 0)
        self.xhi = Parameter(name, 'xhi', 0)
        self.ylow = Parameter(name, 'ylow', 0)
        self.yhi = Parameter(name, 'yhi', 0)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.xlow, self.xhi, self.ylow, self.yhi,
                                  self.ampl))
        self.cache = 0


    def guess(self, dep, *args, **kwargs):
        xlo, xhi = guess_bounds(args[0])
        ylo, yhi = guess_bounds(args[1])
        norm = guess_amplitude2d(dep, *args)
        param_apply_limits(xlo, self.xlow, **kwargs)
        param_apply_limits(xhi, self.xhi, **kwargs)
        param_apply_limits(ylo, self.ylow, **kwargs)
        param_apply_limits(yhi, self.yhi, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.box2d(*args, **kwargs)


class Const2D(Const1D):

    def __init__(self, name='const2d'):
        Const1D.__init__(self, name)
        self.cache = 0

    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.const2d(*args, **kwargs)


class Scale2D(Const2D):

    def __init__(self, name='scale2d'):
        Const2D.__init__(self, name)
        self.integrate=False
        self.cache = 0


class Delta2D(ArithmeticModel):

    def __init__(self, name='delta2d'):
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.xpos, self.ypos, self.ampl))
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
        return _modelfcts.delta2d(*args, **kwargs)

class Gauss2D(ArithmeticModel):

    def __init__(self, name='gauss2d'):
        self.fwhm = Parameter(name, 'fwhm', 10, tinyval, hard_min=tinyval)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999,
                               frozen=True)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians', frozen=True)
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
        return _modelfcts.gauss2d(*args, **kwargs)

class SigmaGauss2D(Gauss2D):

    def __init__(self, name='sigmagauss2d'):
        self.sigma_a = Parameter(name, 'sigma_a', 10, tinyval, hard_min=tinyval)
        self.sigma_b = Parameter(name, 'sigma_b', 10, tinyval, hard_min=tinyval)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.theta = \
            Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, \
                          -2*numpy.pi, 4*numpy.pi, 'radians', frozen=True)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.sigma_a, self.sigma_b, self.xpos,
                                  self.ypos, self.theta, self.ampl))
        self.cache = 0

    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.sigmagauss2d(*args, **kwargs)


class NormGauss2D(ArithmeticModel):

    def __init__(self, name='normgauss2d'):
        self.fwhm = Parameter(name, 'fwhm', 10, tinyval, hard_min=tinyval)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999,
                               frozen=True)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians', frozen=True)
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
        ampl = guess_amplitude2d(dep, *args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)

        # Apply normalization factor to guessed amplitude
        norm = (numpy.pi/_gfactor)*self.fwhm.val*self.fwhm.val*numpy.sqrt(1.0 - (self.ellip.val*self.ellip.val))
        for key in ampl.keys():
            if ampl[key] is not None:
                ampl[key] *= norm
        param_apply_limits(ampl, self.ampl, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.ngauss2d(*args, **kwargs)

    
class Polynom2D(ArithmeticModel):

    def __init__(self, name='polynom2d'):
        self.c = Parameter(name, 'c', 1)
        self.cx1 = Parameter(name, 'cx1', 0)
        self.cx2 = Parameter(name, 'cx2', 0)
        self.cy1 = Parameter(name, 'cy1', 0)
        self.cy2 = Parameter(name, 'cy2', 0)
        self.cx1y1 = Parameter(name, 'cx1y1', 0)
        self.cx1y2 = Parameter(name, 'cx1y2', 0)
        self.cx2y1 = Parameter(name, 'cx2y1', 0)
        self.cx2y2 = Parameter(name, 'cx2y2', 0)
        ArithmeticModel.__init__(self, name,
                                 (self.c, self.cy1, self.cy2, self.cx1,
                                  self.cx1y1, self.cx1y2, self.cx2,
                                  self.cx2y1, self.cx2y2))
        self.cache = 0


    def guess(self, dep, *args, **kwargs):
        x0min = args[0].min()
        x0max = args[0].max()
        x1min = args[1].min()
        x1max = args[1].max()
        ymin = dep.min()
        ymax = dep.max()

        ylo = 0
        if numpy.abs(ymin - 0) > DBL_EPSILON:
            ylo = -ymin
            if ymin < 0: ylo = ymin

        yhi = 0
        if numpy.abs(ymax - 0) > DBL_EPSILON:
            yhi = -ymax
            if ymax > 0: yhi = ymax

        dydx0 = (ymax-ymin)/(x0max-x0min)
        dydx1 = (ymax-ymin)/(x1max-x1min)
        dyd2x0 = (ymax-ymin)/((x0max-x0min)*(x0max-x0min))
        dyd2x1 = (ymax-ymin)/((x1max-x1min)*(x1max-x1min))
        dydx0dx1 = (ymax-ymin)/((x0max-x0min)*(x1max-x1min))
        
        c     = {'val':(ymax+ymin)/2., 'min': ylo, 'max': yhi }
        cx1   = {'val': 0., 'min': -100*dydx0, 'max': 100*dydx0 }
        cy1   = {'val': 0., 'min': -100*dydx1, 'max': 100*dydx1 }
        cx2   = {'val': 0., 'min': -100*dyd2x0, 'max': 100*dyd2x0 }
        cy2   = {'val': 0., 'min': -100*dyd2x1, 'max': 100*dyd2x1 }
        cx1y1 = {'val': 0., 'min': -100*dydx0dx1, 'max': 100*dydx0dx1 }
        c22   = {'val': 0., 'min': ylo, 'max': yhi }
        
        param_apply_limits(c, self.c, **kwargs)
        param_apply_limits(cx1, self.cx1, **kwargs)
        param_apply_limits(cy1, self.cy1, **kwargs)
        param_apply_limits(cx2, self.cx2, **kwargs)
        param_apply_limits(cy2, self.cy2, **kwargs)
        param_apply_limits(cx1y1, self.cx1y1, **kwargs)
        param_apply_limits(c22, self.cx1y2, **kwargs)
        param_apply_limits(c22, self.cx2y1, **kwargs)
        param_apply_limits(c22, self.cx2y2, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.poly2d(*args, **kwargs)


class TableModel(ArithmeticModel):
    def __init__(self, name='tablemodel'):
        # these attributes should remain somewhat private
        # as not to conflict with user defined parameter names
        self.__x = None
        self.__y = None
        self.__filtered_y = None
        self.filename = None
        self.method = linear_interp  # interpolation method
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.ampl,))


    def __setstate__(self, state):
        self.__x = None
        self.__y = state.pop('_y', None)
        self.__filtered_y = state.pop('_filtered_y', None)
        self.filename = state.pop('_file', None)
        ArithmeticModel.__setstate__(self, state)

    def load(self, x, y):
        self.__y = y
        self.__x = x
 
        # Input grid is sorted!
        if x is not None:
            idx = numpy.asarray(x).argsort()
            self.__y = numpy.asarray(y)[idx]
            self.__x = numpy.asarray(x)[idx]


    def get_x(self):
        return self.__x


    def get_y(self):
        return self.__y


    def fold(self, data):
        mask = data.mask
        if self.__x is None and numpy.iterable(mask):
            if len(mask) != len(self.__y):
                raise ModelErr("filtermismatch", 'table model',
                               'data, (%s vs %s)' %
                               (len(self.__y), len(mask)))
            self.__filtered_y = self.__y[mask]

    @modelCacher1d
    def calc(self, p, x0, x1=None, *args, **kwargs):

        if self.__x is not None and self.__y is not None:
            return p[0] * interpolate(x0, self.__x, self.__y, function=self.method)

        elif (self.__filtered_y is not None and
              len(x0) == len(self.__filtered_y)):
            return p[0] * self.__filtered_y

        elif (self.__y is not None and
              len(x0) == len(self.__y)):
            return p[0] * self.__y

        raise ModelErr("filtermismatch", 'table model', 'data, (%s vs %s)' %
                       (len(self.__y), len(x0)))

class UserModel(ArithmeticModel):
    def __init__(self, name='usermodel', pars=None):
        # these attributes should remain somewhat private
        # as not to conflict with user defined parameter names
        self._y = []
        self._file = None
        if pars is None:
            self.ampl = Parameter(name, 'ampl', 1)
            pars = (self.ampl,)
        else:
            for par in pars:
                self.__dict__[par.name] = par
        ArithmeticModel.__init__(self, name, pars)



class Integrator1D(CompositeModel, ArithmeticModel):

    @staticmethod
    def wrapobj(obj):
	if isinstance(obj, ArithmeticModel):
	    return obj
	return ArithmeticFunctionModel(obj)

    def __init__(self, model, *otherargs, **otherkwargs):
        self.model = self.wrapobj(model)
        self.otherargs = otherargs
        self.otherkwargs = otherkwargs
        self._errflag = 0
	CompositeModel.__init__(self,
                                ('integrate1d(%s)' % self.model.name),
                                (self.model,))


    def startup(self):
        self.model.startup()
        self._errflag = 1
        CompositeModel.startup(self)


    def teardown(self):
        self.model.teardown()
        CompositeModel.teardown(self)


    def calc(self, p, xlo, xhi=None, **kwargs):
        if xhi is None:
            raise ModelErr('needsint')

	return _modelfcts.integrate1d(self.model.calc,
                                      p, xlo, xhi, **self.otherkwargs)

class Integrate1D(ArithmeticModel):

    def __init__(self, name='integrate1d'):
        tol = numpy.finfo(float).eps
        self.epsabs = Parameter(name, 'epsabs', tol, alwaysfrozen=True)
        self.epsrel = Parameter(name, 'epsrel', 0, alwaysfrozen=True)
        self.maxeval = Parameter(name, 'maxeval', 10000, alwaysfrozen=True)
        ArithmeticModel.__init__(self, name, (self.epsabs,self.epsrel,
                                              self.maxeval))


    def __call__(self, model):
        return Integrator1D(model,
                            epsabs=self.epsabs.val,
                            epsrel=self.epsrel.val,
                            maxeval=int(self.maxeval.val),
                            logger=warning)
