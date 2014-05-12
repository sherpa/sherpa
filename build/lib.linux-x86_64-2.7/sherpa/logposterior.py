from sherpa.stats import Likelihood, truncation_value, Stat
from sherpa.models import Parameter
from sherpa.utils import NoNewAttributesAfterInit
from sherpa.utils.err import StatErr
from itertools import izip
import numpy


class Prior(Likelihood):

    # Provide a Model-like parameter interface
    def __getattr__(self, name):

        par = self.__dict__.get(name.lower())
        if (par is not None) and isinstance(par, Parameter):
            return par
        raise AttributeError("'%s' object has no attribute '%s'" %
                             (type(self).__name__, name))

    def __setattr__(self, name, val):

        par = getattr(self, name.lower(), None)
        if (par is not None) and isinstance(par, Parameter):
            par.val = val
        else:
            NoNewAttributesAfterInit.__setattr__(self, name, val)


    def __init__(self, statfunc=None, hyperpars={}, pars={}, name='prior'):

        # Posterior hyper-parameters
        self.hyperpars = []
        for key in hyperpars.keys():
            val = hyperpars[key]
            param = Parameter(name, key, val, alwaysfrozen=True)
            self.__dict__[key] = param
            self.hyperpars.append(param)

        # References to parameters in source model
        self.pars = []
        for key in pars.keys():
            self.__dict__[key] = pars[key]
            self.pars.append(pars[key])

        self._statfuncset = False
        self.statfunc = (lambda x: None)
        if statfunc is not None:
            self.statfunc = statfunc
            self._statfuncset = True

        Likelihood.__init__(self, name)


    def __str__(self):

        s = self.name
        hfmt = '\n   %-15s %-6s %12s'
        s += hfmt % ('Param', 'Type', 'Value')
        s += hfmt % ('-'*5, '-'*4, '-'*5)
	for p in self.hyperpars:
            s += ('\n   %-15s %-6s %12g' %
                  (p.fullname,'frozen', p.val))
        return s


    def set_statfunc(self, func):
        self.statfunc = func
        self._statfuncset = True

    def calc_stat(self, data, model, staterror=None, syserror=None,
                  weight=None):
        if not self._statfuncset:
            raise StatErr('nostat', self.name, 'calc_stat()')
        return self.statfunc(self, data, model, staterror, syserror, weight)
