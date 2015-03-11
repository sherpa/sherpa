# 
#  Copyright (C) 2009  Smithsonian Astrophysical Observatory
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
from sherpa.utils import NoNewAttributesAfterInit
from sherpa.utils.err import StatErr
import sherpa.stats._statfcts


__all__ = ('Stat', 'Cash', 'CStat', 'LeastSq',
           'Chi2Gehrels', 'Chi2ConstVar', 'Chi2DataVar', 'Chi2ModVar',
           'Chi2XspecVar', 'Chi2')


from sherpa import get_config
from ConfigParser import ConfigParser

config = ConfigParser()
config.read(get_config())

# truncation_flag indicates whether or not model truncation
# should be performed.  If true, use the truncation_value from
# the config file.
truncation_flag = config.get('statistics','truncate').upper()
truncation_value = float(config.get('statistics','trunc_value'))
if (bool(truncation_flag) is False or truncation_flag == "FALSE" or
    truncation_flag == "NONE" or truncation_flag == "0"):
    truncation_value = -1.0


class Stat(NoNewAttributesAfterInit):

    def __init__(self, name):
        self.name = name
        NoNewAttributesAfterInit.__init__(self)

    def __repr__(self):
        if self.__doc__ is not None:
            return self.__doc__
        return ("<%s statistic instance '%s'>" %
                (type(self).__name__, self.name))

    def calc_staterror(self, data):
        raise NotImplementedError

    def calc_stat(self, data, model, staterror=None, syserror=None,
                  weight=None):
        raise NotImplementedError

class Likelihood(Stat):
    """Maximum likelihood function"""
    def __init__(self, name='likelihood'):
        Stat.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        # Likelihood stats do not have 'errors' associated with them.
        # return 1 to avoid dividing by 0 by some optimization methods.
        return numpy.ones_like(data)


class Cash(Likelihood):
    """Maximum likelihood function"""
    def __init__(self, name='cash'):
        Likelihood.__init__(self, name)

    @staticmethod
    def calc_stat(data, model, staterror=None, syserror=None, weight=None):
        return _statfcts.calc_cash_stat(data, model, staterror, syserror,
                                        weight, truncation_value)


class CStat(Likelihood):
    """Maximum likelihood function (XSPEC style)"""
    def __init__(self, name='cstat'):
        Likelihood.__init__(self, name)

    @staticmethod
    def calc_stat(data, model, staterror=None, syserror=None, weight=None):
        return _statfcts.calc_cstat_stat(data, model, staterror, syserror,
                                         weight, truncation_value)


class Chi2(Stat):
    """Chi Squared"""
    def __init__(self, name='chi2'):
        Stat.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
         raise StatErr('chi2noerr')

    @staticmethod
    def calc_stat(data, model, staterror, syserror=None, weight=None):
        return _statfcts.calc_chi2_stat(data, model, staterror,
                                        syserror, weight, truncation_value)

class LeastSq(Chi2):
    """Least Squared"""
    def __init__(self, name='leastsq'):
        Stat.__init__(self, name)

    @staticmethod
    def calc_staterror(data):
        return numpy.ones_like(data)        

    @staticmethod
    def calc_stat(data, model, staterror, syserror=None, weight=None):
        return _statfcts.calc_lsq_stat(data, model, staterror,
                                       syserror, weight, truncation_value)
    

class Chi2Gehrels(Chi2):
    """Chi Squared with Gehrels variance"""
    def __init__(self, name='chi2gehrels'):
        Chi2.__init__(self, name)

    calc_staterror = _statfcts.calc_chi2gehrels_errors


class Chi2ConstVar(Chi2):
    """Chi Squared with constant variance"""
    def __init__(self, name='chi2constvar'):
        Chi2.__init__(self, name)

    calc_staterror = _statfcts.calc_chi2constvar_errors


class Chi2DataVar(Chi2):
    """Chi Squared with data variance"""
    def __init__(self, name='chi2datavar'):
        Chi2.__init__(self, name)

    calc_staterror = _statfcts.calc_chi2datavar_errors


class Chi2ModVar(Chi2):
    """Chi Squared with model amplitude variance"""
    def __init__(self, name='chi2modvar'):
        Chi2.__init__(self, name)

    # Statistical errors are not used
    @staticmethod
    def calc_staterror(data):
        return numpy.zeros_like(data)

    @staticmethod
    def calc_stat(data, model, staterror, syserror=None, weight=None):
        return _statfcts.calc_chi2modvar_stat(data, model, staterror,
                                              syserror, weight,
                                              truncation_value)


class Chi2XspecVar(Chi2):
    """Chi Squared with data variance (XSPEC style)"""
    def __init__(self, name='chi2xspecvar'):
        Chi2.__init__(self, name)

    calc_staterror = _statfcts.calc_chi2xspecvar_errors


class UserStat(Stat):

    def __init__(self, statfunc=None, errfunc=None, name='userstat'):
        self._statfuncset = False
        self.statfunc = (lambda x: None)
        
        self._staterrfuncset = False
        self.errfunc = (lambda x: None)
        
        if statfunc is not None:
            self.statfunc = statfunc
            self._statfuncset = True

        if errfunc is not None:
            self.errfunc = errfunc
            self._staterrfuncset = True
            

        Stat.__init__(self, name)


    def __getstate__(self):
        state = self.__dict__.copy()
        # Function pointers to methods of the class
        # (of type 'instancemethod') are NOT picklable
        # remove them and restore later with a coord init
        del state['statfunc']
        del state['errfunc']

        return state

    def __setstate__(self, state):
        # Populate the function pointers we deleted at pickle time with
        # no-ops.
        self.__dict__['statfunc']=(lambda x: None)
        self.__dict__['errfunc']=(lambda x: None)
        self.__dict__.update(state)


    def set_statfunc(self, func):
        self.statfunc = func
        self._statfuncset = True


    def set_errfunc(self, func):
        self.errfunc = func
        self._staterrfuncset = True


    def calc_staterror(self, data):
        if not self._staterrfuncset:
            raise StatErr('nostat', self.name, 'calc_staterror()')
        return self.errfunc(data)


    def calc_stat(self, data, model, staterror=None, syserror=None,
                  weight=None):
        if not self._statfuncset:
            raise StatErr('nostat', self.name, 'calc_stat()')
        return self.statfunc(data, model, staterror, syserror, weight)
