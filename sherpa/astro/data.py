# 
#  Copyright (C) 2008  Smithsonian Astrophysical Observatory
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

"""
Classes for storing, inspecting, and manipulating astronomical data sets
"""

from itertools import izip
import os.path
import numpy
from sherpa.data import BaseData, Data1DInt, Data2D, DataND
from sherpa.utils.err import DataErr, ImportErr
from sherpa.utils import SherpaFloat, pad_bounding_box, interpolate, \
    create_expr, parse_expr, bool_cast, rebin, filter_bins
from sherpa.astro.utils import *

import logging
warning = logging.getLogger(__name__).warning

groupstatus = False
try:
    import group as pygroup
    groupstatus = True
except:
    groupstatus = False
    warning('the group module (from the CIAO tools package) is not installed.'
            + '\nDynamic grouping functions will not be available.')


__all__ = ('DataARF', 'DataRMF', 'DataPHA', 'DataIMG', 'DataIMGInt')

def _notice_resp(chans, arf, rmf):
    bin_mask=None

    if rmf is not None and arf is not None:

        bin_mask = rmf.notice(chans)        
        if len(rmf.energ_lo) == len(arf.energ_lo):
            arf.notice(bin_mask)

        # If the response is mis-matched, determine which energy bins in the RMF
        # coorespond to energy bins in the ARF and which are noticed.
        # Propogate the noticed RMF energy bins to the ARF energy  bins.
        elif len(rmf.energ_lo) < len(arf.energ_lo):
            arf_mask = None
            if bin_mask is not None:
                arf_mask = numpy.zeros(len(arf.energ_lo), dtype=bool)
                for ii, val in enumerate(bin_mask):
                    if val:
                        slice = filter_bins( (rmf.energ_lo[ii],), (rmf.energ_hi[ii],),
                                             (arf.energ_lo,arf.energ_hi) ).nonzero()[0]
                        arf_mask[slice]=True
            arf.notice(arf_mask)

    else:
        if rmf is not None:
            bin_mask = rmf.notice(chans)
        if arf is not None:
            arf.notice(bin_mask)


class DataARF(Data1DInt):
    "ARF data set"

    mask = property(BaseData._get_mask, BaseData._set_mask)

    def _get_specresp(self):
        return self._specresp

    def _set_specresp(self, val):
        self._specresp = val
        self._rsp = val

    specresp = property(_get_specresp, _set_specresp)

    def __init__(self, name, energ_lo, energ_hi, specresp, bin_lo=None,
                 bin_hi=None, exposure=None, header=None):
        self._lo = energ_lo
        self._hi = energ_hi
        BaseData.__init__(self)

    def __str__(self):
        # Print the metadata first
        old = self._fields
        ss = old
        try:
            self._fields = filter((lambda x: x!='header'), self._fields)
            ss = BaseData.__str__(self)
        finally:
            self._fields = old
        return ss


    def __setstate__(self, state):
        if not state.has_key('header'):
            self.header=None
        self.__dict__.update(state)
        
        if not state.has_key('_specresp'):
            self.__dict__['_specresp'] = state.get('specresp',None)
            self.__dict__['_rsp'] = state.get('specresp',None)


    def apply_arf(self, src, *args, **kwargs):
        "Fold the source array src through the ARF and return the result"
        model = arf_fold(src, self._rsp)

        # Rebin the high-res source model folded through ARF down to the size
        # the PHA or RMF expects.
        if args != ():
            (arf, rmf) = args
            if rmf != () and len(arf[0]) > len(rmf[0]):
                model = rebin(model, arf[0], arf[1], rmf[0], rmf[1])

        return model


    def notice(self, bin_mask=None):
        self._rsp = self.specresp
        self._lo = self.energ_lo
        self._hi = self.energ_hi
        if bin_mask is not None:
            self._rsp = self.specresp[bin_mask]
            self._lo = self.energ_lo[bin_mask]
            self._hi = self.energ_hi[bin_mask]

    def get_indep(self, filter=False):
        filter=bool_cast(filter)
        return (self._lo, self._hi)

    def get_dep(self, filter=False):
        filter=bool_cast(filter)
        return self._rsp

    def get_xlabel(self):
        return 'Energy (keV)'

    def get_ylabel(self):
        return 'cm^2'

class DataRMF(Data1DInt):
    "RMF data set"

    mask = property(BaseData._get_mask, BaseData._set_mask)

    def __init__(self, name, detchans, energ_lo, energ_hi, n_grp, f_chan,
                 n_chan, matrix, offset=1, e_min=None, e_max=None, header=None):
        self._fch = f_chan;  self._nch = n_chan
        self._grp = n_grp;   self._rsp = matrix
        self._lo = energ_lo; self._hi = energ_hi
        BaseData.__init__(self)

    def __str__(self):
        # Print the metadata first
        old = self._fields
        ss = old
        try:
            self._fields = filter((lambda x: x!='header'), self._fields)
            ss = BaseData.__str__(self)
        finally:
            self._fields = old
        return ss


    def __setstate__(self, state):
        if not state.has_key('header'):
            self.header=None
        self.__dict__.update(state)


    def apply_rmf(self, src, *args, **kwargs):
        "Fold the source array src through the RMF and return the result"

        # Rebin the high-res source model from the PHA down to the size
        # the RMF expects.
        if args != ():
            (rmf, pha) = args
            if pha != () and len(pha[0]) > len(rmf[0]):
                src = rebin(src, pha[0], pha[1], rmf[0], rmf[1])

        if len(src) != len(self._lo):
            raise TypeError("Mismatched filter between ARF and RMF or PHA and RMF")

        return rmf_fold(src, self._grp, self._fch, self._nch, self._rsp,
                        self.detchans, self.offset)

    def notice(self, noticed_chans=None):
        bin_mask=None
        self._fch = self.f_chan;  self._nch = self.n_chan
        self._grp = self.n_grp;   self._rsp = self.matrix
        self._lo = self.energ_lo; self._hi = self.energ_hi
        if noticed_chans is not None:
            (self._grp, self._fch, self._nch, self._rsp,
             bin_mask) = filter_resp(noticed_chans, self.n_grp, self.f_chan,
                                     self.n_chan, self.matrix, self.offset)
            self._lo = self.energ_lo[bin_mask]
            self._hi = self.energ_hi[bin_mask]

        return bin_mask

    #def get_indep(self):
    #    if (self.e_min is not None) and (self.e_max is not None):
    #        return (self.e_min, self.e_max)
    #    channels = numpy.arange(1.0, self.detchans+1.0, 1.0, SherpaFloat)
    #    return (channels - 0.5, channels + 0.5)

    def get_indep(self, filter=False):
        filter=bool_cast(filter)
	return (self._lo, self._hi)

    def get_dep(self, filter=False):
        filter=bool_cast(filter)
        return self.apply_rmf(numpy.ones(self.energ_lo.shape, SherpaFloat))

    def get_xlabel(self):
        if (self.e_min is not None) and (self.e_max is not None):
            return 'Energy (keV)'
        return 'Channel'

    def get_ylabel(self):
        return 'Counts'


class DataPHA(Data1DInt):
    "PHA data set, including any associated instrument and background data"

    mask = property(BaseData._get_mask, BaseData._set_mask)

    def _get_grouped(self):
        return self._grouped
    def _set_grouped(self, val):
        val = bool(val)

        if val and (self.grouping is None):
            raise DataErr('nogrouping', self.name)

        # If grouping status is being changed, we need to reset the mask
        # to be correct size, while still noticing groups within the filter
        if self._grouped != val:
            do_notice = numpy.iterable(self.mask)
            if do_notice:
                old_filter = self.get_filter(val)
                self._grouped = val
                self.ignore()
                for vals in parse_expr(old_filter):
                    self.notice(*vals)
            #self.mask = True

        self._grouped = val
    grouped = property(_get_grouped, _set_grouped, doc='Are the data grouped?')

    def _get_subtracted(self):
        return self._subtracted
    def _set_subtracted(self, val):
        val = bool(val)
        if len(self._backgrounds) == 0:
            raise DataErr('nobkg', self.name)
        self._subtracted = val
    subtracted = property(_get_subtracted, _set_subtracted,
                          doc='Are the background data subtracted?')

    def _get_units(self):
        return self._units
    def _set_units(self, val):
        units = str(val).strip().lower()

        if units == 'bin':
            units = 'channel'

        if units.startswith('chan'):
            self._to_channel   = (lambda x, group=True, response_id=None: x)
            self._from_channel = (lambda x, group=True, response_id=None: x)
            units = 'channel'

        elif units.startswith('ener'):
            self._to_channel   = self._energy_to_channel
            self._from_channel = self._channel_to_energy
            units = 'energy'

        elif units.startswith('wave'):
            self._to_channel   = self._wavelength_to_channel
            self._from_channel = self._channel_to_wavelength
            units = 'wavelength'

        else:
            raise DataErr('bad', 'quantity', val)

        for id in self.background_ids:
            bkg = self.get_background(id)
            if (bkg.get_response() != (None,None) or
                (bkg.bin_lo is not None and bkg.bin_hi is not None)):
                bkg.units = units

        self._units = units

    units = property(_get_units, _set_units, doc='Units of independent axis')

    def _get_rate(self):
        return self._rate
    def _set_rate(self, val):
        self._rate = bool_cast(val)
        for id in self.background_ids:
            self.get_background(id).rate = val

    rate = property(_get_rate, _set_rate, doc='Quantity of y-axis: counts or counts/sec')

    def _get_plot_fac(self):
        return self._plot_fac
    def _set_plot_fac(self, val):
        self._plot_fac = int(val)
        for id in self.background_ids:
            self.get_background(id).plot_fac = val

    plot_fac = property(_get_plot_fac, _set_plot_fac, doc='Number of times to multiply' +
                        ' the y-axis quantity by x-axis bin size')

    def _get_response_ids(self):
        return self._response_ids

    def _set_response_ids(self, ids):
        if not numpy.iterable(ids):
            raise DataErr('idsnotarray', 'response', str(ids))
        keys = self._responses.keys()
        for id in ids:
            if id not in keys:
                raise DataErr('badids', str(id), 'response', str(keys))
        ids = list(ids)
        self._response_ids = ids

    response_ids = property(_get_response_ids, _set_response_ids,
                            doc=('IDs of defined instrument responses ' +
                                 '(ARF/RMF pairs)'))

    def _get_background_ids(self):
        return self._background_ids

    def _set_background_ids(self, ids):
        if not numpy.iterable(ids):
            raise DataErr('idsnotarray', 'background', str(ids))
        keys = self._backgrounds.keys()
        for id in ids:
            if id not in keys:
                raise DataErr('badids', str(id), 'background', str(keys))
        ids = list(ids)
        self._background_ids = ids

    background_ids = property(_get_background_ids, _set_background_ids,
                              doc='IDs of defined background data sets')

    _extra_fields = ('grouped', 'subtracted', 'units', 'rate', 'plot_fac',
                     'response_ids','background_ids')

    def __init__(self, name, channel, counts, staterror=None, syserror=None,
                 bin_lo=None, bin_hi=None, grouping=None, quality=None,
                 exposure=None, backscal=None, areascal=None, header=None):
        self._grouped = (grouping is not None)
        self._original_groups = True
        self._subtracted = False
        self._response_ids = []
        self._background_ids = []
        self._responses   = {}
        self._backgrounds = {}
        self._rate=True
        self._plot_fac=0
        self.units = 'channel'
        self.quality_filter = None
        BaseData.__init__(self)

    def __str__(self):
        # Print the metadata first
        old = self._fields
        ss = old
        try:
            self._fields = filter((lambda x: x!='header'), self._fields)
            ss = BaseData.__str__(self)
        finally:
            self._fields = old
        return ss


    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_to_channel']
        del state['_from_channel']
        return state

    def __setstate__(self, state):
        self._background_ids = state['_background_ids']
        self._backgrounds = state['_backgrounds']
        self._set_units(state['_units'])

        if not state.has_key('header'):
            self.header=None
        self.__dict__.update(state)

 
    primary_response_id = 1


    def set_analysis(self, quantity, type='rate', factor=0):
        self.plot_fac = factor

        type = str(type).strip().lower()
        if not (type.startswith('counts') or type.startswith('rate')):
            raise DataErr("plottype", type, "'rate' or 'counts'")

        self.rate = (type=='rate')

        arf, rmf = self.get_response()
        if rmf is not None and rmf.detchans != len(self.channel):
            raise DataErr("incompatibleresp", rmf.name, self.name)

        if ( (rmf is None and arf is None) and
             (self.bin_lo is None and self.bin_hi is None) and quantity != 'channel'):
            raise DataErr('noinstr', self.name)

        if (rmf is None and arf is not None and quantity != 'channel' and 
            len(arf.energ_lo) != len(self.channel)):
            raise DataErr("incompleteresp", self.name)

        self.units = quantity

    def get_analysis(self):
        return self.units

    def _fix_response_id(self, id):
        if id is None:
            id = self.primary_response_id
        return id

    def get_response(self, id=None):
        id = self._fix_response_id(id)
        return self._responses.get(id, (None, None))

    def set_response(self, arf=None, rmf=None, id=None):
        if (arf is None) and (rmf is None):
            return
        # Multiple responses re-enabled for CIAOX
        id = self._fix_response_id(id)
        self._responses[id] = (arf, rmf)
        ids = self.response_ids[:]
        if id not in ids:
            ids.append(id)
        self.response_ids = ids

    def delete_response(self, id=None):
        id = self._fix_response_id(id)
        self._responses.pop(id, None)
        ids = self.response_ids[:]
        ids.remove(id)
        self.response_ids=ids

    def get_arf(self, id=None):
        return self.get_response(id)[0]

    def get_rmf(self, id=None):
        return self.get_response(id)[1]

    def set_arf(self, arf, id=None):
        self.set_response(arf, self.get_rmf(id), id)

    def set_rmf(self, rmf, id=None):
        self.set_response(self.get_arf(id), rmf, id)


    def get_specresp(self, filter=False):
        filter=bool_cast(filter)
        self.notice_response(False)
        arf,rmf = self.get_response()
        newarf=None

        if arf is not None and rmf is not None:
            specresp = arf.get_dep()
            elo, ehi = arf.get_indep()
            lo, hi = self._get_ebins(group=False)

            newarf = interpolate(lo, elo, specresp)
            newarf[newarf<=0]=1.

            if filter:
                newarf = self.apply_filter(newarf, self._middle)

        return newarf

    # The energy bins can be grouped or ungrouped.  By default,
    # they should be grouped if the data are grouped.  There are
    # certain contexts (e.g., plotting) where we will retrieve the
    # energy bins, and later filter the data; but filtering
    # is automatically followed by grouping.  Grouping the data
    # twice is an error.
    def _get_ebins(self, response_id=None, group=True):
        group=bool_cast(group)
        arf, rmf = self.get_response(response_id)
        if (self.bin_lo is not None) and (self.bin_hi is not None):
            elo = self.bin_lo
            ehi = self.bin_hi
            if (elo[0] > elo[-1]) and (ehi[0] > ehi[-1]):
                elo = self._hc / self.bin_hi
                ehi = self._hc / self.bin_lo
        elif rmf is not None:
            if (rmf.e_min is None) or (rmf.e_max is None):
                raise DataErr('noenergybins', 'RMF')
            elo = rmf.e_min
            ehi = rmf.e_max
        elif arf is not None:
            elo = arf.energ_lo
            ehi = arf.energ_hi
        else:
            elo = self.channel - 0.5
            ehi = self.channel + 0.5

        if self.units == 'channel':
            elo = self.channel - 0.5
            ehi = self.channel + 0.5

        # If the data are grouped, then we should group up
        # the energy bins as well.  E.g., if group 1 is
        # channels 1-5, then the energy boundaries for the
        # *group* should be elo[0], ehi[4].
        if (self.grouped and group):
            elo = self.apply_grouping(elo, self._min)
            ehi = self.apply_grouping(ehi, self._max)

        return (elo, ehi)


    def get_indep(self, filter=True):
        if filter:
            return (self.get_noticed_channels(),)

        return (self.channel,)


    def _get_indep(self, filter=False):
        if (self.bin_lo is not None) and (self.bin_hi is not None):
            elo = self.bin_lo
            ehi = self.bin_hi
            if (elo[0] > elo[-1]) and (ehi[0] > ehi[-1]):
                if self.units == 'wavelength':
                    return (elo, ehi)

                elo = self._hc / self.bin_hi
                ehi = self._hc / self.bin_lo

        else:
            energylist = []
            for id in self.response_ids:
                arf, rmf = self.get_response(id)
                lo = None; hi = None

                if rmf is not None:
                    lo = rmf.energ_lo
                    hi = rmf.energ_hi
                    if filter:
                        lo, hi = rmf.get_indep()

                elif arf is not None:
                    lo = arf.energ_lo
                    hi = arf.energ_hi
                    if filter:
                        lo, hi = arf.get_indep()

                energylist.append((lo,hi))

            if len(energylist) > 1:
                elo, ehi, lookuptable = compile_energy_grid(energylist)
            elif (not energylist or
                  ( len(energylist) == 1 and
                    energylist[0] == (None,None) )):
                raise DataErr('noenergybins', 'Response')
            else:
                elo, ehi = energylist[0]

        lo, hi = elo, ehi
        if self.units == 'wavelength':
            lo = self._hc/ehi
            hi = self._hc/elo

        return (lo,hi)


    def _channel_to_energy(self, val, group=True, response_id=None):
        elo, ehi = self._get_ebins(response_id=response_id, group=group)
        val = numpy.asarray(val).astype(numpy.int_) - 1
        try:
            return (elo[val] + ehi[val]) / 2.0
        except IndexError:
            raise DataErr('invalidchannel', val)

    def _energy_to_channel(self, val):
        elo, ehi = self._get_ebins()

        val = numpy.asarray(val)
        res = []
        for v in val.flat:
            if tuple(numpy.flatnonzero(elo <= v)) == ():
                if elo[0] > elo[-1] and  ehi[0] > ehi[-1]:
                    res.append(SherpaFloat(len(elo)))
                else:
                    res.append(SherpaFloat(1))
            elif tuple(numpy.flatnonzero(ehi > v)) == ():
                if elo[0] > elo[-1] and  ehi[0] > ehi[-1]:
                    res.append(SherpaFloat(1))
                else:
                    res.append(SherpaFloat(len(ehi)))
            elif tuple(numpy.flatnonzero((elo <= v) & (ehi > v) ) + 1) != ():
                res.append(SherpaFloat(
                        numpy.flatnonzero((elo <= v) & (ehi > v) ) + 1))
            elif (elo<=v).argmin() == (ehi>v).argmax():
                res.append(SherpaFloat((elo<=v).argmin()))
            else:
                raise DataErr("energytochannel", v)

        if val.shape == ():
            return res[0]

        return numpy.asarray(res, SherpaFloat)

    _hc = 12.39841874  # nist.gov in [keV-Angstrom]

    def _channel_to_wavelength(self, val, group=True, response_id=None):
        tiny = numpy.finfo(numpy.float32).tiny
        vals = numpy.asarray(self._channel_to_energy(val, group, response_id))
        if vals.shape == ():
            if vals == 0.0:  vals = tiny
        else:
            vals[ vals == 0.0 ] = tiny
        vals = self._hc / vals
        return vals

    def _wavelength_to_channel(self, val):
        tiny = numpy.finfo(numpy.float32).tiny
        vals = numpy.asarray(val)
        if vals.shape == ():
            if vals == 0.0:  vals = tiny
        else:
            vals[ vals == 0.0 ] = tiny
        vals = self._hc / vals
        return self._energy_to_channel(vals)

    default_background_id = 1

    def _fix_background_id(self, id):
        if id is None:
            id = self.default_background_id
        return id

    def get_background(self, id=None):
        id = self._fix_background_id(id)
        return self._backgrounds.get(id)

    def set_background(self, bkg, id=None):
        id = self._fix_background_id(id)
        self._backgrounds[id] = bkg
        ids = self.background_ids[:]
        if id not in ids:
            ids.append(id)
        self.background_ids = ids

    def delete_background(self, id=None):
        id = self._fix_background_id(id)
        self._backgrounds.pop(id, None)
        if len(self._backgrounds) == 0:
            self._subtracted = False            
        ids = self.background_ids[:]
        if id in ids:
            ids.remove(id)
        self.background_ids = ids


    def get_background_scale(self):
        if len(self.background_ids) == 0:
            return None
        return self.sum_background_data(lambda key, bkg: 1.)

    def _check_scale(self, scale, group=True, filter=False):
        if numpy.isscalar(scale) and scale <= 0.0:
            scale = 1.0
        elif numpy.iterable(scale):
            scale = numpy.asarray(scale, dtype=SherpaFloat)
            if group:
                if filter:
                    scale = self.apply_filter(scale, self._middle)
                else:
                    scale = self.apply_grouping(scale, self._middle)

            scale[scale<=0.0] = 1.0
        return scale

    def get_backscal(self, group=True, filter=False):
        backscal = self.backscal
        if backscal is not None:
            backscal = self._check_scale(backscal, group, filter)
        return backscal

    def get_areascal(self, group=True, filter=False):
        areascal = self.areascal
        if areascal is not None:
            areascal = self._check_scale(areascal, group, filter)
        return areascal

    def apply_filter(self, data, groupfunc=numpy.sum):
        """

        Filter the array data, first passing it through apply_grouping()
        (using groupfunc) and then applying the general filters
        
        """
        if (data is None):
            return data
        elif len(data) != len(self.counts):
            counts = numpy.zeros(len(self.counts), dtype=SherpaFloat)
            mask = self.get_mask()
            if mask is not None:
                counts[mask] = numpy.asarray(data, dtype=SherpaFloat)
                data = counts
#            else:
#                raise DataErr('mismatch', "filter", "data array")
        return Data1DInt.apply_filter(self,
                                      self.apply_grouping(data, groupfunc))

    def apply_grouping(self, data, groupfunc=numpy.sum):
        """

        Apply the data set's grouping scheme to the array data,
        combining the grouped data points with groupfunc, and return
        the grouped array.  If the data set has no associated grouping
        scheme or the data are ungrouped, data is returned unaltered.

        """
        if (data is None) or (not self.grouped):
            return data

        groups = self.grouping
        filter = self.quality_filter
        if filter is None:
            return do_group( data, groups, groupfunc.__name__ )

        if (len(data) != len(filter) or
            len(groups) != len(filter)):
            raise DataErr('mismatch', "quality filter", "data array")
        
        filtered_data = numpy.asarray(data)[filter]
        groups = numpy.asarray(groups)[filter]
        grouped_data = do_group( filtered_data, groups, groupfunc.__name__ )
        
        if data is self.channel and groupfunc is self._make_groups:
            return numpy.arange(1, len(grouped_data)+1, dtype=int)

        return grouped_data

    def ignore_bad(self):
        if self.quality is None:
            raise DataErr("noquality", self.name)

        qual_flags = ~numpy.asarray(self.quality, bool)

        if self.grouped and (not self.mask is True):
            self.notice()
            warning('filtering grouped data with quality flags,' +
                    ' previous filters deleted' )

        elif not self.grouped:
            # if ungrouped, create/combine with self.mask
            if not self.mask is True:
                self.mask = self.mask & qual_flags
                return
            else:
                self.mask = qual_flags
                return

        # self.quality_filter used for pre-grouping filter
        self.quality_filter = qual_flags


    def _dynamic_group(self, group_func, *args, **kwargs):

        keys = kwargs.keys()[:]
        for key in keys:
            if kwargs[key] is None:
                kwargs.pop(key)

        old_filter = self.get_filter(group=False)
        do_notice = numpy.iterable(self.mask)

        self.grouping, self.quality = group_func(*args, **kwargs)
        self.group()
        self._original_groups = False

        if do_notice:
            # self.group() above has cleared the filter if applicable
            # No, that just sets a flag.  So manually clear filter
            # here
            self.ignore()
            for vals in parse_expr(old_filter):
                self.notice(*vals)

        #warning('grouping flags have changed, noticing all bins')

    # Have to move this check here; as formerly written, reference
    # to pygroup functions happened *before* checking groupstatus,
    # in _dynamic_group.  So we did not return the intended error
    # message; rather, a NameError was raised stating that pygroup
    # did not exist in global scope (not too clear to the user).
    #
    # The groupstatus check thus has to be done in *each* of the following
    # group functions.

    ## Dynamic grouping functions now automatically impose the
    ## same grouping conditions on *all* associated background data sets.
    ## CIAO 4.5 bug fix, 05/01/2012
    def group_bins(self, num, tabStops=None):
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpNumBins, len(self.channel), num,
                            tabStops=tabStops)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id) 
            if (hasattr(bkg, "group_bins")):
                bkg.group_bins(num, tabStops=tabStops)

    def group_width(self, val, tabStops=None):
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpBinWidth, len(self.channel), val,
                            tabStops=tabStops)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id) 
            if (hasattr(bkg, "group_width")):
                bkg.group_width(val, tabStops=tabStops)

    def group_counts(self, num, maxLength=None, tabStops=None):
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpNumCounts, self.counts, num,
                            maxLength=maxLength, tabStops=tabStops)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id) 
            if (hasattr(bkg, "group_counts")):
                bkg.group_counts(num, maxLength=maxLength, tabStops=tabStops)

    def group_snr(self, snr, maxLength=None, tabStops=None, errorCol=None):
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpSnr, self.counts, snr,
                            maxLength=maxLength, tabStops=tabStops,
                            errorCol=errorCol)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id) 
            if (hasattr(bkg, "group_snr")):
                bkg.group_snr(snr, maxLength=maxLength, tabStops=tabStops, errorCol=errorCol)

    def group_adapt(self, minimum, maxLength=None, tabStops=None):
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpAdaptive, self.counts, minimum,
                            maxLength=maxLength, tabStops=tabStops)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id) 
            if (hasattr(bkg, "group_adapt")):
                bkg.group_adapt(minimum, maxLength=maxLength, tabStops=tabStops)

    def group_adapt_snr(self, minimum, maxLength=None, tabStops=None, errorCol=None):
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpAdaptiveSnr, self.counts, minimum,
                            maxLength=maxLength, tabStops=tabStops,
                            errorCol=errorCol)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id) 
            if (hasattr(bkg, "group_adapt_snr")):
                bkg.group_adapt_snr(minimum, maxLength=maxLength, tabStops=tabStops, errorCol=errorCol)

    def eval_model(self, modelfunc):
        return modelfunc(*self.get_indep(filter=False))

    def eval_model_to_fit(self, modelfunc):
        return self.apply_filter(modelfunc(*self.get_indep(filter=True)))


    def sum_background_data(self,
                            get_bdata_func=(lambda key, bkg: bkg.counts)):
        bdata_list = []

        #for key, bkg in self._backgrounds.items():
        for key in self.background_ids:
            bkg = self.get_background(key)
            bdata = get_bdata_func(key, bkg)

            backscal = bkg.backscal
            if backscal is not None:
                backscal = self._check_scale(backscal, group=False)
                bdata = bdata / backscal

            areascal = bkg.areascal
            if areascal is not None:
                areascal = self._check_scale(areascal, group=False)
                bdata = bdata / areascal

            if bkg.exposure is not None:
                bdata = bdata / bkg.exposure

            bdata_list.append(bdata)

        nbkg = len(bdata_list)
        assert (nbkg > 0)
        if nbkg == 1:
            bkgsum = bdata_list[0]
        else:
            bkgsum = sum(bdata_list)

        backscal = self.backscal
        if backscal is not None:
            backscal = self._check_scale(backscal, group=False)
            bkgsum = backscal * bkgsum

        areascal = self.areascal
        if areascal is not None:
            areascal = self._check_scale(areascal, group=False)
            bkgsum = areascal * bkgsum

        if self.exposure is not None:
            bkgsum = self.exposure * bkgsum

        return bkgsum / SherpaFloat(nbkg)

    def get_dep(self, filter=False):
        # FIXME: Aneta says we need to group *before* subtracting, but that
        # won't work (I think) when backscal is an array
        #if not self.subtracted:
        #    return self.counts
        #return self.counts - self.sum_background_data()
        dep = self.counts
        filter=bool_cast(filter)
        if self.subtracted:
            bkg = self.sum_background_data()
            if len(dep) != len(bkg):
                raise DataErr("subtractlength")
            dep = dep - bkg
        if filter:
            dep = self.apply_filter(dep)
        return dep

    def set_dep(self, val):
        dep = None
        if numpy.iterable(val):
            dep = numpy.asarray(val, SherpaFloat)
        else:
            val = SherpaFloat(val)
            dep = numpy.array([val]*len(self.get_indep()[0]))
        setattr(self, 'counts', dep)

    def get_staterror(self, filter=False, staterrfunc=None):
        staterr = self.staterror
        filter=bool_cast(filter)
        if filter:
            staterr = self.apply_filter(staterr, self._sum_sq)
        else:
            staterr = self.apply_grouping(staterr, self._sum_sq)
            
        if (staterr is None) and (staterrfunc is not None):
            cnts = self.counts
            if filter:
                cnts = self.apply_filter(cnts)
            else:
                cnts = self.apply_grouping(cnts)
            
            staterr = staterrfunc(cnts)

        if (staterr is not None) and self.subtracted:
            bkg_staterr_list = []

            #for bkg in self._backgrounds.values():
            for key in self.background_ids:
                bkg = self.get_background(key)
                berr = bkg.staterror
                if filter:
                    berr = self.apply_filter(berr, self._sum_sq)
                else:
                    berr = self.apply_grouping(berr, self._sum_sq)
                    
                if (berr is None) and (staterrfunc is not None):
                    bkg_cnts = bkg.counts
                    if filter:
                        bkg_cnts = self.apply_filter(bkg_cnts)
                    else:
                        bkg_cnts = self.apply_grouping(bkg_cnts)
                        
                    if (hasattr(staterrfunc,'__name__') and
                        staterrfunc.__name__ == 'calc_chi2datavar_errors' and
                        0.0 in bkg_cnts):
                        mask = (numpy.asarray(bkg_cnts)!=0.0)
                        berr = numpy.zeros(len(bkg_cnts))
                        berr[mask] = staterrfunc(bkg_cnts[mask])
                    else:
                        berr = staterrfunc(bkg_cnts)

                # FIXME: handle this
                # assert (berr is not None)

                # This case appears when the source dataset has an error
                # column and at least one of the background(s) do not.
                # Because the staterr is not None and staterrfunc is, I think
                # we should return None.  This way the user knows to call with
                # staterrfunc next time.
                if berr is None:
                    return None

                bksl = bkg.backscal
                if bksl is not None:
                    bksl = self._check_scale(bksl, filter=filter)
                    berr = berr / bksl
                
                area = bkg.areascal
                if area is not None:
                    area = self._check_scale(area, filter=filter)
                    berr = berr / area

                if bkg.exposure is not None:
                    berr = berr / bkg.exposure

                berr = berr * berr
                bkg_staterr_list.append(berr)

            nbkg = len(bkg_staterr_list)
            assert (nbkg > 0)
            if nbkg == 1:
                bkgsum = bkg_staterr_list[0]
            else:
                bkgsum = sum(bkg_staterr_list)

            bscal = self.backscal
            if bscal is not None:
                bscal = self._check_scale(bscal, filter=filter)
                bkgsum = (bscal * bscal) * bkgsum

            area = self.areascal
            if area is not None:
                area = self._check_scale(area, filter=filter)
                bkgsum = (area * area) * bkgsum

            if self.exposure is not None:
                bkgsum = (self.exposure * self.exposure) * bkgsum

            nbkg = SherpaFloat(nbkg)

            if staterr is not None:
                staterr = staterr*staterr + bkgsum / (nbkg * nbkg)
                staterr = numpy.sqrt(staterr)

        return staterr

    def get_syserror(self, filter=False):
        syserr = self.syserror
        filter=bool_cast(filter)
        if filter:
            syserr = self.apply_filter(syserr, self._sum_sq)
        else:
            syserr = self.apply_grouping(syserr, self._sum_sq)
        return syserr
    
    def get_x(self, filter=False, response_id=None):
        # If we are already in channel space, self._from_channel
        # is always ungrouped.  In any other space, we must
        # disable grouping when calling self._from_channel.
        if self.units != 'channel':
            elo,ehi = self._get_ebins(group=False)
            if len(elo) != len(self.channel):
                raise DataErr("incompleteresp", self.name)
            return self._from_channel(self.channel,group=False,
                                      response_id=response_id)
        else:
            return self._from_channel(self.channel)

    def get_xlabel(self):
        xlabel = self.units.capitalize()
        if self.units == 'energy':
            xlabel += ' (keV)'
        elif self.units == 'wavelength':
            xlabel += ' (Angstrom)'
        #elif self.units == 'channel' and self.grouped:
        #    xlabel = 'Group Number'
        return xlabel


    def _set_initial_quantity(self):
        arf, rmf = self.get_response()

        # Change analysis if ARFs equal or of higher resolution to
        # allow for high-res model evaluation.
        if arf is not None and rmf is None:
            if len(arf.energ_lo) == len(self.channel):
                self.units = 'energy'

        # Only change analysis if RMF matches the parent PHA dataset.
        if rmf is not None:
            if len(self.channel) != len(rmf.e_min):
                raise DataErr("incompatibleresp", rmf.name, self.name)
            self.units = 'energy'


    def _fix_y_units(self, val, filter=False, response_id=None):
        if val is None:
            return val

        filter=bool_cast(filter)
        # make a copy of data for units manipulation
        val = numpy.array(val, dtype=SherpaFloat)

        if self.rate and self.exposure:
            val /= self.exposure
            areascal = self.areascal
            if areascal is not None:
                areascal = self._check_scale(areascal, filter=filter)
                val /= areascal

        if self.grouped or self.rate:

            if self.units != 'channel':
                elo, ehi = self._get_ebins(response_id, group=False)
            else:
                elo, ehi = (self.channel, self.channel+1.)

            if filter:
                # If we apply a filter, make sure that
                # ebins are ungrouped before applying
                # the filter.
                elo = self.apply_filter(elo, self._min)
                ehi = self.apply_filter(ehi, self._max)
            elif self.grouped:
                elo = self.apply_grouping(elo, self._min)
                ehi = self.apply_grouping(ehi, self._max)

            if self.units == 'energy':
                ebin = ehi - elo
            elif self.units == 'wavelength':
                ebin = self._hc/elo - self._hc/ehi
            elif self.units == 'channel':
                ebin = ehi - elo
            else:
                raise DataErr("bad", "quantity", self.units)

            val /= numpy.abs(ebin)

        for ii in range(self.plot_fac):
            val *= self.apply_filter(self.get_x(response_id=response_id),
                                     self._middle)

        return val


    def get_y(self, filter=False, yfunc=None, response_id=None):
        vallist = Data1DInt.get_y(self, yfunc=yfunc)
        filter=bool_cast(filter)
        
        if not isinstance(vallist, tuple):
            vallist = (vallist,)

        newvallist = []

        for val in vallist:
            if filter:
                val = self.apply_filter(val)
            else:
                val = self.apply_grouping(val)
            val = self._fix_y_units(val, filter, response_id)
            newvallist.append(val)

        if len(vallist) == 1:
            vallist = newvallist[0]
        else:
            vallist = tuple(newvallist)

        return vallist

    def get_yerr(self, filter=False, staterrfunc=None, response_id=None):
        filter=bool_cast(filter)
        err = self.get_error(filter, staterrfunc)
        return self._fix_y_units(err, filter, response_id)


    def get_xerr(self, filter=False, response_id=None):
        elo, ehi = self._get_ebins(response_id=response_id)
        filter=bool_cast(filter)
        if filter:
            # If we apply a filter, make sure that
            # ebins are ungrouped before applying
            # the filter.
            elo, ehi = self._get_ebins(response_id, group=False)
            elo = self.apply_filter(elo, self._min)
            ehi = self.apply_filter(ehi, self._max)

        return ehi-elo


    def get_ylabel(self):
        ylabel = 'Counts'

        if self.rate and self.exposure:
            ylabel += '/sec'

        if self.rate or self.grouped:
            if self.units == 'energy':
                ylabel += '/keV'
            elif self.units == 'wavelength':
                ylabel += '/Angstrom'
            elif self.units == 'channel':
                ylabel += '/channel'

        if self.plot_fac:
            ylabel += ' X %s^%s' % (self.units.capitalize(), str(self.plot_fac))

        return ylabel

    @staticmethod
    # Dummy function to tell apply_grouping to construct
    # an array of groups.
    def _make_groups(array):
        pass

    @staticmethod
    def _middle(array):
        array = numpy.asarray(array)
        return (array.min() + array.max()) / 2.0

    @staticmethod
    def _min(array):
        array = numpy.asarray(array)
        return array.min()

    @staticmethod
    def _max(array):
        array = numpy.asarray(array)
        return array.max()

    @staticmethod
    def _sum_sq(array):
        return numpy.sqrt(numpy.sum(array * array))

    def get_noticed_channels(self):
        chans = self.channel
        mask = self.get_mask()
        if mask is not None:
            chans = chans[mask]
        return chans

    def get_mask(self):
        groups = self.grouping
        if self.mask is False:
            return None

        if self.mask is True or not self.grouped:
            if self.quality_filter is not None:
                return self.quality_filter
            elif numpy.iterable(self.mask):
                return self.mask
            return None

        if self.quality_filter is not None:
            groups = groups[self.quality_filter]
        return expand_grouped_mask(self.mask, groups)

    def get_noticed_expr(self):
        chans = self.get_noticed_channels()
        if self.mask is False or len(chans) == 0:
            return 'No noticed channels'
        return create_expr(chans, format='%i')

    def get_filter(self, group=True, format = '%.12f', delim=':'):
        """  
        Integrated values returned are measured from center of bin
        """
        if self.mask is False:
            return 'No noticed bins'

        x = self.get_noticed_channels() # ungrouped noticed channels
        if group:
            # grouped noticed channels
            x = self.apply_filter(self.channel, self._make_groups)

        # convert channels to appropriate quantity if necessary.
        x = self._from_channel(x, group=group)  # knows the units underneath

        if self.units in ('channel',):
            format = '%i'            

        mask = numpy.ones(len(x), dtype=bool)
        if numpy.iterable(self.mask):
            mask = self.mask

        if self.units in ('wavelength',):
            x = x[::-1]
            mask = mask[::-1]
        return create_expr(x, mask, format, delim)

    def get_filter_expr(self):
        return (self.get_filter(format='%.4f', delim='-') +
                ' ' + self.get_xlabel())

    def notice_response(self, notice_resp=True, noticed_chans=None):
        notice_resp=bool_cast(notice_resp)

        if notice_resp and noticed_chans is None:
            noticed_chans = self.get_noticed_channels()

        for id in self.response_ids:
            arf, rmf = self.get_response(id)
            _notice_resp(noticed_chans, arf, rmf)


    def notice(self, lo=None, hi=None, ignore=False, bkg_id=None):
        # If any background IDs are actually given, then impose
        # the filter on those backgrounds *only*, and return.  Do
        # *not* impose filter on data itself.  (Revision possibly
        # this should be done in high-level UI?)  SMD 10/25/12

        filter_background_only = False
        if (bkg_id is not None):
            if (not(numpy.iterable(bkg_id))):
                bkg_id = [bkg_id]
            filter_background_only = True
        else:
            bkg_id = self.background_ids

        # Automatically impose data's filter on background data sets.
        # Units must agree for this to be meaningful, so temporarily
        # make data and background units match. SMD 10/25/12
        for bid in bkg_id:
            bkg = self.get_background(bid)
            old_bkg_units = bkg.units
            bkg.units = self.units
            bkg.notice(lo, hi, ignore)
            bkg.units = old_bkg_units

        # If we're only supposed to filter backgrounds, return
        if (filter_background_only == True):
            return

        # Go on if we are also supposed to filter the source data
        ignore=bool_cast(ignore)
        if lo is None and hi is None:
            self.quality_filter=None
            self.notice_response(False)

        elo,ehi = self._get_ebins()
        if lo is not None and type(lo) != str:
            lo = self._to_channel(lo)
        if hi is not None and type(hi) != str:
            hi = self._to_channel(hi)

        if( ( self.units=="wavelength" and
              elo[0] < elo[-1] and  ehi[0] < ehi[-1] ) or
            ( self.units=="energy" and
              elo[0] > elo[-1] and  ehi[0] > ehi[-1] ) ):
            lo, hi = hi, lo

        # If we are working in channel space, and the data are
        # grouped, we must correct for the fact that bounds expressed
        # expressed in channels must be converted to group number.
        # This is the only set of units for which this must be done;
        # energy and wavelength conversions above already take care of
        # the distinction between grouped and ungrouped.

        if (self.units=="channel" and self.grouped==True):

            if (lo is not None and 
                type(lo) != str and 
                not(lo < self.channel[0])):

                # Find the location of the first channel greater than
                # or equal to lo in self.channel
                # Then find out how many groups there are that contain
                # the channels less than lo, and convert lo from a 
                # channel number to the first group number that has channels 
                # greater than or equal to lo.
                
                lo_index = numpy.where(self.channel >= lo)[0][0]
                lo = len(numpy.where(self.grouping[:lo_index] > -1)[0]) + 1

            if (hi is not None and 
                type(hi) != str and 
                not(hi > self.channel[-1])):

                # Find the location of the first channel greater than
                # or equal to hi in self.channel
                # Then find out how many groups there are that contain
                # the channels less than hi, and convert hi from a 
                # channel number to the first group number that has channels 
                # greater than or equal to hi.
                hi_index = numpy.where(self.channel >= hi)[0][0]
                hi = len(numpy.where(self.grouping[:hi_index] > -1)[0])

                # If the original channel hi starts a new group,
                # increment the group number 
                if (self.grouping[hi_index] > -1):
                    hi = hi + 1

                # If the original channel hi is in a group such that
                # the group has channels greater than original hi,
                # then use the previous group as the highest group included
                # in the filter. Avoid indexing beyond the end of the
                # grouping array.
                if (hi_index + 1 < len(self.grouping)):
                    if (not(self.grouping[hi_index+1] > -1)):
                        hi = hi - 1

        # Don't use the middle of the channel anymore as the
        # grouping function.  That was just plain dumb.
        # So just get back an array of groups 1-N, if grouped
        BaseData.notice(self, (lo,), (hi,),
                        (self.apply_grouping(self.channel,
                                             self._make_groups),),
                        ignore)

    def to_guess(self):
        elo, ehi = self._get_ebins(group=False)
        elo = self.apply_filter(elo, self._min)
        ehi = self.apply_filter(ehi, self._max)
        if self.units=="wavelength":
            lo = self._hc / ehi
            hi = self._hc / elo
            elo = lo; ehi = hi
        cnt = self.get_dep(True)
        arf = self.get_specresp(filter=True)

        y = cnt/(ehi-elo)
        if self.exposure is not None:
            y /= self.exposure           # photons/keV/sec or
                                         # photons/Ang/sec
        #y = cnt/arf/self.exposure
        if arf is not None:
            y /= arf                     # photons/keV/cm^2/sec or 
                                         # photons/Ang/cm^2/sec
        return (y, elo, ehi)

    def to_fit(self, staterrfunc=None):
        return (self.get_dep(True),
                self.get_staterror(True, staterrfunc),
                self.get_syserror(True))

    def to_plot(self, yfunc=None, staterrfunc=None, response_id=None):
        return (self.apply_filter(self.get_x(response_id=response_id), self._middle),
                self.get_y(True, yfunc, response_id=response_id),
                self.get_yerr(True, staterrfunc, response_id=response_id),
                self.get_xerr(True, response_id=response_id),
                self.get_xlabel(),
                self.get_ylabel())

    def group(self):
        "Group the data according to the data set's grouping scheme"
        self.grouped = True

    def ungroup(self):
        "Ungroup the data"
        self.grouped = False

    def subtract(self):
        "Subtract the background data"
        self.subtracted = True

    def unsubtract(self):
        "Remove background subtraction"
        self.subtracted = False

class DataIMG(Data2D):
    "Image data set, including functions for coordinate transformations"

    def _get_coord(self):
        return self._coord
    
    def _set_coord(self, val):
        coord = str(val).strip().lower()

        if coord in ('logical', 'image'):
            coord = 'logical'

        elif coord in ('physical',):
            self._check_physical_transform()
            coord = 'physical'
            
        elif coord in ('world','wcs'):
            self._check_world_transform()
            coord = 'world'

        else:
            raise DataErr('bad', 'coordinates', val)

        self._coord = coord

    coord = property(_get_coord, _set_coord,
                     doc='Coordinate system of independent axes')

    def __init__(self, name, x0, x1, y, shape=None, staterror=None,
                 syserror=None, sky=None, eqpos=None, coord='logical',
                 header=None):
        self._x0 = x0
        self._x1 = x1
        self._region = None
        BaseData.__init__(self)


    def __str__(self):
        # Print the metadata first
        old = self._fields
        ss = old
        try:
            self._fields = filter((lambda x: x!='header'), self._fields)
            ss = BaseData.__str__(self)
        finally:
            self._fields = old
        return ss


    def __getstate__(self):
        state = self.__dict__.copy()
        # Function pointers to methods of the class
        # (of type 'instancemethod') are NOT picklable
        # remove them and restore later with a coord init
        #del state['_get_logical']
        #del state['_get_physical']
        #del state['_get_world']

        # PyRegion objects (of type 'extension') are NOT picklable, yet.
        # preserve the region string and restore later with constructor
        state['_region'] = state['_region'].__str__()
        return state

    def __setstate__(self, state):
        # Populate the function pointers we deleted at pickle time with
        # no-ops.
        #self.__dict__['_get_logical']=(lambda : None)
        #self.__dict__['_get_physical']=(lambda : None)
        #self.__dict__['_get_world']=(lambda : None)

        if not state.has_key('header'):
            self.header=None

        self.__dict__.update(state)

        # _set_coord will correctly define the _get_* WCS function pointers.
        self._set_coord(state['_coord'])
        self._region = Region(self._region)

    def _check_physical_transform(self):
        if self.sky is None:
            raise DataErr('nocoord', self.name, 'physical')

    def _check_world_transform(self):
        if self.eqpos is None:
            raise DataErr('nocoord', self.name, 'world')

    def _logical_to_physical(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_physical_transform()

        # logical -> physical
        x0, x1 = self.sky.apply(x0, x1)

        return (x0, x1)

    def _logical_to_world(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_world_transform()

        # logical -> physical
        if self.sky is not None:
            x0, x1 = self.sky.apply(x0, x1)

        # physical -> world
        x0, x1 = self.eqpos.apply(x0, x1)

        return (x0, x1)

    def _physical_to_logical(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_physical_transform()

        # physical -> logical
        x0, x1 = self.sky.invert(x0, x1)

        return (x0, x1)

    def _physical_to_world(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_world_transform()

        # physical -> world
        x0, x1 = self.eqpos.apply(x0, x1)

        return (x0, x1)

    def _world_to_logical(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_world_transform()

        # world -> physical
        x0, x1 = self.eqpos.invert(x0, x1)

        # physical -> logical
        if self.sky is not None:
            x0, x1 = self.sky.invert(x0, x1)

        return (x0, x1)

    def _world_to_physical(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_world_transform()

        # world -> physical
        x0, x1 = self.eqpos.invert(x0, x1)

        return (x0, x1)

    def get_logical(self):
        coord = self.coord
        x0, x1 = self.get_indep()
        if coord is not 'logical':
            x0 = x0.copy()
            x1 = x1.copy()
            x0, x1 = getattr(self, '_'+coord+'_to_logical')(x0, x1)
        return (x0, x1)

    def get_physical(self):
        coord = self.coord
        x0, x1 = self.get_indep()
        if coord is not 'physical':
            x0 = x0.copy()
            x1 = x1.copy()
            x0, x1 = getattr(self, '_'+coord+'_to_physical')(x0, x1)
        return (x0, x1)


    def get_world(self):
        coord = self.coord
        x0, x1 = self.get_indep()
        if coord is not 'world':
            x0 = x0.copy()
            x1 = x1.copy()
            x0, x1 = getattr(self, '_'+coord+'_to_world')(x0, x1)
        return (x0, x1)


    # For compatibility with old Sherpa keywords
    get_image = get_logical
    get_wcs = get_world

    def set_coord(self, coord):
        coord = str(coord).strip().lower()
        # Destroys original data to conserve memory for big imgs
        good = ('logical','image','physical','world','wcs')
        if coord not in good:
            raise DataErr('badchoices', 'coordinates', coord, ", ".join(good))

        if coord.startswith('wcs'):
            coord = 'world'
        elif coord.startswith('image'):
            coord = 'logical'

        self.x0, self.x1 = getattr(self, 'get_'+coord)()
        self._x0 = self.apply_filter(self.x0)
        self._x1 = self.apply_filter(self.x1)

        self._set_coord(coord)

    def get_filter_expr(self):
        if self._region is not None:
            return str(self._region)
        return ''

    get_filter = get_filter_expr

    def notice2d(self, val=None, ignore=False):
        mask = None
        ignore=bool_cast(ignore)
        if val is not None:
            val = str(val).strip()
            (self._region,
             mask) = region_mask(self._region, val,
                                 self.get_x0(), self.get_x1(),
                                 os.path.isfile(val), ignore)
            mask = numpy.asarray(mask, dtype=numpy.bool_)
        else:
            self._region = None

        if mask is None:
            self.mask = not ignore
            self._region = None

        elif not ignore:
            if self.mask is True:
                self._set_mask(mask)
            else:
                self.mask |= mask
        else:
            mask = ~mask
            if self.mask is False:
                self.mask = mask
            else:
                self.mask &= mask

#        self._x0 = self.apply_filter(self.x0)
#        self._x1 = self.apply_filter(self.x1)

    def get_bounding_mask(self):
        mask = self.mask
        shape = None
        if numpy.iterable(self.mask):
            # create bounding box around noticed image regions
            mask = numpy.array(self.mask).reshape(*self.shape)
            x0_i, x1_i = numpy.where(mask == True)

            x0_lo = x0_i.min()
            x0_hi = x0_i.max()
            x1_lo = x1_i.min()
            x1_hi = x1_i.max()

            shape = mask[x0_lo:x0_hi+1,x1_lo:x1_hi+1].shape
            mask = mask[x0_lo:x0_hi+1,x1_lo:x1_hi+1]

            mask = mask.ravel()
        return mask, shape

    def get_img(self, yfunc=None):
        # FIXME add support for coords to image class -> DS9
        self._check_shape()
        y_img = self.filter_region(self.get_dep(False))
        if yfunc is not None:
            m = self.eval_model_to_fit(yfunc)
            if numpy.iterable(self.mask):
                # if filtered, the calculated model must be padded up
                # to the data size to preserve img shape and WCS coord
                m = pad_bounding_box(m, self.mask)
            y_img = (y_img, self.filter_region(m))

        if yfunc is not None:
            y_img = (y_img[0].reshape(*self.shape),
                     y_img[1].reshape(*self.shape))
        else:
            y_img = y_img.reshape(*self.shape)
	return y_img

    def get_axes(self):
        # FIXME: how to filter an axis when self.mask is size of self.y?
        self._check_shape()

        # dummy placeholders needed b/c img shape may not be square!
        axis0 = numpy.arange(self.shape[1], dtype=float)+1.
        axis1 = numpy.arange(self.shape[0], dtype=float)+1.
        dummy0 = numpy.ones(axis0.size, dtype=float)
        dummy1 = numpy.ones(axis1.size, dtype=float)

        if self.coord == 'physical':
            axis0, dummy = self._logical_to_physical(axis0, dummy0)
            dummy, axis1 = self._logical_to_physical(dummy1, axis1)

        elif self.coord == 'world':
            axis0, dummy = self._logical_to_world(axis0, dummy0)
            dummy, axis1 = self._logical_to_world(dummy1, axis1)

        return (axis0, axis1)

    def get_x0label(self):
        "Return label for first dimension in 2-D view of independent axis/axes"
        if self.coord in ('logical', 'image'):
            return 'x0'
        elif self.coord in ('physical',):
            return 'x0 (pixels)'
        elif self.coord in ('world', 'wcs'):
            return 'RA (deg)'
        else:
            return 'x0'
    
    def get_x1label(self):
        """
        Return label for second dimension in 2-D view of independent axis/axes
        """        
        if self.coord in ('logical', 'image'):
            return 'x1'
        elif self.coord in ('physical',):
            return 'x1 (pixels)'
        elif self.coord in ('world', 'wcs'):
            return 'DEC (deg)'
        else:
            return 'x1'

    def to_contour(self, yfunc=None):
        y = self.filter_region(self.get_dep(False))
        if yfunc is not None:
            m = self.eval_model_to_fit(yfunc)
            if numpy.iterable(self.mask):
                # if filtered, the calculated model must be padded up
                # to the data size to preserve img shape and WCS coord
                m = self.filter_region(pad_bounding_box(m, self.mask))
            y = (y, m)

        return (self.get_x0(),
                self.get_x1(),
                y,
                self.get_x0label(),
                self.get_x1label())

    def filter_region(self, data):
        if data is not None and numpy.iterable(self.mask):
            filter = numpy.ones(len(self.mask), dtype=SherpaFloat)
            filter[~self.mask]=numpy.nan
            return data*filter
        return data

class DataIMGInt(DataIMG):

    def _set_mask(self, val):
        DataND._set_mask(self, val)
        try:
            self._x0lo = self.apply_filter(self.x0lo)
            self._x0hi = self.apply_filter(self.x0hi)
            self._x1lo = self.apply_filter(self.x1lo)
            self._x1hi = self.apply_filter(self.x1hi)
        except DataErr:
            self._x0lo = self.x0lo
            self._x1lo = self.x1lo
            self._x0hi = self.x0hi
            self._x1hi = self.x1hi

    mask = property(DataND._get_mask, _set_mask,
                    doc='Mask array for dependent variable')

    def __init__(self, name, x0lo, x1lo, x0hi, x1hi, y, shape=None, 
                 staterror=None, syserror=None, sky=None, eqpos=None,
                 coord='logical', header=None):
        self._x0lo = x0lo
        self._x1lo = x1lo
        self._x0hi = x0hi
        self._x1hi = x1hi
        self._region = None
        BaseData.__init__(self)

    def set_coord(self, coord):
        coord = str(coord).strip().lower()
        # Destroys original data to conserve memory for big imgs
        good = ('logical','image','physical','world','wcs')
        if coord not in good:
            raise DataErr('bad', 'coordinates', coord)

        if coord.startswith('wcs'):
            coord = 'world'
        elif coord.startswith('image'):
            coord = 'logical'

        self.x0lo, self.x1lo, self.x0hi, self.x1hi = getattr(self, 'get_'+coord)()
        self._x0lo = self.apply_filter(self.x0lo)
        self._x0hi = self.apply_filter(self.x0hi)
        self._x1lo = self.apply_filter(self.x1lo)
        self._x1hi = self.apply_filter(self.x1hi)
        self._set_coord(coord)

    def get_logical(self):
        coord = self.coord
        x0lo, x1lo, x0hi, x1hi = self.get_indep()
        if coord is not 'logical':
            x0lo = x0lo.copy()
            x1lo = x1lo.copy()
            x0lo, x1lo = getattr(self, '_'+coord+'_to_logical')(x0lo, x1lo)

            x0hi = x0hi.copy()
            x1hi = x1hi.copy()
            x0hi, x1hi = getattr(self, '_'+coord+'_to_logical')(x0hi, x1hi)

        return (x0lo, x1lo, x0hi, x1hi)

    def get_physical(self):
        coord = self.coord
        x0lo, x1lo, x0hi, x1hi = self.get_indep()
        if coord is not 'physical':
            x0lo = x0lo.copy()
            x1lo = x1lo.copy()
            x0lo, x1lo = getattr(self, '_'+coord+'_to_physical')(x0lo, x1lo)

            x0hi = x0hi.copy()
            x1hi = x1hi.copy()
            x0hi, x1hi = getattr(self, '_'+coord+'_to_physical')(x0hi, x1hi)

        return (x0lo, x1lo, x0hi, x1hi)

    def get_world(self):
        coord = self.coord
        x0lo, x1lo, x0hi, x1hi = self.get_indep()
        if coord is not 'world':
            x0lo = x0lo.copy()
            x1lo = x1lo.copy()
            x0lo, x1lo = getattr(self, '_'+coord+'_to_world')(x0lo, x1lo)

            x0hi = x0hi.copy()
            x1hi = x1hi.copy()
            x0hi, x1hi = getattr(self, '_'+coord+'_to_world')(x0hi, x1hi)

        return (x0lo, x1lo, x0hi, x1hi)


    # def get_indep(self, filter=False):
    #     x0, x1 = DataIMG.get_indep(self, filter=filter)
        
    #     halfwidth = numpy.array([.5,.5])
    #     if self.coord == 'physical' and self.sky is not None:
    #         halfwidth = numpy.array(self.sky.cdelt)/2.
    #     elif self.coord == 'world' and self.eqpos is not None:
    #         halfwidth = numpy.array(self.eqpos.cdelt)/2.

    #     return (x0-halfwidth[0],x1-halfwidth[1],
    #             x0+halfwidth[0],x1+halfwidth[1])


    def get_indep(self, filter=False):
        filter=bool_cast(filter)
        if filter:
            return (self._x0lo, self._x1lo, self._x0hi, self._x1hi)
	return (self.x0lo, self.x1lo, self.x0hi, self.x1hi)

    def get_x0(self, filter=False):
        indep = self.get_indep(filter)
        return (indep[0] + indep[2]) / 2.0

    def get_x1(self, filter=False):
        indep = self.get_indep(filter)
        return (indep[1] + indep[3]) / 2.0


    def get_axes(self):
        # FIXME: how to filter an axis when self.mask is size of self.y?
        self._check_shape()

        # dummy placeholders needed b/c img shape may not be square!
        axis0lo = numpy.arange(self.shape[1], dtype=float)-0.5
        axis1lo = numpy.arange(self.shape[0], dtype=float)-0.5

        axis0hi = numpy.arange(self.shape[1], dtype=float)+0.5
        axis1hi = numpy.arange(self.shape[0], dtype=float)+0.5

        dummy0 = numpy.ones(axis0lo.size, dtype=float)
        dummy1 = numpy.ones(axis1lo.size, dtype=float)

        if self.coord == 'physical':
            axis0lo, dummy = self._logical_to_physical(axis0lo, dummy0)
            axis0hi, dummy = self._logical_to_physical(axis0hi, dummy0)

            dummy, axis1lo = self._logical_to_physical(dummy1, axis1lo)
            dummy, axis1hi = self._logical_to_physical(dummy1, axis1hi)

        elif self.coord == 'world':
            axis0lo, dummy = self._logical_to_world(axis0lo, dummy0)
            axis0hi, dummy = self._logical_to_world(axis0hi, dummy0)

            dummy, axis1lo = self._logical_to_world(dummy1, axis1lo)
            dummy, axis1hi = self._logical_to_world(dummy1, axis1hi)

        return (axis0lo, axis1lo, axis0hi, axis1hi)
