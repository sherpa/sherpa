#
#  Copyright (C) 2007, 2015, 2016  Smithsonian Astrophysical Observatory
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
Provide I/O routines to Sherpa, focussing on FITS.

This module should not be imported directly if there is no
FITS support present (i.e. neither the crates or astropy
back ends are in use).
"""

from six.moves import zip as izip
from six.moves.configparser import ConfigParser

import logging
import os
import os.path
import sys
import numpy
import sherpa.io
from sherpa.utils.err import IOErr
from sherpa.utils import SherpaFloat
from sherpa.data import Data2D, Data1D, BaseData, Data2DInt
from sherpa.astro.data import DataIMG, DataIMGInt, DataARF, DataRMF, DataPHA
from sherpa import get_config

import importlib

config = ConfigParser()
config.read(get_config())
io_opt = config.get('options', 'io_pkg')
io_opt = str(io_opt).strip().lower()

if io_opt.startswith('pycrates') or io_opt.startswith('crates'):
    io_opt = 'crates_backend'

elif io_opt.startswith('pyfits'):
    io_opt = 'pyfits_backend'

try:
    importlib.import_module('.' + io_opt, package='sherpa.astro.io')
    backend = sys.modules['sherpa.astro.io.' + io_opt]
except ImportError:
    raise ImportError("""Cannot import selected FITS I/O backend {}.
    If you are using CIAO, this is most likely an error and you should contact the CIAO helpdesk.
    If you are using Standalone Sherpa, please install astropy (pyfits should also work, but it's deprecated)."""
                      .format(io_opt))

warning = logging.getLogger(__name__).warning
info    = logging.getLogger(__name__).info


read_table_blocks = backend.read_table_blocks

__all__ = ('read_table', 'read_image', 'read_arf', 'read_rmf', 'read_arrays',
           'read_pha', 'write_image', 'write_pha', 'write_table',
           'pack_table', 'pack_image', 'pack_pha', 'read_table_blocks')


def _is_subclass(t1, t2):
    return isinstance(t1, type) and issubclass(t1, t2) and (t1 is not t2)


def read_arrays(*args):
    """
    read_table( *NumPy_args [, dstype = Data1D ] )

    read_table( *CrateData_args [, dstype = Data1D ] )

    read_table( *NumPy_and_CrateData_args [, dstype = Data1D ] )
    """
    args = list(args)
    if len(args) == 0:
        raise IOErr('noarrays')

    dstype = Data1D
    if _is_subclass(args[-1], BaseData):
        dstype = args.pop()

    args = backend.get_column_data(*args)

    # Determine max number of args for dataset constructor
    sherpa.io._check_args(len(args), dstype)

    return dstype('', *args)


def read_table(arg, ncols=2, colkeys=None, dstype=Data1D):
    """
    read_table( filename [, ncols=2 [, colkeys=None [, dstype=Data1D ]]] )

    read_table( TABLECrate [, ncols=2 [, colkeys=None [, dstype=Data1D ]]] )
    """
    colnames, cols, name, hdr = backend.get_table_data(arg, ncols, colkeys)

    # Determine max number of args for dataset constructor
    sherpa.io._check_args(len(cols), dstype)

    return dstype(name, *cols)


def read_ascii(filename, ncols=2, colkeys=None, dstype=Data1D, **kwargs):
    """
    read_ascii( filename [, ncols=2 [, colkeys=None [, sep=' ' [,
                dstype=Data1D [, comment='#' ]]]]] )
    """
    colnames, cols, name = backend.get_ascii_data(filename, ncols=ncols,
                                                  colkeys=colkeys,
                                                  dstype=dstype, **kwargs)

    # Determine max number of args for dataset constructor
    sherpa.io._check_args(len(cols), dstype)

    return dstype(name, *cols)


def read_image(arg, coord='logical', dstype=DataIMG):
    """
    read_image( filename [, coord='logical' [, dstype=DataIMG ]])

    read_image( IMAGECrate [, coord='logical' [, dstype=DataIMG ]])
    """
    data, filename = backend.get_image_data(arg)
    axlens = data['y'].shape

    x0 = numpy.arange(axlens[1], dtype=SherpaFloat) + 1.
    x1 = numpy.arange(axlens[0], dtype=SherpaFloat) + 1.
    x0, x1 = numpy.meshgrid(x0, x1)
    x0 = x0.ravel()
    x1 = x1.ravel()

    data['y'] = data['y'].ravel()
    data['coord'] = coord
    data['shape'] = axlens

    if issubclass(dstype, DataIMGInt):
        dataset = dstype(filename, x0 - 0.5, x1 - 0.5, x0 + 0.5, x1 + 0.5,
                         **data)
    elif issubclass(dstype, Data2DInt):
        for name in ['coord', 'eqpos', 'sky', 'header']:
            data.pop(name, None)
        dataset = dstype(filename, x0 - 0.5, x1 - 0.5, x0 + 0.5, x1 + 0.5,
                         **data)
    else:
        dataset = dstype(filename, x0, x1, **data)

    return dataset


def read_arf(arg):
    """
    read_arf( filename )

    read_arf( ARFCrate )
    """
    data, filename = backend.get_arf_data(arg)
    return DataARF(filename, **data)


def read_rmf(arg):
    """
    read_rmf( filename )

    read_rmf( RMFCrate )
    """
    data, filename = backend.get_rmf_data(arg)
    return DataRMF(filename, **data)


def _read_ancillary(data, key, label, dname,
                    read_func, output_once=True):
    """Read in the file, if the key in data is set,
    replacing the value with the full path to the file,
    and return the resulf of calling read_func on the
    file. dname is the directory where the input file is,
    and is used to create the full path if it is not an
    absolute path.
    """

    if not(data[key]) or data[key].lower() == 'none':
        return None

    out = None
    try:
        if os.path.dirname(data[key]) == '':
            data[key] = os.path.join(dname, data[key])

        out = read_func(data[key])
        if output_once:
            info('read {} file {}'.format(label, data[key]))

    except:
        if output_once:
            warning(str(sys.exc_info()[1]))

    return out


def read_pha(arg, use_errors=False, use_background=False):
    """
    read_pha( filename [, use_errors=False [, use_background=False]] )

    read_pha( PHACrate [, use_errors=False [, use_background=False]] )
    """
    datasets, filename = backend.get_pha_data(arg,
                                              use_background=use_background)
    phasets = []
    output_once = True
    for data in datasets:
        if not use_errors:
            if data['staterror'] is not None or data['syserror'] is not None:
                if data['staterror'] is None:
                    msg = 'systematic'
                elif data['syserror'] is None:
                    msg = 'statistical'
                    if output_once:
                        wmsg = "systematic errors were not found in " + \
                               "file '{}'".format(filename)
                        warning(wmsg)
                else:
                    msg = 'statistical and systematic'
                if output_once:
                    imsg = msg + " errors were found in file " + \
                           "'{}' \nbut not used; ".format(filename) + \
                           "to use them, re-read with use_errors=True"
                    info(imsg)
                data['staterror'] = None
                data['syserror'] = None

        dname = os.path.dirname(filename)
        albl = 'ARF'
        rlbl = 'RMF'
        if use_background:
            albl = albl + ' (background)'
            rlbl = rlbl + ' (background)'

        arf = _read_ancillary(data, 'arffile', albl, dname, read_arf,
                              output_once)
        rmf = _read_ancillary(data, 'rmffile', rlbl, dname, read_rmf,
                              output_once)

        backgrounds = []

        if data['backfile'] and data['backfile'].lower() != 'none':
            try:
                if os.path.dirname(data['backfile']) == '':
                    data['backfile'] = os.path.join(os.path.dirname(filename),
                                                    data['backfile'])

                bkg_datasets = []
                # Do not read backgrounds of backgrounds
                if not use_background:
                    bkg_datasets = read_pha(data['backfile'], use_errors, True)

                    if output_once:
                        info('read background file {}'.format(
                            data['backfile']))

                if numpy.iterable(bkg_datasets):
                    for bkg_dataset in bkg_datasets:
                        if bkg_dataset.get_response() == (None, None) and \
                           rmf is not None:
                            bkg_dataset.set_response(arf, rmf)
                        backgrounds.append(bkg_dataset)
                else:
                    if bkg_datasets.get_response() == (None, None) and \
                       rmf is not None:
                        bkg_datasets.set_response(arf, rmf)
                    backgrounds.append(bkg_datasets)

            except:
                if output_once:
                    warning(str(sys.exc_info()[1]))

        for bkg_type, bscal_type in izip(('background_up', 'background_down'),
                                         ('backscup', 'backscdn')):
            if data[bkg_type] is not None:
                b = DataPHA(filename,
                            channel=data['channel'],
                            counts=data[bkg_type],
                            bin_lo=data['bin_lo'],
                            bin_hi=data['bin_hi'],
                            grouping=data['grouping'],
                            quality=data['quality'],
                            exposure=data['exposure'],
                            backscal=data[bscal_type],
                            header=data['header'])
                b.set_response(arf, rmf)
                if output_once:
                    info("read {} into a dataset from file {}".format(
                        bkg_type, filename))
                backgrounds.append(b)

        for k in ['backfile', 'arffile', 'rmffile', 'backscup', 'backscdn',
                  'background_up', 'background_down']:
            data.pop(k, None)

        pha = DataPHA(filename, **data)
        pha.set_response(arf, rmf)
        for id, b in enumerate(backgrounds):
            if b.grouping is None:
                b.grouping = pha.grouping
                b.grouped = (b.grouping is not None)
            if b.quality is None:
                b.quality = pha.quality
            pha.set_background(b, id + 1)

        # set units *after* bkgs have been set
        pha._set_initial_quantity()
        phasets.append(pha)
        output_once = False

    if len(phasets) == 1:
        phasets = phasets[0]

    return phasets


def _pack_table(dataset):
    names = dataset._fields
    cols = [(name, getattr(dataset, name)) for name in names]
    data = dict(cols)
    return data, names


def _pack_image(dataset):
    if not(isinstance(dataset, Data2D) or
           isinstance(dataset, DataIMG)):
        raise IOErr('notimage', dataset.name)
    data = {}

    header = {}
    if hasattr(dataset, 'header') and type(dataset.header) is dict:
        header = dataset.header.copy()

    data['pixels'] = numpy.asarray(dataset.get_img())
    data['sky'] = getattr(dataset, 'sky', None)
    data['eqpos'] = getattr(dataset, 'eqpos', None)

    return data, header


def _set_keyword(header, label, value):
    """Extract the string form of the value, and store
    the last path element in the header using the label
    key.
    """

    if value is None:
        return

    name = getattr(value, 'name', 'none')
    if name is not None and name.find('/') != -1:
        name = name.split('/')[-1]

    header[label] = name


def _pack_pha(dataset):
    if not isinstance(dataset, DataPHA):
        raise IOErr('notpha', dataset.name)

    data = {}

    arf, rmf = dataset.get_response()
    bkg = dataset.get_background()

    # Header Keys
    header = {}
    if hasattr(dataset, 'header'):  # and type(dataset.header) is dict:
        header = dataset.header.copy()

    header['EXPOSURE'] = getattr(dataset, 'exposure', 'none')

    _set_keyword(header, 'RESPFILE', rmf)
    _set_keyword(header, 'ANCRFILE', arf)
    _set_keyword(header, 'BACKFILE', bkg)

    # Columns
    col_names = ['channel', 'counts', 'stat_err', 'sys_err',
                 'bin_lo', 'bin_hi', 'grouping', 'quality']

    data['channel'] = getattr(dataset, 'channel', None)
    data['counts'] = getattr(dataset, 'counts', None)
    data['stat_err'] = getattr(dataset, 'staterror', None)
    data['sys_err'] = getattr(dataset, 'syserror', None)
    data['bin_lo'] = getattr(dataset, 'bin_lo', None)
    data['bin_hi'] = getattr(dataset, 'bin_hi', None)
    data['grouping'] = getattr(dataset, 'grouping', None)
    data['quality'] = getattr(dataset, 'quality', None)

    backscal = getattr(dataset, 'backscal', None)
    if backscal is not None:
        if numpy.isscalar(backscal):
            header['BACKSCAL'] = backscal
        else:
            data['backscal'] = backscal
            col_names.append('backscal')

    areascal = getattr(dataset, 'areascal', None)
    if areascal is not None:
        if numpy.isscalar(areascal):
            header['AREASCAL'] = areascal
        else:
            data['areascal'] = areascal
            col_names.append('areascal')

    return data, col_names, header


def write_arrays(filename, args, fields=None, ascii=True, clobber=False):
    backend.set_arrays(filename, args, fields, ascii=ascii, clobber=clobber)


def write_table(filename, dataset, ascii=True, clobber=False):
    data, names = _pack_table(dataset)
    backend.set_table_data(filename, data, names, ascii=ascii, clobber=clobber)


def write_image(filename, dataset, ascii=True, clobber=False):
    data, hdr = _pack_image(dataset)
    backend.set_image_data(filename, data, hdr, ascii=ascii, clobber=clobber )


def write_pha(filename, dataset, ascii=True, clobber=False):
    data, col_names, hdr = _pack_pha(dataset)
    backend.set_pha_data(filename, data, col_names, hdr, ascii=ascii,
                         clobber=clobber)


def pack_table(dataset):
    data, names = _pack_table(dataset)
    return backend.set_table_data('', data, names, packup=True)


def pack_image(dataset):
    data, hdr = _pack_image(dataset)
    return backend.set_image_data('', data, hdr, packup=True)


def pack_pha(dataset):
    data, col_names, hdr = _pack_pha(dataset)
    return backend.set_pha_data('', data, col_names, hdr, packup=True)
