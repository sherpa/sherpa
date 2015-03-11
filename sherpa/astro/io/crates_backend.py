# 
#  Copyright (C) 2011  Smithsonian Astrophysical Observatory
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


from itertools import izip
import os.path
import numpy
from sherpa.utils.err import IOErr
from sherpa.utils import SherpaInt, SherpaUInt, SherpaFloat, is_binary_file
from sherpa.astro.utils import resp_init
from sherpa.astro.io.meta import *
import pycrates

import logging
warning = logging.getLogger(__name__).warning
error = logging.getLogger(__name__).error
info    = logging.getLogger(__name__).info

transformstatus = False
try:
    from sherpa.astro.io.wcs import WCS
    transformstatus = True
except:
    warning('failed to import WCS module; WCS routines will not be ' +
            'available')

__all__ = ('get_table_data', 'get_image_data', 'get_arf_data', 'get_rmf_data',
           'get_pha_data','set_table_data', 'set_image_data', 'set_pha_data',
           'get_column_data', 'get_ascii_data')



def open_crate_dataset(filename, crateType=pycrates.CrateDataset, mode='r'):
    """
    checked for new crates
    """
    dataset = crateType(filename, mode=mode)
    return dataset


def open_crate(filename, crateType=pycrates.CrateDataset, mode='r'):
    """
    checked for new crates
    """
    dataset = open_crate_dataset(filename, crateType, mode)
    current = dataset.get_current_crate()
    data = dataset.get_crate(current)
    return data

def close_crate_dataset(dataset):

    if (hasattr(dataset, "snip")):
        dataset.snip()
    # Else, it's a plain Crate that can just fall out of scope

def get_filename_from_dmsyntax(filename, blockname=None):

    arg = str(filename)
    isbinary = True
    colnames = True
    dmsyn = ''

    if '[' in filename and ']' in filename:
        parts = filename.split('[')
        filename = parts.pop(0)
        if parts:
            dmsyn = parts.pop(0).lower()

    if not is_binary_file(filename):
        isbinary = False
        fd = open(filename, 'r')
        try:
            last=None
            line = fd.readline().strip()
            while len(line) > 0 and line[0] in '#%':
                last = line
                line = fd.readline().strip()
            if (last is not None and
                (len(last.split(' ')) != len(line.split(' ')) )):
                colnames = False
        finally:
            fd.close()

    if blockname is not None:
        arg += "[%s]" % str(blockname).upper()

    if (not isbinary) and (not colnames) and (not 'cols' in dmsyn):
        arg += "[opt colnames=none]"

    return arg


def _try_hdr_key(crate, name, dtype=str):
    """
    checked for new crates
    """
    name = name.upper().strip()
    if not crate.key_exists(name):
        return None

    return _try_key(crate, name, dtype)


def _try_key(crate, name, dtype=str):
    """
    checked for new crates
    """

    if not crate.key_exists(name):
        return None

    key = pycrates.get_key(crate, name)

    if key is None:
        return None

    key = key.value

    if str(key).find('none') != -1:
        return None

    return dtype( key )


def _get_meta_data(crate):
    """
    checked for new crates
    """
    meta = Meta()
    names = pycrates.get_key_names(crate)
    if names is not None and numpy.iterable(names):
        for name in names:
            val = pycrates.get_keyval(crate, name)

            # empty numpy strings are not recognized by load pickle!
            if type(val) is numpy.str_ and val == '':
                val = ''

            meta[name] = val
    return meta


def _set_key(crate, name, val, fix_type=False, dtype=str):
    """
    checked for new crates
    """
    key = pycrates.CrateKey()
    key.name = str(name).strip().upper()

    if fix_type:
        val = dtype(val)

    key.value = val
    pycrates.add_key(crate, key)


def _set_column(crate, name, val):
    """
    checked for new crates
    """
    col = pycrates.CrateData()
    col.name = name.upper()
    col.values = numpy.array(val)
    pycrates.add_col(crate, col)


def _set_pha_column(crate, col, name, val):
    col.name = name.upper()

    if col.load( numpy.asarray(val), True) != pycrates.dmSUCCESS:
        raise IOErr('setcolfailed', col.name, col.get_status_message())

    exec('crate.' + name + ' = col')


def _require_key(crate, name, key, dtype=str):
    """
    checked for new crates
    """
    key = _try_key(crate, key, dtype)
    if key is None:
        raise IOErr('nokeyword', crate.get_filename(), name)
    return key


def _require_hdr_key(crate, name, dtype=str):
    """
    checked for new crates
    """
    key = _try_hdr_key(crate, name, dtype)
    if key is None:
        raise IOErr('nokeyword', crate.get_filename(), name)
    return key


def _try_unit(col, dtype=str):
    """
    checked for new crates
    """
    if col is None:
        return None

    unit = col.unit

    if str(unit) == '':
        return None

    return str(unit)


def _try_col(crate, colname, make_copy=False, fix_type=False, dtype=SherpaFloat):
    """
    checked for new crates
    """
    if not crate.column_exists(colname):
        return None

    col = crate.get_column(colname)

    if col is None:
        return None

    if make_copy:
        # Make a copy if a filename passed in
        data = numpy.array(col.values).ravel()

    else:
        # Use a reference if a crate passed in
        data = numpy.asarray(col.values).ravel()

    if fix_type:
        data = data.astype(dtype)

    return data


def _require_tbl_col(crate, colname, cnames, make_copy=False,
                     fix_type=False, dtype=SherpaFloat):
    """
    checked for new crates
    """
    name = str(colname).strip()
    if not crate.column_exists(name):
        raise IOErr('reqcol', name, cnames)

    #data = pycrates.get_colvals(crate, name)
    col = crate.get_column(name)

    if make_copy:
        # Make a copy if a filename passed in
        data = numpy.array(col.values)
    else:
        # Use a reference if a crate passed in
        data = numpy.asarray(col.values)

    if fix_type:
        data = data.astype(dtype)

    return numpy.column_stack(data)


def _try_key_list(crate, keyname, num, dtype=SherpaFloat, fix_type=False):
    """
    checked for new crates
    """
    if not crate.key_exists(keyname):
        return None

    key = crate.get_key(keyname)

    # Make a copy of the data, since we don't know that pycrates will
    # do something sensible wrt reference counting
    key = numpy.array([key.value]*num)

    if fix_type:
        key = key.astype(dtype)

    return key


def _try_col_list(crate, colname, num, make_copy=False, fix_type=False,
                  dtype=SherpaFloat):
    """
    checked for the new crates
    """
    if not crate.column_exists(colname):
        return numpy.array([None]*num)

    col = crate.get_column(colname)

    if make_copy:
        # Make a copy if a filename passed in
        col = numpy.array(col.values)
    else:
        # Use a reference if a crate passed in
        col = numpy.asarray(col.values)

    if fix_type:
        col = col.astype(dtype)

    return col


def _require_col_list(crate, colname, num, make_copy=False, fix_type=False,
                      dtype=SherpaFloat):
    """
    check for new crates
    """
    if not crate.column_exists(colname):
        raise IOErr('badcol', colname)
    return _try_col_list(crate, colname, num, make_copy, fix_type, dtype)


def _require_col(crate, colname, make_copy=False, fix_type=False, dtype=SherpaFloat):
    """
    checked for new crates
    """
    data = _try_col(crate, colname, make_copy, fix_type, dtype)
    if data is None:
        raise IOErr('badcol', colname)
    return data


def _try_image(crate, make_copy=False, fix_type=False, dtype=SherpaFloat):
    """
    checked for new crates
    """

    if make_copy:
        # Make a copy if a filename passed in
        dat = pycrates.copy_piximgvals(crate).squeeze()
    else:
        # Use a reference if a crate passed in
        dat = pycrates.get_piximgvals(crate).squeeze()

    if fix_type:
        dat = dat.astype(dtype)

    # FITS standard is FORTRAN filling, Sherpa is c-style
#    return dat.reshape(dat.shape[::-1])

    # Crates now returns image in c-style filling
    return dat


def _require_image(crate, make_copy=False, fix_type=False, dtype=SherpaFloat):
    """
    checked for new crates
    """
    dat = _try_image(crate, make_copy, fix_type, dtype)
    if dat is None:
        raise IOErr('badimg', crate.get_filename())
    return dat


def _get_crate_by_blockname(dataset, blkname):
    ncrates = dataset.get_ncrates()
    for ii in range(1, ncrates+1):
        crate = dataset.get_crate(ii)
        if crate.name.strip().lower() == blkname.strip().lower():
            return crate
    return None


## Read Functions ##

def read_table_blocks(arg, make_copy=False):

    filename = ''
    dataset = None
    close_dataset = False
    if isinstance(arg, pycrates.TABLECrate):
        filename = arg.get_filename()
        dataset = arg.get_dataset()
    elif isinstance(arg, pycrates.CrateDataset):
        filename = arg.get_filename()
        dataset = arg
    elif type(arg) in (str, unicode, numpy.str_):
        filename = arg
        dataset = pycrates.CrateDataset(arg)
        close_dataset = True
    else:
        raise IOErr('badfile', arg, "CrateDataset obj")

    # index of block number starts at 1
    blockidx = numpy.arange(dataset.get_ncrates()) + 1

    cols = {}
    hdr = {}

    for ii in blockidx:
        crate = dataset.get_crate(ii)
        hdr[ii] = {}
        names = crate.get_keynames()
        for name in names:
            hdr[ii][name] = crate.get_key_value(name)

        cols[ii] = {}
        # skip over primary
        if crate.name == 'PRIMARY':
            continue

        names = crate.get_colnames()
        for name in names:
            cols[ii][name] = crate.get_column(name).values

    if close_dataset:
        close_crate_dataset(dataset)
    return filename, cols, hdr


def get_header_data( arg, blockname=None, hdrkeys=None ):
    """
    checked for new crates
    """

    filename = ''
    close_dataset = False
    if type(arg) == str:

        #arg = get_filename_from_dmsyntax(arg, blockname)
        arg = get_filename_from_dmsyntax(arg)

        tbl = None
        try:
            tbl = open_crate(arg)
        except Exception, e:
            raise e

        close_dataset = True
        filename = tbl.get_filename()

        # Make a copy of the data, since we don't know that pycrates will
        # do something sensible wrt reference counting
    elif isinstance(arg, pycrates.TABLECrate):
        tbl = arg
        filename = arg.get_filename()
    else:
        raise IOErr('badfile', arg, 'TABLECrate obj')

    # Crates "caches" open files by their filename in memory.  If you try
    # to open a file multiple times (with DM syntax) it corrupts the Crate
    # in memory.  This is a work-around to open the CrateDataset without
    # DM syntax and iterate through the crates looking for the block
    # name that matches.
    if blockname is not None:
        crate = _get_crate_by_blockname(tbl.get_dataset(), blockname)
        tbl = crate or tbl

    hdr = {}
    if hdrkeys is None:
        #hdrkeys = tbl.get_keynames()
        hdrkeys = pycrates.get_key_names(tbl)

    for key in hdrkeys:
        hdr[key] = _require_hdr_key(tbl, key)

    if close_dataset:
        close_crate_dataset(tbl.get_dataset())
    return hdr


def get_column_data( *args ):
    """
    checked for new crates

    get_column_data( *NumPy_args )

    get_column_data( *CrateData_args )

    get_column_data( *NumPy_and_CrateData_args )
    """
    # args is passed as type list
    if len(args) == 0:
        raise IOErr('noarrays')

    cols = []
    for arg in args:
        if isinstance(arg, pycrates.CrateData):
            #vals = arg.get_values()
            vals = arg.values

        elif arg is None or type(arg) in (numpy.ndarray, list, tuple):
            vals = arg
        else:
            raise IOErr('badarray', arg)

        if arg is not None:
            vals = numpy.asarray( vals )
            for col in numpy.column_stack(vals):
                cols.append( col )
        else:
            cols.append( vals )

    return cols

def get_ascii_data(filename, ncols=2, colkeys=None, **kwargs):
    """
    get_table_data( filename [, ncols=2 [, colkeys=None [, **kwargs ]]] )
    """
    return get_table_data( filename, ncols, colkeys )[:3]


def get_table_data( arg, ncols=1, colkeys=None, make_copy=True, fix_type=True,
                    blockname=None, hdrkeys=None):
    """
    get_table_data( filename , ncols=1 [, colkeys=None [, make_copy=True [,
                    fix_type=True [, blockname=None [, hdrkeys=None ]]]]])

    get_table_data( TABLECrate , ncols=1 [, colkeys=None [, make_copy=True [,
                    fix_type=True [, blockname=None [, hdrkeys=None ] ]]]])
    """
    filename = ''
    close_dataset = False
    if type(arg) == str:

        arg = get_filename_from_dmsyntax(arg)
        tbl = open_crate(arg)
        if not isinstance(tbl, pycrates.TABLECrate):
            #######??????????????????????????????????######## dtn
            close_crate_dataset( tbl.get_dataset() )
            #######??????????????????????????????????######## dtn
            raise IOErr('badfile', arg, 'TABLECrate obj')

        filename = tbl.get_filename()
        close_dataset = True

        # Make a copy of the data, since we don't know that pycrates will
        # do something sensible wrt reference counting
    elif isinstance(arg, pycrates.TABLECrate):
        tbl = arg
        filename = arg.get_filename()
        make_copy=False
    else:
        raise IOErr('badfile', arg, 'TABLECrate obj')


    # Crates "caches" open files by their filename in memory.  If you try
    # to open a file multiple times (with DM syntax) it corrupts the Crate
    # in memory.  This is a work-around to open the CrateDataset without
    # DM syntax and iterate through the crates looking for the block
    # name that matches.
    if blockname is not None:
        crate = _get_crate_by_blockname(tbl.get_dataset(), blockname)
        tbl = crate or tbl


    cnames = list(pycrates.get_col_names(tbl, vectors=False, rawonly=True))

    if colkeys is not None:
        colkeys = [str(name).strip() for name in list(colkeys)]

    elif (type(arg) == str and (not os.path.isfile(arg))
          and '[' in arg and ']' in arg):
        colkeys = cnames

    # Try Channel, Counts or X,Y before defaulting to first two table cols
    elif 'CHANNEL' in cnames and 'COUNTS' in cnames:
        colkeys = ['CHANNEL','COUNTS']

    elif 'X' in cnames and 'Y' in cnames:
        colkeys = ['X','Y']

    else:
        colkeys = cnames[:ncols]

    cols = []
    for name in colkeys:
        for col in _require_tbl_col(tbl, name, cnames, make_copy, fix_type):
            cols.append(col)

    hdr={}
    if hdrkeys is not None:
        for key in hdrkeys:
            hdr[key] = _require_hdr_key(tbl, key)

    if close_dataset:
        close_crate_dataset(tbl.get_dataset())
    return colkeys, cols, filename, hdr


def get_image_data(arg, make_copy=True, fix_type=True):
    """
    get_image_data ( filename [, make_copy=True, fix_type=True ])

    get_image_data ( IMAGECrate [, make_copy=True, fix_type=True ])
    """
    filename = ''
    close_dataset = False
    if type(arg) == str:
        img = open_crate(arg)

        if not isinstance(img, pycrates.IMAGECrate):
            #######??????????????????????????????????######## dtn
            close_crate_dataset( img.get_dataset() )
            #######??????????????????????????????????######## dtn
            raise IOErr('badfile', arg, "IMAGECrate obj")

        filename = arg
        close_dataset = True

    elif isinstance(arg, pycrates.IMAGECrate):
        img = arg
        filename = arg.get_filename()
        make_copy=False

    else:
        raise IOErr('badfile', arg, "IMAGECrate obj")

    data = {}

    data['y'] = _require_image(img, make_copy, fix_type)

    sky = None
    skynames = ['SKY', 'sky', 'pos', 'POS']
    names = img.get_axisnames()

    # find the SKY name using the set intersection
    inter = list(set(names) & set(skynames))
    if inter:
        sky = img.get_transform(inter[0])

    wcs = None
    if 'EQPOS' in names:
        wcs = img.get_transform('EQPOS')

    if sky is not None and transformstatus:
        linear = pycrates.WCSTANTransform()
        linear.set_name("LINEAR")
        linear.set_transform_matrix(sky.get_transform_matrix())
        cdelt = numpy.array(linear.get_parameter_value('CDELT'))
        crpix = numpy.array(linear.get_parameter_value('CRPIX'))
        crval = numpy.array(linear.get_parameter_value('CRVAL'))
        data['sky'] = WCS('physical', 'LINEAR', crval, crpix, cdelt)

    if wcs is not None and transformstatus:
        cdelt = numpy.array(wcs.get_parameter_value('CDELT'))
        crpix = numpy.array(wcs.get_parameter_value('CRPIX'))
        crval = numpy.array(wcs.get_parameter_value('CRVAL'))
        crota = SherpaFloat(wcs.get_parameter_value('CROTA'))
        equin = SherpaFloat(wcs.get_parameter_value('EQUINOX'))
        epoch = SherpaFloat(wcs.get_parameter_value('EPOCH'))
        data['eqpos'] = WCS('world', 'WCS', crval, crpix, cdelt,
                            crota, epoch, equin)

    data['header'] = _get_meta_data(img)

    keys = ['MTYPE1','MFORM1','CTYPE1P','CTYPE2P','WCSNAMEP','CDELT1P',
            'CDELT2P','CRPIX1P','CRPIX2P','CRVAL1P','CRVAL2P',
            'MTYPE2','MFORM2','CTYPE1','CTYPE2','CDELT1','CDELT2','CRPIX1',
            'CRPIX2','CRVAL1','CRVAL2','CUNIT1','CUNIT2','EQUINOX']
#            'WCSTY1P', 'WCSTY2P']

    for key in keys:
        try:
            data['header'].pop(key)
        except KeyError:
            pass

    if close_dataset:
        close_crate_dataset(img.get_dataset())
    return data, filename


def get_arf_data(arg, make_copy=True):
    """
    get_arf_data( filename [, make_copy=True ])

    get_arf_data( ARFCrate [, make_copy=True ])
    """
    filename = ''
    close_dataset = False
    if type(arg) == str:
        arf = open_crate(arg)
        if not isinstance(arf, pycrates.TABLECrate):
            #######??????????????????????????????????######## dtn
            close_crate_dataset( arf.get_dataset() )
            #######??????????????????????????????????######## dtn
            raise IOErr('badfile', arg, "ARFCrate obj")
        filename = arg
        close_dataset = True
    elif isinstance(arg, pycrates.TABLECrate):
        arf = arg
        filename = arg.get_filename()
        make_copy=False

    else:
        raise IOErr('badfile', arg, "ARFCrate obj")

    if arf is None or arf.get_colnames() is None:
        raise IOErr('filenotfound', arg)

    data = {}

    if not arf.column_exists('ENERG_LO'):
        raise IOErr('reqcol', 'ENERG_LO', filename)
    if not arf.column_exists('ENERG_HI'):
        raise IOErr('reqcol', 'ENERG_HI', filename)
    if not arf.column_exists('SPECRESP'):
        raise IOErr('reqcol', 'SPECRESP', filename)

    data['energ_lo'] = _require_col(arf, 'ENERG_LO', make_copy, fix_type=True)
    data['energ_hi'] = _require_col(arf, 'ENERG_HI', make_copy, fix_type=True)
    data['specresp'] = _require_col(arf, 'SPECRESP', make_copy, fix_type=True)
    data['bin_lo']   = _try_col(arf, 'BIN_LO', make_copy, fix_type=True)
    data['bin_hi']   = _try_col(arf, 'BIN_HI', make_copy, fix_type=True)
    data['exposure'] = _try_key(arf, 'EXPOSURE', dtype=SherpaFloat)
    data['header']   = _get_meta_data(arf)
    data['header'].pop('EXPOSURE')

    if close_dataset:
        close_crate_dataset(arf.get_dataset())
    return data, filename


def get_rmf_data(arg, make_copy=True):
    """
    get_rmf_data( filename [, make_copy=True ])

    get_rmf_data( RMFCrate [, make_copy=True ])
    """
    filename = ''
    close_dataset = False
    if type(arg) == str:
        rmfdataset = open_crate_dataset(arg, pycrates.rmfcratedataset.RMFCrateDataset)

        #if (isinstance(rmfdataset, (pycrates.TABLECrate, pycrates.IMAGECrate)) or
        #    pycrates.is_pha(rmfdataset) == 1):
        #    raise IOErr('badfile', arg, "RMFCrateDataset obj")

        if pycrates.is_rmf(rmfdataset) != 1:
            raise IOErr('badfile', arg, "RMFCrateDataset obj")

        filename = arg
        close_dataset = True

    elif pycrates.is_rmf(arg) == 1:
        rmfdataset = arg
        filename = arg.get_filename()
        make_copy = False

    else:
        raise IOErr('badfile', arg, "RMFCrateDataset obj")

    # Open the response matrix by extension name, and try using
    # some of the many, many ways people break the OGIP definition
    # of the extension name for the response matrix.
    rmf = _get_crate_by_blockname(rmfdataset, 'MATRIX')

    if rmf is None:
        rmf = _get_crate_by_blockname(rmfdataset, 'SPECRESP MATRIX')

    if rmf is None:
        rmf = _get_crate_by_blockname(rmfdataset, 'AXAF_RMF')

    if rmf is None:
        rmf = _get_crate_by_blockname(rmfdataset, 'RSP_MATRIX')

    if rmf is None:
        try:
            rmf = rmfdataset.get_crate(2)
        except IndexError:
            rmf = None

    if rmf is None or rmf.get_colnames() is None:
        raise IOErr('filenotfound', arg)

    data = {}

    if not rmf.column_exists('ENERG_LO'):
        raise IOErr('reqcol', 'ENERG_LO', filename)

    if not rmf.column_exists('ENERG_HI'):
        raise IOErr('reqcol', 'ENERG_HI', filename)

    # FIXME: this will be a problem now that we have
    # to pass the name of the matrix column

    if not rmf.column_exists('MATRIX'):
        raise IOErr('reqcol', 'MATRIX', filename)

    if not rmf.column_exists('N_GRP'):
        raise IOErr('reqcol', 'N_GRP', filename)

    if not rmf.column_exists('F_CHAN'):
        raise IOErr('reqcol', 'F_CHAN', filename)

    if not rmf.column_exists('N_CHAN'):
        raise IOErr('reqcol', 'N_CHAN', filename)


    data['detchans'] = _require_hdr_key(rmf, 'DETCHANS', SherpaInt)
    data['energ_lo'] = _require_col(rmf, 'ENERG_LO', make_copy, fix_type=True)
    data['energ_hi'] = _require_col(rmf, 'ENERG_HI', make_copy, fix_type=True)
    data['n_grp']    = _require_col(rmf, 'N_GRP', make_copy, dtype=SherpaUInt, fix_type=True)

    f_chan = rmf.get_column('F_CHAN')
    offset = f_chan.get_tlmin()

    fcbuf = _require_col(rmf, 'F_CHAN', make_copy)
    ncbuf = _require_col(rmf, 'N_CHAN', make_copy)

    respbuf = _require_col_list(rmf, 'MATRIX', 1, make_copy)

    #ebounds = None
    #if rmfdataset.get_current_crate() < rmfdataset.get_ncrates():
    #    ebounds = rmfdataset.get_crate(rmfdataset.get_current_crate() + 1)
    ebounds = _get_crate_by_blockname(rmfdataset, 'EBOUNDS')

    if ebounds is None:
        ebounds = rmfdataset.get_crate(3)

    data['header'] = _get_meta_data(rmf)
    data['header'].pop('DETCHANS')

    channel = None
    if ebounds is not None:
        data['e_min'] = _try_col(ebounds, 'E_MIN', make_copy, fix_type=True)
        data['e_max'] = _try_col(ebounds, 'E_MAX', make_copy, fix_type=True)
        if ebounds.column_exists('CHANNEL'):
            channel = ebounds.get_column('CHANNEL')

        # FIXME: do I include the header keywords from ebounds
        # data['header'].update(_get_meta_data(ebounds))


    if offset < 0:
        error("Failed to locate TLMIN keyword for F_CHAN" +
              " column in RMF file '%s'; "  % filename +
              'Update the offset value in the RMF data set to' +
              ' appropriate TLMIN value prior to fitting')

    if offset < 0 and channel is not None:
        offset = channel.get_tlmin()

    # If response is non-OGIP, tlmin is -(max of type), so resort to default
    if not (offset < 0):
        data['offset'] = offset

    #
    # FIXME:
    #
    # Currently, CRATES does something screwy:  If n_grp is zero in a bin,
    # it appends a zero to f_chan, n_chan, and matrix.  I have no idea what
    # the logic behind this is -- why would you add data that you know you
    # don't need?  Although it's easy enough to filter the zeros out of
    # f_chan and n_chan, it's harder for matrix, since zero is a legitimate
    # value there.
    #
    # I think this crazy behavior of CRATES should be changed, but for the
    # moment we'll just punt in this case.  (If we don't, the calculation
    # in rmf_fold() will be trashed.)


    # CRATES does not support variable length arrays, so here we condense
    # the array of tuples into the proper length array

    chan_width = data['n_grp'].max()
    resp_width = 0
    if len(respbuf.shape) > 1:
        resp_width = respbuf.shape[1]

    (data['f_chan'], data['n_chan'],
     data['matrix'] ) = resp_init( data['n_grp'], fcbuf, ncbuf,
                                   chan_width, respbuf.ravel(), resp_width )

    if close_dataset:
        close_crate_dataset(rmfdataset)
    return data, filename


def get_pha_data(arg, make_copy=True, use_background=False):
    """
    get_pha_data( filename [, make_copy=True [, use_background=False]])

    get_pha_data( PHACrate [, make_copy=True [, use_background=False]])
    """
    filename = ''
    close_dataset = False
    if type(arg) == str:
        phadataset = open_crate_dataset(arg, pycrates.phacratedataset.PHACrateDataset)

        if pycrates.is_pha(phadataset) != 1:
            raise IOErr('badfile', arg, "PHACrateDataset obj")

        filename = arg
        close_dataset = True

    elif pycrates.is_pha(arg) == 1:
        phadataset = arg
        filename = arg.get_filename()
        make_copy=False

    else:
        raise IOErr('badfile', arg, "PHACrateDataset obj")

    pha = _get_crate_by_blockname(phadataset, "SPECTRUM")

    if pha is None:
        pha = phadataset.get_crate(phadataset.get_current_crate())
        if (pha.get_key('HDUCLAS1').value == 'SPECTRUM' or
            pha.get_key('HDUCLAS2').value == 'SPECTRUM'):
            pass
        else:
            pha = phadataset.get_crate(1)
            if (pha.get_key('HDUCLAS1').value == 'SPECTRUM' or
                pha.get_key('HDUCLAS2').value == 'SPECTRUM'):
                pass
            else:
                # If background maybe better to go on to next block?
                pha = None

    if use_background:

        # Used to read BKGs found in an additional block of
        # Chandra Level 3 PHA files
        for ii in range(phadataset.get_ncrates()):
            block = phadataset.get_crate(ii+1)
            hduclas2 = block.get_key('HDUCLAS2')
            if hduclas2 is not None and hduclas2.value == 'BKG':
                pha = block


    if pha is None or pha.get_colnames() is None:
        raise IOErr('filenotfound', arg)

    keys = ['BACKFILE','ANCRFILE','RESPFILE',
            'BACKSCAL','AREASCAL','EXPOSURE']

    keys_or_cols = ['BACKSCAL','BACKSCUP','BACKSCDN','AREASCAL']

    datasets = []

    # Calling phadataset.is_pha_type1() is unreliable when
    # both TYPE:I and TYPE:II keywords are in the header.
    # Here, I instead test for a column, SPEC_NUM, that can
    # *only* be present in Type II. SMD 05/15/13
    if _try_col(pha, 'SPEC_NUM') is None:
        data = {}

        # Keywords
        data['exposure'] = _try_key(pha, 'EXPOSURE', SherpaFloat)
        #data['poisserr'] = _try_key(pha, 'POISSERR', bool)
        data['backfile'] = _try_key(pha, 'BACKFILE')
        data['arffile']  = _try_key(pha, 'ANCRFILE')
        data['rmffile']  = _try_key(pha, 'RESPFILE')

        # Keywords or columns
        for name in keys_or_cols:
            key = name.lower()
            data[key] = _try_key(pha, name, SherpaFloat)
            if data[key] is None:
                data[key] = _try_col(pha, name, make_copy)

        data['header'] = _get_meta_data(pha)
        for key in keys:
            try:
                data['header'].pop(key)
            except KeyError:
                pass

        # Columns

        if not pha.column_exists('CHANNEL'):
            raise IOErr('reqcol', 'CHANNEL', filename)

        data['channel'] = _require_col(pha, 'CHANNEL', make_copy, fix_type=True)
        # Make sure channel numbers, not indices
        if int(data['channel'][0]) == 0 or pha.get_column('CHANNEL').get_tlmin() == 0:
            data['channel'] = data['channel']+1

        data['counts'] = None
        if pha.column_exists('COUNTS'):
            data['counts'] = _require_col(pha, 'COUNTS', make_copy, fix_type=True)
        else:
            if not pha.column_exists('RATE'):
                raise IOErr('reqcol', 'COUNTS or RATE', filename)
            data['counts'] = _require_col(pha, 'RATE', make_copy, fix_type=True)*data['exposure']

        data['staterror']       = _try_col(pha, 'STAT_ERR', make_copy)
        data['syserror']        = _try_col(pha, 'SYS_ERR', make_copy)
        data['background_up']   = _try_col(pha, 'BACKGROUND_UP', make_copy, fix_type=True)
        data['background_down'] = _try_col(pha, 'BACKGROUND_DOWN', make_copy, fix_type=True)
        data['bin_lo']          = _try_col(pha, 'BIN_LO', make_copy,fix_type=True)
        data['bin_hi']          = _try_col(pha, 'BIN_HI', make_copy,fix_type=True)
        data['grouping']        = _try_col(pha, 'GROUPING', make_copy)
        data['quality']         = _try_col(pha, 'QUALITY', make_copy)

        datasets.append(data)

    else:
        # Type 2 PHA file support
        data = {}
        num = pha.get_nrows()

        # Keywords
        exposure = _try_key(pha, 'EXPOSURE', SherpaFloat)
        #poisserr = _try_key(pha, 'POISSERR', bool)
        backfile = _try_key(pha, 'BACKFILE')
        arffile  = _try_key(pha, 'ANCRFILE')
        rmffile  = _try_key(pha, 'RESPFILE')

        # Keywords or columns
        backscal = _try_key_list(pha, 'BACKSCAL', num)
        if backscal is None:
            backscal = _try_col_list(pha, 'BACKSCAL', num, make_copy)

        backscup = _try_key_list(pha, 'BACKSCUP', num)
        if backscup is None:
            backscup = _try_col_list(pha, 'BACKSCUP', num, make_copy)

        backscdn = _try_key_list(pha, 'BACKSCDN', num)
        if backscdn is None:
            backscdn = _try_col_list(pha, 'BACKSCDN', num, make_copy)

        areascal = _try_key_list(pha, 'AREASCAL', num)
        if areascal is None:
            areascal = _try_col_list(pha, 'AREASCAL', num, make_copy)

        # Columns

        if not pha.column_exists('CHANNEL'):
            raise IOErr('reqcol', 'CHANNEL', filename)

        channel         = _require_col_list(pha, 'CHANNEL', num, make_copy, fix_type=True)
        # Make sure channel numbers, not indices
        for ii in range(num):
            if int(channel[ii][0]) == 0:
                channel[ii] += 1

        counts = None
        if pha.column_exists('COUNTS'):
            counts      = _require_col_list(pha, 'COUNTS', num, make_copy, fix_type=True)
        else:
            if not pha.column_exists('RATE'):
                raise IOErr('reqcol', 'COUNTS or RATE', filename)
            counts      = _require_col_list(pha, 'RATE', num, make_copy, fix_type=True) * exposure

        staterror       = _try_col_list(pha, 'STAT_ERR', num, make_copy)
        syserror        = _try_col_list(pha, 'SYS_ERR', num, make_copy)
        background_up   = _try_col_list(pha, 'BACKGROUND_UP', num, make_copy, fix_type=True)
        background_down = _try_col_list(pha, 'BACKGROUND_DOWN', num, make_copy, fix_type=True)
        bin_lo          = _try_col_list(pha, 'BIN_LO', num, make_copy, fix_type=True)
        bin_hi          = _try_col_list(pha, 'BIN_HI', num, make_copy, fix_type=True)
        grouping        = _try_col_list(pha, 'GROUPING', num, make_copy)
        quality         = _try_col_list(pha, 'QUALITY', num, make_copy)

        orders = _try_key_list(pha, 'TG_M', num)
        if orders is None:
            orders = _try_col_list(pha, 'TG_M', num, make_copy)

        parts = _try_key_list(pha, 'TG_PART', num)
        if parts is None:
            parts = _try_col_list(pha, 'TG_PART', num, make_copy)

        specnums = _try_col_list(pha, 'SPEC_NUM', num, make_copy)
        srcids   = _try_col_list(pha, 'TG_SRCID', num, make_copy)

        # Iterate over all rows of channels, counts, errors, etc
        # Populate a list of dictionaries containing individual dataset info
        for (bscal, bscup, bscdn, arsc, chan, cnt, staterr, syserr,
             backup, backdown, binlo, binhi, grp, qual, ordr, prt,
             specnum, srcid
             ) in izip(backscal, backscup, backscdn, areascal, channel,
                       counts, staterror, syserror, background_up,
                       background_down, bin_lo, bin_hi, grouping, quality,
                       orders, parts, specnums, srcids):

            data = {}

            data['exposure'] = exposure
            #data['poisserr'] = poisserr
            data['backfile'] = backfile
            data['arffile']  = arffile
            data['rmffile']  = rmffile

            data['backscal'] = bscal
            data['backscup'] = bscup
            data['backscdn'] = bscdn
            data['areascal'] = arsc

            data['channel']         = chan
            data['counts']          = cnt
            data['staterror']       = staterr
            data['syserror']        = syserr
            data['background_up']   = backup
            data['background_down'] = backdown
            data['bin_lo']          = binlo
            data['bin_hi']          = binhi
            data['grouping']        = grp
            data['quality']         = qual
            data['header']            = _get_meta_data(pha)
            data['header']['TG_M']     = ordr
            data['header']['TG_PART']  = prt
            data['header']['SPEC_NUM'] = specnum
            data['header']['TG_SRCID'] = srcid

            for key in keys:
                try:
                    data['header'].pop(key)
                except KeyError:
                    pass

            datasets.append(data)

    if close_dataset:
        close_crate_dataset(phadataset)
    return datasets, filename


#
### Write/Pack Functions ####
#


def set_image_data(filename, data, header, ascii=False, clobber=False,
                   packup=False):

    if not packup and os.path.isfile(filename) and not clobber:
        raise IOErr('filefound', filename)

    img = pycrates.IMAGECrate()

    # Write Image Header Keys
    for key in header.keys():
        if header[key] is None:
            continue
        _set_key(img, key, header[key])

    # Write Image WCS Header Keys
    if data['eqpos'] is not None:
        cdeltw = data['eqpos'].cdelt
        crvalw = data['eqpos'].crval
        crpixw = data['eqpos'].crpix
        equin  = data['eqpos'].equinox

    if data['sky'] is not None:
        cdeltp = data['sky'].cdelt
        crvalp = data['sky'].crval
        crpixp = data['sky'].crpix

        _set_key(img, 'MTYPE1', 'sky     ')
        _set_key(img, 'MFORM1', 'x,y     ')
        _set_key(img, 'CTYPE1P','x       ')
        _set_key(img, 'CTYPE2P','y       ')
        _set_key(img, 'WCSNAMEP','PHYSICAL')
        _set_key(img, 'CDELT1P', cdeltp[0])
        _set_key(img, 'CDELT2P', cdeltp[1])
        _set_key(img, 'CRPIX1P', crpixp[0])
        _set_key(img, 'CRPIX2P', crpixp[1])
        _set_key(img, 'CRVAL1P', crvalp[0])
        _set_key(img, 'CRVAL2P', crvalp[1])

        if data['eqpos'] is not None:
            # Simply the inverse of read transformations in get_image_data
            cdeltw = cdeltw * cdeltp
            crpixw = ((crpixw - crvalp) /  cdeltp + crpixp )

    if data['eqpos'] is not None:
        _set_key(img, 'MTYPE2', 'EQPOS   ')
        _set_key(img, 'MFORM2', 'RA,DEC  ')
        _set_key(img, 'CTYPE1', 'RA---TAN')
        _set_key(img, 'CTYPE2', 'DEC--TAN')
        _set_key(img, 'CDELT1', cdeltw[0])
        _set_key(img, 'CDELT2', cdeltw[1])
        _set_key(img, 'CRPIX1', crpixw[0])
        _set_key(img, 'CRPIX2', crpixw[1])
        _set_key(img, 'CRVAL1', crvalw[0])
        _set_key(img, 'CRVAL2', crvalw[1])
        _set_key(img, 'EQUINOX', equin)

    # Write Image pixel values
    pix_col = pycrates.CrateData()
    pix_col.values = data['pixels']
    #pycrates.add_piximg(img, pix_col)
    img.add_image(pix_col)

    if packup:
        return img

    if ascii and '[' not in filename and ']' not in filename:
        #filename += "[opt kernel=text/simple]"
        raise IOErr('writenoimg')

    #pycrates.write_file(img, filename)
    img.write(filename, clobber=True)
    close_crate_dataset(img.get_dataset())


def set_table_data(filename, data, col_names, hdr=None, hdrnames=None,
                   ascii=False, clobber=False, packup=False):

    if not packup and os.path.isfile(filename) and not clobber:
        raise IOErr("filefound", filename)

    tbl = pycrates.TABLECrate()

    col_names = [name for name in col_names if data[name] != None]
    col_names.remove('name')
    try:
        for name in col_names:
            if data[name] is None:
                continue
            _set_column(tbl, name, data[name])

    finally:
        if packup:
            return tbl

        if ascii and '[' not in filename and ']' not in filename:
            filename += "[opt kernel=text/simple]"

        #pycrates.write_file(tbl, filename)
        tbl.write(filename, clobber=True)
        close_crate_dataset(tbl.get_dataset())


def set_pha_data(filename, data, col_names, header=None,
                 ascii=False, clobber=False, packup=False):

    if not packup and os.path.isfile(filename) and not clobber:
        raise IOErr("filefound", filename)

    phadataset = pycrates.phacratedataset.PHACrateDataset()

    # FIXME: Placeholder for pycrates2 bug
    phadataset.set_rw_mode('rw')


    pha = pycrates.TABLECrate()
    pha.name = "SPECTRUM"
    try:

        # Write header values using CrateKey objects
        for key in header.keys():
            if header[key] is None:
                continue
            _set_key(pha, key, header[key])

        # Write column values using CrateData objects
        for name in col_names:
            if data[name] is None:
                continue
            _set_column(pha, name, data[name])

    finally:
        phadataset.add_crate(pha)

        if packup:
            return phadataset

        if ascii and '[' not in filename and ']' not in filename:
            filename += "[opt kernel=text/simple]"

        #pycrates.write_pha(phadataset, filename)
        phadataset.write(filename, clobber=True)
        close_crate_dataset(phadataset)

def set_arrays(filename, args, fields=None, ascii=True, clobber=False):

    if os.path.isfile(filename) and not clobber:
        raise IOErr("filefound", filename)

    if not numpy.iterable(args) or len(args) == 0:
        raise IOErr('noarrayswrite')

    if not numpy.iterable(args[0]):
        raise IOErr('noarrayswrite')

    size = len(args[0])
    for arg in args:
        if not numpy.iterable(arg):
            raise IOErr('noarrayswrite')
        elif len(arg) != size:
            raise IOErr('arraysnoteq')

    if ascii and '[' not in filename and ']' not in filename:
        filename += "[opt kernel=text/simple]"

    tbl = pycrates.TABLECrate()

    if fields is None:
        fields = ['col%i' % (ii+1) for ii in range(len(args))]

    if len(args) != len(fields):
        raise IOErr('toomanycols', str(len(fields)), str(len(args)))

    for val, name in izip(args, fields):
        _set_column(tbl, name, val)

    pycrates.write_file(tbl, filename, clobber=True)
    close_crate_dataset(tbl.get_dataset())
