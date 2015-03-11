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


import numpy
import pyfits
import os
from itertools import izip
from sherpa.utils.err import IOErr
from sherpa.utils import SherpaInt, SherpaUInt, SherpaFloat, is_binary_file
from sherpa.io import get_ascii_data, write_arrays
from sherpa.astro.io.meta import *

import logging
warning = logging.getLogger(__name__).warning
error = logging.getLogger(__name__).error

transformstatus = False
try:
    from sherpa.astro.io.wcs import WCS
    transformstatus = True
except:
    warning('failed to import WCS module; WCS routines will not be ' +
            'available')

__all__ = ('get_table_data', 'get_image_data', 'get_arf_data', 'get_rmf_data',
           'get_pha_data', 'set_table_data', 'set_image_data', 'set_pha_data',
           'get_column_data', 'get_ascii_data')

try:
    # pyfits-1.3 support
    _VLF = pyfits.NP_pyfits._VLF
except AttributeError:
    # pyfits-2.3.1 support
    try:
        _VLF = pyfits.core._VLF
    except AttributeError:
        # pyfits-3.0 support
        _VLF = pyfits.column._VLF


def _has_hdu(hdulist, id):
    try:
        hdulist[id]
    except (KeyError, IndexError):
        return False
    return True


def _has_key(hdu, name):
    #return hdu.header.has_key(name)
    return name in hdu.header


def _try_key(hdu, name, fix_type=False, dtype=SherpaFloat):
    if _has_key(hdu, name):
        key = hdu.header[name]

        if str(key).find('none') != -1:
            return None

        if fix_type:
            key = dtype(key)
        return key

    return None


def _require_key(hdu, name, fix_type=False, dtype=SherpaFloat):
    key = _try_key(hdu, name, fix_type, dtype)
    if key is None:
        raise IOErr('nokeyword', hdu._file.name, name)
    return key


def _get_meta_data(hdu):
    meta = Meta()
    for key in dict(hdu.header.items()).keys():
        val = hdu.header[key]

        # empty numpy strings are not recognized by load pickle!
        if type(val) is numpy.str_ and val == '':
            val = ''

        meta[key] = val
    return meta


def _get_wcs_key(hdu, key0, key1, fix_type=False, dtype=SherpaFloat):
    if _has_key(hdu, key0) and _has_key(hdu, key1):
        return numpy.array([_try_key(hdu, key0, fix_type, dtype),
                            _try_key(hdu, key1, fix_type, dtype)],dtype)
    return ()

def _set_wcs_key(hdulist, name, val):
    card = pyfits.Card( str(name), val )
    hdulist.append(card)

def _try_col(hdu, name, dtype=SherpaFloat, fix_type=False):
    if name not in hdu.columns.names:
        return None

    col = hdu.data.field(name)

    if isinstance(col, _VLF):
        col = numpy.concatenate([numpy.asarray(row) for row in col])
    else:
        col = numpy.asarray(col).ravel()

    if fix_type:
        col = col.astype(dtype)

    return col

def _try_tbl_col(hdu, name, dtype=SherpaFloat, fix_type=False):
    if name not in hdu.columns.names:
        return (None,)

    col = hdu.data.field(name)

    if isinstance(col, _VLF):
        col = numpy.concatenate([numpy.asarray(row) for row in col])
    else:
        col = numpy.asarray(col)

    if fix_type:
        col = col.astype(dtype)

    return numpy.column_stack(col)

def _try_vec(hdu, name, size=2, dtype=SherpaFloat, fix_type=False):
    if name not in hdu.columns.names:
        return numpy.array([None]*size)

    col = hdu.data.field(name)

    if isinstance(col, _VLF):
        col = numpy.concatenate([numpy.asarray(row) for row in col])
    else:
        col = numpy.asarray(col)

    if fix_type:
        col = col.astype(dtype)

    if col is None:
        return numpy.array([None]*size)

    return col

def _require_col(hdu, name, dtype=SherpaFloat, fix_type=False):
    col = _try_col(hdu, name, dtype, fix_type)
    if col is None:
        raise IOErr('reqcol', name, hdu._file.name)
    return col

def _require_tbl_col(hdu, name, dtype=SherpaFloat, fix_type=False):
    col = _try_tbl_col(hdu, name, dtype, fix_type)
    if len(col) > 0 and col[0]==None:
        raise IOErr('reqcol', name, hdu._file.name)
    return col

def _require_vec(hdu, name, size=2, dtype=SherpaFloat, fix_type=False):
    col = _try_vec(hdu, name, size, dtype, fix_type)
    if None in list(col):
        raise IOErr('reqcol', name, hdu._file.name)
    return col

def _try_col_or_key(hdu, name, dtype=SherpaFloat, fix_type=False):
    col = _try_col(hdu, name, dtype, fix_type)
    if col is not None:
        return col
    return _try_key(hdu, name, fix_type, dtype)

def _try_vec_or_key(hdu, name, size, dtype=SherpaFloat, fix_type=False):
    col = _try_col(hdu, name, dtype, fix_type)
    if col is not None:
        return col
    return numpy.array([_try_key(hdu, name, fix_type, dtype)]*size)



## Read Functions ##

def read_table_blocks(arg, make_copy=False):

    filename = ''
    hdus = None
    if type(arg) is pyfits.HDUList:
        filename = arg[0]._file.name
        hdus = arg
    elif type(arg) in (str, unicode, numpy.str_) and is_binary_file(arg):
        filename = arg
        hdus = pyfits.open(arg)
    else:
        raise IOErr('badfile', arg, "a binary FITS table or a PyFITS.BinTableHDU list")

    cols = {}
    hdr = {}

    for ii, hdu in enumerate(hdus):
        blockidx = ii+1

        hdr[blockidx] = {}
        header = hdu.header
        if header is not None:
            for key in header.keys():
                hdr[blockidx][key] = header[key]

        # skip over primary, hdu.data is None

        cols[blockidx] = {}
        recarray = hdu.data
        if recarray is not None:
            for colname in recarray.names:
                cols[blockidx][colname] = recarray[colname]

    return filename, cols, hdr





def get_header_data( arg, blockname=None, hdrkeys=None ):

    filename = ''
    if type(arg) == str and is_binary_file(arg):
        tbl = pyfits.open(arg)
        filename = arg
    elif ( (type(arg) is pyfits.HDUList) and
           (len(arg) > 0) and
           (arg[0].__class__ is pyfits.PrimaryHDU) ):
        tbl = arg
        filename = tbl[0]._file.name
    else:
        raise IOErr('badfile', arg, "a binary FITS table or a PyFITS.BinTableHDU list")

    hdr={}
    try:
        # Use the first binary table extension we find.  Throw an exception
        # if there aren't any.
        for hdu in tbl:
            if blockname is None:
                if hdu.__class__ is pyfits.BinTableHDU:
                    break
                else:
                    continue
            elif (hdu.name.lower() == str(blockname).strip().lower()):
                break

        else:
            raise IOErr('badext', filename)

        if hdrkeys is not None:
            for key in hdrkeys:
                hdr[key] = _require_key(hdu, key, dtype=str)
        else:
            for key in hdu.header.keys():
                hdr[key] = _require_key(hdu, key, dtype=str)

    finally:
        tbl.close()

    return hdr


def get_column_data( *args ):
    """
    get_column_data( *NumPy_args )
    """
    # args is passed as type list
    if len(args) == 0:
        raise IOErr('noarrays')

    cols = []
    for arg in args:
        if not (type(arg) in (numpy.ndarray, list, tuple) or arg is None):
            raise IOErr('badarray', arg)
        if arg is not None:
            vals = numpy.asarray(arg)
            for col in numpy.column_stack(vals):
                cols.append( col )
        else:
            cols.append( arg )

    return cols

def get_table_data(arg, ncols=1, colkeys=None, make_copy=False, fix_type=False,
                   blockname = None, hdrkeys=None):
    """
    get_table_data( filename , ncols=1 [, colkeys=None [, make_copy=False [, blockname=None [, hdrkeys=None ]]]])

    get_table_data( [PrimaryHDU, BinTableHDU] , ncols=1 [, colkeys=None [, make_copy=False [, blockname=None [, hdrkeys=None ]]]])
    """
    filename = ''
    if type(arg) == str and is_binary_file(arg):
        tbl = pyfits.open(arg)
        filename = arg
    elif ( (type(arg) is pyfits.HDUList) and
           (len(arg) > 0) and
           (arg[0].__class__ is pyfits.PrimaryHDU) ):
        tbl = arg
        filename = tbl[0]._file.name
    else:
        raise IOErr('badfile', arg, "a binary FITS table or a PyFITS.BinTableHDU list")

    try:
        # Use the first binary table extension we find.  Throw an exception
        # if there aren't any.
        for hdu in tbl:
            if blockname is None:
                if hdu.__class__ is pyfits.BinTableHDU:
                    break
                else:
                    continue
            elif (hdu.name.lower() == str(blockname).strip().lower() and
                  hdu.__class__ is pyfits.BinTableHDU):
                break

        else:
            raise IOErr('badext', filename)

        cnames = list(hdu.columns.names)

        if colkeys is not None:
            colkeys = [name.strip().upper() for name in list(colkeys)]
        # Try Channel, Counts or X,Y before defaulting to first two table cols
        elif ('CHANNEL' in cnames) and ('COUNTS' in cnames):
            colkeys = ['CHANNEL','COUNTS']
        elif ('X' in cnames) and ('Y' in cnames):
            colkeys = ['X','Y']
        else:
            colkeys = cnames[:ncols]

        cols = []
        for name in colkeys:
            for col in _require_tbl_col(hdu, name, fix_type=fix_type):
                cols.append(col)

        hdr={}
        if hdrkeys is not None:
            for key in hdrkeys:
                hdr[key] = _require_key(hdu, key)

    finally:
        tbl.close()

    return colkeys, cols, filename, hdr


def get_image_data(arg, make_copy=False):
    """
    get_image_data( filename [, make_copy=False ])

    get_image_data( [PrimaryHDU] [, make_copy=False ])
    """
    filename = ''
    if type(arg) == str and is_binary_file(arg):
        hdu = pyfits.open(arg)
        filename = arg
    elif ( (type(arg) is pyfits.HDUList) and
           (len(arg) > 0 ) and
           (arg[0].__class__ is pyfits.PrimaryHDU) ):
        hdu = arg
        filename = hdu[0]._file.name
    else:
        raise IOErr('badfile', arg, "a binary FITS file or a PyFITS.PrimaryHDU list")

#   FITS uses logical-to-world where we use physical-to-world.
#   For all transforms, update their physical-to-world
#   values from their logical-to-world values.
#   Find the matching physical transform
#      (same axis no, but sub = 'P' )
#   and use it for the update.
#   Physical tfms themselves do not get updated.
#
#  Fill the physical-to-world transform given the
#  logical-to-world and the associated logical-to-physical.
#      W = wv + wd * ( P - wp )
#      P = pv + pd * ( L - pp )
#      W = lv + ld * ( L - lp )
# Then
#      L = pp + ( P - pv ) / pd
# so   W = lv + ld * ( pp + (P-pv)/pd - lp )
#        = lv + ( ld / pd ) * ( P - [ pv +  (lp-pp)*pd ] )
# Hence
#      wv = lv
#      wd = ld / pd
#      wp = pv + ( lp - pp ) * pd

#  EG suppose phys-to-world is
#         W =  1000 + 2.0 * ( P - 4.0 )
#  and we bin and scale to generate a logical-to-phys of
#         P =  20 + 4.0 * ( L - 10 )
#  Then
#         W = 1000 + 2.0 * ( (20-4) - 4 * 10 ) + 2 * 4 $
#

    try:
        data = {}

        img = hdu[0]
        if hdu[0].data is None:
            img = hdu[1]
            if hdu[1].data is None:
                raise IOErr('badimg', '')

        data['y'] = numpy.asarray(img.data)

        cdeltp = _get_wcs_key(img, 'CDELT1P', 'CDELT2P')
        crpixp = _get_wcs_key(img, 'CRPIX1P', 'CRPIX2P')
        crvalp = _get_wcs_key(img, 'CRVAL1P', 'CRVAL2P')
        cdeltw = _get_wcs_key(img, 'CDELT1', 'CDELT2')
        crpixw = _get_wcs_key(img, 'CRPIX1', 'CRPIX2')
        crvalw = _get_wcs_key(img, 'CRVAL1', 'CRVAL2')

        # proper calculation of cdelt wrt PHYSICAL coords
        if (( cdeltw != () ) and ( cdeltp != () ) ):
            cdeltw = cdeltw/cdeltp

        # proper calculation of crpix wrt PHYSICAL coords
        if (( crpixw != () ) and ( crvalp != () ) and
            ( cdeltp != () ) and ( crpixp != () ) ):
            crpixw = crvalp + ( crpixw - crpixp ) * cdeltp

        sky = None
        if(cdeltp != () and crpixp != () and crvalp != () and transformstatus):
            sky = WCS('physical', 'LINEAR', crvalp, crpixp, cdeltp)

        eqpos = None
        if(cdeltw != () and crpixw != () and crvalw != () and transformstatus):
            eqpos = WCS('world', 'WCS', crvalw, crpixw, cdeltw)

        data['sky'] = sky
        data['eqpos'] = eqpos
        data['header'] = _get_meta_data(img)

        keys = ['MTYPE1','MFORM1','CTYPE1P','CTYPE2P','WCSNAMEP','CDELT1P',
                'CDELT2P','CRPIX1P','CRPIX2P','CRVAL1P','CRVAL2P',
                'MTYPE2','MFORM2','CTYPE1','CTYPE2','CDELT1','CDELT2','CRPIX1',
                'CRPIX2','CRVAL1','CRVAL2','CUNIT1','CUNIT2','EQUINOX']

        for key in keys:
            try:
                data['header'].pop(key)
            except KeyError:
                pass

    finally:
        hdu.close()

    return data, filename


def get_arf_data(arg, make_copy=False):
    """
    get_arf_data( filename [, make_copy=False ])

    get_arf_data( [PrimaryHDU, BinTableHDU] [, make_copy=False ])
    """
    filename = ''
    if type(arg) == str:
        arf = pyfits.open(arg)
        filename = arg
    elif ( (type(arg) is pyfits.HDUList) and
           (len(arg) > 0) and
           (arg[0].__class__ is pyfits.PrimaryHDU) ):
        arf = arg
        filename = arf[0]._file.name
    else:
        raise IOErr('badfile', arg, "a binary FITS file or a PyFITS.BinTableHDU list")

    try:
        if _has_hdu(arf, 'SPECRESP'):
            hdu = arf['SPECRESP']
        elif _has_hdu(arf, 'AXAF_ARF'):
            hdu = arf['AXAF_ARF']
        elif (_has_hdu(arf, 1) and
              ((_try_key(arf[1], 'HDUCLAS1') == 'SPECRESP') or
               (_try_key(arf[1], 'HDUCLAS2') == 'SPECRESP'))):
            hdu = arf[1]
        else:
            raise IOErr('notrsp', filename, 'an ARF')

        data = {}

        data['exposure'] = _try_key(hdu, 'EXPOSURE', fix_type=True)

        data['energ_lo'] = _require_col(hdu, 'ENERG_LO', fix_type=True)
        data['energ_hi'] = _require_col(hdu, 'ENERG_HI', fix_type=True)
        data['specresp'] = _require_col(hdu, 'SPECRESP', fix_type=True)
        data['bin_lo']   = _try_col(hdu, 'BIN_LO', fix_type=True)
        data['bin_hi']   = _try_col(hdu, 'BIN_HI', fix_type=True)
        data['header']     = _get_meta_data(hdu)
        data['header'].pop('EXPOSURE')

    finally:
        arf.close()

    return data, filename


def get_rmf_data(arg, make_copy=False):
    """
    get_rmf_data( filename [, make_copy=False ])

    get_rmf_data( [PrimaryHDU, BinTableHDU] [, make_copy=False ])
    """
    filename = ''
    if type(arg) == str:
        rmf = pyfits.open(arg)
        filename = arg
    elif ( (type(arg) is pyfits.HDUList) and
           (len(arg) > 0) and
           (arg[0].__class__ is pyfits.PrimaryHDU) ):
        rmf = arg
        filename = rmf[0]._file.name
    else:
        raise IOErr('badfile', arg, "a binary FITS file or a PyFITS.BinTableHDU list")

    try:
        if _has_hdu(rmf, 'MATRIX'):
            hdu = rmf['MATRIX']
        elif _has_hdu(rmf, 'SPECRESP MATRIX'):
            hdu = rmf['SPECRESP MATRIX']
        elif _has_hdu(rmf, 'AXAF_RMF'):
            hdu = rmf['AXAF_RMF']
        elif (_has_hdu(rmf, 1) and
              ((_try_key(rmf[1], 'HDUCLAS1') == 'RESPONSE') or
               (_try_key(rmf[1], 'HDUCLAS2') == 'RSP_MATRIX'))):
            hdu = rmf[1]
        else:
            raise IOErr('notrsp', filename, 'an RMF')

        data = {}

        data['detchans'] = SherpaUInt(_require_key(hdu, 'DETCHANS'))
        data['energ_lo'] = _require_col(hdu, 'ENERG_LO', fix_type=True)
        data['energ_hi'] = _require_col(hdu, 'ENERG_HI', fix_type=True)
        data['n_grp']    = _require_col(hdu, 'N_GRP', fix_type=True,
                                        dtype=SherpaUInt)
        data['f_chan']   = _require_vec(hdu, 'F_CHAN', fix_type=True,
                                        dtype=SherpaUInt)
        data['n_chan']   = _require_vec(hdu, 'N_CHAN', fix_type=True,
                                        dtype=SherpaUInt)
        # Read MATRIX as-is -- we will flatten it below, because
        # we need to remove all rows corresponding to n_grp[row] == 0
        data['matrix']   = None
        if 'MATRIX' not in hdu.columns.names:
            pass
        else:
            data['matrix']   = hdu.data.field('MATRIX')

        data['header']     = _get_meta_data(hdu)
        data['header'].pop('DETCHANS')

        # Beginning of non-Chandra RMF support
        fchan_col = list(hdu.columns.names).index('F_CHAN') + 1
        tlmin = _try_key(hdu, 'TLMIN'+str(fchan_col), True, SherpaUInt)
        if tlmin is not None:
            data['offset'] = tlmin
        else:
            error("Failed to locate TLMIN keyword for F_CHAN" +
                  " column in RMF file '%s'; "  % filename +
                  'Update the offset value in the RMF data set to' +
                  ' appropriate TLMIN value prior to fitting')

        if _has_hdu(rmf, 'EBOUNDS'):
            hdu = rmf['EBOUNDS']
            data['e_min'] = _try_col(hdu, 'E_MIN', fix_type=True)
            data['e_max'] = _try_col(hdu, 'E_MAX', fix_type=True)

            # Beginning of non-Chandra RMF support
            chan_col = list(hdu.columns.names).index('CHANNEL') + 1
            tlmin = _try_key(hdu, 'TLMIN'+str(chan_col), True, SherpaUInt)
            if tlmin is not None:
                data['offset'] = tlmin

        else:
            data['e_min'] = None
            data['e_max'] = None
    finally:
        rmf.close()

    ### For every row i of the response matrix, such that
    ### n_grp[i] == 0, we need to remove that row from the
    ### n_chan, f_chan, and matrix arrays we are constructing
    ### to be passed up to the DataRMF data structure.
 
    ### This is trivial for n_chan and f_chan.  For the matrix
    ### array this can be more work -- can't just remove all
    ### zeroes, because some rows where n_grp[row] > 0 might
    ### still have zeroes in the matrix.  I add new code first
    ### to deal with the matrix, then simpler code to remove zeroes
    ### from n_chan and f_chan.

    # Read in MATRIX column with structure as-is -- i.e., as an array of 
    # arrays.  Then flatten it, but include only those arrays that come from
    # rows where n_grp[row] > 0.  Zero elements can only be included from
    # rows where n_grp[row] > 0.  SMD 05/23/13

    if not isinstance(data['matrix'], _VLF):
        data['matrix'] = None
        raise IOErr('badfile', filename, " MATRIX column not a variable length field")

    good = (data['n_grp'] > 0)
    data['matrix'] = data['matrix'][good]
    data['matrix'] = numpy.concatenate([numpy.asarray(row) for row in data['matrix']])
    data['matrix'] = data['matrix'].astype(SherpaFloat)

    #Flatten f_chan and n_chan vectors into 1D arrays as crates does
    # according to group
    if( (data['f_chan'].ndim > 1) and ( data['n_chan'].ndim > 1) ):
        f_chan = []; n_chan = [];
        for grp,fch,nch, in izip(data['n_grp'],data['f_chan'],data['n_chan']):
            for i in xrange(grp):
                f_chan.append(fch[i])
                n_chan.append(nch[i])

        data['f_chan'] = numpy.asarray(f_chan, SherpaUInt)
        data['n_chan'] = numpy.asarray(n_chan, SherpaUInt)
    else:
        if( len(data['n_grp']) == len(data['f_chan']) ):
            # filter out groups with zeroes.
            good = (data['n_grp'] > 0)
            data['f_chan'] = data['f_chan'][good]
            data['n_chan'] = data['n_chan'][good]
    
    return data, filename


def get_pha_data(arg, make_copy=False, use_background=False):
    """
    get_pha_data( filename [, make_copy=False [, use_background=False[])

    get_pha_data( [PrimaryHDU, BinTableHDU] [, make_copy=False
                  [, use_background=False]])
    """
    filename = ''
    if type(arg) == str and is_binary_file(arg):
        pha = pyfits.open(arg)
        filename = arg
    elif ( (type(arg) is pyfits.HDUList) and
           (len(arg) > 0) and
           (arg[0].__class__ is pyfits.PrimaryHDU) ):
        pha = arg
        filename = pha[0]._file.name
    else:
        raise IOErr('badfile', arg, "a binary FITS spectrum or a PyFITS.BinTableHDU list")

    try:
        if _has_hdu(pha, 'SPECTRUM'):
            hdu = pha['SPECTRUM']
        elif (_has_hdu(pha, 1) and
              ((_try_key(pha[1], 'HDUCLAS1') == 'SPECTRUM') or
               (_try_key(pha[1], 'HDUCLAS2') == 'SPECTRUM'))):
            hdu = pha[1]
        else:
            raise IOErr('notrsp',  filename, "a PHA spectrum")

        if use_background:
            for block in pha:
                if (_try_key(block, 'HDUCLAS2') == 'BKG'):
                         hdu = block

        keys = ['BACKFILE','ANCRFILE','RESPFILE',
                'BACKSCAL','AREASCAL','EXPOSURE']
        datasets = []

        if _try_col(hdu, 'SPEC_NUM') is None:
            data = {}

            # Keywords
            data['exposure'] = _try_key(hdu, 'EXPOSURE', True, SherpaFloat)
            #data['poisserr'] = _try_key(hdu, 'POISSERR', True, bool)
            data['backfile'] = _try_key(hdu, 'BACKFILE')
            data['arffile']  = _try_key(hdu, 'ANCRFILE')
            data['rmffile']  = _try_key(hdu, 'RESPFILE')

            # Keywords or columns
            data['backscal'] = _try_col_or_key(hdu, 'BACKSCAL', fix_type=True)
            data['backscup'] = _try_col_or_key(hdu, 'BACKSCUP', fix_type=True)
            data['backscdn'] = _try_col_or_key(hdu, 'BACKSCDN', fix_type=True)
            data['areascal'] = _try_col_or_key(hdu, 'AREASCAL', fix_type=True)

            # Columns
            data['channel']         = _require_col(hdu, 'CHANNEL', fix_type=True)
            #Make sure channel numbers not indices
            chan = list(hdu.columns.names).index('CHANNEL') + 1
            tlmin = _try_key(hdu, 'TLMIN'+str(chan), True, SherpaUInt)
            if int(data['channel'][0]) == 0 or ((tlmin is not None) and tlmin == 0):
                data['channel'] = data['channel']+1

            data['counts']      = _try_col(hdu, 'COUNTS', fix_type=True)
            if data['counts'] is None:
                data['counts']  = _require_col(hdu, 'RATE', fix_type=True) * data['exposure']
            data['staterror']       = _try_col(hdu, 'STAT_ERR')
            data['syserror']        = _try_col(hdu, 'SYS_ERR')
            data['background_up']   = _try_col(hdu, 'BACKGROUND_UP', fix_type=True)
            data['background_down'] = _try_col(hdu, 'BACKGROUND_DOWN', fix_type=True)
            data['bin_lo']          = _try_col(hdu, 'BIN_LO', fix_type=True)
            data['bin_hi']          = _try_col(hdu, 'BIN_HI', fix_type=True)
            data['grouping']        = _try_col(hdu, 'GROUPING', SherpaInt)
            data['quality']         = _try_col(hdu, 'QUALITY', SherpaInt)
            data['header']            = _get_meta_data(hdu)
            for key in keys:
                try:
                    data['header'].pop(key)
                except KeyError:
                    pass

            if data['syserror'] is not None:
                # SYS_ERR is the fractional systematic error
                data['syserror'] = data['syserror'] * data['counts']

            datasets.append(data)

        else:
            data = {}
            # Type 2 PHA file support

            specnum = _try_col_or_key(hdu, 'SPEC_NUM')
            num = len(specnum)

            # Keywords
            exposure = _try_key(hdu, 'EXPOSURE', True, SherpaFloat)
            #poisserr = _try_key(hdu, 'POISSERR', True, bool)
            backfile = _try_key(hdu, 'BACKFILE')
            arffile  = _try_key(hdu, 'ANCRFILE')
            rmffile  = _try_key(hdu, 'RESPFILE')

            # Keywords or columns
            backscal = _try_vec_or_key(hdu, 'BACKSCAL', num, fix_type=True)
            backscup = _try_vec_or_key(hdu, 'BACKSCUP', num, fix_type=True)
            backscdn = _try_vec_or_key(hdu, 'BACKSCDN', num, fix_type=True)
            areascal = _try_vec_or_key(hdu, 'AREASCAL', num, fix_type=True)

            # Columns
            channel         = _require_vec(hdu, 'CHANNEL', num, fix_type=True)

            #Make sure channel numbers not indices
            chan = list(hdu.columns.names).index('CHANNEL') + 1
            tlmin = _try_key(hdu, 'TLMIN'+str(chan), True, SherpaUInt)
            
            for ii in range(num):
                if int(channel[ii][0]) == 0:
                    channel[ii] += 1
            #if ((tlmin is not None) and tlmin == 0) or int(channel[0]) == 0:
            #    channel += 1

            counts =  _try_vec(hdu, 'COUNTS', num, fix_type=True)
            if None in counts:
                counts =  _require_vec(hdu, 'RATE', num, fix_type=True) * data['exposure']
            staterror       = _try_vec(hdu, 'STAT_ERR', num)
            syserror        = _try_vec(hdu, 'SYS_ERR', num)
            background_up   = _try_vec(hdu, 'BACKGROUND_UP', num, fix_type=True)
            background_down = _try_vec(hdu, 'BACKGROUND_DOWN', num, fix_type=True)
            bin_lo          = _try_vec(hdu, 'BIN_LO', num, fix_type=True)
            bin_hi          = _try_vec(hdu, 'BIN_HI', num, fix_type=True)
            grouping        = _try_vec(hdu, 'GROUPING', num, SherpaInt)
            quality         = _try_vec(hdu, 'QUALITY', num, SherpaInt)

            orders          = _try_vec(hdu, 'TG_M', num, SherpaInt)
            parts           = _try_vec(hdu, 'TG_PART', num, SherpaInt)
            specnums        = _try_vec(hdu, 'SPEC_NUM', num, SherpaInt)
            srcids          = _try_vec(hdu, 'TG_SRCID', num, SherpaInt)


            # Iterate over all rows of channels, counts, errors, etc
            # Populate a list of dictionaries containing
            # individual dataset info
            for (bscal, bscup, bscdn, arsc, chan, cnt, staterr, syserr,
                 backup, backdown, binlo, binhi, group, qual, ordr, prt,
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
                data['grouping']        = group
                data['quality']         = qual
                data['header']            = _get_meta_data(hdu)
                data['header']['TG_M']     = ordr
                data['header']['TG_PART']  = prt
                data['header']['SPEC_NUM'] = specnum
                data['header']['TG_SRCID'] = srcid

                for key in keys:
                    try:
                        data['header'].pop(key)
                    except KeyError:
                        pass

                if syserr is not None:
                    # SYS_ERR is the fractional systematic error
                    data['syserror'] = syserr * cnt

                datasets.append(data)

    finally:
        pha.close()

    return datasets, filename


#
### Write Functions ###
#


def set_table_data(filename, data, col_names, hdr=None, hdrnames=None,
                   ascii=False, clobber=False, packup=False):

    if not packup and os.path.isfile(filename) and not clobber:
        raise IOErr("filefound", filename)

    col_names = list(col_names)
    col_names.remove("name")
    #hdrlist = pyfits.CardList()

    #for name in ['exposure','backscal', 'areascal']:
    #    hdrlist.append(pyfits.Card(key=name.upper(),
    #                   value=data[name]))

    collist = []
    cols = []
    coldefs = []
    for name in col_names:
        if data[name] is None:
            continue
        col = pyfits.Column(name=name.upper(),
                            format=data[name].dtype.name.upper(),
                            array=data[name])
        cols.append(data[name])
        coldefs.append(name.upper())
        collist.append(col)

    if ascii:
        set_arrays(filename, cols, coldefs, ascii=ascii, clobber=clobber)
        return

    tbl = pyfits.new_table(pyfits.ColDefs(collist))
                           #, header=pyfits.Header(hdrlist))
    tbl.name = 'HISTOGRAM'
    if packup:
        return tbl
    tbl.writeto(filename, clobber=True)


def set_pha_data(filename, data, col_names, header=None,
                 ascii=False, clobber=False, packup=False):

    if not packup and os.path.isfile(filename) and not clobber:
        raise IOErr("filefound", filename)

    hdrlist = pyfits.CardList()

    for key in header.keys():
        if header[key] is None:
            continue
        hdrlist.append(pyfits.Card( str(key.upper()), header[key] ))

    collist = []
    cols = []
    coldefs = []
    for name in col_names:
        if data[name] is None:
            continue
        col = pyfits.Column(name=name.upper(),
                            format=data[name].dtype.name.upper(),
                            array=data[name])
        cols.append(data[name])
        coldefs.append(name.upper())
        collist.append(col)

    if ascii:
        set_arrays(filename, cols, coldefs, ascii=ascii, clobber=clobber)
        return

    pha = pyfits.new_table(pyfits.ColDefs(collist),
                           header=pyfits.Header(hdrlist))
    pha.name = 'SPECTRUM'
    if packup:
        return pha
    pha.writeto(filename, clobber=True)


def set_image_data(filename, data, header, ascii=False, clobber=False,
                   packup=False):

    if not packup and os.path.isfile(filename) and not clobber:
        raise IOErr("filefound", filename)

    if ascii:
        set_arrays(filename, [data['pixels'].ravel()],
                   ascii=ascii, clobber=clobber)
        return

    hdrlist = pyfits.CardList()

    # Write Image Header Keys
    for key in header.keys():
        if header[key] is None:
            continue
        _set_wcs_key(hdrlist, key, header[key])

    # Write Image WCS Header Keys
    if data['eqpos'] is not None:
        cdeltw = data['eqpos'].cdelt
        crpixw = data['eqpos'].crpix
        crvalw = data['eqpos'].crval
        equin  = data['eqpos'].equinox

    if data['sky'] is not None:
        cdeltp = data['sky'].cdelt
        crpixp = data['sky'].crpix
        crvalp = data['sky'].crval

        _set_wcs_key(hdrlist, 'MTYPE1', 'sky     ')
        _set_wcs_key(hdrlist, 'MFORM1', 'x,y     ')
        _set_wcs_key(hdrlist, 'CTYPE1P', 'x      ')
        _set_wcs_key(hdrlist, 'CTYPE2P', 'y      ')
        _set_wcs_key(hdrlist, 'WCSNAMEP','PHYSICAL')
        _set_wcs_key(hdrlist, 'CDELT1P', cdeltp[0])
        _set_wcs_key(hdrlist, 'CDELT2P', cdeltp[1])
        _set_wcs_key(hdrlist, 'CRPIX1P', crpixp[0])
        _set_wcs_key(hdrlist, 'CRPIX2P', crpixp[1])
        _set_wcs_key(hdrlist, 'CRVAL1P', crvalp[0])
        _set_wcs_key(hdrlist, 'CRVAL2P', crvalp[1])

        if data['eqpos'] is not None:
            # Simply the inverse of read transformations in get_image_data
            cdeltw = cdeltw * cdeltp
            crpixw = ((crpixw - crvalp) / cdeltp + crpixp )

    if data['eqpos'] is not None:
        _set_wcs_key(hdrlist, 'MTYPE2', 'EQPOS   ')
        _set_wcs_key(hdrlist, 'MFORM2', 'RA,DEC  ')
        _set_wcs_key(hdrlist, 'CTYPE1', 'RA---TAN')
        _set_wcs_key(hdrlist, 'CTYPE2', 'DEC--TAN')
        _set_wcs_key(hdrlist, 'CDELT1', cdeltw[0])
        _set_wcs_key(hdrlist, 'CDELT2', cdeltw[1])
        _set_wcs_key(hdrlist, 'CRPIX1', crpixw[0])
        _set_wcs_key(hdrlist, 'CRPIX2', crpixw[1])
        _set_wcs_key(hdrlist, 'CRVAL1', crvalw[0])
        _set_wcs_key(hdrlist, 'CRVAL2', crvalw[1])
        _set_wcs_key(hdrlist, 'CUNIT1', 'deg     ')
        _set_wcs_key(hdrlist, 'CUNIT2', 'deg     ')
        _set_wcs_key(hdrlist, 'EQUINOX', equin)

    #
    img = pyfits.PrimaryHDU(data['pixels'], header=pyfits.Header(hdrlist))
    if packup:
        return img
    img.writeto(filename, clobber=True)


def set_arrays(filename, args, fields=None, ascii=True, clobber=False):

    if ascii:
        write_arrays(filename, args, fields, clobber=clobber)
        return

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

    if fields is None:
        fields = ['col%i' % (ii+1) for ii in range(len(args))]

    if len(args) != len(fields):
        raise IOErr("toomanycols", str(len(fields)), str(len(args)))

    cols = []
    for val, name in izip(args, fields):
        col = pyfits.Column(name=name.upper(), format=val.dtype.name.upper(),
                            array=val)
        cols.append(col)

    tbl = pyfits.new_table(pyfits.ColDefs(cols))
    tbl.name = 'TABLE'
    tbl.writeto(filename, clobber=True)
