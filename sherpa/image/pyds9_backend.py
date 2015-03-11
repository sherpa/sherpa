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

from itertools import izip
import numpy
import time
import ds9

from sherpa.utils import get_keyword_defaults, SherpaFloat
from sherpa.utils.err import DS9Err

_target = 'sherpa'

def _get_win():
    return ds9.ds9(_target)    

def doOpen():
    _get_win()
    
def isOpen():

    targets = ds9.ds9_targets()
    if targets is None:
        return False

    if type(targets) in (list,):
        for target in targets:
            if _target in target:
                return True

    return False

def close():
    if isOpen():
        imager = _get_win()
        imager.set("quit")

def delete_frames():
    if not isOpen():
        raise DS9Err('open')

    imager = _get_win()

    try:
        imager.set("frame delete all")
        return imager.set("frame new")
    except:
        raise DS9Err('delframe')
        
def get_region(coord):
    if not isOpen():
        raise DS9Err('open')

    imager = _get_win()

    try:
        regionstr = "regions -format saoimage -strip yes"
        if (coord != ''):
            if (coord != 'image'):
                regionstr = "regions -format ciao -strip yes -system " + str(coord)
            else:
                regionstr = "regions -format saoimage -strip yes -system image"
                                                    
        reg = imager.get(regionstr)
        reg = reg.replace(';','')
        return reg

    except:
        raise DS9Err('retreg')
    
def image(arr, newframe=False, tile=False):
    if not isOpen():
        doOpen()

    imager = _get_win()

    if newframe is True:
        try:
            imager.set("frame new")
            imager.set("frame last")
        except:
            raise DS9Err('newframe')
    try:
        if tile is True:
            imager.set("tile yes")
        else:
            imager.set("tile no")
    except:
        raise DS9Err('settile')
    time.sleep(1)
    try:
        # pyds9 expects shape[::-1] compared to DS9.py
        # therefore transpose the image before sending
        arr = numpy.asarray(arr, dtype=SherpaFloat)
        imager.set_np2arr(arr.T)
    except:
        raise # DS9Err('noimage')

def _set_wcs(keys):
    eqpos, sky, name = keys

    phys = ''
    wcs  = "OBJECT = '%s'\n" % name

    if eqpos is not None:
        wcrpix  = eqpos.crpix
        wcrval  = eqpos.crval
        wcdelt  = eqpos.cdelt

    if sky is not None:
        pcrpix = sky.crpix
        pcrval = sky.crval
        pcdelt = sky.cdelt

        # join together all strings with a '\n' between each
        phys = '\n'.join(["WCSNAMEP = 'PHYSICAL'",
                          "CTYPE1P = 'x       '",
                          'CRVAL1P = %.14E' % pcrval[0],
                          'CRPIX1P = %.14E' % pcrpix[0],
                          'CDELT1P = %.14E' % pcdelt[0],
                          "CTYPE2P = 'y       '",
                          'CRVAL2P = %.14E' % pcrval[1],
                          'CRPIX2P = %.14E' % pcrpix[1],
                          'CDELT2P = %.14E' % pcdelt[1]])

        if eqpos is not None:
            wcdelt = wcdelt * pcdelt
            wcrpix = ((wcrpix - pcrval) /
                      pcdelt + pcrpix )

    if eqpos is not None:
        # join together all strings with a '\n' between each
        wcs = wcs + '\n'.join(["RADECSYS = 'ICRS    '",
                               "CTYPE1  = 'RA---TAN'",
                               'CRVAL1  = %.14E' % wcrval[0],
                               'CRPIX1  = %.14E' % wcrpix[0],
                               'CDELT1  = %.14E' % wcdelt[0],
                               "CTYPE2  = 'DEC--TAN'",
                               'CRVAL2  = %.14E' % wcrval[1],
                               'CRPIX2  = %.14E' % wcrpix[1],
                               'CDELT2  = %.14E' % wcdelt[1]])

    # join the wcs and physical with '\n' between them and at the end
    return ('\n'.join([wcs,phys]) + '\n')

def wcs(keys):

    if not isOpen():
        raise DS9Err('open')

    imager = _get_win()
    info = _set_wcs( keys )

    try:
        # use stdin to pass the WCS info
        imager.set('wcs replace', info)
    except:
        raise DS9Err('setwcs')


def open():
    doOpen()

def set_region(reg, coord):
    if not isOpen():
        raise DS9Err('open')

    imager = _get_win()

    try:
        if (access(reg, R_OK) is True):
            imager.set("regions load " + "'" + reg + "'")
        else:
            # Assume region string has to be in CIAO format
            regions = reg.split(";")
            for region in regions:
                if (region != ''):
                    if (coord != ''):
                        imager.set("regions", str(coord) + ";" + region)
                    else:
                        imager.set("regions", region)
    except:
        raise DS9Err('badreg', str(reg))

def xpaget(arg):
    if not isOpen():
        raise DS9Err('open')
    imager = _get_win()
    return imager.get(arg)

def xpaset(arg, data=None):
    if not isOpen():
        raise DS9Err('open')
    imager = _get_win()
    return imager.set(arg, data)
