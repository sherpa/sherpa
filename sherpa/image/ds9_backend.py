#
#  Copyright (C) 2007, 2016-2017, 2021, 2026
#  Smithsonian Astrophysical Observatory
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

from os import access, R_OK
import time

from sherpa.astro.io.wcs import WCS
from sherpa.utils.err import DS9Err

from . import DS9

imager = DS9.DS9Win(DS9._DefTemplate, False)


# TODO: the except blocks would ideally catch explicit errors; the present
#       catch-anything approach means that we lose potentially-useful
#       information on the type of error.
#

def close():
    if imager.isOpen():
        imager.xpaset("quit")


def delete_frames():
    if not imager.isOpen():
        raise DS9Err('open')
    try:
        imager.xpaset("frame delete all")
        return imager.xpaset("frame new")
    except:
        raise DS9Err('delframe')


def get_region(coord):
    if not imager.isOpen():
        raise DS9Err('open')
    try:
        regionstr = "regions -format saoimage -strip yes"
        if coord != '':
            if coord != 'image':
                regionfmt = 'ciao'
            else:
                regionfmt = 'saoimage'

            regionstr = "regions -format {} ".format(regionfmt) + \
                        "-strip yes -system {}".format(coord)

        regionstr = imager.xpaget(regionstr)
        return regionstr

    except:
        raise DS9Err('retreg')


def image(arr, newframe=False, tile=False):
    if not imager.isOpen():
        imager.doOpen()
    # Create a new frame if the user requested it, *or* if
    # there happen to be no DS9 frames.
    if newframe or imager.xpaget("frame all") == "\n":
        try:
            imager.xpaset("frame new")
            imager.xpaset("frame last")
        except:
            raise DS9Err('newframe')
    try:
        if tile:
            imager.xpaset("tile yes")
        else:
            imager.xpaset("tile no")
    except:
        raise DS9Err('settile')
    time.sleep(1)
    try:
        imager.showArray(arr)
    except:
        raise DS9Err('noimage')


def _set_wcs(keys: tuple[WCS | None, WCS | None, str]) -> str:
    """Convert the settings into a string to send via XPA"""

    eqpos, sky, name = keys

    # DS9 can be very particular about the WCS settings, so attempt to
    # set everything to avoid problems with ds9 preference settings.
    #
    out = [f"OBJECT = '{name}'"]

    if sky is not None:
        pcrpix = sky.crpix
        pcrval = sky.crval
        pcdelt = sky.cdelt

        crval = [f'{pcrval[0]:.14E}', f'{pcrval[1]:.14E}']
        crpix = [f'{pcrpix[0]:.14E}', f'{pcrpix[1]:.14E}']
        cdelt = [f'{pcdelt[0]:.14E}', f'{pcdelt[1]:.14E}']

    else:
        # Ensure we have the identity transform defined.
        crval = ['1.0', '1.0']
        crpix = ['1.0', '1.0']
        cdelt = ['1.0', '1.0']

    out.extend(['WCSAXESP = 2',
                "WCSNAMEP = 'PHYSICAL'",
                "CTYPE1P  = 'x       '",
                "CTYPE2P  = 'y       '",
                f'CRVAL1P  = {crval[0]}',
                f'CRPIX1P  = {crpix[0]}',
                f'CDELT1P  = {cdelt[0]}',
                f'CRVAL2P  = {crval[1]}',
                f'CRPIX2P  = {crpix[1]}',
                f'CDELT2P  = {cdelt[1]}'])

    if eqpos is not None:

        wcrval = eqpos.crval
        if sky is not None:
            wcdelt = eqpos.cdelt * sky.cdelt
            wcrpix = (eqpos.crpix - sky.crval) / sky.cdelt + sky.crpix
        else:
            wcdelt = eqpos.cdelt
            wcrpix = eqpos.crpix

        out.extend(['WCSAXES = 2',
                    "RADECSYS = 'ICRS    '",
                    "CTYPE1  = 'RA---TAN'",
                    "CTYPE2  = 'DEC--TAN'",
                    f'CRVAL1  = {wcrval[0]:.14E}',
                    f'CRPIX1  = {wcrpix[0]:.14E}',
                    f'CDELT1  = {wcdelt[0]:.14E}',
                    f'CRVAL2  = {wcrval[1]:.14E}',
                    f'CRPIX2  = {wcrpix[1]:.14E}',
                    f'CDELT2  = {wcdelt[1]:.14E}'])

    # Ensure the string ends with a new-line
    out.append('')
    return '\n'.join(out)


def wcs(keys: tuple[WCS | None, WCS | None, str]) -> None:

    if not imager.isOpen():
        raise DS9Err('open')

    info = _set_wcs(keys)

    try:
        imager.xpaset('wcs replace', info)
    except:
        raise DS9Err('setwcs')


def open():
    imager.doOpen()


def set_region(reg, coord):
    if not imager.isOpen():
        raise DS9Err('open')
    try:
        # Assume a region file defines everything correctly
        if access(reg, R_OK):
            imager.xpaset("regions load " + "'" + reg + "'")
        else:
            # Assume region string has to be in CIAO format
            regions = reg.split(";")
            for region in regions:
                if region != '':
                    if coord != '':
                        data = str(coord) + ";" + region
                    else:
                        data = region

                    imager.xpaset("regions", data=data)

    except:
        raise DS9Err('badreg', str(reg))


def xpaget(arg):
    if not imager.isOpen():
        raise DS9Err('open')
    return imager.xpaget(arg)


def xpaset(arg, data=None):
    if not imager.isOpen():
        raise DS9Err('open')
    return imager.xpaset(arg, data)
