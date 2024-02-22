#
#  Copyright (C) 2007, 2016, 2017, 2021, 2023, 2024
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

"""Control commumication with the DS9 backend using XPA directly."""

import os
import time

from sherpa.utils.err import DS9Err

from . import DS9

imager = DS9.DS9Win(DS9._DefTemplate, False)


def close():
    """Ensure the DS9 instance is closed."""
    if imager.isOpen():
        imager.xpaset("quit")


def delete_frames():
    """Remove all frames and create a new frame."""

    if not imager.isOpen():
        raise DS9Err('open')

    try:
        imager.xpaset("frame delete all")
        return imager.xpaset("frame new")

    except Exception as exc:
        raise DS9Err('delframe') from exc


def get_region(coord):
    """Return the current region as a string.

    Parameters
    ----------
    coord : str
       The coordinate setting to use. It may be empty.

    Returns
    -------
    region : str
       The current region.

    """

    if not imager.isOpen():
        raise DS9Err('open')

    try:
        regionstr = "regions -format saoimage -strip yes"
        if coord != '':
            if coord != 'image':
                regionfmt = 'ciao'
            else:
                regionfmt = 'saoimage'

            regionstr = f"regions -format {regionfmt} " + \
                        f"-strip yes -system {coord}"

        return imager.xpaget(regionstr)

    except Exception as exc:
        raise DS9Err('retreg') from exc


def image(arr, newframe=False, tile=False):
    """Display the data as an image in DS9.

    Parameters
    ----------
    arr : ndarray
       The pixel data. It is required to be 2D (Y, X) or 3D (Z, Y, X)
       ordering.
    newframe : bool, optional
       Should the image be displayed in a new frame?
    tile : bool, optional
       Should DS9 tiling mode be selected?

    """

    if not imager.isOpen():
        imager.doOpen()

    # Create a new frame if the user requested it, *or* if
    # there happen to be no DS9 frames.
    if newframe or imager.xpaget("frame all") == "\n":
        try:
            imager.xpaset("frame new")
            imager.xpaset("frame last")

        except Exception as exc:
            raise DS9Err('newframe') from exc

    try:
        if tile:
            imager.xpaset("tile yes")
        else:
            imager.xpaset("tile no")

    except Exception as exc:
        raise DS9Err('settile') from exc

    time.sleep(1)
    try:
        imager.showArray(arr)

    except Exception as exc:
        raise DS9Err('noimage') from exc


def _set_wcs(keys):
    """Convert the WCS information into FITS metadata.

    Parameters
    ----------
    keys : triple
       The (eqpos, sky, name) values, where eqpos and sky are None or
       a WCS object. The name field is a string used as the OBJECT
       name.

    Returns
    -------
    header : str
       The FITS metadata to send to represent the WCS data.

    """

    eqpos, sky, name = keys

    phys_keys = []
    wcs_keys = [f"OBJECT = '{name}'\n"]

    if eqpos is not None:
        wcrpix = eqpos.crpix
        wcrval = eqpos.crval
        wcdelt = eqpos.cdelt

    if sky is not None:
        pcrpix = sky.crpix
        pcrval = sky.crval
        pcdelt = sky.cdelt

        phys_keys = ["WCSNAMEP = 'PHYSICAL'",
                     "CTYPE1P = 'x       '",
                     f'CRVAL1P = {pcrval[0]:.14E}',
                     f'CRPIX1P = {pcrpix[0]:.14E}',
                     f'CDELT1P = {pcdelt[0]:.14E}',
                     "CTYPE2P = 'y       '",
                     f'CRVAL2P = {pcrval[1]:.14E}',
                     f'CRPIX2P = {pcrpix[1]:.14E}',
                     f'CDELT2P = {pcdelt[1]:.14E}']

        if eqpos is not None:
            wcdelt = wcdelt * pcdelt
            wcrpix = (wcrpix - pcrval) / pcdelt + pcrpix

    if eqpos is not None:
        wcs_keys += ["RADECSYS = 'ICRS    '",
                     "CTYPE1  = 'RA---TAN'",
                     f'CRVAL1  = {wcrval[0]:.14E}',
                     f'CRPIX1  = {wcrpix[0]:.14E}',
                     f'CDELT1  = {wcdelt[0]:.14E}',
                     "CTYPE2  = 'DEC--TAN'",
                     f'CRVAL2  = {wcrval[1]:.14E}',
                     f'CRPIX2  = {wcrpix[1]:.14E}',
                     f'CDELT2  = {wcdelt[1]:.14E}']

    # Adding an empty string ensures we end with \n
    return '\n'.join(wcs_keys + phys_keys + [""])


def wcs(keys):
    """Send the WCS data to DS9 for the current image.

    Parameters
    ----------
    keys : triple
       The (eqpos, sky, name) values, where eqpos and sky are None or
       a WCS object. The name field is a string used as the OBJECT
       name.

    """

    if not imager.isOpen():
        raise DS9Err('open')

    info = _set_wcs(keys)

    try:
        # use stdin to pass the WCS info
        imager.xpaset('wcs replace', info)

    except Exception as exc:
        raise DS9Err('setwcs') from exc


def open():
    """Start the DS9 instance (if not already started)."""
    imager.doOpen()


def set_region(reg, coord):
    """Send the region to DS9.

    Parameters
    ----------
    reg : str or filename
       The file containing the region or a string containing regions
       separated by semi-colon characters.

    """
    if not imager.isOpen():
        raise DS9Err('open')

    try:
        # Assume a region file defines everything correctly
        if os.access(reg, os.R_OK):
            imager.xpaset(f"regions load '{reg}'")
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

    except Exception as exc:
        raise DS9Err('badreg', str(reg)) from exc


def xpaget(arg):
    """Send a XPA query to DS9.

    Parameter
    ---------
    arg : str
       The XPA query.

    Returns
    -------
    result
       The response from DS9.

    """
    if not imager.isOpen():
        raise DS9Err('open')

    return imager.xpaget(arg)


def xpaset(arg, data=None):
    """Send an XPA command to DS9.

    Parameter
    ---------
    arg : str
       The XPA command.
    data : optional
       Any data the XPA command requires.

    Returns
    -------
    result
       The response from DS9.

    """

    if not imager.isOpen():
        raise DS9Err('open')

    return imager.xpaset(arg, data)
