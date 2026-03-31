#
#  Copyright (C) 2007, 2016, 2017, 2021, 2025
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

import numpy as np

from sherpa.astro.io.wcs import WCS
from sherpa.utils.err import DS9Err

from . import DS9


imager = DS9.DS9Win(template=DS9._DefTemplate, doOpen=False)

name: str = "ds9"


# The except blocks would ideally catch explicit errors; the present
# catch-anything approach means that we lose potentially-useful
# information on the type of error. To reduce the loss of information
# the new exceptions are linked to the original exception.
#
# For now the try/except blocks have been made to catch BaseException
# rather than Exception, since things like "ds9 being killed" is
# probably something we either want to catch or users are already
# relying on.
#

def close() -> None:
    """Stop the image viewer."""
    if imager.isOpen():
        imager.xpaset("quit")


def delete_frames() -> None:
    """Delete all the frames open in the image viewer."""
    if not imager.isOpen():
        raise DS9Err('open')
    try:
        imager.xpaset("frame delete all")
        imager.xpaset("frame new")
    except BaseException as exc:
        raise DS9Err('delframe') from exc


def get_region(coord: str) -> str:
    """Return the region defined in the image viewer.

    Parameters
    ----------
    coord : str
       The name of the coordinate system (the empty string means
       to use the current system).

    Returns
    -------
    region : str
       The region, or regions, or the empty string.

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

    except BaseException as exc:
        raise DS9Err('retreg') from exc


def image(arr: np.ndarray,
          newframe: bool = False,
          tile: bool = False
          ) -> None:
    """Send the data to the image viewer to display.

    Parameters
    ----------
    array
       The pixel values
    newframe
       Should the pixels be displayed in a new frame?
    tile
       Should the display be tiled?

    """
    if not imager.isOpen():
        imager.doOpen()

    try:
        # Create a new frame if the user requested it, *or* if
        # there happen to be no DS9 frames.
        if newframe or imager.xpaget("frame all") == "\n":
            imager.xpaset("frame new")
            imager.xpaset("frame last")

        if tile:
            imager.xpaset("tile yes")
        else:
            imager.xpaset("tile no")

        time.sleep(1)
        imager.showArray(arr)

    except BaseException as exc:
        raise DS9Err('newframe') from exc

def _set_wcs(eqpos: WCS | None, sky: WCS | None, name: str) -> str:
    """Convert the settings into a string to send via XPA"""

    # DS9 can be very particular about the WCS settings, so attempt to
    # set everything to avoid problems with ds9 preference settings.
    #
    out = [f"OBJECT = '{name}'"]

    if eqpos is not None:
        wcrpix = eqpos.crpix
        wcrval = eqpos.crval
        wcdelt = eqpos.cdelt

    if sky is not None:
        pcrpix = sky.crpix
        pcrval = sky.crval
        pcdelt = sky.cdelt

        out.extend(["WCSAXESP = 2",
                    "WCSNAMEP = 'PHYSICAL'",
                    "CTYPE1P = 'x       '",
                    "CTYPE2P = 'y       '",
                    f'CRVAL1P = {pcrval[0]:.14E}',
                    f'CRPIX1P = {pcrpix[0]:.14E}',
                    f'CDELT1P = {pcdelt[0]:.14E}',
                    f'CRVAL2P = {pcrval[1]:.14E}',
                    f'CRPIX2P = {pcrpix[1]:.14E}',
                    f'CDELT2P = {pcdelt[1]:.14E}'])

        if eqpos is not None:
            wcdelt = wcdelt * pcdelt
            wcrpix = (wcrpix - pcrval) / pcdelt + pcrpix

    if eqpos is not None:
        out.extend(["WCSAXES = 2",
                    "RADECSYS = 'ICRS    '",
                    # f"EQUINOX = {eqpos.equinox}",
                    "CTYPE1  = 'RA---TAN'",
                    "CTYPE2  = 'DEC--TAN'",
                    f'CRVAL1  = {wcrval[0]:.14E}',
                    f'CRPIX1  = {wcrpix[0]:.14E}',
                    f'CDELT1  = {wcdelt[0]:.14E}',
                    f'CRVAL2  = {wcrval[1]:.14E}',
                    f'CRPIX2  = {wcrpix[1]:.14E}',
                    f'CDELT2  = {wcdelt[1]:.14E}'])

    # Ensure the string ends with a new-line.
    return ('\n'.join(out) + '\n')


def wcs(keys: tuple[WCS | None, WCS | None, str]) -> None:
    """Send the WCS informatiom to the image viewer.

    Parameters
    ----------
    keys
       The eqpos and sky transforms, and the name of the display.

    """

    if not imager.isOpen():
        raise DS9Err('open')

    info = _set_wcs(keys[0], keys[1], keys[2])

    try:
        # use stdin to pass the WCS info
        imager.xpaset('wcs replace', info)
    except BaseException as exc:
        raise DS9Err('setwcs') from exc


def open() -> None:
    """Start the image viewer."""
    imager.doOpen()


def set_region(reg: str, coord: str) -> None:
    """Set the region to display in the image viewer.

    Parameters
    ----------
    reg : str
       The region to display.
    coord : str
       The name of the coordinate system (the empty string means
       to use the current system).

    """
    if not imager.isOpen():
        raise DS9Err('open')
    try:
        # Assume a region file defines everything correctly
        if access(reg, R_OK):
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

    except BaseException as exc:
        raise DS9Err('badreg', str(reg)) from exc


def xpaget(arg: str) -> str:
    """Query the image viewer via XPA.

    Retrieve the results of a query to the image viewer.

    Parameters
    ----------
    arg : str
       A command to send to the image viewer via XPA.

    Returns
    -------
    returnval : str

    """
    if not imager.isOpen():
        raise DS9Err('open')
    return imager.xpaget(arg)


def xpaset(arg: str, data: str | bytes | None = None) -> None:
    """Send the image viewer a command via XPA.

    Send a command to the image viewer.

    Parameters
    ----------
    arg : str
       A command to send to the image viewer via XPA.
    data : optional
       The data for the command.

    """
    if not imager.isOpen():
        raise DS9Err('open')
    imager.xpaset(arg, data)
