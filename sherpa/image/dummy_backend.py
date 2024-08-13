#
#  Copyright (C) 2009, 2010, 2016, 2024
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

"""Used when DS9 is not available."""

imager = None
"""The DS9 instance (if set)."""


def close():
    """Ensure the DS9 instance is closed."""
    pass


def delete_frames():
    """Remove all frames and create a new frame."""
    pass


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
    pass


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


def wcs(keys):
    """Send the WCS data to DS9 for the current image.

    Parameters
    ----------
    keys : triple
       The (eqpos, sky, name) values, where eqpos and sky are None or
       a WCS object. The name field is a string used as the OBJECT
       name.

    """


def open():
    """Start the DS9 instance (if not already started)."""
    pass


def set_region(reg, coord):
    """Send the region to DS9.

    Parameters
    ----------
    reg : str or filename
       The file containing the region or a string containing regions
       separated by semi-colon characters.

    """
    pass


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
    pass


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
    pass
