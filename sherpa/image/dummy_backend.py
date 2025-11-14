#
#  Copyright (C) 2009, 2010, 2016, 2026
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


imager = None
"""The DS9 window, or None"""


def close():
    """Stop the image viewer."""
    pass


def delete_frames():
    """Delete all the frames open in the image viewer."""
    pass


def get_region(coord):
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
    pass


def image(array,
          newframe=False,
          tile=False
          ):
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
    pass


def wcs(keys):
    """Send the WCS information to the image viewer.

    Parameters
    ----------
    keys
       The eqpos and sky transforms, and the name of the display.

    """
    pass


def open():
    """Start the image viewer."""
    pass


def set_region(reg, coord):
    """Set the region to display in the image viewer.

    Parameters
    ----------
    reg : str
       The region to display.
    coord : str
       The name of the coordinate system (the empty string means
       to use the current system).

    """
    pass


def xpaget(arg):
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
    pass


def xpaset(arg, data=None):
    """Send the image viewer a command via XPA.

    Send a command to the image viewer.

    Parameters
    ----------
    arg : str
       A command to send to the image viewer via XPA.
    data : optional
       The data for the command.

    """
    pass
