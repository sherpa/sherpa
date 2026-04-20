#
#  Copyright (C) 2009,2010,2016,2025
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

import numpy as np

from sherpa.astro.io.wcs import WCS


imager = None
"""The DS9 window, or None"""

# If this file is xxx_backend.py then 'name = "xxx"'.
name: str = "dummy"
"""The name of the backend."""


def close() -> None:
    """Stop the image viewer."""
    pass


def delete_frames() -> None:
    """Delete all the frames open in the image viewer."""
    pass


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
    return ""


def image(array: np.ndarray,
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
    pass


def wcs(keys: tuple[WCS | None, WCS | None, str]) -> None:
    """Send the WCS informatiom to the image viewer.

    Parameters
    ----------
    keys
       The eqpos and sky transforms, and the name of the display.

    """
    pass


def open() -> None:
    """Start the image viewer."""
    pass


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
    pass


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
    return ""


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
    pass
