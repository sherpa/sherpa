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

"""The no-op backend when DS9 is not available.

.. versionchanged:: 4.19.0
   The API has slightly changed to better match the DS9 version
   (types have been added and the arguments sent to the `wcs` routine
   have changed.

"""

import numpy as np

from sherpa.astro.io.wcs import WCS


imager = None
"""The DS9 window, or None"""


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
    array : np.ndarray
       The pixel values
    newframe : bool, optional
       Should the pixels be displayed in a new frame?
    tile : bool, optional
       Should the display be tiled?

    """
    pass


def wcs(eqpos: WCS | None,
        sky: WCS | None,
        name: str
        ) -> None:
    """Send the WCS information to the image viewer.

    .. versionchanged:: 4.19.0
       The arguments have changed.

    Parameters
    ----------
    eqpos, sky : WCS | None
       The transforms, if available.
    name : str
       The name to display.

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
