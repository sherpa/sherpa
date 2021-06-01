#
#  Copyright (C) 2008, 2015, 2016, 2017, 2018, 2019, 2020, 2021
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

"""Classes for storing, inspecting, and manipulating astronomical data sets.

The two types of Astronomical data supported in this module are
two-dimensional images (:py:class:`DataIMG`) and X-ray spectra
(:py:class:`DataPHA`), along with associated response information
(:py:class:`DataARF` and :py:class:`DataRMF`). These objects can be
constructed directly or read from :term:`FITS` files using the
:py:mod:`sherpa.astro.io` routines.

Both types of data extend the capabilities of the
:py:class:`sherpa.data.Data` class:

- using geometric shapes (regions) to filter images;

- support different units for filtering images (logical, physical, and
  :term:`WCS`), depending on the available metadata;

- support different analysis units for filtering and display for
  :term:`PHA` files (channels, energy, and wavelengths);

- dynamically re-bin PHA data to improve the signal to noise (grouping
  and quality);

- and automatically support one or more spectra that define the
  background for the observation (for PHA files) that can then be
  subtracted from the data or a background model fit to them.

Notes
-----

Some functionality depends on the presence of the region and grouping
Sherpa modules, which are optional components of Sherpa.

Notebook support
----------------

The Data objects support the rich display protocol of IPython, with
HTML display of a table of information highlighting the relevant data
and, for some classes, SVG images. Examples can be found at
[AstroNoteBook]_.

References
----------

.. [AstroNoteBook] https://sherpa.readthedocs.io/en/latest/NotebookSupport.html

Examples
--------

Read in a 2D dataset from the file 'clus.fits' and then filter it to
only use those pixels that lie within 45 units from the physical
coordinate 3150,4515:

>>> from sherpa.astro.io import read_image
>>> img = read_image('clus.fits')
>>> img.set_coord('physical')
>>> img.notice2d('circle(3150,4515,45)')

Read in a PHA dataset from the file 'src.pi', subtract the background,
filter to only use the data 0.5 to 7 keV, and re-group the data within
this range to have at least 20 counts per group:

>>> from sherpa.astro.io import read_pha
>>> pha = read_pha('src.pi')
>>> pha.subtract()
>>> pha.set_analysis('energy')
>>> pha.notice(0.5, 7)
>>> pha.group_counts(20, tabStops=~pha.mask)

"""

import os.path
import logging
import warnings

import numpy

from sherpa.data import Data1DInt, Data2D, Data, Data2DInt, Data1D, \
    IntegratedDataSpace2D, _check
from sherpa.models.regrid import EvaluationSpace1D
from sherpa.utils.err import DataErr, ImportErr
from sherpa.utils import SherpaFloat, pad_bounding_box, interpolate, \
    create_expr, parse_expr, bool_cast, rebin, filter_bins
from sherpa.utils import formatting
from sherpa.astro import hc

# There are currently (Sep 2015) no tests that exercise the code that
# uses the compile_energy_grid symbols.
from sherpa.astro.utils import arf_fold, rmf_fold, filter_resp, \
    compile_energy_grid, do_group, expand_grouped_mask

info = logging.getLogger(__name__).info
warning = logging.getLogger(__name__).warning

regstatus = False
try:
    from sherpa.astro.utils._region import Region
    regstatus = True
except ImportError:
    warning('failed to import sherpa.astro.utils._region; Region routines ' +
            'will not be available')

groupstatus = False
try:
    import group as pygroup
    groupstatus = True
except ImportError:
    groupstatus = False
    warning('the group module (from the CIAO tools package) is not ' +
            'installed.\nDynamic grouping functions will not be available.')


__all__ = ('DataARF', 'DataRMF', 'DataPHA', 'DataIMG', 'DataIMGInt', 'DataRosatRMF')


def _notice_resp(chans, arf, rmf):
    bin_mask = None

    if rmf is not None and arf is not None:

        bin_mask = rmf.notice(chans)
        if len(rmf.energ_lo) == len(arf.energ_lo):
            arf.notice(bin_mask)

        # If the response is mis-matched, determine which energy bins in the
        # RMF correspond to energy bins in the ARF and which are noticed.
        # Propogate the noticed RMF energy bins to the ARF energy  bins.
        elif len(rmf.energ_lo) < len(arf.energ_lo):
            arf_mask = None
            if bin_mask is not None:
                arf_mask = numpy.zeros(len(arf.energ_lo), dtype=bool)
                for ii, val in enumerate(bin_mask):
                    if val:
                        los = (rmf.energ_lo[ii],)
                        his = (rmf.energ_hi[ii],)
                        grid = (arf.energ_lo, arf.energ_hi)
                        idx = filter_bins(los, his, grid).nonzero()[0]
                        arf_mask[idx] = True
            arf.notice(arf_mask)

    else:
        if rmf is not None:
            bin_mask = rmf.notice(chans)
        if arf is not None:
            arf.notice(bin_mask)


def display_header(header, key):
    """Return the header value for display by _repr_html

    The value is not displayed if it doesn't exist, is None,
    is empty, or is the string 'NONE'. This is intended for
    PHA responses.

    Parameters
    ----------
    header : dict-like
    key : str
        The key to display

    Returns
    -------
    value : None or value
        The value to display, or None.

    Notes
    -----
    It is not clear if the Meta class is intended to only store
    string values or not. Limited protection is provided in case
    the value stored is not a string.
    """

    try:
        val = header[key]
    except KeyError:
        return None

    # Unclear if this can happen
    if val is None:
        return None

    # The metadata value is not guaranteed to be a string
    try:
        val = val.strip()
        if val in ['', 'NONE']:
            return None
    except AttributeError:
        pass

    return val


def make_metadata(header, items):
    """Create the metadata table.

    Parameters
    ----------
    header : dict-like
        The header.
    items : list of (str, str)
        The keys to display (in order), if set. The first element
        is the key name, and the second is the label in the header
        to display.

    Returns
    -------
    meta : list of (str, str) or None
        The two-element table rows to display. If no rows matched
        return None.

    """

    meta = []
    for key, desc in items:
        val = display_header(header, key)
        if val is None:
            continue

        meta.append((desc, val))

    if len(meta) == 0:
        return None

    return meta


def _extract_fields(obj, stop, summary, open_block=True):
    """Extract the fields up until the stop field.

    Parameters
    ----------
    obj : Data instance
        It has to have a _fields attribute
    stop : str
        The attribute at which to stop (and is not included).
    summary : str
        The label for the details tab.
    open_block : bool, optional
        Is the details tab open or closed?

    Returns
    -------
    html : str
        The HTML for this section.
    """

    meta = []
    for f in obj._fields[1:]:
        if f == stop:
            break

        v = getattr(obj, f)
        if v is None:
            continue

        meta.append((f.upper(), v))

    return formatting.html_section(meta, summary=summary,
                                   open_block=open_block)


def html_pha(pha):
    """HTML representation: PHA"""

    from sherpa.astro.plot import DataPHAPlot, backend

    ls = []

    plotter = DataPHAPlot()
    plotter.prepare(pha)

    try:
        out = backend.as_html_plot(plotter, 'PHA Plot')
    except AttributeError:
        out = None

    if out is None:
        out = _extract_fields(pha, 'grouped', 'PHA Data')

    ls.append(out)

    # Summary properties
    meta = []
    if pha.name is not None and pha.name != '':
        meta.append(('Identifier', pha.name))

    if pha.exposure is not None:
        meta.append(('Exposure', '{:g} s'.format(pha.exposure)))

    meta.append(('Number of bins', len(pha.channel)))

    meta.append(('Channel range', '{} - {}'.format(int(pha.channel[0]),
                                                   int(pha.channel[-1]))))

    # Although assume the counts are integers, do not force this
    cmin = pha.counts.min()
    cmax = pha.counts.max()
    meta.append(('Count range', '{} - {}'.format(cmin, cmax)))

    if pha.background_ids != []:
        if pha.subtracted:
            msg = 'Subtracted'
        else:
            msg = 'Not subtracted'

        meta.append(('Background', msg))

    # Make sure show all groups (not just those that are within
    # the filter applied to the object).
    #
    if pha.grouping is not None:
        if pha.grouped:
            ngrp = pha.apply_grouping(pha.counts).size
            msg = 'Applied ({} groups)'.format(ngrp)
        else:
            msg = 'Not applied'

        meta.append(('Grouping', msg))

    # Should this only be displayed if a filter has been applied?
    #
    fexpr = pha.get_filter_expr()
    bintype = 'groups' if pha.grouped else 'channels'
    nbins = pha.get_dep(filter=True).size
    meta.append(('Using', '{} with {} {}'.format(fexpr, nbins, bintype)))

    ls.append(formatting.html_section(meta, summary='Summary',
                                      open_block=True))

    # TODO:
    #   correction factors

    # Display a subset of header values
    # - maybe don't display the FITLER if NONE
    # - how about RESPFILE / PHAFILE
    meta = make_metadata(pha.header,
                         [('TELESCOP', 'Mission or Satellite'),
                          ('INSTRUME', 'Instrument or Detector'),
                          ('GRATING', 'Grating type'),
                          ('ORDER', 'Diffraction order'),
                          ('FILTER', 'Instrument filter'),
                          ('OBJECT', 'Object'),
                          ('TITLE', 'Program description'),
                          ('DATE-OBS', 'Observation date'),
                          ('CREATOR', 'Program that created the PHA'),
                          ('CHANTYPE', 'The channel type'),
                          ('HDUCLAS2', 'Data stored'),
                          ('HDUCLAS3', 'Data format'),
                          ('HDUCLAS4', 'PHA format')])

    if meta is not None:
        ls.append(formatting.html_section(meta, summary='Metadata'))

    return formatting.html_from_sections(pha, ls)


def _calc_erange(elo, ehi):
    """Create the energy range information.

    Parameters
    ----------
    elo, ehi - NumPy array
        The low and high energy bins, in keV.

    Returns
    -------
    erange : str
        The string representation of the energy range

    """

    # Have we guaranteed the ordering here or not? Assuming
    # NumPy arrays.
    e1 = elo[0]
    e2 = ehi[-1]
    emin, emax = (e1, e2) if e1 <= e2 else (e2, e1)
    erange = '{:g} - {:g} keV'.format(emin, emax)

    # Randomly pick 1% as the cut-off for a constant bin width
    #
    de = numpy.abs(ehi - elo)
    demin = de.min()
    demax = de.max()
    if demin > 0.0:
        dedelta = (demax - demin) / demin
    else:
        dedelta = 1

    if dedelta <= 0.01:
        erange += ', bin size {:g} keV'.format(demax)
    else:
        erange += ', bin size {:g} - {:g} keV'.format(demin, demax)

    return erange


def _calc_wrange(wlo, whi):
    """Create the wavelength range information.

    Parameters
    ----------
    wlo, whi - NumPy array
        The low and high wavelength bins, in Angstroms.

    Returns
    -------
    wrange : str
        The string representation of the wavelength range

    """

    w1 = wlo[0]
    w2 = whi[-1]
    wmin, wmax = (w1, w2) if w1 <= w2 else (w2, w1)
    wrange = '{:g} - {:g} &#8491;'.format(wmin, wmax)

    # Randomly pick 1% as the cut-off for a constant bin width
    #
    dw = numpy.abs(whi - wlo)
    dwmin = dw.min()
    dwmax = dw.max()
    if dwmin > 0.0:
        dwdelta = (dwmax - dwmin) / dwmin
    else:
        dwdelta = 1

    if dwdelta <= 0.01:
        wrange += ', bin size {:g} &#8491;'.format(dwmax)
    else:
        wrange += ', bin size {:g} - {:g} &#8491;'.format(dwmin, dwmax)

    return wrange


def html_arf(arf):
    """HTML representation: ARF"""

    # Unlike the string representation, this provides extra
    # information (e.g. energy range covered). Should it include
    # any filters or masks? How about bin_lo/hi values?
    #
    # It also assumes the units are keV/cm^2 which is not
    # guaranteed.

    from sherpa.astro.plot import ARFPlot, backend

    ls = []

    plotter = ARFPlot()
    plotter.prepare(arf)

    try:
        out = backend.as_html_plot(plotter, 'ARF Plot')
    except AttributeError:
        out = None

    if out is None:
        out = _extract_fields(arf, 'exposure', 'ARF Data')

    ls.append(out)

    # Summary properties
    meta = []
    if arf.name is not None and arf.name != '':
        meta.append(('Identifier', arf.name))

    if arf.exposure is not None:
        meta.append(('Exposure', '{:g} s'.format(arf.exposure)))

    meta.append(('Number of bins', len(arf.specresp)))

    erange = _calc_erange(arf.energ_lo, arf.energ_hi)
    meta.append(('Energy range', erange))

    # repeat for wavelengths (without the energy threshold)
    #
    if arf.bin_lo is not None and arf.bin_hi is not None:
        wrange = _calc_wrange(arf.bin_lo, arf.bin_hi)
        meta.append(('Wavelength range', wrange))

    a1 = numpy.min(arf.specresp)
    a2 = numpy.max(arf.specresp)
    meta.append(('Area range', '{:g} - {:g} cm<sup>2</sup>'.format(a1, a2)))

    ls.append(formatting.html_section(meta, summary='Summary',
                                      open_block=True))

    # Display a subset of header values
    # - maybe don't display the FITLER if NONE
    # - how about RESPFILE / PHAFILE
    meta = make_metadata(arf.header,
                         [('TELESCOP', 'Mission or Satellite'),
                          ('INSTRUME', 'Instrument or Detector'),
                          ('GRATING', 'Grating type'),
                          ('ORDER', 'Diffraction order'),
                          ('TG_M', 'Diffraction order'),
                          ('FILTER', 'Instrument filter'),
                          ('OBJECT', 'Object'),
                          ('TITLE', 'Program description'),
                          ('DATE-OBS', 'Observation date'),
                          ('CREATOR', 'Program that created the ARF')])

    if meta is not None:
        ls.append(formatting.html_section(meta, summary='Metadata'))

    return formatting.html_from_sections(arf, ls)


def html_rmf(rmf):
    """HTML representation: RMF"""

    # See _html_arf for general comments

    ls = []

    svg = simulate_rmf_plot(rmf)
    if svg is not None:
        out = formatting.html_svg(svg, 'RMF Plot')
    else:
        out = _extract_fields(rmf, 'ethresh', 'RMF Data')

    ls.append(out)

    # Summary properties
    meta = []
    if rmf.name is not None and rmf.name != '':
        meta.append(('Identifier', rmf.name))

    meta.append(('Number of channels', rmf.detchans))
    meta.append(('Number of energies', len(rmf.energ_hi)))

    erange = _calc_erange(rmf.energ_lo, rmf.energ_hi)
    if rmf.ethresh is not None and rmf.energ_lo[0] <= rmf.ethresh:
        # Not entirely happy with the wording of this
        erange += ' (minimum threshold of {} was used)'.format(rmf.ethresh)

    meta.append(('Energy range', erange))

    meta.append(('Channel range', '{} - {}'.format(int(rmf.offset),
                                                   int(rmf.offset + rmf.detchans - 1))))

    # Could show the energy range as given by e_min/e_max but
    # is this useful?

    ls.append(formatting.html_section(meta, summary='Summary',
                                      open_block=True))

    # Display a subset of header values
    # - how about PHAFILE
    meta = make_metadata(rmf.header,
                         [('TELESCOP', 'Mission or Satellite'),
                          ('INSTRUME', 'Instrument or Detector'),
                          ('GRATING', 'Grating type'),
                          ('ORDER', 'Diffraction order'),
                          ('FILTER', 'Instrument filter'),
                          ('OBJECT', 'Object'),
                          ('TITLE', 'Program description'),
                          ('DATE-OBS', 'Observation date'),
                          ('CREATOR', 'Program that created the RMF'),
                          ('CHANTYPE', 'The channel type'),
                          ('LO_THRES', 'The minimum probability threshold'),
                          ('HDUCLAS3', 'Matrix contents')])

    if meta is not None:
        ls.append(formatting.html_section(meta, summary='Metadata'))

    return formatting.html_from_sections(rmf, ls)


def html_img(img):
    """HTML representation: IMG

    Special-case of the Data2D handling. It would be nice to re-use
    parts of the superclass behavior.
    """

    ls = []
    dtype = type(img).__name__

    svg = img_plot(img)
    if svg is not None:
        out = formatting.html_svg(svg, '{} Plot'.format(dtype))
        summary = ''
    else:
        # Only add prefix to summary if there's no plot
        summary = '{} '.format(dtype)

        # Summary properties
        #
        meta = []
        if img.name is not None and img.name != '':
            meta.append(('Identifier', img.name))

        # shape is better defined for DataIMG than Data2D
        meta.append(('Shape',
                     ('{1} by {0} pixels'.format(*img.shape))))

        meta.append(('Number of bins', len(img.y)))

        # Rely on the _fields ordering, ending at shape
        for f in img._fields[1:]:
            if f == 'shape':
                break

            meta.append((f.upper(), getattr(img, f)))

        if img.staterror is not None:
            meta.append(('Statistical error', img.staterror))

        if img.syserror is not None:
            meta.append(('Systematic error', img.syserror))

        out = formatting.html_section(meta, summary=summary + 'Data',
                                      open_block=True)

    ls.append(out)

    # Add coordinate-system information. The WCS structure in Sherpa
    # is not really sufficient to identify the transform.
    #
    if img.sky is not None:
        meta = []
        meta.append(('Center pixel (logical)', img.sky.crpix))
        meta.append(('Center pixel (physical)', img.sky.crval))
        meta.append(('Pixel size', img.sky.cdelt))

        ls.append(formatting.html_section(meta,
                                          summary='Coordinates: {}'.format(img.sky.name)))

    if img.eqpos is not None:
        meta = []
        meta.append(('Center pixel (physical)', img.eqpos.crpix))
        # could convert to RA/Dec
        meta.append(('Center pixel (world)', img.eqpos.crval))
        meta.append(('Pixel size', img.eqpos.cdelt))

        meta.append(('Rotation', img.eqpos.crota))
        meta.append(('Epoch', img.eqpos.epoch))
        meta.append(('Equinox', img.eqpos.equinox))

        ls.append(formatting.html_section(meta,
                                          summary='Coordinates: {}'.format(img.eqpos.name)))

    meta = make_metadata(img.header,
                         [('TELESCOP', 'Mission or Satellite'),
                          ('INSTRUME', 'Instrument or Detector'),
                          ('FILTER', 'Instrument filter'),
                          ('OBJECT', 'Object'),
                          ('TITLE', 'Program description'),
                          ('OBSERVER', 'Observer'),
                          ('EXPOSURE', 'Exposure time'),
                          ('DATE-OBS', 'Observation date'),
                          ('CREATOR', 'Program that created the image')])

    if meta is not None:
        ls.append(formatting.html_section(meta, summary='Metadata'))

    return formatting.html_from_sections(img, ls)


def simulate_rmf_plot(rmf):
    """Create a plot which shows the response to monochromatic energies.

    The SVG of the plot is returned if matplotlib is selected as the
    backend. The choice of energies used to create the response to
    monochromatic energies is based on the data range (using log
    scaling).

    """

    from sherpa.models.basic import Delta1D
    from sherpa.plot import backend

    try:
        from matplotlib import pyplot as plt
    except ImportError:
        return None

    # X access
    #
    if rmf.e_min is None:
        x = numpy.arange(rmf.offset, rmf.detchans + rmf.offset)
        xlabel = 'Channel'
    else:
        x = 0.5 * (rmf.e_min + rmf.e_max)
        xlabel = 'Energy (keV)'

    # How many monochromatic lines to use
    #
    nlines = 5

    # for now let's just create log-spaced energies
    #
    elo, ehi = rmf.energ_lo, rmf.energ_hi
    l1 = numpy.log10(elo[0])
    l2 = numpy.log10(ehi[-1])
    dl = (l2 - l1) / (nlines + 1)

    lines = l1 + dl * numpy.arange(1, nlines + 1)
    energies = numpy.power(10, lines)

    mdl = Delta1D()

    def plotfunc():
        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.2)

        for energy in energies:
            mdl.pos = energy
            y = rmf.apply_rmf(mdl(elo, ehi))
            ax.plot(x, y, label='{:.2g} keV'.format(energy))

        # Try to get the legend centered nicely below the plot
        fig.legend(loc='center', ncol=nlines, bbox_to_anchor=(0.0, 0, 1, 0.1))

        ax.set_xlabel(xlabel)
        ax.set_title(rmf.name)
        ax.set_xscale('log')
        ax.set_yscale('log')

        return fig

    try:
        return backend.as_svg(plotfunc)
    except AttributeError:
        return None


def img_plot(img):
    """Display the image.

    The SVG of the plot is returned if matplotlib is selected as the
    backend.

    The eqpos/wcs coordinate system is not used; it uses physical
    instead. This greatly simplifies the plot (no need to handle WCS).

    """

    from sherpa.plot import backend

    try:
        from matplotlib import pyplot as plt
    except ImportError:
        return None

    # Apply filter and coordinate system
    #
    y = img.get_img()

    # extent is left, right, bottom, top and describes the
    # outer-edge of the pixels.
    #
    ny, nx = img.shape
    coord = img.coord
    if coord in ['physical', 'world']:
        x0, y0 = img._logical_to_physical(0.5, 0.5)
        x1, y1 = img._logical_to_physical(nx + 0.5, ny + 0.5)
        extent = (x0, x1, y0, y1)
        lbl = 'physical'
        cdelt = img.sky.cdelt
        aspect = 'equal' if cdelt[1] == cdelt[0] else 'auto'

    else:
        extent = (0.5, nx + 0.5, 0.5, ny + 0.5)
        aspect = 'equal'
        lbl = 'logical'

    # What is the filtered dataset?
    #
    if img.get_filter_expr() != '':
        x0, x1 = img.get_indep(filter=True)

        x0min, x0max = numpy.min(x0), numpy.max(x0)
        x1min, x1max = numpy.min(x1), numpy.max(x1)

        # Should add in half cdelt to padd these, but
        # it looks like it isn't necessary.
        filtered = (x0min, x1min, x0max, x1max)

    else:
        filtered = None

    def plotfunc():
        fig, ax = plt.subplots()

        im = ax.imshow(y, origin='lower', extent=extent, aspect=aspect)
        fig.colorbar(im, ax=ax)

        if filtered is not None:
            ax.set_xlim(filtered[0], filtered[2])
            ax.set_ylim(filtered[1], filtered[3])

        ax.set_xlabel('X ({})'.format(lbl))
        ax.set_ylabel('Y ({})'.format(lbl))
        if img.name is not None and img.name != '':
            ax.set_title(img.name)

        return fig

    try:
        return backend.as_svg(plotfunc)
    except AttributeError:
        return None


class DataOgipResponse(Data1DInt):
    """
    Parent class for OGIP responses, in particular ARF and RMF. This class implements some common validation code that
    inheriting classes can call in their initializers.

    Inheriting classes should override the protected class field `_ui_name` to provide a more specific label for user
    messages.
    """
    _ui_name = "OGIP Response"

    # FIXME For a future time when we'll review this code in a deeper way: we
    # could have better separation of concerns if the initializers of `DataARF`
    # and `DataRMF` did not rely on the `Data` initializer, and if the
    # class hierarchy was better organized (e.g. it looks like children must
    # not call their super's initializer.  Also, I'd expect validation to
    # happen in individual methods rather than in a large one, and nested ifs
    # should be avoided if possible.
    #
    # The shift to creating a warning message instead of raising an
    # error has made this messier.
    #
    def _validate_energy_ranges(self, label, elo, ehi, ethresh):
        """Check the lo/hi values are > 0, handling common error case.

        Several checks are made, to make sure the parameters follow
        the OGIP standard. At present a failed check can result in
        either a warning message being logged, or an error raised.
        It was felt that raising an error in all cases would not be
        helpful to a user, who can't (easily) change the response
        files.

        Parameters
        ----------
        label : str
            The response file identifier.
        elo, ehi : numpy.ndarray
            The input ENERG_LO and ENERG_HI arrays. They are assumed
            to be one-dimensional and have the same number of elements.
        ethresh : None or float, optional
            If None, then elo must be greater than 0. When set, the
            start bin can have a low-energy edge of 0; it is replaced
            by ethresh. If set, ethresh must be greater than 0.
            An error is raised if ethresh is larger than the upper-edge
            of the first bin (only if the lower edge has been replaced).

        Returns
        -------
        elo, ehi : numpy arrays
            The validated energy limits. These can be the input arrays
            or a copy of them. At present the ehi array is the same as
            the input array, but this may change in the future.

        Notes
        -----
        Only some of the constraints provided by the OGIP standard are
        checked here, since there are issues involving numerical effects
        (e.g. when checking that two bins do not overlap), as well as
        uncertainty over what possible  behavior is seen in released
        data products for missions. The current set of checks are:

          - ehi > elo for each bin
          - elo is monotonic (ascending or descending)
          - when emin is set, the lowest value in elo is >= 0,
            otherwise it is > 0.
          - ethresh (if set) is less than the minimum value in ENERG_HI

        """

        rtype = self._ui_name

        if elo.size != ehi.size:
            raise ValueError("The energy arrays must have the same size, not {} and {}" .format(elo.size, ehi.size))

        if ethresh is not None and ethresh <= 0.0:
            raise ValueError("ethresh is None or > 0")

        if (elo >= ehi).any():
            # raise DataErr('ogip-error', rtype, label,
            #               'has at least one bin with ENERG_HI < ENERG_LO')
            wmsg = "The {} '{}' ".format(rtype, label) + \
                   'has at least one bin with ENERG_HI < ENERG_LO'
            warnings.warn(wmsg)

        # if elo is monotonically increasing, all elements will be True
        #                         decreasing,                      False
        #
        # so the sum will be number of elements or 0
        #
        increasing = numpy.diff(elo, n=1) > 0.0
        nincreasing = increasing.sum()
        if nincreasing > 0 and nincreasing != len(increasing):
            # raise DataErr('ogip-error', rtype, label,
            #               'has a non-monotonic ENERG_LO array')
            wmsg = "The {} '{}' ".format(rtype, label) + \
                   'has a non-monotonic ENERG_LO array'
            warnings.warn(wmsg)

        if nincreasing == 0:
            startidx = -1
        else:
            startidx = 0

        e0 = elo[startidx]
        if ethresh is None:
            if e0 <= 0.0:
                raise DataErr('ogip-error', rtype, label,
                              'has an ENERG_LO value <= 0')
        else:
            # TODO: should this equality be replaced by an approximation test?
            if e0 == 0.0:

                if ehi[startidx] <= ethresh:
                    raise DataErr('ogip-error', rtype, label,
                                  'has an ENERG_HI value <= the replacement ' +
                                  'value of {}'.format(ethresh))

                elo = elo.copy()
                elo[startidx] = ethresh
                wmsg = "The minimum ENERG_LO in the " + \
                       "{} '{}' was 0 and has been ".format(rtype, label) + \
                       "replaced by {}".format(ethresh)
                warnings.warn(wmsg)

            elif e0 < 0.0:
                # raise DataErr('ogip-error', rtype, label,
                #               'has an ENERG_LO value < 0')
                wmsg = "The {} '{}' ".format(rtype, label) + \
                       'has an ENERG_LO value < 0'
                warnings.warn(wmsg)

        return elo, ehi

    def _get_data_space(self, filter=False):
        return EvaluationSpace1D(self._lo, self._hi)


class DataARF(DataOgipResponse):
    """ARF data set.

    The ARF format is described in OGIP documents [1]_ and [2]_.

    Parameters
    ----------
    name : str
        The name of the data set; often set to the name of the file
        containing the data.
    energ_lo, energ_hi, specresp : numpy.ndarray
        The values of the ENERG_LO, ENERG_HI, and SPECRESP columns
        for the ARF. The ENERG_HI values must be greater than the
        ENERG_LO values for each bin, and the energy arrays must be
        in increasing or decreasing order.
    bin_lo, bin_hi : array or None, optional
    exposure : number or None, optional
        The exposure time for the ARF, in seconds.
    header : dict or None, optional
    ethresh : number or None, optional
        If set it must be greater than 0 and is the replacement value
        to use if the lowest-energy value is 0.0.

    Raises
    ------
    sherpa.utils.err.DataErr
        This is raised if the energy arrays do not follow some of the
        OGIP standards.

    Notes
    -----
    There is limited checking that the ARF matches the OGIP standard,
    but as there are cases of released data products that do not follow
    the standard, these checks can not cover all cases.

    References
    ----------

    .. [1] "The Calibration Requirements for Spectral Analysis (Definition of RMF and ARF file formats)", https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html

    .. [2] "The Calibration Requirements for Spectral Analysis Addendum: Changes log", https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002a/cal_gen_92_002a.html

    """
    _ui_name = "ARF"
    _fields = ("name", "energ_lo", "energ_hi", "specresp", "bin_lo", "bin_hi", "exposure", "ethresh")

    def _get_specresp(self):
        return self._specresp

    def _set_specresp(self, val):
        self._specresp = val
        self._rsp = val

    specresp = property(_get_specresp, _set_specresp)

    def __init__(self, name, energ_lo, energ_hi, specresp, bin_lo=None,
                 bin_hi=None, exposure=None, header=None, ethresh=None):
        self.specresp = specresp
        self.bin_lo = bin_lo
        self.bin_hi = bin_hi
        self.exposure = exposure
        self.header = {} if header is None else header
        self.ethresh = ethresh
        energ_lo, energ_hi = self._validate_energy_ranges(name, energ_lo, energ_hi, ethresh)
        self._lo, self._hi = energ_lo, energ_hi
        self.energ_lo = energ_lo
        self.energ_hi = energ_hi
        Data1DInt.__init__(self, name, energ_lo, energ_hi, specresp)

    def __str__(self):
        # Print the metadata first
        try:
            ss = Data.__str__(self)
        except:
            ss = self._fields
        return ss

    def _repr_html_(self):
        """Return a HTML (string) representation of the ARF
        """
        return html_arf(self)

    def __setstate__(self, state):
        if 'header' not in state:
            self.header = {}
        self.__dict__.update(state)

        if '_specresp' not in state:
            self.__dict__['_specresp'] = state.get('specresp', None)
            self.__dict__['_rsp'] = state.get('specresp', None)

    def apply_arf(self, src, *args, **kwargs):
        "Fold the source array src through the ARF and return the result"

        # an external function must be called so all ARFs go through
        # a single entry point in order for caching to 'work'
        model = arf_fold(src, self._rsp)

        # Rebin the high-res source model folded through ARF down to the size
        # the PHA or RMF expects.
        if args != ():
            (arf, rmf) = args
            if rmf != () and len(arf[0]) > len(rmf[0]):
                model = rebin(model, arf[0], arf[1], rmf[0], rmf[1])

        return model

    def notice(self, bin_mask=None):
        self._rsp = self.specresp
        self._lo = self.energ_lo
        self._hi = self.energ_hi
        if bin_mask is not None:
            self._rsp = self.specresp[bin_mask]
            self._lo = self.energ_lo[bin_mask]
            self._hi = self.energ_hi[bin_mask]

    def get_indep(self, filter=False):
        return (self._lo, self._hi)

    def get_dep(self, filter=False):
        return self._rsp

    def get_xlabel(self):
        return 'Energy (keV)'

    def get_ylabel(self):
        from sherpa.plot import backend
        return 'cm' + backend.get_latex_for_string('^2')


class DataRMF(DataOgipResponse):
    """RMF data set.

    The RMF format is described in OGIP documents [1]_ and [2]_.

    Parameters
    ----------
    name : str
        The name of the data set; often set to the name of the file
        containing the data.
    detchans : int
    energ_lo, energ_hi : array
        The values of the ENERG_LO, ENERG_HI, and SPECRESP columns
        for the ARF. The ENERG_HI values must be greater than the
        ENERG_LO values for each bin, and the energy arrays must be
        in increasing or decreasing order.
    n_grp, f_chan, n_chan, matrix : array-like
    offset : int, optional
    e_min, e_max : array-like or None, optional
    header : dict or None, optional
    ethresh : number or None, optional
        If set it must be greater than 0 and is the replacement value
        to use if the lowest-energy value is 0.0.

    Notes
    -----
    There is limited checking that the RMF matches the OGIP standard,
    but as there are cases of released data products that do not follow
    the standard, these checks can not cover all cases. If a check fails
    then a warning message is logged.

    References
    ----------

    .. [1] "The Calibration Requirements for Spectral Analysis (Definition of RMF and ARF file formats)", https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html

    .. [2] "The Calibration Requirements for Spectral Analysis Addendum: Changes log", https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002a/cal_gen_92_002a.html

    """
    _ui_name = "RMF"
    _fields = ("name", "detchans", "energ_lo", "energ_hi", "n_grp", "f_chan", "n_chan", "matrix", "offset", "e_min",
               "e_max", "ethresh")

    def __init__(self, name, detchans, energ_lo, energ_hi, n_grp, f_chan,
                 n_chan, matrix, offset=1, e_min=None, e_max=None,
                 header=None, ethresh=None):
        energ_lo, energ_hi = self._validate(name, energ_lo, energ_hi, ethresh)

        if offset < 0:
            raise ValueError("offset must be >=0, not {}".format(offset))
        self.energ_lo = energ_lo
        self.energ_hi = energ_hi
        self.offset = offset
        self.detchans = detchans
        self.e_min = e_min
        self.e_max = e_max
        self.header = {} if header is None else header
        self.n_grp = n_grp
        self.f_chan = f_chan
        self.n_chan = n_chan
        self.matrix = matrix
        self.ethresh = ethresh
        self._fch = f_chan
        self._nch = n_chan
        self._grp = n_grp
        self._rsp = matrix
        self._lo = energ_lo
        self._hi = energ_hi
        Data1DInt.__init__(self, name, energ_lo, energ_hi, matrix)

    def __str__(self):
        # Print the metadata first
        old = self._fields
        ss = old
        try:
            self._fields = tuple(filter((lambda x: x != 'header'),
                                        self._fields))
            ss = Data.__str__(self)
        finally:
            self._fields = old
        return ss

    def _repr_html_(self):
        """Return a HTML (string) representation of the RMF
        """
        return html_rmf(self)

    def __setstate__(self, state):
        if 'header' not in state:
            self.header = {}
        self.__dict__.update(state)

    def _validate(self, name, energy_lo, energy_hi, ethresh):
        """
        Validate energy ranges and, if necessary, make adjustments.
        Subclasses may override this method to perform different validations
        or skip validation altogether.

        Parameters
        ----------
        name : str
            The name/label of the current file
        energy_lo, energ_hi : NumPy array
            The lower bounds of the energy bins. The arrays must have the same size
        ethresh : float
            The lowest energy value

        Returns
        -------
        energy_lo, energ_hi : NumPy array
            The energy values to use for the bin boundaries
        """
        return self._validate_energy_ranges(name, energy_lo, energy_hi, ethresh)

    def apply_rmf(self, src, *args, **kwargs):
        "Fold the source array src through the RMF and return the result"

        # Rebin the high-res source model from the PHA down to the size
        # the RMF expects.
        if args != ():
            (rmf, pha) = args
            if pha != () and len(pha[0]) > len(rmf[0]):
                src = rebin(src, pha[0], pha[1], rmf[0], rmf[1])

        if len(src) != len(self._lo):
            raise TypeError("Mismatched filter between ARF and RMF " +
                            "or PHA and RMF")

        return rmf_fold(src, self._grp, self._fch, self._nch, self._rsp,
                        self.detchans, self.offset)

    def notice(self, noticed_chans=None):
        bin_mask = None
        self._fch = self.f_chan
        self._nch = self.n_chan
        self._grp = self.n_grp
        self._rsp = self.matrix
        self._lo = self.energ_lo
        self._hi = self.energ_hi
        if noticed_chans is not None:
            (self._grp, self._fch, self._nch, self._rsp,
             bin_mask) = filter_resp(noticed_chans, self.n_grp, self.f_chan,
                                     self.n_chan, self.matrix, self.offset)
            self._lo = self.energ_lo[bin_mask]
            self._hi = self.energ_hi[bin_mask]
        return bin_mask

    def get_indep(self, filter=False):
        return (self._lo, self._hi)

    def get_dep(self, filter=False):
        return self.apply_rmf(numpy.ones(self.energ_lo.shape, SherpaFloat))

    def get_xlabel(self):
        if (self.e_min is not None) and (self.e_max is not None):
            return 'Energy (keV)'
        return 'Channel'

    def get_ylabel(self):
        return 'Counts'


# FIXME There are places in the code that explicitly check if an object is an instance of sherpa.astro.data.DataRMF.
# So it's safer to make DataRosatRMF a subclass of the default class, although in principle they should be siblings
# and subclasses of the same superclass.
class DataRosatRMF(DataRMF):
    ui_name = "ROSAT RMF"

    def _validate(self, name, energy_lo, energy_hi, ethresh):
        return energy_lo, energy_hi


class DataPHA(Data1D):
    """PHA data set, including any associated instrument and background data.

    The PHA format is described in an OGIP document [1]_.

    Parameters
    ----------
    name : str
        The name of the data set; often set to the name of the file
        containing the data.
    channel, counts : array of int
        The PHA data.
    staterror, syserror : scalar or array or None, optional
        The statistical and systematic errors for the data, if
        defined.
    bin_lo, bin_hi : array or None, optional
    grouping : array of int or None, optional
    quality : array of int or None, optional
    exposure : number or None, optional
        The exposure time for the PHA data set, in seconds.
    backscal : scalar or array or None, optional
    areascal : scalar or array or None, optional
    header : dict or None, optional

    Attributes
    ----------
    name : str
        Used to store the file name, for data read from a file.
    channel
    counts
    staterror
    syserror
    bin_lo
    bin_hi
    grouping
    quality
    exposure
    backscal
    areascal

    Notes
    -----
    The original data is stored in the attributes - e.g. `counts` - and
    the data-access methods, such as `get_dep` and `get_staterror`,
    provide any necessary data manipulation to handle cases such as:
    background subtraction, filtering, and grouping.

    The handling of the AREASCAl value - whether it is a scalar or
    array - is currently in flux. It is a value that is stored with the
    PHA file, and the OGIP PHA standard ([1]_) describes the observed
    counts being divided by the area scaling before comparison to the
    model. However, this is not valid for Poisson-based statistics, and
    is also not how XSPEC handles AREASCAL ([2]_); the AREASCAL values
    are used to scale the exposure times instead. The aim is to add
    this logic to the instrument models in `sherpa.astro.instrument`,
    such as `sherpa.astro.instrument.RMFModelPHA`. The area scaling still
    has to be applied when calculating the background contribution to
    a spectrum, as well as when calculating the data and model values used
    for plots (following XSPEC so as to avoid sharp discontinuities where
    the area-scaling factor changes strongly).

    References
    ----------

    .. [1] "The OGIP Spectral File Format", https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html

    .. [2] Private communication with Keith Arnaud

    """
    _fields = ('name', 'channel', 'counts', 'staterror', 'syserror', 'bin_lo', 'bin_hi', 'grouping', 'quality',
               'exposure', 'backscal', 'areascal', 'grouped', 'subtracted', 'units', 'rate', 'plot_fac', 'response_ids',
               'background_ids')

    def _get_grouped(self):
        return self._grouped

    def _set_grouped(self, val):
        val = bool(val)

        if val and self.grouping is None:
            raise DataErr('nogrouping', self.name)

        if self._grouped == val:
            return

        # As the grouping status is being changed, we need to reset the mask
        # to be correct size, while still noticing groups within the filter
        #
        if numpy.iterable(self.mask):
            old_filter = self.get_filter(group=val)
            self._grouped = val
            self.ignore()
            for vals in parse_expr(old_filter):
                self.notice(*vals)

        self._grouped = val

    grouped = property(_get_grouped, _set_grouped,
                       doc='Are the data grouped?')

    def _get_subtracted(self):
        return self._subtracted

    def _set_subtracted(self, val):
        val = bool(val)
        if len(self._backgrounds) == 0:
            raise DataErr('nobkg', self.name)
        self._subtracted = val

    subtracted = property(_get_subtracted, _set_subtracted,
                          doc='Are the background data subtracted?')

    def _get_units(self):
        return self._units

    def _set_units(self, val):
        units = str(val).strip().lower()

        if units == 'bin':
            units = 'channel'

        if units.startswith('chan'):
            # Note: the names of these routines appear confusing because of the
            #       way group values are used
            self._to_channel = self._channel_to_group
            self._from_channel = self._group_to_channel
            units = 'channel'

        elif units.startswith('ener'):
            self._to_channel = self._energy_to_channel
            self._from_channel = self._channel_to_energy
            units = 'energy'

        elif units.startswith('wave'):
            self._to_channel = self._wavelength_to_channel
            self._from_channel = self._channel_to_wavelength
            units = 'wavelength'

        else:
            raise DataErr('bad', 'quantity', val)

        for id in self.background_ids:
            bkg = self.get_background(id)
            if bkg.get_response() != (None, None) or \
               (bkg.bin_lo is not None and bkg.bin_hi is not None):
                bkg.units = units

        self._units = units

    units = property(_get_units, _set_units,
                     doc="Units of the independent axis: one of 'channel', 'energy', 'wavelength'.")

    def _get_rate(self):
        return self._rate

    def _set_rate(self, val):
        self._rate = bool_cast(val)
        for id in self.background_ids:
            # TODO: shouldn't this store bool_cast(val) instead?
            self.get_background(id).rate = val

    rate = property(_get_rate, _set_rate,
                    doc='Quantity of y-axis: counts or counts/sec')

    def _get_plot_fac(self):
        return self._plot_fac

    def _set_plot_fac(self, val):
        self._plot_fac = int(val)
        for id in self.background_ids:
            self.get_background(id).plot_fac = val

    plot_fac = property(_get_plot_fac, _set_plot_fac,
                        doc='Number of times to multiply the y-axis ' +
                        'quantity by x-axis bin size')

    def _get_response_ids(self):
        return self._response_ids

    def _set_response_ids(self, ids):
        if not numpy.iterable(ids):
            raise DataErr('idsnotarray', 'response', str(ids))
        keys = self._responses.keys()
        for id in ids:
            if id not in keys:
                raise DataErr('badids', str(id), 'response', str(keys))
        ids = list(ids)
        self._response_ids = ids

    response_ids = property(_get_response_ids, _set_response_ids,
                            doc=('IDs of defined instrument responses ' +
                                 '(ARF/RMF pairs)'))

    def _get_background_ids(self):
        return self._background_ids

    def _set_background_ids(self, ids):
        if not numpy.iterable(ids):
            raise DataErr('idsnotarray', 'background', str(ids))
        keys = self._backgrounds.keys()
        for id in ids:
            if id not in keys:
                raise DataErr('badids', str(id), 'background', str(keys))
        ids = list(ids)
        self._background_ids = ids

    background_ids = property(_get_background_ids, _set_background_ids,
                              doc='IDs of defined background data sets')

    def __init__(self, name, channel, counts, staterror=None, syserror=None,
                 bin_lo=None, bin_hi=None, grouping=None, quality=None,
                 exposure=None, backscal=None, areascal=None, header=None):

        channel = _check(channel)
        counts = _check(counts)

        self.channel = channel
        self.counts = counts
        self.bin_lo = bin_lo
        self.bin_hi = bin_hi
        self.quality = quality
        self.grouping = grouping
        self.exposure = exposure
        self.backscal = backscal
        self.areascal = areascal
        self.header = {} if header is None else header
        self._grouped = (grouping is not None)
        # _original_groups is set False if the grouping is changed via
        # the _dynamic_groups method. This is currently only used by the
        # serialization code (sherpa.astro.ui.serialize) to determine
        # whether to write out the grouping data.
        #
        self._original_groups = True
        self._subtracted = False
        self._response_ids = []
        self._background_ids = []
        self._responses = {}
        self._backgrounds = {}
        self._rate = True
        self._plot_fac = 0
        self.units = 'channel'
        self.quality_filter = None
        Data1D.__init__(self, name, channel, counts, staterror, syserror)

    def __str__(self):
        # Print the metadata first
        old = self._fields
        ss = old
        try:
            self._fields = tuple(filter((lambda x: x != 'header'),
                                        self._fields))
            ss = Data.__str__(self)
        finally:
            self._fields = old
        return ss

    def _repr_html_(self):
        """Return a HTML (string) representation of the PHA
        """
        return html_pha(self)

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_to_channel']
        del state['_from_channel']
        return state

    def __setstate__(self, state):
        self._background_ids = state['_background_ids']
        self._backgrounds = state['_backgrounds']
        self._set_units(state['_units'])

        if 'header' not in state:
            self.header = {}
        self.__dict__.update(state)

    primary_response_id = 1
    """The identifier for the response component when not set."""

    def set_analysis(self, quantity, type='rate', factor=0):
        """Set the units used when fitting and plotting spectral data.

        Parameters
        ----------
        quantity : {'channel', 'energy', 'wavelength'}
            The analysis setting.
        type : {'rate', 'counts'}, optional
            Do plots display a rate or show counts?
        factor : int, optional
           The Y axis of plots is multiplied by Energy^factor or
           Wavelength^factor before display. The default is 0.

        Raises
        ------
        sherpa.utils.err.DatatErr
           If the type argument is invalid, the RMF or ARF has the
           wrong size, or there in no response.

        See Also
        --------
        get_analysis

        Examples
        --------

        >>> pha.set_analysis('energy')

        >>> pha.set_analysis('wave', type='counts', factor=1)
        >>> pha.units
        'wavelength'

        """
        self.plot_fac = factor

        type = str(type).strip().lower()
        if not (type.startswith('counts') or type.startswith('rate')):
            raise DataErr("plottype", type, "'rate' or 'counts'")

        self.rate = (type == 'rate')

        arf, rmf = self.get_response()
        if rmf is not None and rmf.detchans != len(self.channel):
            raise DataErr("incompatibleresp", rmf.name, self.name)

        if (rmf is None and arf is None) and \
           (self.bin_lo is None and self.bin_hi is None) and \
           quantity != 'channel':
            raise DataErr('norsp', self.name)

        if rmf is None and arf is not None and quantity != 'channel' and \
           len(arf.energ_lo) != len(self.channel):
            raise DataErr("incompleteresp", self.name)

        self.units = quantity

    def get_analysis(self):
        """Return the units used when fitting spectral data.

        Returns
        -------
        setting : { 'channel', 'energy', 'wavelength' }
            The analysis setting.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain PHA data.
        sherpa.utils.err.IdentifierErr
           If the `id` argument is not recognized.

        See Also
        --------
        set_analysis

        Examples
        --------

        >>> is_wave = pha.get_analysis() == 'wavelength'

        """
        return self.units

    def _fix_response_id(self, id):
        if id is None:
            id = self.primary_response_id
        return id

    def get_response(self, id=None):
        """Return the response component.

        Parameters
        ----------
        id : int or str, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        Returns
        -------
        arf, rmf: sherpa.astro.data.DataARF,sherpa.astro.data.DataRMF instances or None
           The response, as an ARF and RMF. Either, or both,
           components can be None.

        See Also
        --------
        delete_response, get_arf, get_rmf, set_response

        """
        id = self._fix_response_id(id)
        return self._responses.get(id, (None, None))

    def set_response(self, arf=None, rmf=None, id=None):
        """Add or replace a response component.

        To remove a response use delete_response(), as setting arf and
        rmf to None here does nothing.

        Parameters
        ----------
        arf : sherpa.astro.data.DataARF instance or None, optional
           The ARF to add if any.
        rmf : sherpa.astro.data.DataRMF instance or None, optional
           The RMF to add, if any.
        id : int or str, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        See Also
        --------
        delete_response, get_response, set_arf, set_rmf

        """
        if (arf is None) and (rmf is None):
            return

        id = self._fix_response_id(id)
        self._responses[id] = (arf, rmf)
        ids = self.response_ids[:]
        if id not in ids:
            ids.append(id)
        self.response_ids = ids

    def delete_response(self, id=None):
        """Remove the response component.

        If the response component does not exist then the method
        does nothing.

        Parameters
        ----------
        id : int or str, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        See Also
        --------
        set_response

        """
        id = self._fix_response_id(id)
        self._responses.pop(id, None)
        ids = self.response_ids[:]
        ids.remove(id)
        self.response_ids = ids

    def get_arf(self, id=None):
        """Return the ARF from the response.

        Parameters
        ----------
        id : int or str, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        Returns
        -------
        arf: sherpa.astro.data.DataARF instance or None
           The ARF, if set.

        See Also
        --------
        get_response, get_rmf

        """
        return self.get_response(id)[0]

    def get_rmf(self, id=None):
        """Return the RMF from the response.

        Parameters
        ----------
        id : int or str, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        Returns
        -------
        rmf: sherpa.astro.data.DataRMF instance or None
           The RMF, if set.

        See Also
        --------
        get_arf, get_response

        """
        return self.get_response(id)[1]

    def set_arf(self, arf, id=None):
        """Add or replace the ARF in a response component.

        This replaces the existing ARF of the response, keeping the
        previous RMF (if set). Use the delete_response method to
        remove the response, rather than setting arf to None.

        Parameters
        ----------
        arf : sherpa.astro.data.DataARF instance
           The ARF to add.
        id : int or str, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        See Also
        --------
        delete_response, set_response, set_rmf

        """
        self.set_response(arf, self.get_rmf(id), id)

    def set_rmf(self, rmf, id=None):
        """Add or replace the RMF in a response component.

        This replaces the existing RMF of the response, keeping the
        previous ARF (if set). Use the delete_response method to
        remove the response, rather than setting rmf to None.

        Parameters
        ----------
        rmf : sherpa.astro.data.DataRMF instance
           The RMF to add.
        id : int or str, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        See Also
        --------
        delete_response, set_response, set_arf

        """
        self.set_response(self.get_arf(id), rmf, id)

    def get_specresp(self, filter=False):
        """Return the effective area values for the data set.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the ARF or not. The default is `False`.

        Returns
        -------
        arf : array
           The effective area values for the data set (or background
           component).

        """
        filter = bool_cast(filter)
        self.notice_response(False)
        arf, rmf = self.get_response()
        newarf = None

        if arf is not None and rmf is not None:
            specresp = arf.get_dep()
            elo, ehi = arf.get_indep()
            lo, hi = self._get_ebins(group=False)

            newarf = interpolate(lo, elo, specresp)
            newarf[newarf <= 0] = 1.

            if filter:
                newarf = self.apply_filter(newarf, self._middle)

        return newarf

    def _get_ebins(self, response_id=None, group=True):
        """Return the low and high edges of the independent axis.

        This method is badly named as it will return values in either
        channel or energy units, depending on the units setting and
        the associated response information. When the response
        includes a RMF then it returns the approximation of the
        mapping from channel space to energy - that is the E_MIN and
        E_MAX columns from the RMF EBOUNDS block rather than from the
        ENERG_LO and ENERG_HI columns from the MATRIX block.

        Parameters
        ----------
        response_id : int or None, optional
            The response to use when units are not "channel". The
            default is to use the default response identifier.
        group : bool, optional
            Should the current grouping setting be applied. This is
            only used if the "grouped" attribute is set.

        Returns
        -------
        lo, hi : ndarray
            The low and high edges of each bin, in either channels or
            keV: energy is used unless the units setting is channel or
            there is no associated response. If the group flag is set
            and the data set is grouped then it uses the grouping
            settings, otherwise the data is for each channel. No
            filtering is applied.

        See Also
        --------
        _get_indep

        Examples
        --------

        >>> pha.ungroup()
        >>> pha.units = 'channel'
        >>> clo, chi = pha._get_ebins()
        >>> (clo == pha.channel).all()
        True
        >>> (chi == clo + 1).all()
        True

        >>> pha.units = 'energy'
        >>> elo, ehi = pha._get_ebins()
        >>> elo.size == pha.channel.size
        True
        >>> elo[0:5]
        array([0.00146, 0.0146 , 0.0292 , 0.0438 , 0.0584 ])
        >>> (elo[1:] == ehi[:-1]).all()
        True

        >>> pha.group()
        >>> glo, ghi = pha._get_ebins()
        >>> glo[0:5]
        array([0.00146   , 0.2482    , 0.3066    , 0.46720001, 0.56940001])

        Note that the returned units are energy even if units is set
        to "wavelength":

        >>> pha.units = 'wave'
        >>> wlo, whi = pha._get_ebins()
        >>> (wlo == glo).all()

        """
        group = bool_cast(group)

        if self.units == 'channel':
            elo = self.channel
            ehi = self.channel + 1
        elif (self.bin_lo is not None) and (self.bin_hi is not None):
            elo = self.bin_lo
            ehi = self.bin_hi
            if (elo[0] > elo[-1]) and (ehi[0] > ehi[-1]):
                elo = hc / self.bin_hi
                ehi = hc / self.bin_lo
        else:
            arf, rmf = self.get_response(response_id)
            if rmf is not None:
                if (rmf.e_min is None) or (rmf.e_max is None):
                    raise DataErr('noenergybins', 'RMF')
                elo = rmf.e_min
                ehi = rmf.e_max
            elif arf is not None:
                elo = arf.energ_lo
                ehi = arf.energ_hi
            else:
                elo = self.channel
                ehi = self.channel + 1

        if self.grouped and group:
            elo = self.apply_grouping(elo, self._min)
            ehi = self.apply_grouping(ehi, self._max)

        # apply_grouping applies a quality filter to the output
        # but if we get here then there is no equivalent. This
        # is likely confusing, at best, but we don't have good
        # tests to check what we should be doing.
        #
        return (elo, ehi)

    def get_indep(self, filter=True):
        if filter:
            return (self.get_noticed_channels(),)

        return (self.channel,)

    def _get_indep(self, filter=False):
        """Return the low and high edges of the independent axis.

        Unlike _get_ebins, this returns values in the "native" space
        of the response - i.e. for a RMF, it returns the bounds from
        the MATRIX rather than EBOUNDS extension of the RMF - and not
        the approximation used in _get_ebins.

        Parameters
        ----------
        filter : bool, optional
            It is not clear what this option means.

        Returns
        -------
        lo, hi : ndarray
            The low and high edges of each bin, in either keV or
            Angstroms.

        Raises
        ------
        sherpa.utils.err.DataErr
            The data set does not contain a response.

        See Also
        --------
        _get_ebins

        Notes
        -----
        If the PHA file contains multiple responses then they are
        combined to create the overall grid.

        Examples
        --------

        >>> pha.units = 'energy'
        >>> elo, eho = pha._get_indep()
        >>> elo.shape
        (1090,)
        >>> pha.channel.shape
        (1024,)
        >>> elo[0:5]
        array([0.1 , 0.11, 0.12, 0.13, 0.14])
        >>> ehi[0:5]
        array([0.11      , 0.12      , 0.13      , 0.14      , 0.15000001])
        >>> (elo[1:] == ehi[:-1]).all()
        True

        >>> pha.units = 'wave'
        >>> wlo, who = pha._get_indep()
        >>> wlo[0:4]
        array([112.71289825, 103.32015848,  95.37245534,  88.56013348])
        >>> whi[0:4]
        array([123.98418555, 112.71289825, 103.32015848,  95.37245534])
        >>> (wlo[:-1] == whi[1:]).all()
        True

        """

        if (self.bin_lo is not None) and (self.bin_hi is not None):
            elo = self.bin_lo
            ehi = self.bin_hi
            if (elo[0] > elo[-1]) and (ehi[0] > ehi[-1]):
                if self.units == 'wavelength':
                    return (elo, ehi)

                elo = hc / self.bin_hi
                ehi = hc / self.bin_lo

        else:
            energylist = []
            for id in self.response_ids:
                arf, rmf = self.get_response(id)
                lo = None
                hi = None

                if rmf is not None:
                    lo = rmf.energ_lo
                    hi = rmf.energ_hi
                    if filter:
                        lo, hi = rmf.get_indep()

                elif arf is not None:
                    lo = arf.energ_lo
                    hi = arf.energ_hi
                    if filter:
                        lo, hi = arf.get_indep()

                energylist.append((lo, hi))

            if len(energylist) > 1:
                # TODO: This is only tested by test_eval_multi_xxx and not with
                # actual (i.e. real world) data
                elo, ehi, lookuptable = compile_energy_grid(energylist)
            elif (not energylist or
                  (len(energylist) == 1 and
                      numpy.equal(energylist[0], None).any())):
                raise DataErr('noenergybins', 'Response')
            else:
                elo, ehi = energylist[0]

        lo, hi = elo, ehi
        if self.units == 'wavelength':
            lo = hc / ehi
            hi = hc / elo

        return (lo, hi)

    def _channel_to_group(self, val):
        """Convert channel number to group number.

        For ungrouped data channel and group numbering are the
        same.
        """
        if not self.grouped:
            return val

        # The edge channels of each group.
        #
        lo = self.apply_grouping(self.channel, self._min)
        hi = self.apply_grouping(self.channel, self._max)

        val = numpy.asarray(val).astype(numpy.int_)
        res = []
        for v in val.flat:
            # could follow _energy_to_channel but for now go
            # with something simple
            if v < self.channel[0]:
                ans = self.channel[0]
            elif v > self.channel[-1]:
                ans = self.channel[-1]
            else:
                idx, = numpy.where((v >= lo) & (v <= hi))
                ans = idx[0] + 1

            res.append(ans)

        res = numpy.asarray(res, SherpaFloat)
        if val.shape == ():
            return res[0]

        return res

    def _group_to_channel(self, val, group=True, response_id=None):
        """Convert group number to channel number.

        For ungrouped data channel and group numbering are the
        same. The mid-point of each group is used (rounded down
        if not an integer).
        """

        if not self.grouped or not group:
            return val

        # The middle channel of each group.
        #
        mid = self.apply_grouping(self.channel, self._middle)

        # Convert to an integer (this keeps the channel within
        # the group).
        #
        mid = numpy.floor(mid)

        val = numpy.asarray(val).astype(numpy.int_) - 1
        try:
            return mid[val]
        except IndexError:
            raise DataErr('invalid group number: {}'.format(val))

    def _channel_to_energy(self, val, group=True, response_id=None):
        elo, ehi = self._get_ebins(response_id=response_id, group=group)
        val = numpy.asarray(val).astype(numpy.int_) - 1
        try:
            return (elo[val] + ehi[val]) / 2.0
        except IndexError:
            raise DataErr('invalidchannel', val)

    def _energy_to_channel(self, val):
        elo, ehi = self._get_ebins()

        # special case handling no noticed data (e.g. ignore_bad
        # removes all bins); assume if elo is empty then so is ehi.
        #
        if len(elo) == 0:
            raise DataErr('notmask')

        val = numpy.asarray(val)
        res = []
        for v in val.flat:
            if tuple(numpy.flatnonzero(elo <= v)) == ():
                if elo[0] > elo[-1] and ehi[0] > ehi[-1]:
                    res.append(SherpaFloat(len(elo)))
                else:
                    res.append(SherpaFloat(1))
            elif tuple(numpy.flatnonzero(ehi > v)) == ():
                if elo[0] > elo[-1] and ehi[0] > ehi[-1]:
                    res.append(SherpaFloat(1))
                else:
                    res.append(SherpaFloat(len(ehi)))
            elif tuple(numpy.flatnonzero((elo <= v) & (ehi > v)) + 1) != ():
                res.append(SherpaFloat(
                    numpy.flatnonzero((elo <= v) & (ehi > v)) + 1))
            elif (elo <= v).argmin() == (ehi > v).argmax():
                res.append(SherpaFloat((elo <= v).argmin()))
            else:
                raise DataErr("energytochannel", v)

        if val.shape == ():
            return res[0]

        return numpy.asarray(res, SherpaFloat)

    def _channel_to_wavelength(self, val, group=True, response_id=None):
        tiny = numpy.finfo(numpy.float32).tiny
        vals = numpy.asarray(self._channel_to_energy(val, group, response_id))
        if vals.shape == ():
            if vals == 0.0:
                vals = tiny
        else:
            vals[vals == 0.0] = tiny
        vals = hc / vals
        return vals

    def _wavelength_to_channel(self, val):
        """Convert a wavelength to group or channel number.

        Note that values of 0 or less are replaced by a
        small value (~1e-38).

        """
        tiny = numpy.finfo(numpy.float32).tiny
        vals = numpy.asarray(val)
        if vals.shape == ():
            if vals <= 0.0:
                vals = tiny
        else:
            vals[vals <= 0.0] = tiny
        vals = hc / vals
        return self._energy_to_channel(vals)

    default_background_id = 1
    """The identifier for the background component when not set."""

    def _fix_background_id(self, id):
        if id is None:
            id = self.default_background_id
        return id

    def get_background(self, id=None):
        """Return the background component.

        Parameters
        ----------
        id : int or str, optional
           The identifier of the background component. If it is None
           then the default background identifier is used.

        Returns
        -------
        bkg : sherpa.astro.data.DataPHA instance or None
           The background dataset. If there is no component then None
           is returned.

        See Also
        --------
        delete_background, set_background

        """
        id = self._fix_background_id(id)
        return self._backgrounds.get(id)

    def set_background(self, bkg, id=None):
        """Add or replace a background component.

        If the background has no grouping of quality arrays then they
        are copied from the source region. If the background has no
        response information (ARF or RMF) then the response is copied
        from the source region.

        Parameters
        ----------
        bkg : sherpa.astro.data.DataPHA instance
           The background dataset to add. This object may be changed
           by this method.
        id : int or str, optional
           The identifier of the background component. If it is None
           then the default background identifier is used.

        See Also
        --------
        delete_background, get_background

        """
        id = self._fix_background_id(id)
        self._backgrounds[id] = bkg
        ids = self.background_ids[:]
        if id not in ids:
            ids.append(id)
        self.background_ids = ids

        # Copy over data from the source to the background
        # if its not present in the background:
        #  - background and grouping
        #  - response information (ONLY THE FIRST TERM)
        #
        # The units (only when a response is present), rate, and
        # plot_fac values are always copied.
        #
        if bkg.grouping is None:
            bkg.grouping = self.grouping
            bkg.grouped = bkg.grouping is not None
        if bkg.quality is None:
            bkg.quality = self.quality

        if bkg.get_response() == (None, None):
            bkg.set_response(*self.get_response())

        if bkg.get_response() != (None, None):
            bkg.units = self.units

        bkg.rate = self.rate
        bkg.plot_fac = self.plot_fac

    def delete_background(self, id=None):
        """Remove the background component.

        If the background component does not exist then the method
        does nothing.

        Parameters
        ----------
        id : int or str, optional
           The identifier of the background component. If it is None
           then the default background identifier is used.

        See Also
        --------
        set_background

        Notes
        -----
        If this call removes the last of the background components
        then the subtracted flag is cleared (if set).

        """

        id = self._fix_background_id(id)
        self._backgrounds.pop(id, None)
        if len(self._backgrounds) == 0:
            self._subtracted = False
        ids = self.background_ids[:]
        if id in ids:
            ids.remove(id)
        self.background_ids = ids

    def get_background_scale(self, bkg_id=1, units='counts',
                             group=True, filter=False):
        """Return the correction factor for the background dataset.

        .. versionchanged:: 4.12.2
           The bkg_id, units, group, and filter parameters have been
           added and the routine no-longer calculates the average
           scaling for all the background components but just for the
           given component.

        Parameters
        ----------
        bkg_id : int or str, optional
           The background component to use (the default is 1).
        units : {'counts', 'rate'}, optional
           The correction is applied to a model defined as counts, the
           default, or a rate. The latter should be used when
           calculating the correction factor for adding the background
           data to the source aperture.
        group : bool, optional
            Should the values be grouped to match the data?
        filter : bool, optional
            Should the values be filtered to match the data?

        Returns
        -------
        scale : None, number, or NumPy array
            The scaling factor to correct the background data onto the
            source data set. If bkg_id is not valid then None is
            returned.

        Notes
        -----
        The correction factor when units is 'counts' is::

            scale_exposure * scale_backscal * scale_areascal / nbkg

        where nbkg is the number of background components and
        scale_x is the source value divided by the background
        value for the field x.

        When units is 'rate' the correction is:

            scale_backscal / nbkg

        and it is currently uncertain whether it should include the
        AREASCAL scaling.

        """

        if units not in ['counts', 'rate']:
            raise ValueError("Invalid units argument: {}".format(units))

        if bkg_id not in self.background_ids:
            return None

        nbkg = len(self.background_ids)

        def correct(obj):
            """Correction factor for the object"""
            ans = 1.0

            # Should we set 0 values to 1 at this stage?
            #
            if obj.backscal is not None:
                ans *= self._check_scale(obj.backscal, group=False)

            if obj.areascal is not None and units == 'counts':
                ans *= self._check_scale(obj.areascal, group=False)

            if obj.exposure is not None and units == 'counts':
                ans *= self._check_scale(obj.exposure, group=False)

            return ans

        src = correct(self)
        bkg = correct(self.get_background(bkg_id))
        scale = src / bkg / nbkg
        return self._check_scale(scale, group=group, filter=filter)

    def _check_scale(self, scale, group=True, filter=False):
        """Ensure the scale value is positive and filtered/grouped.

        Parameters
        ----------
        scale : number or numpy array
            The scale factor.
        group : bool, optional
            Is any grouping applied to the data? This is only
            relevant for an array.
        filter : bool, optional
            Is any filter applied? This is only checked if group
            is True.

        Returns
        -------
        scale : number or numpy array
            Negative values are replaced by 1.0.

        """
        if numpy.isscalar(scale) and scale <= 0.0:
            scale = 1.0
        elif numpy.iterable(scale):
            scale = numpy.asarray(scale, dtype=SherpaFloat)
            if group:
                if filter:
                    scale = self.apply_filter(scale, self._middle)
                else:
                    scale = self.apply_grouping(scale, self._middle)

            scale[scale <= 0.0] = 1.0
        return scale

    def get_backscal(self, group=True, filter=False):
        """Return the background scaling of the PHA data set.

        Return the BACKSCAL setting [BSCAL]_ for the PHA data set.

        Parameters
        ----------
        group : bool, optional
            Should the values be grouped to match the data?
        filter : bool, optional
            Should the values be filtered to match the data?

        Returns
        -------
        backscal : number or ndarray
           The BACKSCAL value, which can be a scalar or a 1D array.

        See Also
        --------
        get_areascal, get_background_scale

        Notes
        -----
        The BACKSCAL value can be defined as the ratio of the area of
        the source (or background) extraction region in image pixels
        to the total number of image pixels. The fact that there is no
        ironclad definition for this quantity does not matter so long
        as the value for a source dataset and its associated
        background dataset are defined in the same manner, because
        only the ratio of source and background BACKSCAL values is
        used. It can be a scalar or an array.

        References
        ----------

        .. [BSCAL] "The OGIP Spectral File Format", Arnaud, K. & George, I.
               http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html

        Examples
        --------

        >>> pha.get_backscal()
        7.8504301607718007e-06

        """
        backscal = self.backscal
        if backscal is not None:
            backscal = self._check_scale(backscal, group, filter)
        return backscal

    def get_areascal(self, group=True, filter=False):
        """Return the fractional area factor of the PHA data set.

        Return the AREASCAL setting [ASCAL]_ for the PHA data set.

        Parameters
        ----------
        group : bool, optional
            Should the values be grouped to match the data?
        filter : bool, optional
            Should the values be filtered to match the data?

        Returns
        -------
        areascal : number or ndarray
           The AREASCAL value, which can be a scalar or a 1D array.

        See Also
        --------
        get_backscal, get_background_scale

        Notes
        -----
        The fractional area scale is normally set to 1, with the ARF used
        to scale the model.

        References
        ----------

        .. [ASCAL] "The OGIP Spectral File Format", Arnaud, K. & George, I.
               http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html

        Examples
        --------

        >>> pha.get_areascal()
        1.0

        """
        areascal = self.areascal
        if areascal is not None:
            areascal = self._check_scale(areascal, group, filter)
        return areascal

    def apply_filter(self, data, groupfunc=numpy.sum):
        """

        Filter the array data, first passing it through apply_grouping()
        (using groupfunc) and then applying the general filters

        """
        if data is None:
            return data

        if len(data) != len(self.counts):
            counts = numpy.zeros(len(self.counts), dtype=SherpaFloat)
            mask = self.get_mask()
            if mask is not None:
                counts[mask] = numpy.asarray(data, dtype=SherpaFloat)
                data = counts
            # else:
            #     raise DataErr('mismatch', "filter", "data array")

        return super().apply_filter(self.apply_grouping(data, groupfunc))

    def apply_grouping(self, data, groupfunc=numpy.sum):
        """

        Apply the data set's grouping scheme to the array data,
        combining the grouped data points with groupfunc, and return
        the grouped array.  If the data set has no associated grouping
        scheme or the data are ungrouped, data is returned unaltered.

        """
        if data is None or not self.grouped:
            return data

        groups = self.grouping
        filter = self.quality_filter
        if filter is None:
            return do_group(data, groups, groupfunc.__name__)

        if len(data) != len(filter) or len(groups) != len(filter):
            raise DataErr('mismatch', "quality filter", "data array")

        filtered_data = numpy.asarray(data)[filter]
        groups = numpy.asarray(groups)[filter]
        grouped_data = do_group(filtered_data, groups, groupfunc.__name__)

        if data is self.channel and groupfunc is self._make_groups:
            return numpy.arange(1, len(grouped_data) + 1, dtype=int)

        return grouped_data

    def ignore_bad(self):
        """Exclude channels marked as bad.

        Ignore any bin in the PHA data set which has a quality value
        that is larger than zero.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set has no quality array.

        See Also
        --------
        ignore : Exclude data from the fit.
        notice : Include data in the fit.

        Notes
        -----
        Bins with a non-zero quality setting are not automatically
        excluded when a data set is created.

        If the data set has been grouped, then calling `ignore_bad`
        will remove any filter applied to the data set. If this
        happens a warning message will be displayed.

        """
        if self.quality is None:
            raise DataErr("noquality", self.name)

        qual_flags = ~numpy.asarray(self.quality, bool)

        if self.grouped and (self.mask is not True):
            self.notice()
            warning('filtering grouped data with quality flags,' +
                    ' previous filters deleted')

        elif not self.grouped:
            # if ungrouped, create/combine with self.mask
            if self.mask is not True:
                self.mask = self.mask & qual_flags
                return
            else:
                self.mask = qual_flags
                return

        # self.quality_filter used for pre-grouping filter
        self.quality_filter = qual_flags

    def _dynamic_group(self, group_func, *args, **kwargs):

        keys = list(kwargs.keys())[:]
        for key in keys:
            if kwargs[key] is None:
                kwargs.pop(key)

        old_filter = self.get_filter(group=False)
        do_notice = numpy.iterable(self.mask)

        self.grouping, self.quality = group_func(*args, **kwargs)
        self.group()
        self._original_groups = False

        if do_notice:
            # self.group() above has cleared the filter if applicable
            # No, that just sets a flag.  So manually clear filter
            # here
            self.ignore()
            for vals in parse_expr(old_filter):
                self.notice(*vals)

        # warning('grouping flags have changed, noticing all bins')

    # Have to move this check here; as formerly written, reference
    # to pygroup functions happened *before* checking groupstatus,
    # in _dynamic_group.  So we did not return the intended error
    # message; rather, a NameError was raised stating that pygroup
    # did not exist in global scope (not too clear to the user).
    #
    # The groupstatus check thus has to be done in *each* of the following
    # group functions.

    # # Dynamic grouping functions now automatically impose the
    # # same grouping conditions on *all* associated background data sets.
    # # CIAO 4.5 bug fix, 05/01/2012
    def group_bins(self, num, tabStops=None):
        """Group into a fixed number of bins.

        Combine the data so that there `num` equal-width bins (or
        groups). The binning scheme is applied to all the channels,
        but any existing filter - created by the `ignore` or `notice`
        set of functions - is re-applied after the data has been
        grouped.

        Parameters
        ----------
        num : int
           The number of bins in the grouped data set. Each bin
           will contain the same number of channels.
        tabStops : array of int or bool, optional
           If set, indicate one or more ranges of channels that should
           not be included in the grouped output. The array should
           match the number of channels in the data set and non-zero or
           `True` means that the channel should be ignored from the
           grouping (use 0 or `False` otherwise).

        See Also
        --------
        group_adapt : Adaptively group to a minimum number of counts.
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_counts : Group into a minimum number of counts per bin.
        group_snr : Group into a minimum signal-to-noise ratio.
        group_width : Group into a fixed bin width.

        Notes
        -----
        Since the bin width is an integer number of channels, it is
        likely that some channels will be "left over". This is even
        more likely when the `tabStops` parameter is set. If this
        happens, a warning message will be displayed to the screen and
        the quality value for these channels will be set to 2.

        """
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpNumBins, len(self.channel), num,
                            tabStops=tabStops)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_bins(num, tabStops=tabStops)

    def group_width(self, val, tabStops=None):
        """Group into a fixed bin width.

        Combine the data so that each bin contains `num` channels.
        The binning scheme is applied to all the channels, but any
        existing filter - created by the `ignore` or `notice` set of
        functions - is re-applied after the data has been grouped.

        Parameters
        ----------
        val : int
           The number of channels to combine into a group.
        tabStops : array of int or bool, optional
           If set, indicate one or more ranges of channels that should
           not be included in the grouped output. The array should
           match the number of channels in the data set and non-zero or
           `True` means that the channel should be ignored from the
           grouping (use 0 or `False` otherwise).

        See Also
        --------
        group_adapt : Adaptively group to a minimum number of counts.
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_bins : Group into a fixed number of bins.
        group_counts : Group into a minimum number of counts per bin.
        group_snr : Group into a minimum signal-to-noise ratio.

        Notes
        -----
        Unless the requested bin width is a factor of the number of
        channels (and no `tabStops` parameter is given), then some
        channels will be "left over". If this happens, a warning
        message will be displayed to the screen and the quality value
        for these channels will be set to 2.

        """
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpBinWidth, len(self.channel), val,
                            tabStops=tabStops)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_width(val, tabStops=tabStops)

    def group_counts(self, num, maxLength=None, tabStops=None):
        """Group into a minimum number of counts per bin.

        Combine the data so that each bin contains `num` or more
        counts. The binning scheme is applied to all the channels, but
        any existing filter - created by the `ignore` or `notice` set
        of functions - is re-applied after the data has been grouped.
        The background is *not* included in this calculation; the
        calculation is done on the raw data even if `subtract` has
        been called on this data set.

        Parameters
        ----------
        num : int
           The number of channels to combine into a group.
        maxLength : int, optional
           The maximum number of channels that can be combined into a
           single group.
        tabStops : array of int or bool, optional
           If set, indicate one or more ranges of channels that should
           not be included in the grouped output. The array should
           match the number of channels in the data set and non-zero or
           `True` means that the channel should be ignored from the
           grouping (use 0 or `False` otherwise).

        See Also
        --------
        group_adapt : Adaptively group to a minimum number of counts.
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_bins : Group into a fixed number of bins.
        group_snr : Group into a minimum signal-to-noise ratio.
        group_width : Group into a fixed bin width.

        Notes
        -----
        If channels can not be placed into a "valid" group, then a
        warning message will be displayed to the screen and the
        quality value for these channels will be set to 2.

        """
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpNumCounts, self.counts, num,
                            maxLength=maxLength, tabStops=tabStops)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_counts(num, maxLength=maxLength, tabStops=tabStops)

    # DOC-TODO: see discussion in astro.ui.utils regarding errorCol
    def group_snr(self, snr, maxLength=None, tabStops=None, errorCol=None):
        """Group into a minimum signal-to-noise ratio.

        Combine the data so that each bin has a signal-to-noise ratio
        which exceeds `snr`. The binning scheme is applied to all the
        channels, but any existing filter - created by the `ignore` or
        `notice` set of functions - is re-applied after the data has
        been grouped.  The background is *not* included in this
        calculation; the calculation is done on the raw data even if
        `subtract` has been called on this data set.

        Parameters
        ----------
        snr : number
           The minimum signal-to-noise ratio that must be exceeded
           to form a group of channels.
        maxLength : int, optional
           The maximum number of channels that can be combined into a
           single group.
        tabStops : array of int or bool, optional
           If set, indicate one or more ranges of channels that should
           not be included in the grouped output. The array should
           match the number of channels in the data set and non-zero or
           `True` means that the channel should be ignored from the
           grouping (use 0 or `False` otherwise).
        errorCol : array of num, optional
           If set, the error to use for each channel when calculating
           the signal-to-noise ratio. If not given then Poisson
           statistics is assumed. A warning is displayed for each
           zero-valued error estimate.

        See Also
        --------
        group_adapt : Adaptively group to a minimum number of counts.
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_bins : Group into a fixed number of bins.
        group_counts : Group into a minimum number of counts per bin.
        group_width : Group into a fixed bin width.

        Notes
        -----
        If channels can not be placed into a "valid" group, then a
        warning message will be displayed to the screen and the
        quality value for these channels will be set to 2.

        """
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpSnr, self.counts, snr,
                            maxLength=maxLength, tabStops=tabStops,
                            errorCol=errorCol)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_snr(snr, maxLength=maxLength, tabStops=tabStops,
                          errorCol=errorCol)

    def group_adapt(self, minimum, maxLength=None, tabStops=None):
        """Adaptively group to a minimum number of counts.

        Combine the data so that each bin contains `num` or more
        counts. The difference to `group_counts` is that this
        algorithm starts with the bins with the largest signal, in
        order to avoid over-grouping bright features, rather than at
        the first channel of the data. The adaptive nature means that
        low-count regions between bright features may not end up in
        groups with the minimum number of counts.  The binning scheme
        is applied to all the channels, but any existing filter -
        created by the `ignore` or `notice` set of functions - is
        re-applied after the data has been grouped.

        Parameters
        ----------
        minimum : int
           The number of channels to combine into a group.
        maxLength : int, optional
           The maximum number of channels that can be combined into a
           single group.
        tabStops : array of int or bool, optional
           If set, indicate one or more ranges of channels that should
           not be included in the grouped output. The array should
           match the number of channels in the data set and non-zero or
           `True` means that the channel should be ignored from the
           grouping (use 0 or `False` otherwise).

        See Also
        --------
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_bins : Group into a fixed number of bins.
        group_counts : Group into a minimum number of counts per bin.
        group_snr : Group into a minimum signal-to-noise ratio.
        group_width : Group into a fixed bin width.

        Notes
        -----
        If channels can not be placed into a "valid" group, then a
        warning message will be displayed to the screen and the
        quality value for these channels will be set to 2.

        """
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpAdaptive, self.counts, minimum,
                            maxLength=maxLength, tabStops=tabStops)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_adapt(minimum, maxLength=maxLength,
                            tabStops=tabStops)

    # DOC-TODO: see discussion in astro.ui.utils regarding errorCol
    def group_adapt_snr(self, minimum, maxLength=None, tabStops=None,
                        errorCol=None):
        """Adaptively group to a minimum signal-to-noise ratio.

        Combine the data so that each bin has a signal-to-noise ratio
        which exceeds `minimum`. The difference to `group_snr` is that
        this algorithm starts with the bins with the largest signal,
        in order to avoid over-grouping bright features, rather than
        at the first channel of the data. The adaptive nature means
        that low-count regions between bright features may not end up
        in groups with the minimum number of counts.  The binning
        scheme is applied to all the channels, but any existing filter
        - created by the `ignore` or `notice` set of functions - is
        re-applied after the data has been grouped.

        Parameters
        ----------
        minimum : number
           The minimum signal-to-noise ratio that must be exceeded
           to form a group of channels.
        maxLength : int, optional
           The maximum number of channels that can be combined into a
           single group.
        tabStops : array of int or bool, optional
           If set, indicate one or more ranges of channels that should
           not be included in the grouped output. The array should
           match the number of channels in the data set and non-zero or
           `True` means that the channel should be ignored from the
           grouping (use 0 or `False` otherwise).
        errorCol : array of num, optional
           If set, the error to use for each channel when calculating
           the signal-to-noise ratio. If not given then Poisson
           statistics is assumed. A warning is displayed for each
           zero-valued error estimate.

        See Also
        --------
        group_adapt : Adaptively group to a minimum number of counts.
        group_bins : Group into a fixed number of bins.
        group_counts : Group into a minimum number of counts per bin.
        group_snr : Group into a minimum signal-to-noise ratio.
        group_width : Group into a fixed bin width.

        Notes
        -----
        If channels can not be placed into a "valid" group, then a
        warning message will be displayed to the screen and the
        quality value for these channels will be set to 2.

        """
        if not groupstatus:
            raise ImportErr('importfailed', 'group', 'dynamic grouping')
        self._dynamic_group(pygroup.grpAdaptiveSnr, self.counts, minimum,
                            maxLength=maxLength, tabStops=tabStops,
                            errorCol=errorCol)
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_adapt_snr(minimum, maxLength=maxLength,
                                tabStops=tabStops, errorCol=errorCol)

    def eval_model(self, modelfunc):
        return modelfunc(*self.get_indep(filter=False))

    def eval_model_to_fit(self, modelfunc):
        return self.apply_filter(modelfunc(*self.get_indep(filter=True)))

    def sum_background_data(self,
                            get_bdata_func=(lambda key, bkg: bkg.counts)):
        """Sum up data, applying the background correction value.

        Parameters
        ----------
        get_bdata_func : function, optional
            What data should be used for each background dataset. The
            function takes the background identifier and background
            DataPHA object and returns the data to use. The default is
            to use the counts array of the background dataset.

        Returns
        -------
        value : scalar or NumPy array
            The sum of the data, including any area, background, and
            exposure-time corrections.

        Notes
        -----
        For each associated background, the data is retrieved (via
        the get_bdata_func parameter), and then

          - divided by its BACKSCAL value (if set)
          - divided by its AREASCAL value (if set)
          - divided by its exposure time (if set)

        The individual background components are then summed together,
        and then multiplied by the source BACKSCAL (if set),
        multiplied by the source AREASCAL (if set), and multiplied
        by the source exposure time (if set). The final step is
        to divide by the number of background files used.

        Example
        -------

        Calculate the background counts, per channel, scaled to match
        the source:

        >>> bcounts = src.sum_background_data()

        Calculate the scaling factor that you need to multiply the
        background data to match the source data. In this case the
        background data has been replaced by the value 1 (rather than
        the per-channel values used with the default argument):

        >>> bscale = src.sum_background_data(lambda k, d: 1)

        """

        bdata_list = []

        for key in self.background_ids:
            bkg = self.get_background(key)
            bdata = get_bdata_func(key, bkg)

            backscal = bkg.backscal
            if backscal is not None:
                backscal = self._check_scale(backscal, group=False)
                bdata = bdata / backscal

            areascal = bkg.get_areascal(group=False)
            if areascal is not None:
                bdata = bdata / areascal

            if bkg.exposure is not None:
                bdata = bdata / bkg.exposure

            bdata_list.append(bdata)

        nbkg = len(bdata_list)
        if nbkg == 0:
            # do not have a good id to use for the error message
            raise DataErr('nobkg', self.name)

        if nbkg == 1:
            bkgsum = bdata_list[0]
        else:
            bkgsum = sum(bdata_list)

        backscal = self.backscal
        if backscal is not None:
            backscal = self._check_scale(backscal, group=False)
            bkgsum = backscal * bkgsum

        areascal = self.areascal
        if areascal is not None:
            areascal = self._check_scale(areascal, group=False)
            bkgsum = areascal * bkgsum

        if self.exposure is not None:
            bkgsum = self.exposure * bkgsum

        return bkgsum / SherpaFloat(nbkg)

    def get_dep(self, filter=False):
        # FIXME: Aneta says we need to group *before* subtracting, but that
        # won't work (I think) when backscal is an array
        # if not self.subtracted:
        #     return self.counts
        # return self.counts - self.sum_background_data()

        dep = self.counts
        filter = bool_cast(filter)

        # The area scaling is not applied to the data, since it
        # should be being applied to the model via the *PHA
        # instrument model. Note however that the background
        # contribution does include the source AREASCAL value
        # (in the same way that the source BACKSCAL value
        # is used).
        #
        if self.subtracted:
            bkg = self.sum_background_data()
            if len(dep) != len(bkg):
                raise DataErr("subtractlength")
            dep = dep - bkg

        if filter:
            dep = self.apply_filter(dep)
        return dep

    def set_dep(self, val):
        # QUS: should this "invert" the areascaling to val
        #      to get the stored values?
        #
        #      Otherwise, when areascal /= 1
        #            y1 = d.get_dep()
        #            d.set_dep(y1)
        #            y2 = d.get_dep()
        #            y1 != y2
        #
        # Or perhaps it removes the areascal value in this case?
        # We already have this split in the API when background data
        # is available and is subtracted.
        #
        if numpy.iterable(val):
            dep = numpy.asarray(val, SherpaFloat)
        else:
            val = SherpaFloat(val)
            dep = numpy.array([val] * len(self.get_indep()[0]))
        setattr(self, 'counts', dep)

    def get_staterror(self, filter=False, staterrfunc=None):
        """Return the statistical error.

        The staterror column is used if defined, otherwise the
        function provided by the staterrfunc argument is used to
        calculate the values.

        Parameters
        ----------
        filter : bool, optional
            Should the channel filter be applied to the return values?
        staterrfunc : function reference, optional
            The function to use to calculate the errors if the
            staterror field is None. The function takes one argument,
            the counts (after grouping and filtering), and returns an
            array of values which represents the one-sigma error for each
            element of the input array. This argument is designed to
            work with implementations of the sherpa.stats.Stat.calc_staterror
            method.

        Returns
        -------
        staterror : array or None
            The statistical error. It will be grouped and,
            if filter=True, filtered. The contribution from any
            associated background components will be included if
            the background-subtraction flag is set.

        Notes
        -----
        There is no scaling by the AREASCAL setting, but background
        values are scaled by their AREASCAL settings. It is not at all
        obvious that the current code is doing the right thing, or that
        this is the right approach.

        Examples
        --------

        >>> dy = dset.get_staterror()

        Ensure that there is no pre-defined statistical-error column
        and then use the Chi2DataVar statistic to calculate the errors:

        >>> stat = sherpa.stats.Chi2DataVar()
        >>> dset.set_staterror(None)
        >>> dy = dset.get_staterror(staterrfunc=stat.calc_staterror)

        """
        staterr = self.staterror

        filter = bool_cast(filter)
        if filter:
            staterr = self.apply_filter(staterr, self._sum_sq)
        else:
            staterr = self.apply_grouping(staterr, self._sum_sq)

        # The source AREASCAL is not applied here, but the
        # background term is.
        #
        if (staterr is None) and (staterrfunc is not None):
            cnts = self.counts

            if filter:
                cnts = self.apply_filter(cnts)
            else:
                cnts = self.apply_grouping(cnts)

            staterr = staterrfunc(cnts)

            # Need to apply the area scaling to the calculated
            # errors. Grouping and filtering complicate this; is
            # _middle the best choice here?
            #
            """
            area = self.areascal
            if staterr is not None and area is not None:
                if numpy.isscalar(area):
                    area = numpy.zeros(self.channel.size) + area

                # TODO: replace with _check_scale?
                if filter:
                    area = self.apply_filter(area, self._middle)
                else:
                    area = self.apply_grouping(area, self._middle)
                staterr = staterr / area
            """

        if (staterr is not None) and self.subtracted:
            bkg_staterr_list = []

            # for bkg in self._backgrounds.values():
            for key in self.background_ids:
                bkg = self.get_background(key)
                berr = bkg.staterror
                if filter:
                    berr = self.apply_filter(berr, self._sum_sq)
                else:
                    berr = self.apply_grouping(berr, self._sum_sq)

                if (berr is None) and (staterrfunc is not None):
                    bkg_cnts = bkg.counts
                    if filter:
                        bkg_cnts = self.apply_filter(bkg_cnts)
                    else:
                        bkg_cnts = self.apply_grouping(bkg_cnts)

                    # TODO: shouldn't the following logic be somewhere
                    #       else more general?
                    if hasattr(staterrfunc, '__name__') and \
                       staterrfunc.__name__ == 'calc_chi2datavar_errors' and \
                       0.0 in bkg_cnts:
                        mask = (numpy.asarray(bkg_cnts) != 0.0)
                        berr = numpy.zeros(len(bkg_cnts))
                        berr[mask] = staterrfunc(bkg_cnts[mask])
                    else:
                        berr = staterrfunc(bkg_cnts)

                # FIXME: handle this
                # assert (berr is not None)

                # This case appears when the source dataset has an error
                # column and at least one of the background(s) do not.
                # Because the staterr is not None and staterrfunc is, I think
                # we should return None.  This way the user knows to call with
                # staterrfunc next time.
                if berr is None:
                    return None

                bksl = bkg.backscal
                if bksl is not None:
                    bksl = self._check_scale(bksl, filter=filter)
                    berr = berr / bksl

                # Need to apply filter/grouping of the source dataset
                # to the background areascal, so can not just say
                #   area = bkg.get_areascal(filter=filter)
                #
                area = bkg.areascal
                if area is not None:
                    area = self._check_scale(area, filter=filter)
                    berr = berr / area

                if bkg.exposure is not None:
                    berr = berr / bkg.exposure

                berr = berr * berr
                bkg_staterr_list.append(berr)

            nbkg = len(bkg_staterr_list)
            assert (nbkg > 0)
            if nbkg == 1:
                bkgsum = bkg_staterr_list[0]
            else:
                bkgsum = sum(bkg_staterr_list)

            bscal = self.backscal
            if bscal is not None:
                bscal = self._check_scale(bscal, filter=filter)
                bkgsum = (bscal * bscal) * bkgsum

            # Correct the background counts by the source AREASCAL
            # setting. Is this correct?
            ascal = self.areascal
            if ascal is not None:
                ascal = self._check_scale(ascal, filter=filter)
                bkgsum = (ascal * ascal) * bkgsum

            if self.exposure is not None:
                bkgsum = (self.exposure * self.exposure) * bkgsum

            nbkg = SherpaFloat(nbkg)

            if staterr is not None:
                staterr = staterr * staterr + bkgsum / (nbkg * nbkg)
                staterr = numpy.sqrt(staterr)

        return staterr

    def get_syserror(self, filter=False):
        """Return any systematic error.

        Parameters
        ----------
        filter : bool, optional
            Should the channel filter be applied to the return values?

        Returns
        -------
        syserror : array or None
            The systematic error, if set. It will be grouped and,
            if filter=True, filtered.

        Notes
        -----
        There is no scaling by the AREASCAL setting.
        """
        syserr = self.syserror
        filter = bool_cast(filter)
        if filter:
            syserr = self.apply_filter(syserr, self._sum_sq)
        else:
            syserr = self.apply_grouping(syserr, self._sum_sq)
        return syserr

    def get_x(self, filter=False, response_id=None):
        # We want the full channel grid with no grouping.
        #
        return self._from_channel(self.channel, group=False, response_id=response_id)

    def get_xlabel(self):
        xlabel = self.units.capitalize()
        if self.units == 'energy':
            xlabel += ' (keV)'
        elif self.units == 'wavelength':
            xlabel += ' (Angstrom)'
        # elif self.units == 'channel' and self.grouped:
        #     xlabel = 'Group Number'
        return xlabel

    def _set_initial_quantity(self):
        arf, rmf = self.get_response()

        # Change analysis if ARFs equal or of higher resolution to
        # allow for high-res model evaluation.
        if arf is not None and rmf is None:
            if len(arf.energ_lo) == len(self.channel):
                self.units = 'energy'

        # Only change analysis if RMF matches the parent PHA dataset.
        if rmf is not None:
            if len(self.channel) != len(rmf.e_min):
                raise DataErr("incompatibleresp", rmf.name, self.name)
            self.units = 'energy'

    def _fix_y_units(self, val, filter=False, response_id=None):
        """Rescale the data to match the 'y' axis."""

        if val is None:
            return val

        filter = bool_cast(filter)
        # make a copy of data for units manipulation
        val = numpy.array(val, dtype=SherpaFloat)

        if self.rate and self.exposure is not None:
            val /= self.exposure

        # TODO: It is not clear if the areascal should always be applied,
        #       or only if self.rate is set (since it is being considered
        #       a "correction" to the exposure time, but don't we want
        #       to apply it in plots even if the Y axis is in counts?)
        #
        if self.areascal is not None:
            areascal = self._check_scale(self.areascal, filter=filter)
            val /= areascal

        if self.grouped or self.rate:

            if self.units != 'channel':
                elo, ehi = self._get_ebins(response_id, group=False)
            else:
                elo, ehi = (self.channel, self.channel + 1.)

            if filter:
                # If we apply a filter, make sure that
                # ebins are ungrouped before applying
                # the filter.
                elo = self.apply_filter(elo, self._min)
                ehi = self.apply_filter(ehi, self._max)
            elif self.grouped:
                elo = self.apply_grouping(elo, self._min)
                ehi = self.apply_grouping(ehi, self._max)

            if self.units == 'energy':
                ebin = ehi - elo
            elif self.units == 'wavelength':
                ebin = hc / elo - hc / ehi
            elif self.units == 'channel':
                ebin = ehi - elo
            else:
                raise DataErr("bad", "quantity", self.units)

            val /= numpy.abs(ebin)

        # The final step is to multiply by the X axis self.plot_fac
        # times.
        if self.plot_fac <= 0:
            return val

        scale = self.apply_filter(self.get_x(response_id=response_id),
                                  self._middle)
        for ii in range(self.plot_fac):
            val *= scale

        return val

    def get_y(self, filter=False, yfunc=None, response_id=None, use_evaluation_space=False):
        vallist = Data.get_y(self, yfunc=yfunc)
        filter = bool_cast(filter)

        if not isinstance(vallist, tuple):
            vallist = (vallist,)

        newvallist = []

        for val in vallist:
            if filter:
                val = self.apply_filter(val)
            else:
                val = self.apply_grouping(val)
            val = self._fix_y_units(val, filter, response_id)
            newvallist.append(val)

        if len(vallist) == 1:
            vallist = newvallist[0]
        else:
            vallist = tuple(newvallist)

        return vallist

    def get_yerr(self, filter=False, staterrfunc=None, response_id=None):
        filter = bool_cast(filter)
        err = self.get_error(filter, staterrfunc)
        return self._fix_y_units(err, filter, response_id)

    def get_xerr(self, filter=False, response_id=None):
        elo, ehi = self._get_ebins(response_id=response_id)
        filter = bool_cast(filter)
        if filter:
            # If we apply a filter, make sure that
            # ebins are ungrouped before applying
            # the filter.
            elo, ehi = self._get_ebins(response_id, group=False)
            elo = self.apply_filter(elo, self._min)
            ehi = self.apply_filter(ehi, self._max)

        return ehi - elo

    def get_ylabel(self):
        ylabel = 'Counts'

        if self.rate and self.exposure:
            ylabel += '/sec'

        if self.rate or self.grouped:
            if self.units == 'energy':
                ylabel += '/keV'
            elif self.units == 'wavelength':
                ylabel += '/Angstrom'
            elif self.units == 'channel':
                ylabel += '/channel'

        if self.plot_fac:
            from sherpa.plot import backend
            latex = backend.get_latex_for_string(
                '^{}'.format(self.plot_fac))
            ylabel += ' X {}{}'.format(self.units.capitalize(), latex)

        return ylabel

    @staticmethod
    # Dummy function to tell apply_grouping to construct
    # an array of groups.
    def _make_groups(array):
        pass

    @staticmethod
    def _middle(array):
        array = numpy.asarray(array)
        return (array.min() + array.max()) / 2.0

    @staticmethod
    def _min(array):
        array = numpy.asarray(array)
        return array.min()

    @staticmethod
    def _max(array):
        array = numpy.asarray(array)
        return array.max()

    @staticmethod
    def _sum_sq(array):
        return numpy.sqrt(numpy.sum(array * array))

    def get_noticed_channels(self):
        """Return the noticed channels.

        Returns
        -------
        channels : ndarray
            The noticed channels (this is independent of the
            analysis setting).

        """
        chans = self.channel
        mask = self.get_mask()
        if mask is None:
            return chans

        # This is added to address issue #361
        #
        # If there is a quality filter then the mask may be
        # smaller than the chans array. It is not clear if this
        # is the best location for this. If it is, then are there
        # other locations where this logic is needed?
        #
        if self.quality_filter is not None and \
           self.quality_filter.size != mask.size:
            chans = chans[self.quality_filter]

        return chans[mask]

    def get_mask(self):
        """Returns the (ungrouped) mask.

        Returns
        -------
        mask : ndarray or None
            The mask, in channels, or None.

        """
        groups = self.grouping
        if self.mask is False:
            return None

        if self.mask is True or not self.grouped:
            if self.quality_filter is not None:
                return self.quality_filter
            elif numpy.iterable(self.mask):
                return self.mask
            return None

        if self.quality_filter is not None:
            groups = groups[self.quality_filter]
        return expand_grouped_mask(self.mask, groups)

    def get_noticed_expr(self):
        chans = self.get_noticed_channels()
        if self.mask is False or len(chans) == 0:
            return 'No noticed channels'
        return create_expr(chans, format='%i')

    def get_filter(self, group=True, format='%.12f', delim=':'):
        """Return the data filter as a string.

        For grouped data, or when the analysis setting is not
        channel, filter values refer to the center of the
        channel or group.

        Parameters
        ----------
        group : bool, optional
            Should the filter reflect the grouped data?
        format : str, optional
            The formatting of the numeric values (this is
            ignored for channel units, as a format of "%i"
            is used).
        delim : str, optional
            The string used to mark the low-to-high range.

        Examples
        --------
        For a Chandra non-grating dataset which has been grouped:

        >>> pha.set_analysis('energy')
        >>> pha.notice(0.5, 7)
        >>> pha.get_filter(format=%.4f')
        ''0.5183:8.2198''
        >>> pha.set_analysis('channel')
        >>> pha.get_filter()
        '36:563'

        The default is to show the data range for the grouped
        dataset, which uses the center of each group. If
        the grouping is turned off then the center of the
        start and ending channel of each group is used
        (and so show a larger data range):

        >>> pha.get_filter(format=%.4f')
        '0.5183:8.2198'
        >>> pha.get_filter(group=False, format=%.4f')
        '0.4745:9.8623'

        """
        if self.mask is False:
            return 'No noticed bins'

        if numpy.iterable(self.mask):
            mask = self.mask
        else:
            mask = None

        if group:
            # grouped noticed channels
            #
            x = self.apply_filter(self.channel, self._make_groups)

        else:
            # ungrouped noticed channels
            x = self.get_noticed_channels()

            # We need the "ungrouped" mask array. Need to check
            # issue #361 since get_noticed_channels notes an
            # issue that may be relevant here (so far this
            # doesn't seem to be the case).
            #
            mask = self.get_mask()

            # Safety check for users. Warn, but continue.
            #
            if mask is not None and mask.sum() != x.size:
                warning("There is a mis-match in the ungrouped mask " +
                        "and data ({} vs {})".format(mask.sum(), x.size))

        # Convert channels to appropriate quantity if necessary
        x = self._from_channel(x, group=group)

        if mask is None:
            mask = numpy.ones(len(x), dtype=bool)

        # Ensure the data is in ascending order for create_expr.
        #
        if self.units == 'wavelength':
            x = x[::-1]
            mask = mask[::-1]

        if self.units == 'channel':
            format = '%i'

        return create_expr(x, mask=mask, format=format, delim=delim)

    def get_filter_expr(self):
        return (self.get_filter(format='%.4f', delim='-') +
                ' ' + self.get_xlabel())

    def notice_response(self, notice_resp=True, noticed_chans=None):
        notice_resp = bool_cast(notice_resp)

        if notice_resp and noticed_chans is None:
            noticed_chans = self.get_noticed_channels()

        for id in self.response_ids:
            arf, rmf = self.get_response(id)
            _notice_resp(noticed_chans, arf, rmf)

    def notice(self, lo=None, hi=None, ignore=False, bkg_id=None):
        """Notice or ignore the given range.

        Parameters
        ----------
        lo, hi : number or None, optional
            The range to change. A value of None means the minimum or
            maximum permitted value. The units of lo and hi are set by
            the units field.
        ignore : bool, optional
            Set to True if the range should be ignored. The default is
            to notice the range.
        bkg_id : int or sequence of int or None, optional
            If not None then apply the filter to the given background
            dataset or datasets, otherwise change the object and all
            its background datasets.

        See Also
        --------
        get_filter, get_filter_expr, get_mask

        Notes
        -----
        If no channels have been ignored then a call to `notice` with
        `ignore=False` will select just the `lo` to `hi` range, and
        exclude any channels outside this range. If there has been a
        filter applied then the range `lo` to `hi` will be added to the
        range of noticed data (when `ignore=False`).

        So, for an ungrouped PHA file with 1024 channels:

        >>> pha.units = 'channel'
        >>> pha.get_filter()
        '1:1024'
        >>> pha.notice(20, 200)
        >>> pha.get_filter()
        '20:200'
        >>> pha.notice(300, 500)
        '20:200,300:500'

        """

        ignore = bool_cast(ignore)

        # This condition is checked for in the _data_space.filter call
        # at the end of the method, but it is easier to enforce it
        # here so we do not need to worry about possible type errors
        # when comparing string and number values.
        #
        for val, label in zip([lo, hi], ['lower', 'upper']):
            if isinstance(val, str):
                # match the error seen from other data classes here
                raise DataErr('typecheck', f'{label} bound')

        # If any background IDs are actually given, then impose
        # the filter on those backgrounds *only*, and return.  Do
        # *not* impose filter on data itself.  (Revision possibly
        # this should be done in high-level UI?)  SMD 10/25/12

        filter_background_only = False
        if bkg_id is not None:
            if not numpy.iterable(bkg_id):
                bkg_id = [bkg_id]
            filter_background_only = True
        else:
            bkg_id = self.background_ids

        # Automatically impose data's filter on background data sets.
        # Units must agree for this to be meaningful, so temporarily
        # make data and background units match. SMD 10/25/12
        for bid in bkg_id:
            bkg = self.get_background(bid)
            old_bkg_units = bkg.units
            try:
                bkg.units = self.units
                # If the background is all ignored then bkg.notice will
                # do nothing (other than display an INFO message).
                #
                bkg.notice(lo, hi, ignore)
            finally:
                bkg.units = old_bkg_units

        # If we're only supposed to filter backgrounds, return
        if filter_background_only:
            return

        # Go on if we are also supposed to filter the source data
        if lo is None and hi is None:
            self.quality_filter = None
            self.notice_response(False)

        # When convert_limit is called, any limit that is outside the
        # supported range (e.g. lo = 0.01 when the response is 0.1 to
        # 11 keV) will get clamped to the minumum or maximum value.
        # This is okay when lo is smaller and hi is larger than the
        # supported values but is a problem when they are both
        # smaller or larger than the range, since the result is then
        # that the first or last group would be filtered. It is
        # therefore necessary to check for this condition.
        #
        # This would be easier if _to_channel were to identify that
        # the range had been (or should be) clamped, but that is a
        # larger change than I want to implement, so instead we repeat
        # some of the checks that the _xxx_to_channel methods make to
        # see what the minimum and maximum allowed values are.
        #
        try:
            elo, ehi = self._get_ebins(group=False)
        except DataErr as de:
            info(f"Skipping dataset {self.name}: {de}")
            return

        if self.units == 'energy':
            # It's not guaranteed that channels are in increasing energy
            # but we can certainly assume that they are monotonic.
            # Thus, the min of the first and last entry is the minimum
            # energy of the elo array
            umin = min(elo[0], elo[-1])
            umax = max(ehi[0], ehi[-1])
        elif self.units == 'wavelength':
            umin = min(hc / ehi[[0, -1]])
            umax = max(hc / elo[[0, -1]])
        else:
            # assume channel units
            umin = min(self.channel[0], self.channel[-1])
            umax = max(self.channel[0], self.channel[-1])

        assert umin < umax, (self.units, umin, umax)

        # Finally we can check that the request is not
        # outside the valid range "on the same side".
        #
        if lo is not None and hi is not None:
            assert lo <= hi, (lo, hi)
            if hi <= umin:
                return
            if lo >= umax:
                return

        # We do not want a "all data are masked out" error to cause
        # this to fail; it should just do nothing (as trying to set
        # a noticed range to include masked-out ranges would also
        # be ignored). Note that originally this try/except caught
        # errors like "RMF has no ebounds data", but this is now
        # caught earlier in the call to _get_ebins, so it is possible
        # that the check is no-longer needed.
        #
        # Convert to "group number" (which, for ungrouped data,
        # is just channel number).
        #
        def convert_limit(x):
            """Convert the limit to 'group number'.

            x must not be None or a string. The return value is None if
            the value is invalid (and so should be skipped).
            """

            try:
                return self._to_channel(x)
            except DataErr as de:
                info(f"Skipping dataset {self.name}: {de}")
                return None

        if lo is not None:
            lo = convert_limit(lo)
            if lo is None:
                return

        if hi is not None:
            hi = convert_limit(hi)
            if hi is None:
                return

        # At this point lo and hi are in "group number".
        #
        # When do we need to flip the limits? It is not clear whether
        # we always need lo <= hi.
        #
        if ((self.units == "wavelength" and
             elo[0] < elo[-1] and ehi[0] < ehi[-1]) or
            (self.units == "energy" and
             elo[0] > elo[-1] and ehi[0] > ehi[-1])):
            lo, hi = hi, lo

        # Don't use the middle of the channel anymore as the
        # grouping function.  That was just plain dumb.
        # So just get back an array of groups 1-N, if grouped
        # DATA-NOTE: need to clean this up.
        #
        groups = self.apply_grouping(self.channel,
                                     self._make_groups)
        self._data_space.filter.notice((lo,), (hi,),
                                       (groups,), ignore)

    def to_guess(self):
        elo, ehi = self._get_ebins(group=False)
        elo = self.apply_filter(elo, self._min)
        ehi = self.apply_filter(ehi, self._max)
        if self.units == "wavelength":
            lo = hc / ehi
            hi = hc / elo
            elo = lo
            ehi = hi
        cnt = self.get_dep(True)
        arf = self.get_specresp(filter=True)

        y = cnt / (ehi - elo)
        if self.exposure is not None:
            y /= self.exposure   # photons/keV/sec or photons/Ang/sec
        # y = cnt/arf/self.exposure
        if arf is not None:
            y /= arf  # photons/keV/cm^2/sec or photons/Ang/cm^2/sec
        return (y, elo, ehi)

    def to_fit(self, staterrfunc=None):
        return (self.get_dep(True),
                self.get_staterror(True, staterrfunc),
                self.get_syserror(True))

    def to_plot(self, yfunc=None, staterrfunc=None, response_id=None):
        return (self.apply_filter(self.get_x(response_id=response_id),
                                  self._middle),
                self.get_y(True, yfunc, response_id=response_id),
                self.get_yerr(True, staterrfunc, response_id=response_id),
                self.get_xerr(True, response_id=response_id),
                self.get_xlabel(),
                self.get_ylabel())

    def group(self):
        """Group the data according to the data set's grouping scheme.

        This sets the grouping flag which means that the value of the
        grouping attribute will be used when accessing data values. This
        can be called even if the grouping attribute is empty.

        See Also
        --------
        ungroup
        """
        self.grouped = True

    def ungroup(self):
        """Remove any data grouping.

        This un-sets the grouping flag which means that the grouping
        attribute will not be used when accessing data values.

        See Also
        --------
        group
        """
        self.grouped = False

    def subtract(self):
        """Subtract the background data.

        See Also
        --------
        unsubtract
        """
        self.subtracted = True

    def unsubtract(self):
        """Remove background subtraction.

        See Also
        --------
        subtract
        """
        self.subtracted = False


class DataIMG(Data2D):
    "Image data set, including functions for coordinate transformations"
    _fields = Data2D._fields + ("sky", "eqpos", "coord", "header")

    def _get_coord(self):
        return self._coord

    def _set_coord(self, val):
        coord = str(val).strip().lower()

        if coord in ('logical', 'image'):
            coord = 'logical'

        elif coord in ('physical',):
            self._check_physical_transform()
            coord = 'physical'

        elif coord in ('world', 'wcs'):
            self._check_world_transform()
            coord = 'world'

        else:
            raise DataErr('bad', 'coordinates', val)

        self._coord = coord

    # You should use set_coord rather than changing coord directly,
    # otherwise constraints set in set_coord are not run. This is
    # probably an error in set_coord (i.e. this logic should be
    # moved into _set_coord).
    #
    coord = property(_get_coord, _set_coord,
                     doc='Coordinate system of independent axes')

    def __init__(self, name, x0, x1, y, shape=None, staterror=None,
                 syserror=None, sky=None, eqpos=None, coord='logical',
                 header=None):
        self.sky = sky
        self.eqpos = eqpos
        self.coord = coord
        self.header = {} if header is None else header
        self._region = None
        Data2D.__init__(self, name, x0, x1, y, shape, staterror, syserror)

    def __str__(self):
        # Print the metadata first
        old = self._fields
        ss = old
        try:
            self._fields = tuple(filter((lambda x: x != 'header'),
                                        self._fields))
            ss = Data.__str__(self)
        finally:
            self._fields = old
        return ss

    def _repr_html_(self):
        """Return a HTML (string) representation of the data
        """
        return html_img(self)

    def __getstate__(self):
        state = self.__dict__.copy()
        # Function pointers to methods of the class
        # (of type 'instancemethod') are NOT picklable
        # remove them and restore later with a coord init
        # del state['_get_logical']
        # del state['_get_physical']
        # del state['_get_world']

        # PyRegion objects (of type 'extension') are NOT picklable, yet.
        # preserve the region string and restore later with constructor
        state['_region'] = state['_region'].__str__()
        return state

    def __setstate__(self, state):
        # Populate the function pointers we deleted at pickle time with
        # no-ops.
        # self.__dict__['_get_logical']=(lambda : None)
        # self.__dict__['_get_physical']=(lambda : None)
        # self.__dict__['_get_world']=(lambda : None)

        if 'header' not in state:
            self.header = {}

        self.__dict__.update(state)

        # _set_coord will correctly define the _get_* WCS function pointers.
        self._set_coord(state['_coord'])
        if regstatus:
            self._region = Region(self._region)
        else:
            # An ImportErr could be raised rather than display a
            # warnng, but that would make it harder for the user
            # to extract useful data (e.g. in the case of triggering
            # this when loading a pickled file).
            #
            if self._region is not None and self._region != '':
                warning("Unable to restore region={} as region module is not avaialable.".format(self._region))

            self._region = None

    def _check_physical_transform(self):
        if self.sky is None:
            raise DataErr('nocoord', self.name, 'physical')

    def _check_world_transform(self):
        if self.eqpos is None:
            raise DataErr('nocoord', self.name, 'world')

    def _logical_to_physical(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_physical_transform()

        # logical -> physical
        x0, x1 = self.sky.apply(x0, x1)

        return (x0, x1)

    def _logical_to_world(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_world_transform()

        # logical -> physical
        if self.sky is not None:
            x0, x1 = self.sky.apply(x0, x1)

        # physical -> world
        x0, x1 = self.eqpos.apply(x0, x1)

        return (x0, x1)

    def _physical_to_logical(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_physical_transform()

        # physical -> logical
        x0, x1 = self.sky.invert(x0, x1)

        return (x0, x1)

    def _physical_to_world(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_world_transform()

        # physical -> world
        x0, x1 = self.eqpos.apply(x0, x1)

        return (x0, x1)

    def _world_to_logical(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_world_transform()

        # world -> physical
        x0, x1 = self.eqpos.invert(x0, x1)

        # physical -> logical
        if self.sky is not None:
            x0, x1 = self.sky.invert(x0, x1)

        return (x0, x1)

    def _world_to_physical(self, x0=None, x1=None):
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_world_transform()

        # world -> physical
        x0, x1 = self.eqpos.invert(x0, x1)

        return (x0, x1)

    def get_logical(self):
        coord = self.coord
        x0, x1 = self.get_indep()
        if coord != 'logical':
            x0 = x0.copy()
            x1 = x1.copy()
            x0, x1 = getattr(self, '_' + coord + '_to_logical')(x0, x1)
        return (x0, x1)

    def get_physical(self):
        coord = self.coord
        x0, x1 = self.get_indep()
        if coord != 'physical':
            x0 = x0.copy()
            x1 = x1.copy()
            x0, x1 = getattr(self, '_' + coord + '_to_physical')(x0, x1)
        return (x0, x1)

    def get_world(self):
        coord = self.coord
        x0, x1 = self.get_indep()
        if coord != 'world':
            x0 = x0.copy()
            x1 = x1.copy()
            x0, x1 = getattr(self, '_' + coord + '_to_world')(x0, x1)
        return (x0, x1)

    # For compatibility with old Sherpa keywords
    get_image = get_logical
    get_wcs = get_world

    def set_coord(self, coord):
        coord = str(coord).strip().lower()
        # Destroys original data to conserve memory for big imgs
        good = ('logical', 'image', 'physical', 'world', 'wcs')
        if coord not in good:
            raise DataErr('badchoices', 'coordinates', coord, ", ".join(good))

        if coord.startswith('wcs'):
            coord = 'world'
        elif coord.startswith('image'):
            coord = 'logical'

        func = getattr(self, 'get_' + coord)
        self.indep = func()
        self._set_coord(coord)

    def get_filter_expr(self):
        if self._region is not None:
            return str(self._region)
        return ''

    get_filter = get_filter_expr

    def notice2d(self, val=None, ignore=False):
        """Apply a 2D filter.

        Parameters
        ----------
        val : str or None, optional
            The filter to apply. It can be a region string or a
            filename.
        ignore : bool, optional
            If set then the filter should be ignored, not noticed.

        """

        ignore = bool_cast(ignore)

        # This was originally a bit-more complex, but it has been
        # simplified.
        #
        if val is None:
            self.mask = not ignore
            self._region = None
            return

        if not regstatus:
            raise ImportErr('importfailed', 'region', 'notice2d')

        # Crete the new region
        #
        val = str(val).strip()
        isfile = os.path.isfile(val)
        reg = Region(val, isfile)

        # Calculate the mask for this region as an "included"
        # region.
        #
        mask = reg.mask(self.get_x0(), self.get_x1())
        mask = mask.astype(bool)

        # Apply the new mask to the existing mask.
        #
        if not ignore:
            if self.mask is True:
                self.mask = mask
            else:
                self.mask |= mask
        else:
            # Invert the response from region_mask
            mask = ~mask
            if self.mask is False:
                self.mask = mask
            else:
                self.mask &= mask

        # Create the new region shape.
        #
        if self._region is None:
            if ignore:
                reg.invert()

            self._region = reg
        else:
            self._region = self._region.combine(reg, ignore)

    def get_bounding_mask(self):
        mask = self.mask
        shape = None
        if numpy.iterable(self.mask):
            # create bounding box around noticed image regions
            mask = numpy.array(self.mask).reshape(*self.shape)

            # TODO: should replace 'mask == True' with mask but
            # not sure we have a good set of tests
            x0_i, x1_i = numpy.where(mask == True)

            x0_lo = x0_i.min()
            x0_hi = x0_i.max()
            x1_lo = x1_i.min()
            x1_hi = x1_i.max()

            # TODO: subset mask and then ask its shape
            shape = mask[x0_lo:x0_hi + 1, x1_lo:x1_hi + 1].shape
            mask = mask[x0_lo:x0_hi + 1, x1_lo:x1_hi + 1]

            mask = mask.ravel()
        return mask, shape

    def get_img(self, yfunc=None):
        # FIXME add support for coords to image class -> DS9
        self._check_shape()
        y_img = self.filter_region(self.get_dep(False))
        if yfunc is not None:
            m = self.eval_model_to_fit(yfunc)
            if numpy.iterable(self.mask):
                # if filtered, the calculated model must be padded up
                # to the data size to preserve img shape and WCS coord
                m = pad_bounding_box(m, self.mask)
            y_img = (y_img, self.filter_region(m))

        if yfunc is not None:
            y_img = (y_img[0].reshape(*self.shape),
                     y_img[1].reshape(*self.shape))
        else:
            y_img = y_img.reshape(*self.shape)
        return y_img

    def get_axes(self):
        # FIXME: how to filter an axis when self.mask is size of self.y?
        self._check_shape()

        # dummy placeholders needed b/c img shape may not be square!
        axis0 = numpy.arange(self.shape[1], dtype=float) + 1.
        axis1 = numpy.arange(self.shape[0], dtype=float) + 1.
        dummy0 = numpy.ones(axis0.size, dtype=float)
        dummy1 = numpy.ones(axis1.size, dtype=float)

        if self.coord == 'physical':
            axis0, dummy = self._logical_to_physical(axis0, dummy0)
            dummy, axis1 = self._logical_to_physical(dummy1, axis1)

        elif self.coord == 'world':
            axis0, dummy = self._logical_to_world(axis0, dummy0)
            dummy, axis1 = self._logical_to_world(dummy1, axis1)

        return (axis0, axis1)

    def get_x0label(self):
        "Return label for first dimension in 2-D view of independent axis/axes"
        if self.coord in ('logical', 'image'):
            return 'x0'
        elif self.coord in ('physical',):
            return 'x0 (pixels)'
        elif self.coord in ('world', 'wcs'):
            return 'RA (deg)'
        else:
            return 'x0'

    def get_x1label(self):
        """
        Return label for second dimension in 2-D view of independent axis/axes
        """
        if self.coord in ('logical', 'image'):
            return 'x1'
        elif self.coord in ('physical',):
            return 'x1 (pixels)'
        elif self.coord in ('world', 'wcs'):
            return 'DEC (deg)'
        else:
            return 'x1'

    def to_contour(self, yfunc=None):
        y = self.filter_region(self.get_dep(False))
        if yfunc is not None:
            m = self.eval_model_to_fit(yfunc)
            if numpy.iterable(self.mask):
                # if filtered, the calculated model must be padded up
                # to the data size to preserve img shape and WCS coord
                m = self.filter_region(pad_bounding_box(m, self.mask))
            y = (y, m)

        return (self.get_x0(),
                self.get_x1(),
                y,
                self.get_x0label(),
                self.get_x1label())

    def filter_region(self, data):
        if data is not None and numpy.iterable(self.mask):
            filter = numpy.ones(len(self.mask), dtype=SherpaFloat)
            filter[~self.mask] = numpy.nan
            return data * filter
        return data


class DataIMGInt(DataIMG):
    _fields = Data2DInt._fields + ("sky", "eqpos", "coord")

    def __init__(self, name, x0lo, x1lo, x0hi, x1hi, y, shape=None,
                 staterror=None, syserror=None, sky=None, eqpos=None,
                 coord='logical', header=None):
        self._region = None
        self.sky = sky
        self.eqpos = eqpos
        self.coord = coord
        self.header = {} if header is None else header
        self.shape = shape
        Data.__init__(self, name, (x0lo, x1lo, x0hi, x1hi), y, staterror, syserror)

    def _init_data_space(self, filter, *data):
        return IntegratedDataSpace2D(filter, *data)

    def get_logical(self):
        coord = self.coord
        x0lo, x1lo, x0hi, x1hi = self.get_indep()
        if coord != 'logical':
            x0lo = x0lo.copy()
            x1lo = x1lo.copy()
            convert = getattr(self, '_' + coord + '_to_logical')
            x0lo, x1lo = convert(x0lo, x1lo)

            x0hi = x0hi.copy()
            x1hi = x1hi.copy()
            x0hi, x1hi = convert(x0hi, x1hi)

        return (x0lo, x1lo, x0hi, x1hi)

    def get_physical(self):
        coord = self.coord
        x0lo, x1lo, x0hi, x1hi = self.get_indep()
        if coord != 'physical':
            x0lo = x0lo.copy()
            x1lo = x1lo.copy()
            convert = getattr(self, '_' + coord + '_to_physical')
            x0lo, x1lo = convert(x0lo, x1lo)

            x0hi = x0hi.copy()
            x1hi = x1hi.copy()
            x0hi, x1hi = convert(x0hi, x1hi)

        return (x0lo, x1lo, x0hi, x1hi)

    def get_world(self):
        coord = self.coord
        x0lo, x1lo, x0hi, x1hi = self.get_indep()
        if coord != 'world':
            x0lo = x0lo.copy()
            x1lo = x1lo.copy()
            convert = getattr(self, '_' + coord + '_to_world')
            x0lo, x1lo = convert(x0lo, x1lo)

            x0hi = x0hi.copy()
            x1hi = x1hi.copy()
            x0hi, x1hi = convert(x0hi, x1hi)

        return (x0lo, x1lo, x0hi, x1hi)

    def get_axes(self):
        # FIXME: how to filter an axis when self.mask is size of self.y?
        self._check_shape()

        # dummy placeholders needed b/c img shape may not be square!
        axis0lo = numpy.arange(self.shape[1], dtype=float) - 0.5
        axis1lo = numpy.arange(self.shape[0], dtype=float) - 0.5

        axis0hi = numpy.arange(self.shape[1], dtype=float) + 0.5
        axis1hi = numpy.arange(self.shape[0], dtype=float) + 0.5

        dummy0 = numpy.ones(axis0lo.size, dtype=float)
        dummy1 = numpy.ones(axis1lo.size, dtype=float)

        if self.coord == 'physical':
            axis0lo, dummy = self._logical_to_physical(axis0lo, dummy0)
            axis0hi, dummy = self._logical_to_physical(axis0hi, dummy0)

            dummy, axis1lo = self._logical_to_physical(dummy1, axis1lo)
            dummy, axis1hi = self._logical_to_physical(dummy1, axis1hi)

        elif self.coord == 'world':
            axis0lo, dummy = self._logical_to_world(axis0lo, dummy0)
            axis0hi, dummy = self._logical_to_world(axis0hi, dummy0)

            dummy, axis1lo = self._logical_to_world(dummy1, axis1lo)
            dummy, axis1hi = self._logical_to_world(dummy1, axis1hi)

        return (axis0lo, axis1lo, axis0hi, axis1hi)
