#
#  Copyright (C) 2008, 2015 - 2024
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

Design
------

The `DataPHA` class adds support for grouping data - effectively
reducing the number of values in a data set - and adds an extra way
to filter the data with the quality array. The class extends
`~sherpa.data.Data1D`, since the primary data is channels and
counts, but it also has to act like an integrated data set
(`~sherpa.data.Data1DInt`) in some cases. In an extension to
OGIP support, there is limited support for the ``BIN_LO`` and
``BIN_HI`` fields provided with Chandra grating data.

The `DataIMG` class extends 2D support for "gridded" data, with
multiple possible coordinate systems (e.g. ``logical``, ``physical``,
and ``world``).  Along with this, spatial filters can be applied,
using the CIAO region syntax [REGION]_.

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

.. [OGIP_92_007] "The OGIP Spectral File Format", https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html

.. [OGIP_92_007a] "The OGIP Spectral File Format Addendum: Changes log ", https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007a/ogip_92_007a.html

.. [CAL_92_002] "The Calibration Requirements for Spectral Analysis (Definition of RMF and ARF file formats)", https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html

.. [CAL_92_002a] "The Calibration Requirements for Spectral Analysis Addendum: Changes log", https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002a/cal_gen_92_002a.html

.. [PRIVATE_KA] Private communication with Keith Arnaud

.. [REGION] https://cxc.harvard.edu/ciao/ahelp/dmregions.html

Examples
--------

Read in a 2D dataset from the file 'image2.fits' and then filter it to
only use those pixels that lie within 45 units from the physical
coordinate 3150,4515:

>>> from sherpa.astro.io import read_image
>>> img = read_image(data_dir + 'image2.fits')
>>> img.set_coord('physical')
>>> img.notice2d('circle(3150,4515,45)')

Read in a PHA dataset from the file '3c273.pi', subtract the background,
filter to only use the data 0.5 to 7 keV, and re-group the data within
this range to have at least 20 counts per group:

>>> from sherpa.astro.io import read_pha
>>> pha = read_pha(data_3c273 + '3c273.pi')
>>> pha.subtract()
>>> pha.set_analysis('energy')
>>> pha.notice(0.5, 7)
>>> pha.group_counts(20)

"""

from __future__ import annotations

import logging
import os
from typing import Any, Callable, Literal, Mapping, Optional, Sequence, \
    Union, cast, overload
import warnings

import numpy as np

from sherpa.astro import hc
from sherpa.data import Data1DInt, Data2D, Data, Data1D, \
    IntegratedDataSpace2D, _check
from sherpa.models.regrid import EvaluationSpace1D
from sherpa.stats import Chi2XspecVar
from sherpa.utils import pad_bounding_box, interpolate, \
    create_expr, create_expr_integrated, parse_expr, bool_cast, \
    rebin, filter_bins
from sherpa.utils.err import ArgumentTypeErr, DataErr, ImportErr
from sherpa.utils import formatting
from sherpa.utils.numeric_types import SherpaFloat
from sherpa.utils.types import ArrayType, IdType, ModelFunc, StatErrFunc

# There are currently (Sep 2015) no tests that exercise the code that
# uses the compile_energy_grid symbols.
from sherpa.astro.utils import arf_fold, rmf_fold, filter_resp, \
    compile_energy_grid, do_group, expand_grouped_mask

info = logging.getLogger(__name__).info
warning = logging.getLogger(__name__).warning

regstatus = False
try:
    from sherpa.astro.utils._region import Region  # type: ignore
    regstatus = True
except ImportError:
    warning('failed to import sherpa.astro.utils._region; Region routines ' +
            'will not be available')

groupstatus = False
try:
    import group as pygroup  # type: ignore
    groupstatus = True
except ImportError:
    groupstatus = False
    warning('the group module (from the CIAO tools package) is not ' +
            'installed.\nDynamic grouping functions will not be available.')


__all__ = ('DataARF', 'DataRMF', 'DataPHA', 'DataIMG', 'DataIMGInt', 'DataRosatRMF')


AnalysisType = Literal["channel", "energy", "wavelength"]

RateType = Literal["counts", "rate"]

# can arf/rmf be sent ARF1D/RMF1D too?
#
def _notice_resp(chans: np.ndarray,
                 arf: Optional[DataARF],
                 rmf: Optional[DataRMF]
                 ) -> None:
    """Notice the channel range for the associated responses

    Parameters
    ----------
    chans : ndarray
        The noticed channel values.
    arf : sherpa.astro.data.DataARF or None
    rmf : sherpa.astro.data.DataRMF or None

    """

    if rmf is not None and arf is not None:

        bin_mask = rmf.notice(chans)
        if len(rmf.energ_lo) == len(arf.energ_lo):
            arf.notice(bin_mask)

        # If the response is mis-matched, determine which energy bins in the
        # RMF correspond to energy bins in the ARF and which are noticed.
        # Propagate the noticed RMF energy bins to the ARF energy bins.
        elif len(rmf.energ_lo) < len(arf.energ_lo):
            arf_mask = None
            if bin_mask is not None:
                arf_mask = np.zeros(len(arf.energ_lo), dtype=bool)
                for ii, val in enumerate(bin_mask):
                    if val:
                        los = (rmf.energ_lo[ii],)
                        his = (rmf.energ_hi[ii],)
                        grid = (arf.energ_lo, arf.energ_hi)
                        # TODO: should this set integrated=True?
                        #       we only have one test of this code in
                        #       sherpa/astro/tests/test_astro.py:test_missmatch_arf
                        mask = filter_bins(los, his, grid)
                        if mask is None:
                            raise DataErr("notmask")

                        idx = mask.nonzero()[0]
                        arf_mask[idx] = True

            arf.notice(arf_mask)

    elif rmf is not None:
        rmf.notice(chans)

    elif arf is not None:
        # Is this the correct thing to do?
        arf.notice(None)


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


def _extract_fields(obj, summary):
    """Extract the "column" fields.

    Write out the _fields values (that are not None) for the Data
    object. The _extra_fields values are ignored, as they are assumed
    to be included in separate "metadata" section.

    Parameters
    ----------
    obj : Data instance
    summary : str
        The label for the details tab.

    Returns
    -------
    html : str
        The HTML for this section.

    """

    meta = []
    fields = [f for f in obj._fields if f != 'name']
    for f in fields:
        v = getattr(obj, f)
        if v is None:
            continue

        meta.append((f.upper(), v))

    return formatting.html_section(meta, summary=summary,
                                   open_block=True)


def html_pha(pha):
    """HTML representation: PHA"""

    from sherpa.astro.plot import DataPHAPlot
    from sherpa import plot

    ls = []

    plotter = DataPHAPlot()
    plotter.prepare(pha)

    try:
        out = plot.backend.as_html_plot(plotter, 'PHA Plot')
    except AttributeError:
        out = None

    if out is None:
        out = _extract_fields(pha, 'PHA Data')

    ls.append(out)

    # Summary properties
    meta = []
    if pha.name is not None and pha.name != '':
        meta.append(('Identifier', pha.name))

    if pha.exposure is not None:
        meta.append(('Exposure', f'{pha.exposure:g} s'))

    meta.append(('Number of bins', len(pha.channel)))
    meta.append(('Channel range', f'{int(pha.channel[0])} - {int(pha.channel[-1])}'))

    # Although assume the counts are integers, do not force this
    cmin = pha.counts.min()
    cmax = pha.counts.max()
    meta.append(('Count range', f'{cmin} - {cmax}'))

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
            msg = f'Applied ({ngrp} groups)'
        else:
            msg = 'Not applied'

        meta.append(('Grouping', msg))

    # Should this only be displayed if a filter has been applied?
    #
    fexpr = pha.get_filter_expr()
    bintype = 'groups' if pha.grouped else 'channels'
    nbins = pha.get_dep(filter=True).size
    meta.append(('Using', f'{fexpr} with {nbins} {bintype}'))

    ls.append(formatting.html_section(meta, summary='Summary',
                                      open_block=True))

    # TODO:
    #   correction factors

    # Display a subset of header values
    # - maybe don't display the FILTER if NONE
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
                          ('HDUCLAS4', 'PHA format'),
                          ('XFLT0001', 'XSPEC filter 0001')])

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
    erange = f'{emin:g} - {emax:g} keV'

    # Randomly pick 1% as the cut-off for a constant bin width
    #
    de = np.abs(ehi - elo)
    demin = de.min()
    demax = de.max()
    if demin > 0.0:
        dedelta = (demax - demin) / demin
    else:
        dedelta = 1

    if dedelta <= 0.01:
        erange += f', bin size {demax:g} keV'
    else:
        erange += f', bin size {demin:g} - {demax:g} keV'

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
    wrange = f'{wmin:g} - {wmax:g} &#8491;'

    # Randomly pick 1% as the cut-off for a constant bin width
    #
    dw = np.abs(whi - wlo)
    dwmin = dw.min()
    dwmax = dw.max()
    if dwmin > 0.0:
        dwdelta = (dwmax - dwmin) / dwmin
    else:
        dwdelta = 1

    if dwdelta <= 0.01:
        wrange += f', bin size {dwmax:g} &#8491;'
    else:
        wrange += f', bin size {dwmin:g} - {dwmax:g} &#8491;'

    return wrange


def html_arf(arf):
    """HTML representation: ARF"""

    # Unlike the string representation, this provides extra
    # information (e.g. energy range covered). Should it include
    # any filters or masks? How about bin_lo/hi values?
    #
    # It also assumes the units are keV/cm^2 which is not
    # guaranteed.

    from sherpa.astro.plot import ARFPlot
    from sherpa import plot

    ls = []

    plotter = ARFPlot()
    plotter.prepare(arf)

    try:
        out = plot.backend.as_html_plot(plotter, 'ARF Plot')
    except AttributeError:
        out = None

    if out is None:
        out = _extract_fields(arf, 'ARF Data')

    ls.append(out)

    # Summary properties
    meta = []
    if arf.name is not None and arf.name != '':
        meta.append(('Identifier', arf.name))

    if arf.exposure is not None:
        meta.append(('Exposure', f'{arf.exposure:g} s'))

    meta.append(('Number of bins', len(arf.specresp)))

    erange = _calc_erange(arf.energ_lo, arf.energ_hi)
    meta.append(('Energy range', erange))

    # repeat for wavelengths (without the energy threshold)
    #
    if arf.bin_lo is not None and arf.bin_hi is not None:
        wrange = _calc_wrange(arf.bin_lo, arf.bin_hi)
        meta.append(('Wavelength range', wrange))

    a1 = np.min(arf.specresp)
    a2 = np.max(arf.specresp)
    meta.append(('Area range', f'{a1:g} - {a2:g} cm<sup>2</sup>'))

    ls.append(formatting.html_section(meta, summary='Summary',
                                      open_block=True))

    # Display a subset of header values
    # - maybe don't display the FILTER if NONE
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

    from sherpa.astro.plot import RMFPlot
    from sherpa import plot

    ls = []

    plotter = RMFPlot()
    plotter.prepare(rmf)

    try:
        out = plot.backend.as_html_plot(plotter, 'RMF Plot')
    except AttributeError:
        out = None

    if out is None:
        out = _extract_fields(rmf, 'RMF Data')


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
        erange += f' (minimum threshold of {rmf.ethresh} was used)'

    meta.append(('Energy range', erange))

    meta.append(('Channel range', f'{int(rmf.offset)} - {int(rmf.offset + rmf.detchans - 1)}'))

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

    Special-case of the Data2D handling. It would be nice to reuse
    parts of the superclass behavior.
    """
    from sherpa.astro.plot import DataIMGPlot
    from sherpa import plot
    ls = []
    dtype = type(img).__name__

    plotter = DataIMGPlot()
    plotter.prepare(img)

    try:
        out = plot.backend.as_html_image(plotter, f'{dtype} Plot')
    except AttributeError:
        out = None

    if out is None:
        # Only add prefix to summary if there's no plot
        summary = f'{dtype} '

        # Summary properties
        #
        meta = []
        if img.name is not None and img.name != '':
            meta.append(('Identifier', img.name))

        # shape is better defined for DataIMG than Data2D
        meta.append(('Shape',
                     (f'{img.shape[1]} by {img.shape[0]} pixels')))

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
                                          summary=f'Coordinates: {img.sky.name}'))

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
                                          summary=f'Coordinates: {img.eqpos.name}'))

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
            raise ValueError(f"The energy arrays must have the same size, not {elo.size} and {ehi.size}")

        if ethresh is not None and ethresh <= 0.0:
            raise ValueError("ethresh is None or > 0")

        if (elo >= ehi).any():
            # raise DataErr('ogip-error', rtype, label,
            #               'has at least one bin with ENERG_HI < ENERG_LO')
            wmsg = f"The {rtype} '{label}' " + \
                   'has at least one bin with ENERG_HI < ENERG_LO'
            warnings.warn(wmsg)

        # if elo is monotonically increasing, all elements will be True
        #                         decreasing,                      False
        #
        # so the sum will be number of elements or 0
        #
        increasing = np.diff(elo, n=1) > 0.0
        nincreasing = increasing.sum()
        if nincreasing > 0 and nincreasing != len(increasing):
            # raise DataErr('ogip-error', rtype, label,
            #               'has a non-monotonic ENERG_LO array')
            wmsg = f"The {rtype} '{label}' " + \
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
                                  f'value of {ethresh}')

                elo = elo.copy()
                elo[startidx] = ethresh
                wmsg = "The minimum ENERG_LO in the " + \
                       f"{rtype} '{label}' was 0 and has been " + \
                       f"replaced by {ethresh}"
                warnings.warn(wmsg)

            elif e0 < 0.0:
                # raise DataErr('ogip-error', rtype, label,
                #               'has an ENERG_LO value < 0')
                wmsg = f"The {rtype} '{label}' " + \
                       'has an ENERG_LO value < 0'
                warnings.warn(wmsg)

        return elo, ehi

    def _get_data_space(self, filter=False):
        # TODO: the class has no _lo/_hi attributes so what is this
        #       meant to do?
        return EvaluationSpace1D(self._lo, self._hi)


class DataARF(DataOgipResponse):
    """ARF data set.

    The ARF format is described in OGIP documents [CAL_92_002]_ and
    [CAL_92_002a]_.

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

    """
    _ui_name = "ARF"
    _fields = ("name", "energ_lo", "energ_hi", "specresp", "bin_lo", "bin_hi")
    _extra_fields = ("exposure", "ethresh")

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

        super().__init__(name, energ_lo, energ_hi, specresp)

        self._xlabel = 'Energy (keV)'
        self._ylabel = None

    def _repr_html_(self) -> str:
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
        if args:
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

    def get_ylabel(self, yfunc=None) -> str:
        """The label for the dependent axis.

        Parameters
        ----------
        yfunc
           Unused.

        See Also
        --------
        set_ylabel

        """
        if self._ylabel is not None:
            return self._ylabel

        # This will pick up the latest plot backend.
        #
        from sherpa import plot
        return 'cm' + plot.backend.get_latex_for_string('^2')


class DataRMF(DataOgipResponse):
    """RMF data set.

    The RMF format is described in OGIP documents [CAL_92_002]_ and
    [CAL_92_002a]_.

    Parameters
    ----------
    name : str
        The name of the data set; often set to the name of the file
        containing the data.
    detchans : int
    energ_lo, energ_hi : array
        The values of the ENERG_LO and ENERG_HI columns for the RMF
        (for the MATRIX block). The ENERG_HI values must be greater
        than the ENERG_LO values for each bin, and the energy arrays
        must be in increasing or decreasing order.
    n_grp, f_chan, n_chan, matrix : array-like
    offset : int, optional
        This must be 0 or greater. This maps to the TLMIN value of
        the F_CHAN column (when reading from a FITS file).
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

    """
    _ui_name = "RMF"
    _fields = ("name", "energ_lo", "energ_hi", "n_grp", "f_chan", "n_chan", "matrix", "e_min",
               "e_max")
    _extra_fields = ("detchans", "offset", "ethresh")

    def __init__(self, name, detchans, energ_lo, energ_hi, n_grp, f_chan,
                 n_chan, matrix, offset=1, e_min=None, e_max=None,
                 header=None, ethresh=None):
        energ_lo, energ_hi = self._validate(name, energ_lo, energ_hi, ethresh)

        # Check it's an integer.
        if not float(offset).is_integer():
            raise ValueError(f"offset must be an integer, not {offset}")

        if offset < 0:
            raise ValueError(f"offset must be >=0, not {offset}")

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

        # It is assumed, but not yet required, that the RMF components
        # are set with the __init__ call, and not changed after the
        # fact (this is to avoid having to have a complex system of
        # checks as used by the parent Data class for the independent
        # and dependent axes). It would make sense to make these
        # fields either read-only or add validation when changed, but
        # either approach is a large change, so stick with the simple
        # validation here for now (but allow for them to be set
        # later).
        #
        if (self.e_min is not None) ^ (self.e_max is not None):
            raise DataErr("e_min/max must both be set or empty")

        if self.e_min is not None:
            nelo = len(self.e_min)
            nehi = len(self.e_max)
            if nelo != nehi:
                raise DataErr(f"e_min/max mismatch in size: {nelo} vs {nehi}")

            if nelo != self.detchans:
                raise DataErr(f"detchans mismatch with e_min/max: {self.detchans} vs {nelo}")

        super().__init__(name, energ_lo, energ_hi, matrix)

        # special case the xlabel behaviour
        if self.e_min is None or self.e_max is None:
            self._xlabel = 'Channel'
        else:
            self._xlabel = 'Energy (keV)'

        self._ylabel = 'Counts'

    # Although we have a Data1DInt-like dataset, the dependent axis
    # does not match the independent axes. So we have to remove the
    # checks on y. We could enforce the internal constraints but this
    # is hard to do sensibly, so just do not bother.
    #
    @property
    def y(self):
        """The dependent axis."""
        return self._y

    @y.setter
    def y(self, val):
        self._y = val

    def _repr_html_(self) -> str:
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
        if args:
            (rmf, pha) = args
            if pha != () and len(pha[0]) > len(rmf[0]):
                src = rebin(src, pha[0], pha[1], rmf[0], rmf[1])

        if len(src) != len(self._lo):
            raise TypeError("Mismatched filter between ARF and RMF " +
                            "or PHA and RMF")

        return rmf_fold(src, self._grp, self._fch, self._nch, self._rsp,
                        self.detchans, self.offset)

    @overload
    def notice(self, noticed_chans: None) -> None:
        ...

    @overload
    def notice(self, noticed_chans: np.ndarray) -> np.ndarray:
        ...

    # Note that noticed_chans is not a mask, but the channel values.
    def notice(self, noticed_chans=None):
        """Filter the response to match the requested channels."""

        self._fch = self.f_chan
        self._nch = self.n_chan
        self._grp = self.n_grp
        self._rsp = self.matrix
        self._lo = self.energ_lo
        self._hi = self.energ_hi

        # This could also return here if noticed_chans contains all
        # channels, but that is harder to check.
        #
        if noticed_chans is None:
            return None

        (self._grp, self._fch, self._nch, self._rsp,
         bin_mask) = filter_resp(noticed_chans, self.n_grp, self.f_chan,
                                 self.n_chan, self.matrix, self.offset)
        self._lo = self.energ_lo[bin_mask]
        self._hi = self.energ_hi[bin_mask]
        return bin_mask

    def get_indep(self, filter=False):
        return (self._lo, self._hi)

    def get_dep(self, filter=False):
        return self.apply_rmf(np.ones(self.energ_lo.shape, SherpaFloat))


# FIXME There are places in the code that explicitly check if an object is an instance of sherpa.astro.data.DataRMF.
# So it's safer to make DataRosatRMF a subclass of the default class, although in principle they should be siblings
# and subclasses of the same superclass.
class DataRosatRMF(DataRMF):
    ui_name = "ROSAT RMF"

    def _validate(self, name, energy_lo, energy_hi, ethresh):
        return energy_lo, energy_hi


def validate_wavelength_limits(wlo, whi, emax):
    """Check that the wavelength limits are sensible.

    This is used by DataPHA.notice to ensure that the wavelength
    limits are meaningful. It deals with converting to energy
    and handling 0 limits.

    Parameters
    ----------
    wlo, whi : number or None
        The wavelength limit for the low and high edges, or None.
        We have wlo >= 0 and wlo <= whi.
    emax : number
        The maximum energy of the response (the upper edge of the
        last bin).

    Returns
    -------
    lims : (lo, hi) or None
        The energy filter or None if it doesn't overlap the response.
        The lo and hi values can be None.

    Notes
    -----
    This routine can return a range that doesn't overlap the
    response, but that is handled downstream (i.e. it only returns
    None in certain circumstances, not all cases).

    """

    # As we allow wlo and whi to be 0 we need to handle this here,
    # otherwise we'd have try hc / 0. We can either replace 0 by a
    # value such as numpy.finfo(numpy.float32).tiny which will result
    # in a very-large energy, or we can use domain knowledge - the
    # maximum energy value in the grid and the requested limits. We
    # try the latter, but it's not as simple
    #
    if whi is None:
        lo = None
    elif whi == 0.0:
        lo = -1
    else:
        lo = hc / whi

    if wlo is None:
        hi = None
    elif wlo == 0.0:
        hi = -1
    else:
        hi = hc / wlo

    # If either of the arguments were 0 then we need to ensure they
    # are sensible when combined together. Note that we can't have
    # whi = 0 and wlo > 0 (ie lo < -1 and hi > 0).
    #
    if lo is not None and hi is not None:
        if lo < 0 and hi < 0:
            # Both limits were 0 so we can do nothing
            return None

        if hi < 0:
            # The original query was 0 to x which maps to hc/x to None
            # but we need to know if hc/x is > emax or not
            if lo < emax:
                hi = None
            else:
                # there is no valid filter here
                return None

    elif lo is not None and lo < 0:
        lo = None

    elif hi is not None and hi < 0:
        hi = None

    return lo, hi


def replace_xspecvar_values(src_counts, bkg_counts,
                            staterr, bkg_variances,
                            src_scale, bkg_scales):
    """Handle error/variances when 0's are present using XSPEC var.

    Correct the error for background-subtracted data when using the
    XSPEC variance when at least one of the source or background terms
    has 0 counts in it. This is to address
    https://github.com/sherpa/sherpa/issues/356

    Parameters
    ----------
    src_counts : ndarray
        The number of source counts in each group.
    bkg_counts : list of ndarray
        The number of background counts in each group, for each
        background element.
    staterr : ndarray
        The error terms for the source counts (expected to be the XSPEC
        error function ran on src_counts). This may be changed.
    bkg_variances : list of ndarray
        The variances for the background counts (expected to be the
        square of the XSPEC error function ran on the bkg_counts
        elements). This may be changed. The length of bkg_variances must
        match bkg_counts.
    src_scale : number or ndarray
        The combined exposure, backscal, and areascal for the source.
    bkg_scales : list of number or ndarray
        The combined exposure, backscal, and areascal for the background
        components.

    Notes
    -----
    The chi2xspecvar error returns 1 when given 0 counts, which is
    okay when there's no background, but fails

    - when there's background subtraction, as it doesn't match XSPEC,
      and is addressed here

    - is not at all obvious what to do when using chi2xspecvar but are
      trying to fit the background rather than subtract is,
      and this is *not* addressed here.

    When looking at the simple case of a source bin with a single
    background component then we have to be concerned with

    a) either the source or the background component has 0 counts but
       the other one has positive (non-zero) counts

    b) both source and background have 0 counts.

    For case a, we can just set the error (or variance) value to 0 for
    the element with 0 counts. This essentially drops that element
    from the calculation.

    For case b then, in order to match XSPEC, we have to clear out the
    background component and then replace the source error term with a
    value depending on the value of the relative scaling terms of the
    source and background components.

    """

    # We care about those elements where either
    # - all components (source, background_1 .. n) are 0
    # - some are 0
    #
    counts = np.asarray([src_counts] + bkg_counts)
    numzero = np.sum(counts == 0, axis=0)

    nbkg = len(bkg_variances)
    ncpts = 1 + nbkg
    all_zero = numzero == ncpts
    some_zero = (numzero > 0) & (numzero < ncpts)

    # For the some_zero cases we just set them to 0. Unfortunately
    # identifying them requires another loop as some_zero has
    # collapsed the data.
    #
    idx, = np.asarray(src_counts == 0 & some_zero).nonzero()
    staterr[idx] = 0

    for bcnts, bvar in zip(bkg_counts, bkg_variances):
        idx, = np.asarray(bcnts == 0 & some_zero).nonzero()
        bvar[idx] = 0

    # Do we have any bins where all source and background components
    # are zero? If not, we can return immediately.
    #
    nzero = all_zero.sum()
    if nzero == 0:
        return

    # When both source and background are 0, the replacement
    # depends on how the source and background scaling works.
    # This is based on a discussion with the XSPEC developers
    # and the code in Xspec/sc/XSStat/ChiSquare.cxx, in
    # particular ChiSquare::applyMinVariance. The XSPEC code
    # only deals with a single background component, so here
    # we average the scaling factors from multiple components
    # for the check.
    #
    # The XSPEC code uses a value (for the variance) of the
    # minimum of
    #
    #    1 / src_texp^2
    #    scale^2 / bkg_texp^2
    #
    # where scale = (AREASCAL * BACKSCAL)_src / (AREASCAL * BACKSCAL)_bkg
    #
    # However, XSPEC is calculating a rate and we want counts,
    # so we care about the minimum of (for the variance)
    #
    #   1
    #   scale^2 * src_texp^2 / bkg_texp^2
    #
    # although we actually use the square root of this, as we
    # change the staterr array (to ensure it doesn't get
    # re-scaled when calculating statvar below).
    #

    # The conversion is per-bin, as the scaling factors are per-bin,
    # but they may be specified as a scalar.  In that case, convert
    # into an array to simplify the following. The requirement here is
    # that we have the same channel range in each background component
    # as the source.
    #
    nelem = len(staterr)
    scales = np.zeros((nbkg, nelem))
    for idx, bscale in enumerate(bkg_scales):
        if np.isscalar(bscale):
            scales[idx] = np.ones(nelem) * bscale
        else:
            scales[idx] = bscale

    # Use the average background-scaling factor for the background, as
    # it's not clear what is done in XSPEC for multiple-background
    # cases.
    #
    bkg_scale = scales.sum(axis=0) / nbkg

    s = np.ones(nzero)
    b = (src_scale / bkg_scale)[all_zero]
    combined = np.asarray([s, b])
    minval = np.min(combined, axis=0)

    staterr[all_zero] = minval
    for bvar in bkg_variances:
        bvar[all_zero] = 0


class DataPHA(Data1D):
    """PHA data set, including any associated instrument and background data.

    The PHA format is described in an OGIP document [OGIP_92_007]_ and
    [OGIP_92_007a]_.

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
        The wavelength ranges for the channels. This is intended to support
        Chandra grating spectra.
    grouping : array of int or None, optional
    quality : array of int or None, optional
    exposure : number or None, optional
        The exposure time for the PHA data set, in seconds.
    backscal : scalar or array or None, optional
    areascal : scalar or array or None, optional
    header : dict or None, optional
        If ``None`` the header will be pre-populated with a minimal set of
        keywords that would be found in an OGIP compliant PHA I file.

    Attributes
    ----------
    name : str
        Used to store the file name, for data read from a file.
    exposure

    Notes
    -----
    The original data is stored in the attributes - e.g. `counts` -
    and the data-access methods, such as `get_dep` and
    `get_staterror`, provide any necessary data manipulation to handle
    cases such as: background subtraction, filtering, and grouping.

    There is additional complexity compared to the Data1D case when
    filtering data because:

    * although the data uses channel numbers, users will often want to
      filter the data using derived values (in energy or wavelength
      units, such as 0.5 to 7.0 keV or 16 to 18 Angstroms);

    * although derived from the Data1D case, PHA data is more-properly
      thought about as being an integrated data set, so each channel
      maps to a range of energy or wavelength values;

    * the data is often grouped to improve the signal-to-noise, and so
      requests for values need to determine whether to filter the data
      or not, whether to group the data or not, and how to combine the
      data within each group;

    * and there is also the quality array, which indicates whether or
      not a channel is trust-worthy or not (and so acts as an
      additional filtering term).

    The handling of the AREASCAl value - whether it is a scalar or
    array - is currently in flux. It is a value that is stored with
    the PHA file, and the OGIP PHA standard ([OGIP_92_007]_,
    [OGIP_92_007a]_) describes the observed counts being divided by
    the area scaling before comparison to the model. However, this is
    not valid for Poisson-based statistics, and is also not how XSPEC
    handles AREASCAL ([PRIVATE_KA]_); the AREASCAL values are used to
    scale the exposure times instead. The aim is to add this logic to
    the instrument models in `sherpa.astro.instrument`, such as
    `sherpa.astro.instrument.RMFModelPHA`. The area scaling still has
    to be applied when calculating the background contribution to a
    spectrum, as well as when calculating the data and model values
    used for plots (following XSPEC so as to avoid sharp
    discontinuities where the area-scaling factor changes strongly).

    """
    _fields = ('name', 'channel', 'counts', 'staterror', 'syserror', 'bin_lo', 'bin_hi', 'grouping', 'quality')
    _extra_fields = ('exposure', 'backscal', 'areascal', 'grouped', 'subtracted', 'units', 'rate',
                     'plot_fac', 'response_ids', 'background_ids')

    _related_fields = Data1D._related_fields + ("bin_lo", "bin_hi", "counts", "grouping", "quality",
                                                "backscal", "areascal")

    def __init__(self,
                 name: str,
                 channel: Optional[ArrayType],
                 counts: Optional[ArrayType],
                 staterror: Optional[ArrayType] = None,
                 syserror: Optional[ArrayType] = None,
                 bin_lo: Optional[ArrayType] = None,
                 bin_hi: Optional[ArrayType] = None,
                 grouping: Optional[ArrayType] = None,
                 quality: Optional[ArrayType] = None,
                 exposure=None,
                 backscal=None,
                 areascal=None,
                 header: Optional[Mapping[str, Any]] = None
                 ) -> None:

        # Set the size of the object as soon as we know (it makes it
        # easier to get usable error messages when checking the
        # related fields). This is only done for the channel case to
        # allow the counts field to be set but the channel field
        # unset.
        #
        channel = _check(channel)
        if channel is not None:
            self._size = len(channel)

        counts = _check(counts)

        # Assert types: is there a better way to do this?
        #
        self._bin_lo: Optional[np.ndarray]
        self._bin_hi: Optional[np.ndarray]
        self._grouping: Optional[np.ndarray]
        self._quality: Optional[np.ndarray]
        self._quality_filter: Optional[np.ndarray]

        self.bin_lo = _check(bin_lo)
        self.bin_hi = _check(bin_hi)
        self.quality = _check(quality)
        self.grouping = _check(grouping)
        self.exposure = exposure
        self.backscal = backscal
        self.areascal = areascal
        if header is None:
            header = {"HDUCLASS": "OGIP", "HDUCLAS1": "SPECTRUM",
                      "HDUCLAS2": "TOTAL", "HDUCLAS3": "TYPE:I",
                      "HDUCLAS4": "COUNT", "HDUVERS": "1.2.1",
                      "TELESCOP": "none", "INSTRUME": "none",
                      "FILTER": "none", "POISSERR": True}

        self.header = header
        self._grouped = grouping is not None

        # _original_groups is set False if the grouping is changed via
        # the _dynamic_groups method. This is currently only used by the
        # serialization code (sherpa.astro.ui.serialize) to determine
        # whether to write out the grouping data.
        #
        self._original_groups = True
        self._subtracted = False
        self._response_ids: list[IdType] = []
        self._background_ids: list[IdType] = []

        # Can these responses be ARF1D / RMF1D too? For now use just
        # DataARF/DataRMF. Since the response versions actually just
        # send any call to the underlying Data object it technically
        # is the correct typing, and avoids playing around with
        # self-referential type imports for now.
        #
        self._responses: dict[IdType, tuple[Optional[DataARF],
                                            Optional[DataRMF]]] = {}

        self._backgrounds: dict[IdType, DataPHA] = {}
        self._rate = True
        self._plot_fac = 0
        self.units = "channel"
        self.quality_filter = None
        super().__init__(name, channel, counts, staterror, syserror)

        # special-case the labels
        self._xlabel = None
        self._ylabel = None

    def _get_grouped(self) -> bool:
        return self._grouped

    def _set_grouped(self, val: bool) -> None:
        val = bool(val)

        if val and self.grouping is None:
            raise DataErr('nogrouping', self.name)

        # Short cut if the grouping isn't changing.
        #
        # This could be dangerous, as there may be times we
        # would want to generate things, but it is believed
        # that we have identified all these changes, such
        # as changing the grouping field.
        #
        if val == self._grouped:
            return

        if not np.iterable(self.mask):
            self._grouped = val
            return

        # As the grouping has changed AND there's a partial filter
        # then it is important to re-apply the filter so that the mask
        # matches the new grouping scheme.
        #
        old_filter = self.get_filter()
        self._grouped = val
        self.ignore()
        for vals in parse_expr(old_filter):
            self.notice(*vals)

    grouped = property(_get_grouped, _set_grouped,
                       doc='Are the data grouped?')

    def _get_subtracted(self) -> bool:
        return self._subtracted

    def _set_subtracted(self, val: bool) -> None:
        val = bool(val)
        if len(self._backgrounds) == 0:
            raise DataErr('nobkg', self.name)
        self._subtracted = val

    subtracted = property(_get_subtracted, _set_subtracted,
                          doc='Are the background data subtracted?')

    def _set_units(self, val: str) -> None:
        units = str(val).strip().lower()

        if units == 'bin':
            units = 'channel'

        if units.startswith('chan'):
            units = 'channel'

        elif units.startswith('ener'):
            units = 'energy'

        elif units.startswith('wave'):
            units = 'wavelength'
        else:
            raise DataErr('bad', 'quantity', val)

        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            if bkg.get_response() != (None, None) or \
               (bkg.bin_lo is not None and bkg.bin_hi is not None):
                bkg.units = units

        self._units = units

    def _get_units(self) -> AnalysisType:
        return cast(AnalysisType, self._units)

    units = property(_get_units, _set_units,
                     doc="Units of the independent axis: one of 'channel', 'energy', 'wavelength'.")

    def _get_rate(self) -> bool:
        return self._rate

    def _set_rate(self, val: bool) -> None:
        self._rate = bool_cast(val)
        for bkg_id in self.background_ids:
            self.get_background(bkg_id).rate = val

    rate = property(_get_rate, _set_rate,
                    doc="""Is the Y axis displayed as a rate when plotting data?

When True the y axis is normalised by the exposure time to display
a rate.""")

    def _get_plot_fac(self) -> int:
        return self._plot_fac

    def _set_plot_fac(self, val: int) -> None:
        # I'd prefer to check whether val is an integer, but there may
        # be users who have set the value to 2.0 and it doesn't seem
        # worth breaking that code. We do however want to error out if
        # val=0.5 (say), hence this slightly-contrived check.
        #
        try:
            ival = int(val)
            okay = val == ival
        except (ValueError, TypeError):
            # For when int() can't convert val, which can raise
            # different errors.
            okay = False

        if not okay:
            raise DataErr("bad", "plot_fac setting", val)

        self._plot_fac = ival
        for bkg_id in self.background_ids:
            self.get_background(bkg_id).plot_fac = ival

    plot_fac = property(_get_plot_fac, _set_plot_fac,
                        doc="""How the X axis is used to create the Y axis when plotting data.

The Y axis values are multiplied by X^plot_fac. The default
value of 0 means the X axis is not used in plots. The value
must be an integer.""")

    def _get_response_ids(self) -> list[IdType]:
        return self._response_ids

    def _set_response_ids(self, ids: Sequence[IdType]) -> None:
        if not np.iterable(ids):
            raise DataErr('idsnotarray', 'response', str(ids))

        keys = list(self._responses.keys())
        for resp_id in ids:
            if resp_id not in keys:
                raise DataErr('badids', str(resp_id), 'response', str(keys))

        self._response_ids = list(ids)

    response_ids = property(_get_response_ids, _set_response_ids,
                            doc="""IDs of defined instrument responses (ARF/RMF pairs).

If set, the identifiers must already exist, and any other responses
will be removed. The identifiers can be integers or strings.
""")

    def _get_background_ids(self) -> list[IdType]:
        return self._background_ids

    def _set_background_ids(self, ids: Sequence[IdType]) -> None:
        if not np.iterable(ids):
            raise DataErr('idsnotarray', 'background', str(ids))

        keys = list(self._backgrounds.keys())
        for bkg_id in ids:
            if bkg_id not in keys:
                raise DataErr('badids', str(bkg_id), 'background', str(keys))

        self._background_ids = list(ids)

    background_ids = property(_get_background_ids, _set_background_ids,
                              doc="""IDs of defined background data sets.

If set, the identifiers must already exist, and any other backgrounds
will be removed. The identifiers can be integers or strings.
""")

    @overload
    def _set_related(self,
                     attr: str,
                     val: Optional[ArrayType],
                     check_mask: bool = True,
                     allow_scalar: object = Literal[False],
                     **kwargs
                     ) -> None:
        ...

    @overload
    def _set_related(self,
                     attr: str,
                     val: float,
                     check_mask: bool = True,
                     allow_scalar: object = Literal[True],
                     **kwargs
                     ) -> None:
        ...

    def _set_related(self, attr, val, check_mask=True,
                     allow_scalar=False, **kwargs):
        """Set a field that must match the independent axes size.

        The value can be None, a scalar (if allow_scalar is set), or
        something with the same length as the independent axis. This
        is intended to be used from the property setter.

        """
        if val is None:
            setattr(self, f"_{attr}", None)
            return

        if not np.iterable(val):
            if not allow_scalar:
                raise DataErr("notanarray")

            setattr(self, f"_{attr}", val)
            return

        # At this point val is an iterable, but it can be hard for the
        # type checker to know this, so convert to an ndarray to make
        # sure.
        #
        super()._set_related(attr, np.asarray(val),
                             check_mask=check_mask)

    # Set up the aliases for channel and counts
    #
    @property
    def channel(self) -> Optional[np.ndarray]:
        """The channel array.

        This is the first, and only, element of the indep attribute.
        """
        return self.indep[0]

    @channel.setter
    def channel(self, val: Optional[ArrayType]) -> None:
        # This case is handled by _check but type checkers may not
        # recognize it, so handle separately.
        #
        if val is None:
            self.indep = (None, )
            return

        self.indep = (_check(val), )

    @property
    def counts(self) -> Optional[np.ndarray]:
        """The counts array.

        This is an alias for the y attribute.
        """
        return self.y

    @counts.setter
    def counts(self, val: Optional[np.ndarray]):
        self.y = val

    # Override the mask handling because the mask matches the grouped
    # data length, not the independent axis.
    #
    @Data1D.mask.setter
    def mask(self, val: Union[ArrayType, bool]) -> None:

        # We only need to over-ride the behavior if the data is
        # grouped and val is a sequence (so we test with isscalar
        # rather than iterable, to avoid selecting strings).
        #
        if self.grouped and val is not None and not np.isscalar(val):
            # The assumption is that if the data is grouped then it contains data.
            nexp = len(self.get_y(filter=False))
            if len(val) != nexp:
                raise DataErr("mismatchn", "grouped data", "mask", nexp, len(val))

            self._data_space.filter.mask = val
            return

        # This is a bit messy just to call the original code
        Data1D.mask.fset(self, val)

    # Set up the properties for the related fields
    #
    @property
    def bin_lo(self) -> Optional[np.ndarray]:
        """The lower edge of each channel, in Angstroms, or None.

        The values are expected to be in descending order. This is
        only expected to be set for Chandra grating data.
        """
        return self._bin_lo

    @bin_lo.setter
    def bin_lo(self, val: Optional[ArrayType]) -> None:
        self._set_related("bin_lo", val)

    @property
    def bin_hi(self) -> Optional[np.ndarray]:
        """The upper edge of each channel, in Angstroms, or None.

        The values are expected to be in descending order, with the
        bin_hi value larger than the corresponding bin_lo element.
        This is only expected to be set for Chandra grating data.
        """
        return self._bin_hi

    @bin_hi.setter
    def bin_hi(self, val: Optional[ArrayType]) -> None:
        self._set_related("bin_hi", val)

    @property
    def grouping(self) -> Optional[np.ndarray]:
        """The grouping data.

        A group is indicated by a sequence of flag values starting
        with ``1`` and then ``-1`` for all the channels in the group,
        following [OGIP_92_007]_.  The grouping array must match the
        number of channels and it will be converted to an integer type
        if necessary.

        .. versionchanged:: 4.15.1
           The filter is now re-calculated when the grouping is
           changed. It is suggested that the filter be checked with
           `get_filter` to check it is still sensible. If set to None
           then the group flag is cleared.

        Returns
        -------
        grouping : numpy.ndarray or None

        See Also
        --------
        group, grouped, quality

        """
        return self._grouping

    @grouping.setter
    def grouping(self, val: Optional[ArrayType]) -> None:
        # _set_related checks if it's a scalar value, so we just need
        # to check it's convertible to ndarray.
        #
        if val is not None:
            try:
                val = np.asarray(val, dtype=np.int16)
            except TypeError:
                raise DataErr("notanintarray") from None

        # If this is called within __init__ then we do not want to
        # recreate the filter expression. If the mask is either True
        # or False then we do not want to recreate filter (as a mask
        # of True will get converted to [True] * nbins which then
        # breaks downstream processing, as it looks like notice/ignore
        # has already been called).
        #
        ofilter = None
        if self._NoNewAttributesAfterInit__initialized and np.iterable(self.mask):
            ofilter = self.get_filter()

        self._set_related("grouping", val)

        # If the array has been removed then we need to reset the
        # group flag.
        #
        if val is None and self._NoNewAttributesAfterInit__initialized and self.grouped:
            self.grouped = False

        if ofilter is not None:
            # If the data has been filtered then we want to re-create
            # the filter (which is needed because the mask size may
            # have changed).  This is not ideal, since the actual
            # filter is not necessarily the requested filter (thanks
            # to the width of a group), so this could end up selecting
            # a surprising range.
            #
            self.ignore()
            for vals in parse_expr(ofilter):
                self.notice(*vals)

    @property
    def quality(self) -> Optional[np.ndarray]:
        """The quality data.

        A quality value of 0 indicates a good channel, otherwise
        (values >=1) the channel is considered bad and can be excluded
        using the `ignore_bad` method, as discussed in [OGIP_92_007]_. The
        quality array must match the number of channels and it will be
        converted to an integer type if necessary.

        Returns
        -------
        quality : numpy.ndarray or None

        See Also
        --------
        group, grouping

        """
        return self._quality

    @quality.setter
    def quality(self, val: Optional[ArrayType]) -> None:
        # _set_related checks if it's a scalar value, so we just need
        # to check it's convertible to ndarray.
        #
        if val is not None:
            try:
                val = np.asarray(val, dtype=np.int16)
            except TypeError:
                raise DataErr("notanintarray") from None

        # If the quality changes, should we re-create the filter as we
        # do with grouping? No, because the quality array is currently
        # not automatically applied. Although perhaps we should reset
        # quality_filter?
        #
        self._set_related("quality", val)

    # It is unclear exactly what the quality_filter is meant to
    # represent, so as part of improving support for ignore_bad, add
    # some semantics.
    #
    @property
    def quality_filter(self) -> Optional[np.ndarray]:
        """The applied quality filter.

        If set, this indicates which channels are currently selected
        (a True value) and those that are ignored (False) by a
        `ignore_bad` call.

        Returns
        -------
        quality_filter : numpy.ndarray or None

        See Also
        --------
        ignore_bad, quality

        """
        return self._quality_filter

    # TODO: should this be allowed? Maybe this should only be set
    # internally with ignore_bad?
    #
    @quality_filter.setter
    def quality_filter(self, val: Optional[ArrayType]) -> None:
        if val is None:
            qfilt = None
        else:
            # ensure we store boolean values
            qfilt = np.asarray(val, dtype=bool)
        self._set_related("quality_filter", qfilt)

    @property
    def areascal(self):
        """The area scaling value (can be a scalar or array).

        If this is an array then it must match the length of channel.
        """
        return self._areascal

    @areascal.setter
    def areascal(self, val):
        self._set_related("areascal", val, allow_scalar=True)

    @property
    def backscal(self):
        """The background scaling value (can be a scalar or array).

        If this is an array then it must match the length of channel.
        """
        return self._backscal

    @backscal.setter
    def backscal(self, val):
        self._set_related("backscal", val, allow_scalar=True)

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the PHA
        """
        return html_pha(self)

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self._background_ids = state['_background_ids']
        self._backgrounds = state['_backgrounds']
        self._set_units(state['_units'])

        if 'header' not in state:
            self.header = {}
        self.__dict__.update(state)

    primary_response_id: IdType = 1
    """The identifier for the response component when not set."""

    def set_analysis(self,
                     quantity: str,
                     type: RateType = 'rate',
                     factor: int = 0) -> None:
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

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> pha.set_analysis('energy')

        >>> pha.set_analysis('wave', type='counts', factor=1)
        >>> pha.units
        'wavelength'

        """
        if type not in ["counts", "rate"]:
            raise DataErr("plottype", type, "'rate' or 'counts'")

        try:
            self.plot_fac = factor
        except DataErr:
            # Create a slightly-different error message
            raise DataErr("bad", "factor setting", factor) from None

        self.rate = type == "rate"

        arf, rmf = self.get_response()
        if rmf is not None and rmf.detchans != len(self.channel):
            raise DataErr("incompatibleresp", rmf.name, self.name)

        if (rmf is None and arf is None) and \
           (self.bin_lo is None and self.bin_hi is None) and \
           quantity != "channel":
            raise DataErr("norsp", self.name)

        if rmf is None and arf is not None and quantity != "channel" and \
           len(arf.energ_lo) != len(self.channel):
            raise DataErr("incompleteresp", self.name)

        self.units = quantity

    def get_analysis(self) -> AnalysisType:
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

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> pha.set_analysis("wave")
        >>> is_wave = pha.get_analysis() == 'wavelength'

        """
        return self.units

    def _fix_response_id(self, id: Optional[IdType]) -> IdType:
        if id is not None:
            return id

        return self.primary_response_id

    def get_response(self, id: Optional[IdType] = None
                     ) -> tuple[Optional[DataARF], Optional[DataRMF]]:
        """Return the response component.

        Parameters
        ----------
        id : int, str, or None, optional
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

    def set_response(self,
                     arf: Optional[DataARF] = None,
                     rmf: Optional[DataRMF] = None,
                     id: Optional[IdType] = None
                     ) -> None:
        """Add or replace a response component.

        To remove a response use delete_response(), as setting arf and
        rmf to None here does nothing.

        Parameters
        ----------
        arf : sherpa.astro.data.DataARF instance or None, optional
           The ARF to add if any.
        rmf : sherpa.astro.data.DataRMF instance or None, optional
           The RMF to add, if any.
        id : int, str, or None, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        See Also
        --------
        delete_response, get_response, set_arf, set_rmf

        Notes
        -----
        If the PHA header does not have the TELESCOP, INSTRUME, or
        FILTER header keywords set (or they are set to "none"), then
        they are taken from the ARF or RMF, if they are not set to
        "none". This is to allow simulated data sets to be used with
        external programs, such as XSPEC.

        """
        if (arf is None) and (rmf is None):
            return

        resp_id = self._fix_response_id(id)
        self._responses[resp_id] = (arf, rmf)
        if resp_id not in self.response_ids:
            self.response_ids.append(resp_id)

        # To support simulated data (e.g. issue #1209) we over-write
        # the header TELESCOP/INSTRUME/FILTER settings to match the
        # response if they are at their "default" value (e.g. "none").
        #
        # TODO: should we warn if these values don't match (when not
        # "none")?
        #
        def set_key(key):
            val = self.header.get(key, "none")
            if val != "none":
                return

            if arf is not None:
                val = arf.header.get(key, "none")
                if val != "none":
                    self.header[key] = val
                    return

            if rmf is None:
                return

            val = rmf.header.get(key, "none")
            if val == "none":
                return

            self.header[key] = val

        set_key("TELESCOP")
        set_key("INSTRUME")
        set_key("FILTER")

    def delete_response(self, id: Optional[IdType] = None) -> None:
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
        resp_id = self._fix_response_id(id)
        self._responses.pop(resp_id, None)
        self.response_ids.remove(resp_id)

    def get_arf(self, id: Optional[IdType] = None) -> Optional[DataARF]:
        """Return the ARF from the response.

        Parameters
        ----------
        id : int, str, or None, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        Returns
        -------
        arf: sherpa.astro.data.DataARF instance or None
           The ARF, if set.

        See Also
        --------
        get_response, get_rmf, get_full_responses

        """
        return self.get_response(id)[0]

    def get_rmf(self, id: Optional[IdType] = None) -> Optional[DataRMF]:
        """Return the RMF from the response.

        Parameters
        ----------
        id : int, str, or None, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        Returns
        -------
        rmf: sherpa.astro.data.DataRMF instance or None
           The RMF, if set.

        See Also
        --------
        get_arf, get_response, get_full_responses

        """
        return self.get_response(id)[1]

    def set_arf(self,
                arf: DataARF,
                id: Optional[IdType] = None) -> None:
        """Add or replace the ARF in a response component.

        This replaces the existing ARF of the response, keeping the
        previous RMF (if set). Use the delete_response method to
        remove the response, rather than setting arf to None.

        Parameters
        ----------
        arf : sherpa.astro.data.DataARF instance
           The ARF to add.
        id : int, str, or None, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        See Also
        --------
        delete_response, set_response, set_rmf

        """
        self.set_response(arf, self.get_rmf(id), id)

    def set_rmf(self,
                rmf: DataRMF,
                id: Optional[IdType] = None) -> None:
        """Add or replace the RMF in a response component.

        This replaces the existing RMF of the response, keeping the
        previous ARF (if set). Use the delete_response method to
        remove the response, rather than setting rmf to None.

        Parameters
        ----------
        rmf : sherpa.astro.data.DataRMF instance
           The RMF to add.
        id : int, str, or None, optional
           The identifier of the response component. If it is None
           then the default response identifier is used.

        See Also
        --------
        delete_response, set_response, set_arf

        """
        self.set_response(self.get_arf(id), rmf, id)

    def get_specresp(self,
                     filter: bool = False
                     ) -> Optional[np.ndarray]:
        """Return the effective area values for the data set.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the ARF or not. The default is `False`.

        Returns
        -------
        arf : array or None
           The effective area values for the data set (or background
           component) if set.

        Notes
        -----
        This will return `None` when a RSP file (a combined ARF and
        RMF) is used, rather than separate responses. The relationship
        between RSP, ARF, and RMF is described in
        `OGIP Calibration Memo CAL/GEN/92-002
        <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html>`_
        and
        `OGIP Calibration Memo CAL/GEN/92-002a
        <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002a/cal_gen_92_002a.html>`_.

        """
        self.notice_response(False)
        arf, rmf = self.get_response()

        # It's not clear why we do interpolation below, why we replace
        # with 1 rather than 0, or how it is even meant to work, since
        # the current code returns different values depending on the
        # units setting - see issue #1582
        #
        if arf is None or rmf is None:
            return None

        specresp = arf.get_dep()
        elo, ehi = arf.get_indep()
        lo, hi = self._get_ebins(group=False)

        newarf = interpolate(lo, elo, specresp)
        newarf[newarf <= 0] = 1.

        if bool_cast(filter):
            newarf = self.apply_filter(newarf, self._middle)

        return newarf

    def get_full_response(self, pileup_model=None):
        """Calculate the response for the dataset.

        Unlike `get_response`, which returns a single response, this function
        returns all responses for datasets that have multiple responses set
        and it offers the possibility to include a pile-up model.

        Parameters
        ----------
        pileup_model : None or a `sherpa.astro.models.JDPileup` instance
            If a pileup model shall be included in the return, then it needs
            to be passed in.

        Returns
        -------
        response
           The return value depends on whether an ARF, RMF, or pile up
           model has been associated with the data set.

        See Also
        --------
        get_response, get_arf, get_rmf
        """
        # import is here because sherpa.astro.instrument depends on
        # sherpa.astro.data. Importing here instead of on the top
        # avoids a circular import.
        from sherpa.astro import instrument

        if pileup_model is not None:
            return instrument.PileupResponse1D(self, pileup_model)

        if len(self._responses) > 1:
            return instrument.MultipleResponse1D(self)

        return instrument.Response1D(self)

    def _get_ebins(self,
                   response_id: Optional[IdType] = None,
                   group: bool = True
                   ) -> tuple[np.ndarray, np.ndarray]:
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
        response_id : int, str, or None, optional
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

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
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
        True
        """

        if self.size is None:
            raise DataErr("sizenotset", self.name)

        # This is primarily to make the type checker happy, since it
        # is expected to be called with channels set.
        #
        if self.channel is None:
            raise DataErr(f"data set '{self.name}' has no channel information")

        if self.units == 'channel':
            elo = self.channel
            ehi = self.channel + 1
        elif (self.bin_lo is not None) and (self.bin_hi is not None):
            # TODO: review whether bin_lo/hi should over-ride the response
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

        if self.grouped and bool_cast(group):
            elo = self.apply_grouping(elo, self._min)
            ehi = self.apply_grouping(ehi, self._max)

            if len(elo) == 0:
                raise DataErr('notmask')

        # apply_grouping applies a quality filter to the output
        # but if we get here then there is no equivalent. This
        # is likely confusing, at best, but we don't have good
        # tests to check what we should be doing.
        #
        return (elo, ehi)

    def get_indep(self,
                  filter: bool = True
                  ) -> Union[tuple[np.ndarray, ...],
                             tuple[None, ...]]:

        # short-cut if no data
        if self.size is None:
            return (None,)

        if filter:
            return (self.get_noticed_channels(),)

        return (self.channel,)

    def _get_indep(self,
                   filter: bool = False
                   ) -> tuple[np.ndarray, np.ndarray]:
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

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> pha.units = 'energy'
        >>> elo, ehi = pha._get_indep()
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
        >>> wlo, whi = pha._get_indep()
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
            for resp_id in self.response_ids:
                arf, rmf = self.get_response(resp_id)
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
                      np.equal(energylist[0], None).any())):
                raise DataErr('noenergybins', 'Response')
            else:
                elo, ehi = energylist[0]

        lo, hi = elo, ehi
        if self.units == 'wavelength':
            lo = hc / ehi
            hi = hc / elo

        return (lo, hi)

    default_background_id: IdType = 1
    """The identifier for the background component when not set.

It is an integer or string.
"""

    def _fix_background_id(self, id: Optional[IdType]) -> IdType:
        """Identify the background identifier.

        Parameters
        ----------
        id : int, str, or None
            The background identifier, or None. If None then the
            default_background_id value will be used.

        Returns
        -------
        bkg_id : int or str
            The background identifier.

        """
        if id is not None:
            return id

        return self.default_background_id

    def get_background(self,
                       id: Optional[IdType] = None
                       ) -> Optional[DataPHA]:
        """Return the background component.

        Parameters
        ----------
        id : int, str, or None, optional
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
        bkg_id = self._fix_background_id(id)
        return self._backgrounds.get(bkg_id)

    def set_background(self,
                       bkg: DataPHA,
                       id: Optional[IdType] = None) -> None:
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
        id : int, str, or None, optional
           The identifier of the background component. If it is None
           then the default background identifier is used.

        See Also
        --------
        delete_background, get_background

        Notes
        -----
        If the PHA header does not have the TELESCOP, INSTRUME, or
        FILTER header keywords set (or they are set to "none"), then
        they are taken from the background, if they are not set to
        "none". This is to allow simulated data sets to be used with
        external programs, such as XSPEC.

        """

        if not isinstance(bkg, DataPHA):
            raise ArgumentTypeErr("badarg", "bkg", "a PHA data set")

        # Check that the background matches the source (i.e. self)
        # dataset.  For this use case we require that the source has
        # set up its channel array (ie self.channel cannot be None),
        # as allowing the background dataset to have channels but not
        # the source makes tracking the state harder and there is no
        # obvious use case for that.
        #
        if self.channel is None:
            raise DataErr("The channel field must be set before adding a background")

        if bkg.channel is None:
            raise DataErr("The channel field of the background must be set")

        if len(self.channel) != len(bkg.channel) or \
           np.any(self.channel != bkg.channel):
            raise DataErr("The source and background channels differ")

        bkg_id = self._fix_background_id(id)
        self._backgrounds[bkg_id] = bkg
        if bkg_id not in self.background_ids:
            self.background_ids.append(bkg_id)

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

        # To support simulated data (e.g. issue #1209) we over-write
        # the header TELESCOP/INSTRUME/FILTER settings to match the
        # background if they are at their "default" value (e.g. "none").
        #
        # TODO: should we warn if these values don't match (when not
        # "none")?
        #
        def set_key(key):
            val = self.header.get(key, "none")
            if val != "none":
                return

            val = bkg.header.get(key, "none")
            if val == "none":
                return

            self.header[key] = val

        set_key("TELESCOP")
        set_key("INSTRUME")
        set_key("FILTER")

    def delete_background(self,
                          id: Optional[IdType] = None
                          ) -> None:
        """Remove the background component.

        If the background component does not exist then the method
        does nothing.

        Parameters
        ----------
        id : int, str, or None, optional
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

        bkg_id = self._fix_background_id(id)
        if bkg_id not in self.background_ids:
            return

        self._backgrounds.pop(bkg_id, None)
        if len(self._backgrounds) == 0:
            self._subtracted = False

        self.background_ids.remove(bkg_id)

    def get_background_scale(self,
                             bkg_id: IdType = 1,
                             units='counts',
                             group=True,
                             filter=False
                             ):
        """Return the correction factor for the background dataset.

        .. versionchanged:: 4.12.2
           The bkg_id, units, group, and filter parameters have been
           added and the routine no-longer calculates the average
           scaling for all the background components but just for the
           given component.

        Parameters
        ----------
        bkg_id : int or str, optional
           The background component to use.
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
            raise ValueError(f"Invalid units argument: {units}")

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
        if np.isscalar(scale) and scale <= 0.0:
            return 1.0

        if np.iterable(scale):
            scale = np.asarray(scale, dtype=SherpaFloat)
            if group:
                if filter:
                    scale = self.apply_filter(scale, self._middle)
                else:
                    scale = self.apply_grouping(scale, self._middle)

            scale[scale <= 0.0] = 1.0

        return scale

    def get_backscal(self, group=True, filter=False):
        """Return the background scaling of the PHA data set.

        Return the BACKSCAL setting [OGIP_92_007]_ for the PHA data
        set.

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

        Examples
        --------

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> pha.get_backscal()
        2.5264364698914e-06

        """
        if self.backscal is None:
            return None

        return self._check_scale(self.backscal, group, filter)

    def get_areascal(self, group=True, filter=False):
        """Return the fractional area factor of the PHA data set.

        Return the AREASCAL setting [OGIP_92_007]_ for the PHA data
        set.

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

        Examples
        --------

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> pha.get_areascal()
        1.0

        """
        if self.areascal is None:
            return None

        return self._check_scale(self.areascal, group, filter)

    @overload
    def apply_filter(self,
                     data: None,
                     groupfunc: Callable) -> None:
        ...

    @overload
    def apply_filter(self,
                     data: ArrayType,
                     groupfunc: Callable) -> np.ndarray:
        ...

    def apply_filter(self, data, groupfunc=np.sum):
        """Group and filter the supplied data to match the data set.

        Parameters
        ----------
        data : ndarray or None
            The data to group, which must match either the number of
            channels of the data set or the number of filtered
            channels.
        groupfunc : function reference
            The grouping function. See apply_grouping for the
            supported values.

        Returns
        -------
        result : ndarray or None
            The grouped and filtered data, or None if the input was
            None.

        Raises
        ------
        sherpa.utils.err.DataErr
            If the data size does not match the number of channels.
        ValueError
            If the name of groupfunc is not supported or the data
            does not match the filtered data.

        See Also
        --------
        apply_grouping, ignore, ignore_bad, notice

        Examples
        --------

        Group and filter the counts array with no filter and then
        with a filter:

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> pha.grouped
        True
        >>> pha.notice()
        >>> print(pha.apply_filter(pha.counts))
        [17.  15.  16.  15. ...
        >>> pha.notice(0.5, 7)
        >>> print(pha.apply_filter(pha.counts))
        [15.  16.  15.  18.  ...

        As the previous example but with no grouping:

        >>> pha.ungroup()
        >>> pha.notice()
        >>> pha.apply_filter(pha.counts)[0:5]
        array([0., 0., 0., 0., 0.])
        >>> pha.notice(0.5, 7)
        >>> pha.apply_filter(pha.counts)[0:5]
        array([4., 3., 0., 1., 1.])

        Rather than group the counts, use the channel numbers and
        return the first and last channel number in each of the
        filtered groups (for the first five groups):

        >>> pha.group()
        >>> pha.notice(0.5, 7.0)
        >>> pha.apply_filter(pha.channel, pha._min)[0:5]
        array([33., 40., 45., 49., 52.])
        >>> pha.apply_filter(pha.channel, pha._max)[0:5]
        array([39., 44., 48., 51., 54.])

        Find the approximate energy range of each selected group from
        the RMF EBOUNDS extension:

        >>> rmf = pha.get_rmf()
        >>> elo = pha.apply_filter(rmf.e_min, pha._min)
        >>> ehi = pha.apply_filter(rmf.e_max, pha._max)

        Calculate the grouped data, after filtering, if the counts were
        increased by 2 per channel. Note that in this case the data to
        apply_filter contains the channel counts after applying the
        current filter:

        >>> orig = pha.counts[pha.get_mask()]
        >>> new = orig + 2
        >>> cts = pha.apply_filter(new)

        """
        if data is None:
            return None

        nelem = self.size
        if nelem is None:
            raise DataErr("sizenotset", self.name)

        data = _check(data)
        ndata = len(data)

        # We allow the data to have either (using the un-grouped data)
        #
        # - the size of the data object (all channels)
        # - the size of the filtered object
        #
        # This is unlike the other Data classes, where only the "all
        # channel" case is supported. We need to allow this new
        # behavior to support model evaluation via eval_model_to_fit
        # when using a PHA-based instrument model.
        #
        if ndata != nelem:

            mask = self.get_mask()
            if mask is None:
                # All elements have been filtered out, so error out.
                #
                raise DataErr("notmask")

            nfiltered = mask.sum()
            if nfiltered != ndata:
                # It is hard to get a concise error here: assume that the user
                # would prefer to know about the filtered length.
                #
                raise DataErr("mismatchn", "filtered data", "array",
                              nfiltered, ndata)

            # Create an array for the full channel range and insert
            # the user-values into it.
            #
            temp = np.zeros(nelem, dtype=SherpaFloat)
            temp[mask] = data
            data = temp

        gdata = self.apply_grouping(data, groupfunc)

        # We can not
        #
        #   return super().apply_filter(gdata)
        #
        # because the super-class does not know that gdata is a
        # grouped dataset, and so may be smaller than self.size.  The
        # problem is that the size attribute in the PHA case is poorly
        # defined: does it mean the number of channels or the number
        # of grouped channels? At the moment it means the former, and
        # there are places where we need this behavior. Can we add an
        # "effective size" property?
        #
        return self._data_space.filter.apply(gdata)

    @overload
    def apply_grouping(self,
                       data: None,
                       groupfunc: Callable) -> None:
        ...

    @overload
    def apply_grouping(self,
                       data: ArrayType,
                       groupfunc: Callable) -> np.ndarray:
        ...

    def apply_grouping(self, data, groupfunc=np.sum):
        """Apply the grouping scheme of the data set to the supplied data.

        Parameters
        ----------
        data : ndarray or None
            The data to group, which must match the number of channels
            of the data set.
        groupfunc : function reference
            The grouping function. Note that what matters is the name
            of the function, not its code. The supported function
            names are: "sum", "_sum_sq", "_min", "_max", "_middle",
            and "_make_groups".

        Returns
        -------
        grouped : ndarray or None
            The grouped data, unless the data set is not grouped or
            the input array was None, when the input data is returned.

        Raises
        ------
        sherpa.utils.err.DataErr
            If the data size does not match the number of channels.
        ValueError
            If the name of groupfunc is not supported.

        See Also
        --------
        apply_filter, ignore_bad

        Notes
        -----
        The supported grouping schemes are:

        ============ ======================================================
        Name         Description
        ============ ======================================================
        sum          Sum all the values in the group.
        _min         The minimum value in the group.
        _max         The maximum value in the group.
        _middle      The average of the minimum and maximum values.
        _sum_sq      The square root of the sum of the squared values.
        _make_groups The group number, starting at the first value of data.
        ============ ======================================================

        There are methods of the DataPHA class that can be used for
        all other than "sum" (the default value).

        The grouped data is not filtered unless a quality filter has
        been applied (e.g. by ignore_bad) in which case the quality
        filter will be applied to the result. In general apply_filter
        should be used if the data is to be filtered as well as
        grouped.

        Examples
        --------

        Sum up the counts in each group (note that the data has not
        been filtered so using get_dep with the filter argument set to
        True is generally preferred to using this method):

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> gcounts = pha.apply_grouping(pha.counts)

        The grouping for an unfiltered PHA data set with 1024 channels
        is used to calculate the number of channels in each group, the
        lowest channel number in each group, the highest channel
        number in each group, and the mid-point between the two:

        >>> pha.grouped
        True
        >>> pha.mask
        True
        >>> len(pha.channel)
        1024
        >>> import numpy as np
        >>> dvals = np.arange(1, 1025)
        >>> print(pha.apply_grouping(np.ones(1024)))
        [ 17.   4.  11.  ...
        >>> print(pha.apply_grouping(dvals, pha._min))
        [  1.  18.  22.  ...
        >>> print(pha.apply_grouping(dvals, pha._max))
        [  17.   21.   32.   ...
        >>> print(pha.apply_grouping(dvals, pha._middle))
        [  9.   19.5  27.   ...

        The grouped data is not filtered (unless ignore_bad has been
        used):

        >>> pha.notice()
        >>> v1 = pha.apply_grouping(dvals)
        >>> pha.notice(1.2, 4.5)
        >>> v2 = pha.apply_grouping(dvals)
        >>> np.all(v1 == v2)
        True

        """
        if data is None:
            return None

        nelem = self.size
        if nelem is None:
            raise DataErr("sizenotset", self.name)

        data = _check(data)
        ndata = len(data)
        if ndata != nelem:
            raise DataErr("mismatchn", "data", "array", nelem, ndata)

        # TODO: This should probably apply the quality filter whether
        # grouped or not.
        #
        if not self.grouped:
            return data

        groups = self.grouping
        filter = self.quality_filter
        if filter is None:
            return do_group(data, groups, groupfunc.__name__)

        nfilter = len(filter)
        if len(data) != nfilter or len(groups) != nfilter:
            raise DataErr("mismatchn", "quality filter", "array", nfilter, len(data))

        filtered_data = np.asarray(data)[filter]
        groups = np.asarray(groups)[filter]
        return do_group(filtered_data, groups, groupfunc.__name__)

    def ignore_bad(self) -> None:
        """Exclude channels marked as bad.

        Ignore any bin in the PHA data set which has a quality value
        that is not equal to zero.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set has no quality array.

        See Also
        --------
        ignore, notice

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

        qual_flags = ~np.asarray(self.quality, bool)

        if self.grouped and (self.mask is not True):
            self.notice()
            warning('filtering grouped data with quality flags,' +
                    ' previous filters deleted')

        elif not self.grouped:
            # if ungrouped, create/combine with self.mask
            if self.mask is not True:
                self.mask = self.mask & qual_flags
                return

            self.mask = qual_flags
            return

        # self.quality_filter used for pre-grouping filter
        self.quality_filter = qual_flags

    def _dynamic_group(self,
                       group_func: Union[Callable, str],
                       *args, **kwargs) -> None:
        """Group the data using the given function and arguments.

        In order to support the grouping module being optional this
        routine will error out if it is not loaded (each time it is
        called). To support Sherpa 4.14.0 and earlier group_func can
        be a callable, but it is expected to be a string which is the
        name of the callable from the "group" module (which in this
        module has been renamed to pygroup if it exists). It also allows
        the user the capability of sending in a callable that they have
        written without the need for the group library.

        If group_func is a callable then it must return the grouping
        and quality arrays for the new scheme.

        """

        if not callable(group_func):
            if not groupstatus:
                raise ImportErr('importfailed', 'group', 'dynamic grouping')

            # The assumption is that the symbol exists so it is
            # not worth catching the AttributeError if it does not,
            # because that's a programming error and would have been
            # caught in testing.
            #
            group_func = getattr(pygroup, group_func)

        keys = list(kwargs.keys())[:]
        for key in keys:
            if kwargs[key] is None:
                kwargs.pop(key)

        # This routine only makes sense to be called if channel and
        # counts are set. It is more to provide help to the type
        # checkers than anything.
        #
        if self.channel is None or self.counts is None:
            raise DataErr(f"data set '{self.name}' can not be grouped as channel or counts is not set")

        # If tabstops is given then we want to ensure it is an
        # ndarray.  Really this should be done on args as well, in
        # case the array is sent in as a positional argument, but we
        # always send it in as a keyword argument. An alternative is
        # to do the conversion in the C++ code, but that is
        # significantly harder to orchestrate so this approach has
        # been taken.
        #
        # We also want to ensure that it is the correct size. If a
        # user tries to call with tabStops=~self.mask, which is the
        # "obvious" thing to do, then the call will fail if the data
        # is already grouped, since mask will have less values in it
        # then required. However, we can identify this case and
        # convert the mask into a "per-channel" array (in the same way
        # that get_mask does it).
        #
        if "tabStops" in kwargs:
            ts = np.asarray(kwargs["tabStops"])

            # We only expand the array if it has the correct size
            # (expand_grouped_mask does not enforce length checks for
            # Sherpa ~ 4.15).
            #
            # TODO: this probably doesn't work if we have quality_filter
            # set.
            #
            nts = len(ts)
            nchan = len(self.channel)
            if self.grouped and nts != nchan and \
               np.iterable(self.mask) and len(self.mask) == nts:
                ts = expand_grouped_mask(ts, self.grouping)

            kwargs["tabStops"] = ts

        else:
            # If there is a mask, and it is an array, invert it for
            # the tabStops. Note that
            #
            # a) use get_mask to ensure we have a value for each channel
            # b) get_mask can return None or an array.
            #
            mask = self.get_mask()
            if np.iterable(mask):
                kwargs["tabStops"] = ~mask

        self.grouping, self.quality = group_func(*args, **kwargs)
        self.group()
        self._original_groups = False

    def group_bins(self, num, tabStops=None) -> None:
        """Group into a fixed number of bins.

        Combine the data so that there `num` equal-width bins (or
        groups). The binning scheme is, by default, applied to only
        the noticed data range. It is suggested that filtering is done
        before calling group_bins.

        .. versionchanged:: 4.16.0
           Grouping now defaults to only using the noticed channel
           range.

        Parameters
        ----------
        num : int
           The number of bins in the grouped data set. Each bin
           will contain the same number of channels.
        tabStops : array of int or bool, optional
           If not set then it will be based on the filtering of the
           data set, so that the grouping only uses the filtered
           data. If set, it should be an array of booleans where True
           indicates that the channel should not be used in the
           grouping (this array must match the number of channels in
           the data set).

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
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_bins(num, tabStops=tabStops)

        self._dynamic_group("grpNumBins", len(self.channel), num,
                            tabStops=tabStops)

    def group_width(self, val, tabStops=None) -> None:
        """Group into a fixed bin width.

        Combine the data so that each bin contains `num` channels.
        The binning scheme is, by default, applied to only the noticed
        data range. It is suggested that filtering is done before
        calling group_width.

        .. versionchanged:: 4.16.0
           Grouping now defaults to only using the noticed channel
           range.

        Parameters
        ----------
        val : int
           The number of channels to combine into a group.
        tabStops : array of int or bool, optional
           If not set then it will be based on the filtering of the
           data set, so that the grouping only uses the filtered
           data. If set, it should be an array of booleans where True
           indicates that the channel should not be used in the
           grouping (this array must match the number of channels in
           the data set).

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
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_width(val, tabStops=tabStops)

        self._dynamic_group("grpBinWidth", len(self.channel), val,
                            tabStops=tabStops)

    def group_counts(self, num, maxLength=None, tabStops=None) -> None:
        """Group into a minimum number of counts per bin.

        Combine the data so that each bin contains `num` or more
        counts. The background is *not* included in this calculation;
        the calculation is done on the raw data even if `subtract` has
        been called on this data set. The binning scheme is, by
        default, applied to only the noticed data range. It is
        suggested that filtering is done before calling group_counts.

        .. versionchanged:: 4.16.0
           Grouping now defaults to only using the noticed channel
           range.

        Parameters
        ----------
        num : int
           The number of channels to combine into a group.
        maxLength : int, optional
           The maximum number of channels that can be combined into a
           single group.
        tabStops : array of int or bool, optional
           If not set then it will be based on the filtering of the
           data set, so that the grouping only uses the filtered
           data. If set, it should be an array of booleans where True
           indicates that the channel should not be used in the
           grouping (this array must match the number of channels in
           the data set).

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

        Examples
        --------

        Group by 20 counts within the range 0.5 to 7 keV (this is
        the default behavior for 4.16 and later):

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> pha.set_analysis("energy")
        >>> pha.notice()
        >>> pha.notice(0.5, 7)
        >>> pha.group_counts(20)

        Group by 20 but over the whole channel range, but then
        filtering to the noticed range of 0.5 to 7 keV (this was the
        default behaviour before 4.16):

        >>> pha.group_counts(20, tabStops=[0] * pha.size)

        """
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_counts(num, maxLength=maxLength, tabStops=tabStops)

        self._dynamic_group("grpNumCounts", self.counts, num,
                            maxLength=maxLength, tabStops=tabStops)

    # DOC-TODO: see discussion in astro.ui.utils regarding errorCol
    def group_snr(self, snr, maxLength=None, tabStops=None, errorCol=None) -> None:
        """Group into a minimum signal-to-noise ratio.

        Combine the data so that each bin has a signal-to-noise ratio
        which exceeds `snr`. The background is *not* included in this
        calculation; the calculation is done on the raw data even if
        `subtract` has been called on this data set. The binning
        scheme is, by default, applied to only the noticed data
        range. It is suggested that filtering is done before calling
        group_snr.

        .. versionchanged:: 4.16.0
           Grouping now defaults to only using the noticed channel
           range.

        Parameters
        ----------
        snr : number
           The minimum signal-to-noise ratio that must be exceeded
           to form a group of channels.
        maxLength : int, optional
           The maximum number of channels that can be combined into a
           single group.
        tabStops : array of int or bool, optional
           If not set then it will be based on the filtering of the
           data set, so that the grouping only uses the filtered
           data. If set, it should be an array of booleans where True
           indicates that the channel should not be used in the
           grouping (this array must match the number of channels in
           the data set).
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
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_snr(snr, maxLength=maxLength, tabStops=tabStops,
                          errorCol=errorCol)

        self._dynamic_group("grpSnr", self.counts, snr,
                            maxLength=maxLength, tabStops=tabStops,
                            errorCol=errorCol)

    def group_adapt(self, minimum, maxLength=None, tabStops=None) -> None:
        """Adaptively group to a minimum number of counts.

        Combine the data so that each bin contains `num` or more
        counts. The difference to `group_counts` is that this
        algorithm starts with the bins with the largest signal, in
        order to avoid over-grouping bright features, rather than at
        the first channel of the data. The adaptive nature means that
        low-count regions between bright features may not end up in
        groups with the minimum number of counts. The binning scheme
        is, by default, applied to only the noticed data range. It is
        suggested that filtering is done before calling group_adapt.

        .. versionchanged:: 4.16.0
           Grouping now defaults to only using the noticed channel
           range.

        Parameters
        ----------
        minimum : int
           The number of channels to combine into a group.
        maxLength : int, optional
           The maximum number of channels that can be combined into a
           single group.
        tabStops : array of int or bool, optional
           If not set then it will be based on the filtering of the
           data set, so that the grouping only uses the filtered
           data. If set, it should be an array of booleans where True
           indicates that the channel should not be used in the
           grouping (this array must match the number of channels in
           the data set).

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
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_adapt(minimum, maxLength=maxLength,
                            tabStops=tabStops)

        self._dynamic_group("grpAdaptive", self.counts, minimum,
                            maxLength=maxLength, tabStops=tabStops)

    # DOC-TODO: see discussion in astro.ui.utils regarding errorCol
    def group_adapt_snr(self, minimum, maxLength=None, tabStops=None,
                        errorCol=None) -> None:
        """Adaptively group to a minimum signal-to-noise ratio.

        Combine the data so that each bin has a signal-to-noise ratio
        which exceeds `minimum`. The difference to `group_snr` is that
        this algorithm starts with the bins with the largest signal,
        in order to avoid over-grouping bright features, rather than
        at the first channel of the data. The adaptive nature means
        that low-count regions between bright features may not end up
        in groups with the minimum number of counts. The binning
        scheme is, by default, applied to only the noticed data
        range. It is suggested that filtering is done before calling
        group_adapt_snr.

        .. versionchanged:: 4.16.0
           Grouping now defaults to only using the noticed channel
           range.

        Parameters
        ----------
        minimum : number
           The minimum signal-to-noise ratio that must be exceeded
           to form a group of channels.
        maxLength : int, optional
           The maximum number of channels that can be combined into a
           single group.
        tabStops : array of int or bool, optional
           If not set then it will be based on the filtering of the
           data set, so that the grouping only uses the filtered
           data. If set, it should be an array of booleans where True
           indicates that the channel should not be used in the
           grouping (this array must match the number of channels in
           the data set).
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
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bkg.group_adapt_snr(minimum, maxLength=maxLength,
                                tabStops=tabStops, errorCol=errorCol)

        self._dynamic_group("grpAdaptiveSnr", self.counts, minimum,
                            maxLength=maxLength, tabStops=tabStops,
                            errorCol=errorCol)

    def eval_model_to_fit(self, modelfunc):
        model = super().eval_model_to_fit(modelfunc)
        return self.apply_filter(model)

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

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> bcounts = pha.sum_background_data()

        Calculate the scaling factor that you need to multiply the
        background data to match the source data. In this case the
        background data has been replaced by the value 1 (rather than
        the per-channel values used with the default argument):

        >>> bscale = pha.sum_background_data(lambda k, d: 1)

        """

        bdata_list = []

        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            bdata = get_bdata_func(bkg_id, bkg)

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

    def get_dep(self,
                filter: bool = False
                ) -> Optional[np.ndarray]:
        # FIXME: Aneta says we need to group *before* subtracting, but that
        # won't work (I think) when backscal is an array
        # if not self.subtracted:
        #     return self.counts
        # return self.counts - self.sum_background_data()

        dep = self.counts

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

        if bool_cast(filter):
            dep = self.apply_filter(dep)

        return dep

    # The code used to re-define set_dep, but the only difference
    # from the parent class was that it set the counts attribute and
    # not the y attribute. These are now the same so it is no-longer
    # needed.
    #
    # There was the following comment in the code which we keep here:
    #
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

    def get_staterror(self,
                      filter: bool = False,
                      staterrfunc: Optional[StatErrFunc] = None
                      ) -> Optional[np.ndarray]:
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

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi', use_errors=True)
        >>> dy = pha.get_staterror()

        Ensure that there is no pre-defined statistical-error column
        and then use the Chi2DataVar statistic to calculate the errors:

        >>> from sherpa.stats import Chi2DataVar
        >>> stat = Chi2DataVar()
        >>> pha.staterror = None
        >>> dy = pha.get_staterror(staterrfunc=stat.calc_staterror)

        """

        filter = bool_cast(filter)
        if filter:
            filterfunc = self.apply_filter
        else:
            filterfunc = self.apply_grouping

        # Use local definitions for routines that calculate the error
        # and scaling terms as they use the same logic for the source
        # and background datasets, including making it obvious that
        # the filtering and grouping comes from the "source" object
        # (i.e. the self object).
        #
        def get_error(dataobj):
            """Get the filtered and grouped error.

            The return value is the error array and then either None
            or the counts array that was used to calculate the errors.
            That is, the second term is only set when staterrfunc was
            called (this is only needed when the XSPEC variance
            routine is used but this is not checked for here).

            """

            if dataobj.staterror is None:
                if staterrfunc is None:
                    return None, None

                cnts = filterfunc(dataobj.counts)
                return staterrfunc(cnts), cnts

            return filterfunc(dataobj.staterror,
                              groupfunc=dataobj._sum_sq), None

        def get_scale(dataobj):
            """Calculate the scaling value.

            Returns EXPOSURE * BACKSCAL * AREASCAL (1.0 is used as a
            replacement for any value is not set).

            """

            scale = 1.0
            if dataobj.backscal is not None:
                scale *= self._check_scale(dataobj.backscal, filter=filter)

            if dataobj.areascal is not None:
                scale *= self._check_scale(dataobj.areascal, filter=filter)

            if dataobj.exposure is not None:
                scale *= dataobj.exposure

            return scale

        staterr, src_counts = get_error(self)
        if staterr is None:
            return None

        if not self.subtracted:
            return staterr

        src_scale = get_scale(self)

        # For each background dataset we filter and/or group the
        # errors to match the source dataset, and then apply the
        # various scaling factors (areascal, backscal, and exposure
        # scaling) to correct them to match the source dataset.  The
        # per-group variance is then stored, so that it can be
        # combined to create the overall background contribution,
        # which can then have the "source" side of the scaling values
        # applied to it before being added to the source error term.
        #
        # A complication is if the errors are being calculated with
        # the datavar method. That is, if the staterror attribute is
        # None and staterrfunc is using "XSPEC" errors - in this case
        # Chi2XspecVar.calc_staterror (technically it's really
        # anything that uses
        # sherpa.stats._statfcts.calc_chi2xspecvar_errors but do we
        # want to catch this as well as the more OO version). This
        # requires storing the scaling factors.
        #
        bkg_variances = []
        bkg_counts = []
        bkg_scales = []

        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            berr, bcounts = get_error(bkg)
            if berr is None:
                # We do not know how to generate an error, so
                # return None. An alternative would be to raise an
                # error, since we have an error from the source.
                #
                return None

            bscale = get_scale(bkg)

            # Apply the scaling factors to the variances
            #
            bkg_variances.append(1.0 * berr * berr / (bscale * bscale))

            # These only need to be stored for the XSPEC variance case.
            #
            bkg_scales.append(bscale)
            bkg_counts.append(bcounts)

        # For the variance case we need to handle things carefully
        # when
        #
        #  1. we are calculating errors rather than them being
        #     sent in (i.e. both src_counts and bkg_counts are
        #     not None);
        #
        #  2. and we are using the XSPEC variance method (we check
        #     only for the Chi2XspecVar method version and not the
        #     low-level
        #     sherpa.stats._statfcts.calc_chi2xspecvar_errors
        #     case, as users are not expected to be using this).
        #
        # The reason for point 2 is that calc_chi2xspecvar_errors
        # returns an error of 1 when the number of counts is 0.  This
        # is okay when there's no background subtraction, but once we
        # have to worry about the background bins it gets complicated.
        # For each bin we have src_i and bkg_i_j (for j=1 to N where N
        # is normally 1), and these refer to the number of counts that
        # are in this bin (i.e. that are then passed to staterrfct).
        #
        # If all the values for the bin are > 0 then we can combine
        # the variances (and so can ignore them here), otherwise we have
        # to recalculate the error terms.
        #
        if staterrfunc == Chi2XspecVar.calc_staterror and \
           src_counts is not None and all(n is not None for n in bkg_counts):
            replace_xspecvar_values(src_counts, bkg_counts,
                                    staterr, bkg_variances,
                                    src_scale, bkg_scales)

        nbkg = len(bkg_variances)
        if nbkg == 1:
            bkgvar = bkg_variances[0]
        else:
            bkgvar = sum(bkg_variances) / (nbkg * nbkg)

        statvar = staterr * staterr + bkgvar * src_scale * src_scale
        return np.sqrt(statvar)

    def get_syserror(self,
                     filter: bool = False
                     ) -> Optional[np.ndarray]:
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
        if bool_cast(filter):
            syserr = self.apply_filter(syserr, self._sum_sq)
        else:
            syserr = self.apply_grouping(syserr, self._sum_sq)

        return syserr

    def get_x(self,
              filter: bool = False,
              response_id: Optional[IdType] = None
              ) -> Optional[np.ndarray]:
        if self.channel is None:
            return None

        # We want the full channel grid with no grouping.
        #
        if self.units == "channel":
            return self.channel

        elo, ehi = self._get_ebins(response_id=response_id, group=False)
        emid = (elo + ehi) / 2

        if self.units == "energy":
            return emid

        # The units must be wavelength.
        #
        tiny = np.finfo(np.float32).tiny
        # In case there are any 0-energy bins replace them
        emid[emid == 0.0] = tiny
        return hc / emid

    def get_xlabel(self) -> str:
        """The label for the independent axis.

        If set_xlabel has been called then the label used in that call
        will be used, otherwise the label will change depending on the
        units setting.

        .. versionchanged:: 4.17.0
           The return value depends on if set_xlabel has been called.

        See Also
        --------
        set_analysis, set_xlabel

        """
        if self._xlabel is not None:
            return self._xlabel

        xlabel = self.units.capitalize()
        if self.units == 'energy':
            xlabel += ' (keV)'
        elif self.units == 'wavelength':
            xlabel += ' (Angstrom)'
        # elif self.units == 'channel' and self.grouped:
        #     xlabel = 'Group Number'

        return xlabel

    def _set_initial_quantity(self) -> None:
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

    @overload
    def _fix_y_units(self,
                     val: None,
                     filter: bool = False,
                     response_id: Optional[IdType] = None
                     ) -> None:
        ...

    @overload
    def _fix_y_units(self,
                     val: ArrayType,
                     filter: bool = False,
                     response_id: Optional[IdType] = None
                     ) -> np.ndarray:
        ...

    def _fix_y_units(self, val, filter=False, response_id=None):
        """Rescale the data to match the 'y' axis."""

        if val is None:
            return val

        filter = bool_cast(filter)
        # make a copy of data for units manipulation
        val = np.array(val, dtype=SherpaFloat)

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

        # Should this be a user-callable method?
        #
        def get_bin_edges():
            if self.units != 'channel':
                xlo, xhi = self._get_ebins(response_id, group=False)
            else:
                xlo, xhi = (self.channel, self.channel + 1.)

            if filter:
                # If we apply a filter, make sure that
                # ebins are ungrouped before applying
                # the filter.
                xlo = self.apply_filter(xlo, self._min)
                xhi = self.apply_filter(xhi, self._max)
            elif self.grouped:
                xlo = self.apply_grouping(xlo, self._min)
                xhi = self.apply_grouping(xhi, self._max)

            return xlo, xhi

        if self.grouped or self.rate:
            xlo, xhi = get_bin_edges()
            if self.units == 'wavelength':
                dx = hc / xlo - hc / xhi
            else:  # Could be "energy" or "channel"
                dx = xhi - xlo

            val /= np.abs(dx)

        # The final step is to multiply by the X axis self.plot_fac
        # times.
        if self.plot_fac <= 0:
            return val

        # Get the bin edges so we can calculate the center of
        # each bin. This used to use
        #
        # xvals = self.get_x(response_id=response_id)
        # if filter:
        #     xvals = self.apply_filter(xvals, self._middle)
        # elif self.grouped:
        #     xvals = self.apply_grouping(xvals, self._middle)
        #
        # but this uses the center of each channel and then
        # averages them, which doesn't quite match the following.
        #
        xlo, xhi = get_bin_edges()
        xmid = (xlo + xhi) / 2
        val *= np.power(xmid, self.plot_fac)
        return val

    @overload
    def get_y(self,
              filter: bool,
              yfunc: None,
              response_id: Optional[IdType] = None,
              use_evaluation_space: bool = False
              ) -> np.ndarray:
        ...

    @overload
    def get_y(self,
              filter: bool,
              yfunc: ModelFunc,
              response_id: Optional[IdType] = None,
              use_evaluation_space: bool = False
              ) -> tuple[np.ndarray, np.ndarray]:
        ...

    # Note that use_evaluation_space is unused.
    #
    def get_y(self, filter=False, yfunc=None, response_id=None,
              use_evaluation_space=False):

        vals = Data.get_y(self, yfunc=yfunc)
        vallist = (vals,) if yfunc is None else vals

        filter = bool_cast(filter)
        out = []
        for val in vallist:
            if filter:
                val = self.apply_filter(val)
            else:
                val = self.apply_grouping(val)
            val = self._fix_y_units(val, filter, response_id)
            out.append(val)

        if yfunc is None:
            return out[0]

        return tuple(out)

    def get_yerr(self,
                 filter: bool = False,
                 staterrfunc: Optional[Callable] = None,
                 response_id: Optional[IdType] = None
                 ) -> Optional[np.ndarray]:
        filter = bool_cast(filter)
        err = self.get_error(filter, staterrfunc)
        return self._fix_y_units(err, filter, response_id)

    def get_xerr(self,
                 filter: bool = False,
                 response_id: Optional[IdType] = None
                 ) -> Optional[np.ndarray]:
        """Returns an X "error".

        The error value for the independent axis is not well defined
        in Sherpa.

        .. versionchanged:: 4.16.1
           The return value is now half the bin width instead of the
           full bin width and is now calculated correctly when the
           analysis is set to "wavelength".

        Parameters
        ----------
        filter : bool, optional
           Should the values be filtered to the current notice range?
        response_id : int, str, or None, optional
           What response should be used?

        Returns
        -------
        xerr : ndarray
           The half-width of each bin (or group) in the current
           analysis units.

        """
        if bool_cast(filter):
            # If we apply a filter, make sure that
            # ebins are ungrouped before applying
            # the filter.
            elo, ehi = self._get_ebins(response_id, group=False)
            elo = self.apply_filter(elo, self._min)
            ehi = self.apply_filter(ehi, self._max)

        else:
            try:
                elo, ehi = self._get_ebins(response_id=response_id)
            except DataErr:
                # What should we do here? This indicates that all bins
                # have been marked as bad (and grouping is present).
                #
                return np.asarray([])

        # Issue #748 #1817 noted we should return the half-width
        # Issue #1985 notes that we need to support wavelength.
        #
        if self.units != "wavelength":
            return (ehi - elo) / 2

        dlam = hc / elo - hc / ehi
        return dlam / 2

    def get_ylabel(self, yfunc=None) -> str:
        """The label for the dependent axis.

        Parameters
        ----------
        yfunc
           Unused.

        See Also
        --------
        set_ylabel

        """

        if self._ylabel is not None:
            return self._ylabel

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
            from sherpa import plot
            latex = plot.backend.get_latex_for_string(
                f'^{self.plot_fac}')
            ylabel += f' X {self.units.capitalize()}{latex}'
        return ylabel

    @staticmethod
    # Dummy function to tell apply_grouping to construct
    # an array of groups.
    def _make_groups(array):
        pass

    @staticmethod
    def _middle(array):
        array = np.asarray(array)
        return (array.min() + array.max()) / 2.0

    @staticmethod
    def _min(array):
        array = np.asarray(array)
        return array.min()

    @staticmethod
    def _max(array):
        array = np.asarray(array)
        return array.max()

    @staticmethod
    def _sum_sq(array):
        return np.sqrt(np.sum(array * array))

    def get_noticed_channels(self) -> np.ndarray:
        """Return the noticed channels.

        .. versionchanged:: 4.17.0
           An empty array is now returned if the filter removes all
           channels and an error is now raised if the data is empty.

        Returns
        -------
        channels : ndarray
            The noticed channels (this is independent of the
            analysis setting).

        """

        if self.size is None:
            raise DataErr("sizenotset", self.name)

        chans = self.channel

        # get_mask returns None if all the data has been filtered.
        mask = self.get_mask()
        if mask is None:
            return np.asarray([], dtype=chans.dtype)

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

    def get_mask(self) -> Optional[np.ndarray]:
        """Returns the (ungrouped) mask.

        .. versionchanged:: 4.17.0
           When all channels are selected the routine now returns an
           array rather than None, and it is an error to call it on an
           empty data set (one whose channel field has not been set).

        Returns
        -------
        mask : ndarray or None
            The mask, in channels, or None.

        """

        if self.size is None:
            raise DataErr("sizenotset", self.name)

        # I am not convinced that self.mask will always be an array if
        # quality filtering is in use, so do not assume that for now.
        #
        if self.mask is False:
            return None

        if self.mask is True:
            if self.quality_filter is None:
                return np.ones(self.size, dtype=bool)

            return self.quality_filter.copy()

        if not self.grouped:
            mask = self.mask.copy()
            # This does not apply the quality filter
            return mask

        groups = self.grouping
        if self.quality_filter is not None:
            groups = groups[self.quality_filter]

        return expand_grouped_mask(self.mask, groups)

    def get_noticed_expr(self) -> str:
        """Returns the current set of noticed channels.

        The values returned are always in channels, no matter the
        current analysis setting.

        Returns
        -------
        expr : str
            The noticed channel range as a string of comma-separated
            "low-high" values. As these are channel filters the low
            and high values are inclusive. If all channels have been
            filtered out then "No noticed channels" is returned.

        See Also
        --------
        get_filter, get_noticed_channels

        """
        chans = self.get_noticed_channels()
        if len(chans) == 0:
            return 'No noticed channels'

        return create_expr(chans, format='%i')

    def get_filter(self,
                   group: bool = True,
                   format: str = '%.12f',
                   delim: str = ':'
                   ) -> str:
        """Return the data filter as a string.

        The filter expression depends on the analysis setting.

        .. versionchanged:: 4.17.0
           The routine now always returns the empty string if all the
           data has been filtered out, rather than sometimes returning
           the phrase "No noticed bins".

        .. versionchanged:: 4.14.0
           Prior to 4.14.0 the filter used the mid-point of the bin,
           not its low or high value.

        Parameters
        ----------
        group : bool, optional
            Should the filter reflect the grouped data?
        format : str, optional
            The formatting of the numeric values (this is ignored for
            channel units, which uses ``format="%i"``).
        delim : str, optional
            The string used to mark the low-to-high range.

        Returns
        -------
        expr : str
            The noticed channel range as a string of comma-separated
            ranges, where the low and high values are separated by
            the `delim` string. The units of the ranges are controlled
            by the analysis setting. If all bins have been
            filtered out then the empty string is returned.

        See Also
        --------
        get_noticed_channels, get_noticed_expr

        Examples
        --------
        For a Chandra non-grating dataset which has been grouped:

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> pha.set_analysis('energy')
        >>> pha.notice(0.5, 7)
        >>> pha.get_filter(format='%.4f')
        '0.4672:9.8696'
        >>> pha.set_analysis('channel')
        >>> pha.get_filter()
        '33:676'

        The filter expression shows the first selected channel to the
        last one, and so is independent of whether the data is grouped
        or not:

        >>> pha.set_analysis('energy')
        >>> pha.get_filter(format='%.4f')
        '0.4672:9.8696'
        >>> pha.get_filter(group=False, format='%.4f')
        '0.4672:9.8696'

        Although the `group` argument does not change the output of
        get_filter, the selected range does depend on whether the data
        is grouped or not (unless the groups align with the filter
        edges):

        >>> pha.ungroup()
        >>> pha.notice()
        >>> pha.notice(0.5, 7)
        >>> pha.get_filter(format='%.3f')
        '0.496:7.008'
        >>> pha.group()
        >>> pha.get_filter(format='%.3f')
        '0.467:9.870'

        >>> pha.notice()
        >>> pha.notice(0.5, 6)
        >>> pha.ignore(2.1, 2.2)
        >>> pha.get_filter(format='%.2f', delim='-')
        '0.47-2.09,2.28-6.57'

        >>> pha.ignore()
        >>> pha.get_filter()
        ''

        """

        # We use get_noticed_channels since it includes quality
        # filtering, which the 'self.mask is True' check below does
        # not make.
        #
        chans = self.get_noticed_channels()
        if len(chans) == 0:
            # This catches both self.mask=False and self.mask=[False] * nchan
            return ''

        # Special case all channels are selected.
        #
        if self.mask is True:
            elo, ehi = self._get_ebins(group=False, response_id=None)
            if self.units == 'energy':
                loval = elo[0]
                hival = ehi[-1]
            elif self.units == 'wavelength':
                loval = hc / ehi[-1]
                hival = hc / elo[0]
            else:
                # Assume channel
                loval = self.channel[0]
                hival = self.channel[-1]
                format = '%i'

            # Check for inversion
            if loval > hival:
                loval, hival = hival, loval

            return f"{format % loval}{delim}{format % hival}"

        mask = self.get_mask()
        # Just to check
        assert mask is not None

        # What channels are selected - note that it is possible to
        # get here and have all channels selected (ie all elements
        # of self.mask are True).
        #
        # We handle channel units differently to energy or wavelength,
        # as channels are treated as Data1D does whereas we want to
        # use Data1DInt for the units.
        #
        # We do everything ungrouped as there is no difference now
        # (as of Sherpa 4.14.0).
        #
        if self.units == 'channel':
            return create_expr(chans, mask=mask, format='%i', delim=delim)

        # Unfortunately we don't have a usable API for accessing the
        # energy or wavelength ranges directly.
        #
        xlo, xhi = self._get_ebins(group=False)
        if self.units == 'wavelength':
            xlo, xhi = hc / xhi, hc / xlo

        # Ensure the data is in ascending order for create_expr_integrated.
        #
        if xlo[0] > xlo[-1]:
            xlo = xlo[::-1]
            xhi = xhi[::-1]
            mask = mask[::-1]

        return create_expr_integrated(xlo[mask], xhi[mask], mask=mask, format=format, delim=delim)

    def get_filter_expr(self) -> str:
        return (self.get_filter(format='%.4f', delim='-') +
                ' ' + self.get_xlabel())

    def notice_response(self,
                        notice_resp: bool = True,
                        noticed_chans: Optional[np.ndarray] = None
                        ) -> None:
        notice_resp = bool_cast(notice_resp)

        if notice_resp and noticed_chans is None:
            noticed_chans = self.get_noticed_channels()

        for resp_id in self.response_ids:
            arf, rmf = self.get_response(resp_id)
            _notice_resp(noticed_chans, arf, rmf)

    def notice(self,
               lo: Optional[float] = None,
               hi: Optional[float] = None,
               ignore: bool = False,
               bkg_id: Optional[Union[IdType, Sequence[IdType]]] = None
               ) -> None:
        """Notice or ignore the given range.

        Select or remove a range of channels when analyzing the data.
        The units for the range, that is the values for the lo and hi
        argument, is set by the analysis field.

        .. versionchanged:: 4.14.0
           PHA filtering has been improved to fix a number of corner
           cases which can result in the same filter now selecting one
           or two fewer channels that done in earlier versions of
           Sherpa. The ``lo`` and ``hi`` arguments are now restricted
           based on the units setting.

        Parameters
        ----------
        lo, hi : number or None, optional
            The range to change. A value of None means the minimum or
            maximum permitted value. The units of ``lo`` and ``hi``
            are set by the units field.
        ignore : bool, optional
            Set to True if the range should be ignored. The default is
            to notice the range.
        bkg_id : int, str, None, or sequence of int or str, optional
            If not None then apply the filter to the given background
            dataset or datasets, otherwise change the object and all
            its background datasets.

        See Also
        --------
        get_filter, get_filter_expr, get_mask

        Notes
        -----
        Calling notice with no arguments selects all points in the
        dataset (or, if ``ignore=True``, it will remove all points).

        If no channels have been ignored then a call to `notice` with
        ``ignore=False`` will select just the ``lo`` to ``hi`` range,
        and exclude any channels outside this range. If there has been
        a filter applied then the range ``lo`` to ``hi`` will be added
        to the range of noticed data (when ``ignore=False``). One
        consequence to the above is that if the first call to `notice`
        (with ``ignore=False``) selects a range outside the data set -
        such as a channel range of 2000 to 3000 when the valid range
        is 1 to 1024 - then all points will be ignored.

        When filtering with channel units then:

        - the ``lo`` and ``hi`` arguments, if set, must be integers,
        - and the ``lo`` and ``hi`` values are inclusive.

        For energy and wavelength filters:

        - the ``lo`` and ``hi`` arguments, if set, must be >= 0,
        - and the ``lo`` limit is inclusive but the ``hi`` limit is exclusive.

        Examples
        --------
        So, for an ungrouped PHA file with 1024 channels:

        >>> from sherpa.astro.io import read_pha
        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> pha.ungroup()
        >>> pha.units = 'channel'
        >>> pha.get_filter()
        '1:1024'
        >>> pha.notice(20, 200)
        >>> pha.get_filter()
        '20:200'
        >>> pha.notice(300, 500)
        >>> pha.get_filter()
        '20:200,300:500'

        Calling `notice` with no arguments removes all the filters:

        >>> pha.notice()
        >>> pha.get_filter()
        '1:1024'

        Ignore the first 30 channels (this is the same as calling
        ```pha.ignore(hi=30)``:

        >>> pha.notice(hi=30, ignore=True)
        >>> pha.get_filter()
        '31:1024'

        When using wavelength or energy units the noticed (or ignored)
        range will not always match the requested range because each
        channel has a finite width in these spaces:

        >>> pha = read_pha(data_3c273 + '3c273.pi')
        >>> pha.grouped
        True
        >>> pha.get_analysis()
        'energy'
        >>> pha.notice()
        >>> pha.notice(0.5, 7)
        >>> pha.get_filter(format='%.3f')
        '0.467:9.870'

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

        # Validate input
        #
        if lo is not None and hi is not None and lo > hi:
            raise DataErr('bad', 'hi argument', 'must be >= lo')

        # Ensure the limits are physically meaningful, that is
        # energy and wavelengths are >= 0. Technically it should be
        # > but using 0 is a nice value for a minimum. We do not
        # enforce limits if channels are being used because it's
        # not clear if channels can technically be negative.
        #
        # For channels we just require the numbers are integers.
        #
        if self.units == 'channel':
            if lo is not None and not float(lo).is_integer():
                raise DataErr('bad', 'lo argument', 'must be an integer channel value')
            if hi is not None and not float(hi).is_integer():
                raise DataErr('bad', 'hi argument', 'must be an integer channel value')

        else:
            if lo is not None and lo < 0:
                raise DataErr('bad', 'lo argument', 'must be >= 0')
            if hi is not None and hi < 0:
                raise DataErr('bad', 'hi argument', 'must be >= 0')

        # If any background IDs are actually given, then impose
        # the filter on those backgrounds *only*, and return.  Do
        # *not* impose filter on data itself.
        #
        filter_background_only = False
        if bkg_id is not None:
            # As bkg_id can be
            #   - int or str
            #   - iterable of int or str
            # it's a bit awkward to identify what is meant.
            #
            if not isinstance(bkg_id, str) and np.iterable(bkg_id):
                bkg_ids = bkg_id
            else:
                bkg_ids = [bkg_id]

            filter_background_only = True
        else:
            bkg_ids = self.background_ids

        # Automatically impose data's filter on background data sets.
        # Units must agree for this to be meaningful, so temporarily
        # make data and background units match.
        #
        for bid in bkg_ids:
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
        #
        if filter_background_only:
            return

        # Go on if we are also supposed to filter the source data
        if lo is None and hi is None:
            self.quality_filter = None
            self.notice_response(False)

        # elo and ehi will be in channel (units=channel) or energy
        # (units=energy or units=wavelength).
        #
        try:
            elo, ehi = self._get_ebins(group=self.grouped)
        except DataErr as de:
            info("Skipping dataset %s: %s", self.name, str(de))
            return

        emin = min(elo[[0, -1]])
        emax = max(ehi[[0, -1]])

        # Convert wavelength limits to energy if necessary.
        #
        if self.units == 'wavelength':
            lims = validate_wavelength_limits(lo, hi, emax)
            if lims is None:
                # No useful filter to apply
                return

            lo, hi = lims

        # safety check
        assert lo is None or hi is None or lo <= hi, (lo, hi, self.name)

        if self.units == 'channel' and hi is not None:
            # A channel range lo to hi is read as [lo, hi] rather than
            # [lo, hi), so we increase the upper limit by 1 to
            # work around this, as the filter call checks for < hi
            # and not <= hi.
            #
            hi += 1

        self._data_space.filter.notice((None, lo), (hi, None),
                                       (elo, ehi), ignore=ignore,
                                       integrated=True)

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

    def group(self) -> None:
        """Group the data according to the data set's grouping scheme.

        This sets the grouping flag which means that the value of the
        grouping attribute will be used when accessing data values. This
        can be called even if the grouping attribute is empty.

        .. versionchanged:: 4.15.1
           The grouping status of any background component is now also
           changed.

        See Also
        --------
        ungroup

        """
        self.grouped = True

        # Ensure any backgrounds are also grouped.
        #
        for bkg_id in self.background_ids:
            bkg = self.get_background(bkg_id)
            try:
                bkg.grouped = True
            except DataErr as exc:
                info(str(exc))

    def ungroup(self) -> None:
        """Remove any data grouping.

        This un-sets the grouping flag which means that the grouping
        attribute will not be used when accessing data values.

        .. versionchanged:: 4.15.1
           The grouping status of any background component is now also
           changed.

        See Also
        --------
        group
        """
        self.grouped = False

        # Ensure any backgrounds are also grouped.
        #
        for bkg_id in self.background_ids:
            # Unlike the group case we do not need to worry about this
            # failing.
            bkg = self.get_background(bkg_id)
            bkg.grouped = False

    def subtract(self) -> None:
        """Subtract the background data.

        See Also
        --------
        unsubtract
        """
        self.subtracted = True

    def unsubtract(self) -> None:
        """Remove background subtraction.

        See Also
        --------
        subtract
        """
        self.subtracted = False


class DataIMG(Data2D):
    '''Image data set

    This class builds on `sherpa.data.Data2D` to add support for
    region filters and the ability to describe the data in different coordinate
    systems, such as the "logical" coordinates (pixels of the image),
    the "physical" coordinates (e.g. detector coordinates), and the "world"
    coordinates (e.g. Ra/Dec on the sky).
    While this class can also be used for sparse data, much of the functionality
    added over its parent class will not be useful in that case.

    Parameters
    ----------
    name : string
        name of this dataset
    x0 : array-like
        Independent coordinate values for the first dimension
    x1 : array-like
        Independent coordinate values for the second dimension
    y : array-like
        The values of the dependent observable.
    shape : tuple
        Shape of the data grid (optional). This is
        used return the data as an image e.g. for display,
        but it not needed for fitting and modelling.
    staterror : array-like
        the statistical error associated with the data
    syserror : array-like
        the systematic error associated with the data
    sky : `sherpa.astro.io.wcs.WCS`
        The optional WCS object that converts to the physical coordinate system.
    eqpos : `sherpa.astro.io.wcs.WCS`
        The optional WCS object that converts to the world coordinate system.
    coord : string
        The coordinate system of the data. This can be one of 'logical',
        'physical', or 'world'.
        'image' is an alias for 'logical' and 'wcs' is an alias for 'world'.
    header : dict
        The FITS header associated with the data to store meta data.

    Examples
    --------
    In this example, we create an image data set, but do not specify the
    coordinate system. This means that the coordinates are expressed in
    the logical coordinates of the image, i.e. in pixels::

    >>> from sherpa.astro.data import DataIMG
    >>> import numpy as np
    >>> x1, x0 = np.mgrid[20:30, 5:20]
    >>> datashape = x0.shape
    >>> y = np.sqrt((x0 - 10)**2 + (x1 - 31)**2)
    >>> x0 = x0.flatten()
    >>> x1 = x1.flatten()
    >>> y = y.flatten()
    >>> image = DataIMG("bimage", x0=x0, x1=x1, y=y, shape=datashape)

    Note that in this example, we end up with a "logical" coordinate system
    in ``image`` and no WCS system to convert it to anything else. On the other hand,
    in FITS standard terminology, the "logical" coordinate system is the
    "image", counting pixels starting at 1, while here the ``x0lo``` and ``x1lo``
    actually start at 20 and 5, respectively.
    This behavior works for now, but might be revisited.
    '''

    _extra_fields = ("sky", "eqpos", "coord")

    sky = None
    """The optional WCS object that converts to the physical coordinate system."""

    eqpos = None
    """The optional WCS object that converts to the world coordinate system."""

    @property
    def coord(self):
        """Return the coordinate setting.

        The attribute is one of 'logical', 'physical', or
        'world'. Use `set_coord` to change the setting.

        """
        return self._coord

    # We do not set this to @coord.setter as the attribute should be
    # changed with set_coord when outside the methods of this class.
    #
    def _set_coord(self, val):
        coord = str(val).strip().lower()

        if coord in ('logical', 'image'):
            coord = 'logical'

        elif coord == 'physical':
            self._check_physical_transform()

        elif coord in ('world', 'wcs'):
            self._check_world_transform()
            coord = 'world'

        else:
            raise DataErr('bad', 'coordinates', val)

        self._coord = coord

    def __init__(self, name, x0, x1, y, shape=None, staterror=None,
                 syserror=None, sky=None, eqpos=None, coord='logical',
                 header=None):
        self.name = name  # needed by transform checks
        self.sky = sky
        self.eqpos = eqpos
        self._set_coord(coord)
        self.header = {} if header is None else header
        self._region = None

        # Store the original axes so we can always recreate the other
        # systems without having to worry about numerical differences
        # from switching between systems. This is an explicit decision
        # to go for repeatable behavior at the expense of using more
        # memory. See #1380 for more information.
        #
        self._orig_indep_axis = (self.coord, x0, x1)

        super().__init__(name, x0, x1, y, shape, staterror, syserror)

        # special case the x-axis labels
        self._x0label = None
        self._x1label = None

    def _clear_filter(self):
        if self._region is None:
            return

        # Remove the spatial filter, as if the axis has changed we know
        # it is going to be wrong.
        #
        self.notice2d()
        warning("Region filter has been removed from '%s'",
                self.name)

    def _repr_html_(self) -> str:
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
        # (but correctly handling the "None" case, to avoid #1214)
        #
        if state['_region'] is not None:
            state['_region'] = state['_region'].__str__()
        return state

    def __setstate__(self, state):
        # Populate the function pointers we deleted at pickle time with
        # no-ops.
        # self.__dict__['_get_logical']=(lambda : None)
        # self.__dict__['_get_physical']=(lambda : None)
        # self.__dict__['_get_world']=(lambda : None)

        # Unfortunately we can't re-create the original data used to
        # create the object, but we can fake it using the last-selected
        # coordinate system (which could lead to some issues for the
        # world system - see #1380 - but there's little we can do here).
        # The two-step process is to get around the behavior of the
        # NoNewAttributesAfterInit parent class.
        #
        if "_orig_indep_axis" not in state:
            state["_orig_indep_axis"] = None

        self.__dict__.update(state)

        if 'header' not in state:
            self.header = {}

        if self._orig_indep_axis is None:
            self._orig_indep_axis = (self.coord, self.x0, self.x1)

        # This may check the data is correct, based on the coord setting,
        # but is it worth it? It may catch a case when data is loaded into
        # a system without WCS support.
        #
        self._set_coord(state['_coord'])

        # This used to always use the _region setting to create a
        # Region filter, but it doesn't make sense if the filter is
        # _None, so we now skip this case (it lead to #1214).
        #
        if self._region is None:
            return

        if regstatus:
            # Should we allow '' to be sent to Region?
            if self._region == '':
                self._region = Region()
            else:
                self._region = Region(self._region)
        else:
            # If the region is "" then str() will produce '' so we want
            # double quotes about it.
            warning('Unable to restore region="%s" as region module is not available.',
                    str(self._region))

            self._region = None

    def _check_physical_transform(self):
        """Error out if the physical transform does not exist"""
        if self.sky is None:
            raise DataErr('nocoord', self.name, 'physical')

    def _check_world_transform(self):
        """Error out if the world transform does not exist"""
        if self.eqpos is None:
            raise DataErr('nocoord', self.name, 'world')

    def _logical_to_physical(self, x0=None, x1=None):
        """Convert logical axes to physical"""
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_physical_transform()

        # logical -> physical
        x0, x1 = self.sky.apply(x0, x1)

        return (x0, x1)

    def _logical_to_world(self, x0=None, x1=None):
        """Convert logical axes to world"""
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
        """Convert physical axes to logical"""
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_physical_transform()

        # physical -> logical
        x0, x1 = self.sky.invert(x0, x1)

        return (x0, x1)

    def _physical_to_world(self, x0=None, x1=None):
        """Convert physical axes to world"""
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_world_transform()

        # physical -> world
        x0, x1 = self.eqpos.apply(x0, x1)

        return (x0, x1)

    def _world_to_logical(self, x0=None, x1=None):
        """Convert world axes to logical"""
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
        """Convert world axes to physical"""
        if x0 is None or x1 is None:
            x0, x1 = self.get_indep()

        self._check_shape()
        self._check_world_transform()

        # world -> physical
        x0, x1 = self.eqpos.invert(x0, x1)

        return (x0, x1)

    # Convert from the _orig_indep_axis tuple (coord, x0, x1) to the
    # required data system (if it isn't already set).
    #
    def _get_coordsys(self, coord):
        if self.coord == coord:
            return self.get_indep()

        (base, x0, x1) = self._orig_indep_axis
        x0 = x0.copy()
        x1 = x1.copy()
        if base == coord:
            return (x0, x1)

        conv = getattr(self, f'_{base}_to_{coord}')
        return conv(x0, x1)

    def get_logical(self):
        return self._get_coordsys("logical")

    def get_physical(self):
        return self._get_coordsys("physical")

    def get_world(self):
        return self._get_coordsys("world")

    # For compatibility with old Sherpa keywords
    get_image = get_logical
    get_wcs = get_world

    def set_coord(self, coord):
        """Change the `coord` attribute.

        .. versionchanged:: 4.14.1
           The filter created by `notice2d` is now cleared when the
           coordinate system is changed.

        Parameters
        ----------
        coord : {'logical', 'image', 'physical', 'world', 'wcs'}
            The coordinate system to use. Note that "image" is a
            synonym for "logical" and "wcs" is a synomyn for "world".

        """

        coord = str(coord).strip().lower()
        good = ('logical', 'image', 'physical', 'world', 'wcs')
        if coord not in good:
            raise DataErr('badchoices', 'coordinates', coord, ", ".join(good))

        if coord == 'wcs':
            coord = 'world'
        elif coord == 'image':
            coord = 'logical'

        if coord == self.coord:
            return

        func = getattr(self, f'get_{coord}')
        self.indep = func()
        self._set_coord(coord)

    def get_filter_expr(self):
        if self._region is None:
            return ""

        return str(self._region)

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
        reg = Region(val, os.path.isfile(val))

        # Calculate the mask for this region as an "included" region.
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
                # add an explicit "whole field" constructor to avoid
                # possible issues with stringification of multiple
                # ignores.
                reg = Region('field()').subtract(reg)

            self._region = reg
        elif ignore:
            self._region = self._region.subtract(reg)
        else:
            self._region = self._region.union(reg)

    def get_bounding_mask(self):
        mask = self.mask
        if not np.iterable(self.mask):
            return mask, None

        # create bounding box around noticed image regions
        mask = np.array(self.mask).reshape(*self.shape)

        x0_i, x1_i = np.where(mask)

        x0_lo = x0_i.min()
        x0_hi = x0_i.max()
        x1_lo = x1_i.min()
        x1_hi = x1_i.max()

        mask = mask[x0_lo:x0_hi + 1, x1_lo:x1_hi + 1]
        shape = mask.shape
        mask = mask.ravel()
        return mask, shape

    def get_img(self, yfunc=None):
        # FIXME add support for coords to image class -> DS9
        self._check_shape()
        y_img = self.filter_region(self.get_dep(False))
        y_img = y_img.reshape(*self.shape)
        if yfunc is None:
            return y_img

        m = self.eval_model_to_fit(yfunc)
        if np.iterable(self.mask):
            # if filtered, the calculated model must be padded up
            # to the data size to preserve img shape and WCS coord
            m = pad_bounding_box(m, self.mask)

        return (y_img,
                self.filter_region(m).reshape(*self.shape))

    def get_axes(self):
        # FIXME: how to filter an axis when self.mask is size of self.y?
        self._check_shape()

        # dummy placeholders needed b/c img shape may not be square!
        axis0 = np.arange(self.shape[1], dtype=float) + 1.
        axis1 = np.arange(self.shape[0], dtype=float) + 1.
        if self.coord == 'logical':
            return (axis0, axis1)

        dummy0 = np.ones(axis0.size, dtype=float)
        dummy1 = np.ones(axis1.size, dtype=float)

        if self.coord == 'physical':
            axis0, dummy = self._logical_to_physical(axis0, dummy0)
            dummy, axis1 = self._logical_to_physical(dummy1, axis1)

        else:
            axis0, dummy = self._logical_to_world(axis0, dummy0)
            dummy, axis1 = self._logical_to_world(dummy1, axis1)

        return (axis0, axis1)

    def get_x0label(self) -> str:
        """Return the label for the first independent axis.

        If set_x0label has been called then the label used in that call
        will be used, otherwise the label will change depending on the
        coordinate setting.

        .. versionchanged:: 4.17.0
           The return value depends on if set_x0label has been called.

        See Also
        --------
        get_x1label, set_coord, set_x0label

        """

        if self._x0label is not None:
            return self._x0label
        if self.coord == 'physical':
            return 'x0 (pixels)'
        if self.coord == 'world':
            return 'RA (deg)'
        return 'x0'

    def get_x1label(self) -> str:
        """Return the label for the second independent axis.

        If set_x1label has been called then the label used in that call
        will be used, otherwise the label will change depending on the
        coordinate setting.

        .. versionchanged:: 4.17.0
           The return value depends on if set_x1label has been called.

        See Also
        --------
        get_x0label, set_coord, set_x1label

        """

        if self._x1label is not None:
            return self._x1label
        if self.coord == 'physical':
            return 'x1 (pixels)'
        if self.coord == 'world':
            return 'DEC (deg)'

        return 'x1'

    def to_contour(self, yfunc=None):
        y = self.filter_region(self.get_dep(False))
        if yfunc is not None:
            m = self.eval_model_to_fit(yfunc)
            if np.iterable(self.mask):
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
        if data is None or not np.iterable(self.mask):
            return data

        # We do not want to change the data array hence the
        # explicit copy: this is done via astype to ensure
        # we convert to a type that can accept NaN values.
        #
        out = data.astype(dtype=SherpaFloat, casting="safe", copy=True)
        out[~self.mask] = np.nan
        return out


class DataIMGInt(DataIMG):
    '''Binned image data set

    This class is very similar to `sherpa.data.DataIMG`, but it stores integrated
    data, i.e. the values of the dependent observable are given between bin edges
    as opposed to at bin centers.

    This class builds on `sherpa.data.Data2DInt` to add support for
    region filters and the ability to describe the data in different coordinate
    systems, such as the "logical" coordinates (pixels of the image),
    the "physical" coordinates (e.g. detector coordinates), and the "world"
    coordinates (e.g. Ra/Dec on the sky).
    While this class can also be used for sparse data, much of the functionality
    added over its parent class will not be useful in that case.

    Parameters
    ----------
    name : string
        name of this dataset
    x0lo : array-like
        Lower bounds of the bins in the first dimension of the independent coordinate
    x1lo : array-like
        Lower bound of the bins in the second dimension of the independent coordinate
    x0hi : array-like
        Upper bound of the bins in the first dimension of the independent coordinate
    x1hi : array-like
        Upper bound of the bins in the second dimension of the independent coordinate
    y : array-like
        The values of the dependent observable.
    shape : tuple
        Shape of the data grid (optional). This is
        used return the data as an image e.g. for display,
        but it not needed for fitting and modelling.
    staterror : array-like
        the statistical error associated with the data
    syserror : array-like
        the systematic error associated with the data
    sky : `sherpa.astro.io.wcs.WCS`
        The optional WCS object that converts to the physical coordinate system.
    eqpos : `sherpa.astro.io.wcs.WCS`
        The optional WCS object that converts to the world coordinate system.
    coord : string
        The coordinate system of the data. This can be one of 'logical',
        'physical', or 'world'.
        'image' is an alias for 'logical' and 'wcs' is an alias for 'world'.
    header : dict
        The FITS header associated with the data to store meta data.

    Example
    -------
    In this example, we first generate a 1000 (x,y) points from a 2D Gaussian.
    This could be, e.g., photons observed from a star. In x direction, the center
    of the Gaussian is at 1.2 and in y direction at 0.0. We then use
    `~numpy.histogram2d` to bin the data into a 2D histogram. Note that numpy's
    convention for the ordering of the axis is not the same as Sherpa's.
    We want to use the bin edges when we create the `DataIMGInt` object, but numpy
    gives us just one set of edges. We therefore first pick either the upper or lower
    edges (``x0edges[:-1]`` or ``x0edges[1:]``) and then repeat them to match
    the number of bins in the flattened data array::

      >>> import numpy as np
      >>> x = np.random.normal(size=1000, loc=1.2)
      >>> y = np.random.normal(size=1000)
      >>> xrange = np.arange(-2, 4.1, 0.5)
      >>> yrange = np.arange(-2, 2.1, 0.5)
      >>> hist, x0edges, x1edges = np.histogram2d(y, x, bins=(yrange, xrange))
      >>> x0lo, x1lo = np.meshgrid(x0edges[:-1], x1edges[:-1])
      >>> x0hi, x1hi = np.meshgrid(x0edges[1:], x1edges[1:])
      >>> image = DataIMGInt("binned_image",
      ...                    x0lo=x0lo.flatten(), x1lo=x1lo.flatten(),
      ...                    x0hi=x0hi.flatten(), x1hi=x1hi.flatten(),
      ...                    y=hist.flatten(), shape=hist.shape)

    Note that in this example, we end up with a "logical" coordinate system
    in ``image`` and no WCS system to convert it to anything else. On the other hand,
    in FITS standard terminology, the "logical" coordinate system is the
    "image", counting pixels starting at 1, while here the ``x0lo``` and ``x1lo``
    actually start at -2. This behavior works for now, but might be revisited.
    '''

    def __init__(self, name, x0lo, x1lo, x0hi, x1hi, y, shape=None,
                 staterror=None, syserror=None, sky=None, eqpos=None,
                 coord='logical', header=None):

        # Note: we do not call the superclass here.
        self._region = None
        self.sky = sky
        self.eqpos = eqpos
        self._set_coord(coord)
        self.header = {} if header is None else header
        self.shape = shape

        self._x0label = 'x0'
        self._x1label = 'x1'
        self._ylabel = 'y'

        Data.__init__(self, name, (x0lo, x1lo, x0hi, x1hi), y, staterror, syserror)

    def _init_data_space(self, filter, *data):
        ndata = len(data)
        if ndata != 4:
            raise DataErr("wrongaxiscount", self.name, 4, ndata)

        ds = IntegratedDataSpace2D(filter, *data)
        self._check_data_space(ds)
        return ds

    def get_x0(self, filter=False):
        if self.size is None:
            return None
        indep = self._data_space.get(filter)
        return (indep.x0lo + indep.x0hi) / 2.0

    def get_x1(self, filter=False):
        if self.size is None:
            return None
        indep = self._data_space.get(filter)
        return (indep.x1lo + indep.x1hi) / 2.0

    @property
    def x0lo(self):
        """
        Property kept for compatibility
        """
        return self._data_space.x0lo

    @property
    def x0hi(self):
        """
        Property kept for compatibility
        """
        return self._data_space.x0hi

    @property
    def x1lo(self):
        """
        Property kept for compatibility
        """
        return self._data_space.x1lo

    @property
    def x1hi(self):
        """
        Property kept for compatibility
        """
        return self._data_space.x1hi

    def get_logical(self):
        coord = self.coord
        x0lo, x1lo, x0hi, x1hi = self.get_indep()
        if coord == 'logical':
            return (x0lo, x1lo, x0hi, x1hi)

        x0lo = x0lo.copy()
        x1lo = x1lo.copy()
        convert = getattr(self, f'_{coord}_to_logical')
        x0lo, x1lo = convert(x0lo, x1lo)

        x0hi = x0hi.copy()
        x1hi = x1hi.copy()
        x0hi, x1hi = convert(x0hi, x1hi)

        return (x0lo, x1lo, x0hi, x1hi)

    def get_physical(self):
        coord = self.coord
        x0lo, x1lo, x0hi, x1hi = self.get_indep()
        if coord == 'physical':
            return (x0lo, x1lo, x0hi, x1hi)

        x0lo = x0lo.copy()
        x1lo = x1lo.copy()
        convert = getattr(self, f'_{coord}_to_physical')
        x0lo, x1lo = convert(x0lo, x1lo)

        x0hi = x0hi.copy()
        x1hi = x1hi.copy()
        x0hi, x1hi = convert(x0hi, x1hi)

        return (x0lo, x1lo, x0hi, x1hi)

    def get_world(self):
        coord = self.coord
        x0lo, x1lo, x0hi, x1hi = self.get_indep()
        if coord == 'world':
            return (x0lo, x1lo, x0hi, x1hi)

        x0lo = x0lo.copy()
        x1lo = x1lo.copy()
        convert = getattr(self, f'_{coord}_to_world')
        x0lo, x1lo = convert(x0lo, x1lo)

        x0hi = x0hi.copy()
        x1hi = x1hi.copy()
        x0hi, x1hi = convert(x0hi, x1hi)

        return (x0lo, x1lo, x0hi, x1hi)

    def get_axes(self):
        # FIXME: how to filter an axis when self.mask is size of self.y?
        self._check_shape()

        # dummy placeholders needed b/c img shape may not be square!
        axis0lo = np.arange(self.shape[1], dtype=float) - 0.5
        axis1lo = np.arange(self.shape[0], dtype=float) - 0.5

        axis0hi = np.arange(self.shape[1], dtype=float) + 0.5
        axis1hi = np.arange(self.shape[0], dtype=float) + 0.5

        if self.coord == 'logical':
            return (axis0lo, axis1lo, axis0hi, axis1hi)

        dummy0 = np.ones(axis0lo.size, dtype=float)
        dummy1 = np.ones(axis1lo.size, dtype=float)

        if self.coord == 'physical':
            axis0lo, dummy = self._logical_to_physical(axis0lo, dummy0)
            axis0hi, dummy = self._logical_to_physical(axis0hi, dummy0)

            dummy, axis1lo = self._logical_to_physical(dummy1, axis1lo)
            dummy, axis1hi = self._logical_to_physical(dummy1, axis1hi)

        else:
            axis0lo, dummy = self._logical_to_world(axis0lo, dummy0)
            axis0hi, dummy = self._logical_to_world(axis0hi, dummy0)

            dummy, axis1lo = self._logical_to_world(dummy1, axis1lo)
            dummy, axis1hi = self._logical_to_world(dummy1, axis1hi)

        return (axis0lo, axis1lo, axis0hi, axis1hi)
