# 
#  Copyright (C) 2010  Smithsonian Astrophysical Observatory
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

import os
import logging
import numpy
from itertools import izip
import sherpa.ui.utils
from sherpa.ui.utils import _check_type, _send_to_pager
from sherpa.utils import SherpaInt, SherpaFloat, sao_arange
from sherpa.utils.err import *
from sherpa.data import Data1D
import sherpa.astro.all
import sherpa.astro.plot

from sherpa.ui.utils import loggable

warning = logging.getLogger(__name__).warning
info = logging.getLogger(__name__).info


__all__ = ('Session')


class Session(sherpa.ui.utils.Session):

    ###########################################################################
    # Standard methods
    ###########################################################################


    def __init__(self):
        self._pileup_models     = {}
        self._background_models = {}
        self._background_sources = {}

        self._astrosourceplot= sherpa.astro.plot.SourcePlot()
        self._astrocompsrcplot = sherpa.astro.plot.ComponentSourcePlot()
        self._astrocompmdlplot = sherpa.astro.plot.ComponentModelPlot()
        self._modelhisto = sherpa.astro.plot.ModelHistogram()
        self._bkgmodelhisto = sherpa.astro.plot.BkgModelHistogram()

        self._bkgdataplot = sherpa.astro.plot.BkgDataPlot()
        self._bkgmodelplot = sherpa.astro.plot.BkgModelPlot()
        self._bkgfitplot = sherpa.astro.plot.BkgFitPlot()
        self._bkgchisqrplot = sherpa.astro.plot.BkgChisqrPlot()
        self._bkgdelchiplot = sherpa.astro.plot.BkgDelchiPlot()
        self._bkgresidplot = sherpa.astro.plot.BkgResidPlot()
        self._bkgratioplot = sherpa.astro.plot.BkgRatioPlot()
        self._bkgsourceplot = sherpa.astro.plot.BkgSourcePlot()
        self._arfplot = sherpa.astro.plot.ARFPlot()
        self._orderplot = sherpa.astro.plot.OrderPlot()
        self._energyfluxplot = sherpa.astro.plot.EnergyFluxHistogram()
        self._photonfluxplot = sherpa.astro.plot.PhotonFluxHistogram()

        # This is a new dictionary of XSPEC module settings.  It
        # is meant only to be populated by the save function, so
        # that the user's XSPEC settings can be saved in the pickle
        # file.  Then, restore can peel out settings from the
        # restored _xspec_state variable, and set abundance,
        # cross-section, etc. in the XSPEC module.
        self._xspec_state = None
        
        sherpa.ui.utils.Session.__init__(self)

        self._pyblocxs = sherpa.astro.sim.MCMC()

        self._plot_types['order']=self._orderplot
        self._plot_types['energy']=self._energyfluxplot
        self._plot_types['photon']=self._photonfluxplot
        self._plot_types['astrocompsource']=self._astrocompsrcplot
        self._plot_types['astrocompmodel']=self._astrocompmdlplot

        self._plot_types['astrosource']=self._astrosourceplot
        self._plot_types['astromodel']=self._modelhisto
        self._plot_types['arf']=self._arfplot
        self._plot_types['bkg']=self._bkgdataplot
        self._plot_types['bkgmodel']=self._bkgmodelhisto
        self._plot_types['bkgfit']=self._bkgfitplot
        self._plot_types['bkgsource']=self._bkgsourceplot
        self._plot_types['bkgratio']=self._bkgratioplot
        self._plot_types['bkgresid']=self._bkgresidplot
        self._plot_types['bkgdelchi']=self._bkgdelchiplot
        self._plot_types['bkgchisqr']=self._bkgchisqrplot


    ###########################################################################
    # High-level utilities
    ###########################################################################

    def __setstate__(self, state):
        if not state.has_key('_background_sources'):
            self.__dict__['_background_sources'] = state.pop('_background_models')

        sherpa.ui.utils.Session.__setstate__(self, state)


    ### Ahelp ingest: 2015-04-27 DJB
    def clean(self):
        """Clear out the current Sherpa session.

        The `clean` function removes all data sets and model
        assignments, and restores the default settings for the
        optimisation and fit statistic.

        See Also
        --------
        save : Save the current Sherpa session to a file.
        restore : Load in a Sherpa session from a file.
        save_all : Save the Sherpa session as an ASCII file.

        Examples
        --------

        >>> clean()

        """
        self._pileup_models     = {}
        self._background_models = {}
        self._background_sources = {}

        self._astrosourceplot= sherpa.astro.plot.SourcePlot()
        self._astrocompsrcplot = sherpa.astro.plot.ComponentSourcePlot()
        self._astrocompmdlplot = sherpa.astro.plot.ComponentModelPlot()
        self._modelhisto = sherpa.astro.plot.ModelHistogram()
        self._bkgmodelhisto = sherpa.astro.plot.BkgModelHistogram()

        self._bkgdataplot = sherpa.astro.plot.BkgDataPlot()
        self._bkgmodelplot = sherpa.astro.plot.BkgModelPlot()
        self._bkgfitplot = sherpa.astro.plot.BkgFitPlot()
        self._bkgchisqrplot = sherpa.astro.plot.BkgChisqrPlot()
        self._bkgdelchiplot = sherpa.astro.plot.BkgDelchiPlot()
        self._bkgresidplot = sherpa.astro.plot.BkgResidPlot()
        self._bkgratioplot = sherpa.astro.plot.BkgRatioPlot()
        self._bkgsourceplot = sherpa.astro.plot.BkgSourcePlot()
        self._arfplot = sherpa.astro.plot.ARFPlot()
        self._orderplot = sherpa.astro.plot.OrderPlot()
        self._energyfluxplot = sherpa.astro.plot.EnergyFluxHistogram()
        self._photonfluxplot = sherpa.astro.plot.PhotonFluxHistogram()

        self._xspec_state = None
        
        sherpa.ui.utils.Session.clean(self)

        self._pyblocxs = sherpa.astro.sim.MCMC()

        self._plot_types['order']=self._orderplot
        self._plot_types['energy']=self._energyfluxplot
        self._plot_types['photon']=self._photonfluxplot
        self._plot_types['astrocompsource']=self._astrocompsrcplot
        self._plot_types['astrocompmodel']=self._astrocompmdlplot

        self._plot_types['astrosource']=self._astrosourceplot
        self._plot_types['astromodel']=self._modelhisto
        self._plot_types['arf']=self._arfplot
        self._plot_types['bkg']=self._bkgdataplot
        self._plot_types['bkgmodel']=self._bkgmodelhisto
        self._plot_types['bkgfit']=self._bkgfitplot
        self._plot_types['bkgsource']=self._bkgsourceplot
        self._plot_types['bkgratio']=self._bkgratioplot
        self._plot_types['bkgresid']=self._bkgresidplot
        self._plot_types['bkgdelchi']=self._bkgdelchiplot
        self._plot_types['bkgchisqr']=self._bkgchisqrplot

    ### Ahelp ingest: 2015-04-27 DJB
    # Add ability to save attributes sepcific to the astro package.
    # Save XSPEC module settings that need to be restored.
    def save(self, filename='sherpa.save', clobber=False):
        """Save the current Sherpa session to a file.

        Parameters
        ----------
        filename : str, optional
           The name of the file to write the results to. The default
           is `sherpa.save`.
        clobber : bool, optional
           This flag controls whether an existing file can be
           overwritten (`True`) or if it raises an exception (`False`,
           the default setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        clean : Clear all stored session data.
        restore : Load in a Sherpa session from a file.
        save_all : Save the Sherpa session as an ASCII file.

        Notes
        -----

        The current Sherpa session is saved using the Python `pickle`
        module. The output is a binary file, which may not be portable
        between versions of Sherpa, but is platform independent, and
        contains all the data. This means that files created by `save`
        can be sent to collaborators to share results.

        Examples
        --------

        Save the current session to the file `sherpa.save`.

        >>> save()

        Save the current session to the file `bestfit.sherpa`,
        overwriting any existing version of the file.

        >>> save('bestfit.sherpa', clobber=True)

        """
        if (hasattr(sherpa.astro, "xspec")):
            self._xspec_state = sherpa.astro.xspec.get_xsstate()
        else:
            self._xspec_state = None
        sherpa.ui.utils.Session.save(self, filename, clobber)

    ### Ahelp ingest: 2015-04-27 DJB
    def restore(self, filename='sherpa.save'):
        """Load in a Sherpa session from a file.

        Parameters
        ----------
        filename : str, optional
           The name of the file to read the results from. The default
           is `sherpa.save`.

        Raises
        ------
        IOError
           If `filename` does not exist.

        See Also
        --------
        clean : Clear all stored session data.
        save : Save the current Sherpa session to a file.

        Notes
        -----
        The input to `restore` must have been created with the `save`
        command. This is a binary file, which may not be portable
        between versions of Sherpa, but is platform independent. A
        warning message may be created if a file saved by an older
        (or newer) version of Sherpa is loaded. An example of such
        a message is

        WARNING: Could not determine whether the model is discrete.
        This probably means that you have restored a session saved with a previous version of Sherpa.
        Falling back to assuming that the model is continuous.

        Examples
        --------

        Load in the Sherpa session from `sherpa.save`.

        >>> restore()

        Load in the session from the given file:

        >>> restore('/data/m31/setup.sherpa')

        """
        sherpa.ui.utils.Session.restore(self, filename)
        if (hasattr(sherpa.astro, "xspec")):
            if (self._xspec_state != None):
                sherpa.astro.xspec.set_xsstate(self._xspec_state)
                self._xspec_state = None
    
    def _get_show_data(self, id=None):
        data_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            data = self.get_data(id)

            data_str += 'Data Set: %s\n' % id
            data_str += 'Filter: %s\n' % data.get_filter_expr()
            if isinstance(data, sherpa.astro.data.DataPHA):

                scale = data.get_background_scale()
                if scale is not None and numpy.isscalar(scale):
                    data_str += 'Bkg Scale: %g\n' % float(scale)

                data_str += 'Noticed Channels: %s\n' % data.get_noticed_expr()

            data_str += data.__str__() + '\n\n'

            if isinstance(data, sherpa.astro.data.DataPHA):
                for resp_id in data.response_ids:
                    # ARF or RMF could be None
                    arf, rmf = data.get_response(resp_id)
                    if rmf is not None:
                        data_str += 'RMF Data Set: %s:%s\n' % (id,resp_id)
                        data_str += rmf.__str__() + '\n\n'
                    if arf is not None:
                        data_str += 'ARF Data Set: %s:%s\n' % (id,resp_id)
                        data_str += arf.__str__() + '\n\n'

                data_str += self._get_show_bkg(id)

        return data_str


    def _get_show_bkg(self, id=None, bkg_id=None):
        data_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]

        for id in ids:
            data = self.get_data(id)

            if not isinstance(data, sherpa.astro.data.DataPHA):
                continue

            bkg_ids = data.background_ids
            if bkg_id is not None:
                bkg_ids = [data._fix_background_id(bkg_id)]

            for bkg_id in bkg_ids:
                    bkg = self.get_bkg(id, bkg_id)
                    data_str += 'Background Data Set: %s:%s\n' % (id,bkg_id)
                    data_str += 'Filter: %s\n' % bkg.get_filter_expr()
                    data_str += 'Noticed Channels: %s\n' % bkg.get_noticed_expr()
                    data_str += bkg.__str__() + '\n\n'

                    for bk_rp_id in bkg.response_ids:
                        # ARF or RMF could be None
                        arf, rmf = bkg.get_response(bk_rp_id)
                        if rmf is not None:
                            data_str += ('Background RMF Data Set: %s:%s\n' %
                                         (id,bkg_id))
                            data_str += rmf.__str__() + '\n\n'
                        if arf is not None:
                            data_str += ('Background ARF Data Set: %s:%s\n' %
                                         (id,bkg_id))
                            data_str += arf.__str__() + '\n\n'

        return data_str


    def _get_show_bkg_model(self, id=None, bkg_id=None):
        model_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            if bkg_id is not None:
                bkg_ids = [bkg_id]
            else:
                bkg_ids = self._background_models.get(id, {}).keys()
                bkg_ids.extend(self._background_sources.get(id, {}).keys())
                bkg_ids = list(set(bkg_ids))

            for bkg_id in bkg_ids:
                model_str += 'Background Model: %s:%s\n' % (id, bkg_id)
                model_str += self.get_bkg_model(id, bkg_id).__str__() + '\n\n'

        return model_str


    def _get_show_bkg_source(self, id=None, bkg_id=None):
        model_str = ''
        ids = self.list_data_ids()
        if id is not None:
            ids = [self._fix_id(id)]
        for id in ids:
            if bkg_id is not None:
                bkg_ids = [bkg_id]
            else:
                bkg_ids = self._background_sources.get(id, {}).keys()

            for bkg_id in bkg_ids:
                model_str += 'Background Source: %s:%s\n' % (id, bkg_id)
                model_str += self.get_bkg_source(id, bkg_id).__str__() + '\n\n'

        return model_str


    ### Ahelp ingest: 2015-05-02 DJB
    def show_bkg(self, id=None, bkg_id=None, outfile=None, clobber=False):
        """Show the details of the PHA background data sets.

        This displays information about the background, or
        backgrounds, for the loaded data sets. This includes: any
        filters, the grouping settings, mission-specific header
        keywords, and the details of any associated instrument
        responses files (ARF, RMF).

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then all background data sets
           are displayed.
        bkg_id : int or str, optional
           The background component to display. The default is all
           components.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        list_model_ids : List of all the data sets with a source expression.
        load_bkg : Load the background from a file and add it to a PHA data set.
        show_all : Report the current state of the Sherpa session.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_bkg(id, bkg_id)
        _send_to_pager(all, outfile, clobber)


    ### Ahelp ingest: 2015-05-02 DJB
    def show_bkg_source(self, id=None, bkg_id=None, outfile=None, clobber=False):
        """Display the background model expression for a data set.

        This displays the background model for a data set, that is,
        the expression set by `set_bkg_model` or `set_bkg_source`, as
        well as the parameter values for the model. The
        `show_bkg_model` function displays the model that is fit to
        the data; that is, it includes any instrument responses.

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then all background expressions
           are displayed.
        bkg_id : int or str, optional
           The background component to display. The default is all
           components.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        list_model_ids : List of all the data sets with a source expression.
        set_bkg_model : Set the background model expression for a data set.
        show_all : Report the current state of the Sherpa session.
        show_model : Display the model expression used to fit a data set.
        show_bkg_model : Display the background model expression used to fit a data set.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_bkg_source(id, bkg_id)
        _send_to_pager(all, outfile, clobber)


    ### Ahelp ingest: 2015-05-01 DJB
    def show_bkg_model(self, id=None, bkg_id=None, outfile=None, clobber=False):
        """Display the background model expression used to fit a data set.

        This displays the model used to the the background data set,
        that is, the expression set by `set_bkg_model` or
        `set_bkg_source` combined with any instrumental responses,
        together with the parameter values for the model. The
        `show_bkg_source` function displays just the background model,
        without the instrument components (if any).

        Parameters
        ----------
        id : int or str, optional
           The data set. If not given then all background expressions are
           displayed.
        bkg_id : int or str, optional
           The background component to display. The default is all
           components.
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        list_model_ids : List of all the data sets with a source expression.
        set_bkg_model : Set the background model expression for a data set.
        show_all : Report the current state of the Sherpa session.
        show_model : Display the model expression used to fit a data set.
        show_bkg_source : Display the background model expression for a data set.

        Notes
        -----
        When `outfile` is `None`, the text is displayed via an external
        program to support paging of the information. The program
        used is determined by the `PAGER` environment variable. If
        `PAGER` is not found then '/usr/bin/more' is used.

        """
        all = ''
        all += self._get_show_bkg_model(id, bkg_id)
        _send_to_pager(all, outfile, clobber)

    ###########################################################################
    # Data
    ###########################################################################

    # DOC-NOTE: also in sherpa.utils
    #@loggable(with_id=True, with_name='load_data')
    def dataspace1d(self, start, stop, step=1, numbins=None,
                    id=None, bkg_id=None, dstype=sherpa.data.Data1DInt):
        """Create the independent axis for a 1D data set.

        dataspace1d

        SYNOPSIS
           Populates a blank 1D Sherpa data set by id

        SYNTAX

        Arguments:
           start   -  lower bound of grid

           stop    -  upper bound of grid

           step    -  bin width size
                      default is 1

           numbins -  number of bins desired
                      default is None

           id      -  Sherpa data id
                      defaut is default data id

           bkg_id  -  Sherpa background id
                      defaut is default background id

           dstype  -  Type of data set to use
                      default is Data1DInt

        Returns:
           None

        DESCRIPTION
           Populates a blank 1D Sherpa data set with the specified grid
           by Sherpa data id.  Alternatively, populate a blank PHA background
           data set by bkg_id.
           
           Specifying a dataspace using step size:
           if numbins is None (default) -> numpy.arange(start,stop,step)

           Specifying a dataspace by indicating the number of bins:
           if numbins is not None -> numpy.linspace(start, stop, numbins)

        EXAMPLES
           Blank integrated data set
           
              dataspace1d(0.1,10,0.1)

           Blank non-integrated data set

              dataspace1d(0.1,10,0.1,dstype=Data1D)

           Blank PHA data set

              dataspace1d(0.1,10,0.1,dstype=DataPHA)

           Blank PHA background data set

              dataspace1d(0.1,10,0.1, 1, 1, DataPHA)

        SEE ALSO
           dataspace2d
        """
        # support non-integrated grids with inclusive boundaries
        if dstype in (sherpa.data.Data1D, sherpa.astro.data.DataPHA):
            stop += step

        xlo,xhi,y = sherpa.utils.dataspace1d(start, stop, step=step,
                                             numbins=numbins)
        args = [xlo,xhi,y]
        kwargs={}
        
        if dstype is sherpa.astro.data.DataPHA:
            channel = numpy.arange(1,len(xlo)+1, dtype=float)
            args = [channel, y]
            #kwargs['bin_lo'] = xlo
            #kwargs['bin_hi'] = xhi
        elif dstype is not sherpa.data.Data1DInt:
            args = [xlo, y]

        if bkg_id is not None:
            self._get_pha_data(id).set_background(dstype('bkg_dataspace1d',
                                                         *args, **kwargs),
                                                  bkg_id)
        else:
            self.set_data(id, dstype('dataspace1d', *args, **kwargs))


    # DOC-NOTE: also in sherpa.utils
    def dataspace2d(self, dims, id=None, dstype=sherpa.astro.data.DataIMG):
        """
        dataspace2d

        SYNOPSIS
           Populates a blank 2D Sherpa image data set by data id

        SYNTAX

        Arguments:
           dims    -  array of image dimensions, i.e. [width,height]

           id      -  Sherpa data id
                      defaut is default data id

           dstype  -  Type of data set to use
                      default is DataIMG

        Returns:
           None

        DESCRIPTION
           Populates a blank 2D Sherpa image data set with logical coordinates
           by default and by Sherpa data id.

        SEE ALSO
           dataspace1d
        """
        x0, x1, y, shape = sherpa.utils.dataspace2d(dims)

        dataset=None
        if issubclass(dstype, (sherpa.astro.data.DataIMGInt,
                               sherpa.data.Data2DInt)):
            dataset = dstype('dataspace2d', x0-0.5, x1-0.5, x0+0.5, x1+0.5,
                             y, shape)
        else:
            dataset = dstype('dataspace2d', x0, x1, y, shape)

        self.set_data(id, dataset)


    # DOC-NOTE: also in sherpa.utils
    def unpack_arrays(self, *args):
        """
        unpack_arrays
        
        SYNOPSIS
           Read NumPy arrays into a dataset

        SYNTAX

        Arguments:
           array0     - first NumPy array | first CrateData obj

           ...

           arrayN     - last NumPy array | last CrateData obj

           dstype     - dataset type desired
                        default = Data1D

        Returns:
           Sherpa dataset

        DESCRIPTION
           Read NumPy arrays into a Sherpa dataset or read CrateData objects
           into a Sherpa dataset.  The list can include both NumPy arrays and
           CrateData objects together.

        SEE ALSO
           unpack_pha, unpack_arf, unpack_rmf, unpack_image, unpack_data
        """
        dataset = None
        try:
            dataset = sherpa.astro.io.read_arrays(*args)
        except AttributeError:
            # if the astro backend is not set, fall back on io module version.
            dataset = sherpa.io.read_arrays(*args)
        return dataset

    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: rework the Data type notes section.
    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_arrays(self, id, *args):
        """Create a data set from array values.

        Parameters
        ----------
        id : int or str
           The identifier for the data set to use.
        *args :
           Two or more arrays, followed by the type of data set to
           create.

        See Also
        --------
        copy_data : Copy a data set to a new identifier.
        delete_data : Delete a data set by identifier.
        get_data : Return the data set by identifier.
        load_data : Create a data set from a file.
        set_data : Set a data set.
        unpack_arrays :

        Notes
        -----
        The data type identifier, which defaults to `Data1D`,
        determines the number, and order, of the required inputs.

        `Data1D`
           required fields: x, y
           optional fields: statistical error, systematic error

        `Data1DInt`
           required fields: xlo, xhi, y
           optional fields: statistical error, systematic error

        `Data2D`
           required fields: x0, x1, y
           optional fields: shape, statistical error, systematic error
           The `shape` argument should be a tuple giving the
           size of the data (ny,nx).

        `Data2DInt`
           required fields: x0lo, x1lo, x0hi, x1hi, y
           optional fields: shape, statistical error, systematic error
           The `shape` argument should be a tuple giving the
           size of the data (ny,nx).

        `DataPHA`
           required fields: channel, counts
           optional fields: staterror, syserror, bin_lo, bin_hi,
             grouping, quality

        `DataIMG`
           The arrays should be 1D, not 2D.
           required fields: x0, x1, y
           optional fields: shape, statistical error, systematic error
           The `shape` argument should be a tuple giving the
           size of the data (ny,nx).

        Examples
        --------

        Create a 1D data set with three points:

        >>> load_arrays(1, [10, 12, 15], [4.2, 12.1, 8.4])

        Create a 1D data set, with the identifier 'prof', from the
        arrays `x` (independent axis), `y` (dependent axis), and `dy`
        (statistical error on the dependent axis):

        >>> load_arrays('prof', x, y, dy)

        Explicitly define the type of the data set:

        >>> load_arrays('prof', x, y, dy, Data1D)

        Data set 1 is a histogram, where the bins cover the range
        1-3, 3-5, and 5-7 with values 4, 5, and 9 respectively.

        >>> load_arrays(1, [1,3,5], [3,5,7], [4,5,9], Data1DInt)

        Create an image data set:

        >>> ivals = np.arange(12)
        >>> (y, x) = np.mgrid[0:3, 0:4]
        >>> x = x.flatten()
        >>> y = y.flatten()
        >>> load_arrays('img', x, y, ivals, (3,4), DataIMG)

        """
        self.set_data(id, self.unpack_arrays(*args))

    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: should unpack_ascii be merged into this?
    def unpack_table(self, filename, ncols=2, colkeys=None, dstype=Data1D):
        """Unpack a FITS binary file into a data structure.

        Parameters
        ----------
        filename :
           Identify the file to read: a file name, or a data structure
           representing the data to use, as used by the I/O backend in
           use by Sherpa: a `TABLECrate` for crates, as used by CIAO,
           or a list of AstroPy HDU objects.
        ncols : int, optional
           The number of columns to read in (the first `ncols` columns
           in the file). The meaning of the columns is determined by
           the `dstype` parameter.
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        dstype : optional
           The data class to use. The default is `Data1D`.

        Returns
        -------
        data :
           The class of the returned object is controlled by the
           `dstype` parameter.

        See Also
        --------
        load_table : Load a FITS binary file as a data set.
        set_data : Set a data set.
        unpack_ascii : Unpack an ASCII file into a data structure.

        Examples
        --------

        Read in the first two columns of the file, as the independent
        (X) and dependent (Y) columns of a data set:

        >>> d = unpack_table('sources.fits')

        Read in the first three columns (the third column is taken to
        be the error on the dependent variable):

        >>> d = unpack_table('sources.fits', ncols=3)

        Read in from columns 'RMID' and 'SUR_BRI':

        >>> d = unpack_table('rprof.fits', colkeys=['RMID', 'SUR_BRI'])

        The first three columns are taken to be the two independent
        axes of a two-dimensional data set (`x0` and `x1) and
        the dependent value (`y`):

        >>> d = unpack_table('fields.fits', ncols=3,
                             dstype=sherpa.astro.data.Data2D)

        When using the Crates I/O library, the file name can include
        CIAO Data Model syntax, such as column selection. This can
        also be done using the `colkeys` parameter, as shown above:

        >>> d = unpack_table('rprof.fits[cols rmid,sur_bri,sur_bri_err]',
                             ncols=3)

        """
        return sherpa.astro.io.read_table(filename, ncols, colkeys, dstype)

    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: the field listing really should be somewhere else
    ###           as it's needed in multiple places (ideally in the
    ###           DataX class documentation, but users may not find it)
    ### DOC-TODO: what do the shape arguments for Data2D/Data2DInt mean?
    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_table(self, id, filename=None, ncols=2, colkeys=None,
                   dstype=Data1D):
        """Load a FITS binary file as a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename :
           Identify the file to read: a file name, or a data structure
           representing the data to use, as used by the I/O backend in
           use by Sherpa: a `TABLECrate` for crates, as used by CIAO,
           or a list of AstroPy HDU objects.
        ncols : int, optional
           The number of columns to read in (the first `ncols` columns
           in the file). The meaning of the columns is determined by
           the `dstype` parameter.
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        dstype : optional
           The data class to use. The default is `Data1D`.

        See Also
        --------
        load_arrays : Create a data set from array values.
        load_ascii : Load an ASCII file as a data set.
        load_image : Load an image as a data set.
        set_data : Set a data set.
        unpack_table : Unpack a FITS binary table into a data structure.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The column order for the different data types are as follows,
        where `x` indicates an independent axis and `y` the dependent
        axis:

        `Data1D`
           required fields: x, y
           optional fields: statistical error, systematic error

        `Data1DInt`
           required fields: xlo, xhi, y
           optional fields: statistical error, systematic error

        `Data2D`
           required fields: x0, x1, y
           optional fields: shape, statistical error, systematic error

        `Data2DInt`
           required fields: x0lo, x1lo, x0hi, x1hi, y
           optional fields: shape, statistical error, systematic error

        Examples
        --------

        Read in the first two columns of the file, as the independent
        (X) and dependent (Y) columns of the default data set:

        >>> load_table('sources.fits')

        Read in the first three columns (the third column is taken to
        be the error on the dependent variable):

        >>> load_table('sources.fits', ncols=3)

        Read in from columns 'RMID' and 'SUR_BRI' into data set
        'prof':

        >>> load_table('prof', 'rprof.fits',
                       colkeys=['RMID', 'SUR_BRI'])

        The first three columns are taken to be the two independent
        axes of a two-dimensional data set (`x0` and `x1) and
        the dependent value (`y`):

        >>> load_table('fields.fits', ncols=3,
                       dstype=sherpa.astro.data.Data2D)

        When using the Crates I/O library, the file name can include
        CIAO Data Model syntax, such as column selection. This can
        also be done using the `colkeys` parameter, as shown above:

        >>> load_table('prof',
                       'rprof.fits[cols rmid,sur_bri,sur_bri_err]',
                       ncols=3)

        Read in a data set using Crates:

        >>> cr = pycrates.read_file('table.fits')
        >>> load_table(cr)

        Read in a data set using Crates:

        >>> hdus = astropy.io.fits.open('table.fits')
        >>> load_table(hdus)

        """
        if filename is None:
            id, filename = filename, id
        
        self.set_data(id, self.unpack_table(filename, ncols, colkeys, dstype))
        
    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: should unpack_ascii be merged into unpack_table?
    ### DOC-TODO: I am going to ignore the crates support here as
    ###           it is somewhat meaningless, since the crate could
    ###           have been read from a FITS binary table.
    def unpack_ascii(self, filename, ncols=2, colkeys=None,
                     dstype=Data1D, sep=' ', comment='#'):
        """Unpack an ASCII file into a data structure.
        
        Parameters
        ----------
        filename : str
           The name of the file to read in. Selection of the relevant
           column depends on the I/O library in use (Crates or
           AstroPy).
        ncols : int, optional
           The number of columns to read in (the first `ncols` columns
           in the file). The meaning of the columns is determined by
           the `dstype` parameter.
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        dstype : optional
           The data class to use. The default is `Data1D`.

        Returns
        -------
        data :
           The class of the returned object is controlled by the
           `dstype` parameter.

        See Also
        --------
        load_ascii : Load an ASCII file as a data set.
        set_data : Set a data set.
        unpack_table : Unpack a FITS binary file into a data structure.

        Examples
        --------

        Read in the first two columns of the file, as the independent
        (X) and dependent (Y) columns of a data set:

        >>> d = unpack_ascii('sources.dat')

        Read in the first three columns (the third column is taken to
        be the error on the dependent variable):

        >>> d = unpack_ascii('sources.dat', ncols=3)

        Read in from columns 'col2' and 'col3':

        >>> d = unpack_ascii('tbl.dat', colkeys=['col2', 'col3'])

        The first three columns are taken to be the two independent
        axes of a two-dimensional data set (`x0` and `x1) and
        the dependent value (`y`):

        >>> d = unpack_ascii('fields.dat', ncols=3,
                             dstype=sherpa.astro.data.Data2D)

        When using the Crates I/O library, the file name can include
        CIAO Data Model syntax, such as column selection. This can
        also be done using the `colkeys` parameter, as shown above:

        >>> d = unpack_ascii('tbl.dat[cols rmid,sur_bri,sur_bri_err]',
                             ncols=3)

        """
        return sherpa.astro.io.read_ascii(filename, ncols, colkeys, dstype,
                                          sep=sep, comment=comment)
    
    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: I am going to ignore the crates support here as
    ###           it is somewhat meaningless, since the crate could
    ###           have been read from a FITS binary table.
    ### DOC-TODO: how best to include datastack support?
    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_ascii(self, id, filename=None, ncols=2, colkeys=None,
                   dstype=Data1D, sep=' ', comment='#'):
        """Load an ASCII file as a data set.

        The standard behavior is to create a single data set, but
        multiple data sets can be loaded with this command, as
        described in the `sherpa.astro.datastack` module.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to read in. Selection of the relevant
           column depends on the I/O library in use (Crates or
           AstroPy).
        ncols : int, optional
           The number of columns to read in (the first `ncols` columns
           in the file). The meaning of the columns is determined by
           the `dstype` parameter.
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.
        dstype : optional
           The data class to use. The default is `Data1D`.

        See Also
        --------
        load_table : Load a FITS binary file as a data set.
        load_image : Load an image as a data set.
        set_data : Set a data set.
        unpack_ascii : Unpack an ASCII file into a data structure.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The column order for the different data types are as follows,
        where `x` indicates an independent axis and `y` the dependent
        axis:

        `Data1D`
           required fields: x, y
           optional fields: statistical error, systematic error

        `Data1DInt`
           required fields: xlo, xhi, y
           optional fields: statistical error, systematic error

        `Data2D`
           required fields: x0, x1, y
           optional fields: shape, statistical error, systematic error

        `Data2DInt`
           required fields: x0lo, x1lo, x0hi, x1hi, y
           optional fields: shape, statistical error, systematic error

        Examples
        --------

        Read in the first two columns of the file, as the independent
        (X) and dependent (Y) columns of the default data set:

        >>> load_ascii('sources.dat')

        Read in the first three columns (the third column is taken to
        be the error on the dependent variable):

        >>> load_ascii('sources.dat', ncols=3)

        Read in from columns 'RMID' and 'SUR_BRI' into data set
        'prof':

        >>> load_ascii('prof', 'rprof.dat',
                       colkeys=['RMID', 'SUR_BRI'])

        The first three columns are taken to be the two independent
        axes of a two-dimensional data set (`x0` and `x1) and
        the dependent value (`y`):

        >>> load_ascii('fields.txt', ncols=3,
                       dstype=sherpa.astro.data.Data2D)

        When using the Crates I/O library, the file name can include
        CIAO Data Model syntax, such as column selection. This can
        also be done using the `colkeys` parameter, as shown above:

        >>> load_ascii('prof',
                       'rprof.dat[cols rmid,sur_bri,sur_bri_err]',
                       ncols=3)

        """
        if filename is None:
            id, filename = filename, id
            
        self.set_data(id, self.unpack_ascii(filename, ncols=ncols,
                                            colkeys=colkeys, dstype=dstype,
                                            sep=sep, comment=comment ))
        
    # DOC-NOTE: also in sherpa.utils
    def unpack_data(self, filename, *args, **kwargs):
        """
        unpack_data

        SYNOPSIS
           Read spectrum, table, or ASCII data into a dataset

        SYNTAX

        Arguments:
           filename   - filename and path

        Returns:
           Sherpa dataset

        DESCRIPTION
           Read PHA spectrum data, FITS table data , or tabular data from a
           column-based text file into a Sherpa dataset given a filename
           and path.

        SEE ALSO
           unpack_pha, unpack_arf, unpack_rmf, unpack_image, unpack_data,
           unpack_table, unpack_ascii
        """
        try:
            data = self.unpack_pha(filename, *args, **kwargs)
        except:
            try:
                data = self.unpack_image(filename, *args, **kwargs)
            except:
                try:
                    data = self.unpack_table(filename, *args, **kwargs)
                except:
                    try:
                        data = self.unpack_ascii(filename, *args, **kwargs)
                    except:
                        raise

        return data

    # DOC-NOTE: also in sherpa.utils
    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_data(self, id, filename=None, *args, **kwargs):
        """
        load_data

        SYNOPSIS
           Load spectrum, table, or ASCII data by id

        SYNTAX

        Arguments:
           id         - data id
                        default = default data id

           filename   - filename and path

      Returns:
           None

        DESCRIPTION
           Load PHA spectrum data, FITS table data, or tabular data from a
           column-based text file into a Sherpa dataset given a filename
           and path by data id.
        
        SEE ALSO
           load_pha, load_arf, load_rmf, load_data, load_image,
           load_bkg, load_table, load_ascii
        """
        if filename is None:
            id, filename = filename, id

        data = self.unpack_data(filename, *args, **kwargs)

        if type(data) is list:
            num = len(data)
            if num > 1:
                for id, pha in enumerate(data):
                    self.set_data(id+1, pha )
                info("Multiple data sets have been input: 1-%s" % num)
            else:
                self.set_data(id, data.pop() )
        else:
            self.set_data(id, data)

    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    def unpack_image(self, arg, coord='logical',
                     dstype=sherpa.astro.data.DataIMG):
        """Create an image data structure.

        Parameters
        ----------
        arg :
           Identify the image file: a file name, or a data structure
           representing the data to use, as used by the I/O backend in
           use by Sherpa: an `IMAGECrate` for crates, as used by CIAO,
           or a list of AstroPy HDU objects.
        coord : { 'logical', 'image', 'physical', 'world', 'wcs' }, optional
           Ensure that the image contains the given coordinate system.
        dstype : optional
           The image class to use. The default is `DataIMG`.

        Returns
        -------
        img :
           The class of the returned object is controlled by the
           `dstype` parameter.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the image does not contain the requested coordinate
           system.

        See Also
        --------
        load_image : Load a file as an image data set.
        set_data : Set a data set.

        Examples
        --------

        >>> img1 = unpack_img("img.fits")
        >>> set_data(img1)

        >>> img = unpack_img('img.fits', 'physical')

        Read in an image using Crates:

        >>> cr = pycrates.read_file('broad.img')
        >>> idata = unpack_img(cr)

        Read in an image using AstroPy:

        >>> hdus = astropy.io.fits.open('broad.img')
        >>> idata = unpack_img(hdus)

        """
        return sherpa.astro.io.read_image(arg, coord, dstype)

    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_image(self, id, arg=None, coord='logical',
                     dstype=sherpa.astro.data.DataIMG):
        """Load a file as an image data set.

        load_image

        SYNOPSIS
           Load image data by id

        SYNTAX

        Arguments:
           id         - dataset ID
                        default = default data id
           arg        - filename and path | IMAGECrate obj | PyFITS HDUList obj

           coord      - string keyword identifying coordinate system
                        choices include: logical, image
                                         physical
                                         world, wcs
                        default = logical

           dstype     - Sherpa dataset type (DataIMG, DataIMGInt)
                        default = DataIMG

        Returns:
           None

        DESCRIPTION
           Load image data from a FITS file into a Sherpa dataset given a
           filename by data id or load in image data from a Crate into a Sherpa
           dataset given a IMAGECrate object by data id or read in image data
           from a HDUList into a Sherpa dataset by data id.

        SEE ALSO
           load_pha, load_arf, load_rmf, load_data, load_table,
           load_bkg
        """
        if arg is None:
            id, arg = arg, id
        self.set_data(id, self.unpack_image(arg, coord, dstype))

    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    ### DOC-TODO: what does this return when given a PHA2 file?
    def unpack_pha(self, arg, use_errors=False):
        """Create a PHA data structure.

        Any instrument or background data sets referenced in the
        header of the PHA file - e.g. with the ANCRFILE, RESPFILE,
        and BACKFILE keywords - will also be loaded.

        Parameters
        ----------
        arg :
           Identify the PHA file: a file name, or a data structure
           representing the data to use, as used by the I/O backend in
           use by Sherpa: a `TABLECrate` for crates, as used by CIAO,
           or a list of AstroPy HDU objects.
        use_errors : bool, optional
           If `True` then the statistical errors are taken from the
           input data, rather than calculated by Sherpa from the
           count values. The default is `False`.

        Returns
        -------
        pha : sherpa.astro.data.DataPHA instance

        See Also
        --------
        load_pha : Load a file as a PHA data set.
        pack_pha : Convert a PHA data set into a file structure.
        set_data : Set a data set.

        Examples
        --------

        >>> pha1 = unpack_arf("src1.pi")
        >>> pha2 = unpack_arf("field.pi")
        >>> set_data(1, pha1)
        >>> set_bkg(1, pha2)

        Read in a PHA file using Crates:

        >>> cr = pycrates.read_file("src.fits")
        >>> pha = unpack_pha(cr)

        Read in a PHA file using AstroPy:

        >>> hdus = astropy.io.fits.open("src.fits")
        >>> pha = unpack_pha(hdus)

        """
        use_errors = sherpa.utils.bool_cast(use_errors)
        return sherpa.astro.io.read_pha(arg, use_errors)


    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    ### DOC-TODO: what does this return when given a PHA2 file?
    def unpack_bkg(self, arg, use_errors=False):
        """Create a PHA data structure for a background data set.

        Any instrument information referenced in the header of the PHA
        file - e.g. with the ANCRFILE and RESPFILE, keywords - will
        also be loaded. Unlike `unpack_pha`, background files will not
        be loaded.

        Parameters
        ----------
        arg :
           Identify the PHA file: a file name, or a data structure
           representing the data to use, as used by the I/O backend in
           use by Sherpa: a `TABLECrate` for crates, as used by CIAO,
           or a list of AstroPy HDU objects.
        use_errors : bool, optional
           If `True` then the statistical errors are taken from the
           input data, rather than calculated by Sherpa from the
           count values. The default is `False`.

        Returns
        -------
        pha : sherpa.astro.data.DataPHA instance

        See Also
        --------
        load_bkg : Load the background from a file and add it to a PHA data set.
        set_bkg : Set the background for a data set.

        Examples
        --------

        >>> pha1 = unpack_arf("src1.pi")
        >>> pha2 = unpack_bkg("field.pi")
        >>> set_data(1, pha1)
        >>> set_bkg(1, pha2)

        Read in a PHA file using Crates:

        >>> cr = pycrates.read_file("bg.fits")
        >>> pha = unpack_pha(cr)

        Read in a PHA file using AstroPy:

        >>> hdus = astropy.io.fits.open("bg.fits")
        >>> pha = unpack_pha(hdus)

        """
        use_errors = sherpa.utils.bool_cast(use_errors)
        return sherpa.astro.io.read_pha(arg, use_errors, True)


    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_pha(self, id, arg=None, use_errors=False):
        """Load a PHA data set.

        load_pha

        SYNOPSIS
           Load PHA data by id

        SYNTAX

        Arguments:
           id         - dataset ID
                        default = default data id

           arg        - filename and path | PHACrate obj | PyFITS HDUList obj

           use_errors - flag to use errors
                        default = False

        Returns:
           None

        DESCRIPTION
           Load PHA data from a FITS file or a PHACrate object or a PyFITS
           HDUList object into a Sherpa dataset by data id.

        SEE ALSO
           load_image, load_arf, load_rmf, load_data, load_table,
           load_bkg
        """
        if arg is None:
            id, arg = arg, id

        phasets = self.unpack_pha(arg, use_errors)

        if numpy.iterable(phasets):
            num = len(phasets)
            for id, pha in enumerate(phasets):
                self.set_data(id+1, pha )
            if num > 1:
                info("Multiple data sets have been input: 1-%s" % num)
        else:
            self.set_data(id, phasets)

    def _get_pha_data(self, id):
        data = self.get_data(id)
        if not isinstance(data, sherpa.astro.data.DataPHA):
            raise ArgumentErr('nopha', self._fix_id(id))
        return data

    def _get_img_data(self, id):
        data = self.get_data(id)
        if not isinstance(data, sherpa.astro.data.DataIMG):
            raise ArgumentErr('noimg', self._fix_id(id))
        return data


    # def _read_error(self, filename, *args, **kwargs):
    #     err = None
    #     try:
    #         err = sherpa.astro.io.backend.get_ascii_data(filename, *args,
    #                                                    **kwargs)[1].pop()
    #     except:
    #         try:
    #             err = sherpa.astro.io.backend.get_table_data(filename, *args,
    #                                                        **kwargs)[1].pop()
    #         except:
    #             try:
    #                 err = sherpa.astro.io.read_image(filename, *args, **kwargs)
    #                 err = err.get_dep()
    #             except:
    #                 raise

    #     return err


    def load_filter(self, id, filename=None, bkg_id=None, ignore=False, 
                    ncols=2, *args, **kwargs):
        """
        load_filter

        SYNOPSIS
           Load the dataset filter from file

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filename   - filename with path

           bkg_id     - background data id
                        default = default background data id

           ignore     - non-zero values ignore instead of notice
                        default = False

           ncols      - number of columns to read from
                        default = 2

           colkeys    - column keys
                        default = None

           sep        - separator character
                        default = ' '

           comment    - comment character
                        default = '#'

        Returns:
           None

        DESCRIPTION
           Load the filter for a dataset from file by data id.

        EXAMPLE
           load_filter("data.dat", colkeys=["FILTER"])

        SEE ALSO
            set_filter
        """
        if filename is None:
            id, filename = filename, id

        self.set_filter(id, self._read_user_model(filename, *args, **kwargs)[1],
                        bkg_id=bkg_id, ignore=ignore)


    ### Ahelp ingest: 2015-04-30 DJB
    ### DOC-TODO: does ncols make sense here? (have removed for now)
    ### DOC-TODO: labelling as AstroPy; i.e. assuming conversion
    ###           from PyFITS lands soon.
    ### DOC-TODO: prob. needs a review as the existing ahelp documentation
    ###           talks about 2 cols, but experimentation suggests 1 col.
    def load_grouping(self, id, filename=None, bkg_id=None, *args, **kwargs):
        """Load the grouping scheme from a file and add to a PHA data set.

        This function sets the grouping column but does not
        automatically group the data, since the quality array may also
        need updating. The `group` function will apply the grouping
        information.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file that contains the grouping
           information. This file can be a FITS table or an ASCII
           file. Selection of the relevant column depends on the I/O
           library in use (Crates or AstroPy).
        bkg_id : int or str, optional
           Set if the grouping scheme should be associated with the
           background associated with the data set.
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.

        See Also
        --------
        get_grouping : Return the gouping array for a PHA data set.
        group : Turn on the grouping for a PHA data set.
        load_quality : Load the quality array from a file and add to a PHA data set.
        save_grouping : Save the grouping scheme to a file.
        set_grouping : Apply a set of grouping flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        There is no check made to see if the grouping array contains
        valid data.

        Examples
        --------

        When using Crates as the I/O library, select the grouping
        column from the file `src.pi`, and use it to set the
        values in the default data set:

        >>> load_grouping('src.pi[cols grouping]')

        Use the `colkeys` option to define the column in the input
        file:

        >>> load_grouping('src.pi', colkeys=['grouping'])

        Load the first column in `grp.dat` and use it to populate
        the grouping array of the data set called 'core'.

        >>> load_grouping('core', 'grp.dat')

        Use `group_counts` to calculate a grouping scheme for the
        data set labelled 'src1', save this scheme to the file
        'grp.dat', and then load this scheme in for data set
        'src2'.

        >>> group_counts('src1', 10)
        >>> save_grouping('src1', 'grp.dat')
        >>> load_grouping('src2', 'grp.dat', colkeys=['groups'])

        """
        if filename is None:
            id, filename = filename, id

        self.set_grouping(id,
            self._read_user_model(filename, *args, **kwargs)[1], bkg_id=bkg_id)

    ### Ahelp ingest: 2015-04-30 DJB
    ### DOC-TODO: labelling as AstroPy; i.e. assuming conversion
    ###           from PyFITS lands soon.
    def load_quality(self, id, filename=None, bkg_id=None, *args, **kwargs):
        """Load the quality array from a file and add to a PHA data set.

        This function sets the quality column but does not
        automatically ignore any columns marked as "bad". Use the
        `ignore_bad` function to apply the new quality information.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file that contains the quality
           information. This file can be a FITS table or an ASCII
           file. Selection of the relevant column depends on the I/O
           library in use (Crates or AstroPy).
        bkg_id : int or str, optional
           Set if the quality array should be associated with the
           background associated with the data set.
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.

        See Also
        --------
        get_quality : Return the quality array for a PHA data set.
        ignore_bad : Exclude channels marked as bad in a PHA data set.
        load_grouping : Load the grouping scheme from a file and add to a PHA data set.
        save_quality: Save the quality array to a file.
        set_quality : Apply a set of quality flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        There is no check made to see if the quality array contains
        valid data.

        Examples
        --------

        When using Crates as the I/O library, select the quality
        column from the file `src.pi`, and use it to set the
        values in the default data set:

        >>> load_quality('src.pi[cols quality]')

        Use the `colkeys` option to define the column in the input
        file:

        >>> load_quality('src.pi', colkeys=['quality'])

        Load the first column in `grp.dat` and use it to populate
        the quality array of the data set called 'core'.

        >>> load_quality('core', 'grp.dat')

        """
        if filename is None:
            id, filename = filename, id

        self.set_quality(id,
            self._read_user_model(filename, *args, **kwargs)[1], bkg_id=bkg_id)

    def set_filter(self, id, val=None, bkg_id=None, ignore=False):
        """
        set_filter

        SYNOPSIS
           Set the dataset filter by data id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           ignore     - non-zero values ignore instead of notice
                        default = False

           val        - array of 0s or 1s

        Returns:
           None

        DESCRIPTION
           Set the filter of a dataset by data id.  

        EXAMPLE
           set_filter([0, 1, 1, ...])

        SEE ALSO
           load_filter
        """
        if val is None:
            val, id = id, val

        filter = numpy.asarray(val, dtype=numpy.bool_)

        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        if numpy.iterable(d.mask):
            if len(d.mask) == len(filter):
                if not ignore:
                    d.mask |= filter
                else:
                    d.mask &= ~filter
            else:
                raise sherpa.utils.err.DataErr('mismatch',
                                               len(d.mask), len(filter))
        else:
            if len(d.get_y(False)) == len(filter):
                if not ignore:
                    d.mask = filter
                else:
                    d.mask = ~filter
            else:
                raise sherpa.utils.err.DataErr('mismatch',
                                               len(d.get_y(False)), len(filter))


    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-06 DJB
    ### DOC-NOTE: is ncols really 2 here? Does it make sense?
    def load_staterror(self, id, filename=None, bkg_id=None, *args, **kwargs):
        """Load the statistical errors from a file.

        Read in a column or image from a file and use the values
        as the statistical errors for a data set. This over rides
        the errors calculated by any statistic, such as
        `chi2gehrels` or `chi2datavar`.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to read in. Supported formats depends
           on the I/O library in use (Crates or AstroPy) and the
           type of data set (e.g. 1D or 2D).
        bkg_id : int or str, optional
           Set to identify which background component to set. The
           default value (`None`) means that this is for the source
           component of the data set.
        ncols : int, optional
           The number of columns to read in (the first `ncols` columns
           in the file).
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.

        See Also
        --------
        get_staterror : Return the statistical error on the dependent axis of a data set.
        load_syserror : Load the systematic errors from a file.
        set_staterror : Set the statistical errors on the dependent axis of a data set.
        set_stat : Set the statistical method.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Read in the first column from 'tbl.dat':

        >>> load_staterror('tbl.dat')

        Use the column labelled 'col3'

        >>> load_staterror('tbl.dat', colkeys=['col3'])

        When using the Crates I/O library, the file name can include
        CIAO Data Model syntax, such as column selection:

        >>> load_staterror('tbl.dat[cols col3]')

        Read in the first column from the file 'errors.fits' as the
        statistical errors for the 'core' data set:

        >>> load_staterror('core', 'errors.fits')

        The data set labelled 'img' is loaded from the file
        'image.fits' and the statistical errors from 'err.fits'.
        The dimensions of the two images must be the same.

        >>> load_image('img', 'image.fits')
        >>> load_staterror('img', 'err.fits')

        """
        if filename is None:
            id, filename = filename, id

        self.set_staterror(id,
            self._read_user_model(filename, *args, **kwargs)[1], bkg_id=bkg_id)

    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-06 DJB
    ### DOC-NOTE: is ncols really 2 here? Does it make sense?
    def load_syserror(self, id, filename=None, bkg_id=None, *args, **kwargs):
        """Load the systematic errors from a file.

        Read in a column or image from a file and use the values
        as the systematic errors for a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to read in. Supported formats depends
           on the I/O library in use (Crates or AstroPy) and the
           type of data set (e.g. 1D or 2D).
        bkg_id : int or str, optional
           Set to identify which background component to set. The
           default value (`None`) means that this is for the source
           component of the data set.
        ncols : int, optional
           The number of columns to read in (the first `ncols` columns
           in the file).
        colkeys : array of str, optional
           An array of the column name to read in. The default is
           `None`.
        sep : str, optional
           The separator character. The default is ' '.
        comment : str, optional
           The comment character. The default is '#'.

        See Also
        --------
        get_syserror : Return the systematic error on the dependent axis of a data set.
        load_staterror : Load the statistical errors from a file.
        set_syserror : Set the systematic errors on the dependent axis of a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Read in the first column from 'tbl.dat':

        >>> load_syserror('tbl.dat')

        Use the column labelled 'col3'

        >>> load_syserror('tbl.dat', colkeys=['col3'])

        When using the Crates I/O library, the file name can include
        CIAO Data Model syntax, such as column selection:

        >>> load_syserror('tbl.dat[cols col3]')

        Read in the first column from the file 'errors.fits' as the
        systematic errors for the 'core' data set:

        >>> load_syserror('core', 'errors.fits')

        The data set labelled 'img' is loaded from the file
        'image.fits' and the systematic errors from 'syserr.fits'.
        The dimensions of the two images must be the same.

        >>> load_image('img', 'image.fits')
        >>> load_syserror('img', 'syserr.fits')

        """
        if filename is None:
            id, filename = filename, id

        self.set_syserror(id,
            self._read_user_model(filename, *args, **kwargs)[1], bkg_id=bkg_id)

    def set_dep(self, id, val=None, bkg_id=None):
        """
        set_dep

        SYNOPSIS
           Set the dependent variable of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           val        - dependent variable array or scalar

           bkg_id     - background id
                        default = None

        Returns:
           None

        DESCRIPTION
           Set the dependent variable of a data set by data id or bkg_id.

        EXAMPLE
           set_dep([1,2,3,...])

           set_dep(1,1)

           set_dep([1,2,3,...], bkg_id=1)

        SEE ALSO
           get_dep, get_indep, get_axes
        """
        if val is None:
            val, id = id, val
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)

        dep=None
        if isinstance(d, sherpa.astro.data.DataPHA):
            if numpy.iterable(val):
                dep = numpy.asarray(val, SherpaFloat)
            else:
                val = SherpaFloat(val)
                dep = numpy.array([val]*len(d.channel))
            d.counts = dep
        else:
            if numpy.iterable(val):
                dep = numpy.asarray(val, SherpaFloat)
            else:
                val = SherpaFloat(val)
                dep = numpy.array([val]*len(d.get_indep()[0]))
            d.y = dep


    set_counts = set_dep


    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def set_staterror(self, id, val=None, fractional=False, bkg_id=None):
        """Set the statistical errors on the dependent axis of a data set.

        These values over-ride the errors calculated by any statistic,
        such as `chi2gehrels` or `chi2datavar`.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        val : array or scalar
           The systematic error.
        fractional : bool, optional
           If `False` (the default value), then the `val` parameter is
           the absolute value, otherwise the `val` parameter
           represents the fractional error, so the absolute value is
           calculated as `get_dep() * val` (and `val` must be
           a scalar).
        bkg_id : int or str, optional
           Set to identify which background component to set. The
           default value (`None`) means that this is for the source
           component of the data set.

        See Also
        --------
        load_staterror : Load the statistical errors from a file.
        load_syserror : Load the systematic errors from a file.
        set_syserror : Set the systematic errors on the dependent axis of a data set.
        get_error : Return the errors on the dependent axis of a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `val` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `val` parameters,
        respectively.

        Examples
        --------

        Set the statistical error for the default data set to the value
        in `dys` (a scalar or an array):

        >>> set_staterror(dys)

        Set the statistical error on the `core` data set to be 5% of
        the data values:

        >>> set_staterror('core', 0.05, fractional=True)

        """
        if val is None:
            val, id = id, val
        err=None

        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)

        fractional=sherpa.utils.bool_cast(fractional)
        if numpy.iterable(val):
            err = numpy.asarray(val, SherpaFloat)
        elif val is not None:
            val = SherpaFloat(val)
            if fractional:
                err = val*d.get_dep()
            else:
                err = numpy.array([val]*len(d.get_dep()))
        d.staterror = err


    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def set_syserror(self, id, val=None, fractional=False, bkg_id=None):
        """Set the systematic errors on the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        val : array or scalar
           The systematic error.
        fractional : bool, optional
           If `False` (the default value), then the `val` parameter is
           the absolute value, otherwise the `val` parameter
           represents the fractional error, so the absolute value is
           calculated as `get_dep() * val` (and `val` must be
           a scalar).
        bkg_id : int or str, optional
           Set to identify which background component to set. The
           default value (`None`) means that this is for the source
           component of the data set.

        See Also
        --------
        load_staterror : Set the statistical errors on the dependent axis of a data set.
        load_syserror : Set the systematic errors on the dependent axis of a data set.
        set_staterror : Set the statistical errors on the dependent axis of a data set.
        get_error : Return the errors on the dependent axis of a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `val` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `val` parameters,
        respectively.

        Examples
        --------

        Set the systematic error for the default data set to the value
        in `dys` (a scalar or an array):

        >>> set_syserror(dys)

        Set the systematic error on the `core` data set to be 5% of
        the data values:

        >>> set_syserror('core', 0.05, fractional=True)

        """
        if val is None:
            val, id = id, val
        err=None

        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        fractional=sherpa.utils.bool_cast(fractional)
        
        if numpy.iterable(val):
            err = numpy.asarray(val, SherpaFloat)
        elif val is not None:
            val = SherpaFloat(val)
            if fractional:
                err = val*d.get_dep()
            else:
                err = numpy.array([val]*len(d.get_dep()))
        d.syserror = err


    ### Ahelp ingest: 2015-05-02 DJB
    def set_exposure(self, id, exptime=None, bkg_id=None):
        """Change the exposure time of a PHA data set.

        The exposure time of a PHA data set is taken from the
        EXPTIME keyword in its header, but it can be changed
        once the file has been loaded.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        exptime : num
           The exposure time, in seconds.
        bkg_id : int or str, optional
           Set to identify which background component to set.  The
           default value (`None`) means that this is for the source
           component of the data set.

        See Also
        --------
        get_exposure : Return the exposure time of a PHA data set.
        set_areascal : Change the area scaling of a PHA data set.
        set_backscal : Change the area scaling of a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `exptime` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `exptime` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Increase the exposure time of the default data set by 5 per
        cent.

        >>> etime = get_exposure()
        >>> set_exposure(etime * 1.05)

        Use the EXPOSURE value from the ARF, rather than the EXPTIME
        value from the PHA file, for data set 2:

        >>> set_exposure(2, get_arf(2).exposure)

        Set the exposure time of the second background component
        of the 'jet' data set.

        >>> set_exposure('jet', 12324.45, bkg_id=2)

        """
        if exptime is None:
            exptime, id = id, exptime

        if exptime is not None:
            exptime = SherpaFloat(exptime)

        if bkg_id is not None:
            self.get_bkg(id,bkg_id).exposure = exptime
        else:
            self._get_pha_data(id).exposure = exptime


    def set_backscal(self, id, backscale=None, bkg_id=None):
        """Change the area scaling of a PHA data set.

        set_backscal

        SYNOPSIS
           Set the source or background extraction region areas by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           backscale  - backscal value

           bkg_id     - background data id
                        default = None

        Returns:
           None

        DESCRIPTION
           Set the extraction region areas of a source PHA dataset by data id
           or of a background dataset by bkg_id. Backscale can be defined as the
           ratio of the area of the source (or background) extraction region in
           image pixels to the total number of image pixels. The fact that
           there is no ironclad definition for this quantity does not matter so
           long as the backscale for a source dataset and its associated
           background dataset are defined in the similar manner, because only
           the ratio of source and background backscale is used in analyses.

        EXAMPLE
           set_backscal(2e-6)

           set_backscal(2, 2e-6)

           set_backscal(1, 1e-5, 1)

           set_backscal([1e-6, 1e-6, 1e-6, ...])

        SEE ALSO
           set_exposure, set_areascal
        """
        if backscale is None:
            backscale, id = id, backscale

        if numpy.iterable(backscale):
            backscale = numpy.asarray(backscale)
        elif backscale is not None:
            backscale = SherpaFloat(backscale)

        if bkg_id is not None:
            self.get_bkg(id,bkg_id).backscal = backscale
        else:
            self._get_pha_data(id).backscal = backscale


    def set_areascal(self, id, area=None, bkg_id=None):
        """Change the area scaling of a PHA data set.

        set_areascal

        SYNOPSIS
           Set the source or background fractional area by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           area       - areascal value [pixel]

           bkg_id     - background data id
                        default = None

        Returns:
           None

        DESCRIPTION
           Set the fractional area of a source PHA dataset by data id or of a 
           background data by bkg_id.

        EXAMPLE
           set_areascal(0.75)

           set_areascal(2, 0.75)

           set_areascal(1, 0.75, 1)

        SEE ALSO
           set_backscal, set_exposure
        """
        if area is None:
            area, id = id, area

        if area is not None:
            area = SherpaFloat(area)

        if bkg_id is not None:
            self.get_bkg(id,bkg_id).areascal = area
        else:
            self._get_pha_data(id).areascal = area

    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def get_staterror(self, id=None, filter=False, bkg_id=None):
        """Return the statistical error on the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.
        bkg_id : int or str, optional
           Set if the values returned should be from the given
           background component, instead of the source data set.

        Returns
        -------
        axis : array
           The statistical error for each data point. This may be
           estimated from the data (e.g. with the `chi2gehrels`
           statistic) or have been set explicitly (`set_staterror`).
           For PHA data sets, the return array will match the grouping
           scheme applied to the data set.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        get_syserror : Return the systematic errors on the dependent axis of a data set.
        list_data_ids : List the identifiers for the loaded data sets.
        set_staterror : Set the statistical errors on the dependent axis of a data set.

        Examples
        --------

        If not explicitly given, the statistical errors on a data set
        may be calculated from the data values (the independent axis),
        depending on the chosen statistic:

        >>> load_arrays(1, [10,15,19], [4,5,9])
        >>> set_stat('chi2datavar')
        >>> get_staterror()
        array([ 2.        ,  2.23606798,  3.        ])
        >>> set_stat('chi2gehrels')
        >>> get_staterror()
        array([ 3.17944947,  3.39791576,  4.122499  ])

        If the statistical errors are set - either when the data set
        is created or with a call to `set_errors` - then these values
        will be used, no matter the statistic:

        >>> load_arrays(1, [10,15,19], [4,5,9], [2,3,5])
        >>> set_stat('chi2datavar')
        >>> get_staterror()
        array([2, 3, 5])
        >>> set_stat('chi2gehrels')
        >>> get_staterror()
        array([2, 3, 5])

        """
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        return d.get_staterror(filter, self.get_stat().calc_staterror)


    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def get_syserror(self, id=None, filter=False, bkg_id=None):
        """Return the systematic error on the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.
        bkg_id : int or str, optional
           Set if the values returned should be from the given
           background component, instead of the source data set.

        Returns
        -------
        axis : array
           The systematic error for each data point.

        Raises
        ------
        sherpa.utils.err.DataErr
           If the data set has no systematic errors.
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        get_staterror : Return the statistical errors on the dependent axis of a data set.
        list_data_ids : List the identifiers for the loaded data sets.
        set_syserror : Set the systematic errors on the dependent axis of a data set.

        """
        d = self.get_data(id)
        id = self._fix_id(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        err = d.get_syserror(filter)
        if err is None or not numpy.iterable(err):
            raise sherpa.utils.err.DataErr('nosyserr', id)
        return err


    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def get_error(self, id=None, filter=False, bkg_id=None):
        """Return the errors on the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.
        bkg_id : int or str, optional
           Set if the values returned should be from the given
           background component, instead of the source data set.

        Returns
        -------
        axis : array
           The error for each data point, formed by adding the
           statistical and systematic errors in quadrature.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_dep : Return the dependent axis of a data set.
        get_staterror : Return the statistical errors on the dependent axis of a data set.
        get_syserror : Return the systematic errors on the dependent axis of a data set.
        list_data_ids : List the identifiers for the loaded data sets.

        """
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        return d.get_error(filter, self.get_stat().calc_staterror)


    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def get_indep(self, id=None, filter=False, bkg_id=None):
        """Return the independent axes of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.
        bkg_id : int or str, optional
           Set if the values returned should be from the given
           background component, instead of the source data set.

        Returns
        -------
        axis : tuple of arrays
           The independent axis values. These are the values at which
           the model is evaluated during fitting. The values returned
           depend on the coordinate system in use for the data set (as
           set by `set_coord`). For PHA data sets the value returned
           is always in channels, whatever the `set_analysis` setting
           is, and does not follow any grouping setting for the data
           set.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_axes : Return information about the independent axes of a data set.
        get_dep : Return the dependent axis of a data set.
        list_data_ids : List the identifiers for the loaded data sets.
        set_coord : Set the coordinate system to use for image analysis.

        Examples
        --------

        For a one-dimensional data set, the X values are returned:

        >>> load_arrays(1, [10,15,19], [4,5,9], Data1D)
        >>> get_indep()
        (array([10, 15, 19]),)

        For a 2D data set the X0 and X1 values are returned:

        >>> load_arrays(2, [10,15,12,19], [12,14,10,17], [4,5,9,-2], Data2D)
        >>> get_indep(2)
        (array([10, 15, 12, 19]), array([12, 14, 10, 17]))

        For PHA data sets the return value is in channel units:

        >>> load_pha('spec', 'src.pi')
        >>> set_analysis('spec', 'energy')
        >>> (chans,) = get_indep('spec')
        >>> chans[0:6]
        array([ 1.,  2.,  3.,  4.,  5.,  6.])

        If the `filter` flag is set then the return will be limited to
        the data that is used in the fit:

        >>> notice_id('spec', 0.5, 7)
        >>> (nchans,) = get_indep('spec', filter=True)
        >>> nchans[0:5]
        array([ 35.,  36.,  37.,  38.,  39.])

        For images the pixel coordinates of each pixel are returned,
        as a 1D array.

        >>> load_image('img', 'image.fits')
        >>> (xvals,yvals) = get_indep('img')
        >>> xvals.shape
        (65536,)
        >>> yvals.shape
        (65536,)
        >>> xvals[0:5]
        array([ 1.,  2.,  3.,  4.,  5.])
        >>> yvals[0:5]
        array([ 1.,  1.,  1.,  1.,  1.])

        The coordinate system for image axes is determinated by the
        `set_coord` setting for the data set:

        >>> set_coord('img', 'physical')
        >>> (avals,bvals) = get_indep('img')
        >>> avals[0:5]
        array([  16.5,   48.5,   80.5,  112.5,  144.5])

        """
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        return d.get_indep(filter=filter)


    def get_axes(self, id=None, bkg_id=None):
        """Return information about the independent axes of a data set.

        get_axes

        SYNOPSIS
           Get the alternate grid of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           bkg_id     - background id
                        default = None

        Returns:
           Array of the alternate independent variable

        DESCRIPTION
           Get the data set alternate grid by data id or bkg_id.
           For PHA spectra, this cooresponds to E_MIN,E_MAX.
           for images, this respresents the axes lengths.

        EXAMPLE
           get_axes()

           get_axes(1)

           get_axes(1, 2)

        SEE ALSO

        """
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        if isinstance(d, sherpa.astro.data.DataPHA):
            return d._get_ebins(group=False)
        elif isinstance(d, (sherpa.data.Data2D,
                            sherpa.astro.data.DataIMG)):
            return d.get_axes()
        return d.get_indep()


    # DOC-NOTE: also in sherpa.utils
    ### Ahelp ingest: 2015-05-05 DJB
    def get_dep(self, id=None, filter=False, bkg_id=None):
        """Return the dependent axis of a data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.
        bkg_id : int or str, optional
           Set if the values returned should be from the given
           background component, instead of the source data set.

        Returns
        -------
        axis : array
           The dependent axis values. The model estimate is compared
           to these values during fitting. For PHA data sets, the 
           return array will match the grouping scheme applied to
           the data set.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the data set does not exist.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        get_rate : Return the count rate of a PHA data set.
        list_data_ids : List the identifiers for the loaded data sets.

        Examples
        --------

        >>> load_arrays(1, [10,15,19], [4,5,9], Data1D)
        >>> get_dep()
        array([10, 15, 19])

        >>> load_arrays(2, [10,15,12,19], [12,14,10,17], [4,5,9,-2], Data2D)
        >>> get_dep(2)
        array([4, 5, 9, -2])

        If the `filter` flag is set then the return will be limited to
        the data that is used in the fit:

        >>> notice_id('spec', 0.5, 7)
        >>> yall = get_dep('spec')
        >>> yfilt = get_dep('spec', filter=True)
        >>> yall.size
        1024
        >>> yfilt.size
        446

        For images, the data is returned as a one-dimensional array:

        >>> load_image('img', 'image.fits')
        >>> ivals = get_dep('img')
        >>> ivals.shape
        (65536,)

        """
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        dep = d.get_y(filter)
        if isinstance(d, sherpa.astro.data.DataPHA):
            old = d._rate
            d._rate=False  # return predicted counts, not rate for PHA
            dep = d.get_y(filter)
            d._rate=old
        return dep


    get_counts = get_dep

    ### Ahelp ingest: 2015-05-01 DJB
    def get_rate(self, id=None, filter=False, bkg_id=None):
        """Return the count rate of a PHA data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.
        bkg_id : int or str, optional
           Set if the rate should be taken from the background
           associated with the data set.

        Returns
        -------
        rate : array
           The rate array. The output matches the grouping of the data
           set. The units are controlled by the `set_analysis` setting
           for this data set; that is, the units used in `plot_data`.
           The return array will match the grouping scheme applied to
           the data set.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain PHA data.

        See Also
        --------
        get_dep : Return the data for a data set.
        ignore : Exclude data from the fit.
        notice : Include data in the fit.
        plot_data : Plot the data values.
        set_analysis : Set the units used when fitting and displaying spectral data.

        Examples
        --------

        >>> rate = get_rate()

        The rate of data set 2 will be in units of count/s/Angstrom
        and only cover the range 20 to 22 Angstroms:

        >>> set_analysis(2, 'wave')
        >>> notice_id(2, 20, 22)
        >>> r2 = get_rate(2, filter=True)

        """
        d = self._get_pha_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        old = d._rate
        d._rate=True     # return count rate for PHA
        rate = d.get_y(filter)
        d._rate=old
        return rate

    ### Ahelp ingest: 2015-05-02 DJB
    ### DOC-TODO: how to get the corresponding x bins for this data?
    ###           i.e. what are the X values for these points
    def get_specresp(self, id=None, filter=False, bkg_id=None):
        """Return the effective area values for a PHA data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the ARF or not. The default is `False`.
        bkg_id : int or str, optional
           Set if the ARF should be taken from a background set
           associated with the data set.

        Returns
        -------
        arf : array
           The effective area values for the data set (or background
           component).

        """
        d = self._get_pha_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        return d.get_specresp(filter)


    ### Ahelp ingest: 2015-05-02 DJB
    def get_exposure(self, id=None, bkg_id=None):
        """Return the exposure time of a PHA data set.

        The exposure time of a PHA data set is taken from the
        EXPTIME keyword in its header, but it can be changed
        once the file has been loaded.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        bkg_id : int or str, optional
           Set to identify which background component to use.  The
           default value (`None`) means that the time is for the
           source component of the data set.

        Returns
        -------
        exposure : number
           The exposure time, in seconds.

        See Also
        --------
        get_areascal : Return the area scaling of a PHA data set.
        get_backscal : Return the area scaling of a PHA data set.
        set_exposure : Change the exposure time of a PHA data set.

        """
        if bkg_id is not None:
            return self.get_bkg(id,bkg_id).exposure
        return self._get_pha_data(id).exposure


    def get_backscal(self, id=None, bkg_id=None):
        """
        get_backscal

        SYNOPSIS
           Get the source or background extraction region areas by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           bkg_id     - background data id
                        default = None

        Returns:
           None

        DESCRIPTION
           Get the extraction region areas of a source PHA dataset by data id
           or of a background dataset by bkg_id. Backscale can be defined as the
           ratio of the area of the source (or background) extraction region in
           image pixels to the total number of image pixels. The fact that
           there is no ironclad definition for this quantity does not matter so
           long as the backscale for a source dataset and its associated
           background dataset are defined in the similar manner, because only
           the ratio of source and background backscale is used in analyses.

        EXAMPLE
           get_backscal()

           get_backscal(1)

           get_backscal(1, 2)

        SEE ALSO
           set_exposure, set_areascal,
           get_exposure, get_areascal
        """
        if bkg_id is not None:
            return self.get_bkg(id,bkg_id).backscal
        return self._get_pha_data(id).backscal


    def get_bkg_scale(self, id=None):
        """
        get_bkg_scale

        SYNOPSIS
           Get the PHA background scale factor by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

        Returns:
           None

        DESCRIPTION
           Get the background scaling factor for a source PHA dataset by data id.
           The background source model is scaled by this factor summed in quadrature
           of all revelent backgrounds.

           scale =  BACKSCAL*EXPOSURE/N  *  \sum_i^N 1./BBACKSCAL_i/BEXPOSURE_i
           

        EXAMPLE
           get_bkg_scale()

           get_bkg_scale(1)

        SEE ALSO
           set_exposure, set_areascal,
           get_exposure, get_areascal
        """
        scale = self._get_pha_data(id).get_background_scale()
        if scale is None:
            raise DataErr('nobkg', self._fix_id(id))

        return scale


    def get_areascal(self, id=None, bkg_id=None):
        """Return the area scaling of a PHA data set.

        get_areascal

        SYNOPSIS
           Get the source or background fractional area by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           bkg_id     - background data id
                        default = None

        Returns:
           None

        DESCRIPTION
           Get the fractional area of a source PHA dataset by data id or of a 
           background data by bkg_id.

        EXAMPLE
           get_areascal()

           get_areascal(1)

           get_areascal(1, 2)

        SEE ALSO
           set_backscal, set_exposure,
           get_backscal, get_exposure
        """
        if bkg_id is not None:
            return self.get_bkg(id,bkg_id).areascal
        return self._get_pha_data(id).areascal


    def _save_type(self, objtype, id, filename, bkg_id=None, **kwargs):
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id, bkg_id)

        if type(d) in (sherpa.astro.data.DataIMG, sherpa.astro.data.DataIMGInt,
                       sherpa.data.Data2D, sherpa.data.Data2DInt):
            backup = d.y
            if objtype == 'delchi':
                raise AttributeError("save_delchi() does not apply for images")

            imgtype = getattr(self, 'get_' + objtype + '_image', None)
            if imgtype is None:
                raise AttributeError("'get_%s_image()' not found" % objtype)

            obj = imgtype(id)

            try:
                # write out the array using the source dataset,
                # include WCS etc.
                d.y = obj.y.ravel()

                sherpa.astro.io.write_image(filename, d, **kwargs)
            finally:
                d.y = backup
            return

        funcname = 'get_'
        if bkg_id is not None:
            funcname += 'bkg_'

        plottype = getattr(self, funcname + objtype + '_plot', None)
        if plottype is None:
            raise AttributeError("'%s%s_plot()' not found" % (funcname,
                                                                objtype))

        obj = plottype(id)
        if bkg_id is not None:
            obj = plottype(id, bkg_id=bkg_id)

        args = None
        fields = None

#        if type(d) in (sherpa.data.Data1DInt, sherpa.astro.data.DataPHA):
#            args = [obj.xlo, obj.xhi, obj.y]
#            fields = ["XLO", "XHI", str(objtype).upper()]
        if (type(d) is sherpa.astro.data.DataPHA and
            objtype in ('model', 'source')):
            args = [obj.xlo, obj.xhi, obj.y]
            fields = ["XLO", "XHI", str(objtype).upper()]
        else:
            args = [obj.x, obj.y]
            fields = ["X", str(objtype).upper()]

        sherpa.astro.io.write_arrays(filename, args, fields, **kwargs)

# To fix bug report 13536, save many kinds of data to ASCII by default,
# and let user override if they want FITS (or vice-versa).  The new defaults 
# as of CIAO 4.6 are:
#
# ascii = False (i.e., write to FITS):
#
# save_pha
# save_image
# save_table
#
# ascii = True (i.e., write to ASCII):
# 
# save_arrays
# save_source
# save_model
# save_resid
# save_delchi
# save_filter
# save_staterror
# save_syserror
# save_error
# save_grouping
# save_quality
# save_data
#
# My logic is that in the former group, you are focused on writing to a
# specific kind of file (PHA, image or table) so use FITS by default.
# In the latter group, you are not focused on file type, but on some
# data or model attribute you want to write out and read into some other
# program with ease.  ASCII is probably better for that.
# SMD 05/15/13
#
    def save_arrays(self, filename, args, fields=None, ascii=True,
                    clobber=False):
        """
        save_arrays

        SYNOPSIS
           Write a list of arrays to file as columns

        SYNTAX

        Arguments:
           filename   - filename with path

           args       - list of arrays that correspond to columns

           fields     - list of column headings
                        default = None

           ascii      - boolean indicating use of an ASCII output format
                        default = True

           clobber    - clobber the existing output file
                        default = False

        Returns:
           None

        DESCRIPTION
           Write a list of arrays to file as columns.

        EXAMPLE

           save_arrays("foo.dat", [a,b,c], fields=['a','b','c'])
       

        SEE ALSO
           save_image, save_data, save_table, save_source, save_model,
           save_resid, save_delchi
        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        sherpa.astro.io.write_arrays(filename, args, fields, ascii, clobber)

    def save_source(self, id, filename=None, bkg_id=None, ascii=False,
                    clobber=False):
        """
        save_source

        SYNOPSIS
           Write the unconvolved source model to file

        SYNTAX

        Arguments:
           id         - data id
                        default = default data id

           filename   - filename with path

           bkg_id     - background id
                        default = default background id

           ascii      - boolean indicating use of an ASCII output format
                        default = False

           clobber    - clobber the existing output file
                        default = False

        Returns:
           None

        DESCRIPTION
           Write the unconvolved source model or background source model 
           to file by data id or background id.  NOTE that the source 
           model array written to file ignores the filter.

        EXAMPLE

           save_source("source.dat", ascii=True)

           save_source("source.fits")

           save_source("bkg.dat", ascii=True, bkg_id=1)

           save_source("bkg.fits", bkg_id=1)

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_model,
           save_resid, save_delchi
        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        self._save_type('source', id, filename, ascii=ascii, clobber=clobber,
                        bkg_id=bkg_id)

    def save_model(self, id, filename=None, bkg_id=None, ascii=False,
                   clobber=False):
        """
        save_model

        SYNOPSIS
           Write the convolved source model to file

        SYNTAX

        Arguments:
           id         - data id
                        default = default data id

           filename   - filename with path

           bkg_id     - background id
                        default = default background id

           ascii      - boolean indicating use of an ASCII output format
                        default = False

           clobber    - clobber the existing output file
                        default = False

        Returns:
           None

        DESCRIPTION
           Write the convolved source model or background model to file by
           data id or background id.  NOTE that the source model array 
           written to file respects the filter.  For PHA data, the source 
           model array is noticed and ungrouped.

        EXAMPLE

           save_model("model.dat", ascii=True)

           save_model("model.fits")

           save_model("bkg_model.dat", ascii=True, bkg_id=1)

           save_model("bkg_model.fits", bkg_id=1)

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_resid, save_delchi
        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        self._save_type('model', id, filename, ascii=ascii, clobber=clobber,
                        bkg_id=bkg_id)

    def save_resid(self, id, filename=None, bkg_id=None, ascii=False,
                   clobber=False):
        """
        save_resid

        SYNOPSIS
           Write the data - model residuals to file

        SYNTAX

        Arguments:
           id         - data id
                        default = default data id

           filename   - filename with path

           bkg_id     - background id
                        default = default background id

           ascii      - boolean indicating use of an ASCII output format
                        default = False

           clobber    - clobber the existing output file
                        default = False

        Returns:
           None

        DESCRIPTION
           Write the data - model or bkg - bkg_model residuals to file by data
           or background id.  NOTE that the residuals array written to file 
           respects the filter and/or grouping flags.

        EXAMPLE

           save_resid("resid.dat", ascii=True)

           save_resid("resid.fits")

           save_resid("bkg_resid.dat", ascii=True, bkg_id=1)

           save_resid("bkg_resid.fits", bkg_id=1)

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_delchi
        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        self._save_type('resid', id, filename, ascii=ascii, clobber=clobber,
                        bkg_id=bkg_id)

    def save_delchi(self, id, filename=None, bkg_id=None, ascii=True,
                    clobber=False):
        """
        save_delchi

        SYNOPSIS
           Write the delta chi squared residuals to file

        SYNTAX

        Arguments:
           id         - data id
                        default = default data id

           filename   - filename with path

           bkg_id     - background id
                        default = default background id

           ascii      - boolean indicating use of an ASCII output format
                        default = False

           clobber    - clobber the existing output file
                        default = False

        Returns:
           None

        DESCRIPTION
           Write the source or background delta chi squared residuals to file
           by data id or background id.  NOTE that the delta chi squared 
           residuals array written to file respects the filter and/or grouping
           flags.

        EXAMPLE

           save_delchi("delchi.dat", ascii=True)

           save_delchi("delchi.fits")

           save_delchi("bkg_delchi.dat", ascii=True, bkg_id=1)

           save_delchi("bkg_delchi.fits", bkg_id=1)

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_resid
        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        self._save_type('delchi', id, filename, ascii=ascii, clobber=clobber,
                        bkg_id=bkg_id)


    def save_filter(self, id, filename=None, bkg_id=None, ascii=True,
                    clobber=False):
        """
        save_filter

        SYNOPSIS
           Write the data set filter to file

        SYNTAX

        Arguments:
           id         - data id
                        default = default data id

           filename   - filename with path

           bkg_id     - background data id
                        default = default background data id

           ascii      - boolean indicating use of an ASCII output format
                        default = False

           clobber    - clobber the existing output file
                        default = False

        Returns:
           None

        DESCRIPTION
           Write the data set filter to file by id or background id. NOTE that
           the filter array written to file respects the grouping flags, if any.

        EXAMPLE

           save_filter("filter.fits")

           save_filter("bkgfilter.fits", bkg_id=1)

           save_filter("filter.dat", ascii=True)

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_delchi
        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id, bkg_id)
        id = self._fix_id(id)
        if d.mask is False:
            raise sherpa.utils.err.DataErr('notmask')
        if not numpy.iterable(d.mask):
            raise sherpa.utils.err.DataErr('nomask', id)

        x = d.get_indep(filter=False)[0]
        mask = numpy.asarray(d.mask, numpy.int)
        if isinstance(d, sherpa.astro.data.DataPHA):
            x = d._get_ebins(group=True)[0]

        self.save_arrays(filename, [x, mask], ['X', 'FILTER'],
                         ascii, clobber)


    # DOC-NOTE: also in sherpa.utils with a different interface
    ### Ahelp ingest: 2015-05-06 DJB
    def save_staterror(self, id, filename=None, bkg_id=None, ascii=True,
                    clobber=False):
        """Save the statistical errors to a file.

        If the statistical errors have not been set explicitly, then
        the values calculated by the statistic - such as `chi2gehrels`
        or `chi2datavar` - will be used.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to. The format
           is determined by the `ascii` argument.
        bkg_id : int or str, optional
           Set if the background should be written out rather
           than the source.
        ascii : bool, optional
           If `False` then the data is written to a FITS
           format binary table. The default is `True`. The
           exact format of the output file depends on the
           I/O library in use (Crates or AstroPy).
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        load_staterror : Load the statistical errors from a file.
        save_error : Save the errors to a file.
        save_syserror : Save the systematic errors to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns `X` and `STAT_ERR`.

        Examples
        --------

        Write out the statistical errors from the default data set to the
        file 'errs.dat'.

        >>> save_staterror('errs.dat')

        Over-write the file it it already exists, and take the data
        from the data set "jet":

        >>> save_staterror('jet', 'err.out', clobber=True)

        Write the data out in FITS format:

        >>> save_staterror('staterr.fits', ascii=False)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)        
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id, bkg_id)
        id = self._fix_id(id)

        x = d.get_indep(filter=False)[0]
        if isinstance(d, sherpa.astro.data.DataPHA):
            x = d._get_ebins(group=True)[0]
        err = self.get_staterror(id, filter=False, bkg_id=bkg_id)
        self.save_arrays(filename, [x, err], ['X', 'STAT_ERR'],
                         ascii, clobber)


    # DOC-NOTE: also in sherpa.utils with a different interface
    ### Ahelp ingest: 2015-05-06 DJB
    def save_syserror(self, id, filename=None, bkg_id=None, ascii=True,
                    clobber=False):
        """Save the systematic errors to a file.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to. The format
           is determined by the `ascii` argument.
        bkg_id : int or str, optional
           Set if the background should be written out rather
           than the source.
        ascii : bool, optional
           If `False` then the data is written to a FITS
           format binary table. The default is `True`. The
           exact format of the output file depends on the
           I/O library in use (Crates or AstroPy).
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If the data set does not contain any systematic errors.
        sherpa.utils.err.DataErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        load_syserror : Load the systematic errors from a file.
        save_error : Save the errors to a file.
        save_staterror : Save the statistical errors to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns `X` and `SYS_ERR`.

        Examples
        --------

        Write out the systematic errors from the default data set to the
        file 'errs.dat'.

        >>> save_syserror('errs.dat')

        Over-write the file it it already exists, and take the data
        from the data set "jet":

        >>> save_syserror('jet', 'err.out', clobber=True)

        Write the data out in FITS format:

        >>> save_syserror('syserr.fits', ascii=False)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id, bkg_id)
        id = self._fix_id(id)

        x = d.get_indep(filter=False)[0]
        if isinstance(d, sherpa.astro.data.DataPHA):
            x = d._get_ebins(group=True)[0]
        err = self.get_syserror(id, filter=False, bkg_id=bkg_id)
        self.save_arrays(filename, [x, err], ['X', 'SYS_ERR'],
                         ascii, clobber)


    # DOC-NOTE: also in sherpa.utils with a different interface
    ### Ahelp ingest: 2015-05-06 DJB
    def save_error(self, id, filename=None, bkg_id=None, ascii=True,
                    clobber=False):
        """Save the errors to a file.

        The total errors for a data set are the quadrature combination
        of the statistical and systematic errors. The systematic
        errors can be 0. If the statistical errors have not been set
        explicitly, then the values calculated by the statistic - such
        as `chi2gehrels` or `chi2datavar` - will be used.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to. The format
           is determined by the `ascii` argument.
        bkg_id : int or str, optional
           Set if the background should be written out rather
           than the source.
        ascii : bool, optional
           If `False` then the data is written to a FITS
           format binary table. The default is `True`. The
           exact format of the output file depends on the
           I/O library in use (Crates or AstroPy).
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        load_staterror : Load the statistical errors from a file.
        load_syserror : Load the systematic errors from a file.
        save_data : Save the data to a file.
        save_staterror : Save the statistical errors to a file.
        save_syserror : Save the systematic errors to a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The output file contains the columns `X` and `ERR`.

        Examples
        --------

        Write out the errors from the default data set to the file
        'errs.dat'.

        >>> save_error('errs.dat')

        Over-write the file it it already exists, and take the data
        from the data set "jet":

        >>> save_error('jet', 'err.out', clobber=True)

        Write the data out in FITS format:

        >>> save_error('err.fits', ascii=False)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id, bkg_id)
        id = self._fix_id(id)

        x = d.get_indep(filter=False)[0]
        if isinstance(d, sherpa.astro.data.DataPHA):
            x = d._get_ebins(group=True)[0]
        err = self.get_error(id, filter=False, bkg_id=bkg_id)
        self.save_arrays(filename, [x, err], ['X', 'ERR'],
                         ascii, clobber)


    ### Ahelp ingest: 2015-05-01 DJB
    def save_pha(self, id, filename=None, bkg_id=None, ascii=False, clobber=False):
        """Save a PHA data set to a file.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to. The format
           is determined by the `ascii` argument.
        bkg_id : int or str, optional
           Set if the background should be written out rather
           than the source.
        ascii : bool, optional
           If `False` then the data is written to a FITS
           format binary table. The default is `True`. The
           exact format of the output file depends on the
           I/O library in use (Crates or AstroPy).
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain PHA data.
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        load_pha : Load a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Write out the PHA data from the default data set to the
        file 'src.pi':

        >>> save_pha('src.pi')

        Over-write the file it it already exists, and take the data
        from the data set "jet":

        >>> save_pha('jet', 'out.pi', clobber=True)

        Write the data out as an ASCII file:

        >>> save_pha('pi.dat', ascii=True)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        d = self._get_pha_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id, bkg_id)

        sherpa.astro.io.write_pha(filename, d, ascii, clobber)


    ### Ahelp ingest: 2015-04-30 DJB
    ### DOC-TODO: labelling as AstroPy; i.e. assuming conversion
    ###           from PyFITS lands soon.
    def save_grouping(self, id, filename=None, bkg_id=None, ascii=True, clobber=False):
        """Save the grouping scheme to a file.

        The output is a two-column file, containing the channel and
        grouping columns from the data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to. The format
           is determined by the `ascii` argument.
        bkg_id : int or str, optional
           Set if the grouping array should be taken from the
           background associated with the data set.
        ascii : bool, optional
           If `False` then the data is written to a FITS
           format binary table. The default is `True`. The
           exact format of the output file depends on the
           I/O library in use (Crates or AstroPy).
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        get_grouping : Return the gouping array for a PHA data set.
        load_quality : Load the quality array from a file and add to a PHA data set.
        set_grouping : Apply a set of grouping flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The column names are 'CHANNEL' and 'GROUPS'.

        Examples
        --------

        Save the channel and grouping columns from the default data
        set to the file 'group.dat' as an ASCII file:

        >>> save_grouping('group.dat')

        Over-write the 'grp.fits' file, if it exists, and write
        out the grouping data from the 'jet' data set, as a FITS
        format file:

        >>> save_grouping('jet', 'grp.fits', ascii=False, clobber=True)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        id = self._fix_id(id)
        d = self._get_pha_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id, bkg_id)

        if d.grouping is None or not numpy.iterable(d.grouping):
            raise sherpa.utils.err.DataErr('nogrouping', id)

        sherpa.astro.io.write_arrays(filename, [d.channel, d.grouping],
                                     ['CHANNEL', 'GROUPS'], ascii, clobber)


    ### Ahelp ingest: 2015-04-30 DJB
    ### DOC-TODO: labelling as AstroPy; i.e. assuming conversion
    ###           from PyFITS lands soon.
    def save_quality(self, id, filename=None, bkg_id=None, ascii=True, clobber=False):
        """Save the quality array to a file.

        The output is a two-column file, containing the channel and
        quality columns from the data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        filename : str
           The name of the file to write the array to. The format
           is determined by the `ascii` argument.
        bkg_id : int or str, optional
           Set if the quality array should be taken from the
           background associated with the data set.
        ascii : bool, optional
           If `False` then the data is written to a FITS
           format binary table. The default is `True`. The
           exact format of the output file depends on the
           I/O library in use (Crates or AstroPy).
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `filename` already exists and `clobber` is `False`.

        See Also
        --------
        get_quality : Return the quality array for a PHA data set.
        load_quality : Load the quality array from a file and add to a PHA data set.
        set_quality : Apply a set of quality flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `filename` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `filename` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        The quality column is labelled as 'GROUPS' rather than
        'QUALITY'.

        Examples
        --------

        Save the channel and quality columns from the default data
        set to the file 'quality.dat' as an ASCII file:

        >>> save_quality('quality.dat')

        Over-write the 'qual.fits' file, if it exists, and write
        out the quality array from the 'jet' data set, as a FITS
        format file:

        >>> save_quality('jet', 'qual.fits', ascii=False, clobber=True)

        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        id = self._fix_id(id)
        d = self._get_pha_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id, bkg_id)

        if d.quality is None or not numpy.iterable(d.quality):
            raise sherpa.utils.err.DataErr('noquality', id)

        # BUG: TODO: change 'GROUPS' to 'QUALITY'
        sherpa.astro.io.write_arrays(filename, [d.channel, d.quality],
                                     ['CHANNEL', 'GROUPS'], ascii, clobber)


    def save_image(self, id, filename=None, ascii=False, clobber=False):
        """
        save_image

        SYNOPSIS
           Write image data by id

        SYNTAX

        Arguments:
           id         - dataset ID
                        default = default data id

           filename   - filename with path

           ascii      - boolean indicating use of an ASCII output format
                        default = False

           clobber    - clobber the existing output file
                        default = False

        Returns:
           None

        DESCRIPTION
           Write image data to a FITS file or ASCII file from a Sherpa dataset
           by id.

        EXAMPLE

           save_image(1, "img.fits")

           save_image(1, "img.out", ascii=True)

        SEE ALSO
           save_pha, save_data, save_table
        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        
        sherpa.astro.io.write_image(filename, self.get_data(id),
                                    ascii, clobber)


    def save_table(self, id, filename=None, ascii=False, clobber=False):
        """
        save_table

        SYNOPSIS
           Write tabular data by id

        SYNTAX

        Arguments:
           id         - dataset ID
                        default = default data id

           filename   - filename with path

           ascii      - boolean indicating use of an ASCII output format
                        default = False

           clobber    - clobber the existing output file
                        default = False

        Returns:
           None

        DESCRIPTION
           Write tabular data to a FITS file or column-based text file
           from a Sherpa dataset by id.

        EXAMPLE

           save_table(1, "tbl.fits")

           save_table(1, "tbl.out", ascii=True)

        SEE ALSO
           save_pha, save_data, save_image
        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        
        sherpa.astro.io.write_table(filename, self.get_data(id),
                                    ascii, clobber)


    def save_data(self, id, filename=None, bkg_id=None, ascii=True, clobber=False):
        """
        save_data

        SYNOPSIS
           Write a data set to file by id

        SYNTAX

        Arguments:
           id         - dataset ID
                        default = default data id

           filename   - filename with path

           bkg_id     - background data id
                        default = default background data id

           ascii      - boolean indicating use of an ASCII output format
                        default = False

           clobber    - clobber the existing output file
                        default = False

        Returns:
           None

        DESCRIPTION
           Write data to a FITS file or ASCII file from a Sherpa dataset
           by id.

        EXAMPLE

           save_data(1, "pha.fits")
           
           save_data(1, "img.fits")

           save_data(1, "data.out", ascii=True)

        SEE ALSO
           save_image, save_data, save_table, save_pha
        """
        clobber=sherpa.utils.bool_cast(clobber)
        ascii=sherpa.utils.bool_cast(ascii)
        if filename is None:
            id, filename = filename, id
        _check_type(filename, basestring, 'filename', 'a string')
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id, bkg_id)

        try:
            sherpa.astro.io.write_pha(filename, d, ascii, clobber)
        except IOErr:
            try:
                sherpa.astro.io.write_image(filename, d, ascii, clobber)
            except IOErr:
                try:
                    sherpa.astro.io.write_table(filename, d, ascii, clobber)
                except:
                    try:
                        sherpa.io.write_data(filename, d, clobber)
                    except:
                        raise


    ### Ahelp ingest: 2015-05-02 DJB
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    def pack_pha(self, id=None):
        """Convert a PHA data set into a file structure.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.

        Returns
        -------
        pha
           The return value depends on the I/O library in use.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain PHA data.

        See Also
        --------
        load_pha : Load a file as a PHA data set.
        set_data : Set a data set.
        unpack_pha : Create a PHA data structure.

        """
        return sherpa.astro.io.pack_pha(self._get_pha_data(id))
    
    ### Ahelp ingest: 2015-05-02 DJB
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    def pack_image(self, id=None):
        """Convert a data set into an image structure.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.

        Returns
        -------
        img
           The return value depends on the I/O library in use.

        See Also
        --------
        load_image : Load a file as an image data set.
        set_data : Set a data set.
        unpack_image : Create an image data structure.

        """
        return sherpa.astro.io.pack_image(self.get_data(id))

    ### Ahelp ingest: 2015-05-02 DJB
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    def pack_table(self, id=None):
        """Convert a data set into a table structure.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.

        Returns
        -------
        tbl
           The return value depends on the I/O library in use and the
           type of data (such as `Data1D`, `Data2D`).

        See Also
        --------
        load_table : Load a FITS binary file as a data set.
        set_data : Set a data set.
        unpack_table : Unpack a FITS binary file into a data structure.

        """
        return sherpa.astro.io.pack_table(self.get_data(id))


    #def _check_resp_id(id):
    #    if (id is not None) and (not self._valid_id(id)):
    #        raise ArgumentTypeError('response identifiers must be integers ' +
    #                                'or strings')

    ### Ahelp ingest: 2015-05-02 DJB
    #@loggable(with_id=True)
    def get_arf(self, id=None, resp_id=None, bkg_id=None):
        """Return the ARF associated with a PHA data set.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        resp_id : int or str, optional
           The identifier for the ARF within this data set, if there
           are multiple responses.
        bkg_id : int or str, optional
           Set this to return the given background component.

        Returns
        -------
        arf : sherpa.astro.instrument.ARF1D instance
           This is a reference to the ARF, rather than a copy, so that
           changing the fields of the object will change the values in
           the data set.

        See Also
        --------
        fake_pha : Simulate a PHA data set from a model.
        load_arf : Load an ARF from a file and add it to a PHA data set.
        load_pha : Load a file as a PHA data set.
        set_full_model : Define the convolved model expression for a data set.
        set_arf : Set the ARF for use by a PHA data set.
        set_rmf : Set the RMF for use by a PHA data set.
        unpack_arf : Read in an ARF from a file.

        Examples
        --------

        Return the exposure field of the ARF from the default data
        set:

        >>> get_arf().exposure

        Copy the ARF from the default data set to data set `2`:

        >>> arf1 = get_arf()
        >>> set_arf(2, arf1)

        Retrieve the ARF associated to the second background
        component of the 'core' data set:

        >>> bgarf = get_arf('core', 'bkg.arf', bkg_id=2)

        """
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        arf, rmf = data.get_response(resp_id)
        if arf is None:
            raise IdentifierErr('getitem', 'ARF data set',
                                data._fix_response_id(resp_id),
                                'in PHA data set %s has not been set' %
                                str(self._fix_id(id)))

        if isinstance(arf, sherpa.astro.data.DataARF):
            arf = sherpa.astro.instrument.ARF1D(arf, data, rmf)

        return arf


    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: add an example of a grating/multiple response
    def set_arf(self, id, arf=None, resp_id=None, bkg_id=None):
        """Set the ARF for use by a PHA data set.

        Set the effective area curve for a PHA data set, or its
        background.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        arf :
           An ARF, such as returned by `get_arf` or `unpack_arf`.
        resp_id : int or str, optional
           The identifier for the ARF within this data set, if there
           are multiple responses.
        bkg_id : int or str, optional
           Set this to identify the ARF as being for use with the
           background.

        See Also
        --------
        get_arf : Return the ARF associated with a PHA data set.
        load_arf : Load an ARF from a file and add it to a PHA data set.
        load_pha : Load a file as a PHA data set.
        set_full_model : Define the convolved model expression for a data set.
        set_rmf : Set the RMF for use by a PHA data set.
        unpack_arf : Read in an ARF from a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `arf` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `arf` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        If a PHA data set has an associated ARF - either from when the
        data was loaded or explicitly with the `set_arf` function -
        then the model fit to the data will include the efect of the
        ARF when the model is created with `set_model` or
        `set_source`. In this case the `get_source` function returns
        the user model, and `get_model` the model that is fit to the
        data (i.e. it includes any response information; that is the
        ARF and RMF, if set). To include the ARF explicitly, use
        `set_full_model`.

        Examples
        --------

        Copy the ARF from the default data set to data set `2`:

        >>> arf1 = get_arf()
        >>> set_arf(2, arf1)

        Read in an ARF from the file 'bkg.arf' and set it as the
        ARF for the background model of data set "core":

        >>> arf = unpack_arf('bkg.arf')
        >>> set_arf('core', arf, bkg_id=1)

        """
        if arf is None:
            id, arf = arf, id

        # store only the ARF dataset in the PHA response dict
        if type(arf) in (sherpa.astro.instrument.ARF1D,):
            arf = arf._arf
        _check_type(arf, sherpa.astro.data.DataARF, 'arf', 'an ARF data set')

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.set_arf(arf, resp_id)
        # Set units of source dataset from channel to energy
        if data.units == 'channel':
            data._set_initial_quantity()


    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    def unpack_arf(self, arg):
        """Create an ARF data structure.

        Parameters
        ----------
        arg :
           Identify the ARF: a file name, or a data structure
           representing the data to use, as used by the I/O backend in
           use by Sherpa: a `TABLECrate` for crates, as used by CIAO,
           or a list of AstroPy HDU objects.

        Returns
        -------
        arf : sherpa.astro.instrument.ARF1D instance

        See Also
        --------
        get_arf : Return the ARF associated with a PHA data set.
        load_arf : Load an ARF from a file and add it to a PHA data set.
        load_bkg_arf : Load an ARF from a file and add it to the background of a PHA data set.
        load_multi_arfs : Load multiple ARFs for a PHA data set.
        load_pha : Load a file as a PHA data set.
        load_rmf : Load a RMF from a file and add it to a PHA data set.
        set_full_model : Define the convolved model expression for a data set.

        Examples
        --------

        >>> arf1 = unpack_arf("arf1.fits")
        >>> arf2 = unpack_arf("arf2.fits")

        Read in an ARF using Crates:

        >>> acr = pycrates.read_file("src.arf")
        >>> arf = unpack_arf(acr)

        Read in an ARF using AstroPy:

        >>> hdus = astropy.io.fits.open("src.arf")
        >>> arf = unpack_arf(hdus)

        """
        return sherpa.astro.instrument.ARF1D(sherpa.astro.io.read_arf(arg))

    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: add an example of a grating/multiple response
    ### DOC-TODO: how to describe I/O backend support?
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    #@loggable(with_id=True, with_keyword='arg')
    def load_arf(self, id, arg=None, resp_id=None, bkg_id=None):
        """Load an ARF from a file and add it to a PHA data set.

        Load in the effective area curve for a PHA data set, or its
        background. The `load_bkg_arf` function can be used for
        setting most background ARFs.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        arg :
           Identify the ARF: a file name, or a data structure
           representing the data to use, as used by the I/O backend in
           use by Sherpa: a `TABLECrate` for crates, as used by CIAO,
           or a list of AstroPy HDU objects.
        resp_id : int or str, optional
           The identifier for the ARF within this data set, if there
           are multiple responses.
        bkg_id : int or str, optional
           Set this to identify the ARF as being for use with the
           background.

        See Also
        --------
        get_arf : Return the ARF associated with a PHA data set.
        load_bkg_arf : Load an ARF from a file and add it to the background of a PHA data set.
        load_multi_arfs : Load multiple ARFs for a PHA data set.
        load_pha : Load a file as a PHA data set.
        load_rmf : Load a RMF from a file and add it to a PHA data set.
        set_full_model : Define the convolved model expression for a data set.
        set_arf : Set the ARF for use by a PHA data set.
        unpack_arf : Create an ARF data structure.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `arg` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `arg` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        If a PHA data set has an associated ARF - either from when the
        data was loaded or explicitly with the `set_arf` function -
        then the model fit to the data will include the efect of the
        ARF when the model is created with `set_model` or
        `set_source`. In this case the `get_source` function returns
        the user model, and `get_model` the model that is fit to the
        data (i.e. it includes any response information; that is the
        ARF and RMF, if set). To include the ARF explicitly, use
        `set_full_model`.

        Examples
        --------

        Use the contents of the file 'src.arf' as the ARF for the
        default data set.

        >>> load_arf('src.arf')

        Read in an ARF from the file 'bkg.arf' and set it as the
        ARF for the background model of data set "core":

        >>> load_arf('core', 'bkg.arf', bkg_id=1)

        """
        if arg is None:
            id, arg = arg, id
        self.set_arf(id, self.unpack_arf(arg), resp_id, bkg_id)

    def get_bkg_arf(self, id=None):
        """
        get_bkg_arf

        SYNOPSIS
           Return a bkg ARF dataset by data id using default bkg_id and resp_id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

        Returns:
           Sherpa DataARF dataset

        DESCRIPTION
           Return a dataset containing background ancillary response data
           given a data id using the default background id and default
           response id.

        SEE ALSO
           set_arf, unpack_arf, load_bkg_arf
        """
        bkg_id = self._get_pha_data(id).default_background_id
        resp_id = self._get_pha_data(id).primary_response_id
        return self.get_arf(id, resp_id, bkg_id)

    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: how to describe I/O backend support?
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    #@loggable(with_id=True, with_keyword='arg')
    def load_bkg_arf(self, id, arg=None):
        """Load an ARF from a file and add it to the background of a
        PHA data set.

        Load in the ARF to the background of the given data set. It
        is only for use when there is only one background component,
        and one response, for the source. For multiple backgrounds
        or responses, use `load_arf`.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        arg :
           Identify the ARF: a file name, or a data structure
           representing the data to use, as used by the I/O backend in
           use by Sherpa: a `TABLECrate` for crates, as used by CIAO,
           or a list of AstroPy HDU objects.

        See Also
        --------
        load_arf : Load an ARF from a file and add it to a PHA data set.
        load_bkg_rmf : Load a RMF from a file and add it to the background of a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `arg` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `arg` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Use the contents of the file 'bkg.arf' as the ARF for the
        background of the default data set.

        >>> load_bkg_arf('bkg.arf')

        Set 'core_bkg.arf' as the ARF for the background of data set
        'core':

        >>> load_bkg_arf('core', 'core_bkg.arf')

        """
        if arg is None:
            id, arg = arg, id
        bkg_id = self._get_pha_data(id).default_background_id
        resp_id = self._get_pha_data(id).primary_response_id
        self.set_arf(id, self.unpack_arf(arg), resp_id, bkg_id)

    ### Ahelp ingest: 2015-04-29 DJB
    def load_multi_arfs(self, id, filenames, resp_ids=None):
        """Load multiple ARFs for a PHA data set.

        A grating observation - such as a Chandra HETG data set - may
        require multiple responses. This function lets the multiple
        ARFs for such a data set be loaded with one command. The
        `load_arf` function can instead be used to load them in
        individually.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        filenames : iterable of str
           An array of file names.
        resp_ids : iterable of int or str
           The identifiers for the ARF within this data set.
           The length should match the filenames argument.

        See Also
        --------
        load_arf : Load an ARF from a file and add it to a PHA data set.
        load_multi_rmfs : Load multiple RMFs for a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with two arguments, they are assumed to be
        `filenames` and `resp_ids`, and three positional arguments
        means `id`, `filenames`, and `resp_ids`.

        Examples
        --------

        Load two ARFs into the default data set, using response ids of
        `1` and `2` for `m1.arf` and `p1.arf` respectively:

        >>> arfs = ['m1.arf', 'p1.arf']
        >>> load_multi_arfs(arfs, [1,2])

        Load in the ARFs to the data set with the identifier
        `lowstate`:

        >>> load_multi_arfs('lowstate', arfs, [1,2])

        """
##         if type(filenames) not in (list, tuple):
##             raise ArgumentError('Filenames must be contained in a list')
##         if type(resp_ids) not in (list, tuple):
##             raise ArgumentError('Response IDs must be contained in a list')

        if resp_ids is None:
            id, filenames, resp_ids = resp_ids, id, filenames

        filenames = list(filenames)[:]
        resp_ids = list(resp_ids)[:]
        
        if (len(filenames) != len(resp_ids)):
            raise ArgumentErr('multirsp')

        while (len(filenames) > 0):
            filename = filenames.pop(0)
            resp_id = resp_ids.pop(0)
            self.load_arf(id, filename, resp_id)

    ### Ahelp ingest: 2015-05-02 DJB
    #@loggable(with_id=True)
    def get_rmf(self, id=None, resp_id=None, bkg_id=None):
        """Return the RMF associated with a PHA data set.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        resp_id : int or str, optional
           The identifier for the RMF within this data set, if there
           are multiple responses.
        bkg_id : int or str, optional
           Set this to return the given background component.

        Returns
        -------
        rmf : sherpa.astro.instrument.RMF1D instance
           This is a reference to the RMF, rather than a copy, so that
           changing the fields of the object will change the values in
           the data set.

        See Also
        --------
        fake_pha : Simulate a PHA data set from a model.
        load_pha : Load a file as a PHA data set.
        load_rmf : Load a RMF from a file and add it to a PHA data set.
        set_full_model : Define the convolved model expression for a data set.
        set_arf : Set the ARF for use by a PHA data set.
        set_rmf : Set the RMF for use by a PHA data set.
        unpack_rmf : Read in a RMF from a file.

        Examples
        --------

        Copy the RMF from the default data set to data set `2`:

        >>> rmf1 = get_rmf()
        >>> set_rmf(2, rmf1)

        Retrieve the RMF associated to the second background
        component of the 'core' data set:

        >>> bgrmf = get_rmf('core', 'bkg.rmf', bkg_id=2)

        """
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        arf, rmf = data.get_response(resp_id)
        if rmf is None:
            raise IdentifierErr('getitem', 'RMF data set',
                                data._fix_response_id(resp_id),
                                'in PHA data set %s has not been set' %
                                str(self._fix_id(id)))

        if isinstance(rmf, sherpa.astro.data.DataRMF):
            rmf = sherpa.astro.instrument.RMF1D(rmf, data, arf)

        return rmf


    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: add an example of a grating/multiple response
    def set_rmf(self, id, rmf=None, resp_id=None, bkg_id=None):
        """Set the RMF for use by a PHA data set.

        Set the effective area curve for a PHA data set, or its
        background.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        arf :
           An ARF, such as returned by `get_arf` or `unpack_arf`.
        resp_id : int or str, optional
           The identifier for the ARF within this data set, if there
           are multiple responses.
        bkg_id : int or str, optional
           Set this to identify the ARF as being for use with the
           background.

        See Also
        --------
        get_rmf : Return the RMF associated with a PHA data set.
        load_pha : Load a file as a PHA data set.
        load_rmf : Load a RMF from a file and add it to a PHA data set.
        set_full_model : Define the convolved model expression for a data set.
        set_arf : Set the ARF for use by a PHA data set.
        unpack_rmf : Create a RMF data structure.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `rmf` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `rmf` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        If a PHA data set has an associated RMF - either from when the
        data was loaded or explicitly with the `set_rmf` function -
        then the model fit to the data will include the efect of the
        RMF when the model is created with `set_model` or
        `set_source`. In this case the `get_source` function returns
        the user model, and `get_model` the model that is fit to the
        data (i.e. it includes any response information; that is the
        ARF and RMF, if set). To include the RMF explicitly, use
        `set_full_model`.

        Examples
        --------

        Copy the RMF from the default data set to data set `2`:

        >>> rmf1 = get_rmf()
        >>> set_rmf(2, rmf1)

        Read in a RMF from the file 'bkg.rmf' and set it as the
        RMF for the background model of data set "core":

        >>> rmf = unpack_rmf('bkg.rmf')
        >>> set_rmf('core', rmf, bkg_id=1)

        """
        if rmf is None:
            id, rmf = rmf, id

        # store only the RMF dataset in the PHA response dict
        if type(rmf) in (sherpa.astro.instrument.RMF1D,):
                rmf = rmf._rmf
        _check_type(rmf, sherpa.astro.data.DataRMF, 'rmf', 'an RMF data set')

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.set_rmf(rmf, resp_id)
        # Set units of source dataset from channel to energy
        if data.units == 'channel':
            data._set_initial_quantity()


    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    def unpack_rmf(self, arg):
        """Create a RMF data structure.

        Parameters
        ----------
        arg :
           Identify the RMF: a file name, or a data structure
           representing the data to use, as used by the I/O backend in
           use by Sherpa: a `RMFCrateDataset` for crates, as used by
           CIAO, or a list of AstroPy HDU objects.

        Returns
        -------
        rmf : sherpa.astro.instrument.RMF1D instance

        See Also
        --------
        get_rmf : Return the RMF associated with a PHA data set.
        load_arf : Load a RMF from a file and add it to a PHA data set.
        load_bkg_rmf : Load a RMF from a file and add it to the background of a PHA data set.
        load_multi_rmfs : Load multiple RMFs for a PHA data set.
        load_pha : Load a file as a PHA data set.
        load_rmf : Load a RMF from a file and add it to a PHA data set.
        set_full_model : Define the convolved model expression for a data set.

        Examples
        --------

        >>> rmf1 = unpack_rmf("rmf1.fits")
        >>> rmf2 = unpack_rmf("rmf2.fits")

        Read in a RMF using Crates:

        >>> acr = pycrates.read_rmf("src.rmf")
        >>> rmf = unpack_rmf(acr)

        Read in a RMF using AstroPy:

        >>> hdus = astropy.io.fits.open("src.rmf")
        >>> rmf = unpack_rmf(hdus)

        """
        return sherpa.astro.instrument.RMF1D(sherpa.astro.io.read_rmf(arg))

    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: add an example of a grating/multiple response
    ### DOC-TODO: how to describe I/O backend support?
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    #@loggable(with_id=True, with_keyword='arg')
    def load_rmf(self, id, arg=None, resp_id=None, bkg_id=None):
        """Load a RMF from a file and add it to a PHA data set.

        Load in the redistribution matrix function for a PHA data set,
        or its background. The `load_bkg_rmf` function can be used for
        setting most background RMFs.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        arg :
           Identify the RMF: a file name, or a data structure
           representing the data to use, as used by the I/O
           backend in use by Sherpa: a `RMFCrateDataset` for
           crates, as used by CIAO, or an AstroPy `HDUList` object.
        resp_id : int or str, optional
           The identifier for the RMF within this data set, if there
           are multiple responses.
        bkg_id : int or str, optional
           Set this to identify the RMF as being for use with the
           background.

        See Also
        --------
        get_rmf : Return the RMF associated with a PHA data set.
        load_bkg_rmf : Load a RMF from a file and add it to the background of a PHA data set.
        load_arf : Load an ARF from a file and add it to a PHA data set.
        load_multi_rmfs : Load multiple RMFs for a PHA data set.
        load_pha : Load a file as a PHA data set.
        set_full_model : Define the convolved model expression for a data set.
        set_rmf : Load a RMF from a file and add it to a PHA data set.
        unpack_rmf : Read in a RMF from a file.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `arg` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `arg` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        If a PHA data set has an associated RMF - either from when the
        data was loaded or explicitly with the `set_rmf` function -
        then the model fit to the data will include the efect of the
        RMF when the model is created with `set_model` or
        `set_source`. In this case the `get_source` function returns
        the user model, and `get_model` the model that is fit to the
        data (i.e. it includes any response information; that is the
        ARF and RMF, if set). To include the RMF explicitly, use
        `set_full_model`.

        Examples
        --------

        Use the contents of the file 'src.rmf' as the RMF for the
        default data set.

        >>> load_rmf('src.rmf')

        Read in a RMF from the file 'bkg.rmf' and set it as the
        RMF for the background model of data set "core":

        >>> load_rmf('core', 'bkg.rmf', bkg_id=1)

        """
        if arg is None:
            id, arg = arg, id
        self.set_rmf(id, self.unpack_rmf(arg), resp_id, bkg_id)

    def get_bkg_rmf(self, id=None):
        """
        get_bkg_rmf

        SYNOPSIS
           Return a bkg RMF dataset by data id using default bkg_id and resp_id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

        Returns:
           Sherpa DataRMF dataset

        DESCRIPTION
           Return a dataset containing background response matrix data
           given a data id using the default background id and default
           response id.

        SEE ALSO
           set_rmf, unpack_rmf, load_bkg_rmf
        """
        bkg_id = self._get_pha_data(id).default_background_id
        resp_id = self._get_pha_data(id).primary_response_id
        return self.get_rmf(id, resp_id, bkg_id)

    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: how to describe I/O backend support?
    ### DOC-TODO: labelling as AstroPy HDUList; i.e. assuming conversion
    ###           from PyFITS lands soon.
    #@loggable(with_id=True, with_keyword='arg')
    def load_bkg_rmf(self, id, arg=None):
        """Load a RMF from a file and add it to the background of a
        PHA data set.

        Load in the RMF to the background of the given data set. It
        is only for use when there is only one background component,
        and one response, for the source. For multiple backgrounds
        or responses, use `load_rmf`.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        arg :
           Identify the RMF: a file name, or a data structure
           representing the data to use, as used by the I/O
           backend in use by Sherpa: a `RMFCrateDataset` for
           crates, as used by CIAO, or an AstroPy `HDUList` object.

        See Also
        --------
        load_rmf : Load a RMF from a file and add it to a PHA data set.
        load_bkg_arf : Load an ARF from a file and add it to the background of a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `arg` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `arg` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Use the contents of the file 'bkg.rmf' as the RMF for the
        background of the default data set.

        >>> load_bkg_rmf('bkg.rmf')

        Set 'core_bkg.rmf' as the RMF for the background of data set
        'core':

        >>> load_bkg_arf('core', 'core_bkg.rmf')

        """
        if arg is None:
            id, arg = arg, id
        bkg_id = self._get_pha_data(id).default_background_id
        resp_id = self._get_pha_data(id).primary_response_id
        self.set_rmf(id, self.unpack_rmf(arg), resp_id, bkg_id)

    ### Ahelp ingest: 2015-04-29 DJB
    def load_multi_rmfs(self, id, filenames, resp_ids=None):
        """Load multiple RMFs for a PHA data set.

        A grating observation - such as a Chandra HETG data set - may
        require multiple responses. This function lets the multiple
        RMFs for such a data set be loaded with one command. The
        `load_rmf` function can instead be used to load them in
        individually.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        filenames : iterable of str
           An array of file names.
        resp_ids : iterable of int or str
           The identifiers for the RMF within this data set.
           The length should match the filenames argument.

        See Also
        --------
        load_rmf : Load a RMF from a file and add it to a PHA data set.
        load_multi_arfs : Load multiple ARFs for a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with two arguments, they are assumed to be
        `filenames` and `resp_ids`, and three positional arguments
        means `id`, `filenames`, and `resp_ids`.

        Examples
        --------

        Load two RMFs into the default data set, using response ids of
        `1` and `2` for `m1.rmf` and `p1.rmf` respectively:

        >>> rmfs = ['m1.rmf', 'p1.rmf']
        >>> load_multi_rmfs(rmfs, [1,2])

        Load in the RMFs to the data set with the identifier
        `lowstate`:

        >>> load_multi_rmfs('lowstate', rmfs, [1,2])

        """
##         if type(filenames) not in (list, tuple):
##             raise ArgumentError('Filenames must be contained in a list')
##         if type(resp_ids) not in (list, tuple):
##             raise ArgumentError('Response IDs must be contained in a list')

        if resp_ids is None:
            id, filenames, resp_ids = resp_ids, id, filenames

        filenames = list(filenames)[:]
        resp_ids = list(resp_ids)[:]
        
        if (len(filenames) != len(resp_ids)):
            raise ArgumentErr('multirsp')

        while (len(filenames) > 0):
            filename = filenames.pop(0)
            resp_id = resp_ids.pop(0)
            self.load_rmf(id, filename, resp_id)
        
    def get_bkg(self, id=None, bkg_id=None):
        """Return the background for a PHA data set.

        get_bkg

        SYNOPSIS
           Return an background PHA dataset by data id and bkg_id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           bkg_id    - background id, if multiple bkgs exist
                       default = default background id

        Returns:
           Sherpa DataPHA dataset

        DESCRIPTION
           Return a dataset containing background PHA data
           given a data id and a background id.

        SEE ALSO
           set_bkg, unpack_bkg, load_bkg
        """
        data = self._get_pha_data(id)
        bkg = data.get_background(bkg_id)
        if bkg is None:
            raise IdentifierErr('getitem', 'background data set',
                                data._fix_background_id(bkg_id),
                                'in PHA data set %s has not been set' %
                                str(self._fix_id(id)))
        return bkg

    ### Ahelp ingest: 2015-04-30 DJB
    def set_bkg(self, id, bkg=None, bkg_id=None):
        """Set the background for a PHA data set.

        The background can either be fit with a model - using
        `set_bkg_model` - or removed from the data before fitting,
        using `subtract`.

        Parameters
        ----------
        id : int or str, optional
           The data set to use. If not given then the default
           identifier is used, as returned by `get_default_id`.
        bkg :
           A PHA data set, such as returned by `get_data` or
           `unpack_pha`.
        bkg_id : int or str, optional
           The identifier for this background, which is needed if
           there are multiple background estimates for the source.

        See Also
        --------
        get_bkg : Return the background for a PHA data set.
        load_bkg : Load the background from a file and add it to a PHA data set.
        load_pha : Load a file as a PHA data set.
        set_bkg_model : Set the background model expression for a data set.
        subtract : Subtract the background estimate from a data set.
        unpack_pha : Create a PHA data structure.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `bkg` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `bkg` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Examples
        --------

        Copy the background from the default data set to data set `2`:

        >>> bkg1 = get_bkg()
        >>> set_bkg(2, bkg1)

        Read in the PHA data from the file 'bkg.pi' and set it as the
        second background component of data set "core":

        >>> bkg = unpack_pha('bkg.pi')
        >>> set_bkg('core', bkg, bkg_id=2)

        """
        if bkg is None:
            id, bkg = bkg, id
        data = self._get_pha_data(id)
        _check_type(bkg, sherpa.astro.data.DataPHA, 'bkg', 'a PHA data set')
        data.set_background(bkg, bkg_id)
        if bkg.grouping is None:
            bkg.grouping = data.grouping
            bkg.grouped = (bkg.grouping is not None)
        if bkg.quality is None:
            bkg.quality = data.quality

        if bkg.get_response() == (None,None):
            bkg.set_response(*data.get_response())

        if bkg.get_response() != (None,None):
            bkg.units = data.units

        bkg.rate = data.rate
        bkg.plot_fac = data.plot_fac


    ### Ahelp ingest: 2015-04-29 DJB
    def list_bkg_ids(self, id=None):
        """List all the background identifiers for a data set.

        A PHA data set can contain multiple background datasets, each
        identified by an integer or string. This function returns a
        list of these identifiers for a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set to query. If not given then the default
           identifier is used, as returned by `get_default_id`.

        Returns
        -------
        ids : array of int or str
           The identifiers for the backround data sets for the data
           set. In many cases this will just be `[1]`.

        See Also
        --------
        list_response_ids : List all the response identifiers of a data set.
        load_bkg : Load the background of a PHA data set.

        """
        #return self._get_pha_data(id).background_ids
        return self._get_pha_data(id)._backgrounds.keys()

    ### Ahelp ingest: 2015-04-29 DJB
    def list_response_ids(self, id=None, bkg_id=None):
        """List all the response identifiers of a data set.

        A PHA data set can contain multiple responses, that is,
        pairs of ARF and RMF, each of which has an identifier.
        This function returns a list of these identifiers
        for a data set.

        Parameters
        ----------
        id : int or str, optional
           The data set to query. If not given then the default
           identifier is used, as returned by `get_default_id`.
        bkg_id : int or str, optional
           Set this to identify the background component to query.

        Returns
        -------
        ids : array of int or str
           The identifiers for the response information for the data
           set. In many cases this will just be `[1]`.

        See Also
        --------
        list_bkg_ids : List all the background identifiers for a data set.
        load_arf : Load an ARF from a file and add it to a PHA data set.
        load_rmf : Load a RMF from a file and add it to a PHA data set.

        """
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        #return data.response_ids
        return data._responses.keys()

    ### Ahelp ingest: 2015-04-27 DJB
    ### DOC-TODO: docs need to be added to sherpa.astro.data.set_analysis
    ### DOC-TODO: should the arguments be renamed to better match optional
    ###           nature of the routine (e.g. can call set_analysis('energy'))?
    #@loggable(with_id=True, with_keyword='quantity')
    def set_analysis(self, id, quantity=None, type='rate', factor=0):
        """Set the units used when fitting and displaying spectral data.

        The set_analysis command sets the units for spectral
        analysis. Note that in order to change the units of a data set
        from 'channel' to 'energy' or 'wavelength', the appropriate
        ARF and RMF instrument response files must be loaded for that
        data set. The `type` and `factor` arguments control how
        the data is plotted.

        Parameters
        ----------
        id : int or str
           If only one argument is given then this is taken to be the
           quantity argument (in which case, the change is made to
           all data sets). If multiple arguments are given then this
           is the identifier for the data set to change.

        quantity : { 'channel', 'chan', 'bin', 'energy', 'ener', 'wavelength', 'wave' }
           The units to use for the analysis.
        
        type : { 'rate', 'counts' }, optional
           The units to use on the Y axis of plots. The default
           is 'rate'.

        factor : int, optional
           The Y axis of plots is multiplied by Energy^factor or
           Wavelength^factor befire display. The default is `0`.

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the `id` argument is not recognized.

        See Also
        --------
        get_analysis : Return the analysis setting for a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `quantity` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `quantity` parameters,
        respectively.

        Examples
        --------

        Set all loaded data sets to use wavelength for any future
        fitting or display.

        >>> set_analysis('wave')

        Set the data set with an identifier of `2` to use energy
        units.

        >>> set_analysis(2, 'energy')

        Set data set 1 to use channel units. Plots will use a Y
        axis of count/bin rather than the default count/s/bin.

        >>> set_analysis(1, 'bin', 'counts')

        Set data set 1 to use energy units. Plots of this data set
        will display keV on the X axis and counts keV (i.e.
        counts/keV * keV^2) in the Y axis.

        >>> set_analysis(1, 'energy', 'counts', 2)

        """
        if quantity is None:
            id, quantity = quantity, id

        _check_type(quantity, basestring, 'quantity', 'a string')
        _check_type(type, basestring, 'type', 'a string')

        ids = self.list_data_ids()
        if id is not None:
            ids = [id]

        for id in ids:
            self._get_pha_data(id).set_analysis(quantity, type, factor) 


    ### Ahelp ingest: 2015-04-27 DJB
    ### DOC-TODO: docs need to be added to sherpa.astro.data.get_analysis
    def get_analysis(self, id=None):
        """Return the units used when fitting spectral data.

        Parameters
        ----------
        id : int or str, optional
           The data set to query. If not given then the default
           identifier is used, as returned by `get_default_id`.

        Returns
        -------
        quantity : { 'channel', 'energy', 'wavelength' }

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain PHA data.
        sherpa.utils.err.IdentifierErr
           If the `id` argument is not recognized.

        See Also
        --------
        get_default_id : Return the default data set identifier.
        set_analysis : Change the analysis setting.

        """
        return self._get_pha_data(id).get_analysis()

    ### Ahelp ingest: 2015-04-28 DJB
    ### DOC-TODO: docs need to be added to sherpa.astro.data.set_coord
    ### DOC-TODO: how best to document the wcssubs support?
    #@loggable(with_id=True, with_keyword='coord')
    def set_coord(self, id, coord=None):
        """Set the coordinate system to use for image analysis.

        The default coordinate system - that is, the mapping between
        pixel position and coordinate value, for images (2D data sets)
        is 'logical'. This function can change this setting, so that
        model parameters can be fit using other systems. This setting
        is also used by the `notice2d` and `ignore2d` series of
        commands.

        Parameters
        ----------
        id : int or str
           The data set to change. If not given then the default
           identifier is used, as returned by `get_default_id`.

        coord : { 'logical', 'image', 'physical', 'world', 'wcs' }
           The coordinate system to use. The 'image' option is the
           same as 'logical', and 'wcs' the same as 'world'.  The
           'physical' and 'world' options are only available for
           installations built with the optional `wcssubs` package.

        See Also
        --------
        get_coord : Get the coordinate system used for image analysis.
        guess : Estimate the parameter values and ranges given the loaded data.
        ignore2d : Exclude a spatial region from an image.
        notice2d : Include a spatial region of an image.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `coord` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `coord` parameters,
        respectively.

        Any limits or values already set for model parameters, such as
        those made by `guess`, may need to be changed after changing
        the coordinate system.

        The 'logical' system is one in which the center of the
        lower-left pixel has coordinates `(1,1)` and the center of the
        top-right pixel has coordinates `(nx,ny)`, for a `nx`
        (columns) by `ny` (rows) pixel image. The pixels have a side
        of length `1`, so the first pixel covers the range `x=0.5` to
        `x=1.5` and `y=0.5` to `y=1.5`.

        The 'physical' and 'world' coordinate systems rely on FITS
        World Coordinate System (WCS) standard [1]_. The 'physical'
        system refers to a linear transformation, with possible
        offset, of the 'logical' system. The 'world' system refers to
        the mapping to a celestial coordinate system.

        Support for the 'physical' and 'world' coordinate systems
        relies on the optional `wcssubs` package. This can be checked
        for by

        >>> import sherpa.astro.utils._wcs

        References
        ----------

        .. [1] http://fits.gsfc.nasa.gov/fits_wcs.html

        Examples
        --------

        Change the coordinate system of the default data set to
        the world system ('wcs' is a synonym for 'world').

        >>> set_coord('wcs')

        Change the data set with the id of 'm82' to use the
        physical coordinate system.

        >>> set_coord('m82', 'physical')

        """
        if coord is None:
            id, coord = coord, id

        _check_type(coord, basestring, 'coord', 'a string')
        
        ids = self.list_data_ids()
        if id is not None:
            ids = [id]

	if(len(ids)==0):
	    raise IdentifierErr('nodatasets')

        for id in ids:
           self._get_img_data(id).set_coord(coord)


    ### Ahelp ingest: 2015-04-28 DJB
    ### DOC-TODO: docs need to be added to sherpa.astro.data.get_coord
    def get_coord(self, id=None):
        """Get the coordinate system used for image analysis.

        Parameters
        ----------
        id : int or str, optional
           The data set to query. If not given then the default
           identifier is used, as returned by `get_default_id`.

        Returns
        -------
        coord : { 'logical', 'physical', 'world' }

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain image data.
        sherpa.utils.err.IdentifierErr
           If the `id` argument is not recognized.

        See Also
        --------
        get_default_id : Return the default data set identifier.
        set_coord : Set the coordinate system to use for image analysis.

        """
        return self._get_img_data(id).coord


    def ignore_bad(self, id=None, bkg_id=None):
        """Exclude channels marked as bad in a PHA data set.

        ignore_bad

        SYNOPSIS
           Ignore bins according to quality flags

        SYNTAX

        Arguments:
           id          Data set ID
                       default = default data id

           bkg_id    - background id
                       default = default bkg id

        Returns:
           None

        DESCRIPTION
           Ignore bins according to quality flags

               ignore_bad()

               ignore_bad("src")

               notice(0.1,4.0)
               ignore_bad()
           WARNING: filtering with quality flags, noticing all bins

           Ignore_bad may alter the grouping flags size, so any filters
           in place will be removed in the case where the grouping flags
           size is changed.

        SEE ALSO
           notice2d_id, notice2d_image, ignore2d, ignore2d_id, ignore2d_image,
           notice, ignore, notice_id, ignore_id, notice2d
        """
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.ignore_bad()


    def _notice_warning(self):
        quantities = numpy.asarray([data.get_analysis()
                                    for data in self._data.values()
                                    if isinstance(data, 
                                                  sherpa.astro.data.DataPHA)])

        if len(quantities) > 1 and not (quantities==quantities[0]).all():
            warning("not all PHA datasets have equal analysis quantities")


    def notice(self, lo=None, hi=None, **kwargs):

        if lo is not None or hi is not None:
            self._notice_warning()
        sherpa.ui.utils.Session.notice(self, lo, hi, **kwargs)

    notice.__doc__ = sherpa.ui.utils.Session.notice.__doc__


    def ignore(self, lo=None, hi=None, **kwargs):

        if lo is not None or hi is not None:
            self._notice_warning()
        sherpa.ui.utils.Session.ignore(self, lo, hi, **kwargs)

    ignore.__doc__ = sherpa.ui.utils.Session.ignore.__doc__


    def notice2d(self, val=None):
        """Include a spatial region of an image.

        notice2d

        SYNOPSIS
           Notice a region mask for all Sherpa DataIMG datasets

        SYNTAX

        Arguments:
           val       - filename and path of region file or DM region syntax
                       default = None
        Returns:
           None

        DESCRIPTION
           Notice a region mask for all Sherpa DataIMG datasets using a
           DM region library syntax or a region file.

           Example1: notice2d with region file

	       notice2d( 'region filename' )

           Example2: notice2d with DM region syntax in physical coordinates

               notice2d( 'circle(4071, 4250, 135)' )

        SEE ALSO
           notice2d_id, notice2d_image, ignore2d, ignore2d_id, ignore2d_image,
           notice, ignore, notice_id, ignore_id
        """
        for d in self._data.values():
            _check_type(d, sherpa.astro.data.DataIMG, 'img',
                        'a image data set')
            d.notice2d(val, False)

    def ignore2d(self, val=None):
        """Exclude a spatial region from an image.

        ignore2d

        SYNOPSIS
           Ignore a region mask for all Sherpa DataIMG datasets

        SYNTAX

        Arguments:
           val       - filename and path of region file or DM region syntax
                       default = None

        Returns:
           None

        DESCRIPTION
           Ignore a region mask for all Sherpa DataIMG datasets using a
           DM region library syntax or a region file.

           Example1: ignore2d with region file

               ignore2d( 'region filename' )

           Example2: ignore2d with DM region syntax in physical coordinates

               ignore2d( 'circle(4071, 4250, 135)' )

        SEE ALSO
           notice2d_id, notice2d, notice2d_image, ignore2d_id, ignore2d_image,
           notice, ignore, notice_id, ignore_id
        """
        for d in self._data.values():
            _check_type(d, sherpa.astro.data.DataIMG, 'img',
                        'a image data set')
            d.notice2d(val, True)
    
    def notice2d_id(self, ids, val=None):
        """
        notice2d_id

        SYNOPSIS
           Notice a region mask for specific Sherpa DataIMG datasets

        SYNTAX

        Arguments:
           ids       - list of data ids to apply filter

           val       - filename and path of region file or DM region syntax
                       default = None

        Returns:
           None

        DESCRIPTION
           Notice a region mask for specific Sherpa DataIMG datasets by ids
           using a DM region library syntax or a region file.

           Example1: notice2d_id with region file

               notice2d_id(['foo','bar'], 'region filename' )

           Example2: notice2d_id with DM region syntax in physical coordinates

               notice2d_id([2,5,7], 'circle(4071, 4250, 135)' )

        SEE ALSO
           notice2d, notice2d_image, ignore2d, ignore2d_id, ignore2d_image,
           notice, ignore, notice_id, ignore_id
         """
        if self._valid_id(ids):
            ids = (ids,)
        else:
            try:
                ids = tuple(ids)
            except TypeError:
                _argument_type_error('ids',
                                     'an identifier or list of identifiers')


        for id in ids:
            _check_type(self.get_data(id), sherpa.astro.data.DataIMG,
                        'img', 'a image data set')
            self.get_data(id).notice2d(val, False)
        
    def ignore2d_id(self, ids, val=None):
        """
        ignore2d_id

        SYNOPSIS
           Ignore a region mask for specific Sherpa DataIMG datasets

        SYNTAX

        Arguments:
           ids       - list of data ids to apply filter

           val       - filename and path of region file or DM region syntax
                       default = None

        Returns:
           None

        DESCRIPTION
           Ignore a region mask for specific Sherpa DataIMG datasets by ids
           using a DM region library syntax or a region file.

           Example1: ignore2d_id with region file

               ignore2d_id(['foo','bar'], 'region filename' )

           Example2: ignore2d_id with DM region syntax in physical coordinates

               ignore2d_id([2,5,7], 'circle(4071, 4250, 135)' )

        SEE ALSO
           notice2d, ignore2d, notice2d_id, notice2d_image, ignore2d_image,
           notice, ignore, notice_id, ignore_id
        """
        if self._valid_id(ids):
            ids = (ids,)
        else:
            try:
                ids = tuple(ids)
            except TypeError:
                _argument_type_error('ids',
                                     'an identifier or list of identifiers')
                

        for id in ids:
            _check_type(self.get_data(id), sherpa.astro.data.DataIMG,
                        'img', 'a image data set')
            self.get_data(id).notice2d(val, True)

    def notice2d_image(self, ids=None):
        """
        notice2d_image

        SYNOPSIS
           Get the current DS9 region and notice data points
           within the region.

        SYNTAX

        Arguments:
           ids       - list of data ids to apply filter

        Returns:
           None

        DESCRIPTION
           After the user has drawn regions in a DS9 pane, this
           function can be used to retrieve the region description,
           and notice all pixels within the region.  The region
           description is automatically converted to the data's
           current coordinate system.

        SEE ALSO
           notice2d, ignore2d, notice2d_id, ignore2d_id, ignore2d_image,
           notice, ignore, notice_id, ignore_id
        """
        if (ids == None):
            ids = self._default_id
        if self._valid_id(ids):
            ids = (ids,)
        else:
            try:
                ids = tuple(ids)
            except TypeError:
                _argument_type_error('ids',
                                     'an identifier or list of identifiers')
                

        for id in ids:
            _check_type(self.get_data(id), sherpa.astro.data.DataIMG,
                        'img', 'a image data set')
            coord = self.get_coord(id)
            if (coord == 'logical'):
                coord = 'image'
            if (coord == 'world'):
                coord = 'wcs'
            regions = self.image_getregion(coord).replace(';','')
            self.notice2d_id(id, regions)

    def ignore2d_image(self, ids=None):
        """
        ignore2d_image

        SYNOPSIS
           Get the current DS9 region and ignore data points
           within the region.

        SYNTAX

        Arguments:
           ids       - list of data ids to apply filter

        Returns:
           None

        DESCRIPTION
           After the user has drawn regions in a DS9 pane, this
           function can be used to retrieve the region description,
           and ignore all pixels within the region.  The region
           description is automatically converted to the data's
           current coordinate system.

        SEE ALSO
           notice2d, ignore2d, notice2d_id, notice2d_image, ignore2d_id,
           notice, ignore, notice_id, ignore_id
        """
        if (ids == None):
            ids = self._default_id
        if self._valid_id(ids):
            ids = (ids,)
        else:
            try:
                ids = tuple(ids)
            except TypeError:
                _argument_type_error('ids',
                                     'an identifier or list of identifiers')
                

        for id in ids:
            _check_type(self.get_data(id), sherpa.astro.data.DataIMG,
                        'img', 'a image data set')
            coord = self.get_coord(id)
            if (coord == 'logical'):
                coord = 'image'
            if (coord == 'world'):
                coord = 'wcs'
            regions = self.image_getregion(coord).replace(';','')
            self.ignore2d_id(id, regions)


    #@loggable(with_id=True, with_keyword='arg')
    def load_bkg(self, id, arg=None, use_errors=False, bkg_id=None):
        """Load the background from a file and add it to a PHA data set.

        load_bkg

        SYNOPSIS
           Load background PHA data by id

        SYNTAX

        Arguments:
           id         - dataset ID
                        default = default data id

           arg        - filename and path | PHACrate obj | PyFITS HDUList obj

           use_errors - flag to use errors
                        default = False

           bkg_id     - background id, if multiple bkgs exist
                        default = default background id

        Returns:
           None

        DESCRIPTION
           Load background PHA data from a FITS file or a PHACrate object or a
           PyFITS HDUList object into a Sherpa dataset by data id and
           background id.

        SEE ALSO
           load_image, load_arf, load_rmf, load_data, load_table,
           load_pha
        """
        if arg is None:
            id, arg = arg, id

        bkgsets = self.unpack_bkg(arg, use_errors)

        if numpy.iterable(bkgsets):
            for bkgid, bkg in enumerate(bkgsets):
                self.set_bkg(id, bkg, bkgid+1)
        else:
            self.set_bkg(id, bkgsets, bkg_id)

    ### Ahelp ingest: 2015-04-29 DJB
    #@loggable(with_id=True)
    def group(self, id=None, bkg_id=None):
        """Turn on the grouping for a PHA data set.

        A PHA data set can be grouped either because it contains
        grouping information [1]_, which is automatically applied when
        the data is read in with `load_pha` or `load_data`, or because
        the `group` set of routines has been used to dynamically
        re-group the data. The `ungroup` function removes this
        grouping (however it was created). The `group` function
        re-applies this grouping. The grouping scheme can be
        changed dynamically, using the `group_xxx` series of
        routines.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        bkg_id : int or str, optional
           Set to group the background associated with the data set.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.
        sherpa.utils.err.DataErr
           If the data set is already grouped.

        See Also
        --------
        fit : Fit one or more data sets.
        group_adapt : Adaptively group to a minimum number of counts.
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_bins : Group into a fixed number of bins.
        group_counts : Group into a minimum number of counts per bin.
        group_snr : Group into a minimum signal-to-noise ratio.
        group_width : Group into a fixed bin width.
        set_grouping : Apply a set of grouping flags to a PHA data set.
        set_quality : Apply a set of quality flags to a PHA data set.
        ungroup : Turn off the grouping for a PHA data set.

        Notes
        -----
        PHA data is often grouped to improve the signal to noise of
        the data, by decreasing the number of bins, so that a
        chi-square statistic can be used when fitting the data.  After
        calling `group`, anything that uses the data set - such as a
        plot, fit, or error analysis - will use the grouped data
        values. Models should be re-fit if `group` is called; the
        increase in the signal of the bins may mean that a chi-square
        statistic can now be used.

        The grouping is implemented by separate arrays to the main
        data - the information is stored in the `grouping` and
        `quality` arrays of the PHA data set - so that a data set can
        be grouped and ungrouped many times, without losing
        information. The `group` command does not create this
        information; this is either created by modifying the PHA file
        before it is read in, or by using the `group_xxx` routines
        once the data has been loaded.

        The `grouped` field of a PHA data set is set to `True` when
        the data is grouped.

        References
        ----------

        .. [1] Arnaud., K. & George, I., "The OGIP Spectral File
               Format",
               http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html

        Examples
        --------

        Group the data in the default data set:

        >>> group()
        >>> get_data().grouped
        True

        Group the first background component of the 'core' data set:

        >>> group('core', bkg_id=1)
        >>> get_bkg('core', bkg_id=1).grouped
        True

        The data is fit using the ungrouped data, and then plots of
        the data and best-fit, and the residuals, are created. The
        first plot uses the ungrouped data, and the second plot uses
        the grouped data.

        >>> ungroup()
        >>> fit()
        >>> plot_fit_resid()
        >>> group()
        >>> plot_fit_resid()

        """
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        if bkg_id is None:
            # First, group backgrounds associated with the
            # data set ID; report if background(s) already grouped.
            for bid in data.background_ids:
                try:
                    self.group(id, bid)
                except DataErr, e:
                    info(str(e))

            # Now check if data is already grouped, and send error message
            # if so
            if (data.grouped is True):
                raise DataErr('groupset', 'data set', str(self._fix_id(id)), 'True')
        else:
            # Just grouping one particular background here
            if (data.grouped is True):
                raise DataErr('groupset', 'background', str(self._fix_id(bkg_id)), 'True')

        # If we get here, checks showed data not grouped, so set group flag
        data.group()


    ### Ahelp ingest: 2015-04-29 DJB
    def set_grouping(self, id, val=None, bkg_id=None):
        """Apply a set of grouping flags to a PHA data set.

        A group is indicated by a sequence of flag values starting
        with `1` and then `-1` for all the channels in the group,
        following [1]_.  Setting the grouping column automatically
        turns on the grouping flag for that data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        val : array of int
           This must be an array of grouping values of the same length
           as the data array.
        bkg_id : int or str, optional
           Set to group the background associated with the data set.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.

        See Also
        --------
        fit : Fit one or more data sets.
        get_grouping : Return the grouping flags for a PHA data set.
        group : Turn on the grouping for a PHA data set.
        group_adapt : Adaptively group to a minimum number of counts.
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_bins : Group into a fixed number of bins.
        group_counts : Group into a minimum number of counts per bin.
        group_snr : Group into a minimum signal-to-noise ratio.
        group_width : Group into a fixed bin width.
        load_grouping : Load the grouping scheme from a file and add to a PHA data set.
        set_quality : Apply a set of quality flags to a PHA data set.
        ungroup : Turn off the grouping for a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `val` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `val` parameters,
        respectively.

        References
        ----------

        .. [1] Arnaud., K. & George, I., "The OGIP Spectral File
               Format",
               http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html

        Examples
        --------

        Copy the grouping array from data set 2 into the default data
        set:

        >>> grp = get_data(2).grouping
        >>> set_grouping(grp)

        Copy the grouping from data set "src1" to the source and
        background data sets of "src2":

        >>> grp = get_data("src1").grouping
        >>> set_grouping("src2", grp)
        >>> set_grouping("src2", grp, bkg_id=1)

        """
        if val is None:
            id, val = val, id

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        if val is None:
            data.grouping = None
        else:
            if(type(val) in (numpy.ndarray,) and 
               issubclass(val.dtype.type, numpy.integer)):
                data.grouping = numpy.asarray(val)
            else:
                data.grouping = numpy.asarray(val, SherpaInt)


    ### Ahelp ingest: 2015-04-30 DJB
    def get_grouping(self, id=None, bkg_id=None):
        """Return the grouping array for a PHA data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        bkg_id : int or str, optional
           Set if the grouping flags should be taken from a background
           associated with the data set.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.

        See Also
        --------
        fit : Fit one or more data sets.
        get_quality : Return the quality array for a PHA data set.
        ignore_bad : Exclude channels marked as bad in a PHA data set.
        set_grouping : Apply a set of grouping flags to a PHA data set.

        Examples
        --------

        Copy the grouping array from the default data set to data set
        2:

        >>> grp1 = get_grouping()
        >>> set_grouping(2, grp1)

        Return the grouping array of the background component labelled
        2 for the 'histate' data set:

        >>> grp = get_grouping('histate', bkg_id=2)

        """

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        return data.grouping

    ### Ahelp ingest: 2015-04-29 DJB
    def set_quality(self, id, val=None, bkg_id=None):
        """Apply a set of quality flags to a PHA data set.

        A quality value of 1 or greater indicates a good channel,
        otherwise the channel is considered bad and can be
        excluded using the `ignore_bad` function, as discussed
        in [1]_.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        val : array of int
           This must be an array of quality values of the same length
           as the data array.
        bkg_id : int or str, optional
           Set if the quality values should be associated with the
           background associated with the data set.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.

        See Also
        --------
        fit : Fit one or more data sets.
        get_quality : Return the quality array for a PHA data set.
        ignore_bad : Exclude channels marked as bad in a PHA data set.
        load_quality : Load the quality array from a file and add to a PHA data set.
        set_grouping : Apply a set of grouping flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `val` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `val` parameters,
        respectively.

        References
        ----------

        .. [1] Arnaud., K. & George, I., "The OGIP Spectral File
               Format",
               http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html

        Examples
        --------

        Copy the quality array from data set 2 into the default data
        set, and then ensure that any 'bad' channels are ignored:

        >>> qual = get_data(2).quality
        >>> set_quality(qual)
        >>> ignore_bad()

        Copy the quality array from data set "src1" to the source and
        background data sets of "src2":

        >>> qual = get_data("src1").quality
        >>> set_quality("src2", qual)
        >>> set_quality("src2", qual, bkg_id=1)

        """
        if val is None:
            id, val = val, id

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        if val is None:
            data.quality = None
        else:
            if(type(val) in (numpy.ndarray,) and 
               issubclass(val.dtype.type, numpy.integer)):
                data.quality = numpy.asarray(val)
            else:
                data.quality = numpy.asarray(val, SherpaInt)

    ### DOC TODO: Need to document that routines like get_quality return
    ###           a reference to the data - so can change the data structure
    ###           - and not a copy

    ### DOC-TODO: explain that many of these can be done with
    ###           direct object access
    ###           get_data().exposure [= ...]

    ### Ahelp ingest: 2015-04-30 DJB
    def get_quality(self, id=None, bkg_id=None):
        """Return the quality flags for a PHA data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        bkg_id : int or str, optional
           Set if the quality flags should be taken from a background
           associated with the data set.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.

        See Also
        --------
        fit : Fit one or more data sets.
        get_grouping : Return the grouping array for a PHA data set.
        ignore_bad : Exclude channels marked as bad in a PHA data set.
        set_quality : Apply a set of quality flags to a PHA data set.

        Examples
        --------

        Copy the quality array from the default data set to data set
        2:

        >>> qual1 = get_quality()
        >>> set_quality(2, qual1)

        Return the quality array of the background component labelled
        2 for the 'histate' data set:

        >>> qual = get_quality('histate', bkg_id=2)

        """

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        return data.quality

    ### Ahelp ingest: 2015-04-29 DJB
    #@loggable(with_id=True)
    def ungroup(self, id=None, bkg_id=None):
        """Turn off the grouping for a PHA data set.

        A PHA data set can be grouped either because it contains
        grouping information [1]_, which is automatically applied when
        the data is read in with `load_pha` or `load_data`, or because
        the `group_xxx` set of routines has been used to dynamically
        re-group the data. The `ungroup` function removes this
        grouping (however it was created).

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        bkg_id : int or str, optional
           Set to ungroup the background associated with the data set.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.
        sherpa.utils.err.DataErr
           If the data set is not grouped.

        See Also
        --------
        fit : Fit one or more data sets.
        group : Turn on the grouping for a PHA data set.

        Notes
        -----
        PHA data is often grouped to improve the signal to noise of
        the data, by decreasing the number of bins, so that a
        chi-square statistic can be used when fitting the data.  After
        calling `ungroup`, anything that uses the data set - such as a
        plot, fit, or error analysis - will use the original data
        values. Models should be re-fit if `ungroup` is called; this
        may require a change of statistic depending on the counts per
        channel in the spectrum.

        The grouping is implemented by separate arrays to the main
        data - the information is stored in the `grouping` and
        `quality` arrays of the PHA data set - so that a data set
        can be grouped and ungrouped many times, without losing
        information.

        The `grouped` field of a PHA data set is set to `False` when
        the data is not grouped.

        If subtracting the background estimate from a data set, the
        grouping applied to the source data set is used for both
        source and background data sets.

        References
        ----------

        .. [1] Arnaud., K. & George, I., "The OGIP Spectral File
               Format",
               http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html

        Examples
        --------

        Ungroup the data in the default data set:

        >>> ungroup()
        >>> get_data().grouped
        False

        Ungroup the first background component of the 'core' data set:

        >>> ungroup('core', bkg_id=1)
        >>> get_bkg('core', bkg_id=1).grouped
        False

        """
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        if bkg_id is None:
            # First, ungroup backgrounds associated with the
            # data set ID; report if background(s) already ungrouped.
            for bid in data.background_ids:
                try:
                    self.ungroup(id, bid)
                except DataErr, e:
                    info(str(e))

            # Now check if data is already ungrouped, and send error message
            # if so
            if (data.grouped is False):
                raise DataErr('groupset', 'data set', str(self._fix_id(id)), 'False')
        else:
            if (data.grouped is False):
                # Just ungrouping one particular background here
                raise DataErr('groupset', 'background', str(self._fix_id(bkg_id)), 'False')

        # If we get here, checks showed data grouped, so set ungroup flag
        data.ungroup()

    ### DOC-TODO: need to document somewhere that this ignores existing
    ###           quality flags and how to use tabStops to include
    ###           this information
    ### DOC-TODO: how to set the quality if using tabstops to indicate
    ###           "bad" channels, rather than ones to ignore

    ### Ahelp ingest: 2015-04-30 DJB
    #@loggable(with_id=True, with_keyword='num')
    def group_bins(self, id, num=None, bkg_id=None, tabStops=None):
        """Group into a fixed number of bins.

        Combine the data so that there `num` equal-width bins (or
        groups). The binning scheme is applied to all the channels,
        but any existing filter - created by the `ignore` or `notice`
        set of functions - is re-applied after the data has been
        grouped.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        num : int
           The number of bins in the grouped data set. Each bin
           will contain the same number of channels.
        bkg_id : int or str, optional
           Set to group the background associated with the data set.
           When `bkg_id` is None (which is the default), the
           grouping is applied to all the associated background
           data sets as well as the source data set.
        tabStops : array of int or bool, optional
           If set, indicate one or more ranges of channels that should
           not be included in the grouped output. The array should
           match the number of channels in the data set and non-zero or
           `True` means that the channel should be ignored from the
           grouping (use 0 or `False` otherwise).

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.

        See Also
        --------
        group_adapt : Adaptively group to a minimum number of counts.
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_counts : Group into a minimum number of counts per bin.
        group_snr : Group into a minimum signal-to-noise ratio.
        group_width : Group into a fixed bin width.
        set_grouping : Apply a set of grouping flags to a PHA data set.
        set_quality : Apply a set of quality flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `num` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `num` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Unlike `group`, it is possible to call `group_bins` multiple
        times on the same data set without needing to call `ungroup`.

        Since the bin width is an integer number of channels, it is
        likely that some channels will be "left over". This is even
        more likely when the `tabStops` parameter is set. If this
        happens, a warning message will be displayed to the screen and
        the quality value for these channels will be set to 2. This
        information can be found with the `get_quality` command.

        Examples
        --------

        Group the default data set so that there are 50 bins.

        >>> group_bins(50)

        Group the 'jet' data set to 50 bins and plot the result,
        then re-bin to 100 bins and overplot the data:

        >>> group_bins('jet', 50)
        >>> plot_data('jet')
        >>> group_bins('jet', 100)
        >>> plot_data('jet', overplot=True)

        The grouping is applied to the full data set, and then
        the filter - in this case defined over the range 0.5
        to 8 keV - will be applied. This means that the
        noticed data range will likely contain less than
        50 bins.

        >>> set_analysis('energy')
        >>> notice(0.5, 8)
        >>> group_bins(50)
        >>> plot_data()

        Do not group any channels numbered less than 20 or
        800 or more. Since there are 780 channels to be
        grouped, the width of each bin will be 20 channels
        and there are no "left over" channels:

        >>> notice()
        >>> channels = get_data().channel
        >>> ign = (channels <= 20) | (channels >= 800)
        >>> group_bins(39, tabStops=ign)
        >>> plot_data()

        """
        if num is None:
            id, num = num, id

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_bins(num, tabStops)

    ### Ahelp ingest: 2015-04-30 DJB
    ### DOC-TODO: should num= be renamed val= to better match
    ###           underlying code/differ from group_bins?
    #@loggable(with_id=True, with_keyword='num')
    def group_width(self, id, num=None, bkg_id=None, tabStops=None):
        """Group into a fixed bin width.

        Combine the data so that each bin contains `num` channels.
        The binning scheme is applied to all the channels, but any
        existing filter - created by the `ignore` or `notice` set of
        functions - is re-applied after the data has been grouped.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        num : int
           The number of channels to combine into a group.
        bkg_id : int or str, optional
           Set to group the background associated with the data set.
           When `bkg_id` is None (which is the default), the
           grouping is applied to all the associated background
           data sets as well as the source data set.
        tabStops : array of int or bool, optional
           If set, indicate one or more ranges of channels that should
           not be included in the grouped output. The array should
           match the number of channels in the data set and non-zero or
           `True` means that the channel should be ignored from the
           grouping (use 0 or `False` otherwise).

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.

        See Also
        --------
        group_adapt : Adaptively group to a minimum number of counts.
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_bins : Group into a fixed number of bins.
        group_counts : Group into a minimum number of counts per bin.
        group_snr : Group into a minimum signal-to-noise ratio.
        set_grouping : Apply a set of grouping flags to a PHA data set.
        set_quality : Apply a set of quality flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `num` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `num` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Unlike `group`, it is possible to call `group_width` multiple
        times on the same data set without needing to call `ungroup`.

        Unless the requested bin width is a factor of the number of
        channels (and no `tabStops` parameter is given), then some
        channels will be "left over". If this happens, a warning
        message will be displayed to the screen and the quality value
        for these channels will be set to 2. This information can be
        found with the `get_quality` command.

        Examples
        --------

        Group the default data set so that each bin contains 20
        channels:

        >>> group_width(20)

        Plot two versions of the 'jet' data set: the first uses
        20 channels per group and the second is 50 channels per
        group:

        >>> group_width('jet', 20)
        >>> plot_data('jet')
        >>> group_width('jet', 50)
        >>> plot_data('jet', overplot=True)

        The grouping is applied to the full data set, and then
        the filter - in this case defined over the range 0.5
        to 8 keV - will be applied.

        >>> set_analysis('energy')
        >>> notice(0.5, 8)
        >>> group_width(50)
        >>> plot_data()

        The grouping is not applied to channels 101 to
        149, inclusive:

        >>> notice()
        >>> channels = get_data().channel
        >>> ign = (channels > 100) & (channels < 150)
        >>> group_width(40, tabStops=ign)
        >>> plot_data()

        """
        if num is None:
            id, num = num, id

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_width(num, tabStops)

    ### Ahelp ingest: 2015-04-30 DJB
    #@loggable(with_id=True, with_keyword='num')
    def group_counts(self, id, num=None, bkg_id=None,
                     maxLength=None, tabStops=None):
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
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        num : int
           The number of channels to combine into a group.
        bkg_id : int or str, optional
           Set to group the background associated with the data set.
           When `bkg_id` is None (which is the default), the
           grouping is applied to all the associated background
           data sets as well as the source data set.
        maxLength : int, optional
           The maximum number of channels that can be combined into a
           single group.
        tabStops : array of int or bool, optional
           If set, indicate one or more ranges of channels that should
           not be included in the grouped output. The array should
           match the number of channels in the data set and non-zero or
           `True` means that the channel should be ignored from the
           grouping (use 0 or `False` otherwise).

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.

        See Also
        --------
        group_adapt : Adaptively group to a minimum number of counts.
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_bins : Group into a fixed number of bins.
        group_snr : Group into a minimum signal-to-noise ratio.
        group_width : Group into a fixed bin width.
        set_grouping : Apply a set of grouping flags to a PHA data set.
        set_quality : Apply a set of quality flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `num` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `num` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Unlike `group`, it is possible to call `group_counts` multiple
        times on the same data set without needing to call `ungroup`.

        If channels can not be placed into a "valid" group, then a
        warning message will be displayed to the screen and the
        quality value for these channels will be set to 2. This
        information can be found with the `get_quality` command.

        Examples
        --------

        Group the default data set so that each bin contains at
        least 20 counts:

        >>> group_counts(20)

        Plot two versions of the 'jet' data set: the first uses
        20 counts per group and the second is 50:

        >>> group_counts('jet', 20)
        >>> plot_data('jet')
        >>> group_counts('jet', 50)
        >>> plot_data('jet', overplot=True)

        The grouping is applied to the full data set, and then
        the filter - in this case defined over the range 0.5
        to 8 keV - will be applied.

        >>> set_analysis('energy')
        >>> notice(0.5, 8)
        >>> group_counts(30)
        >>> plot_data()

        If a channel has more than 30 counts then do not group,
        otherwise group channels so that they contain at least 40
        counts. The `group_adapt` and `group_adapt_snr` functions
        provide similar functionality to this example.  A maximum
        length of 10 channels is enforced, to avoid bins getting too
        large when the signal is low.

        >>> notice()
        >>> counts = get_data().counts
        >>> ign = counts > 30
        >>> group_counts(40, tabStops=ign, maxLength=10)

        """
        if num is None:
            id, num = num, id

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_counts(num, maxLength, tabStops)

    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: check the Poisson stats claim; I'm guessing it means
    ###           gaussian (i.e. sqrt(n))
    #@loggable(with_id=True, with_keyword='snr')
    def group_snr(self, id, snr=None, bkg_id=None,
                        maxLength=None, tabStops=None, errorCol=None):
        """Group into a minimum signal-to-noise ratio.

        Combine the data so that each bin has a signal-to-noise ratio
        of at least `snr`. The binning scheme is applied to all the
        channels, but any existing filter - created by the `ignore` or
        `notice` set of functions - is re-applied after the data has
        been grouped.  The background is *not* included in this
        calculation; the calculation is done on the raw data even if
        `subtract` has been called on this data set.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        snr : number
           The minimum signal-to-noise ratio that must be reached
           to form a group of channels.
        bkg_id : int or str, optional
           Set to group the background associated with the data set.
           When `bkg_id` is None (which is the default), the
           grouping is applied to all the associated background
           data sets as well as the source data set.
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

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.

        See Also
        --------
        group_adapt : Adaptively group to a minimum number of counts.
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_bins : Group into a fixed number of bins.
        group_counts : Group into a minimum number of counts per bin.
        group_width : Group into a fixed bin width.
        set_grouping : Apply a set of grouping flags to a PHA data set.
        set_quality : Apply a set of quality flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `snr` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `snr` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Unlike `group`, it is possible to call `group_snr` multiple
        times on the same data set without needing to call `ungroup`.

        If channels can not be placed into a "valid" group, then a
        warning message will be displayed to the screen and the
        quality value for these channels will be set to 2. This
        information can be found with the `get_quality` command.

        Examples
        --------

        Group the default data set so that each bin has a 
        signal-to-noise ratio of at least 5:

        >>> group_snr(20)

        Plot two versions of the 'jet' data set: the first uses
        a signal-to-noise ratio of 3 and the second 5:

        >>> group_snr('jet', 3)
        >>> plot_data('jet')
        >>> group_snr('jet', 5)
        >>> plot_data('jet', overplot=True)

        """
        if snr is None:
            id, snr = snr, id
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_snr(snr, maxLength, tabStops, errorCol)

    ### Ahelp ingest: 2015-05-01 DJB
    #@loggable(with_id=True, with_keyword='min')
    def group_adapt(self, id, min=None, bkg_id=None,
                     maxLength=None, tabStops=None):
        """Adaptively group to a minimum number of counts.

        Combine the data so that each bin contains `min` or more
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
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        min : int
           The number of channels to combine into a group.
        bkg_id : int or str, optional
           Set to group the background associated with the data set.
           When `bkg_id` is None (which is the default), the
           grouping is applied to all the associated background
           data sets as well as the source data set.
        maxLength : int, optional
           The maximum number of channels that can be combined into a
           single group.
        tabStops : array of int or bool, optional
           If set, indicate one or more ranges of channels that should
           not be included in the grouped output. The array should
           match the number of channels in the data set and non-zero or
           `True` means that the channel should be ignored from the
           grouping (use 0 or `False` otherwise).

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.

        See Also
        --------
        group_adapt_snr : Adaptively group to a minimum signal-to-noise ratio.
        group_bins : Group into a fixed number of bins.
        group_counts : Group into a minimum number of counts per bin.
        group_snr : Group into a minimum signal-to-noise ratio.
        group_width : Group into a fixed bin width.
        set_grouping : Apply a set of grouping flags to a PHA data set.
        set_quality : Apply a set of quality flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `min` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `min` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Unlike `group`, it is possible to call `group_adapt` multiple
        times on the same data set without needing to call `ungroup`.

        If channels can not be placed into a "valid" group, then a
        warning message will be displayed to the screen and the
        quality value for these channels will be set to 2. This
        information can be found with the `get_quality` command.

        Examples
        --------

        Group the default data set so that each bin contains at
        least 20 counts:

        >>> group_adapt(20)

        Plot two versions of the 'jet' data set: the first uses
        an adaptive scheme of 20 counts per bin, the second
        the `group_counts` method:

        >>> group_adapt('jet', 20)
        >>> plot_data('jet')
        >>> group_counts('jet', 20)
        >>> plot_data('jet', overplot=True)

        """
        if min is None:
            id, min = min, id
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_adapt(min, maxLength, tabStops)

    ### Ahelp ingest: 2015-05-01 DJB
    ### DOC-TODO: shouldn't this be snr=None rather than min=None
    #@loggable(with_id=True, with_keyword='min')
    def group_adapt_snr(self, id, min=None, bkg_id=None,
                        maxLength=None, tabStops=None, errorCol=None):
        """Adaptively group to a minimum signal-to-noise ratio.

        Combine the data so that each bin has a signal-to-noise ratio
        of at least `num`. The difference to `group_snr` is that this
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
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.
        num : number
           The minimum signal-to-noise ratio that must be reached
           to form a group of channels.
        bkg_id : int or str, optional
           Set to group the background associated with the data set.
           When `bkg_id` is None (which is the default), the
           grouping is applied to all the associated background
           data sets as well as the source data set.
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

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.

        See Also
        --------
        group_adapt : Adaptively group to a minimum number of counts.
        group_bins : Group into a fixed number of bins.
        group_counts : Group into a minimum number of counts per bin.
        group_snr : Group into a minimum signal-to-noise ratio.
        group_width : Group into a fixed bin width.
        set_grouping : Apply a set of grouping flags to a PHA data set.
        set_quality : Apply a set of quality flags to a PHA data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `num` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `num` parameters,
        respectively. The remaining parameters are expected to be
        given as named arguments.

        Unlike `group`, it is possible to call `group_adapt_snr`
        multiple times on the same data set without needing to call
        `ungroup`.

        If channels can not be placed into a "valid" group, then a
        warning message will be displayed to the screen and the
        quality value for these channels will be set to 2. This
        information can be found with the `get_quality` command.

        Examples
        --------

        Group the default data set so that each bin contains 
        a signal-to-noise ratio of at least 5:

        >>> group_adapt_snr(5)

        Plot two versions of the 'jet' data set: the first uses an
        adaptive scheme and the second the non-adaptive version:

        >>> group_adapt_snr('jet', 4)
        >>> plot_data('jet')
        >>> group_snr('jet', 4)
        >>> plot_data('jet', overplot=True)

        """
        if min is None:
            id, min = min, id
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_adapt_snr(min, maxLength, tabStops, errorCol)

    ### Ahelp ingest: 2015-04-28 DJB
    #@loggable(with_id=True)
    def subtract(self, id=None):
        """Subtract the background estimate from a data set.

        The `subtract` function performs a channel-by-channel
        subtraction of the background estimate from the data. After
        this command, anything that uses the data set - such as a
        plot, fit, or error analysis - will use the subtracted
        data. Models should be re-fit if `subtract` is called.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.
        sherpa.utils.err.DataErr
           If the data set is already subtracted.

        See Also
        --------
        fit : Fit one or more data sets.
        unsubtract : Undo any background subtraction for the data set.

        Notes
        -----
        Unlike X-Spec [1]_, Sherpa does not automatically subtract
        the background estimate from the data.

        Background subtraction can only be performed when data and
        background are of the same length.  If the data and background
        are ungrouped, both must have same number of channels.  If
        they are grouped, data and background can start with different
        numbers of channels, but must have the same number of groups
        after grouping.

        The equation for the subtraction is:

           src_counts - bg_counts * (src_exposure * src_backscal)
                                    -----------------------------
                                     (bg_exposure * bg_backscal)

        where src_exposure and bg_exposure are the source and
        background exposure times, and src_backscal and bg_backscal
        are the source and background backscales.  The backscale, read
        from the `BACKSCAL` header keyword of the PHA file [2]_, is
        the ratio of data extraction area to total detector area.

        The `subtracted` field of a dataset is set to `True` when
        the background is subtracted.

        References
        ----------

        .. [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XspecSpectralFitting.html

        .. [2] https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node5.html

        Examples
        --------

        Background subtract the default data set.

        >>> subtract()
        >>> get_data().subtracted
        True

        Remove the background from the data set labelled 'src':

        >>> subtract('src')
        >>> get_data('src').subtracted
        True

        """
        if (self._get_pha_data(id).subtracted is True):
            raise DataErr('subtractset', 'data set', str(self._fix_id(id)), 'True')
        self._get_pha_data(id).subtract()

    ### Ahelp ingest: 2015-04-28 DJB
    #@loggable(with_id=True)
    def unsubtract(self, id=None):
        """Undo any background subtraction for the data set.

        The `unsubtract` function undoes any changes made by
        `subtract`. After this command, anything that uses the data
        set - such as a plot, fit, or error analysis - will use the
        original data values. Models should be re-fit if `subtract` is
        called.

        Parameters
        ----------
        id : int or str, optional
           The identifier for the data set to use. If not given then
           the default identifier is used, as returned by
           `get_default_id`.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain a PHA data set.
        sherpa.utils.err.DataErr
           If the data set does not have its background subtracted.

        See Also
        --------
        fit : Fit one or more data sets.
        subtract : Subtract the background estimate from a data set.

        Notes
        -----
        The `subtracted` field of a PHA data set is set to `False`
        when the background is not subtracted.

        Examples
        --------

        Remove the background subtraction from the default data set.

        >>> subtract()
        >>> get_data().subtracted
        False

        Remove the background subtraction from the data set labelled
        'src':

        >>> subtract('src')
        >>> get_data('src').subtracted
        False

        """
        if (self._get_pha_data(id).subtracted is False):
            raise DataErr('subtractset', 'data set', str(self._fix_id(id)), 'False')
        self._get_pha_data(id).unsubtract()

    def fake_pha(self, id, arf, rmf, exposure, backscal=None, areascal=None,
                 grouping=None, grouped=False, quality=None, bkg=None):
        """Simulate a PHA data set from a model.

        fake_pha

        SYNOPSIS
           Create and fill a Sherpa DataPHA dataset by data id 
           with faked PHA counts using poisson noise.

        SYNTAX

        Arguments:
           id        - data id, if exists overwrites old dataset

           arf       - Sherpa DataARF dataset, defines ancillary response

           rmf       - Sherpa DataRMF dataset, defines response matrix

           exposure  - length of observation in seconds

           backscal  - background scaling factor
                       default = None

           areascal  - area scaling factor
                       default = None

           grouping  - integer array of grouping flags
                       default = None

           grouped   - dataset grouped boolean
                       default = False

           quality   - integer array of quality flags
                       default = None

           bkg       - python DataPHA object defines the background,
                       default = None

        Returns:
           None

        DESCRIPTION
           fake_pha allows for the simulation of spectra given a source model
           and a grid.  The generated counts will contain poisson noise. If the
           data id exists, the dataset's counts will be clobber, if not, a new
           dataset with that data id will be generated.

        SEE ALSO
           save_pha           
        """
        d = sherpa.astro.data.DataPHA('', None, None)
        if self._data.has_key(id):
            d = self._get_pha_data(id)
        else:
            # Make empty header OGIP compliant
            # And add appropriate values to header from input values
            d.header = dict(HDUCLASS="OGIP", HDUCLAS1="SPECTRUM", HDUCLAS2="TOTAL", HDUCLAS3="TYPE:I", HDUCLAS4="COUNT", HDUVERS="1.1.0")
            self.set_data(id, d)

        if rmf is None:
            raise DataErr('normffake', id)

        if(type(rmf) in (str, numpy.string_)):
            if os.path.isfile(rmf):
                rmf = self.unpack_rmf(rmf)
            else:
                raise IOErr("filenotfound", rmf);

        if(arf is not None and type(arf) in (str, numpy.string_)):
            if os.path.isfile(arf):
                arf = self.unpack_arf(arf)
            else:
                raise IOErr("filenotfound", arf);

        d.exposure = exposure

        if backscal is not None:
            d.backscal = backscal

        if areascal is not None:
            d.areascal = areascal

        if quality is not None:
            d.quality = quality

        if d.channel is None:
            d.channel = sao_arange(1, rmf.detchans)

        else:
            if len(d.channel) != rmf.detchans:
                raise DataErr('incompatibleresp', rmf.name, str(id))

        if grouping is not None:
            d.grouping = grouping

        if d.grouping is not None:
            if sherpa.utils.bool_cast(grouped):
                d.group()
            else:
                d.ungroup()

        for resp_id in d.response_ids:
            d.delete_response(resp_id)

        if arf is not None:
            self.set_arf(id, arf)

        self.set_rmf(id, rmf)

        # Update background here.  bkg contains a new background;
        # delete the old background (if any) and add the new background
        # to the simulated data set, BEFORE simulating data, and BEFORE
        # adding scaled background counts to the simulated data.
        if bkg is not None:
            for bkg_id in d.background_ids:
                d.delete_background(bkg_id)
            self.set_bkg(id, bkg)

        # Calculate the source model, and take a Poisson draw based on
        # the source model.  That becomes the simulated data.
        m = self.get_model(id)
        d.counts = sherpa.utils.poisson_noise( d.eval_model(m) )

        # Add in background counts:
        #  -- Scale each background properly given data's
        #     exposure time, BACKSCAL and AREASCAL
        #  -- Take average of scaled backgrounds
        #  -- Take a Poisson draw based on the average scaled background
        #  -- Add that to the simulated data counts
        #
        # Adding background counts is OPTIONAL, only done if user sets
        # "bkg" argument to fake_pha.  The reason is that the user could
        # well set a "source" model that does include a background
        # component.  In that case users should have the option to simulate
        # WITHOUT background counts being added in.
        #
        # If bkg is not None, then backgrounds were previously updated
        # above, so it is OK to use "bkg is not None" as the condition
        # here.
        if bkg is not None:
            nbkg = len(d.background_ids)
            b=0
            for bkg_id in d.background_ids:
                b=b + d.get_background_scale() * d.get_background(bkg_id).counts
                
            if (nbkg > 0):
                b = b / nbkg
                b_poisson = sherpa.utils.poisson_noise(b)
                d.counts = d.counts + b_poisson

        d.name = 'faked'


    ###########################################################################
    # PSF
    ###########################################################################
    #@loggable
    def load_psf(self, modelname, filename_or_model, *args, **kwargs):    
        kernel = filename_or_model
        if isinstance(filename_or_model, basestring):
            try:
                kernel = self._eval_model_expression(filename_or_model)
            except:
                try:
                    kernel = self.unpack_data(filename_or_model,
                                              *args, **kwargs)
                except:
                    raise

        psf = sherpa.astro.instrument.PSFModel(modelname, kernel)
	if isinstance(kernel, sherpa.models.Model): 
	    self.freeze(kernel)
        self._add_model_component(psf)
        self._psf_models.append(psf)


    load_psf.__doc__ = sherpa.ui.utils.Session.load_psf.__doc__


    ###########################################################################
    # Models
    ###########################################################################

    # DOC-NOTE: also in sherpa.utils
    #@loggable(with_id=True, with_keyword='model')
    def set_full_model(self, id, model=None):
        sherpa.ui.utils.Session.set_full_model(self, id, model)

        if model is None:
            id, model = model, id

        data = self.get_data(id)
        if isinstance(data, sherpa.astro.data.DataPHA):
            model = self._get_model(id)

            if data._responses:

                instruments = (sherpa.astro.instrument.RSPModel,
                               sherpa.astro.instrument.RMFModel,
                               sherpa.astro.instrument.ARFModel,
                               sherpa.astro.instrument.MultiResponseSumModel,
                               sherpa.astro.instrument.PileupRMFModel)

                do_warning = True
                #if type(model) in instruments:
                #if isinstance(model, instruments):
                if sherpa.ui.utils._is_subclass(type(model), instruments):
                    do_warning = False
                for part in model:
                    #if type(part) in instruments:
                    #if isinstance(part, instruments):
                    if sherpa.ui.utils._is_subclass(type(part), instruments):
                        do_warning = False
                if do_warning:
                    warning("PHA source model '%s' \ndoes not" %
                            model.name +
                            " have an associated instrument model; " +
                            "consider using \nset_source() instead of" +
                            " set_full_model() to include associated " +
                            "\ninstrument automatically")

    set_full_model.__doc__ = sherpa.ui.utils.Session.set_full_model.__doc__


    def _add_convolution_models(self, id, data, model, is_source):
        model = \
            sherpa.ui.utils.Session._add_convolution_models(self, id, data,
                                                            model, is_source)
        id = self._fix_id(id)
        if (isinstance(data, sherpa.astro.data.DataPHA) and is_source):
            if not data.subtracted:
                bkg_srcs = self._background_sources.get(self._fix_id(id),{})
                if len(bkg_srcs.keys()) != 0:
                    model = (model +
                             sherpa.astro.background.BackgroundSumModel
                             (data, bkg_srcs))

            pileup_model = self._pileup_models.get(self._fix_id(id))
            if pileup_model is not None:
                resp = \
                    sherpa.astro.instrument.PileupResponse1D(data, pileup_model)
                model = resp(model)

            elif len(data._responses) > 1:
                resp = sherpa.astro.instrument.MultipleResponse1D(data)
                model = resp(model)

            else:
                resp = sherpa.astro.instrument.Response1D(data)
                model = resp(model)

        return model

    #@loggable(with_id=True)
    def get_response(self, id=None, bkg_id=None):
        """
        get_response

        SYNOPSIS
           Return a PHA instrument response, multiple PHA instrument response
           or PHA pileup response model by data id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id
           
           bkg_id    - background id
                       default = default bkg_id

        Returns:
           Sherpa response model

        DESCRIPTION
           Retrieve a PHA instrument response, multiple PHA instrument response
           or PHA pileup response model by data id and background id

        SEE ALSO
           get_pileup_model, get_arf, get_rmf
        """
        pha = self._get_pha_data(id)
        if bkg_id is not None:
            pha = self.get_bkg(id, bkg_id)
        resp = None

        pileup_model = self._pileup_models.get(self._fix_id(id))
        if pileup_model is not None:
            resp = sherpa.astro.instrument.PileupResponse1D(pha, pileup_model)
        elif len(pha._responses) > 1:
            resp = sherpa.astro.instrument.MultipleResponse1D(pha)
        else:
            resp = sherpa.astro.instrument.Response1D(pha)

        return resp


    def get_pileup_model(self, id=None):
        """
        get_pileup_model

        SYNOPSIS
           Return a jdpileup model by data id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

        Returns:
           Sherpa jdpileup model

        DESCRIPTION
           Retrieve a previously set jdpileup model by data id

        SEE ALSO
           set_pileup_model, jdpileup
        """
        return self._get_item(id, self._pileup_models, 'pileup model',
                              'has not been set')

    #@loggable(with_id=True, with_keyword='model')
    def set_pileup_model(self, id, model=None):
        """
        set_pileup_model

        SYNOPSIS
           Set a jdpileup model by data id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           model     - Sherpa jdpileup model

        Returns:
           None

        DESCRIPTION
           Put a Sherpa jdpileup model on the stack by data id

        SEE ALSO
           jdpileup, get_pileup_model
        """
        if model is None:
            id, model = model, id
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)
        self._set_item(id, model, self._pileup_models, sherpa.models.Model,
                       'model', 'a model object or model expression string')


    def _get_bkg_model_status(self, id=None, bkg_id=None):
        src = self._background_sources.get(id, {}).get(bkg_id)
        mdl = self._background_models.get(id, {}).get(bkg_id)

        if src is None and mdl is None:
            IdentifierErr('getitem', 'model', id, 'has not been set')

        model = mdl
        is_source = False

        if mdl is None and src is not None:
            is_source = True
            model = src

        return (model, is_source)


    def get_bkg_source(self, id=None, bkg_id=None):
        """
        get_bkg_source

        SYNOPSIS
           Return the background unconvolved model by data id and bkg id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           bkg_id    - bkg id, if multiple bkgs exist
                       default = default bkg id

        Returns:
           Sherpa bkg unconvolved model

        DESCRIPTION
           Retrieve a Sherpa unconvolved background model by data id and
           background id.

        SEE ALSO
           set_bkg_model, delete_bkg_model
        """
        id     = self._fix_id(id)
        bkg_id = self._fix_id(bkg_id)

        model = self._background_sources.get(id, {}).get(bkg_id)
        if model is None:
            raise ModelErr('nobkg', bkg_id, id)

        return model


    def get_bkg_model(self, id=None, bkg_id=None):
        """
        get_bkg_model

        SYNOPSIS
           Return the background convolved model by data id and bkg id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           bkg_id    - bkg id, if multiple bkgs exist
                       default = default bkg id

        Returns:
           Sherpa bkg convolved model

        DESCRIPTION
           Retrieve a Sherpa convolved background model by data id and
           background id.

        SEE ALSO
           set_bkg_model, delete_bkg_model
        """
        id     = self._fix_id(id)
        bkg_id = self._fix_id(bkg_id)
        src, is_source = self._get_bkg_model_status(id, bkg_id)

        if src is None:
            raise ModelErr('nobkg', bkg_id, id)

        data = self._get_pha_data(id)
        bkg = self.get_bkg(id, bkg_id)

        model = src
        if is_source:
            if len(bkg.response_ids)!=0:
                resp = sherpa.astro.instrument.Response1D(bkg)
                model = resp(src)
            else:
                resp = sherpa.astro.instrument.Response1D(data)
                model = resp(src)
        return model

    #@loggable(with_id=True, with_keyword='model')
    def set_bkg_full_model(self, id, model=None, bkg_id=None):
        """
        set_bkg_full_model

        SYNOPSIS
           Set a convolved Sherpa background model by data id
           and bkg id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           model     - Sherpa bkg model

           bkg_id    - bkg id, if multiple bkgs exist
                       default = default bkg id

        Returns:
           None

        DESCRIPTION
           Add a Sherpa background convolved model to the list of
           current background models by data id and background id.

        SEE ALSO
           get_bkg_model, delete_bkg_model, set_bkg_model, 
           set_bkg_source
        """
        if model is None:
            id, model = model, id

        id     = self._fix_id(id)
        bkg_id = self._fix_id(bkg_id)

        if isinstance(model, basestring):
            model = self._eval_model_expression(model)
        _check_type(model, sherpa.models.Model, 'model',
                    'a model object or model expression string')

        self._background_models.setdefault(id, {})[bkg_id] = model

        data = self.get_bkg(id, bkg_id)
        if data.units != 'channel' and data._responses:

            instruments = (sherpa.astro.instrument.RSPModel,
                           sherpa.astro.instrument.RMFModel,
                           sherpa.astro.instrument.ARFModel,
                           sherpa.astro.instrument.MultiResponseSumModel,
                           sherpa.astro.instrument.PileupRMFModel)

            do_warning = True
            #if type(model) in instruments:
            #if isinstance(model, instruments):
            if sherpa.ui.utils._is_subclass(type(model), instruments):
                do_warning = False
            for part in model:
                #if type(part) in instruments:
                #if isinstance(part, instruments):
                if sherpa.ui.utils._is_subclass(type(part), instruments):
                    do_warning = False
            if do_warning:
                self.delete_bkg_model(id,bkg_id)
                raise TypeError("PHA background source model '%s' \n" % model.name +
                                " does not have an associated instrument model;" +
                                " consider using\n set_bkg_source() instead of" +
                                " set_bkg_model() to include associated\n instrument" +
                                " automatically")

        self._runparamprompt(model.pars)

    ### Ahelp ingest: 2015-04-29 DJB
    ### DOC-TODO: should probably explain more about how backgrounds are fit?
    #@loggable(with_id=True, with_keyword='model')
    def set_bkg_model(self, id, model=None, bkg_id=None):
        """Set the background model expression for a data set.

        The background emission can be fit by a model, defined by the
        `set_bkg_model` call, rather than subtracted from the data.
        If the background is subtracted then the background model is
        ignored when fitting the data.

        Parameters
        ----------
        id : int or str, optional
           The data set containing the source expression. If not given
           then the default identifier is used, as returned by
           `get_default_id`.
        model : str or sherpa.models.Model object
           This defines the model used to fit the data. It can be a
           Python expression or a string version of it.
        bkg_id : int or str, optional
           The identifier for the background of the data set, in
           cases where multiple backgrounds are provided.

        See Also
        --------
        delete_model : Delete the model expression from a data set.
        fit : Fit one or more data sets.
        integrate1d : Integrate 1D source expressions.
        set_model : Set the model expression for a data set.
        set_bkg_full_model : Define the convolved background model expression for a data set.
        show_bkg_model : Display the background model expression for a data set.

        Notes
        -----
        The function does not follow the normal Python standards for
        parameter use, since it is designed for easy interactive use.
        When called with a single un-named argument, it is taken to be
        the `model` parameter. If given two un-named arguments, then
        they are interpreted as the `id` and `model` parameters,
        respectively.

        The emission defined by the background model expression is
        included in the fit to the source dataset, scaling by exposure
        time and area size (given by the ratio of the background to
        source `BACKSCAL` values). That is, if `src_model` and
        `bkg_model` represent the source and background model
        expressions set by calls to `set_model` and `set_bkg_model`
        respectively, the source data is fit by `src_model + scale *
        bkg_model`, where `scale` is the scaling factor.

        PHA data sets will automatically apply the instrumental
        response (ARF and RMF) to the background expression. For some
        cases this is not useful - for example, when different
        responses should be applied to different model components - in
        which case `set_bkg_full_model` should be used instead.

        Examples
        --------

        The background is model by a gaussian line (`gauss1d` model
        component called `bline`) together with an absorbed polynomial
        (the `bgnd` component). The absorbing component (`gal`) is
        also used in the source expression.

        >>> set_model(xsphabs.gal*powlaw1d.pl)
        >>> set_bkg_model(gauss1d.bline + gal*polynom1d.bgnd)

        In this example, the default data set has two background
        estimates, so models are set for both components. The same
        model is applied to both, except that the relative
        normalisations are allowed to vary (by inclusion of the
        `scale` component).

        >>> bmodel = xsphabs.gabs * powlaw1d.pl
        >>> set_bkg_model(2, bmodel)
        >>> set_bkg_model(2, bmodel * const1d.scale, bkg_id=2)

        """
        if model is None:
            id, model = model, id

        id     = self._fix_id(id)
        bkg_id = self._fix_id(bkg_id)

        if isinstance(model, basestring):
            model = self._eval_model_expression(model)
        _check_type(model, sherpa.models.Model, 'model',
                    'a model object or model expression string')

        self._background_sources.setdefault(id, {})[bkg_id] = model

        self._runparamprompt(model.pars)

        # Delete any previous model set with set_full_bkg_model()
        bkg_mdl = self._background_models.get(id, {}).pop(bkg_id, None)
        if bkg_mdl is not None:
            warning("Clearing background convolved model\n'%s'\n" %
                    (bkg_mdl.name) + "for dataset %s background %s" %
                    (str(id), str(bkg_id)))


    set_bkg_source = set_bkg_model


    def delete_bkg_model(self, id=None, bkg_id=None):
        """
        delete_bkg_model

        SYNOPSIS
           Remove a bkg model by data id and bkg id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           bkg_id    - bkg id, if multiple bkgs exist
                       default = default bkg id

        Returns:
           None

        DESCRIPTION
           Removes a background model from the stack by data id
           and background id.

        SEE ALSO
           get_bkg_model, set_bkg_model
        """
        id     = self._fix_id(id)
        bkg_id = self._fix_id(bkg_id)
        # remove dependency of having a loaded PHA dataset at the time
        # of bkg model init.
#        bkg_id = self._get_pha_data(id)._fix_background_id(bkg_id)
        self._background_models.get(id, {}).pop(bkg_id, None)
        self._background_sources.get(id, {}).pop(bkg_id, None)


    def _read_user_model(self, filename, *args, **kwargs):
        x = None
        y = None
        try:
            data = self.unpack_ascii(filename, *args, **kwargs)
            x = data.get_x()
            y = data.get_y()

        # we have to check for the case of a *single* column in an ascii file
        # extract the single array from the read and bypass the dataset
        except TypeError:
            y = sherpa.astro.io.backend.get_ascii_data(filename, *args,
                                               **kwargs)[1].pop()
        except:
            try:
                data = self.unpack_table(filename, *args, **kwargs)
                x = data.get_x()
                y = data.get_y()

            # we have to check for the case of a *single* column in a fits table
            # extract the single array from the read and bypass the dataset
            except TypeError:
                y = sherpa.astro.io.backend.get_table_data(filename, *args,
                                                              **kwargs)[1].pop()
            except:
                try:
                    # unpack_data doesn't include a call to try
                    # getting data from image, so try that here.
                    data = self.unpack_image(filename, *args, **kwargs)
                    #x = data.get_x()
                    y = data.get_y()
                except:
                    raise
        return (x,y)

    #@loggable()
    def load_table_model(self, modelname, filename, method=sherpa.utils.linear_interp, *args, **kwargs):
        """
        load_table_model
        
        SYNOPSIS
           Load a table model from file into a Sherpa session
           
        SYNTAX
        
        Arguments:
           modelname  - model label
        
           filename   - file from which table model data are read
        
           method     - interpolation method
                        default = linear {neville, linear}

           args       - optional arguments to pass to data reader

           kwargs     - optional keyword arguments to pass to data reader

        Returns:
           None
           
        DESCRIPTION
           Load data from a file, and put it in a new model.  This
           model can be used in fitting, just as models that containing
           functions can be used.
           
        SEE ALSO
           set_model, load_user_model, add_user_pars        
        """
        tablemodel = sherpa.models.TableModel(modelname)
        # interpolation method
        tablemodel.method = method 
        tablemodel.filename = filename

        try:
            if not sherpa.utils.is_binary_file(filename):
                raise Exception("Not a FITS file")

            read_tbl = sherpa.astro.io.backend.get_table_data
            read_hdr = sherpa.astro.io.backend.get_header_data

            blkname = 'PRIMARY'
            hdrkeys = ['HDUCLAS1', 'REDSHIFT', 'ADDMODEL']
            hdr = read_hdr(filename, blockname=blkname, hdrkeys=hdrkeys)

            addmodel = sherpa.utils.bool_cast(hdr[hdrkeys[2]])
            addredshift = sherpa.utils.bool_cast(hdr[hdrkeys[1]])

            if str(hdr[hdrkeys[0]]).upper() != 'XSPEC TABLE MODEL':
                raise Exception("Not an XSPEC table model")

            XSTableModel = sherpa.astro.xspec.XSTableModel

            blkname = 'PARAMETERS'
            colkeys = ['NAME', 'INITIAL','DELTA','BOTTOM', 'TOP',
                       'MINIMUM', 'MAXIMUM']
            hdrkeys = ['NINTPARM', 'NADDPARM']


            (colnames, cols,
             name, hdr) = read_tbl(filename, colkeys=colkeys, hdrkeys=hdrkeys,
                                       blockname=blkname, fix_type=False)
            nint = int(hdr[hdrkeys[0]])
            tablemodel = XSTableModel(filename, modelname, *cols,
                                      nint=nint, addmodel=addmodel,
                                      addredshift=addredshift)

        except Exception, e:
            #print e, type(e)
            #raise e
            x = None
            y = None
            try:
                x, y = self._read_user_model(filename, *args, **kwargs)
            except:
                # Fall back to reading plain ASCII, if no other
                # more sophisticated I/O backend loaded (such as
                # pyfits or crates) SMD 05/29/13
                data = sherpa.io.read_data(filename, ncols=2)
                x = data.x
                y = data.y
            tablemodel.load(x,y)


        self._tbl_models.append(tablemodel)
        self._add_model_component(tablemodel)

    def load_user_model(self, func, modelname, filename=None, *args, **kwargs):
        """
        load_user_model
        
        SYNOPSIS
           Load a table model from file into a Sherpa session
           
        SYNTAX
        
        Arguments:
           func       - reference to a user model function
           
           modelname  - model label
        
           filename   - file from which table model data are read
                        default = None
        
           args       - optional arguments to pass to data reader

           kwargs     - optional keyword arguments to pass to data reader

        Returns:
           None
           
        DESCRIPTION
           Take a function written by the user, and assign to a new
           user model class.  Instances of the new class can be created,
           and used as models during fits--just as ordinary Sherpa
           models can.  Optionally, data from a file can be attached to
           the model, and used in an arbitrary way by the user model
           function; but data from file is not required, the user model
           can be just a function.  After a user model is created,
           parameters need to be added with the add_user_pars function.
           
        SEE ALSO
           set_model, load_table_model, add_user_pars
        """
        usermodel = sherpa.models.UserModel(modelname)
        usermodel.calc = func
        usermodel._file = filename
        if (filename is not None):
            x, usermodel._y = self._read_user_model(filename, *args, **kwargs)
        self._add_model_component(usermodel)

    ###########################################################################
    # Fitting
    ###########################################################################


    def _add_extra_data_and_models(self, ids, datasets, models, bkg_ids={}):
        for id, d in zip(ids, datasets):
            if isinstance(d, sherpa.astro.data.DataPHA):
                bkg_models = self._background_models.get(id, {})
                bkg_srcs = self._background_sources.get(id, {})
                if d.subtracted:
                    if (bkg_models or bkg_srcs):
                        warning(('data set %r is background-subtracted; ' +
                                 'background models will be ignored') % id)
                elif not (bkg_models or bkg_srcs):
                    if d.background_ids:
                        warning(('data set %r has associated backgrounds, ' +
                                 'but they have not been subtracted, ' +
                                 'nor have background models been set') % id)
                else:
                    bkg_ids[id] = []
                    for bkg_id in d.background_ids:
                        
                        if not (bkg_id in bkg_models or bkg_id in bkg_srcs):
                            raise ModelErr('nobkg', bkg_id, id)

                        bkg = d.get_background(bkg_id)
                        datasets.append(bkg)

                        bkg_data = d
                        if len(bkg.response_ids)!=0:
                            bkg_data = bkg

                        bkg_model = bkg_models.get(bkg_id,None)
                        bkg_src = bkg_srcs.get(bkg_id,None)
                        if (bkg_model is None and bkg_src is not None):
                            resp = sherpa.astro.instrument.Response1D(bkg_data)
                            bkg_model = resp(bkg_src)
                        models.append(bkg_model)
                        bkg_ids[id].append(bkg_id)


    def _prepare_bkg_fit(self, id, otherids=()):

        # prep data ids for fitting
        ids = self._get_fit_ids(id, otherids)

        # Gather up lists of data objects and models to fit
        # to them.  Add to lists *only* if there actually is
        # a model to fit.  E.g., if data sets 1 and 2 exist,
        # but only data set 1 has a model, then "fit all" is
        # understood to mean "fit 1".  If both data sets have
        # models, then "fit all" means "fit 1 and 2 together".
        datasets = []
        models = []
        fit_to_ids = []
        for i in ids:
            
            # get PHA data and associated background models by id
            data = self._get_pha_data(i)
            bkg_models = self._background_models.get(i, {})
            bkg_sources = self._background_sources.get(i, {})

            for bi in data.background_ids:
                mod = None
                ds = self.get_bkg(i, bi)
                if bkg_models.has_key(bi) or bkg_sources.has_key(bi):
                    mod = self.get_bkg_model(i, bi)

                if mod is not None:
                    datasets.append(ds)
                    models.append(mod)

            fit_to_ids.append(i)

        # If no data sets have models assigned to them, stop now.
        if len(models) < 1:
            raise IdentifierErr("nomodels")

        return fit_to_ids, datasets, models


    def _get_bkg_fit(self, id, otherids=(), estmethod=None):

        fit_to_ids, datasets, models = self._prepare_bkg_fit(id, otherids)

        # Do not add backgrounds to backgrounds.
        #self._add_extra_data_and_models(fit_to_ids, datasets, models)

        fit_to_ids = tuple(fit_to_ids)

        f = self._get_fit_obj(datasets, models, estmethod)

        return fit_to_ids, f


    def fit(self, id=None, *otherids, **kwargs):
        """
        fit

        SYNOPSIS
           Perform fitting process using current optimization method and 
           current fit statistic.

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           otherids  - List of other Sherpa data ids

           outfile   - filename and path of parameter value output vs. number
                       of function evaluations
                       default = None

           clobber   - boolean whether to clobber outfile
                       default = False

        Returns:
           Formatted fit results output 

        DESCRIPTION
           Initiate optimization of model parameter values by id(s).

        SEE ALSO
           get_fit_results, conf, proj, covar, show_fit
        """
        kwargs['bkg_only']=False
        self._fit(id, *otherids, **kwargs)


    def fit_bkg(self, id=None, *otherids, **kwargs):
        """
        fit_bkg

        SYNOPSIS
           Perform fitting process on PHA backgrounds using current 
           optimization method and current fit statistic.

        SYNTAX

        Arguments:
           id        - Sherpa background data id
                       default = default background data id

           otherids  - List of other Sherpa background data ids

           outfile   - filename and path of parameter value output vs. number
                       of function evaluations
                       default = None

           clobber   - boolean whether to clobber outfile
                       default = False

        Returns:
           Formatted fit results output 

        DESCRIPTION
           Initiate optimization of model parameter values by background id(s).

        SEE ALSO
           get_fit_results, conf, proj, covar, show_fit
        """
        kwargs['bkg_only']=True
        self._fit(id, *otherids, **kwargs)


    def _fit(self, id=None, *otherids, **kwargs):
        ids = f = None
        fit_bkg=False

        if kwargs.has_key('bkg_only') and kwargs.pop('bkg_only'):
            fit_bkg=True

        # validate the kwds to f.fit() so user typos do not
        # result in regular fit
        #valid_keys = sherpa.utils.get_keyword_names(sherpa.fit.Fit.fit)
        valid_keys = ('outfile', 'clobber', 'filter_nan')
        for key in kwargs.keys():
            if key not in valid_keys:
                raise TypeError("unknown keyword argument: '%s'" % key)

        if fit_bkg:
            ids, f = self._get_bkg_fit(id, otherids)
        else:
            ids, f = self._get_fit(id, otherids)

	if kwargs.has_key('filter_nan') and kwargs.pop('filter_nan'):
		for i in ids:
			self.get_data(i).mask = self.get_data(i).mask & numpy.isfinite(self.get_data(i).get_x())

        res = f.fit(**kwargs)
        res.datasets = ids
        self._fit_results = res
        info(res.format())


    def _get_stat_info(self):

        ids, datasets, models = self._prepare_fit(None)

        extra_ids = {}
        self._add_extra_data_and_models(ids, datasets, models, extra_ids)

        output = []
        nids = len(ids)
        if len(datasets) > 1:
            bkg_datasets = datasets[nids:]
            bkg_models = models[nids:]
            jj = 0
            for id, d, m in izip(ids, datasets[:nids], models[:nids]):
                f = sherpa.fit.Fit(d, m, self._current_stat)

                statinfo = f.calc_stat_info()
                statinfo.name = 'Dataset %s' % (str(id))
                statinfo.ids = (id,)

                output.append(statinfo)

                bkg_ids = extra_ids.get(id, ())
                nbkg_ids = len(bkg_ids)
                idx_lo = jj*nbkg_ids
                idx_hi = idx_lo + nbkg_ids
                for bkg_id, bkg, bkg_mdl in izip(bkg_ids,
                                                 bkg_datasets[idx_lo:idx_hi],
                                                 bkg_models[idx_lo:idx_hi]):

                    bkg_f = sherpa.fit.Fit(bkg, bkg_mdl, self._current_stat)

                    statinfo = bkg_f.calc_stat_info()
                    statinfo.name = ("Background %s for Dataset %s" %
                                     (str(bkg_id), str(id)))
                    statinfo.ids = (id,)
                    statinfo.bkg_ids = (bkg_id,)

                    output.append(statinfo)

                jj += 1


        f = self._get_fit_obj(datasets, models, None)
        statinfo = f.calc_stat_info()
        if len(ids) == 1:
            statinfo.name = 'Dataset %s' % str(ids)
        else:
            statinfo.name = 'Datasets %s' % str(ids).strip("()")
        statinfo.ids = ids
        output.append(statinfo)

        return output


    ###########################################################################
    # Plotting
    ###########################################################################

    def get_model_plot(self, id=None, **kwargs):
        if isinstance(self.get_data(id), sherpa.astro.data.DataPHA):
            self._prepare_plotobj(id, self._modelhisto, **kwargs)
            return self._modelhisto

        self._prepare_plotobj(id, self._modelplot, **kwargs)
        return self._modelplot

    get_model_plot.__doc__ = sherpa.ui.utils.Session.get_model_plot.__doc__

    def get_source_plot(self, id=None, lo=None, hi=None):
        """
        get_source_plot

        SYNOPSIS
           Return a Sherpa source plot

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           lo       - low limit of plot
                      default = None
           
           hi       - high limit of plot
                      default = None

        Returns:
           Sherpa SourcePlot plot

        DESCRIPTION
           The Sherpa source plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              units        - units of grid, read-only

              xlo          - grid array, low bins

              xhi          - grid array, high bins

              flux         - unconvolved counts

              y            - convolved counts

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - None
                 linecolor      - 'red'
                 linestyle      - 1
                 linethickness  - None
                 ratioline      - N/A
                 symbolcolor    - N/A
                 symbolfill     - N/A
                 symbolstyle    - N/A
                 xaxis          - N/A
                 xerrorbars     - False
                 xlog           - False
                 yerrorbars     - False
                 ylog           - False                 

           Functions:

              prepare()
                 calculate the source model and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_source, plot_bkg, plot_arf, get_bkg_plot, get_arf_plot
        """
        ## srcplot obj is possibly reinstantiated depending on data type
        if isinstance(self.get_data(id), sherpa.astro.data.DataPHA):
            return self._prepare_plotobj(id, self._astrosourceplot, lo=lo, hi=hi)
        return self._prepare_plotobj(id, self._sourceplot)


    def get_model_component_plot(self, id, model=None):
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        if isinstance(self.get_data(id), sherpa.astro.data.DataPHA):
            self._prepare_plotobj(id, self._astrocompmdlplot, model=model)
            return self._astrocompmdlplot

        self._prepare_plotobj(id, self._compmdlplot, model=model)
        return self._compmdlplot

    get_model_component_plot.__doc__ = sherpa.ui.utils.Session.get_model_component_plot.__doc__


    def get_source_component_plot(self, id, model=None):
        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        if isinstance(self.get_data(id), sherpa.astro.data.DataPHA):
            self._prepare_plotobj(id, self._astrocompsrcplot, model=model)
            return self._astrocompsrcplot
        elif isinstance(model, sherpa.models.TemplateModel):
            self._prepare_plotobj(id, self._comptmplsrcplot, model=model)
            return self._comptmplsrcplot
        self._prepare_plotobj(id, self._compsrcplot, model=model)
        return self._compsrcplot

    get_source_component_plot.__doc__ = sherpa.ui.utils.Session.get_source_component_plot.__doc__


    def get_order_plot(self, id=None, orders=None):
        """
        get_order_plot

        SYNOPSIS
           Return a Sherpa order plot

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           orders   - array of orders
                      default = all orders

        Returns:
           Sherpa OrderPlot plot

        DESCRIPTION

        SEE ALSO
           plot_order, plot_bkg, plot_arf, get_bkg_plot, get_arf_plot
        """
        self._prepare_plotobj(id, self._orderplot, orders=orders)
        return self._orderplot


    def get_arf_plot(self, id=None, resp_id=None):
        """
        get_arf_plot

        SYNOPSIS
           Return a Sherpa ancillary response plot

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           resp_id  - response id, if multiple response exist
                      default  = default response id

        Returns:
           Sherpa ARFPlot plot

        DESCRIPTION
           The Sherpa ARF plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              units        - units of grid, read-only

              xlo          - grid array, low bins

              xhi          - grid array, high bins

              flux         - unconvolved counts

              y            - convolved counts

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - None
                 linecolor      - 'red'
                 linestyle      - 1
                 linethickness  - None
                 ratioline      - N/A
                 symbolcolor    - N/A
                 symbolfill     - N/A
                 symbolstyle    - N/A
                 xaxis          - N/A
                 xerrorbars     - False
                 xlog           - False
                 yerrorbars     - False
                 ylog           - False                 

           Functions:

              prepare()
                 populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_arf, plot_source, plot_bkg, get_source_plot, get_bkg_plot
        """
        self._prepare_plotobj(id, self._arfplot, resp_id)
        return self._arfplot

    def get_bkg_fit_plot(self, id=None, bkg_id=None):
        """
        get_bkg_fit_plot

        SYNOPSIS
           Return a Sherpa background fit plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                       default = None

        Returns:
           Sherpa BkgFitPlot plot

        DESCRIPTION
           The Sherpa background fit plot object holds a reference to a
           background plot and background model plot instance.

           Attributes:
              bkgdataplot

              bkgmodelplot

        SEE ALSO
           plot_bkg, plot_bkg_model, plot_bkg_fit, plot_bkg_source
        """
        self._prepare_plotobj(id, self._bkgfitplot, bkg_id=bkg_id)
        return self._bkgfitplot


    def get_bkg_model_plot(self, id=None, bkg_id=None):
        """
        get_bkg_model_plot

        SYNOPSIS
           Return a Sherpa background convolved model plot

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

        Returns:
           Sherpa BkgModelPlot object

        DESCRIPTION
           The Sherpa background model plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

           Functions:

              prepare()
                 calculate the model and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           plot_bkg_model, plot_bkg, plot_bkg_fit, plot_bkg_source
        """
        self._prepare_plotobj(id, self._bkgmodelhisto, bkg_id=bkg_id)
        return self._bkgmodelhisto


    def get_bkg_plot(self, id=None, bkg_id=None):
        """
        get_bkg_plot

        SYNOPSIS
           Return a Sherpa background data plot

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

        Returns:
           Sherpa BkgDataPlot object

        DESCRIPTION
           The Sherpa data plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - 'capped'
                 linecolor      - None
                 linestyle      - 0
                 linethickness  - None
                 ratioline      - N/A
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - 2
                 symbolstyle    - 4
                 xaxis          - N/A
                 xerrorbars     - False
                 xlog           - False
                 yerrorbars     - True
                 ylog           - False

           Functions:

              prepare()
                 populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_data, plot_arf, plot_source, plot_bkg, get_data_plot
        """
        self._prepare_plotobj(id, self._bkgdataplot, bkg_id=bkg_id)
        return self._bkgdataplot


    def get_bkg_source_plot(self, id=None, lo=None, hi=None, bkg_id=None):
        """
        get_bkg_source_plot

        SYNOPSIS
           Return a Sherpa background source plot

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           lo       - low limit of plot
                      default = None
           
           hi       - high limit of plot
                      default = None

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

        Returns:
           Sherpa BkgSourcePlot plot

        DESCRIPTION
           The Sherpa source plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              units        - units of grid, read-only

              xlo          - grid array, low bins

              xhi          - grid array, high bins

              flux         - unconvolved counts

              y            - convolved counts

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - None
                 linecolor      - 'red'
                 linestyle      - 1
                 linethickness  - None
                 ratioline      - N/A
                 symbolcolor    - N/A
                 symbolfill     - N/A
                 symbolstyle    - N/A
                 xaxis          - N/A
                 xerrorbars     - False
                 xlog           - False
                 yerrorbars     - False
                 ylog           - False                 

           Functions:

              prepare()
                 calculate the source model and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_bkg, plot_bkg_model, plot_bkg_fit, plot_bkg_source
        """
        self._prepare_plotobj(id, self._bkgsourceplot, bkg_id=bkg_id,
                              lo=lo, hi=hi)
        return self._bkgsourceplot


    def get_bkg_resid_plot(self, id=None, bkg_id=None):
        """
        get_bkg_resid_plot

        SYNOPSIS
           Return a Sherpa background residuals plot

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

        Returns:
           Sherpa BkgResidPlot plot

        DESCRIPTION
           The Sherpa background resid plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - 'line'
                 errthickness   - None
                 linecolor      - None
                 linestyle      - 0
                 linethickness  - None
                 ratioline      - False
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - 3
                 symbolstyle    - 2
                 xaxis          - True
                 xerrorbars     - True
                 xlog           - False
                 yerrorbars     - True
                 ylog           - False

           Functions:

              prepare()
                 calculate the residuals and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           plot_bkg, plot_bkg_model, plot_bkg_fit, plot_bkg_source,
           plot_bkg_ratio, plot_bkg_resid, plot_bkg_delchi
        """
        self._prepare_plotobj(id, self._bkgresidplot, bkg_id=bkg_id)
        return self._bkgresidplot


    def get_bkg_ratio_plot(self, id=None, bkg_id=None):
        """
        get_bkg_ratio_plot

        SYNOPSIS
           Return a Sherpa background ratio plot

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

        Returns:
           Sherpa BkgRatioPlot plot

        DESCRIPTION
           The Sherpa ratio plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - 'line'
                 errthickness   - None
                 linecolor      - None
                 linestyle      - 0
                 linethickness  - None
                 ratioline      - False
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - 3
                 symbolstyle    - 2
                 xaxis          - True
                 xerrorbars     - True
                 xlog           - False
                 yerrorbars     - True
                 ylog           - False

           Functions:

              prepare()
                 calculate the ratio and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           plot_bkg, plot_bkg_model, plot_bkg_fit, plot_bkg_source,
           plot_bkg_ratio, plot_bkg_resid, plot_bkg_delchi
        """
        self._prepare_plotobj(id, self._bkgratioplot, bkg_id=bkg_id)
        return self._bkgratioplot


    def get_bkg_delchi_plot(self, id=None, bkg_id=None):
        """
        get_bkg_delchi_plot

        SYNOPSIS
           Return a Sherpa background delchi plot

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

        Returns:
           Sherpa BkgDelchiPlot plot

        DESCRIPTION
           The Sherpa delchi plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              x            - independent variable array

              y            - dependent variable array

              yerr         - dependent variable uncertainties array

              xerr         - bin size array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - 'line'
                 errthickness   - None
                 linecolor      - None
                 linestyle      - 0
                 linethickness  - None
                 ratioline      - False
                 symbolcolor    - None
                 symbolfill     - True
                 symbolsize     - 3
                 symbolstyle    - 2
                 xaxis          - True
                 xerrorbars     - True
                 xlog           - False
                 yerrorbars     - True
                 ylog           - False

           Functions:

              prepare()
                 calculate the delta chi and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           plot_bkg, plot_bkg_model, plot_bkg_fit, plot_bkg_source,
           plot_bkg_ratio, plot_bkg_resid, plot_bkg_delchi
        """
        self._prepare_plotobj(id, self._bkgdelchiplot, bkg_id=bkg_id)
        return self._bkgdelchiplot


    def get_bkg_chisqr_plot(self, id=None, bkg_id=None):
        """
        get_bkg_chisqr_plot

        SYNOPSIS
           Return a Sherpa background chisqr plot

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

        Returns:
           Sherpa BkgChisqrPlot plot

        DESCRIPTION
           The Sherpa chisqr plot object holds references to various
           plot preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              units        - units of grid, read-only

              x            - x array

              y            - chisqr values array

              plot_prefs   - dictionary of plotting preferences

                 errcolor       - None
                 errstyle       - None
                 linecolor      - 'red'
                 linestyle      - 1
                 linethickness  - None
                 ratioline      - N/A
                 symbolcolor    - N/A
                 symbolfill     - N/A
                 symbolstyle    - N/A
                 xaxis          - N/A
                 xerrorbars     - False
                 xlog           - False
                 yerrorbars     - False
                 ylog           - False                 

           Functions:

              prepare()
                 calculate the chisqr and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization

        SEE ALSO
           plot_bkg, plot_bkg_model, plot_bkg_fit, plot_bkg_source,
           plot_bkg_ratio, plot_bkg_resid, plot_bkg_delchi
        """
        self._prepare_plotobj(id, self._bkgchisqrplot, bkg_id=bkg_id)
        return self._bkgchisqrplot


    def _prepare_energy_flux_plot(self, plot, lo, hi, id, num, bins, correlated, numcores, bkg_id):
        dist = self.sample_energy_flux(lo, hi, id, num, None, correlated, numcores, bkg_id)
        plot.prepare(dist, bins)
        return plot


    def _prepare_photon_flux_plot(self, plot, lo, hi, id, num, bins, correlated, numcores, bkg_id):
        dist = self.sample_photon_flux(lo, hi, id, num, correlated, numcores, bkg_id)
        plot.prepare(dist, bins)
        return plot


    def get_energy_flux_hist(self, lo=None, hi=None, id=None, num=7500, bins=75,
                             correlated=False, numcores=None, bkg_id=None, **kwargs):
        """
        get_energy_flux_hist

        SYNOPSIS
           Return a Sherpa energy flux histogram

        SYNTAX

        Arguments:
           lo          - lower energy bound
                         default = None

           hi          - upper energy bound
                         default = None

           id          - data id
                         default = default data id

           num         - Number of simulations
                         default = 7500

           bins        - Number of bins in the histogram
                         default = 75

           correlated  - Use a multi-variate distribution to sample parameter values
                         default = False

           numcores    - specify the number of cores for parallel processing.
                         All available cores are used by default.
                         default = None

           bkg_id      - Sherpa background id
                         default = default bkg id

           recalc      - Recompute before sending data arrays to visualizer
                         default = True

        Returns:
           Sherpa energy FluxHistogram object

        DESCRIPTION
           The Sherpa FluxHistogram object holds references to various
           histogram preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              units        - units of grid, read-only

              xlo          - grid array, low bins

              xhi          - grid array, high bins

              y            - flux distribution

              histo_prefs  - dictionary of plotting preferences

                 errcolor       - N/A
                 errstyle       - N/A
                 errthickness   - N/A
                 fillcolor      - None
                 fillopacity    - None
                 fillstyle      - None
                 linestyle      - 1
                 linecolor      - 'red'
                 linethickness  - None
                 symbolangle    - N/A
                 symbolcolor    - N/A
                 symbolfill     - N/A
                 symbolsize     - N/A
                 symbolstyle    - N/A
                 xlog           - False
                 yerrorbars     - False
                 ylog           - False

           Functions:

              prepare()
                 calculate the source model and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_energy_flux, get_photon_flux_plot, plot_photon_flux,
           sample_energy_flux, sample_photon_flux
        """
        if sherpa.utils.bool_cast(kwargs.pop('recalc',True)):
            self._prepare_energy_flux_plot(self._energyfluxplot, lo, hi, id, num,
                                           bins, correlated, numcores, bkg_id)
        return self._energyfluxplot


    def get_photon_flux_hist(self, lo=None, hi=None, id=None, num=7500, bins=75,
                             correlated=False, numcores=None, bkg_id=None, **kwargs):
        """
        get_photon_flux_hist

        SYNOPSIS
           Return a Sherpa photon flux histogram

        SYNTAX

        Arguments:
           lo          - lower energy bound
                         default = None

           hi          - upper energy bound
                         default = None

           id          - data id
                         default = default data id

           num         - Number of simulations
                         default = 7500

           bins        - Number of bins in the histogram
                         default = 75

           correlated  - Use a multi-variate distribution to sample parameter values
                         default = False

           numcores    - specify the number of cores for parallel processing.
                         All available cores are used by default.
                         default = None

           bkg_id      - Sherpa background id
                         default = default bkg id

           recalc      - Recompute before sending data arrays to visualizer
                         default = True

        Returns:
           Sherpa photon FluxHistogram object

        DESCRIPTION
           The Sherpa FluxHistogram object holds references to various
           histogram preferences and data arrays.

           Attributes:
              title        - title of plot, read-only

              xlabel       - x axis label, read-only

              ylabel       - y axis label, read-only

              units        - units of grid, read-only

              xlo          - grid array, low bins

              xhi          - grid array, high bins

              y            - flux distribution

              histo_prefs  - dictionary of plotting preferences

                 errcolor       - N/A
                 errstyle       - N/A
                 errthickness   - N/A
                 fillcolor      - None
                 fillopacity    - None
                 fillstyle      - None
                 linestyle      - 1
                 linecolor      - 'red'
                 linethickness  - None
                 symbolangle    - N/A
                 symbolcolor    - N/A
                 symbolfill     - N/A
                 symbolsize     - N/A
                 symbolstyle    - N/A
                 xlog           - False
                 yerrorbars     - False
                 ylog           - False

           Functions:

              prepare()
                 calculate the source model and populate the data arrays

              plot( overplot=False, clearwindow=True )
                 send data arrays to plotter for visualization


        SEE ALSO
           plot_photon_flux, get_energy_flux_plot, plot_energy_flux,
           sample_energy_flux, sample_photon_flux
        """
        if sherpa.utils.bool_cast(kwargs.pop('recalc',True)):
            self._prepare_photon_flux_plot(self._photonfluxplot, lo, hi, id, num,
                                           bins, correlated, numcores, bkg_id)
        return self._photonfluxplot


    def _prepare_plotobj(self, id, plotobj, resp_id=None, bkg_id=None, lo=None,
                         hi=None, orders=None, model=None):
        if isinstance(plotobj, sherpa.astro.plot.BkgFitPlot):
            plotobj.prepare(self._prepare_plotobj(id, self._bkgdataplot,
                                                  bkg_id=bkg_id),
                            self._prepare_plotobj(id, self._bkgmodelplot,
                                                  bkg_id=bkg_id))
        elif isinstance(plotobj, sherpa.plot.FitPlot):
            plotobj.prepare(self._prepare_plotobj(id, self._dataplot),
                            self._prepare_plotobj(id, self._modelplot))
        elif isinstance(plotobj, sherpa.plot.FitContour):
            plotobj.prepare(self._prepare_plotobj(id, self._datacontour),
                            self._prepare_plotobj(id, self._modelcontour))
        else:
            if isinstance(plotobj, sherpa.astro.plot.ARFPlot):
                plotobj.prepare(self._get_pha_data(id).get_arf(resp_id),
                                self._get_pha_data(id))
            elif(isinstance(plotobj, sherpa.plot.ComponentModelPlot) or 
                 isinstance(plotobj, sherpa.plot.ComponentSourcePlot)):
                plotobj.prepare(self.get_data(id), model, self.get_stat())
            elif isinstance(plotobj, sherpa.astro.plot.BkgDataPlot):
                plotobj.prepare(self.get_bkg(id, bkg_id),
                                self.get_stat())
            elif isinstance(plotobj, (sherpa.astro.plot.BkgModelPlot,
                                      sherpa.astro.plot.BkgRatioPlot,
                                      sherpa.astro.plot.BkgResidPlot,
                                      sherpa.astro.plot.BkgDelchiPlot,
                                      sherpa.astro.plot.BkgChisqrPlot,
                                      sherpa.astro.plot.BkgModelHistogram)):
                plotobj.prepare(self.get_bkg(id, bkg_id),
                                self.get_bkg_model(id, bkg_id),
                                self.get_stat())
            elif isinstance(plotobj, sherpa.astro.plot.BkgSourcePlot):
                plotobj.prepare(self.get_bkg(id, bkg_id),
                                self.get_bkg_source(id, bkg_id), lo, hi)
            elif isinstance(plotobj, sherpa.astro.plot.SourcePlot):
                data = self.get_data(id)
                src = self.get_source(id)
                plotobj.prepare(data, src, lo, hi)
            elif isinstance(plotobj, sherpa.plot.SourcePlot):
                data = self.get_data(id)
                src = self.get_source(id)
                plotobj.prepare(data, src)
            elif (isinstance(plotobj, sherpa.plot.PSFPlot) or
                  isinstance(plotobj, sherpa.plot.PSFContour) or
                  isinstance(plotobj, sherpa.plot.PSFKernelPlot) or 
                  isinstance(plotobj, sherpa.plot.PSFKernelContour)):
                plotobj.prepare(self.get_psf(id), self.get_data(id))
            elif(isinstance(plotobj, sherpa.plot.DataPlot) or
                 isinstance(plotobj, sherpa.plot.DataContour)):
                plotobj.prepare(self.get_data(id), self.get_stat())
            elif isinstance(plotobj, sherpa.astro.plot.OrderPlot):
                plotobj.prepare(self._get_pha_data(id),
                                self.get_model(id), orders )
            else:
                # Using _get_fit becomes very complicated using simulfit
                # models and datasets
                #
                #ids, f = self._get_fit(id)
                plotobj.prepare(self.get_data(id), self.get_model(id),
                                self.get_stat())

        return plotobj


    def _set_plot_item(self, plottype, item, value):
        keys = self._plot_types.keys()[:]

        if plottype.strip().lower() != "all":
            if plottype not in keys:
                raise sherpa.utils.err.PlotErr('wrongtype', plottype, str(keys))
            keys = [plottype]

        for key in keys:
            plots = self._plot_types[key]

            # the astro package complicates plotting by using a regular and
            # astro version of model plots, source plots, and
            # component plots.  One is for PHA, the other is everything else.

            # To avoid confusion for the user, when 'model' is passed.  Change
            # both the regular and astro model plot type.  Astro versions are
            # prefixed with 'astro' in the _plot_types key.
            if key in ["model", "source", "compsource", "compmodel"]:
                plots = [plots, self._plot_types["astro"+key]]
            else:
                plots = [plots]

            for plot in plots:
                if sherpa.ui.utils._is_subclass(plot.__class__,
                                                  sherpa.plot.Histogram):
                    plot.histo_prefs[item] = value
                elif sherpa.ui.utils._is_subclass(plot.__class__,
                                                sherpa.plot.Plot):
                    plot.plot_prefs[item] = value


    def plot_model(self, id=None, **kwargs):
        if isinstance(self.get_data(id), sherpa.astro.data.DataPHA):
            self._plot(id, self._modelhisto, **kwargs)
        else:
            self._plot(id, self._modelplot, **kwargs)

    plot_model.__doc__ = sherpa.ui.utils.Session.plot_model.__doc__


    ### Ahelp ingest: 2015-04-29 DJB
    def plot_arf(self, id=None, resp_id=None, **kwargs):
        """Plot the ARF associated with a data set.

        Display the effective area curve from the ARF
        component of a PHA data set.

        Arguments
        ---------
        id : int or str, optional
           The data set with an ARF. If not given then the default
           identifier is used, as returned by `get_default_id`.
        resp_id : int or str, optional
           Which ARF to use in the case that multiple ARFs are
           associated with a data set. The default is `None`,
           which means the first one.
        replot : bool, optional
           Set to `True` to use the values calculated by the last
           call to `plot_data`. The default is `False`.
        overplot : bool, optional
           If `True` then add the data to an exsiting plot, otherwise
           create a new plot. The default is `False`.

        Raises
        ------
        sherpa.utils.err.ArgumentErr
           If the data set does not contain PHA data.

        See Also
        --------
        get_default_id : Return the default data set identifier.
        plot : Create one or more plot types.

        Examples
        --------

        Plot the ARF for the default data set:

        >>> plot_arf()

        Plot the ARF from data set 1 and overplot
        the ARF from data set 2:

        >>> plot_arf(1)
        >>> plot_arf(2, overplot=True)

        Plot the ARFs labelled "arf1" and "arf2" for the
        "src" data set:

        >>> plot_arf("src", "arf1")
        >>> plot_arf("src", "arf2", overplot=True)

        """
        arf = self._get_pha_data(id).get_arf(resp_id)
        if arf is None:
            raise DataErr('noarf', self._fix_id(id))
        self._plot(id, self._arfplot, resp_id, **kwargs)


    def plot_source_component(self, id, model=None, **kwargs):

        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        plotobj = self._compsrcplot
        if isinstance(self.get_data(id), sherpa.astro.data.DataPHA):
            plotobj = self._astrocompsrcplot
        elif isinstance(model, sherpa.models.TemplateModel):
            plotobj = self._comptmplsrcplot

        self._plot(id, plotobj, None, None, None, None, None, model, **kwargs)

    plot_source_component.__doc__ = sherpa.ui.utils.Session.plot_source_component.__doc__

    def plot_model_component(self, id, model=None, **kwargs):

        if model is None:
            id, model = model, id
        self._check_model(model)
        if isinstance(model, basestring):
            model = self._eval_model_expression(model)

        is_source = self._get_model_status(id)[1]
        model = self._add_convolution_models(id, self.get_data(id),
                                             model, is_source)

        plotobj = self._compmdlplot
        if isinstance(self.get_data(id), sherpa.astro.data.DataPHA):
            plotobj = self._astrocompmdlplot

        self._plot(id, plotobj, None, None, None, None, None, model, **kwargs)


    plot_model_component.__doc__ = sherpa.ui.utils.Session.plot_model_component.__doc__

    def plot_source(self, id=None, lo=None, hi=None, **kwargs):
        """
        plot_source

        SYNOPSIS
           Plot unconvolved source model

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           lo       - low limit of plot
                      default = None

           hi       - high limit of plot
                      default = None

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the unconvolved source model in a 1D plot by
           data id.

        SEE ALSO
           plot_model, plot_data, get_source_plot, plot_arf, plot_bkg,
           plot_bkg_source
        """
        if isinstance(self.get_data(id), sherpa.astro.data.DataPHA):
            self._plot(id, self._astrosourceplot, None, None, lo, hi, **kwargs)
        else:
            self._plot(id, self._sourceplot, **kwargs)


    def plot_order(self, id=None, orders=None, **kwargs):
        """
        plot_order

        SYNOPSIS
           Plot convolved source model by multiple response order

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           orders   - list of plot orders
                      default = None

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the convolved source model in a 1D plot by
           data id and multiple response order.

        SEE ALSO
           plot_model, plot_data, get_source_plot, plot_arf, plot_bkg,
           plot_bkg_source
        """
        self._plot(id, self._orderplot, None, None, None, None,orders, **kwargs)


    def plot_bkg(self, id=None, bkg_id=None, **kwargs):
        """
        plot_bkg

        SYNOPSIS
           Plot background counts

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the background counts in a 1D plot by
           data id or background id.

        SEE ALSO
           plot_bkg_model, plot_bkg_fit, plot_bkg_source
        """
        bkg = self.get_bkg(id, bkg_id)
        self._plot(id, self._bkgdataplot, None, bkg_id, **kwargs)


    def plot_bkg_model(self, id=None, bkg_id=None, **kwargs):
        """
        plot_bkg_model

        SYNOPSIS
           Plot background convolved model

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the background convolved model in a 1D plot by
           data id and background id.

        SEE ALSO
           plot_bkg, plot_bkg_fit, plot_bkg_source
        """
        bkg = self.get_bkg(id, bkg_id)
        self._plot(id, self._bkgmodelhisto, None, bkg_id, **kwargs)


    def plot_bkg_resid(self, id=None, bkg_id=None, **kwargs):
        """
        plot_bkg_resid

        SYNOPSIS
           Plot background residuals

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the background residuals (measured background
           counts minus predicted background counts) in a 1D plot by data id and
           background id.

        SEE ALSO
           plot_bkg, plot_bkg_fit, plot_bkg_source
        """
        bkg = self.get_bkg(id, bkg_id)
        self._plot(id, self._bkgresidplot, None, bkg_id, **kwargs)


    def plot_bkg_ratio(self, id=None, bkg_id=None, **kwargs):
        """
        plot_bkg_ratio

        SYNOPSIS
           Plot background ratio

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the background ratio (background measured counts divided
           by background predicted counts) in a 1D plot by data id and
           background id.

        SEE ALSO
           plot_bkg, plot_bkg_fit, plot_bkg_source
        """
        bkg = self.get_bkg(id, bkg_id)
        self._plot(id, self._bkgratioplot, None, bkg_id, **kwargs)


    def plot_bkg_delchi(self, id=None, bkg_id=None, **kwargs):
        """
        plot_bkg_delchi

        SYNOPSIS
           Plot background delta chi

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the background delchi chi (residuals divided by
           background uncertainties) in a 1D plot by data id and background id.

        SEE ALSO
           plot_bkg, plot_bkg_fit, plot_bkg_source
        """
        bkg = self.get_bkg(id, bkg_id)
        self._plot(id, self._bkgdelchiplot, None, bkg_id, **kwargs)


    def plot_bkg_chisqr(self, id=None, bkg_id=None, **kwargs):
        """
        plot_bkg_chisqr

        SYNOPSIS
           Plot background chi squared contributions

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the background chi squared contributions in a 1D plot by
           data id and background id.

        SEE ALSO
           plot_bkg, plot_bkg_fit, plot_bkg_source
        """
        bkg = self.get_bkg(id, bkg_id)
        self._plot(id, self._bkgchisqrplot, None, bkg_id, **kwargs)


    def plot_bkg_fit(self, id=None, bkg_id=None, **kwargs):
        """
        plot_bkg_fit

        SYNOPSIS
           Plot background counts with fitted background model

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the background fit in a 1D plot by
           data id and background id.

        SEE ALSO
           plot_bkg, plot_bkg_model, plot_bkg_source
        """
        bkg = self.get_bkg(id, bkg_id)
        self._plot(id, self._bkgfitplot, None, bkg_id, **kwargs)


    def plot_bkg_source(self, id=None, lo=None, hi=None, bkg_id=None, **kwargs):
        """
        plot_bkg_source

        SYNOPSIS
           Plot the unconvolved background model

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           lo       - low limit of plot
                      default = None
   
           hi       - high limit of plot
                      default = None

           bkg_id   - bkg id, if multiple bkgs exist
                      default = None

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the unconvolved source model in a 1D plot by
           data id and bkg_id.

        SEE ALSO
           plot_bkg, plot_bkg_model, plot_bkg_fit
        """
        bkg = self.get_bkg(id, bkg_id)
        self._plot(id, self._bkgsourceplot, None, bkg_id, lo, hi, **kwargs)


    def plot_energy_flux(self, lo=None, hi=None, id=None, num=7500, bins=75,
                         correlated=False, numcores=None, bkg_id=None, **kwargs):
        """
        plot_energy_flux

        SYNOPSIS
           Send a energy flux distribution to the visualizer

        SYNTAX

        Arguments:
           lo          - lower energy bound
                         default = None

           hi          - upper energy bound
                         default = None

           id          - Sherpa data id
                         default = default data id

           num         - Number of simulations
                         default = 7500

           bins        - Number of bins in the histogram
                         default = 75

           correlated  - Use a multi-variate distribution to sample parameter values
                         default = False

           numcores    - specify the number of cores for parallel processing.
                         All available cores are used by default.
                         default = None

           bkg_id      - Sherpa background id
                         default = default bkg id

           recalc      - Recompute before sending data arrays to visualizer
                         default = True

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a energy flux histogram by Sherpa data id.

        SEE ALSO
           get_energy_flux_plot, get_photon_flux_plot, plot_photon_flux,
           sample_energy_flux, sample_photon_flux
        """
        efplot = self._energyfluxplot
        if sherpa.utils.bool_cast(kwargs.pop('recalc',True)):
            efplot = self._prepare_energy_flux_plot(efplot, lo, hi, id, num,
                                                    bins, correlated, numcores, bkg_id)
        try:
            sherpa.plot.begin()
            efplot.plot(**kwargs)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    def plot_photon_flux(self, lo=None, hi=None, id=None, num=7500, bins=75,
                         correlated=False, numcores=None, bkg_id=None, **kwargs):
        """
        plot_photon_flux

        SYNOPSIS
           Send a photon flux distribution to the visualizer

        SYNTAX

        Arguments:
           lo          - lower energy bound
                         default = None

           hi          - upper energy bound
                         default = None

           id          - Sherpa data id
                         default = default data id

           num         - Number of simulations
                         default = 7500

           bins        - Number of bins in the histogram
                         default = 75

           correlated  - Use a multi-variate distribution to sample parameter values
                         default = False

           numcores    - specify the number of cores for parallel processing.
                         All available cores are used by default.
                         default = None

           bkg_id      - Sherpa background id
                         default = default bkg id

           recalc      - Recompute before sending data arrays to visualizer
                         default = True

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize a photon flux histogram by Sherpa data id.

        SEE ALSO
           get_photon_flux_plot, get_energy_flux_plot, plot_energy_flux,
           sample_energy_flux, sample_photon_flux
        """
        pfplot = self._photonfluxplot
        if sherpa.utils.bool_cast(kwargs.pop('recalc',True)):
            pfplot = self._prepare_photon_flux_plot(pfplot, lo, hi, id, num,
                                                    bins, correlated, numcores, bkg_id)
        try:
            sherpa.plot.begin()
            pfplot.plot(**kwargs)
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    def plot_bkg_fit_resid(self, id=None, bkg_id=None, replot=False,
                           overplot=False, clearwindow=True):
        """
        plot_bkg_fit_resid

        SYNOPSIS
           Send background fit and background residuals plots to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           bkg_id      - bkg id, if multiple bkgs exist
                         default = None

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the background fit plot and background residuals plot in a
           joint plot window by Sherpa data id and bkg_id.

        SEE ALSO
           plot_bkg_resid, plot_bkg_delchi, plot_bkg_ratio, plot_bkg_chisqr,
           plot_bkg_fit, plot_bkg, plot_bkg_model, plot_bkg_fit_delchi
        """
        self._jointplot.reset()
        fp = self._bkgfitplot
        rp = self._bkgresidplot
        if not sherpa.utils.bool_cast(replot):
            fp = self._prepare_plotobj(id, fp, bkg_id=bkg_id)
            rp = self._prepare_plotobj(id, rp, bkg_id=bkg_id)
        try:
            sherpa.plot.begin()            
            self._jointplot.plottop(fp, overplot=overplot,
                                    clearwindow=clearwindow)

            oldval = rp.plot_prefs['xlog']
            if ((self._bkgdataplot.plot_prefs.has_key('xlog') and
                 self._bkgdataplot.plot_prefs['xlog']) or 
                (self._bkgmodelplot.plot_prefs.has_key('xlog') and
                 self._bkgmodelplot.plot_prefs['xlog'])):
                rp.plot_prefs['xlog']=True

            self._jointplot.plotbot(rp, overplot=overplot)

            rp.plot_prefs['xlog'] = oldval
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()


    def plot_bkg_fit_delchi(self, id=None, bkg_id=None, replot=False,
                            overplot=False, clearwindow=True):
        """
        plot_bkg_fit_delchi

        SYNOPSIS
           Send background fit and background delta chi plots to the visualizer

        SYNTAX

        Arguments:
           id          - Sherpa data id
                         default = default data id

           bkg_id      - bkg id, if multiple bkgs exist
                         default = None

           replot      - Send cached data arrays to visualizer
                         default = False

           overplot    - Plot data without clearing previous plot
                         default = False

        Returns:
           None

        DESCRIPTION
           Visualize the fit plot and delta chi plot in a joint plot
           window by Sherpa data id.

        SEE ALSO
           plot_bkg_resid, plot_bkg_delchi, plot_bkg_ratio, plot_bkg_chisqr,
           plot_bkg_fit, plot_bkg, plot_bkg_model, plot_bkg_fit_resid
        """
        self._jointplot.reset()
        fp = self._bkgfitplot
        dp = self._bkgdelchiplot
        if not sherpa.utils.bool_cast(replot):
            fp = self._prepare_plotobj(id, fp, bkg_id=bkg_id)
            dp = self._prepare_plotobj(id, dp, bkg_id=bkg_id)
        try:
            sherpa.plot.begin()           
            self._jointplot.plottop(fp, overplot=overplot,
                                    clearwindow=clearwindow)

            oldval = dp.plot_prefs['xlog']
            if ((self._bkgdataplot.plot_prefs.has_key('xlog') and
                 self._bkgdataplot.plot_prefs['xlog']) or 
                (self._bkgmodelplot.plot_prefs.has_key('xlog') and
                 self._bkgmodelplot.plot_prefs['xlog'])):
                dp.plot_prefs['xlog']=True

            self._jointplot.plotbot(dp, overplot=overplot)

            dp.plot_prefs['xlog'] = oldval
        except:
            sherpa.plot.exceptions()
            raise
        else:
            sherpa.plot.end()

    ###########################################################################
    # Analysis Functions
    ###########################################################################


    def sample_photon_flux(self, lo=None, hi=None, id=None, num=1, scales=None,
                           correlated=False, numcores=None, bkg_id=None):
        """
        sample_photon_flux

        SYNOPSIS
           Get a sample the of photon flux

        SYNTAX

        Arguments:
           lo          - lower energy bound
                         default = None

           hi          - upper energy bound
                         default = None

           id          - Sherpa data id
                         default = default data id

           num         - Number of simulations
                         default = 1

           correlated  - Use a multi-variate distribution to sample parameter values.
                         default = False

           scales      - User supplied scales for the sampling distributions.
                         If correlated is True then scales must be a symmetric
                         and postive semi-definite 2-D array_like of shape 
                         (N,N) where N is the number of free parameters,
                         otherwise scales can be a 1-D array_like, of length N.
                         default = None

           numcores    - specify the number of cores for parallel processing.
                         All available cores are used by default.
                         default = None

           bkg_id      - Sherpa background id
                         default = default bkg id

        Returns:
           array of flux value and parameter values

        DESCRIPTION
           Get a sample of the photon flux at a particular spot in parameter space.

        SEE ALSO
           get_energy_flux_plot, get_photon_flux_plot, plot_photon_flux,
           plot_energy_flux, sample_energy_flux
        """
        ids, fit = self._get_fit(id)
        data = self.get_data(id)
        src  = None
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            src  = self.get_bkg_source(id, bkg_id)
        else:
            src = self.get_source(id)
            
        correlated=sherpa.utils.bool_cast(correlated)

        return sherpa.astro.flux.sample_flux(fit, data, src,
                                             sherpa.astro.utils.calc_photon_flux,
                                             correlated, num, lo, hi, numcores,
                                             scales)


    def sample_energy_flux(self, lo=None, hi=None, id=None, num=1, scales=None,
                           correlated=False, numcores=None, bkg_id=None):
        """
        sample_energy_flux

        SYNOPSIS
           Get a sample the of energy flux

        SYNTAX

        Arguments:
           lo          - lower energy bound
                         default = None

           hi          - upper energy bound
                         default = None

           id          - Sherpa data id
                         default = default data id

           num         - Number of simulations
                         default = 1

           correlated  - Use a multi-variate distribution to sample parameter values.
                         default = False

           scales      - User supplied scales for the sampling distributions.
                         If correlated is True then scales must be a symmetric
                         and postive semi-definite 2-D array_like of shape 
                         (N,N) where N is the number of free parameters,
                         otherwise scales can be a 1-D array_like, of length N.
                         default = None

           numcores    - specify the number of cores for parallel processing.
                         All available cores are used by default.
                         default = None

           bkg_id      - Sherpa background id
                         default = default bkg id

        Returns:
           array of flux value and parameter values

        DESCRIPTION
           Get a sample of the energy flux at a particular spot in parameter space.

        SEE ALSO
           get_energy_flux_plot, get_photon_flux_plot, plot_photon_flux,
           plot_energy_flux, sample_photon_flux
        """
        ids, fit = self._get_fit(id)
        data = self.get_data(id)
        src  = None
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            src  = self.get_bkg_source(id, bkg_id)
        else:
            src = self.get_source(id)
            
        correlated=sherpa.utils.bool_cast(correlated)

        return sherpa.astro.flux.sample_flux(fit, data, src, 
                                             sherpa.astro.utils.calc_energy_flux,
                                             correlated, num, lo, hi, numcores,
                                             scales)

    def sample_flux(self, modelcomponent=None, lo=None, hi=None, id=None,
                     num=1, scales=None, correlated=False,
                     numcores=None, bkg_id=None, Xrays=True, confidence=68):
         """
         sample_flux

         SYNOPSIS
            Get a sample of the parameters with the corresponding flux and a
            flux uncertainty for a model component or a combination of model
            components.


         SYNTAX

         Arguments:
            modelcomponent - a model component or by default the model
                             expression built from the previously defined
                             models.
                             default = None

            lo             - lower energy bound
                             default = None

            hi             - upper energy bound
                             default = None

            id             - Sherpa data id
                             default = default data id

            num            - Number of realization in the sample
                             default = 1

            correlated	   - If True then include a full covariance matrix
			     to set scales for multi-variate distributions,
                             otherwise use only diagonal elements (variances).  
			     default = False

            scales	   - User supplied scales for the sampling
			     distributions.  If correlated is True then scales
                             must be a symmetric and postive semi-definite 2-D
                             array_like of shape (N,N) where N is the number of
                             free parameters, otherwise scales can be a 1-D
                             array_like, of length N.
			     default = None

            numcores       - specify the number of cores for parallel 
			     processing. All available cores are used by
                             default.
                             default = None

            bkg_id         - Sherpa background id
                             default = default bkg_id

            Xrays          - If True then calc_energy_flux used and the
                             returned flux is in units of erg/cm2/s, otherwise
                             the units are not specified and depend the data.
                             default = True

            confidence     - confidence level for the returned flux uncertainty
                             expressed as percentile
			     default = 68

         Returns:
            array of parameter values and a flux value with lower and upper
            bounds.

         DESCRIPTION
	    Get a sample of parameters with a corresponding flux and a
	    flux uncertainty for a model component or a combination of
	    model components. The model components have to be
	    previously defined and used in the fit. The samples are
	    generated from the multi-variate normal distributions with
	    the scales defined by covariance (if at the best fit) or
	    supplied (as "scales"). The flux is calculated for each
	    set of new parameters.  The returned flux value is given
	    by a sample's median with the lower and upper quantiles
	    defined by the confidence level supplied to the function.

	 EXAMPLES


         SEE ALSO
            get_energy_flux_plot, get_photon_flux_plot, plot_photon_flux,
            plot_energy_flux, sample_photon_flux, sample_energy_flux, 
	    calc_energy_flux, calc_photon_flux, plot_cdf, plot_pdf, normal_sample,
	    t_sample, get_draws
         """

         ids, fit = self._get_fit(id)
         data = self.get_data(id)
         src  = None
         if bkg_id is not None:
             data = self.get_bkg(id, bkg_id)
             src  = self.get_bkg_source(id, bkg_id)
         else:
             src = self.get_source(id)
            
         if None == modelcomponent:
             modelcomponent = src
 
         correlated=sherpa.utils.bool_cast(correlated)

         if not isinstance( modelcomponent, sherpa.models.model.Model ):
             raise ArgumentTypeErr( 'badarg', 'modelcomponent', 'a model' )

         if False == Xrays:
             samples = self.calc_energy_flux( lo=lo, hi=hi, id=id,
                                              bkg_id=bkg_id )
         else:
             # num+1 cause sample energy flux is under-reporting its result?
             samples = self.sample_energy_flux( lo=lo, hi=hi, id=id, num=num+1,
                                               scales=scales,
                                               correlated=correlated,
                                               numcores=numcores,
                                               bkg_id=bkg_id )

         return sherpa.astro.flux.calc_sample_flux( id, lo, hi, self, fit, data,
                                                    samples, modelcomponent,
                                                    confidence )

    ### Ahelp ingest: 2015-04-28 DJB
    def eqwidth(self, src, combo, id=None, lo=None, hi=None, bkg_id=None):
        """Calculate the equivalent width of an emission or absorption line.

        Parameters
        ----------
        src :
           The continuum model (this may contain multiple components).
        combo :
           The continuum plus line (absorption or emission) model.
        lo : optional
           The lower limit for the calculation (the units are set by
           `set_analysis` for the data set). The default value (`None`)
           means that the lower range of the data set is used.
        hi : optional
           The upper limit for the calculation (the units are set by
           `set_analysis` for the data set). The default value (`None`)
           means that the upper range of the data set is used.
        id : int or string, optional
           The identifier of the data set to use. The default value
           (`None`) means that the default identifier, as returned by
           `get_default_id`, is used.
        bkg_id : int or string, optional
           The identifier of the background component to use. This
           should only be set when the line to be measured is in the
           background model.

        Returns
        -------
        width : number
           The equivalent width [1]_ in the appropriate units (as given
           by `set_analysis`).

        See Also
        --------
        calc_model_sum : Sum up the fitted model over a pass band.
        calc_source_sum : Calculate the un-convolved model signal.
        get_default_id : Return the default data set identifier.
        set_model : Set the source model expression.

        References
        ----------

        .. [1] http://en.wikipedia.org/wiki/Equivalent_width

        Examples
        --------

        Set a source model (a powerlaw for the continuum and a
        gaussian for the line), fit it, and then evaluate the
        equivalent width of the line. The example assumes that
        this is a PHA data set, with an associated response,
        so that the analysis can be done in wavelength units. 

        >>> set_source(powlaw1d.cont + gauss1d.line)
        >>> set_analysis('wavelength')
        >>> fit()
        >>> eqwidth(cont, cont+line)
        2.1001988282497308

        The calculation is restricted to the range 20 to 20
        Angstroms.

        >>> eqwidth(cont, cont+line, lo=20, hi=24)
        1.9882824973082310

        The calculation is done for the background model of
        data set 2, over the range 0.5 to 2 (the units of this
        are whatever the analysis setting for this data set id).
        
        >>> set_bkg_source(2, const1d.flat + gauss1d.bline)
        >>> eqwidth(flat, flat+bline, id=2, bkg_id=1, lo=0.5, hi=2)
        0.45494599793003426

        """
        data = self.get_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        return sherpa.astro.utils.eqwidth(data, src, combo, lo, hi)


    ### Ahelp ingest: 2015-05-04 DJB
    def calc_photon_flux(self, lo=None, hi=None, id=None, bkg_id=None):
        """Integrate the source model over a pass band.

        Calculate the integral of S(E) over a pass band, where S(E) is
        the spectral model evaluated for each bin (that is, the model
        without any instrumental responses applied to it).

        Parameters
        ----------
        lo : number, optional
           The minimum limit of the band. Use `None`, the default,
           to use the low value of the data set.
        hi : number, optional
           The maximum limit of the band, which must be larger than
           `lo`. Use `None`, the default, to use the upper value of
           the data set.
        id : int or str, optional
           Use the source expression associated with this data set. If
           not given then the default identifier is used, as returned
           by `get_default_id`.
        bkg_id : int or str, optional
           If set, use the model associated with the given background
           component rather than the source model.

        Returns
        -------
        flux :
           The flux from the source model integrated over the given
           band. This represents the flux from the model without any
           instrument response (i.e. the intrinsic flux of the
           source). For X-Spec style models the units will be
           photon/cm^2/s. If `hi` is `None` but `lo` is set then the
           flux density is returned at that point: photon/cm^2/s/keV
           or photon/cm^2/s/Angstrom depending on the analysis
           setting.

        See Also
        --------
        calc_data_sum : Sum up the observed counts over a pass band.
        calc_model_sum : Sum up the fitted model over a pass band.
        calc_energy_flux : Integrate the source model over a pass band.
        calc_source_sum: Sum up the source model over a pass band.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The units of `lo` and `hi` are determined by the analysis
        setting for the data set (e.g. `get_analysis`).

        Any existing filter on the data set - e.g. as created by
        `ignore` or `notice` - is ignored by this function.

        The flux is calculated from the given source model, so if it
        includes an absorbing component then the result will represent
        the absorbed flux. The absorbing component can be removed, or
        set to absorb no photons, to get the un-absorbed flux.

        The units of the answer depend on the model components used in
        the source expression and the axis or axes of the data set.
        It is unlikely to give sensible results for 2D data sets.

        Examples
        --------

        Calculate the photon flux over the ranges 0.5 to 2 and 0.5 to
        7 keV, and compared to the energy fluxes for the same bands:

        >>> set_analysis('energy')
        >>> calc_photon_flux(0.5, 2)
        0.35190275
        >>> calc_photon_flux(0.5, 7)
        0.49050927
        >>> calc_energy_flux(0.5, 2)
        5.7224906878061796e-10
        >>> calc_energy_flux(0.5, 7)
        1.3758131915063825e-09

        Calculate the photon flux density at 0.5 keV for the source
        "core":

        >>> calc_photon_flux(0.5, id="core")
        0.64978176

        """
        
        data = self.get_data(id)
        model= None

        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            model= self.get_bkg_source(id, bkg_id)
        else:
            model = self.get_source(id)
            
        return sherpa.astro.utils.calc_photon_flux(data, model, lo, hi)
    
    ### Ahelp ingest: 2015-05-04 DJB
    def calc_energy_flux(self, lo=None, hi=None, id=None, bkg_id=None):
        """Integrate the source model over a pass band.

        Calculate the integral of E * S(E) over a pass band, where E
        is the energy of the bin and S(E) the spectral model evaluated
        for that bin (that is, the model without any instrumental
        responses applied to it).

        Parameters
        ----------
        lo : number, optional
           The minimum limit of the band. Use `None`, the default,
           to use the low value of the data set.
        hi : number, optional
           The maximum limit of the band, which must be larger than
           `lo`. Use `None`, the default, to use the upper value of
           the data set.
        id : int or str, optional
           Use the source expression associated with this data set. If
           not given then the default identifier is used, as returned
           by `get_default_id`.
        bkg_id : int or str, optional
           If set, use the model associated with the given background
           component rather than the source model.

        Returns
        -------
        flux :
           The flux from the source model integrated over the given
           band. This represents the flux from the model without any
           instrument response (i.e. the intrinsic flux of the
           source). For X-Spec style models the units will be
           erg/cm^2/s. If `hi` is `None` but `lo` is set then the flux
           density is returned at that point: erg/cm^2/s/keV or
           erg/cm^2/s/Angstrom depending on the analysis setting.

        See Also
        --------
        calc_data_sum : Sum up the data values over a pass band.
        calc_model_sum : Sum up the fitted model over a pass band.
        calc_source_sum: Sum up the source model over a pass band.
        calc_photon_flux : Integrate the source model over a pass band.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The units of `lo` and `hi` are determined by the analysis
        setting for the data set (e.g. `get_analysis`).

        Any existing filter on the data set - e.g. as created by
        `ignore` or `notice` - is ignored by this function.

        The flux is calculated from the given source model, so if it
        includes an absorbing component then the result will represent
        the absorbed flux. The absorbing component can be removed, or
        set to absorb no photons, to get the un-absorbed flux.

        The units of the answer depend on the model components used in
        the source expression and the axis or axes of the data set.
        It is unlikely to give sensible results for 2D data sets.

        Examples
        --------

        Calculate the energy flux over the ranges 0.5 to 2 and 0.5 to
        7 keV:

        >>> set_analysis('energy')
        >>> calc_energy_flux(0.5, 2)
        5.7224906878061796e-10
        >>> calc_energy_flux(0.5, 7)
        1.3758131915063825e-09

        Calculate the energy flux density at 0.5 keV for the source
        "core":

        >>> calc_energy_flux(0.5, id="core")
        5.2573786652855304e-10

        """
        data = self.get_data(id)
        model= None

        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            model= self.get_bkg_source(id, bkg_id)
        else:
            model= self.get_source(id)
        return sherpa.astro.utils.calc_energy_flux(data, model, lo, hi)


    ### Ahelp ingest: 2015-05-05 DJB
    ### DOC-TODO: how do lo/hi limits interact with bin edges;
    ###           is it all in or partially in or ...
    def calc_data_sum(self, lo=None, hi=None, id=None, bkg_id=None):
        """Sum up the data values over a pass band.

        Parameters
        ----------
        lo : number, optional
           The minimum limit of the band. Use `None`, the default,
           to use the low value of the data set.
        hi : number, optional
           The maximum limit of the band, which must be larger than
           `lo`. Use `None`, the default, to use the upper value of
           the data set.
        id : int or str, optional
           Use the source expression associated with this data set. If
           not given then the default identifier is used, as returned
           by `get_default_id`.
        bkg_id : int or str, optional
           If set, use the model associated with the given background
           component rather than the source model.

        Returns
        -------
        dsum : number
           The sum of the data values that lie within the given
           limits.  If `hi` is `None` but `lo` is set then the data
           value of the bin containing the `lo` value are returned.
           If a background estimate has been subtracted from the data
           set then the calculation will use the background-subtracted
           values.

        See Also
        --------
        calc_data_sum2d : Sum up the data values of a 2D data set.
        calc_model_sum : Sum up the fitted model over a pass band.
        calc_energy_flux : Integrate the source model over a pass band.
        calc_photon_flux : Integrate the source model over a pass band.
        calc_source_sum: Sum up the source model over a pass band.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The units of `lo` and `hi` are determined by the analysis
        setting for the data set (e.g. `get_analysis`).

        Any existing filter on the data set - e.g. as created by
        `ignore` or `notice` - is ignored by this function.

        If a grouping scheme has been applied to the data set that it
        will be used. This can change the results, since the first and
        last bins of the selected range may extend outside the
        requested range.

        Examples
        --------

        Calculate the number of counts over the ranges 0.5 to 2 and 0.5 to
        7 keV for the default data set, first using the observed signal
        and then, for the 0.5 to 2 keV band - the background-subtraced
        estimate:

        >>> set_analysis('energy')
        >>> calc_data_sum(0.5, 2)
        745.0
        >>> calc_data_sum(0.5, 7)
        60.0
        >>> subtract()
        >>> calc_data_sum(0.5, 2)
        730.9179738207356

        Calculate the data value in the bin containing 0.5 keV for the
        source "core":

        >>> calc_data_sum(0.5, id="core")
        0.0

        """
        data = self.get_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        return sherpa.astro.utils.calc_data_sum(data, lo, hi)
            
    ### Ahelp ingest: 2015-05-05 DJB
    ### DOC-TODO: does lo!=None,hi=None make sense here,
    ###           since this is not an integration but a sum.
    ###           For now I have just not documented this capability.
    ### DOC-TODO: better comparison of calc_source_sum and calc_model_sum
    ###           needed (e.g. integration or results in PHA case?)
    def calc_model_sum(self, lo=None, hi=None, id=None, bkg_id=None):
        """Sum up the fitted model over a pass band.

        Sum up S(E) over a pass band, where S(E) is the model
        evaluated for each bin (that is, the model that is fit to the
        data, so that it includes instrumental responses).

        Parameters
        ----------
        lo : number, optional
           The minimum limit of the band. Use `None`, the default,
           to use the low value of the data set.
        hi : number, optional
           The maximum limit of the band, which must be larger than
           `lo`. Use `None`, the default, to use the upper value of
           the data set.
        id : int or str, optional
           Use the source expression associated with this data set. If
           not given then the default identifier is used, as returned
           by `get_default_id`.
        bkg_id : int or str, optional
           If set, use the model associated with the given background
           component rather than the source model.

        Returns
        -------
        signal : number
           The sum of the model values used to fit the data.

        See Also
        --------
        calc_data_sum : Sum up the observed counts over a pass band.
        calc_energy_flux : Integrate the source model over a pass band.
        calc_photon_flux : Integrate the source model over a pass band.
        calc_source_sum: Sum up the source model over a pass band.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The units of `lo` and `hi` are determined by the analysis
        setting for the data set (e.g. `get_analysis`).

        Any existing filter on the data set - e.g. as created by
        `ignore` or `notice` - is ignored by this function.

        The units of the answer depend on the model components used in
        the source expression and the axis or axes of the data set.
        It is unlikely to give sensible results for 2D data sets.

        Examples
        --------

        Sum up the model over the data range 0.5 to 2 for the default
        data set, and compared to the data over the same range:

        >>> calc_model_sum(0.5, 2)
        404.97796489631639
        >>> calc_data_sum(0.5, 2)
        745.0

        """
        data = self.get_data(id)
        model= None
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            model= self.get_bkg_model(id, bkg_id)
        else:
            model= self.get_model(id)
        return sherpa.astro.utils.calc_model_sum(data, model, lo, hi)

    ### Ahelp ingest: 2015-05-05 DJB
    def calc_data_sum2d(self, reg=None, id=None):
        """Sum up the data values of a 2D data set.

        Parameters
        ----------
        reg : str, optional
           The spatial filter to use. The default, `None`, is to
           use the whole data set.
        id : int or str, optional
           Use the source expression associated with this data set. If
           not given then the default identifier is used, as returned
           by `get_default_id`.

        Returns
        -------
        dsum : number
           The sum of the data values that lie within the given
           region.

        See Also
        --------
        calc_data_sum : Sum up the data values of a data set.
        calc_model_sum2d : Sum up the fitted model for a 2D data set.
        calc_source_sum2d: Sum up the source model for a 2D data set.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The coordinate system of the region filter is determined by
        the coordinate setting for the data set (e.g. `get_coord`).

        Any existing filter on the data set - e.g. as created by
        `ignore2d` or `notice2d` - is ignored by this function.

        Examples
        --------

        The following examples use the data in the default data set
        created with the following calls, which sets the y (data)
        values to be 0 to 11 in a 3 row by 4 column image:

        >>> ivals = np.arange(12)
        >>> (y,x) = np.mgrid[0:3, 0:4]
        >>> y = y.flatten()
        >>> x = x.flatten()
        >>> load_arrays(1, x, y, ivals, (3,4), DataIMG)

        With no argument, the full data set is used:

        >>> calc_data_sum2d()
        66
        >>> ivals.sum()
        66

        Only use pixels that fall within the given spatial filters:

        >>> calc_data_sum2d('circle(2,2,1)')
        36
        >>> calc_data_sum2d('field()-circle(2,2,1)')
        30

        """
        data = self.get_data(id)
        return sherpa.astro.utils.calc_data_sum2d(data, reg)

    ### Ahelp ingest: 2015-05-05 DJB
    ### DOC-TODO: show an example with psf
    ### DOC-TODO: this needs testing as doesn't seem to be working for me
    def calc_model_sum2d(self, reg=None, id=None):
        """Sum up the fitted model for a 2D data set.

        Parameters
        ----------
        reg : str, optional
           The spatial filter to use. The default, `None`, is to
           use the whole data set.
        id : int or str, optional
           Use the source expression associated with this data set. If
           not given then the default identifier is used, as returned
           by `get_default_id`.

        Returns
        -------
        msum : number
           The sum of the model values, as fitted to the data, that
           lie within the given region. This includes any PSF
           included by `set_psf`.

        See Also
        --------
        calc_model_sum : Sum up the fitted model over a pass band.
        calc_model_sum2d : Sum up the fitted model for a 2D data set.
        calc_source_sum2d: Sum up the source model for a 2D data set.
        set_psf : Apply a PSF model to a data set.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The coordinate system of the region filter is determined by
        the coordinate setting for the data set (e.g. `get_coord`).

        Any existing filter on the data set - e.g. as created by
        `ignore2d` or `notice2d` - is ignored by this function.

        Examples
        --------

        The following examples use the data in the default data set
        created with the following calls, which sets the y (data)
        values to be 0 to 11 in a 3 row by 4 column image:

        >>> ivals = np.arange(12)
        >>> (y,x) = np.mgrid[0:3, 0:4]
        >>> y = y.flatten()
        >>> x = x.flatten()
        >>> load_arrays(1, x, y, ivals, (3,4), DataIMG)
        >>> set_source(const2d.bgnd)
        >>> bgnd.c0 = 2

        With no argument, the full data set is used. Since the model
        evaluates to 2 per pixel, and there are 12 pixels in the
        data set, the result is 24:

        >>> calc_model_sum2d()
        24.0

        Only use pixels that fall within the given spatial filters:

        >>> calc_model_sum2d('circle(2,2,1)')
        8.0
        >>> calc_model_sum2d('field()-circle(2,2,1)')
        16.0

        """
        data = self.get_data(id)
        model= self.get_model(id)
        return sherpa.astro.utils.calc_model_sum2d(data, model, reg)

    ### Ahelp ingest: 2015-05-05 DJB
    def calc_source_sum2d(self, reg=None, id=None):
        """Sum up the fitted model for a 2D data set.

        Parameters
        ----------
        reg : str, optional
           The spatial filter to use. The default, `None`, is to
           use the whole data set.
        id : int or str, optional
           Use the source expression associated with this data set. If
           not given then the default identifier is used, as returned
           by `get_default_id`.

        Returns
        -------
        msum : number
           The sum of the model values that lie within the given
           region. This does not include any PSF included by
           `set_psf`.

        See Also
        --------
        calc_model_sum2d : Sum up the fitted model for a 2D data set.
        calc_source_sum : Sum up the model over a pass band.
        set_psf : Apply a PSF model to a data set.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The coordinate system of the region filter is determined by
        the coordinate setting for the data set (e.g. `get_coord`).

        Any existing filter on the data set - e.g. as created by
        `ignore2d` or `notice2d` - is ignored by this function.

        Examples
        --------

        The following examples use the data in the default data set
        created with the following calls, which sets the y (data)
        values to be 0 to 11 in a 3 row by 4 column image:

        >>> ivals = np.arange(12)
        >>> (y,x) = np.mgrid[0:3, 0:4]
        >>> y = y.flatten()
        >>> x = x.flatten()
        >>> load_arrays(1, x, y, ivals, (3,4), DataIMG)
        >>> set_source(const2d.bgnd)
        >>> bgnd.c0 = 2

        With no argument, the full data set is used. Since the model
        evaluates to 2 per pixel, and there are 12 pixels in the
        data set, the result is 24:

        >>> calc_source_sum2d()
        24.0

        Only use pixels that fall within the given spatial filters:

        >>> calc_source_sum2d('circle(2,2,1)')
        8.0
        >>> calc_source_sum2d('field()-circle(2,2,1)')
        16.0

        """
        data = self.get_data(id)
        src= self.get_source(id)
        return sherpa.astro.utils.calc_model_sum2d(data, src, reg)

    ### Ahelp ingest: 2015-05-05 DJB
    ### DOC-TODO: does lo!=None,hi=None make sense here,
    ###           since this is not an integration but a sum.
    ###           For now I have just not documented this capability.
    def calc_source_sum(self, lo=None, hi=None, id=None, bkg_id=None):
        """Sum up the source model over a pass band.

        Sum up S(E) over a pass band, where S(E) is the spectral model
        evaluated for each bin (that is, the model without any
        instrumental responses applied to it).

        Parameters
        ----------
        lo : number, optional
           The minimum limit of the band. Use `None`, the default,
           to use the low value of the data set.
        hi : number, optional
           The maximum limit of the band, which must be larger than
           `lo`. Use `None`, the default, to use the upper value of
           the data set.
        id : int or str, optional
           Use the source expression associated with this data set. If
           not given then the default identifier is used, as returned
           by `get_default_id`.
        bkg_id : int or str, optional
           If set, use the model associated with the given background
           component rather than the source model.

        Returns
        -------
        signal : number
           The source model summed up over the given band. This does
           *not* include the bin width when using histogram-style
           ('integrated' data spaces), such as used with X-Spec
           emission - also known as additive - models.

        See Also
        --------
        calc_data_sum : Sum up the observed counts over a pass band.
        calc_model_sum : Sum up the fitted model over a pass band.
        calc_energy_flux : Integrate the source model over a pass band.
        calc_photon_flux : Integrate the source model over a pass band.
        set_model : Set the source model expression for a data set.

        Notes
        -----
        The units of `lo` and `hi` are determined by the analysis
        setting for the data set (e.g. `get_analysis`).

        Any existing filter on the data set - e.g. as created by
        `ignore` or `notice` - is ignored by this function.

        The units of the answer depend on the model components used in
        the source expression and the axis or axes of the data set.
        It is unlikely to give sensible results for 2D data sets.

        Examples
        --------

        Sum up the model over the data range 0.5 to 2 for the default
        data set:

        >>> calc_source_sum(0.5, 2)
        139.12819041922018

        Compare the output of the `calc_source_sum` and
        `calc_photon_flux` routines. A 1099-bin data space is created,
        with a model which has a value of 1 for each bin. As the bin
        width is constant, at `0.01`, the integrated value, calculated
        by `calc_photon_flux`, is ` one hundredth the value from
        `calc_data_sum`:

        >>> dataspace1d(0.01, 11, 0.01, id="test")
        >>> set_source("test", const1d.bflat)
        >>> bflat.c0 = 1
        >>> calc_source_sum(id="test")
        1099.0
        >>> calc_photon_flux(id="test")
        10.99

        """
        data = self.get_data(id)
        model= None
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            model= self.get_bkg_source(id, bkg_id)
        else:
            model= self.get_source(id)
        return sherpa.astro.utils.calc_source_sum(data, model, lo, hi)

    ### Ahelp ingest: 2015-05-04 DJB
    ### DOC-TODO: no reason can't k-correct wavelength range,
    ###           but need to work out how to identify the units
    def calc_kcorr(self, z, obslo, obshi, restlo=None, resthi=None,
                   id=None, bkg_id=None):
        """Calculate the K correction for a model.

        The K correction ([1]_, [2]_, [3]_, [4]_) is the numeric
        factor applied to measured energy fluxes to convert values in
        an observed energy band to that they are in a rest-frame
        energy band (that is, correct for the change in spectral shape
        between the rest-frame and observed-frame bands). This is
        often used when converting a flux into a luminosity.

        Parameters
        ----------
        z : number or array, >= 0
           The redshift, or redshifts, of the source.
        obslo : number
           The minimum energy of the observed band.
        obshi : number
           The maximum energy of the observed band, which must
           be larger than `obslo`.
        restlo : number or `None`
           The minimum energy of the rest-frame band. If `None` then
           use `obslo`.
        restlo : number or `None`
           The maximum energy of the rest-frame band. It must be
           larger than `restlo`. If `None` then use `obshi`.
        id : int or str, optional
           Use the source expression associated with this data set. If
           not given then the default identifier is used, as returned
           by `get_default_id`.
        bkg_id : int or str, optional
           If set, use the model associated with the given background
           component rather than the source model.

        Returns
        -------
        kz : number or array of numbers

        See Also
        --------
        calc_energy_flux : Integrate the source model over a pass band.
        dataspace1d : Create the independent axis for a 1D data set.

        Notes
        -----
        This is only defined when the analysis is in 'energy' units.

        If the model contains a redshift parameter then it should
        be set to `0`, rather than the source redshift.

        If the source model is at zero redshift, the observed energy
        band is olo to ohi, and the rest frame band is rlo to rhi
        (which need not match the observed band), then the K
        correction at a redshift z can be calculated as::

          frest = calc_energy_flux(rlo, rhi)
          fobs  = calc_energy_flux(olo*(1+z), ohi*(1+z))
          kz    = frest / fobs

        The energy ranges used - rlo to rhi and olo*(1+z) to ohi*(1+z)
        - should be fully covered by the data grid, otherwise the flux
        calculation will be truncated at the grid boundaries, leading
        to incorrect results.

        References
        ----------

        .. [1] "The K correction", Hogg, D.W., et al.
               http://arxiv.org/abs/astro-ph/0210394

        .. [2] Appendix B of Jones et al. 1998, ApJ, vol 495,
               p. 100-114.
               http://adsabs.harvard.edu/abs/1998ApJ...495..100J

        .. [3] "K and evolutionary corrections from UV to IR",
               Poggianti, B.M., A&AS, 1997, vol 122, p. 399-407.
               http://adsabs.harvard.edu/abs/1997A%26AS..122..399P

        .. [4] "Galactic evolution and cosmology - Probing the
               cosmological deceleration parameter", Yoshii, Y. &
               Takahara, F., ApJ, 1988, vol 326, p. 1-18.
               http://adsabs.harvard.edu/abs/1988ApJ...326....1Y

        Examples
        --------

        Calculate the K correction for an X-Spec apec model, with a
        source temperature of 6 keV and abundance of 0.3 solar, for
        the energy band of 0.5 to 2 keV:

        >>> dataspace1d(0.01, 10, 0.01)
        >>> set_source(xsapec.clus)
        >>> clus.kt = 6
        >>> clus.abundanc = 0.3
        >>> calc_kcorr(0.5, 0.5, 2)
        0.82799195070436793

        Calculate the K correction for a range of redshifts (0 to 2)
        using an observed frame of 0. to 2 keV and a rest frame of 0.1
        to 10 keV (the energy grid is set to ensure that it covers the
        full energy range; that is the rest-frame band and the
        observed frame band multiplied by the smallest and largest
        (1+z) terms):

        >>> dataspace1d(0.01, 11, 0.01)
        >>> zs = np.linspace(0, 2, 21)
        >>> ks = calc_kcorr(zs, 0.5, 2, restlo=0.1, resthi=10)

        """
        data = self.get_data(id)
        model= None
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            model= self.get_bkg_source(id, bkg_id)
        else:
            model= self.get_source(id)
            
        return sherpa.astro.utils.calc_kcorr(data, model, z, obslo, obshi,
                                             restlo, resthi)

    ###########################################################################
    # Session Text Save Function
    ###########################################################################

    ### Ahelp ingest: 2015-04-27 DJB
    def save_all(self, outfile=None, clobber=False):
        """Save the information about the current session to a text file.

        This differs to the `save` command in that the output is human
        readable. Three consequences are:

         1. numeric values may not be recorded to their full precision

         2. data sets are not included in the file

         3. some settings and values may not be recorded.

        Parameters
        ----------
        outfile : str, optional
           If not given the results are displayed to the screen,
           otherwise it is taken to be the name of the file to
           write the results to.
        clobber : bool, optional
           If `outfile` is not `None`, then this flag controls
           whether an existing file can be overwritten (`True`)
           or if it raises an exception (`False`, the default
           setting).

        Raises
        ------
        sherpa.utils.err.IOErr
           If `outfile` already exists and `clobber` is `False`.

        See Also
        --------
        save : Save the current Sherpa session to a file.
        restore : Load in a Sherpa session from a file.

        Notes
        -----

        Items which are not saved include:

        - user models

        - any optional keywords to comands such as `load_data`
          or `load_pha`

        - only a subset of Sherpa commands are saved.

        Examples
        --------

        Write the current Sherpa session to the screen:

        >>> save_all()

        Save the session to the file 'fit.sherpa', overwriting
        it if it already exists:

        >>> save_all('fit.sherpa', clobber=True)

        """

        # TODO:
        #
        #    1) Finish RMF, ARF settings for backgrounds    DONE
        #    2) Add PSF models, table models, kernels etc.  DONE
        #    2a) Account for multi-response model           DONE
        #    2b) And background model (set_bkg)             DONE
        #    2c) And pileup (set_pileup_model)              DONE
        #    3) Any way to deal with user models?           SKIP
        #    4) Energy, coord settings for every data set   DONE
        #    5) Filters for each data set                   DONE
        #    6) Group flags for each data set               DONE

        # 7) Save optional keyword arguments for load_data/load_bkg
        #    8) Set model integrate flags                   DONE
        #    9) Subtract flags                              DONE

        # Send output to stdout, or file

        def _send_to_outfile(all, filename=None):
            try:
                if filename is None:
                    print all
                else:
                    outfile = file(filename, 'a')
                    print >> outfile, all
            except:
                raise

        # a helper function to print parameter attributes as we wish
        def _print_par(par):
            linkstr = ""
            if par.link is not None:
                linkstr = "\nlink(%s, %s)\n" % (par.fullname, par.link.fullname)

            unitstr = ""
            if (type(par.units) == str):
                unitstr = "\"%s\"" % par.units
            
            return ((('%s.default_val = %s\n' +
                      '%s.default_min = %s\n' +
                      '%s.default_max = %s\n' +
                      '%s.val     = %s\n' +
                      '%s.min     = %s\n' +
                      '%s.max     = %s\n' +
                      '%s.units   = %s\n' +
                      '%s.frozen  = %s\n' ) %
                     (par.fullname, repr(par.default_val),
                      par.fullname, repr(par.default_min),
                      par.fullname, repr(par.default_max),
                      par.fullname, repr(par.val),
                      par.fullname, repr(par.min),
                      par.fullname, repr(par.max),
                      par.fullname, unitstr,
                      par.fullname, par.frozen)), linkstr)
        
        # Check output file can be written to

        clobber=sherpa.utils.bool_cast(clobber)
        if (type(outfile) != str and outfile != None):
            raise ArgumentTypeErr('badarg', 'string or None')
        if (type(outfile) == str and os.path.isfile(outfile) and not clobber):
            raise IOErr('filefound', outfile)

        # Import numpy
        _send_to_outfile("import numpy", outfile)

        # Save data files

        _send_to_outfile("\n######### Load Data Sets\n", outfile)
        dids = self.list_data_ids()

        cmd_id = ""
        cmd_resp_id = ""
        cmd_bkg_id = ""
    
        for id in dids:
            # But if id is a string, then quote as a string
            # But what about the rest of any possible load_data() options;
            # how do we replicate the optional keywords that were possibly
            # used?  Store them with data object?
            if (type(id) == str):
                cmd_id = "\"%s\"" % id
            else:
                cmd_id = "%s" % id
            cmd = "load_data(%s,\"%s\")" % (cmd_id, self.get_data(id).name)
            _send_to_outfile(cmd, outfile)

            # Set physical or WCS coordinates here if applicable
            # If can't be done, just pass to next
            try:
                _send_to_outfile("\n######### Set Image Coordinates \n", outfile)
                cmd = "set_coord(%s, %s)" % (cmd_id, repr(self.get_coord(id)))
                _send_to_outfile(cmd, outfile)
            except:
                pass
        
            # PHA attributes; group data if applicable
            try:
                # Only store group flags and quality flags if they were changed
                # from flags in the file
                if (self.get_data(id)._original_groups == False):
                    if (self.get_data(id).grouping != None):
                        _send_to_outfile("\n######### Data Group Flags\n", outfile)
                        cmd = "set_grouping(%s, " % cmd_id
                        cmd = cmd + "val=numpy.array(" + repr(self.get_grouping(id).tolist()) + ", numpy." + str(self.get_grouping(id).dtype) + "))"
                        _send_to_outfile(cmd, outfile)
                    if (self.get_data(id).quality != None):
                        _send_to_outfile("\n######### Data Quality Flags\n", outfile)
                        cmd = "set_quality(%s, " % cmd_id
                        cmd = cmd + "val=numpy.array(" + repr(self.get_quality(id).tolist()) + ", numpy." + str(self.get_quality(id).dtype) + "))"
                        _send_to_outfile(cmd, outfile)
                # End check for original groups and quality flags
                if (self.get_data(id).grouped == True):
                    cmd = "if (get_data(%s).grouping != None and get_data(%s).grouped == False):" % (cmd_id, cmd_id)
                    _send_to_outfile(cmd, outfile)
                    _send_to_outfile("\t######### Group Data", outfile)
                    cmd = "\tgroup(%s)" % cmd_id
                    _send_to_outfile(cmd, outfile)                
            except:
                pass

            # Add responses and ARFs, if any
            try:
                _send_to_outfile("\n######### Data Spectral Responses\n", outfile)
                rids = self.list_response_ids(id)
                cmd_resp_id = ""

                for rid in rids:
                    if (type(rid) == str):
                        cmd_resp_id = "\"%s\"" % rid
                    else:
                        cmd_resp_id = "%s" % rid

                    try:
                        arf = self.get_arf(id, rid)
                        cmd = "load_arf(%s,\"%s\",%s)" % (cmd_id, arf.name, cmd_resp_id)
                        _send_to_outfile(cmd, outfile)
                    except:
                        pass

                    try:
                        rmf = self.get_rmf(id, rid)
                        cmd = "load_rmf(%s,\"%s\",%s)" % (cmd_id, rmf.name, cmd_resp_id)
                        _send_to_outfile(cmd, outfile)
                    except:
                        pass
            except:
                pass

            # Check if this data set has associated backgrounds
            try:
                _send_to_outfile("\n######### Load Background Data Sets\n", outfile)
                bids = self.list_bkg_ids(id)
                cmd_bkg_id = ""
                for bid in bids:
                    if (type(bid) == str):
                        cmd_bkg_id = "\"%s\"" % bid
                    else:
                        cmd_bkg_id = "%s" % bid
                    cmd = "load_bkg(%s,\"%s\", bkg_id=%s)" % (cmd_id, self.get_bkg(id,bid).name, cmd_bkg_id)
                    _send_to_outfile(cmd, outfile)

                    # Group data if applicable
                    try:
                        # Only store group flags and quality flags if they were changed
                        # from flags in the file
                        if (self.get_bkg(id, bid)._original_groups == False):
                            if (self.get_bkg(id, bid).grouping != None):
                                _send_to_outfile("\n######### Background Group Flags\n", outfile)
                                cmd = "set_grouping(%s, " % cmd_id
                                cmd = cmd + "val=numpy.array(" + repr(self.get_grouping(id).tolist()) + ", numpy." + str(self.get_grouping(id, bid).dtype) + "), bkg_id=" + cmd_bkg_id + ")"
                                _send_to_outfile(cmd, outfile)
                            if (self.get_bkg(id, bid).quality != None):
                                _send_to_outfile("\n######### Background Quality Flags\n", outfile)
                                cmd = "set_quality(%s, " % cmd_id
                                cmd = cmd + "val=numpy.array(" + repr(self.get_quality(id).tolist()) + ", numpy." + str(self.get_quality(id, bid).dtype) + "), bkg_id=" + cmd_bkg_id + ")"
                                _send_to_outfile(cmd, outfile)
                        # End check for original groups and quality flags
                        if (self.get_bkg(id, bid).grouped == True):
                            cmd = "if (get_bkg(%s,%s).grouping != None and get_bkg(%s,%s).grouped == False):" % (cmd_id, cmd_bkg_id, cmd_id, cmd_bkg_id)
                            _send_to_outfile(cmd, outfile)
                            _send_to_outfile("\t######### Group Background", outfile)
                            cmd = "\tgroup(%s,%s)" % (cmd_id, cmd_bkg_id)
                            _send_to_outfile(cmd, outfile)
                    except:
                        pass
                
                    # Load background response, ARFs if any
                    _send_to_outfile("\n######### Background Spectral Responses\n", outfile)
                    rids = self.list_response_ids(id, bid)
                    cmd_resp_id = ""
                    for rid in rids:
                        if (type(rid) == str):
                            cmd_resp_id = "\"%s\"" % rid
                        else:
                            cmd_resp_id = "%s" % rid

                        try:
                            arf = self.get_arf(id, rid, bid)
                            cmd = "load_arf(%s,\"%s\",%s,%s)" % (cmd_id, arf.name, cmd_resp_id, cmd_bkg_id)
                            _send_to_outfile(cmd, outfile)
                        except:
                            pass

                        try:
                            rmf = self.get_rmf(id, rid, bid)
                            cmd = "load_rmf(%s,\"%s\",%s,%s)" % (cmd_id, rmf.name, cmd_resp_id, cmd_bkg_id)
                            _send_to_outfile(cmd, outfile)
                        except:
                            pass
                                                      
            except:
                pass

            # Set energy units if applicable
            # If can't be done, just pass to next
            try:
                _send_to_outfile("\n######### Set Energy or Wave Units\n", outfile)
                units = self.get_data(id).units
                rate = self.get_data(id).rate
                if (rate == True):
                    rate = "\"rate\""
                else:
                    rate = "\"counts\""
                factor = self.get_data(id).plot_fac
                cmd = "set_analysis(%s, %s, %s, %s)" % (cmd_id,
                                                                repr(units),
                                                                rate,
                                                                repr(factor))
                _send_to_outfile(cmd, outfile)
            except:
                pass

            # Subtract background data if applicable
            try:
                if (self.get_data(id).subtracted == True):
                    cmd = "if (get_data(%s).subtracted == False):" % cmd_id
                    _send_to_outfile(cmd, outfile)
                    _send_to_outfile("\t######### Subtract Background Data", outfile)
                    cmd = "\tsubtract(%s)" % cmd_id
                    _send_to_outfile(cmd, outfile)
            except:
                pass

            # Set filter if applicable
            try:
                _send_to_outfile("\n######### Filter Data\n", outfile)
                if (len(self.get_data(id).get_filter()) > 0):
                    filter = self.get_data(id).get_filter()
                    if (len(self.get_data(id).get_dims()) == 1):
                        cmd = "notice_id(%s,\"%s\")" % (cmd_id, filter)
                        _send_to_outfile(cmd, outfile)
                    if (len(self.get_data(id).get_dims()) == 2):
                        cmd = "notice2d_id(%s,\"%s\")" % (cmd_id, filter)
                        _send_to_outfile(cmd, outfile)
            except:
                pass
        
        _send_to_outfile("", outfile)

        # Save statistic

        _send_to_outfile("\n######### Set Statistic\n", outfile)
        cmd = "set_stat(\"%s\")" % self.get_stat_name()
        _send_to_outfile(cmd, outfile)
        _send_to_outfile("", outfile)

        # Save fitting method

        _send_to_outfile("\n######### Set Fitting Method\n", outfile)
        cmd = "set_method(\"%s\")" % self.get_method_name()
        _send_to_outfile(cmd, outfile)
        _send_to_outfile("", outfile)

        mdict = self.get_method_opt()
        for key in mdict:
            val = mdict.get(key)
            cmd = "set_method_opt(\"%s\", %s)" % (key, val)
            _send_to_outfile(cmd, outfile)
        _send_to_outfile("", outfile)

        # Save iterative fitting method (if any)
        if (self.get_iter_method_name() != 'none'):
            _send_to_outfile("\n######### Set Iterative Fitting Method\n", outfile)
            cmd = "set_iter_method(\"%s\")" % self.get_iter_method_name()
            _send_to_outfile(cmd, outfile)
            _send_to_outfile("", outfile)

            mdict = self.get_iter_method_opt()
            for key in mdict:
                val = mdict.get(key)
                cmd = "set_iter_method_opt(\"%s\", %s)" % (key, val)
                _send_to_outfile(cmd, outfile)
            _send_to_outfile("", outfile)

        # Save all model components

        # Have to call elements in list in reverse order (item at end of
        # list was model first created), so that we get links right, if any

        # To recreate attributes, print out dictionary as ordered pairs,
        # for each parameter

        _send_to_outfile("\n######### Set Model Components and Parameters\n", outfile)
        all_model_components = self.list_model_components()
        all_model_components.reverse()

        # If there are any links between parameters, store link commands here
        # Then, *after* processing all models in the for loop below, send
        # link commands to outfile -- *all* models need to be created before
        # *any* links between parameters can be established.
        linkstr = ""
        for mod in all_model_components:

            # get actual model instance from the name we are given
            # then get model type, and name of this instance.
            mod = eval(mod)
            typename = mod.type
            modelname = mod.name.partition(".")[2]

            # Special cases:

            # account for PSF creation elsewhere (above, with load_data
            # commands);

            # for table models, use "load_table_model", to ensure we
            # add to lists of known models *and* separate list of
            # tabel models;

            # skip user models entirely, as they require importation of
            # user modules, beyond scope of this script.
        
            if (typename != "psfmodel" and typename != "tabelmodel" and
                typename != "usermodel"):
                # Normal case:  create an instance of the model.
                cmd = "eval(\"%s.%s\")" % (typename, modelname)
                _send_to_outfile(cmd, outfile)
            if (typename == "psfmodel"):
                cmd = "load_psf(\"%s\", \"%s\")" % (mod._name, mod.kernel.name)
                _send_to_outfile(cmd, outfile)
                try:
                    psfmod = self.get_psf(id)
                    cmd = "set_psf(%s, %s)" % (cmd_id, psfmod._name)
                    _send_to_outfile(cmd, outfile)
                except:
                    pass
            if (typename == "tablemodel"):
                # Create table model with load_table_model
                cmd = "load_table_model(\"%s\", \"%s\")" % (modelname , mod.filename)
                _send_to_outfile(cmd, outfile)

            if (typename == "convolutionkernel"):
                # Create general convolution kernel with load_conv
                cmd = "load_conv(\"%s\", \"%s\")" % (modelname , mod.kernel.name)
                _send_to_outfile(cmd, outfile)

            if (typename == "usermodel"):
                # Skip user models -- don't create, don't set parameters
                # Go directly to next model in the model component list.
                _send_to_outfile("WARNING: User model not saved, add any user model to save file manually\n", outfile)
                continue

            if (hasattr(mod, "integrate") == True):
                cmd = "%s.integrate = %s" % (modelname, mod.integrate)
                _send_to_outfile(cmd, outfile)
                _send_to_outfile("", outfile)

            from sherpa.models import Parameter
            for par in mod.__dict__.values():
                if (type(par) == Parameter or
                    issubclass(Parameter, type(par)) == True):
                    par_attributes, par_linkstr = _print_par(par)
                    _send_to_outfile(par_attributes, outfile)
                    linkstr = linkstr + par_linkstr
            # If the model is a PSFModel, could have special
            # attributes "size" and "center" -- if so, record them.
            if (typename == "psfmodel"):
                if (hasattr(mod,"size") == True):
                    cmd = "%s.size = %s" % (modelname, repr(mod.size))
                    _send_to_outfile(cmd, outfile)
                    _send_to_outfile("", outfile)
                if (hasattr(mod,"center") == True):
                    cmd = "%s.center = %s" % (modelname, repr(mod.center))
                    _send_to_outfile(cmd, outfile)
                    _send_to_outfile("", outfile)

            

        # If there were any links made between parameters, send those
        # link commands to outfile now; else, linkstr is just an empty string
        _send_to_outfile(linkstr, outfile)
        
        # Save all source, pileup and background models
        
        _send_to_outfile("\n######### Set Source, Pileup and Background Models\n", outfile)
        for id in dids:
            if (type(id) == str):
                cmd_id = "\"%s\"" % id
            else:
                cmd_id = "%s" % id
            
            # If a data set has a source model associated with it,
            # set that here -- try to distinguish cases where
            # source model is different from whole model.
            # If not, just pass
            try:
                the_source = None
                the_full_model = None
                try:
                    the_source = self.get_source(id)
                except:
                    the_source = None
                    pass

                try:
                    the_full_model = self.get_model(id)
                except:
                    the_full_model = None
                    pass

                if (the_source is None and
                    the_full_model is None):
                    cmd=""
                    pass
                elif (the_source is None and
                      the_full_model is not None):
                    cmd = "set_full_model(%s, %s)" % (cmd_id, the_full_model.name)
                elif (the_source is not None and
                      the_full_model is None):
                    cmd = "set_source(%s, %s)" % (cmd_id, self.the_source.name)
                elif (the_source is not None and
                    the_full_model is not None):
                    if (repr(the_source) == repr(the_full_model)):
                        cmd = "set_full_model(%s, %s)" % (cmd_id, the_full_model.name)
                    else:
                        cmd = "set_source(%s, %s)" % (cmd_id, self.the_source.name)
                else:
                    # You can't actually get here
                    cmd=""
                    pass
                _send_to_outfile(cmd, outfile)
                _send_to_outfile("", outfile)
            except:
                pass

            # If any pileup models, try to set them.  If not, just pass.
            try:
                cmd = "set_pileup_model(%s, %s)" % (cmd_id, self.get_pileup_model(id).name)
                _send_to_outfile(cmd, outfile)
            except:
                pass
        
        
            # Set background models (if any) associated with backgrounds
            # tied to this data set -- if none, then pass.  Again, try
            # to distinguish cases where background "source" model is
            # different from whole background model.
            try:
                bids = self.list_bkg_ids(id)
                cmd_bkg_id = ""
                for bid in bids:
                    if (type(bid) == str):
                        cmd_bkg_id = "\"%s\"" % bid
                    else:
                        cmd_bkg_id = "%s" % bid

                    the_bkg_source = None
                    the_bkg_full_model = None
                    try:
                        the_bkg_source = self.get_bkg_source(bid)
                    except:
                        the_bkg_source = None
                        pass

                    try:
                        the_bkg_full_model = self.get_bkg_model(id)
                    except:
                        the_bkg_full_model = None
                        pass

                    if (the_bkg_source is None and
                        the_bkg_full_model is None):
                        cmd=""
                        pass
                    elif (the_bkg_source is None and
                          the_bkg_full_model is not None):
                        cmd = "set_bkg_full_model(%s, %s, bkg_id=%s)" % (cmd_id, the_bkg_full_model.name, cmd_bkg_id)
                    elif (the_bkg_source is not None and
                          the_bkg_full_model is None):
                        cmd = "set_bkg_source(%s, %s, bkg_id=%s)" % (cmd_id, the_bkg_source.name, cmd_bkg_id)
                    elif (the_bkg_source is not None and
                          the_bkg_full_model is not None):
                        if (repr(the_bkg_source) == repr(the_bkg_full_model)):
                            cmd = "set_bkg_full_model(%s, %s, bkg_id=%s)" % (cmd_id, the_bkg_full_model.name, cmd_bkg_id)
                        else:
                            cmd = "set_bkg_source(%s, %s, bkg_id=%s)" % (cmd_id, the_bkg_source.name, cmd_bkg_id)
                    else:
                        # You can't actually get here
                        cmd=""
                        pass
                    _send_to_outfile(cmd, outfile)
                    _send_to_outfile("", outfile)

            except:
                pass

        # Save XSPEC settings if XSPEC module has been loaded.
        if (hasattr(sherpa.astro, "xspec")):
            _send_to_outfile("\n######### XSPEC Module Settings\n", outfile)
            xspec_state = sherpa.astro.xspec.get_xsstate()

            cmd = "set_xschatter(%d)" % xspec_state["chatter"]
            _send_to_outfile(cmd, outfile)
            cmd = "set_xsabund(\"%s\")" % xspec_state["abund"]
            _send_to_outfile(cmd, outfile)
            cmd = "set_xscosmo(%g, %g, %g)" % (xspec_state["cosmo"][0],
                                               xspec_state["cosmo"][1],
                                               xspec_state["cosmo"][2])
            _send_to_outfile(cmd, outfile)
            cmd = "set_xsxsect(\"%s\")" % xspec_state["xsect"]
            _send_to_outfile(cmd, outfile)
            for name in xspec_state["modelstrings"].keys():
                cmd = "set_xsxset(\"%s\", \"%s\")" % (name,
                                                      xspec_state["modelstrings"][name])
                _send_to_outfile(cmd, outfile)






    def save_session(self, outfile=None, clobber=False):

        # Send output to stdout, or file

        def _send_to_outfile(all, filename=None):
            if all is None:
                return
            try:
                if filename is None:
                    print all
                else:
                    outfile = file(filename, 'a')
                    print >> outfile, all
            except:
                raise

        # a helper function to print parameter attributes as we wish
        def _print_par(par):
            linkstr = ""
            if par.link is not None:
                linkstr = "\nlink(%s, %s)\n" % (par.fullname, par.link.fullname)

            unitstr = ""
            if (type(par.units) == str):
                unitstr = "\"%s\"" % par.units
            
            return ((('%s.default_val = %s\n' +
                      '%s.default_min = %s\n' +
                      '%s.default_max = %s\n' +
                      '%s.val     = %s\n' +
                      '%s.min     = %s\n' +
                      '%s.max     = %s\n' +
                      '%s.units   = %s\n' +
                      '%s.frozen  = %s\n' ) %
                     (par.fullname, repr(par.default_val),
                      par.fullname, repr(par.default_min),
                      par.fullname, repr(par.default_max),
                      par.fullname, repr(par.val),
                      par.fullname, repr(par.min),
                      par.fullname, repr(par.max),
                      par.fullname, unitstr,
                      par.fullname, par.frozen)), linkstr)
        
        # Check output file can be written to

        clobber=sherpa.utils.bool_cast(clobber)
        if (type(outfile) != str and outfile != None):
            raise ArgumentTypeErr('badarg', 'string or None')
        if (type(outfile) == str and os.path.isfile(outfile) and not clobber):
            raise IOErr('filefound', outfile)

        # Import numpy
        _send_to_outfile("import numpy", outfile)

        # Save data files

        _send_to_outfile("\n######### Load Data Sets\n", outfile)
        dids = self.list_data_ids()

	def get_logged_call(call_name, id=None):
		if id is not None:
			if self._calls_tracker.has_key(id) and self._calls_tracker[id].has_key(call_name):
				return self._calls_tracker[id][call_name]
		else:
			if self._calls_tracker.has_key(call_name):
				return self._calls_tracker[call_name]
        cmd_id = ""
        cmd_resp_id = ""
        cmd_bkg_id = ""
    
        for id in dids:
            # But if id is a string, then quote as a string
            # But what about the rest of any possible load_data() options;
            # how do we replicate the optional keywords that were possibly
            # used?  Store them with data object?
            if (type(id) == str):
                cmd_id = "\"%s\"" % id
            else:
                cmd_id = "%s" % id

            cmd = get_logged_call('load_data', id)
            _send_to_outfile(cmd, outfile)

            # Set physical or WCS coordinates here if applicable
            # If can't be done, just pass to next
            try:
                _send_to_outfile("\n######### Set Image Coordinates \n", outfile)
                cmd = get_logged_call('set_coord', id)
                _send_to_outfile(cmd, outfile)
            except:
                pass
        
            # PHA attributes; group data if applicable
            try:
                # Only store group flags and quality flags if they were changed
                # from flags in the file
                if (self.get_data(id)._original_groups == False):
                    if (self.get_data(id).grouping != None):
                        _send_to_outfile("\n######### Data Group Flags\n", outfile)
			cmd = get_logged_call('set_grouping')
                        cmd = "set_grouping(%s, " % cmd_id
                        cmd = cmd + "val=numpy.array(" + repr(self.get_grouping(id).tolist()) + ", numpy." + str(self.get_grouping(id).dtype) + "))"
                        _send_to_outfile(cmd, outfile)
                    if (self.get_data(id).quality != None):
                        _send_to_outfile("\n######### Data Quality Flags\n", outfile)
                        cmd = "set_quality(%s, " % cmd_id
                        cmd = cmd + "val=numpy.array(" + repr(self.get_quality(id).tolist()) + ", numpy." + str(self.get_quality(id).dtype) + "))"
                        _send_to_outfile(cmd, outfile)
                # End check for original groups and quality flags
                if (self.get_data(id).grouped == True):
                    cmd = "if (get_data(%s).grouping != None and get_data(%s).grouped == False):" % (cmd_id, cmd_id)
                    _send_to_outfile(cmd, outfile)
                    _send_to_outfile("\t######### Group Data", outfile)
                    cmd = "\tgroup(%s)" % cmd_id
                    _send_to_outfile(cmd, outfile)                
            except:
                pass

            # Add responses and ARFs, if any
            try:
                _send_to_outfile("\n######### Data Spectral Responses\n", outfile)
                rids = self.list_response_ids(id)
                cmd_resp_id = ""

                for rid in rids:
                    if (type(rid) == str):
                        cmd_resp_id = "\"%s\"" % rid
                    else:
                        cmd_resp_id = "%s" % rid

                    try:
                        arf = self.get_arf(id, rid)
                        cmd = "load_arf(%s,\"%s\",%s)" % (cmd_id, arf.name, cmd_resp_id)
                        _send_to_outfile(cmd, outfile)
                    except:
                        pass

                    try:
                        rmf = self.get_rmf(id, rid)
                        cmd = "load_rmf(%s,\"%s\",%s)" % (cmd_id, rmf.name, cmd_resp_id)
                        _send_to_outfile(cmd, outfile)
                    except:
                        pass
            except:
                pass

            # Check if this data set has associated backgrounds
            try:
                _send_to_outfile("\n######### Load Background Data Sets\n", outfile)
                bids = self.list_bkg_ids(id)
                cmd_bkg_id = ""
                for bid in bids:
                    if (type(bid) == str):
                        cmd_bkg_id = "\"%s\"" % bid
                    else:
                        cmd_bkg_id = "%s" % bid
                    cmd = "load_bkg(%s,\"%s\", bkg_id=%s)" % (cmd_id, self.get_bkg(id,bid).name, cmd_bkg_id)
                    _send_to_outfile(cmd, outfile)

                    # Group data if applicable
                    try:
                        # Only store group flags and quality flags if they were changed
                        # from flags in the file
                        if (self.get_bkg(id, bid)._original_groups == False):
                            if (self.get_bkg(id, bid).grouping != None):
                                _send_to_outfile("\n######### Background Group Flags\n", outfile)
                                cmd = "set_grouping(%s, " % cmd_id
                                cmd = cmd + "val=numpy.array(" + repr(self.get_grouping(id).tolist()) + ", numpy." + str(self.get_grouping(id, bid).dtype) + "), bkg_id=" + cmd_bkg_id + ")"
                                _send_to_outfile(cmd, outfile)
                            if (self.get_bkg(id, bid).quality != None):
                                _send_to_outfile("\n######### Background Quality Flags\n", outfile)
                                cmd = "set_quality(%s, " % cmd_id
                                cmd = cmd + "val=numpy.array(" + repr(self.get_quality(id).tolist()) + ", numpy." + str(self.get_quality(id, bid).dtype) + "), bkg_id=" + cmd_bkg_id + ")"
                                _send_to_outfile(cmd, outfile)
                        # End check for original groups and quality flags
                        if (self.get_bkg(id, bid).grouped == True):
                            cmd = "if (get_bkg(%s,%s).grouping != None and get_bkg(%s,%s).grouped == False):" % (cmd_id, cmd_bkg_id, cmd_id, cmd_bkg_id)
                            _send_to_outfile(cmd, outfile)
                            _send_to_outfile("\t######### Group Background", outfile)
                            cmd = "\tgroup(%s,%s)" % (cmd_id, cmd_bkg_id)
                            _send_to_outfile(cmd, outfile)
                    except:
                        pass
                
                    # Load background response, ARFs if any
                    _send_to_outfile("\n######### Background Spectral Responses\n", outfile)
                    rids = self.list_response_ids(id, bid)
                    cmd_resp_id = ""
                    for rid in rids:
                        if (type(rid) == str):
                            cmd_resp_id = "\"%s\"" % rid
                        else:
                            cmd_resp_id = "%s" % rid

                        try:
                            arf = self.get_arf(id, rid, bid)
                            cmd = "load_arf(%s,\"%s\",%s,%s)" % (cmd_id, arf.name, cmd_resp_id, cmd_bkg_id)
                            _send_to_outfile(cmd, outfile)
                        except:
                            pass

                        try:
                            rmf = self.get_rmf(id, rid, bid)
                            cmd = "load_rmf(%s,\"%s\",%s,%s)" % (cmd_id, rmf.name, cmd_resp_id, cmd_bkg_id)
                            _send_to_outfile(cmd, outfile)
                        except:
                            pass
                                                      
            except:
                pass

            # Set energy units if applicable
            # If can't be done, just pass to next
            try:
                _send_to_outfile("\n######### Set Energy or Wave Units\n", outfile)
                units = self.get_data(id).units
                rate = self.get_data(id).rate
                if (rate == True):
                    rate = "\"rate\""
                else:
                    rate = "\"counts\""
                factor = self.get_data(id).plot_fac
                cmd = "set_analysis(%s, %s, %s, %s)" % (cmd_id,
                                                                repr(units),
                                                                rate,
                                                                repr(factor))
                _send_to_outfile(cmd, outfile)
            except:
                pass

            # Subtract background data if applicable
            try:
                if (self.get_data(id).subtracted == True):
                    cmd = "if (get_data(%s).subtracted == False):" % cmd_id
                    _send_to_outfile(cmd, outfile)
                    _send_to_outfile("\t######### Subtract Background Data", outfile)
                    cmd = "\tsubtract(%s)" % cmd_id
                    _send_to_outfile(cmd, outfile)
            except:
                pass

            # Set filter if applicable
            try:
                _send_to_outfile("\n######### Filter Data\n", outfile)
                if (len(self.get_data(id).get_filter()) > 0):
                    filter = self.get_data(id).get_filter()
                    if (len(self.get_data(id).get_dims()) == 1):
                        cmd = "notice_id(%s,\"%s\")" % (cmd_id, filter)
                        _send_to_outfile(cmd, outfile)
                    if (len(self.get_data(id).get_dims()) == 2):
                        cmd = "notice2d_id(%s,\"%s\")" % (cmd_id, filter)
                        _send_to_outfile(cmd, outfile)
            except:
                pass
        
        _send_to_outfile("", outfile)

        # Save statistic

        _send_to_outfile("\n######### Set Statistic\n", outfile)
        cmd = "set_stat(\"%s\")" % self.get_stat_name()
        _send_to_outfile(cmd, outfile)
        _send_to_outfile("", outfile)

        # Save fitting method

        _send_to_outfile("\n######### Set Fitting Method\n", outfile)
        cmd = "set_method(\"%s\")" % self.get_method_name()
        _send_to_outfile(cmd, outfile)
        _send_to_outfile("", outfile)

        mdict = self.get_method_opt()
        for key in mdict:
            val = mdict.get(key)
            cmd = "set_method_opt(\"%s\", %s)" % (key, val)
            _send_to_outfile(cmd, outfile)
        _send_to_outfile("", outfile)

        # Save iterative fitting method (if any)
        if (self.get_iter_method_name() != 'none'):
            _send_to_outfile("\n######### Set Iterative Fitting Method\n", outfile)
            cmd = "set_iter_method(\"%s\")" % self.get_iter_method_name()
            _send_to_outfile(cmd, outfile)
            _send_to_outfile("", outfile)

            mdict = self.get_iter_method_opt()
            for key in mdict:
                val = mdict.get(key)
                cmd = "set_iter_method_opt(\"%s\", %s)" % (key, val)
                _send_to_outfile(cmd, outfile)
            _send_to_outfile("", outfile)

        # Save all model components

        # Have to call elements in list in reverse order (item at end of
        # list was model first created), so that we get links right, if any

        # To recreate attributes, print out dictionary as ordered pairs,
        # for each parameter

        _send_to_outfile("\n######### Set Model Components and Parameters\n", outfile)
        all_model_components = self.list_model_components()
        all_model_components.reverse()

        # If there are any links between parameters, store link commands here
        # Then, *after* processing all models in the for loop below, send
        # link commands to outfile -- *all* models need to be created before
        # *any* links between parameters can be established.
        linkstr = ""
        for mod in all_model_components:

            # get actual model instance from the name we are given
            # then get model type, and name of this instance.
            mod = eval(mod)
            typename = mod.type
            modelname = mod.name.partition(".")[2]

            # Special cases:

            # account for PSF creation elsewhere (above, with load_data
            # commands);

            # for table models, use "load_table_model", to ensure we
            # add to lists of known models *and* separate list of
            # tabel models;

            # skip user models entirely, as they require importation of
            # user modules, beyond scope of this script.
        
            if (typename != "psfmodel" and typename != "tabelmodel" and
                typename != "usermodel"):
                # Normal case:  create an instance of the model.
                cmd = "eval(\"%s.%s\")" % (typename, modelname)
                _send_to_outfile(cmd, outfile)
            if (typename == "psfmodel"):
                cmd = "load_psf(\"%s\", \"%s\")" % (mod._name, mod.kernel.name)
                _send_to_outfile(cmd, outfile)
                try:
                    psfmod = self.get_psf(id)
                    cmd = "set_psf(%s, %s)" % (cmd_id, psfmod._name)
                    _send_to_outfile(cmd, outfile)
                except:
                    pass
            if (typename == "tablemodel"):
                # Create table model with load_table_model
                cmd = "load_table_model(\"%s\", \"%s\")" % (modelname , mod.filename)
                _send_to_outfile(cmd, outfile)

            if (typename == "convolutionkernel"):
                # Create general convolution kernel with load_conv
                cmd = "load_conv(\"%s\", \"%s\")" % (modelname , mod.kernel.name)
                _send_to_outfile(cmd, outfile)

            if (typename == "usermodel"):
                # Skip user models -- don't create, don't set parameters
                # Go directly to next model in the model component list.
                _send_to_outfile("WARNING: User model not saved, add any user model to save file manually\n", outfile)
                continue

            if (hasattr(mod, "integrate") == True):
                cmd = "%s.integrate = %s" % (modelname, mod.integrate)
                _send_to_outfile(cmd, outfile)
                _send_to_outfile("", outfile)

            from sherpa.models import Parameter
            for par in mod.__dict__.values():
                if (type(par) == Parameter or
                    issubclass(Parameter, type(par)) == True):
                    par_attributes, par_linkstr = _print_par(par)
                    _send_to_outfile(par_attributes, outfile)
                    linkstr = linkstr + par_linkstr
            # If the model is a PSFModel, could have special
            # attributes "size" and "center" -- if so, record them.
            if (typename == "psfmodel"):
                if (hasattr(mod,"size") == True):
                    cmd = "%s.size = %s" % (modelname, repr(mod.size))
                    _send_to_outfile(cmd, outfile)
                    _send_to_outfile("", outfile)
                if (hasattr(mod,"center") == True):
                    cmd = "%s.center = %s" % (modelname, repr(mod.center))
                    _send_to_outfile(cmd, outfile)
                    _send_to_outfile("", outfile)

            

        # If there were any links made between parameters, send those
        # link commands to outfile now; else, linkstr is just an empty string
        _send_to_outfile(linkstr, outfile)
        
        # Save all source, pileup and background models
        
        _send_to_outfile("\n######### Set Source, Pileup and Background Models\n", outfile)
        for id in dids:
            if (type(id) == str):
                cmd_id = "\"%s\"" % id
            else:
                cmd_id = "%s" % id
            
            # If a data set has a source model associated with it,
            # set that here -- try to distinguish cases where
            # source model is different from whole model.
            # If not, just pass
            try:
                the_source = None
                the_full_model = None
                try:
                    the_source = self.get_source(id)
                except:
                    the_source = None
                    pass

                try:
                    the_full_model = self.get_model(id)
                except:
                    the_full_model = None
                    pass

                if (the_source is None and
                    the_full_model is None):
                    cmd=""
                    pass
                elif (the_source is None and
                      the_full_model is not None):
                    cmd = "set_full_model(%s, %s)" % (cmd_id, the_full_model.name)
                elif (the_source is not None and
                      the_full_model is None):
                    cmd = "set_source(%s, %s)" % (cmd_id, self.the_source.name)
                elif (the_source is not None and
                    the_full_model is not None):
                    if (repr(the_source) == repr(the_full_model)):
                        cmd = "set_full_model(%s, %s)" % (cmd_id, the_full_model.name)
                    else:
                        cmd = "set_source(%s, %s)" % (cmd_id, self.the_source.name)
                else:
                    # You can't actually get here
                    cmd=""
                    pass
                _send_to_outfile(cmd, outfile)
                _send_to_outfile("", outfile)
            except:
                pass

            # If any pileup models, try to set them.  If not, just pass.
            try:
                cmd = "set_pileup_model(%s, %s)" % (cmd_id, self.get_pileup_model(id).name)
                _send_to_outfile(cmd, outfile)
            except:
                pass
        
        
            # Set background models (if any) associated with backgrounds
            # tied to this data set -- if none, then pass.  Again, try
            # to distinguish cases where background "source" model is
            # different from whole background model.
            try:
                bids = self.list_bkg_ids(id)
                cmd_bkg_id = ""
                for bid in bids:
                    if (type(bid) == str):
                        cmd_bkg_id = "\"%s\"" % bid
                    else:
                        cmd_bkg_id = "%s" % bid

                    the_bkg_source = None
                    the_bkg_full_model = None
                    try:
                        the_bkg_source = self.get_bkg_source(bid)
                    except:
                        the_bkg_source = None
                        pass

                    try:
                        the_bkg_full_model = self.get_bkg_model(id)
                    except:
                        the_bkg_full_model = None
                        pass

                    if (the_bkg_source is None and
                        the_bkg_full_model is None):
                        cmd=""
                        pass
                    elif (the_bkg_source is None and
                          the_bkg_full_model is not None):
                        cmd = "set_bkg_full_model(%s, %s, bkg_id=%s)" % (cmd_id, the_bkg_full_model.name, cmd_bkg_id)
                    elif (the_bkg_source is not None and
                          the_bkg_full_model is None):
                        cmd = "set_bkg_source(%s, %s, bkg_id=%s)" % (cmd_id, the_bkg_source.name, cmd_bkg_id)
                    elif (the_bkg_source is not None and
                          the_bkg_full_model is not None):
                        if (repr(the_bkg_source) == repr(the_bkg_full_model)):
                            cmd = "set_bkg_full_model(%s, %s, bkg_id=%s)" % (cmd_id, the_bkg_full_model.name, cmd_bkg_id)
                        else:
                            cmd = "set_bkg_source(%s, %s, bkg_id=%s)" % (cmd_id, the_bkg_source.name, cmd_bkg_id)
                    else:
                        # You can't actually get here
                        cmd=""
                        pass
                    _send_to_outfile(cmd, outfile)
                    _send_to_outfile("", outfile)

            except:
                pass

        # Save XSPEC settings if XSPEC module has been loaded.
        if (hasattr(sherpa.astro, "xspec")):
            _send_to_outfile("\n######### XSPEC Module Settings\n", outfile)
            xspec_state = sherpa.astro.xspec.get_xsstate()

            cmd = "set_xschatter(%d)" % xspec_state["chatter"]
            _send_to_outfile(cmd, outfile)
            cmd = "set_xsabund(\"%s\")" % xspec_state["abund"]
            _send_to_outfile(cmd, outfile)
            cmd = "set_xscosmo(%g, %g, %g)" % (xspec_state["cosmo"][0],
                                               xspec_state["cosmo"][1],
                                               xspec_state["cosmo"][2])
            _send_to_outfile(cmd, outfile)
            cmd = "set_xsxsect(\"%s\")" % xspec_state["xsect"]
            _send_to_outfile(cmd, outfile)
            for name in xspec_state["modelstrings"].keys():
                cmd = "set_xsxset(\"%s\", \"%s\")" % (name,
                                                      xspec_state["modelstrings"][name])
                _send_to_outfile(cmd, outfile)
