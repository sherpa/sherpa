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


    def clean(self):
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

    # Add ability to save attributes sepcific to the astro package.
    # Save XSPEC module settings that need to be restored.
    def save(self, filename='sherpa.save', clobber=False):
        """
        save

        SYNOPSIS
           Save the current Sherpa session to a pickled file

        SYNTAX

        Arguments:
           filename   - name of save file
                        default = 'sherpa.save'           

           clobber    - clobber the file if it exists
                        default = False

        Returns:
           None

        DESCRIPTION
           Save the current Sherpa session information to a pickled file
           to be restored for future use.

        SEE ALSO
           restore, clean
        """
        if (hasattr(sherpa.astro, "xspec")):
            self._xspec_state = sherpa.astro.xspec.get_xsstate()
        else:
            self._xspec_state = None
        sherpa.ui.utils.Session.save(self, filename, clobber)

    def restore(self, filename='sherpa.save'):
        """
        restore

        SYNOPSIS
           Restore a previous Sherpa session from a pickled file

        SYNTAX

        Arguments:
           filename   - name of saved file
                        default = 'sherpa.save'           

        Returns:
           None

        DESCRIPTION
           Restore previous Sherpa session information from a pickled file
           for continued use.

        SEE ALSO
           save, clean
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


    def show_bkg(self, id=None, bkg_id=None, outfile=None, clobber=False):
        """
        show_bkg

        SYNOPSIS
           Show the Sherpa PHA background data set

        SYNTAX

        Arguments:
           id       - data id
                      default = All available data ids

           bkg_id   - background data id
                      default = All available background data ids per data id

           outfile   - filename to capture the output
                      default = None

           clobber  - overwrite outfile if exists
                      default = False

        Returns:
           None

        DESCRIPTION
           Show the Sherpa PHA background data set by data id and background
           data id

           Examples:
              show_bkg(1, 2)

              show_bkg(outfile="bkg.out", clobber=True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           show_all, show_data, show_bkg_model, show_bkg_source
        """
        all = ''
        all += self._get_show_bkg(id, bkg_id)
        _send_to_pager(all, outfile, clobber)


    def show_bkg_source(self, id=None, bkg_id=None, outfile=None, clobber=False):
        """
        show_bkg_source

        SYNOPSIS
           Show the background source model (unconvolved)

        SYNTAX

        Arguments:
           id       - data id
                      default = All available data ids

           bkg_id   - background data id
                      default = All available background data ids per data id

           outfile  - filename to capture the output
                      default = None

           clobber  - overwrite outfile if exists
                      default = False

        Returns:
           None

        DESCRIPTION
           Show the background unconvolved source model by data id
           and background data id

           Examples:
              show_bkg_source(1, 2)

              show_bkg_source(outfile="bkg_src.out", clobber=True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           show_all, show_data, show_bkg_model, show_bkg
        """
        all = ''
        all += self._get_show_bkg_source(id, bkg_id)
        _send_to_pager(all, outfile, clobber)


    def show_bkg_model(self, id=None, bkg_id=None, outfile=None, clobber=False):
        """
        show_bkg_model

        SYNOPSIS
           Show the background model (convolved)

        SYNTAX

        Arguments:
           id       - data id
                      default = All available data ids

           bkg_id   - background data id
                      default = All available background data ids per data id

           outfile   - filename to capture the output
                      default = None

           clobber  - overwrite outfile if exists
                      default = False

        Returns:
           None

        DESCRIPTION
           Show the background convolved model by data id and
           background data id

           Examples:
              show_bkg_model(1, 2)

              show_bkg_model(outfile="bkg.out", clobber=True)

           The means of paging the text is handled with the PAGER environment
           variable.  If PAGER is not found, '/usr/bin/more' is attempted
           before error.

        SEE ALSO
           show_all, show_data, show_bkg, show_bkg_source
        """
        all = ''
        all += self._get_show_bkg_model(id, bkg_id)
        _send_to_pager(all, outfile, clobber)

    ###########################################################################
    # Data
    ###########################################################################

    #@loggable(with_id=True, with_name='load_data')
    def dataspace1d(self, start, stop, step=1, numbins=None,
                    id=None, bkg_id=None, dstype=sherpa.data.Data1DInt):
        """
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

    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_arrays(self, id, *args):
        """
        load_arrays
        
        SYNOPSIS
           Load NumPy arrays into a dataset

        SYNTAX

        Arguments:
           id         - data id
        
           array0     - first NumPy array | first CrateData obj

           ...

           arrayN     - last NumPy array | last CrateData obj

           dstype     - dataset type desired
                        default = Data1D

        Returns:
           None

        DESCRIPTION
           Load NumPy arrays into a Sherpa dataset by data id or load CrateData
           objects into a Sherpa dataset by data id.  The list can include both
           NumPy arrays and CrateData objects together.

        SEE ALSO
           unpack_pha, unpack_arf, unpack_rmf, unpack_image, unpack_data
        """
        self.set_data(id, self.unpack_arrays(*args))

    def unpack_table(self, filename, ncols=2, colkeys=None, dstype=Data1D):
        """
        unpack_table

        SYNOPSIS
           Read data into a dataset

        SYNTAX

        Arguments:
           filename   - filename and path | TABLECrate obj | PyFITS HDUList obj

        Keyword Arguments:
           ncols      - number of columns
                        default = 2

           colkeys    - list of column names
                      - vector columns return additional arrays
                        default = None

           dstype     - dataset type desired
                        default = Data1D

        Returns:
           Sherpa dataset

        DESCRIPTION
           Read tabular data from a FITS or column-based text file into
           a Sherpa dataset given a filename and path or read in data from a
           Crate into a Sherpa dataset given a TABLECrate object or read in
           data from a HDUList into a Sherpa dataset.

        SEE ALSO
           unpack_pha, unpack_arf, unpack_rmf, unpack_image, unpack_data
        """
        return sherpa.astro.io.read_table(filename, ncols, colkeys, dstype)

    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_table(self, id, filename=None, ncols=2, colkeys=None,
                   dstype=Data1D):
        """
        load_table

        SYNOPSIS
           Load data by id

        SYNTAX

        Arguments:
           id         - data id
                        default = default data id

           filename   - filename and path | TABLECrate obj | PyFITS HDUList obj

        Keyword Arguments:
           ncols      - number of columns
                        default = 2

           colkeys    - list of column names
                      - vector columns return additional arrays
                        default = None

           dstype     - dataset type desired
                        default = Data1D

        Returns:
           None

        DESCRIPTION
           Load tabular data from a FITS or column-based text file into
           a Sherpa dataset given a filename and path by data id or load in
           data from a Crate into a Sherpa dataset given a TABLECrate object
           by data id or load in data from a HDUList into a Sherpa dataset
           by data id.
        
        SEE ALSO
           load_pha, load_arf, load_rmf, load_data, load_image,
           load_bkg
        """
        if filename is None:
            id, filename = filename, id
        
        self.set_data(id, self.unpack_table(filename, ncols, colkeys, dstype))
        
    def unpack_ascii(self, filename, ncols=2, colkeys=None,
                     dstype=Data1D, sep=' ', comment='#'):
        """
        unpack_ascii
        
        SYNOPSIS
           Read ASCII data into a dataset
        
        SYNTAX

        Arguments:
           filename   - filename and path

        Keyword Arguments:
           ncols      - number of columns
                        default = 2

           colkeys    - list of column names
                        default = None

           dstype     - dataset type desired
                        default = Data1D

           sep        - column separating character
                        default = ' '

           comment    - comment character
                        default = '#'

        Returns:
           Sherpa dataset

        DESCRIPTION
           Read tabular data from a column-based text file into a Sherpa
           dataset given a filename and path.

        SEE ALSO
           unpack_pha, unpack_arf, unpack_rmf, unpack_image, unpack_data,
           unpack_table
        """
        return sherpa.astro.io.read_ascii(filename, ncols, colkeys, dstype,
                                          sep=sep, comment=comment)
    
    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_ascii(self, id, filename=None, ncols=2, colkeys=None,
                   dstype=Data1D, sep=' ', comment='#'):
        """
        load_ascii

        SYNOPSIS
           Load ASCII data by id

        SYNTAX

        Arguments:
           id         - data id
                        default = default data id

           filename   - filename and path

        Keyword Arguments:
           ncols      - number of columns
                        default = 2

           colkeys    - list of column names
                        default = None

           dstype     - dataset type desired
                        default = Data1D

           sep        - column separating character
                        default = ' '

           comment    - comment character
                        default = '#'

        Returns:
           None

        DESCRIPTION
           Load tabular data from a column-based text file into a Sherpa
           dataset given a filename and path by data id.
        
        SEE ALSO
           load_pha, load_arf, load_rmf, load_data, load_image,
           load_bkg, load_table
        """
        if filename is None:
            id, filename = filename, id
            
        self.set_data(id, self.unpack_ascii(filename, ncols=ncols,
                                            colkeys=colkeys, dstype=dstype,
                                            sep=sep, comment=comment ))
        
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

    def unpack_image(self, arg, coord='logical',
                     dstype=sherpa.astro.data.DataIMG):
        """
        unpack_image

        SYNOPSIS
           Read image data into a dataset

        SYNTAX

        Arguments:
           arg        - filename and path | IMAGECrate obj | PyFITS HDUList obj
        
           coord      - string keyword identifying coordinate system
                      - choices include: logical, image
                                         physical
                                         world, wcs
                        default = logical

           dstype     - Sherpa dataset type (DataIMG, DataIMGInt)
                        default = DataIMG

        Returns:
           Sherpa DataIMG dataset

        DESCRIPTION
           Read image data from a FITS file into a Sherpa dataset given a
           filename or read in image data from a Crate into a Sherpa dataset
           given a IMAGECrate object or read in image data from a HDUList into
           a Sherpa dataset.

        SEE ALSO
           unpack_pha, unpack_arf, unpack_rmf, unpack_table, unpack_data
        """
        return sherpa.astro.io.read_image(arg, coord, dstype)

    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_image(self, id, arg=None, coord='logical',
                     dstype=sherpa.astro.data.DataIMG):
        """
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

    def unpack_pha(self, arg, use_errors=False):
        """
        unpack_pha

        SYNOPSIS
           Read PHA data into a dataset

        SYNTAX

        Arguments:
           arg        - filename and path | PHACrate obj | PyFITS HDUList obj

           use_errors - flag to use errors
                        default = False

        Returns:
           List or instance of Sherpa DataPHA dataset(s)

        DESCRIPTION
           Read PHA data from a FITS file into a Sherpa dataset given a
           filename or PHACrate object or PyFITS HDUList object.

        SEE ALSO
           unpack_image, unpack_arf, unpack_rmf, unpack_table, unpack_data
        """
        use_errors = sherpa.utils.bool_cast(use_errors)
        return sherpa.astro.io.read_pha(arg, use_errors)


    def unpack_bkg(self, arg, use_errors=False):
        """
        unpack_bkg

        SYNOPSIS
           Read background PHA data into a dataset

        SYNTAX

        Arguments:
           arg        - filename and path | PHACrate obj | PyFITS HDUList obj

           use_errors - flag to use errors
                        default = False

        Returns:
           List or instance of Sherpa background DataPHA dataset(s)

        DESCRIPTION
           Read background PHA data from a FITS file into a Sherpa dataset 
           given a filename or PHACrate object or PyFITS HDUList object.

        SEE ALSO
           unpack_image, unpack_arf, unpack_rmf, unpack_table, unpack_data
        """
        use_errors = sherpa.utils.bool_cast(use_errors)
        return sherpa.astro.io.read_pha(arg, use_errors, True)


    #@loggable(with_id=True, with_keyword='arg', with_name='load_data')
    def load_pha(self, id, arg=None, use_errors=False):
        """
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


    def load_grouping(self, id, filename=None, bkg_id=None, *args, **kwargs):
        """
        load_grouping

        SYNOPSIS
           Load the dataset grouping flags from file

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filename   - filename with path

           bkg_id     - background data id
                        default = default background data id

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
           Load the grouping flags for a dataset from file by data id.

        EXAMPLE
           load_grouping("data.dat", colkeys=["GROUPS"])

        SEE ALSO
            set_grouping
        """
        if filename is None:
            id, filename = filename, id

        self.set_grouping(id,
            self._read_user_model(filename, *args, **kwargs)[1], bkg_id=bkg_id)

    def load_quality(self, id, filename=None, bkg_id=None, *args, **kwargs):
        """
        load_quality

        SYNOPSIS
           Load the dataset quality flags from file

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filename   - filename with path

           bkg_id     - background data id
                        default = default background data id

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
           Load the quality flags for a dataset from file by data id.

        EXAMPLE
           load_quality("data.dat", colkeys=["GROUPS"])

        SEE ALSO
            set_quality
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


    def load_staterror(self, id, filename=None, bkg_id=None, *args, **kwargs):
        """
        load_staterror

        SYNOPSIS
           Load the statistical errors for a dataset from file

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filename   - filename with path

           bkg_id     - background data id
                        default = default background data id

           ncols      - number of columns to read from
                        default = 2

           colkeys    - column keys
                        default = None

        Returns:
           None

        DESCRIPTION
           Load the statistical error for a dataset from file by data id.  
           Users can specify the column name by using the colkeys argument to 
           set the statistical error.

        EXAMPLE
           load_staterror("data.dat", colkeys=["STAT_ERR"])

           load_staterror("bkg.fits", bkg_id=1, colkeys=["STAT_ERR"])

        SEE ALSO
           load_syserror, set_staterror, set_syserror
        """
        if filename is None:
            id, filename = filename, id

        self.set_staterror(id,
            self._read_user_model(filename, *args, **kwargs)[1], bkg_id=bkg_id)

    def load_syserror(self, id, filename=None, bkg_id=None, *args, **kwargs):
        """
        load_syserror

        SYNOPSIS
           Load the systematic errors for a dataset from file

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filename   - filename with path

           bkg_id     - background data id
                        default = default background data id

           ncols      - number of columns to read from
                        default = 2

           colkeys    - column keys
                        default = None

        Returns:
           None

        DESCRIPTION
           Load the systematic error for a dataset from file by data id and by 
           bkg_id.  Users can specify the column name by using the colkeys 
           argument to set the systematic error.

        EXAMPLE
           load_syserror("data.dat", colkeys=["SYS_ERR"])

           load_syserror("bkg.fits", bkg_id=1, colkeys=["SYS_ERR"])

        SEE ALSO
           load_staterror, set_staterror, set_syserror
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


    def set_staterror(self, id, val=None, fractional=False, bkg_id=None):
        """
        set_staterror

        SYNOPSIS
           Set the statistical errors of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           val        - array or scalar error values

           fractional - use fractional portion of dependent array as error,
                        val must be a scalar value
                        default = False

           bkg_id     - background id
                        default = None

        Returns:
           None

        DESCRIPTION
           Set the statistical error of a source or background dataset by data
           id or by bkg_id.  Users can specify the entire error as an array or
           as a single value to be repeated for every bin.  Also, setting the
           fractional argument will use the single value as the fractional
           portion of the dependent array as the error.

        EXAMPLE
           set_staterror([0.040, 0.039, 0.041, ...])

           set_staterror(2, 0.04)

           set_staterror(0.05, fractional=True)

           set_staterror(0.05, bkg_id=1)

        SEE ALSO
           set_syserror, set_exposure, set_backscal, set_areascal
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


    def set_syserror(self, id, val=None, fractional=False, bkg_id=None):
        """
        set_syserror

        SYNOPSIS
           Set the systematic errors of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           val        - array or scalar error values

           fractional - use fractional portion of dependent array as error,
                        val must be a scalar value
                        default = False

           bkg_id     - background id
                        default = None

        Returns:
           None

        DESCRIPTION
           Set the systematic error of a dataset by data id.  Users can specify
           the entire error as an array or as a single value to be repeated for
           every bin.  Also, setting the fractional argument will use the single
           value as the fractional portion of the dependent array as the error.

        EXAMPLE
           set_syserror([0.040, 0.039, 0.041, ...])

           set_syserror(2, 0.04)

           set_syserror(0.05, fractional=True)

           set_syserror(0.05, bkg_id=1)

        SEE ALSO
           set_syserror, set_exposure, set_backscal, set_areascal
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


    def set_exposure(self, id, exptime=None, bkg_id=None):
        """
        set_exposure

        SYNOPSIS
           Set the source or background exposure times by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           exptime    - exposure time

           bkg_id     - background data id
                        default = None

        Returns:
           None

        DESCRIPTION
           Set the exposure time of a source PHA dataset by data id or of a 
           background data by bkg_id.

        EXAMPLE
           set_exposure(10e5)

           set_exposure(2, 10e5)

           set_exposure(1, 10e5, 1)

        SEE ALSO
           set_backscal, set_areascal
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
        """
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
        """
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


    def get_staterror(self, id=None, filter=False, bkg_id=None):
        """
        get_staterror

        SYNOPSIS
           Get the statistical errors of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filter     - apply filter
                        default = False

           bkg_id     - background id
                        default = None

        Returns:
           Statistical error array

        DESCRIPTION
           Get the statistical error of a source or background dataset by data
           id or by bkg_id.

        EXAMPLE
           get_staterror()

           get_staterror(1, True)

           get_staterror(1, 2)

        SEE ALSO
           set_syserror, set_exposure, set_backscal, set_areascal,
           get_syserror, get_exposure, get_backscal, get_areascal
        """
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        return d.get_staterror(filter, self.get_stat().calc_staterror)


    def get_syserror(self, id=None, filter=False, bkg_id=None):
        """
        get_syserror

        SYNOPSIS
           Get the systematic errors of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filter     - apply filter
                        default = False

           bkg_id     - background id
                        default = None

        Returns:
           Systematic error array

        DESCRIPTION
           Get the systematic error of a dataset by data id or bkg_id.

        EXAMPLE
           get_syserror()

           get_syserror(1, True)

           get_syserror(1, bkg_id=2)

        SEE ALSO
           set_syserror, set_exposure, set_backscal, set_areascal,
           get_syserror, get_exposure, get_backscal, get_areascal
        """
        d = self.get_data(id)
        id = self._fix_id(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        err = d.get_syserror(filter)
        if err is None or not numpy.iterable(err):
            raise sherpa.utils.err.DataErr('nosyserr', id)
        return err


    def get_error(self, id=None, filter=False, bkg_id=None):
        """
        get_error

        SYNOPSIS
           Get the total errors of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filter     - apply filter
                        default = False

           bkg_id     - background id
                        default = None

        Returns:
           Total error array

        DESCRIPTION
           Get the total error (statistical + systematic in quadrature) of a
           dataset by data id or bkg_id.

        EXAMPLE
           get_error()

           get_error(1, True)

           get_error(1, bkg_id=2)

        SEE ALSO
           set_syserror, set_exposure, set_backscal, set_areascal,
           get_syserror, get_exposure, get_backscal, get_areascal
        """
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        return d.get_error(filter, self.get_stat().calc_staterror)


    def get_indep(self, id=None, filter=False, bkg_id=None):
        """
        get_indep

        SYNOPSIS
           Get the independent grid of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           bkg_id     - background id
                        default = None

        Returns:
           Array of the independent variable

        DESCRIPTION
           Get the data set independend grid by data id or bkg_id.

        EXAMPLE
           get_indep()

           get_indep(1)

           get_indep(1, 2)

        SEE ALSO

        """
        d = self.get_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        return d.get_indep(filter=filter)


    def get_axes(self, id=None, bkg_id=None):
        """
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


    def get_dep(self, id=None, filter=False, bkg_id=None):
        """
        get_dep

        SYNOPSIS
           Get the dependent variable of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filter     - apply filter
                        default = False

           bkg_id     - background id
                        default = None

        Returns:
           Array of the dependent variable

        DESCRIPTION
           Get the dependent variable array of data set by data id or bkg_id.

        EXAMPLE
           get_dep()

           get_dep(1, True)

           get_dep(1, bkg_id=2)

        SEE ALSO

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

    def get_rate(self, id=None, filter=False, bkg_id=None):
        """
        get_rate

        SYNOPSIS
           Get the measured count rate of a dataset by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filter     - apply filter
                        default = False

           bkg_id     - background id
                        default = None

        Returns:
           count rate array

        DESCRIPTION
           Get the measured count rate of data set or background by data id or
           bkg_id.

        EXAMPLE
           get_rate()

           get_rate(1, True)

           get_rate(1, bkg_id=2)

        SEE ALSO
           get_counts
        """
        d = self._get_pha_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        old = d._rate
        d._rate=True     # return count rate for PHA
        rate = d.get_y(filter)
        d._rate=old
        return rate

    def get_specresp(self, id=None, filter=False, bkg_id=None):
        """
        get_specresp

        SYNOPSIS
           Get the effective area of a PHA spectrum by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           filter     - apply filter
                        default = False

           bkg_id     - background id
                        default = None

        Returns:
           effective area array

        DESCRIPTION
           Get the effective area array of a PHA spectrum by data set or
           background by data id or bkg_id.

        EXAMPLE
           get_specresp()

           get_specresp(1, True)

           get_specresp(1, bkg_id=2)

        SEE ALSO
           get_counts, get_rate
        """
        d = self._get_pha_data(id)
        if bkg_id is not None:
            d = self.get_bkg(id,bkg_id)
        return d.get_specresp(filter)


    def get_exposure(self, id=None, bkg_id=None):
        """
        get_exposure

        SYNOPSIS
           Get the source or background exposure times by id

        SYNTAX

        Arguments:
           id         - session data id
                        default = default data id

           bkg_id     - background data id
                        default = None

        Returns:
           None

        DESCRIPTION
           Get the exposure time of a source PHA dataset by data id or of a 
           background data by bkg_id.

        EXAMPLE
           get_exposure()

           get_exposure(1)

           get_exposure(1, 2)

        SEE ALSO
           set_backscal, set_areascal,
           get_backscal, get_areascal
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
        """
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


    def save_staterror(self, id, filename=None, bkg_id=None, ascii=True,
                    clobber=False):
        """
        save_staterror

        SYNOPSIS
           Write the data set statistical errors to file

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
           Write the data set statistical errors to file by id or background
           id.

        EXAMPLE

           save_staterror("staterror.fits")
           
           save_staterror("bkgstaterr.fits", bkg_id = 1)

           save_staterror("staterror.dat", ascii = True)

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_delchi, save_error, save_syserror
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


    def save_syserror(self, id, filename=None, bkg_id=None, ascii=True,
                    clobber=False):
        """
        save_syserror

        SYNOPSIS
           Write the data set systematic errors to file

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
           Write the data set systematic errors to file by id or background
           id.

        EXAMPLE

           save_syserror("syserror.fits")
           
           save_syserror("bkgsyserr.fits", bkg_id = 1)

           save_syserror("syserror.dat", ascii = True)

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_delchi, save_error, save_staterror
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


    def save_error(self, id, filename=None, bkg_id=None, ascii=True,
                    clobber=False):
        """
        save_error

        SYNOPSIS
           Write the total errors of a data set to file

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
           Write the total errors (statistical + systematic in quadrature) of a
           data set to file by id or background id.

        EXAMPLE

           save_error("error.fits")
           
           save_error("bkgerr.fits", bkg_id = 1)

           save_error("error.dat", ascii = True)

        SEE ALSO
           save_image, save_data, save_table, save_arrays, save_source,
           save_model, save_delchi, save_staterror, save_syserror
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


    def save_pha(self, id, filename=None, bkg_id=None, ascii=False, clobber=False):
        """
        save_pha

        SYNOPSIS
           Write PHA data by id

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
           Write PHA data to a FITS file or ASCII file from a Sherpa dataset
           by id.

        EXAMPLE

           save_pha(1, "pha.fits")

           save_pha(1, "pha.out", ascii=True)

        SEE ALSO
           save_image, save_data, save_table
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


    def save_grouping(self, id, filename=None, bkg_id=None, ascii=True, clobber=False):
        """
        save_grouping

        SYNOPSIS
           Write PHA grouping flags by id

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
           Write PHA grouping flags to a FITS file or ASCII file from a
           Sherpa dataset by id.

        EXAMPLE

           save_grouping(1, "grouping.fits")

           save_grouping(1, "grouping.out", ascii=True)

        SEE ALSO
           save_image, save_data, save_table
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


    def save_quality(self, id, filename=None, bkg_id=None, ascii=True, clobber=False):
        """
        save_quality

        SYNOPSIS
           Write PHA quality flags by id

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
           Write PHA quality flags to a FITS file or ASCII file from a
           Sherpa dataset by id.

        EXAMPLE

           save_quality(1, "quality.fits")

           save_quality(1, "quality.out", ascii=True)

        SEE ALSO
           save_image, save_data, save_table
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


    def pack_pha(self, id=None):
        """
        pack_pha

        SYNOPSIS
           Pack PHA data by id

        SYNTAX

        Arguments:
           id         - dataset ID
                        default = default data id

        Returns:
           PHACrate or PyFITS HDU list

        DESCRIPTION
           Pack up PHA data from a Sherpa dataset by id.

        SEE ALSO
           pack_image, pack_data, pack_table
        """
        return sherpa.astro.io.pack_pha(self._get_pha_data(id))
    
    def pack_image(self, id=None):
        """
        pack_image

        SYNOPSIS
           Pack image data by id

        SYNTAX

        Arguments:
           id         - dataset ID
                        default = default data id

        Returns:
           IMAGECrate or PyFITS HDU list

        DESCRIPTION
           Pack up image data from a Sherpa dataset by id.

        SEE ALSO
           pack_pha, pack_data, pack_table
        """
        return sherpa.astro.io.pack_image(self.get_data(id))

    def pack_table(self, id=None):
        """
        pack_table

        SYNOPSIS
           Pack tabular data by id

        SYNTAX

        Arguments:
           id         - dataset ID
                        default = default data id

        Returns:
           TABLECrate or PyFITS HDU list

        DESCRIPTION
           Pack up tabular data  from a Sherpa dataset by id.

        SEE ALSO
           pack_pha, pack_data, pack_image
        """
        return sherpa.astro.io.pack_table(self.get_data(id))


    #def _check_resp_id(id):
    #    if (id is not None) and (not self._valid_id(id)):
    #        raise ArgumentTypeError('response identifiers must be integers ' +
    #                                'or strings')

    #@loggable(with_id=True)
    def get_arf(self, id=None, resp_id=None, bkg_id=None):
        """
        get_arf

        SYNOPSIS
           Return an ARF1D model by data id and response id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           resp_id   - response id, if multiple responses exist
                       default = default response id

           bkg_id    - background id, if background response(s) exist
                       default = None

        Returns:
           Sherpa ARF1D model

        DESCRIPTION
           Return a ARF1D model containing ancillary response data
           given a data id and a response id.

        SEE ALSO
           set_arf, unpack_arf, load_arf
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


    def set_arf(self, id, arf=None, resp_id=None, bkg_id=None):
        """
        set_arf

        SYNOPSIS
           Set an ARF dataset by data id and response id

        SYNTAX

        Arguments:
           id        - dataset id
                       default = default data id

           arf       - Sherpa DataARF dataset
                       see get_arf for more info

           resp_id   - response id, if multiple responses exist
                       default = default response id

           bkg_id    - background id, if background response(s) exist
                       default = None

        Returns:
           None

        DESCRIPTION
           Set a dataset containing ancillary response data
           by a data id and a response id.

        SEE ALSO
           get_arf, unpack_arf, load_arf
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


    def unpack_arf(self, arg):
        """
        unpack_arf

        SYNOPSIS
           Read ARF data into a dataset

        SYNTAX

        Arguments:
           arg       - filename and path | ARFCrate obj | PyFITS HDUList obj

        Returns:
           Sherpa DataARF dataset

        DESCRIPTION
           Read a FITS file containing ancillary response data given
           a filename or read data from a ARFCrate object into a Sherpa dataset
           or read data from a PyFITS HDUList object into a Sherpa dataset.

        SEE ALSO
           get_arf, set_arf, load_arf
        """
        return sherpa.astro.instrument.ARF1D(sherpa.astro.io.read_arf(arg))

    #@loggable(with_id=True, with_keyword='arg')
    def load_arf(self, id, arg=None, resp_id=None, bkg_id=None):
        """
        load_arf

        SYNOPSIS
           Loads ARF data by id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           arg       - filename with path | ARFCrate obj | PyFITS HDUList obj

           resp_id   - response id, if multiple responses exist
                       default = default response id

           bkg_id    - background id, if background response(s) exist
                       default = None

        Returns:
           None

        DESCRIPTION
           Load a FITS file containing ancillary response data given
           a filename by data id and response id or read data from a ARFCrate
           object into a Sherpa dataset by data id and response id or read data
           from a PyFITS HDUList object into a Sherpa dataset by data id and
           response id.

        SEE ALSO
           get_arf, set_arf, unpack_arf
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

    #@loggable(with_id=True, with_keyword='arg')
    def load_bkg_arf(self, id, arg=None):
        """
        load_bkg_arf

        SYNOPSIS
           Loads bkg ARF data by id and default bkg_id and default resp_id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           arg       - filename with path | ARFCrate obj | PyFITS HDUList obj

        Returns:
           None

        DESCRIPTION
           Load a FITS file containing ancillary response data given a filename
           by data id and response id or read data from a ARFCrate object or
           PyFITS HDUList object into a Sherpa dataset by data id and default
           background id and default response id.

        SEE ALSO
           get_arf, set_arf, unpack_arf
        """
        if arg is None:
            id, arg = arg, id
        bkg_id = self._get_pha_data(id).default_background_id
        resp_id = self._get_pha_data(id).primary_response_id
        self.set_arf(id, self.unpack_arf(arg), resp_id, bkg_id)

    def load_multi_arfs(self, id, filenames, resp_ids=None):
        """
        load_multi_arfs

        SYNOPSIS
           Loads multiple ARF data files by id and resp_id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           filenames - list of ARF files

           resp_ids  - list of response ids

        Returns:
           Load a list of FITS files containing ancillary response data given
           a list of filenames by data id and a list of response ids. 

        DESCRIPTION
           None

        SEE ALSO
           set_arf, get_arf, unpack_arf, load_arf
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

    #@loggable(with_id=True)
    def get_rmf(self, id=None, resp_id=None, bkg_id=None):
        """
        get_rmf

        SYNOPSIS
           Return an RMF1D model by data id and response id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           resp_id   - response id, if multiple responses exist
                       default = default response id

           bkg_id    - background id, if background response(s) exist
                       default = None

        Returns:
           Sherpa RMF1D model

        DESCRIPTION
           Return a RMF1D model containing response matrix data
           given a data id and a response id.

        SEE ALSO
           set_rmf, unpack_rmf, load_rmf
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


    def set_rmf(self, id, rmf=None, resp_id=None, bkg_id=None):
        """
        set_rmf

        SYNOPSIS
           Set an RMF dataset by data id and response id

        SYNTAX

        Arguments:
           id        - dataset id
                       default = default data id

           arf       - Sherpa DataRMF dataset
                       see get_rmf for more info

           resp_id   - response id, if multiple responses exist
                       default = default response id

           bkg_id    - background id, if background response(s) exist
                       default = None

        Returns:
           None

        DESCRIPTION
           Set a dataset containing response matrix data
           by a data id and a response id.

        SEE ALSO
           get_rmf, unpack_rmf, load_rmf
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


    def unpack_rmf(self, arg):
        """
        unpack_rmf

        SYNOPSIS
           Read RMF data into a dataset

        SYNTAX

        Arguments:
           arg       - filename and path | RMFCrate obj | PyFITS HDUList obj

        Returns:
           Sherpa DataRMF dataset

        DESCRIPTION
           Read a FITS file containing response matrix data given
           a filename or read data from a RMFCrate object into a Sherpa dataset
           or read data from a PyFITS HDUList object into a Sherpa dataset.

        SEE ALSO
           get_rmf, set_rmf, load_rmf
        """
        return sherpa.astro.instrument.RMF1D(sherpa.astro.io.read_rmf(arg))

    #@loggable(with_id=True, with_keyword='arg')
    def load_rmf(self, id, arg=None, resp_id=None, bkg_id=None):
        """
        load_rmf

        SYNOPSIS
           Loads RMF data by id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           arg       - filename with path of RMF file

           resp_id   - response id, if multiple responses exist
                       default = default response id

           bkg_id    - background id, if background response(s) exist
                       default = None

        Returns:
           None

        DESCRIPTION
           Load a FITS file containing response matrix data given
           a filename by data id and response id or read data from a RMFCrate
           object into a Sherpa dataset by data id and response id or read data
           from a PyFITS HDUList object into a Sherpa dataset by data id and
           response id.

        SEE ALSO
           get_rmf, set_rmf, unpack_rmf
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

    #@loggable(with_id=True, with_keyword='arg')
    def load_bkg_rmf(self, id, arg=None):
        """
        load_bkg_rmf

        SYNOPSIS
           Loads bkg RMF data by id using default bkg_id and default resp_id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           arg       - filename with path of RMF file

        Returns:
           None

        DESCRIPTION
           Load a FITS file containing response matrix data given
           a filename by data id and response id or read data from a RMFCrate
           object or a PyFITS HDUList object into a Sherpa background dataset
           by data id and default response id and default background id.

        SEE ALSO
           get_bkg_rmf, set_rmf, unpack_rmf
        """
        if arg is None:
            id, arg = arg, id
        bkg_id = self._get_pha_data(id).default_background_id
        resp_id = self._get_pha_data(id).primary_response_id
        self.set_rmf(id, self.unpack_rmf(arg), resp_id, bkg_id)

    def load_multi_rmfs(self, id, filenames, resp_ids=None):
        """
        load_multi_rmfs

        SYNOPSIS
           Loads multiple RMF data files by id and resp_id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           filenames - list of RMF files

           resp_ids   - list of response ids

        Returns:
           Load a list of FITS files containing response data given
           a list of filenames by data id and a list of response ids. 

        DESCRIPTION
           None

        SEE ALSO
           set_rmf, get_rmf, unpack_rmf, load_rmf
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
        """
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

    def set_bkg(self, id, bkg=None, bkg_id=None):
        """
        set_bkg

        SYNOPSIS
           Set a background PHA dataset by data id and 
           background id

        SYNTAX

        Arguments:
           id        - dataset id
                       default = default data id

           bkg       - Sherpa DataPHA dataset
                       see get_bkg for more info

           bkg_id    - background id, if multiple bkgs exist
                       default = default background id

        Returns:
           None

        DESCRIPTION
           Set a dataset containing background PHA data
           by a data id and a background id.

        SEE ALSO
           get_bkg, unpack_bkg, load_bkg
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


    def list_bkg_ids(self, id=None):
        """
        list_bkg_ids

        SYNOPSIS
           List the available Sherpa background ids for a data set by data id

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id

        Returns:
           list of background ids

        DESCRIPTION
           The Sherpa background id ties background data sets to a source data
           set easily referenced by data id.  The id can be a user
           defined string or integer.

        SEE ALSO
           get_bkg, set_bkg
        """
        #return self._get_pha_data(id).background_ids
        return self._get_pha_data(id)._backgrounds.keys()

    def list_response_ids(self, id=None, bkg_id=None):
        """
        list_response_ids

        SYNOPSIS
           List the available Sherpa response ids for a data set by data id

        SYNTAX

        Arguments:
           id        - Sherpa data id
                       default = default data id
        
           bkg_id    - Sherpa background id
                       default = None

        Returns:
           list of response ids

        DESCRIPTION
           The Sherpa response id ties ARF and RMF data sets to a source or
           background data set easily referenced by data id or background id.
           The id can be a user defined string or integer.

        SEE ALSO
           get_arf, set_arf, get_rmf, set_rmf, load_arf, load_rmf
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
    ### DOC-TODO: this is a function that doesn't really match the normal
    ###           Python parameter rules
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
           identifier is used, as set by `set_default_id`.

        Returns
        -------
        quantity : { 'channel', 'energy', 'wavelength' }

        Raises
        ------
        sherpa.utils.err.IdentifierErr
           If the `id` argument is not recognized.

        See Also
        --------
        get_default_id : Return the default data set identifier. 
        set_analysis : Change the analysis setting.
        set_default_id : Change the default data set identifier. 

        """
        return self._get_pha_data(id).get_analysis()

    #@loggable(with_id=True, with_keyword='coord')
    def set_coord(self, id, coord=None):
        """
        set_coord

        SYNOPSIS
           Set the coordinate system by data id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           coord     - string keyword identifying coordinate system

        Returns:
           None

        DESCRIPTION
           Set the coordinate system of a Sherpa DataIMG dataset
           by data id.  Choices include logical, physical, and world
           coordinates.  Alias for logical is image.  Alias for world
           is wcs.

           * 'logical' or 'image'

           * 'physical'

           * 'world' or 'wcs'

        SEE ALSO
           notice2d, notice2d_id, notice2d_image, ignore2d, ignore2d_id,
           ignore2d_image, get_coord
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


    def get_coord(self, id=None):
        """
        get_coord

        SYNOPSIS
           Get the coordinate system by data id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

        Returns:
           Coord

        DESCRIPTION
           Get the coordinate system of a Sherpa DataIMG dataset
           by data id.  Return values include logical, physical, and world
           coordinates.  Alias for logical is image.  Alias for world
           is wcs.

           * 'logical'

           * 'physical'

           * 'world'

        SEE ALSO
           notice2d, notice2d_id, notice2d_image, ignore2d, ignore2d_id,
           ignore2d_image, set_coord
        """
        return self._get_img_data(id).coord


    def ignore_bad(self, id=None, bkg_id=None):
        """
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
        """
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
        """
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
        """
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

    #@loggable(with_id=True)
    def group(self, id=None, bkg_id=None):
        """
        group

        SYNOPSIS
           Turn grouping ON

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           bkg_id    - background id
                       default = default bkg id

        Returns:
           None

        DESCRIPTION
           Set grouping boolean to True in a Sherpa DataPHA
           dataset by data id or background dataset by bkg id
           utilizing native grouping flags.

        SEE ALSO
           set_grouping, ungroup
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


    def set_grouping(self, id, val=None, bkg_id=None):
        """
        set_grouping

        SYNOPSIS
           Apply user defined grouping flags by data id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           val       - properly sized integer array of grouping flags

           bkg_id    - background id
                       default = default bkg id

        Returns:
           None

        DESCRIPTION
           Override native grouping flags (if available) of Sherpa
           DataPHA dataset by data id or background by bkg id to user
           a user defined array of integers.           

        SEE ALSO
           ungroup, group
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


    def get_grouping(self, id=None, bkg_id=None):
        """
        get_grouping

        SYNOPSIS
           Retrieve the grouping flags by data id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           bkg_id    - background id
                       default = default bkg id

        Returns:
           grouping flags array

        DESCRIPTION
           Obtain the native grouping flags (if available) of a Sherpa
           DataPHA dataset by data id or background by bkg id.

        SEE ALSO
           ungroup, group, load_grouping, set_grouping
        """

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        return data.grouping


    def set_quality(self, id, val=None, bkg_id=None):
        """
        set_quality

        SYNOPSIS
           Apply user defined quality flags by data id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           val       - properly sized integer array of quality flags

           bkg_id    - background id
                       default = default bkg id

        Returns:
           None

        DESCRIPTION
           Override native quality flags (if available) of Sherpa
           DataPHA dataset by data id or background by bkg id to user
           a user defined array of integers.           

        SEE ALSO
           ungroup, group
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


    def get_quality(self, id=None, bkg_id=None):
        """
        get_quality

        SYNOPSIS
           Retrieve the quality flags by data id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           bkg_id    - background id
                       default = default bkg id

        Returns:
           quality flags array

        DESCRIPTION
           Obtain the native quality flags (if available) of a Sherpa
           DataPHA dataset by data id or background by bkg id.

        SEE ALSO
           ungroup, group, load_quality, set_quality
        """

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        return data.quality

    #@loggable(with_id=True)
    def ungroup(self, id=None, bkg_id=None):
        """
        ungroup

        SYNOPSIS
           Turn grouping OFF

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           bkg_id    - background id
                       default = default bkg id

        Returns:
           None

        DESCRIPTION
           Set grouping boolean to False in a Sherpa DataPHA
           dataset by data id or background by bkg id utilizing
           native grouping flags.

        SEE ALSO
           set_grouping, group
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

    #@loggable(with_id=True, with_keyword='num')
    def group_bins(self, id, num=None, bkg_id=None, tabStops=None):
        """
        group_bins

        SYNOPSIS
           Create and set grouping flags by number of bins with equal-widths

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           num       - number of groups

           bkg_id    - background id
                       default = default bkg id

           tabStops  - integer array of noticed channels (1 means ignore)
                       default = None

        Returns:
           None

        DESCRIPTION
           Creates and sets grouping flags on a PHA spectrum data set by data ID
           using a number of groups with equal-widths.  Resetting the grouping
           flags clears any filters already in place.

        SEE ALSO
           group_width, group_snr, group_adapt, group_adapt_snr
        """
        if num is None:
            id, num = num, id

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_bins(num, tabStops)

    #@loggable(with_id=True, with_keyword='num')
    def group_width(self, id, num=None, bkg_id=None, tabStops=None):
        """
        group_width

        SYNOPSIS
           Create and set grouping flags by a bin width.

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           num       - bin width

           bkg_id    - background id
                       default = default bkg id

           tabStops  - integer array of noticed channels (1 means ignore)
                       default = None

        Returns:
           None

        DESCRIPTION
           Creates and sets grouping flags on a PHA spectrum data set by data ID
           using a specific bin width.  Resetting the grouping
           flags clears any filters already in place.

        SEE ALSO
           group_bins, group_snr, group_adapt, group_adapt_snr
        """
        if num is None:
            id, num = num, id

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_width(num, tabStops)

    #@loggable(with_id=True, with_keyword='num')
    def group_counts(self, id, num=None, bkg_id=None,
                     maxLength=None, tabStops=None):
        """
        group_counts

        SYNOPSIS
           Create and set grouping flags using minimum number of counts per bin

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           num       - number of counts per bin

           bkg_id    - background id
                       default = default bkg id

           maxLength - number of elements that can be combined into a group
                       default = None

           tabStops  - integer array of noticed channels (1 means ignore)
                       default = None

        Returns:
           None

        DESCRIPTION
           Creates and sets grouping flags on a PHA spectrum data set by data ID
           using a minimum number of counts per bin.  Resetting the grouping
           flags clears any filters already in place.

        SEE ALSO
           group_snr, group_adapt, group_adapt_snr
        """
        if num is None:
            id, num = num, id

        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_counts(num, maxLength, tabStops)

    #@loggable(with_id=True, with_keyword='snr')
    def group_snr(self, id, snr=None, bkg_id=None,
                        maxLength=None, tabStops=None, errorCol=None):
        """
        group_snr

        SYNOPSIS
           Create and set grouping flags so each group has a signal-to-noise
           ratio of at least snr

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           snr       - minimum signal-to-noise ratio

           bkg_id    - background id
                       default = default bkg id

           maxLength - number of elements that can be combined into a group
                       default = None

           tabStops  - integer array of noticed channels (1 means ignore)
                       default = None

           errorCol  - gives the error for each element of the original array
                       default = None

        Returns:
           None

        DESCRIPTION
           Creates and sets grouping flags on a PHA spectrum data set by data ID
           using a minimum signal-to-noise, snr, for each group.  Resetting the
           grouping flags clears any filters already in place.

        SEE ALSO
           group_counts, group_adapt, group_adapt_snr
        """
        if snr is None:
            id, snr = snr, id
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_snr(snr, maxLength, tabStops, errorCol)

    #@loggable(with_id=True, with_keyword='min')
    def group_adapt(self, id, min=None, bkg_id=None,
                     maxLength=None, tabStops=None):
        """
        group_adapt

        SYNOPSIS
           Create and set grouping flags adaptively so that each group contains
           at least min counts.

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           min       - minimum number of counts

           bkg_id    - background id
                       default = default bkg id

           maxLength - number of elements that can be combined into a group
                       default = None

           tabStops  - integer array of noticed channels (1 means ignore)
                       default = None

        Returns:
           None

        DESCRIPTION
           Creates and sets grouping flags adaptively on a PHA spectrum data
           set by data ID using a minimum number of counts for each group.
           Resetting the grouping flags clears any filters already in place.

        SEE ALSO
           group_counts, group_snr, group_adapt_snr
        """
        if min is None:
            id, min = min, id
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_adapt(min, maxLength, tabStops)

    #@loggable(with_id=True, with_keyword='min')
    def group_adapt_snr(self, id, min=None, bkg_id=None,
                        maxLength=None, tabStops=None, errorCol=None):
        """
        group_adapt_snr

        SYNOPSIS
           Create and set grouping flags adaptively so that each group contains
           a signal-to-noise ratio of at least min.

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           min       - minimum number of counts

           bkg_id    - background id
                       default = default bkg id

           maxLength - number of elements that can be combined into a group
                       default = None

           tabStops  - integer array of noticed channels (1 means ignore)
                       default = None

           errorCol  - gives the error for each element of the original array
                       default = None

        Returns:
           None

        DESCRIPTION
           Creates and sets grouping flags adaptively on a PHA spectrum data
           set by data ID using a signal-to-noise ratio of at least min for each
           group.  Resetting the grouping flags clears any filters already in
           place.

        SEE ALSO
           group_counts, group_adapt, group_snr
        """
        if min is None:
            id, min = min, id
        data = self._get_pha_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        data.group_adapt_snr(min, maxLength, tabStops, errorCol)

    #@loggable(with_id=True)
    def subtract(self, id=None):
        """
        subtract

        SYNOPSIS
           Subtract background counts

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id
     
        Returns:
           None

        DESCRIPTION
           Subtract background counts from total counts according
           to the following equation:
           
           Measured = Total  - Back  * Data Exposure * Data Area
           Counts     Counts   Counts  Back Exposure   Back Area

        SEE ALSO
           unsubtract        
        """
        if (self._get_pha_data(id).subtracted is True):
            raise DataErr('subtractset', 'data set', str(self._fix_id(id)), 'True')
        self._get_pha_data(id).subtract()

    #@loggable(with_id=True)
    def unsubtract(self, id=None):
        """
        unsubtract

        SYNOPSIS
           Ignore subtraction of background counts

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

        Returns:
           None

        DESCRIPTION
           Ignore subtraction of background counts total counts
           according to the following equation:
           
           Measured = Total  - Back  * Data Exposure * Data Area
           Counts     Counts   Counts  Back Exposure   Back Area

        SEE ALSO
           unsubtract
        """
        if (self._get_pha_data(id).subtracted is False):
            raise DataErr('subtractset', 'data set', str(self._fix_id(id)), 'False')
        self._get_pha_data(id).unsubtract()

    def fake_pha(self, id, arf, rmf, exposure, backscal=None, areascal=None,
                 grouping=None, grouped=False, quality=None, bkg=None):
        """
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

    #@loggable(with_id=True, with_keyword='model')
    def set_bkg_model(self, id, model=None, bkg_id=None):
        """
        set_bkg_model

        SYNOPSIS
           Set a Sherpa background source model by data id
           and bkg id

        SYNTAX

        Arguments:
           id        - data id
                       default = default data id

           model     - Sherpa bkg source model

           bkg_id    - bkg id, if multiple bkgs exist
                       default = default bkg id

        Returns:
           None

        DESCRIPTION
           Add a Sherpa background source model to the list 
           of current background source models by data id 
           and background id.

        SEE ALSO
           get_bkg_source, delete_bkg_model, set_bkg_full_model,
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


    def plot_arf(self, id=None, resp_id=None, **kwargs):
        """
        plot_arf

        SYNOPSIS
           Plot ancillary response data

        SYNTAX

        Arguments:
           id       - data id
                      default = default data id

           resp_id  - response id, if multiple response exist
                      default  = default response id

           replot   - replot calculated arrays
                      default = False

           overplot - Plot data without clearing previous plot
                      default = False

        Returns:
           None

        DESCRIPTION
           Visualize the ancillary response (effective area) in a 1D plot by
           data id or response id.

        SEE ALSO
           plot_data, plot_source, plot_bkg, plot_model
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

    def eqwidth(self, src, combo, id=None, lo=None, hi=None, bkg_id=None):
        """
        eqwidth

        SYNOPSIS
           Get equivalent width

        SYNTAX

        Arguments:
           src      - continuum, type Sherpa model

           combo    - continuum plus emission line, type Sherpa model

           id       - data id
                      default = default data id
                      
           lo       - lower bin boundry
                      default = None

           hi       - upper bin boundry
                      default = None

           bkg_id   - bkg id
                      default = default bkg id

        Returns:
           eqwidth value

        DESCRIPTION
           Compute the equivalent width of an emission or
           absorption line in a source or background dataset
           by data id or background id.

        SEE ALSO
           calc_model_sum, calc_data_sum, calc_energy_flux, calc_photon_flux,
           calc_source_sum
        
        """
        data = self.get_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)

        return sherpa.astro.utils.eqwidth(data, src, combo, lo, hi)


    def calc_photon_flux(self, lo=None, hi=None, id=None, bkg_id=None):
        """
        calc_photon_flux

        SYNOPSIS
           Get the unconvolved photon flux

        SYNTAX

        Arguments:
           lo       - low limit
                      default = None

           hi       - high limit
                      default = None

           id       - data id,
                      default = default data id

           bkg_id   - bkg id
                      default = default bkg id

        Returns:
           photon flux value

        DESCRIPTION
           Calculate the unconvolved photon flux for a source
           or background dataset by data id or background id.

        SEE ALSO
           calc_energy_flux, eqwidth, calc_data_sum, calc_model_sum,
           calc_source_sum         
        """
        
        data = self.get_data(id)
        model= None

        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            model= self.get_bkg_source(id, bkg_id)
        else:
            model = self.get_source(id)
            
        return sherpa.astro.utils.calc_photon_flux(data, model, lo, hi)
    
    def calc_energy_flux(self, lo=None, hi=None, id=None, bkg_id=None):
        """
        calc_energy_flux

        SYNOPSIS
           Get the unconvolved energy flux

        SYNTAX

        Arguments:
           lo       - low limit
                      default = None

           hi       - high limit
                      default = None

           id       - data id
                      default = default data id

           bkg_id   - bkg id
                      default = default bkg id

        Returns:
           energy flux value

        DESCRIPTION
           Calculates the unconvolved energy flux for a source
           or background dataset by data id or background id.

        SEE ALSO
           calc_photon_flux, eqwidth, calc_data_sum, calc_model_sum,
           calc_source_sum
        """
        data = self.get_data(id)
        model= None

        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            model= self.get_bkg_source(id, bkg_id)
        else:
            model= self.get_source(id)
        return sherpa.astro.utils.calc_energy_flux(data, model, lo, hi)


    def calc_data_sum(self, lo=None, hi=None, id=None, bkg_id=None):
        """
        calc_data_sum

        SYNOPSIS
           Get observed data counts

        SYNTAX

        Arguments:
           lo       - low limit
                      default = None

           hi       - high limit
                      default = None

           id       - data id
                      default = default data id

           bkg_id   - bkg id, if multiple backgrounds exist
                      default = default bkg id

        Returns:
           sum value of observed counts

        DESCRIPTION
           Calculates the sum of observed counts data for a source
           or background dataset by data id or background id.

        SEE ALSO
           calc_model_sum, calc_photon_flux, calc_energy_flux, eqwidth,
           calc_source_sum, calc_data_sum2d, calc_model_sum2d
        """
        data = self.get_data(id)
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
        return sherpa.astro.utils.calc_data_sum(data, lo, hi)
            
    def calc_model_sum(self, lo=None, hi=None, id=None, bkg_id=None):
        """
        calc_model_sum

        SYNOPSIS
           Get the sum of convolved model amplitudes

        SYNTAX

        Arguments:
           lo       - low limit
                      default = None

           hi       - high limit
                      default = None

           id       - dataset ID
                      default = default data id

           bkg_id   - bkg id, if multiple backgrounds exist
                      default = default bkg id

        Returns:
           sum value of convolved model amplitudes

        DESCRIPTION
           Calculates the sum of convolved model amplitudes
           for a source or background dataset by data id or
           background id.

        SEE ALSO
           calc_data_sum, calc_photon_flux, calc_energy_flux, eqwidth,
           calc_source_sum, calc_data_sum2d, calc_model_sum2d
        """
        data = self.get_data(id)
        model= None
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            model= self.get_bkg_model(id, bkg_id)
        else:
            model= self.get_model(id)
        return sherpa.astro.utils.calc_model_sum(data, model, lo, hi)

    def calc_data_sum2d(self, reg=None, id=None):
        """
        calc_data_sum2d

        SYNOPSIS
           Get observed image counts

        SYNTAX

        Arguments:
           reg      - filename and path of region file or DM region syntax
                      default = None

           id       - data id
                      default = default data id

        Returns:
           sum value of observed image counts

        DESCRIPTION
           Calculates the sum of observed counts data for a source
           image by data id

        SEE ALSO
           calc_model_sum, calc_photon_flux, calc_energy_flux, eqwidth,
           calc_source_sum, calc_data_sum, calc_model_sum2d
        """
        data = self.get_data(id)
        return sherpa.astro.utils.calc_data_sum2d(data, reg)

    def calc_model_sum2d(self, reg=None, id=None):
        """
        calc_model_sum2d

        SYNOPSIS
           Get the sum of convolved image model amplitudes

        SYNTAX

        Arguments:
           reg      - filename and path of region file or DM region syntax
                      default = None

           id       - data id
                      default = default data id

        Returns:
           sum value of convolved image model amplitudes

        DESCRIPTION
           Calculates the sum of convolved image model amplitudes
           for a source by data id

        SEE ALSO
           calc_data_sum, calc_photon_flux, calc_energy_flux, eqwidth,
           calc_source_sum, calc_model_sum, calc_data_sum2d
        """
        data = self.get_data(id)
        model= self.get_model(id)
        return sherpa.astro.utils.calc_model_sum2d(data, model, reg)

    def calc_source_sum2d(self, reg=None, id=None):
        """
        calc_source_sum2d

        SYNOPSIS
           Get the sum of unconvolved image model amplitudes

        SYNTAX

        Arguments:
           reg      - filename and path of region file or DM region syntax
                      default = None

           id       - data id
                      default = default data id

        Returns:
           sum value of unconvolved image model amplitudes

        DESCRIPTION
           Calculates the sum of unconvolved image model amplitudes
           for a source by data id

        SEE ALSO
           calc_data_sum, calc_photon_flux, calc_energy_flux, eqwidth,
           calc_source_sum, calc_source_sum, calc_data_sum2d
        """
        data = self.get_data(id)
        src= self.get_source(id)
        return sherpa.astro.utils.calc_model_sum2d(data, src, reg)

    def calc_source_sum(self, lo=None, hi=None, id=None, bkg_id=None):
        """
        calc_source_sum

        SYNOPSIS
           Get the sum of unconvolved model amplitudes

        SYNTAX

        Arguments:
           lo       - low limit
                      default = None

           hi       - high limit
                      default = None

           id       - dataset ID
                      default = default data id

           bkg_id   - bkg id, if multiple backgrounds exist
                      default = default bkg id

        Returns:
           sum value of unconvolved model amplitudes

        DESCRIPTION
           Calculates the sum of unconvolved model amplitudes
           for a source or background dataset by data id or
           background id.

        SEE ALSO
           calc_data_sum, calc_photon_flux, calc_energy_flux, eqwidth,
           calc_model_sum
        """
        data = self.get_data(id)
        model= None
        if bkg_id is not None:
            data = self.get_bkg(id, bkg_id)
            model= self.get_bkg_source(id, bkg_id)
        else:
            model= self.get_source(id)
        return sherpa.astro.utils.calc_source_sum(data, model, lo, hi)


    def calc_kcorr(self, z, obslo, obshi, restlo=None, resthi=None,
                   id=None, bkg_id=None):
        """
        calc_kcorr

        SYNOPSIS
           Calculate the k correction 

        SYNTAX

        Arguments:
           z        - redshift (scalar or array)

           obslo    - observed-frame lower limit

           obshi    - observed-frame upper limit

           restlo   - rest-frame lower limit
                      default = obslo

           resthi   - rest-frame upper limit
                      default = obshi

           id       - dataset ID
                      default = default data id

           bkg_id   - bkg id, if multiple backgrounds exist
                      default = default bkg id

        Returns:
           k correction (scalar or array)

        DESCRIPTION
           Calculates the k correction for a spectral model,
           redshift, and energy range for a source or background
           dataset by data id or background id.

        EXAMPLE
           set_model( xsmekal.clus )
           calc_kcorr(0.5, 0.5, 2)
           1.0733301
           calc_kcorr(0.5, 0.5, 2, 2, 10)
           0.1129745

        SEE ALSO
           calc_data_sum, calc_photon_flux, calc_energy_flux, eqwidth,
           calc_model_sum
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

    def save_all(self, outfile=None, clobber=False):

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
