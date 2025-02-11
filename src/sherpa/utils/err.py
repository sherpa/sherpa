#
#  Copyright (C) 2010, 2016, 2017, 2019, 2020, 2022, 2023
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

"""
Sherpa specific exceptions
"""

__all__ = ('EstErr', 'FitErr', 'SherpaErr', 'ArgumentErr',
           'ArgumentTypeErr', 'IdentifierErr', 'NotImplementedErr',
           'ImportErr', 'ParameterErr', 'DataErr', 'PSFErr', 'InstrumentErr',
           'IOErr', 'ModelErr', 'PlotErr', 'StatErr', 'DS9Err',
           'ConfidenceErr', 'SessionErr')


class SherpaErr(Exception):
    """Base class for all Sherpa exceptions.

    Parameters
    ----------
    edict : dict
        The error dictionary. The keys are the label for the message
        and the value is a format string - using percent formats -
        that is applied.
    args
        The arguments that define the message. If args is not set then
        a generic message is used, otherwise the first element is used
        to identify the format string from edict, to which the
        remaining arguments are applied. If no match exists then the
        first argument is used as the error message.

    """

    def __init__(self, edict, *args):
        if len(args) == 0:
            super().__init__("Generic Error")
            return

        # Could trap any format errors (in case code has used an argument
        # incorrectly).
        #
        key = args[0]
        if key in edict:
            errmsg = edict[key] % args[1:]
        else:
            errmsg = key

        super().__init__(errmsg)


class ArgumentErr(ValueError, SherpaErr):

    dict = {'nosession': "file '%s' does not contain a saved Sherpa session",
            'badmethod': "'%s' is not a valid method",
            'badopt': "'%s' is not a valid option for method %s",
            'badstat': "'%s' is not a valid statistic",
            'notype': "no typename given",
            'noname': "no name given",
            'badtype': "'%s' is not a valid model type",
            'badexpr': 'invalid %s expression: %s',
            'badconf': "'%s' is not a valid confidence limit method",
            'badplottype': "'%s' is not a valid plot type",

            'nopha': 'data set %s does not contain PHA data',
            'noimg': 'data set %s does not contain IMAGE data',
            'multirsp': 'A response ID is required for each file',
            'bad': "Invalid %s: '%s'",
            'badinterval': "interval syntax requires a tuple, 'lo:hi'",
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, ArgumentErr.dict, key, *args)


class ArgumentTypeErr(TypeError, SherpaErr):

    dict = {'badarg': "'%s' must be %s",
            'intstr': 'identifiers must be integers or strings',
            'plotargs': 'not enough arguments to plot()',
            'tempplotbackend': "'%s' is not a backend class, instance, or string name",
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, ArgumentTypeErr.dict, key, *args)


class IdentifierErr(SherpaErr):

    dict = {'badid': "identifier '%s' is a reserved word",
            'getitem': '%s %s %s',
            'nodatasets': "No data sets found",
            'nomodels': "model stack is empty",
            'badidmodel': "'%s' is a model type, not allowed as a model name",
            'badidnative': "'%s' is reserved for the native Python function",
            'nomodelcmpt': "model component '%s' does not exist",
            'noplotbackend': "'%s' is not a valid plotting backend, choose from %s",
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, IdentifierErr.dict, key, *args)


class TypeErr(SherpaErr):

    dict = {'nocomplex': "ds9 cannot handle complex data",
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, TypeErr.dict, key, *args)


class RuntimeErr(SherpaErr):

    dict = {'notonpath': "Could not find %s on your PATH",
            'nowin': 'Could not open ds9 window %r; timeout',
            'cmdfail': "%r failed: %s",
            'only2d3d': "ds9 can only display 2d and 3d arrays",
            'badarr': "Array info not allowed; rejected keywords: %s",
            'badwin': "DS9Win unusable: %s",
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, RuntimeErr.dict, key, *args)


class NotImplementedErr(SherpaErr):

    dict = {'noinstanceallowed': "cannot create %s instances",
            'contourgrids': 'contours on non-uniform grids are not yet supported'
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, NotImplementedErr.dict, key, *args)


class ImportErr(SherpaErr):

    dict = {'importfailed': "failed to import %s module; %s routines are not available",
            'notsupported': '%s support is not enabled in this build of Sherpa'
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, ImportErr.dict, key, *args)


class EstErr(SherpaErr):

    dict = {'noerr4least2': 'cannot estimate confidence limits with %s',
            'nodegfreedom': 'degrees of freedom are zero or lower',
            'rstat>max': 'reduced statistic larger than %s',
            'nocov': 'Covariance matrix could not be resolved',
            'noparameter': '%s is not in a model associated with selected data sets',
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, EstErr.dict, key, *args)


class FitErr(SherpaErr):

    #
    # Make the dictionary of the error a class attribute or self.dict?
    #
    dict = {'statnotforbackgsub': '%s statistics cannot be used with background subtracted data',
            'binhas0': 'zeros found in uncertainties, consider using calculated uncertainties',
            'nobins': 'no noticed bins found in data set',
            'noclobererr': "'%s' exists, and clobber==False",
            'nothawedpar': 'model has no thawed parameters',
            'needchi2': '%s method requires a deviates array; use a chi-square  statistic', }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, FitErr.dict, key, *args)


class ParameterErr(SherpaErr):
    "Error in creating or using a model"

    dict = {'edge': "parameter %s has a %s of %g",
            'noncall': 'attempted to create %s from non-callable object of type %s',
            'alwaysint': '%s model is defined as integrated',
            'filterarray': "filter '%s' is not an array",
            'filtermismatch': "filter '%s' does not match the dimensions %s",
            'alwaysfrozen': 'parameter %s is always frozen and cannot be thawed',
            'frozennolink': 'parameter %s is always frozen and cannot be linked',
            'notlink': 'link value must be a parameter or None',
            'linkcycle': 'requested parameter link creates a cyclic reference',
            'frozen': "parameter '%s' is frozen",
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, ParameterErr.dict, key, *args)


class DataErr(SherpaErr):
    "Error in creating or using a data set"

    dict = {'ismask': "'mask' must be True, False, or a mask array",
            'notmask': 'mask excludes all data',
            'nomask': "data set '%s' has no filter",
            # mismatchn is newer and reports the size as well as the fields
            # (reported as strings to allow a value like "None" to be used)
            'mismatch': 'size mismatch between %s and %s',
            'mismatchn': 'size mismatch between %s and %s: %s vs %s',
            'notanarray': "Array must be a sequence or None",
            'notanintarray': "Array must be a sequence of integers or None",
            'not1darray': "Array must be 1D",
            'wrongaxiscount': "data set '%s' sent wrong tuple size for the independent axis: %d not %d",
            'sizenotset': "The size of '%s' has not been set",
            'typecheck': 'strings not allowed in %s list',
            'wrongdim': "data set '%s' does not contain %d-D data",
            'notimage': "data set '%s' does not contain image data",
            'nodim': "data set '%s' does not have any defined dimensions",  #  TODO: this does not appear to be used
            'zerodatasimulfit':
            "cannot create a %s instance containing no data sets",
            'staterrsimulfit':
            'unable to obtain or estimate statistical errors for all data sets',
            'shape': "data set '%s' does not specify a shape",
            'nogrouping': "data set '%s' does not specify grouping flags",
            'noquality': "data set '%s' does not specify quality flags",
            'nostaterr': "data set '%s' does not specify statistical errors",
            'nosyserr': "data set '%s' does not specify systematic errors",
            'groupset': '%s %s grouping flag is already %s',
            'subtractset': '%s %s subtract flag is already %s',
            'nobkg': "data set '%s' does not have any associated backgrounds",
            'noarf': "data set '%s' does not have an associated ARF",
            'bad': "unknown %s: '%s'",
            'badchoices': "unknown %s: '%s'\nValid options: %s",
            'idsnotarray': "%s ids '%s' does not appear to be an array",
            'badids': "%s is not a valid %s id in %s",
            'noenergybins': "%s does not specify energy bins",
            'invalidchannel': 'invalid channel number: %s',
            'energytochannel': 'Unable to map energy bin %s to channel number',

            'subtractlength': ("subtract can only be used when the input" +
                               " source and background data sets are the" +
                               " same length"),
            'incompleteresp': ('response incomplete for dataset %s,' +
                               ' check the instrument model'),
            'incompatibleresp': "RMF '%s' is incompatible with PHA dataset '%s'",
            'nocoord': "data set '%s' does not contain a %s coordinate system",
            'bkgmodel': 'background %r has no associated model',
            'plottype': "unknown plot type '%s', choose %s",
            'normffake': 'An RMF has not been found or supplied for data set %s',
            'noenerg': 'no energy grid found in PHA response',
            'norsp': 'No instrument response found for dataset %s',
            'nobrsp': 'No instrument response found for dataset %s background %s',
            'ogip-error': "The %s '%s' %s"
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, DataErr.dict, key, *args)


class PSFErr(SherpaErr):

    dict = {'notstr': 'PSF model parameters cannot be strings',
            'nofold': 'PSF model has not been folded',
            'notset': 'PSF kernel has not been set',
            'nopsf': "model '%s' does not have an associated PSF function",
            'mismatch': 'array size mismatch between %s and %s',
            'mismatch_dims': "kernel '%s' and data '%s' do not match: %dD vs %dD",
            'badsize': 'PSF kernel size must be <=  data size, kernel: %s data: %s',
            'ndim': 'PSF model dimension must be <= 2'
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, PSFErr.dict, key, *args)


class InstrumentErr(SherpaErr):

    dict = {'baddata': "data set '%s' %s and cannot be used with a pileup model",
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, InstrumentErr.dict, key, *args)


class IOErr(SherpaErr):

    dict = {'openfailed': '%s',
            'setcolfailed': 'setting column %s has failed with %s',
            'nokeyword': "file '%s' does not have a '%s' keyword",
            'filefound': "file '%s' exists and clobber is not set",
            'filenotfound': "file '%s' not found",
            'badfile': "'%s' is not a filename or %s",
            'noarrs': 'No input array(s) found',
            'noarrayswrite': 'please supply array(s) to write to file',
            'noarrays': 'no arrays found to be loaded',
            'arraysnoteq': 'not all arrays are of equal length',
            'badargs': "data set '%s' takes at least %s args",
            'badarray': "'%s' must be a Numpy array, list, or tuple",
            'notascii': "file '%s' does not appear to be ASCII",
            'noparamcols': 'No parameter columns found in %s',
            'wrongnumcols': "Expected %d columns but found %d",
            'reqcol': "Required column '%s' not found in %s",
            'badcol': "unable to read required column %s from file",
            'badimg': 'unable to read required image %s from file',
            'notimage': "data set '%s' does not contain an image",
            'notpha': "data set '%s' does not contain a PHA spectrum",
            'onecolneedtwo': "Only found 1 column in %s, need at least 2",
            'writenoimg': "writing images in ASCII is not supported",
            'badext': "file '%s' does not contain a binary table extension",
            'notrsp': "file '%s' does not appear to be %s",
            'bad': 'unknown %s: %s',
            'npconv1d': "numpy_convolution for 1D only",
            'start<stop': "start < stop, where start=%s stop=%s",
            'step>0': "step > 0, where step=%s",
            'nobins': "not enough bins, start=%s stop=%s step=%s",
            'dimarr': 'dim must be an array of dimensions',
            'dimdatasp2d': 'dimensions for dataspace2d must be > 1',
            'boundscheck': 'the energy range is not consistent, %g !< %g',
            '>axes': 'Found more than %s axes',
            'z<=0': 'at least one redshift was <= 0',
            'erange': 'the energy range(s) are not 0 <= elo < ehi',
            'energoverlap': 'the model is defined over the energy range %f to %f keV which does not fully overlap the %s energy band of %f to %f keV %s',
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, IOErr.dict, key, *args)


class ModelErr(SherpaErr):
    "Error in creating or using a model"

    dict = {'numthawed': "expected %d thawed parameters, got %d",
            'badinstance': '%s instance cannot be created from another model',
            'noncall': 'attempted to create %s from non-callable object of type %s',
            'needsint': 'A non-overlapping integrated grid is required for model evaluation,\ne.g. [0.1,0.2],[0.2,0.3]',
            'alwaysint': '%s model is defined as integrated',
            'filterarray': "filter '%s' is not an array",
            'filtermismatch': "Mismatch between %s and %s",
            'nobkg': 'background model %s for data set %s has not been set',
            'nogrid': 'There is no grid on which to evaluate the model',
            'needspoint': 'A non-integrated grid is required for model evaluation',
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, ModelErr.dict, key, *args)


class PlotErr(SherpaErr):
    "Error in creating or using a plotting class"

    dict = {'nodataormodel': 'data or model plot is missing',
            'ordercolors': "orders list length '%s' does not match colors list length '%s'",
            'notorder': "'%i' is not a valid order",
            'orderarrfail': "mismatch between model orders and response ids",
            'plotfac': "%s plot not defined for PHA plotting factor '%i'",
            'wrongtype': "Plot type '%s' not found in %s",
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, PlotErr.dict, key, *args)


class StatErr(SherpaErr):

    dict = {'nostat': "User statistic '%s' has no %s function",
            'badstat': '%s not applicable using current statistic: %s',
            'chi2noerr': 'If you select chi2 as the statistic, all datasets must provide a staterror column',
            'usecstat': 'No background data has been supplied. Use cstat',
            'mismatch': 'size mismatch between %s (%d) and %s (%d)',
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, StatErr.dict, key, *args)


class DS9Err(SherpaErr):

    dict = {'open': "Imager not open",
            'delframe': "Could not delete frames",
            'retreg': "Could not return region in CIAO format",
            'newframe': "Could not create new frame",
            'settile': "Could not set tile option",
            'noimage': "Could not display image",
            'setwcs': "Could not replace WCS",
            'setreg': "Could not set CIAO region format",
            'badreg': "Could not use %s as a region or region file",
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, DS9Err.dict, key, *args)


class ConfidenceErr(SherpaErr):

    dict = {'badarg': "%s must be %s",
            'badlimits': 'Bad parameter limits',
            'needlist': 'Please provide a list of %s',
            'frozen': 'Frozen parameter %s cannot be used for %s',
            'thawed': 'Thawed parameter %s not found in %s',
            'badargconf': '%s is inappropriate for confidence limit estimation',
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, ConfidenceErr.dict, key, *args)


class SessionErr(SherpaErr):

    dict = {'nofit': 'no %s has been performed',
            'noaction': "%s has not been performed",
            }

    def __init__(self, key, *args):
        SherpaErr.__init__(self, SessionErr.dict, key, *args)
