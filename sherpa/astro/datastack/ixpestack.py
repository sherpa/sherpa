#
# Copyright (C) 2022
# MIT
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
import re

import sherpa
from sherpa.utils import public
from sherpa.utils.logging import config_logger
from sherpa.astro.datastack.ds import DataStack
from sherpa.astro import ui
from sherpa.astro.models.polarization import PolarizationDependence

logger = config_logger(__name__)

ID_STR = '__ID'

@public
def model_wrapper(func):
    def wrapfunc(self, model):
        """Run a model-setting function for each of the datasets."""
        datasets = self.filter_datasets()
        try:
            # if model is passed as a string
            model = eval(model, globals(), sherpa.astro.ui.__dict__)
        except TypeError:
            pass
        except Exception as exc:
            raise type(exc)('Error converting model "{0}" to a sherpa '
                            'model object: {1}'.format(model, exc))
        for i, dataset in enumerate(datasets):
            id_ = dataset['id']
            logger.info('Setting stack model using {0}() for id={1}'.format(
                func.__name__, id_))
            print('Setting stack model using {0}() for id={1}'.format(
                func.__name__, id_))
            new_model, model_comps = create_stack_model(model, id_)
            
            for comp in model_comps:
                compmodel = model_comps[comp]['model']
                if isinstance(compmodel, PolarizationDependence):
                    if 'XFLT0001' in dataset['data'].header.keys():
                        compmodel.stokes.val = float(dataset['data'].header['XFLT0001'].split(':')[1])
                    else:
                        logger.info(f'No XFLT0001 keyword in dataset {id_}. Assuming Stokes I.')
                        print(f'No XFLT0001 keyword in dataset {id_}. Assuming Stokes I.')
                    # Now link the parameters (except for "stokes" which is by construction the last one
                    # on the list) to the model of the first instance
                    if i != 0:
                        for phere, p0 in zip(compmodel.pars[:-1],
                                         ds.datasets[0]['model_comps'][comp]['model'].pars[:-1]):
                            phere.link = p0
                    
            func(id_, new_model)
            dataset['model_comps'].update(model_comps)
        return None

    wrapfunc.__name__ = func.__name__
    wrapfunc.__doc__ = func.__doc__
    return wrapfunc


@public
def create_stack_model(model, id_):
    """Create any model components for the given data set.
    Parameters
    ----------
    model
       The model expression to use (a string or a Sherpa model
       expression). The names of the model instances will be
       used to create the dataset-specific names, replacing any
       string matching the template with the data set identifier.
    id_ : string or integer
       The identifier for the data set.
    Returns
    -------
    instance, comps
       The instance for the model expression and a dict containing
       any model components.
    See Also
    --------
    set_template_id
    """
    return _create_stack_model(model, id_)


def _create_stack_model(model, model_id, model_components=None):
    model_comps = model_components or {}

    if hasattr(model, 'parts'):
        # Recursively descend through model and create new parts
        # (as needed) corresponding to the stacked model components.
        newparts = []
        for part in model.parts:
            newpart, new_comps = _create_stack_model(
                part, model_id, model_comps)
            model_comps.update(new_comps)
            newparts.append(newpart)

        if hasattr(model, 'op'):
            model = model.op(*newparts)
        elif hasattr(model, 'rmf') and hasattr(model, 'arf'):
            model = sherpa.astro.instrument.RSPModelPHA(rmf=model.rmf,
                                                        model=newparts[0],
                                                        arf=model.arf,
                                                        pha=model.pha)
        elif hasattr(model, 'rmf'):
            model = sherpa.astro.instrument.RMFModelPHA(rmf=model.rmf,
                                                        model=newparts[0],
                                                        pha=model.pha)
        elif hasattr(model, 'arf'):
            model = sherpa.astro.instrument.ARFModelPHA(model=newparts[0],
                                                        arf=model.arf,
                                                        pha=model.pha)
        else:
            raise ValueError(
                "Unexpected composite model {0} (not operator, ARF or RMF)".format(repr(model)))
    else:
        if hasattr(model, 'val'):
            model = model.val
        else:
            try:
                model_type, model_name_ID = model.name.split('.')
            except ValueError:
                raise ValueError(
                    'Model name "{0}" must be in format <model_type>.<name>'.format(model.name))
            if isinstance(model, PolarizationDependence) and (ID_STR not in model_name_ID):
                model_name_ID = model_name_ID + ID_STR
            model_name = re.sub(ID_STR, str(model_id), model_name_ID)
            try:
                model_temp = getattr(getattr(sherpa.astro.ui, model_type),
                                model_name)
            except AttributeError:
                # Must be a user model, so use add_model to put a
                # modelwrapper function into namespace
                sherpa.astro.ui.add_model(model.__class__)
                model_temp = sherpa.astro.ui.create_model_component(
                            model_type, model_name)
                    
            # In these two cases, models shall be unique for each dataset,
            # so we pass on the new instance that we just created.
            if (ID_STR in model_name_ID) or isinstance(model_temp, PolarizationDependence):
                model = model_temp

            model_name_no_ID = re.sub(ID_STR, "", model_name_ID)
            model_comps[model_name_no_ID] = dict(model=model,
                                                 model_name=model_name)

    return model, model_comps


class IXPEStack(DataStack):
            
    def load_filelist(self, id, arg=None, use_errors=True):
        if arg is None:
            id, arg = arg, id

        if id is not None:
            if self._default_instance:
                ui.load_pha(id, arg, use_errors)
                return
            else:
                raise AttributeError(load_error_msg(id))
        try:
            for infile in arg:
                self._load_func(ui.load_pha, infile, use_errors)
        except (NameError, OSError, IOErr):
            self._load_func(ui.load_pha, arg, use_errors)
           
            
    def show_stack(self):
        """Display basic information about the contents of the data stack.
        A brief summary of each data set in the stack is printed to the
        standard output channel. The information displayed depends on the
        type of each data set.
        """
        cols, vallist = self._stack_repr_list()
        print(('{:<8} ' * len(l)).format(*cols))
        for l in vallist:
            print(('{:<8} ' * len(l)).format(*l))
            
    def _stack_repr_list(self):
        cols = ['id', 'name', 'OBJECT', 'OBS_ID', 'TELESCOP', 'INSTRUME', 'LIVETIME', 'MJD-OBS', 'XFLT0001']
        vallist = []
        for dataset in self.filter_datasets():
            values = [dataset['id'], dataset['data'].name]
            hdr = dataset['data'].header if hasattr(dataset['data'], 'header') else {}
            for k in cols[2:]:
                values.append(hdr.get(k, 'N/A'))
            vallist.append(values)
        return cols, vallist
    
    def _repr_html_(self):
        cols, vallist = self._stack_repr_list()
        out = '<table style="width:100%">'
        out = out + '<tr>' + ('<th>{}</th>' * len(cols)).format(*cols) + '</tr>'
        for l in vallist:
            out = out + '<tr>' + ('<td>{}</td>' * len(l)).format(*l) + '</tr>'
        out = out + '</table>'
        return out

    set_source = model_wrapper(ui.set_source)
    set_model = model_wrapper(ui.set_model)
    set_bkg_model = model_wrapper(ui.set_bkg_model)
    set_full_model = model_wrapper(ui.set_full_model)
    set_bkg_full_model = model_wrapper(ui.set_bkg_full_model)
