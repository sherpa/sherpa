#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
from sherpa.models.model import CompositeModel, ArithmeticModel
from sherpa.utils.err import DataErr

__all__ = ('BackgroundSumModel',)


class BackgroundSumModel(CompositeModel, ArithmeticModel):

    def __init__(self, srcdata, bkgmodels):
        self.srcdata = srcdata
        self.bkgmodels = bkgmodels
        scale_factor = self.srcdata.sum_background_data(lambda key, bkg:1)
        bkgnames = [model.name for model in bkgmodels.values()]
        name = '%g * (' % scale_factor + ' + '.join(bkgnames) + ')'
        CompositeModel.__init__(self, name, self.bkgmodels.values())

    def calc(self, p, *args, **kwargs):
        def eval_bkg_model(key, bkg):
            bmodel = self.bkgmodels.get(key)
            if bmodel is None:
                raise DataErr('bkgmodel', key)
            # FIXME: we're not using p here (and therefore assuming that the
            # parameter values have already been updated to match the contents
            # of p)
            return bmodel(*args, **kwargs)

        return self.srcdata.sum_background_data(eval_bkg_model)            
