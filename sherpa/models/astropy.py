try:

    from astropy.modeling.core import _ModelMeta
    from astropy.modeling import models as apy_models
    from astropy.modeling.models import Gaussian1D as _gaussian1d
    HAS_APY=True
except ImportError:
    HAS_APY=False

from sherpa import ui
from sherpa.models import Parameter
from sherpa.models import ArithmeticModel, modelCacher1d


def make_class(model=_gaussian1d):
    class AstropyModelSherpa(ArithmeticModel):
        def __init__(self, *args, **kwargs):
            try:
                self._astropy_model = model(*args,**kwargs)
                if self._astropy_model.name is None:
                  self.name = "apy_%s" % self._astropy_model.__class__.name
                else:
                  self.name = "apy_%s" % self._astropy_model.name
                  self.__class__ = "apy_%s" % self._astropy_model.__class__.name

                for pname in self._astropy_model._param_metrics.keys():
                    setattr(self,pname,Parameter(self.name, pname,))
                ArithmeticModel.__init__(self,self.name,pars=[getattr(self,pn) for pn in self._astropy_model._param_metrics.keys()])
            except:  # This isn't really needed when using object interface
                pass  # This is due to models which require input

        def calc(self, p, *args, **kwargs):
            kwargs['integrate']=bool_cast(self.integrate)
            if len(args)==2:
                return self._astropy_model(args[0], args[1], *p, **kwargs)
            return self._astropy_model(args[0], *p, **kwargs)

    return AstropyModelSherpa
# proabably need to add startup/tardown or/and guess

if HAS_APY:
    #Test that the init of a model uses the correct model and not just the first!
    for apy_mod in _ModelMeta.registry:
        if hasattr(apy_models,apy_mod.name):
            try:
                exec("make_class(apy_models.{0})()".format(apy_mod.name))
            except TypeError:
                pass

    # this loads the models from astropy.modeling.models into sherpa
    for apy_mod in _ModelMeta.registry:
        if hasattr(apy_models,apy_mod.name):
            try:
                exec("apy_{0} = make_class(apy_models.{0})".format(apy_mod.name))
                exec('ui.load_user_model(apy_{0},"apy_{0}")'.format(apy_mod.name))
                exec('ui.add_user_pars("apy_{0}", apy_mod.param_names, apy_mod.parameters)'.format(apy_mod.name))
            except TypeError:
                pass
