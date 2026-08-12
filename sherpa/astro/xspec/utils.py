#
#  Copyright (C) 2017, 2018, 2019, 2022, 2026
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

import re

try:
    from . import _xspec
except ImportError as ie:
    # It would be nicer for the user to say "from None" rather than
    # "from ie" here but it may lose important information (in the
    # case there is some problem initializing the XSPEC library,
    # rather than it just not being built).
    #
    raise ImportError("XSPEC support is not enabled") from ie

from sherpa.astro.utils.xspec import get_version

__all__ = ['ModelMeta', 'include_if', 'version_at_least']


XSPEC_VERSION = get_version(_xspec.get_xsversion())


class ModelMeta(type):
    """Metaclass for XSPEC models.

    If the class has no _calc method then automatically select it from the
    sherpa.astro.xspec._xspec module, selecting the symbol matching the
    _xspec_name attribute.

    .. versionchanged:: 4.19.0
       The library routines are now named after the model name, not the
       function name, which has simplified the logic for this metaclass.

    """
    NOT_COMPILED_FUNCTION_MESSAGE = "Calling an xspec function that was not compiled"

    def __init__(cls, *args, **kwargs):

        # If the _module is set to None then the class is assumed
        # to have its own calc/_calc behaviour.
        #
        if cls._module is not None:
            # Since this is the class we can not use cls.xspec_name as
            # it is a property.
            funcname = cls._xspec_name

            try:
                cls._calc = getattr(cls._module, funcname)
            except (AttributeError, TypeError):
                # The assumption is that this model class is newer
                # than the XSPEC library, or it is a class that the
                # user is not expected to use (e.g. XSModel).
                #
                cls._calc = ModelMeta._not_compiled

        super(ModelMeta, cls).__init__(*args, **kwargs)

    @staticmethod
    def _not_compiled(*args, **kwargs):
        raise AttributeError(ModelMeta.NOT_COMPILED_FUNCTION_MESSAGE)


def equal_or_greater_than(version_string):
    """Compare with the version of the current xspec instance.

    Parameters
    ----------
    version_string : str
        The version to compare against the current XSPEC version. It
        can include an XSPEC patch level, such as "12.12.0c", but the
        patch level is ignored.

    Returns
    -------
    flag : bool
        `True` if the version of xspec is equal or greater than the
        argument, `False` otherwise

    Notes
    -----
    For better or worse the xspec current instance is not cached
    across calls. It probably could be but it just seems safer not to,
    and any overhead insists on models initialization only.

    """
    return XSPEC_VERSION >= get_version(version_string)


class include_if():
    """
    Generic decorator for including xspec models conditionally. It takes a boolean condition as an argument.
    If the boolean condition is not met, then the model is not included, and its function is replaced with a
    dummy function that throws an exception.

    If the model is disabled, then its class's `version_enabled` attribute is set to `False`.
    """
    DISABLED_MODEL_MESSAGE = "Model {} is disabled because of an unmet condition"

    def __init__(self, condition):
        self.condition = condition

    def __call__(self, model_class):
        if not self.condition:
            model_class.version_enabled = False
            model_class._calc = self._disabled(self.get_message(model_class))

        return model_class

    def get_message(self, model_class):
        return self.DISABLED_MODEL_MESSAGE.format(model_class.__name__)

    @staticmethod
    def _disabled(message):
        def wrapped(*args, **kwargs):
            raise AttributeError(message)

        return wrapped


class version_at_least(include_if):
    """
    Decorator which takes a version string as an argument and enables a model only if
    the xspec version detected at runtime is equal or greater than the one provided to the decorator.
    """
    DISABLED_MODEL_MESSAGE = "Model {} is disabled because XSPEC version >= {} is required"

    def __init__(self, version_string):
        include_if.__init__(self, equal_or_greater_than(version_string))
        self.version_string = version_string

    def get_message(self, model_class):
        return self.DISABLED_MODEL_MESSAGE.format(model_class.__name__, self.version_string)
