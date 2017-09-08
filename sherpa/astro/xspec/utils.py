#
#  Copyright (C) 2017  Smithsonian Astrophysical Observatory
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
from . import _xspec

__all__ = ['ModelFunction', 'include_if', 'version_at_least']

NOT_COMPILED_FUNCTION_MESSAGE = "Calling an xspec function that was not compiled"
DISABLED_MODEL_MESSAGE = "Model {} is disabled because of an unmet condition"


class ModelFunction(object):
    """
    This class for xspec model function strings. The string is interpreted as an xspec function
    during the creation of the xspec model python class. This wrapper ensures that Sherpa graciously informs
    the user when they are trying to use an xspec model that is not available in the version they have installed.
    """
    def __init__(self, function_name):
        try:
            self.function = getattr(_xspec, function_name)
        except AttributeError:
            self.function = None

    def __call__(self, *args, **kwargs):
        """
        There may be edge cases where the model meets the condition expressed in the decorator, but the model is not
        included in the xspec build. We need to handle this error case before we just delegate the
        """
        if not self.function:
            raise AttributeError(NOT_COMPILED_FUNCTION_MESSAGE)

        return self.function(*args, **kwargs)


class include_if(object):
    """
    Generic decorator for including xspec models conditionally. It takes a boolean condition as an argument.
    If the boolean condition is not met, then the model is not included, and its function is replaced with a
    dummy function that throws an exception.

    If the model is disabled, then its class's `version_enabled` attribute is set to `False`.
    """
    def __init__(self, condition):
        self.condition = condition

    def __call__(self, model_class):
        def throw(*args, **kwargs):
            raise AttributeError(DISABLED_MODEL_MESSAGE.format(model_class.__name__))

        if not self.condition:
            model_class._calc = throw
            model_class.version_enabled = False

        return model_class


class version_at_least(include_if):
    """
    Decorator which takes a version string as an argument and enables a model only if
    the xspec version detected at runtime is equal or greater than the one provided to the decorator.
    """
    def __init__(self, version_string):
        include_if.__init__(self, _xspec.get_xsversion() >= version_string)
