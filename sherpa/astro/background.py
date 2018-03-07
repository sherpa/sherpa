#
#  Copyright (C) 2007, 2018  Smithsonian Astrophysical Observatory
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

from sherpa.models.model import CompositeModel, ArithmeticModel
from sherpa.utils.err import DataErr

__all__ = ('BackgroundSumModel',)


class BackgroundSumModel(CompositeModel, ArithmeticModel):

    def __init__(self, srcdata, bkgmodels):
        self.srcdata = srcdata
        self.bkgmodels = bkgmodels
        scale_factor = self.srcdata.sum_background_data(lambda key, bkg: 1)
        bkgnames = [model.name for model in bkgmodels.values()]
        bkgstr = ' + '.join(bkgnames)

        # In complex cases scale_factor can be an array, which
        # does not make sense to include here, although if we want
        # to use the model expression to recover the data [*] then there
        # needs to be a way to encode this here. For now, punt
        # and use a "token" - in this case 'bkgscale' - when
        # a scalar value does not work.
        #
        # [*] note that the semantics of the string version here are
        #     unclear; is it purely as a 'decorative' string for the user,
        #     or is it meant to be machine-readable/actionable?
        #
        # Alternatively we could use
        #   name = '{} * ({})'.format(scale_factor, bkgstr)
        #
        # which displays the contents of the array but is not very
        # user-friendly.
        #
        try:
            name = '%g * (' % scale_factor + bkgstr + ')'
        except TypeError:
            name = 'bkgscale * ({})'.format(bkgstr)

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
