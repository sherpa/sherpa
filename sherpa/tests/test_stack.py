#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from sherpa.utils import SherpaTestCase
import stk

import os
_this_dir = os.path.dirname(__file__)

class test_stack(SherpaTestCase):

    def test_build_stack(self):
        def get_name(name):
            return '/'.join((_this_dir, name))

        out = stk.build('@+{}/{}'.format(_this_dir, 'a.lis'))
        self.assertEquals([get_name('a'), get_name('a1'), get_name('a2'),
                           get_name('b'), get_name('b1'), get_name('b2')], out)

