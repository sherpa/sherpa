#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
import sherpa.utils.integration as integration
from sherpa.utils import SherpaTestCase


class test_integration(SherpaTestCase):

    def test_c_api(self):
        self.assert_(hasattr(integration, '_C_API'))
        self.assertEqual(type(integration._C_API).__name__, 'PyCObject')
