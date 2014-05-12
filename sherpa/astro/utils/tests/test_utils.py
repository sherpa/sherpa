#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2008)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from sherpa.utils import SherpaTestCase
from sherpa.astro.utils import is_in

class test_utils(SherpaTestCase):
    
    def setUp(self):
           self.long  = [100,249,400,450,500,601,1024]
           self.short = [100,249,601,1024]

    def test_response_filter_logic(self):

        # outside case
        self.assert_( is_in(self.long, 50, 2400) )

        # lo case
        self.assert_( is_in(self.long, 100, 200) )

        # hi case
        self.assert_( is_in(self.long, 50, 1024) ) 

        # 'hidden' lo case
        self.assert_( is_in(self.long, 250, 2000) )

        # 'hidden' hi case
        self.assert_( is_in(self.long, 50, 250) )

        # 'hidden' interval case w/ noticed channels inside
        self.assert_( is_in(self.long, 250, 600) )

        # 'hidden' interval case w/ *no* noticed channels inside
        self.assert_( not is_in(self.short, 250, 600) )
