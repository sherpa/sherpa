# 
#  Copyright (C) 2008  Smithsonian Astrophysical Observatory
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
