# 
#  Copyright (C) 2010  Smithsonian Astrophysical Observatory
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

import numpy
from sherpa.utils._pykdtree import *
#from sherpa.utils.kdtree import *
from sherpa.utils import SherpaTestCase, SherpaTest


class test_kdtree(SherpaTestCase):


    def test_empty(self):
        nn = KDTree_3Double()
        self.assertEqual(0, nn.size())
        
        actual = nn.find_nearest((2,3,0))
        self.assertTrue(None==actual, "%s != %s"%(str(None), str(actual)))



    def test_get_all(self):
        nn = KDTree_3Double()
        o1 = object()
        nn.add(((1,1,0), id(o1)))
        o2 = object()
        nn.add(((10,10,0), id(o2)))
        o3 = object()
        nn.add(((11,11,0), id(o3)))

        self.assertEqual([((1,1,0), id(o1)), ((10,10,0), id(o2)), ((11,11,0), id(o3))], nn.get_all())
        self.assertEqual(3, nn.size())
        
        nn.remove(((10,10,0), id(o2)))
        self.assertEqual(2, nn.size())        
        self.assertEqual([((1,1,0), id(o1)), ((11,11,0), id(o3))], nn.get_all())



    def test_nearest(self):
        nn = KDTree_3Double()

        nn_id = {}
        
        o1 = object()
        nn.add(((1,1,0), id(o1)))
        nn_id[id(o1)] = o1
        o2 = object()
        nn.add(((10,10,0), id(o2)))
        nn_id[id(o2)] = o2
        o3 = object()
        nn.add(((4.1, 4.1,0), id(o3)))
        nn_id[id(o3)] = o3
        
        expected =  o3
        actual = nn.find_nearest((2.9,2.9,0))[1]
        self.assertTrue(expected==nn_id[actual], "%s != %s"%(str(expected), str(nn_id[actual])))

        expected = o3
        actual = nn.find_nearest((6, 6,0))[1]
        self.assertTrue(expected==nn_id[actual], "%s != %s"%(str(expected), str(nn_id[actual])))

        expected = o2
        actual = nn.find_nearest((8, 9,0))[1]
        self.assertTrue(expected==nn_id[actual], "%s != %s"%(str(expected), str(nn_id[actual])))


    def test_remove(self):
        class C:
            def __init__(self, i):
                self.i = i
                self.next = None

        nn = KDTree_3Double()

        k1, o1 = (1.1,1.1,0), C(7)
        self.assertFalse(nn.remove((k1, id(o1))), "This cannot be removed!")
        nn.add((k1, id(o1)))

        k2, o2 = (1.1,1.1,0), C(7)
        nn.add((k2, id(o2)))

        self.assertEqual(2, nn.size())
        self.assertTrue(nn.remove((k2, id(o2))))
        self.assertEqual(1, nn.size())        
        self.assertFalse(nn.remove((k2, id(o2))))        
        self.assertEqual(1, nn.size())

        nearest = nn.find_nearest(k1)
        self.assertTrue(nearest[1] == id(o1), "%s != %s"%(nearest[1], o1))



    def test_count_within_range(self):
        nn = KDTree_2Double()

        items = [(i,j) for i in range(3) for j in range(3)]

        for p in items:
            nn.add((p, id(p)))


        expected = 1
        res = nn.count_within_range((0,0), 0.0)
        self.assertEqual(expected, res, "Counted %i points instead of %i"%(res, expected))

        expected = 4
        res = nn.count_within_range((0,0), 1.0)
        self.assertEqual(expected, res, "Counted %i points instead of %i"%(res, expected))
        
        expected = 9
        res = nn.count_within_range((0,0), 2.0)
        self.assertEqual(expected, res, "Counted %i points instead of %i"%(res, expected))


    def test_find_within_range(self):
        nn = KDTree_6Double()
 
        nn_id = {}
        
        o1 = object()
        nn.add(((1,1,0,0,0,0), id(o1)))
        nn_id[id(o1)] = o1
        o2 = object()
        nn.add(((10,10,0,0,0,0), id(o2)))
        nn_id[id(o2)] = o2
        o3 = object()
        nn.add(((4.1, 4.1,0,0,0,0), id(o3)))
        nn_id[id(o3)] = o3
        
        expected =  set([long(id(o1)), long(id(o3))])

        actual = set([ident
                      for _coord, ident
                      in nn.find_within_range((2.1,2.1,0,0,0,0), 3.9)])
        self.assertTrue(expected==actual, "%s != %s"%(str(expected), str(actual)))


    def test_find_exact(self):

        nn = KDTree_3Double()
        vector = [(float(i),float(j),float(k)) for i in range(2) for j in range(2) for k in range(2)]

        for p in vector:
            nn.add((p, id(p)))

        expected = set((tuple(vector[1]),))

        actual = set((tuple(nn.find_exact((vector[1], id(vector[1])))[0]),))

        self.assertTrue(expected==actual, "%s != %s"%(str(expected), str(actual)))


if __name__ == '__main__':

    import sherpa.utils as utils
    SherpaTest(utils).test()
