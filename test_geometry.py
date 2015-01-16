import unittest
from geometry import *
from util import *

from sage.all_cmdline import *

# define sage constants
_0 = Integer(0)
_1 = Integer(1)

class TestGeometry(unittest.TestCase):
    def test_symbolic_cone_equality(self):
        c1 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c2 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c3 = SymbolicCone(((_0,_1),(_1,_1)),(_0,_0),(_1,_0))
        c4 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        c5 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        self.assertEqual(c1,c2)
        self.assertNotEqual(c1,c3)
        self.assertNotEqual(c1,c4)
        self.assertNotEqual(c1,c5)
        self.assertNotEqual(c2,c3)
        self.assertNotEqual(c2,c4)
        self.assertNotEqual(c2,c5)
        self.assertNotEqual(c3,c4)
        self.assertNotEqual(c3,c5)
        self.assertEqual(c4,c5)

    def test_symbolic_cone_addition(self):
        c1 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c2 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c3 = SymbolicCone(((_0,_1),(_1,_1)),(_0,_0),(_1,_0))
        c4 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        c5 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        self.assertEqual(c1 + c2, CombinationOfCones({c1: Integer(2)}))
        self.assertEqual(c1 + c3, CombinationOfCones({c1: Integer(1), c3: Integer(1)}))
        self.assertEqual(c1 + c2 + c3 + c4 + c5, CombinationOfCones({c1: Integer(2), c3: Integer(1), c4: Integer(2)}))

    def test_symbolic_cone_multiplication(self):
        c1 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c2 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c3 = SymbolicCone(((_0,_1),(_1,_1)),(_0,_0),(_1,_0))
        c4 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        c5 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        self.assertEqual(c1 * 2, CombinationOfCones({c1: Integer(2)}))
        self.assertEqual(c1 * Integer(2), CombinationOfCones({c1: Integer(2)}))
        self.assertEqual(2 * c1, CombinationOfCones({c1: Integer(2)}))
        self.assertEqual(Integer(2) * c1, CombinationOfCones({c1: Integer(2)}))
        self.assertEqual(2 * c1 + c3 * 4 + (-2) * c5, CombinationOfCones({c1: Integer(2), c3: Integer(4), c4: Integer(-2)}))

    def test_symbolic_cone_flip(self):
        c1 = SymbolicCone(((Integer(1),Integer(-2)),(Integer(-1),Integer(1))), (_0,_0), (_1,_0))
        c1f = SymbolicCone(((Integer(1),Integer(-2)),(Integer(1),Integer(-1))), (_0,_0), (_1,_1))
        self.assertEqual(c1.flip(),CombinationOfCones({c1f:Integer(-1)}))
        self.assertEqual(c1f.flip(),CombinationOfCones({c1f:Integer(1)}))
        c2 = SymbolicCone(((Integer(0),Integer(-1),Integer(1)),(Integer(-0),Integer(2),Integer(-1)),(Integer(0),Integer(-2),Integer(-3))), (_0,_0,_0), (_1,_0,_1))
        c2f = SymbolicCone(((Integer(0),Integer(1),Integer(-1)),(Integer(0),Integer(2),Integer(-1)),(Integer(0),Integer(2),Integer(3))), (_0,_0,_0), (_0,_0,_0))
        self.assertEqual(c2.flip(),CombinationOfCones({c2f:Integer(1)}))
        self.assertEqual(c2f.flip(),CombinationOfCones({c2f:Integer(1)}))

    # Tests of Eliminate Last Coordinate
    # TODO: add many more cases!

    def test_symbolic_cone_eliminate_last_coordinate_simple_intersection(self):
        V = sageify( ((1,0,1), (0,1,-1)) )
        q = sageify( (0,0,0) )
        o = sageify( (0,0) )
        c = SymbolicCone(V,q,o)
        V_r = sageify( ((1,0), (1,1)) )
        q_r = sageify( (0,0) )
        o_r = sageify( (0,0) )
        c_r = SymbolicCone(V_r,q_r,o_r)
        self.assertEqual(c.eliminate_last_coordinate(),CombinationOfCones({c_r: Integer(1)}))
        # c is not modified in the process
        self.assertEqual(c,SymbolicCone(V,q,o))

    def test_symbolic_cone_eliminate_last_coordinate_all_generators_up(self):
        V = sageify( ((1,0,1), (0,1,1)) )
        q = sageify( (0,0,0) )
        o = sageify( (0,0) )
        c = SymbolicCone(V,q,o)
        V_r = sageify( ((1,0), (0,1)) )
        q_r = sageify( (0,0) )
        o_r = sageify( (0,0) )
        c_r = SymbolicCone(V_r,q_r,o_r)
        self.assertEqual(c.eliminate_last_coordinate(),CombinationOfCones({c_r: Integer(1)}))
        # c is not modified in the process
        self.assertEqual(c,SymbolicCone(V,q,o))

    def test_symbolic_cone_eliminate_last_coordinate_simple_equality(self):
        V = sageify( ((1,0,1), (0,1,-1)) )
        q = sageify( (0,0,0) )
        o = sageify( (0,0) )
        c = SymbolicCone(V,q,o)
        V_r = sageify( ((1,1),) )
        q_r = sageify( (0,0) )
        o_r = sageify( (0,) )
        c_r = SymbolicCone(V_r,q_r,o_r)
        self.assertEqual(c.eliminate_last_coordinate(equality=True),CombinationOfCones({c_r: Integer(1)}))
        # c is not modified in the process
        self.assertEqual(c,SymbolicCone(V,q,o))

    def test_symbolic_cone_eliminate_last_coordinate_simple_intersection_two_pieces(self):
        V = sageify( ((1,0,1), (0,1,-1)) )
        q = sageify( (0,0,1) )
        o = sageify( (0,0) )
        c = SymbolicCone(V,q,o)
        V_r1 = sageify( ((1,0), (0,1)) )
        q_r1 = sageify( (0,0) )
        o_r1 = sageify( (0,0) )
        c_r1 = SymbolicCone(V_r1,q_r1,o_r1)
        V_r2 = sageify( ((1,1), (0,1)) )
        q_r2 = sageify( (0,1) )
        o_r2 = sageify( (0,1) )
        c_r2 = SymbolicCone(V_r2,q_r2,o_r2)
        self.assertEqual(c.eliminate_last_coordinate(),CombinationOfCones({c_r1: Integer(1), c_r2: Integer(-1)}))
        # c is not modified in the process
        self.assertEqual(c,SymbolicCone(V,q,o))

    # Tests of fundamental parallelepiped enumeration
    # TODO add many more cases!

    def test_enumerate_fundamental_parallelepiped_standard_basis(self):
        V = sageify( ((1,0), (0,1)) )
        q = sageify( (0,0) )
        o = sageify( (0,0) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (0,0) ]))
        q = sageify( (0,0) )
        o = sageify( (1,0) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (1,0) ]))
        q = sageify( (0,0) )
        o = sageify( (0,1) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (0,1) ]))
        q = sageify( (0,0) )
        o = sageify( (1,1) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (1,1) ]))
        q = sageify( ( Integer(1)/Integer(2), 0) )
        o = sageify( (0,0) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (1,0) ]))
        q = sageify( ( Integer(-1)/Integer(3), Integer(1)/Integer(2) ) )
        o = sageify( (0,0) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (0,1) ]))

    def test_enumerate_fundamental_parallelepiped_simple_2d(self):
        V = sageify( ((4,1), (2,3)) )
        q = sageify( (0,0) )
        o = sageify( (0,0) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (0,0), (1,1), (2,1), (3,1), (2,2), (3,2), (4,2), (3,3), (4,3), (5,3) ]))
        q = sageify( (0,0) )
        o = sageify( (1,0) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (4,1), (1,1), (2,1), (3,1), (2,2), (3,2), (4,2), (3,3), (4,3), (5,3) ]))
        q = sageify( (0,0) )
        o = sageify( (0,1) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (2,3), (1,1), (2,1), (3,1), (2,2), (3,2), (4,2), (3,3), (4,3), (5,3) ]))
        q = sageify( (0,0) )
        o = sageify( (1,1) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (6,4), (1,1), (2,1), (3,1), (2,2), (3,2), (4,2), (3,3), (4,3), (5,3) ]))
        q = sageify( ( Integer(1)/Integer(2) ,0) )
        o = sageify( (0,0) )
        c = SymbolicCone(V,q,o)
        self.assertItemsEqual(c.enumerate_fundamental_parallelepiped(),sageify([ (4,1), (5,2), (2,1), (3,1), (2,2), (3,2), (4,2), (3,3), (4,3), (5,3) ]))


    def test_combination_of_cones_equality(self):
        c1 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c2 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c3 = SymbolicCone(((_0,_1),(_1,_1)),(_0,_0),(_1,_0))
        c4 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        c5 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        C1 = CombinationOfCones({c1:_1,c3:_1})
        C2 = CombinationOfCones({c2:_1,c3:_1})
        C3 = CombinationOfCones({c2:_1,c3:Integer(2)})
        C4 = CombinationOfCones({c2:-Integer(5),c4:Integer(2)})
        C5 = CombinationOfCones({c1:-Integer(5),c5:Integer(2)})
        self.assertEqual(C1,C2)
        self.assertNotEqual(C1,C3)
        self.assertNotEqual(C1,C4)
        self.assertNotEqual(C1,C5)
        self.assertNotEqual(C2,C3)
        self.assertNotEqual(C2,C4)
        self.assertNotEqual(C2,C5)
        self.assertNotEqual(C3,C4)
        self.assertNotEqual(C3,C5)
        self.assertEqual(C4,C5)

    def test_combination_of_cones_in_place_addition(self):
        c1 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c2 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c3 = SymbolicCone(((_0,_1),(_1,_1)),(_0,_0),(_1,_0))
        c4 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        c5 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        C1 = CombinationOfCones()
        self.assertEqual(C1,CombinationOfCones())
        C1 += c1
        self.assertEqual(C1,CombinationOfCones({c1: Integer(1)}))
        self.assertEqual(C1,CombinationOfCones({c2: Integer(1)}))
        C1 += c2
        self.assertEqual(C1,CombinationOfCones({c1: Integer(2)}))
        self.assertEqual(C1,CombinationOfCones({c2: Integer(2)}))
        C1 += c3
        self.assertEqual(C1,CombinationOfCones({c1: Integer(2), c3: Integer(1)}))
        self.assertEqual(C1,CombinationOfCones({c2: Integer(2), c3: Integer(1)}))
        C1 += C1
        self.assertEqual(C1,CombinationOfCones({c1: Integer(4), c3: Integer(2)}))
        self.assertEqual(C1,CombinationOfCones({c2: Integer(4), c3: Integer(2)}))
        C1 += CombinationOfCones({c4: Integer(-2), c3: Integer(-2), c1: Integer(1)})
        self.assertEqual(C1,CombinationOfCones({c1: Integer(5), c4: Integer(-2)}))
        self.assertEqual(C1,CombinationOfCones({c2: Integer(5), c4: Integer(-2)}))
        C1 += CombinationOfCones({c1: Integer(-5)})
        self.assertEqual(C1,CombinationOfCones({c4: Integer(-2)}))

    def test_combination_of_cones_in_place_multiplication(self):
        c1 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c2 = SymbolicCone(((_0,_1),(_1,_1)),(_0,_0),(_1,_0))
        c3 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        C1 = CombinationOfCones({c1: Integer(2), c2: Integer(-3), c3: Integer(1)})
        C1 *= Integer(3)
        self.assertEqual(C1,CombinationOfCones({c1: Integer(6), c2: Integer(-9), c3: Integer(3)}))
        C1 *= Integer(-2)
        self.assertEqual(C1,CombinationOfCones({c1: Integer(-12), c2: Integer(18), c3: Integer(-6)}))
        C1 *= Integer(0)
        self.assertEqual(C1,CombinationOfCones())        

    def test_combination_of_cones_map(self):
        c1 = SymbolicCone(((_0,_1),(_1,_0)),(_0,_0),(_1,_0))
        c2 = SymbolicCone(((_0,_1),(_1,_1)),(_0,_0),(_1,_0))
        c3 = SymbolicCone(((_0,_1),(_1,_1)),(_1,_0),(_0,_0))
        self.assertEqual(
            CombinationOfCones()
                .map_cones(lambda c: CombinationOfCones(c1)),
            CombinationOfCones())
        self.assertEqual(
            CombinationOfCones({c1: Integer(2)})
                .map_cones(lambda c: c1 * 3),
            CombinationOfCones({c1: Integer(6)}))
        self.assertEqual(
            CombinationOfCones({c1: Integer(2), c2: Integer(3)})
                .map_cones(lambda c: c2 * 5 + c3 if c == c1 else c),
            CombinationOfCones({c2: Integer(13), c3: Integer(2)}))


if __name__ == '__main__':
    unittest.main()