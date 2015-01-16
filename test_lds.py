import unittest
from sage.all_cmdline import *

from lds import *
from util import *
from geometry import *

class TestLinearDiophantineSystem(unittest.TestCase):
    def test_simple_2d_ineq(self):
        A = ( (1,-1) )
        b = ( 0, )
        E = [ 0 ]
        lds = LinearDiophantineSystem(A,b,E)
        cones = lds.symbolic_cones()
        V = sageify( ((1,0), (1,1)) )
        q = sageify( (0,0) )
        o = sageify( (0,0) )
        c = SymbolicCone(V,q,o)
        self.assertEqual(cones,CombinationOfCones({c:Integer(1)}))
        pis = lds.fundamental_parallelepipeds()
        self.assertDictEqual(pis,{c: sageify([ (0,0)  ])})
        r = lds.rational_function_string()
        self.assertEqual(r,"1*(z1**0*z2**0)/((1-z1**1*z2**0)*(1-z1**1*z2**1))")

if __name__ == '__main__':
    unittest.main()