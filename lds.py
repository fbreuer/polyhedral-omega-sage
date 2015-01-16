import geometry

from sage.all_cmdline import *

from util import sageify

class LinearDiophantineSystem(object):

    def __init__(self,A,b,E=None):
        """Creates a new LinearDiophantineSystem Ax >=^E b, x >= 0 for a given matrix A, a right-hand side b and a list E of flags indicating whether the corresponding row is an inequality or an equation. 

        The class LinearDiophantineSystem can be used to compute symbolic cone and/or rational function representations of the set of all non-negative integral solutions to the given linear system. The parameters A,b,E cannot be modified once the LinearDiophantineSystem has been constructed. The methods for solving the system, symbolic_ones, fundamental_parallelepipeds, and rational_function_string, each cache the results of their computations so that they can be run several times and the computation is done only once.

        Arguments:
            A: An integer matrix represented as a sequence of rows. Entries can be ints or Sage Integers.
            b: An integer vector, represented as a sequence. Entries can be ints or Sage Integers.
            E: A list with 0-1 entries which encodes with rows are equalities. (1 means equality, 0 means inequality.) If E is None, then all constraints are assumed to be inequalities, i.e., E is constant 0.
        """
        if E == None:
            E = [ 0 for i in xrange(len(A))]
        self._A = A
        self._b = b
        self._E = E
        self._cone = geometry.SymbolicCone.macmahon_cone(A,b)
        self._symbolic_cones = None
        self._fundamental_parallelepipeds = None
        self._rational_function = None

    def symbolic_cones(self,logging=False):
        """Computes a representation of the set of all non-negative integral solutions as a linear combination of symbolic cones.

        Returns a CombinationOfCones representing the result.
        """
        if self._symbolic_cones == None:
            E = list(self._E)
            E.reverse() # reverse, because we want to eliminate the last rows
            self._symbolic_cones = geometry.CombinationOfCones({self._cone: Integer(1)})
            i = 0
            for e in E:
                self._symbolic_cones = self._symbolic_cones.map_cones(lambda c: c.eliminate_last_coordinate(equality=e))
                if logging:
                    i = i + 1
                    d = len(self._A[0])
                    max_cones = binomial(d+i,d)
                    num_cones =len(self._symbolic_cones)
                    print "iteration %d of %d: %d cones (%.2f%% of theoretical maximum %d) " % (i,len(E),num_cones,float(num_cones/max_cones)*100,max_cones)
        return self._symbolic_cones

    def fundamental_parallelepipeds(self):
        """For each of the cones in the solution returned by symbolic_cones, this method enumerates all lattice points in the corresponding fundamental_parallelepipeds.

        Returns a dictionary, mapping cones to the set of lattice points in the fundamental_parallelepipeds. Each set of lattice point is represented as a list of tuples of Sage Integers.
        """
        if self._fundamental_parallelepipeds == None:
            cones = self.symbolic_cones()
            self._fundamental_parallelepipeds = {}
            for cone in cones.keys():
                points = cone.enumerate_fundamental_parallelepiped()
                self._fundamental_parallelepipeds[cone] = points
        return self._fundamental_parallelepipeds

    def rational_function_string(self):
        """Computes the symbolic_cones and the fundamental_parallelepipeds and then turns this data into a string representation of the corresponding rational function.

        Returns the rational function as a string of the form "r1 + r2 + ... + rn" where each ri is a rational function of the form "m * (p) / (q)" where p is a Laurent polynomial expanded with respect to the monomial basis and q is a polynomial written as a product of factors "(1-m)" where m is a multivariate monomial. Variables are "z1", ..., "zn" and exponentiation is written "z1**3", for example.
        """
        if self._rational_function == None:
            def stringify_cone(cone,fundamental_parallelepiped):
                d = cone.ambient_dimension
                V = cone.generators
                def stringify_monomial(z):
                    return "*".join(("z%d**%d" % (i,z[i]) for i in xrange(d)))
                num = "+".join((stringify_monomial(t) for t in fundamental_parallelepiped))
                den = "*".join(( "(1-%s)" % stringify_monomial(c) for c in V))
                return ("(%s)/(%s)" % (num,den))
            rational_functions = []
            for cone, multiplicity in self.symbolic_cones().iteritems():
                pi = self.fundamental_parallelepipeds()[cone]
                rational_functions.append( str(multiplicity) + "*" + stringify_cone(cone,pi) )
            if len(rational_functions) > 0:
                self._rational_function = "+".join(rational_functions)
            else:
                self._rational_function = "0"
        return self._rational_function