import collections
from linear_algebra import *
from util import *
from UserDict import DictMixin

from sage.all_cmdline import *

# define sage constants
_0 = Integer(0)
_1 = Integer(1)


#from sage.rings.integer import Integer

class SymbolicCone(collections.Hashable):
    """The vertex description of a simplicial symbolic cone.

    The vertex decription of a simplicial cone consists of its apex (which may be rational), the matrix of generators (which must be integral and linearly independent) and a vector of boolean flags indicating which faces are open and closed.

    SymbolicCones are immutable and hashable to facilitate using them as keys in a dictionary.

    This class implements intersecting a SymbolicCone with a coordinate half-space, enumerating the fundamental parallelepiped of a SymbolicCone, computing an inequality description of the symbolic cone and flipping a cone "forward".

    Attributes:
        generators: The generators of the cone, given as a tuple of integer vectors. Vectors are tuples of Sage integers. Generators must be linearly independent. Read-only.
        apex: The apex of the cone, given as a tuple of Sage rationals. Read-only.
        openness: A tuple of booleans, indicating which faces of the cone are open. If openness[i] == 1, the coefficient of the i-th generator is strictly bigger than one, so that the opposite face is open. Conversely, if penness[i] == 0, the coefficient of the i-th generator is strictly bigger than one, so that the opposite face is open.
        dimension: Dimension of the cone, that is, the number of generators. Read-only.
        ambient_dimension: The dimension of the ambient space in which the cone lies, that is, the length of each generator. Read-only.
    """

    def __init__(self, generators, apex, openness=None, do_canonicalize=True):
        """Initializes the SymbolicCone with the given generators, apex and openness. 

        If openness is not given (or None), then the cone is closed. Arguments have to be as specified in the class description.
        """
        # TODO: Check for TypeErrors and ValueErrors. Is the performance overhead for these checks acceptable? Also if we verify linear independence of generators?
        self._dimension = len(generators)
        self._ambient_dimension = len(generators[0])
        self._generators = generators
        self._apex = apex
        if openness == None:
            self._openness = tuple([ 0 for i in range(self.dimension)])
        else:
            self._openness = openness
        if do_canonicalize:
            self.canonicalize()

    @property
    def dimension(self):
        return self._dimension
    
    @property
    def ambient_dimension(self):
        return self._ambient_dimension
    
    @property
    def generators(self):
        return self._generators
    
    @property
    def apex(self):
        return self._apex

    @property
    def openness(self):
        return self._openness

    def __hash__(self):
        """Computes hash value of self."""
        return tuple([self._generators,self._apex,self._openness]).__hash__()

    def __eq__(self,other):
        """Two cones are equal if they have equal generators, apex and openness."""
        return self.generators == other.generators and self.apex == other.apex and self.openness == other.openness

    def __rmul__(self,other):
        """Multiplications i * c where i is an integer and c a SymbolicCone."""
        return self.__mul__(other)

    def __mul__(self,other):
        """Multiplication of SymbolicCone with an integer returns a CombinationOfCones.

        Arguments:
            other: An integer, represented as an int or a Sage Integer."""
        if type(other) == Integer:
            return CombinationOfCones({self:other})
        elif type(other) == int:
            return self * Integer(other)
        else:
            return NotImplemented

    def __add__(self,other):
        if type(other) == SymbolicCone:
            C = CombinationOfCones(self)
            C += other
            return C
        else:
            return NotImplemented

    def __repr__(self):
        return "[ Cone: generators = %s, apex = %s, openness = %s ]" % (self.generators.__repr__(), self.apex.__repr__(), self.openness.__repr__())

    def canonicalize(self):
        """Permutes generators of self so that they are in lexicographic order."""
        # We have to permute the generators and the openness vector with the same permutation.
        # To do this, we zip both lists and then permute according to the first argument.
        pairs = zip(self._generators,self._openness)
        pairs.sort(cmp = lambda x,y: lex_cmp(x[0],y[0]))
        gens, opens = zip(*pairs) # zip also does unzip!
        self._generators = tuple(gens)
        self._openness = tuple(opens)

    def generator_matrix(self):
        """Returns a matrix with generators as columns."""
        return column_matrix(ZZ,self._generators)

    def snf_divisors(self):
        """Computes the list of elements on the diagonal of the Smith Normal Form of the matrix of generators."""
        smith_form = self.generator_matrix().smith_form()[0]
        return [smith_form[i][i] for i in xrange(self.dimension)]

    def index(self):
        """Computes the index (the determinant) of self."""
        return prod(self.snf_divisors())

    def flip(self):
        # we have to use lists instead of generators here, because zip does not work for generators
        J = [ not vector_is_forward(v) for v in self.generators ] # 1 = backward, 0 = forward
        Jpm = [ -_1 if Ji else _1 for Ji in J ] # -1 = backward, 1 = forward
        s = prod(Jpm) # s = (-1)**(number of backward generators)
        _o = tuple(( 1 if oi != Ji else 0 for (oi,Ji) in zip(self.openness,J) ))
        _V = tuple(( msv(Jpmi, v) for (Jpmi,v) in zip(Jpm,self.generators) ))
        return s * SymbolicCone(_V,self.apex,_o)

    def eliminate_last_coordinate(self,equality=False):
        """Eliminates the last coordinate of self by intersection.

        More precisely, eliminate_last_coordinate computes a Brion decomposition of self intersected with the half-space in which the last coodinate is non-negative, and projects the result by forgetting the last coordinate. The elements of the Brion decomposition are all simplicial cones. The formula for decomposition assumes that the projection is one-to-one, i.e., the affine hull of self intersects the last-coordinate-equal-zero-hyperplane in codimension one. (This assumption holds, e.g., then the coordinate generated was introduced by MacMahon-lifting.)

        Also supports intersection with the last-coordinate-equal-zero-hyperplane instead of the half-space where the last coordinate is non-negative. To this end, the argument equality has to be set to True.

        Arguments:
            equality: True if we intersect with hyperplane instead of half-space.
        """
        # abbreviate variables
        V = self.generators
        q = self.apex
        o = self.openness
        k = self.dimension
        lambd = len(V[0]) - 1 # index of last coordinate
        w = lambda j: svv(q, msv(-q[lambd]/(V[j][lambd]), V[j])) # q - q[lambd]/(V[j][lambd]) * V[j]
        sgqn = 1 if q[lambd] >= 0 else -1

        def g(i,j):
            if i == j: # TODO: check if this makes sense for e == 1 as well
                return msv(-1,V[j]) if q[lambd] >= 0 else V[j]
            else:
                return msv(sgqn, svv( msv(V[i][lambd], V[j]), msv(-V[j][lambd], V[i]))) # sg(qn) * ( V[i][lambd] * V[j] - V[j][lambd] * V[i] )
        if not equality:
            G = lambda j: prim(( g(i,j)[:-1] for i in xrange(k) ))
        else:
            G = lambda j: prim(( g(i,j)[:-1] for i in xrange(k) if j != i ))
        def _o(j):
            if equality: # drop j-th element
                return o[:j] + o[j+1:]
            else: # replace j-th element with zero
                return o[:j] + (0,) + o[j+1:]
        def Bplus():
            return CombinationOfCones.sum(( SymbolicCone(G(j), w(j)[:-1], _o(j)) for j in xrange(k) if V[j][lambd] > 0 ))
        def Bminus(): # the Bminus in the paper is this Bminus + Cprime
            return CombinationOfCones.sum(( SymbolicCone(G(j), w(j)[:-1], _o(j)) for j in xrange(k) if V[j][lambd] < 0 ))
        def Cprime():
            return CombinationOfCones(SymbolicCone(prim([v[:-1] for v in V]), q[:-1], o))

        n_up = len([Vj for Vj in V if Vj[lambd] > 0])
        n_down = len([Vj for Vj in V if Vj[lambd] < 0])
        n_flat = k - n_up - n_down
        if not equality:
            if q[lambd] < 0:
                B = Bplus()
            else: # if q[lambd] >= 0
                B = Bminus() + Cprime()
            # TODO: optimization: if we have a choice, take decomposition with fewer cones
        else: # if equality
            if q[lambd] < 0 or (q[lambd] == 0 and exists(range(k),lambda j: V[j][lambd] > 0)[0]):
                B = Bplus()
            elif q[lambd] > 0 or (q[lambd] == 0 and exists(range(k),lambda j: V[j][lambd] < 0)[0]):
                B = Bminus()
            else: # q[lambd] == 0 and forall(range(k), lambda j: V[j][lambd] == 0)[0]
                B = Cprime()

        #print "intermediate result before flip", B
        result = B.map_cones(lambda C: C.flip())
        #print "yields result", result
        return result

    def enumerate_fundamental_parallelepiped(self):
        """Returns a list of all integer points in the fundamental parallelepiped of self.

        The integer points in the list are represented as tuples of Sage integers. Their number is self.index(). Note that this list may be exponentially large; its computation may therefore take a long time and consume a lot of memory.
        """
        k = self.dimension
        d = self.ambient_dimension
        # note: in Sage the Smith form is UVW = S, whereas in the draft we use V=USW.
        V = self.generator_matrix()
        q = vector(self.apex,QQ)
        p = vector([0 for i in range(d)],ZZ)
        if k < d:
            # find an integer point in the affine span
            A,b = self.inequality_description_of_affine_hull()
            p = solve_linear_equation_system_over_integers(A,b)
            if not is_integral(p):
                return []
        S,Uinv,Winv = V.smith_form()
        s = [S[i][i] for i in xrange(k)] # we don't need the last d-k 1s: + [1 for i in range(d-k)]
        #T = diagonal_matrix([1/s[i] for i in xrange(k)]).augment(zero_matrix(QQ,k,d-k))
        qhat = Uinv * (q - p)
        _L = CartesianProduct( *[xrange(s[i]) for i in xrange(k)] )
        #WinvT = Winv * T  
        sk = s[k-1]
        sprime = [Integer(sk / si) for si in s]
        Wprime = tuple([tuple([ (Winv[j,i] * sprime[i]) for i in xrange(k)]) for j in xrange(k)])
        qtrans = [ sum([ - Wprime[j][i] * qhat[i]  for i in xrange(k)]) for j in xrange(k) ]
        qfrac = fract_simple(qtrans)
        qint = [ floor(qi) for qi in qtrans ]
        qsummand = tuple((Integer(qi) for qi in sk * q + V * vector(qfrac) )) # this vector is integral
        o = [ (self.openness[j] if qfrac[j] == 0 else 0) for j in xrange(k) ]
        def _transform_integral(z):
            innerRes = []
            j = 0
            for qj in qint: # qint has k entries
                inner = 0
                i = 0
                for zi in z: # z should have k entries
                    inner += Wprime[j][i] * zi # Wprime[i,j] * vj
                    # the following is optional, depending on performance
                    #inner = inner % sk
                    i += 1
                inner += qj
                inner = inner % sk
                # inner has to be modified according to I
                if inner == 0 and o[j]:
                    inner = sk
                innerRes.append(inner)
                j += 1
            outerRes = []
            for l in xrange(d):
                outer = 0
                j = 0
                for innerResi in innerRes:
                    outer += V[l,j] * innerResi
                    j += 1
                outerRes.append(outer) # outerRes is an integral vector
            return tuple(( (ai + bi).divide_knowing_divisible_by(sk) for (ai,bi) in zip(outerRes,qsummand) ))
        result = [ _transform_integral(v) for v in _L ]
        return result

    def inequality_description_of_affine_hull(self):
        k = self.dimension
        d = self.ambient_dimension
        if d == k: 
            # as we assume that V is full rank, this implies that the cone is full dimensional
            # therefore the system describing the affine hull is empty
            return ([],[])
        V = self.generator_matrix()
        q = vector(self.apex,QQ)
        P, L, U = V.LU()
        A = L**(-1)*P**(-1)
        Ap = A[-(d-k):]
        b = Ap * q
        # this last part makes sure that the system Ax=b is integral
        m = Ap.nrows()
        multipliers = [ lcm([ Ap[i,j].denominator() for j in range(d) ] + [b[i].denominator()]) for i in range(m) ]
        S = diagonal_matrix(multipliers)
        result = (S*Ap, S*b)
        return result

    @staticmethod
    def macmahon_cone(A,b):
        """Takes a linear system and returns the corresponding MacMahon cone.

        Arguments:
            A: An integer matrix represented as a tuple of rows. Entries can be ints or Sage Integers.
            b: An integer vector, represented as a tuple. Entries can be ints or Sage Integers.
            E: A 0-1 vector which encodes with rows are equalities. (1 means equality, 0 means inequality.)
        """
        if type(b) != tuple and type(b) != list:
            raise TypeError("Right-hand side b should be a tuple or list of integers, but is %s" % b.__repr__())
        _A = matrix(ZZ,sageify(A))
        _b = vector(ZZ,sageify(b))
        d = _A.ncols()
        
        generators = tuple(( tuple(( vij for vij in column )) for column in identity_matrix(ZZ,d).stack(_A).columns() ))
        apex = tuple( [0 for i in xrange(d)] + [-1 * bi for bi in _b] )
        openness = tuple((0 for i in xrange(d)))

        return SymbolicCone(generators, apex, openness)

    @staticmethod
    def cone_from_homogeneous_system(A):
        """Takes a matrix describing a homogeneous system of linear inequalities and returns the cone they define.

        Arguments:
            A: A matrix that defines a homogeneous linear inequality system Ax >= 0. A is assumed to be square and of full rank. A has to be a Sage matrix.
        """

        # Implementation: The columns of A^(-1) are the generators of the cone define by Ax >= 0. Of course these columns have to be brought in integer form.

        V = A**(-1)
        d = A.dimensions()[0]
        generators = V.columns()
        integer_generators = tuple((prim_v(intify_v(g)) for g in generators))
        apex = tuple((0 for i in range(d)))
        openness = tuple((0 for i in range(d)))

        return SymbolicCone(integer_generators, apex, openness)

class CombinationOfCones(dict):
    def __init__(self,*args,**kwargs):
        if len(args) == 1 and len(kwargs) == 0 and type(args[0]) == SymbolicCone:
            super(CombinationOfCones,self).__init__({args[0]: _1})
        else:
            super(CombinationOfCones,self).__init__(*args,**kwargs)

    def __add__(self,other):
        """Returns a new CombinationOfCones that is the sum of self and another SymbolicCone or CombinationOfCones."""
        if type(other) != SymbolicCone and type(other) != CombinationOfCones:
            return NotImplemented
        result = CombinationOfCones()
        result += self
        result += other
        return result

    def __mul__(self,other):
        """Returns a new CombinationOfCones that is self multiplied with a scalar factor given by other."""
        if type(other) != Integer and type(other) != int:
            return NotImplemented
        result = CombinationOfCones()
        result += self
        result *= other
        return result        

    def __iadd__(self,other):
        """Adds another SymbolicCone or CombinationOfCones to self. Updates self in-place."""
        if type(other) == SymbolicCone:
            self.add_cone(other)
        elif type(other) == CombinationOfCones:
            for cone, multiplicity in other.iteritems():
                self.add_cone(cone,multiplicity)
        else:
            return NotImplemented
        return self

    def __imul__(self,other):
        """Multiplies self with an integer other. Updates self in-place."""
        if type(other) == Integer:
            if other == _0:
                self.clear()
            else:
                for cone in self.keys():
                    self[cone] = self[cone] * other
            return self
        elif type(other) == int:
            self.__imul__(Integer(other))
            return self
        else:
            return NotImplemented

    def __repr__(self):
        r = []
        for cone, multiplicity in self.iteritems():
            r.append("%d * %s " % (multiplicity,cone.__repr__()) )
        return " + ".join(r)

    def add_cone(self,cone,multiplicity=_1):
        """Adds cone to this linear combination of cones.

        Self is modified in place. Multiplicity of added cone may be given. Returns self.

        Arguments:
            cone: The SymbolicCone to be added to self.
            multiplicity: The multiplicity of cone, given as an Integer.
        """
        if multiplicity != _0:
            old_m = self[cone] if self.has_key(cone) else _0
            new_m = old_m + multiplicity
            if new_m == _0 and self.has_key(cone):
                del self[cone]
            else:
                self[cone] = new_m
        return self

    def map_cones(self,f):
        """Applies function f to every cone in self and returns the linear combination of results.

        f is a function taking a SymbolicCone and returning a CombinationOfCones. If self represents a combination of cones a_1 * c_1 + ... + a_k * c_k, then self.map_cones(f) returns the combination a_1 * f(c_1) + ... + a_k * f(c_k). If f(c_i) is itself a linear combination, then the distributive law is applied to flatten the sum and collect terms. map_cones not modify self, but instead returns a new CombinationOfCones instance.

        Arguments:
            f: A function of a single argument, mapping a SymbolicCone to a CombinationOfCones.
        """
        result = CombinationOfCones()
        for cone, multiplicity in self.iteritems():
            # ensure all updates happen in-place through the use of += and *=
            combination = f(cone)
            combination *= multiplicity 
            result += combination
        return result

    @staticmethod
    def sum(gen,first=None):
        """Takes a list or generator of SymbolicCones and CombinationOfCones and returns their sum.

        Arguments:
            gen: A list or generator of SymbolicCones and CombinationOfCones. Both types can be mixed.
            first: An initial CombinationOfCones to which all elements in gen are added. First is modified in-place. If first is omitted a new empty CombinationOfCones is created.
        """
        if first == None:
            first = CombinationOfCones()
        result = first
        for c in gen:
            result += c
        return result
