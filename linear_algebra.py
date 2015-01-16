"""A naive linear algebra library implementing vector and matrix operation on tuples.

Curiously, this is faster than Sage's builting matrix and vector classes in some instances, especially, when dealing with lots of small matrices instead of few big matrices. Also, implements basic operations on lattice vectors such as converting lattice vectors into primitive vectors and checking if vectors are forward wrt. the lexicographic orientation.
"""

from sage.all_cmdline import *  

def mvv(v,w):
    """Compute the scalar product of two vectors v and w."""
    return sum(( vi * wi for (vi,wi) in zip(v,w) ))

def msv(s,v):
    """Multiply all entries of the vector v with the scalar s."""
    return tuple(( s * vi for vi in v ))

def svv(v,w):
    """Sum two vectors v and w."""
    return tuple(( vi + wi for (vi,wi) in zip(v,w) ))

def mmv(m,v):
    """Multiply a matrix m with a vector v.

    Matrices are represented as a tuple of rows."""
    return tuple(( svv(mi,v) for mi in m ))

def vector_is_forward(v):
    """Checks if a vector is forward, i.e., first non-zero coordinate is positive."""
    for vi in v:
        if vi != 0:
            return vi > 0
    return True

def intify_v(v):
    """Takes a rational vector v and returns a positive multiple of v that is integer.

    If the smallest possible multiple of v that is integer is required, apply prim_v afterwards. Return value is a Sage tuple of Integers."""
    denominators = [ vi.denominator() for vi in v ]
    l = lcm(denominators)
    return tuple((Integer(vi * l) for vi in v))

def intify_m(A):
    """Takes a Sage matrix and returns a multiple that is integer."""
    d = A.dimensions()
    denominators = [ A[i,j].denominator() for i in range(d[0]) for j in range(d[1])]
    l = lcm(denominators)
    return l*A

def lex_cmp(v,w):
    """Compare two vectors v, w lexicographically. v, w must be of same length."""
    for i in range(len(v)):
        c = cmp(v[i],w[i])
        if c != 0:
            return c
    return 0

def prim_v(v):
    """Takes an integer vector v and returns the corresponding primitive integer vector.

    The primitive integer vector corresponding to v is the least positive multiple of v that is still an integer vector.

    Arguments:
        v: An integer vector, represented as a tuple of Sage Integers.
    """
    d = abs(gcd(v))
    if d == 1:
        return v
    else:
        return tuple((vi.divide_knowing_divisible_by(d) for vi in v))

def prim(V):
    """Takes a tuple of integer vectors and applies prim_v to each of them."""
    return tuple(( prim_v(v) for v in V ))

def is_integral(v):
    return reduce(lambda a,b: a and b, [vi.is_integral() for vi in v])

def fract_simple(v):
    #res = []
    #for vi in v:
    #    res.append(vi - floor(vi))
    #return res
    return tuple(( vi - floor(vi) for vi in v ))

def solve_linear_equation_system_over_integers(A,b):
    """Takes a Sage matrix A and a Sage vector b and solves the system Ax = b over the integers.

    Returns a solution as a Sage integer vector.

    Arguments:
        A: An integer matrix, given as a Sage object. A should be a matrix over ZZ! If A is a matrix over QQ, but with integral entries, then the Smith form will give unexpeceted results! For convenience, this function coerces matrices over QQ into matrices over ZZ, however.
        b: An integer vector, given as a Sage object.
    """
    A = column_matrix(ZZ,A.columns())
    S, Uinv, Vinv = A.smith_form()
    n = A.ncols()
    m = A.nrows()
    p = Uinv * b
    q = vector([ p[i]/S[i][i] for i in range(m) ] + [ 0 for i in range(n-m)],ZZ)
    r = Vinv * q
    return r

def basis_of_integer_solutions_of_linear_equation_system(A):
    """Takes a Sage matrix A and computes a basis of the set of all integer solutions of Ax = 0, i.e., a basis of the kernel of A, over the integers.

    Assumes that A has more columns than rows. A is allowed to have less than full rank. Returns the solution as a list of Sage integer vectors.

    Arguments:
        A: An integer matrix, given as a Sage object. A should be a matrix over ZZ! If A is a matrix over QQ, but with integral entries, then the Smith form will give unexpeceted results! For convenience, this function coerces matrices over QQ into matrices over ZZ, however.
    """
    A = column_matrix(ZZ,A.columns())
    # A is an m by n matrix
    n = A.ncols()
    m = A.nrows()
    k = A.rank()
    S, Uinv, Vinv = A.smith_form()
    # S is m by n
    # Uinv is m by m
    # Vinv is n by n
    # Univ * A * Vinv == S
    # 
    # The kernel of A is the preimage of 0 under A, denoted preim(A,0).
    # The kernel of a product A * B is ker(A*B) = preim(B, ker(A)).
    # The kernel of U is trivial.
    # The kernel of S is the span of the last n - k unit vectors 0^k + ZZ^(n-k).
    # The preimage of the integer span of a set B of vectors is the integer span of the preimages of the vectors in B.
    # As V is unimodular, the preimage of a vector v is exactly V^(-1) v.
    # Thus ker(A) = ker(U * S * V) =  preim(V,ker(U * S)) = preim(V, preim(S, ker(U)))
    #             = preim(V, preim(S, 0)) = preim(V, ker(S)) = preim(V, 0^k + ZZ^(n-k))
    #             = span of the last n-k columns of V^(-1)
    # The desired basis is given by the last n-k columns of Vinv.
    return Vinv.columns()[k:]

def congruence_sublattice(congruences):
    """Takes a list of linear congruences and returns a basis of the sublattice of all integer points that satisfy all of the congruences.

    The basis is returned as a list of Sage integer vectors.

    Arguments:
        congruences: A list of pairs (a,m) where a is a list of integers a1,...,ak and m is an integer defining the congruence a1 * x1 + ... +  ak * xk = 0 mod m. For all tuples in the list the integer k has to be identical.
    """
    n = len(congruences)
    k = len(congruences[0][0])
    rows = []
    i = 0
    for (a,m) in congruences:
        rows.append( a + [ (m if i == j else 0) for j in range(n)] )
        i = i + 1
    A = matrix(ZZ,rows)
    B = basis_of_integer_solutions_of_linear_equation_system(A)
    return [ b[:k] for b in B ]