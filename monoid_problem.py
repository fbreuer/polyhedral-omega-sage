from geometry import *
from util import *
from linear_algebra import *

def matrix_of_constraints(N):
    divisors = N.divisors()[1:] # divisors strictly greater than 1
    codivisors = [ N.divide_knowing_divisible_by(di) for di in divisors ]
    n = len(divisors) # number of divisors
    a = lambda i,j: ( divisors[i] * ( gcd(divisors[j],codivisors[i])**2 - divisors[j] )) / (24 * gcd(divisors[i], codivisors[i]) * divisors[j])
    A = matrix(QQ,[ [ a(i,j) for j in range(n) ] for i in range(n)] )
    return A

def modular_constraints(N):
    divisors = N.divisors()[1:] # divisors strictly greater than 1
    prime_divisors = N.prime_divisors()
    def f(p,delta):
        # the maximal j such that p^j | delta
        factors = dict(list(delta.factor()))
        if factors.has_key(p):
            return factors[p]
        else:
            return Integer(0)
    constraints = []
    for p in prime_divisors:
        # one constraint for each prime divisor
        a = [ f(p,delta) for delta in divisors ]
        m = 2
        constraints.append((a,m))
    constraints.append( ([delta - 1 for delta in divisors], 24) )
    constraints.append( ([(N.divide_knowing_divisible_by(delta)) - N for delta in divisors], 24) )
    return constraints

def construction(N,omega=False,compute=False):
    Cs, fp, H, GH = None, None, None, None
    print "N = ", N
    print "divisors = ", N.divisors()
    n = len(N.divisors()[1:])
    print "number of divisors great than 1 = ", n
    A = matrix_of_constraints(N)
    print "the constraint matrix A:"
    print A
    congruences = modular_constraints(N)
    print "congruences defining the sublattice:"
    print congruences
    G = matrix(ZZ,congruence_sublattice(congruences)).transpose()
    print "the sublattice transformation:"
    print G
    print "the transformed constraint matrix:"
    print A*G
    if omega:
        print
        print "Omega Construction"
        print "input matrix:"
        Aprime = tupleize(intify_m(A*G))
        print Aprime
        C = SymbolicCone.macmahon_cone(Aprime,tuple((0 for i in range(n))))
        print "macmahon cone:"
        print C
        Cs = CombinationOfCones({C: Integer(1)})
        for i in range(n):
            Cs = Cs.map_cones(lambda c: c.eliminate_last_coordinate())
        print "combination of cones after elimination:"
        print Cs
        print "indices:"
        indices = [ C.index() for C in Cs.keys() ]
        print indices
        print "sum of indices:"
        print sum(indices)
    print
    print "Direct Construction"
    C = SymbolicCone.cone_from_homogeneous_system(A*G)
    print "cone:"
    print C
    print "SNF elementary divisors:"
    print C.snf_divisors()
    print "index:"
    print C.index()
    if compute:
        print
        print "Computation..."
        print "Enumerating fundamental parallelepiped..."
        fp = C.enumerate_fundamental_parallelepiped()
        print "Filtering Hilbert basis..."
        H = hilbert_basis(C,fp)
        print "number of elements of Hilbert basis = ", len(H)
        GH = [ G * vector(ZZ,h) for h in H ]
    print
    return (A,A*G,G,n,Cs,C,fp,H,GH)

def is_zero(v):
    return forall(v,lambda vi: vi == 0)[0]

def is_composite(ws,w):
    for u in ws:
        if not is_zero(u):
            mu = msv(-1,u)
            wmu = svv(mu,w)
            if not is_zero(wmu) and wmu in ws:
                return True
    return False

def fp_hilbert_basis(fp):
    basis = []
    fp = set(fp)
    for v in fp:
        if not is_composite(fp,v):
            if not is_zero(v):
                basis.append(v)
    return basis

def hilbert_basis(cone,fp):
    basis = fp_hilbert_basis(set(fp))
    basis = basis + list(cone.generators)
    return basis

def _construction_barvinok(N):
    r = construction(N)
    A = r[0]
    AG = r[1]
    G = r[2]
    # want x : Ax >= 0
    d = AG.dimensions()[0]
    C = r[5]
    generators = C.generators
    fund_par_offsets = []
    for i in range(d):
        g = generators[i]
        offset = (AG * vector(ZZ,g))[i]
        fund_par_offsets.append(offset)
    sum_of_generators = sum([vector(ZZ,g) for g in generators])

    vec2str = lambda vec: "[%s]" % (",".join([str(vi) for vi in vec]))
    vec2str_ = lambda vec: "{ %s }" % vec2str(vec)
    row2str = lambda vs, row: " + ".join([ "(%d*%s)" % (row[i],vs[i]) for i in range(d) ])

    xvars = [ "x%d" % i for i in range(d) ]
    yvars = [ "y%d" % i for i in range(d) ]
    zvars = [ "z%d" % i for i in range(d) ]
    xvar_str = vec2str(xvars)
    yvar_str = vec2str(yvars)
    zvar_str = vec2str(zvars)
    zero_str = vec2str([0 for i in range(d)])
    generator_set_strs = [ vec2str_(gen) for gen in generators ]

    rowconstraint = lambda i, vs: "0 <= %s < %d" % (row2str(vs,AG.rows()[i]),fund_par_offsets[i])
    xconstraints = " and ".join([ rowconstraint(i,xvars) for i in range(d)])

    Qstr = "Q := { " + xvar_str + ": " + xconstraints + " } - { " + zero_str + " } + " + "+".join(generator_set_strs) + ";\n"

    exists_str = "exists " + ",".join(yvars + zvars) + ":"
    yconstraints = " and ".join([ rowconstraint(i,yvars) for i in range(d)]) + " and " + row2str(yvars,sum_of_generators) + " > 0 "
    zconstraints = " and ".join([ rowconstraint(i,zvars) for i in range(d)]) + " and " + row2str(zvars,sum_of_generators) + " > 0 "
    linkconstraints = " and ".join(["%s = %s + %s" % (xvars[i],yvars[i],zvars[i]) for i in range(d)])

    Sstr = "S := { " + xvar_str + ": " + exists_str + " " + yconstraints + "\n and " + zconstraints + "\n and " + linkconstraints + " };\n"

    Hstr = "H := Q - S;\n"

    Commandstr = "card H;\n"

    return Qstr + Sstr + Hstr + Commandstr

def _compute_barvinok(N,filename):
    _barvinok_str = _construction_barvinok(N)
    print "barvinok problem, isl/iscc format:"
    print _barvinok_str
    print
    print "writing to file ", filename
    f = open(filename,"w")
    f.write(_barvinok_str)
    f.close()


def _construction_normaliz(N):
    r = construction(N)
    A = r[0]
    AG = r[1]
    G = r[2]
    M = AG
    dim_str = "%d\n%d\n" % M.dimensions()
    mat_str = "\n".join([ " ".join([ str(mij) for mij in row ]) for row in M.rows() ])
    return dim_str + mat_str + "\nhyperplanes"   

def _compute_normaliz(N,filename):
    _normaliz_str = _construction_normaliz(N)
    print "normaliz hilbert problem, lattice format:"
    print _normaliz_str
    print
    print "writing to file ", filename
    f = open(filename,"w")
    f.write(_normaliz_str)
    f.close()

def _construction_4ti2(N):
    r = construction(N)
    A = r[0]
    AG = r[1]
    G = r[2]
    M = AG.transpose()
    dim_str = "%d %d\n" % M.dimensions()
    mat_str = "\n".join([ " ".join([ str(mij) for mij in row ]) for row in M.rows() ])
    return dim_str + mat_str

def _compute_4ti2(N,filename):
    _4ti2_str = _construction_4ti2(N)
    print "4ti2 hilbert problem, lattice format:"
    print _4ti2_str
    print
    print "writing to file ", filename
    f = open(filename,"w")
    f.write(_4ti2_str)
    f.close()

def _read_result_4ti2(filename,A,G=None,C=None):
    f = open(filename,"r")
    result = f.read()
    f.close()
    lines = result.splitlines()
    dimensions_line = lines[0].split(" ")
    n_rows = int(dimensions_line[0])
    n_cols = int(dimensions_line[1])
    M = [ [ Integer(e_str) for e_str in row_str.split(" ") ] for row_str in lines[1:] ]
    # the rows of this matrix are of the form zL for L = (A*G).transpose()
    # to get at the elements z, we need to compute (A*G)**(-1) * (zL).T = L.T**(-1) * (zL).T = L.T**(-1) * L.T * z = z
    # In the next step we then have to apply G to z, to map the integer lattice to the congruence sublattice
    # G * (A*G)**(-1) * zL = G * G**(-1) * A**(-1) * zL = A**(-1) * zL
    Ainv = A**(-1)
    result = [ Ainv * vector(ZZ,row) for row in M ]
    print "size of Hilbert basis: %d" % len(result)
    if G != None and C != None:
        V = matrix(ZZ,C.generators)
        d = len(C.generators[0])
        GV = (G*V.T).T
        h = GV**(-1) * vector(ZZ,[Integer(1) for i in range(d)])
        heights = [ float(h * Hi) for Hi in result ]
        max_height = max(heights)
        print "maximal height of a Hilbert basis in the fundamental parallelepiped is %f" % max_height
        print "this is %2.1f%% of the theoretical maximum " % (100.0 * max_height / float(d))
    return result

def compute_order(N,r,debug=False):
    all_divs = N.divisors()
    non_triv_divs = all_divs[1:]
    r0 = - sum(r)
    full_r = [r0] + list(r)
    s = sum([ full_r[i] * all_divs[i] for i in range(len(full_r))])
    if debug:
        if s % 24 != 0:
            print "ERROR: ", s, r
    return s / 24

def report_on_orders(N):
    data = construction(N)
    A = data[0]
    G = data[2]
    C = data[5]
    result = _read_result_4ti2("hilbert_%d.hil" % N,A,G=G,C=C)
    orders = [ compute_order(N,r) for r in result ]
    orders.sort()
    print "max order: ", orders[0]
    print "first ten orders: ", orders[:10]
    print "last ten orders: ", orders[-10:]
    unique_orders = list(set(orders))
    print "number of distinct orders: ", len(unique_orders)
    return (orders,unique_orders)
