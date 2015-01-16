import lds
import geometry

def remove(V,v):
    return [ vi for vi in V if vi != v ]

def rem(L,e):
    """Returns a copy of L that has the first element equal to e removed."""
    result = list(L)
    result.remove(e)
    return result

def lists_of_length(length,V):
    """All lists of length distinct elements from V"""
    if length == 0:
        return [ [] ]
    else:
        return [ [v]+l for v in V for l in lists_of_length(length-1,rem(V,v)) ]

def all_lists(min_length,max_length,V):
    """All lists of length between min_length and max_length (inclusive) from V."""
    return [ l for i in range(min_length,max_length+1) for l in lists_of_length(i,V) ] 

def remove_duplicates(rel,L):
    """Remove all duplicates from list L where duplicates are defined by the equivalence relation rel. Does not modify the list L but returns a copy."""
    result = []
    for x in L:
        if not any(rel(x,y) for y in result): # if there does not exist an element y in result that is in relation to x
            result.append(x)
    return result

def is_cyclic_permutation(a,b):
    """Checks if list a is a cyclic permutation of b."""
    return a in cyclic_permutations(b)

def cyclic_permutations(a):
    """Returns all cyclic permutations of a."""
    return [ a[i:] + a[:i] for i in range(len(a)) ]

def dag_lds(n):
    """Creates the polytope of all DAGs on n vertices."""

    vertices = range(n)
    edges = [ (v,w) for v in vertices for w in vertices if v != w ] # exclude loops

    def weighted_edge_set_to_vector(es):
        """es is a dictionary mapping edges to multiplicities."""
        result = []
        for e in edges:
            if e in es:
                result.append(es[e])
            else:
                result.append(0)
        return result

    A = []
    E = []
    b = []

    for e in edges:
        # the inequality x_e >= 0 is implicit
        # the inequality x_e <= 1 we have to define
        # x_e <= 1 is equivalent to -x_e >= -1
        A.append(weighted_edge_set_to_vector({ e: -1 }))
        b.append(-1)
        E.append(0)

    cycles = all_lists(2,n,vertices)
    # unfortunately, the above is highly redundant, as all cyclic permutations of a cycle yield the same thing!
    cycles = remove_duplicates(is_cyclic_permutation,cycles)
    # print cycles

    for c in cycles:
        # we want sum_{e in C} x_e <= |C| - 1
        # equivalent: sum_{e in C} - x_e >= 1 - |C|
        es = {}
        last_v = c[-1]
        for v in c:
            es[(last_v,v)] = -1
            last_v = v
        A.append(weighted_edge_set_to_vector(es))
        b.append(1 - len(c))
        E.append(0)

    return lds.LinearDiophantineSystem(A,b,E)

def sorted_gen_matrix_from_cone(c):
    gens = c._generators
    lgens = list(gens)
    lgens.sort(cmp=lex_cmp)
    return tuple(lgens)


