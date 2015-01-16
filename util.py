from sage.all_cmdline import *

def sageify(arg):
    """Takes a nested tuple and converts all ints therein into Sage Integers."""
#    print "."
#    print type(arg)
    if type(arg) == tuple:
        return tuple(( sageify(e) for e in arg ))
    elif type(arg) == int:
        return Integer(arg)
    elif type(arg) == list:
        return [ sageify(e) for e in arg ]
    elif type(arg) == Integer:
        return arg
    elif type(arg) == Rational:
        return arg
    #elif type(arg) == Matrix_rational_dense:
    #    return tupleize(arg)
    #elif type(arg) == Matrix_integer_dense:
    #    return tupleize(arg)
    else:
        raise TypeError("Type %s is not supported by sageify." % str(type(arg)))

def tupleize(arg):
    if type(arg) == sage.matrix.matrix_rational_dense.Matrix_rational_dense or type(arg) == sage.matrix.matrix_integer_dense.Matrix_integer_dense:
        return tuple((tupleize(row) for row in arg.rows()))
    elif type(arg) == sage.modules.vector_rational_dense.Vector_rational_dense or tupe(arg) == sage.modules.vector_integer_dense.Vector_integer_dense:
        return tuple((vi for vi in arg))
    else:
        raise TypeError("Type %s is not supported by sageify." % str(type(arg)))