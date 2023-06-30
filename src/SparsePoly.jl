###############################################################################
#
#   SparsePoly.jl : Generic sparse univariate polynomials over rings
#
###############################################################################

###############################################################################
#
#   SparsePolynomialRing constructor
#
###############################################################################

@doc raw"""
    SparsePolynomialRing(R::Ring, s::VarName cached::Bool = true)

Given a base ring `R` and a string `s` specifying how the generator
(variable) should be printed, return a tuple `S, x` representing the new
polynomial ring $T = R[x1, x2, ...]$ and the generator $x$ of the polynomial
ring. By default the parent object `T` will depend only on `R` and `x` and
will be cached. Setting the optional argument `cached` to `false` will
prevent the parent object `T` from being cached.
"""
SparsePolynomialRing(R::Ring, s::VarName)
