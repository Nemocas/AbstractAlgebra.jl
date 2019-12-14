###############################################################################
#
#   LaurentPoly.jl : Generic algorithms for abstract Laurent polynomials
#
###############################################################################

base_ring(p::AbstractAlgebra.LaurentPolyElem) = base_ring(parent(p))
