###############################################################################
#
#   LaurentPoly.jl : Generic algorithms for abstract Laurent polynomials
#
###############################################################################

base_ring(p::LaurentPolyElem) = base_ring(parent(p))


###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, p::LaurentPolynomialRing)
   print(io, "Univariate Laurent Polynomial Ring in ")
   print(io, string(var(p)))
   print(io, " over ")
   print(IOContext(io, :compact => true), base_ring(p))
end
