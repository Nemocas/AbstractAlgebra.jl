###############################################################################
#
#   UnivPoly.jl: Generic universal polynomial ring (variables can be added)
#
###############################################################################

export UniversalPolynomialRing

function show(io::IO, p::UnivPolyRing)
   print(io, "Universal Polynomial Ring")
end

###############################################################################
#
#   UniversalPolynomialRing constructor
#
###############################################################################

function UniversalPolynomialRing(R::Ring; ordering=:lex, cached=true)
   return Generic.UniversalPolynomialRing(R; ordering=ordering, cached=cached)
end

