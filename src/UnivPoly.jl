###############################################################################
#
#   UnivPoly.jl: Generic universal polynomial ring (variables can be added)
#
###############################################################################

export UniversalPolynomialRing

###############################################################################
#
#   Iterators
#
###############################################################################

function coefficients(a::AbstractAlgebra.UniversalPolyRingElem)
   return Generic.UnivPolyCoeffs(a)
end

function exponent_vectors(a::AbstractAlgebra.UniversalPolyRingElem)
   return Generic.UnivPolyExponentVectors(a)
end

function monomials(a::AbstractAlgebra.UniversalPolyRingElem)
   return Generic.UnivPolyMonomials(a)
end

function terms(a::AbstractAlgebra.UniversalPolyRingElem)
   return Generic.UnivPolyTerms(a)
end

###############################################################################
#
#   UniversalPolynomialRing constructor
#
###############################################################################

function UniversalPolynomialRing(R::Ring; ordering=:lex, cached=true)
   return Generic.UniversalPolynomialRing(R; ordering=ordering, cached=cached)
end

