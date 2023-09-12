###############################################################################
#
#   UnivPoly.jl: Generic universal polynomial ring (variables can be added)
#
###############################################################################

###############################################################################
#
#   Iterators
#
###############################################################################

function coefficients(a::UniversalPolyRingElem)
   return Generic.UnivPolyCoeffs(a)
end

function exponent_vectors(a::UniversalPolyRingElem)
   return Generic.UnivPolyExponentVectors(a)
end

function monomials(a::UniversalPolyRingElem)
   return Generic.UnivPolyMonomials(a)
end

function terms(a::UniversalPolyRingElem)
   return Generic.UnivPolyTerms(a)
end

###############################################################################
#
#   UniversalPolynomialRing constructor
#
###############################################################################

function UniversalPolynomialRing(R::Ring; ordering=:lex, cached::Bool=true)
   return Generic.UniversalPolynomialRing(R; ordering=ordering, cached=cached)
end
