###############################################################################
#
#   TotalFraction.jl : Total ring of fractions
#
###############################################################################

export TotalRingOfFractions

###############################################################################
#
#   TotalRingOfFractions constructor
#
###############################################################################

@doc Markdown.doc"""
    TotalRingOfFractions(R::Ring; cached=true)

Return the parent object of the total ring of fractions over the given base
ring $R$, i.e. the localisation of `R` at the complement of the set of zero
divisors.

If `cached == true` (the default), the returned parent object is cached so
that it will always be returned by a call to the constructor when the same
base ring $R$ is supplied.
"""
function TotalRingOfFractions(R::Ring; cached=true)
   return Generic.TotalRingOfFractions(R; cached=cached)
end

