###############################################################################
#
#   RationalFunctionField.jl : Rational function fields
#
###############################################################################

###############################################################################
#
#   RationalFunctionField constructor
#
###############################################################################

function RationalFunctionField(k::Field, s::VarName; cached=true)
   return Generic.RationalFunctionField(k, Symbol(s); cached=cached)
end

function RationalFunctionField(k::Field, s::AbstractVector{<:VarName}; cached=true)
   return Generic.RationalFunctionField(k, [Symbol(si) for si in s]; cached=cached)
end
