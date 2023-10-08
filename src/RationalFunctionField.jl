###############################################################################
#
#   rational_function_field.jl : Rational function fields
#
###############################################################################

###############################################################################
#
#   rational_function_field constructor
#
###############################################################################

function rational_function_field(k::Field, s::VarName; cached::Bool=true)
   return Generic.RationalFunctionField(k, Symbol(s); cached=cached)
end

function rational_function_field(k::Field, s::AbstractVector{<:VarName}; cached::Bool=true)
   return Generic.RationalFunctionField(k, [Symbol(si) for si in s]; cached=cached)
end
