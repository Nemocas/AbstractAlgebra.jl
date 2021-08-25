export RationalFunctionField

function RationalFunctionField(k::Field, s::Symbol; cached=true)
   return Generic.RationalFunctionField(k, s; cached=cached)
end

function RationalFunctionField(k::Field, s::AbstractString; cached=true)
   return Generic.RationalFunctionField(k, Symbol(s); cached=cached)
end



