##############################################################################
#
#  Evaluation
#
##############################################################################

function evaluate(f::Generic.Frac{T}, V::Vector{U}) where {T <: RingElement, U <: RingElem}
    return evaluate(numerator(f), V)//evaluate(denominator(f), V)
end
  
function evaluate(f::Generic.Frac{T}, v::U) where {T <: RingElement, U <: RingElement}
    return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

function evaluate(f::Generic.Frac{T}, v::U) where {T <: PolyElem, U <: Integer}
    return evaluate(numerator(f), fmpz(v))//evaluate(denominator(f), fmpz(v))
end
  
