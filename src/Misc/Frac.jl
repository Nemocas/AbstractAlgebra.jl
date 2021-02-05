##############################################################################
#
#  Derivative
#
##############################################################################

function derivative(f::Generic.Frac{T}, x::T) where {T <: MPolyElem}
    return derivative(f, var_index(x))
end
  
function derivative(f::Generic.Frac{T}, i::Int) where {T <: MPolyElem}
    n = numerator(f)
    d = denominator(f)
    return (derivative(n, i)*d - n*derivative(d, i))//d^2
end

function derivative(f::Generic.Frac{T}) where {T <: PolyElem}
    n = numerator(f)
    d = denominator(f)
    return (derivative(n)*d - n*derivative(d))//d^2
end

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
  