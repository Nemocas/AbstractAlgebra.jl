###############################################################################
#
#  Transcendental functions
#
###############################################################################

function Base.log(a::RelSeriesElem{T}) where T <: FieldElement
    @assert valuation(a) == 0 
    return integral(derivative(a)*inv(a))
end

function Base.exp(a::RelSeriesElem{T}) where T <: FieldElement
    if iszero(a)
      b = parent(a)(1)
      set_precision!(b, precision(a))
      return b
    end
    @assert valuation(a) > 0
    R = base_ring(parent(a))
    x = parent(a)([R(1)], 1, min(2, precision(a)), 0)
    p = precision(a)
    la = [p]
    while la[end] > 1
        push!(la, div(la[end] + 1, 2))
    end
  
    one = parent(a)([R(1)], 1, 2, 0)
  
    n = length(la) - 1
    # x -> x*(1 - log(a) + a) is the recursion
    while n > 0
        set_precision!(x, la[n])
        set_precision!(one, la[n])
        x = x*(one - log(x) + a) # TODO: can be optimized...
        n -= 1 
    end
    return x
end

###############################################################################
#
#  Transcendental functions
#
###############################################################################

@doc Markdown.doc"""
    derivative(f::RelSeriesElem{T}) -> RelSeriesElem

Return the derivative of the power series $f$.
"""
function derivative(f::RelSeriesElem{T}) where T <: RingElement
  g = parent(f)()
  set_precision!(g, precision(f) - 1)
  fit!(g, pol_length(f)) # FIXME: this function is not part of the interface
  v = valuation(f)
  set_valuation!(g, 0)
  if v == 0
    for i = 1:pol_length(f)
      setcoeff!(g, i - 1, (i + v)*polcoeff(f, i))
    end
  else
    for i = 0:pol_length(f)
      setcoeff!(g, i, (i + v)*polcoeff(f, i))
    end
    set_valuation!(g, v - 1)
  end
  renormalize!(g)  
  return g
end

@doc Markdown.doc"""
    integral(f::RelSeriesElem{T}) -> RelSeriesElem

Return the integral of the power series $f$.
"""
function integral(f::RelSeriesElem{T}) where T <: RingElement
    g = parent(f)()
    fit!(g, precision(f) + 1)
    set_precision!(g, precision(f) + 1)
    v = valuation(f)
    set_valuation!(g, v + 1)
    for i = 0:pol_length(f)
        c = polcoeff(f, i)
        if !iszero(c)
            setcoeff!(g, i, divexact(c, i + v + 1))
        end
    end
    renormalize!(g) 
    return g
end
