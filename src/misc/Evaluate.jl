@doc raw"""
    evaluate(E::Expr, vars::Dict)
    evaluate(E::Expr, R::MPolyRing)
    evaluate(E::Expr, R::PolyRing)

Given an rational multivariate expression, evaluate the expression
by replacing the occurring symbols using the dictionary `vars` or
the variables from the polynomial ring `R`.

# Examples

```jldoctest
julia> evaluate(:((x^2 + 1)//(y + 1//2)), Dict(:x => 1, :y => QQ(2)))
4//5

julia> Qx, (x, y) = QQ[:x, :y];

julia> evaluate(:((x^2 + 1)//(y + 1//2)), Qx)
(x^2 + 1)//(y + 1//2)
```
"""
evaluate(::Expr, ::Union{Dict, MPolyRing, PolyRing})

function evaluate(E::Expr, vars::Dict)
  return _eval(E, vars)
end

function evaluate(E::Expr, R::MPolyRing)
  symR = symbols(R) # Vector{Symbol}
  genR = gens(R)
  return evaluate(E, Dict(symR[i] => genR[i] for i in 1:length(symR)))
end

function evaluate(E::Expr, R::PolyRing)
  symR = var(R) # Symbol
  genR = gen(R)
  return evaluate(E, Dict(symR => genR))
end

function _eval(E::Expr, vars)
  E.head != :call && error("Not a rational expression")
  if E.args[1] == :+
    return reduce(+, (_eval(E.args[i], vars) for i in 2:length(E.args)))
  elseif E.args[1] == :*
    return reduce(*, (_eval(E.args[i], vars) for i in 2:length(E.args)))
  elseif E.args[1] == :-
    if length(E.args) == 2
      return -_eval(E.args[2], vars)
    else
      @assert length(E.args) == 3
      return _eval(E.args[2], vars) - _eval(E.args[3], vars)
    end
  elseif E.args[1] == :^
    return _eval(E.args[2], vars)^_eval(E.args[3], vars)
  elseif E.args[1] == ://
    return _eval(E.args[2], vars)//_eval(E.args[3], vars)
  else
    error("Not a rational expression")
  end
end

function _eval(E::Symbol, vars)
  if !haskey(vars, E)
    error("No mapping for symbol :$E")
  end
  return vars[E]
end

function _eval(E::Number, vars)
  return E
end
