###############################################################################
#
#  Power testing and root
#
###############################################################################

function is_power(a::RingElem, n::Int)
    if isone(a) || iszero(a)
        return true, a
    end
    if isone(-a) && isodd(n)
        return true, a
    end
    R = parent(a)
    Rt = AbstractAlgebra.poly_ring(R)
    x = gen(Rt)
    r = roots(x^n - a)
    if length(r) == 0
        return false, a
    else
        return true, r[1]
    end
end

function root(a::RingElem, n::Int)
    fl, b = is_power(a, n)
    fl || error("element does not have a $n-th root")
    return b
end

@doc raw"""
    falling_factorial(x::RingElement, n::Integer)

Return the falling factorial of $x$, i.e. $x(x - 1)(x - 2)\cdots (x - n + 1)$.
If $n < 0$ we throw a `DomainError()`.

# Examples

```jldoctest
julia> R, x = ZZ[:x];

julia> falling_factorial(x, 1)
x

julia> falling_factorial(x, 2)
x^2 - x

julia> falling_factorial(4, 2)
12
```
"""
falling_factorial(x::RingElement, n::Integer) = (-1)^n*rising_factorial(-x, n)

@doc raw"""
    rising_factorial(x::RingElement, n::Integer)

Return the rising factorial of $x$, i.e. $x(x + 1)(x + 2)\cdots (x + n - 1)$.
If $n < 0$ we throw a `DomainError()`.

# Examples

```jldoctest
julia> R, x = ZZ[:x];

julia> rising_factorial(x, 1)
x

julia> rising_factorial(x, 2)
x^2 + x

julia> rising_factorial(4, 2)
20
```
"""
function rising_factorial(x::RingElement, n::Integer)
  n < 0 && throw(DomainError(n, "Argument must be non-negative"))
  n == 0 && return one(x)
  n == 1 && return x
  if x isa Integer
      return reduce(Base.checked_mul, x+i-1 for i in 1:Int(n))
  end
  return prod(x+i-1 for i in 1:Int(n))
end

@doc raw"""
    rising_factorial2(x::RingElement, n::Integer)

Return a tuple containing the rising factorial $x(x + 1)\cdots (x + n - 1)$
and its derivative.
If $n < 0$ we throw a `DomainError()`.

# Examples

```jldoctest
julia> R, x = ZZ[:x];

julia> rising_factorial2(x, 1)
(x, 1)

julia> rising_factorial2(x, 2)
(x^2 + x, 2*x + 1)

julia> rising_factorial2(4, 2)
(20, 9)
```
"""
function rising_factorial2(x::RingElement, n::Integer)
  n < 0 && throw(DomainError(n, "Argument must be non-negative"))
  n == 0 && return (one(x), zero(x))
  n == 1 && return (x, one(x))
  f, F = rising_factorial2(x, Int(n)-1)
  # use product rule:  [(x+n-1)*f(x)]' = (x+n-1)'*f(x) + (x+n-1)*f'(x)
  if x isa Integer
      return Base.checked_mul(x+n-1, f),  f + Base.checked_mul(x+n-1, F)
  end
  return (x+n-1)*f,  f + (x+n-1)*F
end
