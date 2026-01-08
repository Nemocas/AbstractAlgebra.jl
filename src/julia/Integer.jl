###############################################################################
#
#   Integer.jl : Additional AbstractAlgebra functionality for Julia Integer
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

const JuliaZZ = Integers{BigInt}()

const zz = Integers{Int}()

parent(a::T) where T <: Integer = Integers{T}()

elem_type(::Type{Integers{T}}) where T <: Integer = T

parent_type(::Type{T}) where T <: Integer = Integers{T}

base_ring_type(::Type{<:Integers}) = Union{}   # no base ring

is_exact_type(::Type{T}) where T <: Integer = true

is_domain_type(::Type{T}) where T <: Integer = true

base_ring(::Vector{T}) where T <: Integer = T

# ideal interface
ideal_type(::Type{Integers{T}}) where T <: Integer = Generic.Ideal{T}
ideal(R::Integers, V::Vector) = Generic.Ideal(R, V)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(::Integers{T}) where T <: Integer = T(0)

one(::Integers{T}) where T <: Integer = T(1)

is_unit(a::Integer) = a == 1 || a == -1

is_zero_divisor(a::Integer) = is_zero(a)

canonical_unit(a::T) where T <: Integer = a < 0 ? T(-1) : T(1)

characteristic(::Integers) = 0
is_known(::typeof(characteristic), ::Integers) = true

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::Integer; context = nothing)
    return a
end

function show(io::IO, R::Integers{BigInt})
   print(io, "Integers")
end

function show(io::IO, R::Integers{T}) where T
   print(io, "Integers{$T}()")
end

###############################################################################
#
#   Modular arithmetic
#
###############################################################################

function divrem(a::BigInt, b::BigInt)
   r = mod(a, b)
   q = Base.div(a - r, b)
   return q, r
end

function divrem(a::Int, b::Int)
   r = mod(a, b)
   q = Base.div(a - r, b)
   return q, r
end

function divrem(a::BigInt, b::Int)
   r = mod(a, b)
   q = Base.div(a - r, b)
   return q, r
end

function divrem(a::S, b::T) where {S <: Integer, T <: Integer}
   r = mod(a, b)
   q = Base.div(a - r, b)
   return q, r
end

function div(a::S, b::T) where {S <: Integer, T <: Integer}
   r = mod(a, b)
   q = Base.div(a - r, b)
   return q
end

###############################################################################
#
#   Divides
#
###############################################################################

function divides(a::Integer, b::Integer)
   if iszero(b)
      return iszero(a), b
   end
   q, r = divrem(a, b)
   return iszero(r), q
end

@doc raw"""
    is_divisible_by(a::Integer, b::Integer)

Return `true` if $a$ is divisible by $b$, i.e. if there exists $c$ such that
$a = bc$.
"""
function is_divisible_by(a::Integer, b::Integer)
   if iszero(b)
      return iszero(a)
   end
   r = rem(a, b)
   return iszero(r)
end

function is_divisible_by(a::BigInt, b::BigInt)
   if iszero(b)
      return iszero(a)
   end
   return Bool(ccall((:__gmpz_divisible_p, :libgmp), Cint,
                                             (Ref{BigInt}, Ref{BigInt}), a, b))
end

function is_divisible_by(a::BigInt, b::Int)
   if iszero(b)
      return iszero(a)
   end
   return Bool(ccall((:__gmpz_divisible_ui_p, :libgmp), Cint,
                                        (Ref{BigInt}, Int), a, b < 0 ? -b : b))
end

function is_divisible_by(a::BigInt, b::UInt)
   if iszero(b)
      return iszero(a)
   end
   return Bool(ccall((:__gmpz_divisible_ui_p, :libgmp), Cint,
                                                    (Ref{BigInt}, UInt), a, b))
end

@doc raw"""
    is_associated(a::Integer, b::Integer)

Return `true` if $a$ and $b$ are associated, i.e. if there exists a unit $c$ such that
$a = bc$. For integers, this reduces to checking if $a$ and $b$ differ by a factor of $1$ or $-1$.
"""
function is_associated(a::Integer, b::Integer)
   return a == b || a == -b
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::Integer, b::Integer; check::Bool=true)
   if check
      q, r = divrem(a, b)
      @req is_zero(r) "Not an exact division"
   else
      q = div(a, b)
   end
   return q
end

function divexact(a::BigInt, b::BigInt; check::Bool=true)
   iszero(b) && throw(DivideError())
   q = BigInt()
   if check
      r = BigInt()
      ccall((:__gmpz_tdiv_qr, :libgmp), Nothing,
              (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), q, r, a, b)
      @req is_zero(r) "Not an exact division"
   else
      ccall((:__gmpz_divexact, :libgmp), Nothing,
                              (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), q, a, b)
   end
   return q
end

function divexact(a::BigInt, b::Int; check::Bool=true)
   iszero(b) && throw(DivideError())
   q = BigInt()
   sgn = b < 0
   if check
      r = BigInt()
      ccall((:__gmpz_tdiv_qr_ui, :libgmp), Nothing,
           (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, Int), q, r, a, sgn ? -b : b)
      @req is_zero(r) "Not an exact division"
   else
      ccall((:__gmpz_divexact_ui, :libgmp), Nothing,
                           (Ref{BigInt}, Ref{BigInt}, Int), q, a, sgn ? -b : b)
   end
   return sgn ? -q : q
end

function divexact(a::BigInt, b::UInt; check::Bool=true)
   iszero(b) && throw(DivideError())
   q = BigInt()
   if check
      r = BigInt()
      ccall((:__gmpz_tdiv_qr_ui, :libgmp), Nothing,
                     (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, UInt), q, r, a, b)
      @req is_zero(r) "Not an exact division"
   else
      ccall((:__gmpz_divexact_ui, :libgmp), Nothing,
                                     (Ref{BigInt}, Ref{BigInt}, UInt), q, a, b)
   end
   return q
end

###############################################################################
#
#   Inverse
#
###############################################################################

function inv(a::T) where T <: Integer
   if a == 1
      return one(T)
   elseif a == -1
      return -one(T)
   end
   iszero(a) && throw(DivideError())
   throw(ArgumentError("not a unit"))
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcdinv(a::T, b::T) where T <: Integer
   g, s, t = gcdx(a, b)
   return g, s
end

###############################################################################
#
#   Square root
#
###############################################################################

const sqrt_moduli = [3, 5, 7, 8]
const sqrt_residues = [[0, 1], [0, 1, 4], [0, 1, 2, 4], [0, 1, 4]]

@doc raw"""
    sqrt(a::T; check::Bool=true) where T <: Integer

Return the square root of $a$. By default the function will throw an exception
if the input is not square. If `check=false` this test is omitted.
"""
function sqrt(a::T; check::Bool=true) where T <: Integer
   s = isqrt(a)
   (check && s*s != a) && error("Not a square in sqrt")
   return s
end

@doc raw"""
    is_square_with_sqrt(a::T) where T <: Integer

Return `(true, s)` if $a$ is a perfect square, where $s^2 = a$. Otherwise
return `(false, 0)`.
"""
function is_square_with_sqrt(a::T) where T <: Integer
   if a < 0
      return false, zero(T)
   end
   s = isqrt(a)
   if a == s*s
      return true, s
   else
      return false, zero(T)
   end
end

function is_square_with_sqrt(a::BigInt)
   if a < 0
      return false, zero(BigInt)
   end
   for i = 1:length(sqrt_moduli)
      res = mod(a, sqrt_moduli[i])
      if !(res in sqrt_residues[i])
         return false, zero(BigInt)
      end
   end
   z = BigInt()
   r = BigInt()
   ccall((:__gmpz_sqrtrem, :libgmp), Cint,
         (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), z, r, a)
   if iszero(r)
      return true, z
   else
      return false, zero(BigInt)
   end
end

function is_square(a::T) where T <: Integer
   if a < 0
      return false
   end
   s = isqrt(a)
   return a == s*s
end

function is_square(a::BigInt)
   if a < 0
      return false
   end
   return Bool(ccall((:__gmpz_perfect_square_p, :libgmp), Cint,
                     (Ref{BigInt},), a))
end

###############################################################################
#
#   Root
#
###############################################################################

function root(a::BigInt, n::Int; check::Bool=true)
    a < 0 && iseven(n) && throw(DomainError((a, n),
                      "Argument `a` must be positive if exponent `n` is even"))
    n <= 0 && throw(DomainError(n, "Exponent must be positive"))
    z = BigInt()
    exact = Bool(ccall((:__gmpz_root, :libgmp), Cint,
                  (Ref{BigInt}, Ref{BigInt}, Cint), z, a, n))
    check && !exact && error("Not a perfect n-th power (n = $n)")
    return z
end

@doc raw"""
    root(a::T, n::Int; check::Bool=true) where T <: Integer

Return the $n$-th root of $a$. If `check=true` the function will test if the
input was a perfect $n$-th power, otherwise an exception will be raised. We
require $n > 0$.
"""
function root(a::T, n::Int; check::Bool=true) where T <: Integer
   if n == 2
      a < 0 && throw(DomainError((a, n),
                      "Argument `a` must be positive if exponent `n` is even"))
      s = isqrt(a)
      exact = true
      if check
         r = a - s*s
         exact = iszero(r)
         !exact && error("Not a perfect n-th power (n = $n)")
      end
      return s
   else
      return T(root(BigInt(a), n; check=check))
   end
end

# Suppose a = b^n, and let p be a prime divisor of n, i.e. n=m*p. Let q
# be a prime such that p divides q-1. Then by Fermat's little theorem
# either q divides a and hence b, or else b^(q-1) \equiv 1 mod q holds.
# But we also have that a = b^n = b^(m*p) = (b^m)^p mod q. Thus in fact
# we have   a^{(q-1)/p} \equiv 1 mod q.
#
# So below are some tables for various values of p and q and the allowed
# residue classes for a. They could be recomputed via
#
#  [ a for a in 0:q-1 if a == 0 || powermod(a, divexact(q-1,p), q) == 1]

const _p_q_residues = (
    #p => [q1 => residues,  q2 => residues, .... ]
    3 => ( 7 => [0, 1, 6], 13 => [0, 1, 5, 8, 12] ),
    5 => ( 11 => [0, 1, 10], 31 => [0, 1, 5, 6, 25, 26, 30] ),
    7 => ( 29 => [0, 1, 12, 17, 28], 43 => [0, 1, 6, 7, 36, 37, 42] ),
  )

# helper which returns false if a definitely is not an n-th power,
# otherwise return true (indicating we don't know)
function ispower_moduli(a::Integer, n::Int)
   @assert n >= 2
   n == 2 && return true

   # If a is even and an n-th power then it must be divisible by 2^n
   # and hence by 8 (as n >= 3 at this point).
   if iseven(a) && mod(a, 8) != 0
      return false
   end

   for (p, q_residues) in _p_q_residues
     if is_divisible_by(n, p)
        for (q, residues) in q_residues
           if !(mod(a, q) in residues)
              return false
           end
        end
     end
   end

   return true
end

@doc raw"""
    is_power(a::T, n::Int) where T <: Integer

Return `true, q` if $a$ is a perfect $n$-th power with $a = q^n$. Otherwise
return `false, 0`. We require $n > 0$.
"""
function is_power(a::T, n::Int) where T <: Integer
   n <= 0 && throw(DomainError(n, "exponent n must be positive"))
   if n == 1 || iszero(a) || isone(a)
      return (true, a)
   elseif a == -1
      return isodd(n) ? (true, a) : (false, zero(T))
   elseif iseven(n) && a < 0
      return false, zero(T)
   elseif !ispower_moduli(a, n)
      return (false, zero(T))
   end
      
   q = BigInt()
   r = BigInt()
   ccall((:__gmpz_rootrem, :libgmp), Nothing,
                     (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, Int), q, r, a, n)
   return iszero(r) ? (true, T(q)) : (false, zero(T))
end

function iroot(a::BigInt, n::Int)
    a < 0 && iseven(n) && throw(DomainError((a, n),
                      "Argument `a` must be positive if exponent `n` is even"))
    n <= 0 && throw(DomainError(n, "Exponent must be positive"))
    z = BigInt()
    ccall((:__gmpz_root, :libgmp), Cint,
                  (Ref{BigInt}, Ref{BigInt}, Cint), z, a, n)
    return z
end

@doc raw"""
    iroot(a::T, n::Int) where T <: Integer

Return the truncated integer part of the $n$-th root of $a$ (round towards
zero). We require $n > 0$ and also $a \geq 0$ if $n$ is even.
"""
function iroot(a::T, n::Int) where T <: Integer
   if n == 2
       a < 0 && throw(DomainError((a, n),
                      "Argument `a` must be positive if exponent `n` is even"))
       return isqrt(a)
   end
   return T(iroot(BigInt(a), n))
end

###############################################################################
#
#   Exponential
#
###############################################################################

@doc raw"""
    exp(a::T) where T <: Integer

Return $1$ if $a = 0$, otherwise throw an exception. This function is not
generally of use to the user, but is used internally in AbstractAlgebra.jl.
"""
function exp(a::T) where T <: Integer
    a != 0 && throw(DomainError(a, "a must be 0"))
    return T(1)
end

@doc raw"""
    log(a::T) where T <: Integer

Return $0$ if $a = 1$, otherwise throw an exception. This function is not
generally of use to the user, but is used internally in AbstractAlgebra.jl.
"""
function log(a::T) where T <: Integer
    a != 1 && throw(DomainError(a, "a must be 1"))
    return T(0)
end

###############################################################################
#
#   Coprime bases
#
###############################################################################

# Bernstein, "Factoring into coprimes in essentially linear time"
# ppio(a,b) = (c,n) where v_p(c) = v_p(a) if v_p(b) != 0, 0 otherwise
# c*n = a or c = gcd(a, b^infty), n = div(a, c).
# This is used in various Euclidean domains for Chinese remaindering.

@doc raw"""
    ppio(a::T, b::T)

Return a pair $(c,d)$ such that $a=c*d$ and $c = gcd(a, b^\infty)$ if $a\neq 0$,
and $c=b$, $d=0$ if $a=0$.
"""
function ppio(a::T, b::T) where T <: Integer
   iszero(a) && return (b,T(0))
   c = gcd(a, b)
   n = div(a, c)
   g = gcd(c, n)
   while !isone(g)
      c *= g
      n = div(n, g)
      g = gcd(c, n)
   end
   return c, n
end

###############################################################################
#
#   Primality test
#
###############################################################################

function is_probable_prime(x::Integer, reps::Integer=25)
   return ccall((:__gmpz_probab_prime_p, :libgmp), Cint,
                (Ref{BigInt}, Cint), x, reps) != 0
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

# No actual mutation is permitted for Julia types
# See #1077

function addmul!(a::T, b::T, c::T) where T <: Integer # special case, no temporary required
   return a + b*c
end

function neg!(w::Vector{<:Integer})
    return w .*= -1
end

###############################################################################
#
#   Random generation
#
###############################################################################

RandomExtensions.maketype(R::Integers{T}, _) where {T} = T

# define rand(make(ZZ, n:m))
rand(rng::AbstractRNG,
     sp::SamplerTrivial{<:Make2{T, Integers{T}, <:AbstractArray{<:Integer}}}
     ) where {T} =
        sp[][1](rand(rng, sp[][2]))


rand(rng::AbstractRNG, R::Integers, n) = R(rand(rng, n))

rand(R::Integers, n) = rand(Random.default_rng(), R, n)

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::Integers{T}) where {T <: Signed}
  n = T(2)^rand((1,1,1,2,3,10,31,32,33,63,64,65,100))
  return rand(R, -n:n)
end

function ConformanceTests.generate_element(R::Integers{T}) where {T <: Unsigned}
  n = T(2)^rand((1,1,1,2,3,10,31,32,33,63,64,65,100))
  return rand(R, 0:n)
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::Integers{T})() where T <: Integer
   return T(0)
end

function (a::Integers{T})(b::Union{Integer, Rational}) where T <: Integer
   return T(b)
end
