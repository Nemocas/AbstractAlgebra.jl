###############################################################################
#
#   Integer.jl : Additional AbstractAlgebra functionality for Julia Integer
#
###############################################################################

export iroot, ispower, ispower_with_root, root, issquare_with_sqrt,
       isprobable_prime

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

base_ring(a::Integers{T}) where T <: Integer = Union{}

isexact_type(::Type{T}) where T <: Integer = true

isdomain_type(::Type{T}) where T <: Integer = true

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(::Integers{T}) where T <: Integer = T(0)

one(::Integers{T}) where T <: Integer = T(1)

isunit(a::Integer) = a == 1 || a == -1

canonical_unit(a::T) where T <: Integer = a < 0 ? T(-1) : T(1)

characteristic(::Integers{T}) where T <: Integer = 0

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::Integer; context = nothing)
    return a
end

function show(io::IO, R::Integers)
   print(io, "Integers")
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
   if b == 0
      return a == 0, b
   end
   q, r = divrem(a, b)
   return r == 0, q
end

@doc Markdown.doc"""
    isdivisible_by(a::Integer, b::Integer)

Return `true` if $a$ is divisible by $b$, i.e. if there exists $c$ such that
$a = bc$.
"""
function isdivisible_by(a::Integer, b::Integer)
   if b == 0
      return a == 0
   end
   r = rem(a, b)
   return r == 0
end

function isdivisible_by(a::BigInt, b::BigInt)
   if b == 0
      return a == 0
   end
   return Bool(ccall((:__gmpz_divisible_p, :libgmp), Cint,
                                             (Ref{BigInt}, Ref{BigInt}), a, b))
end

function isdivisible_by(a::BigInt, b::Int)
   if b == 0
      return a == 0
   end
   return Bool(ccall((:__gmpz_divisible_ui_p, :libgmp), Cint,
                                        (Ref{BigInt}, Int), a, b < 0 ? -b : b))
end

function isdivisible_by(a::BigInt, b::UInt)
   if b == 0
      return a == 0
   end
   return Bool(ccall((:__gmpz_divisible_ui_p, :libgmp), Cint,
                                                    (Ref{BigInt}, UInt), a, b))
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::Integer, b::Integer; check::Bool=true)
   if check
      q, r = divrem(a, b)
      iszero(r) || throw(ArgumentError("Not an exact division"))
   else
      q = div(a, b)
   end
   return q
end

function divexact(a::BigInt, b::BigInt; check::Bool=true)
   q = BigInt()
   if check
      r = BigInt()
      ccall((:__gmpz_tdiv_qr, :libgmp), Nothing,
              (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), q, r, a, b)
      r != 0 && throw(ArgumentError("Not an exact division"))
   else
      ccall((:__gmpz_divexact, :libgmp), Nothing,
                              (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), q, a, b)
   end
   return q
end

function divexact(a::BigInt, b::Int; check::Bool=true)
   q = BigInt()
   sgn = b < 0
   if check
      r = BigInt()
      ccall((:__gmpz_tdiv_qr_ui, :libgmp), Nothing,
           (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, Int), q, r, a, sgn ? -b : b)
      r != 0 && throw(ArgumentError("Not an exact division"))
   else
      ccall((:__gmpz_divexact_ui, :libgmp), Nothing,
                           (Ref{BigInt}, Ref{BigInt}, Int), q, a, sgn ? -b : b)
   end
   return sgn ? -q : q
end

function divexact(a::BigInt, b::UInt; check::Bool=true)
   q = BigInt()
   if check
      r = BigInt()
      ccall((:__gmpz_tdiv_qr_ui, :libgmp), Nothing,
                     (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, UInt), q, r, a, b)
      r != 0 && throw(ArgumentError("Not an exact division"))
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

sqrt_moduli = [3, 5, 7, 8]
sqrt_residues = [[0, 1], [0, 1, 4], [0, 1, 2, 4], [0, 1, 4]]

@doc Markdown.doc"""
    sqrt(a::T; check::Bool=true) where T <: Integer

Return the integer square root of $a$. By default the function will
throw an exception if the input is not square. If `check=false` this test is
omitted.
"""
function sqrt(a::T; check::Bool=true) where T <: Integer
   s = isqrt(a)
   (check && s*s != a) && error("Not a square in sqrt")
   return s
end

@doc Markdown.doc"""
    issquare_with_sqrt(a::T) where T <: Integer

Return `(true, s)` if $a$ is a perfect square, where $s^2 = a$. Otherwise
return `(false, 0)`.
"""
function issquare_with_sqrt(a::T) where T <: Integer
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

function issquare_with_sqrt(a::BigInt)
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


@doc Markdown.doc"""
    issquare(a::T) where T <: Integer

Return true if $a$ is a square.
"""
function issquare(a::T) where T <: Integer
   if a < 0
      return false
   end
   s = isqrt(a)
   return a == s*s
end

function issquare(a::BigInt)
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

@doc Markdown.doc"""
    root(a::T, n::Int; check::Bool=true) where T <: Integer

Return the $n$-th root of $a$. If `check=true` the function will test if the
input was a perfect $n$-th power, otherwise an exception will be raised. We
require $n > 0$ and also $a \geq 0$ if $n$ is even.
"""
function root(a::T, n::Int; check::Bool=true) where T <: Integer
   if n == 2
      a < 0 && throw(DomainError((a, n),
                      "Argument `a` must be positive if exponent `n` is even"))
      s = isqrt(a)
      exact = true
      if check
         r = a - s*s
         exact = r == 0
         !exact && error("Not a perfect n-th power (n = $n)")
      end
      return s
   else
      return T(root(BigInt(a), n; check=check))
   end
end

moduli3 = [7, 8, 13]
residues3 = [[0, 1, 6], [0, 1, 3, 5, 7], [0, 1, 5, 8, 12]]

moduli5 = [8, 11, 31]
residues5 = [[0, 1, 3, 5, 7], [0, 1, 10], [0, 1, 5, 6, 25, 26, 30]]

moduli7 = [8, 29, 43]
residues7 = [[0, 1, 3, 5, 7], [0, 1, 12, 17, 28], [0, 1, 6, 7, 36, 37, 42]]

function ispower_moduli(a::Integer, n::Int)
   if mod(n, 3) == 0
      for i = 1:length(moduli3)
         if !(mod(a, moduli3[i]) in residues3[i])
            return false
         end
      end
   elseif (n % 5) == 0
      for i = 1:length(moduli5)
         if !(mod(a, moduli5[i]) in residues5[i])
            return false
         end
      end
   elseif (n % 3) == 0
      for i = 1:length(moduli7)
         if !(mod(a, moduli7[i]) in residues7[i])
            return false
         end
      end
   elseif isodd(n)
      if !(mod(a, moduli5[1]) in residues5[1])
            return false
      end
   end
   return true
end

@doc Markdown.doc"""
    ispower_with_root(a::T, n::Int) where T <: Integer

Return `true, q` if $a$ is a perfect $n$-th power with $a = q^n$. Otherwise
return `false, 0`. We require $n > 0$.
"""
function ispower_with_root(a::T, n::Int) where T <: Integer
   n <= 0 && throw(DomainError(n, "exponent n must be positive"))
   if n == 1 || a == 0 || a == 1
      return (true, a)
   elseif a == -1
      return isodd(n) ? (true, a) : (false, zero(T))
   elseif mod(n, 2) == 0 && a < 0
      return false, zero(BigInt)
   elseif !ispower_moduli(a, n)
      return (false, zero(BigInt))
   end
      
   q = BigInt()
   r = BigInt()
   ccall((:__gmpz_rootrem, :libgmp), Nothing,
                     (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, Int), q, r, a, n)
   return iszero(r) ? (true, T(q)) : (false, zero(T))
end

@doc Markdown.doc"""
    ispower(a::T, n::Int) where T <: Integer

Return `true` if $a$ is a perfect $n$-th power, i.e. if there is an integer $b$
such that $a = b^n$. We require $n > 0$.
"""
function ispower(a::T, n::Int) where T <: Integer
   n <= 0 && throw(DomainError(n, "n is not positive"))
   if n == 1 || a == 0 || a == 1
      return true
   elseif a == -1
      return isodd(n) ? true : false
   elseif mod(n, 2) == 0 && a < 0
      return false
   elseif !ispower_moduli(a, n)
      return false
   end
   
   q = BigInt()
   r = BigInt()
   ccall((:__gmpz_rootrem, :libgmp), Nothing,
                     (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, Int), q, r, BigInt(a), n)

   return r == 0
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

@doc Markdown.doc"""
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

@doc Markdown.doc"""
    exp(a::T) where T <: Integer

Return $1$ if $a = 0$, otherwise throw an exception. This function is not
generally of use to the user, but is used internally in AbstractAlgebra.jl.
"""
function exp(a::T) where T <: Integer
    a != 0 && throw(DomainError(a, "a must be 0"))
    return T(1)
end

@doc Markdown.doc"""
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

@doc Markdown.doc"""
    ppio(a::T, b::T)

Split $a$ into $c*d$ where $c = gcd(a, b^\infty)$.
"""
function ppio(a::T, b::T) where T <: Integer
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

function isprobable_prime(x::Integer, reps::Integer=25)
   return ccall((:__gmpz_probab_prime_p, :libgmp), Cint,
                (Ref{BigInt}, Cint), x, reps) != 0
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::Integer)
   return 0
end

function zero!(a::BigInt)
   ccall((:__gmpz_set_si, :libgmp), Nothing, (Ref{BigInt}, Int), a, 0)
   return a
end

function mul!(a::T, b::T, c::T) where T <: Integer
   return b*c
end

function mul!(a::BigInt, b::BigInt, c::BigInt)
   ccall((:__gmpz_mul, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

function add!(a::T, b::T, c::T) where T <: Integer
   return b + c
end

function add!(a::BigInt, b::BigInt, c::BigInt)
   ccall((:__gmpz_add, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

function addeq!(a::T, b::T) where T <: Integer
   return a + b
end

function addeq!(a::BigInt, b::BigInt)
   ccall((:__gmpz_add, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, a, b)
   return a
end

function addmul!(a::T, b::T, c::T, d::T) where T <: Integer
   return a + b*c
end

function addmul!(a::BigInt, b::BigInt, c::BigInt, d::BigInt)
   ccall((:__gmpz_addmul, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

function addmul!(a::T, b::T, c::T) where T <: Integer # special case, no temporary required
   return a + b*c
end

function addmul!(a::BigInt, b::BigInt, c::BigInt) # special case, no temporary required
   ccall((:__gmpz_addmul, :libgmp), Nothing, (Ref{BigInt}, Ref{BigInt}, Ref{BigInt}), a, b, c)
   return a
end

###############################################################################
#
#   Random generation
#
###############################################################################

RandomExtensions.maketype(R::AbstractAlgebra.Integers{T}, _) where {T} = T

# define rand(make(ZZ, n:m))
rand(rng::AbstractRNG,
     sp::SamplerTrivial{<:Make2{T, Integers{T}, <:AbstractArray{<:Integer}}}
     ) where {T} =
        sp[][1](rand(rng, sp[][2]))


rand(rng::AbstractRNG, R::Integers, n) = R(rand(rng, n))

rand(R::Integers, n) = rand(Random.GLOBAL_RNG, R, n)

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

