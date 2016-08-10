###############################################################################
#
#   Poly.jl : Generic multivariate polynomials over rings
#
###############################################################################

export GenMPoly, GenMPolyRing

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type{T, S, N}(::Type{GenMPoly{T, S, N}}) = GenMPolyRing{T, S, N}

elem_type{T <: RingElem, S, N}(::GenMPolyRing{T, S, N}) = GenMPoly{T, S, N}

vars(a::GenMPolyRing) = a.S

function gens{T <:RingElem, S, N}(a::GenMPolyRing{T, S, N})
   if (S == :deglex || S == :degrevlex)
      return [a([base_ring(a)(1)], [Monomial{S, N}(tuple(1, [Int(i == j) for j in 1:a.num_vars]...))])
           for i in 1:a.num_vars]
   else 
      return [a([base_ring(a)(1)], [Monomial{S, N}(tuple([Int(i == j) for j in 1:a.num_vars]...))])
           for i in 1:a.num_vars]
   end
end

###############################################################################
#
#   Monomial operations
#
###############################################################################

zero(::Type{Tuple{Int}}) = (0,)

zero{N}(::Type{NTuple{N, Int}}) = (0, zero(NTuple{N - 1, Int})...)

function +(a::Tuple{Int}, b::Tuple{Int})
    return (getfield(a, 1) + getfield(b, 1),)
end

function +{N}(a::NTuple{N, Int}, b::NTuple{N, Int})
   return (getfield(a, 1) + getfield(b, 1), (Base.tail(a) + Base.tail(b))...)
end

function *(a::Tuple{Int}, n::Int)
   return (getfield(a, 1)*n,)
end

function *{N}(a::NTuple{N, Int}, n::Int)
   return (getfield(a, 1)*n, (Base.tail(a)*n)...)
end

zero{S, N}(::Type{Monomial{S, N}}) = Monomial{S, N}(zero(NTuple{N, Int}))

function +{S, N}(a::Monomial{S, N}, b::Monomial{S, N})
   return Monomial{S, N}(a.exps + b.exps)
end

function isless{N}(a::Monomial{:lex, N}, b::Monomial{:lex, N})
   return a.exps < b.exps
end

function isless{N}(a::Monomial{:deglex, N}, b::Monomial{:deglex, N})
   return a.exps < b.exps
end

function isless{N}(a::Monomial{:revlex, N}, b::Monomial{:revlex, N})
   return reverse(a.exps) < reverse(b.exps)
end

function isless{N}(a::Monomial{:degrevlex, N}, b::Monomial{:degrevlex, N})
   return (getfield(a.exps, 1) < getfield(b.exps, 1)) ||
          (getfield(a.exps, 1) == getfield(b.exps, 1) && 
           reverse(Base.tail(a.exps)) < reverse(Base.tail(b.exps)))
end

function =={S, N}(a::Monomial{S, N}, b::Monomial{S, N})
   return a.exps == b.exps
end

function *{S, N}(a::Monomial{S, N}, n::Int)
   return Monomial{S, N}(a.exps*n)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function coeff(x::GenMPoly, i::Int)
   i < 0 && throw(DomainError())
   return x.coeffs[i + 1]
end

num_vars(x::GenMPoly) = parent(x).num_vars

function normalise(a::GenMPoly, n::Int)
   while n > 0 && iszero(a.coeffs[n]) 
      n -= 1
   end
   return n
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show{T <: RingElem, S, N}(io::IO, x::GenMPoly{T, S, N})
    len = length(x)
    U = [string(x) for x in vars(parent(x))]
    if len == 0
      print(io, base_ring(x)(0))
    else
      for i = 1:len
        c = coeff(x, len - i)
        bracket = needs_parentheses(c)
        if i != 1 && !is_negative(c)
          print(io, "+")
        end
        X = x.exps[len - i + 1].exps
        if !isone(c) && (c != -1 || show_minus_one(typeof(c)))
          if bracket
            print(io, "(")
          end
          show(io, c)
          if bracket
            print(io, ")")
          end
          if c != 1 && !(c == -1 && !show_minus_one(typeof(c))) && X != zero(NTuple{N, Int})
             print(io, "*")
          end
        end
        if c == -1 && !show_minus_one(typeof(c))
          print(io, "-")
        end
        d = (S == :deglex || S == :degrevlex) ? 1 : 0
        if X == zero(NTuple{N, Int})
          if c == 1
             print(io, c)
          elseif c == -1 && !show_minus_one(typeof(c))
             print(io, 1)
          end
        end
        fst = true
        for j = 1:num_vars(x)
          n = X[j + d]
          if n != 0
            if fst
               print(io, U[j])
               fst = false
            else
               print(io, "*", U[j])
            end
            if n != 1
              print(io, "^", n)
            end
          end
        end      
    end
  end
end

function show(io::IO, p::GenMPolyRing)
   const max_vars = 5 # largest number of variables to print
   n = p.num_vars
   print(io, "Multivariate Polynomial Ring in ")
   if n > max_vars
      print(io, p.num_vars)
      print(io, " variables ")
   end
   for i = 1:min(n - 1, max_vars - 1)
      print(io, string(p.S[i]), ", ")
   end
   if n > max_vars
      print(io, "..., ")
   end
   print(io, string(p.S[n]))
   print(io, " over ")
   show(io, base_ring(p))
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function -{T <: RingElem, S, N}(a::GenMPoly{T, S, N})
   r = parent(a)()
   fit!(r, length(a))
   for i = 1:length(a)
      r.coeffs[i] = -a.coeffs[i]
      r.exps[i] = a.exps[i]
   end
   r.length = a.length
   return r
end

function +{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, b::GenMPoly{T, S, N})
   r = parent(a)()
   fit!(r, length(a) + length(b))
   i = 1
   j = 1
   k = 1
   while i <= length(a) && j <= length(b)
      if a.exps[i] < b.exps[j]
         r.coeffs[k] = a.coeffs[i]
         r.exps[k] = a.exps[i]
          i += 1
      elseif a.exps[i] == b.exps[j]
         c = a.coeffs[i] + b.coeffs[j]
         if c != 0
            r.coeffs[k] = c
            r.exps[k] = a.exps[i]
         else
            k -= 1
         end
         i += 1
         j += 1
      else
         r.coeffs[k] = b.coeffs[j]
         r.exps[k] = b.exps[j]
         j += 1
      end
      k += 1
   end
   while i <= length(a)
      r.coeffs[k] = a.coeffs[i]
      r.exps[k] = a.exps[i]
      i += 1
      k += 1
   end
   while j <= length(b)
      r.coeffs[k] = b.coeffs[j]
      r.exps[k] = b.exps[j]
      j += 1
      k += 1
   end
   r.length = k - 1
   return r
end

function -{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, b::GenMPoly{T, S, N})
   r = parent(a)()
   fit!(r, length(a) + length(b))
   i = 1
   j = 1
   k = 1
   while i <= length(a) && j <= length(b)
      if a.exps[i] < b.exps[j]
         r.coeffs[k] = a.coeffs[i]
         r.exps[k] = a.exps[i]
          i += 1
      elseif a.exps[i] == b.exps[j]
         c = a.coeffs[i] - b.coeffs[j]
         if c != 0
            r.coeffs[k] = c
            r.exps[k] = a.exps[i]
         else
            k -= 1
         end
         i += 1
         j += 1
      else
         r.coeffs[k] = -b.coeffs[j]
         r.exps[k] = b.exps[j]
         j += 1
      end
      k += 1
   end
   while i <= length(a)
      r.coeffs[k] = a.coeffs[i]
      r.exps[k] = a.exps[i]
      i += 1
      k += 1
   end
   while j <= length(b)
      r.coeffs[k] = -b.coeffs[j]
      r.exps[k] = b.exps[j]
      j += 1
      k += 1
   end
   r.length = k - 1
   return r
end

function *{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, b::GenMPoly{T, S, N})
   m = length(a)
   n = length(b)
   if m == 0
      return parent(a)()
   end
   Ai = Array(Int, m)
   Ac = Array(T, m*n)
   Bc = Array(T, m*n)
   Ae = Array(Monomial{S, N}, m*n)
   Be = Array(Monomial{S, N}, m*n)
   for i = 1:m
      c = a.coeffs[i]
      d = a.exps[i].exps
      k = 1
      for j = 1:n
         s = Ac[(i-1)*n + k] = c*b.coeffs[j]
         if s != 0
            Ae[(i-1)*n + k] = Monomial{S, N}(d + b.exps[j].exps)
            k += 1
         end
      end
      Ai[i] = k - 1
   end
   while m > 1
      m2 = div(m, 2)
      for t = 1:m2
         i = 1
         j = 1
         k = 1
         s1 = n*(2t - 2)
         s2 = n*(2t - 1)
         r = 2*n*(t - 1)
         while i <= Ai[2t - 1] && j <= Ai[2t]
            if Ae[s1 + i] < Ae[s2 + j]
               Bc[r + k] = Ac[s1 + i]
               Be[r + k] = Ae[s1 + i]
               i += 1
            elseif Ae[s1 + i] == Ae[s2 + j]
               c = Ac[s1 + i] + Ac[s2 + j]
               if c != 0
                  Bc[r + k] = c
                  Be[r + k] = Ae[s1 + i]
               else
                  k -= 1
               end
               i += 1
               j += 1
            else
               Bc[r + k] = Ac[s2 + j]
               Be[r + k] = Ae[s2 + j]
               j += 1
            end
            k += 1
         end
         while i <= Ai[2t - 1]
            Bc[r + k] = Ac[s1 + i]
            Be[r + k] = Ae[s1 + i]
            i += 1
            k += 1
         end
         while j <= Ai[2t]
            Bc[r + k] = Ac[s2 + j]
            Be[r + k] = Ae[s2 + j]
            j += 1
            k += 1
         end
         Ai[t] = k - 1
      end
      if isodd(m)
         s1 = n*(2m2)
         r = 2*n*m2
         for i = 1:Ai[2m2 + 1]
            Bc[r + i] = Ac[s1 + i]
            Be[r + i] = Ae[s1 + i]
         end
         Ai[m2 + 1] = Ai[2m2 + 1]
      end
      n = 2*n
      m = div(m + 1, 2) 
      t1 = Ac
      Ac = Bc
      Bc = t1
      t2 = Ae
      Ae = Be
      Be = t2        
   end
   resize!(Ac, Ai[1])
   resize!(Ae, Ai[1])
   return parent(a)(Ac, Ae)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

function *{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, n::Integer)
   r = parent(a)()
   fit!(r, length(a))
   j = 1
   for i = 1:length(a)
      c = a.coeffs[i]*n
      if c != 0
         r.coeffs[j] = c 
         r.exps[j] = a.exps[i]
         j += 1
      end
   end
   r.length = j - 1
   return r
end

function *{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, n::fmpz)
   r = parent(a)()
   fit!(r, length(a))
   j = 1
   for i = 1:length(a)
      c = a.coeffs[i]*n
      if c != 0
         r.coeffs[j] = c 
         r.exps[j] = a.exps[i]
         j += 1
      end
   end
   r.length = j - 1
   return r
end

function *{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, n::T)
   r = parent(a)()
   fit!(r, length(a))
   j = 1
   for i = 1:length(a)
      c = a.coeffs[i]*n
      if c != 0
         r.coeffs[j] = c 
         r.exps[j] = a.exps[i]
         j += 1
      end
   end
   r.length = j - 1
   return r
end

*{T <: RingElem, S, N}(n::Integer, a::GenMPoly{T, S, N}) = a*n

*{T <: RingElem, S, N}(n::fmpz, a::GenMPoly{T, S, N}) = a*n

*{T <: RingElem, S, N}(n::T, a::GenMPoly{T, S, N}) = a*n

+{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, b::T) = a + parent(a)(b)

+{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, b::Integer) = a + parent(a)(b)

+{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, b::fmpz) = a + parent(a)(b)

-{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, b::T) = a - parent(a)(b)

-{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, b::Integer) = a - parent(a)(b)

-{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, b::fmpz) = a - parent(a)(b)

+{T <: RingElem, S, N}(a::T, b::GenMPoly{T, S, N}) = parent(b)(a) + b

+{T <: RingElem, S, N}(a::Integer, b::GenMPoly{T, S, N}) = parent(b)(a) + b

+{T <: RingElem, S, N}(a::fmpz, b::GenMPoly{T, S, N}) = parent(b)(a) + b

-{T <: RingElem, S, N}(a::T, b::GenMPoly{T, S, N}) = parent(b)(a) - b

-{T <: RingElem, S, N}(a::Integer, b::GenMPoly{T, S, N}) = parent(b)(a) - b

-{T <: RingElem, S, N}(a::fmpz, b::GenMPoly{T, S, N}) = parent(b)(a) - b

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::GenMPoly, b::Int)
   b < 0 && throw(DomainError())
   # special case powers of x for constructing polynomials efficiently
   if length(a) == 0
      return parent(a)()
   elseif length(a) == 1
      return parent(a)([coeff(a, 0)^b], [a.exps[1]*b])
   elseif b == 0
      return parent(a)(1)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit != 0
         z = z*z
         if (UInt(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function fit!{T <: RingElem, S, N}(a::GenMPoly{T, S, N}, n::Int)
   if length(a.coeffs) < n
      resize!(a.coeffs, n)
      resize!(a.exps, n)
   end
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function Base.call{T <: RingElem, S, N}(a::GenMPolyRing{T, S, N}, b::RingElem)
   return a(base_ring(a)(b), a.vars)
end

function Base.call{T <: RingElem, S, N}(a::GenMPolyRing{T, S, N})
   z = GenMPoly{T, S, N}()
   z.parent = a
   return z
end

function Base.call{T <: RingElem, S, N}(a::GenMPolyRing{T, S, N}, b::Integer)
   z = GenMPoly{T, S, N}(base_ring(a)(b))
   z.parent = a
   return z
end

function Base.call{T <: RingElem, S, N}(a::GenMPolyRing{T, S, N}, b::T)
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = GenMPoly{T, S, N}(b)
   z.parent = a
   return z
end

function Base.call{T <: RingElem, S, N}(a::GenMPolyRing{T, S, N}, b::PolyElem{T})
   parent(b) != a && error("Unable to coerce polynomial")
   return b
end

function Base.call{T <: RingElem, S, N}(a::GenMPolyRing{T, S, N}, b::Array{T, 1}, m::Array{Monomial{S, N}, 1})
   if length(b) > 0
      parent(b[1]) != base_ring(a) && error("Unable to coerce to polynomial")
   end
   z = GenMPoly{T, S, N}(b, m)
   z.parent = a
   return z
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

doc"""
    PolynomialRing(R::Ring, s::Array{String, 1}; cached::Bool = true, S::Symbol = :lex)
> Given a base ring `R` and an array of strings `s` specifying how the
> generators (variables) should be printed, return a tuple `S, x1, x2, ...`
> representing the new polynomial ring $T = R[x1, x2, ...]$ and the generators
> $x1, x2, ...$ of the polynomial ring. By default the parent object `T` will
> depend only on `R` and `x1, x2, ...` and will be cached. Setting the optional
> argument `cached` to `false` will prevent the parent object `T` from being
> cached. `S` is a symbol corresponding to the ordering of the polynomial and
> can be one of `:lex`, `:deglex`, `:revlex` or `:degrevlex`.
"""
function PolynomialRing(R::Ring, s::Array{String, 1}; cached::Bool = true, ordering::Symbol = :lex)
   U = [Symbol(x) for x in s]
   T = elem_type(R)
   N = (ordering == :deglex || ordering == :degrevlex) ? length(U) + 1 : length(U)
   parent_obj = GenMPolyRing{T, ordering, N}(R, U, cached)

   return tuple(parent_obj, gens(parent_obj)...)
end