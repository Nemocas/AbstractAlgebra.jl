###############################################################################
#
#   SRPoly.jl : Generic recursive sparse multivariate polynomials over rings
#
###############################################################################

export GenSRPoly, GenSRPolyRing, SRPolynomialRing

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type{T}(::Type{GenSRPoly{T}}) = GenSRPolyRing{T}

elem_type{T <: RingElem}(::GenSRPolyRing{T}) = GenSRPoly{T}

var(a::GenSRPolyRing) = a.S

function gen{T <:RingElem}(a::GenSRPolyRing{T})
   return a([one(base_ring(a))], [UInt(1)])
end

function gens{T <: RingElem}(a::GenSRPolyRing{T})
   A = Array(GenSRPoly, a.num_vars)
   R = a
   for i = 1:a.num_vars
      A[i] = gen(R)
      R = base_ring(R)
   end
   return A
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function coeff(x::GenSRPoly, i::Int)
   i < 0 && throw(DomainError())
   return x.coeffs[i + 1]
end

num_vars(x::GenSRPoly) = parent(x).num_vars

function normalise(a::GenSRPoly, n::Int)
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

function show{T <: RingElem}(io::IO, x::GenSRPoly{T})
    len = length(x)
    U = string(var(parent(x)))
    if len == 0
      print(io, base_ring(x)(0))
    else
      for i = 1:len
        c = coeff(x, len - i)
        bracket = needs_parentheses(c)
        if i != 1 && !is_negative(c)
          print(io, "+")
        end
        X = x.exps[len - i + 1]
        if !isone(c) && (c != -1 || show_minus_one(typeof(c)))
          if bracket
            print(io, "(")
          end
          show(io, c)
          if bracket
            print(io, ")")
          end
          if c != 1 && !(c == -1 && !show_minus_one(typeof(c))) && X != UInt(0)
             print(io, "*")
          end
        end
        if c == -1 && !show_minus_one(typeof(c))
          print(io, "-")
        end
        if X == UInt(0)
          if c == 1
             print(io, c)
          elseif c == -1 && !show_minus_one(typeof(c))
             print(io, 1)
          end
        end
        n = reinterpret(Int, X)
        if n != 0
           print(io, U)
        end
        if n > 1
           print(io, "^", n)
        end
      end
    end
end

function show(io::IO, p::GenSRPolyRing)
   const max_vars = 5 # largest number of variables to print
   n = p.num_vars
   print(io, "Recursive Multivariate Polynomial Ring in ")
   if n > max_vars
      print(io, p.num_vars)
      print(io, " variables ")
   end
   for i = 1:min(n - 1, max_vars - 1)
      print(io, string(p.S), ", ")
      p = base_ring(p)
   end
   if n > max_vars
      print(io, "..., ")
   end
   print(io, string(p.S))
   print(io, " over ")
   show(io, base_ring(p))
end

show_minus_one{T <: RingElem}(::Type{GenSRPoly{T}}) = show_minus_one(T)

needs_parentheses{T <: RingElem}(a::GenSRPoly{T}) = length(a) > 1

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function -{T <: RingElem}(a::GenSRPoly{T})
   r = parent(a)()
   fit!(r, length(a))
   for i = 1:length(a)
      r.coeffs[i] = -a.coeffs[i]
      r.exps[i] = a.exps[i]
   end
   r.length = a.length
   return r
end

function +{T <: RingElem}(a::GenSRPoly{T}, b::GenSRPoly{T})
   par = parent(a)
   r = par()
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

function -{T <: RingElem}(a::GenSRPoly{T}, b::GenSRPoly{T})
   par = parent(a)
   r = par()
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

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function =={T <: RingElem}(a::GenSRPoly{T}, b::Int)
   return length(a) == 0 ? b == 0 : a.length == 1 && 
          a.exps[1] == 0 && a.coeffs[1] == b
end

=={T <: RingElem}(a::Int, b::GenSRPoly{T}) = b == a

function =={T <: RingElem}(a::GenSRPoly{T}, b::fmpz)
   return length(a) == 0 ? b == 0 : a.length == 1 &
          a.exps[1] == 0 && a.coeffs[1] == b
end

=={T <: RingElem}(a::fmpz, b::GenSRPoly{T}) = b == a

function =={T <: RingElem}(a::GenSRPoly{T}, b::T)
   return length(a) == 0 ? b == 0 : a.length == 1 &&
          a.exps[1] == 0 && a.coeffs[1] == b
end

=={T <: RingElem}(a::T, b::GenSRPoly{T}) = b == a

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function fit!{T <: RingElem}(a::GenSRPoly{T}, n::Int)
   if length(a.coeffs) < n
      resize!(a.coeffs, n)
      resize!(a.exps, n)
   end
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

Base.promote_rule{T <: RingElem, V <: Integer}(::Type{GenSRPoly{T}}, ::Type{V}) = GenSRPoly{T}

Base.promote_rule{T <: RingElem}(::Type{GenSRPoly{T}}, ::Type{T}) = GenSRPoly{T}

function promote_rule1{T <: RingElem, U <: RingElem}(::Type{GenSRPoly{T}}, ::Type{GenSRPoly{U}})
   Base.promote_rule(T, GenSRPoly{U}) == T ? GenSRPoly{T} : Union{}
end

function Base.promote_rule{T <: RingElem, U <: RingElem}(::Type{GenSRPoly{T}}, ::Type{U})
   Base.promote_rule(T, U) == T ? GenSRPoly{T} : promote_rule1(U, GenSRPoly{T})
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function Base.call{T <: RingElem}(a::GenSRPolyRing{T}, b::RingElem)
   return a([base_ring(a)(b)], [UInt(0)])
end

function Base.call{T <: RingElem}(a::GenSRPolyRing{T})
   z = GenSRPoly{T}()
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::GenSRPolyRing{T}, b::Integer)
   z = GenSRPoly{T}(base_ring(a)(b))
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::GenSRPolyRing{T}, b::T)
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = GenSRPoly{T}(b)
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::GenSRPolyRing{T}, b::PolyElem{T})
   parent(b) != a && error("Unable to coerce polynomial")
   return b
end

function Base.call{T <: RingElem}(a::GenSRPolyRing{T}, b::Array{T, 1}, m::Array{UInt, 1})
   if length(b) > 0
      parent(b[1]) != base_ring(a) && error("Unable to coerce to polynomial")
   end
   z = GenSRPoly{T}(b, m)
   z.parent = a
   return z
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

doc"""
    SRPolynomialRing(R::Ring, s::Array{String, 1}; cached::Bool = true)
> Given a base ring `R` and an array of strings `s` specifying how the
> generators (variables) should be printed, return a tuple `S, x1, x2, ...`
> representing the new polynomial ring $T = R[x1, x2, ...]$ and the generators
> $x1, x2, ...$ of the polynomial ring. By default the parent object `T` will
> depend only on `R` and `x1, x2, ...` and will be cached. Setting the optional
> argument `cached` to `false` will prevent the parent object `T` from being
> cached.
"""
function SRPolynomialRing(R::Ring, s::Array{ASCIIString{}, 1}; cached::Bool = true)
   U = [Symbol(s[i]) for i in length(s):-1:1]
   T = elem_type(R)

   R1 = R
   T1 = elem_type(R)
   for i = 1:length(U)
      R1 = GenSRPolyRing{T1}(R1, U[i], i, cached)
      T1 = elem_type(R1)
   end

   return tuple(R1, gens(R1)...)
end
