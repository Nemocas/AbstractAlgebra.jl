###############################################################################
#
#   Fraction.jl : fraction fields
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring_type(::Type{<:FracField{T}}) where T<:RingElement = parent_type(T)

function is_domain_type(::Type{T}) where {S <: RingElement, T <: FracElem{S}}
   return is_domain_type(S)
end

function is_exact_type(a::Type{T}) where {S <: RingElement, T <: FracElem{S}}
   return is_exact_type(S)
end

@doc raw"""
    characteristic(R::FracField{T}) where T <: RingElem

Return the characteristic of the given field.
"""
function characteristic(R::FracField{T}) where T <: RingElem
   return characteristic(base_ring(R))
end

@doc raw"""
    vars(a::FracElem{S}) where {S <: MPolyRingElem{<: RingElement}}

Return the variables actually occurring in $a$. Returned variables are elements
of `base_ring(a)`. The variables from the numerator go first.
"""
function vars(a::FracElem{S}) where {S <: MPolyRingElem{<: RingElement}}
   n = numerator(a, false)
   d = denominator(a, false)
   n_vars = vars(n)
   d_vars = vars(d)
   nd_vars = union!(n_vars, d_vars)
   return nd_vars
end

###############################################################################
#
#   Constructors
#
###############################################################################

function //(x::T, y::T) where {T <: RingElem}
   R = parent(x)
   iszero(y) && throw(DivideError())
   g = gcd(x, y)
   z = Generic.FracFieldElem{T}(divexact(x, g), divexact(y, g))
   try
      z.parent = Generic.FracDict[R]
   catch
      z.parent = Generic.fraction_field(R)
   end
   return z
end

//(x::T, y::FracElem{T}) where {T <: RingElem} = parent(y)(x)//y

//(x::FracElem{T}, y::T) where {T <: RingElem} = x//parent(x)(y)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::FracElem, h::UInt)
   b = 0x8a30b0d963237dd5%UInt
   # We canonicalise before hashing
   return xor(b, hash(numerator(a, true), h), hash(denominator(a, true), h), h)
end

# Fall back method for all other fraction types in system
function Base.numerator(a::FracElem, canonicalise::Bool=true)
   return Base.numerator(a) # all other types ignore canonicalise
end

# Fall back method for all other fraction types in system
function Base.denominator(a::FracElem, canonicalise::Bool=true)
   return Base.denominator(a) # all other types ignore canonicalise
end

zero(R::FracField) = R(0)

one(R::FracField) = R(1)

iszero(a::FracElem) = iszero(numerator(a, false))

isone(a::FracElem) = numerator(a, false) == denominator(a, false)

is_unit(a::FracElem) = !iszero(numerator(a, false))

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::FracElem) = a

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function expressify(a::FracElem; context = nothing)
    n = numerator(a, true)
    d = denominator(a, true)
    if isone(d)
        return expressify(n; context)
    else
        return Expr(:call, ://, expressify(n; context), expressify(d; context))
    end
end

@enable_all_show_via_expressify FracElem

function show(io::IO, mime::MIME"text/plain", a::FracField)
  @show_name(io, a)
  @show_special(io, mime, a)

  println(io, "Fraction field")
  io = pretty(io)
  print(io, Indent(), "of ", Lowercase(), base_ring(a))
  print(io, Dedent())
end

function show(io::IO, a::FracField)
  @show_name(io, a)
  @show_special(io, a)
  if is_terse(io)
    print(io, "Fraction field")
  else
    io = pretty(io)
    print(io, "Fraction field of ")
    print(terse(io), Lowercase(), base_ring(a))
  end
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::FracElem)
   return parent(a)(-numerator(a, false), deepcopy(denominator(a, false)))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   d1 = denominator(a, false)
   d2 = denominator(b, false)
   n1 = numerator(a, false)
   n2 = numerator(b, false)
   if d1 == d2
      rnum = n1 + n2
      if isone(d1)
         rden = deepcopy(d1)
      else
         gd = gcd(rnum, d1)
         if isone(gd)
            rden = deepcopy(d1)
         else
            rnum = divexact(rnum, gd)
            rden = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      rnum = n1*d2 + n2
      rden = deepcopy(d2)
   elseif isone(d2)
      rnum = n1 + n2*d1
      rden = deepcopy(d1)
   else
      gd = gcd(d1, d2)
      if isone(gd)
         rnum = n1*d2 + n2*d1
         rden = d1*d2
      else
         q1 = divexact(d1, gd)
         q2 = divexact(d2, gd)
         rnum = q1*n2 + q2*n1
         t = gcd(rnum, gd)
         if isone(t)
            rden = q2*d1
         else
            rnum = divexact(rnum, t)
            gd = divexact(d1, t)
            rden = gd*q2
         end
      end
   end
   return parent(a)(rnum, rden)
end

function -(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   d1 = denominator(a, false)
   d2 = denominator(b, false)
   n1 = numerator(a, false)
   n2 = numerator(b, false)
   if d1 == d2
      rnum = n1 - n2
      if isone(d1)
         rden = deepcopy(d1)
      else
         gd = gcd(rnum, d1)
         if isone(gd)
            rden = deepcopy(d1)
         else
            rnum = divexact(rnum, gd)
            rden = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      rnum = n1*d2 - n2
      rden = deepcopy(d2)
   elseif isone(d2)
      rnum = n1 - n2*d1
      rden = deepcopy(d1)
   else
      gd = gcd(d1, d2)
      if isone(gd)
         rnum = n1*d2 - n2*d1
         rden = d1*d2
      else
         q1 = divexact(d1, gd)
         q2 = divexact(d2, gd)
         rnum = q2*n1 - q1*n2
         t = gcd(rnum, gd)
         if isone(t)
            rden = q2*d1
         else
            rnum = divexact(rnum, t)
            gd = divexact(d1, t)
            rden = gd*q2
         end
      end
   end
   return parent(a)(rnum, rden)
end

function *(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   n1 = numerator(a, false)
   d2 = denominator(b, false)
   n2 = numerator(b, false)
   d1 = denominator(a, false)
   if d1 == d2
      n = n1*n2
      d = d1*d2
   elseif isone(d1)
      gd = gcd(n1, d2)
      if isone(gd)
         n = n1*n2
         d = deepcopy(d2)
      else
         n = divexact(n1, gd)*n2
         d = divexact(d2, gd)
      end
   elseif isone(d2)
      gd = gcd(n2, d1)
      if isone(gd)
         n = n2*n1
         d = deepcopy(d1)
      else
         n = divexact(n2, gd)*n1
         d = divexact(d1, gd)
      end
   else
      g1 = gcd(n1, d2)
      g2 = gcd(n2, d1)
      if !isone(g1)
         n1 = divexact(n1, g1)
         d2 = divexact(d2, g1)
      end
      if !isone(g2)
         n2 = divexact(n2, g2)
         d1 = divexact(d1, g2)
      end
      n = n1*n2
      d = d1*d2
   end
   return parent(a)(n, d)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::FracElem, b::Union{Integer, AbstractFloat})
   c = base_ring(a)(b)
   g = gcd(denominator(a, false), c)
   n = numerator(a, false)*divexact(c, g)
   d = divexact(denominator(a, false), g)
   return parent(a)(n, d)
end

function *(a::FracElem, b::Rational)
   bnum = base_ring(a)(numerator(b))
   bden = base_ring(a)(denominator(b))
   g1 = gcd(denominator(a, false), bnum)
   g2 = gcd(numerator(a, false), bden)
   n = divexact(numerator(a, false), g2)*divexact(bnum, g1)
   d = divexact(denominator(a, false), g1)*divexact(bden, g2)
   return parent(a)(n, d)
end

function *(a::Union{Integer, AbstractFloat}, b::FracElem)
   c = base_ring(b)(a)
   g = gcd(denominator(b, false), c)
   n = numerator(b, false)*divexact(c, g)
   d = divexact(denominator(b, false), g)
   return parent(b)(n, d)
end

function *(a::Rational, b::FracElem)
   anum = base_ring(b)(numerator(a))
   aden = base_ring(b)(denominator(a))
   g1 = gcd(denominator(b, false), anum)
   g2 = gcd(numerator(b, false), aden)
   n = divexact(numerator(b, false), g2)*divexact(anum, g1)
   d = divexact(denominator(b, false), g1)*divexact(aden, g2)
   return parent(b)(n, d)
end

function *(a::FracElem{T}, b::T) where {T <: RingElem}
   g = gcd(denominator(a, false), b)
   n = numerator(a, false)*divexact(b, g)
   d = divexact(denominator(a, false), g)
   return parent(a)(n, d)
end

function *(a::T, b::FracElem{T}) where {T <: RingElem}
   g = gcd(denominator(b, false), a)
   n = numerator(b, false)*divexact(a, g)
   d = divexact(denominator(b, false), g)
   return parent(b)(n, d)
end

function +(a::FracElem, b::Union{Integer, AbstractFloat})
   n = numerator(a, false) + denominator(a, false)*b
   d = denominator(a, false)
   return parent(a)(n, deepcopy(d))
end

+(a::FracElem, b::Rational) = return a + parent(a)(b)

function -(a::FracElem, b::Union{Integer, AbstractFloat})
   n = numerator(a, false) - denominator(a, false)*b
   d = denominator(a, false)
   return parent(a)(n, deepcopy(d))
end

-(a::FracElem, b::Rational) = return a - parent(a)(b)

+(a::Union{Integer, Rational, AbstractFloat}, b::FracElem) = b + a

function -(a::Union{Integer, AbstractFloat}, b::FracElem)
   n = a*denominator(b, false) - numerator(b, false)
   d = denominator(b, false)
   return parent(b)(n, deepcopy(d))
end

-(a::Rational, b::FracElem) = return parent(b)(a) - b

function +(a::FracElem{T}, b::T) where {T <: RingElem}
   n = numerator(a, false) + denominator(a, false)*b
   d = denominator(a, false)
   return parent(a)(n, deepcopy(d))
end

function -(a::FracElem{T}, b::T) where {T <: RingElem}
   n = numerator(a, false) - denominator(a, false)*b
   d = denominator(a, false)
   return parent(a)(n, deepcopy(d))
end

+(a::T, b::FracElem{T}) where {T <: RingElem} = b + a

function -(a::T, b::FracElem{T}) where {T <: RingElem}
   n = a*denominator(b, false) - numerator(b, false)
   d = denominator(b, false)
   return parent(b)(n, deepcopy(d))
end

###############################################################################
#
#   Comparisons
#
###############################################################################

@doc raw"""
    ==(x::FracElem{T}, y::FracElem{T}) where {T <: RingElem}

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(x::FracElem{T}, y::FracElem{T}) where {T <: RingElem}
   b  = check_parent(x, y, false)
   !b && return false

   return (denominator(x, false) == denominator(y, false) &&
           numerator(x, false) == numerator(y, false)) ||
          (denominator(x, true) == denominator(y, true) &&
           numerator(x, true) == numerator(y, true)) ||
          (numerator(x, false)*denominator(y, false) ==
           denominator(x, false)*numerator(y, false))
end

@doc raw"""
    isequal(x::FracElem{T}, y::FracElem{T}) where {T <: RingElem}

Return `true` if $x == y$ exactly, otherwise return `false`. This function is
useful in cases where the numerators and denominators of the fractions are
inexact, e.g. power series. Only if the power series are precisely the same,
to the same precision, are they declared equal by this function.
"""
function isequal(x::FracElem{T}, y::FracElem{T}) where {T <: RingElem}
   if parent(x) != parent(y)
      return false
   end
   return isequal(numerator(x, false)*denominator(y, false),
                  denominator(x, false)*numerator(y, false))
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

@doc raw"""
    ==(x::FracElem, y::Union{Integer, Rational, AbstractFloat})

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
function ==(x::FracElem, y::Union{Integer, AbstractFloat})
   return (isone(denominator(x, false)) && numerator(x, false) == y) ||
          (isone(denominator(x, true)) && numerator(x, true) == y) ||
          (numerator(x, false) == denominator(x, false)*y)
end

function ==(x::FracElem, y::Rational)
   return (numerator(x, false) == numerator(y, false) &&
           denominator(x, false) == denominator(y, false)) ||
           (numerator(x, false)*denominator(y, false) ==
            denominator(x, false)*numerator(y, false))
end

@doc raw"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::FracElem)

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::FracElem) = y == x

@doc raw"""
    ==(x::FracElem{T}, y::T) where {T <: RingElem}

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
function ==(x::FracElem{T}, y::T) where {T <: RingElem}
   return (isone(denominator(x, false)) && numerator(x, false) == y) ||
          (isone(denominator(x, true)) && numerator(x, true) == y) ||
          (numerator(x, false) == denominator(x, false)*y)
end

@doc raw"""
    ==(x::T, y::FracElem{T}) where {T <: RingElem}

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::T, y::FracElem{T}) where {T <: RingElem} = y == x

###############################################################################
#
#   Inversion
#
###############################################################################

@doc raw"""
    Base.inv(a::FracElem)

Return the inverse of the fraction $a$.
"""
function Base.inv(a::FracElem)
   iszero(numerator(a, false)) && throw(NotInvertibleError(a))
   return parent(a)(deepcopy(denominator(a, false)),
                    deepcopy(numerator(a, false)))
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::FracElem{T}, b::FracElem{T}; check::Bool=true) where {T <: RingElem}
   check_parent(a, b)
   n1 = numerator(a, false)
   d2 = denominator(b, false)
   n2 = numerator(b, false)
   d1 = denominator(a, false)
   if d1 == n2
      n = n1*d2
      d = d1*n2
   elseif isone(d1)
      gd = gcd(n1, n2)
      if isone(gd)
         n = n1*d2
         d = deepcopy(n2)
      else
         n = divexact(n1, gd; check=check)*d2
         d = divexact(n2, gd; check=check)
      end
   elseif isone(n2)
      gd = gcd(d2, d1)
      if isone(gd)
         n = d2*n1
         d = deepcopy(d1)
      else
         n = divexact(d2, gd; check=check)*n1
         d = divexact(d1, gd; check=check)
      end
   else
      g1 = gcd(n1, n2)
      g2 = gcd(d2, d1)
      if !isone(g1)
         n1 = divexact(n1, g1; check=check)
         n2 = divexact(n2, g1; check=check)
      end
      if !isone(g2)
         d2 = divexact(d2, g2; check=check)
         d1 = divexact(d1, g2; check=check)
      end
      n = n1*d2
      d = d1*n2
   end
   return parent(a)(n, d)
end

function divides(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   if iszero(a)
     return true, parent(a)()
   end
   if iszero(b)
     return false, parent(a)()
   end
   return true, divexact(a, b)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(a::FracElem, b::Union{Integer, AbstractFloat}; check::Bool=true)
   b == 0 && throw(DivideError())
   c = base_ring(a)(b)
   g = gcd(numerator(a, false), c)
   n = divexact(numerator(a, false), g; check=false)
   d = denominator(a, false)*divexact(c, g; check=false)
   return parent(a)(n, d)
end

function divexact(a::Union{Integer, AbstractFloat}, b::FracElem{T}; check::Bool=true) where T <: RingElem
   iszero(b) && throw(DivideError())
   c = base_ring(b)(a)
   g = gcd(numerator(b, false), c)
   n = denominator(b, false)*divexact(c, g; check=false)
   d = divexact(numerator(b, false), g; check=false)
   return parent(b)(n, d)
end

function divexact(a::FracElem{T}, b::Rational; check::Bool=true) where T <: RingElem
   return divexact(a, parent(a)(b), check=check)
end

function divexact(a::Rational, b::FracElem{T}; check::Bool=true) where T <: RingElem
   return divexact(parent(b)(a), b, check=check)
end

function divexact(a::FracElem{T}, b::T; check::Bool=true) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   g = gcd(numerator(a, false), b)
   n = divexact(numerator(a, false), g; check=false)
   d = denominator(a, false)*divexact(b, g; check=false)
   return parent(a)(n, d)
end

function divexact(a::T, b::FracElem{T}; check::Bool=true) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   g = gcd(numerator(b, false), a)
   n = denominator(b, false)*divexact(a, g; check=false)
   d = divexact(numerator(b, false), g; check=false)
   return parent(b)(n, d)
end

##############################################################################
#
#  Evaluation
#
##############################################################################

function evaluate(f::FracElem{T}, V::Vector{U}) where {T <: RingElement, U <: RingElement}
    return evaluate(numerator(f), V)//evaluate(denominator(f), V)
end
  
function evaluate(f::FracElem{T}, v::U) where {T <: RingElement, U <: RingElement}
    return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

function evaluate(f::FracElem{T}, v::U) where {T <: PolyRingElem, U <: Integer}
    return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

function evaluate(f::FracElem{T}, vars::Vector{Int}, vals::Vector{U}) where {T <: RingElement, U <: RingElement}
     return evaluate(numerator(f), vars, vals)//evaluate(denominator(f), vars, vals)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::FracElem{T}, b::Int) where {T <: RingElem}
   if b < 0
      a = inv(a)
      b = -b
   end
   return parent(a)(numerator(a)^b, denominator(a)^b)
end

##############################################################################
#
#  Derivative
#
##############################################################################

# Return the derivative with respect to `x`.
function derivative(f::FracElem{T}, x::T) where {T <: MPolyRingElem}
    return derivative(f, var_index(x))
end
  
# Return the derivative with respect to the `i`-th variable.
function derivative(f::FracElem{T}, i::Int) where {T <: MPolyRingElem}
    n = numerator(f)
    d = denominator(f)
    return (derivative(n, i)*d - n*derivative(d, i))//d^2
end

function derivative(f::FracElem{T}) where {T <: PolyRingElem}
    n = numerator(f)
    d = denominator(f)
    return (derivative(n)*d - n*derivative(d))//d^2
end

###############################################################################
#
#   Square root
#
###############################################################################

@doc raw"""
    is_square(a::FracElem{T}) where T <: RingElem

Return `true` if $a$ is a square.
"""
function is_square(a::FracElem{T}) where T <: RingElem
   return is_square(numerator(a)) && is_square(denominator(a))
end

@doc raw"""
    Base.sqrt(a::FracElem{T}; check::Bool=true) where T <: RingElem

Return the square root of $a$. By default the function will throw an
exception if the input is not square. If `check=false` this test is omitted.
"""
function Base.sqrt(a::FracElem{T}; check::Bool=true) where T <: RingElem
   return parent(a)(sqrt(numerator(a); check=check), sqrt(denominator(a); check=check))
end

function is_square_with_sqrt(a::FracElem{T}) where T <: RingElem
   S = parent(a)
   f1, s1 = is_square_with_sqrt(numerator(a))
   if !f1
      return false, zero(S)
   end
   f2, s2 = is_square_with_sqrt(denominator(a))
   if !f2
      return false, zero(S)
   end
   return true, s1//s2
end

###############################################################################
#
#   GCD
#
###############################################################################

@doc raw"""
    gcd(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}

Return a greatest common divisor of $a$ and $b$ if one exists. N.B: we define
the GCD of $a/b$ and $c/d$ to be gcd$(ad, bc)/bd$, reduced to lowest terms.
This requires the existence of a greatest common divisor function for the
base ring.
"""
function gcd(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   gbd = gcd(denominator(a, false), denominator(b, false))
   n = gcd(numerator(a, false), numerator(b, false))
   d = divexact(denominator(a, false), gbd)*denominator(b, false)
   u = canonical_unit(n)
   if !iszero(u)
      n = divexact(n, u)
   end
   d = divexact(d, canonical_unit(d))
   return parent(a)(n, d)
end

################################################################################
#
#   Remove and valuation
#
################################################################################

@doc raw"""
    remove(z::FracElem{T}, p::T) where {T <: RingElem}

Return the tuple $n, x$ such that $z = p^nx$ where $x$ has valuation $0$ at
$p$.
"""
function remove(z::FracElem{T}, p) where {T}
   p = convert(T, p)
   iszero(z) && error("Not yet implemented")
   v, d = remove(denominator(z, false), p)
   w, n = remove(numerator(z, false), p)
   return w-v, parent(z)(deepcopy(n), deepcopy(d))
end

@doc raw"""
    valuation(z::FracElem{T}, p::T) where {T <: RingElem}

Return the valuation of $z$ at $p$.
"""
function valuation(z::FracElem{T}, p) where {T}
   p = convert(T, p)
   v, _ = remove(z, p)
   return v
end

###############################################################################
#
#   Random functions
#
###############################################################################

RandomExtensions.maketype(R::FracField, _) = elem_type(R)

function RandomExtensions.make(S::FracField, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      RandomExtensions.Make(S, vs[1]) # forward to default Make constructor
   else
      Make(S, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:RingElement, <:FracField}})
   S, v = sp[][1:end]
   R = base_ring(S)
   n = rand(rng, v)
   d = R()
   while iszero(d)
      d = rand(rng, v)
   end
   return S(n, d)
end

rand(rng::AbstractRNG, S::FracField, v...) =
   rand(rng, make(S, v...))

rand(S::FracField, v...) = rand(GLOBAL_RNG, S, v...)

###############################################################################
#
#   fraction_field constructor
#
###############################################################################

@doc raw"""
    fraction_field(R::Ring; cached::Bool=true)

Return the parent object of the fraction field over the given base ring $R$.
If `cached == true` (the default), the returned parent object is cached so
that it will always be returned by a call to the constructor when the same
base ring $R$ is supplied.
"""
function fraction_field(R::Ring; cached::Bool=true)
   @req !is_trivial(R) "Zero rings are currently not supported as coefficient ring."
   return Generic.fraction_field(R; cached=cached)
end

@doc raw"""
    FactoredFractionField(R::Ring; cached::Bool=true)

Return the parent object of the fraction field over the given base ring $R$,
where the elements are maintained in factored form as much as possible.
"""
function FactoredFractionField(R::Ring; cached::Bool=true)
   return Generic.FactoredFractionField(R; cached=cached)
end

