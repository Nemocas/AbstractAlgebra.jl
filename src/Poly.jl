###############################################################################
#
#   Poly.jl : Univariate polynomials
#
###############################################################################

export PolyCoeffs, PolynomialRing, PolyRing, addmul!, characteristic,
       chebyshev_t, chebyshev_u, coefficients, compose, constant_coefficient,
       content, deflate, deflation, degree, derivative, discriminant, divexact,
       divexact_low, divhigh, divides, evaluate, gcdinv, inflate, integral,
       interpolate, ismonic, issquare, isterm, isterm_recursive,
       map_coefficients, modulus,  monomial_to_newton!, mul_classical,
       mulhigh_n, mul_karatsuba, mul_ks, mullow, mulmod, newton_to_monomial!,
       nvars, polynomial, pow_multinomial, primpart, pseudodivrem, pseudorem,
       remove, resultant, resultant_ducos, resultant_euclidean,
       resultant_lehmer, resultant_subresultant, resultant_sylvester, resx,
       shift_left, shift_right, subst, sylvester_matrix, symbols, tail,
       valuation, var

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring(R::PolyRing{T}) where T <: RingElement = R.base_ring::parent_type(T)

base_ring(a::PolynomialElem) = base_ring(parent(a))

parent(a::PolynomialElem) = a.parent

function isdomain_type(::Type{T}) where {S <: RingElement, T <: PolyElem{S}}
   return isdomain_type(S)
end

function isexact_type(a::Type{T}) where {S <: RingElement, T <: PolyElem{S}}
   return isexact_type(S)
end

@doc Markdown.doc"""
    var(a::PolyRing)

Return the internal name of the generator of the polynomial ring. Note that
this is returned as a `Symbol` not a `String`.
"""
var(a::PolyRing) = a.S

@doc Markdown.doc"""
    symbols(a::PolyRing)

Return an array of the variable names for the polynomial ring. Note that
this is returned as an array of `Symbol` not `String`.
"""
symbols(a::PolyRing) = [a.S]

@doc Markdown.doc"""
    nvars(a::PolyRing)

Return the number of variables of the polynomial ring, which is 1.
"""
nvars(a::PolyRing) = 1

function check_parent(a::PolynomialElem, b::PolynomialElem, throw::Bool = true)
   c = parent(a) != parent(b)
   c && throw && error("Incompatible polynomial rings in polynomial operation")
   return !c
end

characteristic(a::PolyRing) = characteristic(base_ring(a))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::PolyElem, h::UInt)
   b = 0x53dd43cd511044d1%UInt
   for i in 0:length(a) - 1
      b = xor(b, xor(hash(coeff(a, i), h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

length(a::PolynomialElem) = a.length

@doc Markdown.doc"""
    degree(a::PolynomialElem)

Return the degree of the given polynomial. This is defined to be one less
than the length, even for constant polynomials.
"""
degree(a::PolynomialElem) = length(a) - 1

@doc Markdown.doc"""
    modulus(a::PolyElem{T}) where {T <: ResElem}

Return the modulus of the coefficients of the given polynomial.
"""
modulus(a::PolyElem{T}) where {T <: ResElem} = modulus(base_ring(a))

@doc Markdown.doc"""
    leading_coefficient(a::PolynomialElem)

Return the leading coefficient of the given polynomial. This will be the
nonzero coefficient of the term with highest degree unless the polynomial
in the zero polynomial, in which case a zero coefficient is returned.
"""
function leading_coefficient(a::PolynomialElem)
   return length(a) == 0 ? base_ring(a)(0) : coeff(a, length(a) - 1)
end

@doc Markdown.doc"""
    trailing_coefficient(a::PolynomialElem)

Return the trailing coefficient of the given polynomial. This will be the
nonzero coefficient of the term with lowest degree unless the polynomial
is the zero polynomial, in which case a zero coefficient is returned.
"""
function trailing_coefficient(a::PolynomialElem)
   if iszero(a)
      return base_ring(a)(0)
   else
      for i = 1:length(a)
         c = coeff(a, i - 1)
         if !iszero(c)
            return c
         end
      end
      return coeff(a, length(a) - 1)
   end
end

@doc Markdown.doc"""
    constant_coefficient(a::PolynomialElem)

Return the constant coefficient of the given polynomial. If the polynomial is
the zero polynomial, the function will return zero.
"""
function constant_coefficient(a::PolynomialElem)
   if iszero(a)
      return zero(base_ring(a))
   end
   return coeff(a, 0)
end

@doc Markdown.doc"""
    tail(a::PolynomialElem)

Return the tail of the given polynomial, i.e. the polynomial without its
leading term (if any).
"""
function tail(a::PolynomialElem)
   return iszero(a) ? zero(parent(a)) : truncate(a, length(a) - 1)
end

@doc Markdown.doc"""
    zero(R::PolyRing)

Return the zero polynomial in the given polynomial ring.
"""
zero(R::PolyRing) = R(0)

one(R::PolyRing) = R(1)

@doc Markdown.doc"""
    gen(R::PolyRing)

Return the generator of the given polynomial ring.
"""
gen(R::PolyRing) = R([zero(base_ring(R)), one(base_ring(R))])

@doc Markdown.doc"""
    gens(R::PolyRing)

Return an array containing the generator of the given polynomial ring.
"""
gens(R::PolyRing) = [gen(R)]

iszero(a::PolynomialElem) = length(a) == 0

isone(a::PolynomialElem) = length(a) == 1 && isone(coeff(a, 0))

@doc Markdown.doc"""
    isgen(a::PolynomialElem)

Return `true` if the given polynomial is the constant generator of its
polynomial ring, otherwise return `false`.
"""
function isgen(a::PolynomialElem)
    return length(a) == 2 && iszero(coeff(a, 0)) && isone(coeff(a, 1))
end

@doc Markdown.doc"""
    ismonic(a::PolynomialElem)

Return `true` if the given polynomial is monic, i.e. has leading coefficient
equal to one, otherwise return `false`.
"""
function ismonic(a::PolynomialElem)
    return isone(leading_coefficient(a))
end

isunit(a::PolynomialElem) = length(a) == 1 && isunit(coeff(a, 0))

###############################################################################
#
#  Monomial and term
#
###############################################################################

@doc Markdown.doc"""
    isterm(a::PolynomialElem)

Return `true` if the given polynomial has one term.
"""
function isterm(a::PolynomialElem)
   if iszero(a)
      return false
   end
   for i = 1:length(a) - 1
      if !iszero(coeff(a, i - 1))
         return false
      end
   end
   return true
end

isterm_recursive(a::T) where T <: RingElement = true

@doc Markdown.doc"""
    isterm_recursive(a::PolynomialElem)

Return `true` if the given polynomial has one term. This function is
recursive, with all scalar types returning true.
"""
function isterm_recursive(a::PolynomialElem)
   if !isterm_recursive(leading_coefficient(a))
      return false
   end
   for i = 1:length(a) - 1
      if !iszero(coeff(a, i - 1))
         return false
      end
   end
   return true
end

@doc Markdown.doc"""
    ismonomial_recursive(a::PolynomialElem)

Return `true` if the given polynomial is a monomial.
"""
function ismonomial(a::PolynomialElem)
   if !isone(leading_coefficient(a))
      return false
   end
   for i = 1:length(a) - 1
      if !iszero(coeff(a, i - 1))
         return false
      end
   end
   return true
end

ismonomial_recursive(a::T) where T <: RingElement = isone(a)

@doc Markdown.doc"""
    ismonomial_recursive(a::PolynomialElem)

Return `true` if the given polynomial is a monomial. This function is
recursive, with all scalar types returning true.
"""
function ismonomial_recursive(a::PolynomialElem)
   if !ismonomial_recursive(leading_coefficient(a))
      return false
   end
   for i = 1:length(a) - 1
      if !iszero(coeff(a, i - 1))
         return false
      end
   end
   return true
end

###############################################################################
#
#   Similar and zero
#
###############################################################################

function similar(x::PolyElem, R::Ring, var::Symbol=var(parent(x)); cached::Bool=true)
   TT = elem_type(R)
   V = Vector{TT}(undef, 0)
   p = Generic.Poly{TT}(V)
   # Default similar is supposed to return a polynomial
   p.parent = Generic.PolyRing{TT}(R, var, cached)
   p = set_length!(p, 0)
   return p
end

function similar(x::PolyElem, var::Symbol=var(parent(x)); cached::Bool=true)
   return similar(x, base_ring(x), var; cached=cached)
end

function similar(x::PolyElem, R::Ring, var::String; cached::Bool=true)
   return similar(x, R, Symbol(var); cached=cached)
end

function similar(x::PolyElem, var::String; cached::Bool=true)
   return similar(x, base_ring(x), Symbol(var); cached=cached)
end

zero(p::PolyElem, R::Ring, var::Symbol=var(parent(p)); cached::Bool=true) =
   similar(p, R, var; cached=cached)

zero(p::PolyElem, var::Symbol=var(parent(p)); cached::Bool=true) =
   similar(p, base_ring(p), var; cached=cached)

zero(p::PolyElem, R::Ring, var::String; cached::Bool=true) =
   zero(p, R, Symbol(var); cached=cached)

zero(p::PolyElem, var::String; cached::Bool=true) =
   zero(p, base_ring(p), Symbol(var); cached=cached)

###############################################################################
#
#   polynomial constructor
#
###############################################################################

function polynomial(R::Ring, arr::Vector{T}, var::AbstractString="x"; cached::Bool=true) where T
   TT = elem_type(R)
   coeffs = T == Any && length(arr) == 0 ? elem_type(R)[] : map(R, arr)
   p = Generic.Poly{TT}(coeffs)
   # Default is supposed to return a polynomial
   p.parent = Generic.PolyRing{TT}(R, Symbol(var), cached)
   return p
end

###############################################################################
#
#  Iterators
#
###############################################################################

struct PolyCoeffs{T <: RingElement}
   f::T
end
 
function coefficients(f::PolyElem)
   return PolyCoeffs(f)
end
 
function Base.iterate(PC::PolyCoeffs{<:PolyElem}, st::Int = -1)
   st += 1
   if st > degree(PC.f)
       return nothing
   else
       return coeff(PC.f, st), st
   end
end
 
Base.IteratorEltype(M::PolyElem) = Base.HasEltype()

Base.eltype(M::PolyElem{T}) where {T} = T

Base.eltype(M::PolyCoeffs) = Base.eltype(M.f)
 
Base.IteratorSize(M::PolyCoeffs{<:PolyElem}) = Base.HasLength()

Base.length(M::PolyCoeffs{<:PolyElem}) = length(M.f)
 
function Base.lastindex(a::PolyCoeffs{<:PolyElem})
   return degree(a.f)
end
 
function Base.getindex(a::PolyCoeffs{<:PolyElem}, i::Int)
   return coeff(a.f, i)
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::PolynomialElem) = canonical_unit(leading_coefficient(x))

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(@nospecialize(a::Union{PolynomialElem, NCPolyElem}),
   x = var(parent(a)); context = nothing)
   sum = Expr(:call, :+)
   for k in degree(a):-1:0
      c = coeff(a, k)
      if !iszero(c)
         xk = k < 1 ? 1 : k == 1 ? x : Expr(:call, :^, x, k)
         if isone(c)
            push!(sum.args, Expr(:call, :*, xk))
         else
            push!(sum.args, Expr(:call, :*, expressify(c, context = context), xk))
         end
      end
   end
   return sum
end

function Base.show(io::IO, ::MIME"text/plain", a::Union{PolynomialElem, NCPolyElem})
   print(io, obj_to_string(a, context = io))
end

function Base.show(io::IO, a::Union{PolynomialElem, NCPolyElem})
   print(io, obj_to_string(a, context = io))
end

function show(io::IO, p::PolyRing)
   print(io, "Univariate Polynomial Ring in ")
   print(io, string(var(p)))
   print(io, " over ")
   print(IOContext(io, :compact => true), base_ring(p))
end

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(a::PolynomialElem)
   len = length(a)
   z = parent(a)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, -coeff(a, i - 1))
   end
   z = set_length!(z, len)
   return z
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   lenz = max(lena, lenb)
   z = parent(a)()
   fit!(z, lenz)
   i = 1
   while i <= min(lena, lenb)
      z = setcoeff!(z, i - 1, coeff(a, i - 1) + coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      z = setcoeff!(z, i - 1, deepcopy(coeff(a, i - 1)))
      i += 1
   end
   while i <= lenb
      z = setcoeff!(z, i - 1, deepcopy(coeff(b, i - 1)))
      i += 1
   end
   z = set_length!(z, normalise(z, i - 1))
   return z
end

function -(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   lenz = max(lena, lenb)
   z = parent(a)()
   fit!(z, lenz)
   i = 1
   while i <= min(lena, lenb)
      z = setcoeff!(z, i - 1, coeff(a, i - 1) - coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      z = setcoeff!(z, i - 1, deepcopy(coeff(a, i - 1)))
      i += 1
   end
   while i <= lenb
      z = setcoeff!(z, i - 1, -coeff(b, i - 1))
      i += 1
   end
   z = set_length!(z, normalise(z, i - 1))
   return z
end

@doc Markdown.doc"""
    mul_karatsuba(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement

Return $a \times b$ using one non-recursive application the Karatsuba algorithm.
"""
function mul_karatsuba(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   # we assume len(a) != 0 != lenb and parent(a) == parent(b)
   lena = length(a)
   lenb = length(b)
   m = div(max(lena, lenb) + 1, 2)
   if m < lena
      a1 = shift_right(a, m)
      a0 = truncate(a, m)
   else
      return a*truncate(b, m) + shift_left(a*shift_right(b, m), m)
   end
   if a !== b
      if m < lenb
         b1 = shift_right(b, m)
         b0 = truncate(b, m)
      else
         return b*truncate(a, m) + shift_left(b*shift_right(a, m), m)
      end
   else
      b1 = a1
      b0 = a0
   end
   z0 = a0*b0
   z2 = a1*b1
   if a !== b
      z1 = (a1 + a0)*(b1 + b0) - z2 - z0
   else
      s = a1 + a0
      z1 = s*s - z2 - z0
   end
   r = parent(a)()
   fit!(r, lena + lenb - 1)
   for i = 1:length(z0)
      r = setcoeff!(r, i - 1, coeff(z0, i - 1))
   end
   for i = length(z0) + 1:2m
      r = setcoeff!(r, i - 1, base_ring(a)())
   end
   for i = 1:length(z2)
      r = setcoeff!(r, 2m + i - 1, coeff(z2, i - 1))
   end
   for i = 1:length(z1)
      u = coeff(r, i + m - 1)
      u = addeq!(u, coeff(z1, i - 1))
      setcoeff!(r, i + m - 1, u)
   end
   return r
end

function mul_ks(a::PolyElem{T}, b::PolyElem{T}) where {T <: PolyElem}
   lena = length(a)
   lenb = length(b)
   if lena == 0 || lenb == 0
      return parent(a)()
   end
   maxa = 0
   nza = 0
   for i = 1:lena
      lenc = length(coeff(a, i - 1))
      maxa = max(lenc, maxa)
      nza += (lenc == 0 ? 0 : 1)
   end
   if a !== b
      maxb = 0
      nzb = 0
      for i = 1:lenb
         lenc = length(coeff(b, i - 1))
         maxb = max(lenc, maxb)
         nzb += (lenc == 0 ? 0 : 1)
      end
   else
      maxb = maxa
      nzb = nza
   end
   if nza*nzb < 4*max(lena, lenb)
      return mul_classical(a, b)
   end
   m = maxa + maxb - 1
   z = base_ring(base_ring(a))()
   A1 = Array{elem_type(base_ring(base_ring(a)))}(undef, m*lena)
   for i = 1:lena
      c = coeff(a, i - 1)
      for j = 1:length(c)
         A1[(i - 1)*m + j] = coeff(c, j - 1)
      end
      for j = length(c) + 1:m
         A1[(i - 1)*m + j] = z
      end
   end
   ksa = base_ring(a)(A1)
   if a !== b
      A2 = Array{elem_type(base_ring(base_ring(a)))}(undef, m*lenb)
      for i = 1:lenb
         c = coeff(b, i - 1)
         for j = 1:length(c)
            A2[(i - 1)*m + j] = coeff(c, j - 1)
         end
         for j = length(c) + 1:m
            A2[(i - 1)*m + j] = z
         end
      end
      ksb = base_ring(b)(A2)
   else
      ksb = ksa
   end
   p = ksa*ksb
   r = parent(a)()
   lenr = lena + lenb - 1
   fit!(r, lenr)
   for i = 1:lenr
      u = coeff(r, i - 1)
      fit!(u, m)
      for j = 1:m
         u = setcoeff!(u, j - 1, coeff(p, (i - 1)*m + j - 1))
      end
      setcoeff!(r, i - 1, u)
   end
   r = set_length!(r, normalise(r, lenr))
   return r
end

function mul_classical(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   lena = length(a)
   lenb = length(b)
   if lena == 0 || lenb == 0
      return parent(a)()
   end
   R = base_ring(a)
   t = R()
   lenz = lena + lenb - 1
   d = Array{T}(undef, lenz)
   for i = 1:lena
      d[i] = mul_red!(R(), coeff(a, i - 1), coeff(b, 0), false)
   end
   for i = 2:lenb
      d[lena + i - 1] = mul_red!(R(), coeff(a, lena - 1), coeff(b, i - 1), false)
   end
   for i = 1:lena - 1
      for j = 2:lenb
         t = mul_red!(t, coeff(a, i - 1), coeff(b, j - 1), false)
         d[i + j - 1] = addeq!(d[i + j - 1], t)
      end
   end
   for i = 1:lenz
      d[i] = reduce!(d[i])
   end
   z = parent(a)(d)
   z = set_length!(z, normalise(z, lenz))
   return z
end

function *(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   check_parent(a, b)
   return mul_classical(a, b)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::T, b::PolyElem{T}) where {T <: RingElem}
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   return z
end

function *(a::Union{Integer, Rational, AbstractFloat}, b::PolynomialElem)
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   return z
end

*(a::PolyElem{T}, b::T) where {T <: RingElem} = b*a

*(a::PolynomialElem, b::Union{Integer, Rational, AbstractFloat}) = b*a

###############################################################################
#
#   Powering
#
###############################################################################

function pow_multinomial(a::PolyElem{T}, e::Int) where T <: RingElement
   e < 0 && throw(DomainError(e, "exponent must be >= 0"))
   lena = length(a)
   lenz = (lena - 1) * e + 1
   res = Array{T}(undef, lenz)
   for k = 1:lenz
      res[k] = base_ring(a)()
   end
   d = base_ring(a)()
   first = coeff(a, 0)
   res[1] = first ^ e
   for k = 1 : lenz - 1
      u = -k
      for i = 1 : min(k, lena - 1)
         t = coeff(a, i) * res[(k - i) + 1]
         u += e + 1
         res[k + 1] = addeq!(res[k + 1], t * u)
      end
      d = addeq!(d, first)
      res[k + 1] = divexact(res[k + 1], d)
   end
   z = parent(a)(res)
   z = set_length!(z, normalise(z, lenz))
   return z
end

@doc Markdown.doc"""
    ^(a::PolyElem{T}, b::Int) where T <: RingElement

Return $a^b$. We require $b \geq 0$.
"""
function ^(a::PolyElem{T}, b::Int) where T <: RingElement
   b < 0 && throw(DomainError(b, "exponent must be >= 0"))
   # special case powers of x for constructing polynomials efficiently
   R = parent(a)
   if isgen(a)
      z = R()
      fit!(z, b + 1)
      z = setcoeff!(z, b, deepcopy(coeff(a, 1)))
      for i = 1:b
         z = setcoeff!(z, i - 1, deepcopy(coeff(a, 0)))
      end
      z = set_length!(z, b + 1)
      return z
   elseif b == 0
      return one(R)
   elseif length(a) == 0
      return zero(R)
   elseif length(a) == 1
      return R(coeff(a, 0)^b)
   elseif b == 1
      return deepcopy(a)
   else
      if T <: FieldElement && characteristic(base_ring(R)) == 0
         zn = 0
         while iszero(coeff(a, zn))
            zn += 1
         end
         if length(a) - zn < 8 && b > 4
             f = shift_right(a, zn)
             return shift_left(pow_multinomial(f, b), zn*b)
         end
      end
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
#   Comparisons
#
###############################################################################

@doc Markdown.doc"""
    ==(x::PolyElem{T}, y::PolyElem{T}) where T <: RingElement

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(x::PolyElem{T}, y::PolyElem{T}) where T <: RingElement
   b = check_parent(x, y, false)
   !b && return false
   if length(x) != length(y)
      return false
   else
      for i = 1:length(x)
         if coeff(x, i - 1) != coeff(y, i - 1)
            return false
         end
      end
   end
   return true
end

@doc Markdown.doc"""
    isequal(x::PolyElem{T}, y::PolyElem{T}) where T <: RingElement

Return `true` if $x == y$ exactly, otherwise return `false`. This function is
useful in cases where the coefficients of the polynomial are inexact, e.g.
power series. Only if the power series are precisely the same, to the same
precision, are they declared equal by this function.
"""
function isequal(x::PolyElem{T}, y::PolyElem{T}) where T <: RingElement
   if parent(x) != parent(y)
      return false
   end
   if length(x) != length(y)
      return false
   end
   for i = 1:length(x)
      if !isequal(coeff(x, i - 1), coeff(y, i - 1))
         return false
      end
   end
   return true
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

@doc Markdown.doc"""
    ==(x::PolyElem{T}, y::T) where {T <: RingElem}

Return `true` if $x == y$.
"""
==(x::PolyElem{T}, y::T) where T <: RingElem = ((length(x) == 0 && iszero(y))
                        || (length(x) == 1 && coeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::PolynomialElem, y::Union{Integer, Rational, AbstractFloat})

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::PolynomialElem, y::Union{Integer, Rational, AbstractFloat}) = ((length(x) == 0 && iszero(base_ring(x)(y)))
                        || (length(x) == 1 && coeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::T, y::PolyElem{T}) where T <: RingElem = y == x

Return `true` if $x = y$.
"""
==(x::T, y::PolyElem{T}) where T <: RingElem = y == x

@doc Markdown.doc"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::PolyElem)

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::PolyElem) = y == x

###############################################################################
#
#   Approximation
#
###############################################################################

function Base.isapprox(f::PolynomialElem, g::PolynomialElem; atol::Real=sqrt(eps()))
   check_parent(f, g)
   nmin = min(length(f), length(g))
   i = 1
   while i <= nmin
      if !isapprox(coeff(f, i - 1), coeff(g, i - 1); atol=atol)
         return false
      end
      i += 1
   end
   while i <= length(f)
      if !isapprox(coeff(f, i - 1), 0; atol=atol)
         return false
      end
      i += 1
   end
   while i <= length(g)
      if !isapprox(coeff(g, i - 1), 0; atol=atol)
         return false
      end
      i += 1
   end
   return true
end

function Base.isapprox(f::PolynomialElem{T}, g::T; atol::Real=sqrt(eps())) where T
   return isapprox(f, parent(f)(g); atol=atol)
end

function Base.isapprox(f::T, g::PolynomialElem{T}; atol::Real=sqrt(eps())) where T
   return isapprox(parent(g)(f), g; atol=atol)
end

###############################################################################
#
#   Truncation
#
###############################################################################

@doc Markdown.doc"""
    truncate(a::PolynomialElem, n::Int)

Return $a$ truncated to $n$ terms.
"""
function truncate(a::PolynomialElem, n::Int)
   lena = length(a)
   if lena <= n
      return a
   end
   lenz = min(lena, n)
   z = parent(a)()
   fit!(z, lenz)
   for i = 1:lenz
      z = setcoeff!(z, i - 1, coeff(a, i - 1))
   end
   z = set_length!(z, normalise(z, lenz))
   return z
end

@doc Markdown.doc"""
    mullow(a::PolyElem{T}, b::PolyElem{T}, n::Int) where T <: RingElement

Return $a\times b$ truncated to $n$ terms.
"""
function mullow(a::PolyElem{T}, b::PolyElem{T}, n::Int) where T <: RingElement
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   if lena == 0 || lenb == 0
      return zero(parent(a))
   end
   if n < 0
      n = 0
   end
   R = base_ring(a)
   t = R()
   lenz = min(lena + lenb - 1, n)
   d = Array{T}(undef, lenz)
   for i = 1:min(lena, lenz)
      d[i] = mul_red!(R(), coeff(a, i - 1), coeff(b, 0), false)
   end
   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = mul_red!(R(), coeff(a, lena - 1), coeff(b, j - 1), false)
      end
   end
   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            t = mul_red!(t, coeff(a, i - 1), coeff(b, j - 1), false)
            d[i + j - 1] = addeq!(d[i + j - 1], t)
         end
      end
   end
   for i = 1:lenz
      d[i] = reduce!(d[i])
   end
   z = parent(a)(d)
   z = set_length!(z, normalise(z, lenz))
   return z
end

# computes the terms of the product from degree deg(a) + deg(b) down to
# deg(a) + deg(b) - n inclusive, setting the remaining terms to zero
function mulhigh_n(a::PolyElem{T}, b::PolyElem{T}, n::Int) where T <: RingElement
    # if a = sum a_i t^i and b = sum b_j t^j
    # want (i, j) such that i + j >= deg a + deg b - n
    r = parent(a)()
    for i = max(degree(a) - n, 0):degree(a)
        for j = max(degree(a) + degree(b) - n - i, 0):degree(b)
            r = setcoeff!(r, i + j, coeff(r, i + j) + coeff(a, i)*coeff(b, j))
        end
    end
    return r
end

# assuming b divides a (behaviour is undefined otherwise), computes the last n
# terms of the quotient, i.e. computes divexact(a, b) mod x^n
function divexact_low(a::PolyElem{T}, b::PolyElem{T}, n::Int) where T <: RingElement
    r = parent(a)()
    if iszero(b)
       return r
    end
    shift = 0
    for i = 0:degree(b)
       if !iszero(coeff(b, i))
          break
       end
       shift += 1
    end
    if shift != 0
        a = shift_right(a, shift)
        b = shift_right(b, shift)
    end
    a = truncate(a, n)
    b = truncate(b, n)
    for i = 0:n - 1
        flag, q = divides(coeff(a, 0), coeff(b, 0))
        !flag && error("Not an exact division")
        r = setcoeff!(r, i, q)
        a = shift_right(a - q*b, 1)
        b = truncate(b, n - i - 1)
        # truncate both a and b to n - i - 1
    end
    return r
end

# computes the top terms of the quotient of a by b starting with the term of
# degree n0
# if the division is not exact until at least the term of degree n0, an
# exception may be raised
function divhigh(a::PolyElem{T}, b::PolyElem{T}, n0::Int) where T <: RingElement
    r = parent(a)()
    n = degree(a) - degree(b) - n0
    fit!(r, degree(a) - degree(b) + 1)
    a = deepcopy(a)
    da = degree(a)
    R = base_ring(a)
    t = R()
    for i = 0:n
        if da < degree(b)
            break
        end
        # negate quotient so we can use addeq! below
        q = -divexact(coeff(a, da), leading_coefficient(b))
        r = setcoeff!(r, da - degree(b), q)
        da -= 1
        if i != n
            c = coeff(a, da)
            for j = 0:min(degree(b) - 1, i)
                t = mul!(t, coeff(r, da - degree(b) + j + 1),
			                          coeff(b, degree(b) - j - 1))
                c = addeq!(c, t)
            end
            a = setcoeff!(a, da, c)
        end
    end
    if iszero(r)
        return r
    end
    r = set_length!(r, normalise(r, length(r)))
    # negate r to compensate for negation above
    return -r
end

###############################################################################
#
#   Reversal
#
###############################################################################

@doc Markdown.doc"""
    reverse(x::PolynomialElem, len::Int)

Return the reverse of the polynomial $x$, thought of as a polynomial of
the given length (the polynomial will be notionally truncated or padded with
zeroes before the leading term if necessary to match the specified length).
The resulting polynomial is normalised. If `len` is negative we throw a
`DomainError()`.
"""
function reverse(x::PolynomialElem, len::Int)
   len < 0 && throw(DomainError(len, "len must be >= 0"))
   r = parent(x)()
   fit!(r, len)
   for i = 1:len
      z = setcoeff!(r, i - 1, coeff(x, len - i))
   end
   r = set_length!(r, normalise(r, len))
   return r
end

@doc Markdown.doc"""
    reverse(x::PolynomialElem)

Return the reverse of the polynomial $x$, i.e. the leading coefficient
of $x$ becomes the constant coefficient of the result, etc. The resulting
polynomial is normalised.
"""
function reverse(x::PolynomialElem)
   reverse(x, length(x))
end

###############################################################################
#
#   Shifting
#
###############################################################################

@doc Markdown.doc"""
    shift_left(f::PolynomialElem, n::Int)

Return the polynomial $f$ shifted left by $n$ terms, i.e. multiplied by
$x^n$.
"""
function shift_left(f::PolynomialElem, n::Int)
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   if n == 0
      return f
   end
   flen = length(f)
   r = parent(f)()
   fit!(r, flen + n)
   for i = 1:n
      r = setcoeff!(r, i - 1, zero(base_ring(f)))
   end
   for i = 1:flen
      r = setcoeff!(r, i + n - 1, coeff(f, i - 1))
   end
   return r
end

@doc Markdown.doc"""
    shift_right(f::PolynomialElem, n::Int)

Return the polynomial $f$ shifted right by $n$ terms, i.e. divided by
$x^n$.
"""
function shift_right(f::PolynomialElem, n::Int)
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   flen = length(f)
   if n >= flen
      return zero(parent(f))
   end
   if n == 0
      return f
   end
   r = parent(f)()
   fit!(r, flen - n)
   for i = 1:flen - n
      r = setcoeff!(r, i - 1, coeff(f, i + n - 1))
   end
   return r
end

###############################################################################
#
#   Inflation/deflation
#
###############################################################################

@doc Markdown.doc"""
    function deflation(p::PolyElem)

Return a tuple `(shift, defl)` where `shift` is the exponent of the trailing
term of $p$ and `defl` is the gcd of the distance between the exponents of the
nonzero terms of $p$. If $p = 0$, both `shift` and `defl` will be zero.
"""
function deflation(p::PolyElem)
   if iszero(p)
      return 0, 1
   end
   shift = 0
   while iszero(coeff(p, shift))
      shift += 1
   end
   defl = 0
   if shift == degree(p)
      return shift, 1
   end
   for exj in shift + 1:degree(p)
      if !iszero(coeff(p, exj))
         defl = gcd(defl, exj - shift)
         if defl == 1
	    break
	 end
      end
   end
   return shift, defl
end

@doc Markdown.doc"""
    inflate(f::PolyElem, shift::Int64, n::Int64) -> PolyElem

Given a polynomial $f$ in $x$, return $f(x^n)*x^j$, i.e. multiply
all exponents by $n$ and shift $f$ left by $j$.
"""
function inflate(f::PolyElem, j::Int64, n::Int64)
    y = parent(f)()
    for i = 0:degree(f)
        y = setcoeff!(y, n*i + j, coeff(f, i))
    end
    return y
end

@doc Markdown.doc"""
    inflate(f::PolyElem, n::Int64) -> PolyElem

Given a polynomial $f$ in $x$, return $f(x^n)$, i.e. multiply
all exponents by $n$.
"""
inflate(f::PolyElem, n::Int64) = inflate(f, 0, n)

@doc Markdown.doc"""
    deflate(f::PolyElem, shift::Int64, n::Int64) -> PolyElem

Given a polynomial $g$ in $x^n$ such that `f = g(x)*x^{shift}`, write $f$ as
a polynomial in $x$, i.e. divide all exponents of $g$ by $n$.
"""
function deflate(f::PolyElem, j::Int64, n::Int64)
    y = parent(f)()
    for i = 0:div(degree(f) - j, n)
        y = setcoeff!(y, i, coeff(f, n*i + j))
    end
    return y
end

@doc Markdown.doc"""
    deflate(f::PolyElem, n::Int64) -> PolyElem

Given a polynomial $f$ in $x^n$, write it as a polynomial in $x$, i.e. divide
all exponents by $n$.
"""
deflate(f::PolyElem, n::Int64) = deflate(f, 0, n)

@doc Markdown.doc"""
    deflate(x::PolyElem) -> PolyElem, Int

Deflate the polynomial $f$ maximally, i.e. find the largest $n$ s.th.
$f$ can be deflated by $n$, i.e. $f$ is actually a polynomial in $x^n$.
Return $g, n$ where $g$ is the deflation of $f$.
"""
function deflate(f::PolyElem)
   n = 0
   for i = 0:degree(f)
      if coeff(f, i) != 0
         n = gcd(n, i)
         if n == 1
            return f, 1
         end
      end
   end
   return deflate(f, n), n
end

###############################################################################
#
#   Modular arithmetic
#
###############################################################################

@doc Markdown.doc"""
    mulmod(a::PolyElem{T}, b::PolyElem{T}, d::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}

Return $a\times b \pmod{d}$.
"""
function mulmod(a::PolyElem{T}, b::PolyElem{T}, d::PolyElem{T}) where T <: RingElement
   check_parent(a, b)
   check_parent(a, d)
   return mod(a*b, d)
end

@doc Markdown.doc"""
    powermod(a::PolyElem{T}, b::Int, d::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}

Return $a^b \pmod{d}$. There are no restrictions on $b$.
"""
function powermod(a::PolyElem{T}, b::Int, d::PolyElem{T}) where T <: RingElement
   check_parent(a, d)
   if length(a) == 0
      z = zero(parent(a))
   elseif length(a) == 1
      z = parent(a)(coeff(a, 0)^b)
   elseif b == 0
      z = one(parent(a))
   else
      if b < 0
         a = invmod(a, d)
         b = -b
      end
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit !=0
         z = mulmod(z, z, d)
         if (UInt(bit) & b) != 0
            z = mulmod(z, a, d)
         end
         bit >>= 1
      end
   end
   if length(z) >= length(d)
      z = mod(z, d)
   end
   return z
end

@doc Markdown.doc"""
    invmod(a::PolyElem{T}, b::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}

Return $a^{-1} \pmod{d}$.
"""
function invmod(a::PolyElem{T}, b::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}
   check_parent(a, b)
   g, z = gcdinv(a, b)
   if g != 1
      error("Impossible inverse in invmod")
   end
   return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
   check_parent(f, g)
   iszero(g) && throw(DivideError())
   if iszero(f)
      return zero(parent(f))
   end
   lenq = length(f) - length(g) + 1
   d = Array{T}(undef, lenq)
   for i = 1:lenq
      d[i] = zero(base_ring(f))
   end
   x = gen(parent(f))
   leng = length(g)
   while length(f) >= leng
      lenf = length(f)
      q1 = d[lenf - leng + 1] = divexact(coeff(f, lenf - 1), coeff(g, leng - 1))
      f = f - shift_left(q1*g, lenf - leng)
      if length(f) == lenf # inexact case
         f = set_length!(f, normalise(f, lenf - 1))
      end
   end
   length(f) != 0 && throw(ArgumentError("not an exact division"))
   q = parent(f)(d)
   q = set_length!(q, lenq)
   return q
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(a::PolyElem{T}, b::T) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact(coeff(a, i - 1), b))
   end
   z = set_length!(z, length(a))
   return z
end

function divexact(a::PolyElem, b::Union{Integer, Rational, AbstractFloat})
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact(coeff(a, i - 1), b))
   end
   z = set_length!(z, length(a))
   return z
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

@doc Markdown.doc"""
    mod(f::PolyElem{T}, g::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}

Return $f \pmod{g}$.
"""
function mod(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
   check_parent(f, g)
   if length(g) == 0
      throw(DivideError())
   end
   if length(f) >= length(g)
      f = deepcopy(f)
      b = leading_coefficient(g)
      g = inv(b)*g
      c = base_ring(f)()
      while length(f) >= length(g)
         l = -leading_coefficient(f)
         for i = 1:length(g) - 1
            c = mul!(c, coeff(g, i - 1), l)
            u = coeff(f, i + length(f) - length(g) - 1)
            u = addeq!(u, c)
            f = setcoeff!(f, i + length(f) - length(g) - 1, u)
         end
         f = set_length!(f, normalise(f, length(f) - 1))
      end
   end
   return f
end

function rem(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
  return mod(f, g)
end

@doc Markdown.doc"""
    divrem(f::PolyElem{T}, g::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}

Return a tuple $(q, r)$ such that $f = qg + r$ where $q$ is the euclidean
quotient of $f$ by $g$.
"""
function Base.divrem(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
   check_parent(f, g)
   if length(g) == 0
      throw(DivideError())
   end
   if length(f) < length(g)
      return zero(parent(f)), f
   end
   f = deepcopy(f)
   binv = inv(leading_coefficient(g))
   g = divexact(g, leading_coefficient(g))
   qlen = length(f) - length(g) + 1
   q = parent(f)()
   fit!(q, qlen)
   c = base_ring(f)()
   while length(f) >= length(g)
      q1 = leading_coefficient(f)
      l = -q1
      q = setcoeff!(q, length(f) - length(g), q1*binv)
      for i = 1:length(g) - 1
         c = mul!(c, coeff(g, i - 1), l)
         u = coeff(f, i + length(f) - length(g) - 1)
         u = addeq!(u, c)
         f = setcoeff!(f, i + length(f) - length(g) - 1, u)
      end
      f = set_length!(f, normalise(f, length(f) - 1))
   end
   return q, f
end

@doc Markdown.doc"""
    div(f::PolyElem{T}, g::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}

Return the euclidean quotient of $f$ by $g$.
"""
function Base.div(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
   q, r = divrem(f, g)
   return q
end

##############################################################################
#
#  Ad hoc Euclidean division
#
##############################################################################

function Base.div(f::PolyElem{T}, g::T) where T <: Union{FieldElem, ResElem, AbstractFloat, Rational}
   return div(f, parent(f)(g))
end

###############################################################################
#
#   Pseudodivision
#
###############################################################################

@doc Markdown.doc"""
   pseudorem(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement

Return the pseudoremainder of $f$ divided by $g$. If $g = 0$ we throw a
`DivideError()`.
"""
function pseudorem(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
  check_parent(f, g)
  iszero(g) && throw(DivideError())
  if length(f) < length(g)
     return f
  end
  k = length(f) - length(g) + 1
  b = coeff(g, length(g) - 1)
  x = gen(parent(f))
  while length(f) >= length(g)
     f = f*b - shift_left(coeff(f, length(f) - 1)*g, length(f) - length(g))
     k -= 1
  end
  return f*b^k
end

@doc Markdown.doc"""
   pseudodivrem(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement

Return a tuple $(q, r)$ consisting of the pseudoquotient and pseudoremainder
of $f$ divided by $g$. If $g = 0$ we throw a `DivideError()`.
"""
function pseudodivrem(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
  check_parent(f, g)
  iszero(g) && throw(DivideError())
  if length(f) < length(g)
     return zero(parent(f)), f
  end
  lenq = length(f) - length(g) + 1
  k = lenq
  q = parent(f)()
  fit!(q, lenq)
  b = coeff(g, length(g) - 1)
  x = gen(parent(f))
  while length(f) >= length(g)
     for i = length(f) - length(g) + 2:lenq
        q = setcoeff!(q, i - 1, coeff(q, i - 1) * b)
     end
     q = setcoeff!(q, length(f) - length(g), coeff(f, length(f) - 1))
     f = f*b - shift_left(coeff(f, length(f) - 1)*g, length(f) - length(g))
     k -= 1
  end
  while lenq > 0 && iszero(coeff(q, lenq - 1))
     lenq -= 1
  end
  q = set_length!(q, lenq)
  s = b^k
  return q*s, f*s
end

################################################################################
#
#   Remove and valuation
#
################################################################################

#CF TODO: use squaring for fast large valuation

@doc Markdown.doc"""
   remove(z::PolyElem{T}, p::PolyElem{T}) where T <: RingElement

Compute the valuation of $z$ at $p$, that is, the largest $k$ such that
$p^k$ divides $z$. Additionally, $z/p^k$ is returned as well.

See also `valuation`, which only returns the valuation.
"""
function remove(z::PolyElem{T}, p::PolyElem{T}) where T <: RingElement
 check_parent(z, p)
 !isexact_type(T) && error("remove requires an exact ring")
 iszero(z) && error("Not yet implemented")
 (isunit(p) || iszero(p)) && throw(error("Second argument must be a non-zero non-unit"))
 flag, q = divides(z, p)
 if !flag
   return 0, z
 end
 v = 0
 qn = q
 while flag
   q = qn
   flag, qn = divides(q, p)
   v += 1
 end
 return v, q
end

@doc Markdown.doc"""
   remove(z::PolyElem{T}, p::PolyElem{T}) where T <: Union{ResElem, FieldElement}

Compute the valuation of $z$ at $p$, that is, the largest $k$ such that
$p^k$ divides $z$. Additionally, $z/p^k$ is returned as well.

See also `valuation`, which only returns the valuation.
"""
function remove(z::PolyElem{T}, p::PolyElem{T}) where T <: Union{ResElem, FieldElement}
 check_parent(z, p)
 !isexact_type(T) && error("remove requires an exact ring")
 iszero(z) && error("Not yet implemented")
 (isunit(p) || iszero(p)) && throw(error("Second argument must be a non-zero non-unit"))
 q, r = divrem(z, p)
 if !iszero(r)
   return 0, z
 end
 v = 0
 qn = q
 while iszero(r)
   q = qn
   qn, r = divrem(q, p)
   v += 1
 end
 return v, q
end

@doc Markdown.doc"""
   valuation(z::PolyElem{T}, p::PolyElem{T}) where T <: RingElement

Compute the valuation of $z$ at $p$, that is, the largest $k$ such that
$p^k$ divides $z$.

See also `remove`, which also returns $z/p^k$.
"""
function valuation(z::PolyElem{T}, p::PolyElem{T}) where T <: RingElement
 v, _ = remove(z, p)
 return v
end

@doc Markdown.doc"""
   divides(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement

Return a pair consisting of a flag which is set to `true` if $g$ divides
$f$ and `false` otherwise, and a polynomial $h$ such that $f = gh$ if
such a polynomial exists. If not, the value of $h$ is undetermined.
"""
function divides(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
  check_parent(f, g)
  !isexact_type(T) && error("divides requires an exact ring")
  if length(f) == 0
     return true, parent(f)()
  end
  if length(g) == 0
     return false, parent(f)()
  end
  if length(f) < length(g)
     return false, parent(f)()
  end
  f = deepcopy(f)
  g_lead = leading_coefficient(g)
  qlen = length(f) - length(g) + 1
  q = parent(f)()
  fit!(q, qlen)
  c = base_ring(f)()
  while length(f) >= length(g)
     q1 = leading_coefficient(f)
     flag, d = divides(q1, g_lead)
     if !flag
        return false, parent(f)()
     end
     q = setcoeff!(q, length(f) - length(g), d)
     d = -d
     for i = 1:length(g)
        c = mul!(c, coeff(g, i - 1), d)
        u = coeff(f, i + length(f) - length(g) - 1)
        u = addeq!(u, c)
        f = setcoeff!(f, i + length(f) - length(g) - 1, u)
     end
     f = set_length!(f, normalise(f, length(f)))
  end
  return iszero(f), q
end

@doc Markdown.doc"""
   divides(z::PolyElem{T}, x::T) where T <: RingElement

Return a pair consisting of a flag which is set to `true` if $x$ divides
$z$ and `false` otherwise, and a polynomial $y$ such that $z = xy$ if
such a polynomial exists. If not, the value of $y$ is undetermined.
"""
function divides(z::PolyElem{T}, x::T) where T <: RingElement
  parent(x) != base_ring(z) && error("Wrong parents in divides")
  q = parent(z)()
  fit!(q, length(z))
  flag = true
  for i = 1:length(z)
     flag, c = divides(coeff(z, i - 1), x)
     if !flag
        break
     end
     q = setcoeff!(q, i - 1, c)
  end
  q = set_length!(q, flag ? length(z) : 0)
  return flag, q
end

################################################################################
#
#   Square root
#
################################################################################

function sqrt_classical_char2(f::PolyElem{T}) where T <: RingElement
   S = parent(f)
   R = base_ring(f)
   if iszero(f)
      return true, S()
   end
   m = length(f)
   if iseven(m) # square polys have even degree
      return false, S()
   end
   for i = 1:2:m # polynomial must have even exponents
      if !iszero(coeff(f, i))
         return false, S()
      end
   end
   lenq = div(m + 1, 2)
   d = Array{T}(undef, lenq)
   for i = 1:lenq
      c = coeff(f, 2*i - 2)
      if !issquare(c)
         return false, S()
      end
      d[i] = sqrt(c)
   end
   q = S(d)
   q = set_length!(q, lenq)
   return true, q
end

function sqrt_classical(f::PolyElem{T}, check::Bool=true) where T <: RingElement
   S = parent(f)
   R = base_ring(f)
   if characteristic(R) == 2
      return sqrt_classical_char2(f)
   end
   if iszero(f)
      return true, S()
   end
   m = length(f)
   if iseven(m) # square polys have even degree
      return false, S()
   end
   if !issquare(coeff(f, m - 1))
      return false, S()
   end
   lenq = div(m + 1, 2)
   d = Array{T}(undef, lenq)
   d[lenq] = sqrt(coeff(f, m - 1))
   b = -2*d[lenq]
   k = 1
   c = R()
   last = check ? 0 : lenq - 1
   for i = m - 2:-1:last
      qc = -coeff(f, i)
      for j = lenq - k + 1:lenq
         if i - j + 2 >= j && j > 0
            c = mul_red!(c, d[j], d[i - j + 2], false)
            qc = addeq!(qc, c)
            if (j != i - j + 2)
               qc = addeq!(qc, c)
            end
         end
      end
      qc = reduce!(qc)
      if i >= lenq - 1
         flag, d[lenq - k] = divides(qc, b)
         if !flag
            return false, S()
         end
      elseif !iszero(qc)
         return false, S()
      end
      k += 1
   end
   q = S(d)
   q = set_length!(q, lenq)
   return true, q
end

@doc Markdown.doc"""
    Base.sqrt(f::PolyElem{T}, check::Bool=true) where T <: RingElement

Return the square root of $f$ if it is a perfect square, otherwise an
exception is raised. If `check` is set to `false` the function assumes
the input is square and may not fully check this.
"""
function Base.sqrt(f::PolyElem{T}, check::Bool=true) where T <: RingElement
   flag, q = sqrt_classical(f, check)
   !flag && error("Not a square in sqrt")
   return q
end

@doc Markdown.doc"""
    issquare(f::PolyElem{T}) where T <: RingElement

Return `true` if $f$ is a perfect square.
"""
function issquare(f::PolyElem{T}) where T <: RingElement
   flag, q = sqrt_classical(f)
   return flag
end

###############################################################################
#
#   Content, primitive part, GCD and LCM
#
###############################################################################

function term_gcd(a::T, b::T) where T <: RingElement
   return gcd(a, b)
end

function term_content(a::T) where T <: RingElement
   return a
end

function term_gcd(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   d = min(degree(a), degree(b))
   x = gen(parent(a))
   return term_gcd(coeff(a, degree(a)), coeff(b, degree(b)))*x^d
end

function term_content(a::PolyElem{T}) where T <: RingElement
   for i = 1:length(a)
      c = coeff(a, i - 1)
      if !iszero(c)
         g = term_content(c)
         for j = i + 1:length(a)
            c = coeff(a, j - 1)
            if !iszero(c)
               g = term_gcd(g, term_content(c))
            end
         end
         x = gen(parent(a))
         return g*x^(i - 1)
      end
   end
   return parent(a)()
end

@doc Markdown.doc"""
    gcd(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement

Return a greatest common divisor of $a$ and $b$ if it exists.
"""
function gcd(a::PolyElem{T}, b::PolyElem{T}, ignore_content::Bool = false) where T <: RingElement
   check_parent(a, b)
   if length(b) > length(a)
      (a, b) = (b, a)
   end
   if iszero(b)
      if iszero(a)
         return a
      else
         return divexact(a, canonical_unit(leading_coefficient(a)))
      end
   end
   if isone(b)
      return b
   end
   if !ignore_content
      c1 = content(a)
      c2 = content(b)
      a = divexact(a, c1)
      b = divexact(b, c2)
      c = gcd(c1, c2)
   end
   lead_monomial = isterm_recursive(leading_coefficient(a)) ||
                   isterm_recursive(leading_coefficient(b))
   trail_monomial = isterm_recursive(trailing_coefficient(a)) ||
                    isterm_recursive(trailing_coefficient(b))
   lead_a = leading_coefficient(a)
   lead_b = leading_coefficient(b)
   g = one(parent(a))
   h = one(parent(a))
   while true
      d = length(a) - length(b)
      r = pseudorem(a, b)
      if r == 0
         break
      end
      if length(r) == 1
         b = one(parent(a))
         break
      end
      (a, b) = (b, divexact(r, g*h^d))
      g = leading_coefficient(a)
      if d > 1
         h = divexact(g^d, h^(d - 1))
      else
         h = h^(1 - d)*g^d
      end
   end
   if !ignore_content
      if !isterm_recursive(leading_coefficient(b)) &&
         !isterm_recursive(trailing_coefficient(b))
         if lead_monomial # lead term monomial, so content contains rest
            d = divexact(leading_coefficient(b),
	                 term_content(leading_coefficient(b)))
            b = divexact(b, d)
         elseif trail_monomial # trail term is monomial, so ditto
            d = divexact(trailing_coefficient(b),
			 term_content(trailing_coefficient(b)))
            b = divexact(b, d)
         else
            glead = gcd(lead_a, lead_b)
            if isterm_recursive(glead)
               d = divexact(leading_coefficient(b),
			     term_content(leading_coefficient(b)))
               b = divexact(b, d)
            else # last ditched attempt to find easy content
               h = gcd(leading_coefficient(b), glead)
               h = divexact(h, term_content(h))
               flag, q = divides(b, h)
               if flag
                  b = q
               end
            end
         end
      end
      b = divexact(b, content(b))
   end
   b = c*b
   return divexact(b, canonical_unit(leading_coefficient(b)))
end

function gcd(a::PolyElem{T}, b::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}
   check_parent(a, b)
   if length(a) > length(b)
      (a, b) = (b, a)
   end
   if iszero(b)
      if iszero(a)
         return(a)
      else
         d = leading_coefficient(a)
         return divexact(a, d)
      end
   end
   g = gcd(content(a), content(b))
   a = divexact(a, g)
   b = divexact(b, g)
   while !iszero(a)
      (a, b) = (mod(b, a), a)
   end
   b = g*b
   d = leading_coefficient(b)
   return divexact(b, d)
end

@doc Markdown.doc"""
    lcm(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement

Return a least common multiple of $a$ and $b$ if it exists.
"""
function lcm(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   check_parent(a, b)
   return a*divexact(b, gcd(a, b))
end

@doc Markdown.doc"""
    content(a::PolyElem)

Return the content of $a$, i.e. the greatest common divisor of its
coefficients.
"""
function content(a::PolyElem)
   z = base_ring(a)() # normalise first coefficient
   for i = 1:length(a)
      z = gcd(z, coeff(a, i - 1))
   end
   return z
end

@doc Markdown.doc"""
    primpart(a::PolyElem)

Return the primitive part of $a$, i.e. the polynomial divided by its content.
"""
function primpart(a::PolyElem)
   d = content(a)
   if iszero(d)
      return zero(parent(a))
   else
      return divexact(a, d)
   end
end

###############################################################################
#
#   Evaluation/composition
#
###############################################################################

@doc Markdown.doc"""
    evaluate(a::PolyElem, b::T) where T <: RingElement

Evaluate the polynomial expression $a$ at the value $b$ and return the result.
"""
function evaluate(a::PolyElem, b::T) where T <: RingElement
   i = length(a)
   R = base_ring(a)
   S = parent(b)
   if i == 0
      return zero(R) + zero(S)
   end
   if i > 25
      return subst(a, b)
   end
   z = coeff(a, i - 1) * one(S)
   while i > 1
      i -= 1
      z = coeff(a, i - 1) + z*b
      parent(z) # To work around a bug in julia
   end
   return z
end

@doc Markdown.doc"""
    compose(a::PolyElem, b::PolyElem)

Compose the polynomial $a$ with the polynomial $b$ and return the result,
i.e. return $a\circ b$.
"""
function compose(a::PolyElem, b::PolyElem)
   i = length(a)
   R = base_ring(a)
   S = parent(b)
   if i == 0
      return zero(R) + zero(S)
   end
   if i*length(b) > 25
      return subst(a, b)
   end
   z = coeff(a, i - 1) * one(S)
   while i > 1
      i -= 1
      z = z*b + coeff(a, i - 1)
   end
   return z
end

###############################################################################
#
#   Derivative
#
###############################################################################

@doc Markdown.doc"""
    derivative(a::PolynomialElem)

Return the derivative of the polynomial $a$.
"""
function derivative(a::PolynomialElem)
   if iszero(a)
      return zero(parent(a))
   end
   len = length(a)
   z = parent(a)()
   fit!(z, len - 1)
   for i = 1:len - 1
      z = setcoeff!(z, i - 1, i*coeff(a, i))
   end
   z = set_length!(z, normalise(z, len - 1))
   return z
end

###############################################################################
#
#   Integral
#
###############################################################################

@doc Markdown.doc"""
    integral(x::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}

Return the integral of the polynomial $x$.
"""
function integral(x::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}
   len = length(x)
   p = parent(x)()
   fit!(p, len + 1)
   p = setcoeff!(p, 0, zero(base_ring(x)))
   for i = 1:len
      p = setcoeff!(p, i, divexact(coeff(x, i - 1), base_ring(x)(i)))
   end
   len += 1
   while len > 0 && iszero(coeff(p, len - 1)) # FIXME: cannot use normalise here
      len -= 1
   end
   p = set_length!(p, len)
   return p
end

###############################################################################
#
#   Resultant
#
###############################################################################

# Dichotomous Lazard, computes Se0 from Sd0 and Sd1. See the paper,
# "Optimizations of the subresultant algorithm" by Lionel Ducos, J. Pure and
# Appl. Algebra 2000.
function subresultant_lazard(Sd0::PolyElem{T}, Sd1::PolyElem{T}) where T <: RingElement
   n = length(Sd0) - length(Sd1) - 1
   if n == 0
      return Sd1
   end
   x = leading_coefficient(Sd1)
   y = leading_coefficient(Sd0)
   a = 1 << (ndigits(n, base = 2) - 1) # floor(log_2(a))
   c = x
   n = n - a
   while a != 1
      a = a >> 1
      c = divexact(c*c, y)
      if n >= a
         c = divexact(c*x, y)
         n -= a
      end
   end
   return divexact(c*Sd1, y)
end

# Ducos optimised calculation of Se1. See the paper, "Optimizations of the
# subresultant algorithm" by Lionel Ducos, J. Pure and Appl. Algebra 2000.
function subresultant_ducos(A::PolyElem{T}, Sd1::PolyElem{T}, Se0::PolyElem{T}, sd::T) where T <: RingElement
   d1 = length(A)
   e1 = length(Sd1)
   cd1 = leading_coefficient(Sd1)
   se = leading_coefficient(Se0)
   D = parent(A)()
   fit!(D, d1 - 1)
   for j = 0:e1 - 2
      setcoeff!(D, j, se*coeff(A, j))
   end
   D = set_length!(D, normalise(D, e1 - 1))
   Hj = parent(A)()
   fit!(Hj, e1)
   setcoeff!(Hj, e1 - 1, se)
   Hj -= Se0
   D += coeff(A, e1 - 1)*Hj
   for j = e1:d1 - 2
      Hj = shift_left(Hj, 1)
      Hj -= divexact(coeff(Hj, e1 - 1)*Sd1, cd1)
      D += coeff(A, j)*Hj
   end
   D = divexact(D, leading_coefficient(A))
   Hj = shift_left(Hj, 1)
   r = divexact((Hj + D)*cd1 - coeff(Hj, e1 - 1)*Sd1, sd)
   return iseven(d1 - e1) ? -r : r
end

@doc Markdown.doc"""
    resultant_ducos(p::PolyElem{T}, q::PolyElem{T}) where T <: RingElement

Return the resultant of the $p$ and $q$.
"""
function resultant_ducos(p::PolyElem{T}, q::PolyElem{T}) where T <: RingElement
   # See the paper, "Optimizations of the subresultant algorithm" by Lionel
   # Ducos, J. Pure and Appl. Algebra 2000.
   check_parent(p, q)
   if length(p) == 0 || length(q) == 0
      return zero(base_ring(p))
   end
   sgn = 1
   if length(p) < length(q)
      p, q = q, p
      if iseven(length(p)) && iseven(length(q))
         sgn = -sgn
      end
   end
   lp = length(p)
   lq = length(q)
   if lq == 1
      return coeff(q, 0)^(lp - 1)
   end
   c1 = content(p)
   c2 = content(q)
   p = divexact(p, c1)
   q = divexact(q, c2)
   sd = leading_coefficient(q)^(lp - lq)
   Sd0 = parent(p)()
   A = q
   B = pseudorem(p, -A)
   while true
      d1 = length(A)
      e1 = length(B)
      if e1 == 0
         return zero(base_ring(p))
      end
      Sd1 = B
      delta = d1 - e1
      if delta > 1
         if length(Sd0) == 0
            C = divexact(leading_coefficient(B)^(delta - 1)*B, sd^(delta - 1))
         else
            C = subresultant_lazard(Sd0, Sd1)
         end
      else
         C = B
      end
      if e1 == 1
         return coeff(C, 0)*c1^(lq - 1)*c2^(lp - 1)*sgn
      end
      B = subresultant_ducos(A, Sd1, C, sd)
      Sd0 = C
      Sd1 = B
      A = C
      sd = leading_coefficient(A)
   end
end

# details can be found in, "Optimizations of the subresultant algorithm" by
# Lionel Ducos, J. Pure and Appl. Algebra 2000. Note, the resultant is
# the constant coefficient of S_0 (aka S_00 in other sources)
function resultant_subresultant(p::PolyElem{T}, q::PolyElem{T}) where T <: RingElement
   check_parent(p, q)
   if length(p) == 0 || length(q) == 0
      return zero(base_ring(p))
   end
   sgn = 1
   if length(p) < length(q)
      p, q = q, p
      if iseven(length(p)) && iseven(length(q))
         sgn = -sgn
      end
   end
   lp = length(p)
   lq = length(q)
   if lq == 1
      return coeff(q, 0)^(lp - 1)
   end
   s = leading_coefficient(q)^(lp - lq)
   S = parent(p)()
   A = q
   B = pseudorem(p, -q)
   while true
      d1 = length(A)
      e1 = length(B)
      if e1 == 0
         return zero(base_ring(p))
      end
      S = B
      delta = d1 - e1
      if delta > 1
         C = divexact(leading_coefficient(B)^(delta - 1)*B, s^(delta - 1))
         S = C
      else
         C = B
      end
      if e1 == 1
         return coeff(S, 0)*sgn
      end
      B = divexact(pseudorem(A, -B), s^delta*leading_coefficient(A))
      A = C
      s = leading_coefficient(A)
   end
end

function resultant_lehmer(a::PolyElem{T}, b::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}
   local crossover = 40
   R = base_ring(a)
   check_parent(a, b)
   if length(a) == 0 || length(b) == 0
      return zero(base_ring(a))
   end
   sgn = 1
   if length(a) < length(b)
      a, b = b, a
      if iseven(length(a)) && iseven(length(b))
         sgn = -sgn
      end
   end
   lA = lenA = length(a)
   lB = lenB = length(b)
   if lenB == 1
      return coeff(b, 0)^(lA - 1)
   end
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   s = R(1)
   while lenB > crossover/2 + 1
      shift = max(lenA - crossover, 0)
      a = shift_right(A, shift)
      b = shift_right(B, shift)
      u1, v1 = R(1), R(0)
      u2, v2 = R(0), R(1)
      lena = lenA - shift
      lenb = lenB - shift
      if lenb > crossover/2 + 1
         A = truncate(A, shift)
         B = truncate(B, shift)
         while lenb > crossover/2 + 1
            if iseven(lena + shift) && iseven(lenb + shift)
               sgn = -sgn
            end
            (q, b), a = divrem(a, b), b
            u1, u2 = u2, u1 - q*u2
            v1, v2 = v2, v1 - q*v2
            s *= leading_coefficient(a)^(lena - length(b))
            lena = lenb
            lenb = length(b)
         end
         A, B = u1*A + v1*B + shift_left(a, shift), u2*A + v2*B + shift_left(b, shift)
      else
         if iseven(lenA) && iseven(lenB)
               sgn = -sgn
         end
         B, A = mod(A, B), B
         s *= leading_coefficient(A)^(lenA - length(B))
      end
      lenA = length(A)
      lenB = length(B)
      if lenB == 0
         return zero(base_ring(a)), parent(A)(1), parent(A)(1)
      end
   end
   while lenB > 1
      if iseven(lenA) && iseven(lenB)
         sgn = -sgn
      end
      B, A = mod(A, B), B
      s *= leading_coefficient(A)^(lenA - length(B))
      lenA = lenB
      lenB = length(B)
      if lenB == 0
         return zero(base_ring(A))
      end
   end
   s *= leading_coefficient(B)^(lenA - 1)
   return c1^(lB - 1)*c2^(lA - 1)*s*sgn
end

@doc Markdown.doc"""
    sylvester_matrix(p::PolyElem, q::PolyElem)

Return the sylvester matrix of the given polynomials.
"""
function sylvester_matrix(p::PolyElem{T}, q::PolyElem{T}) where T <: RingElement
   check_parent(p, q)
   R = base_ring(p)
   if length(p) == 0 || length(q) == 0
      return zero_matrix(R, 0, 0)
   end
   m = degree(p)
   n = degree(q)
   M = zero_matrix(R, m + n, m + n)
   for i = 1:n
      for j = m:-1:0
         M[i, m - j + i] = coeff(p, j)
      end
   end
   for i = 1:m
      for j = n:-1:0
         M[i + n, n - j + i] = coeff(q, j)
       end
   end
   return M
end

function resultant_sylvester(p::PolyElem{T}, q::PolyElem{T}) where T <: RingElement
   check_parent(p, q)
   R = base_ring(p)
   if length(p) == 0 || length(q) == 0
      return R(0)
   end
   return det_df(sylvester_matrix(p, q))
end

@doc Markdown.doc"""
    resultant(p::PolyElem{T}, q::PolyElem{T}) where T <: RingElement

Return the resultant of the given polynomials.
"""
function resultant(p::PolyElem{T}, q::PolyElem{T}) where T <: RingElement
  R = parent(p)
  if !isexact_type(T)
     return resultant_sylvester(p, q)
  end
  try
     return resultant_ducos(p, q)
  catch
     return resultant_sylvester(p, q)
  end
end

function resultant_euclidean(a::PolyElem{T}, b::PolyElem{T}) where T <: Union{ResElem, FieldElement}
   check_parent(a, b)
   if length(a) == 0 || length(b) == 0
      return zero(base_ring(a))
   end
   sgn = 1
   if length(a) < length(b)
      a, b = b, a
      if iseven(length(a)) && iseven(length(b))
         sgn = -sgn
      end
   end
   la = lena = length(a)
   lb = lenb = length(b)
   if lenb == 1
      return coeff(b, 0)^(la - 1)
   end
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   s = base_ring(A)(1)
   la = lena = length(A)
   lb = lenb = length(B)
   while lenb > 1
      if iseven(lena) && iseven(lenb)
         sgn = -sgn
      end
      B, A = mod(A, B), B
      s *= leading_coefficient(A)^(lena - length(B))
      parent(s) # julia bug
      lena = lenb
      lenb = length(B)
      if lenb == 0
         return zero(base_ring(a))
      end
   end
   s *= leading_coefficient(B)^(lena - 1)
   return c1^(lb - 1)*c2^(la - 1)*s*sgn
end

function resultant(a::PolyElem{T}, b::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}
   try
      return resultant_euclidean(a, b)
   catch
      return resultant_sylvester(a, b)
   end
end

###############################################################################
#
#   Discriminant
#
###############################################################################

@doc Markdown.doc"""
    discriminant(a::PolyElem)

Return the discriminant of the given polynomial.
"""
function discriminant(a::PolyElem)
   d = derivative(a)
   z = resultant(a, d)
   if length(a) - length(d) == 1
      z = divexact(z, leading_coefficient(a))
   else
      z = z*leading_coefficient(a)^(length(a) - length(d) - 2)
   end
   mod4 = (length(a) + 3) % 4 # degree mod 4
   return mod4 == 2 || mod4 == 3 ? -z : z
end

###############################################################################
#
#   RESX
#
###############################################################################

@doc Markdown.doc"""
    resx(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement

Return a tuple $(r, s, t)$ such that $r$ is the resultant of $a$ and $b$ and
such that $r = a\times s + b\times t$.
"""
function resx(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   check_parent(a, b)
   sgn = 1
   swap = false
   if length(a) < length(b)
      a, b = b, a
      swap = true
      if iseven(length(a)) && iseven(length(b))
         sgn = -sgn
      end
   end
   lena = length(a)
   lenb = length(b)
   if lenb == 0
      return zero(base_ring(a)), zero(parent(a)), zero(parent(a))
   end
   (lena <= 1 && lenb <= 1) && error("Constant polynomials in resx")
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   g = one(base_ring(a))
   h = one(base_ring(a))
   u1, u2 = one(parent(a)), zero(parent(a))
   v1, v2 = zero(parent(a)), one(parent(a))
   while lenb > 1
      d = lena - lenb
      if iseven(lena) && iseven(lenb)
         sgn = -sgn
      end
      (Q, B), A = pseudodivrem(A, B), B
      lena = lenb
      lenb = length(B)
      if lenb == 0
         return zero(base_ring(a)), zero(parent(a)), zero(parent(a))
      end
      s = h^d
      B = divexact(B, g*s)
      t = leading_coefficient(A)^(d + 1)
      u2, u1 = divexact(u1*t - Q*u2, g*s), u2
      v2, v1 = divexact(v1*t - Q*v2, g*s), v2
      g = leading_coefficient(A)
      h = divexact(h*g^d, s)
   end
   s = divexact(h*leading_coefficient(B)^(lena - 1), h^(lena - 1))
   res = c1^(length(b) - 1)*c2^(length(a) - 1)*s*sgn
   if length(b) > 1
      u2 *= c1^(length(b) - 2)*c2^(length(a) - 1)*sgn
   else
      u2 *= c2^(length(a) - 1)*sgn
      u2 = divexact(u2, c1)
   end
   if length(a) > 1
      v2 *= c1^(length(b) - 1)*c2^(length(a) - 2)*sgn
   else
      v2 *= c1^(length(b) - 1)*sgn
      v2 = divexact(v2, c2)
   end
   if lena != 2
      if lena > 1
         d1 = leading_coefficient(B)^(lena - 2)
         d2 = h^(lena - 2)
         u2 = divexact(u2*d1, d2)
         v2 = divexact(v2*d1, d2)
      else
         u2 = divexact(u2*h, leading_coefficient(B))
         v2 = divexact(v2*h, leading_coefficient(B))
      end
   end
   if swap
      u2, v2 = v2, u2
   end
   return res, u2, v2
end

###############################################################################
#
#   GCDX
#
###############################################################################

@doc Markdown.doc"""
    gcdx(a::PolyElem{T}, b::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}

Return a tuple $(g, s, t)$ such that $g$ is the greatest common divisor of
$a$ and $b$ and such that $g = a\times s + b\times t$.
"""
function gcdx(a::PolyElem{T}, b::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}
   check_parent(a, b)
   !isexact_type(T) && error("gcdx requires exact Bezout domain")
   if length(a) == 0
      if length(b) == 0
         return zero(parent(a)), zero(parent(a)), zero(parent(a))
      else
         d = leading_coefficient(b)
         return divexact(b, d), zero(parent(a)), divexact(one(parent(a)), d)
      end
   end
   if length(b) == 0
      d = leading_coefficient(a)
      return divexact(a, d), divexact(one(parent(a)), d), zero(parent(a))
   end
   swap = false
   if length(a) < length(b)
      a, b = b, a
      swap = true
   end
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   u1, u2 = parent(a)(inv(c1)), zero(parent(a))
   v1, v2 = zero(parent(a)), parent(a)(inv(c2))
   while length(B) > 0
      (Q, B), A = divrem(A, B), B
      u2, u1 = u1 - Q*u2, u2
      v2, v1 = v1 - Q*v2, v2
   end
   if swap
      u1, v1 = v1, u1
   end
   d = gcd(c1, c2)
   A, u1, v1 = d*A, d*u1, d*v1
   d = leading_coefficient(A)
   return divexact(A, d), divexact(u1, d), divexact(v1, d)
end

@doc Markdown.doc"""
    gcdinv(a::PolyElem{T}, b::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}

Return a tuple $(g, s)$ such that $g$ is the greatest common divisor of $a$
and $b$ and such that $s = a^{-1} \pmod{b}$. This function is useful for
inverting modulo a polynomial and checking that it really was invertible.
"""
function gcdinv(a::PolyElem{T}, b::PolyElem{T}) where {T <: Union{ResElem, FieldElement}}
   check_parent(a, b)
   R = base_ring(a)
   if length(a) == 0
      if length(b) == 0
         return zero(parent(a)), zero(parent(a))
      else
         d = leading_coefficient(b)
         return divexact(b, d), zero(parent(a))
      end
   end
   if length(b) == 0
      d = leading_coefficient(a)
      return divexact(a, d), parent(a)(inv(d))
   end
   if length(a) < length(b)
      a, b = b, a
      u1, u2 = zero(parent(a)), one(parent(a))
   else
      u1, u2 = one(parent(a)), zero(parent(a))
   end
   lena = length(a)
   lenb = length(b)
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   u1 *= inv(c1)
   u2 *= inv(c2)
   while lenb > 0
      d = lena - lenb
      (Q, B), A = divrem(A, B), B
      lena = lenb
      lenb = length(B)
      u2, u1 = u1 - Q*u2, u2
   end
   d = gcd(c1, c2)
   A, u1 = d*A, d*u1
   d = leading_coefficient(A)
   return divexact(A, d), divexact(u1, d)
end

###############################################################################
#
#   Newton representation
#
###############################################################################

@doc Markdown.doc"""
    monomial_to_newton!(P::Array{T, 1}, roots::Array{T, 1}) where T <: RingElement

Converts a polynomial $p$, given as an array of coefficients, in-place
from its coefficients given in the standard monomial basis to the Newton
basis for the roots $r_0, r_1, \ldots, r_{n-2}$. In other words, this
determines output coefficients $c_i$ such that
$$c_0 + c_1(x-r_0) + c_2(x-r_0)(x-r_1) + \ldots + c_{n-1}(x-r_0)(x-r_1)\cdots(x-r_{n-2})$$
is equal to the input polynomial.
"""
function monomial_to_newton!(P::Array{T, 1}, roots::Array{T, 1}) where T <: RingElement
   n = length(roots)
   if n > 0
      R = parent(roots[1])
      t = R()
      for i = 1:n - 1
         for j = n - 1:-1:i
            t = mul!(t, P[j + 1], roots[i])
            P[j] = addeq!(P[j], t)
         end
      end
   end
   return
end

@doc Markdown.doc"""
    newton_to_monomial!(P::Array{T, 1}, roots::Array{T, 1}) where T <: RingElement

Converts a polynomial $p$, given as an array of coefficients, in-place
from its coefficients given in the Newton basis for the roots
$r_0, r_1, \ldots, r_{n-2}$ to the standard monomial basis. In other words,
this evaluates
$$c_0 + c_1(x-r_0) + c_2(x-r_0)(x-r_1) + \ldots + c_{n-1}(x-r_0)(x-r_1)\cdots(x-r_{n-2})$$
where $c_i$ are the input coefficients given by $p$.
"""
function newton_to_monomial!(P::Array{T, 1}, roots::Array{T, 1}) where T <: RingElement
   n = length(roots)
   if n > 0
      R = parent(roots[1])
      t = R()
      for i = n - 1:-1:1
         d = -roots[i]
         for j = i:n - 1
            t = mul!(t, P[j + 1], d)
            P[j] = addeq!(P[j], t)
         end
      end
   end
   return
end

###############################################################################
#
#   Interpolation
#
###############################################################################

@doc Markdown.doc"""
    interpolate(S::PolyRing, x::Array{T, 1}, y::Array{T, 1}) where T <: RingElement

Given two arrays of values $xs$ and $ys$ of the same length $n$, find
the polynomial $f$ in the polynomial ring $R$ of length at most $n$ such that
$f$ has the value $ys$ at the points $xs$. The values in the arrays $xs$ and
$ys$ must belong to the base ring of the polynomial ring $R$. If no such
polynomial exists, an exception is raised.
"""
function interpolate(S::PolyRing, x::Array{T, 1}, y::Array{T, 1}) where T <: RingElement
   length(x) != length(y) && error("Array lengths don't match in interpolate")
   !isdomain_type(T) && error("interpolation requires a domain type")
   n = length(x)
   if n == 0
      return S()
   elseif n == 1
      return S(y[1])
   end
   R = base_ring(S)
   parent(y[1]) != R && error("Polynomial ring does not match inputs")
   P = Array{T}(undef, n)
   for i = 1:n
      P[i] = deepcopy(y[i])
   end
   for i = 2:n
      t = P[i - 1]
      for j = i:n
         p = P[j] - t
         q = x[j] - x[j - i + 1]
         t = P[j]
         P[j] = divexact(p, q) # division is exact over domain (Lipson, 1971)
      end
   end
   newton_to_monomial!(P, x)
   r = S(P)
   r = set_length!(r, normalise(r, n))
   return r
end

function interpolate(S::PolyRing, x::Array{T, 1}, y::Array{T, 1}) where {T <: ResElem}
   length(x) != length(y) && error("Array lengths don't match in interpolate")
   n = length(x)
   if n == 0
      return S()
   elseif n == 1
      return S(y[1])
   end
   R = base_ring(S)
   parent(y[1]) != R && error("Polynomial ring does not match inputs")
   P = Array{T}(undef, n)
   for i = 1:n
      P[i] = deepcopy(y[i])
   end
   for i = 2:n
      t = P[i - 1]
      for j = i:n
         p = P[j] - t
         q = x[j] - x[j - i + 1]
         t = P[j]
         P[j] = p*inv(q) # must have invertible q for now
      end
   end
   newton_to_monomial!(P, x)
   r = S(P)
   r = set_length!(r, normalise(r, n))
   return r
end

################################################################################
#
#  Change base ring
#
################################################################################

function _change_poly_ring(R, Rx, cached)
   P, _ = PolynomialRing(R, string(var(Rx)), cached = cached)
   return P
end

@doc Markdown.doc"""
    change_base_ring(R::Ring, p::PolyElem{<: RingElement}; parent::PolyRing)

Return the polynomial obtained by coercing the non-zero coefficients of `p`
into `R`.

If the optional `parent` keyword is provided, the polynomial will be an
element of `parent`. The caching of the parent object can be controlled
via the `cached` keyword argument.
"""
function change_base_ring(R::Ring, p::PolyElem{T}; cached::Bool = true, parent::PolyRing = _change_poly_ring(R, parent(p), cached)) where T <: RingElement
   return _map(R, p, parent)
end

################################################################################
#
#  Map
#
################################################################################

_make_parent(g, p::PolyElem, cached::Bool) =
   _change_poly_ring(parent(g(zero(base_ring(p)))),
                     parent(p), cached)

@doc Markdown.doc"""
    map_coefficients(f, p::PolyElem{<: RingElement}; cached::Bool=true, parent::PolyRing)

Transform the polynomial `p` by applying `f` on each non-zero coefficient.

If the optional `parent` keyword is provided, the polynomial will be an
element of `parent`. The caching of the parent object can be controlled
via the `cached` keyword argument.
"""
function map_coefficients(g, p::PolyElem{<:RingElement};
                    cached::Bool = true,
                    parent::PolyRing = _make_parent(g, p, cached))
   return _map(g, p, parent)
end

function _map(g, p::PolyElem, Rx)
   R = base_ring(Rx)
   new_coefficients = elem_type(R)[let c = coeff(p, i)
                                     iszero(c) ? zero(R) : R(g(c))
                                   end for i in 0:degree(p)]
   return Rx(new_coefficients)
end

###############################################################################
#
#   Special functions
#
###############################################################################

function chebyshev_t_pair(n::Int, x::PolyElem)
   if n == 0
      return one(parent(x)), x
   elseif n == 1
      return x, one(parent(x))
   elseif n < 0
      a, b = chebyshev_t_pair(1 - n, x)
      return b, a
   elseif iseven(n)
      a, b = chebyshev_t_pair(n >> 1, x)
      return 2*(a*a) - 1, 2*(a*b) - x
   else
      a, b = chebyshev_t_pair((n >> 1) + 1, x)
      return 2*(a*b) - x, 2*(b*b) - 1
   end
end

@doc Markdown.doc"""
    chebyshev_t(n::Int, x::PolyElem)

Return the Chebyshev polynomial of the first kind $T_n(x)$, defined by
$T_n(x) = \cos(n \cos^{-1}(x))$.
"""
function chebyshev_t(n::Int, x::PolyElem)
   if n == 0
      return one(parent(x))
   elseif n == 1
      return x
   elseif n < 0
      return chebyshev_t(-n, x)
   elseif iseven(n)
      a = chebyshev_t(n >> 1, x)
      return 2*(a*a) - 1
   else
      a, b = chebyshev_t_pair((n >> 1) + 1, x)
      return 2*(a*b) - x
   end
end

function chebyshev_u_pair(n::Int, x::PolyElem)
   if n == 0
      return one(parent(x)), zero(parent(x))
   elseif n == 1
      return 2*x, one(parent(x))
   elseif n == -1
      return zero(parent(x)), -one(parent(x))
   elseif n < -1
      a, b = chebyshev_u_pair(-1 - n, x)
      return -b, -a
   elseif iseven(n)
      a, b = chebyshev_u_pair(n >> 1, x)
      return (a+b)*(a - b), 2*b*(a-x*b)
   else
      a, b = chebyshev_u_pair(n >> 1, x)
      return 2*a*(x*a - b), (a + b)*(a - b)
   end
end

@doc Markdown.doc"""
    chebyshev_u(n::Int, x::PolyElem)

Return the Chebyshev polynomial of the first kind $U_n(x)$, defined by
$(n+1) U_n(x) = T'_{n+1}(x)$.
"""
function chebyshev_u(n::Int, x::PolyElem)
   if n == 0
      return one(parent(x))
   elseif n == 1
      return 2*x
   elseif n == -1
      return zero(parent(x))
   elseif n < -1
      return -chebyshev_u(-2 - n, x)
   elseif iseven(n)
      a, b = chebyshev_u_pair(n >> 1, x)
      return (a + b)*(a - b)
   else
      a, b = chebyshev_u_pair(n >> 1, x)
      return 2*a*(x*a - b)
   end
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addmul!(z::PolyElem{T}, x::PolyElem{T}, y::PolyElem{T}, c::PolyElem{T}) where T <: RingElement
   c = mul!(c, x, y)
   z = addeq!(z, c)
   return z
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::PolyRing, dr::UnitRange{Int}, _) = elem_type(S)

RandomExtensions.maketype(S::PolyRing, deg::Int, _) = elem_type(S)

function RandomExtensions.make(S::PolyRing, deg_range::UnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, deg_range, vs[1]) # forward to default Make constructor
   else
      make(S, deg_range, make(R, vs...))
   end
end

function RandomExtensions.make(S::PolyRing, deg::Int, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, deg, vs[1]) # forward to default Make constructor
   else
      make(S, deg, make(R, vs...))
   end
end

# define rand for make(S, deg_range, v)
function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make3{<:RingElement,<:PolyRing,UnitRange{Int}}})
   S, deg_range, v = sp[][1:end]
   R = base_ring(S)
   f = S()
   x = gen(S)
   # degree -1 is zero polynomial
   deg = rand(rng, deg_range)
   if deg == -1
      return f
   end
   for i = 0:deg - 1
      f += rand(rng, v)*x^i
   end
   # ensure leading coefficient is nonzero
   c = R()
   while iszero(c)
      c = rand(rng, v)
   end
   f += c*x^deg
   return f
end

# define rand for make(S, deg, v)
function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make3{<:RingElement,<:PolyRing,Int}})
   S, deg, v = sp[][1:end]
   R = base_ring(S)
   f = S()
   x = gen(S)
   # degree -1 is zero polynomial
   if deg == -1
      return f
   end
   for i = 0:deg - 1
      f += rand(rng, v)*x^i
   end
   # ensure leading coefficient is nonzero
   c = R()
   while iszero(c)
      c = rand(rng, v)
   end
   f += c*x^deg
   return f
end

rand(rng::AbstractRNG, S::PolyRing, deg_range::UnitRange{Int}, v...) =
   rand(rng, make(S, deg_range, v...))

rand(rng::AbstractRNG, S::PolyRing, deg::Int, v...) =
   rand(rng, make(S, deg, v...))

rand(S::PolyRing, degs, v...) = rand(Random.GLOBAL_RNG, S, degs, v...)

###############################################################################
#
#   Polynomial substitution
#
###############################################################################

@doc Markdown.doc"""
    subst(f::PolyElem{T}, a::Any) where T <: RingElement

Evaluate the polynomial $f$ at $a$. Note that $a$ can be anything, whether
a ring element or not.
"""
function subst(f::PolyElem{T}, a::U) where {T <: RingElement, U}
   S = parent(a)
   n = degree(f)
   R = base_ring(f)
   if n < 0
      return zero(S) + zero(R)
   elseif n == 0
      return coeff(f, 0)*S(1)
   elseif n == 1
      return coeff(f, 0)*S(1) + coeff(f, 1)*a
   end
   d1 = isqrt(n)
   d = div(n, d1)

   if (U <: Integer && U != BigInt) ||
      (U <: Rational && U != Rational{BigInt})
      c = zero(R)*zero(U)
      V = typeof(c)
      if U != V
         A = powers(map(parent(c), a), d)
      else
         A = powers(a, d)
      end
   else
      A = powers(a, d)
   end

   s = coeff(f, d1*d)*A[1]
   for j = 1:min(n - d1*d, d - 1)
      c = coeff(f, d1*d + j)
      if !iszero(c)
         s += c*A[j + 1]
      end
   end
   for i = 1:d1
      s *= A[d + 1]
      s += coeff(f, (d1 - i)*d)*A[1]
      for j = 1:min(n - (d1 - i)*d, d - 1)
         c = coeff(f, (d1 - i)*d + j)
         if !iszero(c)
            s += c*A[j + 1]
         end
      end
   end
   return s
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

@doc Markdown.doc"""
    PolynomialRing(R::Ring, s::Union{String, Char, Symbol}; cached::Bool = true)

Given a base ring `R` and string `s` specifying how the generator (variable)
should be printed, return a tuple `S, x` representing the new polynomial
ring $S = R[x]$ and the generator $x$ of the ring. By default the parent
object `S` will depend only on `R` and `x` and will be cached. Setting the
optional argument `cached` to `false` will prevent the parent object `S` from
being cached.
"""
PolynomialRing(R::Ring, s::Union{AbstractString, Char, Symbol}; cached::Bool = true)

function PolynomialRing(R::Ring, s::Symbol; cached::Bool = true)
   return Generic.PolynomialRing(R, s; cached=cached)
end

function PolynomialRing(R::Ring, s::AbstractString; cached::Bool = true)
   return PolynomialRing(R, Symbol(s); cached=cached)
end

function PolynomialRing(R::Ring, s::Char; cached::Bool = true)
   return PolynomialRing(R, Symbol(s); cached=cached)
end

function PolyRing(R::Ring)
   T = elem_type(R)
   return Generic.PolyRing{T}(R, :x, false)
end
