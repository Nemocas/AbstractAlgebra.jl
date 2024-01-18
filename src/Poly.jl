###############################################################################
#
#   Poly.jl : Univariate polynomials
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring_type(::Type{PolyRing{T}}) where {T} = parent_type(T)

base_ring(R::PolyRing{T}) where T <: RingElement = R.base_ring::parent_type(T)

coefficient_ring(R::PolyRing) = base_ring(R)

parent(a::PolynomialElem) = a.parent

dense_poly_type(::Type{T}) where T<:RingElement = Generic.Poly{T}

function is_domain_type(::Type{T}) where {S <: RingElement, T <: PolyRingElem{S}}
   return is_domain_type(S)
end

function is_exact_type(a::Type{T}) where {S <: RingElement, T <: PolyRingElem{S}}
   return is_exact_type(S)
end

@doc raw"""
    var(a::PolyRing)

Return the internal name of the generator of the polynomial ring. Note that
this is returned as a `Symbol` not a `String`.
"""
var(a::PolyRing) = a.S

@doc raw"""
    symbols(a::PolyRing)

Return an array of the variable names for the polynomial ring. Note that
this is returned as an array of `Symbol` not `String`.
"""
symbols(a::PolyRing) = [a.S]

@doc raw"""
    number_of_variables(a::PolyRing)

Return the number of variables of the polynomial ring, which is 1.
"""
number_of_variables(a::PolyRing) = 1

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

function Base.hash(a::PolyRingElem, h::UInt)
   b = 0x53dd43cd511044d1%UInt
   for i in 0:length(a) - 1
      b = xor(b, xor(hash(coeff(a, i), h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

@doc raw"""
    length(a::PolynomialElem)

Return the length of the polynomial. The length of a univariate polynomial is
defined to be the number of coefficients in its dense representation, including
zero coefficients. Thus naturally the zero polynomial has length zero and
additionally for nonzero polynomials the length is one more than the degree.
(Note that the leading coefficient will always be nonzero.)
"""
length(a::PolynomialElem) = a.length

@doc raw"""
    degree(a::PolynomialElem)

Return the degree of the given polynomial. This is defined to be one less
than the length, even for constant polynomials.
"""
degree(a::PolynomialElem) = length(a) - 1


@doc raw"""
    is_constant(a::PolynomialElem)

Return `true` if `a` is a degree zero polynomial or the zero polynomial, i.e.
a constant polynomial.
"""
function is_constant(a::PolynomialElem)
   return length(a) <= 1
end

@doc raw"""
    modulus(a::PolyRingElem{T}) where {T <: ResElem}

Return the modulus of the coefficients of the given polynomial.
"""
modulus(a::PolyRingElem{T}) where {T <: ResElem} = modulus(base_ring(a))

@doc raw"""
    leading_coefficient(a::PolynomialElem)

Return the leading coefficient of the given polynomial. This will be the
nonzero coefficient of the term with highest degree unless the polynomial
in the zero polynomial, in which case a zero coefficient is returned.
"""
function leading_coefficient(a::PolynomialElem)
   return length(a) == 0 ? zero(base_ring(a)) : coeff(a, length(a) - 1)
end

@doc raw"""
    trailing_coefficient(a::PolynomialElem)

Return the trailing coefficient of the given polynomial. This will be the
nonzero coefficient of the term with lowest degree unless the polynomial
is the zero polynomial, in which case a zero coefficient is returned.
"""
function trailing_coefficient(a::PolynomialElem)
   if iszero(a)
      return zero(base_ring(a))
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

@doc raw"""
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

@doc raw"""
    tail(a::PolynomialElem)

Return the tail of the given polynomial, i.e. the polynomial without its
leading term (if any).
"""
function tail(a::PolynomialElem)
   return iszero(a) ? zero(parent(a)) : truncate(a, length(a) - 1)
end

@doc raw"""
    set_coefficient!(c::PolynomialElem{T}, n::Int, a::T) where T <: RingElement
    set_coefficient!(c::PolynomialElem{T}, n::Int, a::U) where {T <: RingElement, U <: Integer}

Set the coefficient of degree $n$ to $a$.
"""
function set_coefficient!(c::PolynomialElem{T}, n::Int, a::T) where T <: RingElement
   return setcoeff!(c, n, a) # merely acts as generic fallback
end

function set_coefficient!(c::PolynomialElem{T}, n::Int, a::U) where {T <: RingElement, U <: Integer}
   return setcoeff!(c, n, base_ring(c)(a)) # merely acts as generic fallback
end

function set_coefficient!(c::PolynomialElem{T}, n::Int, a::T) where T <: Integer
   return setcoeff!(c, n, a) # merely acts as generic fallback
end

@doc raw"""
    zero(R::PolyRing)

Return the zero polynomial in the given polynomial ring.
"""
zero(R::PolyRing) = R(zero(base_ring(R)))

one(R::PolyRing) = R(one(base_ring(R)))

@doc raw"""
    gen(R::PolyRing)

Return the generator of the given polynomial ring.
"""
gen(R::PolyRing) = R([zero(base_ring(R)), one(base_ring(R))])

@doc raw"""
    gens(R::PolyRing)

Return an array containing the generator of the given polynomial ring.
"""
gens(R::PolyRing) = [gen(R)]

iszero(a::PolynomialElem) = length(a) == 0

isone(a::PolynomialElem) = length(a) == 1 && isone(coeff(a, 0))

@doc raw"""
    is_gen(a::PolynomialElem)

Return `true` if the given polynomial is the constant generator of its
polynomial ring, otherwise return `false`.
"""
function is_gen(a::PolynomialElem)
    return length(a) <= 2 && isone(coeff(a, 1)) && iszero(coeff(a, 0))
end

@doc raw"""
    is_monic(a::PolynomialElem)

Return `true` if the given polynomial is monic, i.e. has leading coefficient
equal to one, otherwise return `false`.
"""
function is_monic(a::PolynomialElem)
    return isone(leading_coefficient(a))
end

function is_unit(a::PolynomialElem)
   if length(a) <= 1
      return is_unit(coeff(a, 0))
   elseif is_domain_type(elem_type(coefficient_ring(a)))
      return false
   elseif !is_unit(coeff(a, 0)) || is_unit(coeff(a, length(a) - 1))
      return false
   else
      throw(NotImplementedError(:is_unit, a))
   end
end

is_zero_divisor(a::PolynomialElem) = is_zero_divisor(content(a))

function is_zero_divisor_with_annihilator(a::PolyRingElem{T}) where T <: RingElement
   f, b = is_zero_divisor_with_annihilator(content(a))
   return f, parent(a)(b)
end

###############################################################################
#
#  Monomial and term
#
###############################################################################

@doc raw"""
    is_term(a::PolynomialElem)

Return `true` if the given polynomial has one term.
"""
function is_term(a::PolynomialElem)
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

is_term_recursive(a::T) where T <: RingElement = true

@doc raw"""
    is_term_recursive(a::PolynomialElem)

Return `true` if the given polynomial has one term. This function is
recursive, with all scalar types returning true.
"""
function is_term_recursive(a::PolynomialElem)
   if !is_term_recursive(leading_coefficient(a))
      return false
   end
   for i = 1:length(a) - 1
      if !iszero(coeff(a, i - 1))
         return false
      end
   end
   return true
end

@doc raw"""
    is_monomial(a::PolynomialElem)

Return `true` if the given polynomial is a monomial.
"""
function is_monomial(a::PolynomialElem)
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

is_monomial_recursive(a::T) where T <: RingElement = isone(a)

@doc raw"""
    is_monomial_recursive(a::PolynomialElem)

Return `true` if the given polynomial is a monomial. This function is
recursive, with all scalar types returning true.
"""
function is_monomial_recursive(a::PolynomialElem)
   if !is_monomial_recursive(leading_coefficient(a))
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

function similar(x::PolyRingElem, R::Ring, s::VarName=var(parent(x)); cached::Bool=true)
   TT = elem_type(R)
   V = Vector{TT}(undef, 0)
   p = Generic.Poly{TT}(V)
   # Default similar is supposed to return a polynomial
   if base_ring(x) === R && Symbol(s) == var(parent(x)) && x isa Generic.Poly{TT}
      # steal parent in case it is not cached
      p.parent = parent(x)
   else
      p.parent = Generic.PolyRing{TT}(R, Symbol(s), cached)
   end
   p = set_length!(p, 0)
   return p
end

similar(x::PolyRingElem, var::VarName=var(parent(x)); cached::Bool=true) =
   similar(x, base_ring(x), Symbol(var); cached)

zero(p::PolyRingElem, R::Ring, var::VarName=var(parent(p)); cached::Bool=true) =
   similar(p, R, var; cached=cached)

zero(p::PolyRingElem, var::VarName=var(parent(p)); cached::Bool=true) =
   similar(p, base_ring(p), var; cached=cached)

###############################################################################
#
#   polynomial constructor
#
###############################################################################

function polynomial(R::Ring, arr::Vector{T}, var::VarName=:x; cached::Bool=true) where T
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

@doc raw"""
    exponent_vectors(a::PolyRingElem)

Return an iterator for the exponent vectors of the given polynomial. The
exponent vectors will have length 1 and may correspond to terms with zero
coefficient but will not give exponents higher than the degree.
"""
function exponent_vectors(a::PolyRingElem)
   return Generic.MPolyExponentVectors(a)
end

struct PolyCoeffs{T <: RingElement}
   f::T
end
 
function coefficients(f::PolyRingElem)
   return PolyCoeffs(f)
end
 
function Base.iterate(PC::PolyCoeffs{<:PolyRingElem}, st::Int = -1)
   st += 1
   if st > degree(PC.f)
       return nothing
   else
       return coeff(PC.f, st), st
   end
end

function Base.iterate(PCR::Iterators.Reverse{<:PolyCoeffs{<:PolyRingElem}},
                                                   st::Int = degree(PCR.itr.f) + 1)
   st -= 1
   if st < 0
      return nothing
   else
      return coeff(PCR.itr.f, st), st
   end
end
 
Base.IteratorEltype(M::PolyRingElem) = Base.HasEltype()

Base.eltype(M::PolyRingElem{T}) where {T} = T

Base.eltype(M::PolyCoeffs) = Base.eltype(M.f)
 
Base.eltype(M::Iterators.Reverse{<:PolyCoeffs}) = Base.eltype(M.itr.f)

Base.eltype(M::Iterators.Take{<:PolyCoeffs}) = Base.eltype(M.xs.f)

Base.eltype(M::Iterators.Take{<:Iterators.Reverse{<:PolyCoeffs}}) = Base.eltype(M.xs.itr.f)

Base.IteratorSize(M::PolyCoeffs{<:PolyRingElem}) = Base.HasLength()

Base.length(M::PolyCoeffs{<:PolyRingElem}) = length(M.f)
 
function Base.lastindex(a::PolyCoeffs{<:PolyRingElem})
   return degree(a.f)
end
 
function Base.getindex(a::PolyCoeffs{<:PolyRingElem}, i::Int)
   return coeff(a.f, i)
end

function Base.getindex(a::Iterators.Reverse{<:PolyCoeffs{<:PolyRingElem}}, i::Int)
   return coeff(a.itr.f, degree(a.itr.f) - i)
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

function expressify(@nospecialize(a::Union{PolynomialElem, NCPolyRingElem}),
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

@enable_all_show_via_expressify Union{PolynomialElem, NCPolyRingElem}

function show(io::IO, p::PolyRing)
   if get(io, :supercompact, false)
      print(io, "Univariate polynomial ring")
   else
      io = pretty(io)
      print(io, "Univariate polynomial ring in ", var(p), " over ")
      print(IOContext(io, :supercompact => true), Lowercase(), base_ring(p))
   end
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

function +(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement
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

function -(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement
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

@doc raw"""
    mul_karatsuba(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement

Return $a \times b$ using the Karatsuba algorithm.
"""
function mul_karatsuba(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement
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
   # necessary for finite characteristic
   r = set_length!(r, normalise(r, length(r)))
   return r
end

function mul_ks(a::PolyRingElem{T}, b::PolyRingElem{T}) where {T <: PolyRingElem}
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
   A1 = Vector{elem_type(base_ring(base_ring(a)))}(undef, m*lena)
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
      A2 = Vector{elem_type(base_ring(base_ring(a)))}(undef, m*lenb)
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

function mul_classical(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement
   lena = length(a)
   lenb = length(b)
   if lena == 0 || lenb == 0
      return parent(a)()
   end
   R = base_ring(a)
   t = R()
   lenz = lena + lenb - 1
   d = Vector{T}(undef, lenz)
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

function use_karamul(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement
   return false
end

function *(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement
   check_parent(a, b)
   # karatsuba recurses into * so check lengths are > 1
   if use_karamul(a, b) && length(a) > 1 && length(b) > 1
      return mul_karatsuba(a, b)
   else
      return mul_classical(a, b)
   end
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::T, b::PolyRingElem{T}) where {T <: RingElem}
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

*(a::PolyRingElem{T}, b::T) where {T <: RingElem} = b*a

*(a::PolynomialElem, b::Union{Integer, Rational, AbstractFloat}) = b*a

###############################################################################
#
#   Powering
#
###############################################################################

function pow_multinomial(a::PolyRingElem{T}, e::Int) where T <: RingElement
   e < 0 && throw(DomainError(e, "exponent must be >= 0"))
   lena = length(a)
   lenz = (lena - 1) * e + 1
   res = Vector{T}(undef, lenz)
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
   return z
end

@doc raw"""
    ^(a::PolyRingElem{T}, b::Int) where T <: RingElement

Return $a^b$. We require $b \geq 0$.
"""
function ^(a::PolyRingElem{T}, b::Int) where T <: RingElement
   b < 0 && throw(DomainError(b, "exponent must be >= 0"))
   # special case powers of x for constructing polynomials efficiently
   R = parent(a)
   if is_gen(a)
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

@doc raw"""
    ==(x::PolyRingElem{T}, y::PolyRingElem{T}) where T <: RingElement

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(x::PolyRingElem{T}, y::PolyRingElem{T}) where T <: RingElement
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

@doc raw"""
    isequal(x::PolyRingElem{T}, y::PolyRingElem{T}) where T <: RingElement

Return `true` if $x == y$ exactly, otherwise return `false`. This function is
useful in cases where the coefficients of the polynomial are inexact, e.g.
power series. Only if the power series are precisely the same, to the same
precision, are they declared equal by this function.
"""
function isequal(x::PolyRingElem{T}, y::PolyRingElem{T}) where T <: RingElement
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

@doc raw"""
    ==(x::PolyRingElem{T}, y::T) where {T <: RingElem}

Return `true` if $x == y$.
"""
==(x::PolyRingElem{T}, y::T) where T <: RingElem = ((length(x) == 0 && iszero(y))
                        || (length(x) == 1 && coeff(x, 0) == y))

@doc raw"""
    ==(x::PolynomialElem, y::Union{Integer, Rational, AbstractFloat})

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::PolynomialElem, y::Union{Integer, Rational, AbstractFloat}) = ((length(x) == 0 && iszero(base_ring(x)(y)))
                        || (length(x) == 1 && coeff(x, 0) == y))

@doc raw"""
    ==(x::T, y::PolyRingElem{T}) where T <: RingElem = y == x

Return `true` if $x = y$.
"""
==(x::T, y::PolyRingElem{T}) where T <: RingElem = y == x

@doc raw"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::PolyRingElem)

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::PolyRingElem) = y == x

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

@doc raw"""
    truncate(a::PolynomialElem, n::Int)

Return $a$ truncated to $n$ terms, i.e. the remainder upon division by $x^n$.
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

@doc raw"""
    mullow(a::PolyRingElem{T}, b::PolyRingElem{T}, n::Int) where T <: RingElement

Return $a\times b$ truncated to $n$ terms.
"""
function mullow(a::PolyRingElem{T}, b::PolyRingElem{T}, n::Int) where T <: RingElement
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
   d = Vector{T}(undef, lenz)
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
   return z
end

# computes the terms of the product from degree deg(a) + deg(b) down to
# deg(a) + deg(b) - n inclusive, setting the remaining terms to zero
function mulhigh_n(a::PolyRingElem{T}, b::PolyRingElem{T}, n::Int) where T <: RingElement
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
function divexact_low(a::PolyRingElem{T}, b::PolyRingElem{T}, n::Int) where T <: RingElement
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
        q = divexact(coeff(a, 0), coeff(b, 0); check=false)
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
function divhigh(a::PolyRingElem{T}, b::PolyRingElem{T}, n0::Int) where T <: RingElement
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

@doc raw"""
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

@doc raw"""
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

@doc raw"""
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

@doc raw"""
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

@doc raw"""
    deflation(p::PolyRingElem)

Return a tuple `(shift, defl)` where `shift` is the exponent of the trailing
term of $p$ and `defl` is the gcd of the distance between the exponents of the
nonzero terms of $p$. If $p = 0$, both `shift` and `defl` will be zero.
"""
function deflation(p::PolyRingElem)
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

@doc raw"""
    inflate(f::PolyRingElem, shift::Int64, n::Int64) -> PolyRingElem

Given a polynomial $f$ in $x$, return $f(x^n)*x^j$, i.e. multiply
all exponents by $n$ and shift $f$ left by $j$.
"""
function inflate(f::PolyRingElem, j::Int64, n::Int64)
    y = parent(f)()
    for i = 0:degree(f)
        y = setcoeff!(y, n*i + j, coeff(f, i))
    end
    return y
end

@doc raw"""
    inflate(f::PolyRingElem, n::Int64) -> PolyRingElem

Given a polynomial $f$ in $x$, return $f(x^n)$, i.e. multiply
all exponents by $n$.
"""
inflate(f::PolyRingElem, n::Int64) = inflate(f, 0, n)

@doc raw"""
    deflate(f::PolyRingElem, shift::Int64, n::Int64) -> PolyRingElem

Given a polynomial $g$ in $x^n$ such that `f = g(x)*x^{shift}`, write $f$ as
a polynomial in $x$, i.e. divide all exponents of $g$ by $n$.
"""
function deflate(f::PolyRingElem, j::Int64, n::Int64)
    y = parent(f)()
    for i = 0:div(degree(f) - j, n)
        y = setcoeff!(y, i, coeff(f, n*i + j))
    end
    return y
end

@doc raw"""
    deflate(f::PolyRingElem, n::Int64) -> PolyRingElem

Given a polynomial $f$ in $x^n$, write it as a polynomial in $x$, i.e. divide
all exponents by $n$.
"""
deflate(f::PolyRingElem, n::Int64) = deflate(f, 0, n)

@doc raw"""
    deflate(x::PolyRingElem) -> PolyRingElem, Int

Deflate the polynomial $f$ maximally, i.e. find the largest $n$ s.th.
$f$ can be deflated by $n$, i.e. $f$ is actually a polynomial in $x^n$.
Return $g, n$ where $g$ is the deflation of $f$.
"""
function deflate(f::PolyRingElem)
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

function mulmod(a::PolyRingElem{T}, b::PolyRingElem{T}, d::PolyRingElem{T}) where T <: RingElement
   check_parent(a, b)
   check_parent(a, d)
   return mod(a*b, d)
end

function powermod(a::PolyRingElem{T}, b::Int, d::PolyRingElem{T}) where T <: RingElement
   check_parent(a, d)
   if b == 0
      z = one(parent(a))
   elseif length(a) == 0
      z = zero(parent(a))
   elseif length(a) == 1
      z = parent(a)(coeff(a, 0)^b)
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

function invmod(a::PolyRingElem{T}, b::PolyRingElem{T}) where {T <: Union{ResElem, FieldElement}}
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

function divexact(f::PolyRingElem{T}, g::PolyRingElem{T}; check::Bool=true) where T <: RingElement
   check_parent(f, g)
   iszero(g) && throw(DivideError())
   if iszero(f)
      return zero(parent(f))
   end
   lenq = length(f) - length(g) + 1
   d = Vector{T}(undef, lenq)
   for i = 1:lenq
      d[i] = zero(base_ring(f))
   end
   x = gen(parent(f))
   leng = length(g)
   while length(f) >= leng
      lenf = length(f)
      q1 = d[lenf - leng + 1] = divexact(coeff(f, lenf - 1), coeff(g, leng - 1); check=check)
      f = f - shift_left(q1*g, lenf - leng)
      if length(f) == lenf # inexact case
         f = set_length!(f, normalise(f, lenf - 1))
      end
   end
   check && length(f) != 0 && throw(ArgumentError("not an exact division"))
   q = parent(f)(d)
   q = set_length!(q, lenq)
   return q
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(a::PolyRingElem{T}, b::T; check::Bool=true) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact(coeff(a, i - 1), b; check=check))
   end
   z = set_length!(z, length(a))
   return z
end

function divexact(a::PolyRingElem, b::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact(coeff(a, i - 1), b; check=check))
   end
   z = set_length!(z, length(a))
   return z
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function mod(f::PolyRingElem{T}, g::PolyRingElem{T}) where T <: RingElement
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

function rem(f::PolyRingElem{T}, g::PolyRingElem{T}) where T <: RingElement
  return mod(f, g)
end

function Base.divrem(f::PolyRingElem{T}, g::PolyRingElem{T}) where T <: RingElement
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

function Base.div(f::PolyRingElem{T}, g::PolyRingElem{T}) where T <: RingElement
   q, r = divrem(f, g)
   return q
end

##############################################################################
#
#  Ad hoc Euclidean division
#
##############################################################################

function Base.div(f::PolyRingElem{T}, g::T) where T <: Union{FieldElem, ResElem, AbstractFloat, Rational}
   return div(f, parent(f)(g))
end

###############################################################################
#
#   Pseudodivision
#
###############################################################################

@doc raw"""
    pseudorem(f::PolyRingElem{T}, g::PolyRingElem{T}) where T <: RingElement

Return the pseudoremainder of $f$ divided by $g$. If $g = 0$ we throw a
`DivideError()`.
"""
function pseudorem(f::PolyRingElem{T}, g::PolyRingElem{T}) where T <: RingElement
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

@doc raw"""
    pseudodivrem(f::PolyRingElem{T}, g::PolyRingElem{T}) where T <: RingElement

Return a tuple $(q, r)$ consisting of the pseudoquotient and pseudoremainder
of $f$ divided by $g$. If $g = 0$ we throw a `DivideError()`.
"""
function pseudodivrem(f::PolyRingElem{T}, g::PolyRingElem{T}) where T <: RingElement
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

function remove(z::PolyRingElem{T}, p::PolyRingElem{T}) where T <: RingElement
 check_parent(z, p)
 !is_exact_type(T) && error("remove requires an exact ring")
 iszero(z) && error("Not yet implemented")
 (is_unit(p) || iszero(p)) && error("Second argument must be a non-zero non-unit")
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

function remove(z::PolyRingElem{T}, p::PolyRingElem{T}) where T <: Union{ResElem, FieldElement}
 check_parent(z, p)
 !is_exact_type(T) && error("remove requires an exact ring")
 iszero(z) && error("Not yet implemented")
 (is_unit(p) || iszero(p)) && error("Second argument must be a non-zero non-unit")
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

function divides(f::PolyRingElem{T}, g::PolyRingElem{T}) where T <: RingElement
  check_parent(f, g)
  !is_exact_type(T) && error("divides requires an exact ring")
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

function divides(z::PolyRingElem{T}, x::T) where T <: RingElement
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

function sqrt_classical_char2(f::PolyRingElem{T}; check::Bool=true) where T <: RingElement
   S = parent(f)
   R = base_ring(f)
   if check && iszero(f)
      return true, S()
   end
   m = length(f)
   if check && iseven(m) # square polys have even degree
      return false, S()
   end
   if check
      for i = 1:2:m # polynomial must have even exponents
         if !iszero(coeff(f, i))
            return false, S()
         end
      end
   end
   lenq = div(m + 1, 2)
   d = Vector{T}(undef, lenq)
   for i = 1:lenq
      c = coeff(f, 2*i - 2)
      if check && !is_square(c)
         return false, S()
      end
      d[i] = sqrt(c; check=false)
   end
   q = S(d)
   q = set_length!(q, lenq)
   return true, q
end

function sqrt_classical(f::PolyRingElem{T}; check::Bool=true) where T <: RingElement
   S = parent(f)
   R = base_ring(f)
   if characteristic(R) == 2
      return sqrt_classical_char2(f; check=check)
   end
   if iszero(f)
      return true, S()
   end
   m = length(f)
   if check && iseven(m) # square polys have even degree
      return false, S()
   end
   if check && !is_square(coeff(f, m - 1))
      return false, S()
   end
   lenq = div(m + 1, 2)
   d = Vector{T}(undef, lenq)
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
         if check
            flag, d[lenq - k] = divides(qc, b)
            if !flag
               return false, S()
            end
         else
            d[lenq - k] = divexact(qc, b; check=check)
         end
      elseif check && !iszero(qc)
         return false, S()
      end
      k += 1
   end
   q = S(d)
   q = set_length!(q, lenq)
   return true, q
end

@doc raw"""
    Base.sqrt(f::PolyRingElem{T}; check::Bool=true) where T <: RingElement

Return the square root of $f$. By default the function checks the input is
square and raises an exception if not. If `check=false` this check is omitted.
"""
function Base.sqrt(f::PolyRingElem{T}; check::Bool=true) where T <: RingElement
   flag, q = sqrt_classical(f; check=check)
   check && !flag && error("Not a square in sqrt")
   return q
end

@doc raw"""
    is_square(f::PolyRingElem{T}) where T <: RingElement

Return `true` if $f$ is a perfect square.
"""
function is_square(f::PolyRingElem{T}) where T <: RingElement
   flag, q = sqrt_classical(f)
   return flag
end

function is_square_with_sqrt(f::PolyRingElem{T}) where T <: RingElement
   return sqrt_classical(f, check=true)
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

function term_gcd(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement
   d = min(degree(a), degree(b))
   x = gen(parent(a))
   return term_gcd(coeff(a, degree(a)), coeff(b, degree(b)))*x^d
end

function term_content(a::PolyRingElem{T}) where T <: RingElement
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

function gcd(a::PolyRingElem{T}, b::PolyRingElem{T}, ignore_content::Bool = false) where T <: RingElement
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
   lead_monomial = is_term_recursive(leading_coefficient(a)) ||
                   is_term_recursive(leading_coefficient(b))
   trail_monomial = is_term_recursive(trailing_coefficient(a)) ||
                    is_term_recursive(trailing_coefficient(b))
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
      if !is_term_recursive(leading_coefficient(b)) &&
         !is_term_recursive(trailing_coefficient(b))
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
            if is_term_recursive(glead)
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

# can throw NotInvertibleError
function gcd_basecase(a::PolyRingElem{T}, b::PolyRingElem{T}) where T
   while !iszero(b)
      a, b = b, mod(a, b)
   end
   return iszero(a) ? zero(parent(a)) : divexact(a, leading_coefficient(a))
end

# can throw NotInvertibleError
function gcd_hgcd(a::PolyRingElem{T}, b::PolyRingElem{T}) where T
   while !iszero(b)
      a, b = b, mod(a, b)
      if iszero(b) || hgcd_prefers_basecase(a, b)
         break
      else
         a, b, _, _, _, _, _ = hgcd_recursive(a, b, false)
      end
   end
   return gcd_basecase(a, b)
end

# can throw NotInvertibleError for T <: ResElem
#
# To get a good gcd for ZZModPolyRingElem/zzModPolyRingElem that can throw NotInvertibleError:
#  1. Ensure that flint is only used for polynomial addition and multiplication
#     and is never used for mod or divrem or divexact. Or, if using flint for
#     division, check the invertibility of the leading coefficient first.
#  2. tune .._prefers_..
#
#Cutoffs are currently dictated by the non-exported functions
#`hgcd_prefers_basecase(a, b)`
#`mat22_mul_prefers_classical(a11, a12, a21, a22, b11, b12, b21, b22)`
function gcd(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: Union{ResElem, FieldElement}
   check_parent(a, b)
   if length(a) < length(b)
      (a, b) = (b, a)
   end
   if iszero(b)
      if iszero(a)
         return a
      else
         return divexact(a, leading_coefficient(a))
      end
   end
   if T <: ResElem
      # since we return a monic gcd, this step is not strictly necessary
      a = divexact(a, content(a))
      b = divexact(b, content(b))
   end
   return gcd_hgcd(a, b)
end

function lcm(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement
   check_parent(a, b)
   g = gcd(a, b)
   iszero(g) && return g
   return a*divexact(b, g)
end

@doc raw"""
    content(a::PolyRingElem)

Return the content of $a$, i.e. the greatest common divisor of its
coefficients.
"""
function content(a::PolyRingElem)
   z = base_ring(a)() # normalise first coefficient
   for i = 1:length(a)
      z = gcd(z, coeff(a, i - 1))
   end
   return z
end

@doc raw"""
    primpart(a::PolyRingElem)

Return the primitive part of $a$, i.e. the polynomial divided by its content.
"""
function primpart(a::PolyRingElem)
   d = content(a)
   if iszero(d)
      return zero(parent(a))
   else
      return divexact(a, d)
   end
end

###############################################################################
#
#   Half GCD
#
###############################################################################

#Cutoffs are currently dictated by the non-exported functions
#`hgcd_prefers_basecase(a, b)`
#`mat22_mul_prefers_classical(a11, a12, a21, a22, b11, b12, b21, b22)`
function hgcd_prefers_basecase(a, b)
   return true
end

function mat22_mul_prefers_classical(a11, a12, a21, a22, b11, b12, b21, b22)
   return true
end

function mat22_mul(a11, a12, a21, a22, b11, b12, b21, b22)
   if mat22_mul_prefers_classical(a11, a12, a21, a22, b11, b12, b21, b22)
      return a11*b11 + a12*b21, a11*b12 + a12*b22,
             a21*b11 + a22*b21, a21*b12 + a22*b22
   else
      R = parent(a11)
      T0 = R()
      T1 = R()
      C0 = R()
      C1 = R()
      C2 = R()
      C3 = R()
      T0 = sub!(T0, a11, a21)
      T1 = sub!(T1, b22, b12)
      C2 = mul!(C2, T0, T1)
      T0 = add!(T0, a21, a22)
      T1 = sub!(T1, b12, b11)
      C3 = mul!(C3, T0, T1)
      T0 = sub!(T0, T0, a11)
      T1 = sub!(T1, b22, T1)
      C1 = mul!(C1, T0, T1)
      T0 = sub!(T0, a12, T0)
      C0 = mul!(C0, T0, b22)
      T0 = mul!(T0, a11, b11)
      C1 = add!(C1, T0, C1)
      C2 = add!(C2, C1, C2)
      C1 = add!(C1, C1, C3)
      C3 = add!(C3, C2, C3)
      C1 = add!(C1, C1, C0)
      T1 = sub!(T1, T1, b21)
      C0 = mul!(C0, a22, T1)
      C2 = sub!(C2, C2, C0)
      C0 = mul!(C0, a12, b21)
      C0 = add!(C0, C0, T0)
      return C0, C1, C2, C3
   end
end


# iterative basecase
# may throw NotInvertibleError
function hgcd_basecase(a::PolyRingElem{T}, b::PolyRingElem{T}) where T
   @assert length(a) > length(b)
   R = parent(a)
   s = 1
   A = a; B = b
   m11 = one(R);  m12 = zero(R)
   m21 = zero(R); m22 = one(R)
   n = cld(degree(a), 2)
   while n <= degree(B)
      q, r = divrem(A, B)
      @assert length(r) < length(B)
      s = -s
      A, B = B, r
      m11, m12 = m11*q + m12, m11
      m21, m22 = m21*q + m22, m21
   end
   return A, B, m11, m12, m21, m22, s
end

# Klaus Thull and Chee K. Yap
# "A Unified Approach to HGCD Algorithms for polynomials and integers"
# may throw NotInvertibleError
function hgcd_recursive(
   a::PolyRingElem{T},
   b::PolyRingElem{T},
   want_matrix::Bool = true
) where T

   @assert length(a) > length(b)
   R = parent(a)

   if degree(a) < 5 || hgcd_prefers_basecase(a, b)
      return hgcd_basecase(a, b)
   end

   m = cld(degree(a), 2)
   mp = degree(a) - m
   if degree(b) < m
      return (a, b, one(R), zero(R), zero(R), one(R), 1)
   end

   a0 = shift_right(a, m)
   b0 = shift_right(b, m)

   A0, B0, R11, R12, R21, R22, Rs = hgcd_recursive(a0, b0)

   # the calculation of (ap,bp) = R^-1(A,B) can be optimized using A0 and B0
   #     ap = (R22*a - R12*b)*Rs
   #     bp = (-R21*a + R11*b)*Rs
   ar = truncate(a, m)    # a = a0*x^m + ar
   br = truncate(b, m)    # b = b0*x^m + br
   ap = shift_left(A0, m) + (R22*ar - R12*br)*Rs
   bp = shift_left(B0, m) + (R11*br - R21*ar)*Rs

   if degree(bp) < m
      return (ap, bp, R11, R12, R21, R22, Rs)
   end

   q, r = divrem(ap, bp)
   @assert length(r) < length(bp)
   c = bp
   d = r

   l = degree(c)
   k = 2*m - l
   @assert l - m < cld(mp, 2)

   c0 = shift_right(c, k)
   d0 = shift_right(d, k)
   @assert degree(c0) == 2*(l - m)
   C0, D0, S11, S12, S21, S22, Ss = hgcd_recursive(c0, d0)

   # (A,B) = Q^-1(a,b) = S^-1(c,d) can be optimized as well
   #     A = (Q22*a - Q12*b)*Qs
   #     B = (-Q21*a + Q11*b)*Qs
   cr = truncate(c, k)    # c = c0*x^k + cr
   dr = truncate(d, k)    # d = d0*x^k + dr
   A = shift_left(C0, k) + (S22*cr - S12*dr)*Ss
   B = shift_left(D0, k) + (S11*dr - S21*cr)*Ss

   Qs = -Rs*Ss
   if want_matrix
      (Q11, Q12, Q21, Q22) = mat22_mul(R11*q + R12, R11, R21*q + R22, R21,
                                       S11, S12, S21, S22)
   else
      Q11 = Q12 = Q21 = Q22 = R()
   end

   return A, B, Q11, Q12, Q21, Q22, Qs
end

@doc raw"""
    hgcd(a::PolyRingElem{T}, b::PolyRingElem{T}) where T

Returns the half-GCD of `a` and `b`, that is, a tuple
`(A, B, m11, m12, m21, m22, s::Int)` such that
   1. m11*m22 - m12*m21 = s = +-1
   2. m11*A + m12*B = a
   3. m21*A + m22*B = b
   4. `A` and `B` are consecutive remainders in the sequence for
      division of `a` by `b` satisfying `deg(A) >= cld(deg(a), 2) > deg(B)
   5. `[m11 m12; m21 m22]` is the product of `[q 1; 1 0]` for the quotients `q`
      generated by such a remainder sequence.

Assumes that the input satsifies `degree(a) > degree(b) >= 0`.

Cutoffs are currently dictated by the non-exported functions
`hgcd_prefers_basecase(a, b)`
`mat22_mul_prefers_classical(a11, a12, a21, a22, b11, b12, b21, b22)`

If the base ring of the polynomial ring is not a field, this function may throw
a NotInvertibleError. Otherwise, the output should be valid.
"""
function hgcd(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement
   check_parent(a, b)
   @assert degree(a) > degree(b) >= 0
   return hgcd_recursive(a, b, true)
end

###############################################################################
#
#   Evaluation/composition
#
###############################################################################

@doc raw"""
    evaluate(a::PolyRingElem, b::T) where T <: RingElement

Evaluate the polynomial expression $a$ at the value $b$ and return the result.
"""
function evaluate(a::PolyRingElem, b::T) where T <: RingElement
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

@doc raw"""
    compose(a::PolyRingElem, b::PolyRingElem)

Compose the polynomial $a$ with the polynomial $b$ and return the result,
i.e. return $a\circ b$.
"""
function compose(a::PolyRingElem, b::PolyRingElem)
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

@doc raw"""
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

@doc raw"""
    integral(x::PolyRingElem{T}) where {T <: Union{ResElem, FieldElement}}

Return the integral of the polynomial $x$.
"""
function integral(x::PolyRingElem{T}) where {T <: Union{ResElem, FieldElement}}
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
function subresultant_lazard(Sd0::PolyRingElem{T}, Sd1::PolyRingElem{T}) where T <: RingElement
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
function subresultant_ducos(A::PolyRingElem{T}, Sd1::PolyRingElem{T}, Se0::PolyRingElem{T}, sd::T) where T <: RingElement
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

@doc raw"""
    resultant_ducos(p::PolyRingElem{T}, q::PolyRingElem{T}) where T <: RingElement

Return the resultant of the $p$ and $q$.
"""
function resultant_ducos(p::PolyRingElem{T}, q::PolyRingElem{T}) where T <: RingElement
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
function resultant_subresultant(p::PolyRingElem{T}, q::PolyRingElem{T}) where T <: RingElement
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

function resultant_lehmer(a::PolyRingElem{T}, b::PolyRingElem{T}) where {T <: Union{ResElem, FieldElement}}
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
   s = one(R)
   while lenB > crossover/2 + 1
      shift = max(lenA - crossover, 0)
      a = shift_right(A, shift)
      b = shift_right(B, shift)
      u1, v1 = one(R), zero(R)
      u2, v2 = zero(R), one(R)
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
         return zero(base_ring(a)), one(parent(A)), one(parent(A))
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

@doc raw"""
    sylvester_matrix(p::PolyRingElem, q::PolyRingElem)

Return the sylvester matrix of the given polynomials.
"""
function sylvester_matrix(p::PolyRingElem{T}, q::PolyRingElem{T}) where T <: RingElement
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

function resultant_sylvester(p::PolyRingElem{T}, q::PolyRingElem{T}) where T <: RingElement
   check_parent(p, q)
   R = base_ring(p)
   if length(p) == 0 || length(q) == 0
      return zero(R)
   end
   return det_df(sylvester_matrix(p, q))
end

@doc raw"""
    resultant(p::PolyRingElem{T}, q::PolyRingElem{T}) where T <: RingElement

Return the resultant of the given polynomials.
"""
function resultant(p::PolyRingElem{T}, q::PolyRingElem{T}) where T <: RingElement
  R = parent(p)
  if !is_exact_type(T)
     return resultant_sylvester(p, q)
  end
  try
     return resultant_ducos(p, q)
  catch
     return resultant_sylvester(p, q)
  end
end

function resultant_euclidean(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: Union{ResElem, FieldElement}
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
   s = one(base_ring(A))
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

function resultant(a::PolyRingElem{T}, b::PolyRingElem{T}) where {T <: Union{ResElem, FieldElement}}
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

@doc raw"""
    discriminant(a::PolyRingElem)

Return the discriminant of the given polynomial.
"""
function discriminant(a::PolyRingElem)
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

@doc raw"""
    resx(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement

Return a tuple $(r, s, t)$ such that $r$ is the resultant of $a$ and $b$ and
such that $r = a\times s + b\times t$.
"""
function resx(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: RingElement
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

function gcdx_basecase(a::PolyRingElem{T}, b::PolyRingElem{T}) where T
   @assert length(a) > 0
   @assert length(b) > 0
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

#Cutoffs are currently dictated by the non-exported functions
#`hgcd_prefers_basecase(a, b)`
#`mat22_mul_prefers_classical(a11, a12, a21, a22, b11, b12, b21, b22)`
function gcdx(a::PolyRingElem{T}, b::PolyRingElem{T}) where {T <: Union{ResElem, FieldElement}}
   check_parent(a, b)
   R = parent(a)
   !is_exact_type(T) && error("gcdx requires exact Bezout domain")
   if length(a) == 0
      if length(b) == 0
         return zero(R), zero(R), zero(R)
      else
         d = leading_coefficient(b)
         return divexact(b, d), zero(R), divexact(one(R), d)
      end
   end
   if length(b) == 0
      d = leading_coefficient(a)
      return divexact(a, d), divexact(one(R), d), zero(R)
   end
   if hgcd_prefers_basecase(a, b)
      return gcdx_basecase(a, b)
   else
      g, s = gcdinv_hgcd(a, b)
      return g, s, divexact(g - s*a, b)
   end
end

function gcdinv_basecase(a::PolyRingElem{T}, b::PolyRingElem{T}) where T
   @assert length(a) > 0
   @assert length(b) > 0
   R = parent(a)
   if length(a) < length(b)
      a, b = b, a
      u1, u2 = zero(R), one(R)
   else
      u1, u2 = one(R), zero(R)
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

function gcdinv_hgcd(a::PolyRingElem{T}, b::PolyRingElem{T}) where T
   @assert length(a) > 0
   @assert length(b) > 0
   R = parent(a)
   m21, m22 = zero(R), one(R)
   ms = 1
   if length(a) < length(b)
      a, b = b, a
      m21, m22 = m22, m21
      ms = -ms
   end
   while !iszero(b)
      q, r = divrem(a, b)
      a, b = b, r
      m21, m22 = m21*q + m22, m21
      ms = -ms
      if iszero(b)
         break
      elseif !hgcd_prefers_basecase(a, b)
         a, b, n11, n12, n21, n22, ns = hgcd_recursive(a, b, true)
         m21, m22 = m21*n11 + m22*n21, m21*n12 + m22*n22
         ms *= ns
      end
   end
   l = leading_coefficient(a)
   g = divexact(a, l)
   s = divexact(m22, ms*l)
   return g, s
end

#Cutoffs are currently dictated by the non-exported functions
#`hgcd_prefers_basecase(a, b)`
#`mat22_mul_prefers_classical(a11, a12, a21, a22, b11, b12, b21, b22)`
function gcdinv(a::PolyRingElem{T}, b::PolyRingElem{T}) where T <: Union{ResElem, FieldElement}
   check_parent(a, b)
   R = parent(a)
   if length(a) == 0
      if length(b) == 0
         return zero(R), zero(R)
      else
         d = leading_coefficient(b)
         return divexact(b, d), zero(R)
      end
   end
   if length(b) == 0
      d = leading_coefficient(a)
      return divexact(a, d), R(inv(d))
   end
   if hgcd_prefers_basecase(a, b)
      return gcdinv_basecase(a, b)
   else
      return gcdinv_hgcd(a, b)
   end
end

################################################################################
#
#  Power sums
#
################################################################################

function polynomial_to_power_sums(f::PolyRingElem{T}, n::Int=degree(f)) where T <: FieldElement
    degree(f) < 1 && error("Polynomial has no roots")
    !is_monic(f) && error("Requires monic polynomial")
    iszero(constant_coefficient(f)) && error("Requires nonzero constant coefficient")
    n < 0 && throw(DomainError(n, "number of terms must be non-negative"))
    d = degree(f)
    R = base_ring(f)
    # Beware: converting to power series and derivative do not commute
    dfc = collect(Iterators.take(Iterators.reverse(coefficients(derivative(f))), n + 1))
    A = abs_series(R, dfc, length(dfc), n + 1; cached=false)
    S = parent(A)
    fc = collect(Iterators.take(Iterators.reverse(coefficients(f)), n + 1))
    B = S(fc, length(fc), n + 1)
    L = A*inv(B)
    s = T[coeff(L, i) for i = 1:n]
    return s
end

@doc raw"""
    polynomial_to_power_sums(f::PolyRingElem{T}, n::Int=degree(f)) where T <: RingElement -> Vector{T}

Uses Newton (or Newton-Girard) formulas to compute the first $n$
sums of powers of the roots of $f$ from the coefficients of $f$, starting
with the sum of (first powers of) the roots. The input polynomial must be
monic, at least degree $1$ and have nonzero constant coefficient.
"""
function polynomial_to_power_sums(f::PolyRingElem{T}, n::Int=degree(f)) where T <: RingElement
    # plain vanilla recursion
    degree(f) < 1 && error("Polynomial has no roots")
    !is_monic(f) && error("Requires monic polynomial")
    iszero(constant_coefficient(f)) && error("Requires nonzero constant coefficient")
    n < 0 && throw(DomainError(n, "number of terms must be non-negative"))
    d = degree(f)
    R = base_ring(f)
    if n == 0
       return elem_type(R)[]
    end
    if n == 1
       return [-coeff(f, d - 1)]
    end
    E = T[(-1)^i*coeff(f, d - i) for i = 0:min(d, n)] # elementary symm. polys
    while length(E) <= n
        push!(E, R())
    end
    P = T[]
    push!(P, E[1 + 1])
    for k = 2:n
        push!(P, (-1)^(k - 1)*k*E[k + 1] +
              sum((-1)^(k - 1 + i)*E[k - i + 1]*P[i] for i = 1:k - 1))
    end
    return P
end

@doc raw"""
    power_sums_to_polynomial(P::Vector{T};
                     parent::PolyRing{T}=PolyRing(parent(P[1])) where T <: RingElement -> PolyRingElem{T}

Uses the Newton (or Newton-Girard) identities to obtain the polynomial
with given sums of powers of roots. The list must be nonempty and contain
`degree(f)` entries where $f$ is the polynomial to be recovered. The list
must start with the sum of first powers of the roots.
"""
function power_sums_to_polynomial(P::Vector{T}; 
                           parent::PolyRing{T}=PolyRing(parent(P[1]))) where T <: RingElement
   return power_sums_to_polynomial(P, parent)
end

function power_sums_to_polynomial(P::Vector{T}, Rx::PolyRing{T}) where T <: FieldElement
    d = length(P)
    R = base_ring(Rx)
    s = rel_series(R, P, d, d, 0)
    r = -integral(s)
    r1 = exp(r)
    @assert iszero(valuation(r1))
    return Rx(T[polcoeff(r1, d - i) for i = 0:d])
end

function power_sums_to_polynomial(P::Vector{T}, Rx::PolyRing{T}) where T <: RingElement
    E = T[one(parent(P[1]))]
    R = parent(P[1])
    last_non_zero = 0
    for k = 1:length(P)
        push!(E, divexact(sum((-1)^(i - 1)*E[k - i + 1]*P[i] for i = 1:k), R(k)))
        if E[end] != 0
            last_non_zero = k
        end
    end
    E = E[1:last_non_zero + 1]
    d = length(E) # the length of the resulting polynomial
    for i = 1:div(d, 2)
        E[i], E[d - i + 1] = (-1)^(d - i)*E[d - i + 1], (-1)^(i - 1)*E[i]
    end
    if isodd(d)
        E[div(d + 1, 2)] *= (-1)^div(d, 2)
    end
    return Rx(E)
end

###############################################################################
#
#   Newton representation
#
###############################################################################

@doc raw"""
    monomial_to_newton!(P::Vector{T}, roots::Vector{T}) where T <: RingElement

Converts a polynomial $p$, given as an array of coefficients, in-place
from its coefficients given in the standard monomial basis to the Newton
basis for the roots $r_0, r_1, \ldots, r_{n-2}$. In other words, this
determines output coefficients $c_i$ such that
$$c_0 + c_1(x-r_0) + c_2(x-r_0)(x-r_1) + \ldots + c_{n-1}(x-r_0)(x-r_1)\cdots(x-r_{n-2})$$
is equal to the input polynomial.
"""
function monomial_to_newton!(P::Vector{T}, roots::Vector{T}) where T <: RingElement
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

@doc raw"""
    newton_to_monomial!(P::Vector{T}, roots::Vector{T}) where T <: RingElement

Converts a polynomial $p$, given as an array of coefficients, in-place
from its coefficients given in the Newton basis for the roots
$r_0, r_1, \ldots, r_{n-2}$ to the standard monomial basis. In other words,
this evaluates
$$c_0 + c_1(x-r_0) + c_2(x-r_0)(x-r_1) + \ldots + c_{n-1}(x-r_0)(x-r_1)\cdots(x-r_{n-2})$$
where $c_i$ are the input coefficients given by $p$.
"""
function newton_to_monomial!(P::Vector{T}, roots::Vector{T}) where T <: RingElement
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

@doc raw"""
    interpolate(S::PolyRing, x::Vector{T}, y::Vector{T}) where T <: RingElement

Given two arrays of values $xs$ and $ys$ of the same length $n$, find
the polynomial $f$ in the polynomial ring $R$ of length at most $n$ such that
$f$ has the value $ys$ at the points $xs$. The values in the arrays $xs$ and
$ys$ must belong to the base ring of the polynomial ring $R$. If no such
polynomial exists, an exception is raised.
"""
function interpolate(S::PolyRing, x::Vector{T}, y::Vector{T}) where T <: RingElement
   length(x) != length(y) && error("Array lengths don't match in interpolate")
   !is_domain_type(T) && error("interpolation requires a domain type")
   n = length(x)
   if n == 0
      return S()
   elseif n == 1
      return S(y[1])
   end
   R = base_ring(S)
   parent(y[1]) != R && error("Polynomial ring does not match inputs")
   P = Vector{T}(undef, n)
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
   return r
end

function interpolate(S::PolyRing, x::Vector{T}, y::Vector{T}) where {T <: ResElem}
   length(x) != length(y) && error("Array lengths don't match in interpolate")
   n = length(x)
   if n == 0
      return S()
   elseif n == 1
      return S(y[1])
   end
   R = base_ring(S)
   parent(y[1]) != R && error("Polynomial ring does not match inputs")
   P = Vector{T}(undef, n)
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
   return r
end

################################################################################
#
#  Change base ring
#
################################################################################

function _change_poly_ring(R, Rx, cached)
   P, _ = polynomial_ring(R, var(Rx), cached = cached)
   return P
end

@doc raw"""
    change_base_ring(R::Ring, p::PolyRingElem{<: RingElement}; parent::PolyRing)

Return the polynomial obtained by coercing the non-zero coefficients of `p`
into `R`.

If the optional `parent` keyword is provided, the polynomial will be an
element of `parent`. The caching of the parent object can be controlled
via the `cached` keyword argument.
"""
function change_base_ring(R::Ring, p::PolyRingElem{T}; cached::Bool = true, parent::PolyRing = _change_poly_ring(R, parent(p), cached)) where T <: RingElement
   return _map(R, p, parent)
end

@doc raw"""
    change_coefficient_ring(R::Ring, p::PolyRingElem{<: RingElement}; parent::PolyRing)

Return the polynomial obtained by coercing the non-zero coefficients of `p`
into `R`.

If the optional `parent` keyword is provided, the polynomial will be an
element of `parent`. The caching of the parent object can be controlled
via the `cached` keyword argument.
"""
function change_coefficient_ring(R::Ring, p::PolyRingElem{T}; cached::Bool = true, parent::PolyRing = _change_poly_ring(R, parent(p), cached)) where T <: RingElement
  return change_base_ring(R, p; cached = cached, parent = parent)
end

################################################################################
#
#  Map
#
################################################################################

_make_parent(g::T, p::PolyRingElem, cached::Bool) where T =
   _change_poly_ring(parent(g(zero(base_ring(p)))),
                     parent(p), cached)

@doc raw"""
    map_coefficients(f, p::PolyRingElem{<: RingElement}; cached::Bool=true, parent::PolyRing)

Transform the polynomial `p` by applying `f` on each non-zero coefficient.

If the optional `parent` keyword is provided, the polynomial will be an
element of `parent`. The caching of the parent object can be controlled
via the `cached` keyword argument.
"""
function map_coefficients(g::T, p::PolyRingElem{<:RingElement};
                    cached::Bool = true,
                    parent::PolyRing = _make_parent(g, p, cached)) where T
   return _map(g, p, parent)
end

function _map(g::T, p::PolyRingElem, Rx) where T
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

function chebyshev_t_pair(n::Int, x::PolyRingElem)
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

@doc raw"""
    chebyshev_t(n::Int, x::PolyRingElem)

Return the Chebyshev polynomial of the first kind $T_n(x)$, defined by
$T_n(x) = \cos(n \cos^{-1}(x))$.
"""
function chebyshev_t(n::Int, x::PolyRingElem)
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

function chebyshev_u_pair(n::Int, x::PolyRingElem)
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

@doc raw"""
    chebyshev_u(n::Int, x::PolyRingElem)

Return the Chebyshev polynomial of the first kind $U_n(x)$, defined by
$(n+1) U_n(x) = T'_{n+1}(x)$.
"""
function chebyshev_u(n::Int, x::PolyRingElem)
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

function addmul!(z::PolyRingElem{T}, x::PolyRingElem{T}, y::PolyRingElem{T}, c::PolyRingElem{T}) where T <: RingElement
   c = mul!(c, x, y)
   z = addeq!(z, c)
   return z
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::PolyRing, dr::AbstractUnitRange{Int}, _) = elem_type(S)

RandomExtensions.maketype(S::PolyRing, deg::Int, _) = elem_type(S)

function RandomExtensions.make(S::PolyRing, deg_range::AbstractUnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, deg_range, vs[1]) # forward to default Make constructor
   else
      Make(S, deg_range, make(R, vs...))
   end
end

function RandomExtensions.make(S::PolyRing, deg::Int, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, deg, vs[1]) # forward to default Make constructor
   else
      Make(S, deg, make(R, vs...))
   end
end

# define rand for make(S, deg_range, v)
function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make3{<:RingElement, <:PolyRing, <:AbstractUnitRange{Int}}})
   S, deg_range, v = sp[][1:end]
   R = base_ring(S)
   len = 1 + rand(rng, deg_range)
   if len <= 0
      return zero(S)
   end
   c = elem_type(R)[rand(rng, v) for i in 1:len]
   # ensure leading coefficient is nonzero
   while iszero(c[len])
      c[len] = rand(rng, v)
   end
   return S(c)
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

rand(rng::AbstractRNG, S::PolyRing, deg_range::AbstractUnitRange{Int}, v...) =
   rand(rng, make(S, deg_range, v...))

rand(rng::AbstractRNG, S::PolyRing, deg::Int, v...) =
   rand(rng, make(S, deg, v...))

rand(S::PolyRing, degs, v...) = rand(Random.GLOBAL_RNG, S, degs, v...)

###############################################################################
#
#   Polynomial substitution
#
###############################################################################

@doc raw"""
    subst(f::PolyRingElem{T}, a::Any) where T <: RingElement

Evaluate the polynomial $f$ at $a$. Note that $a$ can be anything, whether
a ring element or not.
"""
function subst(f::PolyRingElem{T}, a::U) where {T <: RingElement, U}
   S = parent(a)
   n = degree(f)
   R = base_ring(f)
   if n < 0
      return zero(S) + zero(R)
   elseif n == 0
      return coeff(f, 0)*one(S)
   elseif n == 1
      return coeff(f, 0)*one(S) + coeff(f, 1)*a
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
