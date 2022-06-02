###############################################################################
#
#   RationalFunctionField.jl : Generic rational function fields
#
###############################################################################

export RationalFunctionField, norm

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{Rat{T, U}}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}} = RationalFunctionField{T, U}

elem_type(::Type{RationalFunctionField{T, U}}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}} = Rat{T, U}

base_ring(a::RationalFunctionField{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}} = a.base_ring::parent_type(T)

base_ring(a::Rat) = base_ring(parent(a))

parent(a::Rat) = a.parent

data(x::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}} = x.d::Union{Frac{U}}

function fraction_field(a::RationalFunctionField{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   return a.fraction_field::Union{FracField{U}}
end

function is_domain_type(::Type{S}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}, S <: Rat{T, U}}
   return true
end

function is_exact_type(a::Type{T}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}, S <: Rat{T, U}}
   return is_exact_type(T)
end

function characteristic(R::RationalFunctionField)
   return characteristic(base_ring(R))
end

function check_parent(a::Rat{T, U}, b::Rat{T, U}, throw::Bool = true) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   fl = parent(a) != parent(b)
   fl && throw && error("Incompatible rings in rationa function field operation")
   return !fl
end

###############################################################################
#
#   Constructors
#
###############################################################################

function //(x::Rat{T, U}, y::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   check_parent(x, y)
   R = parent(x)
   return R(divexact(data(x), data(y)))
end

function //(x::T, y::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   R = parent(y)
   parent(x) != base_ring(y) && error("Incompatible elements")
   return R(divexact(x, data(y)))
end

function //(x::Rat{T, U}, y::T) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   R = parent(x)
   base_ring(x) != parent(y) && error("Incompatible elements")
   return R(divexact(data(x), y))
end

function //(x::Rational{BigInt}, y::Rat{Rational{BigInt}, U}) where U <: Union{PolyElem, MPolyElem}
   R = parent(y)
   return R(divexact(x, data(y)))
end

function //(x::Rat{Rational{BigInt}, U}, y::Rational{BigInt}) where U <: Union{PolyElem, MPolyElem}
   R = parent(x)
   return R(divexact(data(x), y))
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::Rat, h::UInt)
   b = 0x2d122a968560a3c0%UInt
   return xor(b, hash(data(a), h))
end

function Base.numerator(a::Rat, canonicalise::Bool=true)
   return numerator(data(a), canonicalise)
end

function Base.denominator(a::Rat, canonicalise::Bool=true)
   return denominator(data(a), canonicalise)
end

zero(R::RationalFunctionField) = R()

one(R::RationalFunctionField) = R(1)

iszero(a::Rat) = iszero(data(a))

isone(a::Rat) = isone(data(a))

is_unit(a::Rat) = is_unit(data(a))

gen(R::RationalFunctionField) = R(gen(base_ring(R.fraction_field)))

gens(R::RationalFunctionField) = [R(g) for g in gens(base_ring(R.fraction_field))]

function deepcopy_internal(a::Rat, dict::IdDict)
   R = parent(a)
   return R(deepcopy_internal(data(a), dict))
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::Rat) = a

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function AbstractAlgebra.expressify(a::Rat; context = nothing)
   d = data(a)
   return expressify(d; context)
end

function show(io::IO, ::MIME"text/plain", a::Rat)
   print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function show(io::IO, a::Rat)
   print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function show(io::IO, a::RationalFunctionField)
   print(IOContext(io, :compact => true), "Rational function field over ", base_ring(a))
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::Rat)
   R = parent(a)
   return R(-data(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   check_parent(a, b)
   R = parent(a)
   return R(data(a) + data(b))
end

function -(a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   check_parent(a, b)
   R = parent(a)
   return R(data(a) - data(b))
end

function *(a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   check_parent(a, b)
   R = parent(a)
   return R(data(a) * data(b))
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::Rat, b::Union{Integer, Rational, AbstractFloat})
   R = parent(a)
   return R(data(a)*b)
end

function *(a::Union{Integer, Rational, AbstractFloat}, b::Rat)
   R = parent(b)
   return R(a*data(b))
end

function *(a::Rat{T, U}, b::T) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}}
   R = parent(a)
   return R(data(a)*b)
end

function *(a::T, b::Rat{T, U}) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}}
   R = parent(b)
   return R(a*data(b))
end

function +(a::Rat, b::Union{Integer, Rational, AbstractFloat})
   R = parent(a)
   return R(data(a) + b)
end

function +(a::Union{Integer, Rational, AbstractFloat}, b::Rat)
   R = parent(b)
   return R(a + data(b))
end

function +(a::Rat{T, U}, b::T) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}}
   R = parent(a)
   return R(data(a) + b)
end

function +(a::T, b::Rat{T, U}) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}}
   R = parent(b)
   return R(a + data(b))
end

function -(a::Rat, b::Union{Integer, Rational, AbstractFloat})
   R = parent(a)
   return R(data(a) - b)
end

function -(a::Union{Integer, Rational, AbstractFloat}, b::Rat)
   R = parent(b)
   return R(a - data(b))
end

function -(a::Rat{T, U}, b::T) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}}
   R = parent(a)
   return R(data(a) - b)
end

function -(a::T, b::Rat{T, U}) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}}
   R = parent(b)
   return R(a - data(b))
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function ==(a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   check_parent(a, b)
   return data(a) == data(b)
end

function isequal(a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   check_parent(a, b)
   return data(a) == data(b)
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function ==(a::Rat, b::Union{Integer, Rational, AbstractFloat})
   R = parent(a)
   return data(a) == b
end

function ==(a::Union{Integer, Rational, AbstractFloat}, b::Rat)
   R = parent(b)
   return a == data(b)
end

function ==(a::Rat{T, U}, b::T) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}}
   R = parent(a)
   return data(a) == b
end

function ==(a::T, b::Rat{T, U}) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}}
   R = parent(b)
   return a == data(b)
end

function ==(a::Rat{T, U}, b::Poly{T}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   return a == parent(a)(b)
end

function ==(a::Poly{T}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   return parent(b)(a) == b
end

###############################################################################
#
#   Inversion
#
###############################################################################

function Base.inv(a::Rat)
   R = parent(a)
   return R(inv(data(a)))
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::Rat{T, U}, b::Rat{T, U}; check::Bool=true) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   check_parent(a, b)
   R = parent(a)
   return R(divexact(data(a), data(b); check=check))
end

function divides(a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   check_parent(a, b)
   R = parent(a)
   d, q = divides(data(a), data(b))
   return d, R(q)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(a::Rat, b::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
   R = parent(a)
   return R(divexact(data(a), b; check=check))
end

function divexact(a::Union{Integer, Rational, AbstractFloat}, b::Rat; check::Bool=true)
   R = parent(b)
   return R(divexact(a, data(b); check=check))
end

function divexact(a::Rat{T, U}, b::T; check::Bool=true) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}}
   R = parent(a)
   return R(divexact(data(a), b; check=check))
end

function divexact(a::T, b::Rat{T, U}; check::Bool=true) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}}
   R = parent(b)
   return R(divexact(a, data(b); check=check))
end

##############################################################################
#
#  Evaluation
#
##############################################################################

function evaluate(f::Rat{T, U}, v::V) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}, V <: RingElement}
    return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

function evaluate(f::Rat{T, U}, v::Vector{V}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}, V <: RingElement}
   return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::Rat, b::Int)
   R = parent(a)
   return R(data(a)^b)
end

##############################################################################
#
#  Derivative
#
##############################################################################

function derivative(f::Rat)
   R = parent(f)
   return R(derivative(data(f)))
end

function derivative(f::Rat, x::MPolyElem)
   R = parent(f)
   return R(derivative(data(f), x))
end

###############################################################################
#
#   Square root
#
###############################################################################

@doc Markdown.doc"""
    is_square(a::Rat)

Return `true` if $a$ is a square.
"""
function is_square(a::Rat)
   return is_square(data(a))
end

@doc Markdown.doc"""
    Base.sqrt(a::Rat; check::Bool=true)

Return the square root of $a$. By default the function will throw an exception
if the input is not square. If `check=false` this test is omitted.
"""
function Base.sqrt(a::Rat; check::Bool=true)
   R = parent(a)
   return R(sqrt(data(a); check=check))
end

###############################################################################
#
#   GCD
#
###############################################################################

@doc Markdown.doc"""
    gcd(a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}

Return a greatest common divisor of $a$ and $b$ if one exists. N.B: we define
the GCD of $a/b$ and $c/d$ to be gcd$(ad, bc)/bd$, reduced to lowest terms.
"""
function gcd(a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   check_parent(a, b)
   R = parent(a)
   return R(gcd(data(a), data(b)))
end

###############################################################################
#
#   Unsafe operators and functions
#
###############################################################################

function zero!(c::Rat)
   c.d = zero!(data(c))
   return c
end

function mul!(c::Rat{T, U}, a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   c.d = mul!(data(c), data(a), data(b))
   return c
end

function addeq!(a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   a.d = addeq!(data(a), data(b))
   return a
end

function add!(c::Rat{T}, a::Rat{T, U}, b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   c.d = add!(data(c), data(a), data(b))
   return c
end

###############################################################################
#
#   Random functions
#
###############################################################################

RandomExtensions.maketype(R::RationalFunctionField, _) = elem_type(R)

function RandomExtensions.make(S::RationalFunctionField, vs...)
   R = base_ring(fraction_field(S))
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      RandomExtensions.Make(S, vs[1]) # forward to default Make constructor
   else
      make(S, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:FieldElement, <:RationalFunctionField}})
   S, v = sp[][1:end]
   R = base_ring(fraction_field(S))
   n = rand(rng, v)
   d = R()
   while iszero(d)
      d = rand(rng, v)
   end
   return S(fraction_field(S)(n, d))
end

rand(rng::AbstractRNG, S::RationalFunctionField, v...) =
   rand(rng, make(S, v...))

rand(S::RationalFunctionField, v...) = rand(GLOBAL_RNG, S, v...)

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{Rat{T, U}}, ::Type{Rat{T, U}}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}} = Rat{T, U}

promote_rule(::Type{Rat{T, U}}, ::Type{Rat{T, U}}) where {T <: FieldElem, U <: Union{PolyElem, MPolyElem}} = Rat{T, U}

function promote_rule(::Type{Rat{T, U}}, ::Type{V}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}, V <: RingElem}
   promote_rule(Frac{U}, V) === Frac{U} ? Rat{T, U} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::RationalFunctionField{T, U})() where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   K = fraction_field(a)
   z = Rat{T, U}(K())
   z.parent = a
   return z
end

function (a::RationalFunctionField{T, U})(b::Frac{U}) where {T <: FieldElement, U <: Union{PolyElem{T}, MPolyElem{T}}}
   K = fraction_field(a)
   parent(b) != K && error("Unable to coerce rational function")
   z = Rat{T, U}(b)
   z.parent = a
   return z::Rat{T, U}
end

function (a::RationalFunctionField{T, U})(n::U, d::U) where {T <: FieldElement, U <: Union{PolyElem{T}, MPolyElem{T}}}
   R = parent(n)
   g = gcd(n, d)
   if !isone(g)
      n = divexact(n, g)
      d = divexact(d, g)
   end
   r = Frac{U}(n, d)
   try
      r.parent = FracDict[R]
   catch
      r.parent = FractionField(R)
   end
   return a(r)
end

function (a::RationalFunctionField{T, U})(b::Rat{T, U}) where {T <: FieldElement, U <: Union{PolyElem{T}, MPolyElem{T}}}
   parent(b) != a && error("Unable to coerce rational function")
   return b
end

function (a::RationalFunctionField{T, U})(b::Integer) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   K = fraction_field(a)
   z = Rat{T, U}(K(b))
   z.parent = a
   return z
end

function (a::RationalFunctionField{T, U})(b::Rational{<:Integer}) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   K = fraction_field(a)
   z = Rat{T, U}(K(b))
   z.parent = a
   return z
end

function (a::RationalFunctionField)(b::RingElem)
   return a(fraction_field(a)(b))
end

###############################################################################
#
#   RationalFunctionField constructor
#
###############################################################################

function RationalFunctionField(k::Field, s::Symbol; cached=true)
   T = elem_type(k)

   R, x = AbstractAlgebra.PolynomialRing(k, s, cached=cached)

   U = elem_type(R)

   S = FractionField(R)
   g = S(x)
   t = Rat{T, U}(g)

   par_object = RationalFunctionField{T, U}(k, parent(g), s, cached)

   t.parent = par_object

   return par_object, t
end

function RationalFunctionField(k::Field, s::Vector{Symbol}; cached=true)
   T = elem_type(k)

   R, x = AbstractAlgebra.PolynomialRing(k, s, cached=cached)

   U = elem_type(R)

   S = FractionField(R)
   g = [S(xi) for xi in x]
   t = [Rat{T, U}(gi) for gi in g]

   par_object = RationalFunctionField{T, U}(k, parent(g[1]), s, cached)

   for ti in t
      ti.parent = par_object
   end

   return par_object, t
end
