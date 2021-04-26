###############################################################################
#
#   RationalFunctionField.jl : Generic rational function fields
#
###############################################################################

export RationalFunctionField

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{Rat{T}}) where T <: FieldElement = RationalFunctionField{T}

elem_type(::Type{RationalFunctionField{T}}) where T <: FieldElement = Rat{T}

base_ring(a::RationalFunctionField{T}) where T <: FieldElement = a.base_ring::parent_type(T)

base_ring(a::Rat) = base_ring(parent(a))

parent(a::Rat) = a.parent

fraction_field(a::RationalFunctionField{T}) where T <: FieldElement = a.fraction_field

function isdomain_type(::Type{T}) where {S <: FieldElement, T <: Rat{S}}
   return true
end

function isexact_type(a::Type{T}) where {S <: FieldElement, T <: Rat{S}}
   return isexact_type(S)
end

function characteristic(R::RationalFunctionField)
   return characteristic(base_ring(R))
end

function check_parent(a::Rat{T}, b::Rat{T}, throw::Bool = true) where T <: FieldElement
   fl = parent(a) != parent(b)
   fl && throw && error("Incompatible rings in rationa function field operation")
   return !fl
end

###############################################################################
#
#   Constructors
#
###############################################################################

function //(x::Rat{T}, y::Rat{T}) where T <: FieldElement
   check_parent(x, y)
   R = parent(x)
   return R(divexact(data(x), data(y)))
end

function //(x::T, y::Rat{T}) where T <: FieldElement
   R = parent(y)
   parent(x) != base_ring(y) && error("Incompatible elements")
   return R(divexact(x, data(y)))
end

function //(x::Rat{T}, y::T) where T <: FieldElement
   R = parent(x)
   base_ring(x) != parent(y) && error("Incompatible elements")
   return R(divexact(data(x), y))
end

function //(x::Rational{BigInt}, y::Rat{Rational{BigInt}})
   R = parent(y)
   return R(divexact(x, data(y)))
end

function //(x::Rat{Rational{BigInt}}, y::Rational{BigInt})
   R = parent(x)
   return R(divexact(data(x), y))
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

data(x::Rat) = x.d

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

isunit(a::Rat) = isunit(data(a))

function deepcopy_internal(a::Rat, dic::IdDict)
   R = parent(a)
   return R(deepcopy(data(a)))
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
   return expressify(d)
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

function +(a::Rat{T}, b::Rat{T}) where T <: FieldElement
   check_parent(a, b)
   R = parent(a)
   return R(data(a) + data(b))
end

function -(a::Rat{T}, b::Rat{T}) where T <: FieldElement
   check_parent(a, b)
   R = parent(a)
   return R(data(a) - data(b))
end

function *(a::Rat{T}, b::Rat{T}) where T <: FieldElement
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

function *(a::Rat{T}, b::T) where T <: FieldElem
   R = parent(a)
   return R(data(a)*b)
end

function *(a::T, b::Rat{T}) where T <: FieldElem
   R = parent(b)
   return R(a*data(b))
end

function +(a::Rat, b::Union{Integer, Rational, AbstractFloat})
   R = parent(a)
   return R(data(a) + b)
end

function +(a::Union{Integer, Rational, AbstractFloat}, b::Rat)
   R = parent(b)
   return R(a*data(b))
end

function +(a::Rat{T}, b::T) where T <: FieldElem
   R = parent(a)
   return R(data(a) + b)
end

function +(a::T, b::Rat{T}) where T <: FieldElem
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

function -(a::Rat{T}, b::T) where T <: FieldElem
   R = parent(a)
   return R(data(a) - b)
end

function -(a::T, b::Rat{T}) where T <: FieldElem
   R = parent(b)
   return R(a - data(b))
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function ==(a::Rat{T}, b::Rat{T}) where T <: FieldElement
   check_parent(a, b)
   return data(a) == data(b)
end

function isequal(a::Rat{T}, b::Rat{T}) where T <: FieldElement
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

function ==(a::Rat{T}, b::T) where T <: FieldElem
   R = parent(a)
   return data(a) == b
end

function ==(a::T, b::Rat{T}) where T <: FieldElem
   R = parent(b)
   return a == data(b)
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

function divexact(a::Rat{T}, b::Rat{T}) where T <: FieldElement
   check_parent(a, b)
   R = parent(a)
   return R(divexact(data(a), data(b)))
end

function divides(a::Rat{T}, b::Rat{T}) where T <: FieldElement
   return divides(data(a), data(b))
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(a::Rat, b::Union{Integer, Rational, AbstractFloat})
   R = parent(a)
   return R(divexact(data(a), b))
end

function divexact(a::Union{Integer, Rational, AbstractFloat}, b::Rat)
   R = parent(b)
   return R(divexact(a, data(b)))
end

function divexact(a::Rat{T}, b::T) where T <: FieldElem
   R = parent(a)
   return R(divexact(data(a), b))
end

function divexact(a::T, b::Rat{T}) where T <: FieldElem
   R = parent(b)
   return R(divexact(a, data(b)))
end

##############################################################################
#
#  Evaluation
#
##############################################################################

function evaluate(f::Rat{T}, v::U) where {T <: RingElement, U <: RingElement}
    return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

function evaluate(f::Rat{T}, v::U) where {T <: PolyElem, U <: Integer}
    return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::Rat{T}, b::Int) where T <: FieldElement
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

###############################################################################
#
#   Square root
#
###############################################################################

function issquare(a::Rat)
   return issquare(data(a))
end

function Base.sqrt(a::Rat)
   R = parent(a)
   return R(sqrt(data(a)))
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(a::Rat{T}, b::Rat{T}) where T <: FieldElement
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

function mul!(c::Rat{T}, a::Rat{T}, b::Rat{T}) where T <: FieldElement
   c.d = mul!(data(c), data(a), data(b))
   return c
end

function addeq!(a::Rat{T}, b::Rat{T}) where T <: FieldElement
   a.d = addeq!(data(a), data(b))
   return a
end

function add!(c::Rat{T}, a::Rat{T}, b::Rat{T}) where T <: FieldElement
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

promote_rule(::Type{Rat{T}}, ::Type{Rat{T}}) where T <: FieldElement = Rat{T}
promote_rule(::Type{Rat{T}}, ::Type{Rat{T}}) where T <: FieldElem = Rat{T}

function promote_rule(::Type{Rat{T}}, ::Type{U}) where {T <: FieldElement, U <: RingElem}
   promote_rule(T, U) == T ? Rat{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::RationalFunctionField{T})() where T <: FieldElement
   K = fraction_field(a)
   z = Rat{T}(K())
   z.parent = a
   return z
end

function (a::RationalFunctionField{T})(b::Frac{<:PolyElem{T}}) where T <: FieldElement
   K = fraction_field(a)
   parent(b) != K && error("Unable to coerce rational function")
   z = Rat{T}(b)
   z.parent = a
   return z
end

function (a::RationalFunctionField{T})(b::Rat{T}) where T <: FieldElement
   parent(b) != a && error("Unable to coerce rational function")
   return b
end

function (a::RationalFunctionField{T})(b::Integer) where T <: FieldElement
   K = fraction_field(a)
   z = Rat{T}(K(b))
   z.parent = a
   return z
end

function (a::RationalFunctionField{T})(b::Rational{Integer}) where T <: FieldElement
   K = fraction_field(a)
   z = Rat{T}(K(b))
   z.parent = a
   return z
end

###############################################################################
#
#   RationalFunctionField constructor
#
###############################################################################

function RationalFunctionField(k::Field, s::AbstractString; cached=true)
   T = elem_type(k)

   R, x = PolynomialRing(k, s, cached=cached)
   g = x//1
   t = Rat{T}(g)

   par_object = RationalFunctionField{T}(k, parent(g), Symbol(s), cached)

   t.parent = par_object

   return par_object, t
end

