###############################################################################
#
#   RationalFunctionField.jl : Generic rational function fields
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{RationalFunctionFieldElem{T, U}}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}} = RationalFunctionField{T, U}

elem_type(::Type{RationalFunctionField{T, U}}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}} = RationalFunctionFieldElem{T, U}

base_ring_type(::Type{RationalFunctionField{T, U}}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}} = parent_type(T)

base_ring(a::RationalFunctionField{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}} = a.base_ring::parent_type(T)

parent(a::RationalFunctionFieldElem) = a.parent

data(x::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}} = x.d::FracFieldElem{U}

function fraction_field(a::RationalFunctionField{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   return a.fraction_field::Union{FracField{U}}
end

function is_domain_type(::Type{S}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}, S <: RationalFunctionFieldElem{T, U}}
   return true
end

function is_exact_type(a::Type{S}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}, S <: RationalFunctionFieldElem{T, U}}
   return is_exact_type(T)
end

function characteristic(R::RationalFunctionField)
   return characteristic(base_ring(R))
end

function check_parent(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}, throw::Bool = true) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   fl = parent(a) != parent(b)
   fl && throw && error("Incompatible rings in rationa function field operation")
   return !fl
end

###############################################################################
#
#   Constructors
#
###############################################################################

function //(x::RationalFunctionFieldElem{T, U}, y::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   check_parent(x, y)
   R = parent(x)
   return R(divexact(data(x), data(y)))
end

function //(x::T, y::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(y)
   parent(x) != base_ring(y) && error("Incompatible elements")
   return R(divexact(x, data(y)))
end

function //(x::RationalFunctionFieldElem{T, U}, y::T) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(x)
   base_ring(x) != parent(y) && error("Incompatible elements")
   return R(divexact(data(x), y))
end

function //(x::Rational{BigInt}, y::RationalFunctionFieldElem{Rational{BigInt}, U}) where U <: Union{PolyRingElem, MPolyRingElem}
   R = parent(y)
   return R(divexact(x, data(y)))
end

function //(x::RationalFunctionFieldElem{Rational{BigInt}, U}, y::Rational{BigInt}) where U <: Union{PolyRingElem, MPolyRingElem}
   R = parent(x)
   return R(divexact(data(x), y))
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::RationalFunctionFieldElem, h::UInt)
   b = 0x2d122a968560a3c0%UInt
   return xor(b, hash(data(a), h))
end

function Base.numerator(a::RationalFunctionFieldElem, canonicalise::Bool=true)
   return numerator(data(a), canonicalise)
end

function Base.denominator(a::RationalFunctionFieldElem, canonicalise::Bool=true)
   return denominator(data(a), canonicalise)
end

zero(R::RationalFunctionField) = R()

one(R::RationalFunctionField) = R(1)

iszero(a::RationalFunctionFieldElem) = iszero(data(a))

isone(a::RationalFunctionFieldElem) = isone(data(a))

is_unit(a::RationalFunctionFieldElem) = is_unit(data(a))

gen(R::RationalFunctionField) = R(gen(base_ring(R.fraction_field)))

gen(R::RationalFunctionField, i::Int) = R(gen(base_ring(R.fraction_field), i))

gens(R::RationalFunctionField) = R.(gens(base_ring(R.fraction_field)))

number_of_generators(R::RationalFunctionField) = number_of_generators(base_ring(R.fraction_field))

function deepcopy_internal(a::RationalFunctionFieldElem, dict::IdDict)
   R = parent(a)
   return R(deepcopy_internal(data(a), dict))
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::RationalFunctionFieldElem) = a

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function AbstractAlgebra.expressify(a::RationalFunctionFieldElem; context = nothing)
   d = data(a)
   return expressify(d; context)
end

function show(io::IO, ::MIME"text/plain", a::RationalFunctionFieldElem)
   print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function show(io::IO, a::RationalFunctionFieldElem)
   print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function show(io::IO, ::MIME"text/plain", a::RationalFunctionField)
  println(io, "Rational function field")
  io = pretty(io)
  print(io, Indent(), "over ", Lowercase(), base_ring(a))
  print(io, Dedent())
end

function show(io::IO, a::RationalFunctionField)
  if get(io, :supercompact, false)
    # no nested printing
    print(io, "Rational function field")
  else
    io = pretty(io) # we need this to allow printing lowercase
    print(IOContext(io, :supercompact => true),
          "Rational function field over ", Lowercase(), base_ring(a))
  end
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::RationalFunctionFieldElem)
   R = parent(a)
   return R(-data(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   check_parent(a, b)
   R = parent(a)
   return R(data(a) + data(b))
end

function -(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   check_parent(a, b)
   R = parent(a)
   return R(data(a) - data(b))
end

function *(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   check_parent(a, b)
   R = parent(a)
   return R(data(a) * data(b))
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::RationalFunctionFieldElem, b::Union{Integer, Rational, AbstractFloat})
   R = parent(a)
   return R(data(a)*b)
end

function *(a::Union{Integer, Rational, AbstractFloat}, b::RationalFunctionFieldElem)
   R = parent(b)
   return R(a*data(b))
end

function *(a::RationalFunctionFieldElem{T, U}, b::T) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(a)
   return R(data(a)*b)
end

function *(a::T, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(b)
   return R(a*data(b))
end

function +(a::RationalFunctionFieldElem, b::Union{Integer, Rational, AbstractFloat})
   R = parent(a)
   return R(data(a) + b)
end

function +(a::Union{Integer, Rational, AbstractFloat}, b::RationalFunctionFieldElem)
   R = parent(b)
   return R(a + data(b))
end

function +(a::RationalFunctionFieldElem{T, U}, b::T) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(a)
   return R(data(a) + b)
end

function +(a::T, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(b)
   return R(a + data(b))
end

function -(a::RationalFunctionFieldElem, b::Union{Integer, Rational, AbstractFloat})
   R = parent(a)
   return R(data(a) - b)
end

function -(a::Union{Integer, Rational, AbstractFloat}, b::RationalFunctionFieldElem)
   R = parent(b)
   return R(a - data(b))
end

function -(a::RationalFunctionFieldElem{T, U}, b::T) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(a)
   return R(data(a) - b)
end

function -(a::T, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(b)
   return R(a - data(b))
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function ==(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   check_parent(a, b)
   return data(a) == data(b)
end

function isequal(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   check_parent(a, b)
   return data(a) == data(b)
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function ==(a::RationalFunctionFieldElem, b::Union{Integer, Rational, AbstractFloat})
   R = parent(a)
   return data(a) == b
end

function ==(a::Union{Integer, Rational, AbstractFloat}, b::RationalFunctionFieldElem)
   R = parent(b)
   return a == data(b)
end

function ==(a::RationalFunctionFieldElem{T, U}, b::T) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(a)
   return data(a) == b
end

function ==(a::T, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(b)
   return a == data(b)
end

function ==(a::RationalFunctionFieldElem{T, U}, b::Poly{T}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   return a == parent(a)(b)
end

function ==(a::Poly{T}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   return parent(b)(a) == b
end

###############################################################################
#
#   Inversion
#
###############################################################################

function Base.inv(a::RationalFunctionFieldElem)
   R = parent(a)
   return R(inv(data(a)))
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}; check::Bool=true) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   check_parent(a, b)
   R = parent(a)
   return R(divexact(data(a), data(b); check=check))
end

function divides(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
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

function divexact(a::RationalFunctionFieldElem, b::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
   R = parent(a)
   return R(divexact(data(a), b; check=check))
end

function divexact(a::Union{Integer, Rational, AbstractFloat}, b::RationalFunctionFieldElem; check::Bool=true)
   R = parent(b)
   return R(divexact(a, data(b); check=check))
end

function divexact(a::RationalFunctionFieldElem{T, U}, b::T; check::Bool=true) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(a)
   return R(divexact(data(a), b; check=check))
end

function divexact(a::T, b::RationalFunctionFieldElem{T, U}; check::Bool=true) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}}
   R = parent(b)
   return R(divexact(a, data(b); check=check))
end

##############################################################################
#
#  Evaluation
#
##############################################################################

function evaluate(f::RationalFunctionFieldElem{T, U}, v::V) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}, V <: RingElement}
    return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

function evaluate(f::RationalFunctionFieldElem{T, U}, v::Vector{V}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}, V <: RingElement}
   return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::RationalFunctionFieldElem, b::Int)
   R = parent(a)
   return R(data(a)^b)
end

##############################################################################
#
#  Derivative
#
##############################################################################

function derivative(f::RationalFunctionFieldElem)
   R = parent(f)
   return R(derivative(data(f)))
end

function derivative(f::RationalFunctionFieldElem, x::MPolyRingElem)
   R = parent(f)
   return R(derivative(data(f), x))
end

###############################################################################
#
#   Square root
#
###############################################################################

@doc raw"""
    is_square(a::RationalFunctionFieldElem)

Return `true` if $a$ is a square.
"""
function is_square(a::RationalFunctionFieldElem)
   return is_square(data(a))
end

@doc raw"""
    Base.sqrt(a::RationalFunctionFieldElem; check::Bool=true)

Return the square root of $a$. By default the function will throw an exception
if the input is not square. If `check=false` this test is omitted.
"""
function Base.sqrt(a::RationalFunctionFieldElem; check::Bool=true)
   R = parent(a)
   return R(sqrt(data(a); check=check))
end

###############################################################################
#
#   GCD
#
###############################################################################

@doc raw"""
    gcd(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}

Return a greatest common divisor of $a$ and $b$ if one exists. N.B: we define
the GCD of $a/b$ and $c/d$ to be gcd$(ad, bc)/bd$, reduced to lowest terms.
"""
function gcd(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   check_parent(a, b)
   R = parent(a)
   return R(gcd(data(a), data(b)))
end

###############################################################################
#
#   Unsafe operators and functions
#
###############################################################################

function zero!(c::RationalFunctionFieldElem)
   c.d = zero!(data(c))
   return c
end

function mul!(c::RationalFunctionFieldElem{T, U}, a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   c.d = mul!(data(c), data(a), data(b))
   return c
end

function addeq!(a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   a.d = addeq!(data(a), data(b))
   return a
end

function add!(c::RationalFunctionFieldElem{T}, a::RationalFunctionFieldElem{T, U}, b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
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
      Make(S, vs[1]) # forward to default Make constructor
   else
      Make(S, make(R, vs...))
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

promote_rule(::Type{RationalFunctionFieldElem{T, U}}, ::Type{RationalFunctionFieldElem{T, U}}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}} = RationalFunctionFieldElem{T, U}

promote_rule(::Type{RationalFunctionFieldElem{T, U}}, ::Type{RationalFunctionFieldElem{T, U}}) where {T <: FieldElem, U <: Union{PolyRingElem, MPolyRingElem}} = RationalFunctionFieldElem{T, U}

function promote_rule(::Type{RationalFunctionFieldElem{T, U}}, ::Type{V}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}, V <: RingElem}
   promote_rule(FracFieldElem{U}, V) === FracFieldElem{U} ? RationalFunctionFieldElem{T, U} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::RationalFunctionField{T, U})() where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   K = fraction_field(a)
   z = RationalFunctionFieldElem{T, U}(K())
   z.parent = a
   return z
end

function (a::RationalFunctionField{T, U})(b::FracFieldElem{U}) where {T <: FieldElement, U <: Union{PolyRingElem{T}, MPolyRingElem{T}}}
   K = fraction_field(a)
   parent(b) != K && error("Unable to coerce rational function")
   z = RationalFunctionFieldElem{T, U}(b)
   z.parent = a
   return z::RationalFunctionFieldElem{T, U}
end

function (a::RationalFunctionField{T, U})(n::U, d::U) where {T <: FieldElement, U <: Union{PolyRingElem{T}, MPolyRingElem{T}}}
   R = parent(n)
   g = gcd(n, d)
   if !isone(g)
      n = divexact(n, g)
      d = divexact(d, g)
   end
   r = FracFieldElem{U}(n, d)
   try
      r.parent = FracDict[R]
   catch
      r.parent = fraction_field(R)
   end
   return a(r)
end

function (a::RationalFunctionField{T, U})(b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem{T}, MPolyRingElem{T}}}
   parent(b) != a && error("Unable to coerce rational function")
   return b
end

function (a::RationalFunctionField{T, U})(b::Integer) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   K = fraction_field(a)
   z = RationalFunctionFieldElem{T, U}(K(b))
   z.parent = a
   return z
end

function (a::RationalFunctionField{T, U})(b::Rational{<:Integer}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   K = fraction_field(a)
   z = RationalFunctionFieldElem{T, U}(K(b))
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

function rational_function_field(k::Field, s::Symbol; cached::Bool=true)
   T = elem_type(k)

   R, x = AbstractAlgebra.polynomial_ring(k, s, cached=cached)

   U = elem_type(R)

   S = fraction_field(R)
   g = S(x)
   t = RationalFunctionFieldElem{T, U}(g)

   par_object = RationalFunctionField{T, U}(k, parent(g), s, cached)

   t.parent = par_object

   return par_object, t
end

function rational_function_field(k::Field, s::Vector{Symbol}; cached::Bool=true)
   T = elem_type(k)

   R, x = AbstractAlgebra.polynomial_ring(k, s, cached=cached)

   U = elem_type(R)

   S = fraction_field(R)
   g = [S(xi) for xi in x]
   t = [RationalFunctionFieldElem{T, U}(gi) for gi in g]

   par_object = RationalFunctionField{T, U}(k, parent(g[1]), s, cached)

   for ti in t
      ti.parent = par_object
   end

   return par_object, t
end
