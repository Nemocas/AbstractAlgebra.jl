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

function underlying_fraction_field(a::RationalFunctionField{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   return a.fraction_field::FracField{U}
end

function is_domain_type(::Type{S}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}, S <: RationalFunctionFieldElem{T, U}}
   return true
end

function is_exact_type(a::Type{S}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}, S <: RationalFunctionFieldElem{T, U}}
   return is_exact_type(T)
end

function symbols(a::RationalFunctionField)
  S = a.S
  if S isa Symbol
    return [S]
  else
    return S
  end
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

function (R::AbstractAlgebra.PolyRing{T})(x::RationalFunctionFieldElem{T, U}) where {T <: RingElement, U}
   @assert isone(denominator(x))
   y = numerator(x)
   @assert parent(y) === R
   return y
end

# Avoid ambiguity with (::Generic.PolyRing{T})(::RingElement)
function (R::PolyRing{T})(x::RationalFunctionFieldElem{T, U}) where {T <: RingElement, U}
   @assert isone(denominator(x))
   y = numerator(x)
   @assert parent(y) === R
   return y
end

zero(R::RationalFunctionField) = R()

one(R::RationalFunctionField) = R(1)

iszero(a::RationalFunctionFieldElem) = iszero(data(a))

isone(a::RationalFunctionFieldElem) = isone(data(a))

gen(R::RationalFunctionField) = R(gen(base_ring(underlying_fraction_field(R))))

gen(R::RationalFunctionField, i::Int) = R(gen(base_ring(underlying_fraction_field(R)), i))

gens(R::RationalFunctionField) = R.(gens(base_ring(underlying_fraction_field(R))))

number_of_generators(R::RationalFunctionField) = number_of_generators(base_ring(underlying_fraction_field(R)))

function deepcopy_internal(a::RationalFunctionFieldElem, dict::IdDict)
   R = parent(a)
   return R(deepcopy_internal(data(a), dict))
end

function characteristic(R::RationalFunctionField)
   return characteristic(base_ring(R))
end

is_finite(R::RationalFunctionField) = is_finite(base_ring(underlying_fraction_field(R)))

function is_perfect(R::RationalFunctionField)
  if characteristic(R) == 0
    return true
  end
  # char > 0
  if number_of_generators(R) > 0
    return false
  end

  return is_perfect(base_ring(R))
end

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
  @show_name(io, a)
  @show_special(io, a)
  if is_terse(io)
    # no nested printing
    print(io, "Rational function field")
  else
    io = pretty(io) # we need this to allow printing lowercase
    print(terse(io),
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

function *(a::RationalFunctionFieldElem, b::JuliaRingElement)
   R = parent(a)
   return R(data(a)*b)
end

function *(a::JuliaRingElement, b::RationalFunctionFieldElem)
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

function +(a::RationalFunctionFieldElem, b::JuliaRingElement)
   R = parent(a)
   return R(data(a) + b)
end

function +(a::JuliaRingElement, b::RationalFunctionFieldElem)
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

function -(a::RationalFunctionFieldElem, b::JuliaRingElement)
   R = parent(a)
   return R(data(a) - b)
end

function -(a::JuliaRingElement, b::RationalFunctionFieldElem)
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

function ==(a::RationalFunctionFieldElem, b::JuliaRingElement)
   R = parent(a)
   return data(a) == b
end

function ==(a::JuliaRingElement, b::RationalFunctionFieldElem)
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

function divexact(a::RationalFunctionFieldElem, b::JuliaRingElement; check::Bool=true)
   R = parent(a)
   return R(divexact(data(a), b; check=check))
end

function divexact(a::JuliaRingElement, b::RationalFunctionFieldElem; check::Bool=true)
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

function evaluate(f::RationalFunctionFieldElem, v::RingElement)
    return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

function evaluate(f::RationalFunctionFieldElem, v::Vector{<:RingElement})
   return evaluate(numerator(f), v)//evaluate(denominator(f), v)
end

function (a::RationalFunctionFieldElem)(val::RingElement)
   return evaluate(a, val)
end

function (a::RationalFunctionFieldElem)(vals::RingElement...)
   return evaluate(a, [vals...])
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

function is_square(a::RationalFunctionFieldElem)
   return is_square(data(a))
end

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
#   Unsafe functions
#
###############################################################################

zero!(c::RationalFunctionFieldElem) =
  RationalFunctionFieldElem(zero!(data(c)), parent(c))

one!(c::RationalFunctionFieldElem) =
  RationalFunctionFieldElem(one!(data(c)), parent(c))

neg!(c::T, a::T) where {T <: RationalFunctionFieldElem} =
  RationalFunctionFieldElem(neg!(data(c), data(a)), parent(c))

inv!(c::T, a::T) where {T <: RationalFunctionFieldElem} =
  RationalFunctionFieldElem(inv!(data(c), data(a)), parent(c))

add!(c::T, a::T, b::T) where {T <: RationalFunctionFieldElem} =
  RationalFunctionFieldElem(add!(data(c), data(a), data(b)), parent(c))

sub!(c::T, a::T, b::T) where {T <: RationalFunctionFieldElem} =
  RationalFunctionFieldElem(sub!(data(c), data(a), data(b)), parent(c))

mul!(c::T, a::T, b::T) where {T <: RationalFunctionFieldElem} =
  RationalFunctionFieldElem(mul!(data(c), data(a), data(b)), parent(c))

###############################################################################
#
#   Random functions
#
###############################################################################

RandomExtensions.maketype(R::RationalFunctionField, _) = elem_type(R)

function RandomExtensions.make(S::RationalFunctionField, vs...)
   R = base_ring(underlying_fraction_field(S))
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1]) # forward to default Make constructor
   else
      Make(S, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:FieldElement, <:RationalFunctionField}})
   S, v = sp[][1:end]
   R = base_ring(underlying_fraction_field(S))
   n = rand(rng, v)
   d = R()
   while iszero(d)
      d = rand(rng, v)
   end
   return S(underlying_fraction_field(S)(n, d; reduce = true))
end

rand(rng::AbstractRNG, S::RationalFunctionField, v...) =
   rand(rng, make(S, v...))

rand(S::RationalFunctionField, v...) = rand(Random.default_rng(), S, v...)

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::RationalFunctionField{Rational{BigInt}})
  rand(R, 0:3, -3:3)
end

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
#   Parent object call overload
#
###############################################################################

function (a::RationalFunctionField{T, U})() where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   K = underlying_fraction_field(a)
   z = RationalFunctionFieldElem{T, U}(K(), a)
   return z
end

function (a::RationalFunctionField{T, U})(b::FracFieldElem{U}) where {T <: FieldElement, U <: Union{PolyRingElem{T}, MPolyRingElem{T}}}
   K = underlying_fraction_field(a)
   parent(b) != K && error("Unable to coerce rational function")
   z = RationalFunctionFieldElem{T, U}(b, a)
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
   r.parent = get(FracDict, R) do
      return underlying_fraction_field(R)
   end
   return a(r)
end

function (a::RationalFunctionField{T, U})(b::RationalFunctionFieldElem{T, U}) where {T <: FieldElement, U <: Union{PolyRingElem{T}, MPolyRingElem{T}}}
   parent(b) != a && error("Unable to coerce rational function")
   return b
end

function (a::RationalFunctionField{T, U})(b::Integer) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   K = underlying_fraction_field(a)
   z = RationalFunctionFieldElem{T, U}(K(b), a)
   return z
end

function (a::RationalFunctionField{T, U})(b::Rational{<:Integer}) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   K = underlying_fraction_field(a)
   z = RationalFunctionFieldElem{T, U}(K(b), a)
   return z
end

function (a::RationalFunctionField)(b::RingElem)
   return a(underlying_fraction_field(a)(b))
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
   par_object = RationalFunctionField{T, U}(k, parent(g), s, cached)
   t = RationalFunctionFieldElem{T, U}(g, par_object)

   return par_object, t
end

function rational_function_field(k::Field, s::Vector{Symbol}; cached::Bool=true)
   T = elem_type(k)

   R, x = AbstractAlgebra.polynomial_ring(k, s, cached=cached)

   U = elem_type(R)

   S = fraction_field(R)
   g = [S(xi) for xi in x]
   par_object = RationalFunctionField{T, U}(k, S, s, cached)
   t = [RationalFunctionFieldElem{T, U}(gi, par_object) for gi in g]

   return par_object, t
end

################################################################################
#
#  Factoring
#
################################################################################

function _transport_factor(factorfun::F, f, fwd::G, bwd::H) where {F, G, H}
  R = parent(f)
  g = map_coefficients(fwd, f; cached = false)
  facg = factorfun(g)
  D = Dict{typeof(f), Int}()
  for (pg, e) in facg
    D[map_coefficients(bwd, pg; parent = R)] = e
  end
  return Fac(map_coefficients(bwd, unit(facg); parent = R), D)
end

function factor(f::Union{PolyRingElem{T}, MPolyRingElem{T}}) where {T <: RationalFunctionFieldElem}
  R = base_ring(f)
  return _transport_factor(factor, f, data, R)
end

function factor_squarefree(f::Union{PolyRingElem{T}, MPolyRingElem{T}}) where {T <: RationalFunctionFieldElem}
  R = base_ring(f)
  return _transport_factor(factor_squarefree, f, data, R)
end

function is_squarefree(f::Union{PolyRingElem{T}, MPolyRingElem{T}}) where {T <: RationalFunctionFieldElem}
  is_squarefree(map_coefficients(data, f; cached = false))
end

################################################################################
#
#  Generation of finite subset and rank interpolation
#
################################################################################

function evaluation_points(K::RationalFunctionField, n::Int)
   @req n > 0 "n must be positive."
   F = base_ring(K)
   v = evaluation_points(F, n)
   if is_empty(v)
      v = elem_type(K)[]
      base_v = evaluation_points(F, Int(order(F)))
      d = clog(order(F), ZZ(n))
      genKpowers = powers(gen(K), d - 1)
      for coeffs in ProductIterator(base_v, d; inplace=true)
         push!(v, sum(coeffs[j] * genKpowers[j] for j in 1:d))
         if length(v) == n
            break
         end
      end
   end
   return v
end

function rank_interpolation(M::MatrixElem{<: RationalFunctionFieldElem})
   n = nrows(M)
   m = ncols(M)
   if is_zero(n) || is_zero(m)
      return 0
   end
   Kx = base_ring(Generic.underlying_fraction_field(base_ring(M)))
   A = deepcopy(M)
   #All rows/columns of M are multiplied with polynomials such that all entries of M have denominator 1.
   #Then M can be considered as a matrix over a polynomial ring. 
   for i = 1:n
      for j = 1:m
         if !is_one(denominator(A[i, j]))
            if n < m
               multiply_column!(A, denominator(A[i, j]), j)
            else
               multiply_row!(A, denominator(A[i, j]), i)
            end
         end
      end
   end
   return rank_interpolation(matrix(Kx, n, m, [numerator(A[i, j]) for i in 1:n, j in 1:m]))
end

function rank_interpolation_mc(M::MatrixElem{<: RationalFunctionFieldElem}, err::Float64)
   n = nrows(M)
   m = ncols(M)
   if is_zero(n) || is_zero(m)
      return 0
   end
   Kx = base_ring(Generic.underlying_fraction_field(base_ring(M)))
   A = deepcopy(M)
   #All rows/columns of M are multiplied with polynomials such that all entries of M have denominator 1.
   #Then M can be considered as a matrix over some polynomial ring.
   for i = 1:n
      for j = 1:m
         if is_one(denominator(A[i, j]))
            if n < m
               multiply_column!(A, denominator(A[i, j]), j)
            else
               multiply_row!(A, denominator(A[i, j]), i)
            end
         end
      end
   end
   return rank_interpolation_mc(matrix(Kx, n, m, [numerator(A[i, j]) for i in 1:n, j in 1:m]), err)
end
