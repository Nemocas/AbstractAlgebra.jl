export RationalFunctionField

parent_type(::Type{Rat{T}}) where T <: FieldElement = RationalFunctionField{T}

elem_type(::Type{RationalFunctionField{T}}) where T <: FieldElement = Rat{T}

base_ring(a::RationalFunctionField{T}) where T <: FieldElement = a.base_ring::parent_type(T)

base_ring(a::Rat) = base_ring(parent(a))

parent(a::Rat) = a.parent

data(x::Rat{T}) where T <: FieldElement = x.d::Frac{dense_poly_type(T)}

function fraction_field(a::RationalFunctionField{T}) where T <: FieldElement
   return a.fraction_field::FracField{dense_poly_type(T)}
end

function check_parent(a::Rat{T}, b::Rat{T}, throw::Bool = true) where T <: FieldElement
   fl = parent(a) != parent(b)
   fl && throw && error("Incompatible rings in rationa function field operation")
   return !fl
end

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

gen(R::RationalFunctionField) = R(gen(base_ring(R.fraction_field)))

function deepcopy_internal(a::Rat, dict::IdDict)
   R = parent(a)
   return R(deepcopy_internal(data(a), dict))
end

canonical_unit(a::Rat) = a

function -(a::Rat)
   R = parent(a)
   return R(-data(a))
end

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
   return R(a + data(b))
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

function ==(a::Rat{T}, b::Rat{T}) where T <: FieldElement
   check_parent(a, b)
   return data(a) == data(b)
end

function isequal(a::Rat{T}, b::Rat{T}) where T <: FieldElement
   check_parent(a, b)
   return data(a) == data(b)
end

promote_rule(::Type{Rat{T}}, ::Type{Rat{T}}) where T <: FieldElement = Rat{T}

promote_rule(::Type{Rat{T}}, ::Type{Rat{T}}) where T <: FieldElem = Rat{T}

function promote_rule(::Type{Rat{T}}, ::Type{U}) where {T <: FieldElement, U <: RingElem}
   promote_rule(Frac{dense_poly_type(T)}, U) === Frac{dense_poly_type(T)} ? Rat{T} : Union{}
end

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
   return z::Rat{T}
end

function (a::RationalFunctionField{T})(n::S, d::S) where {T <: FieldElement, S <: PolyElem{T}}
   R = parent(n)
   g = gcd(n, d)
   if !isone(g)
      n = divexact(n, g)
      d = divexact(d, g)
   end
   r = Frac{S}(n, d)
   try
      r.parent = FracDict[R]
   catch
      r.parent = FractionField(R)
   end
   return a(r)
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

function (a::RationalFunctionField{T})(b::Rational{<:Integer}) where T <: FieldElement
   K = fraction_field(a)
   z = Rat{T}(K(b))
   z.parent = a
   return z
end

function (a::RationalFunctionField)(b::RingElem)
   return a(fraction_field(a)(b))
end

function RationalFunctionField(k::Field, s::Symbol; cached=true)
   T = elem_type(k)

   R, x = AbstractAlgebra.PolynomialRing(k, s, cached=cached)
   g = x//1
   t = Rat{T}(g)

   par_object = RationalFunctionField{T}(k, parent(g), s, cached)

   t.parent = par_object

   return par_object, t
end
