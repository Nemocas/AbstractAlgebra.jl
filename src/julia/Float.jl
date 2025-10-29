###############################################################################
#
#   Float.jl : Additional AbstractAlgebra functionality for Julia AbstractFloat types
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

const JuliaRealField = Floats{BigFloat}()

const RDF = Floats{Float64}()

parent(a::T) where T <: AbstractFloat = Floats{T}()

elem_type(::Type{Floats{T}}) where T <: AbstractFloat = T

parent_type(::Type{T}) where T <: AbstractFloat = Floats{T}

base_ring_type(::Type{Floats{T}}) where T <: AbstractFloat = Union{}   # no base ring

is_domain_type(::Type{T}) where T <: AbstractFloat = true

is_exact_type(::Type{T}) where T <: AbstractFloat = false

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(::Floats{T}) where T <: AbstractFloat = T(0)

one(::Floats{T}) where T <: AbstractFloat = T(1)

is_unit(a::AbstractFloat) = !is_zero(a)

canonical_unit(a::AbstractFloat) = iszero(a) ? copysign(one(a), a) : a

characteristic(a::Floats{T}) where T <: AbstractFloat = 0

if VERSION < v"1.13.0-DEV.534" # https://github.com/JuliaLang/julia/pull/53677
  is_negative(n::T) where T<:Real = n < zero(T)
  is_positive(n::T) where T<:Real = n > zero(T)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::Floats{BigFloat})
   print(io, "Floats")
end

function show(io::IO, R::Floats{T}) where T
   print(io, "Floats{$T}()")
end

function expressify(a::AbstractFloat; context = nothing)
   if a < 0
      return Expr(:call, :-, string(-a))
   else
      return string(a)
   end
end

###############################################################################
#
#   Exact division
#
###############################################################################

divexact(a::T, b::T; check::Bool=true) where T <: AbstractFloat = a/b

divexact(a::AbstractFloat, b::Union{Signed, Unsigned}; check::Bool=true) = a/b

divexact(a::BigFloat, b::BigInt; check::Bool=true) = a/b

divexact(a::Union{Signed, Unsigned}, b::AbstractFloat; check::Bool=true) = a/b

divexact(a::BigInt, b::BigFloat; check::Bool=true) = a/b

divexact(a::AbstractFloat, b::Rational{T}; check::Bool=true) where T <: Union{Signed, Unsigned} = a/b

divexact(a::BigFloat, b::Rational{BigInt}; check::Bool=true) = a/b

divexact(a::Rational{T}, b::AbstractFloat; check::Bool=true) where T <: Union{Signed, Unsigned} = a/b

divexact(a::Rational{BigInt}, b::BigFloat; check::Bool=true) = a/b

function divides(a::BigFloat, b::BigFloat)
   if iszero(b)
      return false, BigFloat(0)
   else
      return true, divexact(a, b; check=false)
   end
end

function divrem(a::BigFloat, b::BigFloat)
   return a/b, 0
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(a::T, b::T) where T <: AbstractFloat
   if iszero(a) && iszero(b)
      return T(0)
   else
      return T(1)
   end
end

###############################################################################
#
#   Coprime bases
#
###############################################################################

function ppio(a::T, b::T) where T <: AbstractFloat
   if iszero(a) && iszero(b)
      return T(0), T(0)
   else
      return T(1), a
   end
end

###############################################################################
#
#   Square root
#
###############################################################################

function Base.sqrt(a::AbstractFloat; check::Bool=true)
   return Base.sqrt(a)
end

# See generic documentation in NCRings.jl
function is_square(a::AbstractFloat)
   return a >= 0
end

function is_square_with_sqrt(a::T) where T <: AbstractFloat
   if a >= 0
      return true, Base.sqrt(a)
   else
      return false, zero(T)
   end
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

# No actual mutation is permitted for Julia types
# See #1077

function addmul!(a::T, b::T, c::T) where T <: AbstractFloat # special case, no temporary required
   return a + b*c
end

###############################################################################
#
#   Random generation
#
###############################################################################

RandomExtensions.maketype(R::Floats{T}, _) where {T} = T

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{T, Floats{T},
                                         <:AbstractUnitRange{<:Union{AbstractFloat, Int}}}}
              ) where {T}
   R, n = sp[][1:end]
   R(first(n) + rand(rng, Float64)*(last(n) - first(n)))
end


rand(rng::AbstractRNG, R::Floats, n::AbstractUnitRange) = rand(rng, make(R, n))

rand(R::Floats, n) = rand(Random.default_rng(), R, n)

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::Floats{T}) where T
  return rand(T)*rand(-100:100)
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::Floats{T})() where T <: AbstractFloat
   return T(0)
end

function (a::Floats{T})(b::AbstractFloat) where T <: AbstractFloat
   return T(b)
end

function (a::Floats{T})(b::Integer) where T <: AbstractFloat
   return T(b)
end

function (a::Floats{BigFloat})(b::BigInt)
   return BigFloat(b)
end

function (a::Floats{T})(b::Rational{S}) where {S <: Union{Signed, Unsigned}, T <: AbstractFloat}
   return T(b)
end

function (a::Floats{BigFloat})(b::Rational{BigInt})
   return BigFloat(b)
end

