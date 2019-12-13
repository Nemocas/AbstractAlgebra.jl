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

JuliaRealField = Floats{BigFloat}()

RDF = Floats{Float64}()

parent(a::T) where T <: AbstractFloat = Floats{T}()

elem_type(::Type{Floats{T}}) where T <: AbstractFloat = T

parent_type(::Type{T}) where T <: AbstractFloat = Floats{T}

base_ring(a::Floats{T}) where T <: AbstractFloat = Union{}

isdomain_type(::Type{T}) where T <: AbstractFloat = true

isexact_type(::Type{T}) where T <: AbstractFloat = false

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(::Floats{T}) where T <: AbstractFloat = T(0)

one(::Floats{T}) where T <: AbstractFloat = T(1)

isunit(a::AbstractFloat) = a != 0

canonical_unit(a::AbstractFloat) = a

characteristic(a::Floats{T}) where T <: AbstractFloat = 0

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::Floats)
   print(io, "Floats")
end

needs_parentheses(::AbstractFloat) = false

displayed_with_minus_in_front(a::AbstractFloat) = a < 0

show_minus_one(::Type{T}) where T <: AbstractFloat = false

###############################################################################
#
#   Exact division
#
###############################################################################

divexact(a::T, b::T) where T <: AbstractFloat = a/b

divexact(a::AbstractFloat, b::Union{Signed, Unsigned}) = a/b

divexact(a::BigFloat, b::BigInt) = a/b

divexact(a::Union{Signed, Unsigned}, b::AbstractFloat) = a/b

divexact(a::BigInt, b::BigFloat) = a/b

divexact(a::AbstractFloat, b::Rational{T}) where T <: Union{Signed, Unsigned} = a/b

divexact(a::BigFloat, b::Rational{BigInt}) = a/b

divexact(a::Rational{T}, b::AbstractFloat) where T <: Union{Signed, Unsigned} = a/b

divexact(a::Rational{BigInt}, b::BigFloat) = a/b

function divides(a::BigFloat, b::BigFloat)
   if b == 0
      return false, BigFloat(0)
   else
      return true, divexact(a, b)
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
   if a == 0 && b == 0
      return T(0)
   else
      return T(1)
   end
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::AbstractFloat)
   return 0
end

function zero!(a::BigFloat)
   ccall((:mpfr_set_si, :libmpfr), Nothing,
         (Ref{BigFloat}, Int, Int32), a, 0, Base.MPFR.ROUNDING_MODE[])
   return a
end

function mul!(a::T, b::T, c::T) where T <: AbstractFloat
   return b*c
end

function mul!(a::BigFloat, b::BigFloat, c::BigFloat)
   ccall((:mpfr_mul, :libmpfr), Nothing,
         (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
                 a, b, c, Base.MPFR.ROUNDING_MODE[])
   return a
end

function add!(a::T, b::T, c::T) where T <: AbstractFloat
   return b + c
end

function add!(a::BigFloat, b::BigFloat, c::BigFloat)
   ccall((:mpfr_add, :libmpfr), Nothing,
         (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
                 a, b, c, Base.MPFR.ROUNDING_MODE[])
   return a
end

function addeq!(a::T, b::T) where T <: AbstractFloat
   return a + b
end

function addeq!(a::BigFloat, b::BigFloat)
   ccall((:mpfr_add, :libmpfr), Nothing,
         (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
                 a, a, b, Base.MPFR.ROUNDING_MODE[])
   return a
end

function addmul!(a::T, b::T, c::T, d::T) where T <: AbstractFloat
   return a + b*c
end

function addmul!(a::BigFloat, b::BigFloat, c::BigFloat, d::BigFloat)
   ccall((:mpfr_fma, :libmpfr), Nothing,
         (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
                 a, b, c, a, Base.MPFR.ROUNDING_MODE[])
   return a
end

function addmul!(a::T, b::T, c::T) where T <: AbstractFloat # special case, no temporary required
   return a + b*c
end

function addmul!(a::BigFloat, b::BigFloat, c::BigFloat) # special case, no temporary required
   ccall((:mpfr_fma, :libmpfr), Nothing,
         (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
                 a, b, c, a, Base.MPFR.ROUNDING_MODE[])
   return a
end

###############################################################################
#
#   Random generation
#
###############################################################################

function rand(rng::AbstractRNG, R::Floats, n::UnitRange{<:AbstractFloat})
   return R(n.start + rand(rng, Float64)*(n.stop - n.start))
end

function rand(rng::AbstractRNG, R::Floats, n::UnitRange{Int})
   return R(n.start + rand(rng, Float64)*(n.stop - n.start))
end

rand(R::Floats, n) = rand(Random.GLOBAL_RNG, R, n)

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

function (a::Floats{T})(b::Union{Signed, Unsigned}) where T <: AbstractFloat
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
