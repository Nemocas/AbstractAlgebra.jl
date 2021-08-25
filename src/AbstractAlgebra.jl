module AbstractAlgebra

using InteractiveUtils

import Base: inv, gcd, zero, iszero, length, ^, +, -, *, ==, //, !=

export elem_type, parent_type

export RingElem, FieldElem, RingElement

export PolyElem

export PolyRing

export ZZ, QQ

const CacheDictType = Dict

function get_cached!(default::Base.Callable, dict::AbstractDict,
                                             key,
                                             use_cache::Bool)
   return use_cache ? Base.get!(default, dict, key) : default()
end

include("AbstractTypes.jl")

const PolynomialElem{T} = PolyElem{T}

include("julia/JuliaTypes.jl")

include("Poly.jl")
include("RationalFunctionField.jl")
include("Fraction.jl")

include("Generic.jl")

import .Generic: Generic, elem_type, parent_type, fit!, coeff, setcoeff!, normalise,
                 set_length!, zero!, add!, addeq!, mul!

export Generic

include("Rings.jl")

const ZZ = JuliaZZ
const QQ = JuliaQQ

end # module
