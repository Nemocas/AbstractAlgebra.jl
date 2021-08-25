module Generic

import Base: Array, +, -, *, ==, ^, &, |, <<, >>, ~, <=, >=,
             <, >, //, /, !=

import ..AbstractAlgebra: CacheDictType, get_cached!

import ..AbstractAlgebra: Field, FieldElement, Integers,
                          Rationals, Ring, RingElem,
                          RingElement

import ..AbstractAlgebra: base_ring, canonical_unit, denominator,
                          numerator

using ..AbstractAlgebra

include("generic/GenericTypes.jl")

include("generic/Poly.jl")

include("generic/Fraction.jl")

include("generic/RationalFunctionField.jl")

end # generic
