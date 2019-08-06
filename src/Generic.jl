module Generic

import Base: rand

import AbstractAlgebra: Rationals, NCRing, NCRingElem, Ring, RingElem,
       RingElement

using AbstractAlgebra

include("generic/GenericTypes.jl")

include("generic/FreeModule.jl")

include("generic/Submodule.jl")

end # generic
