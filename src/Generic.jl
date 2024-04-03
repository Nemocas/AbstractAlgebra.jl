module Generic

###############################################################################
#
#   Imports and exports
#
###############################################################################

include("generic/imports.jl")

include("generic/exports.jl")

###############################################################################
#
#   All functionality
#
###############################################################################

include("generic/GenericTypes.jl")

include("generic/PermGroups.jl")

include("generic/YoungTabs.jl")

include("generic/Residue.jl")

include("generic/ResidueField.jl")

include("generic/Poly.jl")

include("generic/NCPoly.jl")

include("generic/MPoly.jl")

include("generic/UnivPoly.jl")

include("generic/SparsePoly.jl")

include("generic/LaurentPoly.jl")

include("generic/LaurentMPoly.jl")

include("generic/RelSeries.jl")

include("generic/AbsSeries.jl")

include("generic/AbsMSeries.jl")

include("generic/LaurentSeries.jl")

include("generic/PuiseuxSeries.jl")

include("generic/Matrix.jl")

include("generic/MatRing.jl")

include("generic/FreeAssAlgebra.jl")

include("generic/Fraction.jl")

include("generic/TotalFraction.jl")

include("generic/FactoredFraction.jl")

include("generic/RationalFunctionField.jl")

include("generic/FunctionField.jl")

include("generic/FreeModule.jl")

include("generic/Submodule.jl")

include("generic/QuotientModule.jl")

include("generic/DirectSum.jl")

include("generic/ModuleHomomorphism.jl")

include("generic/Module.jl")

include("generic/InvariantFactorDecomposition.jl")

include("generic/Map.jl")

include("generic/MapWithInverse.jl")

include("generic/MapCache.jl")

include("generic/Ideal.jl")

include("generic/AhoCorasick.jl")

include("generic/FreeAssAlgebraGroebner.jl")

###############################################################################
#
#   Temporary miscellaneous files being moved from Hecke.jl
#
###############################################################################

include("generic/Misc/Poly.jl")
include("generic/Misc/Rings.jl")
include("generic/Misc/Localization.jl")

# TODO/FIXME: deprecate aliases, remove in the future
import ..AbstractAlgebra: @alias
# Deprecated in 0.35.*
Base.@deprecate_binding ResF EuclideanRingResidueFieldElem
Base.@deprecate_binding ResField EuclideanRingResidueField
Base.@deprecate_binding Res EuclideanRingResidueRingElem
Base.@deprecate_binding ResRing EuclideanRingResidueRing
Base.@deprecate_binding ResidueField EuclideanRingResidueField
Base.@deprecate_binding ResidueFieldElem EuclideanRingResidueFieldElem
Base.@deprecate_binding ResidueRing EuclideanRingResidueRing
Base.@deprecate_binding ResidueRingElem EuclideanRingResidueRingElem

end # generic
