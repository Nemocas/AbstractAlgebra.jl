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

include("generic/UniversalRing.jl")

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

include("generic/FreeAssociativeAlgebra.jl")

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

include("generic/FreeAssociativeAlgebraGroebner.jl")

include("generic/PolyRingHom.jl")

###############################################################################
#
#   Temporary miscellaneous files being moved from Hecke.jl
#
###############################################################################

include("generic/Misc/Poly.jl")
include("generic/Misc/Rings.jl")
include("generic/Misc/Localization.jl")

# Deprecated in 0.35.*
Base.@deprecate_binding ResF EuclideanRingResidueFieldElem
Base.@deprecate_binding ResField EuclideanRingResidueField
Base.@deprecate_binding Res EuclideanRingResidueRingElem
Base.@deprecate_binding ResRing EuclideanRingResidueRing
Base.@deprecate_binding ResidueField EuclideanRingResidueField false
Base.@deprecate_binding ResidueFieldElem EuclideanRingResidueFieldElem
Base.@deprecate_binding ResidueRing EuclideanRingResidueRing false
Base.@deprecate_binding ResidueRingElem EuclideanRingResidueRingElem

# Deprecated in 0.41.*
@deprecate hom(M::DirectSumModule{T}, N::DirectSumModule{T}, mp::Vector{ModuleHomomorphism{T}}) where T hom_direct_sum(M, N, mp)
@deprecate hom(A::DirectSumModule{T}, B::DirectSumModule{T}, M::Matrix{<:Map{<:AbstractAlgebra.FPModule{T}, <:AbstractAlgebra.FPModule{T}}}) where {T} hom_direct_sum(A, B, M)

# Deprecated for 0.44
#Base.@deprecate_binding MatSpace AbstractAlgebra.MatSpace false

end # generic
