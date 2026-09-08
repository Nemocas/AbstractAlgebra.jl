###############################################################################
#
#   DirectSumModule.jl : Direct sums of modules
#
###############################################################################

###############################################################################
#
#   DirectSum constructor
#
###############################################################################

@doc raw"""
    direct_sum(m::Vector{<:FPModule{T}}) where T <: RingElement
    direct_sum(vals::FPModule{T}...) where T <: RingElement

Return a tuple $M, f, g$ consisting of $M$ the direct sum of the modules `m`
(supplied as a vector of modules), a vector $f$ of the injections
of the $m[i]$ into $M$ and a vector $g$ of the projections from
$M$ onto the $m[i]$.

# Examples

```jldoctest
julia> F = free_module(ZZ, 5)
Free module of rank 5 over integers

julia> m1 = F(BigInt[4, 7, 8, 2, 6])
(4, 7, 8, 2, 6)

julia> m2 = F(BigInt[9, 7, -2, 2, -4])
(9, 7, -2, 2, -4)

julia> S1, f1 = sub(F, [m1, m2])
(Submodule over integers with 2 generators and no relations, Hom: S1 -> F)

julia> m1 = F(BigInt[3, 1, 7, 7, -7])
(3, 1, 7, 7, -7)

julia> m2 = F(BigInt[-8, 6, 10, -1, 1])
(-8, 6, 10, -1, 1)

julia> S2, f2 = sub(F, [m1, m2])
(Submodule over integers with 2 generators and no relations, Hom: S2 -> F)

julia> m1 = F(BigInt[2, 4, 2, -3, -10])
(2, 4, 2, -3, -10)

julia> m2 = F(BigInt[5, 7, -6, 9, -5])
(5, 7, -6, 9, -5)

julia> S3, f3 = sub(F, [m1, m2])
(Submodule over integers with 2 generators and no relations, Hom: S3 -> F)

julia> D, f = direct_sum(S1, S2, S3)
(DirectSumModule over integers, AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}[Hom: S1 -> D, Hom: S2 -> D, Hom: S3 -> D], AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}[Hom: D -> S1, Hom: D -> S2, Hom: D -> S3])
```
"""
direct_sum(m::Vector{<:FPModule{<:RingElement}}) = Generic.direct_sum(m)
direct_sum(m::FPModule{<:RingElement}...) = Generic.direct_sum([m...])
