###############################################################################
#
#   SNFModule.jl : Invariant factor decomposition of modules
#
###############################################################################

###############################################################################
#
#   SNFModule constructor
#
###############################################################################

@doc raw"""
    invariant_factors(m::FPModule{T}) where T <: RingElement

Return a vector of the invariant factors of the module $M$.

# Examples

```jldoctest; setup = :(import Random; Random.seed!(42))
julia> M = free_module(ZZ, 3)
Free module of rank 3 over integers

julia> m1 = rand(M, -10:10)
(3, -1, 0)

julia> m2 = rand(M, -10:10)
(4, 4, -7)

julia> S, f = sub(M, [m1, m2])
(Submodule over integers with 2 generators and no relations, Hom: S -> M)

julia> Q, g = quo(M, S)
(Quotient module over integers with 2 generators and relations:
[16 -21], Hom: M -> Q)

julia> invariant_factors(Q)
1-element Vector{BigInt}:
 0
```
"""
function invariant_factors(m::FPModule{T}) where T <: RingElement
   return Generic.invariant_factors(m)
end

