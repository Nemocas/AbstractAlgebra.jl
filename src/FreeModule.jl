###############################################################################
#
#   FreeModule.jl : Free modules over rings
#
###############################################################################

###############################################################################
#
#   FreeModule constructor
#
###############################################################################

@doc raw"""
    free_module(R::NCRing, rank::Int; cached::Bool = true)

Return the free module over the ring $R$ with the given rank.
"""
function free_module(R::NCRing, rank::Int; cached::Bool = true)
   return Generic.FreeModule(R, rank; cached=cached)
end

###############################################################################
#
#   VectorSpace constructor
#
###############################################################################

@doc raw"""
    vector_space(R::Field, dim::Int; cached::Bool = true)

Return the vector space over the field $R$ with the given dimension.
"""
function vector_space(R::Field, dim::Int; cached::Bool = true)
   Generic.FreeModule(R, dim; cached=cached)
end

@doc raw"""
    vector_space_dim(M::FPModule)

Return the dimension of the given vector space over `base_ring(M)`.
This method is only supported, when `base_ring(M)` is a field.

# Examples
```jldoctest
julia> M = free_module(QQ, 2)
Vector space of dimension 2 over rationals

julia> vector_space_dim(M)
2

julia> m = M([1, 2])
(1//1, 2//1)

julia> N, = sub(M, [m])
(Subspace over rationals with 1 generator and no relations, Hom: N -> M)

julia> vector_space_dim(N)
1

julia> Q, = quo(M, N)
(Quotient space over rationals with 1 generator and no relations, Hom: M -> Q)

julia> vector_space_dim(Q)
1
```
"""
function vector_space_dim end
