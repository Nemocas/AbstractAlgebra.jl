@doc raw"""
Permutation group ring constructor.

See also: [`perm_group_ring`](@ref).
"""
function perm_group_ring(R::Ring, l::Int, cached::Bool=true)
    Generic.perm_group_ring(R, l, cached)
end
