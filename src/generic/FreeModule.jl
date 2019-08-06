export FreeModule, free_module_elem

parent_type(::Type{free_module_elem{T}}) where T <: RingElement = FreeModule{T}

base_ring(M::FreeModule{T}) where T <: Union{RingElement, NCRingElem} = M.base_ring::parent_type(T)

elem_type(::Type{FreeModule{T}}) where T <: Union{RingElement, NCRingElem} = free_module_elem{T}

rank(M::FreeModule{T}) where T <: Union{RingElement, NCRingElem} = M.rank

ngens(M::FreeModule{T}) where T <: Union{RingElement, NCRingElem} = M.rank

function rand(M::AbstractAlgebra.FPModule{T}, vals...) where T <: RingElement
   R = base_ring(M)
   v = [rand(R, vals...) for i in 1:ngens(M)]
   return M(v)
end

function (M::FreeModule{T})(a::Vector{T}) where T <: Union{RingElement, NCRingElem}
   return free_module_elem{T}(M, a)
end

function FreeModule(R::NCRing, rank::Int)
   T = elem_type(R)
   return FreeModule{T}(R, rank)
end

