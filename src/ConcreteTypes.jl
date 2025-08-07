###############################################################################
#
#   Parent types
#
###############################################################################

# parent type for matrices
struct MatSpace{T <: NCRingElement} <: Module{T}
  base_ring::NCRing
  nrows::Int
  ncols::Int

  function MatSpace{T}(R::NCRing, r::Int, c::Int, cached::Bool = true) where T <: NCRingElement
     # TODO/FIXME: `cached` is ignored and only exists for backwards compatibility
     @assert elem_type(R) === T
     (r < 0 || c < 0) && error("Dimensions must be non-negative")
     return new{T}(R, r, c)
  end
end

# parent type for two-sided ideals
struct DefaultIdealSet{T <: NCRingElement} <: IdealSet{T}
   base_ring::NCRing

   DefaultIdealSet(R::S) where {S <: NCRing}  = new{elem_type(S)}(R)
end
