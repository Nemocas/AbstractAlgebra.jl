###############################################################################
#
#   Mat space
#
###############################################################################

struct MatSpace{T <: NCRingElement} <: Module{T}
  base_ring::NCRing
  nrows::Int
  ncols::Int

  function MatSpace{T}(R::NCRing, r::Int, c::Int, cached::Bool = true) where T <: NCRingElement
     # TODO/FIXME: `cached` is ignored and only exists for backwards compatibility
     @assert isconcretetype(T)
     @assert elem_type(R) === T
     (r < 0 || c < 0) && error("Dimensions must be non-negative")
     return new{T}(R, r, c)
  end
end
