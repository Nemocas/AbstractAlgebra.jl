###############################################################################
#
#   Ideal constructor
#
###############################################################################

function Ideal(R::Ring, V::Vector{T}) where T <: RingElement
   return Generic.Ideal(R, V)
end

function Ideal(R::Ring, v::T...) where T <: RingElement
   return Generic.Ideal(R, v...)
end

###############################################################################
#
#   IdealSet constructor
#
###############################################################################

function IdealSet(R::Ring)
   return Generic.IdealSet(R)
end
