###############################################################################
#
#   Ideal constructor
#
###############################################################################

function Ideal(R::Ring, V::Vector)
   return Generic.Ideal(R, V)
end

function Ideal(R::Ring, v::T...) where T <: RingElement
   return Generic.Ideal(R, v...)
end

function Ideal(R::Ring)
   return Generic.Ideal(R, [])
end

###############################################################################
#
#   IdealSet constructor
#
###############################################################################

function IdealSet(R::Ring)
   return Generic.IdealSet(R)
end
