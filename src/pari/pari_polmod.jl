###########################################################################################
#
#   pari_polmod.jl : Residue rings for polynomials in Pari
#
###########################################################################################

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function residue{S <: PariRing}(a::pari_polmod{S})
   return pari_poly{S}(reinterpret(Ptr{Int}, unsafe_load(a.data, 3)))
end

function modulus{S <: PariRing}(a::pari_polmod{S})
   return pari_poly{S}(reinterpret(Ptr{Int}, unsafe_load(a.data, 2)))
end

###########################################################################################
#
#   Parent object call overloads
#
###########################################################################################

function Base.call{S <: PariRing}(ord::PariPolModRing{S}, b::Ptr{Int})
   z = pari_polmod{S}(b)
   z.parent = ord
   return z
end
