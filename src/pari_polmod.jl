###########################################################################################
#
#   pari_polmod.jl : Maximal orders via Pari nf objects
#
###########################################################################################

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

PariPolModID = Dict{(DataType, Symbol), PariRing}()

type PariPolModRing{S <: PariRing} <: PariRing
   T::Symbol

   function PariPolModRing(t::Symbol)
      try
         return PariPolModID[S, t]
      catch
         return PariPolModID[S, t] = new(t)
      end
   end
end

type pari_polmod{S <: PariRing} <: RingElem
   data::Ptr{Int}
   parent::PariPolModRing{S}

   function pari_polmod(data::Ptr{Int})
      r = new(gclone(data))
      finalizer(r, _pari_polmod_unclone)
      return r
   end
end

_pari_polmod_unclone(a::pari_polmod) = gunclone(a.data)

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
