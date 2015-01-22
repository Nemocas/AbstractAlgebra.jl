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

type PariPolModRing{S <: PariRing, T} <: PariRing

   function PariPolModRing()
      try
         return PariPolModID[S, T]
      catch
         return PariPolModID[S, T] = new()
      end
   end
end

type pari_polmod{S <: PariRing, T} <: RingElem
   data::Ptr{Int}
   parent::PariPolModRing{S, T}

   function pari_polmod(data::Ptr{Int})
      r = new(gclone(data), PariPolModRing{S, T}())
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

function residue{S <: PariRing, T}(a::pari_polmod{S, T})
   return pari_poly{S, T}(reinterpret(Ptr{Int}, unsafe_load(a.data, 3)))
end

function modulus{S <: PariRing, T}(a::pari_polmod{S, T})
   return pari_poly{S, T}(reinterpret(Ptr{Int}, unsafe_load(a.data, 2)))
end

###########################################################################################
#
#   Parent object call overloads
#
###########################################################################################

function Base.call{S <: PariRing, T}(ord::PariPolModRing{S, T}, b::Ptr{Int})
   return pari_polmod{S, T}(b)
end
