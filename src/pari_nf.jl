###########################################################################################
#
#   pari_nf.jl : Pari nf objects
#
###########################################################################################

export PariNumberField

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

PariNumberFieldID = Dict{fmpq_poly, PariRing}()

type PariNumberField{S, T} <: PariRing
   data::Ptr{Int}
   nf::NfNumberField
   
   function PariNumberField(nf::NfNumberField)
      try
         return PariNumberFieldID[nf.pol]
      catch
         av = unsafe_load(avma, 1)
         p = pari(nf.pol)
         d = gclone(ccall((:nfinit, :libpari), Ptr{Int}, 
                           (Ptr{Int}, Int), p.d, 5))
         unsafe_store!(avma, av, 1)
         ord = new(d, nf)
         finalizer(ord, _pari_nf_unclone)
         return PariNumberFieldID[nf.pol] = ord
      end
   end
end

_pari_nf_unclone(a::PariNumberField) = gunclone(a.data)

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function pol{S, T}(nf::PariNumberField{S, T})
   data = reinterpret(Ptr{Int}, unsafe_load(nf.data + sizeof(Int)))
   return pari_poly{PariIntegerRing, T}(data)
end

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, nf::PariNumberField)
   print(io, "Maximal order of Number Field over Rational Field with defining polynomial ")
   print(io, nf.nf.pol)
end




