###########################################################################################
#
#   pari_nf.jl : Maximal orders via Pari nf objects
#
###########################################################################################

export MaximalOrder

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

PariMaximalOrderID = Dict{fmpq_poly, Ring}()

type PariMaximalOrder <: Ring
   nf::Ptr{Int}
   pol::fmpq_poly

   function PariMaximalOrder(pol::fmpq_poly)
      try
         return PariMaximalOrderID[pol]
      catch
         av = unsafe_load(avma, 1)
         p = pari(pol)
         nf = gclone(ccall((:nfinit, :libpari), Ptr{Int}, 
                           (Ptr{Int}, Int), p.d, 5))
         unsafe_store!(avma, av, 1)
         ord = new(nf, pol)
         finalizer(ord, _pari_nf_clear_fn)
         return PariMaximalOrderID[pol] = ord
      end
   end
end

function _pari_nf_clear_fn(a::PariMaximalOrder)
   ccall((:gunclone, :libpari), Void, (Ptr{Int},), a.nf)
end

type PariIdeal <: RingElem
   ideal::Ptr{Int}
   parent::PariMaximalOrder

   function PariIdeal(a::Ptr{Int})
      r = new(a)
      finaliser(r, _pari_ideal_clear_fn)
      return r
   end
end

function _pari_ideal_clear_fn(a::PariIdeal)
   ccall((:gunclone, :libpari), Void, (Ptr{Int},), a.ideal)
end

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, nf::PariMaximalOrder)
   print(io, "Maximal order of Number Field over Rational Field with defining polynomial ")
   print(nf.pol, " with Z-basis ")
   cstr = ccall((:GENtostr, :libpari), Ptr{Uint8}, 
                (Ptr{Int},), reinterpret(Ptr{Int}, unsafe_load(nf.nf + 7*sizeof(Int))))

   print(io, bytestring(cstr))

   ccall((:pari_free, :libpari), Void, (Ptr{Uint8},), cstr)
end

function show(io::IO, id::PariIdeal)
   av = unsafe_load(avma, 1)
   id2 = ccall((:idealtwoelt, :libpari), Ptr{Int}, 
               (Ptr{Int}, Ptr{Int}), id.parent.nf, id.ideal)
   cstr = ccall((:GENtostr, :libpari), Ptr{Uint8}, 
                (Ptr{Int},), id2)
   unsafe_store!(avma, av, 1)
   
   print(io, bytestring(cstr))

   ccall((:pari_free, :libpari), Void, (Ptr{Uint8},), cstr)
end

###########################################################################################
#
#   Binary operations
#
###########################################################################################

function *(a::PariIdeal, b::PariIdeal)
   av = unsafe_load(avma, 1)
   r = gclone(ccall((:idealmul, :libpari), Ptr{Int}, 
                    (Ptr{Int}, Ptr{Int}, Ptr{Int}), a.parent.nf, a.ideal, b.ideal))
   unsafe_store!(avma, av, 1)
   return PariIdeal(r, a.parent)
end

function divexact(a::PariIdeal, b::PariIdeal)
   av = unsafe_load(avma, 1)
   r = gclone(ccall((:idealdiv, :libpari), Ptr{Int}, 
                    (Ptr{Int}, Ptr{Int}, Ptr{Int}), a.parent.nf, a.ideal, b.ideal))
   unsafe_store!(avma, av, 1)
   return PariIdeal(r, a.parent)
end

###########################################################################################
#
#   Parent call object overloads
#
###########################################################################################

function Base.call(ord::PariMaximalOrder, b::fmpq_poly)
   av = unsafe_load(avma, 1)
   ideal = gclone(ccall((:idealhnf, :libpari), Ptr{Int}, 
                 (Ptr{Int}, Ptr{Int}), ord.nf, pari(b).d))
   unsafe_store!(avma, av, 1)
   return PariIdeal(ideal, ord)
end

###########################################################################################
#
#   MaximalOrder constructor
#
###########################################################################################

function MaximalOrder(nf::NfNumberField)
   return PariMaximalOrder(nf.pol)
end

