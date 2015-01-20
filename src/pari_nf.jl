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
         nf = ccall((:nfinit, :libpari), Ptr{Int}, (Ptr{Int}, Int, Int), pari(pol).d, 0, 4)
         r = PariMaximalOrderID[pol] = new(nf, pol)
         return r
      end
   end
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

###########################################################################################
#
#   MaximalOrder constructor
#
###########################################################################################

function MaximalOrder(nf::NfNumberField)
   return PariMaximalOrder(nf.pol)
end

