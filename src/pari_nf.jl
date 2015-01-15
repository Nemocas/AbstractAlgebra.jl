###########################################################################################
#
#   pari_nf.jl : Maximal orders via Pari nf objects
#
###########################################################################################

export PariMaximalOrder

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

PariMaximalOrderID = Dict{fmpz_poly, Ring}()

type PariMaximalOrder <: Ring
   nf::Ptr{Int}
   pol::fmpz_poly

   function PariMaximalOrder(pol::fmpz_poly)
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
   println(io, "Maximal order of Number Field over Rational Field with defining polynomial ", nf.pol)
   cstr = ccall((:GENtostr, :libpari), Ptr{Uint8}, 
                (Ptr{Int},), nf.nf)

   print(io, bytestring(cstr))

   ccall((:pari_free, :libpari), Void, (Ptr{Uint8},), cstr)
end
