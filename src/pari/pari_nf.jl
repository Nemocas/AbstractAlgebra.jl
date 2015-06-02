###############################################################################
#
#   pari_nf.jl : Pari nf objects
#
###############################################################################

export PariNumberField

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function pol(nf::PariNumberField)
   data = reinterpret(Ptr{Int}, unsafe_load(nf.data + sizeof(Int)))
   return pari_poly{pari_int}(data)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, nf::PariNumberField)
   print(io, "Number Field over Rational Field")
   print(io, " with defining polynomial ", nf.nf.pol)
end




