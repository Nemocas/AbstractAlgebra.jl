###############################################################################
#
#   pari_vec.jl : Pari vectors over RingElems
#
###############################################################################

export PariVector, pari_vec

###############################################################################
#
#   Array [] overloading
#
###############################################################################

function getindex(a::pari_vec, n::Int)
   a.parent.base_ring(reinterpret(Ptr{Int}, 
                      unsafe_load(a.data + n*sizeof(Int))))
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, vec::pari_vec)
   print(io, "[")
   len = lg(vec.data) - 1
   for i = 1:len
      print(io, vec.parent.base_ring(reinterpret(Ptr{Int}, 
                                     unsafe_load(vec.data, i + 1))))
      if i != len
         print(io, ", ")
      end
   end
   print(io, "]")
end
