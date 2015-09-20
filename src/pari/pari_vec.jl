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
   _checkbounds(length(a), n) || throw(BoundsError())
   a.parent.base_ring(reinterpret(Ptr{Int}, 
                      unsafe_load(a.data + n*sizeof(Int))))
end

length(vec::pari_vec) = len = lg(vec.data) - 1
   
###############################################################################
#
#   AbstractString{} I/O
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
