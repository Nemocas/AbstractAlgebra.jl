###############################################################################
#
#   PariFactor.jl : Pari factorization
#
###############################################################################

function getindex(a::PariFactor, i::Int)
   i > a.len && throw(IndexError())
   colf = reinterpret(Ptr{Int}, unsafe_load(a.data, 2))
   coln = reinterpret(Ptr{Int}, unsafe_load(a.data, 3))
   f = reinterpret(Ptr{Int}, unsafe_load(colf, i + 1))
   n = fmpz()
   fmpz!(n, reinterpret(Ptr{Int}, unsafe_load(coln, i + 1)))
   return a.parent(f), Int(n)
end

function show(io::IO, a::PariFactor)
   print(io, "[")
   for i = 1:a.len
      print(io, a[i])
      if i != a.len
         print(io, ", ")
      end
   end
   print(io, "]")
end

