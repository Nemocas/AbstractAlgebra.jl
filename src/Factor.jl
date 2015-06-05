#################################################################################
#
#   Factorisation
#
#################################################################################

type Factor{T <: Ring}
   d::PariFactor
   len::Int
   parent::T
end

function getindex(a::Factor{FlintIntegerRing}, i::Int)
   p, n = a.d[i]
   return ZZ(p), n
end

function getindex(a::Factor{FmpzPolyRing}, i::Int)
   p, n = a.d[i]
   return a.parent(p), n
end

function getindex(a::Factor{FmpqPolyRing}, i::Int)
   p, n = a.d[i]
   return a.parent(p), n
end

function show(io::IO, a::Factor)
   print(io, "[")
   for i = 1:a.len
      print(io, a[i])
      if i != a.len
         print(io, ", ")
      end
   end
   print(io, "]")
end

function factor(n::fmpz)
   f = factor(pari(n))
   return Factor{FlintIntegerRing}(f, f.len, ZZ)
end

function factor(g::fmpz_poly)
   h = pari(g)
   f = factor(h)
   return Factor{FmpzPolyRing}(f, f.len, g.parent)
end

function factor(g::fmpq_poly)
   h = pari(g)
   f = factor(h)
   return Factor{FmpqPolyRing}(f, f.len, g.parent)
end

