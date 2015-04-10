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

function getindex(a::Factor{IntegerRing}, i::Int)
   p, n = a.d[i]
   return ZZ(p), n
end

function getindex{S}(a::Factor{FmpzPolyRing{S}}, i::Int)
   p, n = a.d[i]
   return a.parent(p), n
end

function getindex{S}(a::Factor{FmpqPolyRing{S}}, i::Int)
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

function factor(n::BigInt)
   f = factor(pari(n))
   return Factor{IntegerRing}(f, f.len, ZZ)
end

function factor{S}(g::fmpz_poly{S})
   f = factor(pari(g))
   return Factor{FmpzPolyRing{S}}(f, f.len, g.parent)
end
function factor{S}(g::fmpq_poly{S})
   f = factor(pari(g))
   return Factor{FmpqPolyRing{S}}(f, f.len, g.parent)
end

