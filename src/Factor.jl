#################################################################################
#
#   Factorisation
#
#################################################################################

function _Factor(f::PariFactor, R::Ring)
  D = Dict{typeof(zero(R)), Int}()
  for i=1:f.len
    p, n = f[i]
    D[R(p)] = n
  end
  return D
end

function factor(n::fmpz)
   f = factor(pari(n))
   return _Factor(f, FlintZZ)
end

function factor(g::fmpz_poly)
   h = pari(g)
   f = factor(h)
   return _Factor(f, g.parent)
end

function factor(g::fmpq_poly)
   h = pari(g)
   f = factor(h)
   return _Factor(f, g.parent)
end

