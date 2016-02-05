#################################################################################
#
#   Factorisation
#
#################################################################################

function pari_factor_to_dict(f::PariFactor, R::Ring)
  D = Dict{typeof(zero(R)), Int}()
  for i=1:f.len
    p, n = f[i]
    D[R(p)] = n
  end
  return D
end

function factor(n::fmpz)
   f = factor(pari(n))
   return pari_factor_to_dict(f, FlintZZ)
end

function factor(g::fmpz_poly)
   h = pari(g)
   f = factor(h)
   return pari_factor_to_dict(f, g.parent)
end

function factor(g::fmpq_poly)
   h = pari(g)
   f = factor(h)
   return pari_factor_to_dict(f, g.parent)
end

