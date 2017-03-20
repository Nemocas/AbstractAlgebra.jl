###############################################################################
#
#   pari_poly2.jl : functions for Pari t_POL objects over fields
#
###############################################################################

export fmpq_poly!

###############################################################################
#
#   Conversions to/from Poly{fmpq}
#
###############################################################################

function gensize(a::fmpq_poly)
   coeffs = a.coeffs
   total = 0
   for i = 0:a.length - 1
      total += ccall((:fmpz_size, :libflint), Int, 
                     (Ptr{Int},), coeffs + i*sizeof(Int))
   end
   densize = ccall((:fmpz_size, :libflint), Int, (Ptr{Int},), &a.den)
   return total + densize*a.length + 8*a.length + 2
end

function pari!(x::Ptr{Int}, a::fmpq_poly)
   unsafe_store!(x, evaltyp(t_POL) | (a.length + 2), 1) 
   unsafe_store!(x, evalsigne(Int(a.length != 0)) | evalvarn(0), 2)
   s = a.length + 2
   for i = 0:a.length - 1
      z = coeff(a, i)
      if z == 0
         unsafe_store!(x, unsafe_load(gen_0), i + 3)
      else
         unsafe_store!(x, x + sizeof(Int)*s, i + 3)
         s += pari!(x + sizeof(Int)*s, z)
      end
   end
   return s
end

function pari(a::fmpq_poly)
   s = gensize(a)
   g = pari_poly{pari_rat}(s)
   g.parent = PariPolyRing{pari_rat}(PariQQ, var(parent(a)))
   pari!(g.d, a)
   return g
end

function fmpq_poly!(z::fmpq_poly, g::Ptr{Int})
   length = (unsafe_load(g, 1) & LGBITS) - 2
   fit!(z, length)
   z.length = length
   if length == 0
      return
   end
   c = fmpq()
   for i = 0:length - 1
      c_ptr = unsafe_load(g, 3 + i)
      if reinterpret(Ptr{Int}, c_ptr) == reinterpret(Ptr{Int}, unsafe_load(gen_0))
         c = zero(FlintQQ)
      else
         fmpq!(c, reinterpret(Ptr{Int}, c_ptr))
      end
      setcoeff!(z, i, c)
   end
   return z
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (ord::PariPolyRing{pari_rat})(n::Ptr{Int})
   pol = pari_poly{pari_rat}(n)
   pol.parent = PariPolyRing{pari_rat}(PariQQ, var(ord))
   return pol
end

function (a::FmpqPolyRing)(g::pari_poly{pari_rat})
   z = fmpq_poly()
   z.parent = a
   fmpq_poly!(z, g.d)
   return z
end
