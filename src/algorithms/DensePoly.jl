###############################################################################
#
#   DensePoly.jl : Generic algorithms for dense univariate arithmetic
#                  Currently only Karatsuba multiplication is implemented
#
###############################################################################

module DensePoly

import AbstractAlgebra: Ring, elem_type, mul_red!, reduce!, addeq!

# usual polynomial mullow with the "fast" techniques depending on cutoffs
# write the lower Alen coefficients of
#   (B[1] + B[2]*x ... B[Blen]*x^(Blen - 1))*
#   (C[1] + C[2]*x ... C[Clen]*x^(Clen - 1))
function mullow_fast!(
   A::Vector{T}, Alen::Int,
   B::Vector{T}, Blen::Int,
   C::Vector{T}, Clen::Int,
   R::Ring,
   cutoff::Int) where T

   macc(0, A, 1, Alen, B, 1, Blen, C, 1, Clen, R, cutoff)
end

# if dr = 0, write the Alen lower coeffs of the product B*C
#   to A[Aoff + 0], ..., A[Aoff + Alen - 1]
# if dr > 0, do A += B*C, if dr < 0, do A -= B*C
function macc(
   dr::Int,
   A::Vector{T}, Aoff::Int, Alen::Int,
   B::Vector{T}, Boff::Int, Blen::Int,
   C::Vector{T}, Coff::Int, Clen::Int,
   R::Ring,
   cutoff::Int) where T

   @assert elem_type(R) == T
   @assert 0 < Aoff && Aoff + Alen - 1 <= length(A)
   @assert 0 < Boff && Boff + Blen - 1 <= length(B)
   @assert 0 < Coff && Coff + Clen - 1 <= length(C)

   if Blen < 3 || Clen < 3 || Alen < 4
      macc_classical(dr, A, Aoff, Alen, B, Boff, Blen, C, Coff, Clen, R)
   elseif Blen >= 2*Clen
      macc_unbalanced(dr, A, Aoff, Alen, B, Boff, Blen, C, Coff, Clen, R, cutoff)
   elseif Clen >= 2*Blen
      macc_unbalanced(dr, A, Aoff, Alen, C, Coff, Clen, B, Boff, Blen, R, cutoff)
   elseif Blen >= Clen
      macc_balanced(dr, A, Aoff, Alen, B, Boff, Blen, C, Coff, Clen, R, cutoff)
   else
      macc_balanced(dr, A, Aoff, Alen, C, Coff, Clen, B, Boff, Blen, R, cutoff)
   end
end

function macc_classical(
   dr::Int,
   A::Vector{T}, Aoff::Int, Alen::Int,
   B::Vector{T}, Boff::Int, Blen::Int,
   C::Vector{T}, Coff::Int, Clen::Int,
   R::Ring) where T

   @assert elem_type(R) == T

   t = zero(R)
   for i in 0 : Alen - 1
      jstart = max(0, i - (Clen - 1))
      jstop = min(i, Blen - 1)
      if dr > 0
         for j in jstart : jstop
            t = mul_red!(t, B[Boff + j], C[Coff + i - j], false)
            A[Aoff + i] = addeq!(A[Aoff + i], t)
         end
      elseif dr < 0
         z = zero(R)
         for j in jstart : jstop
            t = mul_red!(t, B[Boff + j], C[Coff + i - j], false)
            z = addeq!(z, t)
         end
         A[Aoff + i] -= z
      else
         A[Aoff + i] = zero(R)
         j = jstart
         if j <= jstop
            A[Aoff + i] = mul_red!(A[Aoff + i], B[Boff + j], C[Coff + i - j], false)
         end
         while (j = j + 1) <= jstop
            t = mul_red!(t, B[Boff + j], C[Coff + i - j], false)
            A[Aoff + i] = addeq!(A[Aoff + i], t)
         end
      end
      reduce!(A[Aoff + i])
   end
end

# break B in chunks of size n = Clen
function macc_unbalanced(
   dr::Int,
   A::Vector{T}, Aoff::Int, Alen::Int,
   B::Vector{T}, Boff::Int, Blen::Int,
   C::Vector{T}, Coff::Int, Clen::Int,
   R::Ring,
   cutoff::Int) where T

   @assert elem_type(R) == T
   @assert Blen > Clen > 1

   n = Clen
   if iszero(dr)
      for i in 0 : Alen - 1
         A[Aoff + i] = zero(R)
      end
      dr = 1
   end
   i = 0
   while (l = min(Alen - n*i, 2*n - 1)) > 0
      if (i + 1)*n < Blen
         macc(dr, A, Aoff + n*i, l, B, Boff + n*i, n,
                                          C, Coff, n, R, cutoff)
      else
         macc(dr, A, Aoff + n*i, l, B, Boff + n*i, Blen - n*i,
                                          C, Coff, n, R, cutoff)
         return
      end
      i += 1
   end
end

function macc_balanced(
   dr::Int,
   A::Vector{T}, Aoff::Int, Alen::Int,
   B::Vector{T}, Boff::Int, Blen::Int,
   C::Vector{T}, Coff::Int, Clen::Int,
   R::Ring,
   cutoff::Int) where T

   @assert elem_type(R) == T
   @assert Alen > 2
   @assert Blen >= Clen > 1

   if Alen <= Blen
      # divide trapezoid into square + triangle + trapezoid
      n = div(Alen + 1, 2)
      m = Alen - n + 1
      if n >= Clen
         macc_unbalanced(dr, A, Aoff, Alen, B, Boff, Blen,
                                                  C, Coff, Clen, R, cutoff)
         return
      end
      macc(dr, A, Aoff, Alen, B, Boff, m, C, Coff, n, R, cutoff)
      dr += iszero(dr)
      macc(dr, A, Aoff + m, n - 1, B, Boff + m, n - 1,
                                         C, Coff, n - 1, R, cutoff)
      macc(dr, A, Aoff + n, m - 1, B, Boff, m - 1,
                                  C, Coff + n, min(Clen - n, m - 1), R, cutoff)
   elseif Clen < cutoff
      macc_classical(dr, A, Aoff, Alen, B, Boff, Blen, C, Coff, Clen, R)
   else
      n = div(Clen + 1, 2)
      macc_toom22(dr, A, Aoff, Alen, B, Boff, Blen,
                                           C, Coff, Clen, n, R, cutoff)
   end
end

# write the Alen lower coeffs of the sum B + C
# to A[Aoff + 0], ..., A[Aoff + Alen - 1]
# note! the resulting A might have entries that are identical to those of B or C
function add(
   A::Vector{T}, Aoff::Int, Alen::Int,
   B::Vector{T}, Boff::Int, Blen::Int,
   C::Vector{T}, Coff::Int, Clen::Int,
   R::Ring) where T

   @assert elem_type(R) == T

   for i in 0 : Alen - 1
      if i < Blen
         if i < Clen
            A[Aoff + i] = B[Boff + i] + C[Coff + i]
         else
            A[Aoff + i] = B[Boff + i]
         end
      else
         if i < Clen
            A[Aoff + i] = C[Coff + i]
         else
            A[Aoff + i] = zero(R)
         end
      end
   end
end

function acc(
   dr::Int,
   A::Vector{T}, Aoff::Int,
   B::Vector{T}, Boff::Int, Blen::Int) where T

   @assert !iszero(dr)

   if dr < 0
      for i in 0 : Blen - 1
         A[Aoff + i] -= B[Boff + i]
      end
   else
      for i in 0 : Blen - 1
         A[Aoff + i] = addeq!(A[Aoff + i], B[Boff + i])
      end
   end
end

# (b1*x^n + b0)*(c1*x^n + c0) = p*x^2n + (q - p - r)*x^n + r
# p = b1*c1
# q = (b1 + b0)*(c1 + c0) = u*v
# r = b0*c0
function macc_toom22(
   dr::Int,
   A::Vector{T}, Aoff::Int, Alen::Int,
   B::Vector{T}, Boff::Int, Blen::Int,
   C::Vector{T}, Coff::Int, Clen::Int,
   n::Int,
   R::Ring,
   cutoff::Int) where T

   @assert elem_type(R) == T
   @assert 0 < n < min(Blen, Clen)

   ulen = max(Blen - n, n)
   vlen = max(Clen - n, n)
   rlen = max(1, min(Alen, n + n - 1))
   plen = max(1, min(Alen - n, Blen - n + Clen - n - 1))

   # temp space can be shared a little
   u = r = Vector{T}(undef, max(ulen, rlen))
   v = p = Vector{T}(undef, max(vlen, plen))

   macc(0, r, 1, rlen, B, Boff, n, C, Coff, n, R, cutoff)
   macc(0, p, 1, plen, B, Boff + n, Blen - n, C, Coff + n, Clen - n, R, cutoff)

   if iszero(dr)
      for i in 0 : Alen - 1
         if i < rlen
            A[Aoff + i] = deepcopy(r[1 + i])
         elseif 0 <= i - 2*n < plen
            A[Aoff + i] = deepcopy(p[1 + i - 2*n])
         else
            A[Aoff + i] = zero(R)
         end
      end
      dr = 1
   else
      acc(dr, A, Aoff, r, 1, min(rlen, Alen))
      acc(dr, A, Aoff + 2*n, p, 1, min(plen, Alen - 2*n))
   end

   acc(-dr, A, Aoff + n, p, 1, min(plen, Alen - n))
   acc(-dr, A, Aoff + n, r, 1, min(rlen, Alen - n))

   add(u, 1, ulen, B, Boff + n, Blen - n, B, Boff, n, R)
   add(v, 1, vlen, C, Coff + n, Clen - n, C, Coff, n, R)
   macc(dr, A, Aoff + n, Alen - n, u, 1, ulen, v, 1, vlen, R, cutoff)
end

end
