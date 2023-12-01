"""
Provides generic asymptotically fast matrix methods:
  - mul and mul! using the Strassen scheme
  - solve_tril!
  - lu!
  - solve_triu

Just prefix the function by "Strassen." all 4 functions support a keyword
argument "cutoff" to indicate when the base case should be used.

The speedup depends on the ring and the entry sizes.

#Examples:

```jldoctest; setup = :(using AbstractAlgebra)
julia> m = matrix(ZZ, rand(-10:10, 1000, 1000));

julia> n = similar(m);

julia> mul!(n, m, m);

julia> Strassen.mul!(n, m, m);

julia> Strassen.mul!(n, m, m; cutoff = 100);

```
"""
module Strassen
using AbstractAlgebra
import AbstractAlgebra:Perm

const cutoff = 1500

function mul(A::MatElem{T}, B::MatElem{T}; cutoff::Int = cutoff) where {T}
  C = zero_matrix(base_ring(A), nrows(A), ncols(B))
  mul!(C, A, B; cutoff)
  return C
end

#scheduling copied from the nmod_mat_mul in Flint
function mul!(C::MatElem{T}, A::MatElem{T}, B::MatElem{T}; cutoff::Int = cutoff) where {T}
  sA = size(A)
  sB = size(B)
  sC = size(C)
  a = sA[1]
  b = sA[2]
  c = sB[2]

  @assert a == sC[1] && b == sB[1] && c == sC[2]

  if (a <= cutoff || b <= cutoff || c <= cutoff)
      AbstractAlgebra.mul!(C, A, B)
      return
  end

  anr = div(a, 2)
  anc = div(b, 2)
  bnr = anc
  bnc = div(c, 2)

  #nmod_mat_window_init(A11, A, 0, 0, anr, anc);
  #nmod_mat_window_init(A12, A, 0, anc, anr, 2*anc);
  #nmod_mat_window_init(A21, A, anr, 0, 2*anr, anc);
  #nmod_mat_window_init(A22, A, anr, anc, 2*anr, 2*anc);
  A11 = view(A, 1:anr, 1:anc)
  A12 = view(A, 1:anr, anc+1:2*anc)
  A21 = view(A, anr+1:2*anr, 1:anc)
  A22 = view(A, anr+1:2*anr, anc+1:2*anc)

  #nmod_mat_window_init(B11, B, 0, 0, bnr, bnc);
  #nmod_mat_window_init(B12, B, 0, bnc, bnr, 2*bnc);
  #nmod_mat_window_init(B21, B, bnr, 0, 2*bnr, bnc);
  #nmod_mat_window_init(B22, B, bnr, bnc, 2*bnr, 2*bnc);
  B11 = view(B, 1:bnr, 1:bnc)
  B12 = view(B, 1:bnr, bnc+1:2*bnc)
  B21 = view(B, bnr+1:2*bnr, 1:bnc)
  B22 = view(B, bnr+1:2*bnr, bnc+1:2*bnc)

  #nmod_mat_window_init(C11, C, 0, 0, anr, bnc);
  #nmod_mat_window_init(C12, C, 0, bnc, anr, 2*bnc);
  #nmod_mat_window_init(C21, C, anr, 0, 2*anr, bnc);
  #nmod_mat_window_init(C22, C, anr, bnc, 2*anr, 2*bnc);
  C11 = view(C, 1:anr, 1:bnc)
  C12 = view(C, 1:anr, bnc+1:2*bnc)
  C21 = view(C, anr+1:2*anr, 1:bnc)
  C22 = view(C, anr+1:2*anr, bnc+1:2*bnc)

  #nmod_mat_init(X1, anr, FLINT_MAX(bnc, anc), A->mod.n);
  #nmod_mat_init(X2, anc, bnc, A->mod.n);

  #X1->c = anc;

  #=
      See Jean-Guillaume Dumas, Clement Pernet, Wei Zhou; "Memory
      efficient scheduling of Strassen-Winograd's matrix multiplication
      algorithm"; http://arxiv.org/pdf/0707.2347v3 for reference on the
      used operation scheduling.
  =#

  X1 = A11 - A21
  X2 = B22 - B12
  #nmod_mat_mul(C21, X1, X2);
  mul!(C21, X1, X2; cutoff)

  add!(X1, A21, A22);
  sub!(X2, B12, B11);
  #nmod_mat_mul(C22, X1, X2);
  mul!(C22, X1, X2; cutoff)

  sub!(X1, X1, A11);
  sub!(X2, B22, X2);
  #nmod_mat_mul(C12, X1, X2);
  mul!(C12, X1, X2; cutoff)

  sub!(X1, A12, X1);
  #nmod_mat_mul(C11, X1, B22);
  mul!(C11, X1, B22; cutoff)

  #X1->c = bnc;
  #nmod_mat_mul(X1, A11, B11);
  mul!(X1, A11, B11; cutoff)

  add!(C12, X1, C12);
  add!(C21, C12, C21);
  add!(C12, C12, C22);
  add!(C22, C21, C22);
  add!(C12, C12, C11);
  sub!(X2, X2, B21);
  #nmod_mat_mul(C11, A22, X2);
  mul!(C11, A22, X2; cutoff)

  sub!(C21, C21, C11);

  #nmod_mat_mul(C11, A12, B21);
  mul!(C11, A12, B21; cutoff)

  add!(C11, X1, C11);

  if c > 2*bnc #A by last col of B -> last col of C
      #nmod_mat_window_init(Bc, B, 0, 2*bnc, b, c);
      Bc = view(B, 1:b, 2*bnc+1:c)
      #nmod_mat_window_init(Cc, C, 0, 2*bnc, a, c);
      Cc = view(C, 1:a, 2*bnc+1:c)
      #nmod_mat_mul(Cc, A, Bc);
      AbstractAlgebra.mul!(Cc, A, Bc)
  end

  if a > 2*anr #last row of A by B -> last row of C
      #nmod_mat_window_init(Ar, A, 2*anr, 0, a, b);
      Ar = view(A, 2*anr+1:a, 1:b)
      #nmod_mat_window_init(Cr, C, 2*anr, 0, a, c);
      Cr = view(C, 2*anr+1:a, 1:c)
      #nmod_mat_mul(Cr, Ar, B);
      AbstractAlgebra.mul!(Cr, Ar, B)
  end

  if b > 2*anc # last col of A by last row of B -> C
      #nmod_mat_window_init(Ac, A, 0, 2*anc, 2*anr, b);
      Ac = view(A, 1:2*anr, 2*anc+1:b)
      #nmod_mat_window_init(Br, B, 2*bnr, 0, b, 2*bnc);
      Br = view(B, 2*bnr+1:b, 1:2*bnc)
      #nmod_mat_window_init(Cb, C, 0, 0, 2*anr, 2*bnc);
      Cb = view(C, 1:2*anr, 1:2*bnc)
      #nmod_mat_addmul(Cb, Cb, Ac, Br);
      AbstractAlgebra.mul!(Cb, Ac, Br, true)
  end
end

#solve_tril fast, recursive
# A   *  X  = U
# B C    Y    V
# => X = solve(A, U)
# Y = solve(C, V - B*X)
function solve_tril!(A::MatElem{T}, B::MatElem{T}, C::MatElem{T}, f::Int = 0; cutoff::Int = 2) where T
  if nrows(A) < cutoff || ncols(A) < cutoff
    return AbstractAlgebra.solve_tril!(A, B, C, f)
  end
  n = nrows(B)
  n2 = div(n, 2)
  B11 = view(B, 1:n2, 1:n2)
  B21 = view(B, n2+1:n, 1:n2)
  B22 = view(B, n2+1:n, n2+1:ncols(B))
  X1 = view(A, 1:n2, 1:ncols(A))
  X2 = view(A, n2+1:n, 1:ncols(A))
  C1 = view(C, 1:n2, 1:ncols(A))
  C2 = view(C, n2+1:n, 1:ncols(A))
  solve_tril!(X1, B11, C1, f; cutoff)
  x = B21 * X1  # strassen...
  sub!(X2, C2, x)
  solve_tril!(X2, B22, X2, f; cutoff)
end

function apply!(A::MatElem, P::Perm{Int}; offset::Int = 0)
  n = length(P.d)
  Q = copy(inv(P).d) #the inv is experimentally verified with the other apply
  cnt = 0
  start  = 0
  while cnt < n
    ptr = start = findnext(!iszero, Q, start +1)::Int
    next = Q[start]
    cnt += 1
    while next != start
      swap_rows!(A, ptr, next)
      Q[ptr] = 0
      next = Q[next]
      cnt += 1
    end
    Q[ptr] = 0
  end
end

function apply!(Q::Perm{Int}, P::Perm{Int}; offset::Int = 0)
  n = length(P.d)
  t = zeros(Int, n-offset)
  for i=1:n-offset
    t[i] = Q.d[P.d[i] + offset]
  end
  for i=1:n-offset
    Q.d[i + offset] = t[i]
  end
end

function lu!(P::Perm{Int}, A; cutoff::Int = 300)
  m = nrows(A)

  @assert length(P.d) == m
  n = ncols(A)
  if n < cutoff
    return AbstractAlgebra.lu!(P, A)
  end
  n1 = div(n, 2)
  for i=1:m
    P.d[i] = i
  end
  P1 = AbstractAlgebra.Perm(m)
  A0 = view(A, 1:m, 1:n1)
  r1 = lu!(P1, A0; cutoff)
  @assert r1 == n1
  if r1 > 0 
    apply!(A, P1)
    apply!(P, P1)
  end

  A00 = view(A, 1:r1, 1:r1)
  A10 = view(A, r1+1:m, 1:r1)
  A01 = view(A, 1:r1, n1+1:n)
  A11 = view(A, r1+1:m, n1+1:n)

  if r1 > 0
    #Note: A00 is a view of A0 thus a view of A
    # A0 is lu!, thus implicitly two triangular matrices giving the
    # lu decomosition. solve_tril! looks ONLY at the lower part of A00
    solve_tril!(A01, A00, A01, 1)
    X = A10 * A01
    sub!(A11, A11, X)
  end

  P1 = Perm(nrows(A11))
  r2 = lu!(P1, A11)
  apply!(A, P1, offset = r1)
  apply!(P, P1, offset = r1)

  if (r1 != n1)
    for i=1:m-r1
      for j=1:min(i, r2)
        A[r1+i-1, r1+j-1] = A[r1+i-1, n1+j-1]
        A[r1+i-1, n1+j-1] = 0
      end
    end
  end
  return r1 + r2
end

function solve_triu(T::MatElem, b::MatElem; cutoff::Int = cutoff)
  #b*inv(T), thus solves Tx = b for T upper triangular
  n = ncols(T)
  if n <= cutoff
    R = AbstractAlgebra.solve_triu(T, b)
    return R
  end

  n2 = div(n, 2) + n % 2
  m = nrows(b)
  m2 = div(m, 2) + m % 2

  U = view(b, 1:m2, 1:n2)
  V = view(b, 1:m2, n2+1:n)
  X = view(b, m2+1:m, 1:n2)
  Y = view(b, m2+1:m, n2+1:n)

  A = view(T, 1:n2, 1:n2)
  B = view(T, 1:n2, 1+n2:n)
  C = view(T, 1+n2:n, 1+n2:n)

  S = solve_triu(A, U; cutoff)
  R = solve_triu(A, X; cutoff)

  SS = mul(S, B; cutoff)
  sub!(SS, V, SS)
  SS = solve_triu(C, SS; cutoff)

  RR = mul(R, B; cutoff)
  sub!(RR, Y, RR)
  RR = solve_triu(C, RR; cutoff)

  return [S SS; R RR]
end

end # module
