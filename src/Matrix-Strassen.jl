# FIXME: the doctest below is disabled because it takes way too long,
# especially on macOS, see https://github.com/Nemocas/AbstractAlgebra.jl/issues/2083

"""
Provides generic asymptotically fast matrix methods:
  - `mul` and `mul!` using the Strassen scheme
  - `_solve_tril!`
  - `lu!`
  - `_solve_triu`

Just prefix the function by "Strassen." all 4 functions support a keyword
argument "cutoff" to indicate when the base case should be used.

The speedup depends on the ring and the entry sizes.

# Examples

```julia
julia> m = matrix(ZZ, rand(-10:10, 1000, 1000));

julia> n1 = similar(m); n2 = similar(m); n3 = similar(m);

julia> n1 = mul!(n1, m, m);

julia> n2 = Strassen.mul!(n2, m, m);

julia> n3 = Strassen.mul!(n3, m, m; cutoff = 100);

julia> n1 == n2 == n3
true
```
"""
module Strassen
using AbstractAlgebra
import AbstractAlgebra: Perm

const cutoff = 1500

function mul(A::MatElem{T}, B::MatElem{T}; cutoff::Int = cutoff) where {T}
  C = zero_matrix(base_ring(A), nrows(A), ncols(B))
  C = mul!(C, A, B; cutoff)
  return C
end

#scheduling copied from the nmod_mat_mul in Flint
"""
Fast, recursive, generic matrix multiplication using the Strassen
trick.

`cutoff` indicates when the recursion stops and the base case is called.
"""
function mul!(C::MatElem{T}, A::MatElem{T}, B::MatElem{T}; cutoff::Int = cutoff) where {T}
  sA = size(A)
  sB = size(B)
  sC = size(C)
  a = sA[1]
  b = sA[2]
  c = sB[2]

  @assert a == sC[1] && b == sB[1] && c == sC[2]

  if (a <= cutoff || b <= cutoff || c <= cutoff)
      return AbstractAlgebra.mul!(C, A, B)
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
  C21 = mul!(C21, X1, X2; cutoff)


  X1 = add!(X1, A21, A22);
  X2 = sub!(X2, B12, B11);
  #nmod_mat_mul(C22, X1, X2);
  C22 = mul!(C22, X1, X2; cutoff)

  X1 = sub!(X1, X1, A11);
  X2 = sub!(X2, B22, X2);
  #nmod_mat_mul(C12, X1, X2);
  C12 = mul!(C12, X1, X2; cutoff)

  X1 = sub!(X1, A12, X1);
  #nmod_mat_mul(C11, X1, B22);
  C11 = mul!(C11, X1, B22; cutoff)

  #X1->c = bnc;  #this resizes X1!!!!
  #nmod_mat_mul(X1, A11, B11);

  X1 = zero_matrix(base_ring(A11), nrows(A11), ncols(B11))
  X1 = mul!(X1, A11, B11; cutoff)

  C12 = add!(C12, X1, C12);
  C21 = add!(C21, C12, C21);
  C12 = add!(C12, C12, C22);
  C22 = add!(C22, C21, C22);
  C12 = add!(C12, C12, C11);
  X2 = sub!(X2, X2, B21);
  #nmod_mat_mul(C11, A22, X2);
  C11 = mul!(C11, A22, X2; cutoff)

  C21 = sub!(C21, C21, C11);

  #nmod_mat_mul(C11, A12, B21);
  C11 = mul!(C11, A12, B21; cutoff)

  C11 = add!(C11, X1, C11);

  if c > 2*bnc #A by last col of B -> last col of C
      #nmod_mat_window_init(Bc, B, 0, 2*bnc, b, c);
      Bc = view(B, 1:b, 2*bnc+1:c)
      #nmod_mat_window_init(Cc, C, 0, 2*bnc, a, c);
      Cc = view(C, 1:a, 2*bnc+1:c)
      #nmod_mat_mul(Cc, A, Bc);
      AbstractAlgebra.mul!(Cc, A, Bc) # needs to mutate Cc
  end

  if a > 2*anr #last row of A by B -> last row of C
      #nmod_mat_window_init(Ar, A, 2*anr, 0, a, b);
      Ar = view(A, 2*anr+1:a, 1:b)
      #nmod_mat_window_init(Cr, C, 2*anr, 0, a, c);
      Cr = view(C, 2*anr+1:a, 1:c)
      #nmod_mat_mul(Cr, Ar, B);
      AbstractAlgebra.mul!(Cr, Ar, B) # needs to mutate Cr
  end

  if b > 2*anc # last col of A by last row of B -> C
      #nmod_mat_window_init(Ac, A, 0, 2*anc, 2*anr, b);
      Ac = view(A, 1:2*anr, 2*anc+1:b)
      #nmod_mat_window_init(Br, B, 2*bnr, 0, b, 2*bnc);
      Br = view(B, 2*bnr+1:b, 1:2*bnc)
      #nmod_mat_window_init(Cb, C, 0, 0, 2*anr, 2*bnc);
      Cb = view(C, 1:2*anr, 1:2*bnc)
      #nmod_mat_addmul(Cb, Cb, Ac, Br);
      AbstractAlgebra.addmul!(Cb, Ac, Br) # needs to mutate Cb
  end

  return C
end

#solve_tril fast, recursive
# A   *  X  = U
# B C    Y    V
# => X = solve(A, U)
# Y = solve(C, V - B*X)
function _solve_tril!(A::MatElem{T}, B::MatElem{T}, C::MatElem{T}, f::Int = 0; cutoff::Int = 2) where T
  if nrows(A) < cutoff || ncols(A) < cutoff
    return AbstractAlgebra._solve_tril!(A, B, C, f)
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
  _solve_tril!(X1, B11, C1, f; cutoff)
  x = B21 * X1  # strassen...
  X2 = sub!(X2, C2, x)
  _solve_tril!(X2, B22, X2, f; cutoff)
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

function lu!(P::Perm{Int}, A::MatElem; cutoff::Int = 300)
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
    # lu decomosition. _solve_tril! looks ONLY at the lower part of A00
    _solve_tril!(A01, A00, A01, 1)
    X = A10 * A01
    A11 = sub!(A11, X)
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

function _solve_triu(T::MatElem, b::MatElem; cutoff::Int = cutoff, side::Symbol = :left, unipotent::Bool = false)
  #inv(T)*b, thus solves Tx = b for T upper triangular
  n = ncols(T)
  if n <= cutoff
    R = AbstractAlgebra._solve_triu(T, b; side, unipotent)
    return R
  end
  if side == :left
    return _solve_triu_left(T, b; cutoff, unipotent)
  end
  @assert side == :right
  @assert n == nrows(T) == nrows(b)

  n2 = div(n, 2) + n % 2
  m = ncols(b)
  m2 = div(m, 2) + m % 2
  #=
    b = [U X; V Y]
    T = [A B; 0 C]
    x = [SS RR; S R]

    [0 C] [SS; S] = CS = V
    [0 C] [RR; R] = CR = Y

    [A B] [SS; S] = A SS + B S = U => A SS = U - BS
    [A B] [RR; R] = A RR + B R = U => A RR = X - BR
    
 =#   

  U = view(b, 1:n2, 1:m2)
  X = view(b, 1:n2, m2+1:m)
  V = view(b, n2+1:n, 1:m2)
  Y = view(b, n2+1:n, m2+1:m)

  A = view(T, 1:n2, 1:n2)
  B = view(T, 1:n2, 1+n2:n)
  C = view(T, 1+n2:n, 1+n2:n)

  S = _solve_triu(C, V; cutoff, side, unipotent)
  R = _solve_triu(C, Y; cutoff, side, unipotent)

  SS = mul(B, S; cutoff)
  SS = sub!(SS, U, SS)
  SS = _solve_triu(A, SS; cutoff, side, unipotent)

  RR = mul(B, R; cutoff)
  RR = sub!(RR, X, RR)
  RR = _solve_triu(A, RR; cutoff, side, unipotent)

  return [SS RR; S R]
end

function _solve_triu_left(T::MatElem, b::MatElem; cutoff::Int = cutoff, unipotent::Bool = false)
  #b*inv(T), thus solves xT = b for T upper triangular
  n = ncols(T)
  if n <= cutoff
    R = AbstractAlgebra._solve_triu_left(T, b; unipotent)
    return R
  end
  
  @assert ncols(b) == nrows(T) == n

  n2 = div(n, 2) + n % 2
  m = nrows(b)
  m2 = div(m, 2) + m % 2
  #=
    b = [U X; V Y]
    T = [A B; 0 C]
    x = [S SS; R RR]

    [S SS] [A; 0] = SA = U
    [R RR] [A; 0] = RA = V
    [S SS] [B; C] = SB + SS C = X => SS C = Y - SB
    [R RR] [B; C] = RB + RR C = Y => RR C = Y - RB

 =#   

  U = view(b, 1:m2, 1:n2)
  V = view(b, 1:m2, n2+1:n)
  X = view(b, m2+1:m, 1:n2)
  Y = view(b, m2+1:m, n2+1:n)

  A = view(T, 1:n2, 1:n2)
  B = view(T, 1:n2, 1+n2:n)
  C = view(T, 1+n2:n, 1+n2:n)

  S = _solve_triu_left(A, U; cutoff, unipotent)
  R = _solve_triu_left(A, X; cutoff, unipotent)

  SS = mul(S, B; cutoff)
  SS = sub!(SS, V, SS)
  SS = _solve_triu_left(C, SS; cutoff, unipotent)

  RR = mul(R, B; cutoff)
  RR = sub!(RR, Y, RR)
  RR = _solve_triu_left(C, RR; cutoff, unipotent)
  #THINK: both pairs of solving could be combined: 
  # solve [U; X], A to get S and R...

  return [S SS; R RR]
end

function mul_tt!(A::MatElem, B::MatElem, C::MatElem; cutoff::Int = cutoff)
  #A = BC for upper triagular square matrices A, B, C
  @assert is_upper_triangular(B)
  @assert is_upper_triangular(C)
  n = ncols(A)
  n <= cutoff && return mul!(A, B, C)

  @assert nrows(A) == n
  @assert size(A) == size(B) == size(C)
  
  #A = [a1 a2; 0 a4]
  #B = [b1 b2; 0 b4]
  #C = [c1 c2; 0 c4]
  
  n2 = div(n, 2)

  a1 = view(A, 1:n2, 1:n2)
  a2 = view(A, 1:n2, 1+n2:n)
  a4 = view(A, 1+n2:n, 1+n2:n)

  b1 = view(B, 1:n2, 1:n2)
  b2 = view(B, 1:n2, 1+n2:n)
  b4 = view(B, 1+n2:n, 1+n2:n)

  c1 = view(C, 1:n2, 1:n2)
  c2 = view(C, 1:n2, 1+n2:n)
  c4 = view(C, 1+n2:n, 1+n2:n)

  mul_tt!(a1, b1, c1; cutoff)
  mul_tu!(a2, b1, c2; cutoff)
  x = similar(a2)
  mul_ut!(x,  b2, c4; cutoff)
  add!(a2, a2, x)
  mul_tt!(a4, b4, c4; cutoff)
#  @assert A == B*C
end

function mul_tu!(A::MatElem, B::MatElem, C::MatElem; cutoff::Int = cutoff)
  #A = BC for square matrices A, B, C, B is upper triangular
  @assert is_upper_triangular(B)
  n = nrows(A)
  n <= cutoff && return mul!(A, B, C)

  m = ncols(A)
  k = ncols(B)
  @assert nrows(C) == k
  
  #A = [a1 a2; a3 a4]
  #B = [b1 b2; 0 b4]
  #C = [c1 c2; c3 c4]
  
  n2 = div(n, 2)
  m2 = div(m, 2)
  k2 = div(k, 2)

  a1 = view(A, 1:n2, 1:m2)
  a2 = view(A, 1:n2, 1+m2:m)
  a3 = view(A, 1+n2:n, 1:m2)
  a4 = view(A, 1+n2:n, 1+m2:m)

  b1 = view(B, 1:n2, 1:k2)
  b2 = view(B, 1:n2, 1+k2:k)
  b4 = view(B, 1+n2:n, 1+k2:k)

  c1 = view(C, 1:k2, 1:m2)
  c2 = view(C, 1:k2, 1+m2:m)
  c3 = view(C, 1+k2:k, 1:m2)
  c4 = view(C, 1+k2:k, 1+m2:m)

  mul_tu!(a1, b1, c1; cutoff)
  x = similar(a1)
  mul!(x, b2, c3; cutoff)
  add!(a1, a1, x)

  mul_tu!(a2, b1, c2; cutoff)
  x = similar(a2)
  mul!(x, b2, c4; cutoff)
  add!(a2, a2, x)

  mul_tu!(a3, b4, c3; cutoff)
  mul_tu!(a4, b4, c4; cutoff)
#  @assert A == B*C
end

function mul_ut!(A::MatElem, B::MatElem, C::MatElem; cutoff::Int = cutoff)
  #A = BC for square matrices A, B, C, C is upper triangular
  @assert is_upper_triangular(C)
  n = nrows(A)
  n <= cutoff && return mul!(A, B, C)

  @assert nrows(A) == n
  m = ncols(A)
  k = ncols(B)
  @assert nrows(C) == k
  
  #A = [a1 a2; a3 a4]
  #B = [b1 b2; b3 b4]
  #C = [c1 c2;  0 c4]
  
  n2 = div(n, 2)
  m2 = div(k, 2)
  k2 = div(k, 2)

  a1 = view(A, 1:n2, 1:m2)
  a2 = view(A, 1:n2, 1+m2:m)
  a3 = view(A, 1+n2:n, 1:m2)
  a4 = view(A, 1+n2:n, 1+m2:m)

  b1 = view(B, 1:n2, 1:k2)
  b2 = view(B, 1:n2, 1+k2:k)
  b3 = view(B, 1+n2:n, 1:k2)
  b4 = view(B, 1+n2:n, 1+k2:k)

  c1 = view(C, 1:k2, 1:m2)
  c2 = view(C, 1:k2, 1+m2:m)
  c4 = view(C, 1+k2:k, 1+m2:m)

  mul_ut!(a1, b1, c1; cutoff)

  mul!(a2, b1, c2; cutoff)
  x = similar(a2)
  mul_ut!(x, b2, c4; cutoff)
  add!(a2, a2, x)

  mul_ut!(a3, b3, c1; cutoff)

  mul!(a4, b3, c2; cutoff)
  x = similar(a4)
  mul_ut!(x, b4, c4; cutoff)
  add!(a4, a4, x)
#  @assert A == B*C
end

end # module
