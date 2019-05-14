```@meta
CurrentModule = AbstractAlgebra
```

# Generic matrices

AbstractAlgebra.jl allows the creation of dense matrices over any computable commutative
ring $R$. Generic matrices over a commutative ring are implemented in
`src/generic/Matrix.jl`. Much of the functionality there covers both matrix spaces and
matrix algebras.

Functions specific to generic matrix algebras of $m\times m$ matrices are implemented in 
`src/generic/MatrixAlgebra.jl`.

As well as implementing the entire Matrix interface, including the optional
functionality, there are many additional generic algorithms implemented for matrix
spaces. We describe this functionality below.

All of this generic functionality is part of the Generic submodule of
AbstractAlgebra.jl. This is exported by default, so it is not necessary to qualify names
of functions.

## Types and parent objects

Generic matrices in AbstractAlgebra.jl have type `Generic.MatSpaceElem{T}` for matrices
in a matrix space, or `Generic.MatAlgElem{T}` for matrices in a matrix algebra, where
`T` is the type of elements of the matrix. Internally, generic matrices are implemented
using an object wrapping a Julia two dimensional array, though they are not themselves
Julia arrays. See the file `src/generic/GenericTypes.jl` for details.

For the most part, one doesn't want to work directly with the `MatSpaceElem` type though,
but with an abstract type called `Generic.Mat` which includes `MatSpaceElem` and views
thereof.

Parents of generic matrices (matrix spaces) have type `Generic.MatSpace{T}`. Parents of
matrices in a matrix algebra have type `Generic.MatAlgebra{T}`.

The generic matrix types (matrix spaces) belong to the abstract type
`AbstractAlgebra.MatElem{T}` and the matrix space parent types belong to
`AbstractAlgebra.MatSpace{T}`. Similarly the generic matrix algebra matrix types belong
to the abstract type `AbstractAlgebra.MatAlgElem{T}` and the parent types belong to
 `AbstractAlgebra.MatAlgebra{T}` Note that both
the concrete type of a matrix space parent object and the abstract class it belongs to
have the name `MatElem`, therefore disambiguation is required to specify which is
intended. The same is true for the abstract types for matrix spaces and their elements.

The dimensions and base ring $R$ of a generic matrix are stored in its parent object,
however to allow creation of matrices without first creating the matrix space parent,
generic matrices in Julia do not contain a reference to their parent. They contain the
row and column numbers (or degree, in the case of matrix algebras) and the base ring
on a per matrix basis. The parent object can then be reconstructed from this data on
demand.

## Matrix space constructors

A matrix space in AbstractAlgebra.jl represents a collection of all matrices with
given dimensions and base ring.

In order to construct matrices in AbstractAlgebra.jl, one can first construct the
matrix space itself. This is accomplished with the following constructor. We discuss
creation of matrix algebras separately in a dedicated section elsewhere in the
documentation.

```julia
MatrixSpace(R::Ring, rows::Int, cols::Int; cache::Bool=true)
```

Construct the space of matrices with the given number of rows and columns over the
given base ring. By default such matrix spaces are cached based on the base ring and
numbers of rows and columns. If the optional named parameter `cached` is set to false,
no caching occurs.

Here are some examples of creating matrix spaces and making use of the
resulting parent objects to coerce various elements into the matrix space.

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S = MatrixSpace(R, 3, 3)

A = S()
B = S(12)
C = S(R(11))
```

We also allow matrices over a given base ring to be constructed directly (see the
Matrix interface).

## Matrix element constructors

In addition to coercing elements into a matrix space as above, we provide the
following functions for constructing explicit matrices.

Also see the Matrix interface for a list of other ways to create matrices.

```julia
R[a b c...;...]
```

Create the matrix over the base ring $R$ consisting of the given rows (separated by
semicolons). Each entry is coerced into $R$  automatically. Note that parentheses may
be placed around individual entries if the lists would otherwise be ambiguous, e.g. 
`R[1 2; 2 (-3)]`.

Beware that this syntax does not support the creation of column vectors. See
the notation below for creating those.

```julia
R[a b c...]
```

Create the row vector with entries in $R$ consisting of the given entries
(separated by spaces). Each entry is coerced into $R$ automatically. Note that
parentheses may be placed around individual entries if the list would otherwise
be ambiguous, e.g. `R[1 2 (-3)]`.

```julia
R[a b c...]'
```

Create the column vector with entries in $R$ consisting of the given entries
(separated by spaces). Each entry is coerced into $R$ automatically. Observe
the dash that is used to transpose the row vector notation (for free) to turn
it into a column vector. Note that parentheses may be placed around individual
entries if the list would otherwise be ambiguous, e.g. `R[1 2 (-3)]'`.

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S = MatrixSpace(R, 3, 3)

M = R[t + 1 1; t^2 0]
N = R[t + 1 2 t]
P = R[1 2 t]'
```

## Submatrices

In addition to the functionality described in the Matrix interface for taking
submatrices of a matrix, the following function variant is also available.

```@docs
sub(::MatElem, ::Int, ::Int, ::Int, ::Int)
```

**Examples**

```julia
M = ZZ[1 2 3; 2 3 4]

N = sub(M, 1, 1, 2, 2)
```

## Matrix functionality provided by AbstractAlgebra.jl

Most of the following generic functionality is available for both matrix spaces and
matrix algebras. Exceptions include functions that do not return or accept square
matrices or which cannot specify a parent. Such functions include `solve`, `kernel`,
and `nullspace` which can't be provided for matrix algebras.

For details on functionality that is provided for matrix algebras only, see the dedicated
section of the documentation.

### Basic matrix functionality

As well as the Ring and Matrix interfaces, the following functions are provided to
manipulate matrices and to set and retrieve entries and other basic data associated
with the matrices.

```@docs
dense_matrix_type(::Ring)
```

```@docs
nrows(::MatElem)
```

```@docs
ncols(::MatElem)
```

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S = MatrixSpace(R, 3, 3)

A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])

T = dense_matrix_type(R)
r = nrows(B)
c = ncols(B)
M = A + B
N = 2 + A
M1 = deepcopy(A)
A != B
isone(one(S)) == true
V = A[1:2, :]
W = A^3
Z = divexact(2*A, 2)
```

### Powering

```@docs
powers(::MatElem, ::Int)
```

**Examples**

```julia
M = ZZ[1 2 3; 2 3 4; 4 5 5]

A = powers(M, 4)
```

### Gram matrix

```@docs
gram(::MatElem)
```

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S = MatrixSpace(R, 3, 3)

A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

B = gram(A)
```

### Trace

```@docs
tr(::MatElem)
```

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S = MatrixSpace(R, 3, 3)

A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

b = tr(A)
```

### Content

```@docs
content(::MatElem)
```

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S = MatrixSpace(R, 3, 3)

A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

b = content(A)
```

### Permutation

```@docs
*(::perm, ::MatElem)
```

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S = MatrixSpace(R, 3, 3)
G = PermGroup(3)

A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
P = G([1, 3, 2])

B = P*A
```

### LU factorisation

```@docs
lu{T <: FieldElem}(::MatElem{T}, ::PermGroup)
```

```@docs
fflu{T <: RingElem}(::MatElem{T}, ::PermGroup)
```

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")
K, a = NumberField(x^3 + 3x + 1, "a")
S = MatrixSpace(K, 3, 3)

A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 - 2 a - 1 2a])

r, P, L, U = lu(A)
r, d, P, L, U = fflu(A)
```

### Reduced row-echelon form

```@docs
rref{T <: RingElem}(::MatElem{T})
rref{T <: FieldElem}(::MatElem{T})
```

```@docs
isrref{T <: RingElem}(::MatElem{T})
isrref{T <: FieldElem}(::MatElem{T})
```

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")
K, a = NumberField(x^3 + 3x + 1, "a")
S = MatrixSpace(K, 3, 3)
   
M = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
   
r, A = rref(M)
isrref(A)

R, x = PolynomialRing(ZZ, "x")
S = MatrixSpace(R, 3, 3)

M = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])

r, d, A = rref(M)
isrref(A)
```

### Hermite normal form

```@docs
hnf{T <: RingElem}(::MatElem{T})
hnf_with_transform{T <: RingElem}(::MatElem{T})
```

### Determinant

```@docs
det{T <: RingElem}(::MatElem{T})
det{T <: FieldElem}(::MatElem{T})
```

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")
K, a = NumberField(x^3 + 3x + 1, "a")
S = MatrixSpace(K, 3, 3)
   
A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])

d = det(A)
```

### Rank

```@docs
rank{T <: RingElem}(::MatElem{T})
rank{T <: FieldElem}(::MatElem{T})
```

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")
K, a = NumberField(x^3 + 3x + 1, "a")
S = MatrixSpace(K, 3, 3)
   
A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])

d = rank(A)
```

### Linear solving

```@docs
solve{T <: FieldElem}(::MatElem{T}, ::MatElem{T})
```

```@docs
solve_rational{T <: RingElem}(::MatElem{T}, ::MatElem{T})
```

```@docs
solve_triu{T <: FieldElem}(::MatElem{T}, ::MatElem{T}, ::Bool)
```

```@docs
can_solve_left_reduced_triu{T <: RingElement}(::MatElem{T}, ::MatElem{T})
```

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")
K, a = NumberField(x^3 + 3x + 1, "a")
S = MatrixSpace(K, 3, 3)
U = MatrixSpace(K, 3, 1)

A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
b = U([2a a + 1 (-a - 1)]')

x = solve(A, b)

A = S([a + 1 2a + 3 a^2 + 1; K(0) a^2 - 1 2a; K(0) K(0) a])
b = U([2a a + 1 (-a - 1)]')

x = solve_triu(A, b, false)

R, x = PolynomialRing(ZZ, "x")
S = MatrixSpace(R, 3, 3)
U = MatrixSpace(R, 3, 2)

A = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])
b = U([2x x + 1 (-x - 1); x + 1 (-x) x^2]')

x, d = solve_rational(A, b)

S = MatrixSpace(ZZ, 3, 3)
T = MatrixSpace(ZZ, 3, 1)

A = S([BigInt(2) 3 5; 1 4 7; 9 2 2])   
B = T([BigInt(4), 5, 7])

X, d = solve_rational(A, B)
```

### Inverse

```@docs
inv{T <: RingElem}(::MatElem{T})
inv{T <: FieldElem}(::MatElem{T})
```

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")
K, a = NumberField(x^3 + 3x + 1, "a")
S = MatrixSpace(K, 3, 3)

A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])

X = inv(A)

R, x = PolynomialRing(ZZ, "x")
S = MatrixSpace(R, 3, 3)

A = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])
    
X, d = inv(A)
```

### Nullspace

```@docs
nullspace{T <: RingElem}(::MatElem{T})
nullspace{T <: FieldElem}(::MatElem{T})
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
S = MatrixSpace(R, 4, 4)
   
M = S([-6*x^2+6*x+12 -12*x^2-21*x-15 -15*x^2+21*x+33 -21*x^2-9*x-9;
       -8*x^2+8*x+16 -16*x^2+38*x-20 90*x^2-82*x-44 60*x^2+54*x-34;
       -4*x^2+4*x+8 -8*x^2+13*x-10 35*x^2-31*x-14 22*x^2+21*x-15;
       -10*x^2+10*x+20 -20*x^2+70*x-25 150*x^2-140*x-85 105*x^2+90*x-50])
   
n, N = nullspace(M)
```

### Kernel

```@docs
kernel{T <: RingElem}(::MatElem{T})
left_kernel{T <: RingElem}(::MatElem{T})
right_kernel{T <: RingElem}(::MatElem{T})
```

**Examples**

```julia
S = MatrixSpace(ZZ, 4, 4)

M = S([1 2 0 4;
       2 0 1 1;
       0 1 1 -1;
       2 -1 0 2])

nr, Nr = kernel(M)
nl, Nl = left_kernel(M)
```

### Hessenberg form

```@docs
hessenberg{T <: RingElem}(::MatElem{T})
```

```@docs
ishessenberg{T <: RingElem}(::MatElem{T})
```

**Examples**

```julia
R = ResidueRing(ZZ, 7)
S = MatrixSpace(R, 4, 4)
   
M = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);
       R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])
   
A = hessenberg(M)
ishessenberg(A) == true
```

### Characteristic polynomial

```@docs
charpoly{T <: RingElem}(::Ring, ::MatElem{T})
```

**Examples**

```julia
R = ResidueRing(ZZ, 7)
S = MatrixSpace(R, 4, 4)
T, x = PolynomialRing(R, "x")

M = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);
       R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])
   
A = charpoly(T, M)
```

### Minimal polynomial

```@docs
minpoly{T <: RingElem}(::Ring, ::MatElem{T}, ::Bool)
minpoly{T <: FieldElem}(::Ring, ::MatElem{T}, ::Bool)
```

**Examples**

```julia
R = GF(13)
T, y = PolynomialRing(R, "y")
   
M = R[7 6 1;
      7 7 5;
      8 12 5]

A = minpoly(T, M)
```

### Transforms

```@docs
similarity!{T <: RingElem}(::MatElem{T}, ::Int, ::T)
```

**Examples**

```julia
R = ResidueRing(ZZ, 7)
S = MatrixSpace(R, 4, 4)
   
M = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);
       R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])
   
similarity!(M, 1, R(3))
```

### (Weak) Popov form

AbstractAlgebra.jl provides algorithms for computing the (weak) Popov of a matrix with
entries in a univariate polynomial ring over a field.

```@docs
weak_popov{T <: PolyElem}(::Generic.Mat{T})
weak_popov_with_transform{T <: PolyElem}(::Generic.Mat{T})
popov{T <: PolyElem}(::Generic.Mat{T})
popov_with_transform{T <: PolyElem}(::Generic.Mat{T})
```

