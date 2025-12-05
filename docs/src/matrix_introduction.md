```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```
# Introduction

AbstractAlgebra provides matrix spaces ($m\times n$ matrices) and matrix algebras
($n\times n$ matrices) over a ring. Whilst both types of matrix provide
matrix multiplication for matrices whose dimensions are compatible for
multiplication, only the latter kind of matrices form rings in the system.

Matrix spaces provide a large number of linear algebra operations, including
linear solving, elementary row operations, various canonical forms. The
system also provides characteristic and minimal polynomial computations, LU
decomposition, determinant, matrix inverse, kernel computations.

There is also code for computation of the Hermite and Smith normal forms over
Euclidean domains and Popov form for matrices over polynomial rings over a
field.

Most of this generic functionality is provided for arbitrary matrix types
that satisfy the AbstractAlgebra matrix interface.
