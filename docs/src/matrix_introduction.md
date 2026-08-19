```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```
# Introduction

AbstractAlgebra provides exact linear algebra over a wide range of rings,
including the integers, finite fields, polynomial rings and residue rings.
Unlike Julia's standard linear algebra functionality, matrices are not
restricted to numerical coefficient domains.

Matrices in AbstractAlgebra belong to parent objects. Rectangular
$m \times n$ matrices are elements of matrix spaces, while square
$n \times n$ matrices may additionally be constructed as elements of matrix
algebras, where matrix multiplication gives them the structure of a ring.

The matrix functionality includes standard arithmetic, linear solving,
elementary row and column operations, determinants, inverses, kernels,
characteristic and minimal polynomials, and various decompositions and
canonical forms. Additional specialised algorithms are available over
suitable base rings, such as Hermite and Smith normal forms over Euclidean
domains and Popov forms over polynomial rings.

The following sections describe how to construct matrices, access and modify
their entries, compute with them, and use the parent objects associated to
matrices.
