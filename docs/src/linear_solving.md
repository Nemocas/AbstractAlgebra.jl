```@meta
CurrentModule = AbstractAlgebra.Solve
```

# Linear solving

!!! note

    This functionality is experimental and subject to change.

## Overview of the functionality

The module `AbstractAlgebra.Solve` provides the following four functions for solving linear systems:
* `solve`
* `can_solve`
* `can_solve_with_solution`
* `can_solve_with_solution_and_kernel`

All of these take the same set of arguments, namely:
* a matrix $A$ of type `MatElem{T}`;
* a vector or matrix $B$ of type `Vector{T}` or `MatElem{T}`;
* a keyword argument `side` which can be either `:right` (default) or `:left`.

If `side` is `:right`, the system $Ax = B$ is solved, otherwise the system $xA = B$ is solved.
For matrices defined over a field, the functions internally rely on `rref`.
If the matrices are defined over a ring, the function `hnf_with_transform` is required internally.

The functionality of the functions can be summarized as follows.
* `solve`: return a solution, if it exists, otherwise throw an error.
* `can_solve`: return `true`, if a solution exists, `false` otherwise.
* `can_solve_with_solution`: return `true` and a solution, if this exists, and `false` and an empty vector or matrix otherwise.
* `can_solve_with_solution_and_kernel`: like `can_solve_with_solution` and additionally return a matrix whose columns (respectively rows) give a basis of the kernel of $A$.

## Solving with several right hand sides

Systems $Ax = b_1,\dots, Ax = b_k$ with the same matrix $A$, but several right hand sides $b_i$ can be solved more efficiently, by first initializing a "context object" `C` of type `SolveCtx`.
```@docs
solve_init
```
Now the functions `solve`, `can_solve`, etc. can be used with `C` in place of $A$.
This way the time-consuming part of the solving (i.e. computing a reduced form of $A$) is only done once and the result cached in `C` to be reused.

## Detailed documentation

```@docs
solve
can_solve
can_solve_with_solution
can_solve_with_solution_and_kernel
```
