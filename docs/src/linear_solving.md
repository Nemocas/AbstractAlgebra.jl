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

The functionality of the functions can be summarized as follows.
* `solve`: return a solution, if it exists, otherwise throw an error.
* `can_solve`: return `true`, if a solution exists, `false` otherwise.
* `can_solve_with_solution`: return `true` and a solution, if this exists, and `false` and an empty vector or matrix otherwise.
* `can_solve_with_solution_and_kernel`: like `can_solve_with_solution` and additionally return a matrix whose columns (respectively rows) give a basis of the kernel of $A$.

## Detailed documentation

```@docs
solve
can_solve
can_solve_with_solution
can_solve_with_solution_and_kernel
```
