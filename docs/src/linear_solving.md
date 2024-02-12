```@meta
CurrentModule = AbstractAlgebra.Solve
```

# [Linear solving](@id solving_chapter)

## Overview of the functionality

The module `AbstractAlgebra.Solve` provides the following four functions for solving linear systems:
* `solve`
* `can_solve`
* `can_solve_with_solution`
* `can_solve_with_solution_and_kernel`

All of these take the same set of arguments, namely:
* a matrix $A$ of type `MatElem`;
* a vector or matrix $B$ of type `Vector` or `MatElem`;
* a keyword argument `side` which can be either `:left` (default) or `:right`.

If `side` is `:left`, the system $xA = B$ is solved, otherwise the system $Ax = B$ is solved.

The functionality of the functions can be summarized as follows.
* `solve`: return a solution, if it exists, otherwise throw an error.
* `can_solve`: return `true`, if a solution exists, `false` otherwise.
* `can_solve_with_solution`: return `true` and a solution, if this exists, and `false` and an empty vector or matrix otherwise.
* `can_solve_with_solution_and_kernel`: like `can_solve_with_solution` and additionally return a matrix whose rows (respectively columns) give a basis of the kernel of $A$.

## Solving with several right hand sides

Systems $xA = b_1,\dots, xA = b_k$ with the same matrix $A$, but several right hand sides $b_i$ can be solved more efficiently, by first initializing a "context object" `C`.
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
kernel
```
