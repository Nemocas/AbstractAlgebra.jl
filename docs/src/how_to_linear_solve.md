```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Solving Linear Systems](@id solving_linear_systems)


## Introduction

This page explains how to solve systems of linear equations. It introduces the
basic solving routines, the conventions used for left and right systems, and
related functionality such as kernel computations.


## Solving a left linear system

To solve a system of linear equations, use the function
[`solve`](@ref solve(::Union{MatElem{T}, Solve.SolveCtx{T}}, ::Union{Vector{T}, MatElem{T}}) where T).

!!! note
    By default, `solve` computes a solution `x` of a left system $xA = b$.
    Right linear systems are discussed in the next section.

```jldoctest
julia> A = matrix(QQ, [2 0 0; 0 3 0; 0 0 0])
[2//1   0//1   0//1]
[0//1   3//1   0//1]
[0//1   0//1   0//1]

julia> b = matrix(QQ, 1, 3, [1, 3//2, 0])
[1//1   3//2   0//1]

julia> x = solve(A, b)
[1//2   1//2   0//1]

julia> x * A == b
true
```

If the system is solvable, `solve` returns a particular solution. Otherwise,
it throws an error.

Since the matrix `A` in the above example does not have full rank, the system
may have infinitely many solutions. The full set of solutions is obtained by
adding arbitrary elements of the *left* kernel of `A` to the particular
solution returned by `solve`.

```jldoctest
julia> A = matrix(QQ, [2 0 0; 0 3 0; 0 0 0]);

julia> rank(A)
2

julia> kernel(A)
[0//1   0//1   1//1]
```

See the documentation of [`kernel`](@ref) for more information on computing kernels.


## Left and right linear systems

The function `solve` supports both left systems ($xA=b$) and right systems ($Ax=b$).
The keyword argument `side` determines which convention is used.

* `side = :left` (the default) solves the system `xA = b`.
* `side = :right` solves the system `Ax = b`.

The following example solves a right linear system:

```jldoctest
julia> A = matrix(QQ, [2 0 0; 0 3 0; 0 0 0])
[2//1   0//1   0//1]
[0//1   3//1   0//1]
[0//1   0//1   0//1]

julia> b = matrix(QQ, 3, 1, [1, 3//2, 0])
[1//1]
[3//2]
[0//1]

julia> x = solve(A, b, side = :right)
[1//2]
[1//2]
[0//1]

julia> A * x == b
true

julia> kernel(A, side = :right)
[0//1]
[0//1]
[1//1]
```


## Testing whether a system is solvable

The function [`solve`](@ref solve(::Union{MatElem{T}, Solve.SolveCtx{T}}, ::Union{Vector{T}, MatElem{T}}) where T)
throws an error if the system has no solution. If it is not known in advance
whether a solution exists, one of the following functions may be more
appropriate:

| Function | Description |
| :--- | :--- |
| [`can_solve`](@ref can_solve(::Union{MatElem{T}, Solve.SolveCtx{T}}, ::Union{Vector{T}, MatElem{T}}) where T) | Return whether the system is solvable. |
| [`can_solve_with_solution`](@ref can_solve_with_solution(::Union{MatElem{T}, Solve.SolveCtx{T}}, ::Union{Vector{T}, MatElem{T}}) where T) | Return whether the system is solvable and, if so, a particular solution. |
| [`can_solve_with_solution_and_kernel`](@ref can_solve_with_solution_and_kernel(::Union{MatElem{T}, Solve.SolveCtx{T}}, ::Union{Vector{T}, MatElem{T}}) where T) | Return whether the system is solvable and, if so, a particular solution together with a basis of the corresponding kernel. |

For example:

```jldoctest
julia> A = matrix(QQ, [2 0 0; 0 3 0; 0 0 0]);

julia> b = matrix(QQ, 1, 3, [1, 3//2, 1]);

julia> can_solve(A, b)
false

julia> b2 = matrix(QQ, 1, 3, [1, 3//2, 0]);

julia> can_solve_with_solution(A, b2)
(true, [1//2 1//2 0])

julia> can_solve_with_solution_and_kernel(A, b2)
(true, [1//2 1//2 0], [0 0 1])
```


## Solving multiple linear systems simultaneously

The functions discussed above also accept a matrix as the right-hand side. In
this case, each row (respectively column) represents a separate right-hand
side, and all systems are solved simultaneously.

If several systems with the same coefficient matrix but different right-hand
sides need to be solved one after another, it is more efficient to construct
a solving context first. To construct such a solving context, use the function
[`solve_init`](@ref solve_init(::MatElem)). The resulting context object can
be used in place of the coefficient matrix in calls to `solve`, `can_solve`,
and the related functions. This avoids recomputing the reduced form of the
coefficient matrix each time a new system is solved.


## Reference documentation

For a complete description of the available solving functions and solving
contexts, see the [Linear solving interface](@ref solving_chapter).
