module Solve

import AbstractAlgebra: FieldElement, MatElem

@doc raw"""
    solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement

Return a vector $x$ solving the linear system $Ax = b$, if `side == :right`
(default), or $xA = b$, if `side == :left`.

If no solution exists, an error is raised.

See also [`can_solve_with_solution`](@ref).
"""
function solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement
end

@doc raw"""
    solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement

Return a matrix $x$ solving the linear system $Ax = b$, if `side == :right`
(default), or $xA = b$, if `side == :left`.

If no solution exists, an error is raised.

See also [`can_solve_with_solution`](@ref).
"""
function solve(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement
end

@doc raw"""
    can_solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement
    can_solve(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement

Return `true` if the linear system $Ax = b$ or $xA = b$ with `side == :right`
(default) or `side == :left`, respectively, has a solution and `false` otherwise.

See also [`can_solve_with_solution`](@ref).
"""
can_solve

function can_solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement
end

function can_solve(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement
end

@doc raw"""
    can_solve_with_solution(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement

Return `true` and a vector $x$ solving the linear system $Ax = b$, if such a
solution exists. Return `false` and a "dummy vector", if the system has no solution.

If `side == :left`, the system $xA = b$ is solved.

See also [`solve`](@ref).
"""
function can_solve_with_solution(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement
end

@doc raw"""
    can_solve_with_solution(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement

Return `true` and a matrix $x$ solving the linear system $Ax = b$, if such a
solution exists. Return `false` and a "dummy matrix", if the system has no solution.

If `side == :left`, the system $xA = b$ is solved.

See also [`solve`](@ref).
"""
function can_solve_with_solution(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement
end

@doc raw"""
    can_solve_with_solution_and_kernel(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement

Return `true`, a vector $x$ solving the linear system $Ax = b$, together with a
matrix $K$ giving the kernel of $A$ (i.e. $AK = 0$), if such a solution exists.
Return `false`, a "dummy vector", and the kernel, if the system has no solution.

If `side == :left`, the system $xA = b$ is solved.

See also [`solve`](@ref) and [`kernel`](@ref).
"""
function can_solve_with_solution_and_kernel(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement
end

@doc raw"""
    can_solve_with_solution_and_kernel(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement

Return `true`, a matrix $x$ solving the linear system $Ax = b$, together with a
matrix $K$ giving the kernel of $A$ (i.e. $AK = 0$), if such a solution exists.
Return `false`, a "dummy matrix", and the kernel, if the system has no solution.

If `side == :left`, the system $xA = b$ is solved.

See also [`solve`](@ref) and [`kernel`](@ref).
"""
function can_solve_with_solution_and_kernel(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement
end

end
