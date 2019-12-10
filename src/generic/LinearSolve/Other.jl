
################################################################################
#
#  Inconsistent system error.
#
################################################################################

struct InconsistentLinearSystemError <: Exception
    A::MatElem
    x::MatElem
    b::MatElem
end

function Base.showerror(io::IO, e::InconsistentLinearSystemError)
    A = e.A
    x = e.x
    b = e.b
    residual = A*x - b
    @error "Solve instance `Ax = b` is inconsistent. With: " A x b residual
end

################################################################################
#
#  Checks.
#
################################################################################

function check_solve_instance_is_well_defined(A::MatElem{T}, b::MatElem{T}) where T
    base_ring(A) != base_ring(b) && error("Base rings don't match in solve_lu")
    nrows(A) != nrows(b) && throw(DimensionMismatch("nrows(A) != nrows(b)"))

    if isempty(A) && !isempty(b)
        throw(DimensionMismatch("Solve Ax=b instance with A=$A empty and b=$b non-empty."))
    end
    return true
end

@doc Markdown.doc"""
    check_system_is_consistent(A,x,b,rk)

Checks if `A[rk+1:n, :]*x == b`, where `n` is the number of rows of `A`.
An error is raised otherwise.
"""
function check_system_is_consistent(A, x, b, rk = 0::Int)

    n = nrows(A)
    m = ncols(A)
    nvecs = max(ncols(x), ncols(b))

    # Containers
    t = base_ring(b)()
    sum = base_ring(b)()

    # Check all the dot products.
    for k = 1:nvecs
        for i = rk+1:n
            sum = zero!(sum)
            for j=1:m 
                sum = addmul_delayed_reduction!(sum, A[i, j], x[j, k], t)
            end
            sum = reduce!(sum)

            if sum != b[i,k]
                throw(InconsistentLinearSystemError(A[rk+1:n, :], x, b[rk+1:n, :]))
            end
        end
    end
    return true
end

# Solve the empty instance, assuming it is well-defined.
function _solve_empty(A,b)
    return similar(b, ncols(A), ncols(b))
end

################################################################################
#
#  Solve scaled.
#
################################################################################

# Need to think of a better name than `solve_ff`.

@doc Markdown.doc"""
    solve_scaled_ff(M::MatElem{T}, b::MatElem{T}; den = Val(true)) where {T <: RingElement}

Alias for solve_scaled_fflu.
"""
function solve_scaled_ff(M::MatElem{T}, b::MatElem{T}; den=Val(true)) where {T <: RingElement}
    return solve_fflu(M, b, den=den)
end


"""
    solve_ff(M::MatElem{T}, b::MatElem{T}) where {T <: FieldElement}

Gives a solution to the linear equation `A*x = b` using `solve_fflu`.
"""
function solve_ff(M::MatElem{T}, b::MatElem{T}) where {T <: FieldElement}    
    x, d = solve_fflu(M, b, den=Val(true))
    for i in 1:nrows(x)
        for j in 1:ncols(x)
            x[i, j] = divexact(x[i, j], d)
        end
    end
    return x
end

function solve_with_det(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
    return solve_fflu(M,b; den=Val(:det))
end

function solve_with_det(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: PolyElem}
   x, d = solve_interpolation(M, b)
   return x, d
end

################################################################################
#
#  Solve with interpolation
#
################################################################################

function solve_interpolation(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: PolyElem}
   m = nrows(M)
   h = ncols(b)
   if m == 0
      return b, base_ring(M)()
   end
   R = base_ring(M)
   maxlen = 0
   for i = 1:m
      for j = 1:m
         maxlen = max(maxlen, length(M[i, j]))
      end
   end
   maxlen == 0 && error("Matrix is singular in solve")
   maxlenb = 0
   for i = 1:m
      for j = 1:h
         maxlenb = max(maxlenb, length(b[i, j]))
      end
   end
   # bound from xd = (M*)b where d is the det
   bound = (maxlen - 1)*(m - 1) + max(maxlenb, maxlen)
   tmat = matrix(base_ring(R), 0, 0, elem_type(base_ring(R))[])
   V = Array{typeof(tmat)}(undef, bound)
   d = Array{elem_type(base_ring(R))}(undef, bound)
   y = Array{elem_type(base_ring(R))}(undef, bound)
   bj = Array{elem_type(base_ring(R))}(undef, bound)
   X = similar(tmat, m, m)
   Y = similar(tmat, m, h)
   x = similar(b)
   b2 = div(bound, 2)
   pt1 = base_ring(R)(1 - b2)
   l = 1
   i = 1
   pt = 1
   while l <= bound
      y[l] = base_ring(R)(pt - b2)
      (y[l] == pt1 && pt != 1) && error("Not enough interpolation points in ring")
      for j = 1:m
         for k = 1:m
            X[j, k] = evaluate(M[j, k], y[l])
         end
         for k = 1:h
            Y[j, k] = evaluate(b[j, k], y[l])
         end
      end
      try
         V[l], d[l] = solve_with_det(X, Y)
         l = l + 1
      catch e
          #TODO: One could make an optimization here. If you've found a specialization
          # that is inconsistent, this corresponds to a pole of the solution. This is quite,
          # useful information to know about the solution, and can reduce the size of
          # future solve instances.
         if e isa DivideError
             @warn ("Division error in solve_fflu. This occurs because the diagonal entries "*
                    "of the FFLU form are not sorted in any particular order. Please fix.")
         elseif !(e isa InconsistentLinearSystemError)
            rethrow(e)
         end
         i = i + 1
      end

      # We tested bound evaluation points and either an impossible
      # inverse was encountered, or the matrix was singular for all
      # the values.

      if i > bound && l == 1
         error("Impossible inverse or too many singular matrices in solve_interpolation")
      end

      pt = pt + 1
   end
   for k = 1:h
      for i = 1:m
         for j = 1:bound
            bj[j] = V[j][i, k]
         end
         x[i, k] = interpolate(R, y, bj)
      end
   end
   return x, interpolate(R, y, d)
end

################################################################################
#
#  Solve rational.
#
################################################################################

@doc Markdown.doc"""
    solve_rational(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where T <: RingElement
> Given a non-singular $n\times n$ matrix over a ring and an $n\times m$
> matrix over the same ring, return a tuple $x, d$ consisting of an
> $n\times m$ matrix $x$ and a denominator $d$ such that $Ax = d*b$. The
> denominator will be the determinant of $A$ up to sign. If $A$ is singular an
> exception is raised.
"""
function solve_rational(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where T <: RingElement
   return solve_scaled_ff(M, b, den=Val(true))
end


function solve_rational(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: PolyElem}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != ncols(M) && error("Non-square matrix in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   try
      return solve_interpolation(M, b)
   catch e
      if !isa(e, ErrorException)
         rethrow(e)
      end
      return solve_scaled_ff(M, b)
   end
end


################################################################################
#
#  The solve UI
#
################################################################################

@doc Markdown.doc"""
    solve(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where T
> Given a non-singular $n\times n$ matrix over a field and an $n\times m$
> matrix over the same field, return $x$ an
> $n\times m$ matrix $x$ such that $Ax = b$.
> If $A$ is singular an exception is raised.
"""
function solve(A::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T};
               kernel = Val(false),
               nullspace = Val(false),
               method = nothing,
               scaled = Val(false),
               side = :right
               ) where {T}

    kwds = Dict(:kernel => kernel,
                :nullspace => nullspace,
                :method => method,
                :scaled => scaled,
                :side => side
                )
    
    # Set up solve instance based on parameters, then call the right method, if possible.
    #
    # Ultimately solve is just a UI function designed to parse the solve instance and dispatch
    # to the method that will actually do the work.
    #
    # In the event `solve` cannot determine a method to solve `A*x == b`, it will give a reason
    # and a recommendation/override instructions.

    TRUE_LIST = [Val(true), Val(:true), true]
    FALSE_LIST = [Val(false), Val(:false), false]

    ####
    # Apply the user suggested method.
    if isa(method, Function)
        error("Generic keyword selection not implemented for calling user method.")
        return method(A,b, kwds... )
        
    elseif !isa(method, Nothing)
        typ = typeof(method)
        error("Provided method is not of type Function. typeof(method) = $typ")
    end

    ####
    # Check the side parameter for correctness.
    if side == Val(:right)
        @warn "`side` parameter is a symbol."
        side = :right
    elseif side == Val(:left)
        @warn "`side` parameter is a symbol."
        side = :left
    end

    ####
    # If unspecified, try to choose the right method.
    if T <: FieldElem
        if kernel in TRUE_LIST || nullspace in TRUE_LIST
            kernel = Val(true)
        end
        return solve_lu(A,b, side=side, kernel=kernel)
        
    elseif hasmethod(hnf_with_transform, Tuple{typeof(A)})
        if kernel in TRUE_LIST || nullspace in TRUE_LIST
            kernel = Val(true)
        else
            kernel = Val(false)
        end
        # If you want different behaviour for nullspace rather than kernel, change it
        # here:
        return solve_hnf(A,b, side=side, kernel=kernel)

    elseif scaled in TRUE_LIST
        if kernel in TRUE_LIST
            msg = ("Kernel not supported if no Hermite Normal Form method exists. "*
                   "Either use `nullspace=Val(true)` to get an integral basis for the "*
                   "kernel in the fraction field, or implement `hnf_with_transform`.")
            throw(DomainError(kernel, msg))
        elseif nullspace in TRUE_LIST
            nullspace = Val(true)
        else
            nullspace = Val(false)
        end
        
        return solve_scaled_ff(A,b, nullspace=nullspace)
    end

    # Don't know what to do, throw nice error message.
    return _solve_generic(A,b; kwds...)
end

function _solve_generic(A::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}; kwds...) where T
    display("Unfortunately, solve cannot solve your linear system `A*x = b` with:")
    @info "" typeof(A) typeof(b)
    display("This is because none of the following conditions held: ")
    @info "" (typeof(A) <: Field) hasmethod(hnf_with_transform, Tuple{typeof(A)}) kwds[:scaled]==Val(true)
    display("If it is sufficent to solve `A*x = d*b` for some non-zero `d`, set `scaled=Val(true)`.")
    display("If you control the type `T` of the base ring elements, and have a method in mind, you can write your own solve method and dispatch will take care of the rest.")
    display("finally, the list of available solve functions is: ")
    @info "" solve_lu solve_scaled_ff solve_rational
    display("Please see the documentation for more information.")
    error("No method found for solving `A*x = b` as given.")
end
