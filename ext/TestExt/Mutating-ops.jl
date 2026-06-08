# The following functions should not expect that their input is a `NCRingElem` or similar.
# They should be usable in more general types, that don't even have a `parent/elem` correspondence

@nospecialize

function test_mutating_op_like_zero(f::Function, f!::Function, A)
  a = deepcopy(A)
  a = f!(a)
  @test equality(a, f(A))
end

function test_mutating_op_like_neg(f::Function, f!::Function, A)
  # initialize storage var with different values to check that its value is not used
  for z in [zero(A), deepcopy(A)]
     a = deepcopy(A)
     z = f!(z, a)
     @test equality(z, f(A))
     @test a == A
  end

  a = deepcopy(A)
  a = f!(a)
  @test equality(a, f(A))
end

function test_mutating_op_like_add(f::Function, f!::Function, A, B, T = Any; only3arg::Bool = false)
  # only3arg = don't test f!(a, b)
  @req A isa T || B isa T "Invalid argument types"

  # initialize storage var with different values to check that its value is not used
  storage_values = T[]
  if A isa T
     push!(storage_values, zero(A))
     push!(storage_values, deepcopy(A))
  end
  if B isa T
     push!(storage_values, zero(B))
     push!(storage_values, deepcopy(B))
  end
  for z in storage_values
     a = deepcopy(A)
     b = deepcopy(B)
     z = f!(z, a, b)
     @test equality(z, f(A, B))
     @test a == A
     @test b == B
  end

  if A isa T
     a = deepcopy(A)
     b = deepcopy(B)
     a = f!(a, a, b)
     @test equality(a, f(A, B))
     @test b == B

     if !only3arg
       a = deepcopy(A)
       b = deepcopy(B)
       a = f!(a, b)
       @test equality(a, f(A, B))
       @test b == B
     end
  end

  if B isa T
     a = deepcopy(A)
     b = deepcopy(B)
     b = f!(b, a, b)
     @test equality(b, f(A, B))
     @test a == A
  end

  if A isa T && B isa T
     # `f(B, B)` may fail if `!(A isa T)`, since we call it with different arguments than the intended `f(A, B)` (same for f!)
     a = deepcopy(A)
     b = deepcopy(B)
     a = f!(a, b, b)
     @test equality(a, f(B, B))
     @test b == B

     b = deepcopy(B)
     b = f!(b, b, b)
     @test equality(b, f(B, B))

     if !only3arg
       b = deepcopy(B)
       b = f!(b, b)
       @test equality(b, f(B, B))
     end
  end
end

function test_mutating_op_like_addmul(f::Function, f!_::Function, Z, A, B, T = Any)
  @req Z isa T "Invalid argument types"
  @req A isa T || B isa T "Invalid argument types"

  f!(z, a, b, ::Nothing) = f!_(z, a, b)
  f!(z, a, b, t) = f!_(z, a, b, t)

  # initialize storage var with different values to check that its value is not used
  # and `nothing` for the three-arg dispatch
  storage_values = Union{T,Nothing}[nothing]
  if A isa T
     push!(storage_values, zero(A))
     push!(storage_values, deepcopy(A))
  end
  if B isa T
     push!(storage_values, zero(B))
     push!(storage_values, deepcopy(B))
  end
  for t in storage_values
     z = deepcopy(Z)
     a = deepcopy(A)
     b = deepcopy(B)
     z = f!(z, a, b, t)
     @test equality(z, f(Z, A, B))
     @test a == A
     @test b == B

     if A isa T
        a = deepcopy(A)
        b = deepcopy(B)
        a = f!(a, a, b, t)
        @test equality(a, f(A, A, B))
        @test b == B
     end

     if B isa T
        a = deepcopy(A)
        b = deepcopy(B)
        b = f!(b, a, b, t)
        @test equality(b, f(B, A, B))
        @test a == A
     end

     if A isa T && B isa T
        # `f(B, B)` may fail if `!(A isa T)`, since we call it with different arguments than the intended `f(A, B)` (same for f!)
        a = deepcopy(A)
        b = deepcopy(B)
        a = f!(a, b, b, t)
        @test equality(a, f(A, B, B))
        @test b == B

        b = deepcopy(B)
        b = f!(b, b, b, t)
        @test equality(b, f(B, B, B))
     end
  end
end

@specialize
