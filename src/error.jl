###############################################################################
#
#   Error objects
#
###############################################################################

struct NotImplementedError <: Exception
  head::Symbol
  args::Tuple
end

function NotImplementedError(h::Symbol, args...)
  return NotImplementedError(h, args)
end

function Base.showerror(io::IO, e::NotImplementedError)
  print(io, "function ")
  print(io, e.head)
  println(io, " is not implemented for argument" * "s"^(length(e.args) != 1))
  for a in e.args
    print(io, typeof(a))
    print(io, ": ")
    println(io, a)
  end
end


struct NotInvertibleError{T, S} <: Exception
  data::T   # element that could not be inverted
  mod::S    # ring or modulus with respect to which it could not be inverted
end

function NotInvertibleError(x::T, y::S) where {T <: NCRingElement, S}
  return NotInvertibleError{T, S}(x, y)
end

function NotInvertibleError(x::NCRingElement)
  return NotInvertibleError(x, parent(x))
end

function Base.showerror(io::IO, e::NotInvertibleError{T, S}) where {T, S <: Ring}
  print(io, e.data)
  print(io, " is not invertible in ")
  print(io, e.mod)
end

function Base.showerror(io::IO, e::NotInvertibleError)
  print(io, e.data)
  print(io, " is not invertible modulo ")
  print(io, e.mod)
end


struct ErrorConstrDimMismatch <: Exception
  expect_r::Int
  expect_c::Int
  get_r::Int
  get_c::Int
  get_l::Int

  function ErrorConstrDimMismatch(er::Int, ec::Int, gr::Int, gc::Int)
    e = new(er, ec, gr, gc, -1)
    return e
  end

  function ErrorConstrDimMismatch(er::Int, ec::Int, gl::Int)
    e = new(er, ec, -1, -1, gl)
    return e
  end

  function ErrorConstrDimMismatch(er::Int, ec::Int, a::Matrix{T}) where {T}
    gr, gc = size(a)
    return ErrorConstrDimMismatch(er, ec, gr, gc)
  end

  function ErrorConstrDimMismatch(er::Int, ec::Int, a::Vector{T}) where {T}
    gl = length(a)
    return ErrorConstrDimMismatch(er, ec, gl)
  end
end

function Base.showerror(io::IO, e::ErrorConstrDimMismatch)
  if e.get_l == -1
    print(io, "Expected dimension $(e.expect_r) x $(e.expect_c), ")
    print(io, "got $(e.get_r) x $(e.get_c)")
  else
    print(io, "Expected an array of length $(e.expect_r * e.expect_c), ")
    print(io, "got $(e.get_l)")
  end
end

struct InfiniteDimensionError <: Exception
  msg::String
  check_available::Bool

  @doc """
      InfiniteDimensionError(msg; check_available = false)

  A computation cannot be carried out because of an infinite-dimensional vector
  space. `msg` is a descriptive error message.

  Set `check_available` to `true` to add the line "You may check finiteness with
  `is_finite_dimensional_vector_space`" to the error message.

  # Examples

  ```jldoctest
  julia> throw(InfiniteDimensionError(check_available=true))
  ERROR: Infinite-dimensional vector space
  You may check finiteness with `is_finite_dimensional_vector_space`
  [...]

  julia> throw(InfiniteDimensionError("A custom message"))
  ERROR: A custom message
  [...]
  ```
  """
  function InfiniteDimensionError(msg::String = "Infinite-dimensional vector space";
                                  check_available::Bool = false)
    return new(msg, check_available)
  end
end

function Base.showerror(io::IO, e::InfiniteDimensionError)
  println(io, e.msg)
  if e.check_available
    println(io, "You may check finiteness with `is_finite_dimensional_vector_space`")
  end
end
