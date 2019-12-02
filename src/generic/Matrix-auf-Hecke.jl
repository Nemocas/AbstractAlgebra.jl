################################################################################
#
#  Matrix from Hecke
#
################################################################################
"""
This file contains functions transported from Hecke. It will act as a kind of
quarantine zone until the AbstractAlgebra package maintainers are happy with
the state of the code.
"""

################################################################################
#
#  Dense matrix types
#
################################################################################

dense_matrix_type(::Type{T}) where {T} = Generic.MatSpaceElem{T}

coefficient_type(::Type{Generic.Mat{T}}) where {T} = T

################################################################################
#
#  Matrix constructors
#
################################################################################

function zero_matrix(::Type{MatElem}, R::Ring, n::Int)
  return zero_matrix(R, n)
end

function zero_matrix(::Type{MatElem}, R::Ring, n::Int, m::Int)
  return zero_matrix(R, n, m)
end

function matrix(A::Array{T, 2}) where T <: RingElem
  r, c = size(A)
  (r < 0 || c < 0) && error("Array must be non-empty")
  m = matrix(parent(A[1, 1]), A)
  return m
end

function matrix(A::Array{T, 1}) where T <: RingElem
  return matrix(reshape(A,length(A),1))
end

function scalar_matrix(R::Ring, n::Int, a::RingElement)
  b = R(a)
  z = zero_matrix(R, n, n)
  for i in 1:n
    z[i, i] = b
  end
  return z
end

