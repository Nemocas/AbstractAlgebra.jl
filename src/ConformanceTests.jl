module ConformanceTests

using ..AbstractAlgebra

# This file mostly contains function stubs.
# The actual implementation are in the folder `ext/TestExt/`.


# helper
function equality(a, b)
  if is_exact_type(typeof(a)) && is_exact_type(typeof(b))
     return a == b
  else
     return isapprox(a, b)
  end
end

function equality_up_to_units(a, b)
  iszero(a) && return iszero(b)
  iszero(b) && return iszero(a)
  return divides(a, b)[1] && divides(b, a)[1]
end

const default_adhoc_partner_rings = [
    AbstractAlgebra.Integers{BigInt}(),
    AbstractAlgebra.Integers{Int}(),
    AbstractAlgebra.Integers{UInt}(),
    AbstractAlgebra.Integers{UInt8}(),
  ]

adhoc_partner_rings(R::NCRing) = default_adhoc_partner_rings


#
# add methods for test_elem on ring elements here
#
function test_elem(R::AbstractAlgebra.Integers{T}) where {T <: Signed}
  n = T(2)^rand((1,1,1,2,3,10,31,32,33,63,64,65,100))
  return rand(R, -n:n)
end

function test_elem(R::AbstractAlgebra.Integers{T}) where {T <: Unsigned}
  n = T(2)^rand((1,1,1,2,3,10,31,32,33,63,64,65,100))
  return rand(R, 0:n)
end

function test_elem(R::AbstractAlgebra.Rationals)
  B = base_ring(R)
  n = test_elem(B)
  d = test_elem(B)
  return is_zero(d) ? R(n) : R(n, d)
end

function test_elem(R::AbstractAlgebra.FinField)
  return rand(R)
end

function test_elem(R::AbstractAlgebra.Floats{T}) where T
  return rand(T)*rand(-100:100)
end

function test_elem(Rx::AbstractAlgebra.PolyRing)
  R = base_ring(Rx)
  return Rx(elem_type(R)[test_elem(R) for i in 1:rand(0:6)])
end

function test_elem(Rx::AbstractAlgebra.MPolyRing)
  R = base_ring(Rx)
  num_gens = ngens(Rx)
  iszero(num_gens) && return Rx(test_elem(R))
  len_bound = 8
  exp_bound = rand(1:5)
  len = rand(0:len_bound)
  coeffs = [test_elem(R) for _ in 1:len]
  exps = [[rand(0:exp_bound) for _ in 1:num_gens] for _ in 1:len]
  return Rx(coeffs, exps)
end

function test_elem(S::Union{AbstractAlgebra.MatSpace,
                           AbstractAlgebra.MatRing})
  R = base_ring(S)
  return S(elem_type(R)[test_elem(R) for i in 1:nrows(S), j in 1:ncols(S)])
end

function test_elem(R::AbstractAlgebra.EuclideanRingResidueRing)
  return R(test_elem(base_ring(R)))
end

function test_elem(Rx::AbstractAlgebra.SeriesRing)
  R = base_ring(Rx)
  prec = rand(3:10)
  len = rand(0:prec-1)
  val = rand(0:prec-len)
  # FIXME: constructors don't seem to catch use of negative val
  @assert val >= 0
  A = elem_type(R)[test_elem(R) for i in 1:len]
  if len > 0 && is_zero(A[1])
    A[1] = one(R)
  end
  if elem_type(Rx) <: RelPowerSeriesRingElem
    @assert prec >= len + val
    return Rx(A, len, prec, val)
  else
    @assert prec >= len
    return Rx(A, len, prec)
  end
end

function test_elem(S::AbstractAlgebra.FreeAssociativeAlgebra)
  f = S()
  g = gens(S)
  R = base_ring(S)
  isempty(g) && return S(test_elem(R))
  len_bound = 8
  exp_bound = 6
  for i in 1:rand(0:len_bound)
     f += test_elem(R) * prod(rand(g) for _ in 1:rand(0:exp_bound); init = S(1))
  end
  return f
end


function test_iterate end

# Groups-conformance-tests.jl
function test_Group_interface end
function test_GroupElem_interface end

# Mutating-ops.jl
function test_mutating_op_like_zero end
function test_mutating_op_like_neg end
function test_mutating_op_like_add end
function test_mutating_op_like_addmul end

# Rings-conformance-tests.jl
function test_NCRing_interface end
function test_Ring_interface end
function test_Field_interface end
function test_EuclideanRing_interface end
function test_Poly_interface end
function test_MPoly_interface end
function test_MatSpace_interface end
function test_MatAlgebra_interface end
function test_Ring_interface_recursive end
function test_Field_interface_recursive end


end # module
