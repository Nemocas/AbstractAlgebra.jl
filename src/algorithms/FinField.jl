###############################################################################
#
#   FinField.jl : Generic algorithms for abstract finite fields
#
###############################################################################

###############################################################################
#
#   Iteration
#
###############################################################################

function _absolute_basis(f::FinField)
   basis = [one(f)]
   deg = degree(f)
   if deg > 1
      a = gen(f) # not defined if deg == 1
      for d = 2:deg
         push!(basis, basis[end]*a)
      end
   end
   return basis
end

struct FinFieldIterator{F}
   basis::Vector{F}
   coeffs::Vector{Int64} # coefficients in the base field, assumed to be a prime field;
                         # even for characteristic > typemax(Int64), it's impossible to
                         # exhaust that many scalars with iteration

   function FinFieldIterator(f)
      basis = _absolute_basis(f)
      deg = length(basis)
      BigInt(characteristic(f))^deg == order(f) ||
         throw(ArgumentError("iteration only supported for extension fields over a prime field"))
      new{elem_type(f)}(basis, Int64[])
   end
end

function Base.iterate(f::FinField, st=FinFieldIterator(f))
   basis, coeffs = st.basis, st.coeffs
   deg = length(basis)
   char = characteristic(f)
   elt = zero(f)
   if isempty(coeffs) # "flag" indicating that this is the first iteration
      resize!(coeffs, deg)
      fill!(coeffs, 0)
   else
      allzero = true
      for d=1:deg
         if allzero
            coeffs[d] += 1
            if coeffs[d] == char
               coeffs[d] = 0
            else
               allzero = false
            end
         end
         if !iszero(coeffs[d]) # add! only when necessary
            elt = add!(elt, elt, coeffs[d]*basis[d])
         end
      end
      allzero && return nothing
   end
   elt, st
end

Base.length(f::FinField) = BigInt(order(f))
Base.eltype(::Type{F}) where {F<:FinField} = elem_type(F)

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.test_elem(R::FinField)
  return rand(R)
end
