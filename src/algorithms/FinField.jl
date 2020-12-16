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

struct FinFieldIterator{F}
   basis::Vector{F}
   coeffs::Vector{Int64} # coefficients in the base field, assumed to be a prime field;
                         # even for characteristic > typemax(Int64), it's impossible to
                         # exhaust that many scalars with iteration

   function FinFieldIterator(f)
      deg = degree(f)
      BigInt(characteristic(f))^deg == order(f) ||
         throw(ArgumentError("iteration only supported for extension fields over a prime field"))
      basis = [one(f)]
      if deg > 1
         a = gen(f) # not defined if deg == 1
         for d = 2:deg
            push!(basis, basis[end]*a)
         end
      end
      new{elem_type(f)}(basis, Int64[])
   end
end

function Base.iterate(f::FinField, st=FinFieldIterator(f))
   basis, coeffs = st.basis, st.coeffs
   deg = degree(f)
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

function test_iterate(f::FinField)
   elts = collect(Iterators.take(f, 20))
   @test elts isa Vector{elem_type(f)}
   @test allunique(elts)
   @test length(elts) == min(order(f), 20)
   if order(f) < 100
      elts = collect(f)
      @test elts isa Vector{elem_type(f)}
      @test allunique(elts)
      @test length(elts) == order(f)
   end
   elts
end
