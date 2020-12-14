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
   current::Vector{Int}

   function FinFieldIterator(f)
      deg = degree(f)
      basis = [one(f)]
      if deg > 1
         a = gen(f) # not defined if deg == 1
         for d = 2:deg
            push!(basis, basis[end]*a)
         end
      end
      new{elem_type(f)}(basis, Int[])
   end
end

function Base.iterate(f::FinField, st=FinFieldIterator(f))
   basis, current = st.basis, st.current
   deg = degree(f)
   char = Int(characteristic(f)) # we currently support only "small" characteristic
   elt = zero(f)
   if isempty(current) # "flag" indicating that this is the first iteration
      resize!(current, deg)
      fill!(current, 0)
   else
      allzero = true
      for d=1:deg
         if allzero
            current[d] += 1
            if current[d] == char
               current[d] = 0
            else
               allzero = false
            end
         end
         if !iszero(current[d]) # add! only when necessary
            elt = add!(elt, elt, current[d]*basis[d])
         end
      end
      allzero && return nothing
   end
   elt, st
end

Base.length(f::FinField) = Int(order(f))
Base.eltype(f::FinField) = elem_type(f)

function test_iterate(f::FinField)
   elts = collect(f)
   @test elts isa Vector{elem_type(f)}
   @test allunique(elts)
   elts
end
