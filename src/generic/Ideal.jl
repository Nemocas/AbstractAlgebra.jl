###############################################################################
#
#   Ideal.jl : Generic ideals over Euclidean domains
#
###############################################################################

###############################################################################
#
#   Type and parent functions
#
###############################################################################

base_ring(S::IdealSet) = S.base_ring

base_ring(I::Ideal) = I.base_ring

function parent(I::Ideal)
   R = base_ring(I)
   return IdealSet{elem_type(R)}(R)
end

elem_type(::Type{IdealSet{S}}) where S <: RingElement = Ideal{S}

parent_type(::Type{Ideal{S}}) where S <: RingElement = IdealSet{S}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

gens(I::Ideal) = I.gens

###############################################################################
#
#   Ideal reduction for mulivariates over Euclidean domain
#
###############################################################################

mutable struct lmnode{U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
   poly::U
   up::Union{lmnode{U}, Nothing}   # out (divisible leading monomials)
   next::Union{lmnode{U}, Nothing} # extend current node so out-degree can be > 1
   lm::NTuple{N, Int}   # leading monomial as exponent vector
   lcm::NTuple{N, Int}  # lcm of lm's in tree rooted here, as exponent vector
   in_heap::Bool # whether node is in the heap

   function lmnode{U, N}(p::U) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
      node = new{U, N}(p, nothing, nothing)
      node.lm = Tuple(exponent_vector(p, 1))
      node.lcm = node.lm
      return node::lmnode{U, N}
   end
end

# heap implementation for sorting polys by lm, head = smallest

# Defined in MPoly, repeated here for reference
# heapleft(i::Int) = 2i
# heapright(i::Int) = 2i + 1
# heapparent(i::Int) = div(i, 2)

function lm_precedes(node1::lmnode{U, N}, node2::lmnode{U, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
   for i = 1:N
      if node2.lm[i] > node1.lm[i]
         return false
      end
   end
   return true
end

function lm_lcm(node1::lmnode{U, N}, node2::lmnode{U, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
   return Tuple(max(node1.lcm[i], node2.lcm[i]) for i in 1:N)
end

function heapinsert!(heap::Vector{lmnode{U, N}}, node::lmnode{U, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
   i = n = length(heap) + 1
   if i != 1 && lm_precedes(heap[1], node)
      i = 1
   else
      @inbounds while (j = heapparent(i)) >= 1
         if lm_precedes(heap[j], node)
            i = j
         else
            break
         end
      end
   end
   push!(heap, node)
   while n > i
      heap[n] = heap[heapparent(n)]
      n >>= 1
   end
   heap[i] = node
   node.in_heap = true
   return Nothing
end

function heappop!(heap::Vector{lmnode{U, N}}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
   s = length(heap)
   x = heap[1]
   i = 1
   j = 2
   @inbounds while j < s
      if !lm_precedes(heap[j + 1], heap[j])
         j += 1
      end
      heap[i] = heap[j]
      i = j
      j <<= 1
   end
   n = heap[s]
   j = i >> 1
   @inbounds while i > 1 && lm_precedes(heap[j], n)
      heap[i] = heap[j]
      i = j
      j >>= 1
   end
   heap[i] = heap[s]
   pop!(heap)
   x.in_heap = false
   return x
end

function extract_gens(V::Vector{U}, node::T) where {N, U <: AbstractAlgebra.MPolyElem, T <: lmnode{U, N}}
   node.in_heap = true # mark node as visited
   # depth first
   if node.up != nothing
      V = extract_gens(V, node.up)
   end
   if node.next != nothing
      V = extract_gens(V, node.next)
   end
   push!(V, node.poly)
   return V
end

function extract_gens(B::Vector{T}) where {N, U <: AbstractAlgebra.MPolyElem, T <: lmnode{U, N}}
   V = Vector{U}()
   for node in B
      V = extract_gens(V, node)
   end
   return V
end

function reduce(I::Ideal{U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
   if hasmethod(gcdx, Tuple{T, T})
      V = gens(I)
      # Step 1: compute V a vector of polynomials giving the same basis as I
      #         but for which 1, 2, 3 above hold
      if length(V) > 1
         # make heap
         N = nvars(base_ring(I))
         heap = Vector{lmnode{U, N}}()
         for v in V
            heapinsert!(heap, lmnode{U, N}(v))
         end
         d = heappop!(heap)
         # take gcd of constant polys
         if isconstant(d.poly)
            S = parent(d.poly)
            d0 = constant_coefficient(d.poly)
            while !isempty(heap) && isconstant(heap[1].poly)
               di = constant_coefficient(pop!(heap).poly)
               d0 = gcd(d0, di)
            end
            d.poly = S(d0)
         end
         # canonicalise
         d.poly = divexact(d.poly, canonical_unit(d.poly))
         B = [d]
         # do reduction
         # extract polynomials from B
         V = extract_gens(B)
      end
      return Ideal{U}(base_ring(I), V)
   else
      error("Not implemented")
   end
end

###############################################################################
#
#   Ideal reduction for univariates over Euclidean domain
#
###############################################################################

# Extend the basis V of polynomials satisfying 1-6 below by the
# polynomials in D, all of which have degree at least that of those in
# V and such that the degrees of the polynomials in D is *decreasing*.
function extend_ideal_basis(D::Vector{T}, V::Vector{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   while !isempty(D)
      d = pop!(D)
      V = extend_ideal_basis(d, V)
   end
   return V
end

function reduce_tail(p::T, V::Vector{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   p = divexact(p, canonical_unit(p))
   i = length(V)
   for n = length(p) - 1:-1:1
      while i > 0 && length(V[i]) > n
         i -= 1
      end
      if i != 0
         q = AbstractAlgebra.div(coeff(p, n - 1), leading_coefficient(V[i]))
         p -= q*shift_left(V[i], n - length(V[i]))
      end
   end
   return p
end

# Given a nonempty vector V of polynomials of satisfying 1-6 below and a
# polynomial p whose degree is at least that of all the polynomials in V, add
# p to V and perform reduction steps so that 1-6 still hold.
function extend_ideal_basis(p::T, V::Vector{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   n = length(V)
   lc = leading_coefficient(V[n])
   # check if p can be added without any reduction
   if length(p) > length(V[n]) && !isunit(lc) && ((_, q) = divides(lc, leading_coefficient(p)))[1]
      return vcat(V, [reduce_tail(p, V)])
   end
   # check if p and V[n] are constant
   if isconstant(V[n]) && isconstant(p)
      return [parent(p)(gcd(constant_coefficient(p), constant_coefficient(V[n])))]
   end
   # check if p can replace V[n]
   swap = false
   if length(p) == length(V[n]) && ((_, q) = divides(lc, leading_coefficient(p)))[1]
      s = V[1:n - 1]
      p, V = V[n], vcat(s, [reduce_tail(p, s)])
      swap = true
   end
   # check if leading coefficients divide leading_coefficient of p
   while n >= 1 && (swap || ((_, q) = divides(leading_coefficient(p), leading_coefficient(V[n])))[1])
      p -= q*shift_left(V[n], length(p) - length(V[n]))
      while n >= 1 && length(V[n]) > length(p)
         n -= 1
      end
      swap = false
   end
   if n == 0 # p is smallest polynomial
      if iszero(p) # p was absorbed, yay!
         return V
      end
      return extend_ideal_basis(reverse(V), [divexact(p, canonical_unit(p))])
   end
   if n < length(V) # we made some progress
      return extend_ideal_basis(vcat(reverse(V[n+1:end]), [p]), V[1:n])
   end
   # we made no progress, use gcdx
   n = length(V)
   v = V[n]
   g, s, t = gcdx(leading_coefficient(p), leading_coefficient(v))
   r = s*p + t*shift_left(v, length(p) - length(v)) # r has leading coeff g
   q = divexact(leading_coefficient(p), g)
   p -= q*shift_left(r, length(p) - length(r))
   if length(r) == length(V[n]) # V[n] can be reduced by r and switched
      q = divexact(leading_coefficient(V[n]), g)
      r, V[n] = V[n] - q*r, r
      if n > 1
         V[n] = reduce_tail(V[n], V[1:n-1])
      end
      if length(r) > length(p)
         r, p = p, r
      end
      if length(p) == 0 # both polynomials were absorbed, yay!
         return V
      end
      if length(r) == 0 # one polynomial was absorbed, yay!
         r = p
      else
         # insert p in V
         lenp = length(p)
         n = findfirst(x->length(x) >= lenp, V)
         V = insert!(V, n, p)
      end
   else # length(r) > length(V[n])
      if length(p) == 0 # one polynomial was absorbed, yay
         return vcat(V, [reduce_tail(r, V)])
      end
      V = vcat(V, [r])
      r = p
   end
   lenr = length(r)
   n = findfirst(x->length(x) >= lenr, V)
   if n == 1 # r is the smallest polynomial
      return extend_ideal_basis(reverse(V), [divexact(r, canonical_unit(r))])
   end
   return extend_ideal_basis(vcat(reverse(V[n:end]), [r]), V[1:n - 1])
end
      
# We call an ideal over a polynomial ring over a Euclidean domain reduced if
# 1. There is only one polynomial of each degree in the ideal
# 2. The degree of polynomials in the basis increases
# 3. The leading coefficient of f_i divides that of f_{i-1} for all i
# 4. Only the final polynomial may have leading coefficient that is a unit
# 5. The polynomials are all canonicalised (divided by their canonical_unit)
# 6. The tail of each polynomial is reduced mod the other polynomials in the basis
function reduce(I::Ideal{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   if hasmethod(gcdx, Tuple{U, U})
      V = gens(I)
      # Step 1: compute V a vector of polynomials giving the same basis as I
      #         but for which 1, 2, 3 above hold
      if length(V) > 1
         D = sort(V, by=degree, rev=true)
         d = pop!(D)
         if isconstant(d)
            S = parent(d)
            d0 = constant_coefficient(d)
            while !isempty(D) && isconstant(D[1])
               di = constant_coefficient(pop!(D))
               d0 = gcd(d0, di)
            end
            d = S(d0)
         end
         V = [divexact(d, canonical_unit(d))]
         V = extend_ideal_basis(D, V)
      end
      return Ideal{T}(base_ring(I), V)
   else
      error("Not implemented")
   end
end

###############################################################################
#
#   Ideal reduction in Euclidean domain
#
###############################################################################

function reduce_euclidean(I::Ideal{T}) where T <: RingElement
end

function reduce(I::Ideal{T}) where T <: RingElement
   if hasmethod(gcdx, Tuple{T, T})
      I = reduce_euclidean(I)
   end
end

function reduce(I::Ideal{T}) where {U <: FieldElement, T <: AbstractAlgebra.PolyElem{U}}
   return reduce_euclidean(I)
end

###############################################################################
#
#   Ideal constructor
#
###############################################################################

function Ideal(R::Ring, V::Vector{T}) where T <: RingElement
   I = Ideal{elem_type(R)}(R, filter(!iszero, map(R, V)))
   return reduce(I)
end

function Ideal(R::Ring, v::T...) where T <: RingElement
   I = Ideal{elem_type(R)}(R, filter(!iszero, [map(R, v)...]))
   return reduce(I)
end

###############################################################################
#
#   IdealSet constructor
#
###############################################################################

function IdealSet(R::Ring)
   return IdealSet{elem_type(R)}(R)
end
