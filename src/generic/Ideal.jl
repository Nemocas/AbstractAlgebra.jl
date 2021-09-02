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
   poly::Union{U, Nothing}
   up::Union{lmnode{U}, Nothing}   # out (divisible leading monomials)
   next::Union{lmnode{U}, Nothing} # extend current node so out-degree can be > 1
   lm::NTuple{N, Int}   # leading monomial as exponent vector
   lcm::NTuple{N, Int}  # lcm of lm's in tree rooted here, as exponent vector
   in_heap::Bool # whether node is in the heap
   path::Bool # used for marking paths to divisible nodes

   function lmnode{U, N}(p::Union{U, Nothing}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
      node = new{U, N}(p, nothing, nothing)
      node.lm = Tuple(exponent_vector(p, 1))
      node.lcm = node.lm
      node.path = false
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
      if node2.lm[i] < node1.lm[i]
         return true
      end
      if node2.lm[i] > node1.lm[i]
         return false
      end
   end
   return true
end

function lm_divides(node1::lmnode{U, N}, node2::lmnode{U, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
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

function lm_divides_lcm(node1::lmnode{U, N}, node2::lmnode{U, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
   for i = 1:N
      if node2.lm[i] > node1.lcm[i]
         return false
      end
   end
   return true
end

function lm_is_constant(node::lmnode{U, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
   for i = 1:N
      if node.lm[i] != 0
         return false
      end
   end
   return true
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

# divides1 = leading coeff of d divides leading coeff of b node and ** holds
# divides2 = leading coeff of b node divides leading coeff of d and ** holds
# eq       = leading monomials of b node and d are the same
# ** lm of b node divides that of lm of d
function find_divisor_nodes!(b::T, d::T)  where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N, T <: lmnode{U, N}}
   divides1 = true
   divides2 = false
   eq = false
   b.path = true
   if b.up != nothing
      found = false
      if lm_divides(d, b.up)
         found = true
         divides1, divides2, eq = find_divisor_nodes!(b.up, d)
      end
      while b.next != nothing
         if lm_divides(d, b.next.up)
            found = true
            d1, d2, eq0 = find_divisor_nodes!(b.next.up, d)
            divides1 &= d1
            divides2 |= d2
            eq |= eq0 
         else
            b = b.next
         end
      end
      if !found
         divides1 = divides(leading_coefficient(b.poly), leading_coefficient(d.poly))[1]
         divides2 = divides(leading_coefficient(d.poly), leading_coefficient(b.poly))[1]
         eq = b.lm == d.lm
      end
   else
      divides1 = divides(leading_coefficient(b.poly), leading_coefficient(d.poly))[1]
      divides2 = divides(leading_coefficient(d.poly), leading_coefficient(b.poly))[1]
      eq = b.lm == d.lm
   end
   return divides1, divides2, eq
end

function tree_attach_node(b::T, d::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N, T <: lmnode{U, N}}
   if b.up != nothing
      found = false
      if b.up.path
         found = true
         tree_attach_node(b.up, d)
      end
      b2 = b
      while b2.next != nothing
         b2 = b2.next
         if b2.up.path
            found = true
            tree_attach_node(b2.up, d)
         end
      end
      if !found # we must attach to this node
         b2.next = lmnode{U, N}(nothing)
         b2.next.up = d
      end
   else
      # we must attach to this node
      b.up = d
   end
   b.path = false # clear path flag
   b.lcm = lm_lcm(b, d) # update lcm for current node
end

function reduce_leading_term!(d::T, b::T, q::V) where {V <: RingElement, U <: AbstractAlgebra.MPolyElem{V}, N, T <: lmnode{U, N}}
   infl = [1 for in in 1:N]
   shift = exponent_vector(d.poly, 1) .- exponent_vector(b.poly, 1)
   d.poly -= q*inflate(b.poly, shift, infl)
   if !iszero(d.poly)
      d.lm = Tuple(exponent_vector(d.poly, 1))
      d.lcm = d.lm
   end
end

function tree_reduce_node!(d::T, b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N, T <: lmnode{U, N}}
   reduced = false
   if b.up != nothing
      found = false
      if b.up.path
         found = true
         reduced = tree_reduce_node!(d, b.up)
      end
      b2 = b
      while !reduced && b2.next != nothing
         b2 = b2.next
         if b2.up.path
            found = true
            reduced = tree_reduce_node!(d, b2.up)
         end
      end
      if !reduced && !found # we must check this node
         if ((flag, q) = divides(leading_coefficient(d.poly), leading_coefficient(b.poly)))[1]
            reduce_leading_term!(d, b, q)
            reduced = true
         end
      end
   else
      # we must check this node
      if ((flag, q) = divides(leading_coefficient(d.poly), leading_coefficient(b.poly)))[1]
         reduce_leading_term!(d, b, q)
         reduced = true
      end
   end
   return reduced
end

function tree_clear_path(b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N, T <: lmnode{U, N}}
   if b.up != nothing
      if b.up.path
         tree_clear_path(b.up)
      end
      b2 = b
      while b2.next != nothing
         b2 = b2.next
         if b2.up.path
            tree_clear_path(b2.up)
         end
      end
   end
   b.path = false # clear path flag
end

function tree_remove(b::T, heap::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N, T <: lmnode{U, N}}
   if b.up != nothing
      if !b.up.in_heap
         tree_remove(b.up, heap)
      end
      b.up = nothing
      b2 = b
      while b2.next != nothing
         b2 = b2.next
         if !b2.up.in_heap
            tree_remove(b2.up)
         end
         b2.up = nothing
      end
   end
   b.next = nothing
   heapinsert!(heap, b)
end

function tree_evict(b::T, d::T, heap::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N, T <: lmnode{U, N}}
   node_lcm = Tuple(1 for i in 1:N)
   if b.up != nothing
      if lm_divides_lcm(b.up, d)
         if lm_divides(b.up, d) && !b.up.in_heap
            tree_remove(b.up, heap)
            b.up = nothing
         else
            tree_evict(b.up)
            node_lcm = lcm.(node_lcm, b.up)
         end
      else
         node_lcm = lcm.(node_lcm, b.up)
      end
      b2 = b
      while b2.next != nothing
         if lm_divides_lcm(b2.next.up, d)
            if lm_divides(b2.next.up, d) && !b2.next.up.in_heap
               tree_remove(b2.next.up, heap)
               b2.next.up = nothing
               b2.next = b2.next.next
            else
               tree_evict(b2.next.up)
               node_lcm = lcm.(node_lcm, b2.next.up)
            end
         else
            node_lcm = lcm.(node_lcm, b2.next.up)
         end
         b2 = b2.next
      end
   end
   if b.up == nothing
      if b.next != nothing
         b.up = b.next.up
         b.next = b.next.next
      end
   end
end

function basis_insert(W::Vector{Tuple{Bool, Bool, Bool}}, B::Vector{T}, d::T, heap::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N, T <: lmnode{U, N}}
   node_divides = false
   if lm_is_constant(d)
      # should not happen as nothing in heap can divide something in basis, except if
      # lm is the same, which implies there is nothing but constants in basis, which means
      # this term from heap should have already been merged in reduce function below
      for b in B
         if lm_is_constant(b)
            S = parent(d.poly)
            b.poly = S(gcd(constant_coefficient(b.poly), constant_coefficient(d.poly)))
            node_divides = true
         end
      end
      error("Should not happen")
   else
      for i in 1:length(B)
         b = B[i]
         if lm_divides(d, b) # d is divisible by something in tree b
            node_divides = true
            # first find all maximal divisor nodes and characterise
            W[i] = find_divisor_nodes!(b, d) # mark path to maximal nodes dividing d
         else
            W[i] = true, false, false # just so test below is not interfered with
         end
      end
      w = true, false, false
      for v in W
         w = w[1] & v[1], w[2] | v[2], w[3] | v[3]
      end
      if node_divides
         if w[1] & !w[2] & !w[3]
            # poly can be attached to basis
            d.poly = divexact(d.poly, canonical_unit(d.poly))
            for b in B
               if b.path # d is divisible by at least one node in tree b
                  tree_attach_node(b, d)
               end
            end
         elseif w[2]
            # lt of poly can be reduced
            reduced = false
            for i in 1:length(B)
               if !reduced && W[i][2] # some node of b has leading coeff dividing leading coeff of d
                  # reduce leading coeff of d
                  reduced = tree_reduce_node!(d, B[i])
               end
               if B[i].path
                  tree_clear_path(B[i])
               end
            end
            if !iszero(d.poly)
               # evict terms divisible by new d
               for i in length(B):-1:1
                  if lm_divides_lcm(B[i], d)
                     if lm_divides(B[i], d) && !B[i].in_heap
                        tree_remove(B[i], heap)
                        deleteat!(B, i)
                     else
                        tree_evict(B[i], d, heap)
                     end
                  end
               end
               # insert d on heap
               heapinsert!(heap, d)
            end
         else
            error("Not implemented")
         end
      end
   end
   if !node_divides # d not divisible by the root of any tree
      d.poly = divexact(d.poly, canonical_unit(d.poly))
      push!(B, d) # add to basis
      push!(W, (true, false, false))
   end
   return Nothing
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
         W = [(true, false, false)]
         # do reduction
         while !isempty(heap)
            d = heappop!(heap)
            basis_insert(W, B, d, heap)
         end
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

