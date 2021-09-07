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

mutable struct lmnode{U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
   poly::Union{U, Nothing}
   up::Union{lmnode{U, V, N}, Nothing}   # out (divisible leading monomials)
   next::Union{lmnode{U, V, N}, Nothing} # extend current node so out-degree can be > 1
   lm::NTuple{N, Int}   # leading monomial as exponent vector
   lcm::NTuple{N, Int}  # lcm of lm's in tree rooted here, as exponent vector
   in_heap::Bool # whether node is in the heap
   path::Bool # used for marking paths to divisible nodes
   path2::Bool # mark paths when following existing paths

   function lmnode{U, V, N}(p::Union{U, Nothing}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
      node = new{U, V, N}(p, nothing, nothing)
      if p != nothing && !iszero(p)
         node.lm = Tuple(exponent_vector(p, 1))
         node.lcm = node.lm
      end
      node.path = false
      node.path2 = false
      node.in_heap = false
      return node::lmnode{U, V, N}
   end
end

function show(io::IO, n::lmnode)
   print(io, "Node(", n.poly, ", ", n.up, ", ", n.next, ")")
end

# heap implementation for sorting polys by lm, head = smallest

# Defined in MPoly, repeated here for reference
# heapleft(i::Int) = 2i
# heapright(i::Int) = 2i + 1
# heapparent(i::Int) = div(i, 2)

function lm_precedes(node1::lmnode{U, V, N}, node2::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
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

function lm_divides(node1::lmnode{U, V, N}, node2::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
   for i = 1:N
      if node2.lm[i] > node1.lm[i]
         return false
      end
   end
   return true
end

function lm_divides(f::lmnode{U, V, N}, i::Int, node2::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
   D = exponent_vector(f.poly, i)
   for j = 1:N
      if node2.lm[j] > D[j]
         return false
      end
   end
   return true
end

function lm_lcm(node1::lmnode{U, V, N}, node2::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
   return Tuple(max(node1.lcm[i], node2.lcm[i]) for i in 1:N)
end

function lm_divides_lcm(node1::lmnode{U, V, N}, node2::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
   for i = 1:N
      if node2.lm[i] > node1.lcm[i]
         return false
      end
   end
   return true
end

function lm_is_constant(node::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
   for i = 1:N
      if node.lm[i] != 0
         return false
      end
   end
   return true
end

function heapinsert!(heap::Vector{lmnode{U, V, N}}, node::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
@assert !iszero(node.poly)
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

function heappop!(heap::Vector{lmnode{U, V, N}}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
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
   x.up = nothing
   x.next = nothing
   x.lm = Tuple(exponent_vector(x.poly, 1))
   x.lcm = x.lm
   return x
end

function extract_gens(D::Vector{U}, node::T) where {N, U <: AbstractAlgebra.MPolyElem, V, T <: lmnode{U, V, N}}
   # depth first
   if node.up != nothing
      if !node.up.path
         extract_gens(D, node.up)
      end
   end
   if node.next != nothing
      extract_gens(D, node.next)
   end
   if node.poly != nothing
      push!(D, node.poly)
   end
   node.path = true
   return nothing
end

function extract_gens(B::Vector{T}) where {N, U <: AbstractAlgebra.MPolyElem, V, T <: lmnode{U, V, N}}
   D = Vector{U}()
   for node in B
      extract_gens(D, node)
   end
   return D
end

# divides1 = leading coeff of d divides leading coeff of b node and ** holds
# divides2 = leading coeff of b node divides leading coeff of d and ** holds
# eq       = leading monomials of b node and d are the same
# ** lm of b node divides that of lm of d
function find_divisor_nodes!(b::T, d::T)  where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N, V, T <: lmnode{U, V, N}}
   divides1 = true
   divides2 = false
   eq = false
   b.path = true
   if b.up != nothing
      found = false
      if lm_divides(d, b.up)
         found = true
         if !b.up.path
            divides1, divides2, eq = find_divisor_nodes!(b.up, d)
         end
      end
      b2 = b
      while b2 != nothing && b2.next != nothing
         if lm_divides(d, b2.next.up)
            found = true
            if !b2.next.up.path
               d1, d2, eq0 = find_divisor_nodes!(b2.next.up, d)
               divides1 &= d1
               divides2 |= d2
               eq |= eq0 
            end
         end
         b2 = b2.next
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

function tree_attach_node(b::T, d::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N, V, T <: lmnode{U, V, N}}
   if b.up != nothing
      found = false
      if b.up.path
         found = true
         if !b.up.path2
            tree_attach_node(b.up, d)
         end
      end
      b2 = b
      while b2.next != nothing
         b2 = b2.next
         if b2.up.path
            found = true
            if !b.up.path2
               tree_attach_node(b2.up, d)
            end
         end
      end
      if !found # we must attach to this node
         b2.next = lmnode{U, V, N}(nothing)
         b2.next.up = d
      end
   else
      # we must attach to this node
      b.up = d
   end
   b.path2 = true # mark path as visited
   b.lcm = lm_lcm(b, d) # update lcm for current node
end

function reduce_leading_term!(d::T, b::T, q::W) where {W <: RingElement, U <: AbstractAlgebra.MPolyElem{W}, N, V, T <: lmnode{U, V, N}}
   infl = [1 for in in 1:N]
   shift = exponent_vector(d.poly, 1) .- exponent_vector(b.poly, 1)
   d.poly -= q*inflate(b.poly, shift, infl)
   if !iszero(d.poly)
      d.lm = Tuple(exponent_vector(d.poly, 1))
      d.lcm = d.lm
   end
end

function reduce_leading_term!(d::T, b::U, q::W) where {W <: RingElement, U <: AbstractAlgebra.MPolyElem{W}, V, N, T <: lmnode{U, V, N}}
   infl = [1 for in in 1:N]
   shift = exponent_vector(d.poly, 1) .- exponent_vector(b, 1)
   d.poly -= q*inflate(b, shift, infl)
   if !iszero(d.poly)
      d.lm = Tuple(exponent_vector(d.poly, 1))
      d.lcm = d.lm
   end
end

function tree_reduce_node!(d::T, b::T, B::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   reduced = false
   if b.up != nothing
      found = false
      if b.up.path
         found = true
         if !b.up.path2
            reduced = tree_reduce_node!(d, b.up, B)
         end
      end
      b2 = b
      while !reduced && b2.next != nothing
         b2 = b2.next
         if b2.up.path
            found = true
            if !b2.up.path2
               reduced = tree_reduce_node!(d, b2.up, B)
            end
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
   b.path2 = true
   return reduced
end

function tree_clear_path(b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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
   b.path2 = false
end

function tree_remove(b::T, d::T, heap::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   if b.up != nothing
      if !b.up.in_heap
         tree_remove(b.up, d, heap)
      end
      b.up = nothing
      b2 = b
      while b2.next != nothing
         b2 = b2.next
         if !b2.up.in_heap
            tree_remove(b2.up, d, heap)
         end
         b2.up = nothing
      end
   end
   b.next = nothing
   if b !== d # do not insert b on heap if it is d
      heapinsert!(heap, b)
   end
end

function tree_evict(b::T, d::T, d1::T, heap::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   node_lcm = b.lm
   if b.up != nothing
      if lm_divides_lcm(b.up, d)
         if lm_divides(b.up, d)
            if !b.up.in_heap
               tree_remove(b.up, d1, heap)
            end
            b.up = nothing
         else
            tree_evict(b.up, d, d1, heap)
            node_lcm = max.(node_lcm, b.up.lcm)
         end
      else
         node_lcm = max.(node_lcm, b.up.lcm)
      end
      b2 = b
      while b2 != nothing && b2.next != nothing
         if lm_divides_lcm(b2.next.up, d)
            if lm_divides(b2.next.up, d)
               if !b2.next.up.in_heap
                  tree_remove(b2.next.up, d1, heap)
               end
               b2.next.up = nothing
               b2.next = b2.next.next
            else
               tree_evict(b2.next.up, d, d1, heap)
               node_lcm = max.(node_lcm, b2.next.up.lcm)
               b2 = b2.next
            end
         else
            node_lcm = max.(node_lcm, b2.next.up.lcm)
            b2 = b2.next
         end
      end
   end
   if b.up == nothing
      if b.next != nothing
         b.up = b.next.up
         b.next = b.next.next
      end
   end
   b.lcm = node_lcm
end

function reduce_leading_term!(X::Vector{T}, d::T, b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   eq = d.lm == b.lm
   if eq && ((flag, q) = divides(leading_coefficient(b.poly), leading_coefficient(d.poly)))[1]
      # reduce b by p and mark for eviction
      p = b.poly
      reduce_leading_term!(b, d.poly, q)
      push!(X, d)
      r = lmnode{U, V, N}(b.poly)
      if !iszero(r.poly)
         push!(X, r)
      end
      b.poly = p
      b.lm = Tuple(exponent_vector(b.poly, 1))
      b.lcm = b.lm
      push!(X, b)
   else
      g, s, t = gcdx(leading_coefficient(d.poly), leading_coefficient(b.poly))
      if eq
         p = s*d.poly + t*b.poly # p has leading coeff g dividing lc(d) and lc(b)
      else
         shift = exponent_vector(d.poly, 1) .- exponent_vector(b.poly, 1)
         infl = [1 for i in 1:N]
         p = s*d.poly + t*inflate(b.poly, shift, infl)
      end
      # reduce d and b by p and replace b by p (temporarily)
      q = divexact(leading_coefficient(d.poly), leading_coefficient(p))
      reduce_leading_term!(d, p, q)
      if eq
         q = divexact(leading_coefficient(b.poly), leading_coefficient(p))
         reduce_leading_term!(b, p, q)
         r = lmnode{U, V, N}(b.poly)
         b.poly = p
         if !iszero(b.poly)
            b.lm = Tuple(exponent_vector(b.poly, 1))
            b.lcm = b.lm
         end
         push!(X, d, r, b)
      else
         r = lmnode{U, V, N}(p)
         push!(X, d, r)
      end
   end
end

function tree_reduce_equal_lm!(X::Vector{T}, d::T, b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   reduced = false
   if b.up != nothing
      found = false
      if b.up.path
         found = true
         if !b.up.path2
            reduced = tree_reduce_equal_lm!(X, d, b.up)
         end
      end
      b2 = b
      while !reduced && b2.next != nothing
         b2 = b2.next
         if b2.up.path
            found = true
            if !b2.up.path2
               reduced = tree_reduce_equal_lm!(X, d, b2.up)
            end
         end
      end
      if !reduced && !found # we must check this node
         if (d.lm == b.lm)
            reduce_leading_term!(X, d, b)
            reduced = true
         end
      end
   else
      # we must check this node
      if (d.lm == b.lm)
         reduce_leading_term!(X, d, b)
         reduced = true
      end
   end
   b.path2 = true
   return reduced
end

function tree_reduce_gcdx!(X::Vector{T}, d::T, b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   reduced = false
   if b.up != nothing
      found = false
      if b.up.path
         found = true
         if !b.up.path2
            reduced = tree_reduce_gcdx!(X, d, b.up)
         end
      end
      b2 = b
      while !reduced && b2.next != nothing
         b2 = b2.next
         if b2.up.path
            found = true
            if !b2.up.path2
               reduced = tree_reduce_gcdx!(X, d, b2.up)
            end
         end
      end
      if !reduced && !found && !divides(leading_coefficient(b.poly), leading_coefficient(d.poly))[1]
         # we have reached the right node
         reduce_leading_term!(X, d, b)
         reduced = true
      end
   else
      # we have reached the right node
      if !divides(leading_coefficient(b.poly), leading_coefficient(d.poly))[1]
         reduce_leading_term!(X, d, b)
         reduced = true
      end
   end
   b.path2 = true
   return reduced
end

# Reduce the i-th term of f by the polys in tree b, the root of which divides it
function reduce_tail!(f::T, i::Int, b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   mon = Tuple(exponent_vector(f.poly, i))
   reduced_to_zero = false
   b.path = true
   if b.up != nothing
      if lm_divides(f, i, b.up) && !b.up.path
         reduced_to_zero = reduce_tail!(f, i, b.up)
      end
      b2 = b
      while !reduced_to_zero && b2 != nothing && b2.next != nothing
         if lm_divides(f, i, b2.next.up) && !b2.next.up.path
            reduced_to_zero = reduce_tail!(f, i, b2.next.up)
         end
         b2 = b2.next
      end
   end
   if !reduced_to_zero
      shift = exponent_vector(f.poly, i) .- exponent_vector(b.poly, 1)
      infl = [1 for j in 1:N]
      q = div(coeff(f.poly, i), leading_coefficient(b.poly))
      f.poly -= q*inflate(b.poly, shift, infl)
      reduced_to_zero = length(f.poly) < i || Tuple(exponent_vector(f.poly, i)) != mon
   end
   return reduced_to_zero 
end

function reduce_tail!(f::T, B::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   i = 2
   while i <= length(f.poly)
      reduced_to_zero = false
      for b in B
         if !reduced_to_zero && lm_divides(f, i, b)
            reduced_to_zero = reduce_tail!(f, i, b)
            tree_clear_path(b)
         end
      end
      if !reduced_to_zero
         i += 1
      end
   end
end

function check_tree_path(b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   if b.in_heap
      println(b)
      error("node is in heap")
   end
   if b.next != nothing && b.up == nothing
       println(b)
       error("tree node error")
   end
   if b.up != nothing
       check_tree_path(b.up)
       while b.next != nothing
          check_tree_path(b.next)
          b = b.next
       end
   else
      if b.path
         println(b)
         error("path is set")
      end
      if b.path2
         println(b2)
         error("path2 is set")
      end
   end
end

function extract_nodes(D::Vector{T}, node::T) where {U <: AbstractAlgebra.MPolyElem, V, N, T <: lmnode{U, V, N}}
   # depth first
   if node.up != nothing
      if !node.up.path
         extract_nodes(D, node.up)
      end
   end
   if node.next != nothing
      extract_nodes(D, node.next)
   end
   if node.poly != nothing
      push!(D, node)
   end
   node.path = true
   return nothing
end

function extract_nodes(B::Vector{T}) where {U <: AbstractAlgebra.MPolyElem, V, N, T <: lmnode{U, V, N}}
   D = Vector{T}()
   for node in B
      extract_nodes(D, node)
   end
   for node in B
      tree_clear_path(node)
   end
   return D
end

function basis_check(B::Vector{T}, heap::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   for b in B
      check_tree_path(b)
   end
   D = extract_nodes(B)
   for i = 1:length(D)
      for j = i + 1:length(D)
         if D[i] === D[j]
            println("B = ", B)
            println("D = ", D)
            error("Duplicate entries in basis")
         end
      end
      if iszero(D[i].poly)
         error("poly is zero")
      end
      for k = 1:length(heap)
         if D[i] === heap[k]
            println(D[i])
            error("entry in heap and basis")
         end
      end
   end
   for b in D
      if b.lm != Tuple(exponent_vector(b.poly, 1))
         println(b)
         error("lm field is incorrect")
      end
      if b.up == nothing
         if b.lm != b.lcm
            println(b)
            error("lcm doesn't match lm")
         end
      else
         node_lcm = b.lm
         node_lcm = max.(node_lcm, b.up.lcm)
         if !lm_divides(b.up, b)
            println(b)
            error("lm does not divide up")
         end
         if b.lm == b.up.lm
            if !divides(leading_coefficient(b.poly), leading_coefficient(b.up.poly))[1]
               error("leading coefficients of equal leading monomials don't divide")
               println(b)
            end
         end
         b2 = b
         while b2.next != nothing
            if !lm_divides(b2.next.up, b)
               println(b)
               error("lm does not divide next")
            end
            node_lcm = max.(node_lcm, b2.next.up.lcm)
            b2 = b2.next
         end
         if b.lcm != node_lcm
            println(b)
            println(b.lcm)
            println(node_lcm)
            error("lcm is incorrect")
         end
      end
   end
   for i = 1:length(D)
      for k = 1:length(heap)
         if lm_divides(D[i], heap[k]) && D[i].lm != heap[k].lm
            println(D[i])
            println(heap[k])
            error("entry in heap divides entry in basis")
         end
      end
   end
end

const LMNODE_DEBUG = false

function basis_insert(W::Vector{Tuple{Bool, Bool, Bool}}, X::Vector{T}, B::Vector{T}, d::T, heap::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   node_divides = false
   if lm_is_constant(d)
      for b in B
         if lm_is_constant(b)
            S = parent(d.poly)
            b.poly = S(gcd(constant_coefficient(b.poly), constant_coefficient(d.poly)))
            node_divides = true
         end
      end
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
      if LMNODE_DEBUG
         println("    ***************")
         println("    B = ", B)
         println("    W = ", W)
         println("    d = ", d)
      end
      w = true, false, false
      for v in W
         w = w[1] & v[1], w[2] | v[2], w[3] | v[3]
      end
      if node_divides
         if w[1] & !w[2] & !w[3]
            if LMNODE_DEBUG
               println("CASE 1:")
            end
            # poly can be attached to basis
            for b in B
               if b.path # d is divisible by at least one node in tree b
                  tree_attach_node(b, d)
                  tree_clear_path(b)
               end
            end
            d.poly = divexact(d.poly, canonical_unit(d.poly))
            reduce_tail!(d, B)
         elseif w[2]
            if LMNODE_DEBUG
               println("CASE 2:")
            end
            # lt of poly can be reduced
            reduced = false
            for i in 1:length(B)
               if !reduced && W[i][2] # some node of b has leading coeff dividing leading coeff of d
                  # reduce leading coeff of d
                  reduced = tree_reduce_node!(d, B[i], B)
               end
               if B[i].path
                  tree_clear_path(B[i])
               end
            end
            if !iszero(d.poly)
               # evict terms divisible by new d
               for i in length(B):-1:1
                  if lm_divides_lcm(B[i], d)
                     if lm_divides(B[i], d)
                        if !B[i].in_heap
                           tree_remove(B[i], d, heap)
                        end
                        deleteat!(B, i)
                        deleteat!(W, i)
                     else
                        tree_evict(B[i], d, d, heap)
                     end
                  end
               end
               # insert d on heap
               heapinsert!(heap, d)
            end
         elseif w[3]
            if LMNODE_DEBUG
               println("CASE 3:")
            end
            # lt of poly can be reduced
            reduced = false
            for i in 1:length(B)
               if !reduced && W[i][3] # some node of b has leading coeff dividing leading coeff of d
                  # reduce leading coeff of d
                  reduced = tree_reduce_equal_lm!(X, d, B[i])
               end
               if B[i].path
                  tree_clear_path(B[i])
               end
            end
            d1 = pop!(X)
            # evict terms divisible by old d (it will not be placed back on heap by eviction)
            for i in length(B):-1:1
               if lm_divides_lcm(B[i], d1)
                  if lm_divides(B[i], d1)
                     if !B[i].in_heap
                        tree_remove(B[i], d1, heap)
                     end
                     deleteat!(B, i)
                     deleteat!(W, i)
                  else
                     tree_evict(B[i], d1, d1, heap)
                  end
               end
            end
            while !isempty(X)
               d = pop!(X)
               if !iszero(d.poly)
                  # evict terms divisible by new d
                  for i in length(B):-1:1
                     if lm_divides_lcm(B[i], d)
                        if lm_divides(B[i], d)
                           if !B[i].in_heap
                              tree_remove(B[i], d1, heap)
                           end
                           deleteat!(B, i)
                           deleteat!(W, i)
                        else
                           tree_evict(B[i], d, d1, heap)
                        end
                     end
                  end
                  # insert d on heap
                  heapinsert!(heap, d)
               end
            end
         else
            if LMNODE_DEBUG
               println("CASE 4:")
            end
            # leading term of d can be reduced
            reduced = false
            for i in 1:length(B)
               if !reduced && !W[i][1] && B[i].path # some node of b has leading monomial dividing leading monomial of d
                  # reduce leading coeff of d
                  reduced = tree_reduce_gcdx!(X, d, B[i])
               end
               if B[i].path
                  tree_clear_path(B[i])
               end
            end
            while !isempty(X)
               d = pop!(X)
               if !iszero(d.poly)
                  # evict terms divisible by new d
                  for i in length(B):-1:1
                     if lm_divides_lcm(B[i], d)
                        if lm_divides(B[i], d)
                           if !B[i].in_heap
                              tree_remove(B[i], d, heap)
                           end
                           deleteat!(B, i)
                           deleteat!(W, i)
                        else
                           tree_evict(B[i], d, d, heap)
                        end
                     end
                  end 
                  # insert d on heap
                  heapinsert!(heap, d)
               end
            end
         end
      end
   end
   if !node_divides # d not divisible by the root of any tree
      push!(B, d) # add to basis
      push!(W, (true, false, false))
      d.poly = divexact(d.poly, canonical_unit(d.poly))
      reduce_tail!(d, B)
   end
   if LMNODE_DEBUG
      basis_check(B, heap)
   end
   return Nothing
end

function reduce(I::Ideal{U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
   if hasmethod(gcdx, Tuple{T, T})
      B = gens(I)
      # Step 1: compute B a vector of polynomials giving the same basis as I
      #         but for which 1, 2, 3 above hold
      if length(B) > 1
         # make heap
         V = ordering(parent(B[1]))
         N = nvars(base_ring(I))
         heap = Vector{lmnode{U, V, N}}()
         for v in B
            heapinsert!(heap, lmnode{U, V, N}(v))
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
         B2 = [d]
         W = [(true, false, false)]
         X = lmnode{U, V, N}[]
         # do reduction
         while !isempty(heap)
            d = heappop!(heap)
            basis_insert(W, X, B2, d, heap)
         end
         # extract polynomials from B2
         B = extract_gens(B2)
      end
      return Ideal{U}(base_ring(I), B)
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

