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
#   Heap and nodes
#
###############################################################################

node_num = [0]

mutable struct lmnode{U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
   poly::Union{U, Nothing}
   up::Union{lmnode{U, V, N}, Nothing}   # out (divisible leading monomials)
   next::Union{lmnode{U, V, N}, Nothing} # extend current node so out-degree can be > 1
   equal::Union{lmnode{U, V, N}, Nothing} # to chain nodes with equal lm
   reducer::Union{lmnode{U, V, N}, Nothing} # best reducer found so far for node
   active::Bool # whether polynomial is still actively being reduced
   new_node::Bool # whether the node was added to the lattice since last time
   settled::Int # 0 = newly reduced, 1 = no reducers remove lc, 2 = none reduce lc, 3 = no reducers
   lm::NTuple{N, Int}   # leading monomial as exponent vector
   lcm::NTuple{N, Int}  # lcm of lm's in tree rooted here, as exponent vector
   in_heap::Bool # whether node is in the heap
   path::Bool # used for marking paths to divisible nodes
   path2::Bool # mark paths when following existing paths
   num::Int
   size::Float64

   function lmnode{U, V, N}(p::Union{U, Nothing}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N}
      node = new{U, V, N}(p, nothing, nothing, nothing, nothing, true, true, 0)
      node.size = 0.0
      if p != nothing && !iszero(p)
         if leading_coefficient(p) < 0
            p = -p
         end
         node.lm = Tuple(exponent_vector(p, 1))
         node.lcm = node.lm
         node.size = reducer_size(node)
      end
      node.equal = node
      node.path = false
      node.path2 = false
      node.in_heap = false
      node_num[] += 1
      node.num = node_num[]
      return node::lmnode{U, V, N}
   end
end

function show_inner(io::IO, n::Nothing)
   print(io, "nothing")
end

poly_level = 2

function print_poly(io::IO, p::MPolyElem)
   if poly_level == 1
      print(io, leading_term(p))
      if length(p) > 1
         print(io, "+...")
      end
   else
      print(io, p)
   end
end

function show_inner(io::IO, n::lmnode)
   if !n.path
      print(io, "Node(p", n.num)
      if n.active
         print(io, "(Y)")
      else
         print(io, "(N)")
      end
      print(io, "=")
      print_poly(io, n.poly)
      print(io, ", ")
      show_inner(io, n.up)
      print(io, ", ")
      if n.next == nothing
         print(io, "nothing")
      else
         if n.next.next != nothing
            print(io, "n", n.num, ":[")
         end
         n2 = n.next
         while n2 != nothing
            show_inner(io, n2.up)
            n2 = n2.next
            if n2 != nothing
               print(io, ", ")
            end
         end
         if n.next.next != nothing
            print(io, "]")
         end
      end
      print(io, ", ")
      if n.equal == n
         print(io, "nothing")
      else
         if n.equal.equal != n
            print(io, "e", n.num, ":[")
         end
         n2 = n.equal
         while n2 != n
            print(io, "e", n2.num)
            if n2.active
               print(io, "(Y)")
            else
               print(io, "(N)")
            end
            print(io, "=")
            print_poly(io, n2.poly)
            n2 = n2.equal
            if n2 != n
               print(io, ", ")
            end
         end
         if n.equal.equal != n
            print(io, "]")
         end
      end
      print(io, ")")
   else
      print(io, "N", n.num)
   end
   n.path = true
end

function show(io::IO, n::lmnode)
   show_inner(io, n)
   clear_path(n)
end

# heap implementation for sorting polys by lm, head = smallest

# Defined in MPoly, repeated here for reference
# heapleft(i::Int) = 2i
# heapright(i::Int) = 2i + 1
# heapparent(i::Int) = div(i, 2)

function lm_precedes(node1::lmnode{U, :lex, N}, node2::lmnode{U, :lex, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
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

function lm_precedes(node1::lmnode{U, :deglex, N}, node2::lmnode{U, :deglex, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
   s1 = sum(node1.lm)
   s2 = sum(node2.lm)
   if s2 < s1
      return true
   end
   if s2 > s1
      return false
   end
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

function lm_precedes(node1::lmnode{U, :degrevlex, N}, node2::lmnode{U, :degrevlex, N}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N}
   s1 = sum(node1.lm)
   s2 = sum(node2.lm)
   if s2 < s1
      return true
   end
   if s2 > s1
      return false
   end
   for i = N:-1:1
      if node2.lm[i] > node1.lm[i]
         return true
      end
      if node2.lm[i] < node1.lm[i]
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

###############################################################################
#
#   Ideal reduction for mulivariates over Euclidean domain (old implementation)
#
###############################################################################

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

function tree_attach_node(X::Vector{T}, b::T, d::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, N, V, T <: lmnode{U, V, N}}
   if b.up != nothing
      found = false
      if b.up.path
         found = true
         if !b.up.path2
            tree_attach_node(X, b.up, d)
         end
      end
      b2 = b
      while b2.next != nothing
         b2 = b2.next
         if b2.up.path
            found = true
            if !b.up.path2
               tree_attach_node(X, b2.up, d)
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

const LMNODE_DEBUG = true

function basis_insert_old(W::Vector{Tuple{Bool, Bool, Bool}}, X::Vector{T}, B::Vector{T}, d::T, heap::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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
                  tree_attach_node(X, b, d)
                  tree_clear_path(b)
               end
            end
            d.poly = divexact(d.poly, canonical_unit(d.poly))
            reduce_tail!(d, B)
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

###############################################################################
#
#   Ideal reduction for multivariates over Euclidean domain (new implementation)
#
###############################################################################

function compute_spoly(f::T, g::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   fc = leading_coefficient(f.poly)
   gc = leading_coefficient(g.poly)
   c = lcm(fc, gc)
   llcm = max.(f.lm, g.lm)
   infl = [1 for i in 1:N]
   shiftf = llcm .- exponent_vector(f.poly, 1)
   shiftg = llcm .- exponent_vector(g.poly, 1)
   s = divexact(c, fc)*inflate(f.poly, shiftf, infl) - divexact(c, gc)*inflate(g.poly, shiftg, infl)
   return lmnode{U, V, N}(s)
end

function compute_gpoly(f::T, g::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   fc = leading_coefficient(f.poly)
   gc = leading_coefficient(g.poly)
   _, s, t = gcdx(fc, gc)
   llcm = max.(f.lm, g.lm)
   infl = [1 for i in 1:N]
   shiftf = llcm .- exponent_vector(f.poly, 1)
   shiftg = llcm .- exponent_vector(g.poly, 1)
   g = s*inflate(f.poly, shiftf, infl) + t*inflate(g.poly, shiftg, infl)
   return lmnode{U, V, N}(g)
end

function smod(c::T, h::T) where T <: RingElement
   if h < 0
      h = abs(h)
   end
   r = AbstractAlgebra.mod(c, h)
   if r >= ((h + 1) >> 1)
      r -= h
   end
   return r
end

# heuristic for size of reducer polynomials (smaller is better), used to sort
# potential reducers
function reducer_size(f::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   if f.size != 0.0
      return f.size
   end
   s = 0.0
   # heuristic really punishes polys with lots of terms
   len = length(f.poly)
   j = len
   for i = 1:len
      c = coeff(f.poly, i)
      s += Base.log(ndigits(c; base=2))*j*j
      j -= 1      
   end
   return s
end

function leading_coefficient_size(f::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   return abs(leading_coefficient(f.poly))
end

function reduce_by_reducer(S::Vector{T}, H::Vector{T}, b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   reduced = false
   if b.reducer != nothing
      c = leading_coefficient(b.poly)
      h = leading_coefficient(b.reducer.poly)
      sign = false
      if h < 0
         h = -h
         sign = true
      end
      q, r = AbstractAlgebra.divrem(c, h)
      coeff_divides = false
      if iszero(r)
         coeff_divides = true
      else
         if r >= ((h + 1) >> 1)
            q += 1
         end
      end
      if iszero(q) # need gcd polynomial
         if sign
            h = -h
         end
         g, s, t = gcdx(c, h)
         infl = [1 for in in 1:N]
         shift = exponent_vector(b.poly, 1) .- exponent_vector(b.reducer.poly, 1)
         d = s*b.poly + t*inflate(b.reducer.poly, shift, infl)
         q = divexact(c, g)
         p = b.poly - q*d # reduce b by d (the gcd poly)
         b.poly = d # replace b with gcd poly
         b.size = 0.0 # update size
         b.size = reducer_size(b)
         b.new_node = true
         if length(b.poly) == 1 # special case, must recompute s-polys
            push!(S, b)
         end
         if !iszero(p)
            push!(H, lmnode{U, V, N}(p))
         end
      else
         if sign
            q = -q
         end
         infl = [1 for in in 1:N]
         shift = exponent_vector(b.poly, 1) .- exponent_vector(b.reducer.poly, 1)
         d = b.poly - q*inflate(b.reducer.poly, shift, infl)
         if coeff_divides
            # leading coeff of reducer divides leading coeff of b.poly
            if !iszero(d) # we have a fragment
               push!(H, lmnode{U, V, N}(d))
            end
            b.active = false
         else
            # case 2: leading coeff is reduced
            b.poly = d
            b.size = 0.0 # update size
            b.size = reducer_size(b)
            b.new_node = true
            if length(b.poly) == 1 # special case, must recompute s-polys
               push!(S, b)
            end
         end
      end
      reduced = true
   end
   b.reducer = nothing
   return reduced
end

# reduce nodes in a tree depth first
# X is used to sort equal nodes
# H is used for fragments
# returns true if some reduction actually occurred
function reduce_nodes(S::Vector{T}, H::Vector{T}, b::T, X::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   b.path = true
   # depth first
   reduced = false
   if b.up != nothing
      if !b.up.path
         reduced |= reduce_nodes(S, H, b.up, X)
      end
      b2 = b
      while b2.next != nothing
         b2 = b2.next
         if !b2.up.path
            reduced |= reduce_nodes(S, H, b2.up, X)
         end
      end
   end
   # reduce nodes equal to b
   # first sort nodes by leading coefficient
   push!(X, b)
   b2 = b
   while b2.equal != b
      b2 = b2.equal
      index = searchsortedfirst(X, b2, by=leading_coefficient_size, rev=true) # find index at which to insert x
      insert!(X, index, b2)
   end
   filter!(x->x.active, X)
   i = 1
   len = length(X)
   while i < len
      c1 = leading_coefficient(X[i].poly)
      j = i + 1
      while j <= len && (leading_coefficient(X[j].poly) == c1 ||
                        leading_coefficient(X[j].poly) == -c1)
         if X[j].poly == X[i].poly || X[j] == -X[i].poly
            X[j].active = false
            deleteat!(X, j)
            reduced = true
            j -= 1
            len -= 1
         end
         j += 1
      end
      i += 1
   end
   # do reductions
   for b2 in X
      reduced |= reduce_by_reducer(S, H, b2)
   end
   # empty X
   while !isempty(X)
      pop!(X)
   end
   return reduced
end

# actually performs the reductions which have been attached to nodes by best_reducer
# X is used to sort equal nodes
# H is used for fragments
# returns true if some reduction actually occurred
function reduce_nodes(S::Vector{T}, H::Vector{T}, B::Vector{T}, X::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   reduced = false
   for b in B
      if !b.path
         reduced |= reduce_nodes(S, H, b, X)
      end
   end
   return reduced
end

function find_best_divides(b::T, X::Vector{T}, best::Union{T, Nothing}, best_divides::Bool) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   c = leading_coefficient(b.poly)
   best_size = best == nothing ? 0.0 : reducer_size(best)
   for i = length(X):-1:1
      if X[i] != b # poly can't be reduced by itself
         h = leading_coefficient(X[i].poly)
         if divides(c, h)[1]
            usable = true
            if X[i].lm == b.lm
               x2 = X[i]
               while x2 != nothing
                  if x2 == b
                     usable = false
                  end
                  x2 = x2.reducer
               end
            end
            if usable && ((c != h && c != -h) || X[i].active)
               if best == nothing || !best_divides
                  best = X[i]
                  best_size = reducer_size(X[i])
                  best_divides = true
               else
                  red_size = reducer_size(X[i])
                  if red_size < best_size
                     best = X[i]
                     best_size = red_size
                     best_divides = true
                  end
               end
            end
         end
      end
   end
   return best, best_divides
end

function find_best_reduces(b::T, X::Vector{T}, best::Union{T, Nothing}, best_reduces::Bool) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   c = leading_coefficient(b.poly)
   best_size = best == nothing ? 0.0 : reducer_size(best)
   for i = length(X):-1:1
      if X[i] != b # poly can't be reduced by itself
         h = smod(c, leading_coefficient(X[i].poly))
         if h != c && h != 0 && !divides(leading_coefficient(X[i].poly), c)[1]
            usable = true
            if X[i].lm == b.lm
               x2 = X[i]
               while x2 != nothing
                  if x2 == b
                     usable = false
                  end
                  x2 = x2.reducer
               end
            end
            if usable
               if best == nothing || !best_reduces
                  best = X[i]
                  best_size = reducer_size(X[i])
                  best_reduces = true
               else
                  red_size = reducer_size(X[i])
                  if red_size < best_size
                     best = X[i]
                     best_size = red_size
                     best_reduces = true
                  end
               end
            end
         end
      end
   end
   return best, best_reduces
end

function find_best_gcd(b::T, X::Vector{T}, best::Union{T, Nothing}, best_is_gcd::Bool) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   c = leading_coefficient(b.poly)
   best_size = best == nothing ? 0.0 : reducer_size(best)
   for i = length(X):-1:1
      if X[i] != b # poly can't be reduced by itself
         if !divides(leading_coefficient(X[i].poly), c)[1] &&
            !divides(c, leading_coefficient(X[i].poly))[1]
            usable = true
            if X[i].lm == b.lm
               x2 = X[i]
               while x2 != nothing
                  if x2 == b
                     usable = false
                  end
                  x2 = x2.reducer
               end
            end
            if usable
               if best == nothing || !best_is_gcd
                  best = X[i]
                  best_size = reducer_size(X[i])
                  best_is_gcd = true
               else
                  red_size = reducer_size(X[i])
                  if red_size < best_size
                     best = X[i]
                     best_size = red_size
                     best_is_gcd = true
                  end
               end
            end
         end
      end
   end
   return best, best_is_gcd
end

# used for heuristics, to select one kind of reducer over another preferentially
# current strategy is to prefer reducer which gives degree reduction over
# everything else, but this is probably wrong based on what Singular does
# function will return best reducer in sorted list X for reducing b
function find_best_reducer(b::T, X::Vector{T}, Xnew::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   # prevent nodes that are not active from being reduced
   if !b.active
      return nothing
   end
   c = leading_coefficient(b.poly)
   best = b.reducer
   # 1. check for leading coefficient that divides c
   best_divides = best == nothing ? false : divides(c, leading_coefficient(best.poly))[1]
   if b.settled == 0
      best, best_divides = find_best_divides(b, X, best, best_divides)
   end
   best, best_divides = find_best_divides(b, Xnew, best, best_divides)
   if best_divides
      b.settled = 0
      return best
   else
      b.settled = 1
   end
   # 2. check for leading coefficient that reduces c
   best_reduces = best == nothing ? false : smod(c, leading_coefficient(best.poly)) != c
   if b.settled <= 1
      best, best_reduces = find_best_reduces(b, X, best, best_reduces)
   end
   best, best_reduces = find_best_reduces(b, Xnew, best, best_reduces)
   if best_reduces
      return best
   else
      b.settled = 2
   end
   # 3. check for leading coefficient that is not divisible by c (gcd poly possible)
   best_is_gcd = best == nothing ? false : !divides(leading_coefficient(best.poly), c)[1]
   if b.settled <= 2
      best, best_is_gcd = find_best_gcd(b, X, best, best_is_gcd)
   end
   best, best_is_gcd = find_best_gcd(b, Xnew, best, best_is_gcd)
   if best_is_gcd
      return best
   else
      b.settled = 3 # no reducer found
   end
   return best
end

# attach best reducer to each node in tree
# see best_reducer description below for details
function best_reducer(b::T, X::Vector{T}, Xnew::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   # insert nodes equal to b in X
   b2 = b
   oldlen = length(X)
   newlen = length(Xnew)
   bend = nothing
   while b2 != bend
      if b2.new_node
         push!(Xnew, b2)
      else
         push!(X, b2)
      end
      b2 = b2.equal
      bend = b
   end
   # set best reducer of B for each node equal to b
   b.reducer = find_best_reducer(b, X, Xnew)
   b2 = b.equal
   while b2 != b
      b2.reducer = find_best_reducer(b2, X, Xnew)
      b2 = b2.equal
   end
   # recurse to b.up, b.next etc
   if b.up != nothing
      best_reducer(b.up, X, Xnew)
      b2 = b
      while b2.next != nothing
         b2 = b2.next
         best_reducer(b2.up, X, Xnew)
      end
   end
   # remove nodes equal to b from X and Xnew
   while length(X) != oldlen
      pop!(X)
   end
   while length(Xnew) != newlen
      d = pop!(Xnew)
      d.new_node = false
   end
end

# Setup for one round of reduction of leading terms by others
# Each node will be reduced mod the leading coeff of the best node
# whose leading term divides it and if there are none, by the
# best node whose leading monomial divides it and if none exists
# gcd polynomials will be created if possible
# The reductions are not actually done here; the best reducer
# is just attached to the node
# The best reducers along the current path so far are given by the
# ordered array X unless they are for new nodes since last time, in
# which case they will go in Xnew
function best_reducer(B::Vector{T}, X::Vector{T}, Xnew::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   for b in B
      best_reducer(b, X, Xnew)
   end
end

# insert fragments into basis
function insert_fragments(S::Vector{T}, B::Vector{T}, H::Vector{T}, bound::NTuple{N, Int}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   sort!(H, by=x->sum(max.(x.lm, bound) .- bound), rev=true)
   while !isempty(H)
      if max.(H[end].lm, bound) != bound
         break
      end
      d = pop!(H)
      basis_insert(S, B, d)
   end
end

# insert a node into the given branch of basis
function basis_insert(b::T, d::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   if lm_divides(b, d) # divides in both directions, so equal
      if d.active # avoid inserting polynomials that are a constant multiple of existing ones
         found = false
         c = leading_coefficient(d.poly)
         if length(d.poly) == length(b.poly)
            h = leading_coefficient(b.poly)
            flag, q = divides(c, h)
            if flag
               found = d.poly == q*b.poly
            end
         end
         d2 = b.equal
         while !found && d2 != b
            if length(d.poly) == length(d2.poly)
               h = leading_coefficient(d2.poly)
               flag, q = divides(c, h)
               if flag
                  found = d.poly == q*d2.poly
               end
            end
            d2 = d2.equal
         end
         if !found
            d.equal = b.equal
            b.equal = d
         else
            d.active = false
         end
      end
   else # not equal to current node
      if b.up != nothing
         flag = lm_divides(d, b.up)
         if flag && !b.up.path
            basis_insert(b.up, d)
         end
         b2 = b
         while b2.next != nothing
            b2 = b2.next
            if lm_divides(d, b2.up)
               flag = true
               if !b2.up.path
                  basis_insert(b2.up, d)
               end
            end
         end
         if !flag # not inserted yet
            if b != d && d.active && d.settled != -1
               n = lmnode{U, V, N}(nothing)
               n.up = d
               n.next = b.next
               b.next = n
            end
         end
      elseif d.active && d.settled != -1
         b.up = d
      end
   end
   compute_lcm(b)
   b.path = true
end

# insert a node into the basis
function basis_insert(S::Vector{T}, B::Vector{T}, d::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   insert_links(d, B)
   inserted = false
   for b in B
      if lm_divides(d, b)
         inserted = true
         basis_insert(b, d)
      end
   end
   if inserted == false # was not inserted
      push!(B, d)
   end
   clear_path(B)
   if d.active
      push!(S, d) # store for later computation of spolys
      d.settled = 0
   end
end

function find_base(d::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   if d.up != nothing
      return d
   end
   d2 = d
   while d2.equal != d
      d2 = d2.equal
      if d2.up != nothing
         return d2
      end
   end
   # declare d to be base
   return d
end

# return true if d has a connection up to b
function connected(b::T, d::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   if !lm_divides(b, d)
      return false
   end
   if d.up != nothing
      d2 = d
      while d2 != nothing
         if b == d2.up
            return true
         end
         if b.up == nothing
            b2 = b.equal
            while b2 != b
               if b2 == d2.up
                  return true
               end
               b2 = b2.equal
            end
         end
         d2 = d2.next
      end
      d2 = d
      while d2 != nothing
         if connected(b, d2.up)
            return true
         end
         d2 = d2.next
      end
   end
   return false
end

# compute spolys of everything in tree b with d with exceptions
# given in compute_spolys function below
function compute_spolys(S::Vector{T}, b::T, d::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   # depth first
   if b.up != nothing
      if !b.up.path
         compute_spolys(S, b.up, d)
      end
   end
   b2 = b
   while b2.next != nothing
      b2 = b2.next
      if !b2.up.path
         compute_spolys(S, b2.up, d)
      end
   end
   n = b
   nend = nothing
   while n != nend
      if n.active
         if n != d && !lm_divides(n, d) && !lm_divides(d, n)
            s = compute_spoly(n, d)
            if !divides(s.poly, d.poly)[1] && !divides(s.poly, n.poly)[1]
               push!(S, s)
            end
            g = compute_gpoly(n, d)
            push!(S, g)
            break
         end
         dbase = find_base(d)
         if connected(b, dbase) || connected(dbase, b)
            s = compute_spoly(n, d)
            if !divides(s.poly, d.poly)[1] && !divides(s.poly, n.poly)[1]
               push!(S, s)
            end
            break
         end 
      end
      n = n.equal
      nend = b
   end
   b.path = true
   return nothing
end

# compute spolys of all polys in basis with d, except d itself
# and polys whose lc divides that of d or which that of d divides
# unless directly connected to d
# s-polys are collected in S (not fragments to avoid too many
# s-polys being added in a single round)
# they are later moved to fragments
function compute_spolys(S::Vector{T}, B::Vector{T}, d::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   for b in B
      compute_spolys(S, b, d)
   end
   clear_path(B)
end

# insert link from d to nodes it divides
function insert_links(d::T, b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   remove = false
   if lm_divides(b, d) # we have arrived
      if !lm_divides(d, b) # make sure node is not equal
         if d.settled != -1
            if d.up == nothing
               d.up = b
            else
               n = lmnode{U, V, N}(nothing)
               n.up = b
               n.next = d.next
               d.next = n
            end
         end
         remove = true
      else
         b2 = b
         nend = nothing
         while b2 != nend
            if d.poly == b2.poly || d.poly == -b2.poly
               d.active = false
            end
            b2 = b2.equal
            nend = b
         end
         d.settled = -1 # flag as equal
         d.up = nothing
         d.next = nothing
      end
   else
      if b.up != nothing
         if lm_divides_lcm(b.up, d) && !b.up.path
            insert_links(d, b.up)
         end
         b2 = b
         while b2.next != nothing
            b2 = b2.next
            if lm_divides_lcm(b2.up, d) && !b2.up.path
               insert_links(d, b2.up)
            end   
         end
      end
   end
   b.path = true
   return remove
end

# add links from node d to existing nodes in basis that it divides
function insert_links(d::T, B::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   i = 1
   len = length(B)
   while i <= len
      b = B[i]
      if lm_divides_lcm(b, d)
         if insert_links(d, b)
            deleteat!(B, i) # node is divisible by new node
            len -= 1
            i -= 1
         end
      end
      i += 1
   end
   compute_lcm(d)
   clear_path(B)
end

# (Re)compute the lcm for node b (assuming all nodes above it have correct lcm)
function compute_lcm(b::T)  where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   b.lcm = b.lm
   if b.up != nothing
      b.lcm = max.(b.lcm, b.up.lcm)
   end
   b2 = b
   while b2.next != nothing
      b2 = b2.next
      b.lcm = max.(b.lcm, b2.up.lcm)
   end
end

# clear all path flags in given branch
function clear_path(b::T) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   if b.up != nothing
      if b.up.path
         clear_path(b.up)
      end
      b2 = b
      while b2.next != nothing
         b2 = b2.next
         if b2.up.path
            clear_path(b2.up)
         end
      end
   end
   b.path = false # clear path flag
end

# clear all path flags in basis
function clear_path(B::Vector{T}) where {U <: AbstractAlgebra.MPolyElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   for b in B
      clear_path(b)
   end
end

# extract generators from branch and store in array D
function extract_gens(D::Vector{U}, node::T) where {N, U <: AbstractAlgebra.MPolyElem, V, T <: lmnode{U, V, N}}
   # depth first
   if node.up != nothing
      if !node.up.path
         extract_gens(D, node.up)
      end
   end
   n = node
   while n.next != nothing
      n = n.next
      if !n.up.path
         extract_gens(D, n.up)
      end
   end
   if node.poly != nothing
      if node.active
         push!(D, node.poly)
      end
   end
   n = node
   while n.equal != node
      if n.equal.active
         push!(D, n.equal.poly)
      end
      n = n.equal
   end
   node.path = true
   return nothing
end

# extract all generators from basis into array D which is returned
function extract_gens(B::Vector{T}) where {N, U <: AbstractAlgebra.MPolyElem, V, T <: lmnode{U, V, N}}
   D = Vector{U}()
   for node in B
      extract_gens(D, node)
   end
   clear_path(B)
   return D
end

function generate_spolys(S::Vector{T}, B::Vector{T}, S2::Vector{T}) where {N, U <: AbstractAlgebra.MPolyElem, V, T <: lmnode{U, V, N}}
   i = 1
   len = length(S2)
   while i <= len
      s = S2[i]
      if !s.active
         deleteat!(S2, i)
         i -= 1
         len -= 1
      elseif s.settled >= 3 # only generate spolys for polys without level 3 reducers
         compute_spolys(S, B, s)
         deleteat!(S2, i)
         i -= 1
         len -= 1
      end
      i += 1
   end
end

# main reduction routine
function reduce(I::Ideal{U}) where {T <: RingElement, U <: AbstractAlgebra.MPolyElem{T}}
   node_num[] = 0
   if hasmethod(gcdx, Tuple{T, T})
      B = gens(I)
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
               h = heappop!(heap)
               di = constant_coefficient(h.poly)
               d0 = gcd(d0, di)
            end
            d.poly = S(d0)
         end
         bound = d.lm
         B2 = [d]
         X = lmnode{U, V, N}[] # fragments
         S = lmnode{U, V, N}[] # spolys
         S2 = lmnode{U, V, N}[] # nodes for which we have not computed spolys
         # insert everything in tree
         H = Vector{lmnode{U, V, N}}()
         H2 = Vector{lmnode{U, V, N}}()
         while !isempty(heap)
            d = heappop!(heap)
            bound = max.(bound, d.lm)
            basis_insert(S2, B2, d)
         end
         # do reduction
         G = Vector{U}()
         X = Vector{lmnode{U, V, N}}()
         X2 = Vector{lmnode{U, V, N}}()
         X2new = Vector{lmnode{U, V, N}}()
#println("0: B2 = ", B2)
#println("S = ", S)
#println("S2 = ", S2)
#println("H = ", H)
#println("")
#readline(stdin)
         while true
            reduction_occurs = true
            while reduction_occurs
#println("B2 = ", B2)
#readline(stdin)
               # attach best reducers to nodes
               best_reducer(B2, X2, X2new)
               # do reductions
#println("1: B2 = ", B2)
#readline(stdin)
               reduction_occurs = reduce_nodes(S2, H, B2, X)
               clear_path(B2)
               # move s-polys from S to fragments
#println("2: B2 = ", B2)
#println("")
#readline(stdin)
               while !isempty(S)
                  d = pop!(S)
                  if !iszero(d.poly)
                     push!(H, d)
                  end
               end
               # reduce contents of H if no reduction has occurred
               if !isempty(H)
                  if !reduction_occurs
                     G = extract_gens(B2)
                  end
                  if !isempty(G)
                     while !isempty(H)
                        d = pop!(H)
                        q, p = Base.divrem(d.poly, G)
                        if !iszero(p)
                           pnode = lmnode{U, V, N}(p)
                           push!(H2, pnode)
                        end
                     end
                     H, H2 = H2, H
                  end
               end
               # insert fragments (including s-polys)
               insert_fragments(S2, B2, H, bound)
#println("3: B2 = ", B2)
#println("S = ", S)
#println("S2 = ", S2)
#println("H = ", H)
#println("")
#readline(stdin)
               if !reduction_occurs
                  generate_spolys(S, B2, S2)
               end
#println("4: S = ", S)
#println("4: S2 = ", S2)
#readline(stdin)
            end
#println("B2 = ", B2)
#println("S2 = ", S2)
#println("S = ", S)
#println("")
#readline(stdin)
            if isempty(H)
#println("S2 = ", S2)
#println("S = ", S)               
               if isempty(S) && isempty(S2)
                  break
               end
            else
               bound = max.(bound, H[end].lm)
            end
         end
         # extract polynomials from B2
         B = extract_gens(B2)
      end
      B = [divexact(d, canonical_unit(d)) for d in B]
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

const IDEAL_UNIV_DEBUG = false

function check_d(D::Vector{T}, m::Int) where {T <: AbstractAlgebra.PolyElem{<:RingElement}}
   if !isempty(D)
      len = length(D[1])
      for d in D
         if length(d) > len
            error(m, ": Polynomials not in order in D")
         end
         len = length(d)
      end
   end
end

function check_v(V::Vector{T}, m::Int) where {T <: AbstractAlgebra.PolyElem{<:RingElement}}
   len = length(V[1]) - 1
   c = 2*leading_coefficient(V[1])
   for v in V
      if length(v) <= len
         error(m, ": Polynomials not in order in V")
      end
      len = length(v)
   end
   n = length(V)
   if n >= 2
      if isunit(leading_coefficient(V[n])) && isunit(leading_coefficient(V[n - 1]))
         error(m, ": Two unit leading coefficients")
      end
   end
   for v in V
      if !divides(c, leading_coefficient(v))[1]
         error(m, ": Leading coefficients do not divide in V")
      end
      if divides(leading_coefficient(v), c)[1]
         error(m, ": Associate leading coefficients in V")
      end
      c = leading_coefficient(v)
   end
end

function mysize(f::T) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   s = 0.0
   j = 0
   # heuristic really punishes polys with lots of terms
   for i = 1:length(f)
      c = coeff(f, i - 1)
      if !iszero(c)
         j += 1
         s += Base.log(ndigits(coeff(f, i - 1); base=2))*j*j
      end
   end
   return s
end

function mysize2(f::T) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   return (length(f) << 40) + mysize(f)
end

function mypush!(H::Vector{T}, f::T) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   index = searchsortedfirst(H, f, by=mysize, rev=true) # find index at which to insert x
   insert!(H, index, f)
end

# Extend the basis V of polynomials satisfying 1-6 below by the
# polynomials in D, all of which have degree at least that of those in
# V and such that the degrees of the polynomials in D is *decreasing*.
function extend_ideal_basis(D::Vector{T}, V::Vector{T}, H::Vector{T}, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   sgen = 1
   while !isempty(D)
      d = pop!(D)
      V = extend_ideal_basis(D, d, V, H, res)
      if IDEAL_UNIV_DEBUG
         check_d(D, 1)
         check_v(V, 1)
      end
      if !isempty(H)
         V = insert_fragments(D, V, H, res)
         if IDEAL_UNIV_DEBUG
            check_d(D, 2)
            check_v(V, 2)
         end
      end
   end
   return V
end

function reduce_by_resultant(f::T, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   res12 = ((res + 1) >> 1)
   for i = 1:length(f)
      c = coeff(f, i - 1)
      q, r = AbstractAlgebra.divrem(c, res)
      if r >= res12
         r -= res
      end
      setcoeff!(f, i - 1, r)
   end
   return divexact(f, canonical_unit(f))
end

function reduce_tail(f::T, V::AbstractVector{T}, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   p = iszero(res) ? f : reduce_by_resultant(f, res)
   p = divexact(p, canonical_unit(p))
   i = length(V)
   n = length(p) - 1
   while n > 0
      c = coeff(p, n - 1)
      while n > 0 && iszero(c)
         n -= 1
         if n > 0
            c = coeff(p, n - 1)
         end
      end
      if n > 0
         while i > 0 && length(V[i]) > n
            i -= 1
         end
         if i != 0 && !iszero(c)
            h = leading_coefficient(V[i]) # should be nonnegative
            h < 0 && error("h must be positive")
            q, r = AbstractAlgebra.divrem(c, h)
            if r >= ((h + 1) >> 1)
               q += 1
            end 
            u = shift_left(V[i], n - length(V[i]))
            p -= q*u
            if !iszero(res)
               p = reduce_by_resultant(p, res)
            end
         end
         n -= 1
      end
   end
   return p
end

function reduce(p::T, V::Vector{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   i = length(V)
   n = length(p)
   while n > 0
      c = coeff(p, n - 1)
      while n > 0 && iszero(c)
         n -= 1
         if n > 0
            c = coeff(p, n - 1)
         end
      end
      if n > 0
         while i > 0 && length(V[i]) > n
            i -= 1
         end
         if i != 0 && !iszero(c)
            h = leading_coefficient(V[i]) # should be nonnegative
            q, r = AbstractAlgebra.divrem(c, h)
            if r >= ((h + 1) >> 1)
               q += 1
            end 
            u = shift_left(V[i], n - length(V[i]))
            p -= q*u        
         end
         n = min(n - 1, length(p))
      end
   end
   return p
end

# An S-polynomial of f and g is a linear combination which makes the leading
# coefficients cancel. E.g. if p = 3x^2, q = 2x then S-poly(p, q) = 2*p - 3x*q
function insert_spoly(V::Vector{T}, H::Vector{T}, n::Int, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   p1 = V[n]
   if n > 1
      p2 = V[n - 1]
      c = divexact(leading_coefficient(p2), leading_coefficient(p1))
      s = c*p1 - shift_left(p2, length(p1) - length(p2))
      s = reduce_by_resultant(s, res)
      if !iszero(s)
         mypush!(H, s)
      end
   end
   if n < length(V)
      p2 = V[n + 1]
      c = divexact(leading_coefficient(p1), leading_coefficient(p2))
      s = c*p2 - shift_left(p1, length(p2) - length(p1))
      s = reduce_by_resultant(s, res)
      if !iszero(s)
         mypush!(H, s)
      end
   end
end

function insert_fragments(D::Vector{T}, V::Vector{T}, H::Vector{T}, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   max_len = 0
   for h in H
      max_len = max(max_len, length(h))
   end
   while !isempty(H)
      n = length(V)
      len = length(V[n])
      if max_len >= len # may be fragments that are too long
         V2 = T[]
         i = 1
         j = 1
         max_length = 0
         while i <= length(H)
            if length(H[i]) >= len # fragment is too long
               push!(V2, H[i])
            else
               H[j] = H[i]
               max_length = max(max_length, length(H[j]))
               j += 1
            end
            i += 1
         end
         while length(H) >= j
            pop!(H)
         end
         if length(V2) > 0
            # put long fragments back in D
            sort!(V2, by=mysize2, rev=true)
            append!(D, V2)
         end
      end
      if !isempty(H)
         p = pop!(H)
         while n >= 1 && length(V[n]) > length(p)
            n -= 1
         end
         # Reduce p by leading coefficients that divide
         while n >= 1 && ((_, q) = divides(leading_coefficient(p), leading_coefficient(V[n])))[1]
            best_v = V[n]
            best_size = mysize(V[n])
            best_q = q
            n1 = n - 1
            while n1 >= 1 && ((_, q) = divides(leading_coefficient(p), leading_coefficient(V[n1])))[1]
               new_size = mysize(V[n1])
               if new_size < best_size
                  best_v = V[n1]
                  best_size = new_size
                  best_q = q
               end
               n1 -= 1
            end
            p -= best_q*shift_left(best_v, length(p) - length(best_v))
            p = reduce_by_resultant(p, res)
            while n >= 1 && length(V[n]) > length(p)
               n -= 1
            end
         end
         # Insert p
         if !iszero(p)
            if n >= 1
               v = V[n]
               if length(p) == length(v)
                  if ((flag, q) = divides(leading_coefficient(v), leading_coefficient(p)))[1]
                     # V[n] can be swapped with p and reduced
                     V[n], r = reduce_tail(p, view(V, 1:n - 1), res), V[n]
                     r -= q*p
                     r = reduce_by_resultant(r, res)
                     if !iszero(r)
                        mypush!(H, r)
                     end
                  else # use gcdx
                     g, s, t = gcdx(leading_coefficient(v), leading_coefficient(p))
                     p, r = s*v + t*p, p # p has leading coefficient g dividing that of r (old p) and v
                     V[n] = reduce_tail(p, view(V, 1:n - 1), res)
                     q = divexact(leading_coefficient(r), g)
                     r -= q*p
                     r = reduce_by_resultant(r, res)
                     if !iszero(r)
                        mypush!(H, r)
                     end
                     q = divexact(leading_coefficient(v), g)
                     v -= q*p
                     v = reduce_by_resultant(v, res)
                     if !iszero(v)
                        mypush!(H, v)
                     end
                     p = reduce_by_resultant(p, res)
                  end
               else # length(p) > length(v)
                  if divides(leading_coefficient(v), leading_coefficient(p))[1]
                     # p can be inserted
                     insert!(V, n + 1, reduce_tail(p, view(V, 1:n), res))
                     n += 1
                  else # use gcdx
                     g, s, t = gcdx(leading_coefficient(v), leading_coefficient(p))
                     r = s*shift_left(v, length(p) - length(v)) + t*p # r has leading coeff g dividing that of p
                     q = divexact(leading_coefficient(p), g)
                     p -= r*q
                     p = reduce_by_resultant(p, res)
                     if !iszero(p)
                        mypush!(H, p)
                     end
                     p = r
                     # p can be inserted
                     insert!(V, n + 1, reduce_tail(p, view(V, 1:n), res))
                     n += 1
                  end
               end
            else # p is the smallest polynomial
               p = divexact(p, canonical_unit(p))
               insert!(V, 1, p)
               n = 1
            end
            # adjust entries above p
            orig_n = n
            while n < length(V)
               v = V[n + 1]
               if ((flag, q) = divides(leading_coefficient(v), leading_coefficient(p)))[1]
                  v -= q*shift_left(p, length(v) - length(p))
                  v = reduce_by_resultant(v, res)
                  if !iszero(v)
                     mypush!(H, v)
                  end
                  deleteat!(V, n + 1)
               elseif divides(leading_coefficient(p), leading_coefficient(v))[1]
                  break
               else
                  # use gcdx
                  g, s, t = gcdx(leading_coefficient(v), leading_coefficient(p))
                  p = s*v + t*shift_left(p, length(v) - length(p)) # p has leading coefficient g dividing that of v
                  q = divexact(leading_coefficient(v), g)
                  r, V[n + 1] = v - q*p, reduce_tail(p, view(V, 1:n), res)
                  r = reduce_by_resultant(r, res)
                  if !iszero(r)
                     mypush!(H, r)
                  end
                  if isterm(V[orig_n])
                     V[n + 1] = reduce_tail(V[n + 1], view(V, orig_n:orig_n), res)
                  end
                  n += 1
               end
            end
            for k = orig_n:2:min(n + 1, length(V))
               insert_spoly(V, H, k, res)
            end
            if isterm(V[orig_n])
               while n < length(V)
                  V[n + 1] = reduce_tail(V[n + 1], view(V, orig_n:orig_n), res)
                  n += 1
               end
            end
         end
      end
   end
   return V
end

# Given a nonempty vector V of polynomials of satisfying 1-6 below and a
# polynomial p whose degree is at least that of all the polynomials in V, add
# p to V and perform reduction steps so that 1-6 still hold. Fragments which
# are generated and need to be added go in H.
function extend_ideal_basis(D::Vector{T}, p::T, V::Vector{T}, H::Vector{T}, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   while true
      n = length(V)
      lc = leading_coefficient(V[n])
      # check if p can be added without any reduction
      if length(p) > length(V[n]) && !divides(leading_coefficient(p), lc)[1] && ((_, q) = divides(lc, leading_coefficient(p)))[1]
         p = reduce_tail(p, V, res)
         V = push!(V, p)
         insert_spoly(V, H, n + 1, res)
         return V
      end
      # check if p and V[n] are constant
      if isconstant(V[n]) && isconstant(p)
         return [parent(p)(gcd(constant_coefficient(p), constant_coefficient(V[n])))]
      end
      # check if p can replace V[n]
      swap = false
      if length(p) == length(V[n]) && !isunit(lc) && ((_, q) = divides(lc, leading_coefficient(p)))[1]
         p, V[n] = V[n], reduce_tail(p, view(V, 1:n - 1), res)
         insert_spoly(V, H, n, res)
         swap = true
      end
      # check if leading coefficients divide leading_coefficient of p
      while n >= 1 && (swap || ((_, q) = divides(leading_coefficient(p), leading_coefficient(V[n])))[1])
         best_v = V[n]
         best_size = mysize(V[n])
         best_q = q
         n1 = n - 1
         while n1 >= 1 && ((_, q) = divides(leading_coefficient(p), leading_coefficient(V[n1])))[1]
            new_size = mysize(V[n1])
            if new_size < best_size
               best_v = V[n1]
               best_size = new_size
               best_q = q
            end
            n1 -= 1
         end
         p -= best_q*shift_left(best_v, length(p) - length(best_v))
         p = reduce_by_resultant(p, res)
         while n >= 1 && length(V[n]) > length(p)
            n -= 1
         end
         swap = false
      end
      if iszero(p) # p was absorbed, yay!
         return V
      end
      if n < length(V) # we made some progress, p is a fragment
         push!(H, p)
         return V
      end
      # we made no progress, use gcdx
      n = length(V)
      v = V[n]
      g, s, t = gcdx(leading_coefficient(p), leading_coefficient(v))
      r = s*p + t*shift_left(v, length(p) - length(v)) # r has leading coeff g
      q = divexact(leading_coefficient(p), g)
      p -= q*r
      if length(r) == length(v) # V[n] can be reduced by r and switched
         q = divexact(leading_coefficient(v), g)
         r, V[n] = v - q*r, reduce_tail(r, view(V, 1:n - 1), res)
         insert_spoly(V, H, n, res)
         r = reduce_by_resultant(r, res)
         p = reduce_by_resultant(p, res)
         if length(r) > length(p)
            r, p = p, r
         end
         if length(p) == 0 # both polynomials were absorbed, yay!
            return V
         end
         if length(r) == 0 # one polynomial was absorbed, yay!
            push!(H, p)
            return V
         else
            push!(H, p, r)
            return V
         end
      else # length(r) > length(v)
         r = reduce_by_resultant(r, res)
         p = reduce_by_resultant(p, res)
         if !iszero(p)
            if length(p) >= length(v)
              push!(D, r)
            else 
               push!(H, p)
               p = r
            end
         else
            p = r 
         end
         # p still needs reduction
      end
   end
end
      
# We call an ideal over a polynomial ring over a Euclidean domain reduced if
# 1. There is only one polynomial of each degree in the ideal
# 2. The degree of polynomials in the basis increases
# 3. The leading coefficient of f_i divides that of f_{i-1} for all i
# 4. Only the final polynomial may have leading coefficient that is a unit
# 5. The polynomials are all canonicalised (divided by their canonical_unit)
# 6. The tail of each polynomial is reduced mod the other polynomials in the basis
# 7. All the S-polynomials of pairs in the basis reduce to zero
function reduce(I::Ideal{T}; complete_reduction::Bool=true) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   if hasmethod(gcdx, Tuple{U, U})
      V = gens(I)
      # Compute V a vector of polynomials giving the same basis as I
      # but for which the above hold
      if length(V) > 1
         # compute resultant
         S = parent(V[1])
         D = sort(V, by=mysize2, rev=true)
         res = zero(base_ring(D[1]))
         while !isempty(D) && isconstant(D[end])
            di = constant_coefficient(pop!(D))
            res = gcd(di, res)
         end
         g = one(base_ring(S))
         for i = 1:length(D)
            for j = i + 1:length(D)
               res = gcd(resultant(D[i], D[j]), res)
               if isunit(res)
                  break
               end
            end
         end
         if iszero(res) # remove common gcd of all polys
            g = zero(S)
            for v in D
               g = gcd(v, g)
            end
            for i in 1:length(D)
               D[i] = divexact(D[i], g)
            end
            for i = 1:length(D)
               for j = i + 1:length(D)
                  res = gcd(resultant(D[i], D[j]), res)
                  if isunit(res)
                     break
                  end
               end
            end
         end
         if isunit(res) # everything reduces to 0
            V = [one(S)*g]
         else
            res = divexact(res, canonical_unit(res))
            V = [S(res)]
            B0 = [reduce_by_resultant(p, res) for p in D]
            B0 = filter(!iszero, B0)
            B = sort(B0, by=mysize2, rev=true)
            H = T[]
            V = extend_ideal_basis(B, V, H, res)
            # multiply polys by original gcd
            if !isone(g)
               for i = 1:length(V)
                  V[i] *= g
               end
            end
         end
      elseif length(V) == 1
         d = V[1]
         V = [divexact(d, canonical_unit(d))]
      end
      if complete_reduction
         for i = 2:length(V)
            V[i] = reduce_tail(V[i], V[1:i - 1], zero(base_ring(V[1])))
         end
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

