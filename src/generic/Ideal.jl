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

@doc raw"""
    gens(I::Ideal{T}) where T <: RingElement

Return a list of generators of the ideal `I` in reduced form and canonicalised.
"""
gens(I::Ideal{T}) where T <: RingElement = I.gens

###############################################################################
#
#   Heap and nodes
#
###############################################################################

const IDEAL_MULTIV_DEBUG = false # whether to print debugging information
const poly_level = 1 # 1 = display leading terms only, 2 = display entire polynomials
const print_inactive = false # whether to print polynomials and node nums of inactive nodes
const print_prompt = false # whether to wait for a keypress between debug output

print_node_level = [0] # 0 = print only polynomials in nodes, 1 = print lattice of nodes

node_num = [0] # for enumerating nodes when printing

# The main node structure used in the lattice of leading monomials
# Each node contains a polynomial and a lattice of nodes shows which
# leading monomials are divisible by which others
# The nodes directly above a node (i.e. ones it divides) are given
# by the field `up` and then by `next.up`, `next.next.up`, etc.
# Nodes with polynomials having equal leading monomials are chained
# together via the `equal` field in a closed loop that eventually
# links back to the present node
# Only the nodes on the backbone of the lattice show divisibility
# via the `up` and `next` fields; `equal` nodes have the same
# divisibility relations with other nodes and so the information
# is not duplicated for these
# When reduction occurs, the best reducer node is first attached
# to the node at the `reducer` field
# Nodes are marked as inactive by setting the `active` field to
# `false` if they are no removed from the basis/lattice
# The inactive nodes are retained as potential reducers in the
# lattice, and possibly purged at regular intervals
# One does not have to check if a node can be reduced by a node
# that was already in the lattice before the last lot of reductions
# if the node could not be reduced by anything at the time of the
# last round of reductions, but if a node is new since the last
# round then it is marked as such by setting `new_node` to true
# since it could potentially be a reducer or reduced by anything
# There are various kinds of reduction possible, including having
# the leading coefficient removed and having the leading coefficient
# reduced and nodes that had no reducers of a given kind at a given
# round are marked as such via the `settled` field so that they are
# not checked for this against old nodes again
# The leading monomial of a node is stored for easy access in the
# `lm` field and the least common multiple of all leading monomials
# in the subtree rooted at a given node are also stored in that node
# The `path` flag is used to mark visited nodes when traversing the
# lattice
# Each node is given a unique number so they can be printed easily
# The size of a node with respect the heuristic used to sort
# potential reducers by `size` is cached on the node in the `size`
# field

mutable struct lmnode{U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N}
   poly::Union{U, Nothing}
   up::Union{lmnode{U, V, N}, Nothing}   # out (divisible leading monomials)
   next::Union{lmnode{U, V, N}, Nothing} # extend current node so out-degree can be > 1
   equal::Union{lmnode{U, V, N}, Nothing} # to chain nodes with equal lm
   reducer::Union{lmnode{U, V, N}, Nothing} # best reducer found so far for node
   active::Bool # whether polynomial is still actively being reduced
   new_node::Bool # whether the node was added to the lattice since last time
   settled::Int # 0 = newly reduced, 1 = no reducers remove lc, 2 = none reduce lc, 3 = no reducers
   new_settled::Int # used for temporary when computing settled
   lm::NTuple{N, Int}   # leading monomial as exponent vector
   lcm::NTuple{N, Int}  # lcm of lm's in tree rooted here, as exponent vector
   path::Bool # used for marking paths to divisible nodes
   num::Int # for counting nodes (used in printing)
   size::Float64 # caches the size of the polynomial

   function lmnode{U, V, N}(p::Union{U, Nothing}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N}
      node = new{U, V, N}(p, nothing, nothing, nothing, nothing, true, true, 0, 0)
      node.size = 0.0
      if p != nothing && !iszero(p)
         if leading_coefficient(p) < 0
            node.poly = -node.poly
         end
         node.lm = Tuple(exponent_vector(node.poly, 1))
         node.lcm = node.lm
         node.size = reducer_size(node)
      end
      node.equal = node
      node.path = false
      node_num[] += 1
      node.num = node_num[]
      return node::lmnode{U, V, N}
   end
end

# compute the lcm of the leading monomials of the given nodes as a tuple
function lm_lcm(node1::lmnode{U, V, N}, node2::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N}
   return Tuple(max(node1.lcm[i], node2.lcm[i]) for i in 1:N)
end

# return `true` if the leading monomial of `node2` divides the lcm of the
# leading monomials of the tree rooted at `node1`, i.e. return `true` if
# `node2` could divide some node in the tree rooted at `node1`
function lm_divides_lcm(node1::lmnode{U, V, N}, node2::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N}
   for i = 1:N
      if node2.lm[i] > node1.lcm[i]
         return false
      end
   end
   return true
end

# print "nothing" at up, next nodes which lead nowhere
function show_inner(io::IO, n::Nothing)
   print(io, "nothing")
end

# print a polynomial, either just the leading monomial if the
# `poly_level` debug flag above is 1 and the whole thing
# otherwise
function print_poly(io::IO, p::MPolyRingElem)
   if poly_level == 1
      print(io, leading_term(p))
      if length(p) > 1
         print(io, "+...")
      end
   else
      print(io, p)
   end
end

# print a node and all its fields recursively, excluding
# inactive nodes if the debug flag `print_inactive` above
# is set to `false`
function show_inner(io::IO, n::lmnode)
   if !n.path
      print(io, "Node(")
      if print_inactive || n.active
         print(io, "p", n.num)
         if print_inactive
            if n.active
               print(io, "(Y)")
            else
               print(io, "(N)")
            end
         end
         print(io, "=")
         print_poly(io, n.poly)
      else
         print(io, "_")
      end
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

# print a vector of nodes
function show(io::IO, B::Vector{lmnode{U, V, N}}) where {U, V, N}
   if print_node_level[] == 0
      BB = extract_gens(B)
      print(io, "[")
      for i = 1:length(BB) - 1
         print_poly(io, BB[i])
         print(io, ", ")
      end
      if !isempty(BB)
         print_poly(io, BB[end])
      end
      print(io, "]")
   else
      print(io, "[")
      for i in 1:length(B) - 1
         show_inner(io, B[i])
         print(io, ", ")
      end
      if !isempty(B)
         show_inner(io, B[end])
      end
      print(io, "]")
      clear_path(B)
   end
end

# heap implementation for sorting polys by lm, head = smallest

# Defined in MPoly, repeated here for reference
# heapleft(i::Int) = 2i
# heapright(i::Int) = 2i + 1
# heapparent(i::Int) = div(i, 2)

function lm_precedes(node1::lmnode{U, :lex, N}, node2::lmnode{U, :lex, N}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, N}
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

function lm_precedes(node1::lmnode{U, :deglex, N}, node2::lmnode{U, :deglex, N}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, N}
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

function lm_precedes(node1::lmnode{U, :degrevlex, N}, node2::lmnode{U, :degrevlex, N}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, N}
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

function lm_divides(node1::lmnode{U, V, N}, node2::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N}
   for i = 1:N
      if node2.lm[i] > node1.lm[i]
         return false
      end
   end
   return true
end

function heapinsert!(heap::Vector{lmnode{U, V, N}}, node::lmnode{U, V, N}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N}
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
   return nothing
end
   
function heappop!(heap::Vector{lmnode{U, V, N}}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N}
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
   x.up = nothing
   x.next = nothing
   x.lm = Tuple(exponent_vector(x.poly, 1))
   x.lcm = x.lm
   return x
end

###############################################################################
#
#   Ideal reduction for multivariates over Euclidean domain
#
###############################################################################

# Reduce coefficients of polynomial p starting from term with index `start`
# (numbered from 1). Used for `normal_form` and `tail_reduce`.
function reduce_coefficients(p::U, V::Vector{U}, start::Int) where U <: AbstractAlgebra.MPolyRingElem
   n = start
   len = length(p)
   infl = [1 for i in 1:nvars(parent(p))]
   while n <= len
      for i = 1:length(V)
         c = coeff(p, n)
         mv = exponent_vector(V[i], 1)
         mp = exponent_vector(p, n)
         if max.(mv, mp) == mp # leading monomial divides
            h = leading_coefficient(V[i]) # should be positive
            q, r = AbstractAlgebra.divrem(c, h)
            if !iszero(q)
               shift = mp .- mv
               u = inflate(V[i], shift, infl)
               p -= q*u
               if iszero(r)
                  n -= 1
                  break
               end
            end        
         end
      end
      n += 1
      len = length(p)
   end
   return p
end

# return the polynomial `p` tail reduced by the given basis `V`
function tail_reduce(p::U, V::Vector{U}) where U <: AbstractAlgebra.MPolyRingElem
   len = length(p)
   if len <= 1
      return p
   end
   return reduce_coefficients(p, V, 2)
end

# return the normal form of the polynomial `p` after reduction by
# the given basis `V`
# this is internal, the user facing version being given below
function normal_form(p::U, V::Vector{U}) where U <: AbstractAlgebra.MPolyRingElem
   if iszero(p)
      return zero(parent(p))
   end
   return reduce_coefficients(p, V, 1)
end

# given two nodes, return a node giving the spoly of the two nodes
function compute_spoly(f::T, g::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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

# given two nodes, return a node giving the gpoly of the two nodes
function compute_gpoly(f::T, g::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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

# given two polynomials, return the spoly of the two
function spoly(f::T, g::T) where T <: MPolyRingElem
   n = nvars(parent(f))
   fc = leading_coefficient(f)
   gc = leading_coefficient(g)
   mf = exponent_vector(f, 1)
   mg = exponent_vector(g, 1)
   c = lcm(fc, gc)
   llcm = max.(mf, mg)
   infl = [1 for i in 1:n]
   shiftf = llcm .- mf
   shiftg = llcm .- mg
   s = divexact(c, fc)*inflate(f, shiftf, infl) - divexact(c, gc)*inflate(g, shiftg, infl)
end

# given two polynomials, return the gpoly of the two
function gpoly(f::T, g::T) where T <: MPolyRingElem
   n = nvars(parent(f))
   fc = leading_coefficient(f)
   gc = leading_coefficient(g)
   mf = exponent_vector(f, 1)
   mg = exponent_vector(g, 1)
   _, s, t = gcdx(fc, gc)
   llcm = max.(mf, mg)
   infl = [1 for i in 1:n]
   shiftf = llcm .- mf
   shiftg = llcm .- mg
   g = s*inflate(f, shiftf, infl) + t*inflate(g, shiftg, infl)
end

# used for sorting polynomials in final basis, first by leading monomial and
# then by leading coefficient
function isless_monomial_lc(p::U, q::U) where U <: AbstractAlgebra.MPolyRingElem{<:RingElement}
   plm = leading_monomial(p)
   qlm = leading_monomial(q)
   if plm < qlm
      return true
   end
   if plm == qlm && leading_coefficient(p) < leading_coefficient(q)
      return true
   end
   return false
end

# heuristic for size of reducer polynomials (smaller is better), used to sort
# potential reducers
function reducer_size(f::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   if f.size != 0.0
      return f.size
   end
   s = 0.0
   # heuristic really punishes polys with lots of terms
   len = length(f.poly)
   j = len
   for i = 1:len
      c = coeff(f.poly, i)
      # larger coefficients are also somewhat punished
      s += Base.log(ndigits(c; base=2))*j*j
      j -= 1      
   end
   return s
end

# returns the leading coefficient of a node for use in sorting
function leading_coefficient_size(f::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   return leading_coefficient(f.poly)
end

# reduce the node `b` by the reducer attached to it, inactivating
# `b` if its leading term is removed, and putting any fragments in
# `H` and marking `b` for recomputation of spolys (by placing it in
# `S`) if it has been merely modified
function reduce_by_reducer(S::Vector{T}, H::Vector{T}, b::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   reduced = false
   if b.reducer != nothing
      c = leading_coefficient(b.poly)
      h = leading_coefficient(b.reducer.poly)
      q, r = AbstractAlgebra.divrem(c, h)
      coeff_divides = false
      if iszero(r)
         coeff_divides = true
      end
      if iszero(q) # need gcd polynomial
         sp = compute_spoly(b, b.reducer)
         g, s, t = gcdx(c, h)
         infl = [1 for in in 1:N]
         shift = exponent_vector(b.poly, 1) .- exponent_vector(b.reducer.poly, 1)
         d = s*b.poly + t*inflate(b.reducer.poly, shift, infl)
         q = divexact(c, g)
         p = b.poly - q*d # reduce b by d (the gcd poly)
         if leading_coefficient(d) < 0
            d = -d
         end
         b.poly = d # replace b with gcd poly
         b.size = 0.0 # update size
         b.size = reducer_size(b)
         b.new_node = true
         # must recompute s-polys
         push!(S, b)
         if !iszero(p)
            push!(H, lmnode{U, V, N}(p))
         end
      else
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
            if leading_coefficient(d) < 0
               d = -d
            end
            b.poly = d
            b.size = 0.0 # update size
            b.size = reducer_size(b)
            b.new_node = true
            # must recompute s-polys
            push!(S, b)
         end
      end
      reduced = true
   end
   b.reducer = nothing
   return reduced
end

# reduce nodes in a tree depth first
# X is used as temporary space to sort equal nodes
# H is used for fragments that must be inserted into the basis
# S is used to mark nodes for computation of spolys when they cease
# being reduced
# returns true if some reduction actually occurred
function reduce_nodes(S::Vector{T}, H::Vector{T}, b::T, X::Vector{T}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   b.path = true
   b.settled = b.new_settled
   b.new_settled = 0
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
      b2.settled = b2.new_settled
      b2.new_settled = 0
      index = searchsortedfirst(X, b2, by=leading_coefficient_size, rev=true) # find index at which to insert x
      insert!(X, index, b2)
   end
   # of the active nodes, remove ones which are the same up to units
   filter!(x->x.active, X)
   i = 1
   len = length(X)
   while i < len
      c1 = leading_coefficient(X[i].poly)
      j = i + 1
      while j <= len && leading_coefficient(X[j].poly) == c1
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
   # clear nodes from temporary space X
   while !isempty(X)
      pop!(X)
   end
   return reduced
end

# actually performs the reductions which have been attached to nodes
# by best_reducer
# X is used to sort equal nodes
# H is used for fragments
# S is used to mark nodes for generation of spolys when they stop
# being reduced
# returns true if some reduction actually occurred
function reduce_nodes(S::Vector{T}, H::Vector{T}, B::Vector{T}, X::Vector{T}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   reduced = false
   for b in B
      if !b.path
         reduced |= reduce_nodes(S, H, b, X)
      end
   end
   return reduced
end

# given a node b and a list of potential reducers, find the
# best one which will remove the leading coefficient of `b`
# if it exists
# the current best reducer and whether it would remove
# the leading coefficient are also passed in
function find_best_divides(b::T, X::Vector{T}, best::Union{T, Nothing}, best_divides::Bool) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   c = leading_coefficient(b.poly)
   best_size = best == nothing ? 0.0 : reducer_size(best)
   for i = length(X):-1:1
      # poly can't be reduced by itself or by a node that will be reduced
      # in this round
      if X[i] != b && X[i].reducer == nothing
         h = leading_coefficient(X[i].poly)
         if divides(c, h)[1]
            usable = true
            # make sure a circular chain of reduction can't occur
            if X[i].lm == b.lm
               x2 = X[i]
               while x2 != nothing
                  if x2 == b
                     usable = false
                  end
                  x2 = x2.reducer
               end
            end
            # if there's currently no reducer or this one is better
            # attach the reducer that has been found
            # only reduce by nodes with same lm if active, else it
            # could be removed by an equal node or fragments thereof
            if usable && (c != h || X[i].active)
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

# given a node b and a list of potential reducers, find the
# best one which will reduce the leading coefficient of `b`
# if it exists
# the current best reducer and whether it would reduce
# the leading coefficient are also passed in
function find_best_reduces(b::T, X::Vector{T}, best::Union{T, Nothing}, best_reduces::Bool) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   c = leading_coefficient(b.poly)
   best_size = best == nothing ? 0.0 : reducer_size(best)
   for i = length(X):-1:1
      # poly can't be reduced by itself or a node that will be reduced
      # in this round
      if X[i] != b && X[i].active && X[i].reducer == nothing
         h = AbstractAlgebra.mod(c, leading_coefficient(X[i].poly))
         if h != c && h != 0 && !divides(leading_coefficient(X[i].poly), c)[1]
            usable = true
            # check there are no cycles of reducers
            if X[i].lm == b.lm
               x2 = X[i]
               while x2 != nothing
                  if x2 == b
                     usable = false
                  end
                  x2 = x2.reducer
               end
            end
            # if node currently doesn't have a reducer or this one is
            # better, attach the reducer just found
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

# given a node b and a list of potential reducers, find the
# best one which will give a smaller leading coefficient if
# we gcd with the leading coefficient of `b` if such exists
# the current best reducer and whether it would gcd to reduce
# the leading coefficient are also passed in
function find_best_gcd(b::T, X::Vector{T}, best::Union{T, Nothing}, best_is_gcd::Bool) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   c = leading_coefficient(b.poly)
   best_size = best == nothing ? 0.0 : reducer_size(best)
   for i = length(X):-1:1
       # poly can't be reduced by itself or by node that will be reduced
       # in this round
      if X[i] != b && X[i].active  && X[i].reducer == nothing
         if !divides(leading_coefficient(X[i].poly), c)[1] &&
            !divides(c, leading_coefficient(X[i].poly))[1]
            usable = true
            # check there are no cycles of reducers
            if X[i].lm == b.lm
               x2 = X[i]
               while x2 != nothing
                  if x2 == b
                     usable = false
                  end
                  x2 = x2.reducer
               end
            end
            # if node currently doesn't have a reducer or this one is
            # better, attach the reducer just found
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
# everything else and reduction of leading coefficient after that
# the function takes a list of potential reducers (i.e. nodes from further down
# the tree), with the ones that are new since last round being in `Xnew`
function find_best_reducer(b::T, X::Vector{T}, Xnew::Vector{T}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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
      b.new_settled = 0
      return best
   else
      b.new_settled = 1
   end
   # 2. check for leading coefficient that reduces c
   best_reduces = best == nothing ? false : AbstractAlgebra.mod(c, leading_coefficient(best.poly)) != c
   if b.settled <= 1
      best, best_reduces = find_best_reduces(b, X, best, best_reduces)
   end
   best, best_reduces = find_best_reduces(b, Xnew, best, best_reduces)
   if best_reduces
      return best
   else
      b.new_settled = 2
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
      b.new_settled = 3 # no reducer found
   end
   return best
end

# attach best reducer to each node in tree
# `X` and `Xnew` will be used to make a list of potential
# reducers (i.e. nodes further down in tree), with `Xnew`
# being those that are new since last round
function best_reducer(b::T, X::Vector{T}, Xnew::Vector{T}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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

# Setup for one round of reduction
# Each node will be reduced mod the leading coeff of the best node
# whose leading term divides it and if there are none, by the
# best node whose leading monomial divides it and if none exists
# gcd polynomials will be created if possible
# The reductions are not actually done here; the best reducer
# is just attached to the node
# The potential reducers along the current path so far are given by
# the array `X` unless they are for new nodes since last time, in
# which case they will go in `Xnew`
function best_reducer(B::Vector{T}, X::Vector{T}, Xnew::Vector{T}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   for b in B
      best_reducer(b, X, Xnew)
   end
end

# insert fragments from `H` into basis `B`
# they are sorted by `size` of leading monomial and only
# those below a current bound are inserted
# nodes are marked for generation of spolys when they stop
# being reduced, by placing them in `S`
function insert_fragments(S::Vector{T}, B::Vector{T}, H::Vector{T}, bound::NTuple{N, Int}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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
function basis_insert(b::T, d::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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

# insert a node `d` into the basis `B`
function basis_insert(S::Vector{T}, B::Vector{T}, d::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   # first link `d` to nodes it divides
   insert_links(d, B)
   # now insert it if any existing nodes divide it
   inserted = false
   for b in B
      if lm_divides(d, b)
         inserted = true
         basis_insert(b, d)
      end
   end
   # otherwise just add it to the basis as an independent root
   if inserted == false # was not inserted
      push!(B, d)
   end
   clear_path(B)
   # if the polynomial wasn't a constant multiple of an existing one
   # (which is flagged using the `active` field), mark it for later
   # computation of spolys when it stops being reduced
   if d.active
      push!(S, d) # store for later computation of spolys
      d.settled = 0
   end
end

# give a node in the lattice, find the node with equal lm
# which is on the backbone of the lattice
# if all nodes with equal lm seem to be leaves, just declare
# `d` itself to be the base node (which won't matter where)
# this function is used
function find_base(d::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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

# return true if `d`` has a connection up to `b`, i.e.
# `d` divides `b`
function connected(b::T, d::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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
function compute_spolys(S::Vector{T}, b::T, d::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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
      if n.active && n.settled == 3
         if n != d && !lm_divides(n, d) && !lm_divides(d, n)
            s = compute_spoly(n, d)
            push!(S, s)
            g = compute_gpoly(n, d)
            push!(S, g)
            break
         end
         dbase = find_base(d)
         if connected(b, dbase) || connected(dbase, b)
            s = compute_spoly(n, d)
            push!(S, s)
            break
         end 
      end
      n = n.equal
      nend = b
   end
   b.path = true
   return nothing
end

# compute spolys of all polys in basis `B` with `d`, except `d itself
# and polys whose lc divides that of `d or which that of `d divides
# unless directly connected to `d`
# s-polys are collected in `S` (not in fragments `H` to avoid too many
# s-polys being added in a single round)
# they are later moved to fragments
function compute_spolys(S::Vector{T}, B::Vector{T}, d::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   for b in B
      compute_spolys(S, b, d)
   end
   clear_path(B)
end

# insert link from d to nodes it divides
function insert_links(d::T, b::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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
function insert_links(d::T, B::Vector{T}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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
function compute_lcm(b::T)  where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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

# clear all path flags in given branch, i.e. reset path flags
function clear_path(b::T) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
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

# clear all path flags in basis, i.e. reset path flags
function clear_path(B::Vector{T}) where {U <: AbstractAlgebra.MPolyRingElem{<:RingElement}, V, N, T <: lmnode{U, V, N}}
   for b in B
      clear_path(b)
   end
end

# extract generators (polynomials) from branch and store in
# array `D`
# only polys from active nodes are extracted
function extract_gens(D::Vector{U}, node::T) where {N, U <: AbstractAlgebra.MPolyRingElem, V, T <: lmnode{U, V, N}}
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
function extract_gens(B::Vector{T}) where {N, U <: AbstractAlgebra.MPolyRingElem, V, T <: lmnode{U, V, N}}
   D = Vector{U}()
   for node in B
      extract_gens(D, node)
   end
   clear_path(B)
   return D
end

# for all nodes in `S2` (marked for spoly generation) which are
# still active, generate spolys against all relevant polys in
# basis `B` and place the spolys in `S` (for later insertion
# into fragments)
function generate_spolys(S::Vector{T}, B::Vector{T}, S2::Vector{T}) where {N, U <: AbstractAlgebra.MPolyRingElem, V, T <: lmnode{U, V, N}}
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
function reduce_gens(I::Ideal{U}; complete_reduction::Bool=true) where {T <: RingElement, U <: AbstractAlgebra.MPolyRingElem{T}}
   node_num[] = 0
   B = gens(I)
   # nothing to be done if only one poly
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
      if is_constant(d.poly)
         S = parent(d.poly)
         d0 = constant_coefficient(d.poly)
         while !isempty(heap) && is_constant(heap[1].poly)
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
         push!(H, d)
      end
      # do reduction
      G = Vector{U}()
      X = Vector{lmnode{U, V, N}}()
      X2 = Vector{lmnode{U, V, N}}()
      X2new = Vector{lmnode{U, V, N}}()
      while true
         reduction_occurs = true
         while reduction_occurs
            # attach best reducers to nodes
            best_reducer(B2, X2, X2new)
            # do reductions
            reduction_occurs = reduce_nodes(S2, H, B2, X)
            clear_path(B2)
            # move s-polys from S to fragments
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
                     p = normal_form(d.poly, G)
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
if IDEAL_MULTIV_DEBUG
            print_node_level[] = 1
            println("B2 = ", B2)
            print_node_level[] = 0
            println("S = ", S)
            println("S2 = ", S2)
            println("H = ", H)
            println("")
            if print_prompt
               readline(stdin)
            end
end
            # when no further reduction is possible in lattice
            # generate spolys
            if !reduction_occurs
               generate_spolys(S, B2, S2)
            end
         end
         # if we have no fragments to insert
         if isempty(H)
            # and no spolys to compute or insert
            if isempty(S) && isempty(S2)
               # we are done
               break
            end
         else
            # else if there are fragments, increase
            # bound so next one at least will be inserted
            bound = max.(bound, H[end].lm)
         end
      end
      # get all still active nodes from basis, we are done
      B = extract_gens(B2)
   end
   # canonicalise basis elements
   B = [divexact(d, canonical_unit(d)) for d in B]
   # do tail reduction if requested
   if complete_reduction
      B = [tail_reduce(d, B) for d in B]
   end
   # sort by leading monomial then leading coefficient
   B = sort!(B, lt = isless_monomial_lc)
   return Ideal{U}(base_ring(I), B)
end

###############################################################################
#
#   Ideal reduction for univariates over Euclidean domain
#
###############################################################################

const IDEAL_UNIV_DEBUG = false

# Debugging: check that polynomials are in order of degree
function check_d(D::Vector{T}, m::Int) where {T <: AbstractAlgebra.PolyRingElem{<:RingElement}}
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

# debugging: check polynomials are in order of degree, that there is not more than one
# polynomial with a unit as leading coefficient, that coefficients divide one another
# and that there are not two leading coefficients that are associate
function check_v(V::Vector{T}, m::Int) where {T <: AbstractAlgebra.PolyRingElem{<:RingElement}}
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
      if is_unit(leading_coefficient(V[n])) && is_unit(leading_coefficient(V[n - 1]))
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

# heuristic "size" of polynomials
function mysize(f::T) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
   s = 0.0
   j = 0
   # heuristic really punishes polys with lots of terms
   for i = 1:length(f)
      c = coeff(f, i - 1)
      if !iszero(c)
         j += 1
         # large coeffs are also punished
         s += Base.log(ndigits(coeff(f, i - 1); base=2))*j*j
      end
   end
   return s
end

# heuristic as above, but first orders polynomials with respect to
# degree
function mysize2(f::T) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
   return (length(f) << 40) + mysize(f)
end

# insert polynomial into H in sorted order wrt mysize
function mypush!(H::Vector{T}, f::T) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
   index = searchsortedfirst(H, f, by=mysize, rev=true) # find index at which to insert x
   insert!(H, index, f)
end

# extend the basis `V` of polynomials satisfying 1-6 below by the
# polynomials in `D`, all of which have degree at least that of those in
# `V` and such that the degrees of the polynomials in `D` is *decreasing*.
function extend_ideal_basis(D::Vector{T}, V::Vector{T}, H::Vector{T}, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
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

# reduce a polynomials coefficients by a resultant `res` that
# is in the ideal
function reduce_by_resultant(f::T, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
   f = deepcopy(f)
   len = 0
   for i = 1:length(f)
      c = coeff(f, i - 1)
      q, r = AbstractAlgebra.divrem(c, res)
      setcoeff!(f, i - 1, r)
      if !iszero(r)
         len = i
      end
   end
   if len < length(f)
      f = set_length!(f, len)
   end
   return divexact(f, canonical_unit(f))
end

# do tail reduction of the polynomial `f` by the basis `V`
function reduce_tail(f::T, V::AbstractVector{T}, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
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
            h = leading_coefficient(V[i]) # should be non-negative
            q, r = AbstractAlgebra.divrem(c, h)
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

# return the normal form of `p` by the basis `V`
# this is internal; see below for the user facing version
function normal_form(p::T, V::Vector{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
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
            h = leading_coefficient(V[i]) # should be non-negative
            q, r = AbstractAlgebra.divrem(c, h)
            u = shift_left(V[i], n - length(V[i]))
            p -= q*u       
         end
         n = min(n - 1, length(p))
      end
   end
   return p
end

# an S-polynomial of f and g is a linear combination which makes the leading
# coefficients cancel, e.g. if p = 3x^2, q = 2x then S-poly(p, q) = 2*p - 3x*q
# this function inserts spolys of `v[n]` and `V[n - 1]` and of `V[n]` and
# and `V[n + 1]` (if those exist) into `H`, reducing coefficients by the
# resultant `res` if available
function insert_spoly(V::Vector{T}, H::Vector{T}, n::Int, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
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

# insert fragments from `H` into basis `V`, placing any long
# fragments back in `D` instead
function insert_fragments(D::Vector{T}, V::Vector{T}, H::Vector{T}, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
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
                  if is_term(V[orig_n])
                     V[n + 1] = reduce_tail(V[n + 1], view(V, orig_n:orig_n), res)
                  end
                  n += 1
               end
            end
            for k = orig_n:2:min(n + 1, length(V))
               insert_spoly(V, H, k, res)
            end
            if is_term(V[orig_n])
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

# given a nonempty vector `V` of polynomials of satisfying 1-6 below and a
# polynomial `p` whose degree is at least that of all the polynomials in `V`,
# add `p` to `V` and perform reduction steps so that 1-6 still hold
# fragments which are generated and need to be added go in `H`
function extend_ideal_basis(D::Vector{T}, p::T, V::Vector{T}, H::Vector{T}, res::U) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
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
      if is_constant(V[n]) && is_constant(p)
         return [parent(p)(gcd(constant_coefficient(p), constant_coefficient(V[n])))]
      end
      # check if p can replace V[n]
      swap = false
      if length(p) == length(V[n]) && !is_unit(lc) && ((_, q) = divides(lc, leading_coefficient(p)))[1]
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
      
function remove_constant_multiples(V::Vector{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
   len = length(V)
   for i = 1:len
      j = i + 1
      while j <= len
         if is_constant_multiple(V[j], V[i])
            deleteat!(V, j)
            len -= 1
            j -= 1
         end
         j += 1
      end
   end
   V = collect(reverse(V))
   for i = 1:len
      j = i + 1
      while j <= len
         if is_constant_multiple(V[j], V[i])
            deleteat!(V, j)
            len -= 1
            j -= 1
         end
         j += 1
      end
   end
   return V
end

# if possible, find a resultant (or preferably gcd of resultants)
# which exists in the idea
# this can be used to reduce coefficients of polynomials to keep
# them small throughout the algorithm
function find_resultant_in_ideal(R::Ring, D::Vector{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
   res = zero(R)
   V = similar(D, 0)
   for i = 1:length(D)
      for j = i + 1:length(D)
         r = resultant(D[i], D[j])
         if iszero(r)
            g = gcd(D[i], D[j])
            push!(V, g*resultant(divexact(D[i], g), divexact(D[j], g)))
         else
            res = gcd(r, res)
            if is_unit(res)
               break
            end
         end
      end
      if is_unit(res)
         break
      end
   end
   if !isempty(V)
      V = remove_constant_multiples(V)
      res = gcd(res, find_resultant_in_ideal(R, V))
   end
   return res
end

function is_constant_multiple(f::T, g::T) where T <: AbstractAlgebra.PolyRingElem
   if length(f) == length(g)
      c1 = leading_coefficient(f)
      c2 = leading_coefficient(g)
      if ((flag, q) = divides(c1, c2))[1]
         if f == g*q
            return true
         end
      end
   end
   return false
end

# We call an ideal over a polynomial ring over a Euclidean domain reduced if
# 1. There is only one polynomial of each degree in the ideal
# 2. The degree of polynomials in the basis increases
# 3. The leading coefficient of f_i divides that of f_{i-1} for all i
# 4. Only the final polynomial may have leading coefficient that is a unit
# 5. The polynomials are all canonicalised (divided by their canonical_unit)
# 6. The tail of each polynomial is reduced mod the other polynomials in the basis
# 7. All the S-polynomials of pairs in the basis reduce to zero
# The algorithm proceeds by starting with all polynomials in D which is sorted
# first by degree, then by "size"
# if possible, we find a resultant in the ideal which can be used to reduce
# coefficients
# polynomials are added one at a time to the basis `V`, reducing existing
# elements of the basis
# this may generate fragments which are placed in `H`
# when a new poly is added to the basis, spolys are also generated for all
# the polynomials reduced by the new poly and for the new poly itself
function reduce_gens(I::Ideal{T}; complete_reduction::Bool=true) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
   V = gens(I)
   # Compute V a vector of polynomials giving the same basis as I
   # but for which the above hold
   V = filter(!iszero, V)
   V = remove_constant_multiples(V)
   if length(V) > 1
      # compute resultant
      S = parent(V[1])
      R = base_ring(S)
      D = sort(V, by=mysize2, rev=true)
      res = zero(R)
      while !isempty(D) && is_constant(D[end])
         di = constant_coefficient(pop!(D))
         res = gcd(di, res)
      end
      g = one(R)
      for i = 1:length(D)
         for j = i + 1:length(D)
            res = gcd(resultant(D[i], D[j]), res)
            if is_unit(res)
               break
            end
         end
      end
      if iszero(res) # remove common gcd of all polys
         g = zero(S)
         for v in D
            g = gcd(v, g)
         end
         D2 = Vector{T}()
         r2 = zero(R)
         for i in 1:length(D)
            D[i] = divexact(D[i], g)
            if !is_constant(D[i])
               push!(D2, D[i])
            else
               r2 = gcd(r2, constant_coefficient(D[i]))
            end
         end
         res = find_resultant_in_ideal(R, D2)
         res = gcd(res, r2)
      end
      if is_unit(res) # everything reduces to 0
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
end

###############################################################################
#
#   Normal form
#
###############################################################################

@doc raw"""
    normal_form(p::U, I::Ideal{U}) where {T <: RingElement, U <: Union{AbstractAlgebra.PolyRingElem{T}, AbstractAlgebra.MPolyRingElem{T}}}

Return the normal form of the polynomial `p` with respect to the ideal `I`.
"""
function normal_form(p::U, I::Ideal{U}) where {T <: RingElement, U <: Union{AbstractAlgebra.PolyRingElem{T}, AbstractAlgebra.MPolyRingElem{T}}}
   return normal_form(p, gens(I))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(I::Ideal{T}, J::Ideal{T}) where T <: RingElement
   return gens(I) == gens(J)
end

###############################################################################
#
#   Containment
#
###############################################################################

@doc raw"""
    Base.contains(I::Ideal{T}, J::Ideal{T}) where T <: RingElement

Return `true` if the ideal `J` is contained in the ideal `I`.
"""
function Base.contains(I::Ideal{T}, J::Ideal{T}) where T <: RingElement
   G1 = gens(J)
   G2 = gens(I)
   if isempty(G1)
      return true
   end
   if isempty(G2)
      return false
   end
   return divides(G1[1], G2[1])[1]
end

function Base.contains(I::Ideal{T}, J::Ideal{T}) where {U <: RingElement, T <: Union{AbstractAlgebra.PolyRingElem{U}, AbstractAlgebra.MPolyRingElem{U}}}
   G = gens(J)
   for v in G
      if !iszero(normal_form(v, I))
         return false
      end
   end
   return true
end

###############################################################################
#
#   Intersection
#
###############################################################################

@doc raw"""
    intersect(I::Ideal{T}, J::Ideal{T}) where T <: RingElement

Return the intersection of the ideals `I` and `J`.
"""
function intersect(I::Ideal{T}, J::Ideal{T}) where T <: RingElement
   R = base_ring(I)
   G1 = gens(I)
   G2 = gens(J)
   if isempty(G1)
      return deepcopy(I)
   end
   if isempty(G2)
      return deepcopy(J)
   end
   return Ideal(R, lcm(G1[1], G2[1]))
end

function intersect(I::Ideal{T}, J::Ideal{T}) where {U <: FieldElement, T <: AbstractAlgebra.PolyRingElem{U}}
   R = base_ring(I)
   G1 = gens(I)
   G2 = gens(J)
   if isempty(G1)
      return deepcopy(I)
   end
   if isempty(G2)
      return deepcopy(J)
   end
   return Ideal(R, lcm(G1[1], G2[1]))
end

function intersect(I::Ideal{T}, J::Ideal{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyRingElem{U}}
   if contains(I, J)
      return J
   elseif contains(J, I)
      return I
   end
   S = base_ring(I) # poly ring
   R = base_ring(S) # coefficient ring
   # create ring with additional variable "t" with higher precedence
   tsym = gensym()
   Sup, Supv = AbstractAlgebra.polynomial_ring(R, vcat(tsym, symbols(S)); cached=false, ordering=:degrevlex)
   G1 = gens(I)
   G2 = gens(J)
   ISup = Ideal(Sup, elem_type(Sup)[f(Supv[2:end]...) for f in G1])
   JSup = Ideal(Sup, elem_type(Sup)[f(Supv[2:end]...) for f in G2])
   t = Supv[1]
   # compute t*I + (1 - t)*J
   IntSup = t*ISup + (1 - t)*JSup
   G = filter(x->!in(t, vars(x)), gens(IntSup))
   GInt = elem_type(S)[f(vcat(S(1), gens(S))...) for f in G]
   return Ideal(S, GInt)
end

function intersect(I::Ideal{T}, J::Ideal{T}) where {U <: RingElement, T <: AbstractAlgebra.MPolyRingElem{U}}
   if contains(I, J)
      return J
   elseif contains(J, I)
      return I
   end
   S = base_ring(I) # poly ring
   R = base_ring(S) # coefficient ring
   # create ring with additional variable "t" with higher precedence
   tsym = gensym()
   Sup, Supv = AbstractAlgebra.polynomial_ring(R, vcat(tsym, symbols(S)); cached=false, ordering=ordering(S))
   G1 = gens(I)
   G2 = gens(J)
   ISup = Ideal(Sup, elem_type(Sup)[f(Supv[2:end]...) for f in G1])
   JSup = Ideal(Sup, elem_type(Sup)[f(Supv[2:end]...) for f in G2])
   t = Supv[1]
   # compute t*I + (1 - t)*J
   IntSup = t*ISup + (1 - t)*JSup
   G = filter(x->!in(t, vars(x)), gens(IntSup))
   GInt = elem_type(S)[f(vcat(S(1), gens(S))...) for f in G]
   return Ideal(S, GInt)
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(I::Ideal{T}, J::Ideal{T}) where T <: RingElement
   R = base_ring(I)
   G1 = gens(I)
   G2 = gens(J)
   return Ideal(R, vcat(G1, G2))
end

function *(I::Ideal{T}, J::Ideal{T}) where T <: RingElement
   R = base_ring(I)
   G1 = gens(I)
   G2 = gens(J)
   return Ideal(R, [v*w for v in G1 for w in G2])
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

function *(I::Ideal{T}, p::T) where T <: RingElement
   R = base_ring(I)
   G = gens(I)
   if iszero(p)
      return Ideal(R, T[])
   end
   p = divexact(p, canonical_unit(p))
   return Ideal(R, [v*p for v in G])
end

function *(p::T, I::Ideal{T}) where T <: RingElement
   return I*p
end

function *(I::Ideal{T}, p::S) where {S <: RingElement, T <: RingElement}
   R = base_ring(I)
   G = gens(I)
   if iszero(p*one(R))
      return Ideal(R, T[])
   end
   V = [v*p for v in G]
   V = [divexact(v, canonical_unit(v)) for v in V]
   return Ideal(R, V)
end

function *(p::S, I::Ideal{T}) where {S <: RingElement, T <: RingElement}
   return I*p
end

###############################################################################
#
#   Ideal reduction in Euclidean domain
#
###############################################################################

function reduce_euclidean(I::Ideal{T}) where T <: RingElement
   V = gens(I)
   if length(V) > 1
      v = V[1]
      for i = 2:length(V)
         v = gcd(v, V[i])
      end
      V = [v]
   end
   if !isempty(V)
      if iszero(V[1])
         pop!(V)
      else
         V[1] = divexact(V[1], canonical_unit(V[1]))
      end
   end
   return Ideal{T}(base_ring(I), V)
end

function reduce_gens(I::Ideal{T}) where T <: RingElement
   I = reduce_euclidean(I)
   return I
end

function reduce_gens(I::Ideal{T}) where {U <: FieldElement, T <: AbstractAlgebra.PolyRingElem{U}}
   return reduce_euclidean(I)
end

###############################################################################
#
#   Ideal constructor
#
###############################################################################

function Ideal(R::Ring, V::Vector{T}) where T <: RingElement
   I = Ideal{elem_type(R)}(R, filter(!iszero, map(R, V)))
   return reduce_gens(I)
end

function Ideal(R::Ring, v::T, vs::T...) where T <: RingElement
   I = Ideal(R, filter(!iszero, [R(v), map(R, vs)...]))
   return reduce_gens(I)
end

Ideal(R::Ring) = Ideal(R, elem_type(R)[])
Ideal(R::Ring, V::Vector{Any}) = Ideal(R, elem_type(R)[R(v) for v in V])

###############################################################################
#
#   IdealSet constructor
#
###############################################################################

function IdealSet(R::Ring)
   return IdealSet{elem_type(R)}(R)
end

