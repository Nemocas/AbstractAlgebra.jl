###############################################################################
#
#   PermGroup / perm
#
###############################################################################

const PermID = ObjectIdDict()

type PermGroup <: Group
   n::Int

   function PermGroup(n::Int, cached=true)
      if haskey(PermID, n)
         return PermID[n]::PermGroup
      else
         z = new(n)
         if cached
            PermID[n] = z
         end
         return z
      end
   end
end

type perm <: GroupElem
   d::Array{Int, 1}
   parent::PermGroup

   function perm(n::Int)
      return new(collect(1:n))
   end

   function perm(a::Array{Int, 1})
      return new(a)
   end
end

export PermGroup, perm, parity, elements, cycles

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent_type(::Type{perm}) = PermGroup

elem_type(R::PermGroup) = perm

###############################################################################
#
#   Basic manipulation
#
###############################################################################

doc"""
    parent(a::perm)
>Return the parent of the given permutation group element.

"""
parent(a::perm) = a.parent

function deepcopy_internal(a::perm, dict::ObjectIdDict)
   R = parent(a)
   return R(deepcopy(a.d))
end

function Base.hash(a::perm, h::UInt)
   b = 0x595dee0e71d271d0%UInt
   for i in 1:length(a.d)
      b $= hash(a[i], h) $ h
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

doc"""
    parity(a::perm)
> Return the parity of the given permutation, i.e. the parity of the number of
> transpositions that compose it. The function returns $1$ if the parity is odd
> otherwise it returns $0$.

"""
# TODO: 2x slower than Flint
function parity(a::perm)
   to_visit = trues(a.d)
   parity = length(to_visit)
   k = 1
   while any(to_visit)
      parity -= 1
      k = findnext(to_visit, k)
      to_visit[k] = false
      next = a[k]
      while next != k
         to_visit[next] = false
         next = a[next]
      end
   end
   return parity%2
end

function getindex(a::perm, n::Int)
   return a.d[n]
end

function setindex!(a::perm, d::Int, n::Int)
   a.d[n] = d
end

doc"""
    eye(R::PermGroup)
> Return the identity permutation for the given permutation group.

"""
eye(R::PermGroup) = R()

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::PermGroup)
   print(io, "Permutation group over $(R.n) elements")
end

type DisplayStyle
   format::Symbol
end

const perm_display_style = DisplayStyle(:array)

function setpermstyle(format::Symbol)
   format in (:array, :cycles) || throw("Permutations can be displayed
   only in styles :array or :cycles")
   perm_display_style.format = format
   return  format
end

function show(io::IO, p::perm)
   if perm_display_style.format == :array
      print(io, "[" * join(p.d, ", ") * "]")
   elseif perm_display_style.format == :cycles
      if p == parent(p)()
         print(io, "()")
      else
         print(io, join(["("*join(c, ",")*")" for c in cycles(p) if
            length(c)>1],""))
      end
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

doc"""
    ==(a::perm, b::perm)
> Return `true` if the given permutations are equal, otherwise return `false`.

"""
function ==(a::perm, b::perm)
   parent(a) == parent(b) || return false
   a.d == b.d || return false
   return true
end

###############################################################################
#
#   Binary operators
#
###############################################################################

doc"""
> Return the composition of the two permutations, i.e. $a\circ b$. In other
> words, the permutation corresponding to applying $b$ first, then $a$, is
> returned.

"""
function *(a::perm, b::perm)
   d = similar(a.d)
   for i in 1:length(d)
      d[i] = a[b[i]]
   end
   R = parent(a)
   return R(d)
end

###############################################################################
#
#   Inversion
#
###############################################################################

doc"""
    inv(a::perm)
> Return the inverse of the given permutation, i.e. the permuation $a^{-1}$
> such that $a\circ a^{-1} = a^{-1}\circ a$ is the identity permutation.

"""
function inv(a::perm)
   d = similar(a.d)
   for i in 1:length(a.d)
      d[a[i]] = i
   end
   R = parent(a)
   return R(d)
end

# TODO: can we do that in place??
function inv!(a::perm)
   d = similar(a.d)
   for i in 1:length(a.d)
      d[a[i]] = i
   end
   a.d = d
end

###############################################################################
#
#   Iterating over all permutations
#
###############################################################################

doc"""
   AllPerms(n::Int)
> Returns an iterator over arrays representing all permutations of `1:n`.
> Similar to `Combinatorics.permutations(1:n)`

"""
immutable AllPerms
   n::Int
   all::Int
   AllPerms(n::Int) = new(n, factorial(n))
end

Base.start(A::AllPerms) = (collect(1:A.n), 1, 1, ones(Int, A.n))
Base.next(A::AllPerms, state) = all_perms(state...)
Base.done(A::AllPerms, state) = state[2] > A.all
Base.eltype(::AllPerms) = Vector{Int}
Base.length(A::AllPerms) = factorial(A.n)

function all_permutations(elts, counter, i, c)
   if counter == 1
      return (copy(elts), (elts, counter+1, i, c))
   end
   n = length(elts)
   while i <= n
       if c[i] < i
          if isodd(i)
             elts[1], elts[i] = elts[i], elts[1]
          else
             elts[c[i]], elts[i] = elts[i], elts[c[i]]
          end
          c[i] += 1
          i = 1
          return (copy(elts), (elts, counter+1, i, c))
       else
          c[i] = 1
          i += 1
       end
    end
end

doc"""
    elements(G::PermGroup)
> Returns an iterator over all elements in the group $G$. You may use
> `collect(elements(G))` to get an array with all elements.

"""
elements(G::PermGroup) = (G(p) for p in AllPerms(G.n))

###############################################################################
#
#   Misc
#
###############################################################################

doc"""
    order(G::PermGroup)
> Returns the order of the full permutation group.

"""
order(G::PermGroup) = factorial(G.n)

doc"""
    cycles(a::perm)
> Decomposes permutation into disjoint cycles.

"""
function cycles(a::perm)
   to_visit = trues(a.d)
   cycles = Vector{Vector{Int}}()
   k = 1
   while any(to_visit)
      cycle = Vector{Int}()
      k = findnext(to_visit, k)
      to_visit[k] = false
      push!(cycle, k)
      next = a[k]
      while next != k
         push!(cycle, next)
         to_visit[next] = false
         next = a[next]
      end
      push!(cycles, cycle)
   end
   return cycles
end

doc"""
    matrix_repr(a::perm)
> Return the permutation matrix representing `a` via natural embedding of the
> permutation group into general linear group over ZZ

"""
function matrix_repr(a::perm)
   A = eye(Int, parent(a).n)
   return A[a.d,:]
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (R::PermGroup)()
   z = perm(R.n)
   z.parent = R
   return z
end

function (R::PermGroup)(a::Array{Int, 1})
   length(a) != R.n && error("Unable to coerce to permutation")
   Base.Set(a) != Base.Set(1:length(a)) && error("Unable to coerce to permutation")

   z = perm(a)
   z.parent = R
   return z
end

function (R::PermGroup)(a::perm)
   parent(a) != R && error("Unable to coerce to permutation")
   return a
end

###############################################################################
#
#   PermGroup constructor
#
###############################################################################

# handled by inner constructors
