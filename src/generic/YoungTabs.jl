export Partition, AllParts, YoungTableau, SkewDiagram
export dim, has_left_neighbor, has_bottom_neighbor, inskewdiag, isrimhook,
       hooklength, leglength, matrix_repr, partitionseq

##############################################################################
#
#   Partition type, AbstractVector interface
#
##############################################################################

length(p::Partition) = length(p.part)

size(p::Partition) = size(p.part)

Base.IndexStyle(::Type{Partition}) = Base.IndexLinear()

getindex(p::Partition, i::Integer) = p.part[i]

Base.sum(p::Partition) = p.n
function setindex!(p::Partition, v::Integer, i::Integer)
   prev = sum(p)
   nex = 1
   if i != 1
      prev = p[i-1]
   end
   if i != length(p)
      nex = p[i+1]
   end
   nex <= v <= prev || throw(ArgumentError("Partition must be positive and non-increasing"))
   p.n += v - p.part[i]
   p.part[i] = v
   return p
end

==(p::Partition, m::Partition) = sum(p) == sum(m) && p.part == m.part
hash(p::Partition, h::UInt) = hash(p.part, hash(Partition, h))

##############################################################################
#
#   IO for Partition
#
##############################################################################

function subscriptify(n::Int)
   subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
   return join([Char(subscript_0 + i) for i in reverse(digits(n))])
end

function show(io::IO, p::Partition)
   uniq = unique(p.part)
   mults = [count(i -> i == u, p.part) for u in uniq]
   str = join((string(u)*subscriptify(m) for (u,m) in zip(uniq, mults)))
   print(io, str)
end

##############################################################################
#
#   Iterator interface for Integer AllParts
#
##############################################################################

const _numPartsTable = Dict{Int, Int}(0 => 1, 1 => 1, 2 => 2)
const _numPartsTableBig = Dict{Int, BigInt}()

doc"""
    _numpart(n::Integer)
> Returns the number of all distinct integer partitions of `n`. The function
> uses Euler pentagonal number theorem for recursive formula. For more details
> see OEIS sequence [A000041](http://oeis.org/A000041). Note that
> `_numpart(0) = 1` by convention.
"""
function _numpart(n::Integer)
   if n < 0
      return 0
   elseif n < 395
      n = Int(n)
      lookuptable = _numPartsTable
   else
      n = BigInt(n)
      lookuptable = _numPartsTableBig
   end
   return _numpart(n, lookuptable)
end

function _numpart(n::T, lookuptable::Dict{Int, T}) where T<:Integer
   s = zero(T)
   if !haskey(lookuptable, n)
      for j in 1:floor(T, (1 + sqrt(1+24n))/6)
         p1 = _numpart(n - div(j*(3j-1),2))
         p2 = _numpart(n - div(j*(3j+1),2))
         s += (-1)^(j-1)*(p1 + p2)
      end
      lookuptable[n] = s
   end
   return lookuptable[n]
end

# Implemented following RuleAsc (Algorithm 3.1) from
#    "Generating All Partitions: A Comparison Of Two Encodings"
# by Jerome Kelleher and Barry O’Sullivan, ArXiv:0909.2331

function Base.start(A::AllParts)
   if A.n < 1
      return 0
   elseif A.n == 1
      A.part[1] = 1
      return 0
   else
      A.part .= 0
      A.part[2] = A.n
      return 2
   end
end

Base.next(A::AllParts, state) = nextpart_asc(A.part, state)
Base.done(A::AllParts, state) = state == 1
Base.eltype(::Type{AllParts}) = Partition
length(A::AllParts) = _numpart(A.n)

@inbounds function nextpart_asc(part, k)
   if k == 0
      return Partition(part, false), 1
   end
   y = part[k] - 1
   k -= 1
   x = part[k] + 1
   while x <= y
      part[k] = x
      y -= x
      k += 1
   end
   part[k] = x + y
   return Partition(reverse!(part[1:k]), false), k
end

doc"""
    conj(part::Partition)
> Returns the conjugated partition of `part`, i.e. the partition corresponding
> to the Young tableau of `part` reflected through the main diagonal.

# Examples:
```jldoctest
julia> p = Partition([4,2,1,1,1])
4₁2₁1₃

julia> conj(p)
5₁2₁1₂
```
"""
function Base.conj(part::Partition)
   p = Vector{Int}(maximum(part))
   for i in 1:sum(part)
      for j in length(part):-1:1
         if part[j] >= i
            p[i] = j
            break
         end
      end
   end
   return Partition(p, false)
end

doc"""
    conj(part::Partition, v::Vector)
> Returns the conjugated partition of `part` together with permuted vector `v`.
"""
function Base.conj(part::Partition, v::Vector)
   w = zeros(Int, size(v))

   acc = Vector{Int}(length(part)+1)
   acc[1] = 0
   for i in 1:length(part)
      acc[i+1] = acc[i] + part[i]
   end

   new_idx = 1
   cpart = conj(part)

   for (i, p) in enumerate(cpart)
      for j in 1:p
         w[new_idx] = (v[acc[j]+i])
         new_idx += 1
      end
   end

   return cpart, w
end

##############################################################################
#
#   Partition sequences and Murnaghan-Nakayama formula
#
##############################################################################

doc"""
    partitionseq(lambda::Partition)
> Returns a sequence (as `BitVector`) of `false`s and `true`s constructed from
> `lambda`: tracing the lower contour of the Young Diagram associated to
> `lambda` from left to right a `true` is inserted for every horizontal and
> `false` for every vertical step. The sequence always starts with `true` and
> ends with `false`.
"""
function partitionseq(lambda::Partition)
   seq = trues(maximum(lambda) + length(lambda))
   j = lambda[end]
   for i in (length(lambda)-1):-1:1
      seq[j+1] = false
      j += lambda[i] - lambda[i+1] + 1
   end
   seq[j+1] = false
   return seq
end

partitionseq(v::Vector{T}) where T<:Integer = partitionseq(Partition(v))

doc"""
    partitionseq(seq::BitVector)
> Returns the essential part of the sequence `seq`, i.e. a subsequence starting
> at first `true` and ending at last `false`.
"""
partitionseq(seq::BitVector) = seq[findfirst(seq, true):findlast(seq, false)]

doc"""
    isrimhook(R::BitVector, idx::Int, len::Int)
> `R[idx:idx+len]` forms a rim hook in the Young Diagram of parition
> corresponding to `R` iff `R[idx] == true` and `R[idx+len] == false`.
"""
function isrimhook(R::BitVector, idx::Int, len::Int)
   return (R[idx+len] == false) && (R[idx] == true)
end

doc"""
    MN1inner(R::BitVector, mu::Partition, t::Int, [charvals])
> Returns the value of $\lambda$-th irreducible character on conjugacy class of
> permutations represented by partition `mu`, where `R` is the (binary)
> partition sequence representing $\lambda$. Values already computed are stored
> in `charvals::Dict{Tuple{BitVector,Vector{Int}}, Int}`.
> This is an implementation (with slight modifications) of the
> Murnaghan-Nakayama formula as described in
>
>     Dan Bernstein,
>     "The computational complexity of rules for the character table of Sn"
>     _Journal of Symbolic Computation_, 37(6), 2004, p. 727-748.
"""
function MN1inner(R::BitVector, mu::Partition, t::Integer, charvals)
   if t > length(mu)
      chi = 1
   elseif mu[t] > length(R)
      chi = 0
   else
      chi = 0
      sgn = false

      for j in 1:mu[t]-1
         if R[j] == false
            sgn = !sgn
         end
      end
      for i in 1:length(R)-mu[t]
         if R[i] != R[i+mu[t]-1]
            sgn = !sgn
         end
         if isrimhook(R, i, mu[t])
            R[i], R[i+mu[t]] = R[i+mu[t]], R[i]
            essR = (partitionseq(R), mu[t+1:end])
            if !haskey(charvals, essR)
               charvals[essR] = MN1inner(R, mu, t+1, charvals)
            end
            chi += (-1)^Int(sgn)*charvals[essR]
            R[i], R[i+mu[t]] = R[i+mu[t]], R[i]
         end
      end
   end
   return chi
end

##############################################################################
#
#   YoungTableau type, AbstractVector interface
#
##############################################################################

function YoungTableau(part::Partition, fill::Vector{Int}=collect(1:part.n))
   part.n == length(fill) || throw("Can't fill Young digaram of $part with $fill: different number of elemnets.")
   tab = zeros(Int, length(part), maximum(part))
   k=1
   for (idx, p) in enumerate(part)
      tab[idx, 1:p] = fill[k:k+p-1]
      k += p
   end
   return YoungTableau(part, tab)
end

YoungTableau(p::Vector{T}, fill=collect(1:sum(p))) where T<:Integer = YoungTableau(Partition(p), fill)

size(Y::YoungTableau) = size(Y.tab)

Base.IndexStyle(::Type{T}) where {T <: YoungTableau} = Base.IndexLinear()

getindex(Y::YoungTableau, i::Int) = Y.tab[i]

function ==(Y1::YoungTableau,Y2::YoungTableau)
   Y1.part == Y2.part || return false
   Y1.tab == Y2.tab || return false
   return true
end

hash(Y::YoungTableau, h::UInt) = hash(Y.part, hash(Y.tab, hash(YoungTableau, h)))

doc"""
    conj(Y::YoungTableau)
> Returns the conjugated tableau, i.e. the tableau reflected through the main
> diagonal.
"""
conj(Y::YoungTableau) = YoungTableau(conj(Y.part), transpose(Y.tab))

##############################################################################
#
#   Misc functions for YoungTableaux
#
##############################################################################

rowlen(Y::YoungTableau, i, j) = sum(Y[i, j:end] .> 0)
collen(Y::YoungTableau, i, j) = sum(Y[i:end, j] .> 0)

doc"""
    hooklength(Y::YoungTableau, i, j)
> Returns the hooklength of an element in `Y` at position `(i,j)`. `hooklength`
> will return `0` for `(i,j)` not in the tableau `Y`.
"""
function hooklength(Y::YoungTableau, i, j)
   if Y[i,j] == 0
      return 0
   else
      return rowlen(Y, i, j) + collen(Y, i, j) - 1
   end
end

doc"""
    dim(Y::YoungTableau)
> Returns the dimension of the irreducible representation of
> `PermutationGroup(n)` associated to YoungTableau `Y` of a partition of `n`.
"""
dim(Y::YoungTableau) = dim(BigInt, Y)

function dim(::Type{T}, Y::YoungTableau) where T<:Integer
   n, m = size(Y)
   num = factorial(T(maximum(Y)))
   den = reduce(*, one(T), hooklength(Y,i,j) for i in 1:n, j in 1:m if j <= Y.part[i])
   return divexact(num, den)
end

##############################################################################
#
#   SkewDiagrams
#
##############################################################################

SkewDiagram(lambda::Vector{Int}, mu::Vector{Int}) = SkewDiagram(Partition(lambda), Partition(mu))

/(lambda::Partition, mu::Partition) = SkewDiagram(lambda, mu)

==(xi::SkewDiagram, psi::SkewDiagram) = xi.lam == psi.lam && xi.mu == psi.mu
hash(xi::SkewDiagram, h::UInt) = hash(xi.lam, hash(xi.mu, hash(SkewDiagram, h)))

###############################################################################
#
#   String I/O
#
###############################################################################

doc"""
    matrix_repr(xi::SkewDiagram)
> Returns a binary representation of the diagram `xi`, i.e. a binary array `A`
> where `A[i,j] == 1` if and only if `(i,j)` is in `xi.lam` but not in `xi.mu`.
"""
function matrix_repr(xi::SkewDiagram)
   ydiag = zeros(Int,length(xi.lam), maximum(xi.lam))
   for i in 1:length(xi.mu)
      ydiag[i, xi.mu[i]+1:xi.lam[i]] .= 1
   end
   for i in length(xi.mu)+1:length(xi.lam)
      ydiag[i,1:xi.lam[i]] .= 1
   end
   return ydiag
end

show(io::IO, xi::SkewDiagram) = show(io, MIME("text/plain"), matrix_repr(xi))

##############################################################################
#
#   Misc functions for SkewDiagrams
#
##############################################################################

doc"""
    inskewdiag(xi::SkewDiagram, i::Int, j::Int)
> Checks if box at position `(i,j)` belongs to the skew diagram `xi`.
"""
function inskewdiag(xi::SkewDiagram, i::Int, j::Int)
   if i <= 0 || j <= 0
      return false
   elseif i > length(xi.lam) || j > maximum(xi.lam)
      return false
   elseif length(xi.mu) >= i
      return xi.mu[i] < j ≤ xi.lam[i]
   else
      return j ≤ xi.lam[i]
   end
end

doc"""
    has_left_neighbor(xi::SkewDiagram, i::Int, j::Int)
> Checks if box at position `(i,j)` has neighbour in `xi` to the left.
"""
function has_left_neighbor(xi::SkewDiagram, i::Int, j::Int)
   if j == 1
      return false
   else
      return inskewdiag(xi, i, j) && inskewdiag(xi, i, j-1)
   end
end

doc"""
    has_bottom_neighbor(xi::SkewDiagram, i::Int, j::Int)
> Checks if box at position `(i,j)` has neighbour in `xi` below.
"""
function has_bottom_neighbor(xi::SkewDiagram, i::Int, j::Int)
   if i == length(xi.lam)
      return false
   else
      return inskewdiag(xi, i, j) && inskewdiag(xi, i+1, j)
   end
end

doc"""
    isrimhook(xi::SkewDiagram)
> Checks if `xi` represents a rim-hook diagram, i.e. its diagram is
> edge-connected and contains no $2\times 2$ squares.
"""
function isrimhook(xi::SkewDiagram)
   i = 1
   j = xi.lam[1]
   while i != length(xi.lam) && j != 1
      left = has_left_neighbor(xi, i,j)
      bottom = has_bottom_neighbor(xi, i,j)
      if left && bottom # there is 2×2 square in xi
         return false
      elseif left
         j -= 1
      elseif bottom
         i += 1
      else
         lam_tail = xi.lam[i+1:end]
         mu_tail = zeros(Int, length(lam_tail))
         mu_tail[1:length(xi.mu)-i] = xi.mu[i+1:end]

         if any(lam_tail .- mu_tail .> 0)
            return false # xi is disconnected
         else
            return true # we arrived at the end of xi
         end
      end
   end
   return true
end

doc"""
    leglength(xi::SkewDiagram, check::Bool=true)
> Computes the leglength of a rim-hook `xi`, i.e. the number of rows with
> non-zero entries minus one. If `check` is `false` function will not check
> whether `xi` is actually a rim-hook.
"""
function leglength(xi::SkewDiagram, check::Bool=true)
   if check
      isrimhook(xi) || throw("$xi is not a rimhook. leglength is defined only for rim hooks")
   end
   m = zeros(length(xi.lam))
   m[1:length(xi.mu)] = xi.mu
   return sum((xi.lam .- m) .> 0) - 1
end
