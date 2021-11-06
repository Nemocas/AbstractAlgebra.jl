###############################################################################
#
#   FreeAssAlgebra.jl : free associative algebra R<x1,...,xn>
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

function parent_type(::Type{FreeAssAlgElem{T}}) where T <: RingElement
   return FreeAssAlgebra{T}
end

function elem_type(::Type{FreeAssAlgebra{T}}) where T <: RingElement
   return FreeAssAlgElem{T}
end

function parent(a::FreeAssAlgElem)
   return a.parent
end

function base_ring(a::FreeAssAlgebra)
   return a.base_ring
end

function symbols(a::FreeAssAlgebra)
   return a.S
end

function nvars(a::FreeAssAlgebra)
   return length(a.S)
end

function length(a::FreeAssAlgElem)
   return a.length
end

function check_parent(a::FreeAssAlgElem{T}, b::FreeAssAlgElem{T}, throw::Bool = true) where T <: RingElement
   b = parent(a) != parent(b)
   b & throw && error("Incompatible polynomial rings in polynomial operation")
   return !b
end

###############################################################################
#
# Basic manipulation
#
###############################################################################

function Base.deepcopy_internal(a::FreeAssAlgElem{T}, dict::IdDict) where {T <: RingElement}
   return FreeAssAlgElem{T}(a.parent,
                           deepcopy_internal(a.coeffs, dict),
                           deepcopy_internal(a.exps, dict),
                           a.length)
end

function zero(a::FreeAssAlgebra{T}) where T
   return FreeAssAlgElem{T}(a, T[], Vector{Int}[], 0)
end

function one(a::FreeAssAlgebra{T}) where T
   c = one(base_ring(a))
   !iszero(c) || return zero(a)
   return FreeAssAlgElem{T}(a, [c], [Int[]], 1)
end

function iszero(a::FreeAssAlgElem{T}) where T
   return length(a) == 0
end

function isone(a::FreeAssAlgElem{T}) where T
   if length(a) == 0
      return isone(zero(base_ring(a)))
   else
      return a.length == 1 && isone(a.coeffs[1]) && isempty(a.exps[1])
   end
end

function gen(a::FreeAssAlgebra{T}, i::Int) where T
   0 < i <= nvars(a) || error("variable index out of range")
   c = one(base_ring(a))
   iszero(c) && return zero(a)
   return FreeAssAlgElem{T}(a, T[c], [Int[i]], 1)
end

function gens(a::FreeAssAlgebra{T}) where {T <: RingElement}
   return [gen(a, i) for i in 1:nvars(a)]
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::FreeAssAlgebra{T})(b::T) where T
   iszero(b) && return zero(a)
   return FreeAssAlgElem{T}(a, T[b], [Int[]], 1)
end

function (a::FreeAssAlgebra{T})(b::Integer) where T
   return a(base_ring(a)(b))
end

function (a::FreeAssAlgebra{T})(c::Vector{T}, e::Vector{Vector{Int}}) where T
   for ei in e
      all(i -> (i <= nvars(a)), ei) || error("variable index out of range")
   end
   n = length(c)
   n == length(e) || error("coefficient array and exponent array should have the same length")
   return FreeAssAlgelem{T}(a, copy(c), copy(e), n)
end

###############################################################################
#
# Coefficients, Terms, Etc.
#
###############################################################################

function coeff(a::FreeAssAlgElem, i::Int)
   return a.coeffs[i]
end

function term(a::FreeAssAlgElem{T}, i::Int) where T <: RingElement
   R = parent(a)
   return FreeAssAlgElem{T}(R, a.coeffs[i], [a.exps[i]], 1)
end

function monomial(a::FreeAssAlgElem{T}, i::Int) where T <: RingElement
   R = parent(a)
   return FreeAssAlgElem{T}(R, T[one(R)], [a.exps[i]], 1)
end

function exponent_word(a::FreeAssAlgElem{T}, i::Int) where T <: RingElement
   0 < i <= length(a) || error("index out of range")
   return a.exps[i]
end

function Base.iterate(a::FreeAssAlgExponentWords, state = 0)
   state += 1
   state <= length(a.poly) || return nothing
   return exponent_word(a.poly, state), state
end

function leading_coefficient(a::FreeAssAlgElem{T}) where T
   a.length > 0 || return zero(base_ring(a))
   return a.coeffs[1]
end

function leading_word(a::FreeAssAlgElem{T}) where T
   a.length > 0 || error("element is zero")
   return a.exps[1]
end

function total_degree(a::FreeAssAlgElem{T}) where T
   # currently stored in dexlex
   return length(a) > 0 ? length(leading_word(a)) : -1
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::FreeAssAlgElem{T}, b::FreeAssAlgElem{T}) where T
   fl = check_parent(a, b, false)
   !fl && return false
   return (a.length == b.length) && 
          (view(a.coeffs, 1:a.length) == view(b.coeffs, 1:b.length)) &&
          (view(a.exps, 1:a.length) == view(b.exps, 1:b.length))
end

function word_cmp(a::Vector{Int}, b::Vector{Int})
   if length(a) > length(b)
      return +1
   elseif length(a) < length(b)
      return -1
   else
      # deglex
      for i in 1:length(a)
         if a[i] > b[i]
            return -1
         elseif a[i] < b[i]
            return +1
         end
      end
      return 0
   end
end

function word_gt(a::Vector{Int}, b::Vector{Int})
   return word_cmp(a, b) > 0
end

function sort_terms!(z::FreeAssAlgElem{T}) where T
   p = sortperm(z.exps, lt = word_gt)
   z.coeffs = [z.coeffs[p[i]] for i in 1:length(p)]
   z.exps = [z.exps[p[i]] for i in 1:length(p)]
   return z
end

function combine_like_terms!(z::FreeAssAlgElem{T}) where T
   o = 0
   i = 1
   while i <= z.length
      if o > 0 && word_cmp(z.exps[o], z.exps[i]) == 0
         z.coeffs[o] += z.coeffs[i]
      else
         o += (o < 1 || !iszero(z.coeffs[o]))
         z.exps[o] = z.exps[i]
         z.coeffs[o] = z.coeffs[i]
      end
      i += 1
   end
   o += (o < 1 || !iszero(z.coeffs[o]))
   z.length = o - 1
   return z
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function *(a::FreeAssAlgElem{T}, b::FreeAssAlgElem{T}) where T
   zcoeffs = T[]
   zexps = Vector{Int}[]
   for i in 1:a.length, j in 1:b.length
      push!(zcoeffs, a.coeffs[i]*b.coeffs[j])
      push!(zexps, vcat(a.exps[i], b.exps[j]))
   end
   z = FreeAssAlgElem{T}(parent(a), zcoeffs, zexps, length(zcoeffs))
   return combine_like_terms!(sort_terms!(z))
end

function +(a::FreeAssAlgElem{T}, b::FreeAssAlgElem{T}) where T
   zcoeffs = T[]
   zexps = Vector{Int}[]
   i = j = 1
   while i <= a.length && j <= b.length
      c = word_cmp(a.exps[i], b.exps[j])
      if c < 0
         push!(zcoeffs, b.coeffs[j])
         push!(zexps, b.exps[j])
         j += 1
      elseif c > 0
         push!(zcoeffs, a.coeffs[i])
         push!(zexps, a.exps[i])
         i += 1
      else
         s = a.coeffs[i] + b.coeffs[j]
         if !iszero(s)
            push!(zcoeffs, s)
            push!(zexps, a.exps[i])
         end
         i += 1
         j += 1
      end
   end
   while i <= a.length
      push!(zcoeffs, a.coeffs[i])
      push!(zexps, a.exps[i])
      i += 1
   end
   while j <= b.length
      push!(zcoeffs, b.coeffs[j])
      push!(zexps, b.exps[j])
      j += 1
   end
   return FreeAssAlgElem{T}(parent(a), zcoeffs, zexps, length(zcoeffs))
end

function -(a::FreeAssAlgElem{T}, b::FreeAssAlgElem{T}) where T
   zcoeffs = T[]
   zexps = Vector{Int}[]
   i = j = 1
   while i <= a.length && j <= b.length
      c = word_cmp(a.exps[i], b.exps[j])
      if c < 0
         push!(zcoeffs, -b.coeffs[j])
         push!(zexps, b.exps[j])
         j += 1
      elseif c > 0
         push!(zcoeffs, a.coeffs[i])
         push!(zexps, a.exps[i])
         i += 1
      else
         s = a.coeffs[i] - b.coeffs[j]
         if !iszero(s)
            push!(zcoeffs, s)
            push!(zexps, a.exps[i])
         end
         i += 1
         j += 1
      end
   end
   while i <= a.length
      push!(zcoeffs, a.coeffs[i])
      push!(zexps, a.exps[i])
      i += 1
   end
   while j <= b.length
      push!(zcoeffs, -b.coeffs[j])
      push!(zexps, b.exps[j])
      j += 1
   end
   return FreeAssAlgElem{T}(parent(a), zcoeffs, zexps, length(zcoeffs))
end

function ^(a::FreeAssAlgElem{T}, b::Integer) where T
   if b == 0
      return one(parent(a))
   elseif b == 1
      return deepcopy(a)
   elseif a.length == 1
      if isempty(a.exps[1])
         e = [Int[]]
      else
         b < 0 && throw(NotInvertibleError(a))
         e = [reduce(vcat, [a.exps[1] for i in 1:b])]
      end
      return FreeAssAlgElem{T}(parent(a), [a.coeffs[1]^b], e, 1)
   else
      b < 0 && throw(NotInvertibleError(a))
      return AbstractAlgebra.internal_power(a, b)
   end
end

###############################################################################
#
# Division
#
###############################################################################

# return c*w*a*wp
function mul_term(c::T, w::Vector{Int}, a::FreeAssAlgElem{T}, wp::Vector{Int}) where T
   zcoeffs = isone(c) ? T[a.coeffs[i] for i in 1:a.length] :
                        T[c*a.coeffs[i] for i in 1:a.length]
   zexps = Vector{Int}[vcat(w, a.exps[i], wp) for i in 1:a.length]
   return FreeAssAlgElem{T}(parent(a), zcoeffs, zexps, a.length)
end


# return (true, l, r) with a = l*b*r and length(l) minimal
#     or (false, junk, junk)
function word_divides_leftmost(a::Vector{Int}, b::Vector{Int})
   n = length(b)
   for i in 0:length(a)-n
      match = true
      for j in 1:n
         if b[j] != a[i+j]
            match = false
            break
         end
      end
      if match
         return (true, Int[a[k] for k in 1:i],
                       Int[a[k] for k in 1+i+n:length(a)])
      end
   end
   return (false, Int[], Int[])
end

function word_divides_rightmost(a::Vector{Int}, b::Vector{Int})
   n = length(b)
   for i in length(a)-n:-1:0
      match = true
      for j in 1:n
         if b[j] != a[i+j]
            match = false
            break
         end
      end
      if match
         return (true, Int[a[k] for k in 1:i],
                       Int[a[k] for k in 1+i+n:length(a)])
      end
   end
   return (false, Int[], Int[])
end


function AbstractAlgebra.divexact_left(f::FreeAssAlgElem{T}, g::FreeAssAlgElem{T}; check::Bool = true) where T
   R = parent(f)
   qcoeffs = T[]
   qexps = Vector{Int}[]
   while length(f) > 0
      ok, ml, mr = word_divides_leftmost(f.exps[1], g.exps[1])
      ok && isempty(ml) || error("not exact division")
      qi = divexact(f.coeffs[1], g.coeffs[1])
      push!(qcoeffs, qi)
      push!(qexps, mr)
      f -= mul_term(qi, ml, g, mr)
   end
   return FreeAssAlgElem{T}(R, qcoeffs, qexps, length(qcoeffs))
end

function AbstractAlgebra.divexact_right(f::FreeAssAlgElem{T}, g::FreeAssAlgElem{T}; check::Bool = true) where T
   R = parent(f)
   qcoeffs = T[]
   qexps = Vector{Int}[]
   while length(f) > 0
      ok, ml, mr = word_divides_rightmost(f.exps[1], g.exps[1])
      ok && isempty(mr) || error("not exact division")
      qi = divexact(f.coeffs[1], g.coeffs[1])
      push!(qcoeffs, qi)
      push!(qexps, ml)
      f -= mul_term(qi, ml, g, mr)
   end
   return FreeAssAlgElem{T}(R, qcoeffs, qexps, length(qcoeffs))
end


###############################################################################
#
#   FreeAssociativeAlgebra constructor
#
###############################################################################

function FreeAssociativeAlgebra(R::AbstractAlgebra.Ring, s::Vector{Symbol}; cached::Bool = true)
   parent_obj = FreeAssAlgebra{elem_type(R)}(R, s, cached)
   return (parent_obj, gens(parent_obj))
end

