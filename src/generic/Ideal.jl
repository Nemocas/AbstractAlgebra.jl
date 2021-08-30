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
#   Ideal reduction
#
###############################################################################

function reduce_euclidean(I::Ideal{T}) where T <: RingElement
end

function reduce(I::Ideal{T}) where T <: RingElement
   if hasmethod(gcdx, Tuple{T, T})
      I = reduce_euclidean(I)
   end
end

# Extend the basis V of polynomials satisfying 1, 2, 3, 4 below by the
# polynomials in D, all of which have degree at least that of those in
# V and such that the degrees of the polynomials in D is *decreasing*.
function extend_ideal_basis(D::Vector{T}, V::Vector{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   while !isempty(D)
      d = pop!(D)
@assert length(d) >= length(V[length(V)])
      V = extend_ideal_basis(d, V)
   end
   return V
end

# Given a nonempty vector V of polynomials of satisfying 1, 2, 3, 4 below and a
# polynomial p whose degree is at least that of all the polynomials in V, add
# p to V and perform reduction steps so that 1, 2, 3 and 4 still hold.
function extend_ideal_basis(p::T, V::Vector{T}) where {U <: RingElement, T <: AbstractAlgebra.PolyElem{U}}
   n = length(V)
   lc = leading_coefficient(V[n])
   # check if p can be added without any reduction
   if length(p) > length(V[n]) && !isunit(lc) && ((_, q) = divides(lc, leading_coefficient(p)))[1]
      return vcat(V, [p])
   end
   # check if p and V[n] are constant
   if isconstant(V[n]) && isconstant(p)
      return [parent(p)(gcd(constant_coefficient(p), constant_coefficient(V[n])))]
   end
   # check if p can replace V[n]
   swap = false
   if length(p) == length(V[n]) && ((_, q) = divides(lc, leading_coefficient(p)))[1]
      p, V = V[n], vcat(V[1:n-1], [p])
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
      return extend_ideal_basis(reverse(V), [p])
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
      V = vcat(V, [r])
      if length(p) == 0 # one polynomial was absorbed, yay
         return V
      end
      r = p
   end
   lenr = length(r)
   n = findfirst(x->length(x) >= lenr, V)
   if n == 1 # r is the smallest polynomial
      return extend_ideal_basis(reverse(V), [r])
   end
   return extend_ideal_basis(vcat(reverse(V[n:end]), [r]), V[1:n - 1])
end
      
# We call an ideal over a polynomial ring over a Euclidean domain reduced if
# 1. There is only one polynomial of each degree in the ideal
# 2. The degree of polynomials in the basis increases
# 3. The leading coefficient of f_i divides that of f_{i-1} for all i
# 4. Only the final polynomial may have leading coefficient that is a unit
# 5. The polynomials are all canonicalised (divided by their canonical_unit)
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
         V = [d]
         V = extend_ideal_basis(D, V) # ensure 1, 2, 3, 4 hold
      end
      # deal with 5
      V = map(x->divexact(x, canonical_unit(x)), V)
      return Ideal{T}(base_ring(I), V)
   else
      error("Not implemented")
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

