###############################################################################
#
#   SparsePoly.jl : Generic sparse univariate polynomials over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{SparsePoly{T}}) where {T <: RingElement} = SparsePolyRing{T}

elem_type(::Type{SparsePolyRing{T}}) where {T <: RingElement} = SparsePoly{T}

function is_domain_type(a::Type{SparsePoly{T}}) where T <: RingElement
   return is_domain_type(T)
end

function is_exact_type(a::Type{SparsePoly{T}}) where T <: RingElement
   return is_exact_type(T)
end

var(a::SparsePolyRing) = a.S

function gen(a::SparsePolyRing)
   return a([one(base_ring(a))], [UInt(1)])
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function coeff(x::SparsePoly, i::Int)
   i < 0 && throw(DomainError(i, "cannot get the i-th coefficient with i < 0"))
   return x.coeffs[i + 1]
end

one(R::SparsePolyRing) = R(1)

zero(R::SparsePolyRing) = R(0)

iszero(x::SparsePoly) = length(x) == 0

isone(x::SparsePoly) = x == 1

length(x::SparsePoly) = x.length

leading_coefficient(x::SparsePoly) = x.length == 0 ? base_ring(x)() : x.coeffs[x.length]

trailing_coefficient(x::SparsePoly) = x.length == 0 ? base_ring(x)() : x.coeffs[1]

function normalise(a::SparsePoly, n::Int)
   while n > 0 && iszero(a.coeffs[n])
      n -= 1
   end
   return n
end

base_ring_type(::Type{SparsePolyRing{T}}) where T <: RingElement = parent_type(T)

base_ring(R::SparsePolyRing{T}) where {T <: RingElement} = R.base_ring::parent_type(T)

parent(a::SparsePoly) = a.parent

function Base.deepcopy_internal(a::SparsePoly{T}, dict::IdDict) where {T <: RingElement}
   Re = Base.deepcopy_internal(a.exps, dict)
   Rc = Vector{T}(undef, a.length)
   for i = 1:a.length
      Rc[i] = Base.deepcopy_internal(a.coeffs[i], dict)
   end
   return parent(a)(Rc, Re)
end

function characteristic(a::SparsePolyRing{T}) where T <: RingElement
   return characteristic(base_ring(a))
end

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::SparsePoly, x = var(parent(a)); context = nothing)
   sum = Expr(:call, :+)
   for i in length(a):-1:1 # the polynomials are stored backwards
      c = coeff(a, i-1)    # ???
      e = a.exps[i]
      xe = iszero(e) ? 1 : isone(e) ? x : Expr(:call, :^, x, e)
      if isone(c)
          push!(sum.args, Expr(:call, :*, xe))
      else
          push!(sum.args, Expr(:call, :*, expressify(c, context = context), xe))
      end
   end
   return sum
end

function Base.show(io::IO, ::MIME"text/plain", a::SparsePoly)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function Base.show(io::IO, a::SparsePoly)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function show(io::IO, p::SparsePolyRing)
   print(io, "Sparse Univariate Polynomial Ring in ")
   print(io, string(p.S))
   print(io, " over ")
   print(IOContext(io, :compact => true), base_ring(p))
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function -(a::SparsePoly{T}) where {T <: RingElement}
   r = parent(a)()
   fit!(r, length(a))
   for i = 1:length(a)
      r.coeffs[i] = -a.coeffs[i]
      r.exps[i] = a.exps[i]
   end
   r.length = a.length
   return r
end

function +(a::SparsePoly{T}, b::SparsePoly{T}) where {T <: RingElement}
   par = parent(a)
   r = par()
   fit!(r, length(a) + length(b))
   i = 1
   j = 1
   k = 1
   while i <= length(a) && j <= length(b)
      if a.exps[i] < b.exps[j]
         r.coeffs[k] = a.coeffs[i]
         r.exps[k] = a.exps[i]
         i += 1
      elseif a.exps[i] == b.exps[j]
         c = a.coeffs[i] + b.coeffs[j]
         if c != 0
            r.coeffs[k] = c
            r.exps[k] = a.exps[i]
         else
            k -= 1
         end
         i += 1
         j += 1
      else
         r.coeffs[k] = b.coeffs[j]
         r.exps[k] = b.exps[j]
         j += 1
      end
      k += 1
   end
   while i <= length(a)
      r.coeffs[k] = a.coeffs[i]
      r.exps[k] = a.exps[i]
      i += 1
      k += 1
   end
   while j <= length(b)
      r.coeffs[k] = b.coeffs[j]
      r.exps[k] = b.exps[j]
      j += 1
      k += 1
   end
   r.length = k - 1
   return r
end

function -(a::SparsePoly{T}, b::SparsePoly{T}) where {T <: RingElement}
   par = parent(a)
   r = par()
   fit!(r, length(a) + length(b))
   i = 1
   j = 1
   k = 1
   while i <= length(a) && j <= length(b)
      if a.exps[i] < b.exps[j]
         r.coeffs[k] = a.coeffs[i]
         r.exps[k] = a.exps[i]
         i += 1
      elseif a.exps[i] == b.exps[j]
         c = a.coeffs[i] - b.coeffs[j]
         if c != 0
            r.coeffs[k] = c
            r.exps[k] = a.exps[i]
         else
            k -= 1
         end
         i += 1
         j += 1
      else
         r.coeffs[k] = -b.coeffs[j]
         r.exps[k] = b.exps[j]
         j += 1
      end
      k += 1
   end
   while i <= length(a)
      r.coeffs[k] = a.coeffs[i]
      r.exps[k] = a.exps[i]
      i += 1
      k += 1
   end
   while j <= length(b)
      r.coeffs[k] = -b.coeffs[j]
      r.exps[k] = b.exps[j]
      j += 1
      k += 1
   end
   r.length = k - 1
   return r
end

function *(a::SparsePoly{T}, b::SparsePoly{T}) where {T <: RingElement}
   m = length(a)
   n = length(b)
   if m == 0 || n == 0
      return parent(a)()
   end
   if m < n
      return b*a
   end
   r = a*b.coeffs[1]
   r = shift_left!(r, b.exps[1])
   for i = 2:n
      s = a*b.coeffs[i]
      s = shift_left!(s, b.exps[i])
      r += s
   end
   return r
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

function *(a::SparsePoly, n::Union{Integer, Rational, AbstractFloat})
   r = parent(a)()
   fit!(r, length(a))
   j = 1
   for i = 1:length(a)
      c = a.coeffs[i]*n
      if c != 0
         r.coeffs[j] = c
         r.exps[j] = a.exps[i]
         j += 1
      end
   end
   r.length = j - 1
   return r
end

function *(a::SparsePoly{T}, n::T) where {T <: RingElem}
   r = parent(a)()
   fit!(r, length(a))
   j = 1
   for i = 1:length(a)
      c = a.coeffs[i]*n
      if c != 0
         r.coeffs[j] = c
         r.exps[j] = a.exps[i]
         j += 1
      end
   end
   r.length = j - 1
   return r
end

*(n::T, a::SparsePoly{T}) where {T <: RingElem} = a*n

*(n::Union{Integer, Rational, AbstractFloat}, a::SparsePoly) = a*n

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::SparsePoly{T}, b::SparsePoly{T}) where {T <: RingElement}
   if length(a) != length(b)
      return false
   end
   for i = 1:length(a)
      if a.exps[i] != b.exps[i] || a.coeffs[i] != b.coeffs[i]
         return false
      end
   end
   return true
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function ==(a::SparsePoly, b::Union{Integer, Rational, AbstractFloat})
   return length(a) == 0 ? b == 0 : a.length == 1 &&
          a.exps[1] == 0 && a.coeffs[1] == b
end

==(a::Union{Integer, Rational, AbstractFloat}, b::SparsePoly) = b == a

function ==(a::SparsePoly{T}, b::T) where T <: RingElem
   return length(a) == 0 ? iszero(b) : a.length == 1 &
          a.exps[1] == 0 && a.coeffs[1] == b
end

==(a::T, b::SparsePoly{T}) where T <: RingElem = b == a

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::SparsePoly{T}, b::Int) where {T <: RingElement}
   b < 0 && throw(DomainError(b, "exponent must be >= 0"))
   # special case powers of x for constructing polynomials efficiently
   if length(a) == 0
      return parent(a)()
   elseif length(a) == 1
      return parent(a)([coeff(a, 0)^b], [a.exps[1]*b])
   elseif b == 0
      return one(parent(a))
   elseif b == 1
      return a
   elseif b == 2
      return a*a
   else
      z = a*a
      for i = 3:b
         z *= a
      end
      return z
   end
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divides(a::SparsePoly{T}, b::SparsePoly{T}) where {T <: RingElement}
   d1 = a.exps[a.length]
   d2 = b.exps[b.length] - b.exps[1]
   q_alloc = b.length
   Qe = Vector{UInt}(undef, q_alloc)
   Qc = Vector{T}(undef, q_alloc)
   r = a
   k = 0
   while length(r) > 0
      if r.exps[1] + d2 > d1
         return false, parent(a)()
      end
      flag, c = divides(r.coeffs[1], b.coeffs[1])
      if !flag
         return false, parent(a)()
      end
      s = b*c
      if r.exps[1] < b.exps[1]
         return false, parent(a)()
      end
      d = r.exps[1] - b.exps[1]
      k += 1
      if k > q_alloc
         q_alloc *= 2
         resize!(Qe, q_alloc)
         resize!(Qc, q_alloc)
      end
      Qe[k] = d
      Qc[k] = c
      s = shift_left!(s, d)
      r -= s
   end
   resize!(Qe, k)
   resize!(Qc, k)
   return true, parent(a)(Qc, Qe)
end

function divexact(a::SparsePoly{T}, b::SparsePoly{T}; check::Bool=true) where {T <: RingElement}
   d, q = divides(a, b)
   check && d == false && error("Not an exact division in divexact")
   return q
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divides(a::SparsePoly{T}, b::T) where {T <: RingElem}
   len = a.length
   Qc = Vector{T}(undef, len)
   for i = 1:len
      flag, Qc[i] = divides(a.coeffs[i], b)
      if !flag
         return false, parent(a)()
      end
   end
   return true, parent(a)(Qc, a.exps)
end

function divexact(a::SparsePoly{T}, b::T; check::Bool=true) where {T <: RingElem}
   len = length(a)
   exps = deepcopy(a.exps)
   coeffs = [divexact(a.coeffs[i], b; check=check) for i in 1:len]
   return parent(a)(coeffs, exps)
end

function divexact(a::SparsePoly, b::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
   len = length(a)
   exps = deepcopy(a.exps)
   coeffs = [divexact(a.coeffs[i], b; check=check) for i in 1:len]
   return parent(a)(coeffs, exps)
end

###############################################################################
#
#   Pseudodivision
#
###############################################################################

function pseudorem(a::SparsePoly{T}, b::SparsePoly{T}) where {T <: RingElement}
   par = parent(a)
   R = base_ring(par)
   m = length(a)
   n = length(b)
   n == 0 && throw(DivideError())
   if a.exps[a.length] < b.exps[b.length]
      return deepcopy(a)
   end
   k = reinterpret(Int, a.exps[a.length] - b.exps[b.length]) + 1
   l = leading_coefficient(b)
   while a.length > 0 && a.exps[a.length] >= b.exps[b.length]
      s = leading_coefficient(a)*b
      s = shift_left!(s, a.exps[a.length] - b.exps[b.length])
      a = a*l - s
      k -= 1
   end
   return a*l^k
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(a::SparsePoly{T}, b::S) where {S <: RingElement, T <: RingElement}
   if a.length == 0
      return base_ring(a)()
   end
   r = a.coeffs[a.length]
   for i = 1:a.length - 1
      r *= b^(reinterpret(Int, a.exps[a.length - i + 1] - a.exps[a.length - i]))
      r += a.coeffs[a.length - i]
   end
   if a.exps[1] != 0
      r *= b^(reinterpret(Int, a.exps[1]))
   end
   return r
end

function evaluate(a::SparsePoly{T}, b::Rational{S}) where {S <: Integer, T <: RingElement}
   if a.length == 0
      return base_ring(a)()
   end
   r = a.coeffs[a.length]
   for i = 1:a.length - 1
      r *= b^(reinterpret(Int, a.exps[a.length - i + 1] - a.exps[a.length - i]))
      r += a.coeffs[a.length - i]
   end
   if a.exps[1] != 0
      r *= b^(reinterpret(Int, a.exps[1]))
   end
   return r
end

function evaluate(a::SparsePoly{T}, b::Integer) where {T <: RingElement}
   if a.length == 0
      return base_ring(a)()
   end
   R = base_ring(a)
   r = a.coeffs[a.length]
   for i = 1:a.length - 1
      r *= R(b)^(reinterpret(Int, a.exps[a.length - i + 1] - a.exps[a.length - i]))
      r += a.coeffs[a.length - i]
   end
   if a.exps[1] != 0
      r *= R(b)^(reinterpret(Int, a.exps[1]))
   end
   return r
end

###############################################################################
#
#   GCD, content and primitive part
#
###############################################################################

function gcd(a::SparsePoly{T}, b::SparsePoly{T}, ignore_content::Bool = false) where {T <: RingElement}
   # ensure degree in main variable of a is at least that of b
   if b.exps[b.length] > a.exps[a.length]
      (a, b) = (b, a)
   end
   if iszero(b)
      return deepcopy(a)
   end
   if isone(b)
      return deepcopy(b)
   end
   # compute gcd of contents and divide content out
   if !ignore_content
      c1 = content(a)
      c2 = content(b)
      c = gcd(c1, c2)
      a = divexact(a, c1)
      b = divexact(b, c2)
   end
   # check if we are in the univariate case
   constant_coeffs = true
   for i = 1:a.length
      if !is_constant(a.coeffs[i])
         constant_coeffs = false
         break
      end
   end
   if constant_coeffs
      for i = 1:b.length
         if !is_constant(b.coeffs[i])
            constant_coeffs = false
            break
         end
      end
   end
   # if we are in univariate case, convert to dense, take gcd, convert back
   if constant_coeffs
      # convert polys to univariate dense
      R = AbstractAlgebra.PolyRing(base_ring(base_ring(a)))
      f = R()
      g = R()
      fit!(f, reinterpret(Int, a.exps[a.length] + 1))
      fit!(g, reinterpret(Int, b.exps[b.length] + 1))
      for i = 1:a.length
         f = setcoeff!(f, reinterpret(Int, a.exps[i]), a.coeffs[i].coeffs[1])
      end
      for i = 1:b.length
         g = setcoeff!(g, reinterpret(Int, b.exps[i]), b.coeffs[i].coeffs[1])
      end
      # take gcd of univariate dense polys
      h = gcd(f, g)
      # convert back to sparse polys
      nonzero = 0
      Ac = Vector{T}(undef, 0)
      Ae = zeros(UInt, 0)
      for i = 1:length(h)
         ci = coeff(h, i - 1)
         if ci != 0
            push!(Ac, base_ring(a)(ci))
            push!(Ae, UInt(i - 1))
         end
      end
      r = parent(a)(Ac, Ae)
      if !ignore_content
         return c*r
      else
         return r
      end
   end
   # compute likely degree of gcd
   deg = 0
   # is the lead/trail term a monomial
   lead_monomial = leading_coefficient(a).length == 1 ||
                   leading_coefficient(b).length == 1
   trail_monomial = trailing_coefficient(a).length == 1 ||
                    trailing_coefficient(b).length == 1
   lead_a = leading_coefficient(a)
   lead_b = leading_coefficient(b)
   # psr algorithm
   g = one(base_ring(a))
   h = one(base_ring(a))
   while true
      adeg = reinterpret(Int,  a.exps[a.length])
      bdeg = reinterpret(Int,  b.exps[b.length])
      # optimisation: if degree of b is equal to likely degree of gcd, try
      # exact division, learned from Bernard Parisse
      if bdeg == deg
         flag, q = divides(a, b)
         if flag
            break
         end
      end
      d = reinterpret(Int, adeg - bdeg)
      r = pseudorem(a, b)
      # zero remainder
      if iszero(r)
         break
      end
      # constant remainder
      if r.length == 1 && r.exps[1] == 0
         b = one(parent(a))
         break
      end
      (a, b) = (b, divexact(r, g*h^d))
      g = leading_coefficient(a)
      if d > 1
         h = divexact(g^d, h^(d - 1))
      else
         h = h^(1 - d)*g^d
      end
   end
   # sometimes don't care about content, e.g. when computing likely gcd degree
   if !ignore_content
      # remove content from b as cheaply as possible, as per Bernard Parisse
      if leading_coefficient(b).length != 1 &&
         trailing_coefficient(b).length != 1
         if lead_monomial # lead term monomial, so content contains rest
            d = divexact(leading_coefficient(b),
                         term_content(leading_coefficient(b)))
            b = divexact(b, d)
         elseif trail_monomial # trail term is monomial, so ditto
            d = divexact(trailing_coefficient(b),
                         term_content(trailing_coefficient(b)))
            b = divexact(b, d)
         else
            glead = gcd(lead_a, lead_b)
            if glead.length == 1 # gcd of lead coeffs monomial
               d = divexact(leading_coefficient(b),
                            term_content(leading_coefficient(b)))
               b = divexact(b, d)
            else # last ditched attempt to find easy content
               h = gcd(leading_coefficient(b), glead)
               h = divexact(h, term_content(h))
               flag, q = divides(b, h)
               if flag
                  b = q
               end
            end
         end
      end
      # remove any monomial content
      b1 = term_content(b.coeffs[1])
      for i = 2:b.length
         b1 = gcd(b1, term_content(b.coeffs[i]))
         if isone(b1)
            break
         end
      end
      b = divexact(b, b1)
      # remove any stubborn content and put back actual content
      return c*primpart(b)
   else
      return b
   end
end

function content(a::SparsePoly{T}) where {T <: RingElem}
   for i = 1:length(a)
      if a.coeffs[i].length == 1
         z = term_content(a.coeffs[1])
         for j = 2:length(a)
            if isone(z)
               return z
            end
            z = gcd(z, term_content(a.coeffs[j]))
         end
         return z
      end
   end
   z = base_ring(a)()
   for i = 1:length(a)
      z = gcd(coeff(a, i - 1), z)
   end
   return z
end

function primpart(a::SparsePoly{T}) where {T <: RingElement}
   d = content(a)
   return divexact(a, d)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function fit!(a::SparsePoly{T}, n::Int) where {T <: RingElement}
   if length(a.coeffs) < n
      resize!(a.coeffs, n)
      resize!(a.exps, n)
   end
   return nothing
end

function addmul!(a::SparsePoly{T}, b::SparsePoly{T}, c::SparsePoly{T}, d::SparsePoly{T}) where {T <: RingElement}
   t = b*c
   t += a
   a.coeffs = t.coeffs
   a.exps = t.exps
   a.length = t.length
   return a
end

function mul!(a::SparsePoly{T}, b::SparsePoly{T}, c::SparsePoly{T}) where {T <: RingElement}
   t = b*c
   a.coeffs = t.coeffs
   a.exps = t.exps
   a.length = t.length
   return a
end

function addeq!(a::SparsePoly{T}, b::SparsePoly{T}) where {T <: RingElement}
   t = a + b
   a.coeffs = t.coeffs
   a.exps = t.exps
   a.length = t.length
   return a
end

function add!(a::SparsePoly{T}, b::SparsePoly{T}, c::SparsePoly{T}) where {T <: RingElement}
   t = b + c
   a.coeffs = t.coeffs
   a.exps = t.exps
   a.length = t.length
   return a
end

function zero!(a::SparsePoly{T}) where {T <: RingElement}
   a.length = 0
   return a
end

function shift_left!(a::SparsePoly{T}, n::UInt) where {T <: RingElement}
   for i = 1:length(a)
      a.exps[i] += n
   end
   return a
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{SparsePoly{T}}, ::Type{SparsePoly{T}}) where T <: RingElement = SparsePoly{T}

function promote_rule(::Type{SparsePoly{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? SparsePoly{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::SparsePolyRing{T} where {T <: RingElement})(b::RingElement)
   return a([base_ring(a)(b)], [UInt(0)])
end

function (a::SparsePolyRing{T})() where {T <: RingElement}
   z = SparsePoly{T}()
   z.parent = a
   return z
end

function (a::SparsePolyRing{T})(b::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
   z = SparsePoly{T}(base_ring(a)(b))
   z.parent = a
   return z
end

function (a::SparsePolyRing{T})(b::T) where {T <: RingElement}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = SparsePoly{T}(b)
   z.parent = a
   return z
end

function (a::SparsePolyRing{T})(b::SparsePoly{T}) where {T <: RingElement}
   parent(b) != a && error("Unable to coerce polynomial")
   return b
end

function (a::SparsePolyRing{T})(b::Vector{T}, m::Vector{UInt}) where {T <: RingElement}
   if length(b) > 0
      parent(b[1]) != base_ring(a) && error("Unable to coerce to polynomial")
   end
   z = SparsePoly{T}(b, m)
   z.parent = a
   return z
end

###############################################################################
#
#   SparsePolynomialRing constructor
#
###############################################################################

function SparsePolynomialRing(R::AbstractAlgebra.Ring, s::Symbol; cached::Bool = true)
   T = elem_type(R)

   par = SparsePolyRing{T}(R, s, cached)

   return par, gen(par)
end
