###############################################################################
#
#   residue_field.jl : residue fields (modulo a principal ideal)
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring_type(::Type{ResidueField{T}}) where T <: RingElement = parent_type(T)

base_ring(S::ResidueField{T}) where {T <: RingElement} = S.base_ring::parent_type(T)

parent(a::ResFieldElem) = a.parent

is_domain_type(a::Type{T}) where T <: ResFieldElem = true

function is_exact_type(a::Type{T}) where {S <: RingElement, T <: ResFieldElem{S}}
   return is_exact_type(S)
end

function check_parent(a::ResFieldElem, b::ResFieldElem, throw::Bool = true)
   Ra = parent(a)
   Rb = parent(b)
   if Ra != Rb
      fl = typeof(Ra) == typeof(Rb) && modulus(Ra) == modulus(Rb)
      !fl && throw && error("Incompatible moduli in residue operation")
      #CF: maybe extend to divisibility?
      return fl
   end
   return true
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::ResFieldElem, h::UInt)
   b = 0x539c1c8715c1adc2%UInt
   return xor(b, xor(hash(data(a), h), h))
end

@doc raw"""
    modulus(S::ResidueField)

Return the modulus $a$ of the given residue ring $S = R/(a)$.
"""
function modulus(S::ResidueField)
   return S.modulus
end

@doc raw"""
    modulus(r::ResFieldElem)

Return the modulus $a$ of the residue ring $S = R/(a)$ that the supplied
residue $r$ belongs to.
"""
function modulus(r::ResFieldElem)
   return modulus(parent(r))
end

characteristic(r::ResidueField{T}) where T <: Integer = modulus(r)
is_known(::typeof(characteristic), R::ResidueField{T}) where T <: Integer = true

data(a::ResFieldElem) = a.data

lift(a::ResFieldElem) = data(a)

lift(a::ResFieldElem{Int}) = BigInt(data(a))

zero(R::ResidueField) = R(0)

one(R::ResidueField) = R(1)

iszero(a::ResFieldElem) = iszero(data(a))

isone(a::ResFieldElem) = isone(data(a))

function is_unit(a::ResFieldElem)
   g = gcd(data(a), modulus(a))
   return isone(g)
end

deepcopy_internal(a::ResFieldElem, dict::IdDict) = parent(a)(deepcopy_internal(data(a), dict))

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function expressify(@nospecialize(a::ResFieldElem); context = nothing)
   return expressify(data(a), context = context)
end

@enable_all_show_via_expressify ResFieldElem

function show(io::IO, a::ResidueField)
   @show_name(io, a)
   @show_special(io, a)
   if is_terse(io)
     print(io, "Residue field")
   else
     io = pretty(io)
     print(io, "Residue field of ",)
     print(terse(io), Lowercase(), base_ring(a))
     print(io, " modulo ", modulus(a))
   end
end

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(a::ResFieldElem)
   parent(a)(-data(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) + data(b))
end

function -(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) - data(b))
end

function *(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(data(a) * data(b))
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

*(a::ResFieldElem, b::JuliaRingElement) = parent(a)(data(a) * b)

*(a::ResFieldElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) * b)

*(a::JuliaRingElement, b::ResFieldElem) = parent(b)(a * data(b))

*(a::T, b::ResFieldElem{T}) where {T <: RingElem} = parent(b)(a * data(b))

+(a::ResFieldElem, b::JuliaRingElement) = parent(a)(data(a) + b)

+(a::ResFieldElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) + b)

+(a::JuliaRingElement, b::ResFieldElem) = parent(b)(a + data(b))

+(a::T, b::ResFieldElem{T}) where {T <: RingElem} = parent(b)(a + data(b))

-(a::ResFieldElem, b::JuliaRingElement) = parent(a)(data(a) - b)

-(a::ResFieldElem{T}, b::T) where {T <: RingElem} = parent(a)(data(a) - b)

-(a::JuliaRingElement, b::ResFieldElem) = parent(b)(a - data(b))

-(a::T, b::ResFieldElem{T}) where {T <: RingElem} = parent(b)(a - data(b))

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::ResFieldElem, b::Integer)
   parent(a)(powermod(data(a), b, modulus(a)))
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc raw"""
    ==(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}

Return `true` if $a == b$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   fl = check_parent(a, b, false)
   !fl && return false
   return data(a) == data(b)
end

@doc raw"""
    isequal(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}

Return `true` if $a == b$ exactly, otherwise return `false`. This function is
useful in cases where the data of the residues are inexact, e.g. power series
Only if the power series are precisely the same, to the same precision, are
they declared equal by this function.
"""
function isequal(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return isequal(data(a), data(b))
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function ==(a::ResFieldElem, b::JuliaRingElement)
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

function ==(a::JuliaRingElement, b::ResFieldElem)
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

function ==(a::ResFieldElem{T}, b::T) where {T <: RingElem}
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

function ==(a::T, b::ResFieldElem{T}) where {T <: RingElem}
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc raw"""
    inv(a::ResFieldElem)

Return the inverse of the element $a$ in the residue ring. If an impossible
inverse is encountered, an exception is raised.
"""
function Base.inv(a::ResFieldElem)
   g, ainv = gcdinv(data(a), modulus(a))
   if !isone(g)
      error("Impossible inverse in inv")
   end
   return parent(a)(ainv)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::ResFieldElem{T}, b::ResFieldElem{T}; check::Bool=true) where {T <: RingElement}
   check_parent(a, b)
   fl, q = divides(a, b)
   return q
end

function divides(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   iszero(b) && error("Division by zero in divides")
   return true, a*inv(b)
end

###############################################################################
#
#   GCD
#
###############################################################################

@doc raw"""
    gcd(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}

Return a greatest common divisor of $a$ and $b$ if one exists. This is done
by taking the greatest common divisor of the data associated with the
supplied residues and taking its greatest common divisor with the modulus.
"""
function gcd(a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   check_parent(a, b)
   return parent(a)(gcd(gcd(data(a), modulus(a)), data(b)))
end

###############################################################################
#
#   Square root
#
###############################################################################

function is_square(a::ResFieldElem{T}) where T <: Integer
   if iszero(a)
      return true
   end
   p = modulus(a)
   pm1div2 = div(p - 1, 2)
   return isone(a^pm1div2)
end

function Base.sqrt(a::ResFieldElem{T}; check::Bool=true) where T <: Integer
   U = parent(a)
   p = modulus(a)
   if p == 2 # special case, cannot find a quadratic nonresidue mod 2
      return deepcopy(a)
   end
   # Compute Q, S such that p - 1 = Q*2^S
   Q = p - 1
   S = 0 # power of 2 dividing p - 1
   while iseven(Q)
      Q >>= 1
      S += 1
   end
   # find a quadratic nonresidue z mod p
   z = U(rand(1:p - 1))
   while is_square(z)
      z = U(rand(1:p - 1))
   end
   # set up
   M = S
   c = z^Q
   t = a^Q
   R = a^div(Q + 1, 2)
   # main loop
   while true
      if iszero(t)
         return zero(U)
      end
      if isone(t)
         return R
      end
      u = t
      i = 0
      while i < M
         if isone(u)
            break
         end
         u = u^2
         i += 1
      end
      check && i == M && error("Not a square in sqrt")
      b = c
      for j = 1:M - i - 1
         b = b^2
      end
      M = i
      c = b^2
      t *= c
      R *= b
   end
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::ResFieldElem{T}) where {T <: RingElement}
   a.data = zero!(a.data)
   return a
end

function mul!(c::ResFieldElem{T}, a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   c.data = mod(data(a)*data(b), modulus(a))
   return c
end

function add!(c::ResFieldElem{T}, a::ResFieldElem{T}, b::ResFieldElem{T}) where {T <: RingElement}
   c.data = mod(data(a) + data(b), modulus(a))
   return c
end

###############################################################################
#
#   Random functions
#
###############################################################################

RandomExtensions.maketype(R::ResidueField, _) = elem_type(R)

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:ResFieldElem{T},
                                         <:ResidueField{T}}}
              ) where {T}
   S, v = sp[][1:end]
   S(rand(rng, v))
end

function RandomExtensions.make(S::ResidueField, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1])
   else
      Make(S, make(base_ring(S), vs...))
   end
end

rand(rng::AbstractRNG, S::ResidueField, v...) = rand(rng, make(S, v...))

rand(S::ResidueField, v...) = rand(Random.default_rng(), S, v...)

###############################################################################
#
#   residue_field constructor
#
###############################################################################

@doc raw"""
    residue_field(R::Ring, a::RingElement; cached::Bool = true)

Create the residue ring $R/(a)$ where $a$ is an element of the ring $R$. We
require $a \neq 0$. If `cached == true` (the default) then the resulting
residue ring parent object is cached and returned for any subsequent calls
to the constructor with the same base ring $R$ and element $a$.
"""
function residue_field(R::Ring, a::RingElement; cached::Bool = true)
   @req !is_trivial(R) "Base ring must not be the zero ring."
   iszero(a) && throw(DivideError())
   @req !is_unit(a) "Cannot create a field with one element"
   T = elem_type(R)
   S = EuclideanRingResidueField{T}(R(a), cached)
   return S, Generic.EuclideanRingResidueMap(R, S)
end

@doc raw"""
    quo(::Type{Field}, R::Ring, a::RingElement; cached::Bool = true)

Returns `S, f` where `S = residue_field(R, a)` and `f` is the 
projection map from `R` to `S`. This map is supplied as a map with section
where the section is the lift of an element of the residue field back
to the ring `R`.
"""
function quo(::Type{Field}, R::Ring, a::RingElement; cached::Bool = true)
   S, f = residue_field(R, a; cached = cached)
   return S, f
end
