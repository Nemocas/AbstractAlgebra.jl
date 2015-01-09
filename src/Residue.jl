###########################################################################################
#
#   Residue.jl : generic residue rings (modulo a principal ideal)
#
###########################################################################################

export ResRing, Residue, ResidueRing, inv

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

ModulusDict = Dict{Union(RingElem, BigInt), Ring}()

type ResRing{T <: Union(RingElem, BigInt)} <: Ring
   base_ring::Ring
   modulus::T

   function ResRing(R::Ring, modulus::T)
      return try
         ModulusDict[modulus]
      catch
         ModulusDict[modulus] = new(R, modulus)
      end
   end
end

type Residue{T <: Union(RingElem, BigInt)} <: RingElem
   data::T
   parent::ResRing

   Residue(a::T) = new(a)
end

elem_type{T <: RingElem}(::ResRing{T}) = Residue{T}

base_ring(a::ResRing) = a.base_ring

base_ring(a::Residue) = base_ring(parent(a))

parent(a::Residue) = a.parent

function check_parent(a::Residue, b::Residue)
   parent(a) != parent(b) && error("Incompatible moduli in residue operation")
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function hash(a::Residue)
   h = 0x539c1c8715c1adc2
   return h $ hash(a.data)
end

function modulus(R::ResRing)
   return R.modulus
end

function modulus(a::Residue)
   return modulus(parent(a))
end

zero(R::ResRing) = R(0)

one(R::ResRing) = R(1)

iszero(a::Residue) = iszero(a.data)

isone(a::Residue) = isone(a.data)

function isunit(a::Residue)
   g, ainv = gcdinv(a.data, modulus(a))
   return g == 1
end

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit(a::Residue) = a

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, x::Residue)
   print(io, x.data)
end

function show(io::IO, a::ResRing)
   print(io, "Residue ring of ", base_ring(a), " modulo ", modulus(a))
end

needs_parentheses(x::Residue) = needs_parentheses(x.data)

is_negative(x::Residue) = is_negative(x.data)

show_minus_one{T <: Union(RingElem, BigInt)}(::Type{Residue{T}}) = true

###########################################################################################
#
#   Binary operators
#
###########################################################################################

function +{T <: Union(RingElem, BigInt)}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(a.data + b.data)
end

function +{T <: Union(RingElem, BigInt)}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(a.data + b.data)
end

function +{T <: Union(RingElem, BigInt)}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(a.data + b.data)
end

###########################################################################################
#
#   Ad hoc binary operations
#
###########################################################################################

*{T <: Union(RingElem, BigInt)}(a::Residue{T}, b::Integer) = parent(a)(a.data * b)

*{T <: Union(RingElem, BigInt)}(a::Integer, b::Residue{T}) = parent(b)(a * b.data)

+{T <: Union(RingElem, BigInt)}(a::Residue{T}, b::Integer) = parent(a)(a.data + b)

+{T <: Union(RingElem, BigInt)}(a::Integer, b::Residue{T}) = parent(b)(a + b.data)

-{T <: Union(RingElem, BigInt)}(a::Residue{T}, b::Integer) = parent(a)(a.data - b)

-{T <: Union(RingElem, BigInt)}(a::Integer, b::Residue{T}) = parent(b)(a - b.data)

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{T <: Union(RingElem, BigInt)}(a::Residue{T}, b::Int)
   parent(a)(powmod(a.data, b, modulus(a)))
end

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{T <: RingElem}(a::Residue{T})
   g, ainv = gcdinv(a.data, modulus(a))
   if g != 1
      error("Impossible inverse in inv")
   end
   return parent(a)(ainv)
end

function inv(a::Residue{BigInt})
   ainv = invmod(a.data, modulus(a))
   return parent(a)(ainv)
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{T <: Union(RingElem, BigInt)}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   g, binv = gcdinv(b.data, modulus(b))
   if g != 1
      error("Impossible inverse in divexact")
   end
   return parent(a)(a.data * binv)
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function gcd{T <: Union(RingElem, BigInt)}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(gcd(gcd(a.data, modulus(a)), b.data))
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function mul!{T <: Union(RingElem, BigInt)}(c::Residue{T}, a::Residue{T}, b::Residue{T})
   c.data = mod(a.data*b.data, modulus(a))
end

function addeq!{T <: Union(RingElem, BigInt)}(c::Residue{T}, a::Residue{T})
   c.data = mod(c.data + a.data, modulus(a))
end

###########################################################################################
#
#   Promotion rules
#
###########################################################################################

Base.promote_rule{T <: RingElem}(::Type{Residue{T}}, ::Type{T}) = Residue{T}

Base.promote_rule{T <: RingElem, U <: Integer}(::Type{Residue{T}}, ::Type{U}) = Residue{T}

Base.promote_rule{T <: Integer}(::Type{Residue{BigInt}}, ::Type{T}) = Residue{BigInt}

###########################################################################################
#
#   Parent object call overloading
#
###########################################################################################

function Base.call{T <: RingElem}(a::ResRing{T})
   z = Residue{T}(zero(base_ring(a)))
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::ResRing{T}, b::Integer)
   z = Residue{T}(mod(base_ring(a)(b), modulus(a)))
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::ResRing{T}, b::T)
   check_parent(base_ring(a), b)
   z = Residue{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::ResRing{T}, b::ResRing{T})
   check_parent(a, b)
   return b
end

function Base.call(a::ResRing{BigInt})
   z = Residue{BigInt}(ZZ())
   z.parent = a
   return z
end

function Base.call(a::ResRing{BigInt}, b::Integer)
   z = Residue{BigInt}(mod(ZZ(b), modulus(a)))
   z.parent = a
   return z
end

Base.call(a::ResRing{BigInt}, b::Residue{BigInt}) = b

###########################################################################################
#
#   ResidueRing constructor
#
###########################################################################################

function ResidueRing{T <: RingElem}(R::Ring, el::T)
   parent(el) != R && error("Modulus is not an element of the specified ring")
   el == 0 && throw(DivideError())
   
   return ResRing{T}(R, el)
end

function ResidueRing(R::IntegerRing, el::Integer)
   el == 0 && throw(DivideError())
   
   return ResRing{BigInt}(R, R(el))
end