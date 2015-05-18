###############################################################################
#
#   Residue.jl : generic residue rings (modulo a principal ideal)
#
###############################################################################

export ResidueRing, Residue, inv, modulus

###############################################################################
#
#   Data types and memory management
#
###############################################################################

ModulusDict = Dict{Tuple{Ring, RingElem}, Ring}()

type ResidueRing{T <: RingElem} <: Ring
   base_ring::Ring
   modulus::T

   function ResidueRing(modulus::T)
      return try
         ModulusDict[parent(modulus), modulus]
      catch
         ModulusDict[parent(modulus), modulus] = new(parent(modulus), modulus)
      end
   end
end

type Residue{T <: RingElem} <: RingElem
   data::T
   parent::ResidueRing

   Residue(a::T) = new(a)
end

elem_type{T <: RingElem}(::ResidueRing{T}) = Residue{T}

base_ring(a::ResidueRing) = a.base_ring

base_ring(a::Residue) = base_ring(parent(a))

parent(a::Residue) = a.parent

function check_parent(a::Residue, b::Residue)
   parent(a) != parent(b) && error("Incompatible moduli in residue operation")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function hash(a::Residue)
   h = 0x539c1c8715c1adc2
   return h $ hash(a.data)
end

function modulus(R::ResidueRing)
   return R.modulus
end

function modulus(a::Residue)
   return modulus(parent(a))
end

zero(R::ResidueRing) = R(0)

one(R::ResidueRing) = R(1)

iszero(a::Residue) = iszero(a.data)

isone(a::Residue) = isone(a.data)

function isunit(a::Residue)
   g, ainv = gcdinv(a.data, modulus(a))
   return g == 1
end

deepcopy(a::Residue) = parent(a)(deepcopy(a.data))

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::Residue) = a

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, x::Residue)
   print(io, x.data)
end

function show(io::IO, a::ResidueRing)
   print(io, "Residue ring of ", base_ring(a), " modulo ", modulus(a))
end

needs_parentheses(x::Residue) = needs_parentheses(x.data)

is_negative(x::Residue) = is_negative(x.data)

show_minus_one{T <: RingElem}(::Type{Residue{T}}) = true

###############################################################################
#
#   Unary operations
#
###############################################################################

function -{T <: RingElem}(a::Residue{T})
   parent(a)(-a.data)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +{T <: RingElem}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(a.data + b.data)
end

function -{T <: RingElem}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(a.data - b.data)
end

function *{T <: RingElem}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(a.data * b.data)
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

*{T <: RingElem}(a::Residue{T}, b::Integer) = parent(a)(a.data * b)

*{T <: RingElem}(a::Integer, b::Residue{T}) = parent(b)(a * b.data)

+{T <: RingElem}(a::Residue{T}, b::Integer) = parent(a)(a.data + b)

+{T <: RingElem}(a::Integer, b::Residue{T}) = parent(b)(a + b.data)

-{T <: RingElem}(a::Residue{T}, b::Integer) = parent(a)(a.data - b)

-{T <: RingElem}(a::Integer, b::Residue{T}) = parent(b)(a - b.data)

###############################################################################
#
#   Powering
#
###############################################################################

function ^{T <: RingElem}(a::Residue{T}, b::Int)
   parent(a)(powmod(a.data, b, modulus(a)))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function =={T <: RingElem}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return a.data == b.data
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv{T <: RingElem}(a::Residue{T})
   g, ainv = gcdinv(a.data, modulus(a))
   if g != 1
      error("Impossible inverse in inv")
   end
   return parent(a)(ainv)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact{T <: RingElem}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   g, binv = gcdinv(b.data, modulus(b))
   if g != 1
      error("Impossible inverse in divexact")
   end
   return parent(a)(a.data * binv)
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd{T <: RingElem}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(gcd(gcd(a.data, modulus(a)), b.data))
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function mul!{T <: RingElem}(c::Residue{T}, a::Residue{T}, b::Residue{T})
   c.data = mod(a.data*b.data, modulus(a))
end

function addeq!{T <: RingElem}(c::Residue{T}, a::Residue{T})
   c.data = mod(c.data + a.data, modulus(a))
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

Base.promote_rule{T <: RingElem}(::Type{Residue{T}}, ::Type{T}) = Residue{T}

Base.promote_rule{T <: RingElem, U <: Integer}(::Type{Residue{T}}, ::Type{U}) = Residue{T}

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function Base.call{T <: RingElem}(a::ResidueRing{T}, b::RingElem)
   return a(base_ring(a)(b))
end

function Base.call{T <: RingElem}(a::ResidueRing{T})
   z = Residue{T}(zero(base_ring(a)))
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::ResidueRing{T}, b::Integer)
   z = Residue{T}(mod(base_ring(a)(b), modulus(a)))
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::ResidueRing{T}, b::T)
   base_ring(a) != parent(b) && error("Operation on incompatible objects")
   z = Residue{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::ResidueRing{T}, b::Residue{T})
   a != parent(b) && error("Operation on incompatible objects")
   return b
end

###############################################################################
#
#   ResidueRing constructor
#
###############################################################################

function ResidueRing{T <: RingElem}(R::Ring, el::T)
   parent(el) != R && error("Modulus is not an element of the specified ring")
   el == 0 && throw(DivideError())
   
   R2 = R
   parent_type = Residue{T}
   while base_ring(R2) != None
      R2 = base_ring(R2)
      T2 = elem_type(R2)
      eval(:(Base.promote_rule(::Type{$parent_type}, ::Type{$T2}) = $parent_type))
   end

   return ResidueRing{T}(el)
end

function ResidueRing(R::IntegerRing, el::Integer)
   el == 0 && throw(DivideError())
   
   return ResidueRing{fmpz}(R(el))
end