###############################################################################
#
#   Residue.jl : generic residue rings (modulo a principal ideal)
#
###############################################################################

export ResidueRing, Residue, inv, modulus, data

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type{T}(::Type{Residue{T}}) = ResidueRing{T}

elem_type{T <: RingElem}(::ResidueRing{T}) = Residue{T}

base_ring{T}(a::ResidueRing{T}) = a.base_ring::parent_type(T)

base_ring(a::ResidueElem) = base_ring(parent(a))

parent(a::ResidueElem) = a.parent

function check_parent_type{T <: RingElem}(a::ResidueRing{T}, b::ResidueRing{T})
   # exists only to check types of parents agree
end
   
function check_parent(a::ResidueElem, b::ResidueElem)
   if parent(a) != parent(b)
      check_parent_type(parent(a), parent(b))
      a.hash != b.hash && error("Incompatible moduli in residue operation")
   end
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::ResidueElem, h::UInt)
   b = 0x539c1c8715c1adc2%UInt
   return b $ hash(data(a), h) $ h
end

function modulus(R::ResidueRing)
   return R.modulus
end

function modulus(a::ResidueElem)
   return modulus(parent(a))
end

data(a::Residue) = a.data

zero(R::ResidueRing) = R(0)

one(R::ResidueRing) = R(1)

iszero(a::ResidueElem) = iszero(data(a))

isone(a::ResidueElem) = isone(data(a))

function isunit(a::ResidueElem)
   g, ainv = gcdinv(data(a), modulus(a))
   return g == 1
end

deepcopy(a::ResidueElem) = parent(a)(deepcopy(data(a)))

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::Residue) = a

###############################################################################
#
#   AbstractString{} I/O
#
###############################################################################

function show(io::IO, x::ResidueElem)
   print(io, data(x))
end

function show(io::IO, a::ResidueRing)
   print(io, "Residue ring of ", base_ring(a), " modulo ", modulus(a))
end

needs_parentheses(x::ResidueElem) = needs_parentheses(data(x))

is_negative(x::ResidueElem) = is_negative(data(x))

show_minus_one{T <: RingElem}(::Type{Residue{T}}) = true

###############################################################################
#
#   Unary operations
#
###############################################################################

function -{T <: RingElem}(a::ResidueElem{T})
   parent(a)(-data(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +{T <: RingElem}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(data(a) + data(b))
end

function -{T <: RingElem}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(data(a) - data(b))
end

function *{T <: RingElem}(a::Residue{T}, b::Residue{T})
   check_parent(a, b)
   return parent(a)(data(a) * data(b))
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

*{T <: RingElem}(a::Residue{T}, b::Integer) = parent(a)(data(a) * b)

*{T <: RingElem}(a::Integer, b::Residue{T}) = parent(b)(a * data(b))

+{T <: RingElem}(a::Residue{T}, b::Integer) = parent(a)(data(a) + b)

+{T <: RingElem}(a::Integer, b::Residue{T}) = parent(b)(a + data(b))

-{T <: RingElem}(a::Residue{T}, b::Integer) = parent(a)(data(a) - b)

-{T <: RingElem}(a::Integer, b::Residue{T}) = parent(b)(a - data(b))

###############################################################################
#
#   Powering
#
###############################################################################

function ^{T <: RingElem}(a::ResidueElem{T}, b::Int)
   parent(a)(powmod(data(a), b, modulus(a)))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function =={T <: RingElem}(a::ResidueElem{T}, b::ResidueElem{T})
   check_parent(a, b)
   return data(a) == data(b)
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function =={T <: RingElem}(a::ResidueElem{T}, b::Integer)
   z = base_ring(a)(b)
   return data(a) == mod(z, modulus(a))
end

function =={T <: RingElem}(a::Integer, b::ResidueElem{T})
   z = base_ring(b)(a)
   return data(b) == mod(z, modulus(b))
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv{T <: RingElem}(a::ResidueElem{T})
   g, ainv = gcdinv(data(a), modulus(a))
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

function divexact{T <: RingElem}(a::ResidueElem{T}, b::ResidueElem{T})
   check_parent(a, b)
   g, binv = gcdinv(data(b), modulus(b))
   if g != 1
      error("Impossible inverse in divexact")
   end
   return parent(a)(data(a) * binv)
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd{T <: RingElem}(a::ResidueElem{T}, b::ResidueElem{T})
   check_parent(a, b)
   return parent(a)(gcd(gcd(data(a), modulus(a)), data(b)))
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function mul!{T <: RingElem}(c::Residue{T}, a::Residue{T}, b::Residue{T})
   c.data = mod(data(a)*data(b), modulus(a))
end

function addeq!{T <: RingElem}(c::Residue{T}, a::Residue{T})
   c.data = mod(c.data + data(a), modulus(a))
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

function ResidueRing{T <: RingElem}(R::Ring, el::T; cached=true)
   parent(el) != R && error("Modulus is not an element of the specified ring")
   el == 0 && throw(DivideError())
   
   R2 = R
   par_type = Residue{T}
   sig_table = [x.sig for x in methods(promote_rule)]
   filter!(x -> x <: Tuple{Type{par_type}, Any}, sig_table)
   while base_ring(R2) != Union{}
      R2 = base_ring(R2)
      T2 = elem_type(R2)
      if !in(Tuple{Type{par_type}, Type{T2}}, sig_table)
         eval(:(Base.promote_rule(::Type{$par_type}, ::Type{$T2}) = $par_type))
      end
   end
   return ResidueRing{T}(el, cached)
end

function ResidueRing(R::FlintIntegerRing, el::Integer; cached=true)
   el == 0 && throw(DivideError())
   
   return ResidueRing{fmpz}(R(el), cached)
end
