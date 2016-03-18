###############################################################################
#
#   fq.jl : Flint finite fields
#
###############################################################################

export FlintFiniteField, characteristic, order, fq, FqFiniteField, frobenius,
       pth_root, trace, norm

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent_type(::Type{fq}) = FqFiniteField

elem_type(::FqFiniteField) = fq

base_ring(a::FqFiniteField) = Union{}

base_ring(a::fq) = Union{}

parent(a::fq) = a.parent

function check_parent(a::fq, b::fq) 
   a.parent != b.parent && error("Operations on distinct finite fields not supported")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::fq, h::UInt)
   b = 0xb310fb6ea97e1f1a%UInt
   for i in 1:degree(parent(a)) + 1
         b $= hash(coeff(a, i), h) $ h
         b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

function coeff(x::fq, n::Int)
   n < 0 && throw(DomainError())
   z = fmpz()
   ccall((:fmpz_poly_get_coeff_fmpz, :libflint), Void, 
               (Ptr{fmpz}, Ptr{fq}, Int), &z, &x, n)
   return z
end

function zero(a::FqFiniteField)
   d = a()
   ccall((:fq_zero, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &d, &a)
   return d
end
   
function one(a::FqFiniteField)
   d = a()
   ccall((:fq_one, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &d, &a)
   return d
end
   
function gen(a::FqFiniteField)
   d = a()
   ccall((:fq_gen, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &d, &a)
   return d
end

iszero(a::fq) = ccall((:fq_is_zero, :libflint), Bool, 
                     (Ptr{fq}, Ptr{FqFiniteField}), &a, &a.parent)

isone(a::fq) = ccall((:fq_is_one, :libflint), Bool,
                    (Ptr{fq}, Ptr{FqFiniteField}), &a, &a.parent)

isgen(a::fq) = a == gen(parent(a))

isunit(a::fq) = ccall((:fq_is_invertible, :libflint), Bool, 
                     (Ptr{fq}, Ptr{FqFiniteField}), &a, &a.parent)

function characteristic(a::FqFiniteField)
   d = fmpz()
   ccall((:__fq_ctx_prime, :libflint), Void, 
         (Ptr{fmpz}, Ptr{FqFiniteField}), &d, &a)
   return d
end
   
function order(a::FqFiniteField)
   d = fmpz()
   ccall((:fq_ctx_order, :libflint), Void, 
         (Ptr{fmpz}, Ptr{FqFiniteField}), &d, &a)
   return d
end
   
function degree(a::FqFiniteField)
   return ccall((:fq_ctx_degree, :libflint), Int, (Ptr{FqFiniteField},), &a)
end

function deepcopy(d::fq)
   z = fq(parent(d), d)
   return z
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::fq) = x

###############################################################################
#
#   AbstractString{} I/O
#
###############################################################################

function show(io::IO, x::fq)
   cstr = ccall((:fq_get_str_pretty, :libflint), Ptr{UInt8}, 
                (Ptr{fq}, Ptr{FqFiniteField}), &x, &x.parent)

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{UInt8},), cstr)
end

function show(io::IO, a::FqFiniteField)
   print(io, "Finite field of degree ", degree(a))
   print(io, " over F_", characteristic(a))
end

needs_parentheses(x::fq) = x.length > 1

is_negative(x::fq) = false

show_minus_one(::Type{fq}) = true

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(x::fq)
   z = parent(x)()
   ccall((:fq_neg, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   return z
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(x::fq, y::fq)
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_add, :libflint), Void, 
        (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &y, &y.parent)
   return z
end

function -(x::fq, y::fq)
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_sub, :libflint), Void, 
        (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &y, &y.parent)
   return z
end

function *(x::fq, y::fq)
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_mul, :libflint), Void, 
        (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &y, &y.parent)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::Int, y::fq)
   z = parent(y)()
   ccall((:fq_mul_si, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Int, Ptr{FqFiniteField}), &z, &y, x, &y.parent)
   return z
end

*(x::Integer, y::fq) = fmpz(x)*y

*(x::fq, y::Integer) = y*x

function *(x::fmpz, y::fq)
   z = parent(y)()
   ccall((:fq_mul_fmpz, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{fmpz}, Ptr{FqFiniteField}), 
                                            &z, &y, &x, &y.parent)
   return z
end

*(x::fq, y::fmpz) = y*x

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::fq, y::Int)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = parent(x)()
   ccall((:fq_pow_ui, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Int, Ptr{FqFiniteField}), &z, &x, y, &x.parent)
   return z
end

function ^(x::fq, y::fmpz)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = parent(x)()
   ccall((:fq_pow, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{fmpz}, Ptr{FqFiniteField}),
                                            &z, &x, &y, &x.parent)
   return z
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::fq, y::fq) 
   check_parent(x, y)
   ccall((:fq_equal, :libflint), Bool, 
         (Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &x, &y, &y.parent)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::fq, y::fq)
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   z = parent(y)()
   ccall((:fq_div, :libflint), Void, 
        (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &y, &y.parent)
   return z
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(x::fq)
   iszero(x) && throw(DivideError())
   z = parent(x)()
   ccall((:fq_inv, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   return z
end

###############################################################################
#
#   Special functions
#
###############################################################################

function pth_root(x::fq)
   z = parent(x)()
   ccall((:fq_pth_root, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   return z
end

function trace(x::fq)
   z = fmpz()
   ccall((:fq_trace, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   return parent(x)(z)
end

function norm(x::fq)
   z = fmpz()
   ccall((:fq_norm, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   return parent(x)(z)
end

function frobenius(x::fq, n = 1)
   z = parent(x)()
   ccall((:fq_frobenius, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Int, Ptr{FqFiniteField}), &z, &x, n, &x.parent)
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function mul!(z::fq, x::fq, y::fq)
   ccall((:fq_mul, :libflint), Void, 
        (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &y, &y.parent)
end

function addeq!(z::fq, x::fq)
   ccall((:fq_add, :libflint), Void, 
        (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &z, &x, &x.parent)
end

###############################################################################
#
#   Promotions
#
###############################################################################

Base.promote_rule{T <: Integer}(::Type{fq}, ::Type{T}) = fq

Base.promote_rule(::Type{fq}, ::Type{fmpz}) = fq

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function Base.call(a::FqFiniteField)
   z = fq(a)
   return z
end

Base.call(a::FqFiniteField, b::Integer) = a(fmpz(b))

function Base.call(a::FqFiniteField, b::Int)
   z = fq(a, b)
   return z
end

function Base.call(a::FqFiniteField, b::fmpz)
   z = fq(a, b)
   return z
end

function Base.call(a::FqFiniteField, b::fq)
   parent(b) != a && error("Coercion between finite fields not implemented")
   return b
end

###############################################################################
#
#   FlintFiniteField constructor
#
###############################################################################

function FlintFiniteField(char::fmpz, deg::Int, s::AbstractString{})
   S = symbol(s)
   parent_obj = FqFiniteField(char, deg, S)

   return parent_obj, gen(parent_obj)
end

function FlintFiniteField(char::Integer, deg::Int, s::AbstractString{})
   return FlintFiniteField(fmpz(char), deg, s)
end

function FlintFiniteField(pol::fmpz_mod_poly, s::AbstractString{})
   S = symbol(s)
   parent_obj = FqFiniteField(pol, S)

   return parent_obj, gen(parent_obj)
end
