###############################################################################
#
#   fq.jl : Flint finite fields
#
###############################################################################

export fq, FqFiniteField, frobenius, pth_root, trace, norm

###############################################################################
#
#   Types and memory management
#
###############################################################################

FqFiniteFieldID = Dict{Tuple{fmpz, Int, Symbol}, Field}()

type FqFiniteField <: Field
   p::Int # fmpz
   sparse_modulus::Int
   a::Ptr{Void}
   j::Ptr{Void}
   len::Int
   mod_coeffs::Ptr{Void}
   mod_alloc::Int
   mod_length::Int
   mod_p::Int # fmpz
   inv_coeffs::Ptr{Void}
   inv_alloc::Int
   inv_length::Int
   inv_p::Int # fmpz
   var::Ptr{Void}
   S::Symbol

   function FqFiniteField(char::fmpz, deg::Int, s::Symbol)
      try
         return FqFiniteFieldID[char, deg, s]
      catch
         d = FqFiniteFieldID[char, deg, s] = new()
         finalizer(d, _FqFiniteField_clear_fn)
         ccall((:fq_ctx_init, :libflint), Void, 
               (Ptr{FqFiniteField}, Ptr{fmpz}, Int, Ptr{Uint8}), 
			            &d, &char, deg, bytestring(string(s)))
         d.S = s
         return d
      end
   end
end

function _FqFiniteField_clear_fn(a :: FqFiniteField)
   ccall((:fq_ctx_clear, :libflint), Void, (Ptr{FqFiniteField},), &a)
end

type fq <: FiniteFieldElem
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
   parent::FqFiniteField

   function fq(ctx::FqFiniteField)
      d = new()
      ccall((:fq_init2, :libflint), Void, 
            (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      d.parent = ctx
      return d
   end

   function fq(ctx::FqFiniteField, x::Int)
      d = new()
      ccall((:fq_init2, :libflint), Void, 
            (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      ccall((:fq_set_si, :libflint), Void, 
                (Ptr{fq}, Int, Ptr{FqFiniteField}), &d, x, &ctx)
      d.parent = ctx
      return d
   end

   function fq(ctx::FqFiniteField, x::fmpz)
      d = new()
      ccall((:fq_init2, :libflint), Void, 
            (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      ccall((:fq_set_fmpz, :libflint), Void, 
            (Ptr{fq}, Ptr{fmpz}, Ptr{FqFiniteField}), &d, &x, &ctx)
      d.parent = ctx
      return d
   end

   function fq(ctx::FqFiniteField, x::fq)
      d = new()
      ccall((:fq_init2, :libflint), Void, 
            (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      ccall((:fq_set, :libflint), Void, 
            (Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &d, &x, &ctx)
      d.parent = ctx
      return d
   end
end

function _fq_clear_fn(a::fq)
   ccall((:fq_clear, :libflint), Void, 
         (Ptr{fq}, Ptr{FqFiniteField}), &a, &a.parent)
end

elem_type(::FqFiniteField) = fq

base_ring(a::FqFiniteField) = None

base_ring(a::fq) = None

parent(a::fq) = a.parent

function check_parent(a::fq, b::fq) 
   a.parent != b.parent && error("Operations on distinct finite fields not supported")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function hash(a::fq)
   h = 0xb310fb6ea97e1f1a
   for i in 1:degree(parent(a)) + 1
         h $= hash(coeff(a, i))
         h = (h << 1) | (h >> (sizeof(Int)*8 - 1))
   end
   return h
end

function coeff(x::fq, n::Int)
   n < 0 && throw(DomainError())
   z = ZZ()
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

isunit(a::fq) = ccall((:fq_is_invertible, :libflint), Bool, 
                     (Ptr{fq}, Ptr{FqFiniteField}), &a, &a.parent)

function characteristic(a::FqFiniteField)
   d = ZZ()
   ccall((:__fq_ctx_prime, :libflint), Void, 
         (Ptr{fmpz}, Ptr{FqFiniteField}), &d, &a)
   return d
end
   
function order(a::FqFiniteField)
   d = ZZ()
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
#   String I/O
#
###############################################################################

function show(io::IO, x::fq)
   cstr = ccall((:fq_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{fq}, Ptr{FqFiniteField}), &x, &x.parent)

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
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

*(x::Integer, y::fq) = ZZ(x)*y

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
   z = ZZ()
   ccall((:fq_trace, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   return parent(x)(z)
end

function norm(x::fq)
   z = ZZ()
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

Base.call(a::FqFiniteField, b::Integer) = a(ZZ(b))

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
#   FiniteField constructor
#
###############################################################################

function FiniteField(char::fmpz, deg::Int, s::String)
   S = symbol(s)
   parent_obj = FqFiniteField(char, deg, S)

   return parent_obj, gen(parent_obj) 
end

function FiniteField(char::Integer, deg::Int, s::String)
   return FiniteField(ZZ(char), deg, s)
end
