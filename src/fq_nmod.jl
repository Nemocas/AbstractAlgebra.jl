###############################################################################
#
#   fq_nmod.jl : Flint finite fields
#
###############################################################################

export fq_nmod, FqNmodFiniteField

###############################################################################
#
#   Types and memory management
#
###############################################################################

FqNmodFiniteFieldID = Dict{(fmpz, Int, Symbol), Field}()

type FqNmodFiniteField <: Field
   p :: Int # fmpz
   n :: Int
   ninv :: Int
   norm :: Int
   sparse_modulus :: Int
   a :: Ptr{Void}
   j :: Ptr{Void}
   len :: Int
   mod_coeffs :: Ptr{Void}
   mod_alloc :: Int
   mod_length :: Int
   mod_n :: Int
   mod_ninv :: Int
   mod_norm :: Int
   inv_coeffs :: Ptr{Void}
   inv_alloc :: Int
   inv_length :: Int
   inv_n :: Int
   inv_ninv :: Int
   inv_norm :: Int
   var :: Ptr{Void}
   S::Symbol

   function FqNmodFiniteField(char::fmpz, deg::Int, s::Symbol)
      try
         return FqNmodFiniteFieldID[char, deg, s]
      catch
         d = FqNmodFiniteFieldID[char, deg, s] = new()
         finalizer(d, _FqNmodFiniteField_clear_fn)
         ccall((:fq_nmod_ctx_init, :libflint), Void, 
               (Ptr{FqNmodFiniteField}, Ptr{fmpz}, Int, Ptr{Uint8}), 
			    &d, &char, deg, bytestring(string(s)))
         return d
      end
   end
end

function _FqNmodFiniteField_clear_fn(a :: FqNmodFiniteField)
   ccall((:fq_nmod_ctx_clear, :libflint), Void, (Ptr{FqNmodFiniteField},), &a)
end

type fq_nmod <: FiniteFieldElem
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
   n :: Int
   ninv :: Int
   norm :: Int
   parent::FqNmodFiniteField

   function fq_nmod(ctx::FqNmodFiniteField)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void, 
            (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      return d
   end

   function fq_nmod(ctx::FqNmodFiniteField, x::Int)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void,
            (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      ccall((:fq_nmod_set_si, :libflint), Void, 
                (Ptr{fq_nmod}, Int, Ptr{FqNmodFiniteField}), &d, x, &ctx)
      return d
   end

   function fq_nmod(ctx::FqNmodFiniteField, x::fmpz)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void, 
            (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      ccall((:fq_nmod_set_fmpz, :libflint), Void, 
            (Ptr{fq_nmod}, Ptr{fmpz}, Ptr{FqNmodFiniteField}), &d, &x, &ctx)
      return d
   end

      function fq_nmod(ctx::FqNmodFiniteField, x::fq_nmod)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void, 
            (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      ccall((:fq_nmod_set, :libflint), Void, 
            (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &x, &ctx)
      return d
   end
end

function _fq_nmod_clear_fn(a::fq_nmod)
   ccall((:fq_nmod_clear, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &a, &a.parent)
end

elem_type(::FqNmodFiniteField) = fq_nmod

base_ring(a::FqNmodFiniteField) = None

base_ring(a::fq_nmod) = None

parent(a::fq_nmod) = a.parent

function check_parent(a::fq_nmod, b::fq_nmod) 
   a.parent != b.parent && error("Operations on distinct finite fields not supported")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function hash(a::fq_nmod)
   h = 0x78e5f766c8ace18d
   for i in 1:degree(parent(a)) + 1
         h $= hash(coeff(a, i))
         h = (h << 1) | (h >> (sizeof(Int)*8 - 1))
   end
   return h
end

function coeff(x::fq_nmod, n::Int)
   n < 0 && throw(DomainError())
   return ccall((:nmod_poly_get_coeff_ui, :libflint), Uint, 
                (Ptr{fq_nmod}, Int), &x, n) 
end

function zero(a::FqNmodFiniteField)
   d = a()
   ccall((:fq_nmod_zero, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &a)
   return d
end
   
function one(a::FqNmodFiniteField)
   d = a()
   ccall((:fq_nmod_one, :libflint), Void,
         (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &a)
   return d
end
   
function gen(a::FqNmodFiniteField)
   d = a()
   ccall((:fq_nmod_gen, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &a)
   return d
end

iszero(a::fq_nmod) = ccall((:fq_nmod_is_zero, :libflint), Bool, 
                     (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &a, &a.parent)

isone(a::fq_nmod) = ccall((:fq_nmod_is_one, :libflint), Bool,
                    (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &a, &a.parent)

isunit(a::fq_nmod) = ccall((:fq_nmod_is_invertible, :libflint), Bool, 
                     (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &a, &a.parent)

function characteristic(a::FqNmodFiniteField)
   d = ZZ()
   ccall((:__fq_nmod_ctx_prime, :libflint), Void, 
         (Ptr{fmpz}, Ptr{FqNmodFiniteField}), &d, &a)
   return d
end
   
function order(a::FqNmodFiniteField)
   d = ZZ()
   ccall((:fq_nmod_ctx_order, :libflint), Void, 
         (Ptr{fmpz}, Ptr{FqNmodFiniteField}), &d, &a)
   return d
end
   
function degree(a::FqNmodFiniteField)
   return ccall((:fq_nmod_ctx_degree, :libflint), Int,
                (Ptr{FqNmodFiniteField},), &a)
end

function deepcopy(d::fq_nmod)
   z = fq_nmod(parent(d), d)
   z.parent = parent(d)
   return z
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::fq_nmod) = x

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, x::fq_nmod)
   cstr = ccall((:fq_nmod_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &x, &x.parent)

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
end

function show(io::IO, a::FqNmodFiniteField)
   print(io, "Finite field of degree ", degree(a))
   print(io, " over F_", characteristic(a))
end

needs_parentheses(x::fq_nmod) = x.length > 1

is_negative(x::fq_nmod) = false

show_minus_one(::Type{fq_nmod}) = true

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(x::fq_nmod)
   z = parent(x)()
   ccall((:fq_nmod_neg, :libflint), Void, 
       (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &x.parent)
   return z
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(x::fq_nmod, y::fq_nmod)
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_nmod_add, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), 
                                                       &z, &x, &y, &y.parent)
   return z
end

function -(x::fq_nmod, y::fq_nmod)
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_nmod_sub, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), 
                                                       &z, &x, &y, &y.parent)
   return z
end

function *(x::fq_nmod, y::fq_nmod)
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_nmod_mul, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), 
                                                       &z, &x, &y, &y.parent)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::Int, y::fq_nmod)
   z = parent(y)()
   ccall((:fq_nmod_mul_si, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Int, Ptr{FqNmodFiniteField}), 
                                               &z, &y, x, &y.parent)
   return z
end

*(x::Integer, y::fq_nmod) = ZZ(x)*y

*(x::fq_nmod, y::Integer) = y*x

function *(x::fmpz, y::fq_nmod)
   z = parent(y)()
   ccall((:fq_nmod_mul_fmpz, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fmpz}, Ptr{FqNmodFiniteField}), 
                                                    &z, &y, &x, &y.parent)
   return z
end

*(x::fq_nmod, y::fmpz) = y*x

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::fq_nmod, y::Int)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = parent(x)()
   ccall((:fq_nmod_pow_ui, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Int, Ptr{FqNmodFiniteField}), 
                                               &z, &x, y, &x.parent)
   return z
end

function ^(x::fq_nmod, y::fmpz)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = parent(x)()
   ccall((:fq_nmod_pow, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fmpz}, Ptr{FqNmodFiniteField}),
                                                    &z, &x, &y, &x.parent)
   return z
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::fq_nmod, y::fq_nmod) 
   check_parent(x, y)
   ccall((:fq_nmod_equal, :libflint), Bool, 
       (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &x, &y, &y.parent)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::fq_nmod, y::fq_nmod)
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   z = parent(y)()
   ccall((:fq_nmod_div, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), 
                                                       &z, &x, &y, &y.parent)
   return z
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(x::fq_nmod)
   iszero(x) && throw(DivideError())
   z = parent(x)()
   ccall((:fq_nmod_inv, :libflint), Void, 
       (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &x.parent)
   return z
end

###############################################################################
#
#   Special functions
#
###############################################################################

function pth_root(x::fq_nmod)
   z = parent(x)()
   ccall((:fq_nmod_pth_root, :libflint), Void, 
       (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &x.parent)
   return z
end

function trace(x::fq_nmod)
   z = ZZ()
   ccall((:fq_nmod_trace, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &x.parent)
   return parent(x)(z)
end

function norm(x::fq_nmod)
   z = ZZ()
   ccall((:fq_nmod_norm, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &x.parent)
   return parent(x)(z)
end

function frobenius(x::fq_nmod, n = 1)
   z = parent(x)()
   ccall((:fq_nmod_frobenius, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Int, Ptr{FqNmodFiniteField}), 
                                               &z, &x, n, &x.parent)
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function mul!(z::fq_nmod, x::fq_nmod, y::fq_nmod)
   ccall((:fq_nmod_mul, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), 
                                                       &z, &x, &y, &y.parent)
end

function addeq!(z::fq_nmod, x::fq_nmod)
   ccall((:fq_nmod_add, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), 
                                                       &z, &z, &x, &x.parent)
end

###############################################################################
#
#   Promotions
#
###############################################################################

Base.promote_rule{T <: Integer}(::Type{fq_nmod}, ::Type{T}) = fq_nmod

Base.promote_rule(::Type{fq_nmod}, ::Type{fmpz}) = fq_nmod

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function Base.call(a::FqNmodFiniteField)
   z = fq_nmod(a)
   z.parent = a
   return z
end

Base.call(a::FqNmodFiniteField, b::Integer) = a(ZZ(b))

function Base.call(a::FqNmodFiniteField, b::Int)
   z = fq_nmod(a, b)
   z.parent = a
   return z
end

function Base.call(a::FqNmodFiniteField, b::fmpz)
   z = fq_nmod(a, b)
   z.parent = a
   return z
end

function Base.call(a::FqNmodFiniteField, b::fq_nmod)
   parent(b) != a && error("Coercion between finite fields not implemented")
   return b
end

###############################################################################
#
#   FiniteField constructor
#
###############################################################################

function FiniteField(char::Int, deg::Int, s::String)
   S = symbol(s)
   parent_obj = FqNmodFiniteField(ZZ(char), deg, S)

   return parent_obj, gen(parent_obj) 
end
