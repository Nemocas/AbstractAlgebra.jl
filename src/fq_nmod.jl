###########################################################################################
#
#   fq_nmod.jl : Flint finite fields
#
###########################################################################################

export fq_nmod, FqNmodFiniteField

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

FqNmodFiniteFieldID = Dict{(BigInt, Int, Symbol), Field}()

type FqNmodFiniteField{S} <: Field
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

   function FqNmodFiniteField(char::BigInt, deg::Int)
      try
         return FqNmodFiniteFieldID[char, deg, S]
      catch
         d = FqNmodFiniteFieldID[char, deg, S] = new()
         finalizer(d, _FqNmodFiniteField_clear_fn)
         temp = fmpz_readonly(char)
         ccall((:fq_nmod_ctx_init, :libflint), Void, 
               (Ptr{FqNmodFiniteField}, Ptr{fmpz_readonly}, Int, Ptr{Uint8}), 
			    &d, &temp, deg, bytestring(string(S)))
         return d
      end
   end
end

function _FqNmodFiniteField_clear_fn(a :: FqNmodFiniteField)
   ccall((:fq_nmod_ctx_clear, :libflint), Void, (Ptr{FqNmodFiniteField},), &a)
end

type fq_nmod{S} <: FiniteFieldElem
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
   n :: Int
   ninv :: Int
   norm :: Int
   parent::FqNmodFiniteField{S}

   function fq_nmod(ctx::FqNmodFiniteField)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void, (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      return d
   end

   function fq_nmod(ctx::FqNmodFiniteField, x::Int)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void, (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      ccall((:fq_nmod_set_si, :libflint), Void, 
                (Ptr{fq_nmod}, Int, Ptr{FqNmodFiniteField}), &d, x, &ctx)
      return d
   end

   function fq_nmod(ctx::FqNmodFiniteField, x::BigInt)
      d = new()
      ccall((:fq_nmod_init2, :libflint), Void, (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &ctx)
      finalizer(d, _fq_nmod_clear_fn)
      temp = fmpz_readonly(x)
      ccall((:fq_nmod_set_fmpz, :libflint), Void, 
                (Ptr{fq_nmod}, Ptr{fmpz_readonly}, Ptr{FqNmodFiniteField}), &d, &temp, &ctx)
      return d
   end
end

function _fq_nmod_clear_fn(a::fq_nmod)
   ccall((:fq_nmod_clear, :libflint), Void, (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &a, &a.parent)
end

elem_type{S}(::FqNmodFiniteField{S}) = fq_nmod{S}

base_ring(a::FqNmodFiniteField) = None

base_ring(a::fq_nmod) = None

parent(a::fq_nmod) = a.parent

function check_parent(a::fq_nmod, b::fq_nmod) 
   a.parent != b.parent && error("Operations on distinct finite fields not supported")
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function zero{S}(a::FqNmodFiniteField{S})
   d = a()
   ccall((:fq_nmod_zero, :libflint), Void, (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &a)
   return d
end
   
function one{S}(a::FqNmodFiniteField{S})
   d = a()
   ccall((:fq_nmod_one, :libflint), Void, (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &a)
   return d
end
   
function gen{S}(a::FqNmodFiniteField{S})
   d = a()
   ccall((:fq_nmod_gen, :libflint), Void, (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &d, &a)
   return d
end

iszero(a::fq_nmod) = ccall((:fq_nmod_is_zero, :libflint), Bool, 
                     (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &a, &a.parent)

isone(a::fq_nmod) = ccall((:fq_nmod_is_one, :libflint), Bool,
                    (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &a, &a.parent)

isunit(a::fq_nmod) = ccall((:fq_nmod_is_invertible, :libflint), Bool, 
                     (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &a, &a.parent)

function characteristic(a::FqNmodFiniteField)
   d = fmpz()
   ccall((:__fq_nmod_ctx_prime, :libflint), Void, (Ptr{fmpz}, Ptr{FqNmodFiniteField}), &d, &a)
   r = BigInt(d)
   return r
end
   
function order(a::FqNmodFiniteField)
   d = fmpz()
   ccall((:fq_nmod_ctx_order, :libflint), Void, (Ptr{fmpz}, Ptr{FqNmodFiniteField}), &d, &a)
   r = BigInt(d)
   return r
end
   
function degree(a::FqNmodFiniteField)
   return ccall((:fq_nmod_ctx_degree, :libflint), Int, (Ptr{FqNmodFiniteField},), &a)
end
###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit(x::fq_nmod) = x

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, x::fq_nmod)
   cstr = ccall((:fq_nmod_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &x, &x.parent)

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
end

function show(io::IO, a::FqNmodFiniteField)
   print(io, "Finite field of degree ", degree(a), " over F_", characteristic(a))
end

needs_parentheses(x::fq_nmod) = x.length > 1

is_negative(x::fq_nmod) = false

show_minus_one{S}(::Type{fq_nmod{S}}) = true

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -(x::fq_nmod)
   z = parent(x)()
   ccall((:fq_nmod_neg, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &x.parent)
   return z
end

###########################################################################################
#
#   Binary operations
#
###########################################################################################

function +{S}(x::fq_nmod{S}, y::fq_nmod{S})
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_nmod_add, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &y, &y.parent)
   return z
end

function -{S}(x::fq_nmod{S}, y::fq_nmod{S})
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_nmod_sub, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &y, &y.parent)
   return z
end

function *{S}(x::fq_nmod{S}, y::fq_nmod{S})
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_nmod_mul, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &y, &y.parent)
   return z
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *(x::Int, y::fq_nmod)
   z = parent(y)()
   ccall((:fq_nmod_mul_si, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Int, Ptr{FqNmodFiniteField}), &z, &y, x, &y.parent)
   return z
end

*(x::Integer, y::fq_nmod) = BigInt(x)*y

*(x::fq_nmod, y::Integer) = y*x

function *(x::BigInt, y::fq_nmod)
   z = parent(y)()
   temp = fmpz_readonly(x)
   ccall((:fq_nmod_mul_fmpz, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fmpz_readonly}, Ptr{FqNmodFiniteField}), &z, &y, &temp, &y.parent)
   return z
end

*(x::fq_nmod, y::BigInt) = y*x

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^(x::fq_nmod, y::Int)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = parent(x)()
   ccall((:fq_nmod_pow_ui, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Int, Ptr{FqNmodFiniteField}), &z, &x, y, &x.parent)
   return z
end

function ^(x::fq_nmod, y::BigInt)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = parent(x)()
   temp = fmpz_readonly(y)
   ccall((:fq_nmod_pow, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fmpz_readonly}, Ptr{FqNmodFiniteField}), &z, &x, &temp, &x.parent)
   return z
end

###########################################################################################
#
#   Comparison
#
###########################################################################################

function =={S}(x::fq_nmod{S}, y::fq_nmod{S}) 
   check_parent(x, y)
   ccall((:fq_nmod_equal, :libflint), Bool, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &x, &y, &y.parent)
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{S}(x::fq_nmod{S}, y::fq_nmod{S})
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   z = parent(y)()
   ccall((:fq_nmod_div, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &y, &y.parent)
   return z
end

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv(x::fq_nmod)
   iszero(x) && throw(DivideError())
   z = parent(x)()
   ccall((:fq_nmod_inv, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &x.parent)
   return z
end

###########################################################################################
#
#   Special functions
#
###########################################################################################

function pth_root(x::fq_nmod)
   z = parent(x)()
   ccall((:fq_nmod_pth_root, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &x.parent)
   return z
end

function trace(x::fq_nmod)
   z = fmpz()
   init(z)
   ccall((:fq_nmod_trace, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &x.parent)
   r = BigInt(z)
   clear(z)
   return parent(x)(r)
end

function norm(x::fq_nmod)
   z = fmpz()
   init(z)
   ccall((:fq_nmod_norm, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &x.parent)
   r = BigInt(z)
   clear(z)
   return parent(x)(r)
end

function frobenius(x::fq_nmod, n = 1)
   z = parent(x)()
   ccall((:fq_nmod_frobenius, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Int, Ptr{FqNmodFiniteField}), &z, &x, n, &x.parent)
   return z
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function mul!(z::fq_nmod, x::fq_nmod, y::fq_nmod)
   ccall((:fq_nmod_mul, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &x, &y, &y.parent)
end

function addeq!(z::fq_nmod, x::fq_nmod)
   ccall((:fq_nmod_add, :libflint), Void, 
         (Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{fq_nmod}, Ptr{FqNmodFiniteField}), &z, &z, &x, &x.parent)
end

###########################################################################################
#
#   Promotions
#
###########################################################################################

Base.promote_rule{S, T <: Integer}(::Type{fq_nmod{S}}, ::Type{T}) = fq_nmod{S}

Base.promote_rule{S}(::Type{fq_nmod{S}}, ::Type{BigInt}) = fq_nmod{S}

###########################################################################################
#
#   Parent object call overload
#
###########################################################################################

function Base.call{S}(a::FqNmodFiniteField{S})
   z = fq_nmod{S}(a)
   z.parent = a
   return z
end

Base.call{S}(a::FqNmodFiniteField{S}, b::Integer) = a(BigInt(b))

function Base.call{S}(a::FqNmodFiniteField{S}, b::Int)
   z = fq_nmod{S}(a, b)
   z.parent = a
   return z
end

function Base.call{S}(a::FqNmodFiniteField{S}, b::BigInt)
   z = fq_nmod{S}(a, b)
   z.parent = a
   return z
end

function Base.call{S}(a::FqNmodFiniteField{S}, b::fq_nmod{S})
   parent(b) != a && error("Coercion between finite fields not implemented")
   return b
end

###########################################################################################
#
#   FiniteField constructor
#
###########################################################################################

function FiniteField(char::Int, deg::Int, s::String)
   S = symbol(s)
   parent_obj = FqNmodFiniteField{S}(BigInt(char), deg)

   return parent_obj, gen(parent_obj) 
end

FiniteField(char::Integer, deg::Int, s::String) = FiniteField(BigInt(char), deg, s)