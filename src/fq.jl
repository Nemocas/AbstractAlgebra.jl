###########################################################################################
#
#   fq.jl : Flint finite fields
#
###########################################################################################

export fq, FqFiniteField

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

FqFiniteFieldID = Dict{(BigInt, Int, Symbol), Field}()

type FqFiniteField{S} <: Field
   p :: Int # fmpz
   sparse_modulus :: Int
   a :: Ptr{Void}
   j :: Ptr{Void}
   len :: Int
   mod_coeffs :: Ptr{Void}
   mod_alloc :: Int
   mod_length :: Int
   mod_p :: Int # fmpz
   inv_coeffs :: Ptr{Void}
   inv_alloc :: Int
   inv_length :: Int
   inv_p :: Int # fmpz
   var :: Ptr{Void}

   function FqFiniteField(char::BigInt, deg::Int)
      try
         return FqFiniteFieldID[char, deg, S]
      catch
         d = FqFiniteFieldID[char, deg, S] = new()
         finalizer(d, _FqFiniteField_clear_fn)
         temp = fmpz_readonly(char)
         ccall((:fq_ctx_init, :libflint), Void, 
               (Ptr{FqFiniteField}, Ptr{fmpz_readonly}, Int, Ptr{Uint8}), 
			    &d, &temp, deg, bytestring(string(S)))
         return d
      end
   end
end

function _FqFiniteField_clear_fn(a :: FqFiniteField)
   ccall((:fq_ctx_clear, :libflint), Void, (Ptr{FqFiniteField},), &a)
end

type fq{S} <: FiniteFieldElem
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
   parent::FqFiniteField{S}

   function fq(ctx::FqFiniteField)
      d = new()
      ccall((:fq_init2, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      return d
   end

   function fq(ctx::FqFiniteField, x::Int)
      d = new()
      ccall((:fq_init2, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      ccall((:fq_set_si, :libflint), Void, 
                (Ptr{fq}, Int, Ptr{FqFiniteField}), &d, x, &ctx)
      return d
   end

   function fq(ctx::FqFiniteField, x::BigInt)
      d = new()
      ccall((:fq_init2, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &d, &ctx)
      finalizer(d, _fq_clear_fn)
      temp = fmpz_readonly(x)
      ccall((:fq_set_fmpz, :libflint), Void, 
                (Ptr{fq}, Ptr{fmpz_readonly}, Ptr{FqFiniteField}), &d, &temp, &ctx)
      return d
   end
end

function _fq_clear_fn(a::fq)
   ccall((:fq_clear, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &a, &a.parent)
end

elem_type{S}(::FqFiniteField{S}) = fq{S}

base_ring(a::FqFiniteField) = None

base_ring(a::fq) = None

parent(a::fq) = a.parent

function check_parent(a::fq, b::fq) 
   a.parent != b.parent && error("Operations on distinct finite fields not supported")
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function zero{S}(a::FqFiniteField{S})
   d = a()
   ccall((:fq_zero, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &d, &a)
   return d
end
   
function one{S}(a::FqFiniteField{S})
   d = a()
   ccall((:fq_one, :libflint), Void, (Ptr{fq}, Ptr{FqFiniteField}), &d, &a)
   return d
end
   
function gen{S}(a::FqFiniteField{S})
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
   d = fmpz()
   ccall((:__fq_ctx_prime, :libflint), Void, (Ptr{fmpz}, Ptr{FqFiniteField}), &d, &a)
   r = BigInt(d)
   return r
end
   
function order(a::FqFiniteField)
   d = fmpz()
   ccall((:fq_ctx_order, :libflint), Void, (Ptr{fmpz}, Ptr{FqFiniteField}), &d, &a)
   r = BigInt(d)
   return r
end
   
function degree(a::FqFiniteField)
   return ccall((:fq_ctx_degree, :libflint), Int, (Ptr{FqFiniteField},), &a)
end
###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit(x::fq) = x

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, x::fq)
   cstr = ccall((:fq_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{fq}, Ptr{FqFiniteField}), &x, &x.parent)

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
end

function show(io::IO, a::FqFiniteField)
   print(io, "Finite field of degree ", degree(a), " over F_", characteristic(a))
end

needs_parentheses(x::fq) = x.length > 1

is_negative(x::fq) = false

show_minus_one{S}(::Type{fq{S}}) = true

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -(x::fq)
   z = parent(x)()
   ccall((:fq_neg, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   return z
end

###########################################################################################
#
#   Binary operations
#
###########################################################################################

function +{S}(x::fq{S}, y::fq{S})
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_add, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &y, &y.parent)
   return z
end

function -{S}(x::fq{S}, y::fq{S})
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_sub, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &y, &y.parent)
   return z
end

function *{S}(x::fq{S}, y::fq{S})
   check_parent(x, y)
   z = parent(y)()
   ccall((:fq_mul, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &y, &y.parent)
   return z
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *(x::Int, y::fq)
   z = parent(y)()
   ccall((:fq_mul_si, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Int, Ptr{FqFiniteField}), &z, &y, x, &y.parent)
   return z
end

*(x::Integer, y::fq) = BigInt(x)*y

*(x::fq, y::Integer) = y*x

function *(x::BigInt, y::fq)
   z = parent(y)()
   temp = fmpz_readonly(x)
   ccall((:fq_mul_fmpz, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{fmpz_readonly}, Ptr{FqFiniteField}), &z, &y, &temp, &y.parent)
   return z
end

*(x::fq, y::BigInt) = y*x

###########################################################################################
#
#   Powering
#
###########################################################################################

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

function ^(x::fq, y::BigInt)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = parent(x)()
   temp = fmpz_readonly(y)
   ccall((:fq_pow, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{fmpz_readonly}, Ptr{FqFiniteField}), &z, &x, &temp, &x.parent)
   return z
end

###########################################################################################
#
#   Comparison
#
###########################################################################################

function =={S}(x::fq{S}, y::fq{S}) 
   check_parent(x, y)
   ccall((:fq_equal, :libflint), Bool, 
         (Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &x, &y, &y.parent)
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{S}(x::fq{S}, y::fq{S})
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   z = parent(y)()
   ccall((:fq_div, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &y, &y.parent)
   return z
end

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv(x::fq)
   iszero(x) && throw(DivideError())
   z = parent(x)()
   ccall((:fq_inv, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   return z
end

###########################################################################################
#
#   Special functions
#
###########################################################################################

function pth_root(x::fq)
   z = parent(x)()
   ccall((:fq_pth_root, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   return z
end

function trace(x::fq)
   z = fmpz()
   init(z)
   ccall((:fq_trace, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   r = BigInt(z)
   clear(z)
   return parent(x)(r)
end

function norm(x::fq)
   z = fmpz()
   init(z)
   ccall((:fq_norm, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &x.parent)
   r = BigInt(z)
   clear(z)
   return parent(x)(r)
end

function frobenius(x::fq, n = 1)
   z = parent(x)()
   ccall((:fq_frobenius, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Int, Ptr{FqFiniteField}), &z, &x, n, &x.parent)
   return z
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function mul!(z::fq, x::fq, y::fq)
   ccall((:fq_mul, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &x, &y, &y.parent)
end

function addeq!(z::fq, x::fq)
   ccall((:fq_add, :libflint), Void, 
         (Ptr{fq}, Ptr{fq}, Ptr{fq}, Ptr{FqFiniteField}), &z, &z, &x, &x.parent)
end

###########################################################################################
#
#   Promotions
#
###########################################################################################

Base.promote_rule{S, T <: Integer}(::Type{fq{S}}, ::Type{T}) = fq{S}

Base.promote_rule{S}(::Type{fq{S}}, ::Type{BigInt}) = fq{S}

###########################################################################################
#
#   Parent object call overload
#
###########################################################################################

function Base.call{S}(a::FqFiniteField{S})
   z = fq{S}(a)
   z.parent = a
   return z
end

Base.call{S}(a::FqFiniteField{S}, b::Integer) = a(BigInt(b))

function Base.call{S}(a::FqFiniteField{S}, b::Int)
   z = fq{S}(a, b)
   z.parent = a
   return z
end

function Base.call{S}(a::FqFiniteField{S}, b::BigInt)
   z = fq{S}(a, b)
   z.parent = a
   return z
end

function Base.call{S}(a::FqFiniteField{S}, b::fq{S})
   parent(b) != a && error("Coercion between finite fields not implemented")
   return b
end

###########################################################################################
#
#   FiniteField constructor
#
###########################################################################################

function FiniteField(char::BigInt, deg::Int, s::String)
   S = symbol(s)
   parent_obj = FqFiniteField{S}(char, deg)

   return parent_obj, gen(parent_obj) 
end

FiniteField(char::Integer, deg::Int, s::String) = FiniteField(BigInt(char), deg, s)