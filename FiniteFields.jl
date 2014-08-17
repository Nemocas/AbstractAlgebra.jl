export FiniteField, gen, prime, degree, order, convert, promote_rule

import Base: convert, promote_rule

import Rings: gen

type fq_ctx
   p :: Int # cannot make this a ZZ
   sparse_modulus :: Int
   a :: Ptr{Void} # cannot make Ptr{ZZ}
   j :: Ptr{Void}
   len :: Int
   # this is actually an fmpz_mod_poly_t
   modulus_coeffs :: Ptr{Void}
   modulus_alloc :: Int
   modulus_length :: Int
   modulus_p :: Int # can't make this a ZZ
   # this is actually an fmpz_mod_poly_
   inv_coeffs :: Ptr{Void}
   inv_alloc :: Int
   inv_length :: Int
   inv_p :: Int # can't make this a ZZ
   var :: Ptr{Void}
   function fq_ctx(p::ZZ, deg::Int, var::String) 
      d = new(0, 0, C_NULL, C_NULL, 0, C_NULL, 0, 0, 0, C_NULL, 0, 0, 0, C_NULL)
      finalizer(d, _fq_ctx_clear_fn)
      ccall((:fq_ctx_init, :libflint), Void, (Ptr{fq_ctx}, Ptr{ZZ}, Int, Ptr{Uint8}), &d, &p, deg, bytestring(var))
      return d
   end
end

function _fq_ctx_clear_fn(a :: fq_ctx)
   ccall((:fq_ctx_clear, :libflint), Void, (Ptr{fq_ctx},), &a)
end

type FField{S} <: Field
   # this is really an fmpz_poly
   coeffs :: Ptr{Void}
   alloc :: Int
   length :: Int
   function FField()
      d = new(C_NULL, 0, 0)
      ccall((:fq_init2, :libflint), Void, (Ptr{FField{S}}, Ptr{fq_ctx}), &d, &eval(:($S)))
      finalizer(d, _FField_clear_fn)
      return d
   end
   function FField(x::Int)
      z = FField{S}()
      ccall((:fq_set_si, :libflint), Void, 
                (Ptr{FField{S}}, Int, Ptr{fq_ctx}), 
               &z, x, &eval(:($S)))
      return z
   end
   function FField(x::ZZ)
      z = FField{S}()
      ccall((:fq_set_fmpz, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{ZZ}, Ptr{fq_ctx}), 
               &z, &x, &eval(:($S)))
      return z
   end
   FField(a::FField{S}) = a
end

function _FField_clear_fn{S}(a :: FField{S})
   ccall((:fq_clear, :libflint), Void, (Ptr{FField{S}}, Ptr{fq_ctx}), &a, &eval(:($S)))
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function zero{S}(::Type{FField{S}})
   d = FField{S}()
   ccall((:fq_zero, :libflint), Void, (Ptr{FField{S}}, Ptr{fq_ctx}), &d, &eval(:($S)))
   return d
end
   
function one{S}(::Type{FField{S}})
   d = FField{S}()
   ccall((:fq_one, :libflint), Void, (Ptr{FField{S}}, Ptr{fq_ctx}), &d, &eval(:($S)))
   return d
end
   
function gen{S}(::Type{FField{S}})
   d = FField{S}()
   ccall((:fq_gen, :libflint), Void, (Ptr{FField{S}}, Ptr{fq_ctx}), &d, &eval(:($S)))
   return d
end

function prime{S}(::Type{FField{S}})
   d = ZZ()
   ccall((:__fq_ctx_prime, :libflint), Void, (Ptr{ZZ}, Ptr{fq_ctx}), &d, &eval(:($S)))
   return d
end
   
function order{S}(::Type{FField{S}})
   d = ZZ()
   ccall((:fq_ctx_order, :libflint), Void, (Ptr{ZZ}, Ptr{fq_ctx}), &d, &eval(:($S)))
   return d
end
   
function degree{S}(::Type{FField{S}})
   return ccall((:fq_ctx_degree, :libflint), Int, (Ptr{fq_ctx}, ), &eval(:($S)))
end
   
###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{S}(io::IO, x::FField{S})
   cstr = ccall((:fq_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{FField{S}}, Ptr{fq_ctx}), &x, &eval(:($S)))

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
end

function show{S}(io::IO, ::Type{FField{S}})
   print(io, "Finite field of degree ")
   print(io, degree(FField{S}))
   print(io, " over F_")
   print(io, prime(FField{S}))
end

needs_parentheses{S}(x::FField{S}) = x.length > 1

is_negative{S}(x::FField{S}) = false

show_minus_one{S}(::Type{FField{S}}) = true

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit{S}(x::FField{S}) = x

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -{S}(x::FField{S})
   z = FField{S}()
   ccall((:fq_neg, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{fq_ctx}), 
               &z, &x, &eval(:($S)))
   return z
end

###########################################################################################
#
#   Binary operators and functions
#
###########################################################################################

function +{S}(x::FField{S}, y::FField{S})
   z = FField{S}()
   ccall((:fq_add, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{FField{S}}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($S)))
   return z
end

function -{S}(x::FField{S}, y::FField{S})
   z = FField{S}()
   ccall((:fq_sub, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{FField{S}}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($S)))
   return z
end

function *{S}(x::FField{S}, y::FField{S})
   z = FField{S}()
   ccall((:fq_mul, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{FField{S}}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($S)))
   return z
end

gcd{S}(x::FField{S}, y::FField{S}) = x == 0 && y == 0 ? FField{S}(0) : FField{S}(1)

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function mul!{S}(z::FField{S}, x::FField{S}, y::FField{S})
   ccall((:fq_mul, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{FField{S}}, Ptr{fq_ctx}), 
               &z, &x, &y)
end

function addeq!{S}(z::FField{S}, x::FField{S})
   ccall((:fq_add, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{FField{S}}, Ptr{fq_ctx}), 
               &z, &z, &x, &eval(:($S)))
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{S}(x::Int, y::FField{S})
   z = FField{S}()
   ccall((:fq_mul_si, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Int, Ptr{fq_ctx}), 
               &z, &y, x, &eval(:($S)))
   return z
end

*{S}(x::FField{S}, y::Int) = y*x

function *{S}(x::ZZ, y::FField{S})
   z = FField{S}()
   ccall((:fq_mul_fmpz, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{ZZ}, Ptr{fq_ctx}), 
               &z, &y, &x, &eval(:($S)))
   return z
end

*{S}(x::FField{S}, y::ZZ) = y*x

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{S}(x::FField{S}, y::Int)
   y < 0 && throw(DomainError())
   z = FField{S}()
   ccall((:fq_pow_ui, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Int, Ptr{fq_ctx}), 
               &z, &x, y, &eval(:($S)))
   return z
end

function ^{S}(x::FField{S}, y::ZZ)
   y < 0 && throw(DomainError())
   z = FField{S}()
   ccall((:fq_pow, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{ZZ}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($S)))
   return z
end

###########################################################################################
#
#   Comparison
#
###########################################################################################

=={S}(x::FField{S}, y::FField{S}) = bool(ccall((:fq_equal, :libflint), Cint, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{fq_ctx}), &x, &y, &eval(:($S))))

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{S}(x::FField{S}, y::FField{S})
   y == 0 && throw(DivideError())
   z = FField{S}()
   ccall((:fq_div, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{FField{S}}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($S)))
   return z
end

/{S}(x::FField{S}, y::FField{S}) = divexact(x, y)

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{S}(x::FField{S})
   z = FField{S}()
   ccall((:fq_inv, :libflint), Void, 
                (Ptr{FField{S}}, Ptr{FField{S}}, Ptr{fq_ctx}), 
               &z, &x, &eval(:($S)))
   return z
end

###########################################################################################
#
#   Conversions and promotions
#
###########################################################################################

convert{S}(::Type{FField{S}}, x::Int) = FField{S}(x)

convert{S}(::Type{FField{S}}, x::ZZ) = FField{S}(x)

promote_rule{S}(::Type{FField{S}}, ::Type{Int}) = FField{S}

promote_rule{S}(::Type{FField{S}}, ::Type{ZZ}) = FField{S}

###########################################################################################
#
#   FiniteField constructor
#
###########################################################################################

function FiniteField(p::ZZ, deg::Int, var::String)
   S = gensym("finitefield")
   P = FField{S}
   eval(:($S = fq_ctx($p, $deg, $var)))
   return P, gen(P)
end

FiniteField(p::Int, deg::Int, var::String) = FiniteField(ZZ(p), deg, var)
