export FiniteField, FinFieldElem, gen, characteristic, degree, order, convert, promote_rule, pth_root, 
       trace, norm, frobenius, FinFieldCtx

import Base: convert, promote_rule

import Rings: gen

FinFieldCtx = ObjectIdDict()

type fq_ctx
   p :: Int # cannot make this a ZZ
   sparse_modulus :: Int
   a :: Ptr{Void} # cannot make Ptr{ZZ}
   j :: Ptr{Void}
   len :: Int
   modulus::fmpz_mod_poly_struct
   inv::fmpz_mod_poly_struct
   var :: Ptr{Void}
   function fq_ctx(p::ZZ, deg::Int, var::String) 
      d = new(0, 0, C_NULL, C_NULL, 0, fmpz_mod_poly_struct(), fmpz_mod_poly_struct(), C_NULL)
      finalizer(d, _fq_ctx_clear_fn)
      ccall((:fq_ctx_init, :libflint), Void, (Ptr{fq_ctx}, Ptr{ZZ}, Int, Ptr{Uint8}), &d, &p, deg, bytestring(var))
      return d
   end
end

function _fq_ctx_clear_fn(a :: fq_ctx)
   ccall((:fq_ctx_clear, :libflint), Void, (Ptr{fq_ctx},), &a)
end

type FinFieldElem{S} <: Field
   data::fmpz_poly_struct
   function FinFieldElem()
      d = new(fmpz_poly_struct())
      ccall((:fq_init2, :libflint), Void, (Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), &d, &FinFieldCtx[S])
      finalizer(d, _FinFieldElem_clear_fn)
      return d
   end
   function FinFieldElem(x::Int)
      z = FinFieldElem{S}()
      ccall((:fq_set_si, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Int, Ptr{fq_ctx}), 
               &z, x, &FinFieldCtx[S])
      return z
   end
   function FinFieldElem(x::ZZ)
      z = FinFieldElem{S}()
      ccall((:fq_set_fmpz, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{ZZ}, Ptr{fq_ctx}), 
               &z, &x, &FinFieldCtx[S])
      return z
   end
   FinFieldElem(a::FinFieldElem{S}) = a
end

function _FinFieldElem_clear_fn{S}(a :: FinFieldElem{S})
   ccall((:fq_clear, :libflint), Void, (Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), &a, &FinFieldCtx[S])
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function zero{S}(::Type{FinFieldElem{S}})
   d = FinFieldElem{S}()
   ccall((:fq_zero, :libflint), Void, (Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), &d, &FinFieldCtx[S])
   return d
end
   
function one{S}(::Type{FinFieldElem{S}})
   d = FinFieldElem{S}()
   ccall((:fq_one, :libflint), Void, (Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), &d, &FinFieldCtx[S])
   return d
end
   
function gen{S}(::Type{FinFieldElem{S}})
   d = FinFieldElem{S}()
   ccall((:fq_gen, :libflint), Void, (Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), &d, &FinFieldCtx[S])
   return d
end

iszero{S}(a::FinFieldElem{S}) = bool(ccall((:fq_is_zero, :libflint), Cint, (Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), &a, &FinFieldCtx[S]))

isone{S}(a::FinFieldElem{S}) = bool(ccall((:fq_is_one, :libflint), Cint, (Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), &a, &FinFieldCtx[S]))

isunit{S}(a::FinFieldElem{S}) = bool(ccall((:fq_is_invertible, :libflint), Cint, (Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), &a, &FinFieldCtx[S]))

function characteristic{S}(::Type{FinFieldElem{S}})
   d = ZZ()
   ccall((:__fq_ctx_prime, :libflint), Void, (Ptr{ZZ}, Ptr{fq_ctx}), &d, &FinFieldCtx[S])
   return d
end
   
function order{S}(::Type{FinFieldElem{S}})
   d = ZZ()
   ccall((:fq_ctx_order, :libflint), Void, (Ptr{ZZ}, Ptr{fq_ctx}), &d, &FinFieldCtx[S])
   return d
end
   
function degree{S}(::Type{FinFieldElem{S}})
   return ccall((:fq_ctx_degree, :libflint), Int, (Ptr{fq_ctx}, ), &FinFieldCtx[S])
end
   
###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{S}(io::IO, x::FinFieldElem{S})
   cstr = ccall((:fq_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), &x, &FinFieldCtx[S])

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
end

function show{S}(io::IO, ::Type{FinFieldElem{S}})
   print(io, "Finite field of degree ")
   print(io, degree(FinFieldElem{S}))
   print(io, " over F_")
   print(io, characteristic(FinFieldElem{S}))
end

needs_parentheses{S}(x::FinFieldElem{S}) = x.length > 1

is_negative{S}(x::FinFieldElem{S}) = false

show_minus_one{S}(::Type{FinFieldElem{S}}) = true

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit{S}(x::FinFieldElem{S}) = x

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -{S}(x::FinFieldElem{S})
   z = FinFieldElem{S}()
   ccall((:fq_neg, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &x, &FinFieldCtx[S])
   return z
end

###########################################################################################
#
#   Binary operators and functions
#
###########################################################################################

function +{S}(x::FinFieldElem{S}, y::FinFieldElem{S})
   z = FinFieldElem{S}()
   ccall((:fq_add, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &x, &y, &FinFieldCtx[S])
   return z
end

function -{S}(x::FinFieldElem{S}, y::FinFieldElem{S})
   z = FinFieldElem{S}()
   ccall((:fq_sub, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &x, &y, &FinFieldCtx[S])
   return z
end

function *{S}(x::FinFieldElem{S}, y::FinFieldElem{S})
   z = FinFieldElem{S}()
   ccall((:fq_mul, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &x, &y, &FinFieldCtx[S])
   return z
end

gcd{S}(x::FinFieldElem{S}, y::FinFieldElem{S}) = x == 0 && y == 0 ? FinFieldElem{S}(0) : FinFieldElem{S}(1)

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function mul!{S}(z::FinFieldElem{S}, x::FinFieldElem{S}, y::FinFieldElem{S})
   ccall((:fq_mul, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &x, &y, &FinFieldCtx[S])
end

function addeq!{S}(z::FinFieldElem{S}, x::FinFieldElem{S})
   ccall((:fq_add, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &z, &x, &FinFieldCtx[S])
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{S}(x::Int, y::FinFieldElem{S})
   z = FinFieldElem{S}()
   ccall((:fq_mul_si, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Int, Ptr{fq_ctx}), 
               &z, &y, x, &FinFieldCtx[S])
   return z
end

*{S}(x::FinFieldElem{S}, y::Int) = y*x

function *{S}(x::ZZ, y::FinFieldElem{S})
   z = FinFieldElem{S}()
   ccall((:fq_mul_fmpz, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{ZZ}, Ptr{fq_ctx}), 
               &z, &y, &x, &FinFieldCtx[S])
   return z
end

*{S}(x::FinFieldElem{S}, y::ZZ) = y*x

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{S}(x::FinFieldElem{S}, y::Int)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = FinFieldElem{S}()
   ccall((:fq_pow_ui, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Int, Ptr{fq_ctx}), 
               &z, &x, y, &FinFieldCtx[S])
   return z
end

function ^{S}(x::FinFieldElem{S}, y::ZZ)
   if y < 0
      x = inv(x)
      y = -y
   end
   z = FinFieldElem{S}()
   ccall((:fq_pow, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{ZZ}, Ptr{fq_ctx}), 
               &z, &x, &y, &FinFieldCtx[S])
   return z
end

###########################################################################################
#
#   Comparison
#
###########################################################################################

=={S}(x::FinFieldElem{S}, y::FinFieldElem{S}) = bool(ccall((:fq_equal, :libflint), Cint, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), &x, &y, &FinFieldCtx[S]))

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{S}(x::FinFieldElem{S}, y::FinFieldElem{S})
   y == 0 && throw(DivideError())
   z = FinFieldElem{S}()
   ccall((:fq_div, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &x, &y, &FinFieldCtx[S])
   return z
end

/{S}(x::FinFieldElem{S}, y::FinFieldElem{S}) = divexact(x, y)

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{S}(x::FinFieldElem{S})
   x == 0 && throw(DivideError())
   z = FinFieldElem{S}()
   ccall((:fq_inv, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &x, &FinFieldCtx[S])
   return z
end

###########################################################################################
#
#   Special functions
#
###########################################################################################

function pth_root{S}(x::FinFieldElem{S})
   z = FinFieldElem{S}()
   ccall((:fq_pth_root, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &x, &FinFieldCtx[S])
   return z
end

function trace{S}(x::FinFieldElem{S})
   z = ZZ()
   ccall((:fq_trace, :libflint), Void, 
                (Ptr{ZZ}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &x, &FinFieldCtx[S])
   return FinFieldElem{S}(z)
end

function norm{S}(x::FinFieldElem{S})
   z = ZZ()
   ccall((:fq_norm, :libflint), Void, 
                (Ptr{ZZ}, Ptr{FinFieldElem{S}}, Ptr{fq_ctx}), 
               &z, &x, &FinFieldCtx[S])
   return FinFieldElem{S}(z)
end

function frobenius{S}(x::FinFieldElem{S}, n = 1)
   z = FinFieldElem{S}()
   ccall((:fq_frobenius, :libflint), Void, 
                (Ptr{FinFieldElem{S}}, Ptr{FinFieldElem{S}}, Int, Ptr{fq_ctx}), 
               &z, &x, n, &FinFieldCtx[S])
   return z
end

###########################################################################################
#
#   Conversions and promotions
#
###########################################################################################

convert{S}(::Type{FinFieldElem{S}}, x::Int) = FinFieldElem{S}(x)

convert{S}(::Type{FinFieldElem{S}}, x::ZZ) = FinFieldElem{S}(x)

promote_rule{S}(::Type{FinFieldElem{S}}, ::Type{Int}) = FinFieldElem{S}

promote_rule{S}(::Type{FinFieldElem{S}}, ::Type{ZZ}) = FinFieldElem{S}

###########################################################################################
#
#   FiniteField constructor
#
###########################################################################################

function FiniteField(p::ZZ, deg::Int, var::String)
   deg < 0 && throw(DomainError())
   S = gensym("finitefield")
   P = FinFieldElem{S}
   FinFieldCtx[S] = fq_ctx(p, deg, var)
   return P, gen(P)
end

FiniteField(p::Int, deg::Int, var::String) = FiniteField(ZZ(p), deg, var)
