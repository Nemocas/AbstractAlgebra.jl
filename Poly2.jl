###########################################################################################
#
#   Poly2.jl : Polynomials over fields
#
###########################################################################################    

import Rings: Poly, coeff, isgen, truncate, mullow, 
              divexact, gcd, content, primpart, mod, divrem, evaluate, compose, deriv, 
              resultant, bezout, integral, lcm, reverse, shift_left, shift_right, 
              setcoeff!, mulmod, powmod, length

export coeff, isgen, truncate, mullow, divexact, gcd, content, primpart, mod, divrem,
       evaluate, compose, show, deriv, resultant, bezout, integral, lcm, reverse, 
       shift_left, shift_right, setcoeff!, mulmod, powmod, length

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

function Poly{S}(::Type{Poly{QQ, S}}, a :: Array{QQ, 1})
   z = Poly{QQ, S}(Poly)
   finalizer(z, _fmpq_poly_clear_fn)
   ccall((:fmpq_poly_init2, :libflint), Void, (Ptr{Poly}, Int), &z, length(a))
   for i = 1:length(a)
      ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, (Ptr{Poly}, Int, Ptr{Fraction}),
            &z, i - 1, &a[i])
   end
   return z
end

function Poly{S, T}(::Type{Poly{FinFieldElem{T}, S}}, a :: Array{FinFieldElem{T}, 1})
   z = Poly{FinFieldElem{T}, S}(Poly)
   finalizer(z, _fq_poly_clear_fn)
   ctx = eval(:($T))
   ccall((:fq_poly_init2, :libflint), Void, (Ptr{Poly}, Int, Ptr{fq_ctx}), &z, length(a), &ctx)
   for i = 1:length(a)
      ccall((:fq_poly_set_coeff, :libflint), Void, (Ptr{Poly}, Int, Ptr{FinFieldElem{T}}, Ptr{fq_ctx}),
            &z, i - 1, &a[i], &ctx)
   end
   return z
end

function _fmpq_poly_clear_fn{S}(a :: Poly{QQ, S})
   ccall((:fmpq_poly_clear, :libflint), Void, (Ptr{Poly},), &a)
end

function _fq_poly_clear_fn{T, S}(a :: Poly{FinFieldElem{T}, S})
   ccall((:fq_poly_clear, :libflint), Void, (Ptr{Poly}, Ptr{fq_ctx}), &a, eval(:($T)))
end

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{S}(io::IO, x::Poly{QQ, S})
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fmpq_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{Poly}, Ptr{Uint8}), &x, bytestring(string(S)))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
end

function show{S, T}(io::IO, x::Poly{FinFieldElem{T}, S})
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fq_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{Poly}, Ptr{Uint8}, Ptr{fq_ctx}), &x, bytestring(string(S)), &eval(:($T)))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

length{S}(x::Poly{QQ, S}) = ccall((:fmpq_poly_length, :libflint), Int, (Ptr{Poly},), &x)

length{S, T}(x::Poly{FinFieldElem{T}, S}) = ccall((:fq_poly_length, :libflint), Int, (Ptr{Poly}, Ptr{fq_ctx}), &x, &eval(:($T)))

function coeff{S}(x::Poly{QQ, S}, n::Int)
   z = QQ()
   ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, (Ptr{Fraction}, Ptr{Poly}, Int), &z, &x, n)
   return z
end

function coeff{S, T}(x::Poly{FinFieldElem{T}, S}, n::Int)
   z = FinFieldElem{T}()
   ccall((:fq_poly_get_coeff, :libflint), Void, (Ptr{FinFieldElem{T}}, Ptr{Poly}, Int, Ptr{fq_ctx}), &z, &x, n, &eval(:($T)))
   return z
end

isgen{S}(x::Poly{QQ, S}) = bool(ccall((:fmpq_poly_is_x, :libflint), Int, (Ptr{Poly},), &x))

isgen{S, T}(x::Poly{FinFieldElem{T}, S}) = bool(ccall((:fq_poly_is_x, :libflint), Int, (Ptr{Poly}, Ptr{fq_ctx}), &x, &eval(:($T))))

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -{S}(x::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_neg, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}), 
               &z, &x)
   return z
end

function -{S, T}(x::Poly{FinFieldElem{T}, S})
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_neg, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &z, &x, &eval(:($T)))
   return z
end

###########################################################################################
#
#   Binary operations
#
###########################################################################################

function +{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_add, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &z, &x, &y)
   return z
end

function +{S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S})
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_add, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($T)))
   return z
end

function -{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_sub, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &z, &x, &y)
   return z
end

function -{S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S})
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_sub, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($T)))
   return z
end

function *{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_mul, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &z, &x, &y)
   return z
end

function *{S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S})
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_mul, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($T)))
   return z
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function setcoeff!{S}(z::Poly{QQ, S}, n::Int, x::QQ)
   ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, 
                (Ptr{Poly}, Int, Ptr{Fraction}), 
               &z, n, &x)
end

function setcoeff!{S, T}(z::Poly{FinFieldElem{T}, S}, n::Int, x::FinFieldElem{T})
   ccall((:fq_poly_set_coeff, :libflint), Void, 
                (Ptr{Poly}, Int, Ptr{FinFieldElem{T}}, Ptr{fq_ctx}), 
               &z, n, &x, &eval(:($T)))
end

function mul!{S}(z::Poly{QQ, S}, x::Poly{QQ, S}, y::Poly{QQ, S})
   ccall((:fmpq_poly_mul, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &z, &x, &y)
end

function mul!{S, T}(z::Poly{FinFieldElem{T}, S}, x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S})
   ccall((:fq_poly_mul, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($T)))
end

function addeq!{S}(z::Poly{QQ, S}, x::Poly{QQ, S},)
   ccall((:fmpq_poly_add, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &z, &z, &x)
end

function addeq!{S, T}(z::Poly{FinFieldElem{T}, S}, x::Poly{FinFieldElem{T}, S},)
   ccall((:fq_poly_add, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &z, &z, &x, &eval(:($T)))
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{S}(x::Int, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_scalar_mul_si, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Int), 
               &z, &y, x)
   return z
end

function *{S}(x::ZZ, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_scalar_mul_fmpz, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{ZZ}), 
               &z, &y, &x)
   return z
end

function *{S}(x::QQ, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_scalar_mul_fmpq, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Fraction}), 
               &z, &y, &x)
   return z
end

function *{S, T}(x::FinFieldElem{T}, y::Poly{FinFieldElem{T}, S})
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_scalar_mul_fq, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{FinFieldElem{T}}, Ptr{fq_ctx}), 
               &z, &y, &x, &eval(:($T)))
   return z
end

*{S, T}(x::Int, y::Poly{FinFieldElem{T}, S}) = FinFieldElem{T}(x)*y

*{S, T}(x::ZZ, y::Poly{FinFieldElem{T}, S}) = FinFieldElem{T}(x)*y

function +{S}(x::Poly{QQ, S}, y::Int)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_add_si, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Int), 
               &z, &x, y)
   return z
end

function +{S}(x::Poly{QQ, S}, y::ZZ)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_add_fmpz, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{ZZ}), 
               &z, &x, &y)
   return z
end

function +{S}(x::Poly{QQ, S}, y::QQ)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_add_fmpq, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Fraction}), 
               &z, &x, &y)
   return z
end

function -{S}(x::Poly{QQ, S}, y::Int)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_sub_si, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Int), 
               &z, &x, y)
   return z
end

function -{S}(x::Poly{QQ, S}, y::ZZ)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_sub_fmpz, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{ZZ}), 
               &z, &x, &y)
   return z
end

function -{S}(x::Poly{QQ, S}, y::QQ)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_sub_fmpq, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Fraction}), 
               &z, &x, &y)
   return z
end

function -{S}(x::Int, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_si_sub, :libflint), Void, 
                (Ptr{Poly}, Int, Ptr{Poly}), 
               &z, x, &y)
   return z
end

function -{S}(x::ZZ, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_fmpz_sub, :libflint), Void, 
                (Ptr{Poly}, Ptr{ZZ}, Ptr{Poly}), 
               &z, &x, &y)
   return z
end

function -{S}(x::QQ, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_fmpq_sub, :libflint), Void, 
                (Ptr{Poly}, Ptr{Fraction}, Ptr{Poly}), 
               &z, &x, &y)
   return z
end

+{S}(x::Int, y::Poly{QQ, S}) = y + x

+{S}(x::ZZ, y::Poly{QQ, S}) = y + x

+{S}(x::QQ, y::Poly{QQ, S}) = y + x

*{S}(x::Poly{QQ, S}, y::Int) = y*x

*{S}(x::Poly{QQ, S}, y::ZZ) = y*x

*{S}(x::Poly{QQ, S}, y::QQ) = y*x

###########################################################################################
#
#   Truncation
#
###########################################################################################

function truncate{S}(a::Poly{QQ, S}, n::Int)
   n < 0 && throw(DomainError())
   
   if length(a) <= n
      return a
   end

   z = Poly{QQ, S}()
   ccall((:fmpq_poly_set_trunc, :libflint), Void,
                (Ptr{Poly}, Ptr{Poly}, Int),
               &z, &a, n)

   return z
end

function truncate{S, T}(a::Poly{FinFieldElem{T}, S}, n::Int)
   n < 0 && throw(DomainError())
   
   if length(a) <= n
      return a
   end

   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_set_trunc, :libflint), Void,
                (Ptr{Poly}, Ptr{Poly}, Int, Ptr{fq_ctx}),
               &z, &a, n, &eval(:($T)))

   return z
end

function mullow{S}(x::Poly{QQ, S}, y::Poly{QQ, S}, n::Int)
   n < 0 && throw(DomainError())
   
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_mullow, :libflint), Void,
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Int),
               &z, &x, &y, n)
   return z
end

function mullow{S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S}, n::Int)
   n < 0 && throw(DomainError())
   
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_mullow, :libflint), Void,
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Int, Ptr{fq_ctx}),
               &z, &x, &y, n, &eval(:($T)))
   return z
end

###########################################################################################
#
#   Reversal
#
###########################################################################################

function reverse{S}(x::Poly{QQ, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_reverse, :libflint), Void,
                (Ptr{Poly}, Ptr{Poly}, Int),
               &z, &x, len)
   return z
end

function reverse{S, T}(x::Poly{FinFieldElem{T}, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_reverse, :libflint), Void,
                (Ptr{Poly}, Ptr{Poly}, Int, Ptr{fq_ctx}),
               &z, &x, len, &eval(:($T)))
   return z
end

###########################################################################################
#
#   Shifting
#
###########################################################################################

function shift_left{S}(x::Poly{QQ, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_shift_left, :libflint), Void,
                (Ptr{Poly}, Ptr{Poly}, Int),
               &z, &x, len)
   return z
end

function shift_left{S, T}(x::Poly{FinFieldElem{T}, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_shift_left, :libflint), Void,
                (Ptr{Poly}, Ptr{Poly}, Int, Ptr{fq_ctx}),
               &z, &x, len, &eval(:($T)))
   return z
end

function shift_right{S}(x::Poly{QQ, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_shift_right, :libflint), Void,
                (Ptr{Poly}, Ptr{Poly}, Int),
               &z, &x, len)
   return z
end

function shift_right{S, T}(x::Poly{FinFieldElem{T}, S}, len::Int)
   len < 0 && throw(DomainError())
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_shift_right, :libflint), Void,
                (Ptr{Poly}, Ptr{Poly}, Int, Ptr{fq_ctx}),
               &z, &x, len, &eval(:($T)))
   return z
end

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{S}(x::Poly{QQ, S}, y::Int)
   y < 0 && throw(DomainError())
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_pow, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Int), 
               &z, &x, y)
   return z
end

function ^{S, T}(x::Poly{FinFieldElem{T}, S}, y::Int)
   y < 0 && throw(DomainError())
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_pow, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Int, Ptr{fq_ctx}), 
               &z, &x, y, &eval(:($T)))
   return z
end

###########################################################################################
#
#   Comparisons
#
###########################################################################################

function =={S}(x::Poly{QQ, S}, y::ZZ) 
   if length(x) > 1
      return false
   elseif length(x) == 1 
      z = QQ();
      ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, 
                (Ptr{Fraction}, Ptr{Poly}, Int), 
               &z, &x, 0)
      return num(z) == y && den(z) == 1
   else
      return y == 0
   end 
end

=={S}(x::Poly{QQ, S}, y::Int) = x == ZZ(y)

function =={S}(x::Poly{QQ, S}, y::QQ) 
   if length(x) > 1
      return false
   elseif length(x) == 1 
      z = QQ();
      ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, 
                (Ptr{Fraction}, Ptr{Poly}, Int), 
               &z, &x, 0)
      return z == y
   else
      return y == 0
   end 
end

=={S}(x::QQ, y::Poly{QQ, S}) = y == x

function =={S, T}(x::Poly{FinFieldElem{T}, S}, y::FinFieldElem{T}) 
   return bool(ccall((:fq_poly_equal_fq, :libflint), Cint, 
                (Ptr{Poly}, Ptr{FinFieldElem{T}}, Ptr{fq_ctx}), 
               &x, &y, &eval(:($T))))
end

=={S, T}(x::FinFieldElem{T}, y::Poly{FinFieldElem{T}, S}) = y == x

=={S}(x::Poly{QQ, S}, y::Poly{QQ, S}) = ccall((:fmpq_poly_equal, :libflint), Bool, 
                (Ptr{Poly}, Ptr{Poly}), &x, &y)

=={S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S}) = bool(ccall((:fq_poly_equal, :libflint), Cint, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), &x, &y, &eval(:($T))))

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{S}(x::Poly{QQ, S}, y::QQ)
   y == 0 && throw(DivideError())
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_scalar_div_fmpq, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Fraction}), 
               &z,  &x, &y)
   return z
end

function divexact{S}(x::Poly{QQ, S}, y::ZZ)
   y == 0 && throw(DivideError())
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_scalar_div_fmpz, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{ZZ}), 
               &z,  &x, &y)
   return z
end

function divexact{S}(x::Poly{QQ, S}, y::Int)
   y == 0 && throw(DivideError())
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_scalar_div_si, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Int), 
               &z,  &x, y)
   return z
end

function divexact{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   y == 0 && throw(DivideError())
   if x == 0
      return zero(Poly{QQ, S})
   end
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_div, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &z,  &x, &y)
   return z
end

function divexact{S, T}(x::Poly{FinFieldElem{T}, S}, y::FinFieldElem{T})
   y == 0 && throw(DivideError())
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_scalar_div_fq, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{FinFieldElem{T}}, Ptr{fq_ctx}), 
               &z,  &x, &y, &eval(:($T)))
   return z
end

function divexact{S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S})
   q, r = divrem(x, y)
   return q
end

###########################################################################################
#
#   Euclidean division
#
###########################################################################################

function mod{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   y == 0 && throw(DivideError())
   r = Poly{QQ, S}()
   ccall((:fmpq_poly_rem, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &r, &x, &y)
   return r
end

function mod{S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S})
   y == 0 && throw(DivideError())
   r = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_rem, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &r, &x, &y, &eval(:($T)))
   return r
end

function divrem{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   y == 0 && throw(DivideError())
   q = Poly{QQ, S}()
   r = Poly{QQ, S}()
   ccall((:fmpq_poly_divrem, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &q, &r, &x, &y)
   return q, r
end

function divrem{S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S})
   y == 0 && throw(DivideError())
   q = Poly{FinFieldElem{T}, S}()
   r = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_divrem, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &q, &r, &x, &y, &eval(:($T)))
   return q, r
end

###########################################################################################
#
#   Modular arithmetic
#
###########################################################################################

function mulmod{S, T}(f::Poly{FinFieldElem{T}, S}, g::Poly{FinFieldElem{T}, S}, h::Poly{FinFieldElem{T}, S})
   h == 0 && throw(DivideError())
   r = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_mulmod, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &r, &f, &g, &h, &eval(:($T)))
   return r
end

function powmod{S, T}(f::Poly{FinFieldElem{T}, S}, n::Int, h::Poly{FinFieldElem{T}, S})
   h == 0 && throw(DivideError())
   r = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_powmod_ui_binexp, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Int, Ptr{Poly}, Ptr{fq_ctx}), 
               &r, &f, n, &h, &eval(:($T)))
   return r
end

function powmod{S, T}(f::Poly{FinFieldElem{T}, S}, n::ZZ, h::Poly{FinFieldElem{T}, S})
   h == 0 && throw(DivideError())
   r = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_powmod_fmpz_binexp, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{ZZ}, Ptr{Poly}, Ptr{fq_ctx}), 
               &r, &f, &n, &h, &eval(:($T)))
   return r
end

###########################################################################################
#
#   Content, primitive part, GCD and LCM
#
###########################################################################################

function gcd{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_gcd, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &z, &x, &y)
   return z
end

function gcd{S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S})
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_gcd, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($T)))
   return z
end

function lcm{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_lcm, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &z, &x, &y)
   return z
end

function content{T <: Ring, S}(a::Poly{Fraction{T}, S})
   if length(a) == 0
      return zero(Fraction{T})
   end
   d = den(coeff(a, 0))
   for i = 2:length(a)
      z = den(coeff(a, 0))
      d = d*divexact(z, gcd(z, d))
   end
   n = num(coeff(a, 0))*divexact(d, den(coeff(a, 0)))
   for i = 2:length(a)
      z = num(coeff(a, i - 1))*divexact(d, den(coeff(a, i - 1)))
      n = gcd(z, n)
   end
   return n/d
end

function content{S}(x::Poly{QQ, S})
   z = QQ()
   ccall((:fmpq_poly_content, :libflint), Void, 
                (Ptr{Fraction}, Ptr{Poly}), 
               &z, &x)
   return z
end

function primpart{S}(x::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_primitive_part, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}), 
               &z, &x)
   return z
end

###########################################################################################
#
#   Evaluation/composition
#
###########################################################################################

function evaluate{S}(x::Poly{QQ, S}, y::ZZ)
   z = QQ()
   ccall((:fmpq_poly_evaluate_fmpz, :libflint), Void, 
                (Ptr{Fraction}, Ptr{Poly}, Ptr{ZZ}), 
               &z, &x, &y)
   return z
end

evaluate{S}(x::Poly{QQ, S}, y::Int) = evaluate(x, ZZ(y))

function evaluate{S}(x::Poly{QQ, S}, y::QQ)
   z = QQ()
   ccall((:fmpq_poly_evaluate_fmpq, :libflint), Void, 
                (Ptr{Fraction}, Ptr{Poly}, Ptr{fmpq}), 
               &z, &x, &y)
   return z
end

function evaluate{S, T}(x::Poly{FinFieldElem{T}, S}, y::FinFieldElem{T})
   z = FinFieldElem{T}()
   ccall((:fq_poly_evaluate_fq, :libflint), Void, 
                (Ptr{FinFieldElem{T}}, Ptr{Poly}, Ptr{FinFieldElem{T}}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($T)))
   return z
end

function compose{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_compose, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &z, &x, &y)
   return z
end

function compose{S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S})
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_compose, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &z, &x, &y, &eval(:($T)))
   return z
end

###########################################################################################
#
#   Derivative
#
###########################################################################################

function deriv{S}(x::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_derivative, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}), 
               &z, &x)
   return z
end

function deriv{S, T}(x::Poly{FinFieldElem{T}, S})
   z = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_derivative, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &z, &x, &eval(:($T)))
   return z
end

###########################################################################################
#
#   Integral
#
###########################################################################################

function integral{S}(x::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_integral, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}), 
               &z, &x)
   return z
end

###########################################################################################
#
#   Resultant
#
###########################################################################################

function resultant{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = QQ()
   ccall((:fmpq_poly_resultant, :libflint), Void, 
                (Ptr{Fraction}, Ptr{Poly}, Ptr{Poly}), 
               &z, &x, &y)
   return z
end

###########################################################################################
#
#   Bezout
#
###########################################################################################

function bezout{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   u = Poly{QQ, S}()
   v = Poly{QQ, S}()
   ccall((:fmpq_poly_xgcd, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{Poly}), 
               &z, &u, &v, &x, &y)
   return (z, u, v)
end

function bezout{S, T}(x::Poly{FinFieldElem{T}, S}, y::Poly{FinFieldElem{T}, S})
   z = Poly{FinFieldElem{T}, S}()
   u = Poly{FinFieldElem{T}, S}()
   v = Poly{FinFieldElem{T}, S}()
   ccall((:fq_poly_xgcd, :libflint), Void, 
                (Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{Poly}, Ptr{fq_ctx}), 
               &z, &u, &v, &x, &y, &eval(:($T)))
   return (z, u, v)
end


