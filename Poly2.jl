###########################################################################################
#
#   Poly2.jl : Polynomials over fields
#
###########################################################################################    

import Rings: Poly, fmpq_poly, coeff, isgen, truncate, mullow, divexact, gcd, content,
              primpart, mod, divrem, evaluate, compose, deriv, resultant, bezout, integral,
              lcm

export coeff, isgen, truncate, mullow, divexact, gcd, content, primpart, mod, divrem,
       evaluate, compose, show, deriv, resultant, bezout, integral, lcm

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

function Poly{S}(::Type{Poly{QQ, S}}, a :: Array{QQ, 1})
   z = fmpq_poly(C_NULL, 0, 0, 0)
   ccall((:fmpq_poly_init2, :libflint), Void, (Ptr{fmpq_poly}, Int), &z, length(a))
   for i = 1:length(a)
      ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, (Ptr{fmpq_poly}, Int, Ptr{fmpq}),
            &z, i - 1, &(a[i].data))
   end
   return Poly{QQ, S}(z)
end

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{S}(io::IO, x::Poly{QQ, S})
   if x.data.length == 0
      print(io, "0")
   else
      cstr = ccall((:fmpq_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{fmpq_poly}, Ptr{Uint8}), &(x.data), bytestring(string(S)))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function coeff{S}(x::Poly{QQ, S}, n::Int)
   z = QQ()
   ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq_poly}, Int), &(z.data), &(x.data), n)
   return z
end

isgen{S}(x::Poly{QQ, S}) = bool(ccall((:fmpq_poly_is_x, :libflint), Int, (Ptr{fmpq_poly},), &(x.data)))

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -{S}(x::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_neg, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data))
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
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function -{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function *{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_mul, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function setcoeff!{S}(z::Poly{QQ, S}, n::Int, x::QQ)
   ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Int, Ptr{fmpq}), 
               &(z.data), n, &(x.data))
end

function mul!{S}(z::Poly{QQ, S}, x::Poly{QQ, S}, y::Poly{QQ, S})
   ccall((:fmpq_poly_mul, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data), &(y.data))
end

function addeq!{S}(z::Poly{QQ, S}, x::Poly{QQ, S},)
   ccall((:fmpq_poly_add, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(z.data), &(x.data))
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{S}(x::Int, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_scalar_mul_si, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), 
               &(z.data), &(y.data), x)
   return z
end

function *{S}(x::ZZ, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_scalar_mul_fmpz, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{ZZ}), 
               &(z.data), &(y.data), &x)
   return z
end

function *{S}(x::QQ, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_scalar_mul_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq}), 
               &(z.data), &(y.data), &(x.data))
   return z
end

function +{S}(x::Poly{QQ, S}, y::Int)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_add_si, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), 
               &(z.data), &(x.data), y)
   return z
end

function +{S}(x::Poly{QQ, S}, y::ZZ)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_add_fmpz, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{ZZ}), 
               &(z.data), &(x.data), &y)
   return z
end

function +{S}(x::Poly{QQ, S}, y::QQ)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_add_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function -{S}(x::Poly{QQ, S}, y::Int)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_sub_si, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), 
               &(z.data), &(x.data), y)
   return z
end

function -{S}(x::Poly{QQ, S}, y::ZZ)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_sub_fmpz, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{ZZ}), 
               &(z.data), &(x.data), &y)
   return z
end

function -{S}(x::Poly{QQ, S}, y::QQ)
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_sub_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function -{S}(x::Int, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_si_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Int, Ptr{fmpq_poly}), 
               &(z.data), x, &(y.data))
   return z
end

function -{S}(x::ZZ, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_fmpz_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{ZZ}, Ptr{fmpq_poly}), 
               &(z.data), &x, &(y.data))
   return z
end

function -{S}(x::QQ, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_fmpq_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data), &(y.data))
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
   
   if a.data.length <= n
      return a
   end

   z = Poly{QQ, S}()
   ccall((:fmpq_poly_set_trunc, :libflint), Void,
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int),
               &(z.data), &(a.data), n)

   return z
end

function mullow{S}(x::Poly{QQ, S}, y::Poly{QQ, S}, n::Int)
   n < 0 && throw(DomainError())
   
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_mullow, :libflint), Void,
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int),
               &(z.data), &(x.data), &(y.data), n)
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
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), 
               &(z.data), &(x.data), y)
   return z
end

###########################################################################################
#
#   Comparisons
#
###########################################################################################

function =={S}(x::Poly{QQ, S}, y::ZZ) 
   if x.data.length > 1
      return false
   elseif x.data.length == 1 
      z = QQ();
      ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}, Int), 
               &(z.data), &(x.data), 0)
      return num(z) == y && den(z) == 1
   else
      return y == 0
   end 
end

function =={S}(x::Poly{QQ, S}, y::QQ) 
   if x.data.length > 1
      return false
   elseif x.data.length == 1 
      z = QQ();
      ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}, Int), 
               &(z.data), &(x.data), 0)
      return z == y
   else
      return y == 0
   end 
end

=={S}(x::Poly{QQ, S}, y::Int) = x == ZZ(y)

=={S}(x::Poly{QQ, S}, y::Poly{QQ, S}) = ccall((:fmpq_poly_equal, :libflint), Bool, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &(x.data), &(y.data))

=={S}(x::QQ, y::Poly{QQ, S}) = y == x

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{S}(x::Poly{QQ, S}, y::QQ)
   y == 0 && throw(DivideError())
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_scalar_div_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq}), 
               &(z.data),  &(x.data), &(y.data))
   return z
end

function divexact{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   y == 0 && throw(DivideError())
   if x == 0
      return zero(Poly{QQ, S})
   end
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_div, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data),  &(x.data), &(y.data))
   return z
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
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(r.data), &(x.data), &(y.data))
   return r
end

function divrem{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   y == 0 && throw(DivideError())
   if x == 0
      return zero(Poly{QQ, S})
   end
   q = Poly{QQ, S}()
   r = Poly{QQ, S}()
   ccall((:fmpq_poly_divrem, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(q.data), &(r.data), &(x.data), &(y.data))
   return q, r
end

###########################################################################################
#
#   Content, primitive part, GCD and LCM
#
###########################################################################################

function gcd{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_gcd, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function lcm{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_lcm, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function content{T <: Ring, S}(a::Poly{Fraction{T}, S})
   if a.data.length == 0
      return zero(Fraction{T})
   end
   d = den(coeff(a, 0))
   for i = 2:a.data.length
      z = den(coeff(a, 0))
      d = d*divexact(z, gcd(z, d))
   end
   n = num(coeff(a, 0))*divexact(d, den(coeff(a, 0)))
   for i = 2:a.data.length
      z = num(coeff(a, i - 1))*divexact(d, den(coeff(a, i - 1)))
      n = gcd(z, n)
   end
   return n/d
end

function content{S}(x::Poly{QQ, S})
   z = QQ()
   ccall((:fmpq_poly_content, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data))
   return z
end

function primpart{S}(x::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_primitive_part, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data))
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
                (Ptr{fmpq}, Ptr{fmpq_poly}, Ptr{ZZ}), 
               &(z.data), &(x.data), &y)
   return z
end

evaluate{S}(x::Poly{QQ, S}, y::Int) = evaluate(x, ZZ(y))

function evaluate{S}(x::Poly{QQ, S}, y::QQ)
   z = QQ()
   ccall((:fmpq_poly_evaluate_fmpq, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}, Ptr{fmpq}), 
               &(z.data), &(x.data), &(y.data))
   return z
end

function compose{S}(x::Poly{QQ, S}, y::Poly{QQ, S})
   z = Poly{QQ, S}()
   ccall((:fmpq_poly_compose, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data), &(y.data))
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
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data))
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
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data))
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
                (Ptr{fmpq}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(x.data), &(y.data))
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
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &(z.data), &(u.data), &(v.data), &(x.data), &(y.data))
   return (z, u, v)
end


