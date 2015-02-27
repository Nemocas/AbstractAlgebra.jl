###########################################################################################
#
#   fmpq_poly.jl : Flint polynomials over QQ
#
###########################################################################################

export fmpq_poly, denominator

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

FmpqPolyID = ObjectIdDict()

type FmpqPolyRing{S} <: Ring
   base_ring::Field

   function FmpqPolyRing(R::RationalField)
      return try
         FmpqPolyID[S]
      catch
         FmpqPolyID[S] = new(R)
      end
   end
end

type fmpq_poly{S} <: PolyElem
   coeffs::Ptr{Int}
   den::Int # really an fmpz
   alloc::Int
   length::Int
   parent::FmpqPolyRing{S}

   function fmpq_poly()
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::Array{Rational{BigInt}, 1})
      z = new()
      ccall((:fmpq_poly_init2, :libflint), Void, (Ptr{fmpq_poly}, Int), &z, length(a))
      for i = 1:length(a)
         temp = fmpq_readonly(a[i])
         ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, 
                     (Ptr{fmpq_poly}, Int, Ptr{fmpq_readonly}), &z, i - 1, &temp)
      end
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::Int)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_si, :libflint), Void, (Ptr{fmpq_poly}, Int), &z, a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::BigInt)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_mpz, :libflint), Void, (Ptr{fmpq_poly}, Ptr{BigInt}), &z, &a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::Rational{BigInt})
      z = new()
      temp = fmpq_readonly(a)
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_fmpq, :libflint), Void, (Ptr{fmpq_poly}, Ptr{fmpq_readonly}), &z, &temp)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::fmpz_poly)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_fmpz_poly, :libflint), Void, (Ptr{fmpq_poly}, Ptr{fmpz_poly}), &z, &a)
      return z
   end
end

function _fmpq_poly_clear_fn(a::fmpq_poly)
   ccall((:fmpq_poly_clear, :libflint), Void, (Ptr{fmpq_poly},), &a)
end

elem_type{S}(::FmpqPolyRing{S}) = fmpq_poly{S}

base_ring(a::FmpqPolyRing) = a.base_ring

function denominator(a::fmpq_poly)
   z = fmpz()
   ccall((:fmpq_poly_get_denominator, :libflint), Void,
         (Ptr{fmpz}, Ptr{fmpq_poly}), &z, &a)
   return BigInt(z)
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################    
   
length{S}(x::fmpq_poly{S}) = ccall((:fmpq_poly_length, :libflint), Int, 
                                   (Ptr{fmpq_poly},), &x)

function coeff(x::fmpq_poly, n::Int)
   n < 0 && throw(DomainError())
   temp = fmpq()
   ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, 
               (Ptr{fmpq}, Ptr{fmpq_poly}, Int), &temp, &x, n)
   return Rational(temp)
end

zero(a::FmpqPolyRing) = a(0)

one(a::FmpqPolyRing) = a(1)

gen(a::FmpqPolyRing) = a([zero(base_ring(a)), one(base_ring(a))])

isgen(x::fmpq_poly) = ccall((:fmpq_poly_is_x, :libflint), Bool, (Ptr{fmpq_poly},), &x)

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit(a::fmpq_poly) = canonical_unit(lead(a))

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{S}(io::IO, x::fmpq_poly{S})
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fmpq_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{fmpq_poly{S}}, Ptr{Uint8}), &x, bytestring(string(S)))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
end

function show{S}(io::IO, p::FmpqPolyRing{S})
   print(io, "Univariate Polynomial Ring in ")
   print(io, string(S))
   print(io, " over ")
   show(io, p.base_ring)
end

show_minus_one(::Type{fmpq_poly}) = show_minus_one(RationalField)

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -(x::fmpq_poly)
   z = parent(x)()
   ccall((:fmpq_poly_neg, :libflint), Void, (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x)
   return z
end

###########################################################################################
#
#   Binary operations
#
###########################################################################################

function +{S}(x::fmpq_poly{S}, y::fmpq_poly{S})
   z = parent(x)()
   ccall((:fmpq_poly_add, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly},  Ptr{fmpq_poly}), 
               &z, &x, &y)
   return z
end

function -{S}(x::fmpq_poly{S}, y::fmpq_poly{S})
   z = parent(x)()
   ccall((:fmpq_poly_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly},  Ptr{fmpq_poly}), 
               &z, &x, &y)
   return z
end

function *{S}(x::fmpq_poly{S}, y::fmpq_poly{S})
   z = parent(x)()
   ccall((:fmpq_poly_mul, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly},  Ptr{fmpq_poly}), 
               &z, &x, &y)
   return z
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *(x::Int, y::fmpq_poly)
   z = parent(y)()
   ccall((:fmpq_poly_scalar_mul_si, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &y, x)
   return z
end

function *(x::BigInt, y::fmpq_poly)
   z = parent(y)()
   temp = fmpz_readonly(x)
   ccall((:fmpq_poly_scalar_mul_fmpz, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpz_readonly}), &z, &y, &temp)
   return z
end

function *(x::Rational{BigInt}, y::fmpq_poly)
   z = parent(y)()
   temp = fmpq_readonly(x)
   ccall((:fmpq_poly_scalar_mul_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_readonly}), &z, &y, &temp)
   return z
end

function +(x::fmpq_poly, y::Int)
   z = parent(x)()
   ccall((:fmpq_poly_add_si, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, y)
   return z
end

function +(x::fmpq_poly, y::BigInt)
   z = parent(x)()
   temp = fmpz_readonly(y)
   ccall((:fmpq_poly_add_fmpz, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpz_readonly}), &z, &x, &temp)
   return z
end

function +(x::fmpq_poly, y::Rational{BigInt})
   z = parent(x)()
   temp = fmpq_readonly(y)
   ccall((:fmpq_poly_add_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_readonly}), &z, &x, &temp)
   return z
end

function -(x::fmpq_poly, y::Int)
   z = parent(x)()
   ccall((:fmpq_poly_sub_si, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, y)
   return z
end

function -(x::fmpq_poly, y::BigInt)
   z = parent(x)()
   temp = fmpz_readonly(y)
   ccall((:fmpq_poly_sub_fmpz, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpz_readonly}), &z, &x, &temp)
   return z
end

function -(x::fmpq_poly, y::Rational{BigInt})
   z = parent(x)()
   temp = fmpq_readonly(y)
   ccall((:fmpq_poly_sub_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_readonly}), &z, &x, &temp)
   return z
end

function -(x::Int, y::fmpq_poly)
   z = parent(y)()
   ccall((:fmpq_poly_si_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Int, Ptr{fmpq_poly}), &z, x, &y)
   return z
end

function -(x::BigInt, y::fmpq_poly)
   z = parent(y)()
   temp = fmpz_readonly(x)
   ccall((:fmpq_poly_fmpz_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpz_readonly}, Ptr{fmpq_poly}), &z, &temp, &y)
   return z
end

function -(x::Rational{BigInt}, y::fmpq_poly)
   z = parent(y)()
   temp = fmpq_readonly(x)
   ccall((:fmpq_poly_fmpq_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_readonly}, Ptr{fmpq_poly}), &z, &temp, &y)
   return z
end

+(x::Int, y::fmpq_poly) = y + x

+(x::BigInt, y::fmpq_poly) = y + x

+(x::Rational{BigInt}, y::fmpq_poly) = y + x

*(x::fmpq_poly, y::Int) = y*x

*(x::fmpq_poly, y::BigInt) = y*x

*(x::fmpq_poly, y::Rational{BigInt}) = y*x

+(x::Integer, y::fmpq_poly) = y + BigInt(x)

-(x::Integer, y::fmpq_poly) = BigInt(x) - y

*(x::Integer, y::fmpq_poly) = BigInt(x)*y

+(x::fmpq_poly, y::Integer) = x + BigInt(y)

-(x::fmpq_poly, y::Integer) = x - BigInt(y)

*(x::fmpq_poly, y::Integer) = BigInt(y)*x

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^(x::fmpq_poly, y::Int)
   y < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpq_poly_pow, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), 
               &z, &x, y)
   return z
end

###########################################################################################
#
#   Comparisons
#
###########################################################################################

=={S}(x::fmpq_poly{S}, y::fmpq_poly{S}) = ccall((:fmpq_poly_equal, :libflint), Bool, 
                                       (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &x, &y)

###########################################################################################
#
#   Ad hoc comparisons
#
###########################################################################################

function ==(x::fmpq_poly, y::Rational{BigInt}) 
   if length(x) > 1
      return false
   elseif length(x) == 1 
      z = fmpz()
      temp = fmpq_readonly(y)
      ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, 
                       (Ptr{fmpq}, Ptr{fmpq_poly}, Int), &z, &x, 0)
      return ccall((:fmpq_equal, :libflint), Bool, 
               (Ptr{fmpq}, Ptr{fmpq}, Int), &z, &temp, 0)
   else
      return y == 0
   end 
end

==(x::fmpz_poly, y::Integer) = x == BigInt(y)

###########################################################################################
#
#   Truncation
#
###########################################################################################

function truncate(a::fmpq_poly, n::Int)
   n < 0 && throw(DomainError())
   
   if length(a) <= n
      return a
   end

   z = parent(a)()
   ccall((:fmpq_poly_set_trunc, :libflint), Void,
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &a, n)
   return z
end

function mullow{S}(x::fmpq_poly{S}, y::fmpq_poly{S}, n::Int)
   n < 0 && throw(DomainError())
   
   z = parent(x)()
   ccall((:fmpq_poly_mullow, :libflint), Void,
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, &y, n)
   return z
end

###########################################################################################
#
#   Reversal
#
###########################################################################################

function reverse(x::fmpq_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpq_poly_reverse, :libflint), Void,
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, len)
   return z
end

###########################################################################################
#
#   Shifting
#
###########################################################################################

function shift_left(x::fmpq_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpq_poly_shift_left, :libflint), Void,
      (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, len)
   return z
end

function shift_right(x::fmpq_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpq_poly_shift_right, :libflint), Void,
       (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, len)
   return z
end

###########################################################################################
#
#   Euclidean division
#
###########################################################################################

function mod(x::fmpq_poly, y::fmpq_poly)
   y == 0 && throw(DivideError())
   r = parent(x)()
   ccall((:fmpq_poly_rem, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &r, &x, &y)
   return r
end

function divrem(x::fmpq_poly, y::fmpq_poly)
   y == 0 && throw(DivideError())
   q = parent(x)()
   r = parent(x)()
   ccall((:fmpq_poly_divrem, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &q, &r, &x, &y)
   return q, r
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{S}(x::fmpq_poly{S}, y::fmpq_poly{S})
   y == 0 && throw(DivideError())
   z = parent(x)()
   ccall((:fmpq_poly_div, :libflint), Void, 
            (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x, &y)
   return z
end

###########################################################################################
#
#   Ad hoc exact division
#
###########################################################################################

function divexact(x::fmpq_poly, y::BigInt)
   y == 0 && throw(DivideError())
   z = parent(x)()
   temp = fmpz_readonly(y)
   ccall((:fmpq_poly_scalar_div_fmpz, :libflint), Void, 
          (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpz_readonly}), &z, &x, &temp)
   return z
end

function divexact(x::fmpq_poly, y::Rational{BigInt})
   y == 0 && throw(DivideError())
   z = parent(x)()
   temp = fmpq_readonly(y)
   ccall((:fmpq_poly_scalar_div_fmpq, :libflint), Void, 
          (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_readonly}), &z, &x, &temp)
   return z
end

function divexact(x::fmpq_poly, y::Int)
   y == 0 && throw(DivideError())
   z = parent(x)()
   ccall((:fmpq_poly_scalar_div_si, :libflint), Void, 
                        (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, y)
   return z
end

divexact(x::fmpq_poly, y::Integer) = divexact(x, BigInt(y)) 

###########################################################################################
#
#   Content, primitive part, GCD and LCM
#
###########################################################################################

function gcd{S}(x::fmpq_poly{S}, y::fmpq_poly{S})
   z = parent(x)()
   ccall((:fmpq_poly_gcd, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x, &y)
   return z
end

function content(x::fmpq_poly)
   temp = fmpq()
   ccall((:fmpq_poly_content, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq_poly}), &temp, &x)
   return Rational(temp)
end

function primpart(x::fmpq_poly)
   z = parent(x)()
   ccall((:fmpq_poly_primitive_part, :libflint), Void, 
         (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x)
   return z
end

###########################################################################################
#
#   Evaluation
#
###########################################################################################

function evaluate(x::fmpq_poly, y::BigInt)
   z = fmpq()
   temp = fmpz_readonly(y)
   ccall((:fmpq_poly_evaluate_fmpz, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}, Ptr{fmpz}), &z, &x, &temp)
   return Rational(z)
end

function evaluate(x::fmpq_poly, y::Rational{BigInt})
   z = fmpq()
   temp = fmpq_readonly(y)
   ccall((:fmpq_poly_evaluate_fmpq, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}, Ptr{fmpz_readonly}), &z, &x, &temp)
   return Rational(z)
end

evaluate(x::fmpq_poly, y::Integer) = evaluate(x, BigInt(y))

###########################################################################################
#
#   Composition
#
###########################################################################################

function compose{S}(x::fmpq_poly{S}, y::fmpq_poly{S})
   z = parent(x)()
   ccall((:fmpq_poly_compose, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x, &y)
   return z
end

###########################################################################################
#
#   Derivative
#
###########################################################################################

function derivative(x::fmpq_poly)
   z = parent(x)()
   ccall((:fmpq_poly_derivative, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x)
   return z
end

###########################################################################################
#
#   Resultant
#
###########################################################################################

function resultant{S}(x::fmpq_poly{S}, y::fmpq_poly{S})
   temp = fmpq()
   ccall((:fmpq_poly_resultant, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &temp, &x, &y)
   return Rational(temp)
end

###########################################################################################
#
#   Discriminant
#
###########################################################################################

function discriminant(x::fmpq_poly)
   temp = fmpq()
   ccall((:fmpq_poly_discriminant, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}), &temp, &x)
   return Rational(temp)
end

###########################################################################################
#
#   Bezout
#
###########################################################################################

function bezout{S}(x::fmpq_poly{S}, y::fmpq_poly{S})
   z = parent(x)()
   u = parent(x)()
   v = parent(x)()
   ccall((:fmpq_poly_xgcd, :libflint), Void, 
        (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
            &z, &u, &v, &x, &y)
   return (z, u, v)
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function fit!(z::fmpq_poly, n::Int)
   ccall((:fmpq_poly_fit_length, :libflint), Void, 
                    (Ptr{fmpq_poly}, Int), &z, n)
end

function setcoeff!(z::fmpq_poly, n::Int, x::BigInt)
   temp = fmpz_readonly(x)
   ccall((:fmpq_poly_set_coeff_fmpz, :libflint), Void, 
                    (Ptr{fmpq_poly}, Int, Ptr{fmpz}), &z, n, &temp)
end

function setcoeff!(z::fmpq_poly, n::Int, x::Rational{BigInt})
   temp = fmpq_readonly(x)
   ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, 
                    (Ptr{fmpq_poly}, Int, Ptr{fmpq_readonly}), &z, n, &temp)
end

function mul!{S}(z::fmpq_poly{S}, x::fmpq_poly{S}, y::fmpq_poly{S})
   ccall((:fmpq_poly_mul, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x, &y)
end

function addeq!{S}(z::fmpq_poly{S}, x::fmpq_poly{S})
   ccall((:fmpq_poly_add, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &z, &x)
end

###########################################################################################
#
#   Promotions
#
###########################################################################################

Base.promote_rule{S, T <: Integer}(::Type{fmpq_poly{S}}, ::Type{T}) = fmpq_poly{S}

###########################################################################################
#
#   Parent object call overloads
#
###########################################################################################

function Base.call{S}(a::FmpqPolyRing{S})
   z = fmpq_poly{S}()
   z.parent = a
   return z
end

function Base.call{S}(a::FmpqPolyRing{S}, b::Int)
   z = fmpq_poly{S}(b)
   z.parent = a
   return z
end

function Base.call{S}(a::FmpqPolyRing{S}, b::Integer)
   z = fmpq_poly{S}(BigInt(b))
   z.parent = a
   return z
end

function Base.call{S}(a::FmpqPolyRing{S}, b::BigInt)
   z = fmpq_poly{S}(b)
   z.parent = a
   return z
end

function Base.call{S}(a::FmpqPolyRing{S}, b::Rational{BigInt})
   z = fmpq_poly{S}(b)
   z.parent = a
   return z
end

function Base.call{S}(a::FmpqPolyRing{S}, b::Array{Rational{BigInt}, 1})
   z = fmpq_poly{S}(b)
   z.parent = a
   return z
end

Base.call{S}(a::FmpqPolyRing{S}, b::fmpq_poly{S}) = b

###########################################################################################
#
#   PolynomialRing constructor
#
###########################################################################################

function PolynomialRing(R::RationalField, s::String)
   S = symbol(s)

   parent_obj = FmpqPolyRing{S}(R)
   
   return parent_obj, parent_obj([QQ(0), QQ(1)])
end
