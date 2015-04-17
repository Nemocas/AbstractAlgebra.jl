###############################################################################
#
#   fmpq_poly.jl : Flint polynomials over QQ
#
###############################################################################

export fmpq_poly, denominator

###############################################################################
#
#   Data types and memory management
#
###############################################################################

FmpqPolyID = ObjectIdDict()

type FmpqPolyRing <: Ring
   base_ring::Field
   S::Symbol

   function FmpqPolyRing(R::RationalField, s::Symbol)
      return try
         FmpqPolyID[s]
      catch
         FmpqPolyID[s] = new(R, s)
      end
   end
end

type fmpq_poly <: PolyElem
   coeffs::Ptr{Int}
   den::Int # really an fmpz
   alloc::Int
   length::Int
   parent::FmpqPolyRing

   function fmpq_poly()
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::Array{fmpq, 1})
      z = new()
      ccall((:fmpq_poly_init2, :libflint), Void, 
            (Ptr{fmpq_poly}, Int), &z, length(a))
      for i = 1:length(a)
         ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, 
                     (Ptr{fmpq_poly}, Int, Ptr{fmpq}), &z, i - 1, &a[i])
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

   function fmpq_poly(a::fmpz)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_fmpz, :libflint), Void, 
            (Ptr{fmpq_poly}, Ptr{fmpz}), &z, &a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::fmpq)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_fmpq, :libflint), Void,
            (Ptr{fmpq_poly}, Ptr{fmpq}), &z, &a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end

   function fmpq_poly(a::fmpz_poly)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set_fmpz_poly, :libflint), Void, 
            (Ptr{fmpq_poly}, Ptr{fmpz_poly}), &z, &a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end


   function fmpq_poly(a::fmpq_poly)
      z = new()
      ccall((:fmpq_poly_init, :libflint), Void, (Ptr{fmpq_poly},), &z)
      ccall((:fmpq_poly_set, :libflint), Void, 
            (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &a)
      finalizer(z, _fmpq_poly_clear_fn)
      return z
   end
end

function _fmpq_poly_clear_fn(a::fmpq_poly)
   ccall((:fmpq_poly_clear, :libflint), Void, (Ptr{fmpq_poly},), &a)
end

elem_type(::FmpqPolyRing) = fmpq_poly

base_ring(a::FmpqPolyRing) = a.base_ring

function denominator(a::fmpq_poly)
   z = ZZ()
   ccall((:fmpq_poly_get_denominator, :libflint), Void,
         (Ptr{fmpz}, Ptr{fmpq_poly}), &z, &a)
   return z
end

var(a::FmpqPolyRing) = a.S

###############################################################################
#
#   Basic manipulation
#
###############################################################################    
   
length(x::fmpq_poly) = ccall((:fmpq_poly_length, :libflint), Int, 
                                   (Ptr{fmpq_poly},), &x)

function coeff(x::fmpq_poly, n::Int)
   n < 0 && throw(DomainError())
   z = QQ()
   ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, 
               (Ptr{fmpq}, Ptr{fmpq_poly}, Int), &z, &x, n)
   return z
end

zero(a::FmpqPolyRing) = a(0)

one(a::FmpqPolyRing) = a(1)

gen(a::FmpqPolyRing) = a([zero(base_ring(a)), one(base_ring(a))])

isgen(x::fmpq_poly) = ccall((:fmpq_poly_is_x, :libflint), Bool, 
                            (Ptr{fmpq_poly},), &x)

function deepcopy(a::fmpq_poly)
   z = fmpq_poly(a)
   z.parent = parent(a)
   return z
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::fmpq_poly) = canonical_unit(lead(a))

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, x::fmpq_poly)
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fmpq_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
          (Ptr{fmpq_poly}, Ptr{Uint8}), &x, bytestring(string(var(parent(x)))))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
end

function show(io::IO, p::FmpqPolyRing)
   print(io, "Univariate Polynomial Ring in ")
   print(io, string(var(p)))
   print(io, " over ")
   show(io, p.base_ring)
end

show_minus_one(::Type{fmpq_poly}) = show_minus_one(RationalField)

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(x::fmpq_poly)
   z = parent(x)()
   ccall((:fmpq_poly_neg, :libflint), Void, 
         (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x)
   return z
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   z = parent(x)()
   ccall((:fmpq_poly_add, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly},  Ptr{fmpq_poly}), 
               &z, &x, &y)
   return z
end

function -(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   z = parent(x)()
   ccall((:fmpq_poly_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly},  Ptr{fmpq_poly}), 
               &z, &x, &y)
   return z
end

function *(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   z = parent(x)()
   ccall((:fmpq_poly_mul, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly},  Ptr{fmpq_poly}), 
               &z, &x, &y)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::Int, y::fmpq_poly)
   z = parent(y)()
   ccall((:fmpq_poly_scalar_mul_si, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &y, x)
   return z
end

function *(x::fmpz, y::fmpq_poly)
   z = parent(y)()
   ccall((:fmpq_poly_scalar_mul_fmpz, :libflint), Void, 
         (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpz}), &z, &y, &x)
   return z
end

function *(x::fmpq, y::fmpq_poly)
   z = parent(y)()
   ccall((:fmpq_poly_scalar_mul_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq}), &z, &y, &x)
   return z
end

function +(x::fmpq_poly, y::Int)
   z = parent(x)()
   ccall((:fmpq_poly_add_si, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, y)
   return z
end

function +(x::fmpq_poly, y::fmpz)
   z = parent(x)()
   ccall((:fmpq_poly_add_fmpz, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpz}), &z, &x, &y)
   return z
end

function +(x::fmpq_poly, y::fmpq)
   z = parent(x)()
   ccall((:fmpq_poly_add_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq}), &z, &x, &y)
   return z
end

function -(x::fmpq_poly, y::Int)
   z = parent(x)()
   ccall((:fmpq_poly_sub_si, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, y)
   return z
end

function -(x::fmpq_poly, y::fmpz)
   z = parent(x)()
   ccall((:fmpq_poly_sub_fmpz, :libflint), Void, 
         (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpz}), &z, &x, &y)
   return z
end

function -(x::fmpq_poly, y::fmpq)
   z = parent(x)()
   ccall((:fmpq_poly_sub_fmpq, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq}), &z, &x, &y)
   return z
end

function -(x::Int, y::fmpq_poly)
   z = parent(y)()
   ccall((:fmpq_poly_si_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Int, Ptr{fmpq_poly}), &z, x, &y)
   return z
end

function -(x::fmpz, y::fmpq_poly)
   z = parent(y)()
   ccall((:fmpq_poly_fmpz_sub, :libflint), Void, 
         (Ptr{fmpq_poly}, Ptr{fmpz}, Ptr{fmpq_poly}), &z, &x, &y)
   return z
end

function -(x::fmpq, y::fmpq_poly)
   z = parent(y)()
   ccall((:fmpq_poly_fmpq_sub, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq}, Ptr{fmpq_poly}), &z, &x, &y)
   return z
end

+(x::Int, y::fmpq_poly) = y + x

+(x::fmpz, y::fmpq_poly) = y + x

+(x::fmpq, y::fmpq_poly) = y + x

*(x::fmpq_poly, y::Int) = y*x

*(x::fmpq_poly, y::fmpz) = y*x

*(x::fmpq_poly, y::fmpq) = y*x

+(x::Integer, y::fmpq_poly) = y + ZZ(x)

-(x::Integer, y::fmpq_poly) = ZZ(x) - y

*(x::Integer, y::fmpq_poly) = ZZ(x)*y

+(x::fmpq_poly, y::Integer) = x + ZZ(y)

-(x::fmpq_poly, y::Integer) = x - ZZ(y)

*(x::fmpq_poly, y::Integer) = ZZ(y)*x

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::fmpq_poly, y::Int)
   y < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpq_poly_pow, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), 
               &z, &x, y)
   return z
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function ==(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   return ccall((:fmpq_poly_equal, :libflint), Bool, 
                                      (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &x, &y)
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function ==(x::fmpq_poly, y::fmpq) 
   if length(x) > 1
      return false
   elseif length(x) == 1 
      z = QQ()
      ccall((:fmpq_poly_get_coeff_fmpq, :libflint), Void, 
                       (Ptr{fmpq}, Ptr{fmpq_poly}, Int), &z, &x, 0)
      return ccall((:fmpq_equal, :libflint), Bool, 
               (Ptr{fmpq}, Ptr{fmpq}, Int), &z, &y, 0)
   else
      return y == 0
   end 
end

==(x::fmpz_poly, y::Integer) = x == ZZ(y)

###############################################################################
#
#   Truncation
#
###############################################################################

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

function mullow(x::fmpq_poly, y::fmpq_poly, n::Int)
   check_parent(x, y)
   n < 0 && throw(DomainError())
   
   z = parent(x)()
   ccall((:fmpq_poly_mullow, :libflint), Void,
         (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, &y, n)
   return z
end

###############################################################################
#
#   Reversal
#
###############################################################################

function reverse(x::fmpq_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpq_poly_reverse, :libflint), Void,
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, len)
   return z
end

###############################################################################
#
#   Shifting
#
###############################################################################

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

###############################################################################
#
#   Euclidean division
#
###############################################################################

function mod(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   y == 0 && throw(DivideError())
   r = parent(x)()
   ccall((:fmpq_poly_rem, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &r, &x, &y)
   return r
end

function divrem(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   y == 0 && throw(DivideError())
   q = parent(x)()
   r = parent(x)()
   ccall((:fmpq_poly_divrem, :libflint), Void, 
         (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), 
               &q, &r, &x, &y)
   return q, r
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   y == 0 && throw(DivideError())
   z = parent(x)()
   ccall((:fmpq_poly_div, :libflint), Void, 
            (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x, &y)
   return z
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::fmpq_poly, y::fmpz)
   y == 0 && throw(DivideError())
   z = parent(x)()
   ccall((:fmpq_poly_scalar_div_fmpz, :libflint), Void, 
          (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpz}), &z, &x, &y)
   return z
end

function divexact(x::fmpq_poly, y::fmpq)
   y == 0 && throw(DivideError())
   z = parent(x)()
   ccall((:fmpq_poly_scalar_div_fmpq, :libflint), Void, 
          (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq}), &z, &x, &y)
   return z
end

function divexact(x::fmpq_poly, y::Int)
   y == 0 && throw(DivideError())
   z = parent(x)()
   ccall((:fmpq_poly_scalar_div_si, :libflint), Void, 
                        (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Int), &z, &x, y)
   return z
end

divexact(x::fmpq_poly, y::Integer) = divexact(x, ZZ(y)) 

###############################################################################
#
#   Content, primitive part, GCD and LCM
#
###############################################################################

function gcd(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   z = parent(x)()
   ccall((:fmpq_poly_gcd, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x, &y)
   return z
end

function content(x::fmpq_poly)
   z = QQ()
   ccall((:fmpq_poly_content, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq_poly}), &z, &x)
   return z
end

function primpart(x::fmpq_poly)
   z = parent(x)()
   ccall((:fmpq_poly_primitive_part, :libflint), Void, 
         (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x)
   return z
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(x::fmpq_poly, y::fmpz)
   z = QQ()
   ccall((:fmpq_poly_evaluate_fmpz, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}, Ptr{fmpz}), &z, &x, &y)
   return z
end

function evaluate(x::fmpq_poly, y::fmpq)
   z = QQ()
   ccall((:fmpq_poly_evaluate_fmpq, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}, Ptr{fmpq}), &z, &x, &y)
   return z
end

evaluate(x::fmpq_poly, y::Integer) = evaluate(x, ZZ(y))

###############################################################################
#
#   Composition
#
###############################################################################

function compose(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   z = parent(x)()
   ccall((:fmpq_poly_compose, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x, &y)
   return z
end

###############################################################################
#
#   Derivative
#
###############################################################################

function derivative(x::fmpq_poly)
   z = parent(x)()
   ccall((:fmpq_poly_derivative, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x)
   return z
end

###############################################################################
#
#   Resultant
#
###############################################################################

function resultant(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   z = QQ()
   ccall((:fmpq_poly_resultant, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x, &y)
   return z
end

###############################################################################
#
#   Discriminant
#
###############################################################################

function discriminant(x::fmpq_poly)
   z = QQ()
   ccall((:fmpq_poly_discriminant, :libflint), Void, 
                (Ptr{fmpq}, Ptr{fmpq_poly}), &z, &x)
   return z
end

###############################################################################
#
#   GCDX
#
###############################################################################

function gcdx(x::fmpq_poly, y::fmpq_poly)
   check_parent(x, y)
   z = parent(x)()
   u = parent(x)()
   v = parent(x)()
   ccall((:fmpq_poly_xgcd, :libflint), Void, 
        (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}, 
                                     Ptr{fmpq_poly}), &z, &u, &v, &x, &y)
   return (z, u, v)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function fit!(z::fmpq_poly, n::Int)
   ccall((:fmpq_poly_fit_length, :libflint), Void, 
                    (Ptr{fmpq_poly}, Int), &z, n)
end

function setcoeff!(z::fmpq_poly, n::Int, x::fmpz)
   ccall((:fmpq_poly_set_coeff_fmpz, :libflint), Void, 
                    (Ptr{fmpq_poly}, Int, Ptr{fmpz}), &z, n, &x)
end

function setcoeff!(z::fmpq_poly, n::Int, x::fmpq)
   ccall((:fmpq_poly_set_coeff_fmpq, :libflint), Void, 
                    (Ptr{fmpq_poly}, Int, Ptr{fmpq}), &z, n, &x)
end

function mul!(z::fmpq_poly, x::fmpq_poly, y::fmpq_poly)
   ccall((:fmpq_poly_mul, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &x, &y)
end

function addeq!(z::fmpq_poly, x::fmpq_poly)
   ccall((:fmpq_poly_add, :libflint), Void, 
                (Ptr{fmpq_poly}, Ptr{fmpq_poly}, Ptr{fmpq_poly}), &z, &z, &x)
end

###############################################################################
#
#   Promotions
#
###############################################################################

Base.promote_rule{T <: Integer}(::Type{fmpq_poly}, ::Type{T}) = fmpq_poly

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function Base.call(a::FmpqPolyRing)
   z = fmpq_poly()
   z.parent = a
   return z
end

function Base.call(a::FmpqPolyRing, b::Int)
   z = fmpq_poly(b)
   z.parent = a
   return z
end

function Base.call(a::FmpqPolyRing, b::Integer)
   z = fmpq_poly(ZZ(b))
   z.parent = a
   return z
end

function Base.call(a::FmpqPolyRing, b::fmpz)
   z = fmpq_poly(b)
   z.parent = a
   return z
end

function Base.call(a::FmpqPolyRing, b::fmpq)
   z = fmpq_poly(b)
   z.parent = a
   return z
end

function Base.call(a::FmpqPolyRing, b::Array{fmpq, 1})
   z = fmpq_poly(b)
   z.parent = a
   return z
end

Base.call(a::FmpqPolyRing, b::fmpq_poly) = b

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

function PolynomialRing(R::RationalField, s::String)
   S = symbol(s)

   parent_obj = FmpqPolyRing(R, S)
   
   return parent_obj, parent_obj([QQ(0), QQ(1)])
end
