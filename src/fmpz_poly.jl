###############################################################################
#
#   fmpz_poly.jl : Flint polynomials over ZZ
#
###############################################################################

## FIXME : add clear_readonly and fmpz initialiser/finalizer and use instead of 
##         init/clear
## FIXME : make function names conform to Julia standard
## FIXME : put special polynomials back in
## FIXME : don't use isequal; in Julia it's for objects that hash to the same 
##         value
## FIXME : rename primpart to primitive_part
## FIXME : figure out why length{S}(x::fmpq_poly{S}) requires the {S} when 
##         called from Base.call{S}(a::NfNumberField{S}, pol::fmpq_poly) in 
##         nf.jl
## FIXME : fix needs_parentheses and is_negative in nf.jl
## FIXME : add hashing for all types
## FIXME : canonical_unit for fractions is odd for (1//(x^2+1))//(2//(x+1)) 
##         over rationals
## FIXME : should Fraction only use canonical_unit when printing?

export fmpz_poly

###############################################################################
#
#   Data types and memory management
#
###############################################################################

FmpzPolyID = ObjectIdDict()

type FmpzPolyRing <: Ring
   base_ring::Ring
   S::Symbol

   function FmpzPolyRing(s::Symbol)
      return try
         FmpzPolyID[s]
      catch
         FmpzPolyID[s] = new(ZZ, s)
      end
   end
end

type fmpz_poly <: PolyElem
   coeffs::Ptr{Void}
   alloc::Int
   length::Int
   parent::FmpzPolyRing

   function fmpz_poly()
      z = new()
      ccall((:fmpz_poly_init, :libflint), Void, (Ptr{fmpz_poly},), &z)
      finalizer(z, _fmpz_poly_clear_fn)
      return z
   end

   function fmpz_poly(a::Array{fmpz, 1})
      z = new()
      ccall((:fmpz_poly_init2, :libflint), Void, 
            (Ptr{fmpz_poly}, Int), &z, length(a))
      for i = 1:length(a)
         ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Void, 
                     (Ptr{fmpz_poly}, Int, Ptr{fmpz}), &z, i - 1, &a[i])
      end
      finalizer(z, _fmpz_poly_clear_fn)
      return z
   end

   function fmpz_poly(a::Int)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Void, (Ptr{fmpz_poly},), &z)
      ccall((:fmpz_poly_set_si, :libflint), Void, (Ptr{fmpz_poly}, Int), &z, a)
      finalizer(z, _fmpz_poly_clear_fn)
      return z
   end

   function fmpz_poly(a::fmpz)
      z = new()
      ccall((:fmpz_poly_init, :libflint), Void, (Ptr{fmpz_poly},), &z)
      ccall((:fmpz_poly_set_mpz, :libflint), Void, 
            (Ptr{fmpz_poly}, Ptr{fmpz}), &z, &a)
      finalizer(z, _fmpz_poly_clear_fn)
      return z
   end
end

function _fmpz_poly_clear_fn(a::fmpz_poly)
   ccall((:fmpz_poly_clear, :libflint), Void, (Ptr{fmpz_poly},), &a)
end

elem_type(::FmpzPolyRing) = fmpz_poly

base_ring(a::FmpzPolyRing) = a.base_ring

parent(a::fmpz_poly) = a.parent

var(a::FmpzPolyRing) = a.S

###############################################################################
#
#   Basic manipulation
#
###############################################################################   
   
length(x::fmpz_poly) = ccall((:fmpz_poly_length, :libflint), Int, 
                             (Ptr{fmpz_poly},), &x)

function coeff(x::fmpz_poly, n::Int)
   n < 0 && throw(DomainError())
   temp = ZZ()
   ccall((:fmpz_poly_get_coeff_fmpz, :libflint), Void, 
               (Ptr{fmpz}, Ptr{fmpz_poly}, Int), &temp, &x, n)
   return fmpz(temp)
end

zero(a::FmpzPolyRing) = a(0)

one(a::FmpzPolyRing) = a(1)

gen(a::FmpzPolyRing) = a([zero(base_ring(a)), one(base_ring(a))])

isgen(x::fmpz_poly) = ccall((:fmpz_poly_is_x, :libflint), Bool, 
                            (Ptr{fmpz_poly},), &x)

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::fmpz_poly) = canonical_unit(lead(a))

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, x::fmpz_poly)
   if length(x) == 0
      print(io, "0")
   else
      cstr = ccall((:fmpz_poly_get_str_pretty, :libflint), Ptr{Uint8}, 
          (Ptr{fmpz_poly}, Ptr{Uint8}), &x, bytestring(string(var(parent(x)))))

      print(io, bytestring(cstr))

      ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
   end
end

function show(io::IO, p::FmpzPolyRing)
   print(io, "Univariate Polynomial Ring in ")
   print(io, string(var(p)))
   print(io, " over ")
   show(io, p.base_ring)
end

show_minus_one(::Type{fmpz_poly}) = show_minus_one(fmpz)

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(x::fmpz_poly)
   z = parent(x)()
   ccall((:fmpz_poly_neg, :libflint), Void, 
         (Ptr{fmpz_poly}, Ptr{fmpz_poly}), &z, &x)
   return z
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(x::fmpz_poly, y::fmpz_poly)
   check_parent(x, y)
   z = parent(x)()
   ccall((:fmpz_poly_add, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly},  Ptr{fmpz_poly}), 
               &z, &x, &y)
   return z
end

function -(x::fmpz_poly, y::fmpz_poly)
   check_parent(x, y)
   z = parent(x)()
   ccall((:fmpz_poly_sub, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly},  Ptr{fmpz_poly}), 
               &z, &x, &y)
   return z
end

function *(x::fmpz_poly, y::fmpz_poly)
   check_parent(x, y)
   z = parent(x)()
   ccall((:fmpz_poly_mul, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly},  Ptr{fmpz_poly}), 
               &z, &x, &y)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::Int, y::fmpz_poly)
   z = parent(y)()
   ccall((:fmpz_poly_scalar_mul_si, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), &z, &y, x)
   return z
end

function *(x::fmpz, y::fmpz_poly)
   z = parent(y)()
   ccall((:fmpz_poly_scalar_mul_fmpz, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz}), &z, &y, &x)
   return z
end

function +(x::fmpz_poly, y::Int)
   z = parent(x)()
   ccall((:fmpz_poly_add_si, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), &z, &x, y)
   return z
end

function +(x::fmpz_poly, y::fmpz)
   z = parent(x)()
   ccall((:fmpz_poly_add_fmpz, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz}), &z, &x, &y)
   return z
end

function -(x::fmpz_poly, y::Int)
   z = parent(x)()
   ccall((:fmpz_poly_sub_si, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), &z, &x, y)
   return z
end

function -(x::fmpz_poly, y::fmpz)
   z = parent(x)()
   ccall((:fmpz_poly_sub_fmpz, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz}), &z, &x, &y)
   return z
end

function -(x::Int, y::fmpz_poly)
   z = parent(y)()
   ccall((:fmpz_poly_si_sub, :libflint), Void, 
                (Ptr{fmpz_poly}, Int, Ptr{fmpz_poly}), &z, x, &y)
   return z
end

function -(x::fmpz, y::fmpz_poly)
   z = parent(y)()
   ccall((:fmpz_poly_fmpz_sub, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz}, Ptr{fmpz_poly}), &z, &x, &y)
   return z
end

+(x::Int, y::fmpz_poly) = y + x

+(x::fmpz, y::fmpz_poly) = y + x

*(x::fmpz_poly, y::Int) = y*x

*(x::fmpz_poly, y::fmpz) = y*x

+(x::Integer, y::fmpz_poly) = y + fmpz(x)

-(x::Integer, y::fmpz_poly) = fmpz(x) - y

*(x::Integer, y::fmpz_poly) = fmpz(x)*y

+(x::fmpz_poly, y::Integer) = x + fmpz(y)

-(x::fmpz_poly, y::Integer) = x - fmpz(y)

*(x::fmpz_poly, y::Integer) = fmpz(y)*x

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::fmpz_poly, y::Int)
   y < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpz_poly_pow, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), 
               &z, &x, y)
   return z
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function ==(x::fmpz_poly, y::fmpz_poly)
   check_parent(x, y)
   return ccall((:fmpz_poly_equal, :libflint), Bool, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}), &x, &y)
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function ==(x::fmpz_poly, y::fmpz) 
   if length(x) > 1
      return false
   elseif length(x) == 1 
      z = ZZ()
      ccall((:fmpz_poly_get_coeff_fmpz, :libflint), Void, 
                       (Ptr{fmpz}, Ptr{fmpz_poly}, Int), &z, &x, 0)
      return ccall((:fmpz_equal, :libflint), Bool, 
               (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &y, 0)
   else
      return y == 0
   end 
end

==(x::fmpz_poly, y::Integer) = x == fmpz(y)

###############################################################################
#
#   Truncation
#
###############################################################################

function truncate(a::fmpz_poly, n::Int)
   n < 0 && throw(DomainError())
   
   if length(a) <= n
      return a
   end

   z = parent(a)()
   ccall((:fmpz_poly_set_trunc, :libflint), Void,
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), &z, &a, n)
   return z
end

function mullow(x::fmpz_poly, y::fmpz_poly, n::Int)
   check_parent(x, y)
   n < 0 && throw(DomainError())
   
   z = parent(x)()
   ccall((:fmpz_poly_mullow, :libflint), Void,
         (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), &z, &x, &y, n)
   return z
end

###############################################################################
#
#   Reversal
#
###############################################################################

function reverse(x::fmpz_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpz_poly_reverse, :libflint), Void,
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), &z, &x, len)
   return z
end

###############################################################################
#
#   Shifting
#
###############################################################################

function shift_left(x::fmpz_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpz_poly_shift_left, :libflint), Void,
      (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), &z, &x, len)
   return z
end

function shift_right(x::fmpz_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:fmpz_poly_shift_right, :libflint), Void,
       (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), &z, &x, len)
   return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::fmpz_poly, y::fmpz_poly)
   check_parent(x, y)
   y == 0 && throw(DivideError())
   z = parent(x)()
   ccall((:fmpz_poly_div, :libflint), Void, 
            (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), &z, &x, &y)
   return z
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::fmpz_poly, y::fmpz)
   y == 0 && throw(DivideError())
   z = parent(x)()
   ccall((:fmpz_poly_scalar_divexact_fmpz, :libflint), Void, 
          (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz}), &z, &x, &y)
   return z
end

function divexact(x::fmpz_poly, y::Int)
   y == 0 && throw(DivideError())
   z = parent(x)()
   ccall((:fmpz_poly_scalar_divexact_si, :libflint), Void, 
                        (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Int), &z, &x, y)
   return z
end

divexact(x::fmpz_poly, y::Integer) = divexact(x, fmpz(y)) 

###############################################################################
#
#   Pseudodivision
#
###############################################################################

function pseudorem(x::fmpz_poly, y::fmpz_poly)
   check_parent(x, y)
   y == 0 && throw(DivideError())
   diff = length(x) - length(y)
   r = parent(x)()
   d = Array(Int, 1)
   ccall((:fmpz_poly_pseudo_rem, :libflint), Void, 
     (Ptr{fmpz_poly}, Ptr{Int}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), &r, d, &x, &y)
   if (diff > d[1])
      return lead(y)^(diff - d[1])*r
   else
      return r
   end
end

function pseudodivrem(x::fmpz_poly, y::fmpz_poly)
   check_parent(x, y)
   y == 0 && throw(DivideError())
   diff = length(x) - length(y)
   q = parent(x)()
   r = parent(x)()
   d = Array(Int, 1)
   ccall((:fmpz_poly_pseudo_divrem_divconquer, :libflint), Void, 
    (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{Int}, Ptr{fmpz_poly}, Ptr{fmpz_poly}),
               &q, &r, d, &x, &y)
   if (diff > d[1])
      m = lead(y)^(diff - d[1])
      return m*q, m*r
   else
      return q, r
   end
end

###############################################################################
#
#   Content, primitive part, GCD and LCM
#
###############################################################################

function gcd(x::fmpz_poly, y::fmpz_poly)
   check_parent(x, y)
   z = parent(x)()
   ccall((:fmpz_poly_gcd, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), &z, &x, &y)
   return z
end

function content(x::fmpz_poly)
   temp = ZZ()
   ccall((:fmpz_poly_content, :libflint), Void,
         (Ptr{fmpz}, Ptr{fmpz_poly}), &temp, &x)
   return fmpz(temp)
end

function primpart(x::fmpz_poly)
   z = parent(x)()
   ccall((:fmpz_poly_primitive_part, :libflint), Void, 
         (Ptr{fmpz_poly}, Ptr{fmpz_poly}), &z, &x)
   return z
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(x::fmpz_poly, y::fmpz)
   z = ZZ()
   ccall((:fmpz_poly_evaluate_fmpz, :libflint), Void, 
        (Ptr{fmpz}, Ptr{fmpz_poly}, Ptr{fmpz}), &z, &x, &y)
   return fmpz(z)
end

evaluate(x::fmpz_poly, y::Integer) = evaluate(x, fmpz(y))

###############################################################################
#
#   Composition
#
###############################################################################

function compose(x::fmpz_poly, y::fmpz_poly)
   check_parent(x, y)
   z = parent(x)()
   ccall((:fmpz_poly_compose, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), &z, &x, &y)
   return z
end

###############################################################################
#
#   Derivative
#
###############################################################################

function derivative(x::fmpz_poly)
   z = parent(x)()
   ccall((:fmpz_poly_derivative, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}), &z, &x)
   return z
end

###############################################################################
#
#   Resultant
#
###############################################################################

function resultant(x::fmpz_poly, y::fmpz_poly)
   check_parent(x, y)
   temp = ZZ()
   ccall((:fmpz_poly_resultant, :libflint), Void, 
                (Ptr{fmpz}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), &temp, &x, &y)
   return fmpz(temp)
end

###############################################################################
#
#   Discriminant
#
###############################################################################

function discriminant(x::fmpz_poly)
   temp = ZZ()
   ccall((:fmpz_poly_discriminant, :libflint), Void, 
                (Ptr{fmpz}, Ptr{fmpz_poly}), &temp, &x)
   return fmpz(temp)
end

###############################################################################
#
#   GCDX
#
###############################################################################

function gcdx(a::fmpz_poly, b::fmpz_poly)
   check_parent(x, y)
   lena = length(a)
   lenb = length(b)
   (lena <= 1 || lenb <= 1) && error("Constant polynomial in gcdx")  
   temp = ZZ()
   u = parent(a)()
   v = parent(a)()
   c1 = content(a)
   c2 = content(b)
   x = divexact(a, c1)
   y = divexact(b, c2)
   ccall((:fmpz_poly_xgcd_modular, :libflint), Void, 
   (Ptr{fmpz}, Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), 
            &temp, &u, &v, &x, &y)
   z = fmpz(temp)*c1^(lenb - 1)*c2^(lena - 1)
   u *= c1^(lenb - 2)*c2^(lena - 1)
   v *= c1^(lenb - 1)*c2^(lena - 2)   
   return (z, u, v)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function fit!(z::fmpz_poly, n::Int)
   ccall((:fmpz_poly_fit_length, :libflint), Void, 
                    (Ptr{fmpz_poly}, Int), &z, n)
end

function setcoeff!(z::fmpz_poly, n::Int, x::fmpz)
   ccall((:fmpz_poly_set_coeff_fmpz, :libflint), Void, 
                    (Ptr{fmpz_poly}, Int, Ptr{fmpz}), &z, n, &x)
end

function mul!(z::fmpz_poly, x::fmpz_poly, y::fmpz_poly)
   ccall((:fmpz_poly_mul, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), &z, &x, &y)
end

function addeq!(z::fmpz_poly, x::fmpz_poly)
   ccall((:fmpz_poly_add, :libflint), Void, 
                (Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz_poly}), &z, &z, &x)
end

###############################################################################
#
#   Promotions
#
###############################################################################

Base.promote_rule{T <: Integer}(::Type{fmpz_poly}, ::Type{T}) = fmpz_poly

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function Base.call(a::FmpzPolyRing)
   z = fmpz_poly()
   z.parent = a
   return z
end

function Base.call(a::FmpzPolyRing, b::Int)
   z = fmpz_poly(b)
   z.parent = a
   return z
end

function Base.call(a::FmpzPolyRing, b::Integer)
   z = fmpz_poly(fmpz(b))
   z.parent = a
   return z
end

function Base.call(a::FmpzPolyRing, b::fmpz)
   z = fmpz_poly(b)
   z.parent = a
   return z
end

function Base.call(a::FmpzPolyRing, b::Array{fmpz, 1})
   z = fmpz_poly(b)
   z.parent = a
   return z
end

Base.call(a::FmpzPolyRing, b::fmpz_poly) = b

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

function PolynomialRing(R::IntegerRing, s::String)
   S = symbol(s)

   parent_obj = FmpzPolyRing(S)
   
   return parent_obj, parent_obj([ZZ(0), ZZ(1)])
end
