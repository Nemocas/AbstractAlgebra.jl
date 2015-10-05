###############################################################################
#
#   acb_poly.jl : Polynomials over arb
#
###############################################################################

export AcbPolyRing, acb_poly, strongequal, derivative, integral, evaluate,
       evaluate2, compose, from_roots, evaluate_iter, evaluate_fast, evaluate,
       interpolate_newton, interpolate_barycentric, interpolate_fast, interpolate

###############################################################################
#
#   Basic manipulation
#
###############################################################################    
  
length(x::acb_poly) = ccall((:acb_poly_length, :libarb), Int, 
                                   (Ptr{acb_poly},), &x)

set_length!(x::acb_poly, n::Int) = ccall((:_acb_poly_set_length, :libarb), Void,
                                   (Ptr{acb_poly}, Int), &x, n)

degree(x::acb_poly) = length(x) - 1

function coeff(a::acb_poly, n::Int)
  n < 0 && throw(DomainError())
  t = parent(a).base_ring()
  ccall((:acb_poly_get_coeff_acb, :libarb), Void,
              (Ptr{acb}, Ptr{acb_poly}, Int), &t, &a, n)
  return t
end

zero(a::AcbPolyRing) = a(0)

one(a::AcbPolyRing) = a(1)

function gen(a::AcbPolyRing)
   z = acb_poly()
   ccall((:acb_poly_set_coeff_si, :libarb), Void,
        (Ptr{acb_poly}, Int, Int), &z, 1, 1)
   z.parent = a
   return z
end

# todo: write a C function for this
function isgen(a::acb_poly)
   return strongequal(a, gen(parent(a)))
end

#function iszero(a::acb_poly)
#   return length(a) == 0
#end

#function isone(a::acb_poly)
#   return strongequal(a, one(parent(a)))
#end

function deepcopy(a::acb_poly)
   z = acb_poly(a)
   z.parent = parent(a)
   return z
end

###############################################################################
#
#   AbstractString{} I/O
#
###############################################################################

function show(io::IO, x::AcbPolyRing)
  print(io, "Univariate Polynomial Ring in ")
  print(io, string(var(x)))
  print(io, " over ")
  show(io, x.base_ring)
end

function show(io::IO, f::acb_poly)
  if length(f) == 0
    print(io, "0")
  else
    print(io, "[ ")
    for i in 0:degree(f)-1
      show(io, coeff(f,i))
      print(io, ", ")
    end
    show(coeff(f,degree(f)))
    print(io, " ]")
  end
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function strongequal(x::acb_poly, y::acb_poly)
   return ccall((:acb_poly_equal, :libarb), Bool, 
                                      (Ptr{acb_poly}, Ptr{acb_poly}), &x, &y)
end

function overlaps(x::acb_poly, y::acb_poly)
   return ccall((:acb_poly_overlaps, :libarb), Bool, 
                                      (Ptr{acb_poly}, Ptr{acb_poly}), &x, &y)
end

function contains(x::acb_poly, y::acb_poly)
   return ccall((:acb_poly_contains, :libarb), Bool, 
                                      (Ptr{acb_poly}, Ptr{acb_poly}), &x, &y)
end

function contains(x::acb_poly, y::fmpz_poly)
   return ccall((:acb_poly_contains_fmpz_poly, :libarb), Bool, 
                                      (Ptr{acb_poly}, Ptr{fmpz_poly}), &x, &y)
end

function contains(x::acb_poly, y::fmpq_poly)
   return ccall((:acb_poly_contains_fmpq_poly, :libarb), Bool, 
                                      (Ptr{acb_poly}, Ptr{fmpq_poly}), &x, &y)
end

function ==(x::acb_poly, y::acb_poly)
    if length(x) != length(y)
        return false
    end
    for i = 0:degree(x)
        if !(coeff(x, i) == coeff(y, i))
            return false
        end
    end
    return true
end

function !=(x::acb_poly, y::acb_poly)
    for i = 0:max(degree(x), degree(y))
        if coeff(x, i) != coeff(y, i)
            return true
        end
    end
    return false
end

function unique_integer(x::acb_poly)
  z = FmpzPolyRing(var(parent(x)))()
  unique = ccall((:acb_poly_get_unique_fmpz_poly, :libarb), Int,
    (Ptr{fmpz_poly}, Ptr{acb_poly}), &z, &x)
  return (unique != 0, z)
end

###############################################################################
#
#   Shifting
#
###############################################################################

function shift_left(x::acb_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:acb_poly_shift_left, :libarb), Void,
      (Ptr{acb_poly}, Ptr{acb_poly}, Int), &z, &x, len)
   return z
end

function shift_right(x::acb_poly, len::Int)
   len < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:acb_poly_shift_right, :libarb), Void,
       (Ptr{acb_poly}, Ptr{acb_poly}, Int), &z, &x, len)
   return z
end

################################################################################
#
#  Unary operations
#
################################################################################

function -(x::acb_poly)
  z = parent(x)()
  ccall((:acb_poly_neg, :libarb), Void, (Ptr{acb_poly}, Ptr{acb_poly}), &z, &x)
  return z
end

################################################################################
#
#  Binary operations
#
################################################################################

function +(x::acb_poly, y::acb_poly)
  z = parent(x)()
  ccall((:acb_poly_add, :libarb), Void,
              (Ptr{acb_poly}, Ptr{acb_poly}, Ptr{acb_poly}, Int),
              &z, &x, &y, prec(parent(x)))
  return z
end

function *(x::acb_poly, y::acb_poly)
  z = parent(x)()
  ccall((:acb_poly_mul, :libarb), Void,
              (Ptr{acb_poly}, Ptr{acb_poly}, Ptr{acb_poly}, Int),
              &z, &x, &y, prec(parent(x)))
  return z
end

function -(x::acb_poly, y::acb_poly)
  z = parent(x)()
  ccall((:acb_poly_sub, :libarb), Void,
              (Ptr{acb_poly}, Ptr{acb_poly}, Ptr{acb_poly}, Int),
              &z, &x, &y, prec(parent(x)))
  return z
end

function ^(x::acb_poly, y::Int)
  y < 0 && throw(DomainError())
  z = parent(x)()
  ccall((:acb_poly_pow_ui, :libarb), Void,
              (Ptr{acb_poly}, Ptr{acb_poly}, UInt, Int),
              &z, &x, y, prec(parent(x)))
  return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

#function +(x::acb_poly, y::Union{Int,fmpz,fmpq,Float64,arb,fmpz_poly,fmpq_poly})
#    return x + parent(x)(y)
#end

#function -(x::acb_poly, y::Union{Int,fmpz,fmpq,Float64,arb,fmpz_poly,fmpq_poly})
#    return x - parent(x)(y)
#end

#function *(x::acb_poly, y::Union{Int,fmpz,fmpq,Float64,arb,fmpz_poly,fmpq_poly})
#    return x * parent(x)(y)
#end

#function +(x::Union{Int,fmpz,fmpq,Float64,arb,fmpz_poly,fmpq_poly}, y::acb_poly)
#    return parent(y)(x) + y
#end

#function -(x::Union{Int,fmpz,fmpq,Float64,arb,fmpz_poly,fmpq_poly}, y::acb_poly)
#    return parent(y)(x) - y
#end

#function *(x::Union{Int,fmpz,fmpq,Float64,arb,fmpz_poly,fmpq_poly}, y::acb_poly)
#    return parent(y)(x) * y
#end

###############################################################################
#
#   Scalar division
#
###############################################################################

function divexact(x::acb_poly, y::Union{Int,fmpz,fmpq,Float64,arb})
    return x * inv(base_ring(parent(x))(y))
end

//(x::acb_poly, y::Union{Int,fmpz,fmpq,Float64,arb}) = divexact(x, y)

###############################################################################
#
#   Euclidean division
#
###############################################################################

function divrem(x::acb_poly, y::acb_poly)
   iszero(y) && throw(DivideError())
   q = parent(x)()
   r = parent(x)()
   if (ccall((:acb_poly_divrem, :libarb), Int, 
         (Ptr{acb_poly}, Ptr{acb_poly}, Ptr{acb_poly}, Ptr{acb_poly}, Int), 
               &q, &r, &x, &y, prec(parent(x))) == 1)
      return (q, r)
   else
      throw(DivideError())
   end
end

function mod(x::acb_poly, y::acb_poly)
   return divrem(x, y)[2]
end

function divexact(x::acb_poly, y::acb_poly)
   return divrem(x, y)[1]
end

###############################################################################
#
#   Truncation
#
###############################################################################

function truncate(a::acb_poly, n::Int)
   n < 0 && throw(DomainError())
   if length(a) <= n
      return a
   end
   # todo: implement set_trunc in arb
   z = deepcopy(a)
   ccall((:acb_poly_truncate, :libarb), Void,
                (Ptr{acb_poly}, Int), &z, n)
   return z
end

function mullow(x::acb_poly, y::acb_poly, n::Int)
   n < 0 && throw(DomainError())
   z = parent(x)()
   ccall((:acb_poly_mullow, :libarb), Void,
         (Ptr{acb_poly}, Ptr{acb_poly}, Ptr{acb_poly}, Int, Int),
            &z, &x, &y, n, prec(parent(x)))
   return z
end

###############################################################################
#
#   Reversal
#
###############################################################################

#function reverse(x::acb_poly, len::Int)
#   len < 0 && throw(DomainError())
#   z = parent(x)()
#   ccall((:acb_poly_reverse, :libarb), Void,
#                (Ptr{acb_poly}, Ptr{acb_poly}, Int), &z, &x, len)
#   return z
#end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(x::acb_poly, y::acb)
   z = parent(y)()
   ccall((:acb_poly_evaluate, :libarb), Void, 
                (Ptr{acb}, Ptr{acb_poly}, Ptr{acb}, Int),
                &z, &x, &y, prec(parent(y)))
   return z
end

function evaluate2(x::acb_poly, y::acb)
   z = parent(y)()
   w = parent(y)()
   ccall((:acb_poly_evaluate2, :libarb), Void, 
                (Ptr{acb}, Ptr{acb}, Ptr{acb_poly}, Ptr{acb}, Int),
                &z, &w, &x, &y, prec(parent(y)))
   return z, w
end

function evaluate(x::acb_poly, y::Union{Int,Float64,fmpq,arb})
    return evaluate(x, base_ring(parent(x))(y))
end

function evaluate(x::acb_poly, y::fmpz)
    return evaluate(x, base_ring(parent(x))(y))
end

function evaluate2(x::acb_poly, y::Union{Int,Float64,fmpz,fmpq,arb})
    return evaluate2(x, base_ring(parent(x))(y))
end

###############################################################################
#
#   Composition
#
###############################################################################

function compose(x::acb_poly, y::acb_poly)
   z = parent(x)()
   ccall((:acb_poly_compose, :libarb), Void, 
                (Ptr{acb_poly}, Ptr{acb_poly}, Ptr{acb_poly}, Int),
                &z, &x, &y, prec(parent(x)))
   return z
end

###############################################################################
#
#   Derivative and integral
#
###############################################################################

function derivative(x::acb_poly)
   z = parent(x)()
   ccall((:acb_poly_derivative, :libarb), Void, 
                (Ptr{acb_poly}, Ptr{acb_poly}, Int), &z, &x, prec(parent(x)))
   return z
end

function integral(x::acb_poly)
   z = parent(x)()
   ccall((:acb_poly_integral, :libarb), Void, 
                (Ptr{acb_poly}, Ptr{acb_poly}, Int), &z, &x, prec(parent(x)))
   return z
end

###############################################################################
#
#   Multipoint evaluation and interpolation
#
###############################################################################

function acb_vec(n::Int)
   return ccall((:_acb_vec_init, :libarb), Ptr{acb_struct}, (Int,), n)
end

function acb_vec(b::Array{acb, 1})
   v = ccall((:_acb_vec_init, :libarb), Ptr{acb_struct}, (Int,), length(b))
   for i=1:length(b)
       ccall((:acb_set, :libarb), Void, (Ptr{acb_struct}, Ptr{acb}),
           v + (i-1)*sizeof(acb_struct), &b[i])
   end
   return v
end

function array(R::AcbField, v::Ptr{acb_struct}, n::Int)
   r = Array(acb, n)
   for i=1:n
       r[i] = R()
       ccall((:acb_set, :libarb), Void, (Ptr{acb}, Ptr{acb_struct}),
           &r[i], v + (i-1)*sizeof(acb_struct))
   end
   return r
end

function acb_vec_clear(v::Ptr{acb_struct}, n::Int)
   ccall((:_acb_vec_clear, :libarb), Void, (Ptr{acb_struct}, Int), v, n)
end

function from_roots(R::AcbPolyRing, b::Array{acb, 1})
   z = R()
   tmp = acb_vec(b)
   ccall((:acb_poly_product_roots, :libarb), Void, 
                (Ptr{acb_poly}, Ptr{acb_struct}, Int, Int), &z, tmp, length(b), prec(R))
   acb_vec_clear(tmp, length(b))
   return z
end

function evaluate_iter(x::acb_poly, b::Array{acb, 1})
   return acb[evaluate(x, b[i]) for i=1:length(b)]
end

function evaluate_fast(x::acb_poly, b::Array{acb, 1})
   tmp = acb_vec(b)
   ccall((:acb_poly_evaluate_vec_fast, :libarb), Void, 
                (Ptr{acb_struct}, Ptr{acb_poly}, Ptr{acb_struct}, Int, Int),
            tmp, &x, tmp, length(b), prec(parent(x)))
   res = array(base_ring(parent(x)), tmp, length(b))
   acb_vec_clear(tmp, length(b))
   return res
end

function interpolate_newton(R::AcbPolyRing, xs::Array{acb, 1}, ys::Array{acb, 1})
   length(xs) != length(ys) && error()
   z = R()
   xsv = acb_vec(xs)
   ysv = acb_vec(ys)
   ccall((:acb_poly_interpolate_newton, :libarb), Void, 
                (Ptr{acb_poly}, Ptr{acb_struct}, Ptr{acb_struct}, Int, Int),
            &z, xsv, ysv, length(xs), prec(R))
   acb_vec_clear(xsv, length(xs))
   acb_vec_clear(ysv, length(ys))
   return z
end

function interpolate_barycentric(R::AcbPolyRing, xs::Array{acb, 1}, ys::Array{acb, 1})
   length(xs) != length(ys) && error()
   z = R()
   xsv = acb_vec(xs)
   ysv = acb_vec(ys)
   ccall((:acb_poly_interpolate_barycentric, :libarb), Void, 
                (Ptr{acb_poly}, Ptr{acb_struct}, Ptr{acb_struct}, Int, Int),
            &z, xsv, ysv, length(xs), prec(R))
   acb_vec_clear(xsv, length(xs))
   acb_vec_clear(ysv, length(ys))
   return z
end

function interpolate_fast(R::AcbPolyRing, xs::Array{acb, 1}, ys::Array{acb, 1})
   length(xs) != length(ys) && error()
   z = R()
   xsv = acb_vec(xs)
   ysv = acb_vec(ys)
   ccall((:acb_poly_interpolate_fast, :libarb), Void, 
                (Ptr{acb_poly}, Ptr{acb_struct}, Ptr{acb_struct}, Int, Int),
            &z, xsv, ysv, length(xs), prec(R))
   acb_vec_clear(xsv, length(xs))
   acb_vec_clear(ysv, length(ys))
   return z
end

# todo: cutoffs for fast algorithm
function interpolate(R::AcbPolyRing, xs::Array{acb, 1}, ys::Array{acb, 1})
   return interpolate_newton(R, xs, ys)
end

# todo: cutoffs for fast algorithm
function evaluate(x::acb_poly, b::Array{acb, 1})
   return evaluate_iter(x, b)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function fit!(z::acb_poly, n::Int)
   ccall((:acb_poly_fit_length, :libarb), Void, 
                    (Ptr{acb_poly}, Int), &z, n)
end

function setcoeff!(z::acb_poly, n::Int, x::fmpz)
   ccall((:acb_poly_set_coeff_fmpz, :libarb), Void, 
                    (Ptr{acb_poly}, Int, Ptr{fmpz}), &z, n, &x)
end

function setcoeff!(z::acb_poly, n::Int, x::acb)
   ccall((:acb_poly_set_coeff_acb, :libarb), Void, 
                    (Ptr{acb_poly}, Int, Ptr{acb}), &z, n, &x)
end

function mul!(z::acb_poly, x::acb_poly, y::acb_poly)
   ccall((:acb_poly_mul, :libarb), Void, 
                (Ptr{acb_poly}, Ptr{acb_poly}, Ptr{acb_poly}, Int),
                    &z, &x, &y, prec(parent(z)))
end

function addeq!(z::acb_poly, x::acb_poly)
   ccall((:acb_poly_add, :libarb), Void, 
                (Ptr{acb_poly}, Ptr{acb_poly}, Ptr{acb_poly}, Int),
                    &z, &z, &x, prec(parent(z)))
end

###############################################################################
#
#   Promotions
#
###############################################################################

Base.promote_rule(::Type{acb_poly}, ::Type{Float64}) = acb_poly

Base.promote_rule(::Type{acb_poly}, ::Type{Complex{Float64}}) = acb_poly

Base.promote_rule(::Type{acb_poly}, ::Type{Int}) = acb_poly

Base.promote_rule(::Type{acb_poly}, ::Type{Complex{Int}}) = acb_poly

Base.promote_rule(::Type{acb_poly}, ::Type{fmpz}) = acb_poly

Base.promote_rule(::Type{acb_poly}, ::Type{fmpq}) = acb_poly

Base.promote_rule(::Type{acb_poly}, ::Type{arb}) = acb_poly

Base.promote_rule(::Type{acb_poly}, ::Type{acb}) = acb_poly

Base.promote_rule(::Type{acb_poly}, ::Type{fmpz_poly}) = acb_poly

Base.promote_rule(::Type{acb_poly}, ::Type{fmpq_poly}) = acb_poly

Base.promote_rule(::Type{acb_poly}, ::Type{arb_poly}) = acb_poly

################################################################################
#
#  Parent object call overloads
#
################################################################################

function Base.call(a::AcbPolyRing)
   z = acb_poly()
   z.parent = a
   return z
end

function Base.call(a::AcbPolyRing, b::Union{Int,fmpz,fmpq,Float64,Complex{Float64},Complex{Int},arb,acb})
   z = acb_poly(base_ring(a)(b), a.base_ring.prec)
   z.parent = a
   return z
end

function Base.call(a::AcbPolyRing, b::Array{acb, 1})
   z = acb_poly(b, a.base_ring.prec)
   z.parent = a
   return z
end

function Base.call(a::AcbPolyRing, b::fmpz_poly)
   z = acb_poly(b, a.base_ring.prec)
   z.parent = a
   return z
end

function Base.call(a::AcbPolyRing, b::fmpq_poly)
   z = acb_poly(b, a.base_ring.prec)
   z.parent = a
   return z
end

function Base.call(a::AcbPolyRing, b::arb_poly)
   z = acb_poly(b, a.base_ring.prec)
   z.parent = a
   return z
end

function Base.call(a::AcbPolyRing, b::acb_poly)
   z = acb_poly(b, a.base_ring.prec)
   z.parent = a
   return z
end

################################################################################
#
#  PolynomialRing constructor
#
################################################################################

function PolynomialRing(R::AcbField, s::String)
  S = symbol(s)
  parent_obj = AcbPolyRing(R, S)
  return parent_obj, parent_obj(fmpz_poly([fmpz(0), fmpz(1)]))
end

