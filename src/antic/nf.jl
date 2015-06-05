###############################################################################
#
#   nf.jl : Antic number fields
#
###############################################################################

export NumberField, norm, trace, CyclotomicField, MaximalRealSubfield, add!,
       sub!, mul!

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent(a::nf_elem) = a.parent

elem_type(::NfNumberField) = nf_elem

base_ring(a::NfNumberField) = None

var(a::NfNumberField) = a.S

function check_parent(a::nf_elem, b::nf_elem) 
   a.parent != b.parent && error("Incompatible number field elements")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function hash(a::nf_elem)
   h = 0xc2a44fbe466a1827
   for i in 1:degree(parent(a)) + 1
         h $= hash(coeff(a, i))
         h = (h << 1) | (h >> (sizeof(Int)*8 - 1))
   end
   return h
end

function coeff(x::nf_elem, n::Int)
   n < 0 && throw(DomainError())
   z = fmpq()
   ccall((:nf_elem_get_coeff_fmpq, :libflint), Void, 
     (Ptr{fmpq}, Ptr{nf_elem}, Int, Ptr{NfNumberField}), &z, &x, n, &parent(x))
   return z
end

function gen(a::NfNumberField)
   r = nf_elem(a)
   ccall((:nf_elem_gen, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &a)
   return r
end

function one(a::NfNumberField)
   r = nf_elem(a)
   ccall((:nf_elem_one, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &a)
   return r
end

function zero(a::NfNumberField)
   r = nf_elem(a)
   ccall((:nf_elem_zero, :libflint), Void,
         (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &a)
   return r
end

function isgen(a::nf_elem)
   return ccall((:nf_elem_is_gen, :libflint), Bool, 
                (Ptr{nf_elem}, Ptr{NfNumberField}), &a, &a.parent)
end

function isone(a::nf_elem)
   return ccall((:nf_elem_is_one, :libflint), Bool, 
                (Ptr{nf_elem}, Ptr{NfNumberField}), &a, &a.parent)
end

function iszero(a::nf_elem)
   return ccall((:nf_elem_is_zero, :libflint), Bool, 
                (Ptr{nf_elem}, Ptr{NfNumberField}), &a, &a.parent)
end

degree(a::NfNumberField) = degree(a.pol)

function deepcopy(d::nf_elem)
   z = nf_elem(parent(d), d)
   return z
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, a::NfNumberField)
   print(io, "Number field over Rational Field")
   print(io, " with defining polynomial ", a.pol)
end

function show(io::IO, x::nf_elem)
   cstr = ccall((:nf_elem_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{nf_elem}, Ptr{Uint8}, Ptr{NfNumberField}), 
                 &x, bytestring(string(var(parent(x)))), &parent(x))

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
end

needs_parentheses(::Nemo.nf_elem) = true

is_negative(::nf_elem) = false

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::nf_elem)
   r = a.parent()
   ccall((:nf_elem_neg, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &a.parent)
   return r
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::nf_elem, b::nf_elem)
   check_parent(a, b)
   r = a.parent()
   ccall((:nf_elem_add, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function -(a::nf_elem, b::nf_elem)
   check_parent(a, b)
   r = a.parent()
   ccall((:nf_elem_sub, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function *(a::nf_elem, b::nf_elem)
   check_parent(a, b)
   r = a.parent()
   ccall((:nf_elem_mul, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function +(a::nf_elem, b::Int)
   r = a.parent()
   ccall((:nf_elem_add_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &r, &a, b, &a.parent)
   return r
end

function +(a::nf_elem, b::fmpz)
   r = a.parent()
   ccall((:nf_elem_add_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function +(a::nf_elem, b::fmpq)
   r = a.parent()
   ccall((:nf_elem_add_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function -(a::nf_elem, b::Int)
   r = a.parent()
   ccall((:nf_elem_sub_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &r, &a, b, &a.parent)
   return r
end

function -(a::nf_elem, b::fmpz)
   r = a.parent()
   ccall((:nf_elem_sub_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function -(a::nf_elem, b::fmpq)
   r = a.parent()
   ccall((:nf_elem_sub_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function -(a::Int, b::nf_elem)
   r = a.parent()
   ccall((:nf_elem_si_sub, :libflint), Void,
         (Ptr{nf_elem}, Int, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, a, &b, &a.parent)
   return r
end

function -(a::fmpz, b::nf_elem)
   r = a.parent()
   ccall((:nf_elem_fmpz_sub, :libflint), Void,
         (Ptr{nf_elem}, Ptr{fmpz}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function -(a::fmpq, b::nf_elem)
   r = a.parent()
   ccall((:nf_elem_fmpq_sub, :libflint), Void,
         (Ptr{nf_elem}, Ptr{fmpq}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

+(a::nf_elem, b::Integer) = a + fmpz(b)

-(a::nf_elem, b::Integer) = a - fmpz(b)

-(a::Integer, b::nf_elem) = fmpz(a) - b

+(a::Integer, b::nf_elem) = b + a

+(a::fmpq, b::nf_elem) = b + a

function *(a::nf_elem, b::Int)
   r = a.parent()
   ccall((:nf_elem_scalar_mul_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &r, &a, b, &a.parent)
   return r
end

function *(a::nf_elem, b::fmpz)
   r = a.parent()
   ccall((:nf_elem_scalar_mul_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function *(a::nf_elem, b::fmpq)
   r = a.parent()
   ccall((:nf_elem_scalar_mul_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

*(a::Integer, b::nf_elem) = b * a

*(a::fmpq, b::nf_elem) = b * a

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::nf_elem, n::Int)
   n < 0 && throw(DomainError)
   r = a.parent()
   ccall((:nf_elem_pow, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &r, &a, n, &a.parent)
   return r
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::nf_elem, b::nf_elem)
   check_parent(a, b)
   return ccall((:nf_elem_equal, :libflint), Bool, 
           (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}), &a, &b, &a.parent)
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::nf_elem)
   r = a.parent()
   ccall((:nf_elem_inv, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &a.parent)
   return r
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::nf_elem, b::nf_elem)
   check_parent(a, b)
   r = a.parent()
   ccall((:nf_elem_div, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(a::nf_elem, b::Int)
   r = a.parent()
   ccall((:nf_elem_scalar_div_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &r, &a, b, &a.parent)
   return r
end

function divexact(a::nf_elem, b::fmpz)
   r = a.parent()
   ccall((:nf_elem_scalar_div_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function divexact(a::nf_elem, b::fmpq)
   r = a.parent()
   ccall((:nf_elem_scalar_div_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

###############################################################################
#
#   Norm and trace
#
###############################################################################

function norm(a::nf_elem)
   z = fmpq()
   ccall((:nf_elem_norm, :libflint), Void,
         (Ptr{fmpq}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &z, &a, &a.parent)
   return z
end

function trace(a::nf_elem)
   z = fmpq()
   ccall((:nf_elem_trace, :libflint), Void,
         (Ptr{fmpq}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &z, &a, &a.parent)
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function mul!(z::nf_elem, x::nf_elem, y::nf_elem)
   ccall((:nf_elem_mul, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}), 
                                                  &z, &x, &y, &parent(x))
end

function addeq!(z::nf_elem, x::nf_elem)
   ccall((:nf_elem_add, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}), 
                                                  &z, &z, &x, &parent(x))
end

function add!(a::nf_elem, b::nf_elem, c::nf_elem)
   ccall((:nf_elem_add, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

function mul!(a::nf_elem, b::nf_elem, c::nf_elem)
   ccall((:nf_elem_mul, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

###############################################################################
#
#   Ad hoc unsafe functions
#
###############################################################################

function add!(c::nf_elem, a::nf_elem, b::fmpq)
   ccall((:nf_elem_add_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

function add!(c::nf_elem, a::nf_elem, b::fmpz)
   ccall((:nf_elem_add_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

function add!(c::nf_elem, a::nf_elem, b::Int)
   ccall((:nf_elem_add_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

add!(c::nf_elem, a::nf_elem, b::Integer) = add!(c, a, fmpz(b))

function sub!(c::nf_elem, a::nf_elem, b::fmpq)
   ccall((:nf_elem_sub_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

function sub!(c::nf_elem, a::nf_elem, b::fmpz)
   ccall((:nf_elem_sub_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

function sub!(c::nf_elem, a::nf_elem, b::Int)
   ccall((:nf_elem_sub_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

sub!(c::nf_elem, a::nf_elem, b::Integer) = sub!(c, a, fmpz(b))

function sub!(c::nf_elem, a::fmpq, b::nf_elem)
   ccall((:nf_elem_fmpq_sub, :libflint), Void,
         (Ptr{nf_elem}, Ptr{fmpq}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

function sub!(c::nf_elem, a::fmpz, b::nf_elem)
   ccall((:nf_elem_fmpz_sub, :libflint), Void,
         (Ptr{nf_elem}, Ptr{fmpz}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

function sub!(c::nf_elem, a::Int, b::nf_elem)
   ccall((:nf_elem_si_sub, :libflint), Void,
         (Ptr{nf_elem}, Int, Ptr{nf_elem}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

sub!(c::nf_elem, a::Integer, b::nf_elem) = sub!(c, fmpz(a), b)

function mul!(c::nf_elem, a::nf_elem, b::fmpq)
   ccall((:nf_elem_scalar_mul_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

function mul!(c::nf_elem, a::nf_elem, b::fmpz)
   ccall((:nf_elem_scalar_mul_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz}, Ptr{NfNumberField}),
         &c, &a, &b, &a.parent)
end

function mul!(c::nf_elem, a::nf_elem, b::Int) 
   ccall((:nf_elem_scalar_mul_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &c, &a, b, &a.parent)
end

mul!(c::nf_elem, a::nf_elem, b::Integer) = mul!(c, a, fmpz(b))

###############################################################################
#
#   Promotions
#
###############################################################################

Base.promote_rule{T <: Integer}(::Type{nf_elem}, ::Type{T}) = nf_elem

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function Base.call(a::NfNumberField)
   z = nf_elem(a)
   return z
end

function Base.call(a::NfNumberField, c::Int)
   z = nf_elem(a)
   ccall((:nf_elem_set_si, :libflint), Void, 
         (Ptr{nf_elem}, Int, Ptr{NfNumberField}), &z, c, &a)
   return z
end

function Base.call(a::NfNumberField, c::fmpz)
   z = nf_elem(a)
   ccall((:nf_elem_set_fmpz, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{fmpz}, Ptr{NfNumberField}), &z, &c, &a)
   return z
end

function Base.call(a::NfNumberField, c::fmpq)
   z = nf_elem(a)
   ccall((:nf_elem_set_fmpq, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{fmpq}, Ptr{NfNumberField}), &z, &c, &a)
   return z
end

function Base.call(a::NfNumberField, b::nf_elem)
   parent(b) != a && error("Cannot coerce number field element")
   return b
end

function Base.call(a::NfNumberField, pol::fmpq_poly)
   pol = parent(a.pol)(pol) # check pol has correct parent
   z = nf_elem(a)
   if length(pol) >= length(a.pol)
      pol = mod(pol, a.pol)
   end
   ccall((:nf_elem_set_fmpq_poly, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{fmpq_poly}, Ptr{NfNumberField}), &z, &pol, &a)
   return z
end

function Base.call(a::FmpqPolyRing, b::nf_elem)
   parent(parent(b).pol) != a && error("Cannot coerce from number field to polynomial ring")
   r = a()
   ccall((:fmpq_poly_set, :libflint), Void, 
         (Ptr{fmpq_poly}, Ptr{nf_elem}), &r, &b)
   return r
end

###############################################################################
#
#   NumberField constructor
#
###############################################################################

function NumberField(pol::fmpq_poly, s::String)
   S = symbol(s)
   parent_obj = NfNumberField(pol, S)

   return parent_obj, gen(parent_obj) 
end

function CyclotomicField(n::Int, s::String, t = "\$")
   Zx, x = PolynomialRing(flintZZ, string(gensym()))
   Qx, = PolynomialRing(flintQQ, t)
   f = cyclotomic(n, x)
   return NumberField(Qx(f), s)
end

function MaximalRealSubfield(n::Int, s::String, t = "\$")
   Zx, x = PolynomialRing(flintZZ, string(gensym()))
   Qx, = PolynomialRing(flintQQ, t)
   f = cos_minpoly(n, x)
   return NumberField(Qx(f), s)
end
