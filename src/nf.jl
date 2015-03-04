###########################################################################################
#
#   nf.jl : Antic number fields
#
###########################################################################################

export norm, trace

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

NfNumberFieldID = Dict{(fmpq_poly, Symbol, Symbol), Field}()

type NfNumberField{S, T} <: Field
   pol_coeffs::Ptr{Void}
   pol_den::Int
   pol_alloc::Int
   pol_length::Int
   pinv_dinv::Ptr{Void}
   pinv_n::Int
   pinv_norm::Int
   powers::Ptr{Void}
   powers_len::Int
   traces_coeffs::Ptr{Void}
   traces_den::Int
   traces_alloc::Int
   traces_length::Int
   flag::Uint
   pol::fmpq_poly{T}

   function NfNumberField(pol::fmpq_poly)
      try
         return NfNumberFieldID[pol, S, T]
      catch
         nf = NfNumberFieldID[pol, S, T] = new()
         nf.pol = pol
         ccall((:nf_init, :libflint), Void, 
            (Ptr{NfNumberField}, Ptr{fmpq_poly}), &nf, &pol)
         finalizer(nf, _NfNumberField_clear_fn)
         return nf
      end
   end
end

_NfNumberField_clear_fn(a::NfNumberField) = ccall((:nf_clear, :libflint), Void, 
                       (Ptr{NfNumberField},), &a)

type nf_elem{S, T} <: NumberFieldElem
   elem_coeffs::Ptr{Void}
   elem_den::Int
   elem_alloc::Int
   elem_length::Int
   parent::NfNumberField{S, T}

   function nf_elem(p::NfNumberField)
      r = new()
      ccall((:nf_elem_init, :libflint), Void, (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &p)
      r.parent = p
      finalizer(r, _nf_elem_clear_fn)
      return r
   end
end

_nf_elem_clear_fn(a::nf_elem) = ccall((:nf_elem_clear, :libflint), Void, 
                                      (Ptr{nf_elem}, Ptr{NfNumberField}), &a, &a.parent)

parent(a::nf_elem) = a.parent

elem_type{S, T}(::NfNumberField{S, T}) = nf_elem{S, T}

base_ring{S, T}(a::NfNumberField{S, T}) = None

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function gen{S, T}(a::NfNumberField{S, T})
   r = nf_elem{S, T}(a)
   ccall((:nf_elem_gen, :libflint), Void, (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &a)
   return r
end

function one{S, T}(a::NfNumberField{S, T})
   r = nf_elem{S, T}(a)
   ccall((:nf_elem_one, :libflint), Void, (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &a)
   return r
end

function zero{S, T}(a::NfNumberField{S, T})
   r = nf_elem{S, T}(a)
   ccall((:nf_elem_zero, :libflint), Void, (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &a)
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

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, a::NfNumberField)
   print(io, "Number field over Rational Field with defining polynomial ", a.pol)
end

function show{S, T}(io::IO, x::nf_elem{S, T})
   cstr = ccall((:nf_elem_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{nf_elem{S, T}}, Ptr{Uint8}, Ptr{NfNumberField{S, T}}), 
                 &x, bytestring(string(S)), &x.parent)

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
end

needs_parentheses(::Nemo.nf_elem) = true

is_negative{S, T}(::nf_elem{S, T}) = false

###########################################################################################
#
#   Unary operators
#
###########################################################################################

function -(a::nf_elem)
   r = a.parent()
   ccall((:nf_elem_neg, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &a.parent)
   return r
end

###########################################################################################
#
#   Binary operators
#
###########################################################################################

function +(a::nf_elem, b::nf_elem)
   r = a.parent()
   ccall((:nf_elem_add, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function -(a::nf_elem, b::nf_elem)
   r = a.parent()
   ccall((:nf_elem_sub, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

function *{S, T}(a::nf_elem{S, T}, b::nf_elem{S, T})
   r = a.parent()
   ccall((:nf_elem_mul, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function +(a::nf_elem, b::Int)
   r = a.parent()
   ccall((:nf_elem_add_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &r, &a, b, &a.parent)
   return r
end

function +(a::nf_elem, b::BigInt)
   r = a.parent()
   temp = fmpz_readonly(b)
   ccall((:nf_elem_add_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz_readonly}, Ptr{NfNumberField}),
         &r, &a, &temp, &a.parent)
   return r
end

function +(a::nf_elem, b::Rational{BigInt})
   r = a.parent()
   temp = fmpq_readonly(b)
   ccall((:nf_elem_add_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq_readonly}, Ptr{NfNumberField}),
         &r, &a, &temp, &a.parent)
   return r
end

function -(a::nf_elem, b::Int)
   r = a.parent()
   ccall((:nf_elem_sub_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &r, &a, b, &a.parent)
   return r
end

function -(a::nf_elem, b::BigInt)
   r = a.parent()
   temp = fmpz_readonly(b)
   ccall((:nf_elem_sub_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz_readonly}, Ptr{NfNumberField}),
         &r, &a, &temp, &a.parent)
   return r
end

function -(a::nf_elem, b::Rational{BigInt})
   r = a.parent()
   temp = fmpq_readonly(b)
   ccall((:nf_elem_sub_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq_readonly}, Ptr{NfNumberField}),
         &r, &a, &temp, &a.parent)
   return r
end

function -(a::Int, b::nf_elem)
   r = a.parent()
   ccall((:nf_elem_si_sub, :libflint), Void,
         (Ptr{nf_elem}, Int, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, a, &b, &a.parent)
   return r
end

function -(a::BigInt, b::nf_elem)
   r = a.parent()
   temp = fmpz_readonly(a)
   ccall((:nf_elem_fmpz_sub, :libflint), Void,
         (Ptr{nf_elem}, Ptr{fmpz_readonly}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &temp, &b, &a.parent)
   return r
end

function -(a::Rational{BigInt}, b::nf_elem)
   r = a.parent()
   temp = fmpq_readonly(a)
   ccall((:nf_elem_fmpq_sub, :libflint), Void,
         (Ptr{nf_elem}, Ptr{fmpq_readonlyl}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &temp, &b, &a.parent)
   return r
end

+(a::nf_elem, b::Integer) = a + BigInt(b)

-(a::nf_elem, b::Integer) = a - BigInt(b)

-(a::Integer, b::nf_elem) = BigInt(a) - b

+(a::Integer, b::nf_elem) = b + a

+(a::Rational{BigInt}, b::nf_elem) = b + a

function *(a::nf_elem, b::Int)
   r = a.parent()
   ccall((:nf_elem_scalar_mul_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &r, &a, b, &a.parent)
   return r
end

function *(a::nf_elem, b::BigInt)
   r = a.parent()
   temp = fmpz_readonly(b)
   ccall((:nf_elem_scalar_mul_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz_readonly}, Ptr{NfNumberField}),
         &r, &a, &temp, &a.parent)
   return r
end

function *(a::nf_elem, b::Rational{BigInt})
   r = a.parent()
   temp = fmpq_readonly(b)
   ccall((:nf_elem_scalar_mul_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq_readonly}, Ptr{NfNumberField}),
         &r, &a, &temp, &a.parent)
   return r
end

*(a::Integer, b::nf_elem) = b * a

*(a::Rational{BigInt}, b::nf_elem) = b * a

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^(a::nf_elem, n::Int)
   n < 0 && throw(DomainError)
   r = a.parent()
   ccall((:nf_elem_pow, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &r, &a, n, &a.parent)
   return r
end

###########################################################################################
#
#   Comparison
#
###########################################################################################

function =={S, T}(a::nf_elem{S, T}, b::nf_elem{S, T})
   return ccall((:nf_elem_equal, :libflint), Bool, 
                (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}), &a, &b, &a.parent)
end

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{S, T}(a::nf_elem{S, T})
   r = a.parent()
   ccall((:nf_elem_inv, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &a.parent)
   return r
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact(a::nf_elem, b::nf_elem)
   r = a.parent()
   ccall((:nf_elem_div, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &r, &a, &b, &a.parent)
   return r
end

###########################################################################################
#
#   Ad hoc exact division
#
###########################################################################################

function divexact(a::nf_elem, b::Int)
   r = a.parent()
   ccall((:nf_elem_scalar_div_si, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Int, Ptr{NfNumberField}),
         &r, &a, b, &a.parent)
   return r
end

function divexact(a::nf_elem, b::BigInt)
   r = a.parent()
   temp = fmpz_readonly(b)
   ccall((:nf_elem_scalar_div_fmpz, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpz_readonly}, Ptr{NfNumberField}),
         &r, &a, &temp, &a.parent)
   return r
end

function divexact(a::nf_elem, b::Rational{BigInt})
   r = a.parent()
   temp = fmpq_readonly(b)
   ccall((:nf_elem_scalar_div_fmpq, :libflint), Void,
         (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{fmpq_readonly}, Ptr{NfNumberField}),
         &r, &a, &temp, &a.parent)
   return r
end

###########################################################################################
#
#   Norm and trace
#
###########################################################################################

function norm{S,T}(a::nf_elem{S, T})
   temp = fmpq()
   ccall((:nf_elem_norm, :libflint), Void,
         (Ptr{fmpq}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &temp, &a, &a.parent)
   return Rational(temp)
end

function trace{S, T}(a::nf_elem{S, T})
   temp = fmpq()
   ccall((:nf_elem_trace, :libflint), Void,
         (Ptr{fmpq}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &temp, &a, &a.parent)
   return Rational(temp)
end

###########################################################################################
#
#   Promotions
#
###########################################################################################

Base.promote_rule{S, T, U <: Integer}(::Type{nf_elem{S, T}}, ::Type{U}) = nf_elem{S, T}

###########################################################################################
#
#   Parent object call overloads
#
###########################################################################################

function Base.call{S, T}(a::NfNumberField{S, T})
   z = nf_elem{S, T}(a)
   return z
end

function Base.call{S, T}(a::NfNumberField{S, T}, c::Int)
   z = nf_elem{S, T}(a)
   ccall((:nf_elem_set_si, :libflint), Void, 
         (Ptr{nf_elem}, Int, Ptr{NfNumberField}), &z, c, &a)
   return z
end

function Base.call{S, T}(a::NfNumberField{S, T}, c::BigInt)
   z = nf_elem{S, T}(a)
   temp = fmpz_readonly(c)
   ccall((:nf_elem_set_fmpz, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{fmpz_readonly}, Ptr{NfNumberField}), &z, &temp, &a)
   return z
end

function Base.call{S, T}(a::NfNumberField{S, T}, c::Rational{BigInt})
   z = nf_elem{S, T}(a)
   temp = fmpq_readonly(c)
   ccall((:nf_elem_set_fmpq, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{fmpq_readonly}, Ptr{NfNumberField}), &z, &temp, &a)
   return z
end

Base.call{S, T}(a::NfNumberField{S, T}, b::nf_elem{S, T}) = b

function Base.call{S, T}(a::NfNumberField{S, T}, pol::fmpq_poly)
   pol = parent(a.pol)(pol) # check pol has correct parent
   z = nf_elem{S, T}(a)
   if length(pol) >= length(a.pol)
      pol = mod(pol, a.pol)
   end
   ccall((:nf_elem_set_fmpq_poly, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{fmpq_poly}, Ptr{NfNumberField}), &z, &pol, &a)
   return z
end

function Base.call(a::FmpqPolyRing, b::nf_elem)
   b.parent.pol.parent != a && error("Cannot coerce from number field to polynomial ring")
   r = a()
   ccall((:fmpq_poly_set, :libflint), Void, (Ptr{fmpq_poly}, Ptr{nf_elem}), &r, &b)
   return r
end

###########################################################################################
#
#   NumberField constructor
#
###########################################################################################

function NumberField{T}(pol::fmpq_poly{T}, s::String)
   S = symbol(s)
   parent_obj = NfNumberField{S, T}(pol)

   return parent_obj, gen(parent_obj) 
end


