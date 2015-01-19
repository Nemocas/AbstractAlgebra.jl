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

NfNumberFieldID = Dict{(fmpz_poly, Symbol), Field}()

type NfNumberField{S} <: Field
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
   pol::fmpz_poly

   function NfNumberField(pol::fmpz_poly)
      polq = fmpq_poly{S}(pol)
      polq.parent = FmpqPolyRing{Rational{BigInt}, S}(QQ)
      try
         return NfNumberFieldID[pol, S]
      catch
         nf = NfNumberFieldID[pol, S] = new()
         nf.pol = pol
         ccall((:nf_init, :libflint), Void, 
            (Ptr{NfNumberField}, Ptr{fmpq_poly}), &nf, &polq)
         finalizer(nf, _NfNumberField_clear_fn)
         return nf
      end
   end
end

_NfNumberField_clear_fn(a::NfNumberField) = ccall((:nf_clear, :libflint), Void, 
                       (Ptr{NfNumberField},), &a)

type nf_elem{S} <: NumberFieldElem
   elem_coeffs::Ptr{Void}
   elem_den::Int
   elem_alloc::Int
   elem_length::Int
   parent::NfNumberField{S}

   function nf_elem(p::NfNumberField)
      r = new()
      ccall((:nf_elem_init, :libflint), Void, (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &p)
      r.parent = p
      finalizer(r, _nf_elem_clear_fn)
      return r
   end
end

_nf_elem_clear_fn(a::nf_elem) = ccall((:nf_elem_clear, :libflint), Void, 
                                      (Ptr{nf_elem},), &a)

parent(a::nf_elem) = a.parent

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function gen{S}(a::NfNumberField{S})
   r = nf_elem{S}(a)
   ccall((:nf_elem_gen, :libflint), Void, (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &a)
   return r
end

function one{S}(a::NfNumberField{S})
   r = nf_elem{S}(a)
   ccall((:nf_elem_one, :libflint), Void, (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &a)
   return r
end

function zero{S}(a::NfNumberField{S})
   r = nf_elem{S}(a)
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

function show{S}(io::IO, x::nf_elem{S})
   cstr = ccall((:nf_elem_get_str_pretty, :libflint), Ptr{Uint8}, 
                (Ptr{nf_elem{S}}, Ptr{Uint8}, Ptr{NfNumberField{S}}), 
                 &x, bytestring(string(S)), &x.parent)

   print(io, bytestring(cstr))

   ccall((:flint_free, :libflint), Void, (Ptr{Uint8},), cstr)
end

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

function *(a::nf_elem, b::nf_elem)
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

function =={S}(a::nf_elem{S}, b::nf_elem{S})
   return ccall((:nf_elem_equal, :libflint), Bool, 
                (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}), &a, &b, &a.parent)
end

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv(a::nf_elem)
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

function norm(a::nf_elem)
   temp = fmpq()
   ccall((:nf_elem_norm, :libflint), Void,
         (Ptr{fmpq}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &temp, &a, &a.parent)
   return Rational(temp)
end

function trace(a::nf_elem)
   temp = fmpq()
   ccall((:nf_elem_trace, :libflint), Void,
         (Ptr{fmpq}, Ptr{nf_elem}, Ptr{NfNumberField}),
         &temp, &a, &a.parent)
   return Rational(temp)
end

###########################################################################################
#
#   Parent object call overloads
#
###########################################################################################

function Base.call{S}(a::NfNumberField{S})
   z = nf_elem{S}(a)
   return z
end

###########################################################################################
#
#   NumberField constructor
#
###########################################################################################

function NumberField(pol::fmpz_poly, s::String)
   S = symbol(s)
   parent_obj = NfNumberField{S}(pol)

   return parent_obj, gen(parent_obj) 
end


