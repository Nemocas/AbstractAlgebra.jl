###############################################################################
#
#   AnticTypes.jl : Antic types
#
###############################################################################

###############################################################################
#
#   NfNumberField / nf_elem
#
###############################################################################

NfNumberFieldID = Dict{Tuple{fmpq_poly, Symbol}, Field}()

type NfNumberField <: Field{Antic}
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
   pol::fmpq_poly
   S::Symbol

   function NfNumberField(pol::fmpq_poly, s::Symbol)
      try
         return NfNumberFieldID[pol, s]
      catch
         nf = NfNumberFieldID[pol, s] = new()
         nf.pol = pol
         ccall((:nf_init, :libflint), Void, 
            (Ptr{NfNumberField}, Ptr{fmpq_poly}), &nf, &pol)
         finalizer(nf, _NfNumberField_clear_fn)
         nf.S = s
         return nf
      end
   end
end

function _NfNumberField_clear_fn(a::NfNumberField)
   ccall((:nf_clear, :libflint), Void, (Ptr{NfNumberField},), &a)
end

type nf_elem <: NumberFieldElem
   elem_coeffs::Ptr{Void}
   elem_den::Int
   elem_alloc::Int
   elem_length::Int
   parent::NfNumberField

   function nf_elem(p::NfNumberField)
      r = new()
      ccall((:nf_elem_init, :libflint), Void, 
            (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &p)
      r.parent = p
      finalizer(r, _nf_elem_clear_fn)
      return r
   end

   function nf_elem(p::NfNumberField, a::nf_elem)
      r = new()
      ccall((:nf_elem_init, :libflint), Void, 
            (Ptr{nf_elem}, Ptr{NfNumberField}), &r, &p)
      ccall((:nf_elem_set, :libflint), Void,
            (Ptr{nf_elem}, Ptr{nf_elem}, Ptr{NfNumberField}), &r, &a, &p)
      r.parent = p
      finalizer(r, _nf_elem_clear_fn)
      return r
   end
end

function _nf_elem_clear_fn(a::nf_elem)
   ccall((:nf_elem_clear, :libflint), Void, 
         (Ptr{nf_elem}, Ptr{NfNumberField}), &a, &a.parent)
end
