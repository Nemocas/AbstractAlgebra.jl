###########################################################################################
#
#   nf.jl : Antic number fields
#
###########################################################################################

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

NfNumberFieldID = Dict{(fmpz_poly, Symbol), Field}()

type NfNumberField{S} <: Field
   pol::fmpz_poly

   function NfNumberField(pol::fmpz_poly)
      try
         return NfNumberFieldID[pol, S]
      catch
         nf = NfNumberFieldID[pol, S] = new(pol)
         return nf
      end
   end
end

type nf{S} <: NumberFieldElem
   parent::NfNumberField{S}
end

parent(a::nf) = nf.parent

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, a::NfNumberField)
   print(io, "Number field over Rational Field with defining polynomial ", a.pol)
end

###########################################################################################
#
#   NumberField constructor
#
###########################################################################################

function NumberField(pol::fmpz_poly, s::String)
   S = symbol(s)
   parent_obj = NfNumberField{S}(pol)

   return parent_obj, S # gen(parent_obj) 
end


