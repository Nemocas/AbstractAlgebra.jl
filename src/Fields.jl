###########################################################################################
#
#   Fields.jl : generic fields
#
###########################################################################################

include("Fraction.jl")

include("fmpq.jl")

include("FiniteFields.jl")

//{T <: FieldElem}(a::T, b::T) = divexact(a, b)



