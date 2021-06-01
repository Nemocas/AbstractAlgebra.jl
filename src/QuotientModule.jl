###############################################################################
#
#   QuotientModule.jl : Quotients of modules by submodules
#
###############################################################################

export QuotientModule, QuotientModuleElem, quo

###############################################################################
#
#   QuotientModule constructor
#
###############################################################################

@doc Markdown.doc"""
    quo(m::FPModule{T}, subm::FPModule{T}) where T <: RingElement

Return the quotient `M` of the module `m` by the module `subm` (which must
have been (transitively) constructed as a submodule of `m` or be `m` itself)
along with the canonical quotient map from `m` to `M`.
"""
function quo(m::FPModule{T}, subm::FPModule{T}) where T <: RingElement
   return Generic.quo(m, subm)  
end
