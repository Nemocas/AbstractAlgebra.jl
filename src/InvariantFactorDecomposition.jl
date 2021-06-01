###############################################################################
#
#   SNFModule.jl : Invariant factor decomposition of modules
#
###############################################################################

export SNFModule, SNFModuleElem, invariant_factors


###############################################################################
#
#   SNFModule constructor
#
###############################################################################

@doc Markdown.doc"""
    snf(m::FPModule{T}) where T <: RingElement

Return a pair `M, f` consisting of the invariant factor decomposition $M$ of
the module `m` and a module homomorphism (isomorphisms) $f : M \to m$. The
module `M` is itself a module which can be manipulated as any other module
in the system.
"""
snf(m::FPModule{T}) where T <: RingElement

@doc Markdown.doc"""
    invariant_factors(m::FPModule{T}) where T <: RingElement

Return a vector of the invariant factors of the module $M$.
"""
function invariant_factors(m::FPModule{T}) where T <: RingElement
   return Generic.invariant_factors(m)
end

