module IJuliaExt

using AbstractAlgebra
using IJulia

function Base.showable(mi::MIME"text/latex",
                       x::Union{AbstractAlgebra.Set, AbstractAlgebra.SetElem, AbstractAlgebra.Fac})
  return !IJulia.inited
end

function Base.showable(mi::MIME"text/html",
                       x::Union{AbstractAlgebra.Set, AbstractAlgebra.SetElem, AbstractAlgebra.Fac})
  return !IJulia.inited || AbstractAlgebra.get_html_as_latex()
end

end # module
