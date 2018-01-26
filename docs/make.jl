using Documenter, AbstractAlgebra

makedocs(
         format   = :html,
         sitename = "AbstractAlgebra.jl",
         pages    = [
             "index.md",
             "constructors.md",
             "rings.md",
             "types.md" # Appendix A
         ]
)

deploydocs(
   repo   = "github.com/Nemocas/AbstractAlgebra.jl.git",
   target = "build",
   deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
   make   = nothing
)

