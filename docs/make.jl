using Documenter, AbstractAlgebra

makedocs(
         format   = :html,
         sitename = "AbstractAlgebra.jl",
         pages    = [
             "index.md",
             "about.md",
             "types.md",
             "constructors.md"
         ]
)

deploydocs(
   repo   = "github.com/Nemocas/AbstractAlgebra.jl.git",
   target = "build",
   deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
   make   = nothing
)

