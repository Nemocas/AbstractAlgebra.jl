using Documenter, AbstractAlgebra

makedocs()

deploydocs(
   deps   = Deps.pip("mkdocs", "python-markdown-math"),
   repo = "github.com/Nemocas/AbstractAlgebra.jl.git",
   julia = "0.6"
)

