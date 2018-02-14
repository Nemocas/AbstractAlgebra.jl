using Documenter, AbstractAlgebra

makedocs(
         format   = :html,
         sitename = "AbstractAlgebra.jl",
         pages    = [
             "index.md",
             "constructors.md",
             "rings.md",
             "euclidean.md",
             "polynomial_rings.md",
             "polynomial.md",
             "mpolynomial_rings.md",
             "mpolynomial.md",
             "series_rings.md",
             "series.md",
             "residue_rings.md",
             "residue.md",
             "matrix_spaces.md",
             "matrix.md",
             "types.md" # Appendix A
         ]
)

deploydocs(
   repo   = "github.com/Nemocas/AbstractAlgebra.jl.git",
   target = "build",
   deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
   make   = nothing
)

