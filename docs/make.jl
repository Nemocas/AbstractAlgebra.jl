using Documenter, AbstractAlgebra

makedocs(
         format   = :html,
         sitename = "AbstractAlgebra.jl",
         pages    = [
             "index.md",
             "constructors.md",
             "Rings" => [ "rings.md",
                          "euclidean.md",
                          "integer.md",
                          "polynomial_rings.md",
                          "polynomial.md",
                          "mpolynomial_rings.md",
                          "mpolynomial.md",
                          "series_rings.md",
                          "series.md",
                          "residue_rings.md",
                          "residue.md"],
             "Fields" => [ "fields.md",
                           "fraction_fields.md",
                           "fraction.md",
                           "rational.md",
                           "finfield.md"],
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

