using Documenter, AbstractAlgebra

makedocs(
         format   = Documenter.HTML(),
         sitename = "AbstractAlgebra.jl",
         modules = [AbstractAlgebra],
         clean = true,
         doctest = false,
         pages    = [
             "index.md",
             "constructors.md",
             "Rings" => [ "rings.md",
                          "ncrings.md",
                          "euclidean.md",
                          "integer.md",
                          "polynomial_rings.md",
                          "polynomial.md",
			  "ncpolynomial.md",
                          "mpolynomial_rings.md",
                          "mpolynomial.md",
                          "series_rings.md",
                          "series.md",
                          "puiseux.md",
                          "residue_rings.md",
                          "residue.md"],
             "Fields" => [ "fields.md",
                           "fraction_fields.md",
                           "fraction.md",
                           "rational.md",
                           "finfield.md",
                           "real.md",
                           "numberfield.md"],
             "Groups" => [ "perm.md",
                           "ytabs.md"],
             "Modules" => [ "module.md",
                            "free_module.md",
                            "submodule.md",
                            "quotient_module.md",
                            "module_homomorphism.md"],
             "Matrices" => [ "matrix_spaces.md",
                             "matrix.md",
                             "matrix_algebras.md"],
             "Maps" => [ "map.md",
                         "functional_map.md",
                         "map_cache.md",
                         "map_with_inverse.md"],
             "types.md" # Appendix A
         ]
)

deploydocs(
   repo   = "github.com/Nemocas/AbstractAlgebra.jl.git",
   target = "build",
)

