using Documenter, AbstractAlgebra

DocMeta.setdocmeta!(AbstractAlgebra, :DocTestSetup, AbstractAlgebra.doctestsetup(); recursive = true)

let dir = joinpath(@__DIR__, "src", "assets")
    for n in ("parents_diagram", "elements_diagram")
        svg = joinpath(dir, n * ".svg"); pdf = joinpath(dir, n * ".pdf")
        isfile(svg) || (@warn "SVG missing, skipping" svg; continue)
        (!isfile(pdf) || stat(svg).mtime > stat(pdf).mtime) &&
            (run(`rsvg-convert -f pdf -o $pdf $svg`))
    end
end

makedocs(
         format = [
            Documenter.HTML(;
            size_threshold_warn = 204800,
            size_threshold = 409600,
         ),
         Documenter.LaTeX(platform="tectonic"),
         ],
         sitename = "AbstractAlgebra.jl",
         modules = [AbstractAlgebra],
         clean = true,
         doctest = true,
         checkdocs = :none,
         pages    = [
             "index.md",
             "Fundamental interface of AbstractAlgebra.jl" => [
                 "types.md",
                 "visualizing_types.md",
                 "extending_abstractalgebra.md",
             ],
             "constructors.md",
             "Rings" => [
                 "ring_introduction.md",
                 "ring.md",
                 "integer.md",
                 "total_fraction.md",
                 "polynomial.md",
                 "ncpolynomial.md",
                 "mpolynomial.md",
                 "univpolynomial.md",
                 "laurent_polynomial.md",
                 "laurent_mpolynomial.md",
                 "series.md",
                 "puiseux.md",
                 "mseries.md",
                 "residue.md",
                 "free_associative_algebra.md",
             ],
             "Fields" => [
                 "field_introduction.md",
                 "field.md",
                 "fraction.md",
                 "rational.md",
                 "function_field.md",
                 "finfield.md",
                 "real.md",
             ],
             "Groups" => [
                 "perm.md",
                 "ytabs.md",
             ],
             "Modules" => [
                 "module_introduction.md",
                 "module.md",
                 "free_module.md",
                 "submodule.md",
                 "quotient_module.md",
                 "direct_sum.md",
                 "module_homomorphism.md",
             ],
             "Ideals" => [
                 "ideal.md",
             ],
             "Matrices" => [
                 "matrix_introduction.md",
                 "matrix.md",
                 "matrix_spaces.md",
                 "matrix_implementation.md",
                 "matrix_algebras.md",
             ],
             "Maps" => [
                 "map_introduction.md",
                 "functional_map.md",
                 "map_cache.md",
                 "map_with_inverse.md",
             ],
             "Miscellaneous" => [
                 "misc.md",
                 "assertions.md"
             ],
             "Interfaces" => [
                 "interface_introduction.md",
                 "ring_interface.md",
                 "euclidean_interface.md",
                 "poly_interface.md",
                 "mpoly_interface.md",
                 "series_interface.md",
                 "residue_interface.md",
                 "field_interface.md",
                 "fraction_interface.md",
                 "module_interface.md",
                 "ideal_interface.md",
                 "matrix_interface.md",
                 "map_interface.md",
                 "rand.md",
                 "linear_solving.md",
             ],
         ]
)

deploydocs(
   repo   = "github.com/Nemocas/AbstractAlgebra.jl.git",
   target = "build",
)
