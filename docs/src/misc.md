```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Miscellaneous

## Printing options

AbstractAlgebra supports printing to LaTeX using the MIME type "text/latex". To
enable LaTeX rendering in Jupyter notebooks and query for the current state,
use the following functions:

```@docs
set_html_as_latex
get_html_as_latex
```
