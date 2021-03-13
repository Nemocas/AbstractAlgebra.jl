# Elements and parents

Any implementation with elements and parents should implement the following
interface:

```@docs
parent
elem_type
parent_type
```

Further, if one has a base ring, like polynomials over the integers
$\mathbb{Z}[x]$, then one should implement

```@docs
base_ring
```
