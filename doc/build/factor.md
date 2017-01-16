


<a id='Introduction-1'></a>

## Introduction


Nemo provides a unified interface to handle factorizations using the `Fact` objects. These can only be constructed using the factor function for the respective ring elements. This is best illustrated by an example.


```julia
julia> fac = factor(-6000361807272228723606)
-1 * 2 * 229^100 * 43669^100 * 3

julia> unit(fac)
-1

julia> -6000361807272228723606 == unit(fac) * prod([ p^e for (p, e) in fac])
true

julia> for (p, e) in fac; println("$p $e"); end
2 1
229 3
43669 3
3 1

julia> 229 in fac
true

julia> fac[229]
3
```


<a id='Basic-functionality-1'></a>

## Basic functionality


Objects of type `Fac` are iterable, that is, if `a` is an object of type `Fac`, then `for (p, e) in a` will iterate through all pairs `(p, e)`, where `p` is a factor and `e` the corresponding exponent.


```
getindex(collection, key...)
```

Retrieve the value(s) stored at the given key or index within a collection. The syntax `a[i,j,...]` is converted by the compiler to `getindex(a, i, j, ...)`.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/3c9d75391c72d7c32eea75ff187ce77b2d5effc8/base/docs/helpdb/Base.jl#L466-471' class='documenter-source'>source</a><br>


```
getindex(a::Fac{T}, b::T) -> Int
```

> If `b` is a factor of `a`, the corresponding exponent is returned. Otherwise an error is thrown.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/Factor.jl#L44' class='documenter-source'>source</a><br>


```
length(collection) -> Integer
```

For ordered, indexable collections, the maximum index `i` for which `getindex(collection, i)` is valid. For unordered collections, the number of elements.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/3c9d75391c72d7c32eea75ff187ce77b2d5effc8/base/docs/helpdb/Base.jl#L5900-5905' class='documenter-source'>source</a><br>


```
length(a::Fac) -> Int
```

> Returns the number of factors of `a`, not including the unit.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/Factor.jl#L86' class='documenter-source'>source</a><br>


```
unit(a::Fac{T}) -> T
```

> Returns the unit of the factorization.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/Factor.jl#L20' class='documenter-source'>source</a><br>

