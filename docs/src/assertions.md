```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Assertion and Verbosity Macros

We describe here various macros provided by AbstractAlgebra.

## Verbosity macros
There is a (global) list of symbols called *verbosity scopes* which represent keywords used to
trigger some verbosity macros within the code. Each of these verbosity scopes has its own
integer *verbosity level*, which is set to $0$ by default. A verbosity macro call
must specify the verbosity scope `S` and optionally the trigger level `k` (defaulting to $1$) such that,
if the current verbosity level `l` of `S` is bigger than or equal to `k`, then the
macro triggers a given action.  Inside a module, the function `add_verbosity_scope` must be
called in the `__init__` function of that module.

```@docs
add_verbosity_scope(s::Symbol)
set_verbosity_level(s::Symbol, l::Int)
get_verbosity_level(s::Symbol)
```

### Printings

```@docs
@vprintln
@vprint
```

### Actions

```@docs
@v_do
```

## Assertion macros
There is a list of symbols called *assertion scopes* which represent keywords used to
trigger some particular macros within the codes. Each of these assertion scopes is
associated with an *assertion level*, being set to $0$ by default. An assertion macro
is joined to an assertion scope `S` and a value `k` (set to $1$ by default) such that,
if the current assertion level `l` of `S` is bigger than or equal to `k`, then the
macro triggers an action on the given assertion

```@docs
add_assertion_scope(s::Symbol)
set_assertion_level(s::Symbol, l::Int)
get_assertion_level(s::Symbol)
```

### Check

```@docs
@hassert
```

## Miscellaneous

```@docs
@req
```
