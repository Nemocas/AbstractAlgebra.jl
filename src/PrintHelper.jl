"""
    pluralize(noun::String)

Return the plural form of a given english noun in singular.

This function employs certain heuristics and thus may not always
obtain the correct pluralization. But it works well enough in the
vast majority of cases.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> pluralize("generator")
"generators"

julia> pluralize("variety")
"varieties"
```
"""
function pluralize(noun::String)
  noun == "child" && return "children"

  # matrix -> matrices
  endswith(noun, "ix") && return noun[1:end-2] * "ices"

  # simplex -> simplices
  # vertex -> vertices
  endswith(noun, "ex") && return noun[1:end-2] * "ices"

  # polyhedron -> polyhedra
  endswith(noun, "ron") && return noun[1:end-2] * "a"

  # basis -> bases
  endswith(noun, "sis") && return noun[1:end-3] * "ses"

  # maximum -> maxima
  endswith(noun, "um") && return noun[1:end-2] * "a"

  # patch -> patches
  endswith(noun, "tch") && return noun * "es"

  # family -> families
  # variety -> varieties
  endswith(noun, "y") && !(noun[end-1] in "aeiouy") && return noun[1:end-1] * "ies"

  # otherwise fall back to the default and add 's'
  return noun * "s"
end


# The following code is inspired by https://github.com/TotalVerb/EnglishText.jl/
#
# We want to output a quantity followed by a noun, putting the noun into
# plural form if necessary. We use a dedicated type for this with a `show`
# method. This way we avoid creating a temporary string when this code is used
# for printing to an IO stream.

"""
    ItemQuantity(count::Int, noun::String)
    ItemQuantity(count::Int, noun::String, noun_plural::String)

A helper object which has the sole purpose of being printed. If `count` is `1`
then it is printed as `1 noun`. For any other `count`, the output will be
`count noun_plural`. If `noun_plural` is not given, this defaults to
`plural(noun)`.

The reason we allow specifying an explicit plural form is that
[`plural`](@ref) is not perfect and there may be situations where it does not
produce the correct plural form, thus a fallback alternative is prudent to
have. Ideally, though, please instead improve `pluralize` to handle your needs
correctly.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> ItemQuantity(0, "generator")
0 generators

julia> ItemQuantity(1, "generator")
1 generator

julia> ItemQuantity(2, "generator")
2 generators
```

Here is an example with a custom plural form.
```jldoctest; setup = :(using AbstractAlgebra)
julia> ItemQuantity(0, "ox", "oxen")
0 oxen

julia> ItemQuantity(1, "ox", "oxen")
1 ox

julia> ItemQuantity(2, "ox", "oxen")
2 oxen
```
"""
struct ItemQuantity
    count::Int
    noun::String
    noun_plural::String
    ItemQuantity(count::Int, noun::String) = new(count, noun)
    ItemQuantity(count::Int, noun::String, noun_plural::String) = new(count, noun, noun_plural)
end

function Base.show(io::IO, quantity::ItemQuantity)
    print(io, quantity.count)
    print(io, " ")
    if quantity.count == 1
        print(io, quantity.noun)
    elseif isdefined(quantity, :noun_plural)
        print(io, quantity.noun_plural)
    else
        print(io, pluralize(quantity.noun))
    end
end

"""
    ordinal_number_string(number::Int)

A helper function which returns `number` as a string in ordinal form.

# Examples
```jldoctest; setup = :(using AbstractAlgebra; ordinal_number_string=AbstractAlgebra.ordinal_number_string)
julia> ordinal_number_string(1)
"1st"

julia> ordinal_number_string(2)
"2nd"

julia> ordinal_number_string(3)
"3rd"

julia> ordinal_number_string(4)
"4th"
```
"""
function ordinal_number_string(number::Int)
    number >= 0 || error("number must be non-negative")
    mod100 = mod(number, 100)
    mod10 = mod(number, 10)
    if mod100 in 11:13
        return "$(number)th"
    elseif mod10 == 1
        return "$(number)st"
    elseif mod10 == 2
        return "$(number)nd"
    elseif mod10 == 3
        return "$(number)rd"
    else
        return "$(number)th"
    end
end
