using AbstractAlgebra
using InteractiveUtils: subtypes
using Kroki

###############################################################################
#
#   The hierarchy, read off the loaded package
#
###############################################################################

# `Module{T<:NCRingElement}` is drawn as `Module{T}`; the bounds would only add
# noise to the picture.
function type_label(T)
  vars = String[]
  while T isa UnionAll
    push!(vars, string(T.var.name))
    T = T.body
  end

  isempty(vars) && return string(nameof(T))
  return string(nameof(T), "{", join(vars, ", "), "}")
end

# `IdealElem` is unused; drop this once #2524 has removed it.
const IGNORED = ["IdealElem{T}"]

# Only AbstractAlgebra's own abstract types.  `Generic`'s refinements and the
# concrete types are covered by the prose in visualizing_types.md instead.
is_shown(T) =
  isabstracttype(T) && parentmodule(T) === AbstractAlgebra && !(type_label(T) in IGNORED)

# Short branches first, so the long chains cascade towards one side instead of
# splitting the picture down the middle.
subtypes_shown(T) = sort!(filter(is_shown, subtypes(T)); by = S -> (leaves(S), type_label(S)))

leaves(T) = (S = subtypes_shown(T); isempty(S) ? 1 : sum(leaves, S))

###############################################################################
#
#   PlantUML
#
###############################################################################

const PREAMBLE = """
@startuml
skinparam monochrome true
skinparam defaultFontName Monospaced
skinparam defaultFontSize 16
skinparam objectArrowColor DarkGray
skinparam RoundCorner 15
left to right direction
"""

const EPILOGUE = """
hide members
hide circle

@enduml
"""

function emit_subtypes(io::IO, T)
  for S in subtypes_shown(T)
    println(io, "\"", type_label(T), "\" --> \"", type_label(S), "\"")
    emit_subtypes(io, S)
  end
end

# `extra` carries relations that are not subtyping, so cannot be discovered by
# reflection.
function diagram(roots::Vector; extra::String = "")
  io = IOBuffer()
  print(io, PREAMBLE)

  for T in roots
    println(io)
    emit_subtypes(io, T)
  end

  print(io, extra)
  print(io, EPILOGUE)

  return Diagram(:plantuml, String(take!(io)))
end

parents = diagram([AbstractAlgebra.Set])

# `SetMap` is a hierarchy of its own: its members are not elements, but the
# third parameter `S` of `Map{D, C, S, T}`.
elements = diagram([AbstractAlgebra.SetElem, AbstractAlgebra.SetMap];
                   extra = "\n\"Map{D, C, S, T}\" .. \"SetMap\"\n")


open(joinpath(@__DIR__, "src", "assets", "parents_diagram.svg"), "w") do io
  write(io, sprint(show, "image/svg+xml", parents))
end

open(joinpath(@__DIR__, "src", "assets", "elements_diagram.svg"), "w") do io
  write(io, sprint(show, "image/svg+xml", elements))
end

open(joinpath(@__DIR__, "src", "assets", "parents_diagram.pdf"), "w") do io
  write(io, sprint(show, "application/pdf", parents))
end

open(joinpath(@__DIR__, "src", "assets", "elements_diagram.pdf"), "w") do io
  write(io, sprint(show, "application/pdf", elements))
end
