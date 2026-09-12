###############################################################################
#
#   Type diagrams for the documentation
#
#   Regenerates docs/src/assets/{parents,elements}_diagram.{svg,pdf} from the
#   abstract types AbstractAlgebra actually defines:
#
#       julia --project=docs docs/create_type_diagrams.jl
#
#   Run it after adding, removing or renaming an abstract type.  Every box is
#   sized from its label, so a rename can never squeeze the text again.
#
###############################################################################

using AbstractAlgebra
using InteractiveUtils: subtypes

###############################################################################
#
#   Hierarchy
#
###############################################################################

struct TypeNode
  label::String
  children::Vector{TypeNode}
end

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

# Types deliberately left out.  `Ideal` is declared `<: Set`, yet the
# `parent_type` of a `Generic.Ideal` is an `IdealSet`, so it is neither a parent
# nor plainly an element.  Until that is settled it belongs in no diagram.
const EXCLUDED = ["Ideal{T}"]

# Only AbstractAlgebra's own abstract types.  `Generic`'s refinements and all
# concrete types are covered by the prose instead.
is_diagram_type(T) =
  isabstracttype(T) && parentmodule(T) === AbstractAlgebra && !(type_label(T) in EXCLUDED)

leaf_count(n::TypeNode) = isempty(n.children) ? 1 : sum(leaf_count, n.children)

# Types that read better at the bottom of their sibling block, because a
# separate hierarchy is drawn right underneath them.
const PLACE_LAST = ["Map{D, C, S, T}"]

function build_tree(T)
  children = [build_tree(S) for S in subtypes(T) if is_diagram_type(S)]

  # short branches first, so the long chains cascade towards the bottom right
  sort!(children; by = n -> (n.label in PLACE_LAST, leaf_count(n), n.label))

  return TypeNode(type_label(T), children)
end

###############################################################################
#
#   Layout
#
###############################################################################

const FONT_SIZE = 16.0
const CHAR_WIDTH = 0.6 * FONT_SIZE  # every monospace face, and PDF's Courier, is 0.6em wide

const BOX_HEIGHT = 30.0
const BOX_PADDING = 11.0            # left and right of the label
const CORNER_RADIUS = 8.0

const ROW_GAP = 12.0                # between two stacked boxes
const COLUMN_GAP = 46.0             # between a type and its subtypes
const TREE_GAP = 34.0               # between two independent hierarchies
const MARGIN = 4.0

const ARROW_LENGTH = 9.0
const ARROW_WIDTH = 7.0

struct Box
  label::String
  x::Float64                        # left edge
  y::Float64                        # top edge
  width::Float64
end

# `:subtype` is drawn as an arrow along the tree, `:uses` as a plain dashed
# line between two hierarchies.
struct Edge
  from::Int                         # index into Diagram.boxes
  to::Int
  kind::Symbol
end

struct Diagram
  boxes::Vector{Box}
  edges::Vector{Edge}
  width::Float64
  height::Float64
end

box_width(label::AbstractString) = CHAR_WIDTH * length(label) + 2 * BOX_PADDING

# A hierarchy to draw.  Trees are stacked in the given order; `align_under`
# names a box of an earlier tree whose indentation this one adopts.
struct Tree
  root::TypeNode
  align_under::Union{Nothing,String}
end

Tree(root::TypeNode) = Tree(root, nothing)

mutable struct Placement
  boxes::Vector{Box}
  index::Dict{String,Int}
  cursor::Float64
end

# Subtypes start right behind their supertype, so only a genuinely long chain of
# names widens the picture.  Leaves are stacked top to bottom in depth-first
# order; every other type is centred on the block of its subtypes.  Sibling
# subtrees occupy disjoint bands, so no two boxes can collide.
function place!(p::Placement, node::TypeNode, x::Float64)
  width = box_width(node.label)

  if isempty(node.children)
    y = p.cursor
    p.cursor += BOX_HEIGHT + ROW_GAP
  else
    ys = [place!(p, c, x + width + COLUMN_GAP) for c in node.children]
    y = (first(ys) + last(ys)) / 2
  end

  push!(p.boxes, Box(node.label, x, y, width))
  @assert !haskey(p.index, node.label) "duplicate label $(node.label)"
  p.index[node.label] = length(p.boxes)

  return y
end

function collect_edges!(edges::Vector{Edge}, index::Dict{String,Int}, node::TypeNode)
  for c in node.children
    push!(edges, Edge(index[node.label], index[c.label], :subtype))
    collect_edges!(edges, index, c)
  end
end

function build_diagram(trees::Vector{Tree}; links = Tuple{String,String}[])
  p = Placement(Box[], Dict{String,Int}(), MARGIN)

  for (i, t) in enumerate(trees)
    i > 1 && (p.cursor += TREE_GAP - ROW_GAP)
    x = t.align_under === nothing ? MARGIN : p.boxes[p.index[t.align_under]].x
    place!(p, t.root, x)
  end

  edges = Edge[]
  for t in trees
    collect_edges!(edges, p.index, t.root)
  end
  for (from, to) in links
    push!(edges, Edge(p.index[from], p.index[to], :uses))
  end

  width = maximum(b.x + b.width for b in p.boxes) + MARGIN
  height = p.cursor - ROW_GAP + MARGIN

  return Diagram(p.boxes, edges, width, height)
end

###############################################################################
#
#   Shared geometry of the two renderers
#
###############################################################################

# Trim the trailing `.0` that `string(::Float64)` insists on.
fmt(v::Real) = (r = round(v; digits = 2); isinteger(r) ? string(Int(r)) : string(r))

# A subtype arrow leaves the right edge of the supertype, shares one vertical
# trunk with its siblings, and enters the left edge of each subtype:
#
#     Ring ──┐
#            ├──> PolyRing{T}
#            └──> Field
#
# All subtypes of one supertype share a left edge, so they share the trunk too.
function subtype_elbow(from::Box, to::Box)
  x0 = from.x + from.width
  y0 = from.y + BOX_HEIGHT / 2
  x1 = to.x - ARROW_LENGTH
  y1 = to.y + BOX_HEIGHT / 2
  trunk = to.x - COLUMN_GAP / 2

  return [(x0, y0), (trunk, y0), (trunk, y1), (x1, y1)]
end

# A `:uses` link is a plain vertical line between two equally indented boxes.
function uses_line(from::Box, to::Box)
  x = from.x + min(from.width, to.width) / 2

  return [(x, from.y + BOX_HEIGHT), (x, to.y)]
end

edge_points(d::Diagram, e::Edge) =
  e.kind === :subtype ? subtype_elbow(d.boxes[e.from], d.boxes[e.to]) :
                        uses_line(d.boxes[e.from], d.boxes[e.to])

# Filled triangle at the head of a subtype arrow, pointing right.
function arrow_head(to::Box)
  tip = (to.x, to.y + BOX_HEIGHT / 2)

  return [tip,
          (tip[1] - ARROW_LENGTH, tip[2] - ARROW_WIDTH / 2),
          (tip[1] - ARROW_LENGTH, tip[2] + ARROW_WIDTH / 2)]
end

###############################################################################
#
#   SVG renderer
#
###############################################################################

const BOX_FILL = "#F8F8F8"
const BOX_STROKE = "#383838"
const TEXT_COLOR = "#000000"
const EDGE_COLOR = "#808080"
const FONT_STACK = "ui-monospace, SFMono-Regular, Menlo, Consolas, monospace"

points(pts) = join(("$(fmt(x)),$(fmt(y))" for (x, y) in pts), " ")

function render_svg(d::Diagram)
  io = IOBuffer()

  # The background stays transparent so the diagram sits on the page colour of
  # either documentation theme.
  println(io, """<?xml version="1.0" encoding="UTF-8"?>""")
  println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(fmt(d.width))" """ *
              """height="$(fmt(d.height))" viewBox="0 0 $(fmt(d.width)) $(fmt(d.height))">""")
  println(io, """<g fill="none" stroke="$EDGE_COLOR" stroke-width="1.2">""")

  for e in d.edges
    dash = e.kind === :uses ? """ stroke-dasharray="5 4\"""" : ""
    println(io, """<polyline points="$(points(edge_points(d, e)))"$dash/>""")
  end

  println(io, "</g>")
  println(io, """<g fill="$EDGE_COLOR" stroke="none">""")

  for e in d.edges
    e.kind === :subtype || continue
    println(io, """<polygon points="$(points(arrow_head(d.boxes[e.to])))"/>""")
  end

  println(io, "</g>")

  for b in d.boxes
    println(io, """<rect x="$(fmt(b.x))" y="$(fmt(b.y))" width="$(fmt(b.width))" """ *
                """height="$(fmt(BOX_HEIGHT))" rx="$(fmt(CORNER_RADIUS))" """ *
                """fill="$BOX_FILL" stroke="$BOX_STROKE" stroke-width="1.5"/>""")

    # `textLength` pins the label to the width the box was computed for, so the
    # two never drift apart no matter which monospace face the reader has.
    println(io, """<text x="$(fmt(b.x + BOX_PADDING))" """ *
                """y="$(fmt(b.y + BOX_HEIGHT / 2 + 0.3 * FONT_SIZE))" """ *
                """font-family="$FONT_STACK" font-size="$(fmt(FONT_SIZE))" """ *
                """textLength="$(fmt(b.width - 2 * BOX_PADDING))" lengthAdjust="spacing" """ *
                """fill="$TEXT_COLOR">$(b.label)</text>""")
  end

  println(io, "</svg>")

  return String(take!(io))
end

###############################################################################
#
#   PDF renderer
#
#   A hand-rolled one-page PDF avoids pulling a converter into the docs
#   toolchain.  Its Courier is metrically identical to what the SVG asks for
#   (0.6em per glyph), so both renderings agree on every box.
#
###############################################################################

# PDF's origin is the bottom left corner, the layout's is the top left one.
flip(d::Diagram, y::Real) = d.height - y

rgb(hex::String) = join((fmt(parse(Int, hex[i:i+1]; base = 16) / 255) for i in (2, 4, 6)), " ")

function pdf_polyline(io::IO, pts)
  (x, y) = first(pts)
  println(io, "$(fmt(x)) $(fmt(y)) m")

  for (x, y) in Iterators.drop(pts, 1)
    println(io, "$(fmt(x)) $(fmt(y)) l")
  end
end

function pdf_round_rect(io::IO, x, y, w, h, r)
  k = r * 0.5523                    # circle-to-bezier magic constant
  x0, y0, x1, y1 = x, y, x + w, y + h

  println(io, "$(fmt(x0 + r)) $(fmt(y0)) m")
  println(io, "$(fmt(x1 - r)) $(fmt(y0)) l")
  println(io, "$(fmt(x1 - r + k)) $(fmt(y0)) $(fmt(x1)) $(fmt(y0 + r - k)) $(fmt(x1)) $(fmt(y0 + r)) c")
  println(io, "$(fmt(x1)) $(fmt(y1 - r)) l")
  println(io, "$(fmt(x1)) $(fmt(y1 - r + k)) $(fmt(x1 - r + k)) $(fmt(y1)) $(fmt(x1 - r)) $(fmt(y1)) c")
  println(io, "$(fmt(x0 + r)) $(fmt(y1)) l")
  println(io, "$(fmt(x0 + r - k)) $(fmt(y1)) $(fmt(x0)) $(fmt(y1 - r + k)) $(fmt(x0)) $(fmt(y1 - r)) c")
  println(io, "$(fmt(x0)) $(fmt(y0 + r)) l")
  println(io, "$(fmt(x0)) $(fmt(y0 + r - k)) $(fmt(x0 + r - k)) $(fmt(y0)) $(fmt(x0 + r)) $(fmt(y0)) c")
  println(io, "h")
end

function pdf_content(d::Diagram)
  io = IOBuffer()

  println(io, "$(rgb(EDGE_COLOR)) RG $(rgb(EDGE_COLOR)) rg 1.2 w 1 J 1 j")

  for e in d.edges
    println(io, e.kind === :uses ? "[5 4] 0 d" : "[] 0 d")
    pdf_polyline(io, ((x, flip(d, y)) for (x, y) in edge_points(d, e)))
    println(io, "S")
  end

  println(io, "[] 0 d")

  for e in d.edges
    e.kind === :subtype || continue
    pdf_polyline(io, ((x, flip(d, y)) for (x, y) in arrow_head(d.boxes[e.to])))
    println(io, "f")
  end

  println(io, "$(rgb(BOX_STROKE)) RG $(rgb(BOX_FILL)) rg 1.5 w")

  for b in d.boxes
    pdf_round_rect(io, b.x, flip(d, b.y + BOX_HEIGHT), b.width, BOX_HEIGHT, CORNER_RADIUS)
    println(io, "B")
  end

  println(io, "BT /F1 $(fmt(FONT_SIZE)) Tf $(rgb(TEXT_COLOR)) rg")

  for b in d.boxes
    baseline = flip(d, b.y + BOX_HEIGHT / 2 + 0.3 * FONT_SIZE)
    println(io, "1 0 0 1 $(fmt(b.x + BOX_PADDING)) $(fmt(baseline)) Tm ($(b.label)) Tj")
  end

  println(io, "ET")

  return String(take!(io))
end

function render_pdf(d::Diagram)
  content = pdf_content(d)
  objects = [
    "<< /Type /Catalog /Pages 2 0 R >>",
    "<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
    "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 $(fmt(d.width)) $(fmt(d.height))] " *
      "/Resources << /Font << /F1 5 0 R >> >> /Contents 4 0 R >>",
    "<< /Length $(ncodeunits(content)) >>\nstream\n$(content)endstream",
    "<< /Type /Font /Subtype /Type1 /BaseFont /Courier >>",
  ]

  io = IOBuffer()
  print(io, "%PDF-1.4\n")

  offsets = Int[]
  for (i, obj) in enumerate(objects)
    push!(offsets, position(io))
    print(io, "$i 0 obj\n", obj, "\nendobj\n")
  end

  xref = position(io)
  print(io, "xref\n0 $(length(objects) + 1)\n0000000000 65535 f \n")
  for offset in offsets
    print(io, lpad(offset, 10, '0'), " 00000 n \n")
  end
  print(io, "trailer\n<< /Size $(length(objects) + 1) /Root 1 0 R >>\n")
  print(io, "startxref\n$xref\n%%EOF\n")

  return String(take!(io))
end

###############################################################################
#
#   Output
#
###############################################################################

function write_diagram(name::String, d::Diagram)
  dir = joinpath(@__DIR__, "src", "assets")

  write(joinpath(dir, "$(name).svg"), render_svg(d))
  write(joinpath(dir, "$(name).pdf"), render_pdf(d))

  println("$name: $(length(d.boxes)) types, $(fmt(d.width))x$(fmt(d.height))")
end

parents = build_diagram([Tree(build_tree(AbstractAlgebra.Set))])

# `SetMap` is a hierarchy of its own: its members are not elements but the
# third parameter of `Map{D, C, S, T}`, which the dashed line records.
elements = build_diagram([Tree(build_tree(AbstractAlgebra.SetElem)),
                          Tree(build_tree(AbstractAlgebra.SetMap), "Map{D, C, S, T}")];
                         links = [("Map{D, C, S, T}", "SetMap")])

write_diagram("parents_diagram", parents)
write_diagram("elements_diagram", elements)
