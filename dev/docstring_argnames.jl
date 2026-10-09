#!/usr/bin/env julia
#
# Compare argument names in docstring signatures with those used in the code,
# and optionally rename the code to match.
#
#   julia dev/docstring_argnames.jl [--fix] [path ...]
#
# Paths default to src/. For every docstring whose first block lists
# signatures, each signature is matched (by function name, arity and type
# annotations) against the documented definition and against the undocumented
# definitions of the same function directly following it. Positional argument
# names that differ are reported.
#
# With --fix, definitions are rewritten to use the docstring's names, except
# when a new name already occurs somewhere in the definition. Those are
# reported and left for manual work. Only identifiers are touched; field
# names, keyword names in calls, quoted symbols and comments stay as they are.
#
#   julia dev/docstring_argnames.jl --rename M=A [--type Mat] [path ...]
#
# Renames every positional argument named `M` whose type annotation starts
# with `Mat` (or is a `Union{...}` mentioning it), together with its uses in
# the definition, the signature lines of its docstring, and `M` inside
# backtick or `$...$` spans of that docstring. In Markdown files, lines of
# code blocks that contain `M::Mat...` are rewritten. The same conflict rule
# applies: a definition or docstring already using `A` is skipped.
#
# Requires Julia >= 1.10 (uses Base.JuliaSyntax).

const JS = Base.JuliaSyntax
const SyntaxNode = JS.SyntaxNode

kind_(n) = JS.kind(n)
children(n) = (c = JS.children(n); c === nothing ? () : c)
text(n) = JS.sourcetext(n)
normtext(n) = replace(text(n), r"\s+" => "")
basename_(name) = String(last(split(name, '.')))

# ---- signatures --------------------------------------------------------------

struct Signature
   name::String
   args::Vector{Tuple{Union{Nothing,SyntaxNode},String}}   # (name node, type text)
end

argnames(s::Signature) = [a === nothing ? "" : text(a) for (a, _) in s.args]

# Strip `where`, return-type annotation and `->` wrappers around a call.
function unwrap_call(n)
   while true
      k = kind_(n)
      c = children(n)
      if k == JS.K"where" || k == JS.K"->"
         n = c[1]
      elseif k == JS.K"::" && length(c) == 2 && kind_(c[1]) in (JS.K"call", JS.K"where")
         n = c[1]
      else
         return n
      end
   end
end

# One positional argument: (name node or nothing, type text)
function argument(a)
   k = kind_(a)
   c = children(a)
   if k == JS.K"::"
      return length(c) == 2 ? (argument(c[1])[1], normtext(c[2])) : (nothing, normtext(c[1]))
   elseif k == JS.K"=" || k == JS.K"..."
      return argument(c[1])
   elseif k == JS.K"Identifier"
      return (a, "")
   end
   return (nothing, "")      # destructuring etc.
end

function signature(n)
   n = unwrap_call(n)
   kind_(n) == JS.K"call" || return nothing
   c = children(n)
   isempty(c) && return nothing

   args = Tuple{Union{Nothing,SyntaxNode},String}[]
   for a in c[2:end]
      kind_(a) == JS.K"parameters" && break
      push!(args, argument(a))
   end
   return Signature(text(c[1]), args)
end

# The signature node of a definition (function, short form, macro-wrapped)
function def_signature_node(n)
   k = kind_(n)
   c = children(n)
   if (k == JS.K"function" || k == JS.K"=") && !isempty(c)
      return c[1]
   elseif k == JS.K"macrocall" && length(c) >= 2
      return def_signature_node(c[end])
   end
   return nothing
end

is_definition(n) = (s = def_signature_node(n); s !== nothing && signature(s) !== nothing)

# Type annotations compared without parameters: `MatElem{T}` matches `MatElem`
strip_params(t) = replace(t, r"\{.*\}" => "")

# Same name and arity, and no contradicting type annotations. Returns the
# number of positions where both annotate the same type, and whether the
# docstring annotates every position the code annotates.
function match_score(doc::Signature, code::Signature)
   basename_(doc.name) == basename_(code.name) || return nothing
   length(doc.args) == length(code.args) || return nothing

   score = 0
   complete = true
   for ((_, td), (_, tc)) in zip(doc.args, code.args)
      isempty(tc) && continue
      if isempty(td)
         complete = false
         continue
      end
      strip_params(td) == strip_params(tc) || return nothing
      score += 1
   end
   return score, complete
end

# ---- docstrings --------------------------------------------------------------

# (docstring node, documented node) for `"""...""" def` and `@doc raw"""...""" def`
function doc_pair(n)
   k = kind_(n)
   c = children(n)
   k == JS.K"doc" && length(c) == 2 && return (c[1], c[2])
   k == JS.K"macrocall" && length(c) == 3 && text(c[1]) in ("doc", "@doc") && return (c[2], c[3])
   return nothing
end

# Value of the innermost string node
function string_value(n)
   if kind_(n) == JS.K"string"
      return join(string(c.val) for c in children(n) if c.val isa AbstractString)
   end
   for c in children(n)
      s = string_value(c)
      s === nothing || return s
   end
   return nothing
end

# Signatures in the indented block at the start of a docstring
function doc_signatures(doc)
   lines = String[]
   for l in split(doc, '\n')
      if startswith(l, "    ") && !isempty(strip(l))
         push!(lines, strip(l))
      elseif !isempty(lines)
         break
      end
   end

   sigs = Signature[]
   for l in lines
      s = try signature(JS.parsestmt(SyntaxNode, l)) catch; nothing end
      s === nothing || push!(sigs, s)
   end
   return sigs
end

# ---- scanning ----------------------------------------------------------------

# Best matching docstring signature for a code signature, or nothing. For an
# undocumented method following the docstring, the docstring signature must
# be at least as specific: `f(M::MatElem, rows, cols)` does not describe
# `f(M::MatElem, i::Int, cols::AbstractVector{Int})`.
function best_doc_signature(docsigs, code, documented::Bool)
   best = nothing
   bestscore = -1
   for ds in docsigs
      m = match_score(ds, code)
      m === nothing && continue
      score, complete = m
      (documented || complete) || continue
      score > bestscore && ((best, bestscore) = (ds, score))
   end
   return best
end

function scan!(out, node)
   sibs = children(node)
   for (i, n) in enumerate(sibs)
      if kind_(n) in (JS.K"module", JS.K"block", JS.K"toplevel")
         scan!(out, n)
         continue
      end

      dp = doc_pair(n)
      dp === nothing && continue
      docsigs = doc_signatures(something(string_value(dp[1]), ""))
      isempty(docsigs) && continue

      # the documented definition and the undocumented ones right after it
      defs = Any[dp[2]]
      for m in sibs[i+1:end]
         (doc_pair(m) === nothing && is_definition(m)) || break
         push!(defs, m)
      end

      for d in defs
         signode = def_signature_node(d)
         signode === nothing && continue
         code = signature(signode)
         code === nothing && continue
         doc = best_doc_signature(docsigs, code, d === dp[2])
         doc === nothing && continue

         renames = Pair{SyntaxNode,String}[]     # code name node => doc name
         for ((dn, _), (cn, _)) in zip(doc.args, code.args)
            (dn === nothing || cn === nothing) && continue
            text(dn) == text(cn) && continue
            push!(renames, cn => text(dn))
         end
         isempty(renames) && continue

         push!(out, (def = d, sig = signode, name = code.name, line = JS.source_line(signode),
                     docnames = argnames(doc), codenames = argnames(code), renames = renames))
      end
   end
end

# ---- renaming ----------------------------------------------------------------

function find_identifiers!(out, n, name)
   kind_(n) == JS.K"Identifier" && text(n) == name && push!(out, n)
   for c in children(n)
      find_identifiers!(out, c, name)
   end
end

const CALL_KINDS = (JS.K"call", JS.K"macrocall", JS.K"parameters")

# Identifier nodes named `name` that refer to the variable. Skipped: field
# names (`x.name`), quoted symbols, and keyword names in calls (`f(x; name = 1)`).
# Inside the signature `sig`, an `=` binds a default value, so its left side
# is renamed.
function rename_targets!(out, n, name, sig, in_sig, parent_kind)
   k = kind_(n)
   k == JS.K"quote" && return
   if k == JS.K"Identifier"
      text(n) == name && push!(out, n)
      return
   end

   in_sig = in_sig || n === sig
   is_keyword = k == JS.K"=" && !in_sig && parent_kind in CALL_KINDS
   for (i, ch) in enumerate(children(n))
      k == JS.K"." && i == 2 && continue
      is_keyword && i == 1 && continue
      rename_targets!(out, ch, name, sig, in_sig, k)
   end
end

# Byte edits for one mismatch, or the clashing name if a rename is unsafe
function rename_edits(m)
   edits = Pair{UnitRange{Int},String}[]
   for (node, newname) in m.renames
      clashes = SyntaxNode[]
      find_identifiers!(clashes, m.def, newname)
      isempty(clashes) || return nothing, newname

      targets = SyntaxNode[]
      rename_targets!(targets, m.def, text(node), m.sig, false, nothing)
      for t in targets
         push!(edits, JS.byte_range(t) => newname)
      end
   end
   return edits, nothing
end

function apply_edits(src, edits)
   bytes = Vector{UInt8}(codeunits(src))
   for (r, s) in sort(edits, by = first, rev = true)
      splice!(bytes, r, codeunits(s))
   end
   return String(bytes)
end

# ---- explicit rename mode ----------------------------------------------------

type_matches(t, prefix) = startswith(t, prefix) || (startswith(t, "Union{") && occursin(prefix, t))

# Inline code and math spans of a docstring line
const PROSE_SPANS = r"`[^`\n]*`|\$[^$\n]*\$"

# Rewrite the source text of a docstring: `old` in the leading signature
# lines, and inside prose spans. Returns (new text, prose conflict flag).
function rename_docstring(s, old, new)
   word = Regex("\\b" * old * "\\b")
   newword = Regex("\\b" * new * "\\b")
   lines = split(s, '\n')
   insigs = true
   conflict = false
   for i in 2:length(lines)
      l = lines[i]
      if insigs && startswith(l, "    ") && !isempty(strip(l))
         lines[i] = replace(l, word => new)
         continue
      end
      insigs = false
      conflict |= any(m -> occursin(newword, m.match), eachmatch(PROSE_SPANS, l))
   end
   conflict && return join(lines, '\n'), true
   for i in 2:length(lines)
      lines[i] = replace(lines[i], PROSE_SPANS => m -> replace(m, word => new))
   end
   return join(lines, '\n'), false
end

function rename_src!(f, old, new, prefix)
   src = read(f, String)
   tree = JS.parseall(SyntaxNode, src; filename = f, ignore_errors = true)
   edits = Pair{UnitRange{Int},String}[]
   nskipped = 0
   for (docnode, d) in definitions(tree)
      signode = def_signature_node(d)
      signode === nothing && continue
      code = signature(signode)
      code === nothing && continue
      any(((nm, ty),) -> nm !== nothing && text(nm) == old && type_matches(ty, prefix), code.args) || continue
      line = JS.source_line(signode)

      clashes = SyntaxNode[]
      find_identifiers!(clashes, d, new)
      if !isempty(clashes)
         println("$f:$line SKIPPED, `$new` already used in $(code.name)")
         nskipped += 1
         continue
      end

      targets = SyntaxNode[]
      rename_targets!(targets, d, old, signode, false, nothing)
      for t in targets
         push!(edits, JS.byte_range(t) => new)
      end
      docnode === nothing && continue

      doctext, conflict = rename_docstring(text(docnode), old, new)
      conflict && println("$f:$line docstring of $(code.name) mentions `$new`, prose not rewritten")
      doctext == text(docnode) || push!(edits, JS.byte_range(docnode) => doctext)
   end
   isempty(edits) || write(f, apply_edits(src, edits))
   return length(edits), nskipped
end

# (docstring node or nothing, definition node) for every top-level definition
function definitions(node, out = [])
   for n in children(node)
      if kind_(n) in (JS.K"module", JS.K"block", JS.K"toplevel")
         definitions(n, out)
         continue
      end
      dp = doc_pair(n)
      d = dp === nothing ? n : dp[2]
      is_definition(d) && push!(out, (dp === nothing ? nothing : dp[1], d))
   end
   return out
end

function rename_md!(f, old, new, prefix)
   word = Regex("\\b" * old * "\\b")
   sigline = Regex("\\b" * old * "::" * prefix)
   lines = readlines(f; keep = true)
   n = 0
   for i in eachindex(lines)
      occursin(sigline, lines[i]) || continue
      lines[i] = replace(lines[i], word => new)
      n += 1
   end
   n == 0 || write(f, join(lines))
   return n
end

function rename_mode(files, old, new, prefix)
   nedits = 0
   nskipped = 0
   for f in files
      if endswith(f, ".md")
         nedits += rename_md!(f, old, new, prefix)
         continue
      end
      e, s = rename_src!(f, old, new, prefix)
      nedits += e
      nskipped += s
   end
   println("$nedits edit(s), $nskipped definition(s) skipped")
   exit(0)
end

# ---- main ---------------------------------------------------------------------

function julia_files(paths, exts = (".jl",))
   files = String[]
   for p in paths
      if !isdir(p)
         push!(files, p)
         continue
      end
      for (root, _, fs) in walkdir(p), f in fs
         any(e -> endswith(f, e), exts) && push!(files, joinpath(root, f))
      end
   end
   return sort(files)
end

describe(m) = "$(m.name): doc ($(join(m.docnames, ", "))) vs code ($(join(m.codenames, ", ")))"

function option(args, name)
   i = findfirst(==(name), args)
   i === nothing && return nothing, args
   return args[i+1], vcat(args[1:i-1], args[i+2:end])
end

function main(args)
   rename, args = option(args, "--rename")
   prefix, args = option(args, "--type")
   fix = "--fix" in args
   paths = filter(a -> a != "--fix", args)
   isempty(paths) && (paths = ["src"])

   if rename !== nothing
      old, new = split(rename, "=")
      rename_mode(julia_files(paths, (".jl", ".md")), String(old), String(new), something(prefix, "Mat"))
   end

   nmismatch = 0
   nskipped = 0
   for f in julia_files(paths)
      src = read(f, String)
      tree = JS.parseall(SyntaxNode, src; filename = f, ignore_errors = true)
      mismatches = []
      scan!(mismatches, tree)
      nmismatch += length(mismatches)

      if !fix
         for m in mismatches
            println("$f:$(m.line) ", describe(m))
         end
         continue
      end

      edits = Pair{UnitRange{Int},String}[]
      for m in mismatches
         e, clash = rename_edits(m)
         if e === nothing
            println("$f:$(m.line) SKIPPED, `$clash` already used in ", describe(m))
            nskipped += 1
            continue
         end
         println("$f:$(m.line) $(m.name): renamed ", join(("$(text(a)) => $b" for (a, b) in m.renames), ", "))
         append!(edits, e)
      end
      isempty(edits) || write(f, apply_edits(src, edits))
   end

   fix ? println("$(nmismatch - nskipped) definition(s) renamed, $nskipped skipped") :
         println("$nmismatch mismatch(es)")
   exit(nmismatch == 0 ? 0 : 1)
end

main(ARGS)
