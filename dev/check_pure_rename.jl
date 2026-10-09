#!/usr/bin/env julia
#
# Check that the changes to Julia source files between two git revisions are
# nothing but consistent renamings of arguments and local variables.
#
#   julia dev/check_pure_rename.jl OLDREV NEWREV [file ...]
#
# NEWREV may be `-` for the working tree. Without file arguments, all .jl
# files under src/ that differ between the revisions are checked.
#
# Top-level method definitions are paired by name and argument types. Each
# changed pair is walked in lock-step to infer a symbol substitution. Reported
# as problems: anything that is not such a substitution, substitutions that
# merge two distinct variables, renamed keyword parameters (API change), and
# definitions or other top-level items that were added or removed.
# Exit status is 1 if problems were found.

isdoc(e) = e isa Expr && e.head == :macrocall &&
           e.args[1] in (GlobalRef(Core, Symbol("@doc")), Symbol("@doc"))

striplines(e) = e isa LineNumberNode ? nothing :
                e isa Expr ? Expr(e.head, map(striplines, e.args)...) : e

short(e) = (s = string(e); length(s) > 70 ? first(s, 67) * "..." : s)

# ---- top-level items ---------------------------------------------------------

# Top-level expressions of a file, docstrings and line numbers stripped,
# recursing into `module` blocks.
function toplevel_items(src)
   items = []
   collect_items!(items, Meta.parseall(src))
   return items
end

function collect_items!(items, e)
   for a in e.args
      a isa LineNumberNode && continue
      if a isa Expr && a.head == :module
         collect_items!(items, a.args[3])
         continue
      end
      isdoc(a) && (a = a.args[end])
      push!(items, striplines(a))
   end
end

# ---- definitions -------------------------------------------------------------

# A bare call at top level is a docstring stub such as `nrows(M::MatElem)`.
sigexpr(e) = e.head == :call ? e : e.args[1]

# The call part of a definition's signature (a Symbol for `function f end`)
function signature(e)
   s = sigexpr(e)
   while s isa Expr && s.head in (:where, :(::))
      s = s.args[1]
   end
   return s
end

function is_def(e)
   e isa Expr || return false
   e.head in (:function, :call) && return true
   e.head == :(=) || return false
   s = signature(e)
   return s isa Expr && s.head == :call
end

# `@inline function f ... end` etc.: the definition inside the macro call
defexpr(e) = e isa Expr && e.head == :macrocall && is_def(e.args[end]) ? e.args[end] : e
is_defitem(e) = is_def(defexpr(e))

# Replace argument names and default values by `_`, so that definitions can be
# paired across revisions by function name and argument types only.
function erase_names(a)
   a isa Symbol && return :_
   a isa Expr || return a
   a.head == :(::) && length(a.args) == 2 && return Expr(:(::), :_, a.args[2])
   a.head == :kw && return Expr(:kw, erase_names(a.args[1]), :_)
   a.head == :... && return Expr(:..., erase_names(a.args[1]))
   a.head == :parameters && return Expr(:parameters, map(erase_names, a.args)...)
   return a
end

function pairing_key(e)
   sig = deepcopy(sigexpr(e))
   c = sig
   while c isa Expr && c.head in (:where, :(::))
      c = c.args[1]
   end
   if c isa Expr && c.head == :call
      for i in 2:length(c.args)
         c.args[i] = erase_names(c.args[i])
      end
   end
   return string(sig)
end

# ---- rename inference --------------------------------------------------------

# Walk old and new expression in lock-step, recording in `subst` which symbol
# became which. `ctx` is :sig (positional parameters), :kwparams (keyword
# parameters) or :body.
function unify!(subst, issues, o, n, ctx)
   if o isa Symbol && n isa Symbol
      if get(subst, o, n) != n
         push!(issues, "inconsistent rename of `$o`: `$(subst[o])` and `$n`")
      end
      subst[o] = n
      return
   end

   if o isa Expr && n isa Expr
      if o.head != n.head || length(o.args) != length(n.args)
         push!(issues, "structural change: `$(short(o))` -> `$(short(n))`")
         return
      end

      # `f(x; key = val)`: the keyword name is not a variable
      if o.head == :kw && ctx == :body
         o.args[1] == n.args[1] ||
            push!(issues, "keyword name changed in call: `$(o.args[1])` -> `$(n.args[1])`")
         unify!(subst, issues, o.args[2], n.args[2], ctx)
         return
      end

      if o.head == :kw && ctx == :kwparams && o.args[1] != n.args[1]
         push!(issues, "keyword parameter renamed (API change): `$(o.args[1])` -> `$(n.args[1])`")
      end

      sub = o.head == :parameters && ctx == :sig ? :kwparams : ctx
      for (a, b) in zip(o.args, n.args)
         unify!(subst, issues, a, b, sub)
      end
      return
   end

   if o isa QuoteNode && n isa QuoteNode
      o.value == n.value || push!(issues, "quoted symbol changed: `$(o.value)` -> `$(n.value)`")
      return
   end

   o == n || push!(issues, "literal changed: `$(short(o))` -> `$(short(n))`")
end

function unify_item!(subst, issues, o, n)
   if o isa Expr && n isa Expr && o.head == n.head == :macrocall && length(o.args) == length(n.args)
      for i in 1:length(o.args)-1
         unify!(subst, issues, o.args[i], n.args[i], :body)
      end
      return unify_item!(subst, issues, o.args[end], n.args[end])
   end

   if is_def(o) && is_def(n) && o.head == n.head == :call
      unify!(subst, issues, o, n, :sig)
      return
   end

   if is_def(o) && is_def(n) && length(o.args) == length(n.args) == 2
      unify!(subst, issues, o.args[1], n.args[1], :sig)
      unify!(subst, issues, o.args[2], n.args[2], :body)
      return
   end

   unify!(subst, issues, o, n, :body)
end

# A rename `a -> b` is only safe if `b` was not already in use.
function check_merges!(issues, subst)
   renames = [o => n for (o, n) in subst if o != n]

   for (o, n) in renames
      get(subst, n, nothing) == n &&
         push!(issues, "`$o` -> `$n` merges with the existing name `$n`")
   end

   by_new = Dict{Symbol,Vector{Symbol}}()
   for (o, n) in renames
      push!(get!(by_new, n, Symbol[]), o)
   end
   for (n, os) in by_new
      length(os) > 1 && push!(issues, "`$(join(os, "`, `"))` all renamed to `$n`")
   end

   return sort(renames, by = string)
end

# ---- main ---------------------------------------------------------------------

readrev(rev, f) = rev == "-" ? read(f, String) : read(`git show $rev:$f`, String)

function changed_files(oldrev, newrev)
   cmd = newrev == "-" ? `git diff --name-only $oldrev -- src` :
                         `git diff --name-only $oldrev $newrev -- src`
   return filter(f -> endswith(f, ".jl"), split(read(cmd, String)))
end

function check_file(f, oldrev, newrev)
   old = toplevel_items(readrev(oldrev, f))
   new = toplevel_items(readrev(newrev, f))
   nproblems = 0

   queue = Dict{String,Vector{Any}}()
   for n in filter(is_defitem, new)
      push!(get!(queue, pairing_key(defexpr(n)), []), n)
   end

   for o in filter(is_defitem, old)
      q = get(queue, pairing_key(defexpr(o)), nothing)
      if q === nothing || isempty(q)
         println("$f: removed, or argument types changed:\n    ", short(sigexpr(defexpr(o))))
         nproblems += 1
         continue
      end
      n = popfirst!(q)
      o == n && continue

      subst = Dict{Symbol,Symbol}()
      issues = String[]
      unify_item!(subst, issues, o, n)
      renames = check_merges!(issues, subst)

      println("$f: ", short(sigexpr(defexpr(n))))
      isempty(renames) || println("    renamed: ", join(("$a => $b" for (a, b) in renames), ", "))
      for i in issues
         println("    !! ", i)
      end
      nproblems += !isempty(issues)
   end

   for n in reduce(vcat, values(queue); init = [])
      println("$f: added, or argument types changed:\n    ", short(sigexpr(defexpr(n))))
      nproblems += 1
   end

   others_old = filter(!is_defitem, old)
   others_new = filter(!is_defitem, new)
   for x in setdiff(others_old, others_new)
      println("$f: non-method top-level item changed or removed: ", short(x))
      nproblems += 1
   end
   for x in setdiff(others_new, others_old)
      println("$f: non-method top-level item changed or added: ", short(x))
      nproblems += 1
   end

   return nproblems
end

function main(args)
   if length(args) < 2
      println(stderr, "usage: julia dev/check_pure_rename.jl OLDREV NEWREV|- [file ...]")
      exit(2)
   end
   oldrev, newrev = args[1], args[2]
   files = length(args) > 2 ? args[3:end] : changed_files(oldrev, newrev)

   nproblems = sum(check_file(f, oldrev, newrev) for f in files; init = 0)
   println(nproblems == 0 ? "OK: only pure renames" : "$nproblems problem(s) found")
   exit(nproblems == 0 ? 0 : 1)
end

main(ARGS)
