module PrettyPrinting

using ..AbstractAlgebra

import ..AbstractAlgebra: MatrixElem
import ..AbstractAlgebra: NCRingElem
import ..AbstractAlgebra: RingElem

using Preferences: Preferences, @load_preference, @set_preferences!

import Base: displaysize
import Base: get
import Base: getindex
import Base: haskey
import Base: in
import Base: lock
import Base: pipe_reader
import Base: pipe_writer
import Base: print
import Base: show
import Base: unlock
import Base: write

export @enable_all_show_via_expressify
export @show_name
export @show_special
export @show_special_elem
export Dedent
export IOCustom
export Indent
export Lowercase
export LowercaseOff
export allow_unicode
export canonicalize
export expr_to_latex_string
export expr_to_string
export expressify
export extra_name
export get_current_module
export get_html_as_latex
export get_name
export get_syntactic_sign_abs
export indent_string!
export is_syntactic_one
export is_syntactic_zero
export is_terse
export is_unicode_allowed
export obj_to_latex_string
export obj_to_string
export obj_to_string_wrt_times
export pretty
export print_integer_string
export print_obj
export printer
export set_current_module
export set_html_as_latex
export set_name!
export show_obj
export show_via_expressify
export terse
export with_unicode


# printing is done with respect to the following precedences
# There is no point in using the julia values because we add our own ops
const prec_lowest      = 0
const prec_inf_Equal   = 5     # infix ==
const prec_inf_Plus    = 11    # infix a+b+c
const prec_inf_Minus   = 11    # infix a-b-c
const prec_inf_Times   = 13    # infix a*b*c
const prec_inf_Divide  = 13    # infix a/b/c
const prec_inf_DoubleDivide  = 14    # infix a//b//c
const prec_inf_CenterDot = 15  # non associative * with spaces, \cdot in latex
const prec_pre_Plus    = 20    # prefix +a not used
const prec_pre_Minus   = 21    # prefix -a
const prec_pre_Times   = 22    # prefix *a not used
const prec_pre_Divide  = 23    # prefix /b not used
const prec_inf_Power   = 30    # infix a^b
const prec_post_call   = 99    # a(b) i.e. whether a+b(c) is (a+b)(c) vs a+(b(c))
const prec_post_ref    = 100   # a[b]

const prec_post_FractionBox    = 50    # precedence for a/b in 2d form
const prec_post_SuperscriptBox = 51    # precedence for a^b in 2d form

################################################################################
#
#  Expression to string
#
################################################################################

function expr_to_string(@nospecialize(obj))
   return sprint(show_obj, MIME("text/plain"), obj)
end

function expr_to_latex_string(@nospecialize(obj))
   return sprint(show_obj, MIME("text/latex"), obj)
end

function obj_to_string(@nospecialize(obj); context = nothing)
   return sprint(show_via_expressify, MIME("text/plain"), obj, context = context)
end

# parenthesize obj as if it were a term in a product
function obj_to_string_wrt_times(@nospecialize(obj); context = nothing)
   io = IOBuffer()
   S = AbstractAlgebra.printer(io)
   print_obj(S, MIME("text/plain"),
             canonicalize(expressify(obj, context = context)),
             prec_inf_Times, prec_inf_Times)
   finish(S)
   return String(take!(io))
end

function obj_to_latex_string(@nospecialize(obj); context = nothing)
   return sprint(show_via_expressify, MIME("text/latex"), obj, context = context)
end

function Base.showable(mi::MIME"text/latex", x::Union{RingElem, NCRingElem, MatrixElem})
    return !AbstractAlgebra.is_ijulia_inited()
end

function Base.show(io::IO, mi::MIME"text/latex", x::Union{RingElem, NCRingElem, MatrixElem})
   show_via_expressify(io, mi, x)
end

function Base.showable(mi::MIME"text/html", x::Union{RingElem, NCRingElem, MatrixElem})
   return !AbstractAlgebra.is_ijulia_inited() || AbstractAlgebra.get_html_as_latex()
end

function Base.show(io::IO, mi::MIME"text/html", x::Union{RingElem, NCRingElem, MatrixElem})
   if AbstractAlgebra.get_html_as_latex()
      io = IOContext(io, :size_limit => 1000)
   end
   show_via_expressify(io, mi, x)
end

function show_via_expressify(io::IO, @nospecialize(obj); context = nothing)
   show_via_expressify(io::IO, MIME("text/plain"), obj, context = context)
end

function show_via_expressify(io::IO, mi::MIME, @nospecialize(obj); context = nothing)
   show_obj(io, mi, canonicalize(expressify(obj, context = context)))
end

# the low-level workhorse
function show_obj(io::IO, mi::MIME, obj)
   S = printer(io)
   print_obj(S, mi, obj, prec_lowest, prec_lowest)
   finish(S)
end

const _html_as_latex = Ref{Bool}(false)

@doc raw"""
    get_html_as_latex()

Returns whether MIME type `text/html` is printed as `text/latex`.
"""
get_html_as_latex() = _html_as_latex[]

@doc raw"""
    set_html_as_latex(fl::Bool)

Toggles whether MIME type `text/html` should be printed as `text/latex`. Note
that this is a global option. The return value is the old value.
"""
function set_html_as_latex(fl::Bool)
  old = get_html_as_latex()
  _html_as_latex[] = fl
  return old
end

function show_obj(io::IO, mi::MIME"text/html", obj)
   if _html_as_latex[]
      print(io, "\$")
      show_obj(io, MIME("text/latex"), obj)
      print(io, "\$")
   else
      #print(io, "<font face=\"mono\">")
      show_obj(io, MIME("text/plain"), obj)
      #print(io, "</font>")
   end
end

################################################################################
#
#  Expressify
#
################################################################################

# The definitive guide to how Expr's are printed is the print_obj function
# below. Here is a quick summary of the latex behaviour:
#  * Expr(:call, ://, a, b) will create a fraction box
#  * Expr(:call, :/, a, b) will not create a fraction box
#  * heads :vcat, :hcat, and :row can be used to make arrays
#  * Expr(:latex_form, a, b) uses b for latex output and a otherwise
#  * Expr(:matrix, m) is a hint that m is a matrix, i.e. enclosed in ()
#  * Expr(:text, ...) instructs string arguments to be wrapped in \text{}
#  * Symbol leaves themselves have some special transformations:
#        Symbol("a")          => a
#        Symbol("α")          => {\alpha}
#        Symbol("x1")         => \operatorname{x1}
#        Symbol("xy_1")       => \operatorname{xy}_{1}
#        Symbol("sin")        => \operatorname{sin}
#        Symbol("sin_cos")    => \operatorname{sin\_cos}
#        Symbol("sin_1")      => \operatorname{sin}_{1}
#        Symbol("sin_cos_1")  => \operatorname{sin\_cos}_{1}
#        Symbol("αaβb_1_2")   => \operatorname{{\alpha}a{\beta}b}_{1,2}
#  * various sequential constructs:
#        Expr(:call, a, b, c)    => a(b,c)
#        Expr(:ref, a, b, c)     => a[b,c]
#        Expr(:vcat, a, b)       => [a;b]       (actually vertical in latex)
#        Expr(:vect, a, b)       => [a,b]
#        Expr(:tuple, a, b)      => (a,b)
#        Expr(:list, a, b)       => {a,b}
#        Expr(:series, a, b)     => a,b
#        Expr(:sequence, a, b)   => ab
#        Expr(:row, a, b)        => a b
#        Expr(:hcat, a, b)       => a b


function expressify(@nospecialize(a); context = nothing)
   return sprint(print, a; context = context)::String
end

# Only when AbstractAlgebra.expressify(a::T; context = nothing) has been
# defined may enable_all_show_via_expressify be used.
# AA defines Base.show for "text/latex" and "text/html" for a general set
# of x, but for backward compatibility it is not defined for general x and
# "text/plain" or the mime-less version.
# Rationale: when neither Base.show nor AA.expressify is defined for T, then,
# since expressify calls Base.show for backward compatibility, a definition of
# Base.show in terms of expressify would give a stack overflow.
#
# We have to be careful for jupyter. When jupyter prints an object x through
# IJulia.jl, it collects sprint(show, m, x) for all m in
# [ "text/plain", "text/html", "text/latex" ... ]
# and then picks one according to some precedence. If there is a string for
# either "text/html" or "text/latex", it will *not* use "text/plain". But
# "text/plain" is the only thing which will trigger rendering using a monospace
# font (hence it will look like "code"). This explains why all julia base
# objects are printed using monospace: none of these objects has print methods
# for text/html or text/latex.
#
# The challenge for us is that we always want to define text/latex and
# text/html methods for sprint and friends. This has the disadvantage that if
# get_html_as_latex() == false, then our objects will print their ordinary
# string presentation, but since it is coming from text/html, it will be
# rendered using the "normal" font.
#
# Thus, when IJulia is initialized we make `showable` return false for
# text/latex (always) and text/html (unless get_html_as_latex()), so IJulia
# falls back to text/plain.
#
# The IJulia state is queried via AbstractAlgebra.is_ijulia_inited(); by
# default this is false, and the IJulia extension installs the real check.
#
# Super easy!

macro enable_all_show_via_expressify(T)
  return quote
    function Base.show(io::IO, x::$(esc(T)))
       AbstractAlgebra.show_via_expressify(io, x)
    end

    function Base.show(io::IO, mi::MIME"text/plain", x::$(esc(T)))
       AbstractAlgebra.show_via_expressify(io, mi, x)
    end

    function Base.showable(mi::MIME"text/latex", x::$(esc(T)))
       return !AbstractAlgebra.is_ijulia_inited()
     end

    function Base.show(io::IO, mi::MIME"text/latex", x::$(esc(T)))
       return AbstractAlgebra.show_via_expressify(io, mi, x)
    end

    function Base.showable(mi::MIME"text/html", x::$(esc(T)))
       return !AbstractAlgebra.is_ijulia_inited() || AbstractAlgebra.get_html_as_latex()
    end

    function Base.show(io::IO, mi::MIME"text/html", x::$(esc(T)))
       return AbstractAlgebra.show_via_expressify(io, mi, x)
    end
  end
end

################################################################################
#
#   Canonicalization
#
################################################################################

# Canonicalization performs the following transformations
#   a+(b+c) => a+b+c    sums are flattened
#   a*(b*c) => a*b*c    products are flattened
#   a*-b*c  => -a*b*c   unary minus can move freely through products
#   (-a)/b  => -(a/b)   unary minus can move through the numerator of a quotient
#   a-b     => a + -b   subtraction is turned into addition with unary minus
#   --a     => a
#   0*a     => 0
#   0+a     => a
#   1*a     => a

# since unary minus has a higher precedence than * and /, we maintain
# -a*b*c as (-a)*b*c and not as -(a*b*c) because the former prints without ()

# is obj a call to op with 1 or more arguments
function isaExprOp(@nospecialize(obj), op::Symbol)
   return isa(obj, Expr) &&
          length(obj.args) > 1 &&
          obj.head === :call &&
          obj.args[1] === op
end

# is obj a call to op with len arguments
function isaExprOp(@nospecialize(obj), op::Symbol, len::Int)
   return isa(obj, Expr) &&
          length(obj.args) == len + 1 &&
          obj.head === :call &&
          obj.args[1] === op
end

# syntactic zeros can be removed from sums and turn a product into 0

is_syntactic_zero(obj::Number) = iszero(obj)

is_syntactic_zero(obj) = false

# syntactic ones can be removed from products

is_syntactic_one(obj::Number) = isone(obj)

is_syntactic_one(obj) = false


function get_syntactic_sign_abs(obj::Number)
   return obj < 0 ? (-1, -obj) : (1, obj)
end

function get_syntactic_sign_abs(obj::Expr)
   if obj.head !== :call
      return (1, obj)
   end
   if length(obj.args) == 2 && obj.args[1] === :-
      # unary minus is negative
      (sgn, abs) = get_syntactic_sign_abs(obj.args[2])
      return (-sgn, obj.args[2])
   elseif length(obj.args) > 2 && obj.args[1] === :*
      # product is negative if first term is
      (sgn, abs) = get_syntactic_sign_abs(obj.args[2])
      if sgn > 0
         return (1, obj)
      else
         newobj = Expr(obj.head)
         newobj.args = copy(obj.args)
         newobj.args[2] = abs
         if is_syntactic_one(newobj.args[2])
            deleteat!(newobj.args, 2)
         end
         if length(newobj.args) == 2
            return (sgn, newobj.args[2])
         else
            return (sgn, newobj)
         end
      end
   elseif length(obj.args) == 3 && (obj.args[1] === :/ || obj.args[1] === ://)
      # quotient is negative if numerator is
      (sgn, abs) = get_syntactic_sign_abs(obj.args[2])
      if sgn > 0
         return (1, obj)
      else
         newobj = Expr(obj.head)
         newobj.args = copy(obj.args)
         newobj.args[2] = abs
         return (sgn, newobj)
      end
   else
      return (1, obj)
   end
end

function get_syntactic_sign_abs(obj)
   return (1, obj)
end


function syntactic_neg!(obj::Expr)
   if isaExprOp(obj, :*) || isaExprOp(obj, :/) || isaExprOp(obj, ://)
      obj.args[2] = syntactic_neg!(obj.args[2])
      return obj
   elseif isaExprOp(obj, :-, 1)
      return obj.args[2]
   else
      return Expr(:call, :-, obj)
   end
end

function syntactic_neg!(obj)
   return Expr(:call, :-, obj)
end


# not actually implemented recursively to avoid stack overflow
function flatten_recursive!(ans::Expr, obj::Expr, op::Symbol)
   stack = reverse(obj.args[2:end])
   while !isempty(stack)
      top = pop!(stack)
      if isaExprOp(top, op)
         for i in length(top.args):-1:2
            push!(stack, top.args[i])
         end
      else
         push!(ans.args, top)
      end
   end
end

# op(op(a,b),op(c,d)) => op(a,b,c,d) etc
function flatten_op(obj::Expr, op::Symbol)
   if !isaExprOp(obj, op)
      return obj
   end
   for i in 2:length(obj.args)
      if isaExprOp(obj.args[i], op)
         ans = Expr(:call, op)
         flatten_recursive!(ans, obj, op)
         return ans
      end
   end
   return obj
end

function canonicalizePlusFinal!(obj::Expr)
   @assert obj.head === :call && obj.args[1] === :+
   if length(obj.args) < 2
      return 0
   elseif length(obj.args) == 2
      return obj.args[2]
   end
   obj = flatten_op(obj, :+)
   for i in 3:length(obj.args)
      (sign, abs) = get_syntactic_sign_abs(obj.args[i])
      if sign < 0
         obj.args[i] = Expr(:call, :-, abs)
      else
         obj.args[i] = abs
      end
   end
   return obj
end

function canonicalizePlus(obj::Expr)
   @assert obj.head === :call && obj.args[1] === :+
   if length(obj.args) < 2
      return 0
   elseif length(obj.args) == 2
      return canonicalize(obj.args[2])
   end
   # this flatten is just to try to avoid stack overflows in canonicalize
   obj = flatten_op(obj, :+)
   newobj = Expr(:call, :+)
   for i in 2:length(obj.args)
      t = canonicalize(obj.args[i])
      if !is_syntactic_zero(t)
         push!(newobj.args, t)
      end
   end
   # now do the real flatten
   return canonicalizePlusFinal!(newobj)
end

function canonicalizeMinus(obj::Expr)
   @assert obj.head === :call && obj.args[1] === :-
   if length(obj.args) < 2
      return 0
   elseif length(obj.args) == 2
      return syntactic_neg!(canonicalize(obj.args[2]))
   end
   newobj = Expr(:call, :+)
   for i in 2:length(obj.args)
      t = canonicalize(obj.args[i])
      if !is_syntactic_zero(t)
         push!(newobj.args, i > 2 ? syntactic_neg!(t) : t)
      end
   end
   return canonicalizePlusFinal!(newobj)
end

function canonicalizeTimes(obj::Expr)
   op = obj.args[1]
   @assert obj.head === :call && (op === :* || op === :cdot)
   if length(obj.args) < 2
      return 1
   elseif length(obj.args) == 2
      return canonicalize(obj.args[2])
   end
   if op === :cdot
      return canonicalize_general_recursive(obj)
   end
   obj = flatten_op(obj, op)
   newobj = Expr(:call, op)
   newsign = 1
   for i in 2:length(obj.args)
      t = canonicalize(obj.args[i])
      if is_syntactic_zero(t)
         return 0
      else
         (sign, abs) = get_syntactic_sign_abs(t)
         newsign *= sign
         if !is_syntactic_one(abs)
             push!(newobj.args, abs)
         end
      end
   end
   if length(newobj.args) < 2
      newobj = 1
   elseif length(newobj.args) == 2
      newobj = newobj.args[2]
   else
      newobj = flatten_op(newobj, op)
   end
   if newsign < 0
      newobj = syntactic_neg!(newobj)
   end
   return newobj
end

function canonicalize_general_recursive(obj::Expr)
   newobj = Expr(obj.head)
   for i in obj.args
      push!(newobj.args, canonicalize(i))
   end
   return newobj
end

function canonicalize(obj::Expr)
   if obj.head === :call && !isempty(obj.args)
      if obj.args[1] === :+
         return canonicalizePlus(obj)
      elseif obj.args[1] === :-
         return canonicalizeMinus(obj)
      elseif obj.args[1] === :* || obj.args[1] === :cdot
         return canonicalizeTimes(obj)
      end
   end
   return canonicalize_general_recursive(obj)
end

#fallback
function canonicalize(obj)
   return obj
end

################################################################################
#
#   Printing
#
################################################################################

mutable struct printer{IOT <: IO}
   io::IOT
   array::Vector{String}

   # if terse_level is positive we print a+b instead of a + b
   terse_level::Int
   size_limit_stack::Vector{Int}  # >= 0 for loosely-defined limit before ...
                                  # < 0 for unrestricted output

   function printer(io::IO)
      terse_level = get(io, :terse_level, 0)
      size_limit = get(io, :size_limit, -1)
      return new{typeof(io)}(io, String[], terse_level, Int[size_limit])
   end
end

# TODO since the subexpressions are not changing much, cache the leaf_count
# so that subexpressions don't have bad quadratic behavior
function leaf_count(S::printer, obj::Expr)
   z = 0
   for i in obj.args
      z += leaf_count(S, i)
   end
   return z
end

function leaf_count(S::printer, obj)
   return 1
end

# size_limit is a rough limit on the number of leaves printed
function size_limit(S::printer)
   return S.size_limit_stack[end]
end

function set_size_limit(S::printer, l::Int)
   push!(S.size_limit_stack, l)
end

function restore_size_limit(S::printer)
   pop!(S.size_limit_stack)
end

# maintain the last few printed things for token (un)mashing purposes
function push(S::printer, s::String)
   while length(S.array) > 3
      print(S.io, popfirst!(S.array))
   end
   push!(S.array, s)
end

function push_left_parenthesis(S::printer, ::MIME)
   push(S, "(")
end

function push_right_parenthesis(S::printer, ::MIME)
   push(S, ")")
end

function push_left_parenthesis(S::printer, ::MIME"text/latex")
   push(S, "\\left(")
end

function push_right_parenthesis(S::printer, ::MIME"text/latex")
   push(S, "\\right)")
end

function push_left_bracket(S::printer, ::MIME)
   push(S, "[")
end

function push_right_bracket(S::printer, ::MIME)
   push(S, "]")
end

function push_left_bracket(S::printer, ::MIME"text/latex")
   push(S, "\\left[")
end

function push_right_bracket(S::printer, ::MIME"text/latex")
   push(S, "\\right]")
end

function push_left_curly(S::printer, ::MIME)
   push(S, "{")
end

function push_right_curly(S::printer, ::MIME)
   push(S, "}")
end

function push_left_curly(S::printer, ::MIME"text/latex")
   push(S, "\\left{")
end

function push_right_curly(S::printer, ::MIME"text/latex")
   push(S, "\\right}")
end

function push_elision(S, ::MIME)
   push(S, "...")
end

function push_elision(S, ::MIME"text/latex")
   push(S, "{\\ldots}")
end

function finish(S::printer)
   while !isempty(S.array)
      print(S.io, popfirst!(S.array))
   end
end

# determine the limits for the subexpressions of a given expression
# look at obj.args[offset + 1], ..., obj.args[offset + n]
function child_limits(S::printer, obj::Expr, off::Int, n::Int)
   l = max(1, size_limit(S))

   if n > l
      leaf_counts = Int[leaf_count(S, obj.args[off +
                                 (2*i <= l + 1 ? i : i + n - l)]) for i in 1:l]
      n = l
   else
      leaf_counts = Int[leaf_count(S, obj.args[off + i]) for i in 1:n]
   end

   total_leaf_count = sum(leaf_counts)
   a = Int[] # will be limits for the terms from the start
   b = Int[] # will be limits for the terms from the end
   abtotal = 0
   while length(a) + length(b) < n
      ai = 1 + length(a)
      bi = n - length(b)
      if length(a) < length(b) || (length(a) == length(b) &&
                                   leaf_counts[ai] <= leaf_counts[bi])
         t = max(min(3, leaf_counts[ai]), div(leaf_counts[ai]*l, total_leaf_count))
         if leaf_counts[ai] < l/32
            t = max(t, leaf_counts[ai])
         end
         push!(a, min(t, l))
      else
         t = max(min(3, leaf_counts[bi]), div(leaf_counts[bi]*l, total_leaf_count))
         if leaf_counts[bi] < l/32
            t = max(t, leaf_counts[bi])
         end
         push!(b, min(t, l))
      end
      abtotal += t
      abtotal < l || break
   end
   return a, b
end

function compare_op_string(mi::MIME, op)
   if op === :(==)
      return "="
   else
      return string(op)::String
   end
end

function compare_op_string(mi::MIME"text/latex", op)
   if op === :(==)
      return "="
   elseif op === :(<=)
      return "\\le"
   elseif op === :(>=)
      return "\\ge"
   elseif op === :(!=)
      return "\\neq"
   else
      return string(op)::String
   end
end

function print_comparison(S::printer, mi::MIME, obj::Expr,
                           left::Int, right::Int, prec::Int)
   n = length(obj.args)
   @assert isodd(n) && n > 1
   sep = S.terse_level > 0 ? "" : " "

   needp = prec <= left || prec <= right
   if needp
      left = right = prec_lowest
      push_left_parenthesis(S, mi)
   end

   print_obj(S, mi, obj.args[1], left, prec)
   for i in 2:2:n-3
      push(S, sep*compare_op_string(mi, obj.args[i])*sep)
      print_obj(S, mi, obj.args[i+1], prec, prec)
   end
   push(S, sep*compare_op_string(mi, obj.args[n-1])*sep)
   print_obj(S, mi, obj.args[n], prec, right)

   if needp
      push_right_parenthesis(S, mi)
   end
end


# dir > 0 left associate: +, -, *, /
# dir < 0 right associative: ^
# dir = 0 non associative:
function printGenericInfix(S::printer, mi::MIME, obj::Expr,
                           left::Int, right::Int,
                           op::String, prec::Int, dir::Int)
   n = length(obj.args)
   if n < 3
      print_call_or_ref(S, mi, obj, left, right)
      return
   end

   needp = prec <= left || prec <= right
   if needp
      left = right = prec_lowest
      push_left_parenthesis(S, mi)
   end

   if size_limit(S) < 0
      # printing with no restriction
      print_obj(S, mi, obj.args[2], left, prec - (dir > 0))
      for i in 3:(n - 1)
         push(S, op)
         print_obj(S, mi, obj.args[i], prec, prec)
      end
      push(S, op)
      print_obj(S, mi, obj.args[n], prec - (dir < 0), right)
   elseif size_limit(S) == 0
      # no space to print anything
      push_elision(S, mi)
   else
      n -= 1
      a, b = child_limits(S, obj, 1, n)
      wrote_elision = false
      for i in 1:n

         if i <= length(a)
            set_size_limit(S, a[i])
         elseif n - i + 1 <= length(b)
            set_size_limit(S, b[n - i + 1])
         else
            if !wrote_elision
               i == 1 || push(S, op)
               push_elision(S, mi)
            end
            wrote_elision = true
            continue
         end
         i == 1 || push(S, op)
         print_obj(S, mi, obj.args[i + 1],
                   i == 1 ? left : i == n ? prec - (dir < 0) : prec,
                   i == 1 ? prec - (dir > 0) : i == n ? right : prec)
         restore_size_limit(S)
      end
   end

   if needp
      push_right_parenthesis(S, mi)
   end
end

function printGenericPrefix(S::printer, mi::MIME, obj::Expr,
                            left::Int, right::Int,
                            op::String, prec::Int)
   @assert length(obj.args) == 2
   needp = prec <= right
   if needp
      left = right = prec_lowest
      push_left_parenthesis(S, mi)
   end
   push(S, op)
   print_obj(S, mi, obj.args[2], prec, right)
   if needp
      push_right_parenthesis(S, mi)
   end
end

function printPlusArg(S::printer, mi::MIME, obj::Expr, i::Int,
                      left::Int, right::Int, prec::Int)
   n = length(obj.args)
   if i == 1
      print_obj(S, mi, obj.args[2], left, prec - 1)
   else
      arg = obj.args[i + 1]
      if isaExprOp(arg, :-, 1)
         push(S, S.terse_level > 0 ? "-" : " - ")
         arg = arg.args[2]
         left_prec = prec_inf_Minus
      else
         push(S, S.terse_level > 0 ? "+" : " + ")
         left_prec = prec_inf_Plus
      end
      right_prec = i + 2 > n ? right :
                  isaExprOp(obj.args[i + 2], :-, 1) ? prec_inf_Minus :
                                                      prec_inf_Plus
      print_obj(S, mi, arg, left_prec, right_prec)
   end
end

# special override for plus which seeks the sign of its arguments
function printPlus(S::printer, mi::MIME, obj::Expr,
                   left::Int, right::Int)
   n = length(obj.args)
   @assert n > 0 && obj.head === :call && obj.args[1] === :+
   if n < 2
      print_call_or_ref(S, mi, obj, left, right)
      return
   elseif n == 2
      printGenericPrefix(S, mi, obj, left, right, "+", prec_pre_Plus)
      return
   end

   prec = prec_inf_Plus
   needp = prec <= left || prec <= right
   if needp
      left = right = prec_lowest
      push_left_parenthesis(S, mi)
   end

   n -= 1
   if size_limit(S) < 0
      for i in 1:n
         printPlusArg(S, mi, obj, i, left, right, prec)
      end
   elseif size_limit(S) == 0
      push_elision(S, mi)
   else
      a, b = child_limits(S, obj, 1, n)
      wrote_elision = false
      for i in 1:n
         if i <= length(a)
            set_size_limit(S, a[i])
         elseif n - i + 1 <= length(b)
            set_size_limit(S, b[n - i + 1])
         else
            if !wrote_elision
               i == 1 || push(S, S.terse_level > 0 ? "+" : " + ")
               push_elision(S, mi)
            end
            wrote_elision = true
            continue
         end
         printPlusArg(S, mi, obj, i, left, right, prec)
         restore_size_limit(S)
      end
   end

   if needp
      push_right_parenthesis(S, mi)
   end
end

# special override for unary minus
function printMinus(S::printer, mi::MIME, obj::Expr,
                    left::Int, right::Int)
   n = length(obj.args)
   @assert n > 0 && obj.head === :call && obj.args[1] === :-
   if n < 2
      print_call_or_ref(S, mi, obj, left, right)
   elseif n == 2
      printGenericPrefix(S, mi, obj, left, right, "-", prec_pre_Minus)
   else
      op = S.terse_level > 0 ? "-" : " - "
      printGenericInfix(S, mi, obj, left, right, op, prec_inf_Minus, 1)
   end
end

function printFraction(S::printer, mi::MIME"text/latex",
                       @nospecialize(num), @nospecialize(den),
                       left::Int, right::Int)
   prec = prec_post_FractionBox
   needp = prec <= left || prec <= right
   if needp
      left = right = prec_lowest
      push_left_parenthesis(S, mi)
   end
   push(S, "\\frac{")
   print_obj(S, mi, num, prec_lowest, prec_lowest)
   push(S, "}{")
   print_obj(S, mi, den, prec_lowest, prec_lowest)
   push(S, "}")
   if needp
      push_right_parenthesis(S, mi)
   end
end

function printDivides(S::printer, mi::MIME"text/latex", obj::Expr,
                      left::Int, right::Int)
   n = length(obj.args)
   @assert n > 0 && obj.head === :call && obj.args[1] === ://
   if n != 3
      printGenericInfix(S, mi, obj, left, right, "//", prec_inf_Divide, +1)
   else
      (sgn, abs) = get_syntactic_sign_abs(obj.args[2])
      if sgn < 0
         prec = prec_pre_Minus
         needp = prec <= right
         if needp
            left = right = prec_lowest
            push_left_parenthesis(S, mi)
         end
         push(S, "-")
         printFraction(S, mi, abs, obj.args[3], prec, right)
         if needp
            push_right_parenthesis(S, mi)
         end
      else
         printFraction(S, mi, abs, obj.args[3], left, right)
      end
   end
end

function printPower(S::printer, mi::MIME"text/latex", obj::Expr,
                    left::Int, right::Int)
   n = length(obj.args)
   @assert n > 0 && obj.head === :call && obj.args[1] === :(^)
   if n != 3
      printGenericInfix(S, mi, obj, left, right, "^", prec_inf_Power, +1)
   else
      prec = prec_post_SuperscriptBox
      needp = prec <= left || prec <= right
      if needp
         left = right = prec_lowest
         push_left_parenthesis(S, mi)
      end
      print_obj(S, mi, obj.args[2], left, prec)
      push(S, "^{")
      print_obj(S, mi, obj.args[3], prec_lowest, prec_lowest)
      push(S, "}")
      if needp
         push_right_parenthesis(S, mi)
      end
   end
end

function print_call_or_ref(S::printer, mi::MIME, obj::Expr,
                           left::Int, right::Int)
   n = length(obj.args)
   @assert n > 0 && (obj.head === :call || obj.head === :ref)
   prec = obj.head === :call ? prec_post_call : prec_post_ref

   needp = prec <= left
   if needp
      left = prec_lowest
      push_left_parenthesis(S, mi)
   end

   if size_limit(S) < 0
      print_obj(S, mi, obj.args[1], left, prec)
      obj.head === :call ? push_left_parenthesis(S, mi) : push_left_bracket(S, mi)
      for i in 2:n
         i == 2 || push(S, S.terse_level > 0 ? "," : ", ")
         print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
      end
      obj.head === :call ? push_right_parenthesis(S, mi) : push_right_bracket(S, mi)
   elseif size_limit(S) == 0
      push_elision(S, mi)
   else
      a, b = child_limits(S, obj, 0, n)
      if 1 <= length(a)
         set_size_limit(S, a[1])
         print_obj(S, mi, obj.args[1], left, prec)
         restore_size_limit(S)
      else
         push_elision(S, mi)
      end
      obj.head === :call ? push_left_parenthesis(S, mi) : push_left_bracket(S, mi)
      wrote_elision = false
      for i in 2:n
         if i <= length(a)
            set_size_limit(S, a[i])
         elseif n - i + 1 <= length(b)
            set_size_limit(S, b[n - i + 1])
         else
            if !wrote_elision
               i == 2 || push(S, S.terse_level > 0 ? "," : ", ")
               push_elision(S, mi)
            end
            wrote_elision = true
            continue
         end
         i == 2 || push(S, S.terse_level > 0 ? "," : ", ")
         print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
         restore_size_limit(S)
      end
      obj.head === :call ? push_right_parenthesis(S, mi) : push_right_bracket(S, mi)
   end

   if needp
      push_right_parenthesis(S, mi)
   end
end

function print_tuple_etc(S::printer, mi::MIME, obj::Expr, left::Int, right::Int)
   n = length(obj.args)

   needp = prec_lowest < left || prec_lowest < right
   sep = S.terse_level > 0 ? "," : ", "
   if obj.head === :vcat
      needp = false
      sep = "; "
   elseif obj.head === :vect
      needp = false
   elseif obj.head === :tuple
      needp = false
   elseif obj.head === :list
      needp = false
   elseif obj.head === :series
   elseif obj.head === :sequence
      sep = ""
   elseif obj.head === :row || obj.head === :hcat
      sep = " "
   else
      error("invalid head")
   end

   needp && push_left_parenthesis(S, mi)

   if obj.head === :vcat || obj.head === :vect
      push_left_bracket(S, mi)
   elseif obj.head === :list
      push_left_curly(S, mi)
   elseif obj.head === :tuple
      push_left_parenthesis(S, mi)
   elseif obj.head === :row || obj.head === :hcat
      @assert S.terse_level >= 0
      S.terse_level += 1
   end

   if size_limit(S) < 0
      for i in 1:n
         i == 1 || push(S, sep)
         print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
      end
   elseif size_limit(S) == 0
      push_elision(S, mi)
   else
      a, b = child_limits(S, obj, 0, n)
      wrote_elision = false
      for i in 1:n
         if i <= length(a)
            set_size_limit(S, a[i])
         elseif n - i + 1 <= length(b)
            set_size_limit(S, b[n - i + 1])
         else
            if !wrote_elision
               i == 1 || push(S, sep)
               push_elision(S, mi)
            end
            wrote_elision = true
            continue
         end
         i == 1 || push(S, sep)
         print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
         restore_size_limit(S)
      end
   end

   if obj.head === :vcat || obj.head === :vect
      push_right_bracket(S, mi)
   elseif obj.head === :list
      push_right_curly(S, mi)
   elseif obj.head === :tuple
      push_right_parenthesis(S, mi)
   elseif obj.head === :row || obj.head === :hcat
      S.terse_level -= 1
      @assert S.terse_level >= 0
   end

   needp && push_right_parenthesis(S, mi)
end

function print_vcat(S::printer, mi::MIME"text/latex", obj::Expr,
                    left::Int, right::Int)

   nrows = length(obj.args)
   if nrows < 1
      push(S, "\\text{empty}")
      return
   end

   l = size_limit(S)
   if l == 0
      push_elision(S, mi)
      return
   end

   # print row i iff i <= b || i > nrows - b
   # print col j iff j <= a || j > ncols - a
   if l < 0
      # no limit on printing
      ncols = 1
      for i in 1:nrows
         ei = obj.args[i]
         if isa(ei, Expr) && (ei.head === :hcat || ei.head === :row)
            ncols = max(ncols, length(ei.args))
         end
      end
      b = nrows
      a = ncols
      entry_limit = l
   else
      ncols = 1
      nleaves = 1
      for i in 1:nrows
         ei = obj.args[i]
         if isa(ei, Expr) && (ei.head === :hcat || ei.head === :row)
            ncols = max(ncols, length(ei.args))
         end
         nleaves += leaf_count(S, ei)
      end
      s = sqrt((l + 1)/(nleaves + 1))
      a = max(1, round(Int, 0.5*s*ncols))
      b = max(1, round(Int, 0.5*s*nrows))
      entry_limit = cld(cld(l, min(2*a, ncols)), min(2*b, nrows))
   end

   # figure how many columns are actually going to be displayed
   if 2*a + 1 < ncols
      display_ncols = 2*a + 1
      use_col_elision = true
   else
      display_ncols = ncols
      use_col_elision = false
      a = ncols   # just take a to the max
   end

   push(S, "\\begin{array}")
   push(S, "{" * "c"^display_ncols * "}\n")

   set_size_limit(S, entry_limit)
   wrote_elision = false
   for i in 1:nrows
      ei = obj.args[i]
      if i <= b || i > nrows - b
         i == 1 || push(S, " \\\\\n")
         if isa(ei, Expr) && (ei.head === :hcat || ei.head === :row)
            eilen = length(ei.args)
            for j in 1:min(a, eilen)
               j == 1 || push(S, " & ")
               print_obj(S, mi, ei.args[j], prec_lowest, prec_lowest)
            end
            if use_col_elision && eilen > a
               push(S, " & \\cdots ")
               for j in (eilen - a + 1):eilen
                  push(S, " & ")
                  print_obj(S, mi, ei.args[j], prec_lowest, prec_lowest)
               end
            end
         else
            print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
         end
      elseif !wrote_elision
         i == 1 || push(S, " \\\\\n")
         if use_col_elision
            push(S, "\\vdots " * " & \\vdots "^(a - 1))
            push(S, " & \\ddots ")
            push(S, " & \\vdots "^a)
         else
            push(S, "\\vdots " * " & \\vdots "^(display_ncols - 1))
         end
         wrote_elision = true
      end
   end
   restore_size_limit(S)

   push(S, "\n\\end{array}")
end

# special case so that we don't have to negate integers
function print_integer_string(S::printer, mi::MIME, obj::String,
                              left::Int, right::Int)
   needp = prec_pre_Minus <= right && obj[1] == '-'
   if needp
      push_left_parenthesis(S, mi)
   end
   push(S, obj)
   if needp
      push_right_parenthesis(S, mi)
   end
end

function print_obj(S::printer, mi::MIME, obj::Integer,
                   left::Int, right::Int)
   print_integer_string(S, mi, string(obj), left, right)
end


const _latex_to_string = Dict{String, String}(
  "Α" => "\\alpha", "α" => "\\alpha", "Β" => "\\Beta", "β" => "\\beta", "Γ" =>
  "\\Gamma", "γ" => "\\gamma", "Δ" => "\\Delta", "δ" => "\\delta", "Ε" =>
  "\\Epsilon", "ε" => "\\epsilon", "Ζ" => "\\Zeta", "ζ" => "\\zeta", "Η" =>
  "\\Eta", "η" => "\\eta", "Θ" => "\\Theta", "θ" => "\\theta", "Ι" => "\\Iota",
  "ι" => "\\iota", "Κ" => "\\Kappa", "κ" => "\\kappa", "Λ" => "\\Lambda", "λ"
  => "\\lambda", "Μ" => "\\Mu", "μ" => "\\mu", "Ν" => "\\Nu", "ν" => "\\nu",
  "Ξ" => "\\Xi", "ξ" => "\\xi", "Ο" => "\\Omicron", "ο" => "\\omicron", "Π" =>
  "\\Pi", "π" => "\\pi", "Ρ" => "\\Rho", "ρ" => "\\rho", "Σ" => "\\Sigma",
  "σ" => "\\sigma", "Τ" => "\\Tau", "τ" => "\\tau", "Υ" => "\\Upsilon", "υ" =>
  "\\upsilon", "Φ" => "\\Phi", "φ" => "\\phi", "Χ" => "\\Chi", "χ" => "\\chi",
  "Ψ" => "\\Psi", "ψ" => "\\psi", "Ω" => "\\Omega", "omega" => "\\omega")

function deunicodify(x::String)
   z = ""
   for i in eachindex(x)
      c = string(x[i])
      if haskey(_latex_to_string, c)
         z *= "{" * _latex_to_string[c] * "}"
      else
         z *= c
      end
   end
   return z
end

function underscorify(x::String)
   if length(x) == 1 || occursin(r"\{|\}", x)
      return x
   end
   y = split(x, '_')
   n = length(y)
   while n > 1 && tryparse(Int, y[n]) != nothing
      n -= 1
   end
   if n == 1 && length(y[1]) == 1
      z = string(y[1])
   else
      # at this point we need operatorname and escaped underscores
      z = "\\mathop{\\mathrm{" * join(y[1:n], "\\_") * "}}"
   end
   if n < length(y)
      z = z * "_{" * join(y[n+1:end], ",") * "}"
   end
   return z
end

function print_obj(S::printer, ::MIME, obj::String,
                   left::Int, right::Int)
   push(S, obj)
end

function print_obj(S::printer, ::MIME, obj::Symbol,
                   left::Int, right::Int)
   push(S, string(obj))
end

function print_obj(S::printer, ::MIME"text/latex", obj::String,
                   left::Int, right::Int)
   push(S, deunicodify(obj))
end

function print_obj(S::printer, ::MIME"text/latex", obj::Symbol,
                   left::Int, right::Int)
   push(S, deunicodify(underscorify(string(obj))))
end

function print_obj(S::printer, mi::MIME, obj::Expr,
                   left::Int, right::Int)
   n = length(obj.args)
   if obj.head === :call && n > 0
      if obj.args[1] === :+
         printPlus(S, mi, obj, left, right)
      elseif obj.args[1] === :-
         printMinus(S, mi, obj, left, right)
      elseif obj.args[1] === :*
         printGenericInfix(S, mi, obj, left, right, "*", prec_inf_Times, +1)
      elseif obj.args[1] === :cdot
         printGenericInfix(S, mi, obj, left, right, " * ", prec_inf_CenterDot, 0)
      elseif obj.args[1] === :/
         printGenericInfix(S, mi, obj, left, right, "/", prec_inf_Divide, +1)
      elseif obj.args[1] === ://
         printGenericInfix(S, mi, obj, left, right, "//", prec_inf_DoubleDivide, +1)
      elseif obj.args[1] === :^
         printGenericInfix(S, mi, obj, left, right, "^", prec_inf_Power, -1)
      elseif obj.args[1] === :(==) || obj.args[1] === :(!=) ||
             obj.args[1] === :(>=) || obj.args[1] === :(>) ||
             obj.args[1] === :(<=) || obj.args[1] === :(<)
         o = compare_op_string(mi, obj.args[1])
         o = S.terse_level > 0 ? o : " "*o*" "
         printGenericInfix(S, mi, obj, left, right, o, prec_inf_Equal, 0)
      else
         print_call_or_ref(S, mi, obj, left, right)
      end
   elseif obj.head === :(=) && n == 2
      o = Expr(:call, :(=), obj.args[1], obj.args[2])
      printGenericInfix(S, mi, o, left, right, " = ", prec_inf_Equal, -1)
   elseif obj.head === :comparison && isodd(n) && n > 1
      print_comparison(S, mi, obj, left, right, prec_inf_Equal)
   elseif obj.head === :vcat || obj.head === :vect || obj.head === :tuple ||
          obj.head === :list || obj.head === :series || obj.head === :sequence ||
          obj.head === :row || obj.head === :hcat
      print_tuple_etc(S, mi, obj, left, right)
   elseif obj.head === :ref && n > 0
      print_call_or_ref(S, mi, obj, left, right)
   elseif obj.head === :text && n == 1
      print_obj(S, mi, obj.args[1], left, right)
   elseif obj.head === :latex_form && length(obj.args) >= 1
      print_obj(S, mi, obj.args[1], left, right)
   elseif obj.head === :matrix && n == 1
      print_obj(S, mi, obj.args[1], left, right)
   else
      push(S, "[??? unknown Expr ???]")
   end
end


function print_obj(S::printer, mi::MIME"text/latex", obj::Expr,
                   left::Int, right::Int)
   n = length(obj.args)
   if obj.head === :call && n > 0
      if obj.args[1] === :+
         printPlus(S, mi, obj, left, right)
      elseif obj.args[1] === :-
         printMinus(S, mi, obj, left, right)
      elseif obj.args[1] === :*
         printGenericInfix(S, mi, obj, left, right, " ", prec_inf_Times, +1)
      elseif obj.args[1] === :cdot
         printGenericInfix(S, mi, obj, left, right, " \\cdot ", prec_inf_CenterDot, 0)
      elseif obj.args[1] === :/
         printGenericInfix(S, mi, obj, left, right, "/", prec_inf_Divide, +1)
      elseif obj.args[1] === ://
         printDivides(S, mi, obj, left, right)
      elseif obj.args[1] === :^
         printPower(S, mi, obj, left, right)
      elseif obj.args[1] === :(==) || obj.args[1] === :(!=) ||
             obj.args[1] === :(>=) || obj.args[1] === :(>) ||
             obj.args[1] === :(<=) || obj.args[1] === :(<)
         o = compare_op_string(mi, obj.args[1])
         o = S.terse_level > 0 ? o : " "*o*" "
         printGenericInfix(S, mi, obj, left, right, o, prec_inf_Equal, 0)
      elseif obj.args[1] === :sqrt && length(obj.args) == 2
         needp = prec_inf_Power <= right # courtesy
         needp && push_left_parenthesis(S, mi)
         push(S, "\\sqrt{")
         print_obj(S, mi, obj.args[2], prec_lowest, prec_lowest)
         push(S, "}")
         needp && push_right_parenthesis(S, mi)
      else
         print_call_or_ref(S, mi, obj, left, right)
      end
   elseif obj.head === :(=) && n == 2
      o = Expr(:call, :(=), obj.args[1], obj.args[2])
      printGenericInfix(S, mi, o, left, right, " = ", prec_inf_Equal, -1)
   elseif obj.head === :comparison && isodd(n) && n > 1
      print_comparison(S, mi, obj, left, right, prec_inf_Equal)
   elseif obj.head === :vcat
      print_vcat(S, mi, obj, left, right)
   elseif obj.head === :vect || obj.head === :tuple ||
          obj.head === :list || obj.head === :series || obj.head === :sequence ||
          obj.head === :row || obj.head === :hcat
      print_tuple_etc(S, mi, obj, left, right)
   elseif obj.head === :ref && n > 0
      print_call_or_ref(S, mi, obj, left, right)
   elseif obj.head === :text && n == 1
      push(S, "\\text{")
      print_obj(S, mi, obj.args[1], left, right)
      push(S, "}")
   elseif obj.head === :latex_form && n >= 2
      print_obj(S, mi, obj.args[2], left, right)
   elseif obj.head === :matrix && n == 1
      push_left_parenthesis(S, mi)
      print_obj(S, mi, obj.args[1], left, right)
      push_right_parenthesis(S, mi)
   else
      push(S, "[??? unknown Expr ???]")
   end
end

function print_obj(S::printer, mi::MIME, obj, left::Int, right::Int)
   push(S, "[??? unknown object ???]")
end


###############################################################################
#
#  Macros for fancy printing.
#
###############################################################################

const CurrentModule = Ref(Main)

function set_current_module(m)
  CurrentModule[] = m
end

function get_current_module()
  return CurrentModule[]
end

"""
    set_name!(obj, name::String; override::Bool=true)

Sets the name of the object `obj` to `name`.
This name is used for printing using [`AbstractAlgebra.@show_name`](@ref).
If `override` is `false`, the name is only set if there is no name already set.

This function errors if `obj` does not support attribute storage.
"""
function set_name!(obj, name::String; override::Bool=true)
  override || isnothing(get_attribute(obj, :_name)) || return
  set_attribute!(obj, :_name => name)
end

"""
    set_name!(obj; override::Bool=true)

Sets the name of the object `obj` to the name of a variable in global (`Main` module) namespace
with value bound to the object `obj`, if such a variable exists (see [`AbstractAlgebra.PrettyPrinting.find_name`](@ref)).
This name is used for printing using [`AbstractAlgebra.@show_name`](@ref).
If `override` is `false`, the name is only set if there is no name already set.

This function errors if `obj` does not support attribute storage.
"""
function set_name!(obj; override::Bool=true)
  override || isnothing(get_attribute(obj, :_name)) || return
  sy = find_name(obj)
  isnothing(sy) && return
  set_name!(obj, string(sy); override=true)
end

"""
    extra_name(obj) -> Union{String,Nothing}

May be overloaded to provide a fallback name for the object `obj` in [`AbstractAlgebra.get_name`](@ref).
The default implementation returns `nothing`.
"""
extra_name(obj) = nothing

"""
    find_name(obj, M = Main; all::Bool = false) -> Union{String,Nothing}

Return name of a variable in `M`'s namespace with value bound to the object `obj`,
or `nothing` if no such variable exists.
If `all` is `true`, private and non-exported variables are also searched.

!!! note
    If the object is stored in several variables, the first one will be used,
    but a name returned once is kept until the variable no longer contains this object.

For this to work in doctests, one should call
`AbstractAlgebra.set_current_module(@__MODULE__)` in the `value` argument of
`Documenter.DocMeta.setdocmeta!` and keep the default value of `M = Main` here.

!!! warning
    This function should not be used directly, but rather through [`AbstractAlgebra.get_name`](@ref).
"""
function find_name(obj, M=Main; all::Bool=false)
  is_attribute_storing(obj) || return find_new_name(obj, M; all)

  cached_name = get_attribute(obj, :_cached_name)
  if !isnothing(cached_name)
   cached_name_sy = Symbol(cached_name)
    if M === Main && get_current_module() != Main
      if isdefined(get_current_module(), cached_name_sy) && getproperty(get_current_module(), cached_name_sy) === obj
        return cached_name
      end
    end
    if isdefined(M, cached_name_sy) && getproperty(M, cached_name_sy) === obj
      return cached_name
    end
  end
  name = find_new_name(obj, M; all)
  set_attribute!(obj, :_cached_name => name)
  return name
end

function find_new_name(obj, M=Main; all::Bool=false)
  # in Documenter, the examples are not run in the REPL.
  # in the REPL: A = ... adds A to the global name space (Main....)
  # in Documenter (doctests) all examples are run in their own module
  # which is stored in CurrentModule, hence we need to search there as well
  #furthermore, they are not exported, hence the "all" option
  if M === Main && get_current_module() != Main
    a = find_name(obj, get_current_module(), all=true)
    if !isnothing(a)
      return a
    end
  end
  for a = names(M; all=all)
    a === :ans && continue
    if isdefined(M, a) && getfield(M, a) === obj
      return string(a)
    end
  end
  return nothing
end

"""
    get_name(obj) -> Union{String,Nothing}

Returns the name of the object `obj` if it is set, or `nothing` otherwise.
This function tries to find a name in the following order:
1. The name set by [`AbstractAlgebra.set_name!`](@ref).
2. The name of a variable in global (`Main` module) namespace with value bound to the object `obj` (see [`AbstractAlgebra.PrettyPrinting.find_name`](@ref)).
3. The name returned by [`AbstractAlgebra.extra_name`](@ref).
"""
function get_name(obj)
  if is_attribute_storing(obj)
    name = get_attribute(obj, :_name)
    isnothing(name) || return name
  end

  sy = find_name(obj)
  isnothing(sy) || return string(sy)

  name_maybe = extra_name(obj)::Union{String,Nothing}
  isnothing(name_maybe) || return name_maybe

  return nothing
end

"""
    @show_name(io::IO, obj)

If either `is_terse(io)` is true or property `:compact` is set to `true` for `io`
(see [`IOContext`](https://docs.julialang.org/en/v1/base/io-network/#Base.IOContext)),
print the name [`get_name(obj)`](@ref AbstractAlgebra.get_name) of the object `obj` to
the `io` stream, then return from the current scope. Otherwise, do nothing.

It is supposed to be used at the start of `show` methods as shown in the documentation.
"""
macro show_name(io, obj)
  return :(
    begin
      local i = $(esc(io))
      local o = $(esc(obj))
      if is_terse(i) || get(i, :compact, false)
        name = get_name(o)
        if !isnothing(name)
          if AbstractAlgebra.PrettyPrinting._supports_io_custom(i)
            print(i, LowercaseOff())
          end
          print(i, name)
          return
        end
      end
    end
  )
end

@doc raw"""
    @show_special(io::IO, obj)

If the `obj` has a `show` attribute, this gets called with `io` and `obj` and
returns from the current scope. Otherwise, does nothing.

If `obj` does not have attribute storage available, this macro does nothing.

It is supposed to be used at the start of `show` methods as shown in the documentation.

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> R = @polynomial_ring(QQ, :x; cached=false)
Univariate polynomial ring in x over rationals

julia> AbstractAlgebra.@show_special(stdout, R)

julia> set_attribute!(R, :show, (i,o) -> print(i, "=> The One True Ring <="))

julia> AbstractAlgebra.@show_special(stdout, R)
=> The One True Ring <=

julia> R   # show for R uses @show_special, so we can observe the effect directly
=> The One True Ring <=
```
"""
macro show_special(io, obj)
  return :(
    begin
      local i = $(esc(io))
      local o = $(esc(obj))
      if is_attribute_storing(o)
        s = get_attribute(o, :show)
        if s !== nothing
          s(i, o)
          return
        end
      end
    end
  )
end

@doc raw"""
    @show_special(io::IO, mime, obj)

If the `obj` has a `show` attribute, this gets called with `io`, `mime` and `obj` (if applicable)
and `io` and `obj` otherwise, and returns from the current scope. Otherwise, does nothing.

If `obj` does not have attribute storage available, this macro does nothing.

It is supposed to be used at the start of `show` methods as shown in the documentation.

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> R = @polynomial_ring(QQ, :x; cached=false)
Univariate polynomial ring in x over rationals

julia> AbstractAlgebra.@show_special(stdout, MIME"text/plain"(), R)

julia> myshow(i,o) = print(i, "=> The One True Ring <=");

julia> myshow(i,m,o) = print(i, "=> The One True Ring with mime type $m <=");

julia> set_attribute!(R, :show, myshow)

julia> AbstractAlgebra.@show_special(stdout, MIME"text/plain"(), R)
=> The One True Ring with mime type text/plain <=

julia> R   # show for R uses @show_special, so we can observe the effect directly
=> The One True Ring <=
```
"""
macro show_special(io, mime, obj)
  return :(
    begin
      local i = $(esc(io))
      local m = $(esc(mime))
      local o = $(esc(obj))
      if is_attribute_storing(o)
        s = get_attribute(o, :show)
        if s !== nothing
          if applicable(s, i, m, o)
            s(i, m, o)
          else
            s(i, o)
          end
          return
        end
      end
    end
  )
end

@doc raw"""
    @show_special_elem(io::IO, obj)

If the `parent` of `obj` has a `show_elem` attribute, this gets called with `io` and `obj` and
returns from the current scope. Otherwise, does nothing.

If `parent(obj)` does not have attribute storage available, this macro does nothing.

It is supposed to be used at the start of `show` methods as shown in the documentation.

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> R = @polynomial_ring(QQ, :x; cached=false)
Univariate polynomial ring in x over rationals

julia> AbstractAlgebra.@show_special_elem(stdout, x)

julia> set_attribute!(R, :show_elem, (i,o) -> print(i, "=> $o <="))

julia> AbstractAlgebra.@show_special_elem(stdout, x)
=> x <=

julia> x   # show for x does not uses @show_special_elem, so x prints as before
x
```
"""
macro show_special_elem(io, obj)
  return :(
    begin
      local i = $(esc(io))
      local o = $(esc(obj))
      local p = parent(o)
      if is_attribute_storing(p)
        s = get_attribute(p, :show_elem)
        if s !== nothing
          s(i, o)
          return
        end
      end
    end
  )
end

@doc raw"""
    @show_special_elem(io::IO, mime, obj)

If the `parent` of `obj` has a `show_elem` attribute, this gets called with `io`, `mime` and `obj` (if applicable)
and `io` and `obj` otherwise, and returns from the current scope. Otherwise, does nothing.

If `parent(obj)` does not have attribute storage available, this macro does nothing.

It is supposed to be used at the start of `show` methods as shown in the documentation.

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> R = @polynomial_ring(QQ, :x; cached=false)
Univariate polynomial ring in x over rationals

julia> AbstractAlgebra.@show_special_elem(stdout, MIME"text/plain"(), x)

julia> set_attribute!(R, :show_elem, (i,m,o) -> print(i, "=> $o with mime type $m <="))

julia> AbstractAlgebra.@show_special_elem(stdout, MIME"text/plain"(), x)
=> x with mime type text/plain <=

julia> x   # show for x does not uses @show_special_elem, so x prints as before
x
```
"""
macro show_special_elem(io, mime, obj)
  return :(
    begin
      local i = $(esc(io))
      local m = $(esc(mime))
      local o = $(esc(obj))
      local p = parent(o)
      if is_attribute_storing(p)
        s = get_attribute(p, :show_elem)
        if s !== nothing
          if applicable(s, i, m, o)
            s(i, m, o)
          else
            s(i, o)
          end
          return
        end
      end
    end
  )
end


################################################################################
#
#  Unicode printing
#
################################################################################

const ALLOW_UNICODE_OVERRIDE_VALUE = Ref{Union{Bool,Nothing}}(nothing)

@doc """
    allow_unicode(allowed::Bool; temporary::Bool=false) -> Bool

Set whether unicode characters are allowed in pretty printing and returns the
previous value.
If `temporary` is `true`, then the change is only active for the current worker and session.
Otherwise, the change is permanently saved in the preferences.
A permanent change will always override a temporary change.

This function may behave arbitrarily if called from within the scope of a
`with_unicode` do-block.

See also [`is_unicode_allowed`](@ref) and [`with_unicode`](@ref).
"""
function allow_unicode(allowed::Bool; temporary::Bool=false)
   if temporary
      old_allowed = is_unicode_allowed()
      ALLOW_UNICODE_OVERRIDE_VALUE[] = allowed
      return old_allowed
   else
      old_allowed = is_unicode_allowed()
      @set_preferences!("unicode" => allowed)
      ALLOW_UNICODE_OVERRIDE_VALUE[] = nothing
      return old_allowed
   end
end

@doc """
    is_unicode_allowed() -> Bool

Return whether unicode characters are allowed in pretty printing.

See also [`allow_unicode`](@ref) and [`with_unicode`](@ref).
"""
function is_unicode_allowed()
   override = ALLOW_UNICODE_OVERRIDE_VALUE[]
   !isnothing(override) && return override
   return @load_preference("unicode", default = false)::Bool
end

@doc """
    with_unicode(f::Function, allowed::Bool=true)

Temporarily set whether unicode characters are allowed in pretty printing
during the execution of `f`.
This is useful for e.g. running doctests independently on the user preference.

`with_unicode` is expected to be called in the following way:
```julia
with_unicode([allowed]) do
   # code that should be executed with unicode allowed/disallowed
end
```

See also [`allow_unicode`](@ref) and [`is_unicode_allowed`](@ref).
"""
function with_unicode(f::Function, allowed::Bool=true)
   previous = ALLOW_UNICODE_OVERRIDE_VALUE[]
   ALLOW_UNICODE_OVERRIDE_VALUE[] = allowed
   try
      f()
   finally
      @assert ALLOW_UNICODE_OVERRIDE_VALUE[] == allowed
      ALLOW_UNICODE_OVERRIDE_VALUE[] = previous
   end
end

################################################################################
#
#  IO context with indentation and lowercasing
#
################################################################################

# The following piece of code draws inspiration from
#   https://github.com/KristofferC/IOIndents.jl
# But we do the indentation differently (and more correctly for multiline
# printing)

"""
    Indent

When printed to an `IOCustom` object, increases the indentation level by one.

# Examples

```repl
julia> io = AbstractAlgebra.pretty(stdout);

julia> print(io, AbstractAlgebra.Indent(), "This is indented")
  This is indented
```
"""
struct Indent end

"""
    Dedent

When printed to an `IOCustom` object, decreases the indentation level by one.

# Examples

```repl
julia> io = AbstractAlgebra.pretty(stdout);

julia> print(io, AbstractAlgebra.Indent(), AbstractAlgebra.Dedent(), "This is indented")
This is indented
```
"""
struct Dedent end

"""
    Lowercase

When printed to an `IOCustom` object, the next letter printed will be lowercase.

# Examples

```repl
julia> io = AbstractAlgebra.pretty(stdout);

julia> print(io, AbstractAlgebra.Lowercase(), "Foo")
foo
```
"""
struct Lowercase end

"""
    LowercaseOff

When printed to an `IOCustom` object, the case of the next letter will not be
changed when printed.

# Examples

```repl
julia> io = AbstractAlgebra.pretty(stdout);

julia> print(io, AbstractAlgebra.Lowercase(), AbstractAlgebra.LowercaseOff(), "Foo")
Foo
```
"""
struct LowercaseOff end

Base.show(io::IO, ::Union{Lowercase, LowercaseOff, Indent, Dedent}) = nothing

mutable struct IOCustom{IO_t <: IO} <: Base.AbstractPipe
    io::IO_t
    indent_level::Int
    indented_line::Bool
    indent_str::String
    printed::Int
    lowercasefirst::Bool
    force_newlines::Bool

    function IOCustom{IO_t}(io::IO_t, indent_level::Int,
                            indented_line::Bool, indent_str::String,
                            printed::Int, lowercasefirst::Bool, force_newlines::Bool) where IO_t <: IO
        @assert(!(IO_t <: IOCustom))
        return new(io, indent_level, indented_line, indent_str, printed, lowercasefirst, force_newlines)
    end
end

_unwrap(io::IOCustom) = io
_unwrap(io::IOContext) = io.io

_supports_io_custom(io::IOCustom) = true
_supports_io_custom(io::IOContext) = _supports_io_custom(io.io)
_supports_io_custom(io::Any) = false

indent_string!(io::IO, str::String) = (_unwrap(io).indent_str = str; io)

IOCustom(io::IO, force_newlines = false) = IOCustom{typeof(io)}(io, 0, false, "  ", 0, false, force_newlines)

IOCustom(io::IOCustom, force_newlines = false) = begin io.force_newlines = force_newlines; io; end

in(key_value::Pair, io::IOCustom) = in(key_value, io.io)
haskey(io::IOCustom, key) = haskey(io.io, key)
getindex(io::IOCustom, key) = getindex(io.io, key)
get(io::IOCustom, key, default) = get(io.io, key, default)

function show(_io::IO, io::IOCustom)
    ioi = IOCustom(_io)
    print(ioi, "IOCustom:", Indent())
    print(ioi, "\nIO: "); show(ioi, io.io)
    print(ioi, "\nIndent string: \"", io.indent_str, "\"")
    print(ioi, "\nIndent: ", io.indent_level)
end

pipe_reader(io::IOCustom) = io.io
pipe_writer(io::IOCustom) = io.io
lock(io::IOCustom) = lock(io.io)
unlock(io::IOCustom) = unlock(io.io)

Base.displaysize(io::IOCustom) = displaysize(io.io)

write(io::IO, ::Indent) = (_unwrap(io).indent_level += 1; nothing)
print(io::IO, ::Indent) = write(io, Indent())
write(io::IO, ::Dedent) = (_unwrap(io).indent_level = max(0, _unwrap(io).indent_level - 1); nothing)
print(io::IO, ::Dedent) = write(io, Dedent())
write(io::IO, ::Lowercase) = (_unwrap(io).lowercasefirst = true; nothing)
print(io::IO, ::Lowercase) = write(io, Lowercase())
write(io::IO, ::LowercaseOff) = (_unwrap(io).lowercasefirst = false; nothing)
print(io::IO, ::LowercaseOff) = write(io, LowercaseOff())

write_indent(io::IOCustom) = write(io.io, io.indent_str^io.indent_level)

write(io::IOCustom, chr::Char) = write(io, string(chr)) # Need to catch writing a '\n'

write(io::IOCustom, chr::Symbol) = write(io, string(chr))

_isbuffer(io::IOBuffer) = true
_isbuffer(io::IO) = false
_isbuffer(io::IOContext) = _isbuffer(io.io)
_isbuffer(io::IOCustom) = _isbuffer(io.io)

function _write_line(io::IOCustom, str::Union{String,SubString{String}})
  written::Int = 0

  # We handle io.indent_level == 0 differently for the following reason
  # $(s) will call print(io, s)
  # but if print(io, s) does io = pretty(io), we would insert new lines
  # after some random width (80?), which is no desirable.
  # So we do not insert newlines at all.
  if io.indent_level == 0
    if io.lowercasefirst
      written += write(io.io, lowercasefirst(str))
      io.lowercasefirst = false
    else
      written += write(io.io, str)
    end
    return written
  end

  # If we are writing to an IOBuffer, don't insert
  # artificial newlines
  #
  # Main application are doctests, since they are
  # printed to an IOBuffer for comparisons
  c = _isbuffer(io) && !io.force_newlines ? typemax(Int) : displaysize(io)[2]
  ind = io.indent_level * textwidth(io.indent_str)
  limit = c - ind > 0 ? c - ind : c
  # there might be already something written
  if c - ind - io.printed < 0
    spaceleft = mod(c - ind - io.printed, c)
  else
    spaceleft = c - ind - io.printed
  end
  #@show spaceleft
  if io.lowercasefirst
   str = lowercasefirst(str)
   io.lowercasefirst = false
  end
  # The following code deals with line wrapping of Unicode text, including
  # double-width symbols and more.
  _graphemes = Base.Unicode.graphemes(str)
  firstlen = min(spaceleft, length(_graphemes))
  # make an iterator over valid indices
  firstiter = Base.Iterators.take(_graphemes, firstlen)
  restiter = Base.Iterators.drop(_graphemes, firstlen)
  firststr = join(firstiter)
  width = textwidth(firststr)
  if length(firstiter) == width
    written += write(io.io, firststr)
    io.printed += width
  else
    #firstline is wider than number of graphemes
    partcollect = collect(firstiter)
    printstr = ""
    j = 1
    width = 0
    while width < limit && j <= length(partcollect)
      printstr *= partcollect[j]
      width += textwidth(partcollect[j])
      j += 1
    end
    written += write(io.io, printstr)
    io.printed += width

    #the spillover string
    written += write(io.io, "\n")
    written += write_indent(io)
    printstr = join(collect(firstiter)[j:end])
    written += write(io.io, printstr)
    io.printed = textwidth(printstr)
  end
  it = Iterators.partition(1:length(restiter), limit)
  restcollect = collect(restiter)
  for i in it
    # partitions of the spillover text
    partcollect = restcollect[i]
    partstr = join(partcollect)
    width = textwidth(partstr)
    if width < limit || length(i) == width
      written += write(io.io, "\n")
      written += write_indent(io)
      written += write(io.io, partstr)
      io.printed = width
    else
      # width is more than the number of graphemes
      # we can only ever get double length lines
      # (assuming non standard width can only be 2.)
      # (see https://github.com/alacritty/alacritty/issues/265#issue-199665364 )
      printstr = ""
      j = 1
      while textwidth(printstr) < (limit)
         printstr *= partcollect[j]
         j += 1
         if j > length(partcollect)
            break
         end
      end
      written += write(io.io, "\n")
      written += write_indent(io)
      written += write(io.io, printstr)
      io.printed = textwidth(printstr)
      # print the second part
      # there are at most two parts due to our assumption
      # that no grapheme exceeds double width
      printstr = join(partcollect[j:end])
      written += write(io.io, "\n")
      written += write_indent(io)
      written += write(io.io, printstr)
      io.printed = textwidth(printstr)
    end
  end
  return written
end

function write(io::IOCustom, str::String)
  written::Int = 0
  if str == "\n"
    written = write(io.io, str)
    io.indented_line = false
    io.printed = 0
    return written
  end

  for (i, line) in enumerate(split(str, "\n"))
    if i != 1
      written += write(io.io, "\n")
      io.printed = 0
      io.indented_line = false
    end

    if !io.indented_line
      written += write_indent(io)
      io.indented_line = true
    end
    written += _write_line(io, line)
  end
  if !isempty(str) && last(str) == '\n'
    io.indented_line = false
  else
    io.indented_line = true
  end
  return written
end

# Base.write on an IOContext does not call Base.write on the unwrapped context ...
Base.write(io::IOContext{<: IOCustom}, s::Union{SubString{String}, String}) = Base.write(Base.unwrapcontext(io)[1], s)

# println(io) redirects to print(io, '\n')
Base.write(io::IOContext{<: IOCustom}, s::Char) = Base.write(Base.unwrapcontext(io)[1], s)

Base.write(io::IOContext{<: IOCustom}, s::Symbol) = Base.write(Base.unwrapcontext(io)[1], s)

Base.take!(io::IOCustom) = take!(io.io)

"""
    pretty(io::IO) -> IOCustom

Wrap `io` into an `IOCustom` object.

# Examples

```repl
julia> io = AbstractAlgebra.pretty(stdout);
```
"""
pretty(io::IO; force_newlines = false) = IOCustom(io, force_newlines)

pretty(io::IOContext; force_newlines = false) = io.io isa IOCustom ? io : IOCustom(io, force_newlines)

# helpers for testing the pretty printing modes
repr_detailed(x) = repr(MIME"text/plain"(), x)
repr_oneline(x) = repr(x)
repr_terse(x) = repr(x, context = :supercompact => true)

# for backwards compatibility
# FIXME/TODO: remove these again
detailed(x) = repr_detailed(x)
oneline(x) = repr_oneline(x)
supercompact(x) = repr_terse(x)


#
"""
    terse(io::IO) -> IO

Return a new IO objects derived from `io` for which "terse" printing
mode has been enabled.

See <https://docs.oscar-system.org/stable/DeveloperDocumentation/printing_details/>
for details.

# Examples

```jldoctest
julia> AbstractAlgebra.is_terse(stdout)
false

julia> io = AbstractAlgebra.terse(stdout);

julia> AbstractAlgebra.is_terse(io)
true
```
"""
terse(io::IO) = IOContext(io, :supercompact => true)


"""
    is_terse(io::IO) -> Bool

Test whether "terse" printing mode is enabled for `io`.

See <https://docs.oscar-system.org/stable/DeveloperDocumentation/printing_details/>
for details.

# Examples

```jldoctest
julia> AbstractAlgebra.is_terse(stdout)
false

julia> io = AbstractAlgebra.terse(stdout);

julia> AbstractAlgebra.is_terse(io)
true
```
"""
is_terse(io::IO) = get(io, :supercompact, false)::Bool


end # PrettyPrinting
