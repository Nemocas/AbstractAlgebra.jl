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

function obj_to_latex_string(@nospecialize(obj); context = nothing)
   return sprint(show_via_expressify, MIME("text/latex"), obj, context = context)
end

function Base.show(io::IO, mi::MIME"text/latex", x::Union{RingElem, NCRingElem, MatrixElem})
   show_via_expressify(io, mi, x)
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

global _html_as_latex = Ref(false)

@doc doc"""
    get_html_as_latex()

Returns whether MIME type `text/html` is printed as `text/latex`.
"""
get_html_as_latex() = _html_as_latex[]

@doc doc"""
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


# generic fallback in case expressify is not defined
function expressify(@nospecialize(a); context = nothing)::String
    s = sprint(print, a; context = context)::String
    if needs_parentheses(a)::Bool
        return "(" * s * ")"
    else
        return s
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
          obj.head == :call &&
          obj.args[1] == op
end

# is obj a call to op with len arguments
function isaExprOp(@nospecialize(obj), op::Symbol, len::Int)
   return isa(obj, Expr) &&
          length(obj.args) == len + 1 &&
          obj.head == :call &&
          obj.args[1] == op
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
   if length(obj.args) == 2 && obj.args[1] == :-
      # unary minus is negative
      (sgn, abs) = get_syntactic_sign_abs(obj.args[2])
      return (-sgn, obj.args[2])
   elseif length(obj.args) > 2 && obj.args[1] == :*
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
   elseif length(obj.args) == 3 && (obj.args[1] == :/ || obj.args[1] == ://)
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

# op(op(a,b),op(c,d)) => op(a,b,c,d) ect
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
   @assert obj.head == :call && obj.args[1] == :+
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
   @assert obj.head == :call && obj.args[1] == :+
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
   @assert obj.head == :call && obj.args[1] == :-
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
   @assert obj.head == :call && (op == :* || op == :cdot)
   if length(obj.args) < 2
      return 1
   elseif length(obj.args) == 2
      return canonicalize(obj.args[2])
   end
   if op == :cdot
      return canonicalizeOperation(obj)
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

function canonicalizeOperation(obj::Expr)
   @assert obj.head == :call
   newobj = Expr(:call, obj.args[1])
   for i in 2:length(obj.args)
      push!(newobj.args, canonicalize(obj.args[i]))
   end
   return newobj
end

function canonicalize(obj::Expr)
   if obj.head == :call && !isempty(obj.args)
      if obj.args[1] == :+
         return canonicalizePlus(obj)
      elseif obj.args[1] == :-
         return canonicalizeMinus(obj)
      elseif obj.args[1] == :* || obj.args[1] == :cdot
         return canonicalizeTimes(obj)
      elseif obj.args[1] == :/ || obj.args[1] == :// || obj.args[1] == :^
         return canonicalizeOperation(obj)
      else
         return obj
      end
   elseif obj.head == :vcat || obj.head == :vect || obj.head == :tuple ||
          obj.head == :list || obj.head == :series || obj.head == :sequence ||
          obj.head == :row || obj.head == :hcat || obj.head === :ref ||
          obj.head == :matrix || obj.head == :latex_form
      newobj = Expr(obj.head)
      for i in obj.args
         push!(newobj.args, canonicalize(i))
      end
      return newobj
   end
   return obj
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

# printing is done with respect to the following precedences
# There is no point in using the julia values because we add our own ops
prec_lowest      = 0
prec_inf_Plus    = 11    # infix a+b+c
prec_inf_Minus   = 11    # infix a-b-c
prec_inf_Times   = 13    # infix a*b*c
prec_inf_Divide  = 13    # infix a/b/c
prec_inf_DoubleDivide  = 14    # infix a//b//c
prec_inf_CenterDot = 15  # non associative * with spaces, \cdot in latex
prec_pre_Plus    = 20    # prefix +a not used
prec_pre_Minus   = 21    # prefix -a
prec_pre_Times   = 22    # prefix *a not used
prec_pre_Divide  = 23    # prefix /b not used
prec_inf_Power   = 30    # infix a^b
prec_post_call   = 99    # a(b) i.e. whether a+b(c) is (a+b)(c) vs a+(b(c))
prec_post_ref    = 100   # a[b]

prec_post_FractionBox    = 50    # precedence for a/b in 2d form
prec_post_SuperscriptBox = 51    # precedence for a^b in 2d form

mutable struct printer
   io::IO
   array::Vector{String}
   terse_level::Int
   size_limit_stack::Vector{Int}  # >= 0 for loosely-defined limit before ...
                                  # < 0 for unrestricted output
end

function printer(io::IO)
   terse_level = get(io, :terse, false) ? 1 : 0
   size_limit = get(io, :size_limit, -1)
   return printer(io, String[], terse_level, Int[size_limit])
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

# terse means a+b instead of a + b
function isterse(S::printer)
   return S.terse_level > 0
end

function set_terse(S::printer)
   S.terse_level = max(1, S.terse_level + 1)
end

function restore_terse(S::printer)
   S.terse_level = S.terse_level - 1
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

# dir > 0 left assiciate: +, -, *, /
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
         push(S, isterse(S) ? "-" : " - ")
         arg = arg.args[2]
         left_prec = prec_inf_Minus
      else
         push(S, isterse(S) ? "+" : " + ")
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
               i == 1 || push(S, isterse(S) ? "+" : " + ")
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
      op = isterse(S) ? "-" : " - "
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
   @assert n > 0 && (obj.head == :call || obj.head == :ref)
   prec = obj.head == :call ? prec_post_call : prec_post_ref

   needp = prec <= left
   if needp
      left = prec_lowest
      push_left_parenthesis(S, mi)
   end

   if size_limit(S) < 0
      print_obj(S, mi, obj.args[1], left, prec)
      obj.head == :call ? push_left_parenthesis(S, mi) : push_left_bracket(S, mi)
      for i in 2:n
         i == 2 || push(S, isterse(S) ? "," : ", ")
         print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
      end
      obj.head == :call ? push_right_parenthesis(S, mi) : push_right_bracket(S, mi)
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
      obj.head == :call ? push_left_parenthesis(S, mi) : push_left_bracket(S, mi)
      wrote_elision = false
      for i in 2:n
         if i <= length(a)
            set_size_limit(S, a[i])
         elseif n - i + 1 <= length(b)
            set_size_limit(S, b[n - i + 1])
         else
            if !wrote_elision
               i == 2 || push(S, isterse(S) ? "," : ", ")
               push_elision(S, mi)
            end
            wrote_elision = true
            continue
         end
         i == 2 || push(S, isterse(S) ? "," : ", ")
         print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
         restore_size_limit(S)
      end
      obj.head == :call ? push_right_parenthesis(S, mi) : push_right_bracket(S, mi)
   end

   if needp
      push_right_parenthesis(S, mi)
   end
end

function print_tuple_ect(S::printer, mi::MIME, obj::Expr, left::Int, right::Int)
   n = length(obj.args)

   needp = prec_lowest < left || prec_lowest < right
   sep = isterse(S) ? "," : ", "
   if obj.head == :vcat
      needp = false
      sep = "; "
   elseif obj.head == :vect
      needp = false
   elseif obj.head == :tuple
      needp = false
   elseif obj.head == :list
      needp = false
   elseif obj.head == :series
   elseif obj.head == :sequence
      sep = ""
   elseif obj.head == :row || obj.head == :hcat
      sep = " "
   else
      error("invalid head")
   end

   needp && push_left_parenthesis(S, mi)

   if obj.head == :vcat || obj.head == :vect
      push_left_bracket(S, mi)
   elseif obj.head == :list
      push_left_curly(S, mi)
   elseif obj.head == :tuple
      push_left_parenthesis(S, mi)
   elseif obj.head == :row || obj.head == :hcat
      set_terse(S)
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

   if obj.head == :vcat || obj.head == :vect
      push_right_bracket(S, mi)
   elseif obj.head == :list
      push_right_curly(S, mi)
   elseif obj.head == :tuple
      push_right_parenthesis(S, mi)
   elseif obj.head == :row || obj.head == :hcat
      restore_terse(S)
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
         if isa(ei, Expr) && (ei.head == :hcat || ei.head == :row)
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
         if isa(ei, Expr) && (ei.head == :hcat || ei.head == :row)
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
         if isa(ei, Expr) && (ei.head == :hcat || ei.head == :row)
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
  "\\Pie", "π" => "\\pie", "Ρ" => "\\Rho", "ρ" => "\\rho", "Σ" => "\\Sigma",
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
   # at this point we need operatorname and escaped underscores
   z = "\\operatorname{" * join(y[1:n], "\\_") * "}"
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
      else
         print_call_or_ref(S, mi, obj, left, right)
      end
   elseif obj.head == :vcat || obj.head == :vect || obj.head == :tuple ||
          obj.head == :list || obj.head == :series || obj.head == :sequence ||
          obj.head == :row || obj.head == :hcat
      print_tuple_ect(S, mi, obj, left, right)
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
   elseif obj.head === :vcat
      print_vcat(S, mi, obj, left, right)
   elseif obj.head == :vect || obj.head == :tuple ||
          obj.head == :list || obj.head == :series || obj.head == :sequence ||
          obj.head == :row || obj.head == :hcat
      print_tuple_ect(S, mi, obj, left, right)
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

