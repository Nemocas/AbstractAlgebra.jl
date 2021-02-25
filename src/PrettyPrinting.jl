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
   return sprint(show_via_expressify, MIME("text/plain"), obj)
end

function obj_to_latex_string(@nospecialize(obj); context = nothing)
   return sprint(show_via_expressify, MIME("text/latex"), obj)
end

function Base.show(io::IO, mi::MIME, x::RingElem)
   show_via_expressify(io, mi, x)
end

function show_via_expressify(io::IO, @nospecialize(obj))
   show_via_expressify(io::IO, MIME("text/plain"), obj)
end

function show_via_expressify(io::IO, mi::MIME, @nospecialize(obj))
   show_obj(io, mi, canonicalize(expressify(obj)))
end

# the low-level workhorse
function show_obj(io::IO, mi::MIME, obj)
   S = printer(io)
   print_obj(S, mi, obj, prec_lowest, prec_lowest)
   finish(S)
end

################################################################################
#
#  Expressify
#
################################################################################

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
   # julia supports only binary subtraction
   if length(obj.args) == 3 && isaExprOp(obj.args[3], :-, 1)
      obj.args[3] = obj.args[3].args[2]
      obj.args[1] = :-
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
   elseif obj.head == :vcat || obj.head == :hcat || obj.head == :row
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
prec_post_call   = 99   # a(b) i.e. whether a+b(c) is (a+b)(c) vs a+(b(c))

prec_post_FractionBox    = 50    # precedence for a/b in 2d form
prec_post_SuperscriptBox = 51    # precedence for a^b in 2d form

mutable struct printer
   io::IO
   array::Vector{String}
end

function printer(io::IO)
   return printer(io, String[])
end

function push(S::printer, s::String)
   # maintain the last few printed things for token (un)mashing purposes
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

function finish(S::printer)
   while !isempty(S.array)
      print(S.io, popfirst!(S.array))
   end
end

# dir > 0 left assiciate: +, -, *, /
# dir < 0 right associative: ^
# dir = 0 non associative:
function printGenericInfix(S::printer, mi::MIME, obj::Expr,
                           left::Int, right::Int,
                           op::String, prec::Int, dir::Int)
   n = length(obj.args)
   if n < 3
      printCall(S, mi, obj, left, right)
      return
   end
   needp = prec <= left || prec <= right
   if needp
      left = right = prec_lowest
      push_left_parenthesis(S, mi)
   end
   print_obj(S, mi, obj.args[2], left, prec - (dir > 0))
   for i in 3:(n - 1)
      push(S, op)
      print_obj(S, mi, obj.args[i], prec, prec)
   end
   push(S, op)
   print_obj(S, mi, obj.args[n], prec - (dir < 0), right)
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

# special override for plus which seeks the sign of its arguments 
function printPlus(S::printer, mi::MIME, obj::Expr,
                   left::Int, right::Int)
   n = length(obj.args)
   @assert n > 0 && obj.head === :call && obj.args[1] === :+
   if n < 2
      printCall(S, mi, obj, left, right)
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
   print_obj(S, mi, obj.args[2], left, prec - 1)
   for i in 3:n
      arg = obj.args[i]
      if isaExprOp(arg, :-, 1)
         push(S, " - ")
         arg = arg.args[2]
         left_prec = prec_inf_Minus
      else
         push(S, " + ")
         left_prec = prec_inf_Plus
      end
      right_prec = i + 1 > n ? right :
                  isaExprOp(obj.args[i + 1], :-, 1) ? prec_inf_Minus :
                                                      prec_inf_Plus
      print_obj(S, mi, arg, left_prec, right_prec)
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
      printCall(S, mi, obj, left, right)
   elseif n == 2
      printGenericPrefix(S, mi, obj, left, right, "-", prec_pre_Minus)
   else
      printGenericInfix(S, mi, obj, left, right, " - ", prec_inf_Minus, 1)
   end
end

function printCall(S::printer, mi::MIME, obj::Expr,
                   left::Int, right::Int)
   n = length(obj.args)
   @assert n > 0 && obj.head === :call
   prec = prec_post_call
   needp = prec <= left
   if needp
      left = prec_lowest
      push_left_parenthesis(S, mi)
   end
   print_obj(S, mi, obj.args[1], left, prec)
   push_left_parenthesis(S, mi)
   for i in 2:n
      print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
      if i < n
         push(S, ", ")
      end
   end
   push_right_parenthesis(S, mi)
   if needp
      push_right_parenthesis(S, mi)
   end
end


# What to do about / and //? Just print them both as /.
function printDivides(S::printer, mi::MIME"text/latex", obj::Expr,
                      left::Int, right::Int)
   n = length(obj.args)
   @assert n > 0 && obj.head === :call && (obj.args[1] === :/ || obj.args[1] === ://)
   if n != 3
      printGenericInfix(S, mi, obj, left, right, "/", prec_inf_Divide, +1)
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
      c = x[i]
      y = get(_latex_to_string, string(c), string(c)) 
      z = z * y
      if i < lastindex(x)
         z = z * " "
      end
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
   push(S, deunicodify(string(obj)))
end

function print_obj(S::printer, mi::MIME, obj::Expr,
                   left::Int, right::Int)
   if obj.head === :call && !isempty(obj.args)
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
         printCall(S, mi, obj, left, right)
      end
   elseif obj.head == :vcat
      push(S, "[")
      for i in 1:length(obj.args)
         if i > 1
            push(S, "; ")
         end
         print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
      end
      push(S, "]")
   elseif obj.head == :hcat || obj.head == :row
      for i in 1:length(obj.args)
         if i > 1
            push(S, " ")
         end
         print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
      end
   else
      push(S, "[??? unknown Expr ???]")
   end
end

function print_obj(S::printer, mi::MIME"text/latex", obj::Expr,
                   left::Int, right::Int)
   if obj.head === :call && !isempty(obj.args)
      if obj.args[1] === :+
         printPlus(S, mi, obj, left, right)
      elseif obj.args[1] === :-
         printMinus(S, mi, obj, left, right)
      elseif obj.args[1] === :*
         printGenericInfix(S, mi, obj, left, right, " ", prec_inf_Times, +1)
      elseif obj.args[1] === :cdot
         printGenericInfix(S, mi, obj, left, right, " \\cdot ", prec_inf_CenterDot, 0)
      elseif obj.args[1] === :/ || obj.args[1] === ://
         printDivides(S, mi, obj, left, right)
      elseif obj.args[1] === :^
         printPower(S, mi, obj, left, right)
      elseif obj.args[1] === :sqrt && length(obj.args) == 2
         push(S, "\\sqrt{")
         print_obj(S, mi, obj.args[2], prec_lowest, prec_lowest)
         push(S, "}")
      else
         printCall(S, mi, obj, left, right)
      end
   elseif obj.head == :vcat
      push(S, "\\begin{array}")
      ncols = 1
      for i in 1:length(obj.args)
         if isa(obj.args[i], Expr) && (obj.args[i].head == :hcat ||
                                       obj.args[i].head == :row)
            ncols = max(ncols, length(obj.args[i].args))
         end
      end
      push(S, "{" * "c"^ncols * "}\n")
      for i in 1:length(obj.args)
         if i > 1
             push(S, " \\\\\n")
         end
         print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
      end
      push(S, "\n\\end{array}\n")
   elseif obj.head == :hcat || obj.head == :row
      for i in 1:length(obj.args)
         if i > 1
            push(S, " & ")
         end
         print_obj(S, mi, obj.args[i], prec_lowest, prec_lowest)
      end
   else
      push(S, "[??? unknown Expr ???]")
   end
end

function print_obj(S::printer, mi::MIME, obj, left::Int, right::Int)
   push(S, "[??? unknown object ???]")
end

