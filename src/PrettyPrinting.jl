################################################################################
#
#  Expression to string
#
################################################################################

# Main function for Expr -> 1D String conversion
function expr_to_string(@nospecialize(obj))::String
    S = printer([])
    printExpr(S, obj, prec_lowest, prec_lowest)
    return getstring(S)
end

function expr_to_latex_string(@nospecialize(obj))::String
    S = printer([])
    printLatexExpr(S, obj, prec_lowest, prec_lowest)
    return getstring(S)
end

# Main function for Expr -> 2D String conversion
# currently pagewidth is unused because a proper implemention is too much work
function expr_to_2dstring(@nospecialize(obj), pagewidth::Int)::Vector{String}
    S = printer2D()
    printExpr(S, obj, prec_lowest, prec_lowest)
    rb = getrowbox(S)
    measure!(rb, pagewidth)
    chararr = [fill(' ', rb.sizex) for i in 1:rb.sizey]
    paint!(chararr, rb, 0, 0)
    return map(join, chararr)
end

# Shortcut for Obj -> 1D String conversion
function obj_to_string(@nospecialize(x); context = nothing)
  return expr_to_string(canonicalize(expressify(x, context = context)))
end

function obj_to_latex_string(@nospecialize(x); context = nothing)
  return expr_to_latex_string(canonicalize(expressify(x, context = context)))
end

function Base.show(io::IO, ::MIME"text/latex", x::RingElem)
  S = obj_to_latex_string(x)
  print(io, S)
end

################################################################################
#
#  Expressify
#
################################################################################

# The user should define expressify for all of their types.
# The conversion to expressions is generally clean and does not require
#   displayed_with_minus_in_front, needs_parentheses, ...

# Legacy API for no expressions
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
#   a-b     => a + -b   subtract is turned into addition with unary minus
#   --a     => a        
#   0*a     => 0
#   0+a     => a
#   1*a     => a 

# syntactic zeros can be removed from sums and turn a product into 0

is_syntactic_zero(obj::Number) = iszero(obj)

is_syntactic_zero(obj) = false

is_syntactic_one(obj::Number) = isone(obj)

is_syntactic_one(obj) = false

function get_syntactic_sign_abs(@nospecialize(obj))
    if isa(obj, Expr)
        if obj.head === :call
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
        else
            return (1, obj)
        end
    elseif isa(obj, Number)
        return obj < 0 ? (-1, -obj) : (1, obj)
    else
        return (1, obj)
    end
end

function syntactic_neg(@nospecialize(obj))
    (sgn, abs) = get_syntactic_sign_abs(obj)
    if sgn < 0
        return abs
    else
        return Expr(:call, :-, abs)
    end
end

# is obj a call to op with 1 or more arguments
function isaExprOp(@nospecialize(obj), op::Symbol)
    return isa(obj, Expr) &&
           length(obj.args) > 1 &&
           obj.head == :call &&
           obj.args[1] == op
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
    if length(newobj.args) < 2
        return 0
    elseif length(newobj.args) == 2
        return newobj.args[2]
    else
        return flatten_op(newobj, :+)
    end
end

function canonicalizeMinus(@nospecialize(obj))
    n = length(obj.args)
    @assert obj.head == :call && obj.args[1] == :-

    if n < 2
        return 0
    elseif n == 2
        return syntactic_neg(canonicalize(obj.args[2]))
    end

    newobj = Expr(:call, :+)
    for i in 2:n
        t = canonicalize(obj.args[i])
        if !is_syntactic_zero(t)
            push!(newobj.args, i > 2 ? syntactic_neg(t) : t)
        end
    end

    if length(newobj.args) < 2
        return 0
    elseif length(newobj.args) == 2
        return newobj.args[2]
    else
        return flatten_op(newobj, :+)
    end
end

function canonicalizeTimes(@nospecialize(obj))
    @assert obj.head == :call && obj.args[1] == :*

    if length(obj.args) < 2
        return 1
    elseif length(obj.args) == 2
        return canonicalize(obj.args[2])
    end

    obj = flatten_op(obj, :*)

    newobj = Expr(:call, :*)
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
        newobj = flatten_op(newobj, :*)
    end

    return newsign > 0 ? newobj : Expr(:call, :-, newobj)
end

function canonicalizeDivides(@nospecialize(obj))
    @assert obj.head == :call && (obj.args[1] == :/ || obj.args[1] == ://)
    newobj = Expr(:call, obj.args[1])
    # lets not do anything fancy with division
    for i in 2:length(obj.args)
        push!(newobj.args, canonicalize(obj.args[i]))
    end
    return newobj
end

function canonicalizePower(@nospecialize(obj))
    @assert obj.head == :call && obj.args[1] == :^
    newobj = Expr(:call, :^)
    for i in 2:length(obj.args)
        push!(newobj.args, canonicalize(obj.args[i]))
    end
    return newobj
end

function canonicalize(@nospecialize(obj))
    if !isa(obj, Expr)
        return obj
    end
    if obj.head == :call && !isempty(obj.args)
        if obj.args[1] == :+
            return canonicalizePlus(obj)
        elseif obj.args[1] == :-
            return canonicalizeMinus(obj)
        elseif obj.args[1] == :*
            return canonicalizeTimes(obj)
        elseif obj.args[1] == :/ || obj.args[1] == ://
            return canonicalizeDivides(obj)
        elseif obj.args[1] == :^
            return canonicalizePower(obj)
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

################################################################################
#
#   Printing
#
################################################################################

# printing is done with respect to the following precedences
prec_lowest      = 0
prec_inf_Plus    = 1    # infix a+b+c
prec_inf_Minus   = 2    # infix a-b-c
prec_pre_Plus    = 3    # prefix +a
prec_pre_Minus   = 4    # prefix -a
prec_inf_Times   = 5    # infix a*b*c
prec_inf_Divide  = 6    # infix a/b/c
prec_pre_Times   = 7    # prefix *a not used
prec_pre_Divide  = 8    # prefix /b not used
prec_inf_Power   = 9    # infix a^b

prec_post_SubsuperscriptBox = 10    # precedence for a^b in 2d form

prec_post_call   = 99   # a(b) i.e. whether a+b(c) is (a+b)(c) vs a+(b(c))

mutable struct printer
    array::Vector{String}
end

function push(S::printer, s::String)
    push!(S.array, s)
end

function getstring(S::printer)
    return join(S.array)
end

# prefix operator: unary minus -, ...
function printGenericPrefix(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n == 2
    needp = prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    push(S, op)
    printExpr(S, obj.args[2], prec, right)
    if needp
        push(S, ")")
    end
end

# infix operator: plus +, times *, ...
function printGenericInfix(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n > 2
    # actual condition is prec < left || prec < right, but because we flatten
    # + and *, we should see in the output if they are not flattened.
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[2], left, prec)
    for i in 3 : n - 1
        push(S, op)
        printExpr(S, obj.args[i], prec, prec)
    end
    push(S, op)
    printExpr(S, obj.args[n], prec, right)
    if needp
        push(S, ")")
    end
end

# left associative infix operator: ?
function printGenericInfixLeft(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n === 3
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[2], left, prec - 1)
    push(S, op)
    printExpr(S, obj.args[3], prec, right)
    if needp
        push(S, ")")
    end
end

# right associative infix operator: power ^, ...
function printGenericInfixRight(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n === 3
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[2], left, prec)
    push(S, op)
    printExpr(S, obj.args[3], prec - 1, right)
    if needp
        push(S, ")")
    end    
end

# non associative infix operator minus -, divides /
function printGenericInfixNone(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n > 2
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[2], left, prec - 1)
    for i in 3 : n - 1
        push(S, op)
        printExpr(S, obj.args[i], prec - 1, prec)
    end
    push(S, op)
    printExpr(S, obj.args[n], prec - 1, right)
    if needp
        push(S, ")")
    end
end

# postfix operator: factorial !, ...
function printGenericPostfix(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n == 2
    needp = prec <= left
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[2], prec, right)
    push(S, op)
    if needp
        push(S, ")")
    end
end

function printPlus(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(+)
    if n < 2
        push(S, "??? plus with zero arguments ???")
        return
    end
    prec = prec_inf_Plus
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[2], left, prec)
    for i in 3 : n
        sgn, abs = get_syntactic_sign_abs(obj.args[i])
        push(S, sgn > 0 ? " + " : " - ")
        printExpr(S, abs, prec, i < n ? prec : right)
    end
    if needp
        push(S, ")")
    end
end

function printMinus(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(-)
    if n < 2
        printcall(S, obj, left, right)
        return
    elseif n == 2
        printGenericPrefix(S, obj, left, right, "-", prec_pre_Minus)
        return
    end
    prec = prec_inf_Minus
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[2], left, prec)
    for i in 3 : n
        sgn, abs = get_syntactic_sign_abs(obj.args[i])
        push(S, sgn > 0 ? " - " : " + ")
        printExpr(S, abs, prec, i < n ? prec : right)
    end
    if needp
        push(S, ")")
    end
end



function printTimes(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(*)
    if n < 2
        printcall(S, obj, left, right)
    elseif n == 2
        printExpr(S, obj.args[2], left, right)
    else
        printGenericInfix(S, obj, left, right, "*", prec_inf_Times)
    end
end


function printDivides(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && (obj.args[1] === :/ || obj.args[1] === ://)
    op_str = (obj.args[1] === :/) ? "/" : "//"
    if n < 2
        printcall(S, obj, left, right)
    elseif n == 2
        printGenericPrefix(S, obj, left, right, op_str, prec_pre_Divide)
    elseif n == 3
        printGenericInfixLeft(S, obj, left, right, op_str, prec_inf_Divide)
    else
        printGenericInfixNone(S, obj, left, right, op_str, prec_inf_Divide)
    end
end

function printPower(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(^)
    if n < 3
        printcall(S, obj, left, right)
    elseif n == 3
        printGenericInfixRight(S, obj, left, right, "^", prec_inf_Power)
    else
        printGenericInfixNone(S, obj, left, right, "^", prec_inf_Power)
    end
end

function printcall(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call
    prec = prec_post_call
    needp = prec <= left
    if needp
        left = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[1], left, prec)
    push(S, "(")
    for i in 2:n
        printExpr(S, obj.args[i], prec_lowest, prec_lowest)
        if i < n
            push(S, ", ")
        end
    end
    push(S, ")")
    if needp
        push(S, ")")
    end
end

function printExpr(S::printer, obj, left::Int, right::Int)
    if isa(obj, String)
        push(S, obj)
    elseif isa(obj, Symbol)
        push(S, string(obj))
    elseif isa(obj, Integer)
        if obj < 0
            printGenericPrefix(S, Expr(:call, :-, string(-obj)), left, right, "-", prec_pre_Minus)
        else
            push(S, string(obj))
        end
    elseif isa(obj, Expr)
        if obj.head === :call && !isempty(obj.args)
            if obj.args[1] === :+
                printPlus(S, obj, left, right)
            elseif obj.args[1] === :-
                printMinus(S, obj, left, right)
            elseif obj.args[1] === :*
                printTimes(S, obj, left, right)
            elseif obj.args[1] === :/ || obj.args[1] === ://
                printDivides(S, obj, left, right)
            elseif obj.args[1] === :^
                printPower(S, obj, left, right)
            else
                printcall(S, obj, left, right)
            end
        elseif obj.head == :vcat
            push(S, "[")
            for i in 1:length(obj.args)
                if i > 1
                    push(S, "; ")
                end
                printExpr(S, obj.args[i], 0, 0)
            end
            push(S, "]")
        elseif obj.head == :hcat || obj.head == :row
            for i in 1:length(obj.args)
                if i > 1
                    push(S, " ")
                end
                printExpr(S, obj.args[i], 0, 0)
            end
        else
            push(S, "??? unknown Expr ???")
        end
    else
        push(S, "??? unknown ???")
    end
end

################################################################################
#
#  Latex printing
#
################################################################################

function printLatexGenericPrefix(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n == 2
    needp = prec <= right
    if needp
        left = right = prec_lowest
        push(S, "\\left(")
    end
    push(S, op)
    printLatexExpr(S, obj.args[2], prec, right)
    if needp
        push(S, "\\right)")
    end
end

# infix operator: plus +, times *, ...
function printLatexGenericInfix(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n > 2
    # actual condition is prec < left || prec < right, but because we flatten
    # + and *, we should see in the output if they are not flattened.
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "\\left(")
    end
    printLatexExpr(S, obj.args[2], left, prec)
    for i in 3 : n - 1
        push(S, op)
        printLatexExpr(S, obj.args[i], prec, prec)
    end
    push(S, op)
    printLatexExpr(S, obj.args[n], prec, right)
    if needp
        push(S, "\\right)")
    end
end

# left associative infix operator: ?
function printLatexGenericInfixLeft(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n === 3
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "\\left(")
    end
    printLatexExpr(S, obj.args[2], left, prec - 1)
    push(S, op)
    printLatexExpr(S, obj.args[3], prec, right)
    if needp
        push(S, "\\right)")
    end
end

# right associative infix operator: power ^, ...
function printLatexGenericInfixRight(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n === 3
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "\\left(")
    end
    printLatexExpr(S, obj.args[2], left, prec)
    push(S, op)
    printLatexExpr(S, obj.args[3], prec - 1, right)
    if needp
        push(S, "\\right)")
    end    
end

# non associative infix operator minus -, divides /
function printLatexGenericInfixNone(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n > 2
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "\\left(")
    end
    printLatexExpr(S, obj.args[2], left, prec - 1)
    for i in 3 : n - 1
        push(S, op)
        printLatexExpr(S, obj.args[i], prec - 1, prec)
    end
    push(S, op)
    printLatexExpr(S, obj.args[n], prec - 1, right)
    if needp
        push(S, "\\right)")
    end
end

# postfix operator: factorial !, ...
function printLatexGenericPostfix(S::printer, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n == 2
    needp = prec <= left
    if needp
        left = right = prec_lowest
        push(S, "\\right(")
    end
    printLatexExpr(S, obj.args[2], prec, right)
    push(S, op)
    if needp
        push(S, "\\left)")
    end
end

function printLatexPlus(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(+)
    if n < 2
        push(S, "??? plus with zero arguments ???")
        return
    end
    prec = prec_inf_Plus
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "\\left(")
    end
    printLatexExpr(S, obj.args[2], left, prec)
    for i in 3 : n
        sgn, abs = get_syntactic_sign_abs(obj.args[i])
        push(S, sgn > 0 ? " + " : " - ")
        printLatexExpr(S, abs, prec, i < n ? prec : right)
    end
    if needp
        push(S, "\\right)")
    end
end

function printLatexMinus(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(-)
    if n < 2
        printLatexcall(S, obj, left, right)
        return
    elseif n == 2
        printLatexGenericPrefix(S, obj, left, right, "-", prec_pre_Minus)
        return
    end
    prec = prec_inf_Minus
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "\\left(")
    end
    printLatexExpr(S, obj.args[2], left, prec)
    for i in 3 : n
        sgn, abs = get_syntactic_sign_abs(obj.args[i])
        push(S, sgn > 0 ? " - " : " + ")
        printLatexExpr(S, abs, prec, i < n ? prec : right)
    end
    if needp
        push(S, "\\right)")
    end
end



function printLatexTimes(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(*)
    if n < 2
        printLatexcall(S, obj, left, right)
    elseif n == 2
        printLatexExpr(S, obj.args[2], left, right)
    else
        printLatexGenericInfix(S, obj, left, right, "", prec_inf_Times)
    end
end


function printLatexDivides(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && (obj.args[1] === :/ || obj.args[1] === ://)
    op_str = (obj.args[1] === :/) ? "/" : "//"
    if n < 2
        printLatexcall(S, obj, left, right)
    elseif n == 2
        printLatexGenericPrefix(S, obj, left, right, op_str, prec_pre_Divide)
    elseif n == 3
      prec = prec_inf_Divide
      needp = prec <= left || prec <= right
      if needp
          left = right = prec_lowest
          push(S, "\\left(")
      end
      push(S, "\\frac{")
      printLatexExpr(S, obj.args[2], prec_lowest, prec_lowest)
      push(S, "}{")
      printLatexExpr(S, obj.args[3], prec_lowest, prec_lowest)
      push(S, "}")
      if needp
          push(S, "\\right)")
      end
    else
        printLatexGenericInfixNone(S, obj, left, right, op_str, prec_inf_Divide)
    end
end

function printLatexPower(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(^)
    if n < 3
        printLatexcall(S, obj, left, right)
    elseif n == 3
        prec = prec_inf_Power
        needp = prec <= left || prec <= right
        if needp
            left = right = prec_lowest
            push(S, "\\left(")
        end
        printLatexExpr(S, obj.args[2], left, prec)
        push(S, "^{")
        printLatexExpr(S, obj.args[3], prec_lowest, prec_lowest)
        push(S, "}")
        if needp
            push(S, "\\right)")
        end
    else
        printLatexGenericInfixNone(S, obj, left, right, "^", prec_inf_Power)
    end
end

function printLatexcall(S::printer, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call
    prec = prec_post_call
    needp = prec <= left
    if needp
        left = prec_lowest
        push(S, "(")
    end
    printLatexExpr(S, obj.args[1], left, prec)
    push(S, "(")
    for i in 2:n
        printLatexExpr(S, obj.args[i], prec_lowest, prec_lowest)
        if i < n
            push(S, ", ")
        end
    end
    push(S, ")")
    if needp
        push(S, ")")
    end
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

function printLatexExpr(S::printer, obj, left::Int, right::Int)
    if isa(obj, String)
        push(S, deunicodify(obj))
    elseif isa(obj, Symbol)
        push(S, deunicodify(string(obj)))
    elseif isa(obj, Integer)
        if obj < 0
            printLatexGenericPrefix(S, Expr(:call, :-, string(-obj)), left, right, "-", prec_pre_Minus)
        else
            push(S, string(obj))
        end
    elseif isa(obj, Expr)
        if obj.head === :call && !isempty(obj.args)
            if obj.args[1] === :+
                printLatexPlus(S, obj, left, right)
            elseif obj.args[1] === :-
                printLatexMinus(S, obj, left, right)
            elseif obj.args[1] === :*
                printLatexTimes(S, obj, left, right)
            elseif obj.args[1] === :/ || obj.args[1] === ://
                printLatexDivides(S, obj, left, right)
            elseif obj.args[1] === :^
                printLatexPower(S, obj, left, right)
            else
                printLatexcall(S, obj, left, right)
            end
        elseif obj.head == :vcat
            push(S, "[")
            for i in 1:length(obj.args)
                if i > 1
                    push(S, "; ")
                end
                printLatexExpr(S, obj.args[i], 0, 0)
            end
            push(S, "]")
        elseif obj.head == :hcat || obj.head == :row
            for i in 1:length(obj.args)
                if i > 1
                    push(S, " ")
                end
                printLatexExpr(S, obj.args[i], 0, 0)
            end
        else
            push(S, "??? unknown Expr ???")
        end
    else
        push(S, "??? unknown ???")
    end
end

################################################################################
#
#  2D printing
#
################################################################################

mutable struct ChildElem
    offx::Int
    offy::Int
    data::Any
end

function ChildElem(d)
    return ChildElem(-1, -1, d)
end

mutable struct StringBox
    sizex::Int
    sizey::Int
    basey::Int
    str::String
end

function StringBox(s::String)
    return StringBox(-1, -1, -1, s)
end

mutable struct RowBox
    sizex::Int
    sizey::Int
    basey::Int
    array::Vector{ChildElem}
end

function RowBox()
    return RowBox(-1, -1, -1, ChildElem[])
end

mutable struct BigBox
    sizex::Int
    sizey::Int
    basey::Int
    char::Char
end

function BigBox(c::Char)
    return BigBox(-1, -1, -1, c)
end

function BigBox(s::String)
    return BigBox(-1, -1, -1, s[1])
end

mutable struct SuperscriptBox
    sizex::Int
    sizey::Int
    basey::Int
    sup::ChildElem
end

function SuperscriptBox(d)
    return SuperscriptBox(-1, -1, -1, ChildElem(d))
end

mutable struct GridBox
    sizex::Int
    sizey::Int
    basey::Int
    array::Vector{Vector{ChildElem}}
end

function GridBox()
    return GridBox(-1, -1, -1, [])
end

mutable struct FractionBox
    sizex::Int
    sizey::Int
    basey::Int
    num::ChildElem
    den::ChildElem
end

function FractionBox(d1, d2)
    return FractionBox(-1, -1, -1, ChildElem(d1), ChildElem(d2))
end

mutable struct printer2D
    array::Vector{Any}
end

function printer2D()
    return printer2D([])
end

function push(S::printer2D, s::String)
    push!(S.array, StringBox(s))
end

function push(S::printer2D, d::Union{RowBox, BigBox, GridBox, SuperscriptBox, FractionBox})
    push!(S.array, d)
end

function getrowbox(S::printer2D)
    r = RowBox()
    for i in S.array
        push!(r.array, ChildElem(i))
    end
    return r
end

# prefix operator: unary minus -, ...
function printGenericPrefix(S::printer2D, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n == 2
    needp = prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    push(S, op)
    printExpr(S, obj.args[2], prec, right)
    if needp
        push(S, ")")
    end
end

# infix operator: plus +, times *, ...
function printGenericInfix(S::printer2D, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n > 2
    # actual condition is prec < left || prec < right, but because we flatten
    # + and *, we should see in the output if they are not flattened.
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[2], left, prec)
    for i in 3 : n - 1
        push(S, op)
        printExpr(S, obj.args[i], prec, prec)
    end
    push(S, op)
    printExpr(S, obj.args[n], prec, right)
    if needp
        push(S, ")")
    end
end

# left associative infix operator: ?
function printGenericInfixLeft(S::printer2D, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n === 3
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[2], left, prec - 1)
    push(S, op)
    printExpr(S, obj.args[3], prec, right)
    if needp
        push(S, ")")
    end
end

# right associative infix operator: power ^, ...
function printGenericInfixRight(S::printer2D, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n === 3
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, "(")
    end
    printExpr(S, obj.args[2], left, prec)
    push(S, op)
    printExpr(S, obj.args[3], prec - 1, right)
    if needp
        push(S, ")")
    end    
end

# non associative infix operator minus -, divides /
function printGenericInfixNone(S::printer2D, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n > 2
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, BigBox("("))
    end
    printExpr(S, obj.args[2], left, prec - 1)
    for i in 3 : n - 1
        push(S, op)
        printExpr(S, obj.args[i], prec - 1, prec)
    end
    push(S, op)
    printExpr(S, obj.args[n], prec - 1, right)
    if needp
        push(S, BigBox(")"))
    end
end

# postfix operator: factorial !, ...
function printGenericPostfix(S::printer2D, obj::Expr, left::Int, right::Int, op::String, prec::Int)
    n = length(obj.args)
    @assert n == 2
    needp = prec <= left
    if needp
        left = right = prec_lowest
        push(S, BigBox("("))
    end
    printExpr(S, obj.args[2], prec, right)
    push(S, op)
    if needp
        push(S, BigBox(")"))
    end
end


function printPlus(S::printer2D, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(+)
    if n < 2
        push(S, "??? plus with zero arguments ???")
        return
    end
    prec = prec_inf_Plus
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, BigBox("("))
    end
    printExpr(S, obj.args[2], left, prec)
    for i in 3 : n
        sgn, abs = get_syntactic_sign_abs(obj.args[i])
        push(S, sgn > 0 ? " + " : " - ")
        printExpr(S, abs, prec, i < n ? prec : right)
    end
    if needp
        push(S, BigBox(")"))
    end
end

function printMinus(S::printer2D, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(-)
    if n < 2
        push(S, "??? minus with zero arguments ???")
        return
    elseif n == 2
        printGenericPrefix(S, obj, left, right, "-", prec_pre_Minus)
        return
    end
    prec = prec_inf_Minus
    needp = prec <= left || prec <= right
    if needp
        left = right = prec_lowest
        push(S, BigBox("("))
    end
    printExpr(S, obj.args[2], left, prec)
    for i in 3 : n
        sgn, abs = get_syntactic_sign_abs(obj.args[i])
        push(S, sgn > 0 ? " - " : " + ")
        printExpr(S, abs, prec, i < n ? prec : right)
    end
    if needp
        push(S, BigBox(")"))
    end
end

function printTimes(S::printer2D, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(*)
    if n < 2
        push(S, "??? times with zero arguments ???")
    elseif n == 2
        printExpr(S, obj.args[2], left, right)
    else
        printGenericInfix(S, obj, left, right, "*", prec_inf_Times)
    end
end

function printDivides(S::printer2D, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && (obj.args[1] === :/ || obj.args[1] === ://)
    op_str = (obj.args[1] === :/) ? "/" : "//"
    if n < 2
        push(S, "??? divides with zero arguments ???")
    elseif n == 2
        printGenericPrefix(S, obj, left, right, op_str, prec_pre_Divide)
    elseif n == 3
        S2 = printer2D()
        printExpr(S2, obj.args[2], prec_lowest, prec_lowest)
        S3 = printer2D()
        printExpr(S3, obj.args[3], prec_lowest, prec_lowest)
        prec = prec_inf_Divide
        needp = prec <= left || prec <= right
        if needp
            push(S, BigBox("("))
        end
        push(S, FractionBox(getrowbox(S2), getrowbox(S3)))
        if needp
            push(S, BigBox(")"))
        end
    else
        printGenericInfixNone(S, obj, left, right, op_str, prec_inf_Divide)
    end
end

function printPower(S::printer2D, obj::Expr, left::Int, right::Int)
    n = length(obj.args)
    @assert n > 0 && obj.head === :call && obj.args[1] === :(^)
    if n < 2
        push(S, "??? power with zero arguments ???")
    elseif n == 2
        push(S, "??? power with one argument ???")
    elseif n == 3
        S2 = printer2D()
        printExpr(S2, obj.args[3], prec_lowest, prec_lowest)
        prec = prec_post_SubsuperscriptBox
        needp = prec <= left
        needp = needp || prec <= right  # extra parenthesis to be nice
        if needp
            left = right = prec_lowest
            push(S, BigBox("("))
        end
        printExpr(S, obj.args[2], left, prec)
        push(S, SuperscriptBox(getrowbox(S2)))
        if needp
            push(S, BigBox(")"))
        end
    else
        printGenericInfixNone(S, obj, left, right, "^", prec_inf_Power)
    end
end

#=
function printExpr(S::printer2D, obj::fmpz, left::Int, right::Int)
    if obj < 0
        printGenericPrefix(S, Expr(:call, :-, string(-obj)), left, right, "-", prec_pre_Minus)
    else
        push(S, string(obj))
    end
end
=#

function printExpr(S::printer2D, @nospecialize(obj), left::Int, right::Int)
    if isa(obj, String)
        push(S, obj)
    elseif isa(obj, Symbol)
        push(S, string(obj))
    elseif isa(obj, Integer)
        if obj < 0
            printGenericPrefix(S, Expr(:call, :-, string(-obj)), left, right, "-", prec_pre_Minus)
        else
            push(S, string(obj))
        end
    elseif isa(obj, Expr)
        if obj.head === :call && !isempty(obj.args)
            if obj.args[1] === :+
                printPlus(S, obj, left, right)
            elseif obj.args[1] === :-
                printMinus(S, obj, left, right)
            elseif obj.args[1] === :*
                printTimes(S, obj, left, right)
            elseif obj.args[1] === :/ || obj.args[1] === ://
                printDivides(S, obj, left, right)
            elseif obj.args[1] === :^
                printPower(S, obj, left, right)
            else
                push(S, "??? unknown Expr call ???")
            end
        elseif obj.head == :vcat
            mat = GridBox()
            for i in 1:length(obj.args)
                obji = obj.args[i]
                row = ChildElem[]
                if obji.head == :hcat || obji.head == :row
                    for j in 1:length(obji.args)
                        S2 = printer2D()
                        printExpr(S2, obji.args[j], prec_lowest, prec_lowest)
                        push!(row, ChildElem(getrowbox(S2)))
                    end                    
                else
                    S2 = printer2D()
                    printExpr(S2, obji, prec_lowest, prec_lowest)
                    push!(row, ChildElem(getrowbox(S2)))
                end
                push!(mat.array, row)
            end
            push(S, mat)
        elseif obj.head == :hcat || obj.head == :row
            mat = GridBox()
            obji = obj.args[i]
            row = ChildElem[]
            for j in 1:length(obj.args)
                S2 = printer2D()
                printExpr(S2, obj[j], prec_lowest, prec_lowest)
                    push!(row, ChildElem(getrowbox(S2)))
            end                    
            push!(mat.array, row)
            push(S, mat)
        else
            push(S, "??? unknown Expr ???")
        end
    else
        push(S, "??? unknown ???")
    end
end



function paint!(s::Vector{Vector{Char}}, a::RowBox, offx::Int, offy::Int)
    for i in a.array
        paint!(s, i.data, offx + i.offx, offy + i.offy)
    end
end

function measure!(a::RowBox, w::Int)
    width = 0
    above = 1
    below = 0
    for i in a.array
        measure!(i.data, w)
        i.offx = width
        width += i.data.sizex
        above = max(above, i.data.basey)
        below = max(below, i.data.sizey - i.data.basey)
    end
    for i in a.array
        i.offy = above - i.data.basey
    end
    a.sizey = above + below
    a.basey = above
    a.sizex = width 
end


function paint!(s::Vector{Vector{Char}}, a::StringBox, offx::Int, offy::Int)
    for i in 1:length(a.str)
        s[1 + offy][offx + i] = a.str[i]
    end
end

function measure!(a::StringBox, w::Int)
    a.sizey = 1
    a.basey = 1
    a.sizex = length(a.str)
end

function paint!(s::Vector{Vector{Char}}, a::SuperscriptBox, offx::Int, offy::Int)
    paint!(s, a.sup.data, offx + a.sup.offx, offy + a.sup.offy)
end

function measure!(a::SuperscriptBox, w::Int)
    measure!(a.sup.data, w)
    a.sup.offx = 0
    a.sup.offy = 0
    a.sizex = a.sup.data.sizex
    a.sizey = a.sup.data.sizey + 1
    a.basey = a.sup.data.sizey + 1
end

function paint!(s::Vector{Vector{Char}}, a::FractionBox, offx::Int, offy::Int)
    paint!(s, a.num.data, offx + a.num.offx, offy + a.num.offy)
    paint!(s, a.den.data, offx + a.den.offx, offy + a.den.offy)
    for i in 3:a.sizex
        s[offy + a.basey][offx + i - 1] = '-'
    end
end

function measure!(a::FractionBox, w::Int)
    measure!(a.num.data, w)
    measure!(a.den.data, w)
    a.sizex = 2 + max(a.num.data.sizex, a.den.data.sizex)
    a.sizey = a.num.data.sizey + 1 + a.den.data.sizey
    a.basey = a.num.data.sizey + 1
    a.num.offx = div(a.sizex - a.num.data.sizex, 2)
    a.num.offy = 0
    a.den.offx = div(a.sizex - a.den.data.sizex, 2)
    a.den.offy = a.basey
end

function paint!(s::Vector{Vector{Char}}, a::BigBox, offx::Int, offy::Int)
    s[1 + offy][1 + offx] = a.char
end

function measure!(a::BigBox, w::Int)
    a.sizex = 1
    a.sizey = 1
    a.basey = 1
end


function paint!(s::Vector{Vector{Char}}, a::GridBox, offx::Int, offy::Int)
    for row in a.array
        for i in row
            paint!(s, i.data, offx + i.offx, offy + i.offy)
        end
    end
    for y in 1:a.sizey
        s[offy + y][1 + offx] = '['
        s[offy + y][offx + a.sizex] = ']'
    end
end

function measure!(a::GridBox, w::Int)
    maxncols = 0
    nrows = length(a.array)
    for row in a.array
        maxncols = max(maxncols, length(row))
        for i in row
            measure!(i.data, w)
        end
    end

    colwidth = zeros(Int, maxncols)
    rowabove = zeros(Int, nrows)
    rowbelow = zeros(Int, nrows)
    for i in 1:length(a.array)
        row = a.array[i]
        for j in 1:length(row)
            rowabove[i] = max(rowabove[i], row[j].data.basey)
            rowbelow[i] = max(rowbelow[i], row[j].data.sizey - row[j].data.basey)
            colwidth[j] = max(colwidth[j], row[j].data.sizex)
        end
    end

    coloffset = zeros(Int, maxncols)
    off = 1
    for j in 1:maxncols
        coloffset[j] = off
        extra = 1
        if j < maxncols
            maxheight = 1
            for i in 1:length(a.array)
                row = a.array[i]
                if j + 1 <= length(row)
                    maxheight = max(maxheight, row[j].data.sizey)
                    maxheight = max(maxheight, row[j + 1].data.sizey)
                end
            end
            if maxheight >= 5
                extra += 2
            elseif maxheight >= 2
                extra += 1
            end 
        end
        off += colwidth[j] +  extra
    end
    a.sizex = max(off, 2)

    off = 0
    for i in 1:length(a.array)
        row = a.array[i]
        for j in 1:length(row)
            row[j].offx = coloffset[j]
            row[j].offy = off + rowabove[i] - row[j].data.basey
        end
        extra = 0
        if i < length(a.array)
            maxheight = max(rowabove[i] + rowbelow[i],
                            rowabove[i + 1] + rowbelow[i + 1])
            if maxheight >= 5
                extra += 2
            elseif maxheight >= 2
                extra += 1
            end 
        end
        off += rowabove[i] + rowbelow[i] + extra
    end
    a.sizey = max(off, 1)
    a.basey = div(a.sizey + 2, 2)
end
