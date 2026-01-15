import AbstractAlgebra.PrettyPrinting

@testset "PrettyPrinting: Expression to string" begin

   function just_string(x)
      return AbstractAlgebra.expr_to_string(x)
   end

   @test just_string(-2) == "-2"
   @test just_string(BigInt(2)) == "2"
   @test just_string(:(+())) == "+()"
   @test just_string(:(-())) == "-()"
   @test just_string(:(+(a))) == "+a"
   @test just_string(:(-(a))) == "-a"
   @test just_string(:(*(a))) == "*(a)"
   @test just_string(:(a+b+c+d)) == "a + b + c + d"
   @test just_string(:(a-b-c-d)) == "a - b - c - d"
   @test just_string(:((a+b)+(c+d))) == "a + b + (c + d)"
   @test just_string(:((a+b)-(c+d))) == "a + b - (c + d)"
   @test just_string(:((a+b)-c+d)) == "a + b - c + d"
   @test just_string(:((a-b)-(c-d))) == "a - b - (c - d)"

   @test just_string(:(-a*b*c*d)) == "-a*b*c*d"
   @test just_string(:(-a/b/c/d)) == "-a/b/c/d"
   @test just_string(:((a*b)/(c*d))) == "a*b/(c*d)"
   @test just_string(:(a*(b/(c*d)))) == "a*(b/(c*d))"
   @test just_string(:(a*(b/c)*d)) == "a*(b/c)*d"
   @test just_string(:(a*((b/c)*d))) == "a*(b/c*d)"

   @test just_string(:(a+b(c,d))) == "a + b(c, d)"

   @test just_string(:(-a/b)) == "-a/b"
   @test just_string(:(-a//b)) == "-a//b"
   @test just_string(:(-(a/b))) == "-(a/b)"
   @test just_string(:(-(a//b))) == "-(a//b)"
   @test just_string(:(a^((-b)/c))) == "a^(-b/c)"
   @test just_string(:(a^(-(b/c)))) == "a^-(b/c)"
   @test just_string(:(a^-(b*c^d))) == "a^-(b*c^d)"
   @test just_string(:((a^-b)*(c^d))) == "a^-b*c^d"

   # // has a higher precedence than /
   @test just_string(:(a*b//c)) == "a*b//c"
   @test just_string(:(a*(b//c))) == "a*b//c"
   @test just_string(:((a*b)//c)) == "(a*b)//c"
   @test just_string(:(-a*b//c*d^e)) == "-a*b//c*d^e"
   @test just_string(:((-a)*(b//c)*d^e)) == "-a*b//c*d^e"

   @test just_string(:(a(b, c))) == "a(b, c)"
   @test just_string(:(a[b, c])) == "a[b, c]"
   @test just_string(:([a, b, c])) == "[a, b, c]"
   @test just_string(:((a, b, c))) == "(a, b, c)"
   @test just_string(Expr(:list, :a, :b, :c)) == "{a, b, c}"
   @test just_string(Expr(:series, :a, :b, :c)) == "a, b, c"
   @test just_string(Expr(:row, :a, :b, :c)) == "a b c"
   @test just_string(Expr(:sequence, :a, :b, :c)) == "abc"
   @test just_string(:([a b; c d])) == "[a b; c d]"

   @test just_string(:(if a; b; end;)) isa String
   @test just_string(1.2) isa String

   function canonical_string(x)
      return AbstractAlgebra.expr_to_string(AbstractAlgebra.canonicalize(x))
   end

   @test canonical_string(:(+())) == "0"
   @test canonical_string(:(-())) == "0"
   @test canonical_string(:(*())) == "1"
   @test canonical_string(:(+(a))) == "a"
   @test canonical_string(:(-(a))) == "-a"
   @test canonical_string(:(*(a))) == "a"
   @test canonical_string(:(a*1*1)) == "a"
   @test canonical_string(:(1*1*1)) == "1"
   @test canonical_string(:(a+b*c)) == "a + b*c"
   @test canonical_string(:(a-b*c/d)) == "a - b*c/d"
   @test canonical_string(:(a-b*c/d*0)) == "a"
   @test canonical_string(:(a+(-1)*b*c)) == "a - b*c"
   @test canonical_string(:(a*((-1)*b*c))) == "-a*b*c"
   @test canonical_string(:(a*((b/c)*d))) == "a*(b/c)*d"
   @test canonical_string(:(a+-(b+d)*c)) == "a - (b + d)*c"
   @test canonical_string(:(a+(-b)*c)) == "a - b*c"
   @test canonical_string(:(a+(-b)*(-c))) == "a + b*c"
   @test canonical_string(:(a+((-1)*b*c))) == "a - b*c"
   @test canonical_string(:((-a)/b)) == "-a/b"
   @test canonical_string(:(-(a/b))) == "-a/b"
   @test canonical_string(:(-(a/-b))) == "-a/-b"
   @test canonical_string(:((-a)^b)) == "(-a)^b"
   @test canonical_string(:(sqrt(-a)*-b)) == "-sqrt(-a)*b"
   @test canonical_string(:(a^(-b/c) - ((-a)*b)*c - ((a*b)*(-c*d)))) ==
                                                "a^(-b/c) + a*b*c + a*b*c*d"
   @test canonical_string(:((-a/b)*c+d*(e/f)*g^h + (-i)*j/k^l)) ==
                                              "-a/b*c + d*(e/f)*g^h - i*j/k^l"

   @test canonical_string(:(x^((-1)//2))) == "x^(-1//2)"

   @test canonical_string(:(cdot(cdot(-a,b),cdot(c,d)))) == "(-a * b) * (c * d)"

   @test canonical_string(:([a b; c d])) == "[a b; c d]"
   @test canonical_string(:([-a-b -c+d; -a -b-c+d])) == "[-a-b -c+d; -a -b-c+d]"

   @test canonical_string(:(a + (b == c))) == "a + (b = c)"
   @test canonical_string(:(a + (b == c == d))) == "a + (b = c = d)"
   @test canonical_string(:(a = b = c+d*1)) == "a = b = c + d"
   @test canonical_string(:((a = b) = c)) == "(a = b) = c"
   @test canonical_string(:(a == b == c)) == "a = b = c"
   @test canonical_string(:((a == b) == c)) == "(a = b) = c"
   @test canonical_string(:(a == (b == c) == d)) == "a = (b = c) = d"
   @test canonical_string(:(a == b)) == "a = b"
   @test canonical_string(:(a != b)) == "a != b"
   @test canonical_string(:(a > b))  == "a > b"
   @test canonical_string(:(a >= b)) == "a >= b"
   @test canonical_string(:(a < b))  == "a < b"
   @test canonical_string(:(a <= b)) == "a <= b"

   @test canonical_string(:(if a; b; end;)) isa String
   @test canonical_string(1.2) isa String

   # bugfix: ensure multiple levels of canonicalize produce :(a + b + -c)
   @test canonical_string(:((*)(a) + (*)((*)(b) + -1*c))) == "a + b - c"
   @test canonical_string(:(f((*)(a) + (*)((*)(b) + -1*c)))) == "f(a + b - c)"

   function latex_string(x)
      return AbstractAlgebra.expr_to_latex_string(AbstractAlgebra.canonicalize(x))
   end

   @test latex_string(Symbol("a")) == "a"
   @test latex_string(Symbol("a_1")) == "a_{1}"
   @test latex_string(Symbol("α")) == "{\\alpha}"
   @test latex_string(Symbol("α_1")) == "{\\alpha}_{1}"
   @test latex_string(Symbol("x1")) == "\\mathop{\\mathrm{x1}}"
   @test latex_string(Symbol("xy_1")) == "\\mathop{\\mathrm{xy}}_{1}"
   @test latex_string(Symbol("sin")) == "\\mathop{\\mathrm{sin}}"
   @test latex_string(Symbol("sin_cos")) == "\\mathop{\\mathrm{sin\\_cos}}"
   @test latex_string(Symbol("sin_1")) == "\\mathop{\\mathrm{sin}}_{1}"
   @test latex_string(Symbol("sin_cos_1")) == "\\mathop{\\mathrm{sin\\_cos}}_{1}"
   @test latex_string(Symbol("αaβb_1_2")) == "\\mathop{\\mathrm{{\\alpha}a{\\beta}b}}_{1,2}"

   @test latex_string(:(+())) == "0"
   @test latex_string(:(-())) == "0"
   @test latex_string(:(*())) == "1"
   @test latex_string(:(+(a))) == "a"
   @test latex_string(:(-(a))) == "-a"
   @test latex_string(:(*(a))) == "a"
   @test latex_string(:(a*1*1)) == "a"
   @test latex_string(:(1*1*1)) == "1"
   @test latex_string(:(a+b*c)) == "a + b c"
   @test latex_string(:(a-(b*c)//d)) == "a - \\frac{b c}{d}"
   @test latex_string(:(a-b*c/d*0)) == "a"
   @test latex_string(:(a+(-1)*b*c)) == "a - b c"
   @test latex_string(:(a*((-1)*b*c))) == "-a b c"
   @test latex_string(:(a*((b//c)*d))) == "a \\frac{b}{c} d"
   @test latex_string(:(a+-(b+d)*c)) == "a - \\left(b + d\\right) c"
   @test latex_string(:(a+(-b)*c)) == "a - b c"
   @test latex_string(:(a+(-b)*(-c))) == "a + b c"
   @test latex_string(:(a+((-1)*b*c))) == "a - b c"
   @test latex_string(:((-a)//b)) == "-\\frac{a}{b}"
   @test latex_string(:(-(a//b))) == "-\\frac{a}{b}"
   @test latex_string(:(-(a//-b))) == "-\\frac{a}{-b}"
   @test latex_string(:((-a)^b)) == "\\left(-a\\right)^{b}"
   @test latex_string(:(a+sqrt(b*c))) == "a + \\sqrt{b c}"
   @test latex_string(:(sqrt(a)^b)) == "\\left(\\sqrt{a}\\right)^{b}"
   @test latex_string(:((a//b)^c)) == "\\left(\\frac{a}{b}\\right)^{c}"
   @test latex_string(:(ff(-a)*-b)) == "-\\mathop{\\mathrm{ff}}\\left(-a\\right) b"
   @test latex_string(:(a^(-b//c) - ((-a)*b)*c - ((a*b)*(-c*d)))) ==
                                         "a^{-\\frac{b}{c}} + a b c + a b c d"
   @test latex_string(:((-a//b)*c+d*(e//f)*g^h + ((-i)*j)//k^l)) ==
                  "-\\frac{a}{b} c + d \\frac{e}{f} g^{h} - \\frac{i j}{k^{l}}"

   @test latex_string(:(x^((-1)//2))) == "x^{-\\frac{1}{2}}"

   # / has same precedence as * in julia
   @test latex_string(:(x^(-a*b/c))) == "x^{-a b/c}"
   # // has higher precedence in julia
   @test latex_string(:(x^(-a*b//c))) == "x^{-a \\frac{b}{c}}"

   @test latex_string(:((-a*b/c)^x)) == "\\left(-a b/c\\right)^{x}"
   @test latex_string(:((a*b/c)^x)) == "\\left(a b/c\\right)^{x}"
   @test latex_string(:((-a*b//c)^x)) == "\\left(-a \\frac{b}{c}\\right)^{x}"
   @test latex_string(:((a*b//c)^x)) == "\\left(a \\frac{b}{c}\\right)^{x}"

   @test latex_string(:(cdot(cdot(-a,b),cdot(c,d)))) ==
             "\\left(-a \\cdot b\\right) \\cdot \\left(c \\cdot d\\right)"

   @test latex_string(:(2*α^2-1*α+1)) == "2 {\\alpha}^{2} - {\\alpha} + 1"

   @test latex_string(Expr(:matrix, :([a b; c d]))) ==
           "\\left(\\begin{array}{cc}\na & b \\\\\nc & d\n\\end{array}\\right)"

   @test latex_string(:(a = b = c+d*1)) == "a = b = c + d"
   @test latex_string(:((a = b) = c)) == "\\left(a = b\\right) = c"
   @test latex_string(:(a == b == c)) == "a = b = c"
   @test latex_string(:((a == b) == c)) == "\\left(a = b\\right) = c"
   @test latex_string(:(a == (b == c) == d)) == "a = \\left(b = c\\right) = d"
   @test latex_string(:(a <= b + -2*c <= d != e)) == "a \\le b - 2 c \\le d \\neq e"
   @test latex_string(:(a == b)) == "a = b"
   @test latex_string(:(a != b)) == "a \\neq b"
   @test latex_string(:(a > b))  == "a > b"
   @test latex_string(:(a >= b)) == "a \\ge b"
   @test latex_string(:(a < b))  == "a < b"
   @test latex_string(:(a <= b)) == "a \\le b"

   @test latex_string(:(if a; b; end;)) isa String
   @test latex_string(1.2) isa String

   function limited_latex_string(x)
      x = AbstractAlgebra.canonicalize(x)
      b = IOBuffer()
      c = IOContext(b, :size_limit => 20)
      AbstractAlgebra.show_obj(c, MIME("text/latex"), x)
      return String(take!(b))
   end

   function limited_string(x)
      x = AbstractAlgebra.canonicalize(x)
      b = IOBuffer()
      c = IOContext(b, :size_limit => 20)
      AbstractAlgebra.show_obj(c, MIME("text/plain"), x)
      return String(take!(b))
   end

   e = Expr(:call, :+, [:(x^$i) for i in 1:200]...)
   @test limited_latex_string(e) == "x^{1} + x^{2} + x^{3} + x^{4} +"*
                                    " {\\ldots} + x^{198} + x^{199} + x^{200}"
   @test limited_string(e) == "x^1 + x^2 + x^3 + x^4 +"*
                              " ... + x^198 + x^199 + x^200"

   e = Expr(:call, :*, [:(x^$i) for i in 1:200]...)
   @test limited_latex_string(e) == "x^{1} x^{2} x^{3} x^{4} "*
                                    "{\\ldots} x^{198} x^{199} x^{200}"
   @test limited_string(e) == "x^1*x^2*x^3*x^4*...*x^198*x^199*x^200"

   e = Expr(:call, [:(x^$i) for i in 1:200]...)
   @test limited_latex_string(e) == "\\left(x^{1}\\right)\\left(x^{2}, x^{3},"*
                         " x^{4}, {\\ldots}, x^{198}, x^{199}, x^{200}\\right)"
   @test limited_string(e) == "(x^1)(x^2, x^3, x^4, ..., x^198, x^199, x^200)"

   e = Expr(:ref, [:(x^$i) for i in 1:200]...)
   @test limited_latex_string(e) == "\\left(x^{1}\\right)\\left[x^{2}, "*
                   "x^{3}, x^{4}, {\\ldots}, x^{198}, x^{199}, x^{200}\\right]"
   @test limited_string(e) == "(x^1)[x^2, x^3, x^4, ..., x^198, x^199, x^200]"

   e = Expr(:list, [:(x^$i) for i in 1:200]...)
   @test limited_latex_string(e) == "\\left{x^{1}, x^{2}, x^{3}, x^{4}, "*
                                 "{\\ldots}, x^{198}, x^{199}, x^{200}\\right}"
   @test limited_string(e) == "{x^1, x^2, x^3, x^4, ..., x^198, x^199, x^200}"

   e = Expr(:sequence, Expr(:text, "Polynomial ring in "),
                       :x,
                       Expr(:text, " over "),
                       Expr(:latex_form, "ZZ", "\\mathbb{Z}"))
   e = Expr(:call, :/, e, :(x^2+1))
   @test canonical_string(e) == "(Polynomial ring in x over ZZ)/(x^2 + 1)"
   @test latex_string(e) == "\\left(\\text{Polynomial ring in }x\\"*
                     "text{ over }\\mathbb{Z}\\right)/\\left(x^{2} + 1\\right)"

   e = :([a b;c;d e])
   @test canonical_string(e) == "[a b; c; d e]"
   @test latex_string(e) ==
                  "\\begin{array}{cc}\na & b \\\\\nc \\\\\nd & e\n\\end{array}"

   e = Expr(:matrix, e)
   @test canonical_string(e) == "[a b; c; d e]"
   @test latex_string(e) ==
   "\\left(\\begin{array}{cc}\na & b \\\\\nc \\\\\nd & e\n\\end{array}\\right)"


   R, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"])
   p = a + b^2 + c^3
   @test sprint(show, p) == "a + b^2 + c^3"
   @test sprint(show, p, context = :compact => true) == "a + b^2 + c^3"
   @test sprint(show, p, context = :terse_level => 1) == "a+b^2+c^3"

   Qx, x = QQ["x"]
   AbstractAlgebra.set_html_as_latex(true)
   @test AbstractAlgebra.get_html_as_latex() == true
   @test sprint(show, "text/html", x^12) == "\$x^{12}\$"
   fl = AbstractAlgebra.set_html_as_latex(false)
   @assert fl == true
   @test AbstractAlgebra.get_html_as_latex() == false
   @test sprint(show, "text/html", x^12) == "x^12"

   R,(x,y) = ZZ["x","y"]

   @test latex_string(Expr(:vcat)) == "\\text{empty}"

   m = matrix(R, [x^i*y^j for i in 0:20, j in 0:20])
   e = AbstractAlgebra.expressify(m)
   @test length(latex_string(e)) > 21*21
   @test limited_latex_string(e) == "\\begin{array}{ccc}\n"*
                                    "1 & \\cdots  & y^{20} \\\\\n"*
                                    "\\vdots  & \\ddots  & \\vdots  \\\\\n"*
                                    "x^{20} & \\cdots  & x^{20} y^{20}\n"*
                                    "\\end{array}"

   m = matrix(R, [x^i*y^j for i in 0:20, j in 0:1])
   e = AbstractAlgebra.expressify(m)
   @test length(latex_string(e)) > 21*2
   @test limited_latex_string(e) == "\\begin{array}{cc}\n"*
                                    "1 & y \\\\\n"*
                                    "x & x y \\\\\n"*
                                    "x^{2} & x^{2} y \\\\\n"*
                                    "x^{3} & x^{3} y \\\\\n"*
                                    "\\vdots  & \\vdots  \\\\\n"*
                                    "x^{17} & x^{17} y \\\\\n"*
                                    "x^{18} & x^{18} y \\\\\n"*
                                    "x^{19} & x^{19} y \\\\\n"*
                                    "x^{20} & x^{20} y\n"*
                                    "\\end{array}"

   R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
   @test AbstractAlgebra.obj_to_string_wrt_times(x^2) == "x^2"
   @test AbstractAlgebra.obj_to_string_wrt_times(x*y) == "(x*y)"
   @test AbstractAlgebra.obj_to_string_wrt_times(x + y) == "(x + y)"
end

@testset "PrettyPrinting: Special printing macros" begin
  # TODO

end


@testset "PrettyPrinting: Unicode preferences" begin
  # TODO

end

@testset "PrettyPrinting: Three print modes" begin
  # Test various examples from the Oscar manual

  #
  #
  #
  struct NewRing
    base_ring
  end

  base_ring(R::NewRing) = R.base_ring

  function Base.show(io::IO, ::MIME"text/plain", R::NewRing)
    println(io, "I am a new ring")  # at least one new line is needed
    println(io, "I print with newlines")
    print(io, base_ring(R)) # the last print statement must not add a new line
  end

  function Base.show(io::IO, R::NewRing)
    if PrettyPrinting.is_terse(io)
      # no nested printing
      print(io, "terse printing of newring ")
    else
      # nested printing allowed, preferably terse
      print(io, "one line printing of newring with ")
      print(PrettyPrinting.terse(io), "terse ", base_ring(R))
    end
  end

  R = NewRing(QQ)
  @test PrettyPrinting.repr_detailed(R) ==
        """I am a new ring
           I print with newlines
           Rationals"""
  @test PrettyPrinting.repr_oneline(R) == "one line printing of newring with terse Rationals"
  @test PrettyPrinting.repr_terse(R) == "terse printing of newring "

  #
  #
  #
  struct A{T}
    x::T
  end

  function Base.show(io::IO, a::A)
    io = AbstractAlgebra.pretty(io)
    println(io, "Something of type A")
    print(io, AbstractAlgebra.Indent(), "over ", AbstractAlgebra.Lowercase(), a.x)
    print(io, AbstractAlgebra.Dedent()) # don't forget to undo the indentation!
  end

  struct B
  end

  function Base.show(io::IO, b::B)
    io = AbstractAlgebra.pretty(io)
    print(io, AbstractAlgebra.LowercaseOff(), "Hilbert thing")
  end

  x = A(2)
  y = """
      Something of type A
        over 2"""
  @test PrettyPrinting.repr_detailed(x) == y
  @test PrettyPrinting.repr_oneline(x) == y
  @test PrettyPrinting.repr_terse(x) == y

  x = A(A(2))
  y = """
      Something of type A
        over something of type A
          over 2"""
  @test PrettyPrinting.repr_detailed(x) == y
  @test PrettyPrinting.repr_oneline(x) == y
  @test PrettyPrinting.repr_terse(x) == y

  x = A(B())
  y = """
      Something of type A
        over Hilbert thing"""
  @test PrettyPrinting.repr_detailed(x) == y
  @test PrettyPrinting.repr_oneline(x) == y
  @test PrettyPrinting.repr_terse(x) == y

end

@testset "PrettyPrinting: Intendation and Decapitalization" begin
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  @test io isa AbstractAlgebra.PrettyPrinting.IOCustom
  @test io === AbstractAlgebra.pretty(io)
  print(io, AbstractAlgebra.Indent(), "test")
  println(io, " test")
  print(io, "test")
  println(io, AbstractAlgebra.Indent())
  print(io, "test")
  println(io, AbstractAlgebra.Dedent())
  print(io, "test")
  @test String(take!(io)) == "  test test\n" *
                             "  test\n" *
                             "    test\n" *
                             "  test"

  # Test unicode
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  println(io, "testing unicode")
  print(io, AbstractAlgebra.Indent(), "ŎŚĊĂŖ")
  @test String(take!(io)) == "testing unicode\n" *
                             "  ŎŚĊĂŖ"

  # Test evil unicodes
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  _, c = displaysize(io)
  print(io, AbstractAlgebra.Indent())
  ellipses = String([0xe2, 0x80, 0xa6])
  wedge = String([0xe2, 0x88, 0xa7])
  iacute = String([0xc3, 0xad])
  str = wedge ^25 * ellipses^25 * iacute^50
  print(io, "aa", str)
  @test String(take!(io)) == "  aa∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧" *
                             "…………………………………………………………………" *
                             "íííííííííííííííííííííííííí\n" *
                             "  íííííííííííííííííííííííí"


  # Test string longer than width
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  _, c = displaysize(io)
  print(io, AbstractAlgebra.Indent())
  println(io, "t"^c)
  println(io, "aa", "t"^c)
  print(io, AbstractAlgebra.Indent())
  print(io, "aa", "t"^c)
  @test String(take!(io)) == "  " * "t"^(c - 2) * "\n" *
                             "  tt\n" *
                             "  aa" * "t"^(c - 4) * "\n" *
                             "  tttt\n" *
                             "    aa" * "t"^(c - 6) * "\n" *
                             "    tttttt"

  # Test unicode string longer than width
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  _, c = displaysize(io)
  print(io, AbstractAlgebra.Indent())
  println(io, "Ŏ"^c)
  println(io, "aa", "Ś"^c)
  print(io, AbstractAlgebra.Indent())
  print(io, "aa", "Ŗ"^c)
  @test String(take!(io)) == "  " * "Ŏ"^(c-2) * "\n" *
                             "  ŎŎ" * "\n" *
                             "  aa" * "Ś"^(c-4) * "\n" *
                             "  ŚŚŚŚ" * "\n" *
                             "    aa" * "Ŗ"^(c-6) * "\n" *
                             "    ŖŖŖŖŖŖ"

  # Test evil unicode string much longer than width
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  _, c = displaysize(io)
  ellipses = String([0xe2, 0x80, 0xa6])
  wedge = String([0xe2, 0x88, 0xa7])
  iacute = String([0xc3, 0xad])
  evil_a = String([0x61, 0xcc, 0x81, 0xcc, 0xa7, 0xcc, 0xa7])
  print(io, AbstractAlgebra.Indent())
  println(io, "Ŏ"^c)
  println(io, ellipses^c)
  println(io, "aa", "Ś"^c)
  println(io, "bb", wedge^c)
  print(io, AbstractAlgebra.Indent())
  println(io, "aa", "Ŗ"^c)
  print(io, iacute^c)
  println(io, evil_a^c)
  print(io, evil_a^c)
  @test String(take!(io)) == "  " * "Ŏ"^(c-2) * "\n" *
                             "  ŎŎ" * "\n" *
                             "  " * ellipses^(c-2) * "\n" *
                             "  " * ellipses^2 * "\n" *
                             "  aa" * "Ś"^(c-4) * "\n" *
                             "  ŚŚŚŚ" * "\n" *
                             "  bb" * wedge^(c-4) * "\n" *
                             "  " * wedge^4 * "\n" *
                             "    aa" * "Ŗ"^(c-6) * "\n" *
                             "    ŖŖŖŖŖŖ" * "\n" *
                             "    " * iacute^(c-4) * "\n" *
                             "    " * iacute^4 * evil_a^(c-8) * "\n" *
                             "    " * evil_a^(8) * "\n" *
                             "    " * evil_a^(c-4) * "\n" *
                             "    " * evil_a^4

  # Test graphemes with non standard width
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  _, c = displaysize(io)
  boat = String([0xe2, 0x9b, 0xb5])
  family = String([0xf0, 0x9f, 0x91, 0xaa])
  print(io, AbstractAlgebra.Indent())
  println(io, (boat * family)^40)
  print(io, (boat * family)^40)
  @test String(take!(io)) == "  " * (boat*family)^19 * boat * "\n" *
                             "  " * (family*boat)^19 * family * "\n" *
                             "  " * boat * family * "\n" *
                             "  " * (boat*family)^19 * boat * "\n" *
                             "  " * (family*boat)^19 * family * "\n" *
                             "  " * boat * family

  # Test graphemes with standard and non standard width mixed in
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  _, c = displaysize(io)
  ellipses = String([0xe2, 0x80, 0xa6])
  wedge = String([0xe2, 0x88, 0xa7])
  iacute = String([0xc3, 0xad])
  evil_a = String([0x61, 0xcc, 0x81, 0xcc, 0xa7, 0xcc, 0xa7])
  boat = String([0xe2, 0x9b, 0xb5])
  family = String([0xf0, 0x9f, 0x91, 0xaa])
  print(io, AbstractAlgebra.Indent())
  println(io, "Ŏ"^c)
  println(io, ellipses^c)
  println(io, "aa", "Ś"^c)
  println(io, boat^(3*c))
  println(io, "bb", wedge^c)
  print(io, AbstractAlgebra.Indent())
  println(io, "aa", "Ŗ"^c)
  println(io, family^(3*c))
  println(io, iacute^c)
  println(io, evil_a^c)
  print(io, evil_a^c)
  @test String(take!(io)) == "  " * "Ŏ"^(c-2) * "\n" *
                             "  ŎŎ" * "\n" *
                             "  " * ellipses^(c-2) * "\n" *
                             "  " * ellipses^2 * "\n" *
                             "  aa" * "Ś"^(c-4) * "\n" *
                             "  ŚŚŚŚ" * "\n" *
                             "  " * boat^39 * "\n" *
                             "  " * boat^39 * "\n" *
                             "  " * boat^39 * "\n" *
                             "  " * boat^39 * "\n" *
                             "  " * boat^39 * "\n" *
                             "  " * boat^39 * "\n" *
                             "  " * boat^6 * "\n" *
                             "  bb" * wedge^(c-4) * "\n" *
                             "  " * wedge^4 * "\n" *
                             "    aa" * "Ŗ"^(c-6) * "\n" *
                             "    ŖŖŖŖŖŖ" * "\n" *
                             "    " * family^38 * "\n" *
                             "    " * family^38 * "\n" *
                             "    " * family^38 * "\n" *
                             "    " * family^38 * "\n" *
                             "    " * family^38 * "\n" *
                             "    " * family^38 * "\n" *
                             "    " * family^12 * "\n" *
                             "    " * iacute^(c-4) * "\n" *
                             "    " * iacute^4 *"\n" *
                             "    " * evil_a^(c-4) * "\n" *
                             "    " * evil_a^(4) * "\n" *
                             "    " * evil_a^(c-4) * "\n" *
                             "    " * evil_a^4

  # Test too much indentation
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  _, c = displaysize(io)
  for i in 1:c
    print(io, AbstractAlgebra.Indent())
  end
  print(io, "test")
  @test String(take!(io)) isa String

  # Lowercase/LowercaseOff
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  print(io, AbstractAlgebra.Lowercase(), "Test")
  @test String(take!(io)) == "test"
  print(io, AbstractAlgebra.Lowercase(), AbstractAlgebra.LowercaseOff(), "Test")
  @test String(take!(io)) == "Test"

  # Fix bug for pretty(IOCustom(pretty(io))) forgetting everything
  io = IOBuffer()
  io = AbstractAlgebra.pretty(io, force_newlines = true)
  print(io, AbstractAlgebra.Lowercase())
  io = AbstractAlgebra.pretty(IOContext(io), force_newlines = true)
  print(io, "A")
  @test String(take!(io.io)) == "a"

  # Check that printing to IOBuffer does not introduce newlines
  io = AbstractAlgebra.pretty(IOBuffer())
  _, c = displaysize(io)
  print(io, AbstractAlgebra.Indent(), "a"^(c + 1))
  @test String(take!(io.io)) == "  " * "a"^(c + 1)

  # Fix #1435
  io = IOContext(AbstractAlgebra.pretty(IOBuffer()))
  print(io, AbstractAlgebra.Indent(), sin)
  @test String(take!(io.io)) == "  sin"
  f = x-> 2*x
  io = IOContext(AbstractAlgebra.pretty(IOBuffer()))
  print(io, AbstractAlgebra.Indent(), f)
  @test startswith(String(take!(io.io)), "  ")
end
