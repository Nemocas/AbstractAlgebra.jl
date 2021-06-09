@testset "PrettyPrinting" begin

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

   @test canonical_string(:(if a; b; end;)) isa String
   @test canonical_string(1.2) isa String

   function latex_string(x)
      return AbstractAlgebra.expr_to_latex_string(AbstractAlgebra.canonicalize(x))
   end

   @test latex_string(Symbol("a")) == "a"
   @test latex_string(Symbol("α")) == "{\\alpha}"
   @test latex_string(Symbol("x1")) == "\\operatorname{x1}"
   @test latex_string(Symbol("xy_1")) == "\\operatorname{xy}_{1}"
   @test latex_string(Symbol("sin")) == "\\operatorname{sin}"
   @test latex_string(Symbol("sin_cos")) == "\\operatorname{sin\\_cos}"
   @test latex_string(Symbol("sin_1")) == "\\operatorname{sin}_{1}"
   @test latex_string(Symbol("sin_cos_1")) == "\\operatorname{sin\\_cos}_{1}"
   @test latex_string(Symbol("αaβb_1_2")) == "\\operatorname{{\\alpha}a{\\beta}b}_{1,2}"

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
   @test latex_string(:(ff(-a)*-b)) == "-\\operatorname{ff}\\left(-a\\right) b"
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


   R, (a, b, c) = PolynomialRing(QQ, ["a", "b", "c"])
   p = a + b^2 + c^3
   @test sprint(show, p) == "a + b^2 + c^3"
   @test sprint(show, p, context = :compact => true) == "a + b^2 + c^3"
   @test sprint(show, p, context = :terse => true) == "a+b^2+c^3"

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
end

