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
   @test canonical_string(:(if a; b; end;)) isa String
   @test canonical_string(1.2) isa String

   function latex_string(x)
      return AbstractAlgebra.expr_to_latex_string(AbstractAlgebra.canonicalize(x))
   end

   @test latex_string(:(+())) == "0"
   @test latex_string(:(-())) == "0"
   @test latex_string(:(*())) == "1"
   @test latex_string(:(+(a))) == "a"
   @test latex_string(:(-(a))) == "-a"
   @test latex_string(:(*(a))) == "a"
   @test latex_string(:(a*1*1)) == "a"
   @test latex_string(:(1*1*1)) == "1"
   @test latex_string(:(a+b*c)) == "a + b c"
   @test latex_string(:(a-b*c/d)) == "a - \\frac{b c}{d}"
   @test latex_string(:(a-b*c/d*0)) == "a"
   @test latex_string(:(a+(-1)*b*c)) == "a - b c"
   @test latex_string(:(a*((-1)*b*c))) == "-a b c"
   @test latex_string(:(a*((b/c)*d))) == "a \\frac{b}{c} d"
   @test latex_string(:(a+-(b+d)*c)) == "a - \\left(b + d\\right) c"
   @test latex_string(:(a+(-b)*c)) == "a - b c"
   @test latex_string(:(a+(-b)*(-c))) == "a + b c"
   @test latex_string(:(a+((-1)*b*c))) == "a - b c"
   @test latex_string(:((-a)/b)) == "-\\frac{a}{b}"
   @test latex_string(:(-(a/b))) == "-\\frac{a}{b}"
   @test latex_string(:(-(a/-b))) == "-\\frac{a}{-b}"
   @test latex_string(:((-a)^b)) == "\\left(-a\\right)^{b}"
   @test latex_string(:(a+sqrt(b*c))) == "a + \\sqrt{b c}"
   @test latex_string(:(ff(-a)*-b)) == "-f f\\left(-a\\right) b"
   @test latex_string(:(a^(-b/c) - ((-a)*b)*c - ((a*b)*(-c*d)))) ==
                                         "a^{-\\frac{b}{c}} + a b c + a b c d"                             
   @test latex_string(:((-a/b)*c+d*(e/f)*g^h + (-i)*j/k^l)) ==
                  "-\\frac{a}{b} c + d \\frac{e}{f} g^{h} - \\frac{i j}{k^{l}}"

   @test latex_string(:(x^((-1)//2))) == "x^{-\\frac{1}{2}}"

   # / has same precedence as * in julia
   @test latex_string(:(x^(-a*b/c))) == "x^{-\\frac{a b}{c}}"
   # // has higher precedence in julia
   @test latex_string(:(x^(-a*b//c))) == "x^{-a \\frac{b}{c}}"

   @test latex_string(:((-a*b/c)^x)) == "\\left(-\\frac{a b}{c}\\right)^{x}"
   @test latex_string(:((a*b/c)^x)) == "\\left(\\frac{a b}{c}\\right)^{x}"

   @test latex_string(:(cdot(cdot(-a,b),cdot(c,d)))) ==
             "\\left(-a \\cdot b\\right) \\cdot \\left(c \\cdot d\\right)"

   @test latex_string(:(2*α^2-1*α+1)) == "2 \\alpha^{2} - \\alpha + 1"

   @test latex_string(:([a b; c d])) ==
                        "\\begin{array}{cc}\na & b \\\\\nc & d\n\\end{array}\n"

   @test latex_string(:(if a; b; end;)) isa String
   @test latex_string(1.2) isa String

end

