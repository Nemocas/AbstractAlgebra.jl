var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Nemo",
    "title": "Nemo",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Nemo-1",
    "page": "Nemo",
    "title": "Nemo",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Introduction-1",
    "page": "Nemo",
    "title": "Introduction",
    "category": "section",
    "text": "Nemo is a computer algebra package for the Julia programming language, maintained by William Hart,  Tommy Hofmann, Claus Fieker, Fredrik Johansson with additional code by Oleksandr Motsak and other contributors.http://nemocas.org (Website)\nhttps://github.com/Nemocas/Nemo.jl (Source code)\nhttp://nemocas.github.io/Nemo.jl/latest/ (Online documentation)The features of Nemo so far include:Multiprecision integers and rationals\nIntegers modulo n\np-adic numbers\nFinite fields (prime and non-prime order)\nNumber field arithmetic\nMaximal orders of number fields\nArithmetic of ideals in maximal orders\nArbitrary precision real and complex balls\nUnivariate polynomials and matrices over the above\nGeneric polynomials, power series, fraction fields, residue rings and matrices"
},

{
    "location": "index.html#Installation-1",
    "page": "Nemo",
    "title": "Installation",
    "category": "section",
    "text": "To use Nemo we require Julia 0.4 or higher. Please see http://julialang.org/downloads for instructions on how to obtain julia for your system.At the Julia prompt simply typejulia> Pkg.add(\"Nemo\")\njulia> Pkg.build(\"Nemo\")Alternatively, if you don't want to set Julia up yourself, Julia and Nemo are available on https://cloud.sagemath.com/."
},

{
    "location": "index.html#Quick-start-1",
    "page": "Nemo",
    "title": "Quick start",
    "category": "section",
    "text": "Here are some examples of using Nemo.This example computes recursive univariate polynomials.julia> using Nemo\n\njulia> R, x = PolynomialRing(ZZ, \"x\")\n(Univariate Polynomial Ring in x over Integer Ring,x)\n\njulia> S, y = PolynomialRing(R, \"y\")\n(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integer Ring,y)\n\njulia> T, z = PolynomialRing(S, \"z\")\n(Univariate Polynomial Ring in z over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integer Ring,z)\n\njulia> f = x + y + z + 1\nz+(y+(x+1))\n\njulia> p = f^30; # semicolon supresses output\n\njulia> @time q = p*(p+1);\n  0.325521 seconds (140.64 k allocations: 3.313 MB)Here is an example using generic recursive ring constructions.julia> using Nemo\n\njulia> R, x = FiniteField(7, 11, \"x\")\n(Finite field of degree 11 over F_7,x)\n\njulia> S, y = PolynomialRing(R, \"y\")\n(Univariate Polynomial Ring in y over Finite field of degree 11 over F_7,y)\n\njulia> T = ResidueRing(S, y^3 + 3x*y + 1)\nResidue ring of Univariate Polynomial Ring in y over Finite field of degree 11 over F_7 modulo y^3+(3*x)*y+(1)\n\njulia> U, z = PolynomialRing(T, \"z\")\n(Univariate Polynomial Ring in z over Residue ring of Univariate Polynomial Ring in y over Finite field of degree 11 over F_7 modulo y^3+(3*x)*y+(1),z)\n\njulia> f = (3y^2 + y + x)*z^2 + ((x + 2)*y^2 + x + 1)*z + 4x*y + 3;\n\njulia> g = (7y^2 - y + 2x + 7)*z^2 + (3y^2 + 4x + 1)*z + (2x + 1)*y + 1;\n\njulia> s = f^12;\n\njulia> t = (s + g)^12;\n\njulia> @time resultant(s, t)\n  0.426612 seconds (705.88 k allocations: 52.346 MB, 2.79% gc time)\n(x^10+4*x^8+6*x^7+3*x^6+4*x^5+x^4+6*x^3+5*x^2+x)*y^2+(5*x^10+x^8+4*x^7+3*x^5+5*x^4+3*x^3+x^2+x+6)*y+(2*x^10+6*x^9+5*x^8+5*x^7+x^6+6*x^5+5*x^4+4*x^3+x+3)Here is an example using matrices.julia> using Nemo\n\njulia> M = MatrixSpace(R, 40, 40)();\n\njulia> R, x = PolynomialRing(ZZ, \"x\")\n(Univariate Polynomial Ring in x over Integer Ring,x)\n\njulia> M = MatrixSpace(R, 40, 40)();\n\njulia> for i in 1:40\n          for j in 1:40\n             M[i, j] = R(map(fmpz, rand(-20:20, 3)))\n          end\n       end\n\njulia> @time det(M);\n  0.174888 seconds (268.40 k allocations: 26.537 MB, 4.47% gc time)And here is an example with power series.julia> using Nemo\n\njulia> R, x = QQ[\"x\"]\n(Univariate Polynomial Ring in x over Rational Field,x)\n\njulia> S, t = PowerSeriesRing(R, 100, \"t\")\n(Univariate power series ring in t over Univariate Polynomial Ring in x over Rational Field,t+O(t^101))\n\njulia> u = t + O(t^100)\nt+O(t^100)\n\njulia> @time divexact((u*exp(x*u)), (exp(u)-1));\n  0.042663 seconds (64.01 k allocations: 1.999 MB, 15.40% gc time)"
},

{
    "location": "about.html#",
    "page": "Nemo",
    "title": "Nemo",
    "category": "page",
    "text": ""
},

{
    "location": "about.html#Nemo-1",
    "page": "Nemo",
    "title": "Nemo",
    "category": "section",
    "text": ""
},

{
    "location": "about.html#What-is-Nemo?-1",
    "page": "Nemo",
    "title": "What is Nemo?",
    "category": "section",
    "text": "Nemo is a computer algebra package for the Julia programming language. Our aim is to provide a highly performant computer algebra package coveringCommutative Algebra\nNumber Theory\nGroup TheoryNemo consists of wrappers of specialised C/C++ libraries:Flint    [http://flintlib.org/]\nArb      [http://fredrikj.net/arb/]\nAntic    [https://github.com/wbhart/antic/]\nSingular [https://www.singular.uni-kl.de/]\nPari     [http://pari.math.u-bordeaux.fr/]It will also eventually provide interfaces to interpreted library code from other computer algebra systems such as Gap and Singular.Nemo also provides implementations of generic algorithms and mathematical data structures. So far the fully recursive constructions includeUnivariate polynomial rings\nPower series rings\nResidue rings (modulo principal ideals)\nFraction fields\nMatrices"
},

{
    "location": "about.html#Why-Julia?-1",
    "page": "Nemo",
    "title": "Why Julia?",
    "category": "section",
    "text": "Julia is a sophisticated, modern programming language which is designed to be both performant and flexible. It was written by mathematicians, for mathematicians.The benefits of Julia includeFamiliar imperative syntax\nJIT compilation (provides near native performance, even for highly generic code)\nREPL console (cuts down on development time)\nParametric types (allows for fast generic constructions over other data types)\nPowerful metaprogramming facilities\nOperator overloading\nMultiple dispatch (dispatch on every argument of a function)\nEfficient native C interface (no wrapper overhead)\nExperimental C++ interface\nDynamic type inference\nBuilt-in bignums\nAble to be embedded in C programs\nHigh performance collection types (dictionaries, iterators, arrays, etc.)\nJupyter support (for web based notebooks)The main benefits for Nemo are the parametric type system and JIT compilation. The former allows us to model many mathematical types, e.g. generic polynomial rings over an arbitrary base ring. The latter speeds up the runtime performance, even of highly generic mathematical procedures."
},

{
    "location": "types.html#",
    "page": "Types in Nemo",
    "title": "Types in Nemo",
    "category": "page",
    "text": ""
},

{
    "location": "types.html#Types-in-Nemo-1",
    "page": "Types in Nemo",
    "title": "Types in Nemo",
    "category": "section",
    "text": "On this page we discuss the type hierarchy in Nemo and a concept known as parents. These details are quite technical and should be skipped or skimmed by new users of Julia/Nemo. Types are almost never dealt with directly when scripting Nemo to do mathematical computations. In contrast, Nemo developers will certainly want to know how we model mathematical objects and the rings, fields, groups, etc. that they belong to in Nemo."
},

{
    "location": "types.html#Introduction-1",
    "page": "Types in Nemo",
    "title": "Introduction",
    "category": "section",
    "text": "Julia provides two levels of types that we make use ofabstract types\nconcrete typesConcrete types are just like the usual types everyone is familiar with from C or C++.Abstract types can be thought of as collections of types. They are used when writing generic functions that should work for any type in the given collection.To write a generic function that accepts any type in a given collection of types, we first create an abstract type. Then we create the individual concrete types that belong to that abstract type. A generic function can then be constructed with a type parameter, T say, similar to a template parameter in C++. The main difference is that we can specify which abstract type our type parameter T must belong to.We use the symbol <: in Julia to determine that a given type belongs to a given abstract type. For example the built-in Julia type Int64 for 64 bit machine integers belongs to the Julia abstract type Integer. Thus Int <: Integer returns true.Here is some Julia code illustrating this with a more complex example. We create an abstract type called Shape and two user defined concrete types square and circle belonging to Shape. We then show how to write methods that accept each of the concrete types and then show how to write a generic function for any type T belonging to the abstract type Shape.Note that in the type definitions of square and circle we specify that those types belong to the abstract type Shape using the <: operator.abstract Shape\n\ntype square <: Shape\n   width::Int\n   border_thickness::Int\nend\n\ntype circle <: Shape\n   centre::Tuple{Int, Int}\n   radius::Int\n   border_thickness::Int\nend\n\nfunction area(s::square)\n   return s.width^2\nend\n\nfunction area(s::circle)\n   return pi*s.radius^2\nend\n\nfunction border_thickness{T <: Shape}(s::T)\n   return s.border_thickness\nend\n\ns = square(3, 1)\nc = circle((3, 4), 2, 2)\n\narea(s)\narea(c)\nborder_thickness(s)\nborder_thickness(c)"
},

{
    "location": "types.html#The-abstract-type-hierarchy-in-Nemo-1",
    "page": "Types in Nemo",
    "title": "The abstract type hierarchy in Nemo",
    "category": "section",
    "text": "Abstract types in Julia can also belong to one another in a hierarchy. For example, the Nemo.Field abstract type belongs to the Nemo.Ring abstract type. An object representing a field in Nemo has type belonging to Nemo.Field. But because we define the inclusion Nemo.Field <: Nemo.Ring in Nemo, the type of such an object also automatically belongs to Nemo.Ring. This means that any generic function in Nemo which is designed to work with ring objects will certainly also work with field objects.In Nemo we also distinguish between the elements of a field, say, and the field itself, and similarly for groups and rings and all other kinds of domains in Nemo. For example, we have an object of type GenPolyRing to model a generic polynomial ring, and elements of that polynomial ring would have type GenPoly. In order to model this distinction between elements and the domains they belong to, Nemo has two main branches in its abstract type hierarchy, as shown in the following diagram. One branch consists of the abstract types for the domains available in Nemo and the other branch is for the abstract types for elements of those domains. (Image: alt text)All objects in Nemo, whether they represent rings, fields, groups, sets, etc. on the one hand, or ring elements, field elements, etc. on the other hand, have concrete types that belong to one of the abstract types shown above."
},

{
    "location": "types.html#Why-types-aren't-enough-1",
    "page": "Types in Nemo",
    "title": "Why types aren't enough",
    "category": "section",
    "text": "Naively, one may expect that rings in Nemo can be modeled as types and their elements as objects with the given type. But there are various reasons why this is not a good model.As an example, consider the ring R = mathbbZnmathbbZ for a multiprecision integer n. If we were to model the ring R as a type, then the type would somehow need to contain the modulus n. This is not possible in Julia, and in fact it is not desirable either.Julia dispatches on type, and each time we call a generic function with different types, a new version of the function is compiled at runtime for performance. But this would be a disaster if we were writing a multimodular algorithm, say. In such an algorithm many rings mathbbZnmathbbZ would be needed and every function we use would be recompiled over and over for each different n. This would result  in a huge delay as the compiler is invoked many times.For this reason, the modulus n needs to be attached to the elements of the ring, not to type associated with those elements.But now we have a problem. How do we create new elements of the ring mathbbZnmathbbZ given only the type? Suppose all rings mathbbZnmathbbZ were represented by the same type Zmod say. How would we create a = 3 pmod7? We could not write a = Zmod(3) since the modulus 7 is not contained in the type Zmod.We could of course use the notation a = Zmod(3, 7), but this would make implementation of generic algorithms very difficult, as they would need to distinguish the case where constructors take a single argument, such as a = ZZ(7) and cases where they take a modulus, such as a = Zmod(3, 7).The way we get around this in Nemo is to have special (singleton) objects that act like types, but are really just ordinary Julia objects. These objects, called parent objects can contain extra information, such as the modulus n. In order to create new elements of mathbbZnmathbbZ as above, we overload the call operator for the parent object, making it callable. Making a parent object callable is exactly analogous to writing a constructor for a type.In the following Nemo example, we create the parent object R corresponding to the ring mathbbZ7mathbbZ. We then create a new element a of this ring by calling the parent object R, just as though R were a type with a constructor accepting an Int parameter. R = ResidueRing(ZZ, 7)\na = R(3)This example creates the element a = 3 pmod7. The important point is that unlike a type, a parent object such as R can contain additional information that a type cannot contain, such as the modulus 7 of the ring in this example, or context objects required by C libraries in other examples."
},

{
    "location": "types.html#More-complex-example-of-parent-objects-1",
    "page": "Types in Nemo",
    "title": "More complex example of parent objects",
    "category": "section",
    "text": "Here is some Julia/Nemo code which constructs a polynomial ring over the integers, a polynomial in that ring and then does some introspection to illustrate the various relations between the objects and types.julia> using Nemo\n\njulia> R, x = ZZ[\"x\"]\n(Univariate Polynomial Ring in x over Integer Ring,x)\n\njulia> f = x^2 + 3x + 1\nx^2+3*x+1\n\njulia> typeof(R)\nNemo.FmpzPolyRing\n\njulia> typeof(f)\nNemo.fmpz_poly\n\njulia> parent(f)\nUnivariate Polynomial Ring in x over Integer Ring\n\njulia> typeof(R) <: PolyRing\ntrue\n\njulia> typeof(f) <: PolyElem\ntrue\n\njulia> parent(f) == R\ntrue"
},

{
    "location": "types.html#Concrete-types-in-Nemo-1",
    "page": "Types in Nemo",
    "title": "Concrete types in Nemo",
    "category": "section",
    "text": "Finally we come to all the concrete types in Nemo. These are of two main kinds: those for generic constructions (e.g. generic polynomials over an arbitrary ring) and those for specific implementations, usually provided by a C library (e.g. polynomials over the integers, provided by Flint).Below we give the type of each kind of element available in Nemo. In parentheses we list the types of their corresponding parent objects. Note that these are the types of the element objects and parent objects respectively, not the abstract types to which these types belong, which the reader can easily guess. For example, fmpz belongs to the abstract type RingElem and FlintIntegerRing belongs to Ring. Similarly Poly{T} belongs to PolyElem whereas PolynomialRing{T} belongs to PolyRing. We also have that fmpz_poly belongs to PolyElem and FmpzPolyRing belongs to PolyRing, and so on.All the generic types are parameterised by a type T which is the type of the elements of the ring they are defined over. Generic\nGenPoly{T} (GenPolyRing{T})\nGenRelSeries{T} (GenRelSeriesRing{T})\nGenRes{T} (GenResRing{T})\nGenFrac{T} (GenFracField{T})\nGenMat{T} (GenMatSpace{T})\nFlint\nfmpz (FlintIntegerRing)\nfmpq (FlintRationalField)\nfq_nmod (FqNmodFiniteField)\nfq (FqFiniteField)\npadic (FlintPadicField)\nfmpz_poly (FmpzPolyRing)\nfmpq_poly (FmpqPolyRing)\nnmod_poly (NmodPolyRing)\nfmpz_mod_poly (FmpzModPolyRing)\nfq_poly (FqPolyRing)\nfq_nmod_poly (FqNmodPolyRing)\nfmpz_rel_series (FmpzRelSeriesRing)\nfmpq_rel_series (FmpqRelSeriesRing)\nfmpz_mod_rel_series (FmpzModRelSeriesRing)\nfq_nmod_rel_series (FqNmodRelSeriesRing)\nfq_rel_series (FqRelSeriesRing)\nfmpz_mat (FmpzMatSpace)\nnmod_mat (NmodMatSpace)\nperm (PermGroup)\nAntic\nnf_elem (AnticNumberField)\nArb\narb (ArbField)\nacb (AcbField)"
},

{
    "location": "constructors.html#",
    "page": "Constructing mathematical objects in Nemo",
    "title": "Constructing mathematical objects in Nemo",
    "category": "page",
    "text": ""
},

{
    "location": "constructors.html#Constructing-mathematical-objects-in-Nemo-1",
    "page": "Constructing mathematical objects in Nemo",
    "title": "Constructing mathematical objects in Nemo",
    "category": "section",
    "text": ""
},

{
    "location": "constructors.html#Constructing-objects-in-Julia-1",
    "page": "Constructing mathematical objects in Nemo",
    "title": "Constructing objects in Julia",
    "category": "section",
    "text": "In Julia, one constructs objects of a given type by calling a type constructor. This is simply a function with the same name as the type itself. For example, to construct a BigInt object in Julia, we simply call the BigInt constructor:n = BigInt(\"1234567898765434567898765434567876543456787654567890\")Julia also uses constructors to convert between types. For example, to convert an Int to a BigInt:m = BigInt(123)"
},

{
    "location": "constructors.html#How-we-construct-objects-in-Nemo-1",
    "page": "Constructing mathematical objects in Nemo",
    "title": "How we construct objects in Nemo",
    "category": "section",
    "text": "As we explained in the previous section, Julia types don't contain enough information to properly model the ring of integers modulo n for a multiprecision modulus n. Instead of using types to construct objects, we use special objects that we refer to as parent objects. They behave a lot like Julia types.Consider the following simple example, to create a Flint multiprecision integer:n = ZZ(\"12345678765456787654567890987654567898765678909876567890\")Here ZZ is not a Julia type, but a callable object. However, for most purposes one can think of such a parent object ZZ as though it were a type."
},

{
    "location": "constructors.html#Constructing-parent-objects-1",
    "page": "Constructing mathematical objects in Nemo",
    "title": "Constructing parent objects",
    "category": "section",
    "text": "For more complicated groups, rings, fields, etc., one first needs to construct the parent object before one can use it to construct element objects.Nemo provides a set of functions for constructing such parent objects. For example, to create a parent object for polynomials over the integers, we use the PolynomialRing parent object constructor.R, x = PolynomialRing(ZZ, \"x\")\nf = x^3 + 3x + 1\ng = R(12)In this example, R is the parent object and we use it to convert the Int value 12 to an element of the polynomial ring mathbbZx."
},

{
    "location": "constructors.html#List-of-parent-object-constructors-1",
    "page": "Constructing mathematical objects in Nemo",
    "title": "List of parent object constructors",
    "category": "section",
    "text": "For convenience, we provide a list of all the parent object constructors in Nemo and explain what domains they represent.Mathematics Nemo constructor\nR = mathbbZ R = ZZ\nR = mathbbQ R = QQ\nR = mathbbF_p^n R, a = FiniteField(p, n, \"a\")\nR = mathbbZnmathbbZ R = ResidueRing(ZZ, n)\nS = Rx S, x = PolynomialRing(R, \"x\")\nS = Rx (to precision n) S, x = PowerSeriesRing(R, n, \"x\")\nS = mboxFrac_R S = FractionField(R)\nS = R(f) S = ResidueRing(R, f)\nS = mboxMat_mtimes n(R) S = MatrixSpace(R, m, n)\nS = mathbbQx(f) S, a = NumberField(f, \"a\")\nO = mathcalO_K S = MaximalOrder(K)\nideal I of O = mathcalO_K I = Ideal(O, gens, ...)"
},

]}
