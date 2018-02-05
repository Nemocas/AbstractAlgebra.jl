var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "AbstractAlgebra.jl",
    "title": "AbstractAlgebra.jl",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#AbstractAlgebra.jl-1",
    "page": "AbstractAlgebra.jl",
    "title": "AbstractAlgebra.jl",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Introduction-1",
    "page": "AbstractAlgebra.jl",
    "title": "Introduction",
    "category": "section",
    "text": "AbstractAlgebra.jl is a computer algebra package for the Julia programming language,  maintained by William Hart, Tommy Hofmann, Claus Fieker and Fredrik Johansson and other interested contributors.Source code\nOnline documentationAbstractAlgebra.jl grew out of the Nemo project after a number of requests from the community for the pure Julia part of Nemo to be split off into a separate project. See the Nemo website for more details about Nemo.Nemo website"
},

{
    "location": "index.html#Features-1",
    "page": "AbstractAlgebra.jl",
    "title": "Features",
    "category": "section",
    "text": "The features of AbstractAlgebra.jl include:Use of Julia multiprecision integers and rationals\nFinite fields (prime order)\nNumber fields\nUnivariate polynomials\nMultivariate polynomials\nRelative and absolute power series\nLaurent series\nFraction fields\nResidue rings, including mathbbZnmathbbZ\nMatrices and linear algebraAll implementations are fully recursive and generic, so that one can build matrices over polynomial rings, over a finite field, for example.AbstractAlgebra.jl also provides a set of abstract types for Groups, Rings, Fields, Modules and elements thereof, which allow external types to be made part of the AbstractAlgebra.jl type hierarchy."
},

{
    "location": "index.html#Installation-1",
    "page": "AbstractAlgebra.jl",
    "title": "Installation",
    "category": "section",
    "text": "To use AbstractAlgebra we require Julia 0.6 or higher. Please see http://julialang.org/downloads for instructions on  how to obtain Julia for your system.At the Julia prompt simply typejulia> Pkg.add(\"AbstractAlgebra\")"
},

{
    "location": "index.html#Quick-start-1",
    "page": "AbstractAlgebra.jl",
    "title": "Quick start",
    "category": "section",
    "text": "Here are some examples of using AbstractAlgebra.jl.This example computes recursive univariate polynomials.using AbstractAlgebra\n\nR, (x, y, z) = PolynomialRing(JuliaZZ, [\"x\", \"y\", \"z\"])\n\nf = x + y + z + 1\n\np = f^20;\n\n@time q = p*(p+1);Here is an example using generic recursive ring constructions.using AbstractAlgebra\n\nR = GF(7)\n\nS, y = PolynomialRing(R, \"y\")\n\nT = ResidueRing(S, y^3 + 3y + 1)\n\nU, z = PolynomialRing(T, \"z\")\n\nf = (3y^2 + y + 2)*z^2 + (2*y^2 + 1)*z + 4y + 3;\n\ng = (7y^2 - y + 7)*z^2 + (3y^2 + 1)*z + 2y + 1;\n\ns = f^4;\n\nt = (s + g)^4;\n\n@time resultant(s, t)Here is an example using matrices.using AbstractAlgebra\n\nR, x = PolynomialRing(JuliaZZ, \"x\")\n\nS = MatrixSpace(R, 10, 10)\n\nM = rand(S, 0:3, -10:10);\n\n@time det(M)And here is an example with power series.using AbstractAlgebra\n\nR, x = JuliaQQ[\"x\"]\n\nS, t = PowerSeriesRing(R, 30, \"t\")\n\nu = t + O(t^100)\n\n@time divexact((u*exp(x*u)), (exp(u)-1));"
},

{
    "location": "constructors.html#",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "Constructing mathematical objects in AbstractAlgebra.jl",
    "category": "page",
    "text": ""
},

{
    "location": "constructors.html#Constructing-mathematical-objects-in-AbstractAlgebra.jl-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "Constructing mathematical objects in AbstractAlgebra.jl",
    "category": "section",
    "text": ""
},

{
    "location": "constructors.html#Constructing-objects-in-Julia-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "Constructing objects in Julia",
    "category": "section",
    "text": "In Julia, one constructs objects of a given type by calling a type constructor. This is simply a function with the same name as the type itself. For example, to construct a  BigInt object in Julia, we simply call the BigInt constructor:n = BigInt(\"1234567898765434567898765434567876543456787654567890\")Julia also uses constructors to convert between types. For example, to convert an Int to a BigInt:m = BigInt(123)"
},

{
    "location": "constructors.html#How-we-construct-objects-in-AbstractAlgebra.jl-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "How we construct objects in AbstractAlgebra.jl",
    "category": "section",
    "text": "As we explain in Appendix A, Julia types don't contain enough information to properly model groups, rings, fields, etc. Instead of using types to construct objects, we use special objects that we refer to as parent objects. They behave a lot like Julia types.Consider the following simple example, to create a multiprecision integer:n = JuliaZZ(\"12345678765456787654567890987654567898765678909876567890\")Here JuliaZZ is not a Julia type, but a callable object. However, for most purposes one can think of such a parent object as though it were a type."
},

{
    "location": "constructors.html#Constructing-parent-objects-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "Constructing parent objects",
    "category": "section",
    "text": "For more complicated groups, rings, fields, etc., one first needs to construct the parent object before one can use it to construct element objects.AbstractAlgebra.jl provides a set of functions for constructing such parent objects. For example, to create a parent object for univariate polynomials over the integers, we use the PolynomialRing parent object constructor.R, x = PolynomialRing(JuliaZZ, \"x\")\nf = x^3 + 3x + 1\ng = R(12)In this example, R is the parent object and we use it to convert the Int value 12 to an element of the polynomial ring mathbbZx."
},

{
    "location": "constructors.html#List-of-parent-object-constructors-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "List of parent object constructors",
    "category": "section",
    "text": "For convenience, we provide a list of all the parent object constructors in AbstractAlgebra.jl and explain what mathematical domains they represent.Mathematics AbstractAlgebra.jl constructor\nR = mathbbZ R = JuliaZZ\nR = mathbbQ R = JuliaQQ\nR = mathbbF_p R = GF(p)\nR = mathbbZnmathbbZ R = ResidueRing(JuliaZZ, n)\nS = Rx S, x = PolynomialRing(R, \"x\")\nS = Rx y S, (x, y) = PolynomialRing(R, [\"x\", \"y\"])\nS = Rx (to precision n) S, x = PowerSeriesRing(R, n, \"x\")\nS = R((x)) (to precision n) S, x = LaurentSeriesRing(R, n, \"x\")\nS = K((x)) (to precision n) S, x = LaurentSeriesField(K, n, \"x\")\nS = mboxFrac_R S = FractionField(R)\nS = R(f) S = ResidueRing(R, f)\nS = mboxMat_mtimes n(R) S = MatrixSpace(R, m, n)\nS = mathbbQx(f) S, a = NumberField(f, \"a\")"
},

{
    "location": "rings.html#",
    "page": "Ring interface",
    "title": "Ring interface",
    "category": "page",
    "text": ""
},

{
    "location": "rings.html#Ring-interface-1",
    "page": "Ring interface",
    "title": "Ring interface",
    "category": "section",
    "text": "AbstractAlgebra.jl generic code makes use of a standardised set of functions which it expects to be implemented for all rings. Here we document this interface. All libraries which want to make use of the generic capabilities of AbstractAlgebra.jl must supply all of the required functionality for their rings.In addition to the required functions, there are also optional functions which can be provided for certain types of rings, e.g. GCD domains or fields, etc. If implemented, these allow the generic code to provide additional functionality for those rings, or in some cases, to select more efficient algorithms."
},

{
    "location": "rings.html#Types-1",
    "page": "Ring interface",
    "title": "Types",
    "category": "section",
    "text": "Most rings must supply two types:a type for the parent object (representing the ring itself)\na type for elements of that ringFor example, the generic univariate polynomial type in AbstractAlgebra.jl provides two  types in generic/GenericTypes.jl: Generic.PolyRing{T} for the parent objects\nGeneric.PolyR{T} for the actual polynomialsThe parent type must belong to AbstractAlgebra.Ring and the element type must belong to AbstractAlgebra.RingElem. Of course, the types may belong to these abstract types transitively, e.g. Poly{T} actually belongs to AbstractAlgebra.PolyElem{T} which in turn belongs to AbstractAlgebra.RingElem.For parameterised rings, we advise that the types of both the parent objects and element objects to be parameterised by the types of the elements of the base ring (see the function base_ring below for a definition).There can be variations on this theme: e.g. in some areas of mathematics there is a notion of a coefficient domain, in which case it may make sense to parameterise all types by the type of elements of this coefficient domain. But note that this may have implications for the ad hoc operators one might like to explicitly implement."
},

{
    "location": "rings.html#Parent-object-caches-1",
    "page": "Ring interface",
    "title": "Parent object caches",
    "category": "section",
    "text": "In many cases, it is desirable to have only one object in the system to represent each ring. This means that if the same ring is constructed twice, elements of the two rings will be compatible as far as arithmetic is concerned.In order to facilitate this, global caches of rings are stored in AbstractAlgebra.jl, usually implemented using dictionaries. For example, the Generic.PolyRing parent objects are looked up in a dictionary PolyID to see if they have been previously defined.Whether these global caches are provided or not, depends on both mathematical and algorithmic considerations. E.g. in the case of number fields, it isn't desirable to identify all number fields with the same defining polynomial, as they may be considered with distinct embeddings into one another. In other cases, identifying whether two rings  are the same may be prohibitively expensive. Generally, it may only make sense algorithmically to identify two rings if they were constructed from identical data.If a global cache is provided, it must be optionally possible to construct the parent objects without caching. This is done by passing a boolean value cached to the inner constructor of the parent object. See generic/GenericTypes.jl` for examples of how to construct and handle such caches."
},

{
    "location": "rings.html#Required-functions-for-all-rings-1",
    "page": "Ring interface",
    "title": "Required functions for all rings",
    "category": "section",
    "text": "In the following, we list all the functions that are required to be provided for rings in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.We give this interface for fictitious types MyParent for the type of the ring parent object R and MyElem for the type of the elements of the ring.Note that generic functions in AbstractAlgebra.jl may not rely on the existence of functions that are not documented here. If they do, those functions will only be available for rings that implement that additional functionality, and should be documented as such."
},

{
    "location": "rings.html#Data-type-and-parent-object-methods-1",
    "page": "Ring interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "parent_type(::Type{MyElem})Returns the type of the corresponding parent object for the given element type. For example, parent_type(Generic.Poly{T}) will return Generic.PolyRing{T}.elem_type(::Type{MyParent})Returns the type of the elements of the ring whose parent object has the given type. This is the inverse of the parent_type function, i.e. elem_type(Generic.PolyRing{T}) will return Generic.Poly{T}.base_ring(R::MyParent)Given a parent object R, representing a ring, this function returns the parent object of any base ring that parameterises this ring. For example, the base ring of the ring of polynomials over the integers would be the integer ring.If the ring is not parameterised by another ring, this function must return Union{}.Note that there is a distinction between a base ring and other kinds of parameters. For example, in the ring mathbbZnmathbbZ, the modulus n is a parameter, but the only base ring is mathbbZ. We consider the ring mathbbZnmathbbZ to have been constructed from the base ring mathbbZ by taking its quotient by a (principal) ideal.parent(f::MyElem)Return the parent object of the given element, i.e. return the ring to which the given element belongs.This is usually stored in a field parent in each ring element. (If the parent objects have mutable struct types, the internal overhead here is just an additional machine  pointer stored in each element of the ring.)For some element types it isn't necessary to append the parent object as a field of every element. This is the case when the parent object can be reconstructed just given the type of the elements. For example, this is the case for the ring of integers and in fact for any ring element type that isn't parameterised or generic in any way.isdomain_type(::Type{MyElem})Returns true if every element of the given element type (which may be parameterised or an abstract type) necessarily has a parent that is an integral domain, otherwise if this cannot be guaranteed, the function returns false. For example, if MyElem was the type of elements of generic residue rings of a polynomial ring, the answer to the question would depend on the modulus of the residue  ring. Therefore isdomain_type would have to return false, since we cannot guarantee that we are dealing with elements of an integral domain in general. But if the given element type was for rational integers, the answer would be true, since every rational integer has as parent the ring of rational integers, which is an integral domain.Note that this function depends only on the type of an element and cannot access information about the object itself, or its parent.isexact_type(::Type{MyElem})Returns true if every element of the given type is represented exactly. For example, p-adic numbers, real and complex floating point numbers and power series are not exact, as we can only represent them in general with finite truncations. Similarly polynomials and matrices over inexact element types are themselves inexact.Integers, rationals, finite fields and polynomials and matrices over them are always exact.Note that MyElem may be parameterised or an abstract type, in which case every element of every type represented by MyElem must be exact, otherwise the function must return false.Base.hash(f::MyElem, h::UInt)Return a hash for the object f of type UInt. This is used as a hopefully cheap way to distinguish objects that differ arithmetically. If the object has components, e.g. the coefficients of a polynomial or elements of a matrix, these should be hashed recursively, passing the same parameter h to all levels. Each component should then be xor'd with h before combining the individual component hashes to give the final hash.The hash functions in AbstractAlgebra.jl usually start from some fixed 64 bit hexadecimal  value that has been picked at random by the library author for that type. That is then truncated to fit a UInt (in case the latter is not 64 bits). This ensures that objects that are the same arithmetically (or that have the same components), but have different types (or structures), are unlikely to hash to the same value.deepcopy_internal(f::MyElem, dict::ObjectIdDict)Return a copy of the given element, recursively copying all components of the object.Obviously the parent, if it is stored in the element, should not be copied. The new element should have precisely the same parent as the old object.For types that cannot self-reference themselves anywhere internally, the dict argument may be ignored.In the case that internal self-references are possible, please consult the Julia documentation on how to implement deepcopy_internal."
},

{
    "location": "rings.html#Constructors-1",
    "page": "Ring interface",
    "title": "Constructors",
    "category": "section",
    "text": "Outer constructors for most AbstractAlgebra types are provided by overloading the call syntax for parent objects. If R is a parent object for a given ring we provide the following constructors.(R::MyParent)()Return the zero object of the given ring.(R::MyParent)(a::Integer)Coerce the given integer into the given ring.(R::MyParent)(a::MyElem)If a belongs to the given ring, the function returns it (without making a copy). Otherwise an error is thrown.For parameterised rings we also require a function to coerce from the base ring into the parent ring.(R::MyParent{T})(a::T) where T <: AbstractAlgebra.RingElemCoerce a into the ring R if a belongs to the base ring of R."
},

{
    "location": "rings.html#Basic-manipulation-of-rings-and-elements-1",
    "page": "Ring interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "zero(R::MyParent)Return the zero element of the given ring.one(R::MyParent)Return the multiplicative identity of the given ring.iszero(f::MyElem)Return true if the given element is the zero element of the ring it belongs to.isone(f::MyElem)Return true if the given element is the multiplicative identity of the ring it belongs to."
},

{
    "location": "rings.html#Canonicalisation-1",
    "page": "Ring interface",
    "title": "Canonicalisation",
    "category": "section",
    "text": "canonical_unit(f::MyElem)When fractions are created with two elements of the given type, it is nice to be able to represent them in some kind of canonical form. This is of course not always possible. But for example, fractions of integers can be canonicalised by first removing any common factors of the numerator and denominator, then making the denominator positive.In AbstractAlgebra.jl, the denominator would be made positive by dividing both the numerator and denominator by the canonical unit of the denominator. For a negative denominator, this would be -1.For elements of a field, canonical_unit simply returns the element itself. In general, canonical_unit of an invertible element should be that element. Finally, if a = ub we should have the identity canonical_unit(a) = canonical_unit(u)*canonical_unit(b).For some rings, it is completely impractical to implement this function, in which case it may return 1 in the given ring. The function must however always exist, and always return an element of the ring."
},

{
    "location": "rings.html#String-I/O-1",
    "page": "Ring interface",
    "title": "String I/O",
    "category": "section",
    "text": "show(io::IO, R::MyParent)This should print (to the given IO object), an English description of the parent ring. If the ring is parameterised, it can call the corresponding show function for any rings it depends on.show(io::IO, f::MyElem)This should print a human readable, textual representation of the object (to the given IO object). It can recursively call the corresponding show functions for any of its components.It may be necessary in some cases to print parentheses around components of f or to print signs of components. For these, the following functions will exist for each component or component type.needs_parentheses(f::MyElem)Should returns true if parentheses are needed around this object when printed, e.g. as a coefficient of a polynomial. As an example, non-constant polynomials would need such parentheses if used as coefficients of another polynomial.isnegative(f::MyElem)When printing polynomials, a + sign is usually inserted automatically between terms of the polynomial. However, this is not desirable if the coefficient is negative and that negative sign is already printed when the coefficient is printed.This function must return true if a negative sign would already be prepended when f is printed. This suppresses the automatic printing of a + sign by polynomial printing functions that are printing f as a coefficient of a term.Note that if needs_parentheses returns true for f, then isnegative should always return false for that f, since an automatic + will need to be printed in front of a coefficient that is printed with parentheses.show_minus_one(::Type{MyElem})When printing polynomials, we prefer to print x rather than 1*x if the degree 1 term has coefficient 1. This can be taken care of without any special support.However, we also prefer to print -x rather than -1*x. This requires special support, since -1 in some rings is not printed as -1 (e.g. -1 in mathbbZ3mathbbZ might be printed as 2). In such rings, show_minus_one should return true.If show_minus_one returns true, polynomial printing functions will not print -x for terms of degree 1 with coefficient -1, but will use the printing function of the given type to print the coefficient in that case."
},

{
    "location": "rings.html#Unary-operations-1",
    "page": "Ring interface",
    "title": "Unary operations",
    "category": "section",
    "text": "-(f::MyElem)Returns -f."
},

{
    "location": "rings.html#Binary-operations-1",
    "page": "Ring interface",
    "title": "Binary operations",
    "category": "section",
    "text": "+(f::MyElem, g::MyElem)\n-(f::MyElem, g::MyElem)\n*(f::MyElem, g::MyElem)Returns f + g, f - g or fg, respectively."
},

{
    "location": "rings.html#Comparison-1",
    "page": "Ring interface",
    "title": "Comparison",
    "category": "section",
    "text": "==(f::MyElem, g::MyElem)Returns true if f and g are arithmetically equal. In the case where the two elements are inexact, the function returns true if they agree to the minimum precision of the two.isequal(f::MyElem, g::MyElem)For exact rings, this should return the same thing as == above. For inexact rings, this returns true only if the two elements are arithmetically equal and have the same precision."
},

{
    "location": "rings.html#Powering-1",
    "page": "Ring interface",
    "title": "Powering",
    "category": "section",
    "text": "^(f::MyElem, e::Int)Return f^e. The function should throw a DomainError() if negative exponents don't make sense but are passed to the function."
},

{
    "location": "rings.html#Exact-division-1",
    "page": "Ring interface",
    "title": "Exact division",
    "category": "section",
    "text": "divexact(f::MyElem, g::MyElem)Returns fg, though note that Julia uses / for floating point division. Here we mean exact division in the ring, i.e. return q such that f = gq. A DivideError() should be thrown if g is zero. If no exact quotient exists or an impossible inverse is unavoidably encountered, an error should be thrown."
},

{
    "location": "rings.html#Unsafe-operators-1",
    "page": "Ring interface",
    "title": "Unsafe operators",
    "category": "section",
    "text": "To speed up polynomial and matrix arithmetic, it sometimes makes sense to mutate values in place rather than replace them with a newly created object every time they are modified.For this purpose, certain mutating operators are required. In order to support immutable types (struct in Julia) and systems that don't have in-place operators, all unsafe operators must return the (ostensibly) mutated value. Only the returned value is used in computations, so this lifts the requirement that the unsafe operators actually mutate the value.Note the exclamation point is a convention, which indicates that the object may be mutated in-place.To make use of these functions, one must be certain that no other references are held to the object being mutated, otherwise those values will also be changed!The results of deepcopy and all arithmetic operations, including powering and division can be assumed to be new objects without other references being held, as can objects returned from constructors.Note that R(a) where R is the ring a belongs to, does not create a new value. For this case, use deepcopy(a).zero!(f::MyElem)Set the value f to zero in place. Return the mutated value.mul!(c::MyElem, a::MyElem, b::MyElem)Set c to the value ab in place. Return the mutated value. Aliasing is permitted.add!(c::MyElem, a::MyElem, b::MyElem)Set c to the value a + b in place. Return the mutated value. Aliasing is permitted."
},

{
    "location": "rings.html#Random-generation-1",
    "page": "Ring interface",
    "title": "Random generation",
    "category": "section",
    "text": "The random functions are only used for test code to generate test data. They therefore don't need to provide any guarantees on uniformity, and in fact, test values that are known to be a good source of corner cases can be supplied.rand(R::MyParent, v...)Returns a random element in the given ring of the specified size.There can be as many arguments as is necessary to specify the size of the test example which is being produced."
},

{
    "location": "rings.html#Promotion-rules-1",
    "page": "Ring interface",
    "title": "Promotion rules",
    "category": "section",
    "text": "In order for AbstractAlgebra to be able to automatically coerce up towers of rings, certain promotion rules must be defined. For every ring, one wants to be able to coerce integers into the ring. And for any ring constructed over a base ring, one would like to be able to coerce from the base ring into the ring.The promotion rules look a bit different depending on whether the element type is parameterised or not and whether it is built on a base ring.For ring element types MyElem that are neither parameterised, not built over a base ring, the promotion rules can be defined as follows:promote_rule(::Type{MyElem}, ::Type{T}) where {T <: Integer} = MyElemFor ring element types MyType that aren't parameterised, but which have a base ring with concrete element type T the promotion rules can be defined as follows:promote_rule(::Type{MyElem}, ::Type{U}) where U <: Integer = MyElempromote_rule(::Type{MyElem}, ::Type{T}) = MyElemFor ring element types MyElem{T} that are parameterised by the type of elements of the base ring, the promotion rules can be defined as follows:promote_rule(::Type{MyElem{T}}, ::Type{MyElem{T}}) where T <: RingElement = MyElem{T}\njulia function promote_rule(::Type{MyElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingEle ment}    promote_rule(T, U) == T ? MyElem{T} : Union{} end\n## Required functionality for inexact rings\n\n### Approximation (floating point and ball arithmetic only)\njulia isapprox(f::MyElem, g::MyElem; atol::Real=sqrt(eps()))\nThis is used by test code that uses rings involving floating point or ball arithmetic.\nThe function should return `true` if all components of $f$ and $g$ are equal to\nwithin the square root of the Julia epsilon, since numerical noise may make an exact\ncomparison impossible.\n\nFor parameterised rings over an inexact ring, we also require the following ad hoc\napproximation functionality.\njulia isapprox(f::MyElem{T}, g::T; atol::Real=sqrt(eps())) where T <: AbstractAlgebra.RingElemjulia isapprox(f::T, g::MyElem{T}; atol::Real=sqrt(eps())) where T <: AbstractAlgebra.RingElem\nThese notionally coerce the element of the base ring into the parameterised ring and do\na full comparison.\n\n## Optional functionality\n\nSome functionality is difficult or impossible to implement for all rings in the system.\nIf it is provided, additional functionality or performance may become available. Here\nis a list of all functions that are considered optional and can't be relied on by\ngeneric functions in the AbstractAlgebra Ring interface.\n\nIt may be that no algorithm, or no efficient algorithm is known to implement these\nfunctions. As these functions are optional, they do not need to exist. Julia will\nalready inform the user that the function has not been implemented if it is called but\ndoesn't exist.\n\n### Optional basic manipulation functionality\njulia isunit(f::MyElem)\nReturn `true` if the given element is a unit in the ring it belongs to. \n \n### Optional binary ad hoc operators\n\nBy default, ad hoc operations are handled by AbstractAlgebra.jl if they are not defined\nexplicitly, by coercing both operands into the same ring and then performing the\nrequired operation.\n\nIn some cases, e.g. for matrices, this leads to very inefficient behaviour. In such\ncases, it is advised to implement some of these operators explicitly.\n\nIt can occasionally be worth adding a separate set of ad hoc binary operators for the\ntype `Int`, if this can be done more efficiently than for arbitrary Julia Integer types.\njulia +(f::MyElem, c::Integer) -(f::MyElem, c::Integer) *(f::MyElem, c::Integer)julia +(c::Integer, f::MyElem) -(c::Integer, f::MyElem) *(c::Integer, f::MyElem)\nFor parameterised types, it is also sometimes more performant to provide explicit ad\nhoc operators with elements of the base ring.\njulia +(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem -(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem *(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElemjulia +(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem -(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem *(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem\n### Optional ad hoc comparisons\njulia ==(f::MyElem, c::Integer)julia ==(c::Integer, f::MyElem)julia ==(f::MyElem{T}, c:T) where T <: AbstractAlgebra.RingElemjulia ==(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem\n### Optional ad hoc exact division functions\njulia divexact(a::MyType{T}, b::T) where T <: AbstractAlgebra.RingElemjulia divexact(a::MyType, b::Integer)\n### Optional powering functions\njulia ^(f::MyElem, e::BigInt)\nIn case $f$ cannot explode in size when powered by a very large integer, and it is\npractical to do so, one may provide this function to support powering with `BigInt`\nexponents.\n\n### Optional unsafe operators\njulia addmul!(c::MyElem, a::MyElem, b::MyElem, t::MyElem) ```Set c = c + ab in-place. Return the mutated value. The value t should be a temporary of the same type as a, b and c, which can be used arbitrarily by the implementation to speed up the computation. Aliasing between a, b and c is  permitted."
},

{
    "location": "euclidean.html#",
    "page": "Euclidean ring interface",
    "title": "Euclidean ring interface",
    "category": "page",
    "text": ""
},

{
    "location": "euclidean.html#Euclidean-ring-interface-1",
    "page": "Euclidean ring interface",
    "title": "Euclidean ring interface",
    "category": "section",
    "text": "If a ring provides a meaningful Euclidean structure such that a useful Euclidean remainder can be computed practically, various additional functionality is provided by AbstractAlgebra.jl for those rings. This functionality depends on the following functions existing.mod(f::MyElem, g::MyElem)Returns the Euclidean remainder of f by g. A DivideError() should be thrown if g is zero. An error should be thrown if an impossible inverse is encountered.divrem(f::MyElem, g::MyElem)Returns a pair q, r consisting of the Euclidean quotient and remainder of f by g. A DivideError should be thrown if g is zero. An error should be thrown if an impossible inverse is encountered.div(f::MyElem, g::MyElem)Returns the Euclidean quotient of f by g. A DivideError should be thrown if g is zero. An error should be thrown if an impossible inverse is encountered.mulmod(f::MyElem, g::MyElem, m::MyElem)Returns fg pmodm.powmod(f::MyElem, e::Int, m::MyElem)Returns f^e pmodm.invmod(f::MyElem, m::MyElem)Returns the inverse of f modulo m. If such an inverse doesn't exist, an impossible inverse error should be thrown.divides(f::MyElem, g::MyElem)Returns a pair, flag, q, where flag is set to true if g divides f, in which case the quotient is set to the quotient, or flag is set to false and the quotient is set to zero in the same ring as f and g.remove(f::MyElem, p::MyElem)Returns a pair v, q where p^v is the highest power of p dividing f, and q is the cofactor after f is divided by this power.valuation(f::MyElem, p::MyElem)Returns v where p^v is the highest power of p dividing f.gcd(f::MyElem, g::MyElem)Returns a greatest common divisor of f and g.lcm(f::MyElem, g::MyElem)Returns fggcd(f g) if either f or g is not zero, otherwise it throws a DivideError().gcdx(f::MyElem, g::MyElem)Returns a triple d, s, t such that d = gcd(f g) and d = sf + tg, with s reduced modulo g and t reduced modulo f.gcdinv(f::MyElem, g::MyElem)Returns a tuple d, s such that d = gcd(f g) and s = (fd)^-1 pmodgd. Note that d = 1 iff f is invertible modulo g, in which case s = f^-1 pmodg."
},

{
    "location": "polynomial_rings.html#",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Univariate Polynomial Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "polynomial_rings.html#Univariate-Polynomial-Ring-Interface-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Univariate Polynomial Ring Interface",
    "category": "section",
    "text": "Univariate polynomial rings are supported in AbstractAlgebra, and in addition to the standard Ring interface, numerous additional functions are required to be present for univariate polynomial rings.Univariate polynomial rings over a field are also Euclidean and therefore such rings must implement the Euclidean interface.Since a sparse distributed multivariate format can generally also handle sparse univariate polynomials, the univariate polynomial interface is designed around the assumption that they are dense. This is not a requirement, but it may be easier to use the multivariate interface for sparse univariate types."
},

{
    "location": "polynomial_rings.html#Types-and-parents-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra provides two abstract types for polynomial rings and their elements:PolyRing{T} is the abstract type for univariate polynomial ring parent types\nPolyElem{T} is the abstract type for univariate polynomial typesWe have that PolyRing{T} <: AbstractAlgebra.Ring and PolyElem{T} <: AbstractAlgebra.RingElem.Note that both abstract types are parameterised. The type T should usually be the type of elements of the coefficient ring of the polynomial ring. For example, in the case of mathbbZx the type T would be the type of an integer, e.g. BigInt.If the parent object for such a ring has type MyZX and polynomials in that ring have type MyZXPoly then one would have:MyZX <: PolyRing{BigInt}\nMyZXPoly <: PolyElem{BigInt}Polynomial rings should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Polynomial rings should at least be distinguished based on their base (coefficient) ring. But if they have the same base ring and symbol (for their variable/generator), they should certainly have the same parent object.See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "polynomial_rings.html#Required-functionality-for-univariate-polynomials-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Required functionality for univariate polynomials",
    "category": "section",
    "text": "In addition to the required functionality for the Ring interface (and in the case of polynomials over a field, the Euclidean Ring interface), the Polynomial Ring interface has the following required functions.We suppose that R is a fictitious base ring (coefficient ring) and that S is a univariate polynomial ring over R (i.e. S = Rx) with parent object S of type MyPolyRing{T}. We also assume the polynomials in the ring have type MyPoly{T}, where T is the type of elements of the base (coefficient) ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem."
},

{
    "location": "polynomial_rings.html#Constructors-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Constructors",
    "category": "section",
    "text": "In addition to the standard constructors, the following constructors, taking an array of coefficients, must be available. In each case, the array may have trailing zero entries, but the resulting polynomial should usually be normalised (of length zero, or with nonzero leading coefficient).(S::MyPolyRing{T})(A::Array{T, 1}) where T <: AbstractAlgebra.RingElemCreate the polynomial in the given ring whose degree i coefficient is given by A[i].(S::MyPolyRing{T})(A::Array{U, 1}) where T <: AbstractAlgebra.RingElem, U <: AbstractAlgebra.RingElemCreate the polynomial in the given ring whose degree i coefficient is given by A[i]. The elements of the array are assumed to be able to be coerced into the base ring R.(S::MyPolyRing{T})(A::Array{U, 1}) where T <: AbstractAlgebra.RingElem, U <: IntegerCreate the polynomial in the given ring whose degree i coefficient is given by A[i].It may be desirable to have a additional version of the function that accepts an array of Julia Int values  if this can be done more efficiently.** Examples **S, x = PolynomialRing(JuliaQQ, \"x\")\n\nf = S(Rational{BigInt}[2, 3, 1])\ng = S(BigInt[1, 0, 4])\nh = S([4, 7, 2, 9])"
},

{
    "location": "polynomial_rings.html#Data-type-and-parent-object-methods-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "var(S::MyPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn a Symbol representing the variable (generator) of the polynomial ring. Note that this is a Symbol not a String, though its string value will usually be used when printing polynomials.vars(S::MyPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn the array [s] where s	 is aSymbol` representing the variable of the given polynomial ring. This is provided for uniformity with the multivariate interface, where there is more than one variable, and hence an array of symbols.** Examples **S, x = PolynomialRing(JuliaQQ, \"x\")\n\nvsym = var(S)\nV = vars(S)"
},

{
    "location": "polynomial_rings.html#Basic-manipulation-of-rings-and-elements-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "length(f::MyPoly{T}) where T <: AbstractAlgebra.RingElemReturn the length of the given polynomial. The length of the zero polynomial is defined to be 0, otherwise the length is the degree plus 1. The return value should be of type Int.set_length!(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemThis function must zero any coefficients beyond the requested length n and then set the length of the polynomial to n. This function does not need to normalise the polynomial and is not useful to the user, but is used extensively by the AbstractAlgebra generic functionality.This function mutates the existing polynomial in-place, but does not return the polynomial.coeff(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemReturn the coefficient of the polynomial f of degree n. If n is larger than the degree of the polynomial, it should return zero in the coefficient ring. setcoeff!(f::MyPoly{T}, n::Int, a::T) where T <: AbstractAlgebra.RingElemSet the degree n coefficient of f to a. This mutates the polynomial in-place if possible and returns the mutated polynomial (so that immutable types can also be supported). The function must not assume that the polynomial already has space for n + 1 coefficients. The polynomial must be resized if this is not the case.Note that this function is not required to normalise the polynomial and is not necessarily useful to the user, but is used extensively by the generic functionality in AbstractAlgebra.jl. It is for setting raw coefficients in the representation.normalise(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemGiven a polynomial whose length is currently n, including any leading zero coefficients, return the length of the normalised polynomial (either zero of the length of the polynomial with nonzero leading coefficient). Note that the function does not actually perform the normalisation.fit!(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemEnsure that the polynomial f internally has space for n coefficients. This function must mutate the function in-place if it is mutable. It does not return the mutated polynomial. Immutable types can still be supported by defining this function to do nothing.Some interfaces for C polynomial types automatically manage the internal allocation of polynomials in every function that can be called on them. Explicit adjustment by the generic code in AbstractAlgebra.jl is not required. In such cases, this function can also be defined to do nothing.** Examples **S, x = PolynomialRing(JuliaZZ, \"x\")\n\nf = x^3 + 3x + 1\ng = S(BigInt[1, 2, 0, 1, 0, 0, 0]);\n\nn = length(f)\nc = coeff(f, 1)\nset_length!(g, normalise(g, 7))\ng = setcoeff!(g, 2, BigInt(11))\nfit!(g, 8)\ng = setcoeff!(g, 7, BigInt(4))\n"
},

{
    "location": "polynomial_rings.html#Optional-functionality-for-polynomial-rings-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Optional functionality for polynomial rings",
    "category": "section",
    "text": "Sometimes parts of the Euclidean Ring interface can and should be implemented for polynomials over a ring that is not necessarily a field.When divisibility testing can be implemented for a polynomial ring over a field, it  should be possible to implement the following functions from the Euclidean Ring interface:divides\nremove\nvaluationWhen the given polynomial ring is a GCD domain, with an effective GCD algorithm, it may be possible to implement the following functions:gcd\nlcmPolynomial rings can optionally implement any part of the generic univariate polynomial functionality provided by AbstractAlgebra.jl, using the same interface. Obviously additional functionality can also be added to that provided by AbstractAlgebra.jl on an ad hoc basis."
},

{
    "location": "types.html#",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "Appendix A: Types in AbstractAlgebra.jl",
    "category": "page",
    "text": ""
},

{
    "location": "types.html#Appendix-A:-Types-in-AbstractAlgebra.jl-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "Appendix A: Types in AbstractAlgebra.jl",
    "category": "section",
    "text": "On this page we discuss the abstract type hierarchy in AbstractAlgebra.jl and objects known as parents which contain additional information about groups, rings, fields and modules, etc., that can't be stored in types alone.These details are technical and can be skipped or skimmed by new users of  Julia/AbstractAlgebra.jl. Types are almost never dealt with directly when scripting  AbstractAlgebra.jl to do mathematical computations. In contrast, AbstractAlgebra.jl developers will want to know how we model mathematical objects and their rings, fields, groups, etc."
},

{
    "location": "types.html#The-abstract-type-hierarchy-in-AbstractAlgebra.jl-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "The abstract type hierarchy in AbstractAlgebra.jl",
    "category": "section",
    "text": "Abstract types in Julia can also belong to one another in a hierarchy. We make use of such a hierarchy to organise the kinds of mathematical objects in AbstractAlgebra.jl.For example, the AbstractAlgebra.Field abstract type belongs to the  AbstractAlgebra.Ring abstract type. In practice, this means that any generic function in AbstractAlgebra.jl which is designed to work with ring objects will also work with field objects.In AbstractAlgebra.jl we also distinguish between the elements of a field, say, and the field itself.For example, we have an object of type Generic.PolyRing to model a generic polynomial ring, and elements of that polynomial ring would have type Generic.Poly. For this purpose, we also have a hierarchy of abstract types, such as FieldElem, that the types of element objects can belong to.(Image: alt text)"
},

{
    "location": "types.html#Why-types-aren't-enough-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "Why types aren't enough",
    "category": "section",
    "text": "Naively, one might have expected that rings in AbstractAlgebra.jl could be modeled as types and their elements as objects with the given type. But there are various reasons why this is not a good model.Consider the ring R = mathbbZnmathbbZ for a multiprecision integer n. If we were to model the ring R as a type, then the type would somehow need to contain the modulus n. This is not possible in Julia, and in fact it is not desirable, since the compiler would then recompile all the associated functions every time a different modulus n was used.We could attach the modulus n to the objects representing elements of the ring, rather than their type.But now we cannot create new elements of the ring mathbbZnmathbbZ given only their type, since the type no longer contains the modulus n.Instead, the way we get around this in AbstractAlgebra.jl is to have special (singleton) objects that act like types, but are really just ordinary Julia objects. These objects, called parent objects can contain extra information, such as the modulus n. In order to create new elements of mathbbZnmathbbZ as above, we overload the call operator for the parent object.In the following AbstractAlgebra.jl example, we create the parent object R corresponding to the ring mathbbZ7mathbbZ. We then create a new element a of this ring by calling the parent object R.R = ResidueRing(JuliaZZ, 7)\na = R(3)Here, R is the parent object, containing the modulus 7. So this example creates  the element a = 3 pmod7."
},

{
    "location": "types.html#More-complex-example-of-parent-objects-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "More complex example of parent objects",
    "category": "section",
    "text": "Here is some Julia/AbstractAlgebra.jl code which constructs a polynomial ring over the integers, a polynomial in that ring and then does some introspection to illustrate the various relations between the objects and types.using AbstractAlgebra\n\nR, x = JuliaZZ[\"x\"]\n\nf = x^2 + 3x + 1\n\ntypeof(R) <: PolyRing\n\ntypeof(f) <: PolyElem\n\nparent(f) == R"
},

{
    "location": "types.html#Concrete-types-in-AbstractAlgebra.jl-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "Concrete types in AbstractAlgebra.jl",
    "category": "section",
    "text": "Here we give a list of the concrete types in AbstractAlgebra.jl.In parentheses we put the types of the corresponding parent objects.gfelem{<:Integer} (GFField{<:Integer})We also think of various Julia types as though they were AbstractAlgebra.jl types:BigInt (Integers{BigInt})\nRational{BigInt} (Rationals{BigInt})Then there are various types for generic constructions over a base ring. They are all parameterised by a type T which is the type of the elements of the base ring they are defined over. Generic.Poly{T} (Generic.PolyRing{T})\nGeneric.MPoly{T} (Generic.MPolyRing{T})\nGeneric.RelSeries{T} (Generic.RelSeriesRing{T})\nGeneric.AbsSeries{T} (Generic.AbsSeriesRing{T})\nGeneric.LaurentSeriesRingElem{T} (Generic.LaurentSeriesRing{T})\nGeneric.LaurentSeriesFieldElem{T} (Generic.LaurentSeriesField{T})\nGeneric.Res{T} (Generic.ResRing{T})\nGeneric.Frac{T} (Generic.FracField{T})\nGeneric.Mat{T} (Generic.MatSpace{T})"
},

]}
