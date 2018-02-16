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
    "text": "The features of AbstractAlgebra.jl include:Use of Julia multiprecision integers and rationals\nFinite fields (prime order, naive implementation only)\nNumber fields (naive implementation only)\nUnivariate polynomials\nMultivariate polynomials\nRelative and absolute power series\nLaurent series\nFraction fields\nResidue rings, including mathbbZnmathbbZ\nMatrices and linear algebraAll implementations are fully recursive and generic, so that one can build matrices over polynomial rings, over a finite field, for example.AbstractAlgebra.jl also provides a set of abstract types for Groups, Rings, Fields, Modules and elements thereof, which allow external types to be made part of the AbstractAlgebra.jl type hierarchy."
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
    "text": "Here are some examples of using AbstractAlgebra.jl.This example computes recursive univariate polynomials.using AbstractAlgebra\n\nR, (x, y, z) = PolynomialRing(JuliaZZ, [\"x\", \"y\", \"z\"])\n\nf = x + y + z + 1\n\np = f^20;\n\n@time q = p*(p+1);Here is an example using generic recursive ring constructions.using AbstractAlgebra\n\nR = GF(7)\n\nS, y = PolynomialRing(R, \"y\")\n\nT = ResidueRing(S, y^3 + 3y + 1)\n\nU, z = PolynomialRing(T, \"z\")\n\nf = (3y^2 + y + 2)*z^2 + (2*y^2 + 1)*z + 4y + 3;\n\ng = (7y^2 - y + 7)*z^2 + (3y^2 + 1)*z + 2y + 1;\n\ns = f^4;\n\nt = (s + g)^4;\n\n@time resultant(s, t)Here is an example using matrices.using AbstractAlgebra\n\nR, x = PolynomialRing(JuliaZZ, \"x\")\n\nS = MatrixSpace(R, 10, 10)\n\nM = rand(S, 0:3, -10:10);\n\n@time det(M)And here is an example with power series.using AbstractAlgebra\n\nR, x = QQ[\"x\"]\n\nS, t = PowerSeriesRing(R, 30, \"t\")\n\nu = t + O(t^100)\n\n@time divexact((u*exp(x*u)), (exp(u)-1));"
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
    "text": "As we explain in Appendix A, Julia types don't contain enough information to properly model groups, rings, fields, etc. Instead of using types to construct objects, we use special objects that we refer to as parent objects. They behave a lot like Julia types.Consider the following simple example, to create a multiprecision integer:n = ZZ(\"12345678765456787654567890987654567898765678909876567890\")Here ZZ is not a Julia type, but a callable object. However, for most purposes one can think of such a parent object as though it were a type."
},

{
    "location": "constructors.html#Constructing-parent-objects-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "Constructing parent objects",
    "category": "section",
    "text": "For more complicated groups, rings, fields, etc., one first needs to construct the parent object before one can use it to construct element objects.AbstractAlgebra.jl provides a set of functions for constructing such parent objects. For example, to create a parent object for univariate polynomials over the integers, we use the PolynomialRing parent object constructor.R, x = PolynomialRing(ZZ, \"x\")\nf = x^3 + 3x + 1\ng = R(12)In this example, R is the parent object and we use it to convert the Int value 12 to an element of the polynomial ring mathbbZx."
},

{
    "location": "constructors.html#List-of-parent-object-constructors-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "List of parent object constructors",
    "category": "section",
    "text": "For convenience, we provide a list of all the parent object constructors in AbstractAlgebra.jl and explain what mathematical domains they represent.Mathematics AbstractAlgebra.jl constructor\nR = mathbbZ R = ZZ\nR = mathbbQ R = QQ\nR = mathbbF_p R = GF(p)\nR = mathbbZnmathbbZ R = ResidueRing(ZZ, n)\nS = Rx S, x = PolynomialRing(R, \"x\")\nS = Rx y S, (x, y) = PolynomialRing(R, [\"x\", \"y\"])\nS = Rx (to precision n) S, x = PowerSeriesRing(R, n, \"x\")\nS = R((x)) (to precision n) S, x = LaurentSeriesRing(R, n, \"x\")\nS = K((x)) (to precision n) S, x = LaurentSeriesField(K, n, \"x\")\nS = mboxFrac_R S = FractionField(R)\nS = R(f) S = ResidueRing(R, f)\nS = mboxMat_mtimes n(R) S = MatrixSpace(R, m, n)\nS = mathbbQx(f) S, a = NumberField(f, \"a\")"
},

{
    "location": "rings.html#",
    "page": "Ring Interface",
    "title": "Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "rings.html#Ring-Interface-1",
    "page": "Ring Interface",
    "title": "Ring Interface",
    "category": "section",
    "text": "AbstractAlgebra.jl generic code makes use of a standardised set of functions which it expects to be implemented for all rings. Here we document this interface. All libraries which want to make use of the generic capabilities of AbstractAlgebra.jl must supply all of the required functionality for their rings.In addition to the required functions, there are also optional functions which can be provided for certain types of rings, e.g. GCD domains or fields, etc. If implemented, these allow the generic code to provide additional functionality for those rings, or in some cases, to select more efficient algorithms."
},

{
    "location": "rings.html#Types-1",
    "page": "Ring Interface",
    "title": "Types",
    "category": "section",
    "text": "Most rings must supply two types:a type for the parent object (representing the ring itself)\na type for elements of that ringFor example, the generic univariate polynomial type in AbstractAlgebra.jl provides two  types in generic/GenericTypes.jl: Generic.PolyRing{T} for the parent objects\nGeneric.PolyR{T} for the actual polynomialsThe parent type must belong to AbstractAlgebra.Ring and the element type must belong to AbstractAlgebra.RingElem. Of course, the types may belong to these abstract types transitively, e.g. Poly{T} actually belongs to AbstractAlgebra.PolyElem{T} which in turn belongs to AbstractAlgebra.RingElem.For parameterised rings, we advise that the types of both the parent objects and element objects to be parameterised by the types of the elements of the base ring (see the function base_ring below for a definition).There can be variations on this theme: e.g. in some areas of mathematics there is a notion of a coefficient domain, in which case it may make sense to parameterise all types by the type of elements of this coefficient domain. But note that this may have implications for the ad hoc operators one might like to explicitly implement."
},

{
    "location": "rings.html#Parent-object-caches-1",
    "page": "Ring Interface",
    "title": "Parent object caches",
    "category": "section",
    "text": "In many cases, it is desirable to have only one object in the system to represent each ring. This means that if the same ring is constructed twice, elements of the two rings will be compatible as far as arithmetic is concerned.In order to facilitate this, global caches of rings are stored in AbstractAlgebra.jl, usually implemented using dictionaries. For example, the Generic.PolyRing parent objects are looked up in a dictionary PolyID to see if they have been previously defined.Whether these global caches are provided or not, depends on both mathematical and algorithmic considerations. E.g. in the case of number fields, it isn't desirable to identify all number fields with the same defining polynomial, as they may be considered with distinct embeddings into one another. In other cases, identifying whether two rings  are the same may be prohibitively expensive. Generally, it may only make sense algorithmically to identify two rings if they were constructed from identical data.If a global cache is provided, it must be optionally possible to construct the parent objects without caching. This is done by passing a boolean value cached to the inner constructor of the parent object. See generic/GenericTypes.jl` for examples of how to construct and handle such caches."
},

{
    "location": "rings.html#Required-functions-for-all-rings-1",
    "page": "Ring Interface",
    "title": "Required functions for all rings",
    "category": "section",
    "text": "In the following, we list all the functions that are required to be provided for rings in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.We give this interface for fictitious types MyParent for the type of the ring parent object R and MyElem for the type of the elements of the ring.Note that generic functions in AbstractAlgebra.jl may not rely on the existence of functions that are not documented here. If they do, those functions will only be available for rings that implement that additional functionality, and should be documented as such."
},

{
    "location": "rings.html#Data-type-and-parent-object-methods-1",
    "page": "Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "parent_type(::Type{MyElem})Returns the type of the corresponding parent object for the given element type. For example, parent_type(Generic.Poly{T}) will return Generic.PolyRing{T}.elem_type(::Type{MyParent})Returns the type of the elements of the ring whose parent object has the given type. This is the inverse of the parent_type function, i.e. elem_type(Generic.PolyRing{T}) will return Generic.Poly{T}.base_ring(R::MyParent)Given a parent object R, representing a ring, this function returns the parent object of any base ring that parameterises this ring. For example, the base ring of the ring of polynomials over the integers would be the integer ring.If the ring is not parameterised by another ring, this function must return Union{}.Note that there is a distinction between a base ring and other kinds of parameters. For example, in the ring mathbbZnmathbbZ, the modulus n is a parameter, but the only base ring is mathbbZ. We consider the ring mathbbZnmathbbZ to have been constructed from the base ring mathbbZ by taking its quotient by a (principal) ideal.parent(f::MyElem)Return the parent object of the given element, i.e. return the ring to which the given element belongs.This is usually stored in a field parent in each ring element. (If the parent objects have mutable struct types, the internal overhead here is just an additional machine  pointer stored in each element of the ring.)For some element types it isn't necessary to append the parent object as a field of every element. This is the case when the parent object can be reconstructed just given the type of the elements. For example, this is the case for the ring of integers and in fact for any ring element type that isn't parameterised or generic in any way.isdomain_type(::Type{MyElem})Returns true if every element of the given element type (which may be parameterised or an abstract type) necessarily has a parent that is an integral domain, otherwise if this cannot be guaranteed, the function returns false. For example, if MyElem was the type of elements of generic residue rings of a polynomial ring, the answer to the question would depend on the modulus of the residue  ring. Therefore isdomain_type would have to return false, since we cannot guarantee that we are dealing with elements of an integral domain in general. But if the given element type was for rational integers, the answer would be true, since every rational integer has as parent the ring of rational integers, which is an integral domain.Note that this function depends only on the type of an element and cannot access information about the object itself, or its parent.isexact_type(::Type{MyElem})Returns true if every element of the given type is represented exactly. For example, p-adic numbers, real and complex floating point numbers and power series are not exact, as we can only represent them in general with finite truncations. Similarly polynomials and matrices over inexact element types are themselves inexact.Integers, rationals, finite fields and polynomials and matrices over them are always exact.Note that MyElem may be parameterised or an abstract type, in which case every element of every type represented by MyElem must be exact, otherwise the function must return false.Base.hash(f::MyElem, h::UInt)Return a hash for the object f of type UInt. This is used as a hopefully cheap way to distinguish objects that differ arithmetically. If the object has components, e.g. the coefficients of a polynomial or elements of a matrix, these should be hashed recursively, passing the same parameter h to all levels. Each component should then be xor'd with h before combining the individual component hashes to give the final hash.The hash functions in AbstractAlgebra.jl usually start from some fixed 64 bit hexadecimal  value that has been picked at random by the library author for that type. That is then truncated to fit a UInt (in case the latter is not 64 bits). This ensures that objects that are the same arithmetically (or that have the same components), but have different types (or structures), are unlikely to hash to the same value.deepcopy_internal(f::MyElem, dict::ObjectIdDict)Return a copy of the given element, recursively copying all components of the object.Obviously the parent, if it is stored in the element, should not be copied. The new element should have precisely the same parent as the old object.For types that cannot self-reference themselves anywhere internally, the dict argument may be ignored.In the case that internal self-references are possible, please consult the Julia documentation on how to implement deepcopy_internal."
},

{
    "location": "rings.html#Constructors-1",
    "page": "Ring Interface",
    "title": "Constructors",
    "category": "section",
    "text": "Outer constructors for most AbstractAlgebra types are provided by overloading the call syntax for parent objects. If R is a parent object for a given ring we provide the following constructors.(R::MyParent)()Return the zero object of the given ring.(R::MyParent)(a::Integer)Coerce the given integer into the given ring.(R::MyParent)(a::MyElem)If a belongs to the given ring, the function returns it (without making a copy). Otherwise an error is thrown.For parameterised rings we also require a function to coerce from the base ring into the parent ring.(R::MyParent{T})(a::T) where T <: AbstractAlgebra.RingElemCoerce a into the ring R if a belongs to the base ring of R."
},

{
    "location": "rings.html#Basic-manipulation-of-rings-and-elements-1",
    "page": "Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "zero(R::MyParent)Return the zero element of the given ring.one(R::MyParent)Return the multiplicative identity of the given ring.iszero(f::MyElem)Return true if the given element is the zero element of the ring it belongs to.isone(f::MyElem)Return true if the given element is the multiplicative identity of the ring it belongs to."
},

{
    "location": "rings.html#Canonicalisation-1",
    "page": "Ring Interface",
    "title": "Canonicalisation",
    "category": "section",
    "text": "canonical_unit(f::MyElem)When fractions are created with two elements of the given type, it is nice to be able to represent them in some kind of canonical form. This is of course not always possible. But for example, fractions of integers can be canonicalised by first removing any common factors of the numerator and denominator, then making the denominator positive.In AbstractAlgebra.jl, the denominator would be made positive by dividing both the numerator and denominator by the canonical unit of the denominator. For a negative denominator, this would be -1.For elements of a field, canonical_unit simply returns the element itself. In general, canonical_unit of an invertible element should be that element. Finally, if a = ub we should have the identity canonical_unit(a) = canonical_unit(u)*canonical_unit(b).For some rings, it is completely impractical to implement this function, in which case it may return 1 in the given ring. The function must however always exist, and always return an element of the ring."
},

{
    "location": "rings.html#String-I/O-1",
    "page": "Ring Interface",
    "title": "String I/O",
    "category": "section",
    "text": "show(io::IO, R::MyParent)This should print (to the given IO object), an English description of the parent ring. If the ring is parameterised, it can call the corresponding show function for any rings it depends on.show(io::IO, f::MyElem)This should print a human readable, textual representation of the object (to the given IO object). It can recursively call the corresponding show functions for any of its components.It may be necessary in some cases to print parentheses around components of f or to print signs of components. For these, the following functions will exist for each component or component type.needs_parentheses(f::MyElem)Should returns true if parentheses are needed around this object when printed, e.g. as a coefficient of a polynomial. As an example, non-constant polynomials would need such parentheses if used as coefficients of another polynomial.isnegative(f::MyElem)When printing polynomials, a + sign is usually inserted automatically between terms of the polynomial. However, this is not desirable if the coefficient is negative and that negative sign is already printed when the coefficient is printed.This function must return true if a negative sign would already be prepended when f is printed. This suppresses the automatic printing of a + sign by polynomial printing functions that are printing f as a coefficient of a term.Note that if needs_parentheses returns true for f, then isnegative should always return false for that f, since an automatic + will need to be printed in front of a coefficient that is printed with parentheses.show_minus_one(::Type{MyElem})When printing polynomials, we prefer to print x rather than 1*x if the degree 1 term has coefficient 1. This can be taken care of without any special support.However, we also prefer to print -x rather than -1*x. This requires special support, since -1 in some rings is not printed as -1 (e.g. -1 in mathbbZ3mathbbZ might be printed as 2). In such rings, show_minus_one should return true.If show_minus_one returns true, polynomial printing functions will not print -x for terms of degree 1 with coefficient -1, but will use the printing function of the given type to print the coefficient in that case."
},

{
    "location": "rings.html#Unary-operations-1",
    "page": "Ring Interface",
    "title": "Unary operations",
    "category": "section",
    "text": "-(f::MyElem)Returns -f."
},

{
    "location": "rings.html#Binary-operations-1",
    "page": "Ring Interface",
    "title": "Binary operations",
    "category": "section",
    "text": "+(f::MyElem, g::MyElem)\n-(f::MyElem, g::MyElem)\n*(f::MyElem, g::MyElem)Returns f + g, f - g or fg, respectively."
},

{
    "location": "rings.html#Comparison-1",
    "page": "Ring Interface",
    "title": "Comparison",
    "category": "section",
    "text": "==(f::MyElem, g::MyElem)Returns true if f and g are arithmetically equal. In the case where the two elements are inexact, the function returns true if they agree to the minimum precision of the two.isequal(f::MyElem, g::MyElem)For exact rings, this should return the same thing as == above. For inexact rings, this returns true only if the two elements are arithmetically equal and have the same precision."
},

{
    "location": "rings.html#Powering-1",
    "page": "Ring Interface",
    "title": "Powering",
    "category": "section",
    "text": "^(f::MyElem, e::Int)Return f^e. The function should throw a DomainError() if negative exponents don't make sense but are passed to the function."
},

{
    "location": "rings.html#Exact-division-1",
    "page": "Ring Interface",
    "title": "Exact division",
    "category": "section",
    "text": "divexact(f::MyElem, g::MyElem)Returns fg, though note that Julia uses / for floating point division. Here we mean exact division in the ring, i.e. return q such that f = gq. A DivideError() should be thrown if g is zero. If no exact quotient exists or an impossible inverse is unavoidably encountered, an error should be thrown."
},

{
    "location": "rings.html#Unsafe-operators-1",
    "page": "Ring Interface",
    "title": "Unsafe operators",
    "category": "section",
    "text": "To speed up polynomial and matrix arithmetic, it sometimes makes sense to mutate values in place rather than replace them with a newly created object every time they are modified.For this purpose, certain mutating operators are required. In order to support immutable types (struct in Julia) and systems that don't have in-place operators, all unsafe operators must return the (ostensibly) mutated value. Only the returned value is used in computations, so this lifts the requirement that the unsafe operators actually mutate the value.Note the exclamation point is a convention, which indicates that the object may be mutated in-place.To make use of these functions, one must be certain that no other references are held to the object being mutated, otherwise those values will also be changed!The results of deepcopy and all arithmetic operations, including powering and division can be assumed to be new objects without other references being held, as can objects returned from constructors.Note that R(a) where R is the ring a belongs to, does not create a new value. For this case, use deepcopy(a).zero!(f::MyElem)Set the value f to zero in place. Return the mutated value.mul!(c::MyElem, a::MyElem, b::MyElem)Set c to the value ab in place. Return the mutated value. Aliasing is permitted.add!(c::MyElem, a::MyElem, b::MyElem)Set c to the value a + b in place. Return the mutated value. Aliasing is permitted.addeq!(a::MyElem, b::MyElem)Set a to a + b in place. Return the mutated value. Aliasing is permitted."
},

{
    "location": "rings.html#Random-generation-1",
    "page": "Ring Interface",
    "title": "Random generation",
    "category": "section",
    "text": "The random functions are only used for test code to generate test data. They therefore don't need to provide any guarantees on uniformity, and in fact, test values that are known to be a good source of corner cases can be supplied.rand(R::MyParent, v...)Returns a random element in the given ring of the specified size.There can be as many arguments as is necessary to specify the size of the test example which is being produced."
},

{
    "location": "rings.html#Promotion-rules-1",
    "page": "Ring Interface",
    "title": "Promotion rules",
    "category": "section",
    "text": "In order for AbstractAlgebra to be able to automatically coerce up towers of rings, certain promotion rules must be defined. For every ring, one wants to be able to coerce integers into the ring. And for any ring constructed over a base ring, one would like to be able to coerce from the base ring into the ring.The promotion rules look a bit different depending on whether the element type is parameterised or not and whether it is built on a base ring.For ring element types MyElem that are neither parameterised, not built over a base ring, the promotion rules can be defined as follows:promote_rule(::Type{MyElem}, ::Type{T}) where {T <: Integer} = MyElemFor ring element types MyType that aren't parameterised, but which have a base ring with concrete element type T the promotion rules can be defined as follows:promote_rule(::Type{MyElem}, ::Type{U}) where U <: Integer = MyElempromote_rule(::Type{MyElem}, ::Type{T}) = MyElemFor ring element types MyElem{T} that are parameterised by the type of elements of the base ring, the promotion rules can be defined as follows:promote_rule(::Type{MyElem{T}}, ::Type{MyElem{T}}) where T <: RingElement = MyElem{T}function promote_rule(::Type{MyElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingEle\nment}\n   promote_rule(T, U) == T ? MyElem{T} : Union{}\nend"
},

{
    "location": "rings.html#Required-functionality-for-inexact-rings-1",
    "page": "Ring Interface",
    "title": "Required functionality for inexact rings",
    "category": "section",
    "text": ""
},

{
    "location": "rings.html#Approximation-(floating-point-and-ball-arithmetic-only)-1",
    "page": "Ring Interface",
    "title": "Approximation (floating point and ball arithmetic only)",
    "category": "section",
    "text": "isapprox(f::MyElem, g::MyElem; atol::Real=sqrt(eps()))This is used by test code that uses rings involving floating point or ball arithmetic. The function should return true if all components of f and g are equal to within the square root of the Julia epsilon, since numerical noise may make an exact comparison impossible.For parameterised rings over an inexact ring, we also require the following ad hoc approximation functionality.isapprox(f::MyElem{T}, g::T; atol::Real=sqrt(eps())) where T <: AbstractAlgebra.RingElemisapprox(f::T, g::MyElem{T}; atol::Real=sqrt(eps())) where T <: AbstractAlgebra.RingElemThese notionally coerce the element of the base ring into the parameterised ring and do a full comparison."
},

{
    "location": "rings.html#Optional-functionality-1",
    "page": "Ring Interface",
    "title": "Optional functionality",
    "category": "section",
    "text": "Some functionality is difficult or impossible to implement for all rings in the system. If it is provided, additional functionality or performance may become available. Here is a list of all functions that are considered optional and can't be relied on by generic functions in the AbstractAlgebra Ring interface.It may be that no algorithm, or no efficient algorithm is known to implement these functions. As these functions are optional, they do not need to exist. Julia will already inform the user that the function has not been implemented if it is called but doesn't exist."
},

{
    "location": "rings.html#Optional-basic-manipulation-functionality-1",
    "page": "Ring Interface",
    "title": "Optional basic manipulation functionality",
    "category": "section",
    "text": "isunit(f::MyElem)Return true if the given element is a unit in the ring it belongs to. "
},

{
    "location": "rings.html#Optional-binary-ad-hoc-operators-1",
    "page": "Ring Interface",
    "title": "Optional binary ad hoc operators",
    "category": "section",
    "text": "By default, ad hoc operations are handled by AbstractAlgebra.jl if they are not defined explicitly, by coercing both operands into the same ring and then performing the required operation.In some cases, e.g. for matrices, this leads to very inefficient behaviour. In such cases, it is advised to implement some of these operators explicitly.It can occasionally be worth adding a separate set of ad hoc binary operators for the type Int, if this can be done more efficiently than for arbitrary Julia Integer types.+(f::MyElem, c::Integer)\n-(f::MyElem, c::Integer)\n*(f::MyElem, c::Integer)+(c::Integer, f::MyElem)\n-(c::Integer, f::MyElem)\n*(c::Integer, f::MyElem)For parameterised types, it is also sometimes more performant to provide explicit ad hoc operators with elements of the base ring.+(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem\n-(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem\n*(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem+(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem\n-(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem\n*(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem"
},

{
    "location": "rings.html#Optional-ad-hoc-comparisons-1",
    "page": "Ring Interface",
    "title": "Optional ad hoc comparisons",
    "category": "section",
    "text": "==(f::MyElem, c::Integer)==(c::Integer, f::MyElem)==(f::MyElem{T}, c:T) where T <: AbstractAlgebra.RingElem==(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem"
},

{
    "location": "rings.html#Optional-ad-hoc-exact-division-functions-1",
    "page": "Ring Interface",
    "title": "Optional ad hoc exact division functions",
    "category": "section",
    "text": "divexact(a::MyType{T}, b::T) where T <: AbstractAlgebra.RingElemdivexact(a::MyType, b::Integer)"
},

{
    "location": "rings.html#Optional-powering-functions-1",
    "page": "Ring Interface",
    "title": "Optional powering functions",
    "category": "section",
    "text": "^(f::MyElem, e::BigInt)In case f cannot explode in size when powered by a very large integer, and it is practical to do so, one may provide this function to support powering with BigInt exponents (or for external modules, any other big integer type)."
},

{
    "location": "rings.html#Optional-unsafe-operators-1",
    "page": "Ring Interface",
    "title": "Optional unsafe operators",
    "category": "section",
    "text": "addmul!(c::MyElem, a::MyElem, b::MyElem, t::MyElem)Set c = c + ab in-place. Return the mutated value. The value t should be a temporary of the same type as a, b and c, which can be used arbitrarily by the implementation to speed up the computation. Aliasing between a, b and c is  permitted."
},

{
    "location": "euclidean.html#",
    "page": "Euclidean Ring Interface",
    "title": "Euclidean Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "euclidean.html#Euclidean-Ring-Interface-1",
    "page": "Euclidean Ring Interface",
    "title": "Euclidean Ring Interface",
    "category": "section",
    "text": "If a ring provides a meaningful Euclidean structure such that a useful Euclidean remainder can be computed practically, various additional functionality is provided by AbstractAlgebra.jl for those rings. This functionality depends on the following functions existing.mod(f::MyElem, g::MyElem)Returns the Euclidean remainder of f by g. A DivideError() should be thrown if g is zero. An error should be thrown if an impossible inverse is encountered.divrem(f::MyElem, g::MyElem)Returns a pair q, r consisting of the Euclidean quotient and remainder of f by g. A DivideError should be thrown if g is zero. An error should be thrown if an impossible inverse is encountered.div(f::MyElem, g::MyElem)Returns the Euclidean quotient of f by g. A DivideError should be thrown if g is zero. An error should be thrown if an impossible inverse is encountered.mulmod(f::MyElem, g::MyElem, m::MyElem)Returns fg pmodm.powmod(f::MyElem, e::Int, m::MyElem)Returns f^e pmodm.invmod(f::MyElem, m::MyElem)Returns the inverse of f modulo m. If such an inverse doesn't exist, an impossible inverse error should be thrown.divides(f::MyElem, g::MyElem)Returns a pair, flag, q, where flag is set to true if g divides f, in which case the quotient is set to the quotient, or flag is set to false and the quotient is set to zero in the same ring as f and g.remove(f::MyElem, p::MyElem)Returns a pair v, q where p^v is the highest power of p dividing f, and q is the cofactor after f is divided by this power.valuation(f::MyElem, p::MyElem)Returns v where p^v is the highest power of p dividing f.gcd(f::MyElem, g::MyElem)Returns a greatest common divisor of f and g.lcm(f::MyElem, g::MyElem)Returns fggcd(f g) if either f or g is not zero, otherwise it throws a DivideError().gcdx(f::MyElem, g::MyElem)Returns a triple d, s, t such that d = gcd(f g) and d = sf + tg, with s reduced modulo g and t reduced modulo f.gcdinv(f::MyElem, g::MyElem)Returns a tuple d, s such that d = gcd(f g) and s = (fd)^-1 pmodgd. Note that d = 1 iff f is invertible modulo g, in which case s = f^-1 pmodg."
},

{
    "location": "integer.html#",
    "page": "Integer ring",
    "title": "Integer ring",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "integer.html#Integer-ring-1",
    "page": "Integer ring",
    "title": "Integer ring",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/julia/Integer.jl for making Julia BigInts conform to the AbstractAlgebra.jl Ring interface.In addition to providing a parent object ZZ for Julia BigInts, we implement any additional functionality required by AbstractAlgebra.jl.Because BigInt cannot be directly included in the AbstractAlgebra.jl abstract type hierarchy, we achieve integration of Julia BigInts by introducing a type union, called RingElement, which is a union of AbstractAlgebra.RingElem and a number of Julia types, including BigInt. Everywhere that RingElem is notionally used in AbstractAlgebra.jl, we are in fact using RingElement, with additional care being taken to avoid ambiguities.The details of how this is done are technical, and we refer the reader to the implementation for details. For most intents and purposes, one can think of the Julia BigInt type as belonging to AbstractAlgebra.RingElem.One other technicality is that Julia defines certain functions for BigInt, such as sqrt and exp differently to what AbstractAlgebra.jl requires. To get around this, we redefine these functions internally to AbstractAlgebra.jl, without redefining them for users of AbstractAlgebra.jl. This allows the internals of AbstractAlgebra.jl to function correctly, without broadcasting pirate definitions of already defined Julia functions to the world.To access the internal definitions, one can use AbstractAlgebra.sqrt and AbstractAlgebra.exp, etc."
},

{
    "location": "integer.html#Types-and-parent-objects-1",
    "page": "Integer ring",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Integers have type BigInt, as in Julia itself. We simply supplement the functionality for this type as required for computer algebra.The parent objects of such integers has type Integers{BigInt}.For convenience, we also make Int a part of the AbstractAlgebra.jl type hierarchy and its parent object (accessible as zz) has type Integers{Int}. But we caution that this type is not particularly useful as a model of the integers and may not function as expected within AbstractAlgebra.jl."
},

{
    "location": "integer.html#Integer-constructors-1",
    "page": "Integer ring",
    "title": "Integer constructors",
    "category": "section",
    "text": "In order to construct integers in AbstractAlgebra.jl, one can first construct the integer ring itself. This is accomplished using the following constructor.Integers{BigInt}()This gives the unique object of type Integers{BigInt} representing the ring of integers in AbstractAlgebra.jl.In practice, one simply uses ZZ which is assigned to be the return value of the above constructor. There is no need to call the constructor in practice.Here are some examples of creating the integer ring and making use of the resulting parent object to coerce various elements into the ring.Examplesf = ZZ()\ng = ZZ(123)\nh = ZZ(BigInt(1234))"
},

{
    "location": "integer.html#Basic-ring-functionality-1",
    "page": "Integer ring",
    "title": "Basic ring functionality",
    "category": "section",
    "text": "The integer ring in AbstractAlgebra.jl implements the full Ring interface and the  Euclidean Ring interface.We give some examples of such functionality.Examplesf = ZZ(12)\n\nh = zero(ZZ)\nk = one(ZZ)\nisone(k) == true\niszero(f) == false\nU = base_ring(ZZ)\nV = base_ring(f)\nT = parent(f)\nf == deepcopy(f)\ng = f + 12\nh = powmod(f, 12, ZZ(17))\nflag, q = divides(f, ZZ(3))"
},

{
    "location": "integer.html#Integer-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Integer ring",
    "title": "Integer functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality below supplements that provided by Julia itself for its BigInt type."
},

{
    "location": "integer.html#AbstractAlgebra.Generic.isunit-Tuple{Integer}",
    "page": "Integer ring",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "Method",
    "text": "isunit(a::Integer)\n\nReturn true if a is 1 or -1.\n\n\n\n"
},

{
    "location": "integer.html#Basic-functionality-1",
    "page": "Integer ring",
    "title": "Basic functionality",
    "category": "section",
    "text": "isunit(::Integer)Examplesr = ZZ(-1)\n\nisunit(r) == true"
},

{
    "location": "integer.html#AbstractAlgebra.sqrt-Tuple{BigInt}",
    "page": "Integer ring",
    "title": "AbstractAlgebra.sqrt",
    "category": "Method",
    "text": "sqrt{T <: Integer}(a::T)\n\nReturn the integer square root of a. If a is not a perfect square an error is thrown.\n\n\n\nsqrt{T <: Integer}(a::Rational{T})\n\nReturn the square root of a if it is the square of a rational, otherwise throw an error.\n\n\n\n"
},

{
    "location": "integer.html#AbstractAlgebra.exp-Tuple{BigInt}",
    "page": "Integer ring",
    "title": "AbstractAlgebra.exp",
    "category": "Method",
    "text": "exp{a <: Integer}(a::T)\n\nReturn 1 if a = 0, otherwise throw an exception. This function is not generally of use to the user, but is used internally in AbstractAlgebra.jl.\n\n\n\nexp{T <: Integer}(a::Rational{T})\n\nReturn 1 if a = 0, otherwise throw an exception.\n\n\n\n"
},

{
    "location": "integer.html#Square-root-1",
    "page": "Integer ring",
    "title": "Square root",
    "category": "section",
    "text": "AbstractAlgebra.sqrt(a::BigInt)AbstractAlgebra.exp(a::BigInt)Examplesd = AbstractAlgebra.sqrt(ZZ(36))\nm = AbstractAlgebra.exp(ZZ(0))"
},

{
    "location": "integer.html#AbstractAlgebra.ppio-Tuple{BigInt,BigInt}",
    "page": "Integer ring",
    "title": "AbstractAlgebra.ppio",
    "category": "Method",
    "text": "ppio(a::T, b::T)\n\nSplit a into c*d where c = gcd(a b^infty).\n\n\n\n"
},

{
    "location": "integer.html#Coprime-bases-1",
    "page": "Integer ring",
    "title": "Coprime bases",
    "category": "section",
    "text": "ppio(a::BigInt, b::BigInt)Examplesc, n = ppio(ZZ(12), ZZ(26))"
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
    "text": "In addition to the standard constructors, the following constructors, taking an array of coefficients, must be available.(S::MyPolyRing{T})(A::Array{T, 1}) where T <: AbstractAlgebra.RingElemCreate the polynomial in the given ring whose degree i coefficient is given by A[i].(S::MyPolyRing{T})(A::Array{U, 1}) where T <: AbstractAlgebra.RingElem, U <: AbstractAlgebra.RingElemCreate the polynomial in the given ring whose degree i coefficient is given by A[i]. The elements of the array are assumed to be able to be coerced into the base ring R.(S::MyPolyRing{T})(A::Array{U, 1}) where T <: AbstractAlgebra.RingElem, U <: IntegerCreate the polynomial in the given ring whose degree i coefficient is given by A[i].It may be desirable to have a additional version of the function that accepts an array of Julia Int values  if this can be done more efficiently.ExamplesS, x = PolynomialRing(QQ, \"x\")\n\nf = S(Rational{BigInt}[2, 3, 1])\ng = S(BigInt[1, 0, 4])\nh = S([4, 7, 2, 9])"
},

{
    "location": "polynomial_rings.html#Data-type-and-parent-object-methods-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "var(S::MyPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn a Symbol representing the variable (generator) of the polynomial ring. Note that this is a Symbol not a String, though its string value will usually be used when printing polynomials.vars(S::MyPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn the array [s] where s	 is aSymbol` representing the variable of the given polynomial ring. This is provided for uniformity with the multivariate interface, where there is more than one variable, and hence an array of symbols.ExamplesS, x = PolynomialRing(QQ, \"x\")\n\nvsym = var(S)\nV = vars(S)"
},

{
    "location": "polynomial_rings.html#Basic-manipulation-of-rings-and-elements-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "length(f::MyPoly{T}) where T <: AbstractAlgebra.RingElemReturn the length of the given polynomial. The length of the zero polynomial is defined to be 0, otherwise the length is the degree plus 1. The return value should be of type Int.set_length!(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemThis function must zero any coefficients beyond the requested length n and then set the length of the polynomial to n. This function does not need to normalise the polynomial and is not useful to the user, but is used extensively by the AbstractAlgebra generic functionality.This function mutates the existing polynomial in-place, but does not return the polynomial.coeff(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemReturn the coefficient of the polynomial f of degree n. If n is larger than the degree of the polynomial, it should return zero in the coefficient ring. setcoeff!(f::MyPoly{T}, n::Int, a::T) where T <: AbstractAlgebra.RingElemSet the degree n coefficient of f to a. This mutates the polynomial in-place if possible and returns the mutated polynomial (so that immutable types can also be supported). The function must not assume that the polynomial already has space for n + 1 coefficients. The polynomial must be resized if this is not the case.Note that this function is not required to normalise the polynomial and is not necessarily useful to the user, but is used extensively by the generic functionality in AbstractAlgebra.jl. It is for setting raw coefficients in the representation.normalise(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemGiven a polynomial whose length is currently n, including any leading zero coefficients, return the length of the normalised polynomial (either zero of the length of the polynomial with nonzero leading coefficient). Note that the function does not actually perform the normalisation.fit!(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemEnsure that the polynomial f internally has space for n coefficients. This function must mutate the function in-place if it is mutable. It does not return the mutated polynomial. Immutable types can still be supported by defining this function to do nothing.Some interfaces for C polynomial types automatically manage the internal allocation of polynomials in every function that can be called on them. Explicit adjustment by the generic code in AbstractAlgebra.jl is not required. In such cases, this function can also be defined to do nothing.ExamplesS, x = PolynomialRing(ZZ, \"x\")\n\nf = x^3 + 3x + 1\ng = S(BigInt[1, 2, 0, 1, 0, 0, 0]);\n\nn = length(f)\nc = coeff(f, 1)\nset_length!(g, normalise(g, 7))\ng = setcoeff!(g, 2, BigInt(11))\nfit!(g, 8)\ng = setcoeff!(g, 7, BigInt(4))\n"
},

{
    "location": "polynomial_rings.html#Optional-functionality-for-polynomial-rings-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Optional functionality for polynomial rings",
    "category": "section",
    "text": "Sometimes parts of the Euclidean Ring interface can and should be implemented for polynomials over a ring that is not necessarily a field.When divisibility testing can be implemented for a polynomial ring over a field, it  should be possible to implement the following functions from the Euclidean Ring interface:divides\nremove\nvaluationWhen the given polynomial ring is a GCD domain, with an effective GCD algorithm, it may be possible to implement the following functions:gcd\nlcmPolynomial rings can optionally implement any part of the generic univariate polynomial functionality provided by AbstractAlgebra.jl, using the same interface. Obviously additional functionality can also be added to that provided by AbstractAlgebra.jl on an ad hoc basis."
},

{
    "location": "polynomial.html#",
    "page": "Generic univariate polynomials",
    "title": "Generic univariate polynomials",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "polynomial.html#Generic-univariate-polynomials-1",
    "page": "Generic univariate polynomials",
    "title": "Generic univariate polynomials",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/generic/Poly.jl for generic polynomials over any commutative ring belonging to the AbstractAlgebra abstract type hierarchy.As well as implementing the Univariate Polynomial interface, and relevant parts of the Euclidean Ring interface for polynomials over a field, there are many additional generic algorithms implemented for such polynomial rings. We describe this generic functionality below.All of the generic functionality is part of a submodule of AbstractAlgebra called Generic. This is exported by default so that it is not necessary to qualify the function names with the submodule name."
},

{
    "location": "polynomial.html#Types-and-parent-objects-1",
    "page": "Generic univariate polynomials",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Polynomials implemented using the AbstractAlgebra generics have type Generic.Poly{T} where T is the type of elements of the coefficient ring. Internally they consist of a Julia array of coefficients and some additional fields for length and a parent object, etc. See the file src/generic/GenericTypes.jl for details.Parent objects of such polynomials have type Generic.PolyRing{T}.The string representation of the variable of the polynomial ring, and the base/coefficient ring R is stored in the parent object. The polynomial element types belong to the abstract type AbstractAlgebra.PolyElem{T} and the polynomial ring types belong to the abstract type AbstractAlgebra.PolyRing{T}. This enables one to write generic functions that can accept any AbstractAlgebra polynomial type.Note that both the generic polynomial ring type Generic.PolyRing{T} and the abstract type it belongs to, AbstractAlgebra.PolyRing{T} are both called PolyRing. The  former is a (parameterised) concrete type for a polynomial ring over a given base ring whose elements have type T. The latter is an abstract type representing all polynomial ring types in AbstractAlgebra.jl, whether generic or very specialised (e.g. supplied by a C library)."
},

{
    "location": "polynomial.html#Polynomial-ring-constructors-1",
    "page": "Generic univariate polynomials",
    "title": "Polynomial ring constructors",
    "category": "section",
    "text": "In order to construct polynomials in AbstractAlgebra.jl, one must first construct the polynomial ring itself. This is accomplished with the following constructor.PolynomialRing(R::AbstractAlgebra.Ring, s::AbstractString; cached::Bool = true)Given a base ring R and string s specifying how the generator (variable) should be printed, return a tuple S, x representing the new polynomial ring S = Rx and the generator x of the ring. By default the parent object S will depend only on R and  x and will be cached. Setting the optional argument cached to false will prevent the parent object S from being cached.A shorthand version of this function is provided: given a base ring R, we abbreviate the constructor as follows.R[\"x\"]Here are some examples of creating polynomial rings and making use of the resulting parent objects to coerce various elements into the polynomial ring.ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\nT, z = QQ[\"z\"]\n\nf = R()\ng = R(123)\nh = S(BigInt(1234))\nk = S(x + 1)\nm = T(z + 1)All of the examples here are generic polynomial rings, but specialised implementations of polynomial rings provided by external modules will also usually provide a PolynomialRing constructor to allow creation of their polynomial rings."
},

{
    "location": "polynomial.html#Basic-ring-functionality-1",
    "page": "Generic univariate polynomials",
    "title": "Basic ring functionality",
    "category": "section",
    "text": "Once a polynomial ring is constructed, there are various ways to construct polynomials in that ring.The easiest way is simply using the generator returned by the PolynomialRing constructor and build up the polynomial using basic arithmetic, as described in the Ring interface. The Julia language also has special syntax for the construction of polynomials in terms of a generator, e.g. we can write 2x instead of 2*x.The polynomial rings in AbstractAlgebra.jl implement the full Ring interface. Of course the entire Univariate Polynomial Ring interface is also implemented.We give some examples of such functionality.ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = x^3 + 3x + 21\ng = (x + 1)*y^2 + 2x + 1\n\nh = zero(S)\nk = one(R)\nisone(k) == true\niszero(f) == false\nn = length(g)\nU = base_ring(S)\nV = base_ring(y + 1)\nv = var(S)\nT = parent(y + 1)\ng == deepcopy(g)\nt = divexact(2g, 2)For polynomials over a field, the Euclidean Ring interface is implemented.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\nT, y = PolynomialRing(S, \"y\")\n\nf = (3*x^2 + x + 2)*y + x^2 + 1\ng = (5*x^2 + 2*x + 1)*y^2 + 2x*y + x + 1\nh = (3*x^3 + 2*x^2 + x + 7)*y^5 + 2x*y + 1\n\ninvmod(f, g)\nmulmod(f, g, h)\npowmod(f, 3, h)\nh = mod(f, g)\nq, r = divrem(f, g)\nd = gcd(f*h, g*h)\nk = gcdinv(f, h)\nm = lcm(f, h)\nflag, q = divides(g^2, g)\nvaluation(3g^3, g) == 3\nval, q = remove(5g^3, g)\nr, s, t = gcdx(g, h)Functions in the Euclidean Ring interface are supported over residue rings that are not fields, except that if an impossible inverse is encountered during the computation an error is thrown."
},

{
    "location": "polynomial.html#Polynomial-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic univariate polynomials",
    "title": "Polynomial functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality listed below is automatically provided by AbstractAlgebra.jl for any polynomial module that implements the full Univariate Polynomial Ring interface. This includes AbstractAlgebra.jl's own generic polynomial rings.But if a C library provides all the functionality documented in the Univariate Polynomial Ring interface, then all the functions described here will also be  automatically supplied by AbstractAlgebra.jl for that polynomial type.Of course, modules are free to provide specific implementations of the functions described here, that override the generic implementation."
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.modulus-Union{Tuple{AbstractAlgebra.PolyElem{T}}, Tuple{T}} where T<:AbstractAlgebra.ResElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.modulus",
    "category": "Method",
    "text": "modulus{T <: ResElem}(a::AbstractAlgebra.PolyElem{T})\n\nReturn the modulus of the coefficients of the given polynomial.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.lead-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.lead",
    "category": "Method",
    "text": "lead(x::AbstractAlgebra.PolyElem)\n\nReturn the leading coefficient of the given polynomial. This will be the nonzero coefficient of the term with highest degree unless the polynomial in the zero polynomial, in which case a zero coefficient is returned.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.trail-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.trail",
    "category": "Method",
    "text": "trail(x::AbstractAlgebra.PolyElem)\n\nReturn the trailing coefficient of the given polynomial. This will be the nonzero coefficient of the term with lowest degree unless the polynomial in the zero polynomial, in which case a zero coefficient is returned.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.gen-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.gen",
    "category": "Method",
    "text": "gen{T <: RingElement}(R::AbsSeriesRing{T})\n\nReturn the generator of the power series ring, i.e. x + O(x^n) where n is the precision of the power series ring R.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.isgen-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.isgen",
    "category": "Method",
    "text": "isgen(a::AbstractAlgebra.PolyElem)\n\nReturn true if the given polynomial is the constant generator of its polynomial ring, otherwise return false.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.isunit-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "Method",
    "text": "isunit(a::AbstractAlgebra.PolyElem)\n\nReturn true if the given polynomial is a unit in its polynomial ring, otherwise return false.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.ismonomial-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.ismonomial",
    "category": "Method",
    "text": "ismonomial(a::AbstractAlgebra.PolyElem)\n\nReturn true if the given polynomial is a monomial.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.isterm-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.isterm",
    "category": "Method",
    "text": "isterm(a::AbstractAlgebra.PolyElem)\n\nReturn true if the given polynomial is has one term. This function is recursive, with all scalar types returning true.\n\n\n\n"
},

{
    "location": "polynomial.html#Basic-functionality-1",
    "page": "Generic univariate polynomials",
    "title": "Basic functionality",
    "category": "section",
    "text": "modulus{T <: ResElem}(::PolyElem{T})lead(::PolyElem)\ntrail(::PolyElem)gen(::PolyElem)isgen(::PolyElem)isunit(::PolyElem)ismonomial(::PolyElem)isterm(::PolyElem)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\nT, z = PolynomialRing(QQ, \"z\")\nU = ResidueRing(ZZ, 17)\nV, w = PolynomialRing(U, \"w\")\n\na = zero(S)\nb = one(S)\n\nc = BigInt(1)//2*z^2 + BigInt(1)//3\nd = x*y^2 + (x + 1)*y + 3\n\nf = lead(d)\ny = gen(S)\ng = isgen(w)\nm = isunit(b)\nn = degree(d)\nr = modulus(w)\nisterm(2y^2) == true\nismonomial(y^2) == true"
},

{
    "location": "polynomial.html#Base.truncate-Tuple{AbstractAlgebra.PolyElem,Int64}",
    "page": "Generic univariate polynomials",
    "title": "Base.truncate",
    "category": "Method",
    "text": "truncate(a::AbstractAlgebra.PolyElem, n::Int)\n\nReturn a truncated to n terms.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.mullow-Union{Tuple{AbstractAlgebra.PolyElem{T},AbstractAlgebra.PolyElem{T},Int64}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.mullow",
    "category": "Method",
    "text": "mullow{T <: RingElement}(a::AbstractAlgebra.PolyElem{T}, b::AbstractAlgebra.PolyElem{T}, n::Int)\n\nReturn atimes b truncated to n terms.\n\n\n\n"
},

{
    "location": "polynomial.html#Truncation-1",
    "page": "Generic univariate polynomials",
    "title": "Truncation",
    "category": "section",
    "text": "truncate(::PolyElem, ::Int)mullow{T <: RingElem}(::PolyElem{T}, ::PolyElem{T}, ::Int)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\ng = (x + 1)*y + (x^3 + 2x + 2)\n\nh = truncate(f, 1)\nk = mullow(f, g, 4)"
},

{
    "location": "polynomial.html#Base.reverse-Tuple{AbstractAlgebra.PolyElem,Int64}",
    "page": "Generic univariate polynomials",
    "title": "Base.reverse",
    "category": "Method",
    "text": "reverse(x::AbstractAlgebra.PolyElem, len::Int)\n\nReturn the reverse of the polynomial x, thought of as a polynomial of the given length (the polynomial will be notionally truncated or padded with zeroes before the leading term if necessary to match the specified length). The resulting polynomial is normalised. If len is negative we throw a DomainError().\n\n\n\n"
},

{
    "location": "polynomial.html#Base.reverse-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "Base.reverse",
    "category": "Method",
    "text": "reverse(x::AbstractAlgebra.PolyElem)\n\nReturn the reverse of the polynomial x, i.e. the leading coefficient of x becomes the constant coefficient of the result, etc. The resulting polynomial is normalised.\n\n\n\n"
},

{
    "location": "polynomial.html#Reversal-1",
    "page": "Generic univariate polynomials",
    "title": "Reversal",
    "category": "section",
    "text": "reverse(::PolyElem, ::Int)\nreverse(::PolyElem)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\n\ng = reverse(f, 7)\nh = reverse(f)"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.shift_left-Tuple{AbstractAlgebra.PolyElem,Int64}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.shift_left",
    "category": "Method",
    "text": "shift_left(x::AbstractAlgebra.PolyElem, n::Int)\n\nReturn the polynomial f shifted left by n terms, i.e. multiplied by x^n.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.shift_right-Tuple{AbstractAlgebra.PolyElem,Int64}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.shift_right",
    "category": "Method",
    "text": "shift_right(f::AbstractAlgebra.PolyElem, n::Int)\n\nReturn the polynomial f shifted right by n terms, i.e. divided by x^n.\n\n\n\n"
},

{
    "location": "polynomial.html#Shifting-1",
    "page": "Generic univariate polynomials",
    "title": "Shifting",
    "category": "section",
    "text": "shift_left(::PolyElem, ::Int)shift_right(::PolyElem, ::Int)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\n\ng = shift_left(f, 7)\nh = shift_right(f, 2)"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.pseudorem-Union{Tuple{AbstractAlgebra.PolyElem{T},AbstractAlgebra.PolyElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.pseudorem",
    "category": "Method",
    "text": "pseudorem{T <: RingElement}(f::AbstractAlgebra.PolyElem{T}, g::AbstractAlgebra.PolyElem{T})\n\nReturn the pseudoremainder of a divided by b. If b = 0 we throw a DivideError().\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.pseudodivrem-Union{Tuple{AbstractAlgebra.PolyElem{T},AbstractAlgebra.PolyElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.pseudodivrem",
    "category": "Method",
    "text": "pseudodivrem{T <: RingElement}(f::AbstractAlgebra.PolyElem{T}, g::AbstractAlgebra.PolyElem{T})\n\nReturn a tuple (q r) consisting of the pseudoquotient and pseudoremainder of a divided by b. If b = 0 we throw a DivideError().\n\n\n\n"
},

{
    "location": "polynomial.html#Pseudodivision-1",
    "page": "Generic univariate polynomials",
    "title": "Pseudodivision",
    "category": "section",
    "text": "Given two polynomials a b, pseudodivision computes polynomials q and r with length(r)  length(b) such that L^d a = bq + r where d = length(a) - length(b) + 1 and L is the leading coefficient of b.We call q the pseudoquotient and r the pseudoremainder.pseudorem{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})pseudodivrem{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\ng = (x + 1)*y + (x^3 + 2x + 2)\n\nh = pseudorem(f, g)\nq, r = pseudodivrem(f, g)"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.content-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.content",
    "category": "Method",
    "text": "content(a::AbstractAlgebra.PolyElem)\n\nReturn the content of a, i.e. the greatest common divisor of its coefficients.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.primpart-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.primpart",
    "category": "Method",
    "text": "primpart(a::AbstractAlgebra.PolyElem)\n\nReturn the primitive part of a, i.e. the polynomial divided by its content.\n\n\n\n"
},

{
    "location": "polynomial.html#Content-and-primitive-part-1",
    "page": "Generic univariate polynomials",
    "title": "Content and primitive part",
    "category": "section",
    "text": "content(::PolyElem)primpart(::PolyElem)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nk = x*y^2 + (x + 1)*y + 3\n\nn = content(k)\np = primpart(k*(x^2 + 1))"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.evaluate-Union{Tuple{AbstractAlgebra.PolyElem{T},T}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "Method",
    "text": "evaluate{T <: RingElement}(a::AbstractAlgebra.PolyElem{T}, b::T)\n\nEvaluate the polynomial a at the value b and return the result.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.evaluate-Tuple{AbstractAlgebra.PolyElem,Integer}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "Method",
    "text": "evaluate{T <: RingElement}(a::AbstractAlgebra.PolyElem{T}, b::T)\n\nEvaluate the polynomial a at the value b and return the result.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.compose-Tuple{AbstractAlgebra.PolyElem,AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.compose",
    "category": "Method",
    "text": "compose(a::AbstractAlgebra.PolyElem, b::AbstractAlgebra.PolyElem)\n\nCompose the polynomial a with the polynomial b and return the result, i.e. return acirc b.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.subst-Union{Tuple{AbstractAlgebra.PolyElem{T},Any}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.subst",
    "category": "Method",
    "text": "subst{T <: RingElement}(f::AbstractAlgebra.PolyElem{T}, a::Any)\n\nEvaluate the polynomial f at a. Note that a can be anything, whether a ring element or not.\n\n\n\n"
},

{
    "location": "polynomial.html#Evaluation,-composition-and-substitution-1",
    "page": "Generic univariate polynomials",
    "title": "Evaluation, composition and substitution",
    "category": "section",
    "text": "evaluate{T <: RingElem}(::PolyElem{T}, ::T)\nevaluate(::PolyElem, ::Integer)compose(::PolyElem, ::PolyElem)subst{T <: RingElem}(::PolyElem{T}, ::Any)We also overload the functional notation so that the polynomial f can be evaluated at a by writing f(a). ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n   \nf = x*y^2 + (x + 1)*y + 3\ng = (x + 1)*y + (x^3 + 2x + 2)\nM = R[x + 1 2x; x - 3 2x - 1]\n\nk = evaluate(f, 3)\nm = evaluate(f, x^2 + 2x + 1)\nn = compose(f, g)\np = subst(f, M)\nq = f(M)\nr = f(23)"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.derivative-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.derivative",
    "category": "Method",
    "text": "derivative(a::AbstractAlgebra.PolyElem)\n\nReturn the derivative of the polynomial a.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.integral-Union{Tuple{AbstractAlgebra.PolyElem{T}}, Tuple{T}} where T<:Union{AbstractAlgebra.FieldElem, AbstractAlgebra.ResElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.integral",
    "category": "Method",
    "text": "integral{T <: Union{AbstractAlgebra.ResElem, FieldElement}}(x::AbstractAlgebra.PolyElem{T})\n\nReturn the integral of the polynomial a.\n\n\n\n"
},

{
    "location": "polynomial.html#Derivative-and-integral-1",
    "page": "Generic univariate polynomials",
    "title": "Derivative and integral",
    "category": "section",
    "text": "derivative(::PolyElem)integral{T <: Union{ResElem, FieldElem}}(::PolyElem{T})ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\nT, z = PolynomialRing(QQ, \"z\")\nU = ResidueRing(T, z^3 + 3z + 1)\nV, w = PolynomialRing(U, \"w\")\n\nf = x*y^2 + (x + 1)*y + 3\ng = (z^2 + 2z + 1)*w^2 + (z + 1)*w - 2z + 4\n\nh = derivative(f)\nk = integral(g)   "
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.resultant-Union{Tuple{AbstractAlgebra.PolyElem{T},AbstractAlgebra.PolyElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.resultant",
    "category": "Method",
    "text": "resultant{T <: RingElem}(a::AbstractAlgebra.PolyElem{T}, b::AbstractAlgebra.PolyElem{T})\n\nReturn the resultant of the given polynomials.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.resx-Union{Tuple{AbstractAlgebra.PolyElem{T},AbstractAlgebra.PolyElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.resx",
    "category": "Method",
    "text": "resx{T <: RingElement}(a::AbstractAlgebra.PolyElem{T}, b::AbstractAlgebra.PolyElem{T})\n\nReturn a tuple (r s t) such that r is the resultant of a and b and such that r = atimes s + btimes t.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.discriminant-Tuple{AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.discriminant",
    "category": "Method",
    "text": "discriminant(a::AbstractAlgebra.PolyElem)\n\nReturn the discrimnant of the given polynomial.\n\n\n\n"
},

{
    "location": "polynomial.html#Resultant-and-discriminant-1",
    "page": "Generic univariate polynomials",
    "title": "Resultant and discriminant",
    "category": "section",
    "text": "resultant{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})resx{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})discriminant(a::PolyElem)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = 3x*y^2 + (x + 1)*y + 3\ng = 6(x + 1)*y + (x^3 + 2x + 2)\n\nh = resultant(f, g)\nk = discriminant(f)"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.monomial_to_newton!-Union{Tuple{Array{T,1},Array{T,1}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.monomial_to_newton!",
    "category": "Method",
    "text": "monomial_to_newton!{T <: RingElement}(P::Array{T, 1}, roots::Array{T, 1})\n\nConverts a polynomial p, given as an array of coefficients, in-place from its coefficients given in the standard monomial basis to the Newton basis for the roots r_0 r_1 ldots r_n-2. In other words, this determines output coefficients c_i such that c_0 + c_1(x-r_0) + c_2(x-r_0)(x-r_1) + ldots + c_n-1(x-r_0)(x-r_1)cdots(x-r_n-2) is equal to the input polynomial.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.newton_to_monomial!-Union{Tuple{Array{T,1},Array{T,1}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.newton_to_monomial!",
    "category": "Method",
    "text": "newton_to_monomial!{T <: RingElement}(P::Array{T, 1}, roots::Array{T, 1})\n\nConverts a polynomial p, given as an array of coefficients, in-place from its coefficients given in the Newton basis for the roots r_0 r_1 ldots r_n-2 to the standard monomial basis. In other words, this evaluates c_0 + c_1(x-r_0) + c_2(x-r_0)(x-r_1) + ldots + c_n-1(x-r_0)(x-r_1)cdots(x-r_n-2) where c_i are the input coefficients given by p.\n\n\n\n"
},

{
    "location": "polynomial.html#Newton-representation-1",
    "page": "Generic univariate polynomials",
    "title": "Newton representation",
    "category": "section",
    "text": "monomial_to_newton!{T <: RingElem}(::Array{T, 1}, ::Array{T, 1})newton_to_monomial!{T <: RingElem}(::Array{T, 1}, ::Array{T, 1})ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = 3x*y^2 + (x + 1)*y + 3\ng = deepcopy(f)\nroots = [R(1), R(2), R(3)]\n\nmonomial_to_newton!(g.coeffs, roots)\nnewton_to_monomial!(g.coeffs, roots)"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.interpolate-Union{Tuple{AbstractAlgebra.PolyRing,Array{T,1},Array{T,1}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.interpolate",
    "category": "Method",
    "text": "interpolate{T <: RingElement}(S::AbstractAlgebra.PolyRing, x::Array{T, 1}, y::Array{T, 1})\n\nGiven two arrays of values xs and ys of the same length n, find the polynomial f in the polynomial ring R of length at most n such that f has the value ys at the points xs. The values in the arrays xs and ys must belong to the base ring of the polynomial ring R. If no such polynomial exists, an exception is raised.\n\n\n\n"
},

{
    "location": "polynomial.html#Interpolation-1",
    "page": "Generic univariate polynomials",
    "title": "Interpolation",
    "category": "section",
    "text": "interpolate{T <: RingElem}(::PolyRing, ::Array{T, 1}, ::Array{T, 1})ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nxs = [R(1), R(2), R(3), R(4)]\nys = [R(1), R(4), R(9), R(16)]\n\nf = interpolate(S, xs, ys)"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.chebyshev_t-Tuple{Int64,AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.chebyshev_t",
    "category": "Method",
    "text": "chebyshev_t(n::Int, x::AbstractAlgebra.PolyElem)\n\nReturn the Chebyshev polynomial of the first kind T_n(x), defined by T_n(x) = cos(n cos^-1(x)).\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.chebyshev_u-Tuple{Int64,AbstractAlgebra.PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.chebyshev_u",
    "category": "Method",
    "text": "chebyshev_u(n::Int, x::AbstractAlgebra.PolyElem)\n\nReturn the Chebyshev polynomial of the first kind U_n(x), defined by (n+1) U_n(x) = T_n+1(x).\n\n\n\n"
},

{
    "location": "polynomial.html#Special-functions-1",
    "page": "Generic univariate polynomials",
    "title": "Special functions",
    "category": "section",
    "text": "The following special functions can be computed for any polynomial ring. Typically one uses the generator x of a polynomial ring to get the respective special polynomials expressed in terms of that generator.chebyshev_t(::Int, ::PolyElem)chebyshev_u(::Int, ::PolyElem)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = chebyshev_t(20, y)\ng = chebyshev_u(15, y)"
},

{
    "location": "mpolynomial_rings.html#",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Multvariate Polynomial Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "mpolynomial_rings.html#Multvariate-Polynomial-Ring-Interface-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Multvariate Polynomial Ring Interface",
    "category": "section",
    "text": "Multivariate polynomial rings are supported in AbstractAlgebra.jl, and in addition to the standard Ring interface, numerous additional functions are provided.Unlike other kinds of rings, even complex operations such as GCD depend heavily on the multivariate representation. Therefore AbstractAlgebra.jl cannot provide much in the way of additional functionality to external multivariate implementations.This means that external libraries must be able to implement their multivariate formats in whatever way they see fit. The required interface here should be implemented, even if it is not optimal. But it can be extended, either by implementing one of the optional interfaces, or by extending the required interface in some other way.Naturally, any multivariate polynomial ring implementation provides the full Ring interface, in order to be treated as a ring for the sake of AbstractAlgebra.jl.Considerations which make it impossible for AbstractAlgebra.jl to provide generic functionality on top of an arbitrary multivariate module include:orderings (lexical, degree, weighted, block, arbitrary)\nsparse or dense representation\ndistributed or recursive representation\npacked or unpacked exponents\nexponent bounds (and whether adaptive or not)\nrandom access or iterators\nwhether monomials and polynomials have the same type\nwhether special cache aware data structures such as Geobuckets are used"
},

{
    "location": "mpolynomial_rings.html#Types-and-parents-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra.jl provides two abstract types for multivariate polynomial rings and their elements:MPolyRing{T} is the abstract type for multivariate polynomial ring parent types\nMPolyElem{T} is the abstract type for multivariate polynomial typesWe have that MPolyRing{T} <: AbstractAlgebra.Ring and  MPolyElem{T} <: AbstractAlgebra.RingElem.Note that both abstract types are parameterised. The type T should usually be the type of elements of the coefficient ring of the polynomial ring. For example, in the case of mathbbZx y the type T would be the type of an integer, e.g. BigInt.Multivariate polynomial rings should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Multivariate polynomial rings should at least be distinguished based on their base (coefficient) ring and number of variables. But if they have the same base ring, symbols (for their variables/generators) and ordering, they should certainly have the same parent object.See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "mpolynomial_rings.html#Required-functionality-for-multivariate-polynomials-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Required functionality for multivariate polynomials",
    "category": "section",
    "text": "In addition to the required functionality for the Ring interface, the Multivariate Polynomial interface has the following required functions.We suppose that R is a fictitious base ring (coefficient ring) and that S is a multivariate polynomial ring over R (i.e. S = Rx y ldots) with parent object S of type MyMPolyRing{T}. We also assume the polynomials in the ring have type MyMPoly{T}, where T is the type of elements of the base (coefficient) ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem."
},

{
    "location": "mpolynomial_rings.html#Data-type-and-parent-object-methods-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "vars(S::MyMPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn an array of Symbols representing the variables (generators) of the polynomial ring. Note that these are Symbols not Strings, though their string values will usually be used when printing polynomials.gens(S::MyMPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn an array of all the generators (variables) of the given polynomial ring (as polynomials).The first entry in the array will be the variable with most significance with respect to the ordering.ordering(S::MyMPolyRing{T})Return the ordering of the given polynomial ring as a symbol. Supported values currently include :lex, :deglex and :degrevlex.ExamplesS, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"]; ordering=:deglex)\n\nV = vars(S)\nX = gens(S)\nord = ordering(S)"
},

{
    "location": "mpolynomial_rings.html#Basic-manipulation-of-rings-and-elements-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "length(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn the number of nonzero terms of the given polynomial. The length of the zero polynomial is defined to be 0. The return value should be of type Int.isgen(x::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn true if x is a generator of the polynomial ring.max_degrees(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturns a tuple (B, b) consisting of an array of Ints specifying the highest power of each variable that appears in the given polynomial and b the largest of the values in B.nvars(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn the number of variables of the polynomial ring that f belongs to.isunit(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn true if f is a unit in its parent polynomial ring.isconstant(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn true if f is a constant polynomial. The zero polynomial is considered constant for the purposes of this function.isterm(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn true if f consists of a single term.ismonomial(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn true if f consists of a single term with coefficient 1.ExamplesS, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = x^3*y + 3x*y^2 + 1\n\nn = length(f)\nisgen(y) == true\nB, b = max_degrees(f)\nnvars(f) == 2\nisunit(f) == false\nisconstant(f) == false\nisterm(2x*y) == true\nismonomial(x*y) == false"
},

{
    "location": "mpolynomial_rings.html#Exact-division-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Exact division",
    "category": "section",
    "text": "For any ring that implements exact division, the following can be implemented.divexact(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn the exact quotient of f by g if it exists, otherwise throw an error.divides(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn a tuple (flag, q) where flag is true if g divides f, in which case q will be the exact quotient, or flag is false and q is set to zero.remove(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturns a tuple (v q) such that the highest power of g that divides f is g^v and the cofactor is q.valuation(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturns v such that the highest power of g that divides f is g^v.ExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = 2x^2*y + 2x + y + 1\ng = x^2*y^2 + 1\n\nflag, q = divides(f*g, f)\nd = divexact(f*g, f)\nv, q = remove(f*g^3, g)\nn = valuation(f*g^3, g)"
},

{
    "location": "mpolynomial_rings.html#Euclidean-division-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Euclidean division",
    "category": "section",
    "text": "Although multivariate polynomial rings are not in general Euclidean, it is possible to define a quotient with remainder function that depends on the polynomial ordering in the case that the quotient ring is a field or a Euclidean domain. In the case that a polynomial g divides a polynomial f, the result no longer depends on the ordering and the remainder is zero, with the quotient agreeing with the exact quotient.divrem(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn a tuple (q r) such that f = qg + r, where the coefficients of terms of r whose monomials are divisible by the leading monomial of g are reduced modulo the leading coefficient of g (according to the Euclidean function on the coefficients).Note that the result of this function depends on the ordering of the polynomial ring.div(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemAs per the divrem function, but returning the quotient only. Especially when the quotient happens to be exact, this function can be exceedingly fast.divrem(f::MyMPoly{T}, G::Array{MyMPoly{T}, 1}) where T <: AbstractAlgebra.RingElemAs per the divrem function above, except that each term of r starting with the most significant term, is reduced modulo the leading terms of each of the polynomials in the array G for which the leading monomial is a divisor.A tuple (Q r) is returned from the function, where Q is an array of polynomials of the same length as G, and such that f = r + sum QiGi.The result is again dependent on the ordering in general, but if the polynomials in G are over a field and the reduced generators of a Groebner basis, then the result is unique.ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nf = 2x^2*y + 2x + y + 1\ng = x + y\nh = y + 1\n\nq = div(f, g)\nq, r = divrem(f, g)\nQ, r = divrem(f, [g, h])"
},

{
    "location": "mpolynomial_rings.html#Evaluation-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Evaluation",
    "category": "section",
    "text": "evaluate(f::MyMPoly{T}, A::Array{T, 1}) where T <: AbstractAlgebra.RingElemEvaluate the polynomial f at the values specified by the entries of the array A.evaluate(f::MPoly{T}, A::Array{T, 1}) where T <: IntegerEvaluate the polynomial f at the values specified by the entries of the array A.ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nf = 2x^2*y + 2x + y + 1\n\nm = evaluate(f, Rational{BigInt}[2, 3])\nn = evaluate(f, [2, 3])"
},

{
    "location": "mpolynomial_rings.html#GCD-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "GCD",
    "category": "section",
    "text": "In cases where there is a meaningful Euclidean structure on the coefficient ring, it is possible to compute the GCD of multivariate polynomials.gcd(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn a greatest common divisor of f and g.ExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = 2x^2*y + 2x + y + 1\ng = x^2*y^2 + 1\n\nd = gcd(f*g^2, f^2*g)"
},

{
    "location": "mpolynomial_rings.html#Interface-for-sparse-distributed,-random-access-multivariates-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Interface for sparse distributed, random access multivariates",
    "category": "section",
    "text": "The following additional functions should be implemented by libraries that provide a sparse distributed polynomial format, stored in a representation for which terms can be accessed in constant time (e.g. where arrays are used to store coefficients and exponent vectors)."
},

{
    "location": "mpolynomial_rings.html#Sparse-distributed,-random-access-constructors-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Sparse distributed, random access constructors",
    "category": "section",
    "text": "In addition to the standard constructors, the following constructor, taking arrays of coefficients and exponent vectors, should be provided.(S::MyMPolyRing{T})(A::Array{T, 1}, m::Array{UInt, 2}) where T <: AbstractAlgebra.RingEle\nmCreate the polynomial in the given ring with nonzero coefficients specified by the elements of A and corresponding exponent vectors given by the elements of m. For efficiency reason, the exponents of term i are given by the vector m[:, i] since Julia uses column major two dimensional arrays.For maximum compatibility with external libraries, the coefficient (and term) at index 1 correspond to the most significant term with respect to the polynomial ring ordering.Each exponent vector uses a separate word for each exponent field, the first of which should be any degree or weight, and otherwise should be the exponent for the most significant variable with respect to the ordering. The top bit of each word is reserved to detect overflows.If a full word is not used for exponents, a check should be done to ensure there are no overflows before setting the exponents.A library may also optionally provide an interface that makes use of BigInt (or any other big integer type) for exponents instead of UInt.ExamplesS, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nf = S(Rational{BigInt}[2, 3, 1], UInt[3 2 1; 0 1 0])"
},

{
    "location": "mpolynomial_rings.html#Sparse-distributed,-random-access-basic-manipulation-1",
    "page": "Multvariate Polynomial Ring Interface",
    "title": "Sparse distributed, random access basic manipulation",
    "category": "section",
    "text": "coeff(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemReturn the coefficient of the (n+1)-th term of f. The first term should be the most significant term with respect to the ordering.exponent(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemReturn an array of Ints giving the vector of exponents for the n + 1-th term of f. The first entry of the array should correspond to the exponent of the most significant variable with respect to the ordering.exponent!(A::Array{Int, 1}, f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemAs per exponent, but set the values in the array A rather than allocating an array for this purpose. The array is also returned by the function after being mutated.fit!(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemEnsure that the polynomial f internally has space for n nonzero terms. This function must mutate the function in-place if it is mutable. It does not return the mutated polynomial. Immutable types can still be supported by defining this function to do nothing.ExamplesS, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = x^3*y + 3x*y^2 + 1\n\nc = coeff(f, 1)\nfit!(f, 8)"
},

{
    "location": "mpolynomial.html#",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Generic sparse distributed multivariate polynomials",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "mpolynomial.html#Generic-sparse-distributed-multivariate-polynomials-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Generic sparse distributed multivariate polynomials",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/generic/MPoly.jl for generic sparse distributed multivariate polynomials over any commutative ring belonging to the AbstractAlgebra abstract type hierarchy.This modules implements the Multivariate Polynomial interface, including the sparse distributed, random access part of the interface.All of the generic functionality is part of a submodule of AbstractAlgebra called Generic. This is exported by default so that it is not necessary to qualify the function names with the submodule name.Multivariates are implemented in this module using a Julia array of coefficients and a 2-dimensional Julia array of UInts for the exponent vectors. Note that exponent n is represented by the n-th column of the exponent array, not the n-th row. This is because Julia uses a column major representation."
},

{
    "location": "mpolynomial.html#Types-and-parent-objects-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Multivariate polynomials implemented in AbstractAlgebra.jl have type Generic.MPoly{T} where T is the type of elements of the coefficient ring.The polynomials are implemented using a Julia array of coefficients and a 2-dimensional Julia array of UInts for the exponent vectors. Note that exponent n is represented by the n-th column of the exponent array, not the n-th row. This is because Julia uses a column major representation. See the file src/generic/GenericTypes.jl for details.The top bit of each UInt is reserved for overflow detection.Parent objects of such polynomials have type Generic.MPolyRing{T}.The string representation of the variables of the polynomial ring, the base/coefficient ring R and the ordering are stored in the parent object. The polynomial element types belong to the abstract type AbstractAlgebra.MPolyElem{T} and the polynomial ring types belong to the abstract type AbstractAlgebra.MPolyRing{T}.Note that both the generic polynomial ring type Generic.MPolyRing{T} and the abstract type it belongs to, AbstractAlgebra.MPolyRing{T} are both called MPolyRing. The  former is a (parameterised) concrete type for a polynomial ring over a given base ring whose elements have type T. The latter is an abstract type representing all multivariate polynomial ring types in AbstractAlgebra.jl, whether generic or very specialised (e.g. supplied by a C library)."
},

{
    "location": "mpolynomial.html#Polynomial-ring-constructors-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Polynomial ring constructors",
    "category": "section",
    "text": "In order to construct multivariate polynomials in AbstractAlgebra.jl, one must first construct the polynomial ring itself. This is accomplished with the following constructor.PolynomialRing(R::AbstractAlgebra.Ring, S::Array{String, 1}; cached::Bool = true, ordering::Symbol=:lex)Given a base ring R and and array S of strings specifying how the generators (variables) should be printed, return a tuple S, (x, ...) representing the new polynomial ring S = Rx ldots and a tuple of the generators (x ) of the ring. By default the parent object S will depend only on R and  (x, ...) and will be cached. Setting the optional argument cached to false will prevent the parent object  S from being cached.The optional named argument ordering can be used to specify an ordering. The currently supported options are :lex, :deglex and `:degrevlex	.Here are some examples of creating multivariate polynomial rings and making use of the resulting parent objects to coerce various elements into the polynomial ring.ExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"]; ordering=:deglex)\n\nf = R()\ng = R(123)\nh = R(BigInt(1234))\nk = R(x + 1)\nm = R(x + y + 1)All of the examples here are generic polynomial rings, but specialised implementations of polynomial rings provided by external modules will also usually provide a PolynomialRing constructor to allow creation of their polynomial rings."
},

{
    "location": "series_rings.html#",
    "page": "Series Ring Interface",
    "title": "Series Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "series_rings.html#Series-Ring-Interface-1",
    "page": "Series Ring Interface",
    "title": "Series Ring Interface",
    "category": "section",
    "text": "Univariate power series rings are supported in AbstractAlgebra in a variety of different forms, including absolute and relative precision models and Laurent series.In addition to the standard Ring interface, numerous additional functions are required to be present for power series rings."
},

{
    "location": "series_rings.html#Types-and-parents-1",
    "page": "Series Ring Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra provides two abstract types for power series rings and their elements:SeriesRing{T} is the abstract type for all power series ring parent types\nSeriesElem{T} is the abstract type for all power series typesWe have that SeriesRing{T} <: AbstractAlgebra.Ring and  SeriesElem{T} <: AbstractAlgebra.RingElem.Note that both abstract types are parameterised. The type T should usually be the type of elements of the coefficient ring of the power series ring. For example, in the case of mathbbZx the type T would be the type of an integer, e.g. BigInt.Within the SeriesElem{T} abstract type is the abstract type RelSeriesElem{T} for relative power series, and AbsSeriesElem{T} for absolute power series.Relative series are typically stored with a valuation and a series that is either zero or that has nonzero constant term. Absolute series are stored starting from the constant term, even if it is zero.If the parent object for a relative series ring over the bignum integers has type MySeriesRing and series in that ring have type MySeries then one would have:MySeriesRing <: SeriesRing{BigInt}\nMySeries <: RelSeriesElem{BigInt}Series rings should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Series rings should at least be distinguished based on their base (coefficient) ring. But if they have the same base ring and symbol (for their variable/generator) and same default precision, they should certainly have the same parent object.See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "series_rings.html#Required-functionality-for-series-1",
    "page": "Series Ring Interface",
    "title": "Required functionality for series",
    "category": "section",
    "text": "In addition to the required functionality for the Ring interface the Series Ring interface has the following required functions.We suppose that R is a fictitious base ring (coefficient ring) and that S is a series ring over R (e.g. S = Rx) with parent object S of type MySeriesRing{T}. We also assume the series in the ring have type MySeries{T}, where  T is the type of elements of the base (coefficient) ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem."
},

{
    "location": "series_rings.html#Constructors-1",
    "page": "Series Ring Interface",
    "title": "Constructors",
    "category": "section",
    "text": "In addition to the standard constructors, the following constructors, taking an array of coefficients, must be available.For relative power series and Laurent series we have:(S::MySeriesRing{T})(A::Array{T, 1}, len::Int, prec::Int, val::Int) where T <: AbstractAlgebra.RingElemCreate the series in the given ring whose valuation is val, whose absolute precision is given by prec and the coefficients of which are given by A, starting from the first nonzero term. Only len terms of the array are used, the remaining terms being ignored. The value len cannot exceed the length of the supplied array.It is permitted to have trailing zeros in the array, but it is not needed, even if the precision minus the valuation is bigger than the length of the array.ExamplesS, x = PowerSeriesRing(QQ, 10, \"x\"; model=:capped_relative)\nT, y = LaurentSeriesRing(ZZ, 10, \"y\")\nU, z = LaurentSeriesField(QQ, 10, \"z\")\n \nf = S(Rational{BigInt}[2, 3, 1], 3, 6, 2)\ng = T(BigInt[2, 3, 1], 3, 6, 2)\nh = U(Rational{BigInt}[2, 3, 1], 3, 6, 2)For absolute power series we have:(S::MySeriesRing{T})(A::Array{T, 1}, len::Int, prec::Int) where T <: AbstractAlgebra.RingElemCreate the series in the given ring whose absolute precision is given by prec and the coefficients of which are given by A, starting from the constant term. Only len terms of the array are used, the remaining terms being ignored.Note that len is usually maintained separately of any polynomial that is underlying the power series. This allows for easy trucation of a power series without actually modifying the polynomial underlying it.It is permitted to have trailing zeros in the array, but it is not needed, even if the precision is bigger than the length of the array.ExamplesS, x = PowerSeriesRing(QQ, 10, \"x\"; model=:capped_absolute)\n\nf = S(Rational{BigInt}[0, 2, 3, 1], 4, 6)"
},

{
    "location": "series_rings.html#Data-type-and-parent-object-methods-1",
    "page": "Series Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "var(S::MySeriesRing{T}) where T <: AbstractAlgebra.RingElemReturn a Symbol representing the variable (generator) of the series ring. Note that this is a Symbol not a String, though its string value will usually be used when printing series.max_precision(S::MySeriesRing{T}) where T <: AbstractAlgebra.RingElemReturn the (default) maximum precision of the power series ring. This is the precision that the output of an operation will be if it cannot be represented to full precision (e.g. because it mathematically has infinite precision).This value is usually supplied upon creation of the series ring and stored in the ring. It is independent of the precision which each series in the ring actually has. Those are stored on a per element basis in the actual series elements.ExamplesS, x = PowerSeriesRing(QQ, 10, \"x\")\n\nvsym = var(S)\nmax_precision(S) == 10"
},

{
    "location": "series_rings.html#Basic-manipulation-of-rings-and-elements-1",
    "page": "Series Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "pol_length(f::MySeries{T}) where T <: AbstractAlgebra.RingElemReturn the length of the polynomial underlying the given power series. This is not generally useful to the user, but is used internally.set_length!(f::MySeries{T}, n::Int) where T <: AbstractAlgebra.RingElemThis function sets the effective length of the polynomial underlying the given series. The function doesn't modify the actual polynomial, but simply changes the number of terms of the polynomial which are considered to belong to the power series. The remaining terms are ignored.This function cannot set the length to a value greater than the length of any underlying polynomial.The function mutates the series in-place but does not return the mutated series.precision(f::MySeries{T})Returns the absolute precision of f.set_prec!(f::MySeries{T}, prec::Int)Set the absolute precision of the given series to the given value.This function mutates the series in-place but does not return the mutated series.valuation(f::MySeries{T})Return the valuation of the given series.set_val!(f::MySeries{T}, val::Int)For relative series and Laurent series only, this function alters the valuation of the given series to the given value.The series is mutated in-place but does not return the mutated series.polcoeff(f::MySeries{T}, n::Int) Return the coefficient of degree n of the polynomial underlying the series. If n is larger than the degree of this polynomial, zero is returned. This function is not generally of use to the user but is used internally.setcoeff!(f::MySeries{T}, n::Int, a::T) where T <: AbstractAlgebra.RingElemSet the degree n coefficient of the polynomial underlying f to a. This mutates the polynomial in-place if possible and returns the mutated series (so that immutable types can also be supported). The function must not assume that the polynomial already has space for n + 1 coefficients. The polynomial must be resized if this is not the case.Note that this function is not required to normalise the polynomial and is not necessarily useful to the user, but is used extensively by the generic functionality in AbstractAlgebra.jl. It is for setting raw coefficients in the representation.normalise(f::MySeries{T}, n::Int)Given a series f represented by a polynomial of at least the given length, return the normalised length of the underlying polynomial assuming it has length at most n. This function does not actually normalise the polynomial and is not particularly useful to the user. It is used internally.renormalize!(f::MySeries{T}) where T <: AbstractAlgebra.RingElemGiven a relative series or Laurent series whose underlying polynomial has zero constant term, say as the result of some internal computation, renormalise the series so that the  polynomial has nonzero constant term. The precision and valuation of the series are adjusted to compensate. This function is not intended to be useful to the user, but is  used internally.fit!(f::MySeries{T}, n::Int) where T <: AbstractAlgebra.RingElemEnsure that the polynomial underlying f internally has space for n coefficients. This function must mutate the series in-place if it is mutable. It does not return the mutated series. Immutable types can still be supported by defining this function to do nothing.Some interfaces for C polynomial types automatically manage the internal allocation of polynomials in every function that can be called on them. Explicit adjustment by the generic code in AbstractAlgebra.jl is not required. In such cases, this function can also be defined to do nothing.gen(R::MySeriesRing{T}) where T <: AbstractAlgebra.RingElemReturn the generator x of the series ring.ExamplesS, x = PowerSeriesRing(ZZ, 10, \"x\")\n\nf = 1 + 3x + x^3 + O(x^5)\ng = S(BigInt[1, 2, 0, 1, 0, 0, 0], 4, 10, 3);\n\nn = pol_length(f)\nc = polcoeff(f, 1)\nset_length!(g, 3)\ng = setcoeff!(g, 2, BigInt(11))\nfit!(g, 8)\ng = setcoeff!(g, 7, BigInt(4))\nw = gen(S)\nisgen(w) == true"
},

{
    "location": "series.html#",
    "page": "Generic power series",
    "title": "Generic power series",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "series.html#Generic-power-series-1",
    "page": "Generic power series",
    "title": "Generic power series",
    "category": "section",
    "text": "AbstractAlgebra.jl allows the creation of capped relative and absolute power series over  any computable commutative ring R.Capped relative power series are power series of the form a_jx^j + a_j+1x^j+1 + cdots + a_k-1x^k-1 + O(x^k) where i geq 0, a_i in R and the relative precision k - j is at most equal to some specified precision n.Capped absolute power series are power series of the form a_jx^j + a_j+1x^j+1 + cdots + a_n-1x^n-1 + O(x^n) where j geq 0, a_j in R and the precision n is fixed.There are two implementations of relative series: relative power series, implemented in src/generic/RelSeries.jl and Laurent series, implemented in src/generic/Laurent.jl. Note that there are two implementations for Laurent series, one over rings and one over fields, though in practice most of the implementation uses the same code in both cases.There is a single implementation of absolute series: absolute power series, implemented in src/generic/AbsSeries.jl.As well as implementing the Series Ring interface, the series modules in AbstractAlgebra.jl implement the generic algorithms described below.All of the generic functionality is part of the Generic submodule of AbstractAlgebra.jl. This is exported by default so that it is not necessary to qualify function names."
},

{
    "location": "series.html#Types-and-parent-objects-1",
    "page": "Generic power series",
    "title": "Types and parent objects",
    "category": "section",
    "text": "The types of generic polynomials implemented by AbstractAlgebra.jl are Generic.RelSeries{T}, Generic.AbsSeries{T}, LaurentSeriesRingElem{T} and LaurentSeriesFieldElem{T}.Relative power series elements belong to the abstract type AbstractAlgebra.RelSeriesElem.Laurent series elements belong directly to either AbstractAlgebra.RingElem or AbstractAlgebra.FieldElem since it is more useful to be able to distinguish whether they belong to a ring or field than it is to distinguish that they are relative series.Absolute power series elements belong to AbstractAlgebra.AbsSeriesElem.The parent types for relative and absolute power series, Generic.RelSeriesRing{T}  and Generic.AbsSeriesRing{T} respectively, belong to AbstractAlgebra.SeriesRing{T}.The parent types for Laurent series rings and fields, Generic.LaurentSeriesRing{T} and Generic.LaurentSeriesField{T} respectively, belong directly to  AbstractAlgebra.Ring and AbstractAlgebra.Field respectively.The default precision, string representation of the variable and base ring R of a generic power series are stored in its parent object. "
},

{
    "location": "series.html#Series-ring-constructors-1",
    "page": "Generic power series",
    "title": "Series ring constructors",
    "category": "section",
    "text": "In order to construct series in AbstractAlgebra.jl, one must first construct the ring itself. This is accomplished with any of the following constructors.PowerSeriesRing(R::AbstractAlgebra.Ring, prec_max::Int, s::AbstractString; cached::Bool = true, model=:capped_relative)LaurentSeriesRing(R::AbstractAlgebra.Ring, prec_max::Int, s::AbstractString; cached::Bool = true)LaurentSeriesRing(R::AbstractAlgebra.Field, prec_max::Int, s::AbstractString; cached::Bool = true)Given a base ring R, a maximum precision (relative or absolute, depending on the model) and a string s specifying how the generator (variable) should be printed, return a typle S, x representing the series ring and its generator.By default, S will depend only on S, x and the maximum precision and will be cached. Setting the optional argument cached to false will prevent this.In the case of power series, the optional argument model can be set to either :capped_absolute or capped_relative, depending on which power series model is required.Here are some examples of constructing various kinds of series rings and coercing various elements into those rings.ExamplesR, x = PowerSeriesRing(ZZ, 10, \"x\")\nS, y = PowerSeriesRing(ZZ, 10, \"y\"; model=:capped_absolute)\nT, z = LaurentSeriesRing(ZZ, 10, \"z\")\nU, w = LaurentSeriesField(QQ, 10, \"w\")\n\nf = R()\ng = S(123)\nh = U(BigInt(1234))\nk = T(z + 1)"
},

{
    "location": "series.html#Big-oh-notation-1",
    "page": "Generic power series",
    "title": "Big-oh notation",
    "category": "section",
    "text": "Series elements can be given a precision using the big-oh notation. This is provided by a function of the following form, (or something equivalent for Laurent series):O(x::SeriesElem)ExamplesR, x = PowerSeriesRing(ZZ, 10, \"x\")\nS, y = LaurentSeriesRing(ZZ, 10, \"y\")\n\nf = 1 + 2x + O(x^5)\ng = 2y + 7y^2 + O(y^7)What is happening here in practice is that O(x^n) is creating the series 0 + O(x^n) and the rules for addition of series dictate that if this is added to a series of  greater precision, then the lower of the two precisions must be used.Of course it may be that the precision of the series that O(x^n) is added to is already lower than n, in which case adding O(x^n) has no effect. This is the case if the default precision is too low, since x on its own has the default precision."
},

{
    "location": "series.html#Power-series-models-1",
    "page": "Generic power series",
    "title": "Power series models",
    "category": "section",
    "text": "Capped relative power series have their maximum relative precision capped at some value prec_max. This means that if the leading term of a nonzero power series element is c_ax^a and the precision is b then the power series is of the form  c_ax^a + c_a+1x^a+1 + ldots + O(x^a + b).The zero power series is simply taken to be 0 + O(x^b).The capped relative model has the advantage that power series are stable multiplicatively. In other words, for nonzero power series f and g we have that divexact(f*g), g) == f.However, capped relative power series are not additively stable, i.e. we do not always have (f + g) - g = f.Similar comments apply to Laurent series.On the other hand, capped absolute power series have their absolute precision capped. This means that if the leading term of a nonzero power series element is c_ax^a and the precision is b then the power series is of the form c_ax^a + c_a+1x^a+1 + ldots + O(x^b).Capped absolute series are additively stable, but not necessarily multiplicatively stable.For all models, the maximum precision is also used as a default precision in the case of coercing coefficients into the ring and for any computation where the result could mathematically be given to infinite precision.In all models we say that two power series are equal if they agree up to the minimum absolute precision of the two power series.Thus, for example, x^5 + O(x^10) == 0 + O(x^5), since the minimum absolute precision is 5.During computations, it is possible for power series to lose relative precision due to cancellation. For example if f = x^3 + x^5 + O(x^8) and g = x^3 + x^6 + O(x^8) then f - g = x^5 - x^6 + O(x^8) which now has relative precision 3 instead of relative precision 5.Amongst other things, this means that equality is not transitive. For example x^6 + O(x^11) == 0 + O(x^5) and x^7 + O(x^12) == 0 + O(x^5) but x^6 + O(x^11) neq x^7 + O(x^12).Sometimes it is necessary to compare power series not just for arithmetic equality, as above, but to see if they have precisely the same precision and terms. For this purpose we introduce the isequal function.For example, if f = x^2 + O(x^7) and g = x^2 + O(x^8) and h = 0 + O(x^2) then f == g, f == h and g == h, but isequal(f, g), isequal(f, h) and isequal(g, h) would all return false. However, if k = x^2 + O(x^7) then isequal(f, k) would return true.There are further difficulties if we construct polynomial over power series. For example, consider the polynomial in y over the power series ring in x over the rationals. Normalisation of such polynomials is problematic. For instance, what is the leading coefficient of (0 + O(x^10))y + (1 + O(x^10))?If one takes it to be (0 + O(x^10)) then some functions may not terminate due to the fact that algorithms may require the degree of polynomials to decrease with each iteration. Instead, the degree may remain constant and simply accumulate leading terms which are arithmetically zero but not identically zero.On the other hand, when constructing power series over other power series, if we simply throw away terms which are arithmetically equal to zero, our computations may have different output depending on the order in which the power series are added!One should be aware of these difficulties when working with power series. Power series, as represented on a computer, simply don't satisfy the axioms of a ring. They must be used with care in order to approximate operations in a mathematical power series ring.Simply increasing the precision will not necessarily give a \"more correct\" answer and some computations may not even terminate due to the presence of arithmetic zeroes!An absolute power series ring over a ring R with precision p behaves  very much like the quotient Rx(x^p) of the polynomial ring over R. Therefore one can often treat absolute power series rings as though they were rings. However, this depends on all series being given a precision equal to the specified maximum precision and not a lower precision."
},

{
    "location": "series.html#Basic-ring-functionality-1",
    "page": "Generic power series",
    "title": "Basic ring functionality",
    "category": "section",
    "text": "All power series models provide the functionality described in the Ring and Series Ring interfaces.ExamplesS, x = PowerSeriesRing(ZZ, 10, \"x\")\n\nf = 1 + 3x + x^3 + O(x^10)\ng = 1 + 2x + x^2 + O(x^10)\n\nh = zero(S)\nk = one(S)\nisone(k) == true\niszero(f) == false\nn = pol_length(f)\nc = polcoeff(f, 3)\nU = base_ring(S)\nv = var(S)\nT = parent(x + 1)\ng == deepcopy(g)\nt = divexact(2g, 2)\np = precision(f)"
},

{
    "location": "series.html#Series-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic power series",
    "title": "Series functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality below is automatically provided by AbstractAlgebra.jl for any series module that implements the full Series Ring interface. This includes AbstractAlgebra's own generic series rings.Of course, modules are encouraged to provide specific implementations of the functions described here, that override the generic implementation.Unless otherwise noted, the functions are available for all series models, including Laurent series. We denote this by using the abstract type AbstractAlgebra.RelSeriesElem, even though absolute series and Laurent series types do not belong to this abstract type."
},

{
    "location": "series.html#AbstractAlgebra.Generic.modulus-Union{Tuple{AbstractAlgebra.SeriesElem{T}}, Tuple{T}} where T<:AbstractAlgebra.ResElem",
    "page": "Generic power series",
    "title": "AbstractAlgebra.Generic.modulus",
    "category": "Method",
    "text": "modulus{T <: ResElem}(a::AbstractAlgebra.SeriesElem{T})\n\nReturn the modulus of the coefficients of the given power series.\n\n\n\n"
},

{
    "location": "series.html#AbstractAlgebra.Generic.isgen-Tuple{AbstractAlgebra.RelSeriesElem}",
    "page": "Generic power series",
    "title": "AbstractAlgebra.Generic.isgen",
    "category": "Method",
    "text": "isgen(a::RelSeriesElem)\n\nReturn true if the given power series is arithmetically equal to the generator of its power series ring to its current precision, otherwise return false.\n\n\n\n"
},

{
    "location": "series.html#AbstractAlgebra.Generic.isunit-Tuple{AbstractAlgebra.RelSeriesElem}",
    "page": "Generic power series",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "Method",
    "text": "isunit(a::AbstractAlgebra.RelSeriesElem)\n\nReturn true if the given power series is arithmetically equal to a unit, i.e. is invertible, otherwise return false.\n\n\n\n"
},

{
    "location": "series.html#Basic-functionality-1",
    "page": "Generic power series",
    "title": "Basic functionality",
    "category": "section",
    "text": "coeff(a::AbstractAlgebra.SeriesElem, n::Int)Return the degree n coefficient of the given power series. Note coefficients are numbered from n = 0 for the constant coefficient. If n exceeds the current precision of the power series, the function returns a zero coefficient.For power series types, n must be non-negative. Laurent series do not have this restriction.modulus{T <: ResElem}(::SeriesElem{T})isgen(::RelSeriesElem)isunit(::RelSeriesElem)ExamplesR, t = PowerSeriesRing(QQ, 10, \"t\")\nS, x = PowerSeriesRing(R, 30, \"x\")\n\na = O(x^4)\nb = (t + 3)*x + (t^2 + 1)*x^2 + O(x^4)\n\nk = isgen(gen(R))\nm = isunit(-1 + x + 2x^2)\nn = valuation(a)\np = valuation(b)\nc = coeff(b, 2)"
},

{
    "location": "series.html#AbstractAlgebra.Generic.shift_left-Union{Tuple{AbstractAlgebra.RelSeriesElem{T},Int64}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic power series",
    "title": "AbstractAlgebra.Generic.shift_left",
    "category": "Method",
    "text": "shift_left(x::AbstractAlgebra.RelSeriesElem, n::Int)\n\nReturn the power series f shifted left by n terms, i.e. multiplied by x^n.\n\n\n\n"
},

{
    "location": "series.html#AbstractAlgebra.Generic.shift_right-Union{Tuple{AbstractAlgebra.RelSeriesElem{T},Int64}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic power series",
    "title": "AbstractAlgebra.Generic.shift_right",
    "category": "Method",
    "text": "shift_right(f::AbstractAlgebra.RelSeriesElem, n::Int)\n\nReturn the power series f shifted right by n terms, i.e. divided by x^n.\n\n\n\n"
},

{
    "location": "series.html#Shifting-1",
    "page": "Generic power series",
    "title": "Shifting",
    "category": "section",
    "text": "shift_left{T <: RingElem}(::RelSeriesElem{T}, ::Int)shift_right{T <: RingElem}(::RelSeriesElem{T}, ::Int)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS, x = PowerSeriesRing(R, 30, \"x\")\n\na = 2x + x^3\nb = O(x^4)\nc = 1 + x + 2x^2 + O(x^5)\nd = 2x + x^3 + O(x^4)\n\nf = shift_left(a, 2)\ng = shift_left(b, 2)\nh = shift_right(c, 1)\nk = shift_right(d, 3)"
},

{
    "location": "series.html#Base.truncate-Union{Tuple{AbstractAlgebra.RelSeriesElem{T},Int64}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic power series",
    "title": "Base.truncate",
    "category": "Method",
    "text": "truncate(a::AbstractAlgebra.RelSeriesElem, n::Int)\n\nReturn a truncated to (absolute) precision n.\n\n\n\n"
},

{
    "location": "series.html#Truncation-1",
    "page": "Generic power series",
    "title": "Truncation",
    "category": "section",
    "text": "truncate{T <: RingElem}(::RelSeriesElem{T}, ::Int)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS, x = PowerSeriesRing(R, 30, \"x\")\n\na = 2x + x^3\nb = O(x^4)\nc = 1 + x + 2x^2 + O(x^5)\nd = 2x + x^3 + O(x^4)\n\nf = truncate(a, 3)\ng = truncate(b, 2)\nh = truncate(c, 7)\nk = truncate(d, 5)"
},

{
    "location": "series.html#Base.inv-Tuple{AbstractAlgebra.RelSeriesElem}",
    "page": "Generic power series",
    "title": "Base.inv",
    "category": "Method",
    "text": "inv(a::AbstractAlgebra.RelSeriesElem)\n\nReturn the inverse of the power series a, i.e. 1a.\n\n\n\n"
},

{
    "location": "series.html#Division-1",
    "page": "Generic power series",
    "title": "Division",
    "category": "section",
    "text": "inv(::RelSeriesElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS, x = PowerSeriesRing(R, 30, \"x\")\n\na = 1 + x + 2x^2 + O(x^5)\nb = S(-1)\n\nc = inv(a)\nd = inv(b)"
},

{
    "location": "series.html#Base.exp-Tuple{AbstractAlgebra.RelSeriesElem}",
    "page": "Generic power series",
    "title": "Base.exp",
    "category": "Method",
    "text": "exp(a::AbstractAlgebra.RelSeriesElem)\n\nReturn the exponential of the power series a.\n\n\n\n"
},

{
    "location": "series.html#Base.sqrt-Tuple{AbstractAlgebra.RelSeriesElem}",
    "page": "Generic power series",
    "title": "Base.sqrt",
    "category": "Method",
    "text": "sqrt(a::AbstractAlgebra.RelSeriesElem)\n\nReturn the square root of the power series a.\n\n\n\n"
},

{
    "location": "series.html#Special-functions-1",
    "page": "Generic power series",
    "title": "Special functions",
    "category": "section",
    "text": "Base.exp(a::RelSeriesElem)Base.sqrt(a::RelSeriesElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS, x = PowerSeriesRing(R, 30, \"x\")\nT, z = PowerSeriesRing(QQ, 30, \"z\")\n\na = 1 + z + 3z^2 + O(z^5)\nb = z + 2z^2 + 5z^3 + O(z^5)\n\nc = exp(x + O(x^40))\nd = divexact(x, exp(x + O(x^40)) - 1)\nf = exp(b)\nh = sqrt(a)"
},

{
    "location": "residue_rings.html#",
    "page": "Residue Ring Interface",
    "title": "Residue Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "residue_rings.html#Residue-Ring-Interface-1",
    "page": "Residue Ring Interface",
    "title": "Residue Ring Interface",
    "category": "section",
    "text": "Residue rings (currently a quotient ring modulo a principal ideal) are supported in AbstractAlgebra.jl, at least for Euclidean base rings. In addition to the standard Ring interface, some additional functions are required to be present for residue rings."
},

{
    "location": "residue_rings.html#Types-and-parents-1",
    "page": "Residue Ring Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra provides four abstract types for residue rings and their elements:ResRing{T} is the abstract type for residue ring parent types\nResField{T} is the abstract type for residue rings known to be fields\nResElem{T} is the abstract type for types of elements of residue rings (residues)\nResFieldElem{T} is the abstract type for types of elements of residue fieldsWe have that ResRing{T} <: AbstractAlgebra.Ring and  ResElem{T} <: AbstractAlgebra.RingElem.Note that these abstract types are parameterised. The type T should usually be the type of elements of the base ring of the residue ring/field.If the parent object for a residue ring has type MyResRing and residues in that ring have type MyRes then one would have:MyResRing <: ResRing{BigInt}\nMyRes <: ResElem{BigInt}Residue rings should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Residue rings should at least be distinguished based on their base ring and modulus (the principal ideal one is taking a quotient of the base ring by).See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "residue_rings.html#Required-functionality-for-residue-rings-1",
    "page": "Residue Ring Interface",
    "title": "Required functionality for residue rings",
    "category": "section",
    "text": "In addition to the required functionality for the Ring interface the Residue Ring interface has the following required functions.We suppose that R is a fictitious base ring, m is an element of that ring, and that S is the residue ring (quotient ring) R(m) with parent object S of type MyResRing{T}. We also assume the residues r pmodm in the residue ring have type MyRes{T}, where T is the type of elements of the base ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem."
},

{
    "location": "residue_rings.html#Data-type-and-parent-object-methods-1",
    "page": "Residue Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "modulus(S::MyResRing{T}) where T <: AbstractAlgebra.RingElemReturn the modulus of the given residue ring, i.e. if the residue ring S was specified to be R(m), return m.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\n\nm = modulus(S)"
},

{
    "location": "residue_rings.html#Basic-manipulation-of-rings-and-elements-1",
    "page": "Residue Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "data(f::MyRes{T}) where T <: AbstractAlgebra.RingElemGiven a residue r pmodm, represented as such, return r.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\n\nf = S(x^2 + 2)\n\nd = data(f)"
},

{
    "location": "residue.html#",
    "page": "Generic residue rings",
    "title": "Generic residue rings",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "residue.html#Generic-residue-rings-1",
    "page": "Generic residue rings",
    "title": "Generic residue rings",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/generic/Residue.jl for generic residue rings over any Euclidean domain (in practice most of the functionality is provided for GCD domains that provide a meaningful GCD function) belonging to the AbstractAlgebra.jl abstract type hierarchy.As well as implementing the Residue Ring interface a number of generic algorithms are implemented for residue rings. We describe this generic functionality below.All of the generic functionality is part of a submodule of AbstractAlgebra called Generic. This is exported by default so that it is not necessary to qualify the function names with the submodule name."
},

{
    "location": "residue.html#Types-and-parent-objects-1",
    "page": "Generic residue rings",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Residues implemented using the AbstractAlgebra generics have type Generic.Res{T} or in the case of residue rings that are known to be fields, Generic.ResF{T}, where T is the type of elements of the base ring. See the file src/generic/GenericTypes.jl for details.Parent objects of residue ring elements have type Generic.ResRing{T} and those of residue fields have type GenericResField{T}.The defining modulus of the residue ring is stored in the parent object.The residue element types belong to the abstract type AbstractAlgebra.ResElem{T} or AbstractAlgebra.ResFieldElem{T} in the case of residue fields, and the residue ring types belong to the abstract type AbstractAlgebra.ResRing{T} or AbstractAlgebra.ResField{T} respectively. This enables one to write generic functions that can accept any AbstractAlgebra residue type.Note that both the generic residue ring type Generic.ResRing{T} and the abstract type it belongs to, AbstractAlgebra.ResRing{T} are both called ResRing, and  similarly for the residue field types. In each case, the  former is a (parameterised) concrete type for a residue ring over a given base ring whose elements have type T. The latter is an abstract type representing all residue ring types in  AbstractAlgebra.jl, whether generic or very specialised (e.g. supplied by a C library)."
},

{
    "location": "residue.html#Residue-ring-constructors-1",
    "page": "Generic residue rings",
    "title": "Residue ring constructors",
    "category": "section",
    "text": "In order to construct residues in AbstractAlgebra.jl, one must first construct the resiude ring itself. This is accomplished with one of the following constructors.ResidueRing(R::AbstractAlgebra.Ring, m::AbstractAlgebra.RingElem; cached::Bool = true)ResidueField(R::AbstractAlgebra.Ring, m::AbstractAlgebra.RingElem; cached::Bool = true)Given a base ring R and residue m contained in this ring, return the parent object of the residue ring R(m). By default the parent object S will depend only on R and m and will be cached. Setting the optional argument cached to false will prevent the parent object S from being cached.The ResidueField constructor does the same thing as the ResidueRing constructor, but the resulting object has type belonging to Field rather than Ring, so it can be used anywhere a field is expected in AbstractAlgebra.jl. No check is made for maximality of the ideal generated by m.Here are some examples of creating residue rings and making use of the resulting parent objects to coerce various elements into the residue ring.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\n\nf = S()\ng = S(123)\nh = S(BigInt(1234))\nk = S(x + 1)All of the examples here are generic residue rings, but specialised implementations of residue rings provided by external modules will also usually provide a ResidueRing constructor to allow creation of their residue rings."
},

{
    "location": "residue.html#Basic-ring-functionality-1",
    "page": "Generic residue rings",
    "title": "Basic ring functionality",
    "category": "section",
    "text": "Residue rings in AbstractAlgebra.jl implement the full Ring interface. Of course the entire Residue Ring interface is also implemented.We give some examples of such functionality.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\n\nf = S(x + 1)\n\nh = zero(S)\nk = one(S)\nisone(k) == true\niszero(f) == false\nm = modulus(S)\nU = base_ring(S)\nV = base_ring(f)\nT = parent(f)\nf == deepcopy(f)"
},

{
    "location": "residue.html#Residue-ring-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic residue rings",
    "title": "Residue ring functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality listed below is automatically provided by AbstractAlgebra.jl for any residue ring module that implements the full Residue Ring interface. This includes AbstractAlgebra.jl's own generic residue rings.But if a C library provides all the functionality documented in the Residue Ring interface, then all the functions described here will also be automatically supplied by AbstractAlgebra.jl for that residue ring type.Of course, modules are free to provide specific implementations of the functions described here, that override the generic implementation."
},

{
    "location": "residue.html#AbstractAlgebra.Generic.modulus-Tuple{AbstractAlgebra.ResElem}",
    "page": "Generic residue rings",
    "title": "AbstractAlgebra.Generic.modulus",
    "category": "Method",
    "text": "modulus(R::AbstractAlgebra.ResElem)\n\nReturn the modulus a of the residue ring S = R(a) that the supplied residue r belongs to.\n\n\n\n"
},

{
    "location": "residue.html#AbstractAlgebra.Generic.isunit-Tuple{AbstractAlgebra.ResElem}",
    "page": "Generic residue rings",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "Method",
    "text": "isunit(a::AbstractAlgebra.ResElem)\n\nReturn true if the supplied element a is invertible in the residue ring it belongs to, otherwise return false.\n\n\n\n"
},

{
    "location": "residue.html#Basic-functionality-1",
    "page": "Generic residue rings",
    "title": "Basic functionality",
    "category": "section",
    "text": "modulus(::AbstractAlgebra.ResElem)isunit(::AbstractAlgebra.ResElem)ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\n\nr = S(x + 1)\n\na = modulus(S)\nisunit(r) == true"
},

{
    "location": "residue.html#Base.inv-Tuple{AbstractAlgebra.ResElem}",
    "page": "Generic residue rings",
    "title": "Base.inv",
    "category": "Method",
    "text": "inv(a::AbstractAlgebra.ResElem)\n\nReturn the inverse of the element a in the residue ring. If an impossible inverse is encountered, an exception is raised.\n\n\n\n"
},

{
    "location": "residue.html#Inversion-1",
    "page": "Generic residue rings",
    "title": "Inversion",
    "category": "section",
    "text": "inv(::AbstractAlgebra.ResElem)ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R)\n\nf = S(x + 1)\n\ng = inv(f)"
},

{
    "location": "residue.html#Base.gcd-Union{Tuple{AbstractAlgebra.ResElem{T},AbstractAlgebra.ResElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic residue rings",
    "title": "Base.gcd",
    "category": "Method",
    "text": "gcd{T <: RingElement}(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T})\n\nReturn a greatest common divisor of a and b if one exists. This is done by taking the greatest common divisor of the data associated with the supplied residues and taking its greatest common divisor with the modulus.\n\n\n\n"
},

{
    "location": "residue.html#Greatest-common-divisor-1",
    "page": "Generic residue rings",
    "title": "Greatest common divisor",
    "category": "section",
    "text": "gcd{T <: RingElem}(::ResElem{T}, ::ResElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R)\n\nf = S(x + 1)\ng = S(x^2 + 2x + 1)\n\nh = gcd(f, g)"
},

{
    "location": "fields.html#",
    "page": "Field Interface",
    "title": "Field Interface",
    "category": "page",
    "text": ""
},

{
    "location": "fields.html#Field-Interface-1",
    "page": "Field Interface",
    "title": "Field Interface",
    "category": "section",
    "text": "AbstractAlgebra.jl generic code makes use of a standardised set of functions which it expects to be implemented for all fields. Here we document this interface. All libraries which want to make use of the generic capabilities of AbstractAlgebra.jl must supply all of the required functionality for their fields."
},

{
    "location": "fields.html#Types-1",
    "page": "Field Interface",
    "title": "Types",
    "category": "section",
    "text": "Most fields must supply two types:a type for the parent object (representing the field itself)\na type for elements of that fieldFor example, the generic fraction field type in AbstractAlgebra.jl provides two  types in generic/GenericTypes.jl: Generic.FracField{T} for the parent objects\nGeneric.Frac{T} for the actual fractionsThe parent type must belong to AbstractAlgebra.Field and the element type must belong to AbstractAlgebra.FieldElem. Of course, the types may belong to these abstract types transitively.For parameterised fields, we advise that the types of both the parent objects and element objects to be parameterised by the types of the elements of the base ring.There can be variations on this theme: e.g. in some areas of mathematics there is a notion of a coefficient domain, in which case it may make sense to parameterise all types by the type of elements of this coefficient domain. But note that this may have implications for the ad hoc operators one might like to explicitly implement."
},

{
    "location": "fields.html#Parent-object-caches-1",
    "page": "Field Interface",
    "title": "Parent object caches",
    "category": "section",
    "text": "In many cases, it is desirable to have only one object in the system to represent each field. This means that if the same field is constructed twice, elements of the two fields will be compatible as far as arithmetic is concerned.In order to facilitate this, global caches of fields are stored in AbstractAlgebra.jl, usually implemented using dictionaries. For example, the Generic.FracField parent objects are looked up in a dictionary FracDict to see if they have been previously defined.Whether these global caches are provided or not, depends on both mathematical and algorithmic considerations. E.g. in the case of number fields, it isn't desirable to identify all number fields with the same defining polynomial, as they may be considered with distinct embeddings into one another. In other cases, identifying whether two fields are the same may be prohibitively expensive. Generally, it may only make sense algorithmically to identify two fields if they were constructed from identical data.If a global cache is provided, it must be optionally possible to construct the parent objects without caching. This is done by passing a boolean value cached to the inner constructor of the parent object. See generic/GenericTypes.jl` for examples of how to construct and handle such caches."
},

{
    "location": "fields.html#Required-functions-for-all-fields-1",
    "page": "Field Interface",
    "title": "Required functions for all fields",
    "category": "section",
    "text": "In the following, we list all the functions that are required to be provided for fields in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.We give this interface for fictitious types MyParent for the type of the field parent object R and MyElem for the type of the elements of the field.Note that generic functions in AbstractAlgebra.jl may not rely on the existence of functions that are not documented here. If they do, those functions will only be available for fields that implement that additional functionality, and should be documented as such.In the first place, all fields are rings and therefore any field type must implement all of the Ring interface. The functionality below is in addition to this basic functionality."
},

{
    "location": "fields.html#Data-type-and-parent-object-methods-1",
    "page": "Field Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "characteristic(R::MyParent)Return the characteristic of the field."
},

{
    "location": "fields.html#Basic-manipulation-of-rings-and-elements-1",
    "page": "Field Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "isunit(f::MyElem)Return true if the given element is invertible, i.e. nonzero in the field."
},

{
    "location": "fields.html#Inversion-1",
    "page": "Field Interface",
    "title": "Inversion",
    "category": "section",
    "text": "inv(f::MyElem)Return the inverse of the given element in the field. If f = 0, an error is thrown."
},

{
    "location": "fraction_fields.html#",
    "page": "Fraction Field Interface",
    "title": "Fraction Field Interface",
    "category": "page",
    "text": ""
},

{
    "location": "fraction_fields.html#Fraction-Field-Interface-1",
    "page": "Fraction Field Interface",
    "title": "Fraction Field Interface",
    "category": "section",
    "text": "Fraction fields are supported in AbstractAlgebra.jl, at least for gcd domains. In addition to the standard Ring interface, some additional functions are required to be present for fraction fields."
},

{
    "location": "fraction_fields.html#Types-and-parents-1",
    "page": "Fraction Field Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra provides two abstract types for fraction fields and their elements:FracField{T} is the abstract type for fraction field parent types\nFracElem{T} is the abstract type for types of fractionsWe have that FracField{T} <: AbstractAlgebra.Field and  FracElem{T} <: AbstractAlgebra.FieldElem.Note that both abstract types are parameterised. The type T should usually be the type of elements of the base ring of the fraction field.Fraction fields should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Fraction fields should at least be distinguished based on their base ring.See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "fraction_fields.html#Required-functionality-for-fraction-fields-1",
    "page": "Fraction Field Interface",
    "title": "Required functionality for fraction fields",
    "category": "section",
    "text": "In addition to the required functionality for the Field interface the Fraction Field interface has the following required functions.We suppose that R is a fictitious base ring, and that S is the fraction field with  parent object S of type MyFracField{T}. We also assume the fractions in the field  have type MyFrac{T}, where T is the type of elements of the base ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem."
},

{
    "location": "fraction_fields.html#Constructors-1",
    "page": "Fraction Field Interface",
    "title": "Constructors",
    "category": "section",
    "text": "We provide the following constructors. Note that these constructors don't require construction of the parent object first. This is easier to achieve if the fraction element type doesn't contain a reference to the parent object, but merely contains a reference to the base ring. The parent object can then be constructed on demand.//(x::T, y::T) where T <: AbstractAlgebra.RingElemReturn the fraction xy.//(x::T, y::AbstractAlgebra.FracElem{T}) where T <: AbstractAlgebra.RingElemReturn xy where x is in the base ring of y.//(x::AbstractAlgebra.FracElem{T}, y::T) where T <: AbstractAlgebra.RingElemReturn xy where y is in the base ring of x.ExamplesR, x = PolynomialRing(ZZ, \"x\")\n\nf = (x^2 + x + 1)//(x^3 + 3x + 1)\ng = f//x\nh = x//f"
},

{
    "location": "fraction_fields.html#Basic-manipulation-of-fields-and-elements-1",
    "page": "Fraction Field Interface",
    "title": "Basic manipulation of fields and elements",
    "category": "section",
    "text": "numerator(d::MyFrac{T}) where T <: AbstractAlgebra.RingElemGiven a fraction d = ab return a, where ab is in lowest terms with respect to the canonical_unit and gcd functions on the base ring.denominator(d::MyFrac{T}) where T <: AbstractAlgebra.RingElemGiven a fraction d = ab return b, where ab is in lowest terms with respect to the canonical_unit and gcd functions on the base ring.ExamplesR, x = PolynomialRing(QQ, \"x\")\n\nf = (x^2 + x + 1)//(x^3 + 3x + 1)\n\nn = numerator(f)\nd = denominator(f)"
},

{
    "location": "fraction.html#",
    "page": "Generic fraction fields",
    "title": "Generic fraction fields",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "fraction.html#Generic-fraction-fields-1",
    "page": "Generic fraction fields",
    "title": "Generic fraction fields",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/generic/Fraction.jl for generic fraction fields over any gcd domain belonging to the AbstractAlgebra.jl abstract type hierarchy.As well as implementing the Fraction Field interface a number of generic algorithms are implemented for fraction fields. We describe this generic functionality below.All of the generic functionality is part of a submodule of AbstractAlgebra called Generic. This is exported by default so that it is not necessary to qualify the function names with the submodule name."
},

{
    "location": "fraction.html#Types-and-parent-objects-1",
    "page": "Generic fraction fields",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Fractions implemented using the AbstractAlgebra generics have type Generic.Frac{T} where T is the type of elements of the base ring. See the file src/generic/GenericTypes.jl for details.Parent objects of such fraction elements have type Generic.FracField{T}.The fraction element types belong to the abstract type AbstractAlgebra.FracElem{T} and the fraction field types belong to the abstract type AbstractAlgebra.FracRing{T}. This enables one to write generic functions that can accept any AbstractAlgebra fraction type.Note that both the generic fraction field type Generic.FracField{T} and the abstract type it belongs to, AbstractAlgebra.FracField{T} are both called FracField. The  former is a (parameterised) concrete type for a fraction field over a given base ring whose elements have type T. The latter is an abstract type representing all fraction field types in AbstractAlgebra.jl, whether generic or very specialised (e.g. supplied by a C library)."
},

{
    "location": "fraction.html#Fraction-field-constructors-1",
    "page": "Generic fraction fields",
    "title": "Fraction field constructors",
    "category": "section",
    "text": "In order to construct fractions in AbstractAlgebra.jl, one can first construct the fraction field itself. This is accomplished with the following constructor.FractionField(R::AbstractAlgebra.Ring; cached::Bool = true)Given a base ring R return the parent object of the fraction field of R. By default the parent object S will depend only on R and will be cached. Setting the optional argument cached to false will prevent the parent object S from being cached.Here are some examples of creating fraction fields and making use of the resulting parent objects to coerce various elements into the fraction field.ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS = FractionField(R)\n\nf = S()\ng = S(123)\nh = S(BigInt(1234))\nk = S(x + 1)All of the examples here are generic fraction fields, but specialised implementations of fraction fields provided by external modules will also usually provide a FractionField constructor to allow creation of the fraction fields they provide."
},

{
    "location": "fraction.html#Basic-field-functionality-1",
    "page": "Generic fraction fields",
    "title": "Basic field functionality",
    "category": "section",
    "text": "Fraction fields in AbstractAlgebra.jl implement the full Field interface. Of course the entire Fraction Field interface is also implemented.We give some examples of such functionality.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = FractionField(R)\n\nf = S(x + 1)\ng = (x^2 + x + 1)//(x^3 + 3x + 1)\n\nh = zero(S)\nk = one(S)\nisone(k) == true\niszero(f) == false\nm = characteristic(S)\nU = base_ring(S)\nV = base_ring(f)\nT = parent(f)\nr = deepcopy(f)\nn = numerator(g)\nd = denominator(g)"
},

{
    "location": "fraction.html#Fraction-field-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic fraction fields",
    "title": "Fraction field functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality listed below is automatically provided by AbstractAlgebra.jl for any fraction field module that implements the full Fraction Field interface. This includes AbstractAlgebra.jl's own generic fraction fields.But if a C library provides all the functionality documented in the Fraction Field interface, then all the functions described here will also be automatically supplied by AbstractAlgebra.jl for that fraction field type.Of course, modules are free to provide specific implementations of the functions described here, that override the generic implementation."
},

{
    "location": "fraction.html#Base.gcd-Union{Tuple{AbstractAlgebra.FracElem{T},AbstractAlgebra.FracElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic fraction fields",
    "title": "Base.gcd",
    "category": "Method",
    "text": "gcd{T <: RingElem}(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T})\n\nReturn a greatest common divisor of a and b if one exists. N.B: we define the GCD of ab and cd to be gcd(ad bc)bd, reduced to lowest terms. This requires the existence of a greatest common divisor function for the base ring.\n\n\n\n"
},

{
    "location": "fraction.html#Greatest-common-divisor-1",
    "page": "Generic fraction fields",
    "title": "Greatest common divisor",
    "category": "section",
    "text": "gcd{T <: RingElem}(::FracElem{T}, ::FracElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\n\nf = (x + 1)//(x^3 + 3x + 1)\ng = (x^2 + 2x + 1)//(x^2 + x + 1)\n\nh = gcd(f, g)"
},

{
    "location": "fraction.html#AbstractAlgebra.Generic.remove-Union{Tuple{AbstractAlgebra.FracElem{T},T}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic fraction fields",
    "title": "AbstractAlgebra.Generic.remove",
    "category": "Method",
    "text": "remove{T <: RingElem}(z::AbstractAlgebra.FracElem{T}, p::T)\n\nReturn the tuple n x such that z = p^nx where x has valuation 0 at p.\n\n\n\n"
},

{
    "location": "fraction.html#AbstractAlgebra.Generic.valuation-Union{Tuple{AbstractAlgebra.FracElem{T},T}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic fraction fields",
    "title": "AbstractAlgebra.Generic.valuation",
    "category": "Method",
    "text": "valuation{T <: RingElem}(z::AbstractAlgebra.FracElem{T}, p::T)\n\nReturn the valuation of z at p.\n\n\n\n"
},

{
    "location": "fraction.html#Remove-and-valuation-1",
    "page": "Generic fraction fields",
    "title": "Remove and valuation",
    "category": "section",
    "text": "When working over a Euclidean domain, it is convenient to extend valuations to the fraction field. To facilitate this, we define the following functions.remove{T <: RingElem}(::FracElem{T}, ::T)valuation{T <: RingElem}(::FracElem{T}, ::T)ExamplesR, x = PolynomialRing(ZZ, \"x\")\n\nf = (x + 1)//(x^3 + 3x + 1)\ng = (x^2 + 1)//(x^2 + x + 1)\n\nv, q = remove(f^3*g, x + 1)\nv = valuation(f^3*g, x + 1)"
},

{
    "location": "rational.html#",
    "page": "Rational field",
    "title": "Rational field",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "rational.html#Rational-field-1",
    "page": "Rational field",
    "title": "Rational field",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/julia/Rational.jl for making Julia Rational{BigInt}s conform to the AbstractAlgebra.jl Field interface.In addition to providing a parent object QQ for Julia Rational{BigInt}s, we implement any additional functionality required by AbstractAlgebra.jl.Because Rational{BigInt} cannot be directly included in the AbstractAlgebra.jl abstract type hierarchy, we achieve integration of Julia Rational{BigInt}s by introducing a type union, called FieldElement, which is a union of AbstractAlgebra.FieldElem and a number of Julia types, including Rational{BigInt}. Everywhere that FieldElem is notionally used in AbstractAlgebra.jl, we are in fact using FieldElement, with additional care being taken to avoid ambiguities.The details of how this is done are technical, and we refer the reader to the implementation for details. For most intents and purposes, one can think of the Julia Rational{BigInt} type as belonging to AbstractAlgebra.FieldElem.One other technicality is that Julia defines certain functions for Rational{BigInt}, such as sqrt and exp differently to what AbstractAlgebra.jl requires. To get around this, we redefine these functions internally to AbstractAlgebra.jl, without redefining them for users of AbstractAlgebra.jl. This allows the internals of AbstractAlgebra.jl to function correctly, without broadcasting pirate definitions of already defined Julia functions to the world.To access the internal definitions, one can use AbstractAlgebra.sqrt and AbstractAlgebra.exp, etc."
},

{
    "location": "rational.html#Types-and-parent-objects-1",
    "page": "Rational field",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Rationals have type Rational{BigInt}, as in Julia itself. We simply supplement the functionality for this type as required for computer algebra.The parent objects of such integers has type Rationals{BigInt}.For convenience, we also make Rational{Int} a part of the AbstractAlgebra.jl type hierarchy and its parent object (accessible as qq) has type Rationals{Int}. But we caution that this type is not particularly useful as a model of the rationals and may not function as expected within AbstractAlgebra.jl."
},

{
    "location": "rational.html#Rational-constructors-1",
    "page": "Rational field",
    "title": "Rational constructors",
    "category": "section",
    "text": "In order to construct rationals in AbstractAlgebra.jl, one can first construct the rational field itself. This is accomplished using either of the following constructors.FractionField(R::Integers{BigInt})Rationals{BigInt}()This gives the unique object of type Rationals{BigInt} representing the field of rationals in AbstractAlgebra.jl.In practice, one simply uses QQ which is assigned to be the return value of the above constructor. There is no need to call the constructor in practice.Here are some examples of creating the rational field and making use of the resulting parent object to coerce various elements into the field.Examplesf = QQ()\ng = QQ(123)\nh = QQ(BigInt(1234))\nk = QQ(BigInt(12), BigInt(7))\n\nQQ == FractionField(ZZ)"
},

{
    "location": "rational.html#Basic-field-functionality-1",
    "page": "Rational field",
    "title": "Basic field functionality",
    "category": "section",
    "text": "The rational field in AbstractAlgebra.jl implements the full Field and Fraction Field interfaces.We give some examples of such functionality.Examplesf = QQ(12, 7)\n\nh = zero(QQ)\nk = one(QQ)\nisone(k) == true\niszero(f) == false\nU = base_ring(QQ)\nV = base_ring(f)\nT = parent(f)\nf == deepcopy(f)\ng = f + 12\nr = ZZ(12)//ZZ(7)\nn = numerator(r)"
},

{
    "location": "rational.html#Rational-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Rational field",
    "title": "Rational functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality below supplements that provided by Julia itself for its Rational{BigInt} type."
},

{
    "location": "rational.html#AbstractAlgebra.sqrt-Tuple{Rational{BigInt}}",
    "page": "Rational field",
    "title": "AbstractAlgebra.sqrt",
    "category": "Method",
    "text": "sqrt{T <: Integer}(a::Rational{T})\n\nReturn the square root of a if it is the square of a rational, otherwise throw an error.\n\n\n\n"
},

{
    "location": "rational.html#AbstractAlgebra.exp-Tuple{Rational{BigInt}}",
    "page": "Rational field",
    "title": "AbstractAlgebra.exp",
    "category": "Method",
    "text": "exp{T <: Integer}(a::Rational{T})\n\nReturn 1 if a = 0, otherwise throw an exception.\n\n\n\n"
},

{
    "location": "rational.html#Square-root-1",
    "page": "Rational field",
    "title": "Square root",
    "category": "section",
    "text": "AbstractAlgebra.sqrt(a::Rational{BigInt})AbstractAlgebra.exp(a::Rational{BigInt})Examplesd = AbstractAlgebra.sqrt(ZZ(36)//ZZ(25))\nm = AbstractAlgebra.exp(ZZ(0)//ZZ(1))"
},

{
    "location": "finfield.html#",
    "page": "Finite fields",
    "title": "Finite fields",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "finfield.html#Finite-fields-1",
    "page": "Finite fields",
    "title": "Finite fields",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/julia/GF.jl for finite fields. The module is a naive implementation that supports only fields of degree 1 (prime fields). They are modelled as mathbbZpmathbbZ for p a prime."
},

{
    "location": "finfield.html#Types-and-parent-objects-1",
    "page": "Finite fields",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Finite fields have type GFField{T} where T is either Int or BigInt.Elements of such a finite field have type gfelem{T}."
},

{
    "location": "finfield.html#AbstractAlgebra.GF-Union{Tuple{T}, Tuple{T}} where T<:Integer",
    "page": "Finite fields",
    "title": "AbstractAlgebra.GF",
    "category": "Method",
    "text": "GF{T <: Integer}(p::T)\n\nReturn the finite field mathbbF_p, where p is a prime. The integer p is not checked for primality, but the behaviour of the resulting object is undefined if p is composite.\n\n\n\n"
},

{
    "location": "finfield.html#Finite-field-constructors-1",
    "page": "Finite fields",
    "title": "Finite field constructors",
    "category": "section",
    "text": "In order to construct finite fields in AbstractAlgebra.jl, one must first construct the field itself. This is accomplished with the following constructors.GF(p::T) where T <: IntegerHere are some examples of creating a finite field and making use of the resulting parent object to coerce various elements into the field.ExamplesF = GF(13)\n\ng = F(3)\nh = F(g)"
},

{
    "location": "finfield.html#Basic-field-functionality-1",
    "page": "Finite fields",
    "title": "Basic field functionality",
    "category": "section",
    "text": "The finite field module in AbstractAlgebra.jl implements the full Field interface.We give some examples of such functionality.ExamplesF = GF(13)\n\nh = zero(F)\nk = one(F)\nisone(k) == true\niszero(f) == false\nU = base_ring(F)\nV = base_ring(h)\nT = parent(h)\nh == deepcopy(h)\nh = h + 2\nm = inv(k)"
},

{
    "location": "finfield.html#AbstractAlgebra.Generic.gen-Union{Tuple{AbstractAlgebra.GFField{T}}, Tuple{T}} where T<:Integer",
    "page": "Finite fields",
    "title": "AbstractAlgebra.Generic.gen",
    "category": "Method",
    "text": "gen{T <: Integer}(a::GFField{T})\n\nReturn a generator of the field. Currently this returns 1.\n\n\n\n"
},

{
    "location": "finfield.html#AbstractAlgebra.Generic.order-Tuple{AbstractAlgebra.GFField}",
    "page": "Finite fields",
    "title": "AbstractAlgebra.Generic.order",
    "category": "Method",
    "text": "order(R::GFField)\n\nReturn the order, i.e. the number of element in, the given finite field.\n\n\n\n"
},

{
    "location": "finfield.html#AbstractAlgebra.Generic.degree-Tuple{AbstractAlgebra.GFField}",
    "page": "Finite fields",
    "title": "AbstractAlgebra.Generic.degree",
    "category": "Method",
    "text": "degree(R::GFField)\n\nReturn the degree of the given finite field.\n\n\n\n"
},

{
    "location": "finfield.html#Basic-manipulation-of-fields-and-elements-1",
    "page": "Finite fields",
    "title": "Basic manipulation of fields and elements",
    "category": "section",
    "text": "gen{T <: Integer}(F::GFField{T})order(F::GFField)degree(F::GFField)ExamplesF = GF(13)\n\nd = degree(F)\nn = order(F)\ng = gen(F)"
},

{
    "location": "numberfield.html#",
    "page": "Number fields",
    "title": "Number fields",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "numberfield.html#Number-fields-1",
    "page": "Number fields",
    "title": "Number fields",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a very naive implementation of number fields. This allows arithmetic in algebraic number fields, which are currently modeled as mathbbZx modulo an irreducible polynomial, i.e. as a residue field.In fact, the definition of the number field constructor is currently given in src/generic/ResidueField.jl and no type is defined for a number field. The definition mainly exists for testing purposes. It may later be replaced by a more standard implementation. For a more fully fleshed out number field implementation (based on a very high performance C library), see Nemo.jl."
},

{
    "location": "numberfield.html#Number-field-constructors-1",
    "page": "Number fields",
    "title": "Number field constructors",
    "category": "section",
    "text": "In order to construct number fields in AbstractAlgebra.jl, one must first construct the field itself. This is accomplished with the following constructor.NumberField(f::AbstractAlgebra.Generic.Poly{Rational{BigInt}}, s::AbstractString, t = \"\\$\"; cached = true)Given an irreducible defining polynomial f in mathbbZx, return a tuple (K x) consisting of the number field defined by that polynomial and a generator. The string fields are currently ignored, but are reserved for future use.Currently the generator of the number field prints the same way as the variable in mathbbZx.ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(R, x^3 + 3x + 1, \"a\")\n\nf = a^2 + 2a + 7"
},

{
    "location": "numberfield.html#Basic-field-functionality-1",
    "page": "Number fields",
    "title": "Basic field functionality",
    "category": "section",
    "text": "The number field module in AbstractAlgebra.jl implements the full Field and ResidueRing interfaces."
},

{
    "location": "matrix_spaces.html#",
    "page": "Matrix Interface",
    "title": "Matrix Interface",
    "category": "page",
    "text": ""
},

{
    "location": "matrix_spaces.html#Matrix-Interface-1",
    "page": "Matrix Interface",
    "title": "Matrix Interface",
    "category": "section",
    "text": "Generic matrices are supported in AbstractAlgebra.jl. As the space of mtimes n matrices over a commutative ring is not itself a commutative ring, not all of the Ring interface needs to be implemented for matrices in general.In particular, the following functions do not need to be implemented: isdomain_type, needs_parentheses, isnegative, show_minus_one and divexact. The canonical_unit function should be implemented, but simply needs to return the corresponding value for entry 1 1 (the function is never called on empty matrices).Note that AbstractAlgebra.jl matrices are not the same as Julia matrices. We store a base ring in our matrix and matrices are row major instead of column major in order to support the numerous large C libraries that use this convention.All AbstractAlgebra.jl matrices are assumed to be mutable. This is usually critical to performance."
},

{
    "location": "matrix_spaces.html#Types-and-parents-1",
    "page": "Matrix Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra provides two abstract types for matrix spaces and their elements:MatSpace{T} is the abstract type for matrix space parent types\nMatElem{T} is the abstract type for matrix typesNote that both abstract types are parameterised. The type T should usually be the type of elements of the matrices.Matrix spaces should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Matrix spaces should at least be distinguished based on their base (coefficient) ring and the dimensions of the matrices in the space.See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "matrix_spaces.html#Required-functionality-for-matrices-1",
    "page": "Matrix Interface",
    "title": "Required functionality for matrices",
    "category": "section",
    "text": "In addition to the required (relevant) functionality for the Ring interface (see above), the following functionality is required for the Matrix interface.We suppose that R is a fictitious base ring (coefficient ring) and that S is a space of mtimes n matrices over R with parent object S of type MyMatSpace{T}. We also assume the matrices in the space have type MyMat{T}, where T is the type of elements of the base (element) ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem."
},

{
    "location": "matrix_spaces.html#Constructors-1",
    "page": "Matrix Interface",
    "title": "Constructors",
    "category": "section",
    "text": "In addition to the standard constructors, the following constructors, taking an array of elements, must be available.(S::MyMatSpace{T})(A::Array{T, 2}) where T <: AbstractAlgebra.RingElemCreate the matrix in the given space whose (i j) entry is given by A[i, j].(S::MyMatSpace{T})(A::Array{S, 2}) where {S <: AbstractAlgebra.RingElem, T <: AbstractAlgebra.RingElem}Create the matrix in the given space whose (i j) entry is given by A[i, j], where S is the type of elements that can be coerced into the base ring of the matrix.(S::MyMatSpace{T})(A::Array{S, 1}) where {S <: AbstractAlgebra.RingElem, T <: AbstractAlgebra.RingElem}Create the matrix in the given space of matrices (with dimensions mtimes n say), whose (i j) entry is given by A[i*(n - 1) + j] and where S is the type of elements that can be coerced into the base ring of the matrix.ExamplesS = MatrixSpace(QQ, 2, 3)\n\nM1 = S(Rational{BigInt}[2 3 1; 1 0 4])\nM2 = S(BigInt[2 3 1; 1 0 4])\nM3 = S(BigInt[2, 3, 1, 1, 0, 4])It is also possible to create matrices directly, without first creating the corresponding matrix space (the inner constructor should be called directly). Note that to support this, matrix space parent objects don't contain a reference to their parent. Instead, parents are constructed on-the-fly if requested.matrix(R::Ring, arr::Array{T, 2}) where T <: AbstractAlgebra.RingElemGiven an mtimes n Julia matrix of entries, construct the corresponding AbstractAlgebra.jl matrix over the given ring R, assuming all the entries can be coerced into R.matrix(R::Ring, r::Int, c::Int, A::Array{T, 1}) where T <: AbstractAlgebra.RingElemConstruct the given rtimes c AbstractAlgebra.jl matrix over the ring R whose (i j) entry is given by A[c*(i - 1) + j], assuming that all the entries can be coerced into R.zero_matrix(R::Ring, r::Int, c::Int)Construct the rtimes c AbstractAlgebra.jl zero matrix over the ring R.identity_matrix(R::Ring, n::Int)Construct the ntimes n AbstractAlgebra.jl identity matrix over the ring R.similar(x::MyMat{T}) where T <: AbstractAlgebra.RingElemConstruct the zero matrix with the same dimensions and base ring as the given matrix.similar(x::MyMat{T}, r::Int, c::Int) where T <: AbstractAlgebra.RingElemConstruct the rtimes c zero matrix with the same base ring as the given matrix.ExamplesM = matrix(ZZ, BigInt[3 1 2; 2 0 1])\nN = matrix(ZZ, 3, 2, BigInt[3, 1, 2, 2, 0, 1])\nP = zero_matrix(ZZ, 3, 2)\nQ = identity_matrix(ZZ, 4)\nC = similar(P)\nD = similar(Q, 4, 5)"
},

{
    "location": "matrix_spaces.html#Basic-manipulation-of-matrices-1",
    "page": "Matrix Interface",
    "title": "Basic manipulation of matrices",
    "category": "section",
    "text": "rows(f::MyMat{T}) where T <: AbstractAlgebra.RingElemReturn the number of rows of the given matrix.cols(f::MyMat{T}) where T <: AbstractAlgebra.RingElemReturns the number of columns of the given matrix.getindex(M::MyMat{T}, r::Int, c::Int) where T <: AbstractAlgebra.RingElemReturn the (i j)-th entry of the matrix M.setindex!(M::MyMat{T}, d::T, r::Int, c::Int) where T <: AbstractAlgebra.RingElemSet the (i j)-th entry of the matrix M to d, which is assumed to be in the base ring of the matrix. The matrix must have such an entry and the matrix is mutated in place and not returned from the function.ExamplesM = matrix(ZZ, BigInt[2 3 0; 1 1 1])\n\nm = rows(M)\nn = cols(M)\nM[1, 2] = BigInt(4)\nc = M[1, 1]"
},

{
    "location": "matrix_spaces.html#Transpose-1",
    "page": "Matrix Interface",
    "title": "Transpose",
    "category": "section",
    "text": "transpose(::MyMat{T}) where T <: AbstractAlgebra.RingElemReturn the transpose of the given matrix.The standard Julia tick notation can also be used for transposing a matrix.ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\n\nB = transpose(A)\nC = A'"
},

{
    "location": "matrix_spaces.html#Optional-functionality-for-matrices-1",
    "page": "Matrix Interface",
    "title": "Optional functionality for matrices",
    "category": "section",
    "text": "Especially when wrapping C libraries, some functions are best implemented directly, rather than relying on the generic functionality. The following are all provided by the AbstractAlgebra.jl generic code, but can optionally be implemented directly for performance reasons."
},

{
    "location": "matrix_spaces.html#Optional-constructors-1",
    "page": "Matrix Interface",
    "title": "Optional constructors",
    "category": "section",
    "text": "eye(M::MyMat{T}) where T <: AbstractAlgebra.RingElemConstruct the identity matrix with the same dimensions and base ring as the given matrix.eye(M::MyMat{T}, n::Int) where T <: AbstractAlgebra.RingElemConstruct the ntimes n identity matrix with the same base ring as the given matrix.ExamplesM = matrix(ZZ, BigInt[1 2 3; 4 5 6])\n\nN = eye(M)\nP = eye(M, 2)"
},

{
    "location": "matrix_spaces.html#Optional-submatrices-1",
    "page": "Matrix Interface",
    "title": "Optional submatrices",
    "category": "section",
    "text": "sub(M::MyMat{T}, rows::UnitRange{Int}, cols::UnitRange{Int}) where T <: AbstractAlgebra.RingElemReturn a new matrix with the same entries as the submatrix with the given range of rows and columns.ExamplesM = matrix(ZZ, BigInt[1 2 3; 2 3 4; 3 4 5])\n\nN1 = M[1:2, :]\nN2 = M[:, :]\nN3 = M[2:3, 2:3]"
},

{
    "location": "matrix_spaces.html#Optional-row-swapping-1",
    "page": "Matrix Interface",
    "title": "Optional row swapping",
    "category": "section",
    "text": "swap_rows!(M::MyMat{T}, i::Int, j::Int) where T <: AbstractAlgebra.RingElemSwap the rows of M in place. The function does not return the mutated matrix (since matrices are assumed to be mutable in AbstractAlgebra.jl).ExamplesM = identity_matrix(ZZ, 3)\n\nswap_rows!(M, 1, 2)"
},

{
    "location": "matrix_spaces.html#Optional-concatenation-1",
    "page": "Matrix Interface",
    "title": "Optional concatenation",
    "category": "section",
    "text": "hcat(M::MyMat{T}, N::MyMat{T}) where T <: AbstractAlgebra.RingElemReturn the horizontal concatenation of M and N. It is assumed that the number of rows of M and N are the same.vcat(M::MyMat{T}, N::MyMat{T}) where T <: AbstractAlgebra.RingElemReturn the vertical concatenation of M and N. It is assumed that the number of columns of M and N are the same.ExamplesM = matrix(ZZ, BigInt[1 2 3; 2 3 4; 3 4 5])\nN = matrix(ZZ, BigInt[1 0 1; 0 1 0; 1 0 1])\n\nP = hcat(M, N)\nQ = vcat(M, N)"
},

{
    "location": "matrix.html#",
    "page": "Generic matrices",
    "title": "Generic matrices",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "matrix.html#Generic-matrices-1",
    "page": "Generic matrices",
    "title": "Generic matrices",
    "category": "section",
    "text": "AbstractAlgebra.jl allows the creation of dense matrices over any computable commutative ring R. Generic matrices over a commutative ring are implemented in src/generic/Matrix.jl.As well as implementing the entire Matrix interface, including the optional functionality, there are many additional generic algorithms implemented for matrix spaces. We describe this functionality below.All of this generic functionality is part of the Generic submodule of AbstractAlgebra.jl. This is exported by default, so it is not necessary to qualify names of functions."
},

{
    "location": "matrix.html#Types-and-parent-objects-1",
    "page": "Generic matrices",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Generic matrices in AbstractAlgebra.jl have type Generic.Mat{T} where T is the type of elements of the matrix. Internally, generic matrices are implemented using an object wrapping a Julia two dimensional array, though they are not themselves Julia arrays. See the file src/generic/GenericTypes.jl for details.Parents of generic matrices (matrix spaces) have type Generic.MatSpace{T}.The generic matrix types belong to the abstract type AbstractAlgebra.MatElem{T} and the matrix space parent types belong to AbstractAlgebra.MatSpace{T}. Note that both the concrete type of a matrix space parent object and the abstract class it belongs to have the name MatElem, therefore disambiguation is required to specify which is intended.The dimensions and base ring R of a generic matrix are stored in its parent object, however to allow creation of matrices without first creating the matrix space parent, generic matrices in Julia do not contain a reference to their parent. They contain the row and column numbers and the base ring on a per matrix basis. The parent object can then be reconstructed from this data on demand."
},

{
    "location": "matrix.html#Matrix-space-constructors-1",
    "page": "Generic matrices",
    "title": "Matrix space constructors",
    "category": "section",
    "text": "A matrix space in AbstractAlgebra.jl represents a collection of all matrices with given dimensions and base ring.In order to construct matrices in AbstractAlgebra.jl, one can first constructs the matrix space itself. This is accomplished with the following constructor.MatrixSpace(R::Ring, rows::Int, cols::Int; cache::Bool=true)Construct the space of matrices with the given number of rows and columns over the given base ring. By default such matrix spaces are cached based on the base ring and numbers of rows and columns. If the optional named parameter cached is set to false, no caching occurs.Here are some examples of creating matrix spaces and making use of the resulting parent objects to coerce various elements into the matrix space.ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S()\nB = S(12)\nC = S(R(11))We also allow matrices over a given base ring to be constructed directly (see the Matrix interface)."
},

{
    "location": "matrix.html#Matrix-element-constructors-1",
    "page": "Generic matrices",
    "title": "Matrix element constructors",
    "category": "section",
    "text": "In addition to coercing elements into a matrix space as above, we provide the following functions for constructing explicit matrices.Also see the Matrix interface for a list of other ways to create matrices.R[a b c...;...]Create the matrix over the base ring R consisting of the given rows (separated by semicolons). Each entry is coerced into R  automatically. Note that parentheses may be placed around individual entries if the lists would otherwise be ambiguous, e.g.  R[1 2; 2 (-3)].Beware that this syntax does not support the creation of column vectors. See the notation below for creating those.R[a b c...]Create the row vector with entries in R consisting of the given entries (separated by spaces). Each entry is coerced into R automatically. Note that parentheses may be placed around individual entries if the list would otherwise be ambiguous, e.g. R[1 2 (-3)].R[a b c...]'Create the column vector with entries in R consisting of the given entries (separated by spaces). Each entry is coerced into R automatically. Observe the dash that is used to transpose the row vector notation (for free) to turn it into a column vector. Note that parentheses may be placed around individual entries if the list would otherwise be ambiguous, e.g. R[1 2 (-3)]'.ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nM = R[t + 1 1; t^2 0]\nN = R[t + 1 2 t]\nP = R[1 2 t]'"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.sub-Tuple{AbstractAlgebra.MatElem,Int64,Int64,Int64,Int64}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.sub",
    "category": "Method",
    "text": "sub(M::AbstractAlgebra.MatElem, r1::Int, c1::Int, r2::Int, c2::Int)\n\nReturn a copy of the submatrix of M from (r1 c1) to (r2 c2) inclusive. Note that is the copy is modified, the original matrix is not.\n\n\n\n"
},

{
    "location": "matrix.html#Submatrices-1",
    "page": "Generic matrices",
    "title": "Submatrices",
    "category": "section",
    "text": "In addition to the functionality described in the Matrix interface for taking submatrices of a matrix, the following function variant is also available.sub(::MatElem, ::Int, ::Int, ::Int, ::Int)ExamplesM = ZZ[1 2 3; 2 3 4]\n\nN = sub(M, 1, 1, 2, 2)"
},

{
    "location": "matrix.html#Matrix-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic matrices",
    "title": "Matrix functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": ""
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.rows-Tuple{AbstractAlgebra.MatElem}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.rows",
    "category": "Method",
    "text": "rows(a::AbstractAlgebra.MatElem)\n\nReturn the number of rows of the given matrix.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.cols-Tuple{AbstractAlgebra.MatElem}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.cols",
    "category": "Method",
    "text": "cols(a::AbstractAlgebra.MatElem)\n\nReturn the number of columns of the given matrix.\n\n\n\n"
},

{
    "location": "matrix.html#Basic-matrix-functionality-1",
    "page": "Generic matrices",
    "title": "Basic matrix functionality",
    "category": "section",
    "text": "As well as the Ring and Matrix interfaces, the following functions are provided to manipulate matrices and to set and retrieve entries and other basic data associated with the matrices.rows(::MatElem)cols(::MatElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\nB = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])\n\nr = rows(B)\nc = cols(B)\nM = A + B\nN = 2 + A\nM1 = deepcopy(A)\nA != B\nisone(one(S)) == true\nV = A[1:2, :]\nW = A^3\nZ = divexact(2*A, 2)"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.powers-Tuple{AbstractAlgebra.MatElem,Int64}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.powers",
    "category": "Method",
    "text": "powers{T <: RingElement}(a::AbstractAlgebra.MatElem{T}, d::Int)\n\nReturn an array of matrices M wher Mi + 1 = a^i for i = 0d\n\n\n\n"
},

{
    "location": "matrix.html#Powering-1",
    "page": "Generic matrices",
    "title": "Powering",
    "category": "section",
    "text": "powers(::MatElem, ::Int)ExamplesM = ZZ[1 2 3; 2 3 4; 4 5 5]\n\nA = powers(M, 4)"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.gram-Tuple{AbstractAlgebra.MatElem}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.gram",
    "category": "Method",
    "text": "gram(x::AbstractAlgebra.MatElem)\n\nReturn the Gram matrix of x, i.e. if x is an rtimes c matrix return the rtimes r matrix whose entries i j are the dot products of the i-th and j-th rows, respectively.\n\n\n\n"
},

{
    "location": "matrix.html#Gram-matrix-1",
    "page": "Generic matrices",
    "title": "Gram matrix",
    "category": "section",
    "text": "gram(::MatElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\n\nB = gram(A)"
},

{
    "location": "matrix.html#Base.LinAlg.trace-Tuple{AbstractAlgebra.MatElem}",
    "page": "Generic matrices",
    "title": "Base.LinAlg.trace",
    "category": "Method",
    "text": "trace(x::AbstractAlgebra.MatElem)\n\nReturn the trace of the matrix a, i.e. the sum of the diagonal elements. We require the matrix to be square.\n\n\n\n"
},

{
    "location": "matrix.html#Trace-1",
    "page": "Generic matrices",
    "title": "Trace",
    "category": "section",
    "text": "trace(::MatElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\n\nb = trace(A)"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.content-Tuple{AbstractAlgebra.MatElem}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.content",
    "category": "Method",
    "text": "content(x::AbstractAlgebra.MatElem)\n\nReturn the content of the matrix a, i.e. the greatest common divisor of all its entries, assuming it exists.\n\n\n\n"
},

{
    "location": "matrix.html#Content-1",
    "page": "Generic matrices",
    "title": "Content",
    "category": "section",
    "text": "content(::MatElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\n\nb = content(A)"
},

{
    "location": "matrix.html#Base.:*-Tuple{AbstractAlgebra.perm,AbstractAlgebra.MatElem}",
    "page": "Generic matrices",
    "title": "Base.:*",
    "category": "Method",
    "text": "*(x, y...)\n\nMultiplication operator. x*y*z*... calls this function with all arguments, i.e. *(x, y, z, ...).\n\n\n\n"
},

{
    "location": "matrix.html#Permutation-1",
    "page": "Generic matrices",
    "title": "Permutation",
    "category": "section",
    "text": "*(::perm, ::MatElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\nG = PermGroup(3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\nP = G([1, 3, 2])\n\nB = P*A"
},

{
    "location": "matrix.html#Base.LinAlg.lufact-Union{Tuple{AbstractAlgebra.MatElem{T},AbstractAlgebra.PermGroup}, Tuple{T}} where T<:AbstractAlgebra.FieldElem",
    "page": "Generic matrices",
    "title": "Base.LinAlg.lufact",
    "category": "Method",
    "text": "lufact{T <: FieldElement}(A::AbstractAlgebra.MatElem{T}, P = PermGroup(rows(A)))\n\nReturn a tuple r p L U consisting of the rank of A, a permutation p of A belonging to P, a lower triangular matrix L and an upper triangular matrix U such that p(A) = LU, where p(A) stands for the matrix whose rows are the given permutation p of the rows of A.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.fflu-Union{Tuple{AbstractAlgebra.MatElem{T},AbstractAlgebra.PermGroup}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.fflu",
    "category": "Method",
    "text": "fflu{T <: RingElement}(A::AbstractAlgebra.MatElem{T}, P = PermGroup(rows(A)))\n\nReturn a tuple r d p L U consisting of the rank of A, a denominator d, a permutation p of A belonging to P, a lower triangular matrix L and an upper triangular matrix U such that p(A) = LD^1U, where p(A) stands for the matrix whose rows are the given permutation p of the rows of A and such that D is the diagonal matrix diag(p_1 p_1p_2 ldots p_n-2p_n-1 p_n-1 where the p_i are the inverses of the diagonal entries of U. The denominator d is set to pm mboxdet(S) where S is an appropriate submatrix of A (S = A if A is square) and the sign is decided by the parity of the permutation.\n\n\n\n"
},

{
    "location": "matrix.html#LU-factorisation-1",
    "page": "Generic matrices",
    "title": "LU factorisation",
    "category": "section",
    "text": "lufact{T <: FieldElem}(::MatElem{T}, ::PermGroup)fflu{T <: RingElem}(::MatElem{T}, ::PermGroup)ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\n\nA = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 - 2 a - 1 2a])\n\nr, P, L, U = lufact(A)\nr, d, P, L, U = fflu(A)"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.rref-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.rref",
    "category": "Method",
    "text": "rref{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nReturns a tuple (r d A) consisting of the rank r of M and a denominator d in the base ring of M and a matrix A such that Ad is the reduced row echelon form of M. Note that the denominator is not usually minimal.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.rref-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.FieldElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.rref",
    "category": "Method",
    "text": "rref{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nReturns a tuple (r d A) consisting of the rank r of M and a denominator d in the base ring of M and a matrix A such that Ad is the reduced row echelon form of M. Note that the denominator is not usually minimal.\n\n\n\nrref{T <: FieldElement}(M::AbstractAlgebra.MatElem{T})\n\nReturns a tuple (r A) consisting of the rank r of M and a reduced row echelon form A of M.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.isrref-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.isrref",
    "category": "Method",
    "text": "isrref{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nReturn true if M is in reduced row echelon form, otherwise return false.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.isrref-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.FieldElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.isrref",
    "category": "Method",
    "text": "isrref{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nReturn true if M is in reduced row echelon form, otherwise return false.\n\n\n\nisrref(M::AbstractAlgebra.MatElem{T}) where {T <: FieldElement}\n\nReturn true if M is in reduced row echelon form, otherwise return false.\n\n\n\n"
},

{
    "location": "matrix.html#Reduced-row-echelon-form-1",
    "page": "Generic matrices",
    "title": "Reduced row-echelon form",
    "category": "section",
    "text": "rref{T <: RingElem}(::MatElem{T})\nrref{T <: FieldElem}(::MatElem{T})isrref{T <: RingElem}(::MatElem{T})\nisrref{T <: FieldElem}(::MatElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\n   \nM = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])\n   \nr, A = rref(M)\nisrref(A)\n\nR, x = PolynomialRing(ZZ, \"x\")\nS = MatrixSpace(R, 3, 3)\n\nM = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])\n\nr, d, A = rref(M)\nisrref(A)"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.hnf-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.hnf",
    "category": "Method",
    "text": "hnf{T <: RingElement}(A::Mat{T}) -> Mat{T}\n\nReturn the upper right row Hermite normal form of A.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.hnf_with_trafo-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.hnf_with_trafo",
    "category": "Method",
    "text": "hnf{T <: RingElement}(A::Mat{T}) -> Mat{T}, Mat{T}\n\nReturn the upper right row Hermite normal form H of A together with invertible matrix U such that UA = H.\n\n\n\n"
},

{
    "location": "matrix.html#Hermite-normal-form-1",
    "page": "Generic matrices",
    "title": "Hermite normal form",
    "category": "section",
    "text": "hnf{T <: RingElem}(::MatElem{T})\nhnf_with_trafo{T <: RingElem}(::MatElem{T})"
},

{
    "location": "matrix.html#Base.LinAlg.det-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "Base.LinAlg.det",
    "category": "Method",
    "text": "det(M)\n\nMatrix determinant.\n\nExample\n\njulia> M = [1 0; 2 2]\n2×2 Array{Int64,2}:\n 1  0\n 2  2\n\njulia> det(M)\n2.0\n\n\n\ndet{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nReturn the determinant of the matrix M. We assume M is square.\n\n\n\n"
},

{
    "location": "matrix.html#Base.LinAlg.det-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.FieldElem",
    "page": "Generic matrices",
    "title": "Base.LinAlg.det",
    "category": "Method",
    "text": "det(M)\n\nMatrix determinant.\n\nExample\n\njulia> M = [1 0; 2 2]\n2×2 Array{Int64,2}:\n 1  0\n 2  2\n\njulia> det(M)\n2.0\n\n\n\ndet{T <: FieldElement}(M::AbstractAlgebra.MatElem{T})\n\nReturn the determinant of the matrix M. We assume M is square.\n\n\n\ndet{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nReturn the determinant of the matrix M. We assume M is square.\n\n\n\n"
},

{
    "location": "matrix.html#Determinant-1",
    "page": "Generic matrices",
    "title": "Determinant",
    "category": "section",
    "text": "det{T <: RingElem}(::MatElem{T})\ndet{T <: FieldElem}(::MatElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\n   \nA = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])\n\nd = det(A)"
},

{
    "location": "matrix.html#Base.LinAlg.rank-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "Base.LinAlg.rank",
    "category": "Method",
    "text": "rank{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nReturn the rank of the matrix M.\n\n\n\n"
},

{
    "location": "matrix.html#Base.LinAlg.rank-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.FieldElem",
    "page": "Generic matrices",
    "title": "Base.LinAlg.rank",
    "category": "Method",
    "text": "rank{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nReturn the rank of the matrix M.\n\n\n\nrank{T <: FieldElement}(M::AbstractAlgebra.MatElem{T})\n\nReturn the rank of the matrix M.\n\n\n\n"
},

{
    "location": "matrix.html#Rank-1",
    "page": "Generic matrices",
    "title": "Rank",
    "category": "section",
    "text": "rank{T <: RingElem}(::MatElem{T})\nrank{T <: FieldElem}(::MatElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\n   \nA = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])\n\nd = rank(A)"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.solve-Union{Tuple{AbstractAlgebra.MatElem{T},AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.FieldElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.solve",
    "category": "Method",
    "text": "solve{T <: FieldElement}(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T})\n\nGiven a non-singular ntimes n matrix over a field and an ntimes m matrix over the same field, return x an ntimes m matrix x such that Ax = b. If A is singular an exception is raised.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.solve_rational-Union{Tuple{AbstractAlgebra.MatElem{T},AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.solve_rational",
    "category": "Method",
    "text": "solve_rational{T <: RingElement}(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T})\n\nGiven a non-singular ntimes n matrix over a ring and an ntimes m matrix over the same ring, return a tuple x d consisting of an ntimes m matrix x and a denominator d such that Ax = db. The denominator will be the determinant of A up to sign. If A is singular an exception is raised.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.solve_triu-Union{Tuple{AbstractAlgebra.MatElem{T},AbstractAlgebra.MatElem{T},Bool}, Tuple{T}} where T<:AbstractAlgebra.FieldElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.solve_triu",
    "category": "Method",
    "text": "solve_triu{T <: FieldElement}(U::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}, unit=false)\n\nGiven a non-singular ntimes n matrix over a field which is upper triangular, and an ntimes m matrix over the same field, return an ntimes m matrix x such that Ax = b. If A is singular an exception is raised. If unit is true then U is assumed to have ones on its diagonal, and the diagonal will not be read.\n\n\n\n"
},

{
    "location": "matrix.html#Linear-solving-1",
    "page": "Generic matrices",
    "title": "Linear solving",
    "category": "section",
    "text": "solve{T <: FieldElem}(::MatElem{T}, ::MatElem{T})solve_rational{T <: RingElem}(::MatElem{T}, ::MatElem{T})solve_triu{T <: FieldElem}(::MatElem{T}, ::MatElem{T}, ::Bool)ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\nU = MatrixSpace(K, 3, 1)\n\nA = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])\nb = U([2a a + 1 (-a - 1)]')\n\nx = solve(A, b)\n\nA = S([a + 1 2a + 3 a^2 + 1; K(0) a^2 - 1 2a; K(0) K(0) a])\nb = U([2a a + 1 (-a - 1)]')\n\nx = solve_triu(A, b, false)\n\nR, x = PolynomialRing(ZZ, \"x\")\nS = MatrixSpace(R, 3, 3)\nU = MatrixSpace(R, 3, 2)\n\nA = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])\nb = U([2x x + 1 (-x - 1); x + 1 (-x) x^2]')\n\nx, d = solve_rational(A, b)\n\nS = MatrixSpace(ZZ, 3, 3)\nT = MatrixSpace(ZZ, 3, 1)\n\nA = S([BigInt(2) 3 5; 1 4 7; 9 2 2])   \nB = T([BigInt(4), 5, 7])\n\nX, d = solve_rational(A, B)"
},

{
    "location": "matrix.html#Base.inv-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "Base.inv",
    "category": "Method",
    "text": "inv(M)\n\nMatrix inverse. Computes matrix N such that M * N = I, where I is the identity matrix. Computed by solving the left-division N = M \\ I.\n\nExample\n\njulia> M = [2 5; 1 3]\n2×2 Array{Int64,2}:\n 2  5\n 1  3\n\njulia> N = inv(M)\n2×2 Array{Float64,2}:\n  3.0  -5.0\n -1.0   2.0\n\njulia> M*N == N*M == eye(2)\ntrue\n\n\n\ninv{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nGiven a non-singular ntimes n matrix over a ring the tuple X d consisting of an ntimes n matrix X and a denominator d such that AX = dI_n, where I_n is the ntimes n identity matrix. The denominator will be the determinant of A up to sign. If A is singular an exception is raised.\n\n\n\n"
},

{
    "location": "matrix.html#Base.inv-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.FieldElem",
    "page": "Generic matrices",
    "title": "Base.inv",
    "category": "Method",
    "text": "inv(M)\n\nMatrix inverse. Computes matrix N such that M * N = I, where I is the identity matrix. Computed by solving the left-division N = M \\ I.\n\nExample\n\njulia> M = [2 5; 1 3]\n2×2 Array{Int64,2}:\n 2  5\n 1  3\n\njulia> N = inv(M)\n2×2 Array{Float64,2}:\n  3.0  -5.0\n -1.0   2.0\n\njulia> M*N == N*M == eye(2)\ntrue\n\n\n\ninv{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nGiven a non-singular ntimes n matrix over a ring the tuple X d consisting of an ntimes n matrix X and a denominator d such that AX = dI_n, where I_n is the ntimes n identity matrix. The denominator will be the determinant of A up to sign. If A is singular an exception is raised.\n\n\n\ninv{T <: FieldElement}(M::AbstractAlgebra.MatElem{T})\n\nGiven a non-singular ntimes n matrix over a field, return an ntimes n matrix X such that AX = I_n where I_n is the ntimes n identity matrix. If A is singular an exception is raised.\n\n\n\n"
},

{
    "location": "matrix.html#Inverse-1",
    "page": "Generic matrices",
    "title": "Inverse",
    "category": "section",
    "text": "inv{T <: RingElem}(::MatElem{T})\ninv{T <: FieldElem}(::MatElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\n\nA = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])\n\nX = inv(A)\n\nR, x = PolynomialRing(ZZ, \"x\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])\n    \nX, d = inv(A)"
},

{
    "location": "matrix.html#Base.LinAlg.nullspace-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "Base.LinAlg.nullspace",
    "category": "Method",
    "text": "nullspace(M)\n\nBasis for nullspace of M.\n\nExample\n\njulia> M = [1 0 0; 0 1 0; 0 0 0]\n3×3 Array{Int64,2}:\n 1  0  0\n 0  1  0\n 0  0  0\n\njulia> nullspace(M)\n3×1 Array{Float64,2}:\n 0.0\n 0.0\n 1.0\n\n\n\nnullspace{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nReturns a tuple (nu N) consisting of the nullity nu of M and a basis N (consisting of column vectors) for the right nullspace of M, i.e. such that MN is the zero matrix. If M is an mtimes n matrix N will be an ntimes nu matrix. Note that the nullspace is taken to be the vector space kernel over the fraction field of the base ring if the latter is not a field. In AbstractAlgebra we use the name ``kernel'' for a function to compute an integral kernel.\n\n\n\n"
},

{
    "location": "matrix.html#Base.LinAlg.nullspace-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.FieldElem",
    "page": "Generic matrices",
    "title": "Base.LinAlg.nullspace",
    "category": "Method",
    "text": "nullspace(M)\n\nBasis for nullspace of M.\n\nExample\n\njulia> M = [1 0 0; 0 1 0; 0 0 0]\n3×3 Array{Int64,2}:\n 1  0  0\n 0  1  0\n 0  0  0\n\njulia> nullspace(M)\n3×1 Array{Float64,2}:\n 0.0\n 0.0\n 1.0\n\n\n\nnullspace{T <: RingElement}(M::AbstractAlgebra.MatElem{T})\n\nReturns a tuple (nu N) consisting of the nullity nu of M and a basis N (consisting of column vectors) for the right nullspace of M, i.e. such that MN is the zero matrix. If M is an mtimes n matrix N will be an ntimes nu matrix. Note that the nullspace is taken to be the vector space kernel over the fraction field of the base ring if the latter is not a field. In AbstractAlgebra we use the name ``kernel'' for a function to compute an integral kernel.\n\n\n\nnullspace{T <: FieldElement}(M::AbstractAlgebra.MatElem{T})\n\nReturns a tuple (nu N) consisting of the nullity nu of M and a basis N (consisting of column vectors) for the right nullspace of M, i.e. such that MN is the zero matrix. If M is an mtimes n matrix N will be an ntimes nu matrix. Note that the nullspace is taken to be the vector space kernel over the fraction field of the base ring if the latter is not a field. In Nemo we use the name ``kernel'' for a function to compute an integral kernel.\n\n\n\n"
},

{
    "location": "matrix.html#Nullspace-1",
    "page": "Generic matrices",
    "title": "Nullspace",
    "category": "section",
    "text": "nullspace{T <: RingElem}(::MatElem{T})\nnullspace{T <: FieldElem}(::MatElem{T})ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS = MatrixSpace(R, 4, 4)\n   \nM = S([-6*x^2+6*x+12 -12*x^2-21*x-15 -15*x^2+21*x+33 -21*x^2-9*x-9;\n       -8*x^2+8*x+16 -16*x^2+38*x-20 90*x^2-82*x-44 60*x^2+54*x-34;\n       -4*x^2+4*x+8 -8*x^2+13*x-10 35*x^2-31*x-14 22*x^2+21*x-15;\n       -10*x^2+10*x+20 -20*x^2+70*x-25 150*x^2-140*x-85 105*x^2+90*x-50])\n   \nn, N = nullspace(M)"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.hessenberg-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.hessenberg",
    "category": "Method",
    "text": "hessenberg(A::AbstractAlgebra.MatElem{T}) where {T <: RingElement}\n\nReturns the Hessenberg form of M, i.e. an upper Hessenberg matrix which is similar to M. The upper Hessenberg form has nonzero entries above and on the diagonal and in the diagonal line immediately below the diagonal.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.ishessenberg-Union{Tuple{AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.ishessenberg",
    "category": "Method",
    "text": "ishessenberg{T <: RingElement}(A::AbstractAlgebra.MatElem{T})\n\nReturns true if M is in Hessenberg form, otherwise returns false.\n\n\n\n"
},

{
    "location": "matrix.html#Hessenberg-form-1",
    "page": "Generic matrices",
    "title": "Hessenberg form",
    "category": "section",
    "text": "hessenberg{T <: RingElem}(::MatElem{T})ishessenberg{T <: RingElem}(::MatElem{T})ExamplesR = ResidueRing(ZZ, 7)\nS = MatrixSpace(R, 4, 4)\n   \nM = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);\n       R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])\n   \nA = hessenberg(M)\nishessenberg(A) == true"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.charpoly-Union{Tuple{AbstractAlgebra.Ring,AbstractAlgebra.MatElem{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.charpoly",
    "category": "Method",
    "text": "charpoly{T <: RingElement}(V::Ring, Y::AbstractAlgebra.MatElem{T})\n\nReturns the characteristic polynomial p of the matrix M. The polynomial ring R of the resulting polynomial must be supplied and the matrix is assumed to be square.\n\n\n\n"
},

{
    "location": "matrix.html#Characteristic-polynomial-1",
    "page": "Generic matrices",
    "title": "Characteristic polynomial",
    "category": "section",
    "text": "charpoly{T <: RingElem}(::Ring, ::MatElem{T})ExamplesR = ResidueRing(ZZ, 7)\nS = MatrixSpace(R, 4, 4)\nT, x = PolynomialRing(R, \"x\")\n\nM = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);\n       R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])\n   \nA = charpoly(T, M)"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.minpoly-Union{Tuple{AbstractAlgebra.Ring,AbstractAlgebra.MatElem{T},Bool}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.minpoly",
    "category": "Method",
    "text": "minpoly{T <: RingElement}(S::Ring, M::AbstractAlgebra.MatElem{T}, charpoly_only = false)\n\nReturns the minimal polynomial p of the matrix M. The polynomial ring R of the resulting polynomial must be supplied and the matrix must be square.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.minpoly-Union{Tuple{AbstractAlgebra.Ring,AbstractAlgebra.MatElem{T},Bool}, Tuple{T}} where T<:AbstractAlgebra.FieldElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.minpoly",
    "category": "Method",
    "text": "minpoly{T <: FieldElement}(S::Ring, M::AbstractAlgebra.MatElem{T}, charpoly_only = false)\n\nReturns the minimal polynomial p of the matrix M. The polynomial ring R of the resulting polynomial must be supplied and the matrix must be square.\n\n\n\nminpoly{T <: RingElement}(S::Ring, M::AbstractAlgebra.MatElem{T}, charpoly_only = false)\n\nReturns the minimal polynomial p of the matrix M. The polynomial ring R of the resulting polynomial must be supplied and the matrix must be square.\n\n\n\n"
},

{
    "location": "matrix.html#Minimal-polynomial-1",
    "page": "Generic matrices",
    "title": "Minimal polynomial",
    "category": "section",
    "text": "minpoly{T <: RingElem}(::Ring, ::MatElem{T}, ::Bool)\nminpoly{T <: FieldElem}(::Ring, ::MatElem{T}, ::Bool)ExamplesR = GF(13)\nT, y = PolynomialRing(R, \"y\")\n   \nM = R[7 6 1;\n      7 7 5;\n      8 12 5]\n\nA = minpoly(T, M)"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.similarity!-Union{Tuple{AbstractAlgebra.MatElem{T},Int64,T}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.similarity!",
    "category": "Method",
    "text": "similarity!{T <: RingElement}(A::AbstractAlgebra.MatElem{T}, r::Int, d::T)\n\nApplies a similarity transform to the ntimes n matrix M in-place. Let P be the ntimes n identity matrix that has had all zero entries of row r replaced with d, then the transform applied is equivalent to M = P^-1MP. We require M to be a square matrix. A similarity transform preserves the minimal and characteristic polynomials of a matrix.\n\n\n\n"
},

{
    "location": "matrix.html#Transforms-1",
    "page": "Generic matrices",
    "title": "Transforms",
    "category": "section",
    "text": "similarity!{T <: RingElem}(::MatElem{T}, ::Int, ::T)ExamplesR = ResidueRing(ZZ, 7)\nS = MatrixSpace(R, 4, 4)\n   \nM = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);\n       R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])\n   \nsimilarity!(M, 1, R(3))"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.weak_popov-Union{Tuple{AbstractAlgebra.Generic.Mat{T}}, Tuple{T}} where T<:AbstractAlgebra.PolyElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.weak_popov",
    "category": "Method",
    "text": "weak_popov{T <: PolyElem}(A::Mat{T})\n\nReturn the weak Popov form of A.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.weak_popov_with_trafo-Union{Tuple{AbstractAlgebra.Generic.Mat{T}}, Tuple{T}} where T<:AbstractAlgebra.PolyElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.weak_popov_with_trafo",
    "category": "Method",
    "text": "weak_popov_with_trafo{T <: PolyElem}(A::Mat{T})\n\nCompute a tuple (P U) where P is the weak Popov form of A and U is a transformation matrix so that P = UA.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.popov-Union{Tuple{AbstractAlgebra.Generic.Mat{T}}, Tuple{T}} where T<:AbstractAlgebra.PolyElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.popov",
    "category": "Method",
    "text": "popov{T <: PolyElem}(A::Mat{T})\n\nReturn the Popov form of A.\n\n\n\n"
},

{
    "location": "matrix.html#AbstractAlgebra.Generic.popov_with_trafo-Union{Tuple{AbstractAlgebra.Generic.Mat{T}}, Tuple{T}} where T<:AbstractAlgebra.PolyElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.popov_with_trafo",
    "category": "Method",
    "text": "popov_with_trafo{T <: PolyElem}(A::Mat{T})\n\nCompute a tuple (P U) where P is the Popov form of A and U is a transformation matrix so that P = UA.\n\n\n\n"
},

{
    "location": "matrix.html#(Weak)-Popov-form-1",
    "page": "Generic matrices",
    "title": "(Weak) Popov form",
    "category": "section",
    "text": "AbstractAlgebra.jl provides algorithms for computing the (weak) Popov of a matrix with entries in a univariate polynomial ring over a field.weak_popov{T <: PolyElem}(::Generic.Mat{T})\nweak_popov_with_trafo{T <: PolyElem}(::Generic.Mat{T})\npopov{T <: PolyElem}(::Generic.Mat{T})\npopov_with_trafo{T <: PolyElem}(::Generic.Mat{T})"
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
    "text": "Naively, one might have expected that rings in AbstractAlgebra.jl could be modeled as types and their elements as objects with the given type. But there are various reasons why this is not a good model.Consider the ring R = mathbbZnmathbbZ for a multiprecision integer n. If we were to model the ring R as a type, then the type would somehow need to contain the modulus n. This is not possible in Julia, and in fact it is not desirable, since the compiler would then recompile all the associated functions every time a different modulus n was used.We could attach the modulus n to the objects representing elements of the ring, rather than their type.But now we cannot create new elements of the ring mathbbZnmathbbZ given only their type, since the type no longer contains the modulus n.Instead, the way we get around this in AbstractAlgebra.jl is to have special (singleton) objects that act like types, but are really just ordinary Julia objects. These objects, called parent objects can contain extra information, such as the modulus n. In order to create new elements of mathbbZnmathbbZ as above, we overload the call operator for the parent object.In the following AbstractAlgebra.jl example, we create the parent object R corresponding to the ring mathbbZ7mathbbZ. We then create a new element a of this ring by calling the parent object R.R = ResidueRing(ZZ, 7)\na = R(3)Here, R is the parent object, containing the modulus 7. So this example creates  the element a = 3 pmod7."
},

{
    "location": "types.html#More-complex-example-of-parent-objects-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "More complex example of parent objects",
    "category": "section",
    "text": "Here is some Julia/AbstractAlgebra.jl code which constructs a polynomial ring over the integers, a polynomial in that ring and then does some introspection to illustrate the various relations between the objects and types.using AbstractAlgebra\n\nR, x = ZZ[\"x\"]\n\nf = x^2 + 3x + 1\n\ntypeof(R) <: PolyRing\n\ntypeof(f) <: PolyElem\n\nparent(f) == R"
},

{
    "location": "types.html#Concrete-types-in-AbstractAlgebra.jl-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "Concrete types in AbstractAlgebra.jl",
    "category": "section",
    "text": "Here we give a list of the concrete types in AbstractAlgebra.jl.In parentheses we put the types of the corresponding parent objects.gfelem{<:Integer} (GFField{<:Integer})We also think of various Julia types as though they were AbstractAlgebra.jl types:BigInt (Integers{BigInt})\nRational{BigInt} (Rationals{BigInt})Then there are various types for generic constructions over a base ring. They are all parameterised by a type T which is the type of the elements of the base ring they are defined over. Generic.Poly{T} (Generic.PolyRing{T})\nGeneric.MPoly{T} (Generic.MPolyRing{T})\nGeneric.RelSeries{T} (Generic.RelSeriesRing{T})\nGeneric.AbsSeries{T} (Generic.AbsSeriesRing{T})\nGeneric.LaurentSeriesRingElem{T} (Generic.LaurentSeriesRing{T})\nGeneric.LaurentSeriesFieldElem{T} (Generic.LaurentSeriesField{T})\nGeneric.Res{T} (Generic.ResRing{T})\nGeneric.Frac{T} (Generic.FracField{T})\nGeneric.Mat{T} (Generic.MatSpace{T})"
},

]}
